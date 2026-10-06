"""The database handle and what comes out of it."""

from __future__ import annotations

import ctypes
import numbers
import warnings
from dataclasses import dataclass

import numpy as np

from . import _lib
from ._lib import (
    CInfo,
    CLocation,
    CPlan,
    CStation,
    GF_NCOMP,
    GF_NDP_LOC,
    GF_NDP_MT,
    GF3DError,
    LIBRARY_LOCK,
    check,
    lib,
)
from .sources import CMTSource, ForceSource

__all__ = ["Database", "Result", "Plan", "Location", "Station", "open"]

_ONSET_WARN = 1.0e-3
_INT_MAX = 2**31 - 1
_DOUBLE_P = ctypes.POINTER(ctypes.c_double)
_FLOAT_P = ctypes.POINTER(ctypes.c_float)
_INT_P = ctypes.POINTER(ctypes.c_int)
_NGLL3 = 125


def _ptr(a: np.ndarray):
    return a.ctypes.data_as(_DOUBLE_P)


def _cint(n) -> int:
    """An int clipped to the C int range, so that a silly one stays silly."""
    return min(max(int(n), -_INT_MAX - 1), _INT_MAX)


def _fptr(a: np.ndarray):
    return a.ctypes.data_as(_FLOAT_P)


def _s(b: bytes) -> str:
    return b.decode("utf-8", "replace").strip()


@dataclass(frozen=True)
class Station:
    """One reciprocal station of a database."""

    id: str
    network: str
    station: str
    latitude: float
    longitude: float
    depth_m: float
    hdur: float
    f_cutoff: float
    factor_force_source: float
    time_shift: float

    @classmethod
    def _from_c(cls, c: CStation) -> "Station":
        return cls(
            id=_s(c.id), network=_s(c.network), station=_s(c.station),
            latitude=c.latitude, longitude=c.longitude, depth_m=c.depth_m,
            hdur=c.hdur, f_cutoff=c.f_cutoff,
            factor_force_source=c.factor_force_source, time_shift=c.time_shift,
        )


@dataclass(frozen=True)
class Location:
    """Where a source sits in the mesh."""

    ielem: int
    morton_hex: str
    xi: float
    eta: float
    gamma: float
    xyz: tuple
    xyz_target: tuple
    distance_km: float
    anchor_err: float
    theta: float
    phi: float
    r_surface: float

    @classmethod
    def _from_c(cls, c: CLocation) -> "Location":
        return cls(
            ielem=c.ielem, morton_hex=_s(c.morton_hex),
            xi=c.xi, eta=c.eta, gamma=c.gamma,
            xyz=tuple(c.xyz), xyz_target=tuple(c.xyz_target),
            distance_km=c.distance_km, anchor_err=c.anchor_err,
            theta=c.theta, phi=c.phi, r_surface=c.r_surface,
        )


@dataclass(frozen=True)
class Plan:
    """The output time axis and the source time function conversion.

    ``nt`` is the length of everything an extraction returns, and
    ``t_first`` its first sample, in seconds relative to the centroid time.
    ``guard`` is set when the requested source is narrower than the one the
    database was built with, in which case no width correction is possible
    and the traces carry the database's own duration.
    """

    nt: int
    nt_db: int
    npad: int
    subsample_step: int
    khalf: int
    guard: bool
    kind_stf: int
    dt: float
    dt_sub: float
    t0_db: float
    t0_req: float          # resolved: never the caller's negative sentinel
    t0: float
    t_first: float
    hdur_src: float
    hdur_target: float
    hdur_db: float
    hdur_corr: float
    trunc: float

    @classmethod
    def _from_c(cls, c: CPlan) -> "Plan":
        return cls(
            nt=c.nt, nt_db=c.nt_db, npad=c.npad, subsample_step=c.subsample_step,
            khalf=c.khalf, guard=bool(c.guard), kind_stf=c.kind_stf,
            dt=c.dt, dt_sub=c.dt_sub, t0_db=c.t0_db, t0_req=c.t0_req, t0=c.t0,
            t_first=c.t_first, hdur_src=c.hdur_src, hdur_target=c.hdur_target,
            hdur_db=c.hdur_db, hdur_corr=c.hdur_corr, trunc=c.trunc,
        )

    @property
    def times(self) -> np.ndarray:
        """The axis, as the library builds it."""
        return self.t_first + np.arange(self.nt) * self.dt_sub


@dataclass
class Result:
    """Seismograms, and optionally their partial derivatives.

    ``data`` has shape ``(nstations, 3, nt)`` in metres, components N, E, Z,
    with ``t`` the shared time axis in seconds relative to the centroid
    time. ``dp``, when present, has shape ``(nstations, ndp, 3, nt)`` with
    the partials in the order named by ``dp_names``: the six moment-tensor
    components per dyne-cm, then latitude and longitude per degree, depth
    per km and centroid time per second.

    So ``(result.dp[:, :6] * cmt.tensor[None, :, None, None]).sum(1)``
    reproduces ``result.data``: the partials are exact and linear, not
    finite differences.
    """

    t: np.ndarray
    data: np.ndarray
    onset: np.ndarray
    plan: Plan
    station_ids: list
    location: Location | None = None
    dp: np.ndarray | None = None
    dp_names: list | None = None
    dp_units: list | None = None
    source: object = None

    @property
    def nt(self) -> int:
        return self.data.shape[2]

    @property
    def nstations(self) -> int:
        return self.data.shape[0]

    def trace(self, station, component: str) -> np.ndarray:
        """One trace by station id (or index) and component name."""
        ista = self.station_ids.index(station) if isinstance(station, str) else station
        icomp = "NEZ".index(component.upper())
        return self.data[ista, icomp]

    def partial(self, station, component: str, name: str) -> np.ndarray:
        """One partial derivative by name, e.g. ``"Mrr"`` or ``"dep"``."""
        if self.dp is None:
            raise ValueError("this result carries no partial derivatives")
        ista = self.station_ids.index(station) if isinstance(station, str) else station
        icomp = "NEZ".index(component.upper())
        return self.dp[ista, self.dp_names.index(name), icomp]

    def to_stream(self):
        """An obspy :class:`~obspy.core.stream.Stream` of the seismograms.

        obspy is an optional dependency; this is the only thing that needs
        it. Channels follow the solver's own naming, ``BXN``/``BXE``/``BXZ``,
        and the start time is the centroid time plus ``t[0]`` when the
        source carried an origin time.
        """
        from obspy import Stream, Trace, UTCDateTime  # noqa: PLC0415
        from obspy.core.trace import Stats  # noqa: PLC0415

        centroid = getattr(self.source, "centroid_time", None)
        traces = []
        for ista, sid in enumerate(self.station_ids):
            net, _, sta = sid.partition(".")
            for icomp, comp in enumerate("NEZ"):
                stats = Stats()
                stats.network = net
                stats.station = sta
                stats.channel = f"BX{comp}"
                stats.location = ""
                stats.delta = self.plan.dt_sub
                if centroid is not None:
                    stats.starttime = UTCDateTime(centroid) + float(self.t[0])
                else:
                    stats.starttime = UTCDateTime(0) + float(self.t[0])
                traces.append(Trace(data=self.data[ista, icomp].copy(), header=stats))
        return Stream(traces)


class Database:
    """An open Green function database.

    Use it as a context manager, or call :meth:`close` when done::

        with gf3d.Database("EXAMPLES/.../regional/GFDB") as db:
            r = db.partials(gf3d.CMTSource.read("CMTSOLUTION"))

    Opening reads the metadata only; the 58 MB topography grid and the
    element files are read on demand, so the cost of an open is small and
    the cost of the first extraction is not.

    Two databases may be open at once, but the search tree is rebuilt
    whenever an extraction switches between them, so alternating is slow --
    and two databases of *different planets or topography grids* must not be
    used alternately at all, since specfem's globals hold one set of
    dimensions for the process.

    ``max_elements`` is how many elements the handle keeps in memory between
    extractions, for a caller that comes back to them -- a sampler walking
    around one source, say. Each takes ``info["bytes_per_element"]`` (every
    station's displacement for that element), and the least recently used
    one is dropped to make room; a position inside a kept element reads no
    displacement file. Every handle, whatever ``max_elements``, also keeps
    the coordinates its locates read (3 kB per element) for the
    ``max(10, max_elements)`` elements used most recently, so returning to a
    position still held reads nothing from disk. The numbers are the same
    with or without either, to the bit. ``0``, the default, keeps no
    elements and reads every extraction's element from disk.
    :attr:`cache_stats` says what it did.

    The library is not thread-safe; calls are serialised through one lock.
    Parallel chains belong in separate processes, each with its own cache,
    so budget ``processes * max_elements * bytes_per_element``.
    """

    def __init__(self, path, check_completion: bool = False, max_elements: int = 0):
        self.path = str(path)
        self._handle = ctypes.c_int(0)
        self._info = None

        # checked here: ctypes would silently wrap a Python int into a C int
        if isinstance(max_elements, bool) or not isinstance(max_elements, numbers.Integral):
            raise TypeError(f"max_elements must be an integer, not {type(max_elements).__name__}")
        if max_elements < 0:
            raise ValueError(f"max_elements must not be negative, got {max_elements}")
        # more than the database holds means all of them, as in the library
        self._max_elements = int(min(max_elements, _INT_MAX))

        with LIBRARY_LOCK:
            code = lib.gf3d_open(
                self.path.encode(), 1 if check_completion else 0, self._max_elements,
                ctypes.byref(self._handle),
            )
            check(code, f"opening {self.path}")

    # -- lifetime ---------------------------------------------------------

    @property
    def closed(self) -> bool:
        return self._handle.value == 0

    def close(self) -> None:
        """Close the database, freeing its element cache, and release the
        search tree. Idempotent."""
        if self.closed:
            return
        h, self._handle.value = self._handle.value, 0
        with LIBRARY_LOCK:
            lib.gf3d_close(h)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False

    def __del__(self):
        # interpreter teardown can have unloaded ctypes already
        try:
            self.close()
        except Exception:
            pass

    def __repr__(self):
        if self.closed:
            return f"<gf3d.Database {self.path!r} (closed)>"
        return (
            f"<gf3d.Database {self.path!r}: {self.info['nstations']} stations, "
            f"{self.info['nelem']} elements, {self.info['nt_subsampled']} samples>"
        )

    def _h(self) -> int:
        if self.closed:
            raise GF3DError(_lib.GF_ERR_ARG, "invalid argument", "the database is closed")
        return self._handle.value

    # -- metadata ---------------------------------------------------------

    @property
    def info(self) -> dict:
        """What the database says about itself."""
        if self._info is None:
            c = CInfo()
            with LIBRARY_LOCK:
                check(lib.gf3d_get_info(self._h(), ctypes.byref(c)), "gf3d_get_info")
            d = {f: getattr(c, f) for f, _ in CInfo._fields_ if f != "pad_"}
            for flag in ("topography", "ellipticity", "rotation", "attenuation", "gravity"):
                d[flag] = bool(d[flag])
            self._info = d
        return self._info

    @property
    def max_elements(self) -> int:
        """How many elements the handle may keep, as asked for at open."""
        return self._max_elements

    @property
    def cache_stats(self) -> dict:
        """What the element cache has done since the database was opened.

        ``hits``, ``misses``: extractions whose element was, or was not,
        already in memory; with ``max_elements=0`` every one is a miss.
        ``evictions``: elements dropped to make room. ``n_cached``: elements
        held now. ``files_read``: element files, coordinates or
        displacement, this handle has read by any route. The coordinate
        store has no counters of its own: its reads show only here.
        """
        hits, misses, evictions, files = (ctypes.c_longlong() for _ in range(4))
        n_cached = ctypes.c_int()
        with LIBRARY_LOCK:
            check(
                lib.gf3d_cache_stats(
                    self._h(), ctypes.byref(hits), ctypes.byref(misses),
                    ctypes.byref(evictions), ctypes.byref(n_cached), ctypes.byref(files),
                ),
                "gf3d_cache_stats",
            )
        return {
            "hits": hits.value,
            "misses": misses.value,
            "evictions": evictions.value,
            "n_cached": n_cached.value,
            "files_read": files.value,
        }

    @property
    def stations(self) -> list:
        """Every station, in the library's own order."""
        # resolved before the lock is taken, not inside the loop: self.info
        # may itself have to call the library
        nstations = self.info["nstations"]

        out = []
        c = CStation()
        with LIBRARY_LOCK:
            h = self._h()
            for i in range(nstations):
                check(lib.gf3d_get_station(h, i, ctypes.byref(c)), f"station {i}")
                out.append(Station._from_c(c))
        return out

    @property
    def station_ids(self) -> list:
        """``['NET.STA', ...]``, the first axis of every extracted array."""
        return [s.id for s in self.stations]

    # -- locating and planning -------------------------------------------

    def locate(self, latitude: float, longitude: float, depth_km: float) -> Location:
        """Find the element containing a geographic point."""
        c = CLocation()
        with LIBRARY_LOCK:
            check(
                lib.gf3d_locate(self._h(), latitude, longitude, depth_km, ctypes.byref(c)),
                f"locating {latitude}, {longitude} at {depth_km} km",
            )
        return Location._from_c(c)

    def plan(self, source, t0: float | None = None) -> Plan:
        """The output axis and conversion for a source, without extracting."""
        c = CPlan()
        csrc = source._to_struct()
        with LIBRARY_LOCK:
            check(
                lib.gf3d_get_plan(
                    self._h(), ctypes.byref(csrc), -1.0 if t0 is None else float(t0),
                    ctypes.byref(c),
                ),
                "gf3d_get_plan",
            )
        return Plan._from_c(c)

    # -- extraction -------------------------------------------------------

    def seismograms(self, source, t0: float | None = None) -> Result:
        """Seismograms at every station.

        ``t0`` is the start time in seconds before the centroid time; the
        default is specfem's own rule for a forward run, 1.5 times the half
        duration for a moment tensor.
        """
        return self._extract(source, t0, kind=0)

    def partials(self, source, kind: int = 2, t0: float | None = None) -> Result:
        """Seismograms and their partial derivatives.

        ``kind=1`` gives the six moment-tensor partials, ``kind=2`` those
        plus latitude, longitude, depth and centroid time. A force source
        has none.
        """
        if kind not in (1, 2):
            raise ValueError("kind must be 1 (moment tensor) or 2 (and centroid)")
        return self._extract(source, t0, kind=kind)

    def _extract(self, source, t0, kind: int) -> Result:
        csrc = source._to_struct()
        t0_req = -1.0 if t0 is None else float(t0)

        plan = self.plan(source, t0)
        nt = plan.nt
        nsta = self.info["nstations"]
        ids = self.station_ids

        ndp = 0
        if kind > 0:
            n = ctypes.c_int(0)
            with LIBRARY_LOCK:
                check(lib.gf3d_ndp(kind, ctypes.byref(n)), "gf3d_ndp")
            ndp = n.value

        seis = np.empty((nsta, GF_NCOMP, nt), dtype=np.float64)
        t = np.empty(nt, dtype=np.float64)
        onset = np.empty(nsta, dtype=np.float64)
        dp = np.empty((nsta, ndp, GF_NCOMP, nt), dtype=np.float64) if ndp else None
        cloc = CLocation()

        with LIBRARY_LOCK:
            h = self._h()
            if ndp:
                code = lib.gf3d_partials(
                    h, ctypes.byref(csrc), t0_req, kind, nt, ndp,
                    _ptr(seis), _ptr(dp), _ptr(t), _ptr(onset), ctypes.byref(cloc),
                )
            else:
                code = lib.gf3d_seismograms(
                    h, ctypes.byref(csrc), t0_req, nt,
                    _ptr(seis), _ptr(t), _ptr(onset), ctypes.byref(cloc),
                )
            check(code, "extraction")

        names = units = None
        if ndp:
            names, units = self.partial_names(ndp)

        worst = float(onset.max()) if nsta else 0.0
        if worst > _ONSET_WARN:
            warnings.warn(
                f"the record starts with {worst:.2e} of the trace peak already "
                "present: the source time function conversion is reaching back "
                "past the beginning of the stored record, so the first arrivals "
                "may be contaminated. A longer reciprocal run, or a later t0, "
                "is the cure.",
                RuntimeWarning,
                stacklevel=3,
            )

        return Result(
            t=t, data=seis, onset=onset, plan=plan, station_ids=ids,
            location=Location._from_c(cloc), dp=dp, dp_names=names,
            dp_units=units, source=source,
        )

    # -- one element's raw data --------------------------------------------

    def _station_indices(self, stations) -> np.ndarray:
        """0-based station indices from ids (``'NET.STA'``) and/or ints.

        An id that is not in the database is a ValueError naming it. An int
        is passed on as it is, clipped to the C int range only so that it
        arrives as an out-of-range index rather than as an OverflowError:
        the library does the range check.
        """
        ids = None
        out = []
        for s in stations:
            if isinstance(s, str):
                if ids is None:
                    ids = self.station_ids
                try:
                    out.append(ids.index(s))
                except ValueError:
                    raise ValueError(f"no station {s!r} in this database") from None
            else:
                out.append(_cint(s))
        return np.array(out, dtype=np.intc)

    def element_block(self, ielem: int, stations=None, nt: int | None = None,
                      out: np.ndarray | None = None) -> np.ndarray:
        """One element's raw displacement for a set of stations.

        Returns a C-contiguous float32 array of shape ``(nsel, 3, nt, 375)``:
        ``block[s, a, t, m]`` is the displacement of station ``s`` along the
        force component ``a`` (N, E, Z) at stored sample ``t``, with
        ``m = p + 3*(i + 5*j + 25*k)``, ``p`` the displacement component at
        the source and ``i, j, k`` the 0-based GLL indices. The numbers are
        the stored ones, unchanged: contract them with weights outside the
        library.

        ``ielem`` is 1-based, as :attr:`Location.ielem`. ``stations`` is
        ``None`` for every station, or a sequence of indices (0-based) and/or
        ids (``'NET.STA'``), in any order, repeats allowed. ``nt`` is how many
        samples from the first stored one to read, 1 to ``nt_subsampled``
        (default all); only the chunks that hold them are read. A trace
        contracted from a prefix is the full read's on every sample, but one
        converted with the source time function afterwards differs in its
        last ``plan.khalf`` samples, so read that many more than are kept.

        ``out``, if given, receives the block and is returned: a float32
        C-contiguous ndarray of exactly that shape.

        The element is read from disk on every call. The handle's element
        cache is neither used nor filled (``hits``, ``misses``, ``evictions``
        and ``n_cached`` do not move); ``files_read`` grows by one per
        station.
        """
        if nt is None:
            nt = self.info["nt_subsampled"]
        nt = _cint(nt)

        if stations is None:
            sel = None
            nsel = self.info["nstations"]
        else:
            sel = self._station_indices(stations)
            nsel = len(sel)

        # the library refuses a size it does not accept before the shape
        # below is allocated or indexed with it
        if nsel < 1 or nt < 1:
            shape = None
        else:
            shape = (nsel, GF_NCOMP, nt, GF_NCOMP * _NGLL3)

        if out is not None:
            if not (isinstance(out, np.ndarray) and out.dtype == np.float32
                    and out.flags["C_CONTIGUOUS"] and out.flags["WRITEABLE"]
                    and out.shape == shape):
                raise ValueError(
                    f"out must be a writeable C-contiguous float32 ndarray of shape {shape}"
                )
        elif shape is not None and nt <= self.info["nt_subsampled"]:
            out = np.empty(shape, dtype=np.float32)

        with LIBRARY_LOCK:
            check(
                lib.gf3d_element_block(
                    self._h(), _cint(ielem), nsel,
                    None if sel is None else sel.ctypes.data_as(_INT_P),
                    nt, None if out is None else _fptr(out),
                ),
                "gf3d_element_block",
            )
        return out

    @staticmethod
    def partial_names(ndp: int = GF_NDP_LOC):
        """``(names, units)`` of the first ``ndp`` partials."""
        names, units = [], []
        nbuf = ctypes.create_string_buffer(16)
        ubuf = ctypes.create_string_buffer(16)
        with LIBRARY_LOCK:
            for ip in range(ndp):
                check(
                    lib.gf3d_partial_name(ip, nbuf, len(nbuf), ubuf, len(ubuf)),
                    f"partial {ip}",
                )
                names.append(_s(nbuf.value))
                units.append(_s(ubuf.value))
        return names, units


def open(path, check_completion: bool = False, max_elements: int = 0) -> Database:  # noqa: A001
    """Open a database. The same as calling :class:`Database`."""
    return Database(path, check_completion=check_completion, max_elements=max_elements)
