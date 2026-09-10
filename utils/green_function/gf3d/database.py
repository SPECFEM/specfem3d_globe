"""The database handle and what comes out of it."""

from __future__ import annotations

import ctypes
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
_DOUBLE_P = ctypes.POINTER(ctypes.c_double)


def _ptr(a: np.ndarray):
    return a.ctypes.data_as(_DOUBLE_P)


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
    t0_req: float
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
    """

    def __init__(self, path, check_completion: bool = False):
        self.path = str(path)
        self._handle = ctypes.c_int(0)
        self._info = None

        with LIBRARY_LOCK:
            code = lib.gf3d_open(
                self.path.encode(), 1 if check_completion else 0, ctypes.byref(self._handle)
            )
            check(code, f"opening {self.path}")

    # -- lifetime ---------------------------------------------------------

    @property
    def closed(self) -> bool:
        return self._handle.value == 0

    def close(self) -> None:
        """Close the database and release the search tree. Idempotent."""
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
        return self._extract(source, t0, itypsokern=0)

    def partials(self, source, kind: int = 2, t0: float | None = None) -> Result:
        """Seismograms and their partial derivatives.

        ``kind=1`` gives the six moment-tensor partials, ``kind=2`` those
        plus latitude, longitude, depth and centroid time. A force source
        has none.
        """
        if kind not in (1, 2):
            raise ValueError("kind must be 1 (moment tensor) or 2 (and centroid)")
        return self._extract(source, t0, itypsokern=kind)

    def _extract(self, source, t0, itypsokern: int) -> Result:
        csrc = source._to_struct()
        t0_req = -1.0 if t0 is None else float(t0)

        plan = self.plan(source, t0)
        nt = plan.nt
        nsta = self.info["nstations"]
        ids = self.station_ids

        ndp = 0
        if itypsokern > 0:
            n = ctypes.c_int(0)
            with LIBRARY_LOCK:
                check(lib.gf3d_ndp(itypsokern, ctypes.byref(n)), "gf3d_ndp")
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
                    h, ctypes.byref(csrc), t0_req, itypsokern, nt, ndp,
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


def open(path, check_completion: bool = False) -> Database:  # noqa: A001
    """Open a database. The same as calling :class:`Database`."""
    return Database(path, check_completion=check_completion)
