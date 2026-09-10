"""ctypes binding to lib/libgf3d.so.

This module is the transcription of ``include/gf3d.h`` and nothing else: the
structures, the entry points, the library search, and the translation of a
status code into an exception. Everything with an opinion about how to use
them lives in :mod:`gf3d.database`.

The header is the contract. Two things guard the transcription against
drifting away from it:

* ``gf3d_sizeof`` is called at load time and compared with
  :func:`ctypes.sizeof` for all five structures, so a field added on one
  side and not the other fails loudly at import rather than silently
  returning shifted numbers;
* every function is given an explicit ``argtypes`` and ``restype``, so a
  wrong argument is a ``ctypes.ArgumentError`` at the call, not a
  reinterpreted pointer.

Thread safety
-------------
The library is not thread-safe. It keeps the last error message, the
kd-tree that serves whichever database located most recently, and specfem's
``shared_parameters`` in process-wide storage. Every call made through this
module therefore takes ``LIBRARY_LOCK``. Extraction releases the GIL (it is
a C call), so other Python threads keep running; they just cannot be inside
the library at the same time.
"""

from __future__ import annotations

import ctypes
import ctypes.util
import os
import threading
from pathlib import Path

__all__ = [
    "GF3DError",
    "GF_OK",
    "GF_SRC_FORCE",
    "GF_SRC_CMT",
    "GF_NCOMP",
    "GF_NDP_MT",
    "GF_NDP_LOC",
    "CSource",
    "CInfo",
    "CStation",
    "CLocation",
    "CPlan",
    "LIBRARY_LOCK",
    "lib",
    "library_path",
    "check",
]

# ---------------------------------------------------------------------------
# constants, from include/gf3d.h
# ---------------------------------------------------------------------------

GF3D_STRLEN = 64
GF3D_MORTON_STRLEN = 24

GF_OK = 0
GF_ERR_NO_HDF5 = 1
GF_ERR_NO_PATH = 2
GF_ERR_NO_FILE = 3
GF_ERR_HDF5 = 4
GF_ERR_IO = 5
GF_ERR_FORMAT = 6
GF_ERR_MISMATCH = 7
GF_ERR_INCOMPLETE = 8
GF_ERR_ALLOC = 9
GF_ERR_ARG = 10
GF_ERR_NO_ELEMENT = 11
GF_ERR_GEOMETRY = 12

GF_SRC_FORCE = 1
GF_SRC_CMT = 2

GF_NCOMP = 3
GF_NDP_MT = 6
GF_NDP_LOC = 10

_c_double = ctypes.c_double
_c_int = ctypes.c_int
_c_char = ctypes.c_char


class GF3DError(Exception):
    """A libgf3d call that did not return ``GF_OK``.

    Carries the numeric code, the library's short name for it and the
    message the library left behind, so that a caller can branch on
    ``err.code`` (for instance ``GF_ERR_NO_ELEMENT`` when a trial source
    wanders outside the database during an inversion) rather than parse
    text.
    """

    def __init__(self, code: int, name: str, message: str):
        self.code = code
        self.name = name
        self.message = message
        super().__init__(f"{message} [{name}, code {code}]" if message else f"{name} (code {code})")


# ---------------------------------------------------------------------------
# the structures
# ---------------------------------------------------------------------------


class CSource(ctypes.Structure):
    """``gf3d_source``"""

    _fields_ = [
        ("source_type", _c_int),
        ("force_stf", _c_int),
        ("latitude", _c_double),
        ("longitude", _c_double),
        ("depth_km", _c_double),
        ("hdur", _c_double),
        ("time_shift", _c_double),
        ("moment", _c_double * 6),
        ("force_factor", _c_double),
        ("force_dir", _c_double * 3),
    ]


class CInfo(ctypes.Structure):
    """``gf3d_info``"""

    _fields_ = [
        ("nelem", _c_int),
        ("nstations", _c_int),
        ("nstep", _c_int),
        ("nt_subsampled", _c_int),
        ("subsample_step", _c_int),
        ("ngllx", _c_int),
        ("nglly", _c_int),
        ("ngllz", _c_int),
        ("topography", _c_int),
        ("ellipticity", _c_int),
        ("rotation", _c_int),
        ("attenuation", _c_int),
        ("gravity", _c_int),
        ("pad_", _c_int),
        ("dt", _c_double),
        ("t0", _c_double),
        ("r_planet", _c_double),
        ("rhoav", _c_double),
        ("scale_displ", _c_double),
    ]


class CStation(ctypes.Structure):
    """``gf3d_station``"""

    _fields_ = [
        ("id", _c_char * GF3D_STRLEN),
        ("network", _c_char * GF3D_STRLEN),
        ("station", _c_char * GF3D_STRLEN),
        ("latitude", _c_double),
        ("longitude", _c_double),
        ("depth_m", _c_double),
        ("hdur", _c_double),
        ("f_cutoff", _c_double),
        ("factor_force_source", _c_double),
        ("time_shift", _c_double),
    ]


class CLocation(ctypes.Structure):
    """``gf3d_location``"""

    _fields_ = [
        ("ielem", _c_int),
        ("morton_hex", _c_char * GF3D_MORTON_STRLEN),
        ("xi", _c_double),
        ("eta", _c_double),
        ("gamma", _c_double),
        ("xyz", _c_double * 3),
        ("xyz_target", _c_double * 3),
        ("distance_km", _c_double),
        ("anchor_err", _c_double),
        ("theta", _c_double),
        ("phi", _c_double),
        ("r_surface", _c_double),
    ]


class CPlan(ctypes.Structure):
    """``gf3d_plan``"""

    _fields_ = [
        ("nt", _c_int),
        ("nt_db", _c_int),
        ("npad", _c_int),
        ("subsample_step", _c_int),
        ("khalf", _c_int),
        ("guard", _c_int),
        ("kind_stf", _c_int),
        ("pad_", _c_int),
        ("dt", _c_double),
        ("dt_sub", _c_double),
        ("t0_db", _c_double),
        ("t0_req", _c_double),
        ("t0", _c_double),
        ("t_first", _c_double),
        ("hdur_src", _c_double),
        ("hdur_target", _c_double),
        ("hdur_db", _c_double),
        ("hdur_corr", _c_double),
        ("trunc", _c_double),
    ]


# ---------------------------------------------------------------------------
# finding the library
# ---------------------------------------------------------------------------

_HELP = """
Could not load libgf3d.

The library is built from a configured specfem3d_globe tree:

    ./configure --with-hdf5 HDF5_INC=<dir> HDF5_LIBS=-L<dir>
    make gf3d

which writes lib/libgf3d.so. Point GF3D_LIB at it if it is somewhere this
package cannot guess:

    export GF3D_LIB=/path/to/specfem3d_globe/lib/libgf3d.so

A load that fails naming libgfortran, libifcore, libmpi or libhdf5 means the
library was found but the compiler's or HDF5's runtime is not on
LD_LIBRARY_PATH: load the same modules the tree was built with.
"""


def _candidates():
    env = os.environ.get("GF3D_LIB")
    if env:
        yield Path(env)

    # utils/green_function/gf3d/_lib.py -> the repository root
    root = Path(__file__).resolve().parents[3]
    for name in ("libgf3d.so", "libgf3d.dylib"):
        yield root / "lib" / name

    found = ctypes.util.find_library("gf3d")
    if found:
        yield Path(found)


def _load():
    tried = []
    for cand in _candidates():
        if cand.exists() or not cand.is_absolute():
            try:
                return ctypes.CDLL(str(cand)), str(cand)
            except OSError as exc:  # found but unloadable: say why
                raise OSError(f"{cand}: {exc}\n{_HELP}") from exc
        tried.append(str(cand))
    raise OSError("looked in:\n  " + "\n  ".join(tried) + "\n" + _HELP)


lib, library_path = _load()

# Reentrant on purpose. The lock serialises calls into a library that keeps
# process-wide state, and that is all it is for -- but the wrapper's own
# properties compose (Database.stations needs Database.info, and either may
# have to call the library), and under a plain Lock any such nesting is a
# deadlock. One existed: Database.stations resolved nstations inside the
# lock, so asking a fresh database for its stations before anything had
# read its info hung the interpreter. That call is now hoisted out, and
# this is reentrant so that the next composition of two locked properties
# is merely slow to write rather than fatal to run. An RLock still
# serialises across threads, which is the property that matters.
LIBRARY_LOCK = threading.RLock()


# ---------------------------------------------------------------------------
# signatures
# ---------------------------------------------------------------------------

_c_double_p = ctypes.POINTER(_c_double)
_c_int_p = ctypes.POINTER(_c_int)

lib.gf3d_version.argtypes = [ctypes.c_char_p, _c_int]
lib.gf3d_version.restype = _c_int

lib.gf3d_sizeof.argtypes = [_c_int_p] * 5
lib.gf3d_sizeof.restype = _c_int

lib.gf3d_last_error.argtypes = [ctypes.c_char_p, _c_int]
lib.gf3d_last_error.restype = _c_int

lib.gf3d_error_string.argtypes = [_c_int, ctypes.c_char_p, _c_int]
lib.gf3d_error_string.restype = _c_int

lib.gf3d_open.argtypes = [ctypes.c_char_p, _c_int, _c_int_p]
lib.gf3d_open.restype = _c_int

lib.gf3d_close.argtypes = [_c_int]
lib.gf3d_close.restype = _c_int

lib.gf3d_get_info.argtypes = [_c_int, ctypes.POINTER(CInfo)]
lib.gf3d_get_info.restype = _c_int

lib.gf3d_get_station.argtypes = [_c_int, _c_int, ctypes.POINTER(CStation)]
lib.gf3d_get_station.restype = _c_int

lib.gf3d_locate.argtypes = [_c_int, _c_double, _c_double, _c_double, ctypes.POINTER(CLocation)]
lib.gf3d_locate.restype = _c_int

lib.gf3d_get_plan.argtypes = [_c_int, ctypes.POINTER(CSource), _c_double, ctypes.POINTER(CPlan)]
lib.gf3d_get_plan.restype = _c_int

lib.gf3d_ndp.argtypes = [_c_int, _c_int_p]
lib.gf3d_ndp.restype = _c_int

lib.gf3d_partial_name.argtypes = [_c_int, ctypes.c_char_p, _c_int, ctypes.c_char_p, _c_int]
lib.gf3d_partial_name.restype = _c_int

lib.gf3d_seismograms.argtypes = [
    _c_int,
    ctypes.POINTER(CSource),
    _c_double,
    _c_int,
    _c_double_p,
    _c_double_p,
    _c_double_p,
    ctypes.POINTER(CLocation),
]
lib.gf3d_seismograms.restype = _c_int

lib.gf3d_partials.argtypes = [
    _c_int,
    ctypes.POINTER(CSource),
    _c_double,
    _c_int,
    _c_int,
    _c_int,
    _c_double_p,
    _c_double_p,
    _c_double_p,
    _c_double_p,
    ctypes.POINTER(CLocation),
]
lib.gf3d_partials.restype = _c_int


# ---------------------------------------------------------------------------
# layout check
# ---------------------------------------------------------------------------


def _check_layout():
    sizes = [_c_int() for _ in range(5)]
    lib.gf3d_sizeof(*[ctypes.byref(s) for s in sizes])
    expected = [CSource, CInfo, CStation, CLocation, CPlan]
    for got, cls in zip(sizes, expected):
        if got.value != ctypes.sizeof(cls):
            raise RuntimeError(
                f"{cls.__name__} is {ctypes.sizeof(cls)} bytes here but "
                f"{got.value} in {library_path}: this package and the library "
                "were built from different versions of include/gf3d.h"
            )


_check_layout()


# ---------------------------------------------------------------------------
# errors
# ---------------------------------------------------------------------------


def check(code: int, context: str = "") -> None:
    """Raise :class:`GF3DError` unless ``code`` is ``GF_OK``.

    Must be called with ``LIBRARY_LOCK`` still held: the message it reads is
    process-wide, so another thread's failure would otherwise be reported
    here.
    """
    if code == GF_OK:
        return

    name_buf = ctypes.create_string_buffer(64)
    msg_buf = ctypes.create_string_buffer(512)
    lib.gf3d_error_string(code, name_buf, len(name_buf))
    lib.gf3d_last_error(msg_buf, len(msg_buf))

    message = msg_buf.value.decode("utf-8", "replace").strip()
    if context:
        message = f"{context}: {message}" if message else context
    raise GF3DError(code, name_buf.value.decode("utf-8", "replace").strip(), message)


def version() -> str:
    """The library's own version string."""
    buf = ctypes.create_string_buffer(64)
    with LIBRARY_LOCK:
        lib.gf3d_version(buf, len(buf))
    return buf.value.decode("utf-8", "replace")
