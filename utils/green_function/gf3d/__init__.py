"""Green function extraction from a specfem3d_globe reciprocal database.

A thin ctypes binding over ``lib/libgf3d.so``, which is the same code
``bin/xgf3d`` runs. Nothing here reimplements any part of the extraction:
the point of the package is to get seismograms and their partial
derivatives into memory, for an inversion loop that cannot afford a process
launch and a SAC file per iteration.

    import gf3d

    with gf3d.Database("EXAMPLES/green_function_database/regional/GFDB") as db:
        cmt = gf3d.CMTSource.read(".../validation_data/CMTSOLUTION")
        r = db.partials(cmt)

        r.data          # (nstations, 3, nt) metres, components N, E, Z
        r.t             # (nt,) seconds, t = 0 at the centroid time
        r.dp            # (nstations, 10, 3, nt) partial derivatives
        r.dp_names      # ['Mrr', ..., 'Mtp', 'lat', 'lon', 'dep', 'tim']
        r.station_ids   # ['IU.LVC', 'IU.SDV', 'IU.SJG']

The partials are analytic and linear, so the moment-tensor half satisfies

    (r.dp[:, :6] * cmt.tensor[None, :, None, None]).sum(1)  ==  r.data

to round-off, with the CMTSOLUTION's own numbers, in dyne-cm.

Finding the library
-------------------
``$GF3D_LIB`` if set, else ``<repo>/lib/libgf3d.so`` relative to this file,
else the system search path. Build it with ``./configure --with-hdf5`` and
``make gf3d``.

Threads
-------
The library keeps process-wide state and is not thread-safe; every call
here is serialised on a module-level lock. Extraction is a C call, so other
Python threads keep running while it is inside the library.
"""

from ._lib import (
    GF_NCOMP,
    GF_NDP_LOC,
    GF_NDP_MT,
    GF_OK,
    GF_ERR_ARG,
    GF_ERR_GEOMETRY,
    GF_ERR_INCOMPLETE,
    GF_ERR_NO_ELEMENT,
    GF_ERR_NO_PATH,
    GF3DError,
    library_path,
    version as library_version,
)
from .database import Database, Location, Plan, Result, Station, open
from .sources import CMTSource, ForceSource

__all__ = [
    "Database",
    "open",
    "CMTSource",
    "ForceSource",
    "Result",
    "Plan",
    "Location",
    "Station",
    "GF3DError",
    "library_version",
    "library_path",
    "GF_NCOMP",
    "GF_NDP_MT",
    "GF_NDP_LOC",
    "GF_OK",
    "GF_ERR_ARG",
    "GF_ERR_NO_PATH",
    "GF_ERR_NO_ELEMENT",
    "GF_ERR_GEOMETRY",
    "GF_ERR_INCOMPLETE",
    "__version__",
]

__version__ = "0.1.0"
