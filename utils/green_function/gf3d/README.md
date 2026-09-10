# `gf3d` — Green functions in memory

A ctypes binding to `lib/libgf3d.so`, the Green function extraction library
of `src/gf3d/`. It runs the same code `bin/xgf3d` runs; what it adds is
getting seismograms and their partial derivatives into numpy arrays without
a process launch and a SAC file per iteration, which is what a source
inversion needs.

## Building the library

The package is pure Python, but the library it binds is not:

```bash
./configure --with-hdf5 HDF5_INC=<dir> HDF5_LIBS=-L<dir>
make gf3d
```

That writes `lib/libgf3d.so` (and `lib/libgf3d.a`, `include/gf3d.h`,
`include/gf3d.mod` for C and Fortran callers). The HDF5 library directories
are baked into the shared object as an rpath, so nothing needs
`LD_LIBRARY_PATH` — except the compiler's own runtime, if the tree was
built with a module-provided compiler.

## Using it

```bash
pip install -e utils/green_function          # or set PYTHONPATH
```

```python
import numpy as np
import gf3d

with gf3d.Database("EXAMPLES/green_function_database/regional/GFDB") as db:
    cmt = gf3d.CMTSource.read(".../validation_data/CMTSOLUTION")
    r = db.partials(cmt)

    r.data          # (nstations, 3, nt), metres, components N, E, Z
    r.t             # (nt,), seconds, t = 0 at the centroid time
    r.dp            # (nstations, 10, 3, nt)
    r.dp_names      # ['Mrr','Mtt','Mpp','Mrt','Mrp','Mtp','lat','lon','dep','tim']
    r.dp_units      # ['m/dyne-cm', ..., 'm/deg', 'm/deg', 'm/km', 'm/s']
    r.station_ids   # ['IU.LVC', 'IU.SDV', 'IU.SJG']
```

The partials are analytic and exactly linear in the moment tensor, so

```python
(r.dp[:, :6] * np.asarray(cmt.tensor)[None, :, None, None]).sum(1)
```

reproduces `r.data` to round-off with the CMTSOLUTION's own numbers, in
dyne·cm. That identity is the cheapest check that an inversion is using the
partials the way the library means them.

`db.seismograms(source)` skips the partials. `db.plan(source)` answers what
the output axis will be without extracting anything. `db.locate(lat, lon,
depth_km)` says which element a point falls in. A `ForceSource` works
wherever a `CMTSource` does, except that it has no partials.

`Result.to_stream()` returns an obspy `Stream` if obspy is installed
(`pip install -e 'utils/green_function[obspy]'`); nothing else in the
package needs it.

## Things worth knowing

**`t = 0` is the centroid time**, not the origin time. The CMTSOLUTION's
`time shift` is carried in `cmt.time_shift` and applied only when an
absolute start time is needed, as `to_stream()` does. Putting it into the
trace instead is a mistake that looks right in isolation and is 29 seconds
wrong against real data for the shipped example.

**Errors are exceptions, never a crash.** Every failure — a database that
is not there, a source outside it, a NaN, a closed handle — raises
`GF3DError` with `.code`, `.name` and `.message`. `err.code ==
gf3d.GF_ERR_NO_ELEMENT` is the one an inversion should expect: a trial
source that wandered outside the database.

**The library is not thread-safe** and this package serialises calls on a
module-level lock. The process-wide state is the last error message, the
search tree, and specfem's own parameter module. Extraction is a C call, so
other Python threads keep running; they just wait at the library's door.

**Two databases can be open at once**, but the search tree is rebuilt each
time an extraction switches between them, so alternating is slow. Two
databases of *different planets or topography grids* must not be used
alternately at all: specfem holds one set of grid dimensions per process.

**The onset warning** fires when a trace already carries more than about
1e-3 of its peak at the first sample. It means the source time function
conversion is reaching back past the start of the stored record and the
first arrivals may be contaminated — a longer reciprocal run, or a later
`t0`, is the cure. It is a warning and not an error because whether it
matters depends on the arrival being measured.

**If two `gf3d` packages collide.** The name is shared with the older
standalone GF3D Python package, which reimplemented extraction over h5py.
This one supersedes it and the two cannot live in one environment.

## Relation to the rest of the tree

| | |
|---|---|
| `include/gf3d.h` | the C ABI this binds, and the document of record for units, shapes and index order |
| `src/gf3d/gf3d_capi.F90` | the `bind(C)` façade behind it |
| `src/gf3d/gf3d.F90` | the same library for a Fortran caller (`use gf3d`) |
| `bin/xgf3d` | the command-line tool, for SAC and ASCII output |
| `tests/gf3d/9f.test_gf_python.sh` | checks this package against `xgf3d`'s own output |
