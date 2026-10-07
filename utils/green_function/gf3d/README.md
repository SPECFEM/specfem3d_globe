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

Two configure options make it faster, neither of them a default:

- `--enable-openmp`: a handle that holds its element (`max_elements > 0`)
  computes the stations over `OMP_NUM_THREADS` threads, every number the
  same as on one. Set `OMP_NUM_THREADS` when several processes share a
  node; with it unset, each takes every core.
- `FCFLAGS="-march=x86-64-v3"` (or the machine's own ISA): about a quarter
  faster again under gfortran. It applies to the whole build, and fused
  multiply-adds move the last digits (~1e-13 of a trace). Intel builds
  already get `-xHost` by default.

## Using it

```bash
pip install -e utils/green_function          # or set PYTHONPATH
```

Or from the lock file: `uv sync --project utils/green_function`, with
`--group test` to add what `tests/gf3d` needs (obspy, scipy, h5py).

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

### Doing the contraction yourself

A seismogram is one dot product per sample between an element's stored
displacement and 375 weights that depend only on the source (the
`w_{p,ijk}` of the user manual's Green function chapter). A sampler that evaluates many
sources in one element can take the two ingredients and do the rest itself,
on a GPU for instance:

```python
W = db.weights(cmt, kind=2)              # locates; W.w (375,), W.dw (9, 375), W.scale (nstations,)
u = db.element_block(W.location.ielem)   # float32 (nstations, 3, nt_db, 375), read once per element
x = W.scale[:, None, None] * (u @ W.w)   # (nstations, 3, nt_db), before the source time function

p = db.plan(cmt)
k = gf3d.stf_kernel(p.kind_stf, p.hdur_corr, p.dt_sub, p.trunc)
```

`x` converted with `k` as `stf_kernel`'s docstring writes out is
`db.seismograms(cmt).data` to rounding; `u @ W.dw[c]` likewise gives the
first nine columns of `db.partials(cmt).dp` (the centroid time has no
weight: it shifts the trace). The source time function is the caller's:
`db.plan` gives the library's own kernel parameters, and any others are
allowed.

- `element_block(ielem, stations=..., nt=...)` reads a station subset and
  only the first `nt` samples; the last `plan.khalf` converted samples need
  data past the prefix, so read that many more than you keep. For 142
  stations and 917 samples that is 0.59 GB, ~0.4 s on one core. `out=`
  takes a preallocated (e.g. pinned) float32 buffer.
- It reads from disk on every call and never touches the element cache:
  `cache_stats` counts its files in `files_read` and nothing else. Keep the
  blocks yourself, keyed by `W.location.ielem`.
- `weights` reads no displacement: ~0.08 ms per call once the handle has
  the element's coordinates. A source that moves into another element
  changes `W.location.ielem`; the forward model jumps there, as the mesh's
  strain does.

`EXAMPLES/green_function_database/global/api_demos/python/gf3d_sampler_demo.py`
does all of this on the global example and checks it against `seismograms`
and `partials`.

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

**A sampler should open with `max_elements`.** An extraction reads its
element's displacement for every station from disk, and on a large
database that read is most of the cost: 185 stations × 3725 samples is
3 GB, about 2.2 s per call on one core (4.5 s with all ten partials), of
which the read is 1.8 s. Held in memory, the same calls take 0.5 s and
2.5 s, or 0.02 s and 0.11 s over 32 threads of a library built with
OpenMP. A
handle opened as `gf3d.Database(path, max_elements=N)` keeps the `N`
elements it used most recently in memory and drops the least recently used
one to make room; a position inside a kept element reads no displacement.
Every handle also keeps the coordinates its locates read, 3 kB per element,
for the `max(10, N)` elements used most recently, so returning to a
position still held reads nothing from disk. The numbers are the same to
the bit either way. Each element costs `db.info["bytes_per_element"]`;
`db.cache_stats` reports hits, misses, evictions, the elements held and the
files read. The default, `0`, keeps nothing.

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
