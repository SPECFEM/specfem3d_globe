# Green Function Database Example

This example demonstrates how to build a Green function (GF) database using
SPECFEM3D_GLOBE and validate it against direct forward simulations.

## Prerequisites

- SPECFEM3D_GLOBE compiled with HDF5 and GF database support
  (`GF_DATABASE_ENABLED = .true.` in `Par_file`)
- MPI (e.g., OpenMPI or MPICH)
- [uv](https://docs.astral.sh/uv/) for Python dependency management

## Quick Start

1. **Install Python dependencies**

   From this directory (`EXAMPLES/green_function_database/`):

   ```bash
   uv sync
   ```

2. **Run the regional workflow**

   ```bash
   cd regional
   snakemake -j1
   ```

   This will:
   - Set up the base directory and run the mesher
   - Run reciprocal simulations (N, E, Z force components) for each station
   - Build the GF database manifest (`GFDB/centroids.bin`)
   - Run forward validation simulations (force and CMT)
   - Extract seismograms from the database with `xgf3d` for both validation
     sources and compare them with the forward runs: per-station figures,
     a record section per source type, a summary figure and table, and the
     pass/fail verdict in `validation_output/gf_compare_gate.json`. The
     workflow fails if any station's misfit exceeds `GF_COMPARE_THRESHOLD`;
     the figures survive a failure, on purpose, so it can be diagnosed.

   Extraction by hand, once the database is built:

   ```bash
   bin/xgf3d --seis <GFDB> <CMTSOLUTION|FORCESOLUTION> <outdir> [--format sac|sacan|ascii|all] \
                    [--partials 1|2] [--t0 <s>]
   ```

   writes `NET.STA.BX{N,E,Z}.sem.sac` with the solver's own header rules
   (`--format sac`, the default), the alphanumeric form (`sacan`), or the
   `NET.STA.gf3d.txt` columns the comparison reads (`ascii`; the workflow
   asks for this one explicitly). `--partials 1` adds the six moment-tensor
   partial derivatives of a CMTSOLUTION's seismograms in the same format
   (`NET.STA.BXN.Mrr.sem.sac`, or `NET.STA.partials.txt`), per dyne-cm;
   `--partials 2` adds latitude, longitude, depth and centroid time.

   **Note 1**: The regional workflow is affected by the absorbing boundary conditions. It is important to choose the stations carefully to avoid strong reflections from the boundaries.

3. **Run the global workflow**

   ```bash
   cd global
   snakemake -j1
   ```

   **Note 1**: The global workflow is computationally more expensive but does not have issues with boundary reflections. It is recommended to run the global workflow if you have sufficient computational resources.

   **Note 2 — runtime**: on a node with two H100 GPUs, the complete global
   workflow — mesher, 12 reciprocal simulations (4 stations × 3 components),
   the two forward validation simulations and the comparison — took about
   **40 minutes** with `snakemake -j1`, i.e. one simulation at a time at
   roughly 3 minutes each. Budget for that before starting. A failed or
   interrupted run resumes from its marker files (`mesher.flag`,
   `simulations/<station>/setup.flag`, `done_{N,E,Z}.flag`), so nothing
   already finished is repeated.


## Calling the library yourself

Everything above goes through `xgf3d` and writes files. For an inversion,
what you want instead is the seismograms and their partial derivatives in
memory. The library is callable directly from Python and from Fortran, and
`global/` holds a pair of runnable demonstrations of each — four short
programs against the database the global workflow builds. They are not part
of the workflow: nothing in either `Snakefile` refers to those directories,
so run them by hand once `global/GFDB` exists.

All four need `lib/libgf3d.so`, which comes from `make gf3d` in a tree
configured `--with-hdf5`.

### Python

The `gf3d` package is a ctypes binding, so it needs no compiler:

```python
import gf3d

with gf3d.Database("global/GFDB") as db:
    cmt = gf3d.CMTSource.read("global/validation_data/CMTSOLUTION")
    r = db.partials(cmt)

    r.data        # (nstations, 3, nt) metres, components N, E, Z
    r.t           # (nt,) seconds, t = 0 at the centroid time
    r.dp          # (nstations, 10, 3, nt) partial derivatives
    r.dp_names    # ['Mrr', ..., 'Mtp', 'lat', 'lon', 'dep', 'tim']
```

`global/api_demos/python/` holds two scripts. Each writes one figure and nothing
else, and reports how long every call into the library took.

```bash
cd global/api_demos/python
../../../.venv/bin/python gf3d_api_demo.py       # the partials, and what they are for
../../../.venv/bin/python gf3d_events_demo.py    # several events from one open database
```

`gf3d_api_demo.py` opens the database, extracts seismograms and all ten
partial derivatives, checks the identity the moment-tensor partials satisfy
exactly, and then uses a centroid partial as a derivative: it predicts the
seismogram of a relocated source from a single Taylor term and shows the
error falling fourfold each time the relocation is halved, which is what a
correct first derivative must do and a plausible-looking wrong one will
not.

`gf3d_events_demo.py` is the catalogue case. The database was not built for
one earthquake — `db_base/DATA/GF_LOCATIONS` declares three hypocentres,
and the reciprocal runs stored the strain around all of them — so the
script opens the database once and extracts for every hypocentre it
declares. The timings are the point: tens of milliseconds per event against
the forty minutes of simulation that built the database.

See `utils/green_function/gf3d/README.md` for the package itself.

### Fortran

A Fortran program needs `use gf3d` and nothing else. `global/api_demos/fortran/`
holds the same two demonstrations, without the figures:

```bash
cd global/api_demos/fortran
make                     # or: make FC=ifort
./gf3d_api_demo
./gf3d_events_demo
```

They print the same numbers as their Python counterparts, because it is the
same library underneath.

The link line is worth reading, in `global/api_demos/fortran/Makefile`:

```
$(FC) demo.f90 -I<repo>/include -L<repo>/lib -lgf3d -Wl,-rpath,<repo>/lib
```

and that is all of it. No HDF5 flags and no list of the shared objects the
library reuses: `lib/libgf3d.so` records HDF5 as a dependency of its own and
the `-rpath` lets the loader find it, so a caller does not have to know
where this tree's HDF5 came from. Linking the static `lib/libgf3d.a` works
too and is what `tests/gf3d/` does, but then the HDF5 link line becomes the
caller's problem.

One thing to know: Fortran module files are compiler-specific, so
`include/gf3d.mod` can only be read by the compiler that wrote it — the one
the tree was configured with. The Makefile defaults to `gfortran`; if that
is not what built the library you will get `Cannot open module file
'gf3d.mod'`, and `make FC=<yours>` is the fix.


## Workflow Configuration

The Snakefile accepts configuration overrides via `--config`. Key options:

| Option              | Default       | Description                              |
|---------------------|---------------|------------------------------------------|
| `SPECFEM_DIR`       | `../../..`    | Path to the specfem3d_globe root         |
| `NPROC`             | `4`           | Number of MPI ranks                      |
| `MPIRUN`            | `mpirun`      | MPI launcher command                     |
| `CREATE_VALIDATION` | `True`        | Run validation forward simulations       |
| `GF_COMPARE_THRESHOLD` | `5e-3`     | Relative L2 misfit (station vector) above which the comparison fails |

Example with overrides:

```bash
snakemake -j1 --config NPROC=6 MPIRUN="srun"
```

## Parallelism

Stations can be run in parallel (components within a station are always
sequential). Use `-jN` with the `mpi` resource to control this:

```bash
snakemake -j4 --resources mpi=2
```

This allows up to 4 tasks in parallel, but limits MPI simulations to 2
concurrent runs.

## Cleaning Up

```bash
cd regional
snakemake clean_all    # remove all generated files
```

Or clean specific parts:

```bash
snakemake clean_base          # mesher output and symlinks
snakemake clean_simulations   # station simulation directories
snakemake clean_database      # GF database files
snakemake clean_validation    # force validation directory
snakemake clean_validation_cmt # CMT validation directory
```

## Rebuilding

Changing `GF_SUBSAMPLE_STEP` in `db_base/DATA/Par_file` or the station list in
`db_base/DATA/STATIONS` needs `snakemake clean_database clean_simulations`
first: the database's `mesh_info.h5` is created once and never overwritten
(it would otherwise disagree with the new traces), and a station's
FORCESOLUTION is written only when its `setup.flag` is absent. The mesher
output and the two forward validation runs do not depend on either and can
stay.

After editing the Snakefile, Snakemake 9 may want to redo jobs whose rule
code changed even though their outputs exist; `snakemake --rerun-triggers
mtime` runs only what is missing (put any target *before* that option).
