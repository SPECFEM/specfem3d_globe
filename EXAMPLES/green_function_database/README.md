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
