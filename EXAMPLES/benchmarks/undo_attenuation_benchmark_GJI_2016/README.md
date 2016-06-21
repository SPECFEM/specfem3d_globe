Undo Attenuation Benchmark
==========================

This benchmark reproduces the results of section 3 of the paper "Anelastic
sensitivity kernels with parsimonious storage for adjoint tomography and full
waveform inversion" by D. Komatitsch, Z. Xie, E. Bozdag, E. Sales de Andrade,
D. B. Peter, Q. Liu, and J. Tromp, published in Geophys. J. Int. 2016,
doi: 10.1093/gji/ggw224. The exact revision of the code used to produce the
figure in the paper is 74ffd4e330f281d7964855d1845d1b72bf13ab74.

There are three phases to the benchmark:
  * Exact computation of adjoint kernel using a full dump of the attenuated
    forward field. As this method requires the full field at all time steps,
    the storage requirements are quite large (~16 TiB) and it is only
    implemented for the alpha kernel. On the other hand, since attenuation does
    not need to be reversed, the result is exact.
  * Calculation of the adjoint kernel using a forward field with partial
    physical dispersion only. Since true attenuation is not used, the result is
    not exact, but it avoids the instability with solving the backward
    equations.
  * Calculation of the adjoint kernel using the newly proposed undo attenuation
    method. Via checkpointing, this method achieves a balance between excessive
    storage space and excessive computation, and is mathematically exact.

Building
--------

Before building, you may wish to modify the `DATA/Par_file` to match the
available processors and memory on your cluster.

Running the `./build.sh` in this directory will configure and build three
sub-directories for each of the stages. Running `./build.sh -h` will print out
any options that are allowed. You may also run `<path/to/>build.sh` from some
path outside of this directory to create the build directories elsewhere.

Each phase's sub-directory contains two sub-directories, `step1-forward` and
`step2-adjoint`, which are used for the forward and adjoint simulations,
respectively.

Running the Forward Simulation
------------------------------

The forward simulations for the three phases are completely independent and may
be run simultaneously. They require the standard run of `xmeshfem3D` and
`xspecfem3D`. You may use whatever normal method of running these programs you
would usually use, e.g., a script for a job manager on your cluster.

Preparing the Adjoint Simulations
---------------------------------

After the forward simulation is complete, run the `mkadjsrc.py` script to
create the adjoint sources and the individual directories for the runs. This
script creates two adjoint simulations, one with the source filtered between
50-100 seconds, and the other filtered between 100-200 seconds. Only the latter
is shown in the referenced paper, but you may wish to run the former to see the
effect in that case. For example, running:

	$ ./mkadjsrc.py exact/step1-forward exact/step2-adjoint

will create the directories `exact/step2-adjoint/50-100s` and
`exact/step2-adjoint/100-200s` for the two filter bands. Run `./mkadjsrc.py -h`
to see additional options for plotting the sources.

NOTE: In order to conserve disk space, the `DATABASES_MPI` and `huge_dumps`
directories are symbolic links to the forward simulation directories. However,
since the kernel results are written to the `DATABASES_MPI` directory, you
cannot run both filter bands at the same time. After running one adjoint
simulation, you can move the resulting kernel files aside before running the
second adjoint simulation. Alternatively, you can change the linked directories
into full copies, at the expense of a few TiB of disk space.

Running the Adjoint Simulations
-------------------------------

The adjoint simulations for each phase are again independent, but as noted
above, the two filtering bands within each phase are not. Running the adjoint
simulations should proceed the same as the forward simulation (e.g., using
whatever cluster you have available.)

The resulting kernel files will be in
`<phase>/step2-adjoint/<band>/DATABASES_MPI/proc??????_alpha_kernel.bin`
where `<phase>` is `exact`, `ppd`, or `undo`; `<band>` is `50-100s` or
`100-200s`. You can visualize these files in whatever manner you wish, such as
the scripts in the `utils/Visualization/VTK_ParaView/` directory of this
repository. The figures created for the paper can be reproduced using the code
in the [undo attenuation benchmark
repository](https://gitlab.com/QuLogic/undo_attenuation_benchmark).
