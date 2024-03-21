README
======


This folder contains several examples of how to run forward and adjoint simulations.
Each example is explained in more detail in a readme file contained in the example folders.


The following examples are provided:

- basic usage:

  * `regional_Greece_small/`

  contains an example for a small regional simulation and an event located in southern Greece;
  the example can be run as a small test on a single desktop machine
  (4 CPUs, forward simulation lasts ~5min, kernel simulation lasts ~10min)

  * `regional_MiddleEast/`

  contains an example for a regional simulation for the Middle East region
  (64 CPUs, forward simulation lasts ~7min, kernel simulation lasts ~15min)

  * `global_s362ani_shakemovie/`

  contains an example for a global shakemovie simulation using S362ANI as 3-D background model
  (384 CPUs, forward simulation lasts ~1h 55min)

  * `global_PREM_kernels/`

  contains examples for amplitude and traveltime kernels using PREM as background model;
  both examples focus on how adjoint sources for filtered measurements would be constructed and used.
  (24 CPUs, forward simulations last ~5min, kernel simulations last ~10min)


- PREM benchmark solutions:

  * `benchmarks/`

  contains PREM benchmarks for comparing SEM outputs with normal-mode solutions
  (various resolutions: Vanuatu 384 CPUs/ Bolivia 1536 CPUs)


- seismic interferometry:

  * `noise_examples/`

  contains examples of how to compute noise sensitivity kernels for global and regional simulations
  (various resolutions: global 600 CPUs/ regional 400 CPUs)



