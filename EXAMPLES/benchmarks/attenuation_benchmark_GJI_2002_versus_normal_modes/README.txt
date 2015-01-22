
To run this benchmark with attenuation from our 2002 GJI paper (Dimitri Komatitsch and Jeroen Tromp, Spectral-element simulations of global seismic wave propagation, Part I: Validation, Geophysical Journal International, vol. 149, p. 390-412 (2002), Figure 21) do this:

1/ run the current version of SPECFEM3D_GLOBE using the Par_file, STATIONS and CMTSOLUTION files that you will find in this directory

2/ when the run is finished, copy the CI.PAS.*.sem.ascii seismograms to this directory i.e. replace the current files

3/ run ./convolve_all.csh
(edit that script to change the name of the Fortran compiler if GNU gfortran is not installed on your system)

4/ type "gnuplot plot_Bolivia_earthquake_comparison_with_modes_at_station_PAS.gnu"; you should see three figures, each showing three seismograms that should be almost identical; if so, the benchmark is successful; otherwise there is a problem somewhere. Small differences of about 1% to 5% with the normal-mode reference solution are normal (the normal-mode solution with attenuation itself is not exact, it is only "quasi-analytical" and attenuation is handled based on an approximation in the normal-mode technique)

Dimitri Komatitsch, CNRS Marseille, France, September 2013.

