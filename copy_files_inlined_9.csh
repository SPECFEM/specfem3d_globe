#!/bin/csh

# copy files for the inlined version of the solver

cp -f compute_forces_crust_mantle_inlined_9.f90 compute_forces_crust_mantle.f90
cp -f compute_forces_inner_core_inlined_9.f90 compute_forces_inner_core.f90
cp -f compute_forces_outer_core_inlined_9.f90 compute_forces_outer_core.f90
cp -f specfem3D_inlined.f90 specfem3D.f90

