#=====================================================================
#
#          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
#          --------------------------------------------------
#
#          Main authors: Dimitri Komatitsch and Jeroen Tromp
#                        Princeton University, USA
#             and University of Pau / CNRS / INRIA, France
# (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
#                            August 2013
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#=====================================================================

#######################################

specfem3D_TARGETS = \
	$E/xspecfem3D \
	$(EMPTY_MACRO)

specfem3D_OBJECTS = \
	$O/assemble_MPI_scalar.solver.o \
	$O/assemble_MPI_vector.solver.o \
	$O/comp_source_spectrum.solver.o \
	$O/comp_source_time_function.solver.o \
	$O/compute_adj_source_frechet.solver.o \
	$O/compute_arrays_source.solver.o \
	$O/convert_time.solver.o \
	$O/define_derivation_matrices.solver.o \
	$O/get_backazimuth.solver.o \
	$O/get_cmt.solver.o \
	$O/get_event_info.solver.o \
	$O/multiply_arrays_source.solver.o \
	$O/netlib_specfun_erf.solver.o \
	$O/recompute_jacobian.solver.o \
	$O/write_output_ASCII.solver.o \
	$O/write_output_SAC.solver.o \
	$(EMPTY_MACRO)

# solver objects with statically allocated arrays; dependent upon
# values_from_mesher.h
specfem3D_OBJECTS += \
	$O/specfem3D_par.solverstatic.o \
	$O/check_stability.solverstatic.o \
	$O/compute_add_sources.solverstatic.o \
	$O/compute_boundary_kernel.solverstatic.o \
	$O/compute_coupling.solverstatic.o \
	$O/compute_element.solverstatic.o \
	$O/compute_element_strain.solverstatic.o \
	$O/compute_forces_acoustic_calling_routine.solverstatic.o \
	$O/compute_forces_viscoelastic_calling_routine.solverstatic.o \
	$O/compute_forces_crust_mantle_noDev.solverstatic.o \
	$O/compute_forces_crust_mantle_Dev.solverstatic.o \
	$O/compute_forces_inner_core_noDev.solverstatic.o \
	$O/compute_forces_inner_core_Dev.solverstatic.o \
	$O/compute_forces_outer_core_noDev.solverstatic.o \
	$O/compute_forces_outer_core_Dev.solverstatic.o \
	$O/compute_kernels.solverstatic.o \
	$O/compute_seismograms.solverstatic.o \
	$O/compute_stacey_crust_mantle.solverstatic.o \
	$O/compute_stacey_outer_core.solverstatic.o \
	$O/finalize_simulation.solverstatic.o \
	$O/get_attenuation.solverstatic.o \
	$O/initialize_simulation.solverstatic.o \
	$O/iterate_time.solverstatic.o \
	$O/locate_receivers.solverstatic.o \
	$O/locate_regular_points.solverstatic.o \
	$O/locate_sources.solver.o \
	$O/noise_tomography.solverstatic.o \
	$O/prepare_timerun.solverstatic.o \
	$O/read_arrays_solver.solverstatic.o \
	$O/read_forward_arrays.solverstatic.o \
	$O/read_mesh_databases.solverstatic.o \
	$O/read_topography_bathymetry.solverstatic.o \
	$O/save_forward_arrays.solverstatic.o \
	$O/save_kernels.solverstatic.o \
	$O/save_regular_kernels.solverstatic.o \
	$O/setup_GLL_points.solverstatic.o \
	$O/setup_sources_receivers.solverstatic.o \
	$O/specfem3D.solverstatic.o \
	$O/update_displacement_Newmark.solverstatic.o \
	$O/write_movie_output.solverstatic.o \
	$O/write_movie_volume.solverstatic.o \
	$O/write_movie_surface.solverstatic.o \
	$O/write_seismograms.solver.o \
	$(EMPTY_MACRO)

# These files come from the shared directory
specfem3D_SHARED_OBJECTS = \
	$O/shared_par.shared.o \
	$O/auto_ner.shared.o \
	$O/broadcast_computed_parameters.shared.o \
	$O/calendar.shared.o \
	$O/count_elements.shared.o \
	$O/count_number_of_sources.shared.o \
	$O/count_points.shared.o \
	$O/create_name_database.shared.o \
	$O/define_all_layers.shared.o \
	$O/euler_angles.shared.o \
	$O/exit_mpi.shared.o \
	$O/force_ftz.cc.o \
	$O/get_model_parameters.shared.o \
	$O/get_timestep_and_layers.shared.o \
	$O/get_value_parameters.shared.o \
	$O/gll_library.shared.o \
	$O/hex_nodes.shared.o \
	$O/intgrl.shared.o \
	$O/lagrange_poly.shared.o \
	$O/make_ellipticity.shared.o \
	$O/make_gravity.shared.o \
	$O/model_prem.shared.o \
	$O/model_topo_bathy.shared.o \
	$O/parallel.sharedmpi.o \
	$O/param_reader.cc.o \
	$O/read_compute_parameters.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/spline_routines.shared.o \
	$O/write_c_binary.cc.o \
	$(EMPTY_MACRO)

###
### CUDA
###

cuda_OBJECTS = \
  $O/assemble_MPI_scalar_cuda.cuda.o \
	$O/assemble_MPI_vector_cuda.cuda.o \
  $O/check_fields_cuda.cuda.o \
	$O/compute_add_sources_elastic_cuda.cuda.o \
	$O/compute_coupling_cuda.cuda.o \
	$O/compute_forces_crust_mantle_cuda.cuda.o \
	$O/compute_forces_inner_core_cuda.cuda.o \
	$O/compute_forces_outer_core_cuda.cuda.o \
	$O/compute_kernels_cuda.cuda.o \
	$O/compute_stacey_acoustic_cuda.cuda.o \
	$O/compute_stacey_elastic_cuda.cuda.o \
	$O/initialize_cuda.cuda.o \
	$O/it_update_displacement_cuda.cuda.o \
	$O/noise_tomography_cuda.cuda.o \
	$O/prepare_mesh_constants_cuda.cuda.o \
	$O/transfer_fields_cuda.cuda.o \
	$O/write_seismograms_cuda.cuda.o \
	$O/save_and_compare_cpu_vs_gpu.cudacc.o \
	$(EMPTY_MACRO)

cuda_STUBS = \
	$O/specfem3D_gpu_cuda_method_stubs.cudacc.o \
	$(EMPTY_MACRO)

cuda_DEVICE_OBJ = \
	$O/cuda_device_obj.o \
	$(EMPTY_MACRO)

LIBSPECFEM_SOLVER = $O/libspecfem_solver.a

ifeq ($(CUDA),yes)
XSPECFEM_OBJECTS = $(specfem3D_SHARED_OBJECTS) $(specfem3D_OBJECTS) $(cuda_OBJECTS)
else
XSPECFEM_OBJECTS = $(specfem3D_SHARED_OBJECTS) $(specfem3D_OBJECTS) $(cuda_STUBS)
endif

###
### ADIOS
###

adios_OBJECTS = \
	$O/adios_helpers.shared.o \
	$O/adios_manager.shared.o \
	$O/read_arrays_solver_adios.solverstatic.o \
	$O/read_attenuation_adios.solverstatic.o \
	$O/read_forward_arrays_adios.solverstatic.o \
	$O/read_mesh_databases_adios.solverstatic.o \
	$O/save_forward_arrays_adios.solverstatic.o \
	$O/write_specfem_adios_header.solverstatic.o \
	$(EMPTY_MACRO)

adios_STUBS = \
	$O/adios_method_stubs.shared.o \
	$(EMPTY_MACRO)

# conditional adios linking
ifeq ($(ADIOS),yes)
XSPECFEM_OBJECTS += $(adios_OBJECTS)
else
XSPECFEM_OBJECTS += $(adios_STUBS)
endif

###
### VTK
###

vtk_OBJECTS = \
  $O/visual_vtk.visualcc.o \
	$(EMPTY_MACRO)

vtk_STUBS = \
	$O/visual_vtk_stubs.visualc.o \
	$(EMPTY_MACRO)

# conditional adios linking
ifeq ($(VTK),yes)
XSPECFEM_OBJECTS += $(vtk_OBJECTS)
else
XSPECFEM_OBJECTS += $(vtk_STUBS)
endif


#######################################

####
#### rules for executables
####

ifeq ($(CUDA),yes)
## cuda version

ifeq ($(CUDA5),yes)

## cuda 5 version
${E}/xspecfem3D: $(XSPECFEM_OBJECTS)
	@echo "building xspecfem3D with CUDA 5 support"
	${NVCCLINK} -o $(cuda_DEVICE_OBJ) $(cuda_OBJECTS)
	${FCLINK} -o ${E}/xspecfem3D $(XSPECFEM_OBJECTS) $(cuda_DEVICE_OBJ) $(MPILIBS) $(CUDA_LINK)
else

## cuda 4 version
${E}/xspecfem3D: $(XSPECFEM_OBJECTS)
	@echo "building xspecfem3D with CUDA 4 support"
	${FCLINK} -o ${E}/xspecfem3D $(XSPECFEM_OBJECTS) $(MPILIBS) $(CUDA_LINK)
endif

else

## non-cuda version
${E}/xspecfem3D: $(XSPECFEM_OBJECTS)
	@echo "building xspecfem3D without CUDA support"
## use MPI here
## DK DK add OpenMP compiler flag here if needed
#	${MPIFCCOMPILE_CHECK} -qsmp=omp -o ${E}/xspecfem3D $(XSPECFEM_OBJECTS) $(MPILIBS)
	${MPIFCCOMPILE_CHECK} -o ${E}/xspecfem3D $(XSPECFEM_OBJECTS) $(MPILIBS)

endif

#######################################

## compilation directories
S := ${S_TOP}/src/specfem3D
$(specfem3D_OBJECTS): S = ${S_TOP}/src/specfem3D

####
#### rule for each .o file below
####

###
### specfem3D - optimized flags and dependence on values from mesher here
###
$O/%.solverstatic.o: $S/%.f90 ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.solverstatic.o: $S/%.F90 ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.solverstatic_openmp.o: $S/%.f90 ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h
## DK DK add OpenMP compiler flag here if needed
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -qsmp=omp -o $@ $<

###
### no dependence on values from mesher here
###

$O/%.solver.o: $S/%.f90 ${SETUP}/constants.h
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.solver.o: $S/%.F90 ${SETUP}/constants.h
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<


###
### VTK compilation
###

$O/%.visualcc.o: $S/%.cpp ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(MPI_INC) -o $@ $<

$O/%.visualc.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(MPI_INC) -o $@ $<

