#=====================================================================
#
#          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
#          --------------------------------------------------
#
#          Main authors: Dimitri Komatitsch and Jeroen Tromp
#                        Princeton University, USA
#             and University of Pau / CNRS / INRIA, France
# (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
#                            April 2011
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
	$O/assemble_MPI_central_cube.o \
	$O/assemble_MPI_scalar.o \
	$O/assemble_MPI_vector.o \
	$O/assemble_MPI_central_cube_block.o \
	$O/assemble_MPI_scalar_block.o \
	$O/assemble_MPI_vector_block.o \
	$O/comp_source_spectrum.o \
	$O/comp_source_time_function.o \
	$O/compute_adj_source_frechet.o \
	$O/compute_arrays_source.o \
	$O/convert_time.o \
	$O/create_central_cube_buffers.o \
	$O/define_derivation_matrices.o \
	$O/get_attenuation.o \
	$O/get_backazimuth.o \
	$O/get_cmt.o \
	$O/get_event_info.o \
	$O/locate_receivers.o \
	$O/locate_regular_points.o \
	$O/locate_sources.o \
	$O/multiply_arrays_source.o \
	$O/netlib_specfun_erf.o \
	$O/read_arrays_buffers_solver.o \
	$O/recompute_jacobian.o \
	$O/save_regular_kernels.o \
	$O/write_seismograms.o \
	$O/write_output_ASCII.o \
	$O/write_output_SAC.o \
	$(EMPTY_MACRO)

# solver objects with statically allocated arrays; dependent upon
# values_from_mesher.h
specfem3D_OBJECTS += \
	$O/check_simulation_stability.o \
	$O/compute_add_sources.o \
	$O/compute_boundary_kernel.o \
	$O/compute_coupling.o \
	$O/compute_element.o \
	$O/compute_forces_crust_mantle.o \
	$O/compute_forces_crust_mantle_Dev.o \
	$O/compute_forces_inner_core.o \
	$O/compute_forces_inner_core_Dev.o \
	$O/compute_forces_outer_core.o \
	$O/compute_forces_outer_core_Dev.o \
	$O/compute_kernels.o \
	$O/compute_seismograms.o \
	$O/compute_stacey_crust_mantle.o \
	$O/compute_stacey_outer_core.o \
	$O/fix_non_blocking_flags.o \
	$O/initialize_simulation.o \
	$O/noise_tomography.o \
	$O/prepare_timerun.o \
	$O/read_arrays_solver.o \
	$O/read_forward_arrays.o \
	$O/read_mesh_databases.o \
	$O/save_forward_arrays.o \
	$O/save_kernels.o \
	$O/setup_sources_receivers.o \
	$O/specfem3D.o \
	$O/write_movie_volume.o \
	$O/write_movie_surface.o \
	$(EMPTY_MACRO)

# These files come from the shared directory
specfem3D_SHARED_OBJECTS = \
	$O/auto_ner.o \
	$O/broadcast_compute_parameters.o \
	$O/calendar.o \
	$O/count_number_of_sources.o \
	$O/create_name_database.o \
	$O/euler_angles.o \
	$O/exit_mpi.o \
	$O/force_ftz.o \
	$O/get_model_parameters.o \
	$O/get_value_parameters.o \
	$O/gll_library.o \
	$O/hex_nodes.o \
	$O/intgrl.o \
	$O/lagrange_poly.o \
	$O/make_ellipticity.o \
	$O/make_gravity.o \
	$O/model_prem.o \
	$O/model_topo_bathy.o \
	$O/param_reader.o \
	$O/read_compute_parameters.o \
	$O/read_parameter_file.o \
	$O/read_value_parameters.o \
	$O/reduce.o \
	$O/rthetaphi_xyz.o \
	$O/spline_routines.o \
	$O/write_c_binary.o \
	$(EMPTY_MACRO)


#######################################

####
#### rules for executables
####

${E}/xspecfem3D: $(specfem3D_OBJECTS) $(specfem3D_SHARED_OBJECTS)
## use MPI here
## DK DK add OpenMP compiler flag here if needed
#	${MPIFCCOMPILE_NO_CHECK} -qsmp=omp -o ${E}/xspecfem3D $(specfem3D_OBJECTS) $(specfem3D_SHARED_OBJECTS) $(MPILIBS)
	${MPIFCCOMPILE_NO_CHECK} -o ${E}/xspecfem3D $(specfem3D_OBJECTS) $(specfem3D_SHARED_OBJECTS) $(MPILIBS)

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
$O/compute_add_sources.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/compute_add_sources.f90
	${FCCOMPILE_NO_CHECK} -c -o $O/compute_add_sources.o ${FCFLAGS_f90} $S/compute_add_sources.f90

$O/compute_boundary_kernel.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/compute_boundary_kernel.f90
	${FCCOMPILE_CHECK} -c -o $O/compute_boundary_kernel.o ${FCFLAGS_f90} $S/compute_boundary_kernel.f90

$O/compute_coupling.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/compute_coupling.f90
	${FCCOMPILE_NO_CHECK} -c -o $O/compute_coupling.o ${FCFLAGS_f90} $S/compute_coupling.f90

$O/compute_element.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/compute_element.F90
	${FCCOMPILE_NO_CHECK} -c -o $O/compute_element.o ${FCFLAGS_f90} $S/compute_element.F90

$O/compute_forces_crust_mantle.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/compute_forces_crust_mantle.f90
	${FCCOMPILE_NO_CHECK} -c -o $O/compute_forces_crust_mantle.o ${FCFLAGS_f90} $S/compute_forces_crust_mantle.f90

$O/compute_forces_crust_mantle_Dev.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/compute_forces_crust_mantle_Dev.F90
## DK DK add OpenMP compiler flag here if needed
#	${FCCOMPILE_NO_CHECK} -c -qsmp=omp -o $O/compute_forces_crust_mantle_Dev.o ${FCFLAGS_f90} $S/compute_forces_crust_mantle_Dev.F90
	${FCCOMPILE_NO_CHECK} -c -o $O/compute_forces_crust_mantle_Dev.o ${FCFLAGS_f90} $S/compute_forces_crust_mantle_Dev.F90

$O/compute_forces_outer_core.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/compute_forces_outer_core.f90
	${FCCOMPILE_NO_CHECK} -c -o $O/compute_forces_outer_core.o ${FCFLAGS_f90} $S/compute_forces_outer_core.f90

$O/compute_forces_outer_core_Dev.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/compute_forces_outer_core_Dev.f90
	${FCCOMPILE_NO_CHECK} -c -o $O/compute_forces_outer_core_Dev.o ${FCFLAGS_f90} $S/compute_forces_outer_core_Dev.f90

$O/compute_forces_inner_core.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/compute_forces_inner_core.f90
	${FCCOMPILE_NO_CHECK} -c -o $O/compute_forces_inner_core.o ${FCFLAGS_f90} $S/compute_forces_inner_core.f90

$O/compute_forces_inner_core_Dev.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/compute_forces_inner_core_Dev.F90
	${FCCOMPILE_NO_CHECK} -c -o $O/compute_forces_inner_core_Dev.o ${FCFLAGS_f90} $S/compute_forces_inner_core_Dev.F90

$O/compute_kernels.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/compute_kernels.f90
	${FCCOMPILE_NO_CHECK} -c -o $O/compute_kernels.o ${FCFLAGS_f90} $S/compute_kernels.f90

$O/compute_seismograms.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/compute_seismograms.f90
	${FCCOMPILE_NO_CHECK} -c -o $O/compute_seismograms.o ${FCFLAGS_f90} $S/compute_seismograms.f90

$O/compute_stacey_crust_mantle.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/compute_stacey_crust_mantle.f90
	${FCCOMPILE_NO_CHECK} -c -o $O/compute_stacey_crust_mantle.o ${FCFLAGS_f90} $S/compute_stacey_crust_mantle.f90

$O/compute_stacey_outer_core.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/compute_stacey_outer_core.f90
	${FCCOMPILE_NO_CHECK} -c -o $O/compute_stacey_outer_core.o ${FCFLAGS_f90} $S/compute_stacey_outer_core.f90

$O/read_arrays_solver.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/read_arrays_solver.f90
	${FCCOMPILE_NO_CHECK} -c -o $O/read_arrays_solver.o ${FCFLAGS_f90} $S/read_arrays_solver.f90

$O/read_forward_arrays.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/read_forward_arrays.f90
	${FCCOMPILE_NO_CHECK} -c -o $O/read_forward_arrays.o ${FCFLAGS_f90} $S/read_forward_arrays.f90

$O/save_forward_arrays.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/save_forward_arrays.f90
	${FCCOMPILE_NO_CHECK} -c -o $O/save_forward_arrays.o ${FCFLAGS_f90} $S/save_forward_arrays.f90

$O/save_kernels.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/save_kernels.f90
	${FCCOMPILE_NO_CHECK} -c -o $O/save_kernels.o ${FCFLAGS_f90} $S/save_kernels.f90

$O/save_regular_kernels.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/save_regular_kernels.f90
	${FCCOMPILE_NO_CHECK} -c -o $O/save_regular_kernels.o ${FCFLAGS_f90} $S/save_regular_kernels.f90

###
### specfem3D - regular compilation options here
###
$O/comp_source_time_function.o: $S/comp_source_time_function.f90
	${FCCOMPILE_CHECK} -c -o $O/comp_source_time_function.o ${FCFLAGS_f90} $S/comp_source_time_function.f90

$O/comp_source_spectrum.o: ${SETUP}/constants.h $S/comp_source_spectrum.f90
	${FCCOMPILE_CHECK} -c -o $O/comp_source_spectrum.o ${FCFLAGS_f90} $S/comp_source_spectrum.f90

$O/compute_adj_source_frechet.o: ${SETUP}/constants.h $S/compute_adj_source_frechet.f90
	${FCCOMPILE_CHECK} -c -o $O/compute_adj_source_frechet.o ${FCFLAGS_f90} $S/compute_adj_source_frechet.f90

$O/compute_arrays_source.o: ${SETUP}/constants.h $S/compute_arrays_source.f90
	${FCCOMPILE_CHECK} -c -o $O/compute_arrays_source.o ${FCFLAGS_f90} $S/compute_arrays_source.f90

$O/convert_time.o: $S/convert_time.f90
	${FCCOMPILE_CHECK} -c -o $O/convert_time.o ${FCFLAGS_f90} $S/convert_time.f90

$O/define_derivation_matrices.o: ${SETUP}/constants.h $S/define_derivation_matrices.f90
	${FCCOMPILE_CHECK} -c -o $O/define_derivation_matrices.o ${FCFLAGS_f90} $S/define_derivation_matrices.f90

$O/get_attenuation.o: ${SETUP}/constants.h $S/get_attenuation.f90
	${FCCOMPILE_CHECK} -c -o $O/get_attenuation.o ${FCFLAGS_f90} $S/get_attenuation.f90

$O/get_backazimuth.o: ${SETUP}/constants.h $S/get_backazimuth.f90
	${FCCOMPILE_CHECK} -c -o $O/get_backazimuth.o ${FCFLAGS_f90} $S/get_backazimuth.f90

$O/get_cmt.o: ${SETUP}/constants.h $S/get_cmt.f90
	${FCCOMPILE_CHECK} -c -o $O/get_cmt.o ${FCFLAGS_f90} $S/get_cmt.f90

$O/multiply_arrays_source.o: ${SETUP}/constants.h $S/multiply_arrays_source.f90
	${FCCOMPILE_CHECK} -c -o $O/multiply_arrays_source.o ${FCFLAGS_f90} $S/multiply_arrays_source.f90

$O/netlib_specfun_erf.o: $S/netlib_specfun_erf.f90
	${FCCOMPILE_CHECK} -c -o $O/netlib_specfun_erf.o ${FCFLAGS_f90} $S/netlib_specfun_erf.f90

$O/recompute_jacobian.o: ${SETUP}/constants.h $S/recompute_jacobian.f90
	${FCCOMPILE_CHECK} -c -o $O/recompute_jacobian.o ${FCFLAGS_f90} $S/recompute_jacobian.f90

##
## specfem3D - use MPI here & dependent on values from mesher here
##
$O/assemble_MPI_central_cube.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/assemble_MPI_central_cube.f90
	${MPIFCCOMPILE_NO_CHECK} -c -o $O/assemble_MPI_central_cube.o ${FCFLAGS_f90} $S/assemble_MPI_central_cube.f90

$O/assemble_MPI_vector.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/assemble_MPI_vector.f90
	${MPIFCCOMPILE_NO_CHECK} -c -o $O/assemble_MPI_vector.o ${FCFLAGS_f90} $S/assemble_MPI_vector.f90

$O/assemble_MPI_vector_block.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/assemble_MPI_vector_block.f90
	${MPIFCCOMPILE_NO_CHECK} -c -o $O/assemble_MPI_vector_block.o ${FCFLAGS_f90} $S/assemble_MPI_vector_block.f90

$O/check_simulation_stability.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/check_simulation_stability.f90
	${MPIFCCOMPILE_NO_CHECK} -c -o $O/check_simulation_stability.o ${FCFLAGS_f90} $S/check_simulation_stability.f90

$O/initialize_simulation.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/initialize_simulation.f90
	${MPIFCCOMPILE_NO_CHECK} -c -o $O/initialize_simulation.o ${FCFLAGS_f90} $S/initialize_simulation.f90

$O/noise_tomography.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/noise_tomography.f90
	${MPIFCCOMPILE_NO_CHECK} -c -o $O/noise_tomography.o ${FCFLAGS_f90} $S/noise_tomography.f90

$O/prepare_timerun.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/prepare_timerun.f90
	${MPIFCCOMPILE_NO_CHECK} -c -o $O/prepare_timerun.o ${FCFLAGS_f90} $S/prepare_timerun.f90

$O/read_mesh_databases.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/read_mesh_databases.f90
	${MPIFCCOMPILE_NO_CHECK} -c -o $O/read_mesh_databases.o ${FCFLAGS_f90} $S/read_mesh_databases.f90

$O/setup_sources_receivers.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/setup_sources_receivers.f90
	${MPIFCCOMPILE_NO_CHECK} -c -o $O/setup_sources_receivers.o ${FCFLAGS_f90} $S/setup_sources_receivers.f90

$O/specfem3D.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/specfem3D.F90
	${MPIFCCOMPILE_NO_CHECK} -c -o $O/specfem3D.o ${FCFLAGS_f90} $S/specfem3D.F90

$O/write_movie_surface.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/write_movie_surface.f90
	${MPIFCCOMPILE_NO_CHECK} -c -o $O/write_movie_surface.o ${FCFLAGS_f90} $S/write_movie_surface.f90

$O/write_movie_volume.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/write_movie_volume.f90
	${MPIFCCOMPILE_NO_CHECK} -c -o $O/write_movie_volume.o ${FCFLAGS_f90} $S/write_movie_volume.f90

##
## specfem3D - non-dependent on values from mesher here
##
$O/assemble_MPI_central_cube_block.o: ${SETUP}/constants.h $S/assemble_MPI_central_cube_block.f90
	${MPIFCCOMPILE_NO_CHECK} -c -o $O/assemble_MPI_central_cube_block.o ${FCFLAGS_f90} $S/assemble_MPI_central_cube_block.f90

$O/assemble_MPI_scalar.o: ${SETUP}/constants.h $S/assemble_MPI_scalar.f90
	${MPIFCCOMPILE_NO_CHECK} -c -o $O/assemble_MPI_scalar.o ${FCFLAGS_f90} $S/assemble_MPI_scalar.f90

$O/assemble_MPI_scalar_block.o: ${SETUP}/constants.h $S/assemble_MPI_scalar_block.f90
	${MPIFCCOMPILE_NO_CHECK} -c -o $O/assemble_MPI_scalar_block.o ${FCFLAGS_f90} $S/assemble_MPI_scalar_block.f90

$O/fix_non_blocking_flags.o: ${SETUP}/constants.h $S/fix_non_blocking_flags.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/fix_non_blocking_flags.o ${FCFLAGS_f90} $S/fix_non_blocking_flags.f90

###
### specfem3D - regular MPI compilation options here
###
$O/create_central_cube_buffers.o: ${SETUP}/constants.h $S/create_central_cube_buffers.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/create_central_cube_buffers.o ${FCFLAGS_f90} $S/create_central_cube_buffers.f90

$O/get_event_info.o: ${SETUP}/constants.h $S/get_event_info.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/get_event_info.o ${FCFLAGS_f90} $S/get_event_info.f90

$O/locate_receivers.o: ${SETUP}/constants.h $S/locate_receivers.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/locate_receivers.o ${FCFLAGS_f90} $S/locate_receivers.f90

$O/locate_regular_points.o: ${SETUP}/constants.h $S/locate_regular_points.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/locate_regular_points.o ${FCFLAGS_f90} $S/locate_regular_points.f90

$O/locate_sources.o: ${SETUP}/constants.h $S/locate_sources.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/locate_sources.o ${FCFLAGS_f90} $S/locate_sources.f90

$O/read_arrays_buffers_solver.o: ${SETUP}/constants.h $S/read_arrays_buffers_solver.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/read_arrays_buffers_solver.o ${FCFLAGS_f90} $S/read_arrays_buffers_solver.f90

$O/write_seismograms.o: ${SETUP}/constants.h $S/write_seismograms.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/write_seismograms.o ${FCFLAGS_f90} $S/write_seismograms.f90

$O/write_output_ASCII.o: ${SETUP}/constants.h $S/write_output_ASCII.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/write_output_ASCII.o ${FCFLAGS_f90} $S/write_output_ASCII.f90

$O/write_output_SAC.o: ${SETUP}/constants.h $S/write_output_SAC.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/write_output_SAC.o ${FCFLAGS_f90} $S/write_output_SAC.f90

