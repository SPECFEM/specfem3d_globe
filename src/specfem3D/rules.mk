#=====================================================================
#
#          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
#          --------------------------------------------------
#
#     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
#                        Princeton University, USA
#                and CNRS / University of Marseille, France
#                 (there are currently many more authors!)
# (c) Princeton University and CNRS / University of Marseille, April 2014
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
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
	$O/compute_adj_source_frechet.solver.o \
	$O/convert_time.solver.o \
	$O/define_derivation_matrices.solver.o \
	$O/file_io_threads.cc.o \
	$O/force_ftz.cc.o \
	$O/get_backazimuth.solver.o \
	$O/get_cmt.solver.o \
	$O/get_force.solver.o \
	$O/get_event_info.solver.o \
	$O/make_gravity.solver.o \
	$O/netlib_specfun_erf.solver.o \
	$(EMPTY_MACRO)

# solver objects with statically allocated arrays; dependent upon
# values_from_mesher.h
specfem3D_OBJECTS += \
	$O/asdf_data.solverstatic_module.o \
	$O/comp_source_time_function.solverstatic.o \
	$O/specfem3D_par.solverstatic_module.o \
	$O/write_seismograms.solverstatic.o \
	$O/check_stability.solverstatic.o \
	$O/compute_add_sources.solverstatic.o \
	$O/compute_arrays_source.solverstatic.o \
	$O/compute_boundary_kernel.solverstatic.o \
	$O/compute_coupling.solverstatic.o \
	$O/compute_element.solverstatic.o \
	$O/compute_element_att_memory.solverstatic.o \
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
	$O/iterate_time_undoatt.solverstatic.o \
	$O/locate_receivers.solverstatic.o \
	$O/locate_regular_points.solverstatic.o \
	$O/locate_sources.solverstatic.o \
	$O/multiply_arrays_source.solverstatic.o \
	$O/noise_tomography.solverstatic.o \
	$O/prepare_attenuation.solverstatic.o \
	$O/prepare_elastic_elements.solverstatic.o \
	$O/prepare_gpu.solverstatic.o \
	$O/prepare_gravity.solverstatic.o \
	$O/prepare_movie.solverstatic.o \
	$O/prepare_noise.solverstatic.o \
	$O/prepare_oceans.solverstatic.o \
	$O/prepare_optimized_arrays.solverstatic.o \
	$O/prepare_stacey.solverstatic.o \
	$O/prepare_timerun.solverstatic.o \
	$O/prepare_vtk_window.solverstatic.o \
	$O/prepare_wavefields.solverstatic.o \
	$O/print_stf_file.solverstatic.o \
	$O/read_adjoint_sources.solverstatic.o \
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
	$O/update_displacement_LDDRK.solverstatic.o \
	$O/update_displacement_Newmark.solverstatic.o \
	$O/write_movie_output.solverstatic.o \
	$O/write_movie_volume.solverstatic.o \
	$O/write_movie_surface.solverstatic.o \
	$O/write_output_ASCII.solverstatic.o \
	$O/write_output_SAC.solverstatic.o \
	$(EMPTY_MACRO)

specfem3D_MODULES = \
	$(FC_MODDIR)/asdf_data.$(FC_MODEXT) \
	$(FC_MODDIR)/constants_solver.$(FC_MODEXT) \
	$(FC_MODDIR)/manager_adios_par.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_crustmantle.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_innercore.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_outercore.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_noise.$(FC_MODEXT) \
	$(FC_MODDIR)/specfem_par_movie.$(FC_MODEXT) \
	$(EMPTY_MACRO)

# These files come from the shared directory
specfem3D_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/adios_manager.shared_adios_module.o \
	$O/auto_ner.shared.o \
	$O/binary_c_io.cc.o \
	$O/broadcast_computed_parameters.shared.o \
	$O/calendar.shared.o \
	$O/count_elements.shared.o \
	$O/count_number_of_sources.shared.o \
	$O/count_points.shared.o \
	$O/create_name_database.shared.o \
	$O/define_all_layers.shared.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/get_model_parameters.shared.o \
	$O/get_timestep_and_layers.shared.o \
	$O/gll_library.shared.o \
	$O/hex_nodes.shared.o \
	$O/init_openmp.shared.o \
	$O/intgrl.shared.o \
	$O/lagrange_poly.shared.o \
	$O/make_ellipticity.shared.o \
	$O/model_prem.shared.o \
	$O/model_topo_bathy.shared.o \
	$O/parallel.sharedmpi.o \
	$O/param_reader.cc.o \
	$O/read_compute_parameters.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/recompute_jacobian.shared.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/spline_routines.shared.o \
	$O/write_VTK_file.shared.o \
	$(EMPTY_MACRO)

###
### GPU
###

gpu_specfem3D_STUBS = \
	$O/specfem3D_gpu_method_stubs.cc.o \
	$(EMPTY_MACRO)

ifdef NO_GPU
specfem3D_OBJECTS += $(gpu_specfem3D_STUBS)
else
specfem3D_OBJECTS += $(gpu_specfem3D_OBJECTS)
endif

###
### ADIOS
###

adios_specfem3D_OBJECTS = \
	$O/read_arrays_solver_adios.solverstatic_adios.o \
	$O/read_attenuation_adios.solverstatic_adios.o \
	$O/read_forward_arrays_adios.solverstatic_adios.o \
	$O/read_mesh_databases_adios.solverstatic_adios.o \
	$O/save_forward_arrays_adios.solverstatic_adios.o \
	$O/save_kernels_adios.solverstatic_adios.o \
	$(EMPTY_MACRO)

adios_specfem3D_SHARED_OBJECTS = \
	$O/adios_helpers_definitions.shared_adios_module.o \
	$O/adios_helpers_writers.shared_adios_module.o \
	$O/adios_helpers.shared_adios.o \
	$(EMPTY_MACRO)

adios_specfem3D_SHARED_STUBS = \
	$O/adios_method_stubs.cc.o \
	$(EMPTY_MACRO)

# conditional adios linking
ifeq ($(ADIOS),yes)
specfem3D_OBJECTS += $(adios_specfem3D_OBJECTS)
specfem3D_SHARED_OBJECTS += $(adios_specfem3D_SHARED_OBJECTS)
else
specfem3D_SHARED_OBJECTS += $(adios_specfem3D_SHARED_STUBS)
endif

###
### ASDF
###

asdf_specfem3D_OBJECTS = \
	$O/write_output_ASDF.solverstatic.o \
	$O/read_adjoint_sources_ASDF.solverstatic.o \
	$(EMPTY_MACRO)

asdf_specfem3D_SHARED_OBJECTS = \
	$O/asdf_manager.shared_asdf.o \
	$(EMPTY_MACRO)

asdf_specfem3D_STUBS = \
	$(EMPTY_MACRO)

asdf_specfem3D_SHARED_STUBS = \
	$O/asdf_method_stubs.cc.o \
	$(EMPTY_MACRO)

# conditional asdf linking
ifeq ($(ASDF),yes)
specfem3D_OBJECTS += $(asdf_specfem3D_OBJECTS)
specfem3D_SHARED_OBJECTS += $(asdf_specfem3D_SHARED_OBJECTS)
else
specfem3D_OBJECTS += $(asdf_specfem3D_STUBS)
specfem3D_SHARED_OBJECTS += ${asdf_specfem3D_SHARED_STUBS}
endif

#
# conditional CEM model
ifeq ($(CEM),yes)
specfem3D_OBJECTS += $O/read_write_netcdf.checknetcdf.o
endif

###
### VTK
###

vtk_specfem3D_OBJECTS = \
  $O/visual_vtk.visualcc.o \
	$(EMPTY_MACRO)

vtk_specfem3D_STUBS = \
	$O/visual_vtk_stubs.visualc.o \
	$(EMPTY_MACRO)

# conditional adios linking
ifeq ($(VTK),yes)
specfem3D_OBJECTS += $(vtk_specfem3D_OBJECTS)
else
specfem3D_OBJECTS += $(vtk_specfem3D_STUBS)
endif

###
### LIBXSMM
###
# conditional module
ifeq ($(XSMM),yes)
specfem3D_MODULES += $(FC_MODDIR)/my_libxsmm.$(FC_MODEXT)
endif

#######################################

####
#### rules for executables
####

specfem3D_ALL_OBJECTS = $(specfem3D_OBJECTS) $(specfem3D_SHARED_OBJECTS)

ifeq ($(CUDA),yes)
SPECFEM_LINK_FLAGS = $(LDFLAGS) $(MPILIBS) $(CUDA_LINK) $(LIBS)
else
SPECFEM_LINK_FLAGS = $(LDFLAGS) $(MPILIBS) $(LIBS)
endif

# non-gpu or opencl
ifeq ($(ASDF),yes)
SPECFEM_LINK_FLAGS += $(ASDF_LIBS) -lhdf5hl_fortran -lhdf5_hl -lhdf5 -lstdc++
endif


${E}/xspecfem3D: $(specfem3D_ALL_OBJECTS)
	@echo ""
	@echo "building xspecfem3D $(BUILD_VERSION_TXT)"
	@echo ""


## use MPI here
## DK DK add OpenMP compiler flag here if needed
#	${MPIFCCOMPILE_CHECK} -qsmp=omp -o ${E}/xspecfem3D $(specfem3D_ALL_OBJECTS) $(LDFLAGS) $(MPILIBS) $(LIBS)
	${FCLINK} -o ${E}/xspecfem3D $(specfem3D_ALL_OBJECTS) $(SPECFEM_LINK_FLAGS)
	@echo ""

#######################################

## compilation directories
S := ${S_TOP}/src/specfem3D
$(specfem3D_OBJECTS): S = ${S_TOP}/src/specfem3D

###
### additional dependencies
###

$O/write_output_ASDF.solverstatic.o: $O/asdf_data.solverstatic_module.o
$O/compute_arrays_source.solverstatic.o: $O/write_seismograms.solverstatic.o
$O/iterate_time.solverstatic.o: $O/write_seismograms.solverstatic.o
$O/iterate_time_undoatt.solverstatic.o: $O/write_seismograms.solverstatic.o
$O/locate_receivers.solverstatic.o: $O/write_seismograms.solverstatic.o
$O/read_adjoint_sources.solverstatic.o: $O/write_seismograms.solverstatic.o

$O/specfem3D_par.solverstatic_module.o: $O/adios_manager.shared_adios_module.o

# Version file
$O/initialize_simulation.solverstatic.o: ${SETUP}/version.fh

###
### specfem3D - optimized flags and dependence on values from mesher here
###
$O/%.solverstatic_module.o: $S/%.F90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.solverstatic_module.o: $S/%.f90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.solverstatic.o: $S/%.f90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o $O/specfem3D_par.solverstatic_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.solverstatic.o: $S/%.F90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o $O/specfem3D_par.solverstatic_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.solverstatic_openmp.o: $S/%.f90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o $O/specfem3D_par.solverstatic_module.o
## DK DK add OpenMP compiler flag here if needed
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -qsmp=omp -o $@ $<


$O/%.solverstatic_adios.o: $S/%.f90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o $O/specfem3D_par.solverstatic_module.o $O/adios_helpers.shared_adios.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.solverstatic_adios.o: $S/%.F90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o $O/specfem3D_par.solverstatic_module.o $O/adios_helpers.shared_adios.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

###
### no dependence on values from mesher here
###

$O/%.solver.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.solver.o: $S/%.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.cc.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) -o $@ $<

###
### VTK compilation
###

$O/%.visualcc.o: $S/%.cpp ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(MPI_INCLUDES) -o $@ $<

$O/%.visualc.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(MPI_INCLUDES) -o $@ $<

###
### CEM
###
$O/%.checknetcdf.o: $S/%.f90 $O/shared_par.shared_module.o $O/specfem3D_par.solverstatic_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} $(NETCDF_INCLUDE) -c -o $@ $<
