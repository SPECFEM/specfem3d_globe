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

meshfem3D_TARGETS = \
	$E/xmeshfem3D \
	$(EMPTY_MACRO)

meshfem3D_OBJECTS = \
	$O/meshfem3D_par.check.o \
	$O/meshfem3D.check.o \
	$O/meshfem3D_models.check.o \
	$O/add_missing_nodes.check.o \
	$O/add_topography.check.o \
	$O/add_topography_410_650.check.o \
	$O/add_topography_cmb.check.o \
	$O/add_topography_icb.check.o \
	$O/assemble_MPI_central_cube_mesh.check.o \
	$O/assemble_MPI_scalar_mesh.check.o \
	$O/assemble_MPI_vector_mesh.check.o \
	$O/calc_jacobian.check.o \
	$O/compute_area.check.o \
	$O/compute_coordinates_grid.check.o \
	$O/compute_element_properties.check.o \
	$O/compute_volumes.check.o \
	$O/create_addressing.check.o \
	$O/create_central_cube.check.o \
	$O/create_central_cube_buffers.check.o \
	$O/create_chunk_buffers.check.o \
	$O/create_doubling_elements.check.o \
	$O/create_mass_matrices.check.o \
	$O/create_meshes.check.o \
	$O/create_MPI_interfaces.check.o \
	$O/create_regions_mesh.check.o \
	$O/create_regular_elements.check.o \
	$O/define_superbrick.check.o \
	$O/finalize_mesher.check.o \
	$O/fix_non_blocking_flags.check.o \
	$O/get_absorb.check.o \
	$O/get_ellipticity.check.o \
	$O/get_global.check.o \
	$O/get_jacobian_boundaries.check.o \
	$O/get_jacobian_discontinuities.check.o \
	$O/get_model.check.o \
	$O/get_MPI_1D_buffers.check.o \
	$O/get_MPI_cutplanes_eta.check.o \
	$O/get_MPI_cutplanes_xi.check.o \
	$O/get_MPI_interfaces.check.o \
	$O/get_perm_color.check.o \
	$O/get_shape2D.check.o \
	$O/get_shape3D.check.o \
	$O/heap_sort.check.o \
	$O/initialize_layers.check.o \
	$O/initialize_mesher.check.o \
	$O/lgndr.check.o \
	$O/model_1dref.check.o \
	$O/model_1066a.check.o \
	$O/model_ak135.check.o \
	$O/model_sea1d.check.o \
	$O/model_aniso_inner_core.check.o \
	$O/model_aniso_mantle.check.o \
	$O/model_atten3D_QRFSI12.check.o \
	$O/model_attenuation.check.o \
	$O/model_crust.check.o \
	$O/model_crustmaps.check.o \
	$O/model_eucrust.check.o \
	$O/model_epcrust.check.o \
	$O/model_gapp2.check.o \
	$O/model_gll.check.o \
	$O/model_heterogen_mantle.check.o \
	$O/model_iasp91.check.o \
	$O/model_jp1d.check.o \
	$O/model_jp3d.check.o \
	$O/model_ppm.check.o \
	$O/model_s20rts.check.o \
	$O/model_s40rts.check.o \
	$O/model_s362ani.check.o \
	$O/model_sea99_s.check.o \
	$O/moho_stretching.check.o \
	$O/save_arrays_solver.check.o \
	$O/setup_color_perm.check.o \
	$O/setup_counters.check.o \
	$O/setup_inner_outer.check.o \
	$O/setup_model.check.o \
	$O/setup_MPI_interfaces.check.o \
	$O/sort_array_coordinates.check.o \
	$O/stretching_function.check.o \
	$O/test_MPI_interfaces.check.o \
	$O/write_AVS_DX_global_chunks_data.check.o \
	$O/write_AVS_DX_global_data.check.o \
	$O/write_AVS_DX_global_faces_data.check.o \
	$O/write_AVS_DX_surface_data.check.o \
	$(EMPTY_MACRO)

meshfem3D_MODULES = \
	$(FC_MODDIR)/gapp2_mantle_model_constants.$(FC_MODEXT) \
	$(FC_MODDIR)/meshfem3d_models_par.$(FC_MODEXT) \
	$(FC_MODDIR)/module_ppm.$(FC_MODEXT) \
	$(EMPTY_MACRO)

# These files come from the shared directory
meshfem3D_SHARED_OBJECTS = \
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
	$O/exit_mpi.sharedmpi.o \
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
	$O/memory_eval.shared.o \
	$O/model_prem.shared.o \
	$O/model_topo_bathy.shared.o \
	$O/param_reader.cc.o \
	$O/read_compute_parameters.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/save_header_file.shared.o \
	$O/spline_routines.shared.o \
	$O/write_VTK_file.shared.o \
	$(EMPTY_MACRO)

XMESHFEM_OBJECTS = $(meshfem3D_SHARED_OBJECTS) $(meshfem3D_OBJECTS)

###
### ADIOS
###

adios_OBJECTS = \
  $O/adios_helpers.shared.o \
  $O/adios_manager.shared.o \
	$(EMPTY_MACRO)

adios_STUBS = \
	$O/adios_method_stubs.shared.o \
	$(EMPTY_MACRO)

# conditional adios linking
ifeq ($(ADIOS),yes)
XMESHFEM_OBJECTS += $(adios_OBJECTS)
else
XMESHFEM_OBJECTS += $(adios_STUBS)
endif

#######################################

####
#### rules for executables
####

${E}/xmeshfem3D: $(XMESHFEM_OBJECTS)
## use MPI here
	${MPIFCCOMPILE_CHECK} -o ${E}/xmeshfem3D $(XMESHFEM_OBJECTS) $(MPILIBS)

#######################################

## compilation directories
S := ${S_TOP}/src/meshfem3D
$(meshfem3D_OBJECTS): S = ${S_TOP}/src/meshfem3D

####
#### rule for each .o file below
####

$O/%.check.o: $S/%.f90 ${SETUP}/constants.h
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.check.o: $S/%.F90 ${SETUP}/constants.h
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.checkmpi.o: $S/%.f90 ${SETUP}/constants.h
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.checkmpi.o: $S/%.F90 ${SETUP}/constants.h
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

