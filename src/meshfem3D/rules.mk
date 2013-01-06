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

meshfem3D_TARGETS = \
	$E/xmeshfem3D \
	$(EMPTY_MACRO)

meshfem3D_OBJECTS = \
	$O/meshfem3D.o \
	$O/meshfem3D_models.o \
	$O/add_missing_nodes.o \
	$O/add_topography.o \
	$O/add_topography_410_650.o \
	$O/add_topography_cmb.o \
	$O/add_topography_icb.o \
	$O/calc_jacobian.o \
	$O/compute_coordinates_grid.o \
	$O/compute_element_properties.o \
	$O/create_central_cube.o \
	$O/create_chunk_buffers.o \
	$O/create_doubling_elements.o \
	$O/create_mass_matrices.o \
	$O/create_regions_mesh.o \
	$O/get_perm_color.o \
	$O/create_regular_elements.o \
	$O/define_superbrick.o \
	$O/get_absorb.o \
	$O/get_ellipticity.o \
	$O/get_global.o \
	$O/get_jacobian_boundaries.o \
	$O/get_jacobian_discontinuities.o \
	$O/get_model.o \
	$O/get_MPI_1D_buffers.o \
	$O/get_MPI_cutplanes_eta.o \
	$O/get_MPI_cutplanes_xi.o \
	$O/get_shape2D.o \
	$O/get_shape3D.o \
	$O/lgndr.o \
	$O/model_sea99_s.o \
	$O/model_1066a.o \
	$O/model_ak135.o \
	$O/model_aniso_inner_core.o \
	$O/model_aniso_mantle.o \
	$O/model_atten3D_QRFSI12.o \
	$O/model_attenuation.o \
	$O/model_crust.o \
	$O/model_eucrust.o \
	$O/model_epcrust.o \
	$O/model_crustmaps.o \
	$O/model_gll.o \
	$O/model_heterogen_mantle.o \
	$O/model_iasp91.o \
	$O/model_jp1d.o \
	$O/model_jp3d.o \
	$O/model_ppm.o \
	$O/model_gapp2.o \
	$O/model_1dref.o \
	$O/model_s20rts.o \
	$O/model_s40rts.o \
	$O/model_s362ani.o \
	$O/model_sea1d.o \
	$O/moho_stretching.o \
	$O/save_arrays_solver.o \
	$O/sort_array_coordinates.o \
	$O/stretching_function.o \
	$O/write_AVS_DX_global_chunks_data.o \
	$O/write_AVS_DX_global_data.o \
	$O/write_AVS_DX_global_faces_data.o \
	$O/write_AVS_DX_surface_data.o \
	$(EMPTY_MACRO)

meshfem3D_MODULES = \
	$(FC_MODDIR)/gapp2_mantle_model_constants.$(FC_MODEXT) \
	$(FC_MODDIR)/meshfem3d_models_par.$(FC_MODEXT) \
	$(FC_MODDIR)/module_ppm.$(FC_MODEXT) \
	$(EMPTY_MACRO)

# These files come from the shared directory
meshfem3D_SHARED_OBJECTS = \
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
	$O/memory_eval.o \
	$O/model_prem.o \
	$O/model_topo_bathy.o \
	$O/param_reader.o \
	$O/read_compute_parameters.o \
	$O/read_parameter_file.o \
	$O/read_value_parameters.o \
	$O/reduce.o \
	$O/rthetaphi_xyz.o \
	$O/save_header_file.o \
	$O/spline_routines.o \
	$(EMPTY_MACRO)

#######################################

####
#### rules for executables
####

${E}/xmeshfem3D: $(meshfem3D_OBJECTS) $(meshfem3D_SHARED_OBJECTS)
## use MPI here
	${MPIFCCOMPILE_CHECK} -o ${E}/xmeshfem3D $(meshfem3D_OBJECTS) $(meshfem3D_SHARED_OBJECTS) $(MPILIBS)

#######################################

## compilation directories
S := ${S_TOP}/src/meshfem3D
$(meshfem3D_OBJECTS): S = ${S_TOP}/src/meshfem3D

####
#### rule for each .o file below
####

###
### meshfem3D objects
###

$O/add_missing_nodes.o: ${SETUP}/constants.h $S/add_missing_nodes.f90
	${FCCOMPILE_CHECK} -c -o $O/add_missing_nodes.o ${FCFLAGS_f90} $S/add_missing_nodes.f90

$O/add_topography.o: ${SETUP}/constants.h $S/add_topography.f90
	${FCCOMPILE_CHECK} -c -o $O/add_topography.o ${FCFLAGS_f90} $S/add_topography.f90

$O/add_topography_410_650.o: ${SETUP}/constants.h $S/add_topography_410_650.f90
	${FCCOMPILE_CHECK} -c -o $O/add_topography_410_650.o ${FCFLAGS_f90} $S/add_topography_410_650.f90

$O/add_topography_cmb.o: ${SETUP}/constants.h $S/add_topography_cmb.f90
	${FCCOMPILE_CHECK} -c -o $O/add_topography_cmb.o ${FCFLAGS_f90} $S/add_topography_cmb.f90

$O/add_topography_icb.o: ${SETUP}/constants.h $S/add_topography_icb.f90
	${FCCOMPILE_CHECK} -c -o $O/add_topography_icb.o ${FCFLAGS_f90} $S/add_topography_icb.f90

$O/calc_jacobian.o: ${SETUP}/constants.h $S/calc_jacobian.f90
	${FCCOMPILE_CHECK} -c -o $O/calc_jacobian.o ${FCFLAGS_f90} $S/calc_jacobian.f90

$O/compute_coordinates_grid.o: ${SETUP}/constants.h $S/compute_coordinates_grid.f90
	${FCCOMPILE_CHECK} -c -o $O/compute_coordinates_grid.o ${FCFLAGS_f90} $S/compute_coordinates_grid.f90

$O/compute_element_properties.o: ${SETUP}/constants.h $S/compute_element_properties.f90 $O/meshfem3D_models.o
	${FCCOMPILE_CHECK} -c -o $O/compute_element_properties.o ${FCFLAGS_f90} $S/compute_element_properties.f90

$O/create_central_cube.o: ${SETUP}/constants.h $S/create_central_cube.f90 $O/meshfem3D_models.o
	${FCCOMPILE_CHECK} -c -o $O/create_central_cube.o ${FCFLAGS_f90} $S/create_central_cube.f90

$O/create_doubling_elements.o: ${SETUP}/constants.h $S/create_doubling_elements.f90 $O/meshfem3D_models.o
	${FCCOMPILE_CHECK} -c -o $O/create_doubling_elements.o ${FCFLAGS_f90} $S/create_doubling_elements.f90

$O/create_mass_matrices.o: ${SETUP}/constants.h $S/create_mass_matrices.f90 $O/meshfem3D_models.o
	${FCCOMPILE_CHECK} -c -o $O/create_mass_matrices.o ${FCFLAGS_f90} $S/create_mass_matrices.f90

$O/create_regions_mesh.o: ${SETUP}/constants.h $S/create_regions_mesh.F90 $O/meshfem3D_models.o
	${MPIFCCOMPILE_CHECK} -c -o $O/create_regions_mesh.o ${FCFLAGS_f90} $S/create_regions_mesh.F90

$O/get_perm_color.o: ${SETUP}/constants.h $S/get_perm_color.f90
	${FCCOMPILE_CHECK} -c -o $O/get_perm_color.o ${FCFLAGS_f90} $S/get_perm_color.f90

$O/create_regular_elements.o: ${SETUP}/constants.h $S/create_regular_elements.f90 $O/meshfem3D_models.o
	${FCCOMPILE_CHECK} -c -o $O/create_regular_elements.o ${FCFLAGS_f90} $S/create_regular_elements.f90

$O/define_superbrick.o: ${SETUP}/constants.h $S/define_superbrick.f90
	${FCCOMPILE_CHECK} -c -o $O/define_superbrick.o ${FCFLAGS_f90} $S/define_superbrick.f90

$O/get_absorb.o: ${SETUP}/constants.h $S/get_absorb.f90
	${FCCOMPILE_CHECK} -c -o $O/get_absorb.o ${FCFLAGS_f90} $S/get_absorb.f90

$O/get_ellipticity.o: ${SETUP}/constants.h $S/get_ellipticity.f90
	${FCCOMPILE_CHECK} -c -o $O/get_ellipticity.o ${FCFLAGS_f90} $S/get_ellipticity.f90

$O/get_global.o: ${SETUP}/constants.h $S/get_global.f90
	${FCCOMPILE_CHECK} -c -o $O/get_global.o ${FCFLAGS_f90} $S/get_global.f90

$O/get_jacobian_boundaries.o: ${SETUP}/constants.h $S/get_jacobian_boundaries.f90
	${FCCOMPILE_CHECK} -c -o $O/get_jacobian_boundaries.o ${FCFLAGS_f90} $S/get_jacobian_boundaries.f90

$O/get_jacobian_discontinuities.o: ${SETUP}/constants.h $S/get_jacobian_discontinuities.f90
	${FCCOMPILE_CHECK} -c -o $O/get_jacobian_discontinuities.o ${FCFLAGS_f90} $S/get_jacobian_discontinuities.f90

$O/get_model.o: ${SETUP}/constants.h $S/get_model.f90 $O/meshfem3D_models.o
	${FCCOMPILE_CHECK} -c -o $O/get_model.o ${FCFLAGS_f90} $S/get_model.f90

$O/get_MPI_1D_buffers.o: ${SETUP}/constants.h $S/get_MPI_1D_buffers.f90
	${FCCOMPILE_CHECK} -c -o $O/get_MPI_1D_buffers.o ${FCFLAGS_f90} $S/get_MPI_1D_buffers.f90

$O/get_MPI_cutplanes_xi.o: ${SETUP}/constants.h $S/get_MPI_cutplanes_xi.f90
	${FCCOMPILE_CHECK} -c -o $O/get_MPI_cutplanes_xi.o ${FCFLAGS_f90} $S/get_MPI_cutplanes_xi.f90

$O/get_MPI_cutplanes_eta.o: ${SETUP}/constants.h $S/get_MPI_cutplanes_eta.f90
	${FCCOMPILE_CHECK} -c -o $O/get_MPI_cutplanes_eta.o ${FCFLAGS_f90} $S/get_MPI_cutplanes_eta.f90

$O/get_shape2D.o: ${SETUP}/constants.h $S/get_shape2D.f90
	${FCCOMPILE_CHECK} -c -o $O/get_shape2D.o ${FCFLAGS_f90} $S/get_shape2D.f90

$O/get_shape3D.o: ${SETUP}/constants.h $S/get_shape3D.f90
	${FCCOMPILE_CHECK} -c -o $O/get_shape3D.o ${FCFLAGS_f90} $S/get_shape3D.f90

$O/lgndr.o: ${SETUP}/constants.h $S/lgndr.f90
	${FCCOMPILE_CHECK} -c -o $O/lgndr.o ${FCFLAGS_f90} $S/lgndr.f90

$O/model_iasp91.o: ${SETUP}/constants.h $S/model_iasp91.f90
	${FCCOMPILE_CHECK} -c -o $O/model_iasp91.o ${FCFLAGS_f90} $S/model_iasp91.f90

$O/model_1066a.o: ${SETUP}/constants.h $S/model_1066a.f90
	${FCCOMPILE_CHECK} -c -o $O/model_1066a.o ${FCFLAGS_f90} $S/model_1066a.f90

$O/model_ak135.o: ${SETUP}/constants.h $S/model_ak135.f90
	${FCCOMPILE_CHECK} -c -o $O/model_ak135.o ${FCFLAGS_f90} $S/model_ak135.f90

$O/model_aniso_inner_core.o: ${SETUP}/constants.h $S/model_aniso_inner_core.f90
	${FCCOMPILE_CHECK} -c -o $O/model_aniso_inner_core.o ${FCFLAGS_f90} $S/model_aniso_inner_core.f90

$O/model_1dref.o: ${SETUP}/constants.h $S/model_1dref.f90
	${FCCOMPILE_CHECK} -c -o $O/model_1dref.o ${FCFLAGS_f90} $S/model_1dref.f90

$O/model_jp1d.o: ${SETUP}/constants.h $S/model_jp1d.f90
	${FCCOMPILE_CHECK} -c -o $O/model_jp1d.o ${FCFLAGS_f90} $S/model_jp1d.f90

$O/model_sea1d.o: ${SETUP}/constants.h $S/model_sea1d.f90
	${FCCOMPILE_CHECK} -c -o $O/model_sea1d.o ${FCFLAGS_f90} $S/model_sea1d.f90

$O/moho_stretching.o: ${SETUP}/constants.h $S/moho_stretching.f90 $O/meshfem3D_models.o
	${FCCOMPILE_CHECK} -c -o $O/moho_stretching.o ${FCFLAGS_f90} $S/moho_stretching.f90

$O/save_arrays_solver.o: ${SETUP}/constants.h $S/save_arrays_solver.f90
	${FCCOMPILE_CHECK} -c -o $O/save_arrays_solver.o ${FCFLAGS_f90} $S/save_arrays_solver.f90

$O/sort_array_coordinates.o: ${SETUP}/constants.h $S/sort_array_coordinates.f90
	${FCCOMPILE_CHECK} -c -o $O/sort_array_coordinates.o ${FCFLAGS_f90} $S/sort_array_coordinates.f90

$O/stretching_function.o: ${SETUP}/constants.h $S/stretching_function.f90
	${FCCOMPILE_CHECK} -c -o $O/stretching_function.o ${FCFLAGS_f90} $S/stretching_function.f90

$O/write_AVS_DX_global_faces_data.o: ${SETUP}/constants.h $S/write_AVS_DX_global_faces_data.f90
	${FCCOMPILE_CHECK} -c -o $O/write_AVS_DX_global_faces_data.o ${FCFLAGS_f90} $S/write_AVS_DX_global_faces_data.f90

$O/write_AVS_DX_global_chunks_data.o: ${SETUP}/constants.h $S/write_AVS_DX_global_chunks_data.f90
	${FCCOMPILE_CHECK} -c -o $O/write_AVS_DX_global_chunks_data.o ${FCFLAGS_f90} $S/write_AVS_DX_global_chunks_data.f90

$O/write_AVS_DX_surface_data.o: ${SETUP}/constants.h $S/write_AVS_DX_surface_data.f90
	${FCCOMPILE_CHECK} -c -o $O/write_AVS_DX_surface_data.o ${FCFLAGS_f90} $S/write_AVS_DX_surface_data.f90

$O/write_AVS_DX_global_data.o: ${SETUP}/constants.h $S/write_AVS_DX_global_data.f90
	${FCCOMPILE_CHECK} -c -o $O/write_AVS_DX_global_data.o ${FCFLAGS_f90} $S/write_AVS_DX_global_data.f90

##
## meshfem3D objects with mpi compilation
##
$O/create_chunk_buffers.o: ${SETUP}/constants.h $S/create_chunk_buffers.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/create_chunk_buffers.o ${FCFLAGS_f90} $S/create_chunk_buffers.f90

$O/meshfem3D.o: ${SETUP}/constants.h $S/meshfem3D.f90 $O/meshfem3D_models.o
	${MPIFCCOMPILE_CHECK} -c -o $O/meshfem3D.o ${FCFLAGS_f90} $S/meshfem3D.f90

$O/meshfem3D_models.o: ${SETUP}/constants.h $S/meshfem3D_models.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/meshfem3D_models.o ${FCFLAGS_f90} $S/meshfem3D_models.f90

$O/model_aniso_mantle.o: ${SETUP}/constants.h $S/model_aniso_mantle.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/model_aniso_mantle.o ${FCFLAGS_f90} $S/model_aniso_mantle.f90

$O/model_atten3D_QRFSI12.o: ${SETUP}/constants.h $S/model_atten3D_QRFSI12.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/model_atten3D_QRFSI12.o ${FCFLAGS_f90} $S/model_atten3D_QRFSI12.f90

$O/model_attenuation.o: ${SETUP}/constants.h $S/model_attenuation.f90 $O/model_ak135.o $O/model_1066a.o $O/model_1dref.o
	${MPIFCCOMPILE_CHECK} -c -o $O/model_attenuation.o ${FCFLAGS_f90} $S/model_attenuation.f90

$O/model_crust.o: ${SETUP}/constants.h $S/model_crust.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/model_crust.o ${FCFLAGS_f90} $S/model_crust.f90

$O/model_eucrust.o: ${SETUP}/constants.h $S/model_eucrust.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/model_eucrust.o ${FCFLAGS_f90} $S/model_eucrust.f90

$O/model_epcrust.o: ${SETUP}/constants.h $S/model_epcrust.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/model_epcrust.o ${FCFLAGS_f90} $S/model_epcrust.f90

$O/model_crustmaps.o: ${SETUP}/constants.h $S/model_crustmaps.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/model_crustmaps.o ${FCFLAGS_f90} $S/model_crustmaps.f90

$O/model_gll.o: ${SETUP}/constants.h $S/model_gll.f90 $O/meshfem3D_models.o
	${MPIFCCOMPILE_CHECK} -c -o $O/model_gll.o ${FCFLAGS_f90} $S/model_gll.f90

$O/model_heterogen_mantle.o: ${SETUP}/constants.h $S/model_heterogen_mantle.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/model_heterogen_mantle.o ${FCFLAGS_f90} $S/model_heterogen_mantle.f90

$O/model_jp3d.o: ${SETUP}/constants.h $S/model_jp3d.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/model_jp3d.o ${FCFLAGS_f90} $S/model_jp3d.f90

$O/model_ppm.o: ${SETUP}/constants.h $S/model_ppm.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/model_ppm.o ${FCFLAGS_f90} $S/model_ppm.f90

$O/model_gapp2.o: ${SETUP}/constants.h $S/model_gapp2.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/model_gapp2.o ${FCFLAGS_f90} $S/model_gapp2.f90

$O/model_s20rts.o: ${SETUP}/constants.h $S/model_s20rts.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/model_s20rts.o ${FCFLAGS_f90} $S/model_s20rts.f90

$O/model_s40rts.o: ${SETUP}/constants.h $S/model_s40rts.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/model_s40rts.o ${FCFLAGS_f90} $S/model_s40rts.f90

$O/model_s362ani.o: ${SETUP}/constants.h $S/model_s362ani.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/model_s362ani.o ${FCFLAGS_f90} $S/model_s362ani.f90

$O/model_sea99_s.o: ${SETUP}/constants.h $S/model_sea99_s.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/model_sea99_s.o ${FCFLAGS_f90} $S/model_sea99_s.f90

