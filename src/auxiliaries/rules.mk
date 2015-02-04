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

auxiliaries_TARGETS = \
	$E/xconvolve_source_timefunction \
	$E/xdetect_duplicates_stations_file \
	$E/xcombine_AVS_DX \
	$E/xcombine_paraview_strain_data \
	$E/xcombine_vol_data \
	$E/xcombine_vol_data_vtk \
	$E/xcombine_surf_data \
	$E/xcreate_movie_AVS_DX \
	$E/xcreate_movie_GMT_global \
	$E/xextract_database \
	$(EMPTY_MACRO)

ifeq ($(ADIOS),yes)
auxiliaries_TARGETS += \
	$E/xcombine_vol_data_adios \
	$E/xcombine_vol_data_vtk_adios \
	$(EMPTY_MACRO)
endif

auxiliaries_OBJECTS = \
	$(xconvolve_source_timefunction_OBJECTS) \
	$(xdetect_duplicates_stations_file_OBJECTS) \
	$(xcombine_AVS_DX_OBJECTS) \
	$(xcombine_paraview_strain_data_OBJECTS) \
	$(xcombine_surf_data_OBJECTS) \
	$(xcombine_vol_data_OBJECTS) \
	$(xcombine_vol_data_adios_OBJECTS) \
	$(xcombine_vol_data_vtk_OBJECTS) \
	$(xcombine_vol_data_vtk_adios_OBJECTS) \
	$(xcreate_movie_AVS_DX_OBJECTS) \
	$(xcreate_movie_GMT_global_OBJECTS) \
	$(EMPTY_MACRO)

# These files come from the shared directory
auxiliaries_SHARED_OBJECTS = \
	$(xconvolve_source_timefunction_SHARED_OBJECTS) \
	$(xdetect_duplicates_stations_file_SHARED_OBJECTS) \
	$(xcombine_AVS_DX_SHARED_OBJECTS) \
	$(xcombine_paraview_strain_data_SHARED_OBJECTS) \
	$(xcombine_surf_data_SHARED_OBJECTS) \
	$(xcombine_vol_data_SHARED_OBJECTS) \
	$(xcombine_vol_data_adios_SHARED_OBJECTS) \
	$(xcombine_vol_data_vtk_SHARED_OBJECTS) \
	$(xcombine_vol_data_vtk_adios_SHARED_OBJECTS) \
	$(xcreate_movie_AVS_DX_SHARED_OBJECTS) \
	$(xcreate_movie_GMT_global_SHARED_OBJECTS) \
	$(xextract_database_SHARED_OBJECTS) \
	$(EMPTY_MACRO)

####
#### rules for executables
####

.PHONY: all_aux aux movies

all_aux: $(auxiliaries_TARGETS)

aux: $(auxiliaries_TARGETS)

movies: \
	$E/xcombine_paraview_strain_data \
	$E/xcreate_movie_AVS_DX \
	$E/xcreate_movie_GMT_global \
	$(EMPTY_MACRO)

#######################################

####
#### rules for each program follow
####

#######################################

xconvolve_source_timefunction_OBJECTS = \
	$O/convolve_source_timefunction.aux.o \
	$(EMPTY_MACRO)

xconvolve_source_timefunction_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$(EMPTY_MACRO)

${E}/xconvolve_source_timefunction: $(xconvolve_source_timefunction_OBJECTS) $(xconvolve_source_timefunction_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} -o $@ $+

#######################################

xdetect_duplicates_stations_file_OBJECTS = \
  $O/detect_duplicates_stations_file.aux.o \
  $(EMPTY_MACRO)

xdetect_duplicates_stations_file_SHARED_OBJECTS = \
  $O/shared_par.shared_module.o \
  $(EMPTY_MACRO)

${E}/xdetect_duplicates_stations_file: $(xdetect_duplicates_stations_file_OBJECTS) $(xdetect_duplicates_stations_file_SHARED_OBJECTS)
  ${FCCOMPILE_CHECK} -o $@ $+

#######################################

xcombine_AVS_DX_OBJECTS = \
	$O/combine_AVS_DX.aux.o \
	$O/get_cmt.solver.o \
	$(EMPTY_MACRO)

xcombine_AVS_DX_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/auto_ner.shared.o \
	$O/calendar.shared.o \
	$O/count_elements.shared.o \
	$O/count_number_of_sources.shared.o \
	$O/count_points.shared.o \
	$O/create_serial_name_database.shared.o \
	$O/define_all_layers.shared.o \
	$O/get_model_parameters.shared.o \
	$O/get_timestep_and_layers.shared.o \
	$O/param_reader.cc.o \
	$O/read_compute_parameters.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$(EMPTY_MACRO)

${E}/xcombine_AVS_DX: $(xcombine_AVS_DX_OBJECTS) $(xcombine_AVS_DX_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} -o $@ $+

#######################################

xcombine_paraview_strain_data_OBJECTS = \
	$O/combine_paraview_strain_data.auxsolver.o \
	$(EMPTY_MACRO)

xcombine_paraview_strain_data_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/binary_c_io.cc.o \
	$(EMPTY_MACRO)

${E}/xcombine_paraview_strain_data: $(xcombine_paraview_strain_data_OBJECTS) $(xcombine_paraview_strain_data_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} -o $@ $+

#######################################

xcombine_surf_data_OBJECTS = \
	$O/combine_surf_data.auxsolver.o \
	$(EMPTY_MACRO)

xcombine_surf_data_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/binary_c_io.cc.o \
	$O/flush_system.shared.o \
	$(EMPTY_MACRO)

${E}/xcombine_surf_data: $(xcombine_surf_data_OBJECTS) $(xcombine_surf_data_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} -o $@ $+

#######################################

xcombine_vol_data_OBJECTS = \
	$O/combine_vol_data.auxsolver.o \
	$(EMPTY_MACRO)

xcombine_vol_data_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/binary_c_io.cc.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/intgrl.shared.o \
	$O/make_ellipticity.shared.o \
	$O/model_prem.shared.o \
	$O/parallel.sharedmpi.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/param_reader.cc.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/spline_routines.shared.o \
	$(EMPTY_MACRO)

${E}/xcombine_vol_data: $(xcombine_vol_data_OBJECTS) $(xcombine_vol_data_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)

#######################################

xcombine_vol_data_adios_OBJECTS = \
	$O/combine_vol_data.auxadios.o \
	$O/combine_vol_data_adios_impl.auxmpi.o \
	$(EMPTY_MACRO)

xcombine_vol_data_adios_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/binary_c_io.cc.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/intgrl.shared.o \
	$O/make_ellipticity.shared.o \
	$O/model_prem.shared.o \
	$O/parallel.sharedmpi.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/param_reader.cc.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/spline_routines.shared.o \
	$(EMPTY_MACRO)

$O/combine_vol_data.auxadios.o: $O/combine_vol_data_adios_impl.auxmpi.o

${E}/xcombine_vol_data_adios: $(xcombine_vol_data_adios_OBJECTS) $(xcombine_vol_data_adios_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)

#######################################

xcombine_vol_data_vtk_OBJECTS = \
	$O/combine_vol_data.auxsolver_vtk.o \
	$(EMPTY_MACRO)

xcombine_vol_data_vtk_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/binary_c_io.cc.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/intgrl.shared.o \
	$O/make_ellipticity.shared.o \
	$O/model_prem.shared.o \
	$O/parallel.sharedmpi.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/param_reader.cc.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/spline_routines.shared.o \
	$(EMPTY_MACRO)

${E}/xcombine_vol_data_vtk: $(xcombine_vol_data_vtk_OBJECTS) $(xcombine_vol_data_vtk_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)

#######################################

xcombine_vol_data_vtk_adios_OBJECTS = \
	$O/combine_vol_data.auxadios_vtk.o \
	$O/combine_vol_data_adios_impl.auxmpi.o \
	$(EMPTY_MACRO)

xcombine_vol_data_vtk_adios_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/binary_c_io.cc.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/intgrl.shared.o \
	$O/make_ellipticity.shared.o \
	$O/model_prem.shared.o \
	$O/parallel.sharedmpi.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/param_reader.cc.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/spline_routines.shared.o \
	$(EMPTY_MACRO)

$O/combine_vol_data.auxadios_vtk.o: $O/combine_vol_data_adios_impl.auxmpi.o

${E}/xcombine_vol_data_vtk_adios: $(xcombine_vol_data_vtk_adios_OBJECTS) $(xcombine_vol_data_vtk_adios_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)

#######################################

xcreate_movie_AVS_DX_OBJECTS = \
	$O/create_movie_AVS_DX.aux.o \
	$(EMPTY_MACRO)

xcreate_movie_AVS_DX_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/auto_ner.shared.o \
	$O/count_elements.shared.o \
	$O/count_number_of_sources.shared.o \
	$O/count_points.shared.o \
	$O/define_all_layers.shared.o \
	$O/get_global.shared.o \
	$O/get_model_parameters.shared.o \
	$O/get_timestep_and_layers.shared.o \
	$O/param_reader.cc.o \
	$O/read_compute_parameters.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/sort_array_coordinates.shared.o \
	$(EMPTY_MACRO)

${E}/xcreate_movie_AVS_DX: $(xcreate_movie_AVS_DX_OBJECTS) $(xcreate_movie_AVS_DX_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} -o $@ $+

#######################################

xcreate_movie_GMT_global_OBJECTS = \
	$O/create_movie_GMT_global.aux.o \
	$(EMPTY_MACRO)

xcreate_movie_GMT_global_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/auto_ner.shared.o \
	$O/count_elements.shared.o \
	$O/count_number_of_sources.shared.o \
	$O/count_points.shared.o \
	$O/define_all_layers.shared.o \
	$O/get_model_parameters.shared.o \
	$O/get_timestep_and_layers.shared.o \
	$O/param_reader.cc.o \
	$O/read_compute_parameters.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$(EMPTY_MACRO)

${E}/xcreate_movie_GMT_global: $(xcreate_movie_GMT_global_OBJECTS) $(xcreate_movie_GMT_global_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} -o $@ $+

#######################################

${E}/xextract_database: $(S_TOP)/utils/extract_database/extract_database.f90 ${OUTPUT}/values_from_mesher.h
	${FCCOMPILE_CHECK} -o ${E}/xextract_database ${FCFLAGS_f90} $(S_TOP)/utils/extract_database/extract_database.f90

#######################################

## compilation directories
S := ${S_TOP}/src/auxiliaries
$(auxiliaries_OBJECTS): S := ${S_TOP}/src/auxiliaries

####
#### rule for each .o file below
####

##
## auxiliaries
##
$O/%.aux.o: $S/%.f90 $O/shared_par.shared_module.o ${OUTPUT}/values_from_mesher.h
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.auxsolver.o: $S/%.f90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.auxmpi.o: $S/%.f90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o $O/parallel.sharedmpi.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.auxsolver.o: $S/%.F90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.auxsolver_vtk.o: $S/%.F90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< $(FC_DEFINE)USE_VTK_INSTEAD_OF_MESH

$O/%.auxadios.o: $S/%.F90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o $O/parallel.sharedmpi.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< $(FC_DEFINE)ADIOS_INPUT

$O/%.auxadios_vtk.o: $S/%.F90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o $O/parallel.sharedmpi.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< $(FC_DEFINE)ADIOS_INPUT $(FC_DEFINE)USE_VTK_INSTEAD_OF_MESH
