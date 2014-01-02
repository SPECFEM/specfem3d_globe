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

auxiliaries_TARGETS = \
	$E/xconvolve_source_timefunction \
	$E/xcombine_AVS_DX \
	$E/xcombine_paraview_strain_data \
	$E/xcombine_vol_data \
	$E/xcombine_vol_data_vtk \
	$E/xcombine_surf_data \
	$E/xcreate_movie_AVS_DX \
	$E/xcreate_movie_GMT_global \
	$E/xextract_database \
	$(EMPTY_MACRO)

auxiliaries_OBJECTS = \
	$O/combine_AVS_DX.aux.o \
	$O/combine_paraview_strain_data.auxsolver.o \
	$O/combine_surf_data.auxsolver.o \
	$O/combine_vol_data.auxsolver.o \
	$O/combine_vol_data_vtk.auxsolver.o \
	$O/convolve_source_timefunction.aux.o \
	$O/create_movie_AVS_DX.aux.o \
	$O/create_movie_GMT_global.aux.o \
	$(EMPTY_MACRO)

# These files come from the shared directory
auxiliaries_SHARED_OBJECTS = \
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
	$O/get_value_parameters.shared.o \
	$O/param_reader.cc.o \
	$O/read_compute_parameters.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$(EMPTY_MACRO)


#######################################

####
#### rules for executables
####

all_aux: required $(auxiliaries_TARGETS)

aux: required $(auxiliaries_TARGETS)

${E}/xconvolve_source_timefunction: $O/convolve_source_timefunction.aux.o
	${FCCOMPILE_CHECK} -o ${E}/xconvolve_source_timefunction $O/convolve_source_timefunction.aux.o

${E}/xcombine_AVS_DX: $(auxiliaries_SHARED_OBJECTS) $O/get_cmt.solver.o $O/combine_AVS_DX.aux.o
	${FCCOMPILE_CHECK} -o ${E}/xcombine_AVS_DX $(auxiliaries_SHARED_OBJECTS) $O/get_cmt.solver.o $O/combine_AVS_DX.aux.o

${E}/xcombine_paraview_strain_data: $(auxiliaries_SHARED_OBJECTS) $O/combine_paraview_strain_data.auxsolver.o $O/write_c_binary.cc.o
	${FCCOMPILE_CHECK} -o ${E}/xcombine_paraview_strain_data $(auxiliaries_SHARED_OBJECTS) $O/combine_paraview_strain_data.auxsolver.o $O/write_c_binary.cc.o

${E}/xcombine_vol_data: $(auxiliaries_SHARED_OBJECTS) $O/combine_vol_data.auxsolver.o $O/write_c_binary.cc.o $O/combine_vol_data_shared.aux.o
	${FCCOMPILE_CHECK} -o ${E}/xcombine_vol_data $(auxiliaries_SHARED_OBJECTS) $O/combine_vol_data.auxsolver.o $O/write_c_binary.cc.o $O/combine_vol_data_shared.aux.o

${E}/xcombine_vol_data_vtk: $(auxiliaries_SHARED_OBJECTS) $O/combine_vol_data.auxsolver_vtk.o $O/write_c_binary.cc.o $O/combine_vol_data_shared.aux.o
	${FCCOMPILE_CHECK} -o ${E}/xcombine_vol_data_vtk $(auxiliaries_SHARED_OBJECTS) $O/combine_vol_data.auxsolver_vtk.o $O/write_c_binary.cc.o $O/combine_vol_data_shared.aux.o

${E}/xcombine_surf_data: $(auxiliaries_SHARED_OBJECTS) $O/combine_surf_data.auxsolver.o $O/write_c_binary.cc.o
	${FCCOMPILE_CHECK} -o ${E}/xcombine_surf_data $(auxiliaries_SHARED_OBJECTS) $O/combine_surf_data.auxsolver.o $O/write_c_binary.cc.o

${E}/xcreate_movie_AVS_DX: $(auxiliaries_SHARED_OBJECTS) $O/create_movie_AVS_DX.aux.o
	${FCCOMPILE_CHECK} -o ${E}/xcreate_movie_AVS_DX $(auxiliaries_SHARED_OBJECTS) $O/create_movie_AVS_DX.aux.o

${E}/xcreate_movie_GMT_global: $(auxiliaries_SHARED_OBJECTS) $O/create_movie_GMT_global.aux.o
	${FCCOMPILE_CHECK} -o ${E}/xcreate_movie_GMT_global $(auxiliaries_SHARED_OBJECTS) $O/create_movie_GMT_global.aux.o

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
$O/%.aux.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.auxsolver.o: $S/%.f90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.auxsolver.o: $S/%.F90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.auxsolver_vtk.o: $S/%.F90 ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< -DUSE_VTK_INSTEAD_OF_MESH
