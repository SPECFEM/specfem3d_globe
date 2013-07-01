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

auxiliaries_TARGETS = \
	$E/xcheck_buffers_1D \
	$E/xcheck_buffers_2D \
	$E/xcheck_buffers_corners_chunks \
	$E/xcheck_buffers_faces_chunks \
	$E/xconvolve_source_timefunction \
	$E/xcombine_AVS_DX \
	$E/xcombine_paraview_strain_data \
	$E/xcombine_vol_data \
	$E/xcombine_surf_data \
	$E/xcreate_movie_AVS_DX \
	$E/xcreate_movie_GMT_global \
	$E/xextract_database \
	$(EMPTY_MACRO)

auxiliaries_OBJECTS = \
	$O/check_buffers_1D.o \
	$O/check_buffers_2D.o \
	$O/check_buffers_corners_chunks.o \
	$O/check_buffers_faces_chunks.o \
	$O/combine_AVS_DX.o \
	$O/combine_paraview_strain_data.o \
	$O/combine_surf_data.o \
	$O/combine_vol_data.o \
	$O/convolve_source_timefunction.o \
	$O/create_movie_AVS_DX.o \
	$O/create_movie_GMT_global.o \
	$(EMPTY_MACRO)

# These files come from the shared directory
auxiliaries_SHARED_OBJECTS = \
	$O/auto_ner.o \
	$O/calendar.o \
	$O/create_serial_name_database.o \
	$O/get_cmt.o \
	$O/get_model_parameters.o \
	$O/get_value_parameters.o \
	$O/param_reader.o \
	$O/read_compute_parameters.o \
	$O/read_parameter_file.o \
	$O/read_value_parameters.o \
	$O/reduce.o \
	$O/rthetaphi_xyz.o \
	$(EMPTY_MACRO)


#######################################

####
#### rules for executables
####

${E}/xcheck_buffers_1D: $O/check_buffers_1D.o $(auxiliaries_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} -o ${E}/xcheck_buffers_1D $O/check_buffers_1D.o $(auxiliaries_SHARED_OBJECTS)

${E}/xcheck_buffers_2D: $O/check_buffers_2D.o $(auxiliaries_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} -o ${E}/xcheck_buffers_2D $O/check_buffers_2D.o $(auxiliaries_SHARED_OBJECTS)

${E}/xcheck_buffers_corners_chunks: $O/check_buffers_corners_chunks.o $(auxiliaries_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} -o ${E}/xcheck_buffers_corners_chunks $O/check_buffers_corners_chunks.o $(auxiliaries_SHARED_OBJECTS)

${E}/xcheck_buffers_faces_chunks: $O/check_buffers_faces_chunks.o $(auxiliaries_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} -o ${E}/xcheck_buffers_faces_chunks $O/check_buffers_faces_chunks.o $(auxiliaries_SHARED_OBJECTS)

${E}/xconvolve_source_timefunction: $O/convolve_source_timefunction.o
	${FCCOMPILE_CHECK} -o ${E}/xconvolve_source_timefunction $O/convolve_source_timefunction.o

${E}/xcombine_AVS_DX: $O/combine_AVS_DX.o $(auxiliaries_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} -o ${E}/xcombine_AVS_DX $O/combine_AVS_DX.o $(auxiliaries_SHARED_OBJECTS)

${E}/xcombine_paraview_strain_data: $O/combine_paraview_strain_data.o $O/write_c_binary.o
	${FCCOMPILE_CHECK} -o ${E}/xcombine_paraview_strain_data  $O/combine_paraview_strain_data.o $O/write_c_binary.o

${E}/xcombine_vol_data: $O/combine_vol_data.o $O/write_c_binary.o
	${FCCOMPILE_CHECK} -o ${E}/xcombine_vol_data  $O/combine_vol_data.o $O/write_c_binary.o

${E}/xcombine_surf_data: $O/combine_surf_data.o $O/write_c_binary.o
	${FCCOMPILE_CHECK} -o ${E}/xcombine_surf_data  $O/combine_surf_data.o $O/write_c_binary.o

${E}/xcreate_movie_AVS_DX: $O/create_movie_AVS_DX.o $(auxiliaries_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} -o ${E}/xcreate_movie_AVS_DX $O/create_movie_AVS_DX.o $(auxiliaries_SHARED_OBJECTS)

${E}/xcreate_movie_GMT_global: $O/create_movie_GMT_global.o $(auxiliaries_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} -o ${E}/xcreate_movie_GMT_global $O/create_movie_GMT_global.o $(auxiliaries_SHARED_OBJECTS)

${E}/xextract_database: $(S_TOP)/utils/extract_database/extract_database.f90
	${FCCOMPILE_CHECK} -o ${E}/xextract_database ${FCFLAGS_f90} $(S_TOP)/utils/extract_database/extract_database.f90

#######################################

## compilation directories
S := ${S_TOP}/src/auxiliaries
$(auxiliaries_OBJECTS): S := ${S_TOP}/src/auxiliaries

####
#### rule for each .o file below
####

##
## auxiliary objects
##

$O/check_buffers_1D.o: ${SETUP}/constants.h $S/check_buffers_1D.f90
	${FCCOMPILE_CHECK} -c -o $O/check_buffers_1D.o ${FCFLAGS_f90} $S/check_buffers_1D.f90

$O/check_buffers_2D.o: ${SETUP}/constants.h $S/check_buffers_2D.f90
	${FCCOMPILE_CHECK} -c -o $O/check_buffers_2D.o ${FCFLAGS_f90} $S/check_buffers_2D.f90

$O/check_buffers_corners_chunks.o: ${SETUP}/constants.h $S/check_buffers_corners_chunks.f90
	${FCCOMPILE_CHECK} -c -o $O/check_buffers_corners_chunks.o ${FCFLAGS_f90} $S/check_buffers_corners_chunks.f90

$O/check_buffers_faces_chunks.o: ${SETUP}/constants.h $S/check_buffers_faces_chunks.f90
	${FCCOMPILE_CHECK} -c -o $O/check_buffers_faces_chunks.o ${FCFLAGS_f90} $S/check_buffers_faces_chunks.f90

$O/combine_AVS_DX.o: ${SETUP}/constants.h $S/combine_AVS_DX.f90
	${FCCOMPILE_CHECK} -c -o $O/combine_AVS_DX.o ${FCFLAGS_f90} $S/combine_AVS_DX.f90

$O/combine_paraview_strain_data.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/combine_paraview_strain_data.f90
	${FCCOMPILE_CHECK} -c -o $O/combine_paraview_strain_data.o ${FCFLAGS_f90} $S/combine_paraview_strain_data.f90

$O/combine_surf_data.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/combine_surf_data.f90
	${FCCOMPILE_CHECK} -c -o $O/combine_surf_data.o ${FCFLAGS_f90} $S/combine_surf_data.f90

$O/combine_vol_data.o: ${SETUP}/constants.h ${OUTPUT}/values_from_mesher.h $S/combine_vol_data.f90
	${FCCOMPILE_CHECK} -c -o $O/combine_vol_data.o ${FCFLAGS_f90} $S/combine_vol_data.f90

$O/convolve_source_timefunction.o: $S/convolve_source_timefunction.f90
	${FCCOMPILE_CHECK} -c -o $O/convolve_source_timefunction.o ${FCFLAGS_f90} $S/convolve_source_timefunction.f90

$O/create_movie_AVS_DX.o: ${SETUP}/constants.h $S/create_movie_AVS_DX.f90
	${FCCOMPILE_CHECK} -c -o $O/create_movie_AVS_DX.o ${FCFLAGS_f90} $S/create_movie_AVS_DX.f90

$O/create_movie_GMT_global.o: ${SETUP}/constants.h $S/create_movie_GMT_global.f90
	${FCCOMPILE_CHECK} -c -o $O/create_movie_GMT_global.o ${FCFLAGS_f90} $S/create_movie_GMT_global.f90

