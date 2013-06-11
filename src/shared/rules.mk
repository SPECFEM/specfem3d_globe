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

shared_TARGETS = \
	$(shared_OBJECTS) \
	$(EMPTY_MACRO)

shared_OBJECTS = \
	$O/auto_ner.o \
	$O/broadcast_compute_parameters.o \
	$O/calendar.o \
	$O/count_number_of_sources.o \
	$O/create_name_database.o \
	$O/create_serial_name_database.o \
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
	$O/write_c_binary.o \
	$(EMPTY_MACRO)

#######################################

## compilation directories
S := ${S_TOP}/src/shared
$(shared_OBJECTS): S = ${S_TOP}/src/shared

####
#### rule for each .o file below
####

##
## shared
##
$O/auto_ner.o: ${SETUP}/constants.h $S/auto_ner.f90
	${FCCOMPILE_CHECK} -c -o $O/auto_ner.o ${FCFLAGS_f90} $S/auto_ner.f90

$O/broadcast_compute_parameters.o: ${SETUP}/constants.h $S/broadcast_compute_parameters.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/broadcast_compute_parameters.o ${FCFLAGS_f90} $S/broadcast_compute_parameters.f90

$O/calendar.o: $S/calendar.f90
	${FCCOMPILE_CHECK} -c -o $O/calendar.o ${FCFLAGS_f90} $S/calendar.f90

$O/count_number_of_sources.o: ${SETUP}/constants.h $S/count_number_of_sources.f90
	${FCCOMPILE_CHECK} -c -o $O/count_number_of_sources.o ${FCFLAGS_f90} $S/count_number_of_sources.f90

$O/create_name_database.o: ${SETUP}/constants.h $S/create_name_database.f90
	${FCCOMPILE_CHECK} -c -o $O/create_name_database.o ${FCFLAGS_f90} $S/create_name_database.f90

$O/create_serial_name_database.o: ${SETUP}/constants.h $S/create_serial_name_database.f90
	${FCCOMPILE_CHECK} -c -o $O/create_serial_name_database.o ${FCFLAGS_f90} $S/create_serial_name_database.f90

$O/euler_angles.o: ${SETUP}/constants.h $S/euler_angles.f90
	${FCCOMPILE_CHECK} -c -o $O/euler_angles.o ${FCFLAGS_f90} $S/euler_angles.f90

### C compilation
$O/force_ftz.o: $S/force_ftz.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) -o $O/force_ftz.o $S/force_ftz.c

$O/get_model_parameters.o: ${SETUP}/constants.h $S/get_model_parameters.f90
	${FCCOMPILE_CHECK} -c -o $O/get_model_parameters.o ${FCFLAGS_f90} $S/get_model_parameters.f90

$O/get_value_parameters.o: ${SETUP}/constants.h $S/get_value_parameters.f90
	${FCCOMPILE_CHECK} -c -o $O/get_value_parameters.o ${FCFLAGS_f90} $S/get_value_parameters.f90

$O/gll_library.o: ${SETUP}/constants.h $S/gll_library.f90
	${FCCOMPILE_CHECK} -c -o $O/gll_library.o ${FCFLAGS_f90} $S/gll_library.f90

$O/hex_nodes.o: ${SETUP}/constants.h $S/hex_nodes.f90
	${FCCOMPILE_CHECK} -c -o $O/hex_nodes.o ${FCFLAGS_f90} $S/hex_nodes.f90

$O/intgrl.o: ${SETUP}/constants.h $S/intgrl.f90
	${FCCOMPILE_CHECK} -c -o $O/intgrl.o ${FCFLAGS_f90} $S/intgrl.f90

$O/lagrange_poly.o: ${SETUP}/constants.h $S/lagrange_poly.f90
	${FCCOMPILE_CHECK} -c -o $O/lagrange_poly.o ${FCFLAGS_f90} $S/lagrange_poly.f90

$O/make_ellipticity.o: ${SETUP}/constants.h $S/make_ellipticity.f90
	${FCCOMPILE_CHECK} -c -o $O/make_ellipticity.o ${FCFLAGS_f90} $S/make_ellipticity.f90

$O/make_gravity.o: ${SETUP}/constants.h $S/make_gravity.f90
	${FCCOMPILE_CHECK} -c -o $O/make_gravity.o ${FCFLAGS_f90} $S/make_gravity.f90

$O/memory_eval.o: ${SETUP}/constants.h $S/memory_eval.f90
	${FCCOMPILE_CHECK} -c -o $O/memory_eval.o ${FCFLAGS_f90} $S/memory_eval.f90

$O/model_prem.o: ${SETUP}/constants.h $S/model_prem.f90
	${FCCOMPILE_CHECK} -c -o $O/model_prem.o ${FCFLAGS_f90} $S/model_prem.f90

$O/param_reader.o: $S/param_reader.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) -o $O/param_reader.o $S/param_reader.c

$O/read_compute_parameters.o: ${SETUP}/constants.h $S/read_compute_parameters.F90
	${FCCOMPILE_CHECK} -c -o $O/read_compute_parameters.o ${FCFLAGS_f90} $S/read_compute_parameters.F90

$O/read_parameter_file.o: ${SETUP}/constants.h $S/read_parameter_file.f90
	${FCCOMPILE_CHECK} -c -o $O/read_parameter_file.o ${FCFLAGS_f90} $S/read_parameter_file.f90

$O/read_value_parameters.o: ${SETUP}/constants.h $S/read_value_parameters.f90
	${FCCOMPILE_CHECK} -c -o $O/read_value_parameters.o ${FCFLAGS_f90} $S/read_value_parameters.f90

$O/reduce.o: ${SETUP}/constants.h $S/reduce.f90
	${FCCOMPILE_CHECK} -c -o $O/reduce.o ${FCFLAGS_f90} $S/reduce.f90

$O/rthetaphi_xyz.o: ${SETUP}/constants.h $S/rthetaphi_xyz.f90
	${FCCOMPILE_CHECK} -c -o $O/rthetaphi_xyz.o ${FCFLAGS_f90} $S/rthetaphi_xyz.f90

$O/save_header_file.o: ${SETUP}/constants.h $S/save_header_file.f90
	${FCCOMPILE_CHECK} -c -o $O/save_header_file.o ${FCFLAGS_f90} $S/save_header_file.f90

$O/spline_routines.o: ${SETUP}/constants.h $S/spline_routines.f90
	${FCCOMPILE_CHECK} -c -o $O/spline_routines.o ${FCFLAGS_f90} $S/spline_routines.f90

### C compilation
$O/write_c_binary.o: $S/write_c_binary.c ${SETUP}/config.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $O/write_c_binary.o $S/write_c_binary.c

##
## shared objects with mpi compilation
##
$O/exit_mpi.o: ${SETUP}/constants.h $S/exit_mpi.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/exit_mpi.o ${FCFLAGS_f90} $S/exit_mpi.f90

$O/model_topo_bathy.o: ${SETUP}/constants.h $S/model_topo_bathy.f90
	${MPIFCCOMPILE_CHECK} -c -o $O/model_topo_bathy.o ${FCFLAGS_f90} $S/model_topo_bathy.f90

