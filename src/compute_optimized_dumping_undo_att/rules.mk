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

compute_optimized_dumping_undo_att_TARGETS = \
	$E/xcompute_optimized_dumping_undo_att \
	$(EMPTY_MACRO)

compute_optimized_dumping_undo_att_OBJECTS = \
	$O/compute_optimized_dumping_undo_att.o \
	$(EMPTY_MACRO)

# These files come from the shared directory
compute_optimized_dumping_undo_att_SHARED_OBJECTS = \
	$O/auto_ner.o \
	$O/count_number_of_sources.o \
	$O/euler_angles.o \
	$O/force_ftz.o \
	$O/get_model_parameters.o \
	$O/get_value_parameters.o \
	$O/memory_eval.o \
	$O/param_reader.o \
	$O/read_compute_parameters.o \
	$O/read_parameter_file.o \
	$O/read_value_parameters.o \
	$O/reduce.o \
	$O/rthetaphi_xyz.o \
	$O/save_header_file.o \
	$(EMPTY_MACRO)

#######################################

####
#### rules for executables
####

${E}/xcompute_optimized_dumping_undo_att: $(compute_optimized_dumping_undo_att_OBJECTS) $(compute_optimized_dumping_undo_att_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} -o ${E}/xcompute_optimized_dumping_undo_att $(compute_optimized_dumping_undo_att_OBJECTS) $(compute_optimized_dumping_undo_att_SHARED_OBJECTS)

## uses MPI compiler to link executable instead (usedful for cross-compilation)
#${E}/xcompute_optimized_dumping_undo_att: $(compute_optimized_dumping_undo_att_OBJECTS) $(compute_optimized_dumping_undo_att_SHARED_OBJECTS)
#	${MPIFCCOMPILE_CHECK} -o ${E}/xcompute_optimized_dumping_undo_att $(compute_optimized_dumping_undo_att_OBJECTS) $(compute_optimized_dumping_undo_att_SHARED_OBJECTS)

#######################################

## compilation directories
S := ${S_TOP}/src/compute_optimized_dumping_undo_att
$(compute_optimized_dumping_undo_att_OBJECTS): S := ${S_TOP}/src/compute_optimized_dumping_undo_att

####
#### rule for each .o file below
####

##
## compute_optimized_dumping_undo_att objects
##

$O/compute_optimized_dumping_undo_att.o: $S/compute_optimized_dumping_undo_att.F90
	${FCCOMPILE_CHECK} -c -o $O/compute_optimized_dumping_undo_att.o ${FCFLAGS_f90} $S/compute_optimized_dumping_undo_att.F90

