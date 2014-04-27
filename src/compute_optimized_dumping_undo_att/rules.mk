#=====================================================================
#
#          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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

compute_optimized_dumping_undo_att_TARGETS = \
	$E/xcompute_optimized_dumping_undo_att \
	$(EMPTY_MACRO)

compute_optimized_dumping_undo_att_OBJECTS = \
	$O/compute_optimized_dumping_undo_att.optdump.o \
	$(EMPTY_MACRO)

# These files come from the shared directory
compute_optimized_dumping_undo_att_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/auto_ner.shared.o \
	$O/count_elements.shared.o \
	$O/count_number_of_sources.shared.o \
	$O/count_points.shared.o \
	$O/define_all_layers.shared.o \
	$O/euler_angles.shared.o \
	$O/force_ftz.cc.o \
	$O/get_model_parameters.shared.o \
	$O/get_timestep_and_layers.shared.o \
	$O/get_value_parameters.shared.o \
	$O/memory_eval.shared.o \
	$O/param_reader.cc.o \
	$O/read_compute_parameters.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/save_header_file.shared.o \
	$(EMPTY_MACRO)

#######################################

####
#### rules for executables
####

${E}/xcompute_optimized_dumping_undo_att: $(compute_optimized_dumping_undo_att_SHARED_OBJECTS) $(compute_optimized_dumping_undo_att_OBJECTS)
	${FCCOMPILE_CHECK} -o $@ $+

## uses MPI compiler to link executable instead (useful for cross-compilation)
#${E}/xcompute_optimized_dumping_undo_att: $(compute_optimized_dumping_undo_att_SHARED_OBJECTS) $(compute_optimized_dumping_undo_att_OBJECTS)
#	${MPIFCCOMPILE_CHECK} -o $@ $+

#######################################

## compilation directories
S := ${S_TOP}/src/compute_optimized_dumping_undo_att
$(compute_optimized_dumping_undo_att_OBJECTS): S := ${S_TOP}/src/compute_optimized_dumping_undo_att

####
#### rule for each .o file below
####

$O/%.optdump.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

