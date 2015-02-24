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

create_header_file_TARGETS = \
	$E/xcreate_header_file \
	${OUTPUT}/values_from_mesher.h \
	$(EMPTY_MACRO)

create_header_file_OBJECTS = \
	$O/create_header_file.header.o \
	$(EMPTY_MACRO)

# These files come from the shared directory
create_header_file_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/auto_ner.shared.o \
	$O/count_elements.shared.o \
	$O/count_number_of_sources.shared.o \
	$O/count_points.shared.o \
	$O/define_all_layers.shared.o \
	$O/euler_angles.shared.o \
	$O/get_model_parameters.shared.o \
	$O/get_timestep_and_layers.shared.o \
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
#### Mesher header file
####

${OUTPUT}/values_from_mesher.h: $E/xcreate_header_file $B/DATA/Par_file
	@-rm -f $@
	@echo ""
	@echo "running xcreate_header_file..."
	@echo ""
	$E/xcreate_header_file
	@echo ""
	@test -f $@

#######################################

####
#### rules for executables
####

${E}/xcreate_header_file: $(create_header_file_SHARED_OBJECTS) $(create_header_file_OBJECTS)
	${FCCOMPILE_CHECK} -o $@ $+

## uses MPI compiler to link executable instead (useful for cross-compilation)
#${E}/xcreate_header_file: $(create_header_file_SHARED_OBJECTS) $(create_header_file_OBJECTS)
#	${MPIFCCOMPILE_CHECK} -o $@ $+

#######################################

## compilation directories
S := ${S_TOP}/src/create_header_file
$(create_header_file_OBJECTS): S := ${S_TOP}/src/create_header_file

####
#### rule for each .o file below
####

##
## rule for create_header_file
##
$O/%.header.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

