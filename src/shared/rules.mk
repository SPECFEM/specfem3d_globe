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

shared_TARGETS = \
	$(shared_OBJECTS) \
	$(EMPTY_MACRO)

shared_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/auto_ner.shared.o \
	$O/broadcast_computed_parameters.shared.o \
	$O/calendar.shared.o \
	$O/count_elements.shared.o \
	$O/count_number_of_sources.shared.o \
	$O/count_points.shared.o \
	$O/create_name_database.shared.o \
	$O/create_serial_name_database.shared.o \
	$O/define_all_layers.shared.o \
	$O/euler_angles.shared.o \
	$O/exit_mpi.shared.o \
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
	$O/parallel.sharedmpi.o \
	$O/param_reader.cc.o \
	$O/read_compute_parameters.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/save_header_file.shared.o \
	$O/spline_routines.shared.o \
	$O/write_c_binary.cc.o \
	$O/write_VTK_file.shared.o \
	$(EMPTY_MACRO)

shared_MODULES = \
 	$(FC_MODDIR)/asdf_var.$(FC_MODEXT) \
	$(FC_MODDIR)/constants.$(FC_MODEXT) \
	$(FC_MODDIR)/shared_input_parameters.$(FC_MODEXT) \
	$(FC_MODDIR)/shared_compute_parameters.$(FC_MODEXT) \
	$(FC_MODDIR)/shared_parameters.$(FC_MODEXT) \
	$(EMPTY_MACRO)

adios_shared_OBJECTS = \
	$O/adios_helpers_definitions.shared_adios_module.o \
	$O/adios_helpers_writers.shared_adios_module.o \
	$O/adios_helpers.shared_adios.o \
	$O/adios_manager.shared_adios.o \
	$O/asdf_helpers_definitions.shared_adios_module.o \
	$O/asdf_helpers_writers.shared_adios_module.o \
	$O/asdf_helpers.shared_adios.o \
	$(EMPTY_MACRO)

adios_shared_MODULES = \
	$(FC_MODDIR)/adios_helpers_definitions_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/adios_helpers_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/adios_helpers_writers_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/asdf_helpers_definitions_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/asdf_helpers_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/asdf_helpers_writer_mod.$(FC_MODEXT) \
	$(EMPTY_MACRO)

adios_shared_STUBS = \
	$O/adios_method_stubs.cc.o \
	$(EMPTY_MACRO)

ifeq ($(ADIOS),yes)
shared_OBJECTS += $(adios_shared_OBJECTS)
shared_MODULES += $(adios_shared_MODULES)
else
shared_OBJECTS += $(adios_shared_STUBS)
endif

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

$O/%.shared_module.o: $S/%.f90 ${SETUP}/constants.h
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared.o: $S/%.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.sharedmpi.o: $S/%.f90 $O/shared_par.shared_module.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

## adios

$O/%.shared_adios_module.o: $S/%.f90
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared_adios.o: $S/%.f90 $O/adios_helpers_writers.shared_adios_module.o $O/adios_helpers_definitions.shared_adios_module.o $O/shared_par.shared_module.o $O/asdf_helpers_writers.shared_adios_module.o $O/asdf_helpers_definitions.shared_adios_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared_adios.o: $S/%.F90 $O/adios_helpers_writers.shared_adios_module.o $O/adios_helpers_definitions.shared_adios_module.o $O/shared_par.shared_module.o $O/asdf_helpers_writers.shared_adios_module.o $O/asdf_helpers_definitions.shared_adios_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<


$O/%.cc.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) -o $@ $< 

