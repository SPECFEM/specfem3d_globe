#=====================================================================
#
#                       S p e c f e m 3 D  G l o b e
#                       ----------------------------
#
#     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
#                        Princeton University, USA
#                and CNRS / University of Marseille, France
#                 (there are currently many more authors!)
# (c) Princeton University and CNRS / University of Marseille, April 2014
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
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
	$O/adios_manager.shared_adios_module.o \
	$O/assemble_MPI_scalar.shared.o \
	$O/assemble_MPI_vector.shared.o \
	$O/auto_ner.shared.o \
	$O/binary_c_io.cc.o \
	$O/broadcast_computed_parameters.shared.o \
	$O/calendar.shared.o \
	$O/count_elements.shared.o \
	$O/count_number_of_sources.shared.o \
	$O/count_points.shared.o \
	$O/create_addressing.shared.o \
	$O/create_name_database.shared.o \
	$O/create_serial_name_database.shared.o \
	$O/define_all_layers.shared.o \
	$O/euler_angles.shared.o \
	$O/exit_mpi.shared.o \
	$O/fft.shared.o \
	$O/flush_system.shared.o \
	$O/get_all_eight_slices.shared.o \
	$O/get_global.shared.o \
	$O/model_vpremoon.shared.o \
	$O/get_model_parameters.shared.o \
	$O/get_timestep_and_layers.shared.o \
	$O/gll_library.shared.o \
	$O/heap_sort.shared.o \
	$O/hex_nodes.shared.o \
	$O/init_openmp.shared.o \
	$O/interpolate.shared.o \
	$O/intgrl.shared.o \
	$O/lagrange_poly.shared.o \
	$O/make_ellipticity.shared.o \
	$O/memory_eval.shared.o \
	$O/model_mars_1D.shared.o \
	$O/model_prem.shared.o \
	$O/model_Sohl.shared.o \
	$O/model_topo_bathy.shared.o \
	$O/parallel.sharedmpi.o \
	$O/param_reader.cc.o \
	$O/read_compute_parameters.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/recompute_jacobian.shared.o \
	$O/reduce.shared.o \
	$O/rotate_tensor.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/save_header_file.shared.o \
	$O/search_kdtree.shared.o \
	$O/SIEM_math_library.shared.o \
	$O/smooth_weights_vec.shared.o \
	$O/sort_array_coordinates.shared.o \
	$O/spline_routines.shared.o \
	$O/write_VTK_file.shared.o \
	$O/ylm.shared.o \
	$(EMPTY_MACRO)

shared_MODULES = \
	$(FC_MODDIR)/constants.$(FC_MODEXT) \
	$(FC_MODDIR)/manager_adios.$(FC_MODEXT) \
	$(FC_MODDIR)/model_prem_par.$(FC_MODEXT) \
	$(FC_MODDIR)/model_sohl_par.$(FC_MODEXT) \
	$(FC_MODDIR)/model_mars_1d_par.$(FC_MODEXT) \
	$(FC_MODDIR)/model_vpremoon_par.$(FC_MODEXT) \
	$(FC_MODDIR)/my_mpi.$(FC_MODEXT) \
	$(FC_MODDIR)/shared_input_parameters.$(FC_MODEXT) \
	$(FC_MODDIR)/shared_compute_parameters.$(FC_MODEXT) \
	$(FC_MODDIR)/shared_parameters.$(FC_MODEXT) \
	$(FC_MODDIR)/siem_math_library.$(FC_MODEXT) \
	$(FC_MODDIR)/siem_math_library_mpi.$(FC_MODEXT) \
	$(FC_MODDIR)/siem_gll_library.$(FC_MODEXT) \
	$(FC_MODDIR)/kdtree_search.$(FC_MODEXT) \
	$(EMPTY_MACRO)


##
## ADIOS
##

adios_shared_OBJECTS = \
	$O/adios_helpers_addons.shared_adios_cc.o \
	$O/adios_helpers_definitions.shared_adios.o \
	$O/adios_helpers_readers.shared_adios.o \
	$O/adios_helpers_writers.shared_adios.o \
	$O/adios_helpers.shared_adios.o \
	$(EMPTY_MACRO)

adios_shared_MODULES = \
	$(FC_MODDIR)/adios_helpers_definitions_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/adios_helpers_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/adios_helpers_readers_mod.$(FC_MODEXT) \
	$(FC_MODDIR)/adios_helpers_writers_mod.$(FC_MODEXT) \
	$(EMPTY_MACRO)

adios_shared_STUBS = \
	$O/adios_method_stubs.cc.o \
	$(EMPTY_MACRO)

#adios2_shared_OBJECTS = \
#	$O/adios2_helpers_read.shared_adios2_module.o \
#	$(EMPTY_MACRO)

#adios2_shared_MODULES = \
#	$(FC_MODDIR)/adios2_helpers_read_mod.$(FC_MODEXT) \
#	$(EMPTY_MACRO)

ifeq ($(ADIOS),yes)
shared_OBJECTS += $(adios_shared_OBJECTS)
shared_MODULES += $(adios_shared_MODULES)
else ifeq ($(ADIOS2),yes)
shared_OBJECTS += $(adios_shared_OBJECTS)
shared_MODULES += $(adios_shared_MODULES)
else
shared_OBJECTS += $(adios_shared_STUBS)
endif

##
## ASDF
##

asdf_shared_OBJECTS = \
	$O/asdf_manager.shared_asdf.o \
	$(EMPTY_MACRO)

asdf_shared_STUBS = \
	$O/asdf_method_stubs.cc.o \
	$(EMPTY_MACRO)

ifeq ($(ASDF),yes)
shared_OBJECTS += $(asdf_shared_OBJECTS)
else
shared_OBJECTS += $(asdf_shared_STUBS)
endif

##
## C++ Parallel STL sorting
##

ifeq ($(PARALLEL_STL),yes)
shared_OBJECTS += $O/sort_array_coordinates_c.shared.o
endif


#######################################

## compilation directories
S := ${S_TOP}/src/shared
$(shared_OBJECTS): S = ${S_TOP}/src/shared


####
#### rule for each .o file below
####


ifeq ($(ADIOS),yes)
$O/adios_helpers.shared_adios.o: $O/adios_helpers_readers.shared_adios.o $O/adios_helpers_writers.shared_adios.o $O/adios_helpers_definitions.shared_adios.o
else ifeq ($(ADIOS2),yes)
$O/adios_helpers.shared_adios.o: $O/adios_helpers_readers.shared_adios.o $O/adios_helpers_writers.shared_adios.o $O/adios_helpers_definitions.shared_adios.o
endif

$O/get_model_parameters.shared.o: \
	$O/model_prem.shared.o \
	$O/model_Sohl.shared.o \
	$O/model_vpremoon.shared.o \
	$O/model_mars_1D.shared.o

$O/make_ellipticity.shared.o: \
	$O/model_prem.shared.o \
	$O/model_Sohl.shared.o \
	$O/model_vpremoon.shared.o

##
## shared
##

$O/read_parameter_file.shared.o: $S/read_parameter_file.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< $(ADIOS_DEF) $(ADIOS2_DEF)

$O/%.shared_module.o: $S/%.f90 ${SETUP}/constants.h
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared.o: $S/%.f90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared.o: $S/%.F90 $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.sharedmpi.o: $S/%.f90 $O/shared_par.shared_module.o
	${MPIFCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

## adios

$O/%.shared_adios_module.o: $S/%.f90 $O/shared_par.shared_module.o $O/parallel.sharedmpi.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared_adios_module.o: $S/%.F90 $O/shared_par.shared_module.o $O/parallel.sharedmpi.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared_adios.o: $S/%.f90 $O/adios_manager.shared_adios_module.o $S/adios_helpers_definitions_1d_generic.inc $S/adios_helpers_readers_1d_generic.inc $S/adios_helpers_writers_1d_generic.inc $S/adios_helpers_writers_1d_generic_offset.inc
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared_adios.o: $S/%.F90 $O/adios_manager.shared_adios_module.o $S/adios_helpers_definitions_1d_generic.inc $S/adios_helpers_readers_1d_generic.inc $S/adios_helpers_writers_1d_generic.inc $S/adios_helpers_writers_1d_generic_offset.inc
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.shared_adios_cc.o: $S/%.c ${SETUP}/config.h
	${MPICC} -c $(CPPFLAGS) $(CFLAGS) $(ADIOS2_DEF) $(ADIOS_DEF) -o $@ $<


## ADIOS2
# special two step build process for adios2_helpers_read
# 1. preprocess and then replace all _NL_ string with new lines (-> $O/adios2_helpers_read.f90)
# 2. compile generated f90
#$S/adios2_helpers_read_generated.f90: $S/adios2_helpers_read_template.F90
#	@echo " == Rule Read-generate triggered for $< "
#	@echo " If this auto-generation does not work, you can do it manually using gfortran (version 5 or higher)"
#	@echo "     $$ gfortran -cpp -E -o $S/adios2_helpers_read_generated.f90 $S/adios2_helpers_read_template.F90"
#	@echo "     $$ sed -i " '"s\/ *_NL_\/\\n\/g"'  " adios2_helpers_read_generated.f90"
#	gfortran -cpp -E -o $S/adios2_helpers_read_generated.f90 $S/adios2_helpers_read_template.F90
#	sed -i -e "s/ *_NL_/\n/g" -e "s/##//g" $S/adios2_helpers_read_generated.f90


#$O/adios2_helpers_read.shared_adios2_module.o: $S/adios2_helpers_read_generated.f90 $O/shared_par.shared_module.o
#	@echo " == Rule Read-compile triggered for $< "
#	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $S/adios2_helpers_read_generated.f90 $(ADIOS2_DEF)

#$O/%.shared_adios2_module.o: $S/%.f90 $O/shared_par.shared_module.o
#	@echo " == Rule 1 triggered for $< "
#	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< $(ADIOS2_DEF)

#$O/%.shared_adios2_module.o: $S/%.F90 $O/shared_par.shared_module.o
#	@echo " == Rule 2 triggered for $< "
#	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< $(ADIOS2_DEF)

#$O/%.shared_adios2.o: $S/%.F90 $O/adios2_manager.shared_adios2_module.o  $O/adios2_helpers_read.shared_adios_module.o  $O/shared_par.shared_module.o
#	@echo " == Rule 3 triggered for $< "
#	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< $(ADIOS2_DEF)

## asdf

$O/%.shared_asdf.o: $S/%.f90
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.cc.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) -o $@ $<

## c++ files
$O/%.shared.o: $S/%.cpp ${SETUP}/config.h
	${CXX} -c $(CXXFLAGS) $(PARALLEL_STL_DEF) -o $@ $<

