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

tomography_TARGETS = \
	$E/xinterpolate_model \
	$(EMPTY_MACRO)

ifeq ($(ADIOS),yes)
tomography_TARGETS += \
	$E/xconvert_model_file_adios \
	$(EMPTY_MACRO)
endif

tomography_OBJECTS = \
	$(xinterpolate_model_OBJECTS) \
	$(xconvert_model_file_adios_OBJECTS) \
	$(EMPTY_MACRO)

# These files come from the shared directory
tomography_SHARED_OBJECTS = \
	$(xinterpolate_model_SHARED_OBJECTS) \
	$(xconvert_model_file_adios_SHARED_OBJECTS) \
	$(EMPTY_MACRO)


tomography_MODULES = \
	$(FC_MODDIR)/kdtree_search.$(FC_MODEXT) \
	$(EMPTY_MACRO)

####
#### rules for executables
####

.PHONY: all_tomo tomo

all_tomo: $(tomo_TARGETS)

tomo: $(tomo_TARGETS)

#######################################

####
#### rules for each program follow
####

#######################################


## xinterpolate_model

xinterpolate_model_OBJECTS = \
	$O/interpolate_model.tomo.o \
	$O/interpolate_model_kdtree.tomo.o \
	$(EMPTY_MACRO)

xinterpolate_model_SHARED_OBJECTS = \
	$O/parallel.sharedmpi.o \
	$O/gll_library.shared.o \
	$O/hex_nodes.shared.o \
	$O/lagrange_poly.shared.o \
	$O/recompute_jacobian.solver.o \
	$(EMPTY_MACRO)

# extra dependencies
$O/interpolate_model.tomo.o: $O/interpolate_model_kdtree.tomo.o

${E}/xinterpolate_model: $(xinterpolate_model_OBJECTS) $(xinterpolate_model_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)



## xconvert_model_file_adios

xconvert_model_file_adios_OBJECTS = \
	$O/convert_model_file_adios.tomoadios.o \
	$(EMPTY_MACRO)

xconvert_model_file_adios_SHARED_OBJECTS = \
	$O/parallel.sharedmpi.o \
	$O/adios_helpers_definitions.shared_adios_module.o \
	$O/adios_helpers_writers.shared_adios_module.o \
	$O/adios_helpers.shared_adios.o \
	$O/adios_manager.shared_adios.o \
	$(EMPTY_MACRO)

# extra dependencies
#$O/convert_model_file_adios.tomoadios.o: $O/adios_helpers.shared_adios.o

${E}/xconvert_model_file_adios: $(xconvert_model_file_adios_OBJECTS) $(xconvert_model_file_adios_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)


#######################################

## compilation directories
S := ${S_TOP}/src/tomography
$(tomography_OBJECTS): S := ${S_TOP}/src/tomography

####
#### rule for each .o file below
####

##
## tomography
##
$O/%.tomo.o: $S/%.f90 $O/shared_par.shared_module.o ${OUTPUT}/values_from_mesher.h $O/parallel.sharedmpi.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.tomo.o: $S/%.F90 $O/shared_par.shared_module.o ${OUTPUT}/values_from_mesher.h $O/parallel.sharedmpi.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.tomoadios.o: $S/%.F90 $O/shared_par.shared_module.o ${OUTPUT}/values_from_mesher.h $O/parallel.sharedmpi.o $O/adios_helpers.shared_adios.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< $(FC_DEFINE)ADIOS_INPUT

$O/%.tomoadios.o: $S/%.f90 $O/shared_par.shared_module.o ${OUTPUT}/values_from_mesher.h $O/parallel.sharedmpi.o $O/adios_helpers.shared_adios.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

