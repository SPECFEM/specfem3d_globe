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

tomography_TARGETS = \
	$E/xadd_model_iso \
	$E/xadd_model_tiso \
	$E/xadd_model_tiso_cg \
	$E/xadd_model_tiso_iso \
	$E/xsum_kernels \
	$E/xsum_preconditioned_kernels \
	$(EMPTY_MACRO)


tomography_OBJECTS = \
	$(xadd_model_iso_OBJECTS) \
	$(xadd_model_tiso_OBJECTS) \
	$(xadd_model_tiso_cg_OBJECTS) \
	$(xadd_model_tiso_iso_OBJECTS) \
	$(xsum_kernels_OBJECTS) \
	$(xsum_preconditioned_kernels_OBJECTS) \
	$(EMPTY_MACRO)


# These files come from the shared directory
tomography_SHARED_OBJECTS = \
	$(xadd_model_SHARED_OBJECTS) \
	$(xsum_kernels_SHARED_OBJECTS) \
	$(xsum_preconditioned_kernels_SHARED_OBJECTS) \
	$(EMPTY_MACRO)


tomography_MODULES = \
	$(FC_MODDIR)/tomography_par.$(FC_MODEXT) \
	$(FC_MODDIR)/tomography_kernels_iso.$(FC_MODEXT) \
	$(FC_MODDIR)/tomography_kernels_tiso.$(FC_MODEXT) \
	$(FC_MODDIR)/tomography_kernels_tiso_cg.$(FC_MODEXT) \
	$(FC_MODDIR)/tomography_model_tiso.$(FC_MODEXT) \
	$(FC_MODDIR)/tomography_model_iso.$(FC_MODEXT) \
	$(EMPTY_MACRO)


####
#### ADIOS versions
####

adios_tomography_TARGETS += \
	$E/xsum_kernels_adios \
	$(EMPTY_MACRO)

adios_tomography_OBJECTS += \
	$(xsum_kernels_adios_OBJECTS) \
	$(EMPTY_MACRO)

ifeq ($(ADIOS),yes)
tomography_TARGETS += $(adios_tomography_TARGETS)
tomography_OBJECTS += $(adios_tomography_OBJECTS)
else ifeq ($(ADIOS2),yes)
tomography_TARGETS += $(adios_tomography_TARGETS)
tomography_OBJECTS += $(adios_tomography_OBJECTS)
endif


####
#### rules for executables
####

.PHONY: all_tomo tomo tomography

all_tomo: $(tomography_TARGETS)

tomo: $(tomography_TARGETS)

tomography: $(tomography_TARGETS)

#######################################

####
#### rules for each program follow
####

#######################################

xadd_model_OBJECTS = \
	$O/tomography_par.tomo_module.o \
	$O/compute_kernel_integral.tomo.o \
	$O/get_cg_direction.tomo.o \
	$O/get_sd_direction.tomo.o \
	$O/read_kernels.tomo.o \
	$O/read_kernels_cg.tomo.o \
	$O/read_model.tomo.o \
	$O/read_parameters_tomo.tomo.o \
	$O/write_gradient.tomo.o \
	$O/write_new_model.tomo.o \
	$O/write_new_model_perturbations.tomo.o \
	$(EMPTY_MACRO)

xadd_model_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/specfem3D_par.solverstatic_module.o \
	$O/read_mesh_parameters.solverstatic.o \
	$O/parallel.sharedmpi.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/gll_library.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

##
## xadd_model_iso
##
xadd_model_iso_OBJECTS = \
	$O/add_model_iso.tomo.o \
	$(xadd_model_OBJECTS) \
	$(EMPTY_MACRO)

${E}/xadd_model_iso: $(xadd_model_iso_OBJECTS) $(xadd_model_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)


##
## xadd_model_tiso
##
xadd_model_tiso_OBJECTS = \
	$O/add_model_tiso.tomo.o \
	$(xadd_model_OBJECTS) \
	$(EMPTY_MACRO)

${E}/xadd_model_tiso: $(xadd_model_tiso_OBJECTS) $(xadd_model_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)

##
## xadd_model_tiso_cg
##
xadd_model_tiso_cg_OBJECTS = \
	$O/add_model_tiso_cg.tomo.o \
	$(xadd_model_OBJECTS) \
	$(EMPTY_MACRO)

${E}/xadd_model_tiso_cg: $(xadd_model_tiso_cg_OBJECTS) $(xadd_model_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)


##
## xadd_model_tiso_iso
##
xadd_model_tiso_iso_OBJECTS = \
	$O/add_model_tiso_iso.tomo.o \
	$(xadd_model_OBJECTS) \
	$(EMPTY_MACRO)

${E}/xadd_model_tiso_iso: $(xadd_model_tiso_iso_OBJECTS) $(xadd_model_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)


##
## xsum_kernels
##
xsum_kernels_OBJECTS = \
	$O/tomography_par.tomo_module.o \
	$O/sum_kernels.tomo.o \
	$(EMPTY_MACRO)

xsum_kernels_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/specfem3D_par.solverstatic_module.o \
	$O/read_mesh_parameters.solverstatic.o \
	$O/parallel.sharedmpi.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

${E}/xsum_kernels: $(xsum_kernels_OBJECTS) $(xsum_kernels_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)

##
## xsum_kernels_adios
##
xsum_kernels_adios_OBJECTS = \
	$O/tomography_par.tomo_module.o \
	$O/sum_kernels.tomo_adios.o \
	$(EMPTY_MACRO)

xsum_kernels_adios_SHARED_OBJECTS = \
	$(xsum_kernels_SHARED_OBJECTS) \
	$O/adios_helpers_addons.shared_adios_cc.o \
	$O/adios_helpers_definitions.shared_adios.o \
	$O/adios_helpers_readers.shared_adios.o \
	$O/adios_helpers_writers.shared_adios.o \
	$O/adios_helpers.shared_adios.o \
	$O/adios_manager.shared_adios_module.o \
	$(EMPTY_MACRO)

# extra dependencies
$O/sum_kernels.tomo_adios.o: $O/adios_manager.shared_adios_module.o

${E}/xsum_kernels_adios: $(xsum_kernels_adios_OBJECTS) $(xsum_kernels_adios_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)


##
## xsum_preconditioned_kernels
##
xsum_preconditioned_kernels_OBJECTS = \
	$O/tomography_par.tomo_module.o \
	$O/sum_preconditioned_kernels.tomo.o \
	$(EMPTY_MACRO)

xsum_preconditioned_kernels_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/specfem3D_par.solverstatic_module.o \
	$O/read_mesh_parameters.solverstatic.o \
	$O/parallel.sharedmpi.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

${E}/xsum_preconditioned_kernels: $(xsum_preconditioned_kernels_OBJECTS) $(xsum_preconditioned_kernels_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)

#######################################

## compilation directories
S := ${S_TOP}/src/tomography
$(tomography_OBJECTS): S := ${S_TOP}/src/tomography

###
### Model dependencies
###
$O/tomography_par.tomo_module.o: $O/shared_par.shared_module.o $O/specfem3D_par.solverstatic_module.o


####
#### rule for each .o file below
####

##
## tomography
##

$O/%.tomo_module.o: $S/%.f90 ${SETUP}/constants_tomography.h $O/shared_par.shared_module.o $O/specfem3D_par.solverstatic_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<


$O/%.tomo.o: $S/%.f90 $O/tomography_par.tomo_module.o $O/parallel.sharedmpi.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.tomo.o: $S/%.F90 $O/tomography_par.tomo_module.o $O/parallel.sharedmpi.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

### ADIOS
$O/%.tomo_adios.o: $S/%.f90 $O/tomography_par.tomo_module.o $O/parallel.sharedmpi.o $O/adios_helpers.shared_adios.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.tomo_adios.o: $S/%.F90 $O/tomography_par.tomo_module.o $O/parallel.sharedmpi.o $O/adios_helpers.shared_adios.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< $(FC_DEFINE)USE_ADIOS_INSTEAD_OF_MESH


