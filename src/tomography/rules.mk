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
	$E/xadd_model_iso \
	$E/xadd_model_tiso \
	$E/xadd_model_tiso_cg \
	$E/xadd_model_tiso_iso \
	$E/xdifference_sem \
	$E/xinterpolate_model \
	$E/xsmooth_sem \
	$E/xsum_kernels \
	$E/xsum_preconditioned_kernels \
	$(EMPTY_MACRO)

ifeq ($(ADIOS),yes)
tomography_TARGETS += \
	$E/xconvert_model_file_adios \
	$(EMPTY_MACRO)
endif

tomography_OBJECTS = \
	$(xadd_model_iso_OBJECTS) \
	$(xadd_model_tiso_OBJECTS) \
	$(xadd_model_tiso_cg_OBJECTS) \
	$(xadd_model_tiso_iso_OBJECTS) \
	$(xdifference_sem_OBJECTS) \
	$(xinterpolate_model_OBJECTS) \
	$(xsmooth_sem_OBJECTS) \
	$(xsum_kernels_OBJECTS) \
	$(xsum_preconditioned_kernels_OBJECTS) \
	$(xconvert_model_file_adios_OBJECTS) \
	$(EMPTY_MACRO)

# These files come from the shared directory
tomography_SHARED_OBJECTS = \
	$(xadd_model_SHARED_OBJECTS) \
	$(xdifference_sem_SHARED_OBJECTS) \
	$(xinterpolate_model_SHARED_OBJECTS) \
	$(xsmooth_sem_SHARED_OBJECTS) \
	$(xsum_kernels_SHARED_OBJECTS) \
	$(xsum_preconditioned_kernels_SHARED_OBJECTS) \
	$(xconvert_model_file_adios_SHARED_OBJECTS) \
	$(EMPTY_MACRO)


tomography_MODULES = \
	$(FC_MODDIR)/tomography_par.$(FC_MODEXT) \
	$(FC_MODDIR)/tomography_kernels_iso.$(FC_MODEXT) \
	$(FC_MODDIR)/tomography_kernels_tiso.$(FC_MODEXT) \
	$(FC_MODDIR)/tomography_kernels_tiso_cg.$(FC_MODEXT) \
	$(FC_MODDIR)/tomography_model_tiso.$(FC_MODEXT) \
	$(FC_MODDIR)/tomography_model_iso.$(FC_MODEXT) \
	$(FC_MODDIR)/kdtree_search.$(FC_MODEXT) \
	$(EMPTY_MACRO)

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
	$O/get_gradient_cg.tomo.o \
	$O/get_gradient_steepest.tomo.o \
	$O/read_kernels.tomo.o \
	$O/read_kernels_cg.tomo.o \
	$O/read_model.tomo.o \
	$O/read_parameters_tomo.tomo.o \
	$O/write_gradients.tomo.o \
	$O/write_new_model.tomo.o \
	$O/write_new_model_perturbations.tomo.o \
	$(EMPTY_MACRO)

xadd_model_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/parallel.sharedmpi.o \
	$O/exit_mpi.shared.o \
	$O/gll_library.shared.o \
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
## xconvert_model_file_adios
##
xconvert_model_file_adios_OBJECTS = \
	$O/convert_model_file_adios.tomoadios.o \
	$(EMPTY_MACRO)

xconvert_model_file_adios_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/parallel.sharedmpi.o \
	$O/adios_helpers_definitions.shared_adios_module.o \
	$O/adios_helpers_writers.shared_adios_module.o \
	$O/adios_helpers.shared_adios.o \
	$O/adios_manager.shared_adios.o \
	$(EMPTY_MACRO)

${E}/xconvert_model_file_adios: $(xconvert_model_file_adios_OBJECTS) $(xconvert_model_file_adios_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)


##
## xdifference_sem
##
xdifference_sem_OBJECTS = \
	$O/difference_sem.tomo.o \
	$(EMPTY_MACRO)

xdifference_sem_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$(EMPTY_MACRO)

${E}/xdifference_sem: $(xdifference_sem_OBJECTS) $(xdifference_sem_SHARED_OBJECTS)
	${FCCOMPILE_CHECK} -o $@ $+

##
## xinterpolate_model
##
xinterpolate_model_OBJECTS = \
	$O/interpolate_model.tomo.o \
	$O/interpolate_model_kdtree.tomo.o \
	$(EMPTY_MACRO)

xinterpolate_model_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
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

##
## xsmooth_sem
##
xsmooth_sem_OBJECTS = \
	$O/smooth_sem.tomo.o \
	$(EMPTY_MACRO)

xsmooth_sem_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/parallel.sharedmpi.o \
	$O/exit_mpi.shared.o \
	$O/get_all_eight_slices.shared.o \
	$O/gll_library.shared.o \
	$O/smooth_weights_vec.shared.o \
	$(EMPTY_MACRO)

${E}/xsmooth_sem: $(xsmooth_sem_OBJECTS) $(xsmooth_sem_SHARED_OBJECTS)
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
	$O/parallel.sharedmpi.o \
	$(EMPTY_MACRO)

${E}/xsum_kernels: $(xsum_kernels_OBJECTS) $(xsum_kernels_SHARED_OBJECTS)
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
	$O/parallel.sharedmpi.o \
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
$O/tomography_par.tomo_module.o: $O/shared_par.shared_module.o


####
#### rule for each .o file below
####

##
## tomography
##

$O/%.tomo_module.o: $S/%.f90 ${SETUP}/constants_tomography.h ${OUTPUT}/values_from_mesher.h $O/shared_par.shared_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<


$O/%.tomo.o: $S/%.f90 $O/tomography_par.tomo_module.o $O/parallel.sharedmpi.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.tomo.o: $S/%.F90 $O/tomography_par.tomo_module.o $O/parallel.sharedmpi.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<


$O/%.tomoadios.o: $S/%.F90 $O/tomography_par.tomo_module.o $O/parallel.sharedmpi.o $O/adios_helpers.shared_adios.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< $(FC_DEFINE)ADIOS_INPUT

$O/%.tomoadios.o: $S/%.f90 $O/tomography_par.tomo_module.o $O/parallel.sharedmpi.o $O/adios_helpers.shared_adios.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

