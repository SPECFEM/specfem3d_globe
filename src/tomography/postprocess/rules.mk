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


postprocess_TARGETS = \
	$E/xaddition_sem \
	$E/xclip_sem \
	$E/xcombine_sem \
	$E/xdifference_sem \
	$E/xinterpolate_model \
	$E/xsmooth_sem \
	$(EMPTY_MACRO)
	
ifeq ($(ADIOS),yes)
postprocess_TARGETS += \
	$E/xconvert_model_file_adios \
	$(EMPTY_MACRO)
endif

postprocess_OBJECTS = \
	$(xaddition_sem_OBJECTS) \
	$(xclip_sem_OBJECTS) \
	$(xcombine_sem_OBJECTS) \
	$(xdifference_sem_OBJECTS) \
	$(xinterpolate_model_OBJECTS) \
	$(xsmooth_sem_OBJECTS) \
	$(xconvert_model_file_adios_OBJECTS) \
	$(EMPTY_MACRO)

# These files come from the shared directory
postprocess_SHARED_OBJECTS = \
	$(xaddition_sem_SHARED_OBJECTS) \
	$(xclip_sem_SHARED_OBJECTS) \
	$(xcombine_sem_SHARED_OBJECTS) \
	$(xdifference_sem_SHARED_OBJECTS) \
	$(xinterpolate_model_SHARED_OBJECTS) \
	$(xsmooth_sem_SHARED_OBJECTS) \
	$(xconvert_model_file_adios_SHARED_OBJECTS) \
	$(EMPTY_MACRO)


postprocess_MODULES = \
	$(FC_MODDIR)/postprocess_par.$(FC_MODEXT) \
	$(EMPTY_MACRO)

####
#### rules for executables
####

.PHONY: postprocess

postprocess: $(postprocess_TARGETS)

#######################################

####
#### rules for each program follow
####

#######################################

##
## xconvert_model_file_adios
##
xconvert_model_file_adios_OBJECTS = \
	$O/convert_model_file_adios.postprocess_adios.o \
	$(EMPTY_MACRO)

xconvert_model_file_adios_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/parallel.sharedmpi.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/adios_helpers_definitions.shared_adios_module.o \
	$O/adios_helpers_writers.shared_adios_module.o \
	$O/adios_helpers.shared_adios.o \
	$O/adios_manager.shared_adios.o \
	$(EMPTY_MACRO)

${E}/xconvert_model_file_adios: $(xconvert_model_file_adios_OBJECTS) $(xconvert_model_file_adios_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)

##
## xaddition_sem
##
xaddition_sem_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/addition_sem.postprocess.o \
	$(EMPTY_MACRO)

xaddition_sem_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/parallel.sharedmpi.o \
	$O/exit_mpi.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/specfem3D_par.solverstatic_module.o \
	$(EMPTY_MACRO)

${E}/xaddition_sem: $(xaddition_sem_OBJECTS) $(xaddition_sem_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+

##
## xclip_sem
##
xclip_sem_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/parse_kernel_names.postprocess.o \
	$O/clip_sem.postprocess.o \
	$(EMPTY_MACRO)

xclip_sem_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/parallel.sharedmpi.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/specfem3D_par.solverstatic_module.o \
	$(EMPTY_MACRO)

${E}/xclip_sem: $(xclip_sem_OBJECTS) $(xclip_sem_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)


##
## xcombine_sem
##
xcombine_sem_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/parse_kernel_names.postprocess.o \
	$O/combine_sem.postprocess.o \
	$(EMPTY_MACRO)

xcombine_sem_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/parallel.sharedmpi.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/specfem3D_par.solverstatic_module.o \
	$(EMPTY_MACRO)

${E}/xcombine_sem: $(xcombine_sem_OBJECTS) $(xcombine_sem_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)


##
## xdifference_sem
##
xdifference_sem_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/difference_sem.postprocess.o \
	$(EMPTY_MACRO)

xdifference_sem_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/parallel.sharedmpi.o \
	$O/exit_mpi.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/specfem3D_par.solverstatic_module.o \
	$(EMPTY_MACRO)

${E}/xdifference_sem: $(xdifference_sem_OBJECTS) $(xdifference_sem_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+

##
## xinterpolate_model
##
xinterpolate_model_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/interpolate_model.postprocess.o \
	$(EMPTY_MACRO)

xinterpolate_model_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/parallel.sharedmpi.o \
	$O/gll_library.shared.o \
	$O/heap_sort.shared.o \
	$O/hex_nodes.shared.o \
	$O/lagrange_poly.shared.o \
	$O/recompute_jacobian.solver.o \
	$O/search_kdtree.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/specfem3D_par.solverstatic_module.o \
	$(EMPTY_MACRO)

# extra dependencies
$O/interpolate_model.postprocess.o: $O/search_kdtree.shared.o

${E}/xinterpolate_model: $(xinterpolate_model_OBJECTS) $(xinterpolate_model_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)

##
## xsmooth_sem
##
xsmooth_sem_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/parse_kernel_names.postprocess.o \
	$O/smooth_sem.postprocess.o \
	$(EMPTY_MACRO)

xsmooth_sem_SHARED_OBJECTS = \
	$O/shared_par.shared_module.o \
	$O/parallel.sharedmpi.o \
	$O/exit_mpi.shared.o \
	$O/get_all_eight_slices.shared.o \
	$O/gll_library.shared.o \
	$O/heap_sort.shared.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/search_kdtree.shared.o \
	$O/smooth_weights_vec.shared.o \
	$O/write_VTK_file.shared.o \
	$O/param_reader.cc.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/specfem3D_par.solverstatic_module.o \
	$(EMPTY_MACRO)

# extra dependencies
$O/smooth_sem.postprocess.o: $O/search_kdtree.shared.o

${E}/xsmooth_sem: $(xsmooth_sem_OBJECTS) $(xsmooth_sem_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)


#######################################

## compilation directories
S := ${S_TOP}/src/postprocess
$(postprocess_OBJECTS): S := ${S_TOP}/src/postprocess

###
### module dependencies
###
$O/postprocess_par.postprocess_module.o: $O/shared_par.shared_module.o $O/specfem3D_par.solverstatic_module.o


####
#### rule for each .o file below
####

$O/%.postprocess_module.o: $S/%.f90 ${OUTPUT}/values_from_mesher.h
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.postprocess.o: $S/%.f90 $O/postprocess_par.postprocess_module.o $O/parallel.sharedmpi.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.postprocess.o: $S/%.F90 $O/postprocess_par.postprocess_module.o $O/parallel.sharedmpi.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<


$O/%.postprocess_adios.o: $S/%.F90 $O/postprocess_par.postprocess_module.o $O/parallel.sharedmpi.o $O/adios_helpers.shared_adios.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< $(FC_DEFINE)ADIOS_INPUT

$O/%.postprocess_adios.o: $S/%.f90 $O/postprocess_par.postprocess_module.o $O/parallel.sharedmpi.o $O/adios_helpers.shared_adios.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

