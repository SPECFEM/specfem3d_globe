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

## compilation directories
S := ${S_TOP}/src/tomography/postprocess_sensitivity_kernels
$(tomography/postprocess_sensitivity_kernels_OBJECTS): S := ${S_TOP}/src/tomography/postprocess_sensitivity_kernels

#######################################################################

tomography/postprocess_sensitivity_kernels_TARGETS = \
	$E/xaddition_sem \
	$E/xclip_sem \
	$E/xcombine_sem \
	$E/xdifference_sem \
	$E/xinterpolate_model \
	$E/xcreate_cross_section \
	$E/xsmooth_sem \
	$E/xsmooth_laplacian_sem \
	$(EMPTY_MACRO)

tomography/postprocess_sensitivity_kernels_OBJECTS = \
	$(xaddition_sem_OBJECTS) \
	$(xclip_sem_OBJECTS) \
	$(xcombine_sem_OBJECTS) \
	$(xdifference_sem_OBJECTS) \
	$(xinterpolate_model_OBJECTS) \
	$(xcreate_cross_section_OBJECTS) \
	$(xsmooth_sem_OBJECTS) \
	$(xsmooth_laplacian_sem_OBJECTS) \
	$(xconvert_model_file_adios_OBJECTS) \
	$(EMPTY_MACRO)


## ADIOS
tomography/adios_postprocess_sensitivity_kernels_TARGETS += \
	$E/xconvert_model_file_adios \
	$E/xinterpolate_model_adios \
	$E/xsmooth_sem_adios \
	$E/xsmooth_laplacian_sem_adios \
	$(EMPTY_MACRO)

tomography/adios_postprocess_sensitivity_kernels_OBJECTS += \
	$(xinterpolate_model_adios_OBJECTS) \
	$(xsmooth_sem_adios_OBJECTS) \
	$(xsmooth_laplacian_sem_adios_OBJECTS) \
	$(EMPTY_MACRO)

ifeq ($(ADIOS),yes)
tomography/postprocess_sensitivity_kernels_TARGETS += $(tomography/adios_postprocess_sensitivity_kernels_TARGETS)
tomography/postprocess_sensitivity_kernels_OBJECTS += $(tomography/adios_postprocess_sensitivity_kernels_OBJECTS)
else ifeq ($(ADIOS2),yes)
tomography/postprocess_sensitivity_kernels_TARGETS += $(tomography/adios_postprocess_sensitivity_kernels_TARGETS)
tomography/postprocess_sensitivity_kernels_OBJECTS += $(tomography/adios_postprocess_sensitivity_kernels_OBJECTS)
endif

# These files come from the shared directory
tomography/postprocess_sensitivity_kernels_SHARED_OBJECTS = \
	$(xaddition_sem_SHARED_OBJECTS) \
	$(xclip_sem_SHARED_OBJECTS) \
	$(xcombine_sem_SHARED_OBJECTS) \
	$(xdifference_sem_SHARED_OBJECTS) \
	$(xinterpolate_model_SHARED_OBJECTS) \
	$(xcreate_cross_section_SHARED_OBJECTS) \
	$(xsmooth_sem_SHARED_OBJECTS) \
	$(xsmooth_laplacian_sem_SHARED_OBJECTS) \
	$(xconvert_model_file_adios_SHARED_OBJECTS) \
	$(EMPTY_MACRO)


tomography/postprocess_sensitivity_kernels_MODULES = \
	$(FC_MODDIR)/postprocess_par.$(FC_MODEXT) \
	$(EMPTY_MACRO)

####
#### rules for executables
####

.PHONY: postprocess

postprocess: $(tomography/postprocess_sensitivity_kernels_TARGETS)

postprocess_sensitivity_kernels: postprocess

tomography/postprocess_sensitivity_kernels: postprocess

#######################################

####
#### rules for each program follow
####

#######################################

##
## xconvert_model_file_adios
##
xconvert_model_file_adios_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/convert_model_file_adios.postprocess_adios.o \
	$(EMPTY_MACRO)

xconvert_model_file_adios_SHARED_OBJECTS = \
	$O/specfem3D_par.solverstatic_module.o \
	$O/read_mesh_parameters.solverstatic.o \
	$O/shared_par.shared_module.o \
	$O/parallel.sharedmpi.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/param_reader.cc.o \
	$O/read_value_parameters.shared.o \
	$O/adios_helpers_addons.shared_adios_cc.o \
	$O/adios_helpers_definitions.shared_adios.o \
	$O/adios_helpers_readers.shared_adios.o \
	$O/adios_helpers_writers.shared_adios.o \
	$O/adios_helpers.shared_adios.o \
	$O/adios_manager.shared_adios_module.o \
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
	$O/specfem3D_par.solverstatic_module.o \
	$O/read_mesh_parameters.solverstatic.o \
	$O/shared_par.shared_module.o \
	$O/parallel.sharedmpi.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/param_reader.cc.o \
	$O/read_value_parameters.shared.o \
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
	$O/specfem3D_par.solverstatic_module.o \
	$O/read_mesh_parameters.solverstatic.o \
	$O/shared_par.shared_module.o \
	$O/parallel.sharedmpi.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/param_reader.cc.o \
	$O/read_value_parameters.shared.o \
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
	$O/specfem3D_par.solverstatic_module.o \
	$O/read_mesh_parameters.solverstatic.o \
	$O/shared_par.shared_module.o \
	$O/parallel.sharedmpi.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/param_reader.cc.o \
	$O/read_value_parameters.shared.o \
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
	$O/specfem3D_par.solverstatic_module.o \
	$O/read_mesh_parameters.solverstatic.o \
	$O/shared_par.shared_module.o \
	$O/parallel.sharedmpi.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/param_reader.cc.o \
	$O/read_value_parameters.shared.o \
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
	$O/specfem3D_par.solverstatic_module.o \
	$O/read_mesh_parameters.solverstatic.o \
	$O/shared_par.shared_module.o \
	$O/parallel.sharedmpi.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/gll_library.shared.o \
	$O/heap_sort.shared.o \
	$O/hex_nodes.shared.o \
	$O/interpolate.shared.o \
	$O/lagrange_poly.shared.o \
	$O/recompute_jacobian.shared.o \
	$O/search_kdtree.shared.o \
	$O/param_reader.cc.o \
	$O/read_value_parameters.shared.o \
	$(EMPTY_MACRO)

# extra dependencies
$O/interpolate_model.postprocess.o: $O/search_kdtree.shared.o


${E}/xinterpolate_model: $(xinterpolate_model_OBJECTS) $(xinterpolate_model_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)

##
## xinterpolate_model_adios
##
xinterpolate_model_adios_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/interpolate_model.postprocess_adios.o \
	$(EMPTY_MACRO)

xinterpolate_model_adios_SHARED_OBJECTS = \
	$(xinterpolate_model_SHARED_OBJECTS)

xinterpolate_model_adios_SHARED_OBJECTS += \
	$O/adios_helpers_addons.shared_adios_cc.o \
	$O/adios_helpers_definitions.shared_adios.o \
	$O/adios_helpers_readers.shared_adios.o \
	$O/adios_helpers_writers.shared_adios.o \
	$O/adios_helpers.shared_adios.o \
	$O/adios_manager.shared_adios_module.o \
	$(EMPTY_MACRO)

# extra dependencies
$O/interpolate_model.postprocess_adios.o: $O/search_kdtree.shared.o $O/adios_manager.shared_adios_module.o

${E}/xinterpolate_model_adios: $(xinterpolate_model_adios_OBJECTS) $(xinterpolate_model_adios_SHARED_OBJECTS)
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
	$O/specfem3D_par.solverstatic_module.o \
	$O/read_mesh_parameters.solverstatic.o \
	$O/shared_par.shared_module.o \
	$O/parallel.sharedmpi.o \
	$O/auto_ner.shared.o \
	$O/count_elements.shared.o \
	$O/count_points.shared.o \
	$O/define_all_layers.shared.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/get_all_eight_slices.shared.o \
	$O/get_model_parameters.shared.o \
	$O/get_timestep_and_layers.shared.o \
	$O/gll_library.shared.o \
	$O/heap_sort.shared.o \
	$O/model_mars_1D.shared.o \
	$O/model_vpremoon.shared.o \
	$O/param_reader.cc.o \
	$O/read_compute_parameters.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/search_kdtree.shared.o \
	$O/smooth_weights_vec.shared.o \
	$O/write_VTK_file.shared.o \
	$(EMPTY_MACRO)

###
### GPU
###

gpu_smooth_OBJECTS = \
	$O/check_fields_gpu.o \
	$O/helper_functions_gpu.o \
	$O/initialize_gpu.o \
	$O/smooth_gpu.o \
	$O/transfer_fields_gpu.o \
	$(EMPTY_MACRO)

gpu_smooth_STUBS = \
	$O/specfem3D_gpu_method_stubs.gpu_cc.o \
	$(EMPTY_MACRO)

# kernel files
gpu_smooth_kernels_OBJS := \
	$O/get_maximum_scalar_kernel.o \
	$O/get_maximum_vector_kernel.o \
	$O/smooth_normalize_data_kernel.o \
	$O/smooth_process_kernel.o \
	$(EMPTY_MACRO)

## CUDA
ifeq ($(CUDA),yes)
	# cuda
  # renames endings
	cuda_smooth_kernels_OBJS:=$(subst .o,.cuda-kernel.o,${gpu_smooth_kernels_OBJS})
	cuda_smooth_DEVICE_OBJ =  $O/cuda_smooth_device_obj.o
endif

## HIP
ifeq ($(HIP),yes)
  # defines $(cuda_smooth_kernels_OBJS)
  # renames endings
  cuda_smooth_kernels_OBJS:=$(subst .o,.hip-kernel.o,${gpu_smooth_kernels_OBJS})
endif

ifdef NO_GPU
	gpu_xs_OBJECTS = $(gpu_smooth_STUBS)
else
	gpu_xs_OBJECTS = $(gpu_smooth_OBJECTS)
endif

# substitutes object endings to assign corresponding compilation rule
ifeq ($(GPU_CUDA_AND_OCL),yes)
	# combines both CUDA and OpenCL kernels compilation
  gpu_smooth_OBJECTS:=$(subst .o,.cuda-ocl.o,${gpu_smooth_OBJECTS})
endif

ifneq ($(GPU_CUDA_AND_OCL),yes)
  # OpenCL kernels only
  ifeq ($(OCL), yes)
    gpu_smooth_OBJECTS:=$(subst .o,.ocl.o,${gpu_smooth_OBJECTS})
  endif

  # CUDA kernels only
  ifeq ($(CUDA),yes)
    gpu_smooth_OBJECTS:=$(subst .o,.cuda.o,${gpu_smooth_OBJECTS})
  endif

  # HIP kernels only
  ifeq ($(HIP), yes)
    gpu_smooth_OBJECTS:=$(subst .o,.hip.o,${gpu_smooth_OBJECTS})
  endif
endif

gpu_smooth_OBJECTS += $(cuda_smooth_DEVICE_OBJ) $(cuda_smooth_kernels_OBJS)

##
## compilation
##
xsmooth_sem_SHARED_OBJECTS += $(gpu_xs_OBJECTS)
xsmooth_sem_LIBS = $(MPILIBS)  # $(LDFLAGS) $(MPILIBS) $(LIBS)
xsmooth_sem_LIBS += $(GPU_LINK)

INFO_SMOOTH="building xsmooth_sem $(BUILD_VERSION_TXT)"

# extra dependencies
$O/smooth_sem.postprocess.o: $O/search_kdtree.shared.o

${E}/xsmooth_sem: $(xsmooth_sem_OBJECTS) $(xsmooth_sem_SHARED_OBJECTS) # $(COND_MPI_OBJECTS)
	@echo ""
	@echo $(INFO_SMOOTH)
	@echo ""
	${FCLINK} -o $@ $+ $(xsmooth_sem_LIBS)
	@echo ""

$(cuda_smooth_DEVICE_OBJ): $(subst $(cuda_smooth_DEVICE_OBJ), ,$(gpu_smooth_OBJECTS)) $(cuda_smooth_kernels_OBJS)
	${NVCCLINK} -o $@ $^


##
## xsmooth_sem_adios
##
xsmooth_sem_adios_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/parse_kernel_names.postprocess.o \
	$O/smooth_sem.postprocess_adios.o \
	$(EMPTY_MACRO)

xsmooth_sem_adios_SHARED_OBJECTS = \
	$(xsmooth_sem_SHARED_OBJECTS)

xsmooth_sem_adios_SHARED_OBJECTS += \
	$O/adios_helpers_addons.shared_adios_cc.o \
	$O/adios_helpers_definitions.shared_adios.o \
	$O/adios_helpers_readers.shared_adios.o \
	$O/adios_helpers_writers.shared_adios.o \
	$O/adios_helpers.shared_adios.o \
	$O/adios_manager.shared_adios_module.o \
	$(EMPTY_MACRO)

##
## compilation
##
xsmooth_sem_adios_LIBS = $(MPILIBS)  # $(LDFLAGS) $(MPILIBS) $(LIBS)
xsmooth_sem_adios_LIBS += $(GPU_LINK)

INFO_SMOOTH_ADIOS="building xsmooth_sem_adios $(BUILD_VERSION_TXT)"

# extra dependencies
$O/smooth_sem.postprocess_adios.o: $O/search_kdtree.shared.o

${E}/xsmooth_sem_adios: $(xsmooth_sem_adios_OBJECTS) $(xsmooth_sem_adios_SHARED_OBJECTS)
	@echo ""
	@echo $(INFO_SMOOTH_ADIOS)
	@echo ""
	${FCLINK} -o $@ $+ $(xsmooth_sem_adios_LIBS)
	@echo ""



##
## xsmooth_laplacian_sem
##
xsmooth_laplacian_sem_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/parse_kernel_names.postprocess.o \
	$O/laplacian_smoothing_sem.postprocess.o \
	$(EMPTY_MACRO)

xsmooth_laplacian_sem_SHARED_OBJECTS = \
	$O/specfem3D_par.solverstatic_module.o \
	$O/read_mesh_parameters.solverstatic.o \
	$O/shared_par.shared_module.o \
	$O/parallel.sharedmpi.o \
	$O/assemble_MPI_scalar.shared.o \
	$O/auto_ner.shared.o \
	$O/count_elements.shared.o \
	$O/count_points.shared.o \
	$O/define_all_layers.shared.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/get_model_parameters.shared.o \
	$O/get_timestep_and_layers.shared.o \
	$O/gll_library.shared.o \
	$O/lagrange_poly.shared.o \
	$O/model_mars_1D.shared.o \
	$O/model_prem.shared.o \
	$O/model_vpremoon.shared.o \
	$O/param_reader.cc.o \
	$O/read_compute_parameters.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$(EMPTY_MACRO)

${E}/xsmooth_laplacian_sem: $(xsmooth_laplacian_sem_OBJECTS) $(xsmooth_laplacian_sem_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)


##
## xsmooth_laplacian_sem_adios
##
xsmooth_laplacian_sem_adios_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/parse_kernel_names.postprocess.o \
	$O/laplacian_smoothing_sem.postprocess_adios.o \
	$(EMPTY_MACRO)

xsmooth_laplacian_sem_adios_SHARED_OBJECTS = \
	$(xsmooth_laplacian_sem_SHARED_OBJECTS)

xsmooth_laplacian_sem_adios_SHARED_OBJECTS += \
	$O/adios_helpers_addons.shared_adios_cc.o \
	$O/adios_helpers_definitions.shared_adios.o \
	$O/adios_helpers_readers.shared_adios.o \
	$O/adios_helpers_writers.shared_adios.o \
	$O/adios_helpers.shared_adios.o \
	$O/adios_manager.shared_adios_module.o \
	$(EMPTY_MACRO)

${E}/xsmooth_laplacian_sem_adios: $(xsmooth_laplacian_sem_adios_OBJECTS) $(xsmooth_laplacian_sem_adios_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)


##
## xcreate_cross_section
##
xcreate_cross_section_OBJECTS = \
	$O/postprocess_par.postprocess_module.o \
	$O/create_cross_section.postprocess.o \
	$O/parse_kernel_names.postprocess.o \
	$(EMPTY_MACRO)

xcreate_cross_section_SHARED_OBJECTS = \
	$O/specfem3D_par.solverstatic_module.o \
	$O/read_mesh_parameters.solverstatic.o \
	$O/parallel.sharedmpi.o \
	$O/shared_par.shared_module.o \
	$O/binary_c_io.cc.o \
	$O/auto_ner.shared.o \
	$O/count_elements.shared.o \
	$O/count_points.shared.o \
	$O/create_name_database.shared.o \
	$O/define_all_layers.shared.o \
	$O/exit_mpi.shared.o \
	$O/flush_system.shared.o \
	$O/get_model_parameters.shared.o \
	$O/get_timestep_and_layers.shared.o \
	$O/gll_library.shared.o \
	$O/heap_sort.shared.o \
	$O/hex_nodes.shared.o \
	$O/interpolate.shared.o \
	$O/intgrl.shared.o \
	$O/lagrange_poly.shared.o \
	$O/make_ellipticity.shared.o \
	$O/model_mars_1D.shared.o \
	$O/model_prem.shared.o \
	$O/model_Sohl.shared.o \
	$O/model_topo_bathy.shared.o \
	$O/model_vpremoon.shared.o \
	$O/param_reader.cc.o \
	$O/read_compute_parameters.shared.o \
	$O/read_parameter_file.shared.o \
	$O/read_value_parameters.shared.o \
	$O/recompute_jacobian.shared.o \
	$O/reduce.shared.o \
	$O/rthetaphi_xyz.shared.o \
	$O/search_kdtree.shared.o \
	$O/spline_routines.shared.o \
	$(EMPTY_MACRO)

# extra dependencies
$O/create_cross_section.postprocess.o: $O/search_kdtree.shared.o

${E}/xcreate_cross_section: $(xcreate_cross_section_OBJECTS) $(xcreate_cross_section_SHARED_OBJECTS)
	${MPIFCCOMPILE_CHECK} -o $@ $+ $(MPILIBS)


#######################################

###
### module dependencies
###
$O/postprocess_par.postprocess_module.o: $O/shared_par.shared_module.o
$O/postprocess_par.postprocess_module.o: $O/specfem3D_par.solverstatic_module.o

####
#### rule for each .o file below
####

$O/%.postprocess_module.o: $S/%.f90 $O/shared_par.shared_module.o $O/specfem3D_par.solverstatic_module.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.postprocess.o: $S/%.f90 $O/postprocess_par.postprocess_module.o $O/parallel.sharedmpi.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

$O/%.postprocess.o: $S/%.F90 $O/postprocess_par.postprocess_module.o $O/parallel.sharedmpi.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

### ADIOS
$O/%.postprocess_adios.o: $S/%.F90 $O/postprocess_par.postprocess_module.o $O/parallel.sharedmpi.o $O/adios_helpers.shared_adios.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $< $(FC_DEFINE)USE_ADIOS_INSTEAD_OF_MESH

$O/%.postprocess_adios.o: $S/%.f90 $O/postprocess_par.postprocess_module.o $O/parallel.sharedmpi.o $O/adios_helpers.shared_adios.o
	${FCCOMPILE_CHECK} ${FCFLAGS_f90} -c -o $@ $<

