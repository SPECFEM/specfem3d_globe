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

## compilation directories

S := ${S_TOP}/src/gpu

BOAST_DIR_NAME := kernels.gen
BOAST_DIR := ${S}/${BOAST_DIR_NAME}

#######################################

gpu_specfem3D_OBJECTS = \
	$O/assemble_MPI_scalar_gpu.o \
	$O/assemble_MPI_vector_gpu.o \
	$O/check_fields_gpu.o \
	$O/compute_add_sources_elastic_gpu.o \
	$O/compute_coupling_gpu.o \
	$O/compute_forces_crust_mantle_gpu.o \
	$O/compute_forces_inner_core_gpu.o \
	$O/compute_forces_outer_core_gpu.o \
	$O/compute_kernels_gpu.o \
	$O/compute_seismograms_gpu.o \
	$O/compute_stacey_acoustic_gpu.o \
	$O/compute_stacey_elastic_gpu.o \
	$O/compute_strain_gpu.o \
	$O/helper_functions_gpu.o \
	$O/initialize_gpu.o \
	$O/noise_tomography_gpu.o \
	$O/prepare_mesh_constants_gpu.o \
	$O/smooth_gpu.o \
	$O/transfer_fields_gpu.o \
	$O/update_displacement_gpu.o \
	$O/update_displacement_LDDRK_gpu.o \
	$O/write_seismograms_gpu.o \
	$O/save_and_compare_cpu_vs_gpu.o \
	$(EMPTY_MACRO)

gpu_specfem3D_STUBS = \
	$O/specfem3D_gpu_method_stubs.gpu_cc.o \
	$(EMPTY_MACRO)

# CUDA kernel files
ifeq ($(CUDA),yes)
  cuda_specfem3D_DEVICE_OBJ =  $O/cuda_device_obj.o

  # defines $(cuda_kernels_OBJS)
  include $(BOAST_DIR)/kernel_cuda.mk
endif

ifeq ($(HIP),yes)
  # defines $(cuda_kernels_OBJS)
  include $(BOAST_DIR)/kernel_cuda.mk
  # renames endings
  cuda_kernels_OBJS:=$(subst .cuda-kernel.o,.hip-kernel.o,${cuda_kernels_OBJS})
endif


ifdef NO_GPU
gpu_OBJECTS = $(gpu_specfem3D_STUBS)
else
gpu_OBJECTS = $(gpu_specfem3D_OBJECTS)
endif


#######################################

# substitutes object endings to assign corresponding compilation rule
ifeq ($(GPU_CUDA_AND_OCL),yes)
	# combines both CUDA and OpenCL kernels compilation
  gpu_specfem3D_OBJECTS:=$(subst .o,.cuda-ocl.o,${gpu_specfem3D_OBJECTS})
endif

ifneq ($(GPU_CUDA_AND_OCL),yes)
  # OpenCL kernels only
  ifeq ($(OCL), yes)
    gpu_specfem3D_OBJECTS:=$(subst .o,.ocl.o,${gpu_specfem3D_OBJECTS})
  endif

  # CUDA kernels only
  ifeq ($(CUDA),yes)
    gpu_specfem3D_OBJECTS:=$(subst .o,.cuda.o,${gpu_specfem3D_OBJECTS})
  endif

  # HIP kernels only
  ifeq ($(HIP), yes)
    gpu_specfem3D_OBJECTS:=$(subst .o,.hip.o,${gpu_specfem3D_OBJECTS})
  endif
endif

gpu_specfem3D_OBJECTS += $(cuda_specfem3D_DEVICE_OBJ) $(cuda_kernels_OBJS)

###
### variables
###

ifeq ($(HAS_GPU), yes)
	BUILD_VERSION_TXT := with
endif
SELECTOR_CFLAG :=

## CUDA compilation
NVCC_CFLAGS := ${NVCC_FLAGS} -x cu
ifeq ($(CUDA),yes)
  BUILD_VERSION_TXT += Cuda
  SELECTOR_CFLAG += $(FC_DEFINE)USE_CUDA
  GPU_LINK = $(CUDA_LINK)

  ifeq ($(CUDA5),yes)
    BUILD_VERSION_TXT += (v5)
  endif
  ifeq ($(CUDA6),yes)
    BUILD_VERSION_TXT += (v6)
  endif
  ifeq ($(CUDA7),yes)
    BUILD_VERSION_TXT += (v7)
  endif
  ifeq ($(CUDA8),yes)
    BUILD_VERSION_TXT += (v8)
  endif
  ifeq ($(CUDA9),yes)
	  BUILD_VERSION_TXT += (v9)
  endif
  ifeq ($(CUDA10),yes)
	  BUILD_VERSION_TXT += (v10)
  endif
  ifeq ($(CUDA11),yes)
	  BUILD_VERSION_TXT += (v11)
  endif
  ifeq ($(CUDA12),yes)
	  BUILD_VERSION_TXT += (v12)
  endif
endif

ifeq ($(GPU_CUDA_AND_OCL),yes)
  BUILD_VERSION_TXT += and
endif

## OpenCL compilation
ifeq ($(OCL), yes)
  BUILD_VERSION_TXT += OpenCL
  OCL_CPU_FLAGS += $(OCL_INC)
  SELECTOR_CFLAG += $(FC_DEFINE)USE_OPENCL
  GPU_LINK = $(OCL_LINK)

  ifneq ($(strip $(OCL_GPU_FLAGS)),)
    SELECTOR_CFLAG += -DOCL_GPU_CFLAGS="$(OCL_GPU_FLAGS)"
  endif
  ifeq ($(CUDA),yes)
    GPU_LINK += $(CUDA_LINK)
    NVCC_CFLAGS += $(OCL_CPU_FLAGS)
  endif
endif

## HIP compilation
HIPCC_CFLAGS := ${HIP_CFLAGS}      # ${HIP_CFLAG_ENDING} adds -x hip depending on platform
ifeq ($(HIP), yes)
  BUILD_VERSION_TXT += HIP
  SELECTOR_CFLAG += $(FC_DEFINE)USE_HIP
  ifneq ($(strip $(HIP_GPU_FLAGS)),)
    SELECTOR_CFLAG += -DHIP_GPU_CFLAGS="$(HIP_GPU_FLAGS)"
  endif
  GPU_LINK = $(HIP_LINK)

  # todo: compile hip with nvcc
  #ifeq ($(CUDA),yes)
  #  GPU_LINK += $(CUDA_LINK)
  #  NVCC_CFLAGS += $(HIP_CPU_FLAGS)
  #endif
endif

ifeq ($(HAS_GPU), yes)
	BUILD_VERSION_TXT += support
endif

###
### building rules
###

.PHONY: boast

boast: boast_kernels

###
### boast kernel generation
###

boast_kernels :
	@echo ""
	@echo "building boast kernels: in directory $(BOAST_DIR_NAME)"
	@echo ""
	cd src/gpu/boast ;\
	mkdir -p ../$(BOAST_DIR_NAME);\
	ruby kernels.rb --output-dir ../$(BOAST_DIR_NAME) --elem $(GPU_ELEM_PER_THREAD)
	@echo ""

test_boast_kernels :
	@echo ""
	@echo "building and testing boast kernels: in directory $(BOAST_DIR_NAME)"
	@echo ""
	cd src/gpu/boast ;\
	mkdir -p ../$(BOAST_DIR_NAME);\
	ruby kernels.rb --output-dir ../$(BOAST_DIR_NAME) --check
	@echo ""

#######################################

####
#### rule for each .o file below
####

###
### CUDA compilation
###

# source kernel files
# cuda kernels
ifeq ($(CUDA),yes)
$O/%.cuda-kernel.o: $(BOAST_DIR)/%.cu $S/mesh_constants_gpu.h $S/mesh_constants_cuda.h $(BOAST_DIR)/kernel_proto.cu.h
	$(NVCC) -c $< -o $@ $(NVCC_CFLAGS) -I${SETUP} -I$(BOAST_DIR) $(SELECTOR_CFLAG) -include $(word 2,$^)

$(cuda_specfem3D_DEVICE_OBJ): $(subst $(cuda_specfem3D_DEVICE_OBJ), ,$(gpu_specfem3D_OBJECTS)) $(cuda_kernels_OBJS)
	${NVCCLINK} -o $@ $^
endif

ifeq ($(HIP),yes)
$O/%.hip-kernel.o: $(BOAST_DIR)/%.cpp $S/mesh_constants_gpu.h $S/mesh_constants_hip.h $(BOAST_DIR)/kernel_proto.cu.h
	$(HIPCC) -c $< -o $@ $(HIP_CFLAGS) -I${SETUP} -I$(BOAST_DIR) $(SELECTOR_CFLAG) -include $(word 2,$^)
endif

$O/%.cuda-ocl.o: $O/%.cuda.o
	cd $O && cp $(shell basename $<) $(shell basename $@)

# source files in src/gpu/
$O/%.ocl.o: $S/%.c ${SETUP}/config.h $S/mesh_constants_gpu.h $S/mesh_constants_ocl.h
	${CC} -c $< -o $@ $(OCL_CPU_FLAGS) -I${SETUP} -I$(BOAST_DIR) $(SELECTOR_CFLAG)

$O/%.cuda.o: $S/%.c ${SETUP}/config.h $S/mesh_constants_gpu.h
	$(NVCC) -c $< -o $@ $(NVCC_CFLAGS) -I${SETUP} -I$(BOAST_DIR) $(SELECTOR_CFLAG)

$O/%.hip.o: $S/%.c ${SETUP}/config.h $S/mesh_constants_gpu.h  $S/mesh_constants_hip.h
	${HIPCC} ${HIP_CFLAG_ENDING} -c $< -o $@ $(HIPCC_CFLAGS) -I${SETUP} -I$(BOAST_DIR) $(SELECTOR_CFLAG)

# C version
$O/%.gpu_cc.o: $S/%.c ${SETUP}/config.h
	${CC} -c $(CPPFLAGS) $(CFLAGS) $(MPI_INCLUDES) -o $@ $<

print-%:
	@echo '$*=$($*)'
