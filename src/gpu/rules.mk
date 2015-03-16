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
	$O/compute_stacey_acoustic_gpu.o \
	$O/compute_stacey_elastic_gpu.o \
	$O/compute_strain_gpu.o \
	$O/helper_functions_gpu.o \
	$O/initialize_gpu.o \
	$O/noise_tomography_gpu.o \
	$O/prepare_mesh_constants_gpu.o \
	$O/transfer_fields_gpu.o \
	$O/update_displacement_gpu.o \
	$O/write_seismograms_gpu.o \
	$O/save_and_compare_cpu_vs_gpu.o \
	$(EMPTY_MACRO)

ifeq ($(CUDA),yes)
  cuda_specfem3D_DEVICE_OBJ =  $O/cuda_device_obj.o 
  include $(BOAST_DIR)/kernel_cuda.mk # defines $(cuda_kernels_OBJS)
endif

#######################################

ifeq ($(GPU_CUDA_AND_OCL),yes)
  gpu_specfem3D_OBJECTS:=$(subst .o,.cuda-ocl.o,${gpu_specfem3D_OBJECTS})
endif

ifneq ($(GPU_CUDA_AND_OCL),yes)
  ifeq ($(OCL), yes)
    gpu_specfem3D_OBJECTS:=$(subst .o,.ocl.o,${gpu_specfem3D_OBJECTS})
  endif
  ifeq ($(CUDA),yes)
    gpu_specfem3D_OBJECTS:=$(subst .o,.cuda.o,${gpu_specfem3D_OBJECTS})
  endif
endif

gpu_specfem3D_OBJECTS += $(cuda_specfem3D_DEVICE_OBJ) $(cuda_kernels_OBJS)

###
### variables
###

NVCC_CFLAGS := ${NVCC_FLAGS} -x cu

BUILD_VERSION_TXT := with
SELECTOR_CFLAG :=

ifeq ($(CUDA),yes)
  BUILD_VERSION_TXT += Cuda
  SELECTOR_CFLAG += -DUSE_CUDA

  ifeq ($(CUDA5),yes)
    BUILD_VERSION_TXT += (v5)
  endif
endif

ifeq ($(GPU_CUDA_AND_OCL),yes)
  BUILD_VERSION_TXT += and
endif

ifeq ($(OCL), yes)
  BUILD_VERSION_TXT += OpenCL
  LDFLAGS += $(OCL_LINK)
  OCL_CPU_FLAGS += $(OCL_INC)
  SELECTOR_CFLAG += -DUSE_OPENCL
  ifneq ($(strip $(OCL_GPU_FLAGS)),)
    SELECTOR_CFLAG += -DOCL_GPU_CFLAGS="$(OCL_GPU_FLAGS)"
  endif
  ifeq ($(CUDA),yes)
    CUDA_LINK += $(OCL_LINK)
    NVCC_CFLAGS += $(OCL_CPU_FLAGS)
  endif
endif

BUILD_VERSION_TXT += support

###
### building rules
###


###
### boast kernel generation
###
GPU_ELEM_PER_THREAD=2
boast_kernels :
	@echo ""
	@echo "building boast kernels: in directory $(BOAST_DIR_NAME)"
	@echo ""
	cd src/gpu/boast ;\
	mkdir -p ../$(BOAST_DIR_NAME);\
	ruby kernels.rb --output-dir ../$(BOAST_DIR_NAME) -e $(GPU_ELEM_PER_THREAD)
	@echo ""

test_boast_kernels :
	@echo ""
	@echo "building and testing boast kernels: in directory $(BOAST_DIR_NAME)"
	@echo ""
	cd src/gpu/boast ;\
	mkdir -p ../$(BOAST_DIR_NAME);\
	ruby kernels.rb --output-dir ../$(BOAST_DIR_NAME) --check
	@echo ""

###
### compilation
###

ifeq ($(CUDA),yes)
$O/%.cuda-kernel.o: $(BOAST_DIR)/%.cu $S/mesh_constants_gpu.h $S/mesh_constants_cuda.h
	$(NVCC) -c $< -o $@ $(NVCC_CFLAGS) -I${SETUP} -I$(BOAST_DIR) $(SELECTOR_CFLAG) -include $(word 2,$^)

$(cuda_specfem3D_DEVICE_OBJ): $(subst $(cuda_specfem3D_DEVICE_OBJ), ,$(gpu_specfem3D_OBJECTS)) $(cuda_kernels_OBJS)
	${NVCCLINK} -o $@ $^
endif

$O/%.cuda-ocl.o: $O/%.cuda.o
	cd $O && cp $(shell basename $<) $(shell basename $@)

$O/%.ocl.o: $S/%.c ${SETUP}/config.h $S/mesh_constants_gpu.h $S/mesh_constants_ocl.h
	${CC} -c $< -o $@ $(OCL_CPU_FLAGS) -I${SETUP} -I$(BOAST_DIR) $(SELECTOR_CFLAG)

$O/%.cuda.o: $S/%.c ${SETUP}/config.h $S/mesh_constants_gpu.h
	$(NVCC) -c $< -o $@ $(NVCC_CFLAGS) -I${SETUP} -I$(BOAST_DIR) $(SELECTOR_CFLAG)

print-%:
	@echo '$*=$($*)'
