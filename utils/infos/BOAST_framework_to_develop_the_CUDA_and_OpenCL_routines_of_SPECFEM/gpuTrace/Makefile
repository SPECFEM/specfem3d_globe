################
### Config  ####
################

WITH_OCL = yes
WITH_CUDA = yes

################
### Defaults  ##
################

SO_NAME := ldChecker.so

MPI_INC := -I/usr/include/openmpi-x86_64/ -I/usr/lib/openmpi/include

PY_CFLAGS := $(shell python3-config --includes) -DPYTHON_MOD_PATH=$(shell pwd)
PY_LDFLAGS := $(shell python3-config --libs)

SO_CFLAGS := -fPIC
SO_LDFLAGS := -fPIC -rdynamic -shared  $(PY_LDFLAGS)

CFLAGS := $(SO_CFLAGS) -g -O0 -std=gnu99 -Werror -Wall -pedantic -Wno-format-security $(MPI_INC)

CUDA_CFLAGS := -I/usr/local/cuda-5.5/targets/x86_64-linux/include/ -I/usr/local/cuda-5.5/include

#CUDA_HELPER = cuda_helper_preprocessed.c
CUDA_HELPER = gpu_helper_py.o

################
### Common  ####
################

ifeq ($(WITH_OCL), yes)
GPU_INSTR := $(GPU_INSTR) instr-ocl.o
HELPER_HEADER := $(HELPER_HEADER) ocl_helper.h
endif

ifeq ($(WITH_CUDA), yes)
GPU_INSTR := $(GPU_INSTR) instr-cuda.o
HELPER_HEADER := $(HELPER_HEADER) cuda_helper.h
endif

all : $(SO_NAME)

$(SO_NAME) : $(GPU_INSTR) gpu_helper_py.o ldChecker.o
	gcc -o $@ $^ $(SO_LDFLAGS) -lpython3.3m

gpu_helper_py.o : gpu_helper_py.c $(HELPER_HEADER)
	gcc -o $@ -c $< $(CFLAGS) $(PY_CFLAGS)

ldChecker.o : ldChecker.c ldChecker.h $(HELPER_HEADER)
	gcc -o $@ -c $< $(CFLAGS) $(SO_CFLAGS) $(MPI_INC)

clean : clean-python
	rm -fv *.o *.so

clean-python :
	rm -rf __pycache__

get_env_preload :
	@echo LD_PRELOAD=$(shell pwd)/$(SO_NAME):libpython3.3m.so

################
### OpenCL  ####
################

instr-ocl.o : instr-ocl.c ocl_helper.h ldChecker.h
	gcc -o $@ -c $< $(CFLAGS)

##############
### CUDA  ####
##############

cuda_helper_preprocessed.c : cuda_helper_preprocessor.py $(CUDA_BIN)
	python3 $< $(CUDA_BIN) > $@

instr-cuda.o : instr-cuda.c $(CUDA_HELPER) ldChecker.o
	gcc -o $@ -c $< $(CFLAGS) $(CUDA_CFLAGS)


