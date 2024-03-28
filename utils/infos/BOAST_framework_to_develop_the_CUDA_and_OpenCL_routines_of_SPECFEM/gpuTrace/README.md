gpuTrace
========

In the same spirit as strace and ltrace, gpuTrace will trace the interactions between Cuda and OpenCL applications and the GPU:

candidate_generation<1,32><1,1>(
     unsigned int iteration=<4>
     unsigned int ncols=<5>
     unsigned int nrows=<55>
     unsigned int*frequent_bitmap=<buffer #0 268435520, 128, 0, 0, 0, 268435520, 524
     unsigned int*candidate_bitmap=<buffer #1 268435520, 128, 0, 0, 0, 268435520, 52
     unsigned int*candidate_count=<buffer #3 4, 3, 2, 1, 0, 8, 7, 6, ...>
     ----
     <out> unsigned int*frequent_bitmap=<buffer #0 268435520, 128, 0, 0, 0, 268435520, 524
     <out> unsigned int*candidate_bitmap=<buffer #1 268435520, 524416, 0, 0, 0, 268435520,
     <out> unsigned int*candidate_count=<buffer #3 3, 2, 1, 0, 2, 1, 0, 1, ...>
);

support_count<1,1><1,53>(
     unsigned int ncols=<100>
     unsigned int nrows=<53>
     unsigned int*bitmapOp1=<buffer #4 4294967295, 4294967295, 4294967295, 429496729
     unsigned int*bitmapOp2=<buffer #5 4294967295, 4294967295, 4294967295, 429496729
     unsigned int*bitmapRst=<buffer #6 4294967295, 4294967295, 4294967295, 429496729
     unsigned int*lookup_table=<buffer #7 0, 1, 1, 2, 1, 2, 2, 3, ...>
     unsigned int*support_table=<buffer #8 3043, 3058, 3068, 3024, 3039, 3049, 3006,
     ----
     <out> unsigned int*bitmapOp1=<buffer #4 4294967295, 4294967295, 4294967295, 429496729
     <out> unsigned int*bitmapOp2=<buffer #5 4294967295, 4294967295, 4294967295, 429496729
     <out> unsigned int*bitmapRst=<buffer #6 4294967295, 4294967295, 4294967295, 429496729
     <out> unsigned int*lookup_table=<buffer #7 0, 1, 1, 2, 1, 2, 2, 3, ...>
     <out> unsigned int*support_table=<buffer #8 3032, 3042, 2999, 3057, 3013, 3023, 3038,
);


Building the environment
========================

* if OpenCL or Cuda is not available in your environment, set `Makefile::WITH_OCL` or  `Makefile::WITH_CUDA` to `no`.

* `make all` builds `gpuTrace`'s library

Tracing the execution
=====================

* `gpuTrace` works with a preloaded library, so you need to set the environment variable `LD_PRELOAD` to `ldChecker.so`.

* `make get_env_preload` prints the definition of the variable (LD_PRELOAD=..).
  * with GDB, use the command `set environment LD_PRELOAD=...` before running the application (and use `start` to let GDB load the symbols of the library.)

* run the application with `LD_PRELOAD=... bin/my_app`

Configuration
=============

The configuration of gpuTrace is currently exclusively done at compile time (that should change one day).


GPU Tracing
-----------

* print create and transfer operations
  * PRINT_KERNEL_BUFFER_CREATION
  * PRINT_BUFFER_TRANSFER
  * PRINT_BUFFER_RELEASE
    * PRINT_BUFFER_TRANSFER_FIRST_BYTES_AS_FLOAT (buffers are typeless)

* print kernel execution
  * PRINT_KERNEL_BEFORE_EXEC
  * PRINT_KERNEL_AFTER_EXEC
    * PRINT_KERNEL_AFTER_EXEC_IGNORE_CONST (const are not supposed to change, so keep it to 0)
  * PRINT_KERNEL_NAME_ONLY (lightweight tracing)

* filters for kernel tracing
  * USE_ADVANCED_KERNEL_FILTER
    * ENV__KERNEL_FILTER = "LD_KERNEL_FILTER"
    * environment variable used to select which kernels should be traced.
    * format is "<name>=<exec_counter>:*"
    * which means that only kernels whose name are in the "<name>" list will be traced, and only when they execution counter is equal to "<exec_counter>".
  * otherwise:
    * FILTER_BY_KERNEL_EXEC_CPT
      * KERNEL_EXEC_CPT_LOWER_BOUND
      * KERNEL_EXEC_CPT_UPPER_BOUND
      * only kernel exection with execution counts between lower and upper will be traced
    * FILTER_BY_KERNEL_NAME
      * KERNEL_NAME_FILTER
      * only kernels with this string in their names will be traced
    * filters on name and execution counter filters can be used simultaneously

* BUFFER_ZERO_IS_NULL (OpenCL hack) first buffer created should be handled as NULL pointer
* FORCE_FINISH_KERNEL (OpenCL hack) should be activated if more than one queue is used

* WITH_MPI disable if you're environment doesn't offer MPI headers
  * ONLY_MPI_ROOT_OUTPUT only MPI process with rank == 0 will print messages

* ENABLE_KERNEL_PROFILING reports profiling information at the end of the execution
  * ENABLE_LEAK_DETECTION does refcounting on OpenCL objects

GPU Testing
-----------

Step 1: gather kernel parametere values

* PRINT_KERNEL_BEFORE_EXEC 1
* PRINT_KERNEL_AFTER_EXEC 1
* PRINT_KERNEL_ARG_FULL_BUFFER 1
  * FULL_BUFFER_SIZE_LIMIT 0
* PRINT_KERNEL_PARAMS_TO_FILE 1
* ENABLE_KERNEL_TESTING 0
* ENV__PARAM_FILE_SUBDIR_PREFIX "LD_KERNEL_SUBDIR_PREFIX"

then

* set filters as you need
* give a name to the dataset with environment variable "LD_KERNEL_SUBDIR_PREFIX"
* run the application

finally

* collect the data in "kernel_test_data"

Step 2: test kernel execution

* [disable tracing options]
* disable buffer collection
  * PRINT_KERNEL_PARAMS_TO_FILE 0
  * PRINT_KERNEL_ARG_FULL_BUFFER 0
* enable kernel testing
  * ENABLE_KERNEL_TESTING 1
* run the application
  * all the dataset for all kernels will be tested
  * currently there is no "testing only" option, so I put an `error` quit in the code, for instance in `kernel_executed_event`.