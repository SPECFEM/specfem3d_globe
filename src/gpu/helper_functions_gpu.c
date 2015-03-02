/*
!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================
*/

#include "mesh_constants_gpu.h"


/* ----------------------------------------------------------------------------------------------- */
// OpenCL setup for memset function
/* ----------------------------------------------------------------------------------------------- */

#ifdef USE_OPENCL

// kernel for memset in OpenCL
const char *memset_kern_code[] = { "\
__kernel void memset_uint4(__global int *mem, const int size, __private int val) { \n\
int tid = get_local_id(0); \n\
int bx = (get_group_id(1)) * (get_num_groups(0)) + get_group_id(0); \n\
int i = tid + (bx) * (get_local_size(0)); \n\
//debug \n\
//if (i == 0) { printf(\"memset size = %i value = %i buffer %i \\n\",size,val,mem[0]); } \n\
if (i < size ) { mem[i]=val; } \n\
}" };


/* ----------------------------------------------------------------------------------------------- */

cl_kernel *setup_ocl_memset (int do_setup) {

  static int inited = 0;
  static cl_kernel memset_kern;
  cl_int errcode;

  if (do_setup) {
    if (!inited) {
      // creates openCL kernel
      cl_program memset_program = clCreateProgramWithSource(mocl.context, 1,
                                                            memset_kern_code, 0,
                                                            clck_(&errcode));
      clCheck (clBuildProgram (memset_program, 0, NULL, NULL, NULL, NULL));
      memset_kern = clCreateKernel (memset_program, "memset_uint4", clck_(&errcode));
      inited = 1;
    }
  } else {
    // releases kernel
    if (inited) { clCheck(clReleaseKernel (memset_kern)); }
  }

  return &memset_kern;
}

/* ----------------------------------------------------------------------------------------------- */

void moclEnqueueFillBuffer (cl_mem *buffer, int val, size_t size_byte) {

  // creates/gets OpenCL memset kernel
  cl_kernel *memset_kern = setup_ocl_memset(1);

  // value to fill buffer
  cl_int value = val;

  // gets size as number of integer values
  int size;
  size = size_byte / sizeof(cl_int);

  size_t global_work_size[2];
  size_t local_work_size[2];
  cl_uint idx = 0;

  clCheck (clSetKernelArg (*memset_kern, idx++, sizeof (cl_mem), (void *) buffer));
  clCheck (clSetKernelArg (*memset_kern, idx++, sizeof (cl_int), (void *) &size));
  clCheck (clSetKernelArg (*memset_kern, idx++, sizeof (cl_int), (void *) &value));

  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int) ceil ((double) size / (double) blocksize)) * blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

  local_work_size[0] = blocksize;
  local_work_size[1] = 1;
  global_work_size[0] = num_blocks_x * blocksize;
  global_work_size[1] = num_blocks_y;

  //debug
  //printf("moclEnqueueFillBuffer: size %i value %i - work_size %zu %zu \n",size,value,local_work_size[0],global_work_size[0]);

  clCheck (clEnqueueNDRangeKernel (mocl.command_queue, *memset_kern, 2, NULL,
                                   global_work_size, local_work_size, 0, NULL, NULL));
  // synchronizes
  clFinish (mocl.command_queue);
}

#endif // USE_OPENCL


/*----------------------------------------------------------------------------------------------- */
// GPU helper functions
/*----------------------------------------------------------------------------------------------- */

// copies integer array from CPU host to GPU device
void gpuCreateCopy_todevice_int (gpu_int_mem *d_array_addr_ptr, int *h_array, int size) {

  TRACE ("gpuCreateCopy_todevice_int");

#ifdef USE_OPENCL
  if (run_opencl) {
    cl_int errcode;

    // allocates memory on GPU
    d_array_addr_ptr->ocl = clCreateBuffer (mocl.context, CL_MEM_READ_WRITE,
                                            size * sizeof (int), NULL, clck_(&errcode));

    // copies values onto GPU
    clCheck (clEnqueueWriteBuffer (mocl.command_queue, d_array_addr_ptr->ocl, CL_TRUE, 0,
                                   size*sizeof (int), h_array, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    // allocates memory on GPU
    //
    // note: cudaMalloc uses a double-pointer, such that it can return an error code in case it fails
    //          we thus pass the address to the pointer above (as void double-pointer) to have it
    //          pointing to the correct pointer of the array here
    print_CUDA_error_if_any(cudaMalloc((void**)&d_array_addr_ptr->cuda,size*sizeof(int)),12001);

    // copies values onto GPU
    //
    // note: cudaMemcpy uses the pointer to the array, we thus re-cast the value of
    //          the double-pointer above to have the correct pointer to the array
    print_CUDA_error_if_any(cudaMemcpy((int*) d_array_addr_ptr->cuda,h_array,size*sizeof(int),cudaMemcpyHostToDevice),12002);
  }
#endif
}

/*----------------------------------------------------------------------------------------------- */

// copies real array from CPU host to GPU device
void gpuCreateCopy_todevice_realw (gpu_realw_mem *d_array_addr_ptr, realw *h_array, int size) {

  TRACE ("gpuCreateCopy_todevice_realw");

  // allocates memory on GPU
#ifdef USE_OPENCL
  if (run_opencl) {
    cl_int errcode;

    d_array_addr_ptr->ocl = clCreateBuffer (mocl.context, CL_MEM_READ_WRITE, size * sizeof (realw),
                                            NULL, clck_(&errcode));

    // copies values onto GPU
    clCheck (clEnqueueWriteBuffer (mocl.command_queue, d_array_addr_ptr->ocl, CL_TRUE, 0,
                                   size * sizeof (realw), h_array, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    // allocates memory on GPU
    print_CUDA_error_if_any(cudaMalloc((void**)&d_array_addr_ptr->cuda,size*sizeof(realw)),22001);
    // copies values onto GPU
    print_CUDA_error_if_any(cudaMemcpy((realw*) d_array_addr_ptr->cuda,h_array,size*sizeof(realw),cudaMemcpyHostToDevice),22002);
  }
#endif
}

/*----------------------------------------------------------------------------------------------- */

// copies real array from CPU host to GPU device
void gpuCopy_todevice_realw (gpu_realw_mem *d_array_addr_ptr, realw *h_array, int size) {

  TRACE ("gpuCopy_todevice_realw");

  // copies memory on from CPU to GPU
  // uses blocking copies
#ifdef USE_OPENCL
  if (run_opencl) {
    // copies values onto GPU
    clCheck (clEnqueueWriteBuffer (mocl.command_queue, d_array_addr_ptr->ocl, CL_TRUE, 0, size * sizeof (realw), h_array, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    // copies values onto GPU
    print_CUDA_error_if_any(cudaMemcpy((realw*) d_array_addr_ptr->cuda,h_array,size*sizeof(realw),cudaMemcpyHostToDevice),22003);
  }
#endif
}

/*----------------------------------------------------------------------------------------------- */

// copies double array from CPU host to GPU device
void gpuCopy_todevice_double (gpu_double_mem *d_array_addr_ptr, double *h_array, int size) {

  TRACE ("gpuCopy_todevice_double");

  // copies memory on from CPU to GPU
  // uses blocking copies
#ifdef USE_OPENCL
  if (run_opencl) {
    // copies values onto GPU
    clCheck (clEnqueueWriteBuffer (mocl.command_queue, d_array_addr_ptr->ocl, CL_TRUE, 0, size * sizeof (double), h_array, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    // copies values onto GPU
    print_CUDA_error_if_any(cudaMemcpy((double*) d_array_addr_ptr->cuda,h_array,size*sizeof(double),cudaMemcpyHostToDevice),22003);
  }
#endif
}

/*----------------------------------------------------------------------------------------------- */

// copies integer array from CPU host to GPU device
void gpuCopy_todevice_int (gpu_int_mem *d_array_addr_ptr, int *h_array, int size) {

  TRACE ("gpuCopy_todevice_int");

  // copies memory on from CPU to GPU
  // uses blocking copies
#ifdef USE_OPENCL
  if (run_opencl) {
    // copies values onto GPU
    clCheck (clEnqueueWriteBuffer (mocl.command_queue, d_array_addr_ptr->ocl, CL_TRUE, 0, size * sizeof (int), h_array, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    // copies values onto GPU
    print_CUDA_error_if_any(cudaMemcpy((int*) d_array_addr_ptr->cuda,h_array,size*sizeof(int),cudaMemcpyHostToDevice),22003);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// copies array from GPU to CPU
void gpuCopy_from_device_realw (gpu_realw_mem *d_array_addr_ptr, realw *h_array, int size) {

  TRACE ("gpuCopy_from_device_realw");

  // copies memory from GPU back to CPU
#ifdef USE_OPENCL
  if (run_opencl) {
    // blocking copy
    clCheck (clEnqueueReadBuffer (mocl.command_queue, d_array_addr_ptr->ocl, CL_TRUE, 0, sizeof (realw) * size, h_array, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    // note: cudaMemcpy implicitly synchronizes all other cuda operations
    print_CUDA_error_if_any(cudaMemcpy(h_array,d_array_addr_ptr->cuda, sizeof(realw)*size, cudaMemcpyDeviceToHost),33001);
  }
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// creates real array on GPU
void gpuMalloc_realw (gpu_realw_mem *buffer, int size) {

  TRACE ("gpuMalloc_realw");

  // allocates array on GPU
#ifdef USE_OPENCL
  if (run_opencl) {
    cl_int errcode;
    buffer->ocl = clCreateBuffer (mocl.context, CL_MEM_READ_WRITE, size * sizeof(realw), NULL, clck_(&errcode));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMalloc((void**)&buffer->cuda, size * sizeof(realw)), 44001);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// creates double array on GPU
void gpuMalloc_double (gpu_double_mem *buffer, int size) {

  TRACE ("gpuMalloc_double");

  // allocates array on GPU
#ifdef USE_OPENCL
  if (run_opencl) {
    cl_int errcode;
    buffer->ocl = clCreateBuffer (mocl.context, CL_MEM_READ_WRITE, size * sizeof(double), NULL, clck_(&errcode));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMalloc((void**)&buffer->cuda, size * sizeof(double)), 44002);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// creates double array on GPU
void gpuMalloc_int (gpu_int_mem *buffer, int size) {

  TRACE ("gpuMalloc_int");

  // allocates array on GPU
#ifdef USE_OPENCL
  if (run_opencl) {
    cl_int errcode;
    buffer->ocl = clCreateBuffer (mocl.context, CL_MEM_READ_WRITE, size * sizeof(int), NULL, clck_(&errcode));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMalloc((void**)&buffer->cuda, size * sizeof(int)), 44003);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// creates double array on GPU
void gpuMemset_realw (gpu_realw_mem *buffer, int size, int value) {

  TRACE ("gpuMemset_realw");

  // initializes values for array on GPU
#ifdef USE_OPENCL
  if (run_opencl) {
    moclEnqueueFillBuffer(&buffer->ocl, value, size * sizeof (realw));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemset(buffer->cuda, value, size*sizeof(realw)),44004);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// setup functions
void gpuSetConst (gpu_realw_mem *buffer, size_t size, realw *array) {

  TRACE ("gpuSetConst");

  // allocates array on GPU
#ifdef USE_OPENCL
  if (run_opencl) {
    cl_int errcode;
    buffer->ocl = clCreateBuffer (mocl.context, CL_MEM_READ_ONLY, size * sizeof(realw), NULL, clck_(&errcode));
    clCheck (clEnqueueWriteBuffer (mocl.command_queue, buffer->ocl, CL_TRUE, 0, size * sizeof(realw), array, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMalloc(&buffer->cuda, size * sizeof(realw)), 1400);
    print_CUDA_error_if_any(cudaMemcpy(buffer->cuda, array, size * sizeof(realw), cudaMemcpyHostToDevice),1401);
  }
#endif
}

/*----------------------------------------------------------------------------------------------- */

void gpuFree (void *d_array_addr_ptr) {

  TRACE ("gpuFree");

  // conversion to a generic gpu_mem pointer
  gpu_mem *gpu_ptr = (gpu_mem *) d_array_addr_ptr;

  // frees memory on GPU
#ifdef USE_OPENCL
  if (run_opencl) { clReleaseMemObject(gpu_ptr->ocl);  }
#endif
#ifdef USE_CUDA
  if (run_cuda) { cudaFree(gpu_ptr->cuda); }
#endif
}

/*----------------------------------------------------------------------------------------------- */

/*
void gpuFreeHost (void *d_array_addr_ptr) {

  TRACE ("gpuFree");

  // frees pinned memory on GPU
#ifdef USE_OPENCL
  if (run_opencl) { RELEASE_PINNED_BUFFER_OCL((cl_mem *) d_array_addr_ptr->ocl);  }
#endif
#ifdef USE_CUDA
  if (run_cuda) { cudaFreeHost(d_array_addr_ptr->cuda); }
#endif
}
*/

/* ----------------------------------------------------------------------------------------------- */

void gpuInitialize_buffers(Mesh *mp) {

#ifdef USE_OPENCL
  // sets OpenCL pointers to NULL
  #define INIT_DUMMY_BUFFER(_field_) mp->_field_.ocl = NULL;

  #define GPU_REALW_BUFFER INIT_DUMMY_BUFFER
  #define GPU_INT_BUFFER INIT_DUMMY_BUFFER
  #define GPU_DOUBLE_BUFFER INIT_DUMMY_BUFFER
  #include "gpu_buffer_list.c"
  #undef INIT_DUMMY_BUFFER
#endif
#ifdef USE_CUDA
  // sets CUDA pointers to NULL
  #define INIT_DUMMY_BUFFER(_field_) mp->_field_.cuda = NULL

  #define GPU_REALW_BUFFER INIT_DUMMY_BUFFER
  #define GPU_INT_BUFFER INIT_DUMMY_BUFFER
  #define GPU_DOUBLE_BUFFER INIT_DUMMY_BUFFER
  #include "gpu_buffer_list.c"
  #undef INIT_DUMMY_BUFFER
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// GPU reset
/* ----------------------------------------------------------------------------------------------- */

void gpuReset() {
  // releases previous contexts

  // opencl version
#ifdef USE_OPENCL
  if (run_opencl) clReleaseContext (mocl.context);
#endif

  // cuda version
#ifdef USE_CUDA
  if (run_cuda) {
#if CUDA_VERSION < 4000
    cudaThreadExit();
#else
    cudaDeviceReset();
#endif
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// GPU synchronization
/* ----------------------------------------------------------------------------------------------- */

void gpuSynchronize() {
  // synchronizes device

  // opencl version
#ifdef USE_OPENCL
  if (run_opencl) {
    clFinish (mocl.command_queue);
    if (GPU_ASYNC_COPY) clFinish (mocl.copy_queue);
  }
#endif

  // cuda version
#ifdef USE_CUDA
  if (run_cuda) {
#if CUDA_VERSION < 4000
    cudaThreadSynchronize();
#else
    cudaDeviceSynchronize();
#endif
  }
#endif

}

/*----------------------------------------------------------------------------------------------- */
// OpenCL helper
/*----------------------------------------------------------------------------------------------- */

#ifdef USE_OPENCL

cl_int mocl_errcode;

cl_int clGetLastError () {
  return mocl_errcode;
}

/* The OpenCL Extension Wrangler Library
 * https://code.google.com/p/clew/
 * MIT License
 * */

const char* clewErrorString (cl_int error) {
  static const char* strings[] = {
    // Error Codes
    "CL_SUCCESS"                                  //   0
    , "CL_DEVICE_NOT_FOUND"                         //  -1
    , "CL_DEVICE_NOT_AVAILABLE"                     //  -2
    , "CL_COMPILER_NOT_AVAILABLE"                   //  -3
    , "CL_MEM_OBJECT_ALLOCATION_FAILURE"            //  -4
    , "CL_OUT_OF_RESOURCES"                         //  -5
    , "CL_OUT_OF_HOST_MEMORY"                       //  -6
    , "CL_PROFILING_INFO_NOT_AVAILABLE"             //  -7
    , "CL_MEM_COPY_OVERLAP"                         //  -8
    , "CL_IMAGE_FORMAT_MISMATCH"                    //  -9
    , "CL_IMAGE_FORMAT_NOT_SUPPORTED"               //  -10
    , "CL_BUILD_PROGRAM_FAILURE"                    //  -11
    , "CL_MAP_FAILURE"                              //  -12

    , ""    //  -13
    , ""    //  -14
    , ""    //  -15
    , ""    //  -16
    , ""    //  -17
    , ""    //  -18
    , ""    //  -19

    , ""    //  -20
    , ""    //  -21
    , ""    //  -22
    , ""    //  -23
    , ""    //  -24
    , ""    //  -25
    , ""    //  -26
    , ""    //  -27
    , ""    //  -28
    , ""    //  -29

    , "CL_INVALID_VALUE"                            //  -30
    , "CL_INVALID_DEVICE_TYPE"                      //  -31
    , "CL_INVALID_PLATFORM"                         //  -32
    , "CL_INVALID_DEVICE"                           //  -33
    , "CL_INVALID_CONTEXT"                          //  -34
    , "CL_INVALID_QUEUE_PROPERTIES"                 //  -35
    , "CL_INVALID_COMMAND_QUEUE"                    //  -36
    , "CL_INVALID_HOST_PTR"                         //  -37
    , "CL_INVALID_MEM_OBJECT"                       //  -38
    , "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR"          //  -39
    , "CL_INVALID_IMAGE_SIZE"                       //  -40
    , "CL_INVALID_SAMPLER"                          //  -41
    , "CL_INVALID_BINARY"                           //  -42
    , "CL_INVALID_BUILD_OPTIONS"                    //  -43
    , "CL_INVALID_PROGRAM"                          //  -44
    , "CL_INVALID_PROGRAM_EXECUTABLE"               //  -45
    , "CL_INVALID_KERNEL_NAME"                      //  -46
    , "CL_INVALID_KERNEL_DEFINITION"                //  -47
    , "CL_INVALID_KERNEL"                           //  -48
    , "CL_INVALID_ARG_INDEX"                        //  -49
    , "CL_INVALID_ARG_VALUE"                        //  -50
    , "CL_INVALID_ARG_SIZE"                         //  -51
    , "CL_INVALID_KERNEL_ARGS"                      //  -52
    , "CL_INVALID_WORK_DIMENSION"                   //  -53
    , "CL_INVALID_WORK_GROUP_SIZE"                  //  -54
    , "CL_INVALID_WORK_ITEM_SIZE"                   //  -55
    , "CL_INVALID_GLOBAL_OFFSET"                    //  -56
    , "CL_INVALID_EVENT_WAIT_LIST"                  //  -57
    , "CL_INVALID_EVENT"                            //  -58
    , "CL_INVALID_OPERATION"                        //  -59
    , "CL_INVALID_GL_OBJECT"                        //  -60
    , "CL_INVALID_BUFFER_SIZE"                      //  -61
    , "CL_INVALID_MIP_LEVEL"                        //  -62
    , "CL_INVALID_GLOBAL_WORK_SIZE"                 //  -63
    , "CL_UNKNOWN_ERROR_CODE"
  };

  if (error >= -63 && error <= 0) {
    return strings[-error];
  } else {
    return strings[64];
  }
}

const char *clGetErrorString (cl_int error) {
  return clewErrorString (error);
}

#endif // USE_OPENCL


/* ----------------------------------------------------------------------------------------------- */
// CUDA helper
/* ----------------------------------------------------------------------------------------------- */

#ifdef USE_CUDA

void print_CUDA_error_if_any(cudaError_t err, int num) {
  if (cudaSuccess != err)
  {
    printf("\nCUDA error !!!!! <%s> !!!!! \nat CUDA call error code: # %d\n",cudaGetErrorString(err),num);
    fflush(stdout);

    // outputs error file
    FILE* fp;
    int myrank;
    char filename[BUFSIZ];
#ifdef WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#else
    myrank = 0;
#endif
    sprintf(filename,"OUTPUT_FILES/error_message_%06d.txt",myrank);
    fp = fopen(filename,"a+");
    if (fp != NULL) {
      fprintf(fp,"\nCUDA error !!!!! <%s> !!!!! \nat CUDA call error code: # %d\n",cudaGetErrorString(err),num);
      fclose(fp);
    }

    // stops program
#ifdef WITH_MPI
    MPI_Abort(MPI_COMM_WORLD,1);
#endif
    exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------------------------- */

// Timing helper functions

void start_timing_cuda(cudaEvent_t* start,cudaEvent_t* stop) {
  // creates & starts event
  cudaEventCreate(start);
  cudaEventCreate(stop);
  cudaEventRecord( *start, 0 );
}

/* ----------------------------------------------------------------------------------------------- */

void stop_timing_cuda(cudaEvent_t* start,cudaEvent_t* stop, char* info_str) {
  realw time;
  // stops events
  cudaEventRecord( *stop, 0 );
  cudaEventSynchronize( *stop );
  cudaEventElapsedTime( &time, *start, *stop );
  cudaEventDestroy( *start );
  cudaEventDestroy( *stop );
  // user output
  printf("%s: Execution Time = %f ms\n",info_str,time);
}

#endif // USE_CUDA


/* ----------------------------------------------------------------------------------------------- */
// exit functions
/* ----------------------------------------------------------------------------------------------- */

void exit_on_gpu_error (char *kernel_name) {
  //check to catch errors from previous operations
  // /!\ in opencl, we can't have information about the last ASYNC error
  int error = 0;
  const char *strerr = NULL;

#ifdef USE_OPENCL
  if (run_opencl) {
    clFinish (mocl.command_queue);
    cl_int err = clGetLastError ();

    error = err != CL_SUCCESS;
    strerr = clGetErrorString (err);
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    cudaThreadSynchronize();
    cudaError_t err = cudaGetLastError();
    error = err != cudaSuccess;
    strerr = cudaGetErrorString(err);
  }
#endif

  if (error) {
    fprintf(stderr,"Error after %s: %s\n", kernel_name, strerr);

    // outputs error file
    FILE *fp;
    int myrank;
    char filename[BUFSIZ];
#ifdef WITH_MPI
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
#else
    myrank = 0;
#endif
    sprintf(filename,"OUTPUT_FILES/error_message_%06d.txt",myrank);
    fp = fopen (filename, "a+");
    if (fp != NULL) {
      fprintf (fp, "Error after %s: %s\n", kernel_name, strerr);
      fclose (fp);
    }

    // stops program
#ifdef WITH_MPI
    MPI_Abort (MPI_COMM_WORLD, 1);
#endif
    exit (EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------------------------------- */

void exit_on_error (char *info) {
  printf ("\nERROR: %s\n", info);
  fflush (stdout);

  // outputs error file
  FILE *fp;
  int myrank;
  char filename[BUFSIZ];
#ifdef WITH_MPI
  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
#else
  myrank = 0;
#endif
  sprintf(filename,"OUTPUT_FILES/error_message_%06d.txt",myrank);
  fp = fopen (filename, "a+");
  if (fp != NULL) {
    fprintf (fp, "Error: %s\n", info);
    fclose (fp);
  }

  // stops program
#ifdef WITH_MPI
  MPI_Abort (MPI_COMM_WORLD, 1);
#endif
  exit (EXIT_FAILURE);
}


/*----------------------------------------------------------------------------------------------- */
// additional helper functions
/*----------------------------------------------------------------------------------------------- */

double get_time_val () {
  struct timeval t;
  struct timezone tzp;
  gettimeofday (&t, &tzp);
  return t.tv_sec + t.tv_usec*1e-6;
}

/*----------------------------------------------------------------------------------------------- */

void pause_for_debugger (int pause) {
  if (pause) {
    int myrank;
#ifdef WITH_MPI
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
#else
    myrank = 0;
#endif
    printf ("I'm rank %d\n", myrank);
    int i = 0;
    char hostname[256];
    gethostname (hostname, sizeof (hostname));
    printf ("PID %d on %s:%d ready for attach\n", getpid (), hostname, myrank);

    FILE *file = fopen ("./attach_gdb.txt", "w+");
    if (file != NULL) {
      fprintf (file, "PID %d on %s:%d ready for attach\n", getpid (), hostname, myrank);
      fclose (file);
    }

    fflush (stdout);
    while (0 == i)
      sleep (5);
  }
}

/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (pause_for_debug,
               PAUSE_FOR_DEBUG) () {
  TRACE ("pause_for_debug");

  pause_for_debugger (1);
}

/* ----------------------------------------------------------------------------------------------- */
// MPI synchronization
/* ----------------------------------------------------------------------------------------------- */

void synchronize_mpi () {
#ifdef WITH_MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// kernel setup functions
/* ----------------------------------------------------------------------------------------------- */

void get_blocks_xy (int num_blocks, int *num_blocks_x, int *num_blocks_y) {
  // Initially sets the blocks_x to be the num_blocks, and adds rows as needed (block size limit of 65535).
  // If an additional row is added, the row length is cut in
  // half. If the block count is odd, there will be 1 too many blocks,
  // which must be managed at runtime with an if statement.

  *num_blocks_x = num_blocks;
  *num_blocks_y = 1;

  while (*num_blocks_x > MAXIMUM_GRID_DIM) {
    *num_blocks_x = (int) ceil (*num_blocks_x * 0.5f);
    *num_blocks_y = *num_blocks_y * 2;
  }
}
