/*
!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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

/*----------------------------------------------------------------------------------------------- */
// Helper functions
/*----------------------------------------------------------------------------------------------- */

double get_time () {
  struct timeval t;
  struct timezone tzp;
  gettimeofday (&t, &tzp);
  return t.tv_sec + t.tv_usec*1e-6;
}

/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (pause_for_debug,
               PAUSE_FOR_DEBUG) () {
  TRACE ("pause_for_debug");

  pause_for_debugger (1);
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
#endif

void exit_on_gpu_error (char *kernel_name) {
  //check to catch errors from previous operations
  // /!\ in opencl, we can't have information about the last ASYNC error
  int error = 0;
  const char *strerr = NULL;

#ifdef USE_OPENCL
  if (run_opencl) {
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

    // releases previous contexts

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
    fprintf (fp, "ERROR: %s\n", info);
    fclose (fp);
  }

  // stops program
#ifdef WITH_MPI
  MPI_Abort (MPI_COMM_WORLD, 1);
#endif
  exit (EXIT_FAILURE);
}

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
    sprintf(filename,"../in_out_files/OUTPUT_FILES/error_message_%06d.txt",myrank);
    fp = fopen(filename,"a+");
    if (fp != NULL){
      fprintf(fp,"\nCUDA error !!!!! <%s> !!!!! \nat CUDA call error code: # %d\n",cudaGetErrorString(err),num);
      fclose(fp);
    }

    // stops program
#ifdef WITH_MPI
    MPI_Abort(MPI_COMM_WORLD,1);
#endif
    exit(EXIT_FAILURE);
  }
  return;
}

/* ----------------------------------------------------------------------------------------------- */
// CUDA synchronization
/* ----------------------------------------------------------------------------------------------- */

void synchronize_cuda() {
#if CUDA_VERSION >= 4000
    cudaDeviceSynchronize();
#else
    cudaThreadSynchronize();
#endif
}

#endif

/*----------------------------------------------------------------------------------------------- */

void synchronize_mpi () {
#ifdef WITH_MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif
}

#ifdef USE_CUDA
/*----------------------------------------------------------------------------------------------- */
// Timing helper functions
/* ----------------------------------------------------------------------------------------------- */


void start_timing_cuda(cudaEvent_t* start,cudaEvent_t* stop){
  // creates & starts event
  cudaEventCreate(start);
  cudaEventCreate(stop);
  cudaEventRecord( *start, 0 );
}

/* ----------------------------------------------------------------------------------------------- */

void stop_timing_cuda(cudaEvent_t* start,cudaEvent_t* stop, char* info_str){
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

#endif

/* ----------------------------------------------------------------------------------------------- */
// GPU kernel setup functions
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

/*----------------------------------------------------------------------------------------------- */
// GPU device memory functions
/* ----------------------------------------------------------------------------------------------- */

void get_free_memory (double *free_db, double *used_db, double *total_db) {
#ifdef USE_OPENCL
  if (run_opencl) {
    *free_db = 0;
    *total_db = 0;
    *used_db = 0;
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    // gets memory usage in byte
    size_t free_byte ;
    size_t total_byte ;
    cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;
    if ( cudaSuccess != cuda_status ){
      printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
      exit(EXIT_FAILURE);
    }

    *free_db = (double)free_byte ;
    *total_db = (double)total_byte ;
    *used_db = *total_db - *free_db ;
  }
#endif
}

/*----------------------------------------------------------------------------------------------- */
// Saves GPU memory usage to file

void output_free_memory (int myrank, char *info_str) {
  FILE *fp;
  char filename[BUFSIZ];
  double free_db, used_db, total_db;
  int do_output_info;

  // by default, only master process outputs device info to avoid file cluttering
  do_output_info = 0;
  if( myrank == 0 ){
    do_output_info = 1;
    sprintf(filename,"OUTPUT_FILES/gpu_device_mem_usage.txt");
  }
  // debugging
  if( DEBUG ){
    do_output_info = 1;
    sprintf(filename,"OUTPUT_FILES/gpu_device_mem_usage_proc_%06d.txt",myrank);
  }

  // outputs to file
  if( do_output_info ){

    // gets memory usage
  get_free_memory (&free_db, &used_db, &total_db);

    // file output
  fp = fopen (filename, "a+");
  if (fp != NULL) {
      fprintf(fp,"%d: @%s GPU memory usage: used = %f MB, free = %f MB, total = %f MB\n", myrank, info_str,
              used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);
    fclose (fp);
    }
  }
}

/*----------------------------------------------------------------------------------------------- */

// Fortran-callable version of above method
extern EXTERN_LANG
void FC_FUNC_ (output_free_device_memory,
               OUTPUT_FREE_DEVICE_MEMORY) (int *myrank) {

  TRACE ("output_free_device_memory");

  char info[6];
  sprintf (info, "f %d:", *myrank);
  output_free_memory (*myrank, info);
}

/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (get_free_device_memory,
               get_FREE_DEVICE_MEMORY) (realw *free, realw *used, realw *total) {
  TRACE ("get_free_device_memory");

  double free_db, used_db, total_db;

  get_free_memory (&free_db, &used_db, &total_db);

  // converts to MB
  *free = (realw) free_db / 1024.0 / 1024.0;
  *used = (realw) used_db / 1024.0 / 1024.0;
  *total = (realw) total_db / 1024.0 / 1024.0;
}


/*----------------------------------------------------------------------------------------------- */
// Auxiliary functions
/*----------------------------------------------------------------------------------------------- */

realw get_device_array_maximum_value (Mesh *mp, gpu_realw_mem *d_array, int size) {
  realw max = 0.0f;

  // checks if anything to do
  if (size > 0) {
    realw *h_array = (realw *) calloc(size ,sizeof (realw));

    h_array = (realw *) calloc (size, sizeof (realw));
#ifdef USE_OPENCL
    if (run_opencl) {
      clCheck (clEnqueueReadBuffer (mocl.command_queue, d_array->ocl, CL_TRUE, 0,
                                    sizeof (realw) * size,
                                    h_array, 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      // explicitly wait for cuda kernels to finish
      // (cudaMemcpy implicitly synchronizes all other cuda operations)
      synchronize_cuda();
      print_CUDA_error_if_any(cudaMemcpy(h_array,d_array,sizeof(realw)*size,cudaMemcpyDeviceToHost),33001);
    }
#endif
    // finds maximum value in array
    max = h_array[0];
    int i;
    for (i = 1; i < size; i++) {
      if (abs (h_array[i]) > max)
        max = abs (h_array[i]);
    }
    free (h_array);
  }
  return max;
}


/*-----------------------------------------------------------------------------------------------*/
// scalar arrays (acoustic/fluid outer core)
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (check_norm_acoustic_from_device,
               CHECK_NORM_ACOUSTIC_FROM_DEVICE) (realw *norm,
                                                 long *Mesh_pointer_f,
                                                 int *FORWARD_OR_ADJOINT) {
  TRACE ("check_norm_acoustic_from_device");

  Mesh *mp = (Mesh *) *Mesh_pointer_f;     //get mesh pointer out of Fortran integer container
  realw max;
  gpu_realw_mem d_max;

  max = 0.0f;

  // way 2 b: timing Elapsed time: 1.236916e-03
  // launch simple reduction kernel
  realw *h_max;
  int blocksize = BLOCKSIZE_TRANSFER;

  // outer core
  int size = mp->NGLOB_OUTER_CORE;

  int size_padded = ((int) ceil (((double) size)/ ((double) blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

  h_max = (realw *) calloc (num_blocks_x * num_blocks_y, sizeof (realw));

#ifdef USE_OPENCL
  if (run_opencl) {
    cl_int errcode;
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    d_max.ocl = clCreateBuffer (mocl.context, CL_MEM_READ_WRITE, num_blocks_x * num_blocks_y * sizeof (realw), NULL, &errcode);

    if (*FORWARD_OR_ADJOINT == 1) {
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_displ_outer_core.ocl));
    } else if (*FORWARD_OR_ADJOINT == 3) {
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_displ_outer_core.ocl));
    }

    clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (int), (void *) &size));
    clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (cl_mem), (void *) &d_max.ocl));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.get_maximum_scalar_kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer (mocl.command_queue, d_max.ocl, CL_TRUE, 0,
                                  num_blocks_x * num_blocks_y * sizeof (realw),
                                  h_max, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    cudaMalloc((void**)&d_max.cuda,num_blocks_x*num_blocks_y*sizeof(realw));

    if(*FORWARD_OR_ADJOINT == 1 ){
      get_maximum_scalar_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_displ_outer_core.cuda,size,d_max.cuda);
    }else if(*FORWARD_OR_ADJOINT == 3 ){
      get_maximum_scalar_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_displ_outer_core.cuda,size,d_max.cuda);
    }

    print_CUDA_error_if_any(cudaMemcpy(h_max,d_max.cuda,num_blocks_x*num_blocks_y*sizeof(realw),
                                       cudaMemcpyDeviceToHost),222);
  }
#endif
  // determines max for all blocks
  max = h_max[0];
  int i;
  for (i = 1; i < num_blocks_x * num_blocks_y; i++) {
    if (max < h_max[i])
      max = h_max[i];
  }
#ifdef USE_OPENCL
  if (run_opencl) {
    clReleaseMemObject (d_max.ocl);
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    cudaFree(d_max.cuda);
  }
#endif
  free (h_max);

  // return result
  *norm = max;

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_gpu_error ("after check_norm_acoustic_from_device");
#endif
}

extern EXTERN_LANG
void FC_FUNC_ (check_norm_elastic_from_device,
               CHECK_NORM_ELASTIC_FROM_DEVICE) (realw *norm,
                                                long *Mesh_pointer_f,
                                                int *FORWARD_OR_ADJOINT) {

  TRACE ("check_norm_elastic_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) (*Mesh_pointer_f);

  realw max, max_crust_mantle, max_inner_core;

  int size, size_padded;

  // launch simple reduction kernel
  realw *h_max;
  int blocksize = BLOCKSIZE_TRANSFER;

  // crust_mantle
  max = 0.0f;
  size = mp->NGLOB_CRUST_MANTLE;

  size_padded = ((int) ceil (((double) size) / ((double) blocksize))) * blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (size_padded / blocksize, &num_blocks_x, &num_blocks_y);

  h_max = (realw *) calloc (num_blocks_x * num_blocks_y, sizeof (realw));

  gpu_realw_mem d_max;
#ifdef USE_OPENCL
  cl_int errcode;

  size_t global_work_size[2];
  size_t local_work_size[2];
  cl_uint idx = 0;

  if (run_opencl) {
    d_max.ocl = clCreateBuffer (mocl.context, CL_MEM_READ_WRITE, num_blocks_x*num_blocks_y*sizeof (realw), NULL, &errcode);

    if (*FORWARD_OR_ADJOINT == 1) {
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_displ_crust_mantle.ocl));
    } else if (*FORWARD_OR_ADJOINT == 3) {
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_displ_crust_mantle.ocl));
    }
    clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (int), (void *) &size));
    clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (cl_mem), (void *) &d_max.ocl));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.get_maximum_vector_kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL));

    // copies to CPU
    clCheck (clEnqueueReadBuffer (mocl.command_queue, d_max.ocl, CL_TRUE, 0,
                                  num_blocks_x * num_blocks_y * sizeof (realw),
                                  h_max, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  dim3 grid,threads;

  if (run_cuda) {
    grid = dim3(num_blocks_x,num_blocks_y);
    threads = dim3(blocksize,1,1);

    cudaMalloc((void**)&d_max.cuda,num_blocks_x*num_blocks_y*sizeof(realw));

    if(*FORWARD_OR_ADJOINT == 1 ){
      get_maximum_vector_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_displ_crust_mantle.cuda,size,d_max.cuda);
    }else if(*FORWARD_OR_ADJOINT == 3 ){
      get_maximum_vector_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_displ_crust_mantle.cuda,size,d_max.cuda);
    }
  }
#endif

  // determines max for all blocks
  max = h_max[0];
  int i;
  for (i = 1; i < num_blocks_x * num_blocks_y; i++) {
    if (max < h_max[i])
      max = h_max[i];
  }
  max_crust_mantle = max;

#ifdef USE_OPENCL
  if (run_opencl) {
    clReleaseMemObject (d_max.ocl);
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    cudaFree(d_max.cuda);
  }
#endif
  free (h_max);

  // inner_core
  max = 0.0f;
  size = mp->NGLOB_INNER_CORE;

  size_padded = ((int) ceil (((double) size) / ((double) blocksize))) * blocksize;

  get_blocks_xy (size_padded / blocksize, &num_blocks_x, &num_blocks_y);


  h_max = (realw *) calloc (num_blocks_x * num_blocks_y, sizeof (realw));

#ifdef USE_OPENCL
  if (run_opencl) {
    idx = 0;
    d_max.ocl = clCreateBuffer (mocl.context, CL_MEM_READ_WRITE, num_blocks_x * num_blocks_y * sizeof (realw), NULL, &errcode);

    if (*FORWARD_OR_ADJOINT == 1) {
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_displ_inner_core.ocl));
    } else if (*FORWARD_OR_ADJOINT == 3) {
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_displ_inner_core.ocl));
    } else {
      goto skip_exec;
    }

    clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (int), (void *) &size));
    clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (cl_mem), (void *) &d_max.ocl));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.get_maximum_vector_kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL));

  skip_exec:
    // copies to CPU
    clCheck (clEnqueueReadBuffer (mocl.command_queue, d_max.ocl, CL_TRUE, 0,
                                  num_blocks_x * num_blocks_y * sizeof (realw),
                                  h_max, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    grid = dim3(num_blocks_x,num_blocks_y);
    threads = dim3(blocksize,1,1);

    cudaMalloc((void**)&d_max.cuda,num_blocks_x*num_blocks_y*sizeof(realw));

    if(*FORWARD_OR_ADJOINT == 1 ){
      get_maximum_vector_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_displ_inner_core.cuda,size,d_max.cuda);
    }else if(*FORWARD_OR_ADJOINT == 3 ){
      get_maximum_vector_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_displ_inner_core.cuda,size,d_max.cuda);
    }

    // copies to CPU
    print_CUDA_error_if_any(cudaMemcpy(h_max,d_max.cuda,num_blocks_x*num_blocks_y*sizeof(realw),
                                       cudaMemcpyDeviceToHost),222);
  }
#endif
  // determines max for all blocks
  max = h_max[0];
  for (i = 1; i < num_blocks_x * num_blocks_y; i++) {
    if (max < h_max[i])
      max = h_max[i];
  }
  max_inner_core = max;

#ifdef USE_OPENCL
  if (run_opencl) {
    clReleaseMemObject (d_max.ocl);
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    cudaFree(d_max.cuda);
  }
#endif

  free (h_max);

  // return result
  max = MAX (max_inner_core, max_crust_mantle);
  *norm = max;

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_gpu_error ("after check_norm_elastic_from_device");
#endif
}

/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (check_norm_strain_from_device,
               CHECK_NORM_STRAIN_FROM_DEVICE) (realw *strain_norm,
                                               realw *strain_norm2,
                                               long *Mesh_pointer_f) {

  TRACE ("check_norm_strain_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) (*Mesh_pointer_f);

  realw max, max_eps;
  gpu_realw_mem d_max;
  int num_blocks_x, num_blocks_y;
  int size, size_padded;
#ifdef USE_OPENCL
  cl_int errcode;

  size_t local_work_size[2];
  size_t global_work_size[2];
  cl_uint idx = 0;
#endif
#ifdef USE_CUDA
  dim3 grid,threads;
#endif
  // launch simple reduction kernel
  realw *h_max;
  int blocksize = BLOCKSIZE_TRANSFER;

  // crust_mantle strain arrays
  size = NGLL3 * (mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY);

  size_padded = ((int) ceil (((double) size) / ((double) blocksize))) * blocksize;

  get_blocks_xy (size_padded / blocksize, &num_blocks_x, &num_blocks_y);


  h_max = (realw *) calloc (num_blocks_x * num_blocks_y, sizeof (realw));

  max = 0.0f;

#ifdef USE_OPENCL
  if (run_opencl) {
    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    idx = 0;
    d_max.ocl = clCreateBuffer (mocl.context, CL_MEM_READ_WRITE, num_blocks_x * num_blocks_y * sizeof (realw), NULL, clck_(&errcode));

    // determines max for: eps_trace_over_3_crust_mantle
    clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_eps_trace_over_3_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (int), (void *) &size));
    clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (cl_mem), (void *) &d_max.ocl));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;
    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.get_maximum_scalar_kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer (mocl.command_queue, d_max.ocl, CL_TRUE, 0,
                                  num_blocks_x * num_blocks_y * sizeof (realw),
                                  h_max, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    grid = dim3(num_blocks_x,num_blocks_y);
    threads = dim3(blocksize,1,1);

    cudaMalloc((void**)&d_max.cuda,num_blocks_x*num_blocks_y*sizeof(realw));

    // determines max for: eps_trace_over_3_crust_mantle
    get_maximum_scalar_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_eps_trace_over_3_crust_mantle.cuda,size,d_max.cuda);

    print_CUDA_error_if_any(cudaMemcpy(h_max,d_max.cuda,num_blocks_x*num_blocks_y*sizeof(realw),
                                       cudaMemcpyDeviceToHost),221);
  }
#endif
  max = h_max[0];
  int i;
  for (i = 1; i < num_blocks_x * num_blocks_y; i++) {
    if (max < h_max[i])
      max = h_max[i];
  }
  // strain trace maximum
  *strain_norm = max;

  // frees arrays
#ifdef USE_OPENCL
  if (run_opencl) {
    clReleaseMemObject (d_max.ocl);
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    cudaFree(d_max.cuda);
  }
#endif

  free (h_max);

  // initializes
  // crust_mantle arrays
  size = NGLL3 * mp->NSPEC_CRUST_MANTLE;

  size_padded = ((int) ceil (((double) size) / ((double) blocksize))) * blocksize;

  get_blocks_xy (size_padded / blocksize, &num_blocks_x, &num_blocks_y);


  h_max = (realw *) calloc (num_blocks_x*num_blocks_y, sizeof (realw));
  max_eps = 0.0f;

#ifdef USE_OPENCL
  if (run_opencl) {
    d_max.ocl = clCreateBuffer (mocl.context, CL_MEM_READ_WRITE, num_blocks_x * num_blocks_y * sizeof (realw), NULL, clck_(&errcode));

    int loop;
    gpu_realw_mem d_epsilondev_HH_crust_mantle[] =
      {mp->d_epsilondev_xx_crust_mantle,
       mp->d_epsilondev_yy_crust_mantle,
       mp->d_epsilondev_xy_crust_mantle,
       mp->d_epsilondev_xz_crust_mantle,
       mp->d_epsilondev_yz_crust_mantle
      };
    for (loop = 0; loop < 5; loop++) {
      idx = 0;
      // determines max for: epsilondev_xx_crust_mantle
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (cl_mem), (void *) &d_epsilondev_HH_crust_mantle[loop].ocl));
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (int), (void *) &size));
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (cl_mem), (void *) &d_max.ocl));

      clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.get_maximum_scalar_kernel, 2, NULL,
                                       global_work_size, local_work_size, 0, NULL, NULL));

      clCheck (clEnqueueReadBuffer (mocl.command_queue, d_max.ocl, CL_TRUE, 0,
                                    num_blocks_x * num_blocks_y * sizeof (realw),
                                    h_max, 0, NULL, NULL));

      max = h_max[0];
      for (i = 1; i < num_blocks_x * num_blocks_y; i++) {
        if (max < h_max[i])
          max = h_max[i];
      }
      max_eps = MAX (max_eps, max);
    }
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    cudaMalloc((void**)&d_max.cuda,num_blocks_x*num_blocks_y*sizeof(realw));

    // determines max for: epsilondev_xx_crust_mantle
    get_maximum_scalar_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_epsilondev_xx_crust_mantle.cuda,size,d_max.cuda);

    print_CUDA_error_if_any(cudaMemcpy(h_max,d_max.cuda,num_blocks_x*num_blocks_y*sizeof(realw),
                                       cudaMemcpyDeviceToHost),222);
    max = h_max[0];
    for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
      if( max < h_max[i]) max = h_max[i];
    }
    max_eps = MAX(max_eps,max);

    // determines max for: epsilondev_yy_crust_mantle
    get_maximum_scalar_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_epsilondev_yy_crust_mantle.cuda,size,d_max.cuda);

    print_CUDA_error_if_any(cudaMemcpy(h_max,d_max.cuda,num_blocks_x*num_blocks_y*sizeof(realw),
                                       cudaMemcpyDeviceToHost),223);
    max = h_max[0];
    for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
      if( max < h_max[i]) max = h_max[i];
    }
    max_eps = MAX(max_eps,max);

    // determines max for: epsilondev_xy_crust_mantle
    get_maximum_scalar_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_epsilondev_xy_crust_mantle.cuda,size,d_max.cuda);

    print_CUDA_error_if_any(cudaMemcpy(h_max,d_max.cuda,num_blocks_x*num_blocks_y*sizeof(realw),
                                       cudaMemcpyDeviceToHost),224);
    max = h_max[0];
    for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
      if( max < h_max[i]) max = h_max[i];
    }
    max_eps = MAX(max_eps,max);

    // determines max for: epsilondev_xz_crust_mantle
    get_maximum_scalar_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_epsilondev_xz_crust_mantle.cuda,size,d_max.cuda);

    print_CUDA_error_if_any(cudaMemcpy(h_max,d_max.cuda,num_blocks_x*num_blocks_y*sizeof(realw),
                                       cudaMemcpyDeviceToHost),225);
    max = h_max[0];
    for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
      if( max < h_max[i]) max = h_max[i];
    }
    max_eps = MAX(max_eps,max);

    // determines max for: epsilondev_yz_crust_mantle
    get_maximum_scalar_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_epsilondev_yz_crust_mantle.cuda,size,d_max.cuda);

    print_CUDA_error_if_any(cudaMemcpy(h_max,d_max.cuda,num_blocks_x*num_blocks_y*sizeof(realw),
                                       cudaMemcpyDeviceToHost),226);
    max = h_max[0];
    for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
      if( max < h_max[i]) max = h_max[i];
    }
    max_eps = MAX(max_eps,max);
  }
#endif
  // strain maximum
  *strain_norm2 = max_eps;

  // frees arrays
#ifdef USE_OPENCL
  if (run_opencl) {
    clReleaseMemObject (d_max.ocl);
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    cudaFree(d_max.cuda);
  }
#endif

  free (h_max);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_gpu_error ("after check_norm_strain_from_device");
#endif
}
