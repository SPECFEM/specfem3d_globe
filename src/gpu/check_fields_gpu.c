/*
!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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
    if (cudaSuccess != cuda_status) {
      printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
      exit(EXIT_FAILURE);
    }

    *free_db = (double)free_byte ;
    *total_db = (double)total_byte ;
    *used_db = *total_db - *free_db ;
  }
#endif
#ifdef USE_HIP
  if (run_hip) {
    // gets memory usage in byte
    size_t free_byte ;
    size_t total_byte ;
    hipError_t hip_status = hipMemGetInfo( &free_byte, &total_byte ) ;
    if (hipSuccess != hip_status) {
      printf("Error: hipMemGetInfo fails, %s \n", hipGetErrorString(hip_status) );
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

  // by default, only main process outputs device info to avoid file cluttering
  do_output_info = 0;
  if (myrank == 0) {
    do_output_info = 1;
    sprintf(filename,"OUTPUT_FILES/gpu_device_mem_usage.txt");
  }
  // debugging
  if (DEBUG) {
    do_output_info = 1;
    sprintf(filename,"OUTPUT_FILES/gpu_device_mem_usage_proc_%06d.txt",myrank);
  }

  // outputs to file
  if (do_output_info) {

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

  char info_str[15]; // add extra character for null termination
  int len;

  // safety check to avoid string buffer overflow
  if (*myrank > 99999999) { exit_on_error("Error: rank too large in output_free_device_memory() routine"); }

  len = snprintf (info_str, 15, "rank %8d:", *myrank);
  if (len >= 15){ printf("warning: string length truncated (from %d) in output_free_device_memory() routine\n", len); }

  //debug
  //printf("debug: info ***%s***\n",info_str);

  // writes to output file
  output_free_memory (*myrank, info_str);
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

realw get_device_array_maximum_value (gpu_realw_mem d_array, int size) {

// gets maximum of array on GPU by copying over to CPU and handle it there

  realw *h_array;
  realw max = 0.0f;

  // checks if anything to do
  if (size > 0) {

    // explicitly wait for gpu kernels to finish
    gpuSynchronize();

    h_array = (realw *) calloc (size, sizeof (realw));
    if (h_array == NULL) { exit_on_error("Error allocating h_array array in get_device_array_maximum_value() routine"); }

    // copies values from gpu to cpu array
    gpuCopy_from_device_realw (&d_array, h_array, size);

    // finds maximum value in array
    max = fabs(h_array[0]);
    for (int i = 1; i < size; i++) {
      if (fabs(h_array[i]) > max) max = fabs(h_array[i]);
    }

    // frees temporary array
    free (h_array);
  }
  return max;
}


/*-----------------------------------------------------------------------------------------------*/
// scalar arrays (acoustic/fluid outer core) and vector arrays elastic
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (check_norm_elastic_acoustic_from_device,
               CHECK_NORM_ELASTIC_ACOUSTIC_FROM_DEVICE) (realw *solidnorm,realw *fluidnorm,
                                                         long *Mesh_pointer_f,
                                                         int *FORWARD_OR_ADJOINT) {

  TRACE ("check_norm_elastic_acoustic_from_device");

  // get Mesh from Fortran integer wrapper
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int num_blocks_x, num_blocks_y;
  int size_nonpadded, size_padded;

  // initializes
  *fluidnorm = 0.f;
  *solidnorm = 0.f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in check_norm_elastic_acoustic_from_device() routine");
  }

  //debug
#if (DEBUG_FIELDS == 1)
  {
    printf ("rank %d - check norm max: nglob outer_core %d crust_mantle %d inner_core %d\n",
            mp->myrank, mp->NGLOB_OUTER_CORE,mp->NGLOB_CRUST_MANTLE,mp->NGLOB_INNER_CORE);
    fflush (stdout);
    synchronize_mpi ();
    realw max_d, max_v, max_a;
    max_d = get_device_array_maximum_value(mp->d_displ_outer_core, mp->NGLOB_OUTER_CORE);
    max_v = get_device_array_maximum_value(mp->d_veloc_outer_core, mp->NGLOB_OUTER_CORE);
    max_a = get_device_array_maximum_value(mp->d_accel_outer_core, mp->NGLOB_OUTER_CORE);
    printf ("rank %d - check norm max outer_core displ: %e veloc: %e accel: %e, %i\n", mp->myrank, max_d, max_v, max_a, *FORWARD_OR_ADJOINT);
    fflush (stdout);
    max_d = get_device_array_maximum_value(mp->d_displ_crust_mantle, NDIM * mp->NGLOB_CRUST_MANTLE);
    max_v = get_device_array_maximum_value(mp->d_veloc_crust_mantle, NDIM * mp->NGLOB_CRUST_MANTLE);
    max_a = get_device_array_maximum_value(mp->d_accel_crust_mantle, NDIM * mp->NGLOB_CRUST_MANTLE);
    printf ("rank %d - check norm max crust_mantle displ: %e veloc: %e accel: %e, %i\n", mp->myrank, max_d, max_v, max_a, *FORWARD_OR_ADJOINT);
    fflush (stdout);
    max_d = get_device_array_maximum_value(mp->d_displ_inner_core, NDIM * mp->NGLOB_INNER_CORE);
    max_v = get_device_array_maximum_value(mp->d_veloc_inner_core, NDIM * mp->NGLOB_INNER_CORE);
    max_a = get_device_array_maximum_value(mp->d_accel_inner_core, NDIM * mp->NGLOB_INNER_CORE);
    printf ("rank %d - check norm max inner_core displ: %e veloc: %e accel: %e, %i\n", mp->myrank, max_d, max_v, max_a, *FORWARD_OR_ADJOINT);
    fflush (stdout);
    synchronize_mpi ();
  }
#endif

  // launch simple reduction kernel
  int blocksize = BLOCKSIZE_TRANSFER;
  gpu_realw_mem displ;

  // outer core
  size_nonpadded = mp->NGLOB_OUTER_CORE;

  size_padded = ((int) ceil (((double) size_nonpadded)/ ((double) blocksize)))*blocksize;
  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

  int size_block_fluid = num_blocks_x * num_blocks_y;
  //realw *h_max_fluid = (realw *) calloc (size_block_fluid, sizeof (realw));
  //if (h_max_fluid == NULL) { exit_on_error("Error allocating h_max array in check_norm_elastic_acoustic_from_device() routine"); }

  // cuda graph
#ifdef USE_CUDA_GRAPHS
  if (mp->init_graph_norm){
    // debug: synchronizes first
    //gpuSynchronize(); synchronize_mpi();
    // start capturing
    print_CUDA_error_if_any(cudaStreamBeginCapture(mp->compute_stream),930);
  }
#endif

  // outer core maximum
  if (mp->NGLOB_OUTER_CORE > 0) {
    // sets gpu arrays
    if (*FORWARD_OR_ADJOINT == 1) {
      displ = mp->d_displ_outer_core;
    } else {
      displ = mp->d_b_displ_outer_core;
    }

#ifdef USE_OPENCL
    if (run_opencl) {
      size_t global_work_size[2];
      size_t local_work_size[2];
      cl_uint idx = 0;

      local_work_size[0] = blocksize;
      local_work_size[1] = 1;
      global_work_size[0] = num_blocks_x * blocksize;
      global_work_size[1] = num_blocks_y;

      clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (cl_mem), (void *) &displ.ocl));
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (int), (void *) &size_nonpadded));
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_norm_max.ocl));

      clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.get_maximum_scalar_kernel, 2, NULL,
                                       global_work_size, local_work_size, 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      // graph
#ifdef USE_CUDA_GRAPHS
      if (! mp->use_graph_call_norm){
#endif

      dim3 grid(num_blocks_x,num_blocks_y);
      dim3 threads(blocksize,1,1);

      get_maximum_scalar_kernel<<<grid,threads,0,mp->compute_stream>>>(displ.cuda,size_nonpadded,mp->d_norm_max.cuda);

#ifdef USE_CUDA_GRAPHS
      } // graph
#endif
    }
#endif
#ifdef USE_HIP
    if (run_hip) {
      dim3 grid(num_blocks_x,num_blocks_y);
      dim3 threads(blocksize,1,1);

      hipLaunchKernelGGL(HIP_KERNEL_NAME(get_maximum_scalar_kernel), grid, threads, 0, mp->compute_stream,
                                                                     displ.hip,size_nonpadded,mp->d_norm_max.hip);
    }
#endif
  } //outer_core

  // crust_mantle
  size_nonpadded = mp->NGLOB_CRUST_MANTLE;

  size_padded = ((int) ceil (((double) size_nonpadded) / ((double) blocksize))) * blocksize;
  get_blocks_xy (size_padded / blocksize, &num_blocks_x, &num_blocks_y);

  int size_block_cm = num_blocks_x * num_blocks_y;
  //realw *h_max_cm = (realw *) calloc (size_block_cm, sizeof (realw));
  //if (h_max_cm == NULL) { exit_on_error("Error allocating h_max array in check_norm_elastic_acoustic_from_device() routine"); }

  // sets gpu arrays
  if (*FORWARD_OR_ADJOINT == 1) {
    displ = mp->d_displ_crust_mantle;
  } else {
    displ = mp->d_b_displ_crust_mantle;
  }

  gpu_realw_mem tmp;
  int offset1,offset2;

  INITIALIZE_OFFSET();

  offset1 = size_block_fluid;
  INIT_OFFSET(d_norm_max, offset1);

  tmp = gpuTakeRef(PASS_OFFSET(d_norm_max, offset1));

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (cl_mem), (void *) &displ.ocl));
    clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (int), (void *) &size_nonpadded));
    clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (cl_mem), (void *) &tmp.ocl));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.get_maximum_vector_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    // graph
#ifdef USE_CUDA_GRAPHS
    if (! mp->use_graph_call_norm){
#endif

    dim3 grid = dim3(num_blocks_x,num_blocks_y);
    dim3 threads = dim3(blocksize,1,1);

    get_maximum_vector_kernel<<<grid,threads,0,mp->compute_stream>>>(displ.cuda,size_nonpadded,tmp.cuda);

#ifdef USE_CUDA_GRAPHS
    } // graph
#endif
  }
#endif
#ifdef USE_HIP
  if (run_hip) {
    dim3 grid = dim3(num_blocks_x,num_blocks_y);
    dim3 threads = dim3(blocksize,1,1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(get_maximum_vector_kernel), grid, threads, 0, mp->compute_stream,
                                                                   displ.hip,size_nonpadded,tmp.hip);
  }
#endif

  RELEASE_OFFSET(d_norm_max, offset1);

  // inner_core
  size_nonpadded = mp->NGLOB_INNER_CORE;

  size_padded = ((int) ceil (((double) size_nonpadded) / ((double) blocksize))) * blocksize;
  get_blocks_xy (size_padded / blocksize, &num_blocks_x, &num_blocks_y);

  int size_block_ic = num_blocks_x * num_blocks_y;
  //realw *h_max_ic = (realw *) calloc (size_block_ic, sizeof (realw));
  //if (h_max_ic == NULL) { exit_on_error("Error allocating h_max array for inner core in check_norm_elastic_acoustic_from_device() routine"); }

  if (mp->NGLOB_INNER_CORE > 0) {
    // sets gpu arrays
    if (*FORWARD_OR_ADJOINT == 1) {
      displ = mp->d_displ_inner_core;
    } else {
      displ = mp->d_b_displ_inner_core;
    }

    offset2 = size_block_fluid + size_block_cm;
    INIT_OFFSET(d_norm_max, offset2);

    tmp = gpuTakeRef(PASS_OFFSET(d_norm_max, offset2));

#ifdef USE_OPENCL
    if (run_opencl) {
      size_t global_work_size[2];
      size_t local_work_size[2];
      cl_uint idx = 0;

      local_work_size[0] = blocksize;
      local_work_size[1] = 1;
      global_work_size[0] = num_blocks_x * blocksize;
      global_work_size[1] = num_blocks_y;

      clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (cl_mem), (void *) &displ.ocl));
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (int), (void *) &size_nonpadded));
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (cl_mem), (void *) &tmp.ocl));

      clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.get_maximum_vector_kernel, 2, NULL,
                                       global_work_size, local_work_size, 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      // graph
#ifdef USE_CUDA_GRAPHS
      if (! mp->use_graph_call_norm){
#endif

      dim3 grid = dim3(num_blocks_x,num_blocks_y);
      dim3 threads = dim3(blocksize,1,1);

      get_maximum_vector_kernel<<<grid,threads,0,mp->compute_stream>>>(displ.cuda,size_nonpadded,tmp.cuda);

#ifdef USE_CUDA_GRAPHS
      } // graph
#endif
    }
#endif
#ifdef USE_HIP
    if (run_hip) {
      dim3 grid = dim3(num_blocks_x,num_blocks_y);
      dim3 threads = dim3(blocksize,1,1);

      hipLaunchKernelGGL(HIP_KERNEL_NAME(get_maximum_vector_kernel), grid, threads, 0, mp->compute_stream,
                                                                     displ.hip,size_nonpadded,tmp.hip);
    }
#endif

    RELEASE_OFFSET(d_norm_max, offset2);
  } //inner_core

  // graph
#ifdef USE_CUDA_GRAPHS
  if (! mp->use_graph_call_norm){
#endif

  // note: we use event synchronization to make sure that compute stream has finished before copying data to host.
  //       the same for continuing only after copy to host has finished.
  //
  //       capturing graphs seems only to work with this sort of event synchronization...

  // create event for synchronizing with async copy
  gpuRecordEvent(mp);

  // copies to CPU
  gpuCopy_from_device_realw_asyncEvent(mp, &mp->d_norm_max, mp->h_norm_max, size_block_fluid + size_block_cm + size_block_ic);

  // makes sure copy has finished
  gpuWaitEvent(mp);

#ifdef USE_CUDA_GRAPHS
  } // graph

  // finish creating graph
  if (mp->init_graph_norm){
    // stop capturing
    print_CUDA_error_if_any(cudaStreamEndCapture(mp->compute_stream, &mp->graph_norm),930);

    // get graph info
    size_t numNodes = 0;
    print_CUDA_error_if_any(cudaGraphGetNodes(mp->graph_norm, NULL, &numNodes),931);
    //if (mp->myrank == 0) printf("\nGraph: norm number of nodes = %zu\n",numNodes);

    print_CUDA_error_if_any(cudaGraphInstantiate(&mp->graphExec_norm, mp->graph_norm, NULL, NULL, 0),932);
    //if (mp->myrank == 0) printf("\nGraph: norm instantiated\n");

    // graph is initialized, ready to be called by graph from now on
    mp->init_graph_norm = 0;
    mp->use_graph_call_norm = 1;
  }

  // launches graph instead of separate kernels
  if (mp->use_graph_call_norm){
    // graph
    print_CUDA_error_if_any(cudaGraphLaunch(mp->graphExec_norm, mp->compute_stream),930);
    //if (mp->myrank == 0) printf("\nGraph: norm launch\n");

    // makes sure kernels have finished copies to host
    cudaStreamSynchronize(mp->compute_stream);
  }
#endif

  realw max;
  realw max_crust_mantle, max_inner_core;

  realw *h_norm_max = mp->h_norm_max;

  // determines fluid max for all blocks
  if (mp->NGLOB_OUTER_CORE > 0) {
    max = h_norm_max[0];
    for (int i = 1; i < size_block_fluid; i++) {
      if (max < h_norm_max[i]) max = h_norm_max[i];
    }
  }else{
    max = 0.0f;
  }
  // return result
  *fluidnorm = max;

  // determines elastic max for all blocks crust_mantle
  max = h_norm_max[size_block_fluid];
  for (int i = size_block_fluid+1; i < size_block_fluid + size_block_cm; i++) {
    // sets maximum
    if (max < h_norm_max[i]) max = h_norm_max[i];
  }
  max_crust_mantle = max;

  // determines max for all blocks inner_core
  if (mp->NGLOB_INNER_CORE > 0) {
    max = h_norm_max[size_block_fluid + size_block_cm];
    for (int i = size_block_fluid + size_block_cm + 1; i < size_block_fluid + size_block_cm + size_block_ic; i++) {
      if (max < h_norm_max[i]) max = h_norm_max[i];
    }
  }else{
    max = 0.0f;
  }
  max_inner_core = max;

  // return result
  max = MAX (max_inner_core, max_crust_mantle);
  *solidnorm = max;

  //debug
#if (DEBUG_FIELDS == 1)
  {
    printf ("rank %d - norm elastic: size fluid %d cm %d ic %d\n",mp->myrank,size_block_fluid,size_block_cm,size_block_ic);
    printf ("rank %d - norm elastic: size total used %d - crust_mantle = %e inner_core = %e\n",
            mp->myrank,(size_block_fluid + size_block_cm + size_block_ic),max_crust_mantle,max_inner_core);
    fflush (stdout);
    synchronize_mpi ();
  }
#endif

  // frees arrays
  //free (h_max_fluid);
  //free (h_max_cm);
  //free (h_max_ic);

  GPU_ERROR_CHECKING ("after check_norm_elastic_acoustic_from_device");
}


/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (check_norm_strain_from_device,
               CHECK_NORM_STRAIN_FROM_DEVICE) (realw *strain_norm,
                                               realw *strain_norm2,
                                               long *Mesh_pointer_f) {

  TRACE ("check_norm_strain_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int num_blocks_x, num_blocks_y;
  int size_nonpadded, size_padded;

  // initializes
  *strain_norm = 0.f;
  *strain_norm2 = 0.f;

  // checks if anything to do
  if (! mp->compute_and_store_strain) return;

  //debug
#if (DEBUG_FIELDS == 1)
  {
    printf ("rank %d - norm strain: blocksize_transfer %d nspec_strain_only = %d\n",
            mp->myrank,BLOCKSIZE_TRANSFER,mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY);
    fflush (stdout);
    synchronize_mpi ();
    realw max_1, max_2, max_3;
    max_1 = get_device_array_maximum_value(mp->d_epsilondev_xx_crust_mantle, NGLL3 * mp->NSPEC_CRUST_MANTLE);
    max_2 = get_device_array_maximum_value(mp->d_epsilondev_yy_crust_mantle, NGLL3 * mp->NSPEC_CRUST_MANTLE);
    max_3 = get_device_array_maximum_value(mp->d_epsilondev_xy_crust_mantle, NGLL3 * mp->NSPEC_CRUST_MANTLE);
    printf ("rank %d - norm strain: max xx: %e yy: %e xy: %e\n", mp->myrank, max_1, max_2, max_3);
    fflush (stdout);
    synchronize_mpi ();
  }
#endif

  // launch simple reduction kernel
  const int blocksize = BLOCKSIZE_TRANSFER;

  // crust_mantle strain arrays
  size_nonpadded = NGLL3 * (mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY);

  size_padded = ((int) ceil (((double) size_nonpadded) / ((double) blocksize))) * blocksize;
  get_blocks_xy (size_padded / blocksize, &num_blocks_x, &num_blocks_y);

  int size_block_strain = num_blocks_x * num_blocks_y;
  //realw *h_max_strain = (realw *) calloc (size_block_strain, sizeof (realw));
  //if (h_max_strain == NULL) { exit_on_error("Error allocating h_max array in check_norm_strain_from_device() routine"); }

  // graph
#ifdef USE_CUDA_GRAPHS
  if (mp->init_graph_norm_strain){
    // debug: synchronizes first
    //gpuSynchronize(); synchronize_mpi();
    // start capturing
    print_CUDA_error_if_any(cudaStreamBeginCapture(mp->compute_stream),950);
  }
#endif

  // determines max for: eps_trace_over_3_crust_mantle
  if (mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY > 1){

#ifdef USE_OPENCL
    if (run_opencl) {
      size_t local_work_size[2];
      size_t global_work_size[2];
      cl_uint idx = 0;

      local_work_size[0] = blocksize;
      local_work_size[1] = 1;
      global_work_size[0] = num_blocks_x * blocksize;
      global_work_size[1] = num_blocks_y;

      // reduction kernel
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_eps_trace_over_3_crust_mantle.ocl));
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (int), (void *) &size_nonpadded));
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_norm_strain_max.ocl));

      clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.get_maximum_scalar_kernel, 2, NULL,
                                       global_work_size, local_work_size, 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      // graph
#ifdef USE_CUDA_GRAPHS
      if (! mp->use_graph_call_norm_strain){
#endif

      dim3 grid = dim3(num_blocks_x,num_blocks_y);
      dim3 threads = dim3(blocksize,1,1);

      // reduction kernel
      get_maximum_scalar_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_eps_trace_over_3_crust_mantle.cuda,
                                                                       size_nonpadded,
                                                                       mp->d_norm_strain_max.cuda);

      // graph
#ifdef USE_CUDA_GRAPHS
      }
#endif
    }
#endif
#ifdef USE_HIP
    if (run_hip) {
      dim3 grid = dim3(num_blocks_x,num_blocks_y);
      dim3 threads = dim3(blocksize,1,1);

      // reduction kernel
      hipLaunchKernelGGL(HIP_KERNEL_NAME(get_maximum_scalar_kernel), grid, threads, 0, mp->compute_stream,
                                                                     mp->d_eps_trace_over_3_crust_mantle.hip,
                                                                     size_nonpadded,
                                                                     mp->d_norm_strain_max.hip);
    }
#endif

  } // NSPEC_CRUST_MANTLE_STRAIN_ONLY

  // crust_mantle arrays
  size_nonpadded = NGLL3 * mp->NSPEC_CRUST_MANTLE;

  size_padded = ((int) ceil (((double) size_nonpadded) / ((double) blocksize))) * blocksize;
  get_blocks_xy (size_padded / blocksize, &num_blocks_x, &num_blocks_y);

  int size_block = num_blocks_x * num_blocks_y;
  //realw *h_max = (realw *) calloc (5 * size_block, sizeof (realw));
  //if (h_max == NULL) { exit_on_error("Error allocating h_max array in check_norm_strain_from_device() routine"); }

  // determines max for: epsilondev_xx_crust_mantle,..
  int loop;
  gpu_realw_mem d_epsilondev_HH_crust_mantle[] =
    {mp->d_epsilondev_xx_crust_mantle,
     mp->d_epsilondev_yy_crust_mantle,
     mp->d_epsilondev_xy_crust_mantle,
     mp->d_epsilondev_xz_crust_mantle,
     mp->d_epsilondev_yz_crust_mantle
    };

  INITIALIZE_OFFSET();
  gpu_realw_mem tmp;
  int offset;

  for (loop = 0; loop < 5; loop++) {
    offset = size_block_strain + loop * size_block;
    INIT_OFFSET(d_norm_strain_max, offset);

    tmp = gpuTakeRef(PASS_OFFSET(d_norm_strain_max, offset));

    // gets maximum by gpu reduction kernel
#ifdef USE_OPENCL
    if (run_opencl) {
      size_t local_work_size[2];
      size_t global_work_size[2];
      cl_uint idx = 0;

      local_work_size[0] = blocksize;
      local_work_size[1] = 1;
      global_work_size[0] = num_blocks_x * blocksize;
      global_work_size[1] = num_blocks_y;

      // determines max
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (cl_mem), (void *) &d_epsilondev_HH_crust_mantle[loop].ocl));
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (int), (void *) &size_nonpadded));
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (cl_mem), (void *) &tmp.ocl));

      clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.get_maximum_scalar_kernel, 2, NULL,
                                       global_work_size, local_work_size, 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      // graph
#ifdef USE_CUDA_GRAPHS
      if (! mp->use_graph_call_norm_strain){
#endif

      dim3 grid = dim3(num_blocks_x,num_blocks_y);
      dim3 threads = dim3(blocksize,1,1);

      // determines max for: epsilondev_xx_crust_mantle
      get_maximum_scalar_kernel<<<grid,threads,0,mp->compute_stream>>>(d_epsilondev_HH_crust_mantle[loop].cuda,
                                                                       size_nonpadded,
                                                                       tmp.cuda);

      // graph
#ifdef USE_CUDA_GRAPHS
      }
#endif
    }
#endif
#ifdef USE_HIP
    if (run_hip) {
      dim3 grid = dim3(num_blocks_x,num_blocks_y);
      dim3 threads = dim3(blocksize,1,1);

      // determines max for: epsilondev_xx_crust_mantle
      hipLaunchKernelGGL(HIP_KERNEL_NAME(get_maximum_scalar_kernel), grid, threads, 0, mp->compute_stream,
                                                                     d_epsilondev_HH_crust_mantle[loop].hip,
                                                                     size_nonpadded,
                                                                     tmp.hip);
    }
#endif

    RELEASE_OFFSET(d_norm_strain_max, offset);
  }

  // graph
#ifdef USE_CUDA_GRAPHS
  if (! mp->use_graph_call_norm_strain){
#endif
  // note: we use event synchronization to make sure that compute stream has finished before copying data to host.
  //       the same for continuing only after copy to host has finished.
  //
  //       capturing graphs seems only to work with this sort of event synchronization...

  // create event for synchronizing with async copy
  gpuRecordEvent(mp);

  // copies array to CPU
  gpuCopy_from_device_realw_asyncEvent (mp, &mp->d_norm_strain_max, mp->h_norm_strain_max, size_block_strain + 5 * size_block);

  // makes sure copy has finished
  gpuWaitEvent(mp);

  // graph
#ifdef USE_CUDA_GRAPHS
  }
#endif


#ifdef USE_CUDA
  if (run_cuda) {
    // graph
#ifdef USE_CUDA_GRAPHS
    // finish creating graph
    if (mp->init_graph_norm_strain){
      // stop capturing
      print_CUDA_error_if_any(cudaStreamEndCapture(mp->compute_stream, &mp->graph_norm_strain),950);

      // get graph info
      size_t numNodes = 0;
      print_CUDA_error_if_any(cudaGraphGetNodes(mp->graph_norm_strain, NULL, &numNodes),951);
      //if (mp->myrank == 0) printf("\nGraph: norm strain number of nodes = %zu\n",numNodes);

      print_CUDA_error_if_any(cudaGraphInstantiate(&mp->graphExec_norm_strain, mp->graph_norm_strain, NULL, NULL, 0),952);
      //if (mp->myrank == 0) printf("\nGraph: norm strain instantiated\n");

      // graph is initialized, ready to be called by graph from now on
      mp->init_graph_norm_strain = 0;
      mp->use_graph_call_norm_strain = 1;
    }

    // launches graph instead of separate kernels
    if (mp->use_graph_call_norm_strain){
      // graph
      print_CUDA_error_if_any(cudaGraphLaunch(mp->graphExec_norm_strain, mp->compute_stream),953);
      //if (mp->myrank == 0) printf("\nGraph: norm strain launch \n");

      // makes sure kernels have finished copies to host
      cudaStreamSynchronize(mp->compute_stream);
    }
#endif
  }
#endif

  realw max;
  realw *h_max = mp->h_norm_strain_max;

  // determines maximum
  // strain trace
  if (mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY > 1){
    max = h_max[0];
    for (int i = 1; i < size_block_strain; i++) {
      if (max < h_max[i]) max = h_max[i];
    }
    // strain trace maximum
    *strain_norm = max;
  }

  // strain
  max = h_max[size_block_strain];
  for (int i = size_block_strain + 1; i < size_block_strain + 5 * size_block; i++) {
    if (max < h_max[i]) max = h_max[i];
  }
  //max_eps = MAX (max_eps, max);
  // strain maximum
  *strain_norm2 = max;

  //debug
#if (DEBUG_FIELDS == 1)
  {
    printf ("rank %d - norm strain: size used %d - strain trace = %e strain = %e \n",
            mp->myrank,(size_block_strain + 5 * size_block),*strain_norm,*strain_norm);
    fflush (stdout);
    synchronize_mpi ();
  }
#endif

  // frees arrays
  //free (h_max_strain);
  //free (h_max);

  GPU_ERROR_CHECKING ("after check_norm_strain_from_device");
}
