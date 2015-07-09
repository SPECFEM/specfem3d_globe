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
    int i;
    for (i = 1; i < size; i++) {
      if (fabs(h_array[i]) > max)
        max = fabs(h_array[i]);
    }

    // frees temporary array
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

  // get Mesh from Fortran integer wrapper
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  realw max;
  realw *h_max;

  // initializes
  *norm = 0.f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in check_norm_acoustic_from_device() routine");
  }

  // launch simple reduction kernel
  int blocksize = BLOCKSIZE_TRANSFER;

  // outer core
  int size = mp->NGLOB_OUTER_CORE;

  int size_padded = ((int) ceil (((double) size)/ ((double) blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

  h_max = (realw *) calloc (num_blocks_x * num_blocks_y, sizeof (realw));
  if (h_max == NULL) { exit_on_error("Error allocating h_max array in check_norm_acoustic_from_device() routine"); }

  // sets gpu arrays
  gpu_realw_mem displ;
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
    clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (int), (void *) &size));
    clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_norm_max.ocl));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.get_maximum_scalar_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    get_maximum_scalar_kernel<<<grid,threads,0,mp->compute_stream>>>(displ.cuda,size,mp->d_norm_max.cuda);
  }
#endif

  // copies to CPU
  gpuCopy_from_device_realw (&mp->d_norm_max, h_max, num_blocks_x*num_blocks_y);

  //debug
  if (DEBUG_FIELDS) {
    realw max_d, max_v, max_a;
    max_d = get_device_array_maximum_value(mp->d_displ_outer_core, mp->NGLOB_OUTER_CORE);
    max_v = get_device_array_maximum_value(mp->d_veloc_outer_core, mp->NGLOB_OUTER_CORE);
    max_a = get_device_array_maximum_value(mp->d_accel_outer_core, mp->NGLOB_OUTER_CORE);
    printf ("rank %d - check norm max outer_core displ: %e veloc: %e accel: %e, %i\n", mp->myrank, max_d, max_v, max_a, *FORWARD_OR_ADJOINT);
    fflush (stdout);
    synchronize_mpi ();
  }

  // determines max for all blocks
  max = h_max[0];
  int i;
  for (i = 1; i < num_blocks_x * num_blocks_y; i++) {
    if (max < h_max[i]) max = h_max[i];
  }

  // frees arrays
  free (h_max);

  // return result
  *norm = max;

  GPU_ERROR_CHECKING ("after check_norm_acoustic_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (check_norm_elastic_from_device,
               CHECK_NORM_ELASTIC_FROM_DEVICE) (realw *norm,
                                                long *Mesh_pointer_f,
                                                int *FORWARD_OR_ADJOINT) {

  TRACE ("check_norm_elastic_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  realw max_d, max_v, max_a;
  realw max, max_crust_mantle, max_inner_core;

  int size, size_padded;
  realw *h_max;

  // initializes
  *norm = 0.f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in check_norm_elastic_from_device() routine");
  }

  // launch simple reduction kernel
  int blocksize = BLOCKSIZE_TRANSFER;

  // crust_mantle
  size = mp->NGLOB_CRUST_MANTLE;

  size_padded = ((int) ceil (((double) size) / ((double) blocksize))) * blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (size_padded / blocksize, &num_blocks_x, &num_blocks_y);

  h_max = (realw *) calloc (num_blocks_x * num_blocks_y, sizeof (realw));
  if (h_max == NULL) { exit_on_error("Error allocating h_max array in check_norm_elastic_from_device() routine"); }

  // sets gpu arrays
  gpu_realw_mem displ;
  if (*FORWARD_OR_ADJOINT == 1) {
    displ = mp->d_displ_crust_mantle;
  } else {
    displ = mp->d_b_displ_crust_mantle;
  }

#ifdef USE_OPENCL
  size_t global_work_size[2];
  size_t local_work_size[2];
  cl_uint idx = 0;

  if (run_opencl) {
    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (cl_mem), (void *) &displ.ocl));
    clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (int), (void *) &size));
    clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_norm_max.ocl));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.get_maximum_vector_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  dim3 grid,threads;
  if (run_cuda) {
    grid = dim3(num_blocks_x,num_blocks_y);
    threads = dim3(blocksize,1,1);

    get_maximum_vector_kernel<<<grid,threads,0,mp->compute_stream>>>(displ.cuda,size,mp->d_norm_max.cuda);
  }
#endif

  // copies to CPU
  gpuCopy_from_device_realw (&mp->d_norm_max, h_max, num_blocks_x*num_blocks_y);

  //debug
  if (DEBUG_FIELDS) {
    max_d = get_device_array_maximum_value(mp->d_displ_crust_mantle, NDIM * mp->NGLOB_CRUST_MANTLE);
    max_v = get_device_array_maximum_value(mp->d_veloc_crust_mantle, NDIM * mp->NGLOB_CRUST_MANTLE);
    max_a = get_device_array_maximum_value(mp->d_accel_crust_mantle, NDIM * mp->NGLOB_CRUST_MANTLE);
    printf ("rank %d - check norm max crust_mantle displ: %e veloc: %e accel: %e, %i\n", mp->myrank, max_d, max_v, max_a, *FORWARD_OR_ADJOINT);
    fflush (stdout);
    synchronize_mpi ();
  }

  // determines max for all blocks
  max = h_max[0];
  int i;
  for (i = 1; i < num_blocks_x * num_blocks_y; i++) {
    // sets maximum
    if (max < h_max[i]) max = h_max[i];
  }
  max_crust_mantle = max;

  // frees arrays
  free (h_max);

  // inner_core
  size = mp->NGLOB_INNER_CORE;

  size_padded = ((int) ceil (((double) size) / ((double) blocksize))) * blocksize;

  get_blocks_xy (size_padded / blocksize, &num_blocks_x, &num_blocks_y);

  h_max = (realw *) calloc (num_blocks_x * num_blocks_y, sizeof (realw));
  if (h_max == NULL) { exit_on_error("Error allocating h_max array for inner core in check_norm_elastic_from_device() routine"); }

  // sets gpu arrays
  if (*FORWARD_OR_ADJOINT == 1) {
    displ = mp->d_displ_inner_core;
  } else {
    displ = mp->d_b_displ_inner_core;
  }

#ifdef USE_OPENCL
  if (run_opencl) {
    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    idx = 0;

    clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (cl_mem), (void *) &displ.ocl));
    clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (int), (void *) &size));
    clCheck (clSetKernelArg (mocl.kernels.get_maximum_vector_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_norm_max.ocl));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.get_maximum_vector_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    grid = dim3(num_blocks_x,num_blocks_y);
    threads = dim3(blocksize,1,1);

    get_maximum_vector_kernel<<<grid,threads,0,mp->compute_stream>>>(displ.cuda,size,mp->d_norm_max.cuda);
  }
#endif

  // copies to CPU
  gpuCopy_from_device_realw (&mp->d_norm_max, h_max, num_blocks_x*num_blocks_y);

  //debug
  if (DEBUG_FIELDS) {
    max_d = get_device_array_maximum_value(mp->d_displ_inner_core, NDIM * mp->NGLOB_INNER_CORE);
    max_v = get_device_array_maximum_value(mp->d_veloc_inner_core, NDIM * mp->NGLOB_INNER_CORE);
    max_a = get_device_array_maximum_value(mp->d_accel_inner_core, NDIM * mp->NGLOB_INNER_CORE);
    printf ("rank %d - check norm max inner_core displ: %e veloc: %e accel: %e, %i\n", mp->myrank, max_d, max_v, max_a, *FORWARD_OR_ADJOINT);
    fflush (stdout);
    synchronize_mpi ();
  }

  // determines max for all blocks
  max = h_max[0];
  for (i = 1; i < num_blocks_x * num_blocks_y; i++) {
    if (max < h_max[i]) max = h_max[i];
  }
  max_inner_core = max;

  // frees arrays
  free (h_max);

  //debug
  if (DEBUG_FIELDS) {
    printf ("rank %d - norm elastic: crust_mantle = %e inner_core = %e \n",mp->myrank,max_crust_mantle,max_inner_core);
    fflush (stdout);
    synchronize_mpi ();
  }

  // return result
  max = MAX (max_inner_core, max_crust_mantle);
  *norm = max;

  GPU_ERROR_CHECKING ("after check_norm_elastic_from_device");
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

  realw max, max_eps;
  int num_blocks_x, num_blocks_y;
  int size, size_padded;

  // launch simple reduction kernel
  realw *h_max;
  int blocksize = BLOCKSIZE_TRANSFER;

  // crust_mantle strain arrays
  size = NGLL3 * (mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY);

  size_padded = ((int) ceil (((double) size) / ((double) blocksize))) * blocksize;

  get_blocks_xy (size_padded / blocksize, &num_blocks_x, &num_blocks_y);

  h_max = (realw *) calloc (num_blocks_x * num_blocks_y, sizeof (realw));
  if (h_max == NULL) { exit_on_error("Error allocating h_max array in check_norm_strain_from_device() routine"); }

  // determines max for: eps_trace_over_3_crust_mantle
#ifdef USE_OPENCL
  size_t local_work_size[2];
  size_t global_work_size[2];
  cl_uint idx = 0;
  if (run_opencl) {
    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    // reduction kernel
    idx = 0;
    clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_eps_trace_over_3_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (int), (void *) &size));
    clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_norm_max.ocl));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.get_maximum_scalar_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  dim3 grid,threads;
  if (run_cuda) {
    grid = dim3(num_blocks_x,num_blocks_y);
    threads = dim3(blocksize,1,1);
    // reduction kernel
    get_maximum_scalar_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_eps_trace_over_3_crust_mantle.cuda,size,mp->d_norm_max.cuda);
  }
#endif

  // copies to CPU
  gpuCopy_from_device_realw (&mp->d_norm_max, h_max, num_blocks_x*num_blocks_y);

  // determines maximum
  max = h_max[0];
  int i;
  for (i = 1; i < num_blocks_x * num_blocks_y; i++) {
    if (max < h_max[i])
      max = h_max[i];
  }
  // strain trace maximum
  *strain_norm = max;

  // frees arrays
  free (h_max);

  // initializes
  // crust_mantle arrays
  size = NGLL3 * mp->NSPEC_CRUST_MANTLE;

  size_padded = ((int) ceil (((double) size) / ((double) blocksize))) * blocksize;

  get_blocks_xy (size_padded / blocksize, &num_blocks_x, &num_blocks_y);

  h_max = (realw *) calloc (num_blocks_x * num_blocks_y, sizeof (realw));
  if (h_max == NULL) { exit_on_error("Error allocating h_max array in check_norm_strain_from_device() routine"); }

  max_eps = 0.0f;

  // determines max for: epsilondev_xx_crust_mantle,..
  int loop;
  gpu_realw_mem d_epsilondev_HH_crust_mantle[] =
    {mp->d_epsilondev_xx_crust_mantle,
     mp->d_epsilondev_yy_crust_mantle,
     mp->d_epsilondev_xy_crust_mantle,
     mp->d_epsilondev_xz_crust_mantle,
     mp->d_epsilondev_yz_crust_mantle
    };
  for (loop = 0; loop < 5; loop++) {
    // gets maximum by gpu reduction kernel
#ifdef USE_OPENCL
    if (run_opencl) {
      local_work_size[0] = blocksize;
      local_work_size[1] = 1;
      global_work_size[0] = num_blocks_x * blocksize;
      global_work_size[1] = num_blocks_y;

      // determines max
      idx = 0;
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (cl_mem), (void *) &d_epsilondev_HH_crust_mantle[loop].ocl));
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (int), (void *) &size));
      clCheck (clSetKernelArg (mocl.kernels.get_maximum_scalar_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_norm_max.ocl));

      clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.get_maximum_scalar_kernel, 2, NULL,
                                       global_work_size, local_work_size, 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      grid = dim3(num_blocks_x,num_blocks_y);
      threads = dim3(blocksize,1,1);

      // determines max for: epsilondev_xx_crust_mantle
      get_maximum_scalar_kernel<<<grid,threads,0,mp->compute_stream>>>(d_epsilondev_HH_crust_mantle[loop].cuda,size,mp->d_norm_max.cuda);
    }
#endif
    // copies array to CPU
    gpuCopy_from_device_realw (&mp->d_norm_max, h_max, num_blocks_x * num_blocks_y);

    // determines maximum
    max = h_max[0];
    for (i = 1; i < num_blocks_x * num_blocks_y; i++) {
      if (max < h_max[i]) max = h_max[i];
    }
    max_eps = MAX (max_eps, max);
  }

  // strain maximum
  *strain_norm2 = max_eps;

  // frees arrays
  free (h_max);

  GPU_ERROR_CHECKING ("after check_norm_strain_from_device");
}
