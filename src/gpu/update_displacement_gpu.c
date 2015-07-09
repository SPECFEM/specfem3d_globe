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

extern EXTERN_LANG
void FC_FUNC_ (update_displacement_ic_gpu,
               UPDATE_DISPLACMENT_IC_GPU) (long *Mesh_pointer_f,
                                           realw *deltat_f,
                                           realw *deltatsqover2_f,
                                           realw *deltatover2_f,
                                           int *FORWARD_OR_ADJOINT) {

  TRACE ("update_displacement_ic_gpu");

  //get Mesh from Fortran integer wrapper
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in update_displacement_ic_gpu() routine");
  }

  // inner core
  int size = NDIM * mp->NGLOB_INNER_CORE;

#if (DEBUG_BACKWARD_SIMULATIONS == 1 && DEBUG == 1 ) || DEBUG_FIELDS == 1
  //debug
  realw max_d, max_v, max_a;
  max_d = get_device_array_maximum_value(mp->d_b_displ_inner_core, size);
  max_v = get_device_array_maximum_value(mp->d_b_veloc_inner_core, size);
  max_a = get_device_array_maximum_value(mp->d_b_accel_inner_core, size);
  printf ("rank %d - forward/adjoint: %i, max inner_core displ: %f veloc: %f accel: %f\n", mp->myrank, *FORWARD_OR_ADJOINT, max_d, max_v, max_a);
  fflush (stdout);
  synchronize_mpi ();
#endif

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ( (int)ceil ( ( (double)size)/ ( (double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

  realw deltat = *deltat_f;
  realw deltatsqover2 = *deltatsqover2_f;
  realw deltatover2 = *deltatover2_f;

  // sets gpu arrays
  gpu_realw_mem displ,veloc,accel;
  if (*FORWARD_OR_ADJOINT == 1) {
    displ = mp->d_displ_inner_core;
    veloc = mp->d_veloc_inner_core;
    accel = mp->d_accel_inner_core;
  } else {
    //debug
    DEBUG_BACKWARD_UPDATE ();
    // for backward/reconstructed fields
    displ = mp->d_b_displ_inner_core;
    veloc = mp->d_b_veloc_inner_core;
    accel = mp->d_b_accel_inner_core;
  }

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    //launch kernel
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (cl_mem), (void *) &displ.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (cl_mem), (void *) &veloc.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (int), (void *) &size));
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (realw), (void *) &deltat));
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (realw), (void *) &deltatsqover2));
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (realw), (void *) &deltatover2));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.update_disp_veloc_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    //launch kernel
    update_disp_veloc_kernel<<<grid,threads,0,mp->compute_stream>>>(displ.cuda,
                                                                    veloc.cuda,
                                                                    accel.cuda,
                                                                    size,deltat,deltatsqover2,deltatover2);
  }
#endif

  GPU_ERROR_CHECKING ("update_displacement_ic_gpu");
}

/*----------------------------------------------------------------------------------------------- */
//KERNEL 1
//crust/mantle
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (update_displacement_cm_gpu,
               UPDATE_DISPLACMENT_CM_GPU) (long *Mesh_pointer_f,
                                           realw *deltat_f,
                                           realw *deltatsqover2_f,
                                           realw *deltatover2_f,
                                           int *FORWARD_OR_ADJOINT) {

  TRACE ("update_displacement_cm_gpu");

  //get Mesh from Fortran integer wrapper
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in update_displacement_cm_gpu() routine");
  }

  // crust/mantle
  int size = NDIM * mp->NGLOB_CRUST_MANTLE;

  //debug
#if (DEBUG_BACKWARD_SIMULATIONS == 1 && DEBUG == 1) || DEBUG_FIELDS == 1
  realw max_d, max_v, max_a;
  max_d = get_device_array_maximum_value(mp->d_b_displ_crust_mantle, size);
  max_v = get_device_array_maximum_value(mp->d_b_veloc_crust_mantle, size);
  max_a = get_device_array_maximum_value(mp->d_b_accel_crust_mantle, size);
  printf ("rank %d - forward/adjoint: %i, max crust_mantle displ: %f veloc: %f accel: %f\n", mp->myrank, *FORWARD_OR_ADJOINT, max_d, max_v, max_a);
  fflush (stdout);
  synchronize_mpi ();
#endif

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ( (int)ceil ( ( (double)size)/ ( (double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

  realw deltat = *deltat_f;
  realw deltatsqover2 = *deltatsqover2_f;
  realw deltatover2 = *deltatover2_f;

  // sets gpu arrays
  gpu_realw_mem displ,veloc,accel;
  if (*FORWARD_OR_ADJOINT == 1) {
    displ = mp->d_displ_crust_mantle;
    veloc = mp->d_veloc_crust_mantle;
    accel = mp->d_accel_crust_mantle;
  } else {
    //debug
    DEBUG_BACKWARD_UPDATE ();
    // for backward/reconstructed fields
    displ = mp->d_b_displ_crust_mantle;
    veloc = mp->d_b_veloc_crust_mantle;
    accel = mp->d_b_accel_crust_mantle;
  }

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    //launch kernel
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (cl_mem), (void *) &displ.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (cl_mem), (void *) &veloc.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (int), (void *) &size));
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (realw), (void *) &deltat));
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (realw), (void *) &deltatsqover2));
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (realw), (void *) &deltatover2));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.update_disp_veloc_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    //launch kernel
    update_disp_veloc_kernel<<<grid,threads,0,mp->compute_stream>>>(displ.cuda,
                                                                    veloc.cuda,
                                                                    accel.cuda,
                                                                    size,deltat,deltatsqover2,deltatover2);
  }
#endif

  GPU_ERROR_CHECKING ("update_displacement_cm_gpu");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (update_displacement_oc_gpu,
               UPDATE_DISPLACEMENT_OC_gpu) (long *Mesh_pointer_f,
                                            realw *deltat_f,
                                            realw *deltatsqover2_f,
                                            realw *deltatover2_f,
                                            int *FORWARD_OR_ADJOINT) {

  TRACE ("update_displacement_oc_gpu");

  //get Mesh from Fortran integer wrapper
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in update_displacement_oc_gpu() routine");
  }

  // outer core
  int size = mp->NGLOB_OUTER_CORE;

  //debug
#if (DEBUG_BACKWARD_SIMULATIONS == 1 && DEBUG == 1) || DEBUG_FIELDS == 1
  realw max_d, max_v, max_a;
  max_d = get_device_array_maximum_value(mp->d_b_displ_outer_core, size);
  max_v = get_device_array_maximum_value(mp->d_b_veloc_outer_core, size);
  max_a = get_device_array_maximum_value(mp->d_b_accel_outer_core, size);
  printf ("rank %d - forward/adjoint: %i, max outer_core displ: %f veloc: %f accel: %f\n", mp->myrank, *FORWARD_OR_ADJOINT, max_d, max_v, max_a);
  fflush (stdout);
  synchronize_mpi ();
#endif

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ( (int)ceil ( ( (double)size)/ ( (double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

  realw deltat = *deltat_f;
  realw deltatsqover2 = *deltatsqover2_f;
  realw deltatover2 = *deltatover2_f;

  // sets gpu arrays
  gpu_realw_mem displ,veloc,accel;
  if (*FORWARD_OR_ADJOINT == 1) {
    displ = mp->d_displ_outer_core;
    veloc = mp->d_veloc_outer_core;
    accel = mp->d_accel_outer_core;
  } else {
    //debug
    DEBUG_BACKWARD_UPDATE ();
    // for backward/reconstructed fields
    displ = mp->d_b_displ_outer_core;
    veloc = mp->d_b_veloc_outer_core;
    accel = mp->d_b_accel_outer_core;
  }

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    //launch kernel
    clCheck (clSetKernelArg (mocl.kernels.update_potential_kernel, idx++, sizeof (cl_mem), (void *) &displ.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_potential_kernel, idx++, sizeof (cl_mem), (void *) &veloc.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_potential_kernel, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_potential_kernel, idx++, sizeof (int), (void *) &size));
    clCheck (clSetKernelArg (mocl.kernels.update_potential_kernel, idx++, sizeof (realw), (void *) &deltat));
    clCheck (clSetKernelArg (mocl.kernels.update_potential_kernel, idx++, sizeof (realw), (void *) &deltatsqover2));
    clCheck (clSetKernelArg (mocl.kernels.update_potential_kernel, idx++, sizeof (realw), (void *) &deltatover2));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.update_potential_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    //launch kernel
    update_potential_kernel<<<grid,threads,0,mp->compute_stream>>>(displ.cuda,
                                                                   veloc.cuda,
                                                                   accel.cuda,
                                                                   size,deltat,deltatsqover2,deltatover2);
  }
#endif

  GPU_ERROR_CHECKING ("update_displacement_oc_gpu");
}

/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (multiply_accel_elastic_gpu,
               MULTIPLY_ACCEL_ELASTIC_GPU) (long *Mesh_pointer_f,
                                            int *FORWARD_OR_ADJOINT) {
  TRACE ("multiply_accel_elastic_gpu");

  int size_padded, num_blocks_x, num_blocks_y;

  //get Mesh from Fortran integer wrapper
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in multiply_accel_elastic_gpu() routine");
  }

  // update kernel
  int blocksize = BLOCKSIZE_KERNEL3;

  //multiplies accel with inverse of mass matrix

  //crust/mantle region
  size_padded = ( (int)ceil ( ( (double)mp->NGLOB_CRUST_MANTLE)/ ( (double)blocksize)))*blocksize;

  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

  // sets gpu arrays
  gpu_realw_mem veloc,accel;
  gpu_realw_mem rmassx,rmassy,rmassz;
  realw two_omega_earth;
  if (*FORWARD_OR_ADJOINT == 1) {
    veloc = mp->d_veloc_crust_mantle;
    accel = mp->d_accel_crust_mantle;
    rmassx = mp->d_rmassx_crust_mantle;
    rmassy = mp->d_rmassy_crust_mantle;
    rmassz = mp->d_rmassz_crust_mantle;
    two_omega_earth = mp->two_omega_earth;
  } else {
    //debug
    DEBUG_BACKWARD_UPDATE ();
    // for backward/reconstructed fields
    veloc = mp->d_b_veloc_crust_mantle;
    accel = mp->d_b_accel_crust_mantle;
    rmassx = mp->d_b_rmassx_crust_mantle;
    rmassy = mp->d_b_rmassy_crust_mantle;
    rmassz = mp->d_b_rmassz_crust_mantle;
    two_omega_earth = mp->b_two_omega_earth;
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

    clCheck (clSetKernelArg (mocl.kernels.update_accel_elastic_kernel, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_accel_elastic_kernel, idx++, sizeof (cl_mem), (void *) &veloc.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_accel_elastic_kernel, idx++, sizeof (int), (void *) &mp->NGLOB_CRUST_MANTLE));
    clCheck (clSetKernelArg (mocl.kernels.update_accel_elastic_kernel, idx++, sizeof (realw), (void *) &two_omega_earth));
    clCheck (clSetKernelArg (mocl.kernels.update_accel_elastic_kernel, idx++, sizeof (cl_mem), (void *) &rmassx.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_accel_elastic_kernel, idx++, sizeof (cl_mem), (void *) &rmassy.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_accel_elastic_kernel, idx++, sizeof (cl_mem), (void *) &rmassz.ocl));
    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.update_accel_elastic_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  dim3 grid,threads;

  if (run_cuda) {
    grid = dim3(num_blocks_x,num_blocks_y);
    threads = dim3(blocksize,1,1);

    // launches kernel
    update_accel_elastic_kernel<<< grid, threads,0,mp->compute_stream>>>(accel.cuda,
                                                                         veloc.cuda,
                                                                         mp->NGLOB_CRUST_MANTLE,
                                                                         two_omega_earth,
                                                                         rmassx.cuda,
                                                                         rmassy.cuda,
                                                                         rmassz.cuda);
  }
#endif

  //inner core region
  size_padded = ( (int)ceil ( ( (double)mp->NGLOB_INNER_CORE)/ ( (double)blocksize)))*blocksize;

  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

  // sets gpu arrays
  if (*FORWARD_OR_ADJOINT == 1) {
    veloc = mp->d_veloc_inner_core;
    accel = mp->d_accel_inner_core;
    rmassx = mp->d_rmassx_inner_core;
    rmassy = mp->d_rmassy_inner_core;
    rmassz = mp->d_rmassz_inner_core;
  } else {
    //debug
    DEBUG_BACKWARD_UPDATE ();
    // for backward/reconstructed fields
    veloc = mp->d_b_veloc_inner_core;
    accel = mp->d_b_accel_inner_core;
    rmassx = mp->d_b_rmassx_inner_core;
    rmassy = mp->d_b_rmassy_inner_core;
    rmassz = mp->d_b_rmassz_inner_core;
  }

#ifdef USE_OPENCL
  if (run_opencl) {
    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;
    idx = 0;

    clCheck (clSetKernelArg (mocl.kernels.update_accel_elastic_kernel, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_accel_elastic_kernel, idx++, sizeof (cl_mem), (void *) &veloc.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_accel_elastic_kernel, idx++, sizeof (int), (void *) &mp->NGLOB_INNER_CORE));
    clCheck (clSetKernelArg (mocl.kernels.update_accel_elastic_kernel, idx++, sizeof (realw), (void *) &two_omega_earth));
    clCheck (clSetKernelArg (mocl.kernels.update_accel_elastic_kernel, idx++, sizeof (cl_mem), (void *) &rmassx.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_accel_elastic_kernel, idx++, sizeof (cl_mem), (void *) &rmassy.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_accel_elastic_kernel, idx++, sizeof (cl_mem), (void *) &rmassz.ocl));
    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.update_accel_elastic_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    grid = dim3(num_blocks_x,num_blocks_y);
    threads = dim3(blocksize,1,1);

    // launches kernel
    update_accel_elastic_kernel<<< grid, threads,0,mp->compute_stream>>>(accel.cuda,
                                                                         veloc.cuda,
                                                                         mp->NGLOB_INNER_CORE,
                                                                         two_omega_earth,
                                                                         rmassx.cuda,
                                                                         rmassy.cuda,
                                                                         rmassz.cuda);
  }
#endif

  GPU_ERROR_CHECKING ("after multiply_accel_elastic_gpu");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (update_veloc_elastic_gpu,
               UPDATE_VELOC_ELASTIC_GPU) (long *Mesh_pointer_f,
                                          realw *deltatover2_f,
                                          int *FORWARD_OR_ADJOINT) {

  TRACE ("update_veloc_elastic_gpu");

  int size_padded, num_blocks_x, num_blocks_y;

  //get Mesh from Fortran integer wrapper
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  realw deltatover2 = *deltatover2_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in update_veloc_elastic_gpu() routine");
  }

  // update kernel
  int blocksize = BLOCKSIZE_KERNEL3;

  //updates velocity

  //crust/mantle region
  size_padded = ( (int)ceil ( ( (double)mp->NGLOB_CRUST_MANTLE)/ ( (double)blocksize)))*blocksize;

  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

  // sets gpu arrays
  gpu_realw_mem veloc,accel;
  if (*FORWARD_OR_ADJOINT == 1) {
    veloc = mp->d_veloc_crust_mantle;
    accel = mp->d_accel_crust_mantle;
  } else {
    //debug
    DEBUG_BACKWARD_UPDATE ();
    // for backward/reconstructed fields
    veloc = mp->d_b_veloc_crust_mantle;
    accel = mp->d_b_accel_crust_mantle;
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

    clCheck (clSetKernelArg (mocl.kernels.update_veloc_elastic_kernel, idx++, sizeof (cl_mem), (void *) &veloc.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_veloc_elastic_kernel, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_veloc_elastic_kernel, idx++, sizeof (int), (void *) &mp->NGLOB_CRUST_MANTLE));
    clCheck (clSetKernelArg (mocl.kernels.update_veloc_elastic_kernel, idx++, sizeof (realw), (void *) &deltatover2));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.update_veloc_elastic_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  dim3 grid,threads;

  if (run_cuda) {
    grid = dim3(num_blocks_x,num_blocks_y);
    threads = dim3(blocksize,1,1);

    // launches kernel
    update_veloc_elastic_kernel<<< grid, threads,0,mp->compute_stream>>>(veloc.cuda,
                                                                         accel.cuda,
                                                                         mp->NGLOB_CRUST_MANTLE,
                                                                         deltatover2);
  }
#endif

  //inner core region
  size_padded = ((int) ceil (((double) mp->NGLOB_INNER_CORE) / ((double) blocksize))) * blocksize;

  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

  // sets gpu arrays
  if (*FORWARD_OR_ADJOINT == 1) {
    veloc = mp->d_veloc_inner_core;
    accel = mp->d_accel_inner_core;
  } else {
    //debug
    DEBUG_BACKWARD_UPDATE ();
    // for backward/reconstructed fields
    veloc = mp->d_b_veloc_inner_core;
    accel = mp->d_b_accel_inner_core;
  }

#ifdef USE_OPENCL
  if (run_opencl) {
    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;
    idx = 0;

    clCheck (clSetKernelArg (mocl.kernels.update_veloc_elastic_kernel, idx++, sizeof (cl_mem), (void *) &veloc.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_veloc_elastic_kernel, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_veloc_elastic_kernel, idx++, sizeof (int), (void *) &mp->NGLOB_INNER_CORE));
    clCheck (clSetKernelArg (mocl.kernels.update_veloc_elastic_kernel, idx++, sizeof (realw), (void *) &deltatover2));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.update_veloc_elastic_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    grid = dim3(num_blocks_x,num_blocks_y);
    threads = dim3(blocksize,1,1);

    // launches kernel
    update_veloc_elastic_kernel<<< grid, threads,0,mp->compute_stream>>>(veloc.cuda,
                                                                         accel.cuda,
                                                                         mp->NGLOB_INNER_CORE,
                                                                         deltatover2);
  }
#endif

  GPU_ERROR_CHECKING ("after update_veloc_3_b");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (multiply_accel_acoustic_gpu,
               MULTIPLY_ACCEL_ACOUSTIC_GPU) (long *Mesh_pointer_f,
                                             int *FORWARD_OR_ADJOINT) {
  TRACE ("multiply_accel_acoustic_gpu");

  //get Mesh from Fortran integer wrapper
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in multiply_accel_acoustic_gpu() routine");
  }

  // update kernel
  int blocksize = BLOCKSIZE_KERNEL3;

  // outer core
  int size_padded = ( (int)ceil ( ( (double)mp->NGLOB_OUTER_CORE)/ ( (double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

  // sets gpu arrays
  gpu_realw_mem accel,rmass;
  if (*FORWARD_OR_ADJOINT == 1) {
    accel = mp->d_accel_outer_core;
    rmass = mp->d_rmass_outer_core;
  } else {
    //debug
    DEBUG_BACKWARD_UPDATE ();
    // for backward/reconstructed fields
    accel = mp->d_b_accel_outer_core;
    rmass = mp->d_b_rmass_outer_core;
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

    //multiplies accel with inverse of mass matrix
    clCheck (clSetKernelArg (mocl.kernels.update_accel_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_accel_acoustic_kernel, idx++, sizeof (int), (void *) &mp->NGLOB_OUTER_CORE));
    clCheck (clSetKernelArg (mocl.kernels.update_accel_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &rmass.ocl));
    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.update_accel_acoustic_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // multiplies accel with inverse of mass matrix
    update_accel_acoustic_kernel<<< grid, threads,0,mp->compute_stream>>>(accel.cuda,
                                                                          mp->NGLOB_OUTER_CORE,
                                                                          rmass.cuda);
  }
#endif

  GPU_ERROR_CHECKING ("after multiply_accel_acoustic_gpu");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (update_veloc_acoustic_gpu,
               UPDATE_VELOC_ACOUSTIC_GPU) (long *Mesh_pointer_f,
                                           realw *deltatover2_f,
                                           int *FORWARD_OR_ADJOINT) {

  TRACE ("update_veloc_acoustic_gpu");

  //get Mesh from Fortran integer wrapper
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  realw deltatover2 = *deltatover2_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in update_veloc_acoustic_gpu() routine");
  }

  // update kernel
  int blocksize = BLOCKSIZE_KERNEL3;

  // outer core
  int size_padded = ((int)ceil (((double)mp->NGLOB_OUTER_CORE)/ ( (double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

  // sets gpu arrays
  gpu_realw_mem veloc,accel;
  if (*FORWARD_OR_ADJOINT == 1) {
    veloc = mp->d_veloc_outer_core;
    accel = mp->d_accel_outer_core;
  } else {
    //debug
    DEBUG_BACKWARD_UPDATE ();
    // for backward/reconstructed fields
    veloc = mp->d_b_veloc_outer_core;
    accel = mp->d_b_accel_outer_core;
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

    //updates velocity
    clCheck (clSetKernelArg (mocl.kernels.update_veloc_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &veloc.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_veloc_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_veloc_acoustic_kernel, idx++, sizeof (int), (void *) &mp->NGLOB_OUTER_CORE));
    clCheck (clSetKernelArg (mocl.kernels.update_veloc_acoustic_kernel, idx++, sizeof (realw), (void *) &deltatover2));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.update_veloc_acoustic_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // updates velocity
    update_veloc_acoustic_kernel<<< grid, threads,0,mp->compute_stream>>>(veloc.cuda,
                                                                          accel.cuda,
                                                                          mp->NGLOB_OUTER_CORE,
                                                                          deltatover2);
  }
#endif

  GPU_ERROR_CHECKING ("after update_veloc_acoustic_gpu");
}
