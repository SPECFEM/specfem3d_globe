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

extern EXTERN_LANG
void FC_FUNC_ (compute_stacey_acoustic_gpu,
               COMPUTE_STACEY_ACOUSTIC_GPU) (long *Mesh_pointer_f,
                                             realw *absorb_potential) {
  TRACE ("compute_stacey_acoustic_gpu");
  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // absorbing boundary
  int num_abs_boundary_faces = mp->num_abs_boundary_faces_outer_core;

  // checks if anything to do
  if (num_abs_boundary_faces == 0) return;

  // setup arrays
  gpu_int_mem d_abs_boundary_ispec;
  gpu_int_mem d_abs_boundary_npoin;
  gpu_int_mem d_abs_boundary_ijk;
  gpu_realw_mem d_abs_boundary_jacobian2Dw;
  gpu_realw_mem d_absorb_potential;

  d_abs_boundary_ispec = mp->d_abs_boundary_ispec_outer_core;
  d_abs_boundary_npoin = mp->d_abs_boundary_npoin_outer_core;
  d_abs_boundary_ijk = mp->d_abs_boundary_ijk_outer_core;
  d_abs_boundary_jacobian2Dw = mp->d_abs_boundary_jacobian2Dw_outer_core;

  d_absorb_potential = mp->d_absorb_buffer_outer_core;

  // way 1: Elapsed time: 4.385948e-03
  // > NGLLSQUARE==NGLL2==25, but we handle this inside kernel
  //  int blocksize = 32;

  // way 2: Elapsed time: 4.379034e-03
  // > NGLLSQUARE==NGLL2==25, no further check inside kernel
  int blocksize = NGLL2;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (num_abs_boundary_faces, &num_blocks_x, &num_blocks_y);

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    //daniel debug
    //clCheck (clFinish (mocl.command_queue));
    //printf ("rank %d: stacey a %i, %i save %i num blocks x/y= %i %i nglob %i nspec2D %i nspec %i\n",
    //          mp->myrank,interface_type,num_abs_boundary_faces,mp->save_forward,num_blocks_x,num_blocks_y,
    //          mp->NGLOB_OUTER_CORE,mp->nspec2D_ymin_outer_core,mp->NSPEC_OUTER_CORE);
    //fflush (stdout);
    //synchronize_mpi ();

    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_veloc_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_accel_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (int), (void *) &num_abs_boundary_faces));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &d_abs_boundary_ispec.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &d_abs_boundary_npoin.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &d_abs_boundary_ijk.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &d_abs_boundary_jacobian2Dw.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_vp_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (int), (void *) &mp->save_stacey));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &d_absorb_potential.ocl));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_stacey_acoustic_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));

    //daniel debug
    //clCheck (clFinish (mocl.command_queue));
    //printf ("rank %d: stacey b %i, %i \n", mp->myrank,interface_type,num_abs_boundary_faces);
    //fflush (stdout);
    //synchronize_mpi ();
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    compute_stacey_acoustic_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_veloc_outer_core.cuda,
                                                                          mp->d_accel_outer_core.cuda,
                                                                          num_abs_boundary_faces,
                                                                          d_abs_boundary_ispec.cuda,
                                                                          d_abs_boundary_npoin.cuda,
                                                                          d_abs_boundary_ijk.cuda,
                                                                          d_abs_boundary_jacobian2Dw.cuda,
                                                                          mp->d_ibool_outer_core.cuda,
                                                                          mp->d_vp_outer_core.cuda,
                                                                          mp->save_stacey,
                                                                          d_absorb_potential.cuda);
  }
#endif
#ifdef USE_HIP
  if (run_hip) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(compute_stacey_acoustic_kernel), grid, threads, 0, mp->compute_stream,
                                                                        mp->d_veloc_outer_core.hip,
                                                                        mp->d_accel_outer_core.hip,
                                                                        num_abs_boundary_faces,
                                                                        d_abs_boundary_ispec.hip,
                                                                        d_abs_boundary_npoin.hip,
                                                                        d_abs_boundary_ijk.hip,
                                                                        d_abs_boundary_jacobian2Dw.hip,
                                                                        mp->d_ibool_outer_core.hip,
                                                                        mp->d_vp_outer_core.hip,
                                                                        mp->save_stacey,
                                                                        d_absorb_potential.hip);
  }
#endif

  //  adjoint simulations: stores absorbed wavefield part
  if (mp->save_stacey) {
    // explicitly waits until kernel is finished
    gpuSynchronize();
    // copies array to CPU
    gpuCopy_from_device_realw (&d_absorb_potential, absorb_potential, NGLL2 * num_abs_boundary_faces);
  }

  GPU_ERROR_CHECKING ("compute_stacey_acoustic_kernel");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_stacey_acoustic_backward_gpu,
               COMPUTE_STACEY_ACOUSTIC_BACKWARD_GPU) (long *Mesh_pointer_f,
                                                      realw *absorb_potential) {
  TRACE ("compute_stacey_acoustic_backward_gpu");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // absorbing boundary
  int num_abs_boundary_faces = mp->num_abs_boundary_faces_outer_core;

  // checks if anything to do
  if (num_abs_boundary_faces == 0) return;

  gpu_int_mem d_abs_boundary_ispec;
  gpu_int_mem d_abs_boundary_npoin;
  gpu_int_mem d_abs_boundary_ijk;
  gpu_realw_mem d_absorb_potential;

  d_abs_boundary_ispec = mp->d_abs_boundary_ispec_outer_core;
  d_abs_boundary_npoin = mp->d_abs_boundary_npoin_outer_core;
  d_abs_boundary_ijk = mp->d_abs_boundary_ijk_outer_core;

  d_absorb_potential = mp->d_absorb_buffer_outer_core;

  // way 1: Elapsed time: 4.385948e-03
  // > NGLLSQUARE==NGLL2==25, but we handle this inside kernel
  //  int blocksize = 32;

  // way 2: Elapsed time: 4.379034e-03
  // > NGLLSQUARE==NGLL2==25, no further check inside kernel
  int blocksize = NGLL2;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (num_abs_boundary_faces, &num_blocks_x, &num_blocks_y);

  // adjoint simulations: needs absorbing boundary buffer
  // copies array to GPU
  gpuCopy_todevice_realw (&d_absorb_potential, absorb_potential, NGLL2 * num_abs_boundary_faces);

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_backward_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_accel_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_backward_kernel, idx++, sizeof (cl_mem), (void *) &d_absorb_potential.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_backward_kernel, idx++, sizeof (int), (void *) &num_abs_boundary_faces));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_backward_kernel, idx++, sizeof (cl_mem), (void *) &d_abs_boundary_ispec.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_backward_kernel, idx++, sizeof (cl_mem), (void *) &d_abs_boundary_npoin.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_backward_kernel, idx++, sizeof (cl_mem), (void *) &d_abs_boundary_ijk.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_backward_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_outer_core.ocl));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_stacey_acoustic_backward_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    compute_stacey_acoustic_backward_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_accel_outer_core.cuda,
                                                                                   d_absorb_potential.cuda,
                                                                                   num_abs_boundary_faces,
                                                                                   d_abs_boundary_ispec.cuda,
                                                                                   d_abs_boundary_npoin.cuda,
                                                                                   d_abs_boundary_ijk.cuda,
                                                                                   mp->d_ibool_outer_core.cuda);
  }
#endif
#ifdef USE_HIP
  if (run_hip) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(compute_stacey_acoustic_backward_kernel), grid, threads, 0, mp->compute_stream,
                                                                                 mp->d_b_accel_outer_core.hip,
                                                                                 d_absorb_potential.hip,
                                                                                 num_abs_boundary_faces,
                                                                                 d_abs_boundary_ispec.hip,
                                                                                 d_abs_boundary_npoin.hip,
                                                                                 d_abs_boundary_ijk.hip,
                                                                                 mp->d_ibool_outer_core.hip);
  }
#endif

  GPU_ERROR_CHECKING ("compute_stacey_acoustic_backward_kernel");
}


/*----------------------------------------------------------------------------------------------- */

// undo_attenuation simulation: Stacey for backward/reconstructed wavefield

/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_stacey_acoustic_undoatt_gpu,
               COMPUTE_STACEY_ACOUSTIC_UNDOATT_GPU) (long *Mesh_pointer_f) {

  TRACE ("compute_stacey_acoustic_undoatt_gpu");
  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // checks if anything to do
  if (mp->simulation_type != 3 || mp->save_forward) return;

  // absorbing boundary
  int num_abs_boundary_faces = mp->num_abs_boundary_faces_outer_core;

  // checks if anything to do
  if (num_abs_boundary_faces == 0) return;

  gpu_int_mem d_abs_boundary_ispec;
  gpu_int_mem d_abs_boundary_npoin;
  gpu_int_mem d_abs_boundary_ijk;
  gpu_realw_mem d_abs_boundary_jacobian2Dw;

  d_abs_boundary_ispec = mp->d_abs_boundary_ispec_outer_core;
  d_abs_boundary_npoin = mp->d_abs_boundary_npoin_outer_core;
  d_abs_boundary_ijk = mp->d_abs_boundary_ijk_outer_core;
  d_abs_boundary_jacobian2Dw = mp->d_abs_boundary_jacobian2Dw_outer_core;

  // way 1: Elapsed time: 4.385948e-03
  // > NGLLSQUARE==NGLL2==25, but we handle this inside kernel
  //  int blocksize = 32;

  // way 2: Elapsed time: 4.379034e-03
  // > NGLLSQUARE==NGLL2==25, no further check inside kernel
  int blocksize = NGLL2;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (num_abs_boundary_faces, &num_blocks_x, &num_blocks_y);

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_veloc_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_accel_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (int), (void *) &num_abs_boundary_faces));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &d_abs_boundary_ispec.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &d_abs_boundary_npoin.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &d_abs_boundary_ijk.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &d_abs_boundary_jacobian2Dw.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_vp_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (int), (void *) &mp->save_stacey));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) NULL));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_stacey_acoustic_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    compute_stacey_acoustic_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_veloc_outer_core.cuda,
                                                                          mp->d_b_accel_outer_core.cuda,
                                                                          num_abs_boundary_faces,
                                                                          d_abs_boundary_ispec.cuda,
                                                                          d_abs_boundary_npoin.cuda,
                                                                          d_abs_boundary_ijk.cuda,
                                                                          d_abs_boundary_jacobian2Dw.cuda,
                                                                          mp->d_ibool_outer_core.cuda,
                                                                          mp->d_vp_outer_core.cuda,
                                                                          mp->save_stacey,
                                                                          NULL);
  }
#endif
#ifdef USE_HIP
  if (run_hip) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    hipLaunchKernelGGL(HIP_KERNEL_NAME(compute_stacey_acoustic_kernel), grid, threads, 0, mp->compute_stream,
                                                                        mp->d_b_veloc_outer_core.hip,
                                                                        mp->d_b_accel_outer_core.hip,
                                                                        num_abs_boundary_faces,
                                                                        d_abs_boundary_ispec.hip,
                                                                        d_abs_boundary_npoin.hip,
                                                                        d_abs_boundary_ijk.hip,
                                                                        d_abs_boundary_jacobian2Dw.hip,
                                                                        mp->d_ibool_outer_core.hip,
                                                                        mp->d_vp_outer_core.hip,
                                                                        mp->save_stacey,
                                                                        NULL);
  }
#endif

  GPU_ERROR_CHECKING ("compute_stacey_acoustic_undoatt_gpu");
}
