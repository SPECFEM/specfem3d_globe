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

extern EXTERN_LANG
void FC_FUNC_ (compute_stacey_acoustic_gpu,
               COMPUTE_STACEY_ACOUSTIC_GPU) (long *Mesh_pointer_f,
                                             realw *absorb_potential,
                                             int *itype) {
  TRACE ("compute_stacey_acoustic_gpu");
  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int num_abs_boundary_faces = 0;

  gpu_int_mem d_abs_boundary_ispec;
  gpu_realw_mem d_abs_boundary_jacobian2D;
  gpu_realw_mem d_wgllwgll;
  gpu_realw_mem d_b_absorb_potential;

  // absorbing boundary type
  int interface_type = *itype;
  switch (interface_type) {
  case 4:
    // xmin
    num_abs_boundary_faces = mp->nspec2D_xmin_outer_core;
    d_abs_boundary_ispec = mp->d_ibelm_xmin_outer_core;
    d_abs_boundary_jacobian2D = mp->d_jacobian2D_xmin_outer_core;
    d_b_absorb_potential = mp->d_absorb_xmin_outer_core;
    d_wgllwgll = mp->d_wgllwgll_yz;
    break;

  case 5:
    // xmax
    num_abs_boundary_faces = mp->nspec2D_xmax_outer_core;
    d_abs_boundary_ispec = mp->d_ibelm_xmax_outer_core;
    d_abs_boundary_jacobian2D = mp->d_jacobian2D_xmax_outer_core;
    d_b_absorb_potential = mp->d_absorb_xmax_outer_core;
    d_wgllwgll = mp->d_wgllwgll_yz;
    break;

  case 6:
    // ymin
    num_abs_boundary_faces = mp->nspec2D_ymin_outer_core;
    d_abs_boundary_ispec = mp->d_ibelm_ymin_outer_core;
    d_abs_boundary_jacobian2D = mp->d_jacobian2D_ymin_outer_core;
    d_b_absorb_potential = mp->d_absorb_ymin_outer_core;
    d_wgllwgll = mp->d_wgllwgll_xz;
    break;

  case 7:
    // ymax
    num_abs_boundary_faces = mp->nspec2D_ymax_outer_core;
    d_abs_boundary_ispec = mp->d_ibelm_ymax_outer_core;
    d_abs_boundary_jacobian2D = mp->d_jacobian2D_ymax_outer_core;
    d_b_absorb_potential = mp->d_absorb_ymax_outer_core;
    d_wgllwgll = mp->d_wgllwgll_xz;
    break;

  case 8:
    // zmin
    num_abs_boundary_faces = mp->nspec2D_zmin_outer_core;
    d_abs_boundary_ispec = mp->d_ibelm_bottom_outer_core;
    d_abs_boundary_jacobian2D = mp->d_jacobian2D_bottom_outer_core;
    d_b_absorb_potential = mp->d_absorb_zmin_outer_core;
    d_wgllwgll = mp->d_wgllwgll_xy;
    break;

  default:
    exit_on_error ("compute_stacey_acoustic_gpu: unknown interface type");
  }

  // checks if anything to do
  if (num_abs_boundary_faces == 0) return;

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
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (int), (void *) &interface_type));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (int), (void *) &num_abs_boundary_faces));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &d_abs_boundary_ispec.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_nkmin_xi_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_nkmin_eta_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_njmin_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_njmax_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_nimin_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_nimax_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &d_abs_boundary_jacobian2D.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &d_wgllwgll.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_vp_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (int), (void *) &mp->save_forward));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &d_b_absorb_potential.ocl));

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
                                                     interface_type,
                                                     num_abs_boundary_faces,
                                                     d_abs_boundary_ispec.cuda,
                                                     mp->d_nkmin_xi_outer_core.cuda,
                                                     mp->d_nkmin_eta_outer_core.cuda,
                                                     mp->d_njmin_outer_core.cuda,
                                                     mp->d_njmax_outer_core.cuda,
                                                     mp->d_nimin_outer_core.cuda,
                                                     mp->d_nimax_outer_core.cuda,
                                                     d_abs_boundary_jacobian2D.cuda,
                                                     d_wgllwgll.cuda,
                                                     mp->d_ibool_outer_core.cuda,
                                                     mp->d_vp_outer_core.cuda,
                                                     mp->save_forward,
                                                     d_b_absorb_potential.cuda);
  }
#endif

  //  adjoint simulations: stores absorbed wavefield part
  if (mp->save_forward) {
    // explicitly waits until kernel is finished
    gpuSynchronize();
    // copies array to CPU
    gpuCopy_from_device_realw (&d_b_absorb_potential, absorb_potential, NGLL2 * num_abs_boundary_faces);
  }

  GPU_ERROR_CHECKING ("compute_stacey_acoustic_kernel");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_stacey_acoustic_backward_gpu,
               COMPUTE_STACEY_ACOUSTIC_BACKWARD_GPU) (long *Mesh_pointer_f,
                                                      realw *absorb_potential,
                                                      int *itype) {
  TRACE ("compute_stacey_acoustic_backward_gpu");

  int num_abs_boundary_faces = 0;

  gpu_int_mem d_abs_boundary_ispec;
  gpu_realw_mem d_b_absorb_potential;

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // absorbing boundary type
  int interface_type = *itype;
  switch (interface_type) {
  case 4:
    // xmin
    num_abs_boundary_faces = mp->nspec2D_xmin_outer_core;
    d_abs_boundary_ispec = mp->d_ibelm_xmin_outer_core;
    d_b_absorb_potential = mp->d_absorb_xmin_outer_core;
    break;

  case 5:
    // xmax
    num_abs_boundary_faces = mp->nspec2D_xmax_outer_core;
    d_abs_boundary_ispec = mp->d_ibelm_xmax_outer_core;
    d_b_absorb_potential = mp->d_absorb_xmax_outer_core;
    break;

  case 6:
    // ymin
    num_abs_boundary_faces = mp->nspec2D_ymin_outer_core;
    d_abs_boundary_ispec = mp->d_ibelm_ymin_outer_core;
    d_b_absorb_potential = mp->d_absorb_ymin_outer_core;
    break;

  case 7:
    // ymax
    num_abs_boundary_faces = mp->nspec2D_ymax_outer_core;
    d_abs_boundary_ispec = mp->d_ibelm_ymax_outer_core;
    d_b_absorb_potential = mp->d_absorb_ymax_outer_core;
    break;

  case 8:
    // zmin
    num_abs_boundary_faces = mp->nspec2D_zmin_outer_core;
    d_abs_boundary_ispec = mp->d_ibelm_bottom_outer_core;
    d_b_absorb_potential = mp->d_absorb_zmin_outer_core;
    break;

  default:
    exit_on_error ("compute_stacey_acoustic_backward_gpu: unknown interface type");
  }

  // checks if anything to do
  if (num_abs_boundary_faces == 0) return;

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
  gpuCopy_todevice_realw (&d_b_absorb_potential, absorb_potential, NGLL2 * num_abs_boundary_faces);

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_backward_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_accel_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_backward_kernel, idx++, sizeof (cl_mem), (void *) &d_b_absorb_potential.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_backward_kernel, idx++, sizeof (int), (void *) &interface_type));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_backward_kernel, idx++, sizeof (int), (void *) &num_abs_boundary_faces));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_backward_kernel, idx++, sizeof (cl_mem), (void *) &d_abs_boundary_ispec.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_backward_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_nkmin_xi_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_backward_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_nkmin_eta_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_backward_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_njmin_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_backward_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_njmax_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_backward_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_nimin_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_backward_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_nimax_outer_core.ocl));
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
                                                                                   d_b_absorb_potential.cuda,
                                                                                   interface_type,
                                                                                   num_abs_boundary_faces,
                                                                                   d_abs_boundary_ispec.cuda,
                                                                                   mp->d_nkmin_xi_outer_core.cuda,mp->d_nkmin_eta_outer_core.cuda,
                                                                                   mp->d_njmin_outer_core.cuda,mp->d_njmax_outer_core.cuda,
                                                                                   mp->d_nimin_outer_core.cuda,mp->d_nimax_outer_core.cuda,
                                                                                   mp->d_ibool_outer_core.cuda);
  }
#endif

  GPU_ERROR_CHECKING ("compute_stacey_acoustic_backward_kernel");
}


/*----------------------------------------------------------------------------------------------- */

// undo_attenuation simulation: Stacey for backward/reconstructed wavefield

/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_stacey_acoustic_undoatt_gpu,
               COMPUTE_STACEY_ACOUSTIC_UNDOATT_GPU) (long *Mesh_pointer_f,
                                                     int *itype) {
  TRACE ("compute_stacey_acoustic_undoatt_gpu");
  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int num_abs_boundary_faces = 0;

  gpu_int_mem d_abs_boundary_ispec;
  gpu_realw_mem d_abs_boundary_jacobian2D;
  gpu_realw_mem d_wgllwgll;

  // checks if anything to do
  if (mp->simulation_type /= 3 || mp->save_forward) return;

  // absorbing boundary type
  int interface_type = *itype;
  switch (interface_type) {
  case 4:
    // xmin
    num_abs_boundary_faces = mp->nspec2D_xmin_outer_core;
    d_abs_boundary_ispec = mp->d_ibelm_xmin_outer_core;
    d_abs_boundary_jacobian2D = mp->d_jacobian2D_xmin_outer_core;
    d_wgllwgll = mp->d_wgllwgll_yz;
    break;

  case 5:
    // xmax
    num_abs_boundary_faces = mp->nspec2D_xmax_outer_core;
    d_abs_boundary_ispec = mp->d_ibelm_xmax_outer_core;
    d_abs_boundary_jacobian2D = mp->d_jacobian2D_xmax_outer_core;
    d_wgllwgll = mp->d_wgllwgll_yz;
    break;

  case 6:
    // ymin
    num_abs_boundary_faces = mp->nspec2D_ymin_outer_core;
    d_abs_boundary_ispec = mp->d_ibelm_ymin_outer_core;
    d_abs_boundary_jacobian2D = mp->d_jacobian2D_ymin_outer_core;
    d_wgllwgll = mp->d_wgllwgll_xz;
    break;

  case 7:
    // ymax
    num_abs_boundary_faces = mp->nspec2D_ymax_outer_core;
    d_abs_boundary_ispec = mp->d_ibelm_ymax_outer_core;
    d_abs_boundary_jacobian2D = mp->d_jacobian2D_ymax_outer_core;
    d_wgllwgll = mp->d_wgllwgll_xz;
    break;

  case 8:
    // zmin
    num_abs_boundary_faces = mp->nspec2D_zmin_outer_core;
    d_abs_boundary_ispec = mp->d_ibelm_bottom_outer_core;
    d_abs_boundary_jacobian2D = mp->d_jacobian2D_bottom_outer_core;
    d_wgllwgll = mp->d_wgllwgll_xy;
    break;

  default:
    exit_on_error ("compute_stacey_acoustic_undoatt_gpu: unknown interface type");
  }

  // checks if anything to do
  if (num_abs_boundary_faces == 0) return;

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
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (int), (void *) &interface_type));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (int), (void *) &num_abs_boundary_faces));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &d_abs_boundary_ispec.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_nkmin_xi_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_nkmin_eta_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_njmin_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_njmax_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_nimin_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_nimax_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &d_abs_boundary_jacobian2D.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &d_wgllwgll.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_vp_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_stacey_acoustic_kernel, idx++, sizeof (int), (void *) &mp->save_forward));
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
                                                                          interface_type,
                                                                          num_abs_boundary_faces,
                                                                          d_abs_boundary_ispec.cuda,
                                                                          mp->d_nkmin_xi_outer_core.cuda,
                                                                          mp->d_nkmin_eta_outer_core.cuda,
                                                                          mp->d_njmin_outer_core.cuda,
                                                                          mp->d_njmax_outer_core.cuda,
                                                                          mp->d_nimin_outer_core.cuda,
                                                                          mp->d_nimax_outer_core.cuda,
                                                                          d_abs_boundary_jacobian2D.cuda,
                                                                          d_wgllwgll.cuda,
                                                                          mp->d_ibool_outer_core.cuda,
                                                                          mp->d_vp_outer_core.cuda,
                                                                          mp->save_forward,
                                                                          NULL);
  }
#endif

  GPU_ERROR_CHECKING ("compute_stacey_acoustic_undoatt_gpu");
}
