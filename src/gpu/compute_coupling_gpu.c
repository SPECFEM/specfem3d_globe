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
// ACOUSTIC - ELASTIC coupling
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_coupling_fluid_cmb_gpu,
               COMPUTE_COUPLING_FLUID_CMB_GPU) (long *Mesh_pointer_f,
                                                int *FORWARD_OR_ADJOINT) {

  TRACE ("compute_coupling_fluid_cmb_gpu");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in compute_coupling_fluid_cmb_gpu() routine");
  }

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->nspec2D_top_outer_core, &num_blocks_x, &num_blocks_y);

  // sets gpu arrays
  gpu_realw_mem displ,accel;
  if (*FORWARD_OR_ADJOINT == 1) {
    displ = mp->d_displ_crust_mantle;
    accel = mp->d_accel_outer_core;
  } else {
    // debug
    DEBUG_BACKWARD_COUPLING ();
    // for backward/reconstructed fields
    displ = mp->d_b_displ_crust_mantle;
    accel = mp->d_b_accel_outer_core;
  }

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    // launches GPU kernel
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_CMB_kernel, idx++, sizeof (cl_mem), (void *) &displ.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_CMB_kernel, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_CMB_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_CMB_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibelm_bottom_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_CMB_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_normal_top_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_CMB_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_jacobian2D_top_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_CMB_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_wgllwgll_xy.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_CMB_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_CMB_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibelm_top_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_CMB_kernel, idx++, sizeof (int), (void *) &mp->nspec2D_top_outer_core));

    local_work_size[0] = NGLLX;
    local_work_size[1] = NGLLX;
    global_work_size[0] = num_blocks_x * NGLLX;
    global_work_size[1] = num_blocks_y * NGLLX;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_coupling_fluid_CMB_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(NGLLX,NGLLX,1);

    // launches GPU kernel
    compute_coupling_fluid_CMB_kernel<<<grid,threads,0,mp->compute_stream>>>(displ.cuda,
                                                                             accel.cuda,
                                                                             mp->d_ibool_crust_mantle.cuda,
                                                                             mp->d_ibelm_bottom_crust_mantle.cuda,
                                                                             mp->d_normal_top_outer_core.cuda,
                                                                             mp->d_jacobian2D_top_outer_core.cuda,
                                                                             mp->d_wgllwgll_xy.cuda,
                                                                             mp->d_ibool_outer_core.cuda,
                                                                             mp->d_ibelm_top_outer_core.cuda,
                                                                             mp->nspec2D_top_outer_core);
  }
#endif

  GPU_ERROR_CHECKING ("compute_coupling_fluid_CMB_kernel");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_coupling_fluid_icb_gpu,
               COMPUTE_COUPLING_FLUID_ICB_GPU) (long *Mesh_pointer_f,
                                                int *FORWARD_OR_ADJOINT) {

  TRACE ("compute_coupling_fluid_icb_gpu");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in compute_coupling_fluid_icb_gpu() routine");
  }

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->nspec2D_bottom_outer_core, &num_blocks_x, &num_blocks_y);

  // sets gpu arrays
  gpu_realw_mem displ,accel;
  if (*FORWARD_OR_ADJOINT == 1) {
    displ = mp->d_displ_inner_core;
    accel = mp->d_accel_outer_core;
  } else {
    // debug
    DEBUG_BACKWARD_COUPLING ();
    // for backward/reconstructed fields
    displ = mp->d_b_displ_inner_core;
    accel = mp->d_b_accel_outer_core;
  }

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    // launches GPU kernel
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_ICB_kernel, idx++, sizeof (cl_mem), (void *) &displ.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_ICB_kernel, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_ICB_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_ICB_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibelm_top_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_ICB_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_normal_bottom_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_ICB_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_jacobian2D_bottom_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_ICB_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_wgllwgll_xy.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_ICB_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_ICB_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibelm_bottom_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_fluid_ICB_kernel, idx++, sizeof (int), (void *) &mp->nspec2D_bottom_outer_core));

    local_work_size[0] = NGLLX;
    local_work_size[1] = NGLLX;
    global_work_size[0] = num_blocks_x * NGLLX;
    global_work_size[1] = num_blocks_y * NGLLX;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_coupling_fluid_ICB_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(NGLLX,NGLLX,1);

    // launches GPU kernel
    compute_coupling_fluid_ICB_kernel<<<grid,threads,0,mp->compute_stream>>>(displ.cuda,
                                                                             accel.cuda,
                                                                             mp->d_ibool_inner_core.cuda,
                                                                             mp->d_ibelm_top_inner_core.cuda,
                                                                             mp->d_normal_bottom_outer_core.cuda,
                                                                             mp->d_jacobian2D_bottom_outer_core.cuda,
                                                                             mp->d_wgllwgll_xy.cuda,
                                                                             mp->d_ibool_outer_core.cuda,
                                                                             mp->d_ibelm_bottom_outer_core.cuda,
                                                                             mp->nspec2D_bottom_outer_core);
  }
#endif

  GPU_ERROR_CHECKING ("compute_coupling_fluid_ICB_kernel");
}

/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_coupling_cmb_fluid_gpu,
               COMPUTE_COUPLING_CMB_FLUID_GPU) (long *Mesh_pointer_f,
                                                int *FORWARD_OR_ADJOINT) {

  TRACE ("compute_coupling_cmb_fluid_gpu");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in compute_coupling_cmb_fluid_gpu() routine");
  }

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->nspec2D_bottom_crust_mantle, &num_blocks_x, &num_blocks_y);

  // sets gpu arrays
  gpu_realw_mem displ,accel_cm,accel_oc;
  if (*FORWARD_OR_ADJOINT == 1) {
    displ = mp->d_displ_crust_mantle;
    accel_cm = mp->d_accel_crust_mantle;
    accel_oc = mp->d_accel_outer_core;
  } else {
    // debug
    DEBUG_BACKWARD_COUPLING ();
    // for backward/reconstructed fields
    displ = mp->d_b_displ_crust_mantle;
    accel_cm = mp->d_b_accel_crust_mantle;
    accel_oc = mp->d_b_accel_outer_core;
  }

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    // launches GPU kernel
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_CMB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &displ.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_CMB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &accel_cm.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_CMB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &accel_oc.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_CMB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_CMB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibelm_bottom_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_CMB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_normal_top_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_CMB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_jacobian2D_top_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_CMB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_wgllwgll_xy.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_CMB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_CMB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibelm_top_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_CMB_fluid_kernel, idx++, sizeof (realw), (void *) &mp->RHO_TOP_OC));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_CMB_fluid_kernel, idx++, sizeof (realw), (void *) &mp->minus_g_cmb));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_CMB_fluid_kernel, idx++, sizeof (int), (void *) &mp->gravity));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_CMB_fluid_kernel, idx++, sizeof (int), (void *) &mp->nspec2D_bottom_crust_mantle));

    local_work_size[0] = 5;
    local_work_size[1] = 5;
    global_work_size[0] = num_blocks_x * 5;
    global_work_size[1] = num_blocks_y * 5;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_coupling_CMB_fluid_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(5,5,1);

    // launches GPU kernel
    compute_coupling_CMB_fluid_kernel<<<grid,threads,0,mp->compute_stream>>>(displ.cuda,
                                                                             accel_cm.cuda,
                                                                             accel_oc.cuda,
                                                                             mp->d_ibool_crust_mantle.cuda,
                                                                             mp->d_ibelm_bottom_crust_mantle.cuda,
                                                                             mp->d_normal_top_outer_core.cuda,
                                                                             mp->d_jacobian2D_top_outer_core.cuda,
                                                                             mp->d_wgllwgll_xy.cuda,
                                                                             mp->d_ibool_outer_core.cuda,
                                                                             mp->d_ibelm_top_outer_core.cuda,
                                                                             mp->RHO_TOP_OC,
                                                                             mp->minus_g_cmb,
                                                                             mp->gravity,
                                                                             mp->nspec2D_bottom_crust_mantle);
  }
#endif

  GPU_ERROR_CHECKING ("compute_coupling_CMB_fluid_gpu");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_coupling_icb_fluid_gpu,
               COMPUTE_COUPLING_ICB_FLUID_GPU) (long *Mesh_pointer_f,
                                                int *FORWARD_OR_ADJOINT) {
  TRACE ("compute_coupling_icb_fluid_gpu");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in compute_coupling_icb_fluid_gpu() routine");
  }

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->nspec2D_top_inner_core, &num_blocks_x, &num_blocks_y);

  // sets gpu arrays
  gpu_realw_mem displ,accel_ic,accel_oc;
  if (*FORWARD_OR_ADJOINT == 1) {
    displ = mp->d_displ_inner_core;
    accel_ic = mp->d_accel_inner_core;
    accel_oc = mp->d_accel_outer_core;
  } else {
    // debug
    DEBUG_BACKWARD_COUPLING ();
    // for backward/reconstructed fields
    displ = mp->d_b_displ_inner_core;
    accel_ic = mp->d_b_accel_inner_core;
    accel_oc = mp->d_b_accel_outer_core;
  }

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    // launches GPU kernel
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ICB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &displ.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ICB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &accel_ic.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ICB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &accel_oc.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ICB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ICB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibelm_top_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ICB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_normal_bottom_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ICB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_jacobian2D_bottom_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ICB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_wgllwgll_xy.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ICB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ICB_fluid_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibelm_bottom_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ICB_fluid_kernel, idx++, sizeof (realw), (void *) &mp->RHO_BOTTOM_OC));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ICB_fluid_kernel, idx++, sizeof (realw), (void *) &mp->minus_g_icb));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ICB_fluid_kernel, idx++, sizeof (int), (void *) &mp->gravity));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ICB_fluid_kernel, idx++, sizeof (int), (void *) &mp->nspec2D_top_inner_core));

    local_work_size[0] = NGLLX;
    local_work_size[1] = NGLLX;
    global_work_size[0] = num_blocks_x * NGLLX;
    global_work_size[1] = num_blocks_y * NGLLX;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_coupling_ICB_fluid_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(NGLLX,NGLLX,1);

    // launches GPU kernel
    compute_coupling_ICB_fluid_kernel<<<grid,threads,0,mp->compute_stream>>>(displ.cuda,
                                                                             accel_ic.cuda,
                                                                             accel_oc.cuda,
                                                                             mp->d_ibool_inner_core.cuda,
                                                                             mp->d_ibelm_top_inner_core.cuda,
                                                                             mp->d_normal_bottom_outer_core.cuda,
                                                                             mp->d_jacobian2D_bottom_outer_core.cuda,
                                                                             mp->d_wgllwgll_xy.cuda,
                                                                             mp->d_ibool_outer_core.cuda,
                                                                             mp->d_ibelm_bottom_outer_core.cuda,
                                                                             mp->RHO_BOTTOM_OC,
                                                                             mp->minus_g_icb,
                                                                             mp->gravity,
                                                                             mp->nspec2D_top_inner_core);
  }
#endif

  GPU_ERROR_CHECKING ("compute_coupling_ICB_fluid_gpu");
}

/*----------------------------------------------------------------------------------------------- */
/*OCEANS load coupled on free surface */
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_coupling_ocean_gpu,
               COMPUTE_COUPLING_OCEAN_GPU) (long *Mesh_pointer_f,
                                            int *FORWARD_OR_ADJOINT) {

  TRACE ("compute_coupling_ocean_gpu");
  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in compute_coupling_ocean_gpu() routine");
  }

  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int) ceil (((double) mp->npoin_oceans)/ ((double) blocksize))) * blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

  // sets gpu arrays
  gpu_realw_mem accel,rmassx,rmassy,rmassz;
  if (*FORWARD_OR_ADJOINT == 1) {
    accel = mp->d_accel_crust_mantle;
    rmassx = mp->d_rmassx_crust_mantle;
    rmassy = mp->d_rmassy_crust_mantle;
    rmassz = mp->d_rmassz_crust_mantle;
  } else {
    // debug
    DEBUG_BACKWARD_COUPLING ();
    // for backward/reconstructed fields
    accel = mp->d_b_accel_crust_mantle;
    rmassx = mp->d_b_rmassx_crust_mantle;
    rmassy = mp->d_b_rmassy_crust_mantle;
    rmassz = mp->d_b_rmassz_crust_mantle;
  }

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    // uses corrected mass matrices
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ocean_kernel, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ocean_kernel, idx++, sizeof (cl_mem), (void *) &rmassx.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ocean_kernel, idx++, sizeof (cl_mem), (void *) &rmassy.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ocean_kernel, idx++, sizeof (cl_mem), (void *) &rmassz.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ocean_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_rmass_ocean_load.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ocean_kernel, idx++, sizeof (int), (void *) &mp->npoin_oceans));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ocean_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_ocean_load.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_coupling_ocean_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_normal_ocean_load.ocl));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_coupling_ocean_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));

  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // uses corrected mass matrices
    compute_coupling_ocean_kernel<<<grid,threads,0,mp->compute_stream>>>(accel.cuda,
                                                                         rmassx.cuda,
                                                                         rmassy.cuda,
                                                                         rmassz.cuda,
                                                                         mp->d_rmass_ocean_load.cuda,
                                                                         mp->npoin_oceans,
                                                                         mp->d_ibool_ocean_load.cuda,
                                                                         mp->d_normal_ocean_load.cuda);
  }
#endif

  GPU_ERROR_CHECKING ("compute_coupling_ocean_gpu");
}
