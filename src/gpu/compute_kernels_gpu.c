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

/* ----------------------------------------------------------------------------------------------- */
// attention: compute anisotropic kernels , d_b_displ_crust_mantle
extern EXTERN_LANG
void FC_FUNC_ (compute_kernels_cm_gpu,
               COMPUTE_KERNELS_CM_GPU) (long *Mesh_pointer_f, realw *deltat_f) {

  TRACE ("compute_cm_gpu");
  // debug
  DEBUG_BACKWARD_KERNEL ();

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int blocksize = NGLL3;
  realw deltat = *deltat_f;

  // blocks
  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->NSPEC_CRUST_MANTLE, &num_blocks_x, &num_blocks_y);

  // density kernel
#ifdef USE_OPENCL
  size_t global_work_size[2];
  size_t local_work_size[2];
  cl_uint idx = 0;
  if (run_opencl) {
    // dimensions
    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    // rho kernel
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_accel_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_displ_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_rho_kl_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (int), (void *) &mp->NSPEC_CRUST_MANTLE));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (realw), (void *) &deltat));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_rho_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);
  if (run_cuda) {
    // rho kernel
    compute_rho_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_ibool_crust_mantle.cuda,
                                                              mp->d_accel_crust_mantle.cuda,
                                                              mp->d_b_displ_crust_mantle.cuda,
                                                              mp->d_rho_kl_crust_mantle.cuda,
                                                              mp->NSPEC_CRUST_MANTLE,
                                                              deltat);
  }
#endif

  // wavespeed kernels
  // checks if strain is available
  if (mp->undo_attenuation) {
    // checks strain array size
    if (mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY == 1) {
      exit_on_error("compute_kernels_cm_cuda NSPEC_CRUST_MANTLE_STRAIN_ONLY invalid with undo_att");
    }

#ifdef USE_OPENCL
    if (run_opencl) {
      idx = 0;
      if (!mp->anisotropic_kl) {
        // isotropic kernels
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xx_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_yy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_yz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_eps_trace_over_3_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_beta_kl_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_alpha_kl_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (int), (void *) &mp->NSPEC_CRUST_MANTLE));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (realw), (void *) &deltat));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_displ_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xix_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xiy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xiz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etax_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etay_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etaz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammax_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammay_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammaz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_hprime_xx.ocl));

        clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_iso_undoatt_kernel, 2, NULL,
                                         global_work_size, local_work_size, 0, NULL, NULL));
      } else {
        // anisotropic kernels
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xx_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_yy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_yz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_eps_trace_over_3_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_cijkl_kl_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (int), (void *) &mp->NSPEC_CRUST_MANTLE));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (realw), (void *) &deltat));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_displ_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xix_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xiy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xiz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etax_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etay_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etaz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammax_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammay_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammaz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_hprime_xx.ocl));

        clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_ani_undoatt_kernel, 2, NULL,
                                         global_work_size, local_work_size, 0, NULL, NULL));
      }
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      // computes strain locally based on current backward/reconstructed (b_displ) wavefield
      if (! mp->anisotropic_kl) {
        // isotropic kernels
        compute_iso_undoatt_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_epsilondev_xx_crust_mantle.cuda,
                                                                          mp->d_epsilondev_yy_crust_mantle.cuda,
                                                                          mp->d_epsilondev_xy_crust_mantle.cuda,
                                                                          mp->d_epsilondev_xz_crust_mantle.cuda,
                                                                          mp->d_epsilondev_yz_crust_mantle.cuda,
                                                                          mp->d_eps_trace_over_3_crust_mantle.cuda,
                                                                          mp->d_beta_kl_crust_mantle.cuda,
                                                                          mp->d_alpha_kl_crust_mantle.cuda,
                                                                          mp->NSPEC_CRUST_MANTLE,
                                                                          deltat,
                                                                          mp->d_ibool_crust_mantle.cuda,
                                                                          mp->d_b_displ_crust_mantle.cuda,
                                                                          mp->d_xix_crust_mantle.cuda,mp->d_xiy_crust_mantle.cuda,mp->d_xiz_crust_mantle.cuda,
                                                                          mp->d_etax_crust_mantle.cuda,mp->d_etay_crust_mantle.cuda,mp->d_etaz_crust_mantle.cuda,
                                                                          mp->d_gammax_crust_mantle.cuda,mp->d_gammay_crust_mantle.cuda,mp->d_gammaz_crust_mantle.cuda,
                                                                          mp->d_hprime_xx.cuda);
      } else {
        // anisotropic kernels
        compute_ani_undoatt_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_epsilondev_xx_crust_mantle.cuda,
                                                                          mp->d_epsilondev_yy_crust_mantle.cuda,
                                                                          mp->d_epsilondev_xy_crust_mantle.cuda,
                                                                          mp->d_epsilondev_xz_crust_mantle.cuda,
                                                                          mp->d_epsilondev_yz_crust_mantle.cuda,
                                                                          mp->d_eps_trace_over_3_crust_mantle.cuda,
                                                                          mp->d_cijkl_kl_crust_mantle.cuda,
                                                                          mp->NSPEC_CRUST_MANTLE,
                                                                          deltat,
                                                                          mp->d_ibool_crust_mantle.cuda,
                                                                          mp->d_b_displ_crust_mantle.cuda,
                                                                          mp->d_xix_crust_mantle.cuda,mp->d_xiy_crust_mantle.cuda,mp->d_xiz_crust_mantle.cuda,
                                                                          mp->d_etax_crust_mantle.cuda,mp->d_etay_crust_mantle.cuda,mp->d_etaz_crust_mantle.cuda,
                                                                          mp->d_gammax_crust_mantle.cuda,mp->d_gammay_crust_mantle.cuda,mp->d_gammaz_crust_mantle.cuda,
                                                                          mp->d_hprime_xx.cuda);
      }
    }
#endif
  } else {
    // strain available by compute_forces routine
#ifdef USE_OPENCL
    if (run_opencl) {
      idx = 0;
      if (!mp->anisotropic_kl) {
        // isotropic kernels
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xx_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_yy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_yz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_eps_trace_over_3_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_epsilondev_xx_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_epsilondev_yy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_epsilondev_xy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_epsilondev_xz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_epsilondev_yz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_eps_trace_over_3_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_beta_kl_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_alpha_kl_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (int), (void *) &mp->NSPEC_CRUST_MANTLE));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (realw), (void *) &deltat));

        clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_iso_kernel, 2, NULL,
                                         global_work_size, local_work_size, 0, NULL, NULL));
      } else {
        // anisotropic kernels
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xx_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_yy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_yz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_eps_trace_over_3_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_epsilondev_xx_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_epsilondev_yy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_epsilondev_xy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_epsilondev_xz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_epsilondev_yz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_eps_trace_over_3_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_cijkl_kl_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_kernel, idx++, sizeof (int), (void *) &mp->NSPEC_CRUST_MANTLE));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_kernel, idx++, sizeof (realw), (void *) &deltat));

        clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_ani_kernel, 2, NULL,
                                         global_work_size, local_work_size, 0, NULL, NULL));
      }
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      // takes strain arrays computed from previous compute_forces call
      if (! mp->anisotropic_kl) {
        // isotropic kernels
        compute_iso_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_epsilondev_xx_crust_mantle.cuda,
                                                                  mp->d_epsilondev_yy_crust_mantle.cuda,
                                                                  mp->d_epsilondev_xy_crust_mantle.cuda,
                                                                  mp->d_epsilondev_xz_crust_mantle.cuda,
                                                                  mp->d_epsilondev_yz_crust_mantle.cuda,
                                                                  mp->d_eps_trace_over_3_crust_mantle.cuda,
                                                                  mp->d_b_epsilondev_xx_crust_mantle.cuda,
                                                                  mp->d_b_epsilondev_yy_crust_mantle.cuda,
                                                                  mp->d_b_epsilondev_xy_crust_mantle.cuda,
                                                                  mp->d_b_epsilondev_xz_crust_mantle.cuda,
                                                                  mp->d_b_epsilondev_yz_crust_mantle.cuda,
                                                                  mp->d_b_eps_trace_over_3_crust_mantle.cuda,
                                                                  mp->d_beta_kl_crust_mantle.cuda,
                                                                  mp->d_alpha_kl_crust_mantle.cuda,
                                                                  mp->NSPEC_CRUST_MANTLE,
                                                                  deltat);
      } else {
        // anisotropic kernels
        compute_ani_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_epsilondev_xx_crust_mantle.cuda,
                                                                  mp->d_epsilondev_yy_crust_mantle.cuda,
                                                                  mp->d_epsilondev_xy_crust_mantle.cuda,
                                                                  mp->d_epsilondev_xz_crust_mantle.cuda,
                                                                  mp->d_epsilondev_yz_crust_mantle.cuda,
                                                                  mp->d_eps_trace_over_3_crust_mantle.cuda,
                                                                  mp->d_b_epsilondev_xx_crust_mantle.cuda,
                                                                  mp->d_b_epsilondev_yy_crust_mantle.cuda,
                                                                  mp->d_b_epsilondev_xy_crust_mantle.cuda,
                                                                  mp->d_b_epsilondev_xz_crust_mantle.cuda,
                                                                  mp->d_b_epsilondev_yz_crust_mantle.cuda,
                                                                  mp->d_b_eps_trace_over_3_crust_mantle.cuda,
                                                                  mp->d_cijkl_kl_crust_mantle.cuda,
                                                                  mp->NSPEC_CRUST_MANTLE,
                                                                  deltat);
      }
    }
#endif
  }

  GPU_ERROR_CHECKING ("compute_cm_gpu");
}


/*----------------------------------------------------------------------------------------------- */
// inner_core
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_kernels_ic_gpu,
               COMPUTE_KERNELS_IC_GPU) (long *Mesh_pointer_f, realw *deltat_f) {

  TRACE("compute_kernels_ic_gpu");
  // debug
  DEBUG_BACKWARD_KERNEL();

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int blocksize = NGLL3;
  realw deltat = *deltat_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->NSPEC_INNER_CORE, &num_blocks_x, &num_blocks_y);

  // only isotropic kernels in inner core so far implemented
  // density kernel
#ifdef USE_OPENCL
  size_t global_work_size[2];
  size_t local_work_size[2];
  cl_uint idx = 0;
  if (run_opencl) {
    // dimensions
    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    // rho kernel
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_accel_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_displ_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_rho_kl_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (int), (void *) &mp->NSPEC_INNER_CORE));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (realw), (void *) &deltat));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_rho_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);
  if (run_cuda) {
    // rho kernels
    compute_rho_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_ibool_inner_core.cuda,
                                                              mp->d_accel_inner_core.cuda,
                                                              mp->d_b_displ_inner_core.cuda,
                                                              mp->d_rho_kl_inner_core.cuda,
                                                              mp->NSPEC_INNER_CORE,
                                                              deltat);
  }
#endif

  // wavespeed kernels
  // checks if strain is available
  if (mp->undo_attenuation) {
    // checks strain array size
    if (mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY == 1) {
      exit_on_error("compute_kernels_ic_gpu NSPEC_CRUST_MANTLE_STRAIN_ONLY invalid with undo_att");
    }

    // computes strain locally based on current backward/reconstructed (b_displ) wavefield
    // isotropic kernels (shear, bulk)
#ifdef USE_OPENCL
    if (run_opencl) {
      idx = 0;
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xx_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_yy_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xy_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xz_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_yz_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_eps_trace_over_3_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_beta_kl_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_alpha_kl_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (int), (void *) &mp->NSPEC_INNER_CORE));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (realw), (void *) &deltat));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_displ_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xix_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xiy_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xiz_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etax_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etay_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etaz_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammax_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammay_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammaz_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_undoatt_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_hprime_xx.ocl));

      clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_iso_undoatt_kernel, 2, NULL,
                                       global_work_size, local_work_size, 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      compute_iso_undoatt_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_epsilondev_xx_inner_core.cuda,
                                                                        mp->d_epsilondev_yy_inner_core.cuda,
                                                                        mp->d_epsilondev_xy_inner_core.cuda,
                                                                        mp->d_epsilondev_xz_inner_core.cuda,
                                                                        mp->d_epsilondev_yz_inner_core.cuda,
                                                                        mp->d_eps_trace_over_3_inner_core.cuda,
                                                                        mp->d_beta_kl_inner_core.cuda,
                                                                        mp->d_alpha_kl_inner_core.cuda,
                                                                        mp->NSPEC_INNER_CORE,
                                                                        deltat,
                                                                        mp->d_ibool_inner_core.cuda,
                                                                        mp->d_b_displ_inner_core.cuda,
                                                                        mp->d_xix_inner_core.cuda,
                                                                        mp->d_xiy_inner_core.cuda,
                                                                        mp->d_xiz_inner_core.cuda,
                                                                        mp->d_etax_inner_core.cuda,
                                                                        mp->d_etay_inner_core.cuda,
                                                                        mp->d_etaz_inner_core.cuda,
                                                                        mp->d_gammax_inner_core.cuda,
                                                                        mp->d_gammay_inner_core.cuda,
                                                                        mp->d_gammaz_inner_core.cuda,
                                                                        mp->d_hprime_xx.cuda);
    }
#endif

  }else{
    // takes strain arrays computed from previous compute_forces call

    // isotropic kernels (shear, bulk)
#ifdef USE_OPENCL
    if (run_opencl) {
      idx = 0;
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xx_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_yy_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xy_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xz_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_yz_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_eps_trace_over_3_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_epsilondev_xx_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_epsilondev_yy_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_epsilondev_xy_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_epsilondev_xz_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_epsilondev_yz_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_eps_trace_over_3_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_beta_kl_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_alpha_kl_inner_core.ocl));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (int), (void *) &mp->NSPEC_INNER_CORE));
      clCheck (clSetKernelArg (mocl.kernels.compute_iso_kernel, idx++, sizeof (realw), (void *) &deltat));

      clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_iso_kernel, 2, NULL,
                                       global_work_size, local_work_size, 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      compute_iso_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_epsilondev_xx_inner_core.cuda,
                                                                mp->d_epsilondev_yy_inner_core.cuda,
                                                                mp->d_epsilondev_xy_inner_core.cuda,
                                                                mp->d_epsilondev_xz_inner_core.cuda,
                                                                mp->d_epsilondev_yz_inner_core.cuda,
                                                                mp->d_eps_trace_over_3_inner_core.cuda,
                                                                mp->d_b_epsilondev_xx_inner_core.cuda,
                                                                mp->d_b_epsilondev_yy_inner_core.cuda,
                                                                mp->d_b_epsilondev_xy_inner_core.cuda,
                                                                mp->d_b_epsilondev_xz_inner_core.cuda,
                                                                mp->d_b_epsilondev_yz_inner_core.cuda,
                                                                mp->d_b_eps_trace_over_3_inner_core.cuda,
                                                                mp->d_beta_kl_inner_core.cuda,
                                                                mp->d_alpha_kl_inner_core.cuda,
                                                                mp->NSPEC_INNER_CORE,
                                                                deltat);
    }
#endif
  }

  GPU_ERROR_CHECKING ("compute_ic_gpu");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_kernels_oc_gpu,
               COMPUTE_KERNELS_OC_GPU) (long *Mesh_pointer_f, realw *deltat_f) {

  TRACE ("compute_kernels_oc_gpu");
  // debug
  DEBUG_BACKWARD_KERNEL();

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int blocksize = NGLL3;   // NGLLX*NGLLY*NGLLZ
  realw deltat = *deltat_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->NSPEC_OUTER_CORE, &num_blocks_x, &num_blocks_y);


#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    // dimensions
    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    // acoustic kernels rho/alpha
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_rhostore_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_kappavstore_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_hprime_xx.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xix_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xiy_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xiz_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etax_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etay_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etaz_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammax_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammay_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammaz_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_accel_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_displ_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_accel_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_rho_kl_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_alpha_kl_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (realw), (void *) &deltat));
    clCheck (clSetKernelArg (mocl.kernels.compute_acoustic_kernel, idx++, sizeof (int), (void *) &mp->NSPEC_OUTER_CORE));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_acoustic_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    compute_acoustic_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_ibool_outer_core.cuda,
                                                                   mp->d_rhostore_outer_core.cuda,
                                                                   mp->d_kappavstore_outer_core.cuda,
                                                                   mp->d_hprime_xx.cuda,
                                                                   mp->d_xix_outer_core.cuda,
                                                                   mp->d_xiy_outer_core.cuda,
                                                                   mp->d_xiz_outer_core.cuda,
                                                                   mp->d_etax_outer_core.cuda,
                                                                   mp->d_etay_outer_core.cuda,
                                                                   mp->d_etaz_outer_core.cuda,
                                                                   mp->d_gammax_outer_core.cuda,
                                                                   mp->d_gammay_outer_core.cuda,
                                                                   mp->d_gammaz_outer_core.cuda,
                                                                   mp->d_accel_outer_core.cuda,
                                                                   mp->d_b_displ_outer_core.cuda,
                                                                   mp->d_b_accel_outer_core.cuda,
                                                                   mp->d_rho_kl_outer_core.cuda,
                                                                   mp->d_alpha_kl_outer_core.cuda,
                                                                   deltat,
                                                                   mp->NSPEC_OUTER_CORE);
  }
#endif

  GPU_ERROR_CHECKING ("compute_kernels_oc_kernel");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_kernels_strength_noise_gpu,
               COMPUTE_KERNELS_STRENGTH_NOISE_GPU) (long *Mesh_pointer_f,
                                                    realw *h_noise_surface_movie,
                                                    realw *deltat_f) {

  TRACE ("compute_kernels_strength_noise_gpu");
  // debug
  DEBUG_BACKWARD_KERNEL();

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  realw deltat = *deltat_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->nspec2D_top_crust_mantle, &num_blocks_x, &num_blocks_y);

  // copies surface buffer to GPU
  gpuCopy_todevice_realw (&mp->d_noise_surface_movie, h_noise_surface_movie, NDIM * NGLL2 * mp->nspec2D_top_crust_mantle);

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    // dimensions
    local_work_size[0] = NGLL2;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * NGLL2;
    global_work_size[1] = num_blocks_y;

    // calculates noise strength kernel
    clCheck (clSetKernelArg (mocl.kernels.compute_strength_noise_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_displ_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strength_noise_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibelm_top_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strength_noise_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strength_noise_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_noise_surface_movie.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strength_noise_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_normal_x_noise.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strength_noise_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_normal_y_noise.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strength_noise_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_normal_z_noise.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strength_noise_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_Sigma_kl.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strength_noise_kernel, idx++, sizeof (realw), (void *) &deltat));
    clCheck (clSetKernelArg (mocl.kernels.compute_strength_noise_kernel, idx++, sizeof (int), (void *) &mp->nspec2D_top_crust_mantle));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_strength_noise_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(NGLL2,1,1);

    // calculates noise strength kernel
    compute_strength_noise_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_displ_crust_mantle.cuda,
                                                                         mp->d_ibelm_top_crust_mantle.cuda,
                                                                         mp->d_ibool_crust_mantle.cuda,
                                                                         mp->d_noise_surface_movie.cuda,
                                                                         mp->d_normal_x_noise.cuda,
                                                                         mp->d_normal_y_noise.cuda,
                                                                         mp->d_normal_z_noise.cuda,
                                                                         mp->d_Sigma_kl.cuda,
                                                                         deltat,
                                                                         mp->nspec2D_top_crust_mantle);
  }
#endif

  GPU_ERROR_CHECKING ("compute_strength_noise_kernel_kernel");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_kernels_hess_gpu,
               COMPUTE_KERNELS_HESS_GPU) (long *Mesh_pointer_f,
                                          realw *deltat_f) {
  TRACE ("compute_hess_kernel_gpu");
  // debug
  DEBUG_BACKWARD_KERNEL();

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // checks
  if (! mp->approximate_hess_kl) {
    exit_on_gpu_error ("approximate_hess_kl flag not properly initialized");
  }

  int blocksize = NGLL3;   // NGLLX*NGLLY*NGLLZ
  realw deltat = *deltat_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->NSPEC_CRUST_MANTLE, &num_blocks_x, &num_blocks_y);

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    // dimensions
    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    // approximate hessian
    clCheck (clSetKernelArg (mocl.kernels.compute_hess_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_hess_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_accel_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_hess_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_accel_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_hess_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_hess_kl_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_hess_kernel, idx++, sizeof (realw), (void *) &deltat));
    clCheck (clSetKernelArg (mocl.kernels.compute_hess_kernel, idx++, sizeof (int), (void *) &mp->NSPEC_CRUST_MANTLE));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_hess_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    compute_hess_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_ibool_crust_mantle.cuda,
                                                               mp->d_accel_crust_mantle.cuda,
                                                               mp->d_b_accel_crust_mantle.cuda,
                                                               mp->d_hess_kl_crust_mantle.cuda,
                                                               deltat,
                                                               mp->NSPEC_CRUST_MANTLE);
  }
#endif

  GPU_ERROR_CHECKING ("compute_hess_kernel_gpu");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (resort_array,
               RESORT_ARRAY) (long *Mesh_pointer_f) {

  TRACE ("resort d_cijkl_kl_crust_mantle array");
  // debug

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  if(!(mp->anisotropic_kl))
    return;

  int blocksize = NGLL3;

  // blocks
  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->NSPEC_CRUST_MANTLE, &num_blocks_x, &num_blocks_y);
#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    // dimensions
    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clSetKernelArg (mocl.kernels.resort_array, idx++, sizeof (cl_mem), (void *) &mp->d_cijkl_kl_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.resort_array, idx++, sizeof (int), (void *) &mp->NSPEC_CRUST_MANTLE));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.resort_array, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
      dim3 grid(num_blocks_x,num_blocks_y);
      dim3 threads(blocksize,1,1);
      resort_array<<<grid,threads,0,mp->compute_stream>>>(mp->d_cijkl_kl_crust_mantle.cuda, mp->NSPEC_CRUST_MANTLE);
  }
#endif

  GPU_ERROR_CHECKING ("resort d_cijkl_kl_crust_mantle array");
}
