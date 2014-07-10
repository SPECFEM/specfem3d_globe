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

extern EXTERN_LANG
void FC_FUNC_ (compute_kernels_cm_gpu,
               COMPUTE_KERNELS_CM_GPU) (long *Mesh_pointer, realw *deltat_f) {

  TRACE ("compute_cm_gpu");
  // debug
  DEBUG_BACKWARD_KERNEL ();

  Mesh *mp = (Mesh *) *Mesh_pointer;   //get mesh pointer out of Fortran integer container

  int blocksize = NGLL3;
  realw deltat = *deltat_f;

  // blocks
  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->NSPEC_CRUST_MANTLE, &num_blocks_x, &num_blocks_y);

#ifdef USE_OPENCL
  size_t global_work_size[2];
  size_t local_work_size[2];
  cl_uint idx = 0;

  if (run_opencl) {
    // density kernel
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_accel_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_displ_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_rho_kl_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (int), (void *) &mp->NSPEC_CRUST_MANTLE));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (realw), (void *) &deltat));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_rho_kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  if (run_cuda) {
    // density kernel
    compute_rho_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_ibool_crust_mantle.cuda,
                                                              mp->d_accel_crust_mantle.cuda,
                                                              mp->d_b_displ_crust_mantle.cuda,
                                                              mp->d_rho_kl_crust_mantle.cuda,
                                                              mp->NSPEC_CRUST_MANTLE,
                                                              deltat);
  }
#endif
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
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xx_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_yy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_yz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_eps_trace_over_3_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_beta_kl_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_alpha_kl_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (int), (void *) &mp->NSPEC_CRUST_MANTLE));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (realw), (void *) &deltat));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_displ_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xix_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xiy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xiz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etax_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etay_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etaz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammax_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammay_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammaz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_iso_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_hprime_xx.ocl));


        local_work_size[0] = blocksize;
        local_work_size[1] = 1;
        global_work_size[0] = num_blocks_x * blocksize;
        global_work_size[1] = num_blocks_y;

        clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_iso_undo_att_kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL));

      } else {
        // anisotropic kernels
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xx_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_yy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_xz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_epsilondev_yz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_eps_trace_over_3_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_beta_kl_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_alpha_kl_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (int), (void *) &mp->NSPEC_CRUST_MANTLE));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (realw), (void *) &deltat));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_displ_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xix_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xiy_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xiz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etax_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etay_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etaz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammax_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammay_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammaz_crust_mantle.ocl));
        clCheck (clSetKernelArg (mocl.kernels.compute_ani_undo_att_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_hprime_xx.ocl));

        local_work_size[0] = blocksize;
        local_work_size[1] = 1;
        global_work_size[0] = num_blocks_x * blocksize;
        global_work_size[1] = num_blocks_y;

        clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_ani_undo_att_kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL));
      }
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      // computes strain locally based on current backward/reconstructed (b_displ) wavefield
      if(! mp->anisotropic_kl){
        // isotropic kernels
        compute_iso_undo_att_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_epsilondev_xx_crust_mantle.cuda,
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
        compute_ani_undo_att_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_epsilondev_xx_crust_mantle.cuda,
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
    // computes strain locally based on current backward/reconstructed (b_displ) wavefield
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

        local_work_size[0] = blocksize;
        local_work_size[1] = 1;
        global_work_size[0] = num_blocks_x * blocksize;
        global_work_size[1] = num_blocks_y;

        clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_iso_kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL));

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

        local_work_size[0] = blocksize;
        local_work_size[1] = 1;
        global_work_size[0] = num_blocks_x * blocksize;
        global_work_size[1] = num_blocks_y;

        clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_ani_kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL));
      }
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      // takes strain arrays computed from previous compute_forces call
      if (! mp->anisotropic_kl){
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
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_gpu_error ("compute_cm_gpu");
#endif
}


/*----------------------------------------------------------------------------------------------- */
// inner_core
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_kernels_ic_gpu,
               COMPUTE_KERNELS_IC_GPU) (long *Mesh_pointer, realw *deltat_f) {

  TRACE("compute_kernels_ic_gpu");
  // debug
  DEBUG_BACKWARD_KERNEL();

  Mesh *mp = (Mesh *) *Mesh_pointer;   //get mesh pointer out of Fortran integer container

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
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_accel_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_displ_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_rho_kl_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (int), (void *) &mp->NSPEC_INNER_CORE));
    clCheck (clSetKernelArg (mocl.kernels.compute_rho_kernel, idx++, sizeof (realw), (void *) &deltat));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_rho_kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  if (run_cuda) {
    compute_rho_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_ibool_inner_core.cuda,
                                                              mp->d_accel_inner_core.cuda,
                                                              mp->d_b_displ_inner_core.cuda,
                                                              mp->d_rho_kl_inner_core.cuda,
                                                              mp->NSPEC_INNER_CORE,
                                                              deltat);
  }
#endif

  // checks if strain is available
  if (mp->undo_attenuation) {
    // checks strain array size
    if(mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY == 1) {
      exit_on_error("compute_kernels_cm_cuda NSPEC_CRUST_MANTLE_STRAIN_ONLY invalid with undo_att");
    }
  }

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

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

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
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_gpu_error ("compute_ic_gpu");
#endif
}

extern EXTERN_LANG
void FC_FUNC_ (compute_kernels_oc_gpu,
               COMPUTE_KERNELS_OC_GPU) (long *Mesh_pointer, realw *deltat_f) {

  TRACE ("compute_kernels_oc_gpu");

  Mesh *mp = (Mesh *) *Mesh_pointer;   //get mesh pointer out of Fortran integer container

  int blocksize = NGLL3;   // NGLLX*NGLLY*NGLLZ
  realw deltat = *deltat_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->NSPEC_OUTER_CORE, &num_blocks_x, &num_blocks_y);


#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

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

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_acoustic_kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL));
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
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_gpu_error ("compute_kernels_oc_kernel");
#endif
}

extern EXTERN_LANG
void FC_FUNC_ (compute_kernels_strgth_noise_gpu,
               COMPUTE_KERNELS_STRGTH_NOISE_GPU) (long *Mesh_pointer,
                                                  realw *h_noise_surface_movie,
                                                  realw *deltat_f) {

  TRACE ("compute_kernels_strgth_noise_gpu");

  Mesh *mp = (Mesh *) *Mesh_pointer;   //get mesh pointer out of Fortran integer container

  realw deltat = *deltat_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->nspec2D_top_crust_mantle, &num_blocks_x, &num_blocks_y);

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    // copies surface buffer to GPU
    clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_noise_surface_movie.ocl, CL_FALSE, 0,
                                   NDIM * NGLL2 * mp->nspec2D_top_crust_mantle * sizeof (realw),
                                   h_noise_surface_movie, 0, NULL, NULL));

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

    local_work_size[0] = NGLL2;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * NGLL2;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_strength_noise_kernel,
                                     2, NULL, global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(NGLL2,1,1);

    // copies surface buffer to GPU
    print_CUDA_error_if_any(cudaMemcpy(mp->d_noise_surface_movie.cuda,h_noise_surface_movie,
                                       NDIM*NGLL2*(mp->nspec2D_top_crust_mantle)*sizeof(realw),
                                       cudaMemcpyHostToDevice),90900);

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
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_gpu_error ("compute_strength_noise_kernel_kernel");
#endif
}

extern EXTERN_LANG
void FC_FUNC_ (compute_kernels_hess_gpu,
               COMPUTE_KERNELS_HESS_GPU) (long *Mesh_pointer,
                                          realw *deltat_f) {
  TRACE ("compute_hess_kernel_gpu");

  Mesh *mp = (Mesh *) *Mesh_pointer;   //get mesh pointer out of Fortran integer container

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

    clCheck (clSetKernelArg (mocl.kernels.compute_hess_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_hess_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_accel_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_hess_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_accel_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_hess_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_hess_kl_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_hess_kernel, idx++, sizeof (realw), (void *) &deltat));
    clCheck (clSetKernelArg (mocl.kernels.compute_hess_kernel, idx++, sizeof (int), (void *) &mp->NSPEC_CRUST_MANTLE));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

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
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_gpu_error ("compute_hess_kernel_gpu");
#endif
}
