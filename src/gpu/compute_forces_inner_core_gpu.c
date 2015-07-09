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

#ifdef USE_CUDA
#ifdef USE_TEXTURES_FIELDS
//forward
realw_texture d_displ_ic_tex;
realw_texture d_accel_ic_tex;
//backward/reconstructed
realw_texture d_b_displ_ic_tex;
realw_texture d_b_accel_ic_tex;
#endif
#endif

/* ----------------------------------------------------------------------------------------------- */

void inner_core (int nb_blocks_to_compute, Mesh *mp,
                 int iphase,
                 gpu_int_mem d_ibool,
                 gpu_int_mem d_idoubling,
                 gpu_realw_mem d_xix,
                 gpu_realw_mem d_xiy,
                 gpu_realw_mem d_xiz,
                 gpu_realw_mem d_etax,
                 gpu_realw_mem d_etay,
                 gpu_realw_mem d_etaz,
                 gpu_realw_mem d_gammax,
                 gpu_realw_mem d_gammay,
                 gpu_realw_mem d_gammaz,
                 gpu_realw_mem d_kappav,
                 gpu_realw_mem d_muv,
                 gpu_realw_mem d_epsilondev_xx,
                 gpu_realw_mem d_epsilondev_yy,
                 gpu_realw_mem d_epsilondev_xy,
                 gpu_realw_mem d_epsilondev_xz,
                 gpu_realw_mem d_epsilondev_yz,
                 gpu_realw_mem d_epsilon_trace_over_3,
                 gpu_realw_mem d_one_minus_sum_beta,
                 gpu_realw_mem d_factor_common,
                 gpu_realw_mem d_R_xx,
                 gpu_realw_mem d_R_yy,
                 gpu_realw_mem d_R_xy,
                 gpu_realw_mem d_R_xz,
                 gpu_realw_mem d_R_yz,
                 gpu_realw_mem d_c11store,
                 gpu_realw_mem d_c12store,
                 gpu_realw_mem d_c13store,
                 gpu_realw_mem d_c33store,
                 gpu_realw_mem d_c44store,
                 gpu_realw_mem d_b_epsilondev_xx,
                 gpu_realw_mem d_b_epsilondev_yy,
                 gpu_realw_mem d_b_epsilondev_xy,
                 gpu_realw_mem d_b_epsilondev_xz,
                 gpu_realw_mem d_b_epsilondev_yz,
                 gpu_realw_mem d_b_epsilon_trace_over_3,
                 gpu_realw_mem d_b_R_xx,
                 gpu_realw_mem d_b_R_yy,
                 gpu_realw_mem d_b_R_xy,
                 gpu_realw_mem d_b_R_xz,
                 gpu_realw_mem d_b_R_yz,
                 int FORWARD_OR_ADJOINT) {

  GPU_ERROR_CHECKING ("before kernel inner_core");

  // safety check
  if (FORWARD_OR_ADJOINT != 1 && FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in inner_core() routine");
  }

  // if the grid can handle the number of blocks, we let it be 1D
  // grid_2_x = nb_elem_color;
  // nb_elem_color is just how many blocks we are computing now

  int blocksize = NGLL3_PADDED;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (nb_blocks_to_compute, &num_blocks_x, &num_blocks_y);


#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_kernel *inner_core_kernel_p;
    cl_uint idx = 0;

    if (FORWARD_OR_ADJOINT == 1) {
      inner_core_kernel_p = &mocl.kernels.inner_core_impl_kernel_forward;
    } else {
      // adjoint/kernel simulations
      DEBUG_BACKWARD_FORCES ();
      inner_core_kernel_p = &mocl.kernels.inner_core_impl_kernel_adjoint;
    }

    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (int), (void *) &nb_blocks_to_compute));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_ibool.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_idoubling.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_phase_ispec_inner_inner_core.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (int), (void *) &mp->num_phase_ispec_inner_core));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (int), (void *) &iphase));
    if (FORWARD_OR_ADJOINT == 1) {
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (realw), (void *) &mp->deltat));
    } else {
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (realw), (void *) &mp->b_deltat));
    }

    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (int), (void *) &mp->use_mesh_coloring_gpu));

    if (FORWARD_OR_ADJOINT == 1) {
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_displ_inner_core.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_accel_inner_core.ocl));

    } else {
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_b_displ_inner_core.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_b_accel_inner_core.ocl));
    }

    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_xix.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_xiy.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_xiz.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_etax.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_etay.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_etaz.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_gammax.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_gammay.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_gammaz.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_hprime_xx.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_hprimewgll_xx.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_wgllwgll_xy.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_wgllwgll_xz.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_wgllwgll_yz.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_kappav.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_muv.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (int), (void *) &mp->compute_and_store_strain));

    if (FORWARD_OR_ADJOINT == 1) {
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_epsilondev_xx.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_epsilondev_yy.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_epsilondev_xy.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_epsilondev_xz.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_epsilondev_yz.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_epsilon_trace_over_3.ocl));

    } else {

      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_b_epsilondev_xx.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_b_epsilondev_yy.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_b_epsilondev_xy.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_b_epsilondev_xz.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_b_epsilondev_yz.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_b_epsilon_trace_over_3.ocl));
    }

    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (int), (void *) &mp->attenuation));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (int), (void *) &mp->partial_phys_dispersion_only));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (int), (void *) &mp->use_3d_attenuation_arrays));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_one_minus_sum_beta.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_factor_common.ocl));

    if (FORWARD_OR_ADJOINT == 1) {
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_R_xx.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_R_yy.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_R_xy.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_R_xz.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_R_yz.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_alphaval.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_betaval.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_gammaval.ocl));
    } else {
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_b_R_xx.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_b_R_yy.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_b_R_xy.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_b_R_xz.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_b_R_yz.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_b_alphaval.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_b_betaval.ocl));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_b_gammaval.ocl));
    }

    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (int), (void *) &mp->anisotropic_inner_core));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_c11store.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_c12store.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_c13store.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_c33store.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_c44store.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (int), (void *) &mp->gravity));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_xstore_inner_core.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_ystore_inner_core.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_zstore_inner_core.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_minus_gravity_table.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_minus_deriv_gravity_table.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_density_table.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_wgll_cube.ocl));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (int), (void *) &mp->NSPEC_INNER_CORE_STRAIN_ONLY));
    clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (int), (void *) &mp->NSPEC_INNER_CORE));
#ifdef USE_TEXTURES_FIELDS
    if (FORWARD_OR_ADJOINT == 1) {
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_displ_ic_tex));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_accel_ic_tex));
    } else {
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_b_displ_ic_tex));
      clCheck (clSetKernelArg (*inner_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_b_accel_ic_tex));
    }
#endif
    local_work_size[0] = blocksize / GPU_ELEM_PER_THREAD;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize / GPU_ELEM_PER_THREAD;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, *inner_core_kernel_p, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize / GPU_ELEM_PER_THREAD,1,1);

    if (FORWARD_OR_ADJOINT == 1) {
      inner_core_impl_kernel_forward<<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                            d_ibool.cuda,
                                                                            d_idoubling.cuda,
                                                                            mp->d_phase_ispec_inner_inner_core.cuda,
                                                                            mp->num_phase_ispec_inner_core,
                                                                            iphase,
                                                                            mp->deltat,
                                                                            mp->use_mesh_coloring_gpu,
                                                                            mp->d_displ_inner_core.cuda,
                                                                            mp->d_accel_inner_core.cuda,
                                                                            d_xix.cuda, d_xiy.cuda, d_xiz.cuda,
                                                                            d_etax.cuda, d_etay.cuda, d_etaz.cuda,
                                                                            d_gammax.cuda, d_gammay.cuda, d_gammaz.cuda,
                                                                            mp->d_hprime_xx.cuda,
                                                                            mp->d_hprimewgll_xx.cuda,
                                                                            mp->d_wgllwgll_xy.cuda, mp->d_wgllwgll_xz.cuda, mp->d_wgllwgll_yz.cuda,
                                                                            d_kappav.cuda, d_muv.cuda,
                                                                            mp->compute_and_store_strain,
                                                                            d_epsilondev_xx.cuda,
                                                                            d_epsilondev_yy.cuda,
                                                                            d_epsilondev_xy.cuda,
                                                                            d_epsilondev_xz.cuda,
                                                                            d_epsilondev_yz.cuda,
                                                                            d_epsilon_trace_over_3.cuda,
                                                                            mp->attenuation,
                                                                            mp->partial_phys_dispersion_only,
                                                                            mp->use_3d_attenuation_arrays,
                                                                            d_one_minus_sum_beta.cuda,
                                                                            d_factor_common.cuda,
                                                                            d_R_xx.cuda,d_R_yy.cuda,d_R_xy.cuda,d_R_xz.cuda,d_R_yz.cuda,
                                                                            mp->d_alphaval.cuda,mp->d_betaval.cuda,mp->d_gammaval.cuda,
                                                                            mp->anisotropic_inner_core,
                                                                            d_c11store.cuda,d_c12store.cuda,d_c13store.cuda,
                                                                            d_c33store.cuda,d_c44store.cuda,
                                                                            mp->gravity,
                                                                            mp->d_xstore_inner_core.cuda,mp->d_ystore_inner_core.cuda,mp->d_zstore_inner_core.cuda,
                                                                            mp->d_minus_gravity_table.cuda,
                                                                            mp->d_minus_deriv_gravity_table.cuda,
                                                                            mp->d_density_table.cuda,
                                                                            mp->d_wgll_cube.cuda,
                                                                            mp->NSPEC_INNER_CORE_STRAIN_ONLY,
                                                                            mp->NSPEC_INNER_CORE);
    } else {
      // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
      // debug
      DEBUG_BACKWARD_FORCES();
      inner_core_impl_kernel_adjoint<<< grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                                             d_ibool.cuda,
                                                                             d_idoubling.cuda,
                                                                             mp->d_phase_ispec_inner_inner_core.cuda,
                                                                             mp->num_phase_ispec_inner_core,
                                                                             iphase,
                                                                             mp->b_deltat,
                                                                             mp->use_mesh_coloring_gpu,
                                                                             mp->d_b_displ_inner_core.cuda,
                                                                             mp->d_b_accel_inner_core.cuda,
                                                                             d_xix.cuda, d_xiy.cuda, d_xiz.cuda,
                                                                             d_etax.cuda, d_etay.cuda, d_etaz.cuda,
                                                                             d_gammax.cuda, d_gammay.cuda, d_gammaz.cuda,
                                                                             mp->d_hprime_xx.cuda,
                                                                             mp->d_hprimewgll_xx.cuda,
                                                                             mp->d_wgllwgll_xy.cuda, mp->d_wgllwgll_xz.cuda, mp->d_wgllwgll_yz.cuda,
                                                                             d_kappav.cuda, d_muv.cuda,
                                                                             mp->compute_and_store_strain,
                                                                             d_b_epsilondev_xx.cuda,
                                                                             d_b_epsilondev_yy.cuda,
                                                                             d_b_epsilondev_xy.cuda,
                                                                             d_b_epsilondev_xz.cuda,
                                                                             d_b_epsilondev_yz.cuda,
                                                                             d_b_epsilon_trace_over_3.cuda,
                                                                             mp->attenuation,
                                                                             mp->partial_phys_dispersion_only,
                                                                             mp->use_3d_attenuation_arrays,
                                                                             d_one_minus_sum_beta.cuda,
                                                                             d_factor_common.cuda,
                                                                             d_b_R_xx.cuda,d_b_R_yy.cuda,d_b_R_xy.cuda,d_b_R_xz.cuda,d_b_R_yz.cuda,
                                                                             mp->d_b_alphaval.cuda,mp->d_b_betaval.cuda,mp->d_b_gammaval.cuda,
                                                                             mp->anisotropic_inner_core,
                                                                             d_c11store.cuda,d_c12store.cuda,d_c13store.cuda,
                                                                             d_c33store.cuda,d_c44store.cuda,
                                                                             mp->gravity,
                                                                             mp->d_xstore_inner_core.cuda,mp->d_ystore_inner_core.cuda,mp->d_zstore_inner_core.cuda,
                                                                             mp->d_minus_gravity_table.cuda,
                                                                             mp->d_minus_deriv_gravity_table.cuda,
                                                                             mp->d_density_table.cuda,
                                                                             mp->d_wgll_cube.cuda,
                                                                             mp->NSPEC_INNER_CORE_STRAIN_ONLY,
                                                                             mp->NSPEC_INNER_CORE);
    }
  }
#endif

  GPU_ERROR_CHECKING ("inner_core");
}

/*----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_ (compute_forces_inner_core_gpu,
               COMPUTE_FORCES_INNER_CORE_GPU) (long *Mesh_pointer_f,
                                               int *iphase,
                                               int *FORWARD_OR_ADJOINT_f) {

  TRACE ("compute_forces_inner_core_gpu");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;
  int FORWARD_OR_ADJOINT = *FORWARD_OR_ADJOINT_f;

  // safety check
  if (FORWARD_OR_ADJOINT != 1 && FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in compute_forces_inner_core_gpu() routine");
  }

  // determines number of elements to loop over (inner/outer elements
  int num_elements;
  if (*iphase == 1) {
    num_elements = mp->nspec_outer_inner_core;
  } else {
    num_elements = mp->nspec_inner_inner_core;
  }
  // checks if anything to do
  if (num_elements == 0) return;

  // mesh coloring
  if (mp->use_mesh_coloring_gpu) {

    // note: array offsets require sorted arrays, such that e.g. ibool starts with elastic elements
    //         and followed by acoustic ones.
    //         elastic elements also start with outer than inner element ordering

    int nb_colors, nb_blocks_to_compute;
    int istart;
    int offset, offset_nonpadded;
    int offset_nonpadded_att1, offset_nonpadded_att2, offset_nonpadded_att3;
    int offset_nonpadded_strain;
    int offset_ispec;

    // sets up color loop
    if (mp->NSPEC_INNER_CORE > COLORING_MIN_NSPEC_INNER_CORE) {
      if (*iphase == 1) {
        // outer elements
        nb_colors = mp->num_colors_outer_inner_core;
        istart = 0;

        // array offsets
        offset = 0;
        offset_nonpadded = 0;
        offset_nonpadded_att1 = 0;
        offset_nonpadded_att2 = 0;
        offset_nonpadded_att3 = 0;
        offset_nonpadded_strain = 0;
        offset_ispec = 0;
      } else {
        // inner elements (start after outer elements)
        nb_colors = mp->num_colors_outer_inner_core + mp->num_colors_inner_inner_core;
        istart = mp->num_colors_outer_inner_core;

        // array offsets
        offset = mp->nspec_outer_inner_core * NGLL3_PADDED;
        offset_nonpadded = mp->nspec_outer_inner_core * NGLL3;
        offset_nonpadded_att1 = mp->nspec_outer_inner_core * NGLL3 * N_SLS;
        // for factor_common array
        if (mp->use_3d_attenuation_arrays) {
          offset_nonpadded_att2 = mp->nspec_outer_inner_core * NGLL3;
          offset_nonpadded_att3 = mp->nspec_outer_inner_core * NGLL3 * N_SLS;
        } else {
          offset_nonpadded_att2 = mp->nspec_outer_inner_core;
          offset_nonpadded_att3 = mp->nspec_outer_inner_core * N_SLS;
        }
        // for idoubling array
        offset_ispec = mp->nspec_outer_inner_core;
        // for strain
        if (! mp->NSPEC_INNER_CORE_STRAIN_ONLY == 1) {
          offset_nonpadded_strain = (mp->nspec_outer_inner_core) * NGLL3;
        }
      }
    } else {

      // poor element count, only use 1 color per inner/outer run

      if (*iphase == 1) {
        // outer elements
        nb_colors = 1;
        istart = 0;

        // array offsets
        offset = 0;
        offset_nonpadded = 0;
        offset_nonpadded_att1 = 0;
        offset_nonpadded_att2 = 0;
        offset_nonpadded_att3 = 0;
        offset_nonpadded_strain = 0;
        offset_ispec = 0;
      } else {
        // inner element colors (start after outer elements)
        nb_colors = 1;
        istart = 0;

        // array offsets
        offset = mp->nspec_outer_inner_core * NGLL3_PADDED;
        offset_nonpadded = mp->nspec_outer_inner_core * NGLL3;
        offset_nonpadded_att1 = mp->nspec_outer_inner_core * NGLL3 * N_SLS;
        // for factor_common array
        if (mp->use_3d_attenuation_arrays) {
          offset_nonpadded_att2 = mp->nspec_outer_inner_core * NGLL3;
          offset_nonpadded_att3 = mp->nspec_outer_inner_core * NGLL3 * N_SLS;
        } else {
          offset_nonpadded_att2 = mp->nspec_outer_inner_core;
          offset_nonpadded_att3 = mp->nspec_outer_inner_core * N_SLS;
        }
        // for idoubling array
        offset_ispec = mp->nspec_outer_inner_core;
        // for strain
        if (mp->NSPEC_INNER_CORE_STRAIN_ONLY != 1) {
          offset_nonpadded_strain = mp->nspec_outer_inner_core * NGLL3;
        }
      }
    }


    // loops over colors
    int icolor;
    for (icolor = istart; icolor < nb_colors; icolor++) {

      // gets number of elements for this color
      if (mp->NSPEC_INNER_CORE > COLORING_MIN_NSPEC_INNER_CORE) {
        nb_blocks_to_compute = mp->h_num_elem_colors_inner_core[icolor];
      } else {
        nb_blocks_to_compute = num_elements;
      }

#if ENABLE_VERY_SLOW_ERROR_CHECKING == 1
      // checks
      if (nb_blocks_to_compute <= 0) {
        printf ("Error number of color blocks in inner_core: %d -- color = %d \n",
                nb_blocks_to_compute, icolor);
        exit (EXIT_FAILURE);
      }
#endif
      INITIALIZE_OFFSET();

      INIT_OFFSET(d_ibool_inner_core, offset_nonpadded);
      INIT_OFFSET(d_idoubling_inner_core, offset_ispec);
      INIT_OFFSET(d_xix_inner_core, offset);
      INIT_OFFSET(d_xiy_inner_core, offset);
      INIT_OFFSET(d_xiz_inner_core, offset);
      INIT_OFFSET(d_etax_inner_core, offset);
      INIT_OFFSET(d_etay_inner_core, offset);
      INIT_OFFSET(d_etaz_inner_core, offset);
      INIT_OFFSET(d_gammax_inner_core, offset);
      INIT_OFFSET(d_gammay_inner_core, offset);
      INIT_OFFSET(d_gammaz_inner_core, offset);
      INIT_OFFSET(d_kappavstore_inner_core, offset);
      INIT_OFFSET(d_muvstore_inner_core, offset);
      INIT_OFFSET(d_epsilondev_xx_inner_core, offset_nonpadded);
      INIT_OFFSET(d_epsilondev_yy_inner_core, offset_nonpadded);
      INIT_OFFSET(d_epsilondev_xy_inner_core, offset_nonpadded);
      INIT_OFFSET(d_epsilondev_xz_inner_core, offset_nonpadded);
      INIT_OFFSET(d_epsilondev_yz_inner_core, offset_nonpadded);
      INIT_OFFSET(d_eps_trace_over_3_inner_core, offset_nonpadded_strain);
      INIT_OFFSET(d_one_minus_sum_beta_inner_core, offset_nonpadded_att2);
      INIT_OFFSET(d_factor_common_inner_core, offset_nonpadded_att3);
      INIT_OFFSET(d_R_xx_inner_core, offset_nonpadded_att1);
      INIT_OFFSET(d_R_yy_inner_core, offset_nonpadded_att1);
      INIT_OFFSET(d_R_xy_inner_core, offset_nonpadded_att1);
      INIT_OFFSET(d_R_xz_inner_core, offset_nonpadded_att1);
      INIT_OFFSET(d_R_yz_inner_core, offset_nonpadded_att1);
      INIT_OFFSET(d_c11store_inner_core, offset);
      INIT_OFFSET(d_c12store_inner_core, offset);
      INIT_OFFSET(d_c13store_inner_core, offset);
      INIT_OFFSET(d_c33store_inner_core, offset);
      INIT_OFFSET(d_c44store_inner_core, offset);
      INIT_OFFSET(d_b_epsilondev_xx_inner_core, offset_nonpadded);
      INIT_OFFSET(d_b_epsilondev_yy_inner_core, offset_nonpadded);
      INIT_OFFSET(d_b_epsilondev_xy_inner_core, offset_nonpadded);
      INIT_OFFSET(d_b_epsilondev_xz_inner_core, offset_nonpadded);
      INIT_OFFSET(d_b_epsilondev_yz_inner_core, offset_nonpadded);
      INIT_OFFSET(d_b_eps_trace_over_3_inner_core, offset_nonpadded);
      INIT_OFFSET(d_b_R_xx_inner_core, offset_nonpadded_att1);
      INIT_OFFSET(d_b_R_yy_inner_core, offset_nonpadded_att1);
      INIT_OFFSET(d_b_R_xy_inner_core, offset_nonpadded_att1);
      INIT_OFFSET(d_b_R_xz_inner_core, offset_nonpadded_att1);
      INIT_OFFSET(d_b_R_yz_inner_core, offset_nonpadded_att1);

      inner_core (nb_blocks_to_compute, mp,
                  *iphase,
                  PASS_OFFSET(d_ibool_inner_core, offset_nonpadded),
                  PASS_OFFSET(d_idoubling_inner_core, offset_ispec),
                  PASS_OFFSET(d_xix_inner_core, offset),
                  PASS_OFFSET(d_xiy_inner_core, offset),
                  PASS_OFFSET(d_xiz_inner_core, offset),
                  PASS_OFFSET(d_etax_inner_core, offset),
                  PASS_OFFSET(d_etay_inner_core, offset),
                  PASS_OFFSET(d_etaz_inner_core, offset),
                  PASS_OFFSET(d_gammax_inner_core, offset),
                  PASS_OFFSET(d_gammay_inner_core, offset),
                  PASS_OFFSET(d_gammaz_inner_core, offset),
                  PASS_OFFSET(d_kappavstore_inner_core, offset),
                  PASS_OFFSET(d_muvstore_inner_core, offset),
                  PASS_OFFSET(d_epsilondev_xx_inner_core, offset_nonpadded),
                  PASS_OFFSET(d_epsilondev_yy_inner_core, offset_nonpadded),
                  PASS_OFFSET(d_epsilondev_xy_inner_core, offset_nonpadded),
                  PASS_OFFSET(d_epsilondev_xz_inner_core, offset_nonpadded),
                  PASS_OFFSET(d_epsilondev_yz_inner_core, offset_nonpadded),
                  PASS_OFFSET(d_eps_trace_over_3_inner_core, offset_nonpadded_strain),
                  PASS_OFFSET(d_one_minus_sum_beta_inner_core, offset_nonpadded_att2),
                  PASS_OFFSET(d_factor_common_inner_core, offset_nonpadded_att3),
                  PASS_OFFSET(d_R_xx_inner_core, offset_nonpadded_att1),
                  PASS_OFFSET(d_R_yy_inner_core, offset_nonpadded_att1),
                  PASS_OFFSET(d_R_xy_inner_core, offset_nonpadded_att1),
                  PASS_OFFSET(d_R_xz_inner_core, offset_nonpadded_att1),
                  PASS_OFFSET(d_R_yz_inner_core, offset_nonpadded_att1),
                  PASS_OFFSET(d_c11store_inner_core, offset),
                  PASS_OFFSET(d_c12store_inner_core, offset),
                  PASS_OFFSET(d_c13store_inner_core, offset),
                  PASS_OFFSET(d_c33store_inner_core, offset),
                  PASS_OFFSET(d_c44store_inner_core, offset),
                  PASS_OFFSET(d_b_epsilondev_xx_inner_core, offset_nonpadded),
                  PASS_OFFSET(d_b_epsilondev_yy_inner_core, offset_nonpadded),
                  PASS_OFFSET(d_b_epsilondev_xy_inner_core, offset_nonpadded),
                  PASS_OFFSET(d_b_epsilondev_xz_inner_core, offset_nonpadded),
                  PASS_OFFSET(d_b_epsilondev_yz_inner_core, offset_nonpadded),
                  PASS_OFFSET(d_b_eps_trace_over_3_inner_core, offset_nonpadded),
                  PASS_OFFSET(d_b_R_xx_inner_core, offset_nonpadded_att1),
                  PASS_OFFSET(d_b_R_yy_inner_core, offset_nonpadded_att1),
                  PASS_OFFSET(d_b_R_xy_inner_core, offset_nonpadded_att1),
                  PASS_OFFSET(d_b_R_xz_inner_core, offset_nonpadded_att1),
                  PASS_OFFSET(d_b_R_yz_inner_core, offset_nonpadded_att1),
                  FORWARD_OR_ADJOINT);

      RELEASE_OFFSET(d_ibool_inner_core, offset_nonpadded);
      RELEASE_OFFSET(d_idoubling_inner_core, offset_ispec);
      RELEASE_OFFSET(d_xix_inner_core, offset);
      RELEASE_OFFSET(d_xiy_inner_core, offset);
      RELEASE_OFFSET(d_xiz_inner_core, offset);
      RELEASE_OFFSET(d_etax_inner_core, offset);
      RELEASE_OFFSET(d_etay_inner_core, offset);
      RELEASE_OFFSET(d_etaz_inner_core, offset);
      RELEASE_OFFSET(d_gammax_inner_core, offset);
      RELEASE_OFFSET(d_gammay_inner_core, offset);
      RELEASE_OFFSET(d_gammaz_inner_core, offset);
      RELEASE_OFFSET(d_kappavstore_inner_core, offset);
      RELEASE_OFFSET(d_muvstore_inner_core, offset);
      RELEASE_OFFSET(d_epsilondev_xx_inner_core, offset_nonpadded);
      RELEASE_OFFSET(d_epsilondev_yy_inner_core, offset_nonpadded);
      RELEASE_OFFSET(d_epsilondev_xy_inner_core, offset_nonpadded);
      RELEASE_OFFSET(d_epsilondev_xz_inner_core, offset_nonpadded);
      RELEASE_OFFSET(d_epsilondev_yz_inner_core, offset_nonpadded);
      RELEASE_OFFSET(d_eps_trace_over_3_inner_core, offset_nonpadded_strain);
      RELEASE_OFFSET(d_one_minus_sum_beta_inner_core, offset_nonpadded_att2);
      RELEASE_OFFSET(d_factor_common_inner_core, offset_nonpadded_att3);
      RELEASE_OFFSET(d_R_xx_inner_core, offset_nonpadded_att1);
      RELEASE_OFFSET(d_R_yy_inner_core, offset_nonpadded_att1);
      RELEASE_OFFSET(d_R_xy_inner_core, offset_nonpadded_att1);
      RELEASE_OFFSET(d_R_xz_inner_core, offset_nonpadded_att1);
      RELEASE_OFFSET(d_R_yz_inner_core, offset_nonpadded_att1);
      RELEASE_OFFSET(d_c11store_inner_core, offset);
      RELEASE_OFFSET(d_c12store_inner_core, offset);
      RELEASE_OFFSET(d_c13store_inner_core, offset);
      RELEASE_OFFSET(d_c33store_inner_core, offset);
      RELEASE_OFFSET(d_c44store_inner_core, offset);
      RELEASE_OFFSET(d_b_epsilondev_xx_inner_core, offset_nonpadded);
      RELEASE_OFFSET(d_b_epsilondev_yy_inner_core, offset_nonpadded);
      RELEASE_OFFSET(d_b_epsilondev_xy_inner_core, offset_nonpadded);
      RELEASE_OFFSET(d_b_epsilondev_xz_inner_core, offset_nonpadded);
      RELEASE_OFFSET(d_b_epsilondev_yz_inner_core, offset_nonpadded);
      RELEASE_OFFSET(d_b_eps_trace_over_3_inner_core, offset_nonpadded);
      RELEASE_OFFSET(d_b_R_xx_inner_core, offset_nonpadded_att1);
      RELEASE_OFFSET(d_b_R_yy_inner_core, offset_nonpadded_att1);
      RELEASE_OFFSET(d_b_R_xy_inner_core, offset_nonpadded_att1);
      RELEASE_OFFSET(d_b_R_xz_inner_core, offset_nonpadded_att1);
      RELEASE_OFFSET(d_b_R_yz_inner_core, offset_nonpadded_att1);

      // for padded and aligned arrays
      offset += nb_blocks_to_compute * NGLL3_PADDED;
      // for no-aligned arrays
      offset_nonpadded += nb_blocks_to_compute * NGLL3;
      offset_nonpadded_att1 += nb_blocks_to_compute * NGLL3 * N_SLS;
      // for factor_common array
      if (mp->use_3d_attenuation_arrays) {
        offset_nonpadded_att2 += nb_blocks_to_compute * NGLL3;
        offset_nonpadded_att3 += nb_blocks_to_compute * NGLL3 * N_SLS;
      } else {
        offset_nonpadded_att2 += nb_blocks_to_compute * 1;
        offset_nonpadded_att3 += nb_blocks_to_compute * 1 * N_SLS;
      }
      // for array (ispec)
      offset_ispec += nb_blocks_to_compute;
      // for strain
      if (mp->NSPEC_INNER_CORE_STRAIN_ONLY != 1) {
        offset_nonpadded_strain += nb_blocks_to_compute * NGLL3;
      }
    }   // icolor

  } else {

    // no mesh coloring: uses atomic updates
    inner_core (num_elements, mp,
                *iphase,
                mp->d_ibool_inner_core,
                mp->d_idoubling_inner_core,
                mp->d_xix_inner_core,
                mp->d_xiy_inner_core,
                mp->d_xiz_inner_core,
                mp->d_etax_inner_core,
                mp->d_etay_inner_core,
                mp->d_etaz_inner_core,
                mp->d_gammax_inner_core,
                mp->d_gammay_inner_core,
                mp->d_gammaz_inner_core,
                mp->d_kappavstore_inner_core,
                mp->d_muvstore_inner_core,
                mp->d_epsilondev_xx_inner_core,
                mp->d_epsilondev_yy_inner_core,
                mp->d_epsilondev_xy_inner_core,
                mp->d_epsilondev_xz_inner_core,
                mp->d_epsilondev_yz_inner_core,
                mp->d_eps_trace_over_3_inner_core,
                mp->d_one_minus_sum_beta_inner_core,
                mp->d_factor_common_inner_core,
                mp->d_R_xx_inner_core,
                mp->d_R_yy_inner_core,
                mp->d_R_xy_inner_core,
                mp->d_R_xz_inner_core,
                mp->d_R_yz_inner_core,
                mp->d_c11store_inner_core,
                mp->d_c12store_inner_core,
                mp->d_c13store_inner_core,
                mp->d_c33store_inner_core,
                mp->d_c44store_inner_core,
                mp->d_b_epsilondev_xx_inner_core,
                mp->d_b_epsilondev_yy_inner_core,
                mp->d_b_epsilondev_xy_inner_core,
                mp->d_b_epsilondev_xz_inner_core,
                mp->d_b_epsilondev_yz_inner_core,
                mp->d_b_eps_trace_over_3_inner_core,
                mp->d_b_R_xx_inner_core,
                mp->d_b_R_yy_inner_core,
                mp->d_b_R_xy_inner_core,
                mp->d_b_R_xz_inner_core,
                mp->d_b_R_yz_inner_core,
                FORWARD_OR_ADJOINT);
  }

  GPU_ERROR_CHECKING ("compute_forces_inner_core_ocl");
}
