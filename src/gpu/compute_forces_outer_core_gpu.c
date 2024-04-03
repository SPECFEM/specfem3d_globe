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

#ifdef USE_CUDA
#ifdef USE_TEXTURES_FIELDS
//forward
realw_texture d_displ_oc_tex;
realw_texture d_accel_oc_tex;
//backward/reconstructed
realw_texture d_b_displ_oc_tex;
realw_texture d_b_accel_oc_tex;
#endif
#endif

/* ----------------------------------------------------------------------------------------------- */

void outer_core (int nb_blocks_to_compute, Mesh *mp,
                 int iphase,
                 gpu_int_mem d_ibool,
                 gpu_realw_mem d_xix,
                 gpu_realw_mem d_xiy,
                 gpu_realw_mem d_xiz,
                 gpu_realw_mem d_etax,
                 gpu_realw_mem d_etay,
                 gpu_realw_mem d_etaz,
                 gpu_realw_mem d_gammax,
                 gpu_realw_mem d_gammay,
                 gpu_realw_mem d_gammaz,
                 realw timeval,
                 gpu_realw_mem d_A_array_rotation,
                 gpu_realw_mem d_B_array_rotation,
                 gpu_realw_mem d_b_A_array_rotation,
                 gpu_realw_mem d_b_B_array_rotation,
                 gpu_realw_mem d_A_array_rotation_lddrk,
                 gpu_realw_mem d_B_array_rotation_lddrk,
                 gpu_realw_mem d_b_A_array_rotation_lddrk,
                 gpu_realw_mem d_b_B_array_rotation_lddrk,
                 realw alpha_lddrk,realw beta_lddrk,
                 int FORWARD_OR_ADJOINT) {

  GPU_ERROR_CHECKING ("before outer_core kernel Kernel_2");

  // safety check
  if (FORWARD_OR_ADJOINT != 1 && FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in outer_core() routine");
  }

  // if the grid can handle the number of blocks, we let it be 1D
  // grid_2_x = nb_elem_color;
  // nb_elem_color is just how many blocks we are computing now
  int blocksize = NGLL3_PADDED;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (nb_blocks_to_compute, &num_blocks_x, &num_blocks_y);

  // defines local parameters for forward/adjoint function calls
  realw deltat,two_omega_earth;
  gpu_realw_mem displ,accel;
  gpu_realw_mem A_array_rotation, B_array_rotation;
  gpu_realw_mem A_array_rotation_lddrk, B_array_rotation_lddrk;

  // sets gpu arrays
  if (FORWARD_OR_ADJOINT == 1) {
    displ = mp->d_displ_outer_core;
    accel = mp->d_accel_outer_core;
    deltat = mp->deltat;
    two_omega_earth = mp->two_omega_earth;
    A_array_rotation = d_A_array_rotation;
    B_array_rotation = d_B_array_rotation;
    A_array_rotation_lddrk = d_A_array_rotation_lddrk;
    B_array_rotation_lddrk = d_B_array_rotation_lddrk;
  } else {
    // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
    displ = mp->d_b_displ_outer_core;
    accel = mp->d_b_accel_outer_core;
    deltat = mp->b_deltat;
    two_omega_earth = mp->b_two_omega_earth;
    A_array_rotation = d_b_A_array_rotation;
    B_array_rotation = d_b_B_array_rotation;
    A_array_rotation_lddrk = d_b_A_array_rotation_lddrk;
    B_array_rotation_lddrk = d_b_B_array_rotation_lddrk;
  }

  // kernel timing
  gpu_event start,stop;
  if (GPU_TIMING){ start_timing_gpu(&start,&stop); }

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_kernel *outer_core_kernel_p;
    cl_uint idx = 0;

    // sets kernel function
    if (FORWARD_OR_ADJOINT == 1) {
      outer_core_kernel_p = &mocl.kernels.outer_core_impl_kernel_forward;
    } else {
      // debug
      DEBUG_BACKWARD_FORCES ();
      // adjoint/kernel simulations
      outer_core_kernel_p = &mocl.kernels.outer_core_impl_kernel_adjoint;
    }

    // adds arguments
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (int), (void *) &nb_blocks_to_compute));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_ibool.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_phase_ispec_inner_outer_core.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (int), (void *) &mp->num_phase_ispec_outer_core));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (int), (void *) &iphase));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (int), (void *) &mp->use_mesh_coloring_gpu));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &displ.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_xix.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_xiy.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_xiz.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_etax.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_etay.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_etaz.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_gammax.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_gammay.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &d_gammaz.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_hprime_xx.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_hprimewgll_xx.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_wgllwgll_xy.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_wgllwgll_xz.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_wgllwgll_yz.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (int), (void *) &mp->gravity));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_gravity_pre_store_outer_core.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_wgll_cube.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (int), (void *) &mp->rotation));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (realw), (void *) &timeval));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (realw), (void *) &two_omega_earth));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (realw), (void *) &deltat));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &A_array_rotation.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &B_array_rotation.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &A_array_rotation_lddrk.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &B_array_rotation_lddrk.ocl));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (realw), (void *) &alpha_lddrk));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (realw), (void *) &beta_lddrk));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (int), (void *) &mp->use_lddrk));
    clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (int), (void *) &mp->NSPEC_OUTER_CORE));
#ifdef USE_TEXTURES_FIELDS
    if (FORWARD_OR_ADJOINT == 1) {
      clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_displ_oc_tex));
      clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_accel_oc_tex));
    } else {
      clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_b_displ_oc_tex));
      clCheck (clSetKernelArg (*outer_core_kernel_p, idx++, sizeof (cl_mem), (void *) &mp->d_b_accel_oc_tex));
    }
#endif
    local_work_size[0] = blocksize / GPU_ELEM_PER_THREAD;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize / GPU_ELEM_PER_THREAD;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, *outer_core_kernel_p, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize / GPU_ELEM_PER_THREAD,1,1);

    // defines function pointer to __global__ function (taken from definition in file kernel_proto.cu.h)
    // since forward and adjoint function calls are identical and only the passed arrays change
    outer_core_impl_kernel outer_core_kernel_p;

    // selects function call
    if (FORWARD_OR_ADJOINT == 1) {
      // forward wavefields -> FORWARD_OR_ADJOINT == 1
      outer_core_kernel_p = &outer_core_impl_kernel_forward;
    } else {
      // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
      DEBUG_BACKWARD_FORCES();
      outer_core_kernel_p = &outer_core_impl_kernel_adjoint;
    }

    outer_core_kernel_p<<<grid,threads,0,mp->compute_stream>>>(nb_blocks_to_compute,
                                                               d_ibool.cuda,
                                                               mp->d_phase_ispec_inner_outer_core.cuda,
                                                               mp->num_phase_ispec_outer_core,
                                                               iphase,
                                                               mp->use_mesh_coloring_gpu,
                                                               displ.cuda,
                                                               accel.cuda,
                                                               d_xix.cuda,d_xiy.cuda,d_xiz.cuda,
                                                               d_etax.cuda,d_etay.cuda,d_etaz.cuda,
                                                               d_gammax.cuda,d_gammay.cuda,d_gammaz.cuda,
                                                               mp->d_hprime_xx.cuda,
                                                               mp->d_hprimewgll_xx.cuda,
                                                               mp->d_wgllwgll_xy.cuda,mp->d_wgllwgll_xz.cuda,mp->d_wgllwgll_yz.cuda,
                                                               mp->gravity,
                                                               mp->d_gravity_pre_store_outer_core.cuda,
                                                               mp->d_wgll_cube.cuda,
                                                               mp->rotation,
                                                               timeval,
                                                               two_omega_earth,
                                                               deltat,
                                                               A_array_rotation.cuda,
                                                               B_array_rotation.cuda,
                                                               A_array_rotation_lddrk.cuda,
                                                               B_array_rotation_lddrk.cuda,
                                                               alpha_lddrk,beta_lddrk,
                                                               mp->use_lddrk,
                                                               mp->NSPEC_OUTER_CORE);
  }
#endif
#ifdef USE_HIP
  if (run_hip) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize / GPU_ELEM_PER_THREAD,1,1);

    // defines function pointer to __global__ function (taken from definition in file kernel_proto.cu.h)
    // since forward and adjoint function calls are identical and only the passed arrays change
    outer_core_impl_kernel outer_core_kernel_p;

    // selects function call
    if (FORWARD_OR_ADJOINT == 1) {
      // forward wavefields -> FORWARD_OR_ADJOINT == 1
      outer_core_kernel_p = &outer_core_impl_kernel_forward;
    } else {
      // backward/reconstructed wavefields -> FORWARD_OR_ADJOINT == 3
      DEBUG_BACKWARD_FORCES();
      outer_core_kernel_p = &outer_core_impl_kernel_adjoint;
    }

    hipLaunchKernelGGL( outer_core_kernel_p, grid, threads, 0, mp->compute_stream,
                                                               nb_blocks_to_compute,
                                                               d_ibool.hip,
                                                               mp->d_phase_ispec_inner_outer_core.hip,
                                                               mp->num_phase_ispec_outer_core,
                                                               iphase,
                                                               mp->use_mesh_coloring_gpu,
                                                               displ.hip,
                                                               accel.hip,
                                                               d_xix.hip,d_xiy.hip,d_xiz.hip,
                                                               d_etax.hip,d_etay.hip,d_etaz.hip,
                                                               d_gammax.hip,d_gammay.hip,d_gammaz.hip,
                                                               mp->d_hprime_xx.hip,
                                                               mp->d_hprimewgll_xx.hip,
                                                               mp->d_wgllwgll_xy.hip,mp->d_wgllwgll_xz.hip,mp->d_wgllwgll_yz.hip,
                                                               mp->gravity,
                                                               mp->d_gravity_pre_store_outer_core.hip,
                                                               mp->d_wgll_cube.hip,
                                                               mp->rotation,
                                                               timeval,
                                                               two_omega_earth,
                                                               deltat,
                                                               A_array_rotation.hip,
                                                               B_array_rotation.hip,
                                                               A_array_rotation_lddrk.hip,
                                                               B_array_rotation_lddrk.hip,
                                                               alpha_lddrk,beta_lddrk,
                                                               mp->use_lddrk,
                                                               mp->NSPEC_OUTER_CORE);
  }
#endif

  // kernel timing
  if (GPU_TIMING){ stop_timing_gpu(&start,&stop,"outer_core"); }

  GPU_ERROR_CHECKING ("kernel outer_core");
}

/*----------------------------------------------------------------------------------------------- */

// main compute_forces_outer_core GPU routine

/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_forces_outer_core_gpu,
               COMPUTE_FORCES_OUTER_CORE_GPU) (long *Mesh_pointer_f,
                                               int *iphase,
                                               realw *timeval_f,
                                               realw *alpha_lddrk_f,realw *beta_lddrk_f,
                                               int *FORWARD_OR_ADJOINT_f) {

  TRACE ("compute_forces_outer_core_gpu");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;
  realw timeval = *timeval_f;
  int FORWARD_OR_ADJOINT = *FORWARD_OR_ADJOINT_f;
  realw alpha_lddrk = *alpha_lddrk_f;
  realw beta_lddrk = *beta_lddrk_f;

  // safety check
  if (FORWARD_OR_ADJOINT != 1 && FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in compute_forces_outer_core_gpu() routine");
  }

  // determines number of elements to loop over (inner/outer elements)
  int num_elements;
  if (*iphase == 1) {
    num_elements = mp->nspec_outer_outer_core;
  } else {
    num_elements = mp->nspec_inner_outer_core;
  }

  // checks if anything to do
  if (num_elements == 0) return;

  // mesh coloring
  if (mp->use_mesh_coloring_gpu) {

    // note: array offsets require sorted arrays, such that e.g. ibool starts with elastic elements
    //         and followed by acoustic ones.
    //         acoustic elements also start with outer than inner element ordering

    int nb_colors, nb_blocks_to_compute;
    int istart;
    int offset, offset_nonpadded;

    // sets up color loop
    if (mp->NSPEC_OUTER_CORE > COLORING_MIN_NSPEC_OUTER_CORE) {
      if (*iphase == 1) {
        // outer elements
        nb_colors = mp->num_colors_outer_outer_core;
        istart = 0;

        // array offsets
        offset = 0;
        offset_nonpadded = 0;
      } else {
        // inner element colors (start after outer elements)
        nb_colors = mp->num_colors_outer_outer_core + mp->num_colors_inner_outer_core;
        istart = mp->num_colors_outer_outer_core;

        // array offsets (inner elements start after outer ones)
        offset = mp->nspec_outer_outer_core * NGLL3_PADDED;
        offset_nonpadded = mp->nspec_outer_outer_core * NGLL3;
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
      } else {
        // inner element colors (start after outer elements)
        nb_colors = 1;
        istart = 0;

        // array offsets (inner elements start after outer ones)
        offset = mp->nspec_outer_outer_core * NGLL3_PADDED;
        offset_nonpadded = mp->nspec_outer_outer_core * NGLL3;
      }
    }

    // loops over colors
    int icolor;

    for (icolor = istart; icolor < nb_colors; icolor++) {

      // gets number of elements for this color
      if (mp->NSPEC_OUTER_CORE > COLORING_MIN_NSPEC_OUTER_CORE) {
        nb_blocks_to_compute = mp->h_num_elem_colors_outer_core[icolor];
      } else {
        nb_blocks_to_compute = num_elements;
      }

      // debug
      //printf("debug: compute_forces_outer_core_gpu - iphase %d color %d out of %d, nspec elements %d, offset %d, offset_nonpadded % d, blocks %d\n",
      //       *iphase,icolor,nb_colors,mp->NSPEC_OUTER_CORE,offset,offset_nonpadded,nb_blocks_to_compute);

      INITIALIZE_OFFSET();

      //create cl_mem object _buffer + _ + _offset
      INIT_OFFSET(d_ibool_outer_core, offset_nonpadded);
      INIT_OFFSET(d_xix_outer_core, offset);
      INIT_OFFSET(d_xiy_outer_core, offset);
      INIT_OFFSET(d_xiz_outer_core, offset);
      INIT_OFFSET(d_etax_outer_core, offset);
      INIT_OFFSET(d_etay_outer_core, offset);
      INIT_OFFSET(d_etaz_outer_core, offset);
      INIT_OFFSET(d_gammax_outer_core, offset);
      INIT_OFFSET(d_gammay_outer_core, offset);
      INIT_OFFSET(d_gammaz_outer_core, offset);
      INIT_OFFSET(d_A_array_rotation, offset_nonpadded);
      INIT_OFFSET(d_B_array_rotation, offset_nonpadded);
      INIT_OFFSET(d_b_A_array_rotation, offset_nonpadded);
      INIT_OFFSET(d_b_B_array_rotation, offset_nonpadded);
      INIT_OFFSET(d_A_array_rotation_lddrk, offset_nonpadded);
      INIT_OFFSET(d_B_array_rotation_lddrk, offset_nonpadded);
      INIT_OFFSET(d_b_A_array_rotation_lddrk, offset_nonpadded);
      INIT_OFFSET(d_b_B_array_rotation_lddrk, offset_nonpadded);

      // debug
      //gpu_int_mem my_buffer = mp->d_ibool_outer_core;
      //my_buffer.cuda = my_buffer.cuda + offset_nonpadded;
      //printf("debug: d_ibool_outer_core %p offset_buffer %p my_buffer %p \n",
      //       (void *) mp->d_ibool_outer_core.cuda, (void *) d_ibool_outer_core_offset_nonpadded.cuda, (void*) my_buffer.cuda);

      outer_core (nb_blocks_to_compute, mp,
                  *iphase,
                  PASS_OFFSET(d_ibool_outer_core, offset_nonpadded),
                  PASS_OFFSET(d_xix_outer_core, offset),
                  PASS_OFFSET(d_xiy_outer_core, offset),
                  PASS_OFFSET(d_xiz_outer_core, offset),
                  PASS_OFFSET(d_etax_outer_core, offset),
                  PASS_OFFSET(d_etay_outer_core, offset),
                  PASS_OFFSET(d_etaz_outer_core, offset),
                  PASS_OFFSET(d_gammax_outer_core, offset),
                  PASS_OFFSET(d_gammay_outer_core, offset),
                  PASS_OFFSET(d_gammaz_outer_core, offset),
                  timeval,
                  PASS_OFFSET(d_A_array_rotation, offset_nonpadded),
                  PASS_OFFSET(d_B_array_rotation, offset_nonpadded),
                  PASS_OFFSET(d_b_A_array_rotation, offset_nonpadded),
                  PASS_OFFSET(d_b_B_array_rotation, offset_nonpadded),
                  PASS_OFFSET(d_A_array_rotation_lddrk, offset_nonpadded),
                  PASS_OFFSET(d_B_array_rotation_lddrk, offset_nonpadded),
                  PASS_OFFSET(d_b_A_array_rotation_lddrk, offset_nonpadded),
                  PASS_OFFSET(d_b_B_array_rotation_lddrk, offset_nonpadded),
                  alpha_lddrk,beta_lddrk,
                  FORWARD_OR_ADJOINT);

      RELEASE_OFFSET(d_ibool_outer_core, offset_nonpadded);
      RELEASE_OFFSET(d_xix_outer_core, offset);
      RELEASE_OFFSET(d_xiy_outer_core, offset);
      RELEASE_OFFSET(d_xiz_outer_core, offset);
      RELEASE_OFFSET(d_etax_outer_core, offset);
      RELEASE_OFFSET(d_etay_outer_core, offset);
      RELEASE_OFFSET(d_etaz_outer_core, offset);
      RELEASE_OFFSET(d_gammax_outer_core, offset);
      RELEASE_OFFSET(d_gammay_outer_core, offset);
      RELEASE_OFFSET(d_gammaz_outer_core, offset);
      RELEASE_OFFSET(d_A_array_rotation, offset_nonpadded);
      RELEASE_OFFSET(d_B_array_rotation, offset_nonpadded);
      RELEASE_OFFSET(d_b_A_array_rotation, offset_nonpadded);
      RELEASE_OFFSET(d_b_B_array_rotation, offset_nonpadded);
      RELEASE_OFFSET(d_A_array_rotation_lddrk, offset_nonpadded);
      RELEASE_OFFSET(d_B_array_rotation_lddrk, offset_nonpadded);
      RELEASE_OFFSET(d_b_A_array_rotation_lddrk, offset_nonpadded);
      RELEASE_OFFSET(d_b_B_array_rotation_lddrk, offset_nonpadded);

      // for padded and aligned arrays
      offset += nb_blocks_to_compute * NGLL3_PADDED;
      // for no-aligned arrays
      offset_nonpadded += nb_blocks_to_compute * NGLL3;
    }

  } else {

    // no mesh coloring: uses atomic updates
    outer_core (num_elements, mp,
                *iphase,
                mp->d_ibool_outer_core,
                mp->d_xix_outer_core,
                mp->d_xiy_outer_core,
                mp->d_xiz_outer_core,
                mp->d_etax_outer_core,
                mp->d_etay_outer_core,
                mp->d_etaz_outer_core,
                mp->d_gammax_outer_core,
                mp->d_gammay_outer_core,
                mp->d_gammaz_outer_core,
                timeval,
                mp->d_A_array_rotation,
                mp->d_B_array_rotation,
                mp->d_b_A_array_rotation,
                mp->d_b_B_array_rotation,
                mp->d_A_array_rotation_lddrk,
                mp->d_B_array_rotation_lddrk,
                mp->d_b_A_array_rotation_lddrk,
                mp->d_b_B_array_rotation_lddrk,
                alpha_lddrk,beta_lddrk,
                FORWARD_OR_ADJOINT);
  }

  //double end_time = get_time_val ();
  //printf ("Elapsed time: %e\n", end_time-start_time);
  GPU_ERROR_CHECKING ("compute_forces_outer_core_gpu");
}
