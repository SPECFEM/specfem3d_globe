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

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_strain_gpu,
               COMPUTE_STRAIN_GPU) (long *Mesh_pointer_f, realw *deltat_f, int *FORWARD_OR_ADJOINT) {

  TRACE ("compute_strain_gpu");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int size,size_strain_only;
  int blocksize = NGLL3;

  realw deltat = *deltat_f;

  gpu_realw_mem displ,veloc;
  gpu_realw_mem epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz;
  gpu_realw_mem eps_trace_over_3;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in compute_strain_gpu() routine");
  }

  // checks flag
  if (! mp->undo_attenuation) exit_on_error("Error invalid UNDO_ATTENUATION flag in compute_strain_gpu() routine");

  // crust/mantle
  size = mp->NSPEC_CRUST_MANTLE;
  size_strain_only = mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY;

  // blocks
  int num_blocks_x, num_blocks_y;
  get_blocks_xy (size, &num_blocks_x, &num_blocks_y);

  // sets gpu arrays
  if (*FORWARD_OR_ADJOINT == 1) {
    displ = mp->d_displ_crust_mantle;
    veloc = mp->d_veloc_crust_mantle;
    epsilondev_xx = mp->d_epsilondev_xx_crust_mantle;
    epsilondev_yy = mp->d_epsilondev_yy_crust_mantle;
    epsilondev_xy = mp->d_epsilondev_xy_crust_mantle;
    epsilondev_xz = mp->d_epsilondev_xz_crust_mantle;
    epsilondev_yz = mp->d_epsilondev_yz_crust_mantle;
    eps_trace_over_3 = mp->d_eps_trace_over_3_crust_mantle;
  } else {
    // for backward/reconstructed fields
    displ = mp->d_b_displ_crust_mantle;
    veloc = mp->d_b_veloc_crust_mantle;
    epsilondev_xx = mp->d_b_epsilondev_xx_crust_mantle;
    epsilondev_yy = mp->d_b_epsilondev_yy_crust_mantle;
    epsilondev_xy = mp->d_b_epsilondev_xy_crust_mantle;
    epsilondev_xz = mp->d_b_epsilondev_xz_crust_mantle;
    epsilondev_yz = mp->d_b_epsilondev_yz_crust_mantle;
    eps_trace_over_3 = mp->d_b_eps_trace_over_3_crust_mantle;
  }

#ifdef USE_OPENCL
  size_t global_work_size[2];
  size_t local_work_size[2];
  cl_uint idx;

  if (run_opencl) {
    idx = 0;
    // strain kernel
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &displ.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &veloc.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &epsilondev_xx.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &epsilondev_yy.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &epsilondev_xy.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &epsilondev_xz.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &epsilondev_yz.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &eps_trace_over_3.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (int), (void *) &size));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (int), (void *) &size_strain_only));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (realw), (void *) &deltat));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xix_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xiy_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xiz_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etax_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etay_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etaz_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammax_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammay_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammaz_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_hprime_xx.ocl));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_strain_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  dim3 grid,threads;

  if (run_cuda) {
    grid = dim3(num_blocks_x,num_blocks_y);
    threads = dim3(blocksize,1,1);

    // strain kernel
    compute_strain_kernel<<<grid,threads,0,mp->compute_stream>>>(displ.cuda,
                                                                 veloc.cuda,
                                                                 epsilondev_xx.cuda,
                                                                 epsilondev_yy.cuda,
                                                                 epsilondev_xy.cuda,
                                                                 epsilondev_xz.cuda,
                                                                 epsilondev_yz.cuda,
                                                                 eps_trace_over_3.cuda,
                                                                 size,size_strain_only,
                                                                 deltat,
                                                                 mp->d_ibool_crust_mantle.cuda,
                                                                 mp->d_xix_crust_mantle.cuda,
                                                                 mp->d_xiy_crust_mantle.cuda,
                                                                 mp->d_xiz_crust_mantle.cuda,
                                                                 mp->d_etax_crust_mantle.cuda,
                                                                 mp->d_etay_crust_mantle.cuda,
                                                                 mp->d_etaz_crust_mantle.cuda,
                                                                 mp->d_gammax_crust_mantle.cuda,
                                                                 mp->d_gammay_crust_mantle.cuda,
                                                                 mp->d_gammaz_crust_mantle.cuda,
                                                                 mp->d_hprime_xx.cuda);
  }
#endif

  // inner core
  size = mp->NSPEC_INNER_CORE;
  size_strain_only = mp->NSPEC_INNER_CORE_STRAIN_ONLY;

  // blocks
  get_blocks_xy (size, &num_blocks_x, &num_blocks_y);

  // sets gpu arrays
  if (*FORWARD_OR_ADJOINT == 1) {
    displ = mp->d_displ_inner_core;
    veloc = mp->d_veloc_inner_core;
    epsilondev_xx = mp->d_epsilondev_xx_inner_core;
    epsilondev_yy = mp->d_epsilondev_yy_inner_core;
    epsilondev_xy = mp->d_epsilondev_xy_inner_core;
    epsilondev_xz = mp->d_epsilondev_xz_inner_core;
    epsilondev_yz = mp->d_epsilondev_yz_inner_core;
    eps_trace_over_3 = mp->d_eps_trace_over_3_inner_core;
  } else {
    // for backward/reconstructed fields
    displ = mp->d_b_displ_inner_core;
    veloc = mp->d_b_veloc_inner_core;
    epsilondev_xx = mp->d_b_epsilondev_xx_inner_core;
    epsilondev_yy = mp->d_b_epsilondev_yy_inner_core;
    epsilondev_xy = mp->d_b_epsilondev_xy_inner_core;
    epsilondev_xz = mp->d_b_epsilondev_xz_inner_core;
    epsilondev_yz = mp->d_b_epsilondev_yz_inner_core;
    eps_trace_over_3 = mp->d_b_eps_trace_over_3_inner_core;
  }

#ifdef USE_OPENCL
  if (run_opencl) {
    idx = 0;
    // strain kernel
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &displ.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &veloc.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &epsilondev_xx.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &epsilondev_yy.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &epsilondev_xy.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &epsilondev_xz.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &epsilondev_yz.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &eps_trace_over_3.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (int), (void *) &size));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (int), (void *) &size_strain_only));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (realw), (void *) &deltat));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xix_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xiy_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xiz_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etax_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etay_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etaz_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammax_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammay_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammaz_inner_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_strain_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_hprime_xx.ocl));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_strain_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    grid = dim3(num_blocks_x,num_blocks_y);
    threads = dim3(blocksize,1,1);

    // strain kernel
    compute_strain_kernel<<<grid,threads,0,mp->compute_stream>>>(displ.cuda,
                                                                 veloc.cuda,
                                                                 epsilondev_xx.cuda,
                                                                 epsilondev_yy.cuda,
                                                                 epsilondev_xy.cuda,
                                                                 epsilondev_xz.cuda,
                                                                 epsilondev_yz.cuda,
                                                                 eps_trace_over_3.cuda,
                                                                 size,size_strain_only,
                                                                 deltat,
                                                                 mp->d_ibool_inner_core.cuda,
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

  GPU_ERROR_CHECKING ("compute_strain_gpu");
}

