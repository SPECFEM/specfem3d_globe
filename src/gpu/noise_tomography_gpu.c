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
// noise transfer surface movie
/*----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_ (noise_transfer_surface_to_host,
               NOISE_TRANSFER_SURFACE_TO_HOST) (long *Mesh_pointer_f,
                                                realw *h_noise_surface_movie) {
  TRACE ("noise_transfer_surface_to_host");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->nspec2D_top_crust_mantle, &num_blocks_x, &num_blocks_y);

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    clCheck (clSetKernelArg (mocl.kernels.noise_transfer_surface_to_host_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibelm_top_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.noise_transfer_surface_to_host_kernel, idx++, sizeof (int), (void *) &mp->nspec2D_top_crust_mantle));
    clCheck (clSetKernelArg (mocl.kernels.noise_transfer_surface_to_host_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.noise_transfer_surface_to_host_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_displ_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.noise_transfer_surface_to_host_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_noise_surface_movie.ocl));

    local_work_size[0] = NGLL2;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * NGLL2;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.noise_transfer_surface_to_host_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y,1);
    dim3 threads(NGLL2,1,1);

    noise_transfer_surface_to_host_kernel<<<grid,threads>>>(mp->d_ibelm_top_crust_mantle.cuda,
                                                            mp->nspec2D_top_crust_mantle,
                                                            mp->d_ibool_crust_mantle.cuda,
                                                            mp->d_displ_crust_mantle.cuda,
                                                            mp->d_noise_surface_movie.cuda);
  }
#endif

  // copies noise array to CPU
  gpuCopy_from_device_realw (&mp->d_noise_surface_movie, h_noise_surface_movie, NDIM * NGLL2 * mp->nspec2D_top_crust_mantle);

  GPU_ERROR_CHECKING ("noise_transfer_surface_to_host");
}

/*----------------------------------------------------------------------------------------------- */
// NOISE add source master
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (noise_add_source_master_rec_gpu,
               NOISE_ADD_SOURCE_MASTER_REC_GPU) (long *Mesh_pointer_f,
                                                 int *it_f,
                                                 int *irec_master_noise_f,
                                                 int *islice_selected_rec) {

  TRACE ("noise_add_source_master_rec_cu");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int it = *it_f - 1;   // -1 for Fortran -> C indexing differences
  int irec_master_noise = *irec_master_noise_f-1;

  // adds noise source at master location
  if (mp->myrank == islice_selected_rec[irec_master_noise]) {
#ifdef USE_OPENCL
    if (run_opencl) {
      size_t global_work_size[2];
      size_t local_work_size[2];
      cl_uint idx = 0;

      clCheck (clSetKernelArg (mocl.kernels.noise_add_source_master_rec_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_source_master_rec_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ispec_selected_rec.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_source_master_rec_kernel, idx++, sizeof (int), (void *) &irec_master_noise));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_source_master_rec_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_accel_crust_mantle.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_source_master_rec_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_noise_sourcearray.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_source_master_rec_kernel, idx++, sizeof (int), (void *) &it));

      local_work_size[0] = NGLL3;
      local_work_size[1] = 1;
      global_work_size[0] = 1 * NGLL3;
      global_work_size[1] = 1;

      clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.noise_add_source_master_rec_kernel, 2, NULL,
                                       global_work_size, local_work_size, 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      dim3 grid(1,1,1);
      dim3 threads(NGLL3,1,1);

      noise_add_source_master_rec_kernel<<<grid,threads>>>(mp->d_ibool_crust_mantle.cuda,
                                                           mp->d_ispec_selected_rec.cuda,
                                                           irec_master_noise,
                                                           mp->d_accel_crust_mantle.cuda,
                                                           mp->d_noise_sourcearray.cuda,
                                                           it);
    }
#endif
  }

  GPU_ERROR_CHECKING ("noise_add_source_master_rec_kernel");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (noise_add_surface_movie_gpu,
               NOISE_ADD_SURFACE_MOVIE_GPU) (long *Mesh_pointer_f,
                                             realw *h_noise_surface_movie) {

  TRACE ("noise_add_surface_movie_gpu");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->nspec2D_top_crust_mantle, &num_blocks_x, &num_blocks_y);

#ifdef USE_OPENCL
  size_t global_work_size[2];
  size_t local_work_size[2];
  cl_uint idx = 0;
#endif
#ifdef USE_CUDA
  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(NGLL2,1,1);
#endif

  // copies surface movie to GPU
  gpuCopy_todevice_realw (&mp->d_noise_surface_movie, h_noise_surface_movie, NDIM*NGLL2 *(mp->nspec2D_top_crust_mantle));

  switch (mp->noise_tomography) {
  case 2:
    // adds surface source to forward field
#ifdef USE_OPENCL
    if (run_opencl) {
      idx = 0;
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_accel_crust_mantle.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibelm_top_crust_mantle.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (int), (void *) &mp->nspec2D_top_crust_mantle));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_noise_surface_movie.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_normal_x_noise.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_normal_y_noise.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_normal_z_noise.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_mask_noise.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_jacobian2D_top_crust_mantle.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_wgllwgll_xy.ocl));

      local_work_size[0] = NGLL3;
      local_work_size[1] = 1;
      global_work_size[0] = num_blocks_x * NGLL3;
      global_work_size[1] = num_blocks_y;

      clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.noise_add_surface_movie_kernel, 2, NULL,
                                       global_work_size, local_work_size, 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      noise_add_surface_movie_kernel<<<grid,threads>>>(mp->d_accel_crust_mantle.cuda,
                                                       mp->d_ibool_crust_mantle.cuda,
                                                       mp->d_ibelm_top_crust_mantle.cuda,
                                                       mp->nspec2D_top_crust_mantle,
                                                       mp->d_noise_surface_movie.cuda,
                                                       mp->d_normal_x_noise.cuda,
                                                       mp->d_normal_y_noise.cuda,
                                                       mp->d_normal_z_noise.cuda,
                                                       mp->d_mask_noise.cuda,
                                                       mp->d_jacobian2D_top_crust_mantle.cuda,
                                                       mp->d_wgllwgll_xy.cuda);
    }
#endif
    break;

  case 3:
    // adds surface source to adjoint (backward) field
#ifdef USE_OPENCL
    if (run_opencl) {
      idx = 0;
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_accel_crust_mantle.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibelm_top_crust_mantle.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (int), (void *) &mp->nspec2D_top_crust_mantle));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_noise_surface_movie.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_normal_x_noise.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_normal_y_noise.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_normal_z_noise.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_mask_noise.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_jacobian2D_top_crust_mantle.ocl));
      clCheck (clSetKernelArg (mocl.kernels.noise_add_surface_movie_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_wgllwgll_xy.ocl));

      local_work_size[0] = NGLL3;
      local_work_size[1] = 1;
      global_work_size[0] = num_blocks_x * NGLL3;
      global_work_size[1] = num_blocks_y;

      clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.noise_add_surface_movie_kernel, 2, NULL,
                                       global_work_size, local_work_size, 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      noise_add_surface_movie_kernel<<<grid,threads>>>(mp->d_b_accel_crust_mantle.cuda,
                                                       mp->d_ibool_crust_mantle.cuda,
                                                       mp->d_ibelm_top_crust_mantle.cuda,
                                                       mp->nspec2D_top_crust_mantle,
                                                       mp->d_noise_surface_movie.cuda,
                                                       mp->d_normal_x_noise.cuda,
                                                       mp->d_normal_y_noise.cuda,
                                                       mp->d_normal_z_noise.cuda,
                                                       mp->d_mask_noise.cuda,
                                                       mp->d_jacobian2D_top_crust_mantle.cuda,
                                                       mp->d_wgllwgll_xy.cuda);
    }
#endif
    break;
  }

  GPU_ERROR_CHECKING ("noise_read_add_surface_movie_kernel");
}
