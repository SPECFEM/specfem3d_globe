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
#include "smooth_gpu.h"


/* ----------------------------------------------------------------------------------------------- */

// smoothing routines

/* ----------------------------------------------------------------------------------------------- */

void gpuInitialize_smooth_buffers(Smooth_data *sp) {

  TRACE ("gpuInitialize_smooth_buffers");

#ifdef USE_OPENCL
  // sets OpenCL pointers to NULL
  #define INIT_DUMMY_BUFFER(_field_) sp->_field_.ocl = NULL;
  #define GPU_REALW_BUFFER INIT_DUMMY_BUFFER
#endif
#ifdef USE_CUDA
  // sets CUDA pointers to NULL
  #define INIT_DUMMY_BUFFER(_field_) sp->_field_.cuda = NULL
  #define GPU_REALW_BUFFER INIT_DUMMY_BUFFER
#endif
#ifdef USE_HIP
  // sets HIP pointers to NULL
  #define INIT_DUMMY_BUFFER(_field_) sp->_field_.hip = NULL
  #define GPU_REALW_BUFFER INIT_DUMMY_BUFFER
#endif

  GPU_REALW_BUFFER (x_me);
  GPU_REALW_BUFFER (y_me);
  GPU_REALW_BUFFER (z_me);

  GPU_REALW_BUFFER (x_other);
  GPU_REALW_BUFFER (y_other);
  GPU_REALW_BUFFER (z_other);

  GPU_REALW_BUFFER (integ_factor);
  GPU_REALW_BUFFER (kernel);

  GPU_REALW_BUFFER (data_smooth);
  GPU_REALW_BUFFER (normalisation);

#undef INIT_DUMMY_BUFFER
}

/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(prepare_smooth_gpu,
              PREPARE_SMOOTH_GPU)(long * Container,
                                  realw * xstore_me,
                                  realw * ystore_me,
                                  realw * zstore_me,
                                  realw * sigma_h2,
                                  realw * sigma_v2,
                                  realw * sigma_h3,
                                  realw * sigma_v3,
                                  int * nspec_me,
                                  int * nker){

  TRACE("prepare_smooth_gpu");

  // allocates structure
  Smooth_data* sp = (Smooth_data*) malloc( sizeof(Smooth_data) );

  // checks pointer
  if (! sp) exit_on_error ("Error allocating smooth data pointer");

  // sets fortran pointer
  *Container = (long)sp;

  // initializes gpu array pointers
  gpuInitialize_smooth_buffers(sp);

  // copies arrays to GPU
  gpuCreateCopy_todevice_realw(&sp->x_me, xstore_me, NGLL3*(*nspec_me));
  gpuCreateCopy_todevice_realw(&sp->y_me, ystore_me, NGLL3*(*nspec_me));
  gpuCreateCopy_todevice_realw(&sp->z_me, zstore_me, NGLL3*(*nspec_me));

  // sets variable values
  realw sigma_h2_inv = 1.0f / (*sigma_h2);
  realw sigma_v2_inv = 1.0f / (*sigma_v2);

  realw sigma_h3_sq = (*sigma_h3) * (*sigma_h3);  // squared
  realw sigma_v3_sq = (*sigma_v3) * (*sigma_v3);

  sp->sigma_h2_inv = sigma_h2_inv;
  sp->sigma_v2_inv = sigma_v2_inv;

  sp->h_criterion = sigma_h3_sq;
  sp->v_criterion = sigma_v3_sq;

  sp->nspec_me = *nspec_me;
  sp->nker = *nker;

  // smoothed kernel
  gpuMalloc_realw(&sp->data_smooth, NGLL3*(*nspec_me)*(*nker));
  // initializes values to zero
  gpuMemset_realw(&sp->data_smooth, NGLL3*(*nspec_me)*(*nker), 0);

  // normalization factor
  gpuMalloc_realw(&sp->normalisation, NGLL3*(*nspec_me));
  // initializes values to zero
  gpuMemset_realw(&sp->normalisation, NGLL3*(*nspec_me), 0);

  GPU_ERROR_CHECKING ("prepare_smooth_gpu");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(compute_smooth_gpu,
              COMPUTE_SMOOTH_GPU)(long * smooth_pointer,
                                  realw * data_other,
                                  realw * integ_factor,
                                  realw * xstore_other,
                                  realw * ystore_other,
                                  realw * zstore_other,
                                  const int * nspec_other_f,
                                  const int * use_vector_distance_f){

  TRACE("compute_smooth_gpu");

  Smooth_data *sp = (Smooth_data*) *smooth_pointer;

  int nspec_other = *nspec_other_f;
  int use_vector_distance = *use_vector_distance_f;

  gpuCreateCopy_todevice_realw(&sp->x_other,xstore_other,NGLL3 * nspec_other);
  gpuCreateCopy_todevice_realw(&sp->y_other,ystore_other,NGLL3 * nspec_other);
  gpuCreateCopy_todevice_realw(&sp->z_other,zstore_other,NGLL3 * nspec_other);

  gpuCreateCopy_todevice_realw(&sp->integ_factor,integ_factor,NGLL3 * nspec_other);

  // smoothed kernel
  gpuMalloc_realw(&sp->kernel, NGLL3 * nspec_other);

  //dim3 grid(sp->nspec_me,1);
  //dim3 threads(NGLL3,1,1);

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (sp->nspec_me, &num_blocks_x, &num_blocks_y);

  // loops over kernels
  for (int iker=0; iker < sp->nker; iker++){
    // pointer to corresponding kernel section
    realw * data_p = &data_other[NGLL3 * nspec_other * iker];

    // copy single kernel
    gpuCopy_todevice_realw(&sp->kernel,data_p,NGLL3 * nspec_other);

    // runs gaussian smoothing
#ifdef USE_OPENCL
    if (run_opencl) {
      size_t global_work_size[2];
      size_t local_work_size[2];
      cl_uint idx = 0;

      clCheck (clSetKernelArg (mocl.kernels.smooth_process_kernel, idx++, sizeof (cl_mem), (void *) &sp->x_me.ocl));
      clCheck (clSetKernelArg (mocl.kernels.smooth_process_kernel, idx++, sizeof (cl_mem), (void *) &sp->y_me.ocl));
      clCheck (clSetKernelArg (mocl.kernels.smooth_process_kernel, idx++, sizeof (cl_mem), (void *) &sp->z_me.ocl));
      clCheck (clSetKernelArg (mocl.kernels.smooth_process_kernel, idx++, sizeof (cl_mem), (void *) &sp->x_other.ocl));
      clCheck (clSetKernelArg (mocl.kernels.smooth_process_kernel, idx++, sizeof (cl_mem), (void *) &sp->y_other.ocl));
      clCheck (clSetKernelArg (mocl.kernels.smooth_process_kernel, idx++, sizeof (cl_mem), (void *) &sp->z_other.ocl));
      clCheck (clSetKernelArg (mocl.kernels.smooth_process_kernel, idx++, sizeof (cl_mem), (void *) &sp->kernel.ocl));
      clCheck (clSetKernelArg (mocl.kernels.smooth_process_kernel, idx++, sizeof (realw),  (void *) &sp->sigma_h2_inv));
      clCheck (clSetKernelArg (mocl.kernels.smooth_process_kernel, idx++, sizeof (realw),  (void *) &sp->sigma_v2_inv));
      clCheck (clSetKernelArg (mocl.kernels.smooth_process_kernel, idx++, sizeof (int),    (void *) &iker));
      clCheck (clSetKernelArg (mocl.kernels.smooth_process_kernel, idx++, sizeof (int),    (void *) &sp->nspec_me));
      clCheck (clSetKernelArg (mocl.kernels.smooth_process_kernel, idx++, sizeof (int),    (void *) &nspec_other));
      clCheck (clSetKernelArg (mocl.kernels.smooth_process_kernel, idx++, sizeof (realw),  (void *) &sp->v_criterion));
      clCheck (clSetKernelArg (mocl.kernels.smooth_process_kernel, idx++, sizeof (realw),  (void *) &sp->h_criterion));
      clCheck (clSetKernelArg (mocl.kernels.smooth_process_kernel, idx++, sizeof (cl_mem), (void *) &sp->integ_factor.ocl));
      clCheck (clSetKernelArg (mocl.kernels.smooth_process_kernel, idx++, sizeof (cl_mem), (void *) &sp->data_smooth.ocl));
      clCheck (clSetKernelArg (mocl.kernels.smooth_process_kernel, idx++, sizeof (cl_mem), (void *) &sp->normalisation.ocl));
      clCheck (clSetKernelArg (mocl.kernels.smooth_process_kernel, idx++, sizeof (int),    (void *) &use_vector_distance));

      local_work_size[0] = NGLL3;
      local_work_size[1] = 1;
      global_work_size[0] = num_blocks_x * NGLL3;
      global_work_size[1] = num_blocks_y;

      clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.smooth_process_kernel, 2, NULL,
                                       global_work_size, local_work_size, 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda){
      dim3 grid(num_blocks_x,num_blocks_y);
      dim3 threads(NGLL3,1,1);

      smooth_process_kernel<<<grid,threads>>>(sp->x_me.cuda,
                                              sp->y_me.cuda,
                                              sp->z_me.cuda,
                                              sp->x_other.cuda,
                                              sp->y_other.cuda,
                                              sp->z_other.cuda,
                                              sp->kernel.cuda,
                                              sp->sigma_h2_inv,
                                              sp->sigma_v2_inv,
                                              iker,
                                              sp->nspec_me,
                                              nspec_other,
                                              sp->v_criterion,
                                              sp->h_criterion,
                                              sp->integ_factor.cuda,
                                              sp->data_smooth.cuda,
                                              sp->normalisation.cuda,
                                              use_vector_distance);
    }
#endif
#ifdef USE_HIP
    if (run_hip){
      dim3 grid(num_blocks_x,num_blocks_y);
      dim3 threads(NGLL3,1,1);

      hipLaunchKernelGGL(smooth_process_kernel, dim3(grid), dim3(threads), 0, 0,
                                                 sp->x_me.hip,
                                                 sp->y_me.hip,
                                                 sp->z_me.hip,
                                                 sp->x_other.hip,
                                                 sp->y_other.hip,
                                                 sp->z_other.hip,
                                                 sp->kernel.hip,
                                                 sp->sigma_h2_inv,
                                                 sp->sigma_v2_inv,
                                                 iker,
                                                 sp->nspec_me,
                                                 nspec_other,
                                                 sp->v_criterion,
                                                 sp->h_criterion,
                                                 sp->integ_factor.hip,
                                                 sp->data_smooth.hip,
                                                 sp->normalisation.hip,
                                                 use_vector_distance);
    }
#endif
  }
  gpuSynchronize();

  gpuFree(&sp->kernel);

  gpuFree(&sp->x_other);
  gpuFree(&sp->y_other);
  gpuFree(&sp->z_other);
  gpuFree(&sp->integ_factor);

  GPU_ERROR_CHECKING ("compute_smooth_gpu");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(get_smooth_gpu,
              GET_SMOOTH_GPU)(long * smooth_pointer,
                              realw * data_smooth) {

  TRACE("get_smooth");

  Smooth_data *sp = (Smooth_data*) *smooth_pointer;

  //dim3 grid(sp->nspec_me,1);
  //dim3 threads(NGLL3,1,1);

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (sp->nspec_me, &num_blocks_x, &num_blocks_y);

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    clCheck (clSetKernelArg (mocl.kernels.smooth_normalize_data_kernel, idx++, sizeof (cl_mem), (void *) &sp->data_smooth.ocl));
    clCheck (clSetKernelArg (mocl.kernels.smooth_normalize_data_kernel, idx++, sizeof (cl_mem), (void *) &sp->normalisation.ocl));
    clCheck (clSetKernelArg (mocl.kernels.smooth_normalize_data_kernel, idx++, sizeof (int),    (void *) &sp->nker));
    clCheck (clSetKernelArg (mocl.kernels.smooth_normalize_data_kernel, idx++, sizeof (int),    (void *) &sp->nspec_me));

    local_work_size[0] = NGLL3;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * NGLL3;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.smooth_normalize_data_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
}
#endif
#ifdef USE_CUDA
  if (run_cuda){
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(NGLL3,1,1);

    smooth_normalize_data_kernel<<<grid,threads>>>(sp->data_smooth.cuda,
                                                   sp->normalisation.cuda,
                                                   sp->nker,
                                                   sp->nspec_me);
  }
#endif
#ifdef USE_HIP
  if (run_hip){
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(NGLL3,1,1);

    hipLaunchKernelGGL(smooth_normalize_data_kernel, dim3(grid), dim3(threads), 0, 0,
                                                     sp->data_smooth.hip,
                                                     sp->normalisation.hip,
                                                     sp->nker,sp->nspec_me);
  }
#endif

  gpuCopy_from_device_realw(&sp->data_smooth, data_smooth, NGLL3 * sp->nspec_me * sp->nker);

  gpuFree(&sp->x_me);
  gpuFree(&sp->y_me);
  gpuFree(&sp->z_me);
  gpuFree(&sp->data_smooth);
  gpuFree(&sp->normalisation);

  free(sp);

  GPU_ERROR_CHECKING ("get_smooth_gpu");
}

