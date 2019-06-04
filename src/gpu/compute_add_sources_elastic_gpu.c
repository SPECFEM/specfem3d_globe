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

extern EXTERN_LANG
void FC_FUNC_ (compute_add_sources_gpu,
               COMPUTE_ADD_SOURCES_GPU) (long *Mesh_pointer_f,
                                         int *NSOURCESf,
                                         double *h_stf_pre_compute) {

  TRACE ("compute_add_sources_gpu");

  // get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // checks if anything to do
  if (mp->nsources_local == 0) return;

  int NSOURCES = *NSOURCESf;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (NSOURCES, &num_blocks_x, &num_blocks_y);

  // copies source time function buffer values to GPU
  gpuCopy_todevice_double (&mp->d_stf_pre_compute, h_stf_pre_compute, NSOURCES);

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[3];
    size_t local_work_size[3];
    cl_uint idx = 0;

    // adds source contributions
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_accel_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_sourcearrays.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_stf_pre_compute.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_kernel, idx++, sizeof (int), (void *) &mp->myrank));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_islice_selected_source.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ispec_selected_source.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_kernel, idx++, sizeof (int), (void *) &NSOURCES));

    local_work_size[0] = NGLLX;
    local_work_size[1] = NGLLX;
    local_work_size[2] = NGLLX;
    global_work_size[0] = num_blocks_x * NGLLX;
    global_work_size[1] = num_blocks_y * NGLLX;
    global_work_size[2] = NGLLX;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_add_sources_kernel, 3, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(NGLLX,NGLLX,NGLLX);
    // adds source contributions
    compute_add_sources_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel_crust_mantle.cuda,
                                                                      mp->d_ibool_crust_mantle.cuda,
                                                                      mp->d_sourcearrays.cuda,
                                                                      mp->d_stf_pre_compute.cuda,
                                                                      mp->myrank,
                                                                      mp->d_islice_selected_source.cuda,
                                                                      mp->d_ispec_selected_source.cuda,
                                                                      NSOURCES);
  }
#endif

  GPU_ERROR_CHECKING ("compute_add_sources_gpu");
}

/*----------------------------------------------------------------------------------------------- */
// backward sources
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_add_sources_backward_gpu,
               COMPUTE_ADD_SOURCES_BACKWARD_GPU) (long *Mesh_pointer_f,
                                                  int *NSOURCESf,
                                                  double *h_stf_pre_compute) {
  TRACE ("compute_add_sources_backward_gpu");
  // debug
  DEBUG_BACKWARD_SOURCES ();

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // checks if anything to do
  if (mp->nsources_local == 0) return;

  int NSOURCES = *NSOURCESf;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (NSOURCES, &num_blocks_x, &num_blocks_y);

  // copies source time function buffer values to GPU
  gpuCopy_todevice_double (&mp->d_stf_pre_compute, h_stf_pre_compute, NSOURCES);

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[3];
    size_t local_work_size[3];
    cl_uint idx = 0;

    // adds source contributions
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_b_accel_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_sourcearrays.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_stf_pre_compute.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_kernel, idx++, sizeof (int), (void *) &mp->myrank));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_islice_selected_source.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ispec_selected_source.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_kernel, idx++, sizeof (int), (void *) &NSOURCES));

    local_work_size[0] = NGLLX;
    local_work_size[1] = NGLLX;
    local_work_size[2] = NGLLX;
    global_work_size[0] = num_blocks_x * NGLLX;
    global_work_size[1] = num_blocks_y * NGLLX;
    global_work_size[2] = NGLLX;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_add_sources_kernel, 3, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(NGLLX,NGLLX,NGLLX);

    compute_add_sources_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_accel_crust_mantle.cuda,
                                                                      mp->d_ibool_crust_mantle.cuda,
                                                                      mp->d_sourcearrays.cuda,
                                                                      mp->d_stf_pre_compute.cuda,
                                                                      mp->myrank,
                                                                      mp->d_islice_selected_source.cuda,
                                                                      mp->d_ispec_selected_source.cuda,
                                                                      NSOURCES);
  }
#endif

  GPU_ERROR_CHECKING ("compute_add_sources_backward_gpu");
}


/*----------------------------------------------------------------------------------------------- */
// ADJOINT sources
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (compute_add_sources_adjoint_gpu,
               COMPUTE_ADD_SOURCES_ADJOINT_GPU) (long *Mesh_pointer_f) {

  // adds adjoint sources
  // note: call this routine after transfer_adj_to_device**() to have correct adjoint sourcearrays in array d_source_adjoint

  TRACE("compute_add_sources_adjoint_gpu");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // check if anything to do
  if (mp->nadj_rec_local == 0 ) return;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->nadj_rec_local, &num_blocks_x, &num_blocks_y);

  // the irec_local variable needs to be precomputed (as
  // h_pre_comp..), because normally it is in the loop updating accel,
  // and due to how it's incremented, it cannot be parallelized
#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[3];
    size_t local_work_size[3];
    cl_uint idx = 0;
    cl_event *copy_evt = NULL;
    cl_uint num_evt = 0;

    if (GPU_ASYNC_COPY) {
      // waits for previous copy to finish
      clCheck (clFinish (mocl.copy_queue));
      if (mp->has_last_copy_evt) {
        copy_evt = &mp->last_copy_evt;
        num_evt = 1;
      }
    }

    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_adjoint_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_accel_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_adjoint_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_source_adjoint.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_adjoint_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_hxir_adj.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_adjoint_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_hetar_adj.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_adjoint_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_hgammar_adj.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_adjoint_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_adjoint_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ispec_selected_rec.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_adjoint_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_number_adjsources_global.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_add_sources_adjoint_kernel, idx++, sizeof (int), (void *) &mp->nadj_rec_local));

    local_work_size[0] = NGLLX;
    local_work_size[1] = NGLLX;
    local_work_size[2] = NGLLX;

    global_work_size[0] = num_blocks_x * NGLLX;
    global_work_size[1] = num_blocks_y * NGLLX;
    global_work_size[2] = NGLLX;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_add_sources_adjoint_kernel, 3, NULL,
                                     global_work_size, local_work_size, num_evt, copy_evt, NULL));

    if (GPU_ASYNC_COPY) {
      if (mp->has_last_copy_evt) {
        clCheck (clReleaseEvent (mp->last_copy_evt));
        mp->has_last_copy_evt = 0;
      }
    }
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    // waits for previous transfer_** calls to be finished
    if (GPU_ASYNC_COPY) {
      // waits for asynchronous copy to finish
      cudaStreamSynchronize(mp->copy_stream);
    }

    dim3 grid(num_blocks_x,num_blocks_y,1);
    dim3 threads(NGLLX,NGLLX,NGLLX);

    compute_add_sources_adjoint_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel_crust_mantle.cuda,
                                                                              mp->d_source_adjoint.cuda,
                                                                              mp->d_hxir_adj.cuda,
                                                                              mp->d_hetar_adj.cuda,
                                                                              mp->d_hgammar_adj.cuda,
                                                                              mp->d_ibool_crust_mantle.cuda,
                                                                              mp->d_ispec_selected_rec.cuda,
                                                                              mp->d_number_adjsources_global.cuda,
                                                                              mp->nadj_rec_local);
  }
#endif

  GPU_ERROR_CHECKING("compute_add_sources_adjoint_cuda");
}

/* ----------------------------------------------------------------------------------------------- */
// adjoint memory transfers
/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_adj_to_device,
              TRANSFER_ADJ_TO_DEVICE)(long* Mesh_pointer_f,
                                      int* h_nrec,
                                      realw* h_source_adjoint,
                                      int* h_islice_selected_rec) {

  // transfers adjoint source arrays synchronously to GPU

  TRACE("transfer_adj_to_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // check if anything to do
  if (mp->nadj_rec_local == 0) return;

  // total number of receivers/adjoint sources
  int nrec = *h_nrec;

  int irec,irec_local;

  irec_local = 0;
  for (irec = 0; irec < nrec; irec++) {
    if (mp->myrank == h_islice_selected_rec[irec]) {
      // takes only local sources
        mp->h_source_adjoint[INDEX2(NDIM,0,irec_local)]
              = h_source_adjoint[INDEX2(NDIM,0,irec_local)];

        mp->h_source_adjoint[INDEX2(NDIM,1,irec_local)]
              = h_source_adjoint[INDEX2(NDIM,1,irec_local)];

        mp->h_source_adjoint[INDEX2(NDIM,2,irec_local)]
              = h_source_adjoint[INDEX2(NDIM,2,irec_local)];

      // increases local receivers counter
      irec_local++;
    }
  }
  // check all local sources were added
  if (irec_local != mp->nadj_rec_local) {
    exit_on_error ("irec_local not equal to nadj_rec_local\n");
  }

  // copies extracted array values onto GPU
  gpuCopy_todevice_realw (&mp->d_source_adjoint, mp->h_source_adjoint, mp->nadj_rec_local * NDIM );

  GPU_ERROR_CHECKING ("compute_add_sources_adjoint_gpu");
}


/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(transfer_adj_to_device_async,
              TRANSFER_ADJ_TO_DEVICE_ASYNC)(long *Mesh_pointer_f,
                                            int *h_nrec,
                                            realw *h_source_adjoint,
                                            int *h_islice_selected_rec) {
  // asynchronous transfer for next adjoint source arrays from host to device

  TRACE("transfer_adj_to_device_async");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // check if anything to do
  if (mp->nadj_rec_local == 0) return;

  // checks async-memcpy
  if (! GPU_ASYNC_COPY ) {
    exit_on_error("transfer_adj_to_device_async must be called with GPU_ASYNC_COPY == 1, \
please check mesh_constants_cuda.h");
  }

  // total number of receivers/adjoint sources
  int nrec = *h_nrec;

#ifdef USE_OPENCL
  if (run_opencl) {
    if (mp->has_last_copy_evt) {
      clCheck (clReleaseEvent (mp->last_copy_evt));
      mp->has_last_copy_evt = 0;
    }
    clCheck (clFinish (mocl.copy_queue));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    // waits for previous copy_stream call to be finished
    cudaStreamSynchronize(mp->copy_stream);
  }
#endif

  int irec,irec_local;

  irec_local = 0;
  for (irec = 0; irec < nrec; irec++) {
    if (mp->myrank == h_islice_selected_rec[irec]) {
      // takes only local sources
        mp->h_source_adjoint[INDEX2(NDIM,0,irec_local)]
              = h_source_adjoint[INDEX2(NDIM,0,irec_local)];

        mp->h_source_adjoint[INDEX2(NDIM,1,irec_local)]
              = h_source_adjoint[INDEX2(NDIM,1,irec_local)];

        mp->h_source_adjoint[INDEX2(NDIM,2,irec_local)]
              = h_source_adjoint[INDEX2(NDIM,2,irec_local)];

      // increases local receivers counter
      irec_local++;
    }
  }
  // check all local sources were added
  if (irec_local != mp->nadj_rec_local) exit_on_error("irec_local not equal to nadj_rec_local\n");

  // copies to GPU asynchronuously
#ifdef USE_OPENCL
  if (run_opencl) {
    cl_event *copy_evt = NULL;
    cl_uint num_evt = 0;

    // waits for previous compute_add_sources_adjoint_cuda_kernel() call to be finished
    clCheck (clFinish (mocl.command_queue));

    if (mp->has_last_copy_evt) {
      clCheck (clReleaseEvent (mp->last_copy_evt));
    }

    clCheck (clEnqueueWriteBuffer (mocl.copy_queue, mp->d_source_adjoint.ocl, CL_FALSE, 0,
                                   mp->nadj_rec_local * NDIM * sizeof (realw),
                                   mp->h_source_adjoint, num_evt, copy_evt, &mp->last_copy_evt));
    mp->has_last_copy_evt = 1;
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    // waits for previous compute_add_sources_adjoint_cuda_kernel() call to be finished
    cudaStreamSynchronize(mp->compute_stream);

    // copies extracted array values onto GPU
    // (asynchronous copy to GPU using copy_stream)
    cudaMemcpyAsync(mp->d_source_adjoint.cuda, mp->h_source_adjoint,(mp->nadj_rec_local)*NDIM*sizeof(realw),
                    cudaMemcpyHostToDevice,mp->copy_stream);
  }
#endif
}
