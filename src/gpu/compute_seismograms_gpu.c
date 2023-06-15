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

extern EXTERN_LANG
void FC_FUNC_ (compute_seismograms_gpu,
               COMPUTE_SEISMOGRAMS_GPU) (long *Mesh_pointer,
                                         realw* seismograms,
                                         int* seismo_current_f,
                                         int* it_f,
                                         int* it_end_f,
                                         double* scale_displ_f,
                                         int* nlength_seismogram_f) {
  TRACE ("compute_seismograms_gpu");

  //get Mesh from Fortran integer wrapper
  Mesh *mp = (Mesh *) *Mesh_pointer;

  // checks if anything to do
  if (mp->nrec_local == 0) return;

  int seismo_current = *seismo_current_f - 1 ;  // -1 for Fortran -> C indexing

  int nlength_seismogram = *nlength_seismogram_f;
  int it = *it_f;
  int it_end = *it_end_f;

  realw scale_displ = (realw)(*scale_displ_f);

  int blocksize = NGLL3_PADDED;
  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->nrec_local, &num_blocks_x, &num_blocks_y);

  // seismogram array size
  int size = mp->nrec_local * nlength_seismogram;
  int size_buffer = NDIM * size * sizeof(realw);

  // determines wavefield depending on simulation type
  gpu_realw_mem displ;
  if (mp->simulation_type == 3) {
    // backward/reconstructed wavefield
    displ = mp->d_b_displ_crust_mantle;
  }else{
    // forward/adjoint wavefield
    displ = mp->d_displ_crust_mantle;
  }

  // note: due to subsampling, the last time step it == it_end might not be reached,
  //       but computing seismogram entries might end before.
  //       thus, both checks
  //         it%NTSTEP_BETWEEN_OUTPUT_SEISMOS == 0 || it == it_end
  //       might not be reached. instead we test if the seismogram array is full by
  //         seismo_current == nlength_seismogram - 1
  //       and copy it back whenever.
  //printf("debug: gpu seismo: seismo current/lenght %i/%i - it/it_end %i/%i\n",seismo_current,nlength_seismogram,it,it_end);

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_event kernel_evt = NULL;
    cl_uint idx = 0;

    //prepare field transfer array on device
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (int),    (void *) &mp->nrec_local));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (cl_mem), (void *) &displ.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_hxir.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_hetar.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_hgammar.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_seismograms.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_nu.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ispec_selected_rec.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_number_receiver_global.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (realw),  (void *) &scale_displ));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (int),    (void *) &seismo_current));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_seismograms_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, &kernel_evt));

    // copies array to CPU host
    if (seismo_current == nlength_seismogram - 1 || it == it_end){
      // intermediate copies could be made asynchronous,
      // only the last copy would be a synchronous copy to make sure we have all previous copies done.
      //
      // however, note that asynchronous copies in the copy stream would also require pinned host memory.
      // This is not use here yet for seismograms, thus the copy becomes blocking by default as well.
      // we'll leave the if-case here as a todo for future...
      //if (GPU_ASYNC_COPY && (it != it_end)) {
      //  // waits until kernel is finished before starting async memcpy
      //  clCheck (clFinish (mocl.command_queue));
      //
      //  if (mp->has_last_copy_evt) {
      //    clCheck (clReleaseEvent (mp->last_copy_evt));
      //  }
      //  clCheck (clEnqueueReadBuffer (mocl.copy_queue, mp->d_seismograms.ocl, CL_FALSE, 0,
      //                                size_buffer, seismograms, 1, &kernel_evt, &mp->last_copy_evt));
      //  mp->has_last_copy_evt = 1;
      //} else {
        // blocking copy
        clCheck (clEnqueueReadBuffer (mocl.command_queue, mp->d_seismograms.ocl, CL_TRUE, 0,
                                      size_buffer , seismograms, 0, NULL, NULL));
      //}
    }

    clReleaseEvent (kernel_evt);
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // prepare field transfer array on device
    compute_seismograms_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->nrec_local,
                                                                      displ.cuda,
                                                                      mp->d_ibool_crust_mantle.cuda,
                                                                      mp->d_hxir.cuda,mp->d_hetar.cuda,mp->d_hgammar.cuda,
                                                                      mp->d_seismograms.cuda,
                                                                      mp->d_nu.cuda,
                                                                      mp->d_ispec_selected_rec.cuda,
                                                                      mp->d_number_receiver_global.cuda,
                                                                      scale_displ,
                                                                      seismo_current);
    // copies array to CPU host
    if (seismo_current == nlength_seismogram - 1 || it == it_end){
      // intermediate copies could be made asynchronous,
      // only the last copy would be a synchronous copy to make sure we have all previous copies done.
      //
      // however, note that asynchronous copies in the copy stream would also require pinned host memory.
      // This is not use here yet for seismograms, thus the copy becomes blocking by default as well.
      // we'll leave the if-case here as a todo for future...
      //if (GPU_ASYNC_COPY &&  (it != it_end)) {
      //  // waits until kernel is finished before starting async memcpy
      //  cudaStreamSynchronize(mp->compute_stream);
      //  // copies buffer to CPU
      //  cudaMemcpyAsync(seismograms,mp->d_seismograms.cuda,size_buffer,cudaMemcpyDeviceToHost,mp->copy_stream);
      //} else {
        // synchronous copy
        print_CUDA_error_if_any(cudaMemcpy(seismograms,mp->d_seismograms.cuda,size_buffer,cudaMemcpyDeviceToHost),98000);
      //}
    }
  }
#endif
#ifdef USE_HIP
  if (run_hip) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // prepare field transfer array on device
    hipLaunchKernelGGL(HIP_KERNEL_NAME(compute_seismograms_kernel), grid, threads, 0, mp->compute_stream,
                                                                    mp->nrec_local,
                                                                    displ.hip,
                                                                    mp->d_ibool_crust_mantle.hip,
                                                                    mp->d_hxir.hip,mp->d_hetar.hip,mp->d_hgammar.hip,
                                                                    mp->d_seismograms.hip,
                                                                    mp->d_nu.hip,
                                                                    mp->d_ispec_selected_rec.hip,
                                                                    mp->d_number_receiver_global.hip,
                                                                    scale_displ,
                                                                    seismo_current);

    // copies array to CPU host
    if (seismo_current == nlength_seismogram - 1 || it == it_end){
      // intermediate copies could be made asynchronous,
      // only the last copy would be a synchronous copy to make sure we have all previous copies done.
      //
      // however, note that asynchronous copies in the copy stream would also require pinned host memory.
      // This is not use here yet for seismograms, thus the copy becomes blocking by default as well.
      // we'll leave the if-case here as a todo for future...
      //if (GPU_ASYNC_COPY &&  (it != it_end)) {
      //  // waits until kernel is finished before starting async memcpy
      //  hipStreamSynchronize(mp->compute_stream);
      //  // copies buffer to CPU
      //  hipMemcpyAsync(seismograms,mp->d_seismograms.hip,size_buffer,hipMemcpyDeviceToHost,mp->copy_stream);
      //} else {
        // synchronous copy
        print_HIP_error_if_any(hipMemcpy(seismograms,mp->d_seismograms.hip,size_buffer,hipMemcpyDeviceToHost),98000);
      //}
    }
  }
#endif

  GPU_ERROR_CHECKING ("compute_seismograms_gpu");
}
