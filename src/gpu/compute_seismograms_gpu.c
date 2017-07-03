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

extern EXTERN_LANG
void FC_FUNC_ (compute_seismograms_gpu,
               COMPUTE_SEISMOGRAMS_GPU) (long *Mesh_pointer_f,
                                         realw* seismograms,
                                         int* seismo_currentf,int* itf, double* scale_displf,int* NTSTEP_BETWEEN_OUTPUT_SEISMOSf,int* NSTEPf) {
  TRACE ("compute_seismograms_gpu");

  //get Mesh from Fortran integer wrapper
  Mesh *mp = (Mesh *) *Mesh_pointer_f;
  // checks if anything to do
  if (mp->nrec_local == 0) return;

  int seismo_current = *seismo_currentf - 1 ;
  int NTSTEP_BETWEEN_OUTPUT_SEISMOS = *NTSTEP_BETWEEN_OUTPUT_SEISMOSf;
  int NSTEP = *NSTEPf;
  int it = *itf;
  realw scale_displ = (realw)(*scale_displf);

  int blocksize = NGLL3_PADDED;
  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->nrec_local, &num_blocks_x, &num_blocks_y);

  int size_buffer = 3*mp->nrec_local* sizeof (realw);

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_event kernel_evt = NULL;
    cl_uint idx = 0;

    //prepare field transfer array on device
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (int),    (void *) &mp->nrec_local));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_displ_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_xir.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_etar.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_gammar.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_seismograms.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_nu.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ispec_selected_rec.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_number_receiver_global.ocl));
    clCheck (clSetKernelArg (mocl.kernels.compute_seismograms_kernel, idx++, sizeof (int),    (void *) &scale_displ));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.compute_seismograms_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, &kernel_evt));

    // copies buffer to CPU
    if (GPU_ASYNC_COPY &&  ((seismo_current+1) != NTSTEP_BETWEEN_OUTPUT_SEISMOS) && (it != NSTEP) ) {
      // waits until kernel is finished before starting async memcpy
      clCheck (clFinish (mocl.command_queue));

      if (mp->has_last_copy_evt) {
        clCheck (clReleaseEvent (mp->last_copy_evt));
      }
      clCheck (clEnqueueReadBuffer (mocl.copy_queue, mp->d_seismograms.ocl, CL_FALSE, 0,
                                    size_buffer, seismograms + 3*mp->nrec_local*seismo_current, 1, &kernel_evt, &mp->last_copy_evt));
      mp->has_last_copy_evt = 1;
    } else {
      // blocking copy
      clCheck (clEnqueueReadBuffer (mocl.command_queue, mp->d_seismograms.ocl, CL_TRUE, 0,
                                    size_buffer , seismograms + 3*mp->nrec_local*seismo_current, 0, NULL, NULL));
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
                                                                      mp->d_displ_crust_mantle.cuda,
                                                                      mp->d_ibool_crust_mantle.cuda,
                                                                      mp->d_xir.cuda,mp->d_etar.cuda,mp->d_gammar.cuda,
                                                                      mp->d_seismograms.cuda,
                                                                      mp->d_nu.cuda,
                                                                      mp->d_ispec_selected_rec.cuda,
                                                                      mp->d_number_receiver_global.cuda,
                                                                      scale_displ);
    // copies buffer to CPU
    if (GPU_ASYNC_COPY &&  ((seismo_current+1) != NTSTEP_BETWEEN_OUTPUT_SEISMOS) && (it != NSTEP) ) {
      // waits until kernel is finished before starting async memcpy
      cudaStreamSynchronize(mp->compute_stream);
      // copies buffer to CPU
      cudaMemcpyAsync(seismograms + 3*mp->nrec_local*seismo_current,mp->d_seismograms.cuda,size_buffer,cudaMemcpyDeviceToHost,mp->copy_stream);
    } else {
      // synchronous copy
      print_CUDA_error_if_any(cudaMemcpy(seismograms + 3*mp->nrec_local*seismo_current,mp->d_seismograms.cuda,size_buffer,
                                         cudaMemcpyDeviceToHost),98000);
    }

  }
  
#endif


  GPU_ERROR_CHECKING ("compute_seismograms_gpu");
}
