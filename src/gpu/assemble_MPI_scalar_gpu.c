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

// prepares a device array with with all inter-element edge-nodes -- this
// is followed by a memcpy and MPI operations

/*----------------------------------------------------------------------------------------------- */

// prepares and transfers the inter-element edge-nodes to the host to be MPI'd
// (elements on boundary)

extern EXTERN_LANG
void FC_FUNC_ (transfer_boun_pot_from_device,
               TRANSFER_BOUN_POT_FROM_DEVICE) (long *Mesh_pointer_f,
                                               realw *send_buffer,
                                               int *FORWARD_OR_ADJOINT) {

  TRACE ("transfer_boun_pot_from_device");
  int size_mpi_buffer;

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in transfer_boun_pot_from_device() routine");
  }

  // MPI buffer size
  size_mpi_buffer = mp->max_nibool_interfaces_oc * mp->num_interfaces_outer_core;

  // checks if anything to do
  if (size_mpi_buffer <= 0) return;

  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int) ceil ((double) mp->max_nibool_interfaces_oc / (double) blocksize)) * blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

  // sets gpu arrays
  gpu_realw_mem accel,buffer;
  float *h_buffer;

  if (*FORWARD_OR_ADJOINT == 1) {
    accel = mp->d_accel_outer_core;
    buffer = mp->d_send_accel_buffer_outer_core;
    h_buffer = mp->h_send_accel_buffer_oc;
  } else {
    // backward fields
    // debug
    DEBUG_BACKWARD_ASSEMBLY_OC ();
    accel = mp->d_b_accel_outer_core;
    buffer = mp->d_b_send_accel_buffer_outer_core;
    h_buffer = mp->h_b_send_accel_buffer_oc;
  }


#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;
    cl_event kernel_evt = NULL;

    clCheck (clSetKernelArg (mocl.kernels.prepare_boundary_potential_on_device, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.prepare_boundary_potential_on_device, idx++, sizeof (cl_mem), (void *) &buffer.ocl));
    clCheck (clSetKernelArg (mocl.kernels.prepare_boundary_potential_on_device, idx++, sizeof (int), (void *) &mp->num_interfaces_outer_core));
    clCheck (clSetKernelArg (mocl.kernels.prepare_boundary_potential_on_device, idx++, sizeof (int), (void *) &mp->max_nibool_interfaces_oc));
    clCheck (clSetKernelArg (mocl.kernels.prepare_boundary_potential_on_device, idx++, sizeof (cl_mem), (void *) &mp->d_nibool_interfaces_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.prepare_boundary_potential_on_device, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_interfaces_outer_core.ocl));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.prepare_boundary_potential_on_device, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, &kernel_evt));


    // copies buffer to CPU
    if (GPU_ASYNC_COPY) {
      // waits until kernel is finished before starting async memcpy
      clCheck (clFinish (mocl.command_queue));

      if (mp->has_last_copy_evt) {
        clCheck (clReleaseEvent (mp->last_copy_evt));
      }
      clCheck (clEnqueueReadBuffer (mocl.copy_queue, buffer.ocl, CL_FALSE, 0,
                                    size_mpi_buffer * sizeof (realw), h_buffer, 1, &kernel_evt, &mp->last_copy_evt));
      mp->has_last_copy_evt = 1;
    } else {
      // blocking copy
      clCheck (clEnqueueReadBuffer (mocl.command_queue, buffer.ocl, CL_TRUE, 0,
                                    size_mpi_buffer * sizeof (realw), send_buffer, 0, NULL, NULL));
    }

    clReleaseEvent (kernel_evt);
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    prepare_boundary_potential_on_device<<<grid,threads, 0, mp->compute_stream>>>(accel.cuda,
                                                                                  buffer.cuda,
                                                                                  mp->num_interfaces_outer_core,
                                                                                  mp->max_nibool_interfaces_oc,
                                                                                  mp->d_nibool_interfaces_outer_core.cuda,
                                                                                  mp->d_ibool_interfaces_outer_core.cuda);

    // copies buffer to CPU
    if (GPU_ASYNC_COPY) {
      // waits until kernel is finished before starting async memcpy
      cudaStreamSynchronize(mp->compute_stream);
      // copies buffer to CPU
      cudaMemcpyAsync(h_buffer,buffer.cuda,size_mpi_buffer*sizeof(realw),cudaMemcpyDeviceToHost,mp->copy_stream);
    } else {
      // synchronous copy
      print_CUDA_error_if_any(cudaMemcpy(send_buffer,buffer.cuda,size_mpi_buffer*sizeof(realw),
                                         cudaMemcpyDeviceToHost),98000);
    }
  }
#endif

  GPU_ERROR_CHECKING ("transfer_boun_pot_from_device");
}

/* ----------------------------------------------------------------------------------------------- */
// ASSEMBLY
/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (transfer_asmbl_pot_to_device,
               TRANSFER_ASMBL_POT_TO_DEVICE) (long *Mesh_pointer_f,
                                              realw *buffer_recv_scalar,
                                              int *FORWARD_OR_ADJOINT) {

  TRACE ("transfer_asmbl_pot_to_device");
  int size_mpi_buffer;

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in transfer_asmbl_pot_to_device() routine");
  }

  // buffer size
  size_mpi_buffer = (mp->max_nibool_interfaces_oc)*(mp->num_interfaces_outer_core);

  // checks if anything to do
  if (size_mpi_buffer <= 0 ) return;

  // assembles on GPU
  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int) ceil (((double) mp->max_nibool_interfaces_oc) / ((double) blocksize))) * blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (size_padded / blocksize, &num_blocks_x, &num_blocks_y);

  // sets gpu arrays
  gpu_realw_mem accel,buffer;

  if (*FORWARD_OR_ADJOINT == 1) {
    accel = mp->d_accel_outer_core;
    buffer = mp->d_send_accel_buffer_outer_core;
  } else {
    // backward fields
    // debug
    DEBUG_BACKWARD_ASSEMBLY_OC ();
    accel = mp->d_b_accel_outer_core;
    buffer = mp->d_b_send_accel_buffer_outer_core;
  }

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;
    cl_event *copy_evt = NULL;
    cl_uint num_evt = 0;

    if (GPU_ASYNC_COPY) {
      if (mp->has_last_copy_evt) {
        copy_evt = &mp->last_copy_evt;
        num_evt = 1;
      }
      // waits until previous copy finished
      clCheck (clFinish (mocl.copy_queue));
    }else{
      // copies scalar buffer onto GPU
      clCheck (clEnqueueWriteBuffer (mocl.command_queue, buffer.ocl, CL_TRUE, 0, size_mpi_buffer * sizeof (realw),
                                     buffer_recv_scalar, 0, NULL, NULL));
    }

    //assemble forward/backward field
    clCheck (clSetKernelArg (mocl.kernels.assemble_boundary_potential_on_device, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.assemble_boundary_potential_on_device, idx++, sizeof (cl_mem), (void *) &buffer.ocl));
    clCheck (clSetKernelArg (mocl.kernels.assemble_boundary_potential_on_device, idx++, sizeof (int), (void *) &mp->num_interfaces_outer_core));
    clCheck (clSetKernelArg (mocl.kernels.assemble_boundary_potential_on_device, idx++, sizeof (int), (void *) &mp->max_nibool_interfaces_oc));
    clCheck (clSetKernelArg (mocl.kernels.assemble_boundary_potential_on_device, idx++, sizeof (cl_mem), (void *) &mp->d_nibool_interfaces_outer_core.ocl));
    clCheck (clSetKernelArg (mocl.kernels.assemble_boundary_potential_on_device, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_interfaces_outer_core.ocl));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.assemble_boundary_potential_on_device, 2, NULL,
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
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // asynchronous copy
    if (GPU_ASYNC_COPY) {
      // Wait until previous copy stream finishes. We assemble while other compute kernels execute.
      cudaStreamSynchronize(mp->copy_stream);
    }else{
      // copies scalar buffer onto GPU
      print_CUDA_error_if_any(cudaMemcpy(buffer.cuda, buffer_recv_scalar,size_mpi_buffer*sizeof(realw),
                                         cudaMemcpyHostToDevice),99000);
    }

    //assemble field
    assemble_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(accel.cuda,
                                                                                 buffer.cuda,
                                                                                 mp->num_interfaces_outer_core,
                                                                                 mp->max_nibool_interfaces_oc,
                                                                                 mp->d_nibool_interfaces_outer_core.cuda,
                                                                                 mp->d_ibool_interfaces_outer_core.cuda);

  }
#endif

  //double end_time = get_time_val ();
  //printf ("Elapsed time: %e\n", end_time-start_time);
  GPU_ERROR_CHECKING ("transfer_asmbl_pot_to_device");
}
