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

#include <string.h>

#include "mesh_constants_gpu.h"

/* ----------------------------------------------------------------------------------------------- */
// MPI transfer
/* ----------------------------------------------------------------------------------------------- */

// prepares and transfers the inter-element edge-nodes to the host to be MPI'd
// (elements on boundary)

extern EXTERN_LANG
void FC_FUNC_(transfer_boun_from_device,
              TRANSFER_BOUN_FROM_DEVICE)(long *Mesh_pointer_f,
                                         realw *send_accel_buffer,
                                         int *IREGION,
                                         int *FORWARD_OR_ADJOINT) {
  TRACE("transfer_boun_from_device");

  int size_mpi_buffer = 0;

  int blocksize, size_padded;
  int num_blocks_x, num_blocks_y;

  // gpu arrays
  gpu_realw_mem accel,buffer;
  float *h_buffer;
  int num_interfaces;
  int max_nibool_interfaces;
  gpu_int_mem d_nibool_interfaces;
  gpu_int_mem d_ibool_interfaces;

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in transfer_boun_from_device() routine");
  }

  // sets up arrays for regions
  // crust/mantle region
  switch (*IREGION) {
  case IREGION_CRUST_MANTLE:
    size_mpi_buffer = NDIM*mp->max_nibool_interfaces_cm*mp->num_interfaces_crust_mantle;

    blocksize = BLOCKSIZE_TRANSFER;
    size_padded = ((int) ceil (((double) mp->max_nibool_interfaces_cm) / ((double) blocksize))) * blocksize;

    get_blocks_xy (size_padded / blocksize, &num_blocks_x, &num_blocks_y);

    // sets gpu arrays
    if (*FORWARD_OR_ADJOINT == 1) {
      accel = mp->d_accel_crust_mantle;
      buffer = mp->d_send_accel_buffer_crust_mantle;
      h_buffer = mp->h_send_accel_buffer_cm;
    } else {
      // backward fields
      // debug
      DEBUG_BACKWARD_ASSEMBLY_CM ();
      accel = mp->d_b_accel_crust_mantle;
      buffer = mp->d_b_send_accel_buffer_crust_mantle;
      h_buffer = mp->h_b_send_accel_buffer_cm;
    }
    num_interfaces = mp->num_interfaces_crust_mantle;
    max_nibool_interfaces = mp->max_nibool_interfaces_cm;
    d_nibool_interfaces = mp->d_nibool_interfaces_crust_mantle;
    d_ibool_interfaces = mp->d_ibool_interfaces_crust_mantle;
    break; // IREGION_CRUST_MANTLE

  // inner core region
  case IREGION_INNER_CORE:
    size_mpi_buffer = NDIM*mp->max_nibool_interfaces_ic*mp->num_interfaces_inner_core;

    blocksize = BLOCKSIZE_TRANSFER;
    size_padded = ((int) ceil (((double) mp->max_nibool_interfaces_ic) / ((double) blocksize))) * blocksize;

    get_blocks_xy (size_padded / blocksize, &num_blocks_x, &num_blocks_y);

    // sets gpu arrays
    if (*FORWARD_OR_ADJOINT == 1) {
      accel = mp->d_accel_inner_core;
      buffer = mp->d_send_accel_buffer_inner_core;
      h_buffer = mp->h_send_accel_buffer_ic;
    } else {
      // backward fields
      // debug
      DEBUG_BACKWARD_ASSEMBLY_IC ();
      accel = mp->d_b_accel_inner_core;
      buffer = mp->d_b_send_accel_buffer_inner_core;
      h_buffer = mp->h_b_send_accel_buffer_ic;
    }
    num_interfaces = mp->num_interfaces_inner_core;
    max_nibool_interfaces = mp->max_nibool_interfaces_ic;
    d_nibool_interfaces = mp->d_nibool_interfaces_inner_core;
    d_ibool_interfaces = mp->d_ibool_interfaces_inner_core;
    break; // IREGION_INNER_CORE

  default:
    exit_on_error("Error invalid IREGION in transfer_boun_from_device() routine");
  } // switch (*IREGION)

  // checks if anything to do
  if (size_mpi_buffer == 0 ) return;

  // fills buffer
#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;
    cl_event kernel_evt = NULL;

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clSetKernelArg (mocl.kernels.prepare_boundary_accel_on_device, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.prepare_boundary_accel_on_device, idx++, sizeof (cl_mem), (void *) &buffer.ocl));
    clCheck (clSetKernelArg (mocl.kernels.prepare_boundary_accel_on_device, idx++, sizeof (int), (void *) &num_interfaces));
    clCheck (clSetKernelArg (mocl.kernels.prepare_boundary_accel_on_device, idx++, sizeof (int), (void *) &max_nibool_interfaces));
    clCheck (clSetKernelArg (mocl.kernels.prepare_boundary_accel_on_device, idx++, sizeof (cl_mem), (void *) &d_nibool_interfaces.ocl));
    clCheck (clSetKernelArg (mocl.kernels.prepare_boundary_accel_on_device, idx++, sizeof (cl_mem), (void *) &d_ibool_interfaces.ocl));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.prepare_boundary_accel_on_device, 2, NULL,
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
                                    size_mpi_buffer * sizeof (realw), send_accel_buffer, 0, NULL, NULL));
    }
    clReleaseEvent (kernel_evt);
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    prepare_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(accel.cuda,
                                                                            buffer.cuda,
                                                                            num_interfaces,
                                                                            max_nibool_interfaces,
                                                                            d_nibool_interfaces.cuda,
                                                                            d_ibool_interfaces.cuda);

    // copies buffer to CPU
    if (GPU_ASYNC_COPY) {
      // waits until kernel is finished before starting async memcpy
      cudaStreamSynchronize(mp->compute_stream);
      // copies buffer to CPU
      cudaMemcpyAsync(h_buffer,buffer.cuda,size_mpi_buffer*sizeof(realw),cudaMemcpyDeviceToHost,mp->copy_stream);
    }else{
      // synchronous copy
      print_CUDA_error_if_any(cudaMemcpy(send_accel_buffer,buffer.cuda,size_mpi_buffer*sizeof(realw),
                                         cudaMemcpyDeviceToHost),41000);
    }
  }
#endif

  GPU_ERROR_CHECKING ("transfer_boun_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

// FORWARD_OR_ADJOINT == 1 for accel, and == 3 for b_accel

extern EXTERN_LANG
void FC_FUNC_ (transfer_asmbl_accel_to_device,
               TRANSFER_ASMBL_ACCEL_TO_DEVICE) (long *Mesh_pointer_f,
                                                realw *buffer_recv_vector,
                                                int *IREGION,
                                                int *FORWARD_OR_ADJOINT) {
  TRACE ("transfer_asmbl_accel_to_device");

  int size_mpi_buffer = 0;

  int blocksize, size_padded;
  int num_blocks_x, num_blocks_y;

  // gpu arrays
  gpu_realw_mem accel,buffer;
  int num_interfaces;
  int max_nibool_interfaces;
  gpu_int_mem d_nibool_interfaces;
  gpu_int_mem d_ibool_interfaces;

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in transfer_asmbl_accel_to_device() routine");
  }

  // sets up arrays for regions
  // crust/mantle region
  switch (*IREGION) {
  case IREGION_CRUST_MANTLE:
    size_mpi_buffer = NDIM*(mp->max_nibool_interfaces_cm)*(mp->num_interfaces_crust_mantle);

    // assembles values
    blocksize = BLOCKSIZE_TRANSFER;
    size_padded = ((int) ceil (((double) mp->max_nibool_interfaces_cm) / ((double) blocksize))) * blocksize;

    get_blocks_xy (size_padded / blocksize, &num_blocks_x, &num_blocks_y);

    // sets gpu arrays
    if (*FORWARD_OR_ADJOINT == 1) {
      accel = mp->d_accel_crust_mantle;
      buffer = mp->d_send_accel_buffer_crust_mantle;
    } else {
      // backward fields
      // debug
      DEBUG_BACKWARD_ASSEMBLY_CM ();
      accel = mp->d_b_accel_crust_mantle;
      buffer = mp->d_b_send_accel_buffer_crust_mantle;
    }
    num_interfaces = mp->num_interfaces_crust_mantle;
    max_nibool_interfaces = mp->max_nibool_interfaces_cm;
    d_nibool_interfaces = mp->d_nibool_interfaces_crust_mantle;
    d_ibool_interfaces = mp->d_ibool_interfaces_crust_mantle;
    break; // IREGION_CRUST_MANTLE

  // inner core region
  case IREGION_INNER_CORE:
    size_mpi_buffer = NDIM*(mp->max_nibool_interfaces_ic)*(mp->num_interfaces_inner_core);

    // assembles values
    blocksize = BLOCKSIZE_TRANSFER;
    size_padded = ((int) ceil (((double) mp->max_nibool_interfaces_ic) / ((double) blocksize))) * blocksize;

    get_blocks_xy (size_padded / blocksize, &num_blocks_x, &num_blocks_y);

    // sets gpu arrays
    if (*FORWARD_OR_ADJOINT == 1) {
      accel = mp->d_accel_inner_core;
      buffer = mp->d_send_accel_buffer_inner_core;
    } else {
      // backward fields
      // debug
      DEBUG_BACKWARD_ASSEMBLY_IC ();
      accel = mp->d_b_accel_inner_core;
      buffer = mp->d_b_send_accel_buffer_inner_core;
    }
    num_interfaces = mp->num_interfaces_inner_core;
    max_nibool_interfaces = mp->max_nibool_interfaces_ic;
    d_nibool_interfaces = mp->d_nibool_interfaces_inner_core;
    d_ibool_interfaces = mp->d_ibool_interfaces_inner_core;
    break; // IREGION_INNER_CORE

  default:
    exit_on_error("Error invalid IREGION in transfer_asmbl_accel_to_device() routine");
  } // switch (*IREGION)

  // checks if anything to do
  if (size_mpi_buffer == 0 ) return;

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;
    cl_event *copy_evt = NULL;
    cl_uint num_evt = 0;

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    // gets last copy event
    if (GPU_ASYNC_COPY) {
      if (mp->has_last_copy_evt) {
        copy_evt = &mp->last_copy_evt;
        num_evt = 1;
      }
      // Wait until previous copy stream finishes. We assemble while other compute kernels execute.
      clCheck (clFinish (mocl.copy_queue));
    } else {
      // copies vector buffer values to GPU
      clCheck (clEnqueueWriteBuffer (mocl.command_queue, buffer.ocl, CL_TRUE, 0,
                                     size_mpi_buffer * sizeof (realw), buffer_recv_vector, 0, NULL, NULL));
    }

    //assemble forward/backward accel
    clCheck (clSetKernelArg (mocl.kernels.assemble_boundary_accel_on_device, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.assemble_boundary_accel_on_device, idx++, sizeof (cl_mem), (void *) &buffer.ocl));
    clCheck (clSetKernelArg (mocl.kernels.assemble_boundary_accel_on_device, idx++, sizeof (int), (void *) &num_interfaces));
    clCheck (clSetKernelArg (mocl.kernels.assemble_boundary_accel_on_device, idx++, sizeof (int), (void *) &max_nibool_interfaces));
    clCheck (clSetKernelArg (mocl.kernels.assemble_boundary_accel_on_device, idx++, sizeof (cl_mem), (void *) &d_nibool_interfaces.ocl));
    clCheck (clSetKernelArg (mocl.kernels.assemble_boundary_accel_on_device, idx++, sizeof (cl_mem), (void *) &d_ibool_interfaces.ocl));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.assemble_boundary_accel_on_device, 2, NULL,
                                     global_work_size, local_work_size, num_evt, copy_evt, NULL));

    // releases copy event
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
      // (cudaMemcpy implicitly synchronizes all other cuda operations)
      // copies vector buffer values to GPU
      print_CUDA_error_if_any(cudaMemcpy(buffer.cuda,buffer_recv_vector,size_mpi_buffer*sizeof(realw),
                                         cudaMemcpyHostToDevice),41000);
    }

    //assemble forward/backward accel
    assemble_boundary_accel_on_device<<<grid,threads,0,mp->compute_stream>>>(accel.cuda,
                                                                             buffer.cuda,
                                                                             num_interfaces,
                                                                             max_nibool_interfaces,
                                                                             d_nibool_interfaces.cuda,
                                                                             d_ibool_interfaces.cuda);
  }
#endif

  GPU_ERROR_CHECKING ("transfer_asmbl_accel_to_device");
}


/* ----------------------------------------------------------------------------------------------- */
// Asynchronous memory copy for MPI buffers
/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_buffer_to_device_async,
              TRANSFER_BUFFER_TO_DEVICE_ASYNC)(long* Mesh_pointer_f,
                                               realw* buffer_f,
                                               int* IREGION,
                                               int* FORWARD_OR_ADJOINT) {
  // asynchronous transfer from host to device
  TRACE("transfer_buffer_to_device_async");

  int size_mpi_buffer = 0;

  // gpu arrays
  gpu_realw_mem buffer;
  float *h_buffer;

  // get Mesh from Fortran integer wrapper
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // checks async-memcpy
  if (! GPU_ASYNC_COPY ) {
    exit_on_error("transfer_buffer_to_device_async must be called with GPU_ASYNC_COPY == 1, please check mesh_constants_cuda.h");
  }

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in transfer_buffer_to_device_async() routine");
  }

  // sets up arrays for regions
  switch (*IREGION) {
  case IREGION_CRUST_MANTLE:
    // crust/mantle region
    size_mpi_buffer = NDIM * mp->max_nibool_interfaces_cm * mp->num_interfaces_crust_mantle;
    // sets gpu arrays
    if (*FORWARD_OR_ADJOINT == 1) {
      buffer = mp->d_send_accel_buffer_crust_mantle;
      h_buffer = mp->h_recv_accel_buffer_cm;
    } else {
      // backward fields
      // debug
      DEBUG_BACKWARD_ASSEMBLY_CM ();
      buffer = mp->d_b_send_accel_buffer_crust_mantle;
      h_buffer = mp->h_b_recv_accel_buffer_cm;
    }
    break;

  case IREGION_INNER_CORE:
    // inner core region
    size_mpi_buffer = NDIM * mp->max_nibool_interfaces_ic * mp->num_interfaces_inner_core;
    // sets gpu arrays
    if (*FORWARD_OR_ADJOINT == 1) {
      buffer = mp->d_send_accel_buffer_inner_core;
      h_buffer = mp->h_recv_accel_buffer_ic;
    } else {
      // backward fields
      // debug
      DEBUG_BACKWARD_ASSEMBLY_IC ();
      buffer = mp->d_b_send_accel_buffer_inner_core;
      h_buffer = mp->h_b_recv_accel_buffer_ic;
    }
    break;

  case IREGION_OUTER_CORE:
    // outer core region
    size_mpi_buffer = mp->max_nibool_interfaces_oc * mp->num_interfaces_outer_core;
    // sets gpu arrays
    if (*FORWARD_OR_ADJOINT == 1) {
      buffer = mp->d_send_accel_buffer_outer_core;
      h_buffer = mp->h_recv_accel_buffer_oc;
    } else {
      // backward fields
      // debug
      DEBUG_BACKWARD_ASSEMBLY_IC ();
      buffer = mp->d_b_send_accel_buffer_outer_core;
      h_buffer = mp->h_b_recv_accel_buffer_oc;
    }
    break;

  default:
    exit_on_error("Error invalid IREGION in transfer_buffer_to_device_async() routine");
  } // switch (*IREGION)

  // checks if anything to do
  if (size_mpi_buffer == 0) return;

  // copy on host memory
  memcpy(h_buffer, buffer_f, size_mpi_buffer * sizeof(realw));

  // asynchronous copy to GPU using copy_stream
#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck(clEnqueueWriteBuffer(mocl.copy_queue, buffer.ocl, CL_FALSE, 0,
                                 size_mpi_buffer * sizeof (realw), h_buffer, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    cudaMemcpyAsync(buffer.cuda, h_buffer,size_mpi_buffer*sizeof(realw),cudaMemcpyHostToDevice,mp->copy_stream);
  }
#endif

  GPU_ERROR_CHECKING ("transfer_buffer_to_device_async");
}


/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(sync_copy_from_device,
              SYNC_copy_FROM_DEVICE)(long* Mesh_pointer_f,
                                     int* iphase,
                                     realw* send_buffer,
                                     int* IREGION,
                                     int* FORWARD_OR_ADJOINT) {
  // synchronizes copy stream before copying buffers from pinned memory to CPU host
  TRACE("sync_copy_from_device");

  int size_mpi_buffer = 0;

  // gpu arrays
  float *h_buffer;

  // get Mesh from Fortran integer wrapper
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // checks async-memcpy
  if (! GPU_ASYNC_COPY ) {
    exit_on_error("sync_copy_from_device must be called with GPU_ASYNC_COPY == 1, please check mesh_constants_gpu.h");
  }

  // Wait until async-memcpy of outer elements is finished and start MPI.
  if (*iphase != 2) {
    exit_on_error("sync_copy_from_device must be called for iphase == 2");
  }

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in sync_copy_from_device() routine");
  }

  // sets up arrays for regions
  switch (*IREGION) {
  case IREGION_CRUST_MANTLE:
    // crust/mantle
    size_mpi_buffer = NDIM * mp->max_nibool_interfaces_cm * mp->num_interfaces_crust_mantle;
    // sets host arrays
    if (*FORWARD_OR_ADJOINT == 1) {
      h_buffer = mp->h_send_accel_buffer_cm;
    } else {
      // backward fields
      h_buffer = mp->h_b_send_accel_buffer_cm;
    }
    break;

  case IREGION_INNER_CORE:
    // inner core
    size_mpi_buffer = NDIM * mp->max_nibool_interfaces_ic * mp->num_interfaces_inner_core;
    // sets host arrays
    if (*FORWARD_OR_ADJOINT == 1) {
      h_buffer = mp->h_send_accel_buffer_ic;
    } else {
      // backward fields
      h_buffer = mp->h_b_send_accel_buffer_ic;
    }
    break;

  case IREGION_OUTER_CORE:
    // outer core
    size_mpi_buffer = mp->max_nibool_interfaces_oc * mp->num_interfaces_outer_core;
    // sets host arrays
    if (*FORWARD_OR_ADJOINT == 1) {
      h_buffer = mp->h_send_accel_buffer_oc;
    } else {
      // backward fields
      h_buffer = mp->h_b_send_accel_buffer_oc;
    }
    break;

  default:
    exit_on_error("Error invalid IREGION in sync_copy_from_device() routine");
  } // switch (*IREGION)

  // checks if anything to do
  if (size_mpi_buffer == 0) return;

  // waits for asynchronous copy to finish
#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clFinish (mocl.copy_queue));
    if (mp->has_last_copy_evt) {
      clCheck (clReleaseEvent (mp->last_copy_evt));
      mp->has_last_copy_evt = 0;
    }
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    cudaStreamSynchronize(mp->copy_stream);
  }
#endif

  // There have been problems using the pinned-memory with MPI, so
  // we copy the buffer into a non-pinned region.
  memcpy(send_buffer, h_buffer, size_mpi_buffer * sizeof(realw));

  // memory copy is now finished, so non-blocking MPI send can proceed

  GPU_ERROR_CHECKING ("sync_copy_from_device");
}

