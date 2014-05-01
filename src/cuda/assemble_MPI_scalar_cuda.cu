/*
 !=====================================================================
 !
 !          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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

#include <stdio.h>
#include <cuda.h>
#include <cublas.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"

/* ----------------------------------------------------------------------------------------------- */

// prepares a device array with with all inter-element edge-nodes -- this
// is followed by a memcpy and MPI operations
__global__ void prepare_boundary_potential_on_device(realw* d_potential_dot_dot_acoustic,
                                                     realw* d_send_potential_dot_dot_buffer,
                                                     int num_interfaces,
                                                     int max_nibool_interfaces,
                                                     int* d_nibool_interfaces,
                                                     int* d_ibool_interfaces) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  int iglob,iloc;

  for( int iinterface=0; iinterface < num_interfaces; iinterface++) {
    if(id<d_nibool_interfaces[iinterface]) {

      iloc = id + max_nibool_interfaces*iinterface;
      iglob = d_ibool_interfaces[iloc] - 1;

      // fills buffer
      d_send_potential_dot_dot_buffer[iloc] = d_potential_dot_dot_acoustic[iglob];
    }
  }

}

/* ----------------------------------------------------------------------------------------------- */

// MPI transfer

/* ----------------------------------------------------------------------------------------------- */

// prepares and transfers the inter-element edge-nodes to the host to be MPI'd
// (elements on boundary)

extern "C"
void FC_FUNC_(transfer_boun_pot_from_device,
              TRANSFER_BOUN_POT_FROM_DEVICE)(long* Mesh_pointer_f,
                                             realw* send_buffer,
                                             int* FORWARD_OR_ADJOINT){

  TRACE("transfer_boun_pot_from_device");
  int size_mpi_buffer;

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  // MPI buffer size
  size_mpi_buffer = (mp->max_nibool_interfaces_oc)*(mp->num_interfaces_outer_core);

  // checks if anything to do
  if( size_mpi_buffer <= 0 ) return;

  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int)ceil(((double)(mp->max_nibool_interfaces_oc))/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  if(*FORWARD_OR_ADJOINT == 1) {
    prepare_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel_outer_core,
                                                                                mp->d_send_accel_buffer_outer_core,
                                                                                mp->num_interfaces_outer_core,
                                                                                mp->max_nibool_interfaces_oc,
                                                                                mp->d_nibool_interfaces_outer_core,
                                                                                mp->d_ibool_interfaces_outer_core);

    // copies buffer to CPU
    if( GPU_ASYNC_COPY ){
      // waits until kernel is finished before starting async memcpy
      cudaStreamSynchronize(mp->compute_stream);
      // copies buffer to CPU
      cudaMemcpyAsync(mp->h_send_accel_buffer_oc,mp->d_send_accel_buffer_outer_core,size_mpi_buffer*sizeof(realw),
                      cudaMemcpyDeviceToHost,mp->copy_stream);
    }else{
      // synchronous copy
      print_CUDA_error_if_any(cudaMemcpy(send_buffer,mp->d_send_accel_buffer_outer_core,
                                       size_mpi_buffer*sizeof(realw),
                                       cudaMemcpyDeviceToHost),98000);
    }

  }
  else if(*FORWARD_OR_ADJOINT == 3) {
    // debug
    DEBUG_BACKWARD_ASSEMBLY();

    prepare_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_accel_outer_core,
                                                                                mp->d_b_send_accel_buffer_outer_core,
                                                                                mp->num_interfaces_outer_core,
                                                                                mp->max_nibool_interfaces_oc,
                                                                                mp->d_nibool_interfaces_outer_core,
                                                                                mp->d_ibool_interfaces_outer_core);


    // copies buffer to CPU
    if( GPU_ASYNC_COPY ){
      // waits until kernel is finished before starting async memcpy
      cudaStreamSynchronize(mp->compute_stream);
      // copies buffer to CPU
      cudaMemcpyAsync(mp->h_b_send_accel_buffer_oc,mp->d_b_send_accel_buffer_outer_core,size_mpi_buffer*sizeof(realw),
                      cudaMemcpyDeviceToHost,mp->copy_stream);
    }else{
      // synchronous copy
      print_CUDA_error_if_any(cudaMemcpy(send_buffer,mp->d_b_send_accel_buffer_outer_core,
                                       size_mpi_buffer*sizeof(realw),
                                       cudaMemcpyDeviceToHost),98001);
    }

  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("transfer_boun_pot_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// ASSEMBLY

/* ----------------------------------------------------------------------------------------------- */


__global__ void assemble_boundary_potential_on_device(realw* d_potential_dot_dot_acoustic,
                                                      realw* d_send_potential_dot_dot_buffer,
                                                      int num_interfaces,
                                                      int max_nibool_interfaces,
                                                      int* d_nibool_interfaces,
                                                      int* d_ibool_interfaces) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  int iglob,iloc;

  for( int iinterface=0; iinterface < num_interfaces; iinterface++) {
    if(id<d_nibool_interfaces[iinterface]) {

      iloc = id + max_nibool_interfaces*iinterface;

      iglob = d_ibool_interfaces[iloc]-1;

      // assembles values
      atomicAdd(&d_potential_dot_dot_acoustic[iglob],d_send_potential_dot_dot_buffer[iloc]);
    }
  }
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_asmbl_pot_to_device,
              TRANSFER_ASMBL_POT_TO_DEVICE)(long* Mesh_pointer,
                                            realw* buffer_recv_scalar,
                                            int* FORWARD_OR_ADJOINT) {

  TRACE("transfer_asmbl_pot_to_device");
  int size_mpi_buffer;

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // buffer size
  size_mpi_buffer = (mp->max_nibool_interfaces_oc)*(mp->num_interfaces_outer_core);

  // checks if anything to do
  if( size_mpi_buffer <= 0 ) return;

  // assembles on GPU
  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_oc)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  if(*FORWARD_OR_ADJOINT == 1) {

    // asynchronous copy
    if( GPU_ASYNC_COPY ){
      // Wait until previous copy stream finishes. We assemble while other compute kernels execute.
      cudaStreamSynchronize(mp->copy_stream);
    }else{
      // copies scalar buffer onto GPU
      print_CUDA_error_if_any(cudaMemcpy(mp->d_send_accel_buffer_outer_core, buffer_recv_scalar,size_mpi_buffer*sizeof(realw),
                                         cudaMemcpyHostToDevice),99000);
    }

    //assembles forward field
    assemble_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel_outer_core,
                                                                                 mp->d_send_accel_buffer_outer_core,
                                                                                 mp->num_interfaces_outer_core,
                                                                                 mp->max_nibool_interfaces_oc,
                                                                                 mp->d_nibool_interfaces_outer_core,
                                                                                 mp->d_ibool_interfaces_outer_core);
  }else if(*FORWARD_OR_ADJOINT == 3) {
    // debug
    DEBUG_BACKWARD_ASSEMBLY();

    // asynchronous copy
    if( GPU_ASYNC_COPY ){
      // Wait until previous copy stream finishes. We assemble while other compute kernels execute.
      cudaStreamSynchronize(mp->copy_stream);
    }else{
      // (cudaMemcpy implicitly synchronizes all other cuda operations)
      // copies scalar buffer onto GPU
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_send_accel_buffer_outer_core, buffer_recv_scalar,size_mpi_buffer*sizeof(realw),
                                         cudaMemcpyHostToDevice),99001);
    }

    //assembles reconstructed/backward field
    assemble_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_accel_outer_core,
                                                                                 mp->d_b_send_accel_buffer_outer_core,
                                                                                 mp->num_interfaces_outer_core,
                                                                                 mp->max_nibool_interfaces_oc,
                                                                                 mp->d_nibool_interfaces_outer_core,
                                                                                 mp->d_ibool_interfaces_outer_core);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("transfer_asmbl_pot_to_device");
#endif
}

