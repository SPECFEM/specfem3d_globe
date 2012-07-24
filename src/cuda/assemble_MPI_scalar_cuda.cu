/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  2 . 0
 !               ---------------------------------------
 !
 !          Main authors: Dimitri Komatitsch and Jeroen Tromp
 !    Princeton University, USA and University of Pau / CNRS / INRIA
 ! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
 !                            April 2011
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
  int iinterface=0;

  for( iinterface=0; iinterface < num_interfaces; iinterface++) {
    if(id<d_nibool_interfaces[iinterface]) {
      d_send_potential_dot_dot_buffer[(id + max_nibool_interfaces*iinterface)] =
      d_potential_dot_dot_acoustic[(d_ibool_interfaces[id+max_nibool_interfaces*iinterface]-1)];
    }
  }

}


/* ----------------------------------------------------------------------------------------------- */

// prepares and transfers the inter-element edge-nodes to the host to be MPI'd
extern "C"
void FC_FUNC_(transfer_boun_pot_from_device,
              TRANSFER_BOUN_POT_FROM_DEVICE)(long* Mesh_pointer_f,
                                             realw* send_potential_dot_dot_buffer,
                                             int* FORWARD_OR_ADJOINT){

  TRACE("transfer_boun_pot_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if( mp->num_interfaces_outer_core == 0 ) return;

  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int)ceil(((double)(mp->max_nibool_interfaces_outer_core))/((double)blocksize)))*blocksize;
  int num_blocks_x = size_padded/blocksize;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  if(*FORWARD_OR_ADJOINT == 1) {
    prepare_boundary_potential_on_device<<<grid,threads>>>(mp->d_accel_outer_core,
                                                           mp->d_send_accel_buffer_outer_core,
                                                           mp->num_interfaces_outer_core,
                                                           mp->max_nibool_interfaces_outer_core,
                                                           mp->d_nibool_interfaces_outer_core,
                                                           mp->d_ibool_interfaces_outer_core);
  }
  else if(*FORWARD_OR_ADJOINT == 3) {
    prepare_boundary_potential_on_device<<<grid,threads>>>(mp->d_b_accel_outer_core,
                                                           mp->d_send_accel_buffer_outer_core,
                                                           mp->num_interfaces_outer_core,
                                                           mp->max_nibool_interfaces_outer_core,
                                                           mp->d_nibool_interfaces_outer_core,
                                                           mp->d_ibool_interfaces_outer_core);
  }

  print_CUDA_error_if_any(cudaMemcpy(send_potential_dot_dot_buffer,mp->d_send_accel_buffer_outer_core,
                                     (mp->max_nibool_interfaces_outer_core)*(mp->num_interfaces_outer_core)*sizeof(realw),
                                     cudaMemcpyDeviceToHost),98000);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_boundary_kernel");
#endif
}


/* ----------------------------------------------------------------------------------------------- */


__global__ void assemble_boundary_potential_on_device(realw* d_potential_dot_dot_acoustic,
                                                      realw* d_send_potential_dot_dot_buffer,
                                                      int num_interfaces,
                                                      int max_nibool_interfaces,
                                                      int* d_nibool_interfaces,
                                                      int* d_ibool_interfaces) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  int iinterface=0;

  for( iinterface=0; iinterface < num_interfaces; iinterface++) {
    if(id<d_nibool_interfaces[iinterface]) {

      atomicAdd(&d_potential_dot_dot_acoustic[(d_ibool_interfaces[id+max_nibool_interfaces*iinterface]-1)],
                d_send_potential_dot_dot_buffer[(id + max_nibool_interfaces*iinterface)]);
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

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  //double start_time = get_time();

  // checks if anything to do
  if( mp->num_interfaces_outer_core == 0 ) return;

  // copies scalar buffer onto GPU
  cudaMemcpy(mp->d_send_accel_buffer_outer_core, buffer_recv_scalar,
             (mp->max_nibool_interfaces_outer_core)*(mp->num_interfaces_outer_core)*sizeof(realw),
             cudaMemcpyHostToDevice);

  // assembles on GPU
  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_outer_core)/((double)blocksize)))*blocksize;
  int num_blocks_x = size_padded/blocksize;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  if(*FORWARD_OR_ADJOINT == 1) {
    //assemble forward field
    assemble_boundary_potential_on_device<<<grid,threads>>>(mp->d_accel_outer_core,
                                                            mp->d_send_accel_buffer_outer_core,
                                                            mp->num_interfaces_outer_core,
                                                            mp->max_nibool_interfaces_outer_core,
                                                            mp->d_nibool_interfaces_outer_core,
                                                            mp->d_ibool_interfaces_outer_core);
  }
  else if(*FORWARD_OR_ADJOINT == 3) {
    //assemble reconstructed/backward field
    assemble_boundary_potential_on_device<<<grid,threads>>>(mp->d_b_accel_outer_core,
                                                            mp->d_send_accel_buffer_outer_core,
                                                            mp->num_interfaces_outer_core,
                                                            mp->max_nibool_interfaces_outer_core,
                                                            mp->d_nibool_interfaces_outer_core,
                                                            mp->d_ibool_interfaces_outer_core);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("transfer_asmbl_pot_to_device");
#endif
}

