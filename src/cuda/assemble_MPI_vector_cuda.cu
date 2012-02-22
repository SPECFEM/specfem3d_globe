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

// prepares a device array with all inter-element edge-nodes -- this
// is followed by a memcpy and MPI operations
__global__ void prepare_boundary_accel_on_device(realw* d_accel, realw* d_send_accel_buffer,
                                                 int num_interfaces_ext_mesh,
                                                 int max_nibool_interfaces_ext_mesh,
                                                 int* d_nibool_interfaces_ext_mesh,
                                                 int* d_ibool_interfaces_ext_mesh) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  int iinterface=0;

  for( iinterface=0; iinterface < num_interfaces_ext_mesh; iinterface++) {
    if(id<d_nibool_interfaces_ext_mesh[iinterface]) {
      d_send_accel_buffer[3*(id + max_nibool_interfaces_ext_mesh*iinterface)] =
      d_accel[3*(d_ibool_interfaces_ext_mesh[id+max_nibool_interfaces_ext_mesh*iinterface]-1)];
      d_send_accel_buffer[3*(id + max_nibool_interfaces_ext_mesh*iinterface)+1] =
      d_accel[3*(d_ibool_interfaces_ext_mesh[id+max_nibool_interfaces_ext_mesh*iinterface]-1)+1];
      d_send_accel_buffer[3*(id + max_nibool_interfaces_ext_mesh*iinterface)+2] =
      d_accel[3*(d_ibool_interfaces_ext_mesh[id+max_nibool_interfaces_ext_mesh*iinterface]-1)+2];
    }
  }

}

/* ----------------------------------------------------------------------------------------------- */

// prepares and transfers the inter-element edge-nodes to the host to be MPI'd
// (elements on boundary)
extern "C"
void FC_FUNC_(transfer_boun_accel_from_device,
              TRANSFER_BOUN_ACCEL_FROM_DEVICE)(long* Mesh_pointer_f,
                                                  realw* send_accel_buffer,
                                                  int* IREGION,
                                                  int* FORWARD_OR_ADJOINT){
  TRACE("transfer_boun_accel_from_device");
  int blocksize,size_padded,num_blocks_x,num_blocks_y;

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  // crust/mantle region
  if( *IREGION == IREGION_CRUST_MANTLE ){
    if( mp->num_interfaces_crust_mantle > 0 ){

      blocksize = BLOCKSIZE_TRANSFER;
      size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_crust_mantle)/((double)blocksize)))*blocksize;
      num_blocks_x = size_padded/blocksize;
      num_blocks_y = 1;
      while(num_blocks_x > 65535) {
        num_blocks_x = (int) ceil(num_blocks_x*0.5f);
        num_blocks_y = num_blocks_y*2;
      }
      dim3 grid(num_blocks_x,num_blocks_y);
      dim3 threads(blocksize,1,1);

      if(*FORWARD_OR_ADJOINT == 1) {
        prepare_boundary_accel_on_device<<<grid,threads>>>(mp->d_accel_crust_mantle,
                                                           mp->d_send_accel_buffer_crust_mantle,
                                                           mp->num_interfaces_crust_mantle,
                                                           mp->max_nibool_interfaces_crust_mantle,
                                                           mp->d_nibool_interfaces_crust_mantle,
                                                           mp->d_ibool_interfaces_crust_mantle);
      }
      else if(*FORWARD_OR_ADJOINT == 3) {
        prepare_boundary_accel_on_device<<<grid,threads>>>(mp->d_b_accel_crust_mantle,
                                                           mp->d_send_accel_buffer_crust_mantle,
                                                           mp->num_interfaces_crust_mantle,
                                                           mp->max_nibool_interfaces_crust_mantle,
                                                           mp->d_nibool_interfaces_crust_mantle,
                                                           mp->d_ibool_interfaces_crust_mantle);
      }

      // copies buffer to CPU
      cudaMemcpy(send_accel_buffer,mp->d_send_accel_buffer_crust_mantle,
                 3*mp->max_nibool_interfaces_crust_mantle*mp->num_interfaces_crust_mantle*sizeof(realw),
                 cudaMemcpyDeviceToHost);

    }
  }

  // inner core region
  if( *IREGION == IREGION_INNER_CORE ){
    if( mp->num_interfaces_inner_core > 0 ){

      blocksize = BLOCKSIZE_TRANSFER;
      size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_inner_core)/((double)blocksize)))*blocksize;
      num_blocks_x = size_padded/blocksize;
      num_blocks_y = 1;
      while(num_blocks_x > 65535) {
        num_blocks_x = (int) ceil(num_blocks_x*0.5f);
        num_blocks_y = num_blocks_y*2;
      }
      dim3 grid(num_blocks_x,num_blocks_y);
      dim3 threads(blocksize,1,1);

      if(*FORWARD_OR_ADJOINT == 1) {
        prepare_boundary_accel_on_device<<<grid,threads>>>(mp->d_accel_inner_core,
                                                           mp->d_send_accel_buffer_inner_core,
                                                           mp->num_interfaces_inner_core,
                                                           mp->max_nibool_interfaces_inner_core,
                                                           mp->d_nibool_interfaces_inner_core,
                                                           mp->d_ibool_interfaces_inner_core);
      }
      else if(*FORWARD_OR_ADJOINT == 3) {
        prepare_boundary_accel_on_device<<<grid,threads>>>(mp->d_b_accel_inner_core,
                                                           mp->d_send_accel_buffer_inner_core,
                                                           mp->num_interfaces_inner_core,
                                                           mp->max_nibool_interfaces_inner_core,
                                                           mp->d_nibool_interfaces_inner_core,
                                                           mp->d_ibool_interfaces_inner_core);
      }

      // copies buffer to CPU
      cudaMemcpy(send_accel_buffer,mp->d_send_accel_buffer_inner_core,
                 3*mp->max_nibool_interfaces_inner_core*mp->num_interfaces_inner_core*sizeof(realw),
                 cudaMemcpyDeviceToHost);

    }
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("transfer_boun_accel_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void assemble_boundary_accel_on_device(realw* d_accel, realw* d_send_accel_buffer,
                                                  int num_interfaces_ext_mesh,
                                                  int max_nibool_interfaces_ext_mesh,
                                                  int* d_nibool_interfaces_ext_mesh,
                                                  int* d_ibool_interfaces_ext_mesh) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  int iinterface=0;

  for( iinterface=0; iinterface < num_interfaces_ext_mesh; iinterface++) {
    if(id < d_nibool_interfaces_ext_mesh[iinterface]) {
      atomicAdd(&d_accel[3*(d_ibool_interfaces_ext_mesh[id+max_nibool_interfaces_ext_mesh*iinterface]-1)],
                d_send_accel_buffer[3*(id + max_nibool_interfaces_ext_mesh*iinterface)]);
      atomicAdd(&d_accel[3*(d_ibool_interfaces_ext_mesh[id+max_nibool_interfaces_ext_mesh*iinterface]-1)+1],
                d_send_accel_buffer[3*(id + max_nibool_interfaces_ext_mesh*iinterface)+1]);
      atomicAdd(&d_accel[3*(d_ibool_interfaces_ext_mesh[id+max_nibool_interfaces_ext_mesh*iinterface]-1)+2],
                d_send_accel_buffer[3*(id + max_nibool_interfaces_ext_mesh*iinterface)+2]);
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

// FORWARD_OR_ADJOINT == 1 for accel, and == 3 for b_accel
extern "C"
void FC_FUNC_(transfer_asmbl_accel_to_device,
              TRANSFER_ASMBL_ACCEL_TO_DEVICE)(long* Mesh_pointer,
                                              realw* buffer_recv_vector,
                                              int* IREGION,
                                              int* FORWARD_OR_ADJOINT) {
  TRACE("transfer_asmbl_accel_to_device");
  int blocksize,size_padded,num_blocks_x,num_blocks_y;

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // crust/mantle region
  if( *IREGION == IREGION_CRUST_MANTLE ){
    if( mp->num_interfaces_crust_mantle > 0 ){

      // copies vector buffer values to GPU
      cudaMemcpy(mp->d_send_accel_buffer_crust_mantle, buffer_recv_vector,
                 3*(mp->max_nibool_interfaces_crust_mantle)*(mp->num_interfaces_crust_mantle)*sizeof(realw),
                 cudaMemcpyHostToDevice);

      // assembles values
      blocksize = BLOCKSIZE_TRANSFER;
      size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_crust_mantle)/((double)blocksize)))*blocksize;
      num_blocks_x = size_padded/blocksize;
      num_blocks_y = 1;
      while(num_blocks_x > 65535) {
        num_blocks_x = (int) ceil(num_blocks_x*0.5f);
        num_blocks_y = num_blocks_y*2;
      }

      dim3 grid(num_blocks_x,num_blocks_y);
      dim3 threads(blocksize,1,1);

      if(*FORWARD_OR_ADJOINT == 1) { //assemble forward accel
        assemble_boundary_accel_on_device<<<grid,threads>>>(mp->d_accel_crust_mantle,
                                                            mp->d_send_accel_buffer_crust_mantle,
                                                            mp->num_interfaces_crust_mantle,
                                                            mp->max_nibool_interfaces_crust_mantle,
                                                            mp->d_nibool_interfaces_crust_mantle,
                                                            mp->d_ibool_interfaces_crust_mantle);
      }
      else if(*FORWARD_OR_ADJOINT == 3) { //assemble adjoint accel
        assemble_boundary_accel_on_device<<<grid,threads>>>(mp->d_b_accel_crust_mantle,
                                                            mp->d_send_accel_buffer_crust_mantle,
                                                            mp->num_interfaces_crust_mantle,
                                                            mp->max_nibool_interfaces_crust_mantle,
                                                            mp->d_nibool_interfaces_crust_mantle,
                                                            mp->d_ibool_interfaces_crust_mantle);
      }
    }
  }

  // inner core region
  if( *IREGION == IREGION_INNER_CORE ){
    if( mp->num_interfaces_inner_core > 0 ){
      // copies buffer values to GPU
      cudaMemcpy(mp->d_send_accel_buffer_inner_core, buffer_recv_vector,
                 3*(mp->max_nibool_interfaces_inner_core)*(mp->num_interfaces_inner_core)*sizeof(realw),
                 cudaMemcpyHostToDevice);

      // assembles values
      blocksize = BLOCKSIZE_TRANSFER;
      size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_inner_core)/((double)blocksize)))*blocksize;
      num_blocks_x = size_padded/blocksize;
      num_blocks_y = 1;
      while(num_blocks_x > 65535) {
        num_blocks_x = (int) ceil(num_blocks_x*0.5f);
        num_blocks_y = num_blocks_y*2;
      }

      dim3 grid(num_blocks_x,num_blocks_y);
      dim3 threads(blocksize,1,1);

      if(*FORWARD_OR_ADJOINT == 1) { //assemble forward accel
        assemble_boundary_accel_on_device<<<grid,threads>>>(mp->d_accel_inner_core,
                                                            mp->d_send_accel_buffer_inner_core,
                                                            mp->num_interfaces_inner_core,
                                                            mp->max_nibool_interfaces_inner_core,
                                                            mp->d_nibool_interfaces_inner_core,
                                                            mp->d_ibool_interfaces_inner_core);
      }
      else if(*FORWARD_OR_ADJOINT == 3) { //assemble adjoint accel
        assemble_boundary_accel_on_device<<<grid,threads>>>(mp->d_b_accel_inner_core,
                                                            mp->d_send_accel_buffer_inner_core,
                                                            mp->num_interfaces_inner_core,
                                                            mp->max_nibool_interfaces_inner_core,
                                                            mp->d_nibool_interfaces_inner_core,
                                                            mp->d_ibool_interfaces_inner_core);
      }
    }
  }


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("transfer_asmbl_accel_to_device");
#endif
}
