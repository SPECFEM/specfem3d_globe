/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  2 . 0
 !               ---------------------------------------
 !
 !          Main authors: Dimitri Komatitsch and Jeroen Tromp
 !    Princeton University, USA and University of Pau / CNRS / INRIA
 ! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
 !                            August 2013
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

      // fills buffer
      d_send_accel_buffer[3*iloc] = d_accel[3*iglob];
      d_send_accel_buffer[3*iloc + 1] = d_accel[3*iglob + 1];
      d_send_accel_buffer[3*iloc + 2] = d_accel[3*iglob + 2];
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

  int blocksize,size_padded;
  int num_blocks_x,num_blocks_y;
  dim3 grid,threads;

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  // crust/mantle region
  if( *IREGION == IREGION_CRUST_MANTLE ){
    if( mp->num_interfaces_crust_mantle > 0 ){

      blocksize = BLOCKSIZE_TRANSFER;
      size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_cm)/((double)blocksize)))*blocksize;

      get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

      grid = dim3(num_blocks_x,num_blocks_y);
      threads = dim3(blocksize,1,1);

      if(*FORWARD_OR_ADJOINT == 1) {
        prepare_boundary_accel_on_device<<<grid,threads>>>(mp->d_accel_crust_mantle,
                                                           mp->d_send_accel_buffer_crust_mantle,
                                                           mp->num_interfaces_crust_mantle,
                                                           mp->max_nibool_interfaces_cm,
                                                           mp->d_nibool_interfaces_crust_mantle,
                                                           mp->d_ibool_interfaces_crust_mantle);

        // copies buffer to CPU
        print_CUDA_error_if_any(cudaMemcpy(send_accel_buffer,mp->d_send_accel_buffer_crust_mantle,
                                           NDIM*mp->max_nibool_interfaces_cm*mp->num_interfaces_crust_mantle*sizeof(realw),
                                           cudaMemcpyDeviceToHost),41000);

      }
      else if(*FORWARD_OR_ADJOINT == 3) {
        // debug
        DEBUG_BACKWARD_ASSEMBLY();

        prepare_boundary_accel_on_device<<<grid,threads>>>(mp->d_b_accel_crust_mantle,
                                                           mp->d_b_send_accel_buffer_crust_mantle,
                                                           mp->num_interfaces_crust_mantle,
                                                           mp->max_nibool_interfaces_cm,
                                                           mp->d_nibool_interfaces_crust_mantle,
                                                           mp->d_ibool_interfaces_crust_mantle);
        // copies buffer to CPU
        print_CUDA_error_if_any(cudaMemcpy(send_accel_buffer,mp->d_b_send_accel_buffer_crust_mantle,
                                           NDIM*mp->max_nibool_interfaces_cm*mp->num_interfaces_crust_mantle*sizeof(realw),
                                           cudaMemcpyDeviceToHost),41001);

      }


    }
  }

  // inner core region
  if( *IREGION == IREGION_INNER_CORE ){
    if( mp->num_interfaces_inner_core > 0 ){

      blocksize = BLOCKSIZE_TRANSFER;
      size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ic)/((double)blocksize)))*blocksize;

      get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

      grid = dim3(num_blocks_x,num_blocks_y);
      threads = dim3(blocksize,1,1);

      if(*FORWARD_OR_ADJOINT == 1) {
        prepare_boundary_accel_on_device<<<grid,threads>>>(mp->d_accel_inner_core,
                                                           mp->d_send_accel_buffer_inner_core,
                                                           mp->num_interfaces_inner_core,
                                                           mp->max_nibool_interfaces_ic,
                                                           mp->d_nibool_interfaces_inner_core,
                                                           mp->d_ibool_interfaces_inner_core);

        // copies buffer to CPU
        print_CUDA_error_if_any(cudaMemcpy(send_accel_buffer,mp->d_send_accel_buffer_inner_core,
                                           NDIM*mp->max_nibool_interfaces_ic*mp->num_interfaces_inner_core*sizeof(realw),
                                           cudaMemcpyDeviceToHost),41000);

      }
      else if(*FORWARD_OR_ADJOINT == 3) {
        // debug
        DEBUG_BACKWARD_ASSEMBLY();

        prepare_boundary_accel_on_device<<<grid,threads>>>(mp->d_b_accel_inner_core,
                                                           mp->d_b_send_accel_buffer_inner_core,
                                                           mp->num_interfaces_inner_core,
                                                           mp->max_nibool_interfaces_ic,
                                                           mp->d_nibool_interfaces_inner_core,
                                                           mp->d_ibool_interfaces_inner_core);
        // copies buffer to CPU
        print_CUDA_error_if_any(cudaMemcpy(send_accel_buffer,mp->d_b_send_accel_buffer_inner_core,
                                           NDIM*mp->max_nibool_interfaces_ic*mp->num_interfaces_inner_core*sizeof(realw),
                                           cudaMemcpyDeviceToHost),41001);
      }

    }
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("transfer_boun_accel_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void assemble_boundary_accel_on_device(realw* d_accel, realw* d_send_accel_buffer,
                                                  int num_interfaces,
                                                  int max_nibool_interfaces,
                                                  int* d_nibool_interfaces,
                                                  int* d_ibool_interfaces) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  int iglob,iloc;

  for( int iinterface=0; iinterface < num_interfaces; iinterface++) {
    if(id < d_nibool_interfaces[iinterface]) {

      iloc = id + max_nibool_interfaces*iinterface;
      iglob = d_ibool_interfaces[iloc]-1;

      // assembles acceleration: adds contributions from buffer array
      atomicAdd(&d_accel[3*iglob],d_send_accel_buffer[3*iloc]);
      atomicAdd(&d_accel[3*iglob + 1],d_send_accel_buffer[3*iloc + 1]);
      atomicAdd(&d_accel[3*iglob + 2],d_send_accel_buffer[3*iloc + 2]);
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

  int blocksize,size_padded;
  int num_blocks_x,num_blocks_y;
  dim3 grid,threads;

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // crust/mantle region
  if( *IREGION == IREGION_CRUST_MANTLE ){
    if( mp->num_interfaces_crust_mantle > 0 ){
      // assembles values
      blocksize = BLOCKSIZE_TRANSFER;
      size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_cm)/((double)blocksize)))*blocksize;

      get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

      grid = dim3(num_blocks_x,num_blocks_y);
      threads = dim3(blocksize,1,1);

      if(*FORWARD_OR_ADJOINT == 1) {
        // copies vector buffer values to GPU
        print_CUDA_error_if_any(cudaMemcpy(mp->d_send_accel_buffer_crust_mantle, buffer_recv_vector,
                                           NDIM*(mp->max_nibool_interfaces_cm)*(mp->num_interfaces_crust_mantle)*sizeof(realw),
                                           cudaMemcpyHostToDevice),41000);

        //assemble forward accel
        assemble_boundary_accel_on_device<<<grid,threads>>>(mp->d_accel_crust_mantle,
                                                            mp->d_send_accel_buffer_crust_mantle,
                                                            mp->num_interfaces_crust_mantle,
                                                            mp->max_nibool_interfaces_cm,
                                                            mp->d_nibool_interfaces_crust_mantle,
                                                            mp->d_ibool_interfaces_crust_mantle);
      }
      else if(*FORWARD_OR_ADJOINT == 3) {
        // debug
        DEBUG_BACKWARD_ASSEMBLY();

        // copies vector buffer values to GPU
        print_CUDA_error_if_any(cudaMemcpy(mp->d_b_send_accel_buffer_crust_mantle, buffer_recv_vector,
                                           NDIM*(mp->max_nibool_interfaces_cm)*(mp->num_interfaces_crust_mantle)*sizeof(realw),
                                           cudaMemcpyHostToDevice),41000);

        //assemble adjoint accel
        assemble_boundary_accel_on_device<<<grid,threads>>>(mp->d_b_accel_crust_mantle,
                                                            mp->d_b_send_accel_buffer_crust_mantle,
                                                            mp->num_interfaces_crust_mantle,
                                                            mp->max_nibool_interfaces_cm,
                                                            mp->d_nibool_interfaces_crust_mantle,
                                                            mp->d_ibool_interfaces_crust_mantle);
      }
    }
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("transfer_asmbl_accel_to_device in crust_mantle");
#endif

  // inner core region
  if( *IREGION == IREGION_INNER_CORE ){
    if( mp->num_interfaces_inner_core > 0 ){
      // assembles values
      blocksize = BLOCKSIZE_TRANSFER;
      size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ic)/((double)blocksize)))*blocksize;

      get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

      grid = dim3(num_blocks_x,num_blocks_y);
      threads = dim3(blocksize,1,1);

      if(*FORWARD_OR_ADJOINT == 1) {
        // copies buffer values to GPU
        print_CUDA_error_if_any(cudaMemcpy(mp->d_send_accel_buffer_inner_core, buffer_recv_vector,
                                           NDIM*(mp->max_nibool_interfaces_ic)*(mp->num_interfaces_inner_core)*sizeof(realw),
                                           cudaMemcpyHostToDevice),41001);

        //assemble forward accel
        assemble_boundary_accel_on_device<<<grid,threads>>>(mp->d_accel_inner_core,
                                                            mp->d_send_accel_buffer_inner_core,
                                                            mp->num_interfaces_inner_core,
                                                            mp->max_nibool_interfaces_ic,
                                                            mp->d_nibool_interfaces_inner_core,
                                                            mp->d_ibool_interfaces_inner_core);
      }
      else if(*FORWARD_OR_ADJOINT == 3) {
        // debug
        DEBUG_BACKWARD_ASSEMBLY();

        // copies buffer values to GPU
        print_CUDA_error_if_any(cudaMemcpy(mp->d_b_send_accel_buffer_inner_core, buffer_recv_vector,
                                           NDIM*(mp->max_nibool_interfaces_ic)*(mp->num_interfaces_inner_core)*sizeof(realw),
                                           cudaMemcpyHostToDevice),41001);

        //assemble adjoint accel
        assemble_boundary_accel_on_device<<<grid,threads>>>(mp->d_b_accel_inner_core,
                                                            mp->d_b_send_accel_buffer_inner_core,
                                                            mp->num_interfaces_inner_core,
                                                            mp->max_nibool_interfaces_ic,
                                                            mp->d_nibool_interfaces_inner_core,
                                                            mp->d_ibool_interfaces_inner_core);
      }
    }
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("transfer_asmbl_accel_to_device in inner_core");
#endif
}
