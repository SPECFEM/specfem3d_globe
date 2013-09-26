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

#include "config.h"
#include "mesh_constants_cuda.h"


/* ----------------------------------------------------------------------------------------------- */

// elastic wavefield

// KERNEL 1
/* ----------------------------------------------------------------------------------------------- */


__global__ void UpdateDispVeloc_kernel(realw* displ,
                                       realw* veloc,
                                       realw* accel,
                                       int size,
                                       realw deltat,
                                       realw deltatsqover2,
                                       realw deltatover2) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if(id < size) {
    displ[id] = displ[id] + deltat*veloc[id] + deltatsqover2*accel[id];
    veloc[id] = veloc[id] + deltatover2*accel[id];
    accel[id] = 0.0f; // can do this using memset...not sure if faster
  }
}

/* ----------------------------------------------------------------------------------------------- */

// KERNEL 1
// inner core

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(update_displacement_ic_cuda,
              UPDATE_DISPLACMENT_IC_CUDA)(long* Mesh_pointer_f,
                                          realw* deltat_f,
                                          realw* deltatsqover2_f,
                                          realw* deltatover2_f,
                                          int* FORWARD_OR_ADJOINT) {

TRACE("update_displacement_ic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  int size = NDIM * mp->NGLOB_INNER_CORE;

  // debug
#if DEBUG_BACKWARD_SIMULATIONS == 1
  realw max_d,max_v,max_a;
  max_d = get_device_array_maximum_value(mp->d_b_displ_inner_core, size);
  max_v = get_device_array_maximum_value(mp->d_b_veloc_inner_core, size);
  max_a = get_device_array_maximum_value(mp->d_b_accel_inner_core, size);
  printf("rank %d - max inner_core displ: %f veloc: %f accel: %f\n",mp->myrank,max_d,max_v,max_a);
  fflush(stdout);
  synchronize_mpi();
#endif

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  realw deltat = *deltat_f;
  realw deltatsqover2 = *deltatsqover2_f;
  realw deltatover2 = *deltatover2_f;

  if( *FORWARD_OR_ADJOINT == 1 ){
    //launch kernel
    UpdateDispVeloc_kernel<<<grid,threads>>>(mp->d_displ_inner_core,
                                             mp->d_veloc_inner_core,
                                             mp->d_accel_inner_core,
                                             size,deltat,deltatsqover2,deltatover2);
  }else if( *FORWARD_OR_ADJOINT == 3 ){
    // debug
    DEBUG_BACKWARD_UPDATE();

    // kernel for backward fields
    UpdateDispVeloc_kernel<<<grid,threads>>>(mp->d_b_displ_inner_core,
                                             mp->d_b_veloc_inner_core,
                                             mp->d_b_accel_inner_core,
                                             size,deltat,deltatsqover2,deltatover2);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("update_displacement_ic_cuda");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// KERNEL 1
// crust/mantle

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(update_displacement_cm_cuda,
              UPDATE_DISPLACMENT_CM_CUDA)(long* Mesh_pointer_f,
                                          realw* deltat_f,
                                          realw* deltatsqover2_f,
                                          realw* deltatover2_f,
                                          int* FORWARD_OR_ADJOINT) {

  TRACE("update_displacement_cm_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  int size = NDIM * mp->NGLOB_CRUST_MANTLE;

  // debug
#if DEBUG_BACKWARD_SIMULATIONS == 1
  realw max_d,max_v,max_a;
  max_d = get_device_array_maximum_value(mp->d_b_displ_crust_mantle, size);
  max_v = get_device_array_maximum_value(mp->d_b_veloc_crust_mantle, size);
  max_a = get_device_array_maximum_value(mp->d_b_accel_crust_mantle, size);
  printf("rank %d - max crust_mantle displ: %f veloc: %f accel: %f\n",mp->myrank,max_d,max_v,max_a);
  fflush(stdout);
  synchronize_mpi();
#endif

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  realw deltat = *deltat_f;
  realw deltatsqover2 = *deltatsqover2_f;
  realw deltatover2 = *deltatover2_f;

  if( *FORWARD_OR_ADJOINT == 1 ){
    //launch kernel
    UpdateDispVeloc_kernel<<<grid,threads>>>(mp->d_displ_crust_mantle,
                                             mp->d_veloc_crust_mantle,
                                             mp->d_accel_crust_mantle,
                                             size,deltat,deltatsqover2,deltatover2);
  }else if( *FORWARD_OR_ADJOINT == 3 ){
    // debug
    DEBUG_BACKWARD_UPDATE();

    // kernel for backward fields
    UpdateDispVeloc_kernel<<<grid,threads>>>(mp->d_b_displ_crust_mantle,
                                             mp->d_b_veloc_crust_mantle,
                                             mp->d_b_accel_crust_mantle,
                                             size,deltat,deltatsqover2,deltatover2);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("update_displacement_cm_cuda");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// acoustic wavefield

// KERNEL 1
/* ----------------------------------------------------------------------------------------------- */

__global__ void UpdatePotential_kernel(realw* potential_acoustic,
                                       realw* potential_dot_acoustic,
                                       realw* potential_dot_dot_acoustic,
                                       int size,
                                       realw deltat,
                                       realw deltatsqover2,
                                       realw deltatover2) {
  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  /* because of block and grid sizing problems, there is a small */
  /* amount of buffer at the end of the calculation */
  if(id < size) {
    potential_acoustic[id] = potential_acoustic[id]
                            + deltat*potential_dot_acoustic[id]
                            + deltatsqover2*potential_dot_dot_acoustic[id];

    potential_dot_acoustic[id] = potential_dot_acoustic[id]
                                + deltatover2*potential_dot_dot_acoustic[id];

    potential_dot_dot_acoustic[id] = 0.0f;
  }
}

/* ----------------------------------------------------------------------------------------------- */

// KERNEL 1
// outer core

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(update_displacement_oc_cuda,
              UPDATE_DISPLACEMENT_OC_cuda)(long* Mesh_pointer_f,
                                           realw* deltat_f,
                                           realw* deltatsqover2_f,
                                           realw* deltatover2_f,
                                           int* FORWARD_OR_ADJOINT) {

  TRACE("update_displacement_oc_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  int size = mp->NGLOB_OUTER_CORE;

  // debug
#if DEBUG_BACKWARD_SIMULATIONS == 1
  realw max_d,max_v,max_a;
  max_d = get_device_array_maximum_value(mp->d_b_displ_outer_core, size);
  max_v = get_device_array_maximum_value(mp->d_b_veloc_outer_core, size);
  max_a = get_device_array_maximum_value(mp->d_b_accel_outer_core, size);
  printf("rank %d - max outer_core displ: %f veloc: %f accel: %f\n",mp->myrank,max_d,max_v,max_a);
  fflush(stdout);
  synchronize_mpi();
#endif

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  realw deltat = *deltat_f;
  realw deltatsqover2 = *deltatsqover2_f;
  realw deltatover2 = *deltatover2_f;

  if( *FORWARD_OR_ADJOINT == 1 ){
    //launch kernel
    UpdatePotential_kernel<<<grid,threads>>>(mp->d_displ_outer_core,
                                             mp->d_veloc_outer_core,
                                             mp->d_accel_outer_core,
                                             size,deltat,deltatsqover2,deltatover2);
  }else if( *FORWARD_OR_ADJOINT == 3 ){
    // debug
    DEBUG_BACKWARD_UPDATE();

    UpdatePotential_kernel<<<grid,threads>>>(mp->d_b_displ_outer_core,
                                             mp->d_b_veloc_outer_core,
                                             mp->d_b_accel_outer_core,
                                             size,deltat,deltatsqover2,deltatover2);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("update_displacement_oc_cuda");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// KERNEL 3
//
// elastic domains: crust/mantle and inner core regions

/* ----------------------------------------------------------------------------------------------- */

// unused...
/*
__global__ void kernel_3_cuda_device(realw* veloc,
                                     realw* accel,
                                     int size,
                                     realw deltatover2,
                                     realw two_omega_earth,
                                     realw* rmassx,
                                     realw* rmassy,
                                     realw* rmassz) {
  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if(id < size) {
    // note: update adds rotational acceleration in case two_omega_earth is non-zero
    accel[3*id] = accel[3*id]*rmassx[id] + two_omega_earth*veloc[3*id+1]; // (2,i);
    accel[3*id+1] = accel[3*id+1]*rmassy[id] - two_omega_earth*veloc[3*id]; //(1,i);
    accel[3*id+2] = accel[3*id+2]*rmassz[id];

    veloc[3*id] = veloc[3*id] + deltatover2*accel[3*id];
    veloc[3*id+1] = veloc[3*id+1] + deltatover2*accel[3*id+1];
    veloc[3*id+2] = veloc[3*id+2] + deltatover2*accel[3*id+2];
  }
}
*/

/* ----------------------------------------------------------------------------------------------- */

__global__ void multiply_accel_elastic_cuda_device(realw* accel,
                                                   realw* veloc,
                                                   int size,
                                                   realw two_omega_earth,
                                                   realw* rmassx,
                                                   realw* rmassy,
                                                   realw* rmassz) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if(id < size) {
    // note: update adds rotational acceleration in case two_omega_earth is non-zero
    accel[3*id] = accel[3*id]*rmassx[id] + two_omega_earth*veloc[3*id + 1]; // (2,i);
    accel[3*id + 1] = accel[3*id + 1]*rmassy[id] - two_omega_earth*veloc[3*id]; //(1,i);
    accel[3*id + 2] = accel[3*id + 2]*rmassz[id];
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(multiply_accel_elastic_cuda,
              MULTIPLY_ACCEL_ELASTIC_CUDA)(long* Mesh_pointer,
                                           int* FORWARD_OR_ADJOINT) {
  TRACE("multiply_accel_elastic_cuda");

  int size_padded,num_blocks_x,num_blocks_y;
  dim3 grid,threads;

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int blocksize = BLOCKSIZE_KERNEL3;

  // multiplies accel with inverse of mass matrix

  // crust/mantle region
  size_padded = ((int)ceil(((double)mp->NGLOB_CRUST_MANTLE)/((double)blocksize)))*blocksize;

  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  grid = dim3(num_blocks_x,num_blocks_y);
  threads = dim3(blocksize,1,1);

  if( *FORWARD_OR_ADJOINT == 1 ){
    multiply_accel_elastic_cuda_device<<< grid, threads>>>(mp->d_accel_crust_mantle,
                                                           mp->d_veloc_crust_mantle,
                                                           mp->NGLOB_CRUST_MANTLE,
                                                           mp->two_omega_earth,
                                                           mp->d_rmassx_crust_mantle,
                                                           mp->d_rmassy_crust_mantle,
                                                           mp->d_rmassz_crust_mantle);
  }else if( *FORWARD_OR_ADJOINT == 3 ){
    // debug
    DEBUG_BACKWARD_UPDATE();

    multiply_accel_elastic_cuda_device<<< grid, threads>>>(mp->d_b_accel_crust_mantle,
                                                           mp->d_b_veloc_crust_mantle,
                                                           mp->NGLOB_CRUST_MANTLE,
                                                           mp->b_two_omega_earth,
                                                           mp->d_b_rmassx_crust_mantle,
                                                           mp->d_b_rmassy_crust_mantle,
                                                           mp->d_b_rmassz_crust_mantle);
  }

  // inner core region
  size_padded = ((int)ceil(((double)mp->NGLOB_INNER_CORE)/((double)blocksize)))*blocksize;

  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  grid = dim3(num_blocks_x,num_blocks_y);
  threads = dim3(blocksize,1,1);

  if( *FORWARD_OR_ADJOINT == 1 ){
    multiply_accel_elastic_cuda_device<<< grid, threads>>>(mp->d_accel_inner_core,
                                                           mp->d_veloc_inner_core,
                                                           mp->NGLOB_INNER_CORE,
                                                           mp->two_omega_earth,
                                                           mp->d_rmassx_inner_core,
                                                           mp->d_rmassy_inner_core,
                                                           mp->d_rmassz_inner_core);
  }else if( *FORWARD_OR_ADJOINT == 3 ){
    // debug
    DEBUG_BACKWARD_UPDATE();

    multiply_accel_elastic_cuda_device<<< grid, threads>>>(mp->d_b_accel_inner_core,
                                                           mp->d_b_veloc_inner_core,
                                                           mp->NGLOB_INNER_CORE,
                                                           mp->b_two_omega_earth,
                                                           mp->d_b_rmassx_inner_core,
                                                           mp->d_b_rmassy_inner_core,
                                                           mp->d_b_rmassz_inner_core);
  }


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
  exit_on_cuda_error("after multiply_accel_elastic_cuda");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void update_veloc_elastic_cuda_device(realw* veloc,
                                                 realw* accel,
                                                 int size,
                                                 realw deltatover2) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if(id < size) {
    veloc[3*id] = veloc[3*id] + deltatover2*accel[3*id];
    veloc[3*id + 1] = veloc[3*id + 1] + deltatover2*accel[3*id + 1];
    veloc[3*id + 2] = veloc[3*id + 2] + deltatover2*accel[3*id + 2];
  }
}

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(update_veloc_elastic_cuda,
              UPDATE_VELOC_ELASTIC_CUDA)(long* Mesh_pointer,
                                         realw* deltatover2_f,
                                         int* FORWARD_OR_ADJOINT) {

  TRACE("update_veloc_elastic_cuda");

  int size_padded,num_blocks_x,num_blocks_y;
  dim3 grid,threads;

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  realw deltatover2 = *deltatover2_f;

  int blocksize = BLOCKSIZE_KERNEL3;

  // updates velocity

  // crust/mantle region
  size_padded = ((int)ceil(((double)mp->NGLOB_CRUST_MANTLE)/((double)blocksize)))*blocksize;

  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  grid = dim3(num_blocks_x,num_blocks_y);
  threads = dim3(blocksize,1,1);

  if( *FORWARD_OR_ADJOINT == 1 ){
    update_veloc_elastic_cuda_device<<< grid, threads>>>(mp->d_veloc_crust_mantle,
                                                         mp->d_accel_crust_mantle,
                                                         mp->NGLOB_CRUST_MANTLE,
                                                         deltatover2);
  }else if( *FORWARD_OR_ADJOINT == 3 ){
    update_veloc_elastic_cuda_device<<< grid, threads>>>(mp->d_b_veloc_crust_mantle,
                                                         mp->d_b_accel_crust_mantle,
                                                         mp->NGLOB_CRUST_MANTLE,
                                                         deltatover2);
  }

  // inner core region
  size_padded = ((int)ceil(((double)mp->NGLOB_INNER_CORE)/((double)blocksize)))*blocksize;

  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  grid = dim3(num_blocks_x,num_blocks_y);
  threads = dim3(blocksize,1,1);

  if( *FORWARD_OR_ADJOINT == 1 ){
    update_veloc_elastic_cuda_device<<< grid, threads>>>(mp->d_veloc_inner_core,
                                                         mp->d_accel_inner_core,
                                                         mp->NGLOB_INNER_CORE,
                                                         deltatover2);
  }else if( *FORWARD_OR_ADJOINT == 3 ){
    update_veloc_elastic_cuda_device<<< grid, threads>>>(mp->d_b_veloc_inner_core,
                                                         mp->d_b_accel_inner_core,
                                                         mp->NGLOB_INNER_CORE,
                                                         deltatover2);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
  exit_on_cuda_error("after update_veloc_3_b");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// KERNEL 3
//
// acoustic domains: for fluid outer_core

/* ----------------------------------------------------------------------------------------------- */


__global__ void multiply_accel_acoustic_cuda_device(realw* accel,
                                                    int size,
                                                    realw* rmass) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if(id < size) {
    // multiplies pressure with the inverse of the mass matrix
    accel[id] = accel[id]*rmass[id];
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(multiply_accel_acoustic_cuda,
              MULTIPLY_ACCEL_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                            int* FORWARD_OR_ADJOINT) {
  TRACE("multiply_accel_acoustic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int blocksize = BLOCKSIZE_KERNEL3;

  int size_padded = ((int)ceil(((double)mp->NGLOB_OUTER_CORE)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // multiplies accel with inverse of mass matrix
  if( *FORWARD_OR_ADJOINT == 1 ){
    multiply_accel_acoustic_cuda_device<<< grid, threads>>>(mp->d_accel_outer_core,
                                                            mp->NGLOB_OUTER_CORE,
                                                            mp->d_rmass_outer_core);
  }else if( *FORWARD_OR_ADJOINT == 3 ){
    // debug
    DEBUG_BACKWARD_UPDATE();

    multiply_accel_acoustic_cuda_device<<< grid, threads>>>(mp->d_b_accel_outer_core,
                                                            mp->NGLOB_OUTER_CORE,
                                                            mp->d_b_rmass_outer_core);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after multiply_accel_acoustic_cuda");
#endif
}

/* ----------------------------------------------------------------------------------------------- */


__global__ void update_veloc_acoustic_cuda_device(realw* veloc,
                                                  realw* accel,
                                                  int size,
                                                  realw deltatover2) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if(id < size) {
    // Newmark time scheme: corrector term
    veloc[id] = veloc[id] + deltatover2*accel[id];
  }
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(update_veloc_acoustic_cuda,
              UPDATE_VELOC_ACOUSTIC_CUDA)(long* Mesh_pointer,
                                          realw* deltatover2_f,
                                          int* FORWARD_OR_ADJOINT) {

  TRACE("update_veloc_acoustic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  realw deltatover2 = *deltatover2_f;

  int blocksize = BLOCKSIZE_KERNEL3;

  int size_padded = ((int)ceil(((double)mp->NGLOB_OUTER_CORE)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // updates velocity
  if( *FORWARD_OR_ADJOINT == 1 ){
    update_veloc_acoustic_cuda_device<<< grid, threads>>>(mp->d_veloc_outer_core,
                                                          mp->d_accel_outer_core,
                                                          mp->NGLOB_OUTER_CORE,
                                                          deltatover2);
  }else if( *FORWARD_OR_ADJOINT == 3){
    // debug
    DEBUG_BACKWARD_UPDATE();

    update_veloc_acoustic_cuda_device<<< grid, threads>>>(mp->d_b_veloc_outer_core,
                                                          mp->d_b_accel_outer_core,
                                                          mp->NGLOB_OUTER_CORE,
                                                          deltatover2);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
  exit_on_cuda_error("after update_veloc_acoustic_cuda");
#endif
}

