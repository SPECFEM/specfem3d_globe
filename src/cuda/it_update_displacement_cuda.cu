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

//#include <cublas.h>

#include "config.h"
#include "mesh_constants_cuda.h"


//#define CUBLAS_ERROR(s,n)  if (s != CUBLAS_STATUS_SUCCESS) {  \
//fprintf (stderr, "CUBLAS Memory Write Error @ %d\n",n); \
//exit(EXIT_FAILURE); }


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

  /* because of block and grid sizing problems, there is a small */
  /* amount of buffer at the end of the calculation */
  if(id < size) {
    displ[id] = displ[id] + deltat*veloc[id] + deltatsqover2*accel[id];
    veloc[id] = veloc[id] + deltatover2*accel[id];
    accel[id] = 0; // can do this using memset...not sure if faster
  }
}

/* ----------------------------------------------------------------------------------------------- */

// KERNEL 1
// inner core

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(it_update_displacement_ic_cuda,
              IT_UPDATE_DISPLACMENT_IC_CUDA)(long* Mesh_pointer_f,
                                             realw* deltat_F,
                                             realw* deltatsqover2_F,
                                             realw* deltatover2_F,
                                             realw* b_deltat_F,
                                             realw* b_deltatsqover2_F,
                                             realw* b_deltatover2_F) {

TRACE("it_update_displacement_ic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  int size = NDIM * mp->NGLOB_INNER_CORE;

  realw deltat = *deltat_F;
  realw deltatsqover2 = *deltatsqover2_F;
  realw deltatover2 = *deltatover2_F;

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x = size_padded/blocksize;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  //launch kernel
  UpdateDispVeloc_kernel<<<grid,threads>>>(mp->d_displ_inner_core,
                                           mp->d_veloc_inner_core,
                                           mp->d_accel_inner_core,
                                           size,deltat,deltatsqover2,deltatover2);

  // kernel for backward fields
  if(mp->simulation_type == 3) {
    realw b_deltat = *b_deltat_F;
    realw b_deltatsqover2 = *b_deltatsqover2_F;
    realw b_deltatover2 = *b_deltatover2_F;

    UpdateDispVeloc_kernel<<<grid,threads>>>(mp->d_b_displ_inner_core,
                                             mp->d_b_veloc_inner_core,
                                             mp->d_b_accel_inner_core,
                                             size,b_deltat,b_deltatsqover2,b_deltatover2);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("it_update_displacement_ic_cuda");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// KERNEL 1
// crust/mantle

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(it_update_displacement_cm_cuda,
              IT_UPDATE_DISPLACMENT_CM_CUDA)(long* Mesh_pointer_f,
                                             realw* deltat_F,
                                             realw* deltatsqover2_F,
                                             realw* deltatover2_F,
                                             realw* b_deltat_F,
                                             realw* b_deltatsqover2_F,
                                             realw* b_deltatover2_F) {

  TRACE("it_update_displacement_cm_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  int size = NDIM * mp->NGLOB_CRUST_MANTLE;

  realw deltat = *deltat_F;
  realw deltatsqover2 = *deltatsqover2_F;
  realw deltatover2 = *deltatover2_F;

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x = size_padded/blocksize;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  //launch kernel
  UpdateDispVeloc_kernel<<<grid,threads>>>(mp->d_displ_crust_mantle,
                                           mp->d_veloc_crust_mantle,
                                           mp->d_accel_crust_mantle,
                                           size,deltat,deltatsqover2,deltatover2);

  // kernel for backward fields
  if(mp->simulation_type == 3) {
    realw b_deltat = *b_deltat_F;
    realw b_deltatsqover2 = *b_deltatsqover2_F;
    realw b_deltatover2 = *b_deltatover2_F;

    UpdateDispVeloc_kernel<<<grid,threads>>>(mp->d_b_displ_crust_mantle,
                                             mp->d_b_veloc_crust_mantle,
                                             mp->d_b_accel_crust_mantle,
                                             size,b_deltat,b_deltatsqover2,b_deltatover2);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("it_update_displacement_cm_cuda");
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

    potential_dot_dot_acoustic[id] = 0;
  }
}

/* ----------------------------------------------------------------------------------------------- */

// KERNEL 1
// outer core

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(it_update_displacement_oc_cuda,
              IT_UPDATE_DISPLACEMENT_OC_cuda)(long* Mesh_pointer_f,
                                               realw* deltat_F,
                                               realw* deltatsqover2_F,
                                               realw* deltatover2_F,
                                               realw* b_deltat_F,
                                               realw* b_deltatsqover2_F,
                                               realw* b_deltatover2_F) {

  TRACE("it_update_displacement_oc_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  int size = mp->NGLOB_OUTER_CORE;

  realw deltat = *deltat_F;
  realw deltatsqover2 = *deltatsqover2_F;
  realw deltatover2 = *deltatover2_F;

  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x = size_padded/blocksize;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  //launch kernel
  UpdatePotential_kernel<<<grid,threads>>>(mp->d_displ_outer_core,
                                           mp->d_veloc_outer_core,
                                           mp->d_accel_outer_core,
                                           size,deltat,deltatsqover2,deltatover2);

  if(mp->simulation_type == 3) {
    realw b_deltat = *b_deltat_F;
    realw b_deltatsqover2 = *b_deltatsqover2_F;
    realw b_deltatover2 = *b_deltatover2_F;

    UpdatePotential_kernel<<<grid,threads>>>(mp->d_b_displ_outer_core,
                                             mp->d_b_veloc_outer_core,
                                             mp->d_b_accel_outer_core,
                                             size,b_deltat,b_deltatsqover2,b_deltatover2);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("it_update_displacement_oc_cuda");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// KERNEL 3
//
// crust/mantle and inner core regions

/* ----------------------------------------------------------------------------------------------- */

__global__ void kernel_3_cuda_device(realw* veloc,
                                     realw* accel, int size,
                                     realw deltatover2,
                                     realw* rmassx,
             realw* rmassy,
             realw* rmassz) {
  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  /* because of block and grid sizing problems, there is a small */
  /* amount of buffer at the end of the calculation */
  if(id < size) {
    accel[3*id] = accel[3*id]*rmassx[id];
    accel[3*id+1] = accel[3*id+1]*rmassy[id];
    accel[3*id+2] = accel[3*id+2]*rmassz[id];

    veloc[3*id] = veloc[3*id] + deltatover2*accel[3*id];
    veloc[3*id+1] = veloc[3*id+1] + deltatover2*accel[3*id+1];
    veloc[3*id+2] = veloc[3*id+2] + deltatover2*accel[3*id+2];
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void kernel_3_accel_cuda_device(realw* accel,
                                           int size,
                                           realw* rmassx,
             realw* rmassy,
             realw* rmassz) {
  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  /* because of block and grid sizing problems, there is a small */
  /* amount of buffer at the end of the calculation */
  if(id < size) {
    accel[3*id] = accel[3*id]*rmassx[id];
    accel[3*id+1] = accel[3*id+1]*rmassy[id];
    accel[3*id+2] = accel[3*id+2]*rmassz[id];
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void kernel_3_veloc_cuda_device(realw* veloc,
                                           realw* accel,
                                           int size,
                                           realw deltatover2) {
  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  /* because of block and grid sizing problems, there is a small */
  /* amount of buffer at the end of the calculation */
  if(id < size) {
    veloc[3*id] = veloc[3*id] + deltatover2*accel[3*id];
    veloc[3*id+1] = veloc[3*id+1] + deltatover2*accel[3*id+1];
    veloc[3*id+2] = veloc[3*id+2] + deltatover2*accel[3*id+2];
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(kernel_3_a_cuda,
              KERNEL_3_A_CUDA)(long* Mesh_pointer,
                               realw* deltatover2_F,
                               int* SIMULATION_TYPE_f,
                               realw* b_deltatover2_F,
                               int* OCEANS,
             int* NCHUNKS_VAL) {
  TRACE("kernel_3_a_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int SIMULATION_TYPE = *SIMULATION_TYPE_f;
  realw deltatover2 = *deltatover2_F;
  realw b_deltatover2 = *b_deltatover2_F;

  int blocksize = BLOCKSIZE_KERNEL3;
  int size_padded = ((int)ceil(((double)mp->NGLOB_CRUST_MANTLE)/((double)blocksize)))*blocksize;

  int num_blocks_x = size_padded/blocksize;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // crust/mantle region only
  // check whether we can update accel and veloc, or only accel at this point
  if( *OCEANS == 0 ){
    // updates both, accel and veloc

    if( *NCHUNKS_VAL != 6 && mp->absorbing_conditions){
      kernel_3_cuda_device<<< grid, threads>>>(mp->d_veloc_crust_mantle,
                 mp->d_accel_crust_mantle,
                 mp->NGLOB_CRUST_MANTLE,
                 deltatover2,
                 mp->d_rmassx_crust_mantle,
                 mp->d_rmassy_crust_mantle,
                 mp->d_rmassz_crust_mantle);

      if(SIMULATION_TYPE == 3){
  kernel_3_cuda_device<<< grid, threads>>>(mp->d_b_veloc_crust_mantle,
             mp->d_b_accel_crust_mantle,
             mp->NGLOB_CRUST_MANTLE,
             b_deltatover2,
             mp->d_rmassx_crust_mantle,
             mp->d_rmassy_crust_mantle,
             mp->d_rmassz_crust_mantle);
      }
    }else{
      kernel_3_cuda_device<<< grid, threads>>>(mp->d_veloc_crust_mantle,
                 mp->d_accel_crust_mantle,
                 mp->NGLOB_CRUST_MANTLE,
                 deltatover2,
                 mp->d_rmassz_crust_mantle,
                 mp->d_rmassz_crust_mantle,
                 mp->d_rmassz_crust_mantle);

      if(SIMULATION_TYPE == 3){
  kernel_3_cuda_device<<< grid, threads>>>(mp->d_b_veloc_crust_mantle,
             mp->d_b_accel_crust_mantle,
             mp->NGLOB_CRUST_MANTLE,
             b_deltatover2,
             mp->d_rmassz_crust_mantle,
             mp->d_rmassz_crust_mantle,
             mp->d_rmassz_crust_mantle);
      }
    }

  }else{
    // updates only accel

    if( *NCHUNKS_VAL != 6 && mp->absorbing_conditions){
      kernel_3_accel_cuda_device<<< grid, threads>>>(mp->d_accel_crust_mantle,
                 mp->NGLOB_CRUST_MANTLE,
                 mp->d_rmassx_crust_mantle,
                 mp->d_rmassy_crust_mantle,
                 mp->d_rmassz_crust_mantle);

      if(SIMULATION_TYPE == 3) {
  kernel_3_accel_cuda_device<<< grid, threads>>>(mp->d_b_accel_crust_mantle,
                   mp->NGLOB_CRUST_MANTLE,
                   mp->d_rmassx_crust_mantle,
                   mp->d_rmassy_crust_mantle,
                   mp->d_rmassz_crust_mantle);
      }
    }else{
      kernel_3_accel_cuda_device<<< grid, threads>>>(mp->d_accel_crust_mantle,
                 mp->NGLOB_CRUST_MANTLE,
                 mp->d_rmassz_crust_mantle,
                 mp->d_rmassz_crust_mantle,
                 mp->d_rmassz_crust_mantle);

      if(SIMULATION_TYPE == 3) {
  kernel_3_accel_cuda_device<<< grid, threads>>>(mp->d_b_accel_crust_mantle,
                   mp->NGLOB_CRUST_MANTLE,
                   mp->d_rmassz_crust_mantle,
                   mp->d_rmassz_crust_mantle,
                   mp->d_rmassz_crust_mantle);
      }
    }
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
  exit_on_cuda_error("after kernel_3_a");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(kernel_3_b_cuda,
              KERNEL_3_B_CUDA)(long* Mesh_pointer,
                               realw* deltatover2_F,
                               int* SIMULATION_TYPE_f,
                               realw* b_deltatover2_F,
                               int* OCEANS) {
  TRACE("kernel_3_b_cuda");
  int size_padded,num_blocks_x,num_blocks_y;

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int SIMULATION_TYPE = *SIMULATION_TYPE_f;
  realw deltatover2 = *deltatover2_F;
  realw b_deltatover2 = *b_deltatover2_F;

  int blocksize = BLOCKSIZE_KERNEL3;

  // crust/mantle region
  // in case of ocean loads, we still have to update the velocity for crust/mantle region
  if( *OCEANS ){
    size_padded = ((int)ceil(((double)mp->NGLOB_CRUST_MANTLE)/((double)blocksize)))*blocksize;
    num_blocks_x = size_padded/blocksize;
    num_blocks_y = 1;
    while(num_blocks_x > 65535) {
      num_blocks_x = (int) ceil(num_blocks_x*0.5f);
      num_blocks_y = num_blocks_y*2;
    }
    dim3 grid1(num_blocks_x,num_blocks_y);
    dim3 threads1(blocksize,1,1);

    // updates only veloc at this point
    kernel_3_veloc_cuda_device<<< grid1, threads1>>>(mp->d_veloc_crust_mantle,
                                                     mp->d_accel_crust_mantle,
                                                     mp->NGLOB_CRUST_MANTLE,
                                                     deltatover2);

    if(SIMULATION_TYPE == 3) {
      kernel_3_veloc_cuda_device<<< grid1, threads1>>>(mp->d_b_veloc_crust_mantle,
                                                       mp->d_b_accel_crust_mantle,
                                                       mp->NGLOB_CRUST_MANTLE,
                                                       b_deltatover2);
    }
  }

  // inner core
  size_padded = ((int)ceil(((double)mp->NGLOB_INNER_CORE)/((double)blocksize)))*blocksize;
  num_blocks_x = size_padded/blocksize;
  num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }
  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // updates both, accel and veloc
  kernel_3_cuda_device<<< grid, threads>>>(mp->d_veloc_inner_core,
                                           mp->d_accel_inner_core,
                                           mp->NGLOB_INNER_CORE,
                                           deltatover2,
             mp->d_rmass_inner_core,
             mp->d_rmass_inner_core,
             mp->d_rmass_inner_core);

  if(SIMULATION_TYPE == 3) {
    kernel_3_cuda_device<<< grid, threads>>>(mp->d_b_veloc_inner_core,
                                             mp->d_b_accel_inner_core,
                                             mp->NGLOB_INNER_CORE,
                                             b_deltatover2,
               mp->d_rmass_inner_core,
               mp->d_rmass_inner_core,
               mp->d_rmass_inner_core);
  }


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
  exit_on_cuda_error("after kernel_3_b");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// KERNEL 3
//
// for outer_core

/* ----------------------------------------------------------------------------------------------- */


__global__ void kernel_3_outer_core_cuda_device(realw* veloc,
                                                realw* accel,int size,
                                                realw deltatover2,
                                                realw* rmass) {
  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  /* because of block and grid sizing problems, there is a small */
  /* amount of buffer at the end of the calculation */
  if(id < size) {
    // multiplies pressure with the inverse of the mass matrix
    accel[id] = accel[id]*rmass[id];

    // Newmark time scheme: corrector term
    veloc[id] = veloc[id] + deltatover2*accel[id];
  }
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(kernel_3_outer_core_cuda,
              KERNEL_3_OUTER_CORE_CUDA)(long* Mesh_pointer,
                                        realw* deltatover2_F,
                                        int* SIMULATION_TYPE_f,
                                        realw* b_deltatover2_F) {

  TRACE("kernel_3_outer_core_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

  int SIMULATION_TYPE = *SIMULATION_TYPE_f;
  realw deltatover2 = *deltatover2_F;
  realw b_deltatover2 = *b_deltatover2_F;

  int blocksize = BLOCKSIZE_KERNEL3;
  int size_padded = ((int)ceil(((double)mp->NGLOB_OUTER_CORE)/((double)blocksize)))*blocksize;
  int num_blocks_x = size_padded/blocksize;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }
  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  kernel_3_outer_core_cuda_device<<< grid, threads>>>(mp->d_veloc_outer_core,
                                                      mp->d_accel_outer_core,
                                                      mp->NGLOB_OUTER_CORE,
                                                      deltatover2,mp->d_rmass_outer_core);

  if(SIMULATION_TYPE == 3) {
    kernel_3_outer_core_cuda_device<<< grid, threads>>>(mp->d_b_veloc_outer_core,
                                                        mp->d_b_accel_outer_core,
                                                        mp->NGLOB_OUTER_CORE,
                                                        b_deltatover2,mp->d_rmass_outer_core);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("checking updatedispl_kernel launch...with %dx%d blocks\n",num_blocks_x,num_blocks_y);
  exit_on_cuda_error("after kernel_3_outer_core");
#endif
}

