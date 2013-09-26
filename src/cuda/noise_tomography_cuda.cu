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

//#ifdef WITH_MPI
//#include <mpi.h>
//#endif

#include <sys/types.h>
#include <unistd.h>

#include "config.h"
#include "mesh_constants_cuda.h"


/* ----------------------------------------------------------------------------------------------- */
//unused routines...
/*
extern "C"
void FC_FUNC_(fortranflush,FORTRANFLUSH)(int* rank){
TRACE("fortranflush");

  fflush(stdout);
  fflush(stderr);
  printf("Flushing proc %d!\n",*rank);
}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(fortranprint,FORTRANPRINT)(int* id) {
TRACE("fortranprint");

  int procid;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
#else
  procid = 0;
#endif
  printf("%d: sends msg_id %d\n",procid,*id);
}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(fortranprintf,FORTRANPRINTF)(realw* val) {
TRACE("fortranprintf");

  int procid;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
#else
  procid = 0;
#endif
  printf("%d: sends val %e\n",procid,*val);
}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(fortranprintd,FORTRANPRINTD)(double* val) {
TRACE("fortranprintd");

  int procid;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
#else
  procid = 0;
#endif
  printf("%d: sends val %e\n",procid,*val);
}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
// randomize displ for testing
extern "C"
void FC_FUNC_(make_displ_rand,MAKE_DISPL_RAND)(long* Mesh_pointer_f,realw* h_displ) {
TRACE("make_displ_rand");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper
  // realw* displ_rnd = (realw*)malloc(mp->NGLOB_AB*3*sizeof(realw));
  for(int i=0;i<mp->NGLOB_AB*3;i++) {
    h_displ[i] = rand();
  }
  cudaMemcpy(mp->d_displ,h_displ,mp->NGLOB_AB*3*sizeof(realw),cudaMemcpyHostToDevice);
}
*/

/* ----------------------------------------------------------------------------------------------- */

// noise transfer surface movie

/* ----------------------------------------------------------------------------------------------- */

__global__ void noise_transfer_surface_to_host_kernel(int* ibelm_top,
                                                      int nspec_top,
                                                      int* ibool,
                                                      realw* displ,
                                                      realw* noise_surface_movie) {
  int igll = threadIdx.x;
  int iface = blockIdx.x + blockIdx.y*gridDim.x;

  // int id = tx + blockIdx.x*blockDim.x + blockIdx.y*blockDim.x*gridDim.x;

  if(iface < nspec_top) {
    int ispec = ibelm_top[iface]-1; //-1 for C-based indexing

    int k = NGLLX-1;
    int j = (igll/NGLLX);
    int i = (igll-j*NGLLX);

    int iglob = ibool[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

    noise_surface_movie[INDEX3(NDIM,NGLL2,0,igll,iface)] = displ[iglob*3];
    noise_surface_movie[INDEX3(NDIM,NGLL2,1,igll,iface)] = displ[iglob*3+1];
    noise_surface_movie[INDEX3(NDIM,NGLL2,2,igll,iface)] = displ[iglob*3+2];
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(noise_transfer_surface_to_host,
              NOISE_TRANSFER_SURFACE_TO_HOST)(long* Mesh_pointer_f,
                                              realw* h_noise_surface_movie) {
TRACE("noise_transfer_surface_to_host");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nspec2D_top_crust_mantle,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(NGLL2,1,1);

  noise_transfer_surface_to_host_kernel<<<grid,threads>>>(mp->d_ibelm_top_crust_mantle,
                                                          mp->nspec2D_top_crust_mantle,
                                                          mp->d_ibool_crust_mantle,
                                                          mp->d_displ_crust_mantle,
                                                          mp->d_noise_surface_movie);

  // copies noise array to CPU
  cudaMemcpy(h_noise_surface_movie,mp->d_noise_surface_movie,
             NDIM*NGLL2*(mp->nspec2D_top_crust_mantle)*sizeof(realw),cudaMemcpyDeviceToHost);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("noise_transfer_surface_to_host");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// NOISE add source master

/* ----------------------------------------------------------------------------------------------- */

__global__ void noise_add_source_master_rec_cuda_kernel(int* ibool,
                                                        int* ispec_selected_rec,
                                                        int irec_master_noise,
                                                        realw* accel,
                                                        realw* noise_sourcearray,
                                                        int it) {
  int tx = threadIdx.x;
  int ispec = ispec_selected_rec[irec_master_noise]-1;
  int iglob = ibool[tx + NGLL3*ispec]-1;

  atomicAdd(&accel[iglob*3  ],noise_sourcearray[  3*tx + 3*NGLL3*it]);
  atomicAdd(&accel[iglob*3+1],noise_sourcearray[1+3*tx + 3*NGLL3*it]);
  atomicAdd(&accel[iglob*3+2],noise_sourcearray[2+3*tx + 3*NGLL3*it]);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(noise_add_source_master_rec_cu,
              NOISE_ADD_SOURCE_MASTER_REC_CU)(long* Mesh_pointer_f,
                                              int* it_f,
                                              int* irec_master_noise_f,
                                              int* islice_selected_rec) {

  TRACE("noise_add_source_master_rec_cu");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  int it = *it_f - 1; // -1 for Fortran -> C indexing differences
  int irec_master_noise = *irec_master_noise_f-1;

  dim3 grid(1,1,1);
  dim3 threads(NGLL3,1,1);

  // adds noise source at master location
  if(mp->myrank == islice_selected_rec[irec_master_noise]) {
    noise_add_source_master_rec_cuda_kernel<<<grid,threads>>>(mp->d_ibool_crust_mantle,
                                                              mp->d_ispec_selected_rec,
                                                              irec_master_noise,
                                                              mp->d_accel_crust_mantle,
                                                              mp->d_noise_sourcearray,
                                                              it);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("noise_add_source_master_rec_cuda_kernel");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// noise add surface movie

/* ----------------------------------------------------------------------------------------------- */

__global__ void noise_add_surface_movie_cuda_kernel(realw* accel,
                                                    int* ibool,
                                                    int* ibelm_top,
                                                    int nspec_top,
                                                    realw* noise_surface_movie,
                                                    realw* normal_x_noise,
                                                    realw* normal_y_noise,
                                                    realw* normal_z_noise,
                                                    realw* mask_noise,
                                                    realw* jacobian2D,
                                                    realw* wgllwgll) {

  int igll = threadIdx.x;
  int iface = blockIdx.x + gridDim.x*blockIdx.y; // surface element id

  // when nspec_top > MAXIMUM_GRID_DIM, but mod(nspec_top,2) > 0, we end up with an extra block.
  if(iface < nspec_top) {

    int ispec = ibelm_top[iface]-1;

    int k = NGLLX - 1;
    int j = (igll/NGLLX);
    int i = (igll-j*NGLLX);

    int iglob = ibool[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

    int ipoin = NGLL2*iface + igll;
    realw normal_x = normal_x_noise[ipoin];
    realw normal_y = normal_y_noise[ipoin];
    realw normal_z = normal_z_noise[ipoin];

    realw eta = (noise_surface_movie[INDEX3(NDIM,NGLL2,0,igll,iface)]*normal_x +
                 noise_surface_movie[INDEX3(NDIM,NGLL2,1,igll,iface)]*normal_y +
                 noise_surface_movie[INDEX3(NDIM,NGLL2,2,igll,iface)]*normal_z);

    // weighted jacobian
    realw jacobianw = wgllwgll[k*NGLLX+i]*jacobian2D[igll+NGLL2*iface];

    // note: check error from cuda-memcheck and ddt seems "incorrect", because we
    //          are passing a __constant__ variable pointer around like it was
    //          made using cudaMalloc, which *may* be "incorrect", but produces
    //          correct results.

    // note: global version uses jacobian2D arrays which do not include gll weights wgllwgll,
    //          thus we have to explicitly add: wgllwgll(..) * jacobian2D(..)

    atomicAdd(&accel[iglob*3]  ,eta*mask_noise[ipoin]*normal_x*jacobianw);
    atomicAdd(&accel[iglob*3+1],eta*mask_noise[ipoin]*normal_y*jacobianw);
    atomicAdd(&accel[iglob*3+2],eta*mask_noise[ipoin]*normal_z*jacobianw);

  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(noise_add_surface_movie_cuda,
              NOISE_ADD_SURFACE_MOVIE_CUDA)(long* Mesh_pointer_f,
                                            realw* h_noise_surface_movie) {

  TRACE("noise_add_surface_movie_cuda");


  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nspec2D_top_crust_mantle,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(NGLL2,1,1);

  // copies surface movie to GPU
  cudaMemcpy(mp->d_noise_surface_movie,h_noise_surface_movie,
             NDIM*NGLL2*(mp->nspec2D_top_crust_mantle)*sizeof(realw),cudaMemcpyHostToDevice);

  switch(mp->noise_tomography) {
  case 2:
    // adds surface source to forward field
    noise_add_surface_movie_cuda_kernel<<<grid,threads>>>(mp->d_accel_crust_mantle,
                                                          mp->d_ibool_crust_mantle,
                                                          mp->d_ibelm_top_crust_mantle,
                                                          mp->nspec2D_top_crust_mantle,
                                                          mp->d_noise_surface_movie,
                                                          mp->d_normal_x_noise,
                                                          mp->d_normal_y_noise,
                                                          mp->d_normal_z_noise,
                                                          mp->d_mask_noise,
                                                          mp->d_jacobian2D_top_crust_mantle,
                                                          mp->d_wgllwgll_xy);
    break;

  case 3:
    // adds surface source to adjoint (backward) field
    noise_add_surface_movie_cuda_kernel<<<grid,threads>>>(mp->d_b_accel_crust_mantle,
                                                          mp->d_ibool_crust_mantle,
                                                          mp->d_ibelm_top_crust_mantle,
                                                          mp->nspec2D_top_crust_mantle,
                                                          mp->d_noise_surface_movie,
                                                          mp->d_normal_x_noise,
                                                          mp->d_normal_y_noise,
                                                          mp->d_normal_z_noise,
                                                          mp->d_mask_noise,
                                                          mp->d_jacobian2D_top_crust_mantle,
                                                          mp->d_wgllwgll_xy);
    break;
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("noise_read_add_surface_movie_cuda_kernel");
#endif
}
