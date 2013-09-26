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

#include <sys/types.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"


/* ----------------------------------------------------------------------------------------------- */

// elastic domain sources

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_add_sources_kernel(realw* accel,
                                           int* ibool,
                                           realw* sourcearrays,
                                           double* stf_pre_compute,
                                           int myrank,
                                           int* islice_selected_source,
                                           int* ispec_selected_source,
                                           int NSOURCES) {
  int ispec,iglob;
  realw stf;

  int i = threadIdx.x;
  int j = threadIdx.y;
  int k = threadIdx.z;
  int isource  = blockIdx.x + gridDim.x*blockIdx.y; // bx

  // when NSOURCES > MAXIMUM_GRID_DIM, but mod(nspec_top,2) > 0, we end up with an extra block.
  if(isource < NSOURCES) {
    if(myrank == islice_selected_source[isource]) {

      ispec = ispec_selected_source[isource]-1;

      stf = (realw) stf_pre_compute[isource];
      iglob = ibool[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

      // note: for global version, sourcearrays has dimensions
      //            sourcearrays(NDIM,NGLLX,NGLLY,NGLLZ,NSOURCES)
      atomicAdd(&accel[3*iglob], sourcearrays[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,
                                                     0,i,j,k,isource)]*stf);
      atomicAdd(&accel[3*iglob+1], sourcearrays[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,
                                                       1,i,j,k,isource)]*stf);
      atomicAdd(&accel[3*iglob+2], sourcearrays[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,
                                                       2,i,j,k,isource)]*stf);
    }
  }
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_add_sources_cuda,
              COMPUTE_ADD_SOURCES_CUDA)(long* Mesh_pointer_f,
                                        int* NSOURCESf,
                                        double* h_stf_pre_compute) {

  TRACE("compute_add_sources_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if( mp->nsources_local == 0 ) return;

  int NSOURCES = *NSOURCESf;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(NSOURCES,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLLX,NGLLX,NGLLX);

  // copies source time function buffer values to GPU
  print_CUDA_error_if_any(cudaMemcpy(mp->d_stf_pre_compute,h_stf_pre_compute,
                                     NSOURCES*sizeof(double),cudaMemcpyHostToDevice),71018);

  // adds source contributions
  compute_add_sources_kernel<<<grid,threads>>>(mp->d_accel_crust_mantle,
                                               mp->d_ibool_crust_mantle,
                                               mp->d_sourcearrays,
                                               mp->d_stf_pre_compute,
                                               mp->myrank,
                                               mp->d_islice_selected_source,
                                               mp->d_ispec_selected_source,
                                               NSOURCES);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_add_sources_cuda");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// backward sources

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_add_sources_backward_cuda,
              COMPUTE_ADD_SOURCES_BACKWARD_CUDA)(long* Mesh_pointer_f,
                                                 int* NSOURCESf,
                                                 double* h_stf_pre_compute) {
  TRACE("compute_add_sources_backward_cuda");
  // debug
  DEBUG_BACKWARD_SOURCES();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if( mp->nsources_local == 0 ) return;

  int NSOURCES = *NSOURCESf;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(NSOURCES,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLLX,NGLLX,NGLLX);

  // copies source time function buffer values to GPU
  print_CUDA_error_if_any(cudaMemcpy(mp->d_stf_pre_compute,h_stf_pre_compute,
                                     NSOURCES*sizeof(double),cudaMemcpyHostToDevice),71019);

  compute_add_sources_kernel<<<grid,threads>>>(mp->d_b_accel_crust_mantle,
                                               mp->d_ibool_crust_mantle,
                                               mp->d_sourcearrays,
                                               mp->d_stf_pre_compute,
                                               mp->myrank,
                                               mp->d_islice_selected_source,
                                               mp->d_ispec_selected_source,
                                               NSOURCES);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_add_sources_backward_cuda");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// ADJOINT sources

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_add_sources_adjoint_cuda_kernel(realw* accel,
                                                        int nrec,
                                                        realw* adj_sourcearrays,
                                                        int* ibool,
                                                        int* ispec_selected_rec,
                                                        int* pre_computed_irec,
                                                        int nadj_rec_local) {

  int ispec,iglob;
  int irec,i,j,k;

  int irec_local = blockIdx.x + gridDim.x*blockIdx.y;

  // when nrec > MAXIMUM_GRID_DIM, but mod(nspec_top,2) > 0, we end up with an extra block.
  if(irec_local < nadj_rec_local) {
    irec = pre_computed_irec[irec_local];
    ispec = ispec_selected_rec[irec]-1;

    i = threadIdx.x;
    j = threadIdx.y;
    k = threadIdx.z;
    iglob = ibool[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

    // atomic operations are absolutely necessary for correctness!
    atomicAdd(&accel[3*iglob], adj_sourcearrays[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,
                                                       0,i,j,k,irec_local)]);

    atomicAdd(&accel[3*iglob+1], adj_sourcearrays[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,
                                                         1,i,j,k,irec_local)]);

    atomicAdd(&accel[3*iglob+2], adj_sourcearrays[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,
                                                         2,i,j,k,irec_local)]);
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_add_sources_adjoint_cuda,
              COMPUTE_ADD_SOURCES_ADJOINT_CUDA)(long* Mesh_pointer,
                                                int* nrec,
                                                realw* h_adj_sourcearrays,
                                                int* h_islice_selected_rec,
                                                int* h_ispec_selected_rec,
                                                int* time_index) {

  TRACE("compute_add_sources_adjoint_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // check if anything to do
  if( mp->nadj_rec_local == 0 ) return;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nadj_rec_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y,1);
  dim3 threads(NGLLX,NGLLX,NGLLX);

  // build slice of adj_sourcearrays because full array is *very* large.
  // note: this extracts array values for local adjoint sources at given time step "time_index"
  //          from large adj_sourcearrays array into h_adj_sourcearrays_slice
  int i,j,k;
  int irec_local = 0;
  for(int irec = 0; irec < *nrec; irec++) {
    if(mp->myrank == h_islice_selected_rec[irec]) {
      irec_local++;

      // takes only local sources
      for(k=0;k<NGLLX;k++) {
        for(j=0;j<NGLLX;j++) {
          for(i=0;i<NGLLX;i++) {

            // note: global version uses dimensions
            //          h_adj_sourcearrays(NDIM,NGLLX,NGLLY,NGLLZ,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC)
            mp->h_adj_sourcearrays_slice[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,
                                                0,i,j,k,irec_local-1)]
                    = h_adj_sourcearrays[INDEX6(NDIM,NGLLX,NGLLX,NGLLX,mp->nadj_rec_local,
                                                0,i,j,k,irec_local-1,*time_index-1)];

            mp->h_adj_sourcearrays_slice[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,
                                                1,i,j,k,irec_local-1)]
                    = h_adj_sourcearrays[INDEX6(NDIM,NGLLX,NGLLX,NGLLX,mp->nadj_rec_local,
                                                1,i,j,k,irec_local-1,*time_index-1)];

            mp->h_adj_sourcearrays_slice[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,
                                                2,i,j,k,irec_local-1)]
                    = h_adj_sourcearrays[INDEX6(NDIM,NGLLX,NGLLX,NGLLX,mp->nadj_rec_local,
                                                2,i,j,k,irec_local-1,*time_index-1)];
          }
        }
      }
    }
  }
  // check all local sources were added
  if( irec_local != mp->nadj_rec_local) exit_on_error("irec_local not equal to nadj_rec_local\n");

  // copies extracted array values onto GPU
  print_CUDA_error_if_any(cudaMemcpy(mp->d_adj_sourcearrays, mp->h_adj_sourcearrays_slice,
                                     (mp->nadj_rec_local)*NDIM*NGLL3*sizeof(realw),cudaMemcpyHostToDevice),71000);


  // the irec_local variable needs to be precomputed (as
  // h_pre_comp..), because normally it is in the loop updating accel,
  // and due to how it's incremented, it cannot be parallelized
  compute_add_sources_adjoint_cuda_kernel<<<grid,threads>>>(mp->d_accel_crust_mantle,
                                                            *nrec,
                                                            mp->d_adj_sourcearrays,
                                                            mp->d_ibool_crust_mantle,
                                                            mp->d_ispec_selected_rec,
                                                            mp->d_pre_computed_irec,
                                                            mp->nadj_rec_local);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_add_sources_adjoint_cuda");
#endif
}

