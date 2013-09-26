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

#include "config.h"
#include "mesh_constants_cuda.h"


/* ----------------------------------------------------------------------------------------------- */

// ELASTIC simulations

/* ----------------------------------------------------------------------------------------------- */

__global__ void write_seismograms_transfer_from_device_kernel(int* number_receiver_global,
                                                            int* ispec_selected_rec,
                                                            int* ibool,
                                                            realw* station_seismo_field,
                                                            realw* desired_field,
                                                            int nrec_local) {

// vector fields

  int blockID = blockIdx.x + blockIdx.y*gridDim.x;

  if(blockID < nrec_local) {
    int irec = number_receiver_global[blockID]-1;
    int ispec = ispec_selected_rec[irec]-1;
    int iglob = ibool[threadIdx.x + NGLL3*ispec]-1;

    station_seismo_field[3*NGLL3*blockID + 3*threadIdx.x+0] = desired_field[3*iglob];
    station_seismo_field[3*NGLL3*blockID + 3*threadIdx.x+1] = desired_field[3*iglob+1];
    station_seismo_field[3*NGLL3*blockID + 3*threadIdx.x+2] = desired_field[3*iglob+2];
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void write_seismograms_transfer_scalar_from_device_kernel(int* number_receiver_global,
                                                                     int* ispec_selected_rec,
                                                                     int* ibool,
                                                                     realw* station_seismo_field,
                                                                     realw* desired_field,
                                                                     int nrec_local) {

// scalar fields

  int blockID = blockIdx.x + blockIdx.y*gridDim.x;

  if(blockID < nrec_local) {
    int irec = number_receiver_global[blockID]-1;
    int ispec = ispec_selected_rec[irec]-1;
    int iglob = ibool[threadIdx.x + NGLL3*ispec]-1;

    station_seismo_field[NGLL3*blockID + threadIdx.x] = desired_field[iglob];
  }
}

/* ----------------------------------------------------------------------------------------------- */

void write_seismograms_transfer_from_device(Mesh* mp, realw* d_field,realw* h_field,
                                          int* number_receiver_global,
                                          int* d_ispec_selected,
                                          int* h_ispec_selected,
                                          int* ibool) {

TRACE("write_seismograms_transfer_from_device");

  int irec_local,irec;
  int ispec,iglob,i;

  // checks if anything to do
  if( mp->nrec_local == 0 ) return;

  int blocksize = NGLL3;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nrec_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // prepare field transfer array on device
  write_seismograms_transfer_from_device_kernel<<<grid,threads>>>(mp->d_number_receiver_global,
                                                                  d_ispec_selected,
                                                                  mp->d_ibool_crust_mantle,
                                                                  mp->d_station_seismo_field,
                                                                  d_field,
                                                                  mp->nrec_local);

  // copies array to CPU
  cudaMemcpy(mp->h_station_seismo_field,mp->d_station_seismo_field,
             3*NGLL3*(mp->nrec_local)*sizeof(realw),cudaMemcpyDeviceToHost);

  for(irec_local = 0 ; irec_local < mp->nrec_local; irec_local++) {
    irec = number_receiver_global[irec_local] - 1;
    ispec = h_ispec_selected[irec] - 1;

    for(i = 0; i < NGLL3; i++) {
      iglob = ibool[i+NGLL3*ispec] - 1;
      h_field[0+3*iglob] = mp->h_station_seismo_field[0+3*i+irec_local*NGLL3*3];
      h_field[1+3*iglob] = mp->h_station_seismo_field[1+3*i+irec_local*NGLL3*3];
      h_field[2+3*iglob] = mp->h_station_seismo_field[2+3*i+irec_local*NGLL3*3];
    }
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("write_seismograms_transfer_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

void write_seismograms_transfer_scalar_from_device(Mesh* mp,
                                                   realw* d_field,realw* h_field,
                                                   int* number_receiver_global,
                                                   int* d_ispec_selected,
                                                   int* h_ispec_selected,
                                                   int* ibool) {

  TRACE("write_seismograms_transfer_scalar_from_device");

  int irec_local,irec;
  int ispec,iglob,i;

  // checks if anything to do
  if( mp->nrec_local == 0 ) return;

  int blocksize = NGLL3;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nrec_local,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // prepare field transfer array on device
  write_seismograms_transfer_scalar_from_device_kernel<<<grid,threads>>>(mp->d_number_receiver_global,
                                                                         d_ispec_selected,
                                                                         mp->d_ibool_crust_mantle,
                                                                         mp->d_station_seismo_field,
                                                                         d_field,
                                                                         mp->nrec_local);

  // copies array to CPU
  cudaMemcpy(mp->h_station_seismo_field,mp->d_station_seismo_field,
             NGLL3*(mp->nrec_local)*sizeof(realw),cudaMemcpyDeviceToHost);

  for(irec_local = 0 ; irec_local < mp->nrec_local; irec_local++) {
    irec = number_receiver_global[irec_local] - 1;
    ispec = h_ispec_selected[irec] - 1;
    for(i = 0; i < NGLL3; i++) {
      iglob = ibool[i+NGLL3*ispec] - 1;
      h_field[iglob] = mp->h_station_seismo_field[i+irec_local*NGLL3];
    }
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("write_seismograms_transfer_scalar_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(write_seismograms_transfer_cuda,
              WRITE_SEISMOGRAMS_TRANSFER_CUDA)(long* Mesh_pointer_f,
                                               realw* displ,
                                               realw* b_displ,
                                               realw* eps_trace_over_3,
                                               realw* epsilondev_xx,
                                               realw* epsilondev_yy,
                                               realw* epsilondev_xy,
                                               realw* epsilondev_xz,
                                               realw* epsilondev_yz,
                                               int* number_receiver_global,
                                               int* ispec_selected_rec,
                                               int* ispec_selected_source,
                                               int* ibool) {
TRACE("write_seismograms_transfer_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  // checks if anything to do
  if( mp->nrec_local == 0 ) return;

  switch( mp->simulation_type ){
    case 1:
      write_seismograms_transfer_from_device(mp,mp->d_displ_crust_mantle,
                                             displ,
                                             number_receiver_global,
                                             mp->d_ispec_selected_rec,
                                             ispec_selected_rec, ibool);
      break;

    case 2:
      write_seismograms_transfer_from_device(mp,mp->d_displ_crust_mantle,
                                             displ,
                                             number_receiver_global,
                                             mp->d_ispec_selected_source,
                                             ispec_selected_source, ibool);

      write_seismograms_transfer_scalar_from_device(mp,mp->d_eps_trace_over_3_crust_mantle,
                                                    eps_trace_over_3,
                                                    number_receiver_global,
                                                    mp->d_ispec_selected_source,
                                                    ispec_selected_source, ibool);
      write_seismograms_transfer_scalar_from_device(mp,mp->d_epsilondev_xx_crust_mantle,
                                                    epsilondev_xx,
                                                    number_receiver_global,
                                                    mp->d_ispec_selected_source,
                                                    ispec_selected_source, ibool);
      write_seismograms_transfer_scalar_from_device(mp,mp->d_epsilondev_yy_crust_mantle,
                                                    epsilondev_yy,
                                                    number_receiver_global,
                                                    mp->d_ispec_selected_source,
                                                    ispec_selected_source, ibool);
      write_seismograms_transfer_scalar_from_device(mp,mp->d_epsilondev_xy_crust_mantle,
                                                    epsilondev_xy,
                                                    number_receiver_global,
                                                    mp->d_ispec_selected_source,
                                                    ispec_selected_source, ibool);
      write_seismograms_transfer_scalar_from_device(mp,mp->d_epsilondev_xz_crust_mantle,
                                                    epsilondev_xz,
                                                    number_receiver_global,
                                                    mp->d_ispec_selected_source,
                                                    ispec_selected_source, ibool);
      write_seismograms_transfer_scalar_from_device(mp,mp->d_epsilondev_yz_crust_mantle,
                                                    epsilondev_yz,
                                                    number_receiver_global,
                                                    mp->d_ispec_selected_source,
                                                    ispec_selected_source, ibool);
      break;

    case 3:
      write_seismograms_transfer_from_device(mp,mp->d_b_displ_crust_mantle,b_displ,
                                                 number_receiver_global,
                                                 mp->d_ispec_selected_rec,
                                                 ispec_selected_rec, ibool);
      break;
  }
}

