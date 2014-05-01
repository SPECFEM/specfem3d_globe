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
                                                            realw* d_field,
                                                            int nrec_local) {

// vector fields

  int blockID = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;

  int irec,ispec,iglob;

  if(blockID < nrec_local) {
    irec = number_receiver_global[blockID]-1;
    ispec = ispec_selected_rec[irec]-1;
    iglob = ibool[tx + NGLL3*ispec]-1;

    station_seismo_field[3*NGLL3*blockID + 3*tx+0] = d_field[3*iglob];
    station_seismo_field[3*NGLL3*blockID + 3*tx+1] = d_field[3*iglob+1];
    station_seismo_field[3*NGLL3*blockID + 3*tx+2] = d_field[3*iglob+2];
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void write_seismograms_transfer_strain_from_device_kernel(int* number_receiver_global,
                                                                     int* ispec_selected_rec,
                                                                     int* ibool,
                                                                     realw* station_strain_field,
                                                                     realw* d_field,
                                                                     int nrec_local) {

// scalar fields

  int blockID = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;

  int irec,ispec,iglob;

  if(blockID < nrec_local) {
    irec = number_receiver_global[blockID]-1;
    ispec = ispec_selected_rec[irec]-1;
    iglob = ibool[tx + NGLL3*ispec]-1;

    station_strain_field[NGLL3*blockID + tx] = d_field[iglob];
  }
}

/* ----------------------------------------------------------------------------------------------- */

void write_seismograms_transfer_from_device(Mesh* mp,
                                            realw* d_field,
                                            realw* h_field,
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

  // waits for previous copy call to be finished
  if( GPU_ASYNC_COPY ){
    cudaStreamSynchronize(mp->copy_stream);
  }

  // prepare field transfer array on device
  write_seismograms_transfer_from_device_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_number_receiver_global,
                                                                  d_ispec_selected,
                                                                  mp->d_ibool_crust_mantle,
                                                                  mp->d_station_seismo_field,
                                                                  d_field,
                                                                  mp->nrec_local);

  // copies array to CPU
  if( GPU_ASYNC_COPY ){
    // waits for previous compute call to be finished
    cudaStreamSynchronize(mp->compute_stream);

    // asynchronous copy
    // note: we need to update the host array in a subsequent call to transfer_seismo_from_device_async() routine
    cudaMemcpyAsync(mp->h_station_seismo_field,mp->d_station_seismo_field,
                    3*NGLL3*(mp->nrec_local)*sizeof(realw),
                    cudaMemcpyDeviceToHost,mp->copy_stream);
  }else{
    // synchronous copy
    print_CUDA_error_if_any(cudaMemcpy(mp->h_station_seismo_field,mp->d_station_seismo_field,
                                       3*NGLL3*(mp->nrec_local)*sizeof(realw),cudaMemcpyDeviceToHost),77000);

    // updates array on CPU
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
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("write_seismograms_transfer_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

void write_seismograms_transfer_strain_from_device(Mesh* mp,
                                                   realw* d_field,
                                                   realw* h_field,
                                                   int* number_receiver_global,
                                                   int* d_ispec_selected,
                                                   int* h_ispec_selected,
                                                   int* ibool) {

  TRACE("write_seismograms_transfer_strain_from_device");

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
  write_seismograms_transfer_strain_from_device_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_number_receiver_global,
                                                                         d_ispec_selected,
                                                                         mp->d_ibool_crust_mantle,
                                                                         mp->d_station_strain_field,
                                                                         d_field,
                                                                         mp->nrec_local);

  // copies array to CPU
  // synchronous copy
  print_CUDA_error_if_any(cudaMemcpy(mp->h_station_strain_field,mp->d_station_strain_field,
                                     NGLL3*(mp->nrec_local)*sizeof(realw),cudaMemcpyDeviceToHost),77001);

  // updates host array
  for(irec_local = 0 ; irec_local < mp->nrec_local; irec_local++) {
    irec = number_receiver_global[irec_local] - 1;
    ispec = h_ispec_selected[irec] - 1;
    for(i = 0; i < NGLL3; i++) {
      iglob = ibool[i+NGLL3*ispec] - 1;
      h_field[iglob] = mp->h_station_strain_field[i+irec_local*NGLL3];
    }
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("write_seismograms_transfer_strain_from_device");
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

  // transfers displacement values in receiver elements from GPU to CPU
  switch( mp->simulation_type ){
    case 1:
      // forward simulation
      write_seismograms_transfer_from_device(mp,mp->d_displ_crust_mantle,
                                             displ,
                                             number_receiver_global,
                                             mp->d_ispec_selected_rec,
                                             ispec_selected_rec, ibool);
      break;

    case 2:
      // adjoint simulation
      write_seismograms_transfer_from_device(mp,mp->d_displ_crust_mantle,
                                             displ,
                                             number_receiver_global,
                                             mp->d_ispec_selected_source,
                                             ispec_selected_source, ibool);

      // strain
      write_seismograms_transfer_strain_from_device(mp,mp->d_eps_trace_over_3_crust_mantle,
                                                    eps_trace_over_3,
                                                    number_receiver_global,
                                                    mp->d_ispec_selected_source,
                                                    ispec_selected_source, ibool);
      write_seismograms_transfer_strain_from_device(mp,mp->d_epsilondev_xx_crust_mantle,
                                                    epsilondev_xx,
                                                    number_receiver_global,
                                                    mp->d_ispec_selected_source,
                                                    ispec_selected_source, ibool);
      write_seismograms_transfer_strain_from_device(mp,mp->d_epsilondev_yy_crust_mantle,
                                                    epsilondev_yy,
                                                    number_receiver_global,
                                                    mp->d_ispec_selected_source,
                                                    ispec_selected_source, ibool);
      write_seismograms_transfer_strain_from_device(mp,mp->d_epsilondev_xy_crust_mantle,
                                                    epsilondev_xy,
                                                    number_receiver_global,
                                                    mp->d_ispec_selected_source,
                                                    ispec_selected_source, ibool);
      write_seismograms_transfer_strain_from_device(mp,mp->d_epsilondev_xz_crust_mantle,
                                                    epsilondev_xz,
                                                    number_receiver_global,
                                                    mp->d_ispec_selected_source,
                                                    ispec_selected_source, ibool);
      write_seismograms_transfer_strain_from_device(mp,mp->d_epsilondev_yz_crust_mantle,
                                                    epsilondev_yz,
                                                    number_receiver_global,
                                                    mp->d_ispec_selected_source,
                                                    ispec_selected_source, ibool);

      break;

    case 3:
      // kernel simulation
      write_seismograms_transfer_from_device(mp,mp->d_b_displ_crust_mantle,
                                             b_displ,
                                             number_receiver_global,
                                             mp->d_ispec_selected_rec,
                                             ispec_selected_rec, ibool);
      break;
  }

}


/* ----------------------------------------------------------------------------------------------- */

// data transfer to CPU host

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_seismo_from_device_async,
              TRANSFER_SEISMO_FROM_DEVICE_ASYNC)(long* Mesh_pointer_f,
                                                 realw* displ,
                                                 realw* b_displ,
                                                 int* number_receiver_global,
                                                 int* ispec_selected_rec,
                                                 int* ispec_selected_source,
                                                 int* ibool) {

  TRACE("transfer_seismo_from_device_async");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  int irec,ispec,iglob,i;
  realw* h_field;
  int* h_ispec_selected;

  // checks if anything to do
  if( mp->nrec_local == 0 ) return;

  // checks async-memcpy
  if( GPU_ASYNC_COPY == 0 ){
    exit_on_error("transfer_seismo_from_device_async must be called with GPU_ASYNC_COPY == 1, please check mesh_constants_cuda.h");
  }

  // waits for previous copy call to be finished
  cudaStreamSynchronize(mp->copy_stream);

  // transfers displacements
  // select target array on host
  switch( mp->simulation_type ){
    case 1:
      // forward simulation
      h_field = displ;
      h_ispec_selected = ispec_selected_rec;
      break;

    case 2:
      // adjoint simulation
      h_field = displ;
      h_ispec_selected = ispec_selected_source;
      break;

    case 3:
      // kernel simulation
      h_field = b_displ;
      h_ispec_selected = ispec_selected_rec;
      break;
  }

  // updates corresponding array on CPU
  for(int irec_local = 0 ; irec_local < mp->nrec_local; irec_local++) {
    irec = number_receiver_global[irec_local] - 1;
    ispec = h_ispec_selected[irec] - 1;

    for(i = 0; i < NGLL3; i++) {
      iglob = ibool[i+NGLL3*ispec] - 1;
      h_field[0+3*iglob] = mp->h_station_seismo_field[0+3*i+irec_local*NGLL3*3];
      h_field[1+3*iglob] = mp->h_station_seismo_field[1+3*i+irec_local*NGLL3*3];
      h_field[2+3*iglob] = mp->h_station_seismo_field[2+3*i+irec_local*NGLL3*3];
    }
  }

}


