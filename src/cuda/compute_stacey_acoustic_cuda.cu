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

__global__ void compute_stacey_acoustic_kernel(realw* potential_dot_acoustic,
                                               realw* potential_dot_dot_acoustic,
                                               int interface_type,
                                               int num_abs_boundary_faces,
                                               int* abs_boundary_ispec,
                                               int* nkmin_xi, int* nkmin_eta,
                                               int* njmin, int* njmax,
                                               int* nimin, int* nimax,
                                               realw* abs_boundary_jacobian2D,
                                               realw* wgllwgll,
                                               int* ibool,
                                               realw* vpstore,
                                               int SAVE_FORWARD,
                                               realw* b_absorb_potential) {

  int igll = threadIdx.x;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int i,j,k,iglob,ispec;
  realw sn;
  realw jacobianw,fac1;

  // don't compute points outside NGLLSQUARE==NGLL2==25
  // way 2: no further check needed since blocksize = 25
  if( iface < num_abs_boundary_faces){

  //  if(igll<NGLL2 && iface < num_abs_boundary_faces) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    // determines indices i,j,k depending on absorbing boundary type
    switch( interface_type ){
      case 4:
        // xmin
        if( nkmin_xi[INDEX2(2,0,iface)] == 0 || njmin[INDEX2(2,0,iface)] == 0 ) return;

        i = 0; // index -1
        k = (igll/NGLLX);
        j = (igll-k*NGLLX);

        if( k < nkmin_xi[INDEX2(2,0,iface)]-1 || k > NGLLX-1 ) return;
        if( j < njmin[INDEX2(2,0,iface)]-1 || j > njmax[INDEX2(2,0,iface)]-1 ) return;

        fac1 = wgllwgll[k*NGLLX+j];
        break;

      case 5:
        // xmax
        if( nkmin_xi[INDEX2(2,1,iface)] == 0 || njmin[INDEX2(2,1,iface)] == 0 ) return;

        i = NGLLX-1;
        k = (igll/NGLLX);
        j = (igll-k*NGLLX);

        if( k < nkmin_xi[INDEX2(2,1,iface)]-1 || k > NGLLX-1 ) return;
        if( j < njmin[INDEX2(2,1,iface)]-1 || j > njmax[INDEX2(2,1,iface)]-1 ) return;

        fac1 = wgllwgll[k*NGLLX+j];
        break;

      case 6:
        // ymin
        if( nkmin_eta[INDEX2(2,0,iface)] == 0 || nimin[INDEX2(2,0,iface)] == 0 ) return;

        j = 0;
        k = (igll/NGLLX);
        i = (igll-k*NGLLX);

        if( k < nkmin_eta[INDEX2(2,0,iface)]-1 || k > NGLLX-1 ) return;
        if( i < nimin[INDEX2(2,0,iface)]-1 || i > nimax[INDEX2(2,0,iface)]-1 ) return;

        fac1 = wgllwgll[k*NGLLX+i];
        break;

      case 7:
        // ymax
        if( nkmin_eta[INDEX2(2,1,iface)] == 0 || nimin[INDEX2(2,1,iface)] == 0 ) return;

        j = NGLLX-1;
        k = (igll/NGLLX);
        i = (igll-k*NGLLX);

        if( k < nkmin_eta[INDEX2(2,1,iface)]-1 || k > NGLLX-1 ) return;
        if( i < nimin[INDEX2(2,1,iface)]-1 || i > nimax[INDEX2(2,1,iface)]-1 ) return;

        fac1 = wgllwgll[k*NGLLX+i];
        break;

      case 8:
        // zmin
        k = 0;
        j = (igll/NGLLX);
        i = (igll-j*NGLLX);

        if( j < 0 || j > NGLLX-1 ) return;
        if( i < 0 || i > NGLLX-1 ) return;

        fac1 = wgllwgll[j*NGLLX+i];
        break;

    }

    iglob = ibool[INDEX4(5,5,5,i,j,k,ispec)]-1;

    // determines bulk sound speed
    // velocity
    sn = potential_dot_acoustic[iglob] / vpstore[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] ;

    // gets associated, weighted jacobian
    jacobianw = abs_boundary_jacobian2D[INDEX2(NGLL2,igll,iface)]*fac1;

    // Sommerfeld condition
    atomicAdd(&potential_dot_dot_acoustic[iglob],-sn*jacobianw);

    // adjoint simulations
    if( SAVE_FORWARD ){
      // saves boundary values
      b_absorb_potential[INDEX2(NGLL2,igll,iface)] = sn*jacobianw;
    }

  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_stacey_acoustic_cuda,
              COMPUTE_STACEY_ACOUSTIC_CUDA)(long* Mesh_pointer_f,
                                            realw* absorb_potential,
                                            int* itype) {
TRACE("compute_stacey_acoustic_cuda");
  //double start_time = get_time();

  int num_abs_boundary_faces;
  int* d_abs_boundary_ispec;
  realw* d_abs_boundary_jacobian2D;
  realw* d_wgllwgll;
  realw* d_b_absorb_potential;

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  // absorbing boundary type
  int interface_type = *itype;
  switch( interface_type ){
    case 4:
      // xmin
      num_abs_boundary_faces = mp->nspec2D_xmin_outer_core;
      d_abs_boundary_ispec = mp->d_ibelm_xmin_outer_core;
      d_abs_boundary_jacobian2D = mp->d_jacobian2D_xmin_outer_core;
      d_b_absorb_potential = mp->d_absorb_xmin_outer_core;
      d_wgllwgll = mp->d_wgllwgll_yz;
      break;

    case 5:
      // xmax
      num_abs_boundary_faces = mp->nspec2D_xmax_outer_core;
      d_abs_boundary_ispec = mp->d_ibelm_xmax_outer_core;
      d_abs_boundary_jacobian2D = mp->d_jacobian2D_xmax_outer_core;
      d_b_absorb_potential = mp->d_absorb_xmax_outer_core;
      d_wgllwgll = mp->d_wgllwgll_yz;
      break;

    case 6:
      // ymin
      num_abs_boundary_faces = mp->nspec2D_ymin_outer_core;
      d_abs_boundary_ispec = mp->d_ibelm_ymin_outer_core;
      d_abs_boundary_jacobian2D = mp->d_jacobian2D_ymin_outer_core;
      d_b_absorb_potential = mp->d_absorb_ymin_outer_core;
      d_wgllwgll = mp->d_wgllwgll_xz;
      break;

    case 7:
      // ymax
      num_abs_boundary_faces = mp->nspec2D_ymax_outer_core;
      d_abs_boundary_ispec = mp->d_ibelm_ymax_outer_core;
      d_abs_boundary_jacobian2D = mp->d_jacobian2D_ymax_outer_core;
      d_b_absorb_potential = mp->d_absorb_ymax_outer_core;
      d_wgllwgll = mp->d_wgllwgll_xz;
      break;

    case 8:
      // zmin
      num_abs_boundary_faces = mp->nspec2D_zmin_outer_core;
      d_abs_boundary_ispec = mp->d_ibelm_bottom_outer_core;
      d_abs_boundary_jacobian2D = mp->d_jacobian2D_bottom_outer_core;
      d_b_absorb_potential = mp->d_absorb_zmin_outer_core;
      d_wgllwgll = mp->d_wgllwgll_xy;
      break;

    default:
      exit_on_cuda_error("compute_stacey_acoustic_cuda: unknown interface type");
      break;
  }

  // checks if anything to do
  if( num_abs_boundary_faces == 0 ) return;

  // way 1: Elapsed time: 4.385948e-03
  // > NGLLSQUARE==NGLL2==25, but we handle this inside kernel
  //  int blocksize = 32;

  // way 2: Elapsed time: 4.379034e-03
  // > NGLLSQUARE==NGLL2==25, no further check inside kernel
  int blocksize = NGLL2;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(num_abs_boundary_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  compute_stacey_acoustic_kernel<<<grid,threads>>>(mp->d_veloc_outer_core,
                                                   mp->d_accel_outer_core,
                                                   interface_type,
                                                   num_abs_boundary_faces,
                                                   d_abs_boundary_ispec,
                                                   mp->d_nkmin_xi_outer_core,
                                                   mp->d_nkmin_eta_outer_core,
                                                   mp->d_njmin_outer_core,
                                                   mp->d_njmax_outer_core,
                                                   mp->d_nimin_outer_core,
                                                   mp->d_nimax_outer_core,
                                                   d_abs_boundary_jacobian2D,
                                                   d_wgllwgll,
                                                   mp->d_ibool_outer_core,
                                                   mp->d_vp_outer_core,
                                                   mp->save_forward,
                                                   d_b_absorb_potential);

  //  adjoint simulations: stores absorbed wavefield part
  if( mp->save_forward && num_abs_boundary_faces > 0 ){
    // copies array to CPU
    print_CUDA_error_if_any(cudaMemcpy(absorb_potential,d_b_absorb_potential,
                                       NGLL2*num_abs_boundary_faces*sizeof(realw),
                                       cudaMemcpyDeviceToHost),7701);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_stacey_acoustic_kernel");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_stacey_acoustic_backward_kernel(realw* b_potential_dot_dot_acoustic,
                                                        realw* b_absorb_potential,
                                                        int interface_type,
                                                        int num_abs_boundary_faces,
                                                        int* abs_boundary_ispec,
                                                        int* nkmin_xi, int* nkmin_eta,
                                                        int* njmin, int* njmax,
                                                        int* nimin, int* nimax,
                                                        int* ibool) {

  int igll = threadIdx.x;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int i,j,k,iglob,ispec;

  // don't compute points outside NGLLSQUARE==NGLL2==25
  // way 2: no further check needed since blocksize = 25
  if( iface < num_abs_boundary_faces){

  //  if(igll<NGLL2 && iface < num_abs_boundary_faces) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    // determines indices i,j,k depending on absorbing boundary type
    switch( interface_type ){
      case 4:
        // xmin
        if( nkmin_xi[INDEX2(2,0,iface)] == 0 || njmin[INDEX2(2,0,iface)] == 0 ) return;

        i = 0; // index -1
        k = (igll/NGLLX);
        j = (igll-k*NGLLX);

        if( k < nkmin_xi[INDEX2(2,0,iface)]-1 || k > NGLLX-1 ) return;
        if( j < njmin[INDEX2(2,0,iface)]-1 || j > njmax[INDEX2(2,0,iface)]-1 ) return;

        break;

      case 5:
        // xmax
        if( nkmin_xi[INDEX2(2,1,iface)] == 0 || njmin[INDEX2(2,1,iface)] == 0 ) return;

        i = NGLLX-1;
        k = (igll/NGLLX);
        j = (igll-k*NGLLX);

        if( k < nkmin_xi[INDEX2(2,1,iface)]-1 || k > NGLLX-1 ) return;
        if( j < njmin[INDEX2(2,1,iface)]-1 || j > njmax[INDEX2(2,1,iface)]-1 ) return;

        break;

      case 6:
        // ymin
        if( nkmin_eta[INDEX2(2,0,iface)] == 0 || nimin[INDEX2(2,0,iface)] == 0 ) return;

        j = 0;
        k = (igll/NGLLX);
        i = (igll-k*NGLLX);

        if( k < nkmin_eta[INDEX2(2,0,iface)]-1 || k > NGLLX-1 ) return;
        if( i < nimin[INDEX2(2,0,iface)]-1 || i > nimax[INDEX2(2,0,iface)]-1 ) return;

        break;

      case 7:
        // ymax
        if( nkmin_eta[INDEX2(2,1,iface)] == 0 || nimin[INDEX2(2,1,iface)] == 0 ) return;

        j = NGLLX-1;
        k = (igll/NGLLX);
        i = (igll-k*NGLLX);

        if( k < nkmin_eta[INDEX2(2,1,iface)]-1 || k > NGLLX-1 ) return;
        if( i < nimin[INDEX2(2,1,iface)]-1 || i > nimax[INDEX2(2,1,iface)]-1 ) return;

        break;

      case 8:
        // zmin
        k = 0;
        j = (igll/NGLLX);
        i = (igll-j*NGLLX);

        if( j < 0 || j > NGLLX-1 ) return;
        if( i < 0 || i > NGLLX-1 ) return;

        break;

    }

    iglob = ibool[INDEX4(5,5,5,i,j,k,ispec)]-1;

    // Sommerfeld condition
    atomicAdd(&b_potential_dot_dot_acoustic[iglob],-b_absorb_potential[INDEX2(NGLL2,igll,iface)]);
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_stacey_acoustic_backward_cuda,
              COMPUTE_STACEY_ACOUSTIC_BACKWARD_CUDA)(long* Mesh_pointer_f,
                                                     realw* absorb_potential,
                                                     int* itype) {
TRACE("compute_stacey_acoustic_backward_cuda");
  //double start_time = get_time();

  int num_abs_boundary_faces;
  int* d_abs_boundary_ispec;
  realw* d_b_absorb_potential;

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  // absorbing boundary type
  int interface_type = *itype;
  switch( interface_type ){
    case 4:
      // xmin
      num_abs_boundary_faces = mp->nspec2D_xmin_outer_core;
      d_abs_boundary_ispec = mp->d_ibelm_xmin_outer_core;
      d_b_absorb_potential = mp->d_absorb_xmin_outer_core;
      break;

    case 5:
      // xmax
      num_abs_boundary_faces = mp->nspec2D_xmax_outer_core;
      d_abs_boundary_ispec = mp->d_ibelm_xmax_outer_core;
      d_b_absorb_potential = mp->d_absorb_xmax_outer_core;
      break;

    case 6:
      // ymin
      num_abs_boundary_faces = mp->nspec2D_ymin_outer_core;
      d_abs_boundary_ispec = mp->d_ibelm_ymin_outer_core;
      d_b_absorb_potential = mp->d_absorb_ymin_outer_core;
      break;

    case 7:
      // ymax
      num_abs_boundary_faces = mp->nspec2D_ymax_outer_core;
      d_abs_boundary_ispec = mp->d_ibelm_ymax_outer_core;
      d_b_absorb_potential = mp->d_absorb_ymax_outer_core;
      break;

    case 8:
      // zmin
      num_abs_boundary_faces = mp->nspec2D_zmin_outer_core;
      d_abs_boundary_ispec = mp->d_ibelm_bottom_outer_core;
      d_b_absorb_potential = mp->d_absorb_zmin_outer_core;
      break;

    default:
      exit_on_cuda_error("compute_stacey_acoustic_cuda: unknown interface type");
      break;
  }

  // checks if anything to do
  if( num_abs_boundary_faces == 0 ) return;

  // way 1: Elapsed time: 4.385948e-03
  // > NGLLSQUARE==NGLL2==25, but we handle this inside kernel
  //  int blocksize = 32;

  // way 2: Elapsed time: 4.379034e-03
  // > NGLLSQUARE==NGLL2==25, no further check inside kernel
  int blocksize = NGLL2;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(num_abs_boundary_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  //  adjoint simulations: needs absorbing boundary buffer
  if( num_abs_boundary_faces > 0 ){
    // copies array to GPU
    print_CUDA_error_if_any(cudaMemcpy(d_b_absorb_potential,absorb_potential,
                                       NGLL2*num_abs_boundary_faces*sizeof(realw),
                                       cudaMemcpyHostToDevice),7700);
  }

  compute_stacey_acoustic_backward_kernel<<<grid,threads>>>(mp->d_b_accel_outer_core,
                                                            d_b_absorb_potential,
                                                            interface_type,
                                                            num_abs_boundary_faces,
                                                            d_abs_boundary_ispec,
                                                            mp->d_nkmin_xi_outer_core,mp->d_nkmin_eta_outer_core,
                                                            mp->d_njmin_outer_core,mp->d_njmax_outer_core,
                                                            mp->d_nimin_outer_core,mp->d_nimax_outer_core,
                                                            mp->d_ibool_outer_core);


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_stacey_acoustic_backward_kernel");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// undo_attenuation simulation: stacey for backward/reconstructed wavefield

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_stacey_acoustic_undoatt_cuda,
              COMPUTE_STACEY_ACOUSTIC_UNDOATT_CUDA)(long* Mesh_pointer_f,
                                                    int* itype) {
  TRACE("compute_stacey_acoustic_undoatt_cuda");
  //double start_time = get_time();

  int num_abs_boundary_faces;
  int* d_abs_boundary_ispec;
  realw* d_abs_boundary_jacobian2D;
  realw* d_wgllwgll;
  realw* d_b_absorb_potential = NULL;

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if( mp->simulation_type /= 3 ) return;
  if( mp->save_forward ) return;

  // absorbing boundary type
  int interface_type = *itype;
  switch( interface_type ){
    case 4:
      // xmin
      num_abs_boundary_faces = mp->nspec2D_xmin_outer_core;
      d_abs_boundary_ispec = mp->d_ibelm_xmin_outer_core;
      d_abs_boundary_jacobian2D = mp->d_jacobian2D_xmin_outer_core;
      d_wgllwgll = mp->d_wgllwgll_yz;
      break;

    case 5:
      // xmax
      num_abs_boundary_faces = mp->nspec2D_xmax_outer_core;
      d_abs_boundary_ispec = mp->d_ibelm_xmax_outer_core;
      d_abs_boundary_jacobian2D = mp->d_jacobian2D_xmax_outer_core;
      d_wgllwgll = mp->d_wgllwgll_yz;
      break;

    case 6:
      // ymin
      num_abs_boundary_faces = mp->nspec2D_ymin_outer_core;
      d_abs_boundary_ispec = mp->d_ibelm_ymin_outer_core;
      d_abs_boundary_jacobian2D = mp->d_jacobian2D_ymin_outer_core;
      d_wgllwgll = mp->d_wgllwgll_xz;
      break;

    case 7:
      // ymax
      num_abs_boundary_faces = mp->nspec2D_ymax_outer_core;
      d_abs_boundary_ispec = mp->d_ibelm_ymax_outer_core;
      d_abs_boundary_jacobian2D = mp->d_jacobian2D_ymax_outer_core;
      d_wgllwgll = mp->d_wgllwgll_xz;
      break;

    case 8:
      // zmin
      num_abs_boundary_faces = mp->nspec2D_zmin_outer_core;
      d_abs_boundary_ispec = mp->d_ibelm_bottom_outer_core;
      d_abs_boundary_jacobian2D = mp->d_jacobian2D_bottom_outer_core;
      d_wgllwgll = mp->d_wgllwgll_xy;
      break;

    default:
      exit_on_cuda_error("compute_stacey_acoustic_cuda: unknown interface type");
      break;
  }

  // checks if anything to do
  if( num_abs_boundary_faces == 0 ) return;

  // way 1: Elapsed time: 4.385948e-03
  // > NGLLSQUARE==NGLL2==25, but we handle this inside kernel
  //  int blocksize = 32;

  // way 2: Elapsed time: 4.379034e-03
  // > NGLLSQUARE==NGLL2==25, no further check inside kernel
  int blocksize = NGLL2;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(num_abs_boundary_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  compute_stacey_acoustic_kernel<<<grid,threads>>>(mp->d_b_veloc_outer_core,
                                                   mp->d_b_accel_outer_core,
                                                   interface_type,
                                                   num_abs_boundary_faces,
                                                   d_abs_boundary_ispec,
                                                   mp->d_nkmin_xi_outer_core,
                                                   mp->d_nkmin_eta_outer_core,
                                                   mp->d_njmin_outer_core,
                                                   mp->d_njmax_outer_core,
                                                   mp->d_nimin_outer_core,
                                                   mp->d_nimax_outer_core,
                                                   d_abs_boundary_jacobian2D,
                                                   d_wgllwgll,
                                                   mp->d_ibool_outer_core,
                                                   mp->d_vp_outer_core,
                                                   mp->save_forward,
                                                   d_b_absorb_potential);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_stacey_acoustic_undoatt_cuda");
#endif
}

