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

__global__ void compute_stacey_elastic_kernel(realw* veloc,
                                              realw* accel,
                                              int interface_type,
                                              int num_abs_boundary_faces,
                                              int* abs_boundary_ispec,
                                              int* nkmin_xi, int* nkmin_eta,
                                              int* njmin, int* njmax,
                                              int* nimin, int* nimax,
                                              realw* abs_boundary_normal,
                                              realw* abs_boundary_jacobian2D,
                                              realw* wgllwgll,
                                              int* ibool,
                                              realw* rho_vp,
                                              realw* rho_vs,
                                              int SAVE_FORWARD,
                                              realw* b_absorb_field) {

  int igll = threadIdx.x; // tx
  int iface = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int i,j,k,iglob,ispec;
  realw vx,vy,vz,vn;
  realw nx,ny,nz;
  realw rho_vp_temp,rho_vs_temp;
  realw tx,ty,tz;
  realw jacobianw;
  realw fac1;

  // don't compute surface faces outside of range
  // and don't compute points outside NGLLSQUARE==NGLL2==25
  //if(igll < NGLL2 && iface < num_abs_boundary_faces) {

  // way 2: only check face, no further check needed since blocksize = 25
  if( iface < num_abs_boundary_faces){

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    // determines indices i,j,k depending on absorbing boundary type
    switch( interface_type ){
      case 0:
        // xmin
        if( nkmin_xi[INDEX2(2,0,iface)] == 0 || njmin[INDEX2(2,0,iface)] == 0 ) return;

        i = 0; // index -1
        k = (igll/NGLLX);
        j = (igll-k*NGLLX);

        if( k < nkmin_xi[INDEX2(2,0,iface)]-1 || k > NGLLX-1 ) return;
        if( j < njmin[INDEX2(2,0,iface)]-1 || j > NGLLX-1 ) return;

        fac1 = wgllwgll[k*NGLLX+j];
        break;

      case 1:
        // xmax
        if( nkmin_xi[INDEX2(2,1,iface)] == 0 || njmin[INDEX2(2,1,iface)] == 0 ) return;

        i = NGLLX-1;
        k = (igll/NGLLX);
        j = (igll-k*NGLLX);

        if( k < nkmin_xi[INDEX2(2,1,iface)]-1 || k > NGLLX-1 ) return;
        if( j < njmin[INDEX2(2,1,iface)]-1 || j > njmax[INDEX2(2,1,iface)]-1 ) return;

        fac1 = wgllwgll[k*NGLLX+j];
        break;

      case 2:
        // ymin
        if( nkmin_eta[INDEX2(2,0,iface)] == 0 || nimin[INDEX2(2,0,iface)] == 0 ) return;

        j = 0;
        k = (igll/NGLLX);
        i = (igll-k*NGLLX);

        if( k < nkmin_eta[INDEX2(2,0,iface)]-1 || k > NGLLX-1 ) return;
        if( i < nimin[INDEX2(2,0,iface)]-1 || i > nimax[INDEX2(2,0,iface)]-1 ) return;

        fac1 = wgllwgll[k*NGLLX+i];
        break;

      case 3:
        // ymax
        if( nkmin_eta[INDEX2(2,1,iface)] == 0 || nimin[INDEX2(2,1,iface)] == 0 ) return;

        j = NGLLX-1;
        k = (igll/NGLLX);
        i = (igll-k*NGLLX);

        if( k < nkmin_eta[INDEX2(2,1,iface)]-1 || k > NGLLX-1 ) return;
        if( i < nimin[INDEX2(2,1,iface)]-1 || i > nimax[INDEX2(2,1,iface)]-1 ) return;

        fac1 = wgllwgll[k*NGLLX+i];
        break;
    }

    iglob = ibool[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

    // gets associated velocity
    vx = veloc[iglob*3+0];
    vy = veloc[iglob*3+1];
    vz = veloc[iglob*3+2];

    // gets associated normal
    nx = abs_boundary_normal[INDEX3(NDIM,NGLL2,0,igll,iface)];
    ny = abs_boundary_normal[INDEX3(NDIM,NGLL2,1,igll,iface)];
    nz = abs_boundary_normal[INDEX3(NDIM,NGLL2,2,igll,iface)];

    // // velocity component in normal direction (normal points out of element)
    vn = vx*nx + vy*ny + vz*nz;

    rho_vp_temp = rho_vp[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];
    rho_vs_temp = rho_vs[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];

    tx = rho_vp_temp*vn*nx + rho_vs_temp*(vx-vn*nx);
    ty = rho_vp_temp*vn*ny + rho_vs_temp*(vy-vn*ny);
    tz = rho_vp_temp*vn*nz + rho_vs_temp*(vz-vn*nz);

    jacobianw = abs_boundary_jacobian2D[INDEX2(NGLL2,igll,iface)]*fac1;

    atomicAdd(&accel[iglob*3],-tx*jacobianw);
    atomicAdd(&accel[iglob*3+1],-ty*jacobianw);
    atomicAdd(&accel[iglob*3+2],-tz*jacobianw);

    if( SAVE_FORWARD ){
      b_absorb_field[INDEX3(NDIM,NGLL2,0,igll,iface)] = tx*jacobianw;
      b_absorb_field[INDEX3(NDIM,NGLL2,1,igll,iface)] = ty*jacobianw;
      b_absorb_field[INDEX3(NDIM,NGLL2,2,igll,iface)] = tz*jacobianw;
    } // SIMULATION_TYPE

  } // num_abs_boundary_faces
}

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_stacey_elastic_cuda,
              COMPUTE_STACEY_ELASTIC_CUDA)(long* Mesh_pointer_f,
                                           realw* absorb_field,
                                           int* itype) {

TRACE("compute_stacey_elastic_cuda");

  int num_abs_boundary_faces;
  int* d_abs_boundary_ispec;
  realw* d_abs_boundary_normal;
  realw* d_abs_boundary_jacobian2D;
  realw* d_wgllwgll;
  realw* d_b_absorb_field;

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  // absorbing boundary type
  int interface_type = *itype;
  switch( interface_type ){
    case 0:
      // xmin
      num_abs_boundary_faces = mp->nspec2D_xmin_crust_mantle;
      d_abs_boundary_ispec = mp->d_ibelm_xmin_crust_mantle;
      d_abs_boundary_normal = mp->d_normal_xmin_crust_mantle;
      d_abs_boundary_jacobian2D = mp->d_jacobian2D_xmin_crust_mantle;
      d_b_absorb_field = mp->d_absorb_xmin_crust_mantle;
      d_wgllwgll = mp->d_wgllwgll_yz;
      break;

    case 1:
      // xmax
      num_abs_boundary_faces = mp->nspec2D_xmax_crust_mantle;
      d_abs_boundary_ispec = mp->d_ibelm_xmax_crust_mantle;
      d_abs_boundary_normal = mp->d_normal_xmax_crust_mantle;
      d_abs_boundary_jacobian2D = mp->d_jacobian2D_xmax_crust_mantle;
      d_b_absorb_field = mp->d_absorb_xmax_crust_mantle;
      d_wgllwgll = mp->d_wgllwgll_yz;
      break;

    case 2:
      // ymin
      num_abs_boundary_faces = mp->nspec2D_ymin_crust_mantle;
      d_abs_boundary_ispec = mp->d_ibelm_ymin_crust_mantle;
      d_abs_boundary_normal = mp->d_normal_ymin_crust_mantle;
      d_abs_boundary_jacobian2D = mp->d_jacobian2D_ymin_crust_mantle;
      d_b_absorb_field = mp->d_absorb_ymin_crust_mantle;
      d_wgllwgll = mp->d_wgllwgll_xz;
      break;

    case 3:
      // ymax
      num_abs_boundary_faces = mp->nspec2D_ymax_crust_mantle;
      d_abs_boundary_ispec = mp->d_ibelm_ymax_crust_mantle;
      d_abs_boundary_normal = mp->d_normal_ymax_crust_mantle;
      d_abs_boundary_jacobian2D = mp->d_jacobian2D_ymax_crust_mantle;
      d_b_absorb_field = mp->d_absorb_ymax_crust_mantle;
      d_wgllwgll = mp->d_wgllwgll_xz;
      break;

    default:
      exit_on_cuda_error("compute_stacey_elastic_cuda: unknown interface type");
      break;
  }

  // checks if anything to do
  if( num_abs_boundary_faces == 0 ) return;

  // way 1
  // > NGLLSQUARE==NGLL2==25, but we handle this inside kernel
  //int blocksize = 32;

  // way 2: seems sligthly faster
  // > NGLLSQUARE==NGLL2==25, no further check inside kernel
  int blocksize = NGLL2;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(num_abs_boundary_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // absorbing boundary contributions
  compute_stacey_elastic_kernel<<<grid,threads>>>(mp->d_veloc_crust_mantle,
                                                  mp->d_accel_crust_mantle,
                                                  interface_type,
                                                  num_abs_boundary_faces,
                                                  d_abs_boundary_ispec,
                                                  mp->d_nkmin_xi_crust_mantle,
                                                  mp->d_nkmin_eta_crust_mantle,
                                                  mp->d_njmin_crust_mantle,
                                                  mp->d_njmax_crust_mantle,
                                                  mp->d_nimin_crust_mantle,
                                                  mp->d_nimax_crust_mantle,
                                                  d_abs_boundary_normal,
                                                  d_abs_boundary_jacobian2D,
                                                  d_wgllwgll,
                                                  mp->d_ibool_crust_mantle,
                                                  mp->d_rho_vp_crust_mantle,
                                                  mp->d_rho_vs_crust_mantle,
                                                  mp->save_forward,
                                                  d_b_absorb_field);


  // adjoint simulations: stores absorbed wavefield part
  if(mp->save_forward && num_abs_boundary_faces > 0 ) {
    // copies array to CPU
    print_CUDA_error_if_any(cudaMemcpy(absorb_field,d_b_absorb_field,
                            NDIM*NGLL2*num_abs_boundary_faces*sizeof(realw),cudaMemcpyDeviceToHost),7701);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_stacey_elastic_cuda");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// backward/reconstructed wavefields

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_stacey_elastic_backward_kernel(realw* b_accel,
                                                       realw* b_absorb_field,
                                                       int interface_type,
                                                       int num_abs_boundary_faces,
                                                       int* abs_boundary_ispec,
                                                       int* nkmin_xi, int* nkmin_eta,
                                                       int* njmin, int* njmax,
                                                       int* nimin, int* nimax,
                                                       int* ibool) {

  int igll = threadIdx.x; // tx
  int iface = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int i,j,k,iglob,ispec;

  // don't compute surface faces outside of range
  // and don't compute points outside NGLLSQUARE==NGLL2==25
  //if(igll < NGLL2 && iface < num_abs_boundary_faces) {

  // way 2: only check face, no further check needed since blocksize = 25
  if( iface < num_abs_boundary_faces){

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    // determines indices i,j,k depending on absorbing boundary type
    switch( interface_type ){
      case 0:
        // xmin
        if( nkmin_xi[INDEX2(2,0,iface)] == 0 || njmin[INDEX2(2,0,iface)] == 0 ) return;

        i = 0; // index -1
        k = (igll/NGLLX);
        j = (igll-k*NGLLX);

        if( k < nkmin_xi[INDEX2(2,0,iface)]-1 || k > NGLLX-1 ) return;
        if( j < njmin[INDEX2(2,0,iface)]-1 || j > NGLLX-1 ) return;

        break;

      case 1:
        // xmax
        if( nkmin_xi[INDEX2(2,1,iface)] == 0 || njmin[INDEX2(2,1,iface)] == 0 ) return;

        i = NGLLX-1;
        k = (igll/NGLLX);
        j = (igll-k*NGLLX);

        if( k < nkmin_xi[INDEX2(2,1,iface)]-1 || k > NGLLX-1 ) return;
        if( j < njmin[INDEX2(2,1,iface)]-1 || j > njmax[INDEX2(2,1,iface)]-1 ) return;

        break;

      case 2:
        // ymin
        if( nkmin_eta[INDEX2(2,0,iface)] == 0 || nimin[INDEX2(2,0,iface)] == 0 ) return;

        j = 0;
        k = (igll/NGLLX);
        i = (igll-k*NGLLX);

        if( k < nkmin_eta[INDEX2(2,0,iface)]-1 || k > NGLLX-1 ) return;
        if( i < nimin[INDEX2(2,0,iface)]-1 || i > nimax[INDEX2(2,0,iface)]-1 ) return;

        break;

      case 3:
        // ymax
        if( nkmin_eta[INDEX2(2,1,iface)] == 0 || nimin[INDEX2(2,1,iface)] == 0 ) return;

        j = NGLLX-1;
        k = (igll/NGLLX);
        i = (igll-k*NGLLX);

        if( k < nkmin_eta[INDEX2(2,1,iface)]-1 || k > NGLLX-1 ) return;
        if( i < nimin[INDEX2(2,1,iface)]-1 || i > nimax[INDEX2(2,1,iface)]-1 ) return;

        break;
    }

    iglob = ibool[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

    atomicAdd(&b_accel[iglob*3  ],-b_absorb_field[INDEX3(NDIM,NGLL2,0,igll,iface)]);
    atomicAdd(&b_accel[iglob*3+1],-b_absorb_field[INDEX3(NDIM,NGLL2,1,igll,iface)]);
    atomicAdd(&b_accel[iglob*3+2],-b_absorb_field[INDEX3(NDIM,NGLL2,2,igll,iface)]);

  } // num_abs_boundary_faces
}


/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_stacey_elastic_backward_cuda,
              COMPUTE_STACEY_ELASTIC_BACKWARD_CUDA)(long* Mesh_pointer_f,
                                                    realw* absorb_field,
                                                    int* itype) {

  TRACE("compute_stacey_elastic_backward_cuda");

  int num_abs_boundary_faces;
  int* d_abs_boundary_ispec;
  realw* d_b_absorb_field;

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  // absorbing boundary type
  int interface_type = *itype;
  switch( interface_type ){
    case 0:
      // xmin
      num_abs_boundary_faces = mp->nspec2D_xmin_crust_mantle;
      d_abs_boundary_ispec = mp->d_ibelm_xmin_crust_mantle;
      d_b_absorb_field = mp->d_absorb_xmin_crust_mantle;
      break;

    case 1:
      // xmax
      num_abs_boundary_faces = mp->nspec2D_xmax_crust_mantle;
      d_abs_boundary_ispec = mp->d_ibelm_xmax_crust_mantle;
      d_b_absorb_field = mp->d_absorb_xmax_crust_mantle;
      break;

    case 2:
      // ymin
      num_abs_boundary_faces = mp->nspec2D_ymin_crust_mantle;
      d_abs_boundary_ispec = mp->d_ibelm_ymin_crust_mantle;
      d_b_absorb_field = mp->d_absorb_ymin_crust_mantle;
      break;

    case 3:
      // ymax
      num_abs_boundary_faces = mp->nspec2D_ymax_crust_mantle;
      d_abs_boundary_ispec = mp->d_ibelm_ymax_crust_mantle;
      d_b_absorb_field = mp->d_absorb_ymax_crust_mantle;
      break;

    default:
      exit_on_cuda_error("compute_stacey_elastic_cuda: unknown interface type");
      break;
  }

  // checks if anything to do
  if( num_abs_boundary_faces == 0 ) return;

  // way 1
  // > NGLLSQUARE==NGLL2==25, but we handle this inside kernel
  //int blocksize = 32;

  // way 2: seems sligthly faster
  // > NGLLSQUARE==NGLL2==25, no further check inside kernel
  int blocksize = NGLL2;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(num_abs_boundary_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // adjoint simulations: needs absorbing boundary buffer
  // copies array to GPU
  print_CUDA_error_if_any(cudaMemcpy(d_b_absorb_field,absorb_field,
                                     NDIM*NGLL2*num_abs_boundary_faces*sizeof(realw),cudaMemcpyHostToDevice),7700);

  // absorbing boundary contributions
  compute_stacey_elastic_backward_kernel<<<grid,threads>>>(mp->d_b_accel_crust_mantle,
                                                           d_b_absorb_field,
                                                           interface_type,
                                                           num_abs_boundary_faces,
                                                           d_abs_boundary_ispec,
                                                           mp->d_nkmin_xi_crust_mantle,mp->d_nkmin_eta_crust_mantle,
                                                           mp->d_njmin_crust_mantle,mp->d_njmax_crust_mantle,
                                                           mp->d_nimin_crust_mantle,mp->d_nimax_crust_mantle,
                                                           mp->d_ibool_crust_mantle);


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_stacey_elastic_backward_cuda");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// undo_attenuation simulation: stacey for backward/reconstructed wavefield

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_stacey_elastic_undoatt_cuda,
              COMPUTE_STACEY_ELASTIC_UNDOATT_CUDA)(long* Mesh_pointer_f,
                                                   int* itype) {

  TRACE("compute_stacey_elastic_undoatt_cuda");

  int num_abs_boundary_faces;
  int* d_abs_boundary_ispec;
  realw* d_abs_boundary_normal;
  realw* d_abs_boundary_jacobian2D;
  realw* d_wgllwgll;
  realw* d_b_absorb_field = NULL;

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if( mp->simulation_type /= 3 ) return;
  if( mp->save_forward ) return;

  // absorbing boundary type
  int interface_type = *itype;
  switch( interface_type ){
    case 0:
      // xmin
      num_abs_boundary_faces = mp->nspec2D_xmin_crust_mantle;
      d_abs_boundary_ispec = mp->d_ibelm_xmin_crust_mantle;
      d_abs_boundary_normal = mp->d_normal_xmin_crust_mantle;
      d_abs_boundary_jacobian2D = mp->d_jacobian2D_xmin_crust_mantle;
      d_wgllwgll = mp->d_wgllwgll_yz;
      break;

    case 1:
      // xmax
      num_abs_boundary_faces = mp->nspec2D_xmax_crust_mantle;
      d_abs_boundary_ispec = mp->d_ibelm_xmax_crust_mantle;
      d_abs_boundary_normal = mp->d_normal_xmax_crust_mantle;
      d_abs_boundary_jacobian2D = mp->d_jacobian2D_xmax_crust_mantle;
      d_wgllwgll = mp->d_wgllwgll_yz;
      break;

    case 2:
      // ymin
      num_abs_boundary_faces = mp->nspec2D_ymin_crust_mantle;
      d_abs_boundary_ispec = mp->d_ibelm_ymin_crust_mantle;
      d_abs_boundary_normal = mp->d_normal_ymin_crust_mantle;
      d_abs_boundary_jacobian2D = mp->d_jacobian2D_ymin_crust_mantle;
      d_wgllwgll = mp->d_wgllwgll_xz;
      break;

    case 3:
      // ymax
      num_abs_boundary_faces = mp->nspec2D_ymax_crust_mantle;
      d_abs_boundary_ispec = mp->d_ibelm_ymax_crust_mantle;
      d_abs_boundary_normal = mp->d_normal_ymax_crust_mantle;
      d_abs_boundary_jacobian2D = mp->d_jacobian2D_ymax_crust_mantle;
      d_wgllwgll = mp->d_wgllwgll_xz;
      break;

    default:
      exit_on_cuda_error("compute_stacey_elastic_undoatt_cuda: unknown interface type");
      break;
  }

  // checks if anything to do
  if( num_abs_boundary_faces == 0 ) return;

  // way 1
  // > NGLLSQUARE==NGLL2==25, but we handle this inside kernel
  //int blocksize = 32;

  // way 2: seems sligthly faster
  // > NGLLSQUARE==NGLL2==25, no further check inside kernel
  int blocksize = NGLL2;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(num_abs_boundary_faces,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // absorbing boundary contributions
  compute_stacey_elastic_kernel<<<grid,threads>>>(mp->d_b_veloc_crust_mantle,
                                                  mp->d_b_accel_crust_mantle,
                                                  interface_type,
                                                  num_abs_boundary_faces,
                                                  d_abs_boundary_ispec,
                                                  mp->d_nkmin_xi_crust_mantle,
                                                  mp->d_nkmin_eta_crust_mantle,
                                                  mp->d_njmin_crust_mantle,
                                                  mp->d_njmax_crust_mantle,
                                                  mp->d_nimin_crust_mantle,
                                                  mp->d_nimax_crust_mantle,
                                                  d_abs_boundary_normal,
                                                  d_abs_boundary_jacobian2D,
                                                  d_wgllwgll,
                                                  mp->d_ibool_crust_mantle,
                                                  mp->d_rho_vp_crust_mantle,
                                                  mp->d_rho_vs_crust_mantle,
                                                  mp->save_forward,
                                                  d_b_absorb_field);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_stacey_elastic_undoatt_cuda");
#endif
}

