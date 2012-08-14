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

// ACOUSTIC - ELASTIC coupling

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_coupling_fluid_CMB_kernel(realw* displ_crust_mantle,
                                                  realw* accel_outer_core,
                                                  int* ibool_crust_mantle,
                                                  int* ibelm_bottom_crust_mantle,
                                                  realw* normal_top_outer_core,
                                                  realw* jacobian2D_top_outer_core,
                                                  realw* wgllwgll_xy,
                                                  int* ibool_outer_core,
                                                  int* ibelm_top_outer_core,
                                                  int NSPEC2D_TOP_OC) {

  int i = threadIdx.x;
  int j = threadIdx.y;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int k,k_corresp,iglob_cm,iglob_oc,ispec,ispec_selected;
  realw displ_x,displ_y,displ_z,displ_n;
  realw nx,ny,nz;
  realw weight;

  // for surfaces elements exactly at the top of the outer core (crust mantle bottom)
  if( iface < NSPEC2D_TOP_OC ){

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = ibelm_top_outer_core[iface] - 1;
    ispec_selected = ibelm_bottom_crust_mantle[iface] - 1;

    // only for DOFs exactly on the CMB (top of these elements)
    k = NGLLX - 1;
    // get displacement on the solid side using pointwise matching
    k_corresp = 0;

    // get global point number
    // corresponding points are located at the bottom of the mantle
    iglob_cm = ibool_crust_mantle[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k_corresp,ispec_selected)] - 1;

    // elastic displacement on global point
    displ_x = displ_crust_mantle[iglob_cm*3]; // (1,iglob)
    displ_y = displ_crust_mantle[iglob_cm*3+1]; // (2,iglob)
    displ_z = displ_crust_mantle[iglob_cm*3+2]; // (3,iglob)

    // get normal on the CMB
    nx = normal_top_outer_core[INDEX4(NDIM,NGLLX,NGLLX,0,i,j,iface)]; // (1,i,j,iface)
    ny = normal_top_outer_core[INDEX4(NDIM,NGLLX,NGLLX,1,i,j,iface)]; // (2,i,j,iface)
    nz = normal_top_outer_core[INDEX4(NDIM,NGLLX,NGLLX,2,i,j,iface)]; // (3,i,j,iface)

    // calculates displacement component along normal
    // (normal points outwards of acoustic element)
    displ_n = displ_x*nx + displ_y*ny + displ_z*nz;

    // formulation with generalized potential: gets associated, weighted jacobian
    weight = jacobian2D_top_outer_core[INDEX3(NGLLX,NGLLX,i,j,iface)]*wgllwgll_xy[INDEX2(NGLLX,i,j)];

    // get global point number
    iglob_oc = ibool_outer_core[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] - 1;

    // update fluid acceleration/pressure
    atomicAdd(&accel_outer_core[iglob_oc],+weight*displ_n);

  // }
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_coupling_fluid_ICB_kernel(realw* displ_inner_core,
                                                  realw* accel_outer_core,
                                                  int* ibool_inner_core,
                                                  int* ibelm_top_inner_core,
                                                  realw* normal_bottom_outer_core,
                                                  realw* jacobian2D_bottom_outer_core,
                                                  realw* wgllwgll_xy,
                                                  int* ibool_outer_core,
                                                  int* ibelm_bottom_outer_core,
                                                  int NSPEC2D_BOTTOM_OC) {

  int i = threadIdx.x;
  int j = threadIdx.y;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int k,k_corresp,iglob_ic,iglob_oc,ispec,ispec_selected;
  realw displ_x,displ_y,displ_z,displ_n;
  realw nx,ny,nz;
  realw weight;

  // for surfaces elements exactly at the bottom of the outer core (inner core top)
  if( iface < NSPEC2D_BOTTOM_OC ){

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = ibelm_bottom_outer_core[iface] - 1;
    ispec_selected = ibelm_top_inner_core[iface] - 1;

    // only for DOFs exactly on the ICB (bottom of these elements)
    k = 0;
    // get displacement on the solid side using pointwise matching
    k_corresp = NGLLX - 1;

    // get global point number
    // corresponding points are located at the bottom of the mantle
    iglob_ic = ibool_inner_core[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k_corresp,ispec_selected)] - 1;

    // elastic displacement on global point
    displ_x = displ_inner_core[iglob_ic*3]; // (1,iglob)
    displ_y = displ_inner_core[iglob_ic*3+1]; // (2,iglob)
    displ_z = displ_inner_core[iglob_ic*3+2]; // (3,iglob)

    // get normal on the ICB
    nx = normal_bottom_outer_core[INDEX4(NDIM,NGLLX,NGLLX,0,i,j,iface)]; // (1,i,j,iface)
    ny = normal_bottom_outer_core[INDEX4(NDIM,NGLLX,NGLLX,1,i,j,iface)]; // (2,i,j,iface)
    nz = normal_bottom_outer_core[INDEX4(NDIM,NGLLX,NGLLX,2,i,j,iface)]; // (3,i,j,iface)

    // calculates displacement component along normal
    // (normal points outwards of acoustic element)
    displ_n = displ_x*nx + displ_y*ny + displ_z*nz;

    // formulation with generalized potential: gets associated, weighted jacobian
    weight = jacobian2D_bottom_outer_core[INDEX3(NGLLX,NGLLX,i,j,iface)]*wgllwgll_xy[INDEX2(NGLLX,i,j)];

    // get global point number
    iglob_oc = ibool_outer_core[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] - 1;

    // update fluid acceleration/pressure
    atomicAdd(&accel_outer_core[iglob_oc],+ weight*displ_n);

  // }
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_coupling_fluid_cmb_cuda,
              COMPUTE_COUPLING_FLUID_CMB_CUDA)(long* Mesh_pointer_f) {

  TRACE("compute_coupling_fluid_cmb_cuda");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  int num_blocks_x = mp->nspec2D_top_outer_core;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(5,5,1);

  // launches GPU kernel
  compute_coupling_fluid_CMB_kernel<<<grid,threads>>>(mp->d_displ_crust_mantle,
                                                      mp->d_accel_outer_core,
                                                      mp->d_ibool_crust_mantle,
                                                      mp->d_ibelm_bottom_crust_mantle,
                                                      mp->d_normal_top_outer_core,
                                                      mp->d_jacobian2D_top_outer_core,
                                                      mp->d_wgllwgll_xy,
                                                      mp->d_ibool_outer_core,
                                                      mp->d_ibelm_top_outer_core,
                                                      mp->nspec2D_top_outer_core);

  // adjoint simulations
  if ( mp->simulation_type == 3 ){
    compute_coupling_fluid_CMB_kernel<<<grid,threads>>>(mp->d_b_displ_crust_mantle,
                                                        mp->d_b_accel_outer_core,
                                                        mp->d_ibool_crust_mantle,
                                                        mp->d_ibelm_bottom_crust_mantle,
                                                        mp->d_normal_top_outer_core,
                                                        mp->d_jacobian2D_top_outer_core,
                                                        mp->d_wgllwgll_xy,
                                                        mp->d_ibool_outer_core,
                                                        mp->d_ibelm_top_outer_core,
                                                        mp->nspec2D_top_outer_core);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("compute_coupling_fluid_CMB_kernel");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_coupling_fluid_icb_cuda,
              COMPUTE_COUPLING_FLUID_ICB_CUDA)(long* Mesh_pointer_f) {

  TRACE("compute_coupling_fluid_icb_cuda");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  int num_blocks_x = mp->nspec2D_bottom_outer_core;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(5,5,1);

  // launches GPU kernel
  compute_coupling_fluid_ICB_kernel<<<grid,threads>>>(mp->d_displ_inner_core,
                                                      mp->d_accel_outer_core,
                                                      mp->d_ibool_inner_core,
                                                      mp->d_ibelm_top_inner_core,
                                                      mp->d_normal_bottom_outer_core,
                                                      mp->d_jacobian2D_bottom_outer_core,
                                                      mp->d_wgllwgll_xy,
                                                      mp->d_ibool_outer_core,
                                                      mp->d_ibelm_bottom_outer_core,
                                                      mp->nspec2D_bottom_outer_core);

  // adjoint simulations
  if ( mp->simulation_type == 3 ){
    compute_coupling_fluid_ICB_kernel<<<grid,threads>>>(mp->d_b_displ_inner_core,
                                                        mp->d_b_accel_outer_core,
                                                        mp->d_ibool_inner_core,
                                                        mp->d_ibelm_top_inner_core,
                                                        mp->d_normal_bottom_outer_core,
                                                        mp->d_jacobian2D_bottom_outer_core,
                                                        mp->d_wgllwgll_xy,
                                                        mp->d_ibool_outer_core,
                                                        mp->d_ibelm_bottom_outer_core,
                                                        mp->nspec2D_bottom_outer_core);
  }
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("compute_coupling_fluid_ICB_kernel");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// ELASTIC - ACOUSTIC coupling

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_coupling_CMB_fluid_kernel(realw* displ_crust_mantle,
                                                  realw* accel_crust_mantle,
                                                  realw* accel_outer_core,
                                                  int* ibool_crust_mantle,
                                                  int* ibelm_bottom_crust_mantle,
                                                  realw* normal_top_outer_core,
                                                  realw* jacobian2D_top_outer_core,
                                                  realw* wgllwgll_xy,
                                                  int* ibool_outer_core,
                                                  int* ibelm_top_outer_core,
                                                  double RHO_TOP_OC,
                                                  realw minus_g_cmb,
                                                  int GRAVITY_VAL,
                                                  int NSPEC_BOTTOM_CM) {

  int i = threadIdx.x;
  int j = threadIdx.y;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int k,k_corresp,iglob_cm,iglob_oc,ispec,ispec_selected;
  realw pressure;
  realw nx,ny,nz;
  realw weight;

  // for surfaces elements exactly at the bottom of the crust mantle (outer core top)
  if( iface < NSPEC_BOTTOM_CM ){

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = ibelm_bottom_crust_mantle[iface] - 1;
    ispec_selected = ibelm_top_outer_core[iface] - 1;

    // only for DOFs exactly on the CMB (bottom of these elements)
    k = 0;
    // get velocity potential on the fluid side using pointwise matching
    k_corresp = NGLLX - 1;

    // get normal on the CMB
    nx = normal_top_outer_core[INDEX4(NDIM,NGLLX,NGLLX,0,i,j,iface)]; // (1,i,j,iface)
    ny = normal_top_outer_core[INDEX4(NDIM,NGLLX,NGLLX,1,i,j,iface)]; // (2,i,j,iface)
    nz = normal_top_outer_core[INDEX4(NDIM,NGLLX,NGLLX,2,i,j,iface)]; // (3,i,j,iface)

    // get global point number
    // corresponding points are located at the top of the outer core
    iglob_oc = ibool_outer_core[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k_corresp,ispec_selected)] - 1;
    iglob_cm = ibool_crust_mantle[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] - 1;

    // compute pressure, taking gravity into account
    if( GRAVITY_VAL ){
      pressure = RHO_TOP_OC * ( - accel_outer_core[iglob_oc] + minus_g_cmb
         * (displ_crust_mantle[iglob_cm*3]*nx
            + displ_crust_mantle[iglob_cm*3+1]*ny
            + displ_crust_mantle[iglob_cm*3+2]*nz) );
    }else{
      pressure = - RHO_TOP_OC * accel_outer_core[iglob_oc];
    }

    // formulation with generalized potential: gets associated, weighted jacobian
    weight = jacobian2D_top_outer_core[INDEX3(NGLLX,NGLLX,i,j,iface)]*wgllwgll_xy[INDEX2(NGLLX,i,j)];

    // update fluid acceleration/pressure
    atomicAdd(&accel_crust_mantle[iglob_cm*3],+ weight*nx*pressure);
    atomicAdd(&accel_crust_mantle[iglob_cm*3+1],+ weight*ny*pressure);
    atomicAdd(&accel_crust_mantle[iglob_cm*3+2],+ weight*nz*pressure);

    //  }
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_coupling_ICB_fluid_kernel(realw* displ_inner_core,
                                                  realw* accel_inner_core,
                                                  realw* accel_outer_core,
                                                  int* ibool_inner_core,
                                                  int* ibelm_top_inner_core,
                                                  realw* normal_bottom_outer_core,
                                                  realw* jacobian2D_bottom_outer_core,
                                                  realw* wgllwgll_xy,
                                                  int* ibool_outer_core,
                                                  int* ibelm_bottom_outer_core,
                                                  double RHO_BOTTOM_OC,
                                                  realw minus_g_icb,
                                                  int GRAVITY_VAL,
                                                  int NSPEC_TOP_IC) {

  int i = threadIdx.x;
  int j = threadIdx.y;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int k,k_corresp,iglob_ic,iglob_oc,ispec,ispec_selected;
  realw pressure;
  realw nx,ny,nz;
  realw weight;

  // for surfaces elements exactly at the top of the inner core (outer core bottom)
  if( iface < NSPEC_TOP_IC ){

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = ibelm_top_inner_core[iface] - 1;
    ispec_selected = ibelm_bottom_outer_core[iface] - 1;

    // only for DOFs exactly on the ICB (top of these elements)
    k = NGLLX - 1;
    // get velocity potential on the fluid side using pointwise matching
    k_corresp = 0;

    // get normal on the ICB
    nx = normal_bottom_outer_core[INDEX4(NDIM,NGLLX,NGLLX,0,i,j,iface)]; // (1,i,j,iface)
    ny = normal_bottom_outer_core[INDEX4(NDIM,NGLLX,NGLLX,1,i,j,iface)]; // (2,i,j,iface)
    nz = normal_bottom_outer_core[INDEX4(NDIM,NGLLX,NGLLX,2,i,j,iface)]; // (3,i,j,iface)

    // get global point number
    // corresponding points are located at the bottom of the outer core
    iglob_oc = ibool_outer_core[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k_corresp,ispec_selected)] - 1;
    iglob_ic = ibool_inner_core[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] - 1;

    // compute pressure, taking gravity into account
    if( GRAVITY_VAL ){
      pressure = RHO_BOTTOM_OC * ( - accel_outer_core[iglob_oc] + minus_g_icb
         * (displ_inner_core[iglob_ic*3]*nx
            + displ_inner_core[iglob_ic*3+1]*ny
            + displ_inner_core[iglob_ic*3+2]*nz) );
    }else{
      pressure = - RHO_BOTTOM_OC * accel_outer_core[iglob_oc];
    }

    // formulation with generalized potential: gets associated, weighted jacobian
    weight = jacobian2D_bottom_outer_core[INDEX3(NGLLX,NGLLX,i,j,iface)]*wgllwgll_xy[INDEX2(NGLLX,i,j)];

    // update fluid acceleration/pressure
    atomicAdd(&accel_inner_core[iglob_ic*3],+ weight*nx*pressure);
    atomicAdd(&accel_inner_core[iglob_ic*3+1],+ weight*ny*pressure);
    atomicAdd(&accel_inner_core[iglob_ic*3+2],+ weight*nz*pressure);

    //  }
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_coupling_cmb_fluid_cuda,
              COMPUTE_COUPLING_CMB_FLUID_CUDA)(long* Mesh_pointer_f,
                                               double RHO_TOP_OC,
                                               realw minus_g_cmb,
                                               int GRAVITY_VAL) {

  TRACE("compute_coupling_cmb_fluid_cuda");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  int num_blocks_x = mp->nspec2D_bottom_crust_mantle;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(5,5,1);

  // launches GPU kernel
  compute_coupling_CMB_fluid_kernel<<<grid,threads>>>(mp->d_displ_crust_mantle,
                                                      mp->d_accel_crust_mantle,
                                                      mp->d_accel_outer_core,
                                                      mp->d_ibool_crust_mantle,
                                                      mp->d_ibelm_bottom_crust_mantle,
                                                      mp->d_normal_top_outer_core,
                                                      mp->d_jacobian2D_top_outer_core,
                                                      mp->d_wgllwgll_xy,
                                                      mp->d_ibool_outer_core,
                                                      mp->d_ibelm_top_outer_core,
                                                      RHO_TOP_OC,
                                                      minus_g_cmb,
                                                      GRAVITY_VAL,
                                                      mp->nspec2D_bottom_crust_mantle);

  //  adjoint simulations
  if ( mp->simulation_type == 3 ){
    compute_coupling_CMB_fluid_kernel<<<grid,threads>>>(mp->d_b_displ_crust_mantle,
                                                        mp->d_b_accel_crust_mantle,
                                                        mp->d_b_accel_outer_core,
                                                        mp->d_ibool_crust_mantle,
                                                        mp->d_ibelm_bottom_crust_mantle,
                                                        mp->d_normal_top_outer_core,
                                                        mp->d_jacobian2D_top_outer_core,
                                                        mp->d_wgllwgll_xy,
                                                        mp->d_ibool_outer_core,
                                                        mp->d_ibelm_top_outer_core,
                                                        RHO_TOP_OC,
                                                        minus_g_cmb,
                                                        GRAVITY_VAL,
                                                        mp->nspec2D_bottom_crust_mantle);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("compute_coupling_CMB_fluid_cuda");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_coupling_icb_fluid_cuda,
              COMPUTE_COUPLING_ICB_FLUID_CUDA)(long* Mesh_pointer_f,
                                               double RHO_BOTTOM_OC,
                                               realw minus_g_icb,
                                               int GRAVITY_VAL) {

  TRACE("compute_coupling_icb_fluid_cuda");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  int num_blocks_x = mp->nspec2D_top_inner_core;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(5,5,1);

  // launches GPU kernel
  compute_coupling_ICB_fluid_kernel<<<grid,threads>>>(mp->d_displ_inner_core,
                                                      mp->d_accel_inner_core,
                                                      mp->d_accel_outer_core,
                                                      mp->d_ibool_inner_core,
                                                      mp->d_ibelm_top_inner_core,
                                                      mp->d_normal_bottom_outer_core,
                                                      mp->d_jacobian2D_bottom_outer_core,
                                                      mp->d_wgllwgll_xy,
                                                      mp->d_ibool_outer_core,
                                                      mp->d_ibelm_bottom_outer_core,
                                                      RHO_BOTTOM_OC,
                                                      minus_g_icb,
                                                      GRAVITY_VAL,
                                                      mp->nspec2D_top_inner_core);

  //  adjoint simulations
  if ( mp->simulation_type == 3 ){
    compute_coupling_ICB_fluid_kernel<<<grid,threads>>>(mp->d_b_displ_inner_core,
                                                        mp->d_b_accel_inner_core,
                                                        mp->d_b_accel_outer_core,
                                                        mp->d_ibool_inner_core,
                                                        mp->d_ibelm_top_inner_core,
                                                        mp->d_normal_bottom_outer_core,
                                                        mp->d_jacobian2D_bottom_outer_core,
                                                        mp->d_wgllwgll_xy,
                                                        mp->d_ibool_outer_core,
                                                        mp->d_ibelm_bottom_outer_core,
                                                        RHO_BOTTOM_OC,
                                                        minus_g_icb,
                                                        GRAVITY_VAL,
                                                        mp->nspec2D_top_inner_core);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("compute_coupling_ICB_fluid_cuda");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

/* OCEANS load coupled on free surface */

/* ----------------------------------------------------------------------------------------------- */


__global__ void compute_coupling_ocean_cuda_kernel(realw* accel_crust_mantle,
                                                   realw* rmassx_crust_mantle,
                                                   realw* rmassy_crust_mantle,
                                                   realw* rmassz_crust_mantle,
                                                   realw* rmass_ocean_load,
                                                   realw* normal_top_crust_mantle,
                                                   int* ibool_crust_mantle,
                                                   int* ibelm_top_crust_mantle,
                                                   int* updated_dof_ocean_load,
                                                   int NSPEC2D_TOP_CM) {

  int i = threadIdx.x;
  int j = threadIdx.y;
  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int k,iglob,ispec;
  realw nx,ny,nz;
  realw force_normal_comp;
  realw additional_term_x,additional_term_y,additional_term_z;

  // for surfaces elements exactly at the top of the crust mantle (ocean bottom)
  if( iface < NSPEC2D_TOP_CM ){

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = ibelm_top_crust_mantle[iface] - 1;

    // only for DOFs exactly on the CMB (top of these elements)
    k = NGLLX - 1;

    // get global point number
    iglob = ibool_crust_mantle[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] - 1;

    // only update this global point once

    // daniel: TODO - there might be better ways to implement a mutex like below,
    //            and find a workaround to not use the temporary update array.
    //            As a reminder, atomicExch atomically sets the value of the memory
    //            location stored in address to val and returns the old value.
    //            0 indicates that we still have to do this point

    if( atomicExch(&updated_dof_ocean_load[iglob],1) == 0){

      // get normal on the CMB
      nx = normal_top_crust_mantle[INDEX4(NDIM,NGLLX,NGLLX,0,i,j,iface)]; // (1,i,j,iface)
      ny = normal_top_crust_mantle[INDEX4(NDIM,NGLLX,NGLLX,1,i,j,iface)]; // (2,i,j,iface)
      nz = normal_top_crust_mantle[INDEX4(NDIM,NGLLX,NGLLX,2,i,j,iface)]; // (3,i,j,iface)

      // make updated component of right-hand side
      // we divide by rmass() which is 1 / M
      // we use the total force which includes the Coriolis term above
      force_normal_comp = accel_crust_mantle[iglob*3]*nx / rmassx_crust_mantle[iglob]
                        + accel_crust_mantle[iglob*3+1]*ny / rmassy_crust_mantle[iglob]
                        + accel_crust_mantle[iglob*3+2]*nz / rmassz_crust_mantle[iglob];

      additional_term_x = (rmass_ocean_load[iglob] - rmassx_crust_mantle[iglob]) * force_normal_comp;
      additional_term_y = (rmass_ocean_load[iglob] - rmassy_crust_mantle[iglob]) * force_normal_comp;
      additional_term_z = (rmass_ocean_load[iglob] - rmassz_crust_mantle[iglob]) * force_normal_comp;

      // probably wouldn't need atomicAdd anymore, but just to be sure...
      atomicAdd(&accel_crust_mantle[iglob*3], + additional_term_x * nx);
      atomicAdd(&accel_crust_mantle[iglob*3+1], + additional_term_y * ny);
      atomicAdd(&accel_crust_mantle[iglob*3+2], + additional_term_z * nz);
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_coupling_ocean_cuda,
              COMPUTE_COUPLING_OCEAN_CUDA)(long* Mesh_pointer_f,
                                           int* NCHUNKS_VAL) {

  TRACE("compute_coupling_ocean_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  int num_blocks_x = mp->nspec2D_top_crust_mantle;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(5,5,1);

  // initializes temporary array to zero
  print_CUDA_error_if_any(cudaMemset(mp->d_updated_dof_ocean_load,0,
                                     sizeof(int)*(mp->NGLOB_CRUST_MANTLE_OCEANS)),88501);

  if( *NCHUNKS_VAL != 6 && mp->absorbing_conditions){
    compute_coupling_ocean_cuda_kernel<<<grid,threads>>>(mp->d_accel_crust_mantle,
                                                         mp->d_rmassx_crust_mantle,
                                                         mp->d_rmassy_crust_mantle,
                                                         mp->d_rmassz_crust_mantle,
                                                         mp->d_rmass_ocean_load,
                                                         mp->d_normal_top_crust_mantle,
                                                         mp->d_ibool_crust_mantle,
                                                         mp->d_ibelm_top_crust_mantle,
                                                         mp->d_updated_dof_ocean_load,
                                                         mp->nspec2D_top_crust_mantle);

    // for backward/reconstructed potentials
    if( mp->simulation_type == 3 ) {
      // re-initializes array
      print_CUDA_error_if_any(cudaMemset(mp->d_b_updated_dof_ocean_load,0,
                                         sizeof(int)*(mp->NGLOB_CRUST_MANTLE_OCEANS)),88502);

      compute_coupling_ocean_cuda_kernel<<<grid,threads>>>(mp->d_b_accel_crust_mantle,
                                                           mp->d_rmassx_crust_mantle,
                                                           mp->d_rmassy_crust_mantle,
                                                           mp->d_rmassz_crust_mantle,
                                                           mp->d_rmass_ocean_load,
                                                           mp->d_normal_top_crust_mantle,
                                                           mp->d_ibool_crust_mantle,
                                                           mp->d_ibelm_top_crust_mantle,
                                                           mp->d_b_updated_dof_ocean_load,
                                                           mp->nspec2D_top_crust_mantle);
    }
  }else{
    compute_coupling_ocean_cuda_kernel<<<grid,threads>>>(mp->d_accel_crust_mantle,
                                                         mp->d_rmassz_crust_mantle,
                                                         mp->d_rmassz_crust_mantle,
                                                         mp->d_rmassz_crust_mantle,
                                                         mp->d_rmass_ocean_load,
                                                         mp->d_normal_top_crust_mantle,
                                                         mp->d_ibool_crust_mantle,
                                                         mp->d_ibelm_top_crust_mantle,
                                                         mp->d_updated_dof_ocean_load,
                                                         mp->nspec2D_top_crust_mantle);

    // for backward/reconstructed potentials
    if( mp->simulation_type == 3 ) {
      // re-initializes array
      print_CUDA_error_if_any(cudaMemset(mp->d_b_updated_dof_ocean_load,0,
                                         sizeof(int)*(mp->NGLOB_CRUST_MANTLE_OCEANS)),88502);

      compute_coupling_ocean_cuda_kernel<<<grid,threads>>>(mp->d_b_accel_crust_mantle,
                                                           mp->d_rmassz_crust_mantle,
                                                           mp->d_rmassz_crust_mantle,
                                                           mp->d_rmassz_crust_mantle,
                                                           mp->d_rmass_ocean_load,
                                                           mp->d_normal_top_crust_mantle,
                                                           mp->d_ibool_crust_mantle,
                                                           mp->d_ibelm_top_crust_mantle,
                                                           mp->d_b_updated_dof_ocean_load,
                                                           mp->nspec2D_top_crust_mantle);
    }
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_coupling_ocean_cuda");
#endif
}
