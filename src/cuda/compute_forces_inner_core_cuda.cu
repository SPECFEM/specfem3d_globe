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

#ifdef USE_TEXTURES_FIELDS
texture<realw, cudaTextureType1D, cudaReadModeElementType> d_displ_ic_tex;
texture<realw, cudaTextureType1D, cudaReadModeElementType> d_accel_ic_tex;
#endif

#ifdef USE_TEXTURES_CONSTANTS
texture<realw, cudaTextureType1D, cudaReadModeElementType> d_hprime_xx_ic_tex;
#endif


/* ----------------------------------------------------------------------------------------------- */

// elemental routines

/* ----------------------------------------------------------------------------------------------- */

// updates stress

__device__ void compute_element_ic_att_stress(int tx,int working_element,
                                             realw* R_xx,
                                             realw* R_yy,
                                             realw* R_xy,
                                             realw* R_xz,
                                             realw* R_yz,
                                             realw* sigma_xx,
                                             realw* sigma_yy,
                                             realw* sigma_zz,
                                             realw* sigma_xy,
                                             realw* sigma_xz,
                                             realw* sigma_yz) {

  int offset;
  realw R_xx_val,R_yy_val;

  for(int i_sls = 0; i_sls < N_SLS; i_sls++){
    // index
    // note: index for R_xx,.. here is (i_sls,i,j,k,ispec) and not (i,j,k,ispec,i_sls) as in local version
    //          local version: offset_sls = tx + NGLL3*(working_element + NSPEC*i_sls);
    offset = i_sls + N_SLS*(tx + NGLL3*working_element);

    R_xx_val = R_xx[offset];
    R_yy_val = R_yy[offset];

    *sigma_xx = *sigma_xx - R_xx_val;
    *sigma_yy = *sigma_yy - R_yy_val;
    *sigma_zz = *sigma_zz + R_xx_val + R_yy_val;
    *sigma_xy = *sigma_xy - R_xy[offset];
    *sigma_xz = *sigma_xz - R_xz[offset];
    *sigma_yz = *sigma_yz - R_yz[offset];
  }
}

/* ----------------------------------------------------------------------------------------------- */

// updates R_memory

__device__ void compute_element_ic_att_memory(int tx,int working_element,
                                              realw* d_muv,
                                              realw* factor_common,
                                              realw* alphaval,realw* betaval,realw* gammaval,
                                              realw* R_xx,realw* R_yy,realw* R_xy,realw* R_xz,realw* R_yz,
                                              realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                              realw* epsilondev_xz,realw* epsilondev_yz,
                                              realw epsilondev_xx_loc,realw epsilondev_yy_loc,realw epsilondev_xy_loc,
                                              realw epsilondev_xz_loc,realw epsilondev_yz_loc,
                                              int USE_3D_ATTENUATION_ARRAYS
                                              ){

  int offset;
  realw mul;
  realw alphaval_loc,betaval_loc,gammaval_loc;
  realw factor_loc,Sn,Snp1;

  mul = d_muv[tx + NGLL3_PADDED * working_element];

  // use Runge-Kutta scheme to march in time
  for(int i_sls = 0; i_sls < N_SLS; i_sls++){

    // indices
    // note: index for R_xx,... here is (i_sls,i,j,k,ispec) and not (i,j,k,ispec,i_sls) as in local version
    //          local version: offset_sls = tx + NGLL3*(working_element + NSPEC*i_sls);
    // index for (i_sls,i,j,k,ispec)
    offset = i_sls + N_SLS*(tx + NGLL3*working_element);

    if( USE_3D_ATTENUATION_ARRAYS ){
      factor_loc = mul * factor_common[offset]; //mustore(i,j,k,ispec) * factor_common(i_sls,i,j,k,ispec)
    }else{
      factor_loc = mul * factor_common[i_sls + N_SLS*working_element]; //mustore(i,j,k,ispec) * factor_common(i_sls,1,1,1,ispec)
    }
    alphaval_loc = alphaval[i_sls]; // (i_sls)
    betaval_loc = betaval[i_sls];
    gammaval_loc = gammaval[i_sls];

    // term in xx
    Sn   = factor_loc * epsilondev_xx[tx + NGLL3 * working_element]; //(i,j,k,ispec)
    Snp1   = factor_loc * epsilondev_xx_loc; //(i,j,k)
    R_xx[offset] = alphaval_loc * R_xx[offset] + betaval_loc * Sn + gammaval_loc * Snp1;

    // term in yy
    Sn   = factor_loc * epsilondev_yy[tx + NGLL3 * working_element];
    Snp1   = factor_loc * epsilondev_yy_loc;
    R_yy[offset] = alphaval_loc * R_yy[offset] + betaval_loc * Sn + gammaval_loc * Snp1;
    // term in zz not computed since zero trace

    // term in xy
    Sn   = factor_loc * epsilondev_xy[tx + NGLL3 * working_element];
    Snp1   = factor_loc * epsilondev_xy_loc;
    R_xy[offset] = alphaval_loc * R_xy[offset] + betaval_loc * Sn + gammaval_loc * Snp1;

    // term in xz
    Sn   = factor_loc * epsilondev_xz[tx + NGLL3 * working_element];
    Snp1   = factor_loc * epsilondev_xz_loc;
    R_xz[offset] = alphaval_loc * R_xz[offset] + betaval_loc * Sn + gammaval_loc * Snp1;

    // term in yz
    Sn   = factor_loc * epsilondev_yz[tx + NGLL3 * working_element];
    Snp1   = factor_loc * epsilondev_yz_loc;
    R_yz[offset] = alphaval_loc * R_yz[offset] + betaval_loc * Sn + gammaval_loc * Snp1;
  }
}

/* ----------------------------------------------------------------------------------------------- */

// pre-computes gravity term

__device__ void compute_element_ic_gravity(int tx,int working_element,
                                           int* d_ibool,
                                           realw* d_xstore,realw* d_ystore,realw* d_zstore,
                                           realw* d_minus_gravity_table,
                                           realw* d_minus_deriv_gravity_table,
                                           realw* d_density_table,
                                           realw* wgll_cube,
                                           realw jacobianl,
                                           realw* s_dummyx_loc,
                                           realw* s_dummyy_loc,
                                           realw* s_dummyz_loc,
                                           realw* sigma_xx,
                                           realw* sigma_yy,
                                           realw* sigma_zz,
                                           realw* sigma_xy,
                                           realw* sigma_yx,
                                           realw* sigma_xz,
                                           realw* sigma_zx,
                                           realw* sigma_yz,
                                           realw* sigma_zy,
                                           realw* rho_s_H1,
                                           realw* rho_s_H2,
                                           realw* rho_s_H3){

  realw radius,theta,phi;
  realw cos_theta,sin_theta,cos_phi,sin_phi;
  realw minus_g,minus_dg;
  realw rho;
  realw gxl,gyl,gzl;
  realw minus_g_over_radius,minus_dg_plus_g_over_radius;
  realw cos_theta_sq,sin_theta_sq,cos_phi_sq,sin_phi_sq;
  realw Hxxl,Hyyl,Hzzl,Hxyl,Hxzl,Hyzl;
  realw sx_l,sy_l,sz_l;
  realw factor;

  // R_EARTH_KM is the radius of the bottom of the oceans
  //const realw R_EARTH = 6371000.0f; // in m
  //const realw R_EARTH_KM = 6371.0f; // in km
  // uncomment line below for PREM with oceans
  //const realw R_EARTH = 6368000.0f;
  //const realw R_EARTH_KM = 6368.0f;

  // compute non-symmetric terms for gravity

  // use mesh coordinates to get theta and phi
  // x y z contain r theta phi
  int iglob = d_ibool[working_element*NGLL3 + tx]-1;

  radius = d_xstore[iglob];
  // make sure radius is never zero even for points at center of cube
  // because we later divide by radius
  if(radius < 100.f / (R_EARTH_KM*1000.0f)){ radius = 100.f / (R_EARTH_KM*1000.0f); }

  theta = d_ystore[iglob];
  phi = d_zstore[iglob];

  if( sizeof( theta ) == sizeof( float ) ){
    // float operations
    // sincos function return sinus and cosine for given value
    sincosf(theta, &sin_theta, &cos_theta);
    sincosf(phi, &sin_phi, &cos_phi);
  }else{
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    cos_phi = cos(phi);
    sin_phi = sin(phi);
  }

  // for efficiency replace with lookup table every 100 m in radial direction
  // note: radius in crust mantle should never be zero,
  //          and arrays in C start from 0, thus we need to subtract -1
  int int_radius = rint(radius * R_EARTH_KM * 10.0f ) - 1;
  //make sure we never use below zero for point exactly at the center of the Earth
  if( int_radius < 0 ){int_radius = 0;}

  // get g, rho and dg/dr=dg
  // spherical components of the gravitational acceleration
  // for efficiency replace with lookup table every 100 m in radial direction
  minus_g = d_minus_gravity_table[int_radius];
  minus_dg = d_minus_deriv_gravity_table[int_radius];
  rho = d_density_table[int_radius];

  // Cartesian components of the gravitational acceleration
  gxl = minus_g*sin_theta*cos_phi;
  gyl = minus_g*sin_theta*sin_phi;
  gzl = minus_g*cos_theta;

  // Cartesian components of gradient of gravitational acceleration
  // obtained from spherical components

  minus_g_over_radius = minus_g / radius;
  minus_dg_plus_g_over_radius = minus_dg - minus_g_over_radius;

  cos_theta_sq = cos_theta*cos_theta;
  sin_theta_sq = sin_theta*sin_theta;
  cos_phi_sq = cos_phi*cos_phi;
  sin_phi_sq = sin_phi*sin_phi;

  Hxxl = minus_g_over_radius*(cos_phi_sq*cos_theta_sq + sin_phi_sq) + cos_phi_sq*minus_dg*sin_theta_sq;
  Hyyl = minus_g_over_radius*(cos_phi_sq + cos_theta_sq*sin_phi_sq) + minus_dg*sin_phi_sq*sin_theta_sq;
  Hzzl = cos_theta_sq*minus_dg + minus_g_over_radius*sin_theta_sq;
  Hxyl = cos_phi*minus_dg_plus_g_over_radius*sin_phi*sin_theta_sq;
  Hxzl = cos_phi*cos_theta*minus_dg_plus_g_over_radius*sin_theta;
  Hyzl = cos_theta*minus_dg_plus_g_over_radius*sin_phi*sin_theta;

  // get displacement and multiply by density to compute G tensor
  sx_l = rho * s_dummyx_loc[tx];
  sy_l = rho * s_dummyy_loc[tx];
  sz_l = rho * s_dummyz_loc[tx];

  // compute G tensor from s . g and add to sigma (not symmetric)
  *sigma_xx = *sigma_xx + sy_l*gyl + sz_l*gzl;
  *sigma_yy = *sigma_yy + sx_l*gxl + sz_l*gzl;
  *sigma_zz = *sigma_zz + sx_l*gxl + sy_l*gyl;

  *sigma_xy = *sigma_xy - sx_l * gyl;
  *sigma_yx = *sigma_yx - sy_l * gxl;

  *sigma_xz = *sigma_xz - sx_l * gzl;
  *sigma_zx = *sigma_zx - sz_l * gxl;

  *sigma_yz = *sigma_yz - sy_l * gzl;
  *sigma_zy = *sigma_zy - sz_l * gyl;

  // precompute vector
  factor = jacobianl * wgll_cube[tx];
  *rho_s_H1 = factor * (sx_l * Hxxl + sy_l * Hxyl + sz_l * Hxzl);
  *rho_s_H2 = factor * (sx_l * Hxyl + sy_l * Hyyl + sz_l * Hyzl);
  *rho_s_H3 = factor * (sx_l * Hxzl + sy_l * Hyzl + sz_l * Hzzl);
}


/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2
//
// for inner_core

/* ----------------------------------------------------------------------------------------------- */

__global__ void Kernel_2_inner_core_impl(int nb_blocks_to_compute,
                                         int NGLOB,
                                         int* d_ibool,
                                         int* d_idoubling,
                                         int* d_phase_ispec_inner,
                                         int num_phase_ispec,
                                         int d_iphase,
                                         realw deltat,
                                         int use_mesh_coloring_gpu,
                                         realw* d_displ,
                                         realw* d_veloc,
                                         realw* d_accel,
                                         realw* d_xix, realw* d_xiy, realw* d_xiz,
                                         realw* d_etax, realw* d_etay, realw* d_etaz,
                                         realw* d_gammax, realw* d_gammay, realw* d_gammaz,
                                         realw* d_hprime_xx,
                                         realw* d_hprimewgll_xx,
                                         realw* d_wgllwgll_xy,realw* d_wgllwgll_xz,realw* d_wgllwgll_yz,
                                         realw* d_kappav,
                                         realw* d_muv,
                                         int COMPUTE_AND_STORE_STRAIN,
                                         realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                         realw* epsilondev_xz,realw* epsilondev_yz,
                                         realw* epsilon_trace_over_3,
                                         int ATTENUATION,
                                         int PARTIAL_PHYS_DISPERSION_ONLY,
                                         int USE_3D_ATTENUATION_ARRAYS,
                                         realw* one_minus_sum_beta,realw* factor_common,
                                         realw* R_xx, realw* R_yy, realw* R_xy, realw* R_xz, realw* R_yz,
                                         realw* alphaval,realw* betaval,realw* gammaval,
                                         int ANISOTROPY,
                                         realw* d_c11store,realw* d_c12store,realw* d_c13store,
                                         realw* d_c33store,realw* d_c44store,
                                         int GRAVITY,
                                         realw* d_xstore,realw* d_ystore,realw* d_zstore,
                                         realw* d_minus_gravity_table,
                                         realw* d_minus_deriv_gravity_table,
                                         realw* d_density_table,
                                         realw* wgll_cube,
                                         int NSPEC_INNER_CORE_STRAIN_ONLY,
                                         int NSPEC_INNER_CORE){

  // block id
  int bx = blockIdx.y*gridDim.x+blockIdx.x;
  // thread id
  int tx = threadIdx.x;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int active,offset;
  int iglob;
  int working_element;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;
  realw templ;

  realw fac1,fac2,fac3;
  realw lambdal,mul,lambdalplus2mul,kappal;
  realw mul_iso,mul_aniso;

  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc;
  realw c11,c12,c13,c33,c44;
  realw sum_terms1,sum_terms2,sum_terms3;

  // gravity variables
  realw sigma_yx,sigma_zx,sigma_zy;
  realw rho_s_H1,rho_s_H2,rho_s_H3;

#ifndef MANUALLY_UNROLLED_LOOPS
  int l;
#endif

  __shared__ realw s_dummyx_loc[NGLL3];
  __shared__ realw s_dummyy_loc[NGLL3];
  __shared__ realw s_dummyz_loc[NGLL3];

  __shared__ realw s_tempx1[NGLL3];
  __shared__ realw s_tempx2[NGLL3];
  __shared__ realw s_tempx3[NGLL3];

  __shared__ realw s_tempy1[NGLL3];
  __shared__ realw s_tempy2[NGLL3];
  __shared__ realw s_tempy3[NGLL3];

  __shared__ realw s_tempz1[NGLL3];
  __shared__ realw s_tempz2[NGLL3];
  __shared__ realw s_tempz3[NGLL3];

  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

// use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
// because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses
  active = (tx < NGLL3 && bx < nb_blocks_to_compute) ? 1:0;

// copy from global memory to shared memory
// each thread writes one of the NGLL^3 = 125 data points
  if (active) {

#ifdef USE_MESH_COLORING_GPU
    working_element = bx;
#else
    //mesh coloring
    if( use_mesh_coloring_gpu ){
      working_element = bx;
    }else{
      // iphase-1 and working_element-1 for Fortran->C array conventions
      working_element = d_phase_ispec_inner[bx + num_phase_ispec*(d_iphase-1)]-1;
    }
#endif

    // exclude fictitious elements in central cube
    if( d_idoubling[working_element] == IFLAG_IN_FICTITIOUS_CUBE ){
      active = 0;
    }else{
      // iglob = d_ibool[working_element*NGLL3_PADDED + tx]-1;
      iglob = d_ibool[working_element*NGLL3 + tx]-1;

#ifdef USE_TEXTURES_FIELDS
      s_dummyx_loc[tx] = tex1Dfetch(d_displ_ic_tex, iglob*3);
      s_dummyy_loc[tx] = tex1Dfetch(d_displ_ic_tex, iglob*3 + 1);
      s_dummyz_loc[tx] = tex1Dfetch(d_displ_ic_tex, iglob*3 + 2);
#else
      // changing iglob indexing to match fortran row changes fast style
      s_dummyx_loc[tx] = d_displ[iglob*3];
      s_dummyy_loc[tx] = d_displ[iglob*3 + 1];
      s_dummyz_loc[tx] = d_displ[iglob*3 + 2];
#endif
    }
  }

  if (tx < NGLL2) {
    // hprime
#ifdef USE_TEXTURES_CONSTANTS
    sh_hprime_xx[tx] = tex1Dfetch(d_hprime_xx_ic_tex,tx);
#else
    sh_hprime_xx[tx] = d_hprime_xx[tx];
#endif
    // weighted hprime
    sh_hprimewgll_xx[tx] = d_hprimewgll_xx[tx];
  }

// synchronize all the threads (one thread for each of the NGLL grid points of the
// current spectral element) because we need the whole element to be ready in order
// to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  if (active) {

#ifndef MANUALLY_UNROLLED_LOOPS
    tempx1l = 0.f;
    tempx2l = 0.f;
    tempx3l = 0.f;

    tempy1l = 0.f;
    tempy2l = 0.f;
    tempy3l = 0.f;

    tempz1l = 0.f;
    tempz2l = 0.f;
    tempz3l = 0.f;

    for (l=0;l<NGLLX;l++) {
        fac1 = sh_hprime_xx[l*NGLLX+I];
        offset = K*NGLL2+J*NGLLX+l;
        tempx1l += s_dummyx_loc[offset]*fac1;
        tempy1l += s_dummyy_loc[offset]*fac1;
        tempz1l += s_dummyz_loc[offset]*fac1;

        fac2 = sh_hprime_xx[l*NGLLX+J];
        offset = K*NGLL2+l*NGLLX+I;
        tempx2l += s_dummyx_loc[offset]*fac2;
        tempy2l += s_dummyy_loc[offset]*fac2;
        tempz2l += s_dummyz_loc[offset]*fac2;

        fac3 = sh_hprime_xx[l*NGLLX+K];
        offset = l*NGLL2+J*NGLLX+I;
        tempx3l += s_dummyx_loc[offset]*fac3;
        tempy3l += s_dummyy_loc[offset]*fac3;
        tempz3l += s_dummyz_loc[offset]*fac3;
    }

#else

    tempx1l = s_dummyx_loc[K*NGLL2+J*NGLLX]*sh_hprime_xx[I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+1]*sh_hprime_xx[NGLLX+I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+2]*sh_hprime_xx[2*NGLLX+I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+3]*sh_hprime_xx[3*NGLLX+I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+4]*sh_hprime_xx[4*NGLLX+I];

    tempy1l = s_dummyy_loc[K*NGLL2+J*NGLLX]*sh_hprime_xx[I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+1]*sh_hprime_xx[NGLLX+I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+2]*sh_hprime_xx[2*NGLLX+I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+3]*sh_hprime_xx[3*NGLLX+I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+4]*sh_hprime_xx[4*NGLLX+I];

    tempz1l = s_dummyz_loc[K*NGLL2+J*NGLLX]*sh_hprime_xx[I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+1]*sh_hprime_xx[NGLLX+I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+2]*sh_hprime_xx[2*NGLLX+I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+3]*sh_hprime_xx[3*NGLLX+I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+4]*sh_hprime_xx[4*NGLLX+I];

    tempx2l = s_dummyx_loc[K*NGLL2+I]*sh_hprime_xx[J]
            + s_dummyx_loc[K*NGLL2+NGLLX+I]*sh_hprime_xx[NGLLX+J]
            + s_dummyx_loc[K*NGLL2+2*NGLLX+I]*sh_hprime_xx[2*NGLLX+J]
            + s_dummyx_loc[K*NGLL2+3*NGLLX+I]*sh_hprime_xx[3*NGLLX+J]
            + s_dummyx_loc[K*NGLL2+4*NGLLX+I]*sh_hprime_xx[4*NGLLX+J];

    tempy2l = s_dummyy_loc[K*NGLL2+I]*sh_hprime_xx[J]
            + s_dummyy_loc[K*NGLL2+NGLLX+I]*sh_hprime_xx[NGLLX+J]
            + s_dummyy_loc[K*NGLL2+2*NGLLX+I]*sh_hprime_xx[2*NGLLX+J]
            + s_dummyy_loc[K*NGLL2+3*NGLLX+I]*sh_hprime_xx[3*NGLLX+J]
            + s_dummyy_loc[K*NGLL2+4*NGLLX+I]*sh_hprime_xx[4*NGLLX+J];

    tempz2l = s_dummyz_loc[K*NGLL2+I]*sh_hprime_xx[J]
            + s_dummyz_loc[K*NGLL2+NGLLX+I]*sh_hprime_xx[NGLLX+J]
            + s_dummyz_loc[K*NGLL2+2*NGLLX+I]*sh_hprime_xx[2*NGLLX+J]
            + s_dummyz_loc[K*NGLL2+3*NGLLX+I]*sh_hprime_xx[3*NGLLX+J]
            + s_dummyz_loc[K*NGLL2+4*NGLLX+I]*sh_hprime_xx[4*NGLLX+J];

    tempx3l = s_dummyx_loc[J*NGLLX+I]*sh_hprime_xx[K]
            + s_dummyx_loc[NGLL2+J*NGLLX+I]*sh_hprime_xx[NGLLX+K]
            + s_dummyx_loc[2*NGLL2+J*NGLLX+I]*sh_hprime_xx[2*NGLLX+K]
            + s_dummyx_loc[3*NGLL2+J*NGLLX+I]*sh_hprime_xx[3*NGLLX+K]
            + s_dummyx_loc[4*NGLL2+J*NGLLX+I]*sh_hprime_xx[4*NGLLX+K];

    tempy3l = s_dummyy_loc[J*NGLLX+I]*sh_hprime_xx[K]
            + s_dummyy_loc[NGLL2+J*NGLLX+I]*sh_hprime_xx[NGLLX+K]
            + s_dummyy_loc[2*NGLL2+J*NGLLX+I]*sh_hprime_xx[2*NGLLX+K]
            + s_dummyy_loc[3*NGLL2+J*NGLLX+I]*sh_hprime_xx[3*NGLLX+K]
            + s_dummyy_loc[4*NGLL2+J*NGLLX+I]*sh_hprime_xx[4*NGLLX+K];

    tempz3l = s_dummyz_loc[J*NGLLX+I]*sh_hprime_xx[K]
            + s_dummyz_loc[NGLL2+J*NGLLX+I]*sh_hprime_xx[NGLLX+K]
            + s_dummyz_loc[2*NGLL2+J*NGLLX+I]*sh_hprime_xx[2*NGLLX+K]
            + s_dummyz_loc[3*NGLL2+J*NGLLX+I]*sh_hprime_xx[3*NGLLX+K]
            + s_dummyz_loc[4*NGLL2+J*NGLLX+I]*sh_hprime_xx[4*NGLLX+K];

#endif

// compute derivatives of ux, uy and uz with respect to x, y and z
    offset = working_element*NGLL3_PADDED + tx;

    xixl = d_xix[offset];
    xiyl = d_xiy[offset];
    xizl = d_xiz[offset];
    etaxl = d_etax[offset];
    etayl = d_etay[offset];
    etazl = d_etaz[offset];
    gammaxl = d_gammax[offset];
    gammayl = d_gammay[offset];
    gammazl = d_gammaz[offset];

    duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l;
    duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l;
    duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l;

    duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l;
    duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l;
    duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l;

    duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l;
    duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l;
    duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l;

    // precompute some sums to save CPU time
    duxdxl_plus_duydyl = duxdxl + duydyl;
    duxdxl_plus_duzdzl = duxdxl + duzdzl;
    duydyl_plus_duzdzl = duydyl + duzdzl;
    duxdyl_plus_duydxl = duxdyl + duydxl;
    duzdxl_plus_duxdzl = duzdxl + duxdzl;
    duzdyl_plus_duydzl = duzdyl + duydzl;

    // computes deviatoric strain attenuation and/or for kernel calculations
    if(COMPUTE_AND_STORE_STRAIN) {
      templ = 0.33333333333333333333f * (duxdxl + duydyl + duzdzl); // 1./3. = 0.33333

      // local storage: stresses at this current time step
      epsilondev_xx_loc = duxdxl - templ;
      epsilondev_yy_loc = duydyl - templ;
      epsilondev_xy_loc = 0.5f * duxdyl_plus_duydxl;
      epsilondev_xz_loc = 0.5f * duzdxl_plus_duxdzl;
      epsilondev_yz_loc = 0.5f * duzdyl_plus_duydzl;

      if(NSPEC_INNER_CORE_STRAIN_ONLY == 1) {
        epsilon_trace_over_3[tx] = templ;
      }else{
        epsilon_trace_over_3[tx + working_element*NGLL3] = templ;
      }
    }

    // compute elements with an elastic isotropic rheology
    kappal = d_kappav[offset];
    mul = d_muv[offset];

    // attenuation
    if(ATTENUATION){
      // use unrelaxed parameters if attenuation
      if( USE_3D_ATTENUATION_ARRAYS ){
        mul_iso  = mul * one_minus_sum_beta[tx+working_element*NGLL3]; // (i,j,k,ispec)
        mul_aniso = mul *( one_minus_sum_beta[tx+working_element*NGLL3] - 1.0f );
      }else{
        mul_iso  = mul * one_minus_sum_beta[working_element]; // (1,1,1,ispec)
        mul_aniso = mul *( one_minus_sum_beta[working_element] - 1.0f );
      }
    }else{
      mul_iso = mul;
    }

    // full anisotropic case, stress calculations
    if(ANISOTROPY){

      // elastic tensor for hexagonal symmetry in reduced notation:
      //
      //      c11 c12 c13  0   0        0
      //      c12 c11 c13  0   0        0
      //      c13 c13 c33  0   0        0
      //       0   0   0  c44  0        0
      //       0   0   0   0  c44       0
      //       0   0   0   0   0  (c11-c12)/2
      //
      //       in terms of the A, C, L, N and F of Love (1927):
      //
      //       c11 = A
      //       c12 = A-2N
      //       c13 = F
      //       c33 = C
      //       c44 = L

      c11 = d_c11store[offset];
      c12 = d_c12store[offset];
      c13 = d_c13store[offset];
      c33 = d_c33store[offset];
      c44 = d_c44store[offset];

      // use unrelaxed parameters if attenuation
      if( ATTENUATION){
        c11 = c11 + 1.33333333333333333333f * mul_aniso; // FOUR_THIRDS = 1.33333
        c12 = c12 - 0.66666666666666666666f * mul_aniso; // TWO_THIRDS = 0.66666666666666666666f
        c13 = c13 - 0.66666666666666666666f * mul_aniso;
        c33 = c33 + 1.33333333333333333333f * mul_aniso;
        c44 = c44 + mul_aniso;
      }

      sigma_xx = c11*duxdxl + c12*duydyl + c13*duzdzl;
      sigma_yy = c12*duxdxl + c11*duydyl + c13*duzdzl;
      sigma_zz = c13*duxdxl + c13*duydyl + c33*duzdzl;
      sigma_xy = 0.5f*(c11-c12)*duxdyl_plus_duydxl;
      sigma_xz = c44*duzdxl_plus_duxdzl;
      sigma_yz = c44*duzdyl_plus_duydzl;

    }else{

      // isotropic case

      lambdalplus2mul = kappal + 1.33333333333333333333f * mul_iso;  // 4./3. = 1.3333333
      lambdal = lambdalplus2mul - 2.0f * mul_iso;

      // compute the six components of the stress tensor sigma
      sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
      sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
      sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

      sigma_xy = mul*duxdyl_plus_duydxl;
      sigma_xz = mul*duzdxl_plus_duxdzl;
      sigma_yz = mul*duzdyl_plus_duydzl;
    }

    if(ATTENUATION && ( ! PARTIAL_PHYS_DISPERSION_ONLY ) ){
      // subtracts memory variables if attenuation
      compute_element_ic_att_stress(tx,working_element,
                                    R_xx,R_yy,R_xy,R_xz,R_yz,
                                    &sigma_xx,&sigma_yy,&sigma_zz,&sigma_xy,&sigma_xz,&sigma_yz);
    }

    // define symmetric components (needed for non-symmetric dot product and sigma for gravity)
    sigma_yx = sigma_xy;
    sigma_zx = sigma_xz;
    sigma_zy = sigma_yz;

    // jacobian
    jacobianl = 1.0f / (xixl*(etayl*gammazl-etazl*gammayl)
                        -xiyl*(etaxl*gammazl-etazl*gammaxl)
                        +xizl*(etaxl*gammayl-etayl*gammaxl));

    if( GRAVITY ){
      //  computes non-symmetric terms for gravity
      compute_element_ic_gravity(tx,working_element,
                                 d_ibool,d_xstore,d_ystore,d_zstore,
                                 d_minus_gravity_table,d_minus_deriv_gravity_table,d_density_table,
                                 wgll_cube,jacobianl,
                                 s_dummyx_loc,s_dummyy_loc,s_dummyz_loc,
                                 &sigma_xx,&sigma_yy,&sigma_zz,&sigma_xy,&sigma_yx,
                                 &sigma_xz,&sigma_zx,&sigma_yz,&sigma_zy,
                                 &rho_s_H1,&rho_s_H2,&rho_s_H3);
    }

    // form dot product with test vector, non-symmetric form
    s_tempx1[tx] = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl);
    s_tempy1[tx] = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl);
    s_tempz1[tx] = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl);

    s_tempx2[tx] = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl);
    s_tempy2[tx] = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl);
    s_tempz2[tx] = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl);

    s_tempx3[tx] = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl);
    s_tempy3[tx] = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl);
    s_tempz3[tx] = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl);

  }

// synchronize all the threads (one thread for each of the NGLL grid points of the
// current spectral element) because we need the whole element to be ready in order
// to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  if (active) {

#ifndef MANUALLY_UNROLLED_LOOPS

    tempx1l = 0.f;
    tempy1l = 0.f;
    tempz1l = 0.f;

    tempx2l = 0.f;
    tempy2l = 0.f;
    tempz2l = 0.f;

    tempx3l = 0.f;
    tempy3l = 0.f;
    tempz3l = 0.f;

    for (l=0;l<NGLLX;l++) {

      fac1 = sh_hprimewgll_xx[I*NGLLX+l];
      offset = K*NGLL2+J*NGLLX+l;
      tempx1l += s_tempx1[offset]*fac1;
      tempy1l += s_tempy1[offset]*fac1;
      tempz1l += s_tempz1[offset]*fac1;

      // assumes hprimewgll_xx == hprimewgll_yy == hprimewgll_zz
      fac2 = sh_hprimewgll_xx[J*NGLLX+l];
      offset = K*NGLL2+l*NGLLX+I;
      tempx2l += s_tempx2[offset]*fac2;
      tempy2l += s_tempy2[offset]*fac2;
      tempz2l += s_tempz2[offset]*fac2;

      fac3 = sh_hprimewgll_xx[K*NGLLX+l];
      offset = l*NGLL2+J*NGLLX+I;
      tempx3l += s_tempx3[offset]*fac3;
      tempy3l += s_tempy3[offset]*fac3;
      tempz3l += s_tempz3[offset]*fac3;

    }
#else

    tempx1l = s_tempx1[K*NGLL2+J*NGLLX]*sh_hprimewgll_xx[I*NGLLX]
            + s_tempx1[K*NGLL2+J*NGLLX+1]*sh_hprimewgll_xx[I*NGLLX+1]
            + s_tempx1[K*NGLL2+J*NGLLX+2]*sh_hprimewgll_xx[I*NGLLX+2]
            + s_tempx1[K*NGLL2+J*NGLLX+3]*sh_hprimewgll_xx[I*NGLLX+3]
            + s_tempx1[K*NGLL2+J*NGLLX+4]*sh_hprimewgll_xx[I*NGLLX+4];

    tempy1l = s_tempy1[K*NGLL2+J*NGLLX]*sh_hprimewgll_xx[I*NGLLX]
            + s_tempy1[K*NGLL2+J*NGLLX+1]*sh_hprimewgll_xx[I*NGLLX+1]
            + s_tempy1[K*NGLL2+J*NGLLX+2]*sh_hprimewgll_xx[I*NGLLX+2]
            + s_tempy1[K*NGLL2+J*NGLLX+3]*sh_hprimewgll_xx[I*NGLLX+3]
            + s_tempy1[K*NGLL2+J*NGLLX+4]*sh_hprimewgll_xx[I*NGLLX+4];

    tempz1l = s_tempz1[K*NGLL2+J*NGLLX]*sh_hprimewgll_xx[I*NGLLX]
            + s_tempz1[K*NGLL2+J*NGLLX+1]*sh_hprimewgll_xx[I*NGLLX+1]
            + s_tempz1[K*NGLL2+J*NGLLX+2]*sh_hprimewgll_xx[I*NGLLX+2]
            + s_tempz1[K*NGLL2+J*NGLLX+3]*sh_hprimewgll_xx[I*NGLLX+3]
            + s_tempz1[K*NGLL2+J*NGLLX+4]*sh_hprimewgll_xx[I*NGLLX+4];

    tempx2l = s_tempx2[K*NGLL2+I]*sh_hprimewgll_xx[J*NGLLX]
            + s_tempx2[K*NGLL2+NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+1]
            + s_tempx2[K*NGLL2+2*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+2]
            + s_tempx2[K*NGLL2+3*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+3]
            + s_tempx2[K*NGLL2+4*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+4];

    tempy2l = s_tempy2[K*NGLL2+I]*sh_hprimewgll_xx[J*NGLLX]
            + s_tempy2[K*NGLL2+NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+1]
            + s_tempy2[K*NGLL2+2*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+2]
            + s_tempy2[K*NGLL2+3*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+3]
            + s_tempy2[K*NGLL2+4*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+4];

    tempz2l = s_tempz2[K*NGLL2+I]*sh_hprimewgll_xx[J*NGLLX]
            + s_tempz2[K*NGLL2+NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+1]
            + s_tempz2[K*NGLL2+2*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+2]
            + s_tempz2[K*NGLL2+3*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+3]
            + s_tempz2[K*NGLL2+4*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+4];

    tempx3l = s_tempx3[J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX]
            + s_tempx3[NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+1]
            + s_tempx3[2*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+2]
            + s_tempx3[3*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+3]
            + s_tempx3[4*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+4];

    tempy3l = s_tempy3[J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX]
            + s_tempy3[NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+1]
            + s_tempy3[2*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+2]
            + s_tempy3[3*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+3]
            + s_tempy3[4*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+4];

    tempz3l = s_tempz3[J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX]
            + s_tempz3[NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+1]
            + s_tempz3[2*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+2]
            + s_tempz3[3*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+3]
            + s_tempz3[4*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+4];

#endif

    fac1 = d_wgllwgll_yz[K*NGLLX+J];
    fac2 = d_wgllwgll_xz[K*NGLLX+I];
    fac3 = d_wgllwgll_xy[J*NGLLX+I];

    sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
    sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
    sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

    // adds gravity term
    if( GRAVITY ){
      sum_terms1 += rho_s_H1;
      sum_terms2 += rho_s_H2;
      sum_terms3 += rho_s_H3;
    }


#ifdef USE_MESH_COLORING_GPU
    // no atomic operation needed, colors don't share global points between elements

#ifdef USE_TEXTURES_FIELDS
    d_accel[iglob*3]     = tex1Dfetch(d_accel_ic_tex, iglob*3) + sum_terms1;
    d_accel[iglob*3 + 1] = tex1Dfetch(d_accel_ic_tex, iglob*3 + 1) + sum_terms2;
    d_accel[iglob*3 + 2] = tex1Dfetch(d_accel_ic_tex, iglob*3 + 2) + sum_terms3;
#else
    d_accel[iglob*3]     += sum_terms1;
    d_accel[iglob*3 + 1] += sum_terms2;
    d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

#else // MESH_COLORING

    //mesh coloring
    if( use_mesh_coloring_gpu ){

      if( NSPEC_INNER_CORE > COLORING_MIN_NSPEC_INNER_CORE ){
        // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
        d_accel[iglob*3]     = tex1Dfetch(d_accel_ic_tex, iglob*3) + sum_terms1;
        d_accel[iglob*3 + 1] = tex1Dfetch(d_accel_ic_tex, iglob*3 + 1) + sum_terms2;
        d_accel[iglob*3 + 2] = tex1Dfetch(d_accel_ic_tex, iglob*3 + 2) + sum_terms3;
#else
        d_accel[iglob*3]     += sum_terms1;
        d_accel[iglob*3 + 1] += sum_terms2;
        d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS
      }else{
        // poor element count, only use 1 color per inner/outer run
        // forces atomic operations
        atomicAdd(&d_accel[iglob*3], sum_terms1);
        atomicAdd(&d_accel[iglob*3+1], sum_terms2);
        atomicAdd(&d_accel[iglob*3+2], sum_terms3);
      }

    }else{

      // for testing purposes only: w/out atomic updates
      //d_accel[iglob*3] -= (0.00000001f*tempx1l + 0.00000001f*tempx2l + 0.00000001f*tempx3l);
      //d_accel[iglob*3 + 1] -= (0.00000001f*tempy1l + 0.00000001f*tempy2l + 0.00000001f*tempy3l);
      //d_accel[iglob*3 + 2] -= (0.00000001f*tempz1l + 0.00000001f*tempz2l + 0.00000001f*tempz3l);

      atomicAdd(&d_accel[iglob*3], sum_terms1);
      atomicAdd(&d_accel[iglob*3+1], sum_terms2);
      atomicAdd(&d_accel[iglob*3+2], sum_terms3);

    }
#endif // MESH_COLORING

    // update memory variables based upon the Runge-Kutta scheme
    if( ATTENUATION && ! PARTIAL_PHYS_DISPERSION_ONLY ){
      compute_element_ic_att_memory(tx,working_element,
                                d_muv,
                                factor_common,alphaval,betaval,gammaval,
                                R_xx,R_yy,R_xy,R_xz,R_yz,
                                epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz,
                                epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc,
                                USE_3D_ATTENUATION_ARRAYS);
    }

    // save deviatoric strain for Runge-Kutta scheme
    if( COMPUTE_AND_STORE_STRAIN ){
      // fortran: epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
      epsilondev_xx[tx + working_element*NGLL3] = epsilondev_xx_loc;
      epsilondev_yy[tx + working_element*NGLL3] = epsilondev_yy_loc;
      epsilondev_xy[tx + working_element*NGLL3] = epsilondev_xy_loc;
      epsilondev_xz[tx + working_element*NGLL3] = epsilondev_xz_loc;
      epsilondev_yz[tx + working_element*NGLL3] = epsilondev_yz_loc;
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

void Kernel_2_inner_core(int nb_blocks_to_compute,Mesh* mp,
                         int d_iphase,
                         int* d_ibool,
                         int* d_idoubling,
                         realw* d_xix,realw* d_xiy,realw* d_xiz,
                         realw* d_etax,realw* d_etay,realw* d_etaz,
                         realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                         realw* d_kappav,
                         realw* d_muv,
                         realw* d_epsilondev_xx,
                         realw* d_epsilondev_yy,
                         realw* d_epsilondev_xy,
                         realw* d_epsilondev_xz,
                         realw* d_epsilondev_yz,
                         realw* d_epsilon_trace_over_3,
                         realw* d_one_minus_sum_beta,
                         realw* d_factor_common,
                         realw* d_R_xx,
                         realw* d_R_yy,
                         realw* d_R_xy,
                         realw* d_R_xz,
                         realw* d_R_yz,
                         realw* d_c11store,realw* d_c12store,realw* d_c13store,
                         realw* d_c33store,realw* d_c44store,
                         realw* d_b_epsilondev_xx,
                         realw* d_b_epsilondev_yy,
                         realw* d_b_epsilondev_xy,
                         realw* d_b_epsilondev_xz,
                         realw* d_b_epsilondev_yz,
                         realw* d_b_epsilon_trace_over_3,
                         realw* d_b_R_xx,
                         realw* d_b_R_yy,
                         realw* d_b_R_xy,
                         realw* d_b_R_xz,
                         realw* d_b_R_yz,
                         int FORWARD_OR_ADJOINT){

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("before kernel Kernel_2_inner_core");
#endif

  // if the grid can handle the number of blocks, we let it be 1D
  // grid_2_x = nb_elem_color;
  // nb_elem_color is just how many blocks we are computing now

  int blocksize = NGLL3_PADDED;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(nb_blocks_to_compute,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // Cuda timing
  // cudaEvent_t start, stop;
  // realw time;
  // cudaEventCreate(&start);
  // cudaEventCreate(&stop);
  // cudaEventRecord( start, 0 );
  if( FORWARD_OR_ADJOINT == 1 ){
    Kernel_2_inner_core_impl<<<grid,threads>>>(nb_blocks_to_compute,
                                               mp->NGLOB_INNER_CORE,
                                               d_ibool,
                                               d_idoubling,
                                               mp->d_phase_ispec_inner_inner_core,
                                               mp->num_phase_ispec_inner_core,
                                               d_iphase,
                                               mp->deltat,
                                               mp->use_mesh_coloring_gpu,
                                               mp->d_displ_inner_core,
                                               mp->d_veloc_inner_core,
                                               mp->d_accel_inner_core,
                                               d_xix, d_xiy, d_xiz,
                                               d_etax, d_etay, d_etaz,
                                               d_gammax, d_gammay, d_gammaz,
                                               mp->d_hprime_xx,
                                               mp->d_hprimewgll_xx,
                                               mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                               d_kappav, d_muv,
                                               mp->compute_and_store_strain,
                                               d_epsilondev_xx,
                                               d_epsilondev_yy,
                                               d_epsilondev_xy,
                                               d_epsilondev_xz,
                                               d_epsilondev_yz,
                                               d_epsilon_trace_over_3,
                                               mp->attenuation,
                                               mp->partial_phys_dispersion_only,
                                               mp->use_3d_attenuation_arrays,
                                               d_one_minus_sum_beta,
                                               d_factor_common,
                                               d_R_xx,d_R_yy,d_R_xy,d_R_xz,d_R_yz,
                                               mp->d_alphaval,mp->d_betaval,mp->d_gammaval,
                                               mp->anisotropic_inner_core,
                                               d_c11store,d_c12store,d_c13store,
                                               d_c33store,d_c44store,
                                               mp->gravity,
                                               mp->d_xstore_inner_core,mp->d_ystore_inner_core,mp->d_zstore_inner_core,
                                               mp->d_minus_gravity_table,
                                               mp->d_minus_deriv_gravity_table,
                                               mp->d_density_table,
                                               mp->d_wgll_cube,
                                               mp->NSPEC_INNER_CORE_STRAIN_ONLY,
                                               mp->NSPEC_INNER_CORE);
  }else if( FORWARD_OR_ADJOINT == 3 ){
    // debug
    DEBUG_BACKWARD_FORCES();

    Kernel_2_inner_core_impl<<< grid,threads>>>(nb_blocks_to_compute,
                                                mp->NGLOB_INNER_CORE,
                                                d_ibool,
                                                d_idoubling,
                                                mp->d_phase_ispec_inner_inner_core,
                                                mp->num_phase_ispec_inner_core,
                                                d_iphase,
                                                mp->b_deltat,
                                                mp->use_mesh_coloring_gpu,
                                                mp->d_b_displ_inner_core,
                                                mp->d_b_veloc_inner_core,
                                                mp->d_b_accel_inner_core,
                                                d_xix, d_xiy, d_xiz,
                                                d_etax, d_etay, d_etaz,
                                                d_gammax, d_gammay, d_gammaz,
                                                mp->d_hprime_xx,
                                                mp->d_hprimewgll_xx,
                                                mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                d_kappav, d_muv,
                                                mp->compute_and_store_strain,
                                                d_b_epsilondev_xx,
                                                d_b_epsilondev_yy,
                                                d_b_epsilondev_xy,
                                                d_b_epsilondev_xz,
                                                d_b_epsilondev_yz,
                                                d_b_epsilon_trace_over_3,
                                                mp->attenuation,
                                                mp->partial_phys_dispersion_only,
                                                mp->use_3d_attenuation_arrays,
                                                d_one_minus_sum_beta,
                                                d_factor_common,
                                                d_b_R_xx,d_b_R_yy,d_b_R_xy,d_b_R_xz,d_b_R_yz,
                                                mp->d_b_alphaval,mp->d_b_betaval,mp->d_b_gammaval,
                                                mp->anisotropic_inner_core,
                                                d_c11store,d_c12store,d_c13store,
                                                d_c33store,d_c44store,
                                                mp->gravity,
                                                mp->d_xstore_inner_core,mp->d_ystore_inner_core,mp->d_zstore_inner_core,
                                                mp->d_minus_gravity_table,
                                                mp->d_minus_deriv_gravity_table,
                                                mp->d_density_table,
                                                mp->d_wgll_cube,
                                                mp->NSPEC_INNER_CORE_STRAIN_ONLY,
                                                mp->NSPEC_INNER_CORE);
  }

  // cudaEventRecord( stop, 0 );
  // cudaEventSynchronize( stop );
  // cudaEventElapsedTime( &time, start, stop );
  // cudaEventDestroy( start );
  // cudaEventDestroy( stop );
  // printf("Kernel2 Execution Time: %f ms\n",time);

  /* cudaThreadSynchronize(); */
  /* LOG("Kernel 2 finished"); */
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("Kernel_2_inner_core");
#endif
}

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_forces_inner_core_cuda,
              COMPUTE_FORCES_INNER_CORE_CUDA)(long* Mesh_pointer_f,
                                              int* iphase,
                                              int* FORWARD_OR_ADJOINT_f) {

  TRACE("compute_forces_inner_core_cuda");

  // debug timing
  //printf("Running compute_forces_inner_core_cuda\n");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  int FORWARD_OR_ADJOINT = *FORWARD_OR_ADJOINT_f;

  int num_elements;

  if( *iphase == 1 )
    num_elements = mp->nspec_outer_inner_core;
  else
    num_elements = mp->nspec_inner_inner_core;

  // checks if anything to do
  if( num_elements == 0 ) return;

  // mesh coloring
  if( mp->use_mesh_coloring_gpu ){

    // note: array offsets require sorted arrays, such that e.g. ibool starts with elastic elements
    //         and followed by acoustic ones.
    //         elastic elements also start with outer than inner element ordering

    int nb_colors,nb_blocks_to_compute;
    int istart;
    int offset,offset_nonpadded;
    int offset_nonpadded_att1,offset_nonpadded_att2,offset_nonpadded_att3;
    int offset_nonpadded_strain;
    int offset_ispec;

    // sets up color loop
    if( mp->NSPEC_INNER_CORE > COLORING_MIN_NSPEC_INNER_CORE ){
      if( *iphase == 1 ){
        // outer elements
        nb_colors = mp->num_colors_outer_inner_core;
        istart = 0;

        // array offsets
        offset = 0;
        offset_nonpadded = 0;
        offset_nonpadded_att1 = 0;
        offset_nonpadded_att2 = 0;
        offset_nonpadded_att3 = 0;
        offset_nonpadded_strain = 0;
        offset_ispec = 0;
      }else{
        // inner elements (start after outer elements)
        nb_colors = mp->num_colors_outer_inner_core + mp->num_colors_inner_inner_core;
        istart = mp->num_colors_outer_inner_core;

        // array offsets
        offset = (mp->nspec_outer_inner_core) * NGLL3_PADDED;
        offset_nonpadded = (mp->nspec_outer_inner_core) * NGLL3;
        offset_nonpadded_att1 = (mp->nspec_outer_inner_core) * NGLL3 * N_SLS;
        // for factor_common array
        if( mp->use_3d_attenuation_arrays ){
          offset_nonpadded_att2 = (mp->nspec_outer_inner_core) * NGLL3;
          offset_nonpadded_att3 = (mp->nspec_outer_inner_core) * NGLL3 * N_SLS;
        }else{
          offset_nonpadded_att2 = (mp->nspec_outer_inner_core) * 1;
          offset_nonpadded_att3 = (mp->nspec_outer_inner_core) * 1 * N_SLS;
        }
        // for idoubling array
        offset_ispec = mp->nspec_outer_inner_core;
        // for strain
        if( ! mp->NSPEC_INNER_CORE_STRAIN_ONLY == 1 ){
          offset_nonpadded_strain = (mp->nspec_outer_inner_core) * NGLL3;
        }
      }
    }else{

      // poor element count, only use 1 color per inner/outer run

      if( *iphase == 1 ){
        // outer elements
        nb_colors = 1;
        istart = 0;

        // array offsets
        offset = 0;
        offset_nonpadded = 0;
        offset_nonpadded_att1 = 0;
        offset_nonpadded_att2 = 0;
        offset_nonpadded_att3 = 0;
        offset_nonpadded_strain = 0;
        offset_ispec = 0;
      }else{
        // inner element colors (start after outer elements)
        nb_colors = 1;
        istart = 0;

        // array offsets
        offset = (mp->nspec_outer_inner_core) * NGLL3_PADDED;
        offset_nonpadded = (mp->nspec_outer_inner_core) * NGLL3;
        offset_nonpadded_att1 = (mp->nspec_outer_inner_core) * NGLL3 * N_SLS;
        // for factor_common array
        if( mp->use_3d_attenuation_arrays ){
          offset_nonpadded_att2 = (mp->nspec_outer_inner_core) * NGLL3;
          offset_nonpadded_att3 = (mp->nspec_outer_inner_core) * NGLL3 * N_SLS;
        }else{
          offset_nonpadded_att2 = (mp->nspec_outer_inner_core) * 1;
          offset_nonpadded_att3 = (mp->nspec_outer_inner_core) * 1 * N_SLS;
        }
        // for idoubling array
        offset_ispec = mp->nspec_outer_inner_core;
        // for strain
        if( ! mp->NSPEC_INNER_CORE_STRAIN_ONLY == 1 ){
          offset_nonpadded_strain = (mp->nspec_outer_inner_core) * NGLL3;
        }
      }
    }


    // loops over colors
    for(int icolor = istart; icolor < nb_colors; icolor++){

      // gets number of elements for this color
      if( mp->NSPEC_INNER_CORE > COLORING_MIN_NSPEC_INNER_CORE ){
        nb_blocks_to_compute = mp->h_num_elem_colors_inner_core[icolor];
      }else{
        nb_blocks_to_compute = num_elements;
      }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
      // checks
      if( nb_blocks_to_compute <= 0 ){
        printf("error number of color blocks in inner_core: %d -- color = %d \n",nb_blocks_to_compute,icolor);
        exit(EXIT_FAILURE);
      }
#endif

      Kernel_2_inner_core(nb_blocks_to_compute,mp,
                          *iphase,
                          mp->d_ibool_inner_core + offset_nonpadded,
                          mp->d_idoubling_inner_core + offset_ispec,
                          mp->d_xix_inner_core + offset,mp->d_xiy_inner_core + offset,mp->d_xiz_inner_core + offset,
                          mp->d_etax_inner_core + offset,mp->d_etay_inner_core + offset,mp->d_etaz_inner_core + offset,
                          mp->d_gammax_inner_core + offset,mp->d_gammay_inner_core + offset,mp->d_gammaz_inner_core + offset,
                          mp->d_kappavstore_inner_core + offset,
                          mp->d_muvstore_inner_core + offset,
                          mp->d_epsilondev_xx_inner_core + offset_nonpadded,
                          mp->d_epsilondev_yy_inner_core + offset_nonpadded,
                          mp->d_epsilondev_xy_inner_core + offset_nonpadded,
                          mp->d_epsilondev_xz_inner_core + offset_nonpadded,
                          mp->d_epsilondev_yz_inner_core + offset_nonpadded,
                          mp->d_eps_trace_over_3_inner_core + offset_nonpadded_strain,
                          mp->d_one_minus_sum_beta_inner_core + offset_nonpadded_att2,
                          mp->d_factor_common_inner_core + offset_nonpadded_att3,
                          mp->d_R_xx_inner_core + offset_nonpadded_att1,
                          mp->d_R_yy_inner_core + offset_nonpadded_att1,
                          mp->d_R_xy_inner_core + offset_nonpadded_att1,
                          mp->d_R_xz_inner_core + offset_nonpadded_att1,
                          mp->d_R_yz_inner_core + offset_nonpadded_att1,
                          mp->d_c11store_inner_core + offset,
                          mp->d_c12store_inner_core + offset,
                          mp->d_c13store_inner_core + offset,
                          mp->d_c33store_inner_core + offset,
                          mp->d_c44store_inner_core + offset,
                          mp->d_b_epsilondev_xx_inner_core + offset_nonpadded,
                          mp->d_b_epsilondev_yy_inner_core + offset_nonpadded,
                          mp->d_b_epsilondev_xy_inner_core + offset_nonpadded,
                          mp->d_b_epsilondev_xz_inner_core + offset_nonpadded,
                          mp->d_b_epsilondev_yz_inner_core + offset_nonpadded,
                          mp->d_b_eps_trace_over_3_inner_core + offset_nonpadded,
                          mp->d_b_R_xx_inner_core + offset_nonpadded_att1,
                          mp->d_b_R_yy_inner_core + offset_nonpadded_att1,
                          mp->d_b_R_xy_inner_core + offset_nonpadded_att1,
                          mp->d_b_R_xz_inner_core + offset_nonpadded_att1,
                          mp->d_b_R_yz_inner_core + offset_nonpadded_att1,
                          FORWARD_OR_ADJOINT);

      // for padded and aligned arrays
      offset += nb_blocks_to_compute * NGLL3_PADDED;
      // for no-aligned arrays
      offset_nonpadded += nb_blocks_to_compute * NGLL3;
      offset_nonpadded_att1 += nb_blocks_to_compute * NGLL3 * N_SLS;
      // for factor_common array
      if( mp->use_3d_attenuation_arrays ){
        offset_nonpadded_att2 += nb_blocks_to_compute * NGLL3;
        offset_nonpadded_att3 += nb_blocks_to_compute * NGLL3 * N_SLS;
      }else{
        offset_nonpadded_att2 += nb_blocks_to_compute * 1;
        offset_nonpadded_att3 += nb_blocks_to_compute * 1 * N_SLS;
      }
      // for array(ispec)
      offset_ispec += nb_blocks_to_compute;
      // for strain
      if( ! mp->NSPEC_INNER_CORE_STRAIN_ONLY == 1 ){
        offset_nonpadded_strain += nb_blocks_to_compute * NGLL3;
      }

    } // icolor

  }else{

    // no mesh coloring: uses atomic updates

    Kernel_2_inner_core(num_elements,mp,
                        *iphase,
                        mp->d_ibool_inner_core,
                        mp->d_idoubling_inner_core,
                        mp->d_xix_inner_core,mp->d_xiy_inner_core,mp->d_xiz_inner_core,
                        mp->d_etax_inner_core,mp->d_etay_inner_core,mp->d_etaz_inner_core,
                        mp->d_gammax_inner_core,mp->d_gammay_inner_core,mp->d_gammaz_inner_core,
                        mp->d_kappavstore_inner_core,
                        mp->d_muvstore_inner_core,
                        mp->d_epsilondev_xx_inner_core,
                        mp->d_epsilondev_yy_inner_core,
                        mp->d_epsilondev_xy_inner_core,
                        mp->d_epsilondev_xz_inner_core,
                        mp->d_epsilondev_yz_inner_core,
                        mp->d_eps_trace_over_3_inner_core,
                        mp->d_one_minus_sum_beta_inner_core,
                        mp->d_factor_common_inner_core,
                        mp->d_R_xx_inner_core,
                        mp->d_R_yy_inner_core,
                        mp->d_R_xy_inner_core,
                        mp->d_R_xz_inner_core,
                        mp->d_R_yz_inner_core,
                        mp->d_c11store_inner_core,mp->d_c12store_inner_core,mp->d_c13store_inner_core,
                        mp->d_c33store_inner_core,mp->d_c44store_inner_core,
                        mp->d_b_epsilondev_xx_inner_core,
                        mp->d_b_epsilondev_yy_inner_core,
                        mp->d_b_epsilondev_xy_inner_core,
                        mp->d_b_epsilondev_xz_inner_core,
                        mp->d_b_epsilondev_yz_inner_core,
                        mp->d_b_eps_trace_over_3_inner_core,
                        mp->d_b_R_xx_inner_core,
                        mp->d_b_R_yy_inner_core,
                        mp->d_b_R_xy_inner_core,
                        mp->d_b_R_xz_inner_core,
                        mp->d_b_R_yz_inner_core,
                        FORWARD_OR_ADJOINT);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("compute_forces_inner_core_cuda");
#endif
}

