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


//  cuda constant arrays
__constant__ realw d_hprime_xx[NGLL2];
__constant__ realw d_hprime_yy[NGLL2];
__constant__ realw d_hprime_zz[NGLL2];
__constant__ realw d_hprimewgll_xx[NGLL2];
__constant__ realw d_hprimewgll_yy[NGLL2];
__constant__ realw d_hprimewgll_zz[NGLL2];
__constant__ realw d_wgllwgll_xy[NGLL2];
__constant__ realw d_wgllwgll_xz[NGLL2];
__constant__ realw d_wgllwgll_yz[NGLL2];

__constant__ realw d_wgll_cube[NGLL3]; // needed only for gravity case

/* ----------------------------------------------------------------------------------------------- */

// CONSTANT arrays setup

/* ----------------------------------------------------------------------------------------------- */

/* note:
 constant arrays when used in other compute_forces_***_cuda.cu routines stay zero,
 constant declaration and cudaMemcpyToSymbol would have to be in the same file...

 extern keyword doesn't work for __constant__ declarations.

 also:
 cudaMemcpyToSymbol("deviceCaseParams", caseParams, sizeof(CaseParams));
 ..
 and compile with -arch=sm_20

 see also: http://stackoverflow.com/questions/4008031/how-to-use-cuda-constant-memory-in-a-programmer-pleasant-way
 doesn't seem to work.

 we could keep arrays separated for acoustic and elastic routines...

 workaround:

    for now, we store pointers with cudaGetSymbolAddress() function calls.
    we pass those pointers in all other compute_forces_..() routines

    in this file, we can use the above constant array declarations without need of the pointers.

 */

// constant arrays

void setConst_hprime_xx(realw* array,Mesh* mp)
{

  cudaError_t err = cudaMemcpyToSymbol(d_hprime_xx, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_hprime_xx: %s\n", cudaGetErrorString(err));
    fprintf(stderr, "The problem is maybe -arch sm_13 instead of -arch sm_11 in the Makefile, please doublecheck\n");
    exit(1);
  }

  err = cudaGetSymbolAddress((void**)&(mp->d_hprime_xx),"d_hprime_xx");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_hprime_xx: %s\n", cudaGetErrorString(err));
    exit(1);
  }
}

void setConst_hprime_yy(realw* array,Mesh* mp)
{

  cudaError_t err = cudaMemcpyToSymbol(d_hprime_yy, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_hprime_yy: %s\n", cudaGetErrorString(err));
    fprintf(stderr, "The problem is maybe -arch sm_13 instead of -arch sm_11 in the Makefile, please doublecheck\n");
    exit(1);
  }

  err = cudaGetSymbolAddress((void**)&(mp->d_hprime_yy),"d_hprime_yy");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_hprime_yy: %s\n", cudaGetErrorString(err));
    exit(1);
  }
}

void setConst_hprime_zz(realw* array,Mesh* mp)
{

  cudaError_t err = cudaMemcpyToSymbol(d_hprime_zz, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_hprime_zz: %s\n", cudaGetErrorString(err));
    fprintf(stderr, "The problem is maybe -arch sm_13 instead of -arch sm_11 in the Makefile, please doublecheck\n");
    exit(1);
  }

  err = cudaGetSymbolAddress((void**)&(mp->d_hprime_zz),"d_hprime_zz");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_hprime_zz: %s\n", cudaGetErrorString(err));
    exit(1);
  }
}


void setConst_hprimewgll_xx(realw* array,Mesh* mp)
{
  cudaError_t err = cudaMemcpyToSymbol(d_hprimewgll_xx, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_hprimewgll_xx: %s\n", cudaGetErrorString(err));
    exit(1);
  }

  err = cudaGetSymbolAddress((void**)&(mp->d_hprimewgll_xx),"d_hprimewgll_xx");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_hprimewgll_xx: %s\n", cudaGetErrorString(err));
    exit(1);
  }
}

void setConst_hprimewgll_yy(realw* array,Mesh* mp)
{
  cudaError_t err = cudaMemcpyToSymbol(d_hprimewgll_yy, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_hprimewgll_yy: %s\n", cudaGetErrorString(err));
    exit(1);
  }

  err = cudaGetSymbolAddress((void**)&(mp->d_hprimewgll_yy),"d_hprimewgll_yy");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_hprimewgll_yy: %s\n", cudaGetErrorString(err));
    exit(1);
  }
}

void setConst_hprimewgll_zz(realw* array,Mesh* mp)
{
  cudaError_t err = cudaMemcpyToSymbol(d_hprimewgll_zz, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_hprimewgll_zz: %s\n", cudaGetErrorString(err));
    exit(1);
  }

  err = cudaGetSymbolAddress((void**)&(mp->d_hprimewgll_zz),"d_hprimewgll_zz");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_hprimewgll_zz: %s\n", cudaGetErrorString(err));
    exit(1);
  }
}

void setConst_wgllwgll_xy(realw* array,Mesh* mp)
{
  cudaError_t err = cudaMemcpyToSymbol(d_wgllwgll_xy, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_wgllwgll_xy: %s\n", cudaGetErrorString(err));
    exit(1);
  }
  //mp->d_wgllwgll_xy = d_wgllwgll_xy;
  err = cudaGetSymbolAddress((void**)&(mp->d_wgllwgll_xy),"d_wgllwgll_xy");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_wgllwgll_xy: %s\n", cudaGetErrorString(err));
    exit(1);
  }

}

void setConst_wgllwgll_xz(realw* array,Mesh* mp)
{
  cudaError_t err = cudaMemcpyToSymbol(d_wgllwgll_xz, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in  setConst_wgllwgll_xz: %s\n", cudaGetErrorString(err));
    exit(1);
  }
  //mp->d_wgllwgll_xz = d_wgllwgll_xz;
  err = cudaGetSymbolAddress((void**)&(mp->d_wgllwgll_xz),"d_wgllwgll_xz");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_wgllwgll_xz: %s\n", cudaGetErrorString(err));
    exit(1);
  }

}

void setConst_wgllwgll_yz(realw* array,Mesh* mp)
{
  cudaError_t err = cudaMemcpyToSymbol(d_wgllwgll_yz, array, NGLL2*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_wgllwgll_yz: %s\n", cudaGetErrorString(err));
    exit(1);
  }
  //mp->d_wgllwgll_yz = d_wgllwgll_yz;
  err = cudaGetSymbolAddress((void**)&(mp->d_wgllwgll_yz),"d_wgllwgll_yz");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_wgllwgll_yz: %s\n", cudaGetErrorString(err));
    exit(1);
  }

}

void setConst_wgll_cube(realw* array,Mesh* mp)
{
  cudaError_t err = cudaMemcpyToSymbol(d_wgll_cube, array, NGLL3*sizeof(realw));
  if (err != cudaSuccess)
  {
    fprintf(stderr, "Error in setConst_wgll_cube: %s\n", cudaGetErrorString(err));
    exit(1);
  }
  //mp->d_wgll_cube = d_wgll_cube;
  err = cudaGetSymbolAddress((void**)&(mp->d_wgll_cube),"d_wgll_cube");
  if(err != cudaSuccess) {
    fprintf(stderr, "Error with d_wgll_cube: %s\n", cudaGetErrorString(err));
    exit(1);
  }

}


/* ----------------------------------------------------------------------------------------------- */

// elemental routines

/* ----------------------------------------------------------------------------------------------- */

// updates stress

__device__ void compute_element_cm_att_stress(int tx,int working_element,
                                              realw* R_xx,
                                              realw* R_yy,
                                              realw* R_xy,
                                              realw* R_xz,
                                              realw* R_yz,
                                              reald* sigma_xx,
                                              reald* sigma_yy,
                                              reald* sigma_zz,
                                              reald* sigma_xy,
                                              reald* sigma_xz,
                                              reald* sigma_yz) {

  int i_sls,offset;
  reald R_xx_val,R_yy_val;

  for(i_sls = 0; i_sls < N_SLS; i_sls++){
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
  return;
}

/* ----------------------------------------------------------------------------------------------- */

// updates R_memory

__device__ void compute_element_cm_att_memory(int tx,int working_element,
                                              realw* d_muvstore,
                                              realw* factor_common,
                                              realw* alphaval,realw* betaval,realw* gammaval,
                                              realw* R_xx,realw* R_yy,realw* R_xy,realw* R_xz,realw* R_yz,
                                              realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                              realw* epsilondev_xz,realw* epsilondev_yz,
                                              reald epsilondev_xx_loc,reald epsilondev_yy_loc,reald epsilondev_xy_loc,
                                              reald epsilondev_xz_loc,reald epsilondev_yz_loc,
                                              int ANISOTROPY,
                                              realw* d_c44store
                                              ){

  int i_sls;
  int ijk_ispec;
  int offset_align,offset;
  reald fac;
  reald alphaval_loc,betaval_loc,gammaval_loc;
  reald factor_loc,Sn,Snp1;

  // indices
  offset_align = tx + NGLL3_PADDED * working_element;
  ijk_ispec = tx + NGLL3 * working_element;

  // shear moduli for common factor (only Q_mu attenuation)
  if( ANISOTROPY ){
    fac = d_c44store[offset_align];
  }else{
    fac = d_muvstore[offset_align];
  }

  // use Runge-Kutta scheme to march in time
  for(i_sls = 0; i_sls < N_SLS; i_sls++){
    // indices
    // note: index for R_xx,... here is (i_sls,i,j,k,ispec) and not (i,j,k,ispec,i_sls) as in local version
    //          local version: offset_sls = tx + NGLL3*(working_element + NSPEC*i_sls);
    // index for (i_sls,i,j,k,ispec)
    offset = i_sls + N_SLS*(tx + NGLL3*working_element);
    // index for (i,j,k,ispec,i_sls)
    //offset_sls = tx + NGLL3*(working_element + NSPEC*i_sls);

    // either mustore(i,j,k,ispec) * factor_common(i_sls,i,j,k,ispec)
    // or       factor_common(i_sls,:,:,:,ispec) * c44store(:,:,:,ispec)
    factor_loc = fac * factor_common[offset];

    alphaval_loc = alphaval[i_sls]; // (i_sls)
    betaval_loc = betaval[i_sls];
    gammaval_loc = gammaval[i_sls];

    // term in xx
    Sn   = factor_loc * epsilondev_xx[ijk_ispec]; //(i,j,k,ispec)
    Snp1   = factor_loc * epsilondev_xx_loc; //(i,j,k)
    R_xx[offset] = alphaval_loc * R_xx[offset] + betaval_loc * Sn + gammaval_loc * Snp1;

    // term in yy
    Sn   = factor_loc * epsilondev_yy[ijk_ispec];
    Snp1   = factor_loc * epsilondev_yy_loc;
    R_yy[offset] = alphaval_loc * R_yy[offset] + betaval_loc * Sn + gammaval_loc * Snp1;
    // term in zz not computed since zero trace

    // term in xy
    Sn   = factor_loc * epsilondev_xy[ijk_ispec];
    Snp1   = factor_loc * epsilondev_xy_loc;
    R_xy[offset] = alphaval_loc * R_xy[offset] + betaval_loc * Sn + gammaval_loc * Snp1;

    // term in xz
    Sn   = factor_loc * epsilondev_xz[ijk_ispec];
    Snp1   = factor_loc * epsilondev_xz_loc;
    R_xz[offset] = alphaval_loc * R_xz[offset] + betaval_loc * Sn + gammaval_loc * Snp1;

    // term in yz
    Sn   = factor_loc * epsilondev_yz[ijk_ispec];
    Snp1   = factor_loc * epsilondev_yz_loc;
    R_yz[offset] = alphaval_loc * R_yz[offset] + betaval_loc * Sn + gammaval_loc * Snp1;
  }
  return;
}

/* ----------------------------------------------------------------------------------------------- */

// pre-computes gravity term

__device__ void compute_element_cm_gravity(int tx,int working_element,
                                          int* d_ibool,
                                          realw* d_xstore,realw* d_ystore,realw* d_zstore,
                                          realw* d_minus_gravity_table,
                                          realw* d_minus_deriv_gravity_table,
                                          realw* d_density_table,
                                          realw* wgll_cube,
                                          reald jacobianl,
                                          reald* s_dummyx_loc,
                                          reald* s_dummyy_loc,
                                          reald* s_dummyz_loc,
                                          reald* sigma_xx,
                                          reald* sigma_yy,
                                          reald* sigma_zz,
                                          reald* sigma_xy,
                                          reald* sigma_yx,
                                          reald* sigma_xz,
                                          reald* sigma_zx,
                                          reald* sigma_yz,
                                          reald* sigma_zy,
                                          reald* rho_s_H1,
                                          reald* rho_s_H2,
                                          reald* rho_s_H3){

  reald radius,theta,phi;
  reald cos_theta,sin_theta,cos_phi,sin_phi;
  reald minus_g,minus_dg;
  reald rho;
  reald gxl,gyl,gzl;
  reald minus_g_over_radius,minus_dg_plus_g_over_radius;
  reald cos_theta_sq,sin_theta_sq,cos_phi_sq,sin_phi_sq;
  reald Hxxl,Hyyl,Hzzl,Hxyl,Hxzl,Hyzl;
  reald sx_l,sy_l,sz_l;
  reald factor;

  // R_EARTH_KM is the radius of the bottom of the oceans (radius of Earth in km)
  const reald R_EARTH_KM = 6371.0f;
  // uncomment line below for PREM with oceans
  //const reald R_EARTH_KM = 6368.0f;

  // compute non-symmetric terms for gravity

  // use mesh coordinates to get theta and phi
  // x y z contain r theta phi
  int iglob = d_ibool[working_element*NGLL3 + tx]-1;

  radius = d_xstore[iglob];
  theta = d_ystore[iglob];
  phi = d_zstore[iglob];

  cos_theta = cos(theta);
  sin_theta = sin(theta);
  cos_phi = cos(phi);
  sin_phi = sin(phi);

  // for efficiency replace with lookup table every 100 m in radial direction
  // note: radius in crust mantle should never be zero,
  //          and arrays in C start from 0, thus we need to subtract -1
  int int_radius = rint(radius * R_EARTH_KM * 10.0f ) - 1;

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

  return;
}

/* ----------------------------------------------------------------------------------------------- */

// computes stresses for anisotropic element

__device__ void compute_element_cm_aniso(int offset,
                                         realw* d_c11store,realw* d_c12store,realw* d_c13store,
                                         realw* d_c14store,realw* d_c15store,realw* d_c16store,
                                         realw* d_c22store,realw* d_c23store,realw* d_c24store,
                                         realw* d_c25store,realw* d_c26store,realw* d_c33store,
                                         realw* d_c34store,realw* d_c35store,realw* d_c36store,
                                         realw* d_c44store,realw* d_c45store,realw* d_c46store,
                                         realw* d_c55store,realw* d_c56store,realw* d_c66store,
                                         int ATTENUATION,
                                         reald minus_sum_beta,
                                         reald duxdxl,reald duxdyl,reald duxdzl,
                                         reald duydxl,reald duydyl,reald duydzl,
                                         reald duzdxl,reald duzdyl,reald duzdzl,
                                         reald duxdxl_plus_duydyl,reald duxdxl_plus_duzdzl,reald duydyl_plus_duzdzl,
                                         reald duxdyl_plus_duydxl,reald duzdxl_plus_duxdzl,reald duzdyl_plus_duydzl,
                                         reald* sigma_xx,reald* sigma_yy,reald* sigma_zz,
                                         reald* sigma_xy,reald* sigma_xz,reald* sigma_yz
                                         ){

  reald c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66;
  reald mul;

  c11 = d_c11store[offset];
  c12 = d_c12store[offset];
  c13 = d_c13store[offset];
  c14 = d_c14store[offset];
  c15 = d_c15store[offset];
  c16 = d_c16store[offset];
  c22 = d_c22store[offset];
  c23 = d_c23store[offset];
  c24 = d_c24store[offset];
  c25 = d_c25store[offset];
  c26 = d_c26store[offset];
  c33 = d_c33store[offset];
  c34 = d_c34store[offset];
  c35 = d_c35store[offset];
  c36 = d_c36store[offset];
  c44 = d_c44store[offset];
  c45 = d_c45store[offset];
  c46 = d_c46store[offset];
  c55 = d_c55store[offset];
  c56 = d_c56store[offset];
  c66 = d_c66store[offset];

  // use unrelaxed parameters if attenuation
  if( ATTENUATION){
    mul = c44;
    c11 = c11 + 1.33333333333333333333f * minus_sum_beta * mul;
    c12 = c12 - 0.66666666666666666666f * minus_sum_beta * mul;
    c13 = c13 - 0.66666666666666666666f * minus_sum_beta * mul;
    c22 = c22 + 1.33333333333333333333f * minus_sum_beta * mul;
    c23 = c23 - 0.66666666666666666666f * minus_sum_beta * mul;
    c33 = c33 + 1.33333333333333333333f * minus_sum_beta * mul;
    c44 = c44 + minus_sum_beta * mul;
    c55 = c55 + minus_sum_beta * mul;
    c66 = c66 + minus_sum_beta * mul;
  }

  *sigma_xx = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl +
             c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl;
  *sigma_yy = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl +
             c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl;
  *sigma_zz = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl +
             c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl;
  *sigma_xy = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl +
             c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl;
  *sigma_xz = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl +
             c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl;
  *sigma_yz = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl +
             c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl;
  return;
}

/* ----------------------------------------------------------------------------------------------- */

// computes stresses for isotropic element

__device__ void compute_element_cm_iso(int offset,
                                       realw* d_kappavstore,realw* d_muvstore,
                                       int ATTENUATION,
                                       reald one_minus_sum_beta_use,
                                       reald duxdxl,reald duydyl,reald duzdzl,
                                       reald duxdxl_plus_duydyl,reald duxdxl_plus_duzdzl,reald duydyl_plus_duzdzl,
                                       reald duxdyl_plus_duydxl,reald duzdxl_plus_duxdzl,reald duzdyl_plus_duydzl,
                                       reald* sigma_xx,reald* sigma_yy,reald* sigma_zz,
                                       reald* sigma_xy,reald* sigma_xz,reald* sigma_yz){

  reald lambdal,mul,lambdalplus2mul,kappal;

  // compute elements with an elastic isotropic rheology
  kappal = d_kappavstore[offset];
  mul = d_muvstore[offset];

  // use unrelaxed parameters if attenuation
  if( ATTENUATION ){ mul = mul * one_minus_sum_beta_use; }

  lambdalplus2mul = kappal + 1.33333333333333333333f * mul;  // 4./3. = 1.3333333
  lambdal = lambdalplus2mul - 2.0f * mul;

  // compute the six components of the stress tensor sigma
  *sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl;
  *sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl;
  *sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl;

  *sigma_xy = mul*duxdyl_plus_duydxl;
  *sigma_xz = mul*duzdxl_plus_duxdzl;
  *sigma_yz = mul*duzdyl_plus_duydzl;

}

/* ----------------------------------------------------------------------------------------------- */

// computes stresses for transversely isotropic element

__device__ void compute_element_cm_tiso(int offset,
                                        realw* d_kappavstore,realw* d_muvstore,
                                        realw* d_kappahstore,realw* d_muhstore,realw* d_eta_anisostore,
                                        int ATTENUATION,
                                        reald one_minus_sum_beta_use,
                                        reald duxdxl,reald duxdyl,reald duxdzl,
                                        reald duydxl,reald duydyl,reald duydzl,
                                        reald duzdxl,reald duzdyl,reald duzdzl,
                                        reald duxdxl_plus_duydyl,reald duxdxl_plus_duzdzl,reald duydyl_plus_duzdzl,
                                        reald duxdyl_plus_duydxl,reald duzdxl_plus_duxdzl,reald duzdyl_plus_duydzl,
                                        int iglob,int NGLOB,
                                        realw* d_ystore, realw* d_zstore,
                                        reald* sigma_xx,reald* sigma_yy,reald* sigma_zz,
                                        reald* sigma_xy,reald* sigma_xz,reald* sigma_yz){

  reald kappavl,muvl,kappahl,muhl;
  reald rhovpvsq,rhovphsq,rhovsvsq,rhovshsq,eta_aniso;
  reald costheta,sintheta,cosphi,sinphi;
  reald costhetasq,sinthetasq,cosphisq,sinphisq,costhetafour,sinthetafour,cosphifour,sinphifour;
  reald costwotheta,sintwotheta,costwophi,sintwophi,cosfourtheta,cosfourphi;
  reald costwothetasq,costwophisq,sintwophisq;
  reald etaminone,twoetaminone;
  reald two_eta_aniso,four_eta_aniso,six_eta_aniso;
  reald two_rhovsvsq,two_rhovshsq; // two_rhovpvsq,two_rhovphsq
  reald four_rhovsvsq,four_rhovshsq; // four_rhovpvsq,four_rhovphsq
  reald c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66;

  // cosine and sine function in CUDA only supported for float
  reald theta,phi;

  // use Kappa and mu from transversely isotropic model
  kappavl = d_kappavstore[offset];
  muvl = d_muvstore[offset];

  kappahl = d_kappahstore[offset];
  muhl = d_muhstore[offset];

  // use unrelaxed parameters if attenuation
  // eta does not need to be shifted since it is a ratio
  if( ATTENUATION ){
    muvl = muvl * one_minus_sum_beta_use;
    muhl = muhl * one_minus_sum_beta_use;
  }

  rhovpvsq = kappavl + 1.33333333333333333333f * muvl ; //!!! that is C
  rhovphsq = kappahl + 1.33333333333333333333f * muhl ; //!!! that is A

  rhovsvsq = muvl; // !!! that is L
  rhovshsq = muhl; //!!! that is N

  eta_aniso = d_eta_anisostore[offset]; // !!! that is  F / (A - 2 L)

  // use mesh coordinates to get theta and phi
  //ystore and zstore contain theta and phi
  theta = d_ystore[iglob];
  phi = d_zstore[iglob];

  if( sizeof( theta ) == sizeof( float ) ){
    // float operations

    // daniel: TODO - sincos function
    // or: sincosf(theta, &sintheta, &costheta);
    // or with loss of accuracy:  __sincosf(theta, &sintheta, &costheta);
    // or compile with: -use_fast_math

    costheta = cosf(theta);
    sintheta = sinf(theta);

    cosphi = cosf(phi);
    sinphi = sinf(phi);

    costwotheta = cosf(2.0f * theta);
    sintwotheta = sinf(2.0f * theta);
    costwophi = cosf(2.0f * phi);
    sintwophi = sinf(2.0f * phi);
    cosfourtheta = cosf(4.0f * theta);
    cosfourphi = cosf(4.0f * phi);
  }else{
    // double operations
    costheta = cos(theta);
    sintheta = sin(theta);

    cosphi = cos(phi);
    sinphi = sin(phi);

    costwotheta = cos(2.0f * theta);
    sintwotheta = sin(2.0f * theta);
    costwophi = cos(2.0f * phi);
    sintwophi = sin(2.0f * phi);

    cosfourtheta = cos(4.0f * theta);
    cosfourphi = cos(4.0f * phi);
  }

  costhetasq = costheta * costheta;
  sinthetasq = sintheta * sintheta;
  cosphisq = cosphi * cosphi;
  sinphisq = sinphi * sinphi;

  costhetafour = costhetasq * costhetasq;
  sinthetafour = sinthetasq * sinthetasq;
  cosphifour = cosphisq * cosphisq;
  sinphifour = sinphisq * sinphisq;

  costwothetasq = costwotheta * costwotheta;

  costwophisq = costwophi * costwophi;
  sintwophisq = sintwophi * sintwophi;

  etaminone = eta_aniso - 1.0f;
  twoetaminone = 2.0f * eta_aniso - 1.0f;

  // precompute some products to reduce the CPU time

  two_eta_aniso = 2.0f * eta_aniso;
  four_eta_aniso = 4.0f * eta_aniso;
  six_eta_aniso = 6.0f * eta_aniso;

  //two_rhovpvsq = 2.0f * rhovpvsq;
  //two_rhovphsq = 2.0f * rhovphsq;
  two_rhovsvsq = 2.0f * rhovsvsq;
  two_rhovshsq = 2.0f * rhovshsq;

  //four_rhovpvsq = 4.0f * rhovpvsq;
  //four_rhovphsq = 4.0f * rhovphsq;
  four_rhovsvsq = 4.0f * rhovsvsq;
  four_rhovshsq = 4.0f * rhovshsq;

  // the 21 anisotropic coefficients computed using Mathematica

  c11 = rhovphsq*sinphifour + 2.0f*cosphisq*sinphisq*
        (rhovphsq*costhetasq + (eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)*
        sinthetasq) + cosphifour*
        (rhovphsq*costhetafour + 2.0f*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)*
        costhetasq*sinthetasq + rhovpvsq*sinthetafour);

  c12 = ((rhovphsq - two_rhovshsq)*(3.0f + cosfourphi)*costhetasq)*0.25f -
        four_rhovshsq*cosphisq*costhetasq*sinphisq +
        (rhovphsq*(11.0f + 4.0f*costwotheta + cosfourtheta)*sintwophisq)*0.03125f +
        eta_aniso*(rhovphsq - two_rhovsvsq)*(cosphifour +
        2.0f*cosphisq*costhetasq*sinphisq + sinphifour)*sinthetasq +
        rhovpvsq*cosphisq*sinphisq*sinthetafour -
        rhovsvsq*sintwophisq*sinthetafour;

  c13 = (cosphisq*(rhovphsq + six_eta_aniso*rhovphsq + rhovpvsq - four_rhovsvsq -
        12.0f*eta_aniso*rhovsvsq + (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq -
        four_eta_aniso*rhovsvsq)*cosfourtheta))*0.125f +
        sinphisq*(eta_aniso*(rhovphsq - two_rhovsvsq)*costhetasq +
        (rhovphsq - two_rhovshsq)*sinthetasq);

  c14 = costheta*sinphi*((cosphisq*
        (-rhovphsq + rhovpvsq + four_rhovshsq - four_rhovsvsq +
        (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq -
        four_eta_aniso*rhovsvsq)*costwotheta))*0.5f +
        (etaminone*rhovphsq + 2.0f*(rhovshsq - eta_aniso*rhovsvsq))*sinphisq)* sintheta;

  c15 = cosphi*costheta*((cosphisq* (-rhovphsq + rhovpvsq +
        (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)*
        costwotheta))*0.5f + etaminone*(rhovphsq - two_rhovsvsq)*sinphisq)*sintheta;

  c16 = (cosphi*sinphi*(cosphisq* (-rhovphsq + rhovpvsq +
        (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq -
        four_eta_aniso*rhovsvsq)*costwotheta) +
        2.0f*etaminone*(rhovphsq - two_rhovsvsq)*sinphisq)*sinthetasq)*0.5f;

  c22 = rhovphsq*cosphifour + 2.0f*cosphisq*sinphisq*
        (rhovphsq*costhetasq + (eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)*
        sinthetasq) + sinphifour*
        (rhovphsq*costhetafour + 2.0f*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)*
        costhetasq*sinthetasq + rhovpvsq*sinthetafour);

  c23 = ((rhovphsq + six_eta_aniso*rhovphsq + rhovpvsq - four_rhovsvsq - 12.0f*eta_aniso*rhovsvsq +
        (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)*
        cosfourtheta)*sinphisq)*0.125f +
        cosphisq*(eta_aniso*(rhovphsq - two_rhovsvsq)*costhetasq +
        (rhovphsq - two_rhovshsq)*sinthetasq);

  c24 = costheta*sinphi*(etaminone*(rhovphsq - two_rhovsvsq)*cosphisq +
        ((-rhovphsq + rhovpvsq + (twoetaminone*rhovphsq - rhovpvsq +
        four_rhovsvsq - four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)*0.5f)*sintheta;

  c25 = cosphi*costheta*((etaminone*rhovphsq + 2.0f*(rhovshsq - eta_aniso*rhovsvsq))*
        cosphisq + ((-rhovphsq + rhovpvsq + four_rhovshsq - four_rhovsvsq +
        (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq -
        four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)*0.5f)*sintheta;

  c26 = (cosphi*sinphi*(2.0f*etaminone*(rhovphsq - two_rhovsvsq)*cosphisq +
        (-rhovphsq + rhovpvsq + (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq -
        four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)*sinthetasq)*0.5f;

  c33 = rhovpvsq*costhetafour + 2.0f*(eta_aniso*(rhovphsq - two_rhovsvsq) + two_rhovsvsq)*
        costhetasq*sinthetasq + rhovphsq*sinthetafour;

  c34 = -((rhovphsq - rhovpvsq + (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq
        - four_eta_aniso*rhovsvsq)*costwotheta)*sinphi*sintwotheta)*0.25f;

  c35 = -(cosphi*(rhovphsq - rhovpvsq +
        (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)*
        costwotheta)*sintwotheta)*0.25f;

  c36 = -((rhovphsq - rhovpvsq - four_rhovshsq + four_rhovsvsq +
        (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)*
        costwotheta)*sintwophi*sinthetasq)*0.25f;

  c44 = cosphisq*(rhovsvsq*costhetasq + rhovshsq*sinthetasq) +
        sinphisq*(rhovsvsq*costwothetasq +
        (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + four_eta_aniso*rhovsvsq)*costhetasq* sinthetasq);

  c45 = ((rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq - two_rhovshsq - two_rhovsvsq +
        four_eta_aniso*rhovsvsq + (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq +
        4.0f*etaminone*rhovsvsq)*costwotheta)*sintwophi*sinthetasq)*0.25f;

  c46 = -(cosphi*costheta*((rhovshsq - rhovsvsq)*cosphisq -
        ((rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq - two_rhovshsq - two_rhovsvsq +
        four_eta_aniso*rhovsvsq + (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq +
        four_rhovsvsq - four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)*0.5f)* sintheta);

  c55 = sinphisq*(rhovsvsq*costhetasq + rhovshsq*sinthetasq) +
        cosphisq*(rhovsvsq*costwothetasq +
        (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + four_eta_aniso*rhovsvsq)*costhetasq* sinthetasq);

  c56 = costheta*sinphi*((cosphisq*
        (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq - two_rhovshsq - two_rhovsvsq +
        four_eta_aniso*rhovsvsq + (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq +
        four_rhovsvsq - four_eta_aniso*rhovsvsq)*costwotheta))*0.5f +
        (-rhovshsq + rhovsvsq)*sinphisq)*sintheta;

  c66 = rhovshsq*costwophisq*costhetasq -
        2.0f*(rhovphsq - two_rhovshsq)*cosphisq*costhetasq*sinphisq +
        (rhovphsq*(11.0f + 4.0f*costwotheta + cosfourtheta)*sintwophisq)*0.03125f -
        (rhovsvsq*(-6.0f - 2.0f*cosfourphi + cos(4.0f*phi - 2.0f*theta) - 2.0f*costwotheta +
        cos(2.0f*(2.0f*phi + theta)))*sinthetasq)*0.125f +
        rhovpvsq*cosphisq*sinphisq*sinthetafour -
        (eta_aniso*(rhovphsq - two_rhovsvsq)*sintwophisq*sinthetafour)*0.5f;

  // general expression of stress tensor for full Cijkl with 21 coefficients

  *sigma_xx = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl +
              c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl;

  *sigma_yy = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl +
              c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl;

  *sigma_zz = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl +
              c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl;

  *sigma_xy = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl +
              c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl;

  *sigma_xz = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl +
              c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl;

  *sigma_yz = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl +
              c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl;

  return;
}



/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2
//
// for crust_mantle

/* ----------------------------------------------------------------------------------------------- */

__global__ void Kernel_2_crust_mantle_impl(int nb_blocks_to_compute,
                                          int NGLOB,
                                          int* d_ibool,
                                          int* d_ispec_is_tiso,
                                          int* d_phase_ispec_inner,
                                          int num_phase_ispec,
                                          int d_iphase,
                                          realw d_deltat,
                                          int use_mesh_coloring_gpu,
                                          realw* d_displ,
                                          realw* d_veloc,
                                          realw* d_accel,
                                          realw* d_xix, realw* d_xiy, realw* d_xiz,
                                          realw* d_etax, realw* d_etay, realw* d_etaz,
                                          realw* d_gammax, realw* d_gammay, realw* d_gammaz,
                                          realw* d_kappavstore, realw* d_muvstore,
                                          realw* d_kappahstore, realw* d_muhstore,
                                          realw* d_eta_anisostore,
                                          int COMPUTE_AND_STORE_STRAIN,
                                          realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                          realw* epsilondev_xz,realw* epsilondev_yz,
                                          realw* epsilon_trace_over_3,
                                          int SIMULATION_TYPE,
                                          int ATTENUATION,
                                          int ATTENUATION_NEW,
                                          int USE_ATTENUATION_MIMIC,
                                          realw* one_minus_sum_beta,realw* factor_common,
                                          realw* R_xx, realw* R_yy, realw* R_xy, realw* R_xz, realw* R_yz,
                                          realw* alphaval,realw* betaval,realw* gammaval,
                                          int ANISOTROPY,
                                          realw* d_c11store,realw* d_c12store,realw* d_c13store,
                                          realw* d_c14store,realw* d_c15store,realw* d_c16store,
                                          realw* d_c22store,realw* d_c23store,realw* d_c24store,
                                          realw* d_c25store,realw* d_c26store,realw* d_c33store,
                                          realw* d_c34store,realw* d_c35store,realw* d_c36store,
                                          realw* d_c44store,realw* d_c45store,realw* d_c46store,
                                          realw* d_c55store,realw* d_c56store,realw* d_c66store,
                                          int GRAVITY,
                                          realw* d_xstore,realw* d_ystore,realw* d_zstore,
                                          realw* d_minus_gravity_table,
                                          realw* d_minus_deriv_gravity_table,
                                          realw* d_density_table,
                                          realw* wgll_cube,
                                          int NSPEC_CRUST_MANTLE_STRAIN_ONLY){

  /* int bx = blockIdx.y*blockDim.x+blockIdx.x; //possible bug in original code*/
  int bx = blockIdx.y*gridDim.x+blockIdx.x;
  /* int bx = blockIdx.x; */
  int tx = threadIdx.x;

  //const int NGLLX = 5;
  // const int NGLL2 = 25;
  //const int NGLL3 = NGLL3;
  const int NGLL3_ALIGN = NGLL3_PADDED;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int active,offset;
  int iglob = 0;
  int working_element;

  reald tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  reald xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  reald duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  reald duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  reald duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;

  reald tempx1l_att,tempx2l_att,tempx3l_att,tempy1l_att,tempy2l_att,tempy3l_att,tempz1l_att,tempz2l_att,tempz3l_att;
  reald duxdxl_att,duxdyl_att,duxdzl_att,duydxl_att,duydyl_att,duydzl_att,duzdxl_att,duzdyl_att,duzdzl_att;
  reald duxdyl_plus_duydxl_att,duzdxl_plus_duxdzl_att,duzdyl_plus_duydzl_att;

  reald fac1,fac2,fac3;
  reald minus_sum_beta,one_minus_sum_beta_use;

  reald sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  reald epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc;
  //reald c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66;
  reald sum_terms1,sum_terms2,sum_terms3;

  // gravity variables
  reald sigma_yx,sigma_zx,sigma_zy;
  reald rho_s_H1,rho_s_H2,rho_s_H3;

#ifndef MANUALLY_UNROLLED_LOOPS
    int l;
    realw hp1,hp2,hp3;
#endif

    __shared__ reald s_dummyx_loc[NGLL3];
    __shared__ reald s_dummyy_loc[NGLL3];
    __shared__ reald s_dummyz_loc[NGLL3];

    __shared__ reald s_dummyx_loc_att[NGLL3];
    __shared__ reald s_dummyy_loc_att[NGLL3];
    __shared__ reald s_dummyz_loc_att[NGLL3];

    __shared__ reald s_tempx1[NGLL3];
    __shared__ reald s_tempx2[NGLL3];
    __shared__ reald s_tempx3[NGLL3];
    __shared__ reald s_tempy1[NGLL3];
    __shared__ reald s_tempy2[NGLL3];
    __shared__ reald s_tempy3[NGLL3];
    __shared__ reald s_tempz1[NGLL3];
    __shared__ reald s_tempz2[NGLL3];
    __shared__ reald s_tempz3[NGLL3];

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

      // iglob = d_ibool[working_element*NGLL3_ALIGN + tx]-1;
      iglob = d_ibool[working_element*NGLL3 + tx]-1;

#ifdef USE_TEXTURES
      s_dummyx_loc[tx] = tex1Dfetch(tex_displ, iglob);
      s_dummyy_loc[tx] = tex1Dfetch(tex_displ, iglob + NGLOB);
      s_dummyz_loc[tx] = tex1Dfetch(tex_displ, iglob + 2*NGLOB);
#else
      // changing iglob indexing to match fortran row changes fast style
      s_dummyx_loc[tx] = d_displ[iglob*3];
      s_dummyy_loc[tx] = d_displ[iglob*3 + 1];
      s_dummyz_loc[tx] = d_displ[iglob*3 + 2];
#endif

      if(ATTENUATION){
	if(ATTENUATION_NEW){
	  // takes new routines
	  // use first order Taylor expansion of displacement for local storage of stresses 
	  // at this current time step, to fix attenuation in a consistent way
#ifdef USE_TEXTURES
	  s_dummyx_loc_att[tx] = s_dummyx_loc[tx] + d_deltat * tex1Dfetch(tex_veloc, iglob);
	  s_dummyy_loc_att[tx] = s_dummyy_loc[tx] + d_deltat * tex1Dfetch(tex_veloc, iglob + NGLOB);
	  s_dummyz_loc_att[tx] = s_dummyz_loc[tx] + d_deltat * tex1Dfetch(tex_veloc, iglob + 2*NGLOB);
#else
	  s_dummyx_loc_att[tx] = s_dummyx_loc[tx] + d_deltat * d_veloc[iglob*3];
	  s_dummyy_loc_att[tx] = s_dummyy_loc[tx] + d_deltat * d_veloc[iglob*3 + 1];
	  s_dummyz_loc_att[tx] = s_dummyz_loc[tx] + d_deltat * d_veloc[iglob*3 + 2];
#endif
	}
	else{
	  // takes old routines
	  s_dummyx_loc_att[tx] = s_dummyx_loc[tx];
	  s_dummyy_loc_att[tx] = s_dummyy_loc[tx];
	  s_dummyz_loc_att[tx] = s_dummyz_loc[tx];
	}
      }
    }

// synchronize all the threads (one thread for each of the NGLL grid points of the
// current spectral element) because we need the whole element to be ready in order
// to be able to compute the matrix products along cut planes of the 3D element below
    __syncthreads();

#ifndef MAKE_KERNEL2_BECOME_STUPID_FOR_TESTS

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
          hp1 = d_hprime_xx[l*NGLLX+I];
          offset = K*NGLL2+J*NGLLX+l;
          tempx1l += s_dummyx_loc[offset]*hp1;
          tempy1l += s_dummyy_loc[offset]*hp1;
          tempz1l += s_dummyz_loc[offset]*hp1;

          hp2 = d_hprime_xx[l*NGLLX+J];
          offset = K*NGLL2+l*NGLLX+I;
          tempx2l += s_dummyx_loc[offset]*hp2;
          tempy2l += s_dummyy_loc[offset]*hp2;
          tempz2l += s_dummyz_loc[offset]*hp2;

          hp3 = d_hprime_xx[l*NGLLX+K];
          offset = l*NGLL2+J*NGLLX+I;
          tempx3l += s_dummyx_loc[offset]*hp3;
          tempy3l += s_dummyy_loc[offset]*hp3;
          tempz3l += s_dummyz_loc[offset]*hp3;

      }


      if( ATTENUATION){
	// temporary variables used for fixing attenuation in a consistent way

	tempx1l_att = 0.f;
	tempx2l_att = 0.f;
	tempx3l_att = 0.f;

	tempy1l_att = 0.f;
	tempy2l_att = 0.f;
	tempy3l_att = 0.f;

	tempz1l_att = 0.f;
	tempz2l_att = 0.f;
	tempz3l_att = 0.f;

	for (l=0;l<NGLLX;l++) {
          hp1 = d_hprime_xx[l*NGLLX+I];
          offset = K*NGLL2+J*NGLLX+l;
          tempx1l_att += s_dummyx_loc_att[offset]*hp1;
          tempy1l_att += s_dummyy_loc_att[offset]*hp1;
          tempz1l_att += s_dummyz_loc_att[offset]*hp1;

          hp2 = d_hprime_xx[l*NGLLX+J];
          offset = K*NGLL2+l*NGLLX+I;
          tempx2l_att += s_dummyx_loc_att[offset]*hp2;
          tempy2l_att += s_dummyy_loc_att[offset]*hp2;
          tempz2l_att += s_dummyz_loc_att[offset]*hp2;

          hp3 = d_hprime_xx[l*NGLLX+K];
          offset = l*NGLL2+J*NGLLX+I;
          tempx3l_att += s_dummyx_loc_att[offset]*hp3;
          tempy3l_att += s_dummyy_loc_att[offset]*hp3;
          tempz3l_att += s_dummyz_loc_att[offset]*hp3;

	}
      }
#else

      tempx1l = s_dummyx_loc[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
              + s_dummyx_loc[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
              + s_dummyx_loc[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
              + s_dummyx_loc[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
              + s_dummyx_loc[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

      tempy1l = s_dummyy_loc[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
              + s_dummyy_loc[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
              + s_dummyy_loc[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
              + s_dummyy_loc[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
              + s_dummyy_loc[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

      tempz1l = s_dummyz_loc[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
              + s_dummyz_loc[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
              + s_dummyz_loc[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
              + s_dummyz_loc[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
              + s_dummyz_loc[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

      tempx2l = s_dummyx_loc[K*NGLL2+I]*d_hprime_xx[J]
              + s_dummyx_loc[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
              + s_dummyx_loc[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
              + s_dummyx_loc[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
              + s_dummyx_loc[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

      tempy2l = s_dummyy_loc[K*NGLL2+I]*d_hprime_xx[J]
              + s_dummyy_loc[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
              + s_dummyy_loc[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
              + s_dummyy_loc[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
              + s_dummyy_loc[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

      tempz2l = s_dummyz_loc[K*NGLL2+I]*d_hprime_xx[J]
              + s_dummyz_loc[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
              + s_dummyz_loc[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
              + s_dummyz_loc[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
              + s_dummyz_loc[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

      tempx3l = s_dummyx_loc[J*NGLLX+I]*d_hprime_xx[K]
              + s_dummyx_loc[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
              + s_dummyx_loc[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
              + s_dummyx_loc[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
              + s_dummyx_loc[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

      tempy3l = s_dummyy_loc[J*NGLLX+I]*d_hprime_xx[K]
              + s_dummyy_loc[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
              + s_dummyy_loc[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
              + s_dummyy_loc[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
              + s_dummyy_loc[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

      tempz3l = s_dummyz_loc[J*NGLLX+I]*d_hprime_xx[K]
              + s_dummyz_loc[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
              + s_dummyz_loc[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
              + s_dummyz_loc[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
              + s_dummyz_loc[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

      if( ATTENUATION){
	// temporary variables used for fixing attenuation in a consistent way

	tempx1l_att = s_dummyx_loc_att[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
	  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
	  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
	  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
	  + s_dummyx_loc_att[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

	tempy1l_att = s_dummyy_loc_att[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
	  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
	  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
	  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
	  + s_dummyy_loc_att[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

	tempz1l_att = s_dummyz_loc_att[K*NGLL2+J*NGLLX]*d_hprime_xx[I]
	  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+1]*d_hprime_xx[NGLLX+I]
	  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+2]*d_hprime_xx[2*NGLLX+I]
	  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+3]*d_hprime_xx[3*NGLLX+I]
	  + s_dummyz_loc_att[K*NGLL2+J*NGLLX+4]*d_hprime_xx[4*NGLLX+I];

	tempx2l_att = s_dummyx_loc_att[K*NGLL2+I]*d_hprime_xx[J]
	  + s_dummyx_loc_att[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
	  + s_dummyx_loc_att[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
	  + s_dummyx_loc_att[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
	  + s_dummyx_loc_att[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

	tempy2l_att = s_dummyy_loc_att[K*NGLL2+I]*d_hprime_xx[J]
	  + s_dummyy_loc_att[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
	  + s_dummyy_loc_att[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
	  + s_dummyy_loc_att[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
	  + s_dummyy_loc_att[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

	tempz2l_att = s_dummyz_loc_att[K*NGLL2+I]*d_hprime_xx[J]
	  + s_dummyz_loc_att[K*NGLL2+NGLLX+I]*d_hprime_xx[NGLLX+J]
	  + s_dummyz_loc_att[K*NGLL2+2*NGLLX+I]*d_hprime_xx[2*NGLLX+J]
	  + s_dummyz_loc_att[K*NGLL2+3*NGLLX+I]*d_hprime_xx[3*NGLLX+J]
	  + s_dummyz_loc_att[K*NGLL2+4*NGLLX+I]*d_hprime_xx[4*NGLLX+J];

	tempx3l_att = s_dummyx_loc_att[J*NGLLX+I]*d_hprime_xx[K]
	  + s_dummyx_loc_att[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
	  + s_dummyx_loc_att[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
	  + s_dummyx_loc_att[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
	  + s_dummyx_loc_att[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

	tempy3l_att = s_dummyy_loc_att[J*NGLLX+I]*d_hprime_xx[K]
	  + s_dummyy_loc_att[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
	  + s_dummyy_loc_att[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
	  + s_dummyy_loc_att[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
	  + s_dummyy_loc_att[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];

	tempz3l_att = s_dummyz_loc_att[J*NGLLX+I]*d_hprime_xx[K]
	  + s_dummyz_loc_att[NGLL2+J*NGLLX+I]*d_hprime_xx[NGLLX+K]
	  + s_dummyz_loc_att[2*NGLL2+J*NGLLX+I]*d_hprime_xx[2*NGLLX+K]
	  + s_dummyz_loc_att[3*NGLL2+J*NGLLX+I]*d_hprime_xx[3*NGLLX+K]
	  + s_dummyz_loc_att[4*NGLL2+J*NGLLX+I]*d_hprime_xx[4*NGLLX+K];
      }

#endif

      // compute derivatives of ux, uy and uz with respect to x, y and z
      offset = working_element*NGLL3_ALIGN + tx;

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

      if( ATTENUATION){
	// temporary variables used for fixing attenuation in a consistent way

	duxdxl_att = xixl*tempx1l_att + etaxl*tempx2l_att + gammaxl*tempx3l_att;
	duxdyl_att = xiyl*tempx1l_att + etayl*tempx2l_att + gammayl*tempx3l_att;
	duxdzl_att = xizl*tempx1l_att + etazl*tempx2l_att + gammazl*tempx3l_att;

	duydxl_att = xixl*tempy1l_att + etaxl*tempy2l_att + gammaxl*tempy3l_att;
	duydyl_att = xiyl*tempy1l_att + etayl*tempy2l_att + gammayl*tempy3l_att;
	duydzl_att = xizl*tempy1l_att + etazl*tempy2l_att + gammazl*tempy3l_att;

	duzdxl_att = xixl*tempz1l_att + etaxl*tempz2l_att + gammaxl*tempz3l_att;
	duzdyl_att = xiyl*tempz1l_att + etayl*tempz2l_att + gammayl*tempz3l_att;
	duzdzl_att = xizl*tempz1l_att + etazl*tempz2l_att + gammazl*tempz3l_att;

	// precompute some sums to save CPU time
	duxdyl_plus_duydxl_att = duxdyl_att + duydxl_att;
	duzdxl_plus_duxdzl_att = duzdxl_att + duxdzl_att;
	duzdyl_plus_duydzl_att = duzdyl_att + duydzl_att;

	// computes deviatoric strain attenuation and/or for kernel calculations
	if(COMPUTE_AND_STORE_STRAIN) {
	  realw templ = 0.33333333333333333333f * (duxdxl_att + duydyl_att + duzdzl_att); // 1./3. = 0.33333

	  // local storage: stresses at this current time step
	  epsilondev_xx_loc = duxdxl_att - templ;
	  epsilondev_yy_loc = duydyl_att - templ;
	  epsilondev_xy_loc = 0.5f * duxdyl_plus_duydxl_att;
	  epsilondev_xz_loc = 0.5f * duzdxl_plus_duxdzl_att;
	  epsilondev_yz_loc = 0.5f * duzdyl_plus_duydzl_att;

	  if(NSPEC_CRUST_MANTLE_STRAIN_ONLY == 1) {
	    epsilon_trace_over_3[tx] = templ;
	  }else{
	    epsilon_trace_over_3[tx + working_element*NGLL3] = templ;
	  }
	}
      }else{
	// computes deviatoric strain attenuation and/or for kernel calculations
	if(COMPUTE_AND_STORE_STRAIN) {
	  realw templ = 0.33333333333333333333f * (duxdxl + duydyl + duzdzl); // 1./3. = 0.33333

	  // local storage: stresses at this current time step
	  epsilondev_xx_loc = duxdxl - templ;
	  epsilondev_yy_loc = duydyl - templ;
	  epsilondev_xy_loc = 0.5f * duxdyl_plus_duydxl;
	  epsilondev_xz_loc = 0.5f * duzdxl_plus_duxdzl;
	  epsilondev_yz_loc = 0.5f * duzdyl_plus_duydzl;

	  if(NSPEC_CRUST_MANTLE_STRAIN_ONLY == 1) {
	    epsilon_trace_over_3[tx] = templ;
	  }else{
	    epsilon_trace_over_3[tx + working_element*NGLL3] = templ;
	  }
	}
      }

      // attenuation
      if(ATTENUATION){
        // use unrelaxed parameters if attenuation
        one_minus_sum_beta_use = one_minus_sum_beta[tx+working_element*NGLL3]; // (i,j,k,ispec)
        minus_sum_beta = one_minus_sum_beta_use - 1.0f;
      }

      // computes stresses
      if(ANISOTROPY){
        // full anisotropic case, stress calculations
        compute_element_cm_aniso(offset,
                              d_c11store,d_c12store,d_c13store,d_c14store,d_c15store,d_c16store,d_c22store,
                              d_c23store,d_c24store,d_c25store,d_c26store,d_c33store,d_c34store,d_c35store,
                              d_c36store,d_c44store,d_c45store,d_c46store,d_c55store,d_c56store,d_c66store,
                              ATTENUATION,
                              minus_sum_beta,
                              duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl,
                              duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl,
                              duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl,
                              &sigma_xx,&sigma_yy,&sigma_zz,
                              &sigma_xy,&sigma_xz,&sigma_yz);

      }else{
        if( ! d_ispec_is_tiso[working_element] ){
          // isotropic case
          compute_element_cm_iso(offset,
                              d_kappavstore,d_muvstore,
                              ATTENUATION,
                              one_minus_sum_beta_use,
                              duxdxl,duydyl,duzdzl,
                              duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl,
                              duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl,
                              &sigma_xx,&sigma_yy,&sigma_zz,
                              &sigma_xy,&sigma_xz,&sigma_yz);
        }else{
          // transverse isotropy
          compute_element_cm_tiso(offset,
                                d_kappavstore,d_muvstore,
                                d_kappahstore,d_muhstore,d_eta_anisostore,
                                ATTENUATION,
                                one_minus_sum_beta_use,
                                duxdxl,duxdyl,duxdzl,
                                duydxl,duydyl,duydzl,
                                duzdxl,duzdyl,duzdzl,
                                duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl,
                                duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl,
                                iglob, NGLOB,
                                d_ystore,d_zstore,
                                &sigma_xx,&sigma_yy,&sigma_zz,
                                &sigma_xy,&sigma_xz,&sigma_yz);
        }
      } // ! end of test whether isotropic or anisotropic element


      if(ATTENUATION && (! USE_ATTENUATION_MIMIC ) ){
        // subtracts memory variables if attenuation
        compute_element_cm_att_stress(tx,working_element,
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
        compute_element_cm_gravity(tx,working_element,
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

        fac1 = d_hprimewgll_xx[I*NGLLX+l];
        offset = K*NGLL2+J*NGLLX+l;
        tempx1l += s_tempx1[offset]*fac1;
        tempy1l += s_tempy1[offset]*fac1;
        tempz1l += s_tempz1[offset]*fac1;

        fac2 = d_hprimewgll_yy[J*NGLLX+l];
        offset = K*NGLL2+l*NGLLX+I;
        tempx2l += s_tempx2[offset]*fac2;
        tempy2l += s_tempy2[offset]*fac2;
        tempz2l += s_tempz2[offset]*fac2;

        fac3 = d_hprimewgll_zz[K*NGLLX+l];
        offset = l*NGLL2+J*NGLLX+I;
        tempx3l += s_tempx3[offset]*fac3;
        tempy3l += s_tempy3[offset]*fac3;
        tempz3l += s_tempz3[offset]*fac3;

      }
#else

      tempx1l = s_tempx1[K*NGLL2+J*NGLLX]*d_hprimewgll_xx[I*NGLLX]
              + s_tempx1[K*NGLL2+J*NGLLX+1]*d_hprimewgll_xx[I*NGLLX+1]
              + s_tempx1[K*NGLL2+J*NGLLX+2]*d_hprimewgll_xx[I*NGLLX+2]
              + s_tempx1[K*NGLL2+J*NGLLX+3]*d_hprimewgll_xx[I*NGLLX+3]
              + s_tempx1[K*NGLL2+J*NGLLX+4]*d_hprimewgll_xx[I*NGLLX+4];

      tempy1l = s_tempy1[K*NGLL2+J*NGLLX]*d_hprimewgll_xx[I*NGLLX]
              + s_tempy1[K*NGLL2+J*NGLLX+1]*d_hprimewgll_xx[I*NGLLX+1]
              + s_tempy1[K*NGLL2+J*NGLLX+2]*d_hprimewgll_xx[I*NGLLX+2]
              + s_tempy1[K*NGLL2+J*NGLLX+3]*d_hprimewgll_xx[I*NGLLX+3]
              + s_tempy1[K*NGLL2+J*NGLLX+4]*d_hprimewgll_xx[I*NGLLX+4];

      tempz1l = s_tempz1[K*NGLL2+J*NGLLX]*d_hprimewgll_xx[I*NGLLX]
              + s_tempz1[K*NGLL2+J*NGLLX+1]*d_hprimewgll_xx[I*NGLLX+1]
              + s_tempz1[K*NGLL2+J*NGLLX+2]*d_hprimewgll_xx[I*NGLLX+2]
              + s_tempz1[K*NGLL2+J*NGLLX+3]*d_hprimewgll_xx[I*NGLLX+3]
              + s_tempz1[K*NGLL2+J*NGLLX+4]*d_hprimewgll_xx[I*NGLLX+4];

      tempx2l = s_tempx2[K*NGLL2+I]*d_hprimewgll_yy[J*NGLLX]
              + s_tempx2[K*NGLL2+NGLLX+I]*d_hprimewgll_yy[J*NGLLX+1]
              + s_tempx2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_yy[J*NGLLX+2]
              + s_tempx2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_yy[J*NGLLX+3]
              + s_tempx2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_yy[J*NGLLX+4];

      tempy2l = s_tempy2[K*NGLL2+I]*d_hprimewgll_yy[J*NGLLX]
              + s_tempy2[K*NGLL2+NGLLX+I]*d_hprimewgll_yy[J*NGLLX+1]
              + s_tempy2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_yy[J*NGLLX+2]
              + s_tempy2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_yy[J*NGLLX+3]
              + s_tempy2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_yy[J*NGLLX+4];

      tempz2l = s_tempz2[K*NGLL2+I]*d_hprimewgll_yy[J*NGLLX]
              + s_tempz2[K*NGLL2+NGLLX+I]*d_hprimewgll_yy[J*NGLLX+1]
              + s_tempz2[K*NGLL2+2*NGLLX+I]*d_hprimewgll_yy[J*NGLLX+2]
              + s_tempz2[K*NGLL2+3*NGLLX+I]*d_hprimewgll_yy[J*NGLLX+3]
              + s_tempz2[K*NGLL2+4*NGLLX+I]*d_hprimewgll_yy[J*NGLLX+4];

      tempx3l = s_tempx3[J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX]
              + s_tempx3[NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+1]
              + s_tempx3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+2]
              + s_tempx3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+3]
              + s_tempx3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+4];

      tempy3l = s_tempy3[J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX]
              + s_tempy3[NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+1]
              + s_tempy3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+2]
              + s_tempy3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+3]
              + s_tempy3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+4];

      tempz3l = s_tempz3[J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX]
              + s_tempz3[NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+1]
              + s_tempz3[2*NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+2]
              + s_tempz3[3*NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+3]
              + s_tempz3[4*NGLL2+J*NGLLX+I]*d_hprimewgll_zz[K*NGLLX+4];

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

#ifdef USE_TEXTURES
      d_accel[iglob] = tex1Dfetch(tex_accel, iglob) + sum_terms1 ;
      d_accel[iglob + NGLOB] = tex1Dfetch(tex_accel, iglob + NGLOB) + sum_terms2 ;
      d_accel[iglob + 2*NGLOB] = tex1Dfetch(tex_accel, iglob + 2*NGLOB) + sum_terms3 ;
#else
  /* OLD/To be implemented version that uses coloring to get around race condition. About 1.6x faster */


#ifdef USE_MESH_COLORING_GPU
      // no atomic operation needed, colors don't share global points between elements
      d_accel[iglob*3]     += sum_terms1;
      d_accel[iglob*3 + 1] += sum_terms2;
      d_accel[iglob*3 + 2] += sum_terms3;
#else
      //mesh coloring
      if( use_mesh_coloring_gpu ){

       // no atomic operation needed, colors don't share global points between elements
        d_accel[iglob*3]     += sum_terms1;
        d_accel[iglob*3 + 1] += sum_terms2;
        d_accel[iglob*3 + 2] += sum_terms3;

      }else{

        // for testing purposes only: w/out atomic updates
        //d_accel[iglob*3] -= (0.00000001f*tempx1l + 0.00000001f*tempx2l + 0.00000001f*tempx3l);
        //d_accel[iglob*3 + 1] -= (0.00000001f*tempy1l + 0.00000001f*tempy2l + 0.00000001f*tempy3l);
        //d_accel[iglob*3 + 2] -= (0.00000001f*tempz1l + 0.00000001f*tempz2l + 0.00000001f*tempz3l);

        atomicAdd(&d_accel[iglob*3], sum_terms1);
        atomicAdd(&d_accel[iglob*3+1], sum_terms2);
        atomicAdd(&d_accel[iglob*3+2], sum_terms3);

      }
#endif

#endif

      // update memory variables based upon the Runge-Kutta scheme
      if( ATTENUATION && ( ! USE_ATTENUATION_MIMIC ) ){
        compute_element_cm_att_memory(tx,working_element,
                                  d_muvstore,
                                  factor_common,alphaval,betaval,gammaval,
                                  R_xx,R_yy,R_xy,R_xz,R_yz,
                                  epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz,
                                  epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc,
                                  ANISOTROPY,d_c44store);
      }

      // save deviatoric strain for Runge-Kutta scheme
      if( COMPUTE_AND_STORE_STRAIN ){
        int ijk_ispec = tx + working_element*NGLL3;

        // fortran: epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
        epsilondev_xx[ijk_ispec] = epsilondev_xx_loc;
        epsilondev_yy[ijk_ispec] = epsilondev_yy_loc;
        epsilondev_xy[ijk_ispec] = epsilondev_xy_loc;
        epsilondev_xz[ijk_ispec] = epsilondev_xz_loc;
        epsilondev_yz[ijk_ispec] = epsilondev_yz_loc;
      }

    }

#else  // of #ifndef MAKE_KERNEL2_BECOME_STUPID_FOR_TESTS
    d_accel[iglob] -= 0.00000001f;
    d_accel[iglob + NGLOB] -= 0.00000001f;
    d_accel[iglob + 2*NGLOB] -= 0.00000001f;
#endif // of #ifndef MAKE_KERNEL2_BECOME_STUPID_FOR_TESTS

}
/* ----------------------------------------------------------------------------------------------- */

void Kernel_2_crust_mantle(int nb_blocks_to_compute,Mesh* mp,
                          realw d_deltat,
                          int d_iphase,
                          int* d_ibool,
                          int* d_ispec_is_tiso,
                          realw* d_xix,realw* d_xiy,realw* d_xiz,
                          realw* d_etax,realw* d_etay,realw* d_etaz,
                          realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                          realw* d_kappavstore,realw* d_muvstore,
                          realw* d_kappahstore,realw* d_muhstore,
                          realw* d_eta_anisostore,
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
                          realw* d_c11store,realw* d_c12store,realw* d_c13store,
                          realw* d_c14store,realw* d_c15store,realw* d_c16store,
                          realw* d_c22store,realw* d_c23store,realw* d_c24store,
                          realw* d_c25store,realw* d_c26store,realw* d_c33store,
                          realw* d_c34store,realw* d_c35store,realw* d_c36store,
                          realw* d_c44store,realw* d_c45store,realw* d_c46store,
                          realw* d_c55store,realw* d_c56store,realw* d_c66store
                          ){

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("before kernel Kernel_2_crust_mantle");
#endif

  /* if the grid can handle the number of blocks, we let it be 1D */
  /* grid_2_x = nb_elem_color; */
  /* nb_elem_color is just how many blocks we are computing now */

  int num_blocks_x = nb_blocks_to_compute;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  int blocksize = NGLL3_PADDED;
  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // Cuda timing
  // cudaEvent_t start, stop;
  // realw time;
  // cudaEventCreate(&start);
  // cudaEventCreate(&stop);
  // cudaEventRecord( start, 0 );

  Kernel_2_crust_mantle_impl<<<grid,threads>>>(nb_blocks_to_compute,
                                  mp->NGLOB_CRUST_MANTLE,
                                  d_ibool,
                                  d_ispec_is_tiso,
                                  mp->d_phase_ispec_inner_crust_mantle,
                                  mp->num_phase_ispec_crust_mantle,
                                  d_iphase,
                                  d_deltat,
                                  mp->use_mesh_coloring_gpu,
                                  mp->d_displ_crust_mantle,
                                  mp->d_veloc_crust_mantle,
                                  mp->d_accel_crust_mantle,
                                  d_xix, d_xiy, d_xiz,
                                  d_etax, d_etay, d_etaz,
                                  d_gammax, d_gammay, d_gammaz,
                                  d_kappavstore, d_muvstore,
                                  d_kappahstore, d_muhstore,
                                  d_eta_anisostore,
                                  mp->compute_and_store_strain,
                                  d_epsilondev_xx,d_epsilondev_yy,d_epsilondev_xy,
                                  d_epsilondev_xz,d_epsilondev_yz,
                                  d_epsilon_trace_over_3,
                                  mp->simulation_type,
                                  mp->attenuation,
                                  mp->attenuation_new,
                                  mp->use_attenuation_mimic,
                                  d_one_minus_sum_beta,d_factor_common,
                                  d_R_xx,d_R_yy,d_R_xy,d_R_xz,d_R_yz,
                                  mp->d_alphaval,mp->d_betaval,mp->d_gammaval,
                                  mp->anisotropic_3D_mantle,
                                  d_c11store,d_c12store,d_c13store,
                                  d_c14store,d_c15store,d_c16store,
                                  d_c22store,d_c23store,d_c24store,
                                  d_c25store,d_c26store,d_c33store,
                                  d_c34store,d_c35store,d_c36store,
                                  d_c44store,d_c45store,d_c46store,
                                  d_c55store,d_c56store,d_c66store,
                                  mp->gravity,
                                  mp->d_xstore_crust_mantle,mp->d_ystore_crust_mantle,mp->d_zstore_crust_mantle,
                                  mp->d_minus_gravity_table,
                                  mp->d_minus_deriv_gravity_table,
                                  mp->d_density_table,
                                  mp->d_wgll_cube,
                                  mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY);


  if(mp->simulation_type == 3) {
    Kernel_2_crust_mantle_impl<<< grid,threads>>>(nb_blocks_to_compute,
                                     mp->NGLOB_CRUST_MANTLE,
                                     d_ibool,
                                     d_ispec_is_tiso,
                                     mp->d_phase_ispec_inner_crust_mantle,
                                     mp->num_phase_ispec_crust_mantle,
                                     d_iphase,
                                     d_deltat,
                                     mp->use_mesh_coloring_gpu,
                                     mp->d_b_displ_crust_mantle,
                                     mp->d_b_veloc_crust_mantle,
                                     mp->d_b_accel_crust_mantle,
                                     d_xix, d_xiy, d_xiz,
                                     d_etax, d_etay, d_etaz,
                                     d_gammax, d_gammay, d_gammaz,
                                     d_kappavstore, d_muvstore,
                                     d_kappahstore, d_muhstore,
                                     d_eta_anisostore,
                                     mp->compute_and_store_strain,
                                     d_b_epsilondev_xx,d_b_epsilondev_yy,d_b_epsilondev_xy,
                                     d_b_epsilondev_xz,d_b_epsilondev_yz,
                                     d_b_epsilon_trace_over_3,
                                     mp->simulation_type,
                                     mp->attenuation,
                                     mp->attenuation_new,
                                     mp->use_attenuation_mimic,
                                     d_one_minus_sum_beta,d_factor_common,
                                     d_b_R_xx,d_b_R_yy,d_b_R_xy,d_b_R_xz,d_b_R_yz,
                                     mp->d_b_alphaval,mp->d_b_betaval,mp->d_b_gammaval,
                                     mp->anisotropic_3D_mantle,
                                     d_c11store,d_c12store,d_c13store,
                                     d_c14store,d_c15store,d_c16store,
                                     d_c22store,d_c23store,d_c24store,
                                     d_c25store,d_c26store,d_c33store,
                                     d_c34store,d_c35store,d_c36store,
                                     d_c44store,d_c45store,d_c46store,
                                     d_c55store,d_c56store,d_c66store,
                                     mp->gravity,
                                     mp->d_xstore_crust_mantle,mp->d_ystore_crust_mantle,mp->d_zstore_crust_mantle,
                                     mp->d_minus_gravity_table,
                                     mp->d_minus_deriv_gravity_table,
                                     mp->d_density_table,
                                     mp->d_wgll_cube,
                                     mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY);
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
  exit_on_cuda_error("Kernel_2_crust_mantle");
#endif
}

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_forces_crust_mantle_cuda,
              COMPUTE_FORCES_CRUST_MANTLE_CUDA)(long* Mesh_pointer_f,
                                                realw* deltat,
                                                int* iphase) {

  TRACE("compute_forces_crust_mantle_cuda");

//daniel: debug time
//  printf("Running compute_forces\n");
//  double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  int num_elements;

  if( *iphase == 1 )
    num_elements = mp->nspec_outer_crust_mantle;
  else
    num_elements = mp->nspec_inner_crust_mantle;

  // checks if anything to do
  if( num_elements == 0 ) return;

  // mesh coloring
  if( mp->use_mesh_coloring_gpu ){

    // note: array offsets require sorted arrays, such that e.g. ibool starts with elastic elements
    //         and followed by acoustic ones.
    //         elastic elements also start with outer than inner element ordering

    int nb_colors,nb_blocks_to_compute;
    int istart;
    int color_offset,color_offset_nonpadded,color_offset_nonpadded_att2;
    int color_offset_ispec;

    // sets up color loop
    if( *iphase == 1 ){
      // outer elements
      nb_colors = mp->num_colors_outer_crust_mantle;
      istart = 0;

      // array offsets
      color_offset = 0;
      color_offset_nonpadded = 0;
      color_offset_nonpadded_att2 = 0;
      color_offset_ispec = 0;
    }else{
      // inner elements (start after outer elements)
      nb_colors = mp->num_colors_outer_crust_mantle + mp->num_colors_inner_crust_mantle;
      istart = mp->num_colors_outer_crust_mantle;

      // array offsets
      color_offset = (mp->nspec_outer_crust_mantle) * NGLL3_PADDED;
      color_offset_nonpadded = (mp->nspec_outer_crust_mantle) * NGLL3;
      color_offset_nonpadded_att2 = (mp->nspec_outer_crust_mantle) * NGLL3 * N_SLS;
      color_offset_ispec = mp->nspec_outer_crust_mantle;
    }

    // loops over colors
    for(int icolor = istart; icolor < nb_colors; icolor++){

      nb_blocks_to_compute = mp->h_num_elem_colors_crust_mantle[icolor];

      // checks
      //if( nb_blocks_to_compute <= 0 ){
      //  printf("error number of elastic color blocks: %d -- color = %d \n",nb_blocks_to_compute,icolor);
      //  exit(EXIT_FAILURE);
      //}

      Kernel_2_crust_mantle(nb_blocks_to_compute,mp,
                            *deltat,
                            *iphase,
                            mp->d_ibool_crust_mantle + color_offset_nonpadded,
                            mp->d_ispec_is_tiso_crust_mantle + color_offset_ispec,
                            mp->d_xix_crust_mantle + color_offset,
                            mp->d_xiy_crust_mantle + color_offset,
                            mp->d_xiz_crust_mantle + color_offset,
                            mp->d_etax_crust_mantle + color_offset,
                            mp->d_etay_crust_mantle + color_offset,
                            mp->d_etaz_crust_mantle + color_offset,
                            mp->d_gammax_crust_mantle + color_offset,
                            mp->d_gammay_crust_mantle + color_offset,
                            mp->d_gammaz_crust_mantle + color_offset,
                            mp->d_kappavstore_crust_mantle + color_offset,
                            mp->d_muvstore_crust_mantle + color_offset,
                            mp->d_kappahstore_crust_mantle + color_offset,
                            mp->d_muhstore_crust_mantle + color_offset,
                            mp->d_eta_anisostore_crust_mantle + color_offset,
                            mp->d_epsilondev_xx_crust_mantle + color_offset_nonpadded,
                            mp->d_epsilondev_yy_crust_mantle + color_offset_nonpadded,
                            mp->d_epsilondev_xy_crust_mantle + color_offset_nonpadded,
                            mp->d_epsilondev_xz_crust_mantle + color_offset_nonpadded,
                            mp->d_epsilondev_yz_crust_mantle + color_offset_nonpadded,
                            mp->d_eps_trace_over_3_crust_mantle + color_offset_nonpadded,
                            mp->d_one_minus_sum_beta_crust_mantle + color_offset_nonpadded,
                            mp->d_factor_common_crust_mantle + color_offset_nonpadded_att2,
                            mp->d_R_xx_crust_mantle + color_offset_nonpadded,
                            mp->d_R_yy_crust_mantle + color_offset_nonpadded,
                            mp->d_R_xy_crust_mantle + color_offset_nonpadded,
                            mp->d_R_xz_crust_mantle + color_offset_nonpadded,
                            mp->d_R_yz_crust_mantle + color_offset_nonpadded,
                            mp->d_b_epsilondev_xx_crust_mantle + color_offset_nonpadded,
                            mp->d_b_epsilondev_yy_crust_mantle + color_offset_nonpadded,
                            mp->d_b_epsilondev_xy_crust_mantle + color_offset_nonpadded,
                            mp->d_b_epsilondev_xz_crust_mantle + color_offset_nonpadded,
                            mp->d_b_epsilondev_yz_crust_mantle + color_offset_nonpadded,
                            mp->d_b_eps_trace_over_3_crust_mantle + color_offset_nonpadded,
                            mp->d_b_R_xx_crust_mantle + color_offset_nonpadded,
                            mp->d_b_R_yy_crust_mantle + color_offset_nonpadded,
                            mp->d_b_R_xy_crust_mantle + color_offset_nonpadded,
                            mp->d_b_R_xz_crust_mantle + color_offset_nonpadded,
                            mp->d_b_R_yz_crust_mantle + color_offset_nonpadded,
                            mp->d_c11store_crust_mantle + color_offset,
                            mp->d_c12store_crust_mantle + color_offset,
                            mp->d_c13store_crust_mantle + color_offset,
                            mp->d_c14store_crust_mantle + color_offset,
                            mp->d_c15store_crust_mantle + color_offset,
                            mp->d_c16store_crust_mantle + color_offset,
                            mp->d_c22store_crust_mantle + color_offset,
                            mp->d_c23store_crust_mantle + color_offset,
                            mp->d_c24store_crust_mantle + color_offset,
                            mp->d_c25store_crust_mantle + color_offset,
                            mp->d_c26store_crust_mantle + color_offset,
                            mp->d_c33store_crust_mantle + color_offset,
                            mp->d_c34store_crust_mantle + color_offset,
                            mp->d_c35store_crust_mantle + color_offset,
                            mp->d_c36store_crust_mantle + color_offset,
                            mp->d_c44store_crust_mantle + color_offset,
                            mp->d_c45store_crust_mantle + color_offset,
                            mp->d_c46store_crust_mantle + color_offset,
                            mp->d_c55store_crust_mantle + color_offset,
                            mp->d_c56store_crust_mantle + color_offset,
                            mp->d_c66store_crust_mantle + color_offset
                            );

      // for padded and aligned arrays
      color_offset += nb_blocks_to_compute * NGLL3_PADDED;
      // for no-aligned arrays
      color_offset_nonpadded += nb_blocks_to_compute * NGLL3;
      // for factor_common array
      color_offset_nonpadded_att2 += nb_blocks_to_compute * NGLL3 * N_SLS;
      // for array(ispec)
      color_offset_ispec += nb_blocks_to_compute;
    }

  }else{

    // no mesh coloring: uses atomic updates

    Kernel_2_crust_mantle(num_elements,mp,
                          *deltat,
                          *iphase,
                          mp->d_ibool_crust_mantle,
                          mp->d_ispec_is_tiso_crust_mantle,
                          mp->d_xix_crust_mantle,mp->d_xiy_crust_mantle,mp->d_xiz_crust_mantle,
                          mp->d_etax_crust_mantle,mp->d_etay_crust_mantle,mp->d_etaz_crust_mantle,
                          mp->d_gammax_crust_mantle,mp->d_gammay_crust_mantle,mp->d_gammaz_crust_mantle,
                          mp->d_kappavstore_crust_mantle,mp->d_muvstore_crust_mantle,
                          mp->d_kappahstore_crust_mantle,mp->d_muhstore_crust_mantle,
                          mp->d_eta_anisostore_crust_mantle,
                          mp->d_epsilondev_xx_crust_mantle,
                          mp->d_epsilondev_yy_crust_mantle,
                          mp->d_epsilondev_xy_crust_mantle,
                          mp->d_epsilondev_xz_crust_mantle,
                          mp->d_epsilondev_yz_crust_mantle,
                          mp->d_eps_trace_over_3_crust_mantle,
                          mp->d_one_minus_sum_beta_crust_mantle,
                          mp->d_factor_common_crust_mantle,
                          mp->d_R_xx_crust_mantle,
                          mp->d_R_yy_crust_mantle,
                          mp->d_R_xy_crust_mantle,
                          mp->d_R_xz_crust_mantle,
                          mp->d_R_yz_crust_mantle,
                          mp->d_b_epsilondev_xx_crust_mantle,
                          mp->d_b_epsilondev_yy_crust_mantle,
                          mp->d_b_epsilondev_xy_crust_mantle,
                          mp->d_b_epsilondev_xz_crust_mantle,
                          mp->d_b_epsilondev_yz_crust_mantle,
                          mp->d_b_eps_trace_over_3_crust_mantle,
                          mp->d_b_R_xx_crust_mantle,
                          mp->d_b_R_yy_crust_mantle,
                          mp->d_b_R_xy_crust_mantle,
                          mp->d_b_R_xz_crust_mantle,
                          mp->d_b_R_yz_crust_mantle,
                          mp->d_c11store_crust_mantle,mp->d_c12store_crust_mantle,mp->d_c13store_crust_mantle,
                          mp->d_c14store_crust_mantle,mp->d_c15store_crust_mantle,mp->d_c16store_crust_mantle,
                          mp->d_c22store_crust_mantle,mp->d_c23store_crust_mantle,mp->d_c24store_crust_mantle,
                          mp->d_c25store_crust_mantle,mp->d_c26store_crust_mantle,mp->d_c33store_crust_mantle,
                          mp->d_c34store_crust_mantle,mp->d_c35store_crust_mantle,mp->d_c36store_crust_mantle,
                          mp->d_c44store_crust_mantle,mp->d_c45store_crust_mantle,mp->d_c46store_crust_mantle,
                          mp->d_c55store_crust_mantle,mp->d_c56store_crust_mantle,mp->d_c66store_crust_mantle
                          );
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("compute_forces_crust_mantle_cuda");
#endif
}
