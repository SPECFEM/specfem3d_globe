//note: please do not modify this file manually!
//      this file has been generated automatically by BOAST version 2.1.0
//      by: make boast_kernels

/*
!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

const char * crust_mantle_impl_kernel_forward_program = "\
inline void atomicAdd(volatile __global float *source, const float val) {\n\
  union {\n\
    unsigned int iVal;\n\
    float fVal;\n\
  } res, orig;\n\
  do {\n\
    orig.fVal = *source;\n\
    res.fVal = orig.fVal + val;\n\
  } while (atomic_cmpxchg((volatile __global unsigned int *)source, orig.iVal, res.iVal) != orig.iVal);\n\
}\n\
#ifndef INDEX2\n\
#define INDEX2(isize,i,j) i + isize*j\n\
#endif\n\
#ifndef INDEX3\n\
#define INDEX3(isize,jsize,i,j,k) i + isize*(j + jsize*k)\n\
#endif\n\
#ifndef INDEX4\n\
#define INDEX4(isize,jsize,ksize,i,j,k,x) i + isize*(j + jsize*(k + ksize*x))\n\
#endif\n\
#ifndef INDEX5\n\
#define INDEX5(isize,jsize,ksize,xsize,i,j,k,x,y) i + isize*(j + jsize*(k + ksize*(x + xsize*y)))\n\
#endif\n\
\n\
#ifndef NDIM\n\
#define NDIM 3\n\
#endif\n\
#ifndef NGLLX\n\
#define NGLLX 5\n\
#endif\n\
#ifndef NGLL2\n\
#define NGLL2 25\n\
#endif\n\
#ifndef NGLL3\n\
#define NGLL3 125\n\
#endif\n\
#ifndef NGLL3_PADDED\n\
#define NGLL3_PADDED 128\n\
#endif\n\
#ifndef N_SLS\n\
#define N_SLS 3\n\
#endif\n\
#ifndef IREGION_CRUST_MANTLE\n\
#define IREGION_CRUST_MANTLE 1\n\
#endif\n\
#ifndef IREGION_INNER_CORE\n\
#define IREGION_INNER_CORE 3\n\
#endif\n\
#ifndef IFLAG_IN_FICTITIOUS_CUBE\n\
#define IFLAG_IN_FICTITIOUS_CUBE 11\n\
#endif\n\
#ifndef COLORING_MIN_NSPEC_INNER_CORE\n\
#define COLORING_MIN_NSPEC_INNER_CORE 1000\n\
#endif\n\
#ifndef COLORING_MIN_NSPEC_OUTER_CORE\n\
#define COLORING_MIN_NSPEC_OUTER_CORE 1000\n\
#endif\n\
#ifndef BLOCKSIZE_TRANSFER\n\
#define BLOCKSIZE_TRANSFER 256\n\
#endif\n\
\n\
\n\
\n\
#if __OPENCL_C_VERSION__ && __OPENCL_C_VERSION__ >= 120\n\
static\n\
#endif\n\
void compute_element_cm_att_stress(const int tx, const int working_element, const __global float * R_xx, const __global float * R_yy, const __global float * R_xy, const __global float * R_xz, const __global float * R_yz, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_xz, float * sigma_yz){\n\
  int offset;\n\
  float R_xx_val;\n\
  float R_yy_val;\n\
\n\
  for (int i_sls = 0; i_sls < N_SLS; i_sls += 1) {\n\
\n\
    offset = tx + (NGLL3) * (i_sls + (N_SLS) * (working_element));\n\
\n\
    R_xx_val = R_xx[offset];\n\
    R_yy_val = R_yy[offset];\n\
\n\
    sigma_xx[0] = sigma_xx[0] - (R_xx_val);\n\
    sigma_yy[0] = sigma_yy[0] - (R_yy_val);\n\
    sigma_zz[0] = sigma_zz[0] + R_xx_val + R_yy_val;\n\
    sigma_xy[0] = sigma_xy[0] - (R_xy[offset]);\n\
    sigma_xz[0] = sigma_xz[0] - (R_xz[offset]);\n\
    sigma_yz[0] = sigma_yz[0] - (R_yz[offset]);\n\
  }\n\
}\n\
\n\
#if __OPENCL_C_VERSION__ && __OPENCL_C_VERSION__ >= 120\n\
static\n\
#endif\n\
void compute_element_cm_att_memory(const int tx, const int working_element, const __global float * d_muvstore, const __global float * factor_common, const __global float * alphaval, const __global float * betaval, const __global float * gammaval, __global float * R_xx, __global float * R_yy, __global float * R_xy, __global float * R_xz, __global float * R_yz, const __global float * epsilondev_xx, const __global float * epsilondev_yy, const __global float * epsilondev_xy, const __global float * epsilondev_xz, const __global float * epsilondev_yz, const float epsilondev_xx_loc, const float epsilondev_yy_loc, const float epsilondev_xy_loc, const float epsilondev_xz_loc, const float epsilondev_yz_loc, const int USE_3D_ATTENUATION_ARRAYS){\n\
  int offset;\n\
  float mul;\n\
  float factor_loc;\n\
  float sn;\n\
  float snp1;\n\
  float alphaval_loc;\n\
  float betaval_loc;\n\
  float gammaval_loc;\n\
\n\
  mul = d_muvstore[tx + (NGLL3_PADDED) * (working_element)];\n\
  const int offset_eps = tx + (NGLL3) * (working_element);\n\
\n\
  for (int i_sls = 0; i_sls < N_SLS; i_sls += 1) {\n\
    offset = tx + (NGLL3) * (i_sls + (N_SLS) * (working_element));\n\
    if (USE_3D_ATTENUATION_ARRAYS) {\n\
      factor_loc = (mul) * (factor_common[offset]);\n\
    } else {\n\
      factor_loc = (mul) * (factor_common[i_sls + (N_SLS) * (working_element)]);\n\
    }\n\
\n\
    alphaval_loc = alphaval[i_sls];\n\
    betaval_loc = betaval[i_sls];\n\
    gammaval_loc = gammaval[i_sls];\n\
\n\
    sn = (factor_loc) * (epsilondev_xx[offset_eps]);\n\
    snp1 = (factor_loc) * (epsilondev_xx_loc);\n\
    R_xx[offset] = (alphaval_loc) * (R_xx[offset]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);\n\
    sn = (factor_loc) * (epsilondev_yy[offset_eps]);\n\
    snp1 = (factor_loc) * (epsilondev_yy_loc);\n\
    R_yy[offset] = (alphaval_loc) * (R_yy[offset]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);\n\
    sn = (factor_loc) * (epsilondev_xy[offset_eps]);\n\
    snp1 = (factor_loc) * (epsilondev_xy_loc);\n\
    R_xy[offset] = (alphaval_loc) * (R_xy[offset]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);\n\
    sn = (factor_loc) * (epsilondev_xz[offset_eps]);\n\
    snp1 = (factor_loc) * (epsilondev_xz_loc);\n\
    R_xz[offset] = (alphaval_loc) * (R_xz[offset]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);\n\
    sn = (factor_loc) * (epsilondev_yz[offset_eps]);\n\
    snp1 = (factor_loc) * (epsilondev_yz_loc);\n\
    R_yz[offset] = (alphaval_loc) * (R_yz[offset]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);\n\
  }\n\
}\n\
\n\
#if __OPENCL_C_VERSION__ && __OPENCL_C_VERSION__ >= 120\n\
static\n\
#endif\n\
void compute_element_cm_att_memory_lddrk(const int tx, const int working_element, const __global float * d_muvstore, const __global float * factor_common, const __global float * tau_sigmainvval, __global float * R_xx, __global float * R_yy, __global float * R_xy, __global float * R_xz, __global float * R_yz, __global float * R_xx_lddrk, __global float * R_yy_lddrk, __global float * R_xy_lddrk, __global float * R_xz_lddrk, __global float * R_yz_lddrk, const float alpha_lddrk, const float beta_lddrk, const float deltat, const float epsilondev_xx_loc, const float epsilondev_yy_loc, const float epsilondev_xy_loc, const float epsilondev_xz_loc, const float epsilondev_yz_loc, const int USE_3D_ATTENUATION_ARRAYS){\n\
  int offset;\n\
  float mul;\n\
  float factor_loc;\n\
  float sn;\n\
  float snp1;\n\
  float tau_sigmainv_loc;\n\
\n\
  mul = d_muvstore[tx + (NGLL3_PADDED) * (working_element)];\n\
\n\
  for (int i_sls = 0; i_sls < N_SLS; i_sls += 1) {\n\
    offset = tx + (NGLL3) * (i_sls + (N_SLS) * (working_element));\n\
    if (USE_3D_ATTENUATION_ARRAYS) {\n\
      factor_loc = (mul) * (factor_common[offset]);\n\
    } else {\n\
      factor_loc = (mul) * (factor_common[i_sls + (N_SLS) * (working_element)]);\n\
    }\n\
\n\
    tau_sigmainv_loc = tau_sigmainvval[i_sls];\n\
\n\
    sn = (tau_sigmainv_loc) * (R_xx[offset]);\n\
    snp1 = (factor_loc) * (epsilondev_xx_loc);\n\
    R_xx_lddrk[offset] = (alpha_lddrk) * (R_xx_lddrk[offset]) + (deltat) * (snp1 - (sn));\n\
    R_xx[offset] = R_xx[offset] + (beta_lddrk) * (R_xx_lddrk[offset]);\n\
    sn = (tau_sigmainv_loc) * (R_yy[offset]);\n\
    snp1 = (factor_loc) * (epsilondev_yy_loc);\n\
    R_yy_lddrk[offset] = (alpha_lddrk) * (R_yy_lddrk[offset]) + (deltat) * (snp1 - (sn));\n\
    R_yy[offset] = R_yy[offset] + (beta_lddrk) * (R_yy_lddrk[offset]);\n\
    sn = (tau_sigmainv_loc) * (R_xy[offset]);\n\
    snp1 = (factor_loc) * (epsilondev_xy_loc);\n\
    R_xy_lddrk[offset] = (alpha_lddrk) * (R_xy_lddrk[offset]) + (deltat) * (snp1 - (sn));\n\
    R_xy[offset] = R_xy[offset] + (beta_lddrk) * (R_xy_lddrk[offset]);\n\
    sn = (tau_sigmainv_loc) * (R_xz[offset]);\n\
    snp1 = (factor_loc) * (epsilondev_xz_loc);\n\
    R_xz_lddrk[offset] = (alpha_lddrk) * (R_xz_lddrk[offset]) + (deltat) * (snp1 - (sn));\n\
    R_xz[offset] = R_xz[offset] + (beta_lddrk) * (R_xz_lddrk[offset]);\n\
    sn = (tau_sigmainv_loc) * (R_yz[offset]);\n\
    snp1 = (factor_loc) * (epsilondev_yz_loc);\n\
    R_yz_lddrk[offset] = (alpha_lddrk) * (R_yz_lddrk[offset]) + (deltat) * (snp1 - (sn));\n\
    R_yz[offset] = R_yz[offset] + (beta_lddrk) * (R_yz_lddrk[offset]);\n\
  }\n\
}\n\
\n\
#if __OPENCL_C_VERSION__ && __OPENCL_C_VERSION__ >= 120\n\
static\n\
#endif\n\
void compute_element_cm_gravity(const int tx, const int iglob, const __global float * restrict gravity_pre_store, const __global float * restrict gravity_H, const __global float * restrict wgll_cube, const float jacobianl, const __local float * s_dummyx_loc, const __local float * s_dummyy_loc, const __local float * s_dummyz_loc, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_yx, float * sigma_xz, float * sigma_zx, float * sigma_yz, float * sigma_zy, float * rho_s_H1, float * rho_s_H2, float * rho_s_H3){\n\
  float gxl;\n\
  float gyl;\n\
  float gzl;\n\
  float Hxxl;\n\
  float Hyyl;\n\
  float Hzzl;\n\
  float Hxyl;\n\
  float Hxzl;\n\
  float Hyzl;\n\
  float sx_l;\n\
  float sy_l;\n\
  float sz_l;\n\
  float factor;\n\
\n\
  gxl = gravity_pre_store[0 + (3) * (iglob)];\n\
  gyl = gravity_pre_store[1 + (3) * (iglob)];\n\
  gzl = gravity_pre_store[2 + (3) * (iglob)];\n\
\n\
  Hxxl = gravity_H[0 + (6) * (iglob)];\n\
  Hyyl = gravity_H[1 + (6) * (iglob)];\n\
  Hzzl = gravity_H[2 + (6) * (iglob)];\n\
  Hxyl = gravity_H[3 + (6) * (iglob)];\n\
  Hxzl = gravity_H[4 + (6) * (iglob)];\n\
  Hyzl = gravity_H[5 + (6) * (iglob)];\n\
\n\
  sx_l = s_dummyx_loc[tx];\n\
  sy_l = s_dummyy_loc[tx];\n\
  sz_l = s_dummyz_loc[tx];\n\
\n\
  *(sigma_xx) = *(sigma_xx) + (sy_l) * (gyl) + (sz_l) * (gzl);\n\
  *(sigma_yy) = *(sigma_yy) + (sx_l) * (gxl) + (sz_l) * (gzl);\n\
  *(sigma_zz) = *(sigma_zz) + (sx_l) * (gxl) + (sy_l) * (gyl);\n\
  *(sigma_xy) = *(sigma_xy) - ((sx_l) * (gyl));\n\
  *(sigma_yx) = *(sigma_yx) - ((sy_l) * (gxl));\n\
  *(sigma_xz) = *(sigma_xz) - ((sx_l) * (gzl));\n\
  *(sigma_zx) = *(sigma_zx) - ((sz_l) * (gxl));\n\
  *(sigma_yz) = *(sigma_yz) - ((sy_l) * (gzl));\n\
  *(sigma_zy) = *(sigma_zy) - ((sz_l) * (gyl));\n\
\n\
  factor = (jacobianl) * (wgll_cube[tx]);\n\
  *(rho_s_H1) = (factor) * ((sx_l) * (Hxxl) + (sy_l) * (Hxyl) + (sz_l) * (Hxzl));\n\
  *(rho_s_H2) = (factor) * ((sx_l) * (Hxyl) + (sy_l) * (Hyyl) + (sz_l) * (Hyzl));\n\
  *(rho_s_H3) = (factor) * ((sx_l) * (Hxzl) + (sy_l) * (Hyzl) + (sz_l) * (Hzzl));\n\
}\n\
\n\
#if __OPENCL_C_VERSION__ && __OPENCL_C_VERSION__ >= 120\n\
static\n\
#endif\n\
void compute_element_cm_aniso(const int offset, const __global float * d_c11store, const __global float * d_c12store, const __global float * d_c13store, const __global float * d_c14store, const __global float * d_c15store, const __global float * d_c16store, const __global float * d_c22store, const __global float * d_c23store, const __global float * d_c24store, const __global float * d_c25store, const __global float * d_c26store, const __global float * d_c33store, const __global float * d_c34store, const __global float * d_c35store, const __global float * d_c36store, const __global float * d_c44store, const __global float * d_c45store, const __global float * d_c46store, const __global float * d_c55store, const __global float * d_c56store, const __global float * d_c66store, const float duxdxl, const float duxdyl, const float duxdzl, const float duydxl, const float duydyl, const float duydzl, const float duzdxl, const float duzdyl, const float duzdzl, const float duxdyl_plus_duydxl, const float duzdxl_plus_duxdzl, const float duzdyl_plus_duydzl, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_xz, float * sigma_yz){\n\
  float c11;\n\
  float c12;\n\
  float c13;\n\
  float c14;\n\
  float c15;\n\
  float c16;\n\
  float c22;\n\
  float c23;\n\
  float c24;\n\
  float c25;\n\
  float c26;\n\
  float c33;\n\
  float c34;\n\
  float c35;\n\
  float c36;\n\
  float c44;\n\
  float c45;\n\
  float c46;\n\
  float c55;\n\
  float c56;\n\
  float c66;\n\
\n\
  c11 = d_c11store[offset];\n\
  c12 = d_c12store[offset];\n\
  c13 = d_c13store[offset];\n\
  c14 = d_c14store[offset];\n\
  c15 = d_c15store[offset];\n\
  c16 = d_c16store[offset];\n\
  c22 = d_c22store[offset];\n\
  c23 = d_c23store[offset];\n\
  c24 = d_c24store[offset];\n\
  c25 = d_c25store[offset];\n\
  c26 = d_c26store[offset];\n\
  c33 = d_c33store[offset];\n\
  c34 = d_c34store[offset];\n\
  c35 = d_c35store[offset];\n\
  c36 = d_c36store[offset];\n\
  c44 = d_c44store[offset];\n\
  c45 = d_c45store[offset];\n\
  c46 = d_c46store[offset];\n\
  c55 = d_c55store[offset];\n\
  c56 = d_c56store[offset];\n\
  c66 = d_c66store[offset];\n\
\n\
  *(sigma_xx) = (c11) * (duxdxl) + (c16) * (duxdyl_plus_duydxl) + (c12) * (duydyl) + (c15) * (duzdxl_plus_duxdzl) + (c14) * (duzdyl_plus_duydzl) + (c13) * (duzdzl);\n\
  *(sigma_yy) = (c12) * (duxdxl) + (c26) * (duxdyl_plus_duydxl) + (c22) * (duydyl) + (c25) * (duzdxl_plus_duxdzl) + (c24) * (duzdyl_plus_duydzl) + (c23) * (duzdzl);\n\
  *(sigma_zz) = (c13) * (duxdxl) + (c36) * (duxdyl_plus_duydxl) + (c23) * (duydyl) + (c35) * (duzdxl_plus_duxdzl) + (c34) * (duzdyl_plus_duydzl) + (c33) * (duzdzl);\n\
  *(sigma_xy) = (c16) * (duxdxl) + (c66) * (duxdyl_plus_duydxl) + (c26) * (duydyl) + (c56) * (duzdxl_plus_duxdzl) + (c46) * (duzdyl_plus_duydzl) + (c36) * (duzdzl);\n\
  *(sigma_xz) = (c15) * (duxdxl) + (c56) * (duxdyl_plus_duydxl) + (c25) * (duydyl) + (c55) * (duzdxl_plus_duxdzl) + (c45) * (duzdyl_plus_duydzl) + (c35) * (duzdzl);\n\
  *(sigma_yz) = (c14) * (duxdxl) + (c46) * (duxdyl_plus_duydxl) + (c24) * (duydyl) + (c45) * (duzdxl_plus_duxdzl) + (c44) * (duzdyl_plus_duydzl) + (c34) * (duzdzl);\n\
}\n\
\n\
#if __OPENCL_C_VERSION__ && __OPENCL_C_VERSION__ >= 120\n\
static\n\
#endif\n\
void compute_element_cm_iso(const int offset, const __global float * d_kappavstore, const __global float * d_muvstore, const float duxdxl, const float duydyl, const float duzdzl, const float duxdxl_plus_duydyl, const float duxdxl_plus_duzdzl, const float duydyl_plus_duzdzl, const float duxdyl_plus_duydxl, const float duzdxl_plus_duxdzl, const float duzdyl_plus_duydzl, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_xz, float * sigma_yz){\n\
  float lambdal;\n\
  float mul;\n\
  float lambdalplus2mul;\n\
  float kappal;\n\
\n\
  mul = d_muvstore[offset];\n\
  kappal = d_kappavstore[offset];\n\
\n\
  lambdalplus2mul = kappal + (mul) * (1.3333333333333333f);\n\
  lambdal = lambdalplus2mul - ((mul) * (2.0f));\n\
\n\
  *(sigma_xx) = (lambdalplus2mul) * (duxdxl) + (lambdal) * (duydyl_plus_duzdzl);\n\
  *(sigma_yy) = (lambdalplus2mul) * (duydyl) + (lambdal) * (duxdxl_plus_duzdzl);\n\
  *(sigma_zz) = (lambdalplus2mul) * (duzdzl) + (lambdal) * (duxdxl_plus_duydyl);\n\
  *(sigma_xy) = (mul) * (duxdyl_plus_duydxl);\n\
  *(sigma_xz) = (mul) * (duzdxl_plus_duxdzl);\n\
  *(sigma_yz) = (mul) * (duzdyl_plus_duydzl);\n\
}\n\
\n\
#if __OPENCL_C_VERSION__ && __OPENCL_C_VERSION__ >= 120\n\
static\n\
#endif\n\
void compute_element_cm_tiso(const int offset, const __global float * d_c11store, const __global float * d_c12store, const __global float * d_c13store, const __global float * d_c14store, const __global float * d_c15store, const __global float * d_c16store, const __global float * d_c22store, const __global float * d_c23store, const __global float * d_c24store, const __global float * d_c25store, const __global float * d_c26store, const __global float * d_c33store, const __global float * d_c34store, const __global float * d_c35store, const __global float * d_c36store, const __global float * d_c44store, const __global float * d_c45store, const __global float * d_c46store, const __global float * d_c55store, const __global float * d_c56store, const __global float * d_c66store, const float duxdxl, const float duxdyl, const float duxdzl, const float duydxl, const float duydyl, const float duydzl, const float duzdxl, const float duzdyl, const float duzdzl, const float duxdyl_plus_duydxl, const float duzdxl_plus_duxdzl, const float duzdyl_plus_duydzl, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_xz, float * sigma_yz){\n\
  float c11;\n\
  float c12;\n\
  float c13;\n\
  float c14;\n\
  float c15;\n\
  float c16;\n\
  float c22;\n\
  float c23;\n\
  float c24;\n\
  float c25;\n\
  float c26;\n\
  float c33;\n\
  float c34;\n\
  float c35;\n\
  float c36;\n\
  float c44;\n\
  float c45;\n\
  float c46;\n\
  float c55;\n\
  float c56;\n\
  float c66;\n\
\n\
  c11 = d_c11store[offset];\n\
  c12 = d_c12store[offset];\n\
  c13 = d_c13store[offset];\n\
  c14 = d_c14store[offset];\n\
  c15 = d_c15store[offset];\n\
  c16 = d_c16store[offset];\n\
  c22 = d_c22store[offset];\n\
  c23 = d_c23store[offset];\n\
  c24 = d_c24store[offset];\n\
  c25 = d_c25store[offset];\n\
  c26 = d_c26store[offset];\n\
  c33 = d_c33store[offset];\n\
  c34 = d_c34store[offset];\n\
  c35 = d_c35store[offset];\n\
  c36 = d_c36store[offset];\n\
  c44 = d_c44store[offset];\n\
  c45 = d_c45store[offset];\n\
  c46 = d_c46store[offset];\n\
  c55 = d_c55store[offset];\n\
  c56 = d_c56store[offset];\n\
  c66 = d_c66store[offset];\n\
\n\
  *(sigma_xx) = (c11) * (duxdxl) + (c16) * (duxdyl_plus_duydxl) + (c12) * (duydyl) + (c15) * (duzdxl_plus_duxdzl) + (c14) * (duzdyl_plus_duydzl) + (c13) * (duzdzl);\n\
  *(sigma_yy) = (c12) * (duxdxl) + (c26) * (duxdyl_plus_duydxl) + (c22) * (duydyl) + (c25) * (duzdxl_plus_duxdzl) + (c24) * (duzdyl_plus_duydzl) + (c23) * (duzdzl);\n\
  *(sigma_zz) = (c13) * (duxdxl) + (c36) * (duxdyl_plus_duydxl) + (c23) * (duydyl) + (c35) * (duzdxl_plus_duxdzl) + (c34) * (duzdyl_plus_duydzl) + (c33) * (duzdzl);\n\
  *(sigma_xy) = (c16) * (duxdxl) + (c66) * (duxdyl_plus_duydxl) + (c26) * (duydyl) + (c56) * (duzdxl_plus_duxdzl) + (c46) * (duzdyl_plus_duydzl) + (c36) * (duzdzl);\n\
  *(sigma_xz) = (c15) * (duxdxl) + (c56) * (duxdyl_plus_duydxl) + (c25) * (duydyl) + (c55) * (duzdxl_plus_duxdzl) + (c45) * (duzdyl_plus_duydzl) + (c35) * (duzdzl);\n\
  *(sigma_yz) = (c14) * (duxdxl) + (c46) * (duxdyl_plus_duydxl) + (c24) * (duydyl) + (c45) * (duzdxl_plus_duxdzl) + (c44) * (duzdyl_plus_duydzl) + (c34) * (duzdzl);\n\
}\n\
\n\
\n\
/*----------------------------------------------*/\n\
// main function\n\
/*----------------------------------------------*/\n\
\n\
#ifdef USE_TEXTURES_FIELDS\n\
#ifdef USE_TEXTURES_CONSTANTS\n\
__kernel  void crust_mantle_impl_kernel_forward(const int nb_blocks_to_compute, const __global int * d_ibool, const __global int * d_ispec_is_tiso, const __global int * d_phase_ispec_inner, const int num_phase_ispec, const int d_iphase, const float deltat, const int use_mesh_coloring_gpu, const __global float * restrict d_displ, __global float * d_accel, const __global float * restrict d_xix, const __global float * restrict d_xiy, const __global float * restrict d_xiz, const __global float * restrict d_etax, const __global float * restrict d_etay, const __global float * restrict d_etaz, const __global float * restrict d_gammax, const __global float * restrict d_gammay, const __global float * restrict d_gammaz, const __global float * restrict d_hprime_xx, const __global float * restrict d_hprimewgll_xx, const __global float * restrict d_wgllwgll_xy, const __global float * restrict d_wgllwgll_xz, const __global float * restrict d_wgllwgll_yz, const __global float * restrict d_kappavstore, const __global float * restrict d_muvstore, const __global float * restrict d_kappahstore, const __global float * restrict d_muhstore, const __global float * restrict d_eta_anisostore, const int COMPUTE_AND_STORE_STRAIN, __global float * epsilondev_xx, __global float * epsilondev_yy, __global float * epsilondev_xy, __global float * epsilondev_xz, __global float * epsilondev_yz, __global float * epsilon_trace_over_3, const int ATTENUATION, const int PARTIAL_PHYS_DISPERSION_ONLY, const int USE_3D_ATTENUATION_ARRAYS, const __global float * restrict one_minus_sum_beta, const __global float * restrict factor_common, __global float * R_xx, __global float * R_yy, __global float * R_xy, __global float * R_xz, __global float * R_yz, __global float * R_xx_lddrk, __global float * R_yy_lddrk, __global float * R_xy_lddrk, __global float * R_xz_lddrk, __global float * R_yz_lddrk, const float alpha_lddrk, const float beta_lddrk, const int use_lddrk, const __global float * restrict alphaval, const __global float * restrict betaval, const __global float * restrict gammaval, const __global float * restrict tau_sigmainvval, const __global float * restrict d_c11store, const __global float * restrict d_c12store, const __global float * restrict d_c13store, const __global float * restrict d_c14store, const __global float * restrict d_c15store, const __global float * restrict d_c16store, const __global float * restrict d_c22store, const __global float * restrict d_c23store, const __global float * restrict d_c24store, const __global float * restrict d_c25store, const __global float * restrict d_c26store, const __global float * restrict d_c33store, const __global float * restrict d_c34store, const __global float * restrict d_c35store, const __global float * restrict d_c36store, const __global float * restrict d_c44store, const __global float * restrict d_c45store, const __global float * restrict d_c46store, const __global float * restrict d_c55store, const __global float * restrict d_c56store, const __global float * restrict d_c66store, const __global float * restrict d_rstore, const int GRAVITY, const __global float * restrict d_gravity_pre_store, const __global float * restrict d_gravity_H, const __global float * restrict wgll_cube, const int NSPEC_CRUST_MANTLE_STRAIN_ONLY, __read_only image2d_t d_displ_cm_tex, __read_only image2d_t d_accel_cm_tex, __read_only image2d_t d_hprime_xx_tex, __read_only image2d_t d_hprimewgll_xx_tex){\n\
#else\n\
__kernel  void crust_mantle_impl_kernel_forward(const int nb_blocks_to_compute, const __global int * d_ibool, const __global int * d_ispec_is_tiso, const __global int * d_phase_ispec_inner, const int num_phase_ispec, const int d_iphase, const float deltat, const int use_mesh_coloring_gpu, const __global float * restrict d_displ, __global float * d_accel, const __global float * restrict d_xix, const __global float * restrict d_xiy, const __global float * restrict d_xiz, const __global float * restrict d_etax, const __global float * restrict d_etay, const __global float * restrict d_etaz, const __global float * restrict d_gammax, const __global float * restrict d_gammay, const __global float * restrict d_gammaz, const __global float * restrict d_hprime_xx, const __global float * restrict d_hprimewgll_xx, const __global float * restrict d_wgllwgll_xy, const __global float * restrict d_wgllwgll_xz, const __global float * restrict d_wgllwgll_yz, const __global float * restrict d_kappavstore, const __global float * restrict d_muvstore, const __global float * restrict d_kappahstore, const __global float * restrict d_muhstore, const __global float * restrict d_eta_anisostore, const int COMPUTE_AND_STORE_STRAIN, __global float * epsilondev_xx, __global float * epsilondev_yy, __global float * epsilondev_xy, __global float * epsilondev_xz, __global float * epsilondev_yz, __global float * epsilon_trace_over_3, const int ATTENUATION, const int PARTIAL_PHYS_DISPERSION_ONLY, const int USE_3D_ATTENUATION_ARRAYS, const __global float * restrict one_minus_sum_beta, const __global float * restrict factor_common, __global float * R_xx, __global float * R_yy, __global float * R_xy, __global float * R_xz, __global float * R_yz, __global float * R_xx_lddrk, __global float * R_yy_lddrk, __global float * R_xy_lddrk, __global float * R_xz_lddrk, __global float * R_yz_lddrk, const float alpha_lddrk, const float beta_lddrk, const int use_lddrk, const __global float * restrict alphaval, const __global float * restrict betaval, const __global float * restrict gammaval, const __global float * restrict tau_sigmainvval, const __global float * restrict d_c11store, const __global float * restrict d_c12store, const __global float * restrict d_c13store, const __global float * restrict d_c14store, const __global float * restrict d_c15store, const __global float * restrict d_c16store, const __global float * restrict d_c22store, const __global float * restrict d_c23store, const __global float * restrict d_c24store, const __global float * restrict d_c25store, const __global float * restrict d_c26store, const __global float * restrict d_c33store, const __global float * restrict d_c34store, const __global float * restrict d_c35store, const __global float * restrict d_c36store, const __global float * restrict d_c44store, const __global float * restrict d_c45store, const __global float * restrict d_c46store, const __global float * restrict d_c55store, const __global float * restrict d_c56store, const __global float * restrict d_c66store, const __global float * restrict d_rstore, const int GRAVITY, const __global float * restrict d_gravity_pre_store, const __global float * restrict d_gravity_H, const __global float * restrict wgll_cube, const int NSPEC_CRUST_MANTLE_STRAIN_ONLY, __read_only image2d_t d_displ_cm_tex, __read_only image2d_t d_accel_cm_tex){\n\
#endif\n\
#else\n\
#ifdef USE_TEXTURES_CONSTANTS\n\
__kernel  void crust_mantle_impl_kernel_forward(const int nb_blocks_to_compute, const __global int * d_ibool, const __global int * d_ispec_is_tiso, const __global int * d_phase_ispec_inner, const int num_phase_ispec, const int d_iphase, const float deltat, const int use_mesh_coloring_gpu, const __global float * restrict d_displ, __global float * d_accel, const __global float * restrict d_xix, const __global float * restrict d_xiy, const __global float * restrict d_xiz, const __global float * restrict d_etax, const __global float * restrict d_etay, const __global float * restrict d_etaz, const __global float * restrict d_gammax, const __global float * restrict d_gammay, const __global float * restrict d_gammaz, const __global float * restrict d_hprime_xx, const __global float * restrict d_hprimewgll_xx, const __global float * restrict d_wgllwgll_xy, const __global float * restrict d_wgllwgll_xz, const __global float * restrict d_wgllwgll_yz, const __global float * restrict d_kappavstore, const __global float * restrict d_muvstore, const __global float * restrict d_kappahstore, const __global float * restrict d_muhstore, const __global float * restrict d_eta_anisostore, const int COMPUTE_AND_STORE_STRAIN, __global float * epsilondev_xx, __global float * epsilondev_yy, __global float * epsilondev_xy, __global float * epsilondev_xz, __global float * epsilondev_yz, __global float * epsilon_trace_over_3, const int ATTENUATION, const int PARTIAL_PHYS_DISPERSION_ONLY, const int USE_3D_ATTENUATION_ARRAYS, const __global float * restrict one_minus_sum_beta, const __global float * restrict factor_common, __global float * R_xx, __global float * R_yy, __global float * R_xy, __global float * R_xz, __global float * R_yz, __global float * R_xx_lddrk, __global float * R_yy_lddrk, __global float * R_xy_lddrk, __global float * R_xz_lddrk, __global float * R_yz_lddrk, const float alpha_lddrk, const float beta_lddrk, const int use_lddrk, const __global float * restrict alphaval, const __global float * restrict betaval, const __global float * restrict gammaval, const __global float * restrict tau_sigmainvval, const __global float * restrict d_c11store, const __global float * restrict d_c12store, const __global float * restrict d_c13store, const __global float * restrict d_c14store, const __global float * restrict d_c15store, const __global float * restrict d_c16store, const __global float * restrict d_c22store, const __global float * restrict d_c23store, const __global float * restrict d_c24store, const __global float * restrict d_c25store, const __global float * restrict d_c26store, const __global float * restrict d_c33store, const __global float * restrict d_c34store, const __global float * restrict d_c35store, const __global float * restrict d_c36store, const __global float * restrict d_c44store, const __global float * restrict d_c45store, const __global float * restrict d_c46store, const __global float * restrict d_c55store, const __global float * restrict d_c56store, const __global float * restrict d_c66store, const __global float * restrict d_rstore, const int GRAVITY, const __global float * restrict d_gravity_pre_store, const __global float * restrict d_gravity_H, const __global float * restrict wgll_cube, const int NSPEC_CRUST_MANTLE_STRAIN_ONLY, __read_only image2d_t d_hprime_xx_tex, __read_only image2d_t d_hprimewgll_xx_tex){\n\
#else\n\
__kernel  void crust_mantle_impl_kernel_forward(const int nb_blocks_to_compute, const __global int * d_ibool, const __global int * d_ispec_is_tiso, const __global int * d_phase_ispec_inner, const int num_phase_ispec, const int d_iphase, const float deltat, const int use_mesh_coloring_gpu, const __global float * restrict d_displ, __global float * d_accel, const __global float * restrict d_xix, const __global float * restrict d_xiy, const __global float * restrict d_xiz, const __global float * restrict d_etax, const __global float * restrict d_etay, const __global float * restrict d_etaz, const __global float * restrict d_gammax, const __global float * restrict d_gammay, const __global float * restrict d_gammaz, const __global float * restrict d_hprime_xx, const __global float * restrict d_hprimewgll_xx, const __global float * restrict d_wgllwgll_xy, const __global float * restrict d_wgllwgll_xz, const __global float * restrict d_wgllwgll_yz, const __global float * restrict d_kappavstore, const __global float * restrict d_muvstore, const __global float * restrict d_kappahstore, const __global float * restrict d_muhstore, const __global float * restrict d_eta_anisostore, const int COMPUTE_AND_STORE_STRAIN, __global float * epsilondev_xx, __global float * epsilondev_yy, __global float * epsilondev_xy, __global float * epsilondev_xz, __global float * epsilondev_yz, __global float * epsilon_trace_over_3, const int ATTENUATION, const int PARTIAL_PHYS_DISPERSION_ONLY, const int USE_3D_ATTENUATION_ARRAYS, const __global float * restrict one_minus_sum_beta, const __global float * restrict factor_common, __global float * R_xx, __global float * R_yy, __global float * R_xy, __global float * R_xz, __global float * R_yz, __global float * R_xx_lddrk, __global float * R_yy_lddrk, __global float * R_xy_lddrk, __global float * R_xz_lddrk, __global float * R_yz_lddrk, const float alpha_lddrk, const float beta_lddrk, const int use_lddrk, const __global float * restrict alphaval, const __global float * restrict betaval, const __global float * restrict gammaval, const __global float * restrict tau_sigmainvval, const __global float * restrict d_c11store, const __global float * restrict d_c12store, const __global float * restrict d_c13store, const __global float * restrict d_c14store, const __global float * restrict d_c15store, const __global float * restrict d_c16store, const __global float * restrict d_c22store, const __global float * restrict d_c23store, const __global float * restrict d_c24store, const __global float * restrict d_c25store, const __global float * restrict d_c26store, const __global float * restrict d_c33store, const __global float * restrict d_c34store, const __global float * restrict d_c35store, const __global float * restrict d_c36store, const __global float * restrict d_c44store, const __global float * restrict d_c45store, const __global float * restrict d_c46store, const __global float * restrict d_c55store, const __global float * restrict d_c56store, const __global float * restrict d_c66store, const __global float * restrict d_rstore, const int GRAVITY, const __global float * restrict d_gravity_pre_store, const __global float * restrict d_gravity_H, const __global float * restrict wgll_cube, const int NSPEC_CRUST_MANTLE_STRAIN_ONLY){\n\
#endif\n\
#endif\n\
#ifdef USE_TEXTURES_FIELDS\n\
  const sampler_t sampler_d_displ_cm_tex = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n\
  const sampler_t sampler_d_accel_cm_tex = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n\
#endif\n\
#ifdef USE_TEXTURES_CONSTANTS\n\
  const sampler_t sampler_d_hprime_xx_tex = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n\
  const sampler_t sampler_d_hprimewgll_xx_tex = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n\
#endif\n\
  int K;\n\
  int J;\n\
  int I;\n\
  ushort active_1;\n\
  int offset;\n\
  int iglob_1;\n\
  int working_element;\n\
  float fac1;\n\
  float fac2;\n\
  float fac3;\n\
  float sigma_xx;\n\
  float sigma_xy;\n\
  float sigma_xz;\n\
  float sigma_yx;\n\
  float sigma_yy;\n\
  float sigma_yz;\n\
  float sigma_zx;\n\
  float sigma_zy;\n\
  float sigma_zz;\n\
  float epsilondev_xx_loc_1;\n\
  float epsilondev_yy_loc_1;\n\
  float epsilondev_xy_loc_1;\n\
  float epsilondev_xz_loc_1;\n\
  float epsilondev_yz_loc_1;\n\
  float sum_terms1;\n\
  float sum_terms2;\n\
  float sum_terms3;\n\
  float rho_s_H_1_1;\n\
  float rho_s_H_1_2;\n\
  float rho_s_H_1_3;\n\
\n\
  // shared arrays\n\
  __local float s_dummyx_loc[(NGLL3)];\n\
  __local float s_dummyy_loc[(NGLL3)];\n\
  __local float s_dummyz_loc[(NGLL3)];\n\
\n\
  __local float s_tempx1[(NGLL3)];\n\
  __local float s_tempx2[(NGLL3)];\n\
  __local float s_tempx3[(NGLL3)];\n\
  __local float s_tempy1[(NGLL3)];\n\
  __local float s_tempy2[(NGLL3)];\n\
  __local float s_tempy3[(NGLL3)];\n\
  __local float s_tempz1[(NGLL3)];\n\
  __local float s_tempz2[(NGLL3)];\n\
  __local float s_tempz3[(NGLL3)];\n\
  __local float sh_hprime_xx[(NGLL2)];\n\
  __local float sh_hprimewgll_xx[(NGLL2)];\n\
\n\
  const int bx = (get_group_id(1)) * (get_num_groups(0)) + get_group_id(0);\n\
  if (bx >= nb_blocks_to_compute) {\n\
     return ;\n\
  }\n\
\n\
  const int tx = get_local_id(0);\n\
\n\
  active_1 = (tx < NGLL3 ? 1 : 0);\n\
\n\
  if (active_1) {\n\
#ifdef USE_MESH_COLORING_GPU\n\
    working_element = bx;\n\
#else\n\
    if (use_mesh_coloring_gpu) {\n\
      working_element = bx;\n\
    } else {\n\
      working_element = d_phase_ispec_inner[bx + (num_phase_ispec) * (d_iphase - (1))] - (1);\n\
    }\n\
#endif\n\
    iglob_1 = d_ibool[(working_element) * (NGLL3) + tx] - (1);\n\
#ifdef USE_TEXTURES_FIELDS\n\
    s_dummyx_loc[tx] = as_float(read_imageui(d_displ_cm_tex, sampler_d_displ_cm_tex, int2((iglob_1) * (3) + 0,0)).x);\n\
    s_dummyy_loc[tx] = as_float(read_imageui(d_displ_cm_tex, sampler_d_displ_cm_tex, int2((iglob_1) * (3) + 1,0)).x);\n\
    s_dummyz_loc[tx] = as_float(read_imageui(d_displ_cm_tex, sampler_d_displ_cm_tex, int2((iglob_1) * (3) + 2,0)).x);\n\
#else\n\
    s_dummyx_loc[tx] = d_displ[0 + (3) * (iglob_1)];\n\
    s_dummyy_loc[tx] = d_displ[1 + (3) * (iglob_1)];\n\
    s_dummyz_loc[tx] = d_displ[2 + (3) * (iglob_1)];\n\
#endif\n\
  }\n\
\n\
  if (tx < NGLL2) {\n\
#ifdef USE_TEXTURES_CONSTANTS\n\
    sh_hprime_xx[tx] = as_float(read_imageui(d_hprime_xx_tex, sampler_d_hprime_xx_tex, int2(tx,0)).x);\n\
    sh_hprimewgll_xx[tx] = as_float(read_imageui(d_hprimewgll_xx_tex, sampler_d_hprimewgll_xx_tex, int2(tx,0)).x);\n\
#else\n\
    sh_hprime_xx[tx] = d_hprime_xx[tx];\n\
    sh_hprimewgll_xx[tx] = d_hprimewgll_xx[tx];\n\
#endif\n\
  }\n\
\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
  K = (tx) / (NGLL2);\n\
  J = (tx - ((K) * (NGLL2))) / (NGLLX);\n\
  I = tx - ((K) * (NGLL2)) - ((J) * (NGLLX));\n\
\n\
  if (active_1) {\n\
    float tempx1l;\n\
    float tempx2l;\n\
    float tempx3l;\n\
    float tempy1l;\n\
    float tempy2l;\n\
    float tempy3l;\n\
    float tempz1l;\n\
    float tempz2l;\n\
    float tempz3l;\n\
    tempx1l = 0.0f;\n\
    tempx2l = 0.0f;\n\
    tempx3l = 0.0f;\n\
    tempy1l = 0.0f;\n\
    tempy2l = 0.0f;\n\
    tempy3l = 0.0f;\n\
    tempz1l = 0.0f;\n\
    tempz2l = 0.0f;\n\
    tempz3l = 0.0f;\n\
#ifdef MANUALLY_UNROLLED_LOOPS\n\
    fac1 = sh_hprime_xx[(0) * (NGLLX) + I];\n\
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 0]) * (fac1);\n\
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 0]) * (fac1);\n\
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 0]) * (fac1);\n\
    fac2 = sh_hprime_xx[(0) * (NGLLX) + J];\n\
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (0) * (NGLLX) + I]) * (fac2);\n\
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (0) * (NGLLX) + I]) * (fac2);\n\
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (0) * (NGLLX) + I]) * (fac2);\n\
    fac3 = sh_hprime_xx[(0) * (NGLLX) + K];\n\
    tempx3l = tempx3l + (s_dummyx_loc[(0) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempy3l = tempy3l + (s_dummyy_loc[(0) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempz3l = tempz3l + (s_dummyz_loc[(0) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    fac1 = sh_hprime_xx[(1) * (NGLLX) + I];\n\
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 1]) * (fac1);\n\
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 1]) * (fac1);\n\
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 1]) * (fac1);\n\
    fac2 = sh_hprime_xx[(1) * (NGLLX) + J];\n\
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (1) * (NGLLX) + I]) * (fac2);\n\
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (1) * (NGLLX) + I]) * (fac2);\n\
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (1) * (NGLLX) + I]) * (fac2);\n\
    fac3 = sh_hprime_xx[(1) * (NGLLX) + K];\n\
    tempx3l = tempx3l + (s_dummyx_loc[(1) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempy3l = tempy3l + (s_dummyy_loc[(1) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempz3l = tempz3l + (s_dummyz_loc[(1) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    fac1 = sh_hprime_xx[(2) * (NGLLX) + I];\n\
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 2]) * (fac1);\n\
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 2]) * (fac1);\n\
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 2]) * (fac1);\n\
    fac2 = sh_hprime_xx[(2) * (NGLLX) + J];\n\
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (2) * (NGLLX) + I]) * (fac2);\n\
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (2) * (NGLLX) + I]) * (fac2);\n\
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (2) * (NGLLX) + I]) * (fac2);\n\
    fac3 = sh_hprime_xx[(2) * (NGLLX) + K];\n\
    tempx3l = tempx3l + (s_dummyx_loc[(2) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempy3l = tempy3l + (s_dummyy_loc[(2) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempz3l = tempz3l + (s_dummyz_loc[(2) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    fac1 = sh_hprime_xx[(3) * (NGLLX) + I];\n\
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 3]) * (fac1);\n\
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 3]) * (fac1);\n\
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 3]) * (fac1);\n\
    fac2 = sh_hprime_xx[(3) * (NGLLX) + J];\n\
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (3) * (NGLLX) + I]) * (fac2);\n\
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (3) * (NGLLX) + I]) * (fac2);\n\
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (3) * (NGLLX) + I]) * (fac2);\n\
    fac3 = sh_hprime_xx[(3) * (NGLLX) + K];\n\
    tempx3l = tempx3l + (s_dummyx_loc[(3) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempy3l = tempy3l + (s_dummyy_loc[(3) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempz3l = tempz3l + (s_dummyz_loc[(3) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    fac1 = sh_hprime_xx[(4) * (NGLLX) + I];\n\
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 4]) * (fac1);\n\
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 4]) * (fac1);\n\
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 4]) * (fac1);\n\
    fac2 = sh_hprime_xx[(4) * (NGLLX) + J];\n\
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (4) * (NGLLX) + I]) * (fac2);\n\
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (4) * (NGLLX) + I]) * (fac2);\n\
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (4) * (NGLLX) + I]) * (fac2);\n\
    fac3 = sh_hprime_xx[(4) * (NGLLX) + K];\n\
    tempx3l = tempx3l + (s_dummyx_loc[(4) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempy3l = tempy3l + (s_dummyy_loc[(4) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempz3l = tempz3l + (s_dummyz_loc[(4) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
#else\n\
    for (int l = 0; l <= NGLLX - (1); l += 1) {\n\
      fac1 = sh_hprime_xx[(l) * (NGLLX) + I];\n\
      tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + l]) * (fac1);\n\
      tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + l]) * (fac1);\n\
      tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + l]) * (fac1);\n\
      fac2 = sh_hprime_xx[(l) * (NGLLX) + J];\n\
      tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (l) * (NGLLX) + I]) * (fac2);\n\
      tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (l) * (NGLLX) + I]) * (fac2);\n\
      tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (l) * (NGLLX) + I]) * (fac2);\n\
      fac3 = sh_hprime_xx[(l) * (NGLLX) + K];\n\
      tempx3l = tempx3l + (s_dummyx_loc[(l) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
      tempy3l = tempy3l + (s_dummyy_loc[(l) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
      tempz3l = tempz3l + (s_dummyz_loc[(l) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    }\n\
#endif\n\
\n\
    float xixl;\n\
    float xiyl;\n\
    float xizl;\n\
    float etaxl;\n\
    float etayl;\n\
    float etazl;\n\
    float gammaxl;\n\
    float gammayl;\n\
    float gammazl;\n\
    offset = (working_element) * (NGLL3_PADDED) + tx;\n\
    xixl = d_xix[offset];\n\
    etaxl = d_etax[offset];\n\
    gammaxl = d_gammax[offset];\n\
    xiyl = d_xiy[offset];\n\
    etayl = d_etay[offset];\n\
    gammayl = d_gammay[offset];\n\
    xizl = d_xiz[offset];\n\
    etazl = d_etaz[offset];\n\
    gammazl = d_gammaz[offset];\n\
\n\
    float duxdxl;\n\
    float duxdyl;\n\
    float duxdzl;\n\
    float duydxl;\n\
    float duydyl;\n\
    float duydzl;\n\
    float duzdxl;\n\
    float duzdyl;\n\
    float duzdzl;\n\
    duxdxl = (xixl) * (tempx1l) + (etaxl) * (tempx2l) + (gammaxl) * (tempx3l);\n\
    duxdyl = (xiyl) * (tempx1l) + (etayl) * (tempx2l) + (gammayl) * (tempx3l);\n\
    duxdzl = (xizl) * (tempx1l) + (etazl) * (tempx2l) + (gammazl) * (tempx3l);\n\
    duydxl = (xixl) * (tempy1l) + (etaxl) * (tempy2l) + (gammaxl) * (tempy3l);\n\
    duydyl = (xiyl) * (tempy1l) + (etayl) * (tempy2l) + (gammayl) * (tempy3l);\n\
    duydzl = (xizl) * (tempy1l) + (etazl) * (tempy2l) + (gammazl) * (tempy3l);\n\
    duzdxl = (xixl) * (tempz1l) + (etaxl) * (tempz2l) + (gammaxl) * (tempz3l);\n\
    duzdyl = (xiyl) * (tempz1l) + (etayl) * (tempz2l) + (gammayl) * (tempz3l);\n\
    duzdzl = (xizl) * (tempz1l) + (etazl) * (tempz2l) + (gammazl) * (tempz3l);\n\
\n\
    float duxdxl_plus_duydyl;\n\
    float duxdxl_plus_duzdzl;\n\
    float duydyl_plus_duzdzl;\n\
    float duxdyl_plus_duydxl;\n\
    float duzdxl_plus_duxdzl;\n\
    float duzdyl_plus_duydzl;\n\
    duxdxl_plus_duydyl = duxdxl + duydyl;\n\
    duxdxl_plus_duzdzl = duxdxl + duzdzl;\n\
    duydyl_plus_duzdzl = duydyl + duzdzl;\n\
    duxdyl_plus_duydxl = duxdyl + duydxl;\n\
    duzdxl_plus_duxdzl = duzdxl + duxdzl;\n\
    duzdyl_plus_duydzl = duzdyl + duydzl;\n\
\n\
    if (COMPUTE_AND_STORE_STRAIN) {\n\
      float templ;\n\
      templ = (duxdxl + duydyl + duzdzl) * (0.3333333333333333f);\n\
      epsilondev_xx_loc_1 = duxdxl - (templ);\n\
      epsilondev_yy_loc_1 = duydyl - (templ);\n\
      epsilondev_xy_loc_1 = (duxdyl_plus_duydxl) * (0.5f);\n\
      epsilondev_xz_loc_1 = (duzdxl_plus_duxdzl) * (0.5f);\n\
      epsilondev_yz_loc_1 = (duzdyl_plus_duydzl) * (0.5f);\n\
      if (NSPEC_CRUST_MANTLE_STRAIN_ONLY > 1) {\n\
        epsilon_trace_over_3[tx + (working_element) * (NGLL3)] = templ;\n\
      }\n\
    }\n\
\n\
    if ( !(d_ispec_is_tiso[working_element])) {\n\
      // isotropic elements\n\
      compute_element_cm_iso(offset, d_kappavstore, d_muvstore, duxdxl, duydyl, duzdzl, duxdxl_plus_duydyl, duxdxl_plus_duzdzl, duydyl_plus_duzdzl, duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl,  &sigma_xx,  &sigma_yy,  &sigma_zz,  &sigma_xy,  &sigma_xz,  &sigma_yz);\n\
    } else {\n\
      // transversely isotropic\n\
      compute_element_cm_tiso(offset, d_c11store, d_c12store, d_c13store, d_c14store, d_c15store, d_c16store, d_c22store, d_c23store, d_c24store, d_c25store, d_c26store, d_c33store, d_c34store, d_c35store, d_c36store, d_c44store, d_c45store, d_c46store, d_c55store, d_c56store, d_c66store, duxdxl, duxdyl, duxdzl, duydxl, duydyl, duydzl, duzdxl, duzdyl, duzdzl, duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl,  &sigma_xx,  &sigma_yy,  &sigma_zz,  &sigma_xy,  &sigma_xz,  &sigma_yz);\n\
    }\n\
\n\
    if (ATTENUATION &&  !(PARTIAL_PHYS_DISPERSION_ONLY)) {\n\
      compute_element_cm_att_stress(tx, working_element, R_xx, R_yy, R_xy, R_xz, R_yz,  &sigma_xx,  &sigma_yy,  &sigma_zz,  &sigma_xy,  &sigma_xz,  &sigma_yz);\n\
    }\n\
\n\
    sigma_yx = sigma_xy;\n\
    sigma_zx = sigma_xz;\n\
    sigma_zy = sigma_yz;\n\
\n\
    float jacobianl;\n\
    jacobianl = (1.0f) / ((xixl) * ((etayl) * (gammazl) - ((etazl) * (gammayl))) - ((xiyl) * ((etaxl) * (gammazl) - ((etazl) * (gammaxl)))) + (xizl) * ((etaxl) * (gammayl) - ((etayl) * (gammaxl))));\n\
\n\
    if (GRAVITY) {\n\
      compute_element_cm_gravity(tx, iglob_1, d_gravity_pre_store, d_gravity_H, wgll_cube, jacobianl, s_dummyx_loc, s_dummyy_loc, s_dummyz_loc,  &sigma_xx,  &sigma_yy,  &sigma_zz,  &sigma_xy,  &sigma_yx,  &sigma_xz,  &sigma_zx,  &sigma_yz,  &sigma_zy,  &rho_s_H_1_1,  &rho_s_H_1_2,  &rho_s_H_1_3);\n\
    }\n\
\n\
    s_tempx1[tx] = (jacobianl) * ((sigma_xx) * (xixl) + (sigma_yx) * (xiyl) + (sigma_zx) * (xizl));\n\
    s_tempy1[tx] = (jacobianl) * ((sigma_xy) * (xixl) + (sigma_yy) * (xiyl) + (sigma_zy) * (xizl));\n\
    s_tempz1[tx] = (jacobianl) * ((sigma_xz) * (xixl) + (sigma_yz) * (xiyl) + (sigma_zz) * (xizl));\n\
    s_tempx2[tx] = (jacobianl) * ((sigma_xx) * (etaxl) + (sigma_yx) * (etayl) + (sigma_zx) * (etazl));\n\
    s_tempy2[tx] = (jacobianl) * ((sigma_xy) * (etaxl) + (sigma_yy) * (etayl) + (sigma_zy) * (etazl));\n\
    s_tempz2[tx] = (jacobianl) * ((sigma_xz) * (etaxl) + (sigma_yz) * (etayl) + (sigma_zz) * (etazl));\n\
    s_tempx3[tx] = (jacobianl) * ((sigma_xx) * (gammaxl) + (sigma_yx) * (gammayl) + (sigma_zx) * (gammazl));\n\
    s_tempy3[tx] = (jacobianl) * ((sigma_xy) * (gammaxl) + (sigma_yy) * (gammayl) + (sigma_zy) * (gammazl));\n\
    s_tempz3[tx] = (jacobianl) * ((sigma_xz) * (gammaxl) + (sigma_yz) * (gammayl) + (sigma_zz) * (gammazl));\n\
  }\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
  if (active_1) {\n\
    float tempx1l;\n\
    float tempx2l;\n\
    float tempx3l;\n\
    float tempy1l;\n\
    float tempy2l;\n\
    float tempy3l;\n\
    float tempz1l;\n\
    float tempz2l;\n\
    float tempz3l;\n\
    tempx1l = 0.0f;\n\
    tempx2l = 0.0f;\n\
    tempx3l = 0.0f;\n\
    tempy1l = 0.0f;\n\
    tempy2l = 0.0f;\n\
    tempy3l = 0.0f;\n\
    tempz1l = 0.0f;\n\
    tempz2l = 0.0f;\n\
    tempz3l = 0.0f;\n\
#ifdef MANUALLY_UNROLLED_LOOPS\n\
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 0];\n\
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 0;\n\
    tempx1l = tempx1l + (s_tempx1[offset]) * (fac1);\n\
    tempy1l = tempy1l + (s_tempy1[offset]) * (fac1);\n\
    tempz1l = tempz1l + (s_tempz1[offset]) * (fac1);\n\
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 0];\n\
    offset = (K) * (NGLL2) + (0) * (NGLLX) + I;\n\
    tempx2l = tempx2l + (s_tempx2[offset]) * (fac2);\n\
    tempy2l = tempy2l + (s_tempy2[offset]) * (fac2);\n\
    tempz2l = tempz2l + (s_tempz2[offset]) * (fac2);\n\
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 0];\n\
    offset = (0) * (NGLL2) + (J) * (NGLLX) + I;\n\
    tempx3l = tempx3l + (s_tempx3[offset]) * (fac3);\n\
    tempy3l = tempy3l + (s_tempy3[offset]) * (fac3);\n\
    tempz3l = tempz3l + (s_tempz3[offset]) * (fac3);\n\
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 1];\n\
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 1;\n\
    tempx1l = tempx1l + (s_tempx1[offset]) * (fac1);\n\
    tempy1l = tempy1l + (s_tempy1[offset]) * (fac1);\n\
    tempz1l = tempz1l + (s_tempz1[offset]) * (fac1);\n\
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 1];\n\
    offset = (K) * (NGLL2) + (1) * (NGLLX) + I;\n\
    tempx2l = tempx2l + (s_tempx2[offset]) * (fac2);\n\
    tempy2l = tempy2l + (s_tempy2[offset]) * (fac2);\n\
    tempz2l = tempz2l + (s_tempz2[offset]) * (fac2);\n\
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 1];\n\
    offset = (1) * (NGLL2) + (J) * (NGLLX) + I;\n\
    tempx3l = tempx3l + (s_tempx3[offset]) * (fac3);\n\
    tempy3l = tempy3l + (s_tempy3[offset]) * (fac3);\n\
    tempz3l = tempz3l + (s_tempz3[offset]) * (fac3);\n\
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 2];\n\
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 2;\n\
    tempx1l = tempx1l + (s_tempx1[offset]) * (fac1);\n\
    tempy1l = tempy1l + (s_tempy1[offset]) * (fac1);\n\
    tempz1l = tempz1l + (s_tempz1[offset]) * (fac1);\n\
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 2];\n\
    offset = (K) * (NGLL2) + (2) * (NGLLX) + I;\n\
    tempx2l = tempx2l + (s_tempx2[offset]) * (fac2);\n\
    tempy2l = tempy2l + (s_tempy2[offset]) * (fac2);\n\
    tempz2l = tempz2l + (s_tempz2[offset]) * (fac2);\n\
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 2];\n\
    offset = (2) * (NGLL2) + (J) * (NGLLX) + I;\n\
    tempx3l = tempx3l + (s_tempx3[offset]) * (fac3);\n\
    tempy3l = tempy3l + (s_tempy3[offset]) * (fac3);\n\
    tempz3l = tempz3l + (s_tempz3[offset]) * (fac3);\n\
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 3];\n\
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 3;\n\
    tempx1l = tempx1l + (s_tempx1[offset]) * (fac1);\n\
    tempy1l = tempy1l + (s_tempy1[offset]) * (fac1);\n\
    tempz1l = tempz1l + (s_tempz1[offset]) * (fac1);\n\
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 3];\n\
    offset = (K) * (NGLL2) + (3) * (NGLLX) + I;\n\
    tempx2l = tempx2l + (s_tempx2[offset]) * (fac2);\n\
    tempy2l = tempy2l + (s_tempy2[offset]) * (fac2);\n\
    tempz2l = tempz2l + (s_tempz2[offset]) * (fac2);\n\
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 3];\n\
    offset = (3) * (NGLL2) + (J) * (NGLLX) + I;\n\
    tempx3l = tempx3l + (s_tempx3[offset]) * (fac3);\n\
    tempy3l = tempy3l + (s_tempy3[offset]) * (fac3);\n\
    tempz3l = tempz3l + (s_tempz3[offset]) * (fac3);\n\
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 4];\n\
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 4;\n\
    tempx1l = tempx1l + (s_tempx1[offset]) * (fac1);\n\
    tempy1l = tempy1l + (s_tempy1[offset]) * (fac1);\n\
    tempz1l = tempz1l + (s_tempz1[offset]) * (fac1);\n\
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 4];\n\
    offset = (K) * (NGLL2) + (4) * (NGLLX) + I;\n\
    tempx2l = tempx2l + (s_tempx2[offset]) * (fac2);\n\
    tempy2l = tempy2l + (s_tempy2[offset]) * (fac2);\n\
    tempz2l = tempz2l + (s_tempz2[offset]) * (fac2);\n\
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 4];\n\
    offset = (4) * (NGLL2) + (J) * (NGLLX) + I;\n\
    tempx3l = tempx3l + (s_tempx3[offset]) * (fac3);\n\
    tempy3l = tempy3l + (s_tempy3[offset]) * (fac3);\n\
    tempz3l = tempz3l + (s_tempz3[offset]) * (fac3);\n\
#else\n\
    for (int l = 0; l <= NGLLX - (1); l += 1) {\n\
      fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + l];\n\
      offset = (K) * (NGLL2) + (J) * (NGLLX) + l;\n\
      tempx1l = tempx1l + (s_tempx1[offset]) * (fac1);\n\
      tempy1l = tempy1l + (s_tempy1[offset]) * (fac1);\n\
      tempz1l = tempz1l + (s_tempz1[offset]) * (fac1);\n\
      fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + l];\n\
      offset = (K) * (NGLL2) + (l) * (NGLLX) + I;\n\
      tempx2l = tempx2l + (s_tempx2[offset]) * (fac2);\n\
      tempy2l = tempy2l + (s_tempy2[offset]) * (fac2);\n\
      tempz2l = tempz2l + (s_tempz2[offset]) * (fac2);\n\
      fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + l];\n\
      offset = (l) * (NGLL2) + (J) * (NGLLX) + I;\n\
      tempx3l = tempx3l + (s_tempx3[offset]) * (fac3);\n\
      tempy3l = tempy3l + (s_tempy3[offset]) * (fac3);\n\
      tempz3l = tempz3l + (s_tempz3[offset]) * (fac3);\n\
    }\n\
#endif\n\
\n\
    fac1 = d_wgllwgll_yz[(K) * (NGLLX) + J];\n\
    fac2 = d_wgllwgll_xz[(K) * (NGLLX) + I];\n\
    fac3 = d_wgllwgll_xy[(J) * (NGLLX) + I];\n\
    sum_terms1 =  -((fac1) * (tempx1l) + (fac2) * (tempx2l) + (fac3) * (tempx3l));\n\
    sum_terms2 =  -((fac1) * (tempy1l) + (fac2) * (tempy2l) + (fac3) * (tempy3l));\n\
    sum_terms3 =  -((fac1) * (tempz1l) + (fac2) * (tempz2l) + (fac3) * (tempz3l));\n\
\n\
    if (GRAVITY) {\n\
      sum_terms1 = sum_terms1 + rho_s_H_1_1;\n\
      sum_terms2 = sum_terms2 + rho_s_H_1_2;\n\
      sum_terms3 = sum_terms3 + rho_s_H_1_3;\n\
    }\n\
\n\
#ifdef USE_MESH_COLORING_GPU\n\
#ifdef USE_TEXTURES_FIELDS\n\
    d_accel[0 + (3) * (iglob_1)] = as_float(read_imageui(d_accel_cm_tex, sampler_d_accel_cm_tex, int2((iglob_1) * (3) + 0,0)).x) + sum_terms1;\n\
    d_accel[1 + (3) * (iglob_1)] = as_float(read_imageui(d_accel_cm_tex, sampler_d_accel_cm_tex, int2((iglob_1) * (3) + 1,0)).x) + sum_terms2;\n\
    d_accel[2 + (3) * (iglob_1)] = as_float(read_imageui(d_accel_cm_tex, sampler_d_accel_cm_tex, int2((iglob_1) * (3) + 2,0)).x) + sum_terms3;\n\
#else\n\
    d_accel[0 + (3) * (iglob_1)] = d_accel[0 + (3) * (iglob_1)] + sum_terms1;\n\
    d_accel[1 + (3) * (iglob_1)] = d_accel[1 + (3) * (iglob_1)] + sum_terms2;\n\
    d_accel[2 + (3) * (iglob_1)] = d_accel[2 + (3) * (iglob_1)] + sum_terms3;\n\
#endif\n\
#else\n\
    if (use_mesh_coloring_gpu) {\n\
#ifdef USE_TEXTURES_FIELDS\n\
      d_accel[0 + (3) * (iglob_1)] = as_float(read_imageui(d_accel_cm_tex, sampler_d_accel_cm_tex, int2((iglob_1) * (3) + 0,0)).x) + sum_terms1;\n\
      d_accel[1 + (3) * (iglob_1)] = as_float(read_imageui(d_accel_cm_tex, sampler_d_accel_cm_tex, int2((iglob_1) * (3) + 1,0)).x) + sum_terms2;\n\
      d_accel[2 + (3) * (iglob_1)] = as_float(read_imageui(d_accel_cm_tex, sampler_d_accel_cm_tex, int2((iglob_1) * (3) + 2,0)).x) + sum_terms3;\n\
#else\n\
      d_accel[0 + (3) * (iglob_1)] = d_accel[0 + (3) * (iglob_1)] + sum_terms1;\n\
      d_accel[1 + (3) * (iglob_1)] = d_accel[1 + (3) * (iglob_1)] + sum_terms2;\n\
      d_accel[2 + (3) * (iglob_1)] = d_accel[2 + (3) * (iglob_1)] + sum_terms3;\n\
#endif\n\
    } else {\n\
      atomicAdd(d_accel + (iglob_1) * (3) + 0, sum_terms1);\n\
      atomicAdd(d_accel + (iglob_1) * (3) + 1, sum_terms2);\n\
      atomicAdd(d_accel + (iglob_1) * (3) + 2, sum_terms3);\n\
    }\n\
#endif\n\
\n\
    if (ATTENUATION &&  !(PARTIAL_PHYS_DISPERSION_ONLY)) {\n\
      if ( !(use_lddrk)) {\n\
        compute_element_cm_att_memory(tx, working_element, d_muvstore, factor_common, alphaval, betaval, gammaval, R_xx, R_yy, R_xy, R_xz, R_yz, epsilondev_xx, epsilondev_yy, epsilondev_xy, epsilondev_xz, epsilondev_yz, epsilondev_xx_loc_1, epsilondev_yy_loc_1, epsilondev_xy_loc_1, epsilondev_xz_loc_1, epsilondev_yz_loc_1, USE_3D_ATTENUATION_ARRAYS);\n\
      } else {\n\
        compute_element_cm_att_memory_lddrk(tx, working_element, d_muvstore, factor_common, tau_sigmainvval, R_xx, R_yy, R_xy, R_xz, R_yz, R_xx_lddrk, R_yy_lddrk, R_xy_lddrk, R_xz_lddrk, R_yz_lddrk, alpha_lddrk, beta_lddrk, deltat, epsilondev_xx_loc_1, epsilondev_yy_loc_1, epsilondev_xy_loc_1, epsilondev_xz_loc_1, epsilondev_yz_loc_1, USE_3D_ATTENUATION_ARRAYS);\n\
      }\n\
    }\n\
\n\
    if (COMPUTE_AND_STORE_STRAIN) {\n\
      epsilondev_xx[tx + (working_element) * (NGLL3)] = epsilondev_xx_loc_1;\n\
      epsilondev_yy[tx + (working_element) * (NGLL3)] = epsilondev_yy_loc_1;\n\
      epsilondev_xy[tx + (working_element) * (NGLL3)] = epsilondev_xy_loc_1;\n\
      epsilondev_xz[tx + (working_element) * (NGLL3)] = epsilondev_xz_loc_1;\n\
      epsilondev_yz[tx + (working_element) * (NGLL3)] = epsilondev_yz_loc_1;\n\
    }\n\
  }\n\
}\n\
";

const char * crust_mantle_aniso_impl_kernel_forward_program = "\
inline void atomicAdd(volatile __global float *source, const float val) {\n\
  union {\n\
    unsigned int iVal;\n\
    float fVal;\n\
  } res, orig;\n\
  do {\n\
    orig.fVal = *source;\n\
    res.fVal = orig.fVal + val;\n\
  } while (atomic_cmpxchg((volatile __global unsigned int *)source, orig.iVal, res.iVal) != orig.iVal);\n\
}\n\
#ifndef INDEX2\n\
#define INDEX2(isize,i,j) i + isize*j\n\
#endif\n\
#ifndef INDEX3\n\
#define INDEX3(isize,jsize,i,j,k) i + isize*(j + jsize*k)\n\
#endif\n\
#ifndef INDEX4\n\
#define INDEX4(isize,jsize,ksize,i,j,k,x) i + isize*(j + jsize*(k + ksize*x))\n\
#endif\n\
#ifndef INDEX5\n\
#define INDEX5(isize,jsize,ksize,xsize,i,j,k,x,y) i + isize*(j + jsize*(k + ksize*(x + xsize*y)))\n\
#endif\n\
\n\
#ifndef NDIM\n\
#define NDIM 3\n\
#endif\n\
#ifndef NGLLX\n\
#define NGLLX 5\n\
#endif\n\
#ifndef NGLL2\n\
#define NGLL2 25\n\
#endif\n\
#ifndef NGLL3\n\
#define NGLL3 125\n\
#endif\n\
#ifndef NGLL3_PADDED\n\
#define NGLL3_PADDED 128\n\
#endif\n\
#ifndef N_SLS\n\
#define N_SLS 3\n\
#endif\n\
#ifndef IREGION_CRUST_MANTLE\n\
#define IREGION_CRUST_MANTLE 1\n\
#endif\n\
#ifndef IREGION_INNER_CORE\n\
#define IREGION_INNER_CORE 3\n\
#endif\n\
#ifndef IFLAG_IN_FICTITIOUS_CUBE\n\
#define IFLAG_IN_FICTITIOUS_CUBE 11\n\
#endif\n\
#ifndef COLORING_MIN_NSPEC_INNER_CORE\n\
#define COLORING_MIN_NSPEC_INNER_CORE 1000\n\
#endif\n\
#ifndef COLORING_MIN_NSPEC_OUTER_CORE\n\
#define COLORING_MIN_NSPEC_OUTER_CORE 1000\n\
#endif\n\
#ifndef BLOCKSIZE_TRANSFER\n\
#define BLOCKSIZE_TRANSFER 256\n\
#endif\n\
\n\
\n\
\n\
#if __OPENCL_C_VERSION__ && __OPENCL_C_VERSION__ >= 120\n\
static\n\
#endif\n\
void compute_element_cm_att_stress(const int tx, const int working_element, const __global float * R_xx, const __global float * R_yy, const __global float * R_xy, const __global float * R_xz, const __global float * R_yz, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_xz, float * sigma_yz){\n\
  int offset;\n\
  float R_xx_val;\n\
  float R_yy_val;\n\
\n\
  for (int i_sls = 0; i_sls < N_SLS; i_sls += 1) {\n\
\n\
    offset = tx + (NGLL3) * (i_sls + (N_SLS) * (working_element));\n\
\n\
    R_xx_val = R_xx[offset];\n\
    R_yy_val = R_yy[offset];\n\
\n\
    sigma_xx[0] = sigma_xx[0] - (R_xx_val);\n\
    sigma_yy[0] = sigma_yy[0] - (R_yy_val);\n\
    sigma_zz[0] = sigma_zz[0] + R_xx_val + R_yy_val;\n\
    sigma_xy[0] = sigma_xy[0] - (R_xy[offset]);\n\
    sigma_xz[0] = sigma_xz[0] - (R_xz[offset]);\n\
    sigma_yz[0] = sigma_yz[0] - (R_yz[offset]);\n\
  }\n\
}\n\
\n\
#if __OPENCL_C_VERSION__ && __OPENCL_C_VERSION__ >= 120\n\
static\n\
#endif\n\
void compute_element_cm_att_memory(const int tx, const int working_element, const __global float * d_muvstore, const __global float * factor_common, const __global float * alphaval, const __global float * betaval, const __global float * gammaval, __global float * R_xx, __global float * R_yy, __global float * R_xy, __global float * R_xz, __global float * R_yz, const __global float * epsilondev_xx, const __global float * epsilondev_yy, const __global float * epsilondev_xy, const __global float * epsilondev_xz, const __global float * epsilondev_yz, const float epsilondev_xx_loc, const float epsilondev_yy_loc, const float epsilondev_xy_loc, const float epsilondev_xz_loc, const float epsilondev_yz_loc, const int USE_3D_ATTENUATION_ARRAYS){\n\
  int offset;\n\
  float mul;\n\
  float factor_loc;\n\
  float sn;\n\
  float snp1;\n\
  float alphaval_loc;\n\
  float betaval_loc;\n\
  float gammaval_loc;\n\
\n\
  mul = d_muvstore[tx + (NGLL3_PADDED) * (working_element)];\n\
  const int offset_eps = tx + (NGLL3) * (working_element);\n\
\n\
  for (int i_sls = 0; i_sls < N_SLS; i_sls += 1) {\n\
    offset = tx + (NGLL3) * (i_sls + (N_SLS) * (working_element));\n\
    if (USE_3D_ATTENUATION_ARRAYS) {\n\
      factor_loc = (mul) * (factor_common[offset]);\n\
    } else {\n\
      factor_loc = (mul) * (factor_common[i_sls + (N_SLS) * (working_element)]);\n\
    }\n\
\n\
    alphaval_loc = alphaval[i_sls];\n\
    betaval_loc = betaval[i_sls];\n\
    gammaval_loc = gammaval[i_sls];\n\
\n\
    sn = (factor_loc) * (epsilondev_xx[offset_eps]);\n\
    snp1 = (factor_loc) * (epsilondev_xx_loc);\n\
    R_xx[offset] = (alphaval_loc) * (R_xx[offset]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);\n\
    sn = (factor_loc) * (epsilondev_yy[offset_eps]);\n\
    snp1 = (factor_loc) * (epsilondev_yy_loc);\n\
    R_yy[offset] = (alphaval_loc) * (R_yy[offset]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);\n\
    sn = (factor_loc) * (epsilondev_xy[offset_eps]);\n\
    snp1 = (factor_loc) * (epsilondev_xy_loc);\n\
    R_xy[offset] = (alphaval_loc) * (R_xy[offset]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);\n\
    sn = (factor_loc) * (epsilondev_xz[offset_eps]);\n\
    snp1 = (factor_loc) * (epsilondev_xz_loc);\n\
    R_xz[offset] = (alphaval_loc) * (R_xz[offset]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);\n\
    sn = (factor_loc) * (epsilondev_yz[offset_eps]);\n\
    snp1 = (factor_loc) * (epsilondev_yz_loc);\n\
    R_yz[offset] = (alphaval_loc) * (R_yz[offset]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);\n\
  }\n\
}\n\
\n\
#if __OPENCL_C_VERSION__ && __OPENCL_C_VERSION__ >= 120\n\
static\n\
#endif\n\
void compute_element_cm_att_memory_lddrk(const int tx, const int working_element, const __global float * d_muvstore, const __global float * factor_common, const __global float * tau_sigmainvval, __global float * R_xx, __global float * R_yy, __global float * R_xy, __global float * R_xz, __global float * R_yz, __global float * R_xx_lddrk, __global float * R_yy_lddrk, __global float * R_xy_lddrk, __global float * R_xz_lddrk, __global float * R_yz_lddrk, const float alpha_lddrk, const float beta_lddrk, const float deltat, const float epsilondev_xx_loc, const float epsilondev_yy_loc, const float epsilondev_xy_loc, const float epsilondev_xz_loc, const float epsilondev_yz_loc, const int USE_3D_ATTENUATION_ARRAYS){\n\
  int offset;\n\
  float mul;\n\
  float factor_loc;\n\
  float sn;\n\
  float snp1;\n\
  float tau_sigmainv_loc;\n\
\n\
  mul = d_muvstore[tx + (NGLL3_PADDED) * (working_element)];\n\
\n\
  for (int i_sls = 0; i_sls < N_SLS; i_sls += 1) {\n\
    offset = tx + (NGLL3) * (i_sls + (N_SLS) * (working_element));\n\
    if (USE_3D_ATTENUATION_ARRAYS) {\n\
      factor_loc = (mul) * (factor_common[offset]);\n\
    } else {\n\
      factor_loc = (mul) * (factor_common[i_sls + (N_SLS) * (working_element)]);\n\
    }\n\
\n\
    tau_sigmainv_loc = tau_sigmainvval[i_sls];\n\
\n\
    sn = (tau_sigmainv_loc) * (R_xx[offset]);\n\
    snp1 = (factor_loc) * (epsilondev_xx_loc);\n\
    R_xx_lddrk[offset] = (alpha_lddrk) * (R_xx_lddrk[offset]) + (deltat) * (snp1 - (sn));\n\
    R_xx[offset] = R_xx[offset] + (beta_lddrk) * (R_xx_lddrk[offset]);\n\
    sn = (tau_sigmainv_loc) * (R_yy[offset]);\n\
    snp1 = (factor_loc) * (epsilondev_yy_loc);\n\
    R_yy_lddrk[offset] = (alpha_lddrk) * (R_yy_lddrk[offset]) + (deltat) * (snp1 - (sn));\n\
    R_yy[offset] = R_yy[offset] + (beta_lddrk) * (R_yy_lddrk[offset]);\n\
    sn = (tau_sigmainv_loc) * (R_xy[offset]);\n\
    snp1 = (factor_loc) * (epsilondev_xy_loc);\n\
    R_xy_lddrk[offset] = (alpha_lddrk) * (R_xy_lddrk[offset]) + (deltat) * (snp1 - (sn));\n\
    R_xy[offset] = R_xy[offset] + (beta_lddrk) * (R_xy_lddrk[offset]);\n\
    sn = (tau_sigmainv_loc) * (R_xz[offset]);\n\
    snp1 = (factor_loc) * (epsilondev_xz_loc);\n\
    R_xz_lddrk[offset] = (alpha_lddrk) * (R_xz_lddrk[offset]) + (deltat) * (snp1 - (sn));\n\
    R_xz[offset] = R_xz[offset] + (beta_lddrk) * (R_xz_lddrk[offset]);\n\
    sn = (tau_sigmainv_loc) * (R_yz[offset]);\n\
    snp1 = (factor_loc) * (epsilondev_yz_loc);\n\
    R_yz_lddrk[offset] = (alpha_lddrk) * (R_yz_lddrk[offset]) + (deltat) * (snp1 - (sn));\n\
    R_yz[offset] = R_yz[offset] + (beta_lddrk) * (R_yz_lddrk[offset]);\n\
  }\n\
}\n\
\n\
#if __OPENCL_C_VERSION__ && __OPENCL_C_VERSION__ >= 120\n\
static\n\
#endif\n\
void compute_element_cm_gravity(const int tx, const int iglob, const __global float * restrict gravity_pre_store, const __global float * restrict gravity_H, const __global float * restrict wgll_cube, const float jacobianl, const __local float * s_dummyx_loc, const __local float * s_dummyy_loc, const __local float * s_dummyz_loc, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_yx, float * sigma_xz, float * sigma_zx, float * sigma_yz, float * sigma_zy, float * rho_s_H1, float * rho_s_H2, float * rho_s_H3){\n\
  float gxl;\n\
  float gyl;\n\
  float gzl;\n\
  float Hxxl;\n\
  float Hyyl;\n\
  float Hzzl;\n\
  float Hxyl;\n\
  float Hxzl;\n\
  float Hyzl;\n\
  float sx_l;\n\
  float sy_l;\n\
  float sz_l;\n\
  float factor;\n\
\n\
  gxl = gravity_pre_store[0 + (3) * (iglob)];\n\
  gyl = gravity_pre_store[1 + (3) * (iglob)];\n\
  gzl = gravity_pre_store[2 + (3) * (iglob)];\n\
\n\
  Hxxl = gravity_H[0 + (6) * (iglob)];\n\
  Hyyl = gravity_H[1 + (6) * (iglob)];\n\
  Hzzl = gravity_H[2 + (6) * (iglob)];\n\
  Hxyl = gravity_H[3 + (6) * (iglob)];\n\
  Hxzl = gravity_H[4 + (6) * (iglob)];\n\
  Hyzl = gravity_H[5 + (6) * (iglob)];\n\
\n\
  sx_l = s_dummyx_loc[tx];\n\
  sy_l = s_dummyy_loc[tx];\n\
  sz_l = s_dummyz_loc[tx];\n\
\n\
  *(sigma_xx) = *(sigma_xx) + (sy_l) * (gyl) + (sz_l) * (gzl);\n\
  *(sigma_yy) = *(sigma_yy) + (sx_l) * (gxl) + (sz_l) * (gzl);\n\
  *(sigma_zz) = *(sigma_zz) + (sx_l) * (gxl) + (sy_l) * (gyl);\n\
  *(sigma_xy) = *(sigma_xy) - ((sx_l) * (gyl));\n\
  *(sigma_yx) = *(sigma_yx) - ((sy_l) * (gxl));\n\
  *(sigma_xz) = *(sigma_xz) - ((sx_l) * (gzl));\n\
  *(sigma_zx) = *(sigma_zx) - ((sz_l) * (gxl));\n\
  *(sigma_yz) = *(sigma_yz) - ((sy_l) * (gzl));\n\
  *(sigma_zy) = *(sigma_zy) - ((sz_l) * (gyl));\n\
\n\
  factor = (jacobianl) * (wgll_cube[tx]);\n\
  *(rho_s_H1) = (factor) * ((sx_l) * (Hxxl) + (sy_l) * (Hxyl) + (sz_l) * (Hxzl));\n\
  *(rho_s_H2) = (factor) * ((sx_l) * (Hxyl) + (sy_l) * (Hyyl) + (sz_l) * (Hyzl));\n\
  *(rho_s_H3) = (factor) * ((sx_l) * (Hxzl) + (sy_l) * (Hyzl) + (sz_l) * (Hzzl));\n\
}\n\
\n\
#if __OPENCL_C_VERSION__ && __OPENCL_C_VERSION__ >= 120\n\
static\n\
#endif\n\
void compute_element_cm_aniso(const int offset, const __global float * d_c11store, const __global float * d_c12store, const __global float * d_c13store, const __global float * d_c14store, const __global float * d_c15store, const __global float * d_c16store, const __global float * d_c22store, const __global float * d_c23store, const __global float * d_c24store, const __global float * d_c25store, const __global float * d_c26store, const __global float * d_c33store, const __global float * d_c34store, const __global float * d_c35store, const __global float * d_c36store, const __global float * d_c44store, const __global float * d_c45store, const __global float * d_c46store, const __global float * d_c55store, const __global float * d_c56store, const __global float * d_c66store, const float duxdxl, const float duxdyl, const float duxdzl, const float duydxl, const float duydyl, const float duydzl, const float duzdxl, const float duzdyl, const float duzdzl, const float duxdyl_plus_duydxl, const float duzdxl_plus_duxdzl, const float duzdyl_plus_duydzl, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_xz, float * sigma_yz){\n\
  float c11;\n\
  float c12;\n\
  float c13;\n\
  float c14;\n\
  float c15;\n\
  float c16;\n\
  float c22;\n\
  float c23;\n\
  float c24;\n\
  float c25;\n\
  float c26;\n\
  float c33;\n\
  float c34;\n\
  float c35;\n\
  float c36;\n\
  float c44;\n\
  float c45;\n\
  float c46;\n\
  float c55;\n\
  float c56;\n\
  float c66;\n\
\n\
  c11 = d_c11store[offset];\n\
  c12 = d_c12store[offset];\n\
  c13 = d_c13store[offset];\n\
  c14 = d_c14store[offset];\n\
  c15 = d_c15store[offset];\n\
  c16 = d_c16store[offset];\n\
  c22 = d_c22store[offset];\n\
  c23 = d_c23store[offset];\n\
  c24 = d_c24store[offset];\n\
  c25 = d_c25store[offset];\n\
  c26 = d_c26store[offset];\n\
  c33 = d_c33store[offset];\n\
  c34 = d_c34store[offset];\n\
  c35 = d_c35store[offset];\n\
  c36 = d_c36store[offset];\n\
  c44 = d_c44store[offset];\n\
  c45 = d_c45store[offset];\n\
  c46 = d_c46store[offset];\n\
  c55 = d_c55store[offset];\n\
  c56 = d_c56store[offset];\n\
  c66 = d_c66store[offset];\n\
\n\
  *(sigma_xx) = (c11) * (duxdxl) + (c16) * (duxdyl_plus_duydxl) + (c12) * (duydyl) + (c15) * (duzdxl_plus_duxdzl) + (c14) * (duzdyl_plus_duydzl) + (c13) * (duzdzl);\n\
  *(sigma_yy) = (c12) * (duxdxl) + (c26) * (duxdyl_plus_duydxl) + (c22) * (duydyl) + (c25) * (duzdxl_plus_duxdzl) + (c24) * (duzdyl_plus_duydzl) + (c23) * (duzdzl);\n\
  *(sigma_zz) = (c13) * (duxdxl) + (c36) * (duxdyl_plus_duydxl) + (c23) * (duydyl) + (c35) * (duzdxl_plus_duxdzl) + (c34) * (duzdyl_plus_duydzl) + (c33) * (duzdzl);\n\
  *(sigma_xy) = (c16) * (duxdxl) + (c66) * (duxdyl_plus_duydxl) + (c26) * (duydyl) + (c56) * (duzdxl_plus_duxdzl) + (c46) * (duzdyl_plus_duydzl) + (c36) * (duzdzl);\n\
  *(sigma_xz) = (c15) * (duxdxl) + (c56) * (duxdyl_plus_duydxl) + (c25) * (duydyl) + (c55) * (duzdxl_plus_duxdzl) + (c45) * (duzdyl_plus_duydzl) + (c35) * (duzdzl);\n\
  *(sigma_yz) = (c14) * (duxdxl) + (c46) * (duxdyl_plus_duydxl) + (c24) * (duydyl) + (c45) * (duzdxl_plus_duxdzl) + (c44) * (duzdyl_plus_duydzl) + (c34) * (duzdzl);\n\
}\n\
\n\
#if __OPENCL_C_VERSION__ && __OPENCL_C_VERSION__ >= 120\n\
static\n\
#endif\n\
void compute_element_cm_iso(const int offset, const __global float * d_kappavstore, const __global float * d_muvstore, const float duxdxl, const float duydyl, const float duzdzl, const float duxdxl_plus_duydyl, const float duxdxl_plus_duzdzl, const float duydyl_plus_duzdzl, const float duxdyl_plus_duydxl, const float duzdxl_plus_duxdzl, const float duzdyl_plus_duydzl, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_xz, float * sigma_yz){\n\
  float lambdal;\n\
  float mul;\n\
  float lambdalplus2mul;\n\
  float kappal;\n\
\n\
  mul = d_muvstore[offset];\n\
  kappal = d_kappavstore[offset];\n\
\n\
  lambdalplus2mul = kappal + (mul) * (1.3333333333333333f);\n\
  lambdal = lambdalplus2mul - ((mul) * (2.0f));\n\
\n\
  *(sigma_xx) = (lambdalplus2mul) * (duxdxl) + (lambdal) * (duydyl_plus_duzdzl);\n\
  *(sigma_yy) = (lambdalplus2mul) * (duydyl) + (lambdal) * (duxdxl_plus_duzdzl);\n\
  *(sigma_zz) = (lambdalplus2mul) * (duzdzl) + (lambdal) * (duxdxl_plus_duydyl);\n\
  *(sigma_xy) = (mul) * (duxdyl_plus_duydxl);\n\
  *(sigma_xz) = (mul) * (duzdxl_plus_duxdzl);\n\
  *(sigma_yz) = (mul) * (duzdyl_plus_duydzl);\n\
}\n\
\n\
#if __OPENCL_C_VERSION__ && __OPENCL_C_VERSION__ >= 120\n\
static\n\
#endif\n\
void compute_element_cm_tiso(const int offset, const __global float * d_c11store, const __global float * d_c12store, const __global float * d_c13store, const __global float * d_c14store, const __global float * d_c15store, const __global float * d_c16store, const __global float * d_c22store, const __global float * d_c23store, const __global float * d_c24store, const __global float * d_c25store, const __global float * d_c26store, const __global float * d_c33store, const __global float * d_c34store, const __global float * d_c35store, const __global float * d_c36store, const __global float * d_c44store, const __global float * d_c45store, const __global float * d_c46store, const __global float * d_c55store, const __global float * d_c56store, const __global float * d_c66store, const float duxdxl, const float duxdyl, const float duxdzl, const float duydxl, const float duydyl, const float duydzl, const float duzdxl, const float duzdyl, const float duzdzl, const float duxdyl_plus_duydxl, const float duzdxl_plus_duxdzl, const float duzdyl_plus_duydzl, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_xz, float * sigma_yz){\n\
  float c11;\n\
  float c12;\n\
  float c13;\n\
  float c14;\n\
  float c15;\n\
  float c16;\n\
  float c22;\n\
  float c23;\n\
  float c24;\n\
  float c25;\n\
  float c26;\n\
  float c33;\n\
  float c34;\n\
  float c35;\n\
  float c36;\n\
  float c44;\n\
  float c45;\n\
  float c46;\n\
  float c55;\n\
  float c56;\n\
  float c66;\n\
\n\
  c11 = d_c11store[offset];\n\
  c12 = d_c12store[offset];\n\
  c13 = d_c13store[offset];\n\
  c14 = d_c14store[offset];\n\
  c15 = d_c15store[offset];\n\
  c16 = d_c16store[offset];\n\
  c22 = d_c22store[offset];\n\
  c23 = d_c23store[offset];\n\
  c24 = d_c24store[offset];\n\
  c25 = d_c25store[offset];\n\
  c26 = d_c26store[offset];\n\
  c33 = d_c33store[offset];\n\
  c34 = d_c34store[offset];\n\
  c35 = d_c35store[offset];\n\
  c36 = d_c36store[offset];\n\
  c44 = d_c44store[offset];\n\
  c45 = d_c45store[offset];\n\
  c46 = d_c46store[offset];\n\
  c55 = d_c55store[offset];\n\
  c56 = d_c56store[offset];\n\
  c66 = d_c66store[offset];\n\
\n\
  *(sigma_xx) = (c11) * (duxdxl) + (c16) * (duxdyl_plus_duydxl) + (c12) * (duydyl) + (c15) * (duzdxl_plus_duxdzl) + (c14) * (duzdyl_plus_duydzl) + (c13) * (duzdzl);\n\
  *(sigma_yy) = (c12) * (duxdxl) + (c26) * (duxdyl_plus_duydxl) + (c22) * (duydyl) + (c25) * (duzdxl_plus_duxdzl) + (c24) * (duzdyl_plus_duydzl) + (c23) * (duzdzl);\n\
  *(sigma_zz) = (c13) * (duxdxl) + (c36) * (duxdyl_plus_duydxl) + (c23) * (duydyl) + (c35) * (duzdxl_plus_duxdzl) + (c34) * (duzdyl_plus_duydzl) + (c33) * (duzdzl);\n\
  *(sigma_xy) = (c16) * (duxdxl) + (c66) * (duxdyl_plus_duydxl) + (c26) * (duydyl) + (c56) * (duzdxl_plus_duxdzl) + (c46) * (duzdyl_plus_duydzl) + (c36) * (duzdzl);\n\
  *(sigma_xz) = (c15) * (duxdxl) + (c56) * (duxdyl_plus_duydxl) + (c25) * (duydyl) + (c55) * (duzdxl_plus_duxdzl) + (c45) * (duzdyl_plus_duydzl) + (c35) * (duzdzl);\n\
  *(sigma_yz) = (c14) * (duxdxl) + (c46) * (duxdyl_plus_duydxl) + (c24) * (duydyl) + (c45) * (duzdxl_plus_duxdzl) + (c44) * (duzdyl_plus_duydzl) + (c34) * (duzdzl);\n\
}\n\
\n\
\n\
/*----------------------------------------------*/\n\
// main function\n\
/*----------------------------------------------*/\n\
\n\
#ifdef USE_TEXTURES_FIELDS\n\
#ifdef USE_TEXTURES_CONSTANTS\n\
__kernel  void crust_mantle_aniso_impl_kernel_forward(const int nb_blocks_to_compute, const __global int * d_ibool, const __global int * d_ispec_is_tiso, const __global int * d_phase_ispec_inner, const int num_phase_ispec, const int d_iphase, const float deltat, const int use_mesh_coloring_gpu, const __global float * restrict d_displ, __global float * d_accel, const __global float * restrict d_xix, const __global float * restrict d_xiy, const __global float * restrict d_xiz, const __global float * restrict d_etax, const __global float * restrict d_etay, const __global float * restrict d_etaz, const __global float * restrict d_gammax, const __global float * restrict d_gammay, const __global float * restrict d_gammaz, const __global float * restrict d_hprime_xx, const __global float * restrict d_hprimewgll_xx, const __global float * restrict d_wgllwgll_xy, const __global float * restrict d_wgllwgll_xz, const __global float * restrict d_wgllwgll_yz, const __global float * restrict d_muvstore, const int COMPUTE_AND_STORE_STRAIN, __global float * epsilondev_xx, __global float * epsilondev_yy, __global float * epsilondev_xy, __global float * epsilondev_xz, __global float * epsilondev_yz, __global float * epsilon_trace_over_3, const int ATTENUATION, const int PARTIAL_PHYS_DISPERSION_ONLY, const int USE_3D_ATTENUATION_ARRAYS, const __global float * restrict one_minus_sum_beta, const __global float * restrict factor_common, __global float * R_xx, __global float * R_yy, __global float * R_xy, __global float * R_xz, __global float * R_yz, __global float * R_xx_lddrk, __global float * R_yy_lddrk, __global float * R_xy_lddrk, __global float * R_xz_lddrk, __global float * R_yz_lddrk, const float alpha_lddrk, const float beta_lddrk, const int use_lddrk, const __global float * restrict alphaval, const __global float * restrict betaval, const __global float * restrict gammaval, const __global float * restrict tau_sigmainvval, const __global float * restrict d_c11store, const __global float * restrict d_c12store, const __global float * restrict d_c13store, const __global float * restrict d_c14store, const __global float * restrict d_c15store, const __global float * restrict d_c16store, const __global float * restrict d_c22store, const __global float * restrict d_c23store, const __global float * restrict d_c24store, const __global float * restrict d_c25store, const __global float * restrict d_c26store, const __global float * restrict d_c33store, const __global float * restrict d_c34store, const __global float * restrict d_c35store, const __global float * restrict d_c36store, const __global float * restrict d_c44store, const __global float * restrict d_c45store, const __global float * restrict d_c46store, const __global float * restrict d_c55store, const __global float * restrict d_c56store, const __global float * restrict d_c66store, const int GRAVITY, const __global float * restrict d_gravity_pre_store, const __global float * restrict d_gravity_H, const __global float * restrict wgll_cube, const int NSPEC_CRUST_MANTLE_STRAIN_ONLY, __read_only image2d_t d_displ_cm_tex, __read_only image2d_t d_accel_cm_tex, __read_only image2d_t d_hprime_xx_tex, __read_only image2d_t d_hprimewgll_xx_tex){\n\
#else\n\
__kernel  void crust_mantle_aniso_impl_kernel_forward(const int nb_blocks_to_compute, const __global int * d_ibool, const __global int * d_ispec_is_tiso, const __global int * d_phase_ispec_inner, const int num_phase_ispec, const int d_iphase, const float deltat, const int use_mesh_coloring_gpu, const __global float * restrict d_displ, __global float * d_accel, const __global float * restrict d_xix, const __global float * restrict d_xiy, const __global float * restrict d_xiz, const __global float * restrict d_etax, const __global float * restrict d_etay, const __global float * restrict d_etaz, const __global float * restrict d_gammax, const __global float * restrict d_gammay, const __global float * restrict d_gammaz, const __global float * restrict d_hprime_xx, const __global float * restrict d_hprimewgll_xx, const __global float * restrict d_wgllwgll_xy, const __global float * restrict d_wgllwgll_xz, const __global float * restrict d_wgllwgll_yz, const __global float * restrict d_muvstore, const int COMPUTE_AND_STORE_STRAIN, __global float * epsilondev_xx, __global float * epsilondev_yy, __global float * epsilondev_xy, __global float * epsilondev_xz, __global float * epsilondev_yz, __global float * epsilon_trace_over_3, const int ATTENUATION, const int PARTIAL_PHYS_DISPERSION_ONLY, const int USE_3D_ATTENUATION_ARRAYS, const __global float * restrict one_minus_sum_beta, const __global float * restrict factor_common, __global float * R_xx, __global float * R_yy, __global float * R_xy, __global float * R_xz, __global float * R_yz, __global float * R_xx_lddrk, __global float * R_yy_lddrk, __global float * R_xy_lddrk, __global float * R_xz_lddrk, __global float * R_yz_lddrk, const float alpha_lddrk, const float beta_lddrk, const int use_lddrk, const __global float * restrict alphaval, const __global float * restrict betaval, const __global float * restrict gammaval, const __global float * restrict tau_sigmainvval, const __global float * restrict d_c11store, const __global float * restrict d_c12store, const __global float * restrict d_c13store, const __global float * restrict d_c14store, const __global float * restrict d_c15store, const __global float * restrict d_c16store, const __global float * restrict d_c22store, const __global float * restrict d_c23store, const __global float * restrict d_c24store, const __global float * restrict d_c25store, const __global float * restrict d_c26store, const __global float * restrict d_c33store, const __global float * restrict d_c34store, const __global float * restrict d_c35store, const __global float * restrict d_c36store, const __global float * restrict d_c44store, const __global float * restrict d_c45store, const __global float * restrict d_c46store, const __global float * restrict d_c55store, const __global float * restrict d_c56store, const __global float * restrict d_c66store, const int GRAVITY, const __global float * restrict d_gravity_pre_store, const __global float * restrict d_gravity_H, const __global float * restrict wgll_cube, const int NSPEC_CRUST_MANTLE_STRAIN_ONLY, __read_only image2d_t d_displ_cm_tex, __read_only image2d_t d_accel_cm_tex){\n\
#endif\n\
#else\n\
#ifdef USE_TEXTURES_CONSTANTS\n\
__kernel  void crust_mantle_aniso_impl_kernel_forward(const int nb_blocks_to_compute, const __global int * d_ibool, const __global int * d_ispec_is_tiso, const __global int * d_phase_ispec_inner, const int num_phase_ispec, const int d_iphase, const float deltat, const int use_mesh_coloring_gpu, const __global float * restrict d_displ, __global float * d_accel, const __global float * restrict d_xix, const __global float * restrict d_xiy, const __global float * restrict d_xiz, const __global float * restrict d_etax, const __global float * restrict d_etay, const __global float * restrict d_etaz, const __global float * restrict d_gammax, const __global float * restrict d_gammay, const __global float * restrict d_gammaz, const __global float * restrict d_hprime_xx, const __global float * restrict d_hprimewgll_xx, const __global float * restrict d_wgllwgll_xy, const __global float * restrict d_wgllwgll_xz, const __global float * restrict d_wgllwgll_yz, const __global float * restrict d_muvstore, const int COMPUTE_AND_STORE_STRAIN, __global float * epsilondev_xx, __global float * epsilondev_yy, __global float * epsilondev_xy, __global float * epsilondev_xz, __global float * epsilondev_yz, __global float * epsilon_trace_over_3, const int ATTENUATION, const int PARTIAL_PHYS_DISPERSION_ONLY, const int USE_3D_ATTENUATION_ARRAYS, const __global float * restrict one_minus_sum_beta, const __global float * restrict factor_common, __global float * R_xx, __global float * R_yy, __global float * R_xy, __global float * R_xz, __global float * R_yz, __global float * R_xx_lddrk, __global float * R_yy_lddrk, __global float * R_xy_lddrk, __global float * R_xz_lddrk, __global float * R_yz_lddrk, const float alpha_lddrk, const float beta_lddrk, const int use_lddrk, const __global float * restrict alphaval, const __global float * restrict betaval, const __global float * restrict gammaval, const __global float * restrict tau_sigmainvval, const __global float * restrict d_c11store, const __global float * restrict d_c12store, const __global float * restrict d_c13store, const __global float * restrict d_c14store, const __global float * restrict d_c15store, const __global float * restrict d_c16store, const __global float * restrict d_c22store, const __global float * restrict d_c23store, const __global float * restrict d_c24store, const __global float * restrict d_c25store, const __global float * restrict d_c26store, const __global float * restrict d_c33store, const __global float * restrict d_c34store, const __global float * restrict d_c35store, const __global float * restrict d_c36store, const __global float * restrict d_c44store, const __global float * restrict d_c45store, const __global float * restrict d_c46store, const __global float * restrict d_c55store, const __global float * restrict d_c56store, const __global float * restrict d_c66store, const int GRAVITY, const __global float * restrict d_gravity_pre_store, const __global float * restrict d_gravity_H, const __global float * restrict wgll_cube, const int NSPEC_CRUST_MANTLE_STRAIN_ONLY, __read_only image2d_t d_hprime_xx_tex, __read_only image2d_t d_hprimewgll_xx_tex){\n\
#else\n\
__kernel  void crust_mantle_aniso_impl_kernel_forward(const int nb_blocks_to_compute, const __global int * d_ibool, const __global int * d_ispec_is_tiso, const __global int * d_phase_ispec_inner, const int num_phase_ispec, const int d_iphase, const float deltat, const int use_mesh_coloring_gpu, const __global float * restrict d_displ, __global float * d_accel, const __global float * restrict d_xix, const __global float * restrict d_xiy, const __global float * restrict d_xiz, const __global float * restrict d_etax, const __global float * restrict d_etay, const __global float * restrict d_etaz, const __global float * restrict d_gammax, const __global float * restrict d_gammay, const __global float * restrict d_gammaz, const __global float * restrict d_hprime_xx, const __global float * restrict d_hprimewgll_xx, const __global float * restrict d_wgllwgll_xy, const __global float * restrict d_wgllwgll_xz, const __global float * restrict d_wgllwgll_yz, const __global float * restrict d_muvstore, const int COMPUTE_AND_STORE_STRAIN, __global float * epsilondev_xx, __global float * epsilondev_yy, __global float * epsilondev_xy, __global float * epsilondev_xz, __global float * epsilondev_yz, __global float * epsilon_trace_over_3, const int ATTENUATION, const int PARTIAL_PHYS_DISPERSION_ONLY, const int USE_3D_ATTENUATION_ARRAYS, const __global float * restrict one_minus_sum_beta, const __global float * restrict factor_common, __global float * R_xx, __global float * R_yy, __global float * R_xy, __global float * R_xz, __global float * R_yz, __global float * R_xx_lddrk, __global float * R_yy_lddrk, __global float * R_xy_lddrk, __global float * R_xz_lddrk, __global float * R_yz_lddrk, const float alpha_lddrk, const float beta_lddrk, const int use_lddrk, const __global float * restrict alphaval, const __global float * restrict betaval, const __global float * restrict gammaval, const __global float * restrict tau_sigmainvval, const __global float * restrict d_c11store, const __global float * restrict d_c12store, const __global float * restrict d_c13store, const __global float * restrict d_c14store, const __global float * restrict d_c15store, const __global float * restrict d_c16store, const __global float * restrict d_c22store, const __global float * restrict d_c23store, const __global float * restrict d_c24store, const __global float * restrict d_c25store, const __global float * restrict d_c26store, const __global float * restrict d_c33store, const __global float * restrict d_c34store, const __global float * restrict d_c35store, const __global float * restrict d_c36store, const __global float * restrict d_c44store, const __global float * restrict d_c45store, const __global float * restrict d_c46store, const __global float * restrict d_c55store, const __global float * restrict d_c56store, const __global float * restrict d_c66store, const int GRAVITY, const __global float * restrict d_gravity_pre_store, const __global float * restrict d_gravity_H, const __global float * restrict wgll_cube, const int NSPEC_CRUST_MANTLE_STRAIN_ONLY){\n\
#endif\n\
#endif\n\
#ifdef USE_TEXTURES_FIELDS\n\
  const sampler_t sampler_d_displ_cm_tex = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n\
  const sampler_t sampler_d_accel_cm_tex = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n\
#endif\n\
#ifdef USE_TEXTURES_CONSTANTS\n\
  const sampler_t sampler_d_hprime_xx_tex = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n\
  const sampler_t sampler_d_hprimewgll_xx_tex = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n\
#endif\n\
  int K;\n\
  int J;\n\
  int I;\n\
  ushort active_1;\n\
  int offset;\n\
  int iglob_1;\n\
  int working_element;\n\
  float fac1;\n\
  float fac2;\n\
  float fac3;\n\
  float sigma_xx;\n\
  float sigma_xy;\n\
  float sigma_xz;\n\
  float sigma_yx;\n\
  float sigma_yy;\n\
  float sigma_yz;\n\
  float sigma_zx;\n\
  float sigma_zy;\n\
  float sigma_zz;\n\
  float epsilondev_xx_loc_1;\n\
  float epsilondev_yy_loc_1;\n\
  float epsilondev_xy_loc_1;\n\
  float epsilondev_xz_loc_1;\n\
  float epsilondev_yz_loc_1;\n\
  float sum_terms1;\n\
  float sum_terms2;\n\
  float sum_terms3;\n\
  float rho_s_H_1_1;\n\
  float rho_s_H_1_2;\n\
  float rho_s_H_1_3;\n\
\n\
  // shared arrays\n\
  __local float s_dummyx_loc[(NGLL3)];\n\
  __local float s_dummyy_loc[(NGLL3)];\n\
  __local float s_dummyz_loc[(NGLL3)];\n\
\n\
  __local float s_tempx1[(NGLL3)];\n\
  __local float s_tempx2[(NGLL3)];\n\
  __local float s_tempx3[(NGLL3)];\n\
  __local float s_tempy1[(NGLL3)];\n\
  __local float s_tempy2[(NGLL3)];\n\
  __local float s_tempy3[(NGLL3)];\n\
  __local float s_tempz1[(NGLL3)];\n\
  __local float s_tempz2[(NGLL3)];\n\
  __local float s_tempz3[(NGLL3)];\n\
  __local float sh_hprime_xx[(NGLL2)];\n\
  __local float sh_hprimewgll_xx[(NGLL2)];\n\
\n\
  const int bx = (get_group_id(1)) * (get_num_groups(0)) + get_group_id(0);\n\
  if (bx >= nb_blocks_to_compute) {\n\
     return ;\n\
  }\n\
\n\
  const int tx = get_local_id(0);\n\
\n\
  active_1 = (tx < NGLL3 ? 1 : 0);\n\
\n\
  if (active_1) {\n\
#ifdef USE_MESH_COLORING_GPU\n\
    working_element = bx;\n\
#else\n\
    if (use_mesh_coloring_gpu) {\n\
      working_element = bx;\n\
    } else {\n\
      working_element = d_phase_ispec_inner[bx + (num_phase_ispec) * (d_iphase - (1))] - (1);\n\
    }\n\
#endif\n\
    iglob_1 = d_ibool[(working_element) * (NGLL3) + tx] - (1);\n\
#ifdef USE_TEXTURES_FIELDS\n\
    s_dummyx_loc[tx] = as_float(read_imageui(d_displ_cm_tex, sampler_d_displ_cm_tex, int2((iglob_1) * (3) + 0,0)).x);\n\
    s_dummyy_loc[tx] = as_float(read_imageui(d_displ_cm_tex, sampler_d_displ_cm_tex, int2((iglob_1) * (3) + 1,0)).x);\n\
    s_dummyz_loc[tx] = as_float(read_imageui(d_displ_cm_tex, sampler_d_displ_cm_tex, int2((iglob_1) * (3) + 2,0)).x);\n\
#else\n\
    s_dummyx_loc[tx] = d_displ[0 + (3) * (iglob_1)];\n\
    s_dummyy_loc[tx] = d_displ[1 + (3) * (iglob_1)];\n\
    s_dummyz_loc[tx] = d_displ[2 + (3) * (iglob_1)];\n\
#endif\n\
  }\n\
\n\
  if (tx < NGLL2) {\n\
#ifdef USE_TEXTURES_CONSTANTS\n\
    sh_hprime_xx[tx] = as_float(read_imageui(d_hprime_xx_tex, sampler_d_hprime_xx_tex, int2(tx,0)).x);\n\
    sh_hprimewgll_xx[tx] = as_float(read_imageui(d_hprimewgll_xx_tex, sampler_d_hprimewgll_xx_tex, int2(tx,0)).x);\n\
#else\n\
    sh_hprime_xx[tx] = d_hprime_xx[tx];\n\
    sh_hprimewgll_xx[tx] = d_hprimewgll_xx[tx];\n\
#endif\n\
  }\n\
\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
  K = (tx) / (NGLL2);\n\
  J = (tx - ((K) * (NGLL2))) / (NGLLX);\n\
  I = tx - ((K) * (NGLL2)) - ((J) * (NGLLX));\n\
\n\
  if (active_1) {\n\
    float tempx1l;\n\
    float tempx2l;\n\
    float tempx3l;\n\
    float tempy1l;\n\
    float tempy2l;\n\
    float tempy3l;\n\
    float tempz1l;\n\
    float tempz2l;\n\
    float tempz3l;\n\
    tempx1l = 0.0f;\n\
    tempx2l = 0.0f;\n\
    tempx3l = 0.0f;\n\
    tempy1l = 0.0f;\n\
    tempy2l = 0.0f;\n\
    tempy3l = 0.0f;\n\
    tempz1l = 0.0f;\n\
    tempz2l = 0.0f;\n\
    tempz3l = 0.0f;\n\
#ifdef MANUALLY_UNROLLED_LOOPS\n\
    fac1 = sh_hprime_xx[(0) * (NGLLX) + I];\n\
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 0]) * (fac1);\n\
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 0]) * (fac1);\n\
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 0]) * (fac1);\n\
    fac2 = sh_hprime_xx[(0) * (NGLLX) + J];\n\
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (0) * (NGLLX) + I]) * (fac2);\n\
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (0) * (NGLLX) + I]) * (fac2);\n\
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (0) * (NGLLX) + I]) * (fac2);\n\
    fac3 = sh_hprime_xx[(0) * (NGLLX) + K];\n\
    tempx3l = tempx3l + (s_dummyx_loc[(0) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempy3l = tempy3l + (s_dummyy_loc[(0) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempz3l = tempz3l + (s_dummyz_loc[(0) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    fac1 = sh_hprime_xx[(1) * (NGLLX) + I];\n\
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 1]) * (fac1);\n\
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 1]) * (fac1);\n\
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 1]) * (fac1);\n\
    fac2 = sh_hprime_xx[(1) * (NGLLX) + J];\n\
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (1) * (NGLLX) + I]) * (fac2);\n\
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (1) * (NGLLX) + I]) * (fac2);\n\
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (1) * (NGLLX) + I]) * (fac2);\n\
    fac3 = sh_hprime_xx[(1) * (NGLLX) + K];\n\
    tempx3l = tempx3l + (s_dummyx_loc[(1) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempy3l = tempy3l + (s_dummyy_loc[(1) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempz3l = tempz3l + (s_dummyz_loc[(1) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    fac1 = sh_hprime_xx[(2) * (NGLLX) + I];\n\
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 2]) * (fac1);\n\
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 2]) * (fac1);\n\
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 2]) * (fac1);\n\
    fac2 = sh_hprime_xx[(2) * (NGLLX) + J];\n\
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (2) * (NGLLX) + I]) * (fac2);\n\
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (2) * (NGLLX) + I]) * (fac2);\n\
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (2) * (NGLLX) + I]) * (fac2);\n\
    fac3 = sh_hprime_xx[(2) * (NGLLX) + K];\n\
    tempx3l = tempx3l + (s_dummyx_loc[(2) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempy3l = tempy3l + (s_dummyy_loc[(2) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempz3l = tempz3l + (s_dummyz_loc[(2) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    fac1 = sh_hprime_xx[(3) * (NGLLX) + I];\n\
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 3]) * (fac1);\n\
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 3]) * (fac1);\n\
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 3]) * (fac1);\n\
    fac2 = sh_hprime_xx[(3) * (NGLLX) + J];\n\
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (3) * (NGLLX) + I]) * (fac2);\n\
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (3) * (NGLLX) + I]) * (fac2);\n\
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (3) * (NGLLX) + I]) * (fac2);\n\
    fac3 = sh_hprime_xx[(3) * (NGLLX) + K];\n\
    tempx3l = tempx3l + (s_dummyx_loc[(3) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempy3l = tempy3l + (s_dummyy_loc[(3) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempz3l = tempz3l + (s_dummyz_loc[(3) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    fac1 = sh_hprime_xx[(4) * (NGLLX) + I];\n\
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 4]) * (fac1);\n\
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 4]) * (fac1);\n\
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 4]) * (fac1);\n\
    fac2 = sh_hprime_xx[(4) * (NGLLX) + J];\n\
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (4) * (NGLLX) + I]) * (fac2);\n\
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (4) * (NGLLX) + I]) * (fac2);\n\
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (4) * (NGLLX) + I]) * (fac2);\n\
    fac3 = sh_hprime_xx[(4) * (NGLLX) + K];\n\
    tempx3l = tempx3l + (s_dummyx_loc[(4) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempy3l = tempy3l + (s_dummyy_loc[(4) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempz3l = tempz3l + (s_dummyz_loc[(4) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
#else\n\
    for (int l = 0; l <= NGLLX - (1); l += 1) {\n\
      fac1 = sh_hprime_xx[(l) * (NGLLX) + I];\n\
      tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + l]) * (fac1);\n\
      tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + l]) * (fac1);\n\
      tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + l]) * (fac1);\n\
      fac2 = sh_hprime_xx[(l) * (NGLLX) + J];\n\
      tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (l) * (NGLLX) + I]) * (fac2);\n\
      tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (l) * (NGLLX) + I]) * (fac2);\n\
      tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (l) * (NGLLX) + I]) * (fac2);\n\
      fac3 = sh_hprime_xx[(l) * (NGLLX) + K];\n\
      tempx3l = tempx3l + (s_dummyx_loc[(l) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
      tempy3l = tempy3l + (s_dummyy_loc[(l) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
      tempz3l = tempz3l + (s_dummyz_loc[(l) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    }\n\
#endif\n\
\n\
    float xixl;\n\
    float xiyl;\n\
    float xizl;\n\
    float etaxl;\n\
    float etayl;\n\
    float etazl;\n\
    float gammaxl;\n\
    float gammayl;\n\
    float gammazl;\n\
    offset = (working_element) * (NGLL3_PADDED) + tx;\n\
    xixl = d_xix[offset];\n\
    etaxl = d_etax[offset];\n\
    gammaxl = d_gammax[offset];\n\
    xiyl = d_xiy[offset];\n\
    etayl = d_etay[offset];\n\
    gammayl = d_gammay[offset];\n\
    xizl = d_xiz[offset];\n\
    etazl = d_etaz[offset];\n\
    gammazl = d_gammaz[offset];\n\
\n\
    float duxdxl;\n\
    float duxdyl;\n\
    float duxdzl;\n\
    float duydxl;\n\
    float duydyl;\n\
    float duydzl;\n\
    float duzdxl;\n\
    float duzdyl;\n\
    float duzdzl;\n\
    duxdxl = (xixl) * (tempx1l) + (etaxl) * (tempx2l) + (gammaxl) * (tempx3l);\n\
    duxdyl = (xiyl) * (tempx1l) + (etayl) * (tempx2l) + (gammayl) * (tempx3l);\n\
    duxdzl = (xizl) * (tempx1l) + (etazl) * (tempx2l) + (gammazl) * (tempx3l);\n\
    duydxl = (xixl) * (tempy1l) + (etaxl) * (tempy2l) + (gammaxl) * (tempy3l);\n\
    duydyl = (xiyl) * (tempy1l) + (etayl) * (tempy2l) + (gammayl) * (tempy3l);\n\
    duydzl = (xizl) * (tempy1l) + (etazl) * (tempy2l) + (gammazl) * (tempy3l);\n\
    duzdxl = (xixl) * (tempz1l) + (etaxl) * (tempz2l) + (gammaxl) * (tempz3l);\n\
    duzdyl = (xiyl) * (tempz1l) + (etayl) * (tempz2l) + (gammayl) * (tempz3l);\n\
    duzdzl = (xizl) * (tempz1l) + (etazl) * (tempz2l) + (gammazl) * (tempz3l);\n\
\n\
    float duxdyl_plus_duydxl;\n\
    float duzdxl_plus_duxdzl;\n\
    float duzdyl_plus_duydzl;\n\
    duxdyl_plus_duydxl = duxdyl + duydxl;\n\
    duzdxl_plus_duxdzl = duzdxl + duxdzl;\n\
    duzdyl_plus_duydzl = duzdyl + duydzl;\n\
\n\
    if (COMPUTE_AND_STORE_STRAIN) {\n\
      float templ;\n\
      templ = (duxdxl + duydyl + duzdzl) * (0.3333333333333333f);\n\
      epsilondev_xx_loc_1 = duxdxl - (templ);\n\
      epsilondev_yy_loc_1 = duydyl - (templ);\n\
      epsilondev_xy_loc_1 = (duxdyl_plus_duydxl) * (0.5f);\n\
      epsilondev_xz_loc_1 = (duzdxl_plus_duxdzl) * (0.5f);\n\
      epsilondev_yz_loc_1 = (duzdyl_plus_duydzl) * (0.5f);\n\
      if (NSPEC_CRUST_MANTLE_STRAIN_ONLY > 1) {\n\
        epsilon_trace_over_3[tx + (working_element) * (NGLL3)] = templ;\n\
      }\n\
    }\n\
\n\
    // anisotropic elements\n\
    compute_element_cm_aniso(offset, d_c11store, d_c12store, d_c13store, d_c14store, d_c15store, d_c16store, d_c22store, d_c23store, d_c24store, d_c25store, d_c26store, d_c33store, d_c34store, d_c35store, d_c36store, d_c44store, d_c45store, d_c46store, d_c55store, d_c56store, d_c66store, duxdxl, duxdyl, duxdzl, duydxl, duydyl, duydzl, duzdxl, duzdyl, duzdzl, duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl,  &sigma_xx,  &sigma_yy,  &sigma_zz,  &sigma_xy,  &sigma_xz,  &sigma_yz);\n\
\n\
    if (ATTENUATION &&  !(PARTIAL_PHYS_DISPERSION_ONLY)) {\n\
      compute_element_cm_att_stress(tx, working_element, R_xx, R_yy, R_xy, R_xz, R_yz,  &sigma_xx,  &sigma_yy,  &sigma_zz,  &sigma_xy,  &sigma_xz,  &sigma_yz);\n\
    }\n\
\n\
    sigma_yx = sigma_xy;\n\
    sigma_zx = sigma_xz;\n\
    sigma_zy = sigma_yz;\n\
\n\
    float jacobianl;\n\
    jacobianl = (1.0f) / ((xixl) * ((etayl) * (gammazl) - ((etazl) * (gammayl))) - ((xiyl) * ((etaxl) * (gammazl) - ((etazl) * (gammaxl)))) + (xizl) * ((etaxl) * (gammayl) - ((etayl) * (gammaxl))));\n\
\n\
    if (GRAVITY) {\n\
      compute_element_cm_gravity(tx, iglob_1, d_gravity_pre_store, d_gravity_H, wgll_cube, jacobianl, s_dummyx_loc, s_dummyy_loc, s_dummyz_loc,  &sigma_xx,  &sigma_yy,  &sigma_zz,  &sigma_xy,  &sigma_yx,  &sigma_xz,  &sigma_zx,  &sigma_yz,  &sigma_zy,  &rho_s_H_1_1,  &rho_s_H_1_2,  &rho_s_H_1_3);\n\
    }\n\
\n\
    s_tempx1[tx] = (jacobianl) * ((sigma_xx) * (xixl) + (sigma_yx) * (xiyl) + (sigma_zx) * (xizl));\n\
    s_tempy1[tx] = (jacobianl) * ((sigma_xy) * (xixl) + (sigma_yy) * (xiyl) + (sigma_zy) * (xizl));\n\
    s_tempz1[tx] = (jacobianl) * ((sigma_xz) * (xixl) + (sigma_yz) * (xiyl) + (sigma_zz) * (xizl));\n\
    s_tempx2[tx] = (jacobianl) * ((sigma_xx) * (etaxl) + (sigma_yx) * (etayl) + (sigma_zx) * (etazl));\n\
    s_tempy2[tx] = (jacobianl) * ((sigma_xy) * (etaxl) + (sigma_yy) * (etayl) + (sigma_zy) * (etazl));\n\
    s_tempz2[tx] = (jacobianl) * ((sigma_xz) * (etaxl) + (sigma_yz) * (etayl) + (sigma_zz) * (etazl));\n\
    s_tempx3[tx] = (jacobianl) * ((sigma_xx) * (gammaxl) + (sigma_yx) * (gammayl) + (sigma_zx) * (gammazl));\n\
    s_tempy3[tx] = (jacobianl) * ((sigma_xy) * (gammaxl) + (sigma_yy) * (gammayl) + (sigma_zy) * (gammazl));\n\
    s_tempz3[tx] = (jacobianl) * ((sigma_xz) * (gammaxl) + (sigma_yz) * (gammayl) + (sigma_zz) * (gammazl));\n\
  }\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
  if (active_1) {\n\
    float tempx1l;\n\
    float tempx2l;\n\
    float tempx3l;\n\
    float tempy1l;\n\
    float tempy2l;\n\
    float tempy3l;\n\
    float tempz1l;\n\
    float tempz2l;\n\
    float tempz3l;\n\
    tempx1l = 0.0f;\n\
    tempx2l = 0.0f;\n\
    tempx3l = 0.0f;\n\
    tempy1l = 0.0f;\n\
    tempy2l = 0.0f;\n\
    tempy3l = 0.0f;\n\
    tempz1l = 0.0f;\n\
    tempz2l = 0.0f;\n\
    tempz3l = 0.0f;\n\
#ifdef MANUALLY_UNROLLED_LOOPS\n\
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 0];\n\
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 0;\n\
    tempx1l = tempx1l + (s_tempx1[offset]) * (fac1);\n\
    tempy1l = tempy1l + (s_tempy1[offset]) * (fac1);\n\
    tempz1l = tempz1l + (s_tempz1[offset]) * (fac1);\n\
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 0];\n\
    offset = (K) * (NGLL2) + (0) * (NGLLX) + I;\n\
    tempx2l = tempx2l + (s_tempx2[offset]) * (fac2);\n\
    tempy2l = tempy2l + (s_tempy2[offset]) * (fac2);\n\
    tempz2l = tempz2l + (s_tempz2[offset]) * (fac2);\n\
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 0];\n\
    offset = (0) * (NGLL2) + (J) * (NGLLX) + I;\n\
    tempx3l = tempx3l + (s_tempx3[offset]) * (fac3);\n\
    tempy3l = tempy3l + (s_tempy3[offset]) * (fac3);\n\
    tempz3l = tempz3l + (s_tempz3[offset]) * (fac3);\n\
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 1];\n\
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 1;\n\
    tempx1l = tempx1l + (s_tempx1[offset]) * (fac1);\n\
    tempy1l = tempy1l + (s_tempy1[offset]) * (fac1);\n\
    tempz1l = tempz1l + (s_tempz1[offset]) * (fac1);\n\
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 1];\n\
    offset = (K) * (NGLL2) + (1) * (NGLLX) + I;\n\
    tempx2l = tempx2l + (s_tempx2[offset]) * (fac2);\n\
    tempy2l = tempy2l + (s_tempy2[offset]) * (fac2);\n\
    tempz2l = tempz2l + (s_tempz2[offset]) * (fac2);\n\
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 1];\n\
    offset = (1) * (NGLL2) + (J) * (NGLLX) + I;\n\
    tempx3l = tempx3l + (s_tempx3[offset]) * (fac3);\n\
    tempy3l = tempy3l + (s_tempy3[offset]) * (fac3);\n\
    tempz3l = tempz3l + (s_tempz3[offset]) * (fac3);\n\
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 2];\n\
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 2;\n\
    tempx1l = tempx1l + (s_tempx1[offset]) * (fac1);\n\
    tempy1l = tempy1l + (s_tempy1[offset]) * (fac1);\n\
    tempz1l = tempz1l + (s_tempz1[offset]) * (fac1);\n\
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 2];\n\
    offset = (K) * (NGLL2) + (2) * (NGLLX) + I;\n\
    tempx2l = tempx2l + (s_tempx2[offset]) * (fac2);\n\
    tempy2l = tempy2l + (s_tempy2[offset]) * (fac2);\n\
    tempz2l = tempz2l + (s_tempz2[offset]) * (fac2);\n\
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 2];\n\
    offset = (2) * (NGLL2) + (J) * (NGLLX) + I;\n\
    tempx3l = tempx3l + (s_tempx3[offset]) * (fac3);\n\
    tempy3l = tempy3l + (s_tempy3[offset]) * (fac3);\n\
    tempz3l = tempz3l + (s_tempz3[offset]) * (fac3);\n\
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 3];\n\
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 3;\n\
    tempx1l = tempx1l + (s_tempx1[offset]) * (fac1);\n\
    tempy1l = tempy1l + (s_tempy1[offset]) * (fac1);\n\
    tempz1l = tempz1l + (s_tempz1[offset]) * (fac1);\n\
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 3];\n\
    offset = (K) * (NGLL2) + (3) * (NGLLX) + I;\n\
    tempx2l = tempx2l + (s_tempx2[offset]) * (fac2);\n\
    tempy2l = tempy2l + (s_tempy2[offset]) * (fac2);\n\
    tempz2l = tempz2l + (s_tempz2[offset]) * (fac2);\n\
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 3];\n\
    offset = (3) * (NGLL2) + (J) * (NGLLX) + I;\n\
    tempx3l = tempx3l + (s_tempx3[offset]) * (fac3);\n\
    tempy3l = tempy3l + (s_tempy3[offset]) * (fac3);\n\
    tempz3l = tempz3l + (s_tempz3[offset]) * (fac3);\n\
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 4];\n\
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 4;\n\
    tempx1l = tempx1l + (s_tempx1[offset]) * (fac1);\n\
    tempy1l = tempy1l + (s_tempy1[offset]) * (fac1);\n\
    tempz1l = tempz1l + (s_tempz1[offset]) * (fac1);\n\
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 4];\n\
    offset = (K) * (NGLL2) + (4) * (NGLLX) + I;\n\
    tempx2l = tempx2l + (s_tempx2[offset]) * (fac2);\n\
    tempy2l = tempy2l + (s_tempy2[offset]) * (fac2);\n\
    tempz2l = tempz2l + (s_tempz2[offset]) * (fac2);\n\
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 4];\n\
    offset = (4) * (NGLL2) + (J) * (NGLLX) + I;\n\
    tempx3l = tempx3l + (s_tempx3[offset]) * (fac3);\n\
    tempy3l = tempy3l + (s_tempy3[offset]) * (fac3);\n\
    tempz3l = tempz3l + (s_tempz3[offset]) * (fac3);\n\
#else\n\
    for (int l = 0; l <= NGLLX - (1); l += 1) {\n\
      fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + l];\n\
      offset = (K) * (NGLL2) + (J) * (NGLLX) + l;\n\
      tempx1l = tempx1l + (s_tempx1[offset]) * (fac1);\n\
      tempy1l = tempy1l + (s_tempy1[offset]) * (fac1);\n\
      tempz1l = tempz1l + (s_tempz1[offset]) * (fac1);\n\
      fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + l];\n\
      offset = (K) * (NGLL2) + (l) * (NGLLX) + I;\n\
      tempx2l = tempx2l + (s_tempx2[offset]) * (fac2);\n\
      tempy2l = tempy2l + (s_tempy2[offset]) * (fac2);\n\
      tempz2l = tempz2l + (s_tempz2[offset]) * (fac2);\n\
      fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + l];\n\
      offset = (l) * (NGLL2) + (J) * (NGLLX) + I;\n\
      tempx3l = tempx3l + (s_tempx3[offset]) * (fac3);\n\
      tempy3l = tempy3l + (s_tempy3[offset]) * (fac3);\n\
      tempz3l = tempz3l + (s_tempz3[offset]) * (fac3);\n\
    }\n\
#endif\n\
\n\
    fac1 = d_wgllwgll_yz[(K) * (NGLLX) + J];\n\
    fac2 = d_wgllwgll_xz[(K) * (NGLLX) + I];\n\
    fac3 = d_wgllwgll_xy[(J) * (NGLLX) + I];\n\
    sum_terms1 =  -((fac1) * (tempx1l) + (fac2) * (tempx2l) + (fac3) * (tempx3l));\n\
    sum_terms2 =  -((fac1) * (tempy1l) + (fac2) * (tempy2l) + (fac3) * (tempy3l));\n\
    sum_terms3 =  -((fac1) * (tempz1l) + (fac2) * (tempz2l) + (fac3) * (tempz3l));\n\
\n\
    if (GRAVITY) {\n\
      sum_terms1 = sum_terms1 + rho_s_H_1_1;\n\
      sum_terms2 = sum_terms2 + rho_s_H_1_2;\n\
      sum_terms3 = sum_terms3 + rho_s_H_1_3;\n\
    }\n\
\n\
#ifdef USE_MESH_COLORING_GPU\n\
#ifdef USE_TEXTURES_FIELDS\n\
    d_accel[0 + (3) * (iglob_1)] = as_float(read_imageui(d_accel_cm_tex, sampler_d_accel_cm_tex, int2((iglob_1) * (3) + 0,0)).x) + sum_terms1;\n\
    d_accel[1 + (3) * (iglob_1)] = as_float(read_imageui(d_accel_cm_tex, sampler_d_accel_cm_tex, int2((iglob_1) * (3) + 1,0)).x) + sum_terms2;\n\
    d_accel[2 + (3) * (iglob_1)] = as_float(read_imageui(d_accel_cm_tex, sampler_d_accel_cm_tex, int2((iglob_1) * (3) + 2,0)).x) + sum_terms3;\n\
#else\n\
    d_accel[0 + (3) * (iglob_1)] = d_accel[0 + (3) * (iglob_1)] + sum_terms1;\n\
    d_accel[1 + (3) * (iglob_1)] = d_accel[1 + (3) * (iglob_1)] + sum_terms2;\n\
    d_accel[2 + (3) * (iglob_1)] = d_accel[2 + (3) * (iglob_1)] + sum_terms3;\n\
#endif\n\
#else\n\
    if (use_mesh_coloring_gpu) {\n\
#ifdef USE_TEXTURES_FIELDS\n\
      d_accel[0 + (3) * (iglob_1)] = as_float(read_imageui(d_accel_cm_tex, sampler_d_accel_cm_tex, int2((iglob_1) * (3) + 0,0)).x) + sum_terms1;\n\
      d_accel[1 + (3) * (iglob_1)] = as_float(read_imageui(d_accel_cm_tex, sampler_d_accel_cm_tex, int2((iglob_1) * (3) + 1,0)).x) + sum_terms2;\n\
      d_accel[2 + (3) * (iglob_1)] = as_float(read_imageui(d_accel_cm_tex, sampler_d_accel_cm_tex, int2((iglob_1) * (3) + 2,0)).x) + sum_terms3;\n\
#else\n\
      d_accel[0 + (3) * (iglob_1)] = d_accel[0 + (3) * (iglob_1)] + sum_terms1;\n\
      d_accel[1 + (3) * (iglob_1)] = d_accel[1 + (3) * (iglob_1)] + sum_terms2;\n\
      d_accel[2 + (3) * (iglob_1)] = d_accel[2 + (3) * (iglob_1)] + sum_terms3;\n\
#endif\n\
    } else {\n\
      atomicAdd(d_accel + (iglob_1) * (3) + 0, sum_terms1);\n\
      atomicAdd(d_accel + (iglob_1) * (3) + 1, sum_terms2);\n\
      atomicAdd(d_accel + (iglob_1) * (3) + 2, sum_terms3);\n\
    }\n\
#endif\n\
\n\
    if (ATTENUATION &&  !(PARTIAL_PHYS_DISPERSION_ONLY)) {\n\
      if ( !(use_lddrk)) {\n\
        compute_element_cm_att_memory(tx, working_element, d_muvstore, factor_common, alphaval, betaval, gammaval, R_xx, R_yy, R_xy, R_xz, R_yz, epsilondev_xx, epsilondev_yy, epsilondev_xy, epsilondev_xz, epsilondev_yz, epsilondev_xx_loc_1, epsilondev_yy_loc_1, epsilondev_xy_loc_1, epsilondev_xz_loc_1, epsilondev_yz_loc_1, USE_3D_ATTENUATION_ARRAYS);\n\
      } else {\n\
        compute_element_cm_att_memory_lddrk(tx, working_element, d_muvstore, factor_common, tau_sigmainvval, R_xx, R_yy, R_xy, R_xz, R_yz, R_xx_lddrk, R_yy_lddrk, R_xy_lddrk, R_xz_lddrk, R_yz_lddrk, alpha_lddrk, beta_lddrk, deltat, epsilondev_xx_loc_1, epsilondev_yy_loc_1, epsilondev_xy_loc_1, epsilondev_xz_loc_1, epsilondev_yz_loc_1, USE_3D_ATTENUATION_ARRAYS);\n\
      }\n\
    }\n\
\n\
    if (COMPUTE_AND_STORE_STRAIN) {\n\
      epsilondev_xx[tx + (working_element) * (NGLL3)] = epsilondev_xx_loc_1;\n\
      epsilondev_yy[tx + (working_element) * (NGLL3)] = epsilondev_yy_loc_1;\n\
      epsilondev_xy[tx + (working_element) * (NGLL3)] = epsilondev_xy_loc_1;\n\
      epsilondev_xz[tx + (working_element) * (NGLL3)] = epsilondev_xz_loc_1;\n\
      epsilondev_yz[tx + (working_element) * (NGLL3)] = epsilondev_yz_loc_1;\n\
    }\n\
  }\n\
}\n\
";

