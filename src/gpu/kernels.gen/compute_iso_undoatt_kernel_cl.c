//note: please do not modify this file manually!
//      this file has been generated automatically by BOAST version 2.1.0
//      by: make boast_kernels

/*
!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
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

const char * compute_iso_undoatt_kernel_program = "\
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
#ifndef R_EARTH_KM\n\
#define R_EARTH_KM 6371.0f\n\
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
#if __OPENCL_C_VERSION__ && __OPENCL_C_VERSION__ >= 120\n\
static\n\
#endif\n\
void compute_element_strain_undoatt(const int ispec, const int ijk_ispec, const __global int * d_ibool, const __local float * s_dummyx_loc, const __local float * s_dummyy_loc, const __local float * s_dummyz_loc, const __global float * d_xix, const __global float * d_xiy, const __global float * d_xiz, const __global float * d_etax, const __global float * d_etay, const __global float * d_etaz, const __global float * d_gammax, const __global float * d_gammay, const __global float * d_gammaz, const __local float * sh_hprime_xx, float * epsilondev_loc, float * epsilon_trace_over_3){\n\
  int tx;\n\
  int K;\n\
  int J;\n\
  int I;\n\
  int l;\n\
  int offset;\n\
  float tempx1l;\n\
  float tempx2l;\n\
  float tempx3l;\n\
  float tempy1l;\n\
  float tempy2l;\n\
  float tempy3l;\n\
  float tempz1l;\n\
  float tempz2l;\n\
  float tempz3l;\n\
  float xixl;\n\
  float xiyl;\n\
  float xizl;\n\
  float etaxl;\n\
  float etayl;\n\
  float etazl;\n\
  float gammaxl;\n\
  float gammayl;\n\
  float gammazl;\n\
  float duxdxl;\n\
  float duxdyl;\n\
  float duxdzl;\n\
  float duydxl;\n\
  float duydyl;\n\
  float duydzl;\n\
  float duzdxl;\n\
  float duzdyl;\n\
  float duzdzl;\n\
  float templ;\n\
  float fac1;\n\
  float fac2;\n\
  float fac3;\n\
  tx = get_local_id(0);\n\
  K = (tx) / (NGLL2);\n\
  J = (tx - ((K) * (NGLL2))) / (NGLLX);\n\
  I = tx - ((K) * (NGLL2)) - ((J) * (NGLLX));\n\
  tempx1l = 0.0f;\n\
  tempx2l = 0.0f;\n\
  tempx3l = 0.0f;\n\
  tempy1l = 0.0f;\n\
  tempy2l = 0.0f;\n\
  tempy3l = 0.0f;\n\
  tempz1l = 0.0f;\n\
  tempz2l = 0.0f;\n\
  tempz3l = 0.0f;\n\
  for (l = 0; l <= NGLLX - (1); l += 1) {\n\
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
  offset = (ispec) * (NGLL3_PADDED) + tx;\n\
  xixl = d_xix[offset];\n\
  etaxl = d_etax[offset];\n\
  gammaxl = d_gammax[offset];\n\
  xiyl = d_xiy[offset];\n\
  etayl = d_etay[offset];\n\
  gammayl = d_gammay[offset];\n\
  xizl = d_xiz[offset];\n\
  etazl = d_etaz[offset];\n\
  gammazl = d_gammaz[offset];\n\
  duxdxl = (xixl) * (tempx1l) + (etaxl) * (tempx2l) + (gammaxl) * (tempx3l);\n\
  duxdyl = (xiyl) * (tempx1l) + (etayl) * (tempx2l) + (gammayl) * (tempx3l);\n\
  duxdzl = (xizl) * (tempx1l) + (etazl) * (tempx2l) + (gammazl) * (tempx3l);\n\
  duydxl = (xixl) * (tempy1l) + (etaxl) * (tempy2l) + (gammaxl) * (tempy3l);\n\
  duydyl = (xiyl) * (tempy1l) + (etayl) * (tempy2l) + (gammayl) * (tempy3l);\n\
  duydzl = (xizl) * (tempy1l) + (etazl) * (tempy2l) + (gammazl) * (tempy3l);\n\
  duzdxl = (xixl) * (tempz1l) + (etaxl) * (tempz2l) + (gammaxl) * (tempz3l);\n\
  duzdyl = (xiyl) * (tempz1l) + (etayl) * (tempz2l) + (gammayl) * (tempz3l);\n\
  duzdzl = (xizl) * (tempz1l) + (etazl) * (tempz2l) + (gammazl) * (tempz3l);\n\
  templ = (duxdxl + duydyl + duzdzl) * (0.3333333333333333f);\n\
  epsilondev_loc[0] = duxdxl - (templ);\n\
  epsilondev_loc[1] = duydyl - (templ);\n\
  epsilondev_loc[2] = (duxdyl + duydxl) * (0.5f);\n\
  epsilondev_loc[3] = (duzdxl + duxdzl) * (0.5f);\n\
  epsilondev_loc[4] = (duzdyl + duydzl) * (0.5f);\n\
  *(epsilon_trace_over_3) = templ;\n\
}\n\
__kernel void compute_iso_undoatt_kernel(const __global float * epsilondev_xx, const __global float * epsilondev_yy, const __global float * epsilondev_xy, const __global float * epsilondev_xz, const __global float * epsilondev_yz, const __global float * epsilon_trace_over_3, __global float * mu_kl, __global float * kappa_kl, const int NSPEC, const float deltat, const __global int * d_ibool, const __global float * d_b_displ, const __global float * d_xix, const __global float * d_xiy, const __global float * d_xiz, const __global float * d_etax, const __global float * d_etay, const __global float * d_etaz, const __global float * d_gammax, const __global float * d_gammay, const __global float * d_gammaz, const __global float * d_hprime_xx){\n\
  int ispec;\n\
  int ijk_ispec;\n\
  int tx;\n\
  int iglob;\n\
  float eps_trace_over_3;\n\
  float b_eps_trace_over_3;\n\
  float epsdev[(5)];\n\
  float b_epsdev[(5)];\n\
  __local float s_dummyx_loc[(NGLL3)];\n\
  __local float s_dummyy_loc[(NGLL3)];\n\
  __local float s_dummyz_loc[(NGLL3)];\n\
  __local float sh_hprime_xx[(NGLL2)];\n\
  ispec = get_group_id(0) + (get_group_id(1)) * (get_num_groups(0));\n\
  ijk_ispec = get_local_id(0) + (NGLL3) * (ispec);\n\
  tx = get_local_id(0);\n\
  if (tx < NGLL2) {\n\
    sh_hprime_xx[tx] = d_hprime_xx[tx];\n\
  }\n\
  if (ispec < NSPEC) {\n\
    iglob = d_ibool[ijk_ispec] - (1);\n\
    s_dummyx_loc[tx] = d_b_displ[0 + (3) * (iglob)];\n\
    s_dummyy_loc[tx] = d_b_displ[1 + (3) * (iglob)];\n\
    s_dummyz_loc[tx] = d_b_displ[2 + (3) * (iglob)];\n\
  }\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if (ispec < NSPEC) {\n\
    epsdev[0] = epsilondev_xx[ijk_ispec];\n\
    epsdev[1] = epsilondev_yy[ijk_ispec];\n\
    epsdev[2] = epsilondev_xy[ijk_ispec];\n\
    epsdev[3] = epsilondev_xz[ijk_ispec];\n\
    epsdev[4] = epsilondev_yz[ijk_ispec];\n\
    eps_trace_over_3 = epsilon_trace_over_3[ijk_ispec];\n\
    compute_element_strain_undoatt(ispec, ijk_ispec, d_ibool, s_dummyx_loc, s_dummyy_loc, s_dummyz_loc, d_xix, d_xiy, d_xiz, d_etax, d_etay, d_etaz, d_gammax, d_gammay, d_gammaz, sh_hprime_xx, b_epsdev,  &b_eps_trace_over_3);\n\
    mu_kl[ijk_ispec] = mu_kl[ijk_ispec] + (deltat) * ((epsdev[0]) * (b_epsdev[0]) + (epsdev[1]) * (b_epsdev[1]) + (epsdev[0] + epsdev[1]) * (b_epsdev[0] + b_epsdev[1]) + ((epsdev[2]) * (b_epsdev[2]) + (epsdev[3]) * (b_epsdev[3]) + (epsdev[4]) * (b_epsdev[4])) * (2.0f));\n\
    kappa_kl[ijk_ispec] = kappa_kl[ijk_ispec] + (deltat) * (((eps_trace_over_3) * (b_eps_trace_over_3)) * (9.0f));\n\
  }\n\
}\n\
";
