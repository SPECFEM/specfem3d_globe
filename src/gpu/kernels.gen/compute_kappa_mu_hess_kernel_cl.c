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

const char * compute_kappa_mu_hess_kernel_program = "\
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
#if __OPENCL_C_VERSION__ && __OPENCL_C_VERSION__ >= 120\n\
static\n\
#endif\n\
void compute_vector_gradient_kernel(const int ispec, const __local float * fx, const __local float * fy, const __local float * fz, const __global float * d_xix, const __global float * d_xiy, const __global float * d_xiz, const __global float * d_etax, const __global float * d_etay, const __global float * d_etaz, const __global float * d_gammax, const __global float * d_gammay, const __global float * d_gammaz, const __local float * sh_hprime_xx, float * fgrad_loc){\n\
  int tx;\n\
  int K;\n\
  int J;\n\
  int I;\n\
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
  float dfxdxl;\n\
  float dfxdyl;\n\
  float dfxdzl;\n\
  float dfydxl;\n\
  float dfydyl;\n\
  float dfydzl;\n\
  float dfzdxl;\n\
  float dfzdyl;\n\
  float dfzdzl;\n\
  float fac1;\n\
  float fac2;\n\
  float fac3;\n\
\n\
  tx = get_local_id(0);\n\
  K = (tx) / (NGLL2);\n\
  J = (tx - ((K) * (NGLL2))) / (NGLLX);\n\
  I = tx - ((K) * (NGLL2)) - ((J) * (NGLLX));\n\
\n\
  tempx1l = 0.0f;\n\
  tempx2l = 0.0f;\n\
  tempx3l = 0.0f;\n\
  tempy1l = 0.0f;\n\
  tempy2l = 0.0f;\n\
  tempy3l = 0.0f;\n\
  tempz1l = 0.0f;\n\
  tempz2l = 0.0f;\n\
  tempz3l = 0.0f;\n\
\n\
  for (int l = 0; l < NGLLX; l += 1) {\n\
    fac1 = sh_hprime_xx[(l) * (NGLLX) + I];\n\
    tempx1l = tempx1l + (fx[(K) * (NGLL2) + (J) * (NGLLX) + l]) * (fac1);\n\
    tempy1l = tempy1l + (fy[(K) * (NGLL2) + (J) * (NGLLX) + l]) * (fac1);\n\
    tempz1l = tempz1l + (fz[(K) * (NGLL2) + (J) * (NGLLX) + l]) * (fac1);\n\
    fac2 = sh_hprime_xx[(l) * (NGLLX) + J];\n\
    tempx2l = tempx2l + (fx[(K) * (NGLL2) + (l) * (NGLLX) + I]) * (fac2);\n\
    tempy2l = tempy2l + (fy[(K) * (NGLL2) + (l) * (NGLLX) + I]) * (fac2);\n\
    tempz2l = tempz2l + (fz[(K) * (NGLL2) + (l) * (NGLLX) + I]) * (fac2);\n\
    fac3 = sh_hprime_xx[(l) * (NGLLX) + K];\n\
    tempx3l = tempx3l + (fx[(l) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempy3l = tempy3l + (fy[(l) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
    tempz3l = tempz3l + (fz[(l) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);\n\
  }\n\
\n\
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
\n\
  dfxdxl = (xixl) * (tempx1l) + (etaxl) * (tempx2l) + (gammaxl) * (tempx3l);\n\
  dfxdyl = (xiyl) * (tempx1l) + (etayl) * (tempx2l) + (gammayl) * (tempx3l);\n\
  dfxdzl = (xizl) * (tempx1l) + (etazl) * (tempx2l) + (gammazl) * (tempx3l);\n\
  dfydxl = (xixl) * (tempy1l) + (etaxl) * (tempy2l) + (gammaxl) * (tempy3l);\n\
  dfydyl = (xiyl) * (tempy1l) + (etayl) * (tempy2l) + (gammayl) * (tempy3l);\n\
  dfydzl = (xizl) * (tempy1l) + (etazl) * (tempy2l) + (gammazl) * (tempy3l);\n\
  dfzdxl = (xixl) * (tempz1l) + (etaxl) * (tempz2l) + (gammaxl) * (tempz3l);\n\
  dfzdyl = (xiyl) * (tempz1l) + (etayl) * (tempz2l) + (gammayl) * (tempz3l);\n\
  dfzdzl = (xizl) * (tempz1l) + (etazl) * (tempz2l) + (gammazl) * (tempz3l);\n\
\n\
  fgrad_loc[0] = dfxdxl;\n\
  fgrad_loc[1] = dfxdyl;\n\
  fgrad_loc[2] = dfxdzl;\n\
  fgrad_loc[3] = dfydxl;\n\
  fgrad_loc[4] = dfydyl;\n\
  fgrad_loc[5] = dfydzl;\n\
  fgrad_loc[6] = dfzdxl;\n\
  fgrad_loc[7] = dfzdyl;\n\
  fgrad_loc[8] = dfzdzl;\n\
}\n\
__kernel void compute_kappa_mu_hess_kernel(const __global int * d_ibool, const __global float * d_veloc, const __global float * d_b_veloc, const __global float * d_xix, const __global float * d_xiy, const __global float * d_xiz, const __global float * d_etax, const __global float * d_etay, const __global float * d_etaz, const __global float * d_gammax, const __global float * d_gammay, const __global float * d_gammaz, const __global float * d_hprime_xx, const float deltat, __global float * hess_rho_kl, __global float * hess_kappa_kl, __global float * hess_mu_kl, const int NSPEC, const int USE_SOURCE_RECEIVER_HESSIAN){\n\
  int ispec;\n\
  int ijk_ispec;\n\
  int tx;\n\
  int iglob;\n\
  float vgrad[(9)];\n\
  float b_vgrad[(9)];\n\
  __local float sh_velocx[(NGLL3)];\n\
  __local float sh_velocy[(NGLL3)];\n\
  __local float sh_velocz[(NGLL3)];\n\
  __local float sh_b_velocx[(NGLL3)];\n\
  __local float sh_b_velocy[(NGLL3)];\n\
  __local float sh_b_velocz[(NGLL3)];\n\
  __local float sh_hprime_xx[(NGLL2)];\n\
  float hess_rhol;\n\
  float hess_kappal;\n\
  float hess_mul;\n\
\n\
  ispec = get_group_id(0) + (get_group_id(1)) * (get_num_groups(0));\n\
  ijk_ispec = get_local_id(0) + (NGLL3) * (ispec);\n\
  tx = get_local_id(0);\n\
\n\
  if (tx < NGLL2) {\n\
    sh_hprime_xx[tx] = d_hprime_xx[tx];\n\
  }\n\
\n\
  if (ispec < NSPEC) {\n\
    iglob = d_ibool[ijk_ispec] - (1);\n\
    sh_velocx[tx] = d_veloc[0 + (3) * (iglob)];\n\
    sh_b_velocx[tx] = d_b_veloc[0 + (3) * (iglob)];\n\
    sh_velocy[tx] = d_veloc[1 + (3) * (iglob)];\n\
    sh_b_velocy[tx] = d_b_veloc[1 + (3) * (iglob)];\n\
    sh_velocz[tx] = d_veloc[2 + (3) * (iglob)];\n\
    sh_b_velocz[tx] = d_b_veloc[2 + (3) * (iglob)];\n\
  }\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
  if (ispec < NSPEC) {\n\
    compute_vector_gradient_kernel(ispec, sh_velocx, sh_velocy, sh_velocz, d_xix, d_xiy, d_xiz, d_etax, d_etay, d_etaz, d_gammax, d_gammay, d_gammaz, sh_hprime_xx, vgrad);\n\
\n\
    if (USE_SOURCE_RECEIVER_HESSIAN) {\n\
      compute_vector_gradient_kernel(ispec, sh_b_velocx, sh_b_velocy, sh_b_velocz, d_xix, d_xiy, d_xiz, d_etax, d_etay, d_etaz, d_gammax, d_gammay, d_gammaz, sh_hprime_xx, b_vgrad);\n\
    } else {\n\
      for (int l = 0; l <= 8; l += 1) {\n\
        b_vgrad[l] = vgrad[l];\n\
      }\n\
    }\n\
\n\
    hess_rhol = (vgrad[0]) * (b_vgrad[0]) + (vgrad[4]) * (b_vgrad[4]) + (vgrad[8]) * (b_vgrad[8]);\n\
    hess_kappal = ((3.0f) * (vgrad[0] + vgrad[4] + vgrad[8])) * (b_vgrad[0] + b_vgrad[4] + b_vgrad[8]);\n\
    hess_mul = (vgrad[1] + vgrad[3]) * (b_vgrad[1] + b_vgrad[3]) + (vgrad[2] + vgrad[6]) * (b_vgrad[2] + b_vgrad[6]) + (vgrad[5] + vgrad[7]) * (b_vgrad[5] + b_vgrad[7]) + (((2.0f) * (vgrad[0]) - (vgrad[4]) - (vgrad[8])) * ((2.0f) * (b_vgrad[0]) - (b_vgrad[4]) - (b_vgrad[8])) + ( -(vgrad[0]) + (2.0f) * (vgrad[4]) - (vgrad[8])) * ( -(b_vgrad[0]) + (2.0f) * (b_vgrad[4]) - (b_vgrad[8])) + ( -(vgrad[0]) - (vgrad[4]) + (2.0f) * (vgrad[8])) * ( -(b_vgrad[0]) - (b_vgrad[4]) + (2.0f) * (b_vgrad[8]))) * (0.4444444444444444f);\n\
\n\
    hess_rho_kl[ijk_ispec] = hess_rho_kl[ijk_ispec] + (deltat) * (hess_rhol);\n\
    hess_kappa_kl[ijk_ispec] = hess_kappa_kl[ijk_ispec] + (deltat) * (hess_kappal);\n\
    hess_mu_kl[ijk_ispec] = hess_mu_kl[ijk_ispec] + (deltat) * (hess_mul);\n\
  }\n\
}\n\
";
