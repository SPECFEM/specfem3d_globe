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

#ifndef INDEX2
#define INDEX2(isize,i,j) i + isize*j
#endif
#ifndef INDEX3
#define INDEX3(isize,jsize,i,j,k) i + isize*(j + jsize*k)
#endif
#ifndef INDEX4
#define INDEX4(isize,jsize,ksize,i,j,k,x) i + isize*(j + jsize*(k + ksize*x))
#endif
#ifndef INDEX5
#define INDEX5(isize,jsize,ksize,xsize,i,j,k,x,y) i + isize*(j + jsize*(k + ksize*(x + xsize*y)))
#endif

#ifndef NDIM
#define NDIM 3
#endif
#ifndef NGLLX
#define NGLLX 5
#endif
#ifndef NGLL2
#define NGLL2 25
#endif
#ifndef NGLL3
#define NGLL3 125
#endif
#ifndef NGLL3_PADDED
#define NGLL3_PADDED 128
#endif
#ifndef N_SLS
#define N_SLS 3
#endif
#ifndef IREGION_CRUST_MANTLE
#define IREGION_CRUST_MANTLE 1
#endif
#ifndef IREGION_INNER_CORE
#define IREGION_INNER_CORE 3
#endif
#ifndef IFLAG_IN_FICTITIOUS_CUBE
#define IFLAG_IN_FICTITIOUS_CUBE 11
#endif
#ifndef COLORING_MIN_NSPEC_INNER_CORE
#define COLORING_MIN_NSPEC_INNER_CORE 1000
#endif
#ifndef COLORING_MIN_NSPEC_OUTER_CORE
#define COLORING_MIN_NSPEC_OUTER_CORE 1000
#endif
#ifndef BLOCKSIZE_TRANSFER
#define BLOCKSIZE_TRANSFER 256
#endif

static __device__ void compute_vector_gradient_kernel(const int ispec, const float * fx, const float * fy, const float * fz, const float * d_xix, const float * d_xiy, const float * d_xiz, const float * d_etax, const float * d_etay, const float * d_etaz, const float * d_gammax, const float * d_gammay, const float * d_gammaz, const float * sh_hprime_xx, float * fgrad_loc){
  int tx;
  int K;
  int J;
  int I;
  int offset;
  float tempx1l;
  float tempx2l;
  float tempx3l;
  float tempy1l;
  float tempy2l;
  float tempy3l;
  float tempz1l;
  float tempz2l;
  float tempz3l;
  float xixl;
  float xiyl;
  float xizl;
  float etaxl;
  float etayl;
  float etazl;
  float gammaxl;
  float gammayl;
  float gammazl;
  float dfxdxl;
  float dfxdyl;
  float dfxdzl;
  float dfydxl;
  float dfydyl;
  float dfydzl;
  float dfzdxl;
  float dfzdyl;
  float dfzdzl;
  float fac1;
  float fac2;
  float fac3;

  tx = threadIdx.x;
  K = (tx) / (NGLL2);
  J = (tx - ((K) * (NGLL2))) / (NGLLX);
  I = tx - ((K) * (NGLL2)) - ((J) * (NGLLX));

  tempx1l = 0.0f;
  tempx2l = 0.0f;
  tempx3l = 0.0f;
  tempy1l = 0.0f;
  tempy2l = 0.0f;
  tempy3l = 0.0f;
  tempz1l = 0.0f;
  tempz2l = 0.0f;
  tempz3l = 0.0f;

  for (int l = 0; l < NGLLX; l += 1) {
    fac1 = sh_hprime_xx[(l) * (NGLLX) + I];
    tempx1l = tempx1l + (fx[(K) * (NGLL2) + (J) * (NGLLX) + l]) * (fac1);
    tempy1l = tempy1l + (fy[(K) * (NGLL2) + (J) * (NGLLX) + l]) * (fac1);
    tempz1l = tempz1l + (fz[(K) * (NGLL2) + (J) * (NGLLX) + l]) * (fac1);
    fac2 = sh_hprime_xx[(l) * (NGLLX) + J];
    tempx2l = tempx2l + (fx[(K) * (NGLL2) + (l) * (NGLLX) + I]) * (fac2);
    tempy2l = tempy2l + (fy[(K) * (NGLL2) + (l) * (NGLLX) + I]) * (fac2);
    tempz2l = tempz2l + (fz[(K) * (NGLL2) + (l) * (NGLLX) + I]) * (fac2);
    fac3 = sh_hprime_xx[(l) * (NGLLX) + K];
    tempx3l = tempx3l + (fx[(l) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);
    tempy3l = tempy3l + (fy[(l) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);
    tempz3l = tempz3l + (fz[(l) * (NGLL2) + (J) * (NGLLX) + I]) * (fac3);
  }

  offset = (ispec) * (NGLL3_PADDED) + tx;
  xixl = d_xix[offset];
  etaxl = d_etax[offset];
  gammaxl = d_gammax[offset];
  xiyl = d_xiy[offset];
  etayl = d_etay[offset];
  gammayl = d_gammay[offset];
  xizl = d_xiz[offset];
  etazl = d_etaz[offset];
  gammazl = d_gammaz[offset];

  dfxdxl = (xixl) * (tempx1l) + (etaxl) * (tempx2l) + (gammaxl) * (tempx3l);
  dfxdyl = (xiyl) * (tempx1l) + (etayl) * (tempx2l) + (gammayl) * (tempx3l);
  dfxdzl = (xizl) * (tempx1l) + (etazl) * (tempx2l) + (gammazl) * (tempx3l);
  dfydxl = (xixl) * (tempy1l) + (etaxl) * (tempy2l) + (gammaxl) * (tempy3l);
  dfydyl = (xiyl) * (tempy1l) + (etayl) * (tempy2l) + (gammayl) * (tempy3l);
  dfydzl = (xizl) * (tempy1l) + (etazl) * (tempy2l) + (gammazl) * (tempy3l);
  dfzdxl = (xixl) * (tempz1l) + (etaxl) * (tempz2l) + (gammaxl) * (tempz3l);
  dfzdyl = (xiyl) * (tempz1l) + (etayl) * (tempz2l) + (gammayl) * (tempz3l);
  dfzdzl = (xizl) * (tempz1l) + (etazl) * (tempz2l) + (gammazl) * (tempz3l);

  fgrad_loc[0] = dfxdxl;
  fgrad_loc[1] = dfxdyl;
  fgrad_loc[2] = dfxdzl;
  fgrad_loc[3] = dfydxl;
  fgrad_loc[4] = dfydyl;
  fgrad_loc[5] = dfydzl;
  fgrad_loc[6] = dfzdxl;
  fgrad_loc[7] = dfzdyl;
  fgrad_loc[8] = dfzdzl;
}
__global__ void compute_kappa_mu_hess_kernel(const int * d_ibool, const float * d_veloc, const float * d_b_veloc, const float * d_xix, const float * d_xiy, const float * d_xiz, const float * d_etax, const float * d_etay, const float * d_etaz, const float * d_gammax, const float * d_gammay, const float * d_gammaz, const float * d_hprime_xx, const float deltat, float * hess_rho_kl, float * hess_kappa_kl, float * hess_mu_kl, const int NSPEC, const int USE_SOURCE_RECEIVER_HESSIAN){
  int ispec;
  int ijk_ispec;
  int tx;
  int iglob;
  float vgrad[(9)];
  float b_vgrad[(9)];
  __shared__ float sh_velocx[(NGLL3)];
  __shared__ float sh_velocy[(NGLL3)];
  __shared__ float sh_velocz[(NGLL3)];
  __shared__ float sh_b_velocx[(NGLL3)];
  __shared__ float sh_b_velocy[(NGLL3)];
  __shared__ float sh_b_velocz[(NGLL3)];
  __shared__ float sh_hprime_xx[(NGLL2)];
  float hess_rhol;
  float hess_kappal;
  float hess_mul;

  ispec = blockIdx.x + (blockIdx.y) * (gridDim.x);
  ijk_ispec = threadIdx.x + (NGLL3) * (ispec);
  tx = threadIdx.x;

  if (tx < NGLL2) {
    sh_hprime_xx[tx] = d_hprime_xx[tx];
  }

  if (ispec < NSPEC) {
    iglob = d_ibool[ijk_ispec] - (1);
    sh_velocx[tx] = d_veloc[0 + (3) * (iglob)];
    sh_b_velocx[tx] = d_b_veloc[0 + (3) * (iglob)];
    sh_velocy[tx] = d_veloc[1 + (3) * (iglob)];
    sh_b_velocy[tx] = d_b_veloc[1 + (3) * (iglob)];
    sh_velocz[tx] = d_veloc[2 + (3) * (iglob)];
    sh_b_velocz[tx] = d_b_veloc[2 + (3) * (iglob)];
  }
  __syncthreads();

  if (ispec < NSPEC) {
    compute_vector_gradient_kernel(ispec, sh_velocx, sh_velocy, sh_velocz, d_xix, d_xiy, d_xiz, d_etax, d_etay, d_etaz, d_gammax, d_gammay, d_gammaz, sh_hprime_xx, vgrad);

    if (USE_SOURCE_RECEIVER_HESSIAN) {
      compute_vector_gradient_kernel(ispec, sh_b_velocx, sh_b_velocy, sh_b_velocz, d_xix, d_xiy, d_xiz, d_etax, d_etay, d_etaz, d_gammax, d_gammay, d_gammaz, sh_hprime_xx, b_vgrad);
    } else {
      for (int l = 0; l <= 8; l += 1) {
        b_vgrad[l] = vgrad[l];
      }
    }

    hess_rhol = (vgrad[0]) * (b_vgrad[0]) + (vgrad[4]) * (b_vgrad[4]) + (vgrad[8]) * (b_vgrad[8]);
    hess_kappal = ((3.0f) * (vgrad[0] + vgrad[4] + vgrad[8])) * (b_vgrad[0] + b_vgrad[4] + b_vgrad[8]);
    hess_mul = (vgrad[1] + vgrad[3]) * (b_vgrad[1] + b_vgrad[3]) + (vgrad[2] + vgrad[6]) * (b_vgrad[2] + b_vgrad[6]) + (vgrad[5] + vgrad[7]) * (b_vgrad[5] + b_vgrad[7]) + (((2.0f) * (vgrad[0]) - (vgrad[4]) - (vgrad[8])) * ((2.0f) * (b_vgrad[0]) - (b_vgrad[4]) - (b_vgrad[8])) + ( -(vgrad[0]) + (2.0f) * (vgrad[4]) - (vgrad[8])) * ( -(b_vgrad[0]) + (2.0f) * (b_vgrad[4]) - (b_vgrad[8])) + ( -(vgrad[0]) - (vgrad[4]) + (2.0f) * (vgrad[8])) * ( -(b_vgrad[0]) - (b_vgrad[4]) + (2.0f) * (b_vgrad[8]))) * (0.4444444444444444f);

    hess_rho_kl[ijk_ispec] = hess_rho_kl[ijk_ispec] + (deltat) * (hess_rhol);
    hess_kappa_kl[ijk_ispec] = hess_kappa_kl[ijk_ispec] + (deltat) * (hess_kappal);
    hess_mu_kl[ijk_ispec] = hess_mu_kl[ijk_ispec] + (deltat) * (hess_mul);
  }
}
