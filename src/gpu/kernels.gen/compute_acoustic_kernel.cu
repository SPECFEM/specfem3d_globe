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
#ifndef R_EARTH_KM
#define R_EARTH_KM 6371.0f
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

static __device__ void compute_gradient_kernel(const int ijk, const int ispec, const float * scalar_field, float * vector_field_element, const float * hprime_xx, const float * d_xix, const float * d_xiy, const float * d_xiz, const float * d_etax, const float * d_etay, const float * d_etaz, const float * d_gammax, const float * d_gammay, const float * d_gammaz){
  float temp1l;
  float temp2l;
  float temp3l;
  float hp1;
  float hp2;
  float hp3;
  float xixl;
  float xiyl;
  float xizl;
  float etaxl;
  float etayl;
  float etazl;
  float gammaxl;
  float gammayl;
  float gammazl;
  int l;
  int offset;
  int offset1;
  int offset2;
  int offset3;
  int I;
  int J;
  int K;
  K = (ijk) / (NGLL2);
  J = (ijk - ((K) * (NGLL2))) / (NGLLX);
  I = ijk - ((K) * (NGLL2)) - ((J) * (NGLLX));
  temp1l = 0.0f;
  temp2l = 0.0f;
  temp3l = 0.0f;
  for (l = 0; l <= NGLLX - (1); l += 1) {
    hp1 = hprime_xx[(l) * (NGLLX) + I];
    hp2 = hprime_xx[(l) * (NGLLX) + J];
    hp3 = hprime_xx[(l) * (NGLLX) + K];
    offset1 = (K) * (NGLL2) + (J) * (NGLLX) + l;
    offset2 = (K) * (NGLL2) + (l) * (NGLLX) + I;
    offset3 = (l) * (NGLL2) + (J) * (NGLLX) + I;
    temp1l = temp1l + (scalar_field[offset1]) * (hp1);
    temp2l = temp2l + (scalar_field[offset2]) * (hp2);
    temp3l = temp3l + (scalar_field[offset3]) * (hp3);
  }
  offset = (ispec) * (NGLL3_PADDED) + ijk;
  xixl = d_xix[offset];
  xiyl = d_xiy[offset];
  xizl = d_xiz[offset];
  etaxl = d_etax[offset];
  etayl = d_etay[offset];
  etazl = d_etaz[offset];
  gammaxl = d_gammax[offset];
  gammayl = d_gammay[offset];
  gammazl = d_gammaz[offset];
  vector_field_element[0] = (temp1l) * (xixl) + (temp2l) * (etaxl) + (temp3l) * (gammaxl);
  vector_field_element[1] = (temp1l) * (xiyl) + (temp2l) * (etayl) + (temp3l) * (gammayl);
  vector_field_element[2] = (temp1l) * (xizl) + (temp2l) * (etazl) + (temp3l) * (gammazl);
}
__global__ void compute_acoustic_kernel(const int * ibool, const float * rhostore, const float * kappastore, const float * hprime_xx, const float * d_xix, const float * d_xiy, const float * d_xiz, const float * d_etax, const float * d_etay, const float * d_etaz, const float * d_gammax, const float * d_gammay, const float * d_gammaz, const float * potential_dot_dot_acoustic, const float * b_potential_acoustic, const float * b_potential_dot_dot_acoustic, float * rho_ac_kl, float * kappa_ac_kl, const float deltat, const int NSPEC){
  int ispec;
  int ijk;
  int ijk_ispec;
  int ijk_ispec_padded;
  int iglob;
  float accel_elm[(3)];
  float b_displ_elm[(3)];
  float rhol;
  float kappal;
  float div_displ;
  float b_div_displ;
  __shared__ float scalar_field_displ[(NGLL3)];
  __shared__ float scalar_field_accel[(NGLL3)];
  ispec = blockIdx.x + (blockIdx.y) * (gridDim.x);
  if (ispec < NSPEC) {
    ijk = threadIdx.x;
    ijk_ispec = ijk + (NGLL3) * (ispec);
    ijk_ispec_padded = ijk + (NGLL3_PADDED) * (ispec);
    iglob = ibool[ijk_ispec] - (1);
    scalar_field_displ[ijk] = b_potential_acoustic[iglob];
    scalar_field_accel[ijk] = potential_dot_dot_acoustic[iglob];
    __syncthreads();
    compute_gradient_kernel(ijk, ispec, scalar_field_displ, b_displ_elm, hprime_xx, d_xix, d_xiy, d_xiz, d_etax, d_etay, d_etaz, d_gammax, d_gammay, d_gammaz);
    compute_gradient_kernel(ijk, ispec, scalar_field_accel, accel_elm, hprime_xx, d_xix, d_xiy, d_xiz, d_etax, d_etay, d_etaz, d_gammax, d_gammay, d_gammaz);
    rhol = rhostore[ijk_ispec_padded];
    rho_ac_kl[ijk_ispec] = rho_ac_kl[ijk_ispec] + ((deltat) * (rhol)) * ((accel_elm[0]) * (b_displ_elm[0]) + (accel_elm[1]) * (b_displ_elm[1]) + (accel_elm[2]) * (b_displ_elm[2]));
    kappal = (rhol) / (kappastore[ijk_ispec_padded]);
    div_displ = (kappal) * (potential_dot_dot_acoustic[iglob]);
    b_div_displ = (kappal) * (b_potential_dot_dot_acoustic[iglob]);
    kappa_ac_kl[ijk_ispec] = kappa_ac_kl[ijk_ispec] + ((deltat) * (div_displ)) * (b_div_displ);
  }
}
