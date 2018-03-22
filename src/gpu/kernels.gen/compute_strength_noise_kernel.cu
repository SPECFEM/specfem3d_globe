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

__global__ void compute_strength_noise_kernel(const float * displ, const int * ibelm_top, const int * ibool, const float * noise_surface_movie, const float * normal_x_noise, const float * normal_y_noise, const float * normal_z_noise, float * Sigma_kl, const float deltat, const int nspec_top){
  int iface;
  int ispec;
  int igll;
  int ipoin;
  int i;
  int j;
  int k;
  int iglob;
  float eta;
  iface = blockIdx.x + (blockIdx.y) * (gridDim.x);
  if (iface < nspec_top) {
    ispec = ibelm_top[iface] - (1);
    igll = threadIdx.x;
    ipoin = igll + (NGLL2) * (iface);
    k = NGLLX - (1);
    j = (igll) / (NGLLX);
    i = igll - ((j) * (NGLLX));
    iglob = ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec)] - (1);
    eta = (noise_surface_movie[INDEX3(NDIM, NGLL2, 0, igll, iface)]) * (normal_x_noise[ipoin]) + (noise_surface_movie[INDEX3(NDIM, NGLL2, 1, igll, iface)]) * (normal_y_noise[ipoin]) + (noise_surface_movie[INDEX3(NDIM, NGLL2, 2, igll, iface)]) * (normal_z_noise[ipoin]);
    Sigma_kl[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec)] = Sigma_kl[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec)] + ((deltat) * (eta)) * ((normal_x_noise[ipoin]) * (displ[0 + (3) * (iglob)]) + (normal_y_noise[ipoin]) * (displ[1 + (3) * (iglob)]) + (normal_z_noise[ipoin]) * (displ[2 + (3) * (iglob)]));
  }
}
