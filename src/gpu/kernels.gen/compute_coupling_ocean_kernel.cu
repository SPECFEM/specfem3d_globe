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

__global__ void compute_coupling_ocean_kernel(float * accel_crust_mantle, const float * rmassx_crust_mantle, const float * rmassy_crust_mantle, const float * rmassz_crust_mantle, const float * rmass_ocean_load, const int npoin_ocean_load, const int * ibool_ocean_load, const float * normal_ocean_load){
  int ipoin;
  int iglob;
  float nx;
  float ny;
  float nz;
  float rmass;
  float force_normal_comp;
  float additional_term_x;
  float additional_term_y;
  float additional_term_z;
  ipoin = threadIdx.x + (blockIdx.x) * (blockDim.x) + ((gridDim.x) * (blockDim.x)) * (threadIdx.y + (blockIdx.y) * (blockDim.y));
  if (ipoin < npoin_ocean_load) {
    iglob = ibool_ocean_load[ipoin] - (1);
    nx = normal_ocean_load[INDEX2(NDIM, 0, ipoin)];
    ny = normal_ocean_load[INDEX2(NDIM, 1, ipoin)];
    nz = normal_ocean_load[INDEX2(NDIM, 2, ipoin)];
    force_normal_comp = ((accel_crust_mantle[0 + (3) * (iglob)]) * (nx)) / (rmassx_crust_mantle[iglob]) + ((accel_crust_mantle[1 + (3) * (iglob)]) * (ny)) / (rmassy_crust_mantle[iglob]) + ((accel_crust_mantle[2 + (3) * (iglob)]) * (nz)) / (rmassz_crust_mantle[iglob]);
    rmass = rmass_ocean_load[ipoin];
    additional_term_x = (rmass - (rmassx_crust_mantle[iglob])) * (force_normal_comp);
    additional_term_y = (rmass - (rmassy_crust_mantle[iglob])) * (force_normal_comp);
    additional_term_z = (rmass - (rmassz_crust_mantle[iglob])) * (force_normal_comp);
    accel_crust_mantle[0 + (3) * (iglob)] = accel_crust_mantle[0 + (3) * (iglob)] + (additional_term_x) * (nx);
    accel_crust_mantle[1 + (3) * (iglob)] = accel_crust_mantle[1 + (3) * (iglob)] + (additional_term_y) * (ny);
    accel_crust_mantle[2 + (3) * (iglob)] = accel_crust_mantle[2 + (3) * (iglob)] + (additional_term_z) * (nz);
  }
}
