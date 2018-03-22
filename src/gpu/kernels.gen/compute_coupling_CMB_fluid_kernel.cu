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

__global__ void compute_coupling_CMB_fluid_kernel(const float * displ_crust_mantle, float * accel_crust_mantle, const float * accel_outer_core, const int * ibool_crust_mantle, const int * ibelm_bottom_crust_mantle, const float * normal_top_outer_core, const float * jacobian2D_top_outer_core, const float * wgllwgll_xy, const int * ibool_outer_core, const int * ibelm_top_outer_core, const float RHO_TOP_OC, const float minus_g_cmb, const int GRAVITY, const int NSPEC2D_BOTTOM_CM){
  int i;
  int j;
  int k;
  int iface;
  int k_corresp;
  int iglob_oc;
  int iglob_cm;
  float pressure;
  int ispec;
  int ispec_selected;
  float nx;
  float ny;
  float nz;
  float weight;
  i = threadIdx.x;
  j = threadIdx.y;
  iface = blockIdx.x + (gridDim.x) * (blockIdx.y);
  if (iface < NSPEC2D_BOTTOM_CM) {
    ispec = ibelm_bottom_crust_mantle[iface] - (1);
    ispec_selected = ibelm_top_outer_core[iface] - (1);
    k = 0;
    k_corresp = NGLLX - (1);
    iglob_oc = ibool_outer_core[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k_corresp, ispec_selected)] - (1);
    nx = normal_top_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 0, i, j, iface)];
    ny = normal_top_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 1, i, j, iface)];
    nz = normal_top_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 2, i, j, iface)];
    weight = (jacobian2D_top_outer_core[INDEX3(NGLLX, NGLLX, i, j, iface)]) * (wgllwgll_xy[INDEX2(NGLLX, i, j)]);
    iglob_cm = ibool_crust_mantle[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec)] - (1);
    if (GRAVITY) {
      pressure = (RHO_TOP_OC) * ((minus_g_cmb) * ((displ_crust_mantle[(iglob_cm) * (3)]) * (nx) + (displ_crust_mantle[(iglob_cm) * (3) + 1]) * (ny) + (displ_crust_mantle[(iglob_cm) * (3) + 2]) * (nz)) - (accel_outer_core[iglob_oc]));
    } else {
      pressure = ( -(RHO_TOP_OC)) * (accel_outer_core[iglob_oc]);
    }
    atomicAdd(accel_crust_mantle + (iglob_cm) * (3) + 0, ((weight) * (nx)) * (pressure));
    atomicAdd(accel_crust_mantle + (iglob_cm) * (3) + 1, ((weight) * (ny)) * (pressure));
    atomicAdd(accel_crust_mantle + (iglob_cm) * (3) + 2, ((weight) * (nz)) * (pressure));
  }
}
