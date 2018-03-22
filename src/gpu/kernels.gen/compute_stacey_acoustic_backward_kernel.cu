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

__global__ void compute_stacey_acoustic_backward_kernel(float * b_potential_dot_dot_acoustic, const float * b_absorb_potential, const int interface_type, const int num_abs_boundary_faces, const int * abs_boundary_ispec, const int * nkmin_xi, const int * nkmin_eta, const int * njmin, const int * njmax, const int * nimin, const int * nimax, const int * ibool){
  int igll;
  int iface;
  int i;
  int j;
  int k;
  int iglob;
  int ispec;
  igll = threadIdx.x;
  iface = blockIdx.x + (blockIdx.y) * (gridDim.x);
  if (iface < num_abs_boundary_faces) {
    ispec = abs_boundary_ispec[iface] - (1);
    switch (interface_type) {
      case 4 :
        if (nkmin_xi[INDEX2(2, 0, iface)] == 0 || njmin[INDEX2(2, 0, iface)] == 0) {
           return ;
        }
        i = 0;
        k = (igll) / (NGLLX);
        j = igll - ((k) * (NGLLX));
        if (k < nkmin_xi[INDEX2(2, 0, iface)] - (1) || k > NGLLX - (1)) {
           return ;
        }
        if (j < njmin[INDEX2(2, 0, iface)] - (1) || j > njmax[INDEX2(2, 0, iface)] - (1)) {
           return ;
        }
        break;
      case 5 :
        if (nkmin_xi[INDEX2(2, 1, iface)] == 0 || njmin[INDEX2(2, 1, iface)] == 0) {
           return ;
        }
        i = NGLLX - (1);
        k = (igll) / (NGLLX);
        j = igll - ((k) * (NGLLX));
        if (k < nkmin_xi[INDEX2(2, 1, iface)] - (1) || k > NGLLX - (1)) {
           return ;
        }
        if (j < njmin[INDEX2(2, 1, iface)] - (1) || j > njmax[INDEX2(2, 1, iface)] - (1)) {
           return ;
        }
        break;
      case 6 :
        if (nkmin_eta[INDEX2(2, 0, iface)] == 0 || nimin[INDEX2(2, 0, iface)] == 0) {
           return ;
        }
        j = 0;
        k = (igll) / (NGLLX);
        i = igll - ((k) * (NGLLX));
        if (k < nkmin_eta[INDEX2(2, 0, iface)] - (1) || k > NGLLX - (1)) {
           return ;
        }
        if (i < nimin[INDEX2(2, 0, iface)] - (1) || i > nimax[INDEX2(2, 0, iface)] - (1)) {
           return ;
        }
        break;
      case 7 :
        if (nkmin_eta[INDEX2(2, 1, iface)] == 0 || nimin[INDEX2(2, 1, iface)] == 0) {
           return ;
        }
        j = NGLLX - (1);
        k = (igll) / (NGLLX);
        i = igll - ((k) * (NGLLX));
        if (k < nkmin_eta[INDEX2(2, 1, iface)] - (1) || k > NGLLX - (1)) {
           return ;
        }
        if (i < nimin[INDEX2(2, 1, iface)] - (1) || i > nimax[INDEX2(2, 1, iface)] - (1)) {
           return ;
        }
        break;
      case 8 :
        k = 0;
        j = (igll) / (NGLLX);
        i = igll - ((j) * (NGLLX));
        if (j < 0 || j > NGLLX - (1)) {
           return ;
        }
        if (i < 0 || i > NGLLX - (1)) {
           return ;
        }
        break;
    }
    iglob = ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec)] - (1);
    atomicAdd(b_potential_dot_dot_acoustic + iglob,  -(b_absorb_potential[INDEX2(NGLL2, igll, iface)]));
  }
}
