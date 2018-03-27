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
#define NGLLX false
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

__global__ void compute_seismograms_kernel(const int nrec_local, const float * displ, const int * d_ibool, const float * xir, const float * etar, const float * gammar, float * seismograms, const float * nu, const int * ispec_selected_rec, const int * number_receiver_global, const float scale_displ){
  int ispec;
  int iglob;
  int irec_local;
  int irec;
  int tx;
  float lagrange;
  int i;
  int j;
  int k;
  int l;
  int s;
  __shared__ float sh_dxd[(NGLL3_PADDED)];
  __shared__ float sh_dyd[(NGLL3_PADDED)];
  __shared__ float sh_dzd[(NGLL3_PADDED)];
  tx = threadIdx.x;
  irec_local = blockIdx.x + (gridDim.x) * (blockIdx.y);
  k = (tx) / (NGLL2);
  j = (tx - ((k) * (NGLL2))) / (NGLLX);
  i = tx - ((k) * (NGLL2)) - ((j) * (NGLLX));
  if (irec_local < nrec_local) {
    irec = number_receiver_global[irec_local] - (1);
    ispec = ispec_selected_rec[irec] - (1);
    sh_dxd[tx] = 0;
    sh_dyd[tx] = 0;
    sh_dzd[tx] = 0;
    if (tx < NGLL3) {
      lagrange = ((xir[irec_local + (nrec_local) * (i)]) * (etar[irec_local + (nrec_local) * (j)])) * (gammar[irec_local + (nrec_local) * (k)]);
      iglob = d_ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec)] - (1);
      sh_dxd[tx] = (lagrange) * (displ[(iglob) * (3) + 0]);
      sh_dyd[tx] = (lagrange) * (displ[(iglob) * (3) + 1]);
      sh_dzd[tx] = (lagrange) * (displ[(iglob) * (3) + 2]);
    }
    __syncthreads();
    l = 1;
    s = (l) * (2);
    if (((tx < 0) ^ (s < 0) ? (tx % s) + s : tx % s) == 0) {
      sh_dxd[tx] = sh_dxd[tx] + sh_dxd[tx + l];
      sh_dyd[tx] = sh_dyd[tx] + sh_dyd[tx + l];
      sh_dzd[tx] = sh_dzd[tx] + sh_dzd[tx + l];
    }
    __syncthreads();
    l = (l) * (2);
    s = (l) * (2);
    if (((tx < 0) ^ (s < 0) ? (tx % s) + s : tx % s) == 0) {
      sh_dxd[tx] = sh_dxd[tx] + sh_dxd[tx + l];
      sh_dyd[tx] = sh_dyd[tx] + sh_dyd[tx + l];
      sh_dzd[tx] = sh_dzd[tx] + sh_dzd[tx + l];
    }
    __syncthreads();
    l = (l) * (2);
    s = (l) * (2);
    if (((tx < 0) ^ (s < 0) ? (tx % s) + s : tx % s) == 0) {
      sh_dxd[tx] = sh_dxd[tx] + sh_dxd[tx + l];
      sh_dyd[tx] = sh_dyd[tx] + sh_dyd[tx + l];
      sh_dzd[tx] = sh_dzd[tx] + sh_dzd[tx + l];
    }
    __syncthreads();
    l = (l) * (2);
    s = (l) * (2);
    if (((tx < 0) ^ (s < 0) ? (tx % s) + s : tx % s) == 0) {
      sh_dxd[tx] = sh_dxd[tx] + sh_dxd[tx + l];
      sh_dyd[tx] = sh_dyd[tx] + sh_dyd[tx + l];
      sh_dzd[tx] = sh_dzd[tx] + sh_dzd[tx + l];
    }
    __syncthreads();
    l = (l) * (2);
    s = (l) * (2);
    if (((tx < 0) ^ (s < 0) ? (tx % s) + s : tx % s) == 0) {
      sh_dxd[tx] = sh_dxd[tx] + sh_dxd[tx + l];
      sh_dyd[tx] = sh_dyd[tx] + sh_dyd[tx + l];
      sh_dzd[tx] = sh_dzd[tx] + sh_dzd[tx + l];
    }
    __syncthreads();
    l = (l) * (2);
    s = (l) * (2);
    if (((tx < 0) ^ (s < 0) ? (tx % s) + s : tx % s) == 0) {
      sh_dxd[tx] = sh_dxd[tx] + sh_dxd[tx + l];
      sh_dyd[tx] = sh_dyd[tx] + sh_dyd[tx + l];
      sh_dzd[tx] = sh_dzd[tx] + sh_dzd[tx + l];
    }
    __syncthreads();
    l = (l) * (2);
    s = (l) * (2);
    if (((tx < 0) ^ (s < 0) ? (tx % s) + s : tx % s) == 0) {
      sh_dxd[tx] = sh_dxd[tx] + sh_dxd[tx + l];
      sh_dyd[tx] = sh_dyd[tx] + sh_dyd[tx + l];
      sh_dzd[tx] = sh_dzd[tx] + sh_dzd[tx + l];
    }
    __syncthreads();
    l = (l) * (2);
    if (tx == 0) {
      seismograms[(irec_local) * (3) + 0] = (scale_displ) * ((nu[((irec_local) * (3)) * (3) + 0]) * (sh_dxd[0]) + (nu[((irec_local) * (3) + 1) * (3) + 0]) * (sh_dyd[0]) + (nu[((irec_local) * (3) + 2) * (3) + 0]) * (sh_dzd[0]));
    }
    if (tx == 1) {
      seismograms[(irec_local) * (3) + 1] = (scale_displ) * ((nu[((irec_local) * (3)) * (3) + 1]) * (sh_dxd[0]) + (nu[((irec_local) * (3) + 1) * (3) + 1]) * (sh_dyd[0]) + (nu[((irec_local) * (3) + 2) * (3) + 1]) * (sh_dzd[0]));
    }
    if (tx == 2) {
      seismograms[(irec_local) * (3) + 2] = (scale_displ) * ((nu[((irec_local) * (3)) * (3) + 2]) * (sh_dxd[0]) + (nu[((irec_local) * (3) + 1) * (3) + 2]) * (sh_dyd[0]) + (nu[((irec_local) * (3) + 2) * (3) + 2]) * (sh_dzd[0]));
    }
  }
}
