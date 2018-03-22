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

__global__ void write_seismograms_transfer_from_device_kernel(const int * number_receiver_global, const int * ispec_selected_rec, const int * ibool, float * station_seismo_field, const float * d_field, const int nrec_local){
  int tx;
  int irec;
  int ispec;
  int iglob;
  int blockID;
  blockID = blockIdx.x + (blockIdx.y) * (gridDim.x);
  tx = threadIdx.x;
  if (blockID < nrec_local) {
    irec = number_receiver_global[blockID] - (1);
    ispec = ispec_selected_rec[irec] - (1);
    iglob = ibool[tx + (NGLL3) * (ispec)] - (1);
    station_seismo_field[((NGLL3) * (3)) * (blockID) + (tx) * (3) + 0] = d_field[(iglob) * (3) + 0];
    station_seismo_field[((NGLL3) * (3)) * (blockID) + (tx) * (3) + 1] = d_field[(iglob) * (3) + 1];
    station_seismo_field[((NGLL3) * (3)) * (blockID) + (tx) * (3) + 2] = d_field[(iglob) * (3) + 2];
  }
}
