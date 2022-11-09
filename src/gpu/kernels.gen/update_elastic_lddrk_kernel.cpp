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

__global__ void update_elastic_lddrk_kernel(float * displ, float * veloc, const float * accel, float * displ_lddrk, float * veloc_lddrk, const float alpha_lddrk, const float beta_lddrk, const float deltat, const int size){
  int id;

  id = threadIdx.x + (blockIdx.x) * (blockDim.x) + (blockIdx.y) * ((gridDim.x) * (blockDim.x));

  if (id < size) {
    veloc_lddrk[(id) * (3)] = (alpha_lddrk) * (veloc_lddrk[(id) * (3)]) + (deltat) * (accel[(id) * (3)]);
    veloc_lddrk[(id) * (3) + 1] = (alpha_lddrk) * (veloc_lddrk[(id) * (3) + 1]) + (deltat) * (accel[(id) * (3) + 1]);
    veloc_lddrk[(id) * (3) + 2] = (alpha_lddrk) * (veloc_lddrk[(id) * (3) + 2]) + (deltat) * (accel[(id) * (3) + 2]);
    displ_lddrk[(id) * (3)] = (alpha_lddrk) * (displ_lddrk[(id) * (3)]) + (deltat) * (veloc[(id) * (3)]);
    displ_lddrk[(id) * (3) + 1] = (alpha_lddrk) * (displ_lddrk[(id) * (3) + 1]) + (deltat) * (veloc[(id) * (3) + 1]);
    displ_lddrk[(id) * (3) + 2] = (alpha_lddrk) * (displ_lddrk[(id) * (3) + 2]) + (deltat) * (veloc[(id) * (3) + 2]);
    veloc[(id) * (3)] = veloc[(id) * (3)] + (beta_lddrk) * (veloc_lddrk[(id) * (3)]);
    veloc[(id) * (3) + 1] = veloc[(id) * (3) + 1] + (beta_lddrk) * (veloc_lddrk[(id) * (3) + 1]);
    veloc[(id) * (3) + 2] = veloc[(id) * (3) + 2] + (beta_lddrk) * (veloc_lddrk[(id) * (3) + 2]);
    displ[(id) * (3)] = displ[(id) * (3)] + (beta_lddrk) * (displ_lddrk[(id) * (3)]);
    displ[(id) * (3) + 1] = displ[(id) * (3) + 1] + (beta_lddrk) * (displ_lddrk[(id) * (3) + 1]);
    displ[(id) * (3) + 2] = displ[(id) * (3) + 2] + (beta_lddrk) * (displ_lddrk[(id) * (3) + 2]);
  }
}
