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

__global__ void compute_stacey_elastic_kernel(const float * veloc, float * accel, const int interface_type, const int num_abs_boundary_faces, const int * abs_boundary_ispec, const int * nkmin_xi, const int * nkmin_eta, const int * njmin, const int * njmax, const int * nimin, const int * nimax, const float * abs_boundary_normal, const float * abs_boundary_jacobian2D, const float * wgllwgll, const int * ibool, const float * rho_vp, const float * rho_vs, const int SAVE_STACEY, float * b_absorb_field){
  int igll;
  int iface;
  int i;
  int j;
  int k;
  int iglob;
  int ispec;
  float vx;
  float vy;
  float vz;
  float vn;
  float nx;
  float ny;
  float nz;
  float rho_vp_temp;
  float rho_vs_temp;
  float tx;
  float ty;
  float tz;
  float jacobianw;
  float fac1;
  igll = threadIdx.x;
  iface = blockIdx.x + (blockIdx.y) * (gridDim.x);
  if (iface < num_abs_boundary_faces) {
    ispec = abs_boundary_ispec[iface] - (1);
    switch (interface_type) {
      case 0 :
        if (nkmin_xi[INDEX2(2, 0, iface)] == 0 || njmin[INDEX2(2, 0, iface)] == 0) {
           return ;
        }
        i = 0;
        k = (igll) / (NGLLX);
        j = igll - ((k) * (NGLLX));
        if (k < nkmin_xi[INDEX2(2, 0, iface)] - (1) || k > NGLLX - (1)) {
           return ;
        }
        if (j < njmin[INDEX2(2, 0, iface)] - (1) || j > NGLLX - (1)) {
           return ;
        }
        fac1 = wgllwgll[(k) * (NGLLX) + j];
        break;
      case 1 :
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
        fac1 = wgllwgll[(k) * (NGLLX) + j];
        break;
      case 2 :
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
        fac1 = wgllwgll[(k) * (NGLLX) + i];
        break;
      case 3 :
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
        fac1 = wgllwgll[(k) * (NGLLX) + i];
        break;
    }
    iglob = ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec)] - (1);
    vx = veloc[(iglob) * (3) + 0];
    vy = veloc[(iglob) * (3) + 1];
    vz = veloc[(iglob) * (3) + 2];
    nx = abs_boundary_normal[INDEX3(NDIM, NGLL2, 0, igll, iface)];
    ny = abs_boundary_normal[INDEX3(NDIM, NGLL2, 1, igll, iface)];
    nz = abs_boundary_normal[INDEX3(NDIM, NGLL2, 2, igll, iface)];
    vn = (vx) * (nx) + (vy) * (ny) + (vz) * (nz);
    rho_vp_temp = rho_vp[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec)];
    rho_vs_temp = rho_vs[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec)];
    tx = ((rho_vp_temp) * (vn)) * (nx) + (rho_vs_temp) * (vx - ((vn) * (nx)));
    ty = ((rho_vp_temp) * (vn)) * (ny) + (rho_vs_temp) * (vy - ((vn) * (ny)));
    tz = ((rho_vp_temp) * (vn)) * (nz) + (rho_vs_temp) * (vz - ((vn) * (nz)));
    jacobianw = (abs_boundary_jacobian2D[INDEX2(NGLL2, igll, iface)]) * (fac1);
    atomicAdd(accel + (iglob) * (3) + 0, ( -(tx)) * (jacobianw));
    atomicAdd(accel + (iglob) * (3) + 1, ( -(ty)) * (jacobianw));
    atomicAdd(accel + (iglob) * (3) + 2, ( -(tz)) * (jacobianw));
    if (SAVE_STACEY) {
      b_absorb_field[INDEX3(NDIM, NGLL2, 0, igll, iface)] = (tx) * (jacobianw);
      b_absorb_field[INDEX3(NDIM, NGLL2, 1, igll, iface)] = (ty) * (jacobianw);
      b_absorb_field[INDEX3(NDIM, NGLL2, 2, igll, iface)] = (tz) * (jacobianw);
    }
  }
}
