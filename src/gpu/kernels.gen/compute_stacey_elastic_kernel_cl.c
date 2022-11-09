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

const char * compute_stacey_elastic_kernel_program = "\
inline void atomicAdd(volatile __global float *source, const float val) {\n\
  union {\n\
    unsigned int iVal;\n\
    float fVal;\n\
  } res, orig;\n\
  do {\n\
    orig.fVal = *source;\n\
    res.fVal = orig.fVal + val;\n\
  } while (atomic_cmpxchg((volatile __global unsigned int *)source, orig.iVal, res.iVal) != orig.iVal);\n\
}\n\
#ifndef INDEX2\n\
#define INDEX2(isize,i,j) i + isize*j\n\
#endif\n\
#ifndef INDEX3\n\
#define INDEX3(isize,jsize,i,j,k) i + isize*(j + jsize*k)\n\
#endif\n\
#ifndef INDEX4\n\
#define INDEX4(isize,jsize,ksize,i,j,k,x) i + isize*(j + jsize*(k + ksize*x))\n\
#endif\n\
#ifndef INDEX5\n\
#define INDEX5(isize,jsize,ksize,xsize,i,j,k,x,y) i + isize*(j + jsize*(k + ksize*(x + xsize*y)))\n\
#endif\n\
\n\
#ifndef NDIM\n\
#define NDIM 3\n\
#endif\n\
#ifndef NGLLX\n\
#define NGLLX 5\n\
#endif\n\
#ifndef NGLL2\n\
#define NGLL2 25\n\
#endif\n\
#ifndef NGLL3\n\
#define NGLL3 125\n\
#endif\n\
#ifndef NGLL3_PADDED\n\
#define NGLL3_PADDED 128\n\
#endif\n\
#ifndef N_SLS\n\
#define N_SLS 3\n\
#endif\n\
#ifndef IREGION_CRUST_MANTLE\n\
#define IREGION_CRUST_MANTLE 1\n\
#endif\n\
#ifndef IREGION_INNER_CORE\n\
#define IREGION_INNER_CORE 3\n\
#endif\n\
#ifndef IFLAG_IN_FICTITIOUS_CUBE\n\
#define IFLAG_IN_FICTITIOUS_CUBE 11\n\
#endif\n\
#ifndef COLORING_MIN_NSPEC_INNER_CORE\n\
#define COLORING_MIN_NSPEC_INNER_CORE 1000\n\
#endif\n\
#ifndef COLORING_MIN_NSPEC_OUTER_CORE\n\
#define COLORING_MIN_NSPEC_OUTER_CORE 1000\n\
#endif\n\
#ifndef BLOCKSIZE_TRANSFER\n\
#define BLOCKSIZE_TRANSFER 256\n\
#endif\n\
\n\
__kernel void compute_stacey_elastic_kernel(const __global float * veloc, __global float * accel, const int num_abs_boundary_faces, const __global int * abs_boundary_ispec, const __global int * abs_boundary_npoin, const __global int * abs_boundary_ijk, const __global float * abs_boundary_normal, const __global float * abs_boundary_jacobian2Dw, const __global int * ibool, const __global float * rho_vp, const __global float * rho_vs, const int SAVE_STACEY, __global float * b_absorb_field){\n\
  int npoin;\n\
  int igll;\n\
  int iface;\n\
  int i;\n\
  int j;\n\
  int k;\n\
  int iglob;\n\
  int ispec;\n\
  float vx;\n\
  float vy;\n\
  float vz;\n\
  float vn;\n\
  float nx;\n\
  float ny;\n\
  float nz;\n\
  float rho_vp_temp;\n\
  float rho_vs_temp;\n\
  float tx;\n\
  float ty;\n\
  float tz;\n\
  float weight;\n\
\n\
  igll = get_local_id(0);\n\
  iface = get_group_id(0) + (get_group_id(1)) * (get_num_groups(0));\n\
\n\
  if (iface < num_abs_boundary_faces) {\n\
    npoin = abs_boundary_npoin[iface];\n\
\n\
    if (igll < npoin) {\n\
      ispec = abs_boundary_ispec[iface] - (1);\n\
\n\
      i = abs_boundary_ijk[INDEX3(3, NGLL2, 0, igll, iface)] - (1);\n\
      j = abs_boundary_ijk[INDEX3(3, NGLL2, 1, igll, iface)] - (1);\n\
      k = abs_boundary_ijk[INDEX3(3, NGLL2, 2, igll, iface)] - (1);\n\
\n\
      iglob = ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec)] - (1);\n\
\n\
      vx = veloc[(iglob) * (3) + 0];\n\
      vy = veloc[(iglob) * (3) + 1];\n\
      vz = veloc[(iglob) * (3) + 2];\n\
      nx = abs_boundary_normal[INDEX3(NDIM, NGLL2, 0, igll, iface)];\n\
      ny = abs_boundary_normal[INDEX3(NDIM, NGLL2, 1, igll, iface)];\n\
      nz = abs_boundary_normal[INDEX3(NDIM, NGLL2, 2, igll, iface)];\n\
\n\
      vn = (vx) * (nx) + (vy) * (ny) + (vz) * (nz);\n\
      rho_vp_temp = rho_vp[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec)];\n\
      rho_vs_temp = rho_vs[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec)];\n\
\n\
      tx = ((rho_vp_temp) * (vn)) * (nx) + (rho_vs_temp) * (vx - ((vn) * (nx)));\n\
      ty = ((rho_vp_temp) * (vn)) * (ny) + (rho_vs_temp) * (vy - ((vn) * (ny)));\n\
      tz = ((rho_vp_temp) * (vn)) * (nz) + (rho_vs_temp) * (vz - ((vn) * (nz)));\n\
\n\
      weight = abs_boundary_jacobian2Dw[INDEX2(NGLL2, igll, iface)];\n\
\n\
      atomicAdd(accel + (iglob) * (3) + 0, ( -(tx)) * (weight));\n\
      atomicAdd(accel + (iglob) * (3) + 1, ( -(ty)) * (weight));\n\
      atomicAdd(accel + (iglob) * (3) + 2, ( -(tz)) * (weight));\n\
\n\
      if (SAVE_STACEY) {\n\
        b_absorb_field[INDEX3(NDIM, NGLL2, 0, igll, iface)] = (tx) * (weight);\n\
        b_absorb_field[INDEX3(NDIM, NGLL2, 1, igll, iface)] = (ty) * (weight);\n\
        b_absorb_field[INDEX3(NDIM, NGLL2, 2, igll, iface)] = (tz) * (weight);\n\
      }\n\
    }\n\
  }\n\
}\n\
";
