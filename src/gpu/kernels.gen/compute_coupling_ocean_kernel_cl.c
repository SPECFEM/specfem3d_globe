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

const char * compute_coupling_ocean_kernel_program = "\
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
#ifndef R_EARTH_KM\n\
#define R_EARTH_KM 6371.0f\n\
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
__kernel void compute_coupling_ocean_kernel(__global float * accel_crust_mantle, const __global float * rmassx_crust_mantle, const __global float * rmassy_crust_mantle, const __global float * rmassz_crust_mantle, const __global float * rmass_ocean_load, const int npoin_ocean_load, const __global int * ibool_ocean_load, const __global float * normal_ocean_load){\n\
  int ipoin;\n\
  int iglob;\n\
  float nx;\n\
  float ny;\n\
  float nz;\n\
  float rmass;\n\
  float force_normal_comp;\n\
  float additional_term_x;\n\
  float additional_term_y;\n\
  float additional_term_z;\n\
  ipoin = get_global_id(0) + (get_global_size(0)) * (get_global_id(1));\n\
  if (ipoin < npoin_ocean_load) {\n\
    iglob = ibool_ocean_load[ipoin] - (1);\n\
    nx = normal_ocean_load[INDEX2(NDIM, 0, ipoin)];\n\
    ny = normal_ocean_load[INDEX2(NDIM, 1, ipoin)];\n\
    nz = normal_ocean_load[INDEX2(NDIM, 2, ipoin)];\n\
    force_normal_comp = ((accel_crust_mantle[0 + (3) * (iglob)]) * (nx)) / (rmassx_crust_mantle[iglob]) + ((accel_crust_mantle[1 + (3) * (iglob)]) * (ny)) / (rmassy_crust_mantle[iglob]) + ((accel_crust_mantle[2 + (3) * (iglob)]) * (nz)) / (rmassz_crust_mantle[iglob]);\n\
    rmass = rmass_ocean_load[ipoin];\n\
    additional_term_x = (rmass - (rmassx_crust_mantle[iglob])) * (force_normal_comp);\n\
    additional_term_y = (rmass - (rmassy_crust_mantle[iglob])) * (force_normal_comp);\n\
    additional_term_z = (rmass - (rmassz_crust_mantle[iglob])) * (force_normal_comp);\n\
    accel_crust_mantle[0 + (3) * (iglob)] = accel_crust_mantle[0 + (3) * (iglob)] + (additional_term_x) * (nx);\n\
    accel_crust_mantle[1 + (3) * (iglob)] = accel_crust_mantle[1 + (3) * (iglob)] + (additional_term_y) * (ny);\n\
    accel_crust_mantle[2 + (3) * (iglob)] = accel_crust_mantle[2 + (3) * (iglob)] + (additional_term_z) * (nz);\n\
  }\n\
}\n\
";
