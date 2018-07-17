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

const char * compute_ani_kernel_program = "\
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
#if __OPENCL_C_VERSION__ && __OPENCL_C_VERSION__ >= 120\n\
static\n\
#endif\n\
void compute_strain_product(float * prod, const float eps_trace_over_3, const float * epsdev, const float b_eps_trace_over_3, const float * b_epsdev){\n\
  float eps[(6)];\n\
  float b_eps[(6)];\n\
  eps[0] = epsdev[0] + eps_trace_over_3;\n\
  eps[1] = epsdev[1] + eps_trace_over_3;\n\
  eps[2] =  -(eps[0] + eps[1]) + (eps_trace_over_3) * (3.0f);\n\
  eps[3] = epsdev[4];\n\
  eps[4] = epsdev[3];\n\
  eps[5] = epsdev[2];\n\
  b_eps[0] = b_epsdev[0] + b_eps_trace_over_3;\n\
  b_eps[1] = b_epsdev[1] + b_eps_trace_over_3;\n\
  b_eps[2] =  -(b_eps[0] + b_eps[1]) + (b_eps_trace_over_3) * (3.0f);\n\
  b_eps[3] = b_epsdev[4];\n\
  b_eps[4] = b_epsdev[3];\n\
  b_eps[5] = b_epsdev[2];\n\
  prod[0] = (eps[0]) * (b_eps[0]);\n\
  prod[1] = (eps[0]) * (b_eps[1]);\n\
  prod[1] = prod[1] + (eps[1]) * (b_eps[0]);\n\
  prod[2] = (eps[0]) * (b_eps[2]);\n\
  prod[2] = prod[2] + (eps[2]) * (b_eps[0]);\n\
  prod[3] = (eps[0]) * (b_eps[3]);\n\
  prod[3] = prod[3] + (eps[3]) * (b_eps[0]);\n\
  prod[3] = (prod[3]) * (2.0f);\n\
  prod[4] = (eps[0]) * (b_eps[4]);\n\
  prod[4] = prod[4] + (eps[4]) * (b_eps[0]);\n\
  prod[4] = (prod[4]) * (2.0f);\n\
  prod[5] = (eps[0]) * (b_eps[5]);\n\
  prod[5] = prod[5] + (eps[5]) * (b_eps[0]);\n\
  prod[5] = (prod[5]) * (2.0f);\n\
  prod[6] = (eps[1]) * (b_eps[1]);\n\
  prod[7] = (eps[1]) * (b_eps[2]);\n\
  prod[7] = prod[7] + (eps[2]) * (b_eps[1]);\n\
  prod[8] = (eps[1]) * (b_eps[3]);\n\
  prod[8] = prod[8] + (eps[3]) * (b_eps[1]);\n\
  prod[8] = (prod[8]) * (2.0f);\n\
  prod[9] = (eps[1]) * (b_eps[4]);\n\
  prod[9] = prod[9] + (eps[4]) * (b_eps[1]);\n\
  prod[9] = (prod[9]) * (2.0f);\n\
  prod[10] = (eps[1]) * (b_eps[5]);\n\
  prod[10] = prod[10] + (eps[5]) * (b_eps[1]);\n\
  prod[10] = (prod[10]) * (2.0f);\n\
  prod[11] = (eps[2]) * (b_eps[2]);\n\
  prod[12] = (eps[2]) * (b_eps[3]);\n\
  prod[12] = prod[12] + (eps[3]) * (b_eps[2]);\n\
  prod[12] = (prod[12]) * (2.0f);\n\
  prod[13] = (eps[2]) * (b_eps[4]);\n\
  prod[13] = prod[13] + (eps[4]) * (b_eps[2]);\n\
  prod[13] = (prod[13]) * (2.0f);\n\
  prod[14] = (eps[2]) * (b_eps[5]);\n\
  prod[14] = prod[14] + (eps[5]) * (b_eps[2]);\n\
  prod[14] = (prod[14]) * (2.0f);\n\
  prod[15] = (eps[3]) * (b_eps[3]);\n\
  prod[15] = (prod[15]) * (4.0f);\n\
  prod[16] = (eps[3]) * (b_eps[4]);\n\
  prod[16] = prod[16] + (eps[4]) * (b_eps[3]);\n\
  prod[16] = (prod[16]) * (4.0f);\n\
  prod[17] = (eps[3]) * (b_eps[5]);\n\
  prod[17] = prod[17] + (eps[5]) * (b_eps[3]);\n\
  prod[17] = (prod[17]) * (4.0f);\n\
  prod[18] = (eps[4]) * (b_eps[4]);\n\
  prod[18] = (prod[18]) * (4.0f);\n\
  prod[19] = (eps[4]) * (b_eps[5]);\n\
  prod[19] = prod[19] + (eps[5]) * (b_eps[4]);\n\
  prod[19] = (prod[19]) * (4.0f);\n\
  prod[20] = (eps[5]) * (b_eps[5]);\n\
  prod[20] = (prod[20]) * (4.0f);\n\
}\n\
__kernel void compute_ani_kernel(const __global float * epsilondev_xx, const __global float * epsilondev_yy, const __global float * epsilondev_xy, const __global float * epsilondev_xz, const __global float * epsilondev_yz, const __global float * epsilon_trace_over_3, const __global float * b_epsilondev_xx, const __global float * b_epsilondev_yy, const __global float * b_epsilondev_xy, const __global float * b_epsilondev_xz, const __global float * b_epsilondev_yz, const __global float * b_epsilon_trace_over_3, __global float * cijkl_kl, const int NSPEC, const float deltat){\n\
  int i;\n\
  int ispec;\n\
  int ijk_ispec;\n\
  int offset;\n\
  float eps_trace_over_3;\n\
  float b_eps_trace_over_3;\n\
  float prod[(21)];\n\
  float epsdev[(5)];\n\
  float b_epsdev[(5)];\n\
  ispec = get_group_id(0) + (get_group_id(1)) * (get_num_groups(0));\n\
  if (ispec < NSPEC) {\n\
    ijk_ispec = get_local_id(0) + (NGLL3) * (ispec);\n\
    epsdev[0] = epsilondev_xx[ijk_ispec];\n\
    epsdev[1] = epsilondev_yy[ijk_ispec];\n\
    epsdev[2] = epsilondev_xy[ijk_ispec];\n\
    epsdev[3] = epsilondev_xz[ijk_ispec];\n\
    epsdev[4] = epsilondev_yz[ijk_ispec];\n\
    b_epsdev[0] = b_epsilondev_xx[ijk_ispec];\n\
    b_epsdev[1] = b_epsilondev_yy[ijk_ispec];\n\
    b_epsdev[2] = b_epsilondev_xy[ijk_ispec];\n\
    b_epsdev[3] = b_epsilondev_xz[ijk_ispec];\n\
    b_epsdev[4] = b_epsilondev_yz[ijk_ispec];\n\
    eps_trace_over_3 = epsilon_trace_over_3[ijk_ispec];\n\
    b_eps_trace_over_3 = b_epsilon_trace_over_3[ijk_ispec];\n\
    compute_strain_product(prod, eps_trace_over_3, epsdev, b_eps_trace_over_3, b_epsdev);\n\
    offset = ((ispec) * (NGLL3)) * (21) + get_local_id(0);\n\
    // attention: following array is sorted differently on GPU and CPU, -> use 'resort_array' before copying back to cpu\n\
    for (i = 0; i <= 20; i += 1) {\n\
      cijkl_kl[(i) * (NGLL3) + offset] = cijkl_kl[(i) * (NGLL3) + offset] + (deltat) * (prod[i]);\n\
    }\n\
  }\n\
}\n\
";
