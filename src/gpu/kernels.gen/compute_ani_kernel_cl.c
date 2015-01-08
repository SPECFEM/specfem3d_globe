//note: please do not modify this file manually!
//      this file has been generated automatically by BOAST version 0.99996
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
! the Free Software Foundation; either version 2 of the License, or
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
static void compute_strain_product(float * prod, const float eps_trace_over_3, const float * epsdev, const float b_eps_trace_over_3, const float * b_epsdev){\n\
  float eps[(6)];\n\
  float b_eps[(6)];\n\
  int p;\n\
  int i;\n\
  int j;\n\
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
  p = 0;\n\
  for (i = 0; i <= 5; i += 1) {\n\
    for (j = i; j <= 5; j += 1) {\n\
      prod[p] = (eps[i]) * (b_eps[j]);\n\
      if (j > i) {\n\
        prod[p] = prod[p] + (eps[j]) * (b_eps[i]);\n\
        if (j > 2 && i < 3) {\n\
          prod[p] = (prod[p]) * (2.0f);\n\
        }\n\
      }\n\
      if (i > 2) {\n\
        prod[p] = (prod[p]) * (4.0f);\n\
      }\n\
      p = p + 1;\n\
    }\n\
  }\n\
}\n\
__kernel void compute_ani_kernel(const __global float * epsilondev_xx, const __global float * epsilondev_yy, const __global float * epsilondev_xy, const __global float * epsilondev_xz, const __global float * epsilondev_yz, const __global float * epsilon_trace_over_3, const __global float * b_epsilondev_xx, const __global float * b_epsilondev_yy, const __global float * b_epsilondev_xy, const __global float * b_epsilondev_xz, const __global float * b_epsilondev_yz, const __global float * b_epsilon_trace_over_3, __global float * cijkl_kl, const int NSPEC, const float deltat){\n\
  int i;\n\
  int ispec;\n\
  int ijk_ispec;\n\
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
    for (i = 0; i <= 20; i += 1) {\n\
      cijkl_kl[i + (21) * (ijk_ispec)] = cijkl_kl[i + (21) * (ijk_ispec)] + (deltat) * (prod[i]);\n\
    }\n\
  }\n\
}\n\
";
