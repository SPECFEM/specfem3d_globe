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

const char * compute_seismograms_kernel_program = "\
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
__kernel void compute_seismograms_kernel(const int nrec_local, const __global float * displ, const __global int * d_ibool, const __global float * xir, const __global float * etar, const __global float * gammar, __global float * seismograms, const __global float * nu, const __global int * ispec_selected_rec, const __global int * number_receiver_global, const float scale_displ){\n\
  int ispec;\n\
  int iglob;\n\
  int irec_local;\n\
  int irec;\n\
  int tx;\n\
  float lagrange;\n\
  int i;\n\
  int j;\n\
  int k;\n\
  int l;\n\
  int s;\n\
  __local float sh_dxd[(NGLL3_PADDED)];\n\
  __local float sh_dyd[(NGLL3_PADDED)];\n\
  __local float sh_dzd[(NGLL3_PADDED)];\n\
  tx = get_local_id(0);\n\
  irec_local = get_group_id(0) + (get_num_groups(0)) * (get_group_id(1));\n\
  k = (tx) / (NGLL2);\n\
  j = (tx - ((k) * (NGLL2))) / (NGLLX);\n\
  i = tx - ((k) * (NGLL2)) - ((j) * (NGLLX));\n\
  if (irec_local < nrec_local) {\n\
    irec = number_receiver_global[irec_local] - (1);\n\
    ispec = ispec_selected_rec[irec] - (1);\n\
    sh_dxd[tx] = 0;\n\
    sh_dyd[tx] = 0;\n\
    sh_dzd[tx] = 0;\n\
    if (tx < NGLL3) {\n\
      lagrange = ((xir[irec_local + (nrec_local) * (i)]) * (etar[irec_local + (nrec_local) * (j)])) * (gammar[irec_local + (nrec_local) * (k)]);\n\
      iglob = d_ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec)] - (1);\n\
      sh_dxd[tx] = (lagrange) * (displ[(iglob) * (3) + 0]);\n\
      sh_dyd[tx] = (lagrange) * (displ[(iglob) * (3) + 1]);\n\
      sh_dzd[tx] = (lagrange) * (displ[(iglob) * (3) + 2]);\n\
    }\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    l = 1;\n\
    s = (l) * (2);\n\
    if (((tx < 0) ^ (s < 0) ? (tx % s) + s : tx % s) == 0) {\n\
      sh_dxd[tx] = sh_dxd[tx] + sh_dxd[tx + l];\n\
      sh_dyd[tx] = sh_dyd[tx] + sh_dyd[tx + l];\n\
      sh_dzd[tx] = sh_dzd[tx] + sh_dzd[tx + l];\n\
    }\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    l = (l) * (2);\n\
    s = (l) * (2);\n\
    if (((tx < 0) ^ (s < 0) ? (tx % s) + s : tx % s) == 0) {\n\
      sh_dxd[tx] = sh_dxd[tx] + sh_dxd[tx + l];\n\
      sh_dyd[tx] = sh_dyd[tx] + sh_dyd[tx + l];\n\
      sh_dzd[tx] = sh_dzd[tx] + sh_dzd[tx + l];\n\
    }\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    l = (l) * (2);\n\
    s = (l) * (2);\n\
    if (((tx < 0) ^ (s < 0) ? (tx % s) + s : tx % s) == 0) {\n\
      sh_dxd[tx] = sh_dxd[tx] + sh_dxd[tx + l];\n\
      sh_dyd[tx] = sh_dyd[tx] + sh_dyd[tx + l];\n\
      sh_dzd[tx] = sh_dzd[tx] + sh_dzd[tx + l];\n\
    }\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    l = (l) * (2);\n\
    s = (l) * (2);\n\
    if (((tx < 0) ^ (s < 0) ? (tx % s) + s : tx % s) == 0) {\n\
      sh_dxd[tx] = sh_dxd[tx] + sh_dxd[tx + l];\n\
      sh_dyd[tx] = sh_dyd[tx] + sh_dyd[tx + l];\n\
      sh_dzd[tx] = sh_dzd[tx] + sh_dzd[tx + l];\n\
    }\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    l = (l) * (2);\n\
    s = (l) * (2);\n\
    if (((tx < 0) ^ (s < 0) ? (tx % s) + s : tx % s) == 0) {\n\
      sh_dxd[tx] = sh_dxd[tx] + sh_dxd[tx + l];\n\
      sh_dyd[tx] = sh_dyd[tx] + sh_dyd[tx + l];\n\
      sh_dzd[tx] = sh_dzd[tx] + sh_dzd[tx + l];\n\
    }\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    l = (l) * (2);\n\
    s = (l) * (2);\n\
    if (((tx < 0) ^ (s < 0) ? (tx % s) + s : tx % s) == 0) {\n\
      sh_dxd[tx] = sh_dxd[tx] + sh_dxd[tx + l];\n\
      sh_dyd[tx] = sh_dyd[tx] + sh_dyd[tx + l];\n\
      sh_dzd[tx] = sh_dzd[tx] + sh_dzd[tx + l];\n\
    }\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    l = (l) * (2);\n\
    s = (l) * (2);\n\
    if (((tx < 0) ^ (s < 0) ? (tx % s) + s : tx % s) == 0) {\n\
      sh_dxd[tx] = sh_dxd[tx] + sh_dxd[tx + l];\n\
      sh_dyd[tx] = sh_dyd[tx] + sh_dyd[tx + l];\n\
      sh_dzd[tx] = sh_dzd[tx] + sh_dzd[tx + l];\n\
    }\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    l = (l) * (2);\n\
    if (tx == 0) {\n\
      seismograms[(irec_local) * (3) + 0] = (scale_displ) * ((nu[((irec_local) * (3)) * (3) + 0]) * (sh_dxd[0]) + (nu[((irec_local) * (3) + 1) * (3) + 0]) * (sh_dyd[0]) + (nu[((irec_local) * (3) + 2) * (3) + 0]) * (sh_dzd[0]));\n\
    }\n\
    if (tx == 1) {\n\
      seismograms[(irec_local) * (3) + 1] = (scale_displ) * ((nu[((irec_local) * (3)) * (3) + 1]) * (sh_dxd[0]) + (nu[((irec_local) * (3) + 1) * (3) + 1]) * (sh_dyd[0]) + (nu[((irec_local) * (3) + 2) * (3) + 1]) * (sh_dzd[0]));\n\
    }\n\
    if (tx == 2) {\n\
      seismograms[(irec_local) * (3) + 2] = (scale_displ) * ((nu[((irec_local) * (3)) * (3) + 2]) * (sh_dxd[0]) + (nu[((irec_local) * (3) + 1) * (3) + 2]) * (sh_dyd[0]) + (nu[((irec_local) * (3) + 2) * (3) + 2]) * (sh_dzd[0]));\n\
    }\n\
  }\n\
}\n\
";
