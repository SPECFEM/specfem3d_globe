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

const char * compute_acoustic_kernel_program = "\
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
void compute_gradient_kernel(const int ijk, const int ispec, const __local float * scalar_field, float * vector_field_element, const __global float * hprime_xx, const __global float * d_xix, const __global float * d_xiy, const __global float * d_xiz, const __global float * d_etax, const __global float * d_etay, const __global float * d_etaz, const __global float * d_gammax, const __global float * d_gammay, const __global float * d_gammaz){\n\
  float temp1l;\n\
  float temp2l;\n\
  float temp3l;\n\
  float hp1;\n\
  float hp2;\n\
  float hp3;\n\
  float xixl;\n\
  float xiyl;\n\
  float xizl;\n\
  float etaxl;\n\
  float etayl;\n\
  float etazl;\n\
  float gammaxl;\n\
  float gammayl;\n\
  float gammazl;\n\
  int l;\n\
  int offset;\n\
  int offset1;\n\
  int offset2;\n\
  int offset3;\n\
  int I;\n\
  int J;\n\
  int K;\n\
  K = (ijk) / (NGLL2);\n\
  J = (ijk - ((K) * (NGLL2))) / (NGLLX);\n\
  I = ijk - ((K) * (NGLL2)) - ((J) * (NGLLX));\n\
  temp1l = 0.0f;\n\
  temp2l = 0.0f;\n\
  temp3l = 0.0f;\n\
  for (l = 0; l <= NGLLX - (1); l += 1) {\n\
    hp1 = hprime_xx[(l) * (NGLLX) + I];\n\
    hp2 = hprime_xx[(l) * (NGLLX) + J];\n\
    hp3 = hprime_xx[(l) * (NGLLX) + K];\n\
    offset1 = (K) * (NGLL2) + (J) * (NGLLX) + l;\n\
    offset2 = (K) * (NGLL2) + (l) * (NGLLX) + I;\n\
    offset3 = (l) * (NGLL2) + (J) * (NGLLX) + I;\n\
    temp1l = temp1l + (scalar_field[offset1]) * (hp1);\n\
    temp2l = temp2l + (scalar_field[offset2]) * (hp2);\n\
    temp3l = temp3l + (scalar_field[offset3]) * (hp3);\n\
  }\n\
  offset = (ispec) * (NGLL3_PADDED) + ijk;\n\
  xixl = d_xix[offset];\n\
  xiyl = d_xiy[offset];\n\
  xizl = d_xiz[offset];\n\
  etaxl = d_etax[offset];\n\
  etayl = d_etay[offset];\n\
  etazl = d_etaz[offset];\n\
  gammaxl = d_gammax[offset];\n\
  gammayl = d_gammay[offset];\n\
  gammazl = d_gammaz[offset];\n\
  vector_field_element[0] = (temp1l) * (xixl) + (temp2l) * (etaxl) + (temp3l) * (gammaxl);\n\
  vector_field_element[1] = (temp1l) * (xiyl) + (temp2l) * (etayl) + (temp3l) * (gammayl);\n\
  vector_field_element[2] = (temp1l) * (xizl) + (temp2l) * (etazl) + (temp3l) * (gammazl);\n\
}\n\
__kernel void compute_acoustic_kernel(const __global int * ibool, const __global float * rhostore, const __global float * kappastore, const __global float * hprime_xx, const __global float * d_xix, const __global float * d_xiy, const __global float * d_xiz, const __global float * d_etax, const __global float * d_etay, const __global float * d_etaz, const __global float * d_gammax, const __global float * d_gammay, const __global float * d_gammaz, const __global float * potential_dot_dot_acoustic, const __global float * b_potential_acoustic, const __global float * b_potential_dot_dot_acoustic, __global float * rho_ac_kl, __global float * kappa_ac_kl, const float deltat, const int NSPEC){\n\
  int ispec;\n\
  int ijk;\n\
  int ijk_ispec;\n\
  int ijk_ispec_padded;\n\
  int iglob;\n\
  float accel_elm[(3)];\n\
  float b_displ_elm[(3)];\n\
  float rhol;\n\
  float kappal;\n\
  float div_displ;\n\
  float b_div_displ;\n\
  __local float scalar_field_displ[(NGLL3)];\n\
  __local float scalar_field_accel[(NGLL3)];\n\
  ispec = get_group_id(0) + (get_group_id(1)) * (get_num_groups(0));\n\
  if (ispec < NSPEC) {\n\
    ijk = get_local_id(0);\n\
    ijk_ispec = ijk + (NGLL3) * (ispec);\n\
    ijk_ispec_padded = ijk + (NGLL3_PADDED) * (ispec);\n\
    iglob = ibool[ijk_ispec] - (1);\n\
    scalar_field_displ[ijk] = b_potential_acoustic[iglob];\n\
    scalar_field_accel[ijk] = potential_dot_dot_acoustic[iglob];\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
    compute_gradient_kernel(ijk, ispec, scalar_field_displ, b_displ_elm, hprime_xx, d_xix, d_xiy, d_xiz, d_etax, d_etay, d_etaz, d_gammax, d_gammay, d_gammaz);\n\
    compute_gradient_kernel(ijk, ispec, scalar_field_accel, accel_elm, hprime_xx, d_xix, d_xiy, d_xiz, d_etax, d_etay, d_etaz, d_gammax, d_gammay, d_gammaz);\n\
    rhol = rhostore[ijk_ispec_padded];\n\
    rho_ac_kl[ijk_ispec] = rho_ac_kl[ijk_ispec] + ((deltat) * (rhol)) * ((accel_elm[0]) * (b_displ_elm[0]) + (accel_elm[1]) * (b_displ_elm[1]) + (accel_elm[2]) * (b_displ_elm[2]));\n\
    kappal = (rhol) / (kappastore[ijk_ispec_padded]);\n\
    div_displ = (kappal) * (potential_dot_dot_acoustic[iglob]);\n\
    b_div_displ = (kappal) * (b_potential_dot_dot_acoustic[iglob]);\n\
    kappa_ac_kl[ijk_ispec] = kappa_ac_kl[ijk_ispec] + ((deltat) * (div_displ)) * (b_div_displ);\n\
  }\n\
}\n\
";
