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

const char * smooth_process_kernel_program = "\
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
__kernel void smooth_process_kernel(const __global float * xstore_me, const __global float * ystore_me, const __global float * zstore_me, const __global float * xstore_other, const __global float * ystore_other, const __global float * zstore_other, const __global float * data_other, const float sigma_h2_inv, const float sigma_v2_inv, const int iker, const int nspec_me, const int nspec_other, const float v_criterion, const float h_criterion, const __global float * integ_factor, __global float * data_smooth, __global float * normalisation, const int use_vector_distance){\n\
  int ispec;\n\
  int igll;\n\
  int gll_other;\n\
  int n_loop;\n\
  int ispec_other;\n\
  int ispec_test;\n\
  float x_me;\n\
  float y_me;\n\
  float z_me;\n\
  float x_other;\n\
  float y_other;\n\
  float z_other;\n\
  float center_x;\n\
  float center_y;\n\
  float center_z;\n\
  float vx;\n\
  float vy;\n\
  float vz;\n\
  float alpha;\n\
  float ratio;\n\
  float theta;\n\
  float r0;\n\
  float r1;\n\
  float r0_squared;\n\
  float r1_squared;\n\
  float dist_h;\n\
  float dist_v;\n\
  float val;\n\
  float val_gaussian;\n\
  float coef;\n\
  float normalisation_slice;\n\
  float dat;\n\
\n\
  __local int sh_test[(NGLL3)];\n\
  __local float sh_x_other[(NGLL3)];\n\
  __local float sh_y_other[(NGLL3)];\n\
  __local float sh_z_other[(NGLL3)];\n\
  __local float sh_integ_factor[(NGLL3)];\n\
  __local float sh_data[(NGLL3)];\n\
\n\
  const float PI2 = 9.869604401089358f;\n\
\n\
  ispec = get_group_id(0) + (get_group_id(1)) * (get_num_groups(0));\n\
  igll = get_local_id(0);\n\
\n\
  n_loop = (nspec_other) / (NGLL3) + 1;\n\
  x_me = xstore_me[(NGLL3) * (ispec) + igll];\n\
  y_me = ystore_me[(NGLL3) * (ispec) + igll];\n\
  z_me = zstore_me[(NGLL3) * (ispec) + igll];\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
  dat = 0.0f;\n\
  normalisation_slice = 0.0f;\n\
\n\
  for (int i = 0; i < n_loop; i += 1) {\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
    ispec_other = (NGLL3) * (i) + igll;\n\
    if (ispec_other < nspec_other) {\n\
      center_x = (xstore_other[(ispec_other) * (NGLL3)] + xstore_other[(ispec_other) * (NGLL3) + NGLL3 - (1)]) * (0.5f);\n\
      center_y = (ystore_other[(ispec_other) * (NGLL3)] + ystore_other[(ispec_other) * (NGLL3) + NGLL3 - (1)]) * (0.5f);\n\
      center_z = (zstore_other[(ispec_other) * (NGLL3)] + zstore_other[(ispec_other) * (NGLL3) + NGLL3 - (1)]) * (0.5f);\n\
\n\
      r0_squared = (x_me) * (x_me) + (y_me) * (y_me) + (z_me) * (z_me);\n\
      r1_squared = (center_x) * (center_x) + (center_y) * (center_y) + (center_z) * (center_z);\n\
\n\
      if (use_vector_distance) {\n\
        r0 = sqrt(r0_squared);\n\
        r1 = sqrt(r1_squared);\n\
        dist_v = (r1 - (r0)) * (r1 - (r0));\n\
\n\
        alpha = (r1) / (r0);\n\
        vx = (alpha) * (x_me);\n\
        vy = (alpha) * (y_me);\n\
        vz = (alpha) * (z_me);\n\
        vx = center_x - (vx);\n\
        vy = center_y - (vy);\n\
        vz = center_z - (vz);\n\
        dist_h = (vx) * (vx) + (vy) * (vy) + (vz) * (vz);\n\
      } else {\n\
        alpha = sqrt((r0_squared) * (r1_squared));\n\
        dist_v = r1_squared + r0_squared - ((2.0f) * (alpha));\n\
\n\
        if (alpha > 0.0f) {\n\
          ratio = ((x_me) * (center_x) + (y_me) * (center_y) + (z_me) * (center_z)) / (alpha);\n\
        } else {\n\
          ratio = 1.0f;\n\
        }\n\
        if (ratio >= 1.0f) {\n\
          dist_h = 0.0f;\n\
        } else if (ratio <= -1.0f) {\n\
          dist_h = (r1_squared) * (PI2);\n\
        } else {\n\
          theta = acos(ratio);\n\
          dist_h = (r1_squared) * ((theta) * (theta));\n\
        }\n\
      }\n\
    }\n\
\n\
    sh_test[igll] = (ispec_other >= nspec_other ? 1 : 0);\n\
    sh_test[igll] = (dist_h > h_criterion || sh_test[igll] ? 1 : 0);\n\
    sh_test[igll] = (dist_v > v_criterion || sh_test[igll] ? 1 : 0);\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
    for (int k = 0; k < NGLL3; k += 1) {\n\
      barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
      if (sh_test[k]) {\n\
        continue;\n\
      }\n\
\n\
      ispec_test = (i) * (NGLL3) + k;\n\
      sh_x_other[igll] = xstore_other[(ispec_test) * (NGLL3) + igll];\n\
      sh_y_other[igll] = ystore_other[(ispec_test) * (NGLL3) + igll];\n\
      sh_z_other[igll] = zstore_other[(ispec_test) * (NGLL3) + igll];\n\
\n\
      sh_data[igll] = data_other[(ispec_test) * (NGLL3) + igll];\n\
      sh_integ_factor[igll] = integ_factor[(ispec_test) * (NGLL3) + igll];\n\
\n\
      barrier(CLK_LOCAL_MEM_FENCE);\n\
\n\
      for (int j = 0; j < NGLL3; j += 1) {\n\
        gll_other = ((igll + j < 0) ^ (NGLL3 < 0) ? ((igll + j) % (NGLL3)) + NGLL3 : (igll + j) % (NGLL3));\n\
\n\
        x_other = sh_x_other[gll_other];\n\
        y_other = sh_y_other[gll_other];\n\
        z_other = sh_z_other[gll_other];\n\
\n\
        r0_squared = (x_me) * (x_me) + (y_me) * (y_me) + (z_me) * (z_me);\n\
        r1_squared = (x_other) * (x_other) + (y_other) * (y_other) + (z_other) * (z_other);\n\
\n\
        if (use_vector_distance) {\n\
          r0 = sqrt(r0_squared);\n\
          r1 = sqrt(r1_squared);\n\
          dist_v = (r1 - (r0)) * (r1 - (r0));\n\
\n\
          alpha = (r1) / (r0);\n\
          vx = (alpha) * (x_me);\n\
          vy = (alpha) * (y_me);\n\
          vz = (alpha) * (z_me);\n\
          vx = x_other - (vx);\n\
          vy = y_other - (vy);\n\
          vz = z_other - (vz);\n\
          dist_h = (vx) * (vx) + (vy) * (vy) + (vz) * (vz);\n\
        } else {\n\
          alpha = sqrt((r0_squared) * (r1_squared));\n\
          dist_v = r1_squared + r0_squared - ((2.0f) * (alpha));\n\
\n\
          if (alpha > 0.0f) {\n\
            ratio = ((x_me) * (x_other) + (y_me) * (y_other) + (z_me) * (z_other)) / (alpha);\n\
          } else {\n\
            ratio = 1.0f;\n\
          }\n\
          if (ratio >= 1.0f) {\n\
            dist_h = 0.0f;\n\
          } else if (ratio <= -1.0f) {\n\
            dist_h = (r1_squared) * (PI2);\n\
          } else {\n\
            theta = acos(ratio);\n\
            dist_h = (r1_squared) * ((theta) * (theta));\n\
          }\n\
        }\n\
\n\
        val = ( -(dist_h)) * (sigma_h2_inv) - ((dist_v) * (sigma_v2_inv));\n\
        if (val < -86.0f) {\n\
          val_gaussian = 0.0f;\n\
        } else {\n\
          val_gaussian = exp(val);\n\
        }\n\
        coef = (val_gaussian) * (sh_integ_factor[gll_other]);\n\
        normalisation_slice = normalisation_slice + coef;\n\
        dat = dat + (sh_data[gll_other]) * (coef);\n\
      }\n\
    }\n\
  }\n\
  val = data_smooth[((NGLL3) * (nspec_me)) * (iker) + (NGLL3) * (ispec) + igll] + dat;\n\
  data_smooth[((NGLL3) * (nspec_me)) * (iker) + (NGLL3) * (ispec) + igll] = val;\n\
  val = normalisation[(NGLL3) * (ispec) + igll] + normalisation_slice;\n\
  normalisation[(NGLL3) * (ispec) + igll] = val;\n\
}\n\
";
