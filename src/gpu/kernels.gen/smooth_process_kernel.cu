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

__global__ void smooth_process_kernel(const float * xstore_me, const float * ystore_me, const float * zstore_me, const float * xstore_other, const float * ystore_other, const float * zstore_other, const float * data_other, const float sigma_h2_inv, const float sigma_v2_inv, const int iker, const int nspec_me, const int nspec_other, const float v_criterion, const float h_criterion, const float * integ_factor, float * data_smooth, float * normalisation, const int use_vector_distance){
  int ispec;
  int igll;
  int gll_other;
  int n_loop;
  int ispec_other;
  int ispec_test;
  float x_me;
  float y_me;
  float z_me;
  float x_other;
  float y_other;
  float z_other;
  float center_x;
  float center_y;
  float center_z;
  float vx;
  float vy;
  float vz;
  float alpha;
  float ratio;
  float theta;
  float r0;
  float r1;
  float r0_squared;
  float r1_squared;
  float dist_h;
  float dist_v;
  float val;
  float val_gaussian;
  float coef;
  float normalisation_slice;
  float dat;

  __shared__ int sh_test[(NGLL3)];
  __shared__ float sh_x_other[(NGLL3)];
  __shared__ float sh_y_other[(NGLL3)];
  __shared__ float sh_z_other[(NGLL3)];
  __shared__ float sh_integ_factor[(NGLL3)];
  __shared__ float sh_data[(NGLL3)];

  const float PI2 = 9.869604401089358f;

  ispec = blockIdx.x + (blockIdx.y) * (gridDim.x);
  igll = threadIdx.x;

  n_loop = (nspec_other) / (NGLL3) + 1;
  x_me = xstore_me[(NGLL3) * (ispec) + igll];
  y_me = ystore_me[(NGLL3) * (ispec) + igll];
  z_me = zstore_me[(NGLL3) * (ispec) + igll];
  __syncthreads();

  dat = 0.0f;
  normalisation_slice = 0.0f;

  for (int i = 0; i < n_loop; i += 1) {
    __syncthreads();

    ispec_other = (NGLL3) * (i) + igll;
    if (ispec_other < nspec_other) {
      center_x = (xstore_other[(ispec_other) * (NGLL3)] + xstore_other[(ispec_other) * (NGLL3) + NGLL3 - (1)]) * (0.5f);
      center_y = (ystore_other[(ispec_other) * (NGLL3)] + ystore_other[(ispec_other) * (NGLL3) + NGLL3 - (1)]) * (0.5f);
      center_z = (zstore_other[(ispec_other) * (NGLL3)] + zstore_other[(ispec_other) * (NGLL3) + NGLL3 - (1)]) * (0.5f);

      r0_squared = (x_me) * (x_me) + (y_me) * (y_me) + (z_me) * (z_me);
      r1_squared = (center_x) * (center_x) + (center_y) * (center_y) + (center_z) * (center_z);

      if (use_vector_distance) {
        r0 = sqrt(r0_squared);
        r1 = sqrt(r1_squared);
        dist_v = (r1 - (r0)) * (r1 - (r0));

        alpha = (r1) / (r0);
        vx = (alpha) * (x_me);
        vy = (alpha) * (y_me);
        vz = (alpha) * (z_me);
        vx = center_x - (vx);
        vy = center_y - (vy);
        vz = center_z - (vz);
        dist_h = (vx) * (vx) + (vy) * (vy) + (vz) * (vz);
      } else {
        alpha = sqrt((r0_squared) * (r1_squared));
        dist_v = r1_squared + r0_squared - ((2.0f) * (alpha));

        if (alpha > 0.0f) {
          ratio = ((x_me) * (center_x) + (y_me) * (center_y) + (z_me) * (center_z)) / (alpha);
        } else {
          ratio = 1.0f;
        }
        if (ratio >= 1.0f) {
          dist_h = 0.0f;
        } else if (ratio <= -1.0f) {
          dist_h = (r1_squared) * (PI2);
        } else {
          theta = acos(ratio);
          dist_h = (r1_squared) * ((theta) * (theta));
        }
      }
    }

    sh_test[igll] = (ispec_other >= nspec_other ? 1 : 0);
    sh_test[igll] = (dist_h > h_criterion || sh_test[igll] ? 1 : 0);
    sh_test[igll] = (dist_v > v_criterion || sh_test[igll] ? 1 : 0);
    __syncthreads();

    for (int k = 0; k < NGLL3; k += 1) {
      __syncthreads();

      if (sh_test[k]) {
        continue;
      }

      ispec_test = (i) * (NGLL3) + k;
      sh_x_other[igll] = xstore_other[(ispec_test) * (NGLL3) + igll];
      sh_y_other[igll] = ystore_other[(ispec_test) * (NGLL3) + igll];
      sh_z_other[igll] = zstore_other[(ispec_test) * (NGLL3) + igll];

      sh_data[igll] = data_other[(ispec_test) * (NGLL3) + igll];
      sh_integ_factor[igll] = integ_factor[(ispec_test) * (NGLL3) + igll];

      __syncthreads();

      for (int j = 0; j < NGLL3; j += 1) {
        gll_other = ((igll + j < 0) ^ (NGLL3 < 0) ? ((igll + j) % (NGLL3)) + NGLL3 : (igll + j) % (NGLL3));

        x_other = sh_x_other[gll_other];
        y_other = sh_y_other[gll_other];
        z_other = sh_z_other[gll_other];

        r0_squared = (x_me) * (x_me) + (y_me) * (y_me) + (z_me) * (z_me);
        r1_squared = (x_other) * (x_other) + (y_other) * (y_other) + (z_other) * (z_other);

        if (use_vector_distance) {
          r0 = sqrt(r0_squared);
          r1 = sqrt(r1_squared);
          dist_v = (r1 - (r0)) * (r1 - (r0));

          alpha = (r1) / (r0);
          vx = (alpha) * (x_me);
          vy = (alpha) * (y_me);
          vz = (alpha) * (z_me);
          vx = x_other - (vx);
          vy = y_other - (vy);
          vz = z_other - (vz);
          dist_h = (vx) * (vx) + (vy) * (vy) + (vz) * (vz);
        } else {
          alpha = sqrt((r0_squared) * (r1_squared));
          dist_v = r1_squared + r0_squared - ((2.0f) * (alpha));

          if (alpha > 0.0f) {
            ratio = ((x_me) * (x_other) + (y_me) * (y_other) + (z_me) * (z_other)) / (alpha);
          } else {
            ratio = 1.0f;
          }
          if (ratio >= 1.0f) {
            dist_h = 0.0f;
          } else if (ratio <= -1.0f) {
            dist_h = (r1_squared) * (PI2);
          } else {
            theta = acos(ratio);
            dist_h = (r1_squared) * ((theta) * (theta));
          }
        }

        val = ( -(dist_h)) * (sigma_h2_inv) - ((dist_v) * (sigma_v2_inv));
        if (val < -86.0f) {
          val_gaussian = 0.0f;
        } else {
          val_gaussian = exp(val);
        }
        coef = (val_gaussian) * (sh_integ_factor[gll_other]);
        normalisation_slice = normalisation_slice + coef;
        dat = dat + (sh_data[gll_other]) * (coef);
      }
    }
  }
  val = data_smooth[((NGLL3) * (nspec_me)) * (iker) + (NGLL3) * (ispec) + igll] + dat;
  data_smooth[((NGLL3) * (nspec_me)) * (iker) + (NGLL3) * (ispec) + igll] = val;
  val = normalisation[(NGLL3) * (ispec) + igll] + normalisation_slice;
  normalisation[(NGLL3) * (ispec) + igll] = val;
}
