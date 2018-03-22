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

#ifdef USE_TEXTURES_CONSTANTS
#undef USE_TEXTURES_CONSTANTS
#endif

// fluid rotation

static __device__ void compute_element_oc_rotation(const int tx, const int working_element, const float time, const float two_omega_earth, const float deltat, float * d_A_array_rotation, float * d_B_array_rotation, const float dpotentialdxl, const float dpotentialdyl, float * dpotentialdx_with_rot, float * dpotentialdy_with_rot){
  float two_omega_deltat;
  float cos_two_omega_t;
  float sin_two_omega_t;
  float A_rotation;
  float B_rotation;
  float source_euler_A;
  float source_euler_B;

  // store the source for the Euler scheme for A_rotation and B_rotation
  sincosf((two_omega_earth) * (time),  &sin_two_omega_t,  &cos_two_omega_t);

  // time step deltat of Euler scheme is included in the source
  two_omega_deltat = (deltat) * (two_omega_earth);
  source_euler_A = (two_omega_deltat) * ((cos_two_omega_t) * (dpotentialdyl) + (sin_two_omega_t) * (dpotentialdxl));
  source_euler_B = (two_omega_deltat) * ((sin_two_omega_t) * (dpotentialdyl) - ((cos_two_omega_t) * (dpotentialdxl)));
  A_rotation = d_A_array_rotation[tx + (working_element) * (NGLL3)];
  B_rotation = d_B_array_rotation[tx + (working_element) * (NGLL3)];
  dpotentialdx_with_rot[0] = dpotentialdxl + (A_rotation) * (cos_two_omega_t) + (B_rotation) * (sin_two_omega_t);
  dpotentialdy_with_rot[0] = dpotentialdyl + ( -(A_rotation)) * (sin_two_omega_t) + (B_rotation) * (cos_two_omega_t);

  // updates rotation term with Euler scheme (non-padded offset)
  d_A_array_rotation[tx + (working_element) * (NGLL3)] = d_A_array_rotation[tx + (working_element) * (NGLL3)] + source_euler_A;
  d_B_array_rotation[tx + (working_element) * (NGLL3)] = d_B_array_rotation[tx + (working_element) * (NGLL3)] + source_euler_B;
}

// KERNEL 2
//
// for outer core ( acoustic domain )

__global__ void outer_core_impl_kernel_forward(const int nb_blocks_to_compute, const int * d_ibool, const int * d_phase_ispec_inner, const int num_phase_ispec, const int d_iphase, const int use_mesh_coloring_gpu, const float * __restrict__ d_potential, float * d_potential_dot_dot, const float * __restrict__ d_xix, const float * __restrict__ d_xiy, const float * __restrict__ d_xiz, const float * __restrict__ d_etax, const float * __restrict__ d_etay, const float * __restrict__ d_etaz, const float * __restrict__ d_gammax, const float * __restrict__ d_gammay, const float * __restrict__ d_gammaz, const float * __restrict__ d_hprime_xx, const float * __restrict__ d_hprimewgll_xx, const float * __restrict__ wgllwgll_xy, const float * __restrict__ wgllwgll_xz, const float * __restrict__ wgllwgll_yz, const int GRAVITY, const float * __restrict__ d_rstore, const float * __restrict__ d_d_ln_density_dr_table, const float * __restrict__ d_minus_rho_g_over_kappa_fluid, const float * __restrict__ wgll_cube, const int ROTATION, const float time, const float two_omega_earth, const float deltat, float * d_A_array_rotation, float * d_B_array_rotation, const int NSPEC_OUTER_CORE){
  int bx;
  int tx;
  int K;
  int J;
  int I;
#ifndef MANUALLY_UNROLLED_LOOPS
  int l;
#endif
  unsigned short active_1;
  int offset;
  int iglob_1;
  int working_element;
  float temp1l;
  float temp2l;
  float temp3l;
  float xixl;
  float xiyl;
  float xizl;
  float etaxl;
  float etayl;
  float etazl;
  float gammaxl;
  float gammayl;
  float gammazl;
  float jacobianl;
  float dpotentialdxl;
  float dpotentialdyl;
  float dpotentialdzl;
  float dpotentialdx_with_rot;
  float dpotentialdy_with_rot;
  float sum_terms;
  float gravity_term_1;
  float gxl;
  float gyl;
  float gzl;
  float radius;
  float theta;
  float phi;
  float cos_theta;
  float sin_theta;
  float cos_phi;
  float sin_phi;
  float grad_x_ln_rho;
  float grad_y_ln_rho;
  float grad_z_ln_rho;
  int int_radius;
  __shared__ float s_dummy_loc[(NGLL3)];
  __shared__ float s_temp1[(NGLL3)];
  __shared__ float s_temp2[(NGLL3)];
  __shared__ float s_temp3[(NGLL3)];
  __shared__ float sh_hprime_xx[(NGLL2)];
  __shared__ float sh_hprimewgll_xx[(NGLL2)];
  bx = (blockIdx.y) * (gridDim.x) + blockIdx.x;
  tx = threadIdx.x + ((NGLL3_PADDED) * (0)) / (1);

  // use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  // because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses
  active_1 = (tx < NGLL3 && bx < nb_blocks_to_compute ? 1 : 0);

  // copy from global memory to shared memory
  // each thread writes one of the NGLL^3 = 125 data points
  if (active_1) {
#ifdef USE_MESH_COLORING_GPU
    working_element = bx;
#else
    if (use_mesh_coloring_gpu) {
      working_element = bx;
    } else {
      working_element = d_phase_ispec_inner[bx + (num_phase_ispec) * (d_iphase - (1))] - (1);
    }
#endif
    iglob_1 = d_ibool[(working_element) * (NGLL3) + tx] - (1);
#ifdef USE_TEXTURES_FIELDS
    s_dummy_loc[tx] = tex1Dfetch(d_displ_oc_tex,iglob_1);
#else
    s_dummy_loc[tx] = d_potential[iglob_1];
#endif
  }
  if (tx < NGLL2) {
#ifdef USE_TEXTURES_CONSTANTS
    sh_hprime_xx[tx] = tex1Dfetch(d_hprime_xx_oc_tex,tx);
    sh_hprimewgll_xx[tx] = tex1Dfetch(d_hprimewgll_xx_oc_tex,tx);
#else
    sh_hprime_xx[tx] = d_hprime_xx[tx];
    sh_hprimewgll_xx[tx] = d_hprimewgll_xx[tx];
#endif
  }
  __syncthreads();
  K = (tx) / (NGLL2);
  J = (tx - ((K) * (NGLL2))) / (NGLLX);
  I = tx - ((K) * (NGLL2)) - ((J) * (NGLLX));
  if (active_1) {
    temp1l = 0.0f;
    temp2l = 0.0f;
    temp3l = 0.0f;
#ifdef MANUALLY_UNROLLED_LOOPS
    temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 0]) * (sh_hprime_xx[(0) * (NGLLX) + I]);
    temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (0) * (NGLLX) + I]) * (sh_hprime_xx[(0) * (NGLLX) + J]);
    temp3l = temp3l + (s_dummy_loc[(0) * (NGLL2) + (J) * (NGLLX) + I]) * (sh_hprime_xx[(0) * (NGLLX) + K]);
    temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 1]) * (sh_hprime_xx[(1) * (NGLLX) + I]);
    temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (1) * (NGLLX) + I]) * (sh_hprime_xx[(1) * (NGLLX) + J]);
    temp3l = temp3l + (s_dummy_loc[(1) * (NGLL2) + (J) * (NGLLX) + I]) * (sh_hprime_xx[(1) * (NGLLX) + K]);
    temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 2]) * (sh_hprime_xx[(2) * (NGLLX) + I]);
    temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (2) * (NGLLX) + I]) * (sh_hprime_xx[(2) * (NGLLX) + J]);
    temp3l = temp3l + (s_dummy_loc[(2) * (NGLL2) + (J) * (NGLLX) + I]) * (sh_hprime_xx[(2) * (NGLLX) + K]);
    temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 3]) * (sh_hprime_xx[(3) * (NGLLX) + I]);
    temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (3) * (NGLLX) + I]) * (sh_hprime_xx[(3) * (NGLLX) + J]);
    temp3l = temp3l + (s_dummy_loc[(3) * (NGLL2) + (J) * (NGLLX) + I]) * (sh_hprime_xx[(3) * (NGLLX) + K]);
    temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 4]) * (sh_hprime_xx[(4) * (NGLLX) + I]);
    temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (4) * (NGLLX) + I]) * (sh_hprime_xx[(4) * (NGLLX) + J]);
    temp3l = temp3l + (s_dummy_loc[(4) * (NGLL2) + (J) * (NGLLX) + I]) * (sh_hprime_xx[(4) * (NGLLX) + K]);
#else
    for (l = 0; l <= NGLLX - (1); l += 1) {
      temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + l]) * (sh_hprime_xx[(l) * (NGLLX) + I]);
      temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (l) * (NGLLX) + I]) * (sh_hprime_xx[(l) * (NGLLX) + J]);
      temp3l = temp3l + (s_dummy_loc[(l) * (NGLL2) + (J) * (NGLLX) + I]) * (sh_hprime_xx[(l) * (NGLLX) + K]);
    }
#endif
    offset = (working_element) * (NGLL3_PADDED) + tx;
    xixl = d_xix[offset];
    etaxl = d_etax[offset];
    gammaxl = d_gammax[offset];
    xiyl = d_xiy[offset];
    etayl = d_etay[offset];
    gammayl = d_gammay[offset];
    xizl = d_xiz[offset];
    etazl = d_etaz[offset];
    gammazl = d_gammaz[offset];
    jacobianl = (1.0f) / ((xixl) * ((etayl) * (gammazl) - ((etazl) * (gammayl))) - ((xiyl) * ((etaxl) * (gammazl) - ((etazl) * (gammaxl)))) + (xizl) * ((etaxl) * (gammayl) - ((etayl) * (gammaxl))));
    dpotentialdxl = (xixl) * (temp1l) + (etaxl) * (temp2l) + (gammaxl) * (temp3l);
    dpotentialdyl = (xiyl) * (temp1l) + (etayl) * (temp2l) + (gammayl) * (temp3l);
    dpotentialdzl = (xizl) * (temp1l) + (etazl) * (temp2l) + (gammazl) * (temp3l);
    if (ROTATION) {
      compute_element_oc_rotation(tx, working_element, time, two_omega_earth, deltat, d_A_array_rotation, d_B_array_rotation, dpotentialdxl, dpotentialdyl,  &dpotentialdx_with_rot,  &dpotentialdy_with_rot);
    } else {
      dpotentialdx_with_rot = dpotentialdxl;
      dpotentialdy_with_rot = dpotentialdyl;
    }
    radius = d_rstore[0 + (3) * (iglob_1)];
    theta = d_rstore[1 + (3) * (iglob_1)];
    phi = d_rstore[2 + (3) * (iglob_1)];
    sincosf(theta,  &sin_theta,  &cos_theta);
    sincosf(phi,  &sin_phi,  &cos_phi);
    int_radius = rint(((radius) * (R_EARTH_KM)) * (10.0f)) - (1);
    if ( !(GRAVITY)) {
      grad_x_ln_rho = ((sin_theta) * (cos_phi)) * (d_d_ln_density_dr_table[int_radius]);
      grad_y_ln_rho = ((sin_theta) * (sin_phi)) * (d_d_ln_density_dr_table[int_radius]);
      grad_z_ln_rho = (cos_theta) * (d_d_ln_density_dr_table[int_radius]);
      dpotentialdx_with_rot = dpotentialdx_with_rot + (s_dummy_loc[tx]) * (grad_x_ln_rho);
      dpotentialdy_with_rot = dpotentialdy_with_rot + (s_dummy_loc[tx]) * (grad_y_ln_rho);
      dpotentialdzl = dpotentialdzl + (s_dummy_loc[tx]) * (grad_z_ln_rho);
    } else {
      gxl = (sin_theta) * (cos_phi);
      gyl = (sin_theta) * (sin_phi);
      gzl = cos_theta;
      gravity_term_1 = (((d_minus_rho_g_over_kappa_fluid[int_radius]) * (jacobianl)) * (wgll_cube[tx])) * ((dpotentialdx_with_rot) * (gxl) + (dpotentialdy_with_rot) * (gyl) + (dpotentialdzl) * (gzl));
    }
    s_temp1[tx] = (jacobianl) * ((xixl) * (dpotentialdx_with_rot) + (xiyl) * (dpotentialdy_with_rot) + (xizl) * (dpotentialdzl));
    s_temp2[tx] = (jacobianl) * ((etaxl) * (dpotentialdx_with_rot) + (etayl) * (dpotentialdy_with_rot) + (etazl) * (dpotentialdzl));
    s_temp3[tx] = (jacobianl) * ((gammaxl) * (dpotentialdx_with_rot) + (gammayl) * (dpotentialdy_with_rot) + (gammazl) * (dpotentialdzl));
  }
  __syncthreads();
  if (active_1) {
    temp1l = 0.0f;
    temp2l = 0.0f;
    temp3l = 0.0f;
#ifdef MANUALLY_UNROLLED_LOOPS
    temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + 0]) * (sh_hprimewgll_xx[(I) * (NGLLX) + 0]);
    temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (0) * (NGLLX) + I]) * (sh_hprimewgll_xx[(J) * (NGLLX) + 0]);
    temp3l = temp3l + (s_temp3[(0) * (NGLL2) + (J) * (NGLLX) + I]) * (sh_hprimewgll_xx[(K) * (NGLLX) + 0]);
    temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + 1]) * (sh_hprimewgll_xx[(I) * (NGLLX) + 1]);
    temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (1) * (NGLLX) + I]) * (sh_hprimewgll_xx[(J) * (NGLLX) + 1]);
    temp3l = temp3l + (s_temp3[(1) * (NGLL2) + (J) * (NGLLX) + I]) * (sh_hprimewgll_xx[(K) * (NGLLX) + 1]);
    temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + 2]) * (sh_hprimewgll_xx[(I) * (NGLLX) + 2]);
    temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (2) * (NGLLX) + I]) * (sh_hprimewgll_xx[(J) * (NGLLX) + 2]);
    temp3l = temp3l + (s_temp3[(2) * (NGLL2) + (J) * (NGLLX) + I]) * (sh_hprimewgll_xx[(K) * (NGLLX) + 2]);
    temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + 3]) * (sh_hprimewgll_xx[(I) * (NGLLX) + 3]);
    temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (3) * (NGLLX) + I]) * (sh_hprimewgll_xx[(J) * (NGLLX) + 3]);
    temp3l = temp3l + (s_temp3[(3) * (NGLL2) + (J) * (NGLLX) + I]) * (sh_hprimewgll_xx[(K) * (NGLLX) + 3]);
    temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + 4]) * (sh_hprimewgll_xx[(I) * (NGLLX) + 4]);
    temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (4) * (NGLLX) + I]) * (sh_hprimewgll_xx[(J) * (NGLLX) + 4]);
    temp3l = temp3l + (s_temp3[(4) * (NGLL2) + (J) * (NGLLX) + I]) * (sh_hprimewgll_xx[(K) * (NGLLX) + 4]);
#else
    for (l = 0; l <= NGLLX - (1); l += 1) {
      temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + l]) * (sh_hprimewgll_xx[(I) * (NGLLX) + l]);
      temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (l) * (NGLLX) + I]) * (sh_hprimewgll_xx[(J) * (NGLLX) + l]);
      temp3l = temp3l + (s_temp3[(l) * (NGLL2) + (J) * (NGLLX) + I]) * (sh_hprimewgll_xx[(K) * (NGLLX) + l]);
    }
#endif
    sum_terms =  -((wgllwgll_yz[(K) * (NGLLX) + J]) * (temp1l) + (wgllwgll_xz[(K) * (NGLLX) + I]) * (temp2l) + (wgllwgll_xy[(J) * (NGLLX) + I]) * (temp3l));
    if (GRAVITY) {
      sum_terms = sum_terms + gravity_term_1;
    }
#ifdef USE_MESH_COLORING_GPU
#ifdef USE_TEXTURES_FIELDS
    d_potential_dot_dot[iglob_1] = tex1Dfetch(d_accel_oc_tex,iglob_1) + sum_terms;
#else
    d_potential_dot_dot[iglob_1] = d_potential_dot_dot[iglob_1] + sum_terms;
#endif
#else
    if (use_mesh_coloring_gpu) {
      if (NSPEC_OUTER_CORE > 1000) {
#ifdef USE_TEXTURES_FIELDS
        d_potential_dot_dot[iglob_1] = tex1Dfetch(d_accel_oc_tex,iglob_1) + sum_terms;
#else
        d_potential_dot_dot[iglob_1] = d_potential_dot_dot[iglob_1] + sum_terms;
#endif
      } else {
        atomicAdd(d_potential_dot_dot + iglob_1, sum_terms);
      }
    } else {
      atomicAdd(d_potential_dot_dot + iglob_1, sum_terms);
    }
#endif
  }
}
