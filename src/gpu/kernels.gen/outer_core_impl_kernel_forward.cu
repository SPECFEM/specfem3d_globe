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
static __device__ void compute_element_oc_rotation(const int tx, const int working_element, const float time, const float two_omega_earth, const float deltat, float * d_A_array_rotation, float * d_B_array_rotation, const float dpotentialdxl, const float dpotentialdyl, float * dpotentialdx_with_rot, float * dpotentialdy_with_rot){
  float two_omega_deltat;
  float cos_two_omega_t;
  float sin_two_omega_t;
  float A_rotation;
  float B_rotation;
  float source_euler_A;
  float source_euler_B;
  sincosf((two_omega_earth) * (time),  &sin_two_omega_t,  &cos_two_omega_t);
  two_omega_deltat = (deltat) * (two_omega_earth);
  source_euler_A = (two_omega_deltat) * ((cos_two_omega_t) * (dpotentialdyl) + (sin_two_omega_t) * (dpotentialdxl));
  source_euler_B = (two_omega_deltat) * ((sin_two_omega_t) * (dpotentialdyl) - ((cos_two_omega_t) * (dpotentialdxl)));
  A_rotation = d_A_array_rotation[tx + (working_element) * (NGLL3) - 0];
  B_rotation = d_B_array_rotation[tx + (working_element) * (NGLL3) - 0];
  dpotentialdx_with_rot[0 - 0] = dpotentialdxl + (A_rotation) * (cos_two_omega_t) + (B_rotation) * (sin_two_omega_t);
  dpotentialdy_with_rot[0 - 0] = dpotentialdyl + ( -(A_rotation)) * (sin_two_omega_t) + (B_rotation) * (cos_two_omega_t);
  d_A_array_rotation[tx + (working_element) * (NGLL3) - 0] = d_A_array_rotation[tx + (working_element) * (NGLL3) - 0] + source_euler_A;
  d_B_array_rotation[tx + (working_element) * (NGLL3) - 0] = d_B_array_rotation[tx + (working_element) * (NGLL3) - 0] + source_euler_B;
}
__global__ void outer_core_impl_kernel_forward(const int nb_blocks_to_compute, const int * d_ibool, const int * d_phase_ispec_inner, const int num_phase_ispec, const int d_iphase, const int use_mesh_coloring_gpu, const float * __restrict__ d_potential, float * d_potential_dot_dot, const float * __restrict__ d_xix, const float * __restrict__ d_xiy, const float * __restrict__ d_xiz, const float * __restrict__ d_etax, const float * __restrict__ d_etay, const float * __restrict__ d_etaz, const float * __restrict__ d_gammax, const float * __restrict__ d_gammay, const float * __restrict__ d_gammaz, const float * __restrict__ d_hprime_xx, const float * __restrict__ d_hprimewgll_xx, const float * __restrict__ wgllwgll_xy, const float * __restrict__ wgllwgll_xz, const float * __restrict__ wgllwgll_yz, const int GRAVITY, const float * __restrict__ d_xstore, const float * __restrict__ d_ystore, const float * __restrict__ d_zstore, const float * __restrict__ d_d_ln_density_dr_table, const float * __restrict__ d_minus_rho_g_over_kappa_fluid, const float * __restrict__ wgll_cube, const int ROTATION, const float time, const float two_omega_earth, const float deltat, float * d_A_array_rotation, float * d_B_array_rotation, const int NSPEC_OUTER_CORE){
  int bx;
  int tx;
  int K;
  int J;
  int I;
  unsigned short active;
  int offset;
  int iglob;
  int working_element;
  int l;
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
  float gravity_term;
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
  __shared__ float s_dummy_loc[NGLL3 + 0 - (1) - (0) + 1];
  __shared__ float s_temp1[NGLL3 + 0 - (1) - (0) + 1];
  __shared__ float s_temp2[NGLL3 + 0 - (1) - (0) + 1];
  __shared__ float s_temp3[NGLL3 + 0 - (1) - (0) + 1];
  __shared__ float sh_hprime_xx[NGLL2 + 0 - (1) - (0) + 1];
  __shared__ float sh_hprimewgll_xx[NGLL2 + 0 - (1) - (0) + 1];
  bx = (blockIdx.y) * (gridDim.x) + blockIdx.x;
  tx = threadIdx.x;
  K = (tx) / (NGLL2);
  J = (tx - ((K) * (NGLL2))) / (NGLLX);
  I = tx - ((K) * (NGLL2)) - ((J) * (NGLLX));
  active = (tx < NGLL3 && bx < nb_blocks_to_compute ? 1 : 0);
  if(active){
#ifdef USE_MESH_COLORING_GPU
    working_element = bx;
#else
    if(use_mesh_coloring_gpu){
      working_element = bx;
    } else {
      working_element = d_phase_ispec_inner[bx + (num_phase_ispec) * (d_iphase - (1)) - 0] - (1);
    }
#endif
    iglob = d_ibool[(working_element) * (NGLL3) + tx - 0] - (1);
#ifdef USE_TEXTURES_FIELDS
    s_dummy_loc[tx - 0] = tex1Dfetch(d_displ_oc_tex,iglob);
#else
    s_dummy_loc[tx - 0] = d_potential[iglob - 0];
#endif
  }
  if(tx < NGLL2){
#ifdef USE_TEXTURES_CONSTANTS
    sh_hprime_xx[tx - 0] = tex1Dfetch(d_hprime_xx_oc_tex,tx);
    sh_hprimewgll_xx[tx - 0] = tex1Dfetch(d_hprimewgll_xx_oc_tex,tx);
#else
    sh_hprime_xx[tx - 0] = d_hprime_xx[tx - 0];
    sh_hprimewgll_xx[tx - 0] = d_hprimewgll_xx[tx - 0];
#endif
  }
  __syncthreads();
  if(active){
    temp1l = 0.0f;
    temp2l = 0.0f;
    temp3l = 0.0f;
#ifdef MANUALLY_UNROLLED_LOOPS
    temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 0 - 0]) * (sh_hprime_xx[(0) * (NGLLX) + I - 0]);
    temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (0) * (NGLLX) + I - 0]) * (sh_hprime_xx[(0) * (NGLLX) + J - 0]);
    temp3l = temp3l + (s_dummy_loc[(0) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (sh_hprime_xx[(0) * (NGLLX) + K - 0]);
    temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 1 - 0]) * (sh_hprime_xx[(1) * (NGLLX) + I - 0]);
    temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (1) * (NGLLX) + I - 0]) * (sh_hprime_xx[(1) * (NGLLX) + J - 0]);
    temp3l = temp3l + (s_dummy_loc[(1) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (sh_hprime_xx[(1) * (NGLLX) + K - 0]);
    temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 2 - 0]) * (sh_hprime_xx[(2) * (NGLLX) + I - 0]);
    temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (2) * (NGLLX) + I - 0]) * (sh_hprime_xx[(2) * (NGLLX) + J - 0]);
    temp3l = temp3l + (s_dummy_loc[(2) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (sh_hprime_xx[(2) * (NGLLX) + K - 0]);
    temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 3 - 0]) * (sh_hprime_xx[(3) * (NGLLX) + I - 0]);
    temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (3) * (NGLLX) + I - 0]) * (sh_hprime_xx[(3) * (NGLLX) + J - 0]);
    temp3l = temp3l + (s_dummy_loc[(3) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (sh_hprime_xx[(3) * (NGLLX) + K - 0]);
    temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 4 - 0]) * (sh_hprime_xx[(4) * (NGLLX) + I - 0]);
    temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (4) * (NGLLX) + I - 0]) * (sh_hprime_xx[(4) * (NGLLX) + J - 0]);
    temp3l = temp3l + (s_dummy_loc[(4) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (sh_hprime_xx[(4) * (NGLLX) + K - 0]);
#else
    for(l=0; l<=NGLLX - (1); l+=1){
      temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + l - 0]) * (sh_hprime_xx[(l) * (NGLLX) + I - 0]);
      temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (l) * (NGLLX) + I - 0]) * (sh_hprime_xx[(l) * (NGLLX) + J - 0]);
      temp3l = temp3l + (s_dummy_loc[(l) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (sh_hprime_xx[(l) * (NGLLX) + K - 0]);
    }
#endif
    offset = (working_element) * (NGLL3_PADDED) + tx;
    xixl = d_xix[offset - 0];
    etaxl = d_etax[offset - 0];
    gammaxl = d_gammax[offset - 0];
    xiyl = d_xiy[offset - 0];
    etayl = d_etay[offset - 0];
    gammayl = d_gammay[offset - 0];
    xizl = d_xiz[offset - 0];
    etazl = d_etaz[offset - 0];
    gammazl = d_gammaz[offset - 0];
    jacobianl = (1.0f) / ((xixl) * ((etayl) * (gammazl) - ((etazl) * (gammayl))) - ((xiyl) * ((etaxl) * (gammazl) - ((etazl) * (gammaxl)))) + (xizl) * ((etaxl) * (gammayl) - ((etayl) * (gammaxl))));
    dpotentialdxl = (xixl) * (temp1l) + (etaxl) * (temp2l) + (gammaxl) * (temp3l);
    dpotentialdyl = (xiyl) * (temp1l) + (etayl) * (temp2l) + (gammayl) * (temp3l);
    dpotentialdzl = (xizl) * (temp1l) + (etazl) * (temp2l) + (gammazl) * (temp3l);
    if(ROTATION){
      compute_element_oc_rotation(tx, working_element, time, two_omega_earth, deltat, d_A_array_rotation, d_B_array_rotation, dpotentialdxl, dpotentialdyl,  &dpotentialdx_with_rot,  &dpotentialdy_with_rot);
    } else {
      dpotentialdx_with_rot = dpotentialdxl;
      dpotentialdy_with_rot = dpotentialdyl;
    }
    radius = d_xstore[iglob - 0];
    theta = d_ystore[iglob - 0];
    phi = d_zstore[iglob - 0];
    sincosf(theta,  &sin_theta,  &cos_theta);
    sincosf(phi,  &sin_phi,  &cos_phi);
    int_radius = rint(((radius) * (R_EARTH_KM)) * (10.0f)) - (1);
    if( ! GRAVITY){
      grad_x_ln_rho = ((sin_theta) * (cos_phi)) * (d_d_ln_density_dr_table[int_radius - 0]);
      grad_y_ln_rho = ((sin_theta) * (sin_phi)) * (d_d_ln_density_dr_table[int_radius - 0]);
      grad_z_ln_rho = (cos_theta) * (d_d_ln_density_dr_table[int_radius - 0]);
      dpotentialdx_with_rot = dpotentialdx_with_rot + (s_dummy_loc[tx - 0]) * (grad_x_ln_rho);
      dpotentialdy_with_rot = dpotentialdy_with_rot + (s_dummy_loc[tx - 0]) * (grad_y_ln_rho);
      dpotentialdzl = dpotentialdzl + (s_dummy_loc[tx - 0]) * (grad_z_ln_rho);
    } else {
      gxl = (sin_theta) * (cos_phi);
      gyl = (sin_theta) * (sin_phi);
      gzl = cos_theta;
      gravity_term = (((d_minus_rho_g_over_kappa_fluid[int_radius - 0]) * (jacobianl)) * (wgll_cube[tx - 0])) * ((dpotentialdx_with_rot) * (gxl) + (dpotentialdy_with_rot) * (gyl) + (dpotentialdzl) * (gzl));
    }
    s_temp1[tx - 0] = (jacobianl) * ((xixl) * (dpotentialdx_with_rot) + (xiyl) * (dpotentialdy_with_rot) + (xizl) * (dpotentialdzl));
    s_temp2[tx - 0] = (jacobianl) * ((etaxl) * (dpotentialdx_with_rot) + (etayl) * (dpotentialdy_with_rot) + (etazl) * (dpotentialdzl));
    s_temp3[tx - 0] = (jacobianl) * ((gammaxl) * (dpotentialdx_with_rot) + (gammayl) * (dpotentialdy_with_rot) + (gammazl) * (dpotentialdzl));
  }
  __syncthreads();
  if(active){
    temp1l = 0.0f;
    temp2l = 0.0f;
    temp3l = 0.0f;
#ifdef MANUALLY_UNROLLED_LOOPS
    temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + 0 - 0]) * (sh_hprimewgll_xx[(I) * (NGLLX) + 0 - 0]);
    temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (0) * (NGLLX) + I - 0]) * (sh_hprimewgll_xx[(J) * (NGLLX) + 0 - 0]);
    temp3l = temp3l + (s_temp3[(0) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (sh_hprimewgll_xx[(K) * (NGLLX) + 0 - 0]);
    temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + 1 - 0]) * (sh_hprimewgll_xx[(I) * (NGLLX) + 1 - 0]);
    temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (1) * (NGLLX) + I - 0]) * (sh_hprimewgll_xx[(J) * (NGLLX) + 1 - 0]);
    temp3l = temp3l + (s_temp3[(1) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (sh_hprimewgll_xx[(K) * (NGLLX) + 1 - 0]);
    temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + 2 - 0]) * (sh_hprimewgll_xx[(I) * (NGLLX) + 2 - 0]);
    temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (2) * (NGLLX) + I - 0]) * (sh_hprimewgll_xx[(J) * (NGLLX) + 2 - 0]);
    temp3l = temp3l + (s_temp3[(2) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (sh_hprimewgll_xx[(K) * (NGLLX) + 2 - 0]);
    temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + 3 - 0]) * (sh_hprimewgll_xx[(I) * (NGLLX) + 3 - 0]);
    temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (3) * (NGLLX) + I - 0]) * (sh_hprimewgll_xx[(J) * (NGLLX) + 3 - 0]);
    temp3l = temp3l + (s_temp3[(3) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (sh_hprimewgll_xx[(K) * (NGLLX) + 3 - 0]);
    temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + 4 - 0]) * (sh_hprimewgll_xx[(I) * (NGLLX) + 4 - 0]);
    temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (4) * (NGLLX) + I - 0]) * (sh_hprimewgll_xx[(J) * (NGLLX) + 4 - 0]);
    temp3l = temp3l + (s_temp3[(4) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (sh_hprimewgll_xx[(K) * (NGLLX) + 4 - 0]);
#else
    for(l=0; l<=NGLLX - (1); l+=1){
      temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + l - 0]) * (sh_hprimewgll_xx[(I) * (NGLLX) + l - 0]);
      temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (l) * (NGLLX) + I - 0]) * (sh_hprimewgll_xx[(J) * (NGLLX) + l - 0]);
      temp3l = temp3l + (s_temp3[(l) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (sh_hprimewgll_xx[(K) * (NGLLX) + l - 0]);
    }
#endif
    sum_terms =  -((wgllwgll_yz[(K) * (NGLLX) + J - 0]) * (temp1l) + (wgllwgll_xz[(K) * (NGLLX) + I - 0]) * (temp2l) + (wgllwgll_xy[(J) * (NGLLX) + I - 0]) * (temp3l));
    if(GRAVITY){
      sum_terms = sum_terms + gravity_term;
    }
#ifdef USE_MESH_COLORING_GPU
#ifdef USE_TEXTURES_FIELDS
    d_potential_dot_dot[iglob - 0] = tex1Dfetch(d_accel_oc_tex,iglob) + sum_terms;
#else
    d_potential_dot_dot[iglob - 0] = d_potential_dot_dot[iglob - 0] + sum_terms;
#endif
#else
    if(use_mesh_coloring_gpu){
      if(NSPEC_OUTER_CORE > 1000){
#ifdef USE_TEXTURES_FIELDS
        d_potential_dot_dot[iglob - 0] = tex1Dfetch(d_accel_oc_tex,iglob) + sum_terms;
#else
        d_potential_dot_dot[iglob - 0] = d_potential_dot_dot[iglob - 0] + sum_terms;
#endif
      } else {
        atomicAdd(d_potential_dot_dot + iglob, sum_terms);
      }
    } else {
      atomicAdd(d_potential_dot_dot + iglob, sum_terms);
    }
#endif
  }
}
