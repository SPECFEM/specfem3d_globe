const char * outer_core_impl_kernel_adjoint_program = "\
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
#ifdef USE_TEXTURES_CONSTANTS\n\
#undef USE_TEXTURES_CONSTANTS\n\
#endif\n\
void compute_element_oc_rotation(const int tx, const int working_element, const float time, const float two_omega_earth, const float deltat, __global float * d_A_array_rotation, __global float * d_B_array_rotation, const float dpotentialdxl, const float dpotentialdyl, float * dpotentialdx_with_rot, float * dpotentialdy_with_rot){\n\
  float two_omega_deltat;\n\
  float cos_two_omega_t;\n\
  float sin_two_omega_t;\n\
  float A_rotation;\n\
  float B_rotation;\n\
  float source_euler_A;\n\
  float source_euler_B;\n\
  sin_two_omega_t = sincos((two_omega_earth) * (time),  &cos_two_omega_t);\n\
  two_omega_deltat = (deltat) * (two_omega_earth);\n\
  source_euler_A = (two_omega_deltat) * ((cos_two_omega_t) * (dpotentialdyl) + (sin_two_omega_t) * (dpotentialdxl));\n\
  source_euler_B = (two_omega_deltat) * ((sin_two_omega_t) * (dpotentialdyl) - ((cos_two_omega_t) * (dpotentialdxl)));\n\
  A_rotation = d_A_array_rotation[tx + (working_element) * (NGLL3) - (0)];\n\
  B_rotation = d_B_array_rotation[tx + (working_element) * (NGLL3) - (0)];\n\
  dpotentialdx_with_rot[0 - (0)] = dpotentialdxl + (A_rotation) * (cos_two_omega_t) + (B_rotation) * (sin_two_omega_t);\n\
  dpotentialdy_with_rot[0 - (0)] = dpotentialdyl + ( -(A_rotation)) * (sin_two_omega_t) + (B_rotation) * (cos_two_omega_t);\n\
  d_A_array_rotation[tx + (working_element) * (NGLL3) - (0)] = d_A_array_rotation[tx + (working_element) * (NGLL3) - (0)] + source_euler_A;\n\
  d_B_array_rotation[tx + (working_element) * (NGLL3) - (0)] = d_B_array_rotation[tx + (working_element) * (NGLL3) - (0)] + source_euler_B;\n\
}\n\
__kernel void outer_core_impl_kernel_adjoint(const int nb_blocks_to_compute, const __global int * d_ibool, const __global int * d_phase_ispec_inner, const int num_phase_ispec, const int d_iphase, const int use_mesh_coloring_gpu, const __global float * restrict d_potential, __global float * d_potential_dot_dot, const __global float * restrict d_xix, const __global float * restrict d_xiy, const __global float * restrict d_xiz, const __global float * restrict d_etax, const __global float * restrict d_etay, const __global float * restrict d_etaz, const __global float * restrict d_gammax, const __global float * restrict d_gammay, const __global float * restrict d_gammaz, const __global float * restrict d_hprime_xx, const __global float * restrict d_hprimewgll_xx, const __global float * restrict wgllwgll_xy, const __global float * restrict wgllwgll_xz, const __global float * restrict wgllwgll_yz, const int GRAVITY, const __global float * restrict d_xstore, const __global float * restrict d_ystore, const __global float * restrict d_zstore, const __global float * restrict d_d_ln_density_dr_table, const __global float * restrict d_minus_rho_g_over_kappa_fluid, const __global float * restrict wgll_cube, const int ROTATION, const float time, const float two_omega_earth, const float deltat, __global float * d_A_array_rotation, __global float * d_B_array_rotation, const int NSPEC_OUTER_CORE, __read_only image2d_t d_b_displ_oc_tex, __read_only image2d_t d_b_accel_oc_tex){\n\
#ifdef USE_TEXTURES_FIELDS\n\
  const sampler_t sampler_d_b_displ_oc_tex = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n\
  const sampler_t sampler_d_b_accel_oc_tex = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n\
#endif\n\
#ifdef USE_TEXTURES_CONSTANTS\n\
  const sampler_t sampler_d_hprime_xx_oc_tex = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n\
  const sampler_t sampler_d_hprimewgll_xx_oc_tex = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n\
#endif\n\
  int bx;\n\
  int tx;\n\
  int K;\n\
  int J;\n\
  int I;\n\
  int l;\n\
  ushort active;\n\
  int offset;\n\
  int iglob;\n\
  int working_element;\n\
  float temp1l;\n\
  float temp2l;\n\
  float temp3l;\n\
  float xixl;\n\
  float xiyl;\n\
  float xizl;\n\
  float etaxl;\n\
  float etayl;\n\
  float etazl;\n\
  float gammaxl;\n\
  float gammayl;\n\
  float gammazl;\n\
  float jacobianl;\n\
  float dpotentialdxl;\n\
  float dpotentialdyl;\n\
  float dpotentialdzl;\n\
  float dpotentialdx_with_rot;\n\
  float dpotentialdy_with_rot;\n\
  float sum_terms;\n\
  float gravity_term;\n\
  float gxl;\n\
  float gyl;\n\
  float gzl;\n\
  float radius;\n\
  float theta;\n\
  float phi;\n\
  float cos_theta;\n\
  float sin_theta;\n\
  float cos_phi;\n\
  float sin_phi;\n\
  float grad_x_ln_rho;\n\
  float grad_y_ln_rho;\n\
  float grad_z_ln_rho;\n\
  int int_radius;\n\
  __local float s_dummy_loc[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_temp1[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_temp2[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_temp3[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float sh_hprime_xx[NGLL2 + 0 - (1) - (0) + 1];\n\
  __local float sh_hprimewgll_xx[NGLL2 + 0 - (1) - (0) + 1];\n\
  bx = (get_group_id(1)) * (get_num_groups(0)) + get_group_id(0);\n\
  tx = get_local_id(0);\n\
  K = (tx) / (NGLL2);\n\
  J = (tx - ((K) * (NGLL2))) / (NGLLX);\n\
  I = tx - ((K) * (NGLL2)) - ((J) * (NGLLX));\n\
  active = (tx < NGLL3 && bx < nb_blocks_to_compute ? 1 : 0);\n\
  if(active){\n\
#ifdef USE_MESH_COLORING_GPU\n\
    working_element = bx;\n\
#else\n\
    if(use_mesh_coloring_gpu){\n\
      working_element = bx;\n\
    } else {\n\
      working_element = d_phase_ispec_inner[bx + (num_phase_ispec) * (d_iphase - (1)) - (0)] - (1);\n\
    }\n\
#endif\n\
    iglob = d_ibool[(working_element) * (NGLL3) + tx - (0)] - (1);\n\
#ifdef USE_TEXTURES_FIELDS\n\
    s_dummy_loc[tx - (0)] = as_float(read_imageui(d_b_displ_oc_tex, sampler_d_b_displ_oc_tex, int2(iglob,0)).x);\n\
#else\n\
    s_dummy_loc[tx - (0)] = d_potential[iglob - (0)];\n\
#endif\n\
  }\n\
  if(tx < NGLL2){\n\
#ifdef USE_TEXTURES_CONSTANTS\n\
    sh_hprime_xx[tx - (0)] = as_float(read_imageui(d_hprime_xx_oc_tex, sampler_d_hprime_xx_oc_tex, int2(tx,0)).x);\n\
    sh_hprimewgll_xx[tx - (0)] = as_float(read_imageui(d_hprimewgll_xx_oc_tex, sampler_d_hprimewgll_xx_oc_tex, int2(tx,0)).x);\n\
#else\n\
    sh_hprime_xx[tx - (0)] = d_hprime_xx[tx - (0)];\n\
    sh_hprimewgll_xx[tx - (0)] = d_hprimewgll_xx[tx - (0)];\n\
#endif\n\
  }\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if(active){\n\
    temp1l = 0.0f;\n\
    temp2l = 0.0f;\n\
    temp3l = 0.0f;\n\
#ifdef MANUALLY_UNROLLED_LOOPS\n\
    temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 0 - (0)]) * (sh_hprime_xx[(0) * (NGLLX) + I - (0)]);\n\
    temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (0) * (NGLLX) + I - (0)]) * (sh_hprime_xx[(0) * (NGLLX) + J - (0)]);\n\
    temp3l = temp3l + (s_dummy_loc[(0) * (NGLL2) + (J) * (NGLLX) + I - (0)]) * (sh_hprime_xx[(0) * (NGLLX) + K - (0)]);\n\
    temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 1 - (0)]) * (sh_hprime_xx[(1) * (NGLLX) + I - (0)]);\n\
    temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (1) * (NGLLX) + I - (0)]) * (sh_hprime_xx[(1) * (NGLLX) + J - (0)]);\n\
    temp3l = temp3l + (s_dummy_loc[(1) * (NGLL2) + (J) * (NGLLX) + I - (0)]) * (sh_hprime_xx[(1) * (NGLLX) + K - (0)]);\n\
    temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 2 - (0)]) * (sh_hprime_xx[(2) * (NGLLX) + I - (0)]);\n\
    temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (2) * (NGLLX) + I - (0)]) * (sh_hprime_xx[(2) * (NGLLX) + J - (0)]);\n\
    temp3l = temp3l + (s_dummy_loc[(2) * (NGLL2) + (J) * (NGLLX) + I - (0)]) * (sh_hprime_xx[(2) * (NGLLX) + K - (0)]);\n\
    temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 3 - (0)]) * (sh_hprime_xx[(3) * (NGLLX) + I - (0)]);\n\
    temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (3) * (NGLLX) + I - (0)]) * (sh_hprime_xx[(3) * (NGLLX) + J - (0)]);\n\
    temp3l = temp3l + (s_dummy_loc[(3) * (NGLL2) + (J) * (NGLLX) + I - (0)]) * (sh_hprime_xx[(3) * (NGLLX) + K - (0)]);\n\
    temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 4 - (0)]) * (sh_hprime_xx[(4) * (NGLLX) + I - (0)]);\n\
    temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (4) * (NGLLX) + I - (0)]) * (sh_hprime_xx[(4) * (NGLLX) + J - (0)]);\n\
    temp3l = temp3l + (s_dummy_loc[(4) * (NGLL2) + (J) * (NGLLX) + I - (0)]) * (sh_hprime_xx[(4) * (NGLLX) + K - (0)]);\n\
#else\n\
    for(l=0; l<=NGLLX - (1); l+=1){\n\
      temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + l - (0)]) * (sh_hprime_xx[(l) * (NGLLX) + I - (0)]);\n\
      temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (l) * (NGLLX) + I - (0)]) * (sh_hprime_xx[(l) * (NGLLX) + J - (0)]);\n\
      temp3l = temp3l + (s_dummy_loc[(l) * (NGLL2) + (J) * (NGLLX) + I - (0)]) * (sh_hprime_xx[(l) * (NGLLX) + K - (0)]);\n\
    }\n\
#endif\n\
    offset = (working_element) * (NGLL3_PADDED) + tx;\n\
    xixl = d_xix[offset - (0)];\n\
    etaxl = d_etax[offset - (0)];\n\
    gammaxl = d_gammax[offset - (0)];\n\
    xiyl = d_xiy[offset - (0)];\n\
    etayl = d_etay[offset - (0)];\n\
    gammayl = d_gammay[offset - (0)];\n\
    xizl = d_xiz[offset - (0)];\n\
    etazl = d_etaz[offset - (0)];\n\
    gammazl = d_gammaz[offset - (0)];\n\
    jacobianl = (1.0f) / ((xixl) * ((etayl) * (gammazl) - ((etazl) * (gammayl))) - ((xiyl) * ((etaxl) * (gammazl) - ((etazl) * (gammaxl)))) + (xizl) * ((etaxl) * (gammayl) - ((etayl) * (gammaxl))));\n\
    dpotentialdxl = (xixl) * (temp1l) + (etaxl) * (temp2l) + (gammaxl) * (temp3l);\n\
    dpotentialdyl = (xiyl) * (temp1l) + (etayl) * (temp2l) + (gammayl) * (temp3l);\n\
    dpotentialdzl = (xizl) * (temp1l) + (etazl) * (temp2l) + (gammazl) * (temp3l);\n\
    if(ROTATION){\n\
      compute_element_oc_rotation(tx, working_element, time, two_omega_earth, deltat, d_A_array_rotation, d_B_array_rotation, dpotentialdxl, dpotentialdyl,  &dpotentialdx_with_rot,  &dpotentialdy_with_rot);\n\
    } else {\n\
      dpotentialdx_with_rot = dpotentialdxl;\n\
      dpotentialdy_with_rot = dpotentialdyl;\n\
    }\n\
    radius = d_xstore[iglob - (0)];\n\
    theta = d_ystore[iglob - (0)];\n\
    phi = d_zstore[iglob - (0)];\n\
    sin_theta = sincos(theta,  &cos_theta);\n\
    sin_phi = sincos(phi,  &cos_phi);\n\
    int_radius = rint(((radius) * (R_EARTH_KM)) * (10.0f)) - (1);\n\
    if( ! GRAVITY){\n\
      grad_x_ln_rho = ((sin_theta) * (cos_phi)) * (d_d_ln_density_dr_table[int_radius - (0)]);\n\
      grad_y_ln_rho = ((sin_theta) * (sin_phi)) * (d_d_ln_density_dr_table[int_radius - (0)]);\n\
      grad_z_ln_rho = (cos_theta) * (d_d_ln_density_dr_table[int_radius - (0)]);\n\
      dpotentialdx_with_rot = dpotentialdx_with_rot + (s_dummy_loc[tx - (0)]) * (grad_x_ln_rho);\n\
      dpotentialdy_with_rot = dpotentialdy_with_rot + (s_dummy_loc[tx - (0)]) * (grad_y_ln_rho);\n\
      dpotentialdzl = dpotentialdzl + (s_dummy_loc[tx - (0)]) * (grad_z_ln_rho);\n\
    } else {\n\
      gxl = (sin_theta) * (cos_phi);\n\
      gyl = (sin_theta) * (sin_phi);\n\
      gzl = cos_theta;\n\
      gravity_term = (((d_minus_rho_g_over_kappa_fluid[int_radius - (0)]) * (jacobianl)) * (wgll_cube[tx - (0)])) * ((dpotentialdx_with_rot) * (gxl) + (dpotentialdy_with_rot) * (gyl) + (dpotentialdzl) * (gzl));\n\
    }\n\
    s_temp1[tx - (0)] = (jacobianl) * ((xixl) * (dpotentialdx_with_rot) + (xiyl) * (dpotentialdy_with_rot) + (xizl) * (dpotentialdzl));\n\
    s_temp2[tx - (0)] = (jacobianl) * ((etaxl) * (dpotentialdx_with_rot) + (etayl) * (dpotentialdy_with_rot) + (etazl) * (dpotentialdzl));\n\
    s_temp3[tx - (0)] = (jacobianl) * ((gammaxl) * (dpotentialdx_with_rot) + (gammayl) * (dpotentialdy_with_rot) + (gammazl) * (dpotentialdzl));\n\
  }\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if(active){\n\
    temp1l = 0.0f;\n\
    temp2l = 0.0f;\n\
    temp3l = 0.0f;\n\
#ifdef MANUALLY_UNROLLED_LOOPS\n\
    temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + 0 - (0)]) * (sh_hprimewgll_xx[(I) * (NGLLX) + 0 - (0)]);\n\
    temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (0) * (NGLLX) + I - (0)]) * (sh_hprimewgll_xx[(J) * (NGLLX) + 0 - (0)]);\n\
    temp3l = temp3l + (s_temp3[(0) * (NGLL2) + (J) * (NGLLX) + I - (0)]) * (sh_hprimewgll_xx[(K) * (NGLLX) + 0 - (0)]);\n\
    temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + 1 - (0)]) * (sh_hprimewgll_xx[(I) * (NGLLX) + 1 - (0)]);\n\
    temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (1) * (NGLLX) + I - (0)]) * (sh_hprimewgll_xx[(J) * (NGLLX) + 1 - (0)]);\n\
    temp3l = temp3l + (s_temp3[(1) * (NGLL2) + (J) * (NGLLX) + I - (0)]) * (sh_hprimewgll_xx[(K) * (NGLLX) + 1 - (0)]);\n\
    temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + 2 - (0)]) * (sh_hprimewgll_xx[(I) * (NGLLX) + 2 - (0)]);\n\
    temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (2) * (NGLLX) + I - (0)]) * (sh_hprimewgll_xx[(J) * (NGLLX) + 2 - (0)]);\n\
    temp3l = temp3l + (s_temp3[(2) * (NGLL2) + (J) * (NGLLX) + I - (0)]) * (sh_hprimewgll_xx[(K) * (NGLLX) + 2 - (0)]);\n\
    temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + 3 - (0)]) * (sh_hprimewgll_xx[(I) * (NGLLX) + 3 - (0)]);\n\
    temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (3) * (NGLLX) + I - (0)]) * (sh_hprimewgll_xx[(J) * (NGLLX) + 3 - (0)]);\n\
    temp3l = temp3l + (s_temp3[(3) * (NGLL2) + (J) * (NGLLX) + I - (0)]) * (sh_hprimewgll_xx[(K) * (NGLLX) + 3 - (0)]);\n\
    temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + 4 - (0)]) * (sh_hprimewgll_xx[(I) * (NGLLX) + 4 - (0)]);\n\
    temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (4) * (NGLLX) + I - (0)]) * (sh_hprimewgll_xx[(J) * (NGLLX) + 4 - (0)]);\n\
    temp3l = temp3l + (s_temp3[(4) * (NGLL2) + (J) * (NGLLX) + I - (0)]) * (sh_hprimewgll_xx[(K) * (NGLLX) + 4 - (0)]);\n\
#else\n\
    for(l=0; l<=NGLLX - (1); l+=1){\n\
      temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + l - (0)]) * (sh_hprimewgll_xx[(I) * (NGLLX) + l - (0)]);\n\
      temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (l) * (NGLLX) + I - (0)]) * (sh_hprimewgll_xx[(J) * (NGLLX) + l - (0)]);\n\
      temp3l = temp3l + (s_temp3[(l) * (NGLL2) + (J) * (NGLLX) + I - (0)]) * (sh_hprimewgll_xx[(K) * (NGLLX) + l - (0)]);\n\
    }\n\
#endif\n\
    sum_terms =  -((wgllwgll_yz[(K) * (NGLLX) + J - (0)]) * (temp1l) + (wgllwgll_xz[(K) * (NGLLX) + I - (0)]) * (temp2l) + (wgllwgll_xy[(J) * (NGLLX) + I - (0)]) * (temp3l));\n\
    if(GRAVITY){\n\
      sum_terms = sum_terms + gravity_term;\n\
    }\n\
#ifdef USE_MESH_COLORING_GPU\n\
#ifdef USE_TEXTURES_FIELDS\n\
    d_potential_dot_dot[iglob - (0)] = as_float(read_imageui(d_b_accel_oc_tex, sampler_d_b_accel_oc_tex, int2(iglob,0)).x) + sum_terms;\n\
#else\n\
    d_potential_dot_dot[iglob - (0)] = d_potential_dot_dot[iglob - (0)] + sum_terms;\n\
#endif\n\
#else\n\
    if(use_mesh_coloring_gpu){\n\
      if(NSPEC_OUTER_CORE > 1000){\n\
#ifdef USE_TEXTURES_FIELDS\n\
        d_potential_dot_dot[iglob - (0)] = as_float(read_imageui(d_b_accel_oc_tex, sampler_d_b_accel_oc_tex, int2(iglob,0)).x) + sum_terms;\n\
#else\n\
        d_potential_dot_dot[iglob - (0)] = d_potential_dot_dot[iglob - (0)] + sum_terms;\n\
#endif\n\
      } else {\n\
        atomicAdd(d_potential_dot_dot + iglob, sum_terms);\n\
      }\n\
    } else {\n\
      atomicAdd(d_potential_dot_dot + iglob, sum_terms);\n\
    }\n\
#endif\n\
  }\n\
}\n\
";
