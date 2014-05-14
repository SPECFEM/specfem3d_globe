const char * inner_core_impl_kernel_adjoint_program = "\
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
void compute_element_ic_att_stress(const int tx, const int working_element, const __global float * R_xx, const __global float * R_yy, const __global float * R_xy, const __global float * R_xz, const __global float * R_yz, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_xz, float * sigma_yz){\n\
  int offset;\n\
  int i_sls;\n\
  float R_xx_val;\n\
  float R_yy_val;\n\
  for(i_sls=0; i_sls<=N_SLS - (1); i_sls+=1){\n\
    offset = i_sls + (N_SLS) * (tx + (NGLL3) * (working_element));\n\
    R_xx_val = R_xx[offset - 0];\n\
    R_yy_val = R_yy[offset - 0];\n\
    sigma_xx[0 - 0] = sigma_xx[0 - 0] - (R_xx_val);\n\
    sigma_yy[0 - 0] = sigma_yy[0 - 0] - (R_yy_val);\n\
    sigma_zz[0 - 0] = sigma_zz[0 - 0] + R_xx_val + R_yy_val;\n\
    sigma_xy[0 - 0] = sigma_xy[0 - 0] - (R_xy[offset - 0]);\n\
    sigma_xz[0 - 0] = sigma_xz[0 - 0] - (R_xz[offset - 0]);\n\
    sigma_yz[0 - 0] = sigma_yz[0 - 0] - (R_yz[offset - 0]);\n\
  }\n\
}\n\
void compute_element_ic_gravity(const int tx, const int iglob, const __global float * restrict d_xstore, const __global float * restrict d_ystore, const __global float * restrict d_zstore, const __global float * restrict d_minus_gravity_table, const __global float * restrict d_minus_deriv_gravity_table, const __global float * restrict d_density_table, const __global float * restrict wgll_cube, const float jacobianl, const __local float * s_dummyx_loc, const __local float * s_dummyy_loc, const __local float * s_dummyz_loc, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_yx, float * sigma_xz, float * sigma_zx, float * sigma_yz, float * sigma_zy, float * rho_s_H1, float * rho_s_H2, float * rho_s_H3){\n\
  float radius;\n\
  float theta;\n\
  float phi;\n\
  float cos_theta;\n\
  float sin_theta;\n\
  float cos_phi;\n\
  float sin_phi;\n\
  float cos_theta_sq;\n\
  float sin_theta_sq;\n\
  float cos_phi_sq;\n\
  float sin_phi_sq;\n\
  float minus_g;\n\
  float minus_dg;\n\
  float rho;\n\
  float gxl;\n\
  float gyl;\n\
  float gzl;\n\
  float minus_g_over_radius;\n\
  float minus_dg_plus_g_over_radius;\n\
  float Hxxl;\n\
  float Hyyl;\n\
  float Hzzl;\n\
  float Hxyl;\n\
  float Hxzl;\n\
  float Hyzl;\n\
  float sx_l;\n\
  float sy_l;\n\
  float sz_l;\n\
  float factor;\n\
  int int_radius;\n\
  radius = d_xstore[iglob - 0];\n\
  if(radius < 1.5696123057604773e-05f){\n\
    radius = 1.5696123057604773e-05f;\n\
  }\n\
  theta = d_ystore[iglob - 0];\n\
  phi = d_zstore[iglob - 0];\n\
  sin_theta = sincos(theta,  &cos_theta);\n\
  sin_phi = sincos(phi,  &cos_phi);\n\
  int_radius = rint(((radius) * (6371.0f)) * (10.0f)) - (1);\n\
  if(int_radius < 0){\n\
    int_radius = 0;\n\
  }\n\
  minus_g = d_minus_gravity_table[int_radius - 0];\n\
  minus_dg = d_minus_deriv_gravity_table[int_radius - 0];\n\
  rho = d_density_table[int_radius - 0];\n\
  gxl = ((minus_g) * (sin_theta)) * (cos_phi);\n\
  gyl = ((minus_g) * (sin_theta)) * (sin_phi);\n\
  gzl = (minus_g) * (cos_theta);\n\
  minus_g_over_radius = (minus_g) / (radius);\n\
  minus_dg_plus_g_over_radius = minus_dg - (minus_g_over_radius);\n\
  cos_theta_sq = (cos_theta) * (cos_theta);\n\
  sin_theta_sq = (sin_theta) * (sin_theta);\n\
  cos_phi_sq = (cos_phi) * (cos_phi);\n\
  sin_phi_sq = (sin_phi) * (sin_phi);\n\
  Hxxl = (minus_g_over_radius) * ((cos_phi_sq) * (cos_theta_sq) + sin_phi_sq) + ((cos_phi_sq) * (minus_dg)) * (sin_theta_sq);\n\
  Hyyl = (minus_g_over_radius) * (cos_phi_sq + (cos_theta_sq) * (sin_phi_sq)) + ((minus_dg) * (sin_phi_sq)) * (sin_theta_sq);\n\
  Hzzl = (cos_theta_sq) * (minus_dg) + (minus_g_over_radius) * (sin_theta_sq);\n\
  Hxyl = (((cos_phi) * (minus_dg_plus_g_over_radius)) * (sin_phi)) * (sin_theta_sq);\n\
  Hxzl = (((cos_phi) * (cos_theta)) * (minus_dg_plus_g_over_radius)) * (sin_theta);\n\
  Hyzl = (((cos_theta) * (minus_dg_plus_g_over_radius)) * (sin_phi)) * (sin_theta);\n\
  sx_l = (rho) * (s_dummyx_loc[tx - 0]);\n\
  sy_l = (rho) * (s_dummyy_loc[tx - 0]);\n\
  sz_l = (rho) * (s_dummyz_loc[tx - 0]);\n\
  *(sigma_xx) = *(sigma_xx) + (sy_l) * (gyl) + (sz_l) * (gzl);\n\
  *(sigma_yy) = *(sigma_yy) + (sx_l) * (gxl) + (sz_l) * (gzl);\n\
  *(sigma_zz) = *(sigma_zz) + (sx_l) * (gxl) + (sy_l) * (gyl);\n\
  *(sigma_xy) = *(sigma_xy) - ((sx_l) * (gyl));\n\
  *(sigma_yx) = *(sigma_yx) - ((sy_l) * (gxl));\n\
  *(sigma_xz) = *(sigma_xz) - ((sx_l) * (gzl));\n\
  *(sigma_zx) = *(sigma_zx) - ((sz_l) * (gxl));\n\
  *(sigma_yz) = *(sigma_yz) - ((sy_l) * (gzl));\n\
  *(sigma_zy) = *(sigma_zy) - ((sz_l) * (gyl));\n\
  factor = (jacobianl) * (wgll_cube[tx - 0]);\n\
  rho_s_H1[0 - 0] = (factor) * ((sx_l) * (Hxxl) + (sy_l) * (Hxyl) + (sz_l) * (Hxzl));\n\
  rho_s_H2[0 - 0] = (factor) * ((sx_l) * (Hxyl) + (sy_l) * (Hyyl) + (sz_l) * (Hyzl));\n\
  rho_s_H3[0 - 0] = (factor) * ((sx_l) * (Hxzl) + (sy_l) * (Hyzl) + (sz_l) * (Hzzl));\n\
}\n\
void compute_element_ic_att_memory(const int tx, const int working_element, const __global float * d_muv, const __global float * factor_common, const __global float * alphaval, const __global float * betaval, const __global float * gammaval, __global float * R_xx, __global float * R_yy, __global float * R_xy, __global float * R_xz, __global float * R_yz, const __global float * epsilondev_xx, const __global float * epsilondev_yy, const __global float * epsilondev_xy, const __global float * epsilondev_xz, const __global float * epsilondev_yz, const float epsilondev_xx_loc, const float epsilondev_yy_loc, const float epsilondev_xy_loc, const float epsilondev_xz_loc, const float epsilondev_yz_loc, const int USE_3D_ATTENUATION_ARRAYS){\n\
  int offset;\n\
  int i_sls;\n\
  float mul;\n\
  float alphaval_loc;\n\
  float betaval_loc;\n\
  float gammaval_loc;\n\
  float factor_loc;\n\
  float sn;\n\
  float snp1;\n\
  mul = d_muv[tx + (NGLL3_PADDED) * (working_element) - 0];\n\
  for(i_sls=0; i_sls<=N_SLS - (1); i_sls+=1){\n\
    offset = i_sls + (N_SLS) * (tx + (NGLL3) * (working_element));\n\
    if(USE_3D_ATTENUATION_ARRAYS){\n\
      factor_loc = (mul) * (factor_common[offset - 0]);\n\
    } else {\n\
      factor_loc = (mul) * (factor_common[i_sls + (N_SLS) * (working_element) - 0]);\n\
    }\n\
    alphaval_loc = alphaval[i_sls - 0];\n\
    betaval_loc = betaval[i_sls - 0];\n\
    gammaval_loc = gammaval[i_sls - 0];\n\
    sn = (factor_loc) * (epsilondev_xx[tx + (NGLL3) * (working_element) - 0]);\n\
    snp1 = (factor_loc) * (epsilondev_xx_loc);\n\
    R_xx[offset - 0] = (alphaval_loc) * (R_xx[offset - 0]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);\n\
    sn = (factor_loc) * (epsilondev_yy[tx + (NGLL3) * (working_element) - 0]);\n\
    snp1 = (factor_loc) * (epsilondev_yy_loc);\n\
    R_yy[offset - 0] = (alphaval_loc) * (R_yy[offset - 0]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);\n\
    sn = (factor_loc) * (epsilondev_xy[tx + (NGLL3) * (working_element) - 0]);\n\
    snp1 = (factor_loc) * (epsilondev_xy_loc);\n\
    R_xy[offset - 0] = (alphaval_loc) * (R_xy[offset - 0]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);\n\
    sn = (factor_loc) * (epsilondev_xz[tx + (NGLL3) * (working_element) - 0]);\n\
    snp1 = (factor_loc) * (epsilondev_xz_loc);\n\
    R_xz[offset - 0] = (alphaval_loc) * (R_xz[offset - 0]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);\n\
    sn = (factor_loc) * (epsilondev_yz[tx + (NGLL3) * (working_element) - 0]);\n\
    snp1 = (factor_loc) * (epsilondev_yz_loc);\n\
    R_yz[offset - 0] = (alphaval_loc) * (R_yz[offset - 0]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);\n\
  }\n\
}\n\
__kernel void inner_core_impl_kernel_adjoint(const int nb_blocks_to_compute, const __global int * d_ibool, const __global int * d_idoubling, const __global int * d_phase_ispec_inner, const int num_phase_ispec, const int d_iphase, const float deltat, const int use_mesh_coloring_gpu, const __global float * restrict d_displ, __global float * d_accel, const __global float * restrict d_xix, const __global float * restrict d_xiy, const __global float * restrict d_xiz, const __global float * restrict d_etax, const __global float * restrict d_etay, const __global float * restrict d_etaz, const __global float * restrict d_gammax, const __global float * restrict d_gammay, const __global float * restrict d_gammaz, const __global float * restrict d_hprime_xx, const __global float * restrict d_hprimewgll_xx, const __global float * restrict d_wgllwgll_xy, const __global float * restrict d_wgllwgll_xz, const __global float * restrict d_wgllwgll_yz, const __global float * restrict d_kappavstore, const __global float * restrict d_muvstore, const int COMPUTE_AND_STORE_STRAIN, __global float * epsilondev_xx, __global float * epsilondev_yy, __global float * epsilondev_xy, __global float * epsilondev_xz, __global float * epsilondev_yz, __global float * epsilon_trace_over_3, const int ATTENUATION, const int PARTIAL_PHYS_DISPERSION_ONLY, const int USE_3D_ATTENUATION_ARRAYS, const __global float * restrict one_minus_sum_beta, const __global float * restrict factor_common, __global float * R_xx, __global float * R_yy, __global float * R_xy, __global float * R_xz, __global float * R_yz, const __global float * restrict alphaval, const __global float * restrict betaval, const __global float * restrict gammaval, const int ANISOTROPY, const __global float * restrict d_c11store, const __global float * restrict d_c12store, const __global float * restrict d_c13store, const __global float * restrict d_c33store, const __global float * restrict d_c44store, const int GRAVITY, const __global float * restrict d_xstore, const __global float * restrict d_ystore, const __global float * restrict d_zstore, const __global float * restrict d_minus_gravity_table, const __global float * restrict d_minus_deriv_gravity_table, const __global float * restrict d_density_table, const __global float * restrict wgll_cube, const int NSPEC_INNER_CORE_STRAIN_ONLY, const int NSPEC_INNER_CORE, __read_only image2d_t d_b_displ_ic_tex, __read_only image2d_t d_b_accel_ic_tex){\n\
#ifdef USE_TEXTURES_FIELDS\n\
  const sampler_t sampler_d_b_displ_ic_tex = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n\
  const sampler_t sampler_d_b_accel_ic_tex = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n\
#endif\n\
#ifdef USE_TEXTURES_CONSTANTS\n\
  const sampler_t sampler_d_hprime_xx_ic_tex = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;\n\
#endif\n\
  int bx;\n\
  int tx;\n\
  int K;\n\
  int J;\n\
  int I;\n\
  ushort active;\n\
  int offset;\n\
  int iglob;\n\
  int working_element;\n\
  int l;\n\
  float tempx1l;\n\
  float tempx2l;\n\
  float tempx3l;\n\
  float tempy1l;\n\
  float tempy2l;\n\
  float tempy3l;\n\
  float tempz1l;\n\
  float tempz2l;\n\
  float tempz3l;\n\
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
  float duxdxl;\n\
  float duxdyl;\n\
  float duxdzl;\n\
  float duydxl;\n\
  float duydyl;\n\
  float duydzl;\n\
  float duzdxl;\n\
  float duzdyl;\n\
  float duzdzl;\n\
  float duxdxl_plus_duydyl;\n\
  float duxdxl_plus_duzdzl;\n\
  float duydyl_plus_duzdzl;\n\
  float duxdyl_plus_duydxl;\n\
  float duzdxl_plus_duxdzl;\n\
  float duzdyl_plus_duydzl;\n\
  float templ;\n\
  float fac1;\n\
  float fac2;\n\
  float fac3;\n\
  float lambdal;\n\
  float mul;\n\
  float lambdalplus2mul;\n\
  float kappal;\n\
  float mul_iso;\n\
  float mul_aniso;\n\
  float sigma_xx;\n\
  float sigma_xy;\n\
  float sigma_xz;\n\
  float sigma_yx;\n\
  float sigma_yy;\n\
  float sigma_yz;\n\
  float sigma_zx;\n\
  float sigma_zy;\n\
  float sigma_zz;\n\
  float epsilondev_xx_loc;\n\
  float epsilondev_yy_loc;\n\
  float epsilondev_xy_loc;\n\
  float epsilondev_xz_loc;\n\
  float epsilondev_yz_loc;\n\
  float c11;\n\
  float c12;\n\
  float c13;\n\
  float c33;\n\
  float c44;\n\
  float sum_terms1;\n\
  float sum_terms2;\n\
  float sum_terms3;\n\
  float rho_s_H1;\n\
  float rho_s_H2;\n\
  float rho_s_H3;\n\
  __local float s_dummyx_loc[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_dummyy_loc[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_dummyz_loc[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_tempx1[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_tempx2[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_tempx3[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_tempy1[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_tempy2[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_tempy3[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_tempz1[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_tempz2[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_tempz3[NGLL3 + 0 - (1) - (0) + 1];\n\
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
      working_element = d_phase_ispec_inner[bx + (num_phase_ispec) * (d_iphase - (1)) - 0] - (1);\n\
    }\n\
#endif\n\
    if(d_idoubling[working_element - 0] == IFLAG_IN_FICTITIOUS_CUBE){\n\
      active = 0;\n\
    } else {\n\
      iglob = d_ibool[(working_element) * (NGLL3) + tx - 0] - (1);\n\
#ifdef USE_TEXTURES_FIELDS\n\
      s_dummyx_loc[tx - 0] = as_float(read_imageui(d_b_displ_ic_tex, sampler_d_b_displ_ic_tex, int2((iglob) * (3) + 0,0)).x);\n\
      s_dummyy_loc[tx - 0] = as_float(read_imageui(d_b_displ_ic_tex, sampler_d_b_displ_ic_tex, int2((iglob) * (3) + 1,0)).x);\n\
      s_dummyz_loc[tx - 0] = as_float(read_imageui(d_b_displ_ic_tex, sampler_d_b_displ_ic_tex, int2((iglob) * (3) + 2,0)).x);\n\
#else\n\
      s_dummyx_loc[tx - 0] = d_displ[0 - 0 + (iglob - (0)) * (3)];\n\
      s_dummyy_loc[tx - 0] = d_displ[1 - 0 + (iglob - (0)) * (3)];\n\
      s_dummyz_loc[tx - 0] = d_displ[2 - 0 + (iglob - (0)) * (3)];\n\
#endif\n\
    }\n\
  }\n\
  if(tx < NGLL2){\n\
#ifdef USE_TEXTURES_CONSTANTS\n\
    sh_hprime_xx[tx - 0] = as_float(read_imageui(d_hprime_xx_ic_tex, sampler_d_hprime_xx_ic_tex, int2(tx,0)).x);\n\
    sh_hprimewgll_xx[tx - 0] = as_float(read_imageui(d_hprimewgll_xx_ic_tex, sampler_d_hprimewgll_xx_ic_tex, int2(tx,0)).x);\n\
#else\n\
    sh_hprime_xx[tx - 0] = d_hprime_xx[tx - 0];\n\
    sh_hprimewgll_xx[tx - 0] = d_hprimewgll_xx[tx - 0];\n\
#endif\n\
  }\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if(active){\n\
    tempx1l = 0.0f;\n\
    tempx2l = 0.0f;\n\
    tempx3l = 0.0f;\n\
    tempy1l = 0.0f;\n\
    tempy2l = 0.0f;\n\
    tempy3l = 0.0f;\n\
    tempz1l = 0.0f;\n\
    tempz2l = 0.0f;\n\
    tempz3l = 0.0f;\n\
#ifdef MANUALLY_UNROLLED_LOOPS\n\
    fac1 = sh_hprime_xx[(0) * (NGLLX) + I - 0];\n\
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 0 - 0]) * (fac1);\n\
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 0 - 0]) * (fac1);\n\
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 0 - 0]) * (fac1);\n\
    fac2 = sh_hprime_xx[(0) * (NGLLX) + J - 0];\n\
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (0) * (NGLLX) + I - 0]) * (fac2);\n\
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (0) * (NGLLX) + I - 0]) * (fac2);\n\
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (0) * (NGLLX) + I - 0]) * (fac2);\n\
    fac3 = sh_hprime_xx[(0) * (NGLLX) + K - 0];\n\
    tempx3l = tempx3l + (s_dummyx_loc[(0) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);\n\
    tempy3l = tempy3l + (s_dummyy_loc[(0) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);\n\
    tempz3l = tempz3l + (s_dummyz_loc[(0) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);\n\
    fac1 = sh_hprime_xx[(1) * (NGLLX) + I - 0];\n\
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 1 - 0]) * (fac1);\n\
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 1 - 0]) * (fac1);\n\
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 1 - 0]) * (fac1);\n\
    fac2 = sh_hprime_xx[(1) * (NGLLX) + J - 0];\n\
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (1) * (NGLLX) + I - 0]) * (fac2);\n\
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (1) * (NGLLX) + I - 0]) * (fac2);\n\
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (1) * (NGLLX) + I - 0]) * (fac2);\n\
    fac3 = sh_hprime_xx[(1) * (NGLLX) + K - 0];\n\
    tempx3l = tempx3l + (s_dummyx_loc[(1) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);\n\
    tempy3l = tempy3l + (s_dummyy_loc[(1) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);\n\
    tempz3l = tempz3l + (s_dummyz_loc[(1) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);\n\
    fac1 = sh_hprime_xx[(2) * (NGLLX) + I - 0];\n\
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 2 - 0]) * (fac1);\n\
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 2 - 0]) * (fac1);\n\
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 2 - 0]) * (fac1);\n\
    fac2 = sh_hprime_xx[(2) * (NGLLX) + J - 0];\n\
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (2) * (NGLLX) + I - 0]) * (fac2);\n\
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (2) * (NGLLX) + I - 0]) * (fac2);\n\
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (2) * (NGLLX) + I - 0]) * (fac2);\n\
    fac3 = sh_hprime_xx[(2) * (NGLLX) + K - 0];\n\
    tempx3l = tempx3l + (s_dummyx_loc[(2) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);\n\
    tempy3l = tempy3l + (s_dummyy_loc[(2) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);\n\
    tempz3l = tempz3l + (s_dummyz_loc[(2) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);\n\
    fac1 = sh_hprime_xx[(3) * (NGLLX) + I - 0];\n\
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 3 - 0]) * (fac1);\n\
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 3 - 0]) * (fac1);\n\
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 3 - 0]) * (fac1);\n\
    fac2 = sh_hprime_xx[(3) * (NGLLX) + J - 0];\n\
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (3) * (NGLLX) + I - 0]) * (fac2);\n\
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (3) * (NGLLX) + I - 0]) * (fac2);\n\
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (3) * (NGLLX) + I - 0]) * (fac2);\n\
    fac3 = sh_hprime_xx[(3) * (NGLLX) + K - 0];\n\
    tempx3l = tempx3l + (s_dummyx_loc[(3) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);\n\
    tempy3l = tempy3l + (s_dummyy_loc[(3) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);\n\
    tempz3l = tempz3l + (s_dummyz_loc[(3) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);\n\
    fac1 = sh_hprime_xx[(4) * (NGLLX) + I - 0];\n\
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 4 - 0]) * (fac1);\n\
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 4 - 0]) * (fac1);\n\
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 4 - 0]) * (fac1);\n\
    fac2 = sh_hprime_xx[(4) * (NGLLX) + J - 0];\n\
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (4) * (NGLLX) + I - 0]) * (fac2);\n\
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (4) * (NGLLX) + I - 0]) * (fac2);\n\
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (4) * (NGLLX) + I - 0]) * (fac2);\n\
    fac3 = sh_hprime_xx[(4) * (NGLLX) + K - 0];\n\
    tempx3l = tempx3l + (s_dummyx_loc[(4) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);\n\
    tempy3l = tempy3l + (s_dummyy_loc[(4) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);\n\
    tempz3l = tempz3l + (s_dummyz_loc[(4) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);\n\
#else\n\
    for(l=0; l<=NGLLX - (1); l+=1){\n\
      fac1 = sh_hprime_xx[(l) * (NGLLX) + I - 0];\n\
      tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + l - 0]) * (fac1);\n\
      tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + l - 0]) * (fac1);\n\
      tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + l - 0]) * (fac1);\n\
      fac2 = sh_hprime_xx[(l) * (NGLLX) + J - 0];\n\
      tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (l) * (NGLLX) + I - 0]) * (fac2);\n\
      tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (l) * (NGLLX) + I - 0]) * (fac2);\n\
      tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (l) * (NGLLX) + I - 0]) * (fac2);\n\
      fac3 = sh_hprime_xx[(l) * (NGLLX) + K - 0];\n\
      tempx3l = tempx3l + (s_dummyx_loc[(l) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);\n\
      tempy3l = tempy3l + (s_dummyy_loc[(l) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);\n\
      tempz3l = tempz3l + (s_dummyz_loc[(l) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);\n\
    }\n\
#endif\n\
    offset = (working_element) * (NGLL3_PADDED) + tx;\n\
    xixl = d_xix[offset - 0];\n\
    etaxl = d_etax[offset - 0];\n\
    gammaxl = d_gammax[offset - 0];\n\
    xiyl = d_xiy[offset - 0];\n\
    etayl = d_etay[offset - 0];\n\
    gammayl = d_gammay[offset - 0];\n\
    xizl = d_xiz[offset - 0];\n\
    etazl = d_etaz[offset - 0];\n\
    gammazl = d_gammaz[offset - 0];\n\
    duxdxl = (xixl) * (tempx1l) + (etaxl) * (tempx2l) + (gammaxl) * (tempx3l);\n\
    duxdyl = (xiyl) * (tempx1l) + (etayl) * (tempx2l) + (gammayl) * (tempx3l);\n\
    duxdzl = (xizl) * (tempx1l) + (etazl) * (tempx2l) + (gammazl) * (tempx3l);\n\
    duydxl = (xixl) * (tempy1l) + (etaxl) * (tempy2l) + (gammaxl) * (tempy3l);\n\
    duydyl = (xiyl) * (tempy1l) + (etayl) * (tempy2l) + (gammayl) * (tempy3l);\n\
    duydzl = (xizl) * (tempy1l) + (etazl) * (tempy2l) + (gammazl) * (tempy3l);\n\
    duzdxl = (xixl) * (tempz1l) + (etaxl) * (tempz2l) + (gammaxl) * (tempz3l);\n\
    duzdyl = (xiyl) * (tempz1l) + (etayl) * (tempz2l) + (gammayl) * (tempz3l);\n\
    duzdzl = (xizl) * (tempz1l) + (etazl) * (tempz2l) + (gammazl) * (tempz3l);\n\
    duxdxl_plus_duydyl = duxdxl + duydyl;\n\
    duxdxl_plus_duzdzl = duxdxl + duzdzl;\n\
    duydyl_plus_duzdzl = duydyl + duzdzl;\n\
    duxdyl_plus_duydxl = duxdyl + duydxl;\n\
    duzdxl_plus_duxdzl = duzdxl + duxdzl;\n\
    duzdyl_plus_duydzl = duzdyl + duydzl;\n\
    if(COMPUTE_AND_STORE_STRAIN){\n\
      templ = (duxdxl + duydyl + duzdzl) * (0.3333333333333333f);\n\
      epsilondev_xx_loc = duxdxl - (templ);\n\
      epsilondev_yy_loc = duydyl - (templ);\n\
      epsilondev_xy_loc = (duxdyl_plus_duydxl) * (0.5f);\n\
      epsilondev_xz_loc = (duzdxl_plus_duxdzl) * (0.5f);\n\
      epsilondev_yz_loc = (duzdyl_plus_duydzl) * (0.5f);\n\
      if(NSPEC_INNER_CORE_STRAIN_ONLY == 1){\n\
        epsilon_trace_over_3[tx - 0] = templ;\n\
      } else {\n\
        epsilon_trace_over_3[tx + (working_element) * (NGLL3) - 0] = templ;\n\
      }\n\
    }\n\
    kappal = d_kappavstore[offset - 0];\n\
    mul = d_muvstore[offset - 0];\n\
    if(ATTENUATION){\n\
      if(USE_3D_ATTENUATION_ARRAYS){\n\
        mul_iso = (mul) * (one_minus_sum_beta[tx + (working_element) * (NGLL3) - 0]);\n\
        mul_aniso = (mul) * (one_minus_sum_beta[tx + (working_element) * (NGLL3) - 0] - (1.0f));\n\
      } else {\n\
        mul_iso = (mul) * (one_minus_sum_beta[working_element - 0]);\n\
        mul_aniso = (mul) * (one_minus_sum_beta[working_element - 0] - (1.0f));\n\
      }\n\
    } else {\n\
      mul_iso = mul;\n\
    }\n\
    if(ANISOTROPY){\n\
      c11 = d_c11store[offset - 0];\n\
      c12 = d_c12store[offset - 0];\n\
      c13 = d_c13store[offset - 0];\n\
      c33 = d_c33store[offset - 0];\n\
      c44 = d_c44store[offset - 0];\n\
      if(ATTENUATION){\n\
        c11 = c11 + (mul_aniso) * (1.3333333333333333f);\n\
        c12 = c12 - ((mul_aniso) * (0.6666666666666666f));\n\
        c13 = c13 - ((mul_aniso) * (0.6666666666666666f));\n\
        c33 = c33 + (mul_aniso) * (1.3333333333333333f);\n\
        c44 = c44 + mul_aniso;\n\
      }\n\
      sigma_xx = (c11) * (duxdxl) + (c12) * (duydyl) + (c13) * (duzdzl);\n\
      sigma_yy = (c12) * (duxdxl) + (c11) * (duydyl) + (c13) * (duzdzl);\n\
      sigma_zz = (c13) * (duxdxl) + (c13) * (duydyl) + (c33) * (duzdzl);\n\
      sigma_xy = ((c11 - (c12)) * (duxdyl_plus_duydxl)) * (0.5f);\n\
      sigma_xz = (c44) * (duzdxl_plus_duxdzl);\n\
      sigma_yz = (c44) * (duzdyl_plus_duydzl);\n\
    } else {\n\
      lambdalplus2mul = kappal + (mul_iso) * (1.3333333333333333f);\n\
      lambdal = lambdalplus2mul - ((mul_iso) * (2.0f));\n\
      sigma_xx = (lambdalplus2mul) * (duxdxl) + (lambdal) * (duydyl_plus_duzdzl);\n\
      sigma_yy = (lambdalplus2mul) * (duydyl) + (lambdal) * (duxdxl_plus_duzdzl);\n\
      sigma_zz = (lambdalplus2mul) * (duzdzl) + (lambdal) * (duxdxl_plus_duydyl);\n\
      sigma_xy = (mul) * (duxdyl_plus_duydxl);\n\
      sigma_xz = (mul) * (duzdxl_plus_duxdzl);\n\
      sigma_yz = (mul) * (duzdyl_plus_duydzl);\n\
    }\n\
    if(ATTENUATION &&  ! PARTIAL_PHYS_DISPERSION_ONLY){\n\
      compute_element_ic_att_stress(tx, working_element, R_xx, R_yy, R_xy, R_xz, R_yz,  &sigma_xx,  &sigma_yy,  &sigma_zz,  &sigma_xy,  &sigma_xz,  &sigma_yz);\n\
    }\n\
    sigma_yx = sigma_xy;\n\
    sigma_zx = sigma_xz;\n\
    sigma_zy = sigma_yz;\n\
    jacobianl = (1.0f) / ((xixl) * ((etayl) * (gammazl) - ((etazl) * (gammayl))) - ((xiyl) * ((etaxl) * (gammazl) - ((etazl) * (gammaxl)))) + (xizl) * ((etaxl) * (gammayl) - ((etayl) * (gammaxl))));\n\
    if(GRAVITY){\n\
      compute_element_ic_gravity(tx, iglob, d_xstore, d_ystore, d_zstore, d_minus_gravity_table, d_minus_deriv_gravity_table, d_density_table, wgll_cube, jacobianl, s_dummyx_loc, s_dummyy_loc, s_dummyz_loc,  &sigma_xx,  &sigma_yy,  &sigma_zz,  &sigma_xy,  &sigma_yx,  &sigma_xz,  &sigma_zx,  &sigma_yz,  &sigma_zy,  &rho_s_H1,  &rho_s_H2,  &rho_s_H3);\n\
    }\n\
    s_tempx1[tx - 0] = (jacobianl) * ((sigma_xx) * (xixl) + (sigma_yx) * (xiyl) + (sigma_zx) * (xizl));\n\
    s_tempy1[tx - 0] = (jacobianl) * ((sigma_xy) * (xixl) + (sigma_yy) * (xiyl) + (sigma_zy) * (xizl));\n\
    s_tempz1[tx - 0] = (jacobianl) * ((sigma_xz) * (xixl) + (sigma_yz) * (xiyl) + (sigma_zz) * (xizl));\n\
    s_tempx2[tx - 0] = (jacobianl) * ((sigma_xx) * (etaxl) + (sigma_yx) * (etayl) + (sigma_zx) * (etazl));\n\
    s_tempy2[tx - 0] = (jacobianl) * ((sigma_xy) * (etaxl) + (sigma_yy) * (etayl) + (sigma_zy) * (etazl));\n\
    s_tempz2[tx - 0] = (jacobianl) * ((sigma_xz) * (etaxl) + (sigma_yz) * (etayl) + (sigma_zz) * (etazl));\n\
    s_tempx3[tx - 0] = (jacobianl) * ((sigma_xx) * (gammaxl) + (sigma_yx) * (gammayl) + (sigma_zx) * (gammazl));\n\
    s_tempy3[tx - 0] = (jacobianl) * ((sigma_xy) * (gammaxl) + (sigma_yy) * (gammayl) + (sigma_zy) * (gammazl));\n\
    s_tempz3[tx - 0] = (jacobianl) * ((sigma_xz) * (gammaxl) + (sigma_yz) * (gammayl) + (sigma_zz) * (gammazl));\n\
  }\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if(active){\n\
    tempx1l = 0.0f;\n\
    tempx2l = 0.0f;\n\
    tempx3l = 0.0f;\n\
    tempy1l = 0.0f;\n\
    tempy2l = 0.0f;\n\
    tempy3l = 0.0f;\n\
    tempz1l = 0.0f;\n\
    tempz2l = 0.0f;\n\
    tempz3l = 0.0f;\n\
#ifdef MANUALLY_UNROLLED_LOOPS\n\
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 0 - 0];\n\
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 0;\n\
    tempx1l = tempx1l + (s_tempx1[offset - 0]) * (fac1);\n\
    tempy1l = tempy1l + (s_tempy1[offset - 0]) * (fac1);\n\
    tempz1l = tempz1l + (s_tempz1[offset - 0]) * (fac1);\n\
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 0 - 0];\n\
    offset = (K) * (NGLL2) + (0) * (NGLLX) + I;\n\
    tempx2l = tempx2l + (s_tempx2[offset - 0]) * (fac2);\n\
    tempy2l = tempy2l + (s_tempy2[offset - 0]) * (fac2);\n\
    tempz2l = tempz2l + (s_tempz2[offset - 0]) * (fac2);\n\
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 0 - 0];\n\
    offset = (0) * (NGLL2) + (J) * (NGLLX) + I;\n\
    tempx3l = tempx3l + (s_tempx3[offset - 0]) * (fac3);\n\
    tempy3l = tempy3l + (s_tempy3[offset - 0]) * (fac3);\n\
    tempz3l = tempz3l + (s_tempz3[offset - 0]) * (fac3);\n\
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 1 - 0];\n\
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 1;\n\
    tempx1l = tempx1l + (s_tempx1[offset - 0]) * (fac1);\n\
    tempy1l = tempy1l + (s_tempy1[offset - 0]) * (fac1);\n\
    tempz1l = tempz1l + (s_tempz1[offset - 0]) * (fac1);\n\
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 1 - 0];\n\
    offset = (K) * (NGLL2) + (1) * (NGLLX) + I;\n\
    tempx2l = tempx2l + (s_tempx2[offset - 0]) * (fac2);\n\
    tempy2l = tempy2l + (s_tempy2[offset - 0]) * (fac2);\n\
    tempz2l = tempz2l + (s_tempz2[offset - 0]) * (fac2);\n\
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 1 - 0];\n\
    offset = (1) * (NGLL2) + (J) * (NGLLX) + I;\n\
    tempx3l = tempx3l + (s_tempx3[offset - 0]) * (fac3);\n\
    tempy3l = tempy3l + (s_tempy3[offset - 0]) * (fac3);\n\
    tempz3l = tempz3l + (s_tempz3[offset - 0]) * (fac3);\n\
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 2 - 0];\n\
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 2;\n\
    tempx1l = tempx1l + (s_tempx1[offset - 0]) * (fac1);\n\
    tempy1l = tempy1l + (s_tempy1[offset - 0]) * (fac1);\n\
    tempz1l = tempz1l + (s_tempz1[offset - 0]) * (fac1);\n\
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 2 - 0];\n\
    offset = (K) * (NGLL2) + (2) * (NGLLX) + I;\n\
    tempx2l = tempx2l + (s_tempx2[offset - 0]) * (fac2);\n\
    tempy2l = tempy2l + (s_tempy2[offset - 0]) * (fac2);\n\
    tempz2l = tempz2l + (s_tempz2[offset - 0]) * (fac2);\n\
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 2 - 0];\n\
    offset = (2) * (NGLL2) + (J) * (NGLLX) + I;\n\
    tempx3l = tempx3l + (s_tempx3[offset - 0]) * (fac3);\n\
    tempy3l = tempy3l + (s_tempy3[offset - 0]) * (fac3);\n\
    tempz3l = tempz3l + (s_tempz3[offset - 0]) * (fac3);\n\
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 3 - 0];\n\
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 3;\n\
    tempx1l = tempx1l + (s_tempx1[offset - 0]) * (fac1);\n\
    tempy1l = tempy1l + (s_tempy1[offset - 0]) * (fac1);\n\
    tempz1l = tempz1l + (s_tempz1[offset - 0]) * (fac1);\n\
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 3 - 0];\n\
    offset = (K) * (NGLL2) + (3) * (NGLLX) + I;\n\
    tempx2l = tempx2l + (s_tempx2[offset - 0]) * (fac2);\n\
    tempy2l = tempy2l + (s_tempy2[offset - 0]) * (fac2);\n\
    tempz2l = tempz2l + (s_tempz2[offset - 0]) * (fac2);\n\
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 3 - 0];\n\
    offset = (3) * (NGLL2) + (J) * (NGLLX) + I;\n\
    tempx3l = tempx3l + (s_tempx3[offset - 0]) * (fac3);\n\
    tempy3l = tempy3l + (s_tempy3[offset - 0]) * (fac3);\n\
    tempz3l = tempz3l + (s_tempz3[offset - 0]) * (fac3);\n\
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 4 - 0];\n\
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 4;\n\
    tempx1l = tempx1l + (s_tempx1[offset - 0]) * (fac1);\n\
    tempy1l = tempy1l + (s_tempy1[offset - 0]) * (fac1);\n\
    tempz1l = tempz1l + (s_tempz1[offset - 0]) * (fac1);\n\
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 4 - 0];\n\
    offset = (K) * (NGLL2) + (4) * (NGLLX) + I;\n\
    tempx2l = tempx2l + (s_tempx2[offset - 0]) * (fac2);\n\
    tempy2l = tempy2l + (s_tempy2[offset - 0]) * (fac2);\n\
    tempz2l = tempz2l + (s_tempz2[offset - 0]) * (fac2);\n\
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 4 - 0];\n\
    offset = (4) * (NGLL2) + (J) * (NGLLX) + I;\n\
    tempx3l = tempx3l + (s_tempx3[offset - 0]) * (fac3);\n\
    tempy3l = tempy3l + (s_tempy3[offset - 0]) * (fac3);\n\
    tempz3l = tempz3l + (s_tempz3[offset - 0]) * (fac3);\n\
#else\n\
    for(l=0; l<=NGLLX - (1); l+=1){\n\
      fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + l - 0];\n\
      offset = (K) * (NGLL2) + (J) * (NGLLX) + l;\n\
      tempx1l = tempx1l + (s_tempx1[offset - 0]) * (fac1);\n\
      tempy1l = tempy1l + (s_tempy1[offset - 0]) * (fac1);\n\
      tempz1l = tempz1l + (s_tempz1[offset - 0]) * (fac1);\n\
      fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + l - 0];\n\
      offset = (K) * (NGLL2) + (l) * (NGLLX) + I;\n\
      tempx2l = tempx2l + (s_tempx2[offset - 0]) * (fac2);\n\
      tempy2l = tempy2l + (s_tempy2[offset - 0]) * (fac2);\n\
      tempz2l = tempz2l + (s_tempz2[offset - 0]) * (fac2);\n\
      fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + l - 0];\n\
      offset = (l) * (NGLL2) + (J) * (NGLLX) + I;\n\
      tempx3l = tempx3l + (s_tempx3[offset - 0]) * (fac3);\n\
      tempy3l = tempy3l + (s_tempy3[offset - 0]) * (fac3);\n\
      tempz3l = tempz3l + (s_tempz3[offset - 0]) * (fac3);\n\
    }\n\
#endif\n\
    fac1 = d_wgllwgll_yz[(K) * (NGLLX) + J - 0];\n\
    fac2 = d_wgllwgll_xz[(K) * (NGLLX) + I - 0];\n\
    fac3 = d_wgllwgll_xy[(J) * (NGLLX) + I - 0];\n\
    sum_terms1 =  -((fac1) * (tempx1l) + (fac2) * (tempx2l) + (fac3) * (tempx3l));\n\
    sum_terms2 =  -((fac1) * (tempy1l) + (fac2) * (tempy2l) + (fac3) * (tempy3l));\n\
    sum_terms3 =  -((fac1) * (tempz1l) + (fac2) * (tempz2l) + (fac3) * (tempz3l));\n\
    if(GRAVITY){\n\
      sum_terms1 = sum_terms1 + rho_s_H1;\n\
      sum_terms2 = sum_terms2 + rho_s_H2;\n\
      sum_terms3 = sum_terms3 + rho_s_H3;\n\
    }\n\
#ifdef USE_MESH_COLORING_GPU\n\
#ifdef USE_TEXTURES_FIELDS\n\
    d_accel[0 - 0 + (iglob - (0)) * (3)] = as_float(read_imageui(d_b_accel_ic_tex, sampler_d_b_accel_ic_tex, int2((iglob) * (3) + 0,0)).x) + sum_terms1;\n\
    d_accel[1 - 0 + (iglob - (0)) * (3)] = as_float(read_imageui(d_b_accel_ic_tex, sampler_d_b_accel_ic_tex, int2((iglob) * (3) + 1,0)).x) + sum_terms2;\n\
    d_accel[2 - 0 + (iglob - (0)) * (3)] = as_float(read_imageui(d_b_accel_ic_tex, sampler_d_b_accel_ic_tex, int2((iglob) * (3) + 2,0)).x) + sum_terms3;\n\
#else\n\
    d_accel[0 - 0 + (iglob - (0)) * (3)] = d_accel[0 - 0 + (iglob - (0)) * (3)] + sum_terms1;\n\
    d_accel[1 - 0 + (iglob - (0)) * (3)] = d_accel[1 - 0 + (iglob - (0)) * (3)] + sum_terms2;\n\
    d_accel[2 - 0 + (iglob - (0)) * (3)] = d_accel[2 - 0 + (iglob - (0)) * (3)] + sum_terms3;\n\
#endif\n\
#else\n\
    if(use_mesh_coloring_gpu){\n\
      if(NSPEC_INNER_CORE > 1000){\n\
#ifdef USE_TEXTURES_FIELDS\n\
        d_accel[0 - 0 + (iglob - (0)) * (3)] = as_float(read_imageui(d_b_accel_ic_tex, sampler_d_b_accel_ic_tex, int2((iglob) * (3) + 0,0)).x) + sum_terms1;\n\
        d_accel[1 - 0 + (iglob - (0)) * (3)] = as_float(read_imageui(d_b_accel_ic_tex, sampler_d_b_accel_ic_tex, int2((iglob) * (3) + 1,0)).x) + sum_terms2;\n\
        d_accel[2 - 0 + (iglob - (0)) * (3)] = as_float(read_imageui(d_b_accel_ic_tex, sampler_d_b_accel_ic_tex, int2((iglob) * (3) + 2,0)).x) + sum_terms3;\n\
#else\n\
        d_accel[0 - 0 + (iglob - (0)) * (3)] = d_accel[0 - 0 + (iglob - (0)) * (3)] + sum_terms1;\n\
        d_accel[1 - 0 + (iglob - (0)) * (3)] = d_accel[1 - 0 + (iglob - (0)) * (3)] + sum_terms2;\n\
        d_accel[2 - 0 + (iglob - (0)) * (3)] = d_accel[2 - 0 + (iglob - (0)) * (3)] + sum_terms3;\n\
#endif\n\
      } else {\n\
        atomicAdd(d_accel + (iglob) * (3) + 0, sum_terms1);\n\
        atomicAdd(d_accel + (iglob) * (3) + 1, sum_terms2);\n\
        atomicAdd(d_accel + (iglob) * (3) + 2, sum_terms3);\n\
      }\n\
    } else {\n\
      atomicAdd(d_accel + (iglob) * (3) + 0, sum_terms1);\n\
      atomicAdd(d_accel + (iglob) * (3) + 1, sum_terms2);\n\
      atomicAdd(d_accel + (iglob) * (3) + 2, sum_terms3);\n\
    }\n\
#endif\n\
    if(ATTENUATION &&  ! PARTIAL_PHYS_DISPERSION_ONLY){\n\
      compute_element_ic_att_memory(tx, working_element, d_muvstore, factor_common, alphaval, betaval, gammaval, R_xx, R_yy, R_xy, R_xz, R_yz, epsilondev_xx, epsilondev_yy, epsilondev_xy, epsilondev_xz, epsilondev_yz, epsilondev_xx_loc, epsilondev_yy_loc, epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc, USE_3D_ATTENUATION_ARRAYS);\n\
    }\n\
    if(COMPUTE_AND_STORE_STRAIN){\n\
      epsilondev_xx[tx + (working_element) * (NGLL3) - 0] = epsilondev_xx_loc;\n\
      epsilondev_yy[tx + (working_element) * (NGLL3) - 0] = epsilondev_yy_loc;\n\
      epsilondev_xy[tx + (working_element) * (NGLL3) - 0] = epsilondev_xy_loc;\n\
      epsilondev_xz[tx + (working_element) * (NGLL3) - 0] = epsilondev_xz_loc;\n\
      epsilondev_yz[tx + (working_element) * (NGLL3) - 0] = epsilondev_yz_loc;\n\
    }\n\
  }\n\
}\n\
";
