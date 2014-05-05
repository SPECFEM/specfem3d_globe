#ifndef INDEX2
#define INDEX2(xsize,x,y) x + (y)*xsize
#endif
#ifndef INDEX3
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)
#endif
#ifndef INDEX4
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))
#endif
#ifndef INDEX5
#define INDEX5(xsize,ysize,zsize,isize,x,y,z,i,j) x + xsize*(y + ysize*(z + zsize*(i + isize*(j))))
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
static __device__ void compute_element_ic_att_stress(const int tx, const int working_element, const float * R_xx, const float * R_yy, const float * R_xy, const float * R_xz, const float * R_yz, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_xz, float * sigma_yz){
  int offset;
  int i_sls;
  float R_xx_val;
  float R_yy_val;
  for(i_sls=0; i_sls<=N_SLS - (1); i_sls+=1){
    offset = i_sls + (N_SLS) * (tx + (NGLL3) * (working_element));
    R_xx_val = R_xx[offset - 0];
    R_yy_val = R_yy[offset - 0];
    sigma_xx[0 - 0] = sigma_xx[0 - 0] - (R_xx_val);
    sigma_yy[0 - 0] = sigma_yy[0 - 0] - (R_yy_val);
    sigma_zz[0 - 0] = sigma_zz[0 - 0] + R_xx_val + R_yy_val;
    sigma_xy[0 - 0] = sigma_xy[0 - 0] - (R_xy[offset - 0]);
    sigma_xz[0 - 0] = sigma_xz[0 - 0] - (R_xz[offset - 0]);
    sigma_yz[0 - 0] = sigma_yz[0 - 0] - (R_yz[offset - 0]);
  }
}
static __device__ void compute_element_ic_gravity(const int tx, const int iglob, const float * d_xstore, const float * d_ystore, const float * d_zstore, const float * d_minus_gravity_table, const float * d_minus_deriv_gravity_table, const float * d_density_table, const float * wgll_cube, const float jacobianl, const float * s_dummyx_loc, const float * s_dummyy_loc, const float * s_dummyz_loc, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_yx, float * sigma_xz, float * sigma_zx, float * sigma_yz, float * sigma_zy, float * rho_s_H1, float * rho_s_H2, float * rho_s_H3){
  float radius;
  float theta;
  float phi;
  float cos_theta;
  float sin_theta;
  float cos_phi;
  float sin_phi;
  float cos_theta_sq;
  float sin_theta_sq;
  float cos_phi_sq;
  float sin_phi_sq;
  float minus_g;
  float minus_dg;
  float rho;
  float gxl;
  float gyl;
  float gzl;
  float minus_g_over_radius;
  float minus_dg_plus_g_over_radius;
  float Hxxl;
  float Hyyl;
  float Hzzl;
  float Hxyl;
  float Hxzl;
  float Hyzl;
  float sx_l;
  float sy_l;
  float sz_l;
  float factor;
  int int_radius;
  radius = d_xstore[iglob - 0];
  if(radius < 1.5696123057604773e-05f){
    radius = 1.5696123057604773e-05f;
  }
  theta = d_ystore[iglob - 0];
  phi = d_zstore[iglob - 0];
  sincosf(theta,  &sin_theta,  &cos_theta);
  sincosf(phi,  &sin_phi,  &cos_phi);
  int_radius = rint(((radius) * (6371.0f)) * (10.0f)) - (1);
  if(int_radius < 0){
    int_radius = 0;
  }
  minus_g = d_minus_gravity_table[int_radius - 0];
  minus_dg = d_minus_deriv_gravity_table[int_radius - 0];
  rho = d_density_table[int_radius - 0];
  gxl = ((minus_g) * (sin_theta)) * (cos_phi);
  gyl = ((minus_g) * (sin_theta)) * (sin_phi);
  gzl = (minus_g) * (cos_theta);
  minus_g_over_radius = (minus_g) / (radius);
  minus_dg_plus_g_over_radius = minus_dg - (minus_g_over_radius);
  cos_theta_sq = (cos_theta) * (cos_theta);
  sin_theta_sq = (sin_theta) * (sin_theta);
  cos_phi_sq = (cos_phi) * (cos_phi);
  sin_phi_sq = (sin_phi) * (sin_phi);
  Hxxl = (minus_g_over_radius) * ((cos_phi_sq) * (cos_theta_sq) + sin_phi_sq) + ((cos_phi_sq) * (minus_dg)) * (sin_theta_sq);
  Hyyl = (minus_g_over_radius) * (cos_phi_sq + (cos_theta_sq) * (sin_phi_sq)) + ((minus_dg) * (sin_phi_sq)) * (sin_theta_sq);
  Hzzl = (cos_theta_sq) * (minus_dg) + (minus_g_over_radius) * (sin_theta_sq);
  Hxyl = (((cos_phi) * (minus_dg_plus_g_over_radius)) * (sin_phi)) * (sin_theta_sq);
  Hxzl = (((cos_phi) * (cos_theta)) * (minus_dg_plus_g_over_radius)) * (sin_theta);
  Hyzl = (((cos_theta) * (minus_dg_plus_g_over_radius)) * (sin_phi)) * (sin_theta);
  sx_l = (rho) * (s_dummyx_loc[tx - 0]);
  sy_l = (rho) * (s_dummyy_loc[tx - 0]);
  sz_l = (rho) * (s_dummyz_loc[tx - 0]);
  *(sigma_xx) = *(sigma_xx) + (sy_l) * (gyl) + (sz_l) * (gzl);
  *(sigma_yy) = *(sigma_yy) + (sx_l) * (gxl) + (sz_l) * (gzl);
  *(sigma_zz) = *(sigma_zz) + (sx_l) * (gxl) + (sy_l) * (gyl);
  *(sigma_xy) = *(sigma_xy) - ((sx_l) * (gyl));
  *(sigma_yx) = *(sigma_yx) - ((sy_l) * (gxl));
  *(sigma_xz) = *(sigma_xz) - ((sx_l) * (gzl));
  *(sigma_zx) = *(sigma_zx) - ((sz_l) * (gxl));
  *(sigma_yz) = *(sigma_yz) - ((sy_l) * (gzl));
  *(sigma_zy) = *(sigma_zy) - ((sz_l) * (gyl));
  factor = (jacobianl) * (wgll_cube[tx - 0]);
  rho_s_H1[0 - 0] = (factor) * ((sx_l) * (Hxxl) + (sy_l) * (Hxyl) + (sz_l) * (Hxzl));
  rho_s_H2[0 - 0] = (factor) * ((sx_l) * (Hxyl) + (sy_l) * (Hyyl) + (sz_l) * (Hyzl));
  rho_s_H3[0 - 0] = (factor) * ((sx_l) * (Hxzl) + (sy_l) * (Hyzl) + (sz_l) * (Hzzl));
}
static __device__ void compute_element_ic_att_memory(const int tx, const int working_element, const float * d_muv, const float * factor_common, const float * alphaval, const float * betaval, const float * gammaval, float * R_xx, float * R_yy, float * R_xy, float * R_xz, float * R_yz, const float * epsilondev_xx, const float * epsilondev_yy, const float * epsilondev_xy, const float * epsilondev_xz, const float * epsilondev_yz, const float epsilondev_xx_loc, const float epsilondev_yy_loc, const float epsilondev_xy_loc, const float epsilondev_xz_loc, const float epsilondev_yz_loc, const int USE_3D_ATTENUATION_ARRAYS){
  int offset;
  int i_sls;
  float mul;
  float alphaval_loc;
  float betaval_loc;
  float gammaval_loc;
  float factor_loc;
  float sn;
  float snp1;
  mul = d_muv[tx + (NGLL3_PADDED) * (working_element) - 0];
  for(i_sls=0; i_sls<=N_SLS - (1); i_sls+=1){
    offset = i_sls + (N_SLS) * (tx + (NGLL3) * (working_element));
    if(USE_3D_ATTENUATION_ARRAYS){
      factor_loc = (mul) * (factor_common[offset - 0]);
    } else {
      factor_loc = (mul) * (factor_common[i_sls + (N_SLS) * (working_element) - 0]);
    }
    alphaval_loc = alphaval[i_sls - 0];
    betaval_loc = betaval[i_sls - 0];
    gammaval_loc = gammaval[i_sls - 0];
    sn = (factor_loc) * (epsilondev_xx[tx + (NGLL3) * (working_element) - 0]);
    snp1 = (factor_loc) * (epsilondev_xx_loc);
    R_xx[offset - 0] = (alphaval_loc) * (R_xx[offset - 0]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);
    sn = (factor_loc) * (epsilondev_yy[tx + (NGLL3) * (working_element) - 0]);
    snp1 = (factor_loc) * (epsilondev_yy_loc);
    R_yy[offset - 0] = (alphaval_loc) * (R_yy[offset - 0]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);
    sn = (factor_loc) * (epsilondev_xy[tx + (NGLL3) * (working_element) - 0]);
    snp1 = (factor_loc) * (epsilondev_xy_loc);
    R_xy[offset - 0] = (alphaval_loc) * (R_xy[offset - 0]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);
    sn = (factor_loc) * (epsilondev_xz[tx + (NGLL3) * (working_element) - 0]);
    snp1 = (factor_loc) * (epsilondev_xz_loc);
    R_xz[offset - 0] = (alphaval_loc) * (R_xz[offset - 0]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);
    sn = (factor_loc) * (epsilondev_yz[tx + (NGLL3) * (working_element) - 0]);
    snp1 = (factor_loc) * (epsilondev_yz_loc);
    R_yz[offset - 0] = (alphaval_loc) * (R_yz[offset - 0]) + (betaval_loc) * (sn) + (gammaval_loc) * (snp1);
  }
}
__global__ void inner_core_impl_kernel_adjoint(const int nb_blocks_to_compute, const int * d_ibool, const int * d_idoubling, const int * d_phase_ispec_inner, const int num_phase_ispec, const int d_iphase, const float deltat, const int use_mesh_coloring_gpu, const float * d_displ, float * d_accel, const float * d_xix, const float * d_xiy, const float * d_xiz, const float * d_etax, const float * d_etay, const float * d_etaz, const float * d_gammax, const float * d_gammay, const float * d_gammaz, const float * d_hprime_xx, const float * d_hprimewgll_xx, const float * d_wgllwgll_xy, const float * d_wgllwgll_xz, const float * d_wgllwgll_yz, const float * d_kappavstore, const float * d_muvstore, const int COMPUTE_AND_STORE_STRAIN, float * epsilondev_xx, float * epsilondev_yy, float * epsilondev_xy, float * epsilondev_xz, float * epsilondev_yz, float * epsilon_trace_over_3, const int ATTENUATION, const int PARTIAL_PHYS_DISPERSION_ONLY, const int USE_3D_ATTENUATION_ARRAYS, const float * one_minus_sum_beta, const float * factor_common, float * R_xx, float * R_yy, float * R_xy, float * R_xz, float * R_yz, const float * alphaval, const float * betaval, const float * gammaval, const int ANISOTROPY, const float * d_c11store, const float * d_c12store, const float * d_c13store, const float * d_c33store, const float * d_c44store, const int GRAVITY, const float * d_xstore, const float * d_ystore, const float * d_zstore, const float * d_minus_gravity_table, const float * d_minus_deriv_gravity_table, const float * d_density_table, const float * wgll_cube, const int NSPEC_INNER_CORE_STRAIN_ONLY, const int NSPEC_INNER_CORE){
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
  float tempx1l;
  float tempx2l;
  float tempx3l;
  float tempy1l;
  float tempy2l;
  float tempy3l;
  float tempz1l;
  float tempz2l;
  float tempz3l;
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
  float duxdxl;
  float duxdyl;
  float duxdzl;
  float duydxl;
  float duydyl;
  float duydzl;
  float duzdxl;
  float duzdyl;
  float duzdzl;
  float duxdxl_plus_duydyl;
  float duxdxl_plus_duzdzl;
  float duydyl_plus_duzdzl;
  float duxdyl_plus_duydxl;
  float duzdxl_plus_duxdzl;
  float duzdyl_plus_duydzl;
  float templ;
  float fac1;
  float fac2;
  float fac3;
  float lambdal;
  float mul;
  float lambdalplus2mul;
  float kappal;
  float mul_iso;
  float mul_aniso;
  float sigma_xx;
  float sigma_xy;
  float sigma_xz;
  float sigma_yx;
  float sigma_yy;
  float sigma_yz;
  float sigma_zx;
  float sigma_zy;
  float sigma_zz;
  float epsilondev_xx_loc;
  float epsilondev_yy_loc;
  float epsilondev_xy_loc;
  float epsilondev_xz_loc;
  float epsilondev_yz_loc;
  float c11;
  float c12;
  float c13;
  float c33;
  float c44;
  float sum_terms1;
  float sum_terms2;
  float sum_terms3;
  float rho_s_H1;
  float rho_s_H2;
  float rho_s_H3;
  __shared__ float s_dummyx_loc[NGLL3 + 0 - (1) - (0) + 1];
  __shared__ float s_dummyy_loc[NGLL3 + 0 - (1) - (0) + 1];
  __shared__ float s_dummyz_loc[NGLL3 + 0 - (1) - (0) + 1];
  __shared__ float s_tempx1[NGLL3 + 0 - (1) - (0) + 1];
  __shared__ float s_tempx2[NGLL3 + 0 - (1) - (0) + 1];
  __shared__ float s_tempx3[NGLL3 + 0 - (1) - (0) + 1];
  __shared__ float s_tempy1[NGLL3 + 0 - (1) - (0) + 1];
  __shared__ float s_tempy2[NGLL3 + 0 - (1) - (0) + 1];
  __shared__ float s_tempy3[NGLL3 + 0 - (1) - (0) + 1];
  __shared__ float s_tempz1[NGLL3 + 0 - (1) - (0) + 1];
  __shared__ float s_tempz2[NGLL3 + 0 - (1) - (0) + 1];
  __shared__ float s_tempz3[NGLL3 + 0 - (1) - (0) + 1];
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
    if(d_idoubling[working_element - 0] == IFLAG_IN_FICTITIOUS_CUBE){
      active = 0;
    } else {
      iglob = d_ibool[(working_element) * (NGLL3) + tx - 0] - (1);
#ifdef USE_TEXTURES_FIELDS
      s_dummyx_loc[tx - 0] = tex1Dfetch(d_b_displ_ic_tex,(iglob) * (3) + 0);
      s_dummyy_loc[tx - 0] = tex1Dfetch(d_b_displ_ic_tex,(iglob) * (3) + 1);
      s_dummyz_loc[tx - 0] = tex1Dfetch(d_b_displ_ic_tex,(iglob) * (3) + 2);
#else
      s_dummyx_loc[tx - 0] = d_displ[0 - 0 + (iglob - (0)) * (3)];
      s_dummyy_loc[tx - 0] = d_displ[1 - 0 + (iglob - (0)) * (3)];
      s_dummyz_loc[tx - 0] = d_displ[2 - 0 + (iglob - (0)) * (3)];
#endif
    }
  }
  if(tx < NGLL2){
#ifdef USE_TEXTURES_CONSTANTS
    sh_hprime_xx[tx - 0] = tex1Dfetch(d_hprime_xx_ic_tex,tx);
    sh_hprimewgll_xx[tx - 0] = tex1Dfetch(d_hprimewgll_xx_ic_tex,tx);
#else
    sh_hprime_xx[tx - 0] = d_hprime_xx[tx - 0];
    sh_hprimewgll_xx[tx - 0] = d_hprimewgll_xx[tx - 0];
#endif
  }
  __syncthreads();
  if(active){
    tempx1l = 0.0f;
    tempx2l = 0.0f;
    tempx3l = 0.0f;
    tempy1l = 0.0f;
    tempy2l = 0.0f;
    tempy3l = 0.0f;
    tempz1l = 0.0f;
    tempz2l = 0.0f;
    tempz3l = 0.0f;
#ifdef MANUALLY_UNROLLED_LOOPS
    fac1 = sh_hprime_xx[(0) * (NGLLX) + I - 0];
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 0 - 0]) * (fac1);
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 0 - 0]) * (fac1);
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 0 - 0]) * (fac1);
    fac2 = sh_hprime_xx[(0) * (NGLLX) + J - 0];
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (0) * (NGLLX) + I - 0]) * (fac2);
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (0) * (NGLLX) + I - 0]) * (fac2);
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (0) * (NGLLX) + I - 0]) * (fac2);
    fac3 = sh_hprime_xx[(0) * (NGLLX) + K - 0];
    tempx3l = tempx3l + (s_dummyx_loc[(0) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);
    tempy3l = tempy3l + (s_dummyy_loc[(0) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);
    tempz3l = tempz3l + (s_dummyz_loc[(0) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);
    fac1 = sh_hprime_xx[(1) * (NGLLX) + I - 0];
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 1 - 0]) * (fac1);
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 1 - 0]) * (fac1);
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 1 - 0]) * (fac1);
    fac2 = sh_hprime_xx[(1) * (NGLLX) + J - 0];
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (1) * (NGLLX) + I - 0]) * (fac2);
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (1) * (NGLLX) + I - 0]) * (fac2);
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (1) * (NGLLX) + I - 0]) * (fac2);
    fac3 = sh_hprime_xx[(1) * (NGLLX) + K - 0];
    tempx3l = tempx3l + (s_dummyx_loc[(1) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);
    tempy3l = tempy3l + (s_dummyy_loc[(1) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);
    tempz3l = tempz3l + (s_dummyz_loc[(1) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);
    fac1 = sh_hprime_xx[(2) * (NGLLX) + I - 0];
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 2 - 0]) * (fac1);
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 2 - 0]) * (fac1);
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 2 - 0]) * (fac1);
    fac2 = sh_hprime_xx[(2) * (NGLLX) + J - 0];
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (2) * (NGLLX) + I - 0]) * (fac2);
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (2) * (NGLLX) + I - 0]) * (fac2);
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (2) * (NGLLX) + I - 0]) * (fac2);
    fac3 = sh_hprime_xx[(2) * (NGLLX) + K - 0];
    tempx3l = tempx3l + (s_dummyx_loc[(2) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);
    tempy3l = tempy3l + (s_dummyy_loc[(2) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);
    tempz3l = tempz3l + (s_dummyz_loc[(2) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);
    fac1 = sh_hprime_xx[(3) * (NGLLX) + I - 0];
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 3 - 0]) * (fac1);
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 3 - 0]) * (fac1);
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 3 - 0]) * (fac1);
    fac2 = sh_hprime_xx[(3) * (NGLLX) + J - 0];
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (3) * (NGLLX) + I - 0]) * (fac2);
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (3) * (NGLLX) + I - 0]) * (fac2);
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (3) * (NGLLX) + I - 0]) * (fac2);
    fac3 = sh_hprime_xx[(3) * (NGLLX) + K - 0];
    tempx3l = tempx3l + (s_dummyx_loc[(3) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);
    tempy3l = tempy3l + (s_dummyy_loc[(3) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);
    tempz3l = tempz3l + (s_dummyz_loc[(3) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);
    fac1 = sh_hprime_xx[(4) * (NGLLX) + I - 0];
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + 4 - 0]) * (fac1);
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + 4 - 0]) * (fac1);
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + 4 - 0]) * (fac1);
    fac2 = sh_hprime_xx[(4) * (NGLLX) + J - 0];
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (4) * (NGLLX) + I - 0]) * (fac2);
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (4) * (NGLLX) + I - 0]) * (fac2);
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (4) * (NGLLX) + I - 0]) * (fac2);
    fac3 = sh_hprime_xx[(4) * (NGLLX) + K - 0];
    tempx3l = tempx3l + (s_dummyx_loc[(4) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);
    tempy3l = tempy3l + (s_dummyy_loc[(4) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);
    tempz3l = tempz3l + (s_dummyz_loc[(4) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);
#else
    for(l=0; l<=NGLLX - (1); l+=1){
      fac1 = sh_hprime_xx[(l) * (NGLLX) + I - 0];
      tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + l - 0]) * (fac1);
      tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + l - 0]) * (fac1);
      tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + l - 0]) * (fac1);
      fac2 = sh_hprime_xx[(l) * (NGLLX) + J - 0];
      tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (l) * (NGLLX) + I - 0]) * (fac2);
      tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (l) * (NGLLX) + I - 0]) * (fac2);
      tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (l) * (NGLLX) + I - 0]) * (fac2);
      fac3 = sh_hprime_xx[(l) * (NGLLX) + K - 0];
      tempx3l = tempx3l + (s_dummyx_loc[(l) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);
      tempy3l = tempy3l + (s_dummyy_loc[(l) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);
      tempz3l = tempz3l + (s_dummyz_loc[(l) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (fac3);
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
    duxdxl = (xixl) * (tempx1l) + (etaxl) * (tempx2l) + (gammaxl) * (tempx3l);
    duxdyl = (xiyl) * (tempx1l) + (etayl) * (tempx2l) + (gammayl) * (tempx3l);
    duxdzl = (xizl) * (tempx1l) + (etazl) * (tempx2l) + (gammazl) * (tempx3l);
    duydxl = (xixl) * (tempy1l) + (etaxl) * (tempy2l) + (gammaxl) * (tempy3l);
    duydyl = (xiyl) * (tempy1l) + (etayl) * (tempy2l) + (gammayl) * (tempy3l);
    duydzl = (xizl) * (tempy1l) + (etazl) * (tempy2l) + (gammazl) * (tempy3l);
    duzdxl = (xixl) * (tempz1l) + (etaxl) * (tempz2l) + (gammaxl) * (tempz3l);
    duzdyl = (xiyl) * (tempz1l) + (etayl) * (tempz2l) + (gammayl) * (tempz3l);
    duzdzl = (xizl) * (tempz1l) + (etazl) * (tempz2l) + (gammazl) * (tempz3l);
    duxdxl_plus_duydyl = duxdxl + duydyl;
    duxdxl_plus_duzdzl = duxdxl + duzdzl;
    duydyl_plus_duzdzl = duydyl + duzdzl;
    duxdyl_plus_duydxl = duxdyl + duydxl;
    duzdxl_plus_duxdzl = duzdxl + duxdzl;
    duzdyl_plus_duydzl = duzdyl + duydzl;
    if(COMPUTE_AND_STORE_STRAIN){
      templ = (duxdxl + duydyl + duzdzl) * (0.3333333333333333f);
      epsilondev_xx_loc = duxdxl - (templ);
      epsilondev_yy_loc = duydyl - (templ);
      epsilondev_xy_loc = (duxdyl_plus_duydxl) * (0.5f);
      epsilondev_xz_loc = (duzdxl_plus_duxdzl) * (0.5f);
      epsilondev_yz_loc = (duzdyl_plus_duydzl) * (0.5f);
      if(NSPEC_INNER_CORE_STRAIN_ONLY == 1){
        epsilon_trace_over_3[tx - 0] = templ;
      } else {
        epsilon_trace_over_3[tx + (working_element) * (NGLL3) - 0] = templ;
      }
    }
    kappal = d_kappavstore[offset - 0];
    mul = d_muvstore[offset - 0];
    if(ATTENUATION){
      if(USE_3D_ATTENUATION_ARRAYS){
        mul_iso = (mul) * (one_minus_sum_beta[tx + (working_element) * (NGLL3) - 0]);
        mul_aniso = (mul) * (one_minus_sum_beta[tx + (working_element) * (NGLL3) - 0] - (1.0f));
      } else {
        mul_iso = (mul) * (one_minus_sum_beta[working_element - 0]);
        mul_aniso = (mul) * (one_minus_sum_beta[working_element - 0] - (1.0f));
      }
    } else {
      mul_iso = mul;
    }
    if(ANISOTROPY){
      c11 = d_c11store[offset - 0];
      c12 = d_c12store[offset - 0];
      c13 = d_c13store[offset - 0];
      c33 = d_c33store[offset - 0];
      c44 = d_c44store[offset - 0];
      if(ATTENUATION){
        c11 = c11 + (mul_aniso) * (1.3333333333333333f);
        c12 = c12 - ((mul_aniso) * (0.6666666666666666f));
        c13 = c13 - ((mul_aniso) * (0.6666666666666666f));
        c33 = c33 + (mul_aniso) * (1.3333333333333333f);
        c44 = c44 + mul_aniso;
      }
      sigma_xx = (c11) * (duxdxl) + (c12) * (duydyl) + (c13) * (duzdzl);
      sigma_yy = (c12) * (duxdxl) + (c11) * (duydyl) + (c13) * (duzdzl);
      sigma_zz = (c13) * (duxdxl) + (c13) * (duydyl) + (c33) * (duzdzl);
      sigma_xy = ((c11 - (c12)) * (duxdyl_plus_duydxl)) * (0.5f);
      sigma_xz = (c44) * (duzdxl_plus_duxdzl);
      sigma_yz = (c44) * (duzdyl_plus_duydzl);
    } else {
      lambdalplus2mul = kappal + (mul_iso) * (1.3333333333333333f);
      lambdal = lambdalplus2mul - ((mul_iso) * (2.0f));
      sigma_xx = (lambdalplus2mul) * (duxdxl) + (lambdal) * (duydyl_plus_duzdzl);
      sigma_yy = (lambdalplus2mul) * (duydyl) + (lambdal) * (duxdxl_plus_duzdzl);
      sigma_zz = (lambdalplus2mul) * (duzdzl) + (lambdal) * (duxdxl_plus_duydyl);
      sigma_xy = (mul) * (duxdyl_plus_duydxl);
      sigma_xz = (mul) * (duzdxl_plus_duxdzl);
      sigma_yz = (mul) * (duzdyl_plus_duydzl);
    }
    if(ATTENUATION &&  ! PARTIAL_PHYS_DISPERSION_ONLY){
      compute_element_ic_att_stress(tx, working_element, R_xx, R_yy, R_xy, R_xz, R_yz,  &sigma_xx,  &sigma_yy,  &sigma_zz,  &sigma_xy,  &sigma_xz,  &sigma_yz);
    }
    sigma_yx = sigma_xy;
    sigma_zx = sigma_xz;
    sigma_zy = sigma_yz;
    jacobianl = (1.0f) / ((xixl) * ((etayl) * (gammazl) - ((etazl) * (gammayl))) - ((xiyl) * ((etaxl) * (gammazl) - ((etazl) * (gammaxl)))) + (xizl) * ((etaxl) * (gammayl) - ((etayl) * (gammaxl))));
    if(GRAVITY){
      compute_element_ic_gravity(tx, iglob, d_xstore, d_ystore, d_zstore, d_minus_gravity_table, d_minus_deriv_gravity_table, d_density_table, wgll_cube, jacobianl, s_dummyx_loc, s_dummyy_loc, s_dummyz_loc,  &sigma_xx,  &sigma_yy,  &sigma_zz,  &sigma_xy,  &sigma_yx,  &sigma_xz,  &sigma_zx,  &sigma_yz,  &sigma_zy,  &rho_s_H1,  &rho_s_H2,  &rho_s_H3);
    }
    s_tempx1[tx - 0] = (jacobianl) * ((sigma_xx) * (xixl) + (sigma_yx) * (xiyl) + (sigma_zx) * (xizl));
    s_tempy1[tx - 0] = (jacobianl) * ((sigma_xy) * (xixl) + (sigma_yy) * (xiyl) + (sigma_zy) * (xizl));
    s_tempz1[tx - 0] = (jacobianl) * ((sigma_xz) * (xixl) + (sigma_yz) * (xiyl) + (sigma_zz) * (xizl));
    s_tempx2[tx - 0] = (jacobianl) * ((sigma_xx) * (etaxl) + (sigma_yx) * (etayl) + (sigma_zx) * (etazl));
    s_tempy2[tx - 0] = (jacobianl) * ((sigma_xy) * (etaxl) + (sigma_yy) * (etayl) + (sigma_zy) * (etazl));
    s_tempz2[tx - 0] = (jacobianl) * ((sigma_xz) * (etaxl) + (sigma_yz) * (etayl) + (sigma_zz) * (etazl));
    s_tempx3[tx - 0] = (jacobianl) * ((sigma_xx) * (gammaxl) + (sigma_yx) * (gammayl) + (sigma_zx) * (gammazl));
    s_tempy3[tx - 0] = (jacobianl) * ((sigma_xy) * (gammaxl) + (sigma_yy) * (gammayl) + (sigma_zy) * (gammazl));
    s_tempz3[tx - 0] = (jacobianl) * ((sigma_xz) * (gammaxl) + (sigma_yz) * (gammayl) + (sigma_zz) * (gammazl));
  }
  __syncthreads();
  if(active){
    tempx1l = 0.0f;
    tempx2l = 0.0f;
    tempx3l = 0.0f;
    tempy1l = 0.0f;
    tempy2l = 0.0f;
    tempy3l = 0.0f;
    tempz1l = 0.0f;
    tempz2l = 0.0f;
    tempz3l = 0.0f;
#ifdef MANUALLY_UNROLLED_LOOPS
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 0 - 0];
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 0;
    tempx1l = tempx1l + (s_tempx1[offset - 0]) * (fac1);
    tempy1l = tempy1l + (s_tempy1[offset - 0]) * (fac1);
    tempz1l = tempz1l + (s_tempz1[offset - 0]) * (fac1);
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 0 - 0];
    offset = (K) * (NGLL2) + (0) * (NGLLX) + I;
    tempx2l = tempx2l + (s_tempx2[offset - 0]) * (fac2);
    tempy2l = tempy2l + (s_tempy2[offset - 0]) * (fac2);
    tempz2l = tempz2l + (s_tempz2[offset - 0]) * (fac2);
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 0 - 0];
    offset = (0) * (NGLL2) + (J) * (NGLLX) + I;
    tempx3l = tempx3l + (s_tempx3[offset - 0]) * (fac3);
    tempy3l = tempy3l + (s_tempy3[offset - 0]) * (fac3);
    tempz3l = tempz3l + (s_tempz3[offset - 0]) * (fac3);
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 1 - 0];
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 1;
    tempx1l = tempx1l + (s_tempx1[offset - 0]) * (fac1);
    tempy1l = tempy1l + (s_tempy1[offset - 0]) * (fac1);
    tempz1l = tempz1l + (s_tempz1[offset - 0]) * (fac1);
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 1 - 0];
    offset = (K) * (NGLL2) + (1) * (NGLLX) + I;
    tempx2l = tempx2l + (s_tempx2[offset - 0]) * (fac2);
    tempy2l = tempy2l + (s_tempy2[offset - 0]) * (fac2);
    tempz2l = tempz2l + (s_tempz2[offset - 0]) * (fac2);
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 1 - 0];
    offset = (1) * (NGLL2) + (J) * (NGLLX) + I;
    tempx3l = tempx3l + (s_tempx3[offset - 0]) * (fac3);
    tempy3l = tempy3l + (s_tempy3[offset - 0]) * (fac3);
    tempz3l = tempz3l + (s_tempz3[offset - 0]) * (fac3);
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 2 - 0];
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 2;
    tempx1l = tempx1l + (s_tempx1[offset - 0]) * (fac1);
    tempy1l = tempy1l + (s_tempy1[offset - 0]) * (fac1);
    tempz1l = tempz1l + (s_tempz1[offset - 0]) * (fac1);
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 2 - 0];
    offset = (K) * (NGLL2) + (2) * (NGLLX) + I;
    tempx2l = tempx2l + (s_tempx2[offset - 0]) * (fac2);
    tempy2l = tempy2l + (s_tempy2[offset - 0]) * (fac2);
    tempz2l = tempz2l + (s_tempz2[offset - 0]) * (fac2);
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 2 - 0];
    offset = (2) * (NGLL2) + (J) * (NGLLX) + I;
    tempx3l = tempx3l + (s_tempx3[offset - 0]) * (fac3);
    tempy3l = tempy3l + (s_tempy3[offset - 0]) * (fac3);
    tempz3l = tempz3l + (s_tempz3[offset - 0]) * (fac3);
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 3 - 0];
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 3;
    tempx1l = tempx1l + (s_tempx1[offset - 0]) * (fac1);
    tempy1l = tempy1l + (s_tempy1[offset - 0]) * (fac1);
    tempz1l = tempz1l + (s_tempz1[offset - 0]) * (fac1);
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 3 - 0];
    offset = (K) * (NGLL2) + (3) * (NGLLX) + I;
    tempx2l = tempx2l + (s_tempx2[offset - 0]) * (fac2);
    tempy2l = tempy2l + (s_tempy2[offset - 0]) * (fac2);
    tempz2l = tempz2l + (s_tempz2[offset - 0]) * (fac2);
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 3 - 0];
    offset = (3) * (NGLL2) + (J) * (NGLLX) + I;
    tempx3l = tempx3l + (s_tempx3[offset - 0]) * (fac3);
    tempy3l = tempy3l + (s_tempy3[offset - 0]) * (fac3);
    tempz3l = tempz3l + (s_tempz3[offset - 0]) * (fac3);
    fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + 4 - 0];
    offset = (K) * (NGLL2) + (J) * (NGLLX) + 4;
    tempx1l = tempx1l + (s_tempx1[offset - 0]) * (fac1);
    tempy1l = tempy1l + (s_tempy1[offset - 0]) * (fac1);
    tempz1l = tempz1l + (s_tempz1[offset - 0]) * (fac1);
    fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + 4 - 0];
    offset = (K) * (NGLL2) + (4) * (NGLLX) + I;
    tempx2l = tempx2l + (s_tempx2[offset - 0]) * (fac2);
    tempy2l = tempy2l + (s_tempy2[offset - 0]) * (fac2);
    tempz2l = tempz2l + (s_tempz2[offset - 0]) * (fac2);
    fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + 4 - 0];
    offset = (4) * (NGLL2) + (J) * (NGLLX) + I;
    tempx3l = tempx3l + (s_tempx3[offset - 0]) * (fac3);
    tempy3l = tempy3l + (s_tempy3[offset - 0]) * (fac3);
    tempz3l = tempz3l + (s_tempz3[offset - 0]) * (fac3);
#else
    for(l=0; l<=NGLLX - (1); l+=1){
      fac1 = sh_hprimewgll_xx[(I) * (NGLLX) + l - 0];
      offset = (K) * (NGLL2) + (J) * (NGLLX) + l;
      tempx1l = tempx1l + (s_tempx1[offset - 0]) * (fac1);
      tempy1l = tempy1l + (s_tempy1[offset - 0]) * (fac1);
      tempz1l = tempz1l + (s_tempz1[offset - 0]) * (fac1);
      fac2 = sh_hprimewgll_xx[(J) * (NGLLX) + l - 0];
      offset = (K) * (NGLL2) + (l) * (NGLLX) + I;
      tempx2l = tempx2l + (s_tempx2[offset - 0]) * (fac2);
      tempy2l = tempy2l + (s_tempy2[offset - 0]) * (fac2);
      tempz2l = tempz2l + (s_tempz2[offset - 0]) * (fac2);
      fac3 = sh_hprimewgll_xx[(K) * (NGLLX) + l - 0];
      offset = (l) * (NGLL2) + (J) * (NGLLX) + I;
      tempx3l = tempx3l + (s_tempx3[offset - 0]) * (fac3);
      tempy3l = tempy3l + (s_tempy3[offset - 0]) * (fac3);
      tempz3l = tempz3l + (s_tempz3[offset - 0]) * (fac3);
    }
#endif
    fac1 = d_wgllwgll_yz[(K) * (NGLLX) + J - 0];
    fac2 = d_wgllwgll_xz[(K) * (NGLLX) + I - 0];
    fac3 = d_wgllwgll_xy[(J) * (NGLLX) + I - 0];
    sum_terms1 =  -((fac1) * (tempx1l) + (fac2) * (tempx2l) + (fac3) * (tempx3l));
    sum_terms2 =  -((fac1) * (tempy1l) + (fac2) * (tempy2l) + (fac3) * (tempy3l));
    sum_terms3 =  -((fac1) * (tempz1l) + (fac2) * (tempz2l) + (fac3) * (tempz3l));
    if(GRAVITY){
      sum_terms1 = sum_terms1 + rho_s_H1;
      sum_terms2 = sum_terms2 + rho_s_H2;
      sum_terms3 = sum_terms3 + rho_s_H3;
    }
#ifdef USE_MESH_COLORING_GPU
#ifdef USE_TEXTURES_FIELDS
    d_accel[0 - 0 + (iglob - (0)) * (3)] = tex1Dfetch(d_b_accel_ic_tex,(iglob) * (3) + 0) + sum_terms1;
    d_accel[1 - 0 + (iglob - (0)) * (3)] = tex1Dfetch(d_b_accel_ic_tex,(iglob) * (3) + 1) + sum_terms2;
    d_accel[2 - 0 + (iglob - (0)) * (3)] = tex1Dfetch(d_b_accel_ic_tex,(iglob) * (3) + 2) + sum_terms3;
#else
    d_accel[0 - 0 + (iglob - (0)) * (3)] = d_accel[0 - 0 + (iglob - (0)) * (3)] + sum_terms1;
    d_accel[1 - 0 + (iglob - (0)) * (3)] = d_accel[1 - 0 + (iglob - (0)) * (3)] + sum_terms2;
    d_accel[2 - 0 + (iglob - (0)) * (3)] = d_accel[2 - 0 + (iglob - (0)) * (3)] + sum_terms3;
#endif
#else
    if(use_mesh_coloring_gpu){
      if(NSPEC_INNER_CORE > 1000){
#ifdef USE_TEXTURES_FIELDS
        d_accel[0 - 0 + (iglob - (0)) * (3)] = tex1Dfetch(d_b_accel_ic_tex,(iglob) * (3) + 0) + sum_terms1;
        d_accel[1 - 0 + (iglob - (0)) * (3)] = tex1Dfetch(d_b_accel_ic_tex,(iglob) * (3) + 1) + sum_terms2;
        d_accel[2 - 0 + (iglob - (0)) * (3)] = tex1Dfetch(d_b_accel_ic_tex,(iglob) * (3) + 2) + sum_terms3;
#else
        d_accel[0 - 0 + (iglob - (0)) * (3)] = d_accel[0 - 0 + (iglob - (0)) * (3)] + sum_terms1;
        d_accel[1 - 0 + (iglob - (0)) * (3)] = d_accel[1 - 0 + (iglob - (0)) * (3)] + sum_terms2;
        d_accel[2 - 0 + (iglob - (0)) * (3)] = d_accel[2 - 0 + (iglob - (0)) * (3)] + sum_terms3;
#endif
      } else {
        atomicAdd(d_accel + (iglob) * (3) + 0, sum_terms1);
        atomicAdd(d_accel + (iglob) * (3) + 1, sum_terms2);
        atomicAdd(d_accel + (iglob) * (3) + 2, sum_terms3);
      }
    } else {
      atomicAdd(d_accel + (iglob) * (3) + 0, sum_terms1);
      atomicAdd(d_accel + (iglob) * (3) + 1, sum_terms2);
      atomicAdd(d_accel + (iglob) * (3) + 2, sum_terms3);
    }
#endif
    if(ATTENUATION &&  ! PARTIAL_PHYS_DISPERSION_ONLY){
      compute_element_ic_att_memory(tx, working_element, d_muvstore, factor_common, alphaval, betaval, gammaval, R_xx, R_yy, R_xy, R_xz, R_yz, epsilondev_xx, epsilondev_yy, epsilondev_xy, epsilondev_xz, epsilondev_yz, epsilondev_xx_loc, epsilondev_yy_loc, epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc, USE_3D_ATTENUATION_ARRAYS);
    }
    if(COMPUTE_AND_STORE_STRAIN){
      epsilondev_xx[tx + (working_element) * (NGLL3) - 0] = epsilondev_xx_loc;
      epsilondev_yy[tx + (working_element) * (NGLL3) - 0] = epsilondev_yy_loc;
      epsilondev_xy[tx + (working_element) * (NGLL3) - 0] = epsilondev_xy_loc;
      epsilondev_xz[tx + (working_element) * (NGLL3) - 0] = epsilondev_xz_loc;
      epsilondev_yz[tx + (working_element) * (NGLL3) - 0] = epsilondev_yz_loc;
    }
  }
}
