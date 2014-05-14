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
static __device__ void compute_element_cm_att_stress(const int tx, const int working_element, const float * R_xx, const float * R_yy, const float * R_xy, const float * R_xz, const float * R_yz, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_xz, float * sigma_yz){
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
static __device__ void compute_element_cm_gravity(const int tx, const int iglob, const float * __restrict__ d_xstore, const float * __restrict__ d_ystore, const float * __restrict__ d_zstore, const float * __restrict__ d_minus_gravity_table, const float * __restrict__ d_minus_deriv_gravity_table, const float * __restrict__ d_density_table, const float * __restrict__ wgll_cube, const float jacobianl, const float * s_dummyx_loc, const float * s_dummyy_loc, const float * s_dummyz_loc, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_yx, float * sigma_xz, float * sigma_zx, float * sigma_yz, float * sigma_zy, float * rho_s_H1, float * rho_s_H2, float * rho_s_H3){
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
static __device__ void compute_element_cm_att_memory(const int tx, const int working_element, const float * d_muv, const float * factor_common, const float * alphaval, const float * betaval, const float * gammaval, float * R_xx, float * R_yy, float * R_xy, float * R_xz, float * R_yz, const float * epsilondev_xx, const float * epsilondev_yy, const float * epsilondev_xy, const float * epsilondev_xz, const float * epsilondev_yz, const float epsilondev_xx_loc, const float epsilondev_yy_loc, const float epsilondev_xy_loc, const float epsilondev_xz_loc, const float epsilondev_yz_loc, const float * d_c44store, const int ANISOTROPY, const int USE_3D_ATTENUATION_ARRAYS){
  int offset;
  int i_sls;
  float mul;
  float alphaval_loc;
  float betaval_loc;
  float gammaval_loc;
  float factor_loc;
  float sn;
  float snp1;
  if(ANISOTROPY){
    mul = d_c44store[tx + (NGLL3_PADDED) * (working_element) - 0];
  } else {
    mul = d_muv[tx + (NGLL3_PADDED) * (working_element) - 0];
  }
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
static __device__ void compute_element_cm_aniso(const int offset, const float * d_c11store, const float * d_c12store, const float * d_c13store, const float * d_c14store, const float * d_c15store, const float * d_c16store, const float * d_c22store, const float * d_c23store, const float * d_c24store, const float * d_c25store, const float * d_c26store, const float * d_c33store, const float * d_c34store, const float * d_c35store, const float * d_c36store, const float * d_c44store, const float * d_c45store, const float * d_c46store, const float * d_c55store, const float * d_c56store, const float * d_c66store, const int ATTENUATION, const float one_minus_sum_beta_use, const float duxdxl, const float duxdyl, const float duxdzl, const float duydxl, const float duydyl, const float duydzl, const float duzdxl, const float duzdyl, const float duzdzl, const float duxdyl_plus_duydxl, const float duzdxl_plus_duxdzl, const float duzdyl_plus_duydzl, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_xz, float * sigma_yz){
  float c11;
  float c12;
  float c13;
  float c14;
  float c15;
  float c16;
  float c22;
  float c23;
  float c24;
  float c25;
  float c26;
  float c33;
  float c34;
  float c35;
  float c36;
  float c44;
  float c45;
  float c46;
  float c55;
  float c56;
  float c66;
  float mul;
  float minus_sum_beta;
  c11 = d_c11store[offset - 0];
  c12 = d_c12store[offset - 0];
  c13 = d_c13store[offset - 0];
  c14 = d_c14store[offset - 0];
  c15 = d_c15store[offset - 0];
  c16 = d_c16store[offset - 0];
  c22 = d_c22store[offset - 0];
  c23 = d_c23store[offset - 0];
  c24 = d_c24store[offset - 0];
  c25 = d_c25store[offset - 0];
  c26 = d_c26store[offset - 0];
  c33 = d_c33store[offset - 0];
  c34 = d_c34store[offset - 0];
  c35 = d_c35store[offset - 0];
  c36 = d_c36store[offset - 0];
  c44 = d_c44store[offset - 0];
  c45 = d_c45store[offset - 0];
  c46 = d_c46store[offset - 0];
  c55 = d_c55store[offset - 0];
  c56 = d_c56store[offset - 0];
  c66 = d_c66store[offset - 0];
  if(ATTENUATION){
    minus_sum_beta = one_minus_sum_beta_use - (1.0f);
    mul = (c44) * (minus_sum_beta);
    c11 = c11 + (mul) * (1.3333333333333333f);
    c12 = c12 - ((mul) * (0.6666666666666666f));
    c13 = c13 - ((mul) * (0.6666666666666666f));
    c22 = c22 + (mul) * (1.3333333333333333f);
    c23 = c23 - ((mul) * (0.6666666666666666f));
    c33 = c33 + (mul) * (1.3333333333333333f);
    c44 = c44 + mul;
    c55 = c55 + mul;
    c66 = c66 + mul;
  }
  *(sigma_xx) = (c11) * (duxdxl) + (c16) * (duxdyl_plus_duydxl) + (c12) * (duydyl) + (c15) * (duzdxl_plus_duxdzl) + (c14) * (duzdyl_plus_duydzl) + (c13) * (duzdzl);
  *(sigma_yy) = (c12) * (duxdxl) + (c26) * (duxdyl_plus_duydxl) + (c22) * (duydyl) + (c25) * (duzdxl_plus_duxdzl) + (c24) * (duzdyl_plus_duydzl) + (c23) * (duzdzl);
  *(sigma_zz) = (c13) * (duxdxl) + (c36) * (duxdyl_plus_duydxl) + (c23) * (duydyl) + (c35) * (duzdxl_plus_duxdzl) + (c34) * (duzdyl_plus_duydzl) + (c33) * (duzdzl);
  *(sigma_xy) = (c16) * (duxdxl) + (c66) * (duxdyl_plus_duydxl) + (c26) * (duydyl) + (c56) * (duzdxl_plus_duxdzl) + (c46) * (duzdyl_plus_duydzl) + (c36) * (duzdzl);
  *(sigma_xz) = (c15) * (duxdxl) + (c56) * (duxdyl_plus_duydxl) + (c25) * (duydyl) + (c55) * (duzdxl_plus_duxdzl) + (c45) * (duzdyl_plus_duydzl) + (c35) * (duzdzl);
  *(sigma_yz) = (c14) * (duxdxl) + (c46) * (duxdyl_plus_duydxl) + (c24) * (duydyl) + (c45) * (duzdxl_plus_duxdzl) + (c44) * (duzdyl_plus_duydzl) + (c34) * (duzdzl);
}
static __device__ void compute_element_cm_iso(const int offset, const float * d_kappavstore, const float * d_muvstore, const int ATTENUATION, const float one_minus_sum_beta_use, const float duxdxl, const float duydyl, const float duzdzl, const float duxdxl_plus_duydyl, const float duxdxl_plus_duzdzl, const float duydyl_plus_duzdzl, const float duxdyl_plus_duydxl, const float duzdxl_plus_duxdzl, const float duzdyl_plus_duydzl, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_xz, float * sigma_yz){
  float lambdal;
  float mul;
  float lambdalplus2mul;
  float kappal;
  kappal = d_kappavstore[offset - 0];
  mul = d_muvstore[offset - 0];
  if(ATTENUATION){
    mul = (mul) * (one_minus_sum_beta_use);
  }
  lambdalplus2mul = kappal + (mul) * (1.3333333333333333f);
  lambdal = lambdalplus2mul - ((mul) * (2.0f));
  *(sigma_xx) = (lambdalplus2mul) * (duxdxl) + (lambdal) * (duydyl_plus_duzdzl);
  *(sigma_yy) = (lambdalplus2mul) * (duydyl) + (lambdal) * (duxdxl_plus_duzdzl);
  *(sigma_zz) = (lambdalplus2mul) * (duzdzl) + (lambdal) * (duxdxl_plus_duydyl);
  *(sigma_xy) = (mul) * (duxdyl_plus_duydxl);
  *(sigma_xz) = (mul) * (duzdxl_plus_duxdzl);
  *(sigma_yz) = (mul) * (duzdyl_plus_duydzl);
}
static __device__ void compute_element_cm_tiso(const int offset, const float * d_kappavstore, const float * d_muvstore, const float * d_kappahstore, const float * d_muhstore, const float * d_eta_anisostore, const int ATTENUATION, const float one_minus_sum_beta_use, const float duxdxl, const float duxdyl, const float duxdzl, const float duydxl, const float duydyl, const float duydzl, const float duzdxl, const float duzdyl, const float duzdzl, const float duxdyl_plus_duydxl, const float duzdxl_plus_duxdzl, const float duzdyl_plus_duydzl, const int iglob, const float * d_ystore, const float * d_zstore, float * sigma_xx, float * sigma_yy, float * sigma_zz, float * sigma_xy, float * sigma_xz, float * sigma_yz){
  float kappavl;
  float muvl;
  float kappahl;
  float muhl;
  float rhovpvsq;
  float rhovphsq;
  float rhovsvsq;
  float rhovshsq;
  float eta_aniso;
  float costheta;
  float sintheta;
  float cosphi;
  float sinphi;
  float costhetasq;
  float sinthetasq;
  float cosphisq;
  float sinphisq;
  float costhetafour;
  float sinthetafour;
  float cosphifour;
  float sinphifour;
  float costwotheta;
  float sintwotheta;
  float costwophi;
  float sintwophi;
  float cosfourtheta;
  float cosfourphi;
  float costwothetasq;
  float costwophisq;
  float sintwophisq;
  float etaminone;
  float twoetaminone;
  float two_eta_aniso;
  float four_eta_aniso;
  float six_eta_aniso;
  float two_rhovsvsq;
  float two_rhovshsq;
  float four_rhovsvsq;
  float four_rhovshsq;
  float c11;
  float c12;
  float c13;
  float c14;
  float c15;
  float c16;
  float c22;
  float c23;
  float c24;
  float c25;
  float c26;
  float c33;
  float c34;
  float c35;
  float c36;
  float c44;
  float c45;
  float c46;
  float c55;
  float c56;
  float c66;
  float theta;
  float phi;
  kappavl = d_kappavstore[offset - 0];
  muvl = d_muvstore[offset - 0];
  kappahl = d_kappahstore[offset - 0];
  muhl = d_muhstore[offset - 0];
  if(ATTENUATION){
    muvl = (muvl) * (one_minus_sum_beta_use);
    muhl = (muhl) * (one_minus_sum_beta_use);
  }
  rhovpvsq = kappavl + (muvl) * (1.3333333333333333f);
  rhovphsq = kappahl + (muhl) * (1.3333333333333333f);
  rhovsvsq = muvl;
  rhovshsq = muhl;
  eta_aniso = d_eta_anisostore[offset - 0];
  theta = d_ystore[iglob - 0];
  phi = d_zstore[iglob - 0];
  sincosf(theta,  &sintheta,  &costheta);
  sincosf(phi,  &sinphi,  &cosphi);
  sincosf((theta) * (2.0f),  &sintwotheta,  &costwotheta);
  sincosf((phi) * (2.0f),  &sintwophi,  &costwophi);
  cosfourtheta = cosf((theta) * (4.0f));
  cosfourphi = cosf((phi) * (4.0f));
  costhetasq = (costheta) * (costheta);
  sinthetasq = (sintheta) * (sintheta);
  cosphisq = (cosphi) * (cosphi);
  sinphisq = (sinphi) * (sinphi);
  costhetafour = (costhetasq) * (costhetasq);
  sinthetafour = (sinthetasq) * (sinthetasq);
  cosphifour = (cosphisq) * (cosphisq);
  sinphifour = (sinphisq) * (sinphisq);
  costwothetasq = (costwotheta) * (costwotheta);
  costwophisq = (costwophi) * (costwophi);
  sintwophisq = (sintwophi) * (sintwophi);
  etaminone = eta_aniso - (1.0f);
  twoetaminone = (eta_aniso) * (2.0f) - (1.0f);
  two_eta_aniso = (eta_aniso) * (2.0f);
  four_eta_aniso = (eta_aniso) * (4.0f);
  six_eta_aniso = (eta_aniso) * (6.0f);
  two_rhovsvsq = (rhovsvsq) * (2.0f);
  two_rhovshsq = (rhovshsq) * (2.0f);
  four_rhovsvsq = (rhovsvsq) * (4.0f);
  four_rhovshsq = (rhovshsq) * (4.0f);
  c11 = (rhovphsq) * (sinphifour) + (((cosphisq) * (sinphisq)) * ((rhovphsq) * (costhetasq) + ((eta_aniso) * (rhovphsq) + two_rhovsvsq - ((two_eta_aniso) * (rhovsvsq))) * (sinthetasq))) * (2.0f) + (cosphifour) * ((rhovphsq) * (costhetafour) + ((((eta_aniso) * (rhovphsq) + two_rhovsvsq - ((two_eta_aniso) * (rhovsvsq))) * (costhetasq)) * (sinthetasq)) * (2.0f) + (rhovpvsq) * (sinthetafour));
  c12 = (((rhovphsq - (two_rhovshsq)) * (cosfourphi + 3.0f)) * (costhetasq)) * (0.25f) - ((((four_rhovshsq) * (cosphisq)) * (costhetasq)) * (sinphisq)) + (((rhovphsq) * ((costwotheta) * (4.0f) + cosfourtheta + 11.0f)) * (sintwophisq)) * (0.03125f) + (((eta_aniso) * (rhovphsq - (two_rhovsvsq))) * (cosphifour + (((cosphisq) * (costhetasq)) * (sinphisq)) * (2.0f) + sinphifour)) * (sinthetasq) + (((rhovpvsq) * (cosphisq)) * (sinphisq)) * (sinthetafour) - (((rhovsvsq) * (sintwophisq)) * (sinthetafour));
  c13 = ((cosphisq) * (rhovphsq + (six_eta_aniso) * (rhovphsq) + rhovpvsq - (four_rhovsvsq) - (((eta_aniso) * (rhovsvsq)) * (12.0f)) + ((twoetaminone) * (rhovphsq) - (rhovpvsq) + four_rhovsvsq - ((four_eta_aniso) * (rhovsvsq))) * (cosfourtheta))) * (0.125f) + (sinphisq) * (((eta_aniso) * (rhovphsq - (two_rhovsvsq))) * (costhetasq) + (rhovphsq - (two_rhovshsq)) * (sinthetasq));
  c14 = (((costheta) * (sinphi)) * (((cosphisq) * ( -(rhovphsq) + rhovpvsq + four_rhovshsq - (four_rhovsvsq) + ( -(rhovphsq) + (two_eta_aniso) * (rhovphsq) - (rhovpvsq) + four_rhovsvsq - ((four_eta_aniso) * (rhovsvsq))) * (costwotheta))) * (0.5f) + ((etaminone) * (rhovphsq) + (rhovshsq - ((eta_aniso) * (rhovsvsq))) * (2.0f)) * (sinphisq))) * (sintheta);
  c15 = (((cosphi) * (costheta)) * (((cosphisq) * ( -(rhovphsq) + rhovpvsq + ((twoetaminone) * (rhovphsq) - (rhovpvsq) + four_rhovsvsq - ((four_eta_aniso) * (rhovsvsq))) * (costwotheta))) * (0.5f) + ((etaminone) * (rhovphsq - (two_rhovsvsq))) * (sinphisq))) * (sintheta);
  c16 = ((((cosphi) * (sinphi)) * ((cosphisq) * ( -(rhovphsq) + rhovpvsq + ( -(rhovphsq) + (two_eta_aniso) * (rhovphsq) - (rhovpvsq) + four_rhovsvsq - ((four_eta_aniso) * (rhovsvsq))) * (costwotheta)) + (((etaminone) * (rhovphsq - (two_rhovsvsq))) * (sinphisq)) * (2.0f))) * (sinthetasq)) * (0.5f);
  c22 = (rhovphsq) * (cosphifour) + (((cosphisq) * (sinphisq)) * ((rhovphsq) * (costhetasq) + ((eta_aniso) * (rhovphsq) + two_rhovsvsq - ((two_eta_aniso) * (rhovsvsq))) * (sinthetasq))) * (2.0f) + (sinphifour) * ((rhovphsq) * (costhetafour) + ((((eta_aniso) * (rhovphsq) + two_rhovsvsq - ((two_eta_aniso) * (rhovsvsq))) * (costhetasq)) * (sinthetasq)) * (2.0f) + (rhovpvsq) * (sinthetafour));
  c23 = ((rhovphsq + (six_eta_aniso) * (rhovphsq) + rhovpvsq - (four_rhovsvsq) - (((eta_aniso) * (rhovsvsq)) * (12.0f)) + ((twoetaminone) * (rhovphsq) - (rhovpvsq) + four_rhovsvsq - ((four_eta_aniso) * (rhovsvsq))) * (cosfourtheta)) * (sinphisq)) * (0.125f) + (cosphisq) * (((eta_aniso) * (rhovphsq - (two_rhovsvsq))) * (costhetasq) + (rhovphsq - (two_rhovshsq)) * (sinthetasq));
  c24 = (((costheta) * (sinphi)) * (((etaminone) * (rhovphsq - (two_rhovsvsq))) * (cosphisq) + (( -(rhovphsq) + rhovpvsq + ((twoetaminone) * (rhovphsq) - (rhovpvsq) + four_rhovsvsq - ((four_eta_aniso) * (rhovsvsq))) * (costwotheta)) * (sinphisq)) * (0.5f))) * (sintheta);
  c25 = (((cosphi) * (costheta)) * (((etaminone) * (rhovphsq) + (rhovshsq - ((eta_aniso) * (rhovsvsq))) * (2.0f)) * (cosphisq) + (( -(rhovphsq) + rhovpvsq + four_rhovshsq - (four_rhovsvsq) + ( -(rhovphsq) + (two_eta_aniso) * (rhovphsq) - (rhovpvsq) + four_rhovsvsq - ((four_eta_aniso) * (rhovsvsq))) * (costwotheta)) * (sinphisq)) * (0.5f))) * (sintheta);
  c26 = ((((cosphi) * (sinphi)) * ((((etaminone) * (rhovphsq - (two_rhovsvsq))) * (cosphisq)) * (2.0f) + ( -(rhovphsq) + rhovpvsq + ( -(rhovphsq) + (two_eta_aniso) * (rhovphsq) - (rhovpvsq) + four_rhovsvsq - ((four_eta_aniso) * (rhovsvsq))) * (costwotheta)) * (sinphisq))) * (sinthetasq)) * (0.5f);
  c33 = (rhovpvsq) * (costhetafour) + ((((eta_aniso) * (rhovphsq - (two_rhovsvsq)) + two_rhovsvsq) * (costhetasq)) * (sinthetasq)) * (2.0f) + (rhovphsq) * (sinthetafour);
  c34 = ( -(((rhovphsq - (rhovpvsq) + ((twoetaminone) * (rhovphsq) - (rhovpvsq) + four_rhovsvsq - ((four_eta_aniso) * (rhovsvsq))) * (costwotheta)) * (sinphi)) * (sintwotheta))) * (0.25f);
  c35 = ( -(((cosphi) * (rhovphsq - (rhovpvsq) + ((twoetaminone) * (rhovphsq) - (rhovpvsq) + four_rhovsvsq - ((four_eta_aniso) * (rhovsvsq))) * (costwotheta))) * (sintwotheta))) * (0.25f);
  c36 = ( -(((rhovphsq - (rhovpvsq) - (four_rhovshsq) + four_rhovsvsq + ((twoetaminone) * (rhovphsq) - (rhovpvsq) + four_rhovsvsq - ((four_eta_aniso) * (rhovsvsq))) * (costwotheta)) * (sintwophi)) * (sinthetasq))) * (0.25f);
  c44 = (cosphisq) * ((rhovsvsq) * (costhetasq) + (rhovshsq) * (sinthetasq)) + (sinphisq) * ((rhovsvsq) * (costwothetasq) + ((rhovphsq - ((two_eta_aniso) * (rhovphsq)) + rhovpvsq + (four_eta_aniso) * (rhovsvsq)) * (costhetasq)) * (sinthetasq));
  c45 = (((rhovphsq - ((two_eta_aniso) * (rhovphsq)) + rhovpvsq - (two_rhovshsq) - (two_rhovsvsq) + (four_eta_aniso) * (rhovsvsq) + (rhovphsq - ((two_eta_aniso) * (rhovphsq)) + rhovpvsq + ((etaminone) * (rhovsvsq)) * (4.0f)) * (costwotheta)) * (sintwophi)) * (sinthetasq)) * (0.25f);
  c46 =  -((((cosphi) * (costheta)) * ((rhovshsq - (rhovsvsq)) * (cosphisq) - (((rhovphsq - ((two_eta_aniso) * (rhovphsq)) + rhovpvsq - (two_rhovshsq) - (two_rhovsvsq) + (four_eta_aniso) * (rhovsvsq) + ( -(rhovphsq) + (two_eta_aniso) * (rhovphsq) - (rhovpvsq) + four_rhovsvsq - ((four_eta_aniso) * (rhovsvsq))) * (costwotheta)) * (sinphisq)) * (0.5f)))) * (sintheta));
  c55 = (sinphisq) * ((rhovsvsq) * (costhetasq) + (rhovshsq) * (sinthetasq)) + (cosphisq) * ((rhovsvsq) * (costwothetasq) + ((rhovphsq - ((two_eta_aniso) * (rhovphsq)) + rhovpvsq + (four_eta_aniso) * (rhovsvsq)) * (costhetasq)) * (sinthetasq));
  c56 = (((costheta) * (sinphi)) * (((cosphisq) * (rhovphsq - ((two_eta_aniso) * (rhovphsq)) + rhovpvsq - (two_rhovshsq) - (two_rhovsvsq) + (four_eta_aniso) * (rhovsvsq) + ( -(rhovphsq) + (two_eta_aniso) * (rhovphsq) - (rhovpvsq) + four_rhovsvsq - ((four_eta_aniso) * (rhovsvsq))) * (costwotheta))) * (0.5f) + ( -(rhovshsq) + rhovsvsq) * (sinphisq))) * (sintheta);
  c66 = ((rhovshsq) * (costwophisq)) * (costhetasq) - (((((rhovphsq - (two_rhovshsq)) * (cosphisq)) * (costhetasq)) * (sinphisq)) * (2.0f)) + (((rhovphsq) * ((costwotheta) * (4.0f) + cosfourtheta + 11.0f)) * (sintwophisq)) * (0.03125f) - ((((rhovsvsq) * ((cosfourphi) * (-2.0f) + cos((phi) * (4.0f) - ((theta) * (2.0f))) - ((costwotheta) * (2.0f)) + cos(((phi) * (2.0f) + theta) * (2.0f)) - (6.0f))) * (sinthetasq)) * (0.125f)) + (((rhovpvsq) * (cosphisq)) * (sinphisq)) * (sinthetafour) - (((((eta_aniso) * (rhovphsq - (two_rhovsvsq))) * (sintwophisq)) * (sinthetafour)) * (0.5f));
  *(sigma_xx) = (c11) * (duxdxl) + (c16) * (duxdyl_plus_duydxl) + (c12) * (duydyl) + (c15) * (duzdxl_plus_duxdzl) + (c14) * (duzdyl_plus_duydzl) + (c13) * (duzdzl);
  *(sigma_yy) = (c12) * (duxdxl) + (c26) * (duxdyl_plus_duydxl) + (c22) * (duydyl) + (c25) * (duzdxl_plus_duxdzl) + (c24) * (duzdyl_plus_duydzl) + (c23) * (duzdzl);
  *(sigma_zz) = (c13) * (duxdxl) + (c36) * (duxdyl_plus_duydxl) + (c23) * (duydyl) + (c35) * (duzdxl_plus_duxdzl) + (c34) * (duzdyl_plus_duydzl) + (c33) * (duzdzl);
  *(sigma_xy) = (c16) * (duxdxl) + (c66) * (duxdyl_plus_duydxl) + (c26) * (duydyl) + (c56) * (duzdxl_plus_duxdzl) + (c46) * (duzdyl_plus_duydzl) + (c36) * (duzdzl);
  *(sigma_xz) = (c15) * (duxdxl) + (c56) * (duxdyl_plus_duydxl) + (c25) * (duydyl) + (c55) * (duzdxl_plus_duxdzl) + (c45) * (duzdyl_plus_duydzl) + (c35) * (duzdzl);
  *(sigma_yz) = (c14) * (duxdxl) + (c46) * (duxdyl_plus_duydxl) + (c24) * (duydyl) + (c45) * (duzdxl_plus_duxdzl) + (c44) * (duzdyl_plus_duydzl) + (c34) * (duzdzl);
}
__global__ void crust_mantle_impl_kernel_adjoint(const int nb_blocks_to_compute, const int * d_ibool, const int * d_ispec_is_tiso, const int * d_phase_ispec_inner, const int num_phase_ispec, const int d_iphase, const float deltat, const int use_mesh_coloring_gpu, const float * __restrict__ d_displ, float * d_accel, const float * __restrict__ d_xix, const float * __restrict__ d_xiy, const float * __restrict__ d_xiz, const float * __restrict__ d_etax, const float * __restrict__ d_etay, const float * __restrict__ d_etaz, const float * __restrict__ d_gammax, const float * __restrict__ d_gammay, const float * __restrict__ d_gammaz, const float * __restrict__ d_hprime_xx, const float * __restrict__ d_hprimewgll_xx, const float * __restrict__ d_wgllwgll_xy, const float * __restrict__ d_wgllwgll_xz, const float * __restrict__ d_wgllwgll_yz, const float * __restrict__ d_kappavstore, const float * __restrict__ d_muvstore, const float * __restrict__ d_kappahstore, const float * __restrict__ d_muhstore, const float * __restrict__ d_eta_anisostore, const int COMPUTE_AND_STORE_STRAIN, float * epsilondev_xx, float * epsilondev_yy, float * epsilondev_xy, float * epsilondev_xz, float * epsilondev_yz, float * epsilon_trace_over_3, const int ATTENUATION, const int PARTIAL_PHYS_DISPERSION_ONLY, const int USE_3D_ATTENUATION_ARRAYS, const float * __restrict__ one_minus_sum_beta, const float * __restrict__ factor_common, float * R_xx, float * R_yy, float * R_xy, float * R_xz, float * R_yz, const float * __restrict__ alphaval, const float * __restrict__ betaval, const float * __restrict__ gammaval, const int ANISOTROPY, const float * __restrict__ d_c11store, const float * __restrict__ d_c12store, const float * __restrict__ d_c13store, const float * __restrict__ d_c14store, const float * __restrict__ d_c15store, const float * __restrict__ d_c16store, const float * __restrict__ d_c22store, const float * __restrict__ d_c23store, const float * __restrict__ d_c24store, const float * __restrict__ d_c25store, const float * __restrict__ d_c26store, const float * __restrict__ d_c33store, const float * __restrict__ d_c34store, const float * __restrict__ d_c35store, const float * __restrict__ d_c36store, const float * __restrict__ d_c44store, const float * __restrict__ d_c45store, const float * __restrict__ d_c46store, const float * __restrict__ d_c55store, const float * __restrict__ d_c56store, const float * __restrict__ d_c66store, const int GRAVITY, const float * __restrict__ d_xstore, const float * __restrict__ d_ystore, const float * __restrict__ d_zstore, const float * __restrict__ d_minus_gravity_table, const float * __restrict__ d_minus_deriv_gravity_table, const float * __restrict__ d_density_table, const float * __restrict__ wgll_cube, const int NSPEC_CRUST_MANTLE_STRAIN_ONLY){
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
  float one_minus_sum_beta_use;
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
    iglob = d_ibool[(working_element) * (NGLL3) + tx - 0] - (1);
#ifdef USE_TEXTURES_FIELDS
    s_dummyx_loc[tx - 0] = tex1Dfetch(d_b_displ_cm_tex,(iglob) * (3) + 0);
    s_dummyy_loc[tx - 0] = tex1Dfetch(d_b_displ_cm_tex,(iglob) * (3) + 1);
    s_dummyz_loc[tx - 0] = tex1Dfetch(d_b_displ_cm_tex,(iglob) * (3) + 2);
#else
    s_dummyx_loc[tx - 0] = d_displ[0 - 0 + (iglob - (0)) * (3)];
    s_dummyy_loc[tx - 0] = d_displ[1 - 0 + (iglob - (0)) * (3)];
    s_dummyz_loc[tx - 0] = d_displ[2 - 0 + (iglob - (0)) * (3)];
#endif
  }
  if(tx < NGLL2){
#ifdef USE_TEXTURES_CONSTANTS
    sh_hprime_xx[tx - 0] = tex1Dfetch(d_hprime_xx_tex,tx);
    sh_hprimewgll_xx[tx - 0] = tex1Dfetch(d_hprimewgll_xx_tex,tx);
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
      if(NSPEC_CRUST_MANTLE_STRAIN_ONLY == 1){
        epsilon_trace_over_3[tx - 0] = templ;
      } else {
        epsilon_trace_over_3[tx + (working_element) * (NGLL3) - 0] = templ;
      }
    }
    if(USE_3D_ATTENUATION_ARRAYS){
      one_minus_sum_beta_use = one_minus_sum_beta[tx + (working_element) * (NGLL3) - 0];
    } else {
      one_minus_sum_beta_use = one_minus_sum_beta[working_element - 0];
    }
    if(ANISOTROPY){
      compute_element_cm_aniso(offset, d_c11store, d_c12store, d_c13store, d_c14store, d_c15store, d_c16store, d_c22store, d_c23store, d_c24store, d_c25store, d_c26store, d_c33store, d_c34store, d_c35store, d_c36store, d_c44store, d_c45store, d_c46store, d_c55store, d_c56store, d_c66store, ATTENUATION, one_minus_sum_beta_use, duxdxl, duxdyl, duxdzl, duydxl, duydyl, duydzl, duzdxl, duzdyl, duzdzl, duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl,  &sigma_xx,  &sigma_yy,  &sigma_zz,  &sigma_xy,  &sigma_xz,  &sigma_yz);
    } else {
      if( ! d_ispec_is_tiso[working_element - 0]){
        compute_element_cm_iso(offset, d_kappavstore, d_muvstore, ATTENUATION, one_minus_sum_beta_use, duxdxl, duydyl, duzdzl, duxdxl_plus_duydyl, duxdxl_plus_duzdzl, duydyl_plus_duzdzl, duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl,  &sigma_xx,  &sigma_yy,  &sigma_zz,  &sigma_xy,  &sigma_xz,  &sigma_yz);
      } else {
        compute_element_cm_tiso(offset, d_kappavstore, d_muvstore, d_kappahstore, d_muhstore, d_eta_anisostore, ATTENUATION, one_minus_sum_beta_use, duxdxl, duxdyl, duxdzl, duydxl, duydyl, duydzl, duzdxl, duzdyl, duzdzl, duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl, iglob, d_ystore, d_zstore,  &sigma_xx,  &sigma_yy,  &sigma_zz,  &sigma_xy,  &sigma_xz,  &sigma_yz);
      }
    }
    if(ATTENUATION &&  ! PARTIAL_PHYS_DISPERSION_ONLY){
      compute_element_cm_att_stress(tx, working_element, R_xx, R_yy, R_xy, R_xz, R_yz,  &sigma_xx,  &sigma_yy,  &sigma_zz,  &sigma_xy,  &sigma_xz,  &sigma_yz);
    }
    sigma_yx = sigma_xy;
    sigma_zx = sigma_xz;
    sigma_zy = sigma_yz;
    jacobianl = (1.0f) / ((xixl) * ((etayl) * (gammazl) - ((etazl) * (gammayl))) - ((xiyl) * ((etaxl) * (gammazl) - ((etazl) * (gammaxl)))) + (xizl) * ((etaxl) * (gammayl) - ((etayl) * (gammaxl))));
    if(GRAVITY){
      compute_element_cm_gravity(tx, iglob, d_xstore, d_ystore, d_zstore, d_minus_gravity_table, d_minus_deriv_gravity_table, d_density_table, wgll_cube, jacobianl, s_dummyx_loc, s_dummyy_loc, s_dummyz_loc,  &sigma_xx,  &sigma_yy,  &sigma_zz,  &sigma_xy,  &sigma_yx,  &sigma_xz,  &sigma_zx,  &sigma_yz,  &sigma_zy,  &rho_s_H1,  &rho_s_H2,  &rho_s_H3);
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
    d_accel[0 - 0 + (iglob - (0)) * (3)] = tex1Dfetch(d_b_accel_cm_tex,(iglob) * (3) + 0) + sum_terms1;
    d_accel[1 - 0 + (iglob - (0)) * (3)] = tex1Dfetch(d_b_accel_cm_tex,(iglob) * (3) + 1) + sum_terms2;
    d_accel[2 - 0 + (iglob - (0)) * (3)] = tex1Dfetch(d_b_accel_cm_tex,(iglob) * (3) + 2) + sum_terms3;
#else
    d_accel[0 - 0 + (iglob - (0)) * (3)] = d_accel[0 - 0 + (iglob - (0)) * (3)] + sum_terms1;
    d_accel[1 - 0 + (iglob - (0)) * (3)] = d_accel[1 - 0 + (iglob - (0)) * (3)] + sum_terms2;
    d_accel[2 - 0 + (iglob - (0)) * (3)] = d_accel[2 - 0 + (iglob - (0)) * (3)] + sum_terms3;
#endif
#else
    if(use_mesh_coloring_gpu){
#ifdef USE_TEXTURES_FIELDS
      d_accel[0 - 0 + (iglob - (0)) * (3)] = tex1Dfetch(d_b_accel_cm_tex,(iglob) * (3) + 0) + sum_terms1;
      d_accel[1 - 0 + (iglob - (0)) * (3)] = tex1Dfetch(d_b_accel_cm_tex,(iglob) * (3) + 1) + sum_terms2;
      d_accel[2 - 0 + (iglob - (0)) * (3)] = tex1Dfetch(d_b_accel_cm_tex,(iglob) * (3) + 2) + sum_terms3;
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
#endif
    if(ATTENUATION &&  ! PARTIAL_PHYS_DISPERSION_ONLY){
      compute_element_cm_att_memory(tx, working_element, d_muvstore, factor_common, alphaval, betaval, gammaval, R_xx, R_yy, R_xy, R_xz, R_yz, epsilondev_xx, epsilondev_yy, epsilondev_xy, epsilondev_xz, epsilondev_yz, epsilondev_xx_loc, epsilondev_yy_loc, epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc, d_c44store, ANISOTROPY, USE_3D_ATTENUATION_ARRAYS);
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
