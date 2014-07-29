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
static __device__ void compute_element_strain_undo_att(const int ispec, const int ijk_ispec, const int * d_ibool, const float * s_dummyx_loc, const float * s_dummyy_loc, const float * s_dummyz_loc, const float * d_xix, const float * d_xiy, const float * d_xiz, const float * d_etax, const float * d_etay, const float * d_etaz, const float * d_gammax, const float * d_gammay, const float * d_gammaz, const float * sh_hprime_xx, float * epsilondev_loc, float * epsilon_trace_over_3){
  int tx;
  int K;
  int J;
  int I;
  int l;
  int offset;
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
  float duxdxl;
  float duxdyl;
  float duxdzl;
  float duydxl;
  float duydyl;
  float duydzl;
  float duzdxl;
  float duzdyl;
  float duzdzl;
  float templ;
  float fac1;
  float fac2;
  float fac3;
  tx = threadIdx.x;
  K = (tx) / (NGLL2);
  J = (tx - ((K) * (NGLL2))) / (NGLLX);
  I = tx - ((K) * (NGLL2)) - ((J) * (NGLLX));
  tempx1l = 0.0f;
  tempx2l = 0.0f;
  tempx3l = 0.0f;
  tempy1l = 0.0f;
  tempy2l = 0.0f;
  tempy3l = 0.0f;
  tempz1l = 0.0f;
  tempz2l = 0.0f;
  tempz3l = 0.0f;
  for(l=0; l<=NGLLX - (1); l+=1){
    fac1 = sh_hprime_xx[(l) * (NGLLX) + I - (0)];
    tempx1l = tempx1l + (s_dummyx_loc[(K) * (NGLL2) + (J) * (NGLLX) + l - (0)]) * (fac1);
    tempy1l = tempy1l + (s_dummyy_loc[(K) * (NGLL2) + (J) * (NGLLX) + l - (0)]) * (fac1);
    tempz1l = tempz1l + (s_dummyz_loc[(K) * (NGLL2) + (J) * (NGLLX) + l - (0)]) * (fac1);
    fac2 = sh_hprime_xx[(l) * (NGLLX) + J - (0)];
    tempx2l = tempx2l + (s_dummyx_loc[(K) * (NGLL2) + (l) * (NGLLX) + I - (0)]) * (fac2);
    tempy2l = tempy2l + (s_dummyy_loc[(K) * (NGLL2) + (l) * (NGLLX) + I - (0)]) * (fac2);
    tempz2l = tempz2l + (s_dummyz_loc[(K) * (NGLL2) + (l) * (NGLLX) + I - (0)]) * (fac2);
    fac3 = sh_hprime_xx[(l) * (NGLLX) + K - (0)];
    tempx3l = tempx3l + (s_dummyx_loc[(l) * (NGLL2) + (J) * (NGLLX) + I - (0)]) * (fac3);
    tempy3l = tempy3l + (s_dummyy_loc[(l) * (NGLL2) + (J) * (NGLLX) + I - (0)]) * (fac3);
    tempz3l = tempz3l + (s_dummyz_loc[(l) * (NGLL2) + (J) * (NGLLX) + I - (0)]) * (fac3);
  }
  offset = (ispec) * (NGLL3_PADDED) + tx;
  xixl = d_xix[offset - (0)];
  etaxl = d_etax[offset - (0)];
  gammaxl = d_gammax[offset - (0)];
  xiyl = d_xiy[offset - (0)];
  etayl = d_etay[offset - (0)];
  gammayl = d_gammay[offset - (0)];
  xizl = d_xiz[offset - (0)];
  etazl = d_etaz[offset - (0)];
  gammazl = d_gammaz[offset - (0)];
  duxdxl = (xixl) * (tempx1l) + (etaxl) * (tempx2l) + (gammaxl) * (tempx3l);
  duxdyl = (xiyl) * (tempx1l) + (etayl) * (tempx2l) + (gammayl) * (tempx3l);
  duxdzl = (xizl) * (tempx1l) + (etazl) * (tempx2l) + (gammazl) * (tempx3l);
  duydxl = (xixl) * (tempy1l) + (etaxl) * (tempy2l) + (gammaxl) * (tempy3l);
  duydyl = (xiyl) * (tempy1l) + (etayl) * (tempy2l) + (gammayl) * (tempy3l);
  duydzl = (xizl) * (tempy1l) + (etazl) * (tempy2l) + (gammazl) * (tempy3l);
  duzdxl = (xixl) * (tempz1l) + (etaxl) * (tempz2l) + (gammaxl) * (tempz3l);
  duzdyl = (xiyl) * (tempz1l) + (etayl) * (tempz2l) + (gammayl) * (tempz3l);
  duzdzl = (xizl) * (tempz1l) + (etazl) * (tempz2l) + (gammazl) * (tempz3l);
  templ = (duxdxl + duydyl + duzdzl) * (0.3333333333333333f);
  epsilondev_loc[0 - (0)] = duxdxl - (templ);
  epsilondev_loc[1 - (0)] = duydyl - (templ);
  epsilondev_loc[2 - (0)] = (duxdyl + duydxl) * (0.5f);
  epsilondev_loc[3 - (0)] = (duzdxl + duxdzl) * (0.5f);
  epsilondev_loc[4 - (0)] = (duzdyl + duydzl) * (0.5f);
  *(epsilon_trace_over_3) = templ;
}
static __device__ void compute_strain_product(float * prod, const float eps_trace_over_3, const float * epsdev, const float b_eps_trace_over_3, const float * b_epsdev){
  float eps[6];
  float b_eps[6];
  int p;
  int i;
  int j;
  eps[0 - (0)] = epsdev[0 - (0)] + eps_trace_over_3;
  eps[1 - (0)] = epsdev[1 - (0)] + eps_trace_over_3;
  eps[2 - (0)] =  -(eps[0 - (0)] + eps[1 - (0)]) + (eps_trace_over_3) * (3.0f);
  eps[3 - (0)] = epsdev[4 - (0)];
  eps[4 - (0)] = epsdev[3 - (0)];
  eps[5 - (0)] = epsdev[2 - (0)];
  b_eps[0 - (0)] = b_epsdev[0 - (0)] + b_eps_trace_over_3;
  b_eps[1 - (0)] = b_epsdev[1 - (0)] + b_eps_trace_over_3;
  b_eps[2 - (0)] =  -(b_eps[0 - (0)] + b_eps[1 - (0)]) + (b_eps_trace_over_3) * (3.0f);
  b_eps[3 - (0)] = b_epsdev[4 - (0)];
  b_eps[4 - (0)] = b_epsdev[3 - (0)];
  b_eps[5 - (0)] = b_epsdev[2 - (0)];
  p = 0;
  for(i=0; i<=5; i+=1){
    for(j=0; j<=5; j+=1){
      prod[p - (0)] = (eps[i - (0)]) * (b_eps[j - (0)]);
      if(j > i){
        prod[p - (0)] = prod[p - (0)] + (eps[j - (0)]) * (b_eps[i - (0)]);
        if(j > 2 && i < 3){
          prod[p - (0)] = (prod[p - (0)]) * (2.0f);
        }
        if(i > 2){
          prod[p - (0)] = (prod[p - (0)]) * (4.0f);
        }
        p = p + 1;
      }
    }
  }
}
__global__ void compute_ani_undo_att_kernel(const float * epsilondev_xx, const float * epsilondev_yy, const float * epsilondev_xy, const float * epsilondev_xz, const float * epsilondev_yz, const float * epsilon_trace_over_3, float * cijkl_kl, const int NSPEC, const float deltat, const int * d_ibool, const float * d_b_displ, const float * d_xix, const float * d_xiy, const float * d_xiz, const float * d_etax, const float * d_etay, const float * d_etaz, const float * d_gammax, const float * d_gammay, const float * d_gammaz, const float * d_hprime_xx){
  int ispec;
  int ijk_ispec;
  int tx;
  int iglob;
  float eps_trace_over_3;
  float b_eps_trace_over_3;
  float prod[21];
  int i;
  float epsdev[5];
  float b_epsdev[5];
  __shared__ float s_dummyx_loc[NGLL3 + 0 - (1) - (0) + 1];
  __shared__ float s_dummyy_loc[NGLL3 + 0 - (1) - (0) + 1];
  __shared__ float s_dummyz_loc[NGLL3 + 0 - (1) - (0) + 1];
  __shared__ float sh_hprime_xx[NGLL2 + 0 - (1) - (0) + 1];
  ispec = blockIdx.x + (blockIdx.y) * (gridDim.x);
  ijk_ispec = threadIdx.x + (NGLL3) * (ispec);
  tx = threadIdx.x;
  if(tx < NGLL2){
    sh_hprime_xx[tx - (0)] = d_hprime_xx[tx - (0)];
  }
  if(ispec < NSPEC){
    iglob = d_ibool[ijk_ispec - (0)] - (1);
    s_dummyx_loc[tx - (0)] = d_b_displ[0 - (0) + (iglob - (0)) * (3)];
    s_dummyy_loc[tx - (0)] = d_b_displ[1 - (0) + (iglob - (0)) * (3)];
    s_dummyz_loc[tx - (0)] = d_b_displ[2 - (0) + (iglob - (0)) * (3)];
  }
  __syncthreads();
  if(ispec < NSPEC){
    epsdev[0 - (0)] = epsilondev_xx[ijk_ispec - (0)];
    epsdev[1 - (0)] = epsilondev_yy[ijk_ispec - (0)];
    epsdev[2 - (0)] = epsilondev_xy[ijk_ispec - (0)];
    epsdev[3 - (0)] = epsilondev_xz[ijk_ispec - (0)];
    epsdev[4 - (0)] = epsilondev_yz[ijk_ispec - (0)];
    eps_trace_over_3 = epsilon_trace_over_3[ijk_ispec - (0)];
    compute_element_strain_undo_att(ispec, ijk_ispec, d_ibool, s_dummyx_loc, s_dummyy_loc, s_dummyz_loc, d_xix, d_xiy, d_xiz, d_etax, d_etay, d_etaz, d_gammax, d_gammay, d_gammaz, sh_hprime_xx, b_epsdev,  &b_eps_trace_over_3);
    compute_strain_product(prod, eps_trace_over_3, epsdev, b_eps_trace_over_3, b_epsdev);
    for(i=0; i<=20; i+=1){
      cijkl_kl[i - (0) + (ijk_ispec - (0)) * (21)] = cijkl_kl[i - (0) + (ijk_ispec - (0)) * (21)] + (deltat) * (prod[i - (0)]);
    }
  }
}
