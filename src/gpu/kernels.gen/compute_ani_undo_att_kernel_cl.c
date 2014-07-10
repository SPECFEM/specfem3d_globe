const char * compute_ani_undo_att_kernel_program = "\
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
void compute_element_strain_undo_att(const int ispec, const int ijk_ispec, const __global int * d_ibool, const __local float * s_dummyx_loc, const __local float * s_dummyy_loc, const __local float * s_dummyz_loc, const __global float * d_xix, const __global float * d_xiy, const __global float * d_xiz, const __global float * d_etax, const __global float * d_etay, const __global float * d_etaz, const __global float * d_gammax, const __global float * d_gammay, const __global float * d_gammaz, const __local float * sh_hprime_xx, float * epsilondev_loc, float * epsilon_trace_over_3){\n\
  int tx;\n\
  int K;\n\
  int J;\n\
  int I;\n\
  int l;\n\
  int offset;\n\
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
  float duxdxl;\n\
  float duxdyl;\n\
  float duxdzl;\n\
  float duydxl;\n\
  float duydyl;\n\
  float duydzl;\n\
  float duzdxl;\n\
  float duzdyl;\n\
  float duzdzl;\n\
  float templ;\n\
  float fac1;\n\
  float fac2;\n\
  float fac3;\n\
  tx = get_local_id(0);\n\
  K = (tx) / (NGLL2);\n\
  J = (tx - ((K) * (NGLL2))) / (NGLLX);\n\
  I = tx - ((K) * (NGLL2)) - ((J) * (NGLLX));\n\
  tempx1l = 0.0f;\n\
  tempx2l = 0.0f;\n\
  tempx3l = 0.0f;\n\
  tempy1l = 0.0f;\n\
  tempy2l = 0.0f;\n\
  tempy3l = 0.0f;\n\
  tempz1l = 0.0f;\n\
  tempz2l = 0.0f;\n\
  tempz3l = 0.0f;\n\
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
  offset = (ispec) * (NGLL3_PADDED) + tx;\n\
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
  templ = (duxdxl + duydyl + duzdzl) * (0.3333333333333333f);\n\
  epsilondev_loc[0 - 0] = duxdxl - (templ);\n\
  epsilondev_loc[1 - 0] = duydyl - (templ);\n\
  epsilondev_loc[2 - 0] = (duxdyl + duydxl) * (0.5f);\n\
  epsilondev_loc[3 - 0] = (duzdxl + duxdzl) * (0.5f);\n\
  epsilondev_loc[4 - 0] = (duzdyl + duydzl) * (0.5f);\n\
  *(epsilon_trace_over_3) = templ;\n\
}\n\
void compute_strain_product(float * prod, const float eps_trace_over_3, const float * epsdev, const float b_eps_trace_over_3, const float * b_epsdev){\n\
  float eps[6];\n\
  float b_eps[6];\n\
  int p;\n\
  int i;\n\
  int j;\n\
  eps[0 - 0] = epsdev[0 - 0] + eps_trace_over_3;\n\
  eps[1 - 0] = epsdev[1 - 0] + eps_trace_over_3;\n\
  eps[2 - 0] =  -(eps[0 - 0] + eps[1 - 0]) + (eps_trace_over_3) * (3.0f);\n\
  eps[3 - 0] = epsdev[4 - 0];\n\
  eps[4 - 0] = epsdev[3 - 0];\n\
  eps[5 - 0] = epsdev[2 - 0];\n\
  b_eps[0 - 0] = b_epsdev[0 - 0] + b_eps_trace_over_3;\n\
  b_eps[1 - 0] = b_epsdev[1 - 0] + b_eps_trace_over_3;\n\
  b_eps[2 - 0] =  -(b_eps[0 - 0] + b_eps[1 - 0]) + (b_eps_trace_over_3) * (3.0f);\n\
  b_eps[3 - 0] = b_epsdev[4 - 0];\n\
  b_eps[4 - 0] = b_epsdev[3 - 0];\n\
  b_eps[5 - 0] = b_epsdev[2 - 0];\n\
  p = 0;\n\
  for(i=0; i<=5; i+=1){\n\
    for(j=0; j<=5; j+=1){\n\
      prod[p - 0] = (eps[i - 0]) * (b_eps[j - 0]);\n\
      if(j > i){\n\
        prod[p - 0] = prod[p - 0] + (eps[j - 0]) * (b_eps[i - 0]);\n\
        if(j > 2 && i < 3){\n\
          prod[p - 0] = (prod[p - 0]) * (2.0f);\n\
        }\n\
        if(i > 2){\n\
          prod[p - 0] = (prod[p - 0]) * (4.0f);\n\
        }\n\
        p = p + 1;\n\
      }\n\
    }\n\
  }\n\
}\n\
__kernel void compute_ani_undo_att_kernel(const __global float * epsilondev_xx, const __global float * epsilondev_yy, const __global float * epsilondev_xy, const __global float * epsilondev_xz, const __global float * epsilondev_yz, const __global float * epsilon_trace_over_3, __global float * cijkl_kl, const int NSPEC, const float deltat, const __global int * d_ibool, const __global float * d_b_displ, const __global float * d_xix, const __global float * d_xiy, const __global float * d_xiz, const __global float * d_etax, const __global float * d_etay, const __global float * d_etaz, const __global float * d_gammax, const __global float * d_gammay, const __global float * d_gammaz, const __global float * d_hprime_xx){\n\
  int ispec;\n\
  int ijk_ispec;\n\
  int tx;\n\
  int iglob;\n\
  float eps_trace_over_3;\n\
  float b_eps_trace_over_3;\n\
  float prod[21];\n\
  int i;\n\
  float epsdev[5];\n\
  float b_epsdev[5];\n\
  __local float s_dummyx_loc[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_dummyy_loc[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_dummyz_loc[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float sh_hprime_xx[NGLL2 + 0 - (1) - (0) + 1];\n\
  ispec = get_group_id(0) + (get_group_id(1)) * (get_num_groups(0));\n\
  ijk_ispec = get_local_id(0) + (NGLL3) * (ispec);\n\
  tx = get_local_id(0);\n\
  if(tx < NGLL2){\n\
    sh_hprime_xx[tx - 0] = d_hprime_xx[tx - 0];\n\
  }\n\
  if(ispec < NSPEC){\n\
    iglob = d_ibool[ijk_ispec - 0] - (1);\n\
    s_dummyx_loc[tx - 0] = d_b_displ[0 - 0 + (iglob - (0)) * (3)];\n\
    s_dummyy_loc[tx - 0] = d_b_displ[1 - 0 + (iglob - (0)) * (3)];\n\
    s_dummyz_loc[tx - 0] = d_b_displ[2 - 0 + (iglob - (0)) * (3)];\n\
  }\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if(ispec < NSPEC){\n\
    epsdev[0 - 0] = epsilondev_xx[ijk_ispec - 0];\n\
    epsdev[1 - 0] = epsilondev_yy[ijk_ispec - 0];\n\
    epsdev[2 - 0] = epsilondev_xy[ijk_ispec - 0];\n\
    epsdev[3 - 0] = epsilondev_xz[ijk_ispec - 0];\n\
    epsdev[4 - 0] = epsilondev_yz[ijk_ispec - 0];\n\
    eps_trace_over_3 = epsilon_trace_over_3[ijk_ispec - 0];\n\
    compute_element_strain_undo_att(ispec, ijk_ispec, d_ibool, s_dummyx_loc, s_dummyy_loc, s_dummyz_loc, d_xix, d_xiy, d_xiz, d_etax, d_etay, d_etaz, d_gammax, d_gammay, d_gammaz, sh_hprime_xx, b_epsdev,  &b_eps_trace_over_3);\n\
    compute_strain_product(prod, eps_trace_over_3, epsdev, b_eps_trace_over_3, b_epsdev);\n\
    for(i=0; i<=20; i+=1){\n\
      cijkl_kl[i - 0 + (ijk_ispec - (0)) * (21)] = cijkl_kl[i - 0 + (ijk_ispec - (0)) * (21)] + (deltat) * (prod[i - 0]);\n\
    }\n\
  }\n\
}\n\
";
