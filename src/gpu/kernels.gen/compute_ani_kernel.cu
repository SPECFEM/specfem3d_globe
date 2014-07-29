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
__global__ void compute_ani_kernel(const float * epsilondev_xx, const float * epsilondev_yy, const float * epsilondev_xy, const float * epsilondev_xz, const float * epsilondev_yz, const float * epsilon_trace_over_3, const float * b_epsilondev_xx, const float * b_epsilondev_yy, const float * b_epsilondev_xy, const float * b_epsilondev_xz, const float * b_epsilondev_yz, const float * b_epsilon_trace_over_3, float * cijkl_kl, const int NSPEC, const float deltat){
  int i;
  int ispec;
  int ijk_ispec;
  float eps_trace_over_3;
  float b_eps_trace_over_3;
  float prod[21];
  float epsdev[5];
  float b_epsdev[5];
  ispec = blockIdx.x + (blockIdx.y) * (gridDim.x);
  if(ispec < NSPEC){
    ijk_ispec = threadIdx.x + (NGLL3) * (ispec);
    epsdev[0 - (0)] = epsilondev_xx[ijk_ispec - (0)];
    epsdev[1 - (0)] = epsilondev_yy[ijk_ispec - (0)];
    epsdev[2 - (0)] = epsilondev_xy[ijk_ispec - (0)];
    epsdev[3 - (0)] = epsilondev_xz[ijk_ispec - (0)];
    epsdev[4 - (0)] = epsilondev_yz[ijk_ispec - (0)];
    epsdev[0 - (0)] = b_epsilondev_xx[ijk_ispec - (0)];
    epsdev[1 - (0)] = b_epsilondev_yy[ijk_ispec - (0)];
    epsdev[2 - (0)] = b_epsilondev_xy[ijk_ispec - (0)];
    epsdev[3 - (0)] = b_epsilondev_xz[ijk_ispec - (0)];
    epsdev[4 - (0)] = b_epsilondev_yz[ijk_ispec - (0)];
    eps_trace_over_3 = epsilon_trace_over_3[ijk_ispec - (0)];
    b_eps_trace_over_3 = b_epsilon_trace_over_3[ijk_ispec - (0)];
    compute_strain_product(prod, eps_trace_over_3, epsdev, b_eps_trace_over_3, b_epsdev);
    for(i=0; i<=20; i+=1){
      cijkl_kl[i - (0) + (ijk_ispec - (0)) * (21)] = cijkl_kl[i - (0) + (ijk_ispec - (0)) * (21)] + (deltat) * (prod[i - (0)]);
    }
  }
}
