// from compute_kernels_cuda.cu
#define NGLL3 125

typedef float realw;

__global__ void compute_iso_kernel(realw* epsilondev_xx,
                                   realw* epsilondev_yy,
                                   realw* epsilondev_xy,
                                   realw* epsilondev_xz,
                                   realw* epsilondev_yz,
                                   realw* epsilon_trace_over_3,
                                   realw* b_epsilondev_xx,
                                   realw* b_epsilondev_yy,
                                   realw* b_epsilondev_xy,
                                   realw* b_epsilondev_xz,
                                   realw* b_epsilondev_yz,
                                   realw* b_epsilon_trace_over_3,
                                   realw* mu_kl,
                                   realw* kappa_kl,
                                   int NSPEC,
                                   realw deltat) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;

  // handles case when there is 1 extra block (due to rectangular grid)
  if(ispec < NSPEC) {

    int ijk_ispec = threadIdx.x + NGLL3*ispec;

    // isotropic kernel contributions
    // shear modulus kernel
    mu_kl[ijk_ispec] += deltat * (epsilondev_xx[ijk_ispec]*b_epsilondev_xx[ijk_ispec]+
                                  epsilondev_yy[ijk_ispec]*b_epsilondev_yy[ijk_ispec]+
                                  (epsilondev_xx[ijk_ispec]+epsilondev_yy[ijk_ispec])*
                                    (b_epsilondev_xx[ijk_ispec]+b_epsilondev_yy[ijk_ispec])+
                                    2*(epsilondev_xy[ijk_ispec]*b_epsilondev_xy[ijk_ispec]+
                                       epsilondev_xz[ijk_ispec]*b_epsilondev_xz[ijk_ispec]+
                                       epsilondev_yz[ijk_ispec]*b_epsilondev_yz[ijk_ispec]));

    // bulk modulus kernel
    kappa_kl[ijk_ispec] += deltat * ( 9 * epsilon_trace_over_3[ijk_ispec] * b_epsilon_trace_over_3[ijk_ispec]);
  }
}
