// from compute_kernels_cuda.cu
#define NGLL3 125

typedef float realw;

__device__ void compute_strain_product(realw* prod,
                                       realw eps_trace_over_3,
                                       realw* epsdev,
                                       realw b_eps_trace_over_3,
                                       realw* b_epsdev){

  realw eps[6],b_eps[6];

  // Building of the local matrix of the strain tensor
  // for the adjoint field and the regular backward field

  // note: indices are -1 compared to fortran routine because of fortran -> C array indexing

  // eps11 et eps22
  eps[0] = epsdev[0] + eps_trace_over_3;
  eps[1] = epsdev[1] + eps_trace_over_3;
  //eps33
  eps[2] = - (eps[0] + eps[1]) + 3.0f*eps_trace_over_3;
  //eps23
  eps[3] = epsdev[4];
  //eps13
  eps[4] = epsdev[3];
  //eps12
  eps[5] = epsdev[2];

  b_eps[0] = b_epsdev[0] + b_eps_trace_over_3;
  b_eps[1] = b_epsdev[1] + b_eps_trace_over_3;
  b_eps[2] = - (b_eps[0] + b_eps[1]) + 3.0f*b_eps_trace_over_3;
  b_eps[3] = b_epsdev[4];
  b_eps[4] = b_epsdev[3];
  b_eps[5] = b_epsdev[2];

  // Computing the 21 strain products without assuming eps(i)*b_eps(j) = eps(j)*b_eps(i)
  int p = 0;
  for(int i=0; i<6; i++){
    for(int j=i; j<6; j++){
      prod[p]=eps[i]*b_eps[j];
      if(j>i){
        prod[p]=prod[p]+eps[j]*b_eps[i];
        if(j>2 && i<3){ prod[p] = prod[p]*2.0f;}
      }
      if(i>2){ prod[p]=prod[p]*4.0f;}
      p=p+1;
    }
  }
}

__global__ void compute_ani_kernel(realw* epsilondev_xx,
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
                                   realw* cijkl_kl,
                                   int NSPEC,
                                   realw deltat) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;

  // handles case when there is 1 extra block (due to rectangular grid)
  if(ispec < NSPEC) {

    int ijk_ispec = threadIdx.x + NGLL3*ispec;

    // fully anisotropic kernel contributions
    realw eps_trace_over_3,b_eps_trace_over_3;
    realw prod[21];
    realw epsdev[5];
    realw b_epsdev[5];

    epsdev[0] = epsilondev_xx[ijk_ispec];
    epsdev[1] = epsilondev_yy[ijk_ispec];
    epsdev[2] = epsilondev_xy[ijk_ispec];
    epsdev[3] = epsilondev_xz[ijk_ispec];
    epsdev[4] = epsilondev_yz[ijk_ispec];

    b_epsdev[0] = b_epsilondev_xx[ijk_ispec];
    b_epsdev[1] = b_epsilondev_yy[ijk_ispec];
    b_epsdev[2] = b_epsilondev_xy[ijk_ispec];
    b_epsdev[3] = b_epsilondev_xz[ijk_ispec];
    b_epsdev[4] = b_epsilondev_yz[ijk_ispec];

    eps_trace_over_3 = epsilon_trace_over_3[ijk_ispec];
    b_eps_trace_over_3 = b_epsilon_trace_over_3[ijk_ispec];

    // fully anisotropic kernel contributions
    compute_strain_product(prod,eps_trace_over_3,epsdev,b_eps_trace_over_3,b_epsdev);

    // updates full anisotropic kernel
    for(int i=0;i<21;i++){
      cijkl_kl[i + 21*ijk_ispec] += deltat * prod[i];
    }
  }
}
