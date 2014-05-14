// from compute_kernels_cuda.cu
#define NGLLX 5
#define NGLL2 25
#define NGLL3 125
#define NGLL3_PADDED 128

typedef float realw;

__device__ void compute_element_strain_undo_att(int ispec,int ijk_ispec,
                                                int* d_ibool,
                                                realw* s_dummyx_loc,
                                                realw* s_dummyy_loc,
                                                realw* s_dummyz_loc,
                                                realw* d_xix,realw* d_xiy,realw* d_xiz,
                                                realw* d_etax,realw* d_etay,realw* d_etaz,
                                                realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                                                realw* sh_hprime_xx,
                                                realw* epsilondev_loc,
                                                realw* epsilon_trace_over_3) {


  // thread id == GLL point id
  int tx = threadIdx.x;
  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int offset;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw templ;
  realw fac1,fac2,fac3;

  int l;

// copy from global memory to shared memory
// each thread writes one of the NGLL^3 = 125 data points

  tempx1l = 0.f;
  tempx2l = 0.f;
  tempx3l = 0.f;

  tempy1l = 0.f;
  tempy2l = 0.f;
  tempy3l = 0.f;

  tempz1l = 0.f;
  tempz2l = 0.f;
  tempz3l = 0.f;

  for (l=0;l<NGLLX;l++) {
      fac1 = sh_hprime_xx[l*NGLLX+I];
      tempx1l += s_dummyx_loc[K*NGLL2+J*NGLLX+l]*fac1;
      tempy1l += s_dummyy_loc[K*NGLL2+J*NGLLX+l]*fac1;
      tempz1l += s_dummyz_loc[K*NGLL2+J*NGLLX+l]*fac1;

      fac2 = sh_hprime_xx[l*NGLLX+J];
      tempx2l += s_dummyx_loc[K*NGLL2+l*NGLLX+I]*fac2;
      tempy2l += s_dummyy_loc[K*NGLL2+l*NGLLX+I]*fac2;
      tempz2l += s_dummyz_loc[K*NGLL2+l*NGLLX+I]*fac2;

      fac3 = sh_hprime_xx[l*NGLLX+K];
      tempx3l += s_dummyx_loc[l*NGLL2+J*NGLLX+I]*fac3;
      tempy3l += s_dummyy_loc[l*NGLL2+J*NGLLX+I]*fac3;
      tempz3l += s_dummyz_loc[l*NGLL2+J*NGLLX+I]*fac3;
  }

  // compute derivatives of ux, uy and uz with respect to x, y and z
  offset = ispec*NGLL3_PADDED + tx;

  xixl = d_xix[offset];
  xiyl = d_xiy[offset];
  xizl = d_xiz[offset];
  etaxl = d_etax[offset];
  etayl = d_etay[offset];
  etazl = d_etaz[offset];
  gammaxl = d_gammax[offset];
  gammayl = d_gammay[offset];
  gammazl = d_gammaz[offset];

  duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l;
  duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l;
  duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l;

  duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l;
  duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l;
  duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l;

  duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l;
  duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l;
  duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l;

  // computes deviatoric strain attenuation and/or for kernel calculations
  templ = 0.33333333333333333333f * (duxdxl + duydyl + duzdzl); // 1./3. = 0.33333

  // local storage: stresses at this current time step
  epsilondev_loc[0] = duxdxl - templ;   // xx
  epsilondev_loc[1] = duydyl - templ;   // yy
  epsilondev_loc[2] = 0.5f * ( duxdyl + duydxl ); // xy
  epsilondev_loc[3] = 0.5f * ( duzdxl + duxdzl ); // xz
  epsilondev_loc[4] = 0.5f * ( duzdyl + duydzl ); // yz
  *epsilon_trace_over_3 = templ;
}

__global__ void compute_iso_undo_att_kernel(realw* epsilondev_xx,
                                            realw* epsilondev_yy,
                                            realw* epsilondev_xy,
                                            realw* epsilondev_xz,
                                            realw* epsilondev_yz,
                                            realw* epsilon_trace_over_3,
                                            realw* mu_kl,
                                            realw* kappa_kl,
                                            int NSPEC,
                                            realw deltat,
                                            int* d_ibool,
                                            realw* d_b_displ,
                                            realw* d_xix,
                                            realw* d_xiy,
                                            realw* d_xiz,
                                            realw* d_etax,
                                            realw* d_etay,
                                            realw* d_etaz,
                                            realw* d_gammax,
                                            realw* d_gammay,
                                            realw* d_gammaz,
                                            realw* d_hprime_xx) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ijk_ispec = threadIdx.x + NGLL3*ispec;

  int tx = threadIdx.x;
  int iglob;

  __shared__ realw s_dummyx_loc[NGLL3];
  __shared__ realw s_dummyy_loc[NGLL3];
  __shared__ realw s_dummyz_loc[NGLL3];

  __shared__ realw sh_hprime_xx[NGLL2];

  // loads element displacements
  // all threads load their displacement into shared memory
  if( ispec < NSPEC){
    iglob = d_ibool[ijk_ispec]-1;
    // changing iglob indexing to match fortran row changes fast style
    s_dummyx_loc[tx] = d_b_displ[iglob*3];
    s_dummyy_loc[tx] = d_b_displ[iglob*3 + 1];
    s_dummyz_loc[tx] = d_b_displ[iglob*3 + 2];

    // master thread loads hprime
    if( threadIdx.x == 0 ){
      for(int m=0; m < NGLL2; m++){
        // hprime
        sh_hprime_xx[m] = d_hprime_xx[m];
      }
    }
  }

  // synchronizes threads
  __syncthreads();

  // handles case when there is 1 extra block (due to rectangular grid)
  if(ispec < NSPEC) {

    realw eps_trace_over_3,b_eps_trace_over_3;
    realw epsdev[5];
    realw b_epsdev[5];

    // strain from adjoint wavefield
    epsdev[0] = epsilondev_xx[ijk_ispec];
    epsdev[1] = epsilondev_yy[ijk_ispec];
    epsdev[2] = epsilondev_xy[ijk_ispec];
    epsdev[3] = epsilondev_xz[ijk_ispec];
    epsdev[4] = epsilondev_yz[ijk_ispec];
    eps_trace_over_3 = epsilon_trace_over_3[ijk_ispec];

    // strain from backward/reconstructed forward wavefield
    compute_element_strain_undo_att(ispec,ijk_ispec,
                                    d_ibool,
                                    s_dummyx_loc,s_dummyy_loc,s_dummyz_loc,
                                    d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,
                                    sh_hprime_xx,
                                    b_epsdev,&b_eps_trace_over_3);

    // isotropic kernel contributions
    // shear modulus kernel
    mu_kl[ijk_ispec] += deltat * ( epsdev[0]*b_epsdev[0] + epsdev[1]*b_epsdev[1]
                                   + (epsdev[0]+epsdev[1])*(b_epsdev[0]+b_epsdev[1])
                                   + 2*( epsdev[2]*b_epsdev[2] + epsdev[3]*b_epsdev[3] + epsdev[4]*b_epsdev[4]) );

    // bulk modulus kernel
    kappa_kl[ijk_ispec] += deltat * ( 9 * eps_trace_over_3 * b_eps_trace_over_3);
  }
}
