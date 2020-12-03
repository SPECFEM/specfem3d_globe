// from compute_kernels_cuda.cu
#define NGLL3 125

typedef float realw;

__device__ void compute_gradient(
    const int ispec,
    const realw* fx,
    const realw* fy,
    const realw* fz,
    const realw* d_xix, const realw* d_xiy, const realw* dxiz,
    const realw* d_etax, const realw* d_etay, const realw* d_etaz,
    const realw* d_gammax, const realw* d_gammay, const realw* d_gammaz,
    const realw* sh_hprime_xx,
    realw* f_grad)
{
  int tx = threadIdx.x;
  int K = tx / NGLL2;
  int J = (tx - K*NGLL2) / NGLLX;
  int I = (tx - K*NLGG2 - J*NGLLX);

  int l;
  int offset;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl;
  realw dfxdxl, dfxdyl, dfxdzl, dfydxl, dfydyl, dfydzl, dfzdxl, dfzdyl, dfzdzl;
  realw fac1, fac2, fac3;

  tempx1l = 0.f;
  tempx2l = 0.f;
  tempx3l = 0.f;

  tempy1l = 0.f;
  tempy2l = 0.f;
  tempy3l = 0.f;

  tempz1l = 0.f;
  tempz2l = 0.f;
  tempz3l = 0.f;

  for(l=0; l<NGLLX; ++l) {
    fac1 = sh_hprime_xx[l*NGLLX+I];
    tempx1l += fx[K*NGLL2 + J*NGLLX + l] * fac1;
    tempy1l += fy[K*NGLL2 + J*NGLLX + l] * fac1;
    tempz1l += fz[K*NGLL2 + J*NGLLX + l] * fac1;

    fac2 = sh_hprime_xx[l*NGLLX+J];
    tempx2l += fx[K*NGLL2 + l*NGLLX + I] * fac2;
    tempy2l += fy[K*NGLL2 + l*NGLLX + I] * fac2;
    tempz2l += fz[K*NGLL2 + l*NGLLX + I] * fac2;

    fac2 = sh_hprime_xx[l*NGLLX+K];
    tempx3l += fx[K*NGLL2 + l*NGLLX + I] * fac2;
    tempy3l += fy[K*NGLL2 + l*NGLLX + I] * fac2;
    tempz3l += fz[K*NGLL2 + l*NGLLX + I] * fac2;
  }

  offset = ispec * NGLL3_PADDED + tx;

  xixl = d_xix[offset];
  xiyl = d_xiy[offset];
  xizl = d_xiz[offset];
  etaxl = d_etax[offset];
  etayl = d_etay[offset];
  etazl = d_etaz[offset];
  gammaxl = d_gammax[offset];
  gammayl = d_gammay[offset];
  gammazl = d_gammaz[offset];

  dfxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l;
  dfxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l;
  dfxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l;

  dfydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l;
  dfydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l;
  dfydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l;

  dfzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l;
  dfzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l;
  dfzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l;

  f_grad[0] = dfxdxl; f_grad[1] = dfxdyl; f_grad[2] = dfxdzl;
  f_grad[3] = dfydxl; f_grad[4] = dfydyl; f_grad[5] = dfydzl;
  f_grad[6] = dfzdxl; f_grad[7] = dfzdyl; f_grad[8] = dfzdzl;
}


__global__ void compute_kappa_mu_hess_kernel(int *ibool,
                                             realw* veloc
                                             realw* b_veloc,
                                             realw* d_xis, realw* d_xiy, realw* dxiz,
                                             realw* d_etax,realw* d_etay,realw* d_etaz,
                                             realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                                             realw* d_hprime_xx,
                                             realw deltat,
                                             realw* hess_rho_kl,
                                             realw* hess_kappa_kl,
                                             realw* hess_mu_kl,
                                             int NSPEC_AB,
                                             int USE_SOURCE_RECEIVER_HESSIAN) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;
  int ijk = tx;
  int ijk_ispec = ijk + NGLL3 * ispec;
  int iglob;

  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_velocx[NGLL3];
  __shared__ realw sh_velocy[NGLL3];
  __shared__ realw sh_velocz[NGLL3];
  __shared__ realw sh_b_velocx[NGLL3];
  __shared__ realw sh_b_velocy[NGLL3];
  __shared__ realw sh_b_velocz[NGLL3];

  realw vgrad[9], b_vgrad[9];
  realw hess_rhol, hess_kappal, hess_mul;

  // load hprime into shared memory
  if(tx < NGLL2) {
    sh_hprime_xx[tx] = d_hprime_xx[tx];
  }

  if(ispec < NSPEC_AB) {
    iglob = d_ibool[ijk_ispec] - 1;

    // load velocity (forward) into shared memory
    sh_velocx[ijk] = veloc[iglob*3];
    sh_velocy[ijk] = veloc[iglob*3+1];
    sh_velocz[ijk] = veloc[iglob*3+2];

    // load velocity (adjoint) into shared memory
    sh_b_velocx[ijk] = b_veloc[iglob*3];
    sh_b_velocy[ijk] = b_veloc[iglob*3+1];
    sh_b_velocz[ijk] = b_veloc[iglob*3+2];

  }

  // synchronizes threads
  __syncthreads();

  if(ispec < NSPEC_AB) {

    // compute the gradient of velocity (forward) field
    compute_gradient(
        ispec,
        sh_velocx, sh_velocy, sh_velocz,
        d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,
        sh_hprime_xx,
        vgrad);

    if (USE_SOURCE_RECEIVER_HESSIAN) {
      // compute the gradient of velocity (forward) field
      compute_gradient(ispec,
          sh_b_velocx, sh_b_velocy, sh_b_velocz,
          d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,
          sh_hprime_xx,
          b_vgrad);
    } else {
      // if using source-source hessian, then just copy the forward velocity gradient
      for(int i=0; i<9; ++i) {
        b_vgrad[i] = vgrad[i];
      }
    }

    // velocity gradient is in such order | row-major |:
    // | dVx/dx, dVx/dy, dVx/dz |         |  0, 1, 2  |
    // | dVy/dx, dVy/dy, dVy/dz |         |  3, 4, 5  |
    // | dVz/dx, dVz/dy, dVz/dz |         |  6, 7, 8  |
    hess_rhol = vgrad[0]*b_vgrad[0] + vgrad[4]*b_vgrad[4] + vgrad[8]*b_vgrad[8];
    hess_kappal = (vgrad[0] + vgrad[4] + vgrad[8]) * (b_vgrad[0] + b_vgrad[4] + b_vgrad[8]);
    hess_mul = (vgrad[1] + vgrad[3]) * (b_vgrad[1] + b_vgrad[3])
             + (vgrad[2] + vgrad[6]) * (b_vgrad[2] + b_vgrad[6])
             + (vgrad[5] + vgrad[7]) * (b_vgrad[5] + b_vgrad[7]);

    hess_rho_kl[ijk_ispec] += deltat * hess_rhol;
    hess_kappa_kl[ijk_ispec] += deltat * hess_kappal;
    hess_mu_kl[ijk_ispec] += deltat * hess_mul;
  }
}
