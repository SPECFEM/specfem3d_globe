// from compute_kernels_cuda.cu
#define NGLLX 5
#define NGLL2 25
#define NGLL3 125
#define NGLL3_PADDED 128

typedef float realw;

__device__ void compute_gradient_kernel(int ijk,
                                        int ispec,
                                        realw* scalar_field,
                                        realw* vector_field_element,
                                        realw* hprime_xx,
                                        realw* d_xix,
                                        realw* d_xiy,
                                        realw* d_xiz,
                                        realw* d_etax,
                                        realw* d_etay,
                                        realw* d_etaz,
                                        realw* d_gammax,
                                        realw* d_gammay,
                                        realw* d_gammaz) {

  realw temp1l,temp2l,temp3l;
  realw hp1,hp2,hp3;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl;
  int l,offset,offset1,offset2,offset3;

  int K = (ijk/NGLL2);
  int J = ((ijk-K*NGLL2)/NGLLX);
  int I = (ijk-K*NGLL2-J*NGLLX);

  // derivative along x
  temp1l = 0.f;
  for( l=0; l<NGLLX;l++){
    hp1 = hprime_xx[l*NGLLX+I];
    offset1 = K*NGLL2+J*NGLLX+l;
    temp1l += scalar_field[offset1]*hp1;
  }

  // derivative along y
  temp2l = 0.f;
  for( l=0; l<NGLLX;l++){
    //assumes that hprime_xx = hprime_yy = hprime_zz
    hp2 = hprime_xx[l*NGLLX+J];
    offset2 = K*NGLL2+l*NGLLX+I;
    temp2l += scalar_field[offset2]*hp2;
  }

  // derivative along z
  temp3l = 0.f;
  for( l=0; l<NGLLX;l++){
    //assumes that hprime_xx = hprime_yy = hprime_zz
    hp3 = hprime_xx[l*NGLLX+K];
    offset3 = l*NGLL2+J*NGLLX+I;
    temp3l += scalar_field[offset3]*hp3;
  }

  offset = ispec*NGLL3_PADDED + ijk;

  xixl = d_xix[offset];
  xiyl = d_xiy[offset];
  xizl = d_xiz[offset];
  etaxl = d_etax[offset];
  etayl = d_etay[offset];
  etazl = d_etaz[offset];
  gammaxl = d_gammax[offset];
  gammayl = d_gammay[offset];
  gammazl = d_gammaz[offset];

  // note: global version uses a different potential definition, no need to divide by rho
  //rho_invl = 1.0f / rhol;

  // derivatives of acoustic scalar potential field on GLL points
  vector_field_element[0] = temp1l*xixl + temp2l*etaxl + temp3l*gammaxl;
  vector_field_element[1] = temp1l*xiyl + temp2l*etayl + temp3l*gammayl;
  vector_field_element[2] = temp1l*xizl + temp2l*etazl + temp3l*gammazl;

}

__global__ void compute_acoustic_kernel(int* ibool,
                                        realw* rhostore,
                                        realw* kappastore,
                                        realw* hprime_xx,
                                        realw* d_xix,
                                        realw* d_xiy,
                                        realw* d_xiz,
                                        realw* d_etax,
                                        realw* d_etay,
                                        realw* d_etaz,
                                        realw* d_gammax,
                                        realw* d_gammay,
                                        realw* d_gammaz,
                                        realw* potential_dot_dot_acoustic,
                                        realw* b_potential_acoustic,
                                        realw* b_potential_dot_dot_acoustic,
                                        realw* rho_ac_kl,
                                        realw* kappa_ac_kl,
                                        realw deltat,
                                        int NSPEC) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;

  // handles case when there is 1 extra block (due to rectangular grid)
  if( ispec < NSPEC ){

    int ijk = threadIdx.x;

    // local and global indices
    int ijk_ispec = ijk + NGLL3*ispec;
    int ijk_ispec_padded = ijk + NGLL3_PADDED*ispec;
    int iglob = ibool[ijk_ispec] - 1;

    realw accel_elm[3];
    realw b_displ_elm[3];
    realw rhol,kappal;
    realw div_displ,b_div_displ;

    // shared memory between all threads within this block
    __shared__ realw scalar_field_displ[NGLL3];
    __shared__ realw scalar_field_accel[NGLL3];

    // copy field values
    scalar_field_displ[ijk] = b_potential_acoustic[iglob];
    scalar_field_accel[ijk] = potential_dot_dot_acoustic[iglob];
    __syncthreads();

    // displacement vector from backward field
    compute_gradient_kernel(ijk,ispec,scalar_field_displ,b_displ_elm,
                            hprime_xx,
                            d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz);

    // acceleration vector
    compute_gradient_kernel(ijk,ispec,scalar_field_accel,accel_elm,
                            hprime_xx,
                            d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz);

    // gets material parameter
    rhol = rhostore[ijk_ispec_padded];

    // density kernel
    rho_ac_kl[ijk_ispec] += deltat * rhol * (accel_elm[0]*b_displ_elm[0] +
                                             accel_elm[1]*b_displ_elm[1] +
                                             accel_elm[2]*b_displ_elm[2]);

    // bulk modulus kernel
    kappal = rhol/ kappastore[ijk_ispec_padded];

    div_displ = kappal * potential_dot_dot_acoustic[iglob];
    b_div_displ = kappal * b_potential_dot_dot_acoustic[iglob];

    kappa_ac_kl[ijk_ispec] += deltat * div_displ * b_div_displ;
  }
}
