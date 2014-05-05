// from compute_forces_outer_core_cuda.cu
#define NGLLX 5
#define NGLL2 25
#define NGLL3 125
#define NGLL3_PADDED 128
#define R_EARTH_KM 6371.0f
#define COLORING_MIN_NSPEC_OUTER_CORE 1000

typedef float realw;
typedef float * realw_p;
typedef const float* __restrict__ realw_const_p;

#ifdef USE_TEXTURES_FIELDS
texture<realw, cudaTextureType1D, cudaReadModeElementType> d_displ_oc_tex;
texture<realw, cudaTextureType1D, cudaReadModeElementType> d_accel_oc_tex;
texture<realw, cudaTextureType1D, cudaReadModeElementType> d_b_displ_oc_tex;
texture<realw, cudaTextureType1D, cudaReadModeElementType> d_b_accel_oc_tex;
// templates definitions
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_displ_oc(int x);
template<int FORWARD_OR_ADJOINT> __device__ float texfetch_accel_oc(int x);
// templates for texture fetching
// FORWARD_OR_ADJOINT == 1 <- forward arrays
template<> __device__ float texfetch_displ_oc<1>(int x) { return tex1Dfetch(d_displ_oc_tex, x); }
template<> __device__ float texfetch_accel_oc<1>(int x) { return tex1Dfetch(d_accel_oc_tex, x); }
// FORWARD_OR_ADJOINT == 3 <- backward/reconstructed arrays
template<> __device__ float texfetch_displ_oc<3>(int x) { return tex1Dfetch(d_b_displ_oc_tex, x); }
template<> __device__ float texfetch_accel_oc<3>(int x) { return tex1Dfetch(d_b_accel_oc_tex, x); }
#endif

#ifdef USE_TEXTURES_CONSTANTS
texture<realw, cudaTextureType1D, cudaReadModeElementType> d_hprime_xx_oc_tex;
#endif

__device__ void compute_element_oc_rotation(int tx,int working_element,
                                            realw time,
                                            realw two_omega_earth,
                                            realw deltat,
                                            realw_p d_A_array_rotation,
                                            realw_p d_B_array_rotation,
                                            realw dpotentialdxl,
                                            realw dpotentialdyl,
                                            realw* dpotentialdx_with_rot,
                                            realw* dpotentialdy_with_rot) {

  realw two_omega_deltat,cos_two_omega_t,sin_two_omega_t;
  realw A_rotation,B_rotation;
  realw source_euler_A,source_euler_B;

  // store the source for the Euler scheme for A_rotation and B_rotation
  if( sizeof( sin_two_omega_t ) == sizeof( float ) ){
    // float operations
    // sincos function return sinus and cosine for given value
    sincosf(two_omega_earth*time, &sin_two_omega_t, &cos_two_omega_t);
  }else{
    cos_two_omega_t = cos(two_omega_earth*time);
    sin_two_omega_t = sin(two_omega_earth*time);
  }

  // time step deltat of Euler scheme is included in the source
  two_omega_deltat = deltat * two_omega_earth;

  source_euler_A = two_omega_deltat * (cos_two_omega_t * dpotentialdyl + sin_two_omega_t * dpotentialdxl);
  source_euler_B = two_omega_deltat * (sin_two_omega_t * dpotentialdyl - cos_two_omega_t * dpotentialdxl);

  A_rotation = d_A_array_rotation[tx + working_element*NGLL3]; // (non-padded offset)
  B_rotation = d_B_array_rotation[tx + working_element*NGLL3];

  *dpotentialdx_with_rot = dpotentialdxl + (  A_rotation*cos_two_omega_t + B_rotation*sin_two_omega_t);
  *dpotentialdy_with_rot = dpotentialdyl + (- A_rotation*sin_two_omega_t + B_rotation*cos_two_omega_t);

  // updates rotation term with Euler scheme (non-padded offset)
  d_A_array_rotation[tx + working_element*NGLL3] += source_euler_A;
  d_B_array_rotation[tx + working_element*NGLL3] += source_euler_B;
}



template<int FORWARD_OR_ADJOINT> __global__ void outer_core_impl_kernel(int nb_blocks_to_compute,
                                                                        const int* d_ibool,
                                                                        const int* d_phase_ispec_inner,
                                                                        const int num_phase_ispec,
                                                                        const int d_iphase,
                                                                        const int use_mesh_coloring_gpu,
                                                                        realw_const_p d_potential,
                                                                        realw_p d_potential_dot_dot,
                                                                        realw_const_p d_xix, realw_const_p d_xiy, realw_const_p d_xiz,
                                                                        realw_const_p d_etax, realw_const_p d_etay, realw_const_p d_etaz,
                                                                        realw_const_p d_gammax, realw_const_p d_gammay, realw_const_p d_gammaz,
                                                                        realw_const_p d_hprime_xx,
                                                                        realw_const_p d_hprimewgll_xx,
                                                                        realw_const_p wgllwgll_xy,
                                                                        realw_const_p wgllwgll_xz,
                                                                        realw_const_p wgllwgll_yz,
                                                                        const int GRAVITY,
                                                                        realw_const_p d_xstore, realw_const_p d_ystore, realw_const_p d_zstore,
                                                                        realw_const_p d_d_ln_density_dr_table,
                                                                        realw_const_p d_minus_rho_g_over_kappa_fluid,
                                                                        realw_const_p wgll_cube,
                                                                        const int ROTATION,
                                                                        realw time,
                                                                        realw two_omega_earth,
                                                                        realw deltat,
                                                                        realw_p d_A_array_rotation,
                                                                        realw_p d_B_array_rotation,
                                                                        const int NSPEC_OUTER_CORE){

  // block id == spectral-element id
  int bx = blockIdx.y*gridDim.x+blockIdx.x;
  // thread id == GLL point id
  int tx = threadIdx.x;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  unsigned short int active;
  int iglob,offset;
  int working_element;

  realw temp1l,temp2l,temp3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;

  realw dpotentialdxl,dpotentialdyl,dpotentialdzl;
  realw dpotentialdx_with_rot,dpotentialdy_with_rot;

  realw sum_terms;
  realw gravity_term;
  realw gxl,gyl,gzl;

  realw radius,theta,phi;
  realw cos_theta,sin_theta,cos_phi,sin_phi;
  realw grad_x_ln_rho,grad_y_ln_rho,grad_z_ln_rho;
  int int_radius;
#ifndef MANUALLY_UNROLLED_LOOPS
  int l;
#endif

  __shared__ realw s_dummy_loc[NGLL3];

  __shared__ realw s_temp1[NGLL3];
  __shared__ realw s_temp2[NGLL3];
  __shared__ realw s_temp3[NGLL3];

  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

// use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
// because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses
  active = (tx < NGLL3 && bx < nb_blocks_to_compute) ? 1:0;

// copy from global memory to shared memory
// each thread writes one of the NGLL^3 = 125 data points
  if (active) {

#ifdef USE_MESH_COLORING_GPU
    working_element = bx;
#else
    //mesh coloring
    if( use_mesh_coloring_gpu ){
      working_element = bx;
    }else{
      // iphase-1 and working_element-1 for Fortran->C array conventions
      working_element = d_phase_ispec_inner[bx + num_phase_ispec*(d_iphase-1)]-1;
    }
#endif

    // iglob = d_ibool[working_element*NGLL3_PADDED + tx]-1;
    iglob = d_ibool[working_element*NGLL3 + tx]-1;

#ifdef USE_TEXTURES_FIELDS
    s_dummy_loc[tx] = texfetch_displ_oc<FORWARD_OR_ADJOINT>(iglob);
#else
    // changing iglob indexing to match fortran row changes fast style
    s_dummy_loc[tx] = d_potential[iglob];
#endif
  }

  if (tx < NGLL2) {
    // hprime
    sh_hprime_xx[tx] = d_hprime_xx[tx];
    // weighted hprime
    sh_hprimewgll_xx[tx] = d_hprimewgll_xx[tx];
  }

// synchronize all the threads (one thread for each of the NGLL grid points of the
// current spectral element) because we need the whole element to be ready in order
// to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  if (active) {

#ifndef MANUALLY_UNROLLED_LOOPS
    temp1l = 0.f;
    temp2l = 0.f;
    temp3l = 0.f;
    for (l=0;l<NGLLX;l++) {
      temp1l += s_dummy_loc[K*NGLL2+J*NGLLX+l]*sh_hprime_xx[l*NGLLX+I];
      //assumes that hprime_xx = hprime_yy = hprime_zz
      temp2l += s_dummy_loc[K*NGLL2+l*NGLLX+I]*sh_hprime_xx[l*NGLLX+J];
      temp3l += s_dummy_loc[l*NGLL2+J*NGLLX+I]*sh_hprime_xx[l*NGLLX+K];
    }
#else
    temp1l = s_dummy_loc[K*NGLL2+J*NGLLX]*sh_hprime_xx[I]
            + s_dummy_loc[K*NGLL2+J*NGLLX+1]*sh_hprime_xx[NGLLX+I]
            + s_dummy_loc[K*NGLL2+J*NGLLX+2]*sh_hprime_xx[2*NGLLX+I]
            + s_dummy_loc[K*NGLL2+J*NGLLX+3]*sh_hprime_xx[3*NGLLX+I]
            + s_dummy_loc[K*NGLL2+J*NGLLX+4]*sh_hprime_xx[4*NGLLX+I];

    temp2l = s_dummy_loc[K*NGLL2+I]*sh_hprime_xx[J]
            + s_dummy_loc[K*NGLL2+NGLLX+I]*sh_hprime_xx[NGLLX+J]
            + s_dummy_loc[K*NGLL2+2*NGLLX+I]*sh_hprime_xx[2*NGLLX+J]
            + s_dummy_loc[K*NGLL2+3*NGLLX+I]*sh_hprime_xx[3*NGLLX+J]
            + s_dummy_loc[K*NGLL2+4*NGLLX+I]*sh_hprime_xx[4*NGLLX+J];

    temp3l = s_dummy_loc[J*NGLLX+I]*sh_hprime_xx[K]
            + s_dummy_loc[NGLL2+J*NGLLX+I]*sh_hprime_xx[NGLLX+K]
            + s_dummy_loc[2*NGLL2+J*NGLLX+I]*sh_hprime_xx[2*NGLLX+K]
            + s_dummy_loc[3*NGLL2+J*NGLLX+I]*sh_hprime_xx[3*NGLLX+K]
            + s_dummy_loc[4*NGLL2+J*NGLLX+I]*sh_hprime_xx[4*NGLLX+K];
#endif

    // compute derivatives of ux, uy and uz with respect to x, y and z
    offset = working_element*NGLL3_PADDED + tx;

    xixl = d_xix[offset];
    xiyl = d_xiy[offset];
    xizl = d_xiz[offset];
    etaxl = d_etax[offset];
    etayl = d_etay[offset];
    etazl = d_etaz[offset];
    gammaxl = d_gammax[offset];
    gammayl = d_gammay[offset];
    gammazl = d_gammaz[offset];

    //  compute the jacobian
    jacobianl = 1.f / (xixl*(etayl*gammazl-etazl*gammayl)
                    - xiyl*(etaxl*gammazl-etazl*gammaxl)
                    + xizl*(etaxl*gammayl-etayl*gammaxl));

    // derivatives of potential
    dpotentialdxl = xixl*temp1l + etaxl*temp2l + gammaxl*temp3l;
    dpotentialdyl = xiyl*temp1l + etayl*temp2l + gammayl*temp3l;
    dpotentialdzl = xizl*temp1l + etazl*temp2l + gammazl*temp3l;

    // compute contribution of rotation and add to gradient of potential
    // this term has no Z component
    if( ROTATION ){
      compute_element_oc_rotation(tx,working_element,
                                  time,two_omega_earth,deltat,
                                  d_A_array_rotation,d_B_array_rotation,
                                  dpotentialdxl,dpotentialdyl,
                                  &dpotentialdx_with_rot,&dpotentialdy_with_rot);
    }else{
      dpotentialdx_with_rot = dpotentialdxl;
      dpotentialdy_with_rot = dpotentialdyl;
    }

    // pre-computes gravity terms

    // use mesh coordinates to get theta and phi
    // x y z contain r theta phi
    radius = d_xstore[iglob];
    theta = d_ystore[iglob];
    phi = d_zstore[iglob];

    if( sizeof( theta ) == sizeof( float ) ){
      // float operations
      // sincos function return sinus and cosine for given value
      sincosf(theta, &sin_theta, &cos_theta);
      sincosf(phi, &sin_phi, &cos_phi);
    }else{
      cos_theta = cos(theta);
      sin_theta = sin(theta);
      cos_phi = cos(phi);
      sin_phi = sin(phi);
    }

    // for efficiency replace with lookup table every 100 m in radial direction
    // note: radius in outer core should never be zero,
    //          and arrays in C start from 0, thus we need to subtract -1
    int_radius = rint(radius * R_EARTH_KM * 10.0f ) - 1;

    //debug: checks bounds NRAD_GRAVITY == 70000
    //if( int_radius < 0 || int_radius >= 70000 ){
    //  printf("gravity: in_radius out of bounds %d radius=%e\n",int_radius,radius);
    //}

    // depending on gravity or not, different potential definitions are used
    if( ! GRAVITY ){
      // add (chi/rho)grad(rho) term in no gravity case

      // grad(rho)/rho in Cartesian components
      grad_x_ln_rho = sin_theta * cos_phi * d_d_ln_density_dr_table[int_radius];
      grad_y_ln_rho = sin_theta * sin_phi * d_d_ln_density_dr_table[int_radius];
      grad_z_ln_rho = cos_theta * d_d_ln_density_dr_table[int_radius];

      // adding (chi/rho)grad(rho)
      dpotentialdx_with_rot += s_dummy_loc[tx] * grad_x_ln_rho;
      dpotentialdy_with_rot += s_dummy_loc[tx] * grad_y_ln_rho;
      dpotentialdzl += s_dummy_loc[tx] * grad_z_ln_rho;

    }else{

      // compute divergence of displacement
      // precompute and store gravity term
      //
      // get g, rho and dg/dr=dg
      // spherical components of the gravitational acceleration
      //
      // Cartesian components of the gravitational acceleration
      // integrate and multiply by rho / Kappa
      gxl = sin_theta*cos_phi;
      gyl = sin_theta*sin_phi;
      gzl = cos_theta;

      // uses potential definition: s = grad(chi)
      // gravity term: - rho * g * 1/kappa grad(chi)
      gravity_term = d_minus_rho_g_over_kappa_fluid[int_radius] * jacobianl * wgll_cube[tx] *
                    ( dpotentialdx_with_rot * gxl + dpotentialdy_with_rot * gyl + dpotentialdzl * gzl);

      // divergence of displacement field with gravity on
      // note: these calculations are only considered for SIMULATION_TYPE == 1 .and. SAVE_FORWARD
      //          and one has set MOVIE_VOLUME_TYPE == 4 when MOVIE_VOLUME is .true.;
      //         in case of SIMULATION_TYPE == 3, it gets overwritten by compute_kernels_outer_core()
      //if (NSPEC_OUTER_CORE_ADJOINT /= 1 && MOVIE_VOLUME ){
      //  div_displfluid(i,j,k,ispec) =  d_minus_rho_g_over_kappa_fluid[int_radius] *
      //        (dpotentialdx_with_rot * gxl + dpotentialdy_with_rot * gyl + dpotentialdzl * gzl);
      //}
    }

    // form the dot product with the test vector
    s_temp1[tx] = jacobianl*(xixl*dpotentialdx_with_rot
                             + xiyl*dpotentialdy_with_rot + xizl*dpotentialdzl);
    s_temp2[tx] = jacobianl*(etaxl*dpotentialdx_with_rot
                             + etayl*dpotentialdy_with_rot + etazl*dpotentialdzl);
    s_temp3[tx] = jacobianl*(gammaxl*dpotentialdx_with_rot
                             + gammayl*dpotentialdy_with_rot + gammazl*dpotentialdzl);
  }

// synchronize all the threads (one thread for each of the NGLL grid points of the
// current spectral element) because we need the whole element to be ready in order
// to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  if (active) {

#ifndef MANUALLY_UNROLLED_LOOPS
    temp1l = 0.f;
    temp2l = 0.f;
    temp3l = 0.f;
    for (l=0;l<NGLLX;l++) {
        temp1l += s_temp1[K*NGLL2+J*NGLLX+l]*sh_hprimewgll_xx[I*NGLLX+l];
        //assumes hprimewgll_xx = hprimewgll_yy = hprimewgll_zz
        temp2l += s_temp2[K*NGLL2+l*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+l];
        temp3l += s_temp3[l*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+l];
    }
#else
    temp1l = s_temp1[K*NGLL2+J*NGLLX]*sh_hprimewgll_xx[I*NGLLX]
            + s_temp1[K*NGLL2+J*NGLLX+1]*sh_hprimewgll_xx[I*NGLLX+1]
            + s_temp1[K*NGLL2+J*NGLLX+2]*sh_hprimewgll_xx[I*NGLLX+2]
            + s_temp1[K*NGLL2+J*NGLLX+3]*sh_hprimewgll_xx[I*NGLLX+3]
            + s_temp1[K*NGLL2+J*NGLLX+4]*sh_hprimewgll_xx[I*NGLLX+4];

    temp2l = s_temp2[K*NGLL2+I]*sh_hprimewgll_xx[J*NGLLX]
            + s_temp2[K*NGLL2+NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+1]
            + s_temp2[K*NGLL2+2*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+2]
            + s_temp2[K*NGLL2+3*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+3]
            + s_temp2[K*NGLL2+4*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+4];

    temp3l = s_temp3[J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX]
            + s_temp3[NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+1]
            + s_temp3[2*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+2]
            + s_temp3[3*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+3]
            + s_temp3[4*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+4];
#endif

    sum_terms = - ( wgllwgll_yz[K*NGLLX+J]*temp1l
                  + wgllwgll_xz[K*NGLLX+I]*temp2l
                  + wgllwgll_xy[J*NGLLX+I]*temp3l);

    if( GRAVITY ) sum_terms += gravity_term;

    //iglob = d_ibool[working_element*NGLL3 + tx]-1;

#ifdef USE_MESH_COLORING_GPU
    // no atomic operation needed, colors don't share global points between elements

#ifdef USE_TEXTURES_FIELDS
    d_potential_dot_dot[iglob] = texfetch_accel_oc<FORWARD_OR_ADJOINT>(iglob) + sum_terms;
#else
    d_potential_dot_dot[iglob] += sum_terms;
#endif // USE_TEXTURES_FIELDS

#else // MESH_COLORING

    //mesh coloring
    if( use_mesh_coloring_gpu ){

      if( NSPEC_OUTER_CORE > COLORING_MIN_NSPEC_OUTER_CORE ){
        // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
    d_potential_dot_dot[iglob] = texfetch_accel_oc<FORWARD_OR_ADJOINT>(iglob) + sum_terms;
#else
        d_potential_dot_dot[iglob] += sum_terms;
#endif // USE_TEXTURES_FIELDS
      }else{
        // poor element count, only use 1 color per inner/outer run
        // forces atomic operations
        atomicAdd(&d_potential_dot_dot[iglob],sum_terms);
      }

    }else{

      atomicAdd(&d_potential_dot_dot[iglob],sum_terms);

    }
#endif // MESH_COLORING
  }
}
