/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  2 . 0
 !               ---------------------------------------
 !
 !          Main authors: Dimitri Komatitsch and Jeroen Tromp
 !    Princeton University, USA and University of Pau / CNRS / INRIA
 ! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
 !                            August 2013
 !
 ! This program is free software; you can redistribute it and/or modify
 ! it under the terms of the GNU General Public License as published by
 ! the Free Software Foundation; either version 2 of the License, or
 ! (at your option) any later version.
 !
 ! This program is distributed in the hope that it will be useful,
 ! but WITHOUT ANY WARRANTY; without even the implied warranty of
 ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ! GNU General Public License for more details.
 !
 ! You should have received a copy of the GNU General Public License along
 ! with this program; if not, write to the Free Software Foundation, Inc.,
 ! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 !
 !=====================================================================
 */

#include <stdio.h>

#include <cuda.h>
#include <cublas.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"

#ifdef USE_TEXTURES_FIELDS
texture<realw, cudaTextureType1D, cudaReadModeElementType> d_displ_oc_tex;
texture<realw, cudaTextureType1D, cudaReadModeElementType> d_accel_oc_tex;
#endif

#ifdef USE_TEXTURES_CONSTANTS
texture<realw, cudaTextureType1D, cudaReadModeElementType> d_hprime_xx_oc_tex;
#endif

/* ----------------------------------------------------------------------------------------------- */

// elemental routines

/* ----------------------------------------------------------------------------------------------- */

// fluid rotation

__device__ void compute_element_oc_rotation(int tx,int working_element,
                                            realw time,
                                            realw two_omega_earth,
                                            realw deltat,
                                            realw* d_A_array_rotation,
                                            realw* d_B_array_rotation,
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


/* ----------------------------------------------------------------------------------------------- */

// KERNEL 2
//
// for outer core ( acoustic domain )
/* ----------------------------------------------------------------------------------------------- */


__global__ void Kernel_2_outer_core_impl(int nb_blocks_to_compute,
                                         int NGLOB,
                                         int* d_ibool,
                                         int* d_phase_ispec_inner,
                                         int num_phase_ispec,
                                         int d_iphase,
                                         int use_mesh_coloring_gpu,
                                         realw* d_potential, realw* d_potential_dot_dot,
                                         realw* d_xix, realw* d_xiy, realw* d_xiz,
                                         realw* d_etax, realw* d_etay, realw* d_etaz,
                                         realw* d_gammax, realw* d_gammay, realw* d_gammaz,
                                         realw* d_hprime_xx,
                                         realw* d_hprimewgll_xx,
                                         realw* wgllwgll_xy,realw* wgllwgll_xz,realw* wgllwgll_yz,
                                         int GRAVITY,
                                         realw* d_xstore, realw* d_ystore, realw* d_zstore,
                                         realw* d_d_ln_density_dr_table,
                                         realw* d_minus_rho_g_over_kappa_fluid,
                                         realw* wgll_cube,
                                         int ROTATION,
                                         realw time,
                                         realw two_omega_earth,
                                         realw deltat,
                                         realw* d_A_array_rotation,
                                         realw* d_B_array_rotation,
                                         int NSPEC_OUTER_CORE){

  // block id == spectral-element id
  int bx = blockIdx.y*gridDim.x+blockIdx.x;
  // thread id == GLL point id
  int tx = threadIdx.x;

  // R_EARTH_KM is the radius of the bottom of the oceans (radius of Earth in km)
  //const realw R_EARTH_KM = 6371.0f;
  // uncomment line below for PREM with oceans
  //const realw R_EARTH_KM = 6368.0f;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  int active,offset;
  int iglob;
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
    s_dummy_loc[tx] = tex1Dfetch(d_displ_oc_tex, iglob);
#else
    // changing iglob indexing to match fortran row changes fast style
    s_dummy_loc[tx] = d_potential[iglob];
#endif
  }

  if (tx < NGLL2) {
    // hprime
#ifdef USE_TEXTURES_CONSTANTS
    sh_hprime_xx[tx] = tex1Dfetch(d_hprime_xx_oc_tex,tx);
#else
    sh_hprime_xx[tx] = d_hprime_xx[tx];
#endif
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
    d_potential_dot_dot[iglob] = tex1Dfetch(d_accel_oc_tex, iglob) + sum_terms;
#else
    d_potential_dot_dot[iglob] += sum_terms;
#endif // USE_TEXTURES_FIELDS

#else // MESH_COLORING

    //mesh coloring
    if( use_mesh_coloring_gpu ){

      if( NSPEC_OUTER_CORE > COLORING_MIN_NSPEC_OUTER_CORE ){
        // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
        d_potential_dot_dot[iglob] = tex1Dfetch(d_accel_oc_tex, iglob) + sum_terms;
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

/* ----------------------------------------------------------------------------------------------- */

void Kernel_2_outer_core(int nb_blocks_to_compute, Mesh* mp,
                         int d_iphase,
                         int* d_ibool,
                         realw* d_xix,realw* d_xiy,realw* d_xiz,
                         realw* d_etax,realw* d_etay,realw* d_etaz,
                         realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                         realw time,
                         realw* d_A_array_rotation,
                         realw* d_B_array_rotation,
                         realw* d_b_A_array_rotation,
                         realw* d_b_B_array_rotation,
                         int FORWARD_OR_ADJOINT){

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("before outer_core kernel Kernel_2");
#endif

  // if the grid can handle the number of blocks, we let it be 1D
  // grid_2_x = nb_elem_color;
  // nb_elem_color is just how many blocks we are computing now

  int blocksize = NGLL3_PADDED;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(nb_blocks_to_compute,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // Cuda timing
  // cudaEvent_t start, stop;
  // realw time;
  // cudaEventCreate(&start);
  // cudaEventCreate(&stop);
  // cudaEventRecord( start, 0 );

  if( FORWARD_OR_ADJOINT == 1 ){
    Kernel_2_outer_core_impl<<<grid,threads>>>(nb_blocks_to_compute,
                                                mp->NGLOB_OUTER_CORE,
                                                d_ibool,
                                                mp->d_phase_ispec_inner_outer_core,
                                                mp->num_phase_ispec_outer_core,
                                                d_iphase,
                                                mp->use_mesh_coloring_gpu,
                                                mp->d_displ_outer_core,
                                                mp->d_accel_outer_core,
                                                d_xix, d_xiy, d_xiz,
                                                d_etax, d_etay, d_etaz,
                                                d_gammax, d_gammay, d_gammaz,
                                                mp->d_hprime_xx,
                                                mp->d_hprimewgll_xx,
                                                mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                mp->gravity,
                                                mp->d_xstore_outer_core,mp->d_ystore_outer_core,mp->d_zstore_outer_core,
                                                mp->d_d_ln_density_dr_table,
                                                mp->d_minus_rho_g_over_kappa_fluid,
                                                mp->d_wgll_cube,
                                                mp->rotation,
                                                time,
                                                mp->two_omega_earth,
                                                mp->deltat,
                                                d_A_array_rotation,
                                                d_B_array_rotation,
                                                mp->NSPEC_OUTER_CORE);
  }else if( FORWARD_OR_ADJOINT == 3 ){
    // debug
    DEBUG_BACKWARD_FORCES();

    Kernel_2_outer_core_impl<<<grid,threads>>>(nb_blocks_to_compute,
                                                mp->NGLOB_OUTER_CORE,
                                                d_ibool,
                                                mp->d_phase_ispec_inner_outer_core,
                                                mp->num_phase_ispec_outer_core,
                                                d_iphase,
                                                mp->use_mesh_coloring_gpu,
                                                mp->d_b_displ_outer_core,
                                                mp->d_b_accel_outer_core,
                                                d_xix, d_xiy, d_xiz,
                                                d_etax, d_etay, d_etaz,
                                                d_gammax, d_gammay, d_gammaz,
                                                mp->d_hprime_xx,
                                                mp->d_hprimewgll_xx,
                                                mp->d_wgllwgll_xy, mp->d_wgllwgll_xz, mp->d_wgllwgll_yz,
                                                mp->gravity,
                                                mp->d_xstore_outer_core,mp->d_ystore_outer_core,mp->d_zstore_outer_core,
                                                mp->d_d_ln_density_dr_table,
                                                mp->d_minus_rho_g_over_kappa_fluid,
                                                mp->d_wgll_cube,
                                                mp->rotation,
                                                time,
                                                mp->b_two_omega_earth,
                                                mp->b_deltat,
                                                d_b_A_array_rotation,
                                                d_b_B_array_rotation,
                                                mp->NSPEC_OUTER_CORE);
  }

  // cudaEventRecord( stop, 0 );
  // cudaEventSynchronize( stop );
  // cudaEventElapsedTime( &time, start, stop );
  // cudaEventDestroy( start );
  // cudaEventDestroy( stop );
  // printf("Kernel2 Execution Time: %f ms\n",time);
  /* cudaThreadSynchronize(); */
  /* TRACE("Kernel 2 finished"); */
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //printf("Tried to start with %dx1 blocks\n",nb_blocks_to_compute);
  exit_on_cuda_error("kernel Kernel_2_outer_core");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// main compute_forces_outer_core CUDA routine

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_forces_outer_core_cuda,
              COMPUTE_FORCES_OUTER_CORE_CUDA)(long* Mesh_pointer_f,
                                              int* iphase,
                                              realw* time_f,
                                              int* FORWARD_OR_ADJOINT_f) {

  TRACE("compute_forces_outer_core_cuda");

  // debug timing
  //printf("Running compute_forces_outer_core_cuda\n");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  realw time = *time_f;

  int FORWARD_OR_ADJOINT = *FORWARD_OR_ADJOINT_f;

  int num_elements;

  if( *iphase == 1 )
    num_elements = mp->nspec_outer_outer_core;
  else
    num_elements = mp->nspec_inner_outer_core;

  if( num_elements == 0 ) return;

  // mesh coloring
  if( mp->use_mesh_coloring_gpu){

    // note: array offsets require sorted arrays, such that e.g. ibool starts with elastic elements
    //         and followed by acoustic ones.
    //         acoustic elements also start with outer than inner element ordering

    int nb_colors,nb_blocks_to_compute;
    int istart;
    int offset,offset_nonpadded;

    // sets up color loop
    if( mp->NSPEC_OUTER_CORE > COLORING_MIN_NSPEC_OUTER_CORE ){
      if( *iphase == 1 ){
        // outer elements
        nb_colors = mp->num_colors_outer_outer_core;
        istart = 0;

        // array offsets
        offset = 0;
        offset_nonpadded = 0;
      }else{
        // inner element colors (start after outer elements)
        nb_colors = mp->num_colors_outer_outer_core + mp->num_colors_inner_outer_core;
        istart = mp->num_colors_outer_outer_core;

        // array offsets (inner elements start after outer ones)
        offset = mp->nspec_outer_outer_core * NGLL3_PADDED;
        offset_nonpadded = mp->nspec_outer_outer_core * NGLL3;
      }
    }else{

      // poor element count, only use 1 color per inner/outer run

      if( *iphase == 1 ){
        // outer elements
        nb_colors = 1;
        istart = 0;

        // array offsets
        offset = 0;
        offset_nonpadded = 0;
      }else{
        // inner element colors (start after outer elements)
        nb_colors = 1;
        istart = 0;

        // array offsets (inner elements start after outer ones)
        offset = mp->nspec_outer_outer_core * NGLL3_PADDED;
        offset_nonpadded = mp->nspec_outer_outer_core * NGLL3;
      }
    }

    // loops over colors
    for(int icolor = istart; icolor < nb_colors; icolor++){

      // gets number of elements for this color
      if( mp->NSPEC_OUTER_CORE > COLORING_MIN_NSPEC_OUTER_CORE ){
        nb_blocks_to_compute = mp->h_num_elem_colors_outer_core[icolor];
      }else{
        nb_blocks_to_compute = num_elements;
      }

      Kernel_2_outer_core(nb_blocks_to_compute,mp,
                          *iphase,
                          mp->d_ibool_outer_core + offset_nonpadded,
                          mp->d_xix_outer_core + offset,mp->d_xiy_outer_core + offset,mp->d_xiz_outer_core + offset,
                          mp->d_etax_outer_core + offset,mp->d_etay_outer_core + offset,mp->d_etaz_outer_core + offset,
                          mp->d_gammax_outer_core + offset,mp->d_gammay_outer_core + offset,mp->d_gammaz_outer_core + offset,
                          time,
                          mp->d_A_array_rotation + offset_nonpadded,
                          mp->d_B_array_rotation + offset_nonpadded,
                          mp->d_b_A_array_rotation + offset_nonpadded,
                          mp->d_b_B_array_rotation + offset_nonpadded,
                          FORWARD_OR_ADJOINT);

      // for padded and aligned arrays
      offset += nb_blocks_to_compute * NGLL3_PADDED;
      // for no-aligned arrays
      offset_nonpadded += nb_blocks_to_compute * NGLL3;
    }

  }else{

    // no mesh coloring: uses atomic updates
    Kernel_2_outer_core(num_elements, mp,
                        *iphase,
                        mp->d_ibool_outer_core,
                        mp->d_xix_outer_core,mp->d_xiy_outer_core,mp->d_xiz_outer_core,
                        mp->d_etax_outer_core,mp->d_etay_outer_core,mp->d_etaz_outer_core,
                        mp->d_gammax_outer_core,mp->d_gammay_outer_core,mp->d_gammaz_outer_core,
                        time,
                        mp->d_A_array_rotation,
                        mp->d_B_array_rotation,
                        mp->d_b_A_array_rotation,
                        mp->d_b_B_array_rotation,
                        FORWARD_OR_ADJOINT);

  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("compute_forces_outer_core_cuda");
#endif
}

