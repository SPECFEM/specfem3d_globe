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

#include <sys/types.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"


/* ----------------------------------------------------------------------------------------------- */

// ELASTIC SIMULATIONS

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_kernels_rho_cudakernel(int* ibool,
                                               realw* accel,
                                               realw* b_displ,
                                               realw* rho_kl,
                                               int NSPEC,
                                               realw deltat) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;

  // handles case when there is 1 extra block (due to rectangular grid)
  if(ispec < NSPEC) {

    int ijk_ispec = threadIdx.x + NGLL3*ispec;
    int iglob = ibool[ijk_ispec] - 1 ;

    // density kernel
    rho_kl[ijk_ispec] += deltat * (accel[3*iglob]*b_displ[3*iglob]+
                                   accel[3*iglob+1]*b_displ[3*iglob+1]+
                                   accel[3*iglob+2]*b_displ[3*iglob+2]);
  }
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_kernels_iso_cudakernel(realw* epsilondev_xx,
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

/* ----------------------------------------------------------------------------------------------- */

__device__ void compute_strain_product_cuda(realw* prod,
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

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_kernels_ani_cudakernel(realw* epsilondev_xx,
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
    compute_strain_product_cuda(prod,eps_trace_over_3,epsdev,b_eps_trace_over_3,b_epsdev);

    // updates full anisotropic kernel
    for(int i=0;i<21;i++){
      cijkl_kl[i + 21*ijk_ispec] += deltat * prod[i];
    }
  }
}


/* ----------------------------------------------------------------------------------------------- */


// crust_mantle

extern "C"
void FC_FUNC_(compute_kernels_cm_cuda,
              COMPUTE_KERNELS_CM_CUDA)(long* Mesh_pointer,realw* deltat_f) {

  TRACE("compute_kernels_cm_cuda");
  // debug
  DEBUG_BACKWARD_KERNEL();

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int blocksize = NGLL3;
  realw deltat = *deltat_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->NSPEC_CRUST_MANTLE,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // density kernel
  compute_kernels_rho_cudakernel<<<grid,threads>>>(mp->d_ibool_crust_mantle,
                                                   mp->d_accel_crust_mantle,
                                                   mp->d_b_displ_crust_mantle,
                                                   mp->d_rho_kl_crust_mantle,
                                                   mp->NSPEC_CRUST_MANTLE,
                                                   deltat);

  // checks if strain is available
  if( mp->undo_attenuation ){
    exit_on_error("compute_kernels_cm_cuda not implemented yet for UNDO_ATTENUATION");
  }

  if(! mp->anisotropic_kl){
    // isotropic kernels
    compute_kernels_iso_cudakernel<<<grid,threads>>>(mp->d_epsilondev_xx_crust_mantle,
                                                     mp->d_epsilondev_yy_crust_mantle,
                                                     mp->d_epsilondev_xy_crust_mantle,
                                                     mp->d_epsilondev_xz_crust_mantle,
                                                     mp->d_epsilondev_yz_crust_mantle,
                                                     mp->d_eps_trace_over_3_crust_mantle,
                                                     mp->d_b_epsilondev_xx_crust_mantle,
                                                     mp->d_b_epsilondev_yy_crust_mantle,
                                                     mp->d_b_epsilondev_xy_crust_mantle,
                                                     mp->d_b_epsilondev_xz_crust_mantle,
                                                     mp->d_b_epsilondev_yz_crust_mantle,
                                                     mp->d_b_eps_trace_over_3_crust_mantle,
                                                     mp->d_beta_kl_crust_mantle,
                                                     mp->d_alpha_kl_crust_mantle,
                                                     mp->NSPEC_CRUST_MANTLE,
                                                     deltat);
  }else{
    // anisotropic kernels
    compute_kernels_ani_cudakernel<<<grid,threads>>>(mp->d_epsilondev_xx_crust_mantle,
                                                     mp->d_epsilondev_yy_crust_mantle,
                                                     mp->d_epsilondev_xy_crust_mantle,
                                                     mp->d_epsilondev_xz_crust_mantle,
                                                     mp->d_epsilondev_yz_crust_mantle,
                                                     mp->d_eps_trace_over_3_crust_mantle,
                                                     mp->d_b_epsilondev_xx_crust_mantle,
                                                     mp->d_b_epsilondev_yy_crust_mantle,
                                                     mp->d_b_epsilondev_xy_crust_mantle,
                                                     mp->d_b_epsilondev_xz_crust_mantle,
                                                     mp->d_b_epsilondev_yz_crust_mantle,
                                                     mp->d_b_eps_trace_over_3_crust_mantle,
                                                     mp->d_cijkl_kl_crust_mantle,
                                                     mp->NSPEC_CRUST_MANTLE,
                                                     deltat);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_kernels_cm_cuda");
#endif
}


/* ----------------------------------------------------------------------------------------------- */


// inner_core

extern "C"
void FC_FUNC_(compute_kernels_ic_cuda,
              COMPUTE_KERNELS_IC_CUDA)(long* Mesh_pointer,realw* deltat_f) {

  TRACE("compute_kernels_ic_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int blocksize = NGLL3;
  realw deltat = *deltat_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->NSPEC_INNER_CORE,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // only isotropic kernels in inner core so far implemented
  // density kernel
  compute_kernels_rho_cudakernel<<<grid,threads>>>(mp->d_ibool_inner_core,
                                                   mp->d_accel_inner_core,
                                                   mp->d_b_displ_inner_core,
                                                   mp->d_rho_kl_inner_core,
                                                   mp->NSPEC_INNER_CORE,
                                                   deltat);

  // checks if strain is available
  if( mp->undo_attenuation ){
    exit_on_error("compute_kernels_ic_cuda not implemented yet for UNDO_ATTENUATION");
  }

  // isotropic kernels (shear, bulk)
  compute_kernels_iso_cudakernel<<<grid,threads>>>(mp->d_epsilondev_xx_inner_core,
                                                   mp->d_epsilondev_yy_inner_core,
                                                   mp->d_epsilondev_xy_inner_core,
                                                   mp->d_epsilondev_xz_inner_core,
                                                   mp->d_epsilondev_yz_inner_core,
                                                   mp->d_eps_trace_over_3_inner_core,
                                                   mp->d_b_epsilondev_xx_inner_core,
                                                   mp->d_b_epsilondev_yy_inner_core,
                                                   mp->d_b_epsilondev_xy_inner_core,
                                                   mp->d_b_epsilondev_xz_inner_core,
                                                   mp->d_b_epsilondev_yz_inner_core,
                                                   mp->d_b_eps_trace_over_3_inner_core,
                                                   mp->d_beta_kl_inner_core,
                                                   mp->d_alpha_kl_inner_core,
                                                   mp->NSPEC_INNER_CORE,
                                                   deltat);


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_kernels_ic_cuda");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// ACOUSTIC SIMULATIONS

// for outer core region

/* ----------------------------------------------------------------------------------------------- */


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

/* ----------------------------------------------------------------------------------------------- */


__global__ void compute_kernels_acoustic_kernel(int* ibool,
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

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_kernels_oc_cuda,
              COMPUTE_KERNELS_OC_CUDA)(long* Mesh_pointer,realw* deltat_f) {

TRACE("compute_kernels_oc_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  int blocksize = NGLL3; // NGLLX*NGLLY*NGLLZ
  realw deltat = *deltat_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->NSPEC_OUTER_CORE,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  compute_kernels_acoustic_kernel<<<grid,threads>>>(mp->d_ibool_outer_core,
                                                    mp->d_rhostore_outer_core,
                                                    mp->d_kappavstore_outer_core,
                                                    mp->d_hprime_xx,
                                                    mp->d_xix_outer_core,
                                                    mp->d_xiy_outer_core,
                                                    mp->d_xiz_outer_core,
                                                    mp->d_etax_outer_core,
                                                    mp->d_etay_outer_core,
                                                    mp->d_etaz_outer_core,
                                                    mp->d_gammax_outer_core,
                                                    mp->d_gammay_outer_core,
                                                    mp->d_gammaz_outer_core,
                                                    mp->d_accel_outer_core,
                                                    mp->d_b_displ_outer_core,
                                                    mp->d_b_accel_outer_core,
                                                    mp->d_rho_kl_outer_core,
                                                    mp->d_alpha_kl_outer_core,
                                                    deltat,
                                                    mp->NSPEC_OUTER_CORE);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_kernels_oc_kernel");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// NOISE SIMULATIONS

/* ----------------------------------------------------------------------------------------------- */


__global__ void compute_kernels_strength_noise_cuda_kernel(realw* displ,
                                                           int* ibelm_top,
                                                           int* ibool,
                                                           realw* noise_surface_movie,
                                                           realw* normal_x_noise,
                                                           realw* normal_y_noise,
                                                           realw* normal_z_noise,
                                                           realw* Sigma_kl,
                                                           realw deltat,
                                                           int nspec_top) {
  int iface = blockIdx.x + blockIdx.y*gridDim.x;

  if(iface < nspec_top) {

    int ispec = ibelm_top[iface]-1;
    int igll = threadIdx.x;
    int ipoin = igll + NGLL2*iface;

    int k = NGLLX-1;
    int j = (igll/NGLLX);
    int i = (igll-j*NGLLX);

    int iglob = ibool[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] - 1 ;

    realw eta = ( noise_surface_movie[INDEX3(NDIM,NGLL2,0,igll,iface)]*normal_x_noise[ipoin]+
                 noise_surface_movie[INDEX3(NDIM,NGLL2,1,igll,iface)]*normal_y_noise[ipoin]+
                 noise_surface_movie[INDEX3(NDIM,NGLL2,2,igll,iface)]*normal_z_noise[ipoin]);

    Sigma_kl[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] += deltat*eta*
                                                      (normal_x_noise[ipoin]*displ[3*iglob]+
                                                       normal_y_noise[ipoin]*displ[1+3*iglob]+
                                                       normal_z_noise[ipoin]*displ[2+3*iglob]);
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_kernels_strgth_noise_cu,
              COMPUTE_KERNELS_STRGTH_NOISE_CU)(long* Mesh_pointer,
                                               realw* h_noise_surface_movie,
                                               realw* deltat_f) {

  TRACE("compute_kernels_strgth_noise_cu");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  realw deltat = *deltat_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nspec2D_top_crust_mantle,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLL2,1,1);

  // copies surface buffer to GPU
  print_CUDA_error_if_any(cudaMemcpy(mp->d_noise_surface_movie,h_noise_surface_movie,
                                     NDIM*NGLL2*(mp->nspec2D_top_crust_mantle)*sizeof(realw),
                                     cudaMemcpyHostToDevice),90900);

  // calculates noise strength kernel
  compute_kernels_strength_noise_cuda_kernel<<<grid,threads>>>(mp->d_displ_crust_mantle,
                                                               mp->d_ibelm_top_crust_mantle,
                                                               mp->d_ibool_crust_mantle,
                                                               mp->d_noise_surface_movie,
                                                               mp->d_normal_x_noise,
                                                               mp->d_normal_y_noise,
                                                               mp->d_normal_z_noise,
                                                               mp->d_Sigma_kl,
                                                               deltat,
                                                               mp->nspec2D_top_crust_mantle);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_kernels_strength_noise_cuda_kernel");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// preconditioner (approximate Hessian kernel)

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_kernels_hess_cudakernel(int* ibool,
                                                realw* accel,
                                                realw* b_accel,
                                                realw* hess_kl,
                                                realw deltat,
                                                int NSPEC_AB) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;

  // handles case when there is 1 extra block (due to rectangular grid)
  if(ispec < NSPEC_AB) {

    int ijk = threadIdx.x;
    int ijk_ispec = ijk + NGLL3*ispec;
    int iglob = ibool[ijk_ispec] - 1 ;

    // approximate hessian
    hess_kl[ijk_ispec] += deltat * (accel[3*iglob]*b_accel[3*iglob] +
                                    accel[3*iglob+1]*b_accel[3*iglob+1] +
                                    accel[3*iglob+2]*b_accel[3*iglob+2]);
  }
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_kernels_hess_cuda,
              COMPUTE_KERNELS_HESS_CUDA)(long* Mesh_pointer,
                                         realw* deltat_f) {
  TRACE("compute_kernels_hess_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks
  if( ! mp->approximate_hess_kl ){exit_on_cuda_error("approximate_hess_kl flag not properly initialized");}

  int blocksize = NGLL3; // NGLLX*NGLLY*NGLLZ
  realw deltat = *deltat_f;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->NSPEC_CRUST_MANTLE,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  compute_kernels_hess_cudakernel<<<grid,threads>>>(mp->d_ibool_crust_mantle,
                                                    mp->d_accel_crust_mantle,
                                                    mp->d_b_accel_crust_mantle,
                                                    mp->d_hess_kl_crust_mantle,
                                                    deltat,
                                                    mp->NSPEC_CRUST_MANTLE);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_kernels_hess_cuda");
#endif
}

