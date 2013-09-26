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
#include "prepare_constants_cuda.h"

#ifdef USE_OLDER_CUDA4_GPU
#else
  #ifdef USE_TEXTURES_FIELDS
extern texture<realw, cudaTextureType1D, cudaReadModeElementType> d_displ_cm_tex;
extern texture<realw, cudaTextureType1D, cudaReadModeElementType> d_accel_cm_tex;

extern texture<realw, cudaTextureType1D, cudaReadModeElementType> d_displ_oc_tex;
extern texture<realw, cudaTextureType1D, cudaReadModeElementType> d_accel_oc_tex;

extern texture<realw, cudaTextureType1D, cudaReadModeElementType> d_displ_ic_tex;
extern texture<realw, cudaTextureType1D, cudaReadModeElementType> d_accel_ic_tex;
  #endif

  #ifdef USE_TEXTURES_CONSTANTS
extern texture<realw, cudaTextureType1D, cudaReadModeElementType> d_hprime_xx_cm_tex;
extern texture<realw, cudaTextureType1D, cudaReadModeElementType> d_hprime_xx_oc_tex;
extern texture<realw, cudaTextureType1D, cudaReadModeElementType> d_hprime_xx_ic_tex;
  #endif
#endif

/* ----------------------------------------------------------------------------------------------- */

// helper functions

/* ----------------------------------------------------------------------------------------------- */


// copies integer array from CPU host to GPU device
void copy_todevice_int(void** d_array_addr_ptr,int* h_array,int size){
  TRACE("copy_todevice_int");

  // allocates memory on GPU
  //
  // note: cudaMalloc uses a double-pointer, such that it can return an error code in case it fails
  //          we thus pass the address to the pointer above (as void double-pointer) to have it
  //          pointing to the correct pointer of the array here
  print_CUDA_error_if_any(cudaMalloc((void**)d_array_addr_ptr,size*sizeof(int)),
                          12001);

  // copies values onto GPU
  //
  // note: cudaMemcpy uses the pointer to the array, we thus re-cast the value of
  //          the double-pointer above to have the correct pointer to the array
  print_CUDA_error_if_any(cudaMemcpy((int*) *d_array_addr_ptr,h_array,size*sizeof(int),cudaMemcpyHostToDevice),
                          12002);
}

/* ----------------------------------------------------------------------------------------------- */

// copies integer array from CPU host to GPU device
void copy_todevice_realw(void** d_array_addr_ptr,realw* h_array,int size){
  TRACE("copy_todevice_realw");

  // allocates memory on GPU
  print_CUDA_error_if_any(cudaMalloc((void**)d_array_addr_ptr,size*sizeof(realw)),
                          22001);
  // copies values onto GPU
  print_CUDA_error_if_any(cudaMemcpy((realw*) *d_array_addr_ptr,h_array,size*sizeof(realw),cudaMemcpyHostToDevice),
                          22002);
}


/* ----------------------------------------------------------------------------------------------- */

// GPU preparation

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_constants_device,
              PREPARE_CONSTANTS_DEVICE)(long* Mesh_pointer,
                                        int* myrank_f,
                                        int* h_NGLLX,
                                        realw* h_hprime_xx,realw* h_hprimewgll_xx,
                                        realw* h_wgllwgll_xy,realw* h_wgllwgll_xz,realw* h_wgllwgll_yz,
                                        int* NSOURCES,int* nsources_local,
                                        realw* h_sourcearrays,
                                        int* h_islice_selected_source,int* h_ispec_selected_source,
                                        int* nrec,int* nrec_local, int* nadj_rec_local,
                                        int* h_number_receiver_global,
                                        int* h_islice_selected_rec,int* h_ispec_selected_rec,
                                        int* NSPEC_CRUST_MANTLE, int* NGLOB_CRUST_MANTLE,
                                        int* NSPEC_CRUST_MANTLE_STRAIN_ONLY,
                                        int* NSPEC_OUTER_CORE, int* NGLOB_OUTER_CORE,
                                        int* NSPEC_INNER_CORE, int* NGLOB_INNER_CORE,
                                        int* NSPEC_INNER_CORE_STRAIN_ONLY,
                                        int* SIMULATION_TYPE,int* NOISE_TOMOGRAPHY,
                                        int* SAVE_FORWARD_f,int* ABSORBING_CONDITIONS_f,
                                        int* OCEANS_f,int* GRAVITY_f,
                                        int* ROTATION_f,int* EXACT_MASS_MATRIX_FOR_ROTATION_f,
                                        int* ATTENUATION_f,int* UNDO_ATTENUATION_f,
                                        int* PARTIAL_PHYS_DISPERSION_ONLY_f,int* USE_3D_ATTENUATION_ARRAYS_f,
                                        int* COMPUTE_AND_STORE_STRAIN_f,
                                        int* ANISOTROPIC_3D_MANTLE_f,int* ANISOTROPIC_INNER_CORE_f,
                                        int* SAVE_BOUNDARY_MESH_f,
                                        int* USE_MESH_COLORING_GPU_f,
                                        int* ANISOTROPIC_KL_f,int* APPROXIMATE_HESS_KL_f,
                                        realw* deltat_f,realw* b_deltat_f) {

  TRACE("prepare_constants_device");

  // allocates mesh parameter structure
  Mesh* mp = (Mesh*) malloc( sizeof(Mesh) );
  if (mp == NULL) exit_on_error("error allocating mesh pointer");
  *Mesh_pointer = (long)mp;

  // checks if NGLLX == 5
  if( *h_NGLLX != NGLLX ){
    exit_on_error("NGLLX must be 5 for CUDA devices");
  }

  // sets constant arrays
  setConst_hprime_xx(h_hprime_xx,mp);
  //setConst_hprime_yy(h_hprime_yy,mp); // only needed if NGLLX != NGLLY != NGLLZ
  //setConst_hprime_zz(h_hprime_zz,mp); // only needed if NGLLX != NGLLY != NGLLZ

  setConst_hprimewgll_xx(h_hprimewgll_xx,mp);
  //setConst_hprimewgll_yy(h_hprimewgll_yy,mp); // only needed if NGLLX != NGLLY != NGLLZ
  //setConst_hprimewgll_zz(h_hprimewgll_zz,mp); // only needed if NGLLX != NGLLY != NGLLZ

  setConst_wgllwgll_xy(h_wgllwgll_xy,mp);
  setConst_wgllwgll_xz(h_wgllwgll_xz,mp);
  setConst_wgllwgll_yz(h_wgllwgll_yz,mp);

  // Using texture memory for the hprime-style constants is slower on
  // Fermi generation hardware, but *may* be faster on Kepler
  // generation. We will reevaluate this again, so might as well leave
  // in the code with #USE_TEXTURES_FIELDS not-defined.
  #ifdef USE_TEXTURES_CONSTANTS
  {
    #ifdef USE_OLDER_CUDA4_GPU
      const textureReference* d_hprime_xx_cm_tex_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_hprime_xx_cm_tex_ptr, "d_hprime_xx_cm_tex"), 1101);
      cudaChannelFormatDesc channelDesc1 = cudaCreateChannelDesc<realw>();
      print_CUDA_error_if_any(cudaBindTexture(0, d_hprime_xx_cm_tex_ptr, mp->d_hprime_xx,
                                              &channelDesc1, sizeof(realw)*(NGLL2)), 1102);

      const textureReference* d_hprime_xx_oc_tex_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_hprime_xx_oc_tex_ptr, "d_hprime_xx_oc_tex"), 1103);
      cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<realw>();
      print_CUDA_error_if_any(cudaBindTexture(0, d_hprime_xx_oc_tex_ptr, mp->d_hprime_xx,
                                              &channelDesc2, sizeof(realw)*(NGLL2)), 1104);

      const textureReference* d_hprime_xx_ic_tex_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_hprime_xx_ic_tex_ptr, "d_hprime_xx_ic_tex"), 1105);
      cudaChannelFormatDesc channelDesc3 = cudaCreateChannelDesc<realw>();
      print_CUDA_error_if_any(cudaBindTexture(0, d_hprime_xx_ic_tex_ptr, mp->d_hprime_xx,
                                              &channelDesc3, sizeof(realw)*(NGLL2)), 1106);
    #else
      cudaChannelFormatDesc channelDesc1 = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_hprime_xx_cm_tex, mp->d_hprime_xx,
                                              &channelDesc1, sizeof(realw)*(NGLL2)), 1102);

      cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_hprime_xx_oc_tex, mp->d_hprime_xx,
                                              &channelDesc2, sizeof(realw)*(NGLL2)), 1104);

      cudaChannelFormatDesc channelDesc3 = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_hprime_xx_ic_tex, mp->d_hprime_xx,
                                              &channelDesc3, sizeof(realw)*(NGLL2)), 1106);
    #endif
  }
  #endif


  // sets global parameters
  mp->NSPEC_CRUST_MANTLE = *NSPEC_CRUST_MANTLE;
  mp->NGLOB_CRUST_MANTLE = *NGLOB_CRUST_MANTLE;
  mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY = *NSPEC_CRUST_MANTLE_STRAIN_ONLY;

  mp->NSPEC_OUTER_CORE = *NSPEC_OUTER_CORE;
  mp->NGLOB_OUTER_CORE = *NGLOB_OUTER_CORE;

  mp->NSPEC_INNER_CORE = *NSPEC_INNER_CORE;
  mp->NGLOB_INNER_CORE = *NGLOB_INNER_CORE;
  mp->NSPEC_INNER_CORE_STRAIN_ONLY = *NSPEC_INNER_CORE_STRAIN_ONLY;

  // simulation type
  mp->simulation_type = *SIMULATION_TYPE;
  mp->noise_tomography = *NOISE_TOMOGRAPHY;

  // simulation flags initialization
  mp->save_forward = *SAVE_FORWARD_f;
  mp->absorbing_conditions = *ABSORBING_CONDITIONS_f;
  mp->oceans = *OCEANS_f;
  mp->gravity = *GRAVITY_f;
  mp->rotation = *ROTATION_f;
  mp->exact_mass_matrix_for_rotation = *EXACT_MASS_MATRIX_FOR_ROTATION_f;

  mp->attenuation = *ATTENUATION_f;
  mp->undo_attenuation = *UNDO_ATTENUATION_f;
  mp->partial_phys_dispersion_only = *PARTIAL_PHYS_DISPERSION_ONLY_f;
  mp->use_3d_attenuation_arrays = *USE_3D_ATTENUATION_ARRAYS_f;

  mp->compute_and_store_strain = *COMPUTE_AND_STORE_STRAIN_f;
  mp->anisotropic_3D_mantle = *ANISOTROPIC_3D_MANTLE_f;
  mp->anisotropic_inner_core = *ANISOTROPIC_INNER_CORE_f;
  mp->save_boundary_mesh = *SAVE_BOUNDARY_MESH_f;

  mp->anisotropic_kl = *ANISOTROPIC_KL_f;
  mp->approximate_hess_kl = *APPROXIMATE_HESS_KL_f;

  // mpi process rank
  mp->myrank = *myrank_f;

  // mesh coloring flag
#ifdef USE_MESH_COLORING_GPU
  mp->use_mesh_coloring_gpu = 1;
  if( ! *USE_MESH_COLORING_GPU_f ){exit_on_error("error with USE_MESH_COLORING_GPU constant; please re-compile\n");}
#else
  // mesh coloring
  // note: this here passes the coloring as an option to the kernel routines
  //          the performance seems to be the same if one uses the pre-processing directives above or not
  mp->use_mesh_coloring_gpu = *USE_MESH_COLORING_GPU_f;
#endif

  // sources
  mp->nsources_local = *nsources_local;
  if( mp->simulation_type == 1  || mp->simulation_type == 3 ){
    // not needed in case of pure adjoint simulations (SIMULATION_TYPE == 2)
    copy_todevice_realw((void**)&mp->d_sourcearrays,h_sourcearrays,(*NSOURCES)*NDIM*NGLL3);

    // buffer for source time function values
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_stf_pre_compute,
                                       (*NSOURCES)*sizeof(double)),1303);
  }
  copy_todevice_int((void**)&mp->d_islice_selected_source,h_islice_selected_source,(*NSOURCES));
  copy_todevice_int((void**)&mp->d_ispec_selected_source,h_ispec_selected_source,(*NSOURCES));

  // receiver stations
  // note that:   size(number_receiver_global) = nrec_local
  //                   size(ispec_selected_rec) = nrec
  // number of receiver located in this partition
  mp->nrec_local = *nrec_local;
  if( mp->nrec_local > 0 ){
    copy_todevice_int((void**)&mp->d_number_receiver_global,h_number_receiver_global,mp->nrec_local);

    // for seismograms
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_station_seismo_field),
                                       NDIM*NGLL3*(mp->nrec_local)*sizeof(realw)),4015);

    mp->h_station_seismo_field = (realw*) malloc( NDIM*NGLL3*(mp->nrec_local)*sizeof(realw) );
    if( mp->h_station_seismo_field == NULL) exit_on_error("h_station_seismo_field not allocated \n");

  }
  copy_todevice_int((void**)&mp->d_ispec_selected_rec,h_ispec_selected_rec,(*nrec));

  // receiver adjoint source arrays only used for noise and adjoint simulations
  // adjoint source arrays
  mp->nadj_rec_local = *nadj_rec_local;
  if( mp->nadj_rec_local > 0 ){

    // prepares local irec array:
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_pre_computed_irec,
                                       (mp->nadj_rec_local)*sizeof(int)),6004);

    // the irec_local variable needs to be precomputed (as
    // h_pre_comp..), because normally it is in the loop updating accel,
    // and due to how it's incremented, it cannot be parallelized
    int* h_pre_computed_irec = (int*) malloc( (mp->nadj_rec_local)*sizeof(int) );
    if( h_pre_computed_irec == NULL ) exit_on_error("h_pre_computed_irec not allocated\n");

    int irec_local = 0;
    for(int irec = 0; irec < *nrec; irec++) {
      if(mp->myrank == h_islice_selected_rec[irec]) {
        irec_local++;
        h_pre_computed_irec[irec_local-1] = irec;
      }
    }
    if( irec_local != mp->nadj_rec_local ) exit_on_error("prepare_sim2_or_3_const_device: irec_local not equal\n");
    // copies values onto GPU
    print_CUDA_error_if_any(cudaMemcpy(mp->d_pre_computed_irec,h_pre_computed_irec,
                                       (mp->nadj_rec_local)*sizeof(int),cudaMemcpyHostToDevice),6010);
    free(h_pre_computed_irec);

    // temporary array to prepare extracted source array values
    mp->h_adj_sourcearrays_slice = (realw*) malloc( (mp->nadj_rec_local)*NDIM*NGLL3*sizeof(realw) );
    if( mp->h_adj_sourcearrays_slice == NULL ) exit_on_error("h_adj_sourcearrays_slice not allocated\n");

    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_adj_sourcearrays,
                                       (mp->nadj_rec_local)*NDIM*NGLL3*sizeof(realw)),6003);

  }

  // for rotation and new attenuation
  mp->deltat = *deltat_f;
  if( mp->simulation_type == 3 ){
    mp->b_deltat = *b_deltat_f;
  }
  // initializes for rotational effects
  mp->two_omega_earth = 0.f;
  mp->b_two_omega_earth = 0.f;

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_constants_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// ROTATION simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_rotation_device,
              PREPARE_FIELDS_ROTATION_DEVICE)(long* Mesh_pointer_f,
                                              realw* two_omega_earth_f,
                                              realw* A_array_rotation,
                                              realw* B_array_rotation,
                                              realw* b_two_omega_earth_f,
                                              realw* b_A_array_rotation,
                                              realw* b_B_array_rotation,
                                              int* NSPEC_OUTER_CORE_ROTATION) {

  TRACE("prepare_fields_rotation_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // arrays only needed when rotation is required
  if( ! mp->rotation ){
    exit_on_cuda_error("prepare_fields_rotation_device: rotation flag not properly initialized");
  }
  // checks array size
  if( *NSPEC_OUTER_CORE_ROTATION != mp->NSPEC_OUTER_CORE){
    printf("error prepare_fields_rotation_device: rotation array has wrong size: %d instead of %d\n",
           *NSPEC_OUTER_CORE_ROTATION,mp->NSPEC_OUTER_CORE);
    exit_on_cuda_error("prepare_fields_rotation_device: rotation array has wrong size");
  }

  // rotation arrays (needed only for outer core region)
  mp->two_omega_earth = *two_omega_earth_f;
  copy_todevice_realw((void**)&mp->d_A_array_rotation,A_array_rotation,NGLL3*mp->NSPEC_OUTER_CORE);
  copy_todevice_realw((void**)&mp->d_B_array_rotation,B_array_rotation,NGLL3*mp->NSPEC_OUTER_CORE);

  // backward/reconstructed fields
  if( mp->simulation_type == 3 ){
    mp->b_two_omega_earth = *b_two_omega_earth_f;
    copy_todevice_realw((void**)&mp->d_b_A_array_rotation,b_A_array_rotation,NGLL3*mp->NSPEC_OUTER_CORE);
    copy_todevice_realw((void**)&mp->d_b_B_array_rotation,b_B_array_rotation,NGLL3*mp->NSPEC_OUTER_CORE);
  }
}


/* ----------------------------------------------------------------------------------------------- */

// GRAVITY simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_gravity_device,
              PREPARE_FIELDS_gravity_DEVICE)(long* Mesh_pointer_f,
                                             realw* d_ln_density_dr_table,
                                             realw* minus_rho_g_over_kappa_fluid,
                                             realw* minus_gravity_table,
                                             realw* minus_deriv_gravity_table,
                                             realw* density_table,
                                             realw* h_wgll_cube,
                                             int* NRAD_GRAVITY,
                                             realw* minus_g_icb,
                                             realw* minus_g_cmb,
                                             double* RHO_BOTTOM_OC,
                                             double* RHO_TOP_OC) {

  TRACE("prepare_fields_gravity_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  if( ! mp->gravity ){
    // no gravity case

    // d ln(rho)/dr needed for the no gravity fluid potential
    copy_todevice_realw((void**)&mp->d_d_ln_density_dr_table,d_ln_density_dr_table,(*NRAD_GRAVITY));

  }else{
    // gravity case
    mp->minus_g_icb = *minus_g_icb;
    mp->minus_g_cmb = *minus_g_cmb;

    // sets up gll weights cubed
    setConst_wgll_cube(h_wgll_cube,mp);

    // prepares gravity arrays
    copy_todevice_realw((void**)&mp->d_minus_rho_g_over_kappa_fluid,minus_rho_g_over_kappa_fluid,(*NRAD_GRAVITY));
    copy_todevice_realw((void**)&mp->d_minus_gravity_table,minus_gravity_table,(*NRAD_GRAVITY));
    copy_todevice_realw((void**)&mp->d_minus_deriv_gravity_table,minus_deriv_gravity_table,(*NRAD_GRAVITY));
    copy_todevice_realw((void**)&mp->d_density_table,density_table,(*NRAD_GRAVITY));
  }

  // constants
  mp->RHO_BOTTOM_OC = (realw) *RHO_BOTTOM_OC;
  mp->RHO_TOP_OC = (realw) *RHO_TOP_OC;

}



/* ----------------------------------------------------------------------------------------------- */

// ATTENUATION simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_attenuat_device,
              PREPARE_FIELDS_ATTENUAT_DEVICE)(long* Mesh_pointer_f,
                                                 realw* R_xx_crust_mantle,
                                                 realw* R_yy_crust_mantle,
                                                 realw* R_xy_crust_mantle,
                                                 realw* R_xz_crust_mantle,
                                                 realw* R_yz_crust_mantle,
                                                 realw* b_R_xx_crust_mantle,
                                                 realw* b_R_yy_crust_mantle,
                                                 realw* b_R_xy_crust_mantle,
                                                 realw* b_R_xz_crust_mantle,
                                                 realw* b_R_yz_crust_mantle,
                                                 realw* factor_common_crust_mantle,
                                                 realw* one_minus_sum_beta_crust_mantle,
                                                 realw* R_xx_inner_core,
                                                 realw* R_yy_inner_core,
                                                 realw* R_xy_inner_core,
                                                 realw* R_xz_inner_core,
                                                 realw* R_yz_inner_core,
                                                 realw* b_R_xx_inner_core,
                                                 realw* b_R_yy_inner_core,
                                                 realw* b_R_xy_inner_core,
                                                 realw* b_R_xz_inner_core,
                                                 realw* b_R_yz_inner_core,
                                                 realw* factor_common_inner_core,
                                                 realw* one_minus_sum_beta_inner_core,
                                                 realw* alphaval,realw* betaval,realw* gammaval,
                                                 realw* b_alphaval,realw* b_betaval,realw* b_gammaval) {

  TRACE("prepare_fields_attenuat_device");
  int R_size1,R_size2,R_size3;

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // checks flag
  if( ! mp->attenuation ){ exit_on_cuda_error("prepare_fields_attenuat_device attenuation not properly initialized"); }

  // crust_mantle
  R_size1 = N_SLS*NGLL3*mp->NSPEC_CRUST_MANTLE;
  if( mp->use_3d_attenuation_arrays ){
    R_size2 = NGLL3*mp->NSPEC_CRUST_MANTLE;
    R_size3 = N_SLS*NGLL3*mp->NSPEC_CRUST_MANTLE;
  }else{
    R_size2 = 1*mp->NSPEC_CRUST_MANTLE;
    R_size3 = N_SLS*1*mp->NSPEC_CRUST_MANTLE;
  }

  copy_todevice_realw((void**)&mp->d_one_minus_sum_beta_crust_mantle,one_minus_sum_beta_crust_mantle,R_size2);

  if( ! mp->partial_phys_dispersion_only ){
    // common factor
    copy_todevice_realw((void**)&mp->d_factor_common_crust_mantle,factor_common_crust_mantle,R_size3);
    // memory variables
    copy_todevice_realw((void**)&mp->d_R_xx_crust_mantle,R_xx_crust_mantle,R_size1);
    copy_todevice_realw((void**)&mp->d_R_yy_crust_mantle,R_yy_crust_mantle,R_size1);
    copy_todevice_realw((void**)&mp->d_R_xy_crust_mantle,R_xy_crust_mantle,R_size1);
    copy_todevice_realw((void**)&mp->d_R_xz_crust_mantle,R_xz_crust_mantle,R_size1);
    copy_todevice_realw((void**)&mp->d_R_yz_crust_mantle,R_yz_crust_mantle,R_size1);
  }

  if(mp->simulation_type == 3 ){
    if( ! mp->partial_phys_dispersion_only ){
      // memory variables
      copy_todevice_realw((void**)&mp->d_b_R_xx_crust_mantle,b_R_xx_crust_mantle,R_size1);
      copy_todevice_realw((void**)&mp->d_b_R_yy_crust_mantle,b_R_yy_crust_mantle,R_size1);
      copy_todevice_realw((void**)&mp->d_b_R_xy_crust_mantle,b_R_xy_crust_mantle,R_size1);
      copy_todevice_realw((void**)&mp->d_b_R_xz_crust_mantle,b_R_xz_crust_mantle,R_size1);
      copy_todevice_realw((void**)&mp->d_b_R_yz_crust_mantle,b_R_yz_crust_mantle,R_size1);
    }
  }

  // inner_core
  R_size1 = N_SLS*NGLL3*mp->NSPEC_INNER_CORE;
  if( mp->use_3d_attenuation_arrays ){
    R_size2 = NGLL3*mp->NSPEC_INNER_CORE;
    R_size3 = N_SLS*NGLL3*mp->NSPEC_INNER_CORE;
  }else{
    R_size2 = 1*mp->NSPEC_INNER_CORE;
    R_size3 = N_SLS*1*mp->NSPEC_INNER_CORE;
  }

  copy_todevice_realw((void**)&mp->d_one_minus_sum_beta_inner_core,one_minus_sum_beta_inner_core,R_size2);

  if( ! mp->partial_phys_dispersion_only ){
    // common factor
    copy_todevice_realw((void**)&mp->d_factor_common_inner_core,factor_common_inner_core,R_size3);
    // memory variables
    copy_todevice_realw((void**)&mp->d_R_xx_inner_core,R_xx_inner_core,R_size1);
    copy_todevice_realw((void**)&mp->d_R_yy_inner_core,R_yy_inner_core,R_size1);
    copy_todevice_realw((void**)&mp->d_R_xy_inner_core,R_xy_inner_core,R_size1);
    copy_todevice_realw((void**)&mp->d_R_xz_inner_core,R_xz_inner_core,R_size1);
    copy_todevice_realw((void**)&mp->d_R_yz_inner_core,R_yz_inner_core,R_size1);
  }

  if(mp->simulation_type == 3 ){
    if( ! mp->partial_phys_dispersion_only ){
      // memory variables
      copy_todevice_realw((void**)&mp->d_b_R_xx_inner_core,b_R_xx_inner_core,R_size1);
      copy_todevice_realw((void**)&mp->d_b_R_yy_inner_core,b_R_yy_inner_core,R_size1);
      copy_todevice_realw((void**)&mp->d_b_R_xy_inner_core,b_R_xy_inner_core,R_size1);
      copy_todevice_realw((void**)&mp->d_b_R_xz_inner_core,b_R_xz_inner_core,R_size1);
      copy_todevice_realw((void**)&mp->d_b_R_yz_inner_core,b_R_yz_inner_core,R_size1);
    }
  }

  // alpha,beta,gamma factors
  copy_todevice_realw((void**)&mp->d_alphaval,alphaval,N_SLS);
  copy_todevice_realw((void**)&mp->d_betaval,betaval,N_SLS);
  copy_todevice_realw((void**)&mp->d_gammaval,gammaval,N_SLS);

  if( mp->simulation_type == 3 ){
    // alpha,beta,gamma factors for backward fields
    copy_todevice_realw((void**)&mp->d_b_alphaval,b_alphaval,N_SLS);
    copy_todevice_realw((void**)&mp->d_b_betaval,b_betaval,N_SLS);
    copy_todevice_realw((void**)&mp->d_b_gammaval,b_gammaval,N_SLS);
  }
}

/* ----------------------------------------------------------------------------------------------- */

// STRAIN simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_strain_device,
              PREPARE_FIELDS_STRAIN_DEVICE)(long* Mesh_pointer_f,
                                            realw* epsilondev_xx_crust_mantle,
                                            realw* epsilondev_yy_crust_mantle,
                                            realw* epsilondev_xy_crust_mantle,
                                            realw* epsilondev_xz_crust_mantle,
                                            realw* epsilondev_yz_crust_mantle,
                                            realw* b_epsilondev_xx_crust_mantle,
                                            realw* b_epsilondev_yy_crust_mantle,
                                            realw* b_epsilondev_xy_crust_mantle,
                                            realw* b_epsilondev_xz_crust_mantle,
                                            realw* b_epsilondev_yz_crust_mantle,
                                            realw* eps_trace_over_3_crust_mantle,
                                            realw* b_eps_trace_over_3_crust_mantle,
                                            realw* epsilondev_xx_inner_core,
                                            realw* epsilondev_yy_inner_core,
                                            realw* epsilondev_xy_inner_core,
                                            realw* epsilondev_xz_inner_core,
                                            realw* epsilondev_yz_inner_core,
                                            realw* b_epsilondev_xx_inner_core,
                                            realw* b_epsilondev_yy_inner_core,
                                            realw* b_epsilondev_xy_inner_core,
                                            realw* b_epsilondev_xz_inner_core,
                                            realw* b_epsilondev_yz_inner_core,
                                            realw* eps_trace_over_3_inner_core,
                                            realw* b_eps_trace_over_3_inner_core) {

  TRACE("prepare_fields_strain_device");
  int R_size,size_strain_only;

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // checks flag
  if( ! mp->compute_and_store_strain ){ exit_on_cuda_error("prepare_fields_strain_device strain not properly initialized"); }

  // crust_mantle
  R_size = NGLL3*mp->NSPEC_CRUST_MANTLE;
  copy_todevice_realw((void**)&mp->d_epsilondev_xx_crust_mantle,epsilondev_xx_crust_mantle,R_size);
  copy_todevice_realw((void**)&mp->d_epsilondev_yy_crust_mantle,epsilondev_yy_crust_mantle,R_size);
  copy_todevice_realw((void**)&mp->d_epsilondev_xy_crust_mantle,epsilondev_xy_crust_mantle,R_size);
  copy_todevice_realw((void**)&mp->d_epsilondev_xz_crust_mantle,epsilondev_xz_crust_mantle,R_size);
  copy_todevice_realw((void**)&mp->d_epsilondev_yz_crust_mantle,epsilondev_yz_crust_mantle,R_size);

  // strain
  size_strain_only = NGLL3*(mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY);
  copy_todevice_realw((void**)&mp->d_eps_trace_over_3_crust_mantle,eps_trace_over_3_crust_mantle,size_strain_only);

  // backward/reconstructed fields
  if( mp->simulation_type == 3 ){
    if( mp->undo_attenuation ){
      // strain will be computed locally based on displacement wavefield
      // only uses pointers to already allocated arrays
      mp->d_b_epsilondev_xx_crust_mantle = mp->d_epsilondev_xx_crust_mantle;
      mp->d_b_epsilondev_yy_crust_mantle = mp->d_epsilondev_yy_crust_mantle;
      mp->d_b_epsilondev_xy_crust_mantle = mp->d_epsilondev_xy_crust_mantle;
      mp->d_b_epsilondev_xz_crust_mantle = mp->d_epsilondev_xz_crust_mantle;
      mp->d_b_epsilondev_yz_crust_mantle = mp->d_epsilondev_yz_crust_mantle;
      mp->d_b_eps_trace_over_3_crust_mantle = mp->d_eps_trace_over_3_crust_mantle;
    }else{
      copy_todevice_realw((void**)&mp->d_b_epsilondev_xx_crust_mantle,b_epsilondev_xx_crust_mantle,R_size);
      copy_todevice_realw((void**)&mp->d_b_epsilondev_yy_crust_mantle,b_epsilondev_yy_crust_mantle,R_size);
      copy_todevice_realw((void**)&mp->d_b_epsilondev_xy_crust_mantle,b_epsilondev_xy_crust_mantle,R_size);
      copy_todevice_realw((void**)&mp->d_b_epsilondev_xz_crust_mantle,b_epsilondev_xz_crust_mantle,R_size);
      copy_todevice_realw((void**)&mp->d_b_epsilondev_yz_crust_mantle,b_epsilondev_yz_crust_mantle,R_size);
      //strain
      copy_todevice_realw((void**)&mp->d_b_eps_trace_over_3_crust_mantle,b_eps_trace_over_3_crust_mantle,R_size);
    }
  }

  // inner_core
  R_size = NGLL3*mp->NSPEC_INNER_CORE;
  copy_todevice_realw((void**)&mp->d_epsilondev_xx_inner_core,epsilondev_xx_inner_core,R_size);
  copy_todevice_realw((void**)&mp->d_epsilondev_yy_inner_core,epsilondev_yy_inner_core,R_size);
  copy_todevice_realw((void**)&mp->d_epsilondev_xy_inner_core,epsilondev_xy_inner_core,R_size);
  copy_todevice_realw((void**)&mp->d_epsilondev_xz_inner_core,epsilondev_xz_inner_core,R_size);
  copy_todevice_realw((void**)&mp->d_epsilondev_yz_inner_core,epsilondev_yz_inner_core,R_size);

  // strain
  size_strain_only = NGLL3*(mp->NSPEC_INNER_CORE_STRAIN_ONLY);
  copy_todevice_realw((void**)&mp->d_eps_trace_over_3_inner_core,eps_trace_over_3_inner_core,size_strain_only);

  // backward/reconstructed fields
  if( mp->simulation_type == 3 ){
    if( mp->undo_attenuation ){
      // strain will be computed locally based on displacement wavefield
      // only uses pointers to already allocated arrays
      mp->d_b_epsilondev_xx_inner_core = mp->d_epsilondev_xx_inner_core;
      mp->d_b_epsilondev_yy_inner_core = mp->d_epsilondev_yy_inner_core;
      mp->d_b_epsilondev_xy_inner_core = mp->d_epsilondev_xy_inner_core;
      mp->d_b_epsilondev_xz_inner_core = mp->d_epsilondev_xz_inner_core;
      mp->d_b_epsilondev_yz_inner_core = mp->d_epsilondev_yz_inner_core;
      mp->d_b_eps_trace_over_3_inner_core = mp->d_eps_trace_over_3_inner_core;
    }else{
      copy_todevice_realw((void**)&mp->d_b_epsilondev_xx_inner_core,b_epsilondev_xx_inner_core,R_size);
      copy_todevice_realw((void**)&mp->d_b_epsilondev_yy_inner_core,b_epsilondev_yy_inner_core,R_size);
      copy_todevice_realw((void**)&mp->d_b_epsilondev_xy_inner_core,b_epsilondev_xy_inner_core,R_size);
      copy_todevice_realw((void**)&mp->d_b_epsilondev_xz_inner_core,b_epsilondev_xz_inner_core,R_size);
      copy_todevice_realw((void**)&mp->d_b_epsilondev_yz_inner_core,b_epsilondev_yz_inner_core,R_size);
      // strain
      copy_todevice_realw((void**)&mp->d_b_eps_trace_over_3_inner_core,b_eps_trace_over_3_inner_core,R_size);
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

// STRAIN simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_absorb_device,
              PREPARE_FIELDS_ABSORB_DEVICE)(long* Mesh_pointer_f,
                                            int* nspec2D_xmin_crust_mantle,int* nspec2D_xmax_crust_mantle,
                                            int* nspec2D_ymin_crust_mantle,int* nspec2D_ymax_crust_mantle,
                                            int* NSPEC2DMAX_XMIN_XMAX_CM,int* NSPEC2DMAX_YMIN_YMAX_CM,
                                            int* nimin_crust_mantle,int* nimax_crust_mantle,
                                            int* njmin_crust_mantle,int* njmax_crust_mantle,
                                            int* nkmin_xi_crust_mantle,int* nkmin_eta_crust_mantle,
                                            int* ibelm_xmin_crust_mantle,int* ibelm_xmax_crust_mantle,
                                            int* ibelm_ymin_crust_mantle,int* ibelm_ymax_crust_mantle,
                                            realw* normal_xmin_crust_mantle,realw* normal_xmax_crust_mantle,
                                            realw* normal_ymin_crust_mantle,realw* normal_ymax_crust_mantle,
                                            realw* jacobian2D_xmin_crust_mantle, realw* jacobian2D_xmax_crust_mantle,
                                            realw* jacobian2D_ymin_crust_mantle, realw* jacobian2D_ymax_crust_mantle,
                                            realw* rho_vp_crust_mantle,
                                            realw* rho_vs_crust_mantle,
                                            int* nspec2D_xmin_outer_core,int* nspec2D_xmax_outer_core,
                                            int* nspec2D_ymin_outer_core,int* nspec2D_ymax_outer_core,
                                            int* nspec2D_zmin_outer_core,
                                            int* NSPEC2DMAX_XMIN_XMAX_OC,int* NSPEC2DMAX_YMIN_YMAX_OC,
                                            int* nimin_outer_core,int* nimax_outer_core,
                                            int* njmin_outer_core,int* njmax_outer_core,
                                            int* nkmin_xi_outer_core,int* nkmin_eta_outer_core,
                                            int* ibelm_xmin_outer_core,int* ibelm_xmax_outer_core,
                                            int* ibelm_ymin_outer_core,int* ibelm_ymax_outer_core,
                                            realw* jacobian2D_xmin_outer_core, realw* jacobian2D_xmax_outer_core,
                                            realw* jacobian2D_ymin_outer_core, realw* jacobian2D_ymax_outer_core,
                                            realw* vp_outer_core) {

  TRACE("prepare_fields_absorb_device");
  int size;

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // checks flag
  if( ! mp->absorbing_conditions ){ exit_on_cuda_error("prepare_fields_absorb_device absorbing_conditions not properly initialized"); }

  // crust_mantle
  mp->nspec2D_xmin_crust_mantle = *nspec2D_xmin_crust_mantle;
  mp->nspec2D_xmax_crust_mantle = *nspec2D_xmax_crust_mantle;
  mp->nspec2D_ymin_crust_mantle = *nspec2D_ymin_crust_mantle;
  mp->nspec2D_ymax_crust_mantle = *nspec2D_ymax_crust_mantle;

  // vp & vs
  size = NGLL3*(mp->NSPEC_CRUST_MANTLE);
  copy_todevice_realw((void**)&mp->d_rho_vp_crust_mantle,rho_vp_crust_mantle,size);
  copy_todevice_realw((void**)&mp->d_rho_vs_crust_mantle,rho_vs_crust_mantle,size);

  // ijk index arrays
  copy_todevice_int((void**)&mp->d_nkmin_xi_crust_mantle,nkmin_xi_crust_mantle,2*(*NSPEC2DMAX_XMIN_XMAX_CM));
  copy_todevice_int((void**)&mp->d_nkmin_eta_crust_mantle,nkmin_eta_crust_mantle,2*(*NSPEC2DMAX_YMIN_YMAX_CM));
  copy_todevice_int((void**)&mp->d_njmin_crust_mantle,njmin_crust_mantle,2*(*NSPEC2DMAX_XMIN_XMAX_CM));
  copy_todevice_int((void**)&mp->d_njmax_crust_mantle,njmax_crust_mantle,2*(*NSPEC2DMAX_XMIN_XMAX_CM));
  copy_todevice_int((void**)&mp->d_nimin_crust_mantle,nimin_crust_mantle,2*(*NSPEC2DMAX_YMIN_YMAX_CM));
  copy_todevice_int((void**)&mp->d_nimax_crust_mantle,nimax_crust_mantle,2*(*NSPEC2DMAX_YMIN_YMAX_CM));

  // xmin
  if( mp->nspec2D_xmin_crust_mantle > 0 ){
    copy_todevice_int((void**)&mp->d_ibelm_xmin_crust_mantle,ibelm_xmin_crust_mantle,mp->nspec2D_xmin_crust_mantle);
    copy_todevice_realw((void**)&mp->d_normal_xmin_crust_mantle,normal_xmin_crust_mantle,
                        NDIM*NGLL2*(mp->nspec2D_xmin_crust_mantle));
    copy_todevice_realw((void**)&mp->d_jacobian2D_xmin_crust_mantle,jacobian2D_xmin_crust_mantle,
                        NGLL2*(mp->nspec2D_xmin_crust_mantle));
    // boundary buffer
    if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
      print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_absorb_xmin_crust_mantle,
                              NDIM*NGLL2*(mp->nspec2D_xmin_crust_mantle)*sizeof(realw)),1202);
    }
  }

  // xmax
  if( mp->nspec2D_xmax_crust_mantle > 0 ){
    copy_todevice_int((void**)&mp->d_ibelm_xmax_crust_mantle,ibelm_xmax_crust_mantle,mp->nspec2D_xmax_crust_mantle);
    copy_todevice_realw((void**)&mp->d_normal_xmax_crust_mantle,normal_xmax_crust_mantle,
                        NDIM*NGLL2*(mp->nspec2D_xmax_crust_mantle));
    copy_todevice_realw((void**)&mp->d_jacobian2D_xmax_crust_mantle,jacobian2D_xmax_crust_mantle,
                        NGLL2*(mp->nspec2D_xmax_crust_mantle));
    // boundary buffer
    if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
      print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_absorb_xmax_crust_mantle,
                                         NDIM*NGLL2*(mp->nspec2D_xmax_crust_mantle)*sizeof(realw)),1202);
    }
  }

  // ymin
  if( mp->nspec2D_ymin_crust_mantle > 0 ){
    copy_todevice_int((void**)&mp->d_ibelm_ymin_crust_mantle,ibelm_ymin_crust_mantle,mp->nspec2D_ymin_crust_mantle);
    copy_todevice_realw((void**)&mp->d_normal_ymin_crust_mantle,normal_ymin_crust_mantle,
                        NDIM*NGLL2*(mp->nspec2D_ymin_crust_mantle));
    copy_todevice_realw((void**)&mp->d_jacobian2D_ymin_crust_mantle,jacobian2D_ymin_crust_mantle,
                        NGLL2*(mp->nspec2D_ymin_crust_mantle));
    // boundary buffer
    if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
      print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_absorb_ymin_crust_mantle,
                                         NDIM*NGLL2*(mp->nspec2D_ymin_crust_mantle)*sizeof(realw)),1202);
    }
  }

  // ymax
  if( mp->nspec2D_ymax_crust_mantle > 0 ){
    copy_todevice_int((void**)&mp->d_ibelm_ymax_crust_mantle,ibelm_ymax_crust_mantle,mp->nspec2D_ymax_crust_mantle);
    copy_todevice_realw((void**)&mp->d_normal_ymax_crust_mantle,normal_ymax_crust_mantle,
                        NDIM*NGLL2*(mp->nspec2D_ymax_crust_mantle));
    copy_todevice_realw((void**)&mp->d_jacobian2D_ymax_crust_mantle,jacobian2D_ymax_crust_mantle,
                        NGLL2*(mp->nspec2D_ymax_crust_mantle));
    // boundary buffer
    if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
      print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_absorb_ymax_crust_mantle,
                                         NDIM*NGLL2*(mp->nspec2D_ymax_crust_mantle)*sizeof(realw)),1202);
    }
  }


  // outer_core
  mp->nspec2D_xmin_outer_core = *nspec2D_xmin_outer_core;
  mp->nspec2D_xmax_outer_core = *nspec2D_xmax_outer_core;
  mp->nspec2D_ymin_outer_core = *nspec2D_ymin_outer_core;
  mp->nspec2D_ymax_outer_core = *nspec2D_ymax_outer_core;
  mp->nspec2D_zmin_outer_core = *nspec2D_zmin_outer_core;

  // vp
  size = NGLL3*(mp->NSPEC_OUTER_CORE);
  copy_todevice_realw((void**)&mp->d_vp_outer_core,vp_outer_core,size);

  // ijk index arrays
  copy_todevice_int((void**)&mp->d_nkmin_xi_outer_core,nkmin_xi_outer_core,2*(*NSPEC2DMAX_XMIN_XMAX_OC));
  copy_todevice_int((void**)&mp->d_nkmin_eta_outer_core,nkmin_eta_outer_core,2*(*NSPEC2DMAX_YMIN_YMAX_OC));
  copy_todevice_int((void**)&mp->d_njmin_outer_core,njmin_outer_core,2*(*NSPEC2DMAX_XMIN_XMAX_OC));
  copy_todevice_int((void**)&mp->d_njmax_outer_core,njmax_outer_core,2*(*NSPEC2DMAX_XMIN_XMAX_OC));
  copy_todevice_int((void**)&mp->d_nimin_outer_core,nimin_outer_core,2*(*NSPEC2DMAX_YMIN_YMAX_OC));
  copy_todevice_int((void**)&mp->d_nimax_outer_core,nimax_outer_core,2*(*NSPEC2DMAX_YMIN_YMAX_OC));

  // xmin
  if( mp->nspec2D_xmin_outer_core > 0 ){
    copy_todevice_int((void**)&mp->d_ibelm_xmin_outer_core,ibelm_xmin_outer_core,mp->nspec2D_xmin_outer_core);
    copy_todevice_realw((void**)&mp->d_jacobian2D_xmin_outer_core,jacobian2D_xmin_outer_core,
                        NGLL2*(mp->nspec2D_xmin_outer_core));
    // boundary buffer
    if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
      print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_absorb_xmin_outer_core,
                                         NGLL2*(mp->nspec2D_xmin_outer_core)*sizeof(realw)),1202);
    }
  }

  // xmax
  if( mp->nspec2D_xmax_outer_core > 0 ){
    copy_todevice_int((void**)&mp->d_ibelm_xmax_outer_core,ibelm_xmax_outer_core,mp->nspec2D_xmax_outer_core);
    copy_todevice_realw((void**)&mp->d_jacobian2D_xmax_outer_core,jacobian2D_xmax_outer_core,
                        NGLL2*(mp->nspec2D_xmax_outer_core));
    // boundary buffer
    if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
      print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_absorb_xmax_outer_core,
                                         NGLL2*(mp->nspec2D_xmax_outer_core)*sizeof(realw)),1202);
    }
  }

  // ymin
  if( mp->nspec2D_ymin_outer_core > 0 ){
    copy_todevice_int((void**)&mp->d_ibelm_ymin_outer_core,ibelm_ymin_outer_core,mp->nspec2D_ymin_outer_core);
    copy_todevice_realw((void**)&mp->d_jacobian2D_ymin_outer_core,jacobian2D_ymin_outer_core,
                        NGLL2*(mp->nspec2D_ymin_outer_core));
    // boundary buffer
    if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
      print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_absorb_ymin_outer_core,
                                         NGLL2*(mp->nspec2D_ymin_outer_core)*sizeof(realw)),1202);
    }
  }

  // ymax
  if( mp->nspec2D_ymax_outer_core > 0 ){
    copy_todevice_int((void**)&mp->d_ibelm_ymax_outer_core,ibelm_ymax_outer_core,mp->nspec2D_ymax_outer_core);
    copy_todevice_realw((void**)&mp->d_jacobian2D_ymax_outer_core,jacobian2D_ymax_outer_core,
                        NGLL2*(mp->nspec2D_ymax_outer_core));
    // boundary buffer
    if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
      print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_absorb_ymax_outer_core,
                                         NGLL2*(mp->nspec2D_ymax_outer_core)*sizeof(realw)),1202);
    }
  }

  // zmin
  if( mp->nspec2D_zmin_outer_core > 0 ){
    // note: ibelm_bottom_outer_core and jacobian2D_bottom_outer_core will be allocated
    //          when preparing the outer core
    // boundary buffer
    if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
      print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_absorb_zmin_outer_core,
                                         NGLL2*(mp->nspec2D_zmin_outer_core)*sizeof(realw)),1202);
    }
  }

}

/* ----------------------------------------------------------------------------------------------- */

// MPI interfaces

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_mpi_buffers_device,
              PREPARE_MPI_BUFFERS_DEVICE)(long* Mesh_pointer_f,
                                          int* num_interfaces_crust_mantle,
                                          int* max_nibool_interfaces_cm,
                                          int* nibool_interfaces_crust_mantle,
                                          int* ibool_interfaces_crust_mantle,
                                          int* num_interfaces_inner_core,
                                          int* max_nibool_interfaces_ic,
                                          int* nibool_interfaces_inner_core,
                                          int* ibool_interfaces_inner_core,
                                          int* num_interfaces_outer_core,
                                          int* max_nibool_interfaces_oc,
                                          int* nibool_interfaces_outer_core,
                                          int* ibool_interfaces_outer_core){

  TRACE("prepare_mpi_buffers_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // prepares interprocess-edge exchange information

  // crust/mantle mesh
  mp->num_interfaces_crust_mantle = *num_interfaces_crust_mantle;
  mp->max_nibool_interfaces_cm = *max_nibool_interfaces_cm;
  if( mp->num_interfaces_crust_mantle > 0 ){
    // number of ibool entries array
    copy_todevice_int((void**)&mp->d_nibool_interfaces_crust_mantle,nibool_interfaces_crust_mantle,
                      mp->num_interfaces_crust_mantle);
    // ibool entries (iglob indices) values on interface
    copy_todevice_int((void**)&mp->d_ibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle,
                      (mp->num_interfaces_crust_mantle)*(mp->max_nibool_interfaces_cm));
    // allocates mpi buffer for exchange with cpu
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_send_accel_buffer_crust_mantle),
                                       NDIM*(mp->max_nibool_interfaces_cm)*(mp->num_interfaces_crust_mantle)*sizeof(realw)),4004);
    if( mp->simulation_type == 3){
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_send_accel_buffer_crust_mantle),
                                        NDIM*(mp->max_nibool_interfaces_cm)*(mp->num_interfaces_crust_mantle)*sizeof(realw)),4004);
    }

  }

  // inner core mesh
  mp->num_interfaces_inner_core = *num_interfaces_inner_core;
  mp->max_nibool_interfaces_ic = *max_nibool_interfaces_ic;
  if( mp->num_interfaces_inner_core > 0 ){
    // number of ibool entries array
    copy_todevice_int((void**)&mp->d_nibool_interfaces_inner_core,nibool_interfaces_inner_core,
                      mp->num_interfaces_inner_core);
    // ibool entries (iglob indices) values on interface
    copy_todevice_int((void**)&mp->d_ibool_interfaces_inner_core,ibool_interfaces_inner_core,
                      (mp->num_interfaces_inner_core)*(mp->max_nibool_interfaces_ic));
    // allocates mpi buffer for exchange with cpu
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_send_accel_buffer_inner_core),
                                       NDIM*(mp->max_nibool_interfaces_ic)*(mp->num_interfaces_inner_core)*sizeof(realw)),4004);
    if( mp->simulation_type == 3){
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_send_accel_buffer_inner_core),
                                        NDIM*(mp->max_nibool_interfaces_ic)*(mp->num_interfaces_inner_core)*sizeof(realw)),4004);
    }

  }

  // outer core mesh
  // note: uses only scalar wavefield arrays
  mp->num_interfaces_outer_core = *num_interfaces_outer_core;
  mp->max_nibool_interfaces_oc = *max_nibool_interfaces_oc;
  if( mp->num_interfaces_outer_core > 0 ){
    // number of ibool entries array
    copy_todevice_int((void**)&mp->d_nibool_interfaces_outer_core,nibool_interfaces_outer_core,
                      mp->num_interfaces_outer_core);
    // ibool entries (iglob indices) values on interface
    copy_todevice_int((void**)&mp->d_ibool_interfaces_outer_core,ibool_interfaces_outer_core,
                      (mp->num_interfaces_outer_core)*(mp->max_nibool_interfaces_oc));
    // allocates mpi buffer for exchange with cpu
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_send_accel_buffer_outer_core),
                                       (mp->max_nibool_interfaces_oc)*(mp->num_interfaces_outer_core)*sizeof(realw)),4004);
    if( mp->simulation_type == 3){
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_send_accel_buffer_outer_core),
                                        (mp->max_nibool_interfaces_oc)*(mp->num_interfaces_outer_core)*sizeof(realw)),4004);
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

// for NOISE simulations

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(prepare_fields_noise_device,
              PREPARE_FIELDS_NOISE_DEVICE)(long* Mesh_pointer_f,
                                           int* NSPEC_TOP,
                                           int* NSTEP,
                                           int* h_ibelm_top_crust_mantle,
                                           realw* noise_sourcearray,
                                           realw* normal_x_noise,
                                           realw* normal_y_noise,
                                           realw* normal_z_noise,
                                           realw* mask_noise,
                                           realw* jacobian2D_top_crust_mantle) {

  TRACE("prepare_fields_noise_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // free surface
  mp->nspec2D_top_crust_mantle = *NSPEC_TOP;
  if( mp->nspec2D_top_crust_mantle > 0 ){
    // note: d_ibelm_top_crust_mantle will only be needed for noise computations
    copy_todevice_int((void**)&mp->d_ibelm_top_crust_mantle,h_ibelm_top_crust_mantle,mp->nspec2D_top_crust_mantle);
    // alloc storage for the surface buffer to be copied
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_noise_surface_movie,
                                       NDIM*NGLL2*(mp->nspec2D_top_crust_mantle)*sizeof(realw)),7005);
  }else{
    // for global mesh: each crust/mantle slice should have at top a free surface
    exit_on_cuda_error("prepare_fields_noise_device NSPEC_TOP not properly initialized");
  }


  // prepares noise source array
  if( mp->noise_tomography == 1 ){
    copy_todevice_realw((void**)&mp->d_noise_sourcearray,noise_sourcearray,NDIM*NGLL3*(*NSTEP));
  }

  // prepares noise directions
  if( mp->noise_tomography > 1 ){
    int nface_size = NGLL2*(mp->nspec2D_top_crust_mantle);
    // allocates memory on GPU
    copy_todevice_realw((void**)&mp->d_normal_x_noise,normal_x_noise,nface_size);
    copy_todevice_realw((void**)&mp->d_normal_y_noise,normal_y_noise,nface_size);
    copy_todevice_realw((void**)&mp->d_normal_z_noise,normal_z_noise,nface_size);

    copy_todevice_realw((void**)&mp->d_mask_noise,mask_noise,nface_size);
    copy_todevice_realw((void**)&mp->d_jacobian2D_top_crust_mantle,jacobian2D_top_crust_mantle,nface_size);
  }

  // prepares noise strength kernel
  if( mp->noise_tomography == 3 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_Sigma_kl),
                                       NGLL3*(mp->NSPEC_CRUST_MANTLE)*sizeof(realw)),7401);
    // initializes kernel values to zero
    print_CUDA_error_if_any(cudaMemset(mp->d_Sigma_kl,0,
                                       NGLL3*mp->NSPEC_CRUST_MANTLE*sizeof(realw)),7403);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_fields_noise_device");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// OCEANS

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_oceans_device,
              PREPARE_OCEANS_DEVICE)(long* Mesh_pointer_f,
                                     int* npoin_oceans,
                                     int* h_iglob_ocean_load,
                                     realw* h_rmass_ocean_load_selected,
                                     realw* h_normal_ocean_load) {

  TRACE("prepare_oceans_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // arrays with global points on ocean surface
  mp->npoin_oceans = *npoin_oceans;

  // checks for global partitions, each slice must have a top surface with points on it
  if( mp->npoin_oceans == 0 ){ exit_on_cuda_error("prepare_oceans_device has zero npoin_oceans"); }

  // global point indices
  copy_todevice_int((void**)&mp->d_ibool_ocean_load,h_iglob_ocean_load,mp->npoin_oceans);

  // mass matrix
  copy_todevice_realw((void**)&mp->d_rmass_ocean_load,h_rmass_ocean_load_selected,mp->npoin_oceans);

  // normals
  copy_todevice_realw((void**)&mp->d_normal_ocean_load,h_normal_ocean_load,NDIM*mp->npoin_oceans);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_oceans_device");
#endif
}



/* ----------------------------------------------------------------------------------------------- */

// Earth regions

// CRUST / MANTLE

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_crust_mantle_device,
             PREPARE_CRUST_MANTLE_DEVICE)(long* Mesh_pointer_f,
                                          realw* h_xix, realw* h_xiy, realw* h_xiz,
                                          realw* h_etax, realw* h_etay, realw* h_etaz,
                                          realw* h_gammax, realw* h_gammay, realw* h_gammaz,
                                          realw* h_rho,
                                          realw* h_kappav, realw* h_muv,
                                          realw* h_kappah, realw* h_muh,
                                          realw* h_eta_aniso,
                                          realw* h_rmassx,realw* h_rmassy,realw* h_rmassz,
                                          realw* h_b_rmassx,realw* h_b_rmassy,
                                          int* h_ibool,
                                          realw* h_xstore, realw* h_ystore, realw* h_zstore,
                                          int* h_ispec_is_tiso,
                                          realw *c11store,realw *c12store,realw *c13store,
                                          realw *c14store,realw *c15store,realw *c16store,
                                          realw *c22store,realw *c23store,realw *c24store,
                                          realw *c25store,realw *c26store,realw *c33store,
                                          realw *c34store,realw *c35store,realw *c36store,
                                          realw *c44store,realw *c45store,realw *c46store,
                                          realw *c55store,realw *c56store,realw *c66store,
                                          int* num_phase_ispec,
                                          int* phase_ispec_inner,
                                          int* nspec_outer,
                                          int* nspec_inner,
                                          int* NSPEC2D_BOTTOM_CM,
                                          int* h_ibelm_bottom_crust_mantle,
                                          int* NCHUNKS_VAL,
                                          int* num_colors_outer,
                                          int* num_colors_inner,
                                          int* num_elem_colors) {

  TRACE("prepare_crust_mantle_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  /* Assuming NGLLX=5. Padded is then 128 (5^3+3) */
  int size_padded = NGLL3_PADDED * (mp->NSPEC_CRUST_MANTLE);
  int size_glob = mp->NGLOB_CRUST_MANTLE;

  // mesh
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xix_crust_mantle, size_padded*sizeof(realw)),1001);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xiy_crust_mantle, size_padded*sizeof(realw)),1002);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xiz_crust_mantle, size_padded*sizeof(realw)),1003);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etax_crust_mantle, size_padded*sizeof(realw)),1004);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etay_crust_mantle, size_padded*sizeof(realw)),1005);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etaz_crust_mantle, size_padded*sizeof(realw)),1006);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammax_crust_mantle, size_padded*sizeof(realw)),1007);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammay_crust_mantle, size_padded*sizeof(realw)),1008);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammaz_crust_mantle, size_padded*sizeof(realw)),1009);

  // transfer constant element data with padding
  for(int i=0;i < mp->NSPEC_CRUST_MANTLE;i++) {
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xix_crust_mantle + i*NGLL3_PADDED, &h_xix[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1501);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xiy_crust_mantle+i*NGLL3_PADDED,   &h_xiy[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1502);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xiz_crust_mantle+i*NGLL3_PADDED,   &h_xiz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1503);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etax_crust_mantle+i*NGLL3_PADDED,  &h_etax[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1504);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etay_crust_mantle+i*NGLL3_PADDED,  &h_etay[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1505);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etaz_crust_mantle+i*NGLL3_PADDED,  &h_etaz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1506);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammax_crust_mantle+i*NGLL3_PADDED,&h_gammax[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1507);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammay_crust_mantle+i*NGLL3_PADDED,&h_gammay[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1508);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammaz_crust_mantle+i*NGLL3_PADDED,&h_gammaz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1509);
  }

  // global indexing
  copy_todevice_int((void**)&mp->d_ibool_crust_mantle,h_ibool,NGLL3*(mp->NSPEC_CRUST_MANTLE));

//  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibool_crust_mantle, size_padded*sizeof(int)),1021);
//  print_CUDA_error_if_any(cudaMemcpy(mp->d_ibool_crust_mantle, h_ibool,
//                                     NGLL3*(mp->NSPEC_CRUST_MANTLE)*sizeof(int),cudaMemcpyHostToDevice),1022);

  // transverse isotropic elements
  // only needed if not anisotropic 3D mantle
  if( ! mp->anisotropic_3D_mantle ){
    // no anisotropy

    // transverse isotropy flag
    copy_todevice_int((void**)&mp->d_ispec_is_tiso_crust_mantle,h_ispec_is_tiso,mp->NSPEC_CRUST_MANTLE);

    // kappavstore, kappahstore
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_kappavstore_crust_mantle, size_padded*sizeof(realw)),1010);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_kappahstore_crust_mantle, size_padded*sizeof(realw)),1011);

    // muvstore,muhstore
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_muvstore_crust_mantle, size_padded*sizeof(realw)),1012);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_muhstore_crust_mantle, size_padded*sizeof(realw)),1013);

    // eta_anisostore
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_eta_anisostore_crust_mantle, size_padded*sizeof(realw)),1014);

    // transfer with padding
    for(int i=0;i < mp->NSPEC_CRUST_MANTLE;i++) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_kappavstore_crust_mantle+i*NGLL3_PADDED,&h_kappav[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1510);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_kappahstore_crust_mantle+i*NGLL3_PADDED,&h_kappah[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1511);

      print_CUDA_error_if_any(cudaMemcpy(mp->d_muvstore_crust_mantle+i*NGLL3_PADDED,   &h_muv[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1512);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_muhstore_crust_mantle+i*NGLL3_PADDED,&h_muh[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1513);

      print_CUDA_error_if_any(cudaMemcpy(mp->d_eta_anisostore_crust_mantle+i*NGLL3_PADDED,&h_eta_aniso[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1514);
    }
  }else{
    // anisotropic 3D mantle

    // allocates memory on GPU
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c11store_crust_mantle),
                                       size_padded*sizeof(realw)),4700);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c12store_crust_mantle),
                                       size_padded*sizeof(realw)),4701);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c13store_crust_mantle),
                                       size_padded*sizeof(realw)),4702);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c14store_crust_mantle),
                                       size_padded*sizeof(realw)),4703);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c15store_crust_mantle),
                                       size_padded*sizeof(realw)),4704);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c16store_crust_mantle),
                                       size_padded*sizeof(realw)),4705);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c22store_crust_mantle),
                                       size_padded*sizeof(realw)),4706);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c23store_crust_mantle),
                                       size_padded*sizeof(realw)),4707);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c24store_crust_mantle),
                                       size_padded*sizeof(realw)),4708);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c25store_crust_mantle),
                                       size_padded*sizeof(realw)),4709);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c26store_crust_mantle),
                                       size_padded*sizeof(realw)),4710);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c33store_crust_mantle),
                                       size_padded*sizeof(realw)),4711);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c34store_crust_mantle),
                                       size_padded*sizeof(realw)),4712);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c35store_crust_mantle),
                                       size_padded*sizeof(realw)),4713);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c36store_crust_mantle),
                                       size_padded*sizeof(realw)),4714);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c44store_crust_mantle),
                                       size_padded*sizeof(realw)),4715);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c45store_crust_mantle),
                                       size_padded*sizeof(realw)),4716);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c46store_crust_mantle),
                                       size_padded*sizeof(realw)),4717);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c55store_crust_mantle),
                                       size_padded*sizeof(realw)),4718);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c56store_crust_mantle),
                                       size_padded*sizeof(realw)),4719);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c66store_crust_mantle),
                                       size_padded*sizeof(realw)),4720);

    // transfer constant element data with padding
    for(int i=0;i < mp->NSPEC_CRUST_MANTLE;i++) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c11store_crust_mantle + i*NGLL3_PADDED, &c11store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4800);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c12store_crust_mantle + i*NGLL3_PADDED, &c12store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4801);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c13store_crust_mantle + i*NGLL3_PADDED, &c13store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4802);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c14store_crust_mantle + i*NGLL3_PADDED, &c14store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4803);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c15store_crust_mantle + i*NGLL3_PADDED, &c15store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4804);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c16store_crust_mantle + i*NGLL3_PADDED, &c16store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4805);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c22store_crust_mantle + i*NGLL3_PADDED, &c22store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4806);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c23store_crust_mantle + i*NGLL3_PADDED, &c23store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4807);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c24store_crust_mantle + i*NGLL3_PADDED, &c24store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4808);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c25store_crust_mantle + i*NGLL3_PADDED, &c25store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4809);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c26store_crust_mantle + i*NGLL3_PADDED, &c26store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4810);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c33store_crust_mantle + i*NGLL3_PADDED, &c33store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4811);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c34store_crust_mantle + i*NGLL3_PADDED, &c34store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4812);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c35store_crust_mantle + i*NGLL3_PADDED, &c35store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4813);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c36store_crust_mantle + i*NGLL3_PADDED, &c36store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4814);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c44store_crust_mantle + i*NGLL3_PADDED, &c44store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4815);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c45store_crust_mantle + i*NGLL3_PADDED, &c45store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4816);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c46store_crust_mantle + i*NGLL3_PADDED, &c46store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4817);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c55store_crust_mantle + i*NGLL3_PADDED, &c55store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4818);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c56store_crust_mantle + i*NGLL3_PADDED, &c56store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4819);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c66store_crust_mantle + i*NGLL3_PADDED, &c66store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4820);
    }
  }

  // needed for boundary kernel calculations
  if( mp->simulation_type == 3 && mp->save_boundary_mesh ){
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_rhostore_crust_mantle, size_padded*sizeof(realw)),1010);
    for(int i=0;i < mp->NSPEC_CRUST_MANTLE;i++) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_rhostore_crust_mantle+i*NGLL3_PADDED, &h_rho[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),2106);
    }
  }

  // mesh locations
  // ystore & zstore needed for tiso elements
  copy_todevice_realw((void**)&mp->d_ystore_crust_mantle,h_ystore,size_glob);
  copy_todevice_realw((void**)&mp->d_zstore_crust_mantle,h_zstore,size_glob);


  // xstore only needed when gravity is on
  if( mp->gravity ){
    copy_todevice_realw((void**)&mp->d_xstore_crust_mantle,h_xstore,size_glob);
  }

  // inner/outer elements
  mp->num_phase_ispec_crust_mantle = *num_phase_ispec;
  copy_todevice_int((void**)&mp->d_phase_ispec_inner_crust_mantle,phase_ispec_inner,
                    mp->num_phase_ispec_crust_mantle*2);

  mp->nspec_outer_crust_mantle = *nspec_outer;
  mp->nspec_inner_crust_mantle = *nspec_inner;

  // CMB/fluid outer core coupling
  mp->nspec2D_bottom_crust_mantle = *NSPEC2D_BOTTOM_CM;
  copy_todevice_int((void**)&mp->d_ibelm_bottom_crust_mantle,h_ibelm_bottom_crust_mantle,mp->nspec2D_bottom_crust_mantle);

  // wavefield
  int size = NDIM * mp->NGLOB_CRUST_MANTLE;

  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_displ_crust_mantle),sizeof(realw)*size),4001);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_veloc_crust_mantle),sizeof(realw)*size),4002);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_accel_crust_mantle),sizeof(realw)*size),4003);
  // backward/reconstructed wavefield
  if( mp->simulation_type == 3 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_displ_crust_mantle),sizeof(realw)*size),4011);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_veloc_crust_mantle),sizeof(realw)*size),4012);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_accel_crust_mantle),sizeof(realw)*size),4013);
    // debug
    #if DEBUG_BACKWARD_SIMULATIONS == 1
    //debugging with empty arrays
    print_CUDA_error_if_any(cudaMemset(mp->d_b_displ_crust_mantle,0,sizeof(realw)*size),5111);
    print_CUDA_error_if_any(cudaMemset(mp->d_b_veloc_crust_mantle,0,sizeof(realw)*size),5111);
    print_CUDA_error_if_any(cudaMemset(mp->d_b_accel_crust_mantle,0,sizeof(realw)*size),5111);
    #endif
  }

  #ifdef USE_TEXTURES_FIELDS
  {
    #ifdef USE_OLDER_CUDA4_GPU
      const textureReference* d_displ_cm_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&mp->d_displ_cm_tex_ref_ptr, "d_displ_cm_tex"), 4021);
      cudaChannelFormatDesc channelDesc1 = cudaCreateChannelDesc<realw>();
      print_CUDA_error_if_any(cudaBindTexture(0, mp->d_displ_cm_tex_ref_ptr, mp->d_displ_crust_mantle,
                                              &channelDesc1, sizeof(realw)*size), 4021);

      const textureReference* d_accel_cm_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_accel_cm_tex_ref_ptr, "d_accel_cm_tex"), 4023);
      cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<realw>();
      print_CUDA_error_if_any(cudaBindTexture(0, d_accel_cm_tex_ref_ptr, mp->d_accel_crust_mantle,
                                              &channelDesc2, sizeof(realw)*size), 4023);
    #else
      cudaChannelFormatDesc channelDesc1 = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_displ_cm_tex, mp->d_displ_crust_mantle,
                                              &channelDesc1, sizeof(realw)*size), 4021);

      cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_accel_cm_tex, mp->d_accel_crust_mantle,
                                              &channelDesc2, sizeof(realw)*size), 4023);
    #endif
  }
  #endif


  // mass matrices
  copy_todevice_realw((void**)&mp->d_rmassz_crust_mantle,h_rmassz,size_glob);
  if( (*NCHUNKS_VAL != 6 && mp->absorbing_conditions) || ( mp->rotation && mp->exact_mass_matrix_for_rotation)){
    copy_todevice_realw((void**)&mp->d_rmassx_crust_mantle,h_rmassx,size_glob);
    copy_todevice_realw((void**)&mp->d_rmassy_crust_mantle,h_rmassy,size_glob);
  }else{
    mp->d_rmassx_crust_mantle = mp->d_rmassz_crust_mantle;
    mp->d_rmassy_crust_mantle = mp->d_rmassz_crust_mantle;
  }

  // kernel simulations
  if( mp->simulation_type == 3 ){
    mp->d_b_rmassz_crust_mantle = mp->d_rmassz_crust_mantle;
    if( mp->rotation && mp->exact_mass_matrix_for_rotation ){
      copy_todevice_realw((void**)&mp->d_b_rmassx_crust_mantle,h_b_rmassx,size_glob);
      copy_todevice_realw((void**)&mp->d_b_rmassy_crust_mantle,h_b_rmassy,size_glob);
    }else{
      mp->d_b_rmassx_crust_mantle = mp->d_rmassx_crust_mantle;
      mp->d_b_rmassy_crust_mantle = mp->d_rmassy_crust_mantle;
    }

    // kernels
    size = NGLL3*(mp->NSPEC_CRUST_MANTLE);

    // density kernel
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rho_kl_crust_mantle),
                                       size*sizeof(realw)),5204);
    // initializes kernel values to zero
    print_CUDA_error_if_any(cudaMemset(mp->d_rho_kl_crust_mantle,0,size*sizeof(realw)),5207);

    if( ! mp->anisotropic_kl){
      // isotropic kernels
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_alpha_kl_crust_mantle),
                                         size*sizeof(realw)),5205);
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_beta_kl_crust_mantle),
                                         size*sizeof(realw)),5206);
      print_CUDA_error_if_any(cudaMemset(mp->d_alpha_kl_crust_mantle,0,size*sizeof(realw)),5208);
      print_CUDA_error_if_any(cudaMemset(mp->d_beta_kl_crust_mantle,0,size*sizeof(realw)),5209);
    }else{
      // anisotropic kernels
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_cijkl_kl_crust_mantle),
                                         21*size*sizeof(realw)),5206);
      print_CUDA_error_if_any(cudaMemset(mp->d_cijkl_kl_crust_mantle,0,21*size*sizeof(realw)),5209);
    }

    // preconditioner
    if( mp->approximate_hess_kl ){
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_hess_kl_crust_mantle),
                                         size*sizeof(realw)),3030);
      print_CUDA_error_if_any(cudaMemset(mp->d_hess_kl_crust_mantle,0,size*sizeof(realw)),3031);
    }
  }

  // mesh coloring
  mp->num_colors_outer_crust_mantle = *num_colors_outer;
  mp->num_colors_inner_crust_mantle = *num_colors_inner;
  mp->h_num_elem_colors_crust_mantle = (int*) num_elem_colors;

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_crust_mantle_device");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// OUTER CORE

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_outer_core_device,
              PREPARE_OUTER_CORE_DEVICE)(long* Mesh_pointer_f,
                                         realw* h_xix, realw* h_xiy, realw* h_xiz,
                                         realw* h_etax, realw* h_etay, realw* h_etaz,
                                         realw* h_gammax, realw* h_gammay, realw* h_gammaz,
                                         realw* h_rho, realw* h_kappav,
                                         realw* h_rmass,
                                         int* h_ibool,
                                         realw* h_xstore, realw* h_ystore, realw* h_zstore,
                                         int* num_phase_ispec,
                                         int* phase_ispec_inner,
                                         int* nspec_outer,
                                         int* nspec_inner,
                                         int* NSPEC2D_TOP_OC,
                                         int* NSPEC2D_BOTTOM_OC,
                                         realw* h_normal_top_outer_core,
                                         realw* h_normal_bottom_outer_core,
                                         realw* h_jacobian2D_top_outer_core,
                                         realw* h_jacobian2D_bottom_outer_core,
                                         int* h_ibelm_top_outer_core,
                                         int* h_ibelm_bottom_outer_core,
                                         int* num_colors_outer,
                                         int* num_colors_inner,
                                         int* num_elem_colors) {

  TRACE("prepare_outer_core_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  /* Assuming NGLLX=5. Padded is then 128 (5^3+3) */
  int size_padded = NGLL3_PADDED * (mp->NSPEC_OUTER_CORE);
  int size_glob = mp->NGLOB_OUTER_CORE;

  // mesh
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xix_outer_core, size_padded*sizeof(realw)),1101);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xiy_outer_core, size_padded*sizeof(realw)),1102);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xiz_outer_core, size_padded*sizeof(realw)),1103);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etax_outer_core, size_padded*sizeof(realw)),1104);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etay_outer_core, size_padded*sizeof(realw)),1105);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etaz_outer_core, size_padded*sizeof(realw)),1106);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammax_outer_core, size_padded*sizeof(realw)),1107);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammay_outer_core, size_padded*sizeof(realw)),1108);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammaz_outer_core, size_padded*sizeof(realw)),1109);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_kappavstore_outer_core, size_padded*sizeof(realw)),1110);

  // transfer constant element data with padding
  for(int i=0;i < mp->NSPEC_OUTER_CORE;i++) {
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xix_outer_core + i*NGLL3_PADDED, &h_xix[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1501);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xiy_outer_core+i*NGLL3_PADDED,   &h_xiy[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1502);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xiz_outer_core+i*NGLL3_PADDED,   &h_xiz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1503);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etax_outer_core+i*NGLL3_PADDED,  &h_etax[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1504);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etay_outer_core+i*NGLL3_PADDED,  &h_etay[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1505);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etaz_outer_core+i*NGLL3_PADDED,  &h_etaz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1506);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammax_outer_core+i*NGLL3_PADDED,&h_gammax[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1507);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammay_outer_core+i*NGLL3_PADDED,&h_gammay[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1508);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammaz_outer_core+i*NGLL3_PADDED,&h_gammaz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1509);

    print_CUDA_error_if_any(cudaMemcpy(mp->d_kappavstore_outer_core+i*NGLL3_PADDED,&h_kappav[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1510);
  }

  // needed for kernel calculations
  if( mp->simulation_type == 3 ){
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_rhostore_outer_core, size_padded*sizeof(realw)),1010);
    for(int i=0;i < mp->NSPEC_OUTER_CORE;i++) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_rhostore_outer_core+i*NGLL3_PADDED, &h_rho[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),2106);
    }
  }

  // global indexing
  copy_todevice_int((void**)&mp->d_ibool_outer_core,h_ibool,NGLL3*(mp->NSPEC_OUTER_CORE));

//  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibool_outer_core, size_padded*sizeof(int)),1021);
//  print_CUDA_error_if_any(cudaMemcpy(mp->d_ibool_outer_core, h_ibool,
//                                     NGLL3*(mp->NSPEC_OUTER_CORE)*sizeof(int),cudaMemcpyHostToDevice),1022);

  // mesh locations
  // always needed
  copy_todevice_realw((void**)&mp->d_xstore_outer_core,h_xstore,size_glob);
  copy_todevice_realw((void**)&mp->d_ystore_outer_core,h_ystore,size_glob);
  copy_todevice_realw((void**)&mp->d_zstore_outer_core,h_zstore,size_glob);

  // inner/outer elements
  mp->num_phase_ispec_outer_core = *num_phase_ispec;
  copy_todevice_int((void**)&mp->d_phase_ispec_inner_outer_core,phase_ispec_inner,mp->num_phase_ispec_outer_core*2);

  mp->nspec_outer_outer_core = *nspec_outer;
  mp->nspec_inner_outer_core = *nspec_inner;

  // CMB/ICB coupling
  mp->nspec2D_top_outer_core = *NSPEC2D_TOP_OC;
  mp->nspec2D_bottom_outer_core = *NSPEC2D_BOTTOM_OC;

  copy_todevice_int((void**)&mp->d_ibelm_top_outer_core,h_ibelm_top_outer_core,mp->nspec2D_top_outer_core);
  int size_toc = NGLL2*(mp->nspec2D_top_outer_core);
  copy_todevice_realw((void**)&mp->d_jacobian2D_top_outer_core,h_jacobian2D_top_outer_core,size_toc);
  copy_todevice_realw((void**)&mp->d_normal_top_outer_core,h_normal_top_outer_core,NDIM*size_toc);

  copy_todevice_int((void**)&mp->d_ibelm_bottom_outer_core,h_ibelm_bottom_outer_core,mp->nspec2D_bottom_outer_core);
  int size_boc = NGLL2*(mp->nspec2D_bottom_outer_core);
  copy_todevice_realw((void**)&mp->d_jacobian2D_bottom_outer_core,h_jacobian2D_bottom_outer_core,size_boc);
  copy_todevice_realw((void**)&mp->d_normal_bottom_outer_core,h_normal_bottom_outer_core,NDIM*size_boc);

  // wavefield
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_displ_outer_core),sizeof(realw)*size_glob),5001);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_veloc_outer_core),sizeof(realw)*size_glob),5002);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_accel_outer_core),sizeof(realw)*size_glob),5003);
  // backward/reconstructed wavefield
  if( mp->simulation_type == 3 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_displ_outer_core),sizeof(realw)*size_glob),5011);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_veloc_outer_core),sizeof(realw)*size_glob),5022);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_accel_outer_core),sizeof(realw)*size_glob),5033);
    // debug
    #if DEBUG_BACKWARD_SIMULATIONS == 1
    //debugging with empty arrays
    print_CUDA_error_if_any(cudaMemset(mp->d_b_displ_outer_core,0,sizeof(realw)*size_glob),5111);
    print_CUDA_error_if_any(cudaMemset(mp->d_b_veloc_outer_core,0,sizeof(realw)*size_glob),5111);
    print_CUDA_error_if_any(cudaMemset(mp->d_b_accel_outer_core,0,sizeof(realw)*size_glob),5111);
    #endif
  }

  #ifdef USE_TEXTURES_FIELDS
  {
    #ifdef USE_OLDER_CUDA4_GPU
      const textureReference* d_displ_oc_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&mp->d_displ_oc_tex_ref_ptr, "d_displ_oc_tex"), 5021);
      cudaChannelFormatDesc channelDesc1 = cudaCreateChannelDesc<realw>();
      print_CUDA_error_if_any(cudaBindTexture(0, mp->d_displ_oc_tex_ref_ptr, mp->d_displ_outer_core,
                                              &channelDesc1, sizeof(realw)*size_glob), 5021);

      const textureReference* d_accel_oc_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_accel_oc_tex_ref_ptr, "d_accel_oc_tex"), 5023);
      cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<realw>();
      print_CUDA_error_if_any(cudaBindTexture(0, d_accel_oc_tex_ref_ptr, mp->d_accel_outer_core,
                                              &channelDesc2, sizeof(realw)*size_glob), 5023);
    #else
      cudaChannelFormatDesc channelDesc1 = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_displ_oc_tex, mp->d_displ_outer_core,
                                              &channelDesc1, sizeof(realw)*size_glob), 5021);

      cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_accel_oc_tex, mp->d_accel_outer_core,
                                              &channelDesc2, sizeof(realw)*size_glob), 5023);
    #endif
  }
  #endif

  // mass matrix
  copy_todevice_realw((void**)&mp->d_rmass_outer_core,h_rmass,size_glob);

  // kernel simulations
  if( mp->simulation_type == 3 ){
    // mass matrix
    mp->d_b_rmass_outer_core = mp->d_rmass_outer_core;

    //kernels
    int size = NGLL3*(mp->NSPEC_OUTER_CORE);

    // density kernel
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rho_kl_outer_core),
                                       size*sizeof(realw)),5204);
    print_CUDA_error_if_any(cudaMemset(mp->d_rho_kl_outer_core,0,size*sizeof(realw)),5207);

    // isotropic kernel
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_alpha_kl_outer_core),
                                        size*sizeof(realw)),5205);
    print_CUDA_error_if_any(cudaMemset(mp->d_alpha_kl_outer_core,0,size*sizeof(realw)),5208);
  }

  // mesh coloring
  mp->num_colors_outer_outer_core = *num_colors_outer;
  mp->num_colors_inner_outer_core = *num_colors_inner;
  mp->h_num_elem_colors_outer_core = (int*) num_elem_colors;

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_outer_core_device");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// INNER CORE

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_inner_core_device,
              PREPARE_INNER_CORE_DEVICE)(long* Mesh_pointer_f,
                                         realw* h_xix, realw* h_xiy, realw* h_xiz,
                                         realw* h_etax, realw* h_etay, realw* h_etaz,
                                         realw* h_gammax, realw* h_gammay, realw* h_gammaz,
                                         realw* h_rho, realw* h_kappav, realw* h_muv,
                                         realw* h_rmassx,realw* h_rmassy,realw* h_rmassz,
                                         realw* h_b_rmassx,realw* h_b_rmassy,
                                         int* h_ibool,
                                         realw* h_xstore, realw* h_ystore, realw* h_zstore,
                                         realw *c11store,realw *c12store,realw *c13store,
                                         realw *c33store,realw *c44store,
                                         int* h_idoubling_inner_core,
                                         int* num_phase_ispec,
                                         int* phase_ispec_inner,
                                         int* nspec_outer,
                                         int* nspec_inner,
                                         int* NSPEC2D_TOP_IC,
                                         int* h_ibelm_top_inner_core,
                                         int* num_colors_outer,
                                         int* num_colors_inner,
                                         int* num_elem_colors) {

  TRACE("prepare_inner_core_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  /* Assuming NGLLX=5. Padded is then 128 (5^3+3) */
  int size_padded = NGLL3_PADDED * (mp->NSPEC_INNER_CORE);
  int size_glob = mp->NGLOB_INNER_CORE;

  // mesh
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xix_inner_core, size_padded*sizeof(realw)),1201);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xiy_inner_core, size_padded*sizeof(realw)),1202);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xiz_inner_core, size_padded*sizeof(realw)),1203);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etax_inner_core, size_padded*sizeof(realw)),1204);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etay_inner_core, size_padded*sizeof(realw)),1205);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etaz_inner_core, size_padded*sizeof(realw)),1206);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammax_inner_core, size_padded*sizeof(realw)),1207);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammay_inner_core, size_padded*sizeof(realw)),1208);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammaz_inner_core, size_padded*sizeof(realw)),1209);

  // muvstore needed for attenuatioin also for anisotropic inner core
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_muvstore_inner_core, size_padded*sizeof(realw)),1211);

  // transfer constant element data with padding
  for(int i=0;i < mp->NSPEC_INNER_CORE;i++) {
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xix_inner_core + i*NGLL3_PADDED, &h_xix[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1501);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xiy_inner_core+i*NGLL3_PADDED,   &h_xiy[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1502);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xiz_inner_core+i*NGLL3_PADDED,   &h_xiz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1503);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etax_inner_core+i*NGLL3_PADDED,  &h_etax[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1504);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etay_inner_core+i*NGLL3_PADDED,  &h_etay[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1505);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etaz_inner_core+i*NGLL3_PADDED,  &h_etaz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1506);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammax_inner_core+i*NGLL3_PADDED,&h_gammax[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1507);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammay_inner_core+i*NGLL3_PADDED,&h_gammay[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1508);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammaz_inner_core+i*NGLL3_PADDED,&h_gammaz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1509);

    print_CUDA_error_if_any(cudaMemcpy(mp->d_muvstore_inner_core+i*NGLL3_PADDED,   &h_muv[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1511);
  }

  // anisotropy
  if( ! mp->anisotropic_inner_core ){
    // no anisotropy (uses kappav and muv in inner core)
    // kappavstore needed
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_kappavstore_inner_core, size_padded*sizeof(realw)),1010);
    for(int i=0;i < mp->NSPEC_INNER_CORE;i++) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_kappavstore_inner_core+i*NGLL3_PADDED,&h_kappav[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1510);
    }
  }else{
    // anisotropic inner core
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c11store_inner_core),
                                     size_padded*sizeof(realw)),4700);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c12store_inner_core),
                                     size_padded*sizeof(realw)),4701);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c13store_inner_core),
                                     size_padded*sizeof(realw)),4702);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c33store_inner_core),
                                     size_padded*sizeof(realw)),4703);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c44store_inner_core),
                                     size_padded*sizeof(realw)),4704);

    // transfer constant element data with padding
    for(int i=0;i < mp->NSPEC_INNER_CORE;i++) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c11store_inner_core + i*NGLL3_PADDED, &c11store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4800);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c12store_inner_core + i*NGLL3_PADDED, &c12store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4801);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c13store_inner_core + i*NGLL3_PADDED, &c13store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4802);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c33store_inner_core + i*NGLL3_PADDED, &c33store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4803);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c44store_inner_core + i*NGLL3_PADDED, &c44store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4804);
    }
  }

  // needed for boundary kernel calculations
  if( mp->simulation_type == 3 && mp->save_boundary_mesh ){
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_rhostore_inner_core, size_padded*sizeof(realw)),1010);
    for(int i=0;i < mp->NSPEC_INNER_CORE;i++) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_rhostore_inner_core+i*NGLL3_PADDED, &h_rho[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),2106);
    }
  }

  // global indexing
  copy_todevice_int((void**)&mp->d_ibool_inner_core,h_ibool,NGLL3*(mp->NSPEC_INNER_CORE));

//  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibool_inner_core, size_padded*sizeof(int)),1021);
//  print_CUDA_error_if_any(cudaMemcpy(mp->d_ibool_inner_core, h_ibool,
//                                     NGLL3*(mp->NSPEC_INNER_CORE)*sizeof(int),cudaMemcpyHostToDevice),1022);

  // fictious element flags
  copy_todevice_int((void**)&mp->d_idoubling_inner_core,h_idoubling_inner_core,mp->NSPEC_INNER_CORE);

  // mesh locations
  // only needed when gravity is on
  if( mp->gravity ){
    copy_todevice_realw((void**)&mp->d_xstore_inner_core,h_xstore,size_glob);
    copy_todevice_realw((void**)&mp->d_ystore_inner_core,h_ystore,size_glob);
    copy_todevice_realw((void**)&mp->d_zstore_inner_core,h_zstore,size_glob);
  }

  // inner/outer elements
  mp->num_phase_ispec_inner_core = *num_phase_ispec;
  copy_todevice_int((void**)&mp->d_phase_ispec_inner_inner_core,phase_ispec_inner,mp->num_phase_ispec_inner_core*2);

  mp->nspec_outer_inner_core = *nspec_outer;
  mp->nspec_inner_inner_core = *nspec_inner;

  // boundary elements on top
  mp->nspec2D_top_inner_core = *NSPEC2D_TOP_IC;
  copy_todevice_int((void**)&mp->d_ibelm_top_inner_core,h_ibelm_top_inner_core,mp->nspec2D_top_inner_core);

  // wavefield
  int size = NDIM * mp->NGLOB_INNER_CORE;
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_displ_inner_core),sizeof(realw)*size),6001);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_veloc_inner_core),sizeof(realw)*size),6002);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_accel_inner_core),sizeof(realw)*size),6003);
  // backward/reconstructed wavefield
  if( mp->simulation_type == 3 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_displ_inner_core),sizeof(realw)*size),6011);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_veloc_inner_core),sizeof(realw)*size),6012);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_accel_inner_core),sizeof(realw)*size),6013);
    // debug
    #if DEBUG_BACKWARD_SIMULATIONS == 1
    // debugging with empty arrays
    print_CUDA_error_if_any(cudaMemset(mp->d_b_displ_inner_core,0,sizeof(realw)*size),5111);
    print_CUDA_error_if_any(cudaMemset(mp->d_b_veloc_inner_core,0,sizeof(realw)*size),5111);
    print_CUDA_error_if_any(cudaMemset(mp->d_b_accel_inner_core,0,sizeof(realw)*size),5111);
    #endif
  }

  #ifdef USE_TEXTURES_FIELDS
  {
    #ifdef USE_OLDER_CUDA4_GPU
      const textureReference* d_displ_ic_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&mp->d_displ_ic_tex_ref_ptr, "d_displ_ic_tex"), 6021);
      cudaChannelFormatDesc channelDesc1 = cudaCreateChannelDesc<realw>();
      print_CUDA_error_if_any(cudaBindTexture(0, mp->d_displ_ic_tex_ref_ptr, mp->d_displ_inner_core,
                                              &channelDesc1, sizeof(realw)*size), 6021);

      const textureReference* d_accel_ic_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_accel_ic_tex_ref_ptr, "d_accel_ic_tex"), 6023);
      cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<realw>();
      print_CUDA_error_if_any(cudaBindTexture(0, d_accel_ic_tex_ref_ptr, mp->d_accel_inner_core,
                                              &channelDesc2, sizeof(realw)*size), 6023);
    #else
      cudaChannelFormatDesc channelDesc1 = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_displ_ic_tex, mp->d_displ_inner_core,
                                              &channelDesc1, sizeof(realw)*size), 6021);

      cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_accel_ic_tex, mp->d_accel_inner_core,
                                              &channelDesc2, sizeof(realw)*size), 6023);
    #endif


  }
  #endif

  // mass matrix
  copy_todevice_realw((void**)&mp->d_rmassz_inner_core,h_rmassz,size_glob);
  if( mp->rotation && mp->exact_mass_matrix_for_rotation ){
    copy_todevice_realw((void**)&mp->d_rmassx_inner_core,h_rmassx,size_glob);
    copy_todevice_realw((void**)&mp->d_rmassy_inner_core,h_rmassy,size_glob);
  }else{
    mp->d_rmassx_inner_core = mp->d_rmassz_inner_core;
    mp->d_rmassy_inner_core = mp->d_rmassz_inner_core;
  }

  // kernel simulations
  if( mp->simulation_type == 3 ){
    // mass matrices
    mp->d_b_rmassz_inner_core = mp->d_rmassz_inner_core;
    if( mp->rotation && mp->exact_mass_matrix_for_rotation ){
      copy_todevice_realw((void**)&mp->d_b_rmassx_inner_core,h_b_rmassx,size_glob);
      copy_todevice_realw((void**)&mp->d_b_rmassy_inner_core,h_b_rmassy,size_glob);
    }else{
      mp->d_b_rmassx_inner_core = mp->d_rmassx_inner_core;
      mp->d_b_rmassy_inner_core = mp->d_rmassy_inner_core;
    }

    // kernels
    size = NGLL3*(mp->NSPEC_INNER_CORE);

    // density kernel
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rho_kl_inner_core),
                                       size*sizeof(realw)),5204);
    print_CUDA_error_if_any(cudaMemset(mp->d_rho_kl_inner_core,0,size*sizeof(realw)),5207);

    // isotropic kernel
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_alpha_kl_inner_core),
                                       size*sizeof(realw)),5205);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_beta_kl_inner_core),
                                       size*sizeof(realw)),5205);
    print_CUDA_error_if_any(cudaMemset(mp->d_alpha_kl_inner_core,0,size*sizeof(realw)),5208);
    print_CUDA_error_if_any(cudaMemset(mp->d_beta_kl_inner_core,0,size*sizeof(realw)),5208);
  }

  // mesh coloring
  mp->num_colors_outer_inner_core = *num_colors_outer;
  mp->num_colors_inner_inner_core = *num_colors_inner;
  mp->h_num_elem_colors_inner_core = (int*) num_elem_colors;

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_inner_core_device");
#endif
}



/* ----------------------------------------------------------------------------------------------- */

// cleanup

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_cleanup_device,
              PREPARE_CLEANUP_DEVICE)(long* Mesh_pointer_f,
              int* NCHUNKS_VAL) {

TRACE("prepare_cleanup_device");

  // frees allocated memory arrays
  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // frees memory on GPU

  //------------------------------------------
  // sources
  //------------------------------------------
  if( mp->simulation_type == 1  || mp->simulation_type == 3 ){
    cudaFree(mp->d_sourcearrays);
    cudaFree(mp->d_stf_pre_compute);
  }

  cudaFree(mp->d_islice_selected_source);
  cudaFree(mp->d_ispec_selected_source);

  //------------------------------------------
  // receivers
  //------------------------------------------
  if( mp->nrec_local > 0 ) {
    cudaFree(mp->d_number_receiver_global);
    cudaFree(mp->d_station_seismo_field);
    free(mp->h_station_seismo_field);
  }
  cudaFree(mp->d_ispec_selected_rec);

  if( mp->nadj_rec_local > 0 ){
    cudaFree(mp->d_adj_sourcearrays);
    cudaFree(mp->d_pre_computed_irec);
    free(mp->h_adj_sourcearrays_slice);
  }

  //------------------------------------------
  // rotation arrays
  //------------------------------------------
  if( mp->rotation ){
    cudaFree(mp->d_A_array_rotation);
    cudaFree(mp->d_B_array_rotation);
    if( mp->simulation_type == 3 ){
      cudaFree(mp->d_b_A_array_rotation);
      cudaFree(mp->d_b_B_array_rotation);
    }
  }

  //------------------------------------------
  // gravity arrays
  //------------------------------------------
  if( ! mp->gravity ){
    cudaFree(mp->d_d_ln_density_dr_table);
  }else{
    cudaFree(mp->d_minus_rho_g_over_kappa_fluid);
    cudaFree(mp->d_minus_gravity_table);
    cudaFree(mp->d_minus_deriv_gravity_table);
    cudaFree(mp->d_density_table);
  }

  //------------------------------------------
  // attenuation arrays
  //------------------------------------------
  if( mp->attenuation ){
    cudaFree(mp->d_one_minus_sum_beta_crust_mantle);
    cudaFree(mp->d_one_minus_sum_beta_inner_core);
    if( ! mp->partial_phys_dispersion_only ){
      cudaFree(mp->d_factor_common_crust_mantle);
      cudaFree(mp->d_R_xx_crust_mantle);
      cudaFree(mp->d_R_yy_crust_mantle);
      cudaFree(mp->d_R_xy_crust_mantle);
      cudaFree(mp->d_R_xz_crust_mantle);
      cudaFree(mp->d_R_yz_crust_mantle);
      cudaFree(mp->d_factor_common_inner_core);
      cudaFree(mp->d_R_xx_inner_core);
      cudaFree(mp->d_R_yy_inner_core);
      cudaFree(mp->d_R_xy_inner_core);
      cudaFree(mp->d_R_xz_inner_core);
      cudaFree(mp->d_R_yz_inner_core);
    }
    cudaFree(mp->d_alphaval);
    cudaFree(mp->d_betaval);
    cudaFree(mp->d_gammaval);
    if( mp->simulation_type == 3 ){
      cudaFree(mp->d_b_alphaval);
      cudaFree(mp->d_b_betaval);
      cudaFree(mp->d_b_gammaval);
    }
  }

  //------------------------------------------
  // strain
  //------------------------------------------
  if( mp->compute_and_store_strain ){
    cudaFree(mp->d_epsilondev_xx_crust_mantle);
    cudaFree(mp->d_epsilondev_yy_crust_mantle);
    cudaFree(mp->d_epsilondev_xy_crust_mantle);
    cudaFree(mp->d_epsilondev_xz_crust_mantle);
    cudaFree(mp->d_epsilondev_yz_crust_mantle);

    cudaFree(mp->d_epsilondev_xx_inner_core);
    cudaFree(mp->d_epsilondev_yy_inner_core);
    cudaFree(mp->d_epsilondev_xy_inner_core);
    cudaFree(mp->d_epsilondev_xz_inner_core);
    cudaFree(mp->d_epsilondev_yz_inner_core);

    cudaFree(mp->d_eps_trace_over_3_crust_mantle);
    cudaFree(mp->d_eps_trace_over_3_inner_core);
    if( mp->simulation_type == 3 && ! mp->undo_attenuation ){
      cudaFree(mp->d_b_epsilondev_xx_crust_mantle);
      cudaFree(mp->d_b_epsilondev_yy_crust_mantle);
      cudaFree(mp->d_b_epsilondev_xy_crust_mantle);
      cudaFree(mp->d_b_epsilondev_xz_crust_mantle);
      cudaFree(mp->d_b_epsilondev_yz_crust_mantle);

      cudaFree(mp->d_b_epsilondev_xx_inner_core);
      cudaFree(mp->d_b_epsilondev_yy_inner_core);
      cudaFree(mp->d_b_epsilondev_xy_inner_core);
      cudaFree(mp->d_b_epsilondev_xz_inner_core);
      cudaFree(mp->d_b_epsilondev_yz_inner_core);

      cudaFree(mp->d_b_eps_trace_over_3_crust_mantle);
      cudaFree(mp->d_b_eps_trace_over_3_inner_core);
    }
  }

  //------------------------------------------
  // absorbing boundaries arrays
  //------------------------------------------
  if( mp->absorbing_conditions){
    cudaFree(mp->d_rho_vp_crust_mantle);
    cudaFree(mp->d_rho_vs_crust_mantle);
    cudaFree(mp->d_nkmin_xi_crust_mantle);
    cudaFree(mp->d_nkmin_eta_crust_mantle);
    cudaFree(mp->d_njmin_crust_mantle);
    cudaFree(mp->d_njmax_crust_mantle);
    cudaFree(mp->d_nimin_crust_mantle);
    cudaFree(mp->d_nimax_crust_mantle);
    if( mp->nspec2D_xmin_crust_mantle > 0 ){
      cudaFree(mp->d_ibelm_xmin_crust_mantle);
      cudaFree(mp->d_normal_xmin_crust_mantle);
      cudaFree(mp->d_jacobian2D_xmin_crust_mantle);
      if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
        cudaFree(mp->d_absorb_xmin_crust_mantle);
      }
    }
    if( mp->nspec2D_xmax_crust_mantle > 0 ){
      cudaFree(mp->d_ibelm_xmax_crust_mantle);
      cudaFree(mp->d_normal_xmax_crust_mantle);
      cudaFree(mp->d_jacobian2D_xmax_crust_mantle);
      if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
        cudaFree(mp->d_absorb_xmax_crust_mantle);
      }
    }
    if( mp->nspec2D_ymin_crust_mantle > 0 ){
      cudaFree(mp->d_ibelm_ymin_crust_mantle);
      cudaFree(mp->d_normal_ymin_crust_mantle);
      cudaFree(mp->d_jacobian2D_ymin_crust_mantle);
      if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
        cudaFree(mp->d_absorb_ymin_crust_mantle);
      }
    }
    if( mp->nspec2D_ymax_crust_mantle > 0 ){
      cudaFree(mp->d_ibelm_ymax_crust_mantle);
      cudaFree(mp->d_normal_ymax_crust_mantle);
      cudaFree(mp->d_jacobian2D_ymax_crust_mantle);
      if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
        cudaFree(mp->d_absorb_ymax_crust_mantle);
      }
    }

    cudaFree(mp->d_vp_outer_core);
    cudaFree(mp->d_nkmin_xi_outer_core);
    cudaFree(mp->d_nkmin_eta_outer_core);
    cudaFree(mp->d_njmin_outer_core);
    cudaFree(mp->d_njmax_outer_core);
    cudaFree(mp->d_nimin_outer_core);
    cudaFree(mp->d_nimax_outer_core);
    if( mp->nspec2D_xmin_outer_core > 0 ){
      cudaFree(mp->d_ibelm_xmin_outer_core);
      cudaFree(mp->d_jacobian2D_xmin_outer_core);
      if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
        cudaFree(mp->d_absorb_xmin_outer_core);
      }
    }
    if( mp->nspec2D_xmax_outer_core > 0 ){
      cudaFree(mp->d_ibelm_xmax_outer_core);
      cudaFree(mp->d_jacobian2D_xmax_outer_core);
      if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
        cudaFree(mp->d_absorb_xmax_outer_core);
      }
    }
    if( mp->nspec2D_ymin_outer_core > 0 ){
      cudaFree(mp->d_ibelm_ymin_outer_core);
      cudaFree(mp->d_jacobian2D_ymin_outer_core);
      if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
        cudaFree(mp->d_absorb_ymin_outer_core);
      }
    }
    if( mp->nspec2D_ymax_outer_core > 0 ){
      cudaFree(mp->d_ibelm_ymax_outer_core);
      cudaFree(mp->d_jacobian2D_ymax_outer_core);
      if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
        cudaFree(mp->d_absorb_ymax_outer_core);
      }
    }
    if( mp->nspec2D_zmin_outer_core > 0 ){
      if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
        cudaFree(mp->d_absorb_zmin_outer_core);
      }
    }

  }

  //------------------------------------------
  // mpi buffers
  //------------------------------------------
  if( mp->num_interfaces_crust_mantle > 0 ){
    cudaFree(mp->d_nibool_interfaces_crust_mantle);
    cudaFree(mp->d_ibool_interfaces_crust_mantle);
    cudaFree(mp->d_send_accel_buffer_crust_mantle);
    if( mp->simulation_type == 3 ) cudaFree(mp->d_b_send_accel_buffer_crust_mantle);
  }
  if( mp->num_interfaces_inner_core > 0 ){
    cudaFree(mp->d_nibool_interfaces_inner_core);
    cudaFree(mp->d_ibool_interfaces_inner_core);
    cudaFree(mp->d_send_accel_buffer_inner_core);
    if( mp->simulation_type == 3 ) cudaFree(mp->d_b_send_accel_buffer_inner_core);
  }
  if( mp->num_interfaces_outer_core > 0 ){
    cudaFree(mp->d_nibool_interfaces_outer_core);
    cudaFree(mp->d_ibool_interfaces_outer_core);
    cudaFree(mp->d_send_accel_buffer_outer_core);
    if( mp->simulation_type == 3 ) cudaFree(mp->d_b_send_accel_buffer_outer_core);
  }

  //------------------------------------------
  // NOISE arrays
  //------------------------------------------
  if( mp->noise_tomography > 0 ){
    cudaFree(mp->d_ibelm_top_crust_mantle);
    cudaFree(mp->d_noise_surface_movie);
    if( mp->noise_tomography == 1 ) cudaFree(mp->d_noise_sourcearray);
    if( mp->noise_tomography > 1 ){
      cudaFree(mp->d_normal_x_noise);
      cudaFree(mp->d_normal_y_noise);
      cudaFree(mp->d_normal_z_noise);
      cudaFree(mp->d_mask_noise);
      cudaFree(mp->d_jacobian2D_top_crust_mantle);
    }
    if( mp->noise_tomography == 3 ) cudaFree(mp->d_Sigma_kl);
  }

  //------------------------------------------
  // crust_mantle
  //------------------------------------------
  cudaFree(mp->d_xix_crust_mantle);
  cudaFree(mp->d_xiy_crust_mantle);
  cudaFree(mp->d_xiz_crust_mantle);
  cudaFree(mp->d_etax_crust_mantle);
  cudaFree(mp->d_etay_crust_mantle);
  cudaFree(mp->d_etaz_crust_mantle);
  cudaFree(mp->d_gammax_crust_mantle);
  cudaFree(mp->d_gammay_crust_mantle);
  cudaFree(mp->d_gammaz_crust_mantle);

  cudaFree(mp->d_muvstore_crust_mantle);
  cudaFree(mp->d_ibool_crust_mantle);

  if( ! mp->anisotropic_3D_mantle ){
    cudaFree(mp->d_kappavstore_crust_mantle);
    cudaFree(mp->d_kappahstore_crust_mantle);
    cudaFree(mp->d_muhstore_crust_mantle);
    cudaFree(mp->d_eta_anisostore_crust_mantle);
    cudaFree(mp->d_ispec_is_tiso_crust_mantle);
  }else{
    cudaFree(mp->d_c11store_crust_mantle);
    cudaFree(mp->d_c12store_crust_mantle);
    cudaFree(mp->d_c13store_crust_mantle);
    cudaFree(mp->d_c14store_crust_mantle);
    cudaFree(mp->d_c15store_crust_mantle);
    cudaFree(mp->d_c16store_crust_mantle);
    cudaFree(mp->d_c22store_crust_mantle);
    cudaFree(mp->d_c23store_crust_mantle);
    cudaFree(mp->d_c24store_crust_mantle);
    cudaFree(mp->d_c25store_crust_mantle);
    cudaFree(mp->d_c26store_crust_mantle);
    cudaFree(mp->d_c33store_crust_mantle);
    cudaFree(mp->d_c34store_crust_mantle);
    cudaFree(mp->d_c35store_crust_mantle);
    cudaFree(mp->d_c36store_crust_mantle);
    cudaFree(mp->d_c44store_crust_mantle);
    cudaFree(mp->d_c45store_crust_mantle);
    cudaFree(mp->d_c46store_crust_mantle);
    cudaFree(mp->d_c55store_crust_mantle);
    cudaFree(mp->d_c56store_crust_mantle);
    cudaFree(mp->d_c66store_crust_mantle);
  }

  if( mp->simulation_type == 3 && mp->save_boundary_mesh ){
    cudaFree(mp->d_rhostore_crust_mantle);
  }

  cudaFree(mp->d_ystore_crust_mantle);
  cudaFree(mp->d_zstore_crust_mantle);
  if( mp->gravity ){
    cudaFree(mp->d_xstore_crust_mantle);
  }

  cudaFree(mp->d_phase_ispec_inner_crust_mantle);
  cudaFree(mp->d_ibelm_bottom_crust_mantle);

  cudaFree(mp->d_displ_crust_mantle);
  cudaFree(mp->d_veloc_crust_mantle);
  cudaFree(mp->d_accel_crust_mantle);
  if( mp->simulation_type == 3 ){
    cudaFree(mp->d_b_displ_crust_mantle);
    cudaFree(mp->d_b_veloc_crust_mantle);
    cudaFree(mp->d_b_accel_crust_mantle);
    cudaFree(mp->d_rho_kl_crust_mantle);
    if(mp->anisotropic_kl){
      cudaFree(mp->d_cijkl_kl_crust_mantle);
    }else{
      cudaFree(mp->d_alpha_kl_crust_mantle);
      cudaFree(mp->d_beta_kl_crust_mantle);
    }
    if(mp->approximate_hess_kl){ cudaFree(mp->d_hess_kl_crust_mantle);}
  }
  // mass matrix
  if( *NCHUNKS_VAL != 6 && mp->absorbing_conditions){
    cudaFree(mp->d_rmassx_crust_mantle);
    cudaFree(mp->d_rmassy_crust_mantle);
  }
  cudaFree(mp->d_rmassz_crust_mantle);

  //------------------------------------------
  // outer_core
  //------------------------------------------
  cudaFree(mp->d_xix_outer_core);
  cudaFree(mp->d_xiy_outer_core);
  cudaFree(mp->d_xiz_outer_core);
  cudaFree(mp->d_etax_outer_core);
  cudaFree(mp->d_etay_outer_core);
  cudaFree(mp->d_etaz_outer_core);
  cudaFree(mp->d_gammax_outer_core);
  cudaFree(mp->d_gammay_outer_core);
  cudaFree(mp->d_gammaz_outer_core);

  cudaFree(mp->d_kappavstore_outer_core);
  if( mp->simulation_type == 3 ){
    cudaFree(mp->d_rhostore_outer_core);
  }

  cudaFree(mp->d_xstore_outer_core);
  cudaFree(mp->d_ystore_outer_core);
  cudaFree(mp->d_zstore_outer_core);

  cudaFree(mp->d_ibool_outer_core);
  cudaFree(mp->d_phase_ispec_inner_outer_core);

  cudaFree(mp->d_ibelm_top_outer_core);
  cudaFree(mp->d_jacobian2D_top_outer_core);
  cudaFree(mp->d_normal_top_outer_core);

  cudaFree(mp->d_ibelm_bottom_outer_core);
  cudaFree(mp->d_normal_bottom_outer_core);
  cudaFree(mp->d_jacobian2D_bottom_outer_core);

  cudaFree(mp->d_displ_outer_core);
  cudaFree(mp->d_veloc_outer_core);
  cudaFree(mp->d_accel_outer_core);
  if( mp->simulation_type == 3 ){
    cudaFree(mp->d_b_displ_outer_core);
    cudaFree(mp->d_b_veloc_outer_core);
    cudaFree(mp->d_b_accel_outer_core);
    cudaFree(mp->d_rho_kl_outer_core);
    cudaFree(mp->d_alpha_kl_outer_core);
  }
  // mass matrix
  cudaFree(mp->d_rmass_outer_core);

  //------------------------------------------
  // inner_core
  //------------------------------------------
  cudaFree(mp->d_xix_inner_core);
  cudaFree(mp->d_xiy_inner_core);
  cudaFree(mp->d_xiz_inner_core);
  cudaFree(mp->d_etax_inner_core);
  cudaFree(mp->d_etay_inner_core);
  cudaFree(mp->d_etaz_inner_core);
  cudaFree(mp->d_gammax_inner_core);
  cudaFree(mp->d_gammay_inner_core);
  cudaFree(mp->d_gammaz_inner_core);

  cudaFree(mp->d_muvstore_inner_core);
  cudaFree(mp->d_ibool_inner_core);

  // gravity
  if( mp->gravity ){
    cudaFree(mp->d_xstore_inner_core);
    cudaFree(mp->d_ystore_inner_core);
    cudaFree(mp->d_zstore_inner_core);
  }

  cudaFree(mp->d_ibelm_top_inner_core);

  if( ! mp->anisotropic_inner_core ){
    cudaFree(mp->d_kappavstore_inner_core);
  }else{
    cudaFree(mp->d_c11store_inner_core);
    cudaFree(mp->d_c12store_inner_core);
    cudaFree(mp->d_c13store_inner_core);
    cudaFree(mp->d_c33store_inner_core);
    cudaFree(mp->d_c44store_inner_core);
  }

  if( mp->simulation_type == 3 && mp->save_boundary_mesh ){
    cudaFree(mp->d_rhostore_inner_core);
  }
  cudaFree(mp->d_idoubling_inner_core);
  if( mp->gravity ){
    cudaFree(mp->d_xstore_inner_core);
    cudaFree(mp->d_ystore_inner_core);
    cudaFree(mp->d_zstore_inner_core);
  }
  cudaFree(mp->d_phase_ispec_inner_inner_core);

  cudaFree(mp->d_displ_inner_core);
  cudaFree(mp->d_veloc_inner_core);
  cudaFree(mp->d_accel_inner_core);
  if( mp->simulation_type == 3 ) {
    cudaFree(mp->d_b_displ_inner_core);
    cudaFree(mp->d_b_veloc_inner_core);
    cudaFree(mp->d_b_accel_inner_core);

    cudaFree(mp->d_rho_kl_inner_core);
    cudaFree(mp->d_alpha_kl_inner_core);
    cudaFree(mp->d_beta_kl_inner_core);
  }
  // mass matrix
  cudaFree(mp->d_rmassx_inner_core);
  cudaFree(mp->d_rmassy_inner_core);
  cudaFree(mp->d_rmassz_inner_core);

  // oceans
  if( mp->oceans ){
    cudaFree(mp->d_rmass_ocean_load);
    cudaFree(mp->d_ibool_ocean_load);
    cudaFree(mp->d_normal_ocean_load);
  }

  // releases previous contexts
#if CUDA_VERSION < 4000
  cudaThreadSynchronize();
  cudaThreadExit();
#else
  cudaDeviceSynchronize();
  cudaDeviceReset();
#endif

  // mesh pointer - not needed anymore
  free(mp);
}

