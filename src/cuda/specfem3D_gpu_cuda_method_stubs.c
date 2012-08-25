/* 
!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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
#include <stdlib.h>

#include "config.h"

typedef float realw;

 

//
// src/cuda/assemble_MPI_scalar_cuda.cu
//

void FC_FUNC_(transfer_boun_pot_from_device,
              TRANSFER_BOUN_POT_FROM_DEVICE)(long* Mesh_pointer_f,
                                             realw* send_potential_dot_dot_buffer,
                                             int* FORWARD_OR_ADJOINT){} 

void FC_FUNC_(transfer_asmbl_pot_to_device,
              TRANSFER_ASMBL_POT_TO_DEVICE)(long* Mesh_pointer,
                                            realw* buffer_recv_scalar,
                                            int* FORWARD_OR_ADJOINT) {} 


//
// src/cuda/assemble_MPI_vector_cuda.cu
//

void FC_FUNC_(transfer_boun_accel_from_device,
              TRANSFER_BOUN_ACCEL_FROM_DEVICE)(long* Mesh_pointer_f,
                                                  realw* send_accel_buffer,
                                                  int* IREGION,
                                                  int* FORWARD_OR_ADJOINT){} 

void FC_FUNC_(transfer_asmbl_accel_to_device,
              TRANSFER_ASMBL_ACCEL_TO_DEVICE)(long* Mesh_pointer,
                                              realw* buffer_recv_vector,
                                              int* IREGION,
                                              int* FORWARD_OR_ADJOINT) {} 


//
// src/cuda/check_fields_cuda.cu
//

void FC_FUNC_(pause_for_debug,
              PAUSE_FOR_DEBUG)() {} 

void FC_FUNC_(output_free_device_memory,
              OUTPUT_FREE_DEVICE_MEMORY)(int* myrank) {} 

void FC_FUNC_(get_free_device_memory,
              get_FREE_DEVICE_MEMORY)(realw* free, realw* used, realw* total ) {} 

void FC_FUNC_(check_norm_acoustic_from_device,
              CHECK_NORM_ACOUSTIC_FROM_DEVICE)(realw* norm,
                                                  long* Mesh_pointer_f,
                                                  int* SIMULATION_TYPE) {} 

void FC_FUNC_(check_norm_elastic_from_device,
              CHECK_NORM_ELASTIC_FROM_DEVICE)(realw* norm,
                                              long* Mesh_pointer_f,
                                              int* SIMULATION_TYPE) {} 

void FC_FUNC_(check_norm_strain_from_device,
              CHECK_NORM_STRAIN_FROM_DEVICE)(realw* strain_norm,
                                             realw* strain_norm2,
                                             long* Mesh_pointer_f) {} 

 void FC_FUNC_(check_max_norm_displ_gpu,
 CHECK_MAX_NORM_DISPL_GPU)(int* size, realw* displ,long* Mesh_pointer_f,int* announceID) {} 

 void FC_FUNC_(check_max_norm_vector,
 CHECK_MAX_NORM_VECTOR)(int* size, realw* vector1, int* announceID) {} 

 void FC_FUNC_(check_max_norm_displ,
 CHECK_MAX_NORM_DISPL)(int* size, realw* displ, int* announceID) {} 

 void FC_FUNC_(check_max_norm_b_displ_gpu,
 CHECK_MAX_NORM_B_DISPL_GPU)(int* size, realw* b_displ,long* Mesh_pointer_f,int* announceID) {} 

 void FC_FUNC_(check_max_norm_b_accel_gpu,
 CHECK_MAX_NORM_B_ACCEL_GPU)(int* size, realw* b_accel,long* Mesh_pointer_f,int* announceID) {} 

 void FC_FUNC_(check_max_norm_b_veloc_gpu,
 CHECK_MAX_NORM_B_VELOC_GPU)(int* size, realw* b_veloc,long* Mesh_pointer_f,int* announceID) {} 

 void FC_FUNC_(check_max_norm_b_displ,
 CHECK_MAX_NORM_B_DISPL)(int* size, realw* b_displ,int* announceID) {} 

 void FC_FUNC_(check_max_norm_b_accel,
 CHECK_MAX_NORM_B_ACCEL)(int* size, realw* b_accel,int* announceID) {} 

 void FC_FUNC_(check_error_vectors,
 CHECK_ERROR_VECTORS)(int* sizef, realw* vector1,realw* vector2) {} 

 void FC_FUNC_(get_max_accel,
 GET_MAX_ACCEL)(int* itf,int* sizef,long* Mesh_pointer) {} 


//
// src/cuda/compute_add_sources_elastic_cuda.cu
//

void FC_FUNC_(compute_add_sources_el_cuda,
              COMPUTE_ADD_SOURCES_EL_CUDA)(long* Mesh_pointer_f,
                                           int* NSOURCESf,
                                           double* h_stf_pre_compute) {} 

void FC_FUNC_(compute_add_sources_el_s3_cuda,
              COMPUTE_ADD_SOURCES_EL_S3_CUDA)(long* Mesh_pointer_f,
                                              int* NSOURCESf,
                                              double* h_stf_pre_compute) {} 

void FC_FUNC_(add_sources_el_sim_type_2_or_3,
              ADD_SOURCES_EL_SIM_TYPE_2_OR_3)(long* Mesh_pointer,
                                              int* nrec,
                                              realw* h_adj_sourcearrays,
                                              int* h_islice_selected_rec,
                                              int* h_ispec_selected_rec,
                                              int* time_index) {} 


//
// src/cuda/compute_coupling_cuda.cu
//

void FC_FUNC_(compute_coupling_fluid_cmb_cuda,
              COMPUTE_COUPLING_FLUID_CMB_CUDA)(long* Mesh_pointer_f) {} 

void FC_FUNC_(compute_coupling_fluid_icb_cuda,
              COMPUTE_COUPLING_FLUID_ICB_CUDA)(long* Mesh_pointer_f) {} 

void FC_FUNC_(compute_coupling_cmb_fluid_cuda,
              COMPUTE_COUPLING_CMB_FLUID_CUDA)(long* Mesh_pointer_f) {} 

void FC_FUNC_(compute_coupling_icb_fluid_cuda,
              COMPUTE_COUPLING_ICB_FLUID_CUDA)(long* Mesh_pointer_f) {} 

void FC_FUNC_(compute_coupling_ocean_cuda,
              COMPUTE_COUPLING_OCEAN_CUDA)(long* Mesh_pointer_f,
                                           int* NCHUNKS_VAL) {} 


//
// src/cuda/compute_forces_crust_mantle_cuda.cu
//

void FC_FUNC_(compute_forces_crust_mantle_cuda,
              COMPUTE_FORCES_CRUST_MANTLE_CUDA)(long* Mesh_pointer_f,
                                                int* iphase) {} 


//
// src/cuda/compute_forces_inner_core_cuda.cu
//

void FC_FUNC_(compute_forces_inner_core_cuda,
              COMPUTE_FORCES_INNER_CORE_CUDA)(long* Mesh_pointer_f,
                                              int* iphase) {} 


//
// src/cuda/compute_forces_outer_core_cuda.cu
//

void FC_FUNC_(compute_forces_outer_core_cuda,
              COMPUTE_FORCES_OUTER_CORE_CUDA)(long* Mesh_pointer_f,
                                              int* iphase,
                                              realw* time_f,
                                              realw* b_time_f) {} 


//
// src/cuda/compute_kernels_cuda.cu
//

void FC_FUNC_(compute_kernels_cm_cuda,
              COMPUTE_KERNELS_CM_CUDA)(long* Mesh_pointer,realw* deltat_f) {} 

void FC_FUNC_(compute_kernels_ic_cuda,
              COMPUTE_KERNELS_IC_CUDA)(long* Mesh_pointer,realw* deltat_f) {} 

void FC_FUNC_(compute_kernels_oc_cuda,
              COMPUTE_KERNELS_OC_CUDA)(long* Mesh_pointer,realw* deltat_f) {} 

void FC_FUNC_(compute_kernels_strgth_noise_cu,
              COMPUTE_KERNELS_STRGTH_NOISE_CU)(long* Mesh_pointer,
                                               realw* h_noise_surface_movie,
                                               realw* deltat_f) {} 

void FC_FUNC_(compute_kernels_hess_cuda,
              COMPUTE_KERNELS_HESS_CUDA)(long* Mesh_pointer,
                                         realw* deltat_f) {} 


//
// src/cuda/compute_stacey_acoustic_cuda.cu
//

void FC_FUNC_(compute_stacey_acoustic_cuda,
              COMPUTE_STACEY_ACOUSTIC_CUDA)(long* Mesh_pointer_f,
                                            realw* absorb_potential,
                                            int* itype) {} 


//
// src/cuda/compute_stacey_elastic_cuda.cu
//

void FC_FUNC_(compute_stacey_elastic_cuda,
              COMPUTE_STACEY_ELASTIC_CUDA)(long* Mesh_pointer_f,
                                                realw* absorb_field,
                                                int* itype) {} 


//
// src/cuda/initialize_cuda.cu
//

void FC_FUNC_(initialize_cuda_device,
              INITIALIZE_CUDA_DEVICE)(int* myrank_f,int* ncuda_devices) { 
 fprintf(stderr,"ERROR: GPU_MODE enabled without GPU/CUDA Support. To enable GPU support, reconfigure with --with-cuda flag.\n");
 exit(1);
} 


//
// src/cuda/it_update_displacement_cuda.cu
//

void FC_FUNC_(it_update_displacement_ic_cuda,
              IT_UPDATE_DISPLACMENT_IC_CUDA)(long* Mesh_pointer_f,
                                             realw* deltat_F,
                                             realw* deltatsqover2_F,
                                             realw* deltatover2_F,
                                             realw* b_deltat_F,
                                             realw* b_deltatsqover2_F,
                                             realw* b_deltatover2_F) {} 

void FC_FUNC_(it_update_displacement_cm_cuda,
              IT_UPDATE_DISPLACMENT_CM_CUDA)(long* Mesh_pointer_f,
                                             realw* deltat_F,
                                             realw* deltatsqover2_F,
                                             realw* deltatover2_F,
                                             realw* b_deltat_F,
                                             realw* b_deltatsqover2_F,
                                             realw* b_deltatover2_F) {} 

void FC_FUNC_(it_update_displacement_oc_cuda,
              IT_UPDATE_DISPLACEMENT_OC_cuda)(long* Mesh_pointer_f,
                                               realw* deltat_F,
                                               realw* deltatsqover2_F,
                                               realw* deltatover2_F,
                                               realw* b_deltat_F,
                                               realw* b_deltatsqover2_F,
                                               realw* b_deltatover2_F) {} 

void FC_FUNC_(kernel_3_a_cuda,
              KERNEL_3_A_CUDA)(long* Mesh_pointer,
                               realw* deltatover2_F,
                               int* SIMULATION_TYPE_f,
                               realw* b_deltatover2_F,
                               int* NCHUNKS_VAL) {} 

void FC_FUNC_(kernel_3_b_cuda,
              KERNEL_3_B_CUDA)(long* Mesh_pointer,
                               realw* deltatover2_F,
                               int* SIMULATION_TYPE_f,
                               realw* b_deltatover2_F) {} 

void FC_FUNC_(kernel_3_outer_core_cuda,
              KERNEL_3_OUTER_CORE_CUDA)(long* Mesh_pointer,
                                        realw* deltatover2_F,
                                        int* SIMULATION_TYPE_f,
                                        realw* b_deltatover2_F) {} 


//
// src/cuda/noise_tomography_cuda.cu
//

void FC_FUNC_(fortranflush,FORTRANFLUSH)(int* rank){} 

void FC_FUNC_(fortranprint,FORTRANPRINT)(int* id) {} 

void FC_FUNC_(fortranprintf,FORTRANPRINTF)(realw* val) {} 

void FC_FUNC_(fortranprintd,FORTRANPRINTD)(double* val) {} 

void FC_FUNC_(make_displ_rand,MAKE_DISPL_RAND)(long* Mesh_pointer_f,realw* h_displ) {} 

void FC_FUNC_(noise_transfer_surface_to_host,
              NOISE_TRANSFER_SURFACE_TO_HOST)(long* Mesh_pointer_f,
                                              realw* h_noise_surface_movie) {} 

void FC_FUNC_(noise_add_source_master_rec_cu,
              NOISE_ADD_SOURCE_MASTER_REC_CU)(long* Mesh_pointer_f,
                                              int* it_f,
                                              int* irec_master_noise_f,
                                              int* islice_selected_rec) {} 

void FC_FUNC_(noise_add_surface_movie_cuda,
              NOISE_ADD_SURFACE_MOVIE_CUDA)(long* Mesh_pointer_f,
                                            realw* h_noise_surface_movie) {} 


//
// src/cuda/prepare_mesh_constants_cuda.cu
//

void FC_FUNC_(prepare_constants_device,
              PREPARE_CONSTANTS_DEVICE)(long* Mesh_pointer,
                                        int* myrank_f,
                                        int* h_NGLLX,
                                        realw* h_hprime_xx,
                                        realw* h_hprimewgll_xx,
                                        realw* h_wgllwgll_xy,realw* h_wgllwgll_xz,realw* h_wgllwgll_yz,
                                        int* NSOURCES,int* nsources_local,
                                        realw* h_sourcearrays,
                                        int* h_islice_selected_source,int* h_ispec_selected_source,
                                        int* h_number_receiver_global,
                                        int* h_islice_selected_rec,int* h_ispec_selected_rec,
                                        int* nrec,int* nrec_local, int* nadj_rec_local,
                                        int* NSPEC_CRUST_MANTLE, int* NGLOB_CRUST_MANTLE,
                                        int* NSPEC_CRUST_MANTLE_STRAIN_ONLY,
                                        int* NSPEC_OUTER_CORE, int* NGLOB_OUTER_CORE,
                                        int* NSPEC_INNER_CORE, int* NGLOB_INNER_CORE,
                                        int* NSPEC_INNER_CORE_STRAIN_ONLY,
                                        int* SIMULATION_TYPE,
                                        int* NOISE_TOMOGRAPHY,
                                        int* SAVE_FORWARD_f,
                                        int* ABSORBING_CONDITIONS_f,
                                        int* OCEANS_f,
                                        int* GRAVITY_f,
                                        int* ROTATION_f,
                                        int* ATTENUATION_f,
                                        int* ATTENUATION_NEW_f,
                                        int* USE_ATTENUATION_MIMIC_f,
                                        int* ATTENUATION_3D_VAL_f,
                                        int* COMPUTE_AND_STORE_STRAIN_f,
                                        int* ANISOTROPIC_3D_MANTLE_f,
                                        int* ANISOTROPIC_INNER_CORE_f,
                                        int* SAVE_BOUNDARY_MESH_f,
                                        int* USE_MESH_COLORING_GPU_f,
                                        int* ANISOTROPIC_KL_f,
                                        int* APPROXIMATE_HESS_KL_f,
                                        realw* deltat_f,
                                        realw* b_deltat_f) {} 

void FC_FUNC_(prepare_fields_rotation_device,
              PREPARE_FIELDS_ROTATION_DEVICE)(long* Mesh_pointer_f,
                                              realw* two_omega_earth_f,
                                              realw* A_array_rotation,
                                              realw* B_array_rotation,
                                              realw* b_two_omega_earth_f,
                                              realw* b_A_array_rotation,
                                              realw* b_B_array_rotation,
                                              int* NSPEC_OUTER_CORE_ROTATION) {} 

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
                                             double* RHO_TOP_OC) {} 

void FC_FUNC_(prepare_fields_attenuat_device,
              PREPARE_FIELDS_ATTENUAT_DEVICE)(long* Mesh_pointer_f,
                                                 realw* R_xx_crust_mantle,
                                                 realw* R_yy_crust_mantle,
                                                 realw* R_xy_crust_mantle,
                                                 realw* R_xz_crust_mantle,
                                                 realw* R_yz_crust_mantle,
                                                 realw* factor_common_crust_mantle,
                                                 realw* one_minus_sum_beta_crust_mantle,
                                                 realw* R_xx_inner_core,
                                                 realw* R_yy_inner_core,
                                                 realw* R_xy_inner_core,
                                                 realw* R_xz_inner_core,
                                                 realw* R_yz_inner_core,
                                                 realw* factor_common_inner_core,
                                                 realw* one_minus_sum_beta_inner_core,
                                                 realw* alphaval,realw* betaval,realw* gammaval,
                                                 realw* b_alphaval,realw* b_betaval,realw* b_gammaval) {} 

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
                                            realw* b_eps_trace_over_3_inner_core) {} 

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
                                            realw* vp_outer_core) {} 

void FC_FUNC_(prepare_mpi_buffers_device,
              PREPARE_MPI_BUFFERS_DEVICE)(long* Mesh_pointer_f,
                                          int* num_interfaces_crust_mantle,
                                          int* max_nibool_interfaces_crust_mantle,
                                          int* nibool_interfaces_crust_mantle,
                                          int* ibool_interfaces_crust_mantle,
                                          int* num_interfaces_inner_core,
                                          int* max_nibool_interfaces_inner_core,
                                          int* nibool_interfaces_inner_core,
                                          int* ibool_interfaces_inner_core,
                                          int* num_interfaces_outer_core,
                                          int* max_nibool_interfaces_outer_core,
                                          int* nibool_interfaces_outer_core,
                                          int* ibool_interfaces_outer_core){} 

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
                                           realw* jacobian2D_top_crust_mantle) {} 

void FC_FUNC_(prepare_oceans_device,
              PREPARE_OCEANS_DEVICE)(long* Mesh_pointer_f,
                                     int* npoin_oceans,
                                     int* h_iglob_ocean_load,
                                     realw* h_rmass_ocean_load_selected,
                                     realw* h_normal_ocean_load) {} 

void FC_FUNC_(prepare_crust_mantle_device,
             PREPARE_CRUST_MANTLE_DEVICE)(long* Mesh_pointer_f,
                                          realw* h_xix, realw* h_xiy, realw* h_xiz,
                                          realw* h_etax, realw* h_etay, realw* h_etaz,
                                          realw* h_gammax, realw* h_gammay, realw* h_gammaz,
                                          realw* h_rho,
                                          realw* h_kappav, realw* h_muv,
                                          realw* h_kappah, realw* h_muh,
                                          realw* h_eta_aniso,
                                          realw* h_rmassx,
                                          realw* h_rmassy,
                                          realw* h_rmassz,
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
                                          int* num_elem_colors) {} 

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
                                         int* num_elem_colors) {} 

void FC_FUNC_(prepare_inner_core_device,
              PREPARE_INNER_CORE_DEVICE)(long* Mesh_pointer_f,
                                         realw* h_xix, realw* h_xiy, realw* h_xiz,
                                         realw* h_etax, realw* h_etay, realw* h_etaz,
                                         realw* h_gammax, realw* h_gammay, realw* h_gammaz,
                                         realw* h_rho, realw* h_kappav, realw* h_muv,
                                         realw* h_rmass,
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
                                         int* num_elem_colors) {} 

void FC_FUNC_(prepare_cleanup_device,
              PREPARE_CLEANUP_DEVICE)(long* Mesh_pointer_f,
              int* NCHUNKS_VAL) {} 


//
// src/cuda/transfer_fields_cuda.cu
//

void FC_FUNC_(transfer_fields_cm_to_device,
              TRANSFER_FIELDS_CM_TO_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_fields_ic_to_device,
              TRANSFER_FIELDS_IC_TO_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_fields_oc_to_device,
              TRANSFER_FIELDS_OC_TO_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_b_fields_cm_to_device,
              TRANSFER_FIELDS_B_CM_TO_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                              long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_b_fields_ic_to_device,
              TRANSFER_FIELDS_B_IC_TO_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                              long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_b_fields_oc_to_device,
              TRANSFER_FIELDS_B_OC_TO_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                              long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_fields_cm_from_device,
              TRANSFER_FIELDS_CM_FROM_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_fields_ic_from_device,
              TRANSFER_FIELDS_IC_FROM_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_fields_oc_from_device,
              TRANSFER_FIELDS_OC_FROM_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_b_fields_cm_from_device,
              TRANSFER_B_FIELDS_CM_FROM_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                                long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_b_fields_ic_from_device,
              TRANSFER_B_FIELDS_IC_FROM_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                                long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_b_fields_oc_from_device,
              TRANSFER_B_FIELDS_OC_FROM_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                                long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_accel_cm_to_device,
              TRANSFER_ACCEL_CM_TO_DEVICE)(int* size, realw* accel,long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_displ_cm_from_device,
              TRANSFER_DISPL_CM_FROM_DEVICE)(int* size, realw* displ, long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_b_displ_cm_from_device,
              TRANSFER_B_DISPL_CM_FROM_DEVICE)(int* size, realw* displ, long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_displ_ic_from_device,
              TRANSFER_DISPL_IC_FROM_DEVICE)(int* size, realw* displ, long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_b_displ_ic_from_device,
              TRANSFER_B_DISPL_IC_FROM_DEVICE)(int* size, realw* displ, long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_displ_oc_from_device,
              TRANSFER_DISPL_OC_FROM_DEVICE)(int* size, realw* displ, long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_b_displ_oc_from_device,
              TRANSFER_B_DISPL_OC_FROM_DEVICE)(int* size, realw* displ, long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_veloc_cm_from_device,
              TRANSFER_DISPL_CM_FROM_DEVICE)(int* size, realw* veloc, long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_accel_cm_from_device,
              TRANSFER_ACCEL_CM_FROM_DEVICE)(int* size, realw* accel,long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_b_accel_cm_from_device,
              TRANSFER_B_ACCEL_CM_FROM_DEVICE)(int* size, realw* b_accel,long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_accel_ic_from_device,
              TRANSFER_ACCEL_IC_FROM_DEVICE)(int* size, realw* accel,long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_accel_oc_from_device,
              TRANSFER_ACCEL_OC_FROM_DEVICE)(int* size, realw* accel,long* Mesh_pointer_f) {} 

void FC_FUNC_(transfer_strain_cm_from_device,
              TRANSFER_STRAIN_CM_FROM_DEVICE)(long* Mesh_pointer,
                                                  realw* eps_trace_over_3,
                                                  realw* epsilondev_xx,
                                                  realw* epsilondev_yy,
                                                  realw* epsilondev_xy,
                                                  realw* epsilondev_xz,
                                                  realw* epsilondev_yz) {} 

void FC_FUNC_(transfer_b_strain_cm_to_device,
              TRANSFER_B_STRAIN_CM_TO_DEVICE)(long* Mesh_pointer,
                                              realw* epsilondev_xx,
                                              realw* epsilondev_yy,
                                              realw* epsilondev_xy,
                                              realw* epsilondev_xz,
                                              realw* epsilondev_yz) {} 

void FC_FUNC_(transfer_strain_ic_from_device,
              TRANSFER_STRAIN_IC_FROM_DEVICE)(long* Mesh_pointer,
                                              realw* eps_trace_over_3,
                                              realw* epsilondev_xx,
                                              realw* epsilondev_yy,
                                              realw* epsilondev_xy,
                                              realw* epsilondev_xz,
                                              realw* epsilondev_yz) {} 

void FC_FUNC_(transfer_b_strain_ic_to_device,
              TRANSFER_B_STRAIN_IC_TO_DEVICE)(long* Mesh_pointer,
                                              realw* epsilondev_xx,
                                              realw* epsilondev_yy,
                                              realw* epsilondev_xy,
                                              realw* epsilondev_xz,
                                              realw* epsilondev_yz) {} 

void FC_FUNC_(transfer_rotation_from_device,
              TRANSFER_ROTATION_FROM_DEVICE)(long* Mesh_pointer,
                                             realw* A_array_rotation,
                                             realw* B_array_rotation) {} 

void FC_FUNC_(transfer_b_rotation_to_device,
              TRANSFER_B_ROTATION_TO_DEVICE)(long* Mesh_pointer,
                                              realw* A_array_rotation,
                                              realw* B_array_rotation) {} 

void FC_FUNC_(transfer_kernels_cm_to_host,
              TRANSFER_KERNELS_CM_TO_HOST)(long* Mesh_pointer,
                                           realw* h_rho_kl,
                                           realw* h_alpha_kl,
                                           realw* h_beta_kl,
                                           realw* h_cijkl_kl,
                                           int* NSPEC) {} 

void FC_FUNC_(transfer_kernels_ic_to_host,
              TRANSFER_KERNELS_IC_TO_HOST)(long* Mesh_pointer,
                                                    realw* h_rho_kl,
                                                    realw* h_alpha_kl,
                                                    realw* h_beta_kl,
                                                    int* NSPEC) {} 

void FC_FUNC_(transfer_kernels_oc_to_host,
              TRANSFER_KERNELS_OC_TO_HOST)(long* Mesh_pointer,
                                           realw* h_rho_kl,
                                           realw* h_alpha_kl,
                                           int* NSPEC) {} 

void FC_FUNC_(transfer_kernels_noise_to_host,
              TRANSFER_KERNELS_NOISE_TO_HOST)(long* Mesh_pointer,
                                              realw* h_Sigma_kl,
                                              int* NSPEC) {} 

void FC_FUNC_(transfer_kernels_hess_cm_tohost,
              TRANSFER_KERNELS_HESS_CM_TOHOST)(long* Mesh_pointer,
                                              realw* h_hess_kl,
                                              int* NSPEC) {} 


//
// src/cuda/write_seismograms_cuda.cu
//

void FC_FUNC_(write_seismograms_transfer_cuda,
              WRITE_SEISMOGRAMS_TRANSFER_CUDA)(long* Mesh_pointer_f,
                                               realw* displ,
                                               realw* b_displ,
                                               realw* eps_trace_over_3,
                                               realw* epsilondev_xx,
                                               realw* epsilondev_yy,
                                               realw* epsilondev_xy,
                                               realw* epsilondev_xz,
                                               realw* epsilondev_yz,
                                               int* number_receiver_global,
                                               int* ispec_selected_rec,
                                               int* ispec_selected_source,
                                               int* ibool) {} 

