/*
!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

// ------------------------------------------------------------------   //
// crust_mantle
// ------------------------------------------------------------------   //
// interpolators
GPU_REALW_BUFFER (d_xix_crust_mantle);
GPU_REALW_BUFFER (d_xiy_crust_mantle);
GPU_REALW_BUFFER (d_xiz_crust_mantle);
GPU_REALW_BUFFER (d_etax_crust_mantle);
GPU_REALW_BUFFER (d_etay_crust_mantle);
GPU_REALW_BUFFER (d_etaz_crust_mantle);
GPU_REALW_BUFFER (d_gammax_crust_mantle);
GPU_REALW_BUFFER (d_gammay_crust_mantle);
GPU_REALW_BUFFER (d_gammaz_crust_mantle);

// model parameters
GPU_REALW_BUFFER (d_rhostore_crust_mantle);
GPU_REALW_BUFFER (d_kappavstore_crust_mantle);
GPU_REALW_BUFFER (d_muvstore_crust_mantle);
GPU_REALW_BUFFER (d_kappahstore_crust_mantle);
GPU_REALW_BUFFER (d_muhstore_crust_mantle);
GPU_REALW_BUFFER (d_eta_anisostore_crust_mantle);

// mass matrices
GPU_REALW_BUFFER (d_rmassx_crust_mantle);
GPU_REALW_BUFFER (d_rmassy_crust_mantle);
GPU_REALW_BUFFER (d_rmassz_crust_mantle);
GPU_REALW_BUFFER (d_b_rmassx_crust_mantle);
GPU_REALW_BUFFER (d_b_rmassy_crust_mantle);
GPU_REALW_BUFFER (d_b_rmassz_crust_mantle);

// global indexing
GPU_INT_BUFFER (d_ibool_crust_mantle);
GPU_INT_BUFFER (d_ispec_is_tiso_crust_mantle);

// mesh locations
GPU_REALW_BUFFER (d_rstore_crust_mantle);

// anisotropic 3D mantle
GPU_REALW_BUFFER (d_c11store_crust_mantle);
GPU_REALW_BUFFER (d_c12store_crust_mantle);
GPU_REALW_BUFFER (d_c13store_crust_mantle);
GPU_REALW_BUFFER (d_c14store_crust_mantle);
GPU_REALW_BUFFER (d_c15store_crust_mantle);
GPU_REALW_BUFFER (d_c16store_crust_mantle);
GPU_REALW_BUFFER (d_c22store_crust_mantle);
GPU_REALW_BUFFER (d_c23store_crust_mantle);
GPU_REALW_BUFFER (d_c24store_crust_mantle);
GPU_REALW_BUFFER (d_c25store_crust_mantle);
GPU_REALW_BUFFER (d_c26store_crust_mantle);
GPU_REALW_BUFFER (d_c33store_crust_mantle);
GPU_REALW_BUFFER (d_c34store_crust_mantle);
GPU_REALW_BUFFER (d_c35store_crust_mantle);
GPU_REALW_BUFFER (d_c36store_crust_mantle);
GPU_REALW_BUFFER (d_c44store_crust_mantle);
GPU_REALW_BUFFER (d_c45store_crust_mantle);
GPU_REALW_BUFFER (d_c46store_crust_mantle);
GPU_REALW_BUFFER (d_c55store_crust_mantle);
GPU_REALW_BUFFER (d_c56store_crust_mantle);
GPU_REALW_BUFFER (d_c66store_crust_mantle);

// wavefields
GPU_REALW_BUFFER (d_displ_crust_mantle);
GPU_REALW_BUFFER (d_veloc_crust_mantle);
GPU_REALW_BUFFER (d_accel_crust_mantle);
GPU_REALW_BUFFER (d_b_displ_crust_mantle);
GPU_REALW_BUFFER (d_b_veloc_crust_mantle);
GPU_REALW_BUFFER (d_b_accel_crust_mantle);

// attenuation
GPU_REALW_BUFFER (d_R_xx_crust_mantle);
GPU_REALW_BUFFER (d_R_yy_crust_mantle);
GPU_REALW_BUFFER (d_R_xy_crust_mantle);
GPU_REALW_BUFFER (d_R_xz_crust_mantle);
GPU_REALW_BUFFER (d_R_yz_crust_mantle);
GPU_REALW_BUFFER (d_b_R_xx_crust_mantle);
GPU_REALW_BUFFER (d_b_R_yy_crust_mantle);
GPU_REALW_BUFFER (d_b_R_xy_crust_mantle);
GPU_REALW_BUFFER (d_b_R_xz_crust_mantle);
GPU_REALW_BUFFER (d_b_R_yz_crust_mantle);
GPU_REALW_BUFFER (d_factor_common_crust_mantle);
GPU_REALW_BUFFER (d_one_minus_sum_beta_crust_mantle);
GPU_REALW_BUFFER (d_epsilondev_xx_crust_mantle);
GPU_REALW_BUFFER (d_epsilondev_yy_crust_mantle);
GPU_REALW_BUFFER (d_epsilondev_xy_crust_mantle);
GPU_REALW_BUFFER (d_epsilondev_xz_crust_mantle);
GPU_REALW_BUFFER (d_epsilondev_yz_crust_mantle);
GPU_REALW_BUFFER (d_b_epsilondev_xx_crust_mantle);
GPU_REALW_BUFFER (d_b_epsilondev_yy_crust_mantle);
GPU_REALW_BUFFER (d_b_epsilondev_xy_crust_mantle);
GPU_REALW_BUFFER (d_b_epsilondev_xz_crust_mantle);
GPU_REALW_BUFFER (d_b_epsilondev_yz_crust_mantle);
GPU_REALW_BUFFER (d_eps_trace_over_3_crust_mantle);
GPU_REALW_BUFFER (d_b_eps_trace_over_3_crust_mantle);

// kernels
GPU_REALW_BUFFER (d_rho_kl_crust_mantle);
GPU_REALW_BUFFER (d_alpha_kl_crust_mantle);
GPU_REALW_BUFFER (d_beta_kl_crust_mantle);
GPU_REALW_BUFFER (d_cijkl_kl_crust_mantle);
GPU_REALW_BUFFER (d_hess_kl_crust_mantle);
GPU_REALW_BUFFER (d_hess_rho_kl_crust_mantle);
GPU_REALW_BUFFER (d_hess_kappa_kl_crust_mantle);
GPU_REALW_BUFFER (d_hess_mu_kl_crust_mantle);

// inner / outer elements
GPU_INT_BUFFER (d_phase_ispec_inner_crust_mantle);
GPU_INT_BUFFER (d_ibelm_bottom_crust_mantle);

// ------------------------------------------------------------------   //
// outer_core
// ------------------------------------------------------------------   //
// interpolators
GPU_REALW_BUFFER (d_xix_outer_core);
GPU_REALW_BUFFER (d_xiy_outer_core);
GPU_REALW_BUFFER (d_xiz_outer_core);
GPU_REALW_BUFFER (d_etax_outer_core);
GPU_REALW_BUFFER (d_etay_outer_core);
GPU_REALW_BUFFER (d_etaz_outer_core);
GPU_REALW_BUFFER (d_gammax_outer_core);
GPU_REALW_BUFFER (d_gammay_outer_core);
GPU_REALW_BUFFER (d_gammaz_outer_core);

// model parameters
GPU_REALW_BUFFER (d_rhostore_outer_core);
GPU_REALW_BUFFER (d_kappavstore_outer_core);
GPU_REALW_BUFFER (d_rmass_outer_core);
GPU_REALW_BUFFER (d_b_rmass_outer_core);

// global indexing
GPU_INT_BUFFER (d_ibool_outer_core);

// wavefields
GPU_REALW_BUFFER (d_displ_outer_core);
GPU_REALW_BUFFER (d_veloc_outer_core);
GPU_REALW_BUFFER (d_accel_outer_core);
GPU_REALW_BUFFER (d_b_displ_outer_core);
GPU_REALW_BUFFER (d_b_veloc_outer_core);
GPU_REALW_BUFFER (d_b_accel_outer_core);

// kernels
GPU_REALW_BUFFER (d_rho_kl_outer_core);
GPU_REALW_BUFFER (d_alpha_kl_outer_core);

// inner / outer elements
GPU_INT_BUFFER (d_phase_ispec_inner_outer_core);

GPU_INT_BUFFER (d_ibelm_top_outer_core);
GPU_INT_BUFFER (d_ibelm_bottom_outer_core);

// normals definitions for coupling regions
GPU_REALW_BUFFER (d_normal_top_outer_core);
GPU_REALW_BUFFER (d_normal_bottom_outer_core);

// Jacobian definitions
GPU_REALW_BUFFER (d_jacobian2D_top_outer_core);
GPU_REALW_BUFFER (d_jacobian2D_bottom_outer_core);

// ------------------------------------------------------------------   //
// inner_core
// ------------------------------------------------------------------   //
// interpolators
GPU_REALW_BUFFER (d_xix_inner_core);
GPU_REALW_BUFFER (d_xiy_inner_core);
GPU_REALW_BUFFER (d_xiz_inner_core);
GPU_REALW_BUFFER (d_etax_inner_core);
GPU_REALW_BUFFER (d_etay_inner_core);
GPU_REALW_BUFFER (d_etaz_inner_core);
GPU_REALW_BUFFER (d_gammax_inner_core);
GPU_REALW_BUFFER (d_gammay_inner_core);
GPU_REALW_BUFFER (d_gammaz_inner_core);

// model parameters
GPU_REALW_BUFFER (d_rhostore_inner_core);
GPU_REALW_BUFFER (d_kappavstore_inner_core);
GPU_REALW_BUFFER (d_muvstore_inner_core);

GPU_REALW_BUFFER (d_rmassx_inner_core);
GPU_REALW_BUFFER (d_rmassy_inner_core);
GPU_REALW_BUFFER (d_rmassz_inner_core);
GPU_REALW_BUFFER (d_b_rmassx_inner_core);
GPU_REALW_BUFFER (d_b_rmassy_inner_core);
GPU_REALW_BUFFER (d_b_rmassz_inner_core);

// global indexing
GPU_INT_BUFFER (d_ibool_inner_core);
GPU_INT_BUFFER (d_idoubling_inner_core);

// anisotropic 3D mantle
GPU_REALW_BUFFER (d_c11store_inner_core);
GPU_REALW_BUFFER (d_c12store_inner_core);
GPU_REALW_BUFFER (d_c13store_inner_core);
GPU_REALW_BUFFER (d_c33store_inner_core);
GPU_REALW_BUFFER (d_c44store_inner_core);

// wavefields
GPU_REALW_BUFFER (d_displ_inner_core);
GPU_REALW_BUFFER (d_veloc_inner_core);
GPU_REALW_BUFFER (d_accel_inner_core);
GPU_REALW_BUFFER (d_b_displ_inner_core);
GPU_REALW_BUFFER (d_b_veloc_inner_core);
GPU_REALW_BUFFER (d_b_accel_inner_core);

// attenuation
GPU_REALW_BUFFER (d_R_xx_inner_core);
GPU_REALW_BUFFER (d_R_yy_inner_core);
GPU_REALW_BUFFER (d_R_xy_inner_core);
GPU_REALW_BUFFER (d_R_xz_inner_core);
GPU_REALW_BUFFER (d_R_yz_inner_core);
GPU_REALW_BUFFER (d_b_R_xx_inner_core);
GPU_REALW_BUFFER (d_b_R_yy_inner_core);
GPU_REALW_BUFFER (d_b_R_xy_inner_core);
GPU_REALW_BUFFER (d_b_R_xz_inner_core);
GPU_REALW_BUFFER (d_b_R_yz_inner_core);
GPU_REALW_BUFFER (d_factor_common_inner_core);
GPU_REALW_BUFFER (d_one_minus_sum_beta_inner_core);
GPU_REALW_BUFFER (d_epsilondev_xx_inner_core);
GPU_REALW_BUFFER (d_epsilondev_yy_inner_core);
GPU_REALW_BUFFER (d_epsilondev_xy_inner_core);
GPU_REALW_BUFFER (d_epsilondev_xz_inner_core);
GPU_REALW_BUFFER (d_epsilondev_yz_inner_core);
GPU_REALW_BUFFER (d_b_epsilondev_xx_inner_core);
GPU_REALW_BUFFER (d_b_epsilondev_yy_inner_core);
GPU_REALW_BUFFER (d_b_epsilondev_xy_inner_core);
GPU_REALW_BUFFER (d_b_epsilondev_xz_inner_core);
GPU_REALW_BUFFER (d_b_epsilondev_yz_inner_core);
GPU_REALW_BUFFER (d_eps_trace_over_3_inner_core);
GPU_REALW_BUFFER (d_b_eps_trace_over_3_inner_core);

// kernels
GPU_REALW_BUFFER (d_rho_kl_inner_core);
GPU_REALW_BUFFER (d_alpha_kl_inner_core);
GPU_REALW_BUFFER (d_beta_kl_inner_core);

// inner / outer elements
GPU_INT_BUFFER (d_phase_ispec_inner_inner_core);
GPU_INT_BUFFER (d_ibelm_top_inner_core);

// ------------------------------------------------------------------   //
// oceans
// ------------------------------------------------------------------   //
GPU_INT_BUFFER (d_ibool_ocean_load);
GPU_REALW_BUFFER (d_rmass_ocean_load);
GPU_REALW_BUFFER (d_normal_ocean_load);

// ------------------------------------------------------------------   //
// attenuation
// ------------------------------------------------------------------   //
GPU_REALW_BUFFER (d_alphaval);
GPU_REALW_BUFFER (d_betaval);
GPU_REALW_BUFFER (d_gammaval);
GPU_REALW_BUFFER (d_b_alphaval);
GPU_REALW_BUFFER (d_b_betaval);
GPU_REALW_BUFFER (d_b_gammaval);

// ------------------------------------------------------------------   //
// GLL points & weights
// ------------------------------------------------------------------   //
GPU_REALW_BUFFER (d_hprime_xx);
GPU_REALW_BUFFER (d_hprimewgll_xx);
GPU_REALW_BUFFER (d_wgllwgll_xy);
GPU_REALW_BUFFER (d_wgllwgll_xz);
GPU_REALW_BUFFER (d_wgllwgll_yz);
GPU_REALW_BUFFER (d_wgll_cube);

// ------------------------------------------------------------------   //
// gravity
// ------------------------------------------------------------------   //
GPU_REALW_BUFFER (d_gravity_pre_store_outer_core);
GPU_REALW_BUFFER (d_gravity_pre_store_crust_mantle);
GPU_REALW_BUFFER (d_gravity_pre_store_inner_core);
GPU_REALW_BUFFER (d_gravity_H_crust_mantle);
GPU_REALW_BUFFER (d_gravity_H_inner_core);

// ------------------------------------------------------------------   //
// rotation
// ------------------------------------------------------------------   //
GPU_REALW_BUFFER (d_A_array_rotation);
GPU_REALW_BUFFER (d_B_array_rotation);
GPU_REALW_BUFFER (d_b_A_array_rotation);
GPU_REALW_BUFFER (d_b_B_array_rotation);

// ------------------------------------------------------------------   //
// sources
// ------------------------------------------------------------------   //
// full NSOURCES arrays not needed anymore...
//GPU_REALW_BUFFER (d_sourcearrays);
//GPU_DOUBLE_BUFFER (d_stf_pre_compute);
//GPU_INT_BUFFER (d_islice_selected_source);
// for pure adjoint simulation cases
GPU_INT_BUFFER (d_ispec_selected_source);

// local arrays
GPU_INT_BUFFER (d_ispec_selected_source_local);
GPU_REALW_BUFFER (d_sourcearrays_local);
GPU_REALW_BUFFER (d_stf_local);
GPU_REALW_BUFFER (d_b_stf_local);

// ------------------------------------------------------------------   //
// receivers
// ------------------------------------------------------------------   //
GPU_INT_BUFFER (d_number_receiver_global);
GPU_INT_BUFFER (d_ispec_selected_rec);
GPU_INT_BUFFER (d_islice_selected_rec);

GPU_REALW_BUFFER (d_station_seismo_field);
GPU_REALW_BUFFER (d_station_strain_field);
GPU_REALW_BUFFER (d_nu);
GPU_REALW_BUFFER (d_seismograms);

GPU_REALW_BUFFER (d_hxir);
GPU_REALW_BUFFER (d_hetar);
GPU_REALW_BUFFER (d_hgammar);

// adjoint sources
GPU_INT_BUFFER (d_number_adjsources_global);
GPU_REALW_BUFFER (d_stf_array_adjoint);
GPU_REALW_BUFFER (d_hxir_adj);
GPU_REALW_BUFFER (d_hetar_adj);
GPU_REALW_BUFFER (d_hgammar_adj);

GPU_REALW_BUFFER (d_norm_max);
GPU_REALW_BUFFER (d_norm_strain_max);

// ------------------------------------------------------------------   //
// assembly
// ------------------------------------------------------------------   //
GPU_INT_BUFFER (d_nibool_interfaces_crust_mantle);
GPU_INT_BUFFER (d_ibool_interfaces_crust_mantle);
GPU_REALW_BUFFER (d_send_accel_buffer_crust_mantle);
GPU_REALW_BUFFER (d_b_send_accel_buffer_crust_mantle);
GPU_INT_BUFFER (d_nibool_interfaces_inner_core);
GPU_INT_BUFFER (d_ibool_interfaces_inner_core);
GPU_REALW_BUFFER (d_send_accel_buffer_inner_core);
GPU_REALW_BUFFER (d_b_send_accel_buffer_inner_core);
GPU_INT_BUFFER (d_nibool_interfaces_outer_core);
GPU_INT_BUFFER (d_ibool_interfaces_outer_core);
GPU_REALW_BUFFER (d_send_accel_buffer_outer_core);
GPU_REALW_BUFFER (d_b_send_accel_buffer_outer_core);

// ------------------------------------------------------------------   //
// absorbing boundaries
// ------------------------------------------------------------------   //
// crust/mantle
GPU_INT_BUFFER (d_abs_boundary_ispec_crust_mantle);
GPU_INT_BUFFER (d_abs_boundary_npoin_crust_mantle);
GPU_INT_BUFFER (d_abs_boundary_ijk_crust_mantle);
GPU_REALW_BUFFER (d_abs_boundary_jacobian2Dw_crust_mantle);
GPU_REALW_BUFFER (d_abs_boundary_normal_crust_mantle);
GPU_REALW_BUFFER (d_absorb_buffer_crust_mantle);
GPU_REALW_BUFFER (d_rho_vp_crust_mantle);
GPU_REALW_BUFFER (d_rho_vs_crust_mantle);
// outer core
GPU_INT_BUFFER (d_abs_boundary_ispec_outer_core);
GPU_INT_BUFFER (d_abs_boundary_npoin_outer_core);
GPU_INT_BUFFER (d_abs_boundary_ijk_outer_core);
GPU_REALW_BUFFER (d_abs_boundary_jacobian2Dw_outer_core);
GPU_REALW_BUFFER (d_absorb_buffer_outer_core);
GPU_REALW_BUFFER (d_vp_outer_core);

// ------------------------------------------------------------------   //
// noise tomography
// ------------------------------------------------------------------   //
GPU_REALW_BUFFER (d_noise_surface_movie);
GPU_REALW_BUFFER (d_noise_sourcearray);
GPU_REALW_BUFFER (d_normal_x_noise);
GPU_REALW_BUFFER (d_normal_y_noise);
GPU_REALW_BUFFER (d_normal_z_noise);
GPU_REALW_BUFFER (d_mask_noise);
GPU_INT_BUFFER (d_ibelm_top_crust_mantle);
GPU_REALW_BUFFER (d_jacobian2D_top_crust_mantle);
GPU_REALW_BUFFER (d_Sigma_kl);

// ------------------------------------------------------------------   //
// LDDRK
// ------------------------------------------------------------------   //
GPU_REALW_BUFFER (d_tau_sigmainvval);
// wavefields intermediate
GPU_REALW_BUFFER (d_displ_crust_mantle_lddrk);
GPU_REALW_BUFFER (d_veloc_crust_mantle_lddrk);
GPU_REALW_BUFFER (d_b_displ_crust_mantle_lddrk);
GPU_REALW_BUFFER (d_b_veloc_crust_mantle_lddrk);
GPU_REALW_BUFFER (d_displ_outer_core_lddrk);
GPU_REALW_BUFFER (d_veloc_outer_core_lddrk);
GPU_REALW_BUFFER (d_b_displ_outer_core_lddrk);
GPU_REALW_BUFFER (d_b_veloc_outer_core_lddrk);
GPU_REALW_BUFFER (d_displ_inner_core_lddrk);
GPU_REALW_BUFFER (d_veloc_inner_core_lddrk);
GPU_REALW_BUFFER (d_b_displ_inner_core_lddrk);
GPU_REALW_BUFFER (d_b_veloc_inner_core_lddrk);
// rotation
GPU_REALW_BUFFER (d_A_array_rotation_lddrk);
GPU_REALW_BUFFER (d_B_array_rotation_lddrk);
GPU_REALW_BUFFER (d_b_A_array_rotation_lddrk);
GPU_REALW_BUFFER (d_b_B_array_rotation_lddrk);
// attenuation
GPU_REALW_BUFFER (d_R_xx_crust_mantle_lddrk);
GPU_REALW_BUFFER (d_R_yy_crust_mantle_lddrk);
GPU_REALW_BUFFER (d_R_xy_crust_mantle_lddrk);
GPU_REALW_BUFFER (d_R_xz_crust_mantle_lddrk);
GPU_REALW_BUFFER (d_R_yz_crust_mantle_lddrk);
GPU_REALW_BUFFER (d_b_R_xx_crust_mantle_lddrk);
GPU_REALW_BUFFER (d_b_R_yy_crust_mantle_lddrk);
GPU_REALW_BUFFER (d_b_R_xy_crust_mantle_lddrk);
GPU_REALW_BUFFER (d_b_R_xz_crust_mantle_lddrk);
GPU_REALW_BUFFER (d_b_R_yz_crust_mantle_lddrk);
GPU_REALW_BUFFER (d_R_xx_inner_core_lddrk);
GPU_REALW_BUFFER (d_R_yy_inner_core_lddrk);
GPU_REALW_BUFFER (d_R_xy_inner_core_lddrk);
GPU_REALW_BUFFER (d_R_xz_inner_core_lddrk);
GPU_REALW_BUFFER (d_R_yz_inner_core_lddrk);
GPU_REALW_BUFFER (d_b_R_xx_inner_core_lddrk);
GPU_REALW_BUFFER (d_b_R_yy_inner_core_lddrk);
GPU_REALW_BUFFER (d_b_R_xy_inner_core_lddrk);
GPU_REALW_BUFFER (d_b_R_xz_inner_core_lddrk);
GPU_REALW_BUFFER (d_b_R_yz_inner_core_lddrk);
