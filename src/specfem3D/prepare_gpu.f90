!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
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

  subroutine prepare_GPU()

! initializes and copies arrays onto GPU device

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_noise
  use specfem_par_movie

  implicit none

  ! local parameters
  integer :: ier
  real :: free_mb,used_mb,total_mb
  ! dummy CUSTOM_REAL variables to convert from double precision
  real(kind=CUSTOM_REAL),dimension(:),allocatable:: &
    cr_d_ln_density_dr_table,cr_minus_rho_g_over_kappa_fluid, &
    cr_minus_gravity_table,cr_minus_deriv_gravity_table, &
    cr_density_table
  logical :: USE_3D_ATTENUATION_ARRAYS
  real(kind=CUSTOM_REAL) :: dummy

  ! checks if anything to do
  if (.not. GPU_MODE) return

  ! user output
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) "preparing fields and constants on GPU devices"
    call flush_IMAIN()
  endif
  call synchronize_all()

  if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
    USE_3D_ATTENUATION_ARRAYS = .true.
  else
    USE_3D_ATTENUATION_ARRAYS = .false.
  endif

  ! evaluates memory required
  call memory_eval_gpu()

  ! prepares general fields on GPU
  call prepare_constants_device(Mesh_pointer,myrank,NGLLX, &
                                hprime_xx,hprimewgll_xx, &
                                wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                                NSOURCES, nsources_local, &
                                sourcearrays, &
                                islice_selected_source,ispec_selected_source, &
                                nrec, nrec_local, nadj_rec_local, &
                                number_receiver_global, &
                                islice_selected_rec,ispec_selected_rec, &
                                NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
                                NSPEC_CRUST_MANTLE_STRAIN_ONLY, &
                                NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
                                NSPEC_OUTER_CORE,NGLOB_OUTER_CORE, &
                                NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
                                NSPEC_INNER_CORE_STRAIN_ONLY, &
                                SIMULATION_TYPE,NOISE_TOMOGRAPHY, &
                                SAVE_FORWARD,ABSORBING_CONDITIONS, &
                                OCEANS_VAL,GRAVITY_VAL, &
                                ROTATION_VAL,EXACT_MASS_MATRIX_FOR_ROTATION_VAL, &
                                ATTENUATION_VAL,UNDO_ATTENUATION, &
                                PARTIAL_PHYS_DISPERSION_ONLY,USE_3D_ATTENUATION_ARRAYS, &
                                COMPUTE_AND_STORE_STRAIN, &
                                ANISOTROPIC_3D_MANTLE_VAL,ANISOTROPIC_INNER_CORE_VAL, &
                                SAVE_BOUNDARY_MESH, &
                                USE_MESH_COLORING_GPU, &
                                ANISOTROPIC_KL,APPROXIMATE_HESS_KL, &
                                deltat,b_deltat, &
                                GPU_ASYNC_COPY, &
                                hxir_store,hetar_store,hgammar_store,nu)
  call synchronize_all()

  ! prepares rotation arrays
  if (ROTATION_VAL) then
    if (myrank == 0) then
      write(IMAIN,*) "  loading rotation arrays"
      call flush_IMAIN()
    endif
    call synchronize_all()

    if (SIMULATION_TYPE == 3) then
      call prepare_fields_rotation_device(Mesh_pointer, &
                                          two_omega_earth, &
                                          A_array_rotation,B_array_rotation, &
                                          b_two_omega_earth, &
                                          b_A_array_rotation,b_B_array_rotation, &
                                          NSPEC_OUTER_CORE_ROTATION)
    else
      call prepare_fields_rotation_device(Mesh_pointer, &
                                          two_omega_earth, &
                                          A_array_rotation,B_array_rotation, &
                                          dummy, &
                                          dummy,dummy, &
                                          NSPEC_OUTER_CORE_ROTATION)
    endif
  endif
  call synchronize_all()

  ! prepares arrays related to gravity
  ! note: GPU will use only single-precision (or double precision) for all calculations
  !          using implicit conversion
  if (myrank == 0) then
    write(IMAIN,*) "  loading non-gravity/gravity arrays"
    call flush_IMAIN()
  endif
  call synchronize_all()

  allocate(cr_d_ln_density_dr_table(NRAD_GRAVITY), &
           cr_minus_rho_g_over_kappa_fluid(NRAD_GRAVITY), &
           cr_minus_gravity_table(NRAD_GRAVITY), &
           cr_minus_deriv_gravity_table(NRAD_GRAVITY), &
           cr_density_table(NRAD_GRAVITY), &
           stat=ier)
  if (ier /= 0 ) stop 'Error allocating cr_minus_rho_g_over_kappa_fluid, etc...'

  ! d_ln_density_dr_table needed for no gravity case
  cr_d_ln_density_dr_table(:) = real(d_ln_density_dr_table(:), kind=CUSTOM_REAL)
  ! these are needed for gravity cases only
  cr_minus_rho_g_over_kappa_fluid(:) = real(minus_rho_g_over_kappa_fluid(:), kind=CUSTOM_REAL)
  cr_minus_gravity_table(:) = real(minus_gravity_table(:), kind=CUSTOM_REAL)
  cr_minus_deriv_gravity_table(:) = real(minus_deriv_gravity_table(:), kind=CUSTOM_REAL)
  cr_density_table(:) = real(density_table(:), kind=CUSTOM_REAL)

  ! prepares on GPU
  call prepare_fields_gravity_device(Mesh_pointer, &
                                     cr_d_ln_density_dr_table, &
                                     cr_minus_rho_g_over_kappa_fluid, &
                                     cr_minus_gravity_table, &
                                     cr_minus_deriv_gravity_table, &
                                     cr_density_table, &
                                     wgll_cube, &
                                     NRAD_GRAVITY, &
                                     minus_g_icb,minus_g_cmb, &
                                     RHO_BOTTOM_OC,RHO_TOP_OC)

  deallocate(cr_d_ln_density_dr_table,cr_minus_rho_g_over_kappa_fluid, &
             cr_minus_gravity_table,cr_minus_deriv_gravity_table, &
             cr_density_table)
  call synchronize_all()

  ! prepares attenuation arrays
  if (ATTENUATION_VAL) then
    if (myrank == 0) then
      write(IMAIN,*) "  loading attenuation"
      call flush_IMAIN()
    endif
    call synchronize_all()

    if (SIMULATION_TYPE == 3) then
      call prepare_fields_attenuat_device(Mesh_pointer, &
                                          R_xx_crust_mantle,R_yy_crust_mantle,R_xy_crust_mantle, &
                                          R_xz_crust_mantle,R_yz_crust_mantle, &
                                          b_R_xx_crust_mantle,b_R_yy_crust_mantle,b_R_xy_crust_mantle, &
                                          b_R_xz_crust_mantle,b_R_yz_crust_mantle, &
                                          factor_common_crust_mantle, &
                                          one_minus_sum_beta_crust_mantle, &
                                          R_xx_inner_core,R_yy_inner_core,R_xy_inner_core, &
                                          R_xz_inner_core,R_yz_inner_core, &
                                          b_R_xx_inner_core,b_R_yy_inner_core,b_R_xy_inner_core, &
                                          b_R_xz_inner_core,b_R_yz_inner_core, &
                                          factor_common_inner_core, &
                                          one_minus_sum_beta_inner_core, &
                                          alphaval,betaval,gammaval, &
                                          b_alphaval,b_betaval,b_gammaval)
    else
      call prepare_fields_attenuat_device(Mesh_pointer, &
                                          R_xx_crust_mantle,R_yy_crust_mantle,R_xy_crust_mantle, &
                                          R_xz_crust_mantle,R_yz_crust_mantle, &
                                          dummy,dummy,dummy, &
                                          dummy,dummy, &
                                          factor_common_crust_mantle, &
                                          one_minus_sum_beta_crust_mantle, &
                                          R_xx_inner_core,R_yy_inner_core,R_xy_inner_core, &
                                          R_xz_inner_core,R_yz_inner_core, &
                                          dummy,dummy,dummy, &
                                          dummy,dummy, &
                                          factor_common_inner_core, &
                                          one_minus_sum_beta_inner_core, &
                                          alphaval,betaval,gammaval, &
                                          dummy,dummy,dummy)

    endif
  endif
  call synchronize_all()


  ! prepares attenuation arrays
  if (COMPUTE_AND_STORE_STRAIN) then
    if (myrank == 0) then
      write(IMAIN,*) "  loading strain"
      call flush_IMAIN()
    endif
    call synchronize_all()

    if (SIMULATION_TYPE == 3) then
      call prepare_fields_strain_device(Mesh_pointer, &
                                        epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                                        epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                                        b_epsilondev_xx_crust_mantle,b_epsilondev_yy_crust_mantle,b_epsilondev_xy_crust_mantle, &
                                        b_epsilondev_xz_crust_mantle,b_epsilondev_yz_crust_mantle, &
                                        eps_trace_over_3_crust_mantle, &
                                        b_eps_trace_over_3_crust_mantle, &
                                        epsilondev_xx_inner_core,epsilondev_yy_inner_core,epsilondev_xy_inner_core, &
                                        epsilondev_xz_inner_core,epsilondev_yz_inner_core, &
                                        b_epsilondev_xx_inner_core,b_epsilondev_yy_inner_core,b_epsilondev_xy_inner_core, &
                                        b_epsilondev_xz_inner_core,b_epsilondev_yz_inner_core, &
                                        eps_trace_over_3_inner_core, &
                                        b_eps_trace_over_3_inner_core)
    else
      call prepare_fields_strain_device(Mesh_pointer, &
                                        epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                                        epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                                        dummy,dummy,dummy, &
                                        dummy,dummy, &
                                        eps_trace_over_3_crust_mantle, &
                                        dummy, &
                                        epsilondev_xx_inner_core,epsilondev_yy_inner_core,epsilondev_xy_inner_core, &
                                        epsilondev_xz_inner_core,epsilondev_yz_inner_core, &
                                        dummy,dummy,dummy, &
                                        dummy,dummy, &
                                        eps_trace_over_3_inner_core, &
                                        dummy)

    endif
  endif
  call synchronize_all()

  ! prepares absorbing arrays
  if (NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) then
    if (myrank == 0) then
      write(IMAIN,*) "  loading absorbing boundaries"
      call flush_IMAIN()
    endif
    call synchronize_all()

    call prepare_fields_absorb_device(Mesh_pointer, &
                                      nspec2D_xmin_crust_mantle,nspec2D_xmax_crust_mantle, &
                                      nspec2D_ymin_crust_mantle,nspec2D_ymax_crust_mantle, &
                                      NSPEC2DMAX_XMIN_XMAX_CM,NSPEC2DMAX_YMIN_YMAX_CM, &
                                      nimin_crust_mantle,nimax_crust_mantle, &
                                      njmin_crust_mantle,njmax_crust_mantle, &
                                      nkmin_xi_crust_mantle,nkmin_eta_crust_mantle, &
                                      ibelm_xmin_crust_mantle,ibelm_xmax_crust_mantle, &
                                      ibelm_ymin_crust_mantle,ibelm_ymax_crust_mantle, &
                                      normal_xmin_crust_mantle,normal_xmax_crust_mantle, &
                                      normal_ymin_crust_mantle,normal_ymax_crust_mantle, &
                                      jacobian2D_xmin_crust_mantle,jacobian2D_xmax_crust_mantle, &
                                      jacobian2D_ymin_crust_mantle,jacobian2D_ymax_crust_mantle, &
                                      rho_vp_crust_mantle,rho_vs_crust_mantle, &
                                      nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
                                      nspec2D_ymin_outer_core,nspec2D_ymax_outer_core, &
                                      nspec2D_zmin_outer_core, &
                                      NSPEC2DMAX_XMIN_XMAX_OC,NSPEC2DMAX_YMIN_YMAX_OC, &
                                      nimin_outer_core,nimax_outer_core, &
                                      njmin_outer_core,njmax_outer_core, &
                                      nkmin_xi_outer_core,nkmin_eta_outer_core, &
                                      ibelm_xmin_outer_core,ibelm_xmax_outer_core, &
                                      ibelm_ymin_outer_core,ibelm_ymax_outer_core, &
                                      jacobian2D_xmin_outer_core,jacobian2D_xmax_outer_core, &
                                      jacobian2D_ymin_outer_core,jacobian2D_ymax_outer_core, &
                                      vp_outer_core)

  endif
  call synchronize_all()

  ! prepares MPI interfaces
  if (myrank == 0) then
    write(IMAIN,*) "  loading MPI interfaces"
    call flush_IMAIN()
  endif
  call synchronize_all()

  call prepare_mpi_buffers_device(Mesh_pointer, &
                                  num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                                  nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
                                  num_interfaces_inner_core,max_nibool_interfaces_ic, &
                                  nibool_interfaces_inner_core,ibool_interfaces_inner_core, &
                                  num_interfaces_outer_core,max_nibool_interfaces_oc, &
                                  nibool_interfaces_outer_core,ibool_interfaces_outer_core)
  call synchronize_all()

  ! prepares fields on GPU for noise simulations
  if (NOISE_TOMOGRAPHY > 0) then
    if (myrank == 0) then
      write(IMAIN,*) "  loading noise arrays"
      call flush_IMAIN()
    endif
    call synchronize_all()

    call prepare_fields_noise_device(Mesh_pointer,NSPEC_TOP, &
                                     NSTEP, &
                                     ibelm_top_crust_mantle, &
                                     noise_sourcearray, &
                                     normal_x_noise,normal_y_noise,normal_z_noise, &
                                     mask_noise,jacobian2D_top_crust_mantle)

  endif
  call synchronize_all()

  ! prepares oceans arrays
  if (OCEANS_VAL) then
    if (myrank == 0) then
      write(IMAIN,*) "  loading oceans arrays"
      call flush_IMAIN()
    endif
    call synchronize_all()

    call prepare_oceans_device(Mesh_pointer,npoin_oceans, &
                              ibool_ocean_load, &
                              rmass_ocean_load_selected, &
                              normal_ocean_load)
  endif
  call synchronize_all()

  ! prepares LDDRK arrays
  if (USE_LDDRK) then
    if (myrank == 0) then
      write(IMAIN,*) "  loading LDDRK arrays"
      call flush_IMAIN()
    endif
    call synchronize_all()

    stop 'prepare_lddrk_device not implemented yet'
    !call prepare_lddrk_device(Mesh_pointer)
  endif

  ! crust/mantle region
  if (myrank == 0) then
    write(IMAIN,*) "  loading crust/mantle region"
    call flush_IMAIN()
  endif
  call synchronize_all()

  if (SIMULATION_TYPE == 3) then
    call prepare_crust_mantle_device(Mesh_pointer, &
                                     xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                     etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                     gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                     rhostore_crust_mantle, &
                                     kappavstore_crust_mantle,muvstore_crust_mantle, &
                                     kappahstore_crust_mantle,muhstore_crust_mantle, &
                                     eta_anisostore_crust_mantle, &
                                     rmassx_crust_mantle,rmassy_crust_mantle,rmassz_crust_mantle, &
                                     b_rmassx_crust_mantle,b_rmassy_crust_mantle, &
                                     ibool_crust_mantle, &
                                     rstore_crust_mantle, &
                                     ispec_is_tiso_crust_mantle, &
                                     c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
                                     c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
                                     c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
                                     c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
                                     c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
                                     c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                                     c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                                     num_phase_ispec_crust_mantle,phase_ispec_inner_crust_mantle, &
                                     nspec_outer_crust_mantle,nspec_inner_crust_mantle, &
                                     NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE), &
                                     ibelm_bottom_crust_mantle, &
                                     NCHUNKS_VAL, &
                                     num_colors_outer_crust_mantle,num_colors_inner_crust_mantle, &
                                     num_elem_colors_crust_mantle)
  else
    call prepare_crust_mantle_device(Mesh_pointer, &
                                     xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                     etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                     gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                     rhostore_crust_mantle, &
                                     kappavstore_crust_mantle,muvstore_crust_mantle, &
                                     kappahstore_crust_mantle,muhstore_crust_mantle, &
                                     eta_anisostore_crust_mantle, &
                                     rmassx_crust_mantle,rmassy_crust_mantle,rmassz_crust_mantle, &
                                     dummy,dummy, &
                                     ibool_crust_mantle, &
                                     rstore_crust_mantle, &
                                     ispec_is_tiso_crust_mantle, &
                                     c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
                                     c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
                                     c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
                                     c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
                                     c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
                                     c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                                     c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                                     num_phase_ispec_crust_mantle,phase_ispec_inner_crust_mantle, &
                                     nspec_outer_crust_mantle,nspec_inner_crust_mantle, &
                                     NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE), &
                                     ibelm_bottom_crust_mantle, &
                                     NCHUNKS_VAL, &
                                     num_colors_outer_crust_mantle,num_colors_inner_crust_mantle, &
                                     num_elem_colors_crust_mantle)
  endif
  call synchronize_all()

  ! outer core region
  if (myrank == 0) then
    write(IMAIN,*) "  loading outer core region"
    call flush_IMAIN()
  endif
  call synchronize_all()

  call prepare_outer_core_device(Mesh_pointer, &
                                 xix_outer_core,xiy_outer_core,xiz_outer_core, &
                                 etax_outer_core,etay_outer_core,etaz_outer_core, &
                                 gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
                                 rhostore_outer_core,kappavstore_outer_core, &
                                 rmass_outer_core, &
                                 ibool_outer_core, &
                                 rstore_outer_core, &
                                 num_phase_ispec_outer_core,phase_ispec_inner_outer_core, &
                                 nspec_outer_outer_core,nspec_inner_outer_core, &
                                 NSPEC2D_TOP(IREGION_OUTER_CORE), &
                                 NSPEC2D_BOTTOM(IREGION_OUTER_CORE), &
                                 normal_top_outer_core, &
                                 normal_bottom_outer_core, &
                                 jacobian2D_top_outer_core, &
                                 jacobian2D_bottom_outer_core, &
                                 ibelm_top_outer_core, &
                                 ibelm_bottom_outer_core, &
                                 num_colors_outer_outer_core,num_colors_inner_outer_core, &
                                 num_elem_colors_outer_core)
  call synchronize_all()

  ! inner core region
  if (myrank == 0) then
    write(IMAIN,*) "  loading inner core region"
    call flush_IMAIN()
  endif
  call synchronize_all()

  if (SIMULATION_TYPE == 3) then
    call prepare_inner_core_device(Mesh_pointer, &
                                   xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                   etax_inner_core,etay_inner_core,etaz_inner_core, &
                                   gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                                   rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core, &
                                   rmassx_inner_core,rmassy_inner_core,rmassz_inner_core, &
                                   b_rmassx_inner_core,b_rmassy_inner_core, &
                                   ibool_inner_core, &
                                   rstore_inner_core, &
                                   c11store_inner_core,c12store_inner_core,c13store_inner_core, &
                                   c33store_inner_core,c44store_inner_core, &
                                   idoubling_inner_core, &
                                   num_phase_ispec_inner_core,phase_ispec_inner_inner_core, &
                                   nspec_outer_inner_core,nspec_inner_inner_core, &
                                   NSPEC2D_TOP(IREGION_INNER_CORE), &
                                   ibelm_top_inner_core, &
                                   num_colors_outer_inner_core,num_colors_inner_inner_core, &
                                   num_elem_colors_inner_core)
  else
    call prepare_inner_core_device(Mesh_pointer, &
                                   xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                   etax_inner_core,etay_inner_core,etaz_inner_core, &
                                   gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                                   rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core, &
                                   rmassx_inner_core,rmassy_inner_core,rmassz_inner_core, &
                                   dummy,dummy, &
                                   ibool_inner_core, &
                                   rstore_inner_core, &
                                   c11store_inner_core,c12store_inner_core,c13store_inner_core, &
                                   c33store_inner_core,c44store_inner_core, &
                                   idoubling_inner_core, &
                                   num_phase_ispec_inner_core,phase_ispec_inner_inner_core, &
                                   nspec_outer_inner_core,nspec_inner_inner_core, &
                                   NSPEC2D_TOP(IREGION_INNER_CORE), &
                                   ibelm_top_inner_core, &
                                   num_colors_outer_inner_core,num_colors_inner_inner_core, &
                                   num_elem_colors_inner_core)
  endif
  call synchronize_all()

  ! transfer forward and backward fields to device with initial values
  if (myrank == 0) then
    write(IMAIN,*) "  transferring initial wavefield"
    call flush_IMAIN()
  endif
  call synchronize_all()

  call transfer_fields_cm_to_device(NDIM*NGLOB_CRUST_MANTLE,displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle, &
                                    Mesh_pointer)

  call transfer_fields_ic_to_device(NDIM*NGLOB_INNER_CORE,displ_inner_core,veloc_inner_core,accel_inner_core, &
                                    Mesh_pointer)

  call transfer_fields_oc_to_device(NGLOB_OUTER_CORE,displ_outer_core,veloc_outer_core,accel_outer_core, &
                                    Mesh_pointer)

  if (SIMULATION_TYPE == 3) then
    call transfer_b_fields_cm_to_device(NDIM*NGLOB_CRUST_MANTLE,b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle, &
                                        Mesh_pointer)

    call transfer_b_fields_ic_to_device(NDIM*NGLOB_INNER_CORE,b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core, &
                                        Mesh_pointer)

    call transfer_b_fields_oc_to_device(NGLOB_OUTER_CORE,b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core, &
                                        Mesh_pointer)
  endif

  ! outputs GPU usage to files for all processes
  call output_free_device_memory(myrank)

  ! outputs usage for main process
  if (myrank == 0) then
    ! gets memory usage for main process
    call get_free_device_memory(free_mb,used_mb,total_mb)

    ! outputs info
    if (total_mb /= 0) then
       write(IMAIN,*)
       write(IMAIN,*) "  GPU usage: free  =",free_mb," MB",nint(free_mb/total_mb*100.0),"%"
       write(IMAIN,*) "             used  =",used_mb," MB",nint(used_mb/total_mb*100.0),"%"
       write(IMAIN,*) "             total =",total_mb," MB",nint(total_mb/total_mb*100.0),"%"
       write(IMAIN,*)
    else
       write(IMAIN,*)
       write(IMAIN,*) "  GPU usage: not available."
       write(IMAIN,*)
    endif
    call flush_IMAIN()
  endif

  ! synchronizes processes
  call synchronize_all()

  contains

    subroutine memory_eval_gpu()

    implicit none

    ! local parameters
    double precision :: memory_size
    integer,parameter :: NGLL2 = 25
    integer,parameter :: NGLL3 = 125
    integer,parameter :: NGLL3_PADDED = 128
    integer :: NSPEC_AB,NGLOB_AB
    integer :: NGLL_ATT,NSPEC_2

    memory_size = 0.d0

    ! crust/mantle + inner core
    NSPEC_AB = NSPEC_CRUST_MANTLE + NSPEC_INNER_CORE
    NGLOB_AB = NGLOB_CRUST_MANTLE + NGLOB_INNER_CORE

    ! add size of each set of arrays multiplied by the number of such arrays
    ! d_hprime_xx,d_hprimewgll_xx
    memory_size = memory_size + 2.d0 * NGLL2 * dble(CUSTOM_REAL)

    ! sources
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      ! d_sourcearrays
      memory_size = memory_size + NGLL3 * NSOURCES * NDIM * dble(CUSTOM_REAL)
    endif
    ! d_islice_selected_source,d_ispec_selected_source
    memory_size = memory_size + 2.0 * NSOURCES * dble(SIZE_INTEGER)

    ! receivers
    !d_number_receiver_global
    memory_size = memory_size + nrec_local * dble(SIZE_INTEGER)
    ! d_station_seismo_field
    memory_size = memory_size + NDIM * NGLL3 * nrec_local * dble(CUSTOM_REAL)
    ! d_station_strain_field
    if (SIMULATION_TYPE == 2) then
      memory_size = memory_size + NGLL3 * nrec_local * dble(SIZE_INTEGER)
    endif
    ! d_ispec_selected_rec
    memory_size = memory_size + nrec * dble(SIZE_INTEGER)

    ! d_adj_source_adjoint
    memory_size = memory_size + NDIM *  nadj_rec_local * dble(CUSTOM_REAL)

    ! rotation
    if (ROTATION_VAL) then
      ! d_A_array_rotation,..
      memory_size = memory_size + 2.d0 * NGLL3 * NSPEC_OUTER_CORE * dble(CUSTOM_REAL)
    endif

    ! gravity
    if (GRAVITY_VAL) then
      ! d_minus_rho_g_over_kappa_fluid,..
      memory_size = memory_size + 4.d0 * NRAD_GRAVITY * dble(CUSTOM_REAL)
    else
      ! d_d_ln_density_dr_table
      memory_size = memory_size + NRAD_GRAVITY * dble(CUSTOM_REAL)
    endif

    ! attenuation
    if (ATTENUATION_VAL) then
      if (USE_3D_ATTENUATION_ARRAYS) then
        NGLL_ATT = NGLL3
      else
        NGLL_ATT = 1
      endif
      ! d_one_minus_sum_beta_crust_mantle,..
      memory_size = memory_size + NGLL_ATT * NSPEC_AB * dble(CUSTOM_REAL)

      if (.not. PARTIAL_PHYS_DISPERSION_ONLY) then
        ! d_factor_common_crust_mantle,..
        memory_size = memory_size + N_SLS * NGLL_ATT * NSPEC_AB * dble(CUSTOM_REAL)
        ! d_R_xx_crust_mantle,..
        memory_size = memory_size + 5.d0 * N_SLS * NGLL3 * NSPEC_AB * dble(CUSTOM_REAL)
      endif
      ! alphaval,..
      memory_size = memory_size + 3.d0 * N_SLS * dble(CUSTOM_REAL)
    endif

    ! strains
    if (COMPUTE_AND_STORE_STRAIN) then
      ! d_epsilondev_xx_crust_mantle,..
      memory_size = memory_size + 5.d0 * NGLL3 * NSPEC_AB * dble(CUSTOM_REAL)
      ! d_eps_trace_over_3_crust_mantle,..
      NSPEC_2 = NSPEC_CRUST_MANTLE_STRAIN_ONLY + NSPEC_INNER_CORE_STRAIN_ONLY
      memory_size = memory_size + NGLL3 * NSPEC_2 * dble(CUSTOM_REAL)
    endif

    ! absorbing boundaries
    if (NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) then
      ! d_rho_vp_crust_mantle,..
      memory_size = memory_size + 2.d0 * NGLL3 * NSPEC_CRUST_MANTLE * dble(CUSTOM_REAL)
      ! d_nkmin_xi_crust_mantle,..,d_nkmin_eta_crust_mantle,..
      NSPEC_2 = NSPEC2DMAX_XMIN_XMAX_CM + NSPEC2DMAX_YMIN_YMAX_CM
      memory_size = memory_size + 6.d0 * NSPEC_2 * dble(SIZE_INTEGER)
      ! d_ibelm_xmin_crust_mantle,..
      NSPEC_2 = nspec2D_xmin_crust_mantle + nspec2D_xmax_crust_mantle
      memory_size = memory_size + NSPEC_2 * dble(SIZE_INTEGER)
      ! d_normal_xmax_crust_mantle,..
      memory_size = memory_size + NDIM * NGLL2 * NSPEC_2 * dble(CUSTOM_REAL)
      ! d_jacobian2D_xmax_crust_mantle,..
      memory_size = memory_size + NGLL2 * NSPEC_2 * dble(CUSTOM_REAL)
      ! d_ibelm_ymin_crust_mantle
      NSPEC_2 = nspec2D_ymin_crust_mantle + nspec2D_ymax_crust_mantle
      memory_size = memory_size + NSPEC_2 * dble(SIZE_INTEGER)
      ! d_normal_ymin_crust_mantle,..
      memory_size = memory_size + NDIM * NGLL2 * NSPEC_2 * dble(CUSTOM_REAL)
      ! d_jacobian2D_ymin_crust_mantle,..
      memory_size = memory_size + NGLL2 * NSPEC_2 * dble(CUSTOM_REAL)
      ! d_vp_outer_core
      memory_size = memory_size + NGLL3 * NSPEC_OUTER_CORE * dble(CUSTOM_REAL)
      ! d_nkmin_xi_outer_core,..,d_nkmin_eta_outer_core,..
      NSPEC_2 = NSPEC2DMAX_XMIN_XMAX_OC + NSPEC2DMAX_YMIN_YMAX_OC
      memory_size = memory_size + 6.d0 * NSPEC_2 * dble(SIZE_INTEGER)
    endif

    ! MPI buffers
    ! d_ibool_interfaces_crust_mantle
    memory_size = memory_size + num_interfaces_crust_mantle * max_nibool_interfaces_cm * dble(SIZE_INTEGER)
    ! d_send_accel_buffer_crust_mantle
    memory_size = memory_size + NDIM * num_interfaces_crust_mantle * max_nibool_interfaces_cm * dble(CUSTOM_REAL)
    ! d_ibool_interfaces_inner_core
    memory_size = memory_size + num_interfaces_inner_core * max_nibool_interfaces_ic * dble(SIZE_INTEGER)
    ! d_send_accel_buffer_inner_core
    memory_size = memory_size + NDIM * num_interfaces_inner_core * max_nibool_interfaces_ic * dble(CUSTOM_REAL)
    ! d_ibool_interfaces_outer_core
    memory_size = memory_size + num_interfaces_outer_core * max_nibool_interfaces_oc * dble(SIZE_INTEGER)
    ! d_send_accel_buffer_outer_core
    memory_size = memory_size + num_interfaces_outer_core * max_nibool_interfaces_oc * dble(CUSTOM_REAL)

    ! noise
    if (NOISE_TOMOGRAPHY > 0) then
      ! d_noise_surface_movie
      memory_size = memory_size + NDIM * NGLL2 * NSPEC_TOP * dble(CUSTOM_REAL)
      ! d_noise_sourcearray
      if (NOISE_TOMOGRAPHY == 1) then
        memory_size = memory_size + NDIM * NGLL3 * NSTEP * dble(CUSTOM_REAL)
      endif
      ! d_normal_x_noise,..
      if (NOISE_TOMOGRAPHY == 1) then
        memory_size = memory_size + 4.d0 * NGLL2 * NSPEC_TOP * dble(CUSTOM_REAL)
      endif
    endif

    ! oceans
    if (OCEANS_VAL) then
      ! d_ibool_ocean_load
      memory_size = memory_size + npoin_oceans * dble(SIZE_INTEGER)
      ! d_rmass_ocean_load,..
      memory_size = memory_size + 4.d0 * npoin_oceans * dble(CUSTOM_REAL)
    endif

    ! crust/mantle + inner core
    ! padded xix,..gammaz
    memory_size = memory_size + 9.d0 * NGLL3_PADDED * NSPEC_AB * dble(CUSTOM_REAL)
    ! ibool
    memory_size = memory_size + NGLL3 * NSPEC_AB * dble(SIZE_INTEGER)
    ! crust/mantle
    if (.not. ANISOTROPIC_3D_MANTLE_VAL) then
      ! padded kappav,kappah,..
      memory_size = memory_size + 5.d0 * NGLL3_PADDED * NSPEC_CRUST_MANTLE * dble(CUSTOM_REAL)
    else
      ! padded c11,..
      memory_size = memory_size + 21.d0 * NGLL3_PADDED * NSPEC_CRUST_MANTLE * dble(CUSTOM_REAL)
    endif
    ! rstore
    memory_size = memory_size + 3.d0 * NGLOB_CRUST_MANTLE * dble(CUSTOM_REAL)
    ! inner core
    ! padded muv
    memory_size = memory_size + NGLL3_PADDED * NSPEC_INNER_CORE * dble(CUSTOM_REAL)
    if (.not. ANISOTROPIC_INNER_CORE_VAL) then
      ! padded kappav
      memory_size = memory_size + NGLL3_PADDED * NSPEC_INNER_CORE * dble(CUSTOM_REAL)
    else
      ! padded c11,..
      memory_size = memory_size + 5.d0 * NGLL3_PADDED * NSPEC_INNER_CORE * dble(CUSTOM_REAL)
    endif
    ! rstore
    if (GRAVITY_VAL) then
      memory_size = memory_size + 3.d0 * NGLOB_INNER_CORE * dble(CUSTOM_REAL)
    endif
    ! d_phase_ispec_inner_crust_mantle
    memory_size = memory_size + 2 * num_phase_ispec_crust_mantle * dble(SIZE_INTEGER)
    ! d_displ,..
    memory_size = memory_size + 3.d0 * NDIM * NGLOB_AB * dble(CUSTOM_REAL)
    ! crust/mantle
    ! d_rmassz
    memory_size = memory_size + NGLOB_AB * dble(CUSTOM_REAL)
    ! d_rmassx,..
    if ((NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) .or. (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL)) then
      memory_size = memory_size + 2.d0 * NGLOB_CRUST_MANTLE * dble(CUSTOM_REAL)
    endif
    ! inner core
    if ((ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL)) then
      ! d_rmassx,..
      memory_size = memory_size + 2.d0 * NGLOB_INNER_CORE * dble(CUSTOM_REAL)
    endif

    ! outer core
    ! padded d_xix_outer_core,..
    memory_size = memory_size + 9.d0 * NGLL3_PADDED * NSPEC_OUTER_CORE * dble(CUSTOM_REAL)
    ! padded d_kappav
    memory_size = memory_size + NGLL3_PADDED * NSPEC_OUTER_CORE * dble(CUSTOM_REAL)
    ! ibool
    memory_size = memory_size + NGLL3 * NSPEC_OUTER_CORE * dble(SIZE_INTEGER)
    ! d_rstore_outer_core,..
    memory_size = memory_size + 3.d0 * NGLOB_OUTER_CORE * dble(CUSTOM_REAL)
    ! d_phase_ispec_inner_outer_core
    memory_size = memory_size + 2.d0 * num_phase_ispec_outer_core * dble(SIZE_INTEGER)
    ! d_displ_outer,..
    memory_size = memory_size + 3.d0 * NGLOB_OUTER_CORE * dble(CUSTOM_REAL)
    ! d_rmass_outer_core
    memory_size = memory_size + NGLOB_OUTER_CORE * dble(CUSTOM_REAL)

    ! poor estimate for kernel simulations...
    if (SIMULATION_TYPE == 3) memory_size = 2.d0 * memory_size

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) "  minimum memory requested     : ",memory_size / 1024. / 1024.,"MB per process"
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    end subroutine memory_eval_gpu

  end subroutine prepare_GPU
