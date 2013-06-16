
!! DK DK for kernel calculations


! kernel calculations
  if (SIMULATION_TYPE == 3) then
    ! crust mantle
    call compute_kernels_crust_mantle(ibool_crust_mantle, &
                          rho_kl_crust_mantle,beta_kl_crust_mantle, &
                          alpha_kl_crust_mantle,cijkl_kl_crust_mantle, &
                          accel_crust_mantle,b_displ_crust_mantle, &
                          deltat,displ_crust_mantle,hprime_xx,hprime_xxT,&
                          xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle,&
                          etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
                          gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle)

    ! outer core
    call compute_kernels_outer_core(ibool_outer_core, &
                        xix_outer_core,xiy_outer_core,xiz_outer_core, &
                        etax_outer_core,etay_outer_core,etaz_outer_core, &
                        gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        displ_outer_core,accel_outer_core, &
                        b_displ_outer_core,b_accel_outer_core, &
                        vector_accel_outer_core,vector_displ_outer_core, &
                        b_vector_displ_outer_core, &
                        div_displ_outer_core,b_div_displ_outer_core, &
                        rhostore_outer_core,kappavstore_outer_core, &
                        rho_kl_outer_core,alpha_kl_outer_core, &
                        deviatoric_outercore,nspec_beta_kl_outer_core,beta_kl_outer_core, &
                        deltat)

    ! inner core
    call compute_kernels_inner_core(ibool_inner_core, &
                          rho_kl_inner_core,beta_kl_inner_core, &
                          alpha_kl_inner_core, &
                          accel_inner_core,b_displ_inner_core, &
                          deltat,displ_inner_core,hprime_xx,hprime_xxT,&
                          xix_inner_core,xiy_inner_core,xiz_inner_core,&
                          etax_inner_core,etay_inner_core,etaz_inner_core,&
                          gammax_inner_core,gammay_inner_core,gammaz_inner_core)

    ! NOISE TOMOGRAPHY --- source strength kernel
    if (NOISE_TOMOGRAPHY == 3)  &
       call compute_kernels_strength_noise(nmovie_points,ibool_crust_mantle, &
                          Sigma_kl_crust_mantle,displ_crust_mantle,deltat,it, &
                          normal_x_noise,normal_y_noise,normal_z_noise, &
                          NSPEC2D_TOP(IREGION_CRUST_MANTLE),noise_surface_movie, &
                          ibelm_top_crust_mantle)

    ! --- boundary kernels ------
    if (SAVE_BOUNDARY_MESH) then
      fluid_solid_boundary = .false.
      iregion_code = IREGION_CRUST_MANTLE

      ! Moho
      if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
        call compute_boundary_kernel(displ_crust_mantle,accel_crust_mantle, &
                 b_displ_crust_mantle,nspec_crust_mantle,iregion_code, &
                 ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
                 xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                 etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
                 gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,hprime_xx,hprime_yy,hprime_zz, &
                 rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
                 kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
                 c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
                 c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
                 c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
                 c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
                 c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                 c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                 k_top,ibelm_moho_top,normal_moho,moho_kl_top,fluid_solid_boundary,NSPEC2D_MOHO)

        call compute_boundary_kernel(displ_crust_mantle,accel_crust_mantle, &
                 b_displ_crust_mantle,nspec_crust_mantle,iregion_code, &
                 ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
                 xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                 etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
                 gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,hprime_xx,hprime_yy,hprime_zz, &
                 rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
                 kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
                 c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
                 c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
                 c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
                 c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
                 c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                 c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                 k_bot,ibelm_moho_bot,normal_moho,moho_kl_bot,fluid_solid_boundary,NSPEC2D_MOHO)

        moho_kl = moho_kl + (moho_kl_top - moho_kl_bot) * deltat
      endif

      ! d400
      call compute_boundary_kernel(displ_crust_mantle,accel_crust_mantle, &
                 b_displ_crust_mantle,nspec_crust_mantle,iregion_code, &
                 ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
                 xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                 etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
                 gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,hprime_xx,hprime_yy,hprime_zz, &
                 rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
                 kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
                 c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
                 c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
                 c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
                 c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
                 c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                 c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                 k_top,ibelm_400_top,normal_400,d400_kl_top,fluid_solid_boundary,NSPEC2D_400)

      call compute_boundary_kernel(displ_crust_mantle,accel_crust_mantle, &
                 b_displ_crust_mantle,nspec_crust_mantle,iregion_code, &
                 ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
                 xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                 etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
                 gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,hprime_xx,hprime_yy,hprime_zz, &
                 rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
                 kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
                 c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
                 c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
                 c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
                 c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
                 c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                 c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                 k_bot,ibelm_400_bot,normal_400,d400_kl_bot,fluid_solid_boundary,NSPEC2D_400)

      d400_kl = d400_kl + (d400_kl_top - d400_kl_bot) * deltat

      ! d670
      call compute_boundary_kernel(displ_crust_mantle,accel_crust_mantle, &
                 b_displ_crust_mantle,nspec_crust_mantle,iregion_code, &
                 ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
                 xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                 etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
                 gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,hprime_xx,hprime_yy,hprime_zz, &
                 rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
                 kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
                 c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
                 c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
                 c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
                 c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
                 c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                 c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                 k_top,ibelm_670_top,normal_670,d670_kl_top,fluid_solid_boundary,NSPEC2D_670)

      call compute_boundary_kernel(displ_crust_mantle,accel_crust_mantle, &
                 b_displ_crust_mantle,nspec_crust_mantle,iregion_code, &
                 ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
                 xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                 etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
                 gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,hprime_xx,hprime_yy,hprime_zz, &
                 rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
                 kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
                 c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
                 c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
                 c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
                 c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
                 c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                 c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                 k_bot,ibelm_670_bot,normal_670,d670_kl_bot,fluid_solid_boundary,NSPEC2D_670)

      d670_kl = d670_kl + (d670_kl_top - d670_kl_bot) * deltat

      ! CMB
      fluid_solid_boundary = .true.
      iregion_code = IREGION_CRUST_MANTLE
      call compute_boundary_kernel(displ_crust_mantle,accel_crust_mantle, &
                 b_displ_crust_mantle,nspec_crust_mantle,iregion_code, &
                 ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
                 xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                 etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
                 gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,hprime_xx,hprime_yy,hprime_zz, &
                 rhostore_crust_mantle,kappavstore_crust_mantle, muvstore_crust_mantle, &
                 kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
                 c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
                 c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
                 c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
                 c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
                 c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                 c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                 k_top,ibelm_bottom_crust_mantle,normal_top_outer_core, &
                 cmb_kl_top,fluid_solid_boundary,NSPEC2D_CMB)

      iregion_code = IREGION_OUTER_CORE
      call compute_boundary_kernel(vector_displ_outer_core,vector_accel_outer_core, &
                 b_vector_displ_outer_core,nspec_outer_core, &
                 iregion_code,ystore_outer_core,zstore_outer_core,ibool_outer_core,ispec_is_tiso_outer_core, &
                 xix_outer_core,xiy_outer_core,xiz_outer_core, &
                 etax_outer_core,etay_outer_core,etaz_outer_core,&
                 gammax_outer_core,gammay_outer_core,gammaz_outer_core,hprime_xx,hprime_yy,hprime_zz, &
                 rhostore_outer_core,kappavstore_outer_core,dummy_array, &
                 dummy_array,dummy_array,dummy_array, &
                 dummy_array,dummy_array,dummy_array,dummy_array, &
                 dummy_array,dummy_array,dummy_array, &
                 dummy_array,dummy_array,dummy_array,dummy_array, &
                 dummy_array,dummy_array,dummy_array, &
                 dummy_array,dummy_array,dummy_array,dummy_array, &
                 dummy_array,dummy_array,dummy_array, &
                 k_bot,ibelm_top_outer_core,normal_top_outer_core, &
                 cmb_kl_bot,fluid_solid_boundary,NSPEC2D_CMB)

      cmb_kl = cmb_kl + (cmb_kl_top - cmb_kl_bot) * deltat

      ! ICB
      fluid_solid_boundary = .true.
      call compute_boundary_kernel(vector_displ_outer_core,vector_accel_outer_core, &
                 b_vector_displ_outer_core,nspec_outer_core, &
                 iregion_code,ystore_outer_core,zstore_outer_core,ibool_outer_core,ispec_is_tiso_outer_core, &
                 xix_outer_core,xiy_outer_core,xiz_outer_core, &
                 etax_outer_core,etay_outer_core,etaz_outer_core,&
                 gammax_outer_core,gammay_outer_core,gammaz_outer_core,hprime_xx,hprime_yy,hprime_zz, &
                 rhostore_outer_core,kappavstore_outer_core,dummy_array, &
                 dummy_array,dummy_array,dummy_array, &
                 dummy_array,dummy_array,dummy_array,dummy_array, &
                 dummy_array,dummy_array,dummy_array, &
                 dummy_array,dummy_array,dummy_array,dummy_array, &
                 dummy_array,dummy_array,dummy_array, &
                 dummy_array,dummy_array,dummy_array,dummy_array, &
                 dummy_array,dummy_array,dummy_array, &
                 k_top,ibelm_bottom_outer_core,normal_bottom_outer_core, &
                 icb_kl_top,fluid_solid_boundary,NSPEC2D_ICB)

      iregion_code = IREGION_INNER_CORE
      call compute_boundary_kernel(displ_inner_core,accel_inner_core, &
                 b_displ_inner_core,nspec_inner_core,iregion_code, &
                 ystore_inner_core,zstore_inner_core,ibool_inner_core,ispec_is_tiso_inner_core, &
                 xix_inner_core,xiy_inner_core,xiz_inner_core, &
                 etax_inner_core,etay_inner_core,etaz_inner_core,&
                 gammax_inner_core,gammay_inner_core,gammaz_inner_core,hprime_xx,hprime_yy,hprime_zz, &
                 rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core, &
                 dummy_array,dummy_array,dummy_array, &
                 c11store_inner_core,c12store_inner_core,c13store_inner_core,dummy_array, &
                 dummy_array,dummy_array,dummy_array, &
                 dummy_array,dummy_array,dummy_array,dummy_array, &
                 c33store_inner_core,dummy_array,dummy_array, &
                 dummy_array,c44store_inner_core,dummy_array,dummy_array, &
                 dummy_array,dummy_array,dummy_array, &
                 k_bot,ibelm_top_inner_core,normal_bottom_outer_core, &
                 icb_kl_bot,fluid_solid_boundary,NSPEC2D_ICB)

      icb_kl = icb_kl + (icb_kl_top - icb_kl_bot) * deltat
    endif

    ! approximate hessian
    if( APPROXIMATE_HESS_KL ) then
      call compute_kernels_hessian(ibool_crust_mantle, &
                          hess_kl_crust_mantle,&
                          accel_crust_mantle,b_accel_crust_mantle, &
                          deltat)
    endif

  endif ! end of if computing kernels

!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------

