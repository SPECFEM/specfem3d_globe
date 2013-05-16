
!! DK DK
!! DK DK this first part handles the cases SIMULATION_TYPE == 1 and SIMULATION_TYPE == 2
!! DK DK it also handles the cases NOISE_TOMOGRAPHY == 1 and NOISE_TOMOGRAPHY == 2
!! DK DK

    ! Newmark time scheme update

    ! mantle
    do i=1,NGLOB_CRUST_MANTLE
      displ_crust_mantle(:,i) = displ_crust_mantle(:,i) &
        + deltat*veloc_crust_mantle(:,i) + deltatsqover2*accel_crust_mantle(:,i)
      veloc_crust_mantle(:,i) = veloc_crust_mantle(:,i) &
        + deltatover2*accel_crust_mantle(:,i)
      accel_crust_mantle(:,i) = 0._CUSTOM_REAL
    enddo
    ! outer core
    do i=1,NGLOB_OUTER_CORE
      displ_outer_core(i) = displ_outer_core(i) &
        + deltat*veloc_outer_core(i) + deltatsqover2*accel_outer_core(i)
      veloc_outer_core(i) = veloc_outer_core(i) &
        + deltatover2*accel_outer_core(i)
      accel_outer_core(i) = 0._CUSTOM_REAL
    enddo
    ! inner core
    do i=1,NGLOB_INNER_CORE
      displ_inner_core(:,i) = displ_inner_core(:,i) &
        + deltat*veloc_inner_core(:,i) + deltatsqover2*accel_inner_core(:,i)
      veloc_inner_core(:,i) = veloc_inner_core(:,i) &
        + deltatover2*accel_inner_core(:,i)
      accel_inner_core(:,i) = 0._CUSTOM_REAL
    enddo

!   ! backward field
!   if (SIMULATION_TYPE == 3) then
!     ! mantle
!     do i=1,NGLOB_CRUST_MANTLE
!       b_displ_crust_mantle(:,i) = b_displ_crust_mantle(:,i) &
!         + b_deltat*b_veloc_crust_mantle(:,i) + b_deltatsqover2*b_accel_crust_mantle(:,i)
!       b_veloc_crust_mantle(:,i) = b_veloc_crust_mantle(:,i) &
!         + b_deltatover2*b_accel_crust_mantle(:,i)
!       b_accel_crust_mantle(:,i) = 0._CUSTOM_REAL
!     enddo
!     ! outer core
!     do i=1,NGLOB_OUTER_CORE
!       b_displ_outer_core(i) = b_displ_outer_core(i) &
!         + b_deltat*b_veloc_outer_core(i) + b_deltatsqover2*b_accel_outer_core(i)
!       b_veloc_outer_core(i) = b_veloc_outer_core(i) &
!         + b_deltatover2*b_accel_outer_core(i)
!       b_accel_outer_core(i) = 0._CUSTOM_REAL
!     enddo
!     ! inner core
!     do i=1,NGLOB_INNER_CORE
!       b_displ_inner_core(:,i) = b_displ_inner_core(:,i) &
!         + b_deltat*b_veloc_inner_core(:,i) + b_deltatsqover2*b_accel_inner_core(:,i)
!       b_veloc_inner_core(:,i) = b_veloc_inner_core(:,i) &
!         + b_deltatover2*b_accel_inner_core(:,i)
!       b_accel_inner_core(:,i) = 0._CUSTOM_REAL
!     enddo
!   endif ! SIMULATION_TYPE == 3

    ! integral of strain for adjoint movie volume
    if(MOVIE_VOLUME .and. (MOVIE_VOLUME_TYPE == 2 .or. MOVIE_VOLUME_TYPE == 3) ) then
! do *NOT* use array syntax for that loop, otherwise you will get a compiler error when MOVIE_VOLUME is off
! because the shape of the arrays will not match (due to some arrays purposely declared with a dummy size of 1)
      do ispec = 1,NSPEC_CRUST_MANTLE
        Iepsilondev_crust_mantle(:,:,:,:,ispec) = Iepsilondev_crust_mantle(:,:,:,:,ispec)  &
                                              + deltat*epsilondev_crust_mantle(:,:,:,:,ispec)
        Ieps_trace_over_3_crust_mantle(:,:,:,ispec) = Ieps_trace_over_3_crust_mantle(:,:,:,ispec) &
                                              + deltat*eps_trace_over_3_crust_mantle(:,:,:,ispec)
      enddo
    endif

    ! compute the maximum of the norm of the displacement
    ! in all the slices using an MPI reduction
    ! and output timestamp file to check that simulation is running fine
    if(mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == it_begin+4 .or. it == it_end) then
      call check_simulation_stability(it,displ_crust_mantle,displ_inner_core,displ_outer_core, &
                          eps_trace_over_3_crust_mantle,epsilondev_crust_mantle, &
                          SIMULATION_TYPE,OUTPUT_FILES,time_start,DT,t0,NSTEP, &
                          it_begin,it_end,NUMBER_OF_THIS_RUN,NUMBER_OF_RUNS,myrank)
!     if (SIMULATION_TYPE == 3) then
!       call check_simulation_stability(it,b_displ_crust_mantle,b_displ_inner_core,b_displ_outer_core, &
!                         eps_trace_over_3_crust_mantle,epsilondev_crust_mantle, &
!                         SIMULATION_TYPE,OUTPUT_FILES,time_start,DT,t0,NSTEP, &
!                         it_begin,it_end,NUMBER_OF_THIS_RUN,NUMBER_OF_RUNS,myrank)
!     endif
    endif

    ! ****************************************************
    !   big loop over all spectral elements in the fluid
    ! ****************************************************

    ! compute internal forces in the fluid region
    if(CUSTOM_REAL == SIZE_REAL) then
      time = sngl((dble(it-1)*DT-t0)*scale_t_inv)
    else
      time = (dble(it-1)*DT-t0)*scale_t_inv
    endif

    iphase = 0 ! do not start any non blocking communications at this stage
    icall = 1  ! compute all the outer elements first in the case of non blocking MPI

    if( USE_DEVILLE_PRODUCTS_VAL ) then
      ! uses Deville et al. (2002) routine
      call compute_forces_outer_core_Dev(time,deltat,two_omega_earth, &
           A_array_rotation,B_array_rotation,d_ln_density_dr_table, &
           minus_rho_g_over_kappa_fluid, &
           displ_outer_core,accel_outer_core,div_displ_outer_core, &
           xstore_outer_core,ystore_outer_core,zstore_outer_core, &
           xix_outer_core,xiy_outer_core,xiz_outer_core, &
           etax_outer_core,etay_outer_core,etaz_outer_core, &
           gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
          is_on_a_slice_edge_outer_core, &
          myrank,iproc_xi,iproc_eta,ichunk,addressing, &
          iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
          npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
          iboolfaces_outer_core,iboolcorner_outer_core, &
          iprocfrom_faces,iprocto_faces, &
          iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
          buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
          buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar,iphase,icall, &
           hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
           wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
           ibool_outer_core,MOVIE_VOLUME)
    else
      ! div_displ_outer_core is initialized to zero in the following subroutine.
      call compute_forces_outer_core(time,deltat,two_omega_earth, &
           A_array_rotation,B_array_rotation,d_ln_density_dr_table, &
           minus_rho_g_over_kappa_fluid, &
           displ_outer_core,accel_outer_core,div_displ_outer_core, &
           xstore_outer_core,ystore_outer_core,zstore_outer_core, &
           xix_outer_core,xiy_outer_core,xiz_outer_core, &
           etax_outer_core,etay_outer_core,etaz_outer_core, &
           gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
          is_on_a_slice_edge_outer_core, &
          myrank,iproc_xi,iproc_eta,ichunk,addressing, &
          iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
          npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
          iboolfaces_outer_core,iboolcorner_outer_core, &
          iprocfrom_faces,iprocto_faces, &
          iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
          buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
          buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar,iphase,icall, &
           hprime_xx,hprime_yy,hprime_zz,hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
           wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
           ibool_outer_core,MOVIE_VOLUME)
    endif

!   if (SIMULATION_TYPE == 3) then
!     ! note on backward/reconstructed wavefields:
!     !       time for b_displ( it=1 ) corresponds to (NSTEP - 1)*DT - t0  (after Newmark scheme...)
!     !       as we start with saved wavefields b_displ( 1 ) <-> displ( NSTEP ) which correspond
!     !       to a time (NSTEP - (it-1) - 1)*DT - t0
!     !       for reconstructing the rotational contributions
!     if(CUSTOM_REAL == SIZE_REAL) then
!       time = sngl((dble(NSTEP-it)*DT-t0)*scale_t_inv)
!     else
!       time = (dble(NSTEP-it)*DT-t0)*scale_t_inv
!     endif

!     b_iphase = 0 ! do not start any non blocking communications at this stage
!     b_icall = 1  ! compute all the outer elements first in the case of non blocking MPI

!     if( USE_DEVILLE_PRODUCTS_VAL ) then
!       ! uses Deville et al. (2002) routine
!       call compute_forces_outer_core_Dev(time,b_deltat,b_two_omega_earth, &
!          b_A_array_rotation,b_B_array_rotation,d_ln_density_dr_table, &
!          minus_rho_g_over_kappa_fluid, &
!          b_displ_outer_core,b_accel_outer_core,b_div_displ_outer_core, &
!          xstore_outer_core,ystore_outer_core,zstore_outer_core, &
!          xix_outer_core,xiy_outer_core,xiz_outer_core, &
!          etax_outer_core,etay_outer_core,etaz_outer_core, &
!          gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
!         is_on_a_slice_edge_outer_core, &
!         myrank,iproc_xi,iproc_eta,ichunk,addressing, &
!         iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
!         npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
!         iboolfaces_outer_core,iboolcorner_outer_core, &
!         iprocfrom_faces,iprocto_faces, &
!         iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
!         b_buffer_send_faces,b_buffer_received_faces,npoin2D_max_all_CM_IC, &
!         b_buffer_send_chunkcorn_scalar,b_buffer_recv_chunkcorn_scalar,b_iphase,b_icall, &
!          hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
!          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
!          ibool_outer_core,MOVIE_VOLUME)
!     else
!       call compute_forces_outer_core(time,b_deltat,b_two_omega_earth, &
!          b_A_array_rotation,b_B_array_rotation,d_ln_density_dr_table, &
!          minus_rho_g_over_kappa_fluid, &
!          b_displ_outer_core,b_accel_outer_core,b_div_displ_outer_core, &
!          xstore_outer_core,ystore_outer_core,zstore_outer_core, &
!          xix_outer_core,xiy_outer_core,xiz_outer_core, &
!          etax_outer_core,etay_outer_core,etaz_outer_core, &
!          gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
!         is_on_a_slice_edge_outer_core, &
!         myrank,iproc_xi,iproc_eta,ichunk,addressing, &
!         iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
!         npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
!         iboolfaces_outer_core,iboolcorner_outer_core, &
!         iprocfrom_faces,iprocto_faces, &
!         iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
!         b_buffer_send_faces,b_buffer_received_faces,npoin2D_max_all_CM_IC, &
!         b_buffer_send_chunkcorn_scalar,b_buffer_recv_chunkcorn_scalar,b_iphase,b_icall, &
!          hprime_xx,hprime_yy,hprime_zz,hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
!          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
!          ibool_outer_core,MOVIE_VOLUME)
!     endif
!   endif

    ! Stacey absorbing boundaries
    if(NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) then
      call compute_stacey_outer_core_forward(ichunk,SAVE_FORWARD, &
                              it,ibool_outer_core, &
                              veloc_outer_core,accel_outer_core, &
                              vp_outer_core,wgllwgll_xz,wgllwgll_yz,wgllwgll_xy, &
                              jacobian2D_bottom_outer_core, &
                              jacobian2D_xmin_outer_core,jacobian2D_xmax_outer_core, &
                              jacobian2D_ymin_outer_core,jacobian2D_ymax_outer_core, &
                              ibelm_bottom_outer_core, &
                              ibelm_xmin_outer_core,ibelm_xmax_outer_core, &
                              ibelm_ymin_outer_core,ibelm_ymax_outer_core, &
                              nimin_outer_core,nimax_outer_core, &
                              njmin_outer_core,njmax_outer_core, &
                              nkmin_xi_outer_core,nkmin_eta_outer_core, &
                              NSPEC2D_BOTTOM, &
                              nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
                              nspec2D_ymin_outer_core,nspec2D_ymax_outer_core, &
                              reclen_zmin, &
                              reclen_xmin_outer_core,reclen_xmax_outer_core, &
                              reclen_ymin_outer_core,reclen_ymax_outer_core, &
                              nabs_zmin_oc, &
                              nabs_xmin_oc,nabs_xmax_oc,nabs_ymin_oc,nabs_ymax_oc, &
                              absorb_zmin_outer_core, &
                              absorb_xmin_outer_core,absorb_xmax_outer_core, &
                              absorb_ymin_outer_core,absorb_ymax_outer_core)
!     if (SIMULATION_TYPE == 3) then
!       call compute_stacey_outer_core_backward(ichunk, &
!                             NSTEP,it,ibool_outer_core, &
!                             b_accel_outer_core, &
!                             ibelm_bottom_outer_core, &
!                             ibelm_xmin_outer_core,ibelm_xmax_outer_core, &
!                             ibelm_ymin_outer_core,ibelm_ymax_outer_core, &
!                             nimin_outer_core,nimax_outer_core, &
!                             njmin_outer_core,njmax_outer_core, &
!                             nkmin_xi_outer_core,nkmin_eta_outer_core, &
!                             NSPEC2D_BOTTOM, &
!                             nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
!                             nspec2D_ymin_outer_core,nspec2D_ymax_outer_core, &
!                             reclen_zmin, &
!                             reclen_xmin_outer_core,reclen_xmax_outer_core, &
!                             reclen_ymin_outer_core,reclen_ymax_outer_core, &
!                             nabs_zmin_oc, &
!                             nabs_xmin_oc,nabs_xmax_oc,nabs_ymin_oc,nabs_ymax_oc, &
!                             absorb_zmin_outer_core, &
!                             absorb_xmin_outer_core,absorb_xmax_outer_core, &
!                             absorb_ymin_outer_core,absorb_ymax_outer_core)
!!      call compute_stacey_outer_core(ichunk,SIMULATION_TYPE,SAVE_FORWARD, &
!!                            NSTEP,it,ibool_outer_core, &
!!                            veloc_outer_core,accel_outer_core,b_accel_outer_core, &
!!                            vp_outer_core,wgllwgll_xz,wgllwgll_yz,wgllwgll_xy, &
!!                            jacobian2D_bottom_outer_core, &
!!                            jacobian2D_xmin_outer_core,jacobian2D_xmax_outer_core, &
!!                            jacobian2D_ymin_outer_core,jacobian2D_ymax_outer_core, &
!!                            ibelm_bottom_outer_core, &
!!                            ibelm_xmin_outer_core,ibelm_xmax_outer_core, &
!!                            ibelm_ymin_outer_core,ibelm_ymax_outer_core, &
!!                            nimin_outer_core,nimax_outer_core, &
!!                            njmin_outer_core,njmax_outer_core, &
!!                            nkmin_xi_outer_core,nkmin_eta_outer_core, &
!!                            NSPEC2D_BOTTOM, &
!!                            nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
!!                            nspec2D_ymin_outer_core,nspec2D_ymax_outer_core, &
!!                            reclen_zmin, &
!!                            reclen_xmin_outer_core,reclen_xmax_outer_core, &
!!                            reclen_ymin_outer_core,reclen_ymax_outer_core, &
!!                            nabs_zmin_oc, &
!!                            nabs_xmin_oc,nabs_xmax_oc,nabs_ymin_oc,nabs_ymax_oc, &
!!                            absorb_zmin_outer_core, &
!!                            absorb_xmin_outer_core,absorb_xmax_outer_core, &
!!                            absorb_ymin_outer_core,absorb_ymax_outer_core)
!     endif
    endif ! Stacey conditions


    ! ****************************************************
    ! **********  add matching with solid part  **********
    ! ****************************************************

    ! only for elements in first matching layer in the fluid

    !---
    !--- couple with mantle at the top of the outer core
    !---
    if(ACTUALLY_COUPLE_FLUID_CMB) then
      call compute_coupling_fluid_CMB(displ_crust_mantle, &
                            ibool_crust_mantle,ibelm_bottom_crust_mantle,  &
                            accel_outer_core, &
                            normal_top_outer_core,jacobian2D_top_outer_core, &
                            wgllwgll_xy,ibool_outer_core,ibelm_top_outer_core, &
                            NSPEC2D_TOP(IREGION_OUTER_CORE))
!     if (SIMULATION_TYPE == 3) then
!       call compute_coupling_fluid_CMB(b_displ_crust_mantle, &
!                           ibool_crust_mantle,ibelm_bottom_crust_mantle,  &
!                           b_accel_outer_core, &
!                           normal_top_outer_core,jacobian2D_top_outer_core, &
!                           wgllwgll_xy,ibool_outer_core,ibelm_top_outer_core, &
!                           NSPEC2D_TOP(IREGION_OUTER_CORE))
!     endif
    endif

    !---
    !--- couple with inner core at the bottom of the outer core
    !---
    if(ACTUALLY_COUPLE_FLUID_ICB) then
      call compute_coupling_fluid_ICB(displ_inner_core, &
                            ibool_inner_core,ibelm_top_inner_core,  &
                            accel_outer_core, &
                            normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
                            wgllwgll_xy,ibool_outer_core,ibelm_bottom_outer_core, &
                            NSPEC2D_BOTTOM(IREGION_OUTER_CORE))
!     if (SIMULATION_TYPE == 3) then
!       call compute_coupling_fluid_ICB(b_displ_inner_core, &
!                           ibool_inner_core,ibelm_top_inner_core,  &
!                           b_accel_outer_core, &
!                           normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
!                           wgllwgll_xy,ibool_outer_core,ibelm_bottom_outer_core, &
!                           NSPEC2D_BOTTOM(IREGION_OUTER_CORE))
!     endif
    endif

    ! assemble all the contributions between slices using MPI

    ! outer core
      iphase = 1 ! start the non blocking communications
      call assemble_MPI_scalar(myrank,accel_outer_core,NGLOB_OUTER_CORE, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
            npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
            iboolfaces_outer_core,iboolcorner_outer_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_OUTER_CORE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_OUTER_CORE),NGLOB2DMAX_YMIN_YMAX(IREGION_OUTER_CORE), &
            NGLOB2DMAX_XY_OC_VAL,NCHUNKS_VAL,iphase)

      icall = 2 ! now compute all the inner elements in the case of non blocking MPI

      if( USE_DEVILLE_PRODUCTS_VAL ) then
        ! uses Deville et al. (2002) routine
        call compute_forces_outer_core_Dev(time,deltat,two_omega_earth, &
           A_array_rotation,B_array_rotation,d_ln_density_dr_table, &
           minus_rho_g_over_kappa_fluid,displ_outer_core,accel_outer_core,div_displ_outer_core, &
           xstore_outer_core,ystore_outer_core,zstore_outer_core, &
           xix_outer_core,xiy_outer_core,xiz_outer_core, &
           etax_outer_core,etay_outer_core,etaz_outer_core, &
           gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
          is_on_a_slice_edge_outer_core, &
          myrank,iproc_xi,iproc_eta,ichunk,addressing, &
          iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
          npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
          iboolfaces_outer_core,iboolcorner_outer_core, &
          iprocfrom_faces,iprocto_faces, &
          iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
          buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
          buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar,iphase,icall, &
           hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
           wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
           ibool_outer_core,MOVIE_VOLUME)
      else
        ! div_displ_outer_core is initialized to zero in the following subroutine.
        call compute_forces_outer_core(time,deltat,two_omega_earth, &
           A_array_rotation,B_array_rotation,d_ln_density_dr_table, &
           minus_rho_g_over_kappa_fluid,displ_outer_core,accel_outer_core,div_displ_outer_core, &
           xstore_outer_core,ystore_outer_core,zstore_outer_core, &
           xix_outer_core,xiy_outer_core,xiz_outer_core, &
           etax_outer_core,etay_outer_core,etaz_outer_core, &
           gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
          is_on_a_slice_edge_outer_core, &
          myrank,iproc_xi,iproc_eta,ichunk,addressing, &
          iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
          npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
          iboolfaces_outer_core,iboolcorner_outer_core, &
          iprocfrom_faces,iprocto_faces, &
          iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
          buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
          buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar,iphase,icall, &
           hprime_xx,hprime_yy,hprime_zz,hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
           wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
           ibool_outer_core,MOVIE_VOLUME)
      endif

      do while (iphase <= 7) ! make sure the last communications are finished and processed
        call assemble_MPI_scalar(myrank,accel_outer_core,NGLOB_OUTER_CORE, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
            npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
            iboolfaces_outer_core,iboolcorner_outer_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_OUTER_CORE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_OUTER_CORE),NGLOB2DMAX_YMIN_YMAX(IREGION_OUTER_CORE), &
            NGLOB2DMAX_XY_OC_VAL,NCHUNKS_VAL,iphase)
      enddo

    ! multiply by the inverse of the mass matrix and update velocity
    do i=1,NGLOB_OUTER_CORE
      accel_outer_core(i) = accel_outer_core(i)*rmass_outer_core(i)
      veloc_outer_core(i) = veloc_outer_core(i) + deltatover2*accel_outer_core(i)
    enddo

! ------------------- new non blocking implementation -------------------

!   if (SIMULATION_TYPE == 3) then

!   ! outer core
!       b_iphase = 1 ! start the non blocking communications
!       call assemble_MPI_scalar(myrank,b_accel_outer_core,NGLOB_OUTER_CORE, &
!           iproc_xi,iproc_eta,ichunk,addressing, &
!           iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
!           npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
!           iboolfaces_outer_core,iboolcorner_outer_core, &
!           iprocfrom_faces,iprocto_faces, &
!           iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
!           b_buffer_send_faces,b_buffer_received_faces,npoin2D_max_all_CM_IC, &
!           b_buffer_send_chunkcorn_scalar,b_buffer_recv_chunkcorn_scalar, &
!           NUMMSGS_FACES,NCORNERSCHUNKS, &
!           NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_OUTER_CORE), &
!           NGLOB2DMAX_XMIN_XMAX(IREGION_OUTER_CORE),NGLOB2DMAX_YMIN_YMAX(IREGION_OUTER_CORE), &
!           NGLOB2DMAX_XY_OC_VAL,NCHUNKS_VAL,b_iphase)

!       b_icall = 2 ! now compute all the inner elements in the case of non blocking MPI

!       if( USE_DEVILLE_PRODUCTS_VAL ) then
!         ! uses Deville et al. (2002) routine
!         call compute_forces_outer_core_Dev(time,b_deltat,b_two_omega_earth, &
!          b_A_array_rotation,b_B_array_rotation,d_ln_density_dr_table, &
!          minus_rho_g_over_kappa_fluid, &
!          b_displ_outer_core,b_accel_outer_core,b_div_displ_outer_core, &
!          xstore_outer_core,ystore_outer_core,zstore_outer_core, &
!          xix_outer_core,xiy_outer_core,xiz_outer_core, &
!          etax_outer_core,etay_outer_core,etaz_outer_core, &
!          gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
!         is_on_a_slice_edge_outer_core, &
!         myrank,iproc_xi,iproc_eta,ichunk,addressing, &
!         iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
!         npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
!         iboolfaces_outer_core,iboolcorner_outer_core, &
!         iprocfrom_faces,iprocto_faces, &
!         iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
!         b_buffer_send_faces,b_buffer_received_faces,npoin2D_max_all_CM_IC, &
!         b_buffer_send_chunkcorn_scalar,b_buffer_recv_chunkcorn_scalar,b_iphase,b_icall, &
!          hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
!          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
!          ibool_outer_core,MOVIE_VOLUME)
!       else
!         ! div_displ_outer_core is initialized to zero in the following subroutine.
!         call compute_forces_outer_core(time,b_deltat,b_two_omega_earth, &
!          b_A_array_rotation,b_B_array_rotation,d_ln_density_dr_table, &
!          minus_rho_g_over_kappa_fluid, &
!          b_displ_outer_core,b_accel_outer_core,b_div_displ_outer_core, &
!          xstore_outer_core,ystore_outer_core,zstore_outer_core, &
!          xix_outer_core,xiy_outer_core,xiz_outer_core, &
!          etax_outer_core,etay_outer_core,etaz_outer_core, &
!          gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
!         is_on_a_slice_edge_outer_core, &
!         myrank,iproc_xi,iproc_eta,ichunk,addressing, &
!         iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
!         npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
!         iboolfaces_outer_core,iboolcorner_outer_core, &
!         iprocfrom_faces,iprocto_faces, &
!         iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
!         b_buffer_send_faces,b_buffer_received_faces,npoin2D_max_all_CM_IC, &
!         b_buffer_send_chunkcorn_scalar,b_buffer_recv_chunkcorn_scalar,b_iphase,b_icall, &
!          hprime_xx,hprime_yy,hprime_zz,hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
!          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
!          ibool_outer_core,MOVIE_VOLUME)
!       endif

!       do while (b_iphase <= 7) ! make sure the last communications are finished and processed
!         call assemble_MPI_scalar(myrank,b_accel_outer_core,NGLOB_OUTER_CORE, &
!           iproc_xi,iproc_eta,ichunk,addressing, &
!           iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
!           npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
!           iboolfaces_outer_core,iboolcorner_outer_core, &
!           iprocfrom_faces,iprocto_faces, &
!           iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
!           b_buffer_send_faces,b_buffer_received_faces,npoin2D_max_all_CM_IC, &
!           b_buffer_send_chunkcorn_scalar,b_buffer_recv_chunkcorn_scalar, &
!           NUMMSGS_FACES,NCORNERSCHUNKS, &
!           NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_OUTER_CORE), &
!           NGLOB2DMAX_XMIN_XMAX(IREGION_OUTER_CORE),NGLOB2DMAX_YMIN_YMAX(IREGION_OUTER_CORE), &
!           NGLOB2DMAX_XY_OC_VAL,NCHUNKS_VAL,b_iphase)
!       enddo

! ------------------- new non blocking implementation -------------------

!     ! Newmark time scheme - corrector for fluid parts
!     do i=1,NGLOB_OUTER_CORE
!       b_accel_outer_core(i) = b_accel_outer_core(i)*rmass_outer_core(i)
!       b_veloc_outer_core(i) = b_veloc_outer_core(i) + b_deltatover2*b_accel_outer_core(i)
!     enddo

!   endif

    ! ****************************************************
    !   big loop over all spectral elements in the solid
    ! ****************************************************

    ! compute internal forces in the solid regions

    ! for anisotropy and gravity, x y and z contain r theta and phi

    iphase = 0 ! do not start any non blocking communications at this stage
    iphase_CC = 0 ! do not start any non blocking communications at this stage
    icall = 1  ! compute all the outer elements first in the case of non blocking MPI

    if( USE_DEVILLE_PRODUCTS_VAL ) then
      call compute_forces_crust_mantle_Dev(minus_gravity_table,density_table,minus_deriv_gravity_table, &
          displ_crust_mantle,accel_crust_mantle, &
          xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
          xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
          etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
          gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
            is_on_a_slice_edge_crust_mantle,icall, &
            accel_inner_core,ibool_inner_core,idoubling_inner_core, &
            myrank,iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector,iphase, &
            nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
            npoin2D_cube_from_slices,buffer_all_cube_from_slices,buffer_slices,ibool_central_cube, &
            receiver_cube_from_slices,ibelm_bottom_inner_core,NSPEC2D_BOTTOM_IC,INCLUDE_CENTRAL_CUBE,iphase_CC, &
          hprime_xx,hprime_xxT, &
          hprimewgll_xx,hprimewgll_xxT, &
          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
          kappavstore_crust_mantle,kappahstore_crust_mantle,muvstore_crust_mantle, &
          muhstore_crust_mantle,eta_anisostore_crust_mantle, &
          c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
          c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
          c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
          c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
          c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
          c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
          c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
          ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
          R_memory_crust_mantle,epsilondev_crust_mantle, &
          eps_trace_over_3_crust_mantle,one_minus_sum_beta_crust_mantle, &
          alphaval,betaval,gammaval,factor_common_crust_mantle, &
          size(factor_common_crust_mantle,2), size(factor_common_crust_mantle,3), &
          size(factor_common_crust_mantle,4), size(factor_common_crust_mantle,5) )
    else
      call compute_forces_crust_mantle(minus_gravity_table,density_table,minus_deriv_gravity_table, &
          displ_crust_mantle,accel_crust_mantle, &
          xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
          xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
          etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
          gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
            is_on_a_slice_edge_crust_mantle,icall, &
            accel_inner_core,ibool_inner_core,idoubling_inner_core, &
            myrank,iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector,iphase, &
            nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
            npoin2D_cube_from_slices,buffer_all_cube_from_slices,buffer_slices,ibool_central_cube, &
            receiver_cube_from_slices,ibelm_bottom_inner_core,NSPEC2D_BOTTOM_IC,INCLUDE_CENTRAL_CUBE,iphase_CC, &
          hprime_xx,hprime_yy,hprime_zz, &
          hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
          kappavstore_crust_mantle,kappahstore_crust_mantle,muvstore_crust_mantle, &
          muhstore_crust_mantle,eta_anisostore_crust_mantle, &
          c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
          c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
          c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
          c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
          c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
          c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
          c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
          ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
          R_memory_crust_mantle,epsilondev_crust_mantle, &
          eps_trace_over_3_crust_mantle,one_minus_sum_beta_crust_mantle, &
          alphaval,betaval,gammaval,factor_common_crust_mantle, &
          size(factor_common_crust_mantle,2), size(factor_common_crust_mantle,3), &
          size(factor_common_crust_mantle,4), size(factor_common_crust_mantle,5) )
    endif

!   if (SIMULATION_TYPE == 3 ) then

!     b_iphase = 0 ! do not start any non blocking communications at this stage
!     b_iphase_CC = 0 ! do not start any non blocking communications at this stage
!     b_icall = 1  ! compute all the outer elements first in the case of non blocking MPI

!   ! for anisotropy and gravity, x y and z contain r theta and phi
!     if( USE_DEVILLE_PRODUCTS_VAL ) then
!       call compute_forces_crust_mantle_Dev(minus_gravity_table,density_table,minus_deriv_gravity_table, &
!         b_displ_crust_mantle,b_accel_crust_mantle, &
!         xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
!         xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
!         etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
!         gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
!           is_on_a_slice_edge_crust_mantle,b_icall, &
!           b_accel_inner_core,ibool_inner_core,idoubling_inner_core, &
!           myrank,iproc_xi,iproc_eta,ichunk,addressing, &
!           iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
!           npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
!           iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
!           iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
!           npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
!           iboolfaces_inner_core,iboolcorner_inner_core, &
!           iprocfrom_faces,iprocto_faces, &
!           iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
!           b_buffer_send_faces,b_buffer_received_faces,npoin2D_max_all_CM_IC, &
!           b_buffer_send_chunkcorn_vector,b_buffer_recv_chunkcorn_vector,b_iphase, &
!           nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
!           npoin2D_cube_from_slices,b_buffer_all_cube_from_slices,b_buffer_slices,ibool_central_cube, &
!           receiver_cube_from_slices,ibelm_bottom_inner_core,NSPEC2D_BOTTOM_IC,INCLUDE_CENTRAL_CUBE,b_iphase_CC, &
!         hprime_xx,hprime_xxT, &
!         hprimewgll_xx,hprimewgll_xxT, &
!         wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
!         kappavstore_crust_mantle,kappahstore_crust_mantle,muvstore_crust_mantle, &
!         muhstore_crust_mantle,eta_anisostore_crust_mantle, &
!         c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
!         c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
!         c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
!         c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
!         c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
!         c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
!         c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
!         ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
!         b_R_memory_crust_mantle,b_epsilondev_crust_mantle, &
!         b_eps_trace_over_3_crust_mantle,one_minus_sum_beta_crust_mantle, &
!         b_alphaval,b_betaval,b_gammaval,factor_common_crust_mantle, &
!         size(factor_common_crust_mantle,2), size(factor_common_crust_mantle,3), &
!         size(factor_common_crust_mantle,4), size(factor_common_crust_mantle,5) )
!     else
!       call compute_forces_crust_mantle(minus_gravity_table,density_table,minus_deriv_gravity_table, &
!         b_displ_crust_mantle,b_accel_crust_mantle, &
!         xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
!         xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
!         etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
!         gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
!           is_on_a_slice_edge_crust_mantle,b_icall, &
!           b_accel_inner_core,ibool_inner_core,idoubling_inner_core, &
!           myrank,iproc_xi,iproc_eta,ichunk,addressing, &
!           iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
!           npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
!           iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
!           iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
!           npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
!           iboolfaces_inner_core,iboolcorner_inner_core, &
!           iprocfrom_faces,iprocto_faces, &
!           iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
!           b_buffer_send_faces,b_buffer_received_faces,npoin2D_max_all_CM_IC, &
!           b_buffer_send_chunkcorn_vector,b_buffer_recv_chunkcorn_vector,b_iphase, &
!           nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
!           npoin2D_cube_from_slices,b_buffer_all_cube_from_slices,b_buffer_slices,ibool_central_cube, &
!           receiver_cube_from_slices,ibelm_bottom_inner_core,NSPEC2D_BOTTOM_IC,INCLUDE_CENTRAL_CUBE,b_iphase_CC, &
!         hprime_xx,hprime_yy,hprime_zz, &
!         hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
!         wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
!         kappavstore_crust_mantle,kappahstore_crust_mantle,muvstore_crust_mantle, &
!         muhstore_crust_mantle,eta_anisostore_crust_mantle, &
!         c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
!         c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
!         c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
!         c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
!         c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
!         c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
!         c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
!         ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
!         b_R_memory_crust_mantle,b_epsilondev_crust_mantle, &
!         b_eps_trace_over_3_crust_mantle,one_minus_sum_beta_crust_mantle, &
!         b_alphaval,b_betaval,b_gammaval,factor_common_crust_mantle, &
!         size(factor_common_crust_mantle,2), size(factor_common_crust_mantle,3), &
!         size(factor_common_crust_mantle,4), size(factor_common_crust_mantle,5) )
!     endif
!   endif

    ! Deville routine
    if( USE_DEVILLE_PRODUCTS_VAL ) then
      call compute_forces_inner_core_Dev(minus_gravity_table,density_table,minus_deriv_gravity_table, &
          displ_inner_core,accel_inner_core, &
          xstore_inner_core,ystore_inner_core,zstore_inner_core, &
          xix_inner_core,xiy_inner_core,xiz_inner_core, &
          etax_inner_core,etay_inner_core,etaz_inner_core, &
          gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
            is_on_a_slice_edge_inner_core,icall, &
            accel_crust_mantle,ibool_inner_core,idoubling_inner_core, &
            myrank,iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector,iphase, &
            nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
            npoin2D_cube_from_slices,buffer_all_cube_from_slices,buffer_slices,ibool_central_cube, &
            receiver_cube_from_slices,ibelm_bottom_inner_core,NSPEC2D_BOTTOM_IC,INCLUDE_CENTRAL_CUBE,iphase_CC, &
          hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
          kappavstore_inner_core,muvstore_inner_core,ibool_inner_core,idoubling_inner_core, &
          c11store_inner_core,c33store_inner_core,c12store_inner_core, &
          c13store_inner_core,c44store_inner_core, &
          R_memory_inner_core,epsilondev_inner_core, eps_trace_over_3_inner_core,&
          one_minus_sum_beta_inner_core, &
          alphaval,betaval,gammaval, &
          factor_common_inner_core, &
          size(factor_common_inner_core,2), size(factor_common_inner_core,3), &
          size(factor_common_inner_core,4), size(factor_common_inner_core,5) )
    else
      call compute_forces_inner_core(minus_gravity_table,density_table,minus_deriv_gravity_table, &
          displ_inner_core,accel_inner_core, &
          xstore_inner_core,ystore_inner_core,zstore_inner_core, &
          xix_inner_core,xiy_inner_core,xiz_inner_core, &
          etax_inner_core,etay_inner_core,etaz_inner_core, &
          gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
            is_on_a_slice_edge_inner_core,icall, &
            accel_crust_mantle,ibool_inner_core,idoubling_inner_core, &
            myrank,iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector,iphase, &
            nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
            npoin2D_cube_from_slices,buffer_all_cube_from_slices,buffer_slices,ibool_central_cube, &
            receiver_cube_from_slices,ibelm_bottom_inner_core,NSPEC2D_BOTTOM_IC,INCLUDE_CENTRAL_CUBE,iphase_CC, &
          hprime_xx,hprime_yy,hprime_zz,hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
          kappavstore_inner_core,muvstore_inner_core,ibool_inner_core,idoubling_inner_core, &
          c11store_inner_core,c33store_inner_core,c12store_inner_core, &
          c13store_inner_core,c44store_inner_core, &
          R_memory_inner_core,epsilondev_inner_core, eps_trace_over_3_inner_core,&
          one_minus_sum_beta_inner_core, &
          alphaval,betaval,gammaval, &
          factor_common_inner_core, &
          size(factor_common_inner_core,2), size(factor_common_inner_core,3), &
          size(factor_common_inner_core,4), size(factor_common_inner_core,5) )
    endif

!   if (SIMULATION_TYPE == 3) then
!     if( USE_DEVILLE_PRODUCTS_VAL ) then
!       call compute_forces_inner_core_Dev(minus_gravity_table,density_table,minus_deriv_gravity_table, &
!         b_displ_inner_core,b_accel_inner_core, &
!         xstore_inner_core,ystore_inner_core,zstore_inner_core, &
!         xix_inner_core,xiy_inner_core,xiz_inner_core, &
!         etax_inner_core,etay_inner_core,etaz_inner_core, &
!         gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
!           is_on_a_slice_edge_inner_core,b_icall, &
!           b_accel_crust_mantle,ibool_inner_core,idoubling_inner_core, &
!           myrank,iproc_xi,iproc_eta,ichunk,addressing, &
!           iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
!           npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
!           iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
!           iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
!           npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
!           iboolfaces_inner_core,iboolcorner_inner_core, &
!           iprocfrom_faces,iprocto_faces, &
!           iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
!           b_buffer_send_faces,b_buffer_received_faces,npoin2D_max_all_CM_IC, &
!           b_buffer_send_chunkcorn_vector,b_buffer_recv_chunkcorn_vector,b_iphase, &
!           nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
!           npoin2D_cube_from_slices,b_buffer_all_cube_from_slices,b_buffer_slices,ibool_central_cube, &
!           receiver_cube_from_slices,ibelm_bottom_inner_core,NSPEC2D_BOTTOM_IC,INCLUDE_CENTRAL_CUBE,b_iphase_CC, &
!         hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
!         wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
!         kappavstore_inner_core,muvstore_inner_core,ibool_inner_core,idoubling_inner_core, &
!         c11store_inner_core,c33store_inner_core,c12store_inner_core, &
!         c13store_inner_core,c44store_inner_core, &
!         b_R_memory_inner_core,b_epsilondev_inner_core, b_eps_trace_over_3_inner_core,&
!         one_minus_sum_beta_inner_core, &
!         b_alphaval,b_betaval,b_gammaval, &
!         factor_common_inner_core, &
!         size(factor_common_inner_core,2), size(factor_common_inner_core,3), &
!         size(factor_common_inner_core,4), size(factor_common_inner_core,5) )
!     else
!       call compute_forces_inner_core(minus_gravity_table,density_table,minus_deriv_gravity_table, &
!         b_displ_inner_core,b_accel_inner_core, &
!         xstore_inner_core,ystore_inner_core,zstore_inner_core, &
!         xix_inner_core,xiy_inner_core,xiz_inner_core, &
!         etax_inner_core,etay_inner_core,etaz_inner_core, &
!         gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
!           is_on_a_slice_edge_inner_core,b_icall, &
!           b_accel_crust_mantle,ibool_inner_core,idoubling_inner_core, &
!           myrank,iproc_xi,iproc_eta,ichunk,addressing, &
!           iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
!           npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
!           iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
!           iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
!           npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
!           iboolfaces_inner_core,iboolcorner_inner_core, &
!           iprocfrom_faces,iprocto_faces, &
!           iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
!           b_buffer_send_faces,b_buffer_received_faces,npoin2D_max_all_CM_IC, &
!           b_buffer_send_chunkcorn_vector,b_buffer_recv_chunkcorn_vector,b_iphase, &
!           nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
!           npoin2D_cube_from_slices,b_buffer_all_cube_from_slices,b_buffer_slices,ibool_central_cube, &
!           receiver_cube_from_slices,ibelm_bottom_inner_core,NSPEC2D_BOTTOM_IC,INCLUDE_CENTRAL_CUBE,b_iphase_CC, &
!         hprime_xx,hprime_yy,hprime_zz,hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
!         wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
!         kappavstore_inner_core,muvstore_inner_core,ibool_inner_core,idoubling_inner_core, &
!         c11store_inner_core,c33store_inner_core,c12store_inner_core, &
!         c13store_inner_core,c44store_inner_core, &
!         b_R_memory_inner_core,b_epsilondev_inner_core, b_eps_trace_over_3_inner_core,&
!         one_minus_sum_beta_inner_core, &
!         b_alphaval,b_betaval,b_gammaval, &
!         factor_common_inner_core, &
!         size(factor_common_inner_core,2), size(factor_common_inner_core,3), &
!         size(factor_common_inner_core,4), size(factor_common_inner_core,5) )
!     endif
!   endif

    ! Stacey
    if(NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) then

      call compute_stacey_crust_mantle_forward(ichunk, &
                              it,SAVE_FORWARD,ibool_crust_mantle, &
                              veloc_crust_mantle,accel_crust_mantle, &
                              jacobian2D_xmin_crust_mantle,jacobian2D_xmax_crust_mantle, &
                              jacobian2D_ymin_crust_mantle,jacobian2D_ymax_crust_mantle, &
                              wgllwgll_xz,wgllwgll_yz, &
                              normal_xmin_crust_mantle,normal_xmax_crust_mantle, &
                              normal_ymin_crust_mantle,normal_ymax_crust_mantle, &
                              rho_vp_crust_mantle,rho_vs_crust_mantle, &
                              ibelm_xmin_crust_mantle,ibelm_xmax_crust_mantle, &
                              ibelm_ymin_crust_mantle,ibelm_ymax_crust_mantle, &
                              nimin_crust_mantle,nimax_crust_mantle, &
                              njmin_crust_mantle,njmax_crust_mantle, &
                              nkmin_xi_crust_mantle,nkmin_eta_crust_mantle, &
                              nspec2D_xmin_crust_mantle,nspec2D_xmax_crust_mantle, &
                              nspec2D_ymin_crust_mantle,nspec2D_ymax_crust_mantle, &
                              reclen_xmin_crust_mantle,reclen_xmax_crust_mantle, &
                              reclen_ymin_crust_mantle,reclen_ymax_crust_mantle, &
                              nabs_xmin_cm,nabs_xmax_cm,nabs_ymin_cm,nabs_ymax_cm, &
                              absorb_xmin_crust_mantle5,absorb_xmax_crust_mantle5, &
                              absorb_ymin_crust_mantle5,absorb_ymax_crust_mantle5)

!     if(SIMULATION_TYPE == 3) then
!       call compute_stacey_crust_mantle_backward(ichunk, &
!                             NSTEP,it,ibool_crust_mantle, &
!                             b_accel_crust_mantle, &
!                             ibelm_xmin_crust_mantle,ibelm_xmax_crust_mantle, &
!                             ibelm_ymin_crust_mantle,ibelm_ymax_crust_mantle, &
!                             nimin_crust_mantle,nimax_crust_mantle, &
!                             njmin_crust_mantle,njmax_crust_mantle, &
!                             nkmin_xi_crust_mantle,nkmin_eta_crust_mantle, &
!                             nspec2D_xmin_crust_mantle,nspec2D_xmax_crust_mantle, &
!                             nspec2D_ymin_crust_mantle,nspec2D_ymax_crust_mantle, &
!                             reclen_xmin_crust_mantle,reclen_xmax_crust_mantle, &
!                             reclen_ymin_crust_mantle,reclen_ymax_crust_mantle, &
!                             nabs_xmin_cm,nabs_xmax_cm,nabs_ymin_cm,nabs_ymax_cm, &
!                             absorb_xmin_crust_mantle5,absorb_xmax_crust_mantle5, &
!                             absorb_ymin_crust_mantle5,absorb_ymax_crust_mantle5)
!!      call compute_stacey_crust_mantle(ichunk, &
!!                            NSTEP,it,SAVE_FORWARD,ibool_crust_mantle, &
!!                            veloc_crust_mantle,b_accel_crust_mantle, &
!!                            jacobian2D_xmin_crust_mantle,jacobian2D_xmax_crust_mantle, &
!!                            jacobian2D_ymin_crust_mantle,jacobian2D_ymax_crust_mantle, &
!!                            wgllwgll_xz,wgllwgll_yz, &
!!                            normal_xmin_crust_mantle,normal_xmax_crust_mantle, &
!!                            normal_ymin_crust_mantle,normal_ymax_crust_mantle, &
!!                            rho_vp_crust_mantle,rho_vs_crust_mantle, &
!!                            ibelm_xmin_crust_mantle,ibelm_xmax_crust_mantle, &
!!                            ibelm_ymin_crust_mantle,ibelm_ymax_crust_mantle, &
!!                            nimin_crust_mantle,nimax_crust_mantle, &
!!                            njmin_crust_mantle,njmax_crust_mantle, &
!!                            nkmin_xi_crust_mantle,nkmin_eta_crust_mantle, &
!!                            nspec2D_xmin_crust_mantle,nspec2D_xmax_crust_mantle, &
!!                            nspec2D_ymin_crust_mantle,nspec2D_ymax_crust_mantle, &
!!                            reclen_xmin_crust_mantle,reclen_xmax_crust_mantle, &
!!                            reclen_ymin_crust_mantle,reclen_ymax_crust_mantle, &
!!                            nabs_xmin_cm,nabs_xmax_cm,nabs_ymin_cm,nabs_ymax_cm, &
!!                            absorb_xmin_crust_mantle5,absorb_xmax_crust_mantle5, &
!!                            absorb_ymin_crust_mantle5,absorb_ymax_crust_mantle5)
!     endif

    endif ! Stacey conditions

    ! add the sources
    if (SIMULATION_TYPE == 1) &
      call compute_add_sources(myrank,NSOURCES, &
                                accel_crust_mantle,sourcearrays, &
                                DT,t0,tshift_cmt,hdur_gaussian,ibool_crust_mantle, &
                                islice_selected_source,ispec_selected_source,it, &
                                hdur,xi_source,eta_source,gamma_source,nu_source)

    ! add adjoint sources only if adjoint simulation is performed for source inversion only
    if (SIMULATION_TYPE == 2) then
      if( nadj_rec_local > 0 ) &
        call compute_add_sources_adjoint(myrank,nrec, &
                                nadj_rec_local,NSTEP,NTSTEP_BETWEEN_READ_ADJSRC, &
                                accel_crust_mantle,adj_sourcearrays, &
                                nu,xi_receiver,eta_receiver,gamma_receiver, &
                                xigll,yigll,zigll,ibool_crust_mantle, &
                                islice_selected_rec,ispec_selected_rec, &
                                NSTEP_SUB_ADJ,iadjsrc_len,iadjsrc,iadj_vec, &
                                it,it_begin,station_name,network_name,DT)
    endif

!   ! add adjoint sources and add sources for backward/reconstructed wavefield
!   if (SIMULATION_TYPE == 3) then
!     if( nadj_rec_local > 0 ) &
!       call compute_add_sources_adjoint(myrank,nrec, &
!                               nadj_rec_local,NSTEP,NTSTEP_BETWEEN_READ_ADJSRC, &
!                               accel_crust_mantle,adj_sourcearrays, &
!                               nu,xi_receiver,eta_receiver,gamma_receiver, &
!                               xigll,yigll,zigll,ibool_crust_mantle, &
!                               islice_selected_rec,ispec_selected_rec, &
!                               NSTEP_SUB_ADJ,iadjsrc_len,iadjsrc,iadj_vec, &
!                               it,it_begin,station_name,network_name,DT)
!     call compute_add_sources_backward(myrank,NSOURCES,NSTEP, &
!                               b_accel_crust_mantle,sourcearrays, &
!                               DT,t0,tshift_cmt,hdur_gaussian,ibool_crust_mantle, &
!                               islice_selected_source,ispec_selected_source,it, &
!                               hdur,xi_source,eta_source,gamma_source,nu_source)
!   endif

    ! NOISE_TOMOGRAPHY
    if ( NOISE_TOMOGRAPHY == 1 ) then
       ! the first step of noise tomography is to use |S(\omega)|^2 as a point force source at one of the receivers.
       ! hence, instead of a moment tensor 'sourcearrays', a 'noise_sourcearray' for a point force is needed.
       ! furthermore, the CMTSOLUTION needs to be zero, i.e., no earthquakes.
       ! now this must be manually set in DATA/CMTSOLUTION, by USERS.
       call add_source_master_rec_noise(myrank,nrec, &
                                NSTEP,accel_crust_mantle,noise_sourcearray, &
                                ibool_crust_mantle,islice_selected_rec,ispec_selected_rec, &
                                it,irec_master_noise)
    else if ( NOISE_TOMOGRAPHY == 2 ) then
       ! second step of noise tomography, i.e., read the surface movie saved at every timestep
       ! use the movie to drive the ensemble forward wavefield
       call noise_read_add_surface_movie(nmovie_points,accel_crust_mantle, &
                              normal_x_noise,normal_y_noise,normal_z_noise,mask_noise, &
                              ibelm_top_crust_mantle,ibool_crust_mantle, &
                              NSPEC2D_TOP(IREGION_CRUST_MANTLE),noise_surface_movie, &
                              NSTEP-it+1,jacobian2D_top_crust_mantle,wgllwgll_xy)
        ! be careful, since ensemble forward sources are reversals of generating wavefield "eta"
        ! hence the "NSTEP-it+1", i.e., start reading from the last timestep
        ! note the ensemble forward sources are generally distributed on the surface of the earth
        ! that's to say, the ensemble forward source is kind of a surface force density, not a body force density
        ! therefore, we must add it here, before applying the inverse of mass matrix
!   else if ( NOISE_TOMOGRAPHY == 3 ) then
!       ! third step of noise tomography, i.e., read the surface movie saved at every timestep
!       ! use the movie to reconstruct the ensemble forward wavefield
!       ! the ensemble adjoint wavefield is done as usual
!       ! note instead of "NSTEP-it+1", now we us "it", since reconstruction is a reversal of reversal
!       call noise_read_add_surface_movie(nmovie_points,b_accel_crust_mantle, &
!                             normal_x_noise,normal_y_noise,normal_z_noise,mask_noise, &
!                             ibelm_top_crust_mantle,ibool_crust_mantle, &
!                             NSPEC2D_TOP(IREGION_CRUST_MANTLE),noise_surface_movie, &
!                             it,jacobian2D_top_crust_mantle,wgllwgll_xy)
    endif

    ! ****************************************************
    ! **********  add matching with fluid part  **********
    ! ****************************************************

    ! only for elements in first matching layer in the solid

    !---
    !--- couple with outer core at the bottom of the mantle
    !---
    if(ACTUALLY_COUPLE_FLUID_CMB) then
      call compute_coupling_CMB_fluid(displ_crust_mantle, &
                            accel_crust_mantle, &
                            ibool_crust_mantle,ibelm_bottom_crust_mantle,  &
                            accel_outer_core, &
                            normal_top_outer_core,jacobian2D_top_outer_core, &
                            wgllwgll_xy,ibool_outer_core,ibelm_top_outer_core, &
                            RHO_TOP_OC,minus_g_cmb, &
                            NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE))
!     if(SIMULATION_TYPE == 3) then
!       call compute_coupling_CMB_fluid(b_displ_crust_mantle, &
!                           b_accel_crust_mantle, &
!                           ibool_crust_mantle,ibelm_bottom_crust_mantle,  &
!                           b_accel_outer_core, &
!                           normal_top_outer_core,jacobian2D_top_outer_core, &
!                           wgllwgll_xy,ibool_outer_core,ibelm_top_outer_core, &
!                           RHO_TOP_OC,minus_g_cmb, &
!                           NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE))
!     endif
    endif

    !---
    !--- couple with outer core at the top of the inner core
    !---
    if(ACTUALLY_COUPLE_FLUID_ICB) then
      call compute_coupling_ICB_fluid(displ_inner_core, &
                            accel_inner_core, &
                            ibool_inner_core,ibelm_top_inner_core,  &
                            accel_outer_core, &
                            normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
                            wgllwgll_xy,ibool_outer_core,ibelm_bottom_outer_core, &
                            RHO_BOTTOM_OC,minus_g_icb, &
                            NSPEC2D_TOP(IREGION_INNER_CORE))
!     if(SIMULATION_TYPE == 3) then
!       call compute_coupling_ICB_fluid(b_displ_inner_core, &
!                           b_accel_inner_core, &
!                           ibool_inner_core,ibelm_top_inner_core,  &
!                           b_accel_outer_core, &
!                           normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
!                           wgllwgll_xy,ibool_outer_core,ibelm_bottom_outer_core, &
!                           RHO_BOTTOM_OC,minus_g_icb, &
!                           NSPEC2D_TOP(IREGION_INNER_CORE))
!     endif
    endif

    ! assemble all the contributions between slices using MPI

! assemble all the contributions between slices using MPI
! crust/mantle and inner core handled in the same call
! in order to reduce the number of MPI messages by 2

      iphase = 1 ! initialize the non blocking communication counter
      iphase_CC = 1 ! initialize the non blocking communication counter for the central cube

! start the non blocking communications
      call assemble_MPI_vector(myrank,accel_crust_mantle,accel_inner_core, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector, &
            NUMMSGS_FACES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_CRUST_MANTLE), &
            NGLOB1D_RADIAL(IREGION_INNER_CORE),NCHUNKS_VAL,iphase)

      icall = 2 ! now compute all the inner elements in the case of non blocking MPI

      ! compute internal forces in the solid regions

      ! for anisotropy and gravity, x y and z contain r theta and phi

      if( USE_DEVILLE_PRODUCTS_VAL ) then
        call compute_forces_crust_mantle_Dev(minus_gravity_table,density_table,minus_deriv_gravity_table, &
          displ_crust_mantle,accel_crust_mantle, &
          xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
          xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
          etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
          gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
            is_on_a_slice_edge_crust_mantle,icall, &
            accel_inner_core,ibool_inner_core,idoubling_inner_core, &
            myrank,iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector,iphase, &
            nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
            npoin2D_cube_from_slices,buffer_all_cube_from_slices,buffer_slices,ibool_central_cube, &
            receiver_cube_from_slices,ibelm_bottom_inner_core,NSPEC2D_BOTTOM_IC,INCLUDE_CENTRAL_CUBE,iphase_CC, &
          hprime_xx,hprime_xxT, &
          hprimewgll_xx,hprimewgll_xxT, &
          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
          kappavstore_crust_mantle,kappahstore_crust_mantle,muvstore_crust_mantle, &
          muhstore_crust_mantle,eta_anisostore_crust_mantle, &
          c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
          c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
          c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
          c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
          c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
          c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
          c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
          ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
          R_memory_crust_mantle,epsilondev_crust_mantle, &
          eps_trace_over_3_crust_mantle,one_minus_sum_beta_crust_mantle, &
          alphaval,betaval,gammaval,factor_common_crust_mantle, &
          size(factor_common_crust_mantle,2), size(factor_common_crust_mantle,3), &
          size(factor_common_crust_mantle,4), size(factor_common_crust_mantle,5) )
      else
        call compute_forces_crust_mantle(minus_gravity_table,density_table,minus_deriv_gravity_table, &
          displ_crust_mantle,accel_crust_mantle, &
          xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
          xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
          etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
          gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
            is_on_a_slice_edge_crust_mantle,icall, &
            accel_inner_core,ibool_inner_core,idoubling_inner_core, &
            myrank,iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector,iphase, &
            nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
            npoin2D_cube_from_slices,buffer_all_cube_from_slices,buffer_slices,ibool_central_cube, &
            receiver_cube_from_slices,ibelm_bottom_inner_core,NSPEC2D_BOTTOM_IC,INCLUDE_CENTRAL_CUBE,iphase_CC, &
          hprime_xx,hprime_yy,hprime_zz, &
          hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
          kappavstore_crust_mantle,kappahstore_crust_mantle,muvstore_crust_mantle, &
          muhstore_crust_mantle,eta_anisostore_crust_mantle, &
          c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
          c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
          c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
          c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
          c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
          c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
          c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
          ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
          R_memory_crust_mantle,epsilondev_crust_mantle, &
          eps_trace_over_3_crust_mantle,one_minus_sum_beta_crust_mantle, &
          alphaval,betaval,gammaval,factor_common_crust_mantle, &
          size(factor_common_crust_mantle,2), size(factor_common_crust_mantle,3), &
          size(factor_common_crust_mantle,4), size(factor_common_crust_mantle,5) )
      endif

      ! Deville routine
      if( USE_DEVILLE_PRODUCTS_VAL ) then
        call compute_forces_inner_core_Dev(minus_gravity_table,density_table,minus_deriv_gravity_table, &
          displ_inner_core,accel_inner_core, &
          xstore_inner_core,ystore_inner_core,zstore_inner_core, &
          xix_inner_core,xiy_inner_core,xiz_inner_core, &
          etax_inner_core,etay_inner_core,etaz_inner_core, &
          gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
            is_on_a_slice_edge_inner_core,icall, &
            accel_crust_mantle,ibool_inner_core,idoubling_inner_core, &
            myrank,iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector,iphase, &
            nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
            npoin2D_cube_from_slices,buffer_all_cube_from_slices,buffer_slices,ibool_central_cube, &
            receiver_cube_from_slices,ibelm_bottom_inner_core,NSPEC2D_BOTTOM_IC,INCLUDE_CENTRAL_CUBE,iphase_CC, &
          hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
          kappavstore_inner_core,muvstore_inner_core,ibool_inner_core,idoubling_inner_core, &
          c11store_inner_core,c33store_inner_core,c12store_inner_core, &
          c13store_inner_core,c44store_inner_core, &
          R_memory_inner_core,epsilondev_inner_core, eps_trace_over_3_inner_core,&
          one_minus_sum_beta_inner_core, &
          alphaval,betaval,gammaval, &
          factor_common_inner_core, &
          size(factor_common_inner_core,2), size(factor_common_inner_core,3), &
          size(factor_common_inner_core,4), size(factor_common_inner_core,5) )
      else
        call compute_forces_inner_core(minus_gravity_table,density_table,minus_deriv_gravity_table, &
          displ_inner_core,accel_inner_core, &
          xstore_inner_core,ystore_inner_core,zstore_inner_core, &
          xix_inner_core,xiy_inner_core,xiz_inner_core, &
          etax_inner_core,etay_inner_core,etaz_inner_core, &
          gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
            is_on_a_slice_edge_inner_core,icall, &
            accel_crust_mantle,ibool_inner_core,idoubling_inner_core, &
            myrank,iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector,iphase, &
            nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
            npoin2D_cube_from_slices,buffer_all_cube_from_slices,buffer_slices,ibool_central_cube, &
            receiver_cube_from_slices,ibelm_bottom_inner_core,NSPEC2D_BOTTOM_IC,INCLUDE_CENTRAL_CUBE,iphase_CC, &
          hprime_xx,hprime_yy,hprime_zz,hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
          kappavstore_inner_core,muvstore_inner_core,ibool_inner_core,idoubling_inner_core, &
          c11store_inner_core,c33store_inner_core,c12store_inner_core, &
          c13store_inner_core,c44store_inner_core, &
          R_memory_inner_core,epsilondev_inner_core, eps_trace_over_3_inner_core,&
          one_minus_sum_beta_inner_core, &
          alphaval,betaval,gammaval, &
          factor_common_inner_core, &
          size(factor_common_inner_core,2), size(factor_common_inner_core,3), &
          size(factor_common_inner_core,4), size(factor_common_inner_core,5) )
      endif

! assemble all the contributions between slices using MPI
! crust/mantle and inner core handled in the same call
! in order to reduce the number of MPI messages by 2
      do while (iphase <= 7) ! make sure the last communications are finished and processed
        call assemble_MPI_vector(myrank,accel_crust_mantle,accel_inner_core, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector, &
            NUMMSGS_FACES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_CRUST_MANTLE), &
            NGLOB1D_RADIAL(IREGION_INNER_CORE),NCHUNKS_VAL,iphase)
      enddo

    !---
    !---  use buffers to assemble forces with the central cube
    !---

    if(INCLUDE_CENTRAL_CUBE) then
        do while (iphase_CC <= 4) ! make sure the last communications are finished and processed
          call assemble_MPI_central_cube(ichunk,nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
            npoin2D_cube_from_slices,buffer_all_cube_from_slices,buffer_slices,ibool_central_cube, &
            receiver_cube_from_slices,ibool_inner_core,idoubling_inner_core, &
            ibelm_bottom_inner_core,NSPEC2D_BOTTOM(IREGION_INNER_CORE),accel_inner_core,NDIM,iphase_CC)
        enddo
    endif   ! end of assembling forces with the central cube

    if(NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) then

       do i=1,NGLOB_CRUST_MANTLE
          accel_crust_mantle(1,i) = accel_crust_mantle(1,i)*rmassx_crust_mantle(i) &
               + two_omega_earth*veloc_crust_mantle(2,i)
          accel_crust_mantle(2,i) = accel_crust_mantle(2,i)*rmassy_crust_mantle(i) &
               - two_omega_earth*veloc_crust_mantle(1,i)
          accel_crust_mantle(3,i) = accel_crust_mantle(3,i)*rmassz_crust_mantle(i)
       enddo

    else

       do i=1,NGLOB_CRUST_MANTLE
          accel_crust_mantle(1,i) = accel_crust_mantle(1,i)*rmassz_crust_mantle(i) &
               + two_omega_earth*veloc_crust_mantle(2,i)
          accel_crust_mantle(2,i) = accel_crust_mantle(2,i)*rmassz_crust_mantle(i) &
               - two_omega_earth*veloc_crust_mantle(1,i)
          accel_crust_mantle(3,i) = accel_crust_mantle(3,i)*rmassz_crust_mantle(i)
       enddo

    endif

! ------------------- new non blocking implementation -------------------

!   if (SIMULATION_TYPE == 3) then

!     ! assemble all the contributions between slices using MPI

! assemble all the contributions between slices using MPI
! crust/mantle and inner core handled in the same call
! in order to reduce the number of MPI messages by 2

!       b_iphase = 1 ! initialize the non blocking communication counter
!       b_iphase_CC = 1 ! initialize the non blocking communication counter for the central cube

! start the non blocking communications
!       call assemble_MPI_vector(myrank,b_accel_crust_mantle,b_accel_inner_core, &
!           iproc_xi,iproc_eta,ichunk,addressing, &
!           iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
!           npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
!           iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
!           iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
!           npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
!           iboolfaces_inner_core,iboolcorner_inner_core, &
!           iprocfrom_faces,iprocto_faces, &
!           iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
!           b_buffer_send_faces,b_buffer_received_faces,npoin2D_max_all_CM_IC, &
!           b_buffer_send_chunkcorn_vector,b_buffer_recv_chunkcorn_vector, &
!           NUMMSGS_FACES,NCORNERSCHUNKS, &
!           NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_CRUST_MANTLE), &
!           NGLOB1D_RADIAL(IREGION_INNER_CORE),NCHUNKS_VAL,b_iphase)

!       b_icall = 2 ! now compute all the inner elements in the case of non blocking MPI

!       ! compute internal forces in the solid regions

!       ! for anisotropy and gravity, x y and z contain r theta and phi

!       if( USE_DEVILLE_PRODUCTS_VAL ) then
!         call compute_forces_crust_mantle_Dev(minus_gravity_table,density_table,minus_deriv_gravity_table, &
!         b_displ_crust_mantle,b_accel_crust_mantle, &
!         xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
!         xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
!         etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
!         gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
!           is_on_a_slice_edge_crust_mantle,b_icall, &
!           b_accel_inner_core,ibool_inner_core,idoubling_inner_core, &
!           myrank,iproc_xi,iproc_eta,ichunk,addressing, &
!           iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
!           npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
!           iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
!           iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
!           npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
!           iboolfaces_inner_core,iboolcorner_inner_core, &
!           iprocfrom_faces,iprocto_faces, &
!           iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
!           b_buffer_send_faces,b_buffer_received_faces,npoin2D_max_all_CM_IC, &
!           b_buffer_send_chunkcorn_vector,b_buffer_recv_chunkcorn_vector,b_iphase, &
!           nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
!           npoin2D_cube_from_slices,b_buffer_all_cube_from_slices,b_buffer_slices,ibool_central_cube, &
!           receiver_cube_from_slices,ibelm_bottom_inner_core,NSPEC2D_BOTTOM_IC,INCLUDE_CENTRAL_CUBE,b_iphase_CC, &
!         hprime_xx,hprime_xxT, &
!         hprimewgll_xx,hprimewgll_xxT, &
!         wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
!         kappavstore_crust_mantle,kappahstore_crust_mantle,muvstore_crust_mantle, &
!         muhstore_crust_mantle,eta_anisostore_crust_mantle, &
!         c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
!         c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
!         c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
!         c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
!         c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
!         c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
!         c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
!         ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
!         b_R_memory_crust_mantle,b_epsilondev_crust_mantle, &
!         b_eps_trace_over_3_crust_mantle,one_minus_sum_beta_crust_mantle, &
!         b_alphaval,b_betaval,b_gammaval,factor_common_crust_mantle, &
!         size(factor_common_crust_mantle,2), size(factor_common_crust_mantle,3), &
!         size(factor_common_crust_mantle,4), size(factor_common_crust_mantle,5) )
!       else
!         call compute_forces_crust_mantle(minus_gravity_table,density_table,minus_deriv_gravity_table, &
!         b_displ_crust_mantle,b_accel_crust_mantle, &
!         xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
!         xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
!         etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
!         gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
!           is_on_a_slice_edge_crust_mantle,b_icall, &
!           b_accel_inner_core,ibool_inner_core,idoubling_inner_core, &
!           myrank,iproc_xi,iproc_eta,ichunk,addressing, &
!           iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
!           npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
!           iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
!           iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
!           npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
!           iboolfaces_inner_core,iboolcorner_inner_core, &
!           iprocfrom_faces,iprocto_faces, &
!           iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
!           b_buffer_send_faces,b_buffer_received_faces,npoin2D_max_all_CM_IC, &
!           b_buffer_send_chunkcorn_vector,b_buffer_recv_chunkcorn_vector,b_iphase, &
!           nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
!           npoin2D_cube_from_slices,b_buffer_all_cube_from_slices,b_buffer_slices,ibool_central_cube, &
!           receiver_cube_from_slices,ibelm_bottom_inner_core,NSPEC2D_BOTTOM_IC,INCLUDE_CENTRAL_CUBE,b_iphase_CC, &
!         hprime_xx,hprime_yy,hprime_zz, &
!         hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
!         wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
!         kappavstore_crust_mantle,kappahstore_crust_mantle,muvstore_crust_mantle, &
!         muhstore_crust_mantle,eta_anisostore_crust_mantle, &
!         c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
!         c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
!         c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
!         c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
!         c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
!         c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
!         c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
!         ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
!         b_R_memory_crust_mantle,b_epsilondev_crust_mantle, &
!         b_eps_trace_over_3_crust_mantle,one_minus_sum_beta_crust_mantle, &
!         b_alphaval,b_betaval,b_gammaval,factor_common_crust_mantle, &
!         size(factor_common_crust_mantle,2), size(factor_common_crust_mantle,3), &
!         size(factor_common_crust_mantle,4), size(factor_common_crust_mantle,5) )
!       endif

!       ! Deville routine
!       if( USE_DEVILLE_PRODUCTS_VAL ) then
!         call compute_forces_inner_core_Dev(minus_gravity_table,density_table,minus_deriv_gravity_table, &
!         b_displ_inner_core,b_accel_inner_core, &
!         xstore_inner_core,ystore_inner_core,zstore_inner_core, &
!         xix_inner_core,xiy_inner_core,xiz_inner_core, &
!         etax_inner_core,etay_inner_core,etaz_inner_core, &
!         gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
!           is_on_a_slice_edge_inner_core,b_icall, &
!           b_accel_crust_mantle,ibool_inner_core,idoubling_inner_core, &
!           myrank,iproc_xi,iproc_eta,ichunk,addressing, &
!           iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
!           npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
!           iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
!           iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
!           npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
!           iboolfaces_inner_core,iboolcorner_inner_core, &
!           iprocfrom_faces,iprocto_faces, &
!           iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
!           b_buffer_send_faces,b_buffer_received_faces,npoin2D_max_all_CM_IC, &
!           b_buffer_send_chunkcorn_vector,b_buffer_recv_chunkcorn_vector,b_iphase, &
!           nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
!           npoin2D_cube_from_slices,b_buffer_all_cube_from_slices,b_buffer_slices,ibool_central_cube, &
!           receiver_cube_from_slices,ibelm_bottom_inner_core,NSPEC2D_BOTTOM_IC,INCLUDE_CENTRAL_CUBE,b_iphase_CC, &
!         hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
!         wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
!         kappavstore_inner_core,muvstore_inner_core,ibool_inner_core,idoubling_inner_core, &
!         c11store_inner_core,c33store_inner_core,c12store_inner_core, &
!         c13store_inner_core,c44store_inner_core, &
!         b_R_memory_inner_core,b_epsilondev_inner_core, b_eps_trace_over_3_inner_core,&
!         one_minus_sum_beta_inner_core, &
!         b_alphaval,b_betaval,b_gammaval, &
!         factor_common_inner_core, &
!         size(factor_common_inner_core,2), size(factor_common_inner_core,3), &
!         size(factor_common_inner_core,4), size(factor_common_inner_core,5) )
!       else
!         call compute_forces_inner_core(minus_gravity_table,density_table,minus_deriv_gravity_table, &
!         b_displ_inner_core,b_accel_inner_core, &
!         xstore_inner_core,ystore_inner_core,zstore_inner_core, &
!         xix_inner_core,xiy_inner_core,xiz_inner_core, &
!         etax_inner_core,etay_inner_core,etaz_inner_core, &
!         gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
!           is_on_a_slice_edge_inner_core,b_icall, &
!           b_accel_crust_mantle,ibool_inner_core,idoubling_inner_core, &
!           myrank,iproc_xi,iproc_eta,ichunk,addressing, &
!           iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
!           npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
!           iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
!           iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
!           npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
!           iboolfaces_inner_core,iboolcorner_inner_core, &
!           iprocfrom_faces,iprocto_faces, &
!           iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
!           b_buffer_send_faces,b_buffer_received_faces,npoin2D_max_all_CM_IC, &
!           b_buffer_send_chunkcorn_vector,b_buffer_recv_chunkcorn_vector,b_iphase, &
!           nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
!           npoin2D_cube_from_slices,b_buffer_all_cube_from_slices,b_buffer_slices,ibool_central_cube, &
!           receiver_cube_from_slices,ibelm_bottom_inner_core,NSPEC2D_BOTTOM_IC,INCLUDE_CENTRAL_CUBE,b_iphase_CC, &
!         hprime_xx,hprime_yy,hprime_zz,hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
!         wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
!         kappavstore_inner_core,muvstore_inner_core,ibool_inner_core,idoubling_inner_core, &
!         c11store_inner_core,c33store_inner_core,c12store_inner_core, &
!         c13store_inner_core,c44store_inner_core, &
!         b_R_memory_inner_core,b_epsilondev_inner_core, b_eps_trace_over_3_inner_core,&
!         one_minus_sum_beta_inner_core, &
!         b_alphaval,b_betaval,b_gammaval, &
!         factor_common_inner_core, &
!         size(factor_common_inner_core,2), size(factor_common_inner_core,3), &
!         size(factor_common_inner_core,4), size(factor_common_inner_core,5) )
!       endif

! assemble all the contributions between slices using MPI
! crust/mantle and inner core handled in the same call
! in order to reduce the number of MPI messages by 2
!       do while (b_iphase <= 7) ! make sure the last communications are finished and processed
!         call assemble_MPI_vector(myrank,b_accel_crust_mantle,b_accel_inner_core, &
!           iproc_xi,iproc_eta,ichunk,addressing, &
!           iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
!           npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
!           iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
!           iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
!           npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
!           iboolfaces_inner_core,iboolcorner_inner_core, &
!           iprocfrom_faces,iprocto_faces, &
!           iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
!           b_buffer_send_faces,b_buffer_received_faces,npoin2D_max_all_CM_IC, &
!           b_buffer_send_chunkcorn_vector,b_buffer_recv_chunkcorn_vector, &
!           NUMMSGS_FACES,NCORNERSCHUNKS, &
!           NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_CRUST_MANTLE), &
!           NGLOB1D_RADIAL(IREGION_INNER_CORE),NCHUNKS_VAL,b_iphase)
!       enddo

!     !---
!     !---  use buffers to assemble forces with the central cube
!     !---

!     if(INCLUDE_CENTRAL_CUBE) then
!         do while (b_iphase_CC <= 4) ! make sure the last communications are finished and processed
!           call assemble_MPI_central_cube(ichunk,nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
!             npoin2D_cube_from_slices,b_buffer_all_cube_from_slices,b_buffer_slices,ibool_central_cube, &
!             receiver_cube_from_slices,ibool_inner_core,idoubling_inner_core, &
!             ibelm_bottom_inner_core,NSPEC2D_BOTTOM(IREGION_INNER_CORE),b_accel_inner_core,NDIM,b_iphase_CC)
!         enddo
!     endif   ! end of assembling forces with the central cube

! ------------------- new non blocking implementation -------------------

!     if(NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) then

!        do i=1,NGLOB_CRUST_MANTLE
!           b_accel_crust_mantle(1,i) = b_accel_crust_mantle(1,i)*rmassx_crust_mantle(i) &
!                + b_two_omega_earth*b_veloc_crust_mantle(2,i)
!           b_accel_crust_mantle(2,i) = b_accel_crust_mantle(2,i)*rmassy_crust_mantle(i) &
!                - b_two_omega_earth*b_veloc_crust_mantle(1,i)
!           b_accel_crust_mantle(3,i) = b_accel_crust_mantle(3,i)*rmassz_crust_mantle(i)
!        enddo

!     else

!        do i=1,NGLOB_CRUST_MANTLE
!           b_accel_crust_mantle(1,i) = b_accel_crust_mantle(1,i)*rmassz_crust_mantle(i) &
!                + b_two_omega_earth*b_veloc_crust_mantle(2,i)
!           b_accel_crust_mantle(2,i) = b_accel_crust_mantle(2,i)*rmassz_crust_mantle(i) &
!                - b_two_omega_earth*b_veloc_crust_mantle(1,i)
!           b_accel_crust_mantle(3,i) = b_accel_crust_mantle(3,i)*rmassz_crust_mantle(i)
!        enddo

!     endif

!  endif ! SIMULATION_TYPE == 3

    ! couples ocean with crust mantle
   if(OCEANS_VAL) then
     call compute_coupling_ocean(accel_crust_mantle, &
                                   rmassx_crust_mantle,rmassy_crust_mantle,rmassz_crust_mantle, &
                                   rmass_ocean_load,normal_top_crust_mantle, &
                                   ibool_crust_mantle,ibelm_top_crust_mantle, &
                                   updated_dof_ocean_load,NGLOB_XY, &
                                   NSPEC2D_TOP(IREGION_CRUST_MANTLE), &
                                   ABSORBING_CONDITIONS)
!    if(SIMULATION_TYPE == 3) then
!      call compute_coupling_ocean(b_accel_crust_mantle, &
!                                  rmassx_crust_mantle,rmassy_crust_mantle,rmassz_crust_mantle, &
!                                  rmass_ocean_load,normal_top_crust_mantle, &
!                                  ibool_crust_mantle,ibelm_top_crust_mantle, &
!                                  updated_dof_ocean_load,NGLOB_XY, &
!                                  NSPEC2D_TOP(IREGION_CRUST_MANTLE), &
!                                  ABSORBING_CONDITIONS)
!    endif
   endif

!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------

    ! Newmark time scheme - corrector for elastic parts

    ! mantle
    do i=1,NGLOB_CRUST_MANTLE
      veloc_crust_mantle(:,i) = veloc_crust_mantle(:,i) + deltatover2*accel_crust_mantle(:,i)
    enddo
    ! inner core
    do i=1,NGLOB_INNER_CORE
      accel_inner_core(1,i) = accel_inner_core(1,i)*rmass_inner_core(i) &
             + two_omega_earth*veloc_inner_core(2,i)
      accel_inner_core(2,i) = accel_inner_core(2,i)*rmass_inner_core(i) &
             - two_omega_earth*veloc_inner_core(1,i)
      accel_inner_core(3,i) = accel_inner_core(3,i)*rmass_inner_core(i)

      veloc_inner_core(:,i) = veloc_inner_core(:,i) + deltatover2*accel_inner_core(:,i)
    enddo

!   if (SIMULATION_TYPE == 3) then
!     ! mantle
!     do i=1,NGLOB_CRUST_MANTLE
!       b_veloc_crust_mantle(:,i) = b_veloc_crust_mantle(:,i) + b_deltatover2*b_accel_crust_mantle(:,i)
!     enddo
!     ! inner core
!     do i=1,NGLOB_INNER_CORE
!       b_accel_inner_core(1,i) = b_accel_inner_core(1,i)*rmass_inner_core(i) &
!        + b_two_omega_earth*b_veloc_inner_core(2,i)
!       b_accel_inner_core(2,i) = b_accel_inner_core(2,i)*rmass_inner_core(i) &
!        - b_two_omega_earth*b_veloc_inner_core(1,i)
!       b_accel_inner_core(3,i) = b_accel_inner_core(3,i)*rmass_inner_core(i)
!       b_veloc_inner_core(:,i) = b_veloc_inner_core(:,i) + b_deltatover2*b_accel_inner_core(:,i)
!     enddo
!   endif ! SIMULATION_TYPE == 3

!   ! restores last time snapshot saved for backward/reconstruction of wavefields
!   ! note: this is done here after the Newmark time scheme, otherwise the indexing for sources
!   !          and adjoint sources will become more complicated
!   !          that is, index it for adjoint sources will match index NSTEP - 1 for backward/reconstructed wavefields
!   if(SIMULATION_TYPE == 3 .and. it == 1) then
!     call read_forward_arrays(myrank, &
!                   b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle, &
!                   b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core, &
!                   b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core, &
!                   b_R_memory_crust_mantle,b_R_memory_inner_core, &
!                   b_epsilondev_crust_mantle,b_epsilondev_inner_core, &
!                   b_A_array_rotation,b_B_array_rotation,LOCAL_PATH)
!   endif

! write the seismograms with time shift

! store the seismograms only if there is at least one receiver located in this slice
  if (nrec_local > 0) then
    if (SIMULATION_TYPE == 1) then
      call compute_seismograms(nrec_local,nrec,displ_crust_mantle, &
                                nu,hxir_store,hetar_store,hgammar_store, &
                                scale_displ,ibool_crust_mantle, &
                                ispec_selected_rec,number_receiver_global, &
                                seismo_current,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
                                seismograms)

    else if (SIMULATION_TYPE == 2) then
      call compute_seismograms_adjoint(NSOURCES,nrec_local,displ_crust_mantle, &
                    eps_trace_over_3_crust_mantle,epsilondev_crust_mantle, &
                    nu_source,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                    hxir_store,hetar_store,hgammar_store, &
                    hpxir_store,hpetar_store,hpgammar_store, &
                    tshift_cmt,hdur_gaussian,DT,t0,scale_displ, &
                    hprime_xx,hprime_yy,hprime_zz, &
                    xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                    etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                    gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                    moment_der,sloc_der,stshift_der,shdur_der, &
                    NTSTEP_BETWEEN_OUTPUT_SEISMOS,seismograms,deltat, &
                    ibool_crust_mantle,ispec_selected_source,number_receiver_global, &
                    NSTEP,it,nit_written)

!   else if (SIMULATION_TYPE == 3) then
!     call compute_seismograms_backward(nrec_local,nrec,b_displ_crust_mantle, &
!                               nu,hxir_store,hetar_store,hgammar_store, &
!                               scale_displ,ibool_crust_mantle, &
!                               ispec_selected_rec,number_receiver_global, &
!                               seismo_current,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
!                               seismograms)

    endif
  endif ! nrec_local

  ! write the current or final seismograms
  if(seismo_current == NTSTEP_BETWEEN_OUTPUT_SEISMOS .or. it == it_end) then
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      call write_seismograms(myrank,seismograms,number_receiver_global,station_name, &
            network_name,stlat,stlon,stele,stbur, &
            nrec,nrec_local,ANGULAR_WIDTH_XI_IN_DEGREES,NEX_XI,DT,t0,it_end, &
            yr_SAC,jda_SAC,ho_SAC,mi_SAC,sec_SAC,t_cmt_SAC,t_shift_SAC, &
            elat_SAC,elon_SAC,depth_SAC,event_name_SAC,cmt_lat_SAC,cmt_lon_SAC,&
            cmt_depth_SAC,cmt_hdur_SAC,NPROCTOT_VAL, &
            OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM, &
            OUTPUT_SEISMOS_SAC_BINARY,ROTATE_SEISMOGRAMS_RT,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
            seismo_offset,seismo_current,WRITE_SEISMOGRAMS_BY_MASTER, &
            SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,MODEL)
      if(myrank==0) then
        write(IMAIN,*)
        write(IMAIN,*) ' Total number of time steps written: ', it-it_begin+1
        write(IMAIN,*)
      endif
    else ! case of SIMULATION_TYPE == 2
      if( nrec_local > 0 ) &
        call write_adj_seismograms(seismograms,number_receiver_global, &
                                  nrec_local,it,nit_written,DT, &
                                  NSTEP,NTSTEP_BETWEEN_OUTPUT_SEISMOS,t0,LOCAL_PATH)
        nit_written = it
    endif
    seismo_offset = seismo_offset + seismo_current
    seismo_current = 0
  endif

!
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!

! kernel calculations
! if (SIMULATION_TYPE == 3) then
!   ! crust mantle
!   call compute_kernels_crust_mantle(ibool_crust_mantle, &
!                         rho_kl_crust_mantle,beta_kl_crust_mantle, &
!                         alpha_kl_crust_mantle,cijkl_kl_crust_mantle, &
!                         accel_crust_mantle,b_displ_crust_mantle, &
!                         epsilondev_crust_mantle,b_epsilondev_crust_mantle, &
!                         eps_trace_over_3_crust_mantle,b_eps_trace_over_3_crust_mantle, &
!                         deltat)

!   ! outer core
!   call compute_kernels_outer_core(ibool_outer_core, &
!                       xix_outer_core,xiy_outer_core,xiz_outer_core, &
!                       etax_outer_core,etay_outer_core,etaz_outer_core, &
!                       gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
!                       hprime_xx,hprime_yy,hprime_zz, &
!                       displ_outer_core,accel_outer_core, &
!                       b_displ_outer_core,b_accel_outer_core, &
!                       vector_accel_outer_core,vector_displ_outer_core, &
!                       b_vector_displ_outer_core, &
!                       div_displ_outer_core,b_div_displ_outer_core, &
!                       rhostore_outer_core,kappavstore_outer_core, &
!                       rho_kl_outer_core,alpha_kl_outer_core, &
!                       deviatoric_outercore,nspec_beta_kl_outer_core,beta_kl_outer_core, &
!                       deltat)

!   ! inner core
!   call compute_kernels_inner_core(ibool_inner_core, &
!                         rho_kl_inner_core,beta_kl_inner_core, &
!                         alpha_kl_inner_core, &
!                         accel_inner_core,b_displ_inner_core, &
!                         epsilondev_inner_core,b_epsilondev_inner_core, &
!                         eps_trace_over_3_inner_core,b_eps_trace_over_3_inner_core, &
!                         deltat)

!   ! NOISE TOMOGRAPHY --- source strength kernel
!   if (NOISE_TOMOGRAPHY == 3)  &
!      call compute_kernels_strength_noise(nmovie_points,ibool_crust_mantle, &
!                         Sigma_kl_crust_mantle,displ_crust_mantle,deltat,it, &
!                         normal_x_noise,normal_y_noise,normal_z_noise, &
!                         NSPEC2D_TOP(IREGION_CRUST_MANTLE),noise_surface_movie, &
!                         ibelm_top_crust_mantle)

!   ! --- boundary kernels ------
!   if (SAVE_BOUNDARY_MESH) then
!     fluid_solid_boundary = .false.
!     iregion_code = IREGION_CRUST_MANTLE

!     ! Moho
!     if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
!       call compute_boundary_kernel(displ_crust_mantle,accel_crust_mantle, &
!                b_displ_crust_mantle,nspec_crust_mantle,iregion_code, &
!                ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
!                xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
!                etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
!                gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,hprime_xx,hprime_yy,hprime_zz, &
!                rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
!                kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
!                c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
!                c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
!                c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
!                c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
!                c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
!                c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
!                k_top,ibelm_moho_top,normal_moho,moho_kl_top,fluid_solid_boundary,NSPEC2D_MOHO)

!       call compute_boundary_kernel(displ_crust_mantle,accel_crust_mantle, &
!                b_displ_crust_mantle,nspec_crust_mantle,iregion_code, &
!                ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
!                xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
!                etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
!                gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,hprime_xx,hprime_yy,hprime_zz, &
!                rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
!                kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
!                c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
!                c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
!                c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
!                c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
!                c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
!                c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
!                k_bot,ibelm_moho_bot,normal_moho,moho_kl_bot,fluid_solid_boundary,NSPEC2D_MOHO)

!       moho_kl = moho_kl + (moho_kl_top - moho_kl_bot) * deltat
!     endif

!     ! d400
!     call compute_boundary_kernel(displ_crust_mantle,accel_crust_mantle, &
!                b_displ_crust_mantle,nspec_crust_mantle,iregion_code, &
!                ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
!                xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
!                etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
!                gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,hprime_xx,hprime_yy,hprime_zz, &
!                rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
!                kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
!                c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
!                c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
!                c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
!                c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
!                c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
!                c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
!                k_top,ibelm_400_top,normal_400,d400_kl_top,fluid_solid_boundary,NSPEC2D_400)

!     call compute_boundary_kernel(displ_crust_mantle,accel_crust_mantle, &
!                b_displ_crust_mantle,nspec_crust_mantle,iregion_code, &
!                ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
!                xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
!                etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
!                gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,hprime_xx,hprime_yy,hprime_zz, &
!                rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
!                kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
!                c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
!                c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
!                c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
!                c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
!                c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
!                c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
!                k_bot,ibelm_400_bot,normal_400,d400_kl_bot,fluid_solid_boundary,NSPEC2D_400)

!     d400_kl = d400_kl + (d400_kl_top - d400_kl_bot) * deltat

!     ! d670
!     call compute_boundary_kernel(displ_crust_mantle,accel_crust_mantle, &
!                b_displ_crust_mantle,nspec_crust_mantle,iregion_code, &
!                ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
!                xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
!                etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
!                gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,hprime_xx,hprime_yy,hprime_zz, &
!                rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
!                kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
!                c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
!                c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
!                c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
!                c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
!                c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
!                c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
!                k_top,ibelm_670_top,normal_670,d670_kl_top,fluid_solid_boundary,NSPEC2D_670)

!     call compute_boundary_kernel(displ_crust_mantle,accel_crust_mantle, &
!                b_displ_crust_mantle,nspec_crust_mantle,iregion_code, &
!                ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
!                xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
!                etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
!                gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,hprime_xx,hprime_yy,hprime_zz, &
!                rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
!                kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
!                c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
!                c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
!                c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
!                c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
!                c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
!                c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
!                k_bot,ibelm_670_bot,normal_670,d670_kl_bot,fluid_solid_boundary,NSPEC2D_670)

!     d670_kl = d670_kl + (d670_kl_top - d670_kl_bot) * deltat

!     ! CMB
!     fluid_solid_boundary = .true.
!     iregion_code = IREGION_CRUST_MANTLE
!     call compute_boundary_kernel(displ_crust_mantle,accel_crust_mantle, &
!                b_displ_crust_mantle,nspec_crust_mantle,iregion_code, &
!                ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
!                xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
!                etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
!                gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,hprime_xx,hprime_yy,hprime_zz, &
!                rhostore_crust_mantle,kappavstore_crust_mantle, muvstore_crust_mantle, &
!                kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
!                c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
!                c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
!                c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
!                c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
!                c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
!                c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
!                k_top,ibelm_bottom_crust_mantle,normal_top_outer_core, &
!                cmb_kl_top,fluid_solid_boundary,NSPEC2D_CMB)

!     iregion_code = IREGION_OUTER_CORE
!     call compute_boundary_kernel(vector_displ_outer_core,vector_accel_outer_core, &
!                b_vector_displ_outer_core,nspec_outer_core, &
!                iregion_code,ystore_outer_core,zstore_outer_core,ibool_outer_core,ispec_is_tiso_outer_core, &
!                xix_outer_core,xiy_outer_core,xiz_outer_core, &
!                etax_outer_core,etay_outer_core,etaz_outer_core,&
!                gammax_outer_core,gammay_outer_core,gammaz_outer_core,hprime_xx,hprime_yy,hprime_zz, &
!                rhostore_outer_core,kappavstore_outer_core,dummy_array, &
!                dummy_array,dummy_array,dummy_array, &
!                dummy_array,dummy_array,dummy_array,dummy_array, &
!                dummy_array,dummy_array,dummy_array, &
!                dummy_array,dummy_array,dummy_array,dummy_array, &
!                dummy_array,dummy_array,dummy_array, &
!                dummy_array,dummy_array,dummy_array,dummy_array, &
!                dummy_array,dummy_array,dummy_array, &
!                k_bot,ibelm_top_outer_core,normal_top_outer_core, &
!                cmb_kl_bot,fluid_solid_boundary,NSPEC2D_CMB)

!     cmb_kl = cmb_kl + (cmb_kl_top - cmb_kl_bot) * deltat

!     ! ICB
!     fluid_solid_boundary = .true.
!     call compute_boundary_kernel(vector_displ_outer_core,vector_accel_outer_core, &
!                b_vector_displ_outer_core,nspec_outer_core, &
!                iregion_code,ystore_outer_core,zstore_outer_core,ibool_outer_core,ispec_is_tiso_outer_core, &
!                xix_outer_core,xiy_outer_core,xiz_outer_core, &
!                etax_outer_core,etay_outer_core,etaz_outer_core,&
!                gammax_outer_core,gammay_outer_core,gammaz_outer_core,hprime_xx,hprime_yy,hprime_zz, &
!                rhostore_outer_core,kappavstore_outer_core,dummy_array, &
!                dummy_array,dummy_array,dummy_array, &
!                dummy_array,dummy_array,dummy_array,dummy_array, &
!                dummy_array,dummy_array,dummy_array, &
!                dummy_array,dummy_array,dummy_array,dummy_array, &
!                dummy_array,dummy_array,dummy_array, &
!                dummy_array,dummy_array,dummy_array,dummy_array, &
!                dummy_array,dummy_array,dummy_array, &
!                k_top,ibelm_bottom_outer_core,normal_bottom_outer_core, &
!                icb_kl_top,fluid_solid_boundary,NSPEC2D_ICB)

!     iregion_code = IREGION_INNER_CORE
!     call compute_boundary_kernel(displ_inner_core,accel_inner_core, &
!                b_displ_inner_core,nspec_inner_core,iregion_code, &
!                ystore_inner_core,zstore_inner_core,ibool_inner_core,ispec_is_tiso_inner_core, &
!                xix_inner_core,xiy_inner_core,xiz_inner_core, &
!                etax_inner_core,etay_inner_core,etaz_inner_core,&
!                gammax_inner_core,gammay_inner_core,gammaz_inner_core,hprime_xx,hprime_yy,hprime_zz, &
!                rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core, &
!                dummy_array,dummy_array,dummy_array, &
!                c11store_inner_core,c12store_inner_core,c13store_inner_core,dummy_array, &
!                dummy_array,dummy_array,dummy_array, &
!                dummy_array,dummy_array,dummy_array,dummy_array, &
!                c33store_inner_core,dummy_array,dummy_array, &
!                dummy_array,c44store_inner_core,dummy_array,dummy_array, &
!                dummy_array,dummy_array,dummy_array, &
!                k_bot,ibelm_top_inner_core,normal_bottom_outer_core, &
!                icb_kl_bot,fluid_solid_boundary,NSPEC2D_ICB)

!     icb_kl = icb_kl + (icb_kl_top - icb_kl_bot) * deltat
!   endif

!   ! approximate hessian
!   if( APPROXIMATE_HESS_KL ) then
!     call compute_kernels_hessian(ibool_crust_mantle, &
!                         hess_kl_crust_mantle,&
!                         accel_crust_mantle,b_accel_crust_mantle, &
!                         deltat)
!   endif

! endif ! end of if computing kernels

!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------

  ! first step of noise tomography, i.e., save a surface movie at every time step
  ! modified from the subroutine 'write_movie_surface'
  if ( NOISE_TOMOGRAPHY == 1 ) then
        call noise_save_surface_movie(displ_crust_mantle, &
                            ibelm_top_crust_mantle,ibool_crust_mantle, &
                            NSPEC2D_TOP(IREGION_CRUST_MANTLE),noise_surface_movie,it)
  endif

  ! save movie on surface
  if( MOVIE_SURFACE ) then
    if( mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then
      ! save velocity here to avoid static offset on displacement for movies
      call write_movie_surface(myrank,nmovie_points,scale_veloc,veloc_crust_mantle, &
                    scale_displ,displ_crust_mantle, &
                    xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                    store_val_x,store_val_y,store_val_z, &
                    store_val_x_all,store_val_y_all,store_val_z_all, &
                    store_val_ux,store_val_uy,store_val_uz, &
                    store_val_ux_all,store_val_uy_all,store_val_uz_all, &
                    ibelm_top_crust_mantle,ibool_crust_mantle, &
                    NSPEC2D_TOP(IREGION_CRUST_MANTLE), &
                    NIT,it,OUTPUT_FILES,MOVIE_VOLUME_TYPE)
    endif
  endif


  ! save movie in full 3D mesh
  if(MOVIE_VOLUME ) then
    if( mod(it-MOVIE_START,NTSTEP_BETWEEN_FRAMES) == 0  &
      .and. it >= MOVIE_START .and. it <= MOVIE_STOP) then

      if (MOVIE_VOLUME_TYPE == 1) then  ! output strains
        call  write_movie_volume_strains(myrank,npoints_3dmovie, &
                    LOCAL_PATH,MOVIE_VOLUME_TYPE,MOVIE_COARSE, &
                    it,eps_trace_over_3_crust_mantle,epsilondev_crust_mantle, &
                    muvstore_crust_mantle_3dmovie, &
                    mask_3dmovie,nu_3dmovie)

      else if (MOVIE_VOLUME_TYPE == 2 .or. MOVIE_VOLUME_TYPE == 3) then
        ! output the Time Integral of Strain, or \mu*TIS
        call  write_movie_volume_strains(myrank,npoints_3dmovie, &
                    LOCAL_PATH,MOVIE_VOLUME_TYPE,MOVIE_COARSE, &
                    it,Ieps_trace_over_3_crust_mantle,Iepsilondev_crust_mantle, &
                    muvstore_crust_mantle_3dmovie, &
                    mask_3dmovie,nu_3dmovie)

      else if (MOVIE_VOLUME_TYPE == 4) then ! output divergence and curl in whole volume
        call write_movie_volume_divcurl(myrank,it,eps_trace_over_3_crust_mantle,&
                        div_displ_outer_core, &
                        accel_outer_core,kappavstore_outer_core,rhostore_outer_core,ibool_outer_core, &
                        eps_trace_over_3_inner_core, &
                        epsilondev_crust_mantle,epsilondev_inner_core, &
                        LOCAL_PATH, &
                        displ_crust_mantle,displ_inner_core,displ_outer_core, &
                        veloc_crust_mantle,veloc_inner_core,veloc_outer_core, &
                        accel_crust_mantle,accel_inner_core, &
                        ibool_crust_mantle,ibool_inner_core)

      else if (MOVIE_VOLUME_TYPE == 5) then ! output displacement
        scalingval = scale_displ
        call write_movie_volume_vector(myrank,it,npoints_3dmovie, &
                    LOCAL_PATH,MOVIE_VOLUME_TYPE, &
                    MOVIE_COARSE,ibool_crust_mantle,displ_crust_mantle, &
                    scalingval,mask_3dmovie,nu_3dmovie)

      else if (MOVIE_VOLUME_TYPE == 6) then ! output velocity
        scalingval = scale_veloc
        call write_movie_volume_vector(myrank,it,npoints_3dmovie, &
                    LOCAL_PATH,MOVIE_VOLUME_TYPE, &
                    MOVIE_COARSE,ibool_crust_mantle,veloc_crust_mantle, &
                    scalingval,mask_3dmovie,nu_3dmovie)

      else

        call exit_MPI(myrank, 'MOVIE_VOLUME_TYPE has to be 1,2,3,4,5 or 6')

      endif ! MOVIE_VOLUME_TYPE
    endif
  endif ! MOVIE_VOLUME

