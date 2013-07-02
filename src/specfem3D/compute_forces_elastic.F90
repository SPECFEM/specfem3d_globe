!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
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

  subroutine compute_forces_elastic()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore,only: accel_outer_core,b_accel_outer_core, &
                                  normal_top_outer_core,jacobian2D_top_outer_core, &
                                  normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
                                  ibelm_top_outer_core,ibelm_bottom_outer_core, &
                                  ibool_outer_core
  use specfem_par_movie
  implicit none

  ! local parameters
  ! non blocking MPI
  ! iphase: iphase = 1 is for computing outer elements in the crust_mantle and inner_core regions,
  !              iphase = 2 is for computing inner elements (former icall parameter)
  integer :: iphase
  logical :: phase_is_inner


!daniel: att - debug
!  integer :: iglob
!  logical,parameter :: DEBUG = .false.
!  if( DEBUG ) then
!    iglob = ibool_crust_mantle(1,1,1,100)
!    if( SIMULATION_TYPE == 1) then
!      if( it == NSTEP .and. myrank == 0 ) then
!        print*,'last step',it,'Rxx:',R_xx_crust_mantle(1,1,1,1,100),epsilondev_xx_crust_mantle(1,1,1,100), &
!          displ_crust_mantle(1,iglob),displ_crust_mantle(2,iglob),displ_crust_mantle(3,iglob)
!      endif
!      if( it == NSTEP-1 .and. myrank == 0 ) then
!        print*,'second last step',it,'Rxx:',R_xx_crust_mantle(1,1,1,1,100),epsilondev_xx_crust_mantle(1,1,1,100), &
!          displ_crust_mantle(1,iglob),displ_crust_mantle(2,iglob),displ_crust_mantle(3,iglob)
!      endif
!      if( it == NSTEP-2 .and. myrank == 0 ) then
!        print*,'third last step',it,'Rxx:',R_xx_crust_mantle(1,1,1,1,100),epsilondev_xx_crust_mantle(1,1,1,100), &
!          displ_crust_mantle(1,iglob),displ_crust_mantle(2,iglob),displ_crust_mantle(3,iglob)
!      endif
!    else if( SIMULATION_TYPE == 3 ) then
!      if( it == 1 .and. myrank == 0 ) then
!        print*,'first step',it,'Rxx:',b_R_xx_crust_mantle(1,1,1,1,100),b_epsilondev_xx_crust_mantle(1,1,1,100), &
!          b_displ_crust_mantle(1,iglob),b_displ_crust_mantle(2,iglob),b_displ_crust_mantle(3,iglob)
!      endif
!      if( it == 2 .and. myrank == 0 ) then
!        print*,'second step',it,'Rxx:',b_R_xx_crust_mantle(1,1,1,1,100),b_epsilondev_xx_crust_mantle(1,1,1,100), &
!          b_displ_crust_mantle(1,iglob),b_displ_crust_mantle(2,iglob),b_displ_crust_mantle(3,iglob)
!      endif
!      if( it == 3 .and. myrank == 0 ) then
!        print*,'third step',it,'Rxx:',b_R_xx_crust_mantle(1,1,1,1,100),b_epsilondev_xx_crust_mantle(1,1,1,100), &
!          b_displ_crust_mantle(1,iglob),b_displ_crust_mantle(2,iglob),b_displ_crust_mantle(3,iglob)
!      endif
!    endif
!  endif


  ! ****************************************************
  !   big loop over all spectral elements in the solid
  ! ****************************************************

  ! compute internal forces in the solid regions

  ! for anisotropy and gravity, x y and z contain r theta and phi

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase = 1,2

    ! first, iphase == 1 for points on MPI interfaces (thus outer elements)
    ! second, iphase == 2 for points purely inside partition (thus inner elements)
    !
    ! compute all the outer elements first, then sends out non blocking MPI communication
    ! and continues computing inner elements (overlapping)
    if( iphase == 1 ) then
      phase_is_inner = .false.
    else
      phase_is_inner = .true.
    endif


    if( .NOT. GPU_MODE ) then
       ! on CPU

       ! compute internal forces in the solid regions
       ! note: for anisotropy and gravity, x y and z contain r theta and phi
       if( USE_DEVILLE_PRODUCTS_VAL ) then
          ! uses Deville (2002) optimizations
          ! crust/mantle region
          call compute_forces_crust_mantle_Dev( NSPEC_CRUST_MANTLE_STR_OR_ATT,NGLOB_CRUST_MANTLE, &
               NSPEC_CRUST_MANTLE_ATTENUAT, &
               deltat, &
               displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle, &
               phase_is_inner, &
               R_xx_crust_mantle,R_yy_crust_mantle,R_xy_crust_mantle, &
               R_xz_crust_mantle,R_yz_crust_mantle, &
               epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
               epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
               eps_trace_over_3_crust_mantle, &
               alphaval,betaval,gammaval,factor_common_crust_mantle, &
               size(factor_common_crust_mantle,2), size(factor_common_crust_mantle,3), &
               size(factor_common_crust_mantle,4), size(factor_common_crust_mantle,5), .false. )
          ! inner core region
          call compute_forces_inner_core_Dev( NSPEC_INNER_CORE_STR_OR_ATT,NGLOB_INNER_CORE, &
               NSPEC_INNER_CORE_ATTENUATION, &
               deltat, &
               displ_inner_core,veloc_inner_core,accel_inner_core, &
               phase_is_inner, &
               R_xx_inner_core,R_yy_inner_core,R_xy_inner_core,R_xz_inner_core,R_yz_inner_core, &
               epsilondev_xx_inner_core,epsilondev_yy_inner_core,epsilondev_xy_inner_core, &
               epsilondev_xz_inner_core,epsilondev_yz_inner_core, &
               eps_trace_over_3_inner_core,&
               alphaval,betaval,gammaval, &
               factor_common_inner_core, &
               size(factor_common_inner_core,2), size(factor_common_inner_core,3), &
               size(factor_common_inner_core,4), size(factor_common_inner_core,5), .false. )

       else
          ! no Deville optimization
          ! crust/mantle region
          call compute_forces_crust_mantle(  NSPEC_CRUST_MANTLE_STR_OR_ATT,NGLOB_CRUST_MANTLE, &
               NSPEC_CRUST_MANTLE_ATTENUAT, &
               deltat, &
               displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle, &
               phase_is_inner, &
               R_xx_crust_mantle,R_yy_crust_mantle,R_xy_crust_mantle, &
               R_xz_crust_mantle,R_yz_crust_mantle, &
               epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
               epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
               eps_trace_over_3_crust_mantle, &
               alphaval,betaval,gammaval,factor_common_crust_mantle, &
               size(factor_common_crust_mantle,2), size(factor_common_crust_mantle,3), &
               size(factor_common_crust_mantle,4), size(factor_common_crust_mantle,5) )
          ! inner core region
          call compute_forces_inner_core( NSPEC_INNER_CORE_STR_OR_ATT,NGLOB_INNER_CORE, &
               NSPEC_INNER_CORE_ATTENUATION, &
               deltat, &
               displ_inner_core,veloc_inner_core,accel_inner_core, &
               phase_is_inner, &
               R_xx_inner_core,R_yy_inner_core,R_xy_inner_core,R_xz_inner_core,R_yz_inner_core, &
               epsilondev_xx_inner_core,epsilondev_yy_inner_core,epsilondev_xy_inner_core, &
               epsilondev_xz_inner_core,epsilondev_yz_inner_core, &
               eps_trace_over_3_inner_core,&
               alphaval,betaval,gammaval, &
               factor_common_inner_core, &
               size(factor_common_inner_core,2), size(factor_common_inner_core,3), &
               size(factor_common_inner_core,4), size(factor_common_inner_core,5) )

       endif

       ! adjoint / kernel runs
       if (SIMULATION_TYPE == 3 ) then
          if( USE_DEVILLE_PRODUCTS_VAL ) then
             ! uses Deville (2002) optimizations
             ! crust/mantle region
             call compute_forces_crust_mantle_Dev( NSPEC_CRUST_MANTLE_ADJOINT,NGLOB_CRUST_MANTLE_ADJOINT, &
                  NSPEC_CRUST_MANTLE_STR_AND_ATT, &
                  b_deltat, &
                  b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle, &
                  phase_is_inner, &
                  b_R_xx_crust_mantle,b_R_yy_crust_mantle,b_R_xy_crust_mantle, &
                  b_R_xz_crust_mantle,b_R_yz_crust_mantle, &
                  b_epsilondev_xx_crust_mantle,b_epsilondev_yy_crust_mantle,&
                  b_epsilondev_xy_crust_mantle, &
                  b_epsilondev_xz_crust_mantle,b_epsilondev_yz_crust_mantle, &
                  b_eps_trace_over_3_crust_mantle, &
                  b_alphaval,b_betaval,b_gammaval,factor_common_crust_mantle, &
                  size(factor_common_crust_mantle,2), size(factor_common_crust_mantle,3), &
                  size(factor_common_crust_mantle,4), size(factor_common_crust_mantle,5), .true. )
             ! inner core region
             call compute_forces_inner_core_Dev( NSPEC_INNER_CORE_ADJOINT,NGLOB_INNER_CORE_ADJOINT, &
                  NSPEC_INNER_CORE_STR_AND_ATT, &
                  b_deltat, &
                  b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core, &
                  phase_is_inner, &
                  b_R_xx_inner_core,b_R_yy_inner_core,b_R_xy_inner_core, &
                  b_R_xz_inner_core,b_R_yz_inner_core, &
                  b_epsilondev_xx_inner_core,b_epsilondev_yy_inner_core,b_epsilondev_xy_inner_core, &
                  b_epsilondev_xz_inner_core,b_epsilondev_yz_inner_core, &
                  b_eps_trace_over_3_inner_core,&
                  b_alphaval,b_betaval,b_gammaval, &
                  factor_common_inner_core, &
                  size(factor_common_inner_core,2), size(factor_common_inner_core,3), &
                  size(factor_common_inner_core,4), size(factor_common_inner_core,5), .true. )

          else
             ! no Deville optimization
             ! crust/mantle region
             call compute_forces_crust_mantle( NSPEC_CRUST_MANTLE_ADJOINT,NGLOB_CRUST_MANTLE_ADJOINT, &
                  NSPEC_CRUST_MANTLE_STR_AND_ATT, &
                  b_deltat, &
                  b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle, &
                  phase_is_inner, &
                  b_R_xx_crust_mantle,b_R_yy_crust_mantle,b_R_xy_crust_mantle, &
                  b_R_xz_crust_mantle,b_R_yz_crust_mantle, &
                  b_epsilondev_xx_crust_mantle,b_epsilondev_yy_crust_mantle,&
                  b_epsilondev_xy_crust_mantle, &
                  b_epsilondev_xz_crust_mantle,b_epsilondev_yz_crust_mantle, &
                  b_eps_trace_over_3_crust_mantle, &
                  b_alphaval,b_betaval,b_gammaval,factor_common_crust_mantle, &
                  size(factor_common_crust_mantle,2), size(factor_common_crust_mantle,3), &
                  size(factor_common_crust_mantle,4), size(factor_common_crust_mantle,5) )

             ! inner core region
             call compute_forces_inner_core( NSPEC_INNER_CORE_ADJOINT,NGLOB_INNER_CORE_ADJOINT, &
                  NSPEC_INNER_CORE_STR_AND_ATT, &
                  b_deltat, &
                  b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core, &
                  phase_is_inner, &
                  b_R_xx_inner_core,b_R_yy_inner_core,b_R_xy_inner_core, &
                  b_R_xz_inner_core,b_R_yz_inner_core, &
                  b_epsilondev_xx_inner_core,b_epsilondev_yy_inner_core,b_epsilondev_xy_inner_core, &
                  b_epsilondev_xz_inner_core,b_epsilondev_yz_inner_core, &
                  b_eps_trace_over_3_inner_core,&
                  b_alphaval,b_betaval,b_gammaval, &
                  factor_common_inner_core, &
                  size(factor_common_inner_core,2), size(factor_common_inner_core,3), &
                  size(factor_common_inner_core,4), size(factor_common_inner_core,5) )
          endif
       endif !SIMULATION_TYPE == 3

    else
       ! on GPU
       ! contains both forward SIM_TYPE==1 and backward SIM_TYPE==3 simulations
       ! for crust/mantle
       call compute_forces_crust_mantle_cuda(Mesh_pointer,iphase)
       ! for inner core
       call compute_forces_inner_core_cuda(Mesh_pointer,iphase)
    endif ! GPU_MODE


    ! computes additional contributions to acceleration field
    if( iphase == 1 ) then

       ! absorbing boundaries
       ! Stacey
       if(NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) call compute_stacey_crust_mantle()

       ! add the sources

       ! add adjoint sources
       if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
          if( nadj_rec_local > 0 ) call compute_add_sources_adjoint()
       endif

       ! add the sources
       select case( NOISE_TOMOGRAPHY )
       case( 0 )
          ! regular forward or backward simulation, no noise tomography simulation
          ! adds sources for forward simulation
          if (SIMULATION_TYPE == 1 .and. nsources_local > 0) &
            call compute_add_sources()
          ! add sources for backward/reconstructed wavefield
          if (SIMULATION_TYPE == 3 .and. nsources_local > 0) &
            call compute_add_sources_backward()

       case( 1 )
          ! the first step of noise tomography is to use |S(\omega)|^2 as a point force source at one of the receivers.
          ! hence, instead of a moment tensor 'sourcearrays', a 'noise_sourcearray' for a point force is needed.
          ! furthermore, the CMTSOLUTION needs to be zero, i.e., no earthquakes.
          ! now this must be manually set in DATA/CMTSOLUTION, by USERS.
          call noise_add_source_master_rec()

       case( 2 )
          ! second step of noise tomography, i.e., read the surface movie saved at every timestep
          ! use the movie to drive the ensemble forward wavefield
          call noise_read_add_surface_movie(NGLOB_CRUST_MANTLE,accel_crust_mantle,NSTEP-it+1)
          ! be careful, since ensemble forward sources are reversals of generating wavefield "eta"
          ! hence the "NSTEP-it+1", i.e., start reading from the last timestep
          ! note the ensemble forward sources are generally distributed on the surface of the earth
          ! that's to say, the ensemble forward source is kind of a surface force density, not a body force density
          ! therefore, we must add it here, before applying the inverse of mass matrix

       case( 3 )
          ! third step of noise tomography, i.e., read the surface movie saved at every timestep
          ! use the movie to reconstruct the ensemble forward wavefield
          ! the ensemble adjoint wavefield is done as usual
          ! note instead of "NSTEP-it+1", now we us "it", since reconstruction is a reversal of reversal
          call noise_read_add_surface_movie(NGLOB_CRUST_MANTLE_ADJOINT,b_accel_crust_mantle,it)

       end select


       ! ****************************************************
       ! **********  add matching with fluid part  **********
       ! ****************************************************
       ! only for elements in first matching layer in the solid
       if( .not. GPU_MODE ) then
          ! on CPU
          !---
          !--- couple with outer core at the bottom of the mantle
          !---
          if(ACTUALLY_COUPLE_FLUID_CMB) &
            call compute_coupling_CMB_fluid(displ_crust_mantle,b_displ_crust_mantle, &
                                           accel_crust_mantle,b_accel_crust_mantle, &
                                           ibool_crust_mantle,ibelm_bottom_crust_mantle,  &
                                           accel_outer_core,b_accel_outer_core, &
                                           normal_top_outer_core,jacobian2D_top_outer_core, &
                                           wgllwgll_xy,ibool_outer_core,ibelm_top_outer_core, &
                                           RHO_TOP_OC,minus_g_cmb, &
                                           SIMULATION_TYPE,NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE))

          !---
          !--- couple with outer core at the top of the inner core
          !---
          if(ACTUALLY_COUPLE_FLUID_ICB) &
            call compute_coupling_ICB_fluid(displ_inner_core,b_displ_inner_core, &
                                           accel_inner_core,b_accel_inner_core, &
                                           ibool_inner_core,ibelm_top_inner_core,  &
                                           accel_outer_core,b_accel_outer_core, &
                                           normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
                                           wgllwgll_xy,ibool_outer_core,ibelm_bottom_outer_core, &
                                           RHO_BOTTOM_OC,minus_g_icb, &
                                           SIMULATION_TYPE,NSPEC2D_TOP(IREGION_INNER_CORE))

       else
          ! on GPU
          !---
          !--- couple with outer core at the bottom of the mantle
          !---
          if( ACTUALLY_COUPLE_FLUID_CMB ) &
               call compute_coupling_cmb_fluid_cuda(Mesh_pointer)
          !---
          !--- couple with outer core at the top of the inner core
          !---
          if( ACTUALLY_COUPLE_FLUID_ICB ) &
               call compute_coupling_icb_fluid_cuda(Mesh_pointer)

       endif
    endif ! iphase == 1

    ! assemble all the contributions between slices using MPI

    ! crust/mantle and inner core handled in the same call
    ! in order to reduce the number of MPI messages by 2

    if( iphase == 1 ) then
      ! sends out MPI interface data
      if(.NOT. GPU_MODE) then
        ! on CPU
        ! sends accel values to corresponding MPI interface neighbors
        ! crust mantle
        call assemble_MPI_vector_s(NPROCTOT_VAL,NGLOB_CRUST_MANTLE, &
                      accel_crust_mantle, &
                      buffer_send_vector_crust_mantle,buffer_recv_vector_crust_mantle, &
                      num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                      nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle,&
                      my_neighbours_crust_mantle, &
                      request_send_vector_cm,request_recv_vector_cm)
        ! inner core
        call assemble_MPI_vector_s(NPROCTOT_VAL,NGLOB_INNER_CORE, &
                      accel_inner_core, &
                      buffer_send_vector_inner_core,buffer_recv_vector_inner_core, &
                      num_interfaces_inner_core,max_nibool_interfaces_ic, &
                      nibool_interfaces_inner_core,ibool_interfaces_inner_core,&
                      my_neighbours_inner_core, &
                      request_send_vector_ic,request_recv_vector_ic)
      else
        ! on GPU
        ! crust mantle
        call assemble_MPI_vector_send_cuda(NPROCTOT_VAL, &
                      buffer_send_vector_crust_mantle,buffer_recv_vector_crust_mantle, &
                      num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                      nibool_interfaces_crust_mantle,&
                      my_neighbours_crust_mantle, &
                      request_send_vector_cm,request_recv_vector_cm, &
                      IREGION_CRUST_MANTLE, &
                      1) ! <-- 1 == fwd accel
        ! inner core
        call assemble_MPI_vector_send_cuda(NPROCTOT_VAL, &
                      buffer_send_vector_inner_core,buffer_recv_vector_inner_core, &
                      num_interfaces_inner_core,max_nibool_interfaces_ic, &
                      nibool_interfaces_inner_core,&
                      my_neighbours_inner_core, &
                      request_send_vector_ic,request_recv_vector_ic, &
                      IREGION_INNER_CORE, &
                      1)
      endif ! GPU_MODE

      ! adjoint / kernel runs
      if (SIMULATION_TYPE == 3) then
        if(.NOT. GPU_MODE) then
          ! on CPU
          ! sends accel values to corresponding MPI interface neighbors
          ! crust mantle
          call assemble_MPI_vector_s(NPROCTOT_VAL,NGLOB_CRUST_MANTLE, &
                        b_accel_crust_mantle, &
                        b_buffer_send_vector_cm,b_buffer_recv_vector_cm, &
                        num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                        nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle,&
                        my_neighbours_crust_mantle, &
                        b_request_send_vector_cm,b_request_recv_vector_cm)
          ! inner core
          call assemble_MPI_vector_s(NPROCTOT_VAL,NGLOB_INNER_CORE, &
                        b_accel_inner_core, &
                        b_buffer_send_vector_inner_core,b_buffer_recv_vector_inner_core, &
                        num_interfaces_inner_core,max_nibool_interfaces_ic, &
                        nibool_interfaces_inner_core,ibool_interfaces_inner_core,&
                        my_neighbours_inner_core, &
                        b_request_send_vector_ic,b_request_recv_vector_ic)
        else
          ! on GPU
          ! crust mantle
          call assemble_MPI_vector_send_cuda(NPROCTOT_VAL, &
                      b_buffer_send_vector_cm,b_buffer_recv_vector_cm, &
                      num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                      nibool_interfaces_crust_mantle,&
                      my_neighbours_crust_mantle, &
                      b_request_send_vector_cm,b_request_recv_vector_cm, &
                      IREGION_CRUST_MANTLE, &
                      3) ! <-- 3 == adjoint b_accel
          ! inner core
          call assemble_MPI_vector_send_cuda(NPROCTOT_VAL, &
                      b_buffer_send_vector_inner_core,b_buffer_recv_vector_inner_core, &
                      num_interfaces_inner_core,max_nibool_interfaces_ic, &
                      nibool_interfaces_inner_core,&
                      my_neighbours_inner_core, &
                      b_request_send_vector_ic,b_request_recv_vector_ic, &
                      IREGION_INNER_CORE, &
                      3)
        endif ! GPU
      endif ! SIMULATION_TYPE == 3

    else
      ! waits for send/receive requests to be completed and assembles values
      if(.NOT. GPU_MODE) then
        ! on CPU
        ! crust mantle
        call assemble_MPI_vector_w(NPROCTOT_VAL,NGLOB_CRUST_MANTLE, &
                              accel_crust_mantle, &
                              buffer_recv_vector_crust_mantle,num_interfaces_crust_mantle,&
                              max_nibool_interfaces_cm, &
                              nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
                              request_send_vector_cm,request_recv_vector_cm)
        ! inner core
        call assemble_MPI_vector_w(NPROCTOT_VAL,NGLOB_INNER_CORE, &
                              accel_inner_core, &
                              buffer_recv_vector_inner_core,num_interfaces_inner_core,&
                              max_nibool_interfaces_ic, &
                              nibool_interfaces_inner_core,ibool_interfaces_inner_core, &
                              request_send_vector_ic,request_recv_vector_ic)
      else
        ! on GPU
        ! crust mantle
        call assemble_MPI_vector_write_cuda(NPROCTOT_VAL, &
                            buffer_recv_vector_crust_mantle, &
                            num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                            request_send_vector_cm,request_recv_vector_cm, &
                            IREGION_CRUST_MANTLE, &
                            1) ! <-- 1 == fwd accel
        ! inner core
        call assemble_MPI_vector_write_cuda(NPROCTOT_VAL, &
                            buffer_recv_vector_inner_core, &
                            num_interfaces_inner_core,max_nibool_interfaces_ic, &
                            request_send_vector_ic,request_recv_vector_ic, &
                            IREGION_INNER_CORE, &
                            1)
      endif


      ! adjoint / kernel runs
      if (SIMULATION_TYPE == 3) then
        ! waits for send/receive requests to be completed and assembles values
        if(.NOT. GPU_MODE) then
          ! on CPU
          ! crust mantle
          call assemble_MPI_vector_w(NPROCTOT_VAL,NGLOB_CRUST_MANTLE, &
                              b_accel_crust_mantle, &
                              b_buffer_recv_vector_cm,num_interfaces_crust_mantle,&
                              max_nibool_interfaces_cm, &
                              nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
                              b_request_send_vector_cm,b_request_recv_vector_cm)
          ! inner core
          call assemble_MPI_vector_w(NPROCTOT_VAL,NGLOB_INNER_CORE, &
                              b_accel_inner_core, &
                              b_buffer_recv_vector_inner_core,num_interfaces_inner_core,&
                              max_nibool_interfaces_ic, &
                              nibool_interfaces_inner_core,ibool_interfaces_inner_core, &
                              b_request_send_vector_ic,b_request_recv_vector_ic)

        else
          ! on GPU
          ! crust mantle
          call assemble_MPI_vector_write_cuda(NPROCTOT_VAL, &
                            b_buffer_recv_vector_cm, &
                            num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                            b_request_send_vector_cm,b_request_recv_vector_cm, &
                            IREGION_CRUST_MANTLE, &
                            3) ! <-- 3 == adjoint b_accel
          ! inner core
          call assemble_MPI_vector_write_cuda(NPROCTOT_VAL,&
                            b_buffer_recv_vector_inner_core, &
                            num_interfaces_inner_core,max_nibool_interfaces_ic, &
                            b_request_send_vector_ic,b_request_recv_vector_ic, &
                            IREGION_INNER_CORE, &
                            3)
        endif
      endif ! SIMULATION_TYPE == 3
    endif ! iphase == 1

  enddo ! iphase

  ! updates (only) acceleration w/ rotation in the crust/mantle region (touches oceans)
  if(.NOT. GPU_MODE) then
     ! on CPU
     call compute_forces_el_update_accel(NGLOB_CRUST_MANTLE,NGLOB_XY_CM,veloc_crust_mantle,accel_crust_mantle, &
                                        two_omega_earth, &
                                        rmassx_crust_mantle,rmassy_crust_mantle,rmassz_crust_mantle, &
                                        NCHUNKS_VAL,ABSORBING_CONDITIONS)
     ! adjoint / kernel runs
     if (SIMULATION_TYPE == 3) &
          call compute_forces_el_update_accel(NGLOB_CRUST_MANTLE_ADJOINT,NGLOB_XY_CM,b_veloc_crust_mantle,b_accel_crust_mantle, &
                                             b_two_omega_earth, &
                                             rmassx_crust_mantle,rmassy_crust_mantle,rmassz_crust_mantle, &
                                             NCHUNKS_VAL,ABSORBING_CONDITIONS)
  else
     ! on GPU
     call kernel_3_a_cuda(Mesh_pointer, &
                         deltatover2,SIMULATION_TYPE,b_deltatover2,NCHUNKS_VAL)
  endif

  ! couples ocean with crust mantle
  ! (updates acceleration with ocean load approximation)
  if( OCEANS_VAL ) then
    if(.NOT. GPU_MODE) then
      ! on CPU
      call compute_coupling_ocean(accel_crust_mantle,b_accel_crust_mantle, &
             rmassx_crust_mantle, rmassy_crust_mantle, rmassz_crust_mantle, &
             rmass_ocean_load,normal_top_crust_mantle, &
             ibool_crust_mantle,ibelm_top_crust_mantle, &
             updated_dof_ocean_load,NGLOB_XY_CM, &
             SIMULATION_TYPE,NSPEC2D_TOP(IREGION_CRUST_MANTLE), &
             ABSORBING_CONDITIONS)

    else
      ! on GPU
      call compute_coupling_ocean_cuda(Mesh_pointer,NCHUNKS_VAL)
    endif
  endif

  ! Newmark time scheme:
  ! corrector terms for elastic parts
  ! (updates velocity)
  if(.NOT. GPU_MODE ) then
    ! on CPU
    call compute_forces_el_update_veloc(NGLOB_CRUST_MANTLE,veloc_crust_mantle,accel_crust_mantle, &
                                       NGLOB_INNER_CORE,veloc_inner_core,accel_inner_core, &
                                       deltatover2,two_omega_earth,rmass_inner_core)
    ! adjoint / kernel runs
    if (SIMULATION_TYPE == 3) &
      call compute_forces_el_update_veloc(NGLOB_CRUST_MANTLE_ADJOINT,b_veloc_crust_mantle,b_accel_crust_mantle, &
                                         NGLOB_INNER_CORE_ADJOINT,b_veloc_inner_core,b_accel_inner_core, &
                                         b_deltatover2,b_two_omega_earth,rmass_inner_core)
  else
    ! on GPU
    call kernel_3_b_cuda(Mesh_pointer, &
                        deltatover2,SIMULATION_TYPE,b_deltatover2)
  endif


!daniel: att - debug
!  if( DEBUG ) then
!    if( SIMULATION_TYPE == 1) then
!      if( it > NSTEP - 1000 .and. myrank == 0 ) then
!        print*,'it',it,'Rxx:',R_xx_crust_mantle(1,1,1,1,100),epsilondev_xx_crust_mantle(1,1,1,100)
!      endif
!    else if( SIMULATION_TYPE == 3 ) then
!      if( it <= 1000 .and. myrank == 0 ) then
!        print*,'it',it,'Rxx:',b_R_xx_crust_mantle(1,1,1,1,100),b_epsilondev_xx_crust_mantle(1,1,1,100)
!      endif
!    endif
!  endif

  end subroutine compute_forces_elastic


!=====================================================================

  subroutine compute_forces_el_update_accel(NGLOB,NGLOB_XY,veloc_crust_mantle,accel_crust_mantle, &
                                           two_omega_earth, &
                                           rmassx_crust_mantle,rmassy_crust_mantle,rmassz_crust_mantle, &
                                           NCHUNKS_VAL,ABSORBING_CONDITIONS)

  use constants_solver,only: CUSTOM_REAL,NDIM

  implicit none

  integer :: NGLOB,NGLOB_XY,NCHUNKS_VAL

  ! velocity & acceleration
  ! crust/mantle region
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: veloc_crust_mantle,accel_crust_mantle

  ! mass matrices
  !
  ! in the case of stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
  ! on Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be obsolete
  real(kind=CUSTOM_REAL), dimension(NGLOB_XY) :: rmassx_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_XY) :: rmassy_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB)    :: rmassz_crust_mantle

  real(kind=CUSTOM_REAL) :: two_omega_earth

  logical :: ABSORBING_CONDITIONS

  ! local parameters
  integer :: i

  ! updates acceleration w/ rotation in crust/mantle region only

  if(NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) then

     do i=1,NGLOB
        accel_crust_mantle(1,i) = accel_crust_mantle(1,i)*rmassx_crust_mantle(i) &
             + two_omega_earth*veloc_crust_mantle(2,i)
        accel_crust_mantle(2,i) = accel_crust_mantle(2,i)*rmassy_crust_mantle(i) &
             - two_omega_earth*veloc_crust_mantle(1,i)
        accel_crust_mantle(3,i) = accel_crust_mantle(3,i)*rmassz_crust_mantle(i)
     enddo

  else

     do i=1,NGLOB
        accel_crust_mantle(1,i) = accel_crust_mantle(1,i)*rmassz_crust_mantle(i) &
             + two_omega_earth*veloc_crust_mantle(2,i)
        accel_crust_mantle(2,i) = accel_crust_mantle(2,i)*rmassz_crust_mantle(i) &
             - two_omega_earth*veloc_crust_mantle(1,i)
        accel_crust_mantle(3,i) = accel_crust_mantle(3,i)*rmassz_crust_mantle(i)
     enddo

  endif

  end subroutine compute_forces_el_update_accel


!=====================================================================

  subroutine compute_forces_el_update_veloc(NGLOB_CM,veloc_crust_mantle,accel_crust_mantle, &
                                            NGLOB_IC,veloc_inner_core,accel_inner_core, &
                                            deltatover2,two_omega_earth,rmass_inner_core)

  use constants_solver,only: CUSTOM_REAL,NDIM

  implicit none

  integer :: NGLOB_CM,NGLOB_IC

  ! acceleration & velocity
  ! crust/mantle region
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CM) :: veloc_crust_mantle,accel_crust_mantle
  ! inner core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_IC) :: veloc_inner_core,accel_inner_core

  ! mass matrix
  real(kind=CUSTOM_REAL), dimension(NGLOB_IC) :: rmass_inner_core

  real(kind=CUSTOM_REAL) :: deltatover2,two_omega_earth

  ! local parameters
  integer :: i

  ! Newmark time scheme:
  !
  ! note:
  !   - crust/mantle region
  !         needs only velocity corrector terms
  !         (acceleration already updated before)
  !   - inner core region
  !         needs both, acceleration update & velocity corrector terms

  ! mantle
  do i=1,NGLOB_CM
    veloc_crust_mantle(:,i) = veloc_crust_mantle(:,i) + deltatover2*accel_crust_mantle(:,i)
  enddo
  ! inner core
  do i=1,NGLOB_IC
    accel_inner_core(1,i) = accel_inner_core(1,i)*rmass_inner_core(i) &
           + two_omega_earth*veloc_inner_core(2,i)
    accel_inner_core(2,i) = accel_inner_core(2,i)*rmass_inner_core(i) &
           - two_omega_earth*veloc_inner_core(1,i)
    accel_inner_core(3,i) = accel_inner_core(3,i)*rmass_inner_core(i)

    veloc_inner_core(:,i) = veloc_inner_core(:,i) + deltatover2*accel_inner_core(:,i)
  enddo

  end subroutine compute_forces_el_update_veloc
