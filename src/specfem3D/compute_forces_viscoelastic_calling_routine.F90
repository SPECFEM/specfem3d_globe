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

  subroutine compute_forces_viscoelastic()

! elastic domains for forward or adjoint simulations (SIMULATION_TYPE == 1 or 2 )

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_movie

  implicit none

! note: instead of using size(factor_common_crust_mantle,5), we use the ATT4_VAL/ATT5_VAL in the function calls;
!       this is due to the size(..) function returning either integer(kind=4) or integer(kind=8)
!       depending on compiler flags (-mcmedium), leading to unpredictable results

  ! local parameters
  ! non blocking MPI
  ! iphase: iphase = 1 is for computing outer elements in the crust_mantle and inner_core regions,
  !              iphase = 2 is for computing inner elements (former icall parameter)
  integer :: iphase

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
    !
    ! compute internal forces in the solid regions
    ! note: for anisotropy and gravity, x y and z contain r theta and phi
    if (.not. GPU_MODE) then
      ! on CPU
      ! crust/mantle region
      call compute_forces_crust_mantle(NSPEC_CRUST_MANTLE_STR_OR_ATT,NGLOB_CRUST_MANTLE, &
                                       NSPEC_CRUST_MANTLE_ATTENUATION, &
                                       deltat, &
                                       displ_crust_mantle, &
                                       accel_crust_mantle, &
                                       iphase, &
                                       R_xx_crust_mantle,R_yy_crust_mantle,R_xy_crust_mantle, &
                                       R_xz_crust_mantle,R_yz_crust_mantle, &
                                       R_xx_crust_mantle_lddrk,R_yy_crust_mantle_lddrk,R_xy_crust_mantle_lddrk, &
                                       R_xz_crust_mantle_lddrk,R_yz_crust_mantle_lddrk, &
                                       epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                                       epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                                       eps_trace_over_3_crust_mantle, &
                                       alphaval,betaval,gammaval)

      ! inner core region
      call compute_forces_inner_core(NSPEC_INNER_CORE_STR_OR_ATT,NGLOB_INNER_CORE, &
                                     NSPEC_INNER_CORE_ATTENUATION, &
                                     deltat, &
                                     displ_inner_core, &
                                     accel_inner_core, &
                                     iphase, &
                                     R_xx_inner_core,R_yy_inner_core,R_xy_inner_core, &
                                     R_xz_inner_core,R_yz_inner_core, &
                                     R_xx_inner_core_lddrk,R_yy_inner_core_lddrk,R_xy_inner_core_lddrk, &
                                     R_xz_inner_core_lddrk,R_yz_inner_core_lddrk, &
                                     epsilondev_xx_inner_core,epsilondev_yy_inner_core,epsilondev_xy_inner_core, &
                                     epsilondev_xz_inner_core,epsilondev_yz_inner_core, &
                                     eps_trace_over_3_inner_core, &
                                     alphaval,betaval,gammaval)
    else
      ! on GPU
      ! contains forward FORWARD_OR_ADJOINT == 1
      ! for crust/mantle
      call compute_forces_crust_mantle_gpu(Mesh_pointer,iphase,ALPHA_LDDRK(istage),BETA_LDDRK(istage),1)

      ! initiates asynchronous MPI transfer
      if (NPROCTOT_VAL > 1) then
        if (GPU_ASYNC_COPY .and. iphase == 2) then
          ! crust/mantle region
          ! wait for asynchronous copy to finish
          call sync_copy_from_device(Mesh_pointer,iphase,buffer_send_vector_crust_mantle,IREGION_CRUST_MANTLE,1)
          ! sends MPI buffers
          call assemble_MPI_vector_send_gpu(NPROCTOT_VAL, &
                        buffer_send_vector_crust_mantle,buffer_recv_vector_crust_mantle, &
                        num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                        nibool_interfaces_crust_mantle, &
                        my_neighbors_crust_mantle, &
                        request_send_vector_cm,request_recv_vector_cm)
        endif
      endif

      ! for inner core
      call compute_forces_inner_core_gpu(Mesh_pointer,iphase,ALPHA_LDDRK(istage),BETA_LDDRK(istage),1)

      ! initiates asynchronous MPI transfer
      if (NPROCTOT_VAL > 1) then
        if (GPU_ASYNC_COPY .and. iphase == 2) then
          ! inner core region
          ! wait for asynchronous copy to finish
          call sync_copy_from_device(Mesh_pointer,iphase,buffer_send_vector_inner_core,IREGION_INNER_CORE,1)
          ! sends MPI buffers
          call assemble_MPI_vector_send_gpu(NPROCTOT_VAL, &
                        buffer_send_vector_inner_core,buffer_recv_vector_inner_core, &
                        num_interfaces_inner_core,max_nibool_interfaces_ic, &
                        nibool_interfaces_inner_core, &
                        my_neighbors_inner_core, &
                        request_send_vector_ic,request_recv_vector_ic)
        endif
      endif
    endif ! GPU_MODE

    ! computes additional contributions to acceleration field
    if (iphase == 1) then
       ! during phase for outer elements

       ! absorbing boundaries
       ! Stacey
       if (NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) call compute_stacey_cm_forward()

       ! add the sources

       ! add adjoint sources
       ! note: this must remain here even when SIMULATION_TYPE == 3 because it applies to array
       !       accel_crust_mantle rather than b_accel_crust_mantle
       if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
          if (nadj_rec_local > 0) call compute_add_sources_adjoint()
       endif

       ! add the sources
       select case (NOISE_TOMOGRAPHY)
       case (0)
          ! regular forward or backward simulation, no noise tomography simulation
          ! adds sources for forward simulation
          if (SIMULATION_TYPE == 1 .and. nsources_local > 0) &
            call compute_add_sources()

       case (1)
          ! the first step of noise tomography is to use |S(\omega)|^2 as a point force source at one of the receivers.
          ! hence, instead of a moment tensor 'sourcearrays', a 'noise_sourcearray' for a point force is needed.
          ! furthermore, the CMTSOLUTION needs to be zero, i.e., no earthquakes.
          call noise_add_source_main_rec()

       case (2)
          ! second step of noise tomography, i.e., read the surface movie saved at every timestep
          ! use the movie to drive the ensemble forward wavefield
          call noise_read_add_surface_movie(NGLOB_CRUST_MANTLE,accel_crust_mantle,NSTEP-it+1)
          ! be careful, since ensemble forward sources are reversals of generating wavefield "eta"
          ! hence the "NSTEP-it+1", i.e., start reading from the last timestep
          ! note the ensemble forward sources are generally distributed on the surface of the earth
          ! that's to say, the ensemble forward source is kind of a surface force density, not a body force density
          ! therefore, we must add it here, before applying the inverse of mass matrix
       end select

       ! ****************************************************
       ! **********  add matching with fluid part  **********
       ! ****************************************************
       call compute_coupling_solid()

    endif ! iphase == 1

    ! assemble all the contributions between slices using MPI

    ! crust/mantle and inner core handled in the same call
    ! in order to reduce the number of MPI messages by 2

    if (iphase == 1) then
      ! sends out MPI interface data
      if (NPROCTOT_VAL > 1) then
        if (.not. GPU_MODE) then
          ! on CPU
          ! sends accel values to corresponding MPI interface neighbors
          ! crust mantle
          call assemble_MPI_vector_s(NPROCTOT_VAL,NGLOB_CRUST_MANTLE, &
                        accel_crust_mantle, &
                        buffer_send_vector_crust_mantle,buffer_recv_vector_crust_mantle, &
                        num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                        nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
                        my_neighbors_crust_mantle, &
                        request_send_vector_cm,request_recv_vector_cm)
          ! inner core
          call assemble_MPI_vector_s(NPROCTOT_VAL,NGLOB_INNER_CORE, &
                        accel_inner_core, &
                        buffer_send_vector_inner_core,buffer_recv_vector_inner_core, &
                        num_interfaces_inner_core,max_nibool_interfaces_ic, &
                        nibool_interfaces_inner_core,ibool_interfaces_inner_core, &
                        my_neighbors_inner_core, &
                        request_send_vector_ic,request_recv_vector_ic)
        else
          ! on GPU
          ! sends accel values to corresponding MPI interface neighbors

          ! preparation of the contribution between partitions using MPI
          ! transfers MPI buffers to CPU
          ! note: in case of asynchronous copy, this transfers boundary region to host asynchronously. The
          !       MPI-send is done after compute_forces_viscoelastic_gpu,
          !       once the inner element kernels are launched, and the memcpy has finished.
          call transfer_boun_from_device(Mesh_pointer, &
                                         buffer_send_vector_crust_mantle, &
                                         IREGION_CRUST_MANTLE,1)
          call transfer_boun_from_device(Mesh_pointer, &
                                         buffer_send_vector_inner_core, &
                                         IREGION_INNER_CORE,1)

          if (.not. GPU_ASYNC_COPY) then
            ! for synchronous transfers, sending over MPI can directly proceed
            ! crust mantle
            call assemble_MPI_vector_send_gpu(NPROCTOT_VAL, &
                          buffer_send_vector_crust_mantle,buffer_recv_vector_crust_mantle, &
                          num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                          nibool_interfaces_crust_mantle, &
                          my_neighbors_crust_mantle, &
                          request_send_vector_cm,request_recv_vector_cm)
            ! inner core
            call assemble_MPI_vector_send_gpu(NPROCTOT_VAL, &
                          buffer_send_vector_inner_core,buffer_recv_vector_inner_core, &
                          num_interfaces_inner_core,max_nibool_interfaces_ic, &
                          nibool_interfaces_inner_core, &
                          my_neighbors_inner_core, &
                          request_send_vector_ic,request_recv_vector_ic)
          endif
        endif ! GPU_MODE
      endif
    else
      ! waits for send/receive requests to be completed and assembles values
      if (NPROCTOT_VAL > 1) then
        if (.not. GPU_MODE) then
          ! on CPU
          ! crust mantle
          call assemble_MPI_vector_w(NPROCTOT_VAL,NGLOB_CRUST_MANTLE, &
                                accel_crust_mantle, &
                                buffer_recv_vector_crust_mantle,num_interfaces_crust_mantle, &
                                max_nibool_interfaces_cm, &
                                nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
                                my_neighbors_crust_mantle, &
                                request_send_vector_cm,request_recv_vector_cm)
          ! inner core
          call assemble_MPI_vector_w(NPROCTOT_VAL,NGLOB_INNER_CORE, &
                                accel_inner_core, &
                                buffer_recv_vector_inner_core,num_interfaces_inner_core, &
                                max_nibool_interfaces_ic, &
                                nibool_interfaces_inner_core,ibool_interfaces_inner_core, &
                                my_neighbors_inner_core, &
                                request_send_vector_ic,request_recv_vector_ic)
        else
          ! on GPU
          if (GPU_ASYNC_COPY) then
            ! while inner elements compute "Kernel_2", we wait for MPI to
            ! finish and transfer the boundary terms to the device asynchronously
            !
            ! transfers MPI buffers onto GPU
            ! crust/mantle region
            call transfer_boundary_to_device(Mesh_pointer,NPROCTOT_VAL,buffer_recv_vector_crust_mantle, &
                                             num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                                             request_recv_vector_cm, &
                                             IREGION_CRUST_MANTLE,1)
            ! inner core region
            call transfer_boundary_to_device(Mesh_pointer,NPROCTOT_VAL,buffer_recv_vector_inner_core, &
                                             num_interfaces_inner_core,max_nibool_interfaces_ic, &
                                             request_recv_vector_ic, &
                                             IREGION_INNER_CORE,1)
          endif

          ! waits for MPI send/receive requests to be completed and assembles values
          ! crust mantle
          call assemble_MPI_vector_write_gpu(Mesh_pointer,NPROCTOT_VAL, &
                              buffer_recv_vector_crust_mantle, &
                              num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                              request_send_vector_cm,request_recv_vector_cm, &
                              IREGION_CRUST_MANTLE,1) ! -- 1 == fwd accel
          ! inner core
          call assemble_MPI_vector_write_gpu(Mesh_pointer,NPROCTOT_VAL, &
                              buffer_recv_vector_inner_core, &
                              num_interfaces_inner_core,max_nibool_interfaces_ic, &
                              request_send_vector_ic,request_recv_vector_ic, &
                              IREGION_INNER_CORE,1)
        endif
      endif
    endif ! iphase == 1

  enddo ! iphase

  ! updates (only) acceleration w/ rotation in the crust/mantle and inner core region
  if (.not. GPU_MODE) then
    ! on CPU
    ! crust/mantle and inner core
    call multiply_accel_elastic(two_omega_earth, &
                                NGLOB_CRUST_MANTLE,veloc_crust_mantle,accel_crust_mantle, &
                                rmassx_crust_mantle,rmassy_crust_mantle,rmassz_crust_mantle, &
                                NGLOB_INNER_CORE,veloc_inner_core,accel_inner_core, &
                                rmassx_inner_core,rmassy_inner_core,rmassz_inner_core)
  else
    ! on GPU
    ! includes FORWARD_OR_ADJOINT == 1
    call multiply_accel_elastic_gpu(Mesh_pointer,1)
  endif

  ! couples ocean with crust mantle
  ! (updates acceleration with ocean load approximation)
  if (OCEANS_VAL) then
    if (.not. GPU_MODE) then
      ! on CPU
      call compute_coupling_ocean(NGLOB_CRUST_MANTLE,accel_crust_mantle, &
                                  rmassx_crust_mantle,rmassy_crust_mantle,rmassz_crust_mantle, &
                                  npoin_oceans,ibool_ocean_load,rmass_ocean_load_selected,normal_ocean_load)
    else
      ! on GPU
      call compute_coupling_ocean_gpu(Mesh_pointer,1) ! -- 1 == forward arrays
    endif
  endif

  ! time scheme update
  if (USE_LDDRK) then
    ! Runge-Kutta scheme
    call update_veloc_elastic_lddrk()
  else
    ! Newmark time scheme
    call update_veloc_elastic_newmark()
  endif

  end subroutine compute_forces_viscoelastic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_forces_viscoelastic_backward()

! backward/reconstructed wavefields only

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_movie

  implicit none

  ! local parameters
  ! non blocking MPI
  ! iphase: iphase = 1 is for computing outer elements in the crust_mantle and inner_core regions,
  !              iphase = 2 is for computing inner elements (former icall parameter)
  integer :: iphase

  ! checks
  if (SIMULATION_TYPE /= 3 ) return

!daniel debug: att - debug
!  integer :: iglob
!  logical,parameter :: DEBUG = .false.
!  if (DEBUG) then
!    iglob = ibool_crust_mantle(1,1,1,100)
!    if (SIMULATION_TYPE == 1) then
!      if (it == NSTEP .and. myrank == 0) then
!        print *,'last step',it,'Rxx:',R_xx_crust_mantle(1,1,1,1,100),epsilondev_xx_crust_mantle(1,1,1,100), &
!          displ_crust_mantle(1,iglob),displ_crust_mantle(2,iglob),displ_crust_mantle(3,iglob)
!      endif
!      if (it == NSTEP-1 .and. myrank == 0) then
!        print *,'second last step',it,'Rxx:',R_xx_crust_mantle(1,1,1,1,100),epsilondev_xx_crust_mantle(1,1,1,100), &
!          displ_crust_mantle(1,iglob),displ_crust_mantle(2,iglob),displ_crust_mantle(3,iglob)
!      endif
!      if (it == NSTEP-2 .and. myrank == 0) then
!        print *,'third last step',it,'Rxx:',R_xx_crust_mantle(1,1,1,1,100),epsilondev_xx_crust_mantle(1,1,1,100), &
!          displ_crust_mantle(1,iglob),displ_crust_mantle(2,iglob),displ_crust_mantle(3,iglob)
!      endif
!    else if (SIMULATION_TYPE == 3) then
!      if (it == 1 .and. myrank == 0) then
!        print *,'first step',it,'Rxx:',b_R_xx_crust_mantle(1,1,1,1,100),b_epsilondev_xx_crust_mantle(1,1,1,100), &
!          b_displ_crust_mantle(1,iglob),b_displ_crust_mantle(2,iglob),b_displ_crust_mantle(3,iglob)
!      endif
!      if (it == 2 .and. myrank == 0) then
!        print *,'second step',it,'Rxx:',b_R_xx_crust_mantle(1,1,1,1,100),b_epsilondev_xx_crust_mantle(1,1,1,100), &
!          b_displ_crust_mantle(1,iglob),b_displ_crust_mantle(2,iglob),b_displ_crust_mantle(3,iglob)
!      endif
!      if (it == 3 .and. myrank == 0) then
!        print *,'third step',it,'Rxx:',b_R_xx_crust_mantle(1,1,1,1,100),b_epsilondev_xx_crust_mantle(1,1,1,100), &
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
    !
    ! compute internal forces in the solid regions
    ! note: for anisotropy and gravity, x y and z contain r theta and phi
    if (.not. GPU_MODE) then
      ! on CPU
      ! adjoint / kernel runs
      ! crust/mantle region
      call compute_forces_crust_mantle(NSPEC_CRUST_MANTLE_ADJOINT,NGLOB_CRUST_MANTLE_ADJOINT, &
                                       NSPEC_CRUST_MANTLE_STR_AND_ATT, &
                                       b_deltat, &
                                       b_displ_crust_mantle, &
                                       b_accel_crust_mantle, &
                                       iphase, &
                                       b_R_xx_crust_mantle,b_R_yy_crust_mantle,b_R_xy_crust_mantle, &
                                       b_R_xz_crust_mantle,b_R_yz_crust_mantle, &
                                       b_R_xx_crust_mantle_lddrk,b_R_yy_crust_mantle_lddrk,b_R_xy_crust_mantle_lddrk, &
                                       b_R_xz_crust_mantle_lddrk,b_R_yz_crust_mantle_lddrk, &
                                       b_epsilondev_xx_crust_mantle,b_epsilondev_yy_crust_mantle, &
                                       b_epsilondev_xy_crust_mantle, &
                                       b_epsilondev_xz_crust_mantle,b_epsilondev_yz_crust_mantle, &
                                       b_eps_trace_over_3_crust_mantle, &
                                       b_alphaval,b_betaval,b_gammaval)
      ! inner core region
      call compute_forces_inner_core(NSPEC_INNER_CORE_ADJOINT,NGLOB_INNER_CORE_ADJOINT, &
                                     NSPEC_INNER_CORE_STR_AND_ATT, &
                                     b_deltat, &
                                     b_displ_inner_core, &
                                     b_accel_inner_core, &
                                     iphase, &
                                     b_R_xx_inner_core,b_R_yy_inner_core,b_R_xy_inner_core, &
                                     b_R_xz_inner_core,b_R_yz_inner_core, &
                                     b_R_xx_inner_core_lddrk,b_R_yy_inner_core_lddrk,b_R_xy_inner_core_lddrk, &
                                     b_R_xz_inner_core_lddrk,b_R_yz_inner_core_lddrk, &
                                     b_epsilondev_xx_inner_core,b_epsilondev_yy_inner_core,b_epsilondev_xy_inner_core, &
                                     b_epsilondev_xz_inner_core,b_epsilondev_yz_inner_core, &
                                     b_eps_trace_over_3_inner_core, &
                                     b_alphaval,b_betaval,b_gammaval)
    else
      ! on GPU
      ! contains forward FORWARD_OR_ADJOINT == 3
      ! for crust/mantle
      call compute_forces_crust_mantle_gpu(Mesh_pointer,iphase,ALPHA_LDDRK(istage),BETA_LDDRK(istage),3)

      ! initiates asynchronous MPI transfer
      if (GPU_ASYNC_COPY .and. iphase == 2) then
        ! crust/mantle region
        ! wait for asynchronous copy to finish
        call sync_copy_from_device(Mesh_pointer,iphase,b_buffer_send_vector_cm,IREGION_CRUST_MANTLE,3)
        ! sends MPI buffers
        call assemble_MPI_vector_send_gpu(NPROCTOT_VAL, &
                      b_buffer_send_vector_cm,b_buffer_recv_vector_cm, &
                      num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                      nibool_interfaces_crust_mantle, &
                      my_neighbors_crust_mantle, &
                      b_request_send_vector_cm,b_request_recv_vector_cm)
      endif

      ! for inner core
      call compute_forces_inner_core_gpu(Mesh_pointer,iphase,ALPHA_LDDRK(istage),BETA_LDDRK(istage),3)

      ! initiates asynchronous MPI transfer
      if (GPU_ASYNC_COPY .and. iphase == 2) then
        ! inner core region
        ! wait for asynchronous copy to finish
        call sync_copy_from_device(Mesh_pointer,iphase,b_buffer_send_vector_inner_core,IREGION_INNER_CORE,3)

        ! sends MPI buffers
        call assemble_MPI_vector_send_gpu(NPROCTOT_VAL, &
                      b_buffer_send_vector_inner_core,b_buffer_recv_vector_inner_core, &
                      num_interfaces_inner_core,max_nibool_interfaces_ic, &
                      nibool_interfaces_inner_core, &
                      my_neighbors_inner_core, &
                      b_request_send_vector_ic,b_request_recv_vector_ic)
      endif

    endif ! GPU_MODE

    ! computes additional contributions to acceleration field
    if (iphase == 1) then

      ! absorbing boundaries
      ! Stacey
      if (NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) then
        if (UNDO_ATTENUATION) then
          call compute_stacey_cm_backward_undoatt()
        else
          call compute_stacey_cm_backward()
        endif
      endif

      ! add the sources
      select case (NOISE_TOMOGRAPHY)
      case (0)
        ! add sources for backward/reconstructed wavefield
        if (nsources_local > 0 ) &
          call compute_add_sources_backward()

      case (3)
        ! third step of noise tomography, i.e., read the surface movie saved at every timestep
        ! use the movie to reconstruct the ensemble forward wavefield
        ! the ensemble adjoint wavefield is done as usual
        !
        ! note instead of "NSTEP-it+1", now we use "it", since reconstruction is a reversal of reversal
        call noise_read_add_surface_movie(NGLOB_CRUST_MANTLE_ADJOINT,b_accel_crust_mantle,it)

      end select

      ! ****************************************************
      ! **********  add matching with fluid part  **********
      ! ****************************************************
      call compute_coupling_solid_backward()

    endif ! iphase == 1

    ! assemble all the contributions between slices using MPI

    ! crust/mantle and inner core handled in the same call
    ! in order to reduce the number of MPI messages by 2

    if (iphase == 1) then
      ! sends out MPI interface data
      ! adjoint / kernel runs
      if (.not. GPU_MODE) then
        ! on CPU
        ! sends accel values to corresponding MPI interface neighbors
        ! crust mantle
        call assemble_MPI_vector_s(NPROCTOT_VAL,NGLOB_CRUST_MANTLE, &
                      b_accel_crust_mantle, &
                      b_buffer_send_vector_cm,b_buffer_recv_vector_cm, &
                      num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                      nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
                      my_neighbors_crust_mantle, &
                      b_request_send_vector_cm,b_request_recv_vector_cm)
        ! inner core
        call assemble_MPI_vector_s(NPROCTOT_VAL,NGLOB_INNER_CORE, &
                      b_accel_inner_core, &
                      b_buffer_send_vector_inner_core,b_buffer_recv_vector_inner_core, &
                      num_interfaces_inner_core,max_nibool_interfaces_ic, &
                      nibool_interfaces_inner_core,ibool_interfaces_inner_core, &
                      my_neighbors_inner_core, &
                      b_request_send_vector_ic,b_request_recv_vector_ic)
      else
        ! on GPU
        ! sends accel values to corresponding MPI interface neighbors

        ! preparation of the contribution between partitions using MPI
        ! transfers MPI buffers to CPU
        ! note: in case of asynchronous copy, this transfers boundary region to host asynchronously. The
        !       MPI-send is done after compute_forces_viscoelastic_gpu,
        !       once the inner element kernels are launched, and the memcpy has finished.
        call transfer_boun_from_device(Mesh_pointer, &
                                       b_buffer_send_vector_cm, &
                                       IREGION_CRUST_MANTLE,3)
        call transfer_boun_from_device(Mesh_pointer, &
                                       b_buffer_send_vector_inner_core, &
                                       IREGION_INNER_CORE,3)

        if (.not. GPU_ASYNC_COPY) then
          ! for synchronous transfers, sending over MPI can directly proceed
          ! crust mantle
          call assemble_MPI_vector_send_gpu(NPROCTOT_VAL, &
                      b_buffer_send_vector_cm,b_buffer_recv_vector_cm, &
                      num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                      nibool_interfaces_crust_mantle, &
                      my_neighbors_crust_mantle, &
                      b_request_send_vector_cm,b_request_recv_vector_cm)
          ! inner core
          call assemble_MPI_vector_send_gpu(NPROCTOT_VAL, &
                      b_buffer_send_vector_inner_core,b_buffer_recv_vector_inner_core, &
                      num_interfaces_inner_core,max_nibool_interfaces_ic, &
                      nibool_interfaces_inner_core, &
                      my_neighbors_inner_core, &
                      b_request_send_vector_ic,b_request_recv_vector_ic)
        endif
      endif ! GPU
    else
      ! adjoint / kernel runs
      ! waits for send/receive requests to be completed and assembles values
      if (.not. GPU_MODE) then
        ! on CPU
        ! crust mantle
        call assemble_MPI_vector_w(NPROCTOT_VAL,NGLOB_CRUST_MANTLE, &
                            b_accel_crust_mantle, &
                            b_buffer_recv_vector_cm,num_interfaces_crust_mantle, &
                            max_nibool_interfaces_cm, &
                            nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
                            my_neighbors_crust_mantle, &
                            b_request_send_vector_cm,b_request_recv_vector_cm)
        ! inner core
        call assemble_MPI_vector_w(NPROCTOT_VAL,NGLOB_INNER_CORE, &
                            b_accel_inner_core, &
                            b_buffer_recv_vector_inner_core,num_interfaces_inner_core, &
                            max_nibool_interfaces_ic, &
                            nibool_interfaces_inner_core,ibool_interfaces_inner_core, &
                            my_neighbors_inner_core, &
                            b_request_send_vector_ic,b_request_recv_vector_ic)

      else
        ! on GPU
        if (GPU_ASYNC_COPY) then
          ! while inner elements compute "Kernel_2", we wait for MPI to
          ! finish and transfer the boundary terms to the device asynchronously
          ! wait for asynchronous copy to finish
          !
          ! transfers MPI buffers onto GPU
          ! crust/mantle region
          call transfer_boundary_to_device(Mesh_pointer,NPROCTOT_VAL,b_buffer_recv_vector_cm, &
                                           num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                                           b_request_recv_vector_cm, &
                                           IREGION_CRUST_MANTLE,3)
          ! inner core region
          call transfer_boundary_to_device(Mesh_pointer,NPROCTOT_VAL,b_buffer_recv_vector_inner_core, &
                                           num_interfaces_inner_core,max_nibool_interfaces_ic, &
                                           b_request_recv_vector_ic, &
                                           IREGION_INNER_CORE,3)
        endif

        ! waits for MPI send/receive requests to be completed and assembles values
        ! crust mantle
        call assemble_MPI_vector_write_gpu(Mesh_pointer,NPROCTOT_VAL, &
                          b_buffer_recv_vector_cm, &
                          num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                          b_request_send_vector_cm,b_request_recv_vector_cm, &
                          IREGION_CRUST_MANTLE,3) ! -- 3 == adjoint b_accel
        ! inner core
        call assemble_MPI_vector_write_gpu(Mesh_pointer,NPROCTOT_VAL, &
                          b_buffer_recv_vector_inner_core, &
                          num_interfaces_inner_core,max_nibool_interfaces_ic, &
                          b_request_send_vector_ic,b_request_recv_vector_ic, &
                          IREGION_INNER_CORE,3)
      endif
    endif ! iphase == 1

  enddo ! iphase

  ! updates (only) acceleration w/ rotation in the crust/mantle and inner core region
  if (.not. GPU_MODE) then
    ! on CPU
    ! adjoint / kernel runs
    ! crust/mantle and inner core region
    call multiply_accel_elastic(b_two_omega_earth, &
                                NGLOB_CRUST_MANTLE_ADJOINT,b_veloc_crust_mantle,b_accel_crust_mantle, &
                                b_rmassx_crust_mantle,b_rmassy_crust_mantle,b_rmassz_crust_mantle, &
                                NGLOB_INNER_CORE,b_veloc_inner_core,b_accel_inner_core, &
                                b_rmassx_inner_core,b_rmassy_inner_core,b_rmassz_inner_core)
  else
     ! on GPU
     ! includes FORWARD_OR_ADJOINT == 3
     call multiply_accel_elastic_gpu(Mesh_pointer,3)
  endif

  ! couples ocean with crust mantle
  ! (updates acceleration with ocean load approximation)
  if (OCEANS_VAL) then
    if (.not. GPU_MODE) then
      ! on CPU
      call compute_coupling_ocean(NGLOB_CRUST_MANTLE_ADJOINT,b_accel_crust_mantle, &
                                  b_rmassx_crust_mantle,b_rmassy_crust_mantle,b_rmassz_crust_mantle, &
                                  npoin_oceans,ibool_ocean_load,rmass_ocean_load_selected,normal_ocean_load)
    else
      ! on GPU
      call compute_coupling_ocean_gpu(Mesh_pointer,3) ! -- 3 == backward/reconstructed arrays
    endif
  endif

  ! time scheme update
  if (USE_LDDRK) then
    ! Runge-Kutta scheme
    call update_veloc_elastic_lddrk_backward()
  else
    ! Newmark time scheme
    call update_veloc_elastic_newmark_backward()
  endif

!daniel debug: att - debug
!  if (DEBUG) then
!    if (SIMULATION_TYPE == 1) then
!      if (it > NSTEP - 1000 .and. myrank == 0) then
!        print *,'it',it,'Rxx:',R_xx_crust_mantle(1,1,1,1,100),epsilondev_xx_crust_mantle(1,1,1,100)
!      endif
!    else if (SIMULATION_TYPE == 3) then
!      if (it <= 1000 .and. myrank == 0) then
!        print *,'it',it,'Rxx:',b_R_xx_crust_mantle(1,1,1,1,100),b_epsilondev_xx_crust_mantle(1,1,1,100)
!      endif
!    endif
!  endif

  end subroutine compute_forces_viscoelastic_backward


!
!-------------------------------------------------------------------------------------------------
!


  subroutine compute_forces_crust_mantle( NSPEC_STR_OR_ATT,NGLOB,NSPEC_ATT, &
                                          deltat, &
                                          displ_crust_mantle, &
                                          accel_crust_mantle, &
                                          iphase, &
                                          R_xx,R_yy,R_xy,R_xz,R_yz, &
                                          R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                          epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                          epsilondev_xz,epsilondev_yz, &
                                          epsilon_trace_over_3, &
                                          alphaval,betaval,gammaval)

! wrapper function, decides about Deville optimization
!
! (left in this file to let compiler decide about inlining)

  use constants_solver, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ,USE_DEVILLE_PRODUCTS_VAL, &
    ATT4_VAL,N_SLS,NSPEC_CRUST_MANTLE_STRAIN_ONLY

  ! note: passes sum_terms array as subroutine argument which will help for performance (better than use-statement)
  use specfem_par_crustmantle, only: sum_terms_crust_mantle,factor_common_crust_mantle

  implicit none

  integer,intent(in) :: NSPEC_STR_OR_ATT,NGLOB,NSPEC_ATT

  ! time step
  real(kind=CUSTOM_REAL),intent(in) :: deltat

  ! displacement, velocity and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB),intent(in) :: displ_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB),intent(inout) :: accel_crust_mantle

  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATT),intent(inout) :: R_xx,R_yy,R_xy,R_xz,R_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATT),intent(inout) :: &
    R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STR_OR_ATT),intent(inout) :: &
    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STRAIN_ONLY),intent(inout) :: epsilon_trace_over_3

  ! [alpha,beta,gamma]val reduced to N_SLS and factor_common to N_SLS*NUM_NODES
  real(kind=CUSTOM_REAL), dimension(N_SLS),intent(in) :: alphaval,betaval,gammaval

  ! inner/outer element run flag
  integer,intent(in) :: iphase

  if (USE_DEVILLE_PRODUCTS_VAL) then
    ! uses Deville (2002) optimizations
    call compute_forces_crust_mantle_Dev(NSPEC_STR_OR_ATT,NGLOB,NSPEC_ATT, &
                                         deltat, &
                                         displ_crust_mantle, &
                                         accel_crust_mantle, &
                                         iphase, &
                                         R_xx,R_yy,R_xy,R_xz,R_yz, &
                                         R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                         epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
                                         epsilon_trace_over_3, &
                                         alphaval,betaval,gammaval, &
                                         factor_common_crust_mantle,ATT4_VAL,sum_terms_crust_mantle)
  else
    ! no Deville optimization
    call compute_forces_crust_mantle_noDev(NSPEC_STR_OR_ATT,NGLOB,NSPEC_ATT, &
                                           deltat, &
                                           displ_crust_mantle, &
                                           accel_crust_mantle, &
                                           iphase, &
                                           R_xx,R_yy,R_xy,R_xz,R_yz, &
                                           R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                           epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
                                           epsilon_trace_over_3, &
                                           alphaval,betaval,gammaval, &
                                           factor_common_crust_mantle,ATT4_VAL)
  endif

  end subroutine compute_forces_crust_mantle

!
!-------------------------------------------------------------------------------------------------
!


  subroutine compute_forces_inner_core( NSPEC_STR_OR_ATT,NGLOB,NSPEC_ATT, &
                                        deltat, &
                                        displ_inner_core, &
                                        accel_inner_core, &
                                        iphase, &
                                        R_xx,R_yy,R_xy,R_xz,R_yz, &
                                        R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                        epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                        epsilondev_xz,epsilondev_yz, &
                                        epsilon_trace_over_3, &
                                        alphaval,betaval,gammaval)

! wrapper function, decides about Deville optimization
!
! (left in this file to let compiler decide about inlining)

  use constants_solver, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ,USE_DEVILLE_PRODUCTS_VAL, &
    ATT5_VAL,N_SLS,NSPEC_INNER_CORE_STRAIN_ONLY

  ! note: passes sum_terms array as subroutine argument which will help for performance (better than use-statement)
  use specfem_par_innercore, only: sum_terms_inner_core,factor_common_inner_core

  implicit none

  integer,intent(in) :: NSPEC_STR_OR_ATT,NGLOB,NSPEC_ATT

  ! time step
  real(kind=CUSTOM_REAL),intent(in) :: deltat

  ! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB),intent(in) :: displ_inner_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB),intent(inout) :: accel_inner_core

  ! for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  ! variable lengths for factor_common and one_minus_sum_beta
  real(kind=CUSTOM_REAL), dimension(N_SLS),intent(in) :: alphaval,betaval,gammaval

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATT),intent(inout) :: R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATT),intent(inout) :: &
    R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STR_OR_ATT),intent(inout) :: &
    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STRAIN_ONLY),intent(inout) :: epsilon_trace_over_3

  ! inner/outer element run flag
  integer,intent(in) :: iphase

  if (USE_DEVILLE_PRODUCTS_VAL) then
    ! uses Deville (2002) optimizations
    call compute_forces_inner_core_Dev(NSPEC_STR_OR_ATT,NGLOB,NSPEC_ATT, &
                                       deltat, &
                                       displ_inner_core, &
                                       accel_inner_core, &
                                       iphase, &
                                       R_xx,R_yy,R_xy,R_xz,R_yz, &
                                       R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                       epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
                                       epsilon_trace_over_3, &
                                       alphaval,betaval,gammaval, &
                                       factor_common_inner_core,ATT5_VAL,sum_terms_inner_core)
  else
    ! no Deville optimization
    call compute_forces_inner_core_noDev(NSPEC_STR_OR_ATT,NGLOB,NSPEC_ATT, &
                                         deltat, &
                                         displ_inner_core, &
                                         accel_inner_core, &
                                         iphase, &
                                         R_xx,R_yy,R_xy,R_xz,R_yz, &
                                         R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                         epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
                                         epsilon_trace_over_3, &
                                         alphaval,betaval,gammaval, &
                                         factor_common_inner_core,ATT5_VAL)
  endif

  end subroutine compute_forces_inner_core

