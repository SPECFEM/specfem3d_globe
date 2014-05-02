!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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


  subroutine iterate_time_undoatt()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_movie
  use write_seismograms_mod, only: write_seismograms
  implicit none

  ! local parameters
  integer :: it_temp,seismo_current_temp
  integer :: i,j,ier
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_displ_crust_mantle_store_buffer
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_displ_inner_core_store_buffer
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_displ_outer_core_store_buffer,b_accel_outer_core_store_buffer

  ! timing
  double precision, external :: wtime

!
!   s t a r t   t i m e   i t e r a t i o n s
!

  ! checks
  if( .not. UNDO_ATTENUATION ) return

  ! allocates buffers
  if( SIMULATION_TYPE == 3 ) then
    !! DK DK to Daniel, July 2013: in the case of GPU_MODE it will be *crucial* to leave these arrays on the host
    !! i.e. on the CPU, in order to be able to use all the (unused) memory of the host for them, since they are
    !! (purposely) huge and designed to use almost all the memory available (by carefully optimizing the
    !! value of NT_DUMP_ATTENUATION); when writing to these buffers, it will then be OK to use non-blocking writes
    !! from the device to the host, since we do not reuse the content of these buffers until much later, in a second part
    !! of the run; however when reading back from these buffers, the reads from host to device should then be blocking
    !! because we then use the values read immediately (one value at a time, but to reduce the total number of reads
    !! across the PCI-Express bus we could / should consider reading them 10 by 10 for instance (?) if that fits
    !! in the memory of the GPU
    allocate(b_displ_crust_mantle_store_buffer(NDIM,NGLOB_CRUST_MANTLE_ADJOINT,NT_DUMP_ATTENUATION),stat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error allocating b_displ_crust_mantle_store_buffer')
    allocate(b_displ_outer_core_store_buffer(NGLOB_OUTER_CORE_ADJOINT,NT_DUMP_ATTENUATION),stat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error allocating b_displ_outer_core_store_buffer')
    allocate(b_accel_outer_core_store_buffer(NGLOB_OUTER_CORE_ADJOINT,NT_DUMP_ATTENUATION),stat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error allocating b_displ_outer_core_store_buffer')
    allocate(b_displ_inner_core_store_buffer(NDIM,NGLOB_INNER_CORE_ADJOINT,NT_DUMP_ATTENUATION),stat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error allocating b_displ_inner_core_store_buffer')
  endif

  ! synchronize all processes to make sure everybody is ready to start time loop
  call synchronize_all()

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Starting time iteration loop in undoing attenuation...'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! create an empty file to monitor the start of the simulation
  if(myrank == 0) then
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/starttimeloop_undoatt.txt',status='unknown',action='write')
    write(IOUT,*) 'hello, starting time loop'
    close(IOUT)
  endif

  ! initialize variables for writing seismograms
  seismo_offset = it_begin-1
  seismo_current = 0

  ! get MPI starting time
  time_start = wtime()

  ! *********************************************************
  ! ************* MAIN LOOP OVER THE TIME STEPS *************
  ! *********************************************************

  it = 0
  do iteration_on_subset = 1, NSTEP / NT_DUMP_ATTENUATION

    ! wavefield storage
    if( SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
      ! saves forward wavefields
      call save_forward_arrays_undoatt()

    else if( SIMULATION_TYPE == 3 ) then
      ! reads in last stored forward wavefield
      call read_forward_arrays_undoatt()

      ! note: after reading the restart files of displacement back from disk, recompute the strain from displacement;
      !       this is better than storing the strain to disk as well, which would drastically increase I/O volume
      ! computes strain based on current backward/reconstructed wavefield
      if(COMPUTE_AND_STORE_STRAIN) call itu_compute_strain_att_backward()

      ! intermediate storage of it and seismo_current positions
      it_temp = it
      seismo_current_temp = seismo_current
    endif

    ! time loop within this iteration subset
    select case( SIMULATION_TYPE )
    case( 1, 2 )
      ! forward and adjoint simulations

      ! subset loop
      do it_of_this_subset = 1, NT_DUMP_ATTENUATION

        it = it + 1

        ! simulation status output and stability check
        if( mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == it_begin + 4 .or. it == it_end ) then
          call check_stability()
        endif

        do istage = 1, NSTAGE_TIME_SCHEME ! is equal to 1 if Newmark because only one stage then

          if(USE_LDDRK)then
            ! update displacement using Runge-Kutta time scheme
            call update_displacement_lddrk()
          else
            ! update displacement using Newmark time scheme
            call update_displacement_Newmark()
          endif

          ! acoustic solver for outer core
          ! (needs to be done first, before elastic one)
          call compute_forces_acoustic()

          ! elastic solver for crust/mantle and inner core
          call compute_forces_viscoelastic()

        enddo ! istage

        ! write the seismograms with time shift
        if( nrec_local > 0 .or. ( WRITE_SEISMOGRAMS_BY_MASTER .and. myrank == 0 ) ) then
          call write_seismograms()
        endif

        ! outputs movie files
        call write_movie_output()

        ! first step of noise tomography, i.e., save a surface movie at every time step
        ! modified from the subroutine 'write_movie_surface'
        if( NOISE_TOMOGRAPHY == 1 ) then
          call noise_save_surface_movie()
        endif

        ! updates VTK window
        if( VTK_MODE ) then
          call it_update_vtkwindow()
        endif

      enddo ! subset loop

    case( 3 )
      ! kernel simulations

      ! reconstructs forward wavefield based on last stored wavefield data
      !
      ! note: we step forward in time here, starting from last snapshot.
      !       the newly computed, reconstructed forward wavefields (b_displ_..) get stored in buffers.

      ! subset loop
      do it_of_this_subset = 1, NT_DUMP_ATTENUATION

        it = it + 1

        ! simulation status output and stability check
        if( mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == it_begin + 4 .or. it == it_end ) then
          call check_stability_backward()
        endif

        do istage = 1, NSTAGE_TIME_SCHEME ! is equal to 1 if Newmark because only one stage then

          if(USE_LDDRK)then
            ! update displacement using Runge-Kutta time scheme
            call update_displacement_lddrk_backward()
          else
            ! update displacement using Newmark time scheme
            call update_displacement_Newmark_backward()
          endif

          ! acoustic solver for outer core
          ! (needs to be done first, before elastic one)
          call compute_forces_acoustic_backward()

          ! elastic solver for crust/mantle and inner core
          call compute_forces_viscoelastic_backward()

        enddo ! istage

        ! stores wavefield in buffers
        b_displ_crust_mantle_store_buffer(:,:,it_of_this_subset) = b_displ_crust_mantle(:,:)
        b_displ_outer_core_store_buffer(:,it_of_this_subset) = b_displ_outer_core(:)
        b_accel_outer_core_store_buffer(:,it_of_this_subset) = b_accel_outer_core(:)
        b_displ_inner_core_store_buffer(:,:,it_of_this_subset) = b_displ_inner_core(:,:)

      enddo ! subset loop

      it = it_temp
      seismo_current = seismo_current_temp

      ! computes strain based on current adjoint wavefield
      if(COMPUTE_AND_STORE_STRAIN) call itu_compute_strain_att()

      ! adjoint wavefield simulation
      do it_of_this_subset = 1, NT_DUMP_ATTENUATION

        ! reads backward/reconstructed wavefield from buffers
        ! note: uses wavefield at corresponding time (NSTEP - it + 1 ), i.e. we have now time-reversed wavefields
        ! crust/mantle
        do i = 1, NDIM
          do j =1,NGLOB_CRUST_MANTLE_ADJOINT
            b_displ_crust_mantle(i,j) = b_displ_crust_mantle_store_buffer(i,j,NT_DUMP_ATTENUATION-it_of_this_subset+1)
          enddo
        enddo
        ! outer core
        do j =1,NGLOB_OUTER_CORE_ADJOINT
            b_displ_outer_core(j) = b_displ_outer_core_store_buffer(j,NT_DUMP_ATTENUATION-it_of_this_subset+1)
            b_accel_outer_core(j) = b_accel_outer_core_store_buffer(j,NT_DUMP_ATTENUATION-it_of_this_subset+1)
        enddo
        ! inner core
        do i = 1, NDIM
          do j =1,NGLOB_INNER_CORE_ADJOINT
            b_displ_inner_core(i,j) = b_displ_inner_core_store_buffer(i,j,NT_DUMP_ATTENUATION-it_of_this_subset+1)
          enddo
        enddo

        it = it + 1

        ! simulation status output and stability check
        if( mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == it_begin + 4 .or. it == it_end ) then
          call check_stability()
        endif

        do istage = 1, NSTAGE_TIME_SCHEME ! is equal to 1 if Newmark because only one stage then

          if(USE_LDDRK)then
            ! update displacement using Runge-Kutta time scheme
            call update_displacement_lddrk()
          else
            ! update displacement using Newmark time scheme
            call update_displacement_Newmark()
          endif

          ! acoustic solver for outer core
          ! (needs to be done first, before elastic one)
          call compute_forces_acoustic()

          ! elastic solver for crust/mantle and inner core
          call compute_forces_viscoelastic()

        enddo ! istage

        ! write the seismograms with time shift
        if( nrec_local > 0 .or. ( WRITE_SEISMOGRAMS_BY_MASTER .and. myrank == 0 ) ) then
          call write_seismograms()
        endif

        ! kernel computation
        ! adjoint simulations: kernels
        call compute_kernels()

      enddo ! subset loop

    end select ! SIMULATION_TYPE

  enddo   ! end of main time loop

  !
  !---- end of time iteration loop
  !

  ! frees undo_attenuation buffers
  if( SIMULATION_TYPE == 3 ) then
    deallocate(b_displ_crust_mantle_store_buffer, &
               b_displ_outer_core_store_buffer, &
               b_accel_outer_core_store_buffer, &
               b_displ_inner_core_store_buffer)
  endif

  call print_elapsed_time()

  ! Transfer fields from GPU card to host for further analysis
  if(GPU_MODE) call it_transfer_from_GPU()

  end subroutine iterate_time_undoatt


!
!--------------------------------------------------------------------------------------------
!
! strain for whole domain crust/mantle and inner core
!
!--------------------------------------------------------------------------------------------
!

  subroutine itu_compute_strain_att()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore

  implicit none
  ! local parameters
  integer :: ispec

  ! checks
  if( USE_DEVILLE_PRODUCTS_VAL ) then

    ! inner core
    do ispec = 1, NSPEC_INNER_CORE
      call compute_element_strain_att_Dev(ispec,NGLOB_INNER_CORE,NSPEC_INNER_CORE, &
                                          displ_inner_core,veloc_inner_core,0._CUSTOM_REAL, &
                                          ibool_inner_core, &
                                          hprime_xx,hprime_xxT, &
                                          xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                          etax_inner_core,etay_inner_core,etaz_inner_core, &
                                          gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                                          epsilondev_xx_inner_core, &
                                          epsilondev_yy_inner_core, &
                                          epsilondev_xy_inner_core, &
                                          epsilondev_xz_inner_core, &
                                          epsilondev_yz_inner_core, &
                                          NSPEC_INNER_CORE_STRAIN_ONLY,eps_trace_over_3_inner_core)
    enddo
    ! crust mantle
    do ispec = 1, NSPEC_crust_mantle
      call compute_element_strain_att_Dev(ispec,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
                                          displ_crust_mantle,veloc_crust_mantle,0._CUSTOM_REAL, &
                                          ibool_crust_mantle, &
                                          hprime_xx,hprime_xxT, &
                                          xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                          etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                          gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                          epsilondev_xx_crust_mantle, &
                                          epsilondev_yy_crust_mantle, &
                                          epsilondev_xy_crust_mantle, &
                                          epsilondev_xz_crust_mantle, &
                                          epsilondev_yz_crust_mantle, &
                                          NSPEC_CRUST_MANTLE_STRAIN_ONLY,eps_trace_over_3_crust_mantle)
    enddo

  else

    ! inner core
    do ispec = 1, NSPEC_INNER_CORE
      call compute_element_strain_att_noDev(ispec,NGLOB_INNER_CORE,NSPEC_INNER_CORE, &
                                            displ_inner_core,veloc_inner_core,0._CUSTOM_REAL, &
                                            ibool_inner_core, &
                                            hprime_xx,hprime_yy,hprime_zz, &
                                            xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                            etax_inner_core,etay_inner_core,etaz_inner_core, &
                                            gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                                            epsilondev_xx_inner_core, &
                                            epsilondev_yy_inner_core, &
                                            epsilondev_xy_inner_core, &
                                            epsilondev_xz_inner_core, &
                                            epsilondev_yz_inner_core, &
                                            NSPEC_INNER_CORE_STRAIN_ONLY,eps_trace_over_3_inner_core)
    enddo
    ! crust mantle
    do ispec = 1, NSPEC_CRUST_MANTLE
      call compute_element_strain_att_noDev(ispec,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
                                            displ_crust_mantle,veloc_crust_mantle,0._CUSTOM_REAL, &
                                            ibool_crust_mantle, &
                                            hprime_xx,hprime_yy,hprime_zz, &
                                            xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                            etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                            gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                            epsilondev_xx_crust_mantle, &
                                            epsilondev_yy_crust_mantle, &
                                            epsilondev_xy_crust_mantle, &
                                            epsilondev_xz_crust_mantle, &
                                            epsilondev_yz_crust_mantle, &
                                            NSPEC_CRUST_MANTLE_STRAIN_ONLY,eps_trace_over_3_crust_mantle)
    enddo
  endif

  end subroutine itu_compute_strain_att

!
!--------------------------------------------------------------------------------------------
!

  subroutine itu_compute_strain_att_backward()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore

  implicit none
  ! local parameters
  integer :: ispec

  ! checks
  if( USE_DEVILLE_PRODUCTS_VAL ) then

    ! inner core
    do ispec = 1, NSPEC_INNER_CORE
      call compute_element_strain_att_Dev(ispec,NGLOB_INNER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
                                          b_displ_inner_core,b_veloc_inner_core,0._CUSTOM_REAL, &
                                          ibool_inner_core, &
                                          hprime_xx,hprime_xxT, &
                                          xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                          etax_inner_core,etay_inner_core,etaz_inner_core, &
                                          gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                                          b_epsilondev_xx_inner_core, &
                                          b_epsilondev_yy_inner_core, &
                                          b_epsilondev_xy_inner_core, &
                                          b_epsilondev_xz_inner_core, &
                                          b_epsilondev_yz_inner_core, &
                                          NSPEC_INNER_CORE_STRAIN_ONLY,b_eps_trace_over_3_inner_core)
    enddo

    ! crust mantle
    do ispec = 1, NSPEC_CRUST_MANTLE
      call compute_element_strain_att_Dev(ispec,NGLOB_CRUST_MANTLE_ADJOINT,NSPEC_CRUST_MANTLE_ADJOINT, &
                                          b_displ_crust_mantle,b_veloc_crust_mantle,0._CUSTOM_REAL, &
                                          ibool_crust_mantle, &
                                          hprime_xx,hprime_xxT, &
                                          xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                          etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                          gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                          b_epsilondev_xx_crust_mantle, &
                                          b_epsilondev_yy_crust_mantle, &
                                          b_epsilondev_xy_crust_mantle, &
                                          b_epsilondev_xz_crust_mantle, &
                                          b_epsilondev_yz_crust_mantle, &
                                          NSPEC_CRUST_MANTLE_STRAIN_ONLY,b_eps_trace_over_3_crust_mantle)
    enddo

  else

    ! inner core
    do ispec = 1, NSPEC_INNER_CORE
      call compute_element_strain_att_noDev(ispec,NGLOB_INNER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
                                            b_displ_inner_core,b_veloc_inner_core,0._CUSTOM_REAL, &
                                            ibool_inner_core, &
                                            hprime_xx,hprime_yy,hprime_zz, &
                                            xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                            etax_inner_core,etay_inner_core,etaz_inner_core, &
                                            gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                                            b_epsilondev_xx_inner_core, &
                                            b_epsilondev_yy_inner_core, &
                                            b_epsilondev_xy_inner_core, &
                                            b_epsilondev_xz_inner_core, &
                                            b_epsilondev_yz_inner_core, &
                                            NSPEC_INNER_CORE_STRAIN_ONLY,b_eps_trace_over_3_inner_core)
    enddo
    ! crust mantle
    do ispec = 1, NSPEC_crust_mantle
      call compute_element_strain_att_noDev(ispec,NGLOB_CRUST_MANTLE_ADJOINT,NSPEC_CRUST_MANTLE_ADJOINT, &
                                            b_displ_crust_mantle,b_veloc_crust_mantle,0._CUSTOM_REAL, &
                                            ibool_crust_mantle, &
                                            hprime_xx,hprime_yy,hprime_zz, &
                                            xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                            etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                            gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                            b_epsilondev_xx_crust_mantle, &
                                            b_epsilondev_yy_crust_mantle, &
                                            b_epsilondev_xy_crust_mantle, &
                                            b_epsilondev_xz_crust_mantle, &
                                            b_epsilondev_yz_crust_mantle, &
                                            NSPEC_CRUST_MANTLE_STRAIN_ONLY,b_eps_trace_over_3_crust_mantle)
    enddo
  endif

  end subroutine itu_compute_strain_att_backward

