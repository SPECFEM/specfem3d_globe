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
!
! United States and French Government Sponsorship Acknowledged.

  subroutine iterate_time()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_movie
  use write_seismograms_mod, only: write_seismograms
  implicit none

  ! timing
  double precision, external :: wtime

  ! for EXACT_UNDOING_TO_DISK
  integer :: ispec,iglob,i,j,k,counter,record_length
  real(kind=CUSTOM_REAL) :: radius
  integer, dimension(:), allocatable :: integer_mask_ibool_exact_undo
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: buffer_for_disk
  character(len=MAX_STRING_LEN) outputname

  !----  create a Gnuplot script to display the energy curve in log scale
  if (OUTPUT_ENERGY .and. myrank == 0) then
    open(unit=IOUT_ENERGY,file=trim(OUTPUT_FILES)//'plot_energy.gnu',status='unknown',action='write')
    write(IOUT_ENERGY,*) 'set term wxt'
    write(IOUT_ENERGY,*) '#set term postscript landscape color solid "Helvetica" 22'
    write(IOUT_ENERGY,*) '#set output "energy.ps"'
    write(IOUT_ENERGY,*) 'set logscale y'
    write(IOUT_ENERGY,*) 'set xlabel "Time step number"'
    write(IOUT_ENERGY,*) 'set ylabel "Energy (J)"'
    write(IOUT_ENERGY,'(a152)') '#plot "energy.dat" us 1:2 t ''Kinetic Energy'' w l lc 1, "energy.dat" us 1:3 &
                         &t ''Potential Energy'' w l lc 2, "energy.dat" us 1:4 t ''Total Energy'' w l lc 4'
    write(IOUT_ENERGY,*) '#pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) '#plot "energy.dat" us 1:2 t ''Kinetic Energy'' w l lc 1'
    write(IOUT_ENERGY,*) '#pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) '#plot "energy.dat" us 1:3 t ''Potential Energy'' w l lc 2'
    write(IOUT_ENERGY,*) '#pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) 'plot "energy.dat" us 1:4 t ''Total Energy'' w l lc 4'
    write(IOUT_ENERGY,*) 'pause -1 "Hit any key..."'
    close(IOUT_ENERGY)
  endif

  ! open the file in which we will store the energy curve
  if (OUTPUT_ENERGY .and. myrank == 0) &
    open(unit=IOUT_ENERGY,file=trim(OUTPUT_FILES)//'energy.dat',status='unknown',action='write')

!
!   s t a r t   t i m e   i t e r a t i o n s
!

  ! synchronize all processes to make sure everybody is ready to start time loop
  call synchronize_all()
  if (myrank == 0) write(IMAIN,*) 'All processes are synchronized before time loop'

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Starting time iteration loop...'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! create an empty file to monitor the start of the simulation
  if (myrank == 0) then
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/starttimeloop.txt',status='unknown',action='write')
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

  if (EXACT_UNDOING_TO_DISK) then

    if (GPU_MODE) call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK not supported for GPUs')

    if (UNDO_ATTENUATION) &
      call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK needs UNDO_ATTENUATION to be off because it computes the kernel directly instead')

    if (SIMULATION_TYPE == 1 .and. .not. SAVE_FORWARD) &
      call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK requires SAVE_FORWARD if SIMULATION_TYPE == 1')

    if (ANISOTROPIC_KL) call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK requires ANISOTROPIC_KL to be turned off')

!! DK DK determine the largest value of iglob that we need to save to disk,
!! DK DK since we save the upper part of the mesh only in the case of surface-wave kernels
    ! crust_mantle
    allocate(integer_mask_ibool_exact_undo(NGLOB_CRUST_MANTLE))
    integer_mask_ibool_exact_undo(:) = -1

    counter = 0
    do ispec = 1, NSPEC_CRUST_MANTLE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_crust_mantle(i,j,k,ispec)
            ! xstore ystore zstore have previously been converted to r theta phi, thus xstore now stores the radius
            radius = xstore_crust_mantle(iglob) ! <- radius r (normalized)
            ! save that element only if it is in the upper part of the mesh
            if (radius >= R670 / R_EARTH) then
              ! if this point has not yet been found before
              if (integer_mask_ibool_exact_undo(iglob) == -1) then
                ! create a new unique point
                counter = counter + 1
                integer_mask_ibool_exact_undo(iglob) = counter
              endif
            endif
          enddo
        enddo
      enddo
    enddo

    ! allocate the buffer used to dump a single time step
    allocate(buffer_for_disk(counter))

    ! open the file in which we will dump all the time steps (in a single file)
    write(outputname,"('huge_dumps/proc',i6.6,'_huge_dump_of_all_time_steps.bin')") myrank
    inquire(iolength=record_length) buffer_for_disk
    ! we write to or read from the file depending on the simulation type
    if (SIMULATION_TYPE == 1) then
      open(file=outputname, unit=IFILE_FOR_EXACT_UNDOING, action='write', status='unknown', &
                      form='unformatted', access='direct', recl=record_length)
    else if (SIMULATION_TYPE == 3) then
      open(file=outputname, unit=IFILE_FOR_EXACT_UNDOING, action='read', status='old', &
                      form='unformatted', access='direct', recl=record_length)
    else
      call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK can only be used with SIMULATION_TYPE == 1 or SIMULATION_TYPE == 3')
    endif

  endif ! of if (EXACT_UNDOING_TO_DISK)

  do it = it_begin,it_end

    ! simulation status output and stability check
    if (mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == it_begin + 4 .or. it == it_end) then
      call check_stability()
    endif

    do istage = 1, NSTAGE_TIME_SCHEME ! is equal to 1 if Newmark because only one stage then

      if (USE_LDDRK) then
        ! update displacement using Runge-Kutta time scheme
        call update_displ_lddrk()
      else
        ! update displacement using Newmark time scheme
        call update_displ_Newmark()
      endif

      ! acoustic solver for outer core
      ! (needs to be done first, before elastic one)
      call compute_forces_acoustic()

      ! elastic solver for crust/mantle and inner core
      call compute_forces_viscoelastic()

    enddo ! end of very big external loop on istage for all the stages of the LDDRK time scheme (only one stage if Newmark)

    ! save the forward run to disk for the alpha kernel only
    if (EXACT_UNDOING_TO_DISK .and. SIMULATION_TYPE == 1) then
      do ispec = 1, NSPEC_CRUST_MANTLE
        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX
              iglob = ibool_crust_mantle(i,j,k,ispec)
              if (integer_mask_ibool_exact_undo(iglob) /= -1) &
                buffer_for_disk(integer_mask_ibool_exact_undo(iglob)) = eps_trace_over_3_crust_mantle(i,j,k,ispec)
            enddo
          enddo
        enddo
      enddo
      write(IFILE_FOR_EXACT_UNDOING,rec=it) buffer_for_disk
    endif

    ! kernel simulations (forward and adjoint wavefields)
    if (SIMULATION_TYPE == 3) then

      if (.not. EXACT_UNDOING_TO_DISK) then
      ! note: we step back in time (using time steps - DT ), i.e. wavefields b_displ_..() are time-reversed here

      ! reconstructs forward wavefields based on last stored wavefield data

      ! note: NSTAGE_TIME_SCHEME is equal to 1 if Newmark because only one stage then
      do istage = 1, NSTAGE_TIME_SCHEME

        if (USE_LDDRK) then
          ! update displacement using Runge-Kutta time scheme
          call update_displ_lddrk_backward()
        else
          ! update displacement using Newmark time scheme
          call update_displ_Newmark_backward()
        endif

        ! acoustic solver for outer core
        ! (needs to be done first, before elastic one)
        call compute_forces_acoustic_backward()

        ! elastic solver for crust/mantle and inner core
        call compute_forces_viscoelastic_backward()

      enddo

      ! restores last time snapshot saved for backward/reconstruction of wavefields
      ! note: this is done here after the Newmark time scheme, otherwise the indexing for sources
      !          and adjoint sources will become more complicated
      !          that is, index it for adjoint sources will match index NSTEP - 1 for backward/reconstructed wavefields
      if (it == 1) then
        call read_forward_arrays()
      endif

      else ! of if (.not. EXACT_UNDOING_TO_DISK)

        ! read the forward run from disk for the alpha kernel only
        ! here we time revert the forward run by reading time step NSTEP - it + 1
        read(IFILE_FOR_EXACT_UNDOING,rec=NSTEP-it+1) buffer_for_disk
        do ispec = 1, NSPEC_CRUST_MANTLE
          do k = 1, NGLLZ
            do j = 1, NGLLY
              do i = 1, NGLLX
                iglob = ibool_crust_mantle(i,j,k,ispec)
                if (integer_mask_ibool_exact_undo(iglob) /= -1) &
                  b_eps_trace_over_3_crust_mantle(i,j,k,ispec) = buffer_for_disk(integer_mask_ibool_exact_undo(iglob))
              enddo
            enddo
          enddo
        enddo

      endif ! of if (.not. EXACT_UNDOING_TO_DISK)

    endif ! kernel simulations

    ! calculating gravity field at current timestep
    if (GRAVITY_SIMULATION) call gravity_timeseries()

    ! write the seismograms with time shift (GPU_MODE transfer included)
    call write_seismograms()

    ! adjoint simulations: kernels
    if (SIMULATION_TYPE == 3) then
      call compute_kernels()
    endif

    ! outputs movie files
    call write_movie_output()

    ! first step of noise tomography, i.e., save a surface movie at every time step
    ! modified from the subroutine 'write_movie_surface'
    if (NOISE_TOMOGRAPHY == 1) then
      call noise_save_surface_movie()
    endif

    ! updates VTK window
    if (VTK_MODE) then
      call it_update_vtkwindow()
    endif

  !
  !---- end of time iteration loop
  !
  enddo   ! end of main time loop

  ! close the huge file that contains a dump of all the time steps to disk
  if (EXACT_UNDOING_TO_DISK) close(IFILE_FOR_EXACT_UNDOING)

  call it_print_elapsed_time()

  ! Transfer fields from GPU card to host for further analysis
  if (GPU_MODE) call it_transfer_from_GPU()

!----  close energy file
  if (OUTPUT_ENERGY .and. myrank == 0) close(IOUT_ENERGY)

  end subroutine iterate_time

!
!-------------------------------------------------------------------------------------------------
!

  subroutine it_transfer_from_GPU()

! transfers fields on GPU back onto CPU

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! to store forward wave fields
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then

    call transfer_fields_cm_from_device(NDIM*NGLOB_CRUST_MANTLE, &
                                    displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle, &
                                    Mesh_pointer)
    call transfer_fields_ic_from_device(NDIM*NGLOB_INNER_CORE, &
                                    displ_inner_core,veloc_inner_core,accel_inner_core, &
                                    Mesh_pointer)
    call transfer_fields_oc_from_device(NGLOB_OUTER_CORE, &
                                    displ_outer_core,veloc_outer_core,accel_outer_core, &
                                    Mesh_pointer)

    call transfer_strain_cm_from_device(Mesh_pointer,eps_trace_over_3_crust_mantle, &
                                    epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle, &
                                    epsilondev_xy_crust_mantle,epsilondev_xz_crust_mantle, &
                                    epsilondev_yz_crust_mantle)
    call transfer_strain_ic_from_device(Mesh_pointer,eps_trace_over_3_inner_core, &
                                    epsilondev_xx_inner_core,epsilondev_yy_inner_core, &
                                    epsilondev_xy_inner_core,epsilondev_xz_inner_core, &
                                    epsilondev_yz_inner_core)

    if (ROTATION_VAL) then
      call transfer_rotation_from_device(Mesh_pointer,A_array_rotation,B_array_rotation)
    endif

  else if (SIMULATION_TYPE == 3) then

    ! note: for kernel simulations (SIMULATION_TYPE == 3), attenuation is
    !          only mimicking effects on phase shifts, but not on amplitudes.
    !          flag PARTIAL_PHYS_DISPERSION_ONLY will have to be set to true in this case.
    !
    ! arrays b_R_xx, ... are not used when PARTIAL_PHYS_DISPERSION_ONLY is set,
    ! therefore no need to transfer arrays from GPU to CPU
    !if (ATTENUATION) then
    !endif

    ! to store kernels
    ! inner core
    call transfer_kernels_ic_to_host(Mesh_pointer, &
                                     rho_kl_inner_core, &
                                     alpha_kl_inner_core, &
                                     beta_kl_inner_core,NSPEC_INNER_CORE)
    ! outer core
    call transfer_kernels_oc_to_host(Mesh_pointer, &
                                     rho_kl_outer_core,&
                                     alpha_kl_outer_core,NSPEC_OUTER_CORE)
    ! crust/mantle
    call transfer_kernels_cm_to_host(Mesh_pointer, &
                                     rho_kl_crust_mantle,alpha_kl_crust_mantle,beta_kl_crust_mantle, &
                                     NSPEC_CRUST_MANTLE)
    ! full anisotropic kernel
    if (ANISOTROPIC_KL) then
      call transfer_kernels_ani_cm_to_host(Mesh_pointer,cijkl_kl_crust_mantle,NSPEC_CRUST_MANTLE)
    endif

    ! specific noise strength kernel
    if (NOISE_TOMOGRAPHY == 3) then
      call transfer_kernels_noise_to_host(Mesh_pointer,Sigma_kl_crust_mantle,NSPEC_CRUST_MANTLE)
    endif

    ! approximative hessian for preconditioning kernels
    if (APPROXIMATE_HESS_KL) then
      call transfer_kernels_hess_cm_tohost(Mesh_pointer,hess_kl_crust_mantle,NSPEC_CRUST_MANTLE)
    endif

  endif

  ! from here on, no gpu data is needed anymore
  ! frees allocated memory on GPU
  call prepare_cleanup_device(Mesh_pointer,NCHUNKS_VAL)

  end subroutine it_transfer_from_GPU

!
!-------------------------------------------------------------------------------------------------
!

  subroutine it_update_vtkwindow()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_movie

  implicit none

  real :: currenttime
  integer :: iglob,inum,data_size
  real, dimension(1) :: dummy

  ! VTK rendering at frame interval
  if (mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then

    ! user output
    !if (myrank == 0 ) print*,"  VTK rendering..."

    ! updates time
    currenttime = sngl((it-1)*DT-t0)

    ! transfers fields from GPU to host
    if (GPU_MODE) then
      !if (myrank == 0 ) print*,"  VTK: transferring velocity from GPU"
      call transfer_veloc_cm_from_device(NDIM*NGLOB_CRUST_MANTLE,veloc_crust_mantle,Mesh_pointer)
    endif

    ! updates wavefield
    !if (myrank == 0 ) print*,"  VTK: it = ",it," out of ",it_end," - norm of velocity field"
    inum = 0
    vtkdata(:) = 0.0
    do iglob = 1,NGLOB_CRUST_MANTLE
      if (vtkmask(iglob) .eqv. .true.) then
        inum = inum + 1
        ! stores norm of velocity vector
        vtkdata(inum) = sqrt(veloc_crust_mantle(1,iglob)**2 &
                           + veloc_crust_mantle(2,iglob)**2 &
                           + veloc_crust_mantle(3,iglob)**2)
      endif
    enddo

    ! updates for multiple MPI process
    if (NPROCTOT_VAL > 1) then
      data_size = size(vtkdata)
      if (myrank == 0) then
        ! gather data
        call gatherv_all_r(vtkdata,data_size,&
                            vtkdata_all,vtkdata_points_all,vtkdata_offset_all, &
                            vtkdata_numpoints_all,NPROCTOT_VAL)
        ! updates VTK window
        call visualize_vtkdata(it,currenttime,vtkdata_all)
      else
        ! all other process just send data
        call gatherv_all_r(vtkdata,data_size,&
                            dummy,vtkdata_points_all,vtkdata_offset_all, &
                            1,NPROCTOT_VAL)
      endif
    else
      ! serial run
      ! updates VTK window
      call visualize_vtkdata(it,currenttime,vtkdata)
    endif

  endif

  end subroutine it_update_vtkwindow

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gravity_timeseries()

  implicit none

  stop 'gravity_timeseries() not implemented in this code yet'

  end subroutine gravity_timeseries

