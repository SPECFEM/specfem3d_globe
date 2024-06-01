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
!
! United States and French Government Sponsorship Acknowledged.

  subroutine iterate_time()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_movie

  implicit none

  ! timing
  double precision, external :: wtime

  ! for EXACT_UNDOING_TO_DISK
  integer :: ispec,iglob,i,j,k

  ! energy curve outputs
  if (OUTPUT_ENERGY) call it_open_energy_curve_file()

!
!   s t a r t   t i m e   i t e r a t i o n s
!

  ! synchronize all processes to make sure everybody is ready to start time loop
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) 'All processes are synchronized before time loop'
    write(IMAIN,*)
    write(IMAIN,*) 'Starting time iteration loop...'
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! create an empty file to monitor the start of the simulation
  if (myrank == 0) then
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/starttimeloop.txt',status='unknown',action='write')
    write(IOUT,*) 'hello, starting time loop'
    close(IOUT)
  endif

  ! synchronizes GPU kernels
  if (GPU_MODE) call gpu_synchronize()

  ! initialize variables for writing seismograms
  seismo_offset = it_begin-1
  seismo_current = 0

  ! get MPI starting time
  time_start = wtime()

  ! *********************************************************
  ! ************* MAIN LOOP OVER THE TIME STEPS *************
  ! *********************************************************

  if (EXACT_UNDOING_TO_DISK) call setup_exact_undoing_to_disk()

  ! time loop
  do it = it_begin,it_end

    ! simulation status output and stability check
    if (mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == it_begin + 4 .or. it == it_end) then
      call check_stability()
      if (I_am_running_on_a_slow_node) goto 100
    endif

    do istage = 1, NSTAGE_TIME_SCHEME ! is equal to 1 if Newmark because only one stage then

      if (USE_LDDRK) then
        ! update displacement using Runge-Kutta time scheme
        call update_displ_lddrk()
      else
        ! update displacement using Newmark time scheme
        call update_displ_Newmark()
      endif

      ! update Poisson's load and solve Poisson's equations
      if (FULL_GRAVITY) call SIEM_solve_poisson()

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
                buffer_for_disk(integer_mask_ibool_exact_undo(iglob),it_of_this_exact_subset) = &
                  eps_trace_over_3_crust_mantle(i,j,k,ispec)
            enddo
          enddo
        enddo
      enddo
      if (it_of_this_exact_subset == it_exact_subset_end) then
        do it_of_this_exact_subset = 1, it_exact_subset_end
          write(IFILE_FOR_EXACT_UNDOING,rec=it_exact_subset_offset+it_of_this_exact_subset) &
            buffer_for_disk(:,it_of_this_exact_subset)
        enddo
        it_of_this_exact_subset = 1
        it_exact_subset_offset = it_exact_subset_offset + it_exact_subset_end
        it_exact_subset_end = min(NSTEP_FOR_EXACT_UNDOING, it_end - it_exact_subset_offset)
      else
        it_of_this_exact_subset = it_of_this_exact_subset + 1
      endif
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
        it_of_this_exact_subset = it_of_this_exact_subset + 1
        if (it_of_this_exact_subset > it_exact_subset_end) then
          it_exact_subset_offset = it_exact_subset_offset + it_exact_subset_end
          it_exact_subset_end = min(NSTEP_FOR_EXACT_UNDOING, it_end - it_exact_subset_offset)
          do it_of_this_exact_subset = 1, it_exact_subset_end
            ! here we time revert the forward run by reading time step NSTEP - it + 1
            ! but here, it == it_exact_subset_offset + it_of_this_exact_subset
            read(IFILE_FOR_EXACT_UNDOING,rec=NSTEP-it_exact_subset_offset-it_of_this_exact_subset+1) &
              buffer_for_disk(:,it_of_this_exact_subset)
          enddo
          it_of_this_exact_subset = 1
        endif

        do ispec = 1, NSPEC_CRUST_MANTLE
          do k = 1, NGLLZ
            do j = 1, NGLLY
              do i = 1, NGLLX
                iglob = ibool_crust_mantle(i,j,k,ispec)
                if (integer_mask_ibool_exact_undo(iglob) /= -1) then
                  b_eps_trace_over_3_crust_mantle(i,j,k,ispec) = &
                    buffer_for_disk(integer_mask_ibool_exact_undo(iglob),it_of_this_exact_subset)
                endif
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
    ! attention: for GPU_MODE and ANISOTROPIC_KL it is necessary to use resort_array (see lines 265-268)
    if (SIMULATION_TYPE == 3) then
      call compute_kernels()
    endif

    ! outputs movie files
    if (MOVIE_SURFACE .or. MOVIE_VOLUME) call write_movie_output()

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

 100 continue


  if (SIMULATION_TYPE == 3 .and. GPU_MODE) then
    ! attention: cijkl_kl_crust_mantle is sorted differently on GPU and CPU
    call resort_array(Mesh_pointer)
  endif

  ! full gravity
  if (SIMULATION_TYPE == 3 .and. FULL_GRAVITY) then
    ! calculate the gravity kernels (convolution) using SIEM
    call SIEM_compute_gravity_kernels()
  endif

  ! close the huge file that contains a dump of all the time steps to disk
  if (EXACT_UNDOING_TO_DISK) call finish_exact_undoing_to_disk()

  ! user output of runtime
  call print_elapsed_time()

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
  use specfem_par_noise

  implicit none

  ! to store forward wave fields
  if (SIMULATION_TYPE == 1) then
    if (SAVE_FORWARD .or. (NUMBER_OF_RUNS > 1 .and. NUMBER_OF_THIS_RUN < NUMBER_OF_RUNS)) then
      ! wavefield
      call transfer_fields_cm_from_device(NDIM*NGLOB_CRUST_MANTLE, &
                                          displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle,Mesh_pointer)
      call transfer_fields_ic_from_device(NDIM*NGLOB_INNER_CORE, &
                                          displ_inner_core,veloc_inner_core,accel_inner_core,Mesh_pointer)
      call transfer_fields_oc_from_device(NGLOB_OUTER_CORE, &
                                          displ_outer_core,veloc_outer_core,accel_outer_core,Mesh_pointer)
      ! strain
      call transfer_strain_cm_from_device(Mesh_pointer,eps_trace_over_3_crust_mantle, &
                                          epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle, &
                                          epsilondev_xy_crust_mantle,epsilondev_xz_crust_mantle, &
                                          epsilondev_yz_crust_mantle)
      call transfer_strain_ic_from_device(Mesh_pointer,eps_trace_over_3_inner_core, &
                                          epsilondev_xx_inner_core,epsilondev_yy_inner_core, &
                                          epsilondev_xy_inner_core,epsilondev_xz_inner_core, &
                                          epsilondev_yz_inner_core)
      ! rotation
      if (ROTATION_VAL) then
        call transfer_rotation_from_device(Mesh_pointer,A_array_rotation,B_array_rotation)
      endif

      ! attenuation memory variables
      if (ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL) then
        call transfer_rmemory_cm_from_device(Mesh_pointer,R_xx_crust_mantle,R_yy_crust_mantle,R_xy_crust_mantle, &
                                             R_xz_crust_mantle,R_yz_crust_mantle)
        call transfer_rmemory_ic_from_device(Mesh_pointer,R_xx_inner_core,R_yy_inner_core,R_xy_inner_core, &
                                             R_xz_inner_core,R_yz_inner_core)
      endif
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
      call transfer_kernels_noise_to_host(Mesh_pointer,sigma_kl_crust_mantle,NSPEC_CRUST_MANTLE)
    endif

    ! approximative Hessian for preconditioning kernels
    if (APPROXIMATE_HESS_KL) then
      call transfer_kernels_hess_cm_tohost(Mesh_pointer, &
                                           hess_kl_crust_mantle, &
                                           hess_rho_kl_crust_mantle, &
                                           hess_kappa_kl_crust_mantle, &
                                           hess_mu_kl_crust_mantle, &
                                           NSPEC_CRUST_MANTLE)
    endif

    ! outer core
    if (SAVE_KERNELS_OC) then
      call transfer_kernels_oc_to_host(Mesh_pointer, &
                                       rho_kl_outer_core, &
                                       alpha_kl_outer_core,NSPEC_OUTER_CORE)
    endif

    ! inner core
    if (SAVE_KERNELS_IC) then
      call transfer_kernels_ic_to_host(Mesh_pointer, &
                                       rho_kl_inner_core, &
                                       alpha_kl_inner_core, &
                                       beta_kl_inner_core,NSPEC_INNER_CORE)
    endif
  endif

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
    !if (myrank == 0 ) print *,"  VTK rendering..."

    ! updates time
    currenttime = sngl((it-1)*DT-t0)

    ! transfers fields from GPU to host
    if (GPU_MODE) then
      !if (myrank == 0 ) print *,"  VTK: transferring velocity from GPU"
      call transfer_veloc_cm_from_device(NDIM*NGLOB_CRUST_MANTLE,veloc_crust_mantle,Mesh_pointer)
    endif

    ! updates wavefield
    !if (myrank == 0 ) print *,"  VTK: it = ",it," out of ",it_end," - norm of velocity field"
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
        call gatherv_all_r(vtkdata,data_size, &
                            vtkdata_all,vtkdata_points_all,vtkdata_offset_all, &
                            vtkdata_numpoints_all,NPROCTOT_VAL)
        ! updates VTK window
        call visualize_vtkdata(it,currenttime,vtkdata_all)
      else
        ! all other process just send data
        call gatherv_all_r(vtkdata,data_size, &
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


!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_exact_undoing_to_disk()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  integer :: ispec,iglob,i,j,k
  integer :: counter,record_length
  real(kind=CUSTOM_REAL) :: radius
  character(len=MAX_STRING_LEN) :: outputname

  ! checks if anything to do
  if (.not. EXACT_UNDOING_TO_DISK) return

  ! checks flags
  if (GPU_MODE) call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK not supported for GPUs')

  if (UNDO_ATTENUATION) &
    call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK needs UNDO_ATTENUATION to be off because it computes the kernel directly instead')

  if (SIMULATION_TYPE == 1 .and. .not. SAVE_FORWARD) &
    call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK requires SAVE_FORWARD if SIMULATION_TYPE == 1')

  if (ANISOTROPIC_KL) call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK requires ANISOTROPIC_KL to be turned off')

  if (SIMULATION_TYPE /= 1 .and. SIMULATION_TYPE /= 3) &
    call exit_MPI(myrank,'EXACT_UNDOING_TO_DISK can only be used with SIMULATION_TYPE == 1 or SIMULATION_TYPE == 3')

  ! determine the largest value of iglob that we need to save to disk,
  ! since we save the upper part of the mesh only in the case of surface-wave kernels
  ! crust_mantle
  allocate(integer_mask_ibool_exact_undo(NGLOB_CRUST_MANTLE))
  integer_mask_ibool_exact_undo(:) = -1

  counter = 0
  do ispec = 1, NSPEC_CRUST_MANTLE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          ! xstore ystore zstore have previously been converted to r theta phi in rstore
          radius = rstore_crust_mantle(1,iglob) ! radius r (normalized)
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
  allocate(buffer_for_disk(counter,NSTEP_FOR_EXACT_UNDOING))

  ! open the file in which we will dump all the time steps (in a single file)
  write(outputname,"('huge_dumps/proc',i6.6,'_huge_dump_of_all_time_steps.bin')") myrank
  inquire(iolength=record_length) buffer_for_disk(:,1)
  ! we write to or read from the file depending on the simulation type
  if (SIMULATION_TYPE == 1) then
    open(file=outputname, unit=IFILE_FOR_EXACT_UNDOING, action='write', status='unknown', &
                    form='unformatted', access='direct', recl=record_length)
  else if (SIMULATION_TYPE == 3) then
    open(file=outputname, unit=IFILE_FOR_EXACT_UNDOING, action='read', status='old', &
                    form='unformatted', access='direct', recl=record_length)
  endif

  if (SIMULATION_TYPE == 1) then
    it_of_this_exact_subset = 1
    it_exact_subset_offset = it_begin - 1
    it_exact_subset_end = min(NSTEP_FOR_EXACT_UNDOING, it_end - it_begin + 1)
  else if (SIMULATION_TYPE == 3) then
    ! Trigger a read at the start of the loop
    it_of_this_exact_subset = 0
    it_exact_subset_offset = it_begin - 1
    it_exact_subset_end = 0
  endif

  end subroutine setup_exact_undoing_to_disk

!
!-------------------------------------------------------------------------------------------------
!

  subroutine finish_exact_undoing_to_disk()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! checks if anything to do
  if (.not. EXACT_UNDOING_TO_DISK) return

  ! frees memory
  deallocate(integer_mask_ibool_exact_undo)
  deallocate(buffer_for_disk)

  ! close the huge file that contains a dump of all the time steps to disk
  close(IFILE_FOR_EXACT_UNDOING)

  end subroutine finish_exact_undoing_to_disk

!
!-------------------------------------------------------------------------------------------------
!

  subroutine it_open_energy_curve_file()

  use specfem_par

  implicit none

  ! checks if anything to do
  if (.not. OUTPUT_ENERGY) return

  !----  create a Gnuplot script to display the energy curve in log scale
  if (myrank == 0) then
    open(unit=IOUT_ENERGY,file=trim(OUTPUT_FILES)//'plot_energy.gnu',status='unknown',action='write')
    write(IOUT_ENERGY,*) 'set term wxt'
    write(IOUT_ENERGY,*) '#set term postscript landscape color solid "Helvetica" 22'
    write(IOUT_ENERGY,*) '#set output "energy.ps"'
    write(IOUT_ENERGY,*) 'set logscale y'
    write(IOUT_ENERGY,*) 'set xlabel "Time step number"'
    write(IOUT_ENERGY,*) 'set ylabel "Energy (J)"'
    write(IOUT_ENERGY,'(a152)') '#plot "energy.dat" us 1:2 t "Kinetic Energy" w l lc 1, "energy.dat" us 1:3 &
                         &t "Potential Energy" w l lc 2, "energy.dat" us 1:4 t "Total Energy" w l lc 4'
    write(IOUT_ENERGY,*) '#pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) '#plot "energy.dat" us 1:2 t "Kinetic Energy" w l lc 1'
    write(IOUT_ENERGY,*) '#pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) '#plot "energy.dat" us 1:3 t "Potential Energy" w l lc 2'
    write(IOUT_ENERGY,*) '#pause -1 "Hit any key..."'
    write(IOUT_ENERGY,*) 'plot "energy.dat" us 1:4 t "Total Energy" w l lc 4'
    write(IOUT_ENERGY,*) 'pause -1 "Hit any key..."'
    close(IOUT_ENERGY)
  endif

  ! open the file in which we will store the energy curve
  if (myrank == 0) &
    open(unit=IOUT_ENERGY,file=trim(OUTPUT_FILES)//'energy.dat',status='unknown',action='write')

  end subroutine it_open_energy_curve_file
