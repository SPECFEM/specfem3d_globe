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

  subroutine iterate_time_undoatt()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_noise
  use specfem_par_movie

  implicit none

  ! local parameters
  integer :: it_temp,seismo_current_temp
  integer :: ier
  integer :: buffer_size, it_of_buffer, ntstep_kl
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_displ_cm_store_buffer
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_displ_ic_store_buffer
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_displ_oc_store_buffer,b_accel_oc_store_buffer
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: b_noise_surface_movie_buffer
  double precision :: sizeval
  ! timing
  double precision, external :: wtime

  ! number of buffered snapshot
  ntstep_kl = max(1, NTSTEP_BETWEEN_COMPUTE_KERNELS)
  buffer_size = ceiling(dble(NT_DUMP_ATTENUATION) / ntstep_kl)

  ! energy curve outputs
  if (OUTPUT_ENERGY) call it_open_energy_curve_file()

!
!   s t a r t   t i m e   i t e r a t i o n s
!

  ! checks
  if (.not. UNDO_ATTENUATION) return

  ! checks with undo_attenuation
  if (UNDO_ATTENUATION) then
    ! note: NSTEP must not be a multiple of NT_DUMP_ATTENUATION, but should be equal or larger
    ! makes sure buffer size is not too big for total time length
    if (NSTEP < NT_DUMP_ATTENUATION) then
      print *,'Error undoing attenuation: time steps ',NSTEP,' smaller than buffer size ',NT_DUMP_ATTENUATION
      print *,'Please recompile the solver with your updated parameter set in Par_file.'
      call exit_MPI(myrank,'Error undoing attenuation: number of time steps are too small, please increase record length!')
    endif
  endif

  ! number of time subsets for time loop
  if (NSTEP_STEADY_STATE > 0) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'steady state simulation:'
      write(IMAIN,*) '  total number of time steps                       = ',NSTEP
      write(IMAIN,*) '  number of steady state time steps                = ',NSTEP_STEADY_STATE
      write(IMAIN,*) '  number of transient time steps                   = ',NSTEP-NSTEP_STEADY_STATE
      write(IMAIN,*)
      call flush_IMAIN()
    endif
    ! subsets for steady state only
    NSUBSET_ITERATIONS = ceiling( dble(NSTEP_STEADY_STATE)/dble(NT_DUMP_ATTENUATION) )
  else
    NSUBSET_ITERATIONS = ceiling( dble(NSTEP)/dble(NT_DUMP_ATTENUATION) )
  endif

  ! checks
  if (NSUBSET_ITERATIONS <= 0) call exit_MPI(myrank,'Error invalid number of time subsets for undoing attenuation')

  ! user output
  if (SAVE_FORWARD .or. SIMULATION_TYPE == 3) then
    if (myrank == 0) then
      write(IMAIN,*) 'undoing attenuation:'
      write(IMAIN,*) '  total number of time subsets                     = ',NSUBSET_ITERATIONS
      write(IMAIN,*) '  wavefield snapshots at every NT_DUMP_ATTENUATION = ',NT_DUMP_ATTENUATION
      write(IMAIN,*) '  number of buffered snapshots at each subset      = ',buffer_size
      write(IMAIN,*)
      call flush_IMAIN()

      ! file size per snapshot for save/read routines
      write(IMAIN,*) '  estimated snapshot file storage:'
      ! see fields stored in save_forward_arrays_undoatt() routine in file save_forward_arrays.f90
      ! displ_crust_mantle + veloc_crust_mantle + accel_crust_mantle
      sizeval = 3.d0 * dble(NDIM) * dble(NGLOB_CRUST_MANTLE) * dble(CUSTOM_REAL)
      ! displ_inner_core + veloc_inner_core + accel_inner_core
      sizeval = sizeval + 3.d0 * dble(NDIM) * dble(NGLOB_INNER_CORE) * dble(CUSTOM_REAL)
      ! displ_outer_core + veloc_outer_core + accel_outer_core
      sizeval = sizeval + 3.d0 * dble(NGLOB_OUTER_CORE) * dble(CUSTOM_REAL)
      if (ROTATION_VAL) then
        ! A_array_rotation + B_array_rotation
        sizeval = sizeval + 2.d0 * dble(NGLLX*NGLLY*NGLLZ) * dble(NSPEC_OUTER_CORE_ROTATION) * dble(CUSTOM_REAL)
      endif
      if (ATTENUATION_VAL) then
        ! R_xx_crust_mantle + R_yy_crust_mantle + R_xy_crust_mantle + R_xz_crust_mantle + R_yz_crust_mantle
        sizeval = sizeval + 5.d0 * dble(NGLLX*NGLLY*NGLLZ*N_SLS) * dble(NSPEC_CRUST_MANTLE_ATTENUATION) * dble(CUSTOM_REAL)
        ! R_xx_inner_core + R_yy_inner_core + R_xy_inner_core + R_xz_inner_core + R_yz_inner_core
        sizeval = sizeval + 5.d0 * dble(NGLLX*NGLLY*NGLLZ*N_SLS) * dble(NSPEC_INNER_CORE_ATTENUATION) * dble(CUSTOM_REAL)
      endif
      ! in MB
      sizeval = sizeval / 1024.d0 / 1024.d0
      write(IMAIN,*) '  size of single time snapshot subset per slice         = ', sngl(sizeval),'MB'
      ! for all processes
      sizeval = sizeval * dble(NPROCTOT)
      write(IMAIN,*) '  size of single time snapshot subset for all processes = ', sngl(sizeval),'MB'
      write(IMAIN,*) '                                                        = ', sngl(sizeval/1024.d0),'GB'
      write(IMAIN,*)
      ! for all processes and all subsets
      sizeval = sizeval * dble(NSUBSET_ITERATIONS)
      write(IMAIN,*) '  total size of all time snapshots subsets              = ', sngl(sizeval),'MB'
      write(IMAIN,*) '                                                        = ', sngl(sizeval/1024.d0),'GB'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! allocates buffers
  if (SIMULATION_TYPE == 3) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  wavefield buffers:'
      ! crust/mantle
      ! buffer(NDIM,NGLOB_CRUST_MANTLE_ADJOINT,NT_DUMP_ATTENUATION) in MB
      sizeval = dble(NDIM) * dble(NGLOB_CRUST_MANTLE_ADJOINT) * dble(buffer_size) * dble(CUSTOM_REAL) / 1024.d0 / 1024.d0
      write(IMAIN,*) '  size of crust/mantle wavefield buffer per slice = ', sngl(sizeval),'MB'
      ! outer core
      ! buffer(NGLOB_OUTER_CORE_ADJOINT,NT_DUMP_ATTENUATION) in MB
      sizeval = 2.d0 * dble(NGLOB_OUTER_CORE_ADJOINT) * dble(buffer_size) * dble(CUSTOM_REAL) / 1024.d0 / 1024.d0
      write(IMAIN,*) '  size of   outer core wavefield buffer per slice = ', sngl(sizeval),'MB'
      ! inner core
      ! buffer(NDIM,NGLOB_CRUST_MANTLE_ADJOINT,NT_DUMP_ATTENUATION) in MB
      sizeval = dble(NDIM) * dble(NGLOB_INNER_CORE_ADJOINT) * dble(buffer_size) * dble(CUSTOM_REAL) / 1024.d0 / 1024.d0
      write(IMAIN,*) '  size of   inner core wavefield buffer per slice = ', sngl(sizeval),'MB'
      write(IMAIN,*)
      call flush_IMAIN()

      ! total size of buffers
      ! b_displ_cm_store_buffer
      sizeval = dble(NDIM) * dble(NGLOB_CRUST_MANTLE_ADJOINT) * dble(buffer_size) * dble(CUSTOM_REAL)
      ! b_displ_oc_store_buffer + b_accel_oc_store_buffer
      sizeval = sizeval + 2.d0 * dble(NGLOB_OUTER_CORE_ADJOINT) * dble(buffer_size) * dble(CUSTOM_REAL)
      ! b_displ_ic_store_buffer
      sizeval = sizeval + dble(NDIM) * dble(NGLOB_INNER_CORE_ADJOINT) * dble(buffer_size) * dble(CUSTOM_REAL)
      ! in MB
      sizeval = sizeval / 1024.d0 / 1024.d0
      write(IMAIN,*) '  total size of wavefield buffers per slice       = ', sngl(sizeval),'MB'
      write(IMAIN,*) '                                                  = ', sngl(sizeval/1024.d0),'GB'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    !! DK DK to Daniel, July 2013: in the case of GPU_MODE it would probably be better to leave these arrays on the host
    !! i.e. on the CPU, in order to be able to use all the (unused) memory of the host for them, since they are
    !! (purposely) huge and designed to use almost all the memory available (by carefully optimizing the
    !! value of NT_DUMP_ATTENUATION); when writing to these buffers, it will then be OK to use non-blocking writes
    !! from the device to the host, since we do not reuse the content of these buffers until much later, in a second part
    !! of the run; however when reading back from these buffers, the reads from host to device should then be blocking
    !! because we then use the values read immediately (one value at a time, but to reduce the total number of reads
    !! across the PCI-Express bus we could / should consider reading them 10 by 10 for instance (?) if that fits
    !! in the memory of the GPU
    allocate(b_displ_cm_store_buffer(NDIM,NGLOB_CRUST_MANTLE_ADJOINT,buffer_size),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating b_displ_cm_store_buffer')
    b_displ_cm_store_buffer(:,:,:) = 0.0_CUSTOM_REAL

    allocate(b_displ_oc_store_buffer(NGLOB_OUTER_CORE_ADJOINT,buffer_size),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating b_displ_oc_store_buffer')
    b_displ_oc_store_buffer(:,:) = 0.0_CUSTOM_REAL

    allocate(b_accel_oc_store_buffer(NGLOB_OUTER_CORE_ADJOINT,buffer_size),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating b_accel_oc_store_buffer')
    b_accel_oc_store_buffer(:,:) = 0.0_CUSTOM_REAL

    allocate(b_displ_ic_store_buffer(NDIM,NGLOB_INNER_CORE_ADJOINT,buffer_size),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating b_displ_ic_store_buffer')
    b_displ_ic_store_buffer(:,:,:) = 0.0_CUSTOM_REAL

    ! noise kernel for source strength (sigma_kernel) needs buffer for reconstructed noise_surface_movie array,
    ! otherwise we need file i/o which will considerably slow down performance
    if (NOISE_TOMOGRAPHY == 3) then
      allocate(b_noise_surface_movie_buffer(NDIM,NGLLX,NGLLY,NSPEC_TOP,buffer_size),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating b_noise_surface_movie_buffer')
      b_noise_surface_movie_buffer(:,:,:,:,:) = 0.0_CUSTOM_REAL
    endif
  endif

  ! shift execution of simultaneous events to avoid a high peak bandwidth for snapshots file I/O
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1) call prepare_simultaneous_event_execution_shift_undoatt()

  ! synchronize all processes to make sure everybody is ready to start time loop
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) 'All processes are synchronized before time loop (undoatt)'
    write(IMAIN,*)
    write(IMAIN,*) 'Starting time iteration loop (undoatt)...'
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! create an empty file to monitor the start of the simulation
  if (myrank == 0) then
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/starttimeloop.txt',status='unknown',action='write')
    write(IOUT,*) 'hello, starting time loop in iterate_time_undoatt() routine'
    close(IOUT)
  endif

  ! synchronizes GPU kernels
  if (GPU_MODE) call gpu_synchronize()

  ! initializes variables for writing seismograms
  seismo_offset = it_begin-1
  seismo_current = 0

  ! initializes time increments
  it = 0

  ! get MPI starting time
  time_start = wtime()

  ! *********************************************************
  ! ************* MAIN LOOP OVER THE TIME STEPS *************
  ! *********************************************************

  if (EXACT_UNDOING_TO_DISK) call setup_exact_undoing_to_disk()

  ! transient period simulation
  if (NSTEP_STEADY_STATE > 0) then
    do it = it_begin,NSTEP-NSTEP_STEADY_STATE

      ! simulation status output and stability check
      if (mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == it_begin + 4 .or. it == it_end) then
        call check_stability()
        if (I_am_running_on_a_slow_node) goto 200
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

      ! write the seismograms with time shift (GPU_MODE transfer included)
      call write_seismograms()

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
    enddo   ! end of main time loop of transient state

    200 continue

    it = it - 1
  endif

  ! loops over time subsets
  do iteration_on_subset = 1, NSUBSET_ITERATIONS

    ! wavefield storage
    if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
      ! saves forward wavefields
      call save_forward_arrays_undoatt()

    else if (SIMULATION_TYPE == 3) then
      ! reads in last stored forward wavefield
      call read_forward_arrays_undoatt()

      ! note: after reading the restart files of displacement back from disk, recompute the strain from displacement;
      !       this is better than storing the strain to disk as well, which would drastically increase I/O volume
      ! computes strain based on current backward/reconstructed wavefield
      if (COMPUTE_AND_STORE_STRAIN) call compute_strain_att_backward()
    endif

    ! time loop within this iteration subset
    select case (SIMULATION_TYPE)
    case (1, 2)
      ! forward and adjoint simulations

      ! increment end of this subset
      if (iteration_on_subset < NSUBSET_ITERATIONS) then
        ! takes full length of subset
        it_subset_end = NT_DUMP_ATTENUATION
      else
        ! loops over remaining steps in last subset
        if (NSTEP_STEADY_STATE > 0) then
          it_subset_end = NSTEP_STEADY_STATE - (iteration_on_subset-1)*NT_DUMP_ATTENUATION
        else
          it_subset_end = NSTEP - (iteration_on_subset-1)*NT_DUMP_ATTENUATION
        endif
      endif
      ! checks end index
      if (it_subset_end > NT_DUMP_ATTENUATION) &
        call exit_MPI(myrank,'Error invalid buffer index for undoing attenuation')

      ! subset loop
      do it_of_this_subset = 1, it_subset_end

        it = it + 1

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

        enddo ! istage

        ! write the seismograms with time shift
        call write_seismograms()

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

      enddo ! subset loop

    case (3)
      ! kernel simulations

      ! intermediate storage of it and seismo_current positions
      it_temp = it
      it_of_buffer = 0
      seismo_current_temp = seismo_current

      ! increment end of this subset
      if (iteration_on_subset == 1) then
        ! loops over remaining steps in last forward subset
        if (NSTEP_STEADY_STATE > 0) then
          it_subset_end = NSTEP_STEADY_STATE - (NSUBSET_ITERATIONS-1)*NT_DUMP_ATTENUATION
        else
          it_subset_end = NSTEP - (NSUBSET_ITERATIONS-1)*NT_DUMP_ATTENUATION
        endif
      else
        ! takes full length of subset
        it_subset_end = NT_DUMP_ATTENUATION
      endif
      ! checks end index
      if (it_subset_end > NT_DUMP_ATTENUATION) &
        call exit_MPI(myrank,'Error invalid buffer index for undoing attenuation')

      ! reconstructs forward wavefields based on last stored wavefield data

      ! note: we step forward in time here, starting from last snapshot.
      !       the newly computed, reconstructed forward wavefields (b_displ_..) get stored in buffers.

      ! subset loop
      do it_of_this_subset = 1, it_subset_end

        it = it + 1

        ! simulation status output and stability check
        if (mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == it_begin + 4 .or. it == it_end) then
          call check_stability_backward()
        endif

        do istage = 1, NSTAGE_TIME_SCHEME ! is equal to 1 if Newmark because only one stage then

          if (USE_LDDRK) then
            ! update displacement using Runge-Kutta time scheme
            call update_displ_lddrk_backward()
          else
            ! update displacement using Newmark time scheme
            call update_displ_Newmark_backward()
          endif

          ! update Poisson's load and solve Poisson's equations
          if (FULL_GRAVITY) call SIEM_solve_poisson_backward()

          ! acoustic solver for outer core
          ! (needs to be done first, before elastic one)
          call compute_forces_acoustic_backward()

          ! elastic solver for crust/mantle and inner core
          call compute_forces_viscoelastic_backward()

        enddo ! istage

        ! transfers wavefields from GPU to CPU for buffering
        if (mod(it_temp+it_subset_end-it_of_this_subset+1, ntstep_kl) == 0) then

          it_of_buffer = it_of_buffer + 1

          if (GPU_MODE) then
#if defined(USE_CUDA) || defined(USE_HIP) || defined(USE_OPENCL)
            if (it_of_buffer >= 2) then
              call unregister_host_array(b_displ_cm_store_buffer(:,:, it_of_buffer-1))
            endif
            call register_host_array(NDIM*NGLOB_CRUST_MANTLE_ADJOINT, b_displ_cm_store_buffer(:,:, it_of_buffer))
#endif
            ! daniel debug: check if these transfers could be made async to overlap
            call transfer_ofs_b_displ_cm_from_device(NDIM*NGLOB_CRUST_MANTLE_ADJOINT,it_of_buffer, &
                                                    b_displ_cm_store_buffer,Mesh_pointer)
            call transfer_ofs_b_displ_ic_from_device(NDIM*NGLOB_INNER_CORE_ADJOINT,it_of_buffer, &
                                                    b_displ_ic_store_buffer,Mesh_pointer)
            call transfer_ofs_b_displ_oc_from_device(NGLOB_OUTER_CORE_ADJOINT,it_of_buffer, &
                                                    b_displ_oc_store_buffer,Mesh_pointer)
            call transfer_ofs_b_accel_oc_from_device(NGLOB_OUTER_CORE_ADJOINT,it_of_buffer, &
                                                    b_accel_oc_store_buffer,Mesh_pointer)
          else
            ! stores wavefield in buffers
  ! only the displacement needs to be stored in memory buffers in order to compute the sensitivity kernels,
  ! not the memory variables R_ij, because the sensitivity kernel calculations only involve the displacement
  ! and the strain, not the stress, and the strain can be recomputed on the fly by computing the gradient
  ! of the displacement read back from the memory buffers (see also https://github.com/SPECFEM/specfem3d_globe/issues/194)
            b_displ_cm_store_buffer(:,:,it_of_buffer) = b_displ_crust_mantle(:,:)
            b_displ_oc_store_buffer(:,it_of_buffer) = b_displ_outer_core(:)
            b_accel_oc_store_buffer(:,it_of_buffer) = b_accel_outer_core(:)
            b_displ_ic_store_buffer(:,:,it_of_buffer) = b_displ_inner_core(:,:)
          endif

          ! for noise kernel
          if (NOISE_TOMOGRAPHY == 3) then
            b_noise_surface_movie_buffer(:,:,:,:,it_of_buffer) = noise_surface_movie(:,:,:,:)
          endif

        endif
      enddo ! subset loop

      ! resets current it and seismo_current positions
      it = it_temp
      seismo_current = seismo_current_temp

      ! computes strain based on current adjoint wavefield
      if (COMPUTE_AND_STORE_STRAIN) call compute_strain_att()

      ! adjoint wavefield simulation
      do it_of_this_subset = 1, it_subset_end

        it = it + 1

        if (mod(it, ntstep_kl) == 0) then

          ! reads backward/reconstructed wavefield from buffers
          ! note: uses wavefield at corresponding time (NSTEP - it + 1 ), i.e. we have now time-reversed wavefields
          ! crust/mantle
          ! transfers wavefields from CPU to GPU
          if (GPU_MODE) then
            ! daniel debug: check if these transfers could be made async to overlap
            call transfer_ofs_b_displ_cm_to_device(NDIM*NGLOB_CRUST_MANTLE_ADJOINT,it_of_buffer, &
                                                  b_displ_cm_store_buffer,Mesh_pointer)
            call transfer_ofs_b_displ_ic_to_device(NDIM*NGLOB_INNER_CORE_ADJOINT,it_of_buffer, &
                                                  b_displ_ic_store_buffer,Mesh_pointer)
            call transfer_ofs_b_displ_oc_to_device(NGLOB_OUTER_CORE_ADJOINT,it_of_buffer, &
                                                  b_displ_oc_store_buffer,Mesh_pointer)
            call transfer_ofs_b_accel_oc_to_device(NGLOB_OUTER_CORE_ADJOINT,it_of_buffer, &
                                                  b_accel_oc_store_buffer,Mesh_pointer)
          else
  ! only the displacement needs to be stored in memory buffers in order to compute the sensitivity kernels,
  ! not the memory variables R_ij, because the sensitivity kernel calculations only involve the displacement
  ! and the strain, not the stress, and the strain can be recomputed on the fly by computing the gradient
  ! of the displacement read back from the memory buffers (see also https://github.com/SPECFEM/specfem3d_globe/issues/194)
            b_displ_crust_mantle(:,:) = b_displ_cm_store_buffer(:,:,it_of_buffer)
            b_displ_outer_core(:) = b_displ_oc_store_buffer(:,it_of_buffer)
            b_accel_outer_core(:) = b_accel_oc_store_buffer(:,it_of_buffer)
            b_displ_inner_core(:,:) = b_displ_ic_store_buffer(:,:,it_of_buffer)
          endif

          ! for noise kernel
          if (NOISE_TOMOGRAPHY == 3) then
            noise_surface_movie(:,:,:,:) = b_noise_surface_movie_buffer(:,:,:,:,it_of_buffer)
          endif

          it_of_buffer = it_of_buffer - 1

        endif

        ! simulation status output and stability check
        if (mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == it_begin + 4 .or. it == it_end) then
          call check_stability()
          if (I_am_running_on_a_slow_node) goto 100
        endif

        ! computes adjoint wavefield
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

        enddo ! istage

        ! write the seismograms with time shift
        call write_seismograms()

        ! kernel computation
        ! adjoint simulations: kernels
        ! attention: for GPU_MODE and ANISOTROPIC_KL it is necessary to use resort_array (see lines 442-445)
        if (mod(it, ntstep_kl) == 0) then
#if defined(USE_CUDA) || defined(USE_HIP) || defined(USE_OPENCL)
          if (GPU_MODE) then
            call unregister_host_array(b_displ_cm_store_buffer(:,:, it_of_buffer+1))
            if (it_of_buffer > 0) then
              call register_host_array(NDIM*NGLOB_CRUST_MANTLE_ADJOINT, &
                                      b_displ_cm_store_buffer(:,:, it_of_buffer))
            endif
          endif
#endif
          call compute_kernels()
        endif

      enddo ! subset loop

    end select ! SIMULATION_TYPE

  !
  !---- end of time iteration loop
  !
  enddo   ! end of main time loop

 100 continue


  if (SIMULATION_TYPE == 3 .and. GPU_MODE) then
    ! attention: cijkl_kl_crust_mantle is sorted differently on GPU and CPU
    call resort_array(Mesh_pointer)
  endif

  ! frees undo_attenuation buffers
  if (SIMULATION_TYPE == 3) then
    deallocate(b_displ_cm_store_buffer, &
               b_displ_oc_store_buffer, &
               b_accel_oc_store_buffer, &
               b_displ_ic_store_buffer)
    ! noise simulations
    if (NOISE_TOMOGRAPHY == 3) then
      deallocate(b_noise_surface_movie_buffer)
    endif
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

  ! closes energy file
  if (OUTPUT_ENERGY .and. myrank == 0) close(IOUT_ENERGY)

  ! safety check of last time loop increment
  if (it /= it_end) then
    print *,'Error time increments: it_end = ',it_end,' and last it = ',it,' do not match!'
    call exit_MPI(myrank,'Error invalid time increment ending')
  endif

  end subroutine iterate_time_undoatt
