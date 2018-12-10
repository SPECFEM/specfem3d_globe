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

  subroutine setup_sources_receivers()

  use specfem_par
  implicit none

  ! locates sources and determines simulation start time t0
  call setup_sources()

  ! reads in stations file and locates receivers
  call setup_receivers()

  ! write source and receiver VTK files for Paraview
  call setup_sources_receivers_VTKfile()

  ! pre-compute source arrays
  call setup_sources_precompute_arrays()

  ! pre-compute receiver interpolation factors
  call setup_receivers_precompute_intp()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Total number of samples for seismograms = ',NSTEP
    write(IMAIN,*)
    if (NSOURCES > 1) write(IMAIN,*) 'Using ',NSOURCES,' point sources'
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! frees arrays
  deallocate(theta_source,phi_source)

  ! topography array no more needed
  if (TOPOGRAPHY) then
    if (allocated(ibathy_topo) ) deallocate(ibathy_topo)
  endif

  end subroutine setup_sources_receivers

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_movie
  implicit none

  ! local parameters
  double precision :: min_tshift_src_original
  integer :: isource,ier
  character(len=MAX_STRING_LEN) :: filename

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'sources:',NSOURCES
    call flush_IMAIN()
  endif

  ! allocate arrays for source
  allocate(islice_selected_source(NSOURCES), &
           ispec_selected_source(NSOURCES), &
           Mxx(NSOURCES), &
           Myy(NSOURCES), &
           Mzz(NSOURCES), &
           Mxy(NSOURCES), &
           Mxz(NSOURCES), &
           Myz(NSOURCES),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating source arrays')

  allocate(xi_source(NSOURCES), &
           eta_source(NSOURCES), &
           gamma_source(NSOURCES),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating source arrays')

  allocate(tshift_src(NSOURCES), &
           hdur(NSOURCES), &
           hdur_Gaussian(NSOURCES),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating source arrays')

  allocate(theta_source(NSOURCES), &
           phi_source(NSOURCES),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating source arrays')

  allocate(nu_source(NDIM,NDIM,NSOURCES),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating source arrays')

  if (USE_FORCE_POINT_SOURCE) then
    allocate(force_stf(NSOURCES),factor_force_source(NSOURCES), &
             comp_dir_vect_source_E(NSOURCES), &
             comp_dir_vect_source_N(NSOURCES), &
             comp_dir_vect_source_Z_UP(NSOURCES),stat=ier)
    if (ier /= 0) stop 'error allocating arrays for force point sources'
  endif

  ! sources
  ! BS BS moved open statement and writing of first lines into sr.vtk before the
  ! call to locate_sources, where further write statements to that file follow
  if (myrank == 0) then
  ! write source and receiver VTK files for Paraview
    filename = trim(OUTPUT_FILES)//'/sr_tmp.vtk'
    open(IOUT_VTK,file=trim(filename),status='unknown',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening temporary file sr_temp.vtk')
    write(IOUT_VTK,'(a)') '# vtk DataFile Version 2.0'
    write(IOUT_VTK,'(a)') 'Source and Receiver VTK file'
    write(IOUT_VTK,'(a)') 'ASCII'
    write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
    !  LQY -- won't be able to know NSOURCES+nrec at this point...
    write(IOUT_VTK, '(a,i6,a)') 'POINTS ', NSOURCES, ' float'
    ! closing file, rest of information will be appended later on
    close(IOUT_VTK)
  endif

  ! locate sources in the mesh
  call locate_sources(NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE,ibool_crust_mantle, &
                     xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                     ELLIPTICITY_VAL,min_tshift_src_original)

  ! determines onset time
  call setup_stf_constants(min_tshift_src_original)

  ! count number of sources located in this slice
  nsources_local = 0
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    do isource = 1,NSOURCES
      if (myrank == islice_selected_source(isource)) nsources_local = nsources_local + 1
    enddo
  endif

  ! determines number of times steps for simulation
  call setup_timesteps()

  ! prints source time functions and spectrum to output files
  if (PRINT_SOURCE_TIME_FUNCTION) call print_stf_file()

  ! get information about event name and location
  ! (e.g. needed for SAC seismograms)

  ! The following line is added for get_event_info subroutine.
  ! Because the way NSOURCES_SAC was declared has been changed.
  ! The rest of the changes in this program is just the updates of the subroutines that
  ! I did changes, e.g., adding/removing parameters. by Ebru Bozdag
  call get_event_info_parallel(yr_SAC,jda_SAC,mo_SAC, da_SAC, ho_SAC,mi_SAC,sec_SAC, &
                               event_name_SAC,t_cmt_SAC,t_shift_SAC, &
                               elat_SAC,elon_SAC,depth_SAC,mb_SAC,ms_SAC,cmt_lat_SAC, &
                               cmt_lon_SAC,cmt_depth_SAC,cmt_hdur_SAC,NSOURCES, &
                               Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)

  ! noise simulations ignore the CMTSOLUTIONS sources but employ a noise-spectrum source S_squared instead
  ! checks if anything to do for noise simulations
  if (NOISE_TOMOGRAPHY /= 0) then
    if (myrank == 0) then
      write(IMAIN,*) 'noise simulation will ignore CMT sources'
    endif
  endif

  end subroutine setup_sources

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_stf_constants(min_tshift_src_original)

  use specfem_par
  use specfem_par_movie
  implicit none

  double precision,intent(in) :: min_tshift_src_original

  ! local parameters
  integer :: isource

  ! makes smaller hdur for movies
  logical,parameter :: USE_SMALLER_HDUR_MOVIE = .true.

  if (abs(minval(tshift_src)) > TINYVAL) &
    call exit_MPI(myrank,'one tshift_src must be zero, others must be positive')

  ! filter source time function by Gaussian with hdur = HDUR_MOVIE when writing movies or shakemaps
  if (MOVIE_SURFACE .or. MOVIE_VOLUME) then
    ! smaller hdur_movie will do
    if (USE_SMALLER_HDUR_MOVIE) then
      ! hdur_movie gets assigned an automatic value based on the simulation resolution
      ! this will make that a bit smaller to have a higher-frequency movie output
      HDUR_MOVIE = 0.5 * HDUR_MOVIE
    endif

    ! new hdur for simulation
    hdur = sqrt(hdur**2 + HDUR_MOVIE**2)
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Each source is being convolved with HDUR_MOVIE = ',HDUR_MOVIE
      write(IMAIN,*)
    endif
  endif

  ! convert the half duration for triangle STF to the one for Gaussian STF
  hdur_Gaussian(:) = hdur(:)/SOURCE_DECAY_MIMIC_TRIANGLE

  ! define t0 as the earliest start time
  if (USE_FORCE_POINT_SOURCE) then
    ! point force sources
    ! (might start depending on the frequency given by hdur)
    ! note: point force sources will give the dominant frequency in hdur, thus the main period is 1/hdur.
    !       also, these sources might use a Ricker source time function instead of a Gaussian.
    !       For a Ricker source time function, a start time ~1.2 * main_period is a good choice.
    t0 = 0.d0
    do isource = 1,NSOURCES
      select case(force_stf(isource))
      case (0)
        ! Gaussian source time function
        t0 = min(t0,1.5d0 * (tshift_src(isource) - hdur(isource)))
      case (1)
        ! Ricker source time function
        t0 = min(t0,1.2d0 * (tshift_src(isource) - 1.0d0/hdur(isource)))
      case (2)
        ! Heaviside
        t0 = min(t0,1.5d0 * (tshift_src(isource) - hdur(isource)))
      case default
        stop 'unsupported force_stf value!'
      end select
    enddo
    ! start time defined as positive value, will be subtracted
    t0 = - t0
  else
    ! moment tensors
    ! (based on Heaviside functions)
    t0 = - 1.5d0 * minval( tshift_src(:) - hdur(:) )
  endif

  ! uses an external file for source time function, which starts at time 0.0
  if (EXTERNAL_SOURCE_TIME_FUNCTION) then
    hdur(:) = 0.d0
    t0      = 0.d0
  endif

  ! checks if user set USER_T0 to fix simulation start time
  ! note: USER_T0 has to be positive
  if (USER_T0 > 0.d0) then
    ! user cares about origin time and time shifts of the CMTSOLUTION
    ! and wants to fix simulation start time to a constant start time
    ! time 0 on time axis will correspond to given origin time

    ! notifies user
    if (myrank == 0) then
      write(IMAIN,*) 'USER_T0: ',USER_T0
      write(IMAIN,*) 't0: ',t0,'min_tshift_src_original: ',min_tshift_src_original
      write(IMAIN,*)
    endif

    ! checks if automatically set t0 is too small
    ! note: min_tshift_src_original can be a positive or negative time shift (minimum from all tshift)
    if (t0 <= USER_T0 + min_tshift_src_original) then
      ! by default, tshift_src(:) holds relative time shifts with a minimum time shift set to zero
      ! re-adds (minimum) original time shift such that sources will kick in
      ! according to their absolute time shift
      tshift_src(:) = tshift_src(:) + min_tshift_src_original

      ! sets new simulation start time such that
      ! simulation starts at t = - t0 = - USER_T0
      t0 = USER_T0

      ! notifies user
      if (myrank == 0) then
        write(IMAIN,*) '  set new simulation start time: ', - t0
        write(IMAIN,*)
      endif
    else
      ! start time needs to be at least t0 for numerical stability
      ! notifies user
      if (myrank == 0) then
        write(IMAIN,*) 'Error: USER_T0 is too small'
        write(IMAIN,*) '       must make one of three adjustments:'
        write(IMAIN,*) '       - increase USER_T0 to be at least: ',t0-min_tshift_src_original
        write(IMAIN,*) '       - decrease time shift in CMTSOLUTION file'
        write(IMAIN,*) '       - decrease hdur in CMTSOLUTION file'
        call flush_IMAIN()
      endif
      call exit_mpi(myrank,'Error USER_T0 is set but too small')
    endif
  else if (USER_T0 < 0.d0) then
    if (myrank == 0) then
      write(IMAIN,*) 'Error: USER_T0 is negative, must be set zero or positive!'
    endif
    call exit_mpi(myrank,'Error negative USER_T0 parameter in constants.h')
  endif

  end subroutine setup_stf_constants

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_timesteps()

  use specfem_par
  implicit none

  ! local parameters
  logical :: is_initial_guess

  ! checks if set by initial guess from read_compute_parameters() routine
  if (NTSTEP_BETWEEN_OUTPUT_SEISMOS == NSTEP) then
    is_initial_guess = .true.
  else
    is_initial_guess = .false.
  endif

  ! from initial guess in read_compute_parameters:
  !    compute total number of time steps, rounded to next multiple of 100
  !    NSTEP = 100 * (int(RECORD_LENGTH_IN_MINUTES * 60.d0 / (100.d0*DT)) + 1)
  !
  ! adds initial t0 time to update number of time steps and reach full record length
  if (abs(t0) > 0.d0) then
    ! note: for zero length, nstep has minimal of 5 timesteps for testing
    !       we won't extend this
    !
    ! careful: do not use RECORD_LENGTH_IN_MINUTES here, as it is only read by the master process
    !          when reading the parameter file, but it is not broadcasted to all other processes
    !          NSTEP gets broadcasted, so we work with this values
    if (NSTEP /= 5) then
      ! extend by bulk of 100 steps to account for half-duration rise time
      NSTEP = NSTEP + 100 * (int( abs(t0) / (100.d0*DT)) + 1)
    endif
  endif

  ! if doing benchmark runs to measure scaling of the code for a limited number of time steps only
  if (DO_BENCHMARK_RUN_ONLY) NSTEP = NSTEP_FOR_BENCHMARK

  ! checks length for symmetry in case of noise simulations
  if (NOISE_TOMOGRAPHY /= 0) then
    if (mod(NSTEP+1,2) /= 0) then
      print *,'Error noise simulation: invalid time steps = ',NSTEP,', NSTEP + 1 must be a multiple of 2 due to branch symmetry'
      call exit_MPI(myrank,'Error noise simulation: number of timesteps must be symmetric, due to +/- branches')
    endif
  endif

  ! time loop increments end
  it_end = NSTEP

  ! subsets used to save seismograms must not be larger than the whole time series,
  ! otherwise we waste memory
  if (NTSTEP_BETWEEN_OUTPUT_SEISMOS > NSTEP .or. is_initial_guess) NTSTEP_BETWEEN_OUTPUT_SEISMOS = NSTEP

  ! subsets used to save adjoint sources must not be larger than the whole time series,
  ! otherwise we waste memory
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    if (NTSTEP_BETWEEN_READ_ADJSRC > NSTEP) NTSTEP_BETWEEN_READ_ADJSRC = NSTEP
  endif

  ! buffering with undo_attenuation
  NT_DUMP_ATTENUATION = NT_DUMP_ATTENUATION_VAL
  if (UNDO_ATTENUATION) then
    ! makes sure buffer size is not too big for total time length
    !
    ! note: NSTEP must not be a multiple of NT_DUMP_ATTENUATION.
    !       the value from the header file NT_DUMP_ATTENUATION_VAL gives the optimal (maximum) number of time steps for buffering
    if (NSTEP < NT_DUMP_ATTENUATION) NT_DUMP_ATTENUATION = NSTEP
  endif

  ! debug
  !if (myrank == 0 ) print *,'setup time steps = ',NSTEP,' t0 = ',t0,' DT = ',DT

  end subroutine setup_timesteps

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_receivers()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

! check for imbalance of distribution of receivers or of adjoint sources
  logical, parameter :: CHECK_FOR_IMBALANCE = .false.

  ! local parameters
  integer :: irec,isource,nrec_tot_found,i
  integer :: nrec_simulation
  integer :: nadj_files_found,nadj_files_found_tot
  integer :: ier
  integer,dimension(0:NPROCTOT_VAL-1) :: tmp_rec_local_all
  integer :: maxrec,maxproc(1)
  double precision :: sizeval

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'receivers:'
    call flush_IMAIN()
  endif

  ! allocate memory for receiver arrays
  allocate(islice_selected_rec(nrec), &
           ispec_selected_rec(nrec), &
           xi_receiver(nrec), &
           eta_receiver(nrec), &
           gamma_receiver(nrec),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating receiver arrays')

  allocate(station_name(nrec), &
           network_name(nrec), &
           stlat(nrec), &
           stlon(nrec), &
           stele(nrec), &
           stbur(nrec),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating receiver arrays')

  allocate(nu(NDIM,NDIM,nrec),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating receiver arrays')

  !  receivers
  if (myrank == 0) then
    write(IMAIN,*)
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      write(IMAIN,*) 'Total number of receivers = ', nrec
    else
      write(IMAIN,*) 'Total number of adjoint sources = ', nrec
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! locate receivers in the crust in the mesh
  call locate_receivers(NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE,ibool_crust_mantle, &
                        xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                        yr_SAC,jda_SAC,ho_SAC,mi_SAC,sec_SAC, &
                        theta_source(1),phi_source(1) )

  ! count number of receivers located in this slice
  nrec_local = 0
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    ! note: for 1-chunk simulations, nrec is now the actual number of receivers found in this chunk
    !       (excludes stations located outside of chunk)
    nrec_simulation = nrec
    do irec = 1,nrec
      if (myrank == islice_selected_rec(irec)) nrec_local = nrec_local + 1
    enddo
  else
    nrec_simulation = NSOURCES
    do isource = 1, NSOURCES
      if (myrank == islice_selected_source(isource)) nrec_local = nrec_local + 1
    enddo
  endif

  ! counter for adjoint receiver stations in local slice, used to allocate adjoint source arrays
  nadj_rec_local = 0

  ! counts receivers for adjoint simulations
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! temporary counter to check if any files are found at all
    nadj_files_found = 0
    do irec = 1,nrec
      ! checks if slice is valid
      if (islice_selected_rec(irec) < 0 .or. islice_selected_rec(irec) > NPROCTOT_VAL-1) &
        call exit_MPI(myrank,'something is wrong with the source slice number in adjoint simulation')

      ! adjoint receiver station in this process slice
      if (myrank == islice_selected_rec(irec)) then
        ! updates counter
        nadj_rec_local = nadj_rec_local + 1

        ! checks **net**.**sta**.**MX**.adj files for correct number of time steps
        if (READ_ADJSRC_ASDF) then
          call check_adjoint_sources_asdf(irec,nadj_files_found)
        else
          call check_adjoint_sources(irec,nadj_files_found)
        endif
      endif
    enddo

    ! checks if any adjoint source files found at all
    call sum_all_i(nadj_files_found,nadj_files_found_tot)
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '    ',nadj_files_found_tot,' adjoint component traces found in all slices'
      if (nadj_files_found_tot == 0) &
        call exit_MPI(myrank,'no adjoint traces found, please check adjoint sources in directory SEM/')
    endif
  endif

  ! check that the sum of the number of receivers in each slice is nrec
  call sum_all_i(nrec_local,nrec_tot_found)
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'found a total of ',nrec_tot_found,' receivers in all slices'
    ! checks for number of receivers
    ! note: for 1-chunk simulations, nrec_simulations is the number of receivers/sources found in this chunk
    if (nrec_tot_found /= nrec_simulation) then
      call exit_MPI(myrank,'problem when dispatching the receivers')
    else
      write(IMAIN,*) 'this total is okay'
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! check that the sum of the number of receivers in each slice is nrec
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    if (myrank == 0 .and. nrec_tot_found /= nrec) &
      call exit_MPI(myrank,'total number of receivers is incorrect')
  endif
  call synchronize_all()

  ! statistics about allocation memory for seismograms & adj_sourcearrays
  ! user output info
  ! sources
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    ! user output
    if (myrank == 0) then
      ! note: all process allocate the full sourcearrays array
      ! sourcearrays(NDIM,NGLLX,NGLLY,NGLLZ,NSOURCES)
      sizeval = dble(NSOURCES) * dble(NDIM * NGLLX * NGLLY * NGLLZ * CUSTOM_REAL / 1024. / 1024. )
      ! outputs info
      write(IMAIN,*) 'source arrays:'
      write(IMAIN,*) '  number of sources is ',NSOURCES
      write(IMAIN,*) '  size of source array                 = ', sngl(sizeval),'MB'
      write(IMAIN,*) '                                       = ', sngl(sizeval/1024.d0),'GB'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! seismograms
  ! gather from slaves on master
  tmp_rec_local_all(:) = 0
  tmp_rec_local_all(0) = nrec_local
  if (NPROCTOT_VAL > 1) then
    call gather_all_singlei(nrec_local,tmp_rec_local_all,NPROCTOT_VAL)
  endif
  ! user output
  if (myrank == 0) then
    ! determines maximum number of local receivers and corresponding rank
    maxrec = maxval(tmp_rec_local_all(:))
    ! note: MAXLOC will determine the lower bound index as '1'.
    maxproc = maxloc(tmp_rec_local_all(:)) - 1
    ! seismograms array size in MB
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      ! seismograms need seismograms(NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS)
      sizeval = dble(maxrec) * dble(NDIM * NTSTEP_BETWEEN_OUTPUT_SEISMOS * CUSTOM_REAL / 1024. / 1024. )
    else
      ! adjoint seismograms need seismograms(NDIM*NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS)
      sizeval = dble(maxrec) * dble(NDIM * NDIM * NTSTEP_BETWEEN_OUTPUT_SEISMOS * CUSTOM_REAL / 1024. / 1024. )
    endif
    ! outputs info
    write(IMAIN,*) 'seismograms:'
    if (WRITE_SEISMOGRAMS_BY_MASTER) then
      write(IMAIN,*) '  seismograms written by master process only'
    else
      write(IMAIN,*) '  seismograms written by all processes'
    endif
    write(IMAIN,*) '  writing out seismograms at every NTSTEP_BETWEEN_OUTPUT_SEISMOS = ',NTSTEP_BETWEEN_OUTPUT_SEISMOS
    write(IMAIN,*) '  maximum number of local receivers is ',maxrec,' in slice ',maxproc(1)
    write(IMAIN,*) '  size of maximum seismogram array       = ', sngl(sizeval),'MB'
    write(IMAIN,*) '                                         = ', sngl(sizeval/1024.d0),'GB'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! adjoint sources
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! gather from slaves on master
    tmp_rec_local_all(:) = 0
    tmp_rec_local_all(0) = nadj_rec_local
    if (NPROCTOT_VAL > 1) then
      call gather_all_singlei(nadj_rec_local,tmp_rec_local_all,NPROCTOT_VAL)
    endif
    ! user output
    if (myrank == 0) then
      ! determines maximum number of local receivers and corresponding rank
      maxrec = maxval(tmp_rec_local_all(:))
      ! note: MAXLOC will determine the lower bound index as '1'.
      maxproc = maxloc(tmp_rec_local_all(:)) - 1
      !do i = 1, NPROCTOT_VAL
      !  if (tmp_rec_local_all(i) > maxrec) then
      !    maxrec = tmp_rec_local_all(i)
      !    maxproc = i-1
      !  endif
      !enddo
      ! source_adjoint size in MB
      ! source_adjoint(NDIM,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC)
      sizeval = dble(maxrec) * dble(NDIM * NTSTEP_BETWEEN_READ_ADJSRC * CUSTOM_REAL / 1024. / 1024. )
      ! note: in case IO_ASYNC_COPY is set, and depending of NSTEP_SUB_ADJ,
      !       this memory requirement might double.
      !       at this point, NSTEP_SUB_ADJ is not set yet...
      ! outputs info
      write(IMAIN,*) 'adjoint source arrays:'
      write(IMAIN,*) '  reading adjoint sources at every NTSTEP_BETWEEN_READ_ADJSRC = ',NTSTEP_BETWEEN_READ_ADJSRC
      if (IO_ASYNC_COPY) then
        write(IMAIN,*) '  using asynchronous buffer for file I/O of adjoint sources'
      endif
      write(IMAIN,*) '  maximum number of local adjoint sources is ',maxrec,' in slice ',maxproc(1)
      write(IMAIN,*) '  size of maximum adjoint source array = ', sngl(sizeval),'MB'
      write(IMAIN,*) '                                       = ', sngl(sizeval/1024.d0),'GB'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

! check for imbalance of distribution of receivers or of adjoint sources
  if (CHECK_FOR_IMBALANCE .and. NPROCTOT_VAL > 1) then
    call gather_all_singlei(nrec_local,tmp_rec_local_all,NPROCTOT_VAL)
    if (myrank == 0) then
      open(unit=9977,file='imbalance_of_nrec_local.dat',status='unknown')
      do i = 0,NPROCTOT_VAL-1
        write(9977,*) i,tmp_rec_local_all(i)
      enddo
      close(9977)
    endif

    call gather_all_singlei(nadj_rec_local,tmp_rec_local_all,NPROCTOT_VAL)
    if (myrank == 0) then
      open(unit=9977,file='imbalance_of_nadj_rec_local.dat',status='unknown')
      do i = 0,NPROCTOT_VAL-1
        write(9977,*) i,tmp_rec_local_all(i)
      enddo
      close(9977)
    endif

    if (myrank == 0) then
      open(unit=9977,file='plot_imbalance_histogram.gnu',status='unknown')
      write(9977,*) '#set terminal x11'
      write(9977,*) 'set terminal wxt'
      write(9977,*) '#set terminal gif'
      write(9977,*) '#set output "imbalance_histogram.gif"'
      write(9977,*)
      write(9977,*) 'set xrange [1:',NPROCTOT_VAL,']'
      write(9977,*) '#set xtics 0,0.1,1'
      write(9977,*) 'set boxwidth 1.'
      write(9977,*) 'set xlabel "Mesh slice number"'
      write(9977,*)
      write(9977,*) 'set ylabel "Number of receivers in that mesh slice"'
      write(9977,*) 'plot "imbalance_of_nrec_local.dat" with boxes'
      write(9977,*) 'pause -1 "hit any key..."'
      write(9977,*)
      write(9977,*) 'set ylabel "Number of adjoint sources in that mesh slice"'
      write(9977,*) 'plot "imbalance_of_nadj_rec_local.dat" with boxes'
      write(9977,*) 'pause -1 "hit any key..."'
      close(9977)
    endif

    call synchronize_all()
  endif

  end subroutine setup_receivers

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources_receivers_VTKfile()

  use specfem_par, only: myrank,OUTPUT_FILES,NSOURCES,nrec,MAX_STRING_LEN
  implicit none

  ! local parameters
  character(len=MAX_STRING_LEN) :: filename,filename_new
  character(len=MAX_STRING_LEN) :: command

  ! user output
  if (myrank == 0) then

    ! finishes VTK file
    !  we should know NSOURCES+nrec at this point...
    ! creates source/receiver location file
    filename = trim(OUTPUT_FILES)//'/sr_tmp.vtk'
    filename_new = trim(OUTPUT_FILES)//'/sr.vtk'
    write(command, &
  "('sed -e ',a1,'s/POINTS.*/POINTS',i6,' float/',a1,'<',a,'>',a)")&
      "'",NSOURCES + nrec,"'",trim(filename),trim(filename_new)

    ! note: this system() routine is non-standard Fortran
    call system_command(command)

    ! only extract receiver locations and remove temporary file
    filename_new = trim(OUTPUT_FILES)//'/receiver.vtk'
    write(command, &
  "('awk ',a1,'{if (NR < 5) print $0;if (NR == 6)&
   &print ',a1,'POINTS',i6,' float',a1,';if (NR > 5+',i6,')print $0}',a1,'<',a,'>',a)")&
      "'",'"',nrec,'"',NSOURCES,"'",trim(filename),trim(filename_new)

    ! note: this system() routine is non-standard Fortran
    call system_command(command)

    ! only extract source locations and remove temporary file
    filename_new = trim(OUTPUT_FILES)//'/source.vtk'
    write(command, &
  "('awk ',a1,'{if (NR < 6 + ',i6,') print $0}END{print}',a1,'<',a,'>',a,'; rm -f ',a)")&
      "'",NSOURCES,"'",trim(filename),trim(filename_new),trim(filename)

    ! note: this system() routine is non-standard Fortran
    call system_command(command)

  endif

  end subroutine setup_sources_receivers_VTKfile

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources_precompute_arrays()

  use specfem_par
  use specfem_par_crustmantle
  implicit none

  ! local parameters
  integer :: ier
  integer(kind=8) :: arraysize

  ! allocates source arrays
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    ! source interpolated on all GLL points in source element
    allocate(sourcearrays(NDIM,NGLLX,NGLLY,NGLLZ,NSOURCES),stat=ier)
    if (ier /= 0 ) then
      print *,'Error rank ',myrank,': allocating sourcearrays failed! number of sources = ',NSOURCES
      call exit_MPI(myrank,'Error allocating sourcearrays')
    endif
    ! initializes
    sourcearrays(:,:,:,:,:) = 0._CUSTOM_REAL

    ! stores source arrays
    call setup_sources_receivers_srcarr()
  endif

  ! adjoint source arrays
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! initializes adjoint source buffer
    ! reverse indexing
    allocate(iadj_vec(NSTEP),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating iadj_vec')
    ! initializes iadj_vec
    do it = 1,NSTEP
       iadj_vec(it) = NSTEP-it+1  ! default is for reversing entire record, e.g. 3000,2999,..,1
    enddo

    ! total number of adjoint source blocks to read in
    NSTEP_SUB_ADJ = ceiling( dble(NSTEP)/dble(NTSTEP_BETWEEN_READ_ADJSRC) )

    if (nadj_rec_local > 0) then
      allocate(source_adjoint(NDIM,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC), &
               stat=ier)
      if (ier /= 0 ) then
        print *,'Error rank ',myrank,': allocating source_adjoint failed! Please check your memory usage...'
        print *,'  failed number of local adjoint sources = ',nadj_rec_local,' steps = ',NTSTEP_BETWEEN_READ_ADJSRC
        call exit_MPI(myrank,'Error allocating adjoint sourcearrays')
      endif

      ! additional buffer for asynchronous file i/o
      if (IO_ASYNC_COPY .and. NSTEP_SUB_ADJ > 1) then
        ! allocates read buffer
        allocate(buffer_source_adjoint(NDIM,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC), &
                 stat=ier)
        if (ier /= 0 ) call exit_MPI(myrank,'Error allocating array buffer_source_adjoint')

        ! array size in bytes (note: the multiplication is split into two line to avoid integer-overflow)
        arraysize = NDIM *  CUSTOM_REAL
        arraysize = arraysize * nadj_rec_local * NTSTEP_BETWEEN_READ_ADJSRC

        ! debug
        !print *,'buffer_sourcearrays: size = ',arraysize,' Bytes = ',arraysize/1024./1024.,'MB'

        ! initializes io thread
        call prepare_adj_io_thread(buffer_source_adjoint,arraysize,nadj_rec_local)
      endif

      ! allocate indexing arrays
      allocate(iadjsrc(NSTEP_SUB_ADJ,2), &
               iadjsrc_len(NSTEP_SUB_ADJ),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating adjoint indexing arrays')

      ! initializes iadjsrc, iadjsrc_len and iadj_vec
      call setup_sources_receivers_adjindx(NSTEP,NSTEP_SUB_ADJ, &
                                           NTSTEP_BETWEEN_READ_ADJSRC, &
                                           iadjsrc,iadjsrc_len,iadj_vec)
    endif
  endif

  end subroutine setup_sources_precompute_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources_receivers_srcarr()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  integer :: isource,i,j,k,ispec !,iglob

  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearray
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd

  double precision :: xi,eta,gamma
  double precision :: hlagrange
  double precision :: norm

  do isource = 1,NSOURCES

    ! initializes
    sourcearray(:,:,:,:) = 0._CUSTOM_REAL

    !   check that the source slice number is okay
    if (islice_selected_source(isource) < 0 .or. islice_selected_source(isource) > NPROCTOT_VAL-1) &
      call exit_MPI(myrank,'Error: source slice number invalid')

    !   compute source arrays in source slice
    if (myrank == islice_selected_source(isource)) then

      ! element id which holds source
      ispec = ispec_selected_source(isource)

      ! checks bounds
      if (ispec < 1 .or. ispec > NSPEC_CRUST_MANTLE ) &
        call exit_MPI(myrank,'Error: source ispec number invalid')

      ! gets source location
      xi = xi_source(isource)
      eta = eta_source(isource)
      gamma = gamma_source(isource)

!      ! pre-computes source contribution on GLL points
!      call compute_arrays_source(sourcearray,xi,eta,gamma, &
!                          Mxx(isource),Myy(isource),Mzz(isource),Mxy(isource),Mxz(isource),Myz(isource), &
!                          xix_crust_mantle(:,:,:,ispec),xiy_crust_mantle(:,:,:,ispec),xiz_crust_mantle(:,:,:,ispec), &
!                          etax_crust_mantle(:,:,:,ispec),etay_crust_mantle(:,:,:,ispec),etaz_crust_mantle(:,:,:,ispec), &
!                          gammax_crust_mantle(:,:,:,ispec),gammay_crust_mantle(:,:,:,ispec),gammaz_crust_mantle(:,:,:,ispec), &
!                          xigll,yigll,zigll)
!
!      ! point forces, initializes sourcearray, used for simplified CUDA routines
!    !-------------POINT FORCE-----------------------------------------------
!      if (USE_FORCE_POINT_SOURCE) then
!        ! note: for use_force_point_source xi/eta/gamma are in the range [1,NGLL*]
!        iglob = ibool_crust_mantle(nint(xi),nint(eta),nint(gamma),ispec)
!
!        ! sets sourcearrays
!        do k = 1,NGLLZ
!          do j = 1,NGLLY
!            do i = 1,NGLLX
!              if (ibool_crust_mantle(i,j,k,ispec) == iglob) then
!                ! elastic source components
!                sourcearray(:,i,j,k) = nu_source(COMPONENT_FORCE_SOURCE,:,isource)
!              endif
!            enddo
!          enddo
!        enddo
!      endif
!    !-------------POINT FORCE-----------------------------------------------
!
!      ! stores source excitations
!      sourcearrays(:,:,:,:,isource) = sourcearray(:,:,:,:)
!    endif

      ! compute Lagrange polynomials at the source location
      call lagrange_any(xi,NGLLX,xigll,hxis,hpxis)
      call lagrange_any(eta,NGLLY,yigll,hetas,hpetas)
      call lagrange_any(gamma,NGLLZ,zigll,hgammas,hpgammas)

      if (USE_FORCE_POINT_SOURCE) then ! use of FORCESOLUTION files

        ! note: for use_force_point_source xi/eta/gamma are also in the range [-1,1], for exact positioning

        ! initializes source array
        sourcearrayd(:,:,:,:) = 0.0d0

        ! calculates source array for interpolated location
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX
              hlagrange = hxis(i) * hetas(j) * hgammas(k)

              ! elastic source
              norm = sqrt( comp_dir_vect_source_E(isource)**2 &
                         + comp_dir_vect_source_N(isource)**2 &
                         + comp_dir_vect_source_Z_UP(isource)**2 )

              ! checks norm of component vector
              if (norm < TINYVAL) then
                call exit_MPI(myrank,'error force point source: component vector has (almost) zero norm')
              endif

              ! normalizes vector
              comp_dir_vect_source_E(isource) = comp_dir_vect_source_E(isource) / norm
              comp_dir_vect_source_N(isource) = comp_dir_vect_source_N(isource) / norm
              comp_dir_vect_source_Z_UP(isource) = comp_dir_vect_source_Z_UP(isource) / norm

              ! we use a tilted force defined by its magnitude and the projections
              ! of an arbitrary (non-unitary) direction vector on the E/N/Z_UP basis
              !
              ! note: nu_source(iorientation,:,isource) is the rotation matrix from ECEF to local N-E-UP
              !       (defined in src/specfem3D/locate_sources.f90)
              sourcearrayd(:,i,j,k) = factor_force_source(isource) * hlagrange * &
                                      ( nu_source(1,:,isource) * comp_dir_vect_source_N(isource) + &
                                        nu_source(2,:,isource) * comp_dir_vect_source_E(isource) + &
                                        nu_source(3,:,isource) * comp_dir_vect_source_Z_UP(isource) )
            enddo
          enddo
        enddo

        ! distinguish between single and double precision for reals
        sourcearray(:,:,:,:) = real(sourcearrayd(:,:,:,:),kind=CUSTOM_REAL)

      else ! use of CMTSOLUTION files

        call compute_arrays_source(sourcearray,xi,eta,gamma, &
                          Mxx(isource),Myy(isource),Mzz(isource),Mxy(isource), &
                          Mxz(isource),Myz(isource), &
                          xix_crust_mantle(:,:,:,ispec), &
                          xiy_crust_mantle(:,:,:,ispec), &
                          xiz_crust_mantle(:,:,:,ispec), &
                          etax_crust_mantle(:,:,:,ispec), &
                          etay_crust_mantle(:,:,:,ispec), &
                          etaz_crust_mantle(:,:,:,ispec), &
                          gammax_crust_mantle(:,:,:,ispec), &
                          gammay_crust_mantle(:,:,:,ispec), &
                          gammaz_crust_mantle(:,:,:,ispec), &
                          xigll,yigll,zigll)

      endif

      ! stores source excitations
      sourcearrays(:,:,:,:,isource) = sourcearray(:,:,:,:)

    endif
  enddo

  end subroutine setup_sources_receivers_srcarr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources_receivers_adjindx(NSTEP,NSTEP_SUB_ADJ, &
                                             NTSTEP_BETWEEN_READ_ADJSRC, &
                                             iadjsrc,iadjsrc_len,iadj_vec)

  use constants

  implicit none

  integer NSTEP,NSTEP_SUB_ADJ,NTSTEP_BETWEEN_READ_ADJSRC

  integer, dimension(NSTEP_SUB_ADJ,2) :: iadjsrc ! to read input in chunks
  integer, dimension(NSTEP_SUB_ADJ) :: iadjsrc_len
  integer, dimension(NSTEP) :: iadj_vec

  ! local parameters
  integer :: iadj_block,it,it_sub_adj
  integer :: istart,iend

  ! initializes
  iadjsrc(:,:) = 0
  iadjsrc_len(:) = 0

  ! setting up chunks of NTSTEP_BETWEEN_READ_ADJSRC to read adjoint source traces
  ! i.e. as an example: total length NSTEP = 3000, chunk length NTSTEP_BETWEEN_READ_ADJSRC= 1000
  !                                then it will set first block from 2001 to 3000,
  !                                second block from 1001 to 2000 and so on...
  !
  ! see routine: compute_arrays_source_adjoint()
  !                     how we read in the adjoint source trace in blocks/chunk sizes
  !
  ! see routine: compute_add_sources_adjoint()
  !                     how the adjoint source is added to the (adjoint) acceleration field
  !counts blocks
  ! block number
  ! e.g. increases from 1 (case it=1-1000), 2 (case it=1001-2000) to 3 (case it=2001-3000)
  it_sub_adj = 0
  iadj_block = 1
  do it = 1,NSTEP
    ! we are at the edge of a block
    if (mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC) == 0) then
      ! sets it_sub_adj subset number
      it_sub_adj = iadj_block

      ! block start time ( e.g. 2001)
      istart = NSTEP-it_sub_adj*NTSTEP_BETWEEN_READ_ADJSRC+1
      ! final adj src array
      ! e.g. will be from 1000 to 1, but doesn't go below 1 in cases where NSTEP isn't
      ! a multiple of NTSTEP_BETWEEN_READ_ADJSRC
      if (istart < 1 ) istart = 1

      ! block end time (e.g. 3000)
      iend = NSTEP-(it_sub_adj-1)*NTSTEP_BETWEEN_READ_ADJSRC

      iadjsrc(iadj_block,1) = istart
      iadjsrc(iadj_block,2) = iend

      ! actual block length
      iadjsrc_len(iadj_block) = iend - istart + 1

      ! increases block number
      iadj_block = iadj_block + 1
    endif

    ! time stepping for adjoint sources:
    ! adjoint time step that corresponds to time step in simulation (it).
    ! note, that adjoint source has to be time-reversed with respect to the forward wavefield
    ! e.g.: first block 1 has iadjsrc_len = 1000 with start at 2001 and end at 3000
    !         so iadj_vec(1) = 1000 - 0, iadj_vec(2) = 1000 - 1, ..., to iadj_vec(1000) = 1000 - 999 = 1
    !         then for block 2, iadjsrc_len = 1000 with start at 1001 and end at 2000
    !         so iadj_vec(1001) = 1000 - 0, iadj_vec(1002) = 1000 - 1, .. and so on again down to 1
    !         then block 3 and your guess is right now... iadj_vec(2001) to iadj_vec(3000) is 1000 down to 1. :)
    iadj_vec(it) = iadjsrc_len(it_sub_adj) - mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC)

    ! checks that index is non-negative
    if (iadj_vec(it) < 1 ) iadj_vec(it) = 1
  enddo

  end subroutine setup_sources_receivers_adjindx

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_receivers_precompute_intp()

  use specfem_par
  implicit none

  ! local parameters
  integer :: ier
  integer :: nadj_hprec_local

  ! define local to global receiver numbering mapping
  ! needs to be allocated for subroutine calls (even if nrec_local == 0)
  allocate(number_receiver_global(nrec_local),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating global receiver numbering')

  ! allocates receiver interpolators
  if (nrec_local > 0) then
    ! allocates Lagrange interpolators for receivers
    allocate(hxir_store(nrec_local,NGLLX), &
             hetar_store(nrec_local,NGLLY), &
             hgammar_store(nrec_local,NGLLZ),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating receiver interpolators')

    allocate(hlagrange_store(NGLLX, NGLLY, NGLLZ, nrec_local), stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating array hlagrange_store')

    ! defines and stores Lagrange interpolators at all the receivers
    if (SIMULATION_TYPE == 2) then
      nadj_hprec_local = nrec_local
    else
      nadj_hprec_local = 1
    endif
    allocate(hpxir_store(nadj_hprec_local,NGLLX), &
             hpetar_store(nadj_hprec_local,NGLLY), &
             hpgammar_store(nadj_hprec_local,NGLLZ),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating derivative interpolators')

    ! stores interpolators for receiver positions
    call setup_sources_receivers_intp(NSOURCES, &
                      islice_selected_source, &
                      xi_source,eta_source,gamma_source, &
                      xigll,yigll,zigll, &
                      SIMULATION_TYPE,nrec,nrec_local, &
                      islice_selected_rec,number_receiver_global, &
                      xi_receiver,eta_receiver,gamma_receiver, &
                      hxir_store,hetar_store,hgammar_store, hlagrange_store, &
                      nadj_hprec_local,hpxir_store,hpetar_store,hpgammar_store)

    ! allocates seismogram array
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      allocate(seismograms(NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS),stat=ier)
      if (ier /= 0) stop 'Error while allocating seismograms'
    else
      ! adjoint seismograms
      allocate(seismograms(NDIM*NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS),stat=ier)
      if (ier /= 0) stop 'Error while allocating adjoint seismograms'

      ! allocates Frechet derivatives array
      allocate(moment_der(NDIM,NDIM,nrec_local), &
               sloc_der(NDIM,nrec_local), &
               stshift_der(nrec_local), &
               shdur_der(nrec_local),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating Frechet derivatives arrays')

      moment_der(:,:,:) = 0._CUSTOM_REAL
      sloc_der(:,:) = 0._CUSTOM_REAL
      stshift_der(:) = 0._CUSTOM_REAL
      shdur_der(:) = 0._CUSTOM_REAL
    endif
    ! initializes seismograms
    seismograms(:,:,:) = 0._CUSTOM_REAL
    ! adjoint seismograms
    it_adj_written = 0
  else
    ! allocates dummy array since we need it to pass as argument e.g. in write_seismograms() routine
    ! note: nrec_local is zero, Fortran 90/95 should allow zero-sized array allocation...
    allocate(seismograms(NDIM,0,NTSTEP_BETWEEN_OUTPUT_SEISMOS),stat=ier)
    if (ier /= 0) stop 'Error while allocating zero seismograms'
    ! dummy allocation
    allocate(hxir_store(1,1), &
             hetar_store(1,1), &
             hgammar_store(1,1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating dummy receiver interpolators')
  endif

  end subroutine setup_receivers_precompute_intp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources_receivers_intp(NSOURCES, &
                      islice_selected_source, &
                      xi_source,eta_source,gamma_source, &
                      xigll,yigll,zigll, &
                      SIMULATION_TYPE,nrec,nrec_local, &
                      islice_selected_rec,number_receiver_global, &
                      xi_receiver,eta_receiver,gamma_receiver, &
                      hxir_store,hetar_store,hgammar_store, hlagrange_store, &
                      nadj_hprec_local,hpxir_store,hpetar_store,hpgammar_store)

  use constants

  implicit none

  integer :: NSOURCES

  integer, dimension(NSOURCES) :: islice_selected_source

  double precision, dimension(NSOURCES) :: xi_source,eta_source,gamma_source
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll


  integer :: SIMULATION_TYPE

  integer :: nrec,nrec_local
  integer, dimension(nrec) :: islice_selected_rec
  integer, dimension(nrec_local) :: number_receiver_global
  double precision, dimension(nrec) :: xi_receiver,eta_receiver,gamma_receiver

  double precision, dimension(nrec_local,NGLLX) :: hxir_store
  double precision, dimension(nrec_local,NGLLY) :: hetar_store
  double precision, dimension(nrec_local,NGLLZ) :: hgammar_store
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nrec_local) :: hlagrange_store

  integer :: nadj_hprec_local
  double precision, dimension(nadj_hprec_local,NGLLX) :: hpxir_store
  double precision, dimension(nadj_hprec_local,NGLLY) :: hpetar_store
  double precision, dimension(nadj_hprec_local,NGLLZ) :: hpgammar_store


  ! local parameters
  integer :: isource,irec,irec_local, i, j, k
  double precision, dimension(NGLLX) :: hxir,hpxir
  double precision, dimension(NGLLY) :: hpetar,hetar
  double precision, dimension(NGLLZ) :: hgammar,hpgammar


  ! select local receivers

  ! define local to global receiver numbering mapping
  irec_local = 0
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    do irec = 1,nrec
      if (myrank == islice_selected_rec(irec)) then
        irec_local = irec_local + 1
        ! checks counter
        if (irec_local > nrec_local) call exit_MPI(myrank,'Error receiver interpolators: irec_local exceeds bounds')
        ! stores local to global receiver ids
        number_receiver_global(irec_local) = irec
      endif
    enddo
  else
    do isource = 1,NSOURCES
      if (myrank == islice_selected_source(isource)) then
        irec_local = irec_local + 1
        ! checks counter
        if (irec_local > nrec_local) call exit_MPI(myrank,'Error adjoint source interpolators: irec_local exceeds bounds')
        ! stores local to global receiver/source ids
        number_receiver_global(irec_local) = isource
      endif
    enddo
  endif
  ! checks if all local receivers have been found
  if (irec_local /= nrec_local) call exit_MPI(myrank,'Error number of local receivers do not match')

  ! define and store Lagrange interpolators at all the receivers
  do irec_local = 1,nrec_local
    irec = number_receiver_global(irec_local)

    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      ! receiver positions
      call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
      call lagrange_any(eta_receiver(irec),NGLLY,yigll,hetar,hpetar)
      call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
    else
      ! source positions
      call lagrange_any(xi_source(irec),NGLLX,xigll,hxir,hpxir)
      call lagrange_any(eta_source(irec),NGLLY,yigll,hetar,hpetar)
      call lagrange_any(gamma_source(irec),NGLLZ,zigll,hgammar,hpgammar)
    endif

    ! stores interpolators
    hxir_store(irec_local,:) = hxir(:)
    hetar_store(irec_local,:) = hetar(:)
    hgammar_store(irec_local,:) = hgammar(:)

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          hlagrange_store(i,j,k,irec_local) = hxir_store(irec_local,i)*hetar_store(irec_local,j)*hgammar_store(irec_local,k)
        enddo
      enddo
    enddo

    ! stores derivatives
    if (SIMULATION_TYPE == 2) then
      hpxir_store(irec_local,:) = hpxir(:)
      hpetar_store(irec_local,:) = hpetar(:)
      hpgammar_store(irec_local,:) = hpgammar(:)
    endif
  enddo

  end subroutine setup_sources_receivers_intp

