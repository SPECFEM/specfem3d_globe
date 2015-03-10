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
  double precision :: min_tshift_cmt_original
  integer :: isource
  character(len=MAX_STRING_LEN) :: filename
  integer :: ier

  ! makes smaller hdur for movies
  logical,parameter :: USE_SMALLER_HDUR_MOVIE = .true.

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'sources:'
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

  allocate(tshift_cmt(NSOURCES), &
           hdur(NSOURCES), &
           hdur_gaussian(NSOURCES),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating source arrays')

  allocate(theta_source(NSOURCES), &
           phi_source(NSOURCES),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating source arrays')

  allocate(nu_source(NDIM,NDIM,NSOURCES),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating source arrays')

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
                     ELLIPTICITY_VAL,min_tshift_cmt_original)

  if (abs(minval(tshift_cmt)) > TINYVAL) &
    call exit_MPI(myrank,'one tshift_cmt must be zero, others must be positive')

  ! count number of sources located in this slice
  nsources_local = 0
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    do isource = 1,NSOURCES
      if (myrank == islice_selected_source(isource)) nsources_local = nsources_local + 1
    enddo
  endif

  ! filter source time function by Gaussian with hdur = HDUR_MOVIE when writing movies or shakemaps
  if (MOVIE_SURFACE .or. MOVIE_VOLUME) then
    ! smaller hdur_movie will do
    if (USE_SMALLER_HDUR_MOVIE) then
      ! hdur_movie gets assigned an automatic value based on the simulation resolution
      ! this will make that a bit smaller to have a higher-frequency movie output
      HDUR_MOVIE = 0.5* HDUR_MOVIE
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
  hdur_gaussian(:) = hdur(:)/SOURCE_DECAY_MIMIC_TRIANGLE

  ! define t0 as the earliest start time
  t0 = - 1.5d0*minval( tshift_cmt(:) - hdur(:) )

  ! uses an external file for source time function, which starts at time 0.0
  if (EXTERNAL_SOURCE_TIME_FUNCTION) then
    hdur(:) = 0._CUSTOM_REAL
    t0      = 0.d0
  endif

  ! point force sources will start depending on the frequency given by hdur
  if (USE_FORCE_POINT_SOURCE) then
    ! note: point force sources will give the dominant frequency in hdur,
    !          thus the main period is 1/hdur.
    !          also, these sources use a Ricker source time function instead of a Gaussian.
    !          for a Ricker source time function, a start time ~1.2 * main_period is a good choice
    t0 = - 1.2d0 * minval(tshift_cmt(:) - 1.0d0/hdur(:))
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
      write(IMAIN,*) 't0: ',t0,'min_tshift_cmt_original: ',min_tshift_cmt_original
      write(IMAIN,*)
    endif

    ! checks if automatically set t0 is too small
    ! note: min_tshift_cmt_original can be a positive or negative time shift (minimum from all tshift)
    if (t0 <= USER_T0 + min_tshift_cmt_original) then
      ! by default, tshift_cmt(:) holds relative time shifts with a minimum time shift set to zero
      ! re-adds (minimum) original time shift such that sources will kick in
      ! according to their absolute time shift
      tshift_cmt(:) = tshift_cmt(:) + min_tshift_cmt_original

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
        write(IMAIN,*) '       - increase USER_T0 to be at least: ',t0-min_tshift_cmt_original
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

  ! determines number of times steps for simulation
  call setup_timesteps()

  ! prints source time functions to output files
  if (PRINT_SOURCE_TIME_FUNCTION .and. myrank == 0) then
    do isource = 1,NSOURCES
        ! print source time function and spectrum
         call print_stf(NSOURCES,isource,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                        tshift_cmt,hdur,min_tshift_cmt_original,NSTEP,DT)
    enddo
  endif

  ! get information about event name and location
  ! (e.g. needed for SAC seismograms)

  ! The following line is added for get_event_info subroutine.
  ! Because the way NSOURCES_SAC was declared has been changed.
  ! The rest of the changes in this program is just the updates of the subroutines that
  ! I did changes, e.g., adding/removing parameters. by Ebru Bozdag
  call get_event_info_parallel(myrank,yr_SAC,jda_SAC,ho_SAC,mi_SAC,sec_SAC,&
                              event_name_SAC,t_cmt_SAC,t_shift_SAC, &
                              elat_SAC,elon_SAC,depth_SAC,mb_SAC,cmt_lat_SAC,&
                              cmt_lon_SAC,cmt_depth_SAC,cmt_hdur_SAC,NSOURCES)

  end subroutine setup_sources

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
  NSTEP = NSTEP + 100 * (int( abs(t0) / (100.d0*DT)) + 1)

  ! if doing benchmark runs to measure scaling of the code for a limited number of time steps only
  if (DO_BENCHMARK_RUN_ONLY) NSTEP = NSTEP_FOR_BENCHMARK

  ! noise tomography:
  ! time steps needs to be doubled, due to +/- branches
  if (NOISE_TOMOGRAPHY /= 0 )   NSTEP = 2*NSTEP-1

!! DK DK make sure NSTEP is a multiple of NT_DUMP_ATTENUATION
  if (UNDO_ATTENUATION .and. mod(NSTEP,NT_DUMP_ATTENUATION) /= 0) then
    NSTEP = (NSTEP/NT_DUMP_ATTENUATION + 1)*NT_DUMP_ATTENUATION
  endif
  it_end = NSTEP

  ! subsets used to save seismograms must not be larger than the whole time series,
  ! otherwise we waste memory
  if (NTSTEP_BETWEEN_OUTPUT_SEISMOS > NSTEP .or. is_initial_guess) NTSTEP_BETWEEN_OUTPUT_SEISMOS = NSTEP

  ! re-checks output steps?
  !if (OUTPUT_SEISMOS_SAC_ALPHANUM .and. (mod(NTSTEP_BETWEEN_OUTPUT_SEISMOS,5) /= 0)) &
  !  stop 'if OUTPUT_SEISMOS_SAC_ALPHANUM = .true. then NTSTEP_BETWEEN_OUTPUT_SEISMOS must be a multiple of 5, check the Par_file'

  ! subsets used to save adjoint sources must not be larger than the whole time series,
  ! otherwise we waste memory
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    if (NTSTEP_BETWEEN_READ_ADJSRC > NSTEP) NTSTEP_BETWEEN_READ_ADJSRC = NSTEP
  endif

  end subroutine setup_timesteps

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_receivers()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  integer :: irec,isource,nrec_tot_found
  integer :: nrec_simulation
  integer :: nadj_files_found,nadj_files_found_tot
  integer :: ier
  integer,dimension(:),allocatable :: tmp_rec_local_all
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
        call check_adjoint_sources(irec,nadj_files_found)
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
  ! gathers info about receivers on master
  if (myrank == 0) then
    ! only master process needs full arrays allocated
    allocate(tmp_rec_local_all(NPROCTOT_VAL),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_rec_local_all')
  else
    ! dummy arrays
    allocate(tmp_rec_local_all(1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary array tmp_rec_local_all')
  endif

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
  tmp_rec_local_all(1) = nrec_local
  if (NPROCTOT_VAL > 1) then
    call gather_all_singlei(nrec_local,tmp_rec_local_all,NPROCTOT_VAL)
  endif
  ! user output
  if (myrank == 0) then
    ! determines maximum number of local receivers and corresponding rank
    maxrec = maxval(tmp_rec_local_all(:))
    maxproc = maxloc(tmp_rec_local_all(:))
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
    tmp_rec_local_all(1) = nadj_rec_local
    if (NPROCTOT_VAL > 1) then
      call gather_all_singlei(nadj_rec_local,tmp_rec_local_all,NPROCTOT_VAL)
    endif
    ! user output
    if (myrank == 0) then
      ! determines maximum number of local receivers and corresponding rank
      maxrec = maxval(tmp_rec_local_all(:))
      maxproc = maxloc(tmp_rec_local_all(:))
      !do i = 1, NPROCTOT_VAL
      !  if (tmp_rec_local_all(i) > maxrec) then
      !    maxrec = tmp_rec_local_all(i)
      !    maxproc = i-1
      !  endif
      !enddo
      ! adj_sourcearrays size in MB
      ! adj_sourcearrays(NDIM,NGLLX,NGLLY,NGLLZ,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC)
      sizeval = dble(maxrec) * dble(NDIM * NGLLX * NGLLY * NGLLZ * NTSTEP_BETWEEN_READ_ADJSRC * CUSTOM_REAL / 1024. / 1024. )
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

  deallocate(tmp_rec_local_all)

  end subroutine setup_receivers

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources_receivers_VTKfile()

  use specfem_par,only: myrank,OUTPUT_FILES,NSOURCES,nrec,MAX_STRING_LEN
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
  "('sed -e ',a1,'s/POINTS.*/POINTS',i6,' float/',a1,' < ',a,' > ',a)")&
      "'",NSOURCES + nrec,"'",trim(filename),trim(filename_new)

    ! note: this system() routine is non-standard fortran
    call system_command(command)

    ! only extract receiver locations and remove temporary file
    filename_new = trim(OUTPUT_FILES)//'/receiver.vtk'
    write(command, &
  "('awk ',a1,'{if (NR<5) print $0;if (NR==6)print ',a1,'POINTS',i6,' float',a1,';if (NR>5+',i6,')print $0}',a1,' < ',a,' > ',a)")&
      "'",'"',nrec,'"',NSOURCES,"'",trim(filename),trim(filename_new)

    ! note: this system() routine is non-standard fortran
    call system_command(command)

    ! only extract source locations and remove temporary file
    filename_new = trim(OUTPUT_FILES)//'/source.vtk'
    write(command, &
  "('awk ',a1,'{if (NR< 6 + ',i6,') print $0}END{print}',a1,' < ',a,' > ',a,'; rm -f ',a)")&
      "'",NSOURCES,"'",trim(filename),trim(filename_new),trim(filename)

    ! note: this system() routine is non-standard fortran
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
  if (SIMULATION_TYPE == 1  .or. SIMULATION_TYPE == 3) then
    ! source interpolated on all GLL points in source element
    allocate(sourcearrays(NDIM,NGLLX,NGLLY,NGLLZ,NSOURCES),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating sourcearrays')
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
      ! allocate adjoint source arrays
      allocate(adj_sourcearrays(NDIM,NGLLX,NGLLY,NGLLZ,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC), &
               stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating adjoint sourcearrays')
      adj_sourcearrays(:,:,:,:,:,:) = 0._CUSTOM_REAL

      ! additional buffer for asynchronous file i/o
      if (IO_ASYNC_COPY .and. NSTEP_SUB_ADJ > 1) then
        ! allocates read buffer
        allocate(buffer_sourcearrays(NDIM,NGLLX,NGLLY,NGLLZ,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC), &
                 stat=ier)
        if (ier /= 0 ) call exit_MPI(myrank,'Error allocating array buffer_sourcearrays')

        ! array size in bytes (note: the multiplication is split into two line to avoid integer-overflow)
        arraysize = NDIM * NGLLX * NGLLY * NGLLZ * CUSTOM_REAL
        arraysize = arraysize * nadj_rec_local * NTSTEP_BETWEEN_READ_ADJSRC

        ! debug
        !print*,'buffer_sourcearrays: size = ',arraysize,' Bytes = ',arraysize/1024./1024.,'MB'

        ! initializes io thread
        call prepare_adj_io_thread(buffer_sourcearrays,arraysize,nadj_rec_local)
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
  integer :: isource,iglob,i,j,k,ispec
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearray
  double precision :: xi,eta,gamma

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

      ! pre-computes source contribution on GLL points
      call compute_arrays_source(sourcearray,xi,eta,gamma, &
                          Mxx(isource),Myy(isource),Mzz(isource),Mxy(isource),Mxz(isource),Myz(isource), &
                          xix_crust_mantle(:,:,:,ispec),xiy_crust_mantle(:,:,:,ispec),xiz_crust_mantle(:,:,:,ispec), &
                          etax_crust_mantle(:,:,:,ispec),etay_crust_mantle(:,:,:,ispec),etaz_crust_mantle(:,:,:,ispec), &
                          gammax_crust_mantle(:,:,:,ispec),gammay_crust_mantle(:,:,:,ispec),gammaz_crust_mantle(:,:,:,ispec), &
                          xigll,yigll,zigll)

      ! point forces, initializes sourcearray, used for simplified CUDA routines
      if (USE_FORCE_POINT_SOURCE) then
        ! note: for use_force_point_source xi/eta/gamma are in the range [1,NGLL*]
        iglob = ibool_crust_mantle(nint(xi),nint(eta),nint(gamma),ispec)

        ! sets sourcearrays
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              if (ibool_crust_mantle(i,j,k,ispec) == iglob) then
                ! elastic source components
                sourcearray(:,i,j,k) = nu_source(COMPONENT_FORCE_SOURCE,:,isource)
              endif
            enddo
          enddo
        enddo
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
    call setup_sources_receivers_intp(NSOURCES,myrank, &
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
      allocate(moment_der(NDIM,NDIM,nrec_local),sloc_der(NDIM,nrec_local), &
               stshift_der(nrec_local),shdur_der(nrec_local),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating Frechet derivatives arrays')

      moment_der(:,:,:) = 0._CUSTOM_REAL
      sloc_der(:,:) = 0._CUSTOM_REAL
      stshift_der(:) = 0._CUSTOM_REAL
      shdur_der(:) = 0._CUSTOM_REAL
    endif
    ! initializes seismograms
    seismograms(:,:,:) = 0._CUSTOM_REAL
    nit_written = 0
  else
    ! allocates dummy array since we need it to pass as argument e.g. in write_seismograms() routine
    ! note: nrec_local is zero, fortran 90/95 should allow zero-sized array allocation...
    allocate(seismograms(NDIM,0,NTSTEP_BETWEEN_OUTPUT_SEISMOS),stat=ier)
    if (ier /= 0) stop 'Error while allocating zero seismograms'
  endif

  end subroutine setup_receivers_precompute_intp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources_receivers_intp(NSOURCES,myrank, &
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

  integer :: NSOURCES,myrank

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

