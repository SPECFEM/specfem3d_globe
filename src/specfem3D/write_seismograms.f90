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

  subroutine write_seismograms()

  use constants, only: IMAIN,OUTPUT_ADJOINT_WAVEFIELD_SEISMOGRAMS

  use constants_solver, only: NGLOB_CRUST_MANTLE,NGLOB_CRUST_MANTLE_ADJOINT

  use specfem_par, only: myrank,Mesh_pointer,GPU_MODE,GPU_ASYNC_COPY,SIMULATION_TYPE, &
    nrec_local,number_receiver_global,ispec_selected_rec,ispec_selected_source, &
    it,it_end, &
    seismo_current,seismo_offset, &
    seismograms, &
    nlength_seismogram, &
    NTSTEP_BETWEEN_OUTPUT_SEISMOS,NTSTEP_BETWEEN_OUTPUT_SAMPLE, &
    do_save_seismograms, &
    WRITE_SEISMOGRAMS_BY_MAIN,OUTPUT_SEISMOS_ASDF, &
    SAVE_SEISMOGRAMS_IN_ADJOINT_RUN,SAVE_SEISMOGRAMS_STRAIN, &
    moment_der,sloc_der,shdur_der,stshift_der, &
    scale_displ,NSTEP

  use specfem_par_crustmantle, only: displ_crust_mantle,b_displ_crust_mantle, &
    eps_trace_over_3_crust_mantle,epsilondev_xx_crust_mantle,epsilondev_xy_crust_mantle,epsilondev_xz_crust_mantle, &
    epsilondev_yy_crust_mantle,epsilondev_yz_crust_mantle, &
    ibool_crust_mantle

  implicit none

  ! local parameters
  ! timing
  double precision, external :: wtime
  double precision :: write_time_begin,write_time

  ! checks if anything to do
  if (.not. do_save_seismograms) return

  ! checks subsampling recurrence
  if (mod(it-1,NTSTEP_BETWEEN_OUTPUT_SAMPLE) == 0) then

    ! update position in seismograms
    seismo_current = seismo_current + 1

    ! check for edge effects
    if (seismo_current < 1 .or. seismo_current > nlength_seismogram) &
      call exit_mpi(myrank,'Error: seismo_current out of bounds in recording of seismograms')

    ! compute & store the seismograms only if there is at least one receiver located in this slice
    if (nrec_local > 0) then
      ! gets resulting array values onto CPU
      if (GPU_MODE) then
        ! for forward and kernel simulations, seismograms are computed by the GPU, thus no need to transfer the wavefield
        ! gets field values from GPU
        if (SIMULATION_TYPE == 2 .or. SAVE_SEISMOGRAMS_STRAIN) then
          ! this transfers fields only in elements with stations for efficiency
          call write_seismograms_transfer_gpu(Mesh_pointer, &
                                              displ_crust_mantle,b_displ_crust_mantle, &
                                              eps_trace_over_3_crust_mantle, &
                                              epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                                              epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                                              number_receiver_global, &
                                              ispec_selected_rec,ispec_selected_source, &
                                              ibool_crust_mantle)
          ! synchronizes field values from GPU
          if (GPU_ASYNC_COPY) then
            call transfer_seismo_from_device_async(Mesh_pointer, &
                                                   displ_crust_mantle,b_displ_crust_mantle, &
                                                   number_receiver_global,ispec_selected_rec,ispec_selected_source, &
                                                   ibool_crust_mantle)
          endif
        endif
      endif ! GPU_MODE

      ! computes traces at interpolated receiver locations
      select case (SIMULATION_TYPE)
      case (1)
        ! forward run
        if (.not. GPU_MODE) then
          ! on CPU
          call compute_seismograms(NGLOB_CRUST_MANTLE,displ_crust_mantle,seismo_current,seismograms)
        else
          ! on GPU
          call compute_seismograms_gpu(Mesh_pointer,seismograms,seismo_current,it,it_end, &
                                       scale_displ,nlength_seismogram,NSTEP)

        endif
      case (2)
        ! adjoint run
        call compute_seismograms_adjoint(displ_crust_mantle, &
                                         eps_trace_over_3_crust_mantle, &
                                         epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                                         epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                                         moment_der,sloc_der,stshift_der,shdur_der, &
                                         seismograms)
      case (3)
        ! kernel run
        if (.not. GPU_MODE) then
          ! on CPU
          if (.not. ( SIMULATION_TYPE == 3 .and. (.not. SAVE_SEISMOGRAMS_IN_ADJOINT_RUN)) ) then
            if (OUTPUT_ADJOINT_WAVEFIELD_SEISMOGRAMS) then
              ! uncomment to output adjoint wavefield instead for seismogram output
              call compute_seismograms(NGLOB_CRUST_MANTLE_ADJOINT,displ_crust_mantle,seismo_current,seismograms)
            else
              ! default, backward reconstructed wavefield seismos
              call compute_seismograms(NGLOB_CRUST_MANTLE_ADJOINT,b_displ_crust_mantle,seismo_current,seismograms)
            endif
          endif
        else
          ! on GPU
          if (.not. ( SIMULATION_TYPE == 3 .and. (.not. SAVE_SEISMOGRAMS_IN_ADJOINT_RUN)) ) then
            call compute_seismograms_gpu(Mesh_pointer,seismograms,seismo_current,it,it_end, &
                                         scale_displ,nlength_seismogram,NSTEP)
          endif
        endif
      end select

      ! strain seismograms
      if (SAVE_SEISMOGRAMS_STRAIN) then
        select case (SIMULATION_TYPE)
        case (1)
          ! forward run
          call compute_seismograms_strain(NGLOB_CRUST_MANTLE,displ_crust_mantle)
        case (3)
          ! kernel run
          call compute_seismograms_strain(NGLOB_CRUST_MANTLE,b_displ_crust_mantle)
        end select
      endif

    endif ! nrec_local
  endif

  ! write the current or final seismograms
  !
  ! for example:
  !   NTSTEP_BETWEEN_OUTPUT_SEISMOS = 10 , NTSTEP_BETWEEN_OUTPUT_SAMPLE = 5
  !       -> 10 / 5 == 2 steps (nlength_seismogram)
  !          we will store samples at it==1,it==6 and output at it==10
  !   NTSTEP_BETWEEN_OUTPUT_SEISMOS = 10 , NTSTEP_BETWEEN_OUTPUT_SAMPLE = 1
  !       -> 10 / 1 == 10 steps
  !          we will store samples at it==1,2,3,..,10 and output at it==10
  !
  !   that is for NTSTEP_BETWEEN_OUTPUT_SEISMOS = 10 , NTSTEP_BETWEEN_OUTPUT_SAMPLE = 5,
  !   the seismograms have samples for the wavefield at mod((it-1),NTSTEP_BETWEEN_OUTPUT_SAMPLE) == 0,
  !   which is at it==1,it==6. for writing the seismograms however,
  !   we would write out at mod(it,NTSTEP_BETWEEN_OUTPUT_SEISMOS) == 0, which is at it==10
  !
  if (mod(it,NTSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == it_end) then
    ! timing
    write_time_begin = wtime()

    ! checks if anything to do
    ! note: ASDF uses parallel hdf5 that defines the MPI communicator group that the solver is
    !       run with. this means every processor in the group is needed for write_seismograms
    if (nrec_local > 0 .or. ( WRITE_SEISMOGRAMS_BY_MAIN .and. myrank == 0 ) .or. OUTPUT_SEISMOS_ASDF) then
      ! writes out seismogram files
      select case (SIMULATION_TYPE)
      case (1,3)
        ! forward/reconstructed wavefields
        if (.not. ( SIMULATION_TYPE == 3 .and. (.not. SAVE_SEISMOGRAMS_IN_ADJOINT_RUN) ) ) &
          call write_seismograms_to_file()
        if (SAVE_SEISMOGRAMS_STRAIN) &
          call write_seismograms_strain()
      case (2)
        ! adjoint wavefield
        call write_adj_seismograms()
      end select
    endif

    ! synchronizes processes (waits for all processes to finish writing)
    call synchronize_all()

    ! resets current seismogram position
    seismo_offset = seismo_offset + seismo_current
    seismo_current = 0

    ! user output
    if (myrank == 0) then
      ! timing
      write_time = wtime() - write_time_begin
      ! output
      write(IMAIN,*) 'Total number of time steps written: ', seismo_offset
      if (WRITE_SEISMOGRAMS_BY_MAIN) then
        write(IMAIN,*) 'Writing the seismograms by main proc alone took ',sngl(write_time),' seconds'
      else
        write(IMAIN,*) 'Writing the seismograms in parallel took ',sngl(write_time),' seconds'
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  end subroutine write_seismograms

!
!-------------------------------------------------------------------------------------------------
!

! write seismograms to files
  subroutine write_seismograms_to_file()

  use constants_solver, only: MAX_STRING_LEN,CUSTOM_REAL,NDIM,IMAIN,IOUT,itag

  use specfem_par, only: &
    NPROCTOT_VAL,myrank,nrec,nrec_local, &
    number_receiver_global, &
    seismograms, &
    nlength_seismogram, &
    seismo_offset,seismo_current, &
    islice_num_rec_local, &
    station_name,network_name, &
    OUTPUT_SEISMOS_ASCII_TEXT, &
    OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
    OUTPUT_SEISMOS_ASDF, &
    OUTPUT_SEISMOS_3D_ARRAY, &
    SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE, &
    OUTPUT_FILES, &
    WRITE_SEISMOGRAMS_BY_MAIN, &
    DT,NSTEP,NTSTEP_BETWEEN_OUTPUT_SAMPLE

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: one_seismogram

  integer :: iproc,sender,irec_local,irec,ier,receiver
  integer :: nrec_local_received
  integer :: total_seismos,length_network_name,length_station_name
  character(len=MAX_STRING_LEN) :: sisname,staname

  ! allocates single station seismogram
  allocate(one_seismogram(NDIM,nlength_seismogram),stat=ier)
  if (ier /= 0) call exit_mpi(myrank,'Error while allocating one temporary seismogram')
  one_seismogram(:,:) = 0.0_CUSTOM_REAL

  ! ASDF format
  if (OUTPUT_SEISMOS_ASDF) then
    ! The writing of seismograms by the main proc is handled within write_asdf()
    do irec_local = 1,nrec_local

      ! get global number of that receiver
      irec = number_receiver_global(irec_local)

      one_seismogram(:,:) = seismograms(:,irec_local,:)

      ! write this seismogram
      ! note: ASDF data structure is given in module
      !       stores all traces into ASDF container in case
      call write_one_seismogram(one_seismogram,irec,irec_local,.true.)
    enddo

    ! writes ASDF file output
    call write_asdf()

    call synchronize_all()

    ! deallocate the container
    call close_asdf_data()
  endif

  ! write 3D seismogram array
  if (OUTPUT_SEISMOS_3D_ARRAY) then
    write(sisname,'(A,I5.5)') '/array_seismograms_node_',myrank
    write(staname,'(A,I5.5)') '/array_stations_node_',myrank
    if (seismo_offset == 0) then
      open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.bin',status='unknown',form='unformatted',action='write')
    else
      open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.bin',status='old', &
            form='unformatted',position='append',action='write')
    endif
    write(IOUT) seismograms(:,:,1:seismo_current:NTSTEP_BETWEEN_OUTPUT_SAMPLE)
    close(IOUT)
    ! save list of stations in current processor
    if (seismo_offset == 0) then
      if (myrank == 0) then
        open(unit=IOUT,file=trim(OUTPUT_FILES)//'seismogram_stats.txt',status='unknown',form='formatted',action='write')
        write(IOUT,*) 'DT0    =', DT
        write(IOUT,*) 'NSTEP0 =', NSTEP
        write(IOUT,*) 'DT     =', DT * NTSTEP_BETWEEN_OUTPUT_SAMPLE
        write(IOUT,*) 'NSTEP  =', ceiling(real(seismo_current) / NTSTEP_BETWEEN_OUTPUT_SAMPLE)
        close(IOUT)
      endif
      open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(staname)//'.txt',status='unknown',form='formatted',action='write')
    else
      if (myrank == 0) then
        open(unit=IOUT,file=trim(OUTPUT_FILES)//'seismogram_stats.txt',status='old', &
            form='formatted',position='append',action='write')
        write(IOUT,*) 'NSTEP =', ceiling(real(seismo_current) / NTSTEP_BETWEEN_OUTPUT_SAMPLE)
        close(IOUT)
      endif
      open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(staname)//'.txt',status='old', &
            form='formatted',position='append',action='write')
    endif
    do irec_local = 1,nrec_local
      ! get global number of that receiver
      irec = number_receiver_global(irec_local)
      length_station_name = len_trim(station_name(irec))
      length_network_name = len_trim(network_name(irec))

      write(sisname,"('# ',a,'.',a)") network_name(irec)(1:length_network_name), &
                   station_name(irec)(1:length_station_name)
      write(IOUT,*) sisname(1:len_trim(sisname))
    enddo
    close(IOUT)

  else if (NTSTEP_BETWEEN_OUTPUT_SAMPLE > 1) then
    ! save original DT and NSTEP for adjoint simulation
    if (myrank == 0 .and. seismo_offset == 0) then
      open(unit=IOUT,file=trim(OUTPUT_FILES)//'seismogram_stats.txt',status='unknown',form='formatted',action='write')
      write(IOUT,*) 'DT0     =', DT
      write(IOUT,*) 'NSTEP0  =', NSTEP
      close(IOUT)
    endif
  endif

  ! ASCII / SAC format
  if (OUTPUT_SEISMOS_ASCII_TEXT .or. OUTPUT_SEISMOS_SAC_ALPHANUM .or. OUTPUT_SEISMOS_SAC_BINARY) then

    ! write out seismograms: all processes write their local seismograms themselves
    if (.not. WRITE_SEISMOGRAMS_BY_MAIN) then

      ! all the processes write their local seismograms themselves
      if (SAVE_ALL_SEISMOS_IN_ONE_FILE .and. OUTPUT_SEISMOS_ASCII_TEXT) then
        write(sisname,'(A,I5.5)') '/all_seismograms_node_',myrank

        if (USE_BINARY_FOR_LARGE_FILE) then
          if (seismo_offset == 0) then
            open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.bin',status='unknown',form='unformatted',action='write')
          else
            open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.bin',status='old', &
                 form='unformatted',position='append',action='write')
          endif
        else
          if (seismo_offset == 0) then
            open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.ascii',status='unknown',form='formatted',action='write')
          else
            open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.ascii',status='old', &
                 form='formatted',position='append',action='write')
          endif
        endif
      endif

      ! loop on all the local receivers
      do irec_local = 1,nrec_local

        ! get global number of that receiver
        irec = number_receiver_global(irec_local)

        one_seismogram(:,:) = seismograms(:,irec_local,:)

        ! write this seismogram
        ! note: ASDF data structure is given in module
        !       stores all traces into ASDF container in case
        call write_one_seismogram(one_seismogram,irec,irec_local,.false.)
      enddo

      ! create one large file instead of one small file per station to avoid file system overload
      if (SAVE_ALL_SEISMOS_IN_ONE_FILE .and. OUTPUT_SEISMOS_ASCII_TEXT) close(IOUT)

    else
      ! WRITE_SEISMOGRAMS_BY_MAIN
      ! only the main process does the writing of seismograms

      ! opens file for single file output
      if (SAVE_ALL_SEISMOS_IN_ONE_FILE .and. OUTPUT_SEISMOS_ASCII_TEXT) then
        if (myrank == 0) then
          ! create one large file instead of one small file per station to avoid file system overload
          write(sisname,'(A)') '/all_seismograms'

          if (USE_BINARY_FOR_LARGE_FILE) then
            if (seismo_offset == 0) then
              open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.bin',status='unknown', &
                   form='unformatted',action='write')
            else
              open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.bin',status='old', &
                   form='unformatted',position='append',action='write')
            endif
          else
            if (seismo_offset == 0) then
              open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.ascii',status='unknown', &
                   form='formatted',action='write')
            else
              open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.ascii',status='old', &
                   form='formatted',position='append',action='write')
            endif
          endif
        endif ! myrank
      endif

      ! collects the data from all other processes
      if (myrank == 0) then
        ! on the main, gather all the seismograms
        total_seismos = 0

        ! loop on all the slices
        do iproc = 0,NPROCTOT_VAL-1

          ! communicates only with processes that contain local receivers (to minimize MPI chatter)
          if (islice_num_rec_local(iproc) == 0 ) cycle

          ! receive except from proc 0, which is me and therefore I already have this value
          sender = iproc
          if (iproc == 0) then
            ! main is current slice
            nrec_local_received = nrec_local
          else
            ! receives info from secondary processes
            call recv_singlei(nrec_local_received,sender,itag)
            if (nrec_local_received <= 0) call exit_MPI(myrank,'Error while receiving local number of receivers')
          endif

          if (nrec_local_received > 0) then
            do irec_local = 1,nrec_local_received
              ! init trace
              one_seismogram(:,:) = 0._CUSTOM_REAL

              ! receive except from proc 0, which is myself and therefore I already have these values
              if (iproc == 0) then
                ! get global number of that receiver
                irec = number_receiver_global(irec_local)
                one_seismogram(:,:) = seismograms(:,irec_local,:)
              else
                ! receives info from secondary processes
                call recv_singlei(irec,sender,itag)
                if (irec < 1 .or. irec > nrec) call exit_MPI(myrank,'Error while receiving global receiver number')
                call recv_cr(one_seismogram,NDIM*seismo_current,sender,itag)
              endif

              ! write this seismogram
              call write_one_seismogram(one_seismogram,irec,irec_local,.false.)

              ! counts seismos written
              total_seismos = total_seismos + 1
            enddo
          endif
        enddo

      else
        ! on the nodes, send the seismograms to the main
        receiver = 0
        ! only sends if this slice contains receiver stations
        if (nrec_local > 0) then
          call send_singlei(nrec_local,receiver,itag)
          do irec_local = 1,nrec_local
            ! get global number of that receiver
            irec = number_receiver_global(irec_local)
            call send_singlei(irec,receiver,itag)

            one_seismogram(:,:) = seismograms(:,irec_local,:)
            call send_cr(one_seismogram,NDIM*seismo_current,receiver,itag)
          enddo
        endif
      endif

      ! only main process
      if (myrank == 0) then
        ! output info
        write(IMAIN,*) 'Total number of receivers saved is ',total_seismos,' out of ',nrec
        call flush_IMAIN()

        ! checks
        if (total_seismos /= nrec) call exit_MPI(myrank,'incorrect total number of receivers saved')

        ! create one large file instead of one small file per station to avoid file system overload
        if (SAVE_ALL_SEISMOS_IN_ONE_FILE) close(IOUT)
      endif

    endif ! WRITE_SEISMOGRAMS_BY_MAIN

  endif ! ASCII / SAC format

  deallocate(one_seismogram)

  end subroutine write_seismograms_to_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_one_seismogram(one_seismogram,irec,irec_local,is_for_asdf)

  use constants_solver, only: MAX_STRING_LEN,CUSTOM_REAL,NDIM,DEGREES_TO_RADIANS, &
    MAX_LENGTH_STATION_NAME,MAX_LENGTH_NETWORK_NAME

  use specfem_par, only: &
    myrank, &
    station_name,network_name,stlat,stlon, &
    DT, &
    seismo_current, &
    OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_ASDF, &
    OUTPUT_SEISMOS_SAC_BINARY,ROTATE_SEISMOGRAMS_RT, &
    nlength_seismogram

  use specfem_par, only: &
    cmt_lat => cmt_lat_SAC,cmt_lon => cmt_lon_SAC

  implicit none

  integer,intent(in) :: irec,irec_local
  logical,intent(in) :: is_for_asdf
  real(kind=CUSTOM_REAL), dimension(NDIM,nlength_seismogram),intent(in) :: one_seismogram

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(5,nlength_seismogram) :: seismogram_tmp
  integer :: iorientation,length_station_name,length_network_name

  character(len=4) :: chn
  character(len=MAX_STRING_LEN) :: sisname,sisname_big_file
  character(len=2) :: bic

  ! variables used for calculation of backazimuth and
  ! rotation of components if ROTATE_SEISMOGRAMS=.true.
  integer :: ior_start,ior_end
  double precision :: backaz
  double precision :: phi
  real(kind=CUSTOM_REAL) :: cphi,sphi
  integer :: isample

  ! initializes
  seismogram_tmp(:,:) = 0.0_CUSTOM_REAL

  ! get band code
  call band_instrument_code(DT,bic)

  if (ROTATE_SEISMOGRAMS_RT) then ! iorientation 1=N,2=E,3=Z,4=R,5=T
    ior_start = 3    ! starting from Z
    ior_end   = 5    ! ending with T => ZRT
  else
    ior_start = 1    ! starting from N
    ior_end   = 3    ! ending with Z => NEZ
  endif

  do iorientation = ior_start,ior_end    ! changes according to ROTATE_SEISMOGRAMS_RT

    select case (iorientation)
    case (1)
      !chn = 'LHN'
      chn = bic(1:2)//'N'
    case (2)
      !chn = 'LHE'
      chn = bic(1:2)//'E'
    case (3)
      !chn = 'LHZ'
      chn = bic(1:2)//'Z'
    case (4)
      !chn = 'LHR'
      chn = bic(1:2)//'R'
    case (5)
      !chn = 'LHT'
      chn = bic(1:2)//'T'
    case default
      call exit_MPI(myrank,'incorrect channel value')
    end select

    if (iorientation == 4 .or. iorientation == 5) then
      ! calculate backazimuth needed to rotate East and North
      ! components to Radial and Transverse components
      ! (back-azimuth returned in degrees between [0,360])
      call get_backazimuth(cmt_lat,cmt_lon,stlat(irec),stlon(irec),backaz)

      phi = backaz

      ! back azimuth is the incoming direction of a raypath to a receiving station, i.e. the angle of
      ! the incoming wave front arriving at the station measured between north and the direction to the epicenter in degrees.
      ! (north corresponds to zero degrees)

      ! rotation angle phi takes opposite direction; to have radial direction pointing in outgoing direction
      if (phi > 180.d0) then
         phi = phi-180.d0
      else if (phi < 180.d0) then
         phi = phi+180.d0
      else if (phi == 180.d0) then
         phi = 0.d0
      endif

      cphi = cos(phi*DEGREES_TO_RADIANS)
      sphi = sin(phi*DEGREES_TO_RADIANS)

      ! do the rotation of the components and put result in
      ! new variable seismogram_tmp
      if (iorientation == 4) then ! radial component
        do isample = 1,seismo_current
          seismogram_tmp(iorientation,isample) = &
               cphi * one_seismogram(1,isample) + sphi * one_seismogram(2,isample)
        enddo
      else if (iorientation == 5) then ! transverse component
        do isample = 1,seismo_current
          seismogram_tmp(iorientation,isample) = &
            -1*sphi * one_seismogram(1,isample) + cphi * one_seismogram(2,isample)
        enddo
      endif
    else
      ! keep NEZ components
      do isample = 1,seismo_current
        seismogram_tmp(iorientation,isample) = one_seismogram(iorientation,isample)
      enddo
    endif

    ! create the name of the seismogram file for each slice
    ! file name includes the name of the station and the network
    length_station_name = len_trim(station_name(irec))
    length_network_name = len_trim(network_name(irec))

    ! check that length conforms to standard
    if (length_station_name < 1 .or. length_station_name > MAX_LENGTH_STATION_NAME) &
           call exit_MPI(myrank,'wrong length of station name')

    if (length_network_name < 1 .or. length_network_name > MAX_LENGTH_NETWORK_NAME) &
           call exit_MPI(myrank,'wrong length of network name')

    ! create the name of the seismogram file using the station name and network name
    ! using format: **net**.**sta**.channel
    write(sisname,"('/',a,'.',a,'.',a3,'.sem')") network_name(irec)(1:length_network_name), &
                   station_name(irec)(1:length_station_name),chn

    ! create this name also for the text line added to the unique big seismogram file
    write(sisname_big_file,"(a,'.',a,'.',a3,'.sem')") network_name(irec)(1:length_network_name), &
                   station_name(irec)(1:length_station_name),chn

    if (is_for_asdf) then
      ! ASDF output format
      if (OUTPUT_SEISMOS_ASDF) then
        call store_asdf_data(seismogram_tmp,irec_local,irec,chn,iorientation)
      endif
    else
      ! SAC output format
      if (OUTPUT_SEISMOS_SAC_ALPHANUM .or. OUTPUT_SEISMOS_SAC_BINARY ) then
        call write_output_SAC(seismogram_tmp,irec,iorientation,sisname,chn,phi)
      endif

      ! ASCII output format
      if (OUTPUT_SEISMOS_ASCII_TEXT) then
        call write_output_ASCII(seismogram_tmp,iorientation,sisname,sisname_big_file)
      endif
    endif

  enddo ! do iorientation

  end subroutine write_one_seismogram

!
!-------------------------------------------------------------------------------------------------
!

! write adjoint seismograms to text files

  subroutine write_adj_seismograms()

  use constants, only: MAX_STRING_LEN,CUSTOM_REAL,IOUT

  use specfem_par, only: NTSTEP_BETWEEN_OUTPUT_SEISMOS,NTSTEP_BETWEEN_OUTPUT_SAMPLE, &
    DT,t0,OUTPUT_FILES, &
    seismograms,number_receiver_global,nrec_local, &
    it,seismo_current,seismo_offset, &
    myrank,WRITE_SEISMOGRAMS_BY_MAIN

  implicit none

  ! local parameters
  integer :: irec,irec_local,ier
  integer :: iorientation,isample,it_tmp
  real(kind=CUSTOM_REAL) :: time_t

  character(len=4) :: chn
  character(len=MAX_STRING_LEN) :: sisname
  character(len=2) :: bic

  ! for adjoint simulations, source locations become the "receivers" for storing seismograms

  ! safety check
  if (WRITE_SEISMOGRAMS_BY_MAIN) &
    call exit_MPI(myrank,'Error write_adj_seismograms() needs WRITE_SEISMOGRAMS_BY_MAIN turned off')

  ! checks if anything to do
  if (nrec_local <= 0 ) return

  ! get band code
  call band_instrument_code(DT,bic)

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)

    do iorientation = 1,9
      if (iorientation == 1) then
       chn = 'SNN'
      else if (iorientation == 2) then
       chn = 'SEE'
      else if (iorientation == 3) then
       chn = 'SZZ'
      else if (iorientation == 4) then
       chn = 'SNE'
      else if (iorientation == 5) then
       chn = 'SNZ'
      else if (iorientation == 6) then
       chn = 'SEZ'
      else if (iorientation == 7) then
       !chn = 'LHN'
       chn = bic(1:2)//'N'
      else if (iorientation == 8) then
       chn = bic(1:2)//'E'
      else if (iorientation == 9) then
       chn = bic(1:2)//'Z'
      endif

      ! create the name of the seismogram file for each slice
      ! file name includes the name of the station, the network and the component
      ! for example: NT.S000001.MXN.sem.ascii
      write(sisname,"(a3,'.',a1,i6.6,'.',a3,'.sem.ascii')") '/NT','S',irec,chn

      ! save seismograms in text format
      if (it <= NTSTEP_BETWEEN_OUTPUT_SEISMOS) then
        ! open new file
        open(unit=IOUT,file=trim(OUTPUT_FILES)//sisname(1:len_trim(sisname)), &
              status='unknown',action='write',iostat=ier)
      else
        ! for it > NTSTEP_BETWEEN_OUTPUT_SEISMOS
        ! append to existing file
        open(unit=IOUT,file=trim(OUTPUT_FILES)//sisname(1:len_trim(sisname)), &
              status='old',position='append',action='write',iostat=ier)
      endif
      if (ier /= 0) call exit_mpi(myrank,'Error opening file: '//trim(OUTPUT_FILES)//trim(sisname))

      ! make sure we never write more than the maximum number of time steps
      do isample = 1,seismo_current
        ! current time
        ! current time increment
        it_tmp = seismo_offset + isample

        ! subtract onset time to make sure travel time is correct
        ! distinguish between single and double precision for reals
        time_t = real(dble((it_tmp-1) * NTSTEP_BETWEEN_OUTPUT_SAMPLE) * DT - t0,kind=CUSTOM_REAL)

        ! output
        write(IOUT,*) time_t,seismograms(iorientation,irec_local,isample)
      enddo

      close(IOUT)
    enddo
  enddo

  end subroutine write_adj_seismograms

!
!-------------------------------------------------------------------------------------------------
!

! write strain seismograms to text files

  subroutine write_seismograms_strain()

  use constants, only: MAX_STRING_LEN,CUSTOM_REAL,IOUT,myrank

  use specfem_par, only: NSTEP,DT,t0,seismo_current,seismo_offset, &
    OUTPUT_FILES,WRITE_SEISMOGRAMS_BY_MAIN,SIMULATION_TYPE, &
    number_receiver_global,nrec_local,network_name,station_name, &
    seismograms_eps

  implicit none

  ! local parameters
  integer :: irec,irec_local,ier,it_tmp
  integer :: iorientation,isample
  double precision :: value
  double precision :: timeval

  character(len=3) :: chn
  character(len=MAX_STRING_LEN) :: sisname

  ! safety check
  if (WRITE_SEISMOGRAMS_BY_MAIN) &
    call exit_MPI(myrank,'Error write_seismograms_strain() needs WRITE_SEISMOGRAMS_BY_MAIN turned off')

  ! checks if anything to do
  if (nrec_local <= 0 ) return

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)

    ! strain orientation in local, radial direction (N-E-UP)
    do iorientation = 1,6

      select case (iorientation)
      case (1)
       chn = 'SNN'
      case (2)
       chn = 'SEE'
      case (3)
       chn = 'SZZ'
      case (4)
       chn = 'SNE'
      case (5)
       chn = 'SNZ'
      case (6)
       chn = 'SEZ'
      case default
        call exit_MPI(myrank,'incorrect channel value in write_seismograms_strain()')
      end select

      ! strain seismograms will have channel-code S**
      !
      ! create the name of the strain seismogram file using the station name and network name
      ! using format: **net**.**sta**.channel
      ! for example: IU.KONO.SNN.sem.ascii
      write(sisname,"('/',a,'.',a,'.',a3,'.sem.ascii')") trim(network_name(irec)),trim(station_name(irec)),chn

      ! save seismograms in text format
      if (seismo_offset == 0) then
        !open new file
        open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname),status='unknown',action='write',iostat=ier)
      else
        ! for it > NTSTEP_BETWEEN_OUTPUT_SEISMOS
        !append to existing file
        open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname),status='old',position='append',action='write',iostat=ier)
      endif
      if (ier /= 0) call exit_mpi(myrank,'Error opening file: '//trim(OUTPUT_FILES)//trim(sisname))

      ! subtract half duration of the source to make sure travel time is correct
      do isample = 1,seismo_current
        ! seismogram value
        value = dble(seismograms_eps(iorientation,irec_local,isample))

        ! current time increment
        it_tmp = seismo_offset + isample

        ! current time
        if (SIMULATION_TYPE == 3) then
          timeval = dble(NSTEP-it_tmp)*DT - t0
        else
          timeval = dble(it_tmp-1)*DT - t0
        endif

        ! writes out to file
        ! distinguish between single and double precision for reals
        write(IOUT,*) real(timeval, kind=CUSTOM_REAL), real(value, kind=CUSTOM_REAL)
      enddo
      close(IOUT)

    enddo
  enddo

  end subroutine write_seismograms_strain

!
!-------------------------------------------------------------------------------------------------
!

  subroutine band_instrument_code(DT,bic)

! This subroutine is to choose the appropriate band and instrument codes for channel names of seismograms
! based on the IRIS convention (first two letters of channel codes which were LH(Z/E/N) previously).
! For consistency with observed data, we now use the IRIS convention for band codes (first letter in channel codes)of
! SEM seismograms governed by their sampling rate.
! Instrument code (second letter in channel codes) is fixed to "X" which is assigned by IRIS for synthetic seismograms.
! See the manual for further explanations!
! Ebru, November 2010

  implicit none

  double precision,intent(in) :: DT
  character(len=2),intent(out) :: bic

  bic = ''

  if (1.0d0 <= DT)  bic = 'LX'
  if (0.1d0 < DT .and. DT < 1.0d0) bic = 'MX'
  if (0.0125d0 < DT .and. DT <= 0.1d0) bic = 'BX'
  if (0.004d0 < DT .and. DT <= 0.0125d0) bic = 'HX'
  if (0.001d0 < DT .and. DT <= 0.004d0) bic = 'CX'
  if (DT <= 0.001d0) bic = 'FX'

  end subroutine band_instrument_code

