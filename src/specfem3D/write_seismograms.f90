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

  use constants_solver, only: NGLOB_CRUST_MANTLE,NGLOB_CRUST_MANTLE_ADJOINT,FULL_GRAVITY_VAL,ROTATION_VAL

  use specfem_par, only: myrank,Mesh_pointer,GPU_MODE,GPU_ASYNC_COPY,SIMULATION_TYPE, &
    nrec_local,number_receiver_global,ispec_selected_rec,ispec_selected_source, &
    it,it_end, &
    seismo_current,seismo_offset, &
    seismograms, &
    nlength_seismogram, &
    NTSTEP_BETWEEN_OUTPUT_SEISMOS,NTSTEP_BETWEEN_OUTPUT_SAMPLE, &
    do_save_seismograms, &
    WRITE_SEISMOGRAMS_BY_MAIN,OUTPUT_SEISMOS_ASDF, &
    SAVE_SEISMOGRAMS_STRAIN, &
    moment_der,sloc_der,shdur_der,stshift_der, &
    scale_displ

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
          call compute_seismograms_gpu(Mesh_pointer,seismograms,seismo_current,it,it_end,scale_displ,nlength_seismogram)
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
          if (OUTPUT_ADJOINT_WAVEFIELD_SEISMOGRAMS) then
            ! uncomment to output adjoint wavefield instead for seismogram output
            call compute_seismograms(NGLOB_CRUST_MANTLE_ADJOINT,displ_crust_mantle,seismo_current,seismograms)
          else
            ! default, backward reconstructed wavefield seismos
            call compute_seismograms(NGLOB_CRUST_MANTLE_ADJOINT,b_displ_crust_mantle,seismo_current,seismograms)
          endif
        else
          ! on GPU
          call compute_seismograms_gpu(Mesh_pointer,seismograms,seismo_current,it,it_end,scale_displ,nlength_seismogram)
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

      ! full gravity seismograms
      if (FULL_GRAVITY_VAL) call SIEM_compute_seismos()

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
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'Writing the seismograms'
      call flush_IMAIN()
    endif

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
        if (do_save_seismograms) then
          ! displacement (seismograms)
          call write_seismograms_to_file(1)
          if (FULL_GRAVITY_VAL) then
            call write_seismograms_to_file(5) ! seismograms_phi
            call write_seismograms_to_file(6) ! seismograms_pgrav
            call write_seismograms_to_file(7) ! seismograms_grav
            if (ROTATION_VAL) call write_seismograms_to_file(8) ! seismograms_corio
          endif
        endif
        if (SAVE_SEISMOGRAMS_STRAIN) then
          ! strain
          call write_seismograms_strain(1)
          if (FULL_GRAVITY_VAL) call write_seismograms_strain(2)  ! gravity strain
        endif
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
  subroutine write_seismograms_to_file(istore)

  use constants, only: MAX_STRING_LEN,CUSTOM_REAL,NDIM,IMAIN,IOUT,itag

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

  integer,intent(in) :: istore

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: one_seismogram

  integer :: iproc,sender,irec_local,irec,ier,receiver
  integer :: nrec_local_received
  integer :: total_seismos,length_network_name,length_station_name

  character(len=MAX_STRING_LEN) :: sisname,staname
  character(len=8) :: component

  ! saves displacement or potential
  component = ''
  select case (istore)
  case (1)
    !component = 'd'   ! displacement (seismograms array holds displacement)
    ! for backward compatibility with older globe versions:
    ! we leave this blank to ignore adding a type indicator to the seismogram names
    component = ''
  case (2)
    component = 'v'   ! velocity - not used yet...
  case (3)
    component = 'a'   ! acceleration - not used yet...
  case (4)
    component = 'p'   ! pressure - not used yet...
  case (5)
    component = 'G'   ! gravity potential (seismograms_phi)
  case (6)
    component = 'C.PGRAV'  ! gravity perturbation (seismograms_pgrav)
  case (7)
    component = 'C.GRAV'  ! free-air gravity change (seismograms_grav)
  case (8)
    component = 'C.CORIO'  ! Coriolis acceleration (seismograms_corio)
  case default
    call exit_MPI(myrank,'wrong component to save for seismograms')
  end select

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Component: .sem '//trim(component)
    call flush_IMAIN()
  endif

  ! allocates single station seismogram
  allocate(one_seismogram(NDIM,nlength_seismogram),stat=ier)
  if (ier /= 0) call exit_mpi(myrank,'Error while allocating one temporary seismogram')
  one_seismogram(:,:) = 0.0_CUSTOM_REAL

  ! ASDF format
  if (OUTPUT_SEISMOS_ASDF .and. istore == 1) then     ! only supported for seismogram output
    ! The writing of seismograms by the main proc is handled within write_asdf()
    do irec_local = 1,nrec_local

      ! get global number of that receiver
      irec = number_receiver_global(irec_local)

      one_seismogram(:,:) = seismograms(:,irec_local,:)

      ! write this seismogram
      ! note: ASDF data structure is given in module
      !       stores all traces into ASDF container in case
      call write_one_seismogram(one_seismogram,irec,irec_local,.true.,component,istore)
    enddo

    ! writes ASDF file output
    call write_asdf()

    call synchronize_all()

    ! deallocate the container
    call close_asdf_data()
  endif

  ! write 3D seismogram array
  if (OUTPUT_SEISMOS_3D_ARRAY .and. istore == 1) then   ! only supported for seismogram output
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

  else
    if (NTSTEP_BETWEEN_OUTPUT_SAMPLE > 1 .and. istore == 1) then
      ! save original DT and NSTEP for adjoint simulation
      if (myrank == 0 .and. seismo_offset == 0) then
        open(unit=IOUT,file=trim(OUTPUT_FILES)//'seismogram_stats.txt',status='unknown',form='formatted',action='write')
        write(IOUT,*) 'DT0     =', DT
        write(IOUT,*) 'NSTEP0  =', NSTEP
        close(IOUT)
      endif
    endif
  endif

  ! ASCII / SAC format
  if (OUTPUT_SEISMOS_ASCII_TEXT .or. OUTPUT_SEISMOS_SAC_ALPHANUM .or. OUTPUT_SEISMOS_SAC_BINARY) then

    ! write out seismograms: all processes write their local seismograms themselves
    if (.not. WRITE_SEISMOGRAMS_BY_MAIN) then

      ! all the processes write their local seismograms themselves
      if (SAVE_ALL_SEISMOS_IN_ONE_FILE .and. OUTPUT_SEISMOS_ASCII_TEXT) then
        write(sisname,'(A,I5.5)') '/all_seismograms_'//trim(component)//'_node_',myrank

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
        ! gets trace
        call get_single_trace(istore,irec_local,one_seismogram)

        ! get global number of that receiver
        irec = number_receiver_global(irec_local)

        ! write this seismogram
        ! note: ASDF data structure is given in module
        !       stores all traces into ASDF container in case
        call write_one_seismogram(one_seismogram,irec,irec_local,.false.,component,istore)
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
          write(sisname,'(A)') '/all_seismograms_'//trim(component)//'_main'

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
          if (islice_num_rec_local(iproc) == 0) cycle

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
                ! gets trace
                call get_single_trace(istore,irec_local,one_seismogram)
              else
                ! receives info from secondary processes
                call recv_singlei(irec,sender,itag)
                if (irec < 1 .or. irec > nrec) call exit_MPI(myrank,'Error while receiving global receiver number')
                ! gets trace from secondary processes
                call recv_cr(one_seismogram,NDIM*seismo_current,sender,itag)
              endif

              ! write this seismogram
              call write_one_seismogram(one_seismogram,irec,irec_local,.false.,component,istore)

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
            ! sends receiver number
            call send_singlei(irec,receiver,itag)

            ! gets trace
            call get_single_trace(istore,irec_local,one_seismogram)
            ! sends traces
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

contains

  subroutine get_single_trace(istore,irec_local,one_seismogram)

  use constants, only: NDIM,CUSTOM_REAL
  use specfem_par, only: seismograms,nlength_seismogram,seismo_current
  use specfem_par_full_gravity, only: seismograms_phi,seismograms_pgrav,seismograms_grav,seismograms_corio

  implicit none

  integer,intent(in) :: istore,irec_local
  real(kind=CUSTOM_REAL),intent(out) :: one_seismogram(NDIM,nlength_seismogram)

  ! local parameters
  integer :: i

  ! init trace
  one_seismogram(:,:) = 0.0_CUSTOM_REAL

  select case (istore)
  case (1)
    ! displacement
    !one_seismogram(:,:) = seismograms(:,irec_local,:)
    do i = 1,seismo_current
      one_seismogram(:,i) = seismograms(:,irec_local,i)
    enddo
  case (2)
    ! velocity
    !do i = 1,seismo_current
    !  one_seismogram(:,i) = seismograms_v(:,irec_local,i)
    !enddo
    ! not used yet...
    continue
  case (3)
    ! acceleration
    !do i = 1,seismo_current
    !  one_seismogram(:,i) = seismograms_a(:,irec_local,i)
    !enddo
    ! not used yet...
    continue
  case (4)
    ! pressure
    !do i = 1,seismo_current
    !  one_seismogram(1,i) = seismograms_p(1,irec_local,i) ! single component
    !enddo
    ! not used yet...
    continue
  case (5)
    ! gravity potential (seismograms_phi)
    do i = 1,seismo_current
      one_seismogram(1,i) = seismograms_phi(1,irec_local,i)  ! potential (single component)
    enddo
  case (6)
    ! gravity perturbation (seismograms_pgrav)
    do i = 1,seismo_current
      one_seismogram(:,i) = seismograms_pgrav(:,irec_local,i)
    enddo
  case (7)
    ! free-air gravity (seismograms_grav)
    do i = 1,seismo_current
      one_seismogram(:,i) = seismograms_grav(:,irec_local,i)
    enddo
  case (8)
    ! Coriolis acceleration (seismograms_corio)
    do i = 1,seismo_current
      one_seismogram(:,i) = seismograms_corio(:,irec_local,i)
    enddo
  case default
    call exit_MPI(myrank,'wrong istore component to get single trace')
  end select

  end subroutine get_single_trace

  end subroutine write_seismograms_to_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_one_seismogram(one_seismogram,irec,irec_local,is_for_asdf,component,istore)

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

  integer,intent(in) :: istore
  character(len=8),intent(in) :: component

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

  if (ROTATE_SEISMOGRAMS_RT) then ! iorientation 1 = N,2 = E,3 = Z,4 = R,5 = T
    ior_start = 3    ! starting from Z
    ior_end   = 5    ! ending with T => ZRT
  else
    ior_start = 1    ! starting from N
    ior_end   = 3    ! ending with Z => NEZ
  endif

  ! single component only for pressure & gravitational potential
  if (istore == 4 .or. istore == 5) then
    ior_start = 1
    ior_end = 1
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

    ! single component only for pressure & gravitational potential
    if (istore == 4 .or. istore == 5) then
      chn = bic(1:2)//'P'
    endif

    ! backazimuth rotation
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

      cphi = real(cos(phi*DEGREES_TO_RADIANS),kind=CUSTOM_REAL)
      sphi = real(sin(phi*DEGREES_TO_RADIANS),kind=CUSTOM_REAL)

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

    ! full gravity seismos add an additional component indicator to the name
    if (istore > 4) then
      ! using format: **net**.**sta**.channel.C.GRAV.sem.ascii
      write(sisname,"('/',a,'.',a,'.',a3,'.',a,'.sem')") network_name(irec)(1:length_network_name), &
                     station_name(irec)(1:length_station_name),chn,trim(component)

      ! create this name also for the text line added to the unique big seismogram file
      write(sisname_big_file,"(a,'.',a,'.',a3,'.',a,'.sem')") network_name(irec)(1:length_network_name), &
                     station_name(irec)(1:length_station_name),chn,trim(component)
    endif


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

  subroutine write_seismograms_strain(istore)

  use constants, only: MAX_STRING_LEN,CUSTOM_REAL,IOUT,IMAIN,myrank,NDIM,DEGREES_TO_RADIANS

  use specfem_par, only: NSTEP,DT,t0,seismo_current,seismo_offset, &
    OUTPUT_FILES,WRITE_SEISMOGRAMS_BY_MAIN,ROTATE_SEISMOGRAMS_RT,SIMULATION_TYPE, &
    number_receiver_global,nrec_local,network_name,station_name, &
    stlat,stlon, &
    seismograms_eps

  use specfem_par_full_gravity, only: seismograms_Hgrav

  use specfem_par, only: &
    cmt_lat => cmt_lat_SAC,cmt_lon => cmt_lon_SAC

  implicit none

  integer,intent(in) :: istore

  ! local parameters
  integer :: irec,irec_local,ier,it_tmp
  integer :: iorientation,isample
  double precision :: value
  double precision :: timeval

  character(len=3) :: chn
  character(len=MAX_STRING_LEN) :: sisname

  ! full gravity strain
  integer :: idim,jdim
  integer :: ior_start,ior_end
  real(kind=CUSTOM_REAL),dimension(NDIM,NDIM) :: one_matrix
  real(kind=CUSTOM_REAL),dimension(NDIM,NDIM) :: Qmat,QmatT
  double precision :: backaz,phi
  real(kind=CUSTOM_REAL) :: cphi,sphi

  ! safety check
  if (WRITE_SEISMOGRAMS_BY_MAIN) &
    call exit_MPI(myrank,'Error write_seismograms_strain() needs WRITE_SEISMOGRAMS_BY_MAIN turned off')

  ! user output
  if (myrank == 0) then
    if (istore == 1) &
      write(IMAIN,*) 'Strain        : .S**'
    if (istore == 2) &
      write(IMAIN,*) 'Gravity strain: .H**.PGRAV'
    call flush_IMAIN()
  endif

  ! checks if anything to do
  if (nrec_local <= 0 ) return

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)

    ! gravity strain rotation
    if (istore == 2 .and. ROTATE_SEISMOGRAMS_RT) then
      ! edit the channel names:
      ! orientation 1=N,2=E,3=Z
      ! note here that in write_seismograms the order ZRT is used, but because the original matrix is defined in NEZ order,
      ! which we will rotate, its easier to use RTZ here
      call get_backazimuth(cmt_lat,cmt_lon,stlat(irec),stlon(irec),backaz)

      ! rotation angle phi takes opposite direction; to have radial direction pointing in outgoing direction
      phi = backaz
      if (phi > 180.d0) then
         phi = phi-180.d0
      else if (phi < 180.d0) then
         phi = phi+180.d0
      else if (phi == 180.d0) then
         phi = 0.d0
      endif

      cphi = real(cos(phi*DEGREES_TO_RADIANS),kind=CUSTOM_REAL)
      sphi = real(sin(phi*DEGREES_TO_RADIANS),kind=CUSTOM_REAL)

      ! Create rotation matrix and transpose:
      Qmat(1,1) = cphi
      Qmat(1,2) = sphi
      Qmat(1,3) = 0.0_CUSTOM_REAL
      Qmat(2,1) = -sphi
      Qmat(2,2) = cphi
      Qmat(2,3) = 0.0_CUSTOM_REAL
      Qmat(3,1) = 0.0_CUSTOM_REAL
      Qmat(3,2) = 0.0_CUSTOM_REAL
      Qmat(3,3) = 1.0_CUSTOM_REAL

      QmatT = transpose(Qmat)

      ! Rotate the matrix for each timestep Q M Q^T :
      do isample = 1,seismo_current
        ! gets strain matrix
        one_matrix(:,:) = seismograms_Hgrav(:,:,irec_local,isample)
        ! rotates
        one_matrix(:,:) = matmul(Qmat, one_matrix)
        one_matrix(:,:) = matmul(one_matrix, QmatT)
        ! stores
        seismograms_Hgrav(:,:,irec_local,isample) = one_matrix(:,:)
      enddo
    endif

    ! strain orientation in local, radial direction (N-E-UP)
    ! default - symmetric strain tensor (w/ 6 components)
    ior_start = 1
    ior_end = 6

    ! asymmetric gravity strain tensor (all 3x3 components might be different)
    if (istore == 2) ior_end = 9

    do iorientation = ior_start,ior_end

      ! determines channel name
      chn = ''
      select case (istore)
      case (1)
        ! default
        ! elastic strain
        ! see compute_seismograms_strain() routine for orientations
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

      case (2)
        ! gravity strain
        ! see SIEM_compute_seismograms_pgrav() routine for orientations
        ! global -> local (n-e-up)
        ! eps_xx -> eps_nn
        ! eps_xy -> eps_ne
        ! eps_xz -> eps_nz
        ! eps_yx -> eps_en
        ! eps_yy -> eps_ee
        ! eps_yz -> eps_ez
        ! eps_zx -> eps_zn
        ! eps_zy -> eps_ze
        ! eps_zz -> eps_zz (z in radial direction up)
        select case (iorientation)
        case (1)
         chn = 'HNN'
         if (ROTATE_SEISMOGRAMS_RT) chn = 'HRR'
        case (2)
         chn = 'HNE'
         if (ROTATE_SEISMOGRAMS_RT) chn = 'HRT'
        case (3)
         chn = 'HNZ'
         if (ROTATE_SEISMOGRAMS_RT) chn = 'HRZ'
        case (4)
         chn = 'HEN'
         if (ROTATE_SEISMOGRAMS_RT) chn = 'HTR'
        case (5)
         chn = 'HEE'
         if (ROTATE_SEISMOGRAMS_RT) chn = 'HTT'
        case (6)
         chn = 'HEZ'
         if (ROTATE_SEISMOGRAMS_RT) chn = 'HTZ'
        case (7)
         chn = 'HZN'
         if (ROTATE_SEISMOGRAMS_RT) chn = 'HZR'
        case (8)
         chn = 'HZE'
         if (ROTATE_SEISMOGRAMS_RT) chn = 'HZT'
        case (9)
         chn = 'HZZ'
         if (ROTATE_SEISMOGRAMS_RT) chn = 'HZZ'
        case default
          call exit_MPI(myrank,'incorrect channel value in write_seismograms_strain()')
        end select

      case default
        call exit_MPI(myrank,'wrong component to save for strains')
      end select

      ! strain seismograms will have channel-code S**
      !
      ! create the name of the strain seismogram file using the station name and network name
      ! using format: **net**.**sta**.channel
      ! for example: IU.KONO.SNN.sem.ascii
      write(sisname,"('/',a,'.',a,'.',a3,'.sem.ascii')") trim(network_name(irec)),trim(station_name(irec)),chn

      ! full gravity strain has channel-code H**
      if (istore == 2) then
        ! for example: IU.KONO.HNN.PGRAV.sem.ascii
        write(sisname,"('/',a,'.',a,'.',a3,'.PGRAV','.sem.ascii')") trim(network_name(irec)),trim(station_name(irec)),chn
      endif

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

      ! make sure we never write more than the maximum number of time steps
      do isample = 1,seismo_current
        ! initializes
        value = 0.d0

        ! gets strain value
        select case (istore)
        case (1)
          value = dble(seismograms_eps(iorientation,irec_local,isample))
        case (2)
          ! example: iorientation == 1 -> idim = 1, jdim = 1
          !          iorientation == 2 -> idim = 1, jdim = 2
          !          iorientation == 3 -> idim = 1, jdim = 3
          !          iorientation == 4 -> idim = 2, jdim = 1
          !          ..
          idim = int((iorientation-1)/3) + 1
          jdim = iorientation - (idim-1) * 3
          value = dble(seismograms_Hgrav(idim,jdim,irec_local,isample))
        end select

        ! current time increment
        it_tmp = seismo_offset + isample

        ! subtract half duration of the source to make sure travel time is correct
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

