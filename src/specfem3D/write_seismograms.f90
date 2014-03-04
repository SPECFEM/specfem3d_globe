!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            August 2013
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

  subroutine write_seismograms()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! note: this routine gets called if( nrec_local > 0 .or. ( WRITE_SEISMOGRAMS_BY_MASTER .and. myrank == 0 ) )

  ! update position in seismograms
  seismo_current = seismo_current + 1

  ! compute & store the seismograms only if there is at least one receiver located in this slice
  if (nrec_local > 0) then

    ! gets resulting array values onto CPU
    if( GPU_MODE ) then
      ! gets field values from GPU
      ! this transfers fields only in elements with stations for efficiency
      call write_seismograms_transfer_cuda(Mesh_pointer, &
                                displ_crust_mantle,b_displ_crust_mantle, &
                                eps_trace_over_3_crust_mantle, &
                                epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                                epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                                number_receiver_global, &
                                ispec_selected_rec,ispec_selected_source, &
                                ibool_crust_mantle)

      ! synchronizes field values from GPU
      if( GPU_ASYNC_COPY) then
        call transfer_seismo_from_device_async(Mesh_pointer, &
                                               displ_crust_mantle,b_displ_crust_mantle, &
                                               number_receiver_global,ispec_selected_rec,ispec_selected_source, &
                                               ibool_crust_mantle)
      endif

    endif

    ! computes traces at interpolated receiver locations
    select case( SIMULATION_TYPE )
    case( 1 )
      call compute_seismograms(NGLOB_CRUST_MANTLE,displ_crust_mantle, &
                               seismo_current,seismograms)
    case( 2 )
      call compute_seismograms_adjoint(displ_crust_mantle, &
                                       eps_trace_over_3_crust_mantle, &
                                       epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                                       epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                                       nit_written, &
                                       moment_der,sloc_der,stshift_der,shdur_der, &
                                       seismograms)
    case( 3 )
      call compute_seismograms(NGLOB_CRUST_MANTLE_ADJOINT,b_displ_crust_mantle, &
                               seismo_current,seismograms)
    end select

  endif ! nrec_local

  ! write the current or final seismograms
  if(seismo_current == NTSTEP_BETWEEN_OUTPUT_SEISMOS .or. it == it_end) then

    ! writes out seismogram files
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then

      call write_seismograms_to_file()

      ! user output
      if(myrank==0) then
        write(IMAIN,*)
        write(IMAIN,*) ' Total number of time steps written: ', it-it_begin+1
        write(IMAIN,*)
        call flush_IMAIN()
      endif
    else if( SIMULATION_TYPE == 2 ) then
      if( nrec_local > 0 ) &
        call write_adj_seismograms(nit_written)
      nit_written = it
    endif

    ! resets current seismogram position
    seismo_offset = seismo_offset + seismo_current
    seismo_current = 0

  endif

  end subroutine write_seismograms

!
!-------------------------------------------------------------------------------------------------
!

! write seismograms to files
  subroutine write_seismograms_to_file()

  use asdf_data
  use constants_solver
  use specfem_par,only: &
          NPROCTOT_VAL,myrank,nrec,nrec_local, &
          number_receiver_global,seismograms, &
          islice_selected_rec, &
          seismo_offset,seismo_current, &
          OUTPUT_SEISMOS_ASCII_TEXT, &
          OUTPUT_SEISMOS_ASDF, DT, &
          NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE, &
          OUTPUT_FILES, &
          WRITE_SEISMOGRAMS_BY_MASTER

  implicit none

  ! local parameters
  double precision :: write_time_begin,write_time
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: one_seismogram

  integer :: iproc,sender,irec_local,iorientation,irec,ier,receiver
  integer :: nrec_local_received
  integer :: total_seismos,total_seismos_local
  integer,dimension(:),allocatable:: islice_num_rec_local
  character(len=256) :: sisname
  ! timing
  double precision, external :: wtime
  ! todo: only needed for asdf output but I am passing this around
  type(asdf_event) :: asdf_container

  ! allocates single station seismogram
  allocate(one_seismogram(NDIM,NTSTEP_BETWEEN_OUTPUT_SEISMOS),stat=ier)
  if(ier /= 0) call exit_mpi(myrank,'error while allocating one temporary seismogram')

  ! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

  ! writes out seismograms
  if(.not. WRITE_SEISMOGRAMS_BY_MASTER) then

    ! all the processes write their local seismograms themselves

    write_time_begin = wtime()

    if(OUTPUT_SEISMOS_ASCII_TEXT .and. SAVE_ALL_SEISMOS_IN_ONE_FILE) then
      write(sisname,'(A,I5.5)') '/all_seismograms_node_',myrank

      if(USE_BINARY_FOR_LARGE_FILE) then
        if (seismo_offset==0) then
          open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.bin',status='unknown',form='unformatted',action='write')
        else
          open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.bin',status='old',&
               form='unformatted',position='append',action='write')
        endif
      else
        if (seismo_offset==0) then
          open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.ascii',status='unknown',form='formatted',action='write')
        else
          open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.ascii',status='old',&
               form='formatted',position='append',action='write')
        endif
      endif
    endif

   ! initializes the asdf data structure by allocating arrays
   if (OUTPUT_SEISMOS_ASDF) then
      total_seismos_local = 0
      do irec_local = 1, nrec_local
        do iorientation = 1, 3
          total_seismos_local=total_seismos_local+1
        enddo
      enddo
      call init_asdf_data(asdf_container, total_seismos_local)
    endif

    total_seismos_local = 0
    ! loop on all the local receivers
    do irec_local = 1,nrec_local

      ! get global number of that receiver
      irec = number_receiver_global(irec_local)

      total_seismos_local = total_seismos_local + 1

      one_seismogram = seismograms(:,irec_local,:)

      ! write this seismogram
      ! asdf data structure is passed as an argument
      call write_one_seismogram(one_seismogram,irec,irec_local,asdf_container)
    enddo

    ! write asdf container to the file
    if (OUTPUT_SEISMOS_ASDF) then
      call synchronize_all()
      call write_asdf(asdf_container)
      ! deallocate the contanier
      call close_asdf_data(asdf_container, total_seismos_local)
    endif

    ! create one large file instead of one small file per station to avoid file system overload
    if(OUTPUT_SEISMOS_ASCII_TEXT .and. SAVE_ALL_SEISMOS_IN_ONE_FILE) close(IOUT)

    !if(total_seismos_local/= nrec_local) call exit_MPI(myrank,'incorrect total number of receivers saved')

    ! user output
    if(myrank == 0) then
      write_time = wtime() - write_time_begin
      write(IMAIN,*)
      write(IMAIN,*) 'Writing the seismograms in parallel took ',write_time,' seconds'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

  else ! WRITE_SEISMOGRAMS_BY_MASTER

    ! now only the master process does the writing of seismograms and
    ! collects the data from all other processes

    write_time_begin = wtime()

    if(myrank == 0) then ! on the master, gather all the seismograms

      ! create one large file instead of one small file per station to avoid file system overload
      if(OUTPUT_SEISMOS_ASCII_TEXT .and. SAVE_ALL_SEISMOS_IN_ONE_FILE) then
         write(sisname,'(A)') '/all_seismograms'

       if(USE_BINARY_FOR_LARGE_FILE) then
         if (seismo_offset==0) then
           open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.bin',status='unknown',form='unformatted',action='write')
         else
           open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.bin',status='old',&
                form='unformatted',position='append',action='write')
         endif
       else
         if (seismo_offset==0) then
           open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.ascii',status='unknown',form='formatted',action='write')
         else
           open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.ascii',status='old',&
                form='formatted',position='append',action='write')
         endif
       endif

      endif

      ! counts number of local receivers for each slice
      allocate(islice_num_rec_local(0:NPROCTOT_VAL-1),stat=ier)
      if( ier /= 0 ) call exit_mpi(myrank,'error allocating islice_num_rec_local')

      islice_num_rec_local(:) = 0
      do irec = 1,nrec
        iproc = islice_selected_rec(irec)
        ! checks iproc value
        if( iproc < 0 .or. iproc >= NPROCTOT_VAL ) then
          print*,'error :',myrank,'iproc = ',iproc,'NPROCTOT = ',NPROCTOT_VAL
          call exit_mpi(myrank,'error iproc in islice_selected_rec')
        endif
        ! sums number of receivers for each slice
        islice_num_rec_local(iproc) = islice_num_rec_local(iproc) + 1
      enddo

      total_seismos = 0

      ! loop on all the slices
      do iproc = 0,NPROCTOT_VAL-1

       ! communicates only with processes which contain local receivers (to minimize mpi chatter)
       if( islice_num_rec_local(iproc) == 0 ) cycle

       ! receive except from proc 0, which is me and therefore I already have this value
       sender = iproc
       if(iproc == 0) then
         ! master is current slice
         nrec_local_received = nrec_local
       else
         ! receives info from slave processes
         call recv_singlei(nrec_local_received,sender,itag)
         if(nrec_local_received <= 0) call exit_MPI(myrank,'error while receiving local number of receivers')
       endif
       if (nrec_local_received > 0) then
         do irec_local = 1,nrec_local_received
           ! receive except from proc 0, which is myself and therefore I already have these values
           if(iproc == 0) then
             ! get global number of that receiver
             irec = number_receiver_global(irec_local)
             one_seismogram(:,:) = seismograms(:,irec_local,:)
           else
             ! receives info from slave processes
             call recv_singlei(irec,sender,itag)
             if(irec < 1 .or. irec > nrec) call exit_MPI(myrank,'error while receiving global receiver number')
             call recv_cr(one_seismogram,NDIM*seismo_current,sender,itag)
           endif

           total_seismos = total_seismos + 1
           ! write this seismogram
           call write_one_seismogram(one_seismogram,irec,irec_local)

         enddo
       endif
      enddo
      deallocate(islice_num_rec_local)

      write(IMAIN,*)
      write(IMAIN,*) 'Total number of receivers saved is ',total_seismos,' out of ',nrec
      write(IMAIN,*)

      if(total_seismos /= nrec) call exit_MPI(myrank,'incorrect total number of receivers saved')

      ! create one large file instead of one small file per station to avoid file system overload
      if(SAVE_ALL_SEISMOS_IN_ONE_FILE) close(IOUT)

    else  ! on the nodes, send the seismograms to the master
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


    if(myrank == 0) then
      write_time  = wtime() - write_time_begin
      write(IMAIN,*)
      write(IMAIN,*) 'Writing the seismograms by master proc alone took ',write_time,' seconds'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

  endif ! WRITE_SEISMOGRAMS_BY_MASTER

  deallocate(one_seismogram)

  end subroutine write_seismograms_to_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_one_seismogram(one_seismogram,irec,irec_local,asdf_container)

  use asdf_data
  use constants_solver
  use specfem_par,only: &
          myrank, &
          station_name,network_name,stlat,stlon, &
          DT, &
          seismo_current, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_ASDF,&
          OUTPUT_SEISMOS_SAC_BINARY,ROTATE_SEISMOGRAMS_RT,NTSTEP_BETWEEN_OUTPUT_SEISMOS

  use specfem_par,only: &
          cmt_lat=>cmt_lat_SAC,cmt_lon=>cmt_lon_SAC

  implicit none

  integer :: irec,irec_local
  real(kind=CUSTOM_REAL), dimension(NDIM,NTSTEP_BETWEEN_OUTPUT_SEISMOS) :: one_seismogram

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(5,NTSTEP_BETWEEN_OUTPUT_SEISMOS) :: seismogram_tmp
  integer :: iorientation,length_station_name,length_network_name

  character(len=4) :: chn
  character(len=256) :: sisname,sisname_big_file
  character(len=2) :: bic

  ! variables used for calculation of backazimuth and
  ! rotation of components if ROTATE_SEISMOGRAMS=.true.
  integer :: ior_start,ior_end
  double precision :: backaz
  double precision :: phi
  real(kind=CUSTOM_REAL) :: cphi,sphi
  integer :: isample
  type(asdf_event) :: asdf_container

  ! initializes
  seismogram_tmp(:,:) = 0.0_CUSTOM_REAL

  ! get band code
  call band_instrument_code(DT,bic)

  if (ROTATE_SEISMOGRAMS_RT) then ! iorientation 1=N,2=E,3=Z,4=R,5=T
    ior_start=3    ! starting from Z
    ior_end  =5    ! ending with T => ZRT
  else
    ior_start=1    ! starting from N
    ior_end  =3    ! ending with Z => NEZ
  endif

  do iorientation = ior_start,ior_end      ! BS BS changed according to ROTATE_SEISMOGRAMS_RT

    select case( iorientation )
    case( 1 )
      !chn = 'LHN'
      chn = bic(1:2)//'N'
    case( 2 )
      !chn = 'LHE'
      chn = bic(1:2)//'E'
    case( 3 )
      !chn = 'LHZ'
      chn = bic(1:2)//'Z'
    case( 4 )
      !chn = 'LHR'
      chn = bic(1:2)//'R'
    case( 5 )
      !chn = 'LHT'
      chn = bic(1:2)//'T'
    case default
      call exit_MPI(myrank,'incorrect channel value')
    end select

    if (iorientation == 4 .or. iorientation == 5) then        ! LMU BS BS

      ! BS BS calculate backazimuth needed to rotate East and North
      ! components to Radial and Transverse components
      !  call get_backazimuth(elat,elon,stlat(irec),stlon(irec),backaz)
      call get_backazimuth(cmt_lat,cmt_lon,stlat(irec),stlon(irec),backaz)

      phi = backaz
      if (phi>180.d0) then
         phi = phi-180.d0
      else if (phi<180.d0) then
         phi = phi+180.d0
      else if (phi==180.d0) then
         phi = backaz
      endif

      cphi=cos(phi*DEGREES_TO_RADIANS)
      sphi=sin(phi*DEGREES_TO_RADIANS)

      ! BS BS do the rotation of the components and put result in
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

    else ! keep NEZ components
      do isample = 1,seismo_current
        seismogram_tmp(iorientation,isample) = one_seismogram(iorientation,isample)
      enddo

    endif

    ! create the name of the seismogram file for each slice
    ! file name includes the name of the station and the network
    length_station_name = len_trim(station_name(irec))
    length_network_name = len_trim(network_name(irec))

    ! check that length conforms to standard
    if(length_station_name < 1 .or. length_station_name > MAX_LENGTH_STATION_NAME) &
           call exit_MPI(myrank,'wrong length of station name')

    if(length_network_name < 1 .or. length_network_name > MAX_LENGTH_NETWORK_NAME) &
           call exit_MPI(myrank,'wrong length of network name')

    ! create the name of the seismogram file using the station name and network name
    write(sisname,"('/',a,'.',a,'.',a3,'.sem')") station_name(irec)(1:length_station_name), &
                   network_name(irec)(1:length_network_name),chn

    ! create this name also for the text line added to the unique big seismogram file
    write(sisname_big_file,"(a,'.',a,'.',a3,'.sem')") station_name(irec)(1:length_station_name), &
                   network_name(irec)(1:length_network_name),chn

    ! SAC output format
    if( OUTPUT_SEISMOS_SAC_ALPHANUM .or. OUTPUT_SEISMOS_SAC_BINARY ) &
      call write_output_SAC(seismogram_tmp,irec,iorientation,sisname,chn,phi)

    ! ASCII output format
    if(OUTPUT_SEISMOS_ASCII_TEXT) &
      call write_output_ASCII(seismogram_tmp,iorientation,sisname,sisname_big_file)

    ! ASDF output format
    if(OUTPUT_SEISMOS_ASDF) &
      call store_asdf_data(asdf_container,seismogram_tmp,irec_local,irec,chn,iorientation)

  enddo ! do iorientation

  end subroutine write_one_seismogram

!
!-------------------------------------------------------------------------------------------------
!

! write adjoint seismograms to text files

  subroutine write_adj_seismograms(nit_written)

  use constants
  use specfem_par,only: NSTEP,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
    DT,t0,LOCAL_TMP_PATH, &
    seismograms,number_receiver_global,nrec_local, &
    it

  implicit none

  integer :: nit_written

  ! local parameters
  integer :: irec,irec_local
  integer :: iorientation,isample

  character(len=4) :: chn
  character(len=150) :: clean_LOCAL_PATH,final_LOCAL_PATH
  character(len=256) :: sisname
  character(len=2) :: bic

  call band_instrument_code(DT,bic)

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)

    do iorientation = 1,9
      if(iorientation == 1) then
       chn = 'SNN'
      else if(iorientation == 2) then
       chn = 'SEE'
      else if(iorientation == 3) then
       chn = 'SZZ'
      else if(iorientation == 4) then
       chn = 'SNE'
      else if(iorientation == 5) then
       chn = 'SNZ'
      else if(iorientation == 6) then
       chn = 'SEZ'
      else if(iorientation == 7) then
       !chn = 'LHN'
       chn = bic(1:2)//'N'
      else if(iorientation == 8) then
       chn = bic(1:2)//'E'
      else if(iorientation == 9) then
       chn = bic(1:2)//'Z'
      endif


      ! create the name of the seismogram file for each slice
      ! file name includes the name of the station, the network and the component
      write(sisname,"(a,i6.6,'.',a,'.',a3,'.sem')") 'S',irec,'NT',chn

      ! suppress white spaces if any
      clean_LOCAL_PATH = adjustl(LOCAL_TMP_PATH)

      ! create full final local path
      final_LOCAL_PATH = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) // '/'

      ! save seismograms in text format with no subsampling.
      ! Because we do not subsample the output, this can result in large files
      ! if the simulation uses many time steps. However, subsampling the output
      ! here would result in a loss of accuracy when one later convolves
      ! the results with the source time function
      if(it <= NTSTEP_BETWEEN_OUTPUT_SEISMOS) then
        !open new file
        open(unit=IOUT,file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)),&
              status='unknown',action='write')
      else if(it > NTSTEP_BETWEEN_OUTPUT_SEISMOS) then
        !append to existing file
        open(unit=IOUT,file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)),&
              status='old',position='append',action='write')
      endif
      ! make sure we never write more than the maximum number of time steps
      ! subtract half duration of the source to make sure travel time is correct
      do isample = nit_written+1,min(it,NSTEP)
        ! distinguish between single and double precision for reals
        if(CUSTOM_REAL == SIZE_REAL) then
          write(IOUT,*) sngl(dble(isample-1)*DT - t0),' ',seismograms(iorientation,irec_local,isample-nit_written)
        else
          write(IOUT,*) dble(isample-1)*DT - t0,' ',seismograms(iorientation,irec_local,isample-nit_written)
        endif
      enddo

      close(IOUT)
    enddo
  enddo

  end subroutine write_adj_seismograms

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
  double precision :: DT
  character(len=2) :: bic

  if (DT >= 1.0d0)  bic = 'LX'
  if (DT < 1.0d0 .and. DT > 0.1d0) bic = 'MX'
  if (DT <= 0.1d0 .and. DT > 0.0125d0) bic = 'BX'
  if (DT <= 0.0125d0 .and. DT > 0.004d0) bic = 'HX'
  if (DT <= 0.004d0 .and. DT > 0.001d0) bic = 'CX'
  if (DT <= 0.001d0) bic = 'FX'

 end subroutine band_instrument_code
