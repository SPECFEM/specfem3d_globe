!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2011
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

! write seismograms to files
  subroutine write_seismograms(myrank,seismograms,number_receiver_global,station_name, &
            network_name,stlat,stlon,stele,stbur, &
            nrec,nrec_local,ANGULAR_WIDTH_XI_IN_DEGREES,NEX_XI,DT,hdur,it_end, &
            yr,jda,ho,mi,sec,tshift_cmt,t_shift, &
            elat,elon,depth,event_name,cmt_lat,cmt_lon, &
            cmt_depth,cmt_hdur,NPROCTOT, &
            OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM, &
            OUTPUT_SEISMOS_SAC_BINARY,ROTATE_SEISMOGRAMS_RT,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
            seismo_offset,seismo_current,WRITE_SEISMOGRAMS_BY_MASTER,&
            SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE)

 implicit none

! standard include of the MPI library
 include 'mpif.h'

 include "constants.h"
 include "precision.h"

! parameters
 integer nrec,nrec_local,myrank,it_end,NPROCTOT,NEX_XI !,NSOURCES
 character(len=256) sisname

 integer :: seismo_offset, seismo_current, NTSTEP_BETWEEN_OUTPUT_SEISMOS
 integer, dimension(nrec_local) :: number_receiver_global

 real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS) :: seismograms
 double precision hdur,DT,ANGULAR_WIDTH_XI_IN_DEGREES

 character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
 character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name
 double precision tshift_cmt,t_shift,elat,elon,depth
 double precision cmt_lat,cmt_lon,cmt_depth,cmt_hdur
 double precision, dimension(nrec) :: stlat,stlon,stele,stbur
 integer yr,jda,ho,mi
 double precision sec
 !real mb
! character(len=12) ename
 character(len=20) event_name

! variables
 integer :: iproc,sender,irec_local,irec,ier,receiver,nrec_local_received,nrec_tot_found
 integer :: total_seismos,total_seismos_local
 double precision :: write_time_begin,write_time

 real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: one_seismogram

 integer msg_status(MPI_STATUS_SIZE)

 character(len=150) OUTPUT_FILES

! new flags to decide on seismogram type BS BS 06/2007
  logical OUTPUT_SEISMOS_ASCII_TEXT, OUTPUT_SEISMOS_SAC_ALPHANUM, &
          OUTPUT_SEISMOS_SAC_BINARY
! flag whether seismograms are ouput for North-East-Z component or Radial-Transverse-Z
  logical ROTATE_SEISMOGRAMS_RT

! flag to decide if seismograms are written by master proc only or
! by all processes in parallel (doing the later may create problems on some
! file systems)
  logical WRITE_SEISMOGRAMS_BY_MASTER

! save all seismograms in one large combined file instead of one file per seismogram
! to avoid overloading shared non-local file systems such as GPFS for instance
  logical SAVE_ALL_SEISMOS_IN_ONE_FILE
  logical USE_BINARY_FOR_LARGE_FILE

  allocate(one_seismogram(NDIM,NTSTEP_BETWEEN_OUTPUT_SEISMOS),stat=ier)
  if(ier /= 0) stop 'error while allocating one temporary seismogram'

  ! check that the sum of the number of receivers in each slice is nrec
  call MPI_REDUCE(nrec_local,nrec_tot_found,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ier)
  if(myrank == 0 .and. nrec_tot_found /= nrec) &
      call exit_MPI(myrank,'total number of receivers is incorrect')

  ! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

  ! all the processes write their local seismograms themselves
  if(.not. WRITE_SEISMOGRAMS_BY_MASTER) then

    write_time_begin = MPI_WTIME()

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

    total_seismos_local = 0

    ! loop on all the local receivers
    do irec_local = 1,nrec_local

      ! get global number of that receiver
      irec = number_receiver_global(irec_local)

      total_seismos_local = total_seismos_local + 1

      one_seismogram = seismograms(:,irec_local,:)

      ! write this seismogram
      call write_one_seismogram(one_seismogram,irec, &
                             station_name,network_name,stlat,stlon,stele,stbur,nrec, &
                             ANGULAR_WIDTH_XI_IN_DEGREES,NEX_XI,DT,hdur,it_end, &
                             yr,jda,ho,mi,sec,tshift_cmt,t_shift, &
                             elat,elon,depth,event_name,cmt_lat, &
                             cmt_lon,cmt_depth,cmt_hdur,OUTPUT_FILES, &
                             OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM, &
                             OUTPUT_SEISMOS_SAC_BINARY,ROTATE_SEISMOGRAMS_RT, &
                             NTSTEP_BETWEEN_OUTPUT_SEISMOS,seismo_offset,seismo_current, &
                             SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,myrank)

    enddo

    ! create one large file instead of one small file per station to avoid file system overload
    if(OUTPUT_SEISMOS_ASCII_TEXT .and. SAVE_ALL_SEISMOS_IN_ONE_FILE) close(IOUT)

    if(total_seismos_local/= nrec_local) call exit_MPI(myrank,'incorrect total number of receivers saved')

    write_time = MPI_WTIME() - write_time_begin

    if(myrank == 0) then
     write(IMAIN,*)
     write(IMAIN,*) 'Writing the seismograms in parallel took ',write_time,' seconds'
     write(IMAIN,*)
    endif

  ! now only the master process does the writing of seismograms and
  ! collects the data from all other processes
  else ! WRITE_SEISMOGRAMS_BY_MASTER

    write_time_begin = MPI_WTIME()

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

       total_seismos = 0

       ! loop on all the slices
       do iproc = 0,NPROCTOT-1

         ! receive except from proc 0, which is me and therefore I already have this value
         sender = iproc
         if(iproc /= 0) then
           call MPI_RECV(nrec_local_received,1,MPI_INTEGER,sender,itag,MPI_COMM_WORLD,msg_status,ier)
           if(nrec_local_received < 0) call exit_MPI(myrank,'error while receiving local number of receivers')
         else
           nrec_local_received = nrec_local
         endif
         if (nrec_local_received > 0) then
           do irec_local = 1,nrec_local_received
             ! receive except from proc 0, which is myself and therefore I already have these values
             if(iproc == 0) then
               ! get global number of that receiver
               irec = number_receiver_global(irec_local)
               one_seismogram(:,:) = seismograms(:,irec_local,:)
             else
               call MPI_RECV(irec,1,MPI_INTEGER,sender,itag,MPI_COMM_WORLD,msg_status,ier)
               if(irec < 1 .or. irec > nrec) call exit_MPI(myrank,'error while receiving global receiver number')
               call MPI_RECV(one_seismogram,NDIM*seismo_current,CUSTOM_MPI_TYPE,sender,itag,MPI_COMM_WORLD,msg_status,ier)
             endif

             total_seismos = total_seismos + 1
             ! write this seismogram
             call write_one_seismogram(one_seismogram,irec, &
                                       station_name,network_name,stlat,stlon,stele,stbur,nrec, &
                                       ANGULAR_WIDTH_XI_IN_DEGREES,NEX_XI,DT,hdur,it_end, &
                                       yr,jda,ho,mi,sec,tshift_cmt,t_shift, &
                                       elat,elon,depth,event_name,cmt_lat, &
                                       cmt_lon,cmt_depth,cmt_hdur,OUTPUT_FILES, &
                                       OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM, &
                                       OUTPUT_SEISMOS_SAC_BINARY,ROTATE_SEISMOGRAMS_RT, &
                                       NTSTEP_BETWEEN_OUTPUT_SEISMOS,seismo_offset,seismo_current,&
                                       SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,myrank)
           enddo
         endif
       enddo

       write(IMAIN,*)
       write(IMAIN,*) 'Total number of receivers saved is ',total_seismos,' out of ',nrec
       write(IMAIN,*)

       if(total_seismos /= nrec) call exit_MPI(myrank,'incorrect total number of receivers saved')

       ! create one large file instead of one small file per station to avoid file system overload
       if(SAVE_ALL_SEISMOS_IN_ONE_FILE) close(IOUT)

    else  ! on the nodes, send the seismograms to the master
       receiver = 0
       call MPI_SEND(nrec_local,1,MPI_INTEGER,receiver,itag,MPI_COMM_WORLD,ier)
       if (nrec_local > 0) then
         do irec_local = 1,nrec_local
           ! get global number of that receiver
           irec = number_receiver_global(irec_local)
           call MPI_SEND(irec,1,MPI_INTEGER,receiver,itag,MPI_COMM_WORLD,ier)
           one_seismogram(:,:) = seismograms(:,irec_local,:)
           call MPI_SEND(one_seismogram,NDIM*seismo_current,CUSTOM_MPI_TYPE,receiver,itag,MPI_COMM_WORLD,ier)
         enddo
       endif
    endif

    write_time  = MPI_WTIME() - write_time_begin

    if(myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Writing the seismograms by master proc alone took ',write_time,' seconds'
      write(IMAIN,*)
    endif

  endif ! WRITE_SEISMOGRAMS_BY_MASTER

  deallocate(one_seismogram)

  end subroutine write_seismograms

!=====================================================================

  subroutine write_one_seismogram(one_seismogram,irec, &
              station_name,network_name,stlat,stlon,stele,stbur,nrec, &
              ANGULAR_WIDTH_XI_IN_DEGREES,NEX_XI,DT,hdur,it_end, &
              yr,jda,ho,mi,sec,tshift_cmt,t_shift,&
              elat,elon,depth,event_name,cmt_lat,cmt_lon,cmt_depth,cmt_hdur, &
              OUTPUT_FILES, &
              OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM, &
              OUTPUT_SEISMOS_SAC_BINARY,ROTATE_SEISMOGRAMS_RT, &
              NTSTEP_BETWEEN_OUTPUT_SEISMOS,seismo_offset,seismo_current, &
              SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,myrank)

  implicit none

  include "constants.h"

  integer nrec,it_end,NEX_XI

  integer :: seismo_offset, seismo_current, NTSTEP_BETWEEN_OUTPUT_SEISMOS

  real(kind=CUSTOM_REAL), dimension(NDIM,NTSTEP_BETWEEN_OUTPUT_SEISMOS) :: one_seismogram

  real(kind=CUSTOM_REAL), dimension(5,NTSTEP_BETWEEN_OUTPUT_SEISMOS) :: seismogram_tmp

  integer myrank
  double precision hdur,DT,ANGULAR_WIDTH_XI_IN_DEGREES

  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name

  integer irec,length_station_name,length_network_name
  integer iorientation

  character(len=4) chn
  character(len=256) sisname,sisname_big_file
  character(len=150) OUTPUT_FILES

  ! section added for SAC
  double precision tshift_cmt,t_shift,elat,elon,depth
  double precision cmt_lat,cmt_lon,cmt_depth,cmt_hdur

  double precision, dimension(nrec) :: stlat,stlon,stele,stbur

  ! variables for SAC header fields
  integer yr,jda,ho,mi
  double precision sec
  character(len=20) event_name

  ! flags to determine seismogram type
  logical OUTPUT_SEISMOS_ASCII_TEXT, OUTPUT_SEISMOS_SAC_ALPHANUM, &
          OUTPUT_SEISMOS_SAC_BINARY
  ! flag whether seismograms are ouput for North-East-Z component or Radial-Transverse-Z
  logical ROTATE_SEISMOGRAMS_RT

  ! save all seismograms in one large combined file instead of one file per seismogram
  ! to avoid overloading shared non-local file systems such as GPFS for instance
  logical SAVE_ALL_SEISMOS_IN_ONE_FILE
  logical USE_BINARY_FOR_LARGE_FILE

! local parameters
  character(len=2) bic
  ! variables used for calculation of backazimuth and
  ! rotation of components if ROTATE_SEISMOGRAMS=.true.
  integer ior_start,ior_end
  double precision backaz
  real(kind=CUSTOM_REAL) phi,cphi,sphi
  integer isample

  !----------------------------------------------------------------

  call band_instrument_code(DT,bic)
  if (ROTATE_SEISMOGRAMS_RT) then ! iorientation 1=N,2=E,3=Z,4=R,5=T
    ior_start=3    ! starting from Z
    ior_end  =5    ! ending with T => ZRT
  else
    ior_start=1    ! starting from N
    ior_end  =3    ! ending with Z => NEZ
  endif

    !do iorientation = 1,NDIM
    !do iorientation = 1,5                   ! BS BS changed from 3 (NEZ) to 5 (NEZRT) components
  do iorientation = ior_start,ior_end      ! BS BS changed according to ROTATE_SEISMOGRAMS_RT

    if(iorientation == 1) then
      !chn = 'LHN'
      chn = bic(1:2)//'N'
    else if(iorientation == 2) then
      !chn = 'LHE'
      chn = bic(1:2)//'E'
    else if(iorientation == 3) then
      !chn = 'LHZ'
      chn = bic(1:2)//'Z'
    else if(iorientation == 4) then
      !chn = 'LHR'
      chn = bic(1:2)//'R'
    else if(iorientation == 5) then
      !chn = 'LHT'
      chn = bic(1:2)//'T'
    else
      call exit_MPI(myrank,'incorrect channel value')
    endif

    if (iorientation == 4 .or. iorientation == 5) then        ! LMU BS BS

      ! BS BS calculate backazimuth needed to rotate East and North
      ! components to Radial and Transverse components
      !  call get_backazimuth(elat,elon,stlat(irec),stlon(irec),backaz)
      call get_backazimuth(cmt_lat,cmt_lon,stlat(irec),stlon(irec),backaz)

      phi = backaz
      if (phi>180.) then
         phi = phi-180.
      elseif (phi<180.) then
         phi = phi+180.
      elseif (phi==180.) then
         phi = backaz
      endif

      cphi=cos(phi*pi/180)
      sphi=sin(phi*pi/180)

      ! BS BS do the rotation of the components and put result in
      ! new variable seismogram_tmp
      if (iorientation == 4) then ! radial component
         do isample = 1,seismo_current
            seismogram_tmp(iorientation,isample) = &
               cphi * one_seismogram(1,isample) + sphi * one_seismogram(2,isample)
         enddo
      elseif (iorientation == 5) then ! transverse component
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
    if (OUTPUT_SEISMOS_SAC_ALPHANUM .or. OUTPUT_SEISMOS_SAC_BINARY) then

      call write_output_SAC(seismogram_tmp,irec, &
              station_name,network_name,stlat,stlon,stele,stbur,nrec, &
              ANGULAR_WIDTH_XI_IN_DEGREES,NEX_XI,DT,hdur,it_end, &
              yr,jda,ho,mi,sec,tshift_cmt,t_shift,&
              elat,elon,depth,event_name,cmt_lat,cmt_lon,cmt_depth,cmt_hdur, &
              OUTPUT_FILES, &
              OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
              NTSTEP_BETWEEN_OUTPUT_SEISMOS,seismo_offset,seismo_current, &
              iorientation,phi,chn,sisname)

    endif ! OUTPUT_SEISMOS_SAC_ALPHANUM .or. OUTPUT_SEISMOS_SAC_BINARY

    ! ASCII output format
    if(OUTPUT_SEISMOS_ASCII_TEXT) then

      call write_output_ASCII(seismogram_tmp, &
              DT,hdur,OUTPUT_FILES, &
              NTSTEP_BETWEEN_OUTPUT_SEISMOS,seismo_offset,seismo_current, &
              SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,myrank, &
              iorientation,sisname,sisname_big_file)

    endif  ! OUTPUT_SEISMOS_ASCII_TEXT

  enddo ! do iorientation

  end subroutine write_one_seismogram

!=====================================================================

! write adjoint seismograms to text files

 subroutine write_adj_seismograms(seismograms,number_receiver_global, &
              nrec_local,it,nit_written,DT,NSTEP, &
              NTSTEP_BETWEEN_OUTPUT_SEISMOS,hdur,LOCAL_PATH)

 implicit none

 include "constants.h"

 integer nrec_local,NSTEP,NTSTEP_BETWEEN_OUTPUT_SEISMOS,it,nit_written
 integer, dimension(nrec_local) :: number_receiver_global
 real(kind=CUSTOM_REAL), dimension(9,nrec_local,NSTEP) :: seismograms
 double precision hdur,DT
 character(len=150) LOCAL_PATH

 integer irec,irec_local
 integer iorientation,isample

 character(len=4) chn
 character(len=150) clean_LOCAL_PATH,final_LOCAL_PATH
 character(len=256) sisname
 character(len=2) bic

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
   clean_LOCAL_PATH = adjustl(LOCAL_PATH)

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
         write(IOUT,*) sngl(dble(isample-1)*DT - hdur),' ',seismograms(iorientation,irec_local,isample-nit_written)
       else
         write(IOUT,*) dble(isample-1)*DT - hdur,' ',seismograms(iorientation,irec_local,isample-nit_written)
       endif
     enddo

     close(IOUT)

     enddo

 enddo

 end subroutine write_adj_seismograms

!=====================================================================

 subroutine band_instrument_code(DT,bic)
  ! This subroutine is to choose the appropriate band and instrument codes for channel names of seismograms
  ! based on the IRIS convention (first two letters of channel codes which were LH(Z/E/N) previously).
  ! For consistency with observed data, we now use the IRIS convention for band codes (first letter in channel codes)of
  ! SEM seismograms governed by their sampling rate.
  ! Instrument code (second letter in channel codes) is fixed to "X" which is assigned by IRIS for synthetic seismograms.
  ! See the manual for further explanations!
  ! Ebru, November 2010
  implicit none
  double precision DT
  character(len=2) bic

  if (DT .ge. 1.0d0)  bic = 'LX'
  if (DT .lt. 1.0d0 .and. DT .gt. 0.1d0) bic = 'MX'
  if (DT .le. 0.1d0 .and. DT .gt. 0.0125d0) bic = 'BX'
  if (DT .le. 0.0125d0 .and. DT .gt. 0.004d0) bic = 'HX'
  if (DT .le. 0.004d0 .and. DT .gt. 0.001d0) bic = 'CX'
  if (DT .le. 0.001d0) bic = 'FX'

 end subroutine band_instrument_code
