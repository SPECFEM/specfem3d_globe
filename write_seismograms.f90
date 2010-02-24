!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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
            nrec,nrec_local,DT,hdur,it_end, &
            yr,jda,ho,mi,sec,t_cmt,t_shift, &
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
 integer nrec,nrec_local,myrank,it_end,NPROCTOT !,NSOURCES
 character(len=256) sisname

 integer :: seismo_offset, seismo_current, NTSTEP_BETWEEN_OUTPUT_SEISMOS
 integer, dimension(nrec_local) :: number_receiver_global

 real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS) :: seismograms
 double precision hdur,DT

 character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
 character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name
 double precision t_cmt,t_shift,elat,elon,depth
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
                             DT,hdur,it_end, &
                             yr,jda,ho,mi,sec,t_cmt,t_shift, &
                             elat,elon,depth,event_name,cmt_lat, &
                             cmt_lon,cmt_depth,cmt_hdur,OUTPUT_FILES, &
                             OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM, &
                             OUTPUT_SEISMOS_SAC_BINARY,ROTATE_SEISMOGRAMS_RT, &
                             NTSTEP_BETWEEN_OUTPUT_SEISMOS,seismo_offset,seismo_current, &
                             SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,myrank)

   enddo

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
                                       DT,hdur,it_end, &
                                       yr,jda,ho,mi,sec,t_cmt,t_shift, &
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

!
!----
!

  subroutine write_one_seismogram(one_seismogram,irec, &
              station_name,network_name,stlat,stlon,stele,stbur,nrec, &
              DT,hdur,it_end, &
              yr,jda,ho,mi,sec,t_cmt,t_shift,&
              elat,elon,depth,event_name,cmt_lat,cmt_lon,cmt_depth,cmt_hdur, &
              OUTPUT_FILES, &
              OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM, &
              OUTPUT_SEISMOS_SAC_BINARY,ROTATE_SEISMOGRAMS_RT, &
              NTSTEP_BETWEEN_OUTPUT_SEISMOS,seismo_offset,seismo_current, &
              SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,myrank)

! SAC headers have new format
! by Ebru

  implicit none

  include "constants.h"

  integer nrec,it_end

  integer :: seismo_offset, seismo_current, NTSTEP_BETWEEN_OUTPUT_SEISMOS

  real(kind=CUSTOM_REAL), dimension(NDIM,NTSTEP_BETWEEN_OUTPUT_SEISMOS) :: one_seismogram

  real(kind=CUSTOM_REAL), dimension(5,NTSTEP_BETWEEN_OUTPUT_SEISMOS) :: seismogram_tmp

  integer myrank
  double precision hdur,DT

  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name

  integer irec,length_station_name,length_network_name
  integer iorientation,isample
  double precision value

  character(len=4) chn
  character(len=150) sisname,sisname_big_file
  character(len=150) OUTPUT_FILES

  ! section added for SAC
  !integer NSOURCES

  double precision t_cmt,t_shift,elat,elon,depth
  double precision cmt_lat,cmt_lon,cmt_depth,cmt_hdur
  double precision value1,value2, value3,value4,value5

  double precision, dimension(nrec) :: stlat,stlon,stele,stbur

  character(len=256) sisname_2

  ! variables for SAC header fields
  integer yr,jda,ho,mi,time_sec
  double precision sec
  !real mb
!  character(len=12) ename
  character(len=20) event_name
  logical, external :: is_leap_year

  real DELTA
  real DEPMIN
  real DEPMAX
  real SCALE_F
  real ODELTA
  real B,E,O,A
  real STLA,STLO,STEL,STDP
  real EVLA,EVLO,EVEL,EVDP
  real MAG,DIST,AZ,BAZ,GCARC
  real DEPMEN
  real USER0 !,USER1,USER2,USER3,USER4
  real CMPAZ,CMPINC

  integer NZYEAR,NZJDAY,NZHOUR,NZMIN,NZSEC
  integer NZMSEC,NVHDR,NORID,NEVID
  ! NUMBER of POINTS:
  integer NPTS
  integer IFTYPE,IMAGTYP
  integer IDEP
  integer IZTYPE
  integer IEVTYP
  integer IQUAL
  integer ISYNTH
  ! permission flags:
  integer LEVEN
  integer LPSPOL
  integer LOVROK
  integer LCALDA

  character(len=8) KSTNM
  character(len=16) KEVNM
  character(len=8) KCMPNM
  character(len=8) KNETWK
  character(len=8) KUSER0,KUSER1,KUSER2
  character(len=8), parameter :: str_undef='-12345  '

  real UNUSED   ! header fields unused by SAC
  real undef    ! undefined values
  real INTERNAL ! SAC internal variables, always leave undefined
  real BYSAC
  ! end SAC header variables

  ! flags to determine seismogram type
  logical OUTPUT_SEISMOS_ASCII_TEXT, OUTPUT_SEISMOS_SAC_ALPHANUM, &
          OUTPUT_SEISMOS_SAC_BINARY
  ! flag whether seismograms are ouput for North-East-Z component or Radial-Transverse-Z
  logical ROTATE_SEISMOGRAMS_RT

  ! save all seismograms in one large combined file instead of one file per seismogram
  ! to avoid overloading shared non-local file systems such as GPFS for instance
  logical SAVE_ALL_SEISMOS_IN_ONE_FILE
  logical USE_BINARY_FOR_LARGE_FILE

  ! variables used for calculation of backazimuth and
  ! rotation of components if ROTATE_SEISMOGRAMS=.true.

  integer ior_start,ior_end
  double precision backaz
  real(kind=CUSTOM_REAL) phi,cphi,sphi
  !----------------------------------------------------------------

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
      chn = 'LHN'
    else if(iorientation == 2) then
      chn = 'LHE'
    else if(iorientation == 3) then
      chn = 'LHZ'
    else if(iorientation == 4) then
      chn = 'LHR'
    else if(iorientation == 5) then
      chn = 'LHT'
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

    if (OUTPUT_SEISMOS_SAC_ALPHANUM .or. OUTPUT_SEISMOS_SAC_BINARY) then

!######################## SAC Alphanumeric Seismos ############################
!
! written by Markus Treml and Bernhard Schuberth, Dept. for Earth and Environ-
! mental Sciences, Ludwig-Maximilians-University Munich, Germany
!
! some words about SAC timing:
!==============================
!
!NPTS,DELTA,B,E:
! These define the timing of the seismogram. E is calculated by sac. So, say
! you have 100 NPTS, a DELTA of 0.5, and set B to 0, E should be 50.
! Likewise setting B to -50 gives an E of 0.  Cutting basically cuts out points
! between the two times you designate based on these values.
!KZTIME and KZDATE:
! Now things get funky.  KZTIME defines the exact time that the trace begins
! at. It has no affect on timing per se.  You'll really notice its effect if
! you read in two traces from different dates.

! Reference markers, (e.g. the o-marker) are not defined relative to this time,
! but rather to the begin time (B) of the seismo, so if you adjust B, you also
! need to adjust KZTIME to match.  l would suggest experimenting with this until
! you understand it. It is a little non-intuitive until you see it for yourself.
!
!-----------------------------------------------------------------------------
!
! This file is essentially the alphanumeric equivalent of the SAC binary data
! file. The header section is stored on the first 30 cards. This is followed
! by one or two data sections. The data is in 5G15.7 format.
!----------------------------------------------------------------------

      ! define certain default values

      ! unused or undefined values are set to '-12345.00'
      UNUSED   = -12345.00 ! header fields unused by SAC
      undef    = -12345.00 ! undefined values
      INTERNAL = -12345.00 ! SAC internal variables, always left undefined
      BYSAC    = -12345.00 ! values calculated by SAC from other variables
      !
      DELTA  = DT          ! [REQUIRED]
      DEPMIN = BYSAC
      DEPMAX = BYSAC
      DEPMEN = BYSAC
      SCALE_F= 1000000000  ! factor for y-value, set to 10e9, so that values are in nm
      ODELTA = undef       ! increment from delta

      B      = sngl((seismo_offset)*DT-hdur + t_cmt) ! [REQUIRED]
      E      = BYSAC       ! [REQUIRED]
      O      = 0  !
      A      = undef  !###
      !station values:
      STLA = stlat(irec)
      STLO = stlon(irec)
      STEL = stele(irec)
      STDP = stbur(irec)

      !event values (hypocenter):
      ! note: this writes out the CMT location, which might be different 
      ! to the event location given in the first, PDE line
      EVLA   = cmt_lat
      EVLO   = cmt_lon
      EVEL   = undef  !not defined
      EVDP   = cmt_depth


      ! by Ebru
      ! SAC headers will have new format
      USER0  = cmt_hdur !half duration from CMT file if not changed to hdur=0.d0 (point source)
      
      ! we remove any PDE information, since the simulation could also start
      ! with a "pure" CMT solution, without having any PDE infos
      !
      !USER1  = t_shift !time shift between PDE and CMT solutions 
      !PDE location values (different from CMT location, usually):
      !USER2  = depth !PDE depth
      !USER3  = elat !PDE event latitude
      !USER4  = elon !PDE event longitude
      !  
      !cmt location values (different from hypocenter location, usually):
      ! USER0  = cmt_lat
      ! USER1  = cmt_lon
      !USER0  = elat
      !USER1  = elon
      !USER2  = depth
      !USER3  = cmt_hdur !half duration from CMT if not changed to hdur=0.d0 (point source)

      ! just to avoid compiler warning
      value1 = elat
      value1 = elon
      value1 = depth
        

      ! it is not clear, which magnitude to write out: 
      ! should it be 
      !   body-wave-magnitude (Mb), surface-wave-magnitude (Ms), moment magnitude (Mw)
      !   or leave magnitude and use scalar moment (M0, but calculated by which convention, Harvard?)
      !
      ! it's confusing, and as a result, we will omit it.
      ! by Ebru      
      MAG    = undef
      IMAGTYP= undef

      !MAG    = mb    !
      !IMAGTYP= 52    ! 52 = Mb? 55 = Mw!

      DIST   = BYSAC ! cause
      AZ     = BYSAC ! LCALDA
      BAZ    = BYSAC ! is
      GCARC  = BYSAC ! TRUE

      ! instrument orientation
      if(iorientation == 1) then !N
        CMPAZ  = 0.00
        CMPINC =90.00
      else if(iorientation == 2) then !E
        CMPAZ  =90.00
        CMPINC =90.00
      else if(iorientation == 3) then !Z
        CMPAZ  = 0.00
        CMPINC = 0.00
      else if(iorientation == 4) then !R
        CMPAZ = modulo(phi,360.) ! phi is calculated above (see call distaz())
        CMPINC =90.00
      else if(iorientation == 5) then !T
        CMPAZ = modulo(phi+90.,360.) ! phi is calculated above (see call distaz())
        CMPINC =90.00
      endif
      !----------------end format G15.7--------

      ! date and time:
      NZYEAR =yr
      NZJDAY =jda
      NZHOUR =ho
      NZMIN  =mi

      ! adds time-shift to get the CMT time in the headers as origin time of events  
      ! by Ebru      
      NZSEC  =int(sec+t_shift)
      NZMSEC =int((sec+t_shift-int(sec+t_shift))*1000)

      !NZSEC  =int(sec)
      !NZMSEC =int((sec-int(sec))*1000)

      ! Adjust event time and date after t_shift is added 
      if (NZSEC >= 60) then
       time_sec = jda*24*3600 + ho*3600 + mi*60 + int(sec+t_shift)
       NZJDAY   = int(time_sec/(24*3600))
       NZHOUR   = int(mod(time_sec,24)/3600)
       NZMIN    = int(mod(time_sec,3600)/60)
       NZSEC    = mod(time_sec,60)
       if (NZJDAY  > 365 .and. .not. is_leap_year(NZYEAR)) then
          NZJDAY = mod(NZJDAY,365)
          NZYEAR = yr + 1
       elseif (NZJDAY  > 366 .and. is_leap_year(NZYEAR)) then
          NZJDAY = mod(NZJDAY,366)
          NZYEAR = yr + 1
       elseif (NZJDAY == 366 .and. is_leap_year(NZYEAR)) then
          NZJDAY = 366
       endif
      endif


      NVHDR=6 ! SAC header version number. Current is 6

      ! CSS3.0 variables:
      NORID =int(undef) !origin ID
      NEVID =int(undef) !event  ID
      !NWVID =undef !waveform ID

      ! NUMBER of POINTS:
      NPTS = it_end-seismo_offset ! [REQUIRED]
      ! event type
      IFTYPE = 1 ! 1=ITIME, i.e. seismogram  [REQUIRED] # numbering system is
      IDEP   = 6 ! 6: displ/nm                          # quite strange, best

      IZTYPE = 11 !=origint reference time equivalent ! # by chnhdr and write
      IEVTYP = 40 !event type, 40: Earthquake           # alpha and check
      IQUAL  = int(undef) ! quality
      ISYNTH = int(undef) ! 1 real data, 2...n synth. flag
      ! permission flags:
      LEVEN =1 ! evenly spaced data [REQUIRED]
      LPSPOL=1 ! ? pos. polarity of components (has to be TRUE for LCALDA=1)
      LOVROK=1 ! 1: OK to overwrite file on disk
      LCALDA=1 ! 1: calculate DIST, AZ, BAZ, and GCARC, 0: do nothing
      ! ------------------end format 5I10---------
      !
      !----------------------------------
      KSTNM  = station_name(irec)(1:8) ! A8

      ! writes out event id as event name
      ! by Ebru
      KEVNM  = event_name(1:len_trim(event_name)) ! A16
      
      !if (NSOURCES == 1) then
      !  KEVNM  = ename(1:len_trim(ename))//'_syn'! A16
      !else
      !  KEVNM  = ename(1:len_trim(ename))//'_sFS'! A16
      !endif

      KCMPNM = chn(1:3)           ! 3A8
      KNETWK = network_name(irec) !  A6

      ! indicates SEM synthetics
      ! by Ebru
      KUSER0 = 'SEM'          !  A8
      KUSER1 = 'v5.0.0 ' 
      KUSER2 = 'Tiger' ! aka. awesome (princeton) tiger version :)
      
      !KUSER0 = 'PDE_LAT_'          !  A8
      !KUSER1 = 'PDE_LON_'          !  A8
      !KUSER2 = 'PDEDEPTH'          !  A8
      !----------------------------------

      if (OUTPUT_SEISMOS_SAC_ALPHANUM) then

        ! add .sacan (sac alphanumeric) extension to seismogram file name for SAC seismograms
        write(sisname_2,"('/',a,'.sacan')") trim(sisname)
        if (seismo_offset == 0) then
          open(unit=IOUT_SAC,file=trim(OUTPUT_FILES)//trim(sisname_2),&
            status='unknown',action='write')
        else
          open(unit=IOUT_SAC,file=trim(OUTPUT_FILES)//trim(sisname_2),&
            status='old', position='append',action='write')
        endif

! Formats of alphanumerical SAC header fields
510 format(5G15.7,5G15.7,5G15.7,5G15.7,5G15.7)
520 format(5I10,5I10,5I10,5I10,5I10)
530 format(A8,A16)
540 format(A8,A8,A8)


        if (seismo_offset == 0) then
          !
          ! now write actual header:
          ! ------------------------
          !
          ! real variables:
          !                                 DELTA     DEPMIN   DEPMAX   SCALE   ODELTA
          !                                 B         E        O        A       INTERNAL
          !                                 T0        T1       T2       T3      T4
          !                                 T5        T6       T7       T8      T9
          !                                 F         RESP0    RESP1    RESP2   RESP3
          !                                 RESP4     RESP5    RESP6    RESP7   RESP8
          !                                 RESP9     STLA     STLO     STEL    STDP
          !                                 EVLA      EVLO     EVEL     EVDP    MAG
          !                                 USER0     USER1    USER2    USER3   USER4
          !                                 USER5     USER6    USER7    USER8   USER9
          !                                 DIST      AZ       BAZ      GCARC   INTERNAL
          !                                 INTERNAL  DEPMEN   CMPAZ    CMPINC  XMINIMUM
          !                                 XMAXIMUM  YMINIMUM YMAXIMUM ADJTM   UNUSED
          !
          write(IOUT_SAC,510) DELTA,    DEPMIN,  DEPMAX,  SCALE_F,  ODELTA
          write(IOUT_SAC,510) B,        E,       O,       A,      INTERNAL
          write(IOUT_SAC,510) undef,    undef,   undef,   undef,  undef
          write(IOUT_SAC,510) undef,    undef,   undef,   undef,  undef
          write(IOUT_SAC,510) undef,    undef,   undef,   undef,  undef
          write(IOUT_SAC,510) undef,    undef,   undef,   undef,  undef
          write(IOUT_SAC,510) undef,    STLA,    STLO,    STEL,   STDP
          write(IOUT_SAC,510) EVLA,     EVLO,    EVEL,    EVDP,   MAG
          write(IOUT_SAC,510) USER0,    undef,   undef,   undef,  undef
          !write(IOUT_SAC,510) USER0,    USER1,   USER2,   USER3,  USER4
          write(IOUT_SAC,510) undef,    undef,   undef,   undef,  undef
          write(IOUT_SAC,510) DIST,     AZ,      BAZ,     GCARC,  INTERNAL
          write(IOUT_SAC,510) INTERNAL, DEPMEN,  CMPAZ,   CMPINC, undef
          write(IOUT_SAC,510) undef,    undef,   undef,   undef,  undef
          write(IOUT_SAC,510) UNUSED,   UNUSED,  UNUSED,  UNUSED, UNUSED
          !
          ! integer variables:
          !                                 NSPTS, NWFID, NXSIZE, NYSIZE, UNUSED
          !                                                                    IINST
          !                                 ISTREG IEVREG IEVTYP IQUAL ISYNTH
          !                                 IMAGTYP, IMAGSRC, UNUSED, UNUSED, UNUSED
          !
          write(IOUT_SAC,520) NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC
          write(IOUT_SAC,520) NZMSEC, NVHDR, NORID, NEVID, NPTS
          write(IOUT_SAC,520) int(undef),int(undef),int(undef),int(undef),int(undef)
          write(IOUT_SAC,520) IFTYPE, IDEP, IZTYPE, int(UNUSED), int(undef)
          write(IOUT_SAC,520) int(undef),int(undef),IEVTYP, int(undef), ISYNTH
          write(IOUT_SAC,520) IMAGTYP,int(undef),int(undef),int(undef),int(undef)
          write(IOUT_SAC,520) int(UNUSED), int(UNUSED), int(UNUSED), int(UNUSED), int(UNUSED)
          write(IOUT_SAC,520) LEVEN, LPSPOL, LOVROK, LCALDA, int(UNUSED)
          write(IOUT_SAC,530) KSTNM, KEVNM
          !
          ! character variables:
          !
          !                                   KHOLE    KO       KA
          !                                   KT0      KT1      KT2
          !                                   KT3      KT4      KT5
          !                                   KT6      KT7      KT8
          !                                   KT9      KF       KUSER0
          !                                   KUSER1     KUSER2       KCMPNM
          !                                   KNETWK   KDATRD   KINST
          !
          write(IOUT_SAC,540) '-12345  ','-12345  ','-12345  '
          write(IOUT_SAC,540) '-12345  ','-12345  ','-12345  '
          write(IOUT_SAC,540) '-12345  ','-12345  ','-12345  '
          write(IOUT_SAC,540) '-12345  ','-12345  ','-12345  '
          write(IOUT_SAC,540) '-12345  ','-12345  ',KUSER0
          write(IOUT_SAC,540)   KUSER1, KUSER2, KCMPNM
          write(IOUT_SAC,540)   KNETWK,'-12345  ','-12345  '
        endif

        ! now write data - with five values per row:
        ! ---------------

        do isample = 1+5,seismo_current+1,5

          value1 = dble(seismogram_tmp(iorientation,isample-5))
          value2 = dble(seismogram_tmp(iorientation,isample-4))
          value3 = dble(seismogram_tmp(iorientation,isample-3))
          value4 = dble(seismogram_tmp(iorientation,isample-2))
          value5 = dble(seismogram_tmp(iorientation,isample-1))

          write(IOUT_SAC,510) sngl(value1),sngl(value2),sngl(value3),sngl(value4),sngl(value5)

        enddo

        close(IOUT_SAC)

      endif ! OUTPUT_SEISMOS_SAC_ALPHANUM

      ! For explaination on values set, see above (SAC ASCII)
      if (OUTPUT_SEISMOS_SAC_BINARY) then

        ! add .sac (sac binary) extension to seismogram file name for SAC seismograms
        write(sisname_2,"('/',a,'.sac')") trim(sisname)

        ! open binary file
        if (seismo_offset == 0) then
          call open_file_create(trim(OUTPUT_FILES)//trim(sisname_2)//char(0))
        else
          call open_file_append(trim(OUTPUT_FILES)//trim(sisname_2)//char(0))
        endif

        if (seismo_offset == 0) then
          ! write header variables

          ! write single precision header variables 1:70
          call write_real(DELTA)         !(1)
          call write_real(DEPMIN)        !(2)
          call write_real(DEPMAX)        !(3)
          call write_real(SCALE_F)       !(4)
          call write_real(ODELTA)        !(5)
          call write_real(B)             !(6)
          call write_real(E)             !(7)
          call write_real(O)             !(8)
          call write_real(A)             !(9)
          call write_real(INTERNAL)      !(10)
          call write_real(undef)          !(11)T0
          call write_real(undef)          !(12)T1
          call write_real(undef)          !(13)T2
          call write_real(undef)          !(14)T3
          call write_real(undef)          !(15)T4
          call write_real(undef)          !(16)T5
          call write_real(undef)          !(17)T6
          call write_real(undef)          !(18)T7
          call write_real(undef)          !(19)T8
          call write_real(undef)          !(20)T9
          call write_real(undef)          !(21)F
          call write_real(undef)          !(22)RESP0
          call write_real(undef)          !(23)RESP1
          call write_real(undef)          !(24)RESP2
          call write_real(undef)          !(25)RESP3
          call write_real(undef)          !(26)RESP4
          call write_real(undef)          !(27)RESP5
          call write_real(undef)          !(28)RESP6
          call write_real(undef)          !(29)RESP7
          call write_real(undef)          !(30)RESP8
          call write_real(undef)          !(31)RESP9
          call write_real(STLA)          !(32)
          call write_real(STLO)          !(33)
          call write_real(STEL)          !(34)
          call write_real(STDP)          !(35)
          call write_real(EVLA)          !(36)
          call write_real(EVLO)          !(37)
          call write_real(EVEL)          !(38)
          call write_real(EVDP)          !(39)
          call write_real(MAG)           !(40)
          call write_real(USER0)         !(41)USER0
          call write_real(undef)         !(42)USER1
          call write_real(undef)         !(43)USER2
          call write_real(undef)         !(44)USER3
          call write_real(undef)          !(45)USER4
          call write_real(undef)          !(46)USER5
          call write_real(undef)          !(47)USER6
          call write_real(undef)          !(48)USER7
          call write_real(undef)          !(49)USER8
          call write_real(undef)          !(50)USER9
          call write_real(DIST)          !(51)
          call write_real(AZ)            !(52)
          call write_real(BAZ)           !(53)
          call write_real(GCARC)         !(54)
          call write_real(INTERNAL)      !(55)
          call write_real(INTERNAL)      !(56)
          call write_real(DEPMEN)        !(57)
          call write_real(CMPAZ)         !(58)
          call write_real(CMPINC)        !(59)
          call write_real(undef)          !(60)XMINIMUM
          call write_real(undef)          !(61)XMAXIMUM
          call write_real(undef)          !(62)YMINIMUM
          call write_real(undef)          !(63)YMAXIMUM
          call write_real(undef)          !(64)
          call write_real(undef)          !(65)
          call write_real(undef)          !(66)
          call write_real(undef)          !(67)
          call write_real(undef)          !(68)
          call write_real(undef)          !(69)
          call write_real(undef)          !(70)

          ! write integer header variables 71:105
          call write_integer(NZYEAR)        !(71)
          call write_integer(NZJDAY)        !(72)
          call write_integer(NZHOUR)        !(73)
          call write_integer(NZMIN)         !(74)
          call write_integer(NZSEC)         !(75)
          call write_integer(NZMSEC)        !(76)
          call write_integer(NVHDR)         !(77)
          call write_integer(NORID)         !(78)
          call write_integer(NEVID)         !(79)
          call write_integer(NPTS)          !(80)
          call write_integer(int(undef))     !(81)UNUSED
          call write_integer(int(undef))     !(82)NWFID
          call write_integer(int(undef))     !(83)NXSIZE
          call write_integer(int(undef))     !(84)NYSIZE
          call write_integer(int(undef))     !(85)UNUSED
          call write_integer(IFTYPE)        !(86)
          call write_integer(IDEP)          !(87)
          call write_integer(IZTYPE)        !(88)
          call write_integer(int(undef))     !(89)UNUSED
          call write_integer(int(undef))     !(90)IINST
          call write_integer(int(undef))     !(91)ISTREG
          call write_integer(int(undef))     !(92)IEVREG
          call write_integer(IEVTYP)        !(93)
          call write_integer(int(undef))     !(94)IQUAL
          call write_integer(ISYNTH)        !(95)
          call write_integer(IMAGTYP)       !(96)
          call write_integer(int(undef))     !(97)IMAGSRC
          call write_integer(int(UNUSED))   !(98)
          call write_integer(int(UNUSED))   !(99)
          call write_integer(int(UNUSED))   !(100)
          call write_integer(int(UNUSED))   !(101)
          call write_integer(int(UNUSED))   !(102)
          call write_integer(int(UNUSED))   !(103)
          call write_integer(int(UNUSED))   !(104)
          call write_integer(int(UNUSED))   !(105)

          ! write logical header variables 106:110
          call write_integer(LEVEN)         !(106)
          call write_integer(LPSPOL)        !(107)
          call write_integer(LOVROK)        !(108)
          call write_integer(LCALDA)        !(109)
          call write_integer(int(UNUSED))   !(110)


          ! write character header variables 111:302
          call write_character(KSTNM,8)         !(111:118)
          call write_character(KEVNM,16)         !(119:134)
          call write_character(str_undef,8)      !(135:142)KHOLE
          call write_character(str_undef,8)      !(143:150)KO
          call write_character(str_undef,8)      !(151:158)KA
          call write_character(str_undef,8)      !(159:166)KT0
          call write_character(str_undef,8)      !(167:174)KT1
          call write_character(str_undef,8)      !(175:182)KT2
          call write_character(str_undef,8)      !(183:190)KT3
          call write_character(str_undef,8)      !(191:198)KT4
          call write_character(str_undef,8)      !(199:206)KT5
          call write_character(str_undef,8)      !(207:214)KT6
          call write_character(str_undef,8)      !(215:222)KT7
          call write_character(str_undef,8)      !(223:230)KT8
          call write_character(str_undef,8)      !(231:238)KT9
          call write_character(str_undef,8)      !(239:246)KF
          call write_character(KUSER0,8)        !(247:254)
          call write_character(KUSER1,8)        !(255:262)
          call write_character(KUSER2,8)        !(263:270)
          call write_character(KCMPNM,8)        !(271:278)
          call write_character(KNETWK,8)        !(279:286)
          call write_character(str_undef,8)      !(287:294)KDATRD
          call write_character(str_undef,8)      !(295:302)KINST

        endif

        ! now write SAC time series to file
        ! BS BS write whole time series at once (hope to increase I/O performance
        ! compared to using a loop on it)

        if (CUSTOM_REAL == SIZE_REAL) then
          call write_n_real(seismogram_tmp(iorientation,1:seismo_current),seismo_current)
        elseif (CUSTOM_REAL == SIZE_DOUBLE) then
          call write_n_real(real(seismogram_tmp(iorientation,1:seismo_current)),seismo_current)
        endif

        call close_file()


      endif ! OUTPUT_SEISMOS_SAC_BINARY

!#################### end SAC Alphanumeric Seismos ############################

    endif ! OUTPUT_SEISMOS_SAC_ALPHANUM .or. OUTPUT_SEISMOS_SAC_BINARY

    if(OUTPUT_SEISMOS_ASCII_TEXT) then

! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function

      ! add .ascii extension to seismogram file name for ASCII seismograms
      write(sisname_2,"('/',a,'.ascii')") trim(sisname)

      ! create one large file instead of one small file per station to avoid file system overload
      if(SAVE_ALL_SEISMOS_IN_ONE_FILE) then
        if(USE_BINARY_FOR_LARGE_FILE) then
          write(IOUT) sisname_big_file
        else
          write(IOUT,*) sisname_big_file(1:len_trim(sisname_big_file))
        endif
      else
        if (seismo_offset==0) then
          open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname_2),status='unknown',action='write')
        else
          open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname_2),status='old',position='append',action='write')
        endif
      endif

      ! subtract half duration of the source to make sure travel time is correct
      do isample = 1,seismo_current
        value = dble(seismogram_tmp(iorientation,isample))

        if(SAVE_ALL_SEISMOS_IN_ONE_FILE .and. USE_BINARY_FOR_LARGE_FILE) then
          ! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            write(IOUT) sngl(dble(seismo_offset+isample-1)*DT - hdur),sngl(value)
          else
            write(IOUT) dble(seismo_offset+isample-1)*DT - hdur,value
          endif
        else
          ! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            write(IOUT,*) sngl(dble(seismo_offset+isample-1)*DT - hdur),' ',sngl(value)
          else
            write(IOUT,*) dble(seismo_offset+isample-1)*DT - hdur,' ',value
          endif
        endif
      enddo

      if(.not. SAVE_ALL_SEISMOS_IN_ONE_FILE) close(IOUT)

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
 character(len=150) sisname,clean_LOCAL_PATH,final_LOCAL_PATH

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
       chn = 'LHN'
     else if(iorientation == 8) then
       chn = 'LHE'
     else if(iorientation == 9) then
       chn = 'LHZ'
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
