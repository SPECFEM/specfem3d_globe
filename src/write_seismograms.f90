!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
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
            network_name,stlat,stlon,stele,nrec,nrec_local,DT,hdur,it_end, &
            yr,jda,ho,mi,sec,t_cmt, &
            elat,elon,depth,mb,ename,cmt_lat,cmt_lon, &
            cmt_depth,cmt_hdur,NSOURCES,NPROCTOT, &
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
 integer nrec,nrec_local,myrank,it_end,NPROCTOT,NSOURCES
 character(len=256) sisname

 integer :: seismo_offset, seismo_current, NTSTEP_BETWEEN_OUTPUT_SEISMOS
 integer, dimension(nrec_local) :: number_receiver_global

 real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS) :: seismograms
 double precision hdur,DT

 character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
 character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name
 double precision t_cmt,elat,elon,depth
 double precision cmt_lat,cmt_lon,cmt_depth,cmt_hdur
 double precision, dimension(nrec) :: stlat,stlon,stele
 integer yr,jda,ho,mi
 double precision sec
 real mb
 character(12) ename

! variables
 integer :: iproc,sender,irec_local,irec,ier,receiver,nrec_local_received,nrec_tot_found
 integer :: total_seismos,total_seismos_local
 double precision :: write_time_begin,write_time

! allocate this automatic array in the memory stack to avoid memory fragmentation with "allocate()"
 real(kind=CUSTOM_REAL), dimension(NDIM,NTSTEP_BETWEEN_OUTPUT_SEISMOS) :: one_seismogram

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
          open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.bin', &
                   status='old',form='unformatted',position='append',action='write')
        endif
      else
        if (seismo_offset==0) then
          open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.ascii',status='unknown',form='formatted',action='write')
        else
          open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.ascii', &
                   status='old',form='formatted',position='append',action='write')
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
                             station_name,network_name,stlat,stlon,stele,nrec, &
                             DT,hdur,it_end, &
                             yr,jda,ho,mi,sec,t_cmt,elat,elon,depth,mb,ename,cmt_lat, &
                             cmt_lon,cmt_depth,cmt_hdur,NSOURCES,OUTPUT_FILES, &
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
             open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.bin', &
                   status='old',form='unformatted',position='append',action='write')
           endif
         else
           if (seismo_offset==0) then
             open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.ascii',status='unknown',form='formatted',action='write')
           else
             open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname)//'.ascii', &
                   status='old',form='formatted',position='append',action='write')
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
                                       station_name,network_name,stlat,stlon,stele,nrec, &
                                       DT,hdur,it_end, &
                                       yr,jda,ho,mi,sec,t_cmt,elat,elon,depth,mb,ename,cmt_lat, &
                                       cmt_lon,cmt_depth,cmt_hdur,NSOURCES,OUTPUT_FILES, &
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

  end subroutine write_seismograms

!
!----
!

  subroutine write_one_seismogram(one_seismogram,irec, &
              station_name,network_name,stlat,stlon,stele,nrec, &
              DT,hdur,it_end, &
              yr,jda,ho,mi,sec,t_cmt,elat,elon,depth,mb,ename,cmt_lat,cmt_lon,cmt_depth,cmt_hdur,NSOURCES, &
              OUTPUT_FILES, &
              OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM, &
              OUTPUT_SEISMOS_SAC_BINARY,ROTATE_SEISMOGRAMS_RT, &
              NTSTEP_BETWEEN_OUTPUT_SEISMOS,seismo_offset,seismo_current, &
              SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,myrank)

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
  integer NSOURCES

  double precision t_cmt,elat,elon,depth
  double precision cmt_lat,cmt_lon,cmt_depth,cmt_hdur

  double precision, dimension(nrec) :: stlat,stlon,stele

  character(len=256) sisname_2

! variables for SAC header fields
  integer yr,jda,ho,mi
  double precision sec
  real mb
  character(12) ename

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
  real USER0,USER1,USER2,USER3
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

  character(8) KSTNM
  character(16) KEVNM
  character(8) KCMPNM
  character(8) KNETWK
  character(8) KUSER0,KUSER1,KUSER2
  character(8), parameter :: str_undef='-12345  '

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

          if (backaz>180.) then
             phi=backaz-180.
          elseif (backaz<180.) then
             phi=backaz+180.
          elseif (backaz==180.) then
             phi=backaz
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
     write(sisname,"('/',a,'.',a,'.',a3,'.semd')") station_name(irec)(1:length_station_name), &
                   network_name(irec)(1:length_network_name),chn

! create this name also for the text line added to the unique big seismogram file
     write(sisname_big_file,"(a,'.',a,'.',a3,'.semd')") station_name(irec)(1:length_station_name), &
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
 O      = undef  !###
 A      = undef  !###
!station values:
 STLA = stlat(irec)
 STLO = stlon(irec)
 STEL = stele(irec)
 STDP = undef    !stdep(irec)
!event values (hypocenter):
 EVLA   = elat
 EVLO   = elon
 EVEL   = undef  !not defined
 EVDP   = depth

!cmt location values (different from hypocenter location, usually):
 USER0  = cmt_lat
 USER1  = cmt_lon
 USER2  = cmt_depth

 USER3  = cmt_hdur !half duration from CMT if not changed to hdur=0.d0 (point source)

 MAG    = mb    !
 IMAGTYP= 52    ! 52 = Mb? 55 = Mw!

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
 NZSEC  =int(sec)
 NZMSEC =int((sec-int(sec))*1000)

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

 if (NSOURCES == 1) then
   KEVNM  = ename(1:len_trim(ename))//'_syn'! A16
 else
   KEVNM  = ename(1:len_trim(ename))//'_sFS'! A16
 endif

!----------------------------------
 KCMPNM = chn(3:3)           ! 3A8
 KNETWK = network_name(irec) !  A6

 KUSER0 = 'CMT_LAT_'          !  A8
 KUSER1 = 'CMT_LON_'          !  A8
 KUSER2 = 'CMTDEPTH'          !  A8
!----------------------------------

  if (OUTPUT_SEISMOS_SAC_ALPHANUM) then

  endif ! OUTPUT_SEISMOS_SAC_ALPHANUM

! For explaination on values set, see above (SAC ASCII)
  if (OUTPUT_SEISMOS_SAC_BINARY) then

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

