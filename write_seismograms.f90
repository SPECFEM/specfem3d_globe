!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                    and University of Pau, France
! (c) California Institute of Technology and University of Pau, April 2007
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

! write seismograms to text files
subroutine write_seismograms(myrank,seismograms,number_receiver_global, &
              station_name,network_name,stlat,stlon,stele,nrec,nrec_local, &
              DT,NSTEP,hdur,it_begin,it_end, &
 yr,jda,ho,mi,sec,t_cmt,elat,elon,depth,mb,ename,cmt_lat,cmt_lon,cmt_depth,cmt_hdur,NSOURCES,&
 NPROCTOT,FINAL)

 implicit none

! standard include of the MPI library
 include 'mpif.h'

 include "constants.h"
 include "precision.h"

! parameters
 integer nrec,nrec_local,NSTEP,myrank,it_begin,it_end,NPROCTOT,NSOURCES
 integer, dimension(nrec_local) :: number_receiver_global
 real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NSTEP) :: seismograms
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
 logical FINAL

! variables
 integer :: iproc,sender,irec_local,irec,total_seismos,ier,receiver
 real(kind=CUSTOM_REAL), dimension(NDIM,NSTEP) :: one_seismogram
 integer msg_status(MPI_STATUS_SIZE)
 character(len=150) OUTPUT_FILES

  if(myrank == 0) then ! on the master, gather all the seismograms
    ! get the base pathname for output files
    call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

    total_seismos = 0
    ! receive information from all the slices
    do iproc = 0,NPROCTOT-1
      ! receive except from proc 0, which is me and therefore I already have this value
      sender = iproc
      if(iproc /= 0) call MPI_RECV(nrec_local,1,MPI_INTEGER,sender,itag,MPI_COMM_WORLD,msg_status,ier)
      if (nrec_local > 0) then
        do irec_local = 1,nrec_local
          ! receive except from proc 0, which is myself and therefore I already have these values
          if(iproc == 0) then
            ! get global number of that receiver
            irec = number_receiver_global(irec_local)
            one_seismogram(:,:) = seismograms(:,irec_local,:)
          else
            call MPI_RECV(irec,1,MPI_INTEGER,sender,itag,MPI_COMM_WORLD,msg_status,ier)
            call MPI_RECV(one_seismogram,NDIM*NSTEP,CUSTOM_MPI_TYPE,sender,itag,MPI_COMM_WORLD,msg_status,ier)
          endif
          total_seismos = total_seismos + 1
          call write_one_seismogram(one_seismogram,irec, &
                                    station_name,network_name,stlat,stlon,stele,nrec, &
                                    DT,NSTEP,hdur,it_begin,it_end, &
                                    yr,jda,ho,mi,sec,t_cmt,elat,elon,depth,mb,ename,cmt_lat, &
                                    cmt_lon,cmt_depth,cmt_hdur,NSOURCES,OUTPUT_FILES)
        enddo
      endif
    enddo
    if (FINAL) then
      write(IMAIN,*)
      write(IMAIN,*) 'Total number of receivers saved is ',total_seismos,' out of ',nrec
      write(IMAIN,*)
      if(total_seismos /= nrec) call exit_MPI(myrank, 'incorrect total number of receivers saved')
    endif
  else  ! on the nodes, send the seismograms to the master
    receiver = 0
    call MPI_SEND(nrec_local,1,MPI_INTEGER,receiver,itag,MPI_COMM_WORLD,ier)
    if (nrec_local > 0) then
      do irec_local = 1,nrec_local
        ! get global number of that receiver
        irec = number_receiver_global(irec_local)
        call MPI_SEND(irec,1,MPI_INTEGER,receiver,itag,MPI_COMM_WORLD,ier)
        one_seismogram(:,:) = seismograms(:,irec_local,:)
        call MPI_SEND(one_seismogram,NDIM*NSTEP,CUSTOM_MPI_TYPE,receiver,itag,MPI_COMM_WORLD,ier)
      enddo
    endif
  endif

end subroutine write_seismograms


 subroutine write_one_seismogram(one_seismogram,irec, &
              station_name,network_name,stlat,stlon,stele,nrec, &
              DT,NSTEP,hdur,it_begin,it_end, &
 yr,jda,ho,mi,sec,t_cmt,elat,elon,depth,mb,ename,cmt_lat,cmt_lon,cmt_depth,cmt_hdur,NSOURCES, &
 OUTPUT_FILES)

 implicit none

  include "constants.h"

  integer nrec,NSTEP,it_begin,it_end
  real(kind=CUSTOM_REAL), dimension(NDIM,NSTEP) :: one_seismogram
  double precision hdur,DT

  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name

  integer irec,length_station_name,length_network_name
  integer iorientation,isample
  double precision value

  character(len=4) chn
  character(len=150) sisname
  character(len=150) OUTPUT_FILES

! BS BS begin section added for SAC
  integer NSOURCES
  integer write_counter

  double precision t_cmt,elat,elon,depth
  double precision cmt_lat,cmt_lon,cmt_depth,cmt_hdur
  double precision value1,value2, value3,value4,value5

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

  real UNUSED   ! header fields unused by SAC
  real undef    ! undefined values
  real INTERNAL ! SAC internal variables, always leave undefined
  real BYSAC
! end BS BS SAC header variables

!----------------------------------------------------------------

   do iorientation = 1,NDIM

     if(iorientation == 1) then
       chn = 'LHN'
     else if(iorientation == 2) then
       chn = 'LHE'
     else if(iorientation == 3) then
       chn = 'LHZ'
     else
       stop 'incorrect channel value'
     endif

! create the name of the seismogram file for each slice
! file name includes the name of the station and the network
     length_station_name = len_trim(station_name(irec))
     length_network_name = len_trim(network_name(irec))

! check that length conforms to standard
     if(length_station_name < 1 .or. length_station_name > MAX_LENGTH_STATION_NAME) &
          stop 'wrong length of station name'

     if(length_network_name < 1 .or. length_network_name > MAX_LENGTH_NETWORK_NAME) &
          stop 'wrong length of network name'

! create the name of the seismogram file using the station name and network name
     write(sisname,"('/',a,'.',a,'.',a3,'.semd')") station_name(irec)(1:length_station_name), &
                   network_name(irec)(1:length_network_name),chn

! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function
   open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname),status='unknown')

! subtract half duration of the source to make sure travel time is correct
     do isample = it_begin,it_end
       value = dble(one_seismogram(iorientation,isample))
! distinguish between single and double precision for reals
       if(CUSTOM_REAL == SIZE_REAL) then
         write(IOUT,*) sngl(dble(isample-1)*DT - hdur),' ',sngl(value)
       else
         write(IOUT,*) dble(isample-1)*DT - hdur,' ',value
       endif
     enddo
     close(IOUT)

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

! add .sac extension to seismogram file name for SAC seismograms
 write(sisname_2,"('/',a,'.sac')") trim(sisname)
 open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname_2),status='unknown')

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
 B      = sngl((it_begin -1)*DT-hdur + t_cmt) ! [REQUIRED]
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
 NPTS = it_end-it_begin + 1 ! [REQUIRED]
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
 KSTNM  = station_name(irec) ! A8

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

! Formats of alphanumerical SAC header fields
510 format(5G15.7,5G15.7,5G15.7,5G15.7,5G15.7)
520 format(5I10,5I10,5I10,5I10,5I10)
530 format(A8,A16)
540 format(A8,A8,A8)
!
! now write actual header:
! ------------------------
!
! real variables:
!
 write(IOUT,510) DELTA,    DEPMIN,  DEPMAX,  SCALE_F,  ODELTA
!                                 DELTA     DEPMIN   DEPMAX   SCALE   ODELTA
 write(IOUT,510) B,        E,       O,       A,      INTERNAL
!                                 B         E        O        A       INTERNAL
 write(IOUT,510) undef,    undef,   undef,   undef,  undef
!                                 T0        T1       T2       T3      T4
 write(IOUT,510) undef,    undef,   undef,   undef,  undef
!                                 T5        T6       T7       T8      T9
 write(IOUT,510) undef,    undef,   undef,   undef,  undef
!                                 F         RESP0    RESP1    RESP2   RESP3
 write(IOUT,510) undef,    undef,   undef,   undef,  undef
!                                 RESP4     RESP5    RESP6    RESP7   RESP8
 write(IOUT,510) undef,    STLA,    STLO,    STEL,   STDP
!                                 RESP9     STLA     STLO     STEL    STDP
 write(IOUT,510) EVLA,     EVLO,    EVEL,    EVDP,   MAG
!                                 EVLA      EVLO     EVEL     EVDP    MAG
 write(IOUT,510) USER0,    USER1,   USER2,   USER3,  undef
!                                 USER0     USER1    USER2    USER3   USER4
 write(IOUT,510) undef,    undef,   undef,   undef,  undef
!                                 USER5     USER6    USER7    USER8   USER9
 write(IOUT,510) DIST,     AZ,      BAZ,     GCARC,  INTERNAL
!                                 DIST      AZ       BAZ      GCARC   INTERNAL
 write(IOUT,510) INTERNAL, DEPMEN,  CMPAZ,   CMPINC, undef
!                                 INTERNAL  DEPMEN   CMPAZ    CMPINC  XMINIMUM
 write(IOUT,510) undef,    undef,   undef,   undef,  undef
!                                 XMAXIMUM  YMINIMUM YMAXIMUM ADJTM   UNUSED
 write(IOUT,510) UNUSED,   UNUSED,  UNUSED,  UNUSED, UNUSED
!
! integer variables:
!
 write(IOUT,520) NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC
 write(IOUT,520) NZMSEC, NVHDR, NORID, NEVID, NPTS
 write(IOUT,520) int(undef),int(undef),int(undef),int(undef),int(undef)
!                                 NSPTS, NWFID, NXSIZE, NYSIZE, UNUSED
 write(IOUT,520) IFTYPE, IDEP, IZTYPE, int(UNUSED), int(undef)
!                                                                    IINST
 write(IOUT,520) int(undef),int(undef),IEVTYP, int(undef), ISYNTH
!                                 ISTREG IEVREG IEVTYP IQUAL ISYNTH
 write(IOUT,520) IMAGTYP,int(undef),int(undef),int(undef),int(undef)
!                                 IMAGTYP, IMAGSRC, UNUSED, UNUSED, UNUSED
 write(IOUT,520) int(UNUSED), int(UNUSED), int(UNUSED), int(UNUSED), int(UNUSED)
 write(IOUT,520) LEVEN, LPSPOL, LOVROK, LCALDA, int(UNUSED)
 write(IOUT,530) KSTNM, KEVNM
!
! character variables:
!
 write(IOUT,540) '-12345  ','-12345  ','-12345  '
!                                   KHOLE    KO       KA
 write(IOUT,540) '-12345  ','-12345  ','-12345  '
!                                   KT0      KT1      KT2
 write(IOUT,540) '-12345  ','-12345  ','-12345  '
!                                   KT3      KT4      KT5
 write(IOUT,540) '-12345  ','-12345  ','-12345  '
!                                   KT6      KT7      KT8
 write(IOUT,540) '-12345  ','-12345  ',KUSER0
!                                   KT9      KF       KUSER0
 write(IOUT,540)   KUSER1, KUSER2, KCMPNM
!                                   KUSER1     KUSER2       KCMPNM
 write(IOUT,540)   KNETWK,'-12345  ','-12345  '
!                                   KNETWK   KDATRD   KINST

! now write data - with five values per row:
! ---------------
 write_counter = 0

 do isample = it_begin+5,it_end+1,5

   value1 = dble(one_seismogram(iorientation,isample-5))
   value2 = dble(one_seismogram(iorientation,isample-4))
   value3 = dble(one_seismogram(iorientation,isample-3))
   value4 = dble(one_seismogram(iorientation,isample-2))
   value5 = dble(one_seismogram(iorientation,isample-1))

   write(IOUT,510) sngl(value1),sngl(value2),sngl(value3),sngl(value4),sngl(value5)

   write_counter=write_counter+1

 enddo

 close(IOUT)

!#################### end SAC Alphanumeric Seismos ############################

     enddo ! do iorientation = 1,3

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
           status='unknown')
   else if(it > NTSTEP_BETWEEN_OUTPUT_SEISMOS) then
      !append to existing file
      open(unit=IOUT,file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)),&
           status='old',position='append')
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
