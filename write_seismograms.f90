!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 5
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology July 2004
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
               DT,NSTEP,hdur,LOCAL_PATH,it_begin,it_end)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"

  integer nrec,nrec_local,NSTEP,myrank,it_begin,it_end
  integer, dimension(nrec_local) :: number_receiver_global
  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NSTEP) :: seismograms
  double precision hdur,DT
  character(len=150) LOCAL_PATH

  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name

  integer irec,irec_local,length_station_name,length_network_name
  integer iorientation,isample
  double precision value

  character(len=4) chn
  character(len=150) sisname

! BS BS begin section added for SAC
  double precision, dimension(nrec) :: stlat,stlon,stele
  double precision t_cmt,elat,elon,depth

  double precision cmt_lat,cmt_lon,cmt_depth,cmt_hdur

  double precision value1,value2, value3,value4,value5

  integer write_counter,i,ier

  integer NSOURCES

  integer, parameter :: IOUT_SAC=44
  character(len=256) sisname_2

! variables for SAC header fields
  integer yr,jda, ho, mi
  double precision sec
  real mb
  integer, parameter :: LENGTH_REGION_NAME = 150
  character(len=LENGTH_REGION_NAME) region

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
  real ndef     ! not defined values
  real INTERNAL ! SAC internal variables, always leave undefined
  real BYSAC
! end BS BS SAC header variables

!----------------------------------------------------------------

  do irec_local = 1,nrec_local

! get global number of that receiver
    irec = number_receiver_global(irec_local)

    do iorientation = 1,3

      if(iorientation == 1) then
        chn = 'LHN'
      else if(iorientation == 2) then
        chn = 'LHE'
      else if(iorientation == 3) then
        chn = 'LHZ'
      else
        call exit_MPI(myrank,'incorrect channel value')
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

      write(sisname,"(a,'.',a,'.',a3,'.semd')") station_name(irec)(1:length_station_name), &
                    network_name(irec)(1:length_network_name),chn

! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function
    open(unit=IOUT,file=trim(LOCAL_PATH)//'/'//trim(sisname),status='unknown')

! subtract half duration of the source to make sure travel time is correct
      do isample = it_begin,it_end
        value = dble(seismograms(iorientation,irec_local,isample))
! distinguish between single and double precision for reals
        if(CUSTOM_REAL == SIZE_REAL) then
          write(IOUT,*) sngl(dble(isample-1)*DT - hdur),' ',sngl(value)
        else
          write(IOUT,*) dble(isample-1)*DT - hdur,' ',value
        endif
      enddo

      close(IOUT)

! BS BS beginning of section added (SAC output)

! get event information for SAC header
  if(myrank == 0) call get_event_info(yr,jda,ho,mi,sec,t_cmt,elat,elon,depth,mb,region, &
                        cmt_lat,cmt_lon,cmt_depth,cmt_hdur,NSOURCES,LENGTH_REGION_NAME)
! broadcast the information read on the master to the nodes
  call MPI_BCAST(yr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(jda,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(ho,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(mi,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(NSOURCES,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  call MPI_BCAST(sec,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(t_cmt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(elat,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(elon,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(depth,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(cmt_lat,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(cmt_lon,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(cmt_depth,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(cmt_hdur,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

  call MPI_BCAST(region,LENGTH_REGION_NAME,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

  write(ename(1:12),'(a12)') region(1:12)

  do i=1,len_trim(ename)
    if (ename(i:i)==' ') ename(i:i)='_'
  enddo

  write(sisname_2,"(a,'.sac')") trim(sisname)
  open(unit=IOUT_SAC,file=trim(LOCAL_PATH)//'/'//trim(sisname_2),status='unknown')

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
!
!define certain default values
!
!  unused of undefined values are set to '-12345.00'
  UNUSED   = -12345.00 ! header fields unused by SAC
  ndef     = -12345.00 ! not defined values
  INTERNAL = -12345.00 ! SAC internal variables, always left undefined
  BYSAC    = -12345.00 ! values calculated by SAC from other variables
!
  DELTA  = DT          ! [REQUIRED]
  DEPMIN = BYSAC
  DEPMAX = BYSAC
  DEPMEN = BYSAC
  SCALE_F= 1000000000  ! factor for y-value, set to 10e9, so that values are in nm
  ODELTA = ndef        ! increment from delta
  B      = sngl((it_begin -1)*DT-hdur + t_cmt) ! [REQUIRED]
  E      = BYSAC       ! [REQUIRED]
  O      = ndef  !###
  A      = ndef  !###
!station values:
  STLA = stlat(irec)
  STLO = stlon(irec)
  STEL = stele(irec)
  STDP = ndef    !stdep(irec)
!event values (hypocenter):
  EVLA   = elat
  EVLO   = elon
  EVEL   = ndef  !not defined
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
  NORID =int(ndef) !origin ID
  NEVID =int(ndef) !event  ID
!NWVID =ndef !waveform ID

! NUMBER of POINTS:
  NPTS = it_end-it_begin + 1 ! [REQUIRED]
! event type
  IFTYPE = 1 ! 1=ITIME, i.e. seismogram  [REQUIRED] # numbering system is
  IDEP   = 6 ! 6: displ/nm                          # quite strange, best
  IZTYPE = 11 !=origint reference time equivalent ! # by chnhdr and write
  IEVTYP = 40 !event type, 40: Earthquake           # alpha and check
  IQUAL  = int(ndef) ! quality
  ISYNTH = int(ndef) ! 1 real data, 2...n synth. flag
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
  write(IOUT_SAC,510) DELTA,    DEPMIN,  DEPMAX,  SCALE_F,  ODELTA
!                                 DELTA     DEPMIN   DEPMAX   SCALE   ODELTA
  write(IOUT_SAC,510) B,        E,       O,       A,      INTERNAL
!                                 B         E        O        A       INTERNAL
  write(IOUT_SAC,510) ndef,     ndef,    ndef,    ndef,   ndef
!                                 T0        T1       T2       T3      T4
  write(IOUT_SAC,510) ndef,     ndef,    ndef,    ndef,   ndef
!                                 T5        T6       T7       T8      T9
  write(IOUT_SAC,510) ndef,     ndef,    ndef,    ndef,   ndef
!                                 F         RESP0    RESP1    RESP2   RESP3
  write(IOUT_SAC,510) ndef,     ndef,    ndef,    ndef,   ndef
!                                 RESP4     RESP5    RESP6    RESP7   RESP8
  write(IOUT_SAC,510) ndef,     STLA,    STLO,    STEL,   STDP
!                                 RESP9     STLA     STLO     STEL    STDP
  write(IOUT_SAC,510) EVLA,     EVLO,    EVEL,    EVDP,   MAG
!                                 EVLA      EVLO     EVEL     EVDP    MAG
  write(IOUT_SAC,510) USER0,    USER1,   USER2,   USER3,   ndef
!                                 USER0     USER1    USER2    USER3   USER4
  write(IOUT_SAC,510) ndef,     ndef,    ndef,    ndef,   ndef
!                                 USER5     USER6    USER7    USER8   USER9
  write(IOUT_SAC,510) DIST,     AZ,      BAZ,     GCARC,  INTERNAL
!                                 DIST      AZ       BAZ      GCARC   INTERNAL
  write(IOUT_SAC,510) INTERNAL, DEPMEN,  CMPAZ,   CMPINC, ndef
!                                 INTERNAL  DEPMEN   CMPAZ    CMPINC  XMINIMUM
  write(IOUT_SAC,510) ndef,     ndef,    ndef,    ndef,   ndef
!                                 XMAXIMUM  YMINIMUM YMAXIMUM ADJTM   UNUSED
  write(IOUT_SAC,510) UNUSED,   UNUSED,  UNUSED,  UNUSED, UNUSED
!
! integer variables:
!
  write(IOUT_SAC,520) NZYEAR, NZJDAY, NZHOUR, NZMIN, NZSEC
  write(IOUT_SAC,520) NZMSEC, NVHDR, NORID, NEVID, NPTS
  write(IOUT_SAC,520) int(ndef),int(ndef),int(ndef),int(ndef),int(ndef)
!                                 NSPTS, NWFID, NXSIZE, NYSIZE, UNUSED
  write(IOUT_SAC,520) IFTYPE, IDEP, IZTYPE, int(UNUSED), int(ndef)
!                                                                    IINST
  write(IOUT_SAC,520) int(ndef),int(ndef),IEVTYP, int(ndef), ISYNTH
!                                 ISTREG IEVREG IEVTYP IQUAL ISYNTH
  write(IOUT_SAC,520) IMAGTYP,int(ndef),int(ndef),int(ndef),int(ndef)
!                                 IMAGTYP, IMAGSRC, UNUSED, UNUSED, UNUSED
  write(IOUT_SAC,520) int(UNUSED), int(UNUSED), int(UNUSED), int(UNUSED), int(UNUSED)
  write(IOUT_SAC,520) LEVEN, LPSPOL, LOVROK, LCALDA, int(UNUSED)
  write(IOUT_SAC,530) KSTNM, KEVNM
!
! character variables:
!
  write(IOUT_SAC,540) '-12345  ','-12345  ','-12345  '
!                                   KHOLE    KO       KA
  write(IOUT_SAC,540) '-12345  ','-12345  ','-12345  '
!                                   KT0      KT1      KT2
  write(IOUT_SAC,540) '-12345  ','-12345  ','-12345  '
!                                   KT3      KT4      KT5
  write(IOUT_SAC,540) '-12345  ','-12345  ','-12345  '
!                                   KT6      KT7      KT8
  write(IOUT_SAC,540) '-12345  ','-12345  ',KUSER0
!                                   KT9      KF       KUSER0
  write(IOUT_SAC,540)   KUSER1, KUSER2, KCMPNM
!                                   KUSER1     KUSER2       KCMPNM
  write(IOUT_SAC,540)   KNETWK,'-12345  ','-12345  '
!                                   KNETWK   KDATRD   KINST

! now write data - with five values per row:
! ---------------
  write_counter = 0

  do isample = it_begin+5,it_end+1,5

    value1 = dble(seismograms(iorientation,irec_local,isample-5))
    value2 = dble(seismograms(iorientation,irec_local,isample-4))
    value3 = dble(seismograms(iorientation,irec_local,isample-3))
    value4 = dble(seismograms(iorientation,irec_local,isample-2))
    value5 = dble(seismograms(iorientation,irec_local,isample-1))

    write(IOUT_SAC,510) sngl(value1),sngl(value2),sngl(value3),sngl(value4),sngl(value5)

    write_counter=write_counter+1

  enddo

!#################### end SAC Alphanumeric Seismos ############################

  close(IOUT_SAC)

! BS BS end of section added (SAC output)

      enddo ! do iorientation = 1,3

  enddo ! do irec_local = 1,nrec_local

  end subroutine write_seismograms

!=====================================================================

! write adjoint seismograms to text files

  subroutine write_adj_seismograms(seismograms,number_receiver_global, &
               nrec_local,it,DT,NSTEP,hdur,LOCAL_PATH)

  implicit none

  include "constants.h"

  integer nrec_local,NSTEP,it
  integer, dimension(nrec_local) :: number_receiver_global
  real(kind=CUSTOM_REAL), dimension(6,nrec_local,NSTEP) :: seismograms
  double precision hdur,DT
  character(len=150) LOCAL_PATH

  integer irec,irec_local
  integer iorientation,isample

  character(len=4) chn
  character(len=150) sisname,clean_LOCAL_PATH,final_LOCAL_PATH

  do irec_local = 1,nrec_local

! get global number of that receiver
    irec = number_receiver_global(irec_local)

    do iorientation = 1,6

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
      endif

! create the name of the seismogram file for each slice
! file name includes the name of the station, the network and the component

      write(sisname,"(a,i5.5,'.',a,'.',a3,'.sem')") 'S',irec,'NT',chn

! suppress white spaces if any
    clean_LOCAL_PATH = adjustl(LOCAL_PATH)

! create full final local path
    final_LOCAL_PATH = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) // '/'

! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function
      open(unit=IOUT,file=final_LOCAL_PATH(1:len_trim(final_LOCAL_PATH))//sisname(1:len_trim(sisname)),status='unknown')
! make sure we never write more than the maximum number of time steps
! subtract half duration of the source to make sure travel time is correct
      do isample = 1,min(it,NSTEP)
! distinguish between single and double precision for reals
        if(CUSTOM_REAL == SIZE_REAL) then
          write(IOUT,*) sngl(dble(isample-1)*DT - hdur),' ',seismograms(iorientation,irec_local,isample)
        else
          write(IOUT,*) dble(isample-1)*DT - hdur,' ',seismograms(iorientation,irec_local,isample)
        endif
      enddo

      close(IOUT)

      enddo

  enddo

  end subroutine write_adj_seismograms

!=====================================================================

! BS BS beginning of additional routine for SAC output

  subroutine get_event_info(yr,jda,ho,mi,sec,t_cmt,elat,elon,depth,mb,region,&
                            cmt_lat,cmt_lon,cmt_depth,cmt_hdur,NSOURCES,LENGTH_REGION_NAME)

! written by Bernhard Schuberth

! This subroutine reads the first line of the DATA/CMTSOLUTION file
! and extracts event information needed for SAC or PITSA headers

  implicit none

  include "constants.h"

!--- arguments of the subroutine below

  integer, intent(out) :: NSOURCES,yr,jda,ho,mi

  real, intent(out) :: mb

  double precision, intent(out) :: sec,t_cmt,elat,elon,depth,cmt_lat,cmt_lon,cmt_depth,cmt_hdur

  integer, intent(in) :: LENGTH_REGION_NAME
  character(len=LENGTH_REGION_NAME), intent(out) :: region ! event name for SAC header

!--- local variables here

  integer ios,icounter,mo,da,julian_day

  real ms

  character(len=5) datasource
  character(len=150) string,dummystring,CMTSOLUTION

!
!---- read hypocenter info
!
  call get_value_string(CMTSOLUTION, 'solver.CMTSOLUTION','DATA/CMTSOLUTION')

  open(unit=821,file=CMTSOLUTION,iostat=ios,status='old',action='read')
  if(ios /= 0) stop 'error opening CMTSOLUTION file (in get_event_info)'

  icounter = 0
  do while(ios == 0)
    read(821,"(a)",iostat=ios) dummystring
    if(ios == 0) icounter = icounter + 1
  enddo
  close(821)
  if(mod(icounter,NLINES_PER_CMTSOLUTION_SOURCE) /= 0) &
    stop 'total number of lines in CMTSOLUTION file should be a multiple of NLINES_PER_CMTSOLUTION_SOURCE'
  NSOURCES = icounter / NLINES_PER_CMTSOLUTION_SOURCE
  if(NSOURCES < 1) stop 'need at least one source in CMTSOLUTION file'

  open(unit=821,file=CMTSOLUTION,status='old',action='read')

  ! example header line of CMTSOLUTION file
  !PDE 2003 09 25 19 50 08.93  41.78  144.08  18.0 7.9 8.0 Hokkaido, Japan
  !event_id, date,origin time,latitude,longitude,depth, mb, MS, region

  ! read header with event information
    read(821,*) datasource,yr,mo,da,ho,mi,sec,elat,elon,depth,mb,ms,region

    jda=julian_day(yr,mo,da)

  ! ignore line with event name
    read(821,"(a)") string

  ! read time shift
    read(821,"(a)") string
    read(string(12:len_trim(string)),*) t_cmt

  if (NSOURCES == 1) then

  ! read half duration
    read(821,"(a)") string
    read(string(15:len_trim(string)),*) cmt_hdur

  ! read latitude
    read(821,"(a)") string
    read(string(10:len_trim(string)),*) cmt_lat

  ! read longitude
    read(821,"(a)") string
    read(string(11:len_trim(string)),*) cmt_lon

  ! read depth
    read(821,"(a)") string
    read(string(7:len_trim(string)),*) cmt_depth

  else

    cmt_hdur=-1e8
    cmt_lat=-1e8
    cmt_lon=-1e8
    cmt_depth=-1e8

  endif

  close(821)

  end subroutine get_event_info

