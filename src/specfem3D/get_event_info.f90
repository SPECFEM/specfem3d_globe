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

! get information about event name and location for SAC seismograms: MPI version by Dimitri Komatitsch

! Instead of using region names as event names,
! event names given in the second row of CMT files will be used.
! Thus, I removed old parameters ename, region, LENGTH_REGION_NAME and added event_name!!!!!!!
! Also, t_shift is added as a new parameter to be written on sac headers!
! by Ebru Bozdag

  subroutine get_event_info_parallel(yr,jda,mo,da,ho,mi,sec, &
                                     event_name,tshift_src,t_shift, &
                                     elat,elon,depth,mb,ms,cmt_lat, &
                                     cmt_lon,cmt_depth,cmt_hdur,NSOURCES, &
                                     Mrr, Mtt, Mpp, Mrt, Mrp, Mtp)

  use constants, only: myrank

  implicit none

!--- input or output arguments of the subroutine below

  integer, intent(in) :: NSOURCES ! must be given

  integer, intent(out) :: yr,mo,da,jda,ho,mi
  double precision, intent(out) :: sec
  real, intent(out) :: mb,ms
  double precision, intent(out) :: tshift_src,elat,elon,depth,cmt_lat,cmt_lon,cmt_depth,cmt_hdur
  double precision, intent(out) :: Mrr, Mtt, Mpp, Mrt, Mrp, Mtp
  double precision, intent(out) :: t_shift

  character(len=20), intent(out) :: event_name

  ! get event information for SAC header on the master
  if (myrank == 0) then

    ! note: mb as (body wave) moment magnitude is not used any further,
    !       see comment in write_output_SAC() routine

    call get_event_info_serial(yr,jda,mo,da,ho,mi,sec,event_name,tshift_src,t_shift, &
                               elat,elon,depth,mb,ms, &
                               cmt_lat,cmt_lon,cmt_depth,cmt_hdur,NSOURCES, &
                               Mrr, Mtt, Mpp, Mrt, Mrp, Mtp)

    ! create the event name
    !write(ename(1:12),'(a12)') region(1:12)

    ! replace white spaces with underscores in event name
    !do i = 1,len_trim(ename)
    !  if (ename(i:i) == ' ') ename(i:i) = '_'
    !enddo

  endif

  ! broadcast the information read on the master to the nodes
  call bcast_all_singlei(yr)
  call bcast_all_singlei(jda)
  call bcast_all_singlei(mo)
  call bcast_all_singlei(da)
  call bcast_all_singlei(ho)
  call bcast_all_singlei(mi)

  call bcast_all_singledp(sec)

  call bcast_all_singledp(tshift_src)
  call bcast_all_singledp(t_shift)

  ! event location given on first, PDE line
  call bcast_all_singledp(elat)
  call bcast_all_singledp(elon)
  call bcast_all_singledp(depth)

  ! body and surfaace wave magnitude given on first line
  call bcast_all_singler(mb)
  call bcast_all_singler(ms)

  ! cmt location given in CMT file
  call bcast_all_singledp(cmt_lat)
  call bcast_all_singledp(cmt_lon)
  call bcast_all_singledp(cmt_depth)
  call bcast_all_singledp(cmt_hdur)

  ! moment tensor given in CMT file
  call bcast_all_singledp(Mrr)
  call bcast_all_singledp(Mtt)
  call bcast_all_singledp(Mpp)
  call bcast_all_singledp(Mrt)
  call bcast_all_singledp(Mrp)
  call bcast_all_singledp(Mtp)

  call bcast_all_ch(event_name,20)

  end subroutine get_event_info_parallel

!=====================================================================

! get information about event name and location for SAC seismograms: MPI version by Bernhard Schuberth
! This subroutine reads the first line of the DATA/CMTSOLUTION file
! and extracts event information needed for SAC or PITSA headers

! This subroutine has been modified to read full CMTSOLUTION file particularly for multiple-source cases.
! Time-shifts of all sources can be read and the minimum t_shift is taken to be written in sac headers!
! by Ebru

  subroutine get_event_info_serial(yr,jda,mo,da,ho,mi,sec,event_name,tshift_src,t_shift, &
                                   elat_pde,elon_pde,depth_pde,mb,ms, &
                                   cmt_lat,cmt_lon,cmt_depth,cmt_hdur,NSOURCES, &
                                   Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)

  use constants
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,USE_FORCE_POINT_SOURCE

  implicit none

!--- arguments of the subroutine below

  integer, intent(in) :: NSOURCES

  integer, intent(out) :: yr,jda,mo,da,ho,mi
  double precision, intent(out) :: sec
  double precision, intent(out) :: tshift_src,t_shift
  double precision, intent(out) :: elat_pde,elon_pde,depth_pde
  real, intent(out) :: mb, ms
  double precision, intent(out) :: cmt_lat,cmt_lon,cmt_depth,cmt_hdur
  double precision, intent(out) :: Mrr, Mtt, Mpp, Mrt, Mrp, Mtp

  character(len=20), intent(out) :: event_name ! event name for SAC header

  ! local parameters
  integer :: ios,julian_day
  integer :: isource,idummy
  integer :: i,istart,iend,ier
  double precision, dimension(NSOURCES) :: t_s,hdur,lat,lon,depth
  character(len=20), dimension(NSOURCES) :: e_n
  character(len=5) :: datasource
  character(len=256) :: string

  character(len=MAX_STRING_LEN) :: SOLUTION_FILE, path_to_add

  ! initializes
  yr = 0
  jda = 0
  mo = 0
  da = 0
  ho = 0
  mi = 0
  sec = 0.d0

  event_name = ''

  tshift_src = 0.d0  ! t_cmt_SAC
  t_shift = 0.d0

  elat_pde = 0.d0
  elon_pde = 0.d0
  depth_pde = 0.d0

  mb = 0.0 ! body-wave magnitude
  ms = 0.0 ! surface-wave magnitude

  cmt_lat = -1e8
  cmt_lon = -1e8
  cmt_depth = -1e8
  cmt_hdur = -1e8

  Mrr = 0.d0
  Mtt = 0.d0
  Mpp = 0.d0
  Mrt = 0.d0
  Mrp = 0.d0
  Mtp = 0.d0

  ! reads in source file
  if (USE_FORCE_POINT_SOURCE) then
    ! force solution has no time information
    SOLUTION_FILE = 'DATA/FORCESOLUTION'

    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
      write(path_to_add,"('run',i4.4,'/')") mygroup + 1
      SOLUTION_FILE = path_to_add(1:len_trim(path_to_add))//SOLUTION_FILE(1:len_trim(SOLUTION_FILE))
    endif

    open(unit=IIN,file=trim(SOLUTION_FILE),status='old',action='read',iostat=ios)
    if (ios /= 0) stop 'Error opening DATA/FORCESOLUTION file (in get_event_info_serial)'

    ! read source number isource
    do isource = 1,NSOURCES

      ! example header line of FORCESOLUTION file
      !FORCE 001

      read(IIN,"(a)") string
      ! skips empty lines
      do while (len_trim(string) == 0)
        read(IIN,"(a)") string
      enddo

      ! read header with event information
      read(string,"(a6,i4)") e_n(isource),idummy

      ! read time shift
      read(IIN,"(a)") string
      read(string(12:len_trim(string)),*) t_s(isource)

      read(IIN,"(a)") string
      read(string(15:len_trim(string)),*) hdur(isource)

      ! read latitude
      read(IIN,"(a)") string
      read(string(10:len_trim(string)),*) lat(isource)

      ! read longitude
      read(IIN,"(a)") string
      read(string(11:len_trim(string)),*) lon(isource)

      ! read depth
      read(IIN,"(a)") string
      read(string(7:len_trim(string)),*) depth(isource)

      ! source time function
      read(IIN,"(a)") string
      read(string(22:len_trim(string)),*) idummy ! force_stf(isource)

      ! read magnitude
      read(IIN,"(a)") string
      read(string(21:len_trim(string)),*) mb ! factor_force_source(isource)

      ! read direction vector's East component
      read(IIN,"(a)") string
      read(string(29:len_trim(string)),*) Mpp ! Mcomp_dir_vect_source_E(isource)

      ! read direction vector's North component
      read(IIN,"(a)") string
      read(string(29:len_trim(string)),*) Mtt ! comp_dir_vect_source_N(isource)

      ! read direction vector's vertical component
      read(IIN,"(a)") string
      read(string(32:len_trim(string)),*) Mrr ! comp_dir_vect_source_Z_UP(isource)

    enddo

    close(IIN)

    ! no magnitudes, but force strengths
    ms = mb

  else
    ! CMT sources
    !---- read hypocenter info
    SOLUTION_FILE = 'DATA/CMTSOLUTION'

    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
      write(path_to_add,"('run',i4.4,'/')") mygroup + 1
      SOLUTION_FILE = path_to_add(1:len_trim(path_to_add))//SOLUTION_FILE(1:len_trim(SOLUTION_FILE))
    endif

    open(unit=IIN,file=trim(SOLUTION_FILE),status='old',action='read',iostat=ios)
    if (ios /= 0) stop 'Error opening DATA/CMTSOLUTION file (in get_event_info_serial)'

    ! example header line of CMTSOLUTION file
    !PDE 2003 09 25 19 50 08.93  41.78  144.08  18.0 7.9 8.0 Hokkaido, Japan
    ! which is: event_id, date,origin time,latitude,longitude,depth, mb, MS, region

    ! read source number isource
    do isource = 1,NSOURCES

      ! gets header line
      read(IIN,"(a256)",iostat=ier) string
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading header line in source ',isource
        stop 'Error reading header line in station in CMTSOLUTION file'
      endif

      ! reads header line with event information (free format)
      ! gets rid of the first datasource qualifyer string which can have variable length, like:
      ! "PDE 2014 9 3 .."
      ! " PDEQ2014 9 3 .."
      ! " MLI   1971   1   1 .."
      ! note: globalcmt.org solutions might have missing spaces after datasource qualifier
      !
      ! reads in year,month,day,hour,minutes,seconds
      istart = 1
      iend = len_trim(string)
      ! determines where first number starts
      do i = 1,len_trim(string)
        if (is_numeric(string(i:i))) then
          istart = i
          exit
        endif
      enddo
      if ( istart >= iend ) stop 'Error determining datasource length in header line in CMTSOLUTION file'
      if ( istart <= 1 ) stop 'Error determining datasource length in header line in CMTSOLUTION file'

      ! debug
      !print *,'line ----',string(istart:iend),'----'

      ! read header with event information
      read(string(1:istart-1),*) datasource
      read(string(istart:iend),*) yr,mo,da,ho,mi,sec,elat_pde,elon_pde,depth_pde,mb,ms

      jda = julian_day(yr,mo,da)

      ! read line with event name
      read(IIN,"(a)") string
      read(string(12:len_trim(string)),*) e_n(isource)

      ! read time shift
      read(IIN,"(a)") string
      read(string(12:len_trim(string)),*) t_s(isource)

      ! read half duration
      read(IIN,"(a)") string
      read(string(15:len_trim(string)),*) hdur(isource)

      ! read latitude
      read(IIN,"(a)") string
      read(string(10:len_trim(string)),*) lat(isource)

      ! read longitude
      read(IIN,"(a)") string
      read(string(11:len_trim(string)),*) lon(isource)

      ! read depth
      read(IIN,"(a)") string
      read(string(7:len_trim(string)),*) depth(isource)

      ! read the last 6 lines with moment tensor info
      read(IIN,"(a)") string
      read(string(5:len_trim(string)),*) Mrr
      read(IIN,"(a)") string
      read(string(5:len_trim(string)),*) Mtt
      read(IIN,"(a)") string
      read(string(5:len_trim(string)),*) Mpp
      read(IIN,"(a)") string
      read(string(5:len_trim(string)),*) Mrt
      read(IIN,"(a)") string
      read(string(5:len_trim(string)),*) Mrp
      read(IIN,"(a)") string
      read(string(5:len_trim(string)),*) Mtp

    enddo

    close(IIN)

  endif

  ! sets tshift_src to zero
  tshift_src = 0.d0

  ! takes first event id as event_name
  event_name = e_n(1)

  ! sets cmt info
  if (NSOURCES == 1) then
    cmt_lat = lat(1)
    cmt_lon = lon(1)
    cmt_depth = depth(1)
    cmt_hdur = hdur(1)
    t_shift = t_s(1)
  else
    ! takes minimum time shift of all given sources
    t_shift = minval(t_s(1:NSOURCES))
  endif

contains

  !--------------------------------------------------------------

  logical function is_numeric(char)

  ! returns .true. if input character is a number

  implicit none
  character(len=1), intent(in) :: char

  is_numeric = .false.

  if ( index('0123456789', char) /= 0) then
    is_numeric = .true.
  endif

  end function

  end subroutine get_event_info_serial

