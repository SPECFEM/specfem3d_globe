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

! get information about event name and location for SAC seismograms: MPI version by Dimitri Komatitsch

! Instead of using region names as event names,
! event names given in the second row of CMT files will be used.
! Thus, I removed old parameters ename, region, LENGTH_REGION_NAME and added event_name!!!!!!!
! Also, t_shift is added as a new parameter to be written on sac headers!
! by Ebru Bozdag

  subroutine get_event_info_parallel(myrank,yr,jda,ho,mi,sec,&
                                    event_name,tshift_cmt,t_shift, &
                                    elat,elon,depth,mb,cmt_lat, &
                                    cmt_lon,cmt_depth,cmt_hdur,NSOURCES)

  use constants

  implicit none

!--- input or output arguments of the subroutine below

  integer, intent(in) :: myrank
  integer, intent(in) :: NSOURCES ! must be given

  integer, intent(out) :: yr,jda,ho,mi
  double precision, intent(out) :: sec
  real, intent(out) :: mb
  double precision, intent(out) :: tshift_cmt,elat,elon,depth,cmt_lat,cmt_lon,cmt_depth,cmt_hdur
  double precision, intent(out) :: t_shift

  character(len=20), intent(out) :: event_name

  ! get event information for SAC header on the master
  if (myrank == 0) then

    ! note: mb as (body wave) moment magnitude is not used any further,
    !       see comment in write_output_SAC() routine

    call get_event_info_serial(yr,jda,ho,mi,sec,event_name,tshift_cmt,t_shift, &
                        elat,elon,depth,mb, &
                        cmt_lat,cmt_lon,cmt_depth,cmt_hdur,NSOURCES)

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
  call bcast_all_singlei(ho)
  call bcast_all_singlei(mi)

  call bcast_all_singledp(sec)

  call bcast_all_singledp(tshift_cmt)
  call bcast_all_singledp(t_shift)

  ! event location given on first, PDE line
  call bcast_all_singledp(elat)
  call bcast_all_singledp(elon)
  call bcast_all_singledp(depth)

  ! cmt location given in CMT file
  call bcast_all_singledp(cmt_lat)
  call bcast_all_singledp(cmt_lon)
  call bcast_all_singledp(cmt_depth)
  call bcast_all_singledp(cmt_hdur)

  call bcast_all_ch(event_name,20)

  end subroutine get_event_info_parallel

!=====================================================================

! get information about event name and location for SAC seismograms: MPI version by Bernhard Schuberth
! This subroutine reads the first line of the DATA/CMTSOLUTION file
! and extracts event information needed for SAC or PITSA headers

! This subroutine has been modified to read full CMTSOLUTION file particularly for multiple-source cases.
! Time-shifts of all sources can be read and the minimum t_shift is taken to be written in sac headers!
! by Ebru

  subroutine get_event_info_serial(yr,jda,ho,mi,sec,event_name,tshift_cmt,t_shift,&
                            elat_pde,elon_pde,depth_pde,mb,&
                            cmt_lat,cmt_lon,cmt_depth,cmt_hdur,NSOURCES)

  use constants
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS

  implicit none

!--- arguments of the subroutine below

  integer, intent(in) :: NSOURCES

  integer, intent(out) :: yr,jda,ho,mi
  double precision, intent(out) :: sec
  double precision, intent(out) :: tshift_cmt,t_shift
  double precision, intent(out) :: elat_pde,elon_pde,depth_pde
  real, intent(out) :: mb
  double precision, intent(out) :: cmt_lat,cmt_lon,cmt_depth,cmt_hdur

  character(len=20), intent(out) :: event_name ! event name for SAC header

  ! local parameters
  integer :: ios,mo,da,julian_day
  integer :: isource
  double precision, dimension(NSOURCES) :: t_s,hdur,lat,lon,depth
  character(len=20), dimension(NSOURCES) :: e_n
  real :: ms
  character(len=5) :: datasource
  character(len=256) :: string

  character(len=MAX_STRING_LEN) :: CMTSOLUTION, path_to_add

!
!---- read hypocenter info
  CMTSOLUTION = 'DATA/CMTSOLUTION'

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    CMTSOLUTION=path_to_add(1:len_trim(path_to_add))//CMTSOLUTION(1:len_trim(CMTSOLUTION))
  endif
!
  open(unit=IIN,file=trim(CMTSOLUTION),status='old',action='read',iostat=ios)
  if (ios /= 0) stop 'Error opening DATA/CMTSOLUTION file (in get_event_info_serial)'

  ! example header line of CMTSOLUTION file
  !PDE 2003 09 25 19 50 08.93  41.78  144.08  18.0 7.9 8.0 Hokkaido, Japan
  ! which is: event_id, date,origin time,latitude,longitude,depth, mb, MS, region

  ! read source number isource
  do isource = 1,NSOURCES

    ! read header with event information
    read(IIN,*) datasource,yr,mo,da,ho,mi,sec,elat_pde,elon_pde,depth_pde,mb,ms
    jda=julian_day(yr,mo,da)

    ! ignore line with event name
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

    ! ignore the last 6 lines with moment tensor info
    read(IIN,"(a)") string
    read(IIN,"(a)") string
    read(IIN,"(a)") string
    read(IIN,"(a)") string
    read(IIN,"(a)") string
    read(IIN,"(a)") string
  enddo

  ! sets tshift_cmt to zero
  tshift_cmt = 0.

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
    cmt_lat = -1e8
    cmt_lon = -1e8
    cmt_depth = -1e8
    cmt_hdur = -1e8
    ! takes minimum time shift of all given sources
    t_shift = minval(t_s(1:NSOURCES))
  endif

  close(IIN)

  end subroutine get_event_info_serial

