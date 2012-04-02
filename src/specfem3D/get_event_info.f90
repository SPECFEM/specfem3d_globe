!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"

!--- input or output arguments of the subroutine below

  integer, intent(in) :: myrank

  integer, intent(out) :: yr,jda,ho,mi
  real, intent(out) :: mb
  double precision, intent(out) :: tshift_cmt,elat,elon,depth,cmt_lat,cmt_lon,cmt_depth,cmt_hdur,sec

  !character(len=12), intent(out) :: ename

  integer, intent(in) :: NSOURCES ! must be given
  double precision, intent(out) :: t_shift
  character(len=20), intent(out) :: event_name



!--- local variables below

  integer ier

  !integer, parameter :: LENGTH_REGION_NAME = 150
  !character(len=LENGTH_REGION_NAME) region

! get event information for SAC header on the master
  if(myrank == 0) then

    call get_event_info_serial(yr,jda,ho,mi,sec,event_name,tshift_cmt,t_shift, &
                        elat,elon,depth,mb, &
                        cmt_lat,cmt_lon,cmt_depth,cmt_hdur,NSOURCES)

    ! create the event name
    !write(ename(1:12),'(a12)') region(1:12)

    ! replace white spaces with underscores in event name
    !do i=1,len_trim(ename)
    !  if (ename(i:i) == ' ') ename(i:i) = '_'
    !enddo

  endif

! broadcast the information read on the master to the nodes
  call MPI_BCAST(yr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(jda,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(ho,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(mi,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(sec,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

  call MPI_BCAST(NSOURCES,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  call MPI_BCAST(tshift_cmt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(t_shift,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

  ! event location given on first, PDE line
  call MPI_BCAST(elat,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(elon,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(depth,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

  ! cmt location given in CMT file
  call MPI_BCAST(cmt_lat,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(cmt_lon,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(cmt_depth,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(cmt_hdur,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

  !call MPI_BCAST(ename,12,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(event_name,20,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

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

  implicit none

  include "constants.h"

!--- arguments of the subroutine below

  integer, intent(out) :: yr,jda,ho,mi

  real, intent(out) :: mb

  double precision, intent(out) :: sec,tshift_cmt,t_shift
  double precision, intent(out) :: elat_pde,elon_pde,depth_pde,cmt_lat,cmt_lon,cmt_depth,cmt_hdur

  !integer, intent(in) :: LENGTH_REGION_NAME
  !character(len=LENGTH_REGION_NAME), intent(out) :: region ! event name for SAC header

  character(len=20), intent(out) :: event_name ! event name for SAC header

  integer, intent(in) :: NSOURCES

!--- local variables here

  integer ios,mo,da,julian_day
  integer isource

  double precision, dimension(NSOURCES) :: t_s,hdur,lat,lon,depth
  character(len=20), dimension(NSOURCES) :: e_n

  real ms

  character(len=5) datasource
  character(len=150) string,CMTSOLUTION


!
!---- read hypocenter info
!
  call get_value_string(CMTSOLUTION, 'solver.CMTSOLUTION','DATA/CMTSOLUTION')

  open(unit=821,file=CMTSOLUTION,iostat=ios,status='old',action='read')
  if(ios /= 0) stop 'error opening CMTSOLUTION file (in get_event_info_serial)'

  ! example header line of CMTSOLUTION file
  !PDE 2003 09 25 19 50 08.93  41.78  144.08  18.0 7.9 8.0 Hokkaido, Japan
  ! which is: event_id, date,origin time,latitude,longitude,depth, mb, MS, region

  ! read source number isource
  do isource=1,NSOURCES

    ! read header with event information
    read(821,*) datasource,yr,mo,da,ho,mi,sec,elat_pde,elon_pde,depth_pde,mb,ms
    jda=julian_day(yr,mo,da)

    ! ignore line with event name
    read(821,"(a)") string
    read(string(12:len_trim(string)),*) e_n(isource)

    ! read time shift
    read(821,"(a)") string
    read(string(12:len_trim(string)),*) t_s(isource)

    ! read half duration
    read(821,"(a)") string
    read(string(15:len_trim(string)),*) hdur(isource)

    ! read latitude
    read(821,"(a)") string
    read(string(10:len_trim(string)),*) lat(isource)

    ! read longitude
    read(821,"(a)") string
    read(string(11:len_trim(string)),*) lon(isource)

    ! read depth
    read(821,"(a)") string
    read(string(7:len_trim(string)),*) depth(isource)

    ! ignore the last 6 lines with moment tensor info
    read(821,"(a)") string
    read(821,"(a)") string
    read(821,"(a)") string
    read(821,"(a)") string
    read(821,"(a)") string
    read(821,"(a)") string
  enddo
  ! sets tshift_cmt to zero
  tshift_cmt = 0.

  ! takes first event id as event_name
  event_name = e_n(1)

  ! sets cmt infos
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

  close(821)



!  ! read header with event information
!  read(821,*) datasource,yr,mo,da,ho,mi,sec,elat,elon,depth,mb,ms,region
!
!  jda=julian_day(yr,mo,da)
!
!  ! ignore line with event name
!  read(821,"(a)") string
!
!  ! read time shift
!  read(821,"(a)") string
!  read(string(12:len_trim(string)),*) tshift_cmt
!
!  if (NSOURCES == 1) then
!
!  ! read half duration
!    read(821,"(a)") string
!    read(string(15:len_trim(string)),*) cmt_hdur
!
!  ! read latitude
!    read(821,"(a)") string
!    read(string(10:len_trim(string)),*) cmt_lat
!
!  ! read longitude
!    read(821,"(a)") string
!    read(string(11:len_trim(string)),*) cmt_lon
!
!  ! read depth
!    read(821,"(a)") string
!    read(string(7:len_trim(string)),*) cmt_depth
!
!  else
!
!    cmt_hdur=-1e8
!    cmt_lat=-1e8
!    cmt_lon=-1e8
!    cmt_depth=-1e8
!
!  endif
!
!  close(821)

  end subroutine get_event_info_serial

