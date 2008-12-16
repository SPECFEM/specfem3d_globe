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

! get information about event name and location for SAC seismograms: MPI version by Dimitri Komatitsch

  subroutine get_event_info_parallel(myrank,yr,jda,ho,mi,sec,t_cmt, &
                 elat,elon,depth,mb,ename,cmt_lat,cmt_lon,cmt_depth,cmt_hdur,NSOURCES)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"

!--- input or output arguments of the subroutine below

  integer, intent(in) :: myrank

  integer, intent(out) :: NSOURCES,yr,jda,ho,mi
  real, intent(out) :: mb
  double precision, intent(out) :: t_cmt,elat,elon,depth,cmt_lat,cmt_lon,cmt_depth,cmt_hdur,sec
  character(len=12), intent(out) :: ename

!--- local variables below

  integer i,ier

  integer, parameter :: LENGTH_REGION_NAME = 150
  character(len=LENGTH_REGION_NAME) region

! get event information for SAC header on the master
  if(myrank == 0) then

    call get_event_info_serial(yr,jda,ho,mi,sec,t_cmt,elat,elon,depth,mb,region, &
                        cmt_lat,cmt_lon,cmt_depth,cmt_hdur,NSOURCES,LENGTH_REGION_NAME)

! create the event name
    write(ename(1:12),'(a12)') region(1:12)

! replace white spaces with underscores in event name
    do i=1,len_trim(ename)
      if (ename(i:i) == ' ') ename(i:i) = '_'
    enddo

  endif

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

  call MPI_BCAST(ename,12,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

  end subroutine get_event_info_parallel

!=====================================================================

! get information about event name and location for SAC seismograms: MPI version by Bernhard Schuberth
! This subroutine reads the first line of the DATA/CMTSOLUTION file
! and extracts event information needed for SAC or PITSA headers

  subroutine get_event_info_serial(yr,jda,ho,mi,sec,t_cmt,elat,elon,depth,mb,region,&
                            cmt_lat,cmt_lon,cmt_depth,cmt_hdur,NSOURCES,LENGTH_REGION_NAME)

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
  if(ios /= 0) stop 'error opening CMTSOLUTION file (in get_event_info_serial)'

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
    read(821,"(a4,i5,i3,i3,i3,i3,f6.2,a)") datasource,yr,mo,da,ho,mi,sec,dummystring
    read(dummystring,*) elat,elon,depth,mb,ms,region

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

  end subroutine get_event_info_serial

