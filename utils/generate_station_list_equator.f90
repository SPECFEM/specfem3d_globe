!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
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

! program to generate the STATIONS file to include stations
! every half-degree along a great circle

! Dimitri Komatitsch, University of Pau, France, August 2007

  program generate_station_list

  implicit none

  include "constants.h"

  integer :: nrec,ilongitude

  double precision :: stlat,stlon

  character(len=MAX_LENGTH_STATION_NAME) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME) :: network_name

  print *
  print *,'generating the DATA/STATIONS file to include stations'
  print *,'every half-degree along a great circle'
  print *

  open(unit=IOUT,file='DATA/STATIONS',status='unknown')

! there is a station every half degree
  nrec = 720

  write(IOUT,*) nrec

! select the great circle along the equator
  stlat = 0.d0

  network_name = 'DK'

  do ilongitude = 0,359

! station on integer value of degree YYY
    stlon = ilongitude
    write(station_name,"('S',i3.3,'_00')") ilongitude
    write(IOUT,*) trim(station_name),' ',trim(network_name),' ',sngl(stlat),' ',sngl(stlon),' 0 0'

! station on half value of degree YYY (i.e. in degree YYY + 0.5)
    stlon = ilongitude + 0.5d0
    write(station_name,"('S',i3.3,'_50')") ilongitude
    write(IOUT,*) trim(station_name),' ',trim(network_name),' ',sngl(stlat),' ',sngl(stlon),' 0 0'

  enddo

  close(IOUT)

  end program generate_station_list

