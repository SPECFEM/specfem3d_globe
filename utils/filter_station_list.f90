!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
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

! program to filter the STATIONS file to include stations located in a given region only
! Dimitri Komatitsch, University of Pau, France, May 2006

  program station_filter

  implicit none

  include '../constants.h'

! input
  character(len=100),parameter :: filename='STATIONS',filtered_filename='STATIONS_FILTERED'
  double precision LATITUDE_MIN,LATITUDE_MAX,LONGITUDE_MIN,LONGITUDE_MAX

! output
  integer irec, nrec, nrec_filtered, ios

  double precision stlat,stlon,stele,stbur
  character(len=MAX_LENGTH_STATION_NAME) station_name
  character(len=MAX_LENGTH_NETWORK_NAME) network_name

  print *
  print *,'program to filter the STATIONS file to include stations located in a given region only'
  print *

  print *,'enter LONGITUDE_MIN:'
  read(*,*) LONGITUDE_MIN
  print *

  print *,'enter LONGITUDE_MAX:'
  read(*,*) LONGITUDE_MAX
  print *

  if(LONGITUDE_MIN >= LONGITUDE_MAX) stop 'incorrect longitude range given'

  print *,'enter LATITUDE_MIN:'
  read(*,*) LATITUDE_MIN
  print *

  print *,'enter LATITUDE_MAX:'
  read(*,*) LATITUDE_MAX
  print *

  if(LATITUDE_MIN >= LATITUDE_MAX) stop 'incorrect latitude range given'

  nrec_filtered = 0

  open(unit=IIN, file=trim(filename), status = 'old', iostat = ios)
  if (ios /= 0) stop 'Input STATIONS file does not exist, exiting'
  read(IIN, *) nrec
  do irec = 1, nrec
    read(IIN, *) station_name, network_name, stlat, stlon, stele, stbur
    if(stlat >= LATITUDE_MIN .and. stlat <= LATITUDE_MAX .and. stlon >= LONGITUDE_MIN .and. stlon <= LONGITUDE_MAX) &
          nrec_filtered = nrec_filtered + 1
  enddo
  close(IIN)

    open(unit=IIN,file=trim(filename),status='old')
    open(unit=IOUT,file=trim(filtered_filename),status='unknown')
    read(IIN,*) nrec
    write(IOUT,*) nrec_filtered
    do irec = 1,nrec
      read(IIN,*) station_name,network_name,stlat,stlon,stele,stbur
      if(stlat >= LATITUDE_MIN .and. stlat <= LATITUDE_MAX .and. stlon >= LONGITUDE_MIN .and. stlon <= LONGITUDE_MAX) &
            write(IOUT,*) trim(station_name),' ',trim(network_name),' ',sngl(stlat),' ',sngl(stlon), ' ',sngl(stele), ' ',sngl(stbur)
    enddo
    close(IIN)
    close(IOUT)

    print *
    print *,'there are ',nrec,' stations in file ', trim(filename)
    print *,'saving ',nrec_filtered,' stations inside the model in file ', trim(filtered_filename)
    print *,'excluding ',nrec - nrec_filtered,' stations located outside the model'
    print *

  end program station_filter

