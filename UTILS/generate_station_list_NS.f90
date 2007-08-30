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

! program to generate the STATIONS file to include stations
! every half-degree along a great circle

! Dimitri Komatitsch, University of Pau, France, August 2007

  program generate_station_list

  implicit none

  include "constants.h"

  integer nrec,ilatitude

  double precision stlat

  character(len=MAX_LENGTH_STATION_NAME) station_name
  character(len=MAX_LENGTH_NETWORK_NAME) network_name

  print *
  print *,'generating the DATA/STATIONS file to include stations'
  print *,'every half-degree along a great circle'
  print *

  network_name = 'DK'

! there is a station every half degree
  nrec = 720

  open(unit=IOUT,file='DATA/STATIONS',status='unknown')

  write(IOUT,*) nrec

! first create half the great circle along the Greenwich meridian (longitude = 0)

  do ilatitude = -90,+89

    stlat = ilatitude
    if(ilatitude < 0) then
      write(station_name,"('FRONT_MINUS_',i2.2,'_00')") int(abs(stlat))
    else
      write(station_name,"('FRONT_PLUS_',i2.2,'_00')") int(abs(stlat))
    endif
    write(IOUT,*) trim(station_name),' ',trim(network_name),' ',sngl(stlat),' 0 0 0'

    stlat = ilatitude + 0.5d0
    if(stlat < 0) then
      write(station_name,"('FRONT_MINUS_',i2.2,'_50')") int(abs(stlat))
    else
      write(station_name,"('FRONT_PLUS_',i2.2,'_50')") int(abs(stlat))
    endif
    write(IOUT,*) trim(station_name),' ',trim(network_name),' ',sngl(stlat),' 0 0 0'

  enddo

! then create half the great circle along longitude = 180

  do ilatitude = -90,+89

    stlat = ilatitude
    if(ilatitude < 0) then
      write(station_name,"('BACK_MINUS_',i2.2,'_00')") int(abs(stlat))
    else
      write(station_name,"('BACK_PLUS_',i2.2,'_00')") int(abs(stlat))
    endif
    write(IOUT,*) trim(station_name),' ',trim(network_name),' ',sngl(stlat),' 180 0 0'

    stlat = ilatitude + 0.5d0
    if(stlat < 0) then
      write(station_name,"('BACK_MINUS_',i2.2,'_50')") int(abs(stlat))
    else
      write(station_name,"('BACK_PLUS_',i2.2,'_50')") int(abs(stlat))
    endif
    write(IOUT,*) trim(station_name),' ',trim(network_name),' ',sngl(stlat),' 180 0 0'

  enddo

  close(IOUT)

  end program generate_station_list

