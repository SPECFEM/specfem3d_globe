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
! every 0.1 degree along one-fourth of a great circle

! Dimitri Komatitsch, University of Pau, France, October 2007

  program generate_station_list

  implicit none

  include "constants.h"

  integer nrec,ilatitude,idecimal_degree,upper_bound

  double precision stlat,decimal_degree

  character(len=MAX_LENGTH_STATION_NAME) station_name
  character(len=MAX_LENGTH_NETWORK_NAME) network_name

  print *
  print *,'generating the DATA/STATIONS file to include stations'
  print *,'every 0.1 degree along one-fourth of a great circle'
  print *

  network_name = 'DK'

  nrec = 451

  open(unit=IOUT,file='DATA/STATIONS',status='unknown')

  write(IOUT,*) nrec

! generate one-fourth of the great circle along the Greenwich meridian (longitude = 0)
  do ilatitude = -90,0

! add 1/10th of a degree to generate stations every 0.1 degree between -90 and -50
! otherwise only every degree
    if(ilatitude < -50) then
      upper_bound = 9
    else
      upper_bound = 0
    endif

    do idecimal_degree = 0,upper_bound

      stlat = ilatitude + idecimal_degree / 10.d0
      write(station_name,"('FRONT_MINUS_',i2.2,'_',i2.2)") int(abs(stlat)),nint(100. * (abs(stlat) - int(abs(stlat))))
      write(IOUT,*) trim(station_name),' ',trim(network_name),' ',sngl(stlat),' 0 0 0'

    enddo

  enddo

  close(IOUT)

  end program generate_station_list

