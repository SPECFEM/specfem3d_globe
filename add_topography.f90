!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology September 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine add_topography(myrank,xelm,yelm,zelm,ibathy_topo)

  implicit none

  include "constants.h"

  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)

  integer myrank

! use integer array to store values
  integer ibathy_topo(NX_BATHY,NY_BATHY)

  integer ia

  double precision lat,lon,elevation
  double precision r,theta,phi,colat
  double precision gamma

! we loop on all the points of the element
  do ia = 1,NGNOD

! convert to r theta phi
! slightly move points to avoid roundoff problem when exactly on the polar axis
    call xyz_2_rthetaphi_dble(xelm(ia),yelm(ia),zelm(ia),r,theta,phi)
    theta = theta + 0.0000001d0
    phi = phi + 0.0000001d0
    call reduce(theta,phi)

! convert the geocentric colatitude to a geographic colatitude
  colat = PI/2.0d0 - datan(1.006760466d0*dcos(theta)/dmax1(TINYVAL,dsin(theta)))

! get geographic latitude and longitude in degrees
  lat = 90.0d0 - colat*180.0d0/PI
  lon = phi*180.0d0/PI
  elevation = 0.d0

! compute elevation at current point
  call get_topo_bathy(lat,lon,elevation,ibathy_topo)

! non-dimensionalize the elevation, which is in meters
  elevation = elevation / R_EARTH

! stretching topography between d220 and the surface
  gamma = (r - R220/R_EARTH) / (R_UNIT_SPHERE - R220/R_EARTH)

! add elevation to all the points of that element
! also make sure gamma makes sense
  if(gamma < -0.02 .or. gamma > 1.02) call exit_MPI(myrank,'incorrect value of gamma for topography')

  xelm(ia) = xelm(ia)*(ONE + gamma * elevation / r)
  yelm(ia) = yelm(ia)*(ONE + gamma * elevation / r)
  zelm(ia) = zelm(ia)*(ONE + gamma * elevation / r)

  enddo

  end subroutine add_topography

