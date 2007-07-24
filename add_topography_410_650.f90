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

  subroutine add_topography_410_650(myrank,xelm,yelm,zelm,R220,R400,R670,R771)

  implicit none

  include "constants.h"

  integer myrank

  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)

  double precision R220,R400,R670,R771

  integer ia

  real xcolat,xlon
  real topo410out,topo650out
  double precision topo410,topo650

  double precision r,theta,phi
  double precision gamma

! we loop on all the points of the element
  do ia = 1,NGNOD

! convert to r theta phi
    call xyz_2_rthetaphi_dble(xelm(ia),yelm(ia),zelm(ia),r,theta,phi)
    call reduce(theta,phi)

! get colatitude and longitude in degrees
    xcolat = sngl(theta*180.0d0/PI)
    xlon = sngl(phi*180.0d0/PI)

! compute topography on 410 and 650 at current point
    call subtopo(xcolat,xlon,topo410out,topo650out)

! non-dimensionalize the topography, which is in km
! positive for a depression, so change the sign for a perturbation in radius
    topo410 = -dble(topo410out) / R_EARTH_KM
    topo650 = -dble(topo650out) / R_EARTH_KM

    if(r >= R400/R_EARTH .and. r <= R220/R_EARTH) then
! stretching between R220 and R400
      gamma = (R220/R_EARTH - r) / (R220/R_EARTH - R400/R_EARTH)
      xelm(ia) = xelm(ia)*(ONE + gamma * topo410 / r)
      yelm(ia) = yelm(ia)*(ONE + gamma * topo410 / r)
      zelm(ia) = zelm(ia)*(ONE + gamma * topo410 / r)
    elseif(r>= R771/R_EARTH .and. r <= R670/R_EARTH) then
! stretching between R771 and R670
      gamma = (r - R771/R_EARTH) / (R670/R_EARTH - R771/R_EARTH)
      xelm(ia) = xelm(ia)*(ONE + gamma * topo650 / r)
      yelm(ia) = yelm(ia)*(ONE + gamma * topo650 / r)
      zelm(ia) = zelm(ia)*(ONE + gamma * topo650 / r)
    elseif(r > R670/R_EARTH .and. r < R400/R_EARTH) then
! stretching between R670 and R400
      gamma = (R400/R_EARTH - r) / (R400/R_EARTH - R670/R_EARTH)
      xelm(ia) = xelm(ia)*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
      yelm(ia) = yelm(ia)*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
      zelm(ia) = zelm(ia)*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
    endif
    if(gamma < -0.0001 .or. gamma > 1.0001) call exit_MPI(myrank,'incorrect value of gamma for 410-650 topography')

  enddo

  end subroutine add_topography_410_650

