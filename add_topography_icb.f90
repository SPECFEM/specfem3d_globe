!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                    and University of Pau, France
! (c) California Institute of Technology and University of Pau, October 2007
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine add_topography_icb(myrank,xelm,yelm,zelm,RICB,RCMB)

  implicit none

  include "constants.h"

  integer myrank

  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)

  double precision RICB,RCMB

  integer ia

  double precision topoicb

  double precision r,theta,phi
  double precision gamma

! we loop on all the points of the element
  do ia = 1,NGNOD

! convert to r theta phi
    call xyz_2_rthetaphi_dble(xelm(ia),yelm(ia),zelm(ia),r,theta,phi)
    call reduce(theta,phi)

! compute topography on ICB; the routine subtopo_icb needs to be supplied by the user
!    call subtopo_icb(theta,phi,topoicb)
    topoicb = 0.0d0

! non-dimensionalize the topography, which is in km
! positive for a depression, so change the sign for a perturbation in radius
    topoicb = -topoicb / R_EARTH_KM

    gamma = 0.0d0
    if(r > 0.0d0 .and. r <= RICB/R_EARTH) then
! stretching between center and RICB
      gamma = r/(RICB/R_EARTH)
    elseif(r>= RICB/R_EARTH .and. r <= RCMB/R_EARTH) then
! stretching between RICB and RCMB
      gamma = (r - RCMB/R_EARTH) / (RICB/R_EARTH - RCMB/R_EARTH)
    endif
    if(gamma < -0.0001 .or. gamma > 1.0001) call exit_MPI(myrank,'incorrect value of gamma for CMB topography')

    xelm(ia) = xelm(ia)*(ONE + gamma * topoicb / r)
    yelm(ia) = yelm(ia)*(ONE + gamma * topoicb / r)
    zelm(ia) = zelm(ia)*(ONE + gamma * topoicb / r)

  enddo

  end subroutine add_topography_icb

