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

  subroutine get_ellipticity(xelm,yelm,zelm,nspl,rspl,espl,espl2)

  implicit none

  include "constants.h"

  integer nspl
  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)
  double precision rspl(NR),espl(NR),espl2(NR)

  integer ia

  double precision ell
  double precision r,theta,phi,factor
  double precision cost,p20

  do ia=1,NGNOD

  call xyz_2_rthetaphi_dble(xelm(ia),yelm(ia),zelm(ia),r,theta,phi)

  cost=dcos(theta)
  p20=0.5d0*(3.0d0*cost*cost-1.0d0)

! get ellipticity using spline evaluation
  call splint(rspl,espl,espl2,nspl,r,ell)

  factor=ONE-(TWO/3.0d0)*ell*p20

  xelm(ia)=xelm(ia)*factor
  yelm(ia)=yelm(ia)*factor
  zelm(ia)=zelm(ia)*factor

  enddo

  end subroutine get_ellipticity

