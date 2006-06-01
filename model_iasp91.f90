!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 5
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine model_iasp91(myrank,x,rho,vp,vs,Qkappa,Qmu,idoubling, &
      check_doubling_flag,RICB,RCMB,RTOPDDOUBLEPRIME,R670,R220,R771,R400,RMOHO)

  implicit none

  include "constants.h"

! given a normalized radius x, gives the non-dimesionalized density rho,
! speeds vp and vs, and the quality factors Qkappa and Qmu

  logical check_doubling_flag

  integer idoubling,myrank

  double precision x,rho,vp,vs,Qkappa,Qmu,RICB,RCMB,RTOPDDOUBLEPRIME,R670,R220,R771,R400,RMOHO

  double precision r,scaleval

  double precision, parameter :: r120 = 6251.d3,r020 = 6351.d3
  double precision x1,x2

! compute real physical radius in meters
  r = x * R_EARTH

  x1 = r120 / R_EARTH
  x2 = RMOHO / R_EARTH

! check flags to make sure we correctly honor the discontinuities
! we use strict inequalities since r has been slighly changed in mesher

  if(check_doubling_flag) then

!
!--- inner core
!
  if(r >= 0.d0 .and. r < RICB) then
    if(idoubling /= IFLAG_TOP_INNER_CORE .and. &
       idoubling /= IFLAG_INNER_CORE_NORMAL .and. &
       idoubling /= IFLAG_IN_CENTRAL_CUBE .and. &
       idoubling /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
       idoubling /= IFLAG_TOP_CENTRAL_CUBE .and. &
       idoubling /= IFLAG_IN_FICTITIOUS_CUBE) &
         call exit_MPI(myrank,'wrong doubling flag for inner core point')
!
!--- outer core
!
  else if(r > RICB .and. r < RCMB) then
    if(idoubling /= IFLAG_OUTER_CORE_NORMAL .and. idoubling /= IFLAG_TOP_OUTER_CORE &
    .and. idoubling /= IFLAG_BOTTOM_OUTER_CORE) &
      call exit_MPI(myrank,'wrong doubling flag for outer core point')
!
!--- D" at the base of the mantle
!
  else if(r > RCMB .and. r < RTOPDDOUBLEPRIME) then
    if(idoubling /= IFLAG_BOTTOM_MANTLE .and. idoubling /= IFLAG_MANTLE_NORMAL) &
      call exit_MPI(myrank,'wrong doubling flag for D" point')
!
!--- mantle: from top of D" to d670
!
  else if(r > RTOPDDOUBLEPRIME .and. r < R670) then
    if(idoubling /= IFLAG_MANTLE_NORMAL .and. idoubling /= IFLAG_DOUBLING_670) &
      call exit_MPI(myrank,'wrong doubling flag for top D" -> d670 point')

!
!--- mantle: from d670 to d220
!
  else if(r > R670 .and. r < R220) then
    if(idoubling /= IFLAG_670_220) &
      call exit_MPI(myrank,'wrong doubling flag for d670 -> d220 point')

!
!--- mantle and crust: from d220 to MOHO and then to surface
!
  else if(r > R220) then
    if(idoubling /= IFLAG_220_MOHO .and. idoubling /= IFLAG_CRUST) &
      call exit_MPI(myrank,'wrong doubling flag for d220 -> Moho -> surface point')

  endif

  endif

!
!--- inner core
!
  if(r >= 0.d0 .and. r <= RICB) then
    rho=13.0885d0-8.8381d0*x*x
    vp=11.24094-4.09689*x**2
    vs=3.56454-3.45241*x**2
    Qmu=84.6d0
    Qkappa=1327.7d0
!
!--- outer core
!
  else if(r > RICB .and. r <= RCMB) then
    rho=12.5815d0-1.2638d0*x-3.6426d0*x*x-5.5281d0*x*x*x
    vp=10.03904+3.75665*x-13.67046*x**2
    vs=0.0d0
    Qmu=0.0d0
    Qkappa=57827.0d0
!
!--- D" at the base of the mantle
!
  else if(r > RCMB .and. r <= RTOPDDOUBLEPRIME) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vp=14.49470-1.47089*x
    vs=8.16616-1.58206*x
    Qmu=312.0d0
    Qkappa=57827.0d0
!
!--- mantle: from top of D" to d670
!
  else if(r > RTOPDDOUBLEPRIME .and. r <= R771) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vp=25.1486-41.1538*x+51.9932*x**2-26.6083*x**3
    vs=12.9303-21.2590*x+27.8988*x**2-14.1080*x**3
    Qmu=312.0d0
    Qkappa=57827.0d0
  else if(r > R771 .and. r <= R670) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vp=25.96984-16.93412*x
    vs=20.76890-16.53147*x
    Qmu=312.0d0
    Qkappa=57827.0d0
!
!--- mantle: above d670
!
  else if(r > R670 .and. r <= R400) then
    rho=5.3197d0-1.4836d0*x
    vp=29.38896-21.40656*x
    vs=17.70732-13.50652*x
    Qmu=143.0d0
    Qkappa=57827.0d0
  else if(r > R400 .and. r <= R220) then
    rho=7.1089d0-3.8045d0*x
    vp=30.78765-23.25415*x
    vs=15.24213-11.08552*x
    Qmu=143.0d0
    Qkappa=57827.0d0

! from Sebastien Chevrot: for the IASP91 model
! Depth        R                Vp                    Vs
! 0-20       6351-6371         5.80                  3.36
! 20-35      6336-6351         6.50                  3.75
! 35-120     6251-6336   8.78541-0.74953 x       6.706231-2.248585 x
! with x = r / 6371

  else if(r > R220 .and. r <= r120) then
    rho=2.6910d0+0.6924d0*x
    vp=25.41389-17.69722*x
    vs=5.75020-1.27420*x
    Qmu=80.0d0
    Qkappa=57827.0d0

  else if(r > r120 .and. r <= RMOHO) then
      vp = 8.78541d0-0.74953d0*x
      vs = 6.706231d0-2.248585d0*x
      rho = 3.3713d0 + (3.3198d0-3.3713d0)*(x-x1)/(x2-x1)
      if(rho < 3.30d0 .or. rho > 3.38d0) stop 'incorrect density computed for IASP91'
      Qmu=600.0d0
      Qkappa=57827.0d0

  else if(r > RMOHO .and. r <= r020) then
      vp = 6.5d0
      vs = 3.75d0
      rho = 2.92d0
      Qmu=600.0d0
      Qkappa=57827.0d0
  else
      vp = 5.8d0
      vs = 3.36d0
      rho = 2.72d0
      Qmu=600.0d0
      Qkappa=57827.0d0
  endif

! non-dimensionalize
! time scaling (s^{-1}) is done with scaleval
  scaleval=dsqrt(PI*GRAV*RHOAV)
  rho=rho*1000.0d0/RHOAV
  vp=vp*1000.0d0/(R_EARTH*scaleval)
  vs=vs*1000.0d0/(R_EARTH*scaleval)

  end subroutine model_iasp91

