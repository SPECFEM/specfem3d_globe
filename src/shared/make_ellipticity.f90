!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

  subroutine make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)

! creates a spline for the ellipticity profile in PREM
! radius and density are non-dimensional

  use constants, only: NR,TWO_PI,PI,GRAV,RHOAV,HOURS_PER_DAY,SECONDS_PER_HOUR,R_UNIT_SPHERE

  implicit none

  integer :: nspl

  double precision,dimension(NR) :: rspl,espl,espl2

  logical :: ONE_CRUST

  ! local parameters
  integer :: i
  ! radii
  double precision :: ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R220,R400,R600,R670, &
                   R771,RTOPDDOUBLEPRIME,RCMB,RICB
  double precision :: r_icb,r_cmb,r_topddoubleprime,r_771,r_670,r_600
  double precision :: r_400,r_220,r_80,r_moho,r_middle_crust,r_ocean,r_0

  double precision,dimension(NR) :: r,rho,epsilonval,eta
  double precision,dimension(NR) :: radau,k
  double precision,dimension(NR) :: s1,s2,s3

  double precision :: z,g_a,bom,exponentval,i_rho,i_radau
  double precision :: yp1,ypn

  ! radius of the Earth for gravity calculation
  double precision, parameter :: R_EARTH_ELLIPTICITY = 6371000.d0
  ! radius of the ocean floor for gravity calculation
  double precision, parameter :: ROCEAN_ELLIPTICITY = 6368000.d0

  ! PREM
  ROCEAN = 6368000.d0
  RMIDDLE_CRUST = 6356000.d0
  RMOHO = 6346600.d0
  R80  = 6291000.d0
  R220 = 6151000.d0
  R400 = 5971000.d0
  R600 = 5771000.d0
  R670 = 5701000.d0
  R771 = 5600000.d0
  RTOPDDOUBLEPRIME = 3630000.d0
  RCMB = 3480000.d0
  RICB = 1221000.d0

  ! non-dimensionalize
  r_icb = RICB/R_EARTH_ELLIPTICITY
  r_cmb = RCMB/R_EARTH_ELLIPTICITY
  r_topddoubleprime = RTOPDDOUBLEPRIME/R_EARTH_ELLIPTICITY
  r_771 = R771/R_EARTH_ELLIPTICITY
  r_670 = R670/R_EARTH_ELLIPTICITY
  r_600 = R600/R_EARTH_ELLIPTICITY
  r_400 = R400/R_EARTH_ELLIPTICITY
  r_220 = R220/R_EARTH_ELLIPTICITY
  r_80 = R80/R_EARTH_ELLIPTICITY
  r_moho = RMOHO/R_EARTH_ELLIPTICITY
  r_middle_crust = RMIDDLE_CRUST/R_EARTH_ELLIPTICITY
  r_ocean = ROCEAN_ELLIPTICITY/R_EARTH_ELLIPTICITY
  r_0 = 1.d0

  do i=1,163
    r(i) = r_icb*dble(i-1)/dble(162)
  enddo
  do i=164,323
    r(i) = r_icb+(r_cmb-r_icb)*dble(i-164)/dble(159)
  enddo
  do i=324,336
    r(i) = r_cmb+(r_topddoubleprime-r_cmb)*dble(i-324)/dble(12)
  enddo
  do i=337,517
    r(i) = r_topddoubleprime+(r_771-r_topddoubleprime)*dble(i-337)/dble(180)
  enddo
  do i=518,530
    r(i) = r_771+(r_670-r_771)*dble(i-518)/dble(12)
  enddo
  do i=531,540
    r(i) = r_670+(r_600-r_670)*dble(i-531)/dble(9)
  enddo
  do i=541,565
    r(i) = r_600+(r_400-r_600)*dble(i-541)/dble(24)
  enddo
  do i=566,590
    r(i) = r_400+(r_220-r_400)*dble(i-566)/dble(24)
  enddo
  do i=591,609
    r(i) = r_220+(r_80-r_220)*dble(i-591)/dble(18)
  enddo
  do i=610,619
    r(i) = r_80+(r_moho-r_80)*dble(i-610)/dble(9)
  enddo
  do i=620,626
    r(i) = r_moho+(r_middle_crust-r_moho)*dble(i-620)/dble(6)
  enddo
  do i=627,633
    r(i) = r_middle_crust+(r_ocean-r_middle_crust)*dble(i-627)/dble(6)
  enddo
  do i=634,NR
    r(i) = r_ocean+(r_0-r_ocean)*dble(i-634)/dble(6)
  enddo

  ! use PREM to get the density profile for ellipticity (fine for other 1D reference models)
  do i = 1,NR
    call prem_density(r(i),rho(i),ONE_CRUST,RICB,RCMB,RTOPDDOUBLEPRIME, &
                      R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)
    radau(i) = rho(i)*r(i)*r(i)
  enddo

  eta(1) = 0.0d0

  k(1) = 0.0d0

  do i=2,NR
    call intgrl(i_rho,r,1,i,rho,s1,s2,s3)
! Radau approximation of Clairaut's equation for first-order terms of ellipticity, see e.g. Jeffreys H.,
! The figures of rotating planets, Mon. Not. R. astr. Soc., vol. 113, p. 97-105 (1953).
! The Radau approximation is mentioned on page 97.
! For more details see Section 14.1.2 in Dahlen and Tromp (1998)
! (see also in file ellipticity_equations_from_Dahlen_Tromp_1998.pdf in the "doc" directory of the code).
    call intgrl(i_radau,r,1,i,radau,s1,s2,s3)
    z=(2.0d0/3.0d0)*i_radau/(i_rho*r(i)*r(i))
    ! this comes from equation (14.19) in Dahlen and Tromp (1998)
    eta(i)=(25.0d0/4.0d0)*((1.0d0-(3.0d0/2.0d0)*z)**2.0d0)-1.0d0
    k(i)=eta(i)/(r(i)**3.0d0)
  enddo

  bom = TWO_PI/(HOURS_PER_DAY*SECONDS_PER_HOUR)

  ! non-dimensionalized value
  bom = bom/sqrt(PI*GRAV*RHOAV)

  g_a = 4.0d0*i_rho
  ! this is the equation right above (14.21) in Dahlen and Tromp (1998)
  epsilonval(NR) = (5.0d0/2.d0)*(bom**2.0d0)*R_UNIT_SPHERE / (g_a * (eta(NR)+2.0d0))

  do i = 1,NR-1
    call intgrl(exponentval,r,i,NR,k,s1,s2,s3)
    epsilonval(i) = epsilonval(NR)*exp(-exponentval)
  enddo

  ! initializes spline coefficients
  rspl(:) = 0.d0
  espl(:) = 0.d0
  espl2(:) = 0.d0

  ! get ready to spline epsilonval
  nspl = 1
  rspl(1) = r(1)
  espl(1) = epsilonval(1)
  do i = 2,NR
    if (r(i) /= r(i-1)) then
      nspl = nspl+1
      rspl(nspl) = r(i)
      espl(nspl) = epsilonval(i)
    endif
  enddo

  ! spline epsilonval
  yp1 = 0.0d0
  ypn = (5.0d0/2.0d0)*(bom**2)/g_a-2.0d0*epsilonval(NR)
  call spline_construction(rspl,espl,nspl,yp1,ypn,espl2)

  end subroutine make_ellipticity

!
!-----------------------------------------------------------------
!

  subroutine revert_ellipticity(x,y,z,nspl,rspl,espl,espl2)

! this routine to revert ellipticity and go back to a spherical Earth
! is currently used by src/auxiliaries/combine_vol_data.F90 only

  use constants

  implicit none

  real(kind=CUSTOM_REAL) :: x,y,z
  integer nspl
  double precision rspl(NR),espl(NR),espl2(NR)
  double precision x1,y1,z1

  double precision ell
  double precision r,theta,phi,factor
  double precision cost,p20

  ! gets spherical coordinates
  x1 = x
  y1 = y
  z1 = z
  call xyz_2_rthetaphi_dble(x1,y1,z1,r,theta,phi)

  cost = dcos(theta)
! this is the Legendre polynomial of degree two, P2(cos(theta)), see the discussion above eq (14.4) in Dahlen and Tromp (1998)
  p20 = 0.5d0*(3.0d0*cost*cost-1.0d0)

  ! get ellipticity using spline evaluation
  call spline_evaluation(rspl,espl,espl2,nspl,r,ell)

! this is eq (14.4) in Dahlen and Tromp (1998)
  factor = ONE-(TWO/3.0d0)*ell*p20

  ! removes ellipticity factor
  x = x / factor
  y = y / factor
  z = z / factor

  end subroutine revert_ellipticity

