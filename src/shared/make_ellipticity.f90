!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
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


  subroutine make_ellipticity(nspl,rspl,ellipicity_spline,ellipicity_spline2)

! creates a spline for the ellipticity profile in PREM
! radius and density are non-dimensional
! or
! for mars, creates a spline for the ellipticity profile in Mars (Sohl $ Spohn)

  use constants, only: NR_DENSITY

  implicit none

  integer,intent(inout) :: nspl
  double precision,dimension(NR_DENSITY),intent(inout) :: rspl,ellipicity_spline,ellipicity_spline2

  ! dummy parameter
  double precision,dimension(NR_DENSITY) :: dummy_eta,dummy_eta2

  ! get ellipticity splines
  call make_ellipticity_epsilon_eta(nspl,rspl,ellipicity_spline,ellipicity_spline2,dummy_eta,dummy_eta2)

  end subroutine make_ellipticity

!
!-----------------------------------------------------------------
!

  subroutine make_ellipticity2(nspl,rspl,ellipicity_spline,ellipicity_spline2,eta_spline,eta_spline2)

! creates a spline for the ellipticity profile in PREM
! radius and density are non-dimensional
! or
! for mars, creates a spline for the ellipticity profile in Mars (Sohl $ Spohn)


  use constants, only: NR_DENSITY

  implicit none

  integer,intent(inout) :: nspl
  double precision,dimension(NR_DENSITY),intent(inout) :: rspl,ellipicity_spline,ellipicity_spline2
  double precision,dimension(NR_DENSITY),intent(inout) :: eta_spline,eta_spline2

  ! get both ellipticity and eta splines
  call make_ellipticity_epsilon_eta(nspl,rspl,ellipicity_spline,ellipicity_spline2,eta_spline,eta_spline2)

  end subroutine make_ellipticity2

!
!-----------------------------------------------------------------
!

  subroutine make_ellipticity_epsilon_eta(nspl,rspl,ellipicity_spline,ellipicity_spline2,eta_spline,eta_spline2)

! creates a spline for the ellipticity profile in PREM
! radius and density are non-dimensional
! or
! for mars, creates a spline for the ellipticity profile in Mars (Sohl $ Spohn)


  use constants, only: NR_DENSITY,TWO_PI,PI,GRAV,R_UNIT_SPHERE,myrank

  use shared_parameters, only: PLANET_TYPE,IPLANET_EARTH,IPLANET_MARS,IPLANET_MOON, &
    R_PLANET,R_PLANET_KM,RHOAV, &
    HOURS_PER_DAY,SECONDS_PER_HOUR

  ! reference models
  use model_prem_par
  use model_sohl_par
  use model_vpremoon_par

  implicit none

  integer,intent(inout) :: nspl
  double precision,dimension(NR_DENSITY),intent(inout) :: rspl,ellipicity_spline,ellipicity_spline2
  double precision,dimension(NR_DENSITY),intent(inout) :: eta_spline,eta_spline2

  ! local parameters
  integer :: i
  ! radii
  double precision :: ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R220,R400,R600,R670, &
                      R771,RTOPDDOUBLEPRIME,RCMB,RICB,RSURFACE
  double precision :: r_icb,r_cmb,r_topddoubleprime,r_771,r_670,r_600
  double precision :: r_400,r_220,r_80,r_moho,r_middle_crust,r_ocean,r_0
  double precision :: SOHL_RMOHO,SOHL_R80,SOHL_R220,SOHL_R400,SOHL_R600,SOHL_R670, &
                      SOHL_R771,SOHL_RTOPDDOUBLEPRIME,SOHL_RCMB

  double precision,dimension(NR_DENSITY) :: r,rho,epsilonval,eta
  double precision,dimension(NR_DENSITY) :: radau,k
  double precision,dimension(NR_DENSITY) :: s1,s2,s3

  double precision :: z,g_a,bom,exponentval,integral_rho,integral_radau
  double precision :: yp1,ypn

  ! for eta splines
  double precision :: doteps,ddoteps,epsinv

  ! debugging
  logical, parameter :: DEBUG = .false.
  double precision :: ell,radius

! note: ellipicity is implemented by Radau's approximation to Clairaut's equation.
!
!       details can be found in: doc/notes/ellipticity_equations_from_Dahlen_Tromp_1998.pdf
!
! for calculating the ellipicity factor epsilon(r) anywhere on Earth, one needs to integrate over Earth's density profile.
! we use PREM's density model for that (including the ocean layer on top).
!
! as a todo in future: this PREM density profile might be slightly off for different Earth models,
!                      please check the effect.

  ! Earth
  ! PREM radius of the Earth for gravity calculation
  double precision, parameter :: R_EARTH_ELLIPTICITY = 6371000.d0
  ! PREM radius of the ocean floor for gravity calculation
  double precision, parameter :: ROCEAN_ELLIPTICITY = PREM_ROCEAN

  ! initializes
  RSURFACE = 0.d0
  ROCEAN = 0.d0
  RMIDDLE_CRUST = 0.d0
  RMOHO = 0.d0
  R80  = 0.d0
  R220 = 0.d0
  R400 = 0.d0
  R600 = 0.d0
  R670 = 0.d0
  R771 = 0.d0
  RTOPDDOUBLEPRIME = 0.d0
  RCMB = 0.d0
  RICB = 0.d0

  ! selects radii
  select case (PLANET_TYPE)
  case (IPLANET_EARTH)
    ! Earth
    ! default PREM
    ! radius of the planet for gravity calculation
    RSURFACE = R_EARTH_ELLIPTICITY  ! physical surface (Earth: 6371000, ..)
    ROCEAN = ROCEAN_ELLIPTICITY
    RMIDDLE_CRUST = PREM_RMIDDLE_CRUST
    RMOHO = PREM_RMOHO
    R80  = PREM_R80
    R220 = PREM_R220
    R400 = PREM_R400
    R600 = PREM_R600
    R670 = PREM_R670
    R771 = PREM_R771
    RTOPDDOUBLEPRIME = PREM_RTOPDDOUBLEPRIME
    RCMB = PREM_RCMB
    RICB = PREM_RICB

  case (IPLANET_MARS)
    ! Mars
    ! default Sohn & Spohn Model
    ! gets corresponding Sohl & Spoon model radii
    call get_model_Sohl_radii(SOHL_RMOHO,SOHL_R80,SOHL_R220,SOHL_R400,SOHL_R600,SOHL_R670, &
                              SOHL_R771,SOHL_RTOPDDOUBLEPRIME,SOHL_RCMB)
    ! sets radii for density integration
    RSURFACE = R_PLANET
    ROCEAN = SOHL_ROCEAN
    RMIDDLE_CRUST = SOHL_RMIDDLE_CRUST
    RMOHO = SOHL_RMOHO
    R80  = SOHL_R80
    R220 = SOHL_R220
    R400 = SOHL_R400
    R600 = SOHL_R600
    R670 = SOHL_R670
    R771 = SOHL_R771
    RTOPDDOUBLEPRIME = SOHL_RTOPDDOUBLEPRIME
    RCMB = SOHL_RCMB
    RICB = SOHL_RICB

  case (IPLANET_MOON)
    ! Moon
    ! default VPREMOON
    RSURFACE = R_PLANET
    ROCEAN = VPREMOON_ROCEAN
    RMIDDLE_CRUST = VPREMOON_RMIDDLE_CRUST
    RMOHO = VPREMOON_RMOHO
    R80  = VPREMOON_R80
    R220 = VPREMOON_R220
    R400 = VPREMOON_R400
    R600 = VPREMOON_R600
    R670 = VPREMOON_R670
    R771 = VPREMOON_R771
    RTOPDDOUBLEPRIME = VPREMOON_RTOPDDOUBLEPRIME
    RCMB = VPREMOON_RCMB
    RICB = VPREMOON_RICB

  case default
    call exit_MPI(myrank,'Invalid planet, ellipticity not implemented yet')
  end select

  ! non-dimensionalize
  r_icb = RICB / RSURFACE
  r_cmb = RCMB / RSURFACE
  r_topddoubleprime = RTOPDDOUBLEPRIME / RSURFACE
  r_771 = R771 / RSURFACE
  r_670 = R670 / RSURFACE
  r_600 = R600 / RSURFACE
  r_400 = R400 / RSURFACE
  r_220 = R220 / RSURFACE
  r_80 = R80 / RSURFACE
  r_moho = RMOHO / RSURFACE
  r_middle_crust = RMIDDLE_CRUST / RSURFACE
  r_ocean = ROCEAN / RSURFACE
  r_0 = 1.d0

  ! sets sampling points in different layers
  ! inner core
  do i = 1,163
    r(i) = r_icb*dble(i-1)/dble(162)
  enddo
  ! outer core
  do i = 164,323
    r(i) = r_icb+(r_cmb-r_icb)*dble(i-164)/dble(159)
  enddo
  ! D''
  do i = 324,336
    r(i) = r_cmb+(r_topddoubleprime-r_cmb)*dble(i-324)/dble(12)
  enddo
  ! D'' to 771
  do i = 337,517
    r(i) = r_topddoubleprime+(r_771-r_topddoubleprime)*dble(i-337)/dble(180)
  enddo
  ! 771 to 670
  do i = 518,530
    r(i) = r_771+(r_670-r_771)*dble(i-518)/dble(12)
  enddo
  ! 670 to 600
  do i = 531,540
    r(i) = r_670+(r_600-r_670)*dble(i-531)/dble(9)
  enddo
  ! 600 to 400
  do i = 541,565
    r(i) = r_600+(r_400-r_600)*dble(i-541)/dble(24)
  enddo
  ! 400 to 220
  do i = 566,590
    r(i) = r_400+(r_220-r_400)*dble(i-566)/dble(24)
  enddo
  ! 220 to 80
  do i = 591,609
    r(i) = r_220+(r_80-r_220)*dble(i-591)/dble(18)
  enddo
  ! 80 to Moho
  do i = 610,619
    r(i) = r_80+(r_moho-r_80)*dble(i-610)/dble(9)
  enddo
  ! Moho to middle crust
  do i = 620,626
    r(i) = r_moho+(r_middle_crust-r_moho)*dble(i-620)/dble(6)
  enddo
  ! middle crust to ocean
  do i = 627,633
    r(i) = r_middle_crust+(r_ocean-r_middle_crust)*dble(i-627)/dble(6)
  enddo
  ! ocean
  do i = 634,NR_DENSITY   ! NR_DENSITY = 640
    r(i) = r_ocean+(r_0-r_ocean)*dble(i-634)/dble(6)
  enddo

  ! density profile
  select case(PLANET_TYPE)
  case (IPLANET_EARTH)
    ! Earth
    ! use PREM to get the density profile for ellipticity (fine for other 1D reference models)
    do i = 1,NR_DENSITY
      call prem_density(r(i),rho(i))
      radau(i) = rho(i)*r(i)*r(i)
    enddo

  case (IPLANET_MARS)
    ! Mars
    ! Sohn & Spohn Model A
    ! No Ocean
    do i = 627,NR_DENSITY ! NR_DENSITY = 640
      r(i) = r_middle_crust+(r_0-r_middle_crust)*dble(i-627)/dble(12)
    enddo
    ! use Sohl & Spohn model (1997) to get the density profile for ellipticity.
    do i = 1,NR_DENSITY
      call Sohl_density(r(i),rho(i))
      radau(i) = rho(i)*r(i)*r(i)
    enddo

  case (IPLANET_MOON)
    ! Moon
    ! default VPREMOON
    ! No Ocean
    do i = 627,NR_DENSITY ! NR_DENSITY = 640
      r(i) = r_middle_crust+(r_0-r_middle_crust)*dble(i-627)/dble(12)
    enddo
    do i = 1,NR_DENSITY
      call model_vpremoon_density(r(i),rho(i))
      radau(i) = rho(i)*r(i)*r(i)
    enddo

  case default
    call exit_MPI(myrank,'Invalid planet, ellipticity not implemented yet')
  end select

  eta(1) = 0.0d0
  k(1) = 0.0d0

  do i = 2,NR_DENSITY
    call intgrl(integral_rho,r,1,i,rho,s1,s2,s3)

! Radau approximation of Clairaut's equation for first-order terms of ellipticity, see e.g. Jeffreys H.,
! The figures of rotating planets, Mon. Not. R. astr. Soc., vol. 113, p. 97-105 (1953).
! The Radau approximation is mentioned on page 97.
! For more details see Section 14.1.2 in Dahlen and Tromp (1998)
! (see also in file ellipticity_equations_from_Dahlen_Tromp_1998.pdf in the "doc" directory of the code).
    call intgrl(integral_radau,r,1,i,radau,s1,s2,s3)

    z = (2.0d0/3.0d0) * integral_radau / (integral_rho*r(i)*r(i))

    ! this comes from equation (14.19) in Dahlen and Tromp (1998)
    eta(i) = (25.0d0/4.0d0)*((1.0d0-(3.0d0/2.0d0)*z)**2.0d0)-1.0d0
    k(i) = eta(i)/(r(i)**3.0d0)
  enddo

  ! day rotation
  bom = TWO_PI/(HOURS_PER_DAY*SECONDS_PER_HOUR)

  ! non-dimensionalized value
  bom = bom/sqrt(PI*GRAV*RHOAV)

  g_a = 4.0d0 * integral_rho
  ! this is the equation right above (14.21) in Dahlen and Tromp (1998)
  epsilonval(NR_DENSITY) = (5.0d0/2.d0)*(bom**2.0d0)*R_UNIT_SPHERE / (g_a * (eta(NR_DENSITY)+2.0d0))

  do i = 1,NR_DENSITY-1
    call intgrl(exponentval,r,i,NR_DENSITY,k,s1,s2,s3)
    epsilonval(i) = epsilonval(NR_DENSITY)*exp(-exponentval)
  enddo

  ! initializes spline coefficients
  rspl(:) = 0.d0
  ellipicity_spline(:) = 0.d0
  ellipicity_spline2(:) = 0.d0

  ! get ready to spline epsilonval
  nspl = 1
  rspl(1) = r(1)
  ellipicity_spline(1) = epsilonval(1)
  do i = 2,NR_DENSITY
    if (r(i) /= r(i-1)) then
      nspl = nspl+1
      rspl(nspl) = r(i)
      ellipicity_spline(nspl) = epsilonval(i)
    endif
  enddo

  ! spline epsilonval
  yp1 = 0.0d0
  ypn = (5.0d0/2.0d0)*(bom**2)/g_a-2.0d0*epsilonval(NR_DENSITY)

  call spline_construction(rspl,ellipicity_spline,nspl,yp1,ypn,ellipicity_spline2)

  ! debug
  if (DEBUG) then
    if (myrank == 0) then
      print *,'debug: make ellipticity'
      print *,'debug: number of splines = ',nspl
      print *,'#r(km) #ellipticity_factor #radau #k #r(non-dim) #i'
      do i = 1,NR_DENSITY
        radius = r(i)
        ! get ellipticity using spline evaluation
        call spline_evaluation(rspl,ellipicity_spline,ellipicity_spline2,nspl,radius,ell)
        print *,radius*R_PLANET_KM,ell,radau(i),k(i),radius,i
      enddo
      print *
    endif
  endif

  ! spline evaluation for eta values
  rspl(:) = 0.d0
  eta_spline(:) = 0.d0
  eta_spline2(:) = 0.d0

  ! get ready to spline eta
  nspl = 1
  rspl(1) = r(1)
  eta_spline(1) = eta(1)
  do i = 2,NR_DENSITY
    if (r(i) /= r(i-1)) then
      nspl = nspl+1
      rspl(nspl) = r(i)
      eta_spline(nspl) = eta(i)
    endif
  enddo

  ! spline eta values
  epsinv = 1.0d0/epsilonval(NR_DENSITY)
  doteps = eta(NR_DENSITY)*epsilonval(NR_DENSITY)
  ddoteps = 6.0d0*epsilonval(NR_DENSITY) - 2.0d0*4.0d0*rho(NR_DENSITY)*(doteps+epsilonval(NR_DENSITY))/g_a

  yp1 = 0.0d0
  ypn = doteps*epsinv-0.5d0*doteps*doteps*epsinv*epsinv+ddoteps*epsinv

  call spline_construction(rspl,eta_spline,nspl,yp1,ypn,eta_spline2)

  end subroutine make_ellipticity_epsilon_eta

!
!-----------------------------------------------------------------
!

  subroutine add_ellipticity(x,y,z,nspl,rspl,ellipicity_spline,ellipicity_spline2)

! adds ellipticity factor to position x/y/z

  use constants

  implicit none

  double precision, intent(inout) :: x,y,z
  integer,intent(in) :: nspl
  double precision,intent(in) :: rspl(NR_DENSITY),ellipicity_spline(NR_DENSITY),ellipicity_spline2(NR_DENSITY)

  ! local parameters
  double precision :: ell
  double precision :: r,theta,phi,factor
  double precision :: cost,p20

  ! gets spherical coordinates
  call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)

  cost = dcos(theta)

  ! this is the Legendre polynomial of degree two, P2(cos(theta)),
  ! see the discussion above eq (14.4) in Dahlen and Tromp (1998)
  p20 = 0.5d0*(3.0d0*cost*cost-1.0d0)

  ! get ellipticity using spline evaluation
  call spline_evaluation(rspl,ellipicity_spline,ellipicity_spline2,nspl,r,ell)

  ! this is eq (14.4) in Dahlen and Tromp (1998)
  factor = ONE - (TWO/3.0d0)*ell*p20

  ! applies ellipticity
  x = x * factor
  y = y * factor
  z = z * factor

  end subroutine add_ellipticity

!
!-----------------------------------------------------------------
!

  subroutine add_ellipticity_rtheta(r,theta,nspl,rspl,ellipicity_spline,ellipicity_spline2)

! adds ellipticity factor to radius r

  use constants

  implicit none

  double precision, intent(inout) :: r
  double precision, intent(in) :: theta
  integer,intent(in) :: nspl
  double precision,intent(in) :: rspl(NR_DENSITY),ellipicity_spline(NR_DENSITY),ellipicity_spline2(NR_DENSITY)

  ! local parameters
  double precision :: ell
  double precision :: factor
  double precision :: cost,p20

  ! pre-calculates coefficient
  cost = dcos(theta)

  ! this is the Legendre polynomial of degree two, P2(cos(theta)),
  ! see the discussion above eq (14.4) in Dahlen and Tromp (1998)
  p20 = 0.5d0*(3.0d0*cost*cost-1.0d0)

  ! get ellipticity using spline evaluation
  call spline_evaluation(rspl,ellipicity_spline,ellipicity_spline2,nspl,r,ell)

  ! this is eq (14.4) in Dahlen and Tromp (1998)
  factor = ONE - (TWO/3.0d0)*ell*p20

  ! applies ellipticity
  r = r * factor

  end subroutine add_ellipticity_rtheta

!
!-----------------------------------------------------------------
!

  subroutine revert_ellipticity_cr(x,y,z,nspl,rspl,ellipicity_spline,ellipicity_spline2)

! routine to revert ellipticity and go back to a spherical Earth

  use constants

  implicit none

  real(kind=CUSTOM_REAL),intent(inout) :: x,y,z
  integer,intent(in) :: nspl
  double precision,intent(in) :: rspl(NR_DENSITY),ellipicity_spline(NR_DENSITY),ellipicity_spline2(NR_DENSITY)

  ! local parameters
  double precision :: x1,y1,z1

  ! converts to double
  x1 = x
  y1 = y
  z1 = z

  ! revert ellipticity
  call revert_ellipticity(x1,y1,z1,nspl,rspl,ellipicity_spline,ellipicity_spline2)

  ! converts to custom real
  x = real(x1,kind=CUSTOM_REAL)
  y = real(y1,kind=CUSTOM_REAL)
  z = real(z1,kind=CUSTOM_REAL)

  end subroutine revert_ellipticity_cr

!
!-----------------------------------------------------------------
!

  subroutine revert_ellipticity(x,y,z,nspl,rspl,ellipicity_spline,ellipicity_spline2)

! routine to revert ellipticity and go back to a spherical Earth

! (in double precision)

! note: going from an elliptical reference point x/y/z back to a spherical reference point x_s/y_s/z_s
!       is done here only in an approximate way as we use the stretch factor ell evaluated from x/y/z.
!
!       for an exact projection back to the spherical position, we would need to get the factor ell evaluated
!       from the the original position x_s/y_s/z_s and radius r_s, and not from r as done below.
!
!       one could try to do this iteratively, i.e., once we got the first spherical position, calculate again
!       the elliptical position and compare the input with the newly evaluated, take the difference, correct the
!       the stretch factor accordingly and project back to an updated spherical position - and so on.
!
!       this will take more computational time, so we return only the first approximate spherical position here...

  use constants

  implicit none

  double precision,intent(inout) :: x,y,z
  integer,intent(in) :: nspl
  double precision,intent(in) :: rspl(NR_DENSITY),ellipicity_spline(NR_DENSITY),ellipicity_spline2(NR_DENSITY)

  ! local parameters
  double precision :: ell
  double precision :: r,theta,phi,factor
  double precision :: cost,p20

  ! gets spherical coordinates
  call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)

  cost = dcos(theta)

  ! this is the Legendre polynomial of degree two, P2(cos(theta)),
  ! see the discussion above eq (14.4) in Dahlen and Tromp (1998)
  p20 = 0.5d0*(3.0d0*cost*cost-1.0d0)

  ! get ellipticity using spline evaluation
  call spline_evaluation(rspl,ellipicity_spline,ellipicity_spline2,nspl,r,ell)

  ! this is eq (14.4) in Dahlen and Tromp (1998)
  factor = ONE - (TWO/3.0d0)*ell*p20

  ! removes ellipticity factor from x/y/z position
  ! (assuming that x/y/z position has been stretched before to account for ellipticity,
  !  this puts it back, approximately, to a spherical position)
  x = x / factor
  y = y / factor
  z = z / factor

  end subroutine revert_ellipticity

!
!-----------------------------------------------------------------
!

  subroutine revert_ellipticity_rtheta(r,theta,nspl,rspl,ellipicity_spline,ellipicity_spline2)

! routine to revert ellipticity and go back to a spherical Earth

  use constants

  implicit none

  double precision,intent(inout) :: r
  double precision,intent(in) :: theta

  integer,intent(in) :: nspl
  double precision,intent(in) :: rspl(NR_DENSITY),ellipicity_spline(NR_DENSITY),ellipicity_spline2(NR_DENSITY)

  ! local parameters
  double precision :: ell
  double precision :: factor
  double precision :: cost,p20

  ! P20
  cost = dcos(theta)
  ! this is the Legendre polynomial of degree two, P2(cos(theta)), see the discussion above eq (14.4) in Dahlen and Tromp (1998)
  p20 = 0.5d0*(3.0d0*cost*cost-1.0d0)

  ! get ellipticity using spline evaluation
  call spline_evaluation(rspl,ellipicity_spline,ellipicity_spline2,nspl,r,ell)

  ! this is eq (14.4) in Dahlen and Tromp (1998)
  factor = ONE - (TWO/3.0d0)*ell*p20

  ! removes ellipticity factor
  r = r / factor

  end subroutine revert_ellipticity_rtheta


