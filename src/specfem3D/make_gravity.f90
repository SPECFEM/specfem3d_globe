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

  subroutine make_gravity(nspl,rspl,gravity_spline,gravity_spline2)

! creates a spline for the gravity profile in PREM
! radius and density are non-dimensional

  use constants, only: NR_DENSITY,myrank
  use shared_parameters, only: PLANET_TYPE,IPLANET_EARTH,IPLANET_MARS,IPLANET_MOON,R_PLANET,R_PLANET_KM

  ! reference models
  use model_prem_par
  use model_sohl_par
  use model_vpremoon_par

  implicit none

  integer,intent(out) :: nspl

  double precision,intent(inout) :: rspl(NR_DENSITY),gravity_spline(NR_DENSITY),gravity_spline2(NR_DENSITY)

  ! local parameters
  integer :: i
  double precision :: ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R220,R400,R600,R670, &
                      R771,RTOPDDOUBLEPRIME,RCMB,RICB,RSURFACE
  double precision :: r_icb,r_cmb,r_topddoubleprime,r_771,r_670,r_600
  double precision :: r_400,r_220,r_80,r_moho,r_middle_crust,r_ocean,r_0
  double precision :: r(NR_DENSITY),rho(NR_DENSITY),grav(NR_DENSITY),integral_rho
  double precision :: s1(NR_DENSITY),s2(NR_DENSITY),s3(NR_DENSITY)
  double precision :: yp1,ypn
  double precision :: SOHL_RMOHO,SOHL_R80,SOHL_R220,SOHL_R400,SOHL_R600,SOHL_R670, &
                      SOHL_R771,SOHL_RTOPDDOUBLEPRIME,SOHL_RCMB

  ! Earth
  ! PREM radius of the Earth for gravity calculation
  double precision, parameter :: R_EARTH_GRAVITY = 6371000.d0
  ! PREM radius of the ocean floor for gravity calculation
  double precision, parameter :: ROCEAN_GRAVITY = PREM_ROCEAN

  ! debugging
  logical, parameter :: DEBUG = .false.
  double precision :: gval,radius

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
    RSURFACE = R_EARTH_GRAVITY  ! physical surface (Earth: 6371000, ..)
    ROCEAN = ROCEAN_GRAVITY
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
    ! gets corresponding Sohl&Spoon model radii
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

! note: for Mars
!       discretize the radius point to fit the spline, inherited from PREM
! To do: (may not be necessary)
!        Find number of points at each layer to have better integrations

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
  do i = 634,NR_DENSITY ! NR_DENSITY = 640
    r(i) = r_ocean+(r_0-r_ocean)*dble(i-634)/dble(6)
  enddo

  ! density profile
  select case(PLANET_TYPE)
  case (IPLANET_EARTH)
    ! Earth
    ! use PREM to get the density profile for ellipticity (fine for other 1D reference models)
    do i = 1,NR_DENSITY
      call prem_density(r(i),rho(i))
    enddo

  case (IPLANET_MARS)
    ! Mars
    ! Sohn & Spohn Model A
    ! No Ocean
    do i = 627,NR_DENSITY
      r(i) = r_middle_crust+(r_0-r_middle_crust)*dble(i-627)/dble(12)
    enddo
    ! use Sohl & Spohn model to get the density profile for gravity
    do i = 1,NR_DENSITY
      call Sohl_density(r(i),rho(i))
    enddo

  case (IPLANET_MOON)
    ! Moon
    ! No Ocean
    do i = 627,NR_DENSITY
      r(i) = r_middle_crust+(r_0-r_middle_crust)*dble(i-627)/dble(12)
    enddo
    ! use VPREMOON model to get the density profile for gravity
    do i = 1,NR_DENSITY
      call model_vpremoon_density(r(i),rho(i))
    enddo

  case default
    call exit_MPI(myrank,'Invalid planet, gravity not implemented yet')
  end select

  grav(1) = 0.0d0
  do i = 2,NR_DENSITY
    call intgrl(integral_rho,r,1,i,rho,s1,s2,s3)
    grav(i) = 4.0d0 * integral_rho / (r(i)*r(i))
  enddo

!
! get ready to spline g
!
  nspl = 1
  rspl(1) = r(1)
  gravity_spline(1) = grav(1)
  do i = 2,NR_DENSITY
    if (r(i) /= r(i-1)) then
      nspl = nspl+1
      rspl(nspl) = r(i)
      gravity_spline(nspl) = grav(i)
    endif
  enddo

  yp1 = (4.0d0/3.0d0) * rho(1)
  ypn = 4.0d0*rho(NR_DENSITY) - 2.0d0*grav(NR_DENSITY)/r(NR_DENSITY)

  call spline_construction(rspl,gravity_spline,nspl,yp1,ypn,gravity_spline2)

  ! debug
  if (DEBUG) then
    if (myrank == 0) then
      print *,'debug: make gravity'
      print *,'debug: gravity number of splines = ',nspl
      print *,'debug: gravity yp1,ypn = ',yp1,ypn
      print *,'#r(km) #g #grav(table) #r(non-dim) #i'
      do i = 1,NR_DENSITY
        radius = r(i)
        ! get ellipticity using spline evaluation
        call spline_evaluation(rspl,gravity_spline,gravity_spline2,nspl,radius,gval)
        print *,radius*R_PLANET_KM,gval,grav(i),radius,i
      enddo
      print *
    endif
  endif

  end subroutine make_gravity

