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
! the Free Software Foundation; either version 2 of the License, or
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

  subroutine add_topography_410_650(myrank,xelm,yelm,zelm)

  use constants
  use meshfem3D_par,only: R220,R400,R670,R771

  implicit none

  integer :: myrank

  double precision :: xelm(NGNOD)
  double precision :: yelm(NGNOD)
  double precision :: zelm(NGNOD)

  integer :: ia

  real(kind=4) :: xcolat,xlon
  real(kind=4) :: topo410out,topo650out
  double precision :: topo410,topo650

  double precision :: r,lat,lon,theta,phi
  double precision :: gamma
  double precision :: x,y,z

  !statistics
  logical,parameter :: DEBUG_STATISTICS = .false.
  real(kind=CUSTOM_REAL),save :: min_410 = HUGEVAL,max_410 = - HUGEVAL
  real(kind=CUSTOM_REAL),save :: min_650 = HUGEVAL,max_650 = - HUGEVAL
  real(kind=CUSTOM_REAL) :: min_410_all,max_410_all
  real(kind=CUSTOM_REAL) :: min_650_all,max_650_all

! note: adding topography to 410 and 660 strongly affects PcP, PKiKP, etc. phases,
!       we leave it in and check whether the stretching makes simulation unstable
!
! topography perturbations
!   410-km: minimum / maximum = -13.48 km / + 13.24 km
!   650-km: minimum / maximum = -14.34 km / + 19.19 km

! we loop on all the points of the element
  do ia = 1,NGNOD

    x = xelm(ia)
    y = yelm(ia)
    z = zelm(ia)

    if (USE_OLD_VERSION_5_1_5_FORMAT) then
      ! convert to r theta phi
      call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)
      call reduce(theta,phi)
      ! get colatitude and longitude in degrees
      xcolat = sngl(theta*RADIANS_TO_DEGREES)
      xlon = sngl(phi*RADIANS_TO_DEGREES)
    else
      ! note: the topography on 410 and 650 is given in geographic colat/lon,
      !       thus we need to convert geocentric colatitude to geographic colatitudes
      !
      ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
      call xyz_2_rlatlon_dble(x,y,z,r,lat,lon)

      ! converts lat/lon to (real) colat/lon in degrees
      xcolat = sngl(90.d0 - lat)     ! colatitude
      xlon = sngl(lon)
    endif

    ! stretching occurs between 220 and 770
    if (r > R220/R_EARTH .or. r < R771/R_EARTH) cycle

    ! compute topography on 410 and 650 at current point
    call model_s362ani_subtopo(xcolat,xlon,topo410out,topo650out)

    ! min/max statistics
    if (DEBUG_STATISTICS) then
      if ( topo410out < min_410 ) min_410 = topo410out
      if ( topo410out > max_410 ) max_410 = topo410out
      if ( topo650out < min_650 ) min_650 = topo650out
      if ( topo650out > max_650 ) max_650 = topo650out
      ! debug
      !print*,'topo410 / topo650: ',r,xcolat,xlon,topo410out,topo650out
    endif

    ! non-dimensionalize the topography, which is in km
    ! positive for a depression, so change the sign for a perturbation in radius
    topo410 = -dble(topo410out) / R_EARTH_KM
    topo650 = -dble(topo650out) / R_EARTH_KM

    gamma = 0.d0
    if (r >= R400/R_EARTH .and. r <= R220/R_EARTH) then
      ! stretching between R220 and R400
      gamma = (R220/R_EARTH - r) / (R220/R_EARTH - R400/R_EARTH)
      xelm(ia) = x*(ONE + gamma * topo410 / r)
      yelm(ia) = y*(ONE + gamma * topo410 / r)
      zelm(ia) = z*(ONE + gamma * topo410 / r)
    else if (r>= R771/R_EARTH .and. r <= R670/R_EARTH) then
      ! stretching between R771 and R670
      gamma = (r - R771/R_EARTH) / (R670/R_EARTH - R771/R_EARTH)
      xelm(ia) = x*(ONE + gamma * topo650 / r)
      yelm(ia) = y*(ONE + gamma * topo650 / r)
      zelm(ia) = z*(ONE + gamma * topo650 / r)
    else if (r > R670/R_EARTH .and. r < R400/R_EARTH) then
      ! stretching between R670 and R400
      gamma = (R400/R_EARTH - r) / (R400/R_EARTH - R670/R_EARTH)
      xelm(ia) = x*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
      yelm(ia) = y*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
      zelm(ia) = z*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
    endif
    if (gamma < -0.0001 .or. gamma > 1.0001) call exit_MPI(myrank,'incorrect value of gamma for 410-650 topography')

  enddo

  ! debug
  if (DEBUG_STATISTICS) then
    ! collects min/max on master
    call min_all_cr(min_410,min_410_all)
    call max_all_cr(max_410,max_410_all)
    call min_all_cr(min_650,min_650_all)
    call max_all_cr(max_650,max_650_all)
    if (myrank == 0) then
      if (r <= R220/R_EARTH .and. r >= R771/R_EARTH) then
        print*,'add_topography_410_650: min/max_410 = ',min_410_all,max_410_all,'min/max_650 = ',min_650_all,max_650_all
      endif
    endif
    !if (r <= R220/R_EARTH .and. r >= R771/R_EARTH) then
    !  print*,myrank,'add_topography_410_650: min/max_410 = ',min_410,max_410,'min/max_650 = ',min_650,max_650
    !  print*,myrank,'add_topography_410_650: depth = ',(1.d0 - r)*R_EARTH_KM,' 410-km = ',topo410out,' 650-km = ',topo650out
    !endif
  endif

  end subroutine add_topography_410_650

!
!-------------------------------------------------------------------------------------------------
!

  !> Hejun
  ! use GLL points to capture 410_650 topography
  ! JAN08, 2010
  subroutine add_topography_410_650_gll(myrank,xstore,ystore,zstore,ispec,nspec)

  use constants
  use meshfem3D_par,only: R220,R400,R670,R771

  implicit none

  integer myrank
  integer:: ispec,nspec
  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec):: xstore,ystore,zstore

  integer i,j,k

  real(kind=4) xcolat,xlon
  real(kind=4) topo410out,topo650out
  double precision topo410,topo650

  double precision :: r,lat,lon,theta,phi
  double precision :: gamma
  double precision :: x,y,z

! note: adding topography to 410 and 660 strongly affects PcP, PKiKP, etc. phases,
!       we leave it in and check whether the stretching makes simulation unstable

  ! we loop on all GLL points of the element
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        x = xstore(i,j,k,ispec)
        y = ystore(i,j,k,ispec)
        z = zstore(i,j,k,ispec)

        if (USE_OLD_VERSION_5_1_5_FORMAT) then
          ! convert to r theta phi
          call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)
          call reduce(theta,phi)
          ! get colatitude and longitude in degrees
          xcolat = sngl(theta*RADIANS_TO_DEGREES)
          xlon = sngl(phi*RADIANS_TO_DEGREES)
        else
          ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
          call xyz_2_rlatlon_dble(x,y,z,r,lat,lon)

          ! converts lat/lon to (real) colat/lon in degrees
          xcolat = sngl(90.d0 - lat)     ! colatitude
          xlon = sngl(lon)
        endif

        ! stretching occurs between 220 and 770
        if (r > R220/R_EARTH .or. r < R771/R_EARTH) cycle

        ! compute topography on 410 and 650 at current point
        call model_s362ani_subtopo(xcolat,xlon,topo410out,topo650out)

        ! non-dimensionalize the topography, which is in km
        ! positive for a depression, so change the sign for a perturbation in radius
        topo410 = -dble(topo410out) / R_EARTH_KM
        topo650 = -dble(topo650out) / R_EARTH_KM

        gamma = 0.d0
        if (r >= R400/R_EARTH .and. r <= R220/R_EARTH) then
        ! stretching between R220 and R400
                gamma = (R220/R_EARTH - r) / (R220/R_EARTH - R400/R_EARTH)
                xstore(i,j,k,ispec) = x*(ONE + gamma * topo410 / r)
                ystore(i,j,k,ispec) = y*(ONE + gamma * topo410 / r)
                zstore(i,j,k,ispec) = z*(ONE + gamma * topo410 / r)
        else if (r>= R771/R_EARTH .and. r <= R670/R_EARTH) then
        ! stretching between R771 and R670
                gamma = (r - R771/R_EARTH) / (R670/R_EARTH - R771/R_EARTH)
                xstore(i,j,k,ispec) = x*(ONE + gamma * topo650 / r)
                ystore(i,j,k,ispec) = y*(ONE + gamma * topo650 / r)
                zstore(i,j,k,ispec) = z*(ONE + gamma * topo650 / r)
        else if (r > R670/R_EARTH .and. r < R400/R_EARTH) then
        ! stretching between R670 and R400
                gamma = (R400/R_EARTH - r) / (R400/R_EARTH - R670/R_EARTH)
                xstore(i,j,k,ispec) = x*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
                ystore(i,j,k,ispec) = y*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
                zstore(i,j,k,ispec) = z*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
        endif
        if (gamma < -0.0001 .or. gamma > 1.0001) call exit_MPI(myrank,'incorrect value of gamma for 410-650 topography')

      enddo
    enddo
  enddo

  end subroutine add_topography_410_650_gll
