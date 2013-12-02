!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            August 2013
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

  integer myrank

  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)

  integer ia

  real(kind=4) xcolat,xlon
  real(kind=4) topo410out,topo650out
  double precision topo410,topo650

  double precision :: r,lat,lon
  double precision :: gamma
  double precision :: x,y,z

!! DK DK added this safety test for now
  stop 'there is a currently a bug in this routine, it makes the mesher crash'

! we loop on all the points of the element
  do ia = 1,NGNOD

    x = xelm(ia)
    y = xelm(ia)
    z = xelm(ia)

    ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
    call xyz_2_rlatlon_dble(x,y,z,r,lat,lon)

    ! converts lat/lon to (real) colat/lon in degrees
    xcolat = sngl(90.d0 - lat)     ! colatitude
    xlon = sngl(lon)

    ! compute topography on 410 and 650 at current point
    call model_s362ani_subtopo(xcolat,xlon,topo410out,topo650out)

    ! non-dimensionalize the topography, which is in km
    ! positive for a depression, so change the sign for a perturbation in radius
    topo410 = -dble(topo410out) / R_EARTH_KM
    topo650 = -dble(topo650out) / R_EARTH_KM

    gamma = 0.d0
    if(r >= R400/R_EARTH .and. r <= R220/R_EARTH) then
      ! stretching between R220 and R400
      gamma = (R220/R_EARTH - r) / (R220/R_EARTH - R400/R_EARTH)
      xelm(ia) = x*(ONE + gamma * topo410 / r)
      yelm(ia) = y*(ONE + gamma * topo410 / r)
      zelm(ia) = z*(ONE + gamma * topo410 / r)
    else if(r>= R771/R_EARTH .and. r <= R670/R_EARTH) then
      ! stretching between R771 and R670
      gamma = (r - R771/R_EARTH) / (R670/R_EARTH - R771/R_EARTH)
      xelm(ia) = x*(ONE + gamma * topo650 / r)
      yelm(ia) = y*(ONE + gamma * topo650 / r)
      zelm(ia) = z*(ONE + gamma * topo650 / r)
    else if(r > R670/R_EARTH .and. r < R400/R_EARTH) then
      ! stretching between R670 and R400
      gamma = (R400/R_EARTH - r) / (R400/R_EARTH - R670/R_EARTH)
      xelm(ia) = x*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
      yelm(ia) = y*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
      zelm(ia) = z*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
    endif
    if(gamma < -0.0001 .or. gamma > 1.0001) call exit_MPI(myrank,'incorrect value of gamma for 410-650 topography')

  enddo

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

  double precision :: r,lat,lon
  double precision :: gamma
  double precision :: x,y,z

!! DK DK added this safety test for now
  stop 'there is a currently a bug in this routine, it makes the mesher crash'

  ! we loop on all GLL points of the element
  do k = 1,NGLLZ
     do j = 1,NGLLY
        do i = 1,NGLLX

          x = xstore(i,j,k,ispec)
          y = ystore(i,j,k,ispec)
          z = zstore(i,j,k,ispec)

          ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
          call xyz_2_rlatlon_dble(x,y,z,r,lat,lon)

          ! converts lat/lon to (real) colat/lon in degrees
          xcolat = sngl(90.d0 - lat)     ! colatitude
          xlon = sngl(lon)

          ! compute topography on 410 and 650 at current point
          call model_s362ani_subtopo(xcolat,xlon,topo410out,topo650out)

          ! non-dimensionalize the topography, which is in km
          ! positive for a depression, so change the sign for a perturbation in radius
          topo410 = -dble(topo410out) / R_EARTH_KM
          topo650 = -dble(topo650out) / R_EARTH_KM

          gamma = 0.d0
          if(r >= R400/R_EARTH .and. r <= R220/R_EARTH) then
          ! stretching between R220 and R400
                  gamma = (R220/R_EARTH - r) / (R220/R_EARTH - R400/R_EARTH)
                  xstore(i,j,k,ispec) = x*(ONE + gamma * topo410 / r)
                  ystore(i,j,k,ispec) = y*(ONE + gamma * topo410 / r)
                  zstore(i,j,k,ispec) = z*(ONE + gamma * topo410 / r)
          else if(r>= R771/R_EARTH .and. r <= R670/R_EARTH) then
          ! stretching between R771 and R670
                  gamma = (r - R771/R_EARTH) / (R670/R_EARTH - R771/R_EARTH)
                  xstore(i,j,k,ispec) = x*(ONE + gamma * topo650 / r)
                  ystore(i,j,k,ispec) = y*(ONE + gamma * topo650 / r)
                  zstore(i,j,k,ispec) = z*(ONE + gamma * topo650 / r)
          else if(r > R670/R_EARTH .and. r < R400/R_EARTH) then
          ! stretching between R670 and R400
                  gamma = (R400/R_EARTH - r) / (R400/R_EARTH - R670/R_EARTH)
                  xstore(i,j,k,ispec) = x*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
                  ystore(i,j,k,ispec) = y*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
                  zstore(i,j,k,ispec) = z*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
          endif
          if(gamma < -0.0001 .or. gamma > 1.0001) call exit_MPI(myrank,'incorrect value of gamma for 410-650 topography')

        enddo
     enddo
  enddo

  end subroutine add_topography_410_650_gll
