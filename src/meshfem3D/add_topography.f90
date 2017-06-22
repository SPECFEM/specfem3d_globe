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

  subroutine add_topography(xelm,yelm,zelm,ibathy_topo)

  use constants, only: myrank, &
    NGNOD,NX_BATHY,NY_BATHY,R_EARTH,R_UNIT_SPHERE, &
    PI_OVER_TWO,RADIANS_TO_DEGREES,TINYVAL,ONE

  use meshfem3D_par, only: R220

  implicit none

  double precision,dimension(NGNOD) :: xelm,yelm,zelm

  ! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

  ! local parameters
  double precision :: r,lat,lon,elevation
  double precision :: x,y,z
  double precision :: gamma

  integer :: ia

  ! we loop on all the points of the element
  do ia = 1,NGNOD

    x = xelm(ia)
    y = yelm(ia)
    z = zelm(ia)

    ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
    call xyz_2_rlatlon_dble(x,y,z,r,lat,lon)

    ! compute elevation at current point
    call get_topo_bathy(lat,lon,elevation,ibathy_topo)

    ! non-dimensionalize the elevation, which is in meters
    elevation = elevation / R_EARTH

    ! stretching topography between d220 and the surface
    gamma = (r - R220/R_EARTH) / (R_UNIT_SPHERE - R220/R_EARTH)

    ! add elevation to all the points of that element
    ! also make sure gamma makes sense
    if (gamma < -0.02 .or. gamma > 1.02) call exit_MPI(myrank,'incorrect value of gamma for topography')

    xelm(ia) = x*(ONE + gamma * elevation / r)
    yelm(ia) = y*(ONE + gamma * elevation / r)
    zelm(ia) = z*(ONE + gamma * elevation / r)

  enddo

  end subroutine add_topography

!
!-------------------------------------------------------------------------------------------------
!

  !> Hejun
  ! This subroutine uses GLL points to capture topography variation rather
  ! than using control nodes
  ! Hejun Zhu, OCT16, 2009

  subroutine add_topography_gll(xstore,ystore,zstore,ispec,nspec, &
                                ibathy_topo)

  use constants
  use meshfem3D_par, only: R220

  implicit none

  ! input parameters
  integer:: ispec,nspec

  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec):: xstore,ystore,zstore

  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

  ! local parameters used in this subroutine
  integer :: i,j,k
  double precision :: r,lat,lon,elevation,gamma
  double precision :: x,y,z

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        x = xstore(i,j,k,ispec)
        y = ystore(i,j,k,ispec)
        z = zstore(i,j,k,ispec)

        ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
        call xyz_2_rlatlon_dble(x,y,z,r,lat,lon)

        ! compute elevation at current point
        call get_topo_bathy(lat,lon,elevation,ibathy_topo)

        ! non-dimensionalize the elevation, which is in meters
        elevation = elevation / R_EARTH

        ! stretching topography between d220 and the surface
        gamma = (r - R220/R_EARTH) / (R_UNIT_SPHERE - R220/R_EARTH)

        ! add elevation to all the points of that element
        ! also make sure factor makes sense
        if (gamma < -0.02 .or. gamma > 1.02) then
          call exit_MPI(myrank,'incorrect value of factor for topography GLL points')
        endif

        ! since not all GLL points are exactly at R220, use a small
        ! tolerance for R220 detection
        if (abs(gamma) < SMALLVAL) then
          gamma = 0.d0
        endif

        xstore(i,j,k,ispec) = x*(ONE + gamma * elevation / r)
        ystore(i,j,k,ispec) = y*(ONE + gamma * elevation / r)
        zstore(i,j,k,ispec) = z*(ONE + gamma * elevation / r)
      enddo
    enddo
  enddo

  end subroutine add_topography_gll
