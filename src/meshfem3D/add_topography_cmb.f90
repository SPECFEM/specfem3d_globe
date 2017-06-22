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

  subroutine add_topography_cmb(xelm,yelm,zelm)

! this is only a placeholder function, which is not used yet...user must supply the subtopo_cmb() routine

  use constants
  use meshfem3D_par, only: RTOPDDOUBLEPRIME,RCMB

  implicit none

  double precision,intent(inout) :: xelm(NGNOD)
  double precision,intent(inout) :: yelm(NGNOD)
  double precision,intent(inout) :: zelm(NGNOD)

  ! local parameters
  integer :: ia

  double precision :: r_start,topocmb

  double precision :: r,lat,lon
  double precision :: x,y,z
  double precision :: gamma

! we loop on all the points of the element
  do ia = 1,NGNOD

    x = xelm(ia)
    y = yelm(ia)
    z = zelm(ia)

    ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
    call xyz_2_rlatlon_dble(x,y,z,r,lat,lon)

    ! compute topography on CMB; routine subtopo_cmb needs to be supplied by the user
    ! (see for example routine subtopo_sh_cmb() in model_full_sh.f90)
    !   call subtopo_cmb(lat,lon,topocmb)
    ! until then, no topography pertubations...
    topocmb = 0.0d0

    ! non-dimensionalize the topography, which is in km
    ! positive for a depression, so change the sign for a perturbation in radius
    topocmb = -topocmb / R_EARTH_KM

    ! start stretching a distance RTOPDDOUBLEPRIME - RCMB below the CMB
    ! and finish at RTOPDDOUBLEPRIME of D_double_prime
    r_start = (RCMB - (RTOPDDOUBLEPRIME - RCMB))/R_EARTH
    gamma = 0.0d0
    if (r >= RCMB/R_EARTH .and. r <= RTOPDDOUBLEPRIME/R_EARTH) then
      ! stretching between RCMB and RTOPDDOUBLEPRIME
      gamma = (RTOPDDOUBLEPRIME/R_EARTH - r) / (RTOPDDOUBLEPRIME/R_EARTH - RCMB/R_EARTH)
    else if (r >= r_start .and. r <= RCMB/R_EARTH) then
      ! stretching between r_start and RCMB
      gamma = (r - r_start) / (RCMB/R_EARTH - r_start)
    endif
    if (gamma < -0.0001 .or. gamma > 1.0001) call exit_MPI(myrank,'incorrect value of gamma for CMB topography')

    xelm(ia) = x*(ONE + gamma * topocmb / r)
    yelm(ia) = y*(ONE + gamma * topocmb / r)
    zelm(ia) = z*(ONE + gamma * topocmb / r)

  enddo

  end subroutine add_topography_cmb

