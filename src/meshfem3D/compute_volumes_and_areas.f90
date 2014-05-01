!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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

  subroutine compute_volumes_and_areas(myrank,NCHUNKS,iregion_code,nspec,wxgll,wygll,wzgll,xixstore,xiystore,xizstore, &
                            etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
                            NSPEC2D_BOTTOM,jacobian2D_bottom,NSPEC2D_TOP,jacobian2D_top,idoubling, &
                            volume_total,RCMB,RICB,R_CENTRAL_CUBE)

  use constants

  use meshfem3D_models_par

  implicit none

  integer :: nspec
  double precision :: wxgll(NGLLX),wygll(NGLLY),wzgll(NGLLZ)

  integer :: myrank,NCHUNKS,iregion_code

  double precision :: volume_total
  double precision :: RCMB,RICB,R_CENTRAL_CUBE

  integer,dimension(nspec) :: idoubling

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

  integer :: NSPEC2D_BOTTOM,NSPEC2D_TOP
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM) :: jacobian2D_bottom
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP) :: jacobian2D_top

  ! local parameters
  double precision :: volume_local,area_local_bottom,area_local_top
  double precision :: volume_total_region,area_total_bottom,area_total_top
  double precision :: weight
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  integer :: i,j,k,ispec

  ! initializes
  volume_local = ZERO
  area_local_bottom = ZERO
  area_local_top = ZERO

  ! calculates volume of all elements in mesh
  do ispec = 1,nspec

    ! suppress fictitious elements in central cube
    if(idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          weight = wxgll(i)*wygll(j)*wzgll(k)

          ! compute the Jacobian
          xixl = xixstore(i,j,k,ispec)
          xiyl = xiystore(i,j,k,ispec)
          xizl = xizstore(i,j,k,ispec)
          etaxl = etaxstore(i,j,k,ispec)
          etayl = etaystore(i,j,k,ispec)
          etazl = etazstore(i,j,k,ispec)
          gammaxl = gammaxstore(i,j,k,ispec)
          gammayl = gammaystore(i,j,k,ispec)
          gammazl = gammazstore(i,j,k,ispec)

          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))

          volume_local = volume_local + dble(jacobianl)*weight

        enddo
      enddo
    enddo
  enddo

  ! area of bottom surface
  do ispec = 1,NSPEC2D_BOTTOM
    do i=1,NGLLX
      do j=1,NGLLY
        weight=wxgll(i)*wygll(j)
        area_local_bottom = area_local_bottom + dble(jacobian2D_bottom(i,j,ispec))*weight
      enddo
    enddo
  enddo

  ! area of top surface
  do ispec = 1,NSPEC2D_TOP
    do i=1,NGLLX
      do j=1,NGLLY
        weight=wxgll(i)*wygll(j)
        area_local_top = area_local_top + dble(jacobian2D_top(i,j,ispec))*weight
      enddo
    enddo
  enddo

  ! use an MPI reduction to compute the total area and volume
  volume_total_region = ZERO
  area_total_bottom   = ZERO
  area_total_top   = ZERO

  call sum_all_dp(area_local_bottom,area_total_bottom)
  call sum_all_dp(area_local_top,area_total_top)
  call sum_all_dp(volume_local,volume_total_region)

  if(myrank == 0) then
    !   sum volume over all the regions
    volume_total = volume_total + volume_total_region

    !   check volume of chunk, and bottom and top area
    write(IMAIN,*)
    write(IMAIN,*) '   calculated top area: ',area_total_top

    ! compare to exact theoretical value
    if(NCHUNKS == 6 .and. .not. TOPOGRAPHY) then
      select case(iregion_code)
        case(IREGION_CRUST_MANTLE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*R_UNIT_SPHERE**2
        case(IREGION_OUTER_CORE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RCMB/R_EARTH)**2
        case(IREGION_INNER_CORE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RICB/R_EARTH)**2
        case default
          call exit_MPI(myrank,'incorrect region code')
      end select
    endif

    write(IMAIN,*) 'calculated bottom area: ',area_total_bottom

    ! compare to exact theoretical value
    if(NCHUNKS == 6 .and. .not. TOPOGRAPHY) then
      select case(iregion_code)
        case(IREGION_CRUST_MANTLE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RCMB/R_EARTH)**2
        case(IREGION_OUTER_CORE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RICB/R_EARTH)**2
        case(IREGION_INNER_CORE)
          write(IMAIN,*) '            more or less similar area (central cube): ', &
                                           dble(NCHUNKS)*(2.*(R_CENTRAL_CUBE / R_EARTH)/sqrt(3.))**2
        case default
          call exit_MPI(myrank,'incorrect region code')
      end select
    endif
    call flush_IMAIN()

  endif

  end subroutine compute_volumes_and_areas

!=====================================================================

  ! compute Earth mass of that part of the slice and then total Earth mass

  subroutine compute_Earth_mass(myrank,Earth_mass_total, &
                            nspec,wxgll,wygll,wzgll,xixstore,xiystore,xizstore, &
                            etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,rhostore,idoubling)

  use constants

  implicit none

  double precision :: Earth_mass_total

  integer :: myrank
  integer :: nspec
  double precision :: wxgll(NGLLX),wygll(NGLLY),wzgll(NGLLZ)

  integer,dimension(nspec) :: idoubling

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec) :: rhostore

  ! local parameters
  double precision :: weight
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  integer :: i,j,k,ispec
  double precision :: Earth_mass_local,Earth_mass_total_region

  ! take into account the fact that the density and the radius of the Earth have previously been non-dimensionalized
  double precision, parameter :: non_dimensionalizing_factor = RHOAV*R_EARTH**3

  ! initializes
  Earth_mass_local = ZERO

  ! calculates volume of all elements in mesh
  do ispec = 1,nspec

    ! suppress fictitious elements in central cube
    if(idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          weight = wxgll(i)*wygll(j)*wzgll(k)

          ! compute the Jacobian
          xixl = xixstore(i,j,k,ispec)
          xiyl = xiystore(i,j,k,ispec)
          xizl = xizstore(i,j,k,ispec)
          etaxl = etaxstore(i,j,k,ispec)
          etayl = etaystore(i,j,k,ispec)
          etazl = etazstore(i,j,k,ispec)
          gammaxl = gammaxstore(i,j,k,ispec)
          gammayl = gammaystore(i,j,k,ispec)
          gammazl = gammazstore(i,j,k,ispec)

          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))

          Earth_mass_local = Earth_mass_local + dble(jacobianl)*rhostore(i,j,k,ispec)*weight

        enddo
      enddo
    enddo
  enddo

  ! take into account the fact that the density and the radius of the Earth have previously been non-dimensionalized
  Earth_mass_local = Earth_mass_local * non_dimensionalizing_factor

  ! use an MPI reduction to compute the total Earth mass
  Earth_mass_total_region = ZERO
  call sum_all_dp(Earth_mass_local,Earth_mass_total_region)

  !   sum volume over all the regions
  if(myrank == 0) Earth_mass_total = Earth_mass_total + Earth_mass_total_region

  end subroutine compute_Earth_mass

!=====================================================================

  ! compute Roland_Sylvain integrals of that part of the slice, and then total integrals for the whole Earth

  subroutine compute_Roland_Sylvain_integr(myrank,Roland_Sylvain_integr_total, &
                            nspec,wxgll,wygll,wzgll,xstore,ystore,zstore,xixstore,xiystore,xizstore, &
                            etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,rhostore,idoubling)

  use constants

  implicit none

  double precision, dimension(9) :: Roland_Sylvain_integr_total

  integer :: myrank
  integer :: nspec
  double precision :: wxgll(NGLLX),wygll(NGLLY),wzgll(NGLLZ)

  integer,dimension(nspec) :: idoubling

  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec) :: rhostore

  ! local parameters
  double precision :: weight
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  double precision :: jacobianl
  integer :: i,j,k,ispec
  double precision :: xval,yval,zval
  double precision :: xval_squared,yval_squared,zval_squared
  double precision :: x_meshpoint,y_meshpoint,z_meshpoint,rho_meshpoint
  double precision :: distance,distance_squared,distance_cubed,distance_fifth_power, &
                      three_over_distance_squared,one_over_distance_cubed,three_over_distance_fifth_power
  double precision :: common_multiplying_factor

  double precision, dimension(9) :: Roland_Sylvain_int_local,Roland_Sylvain_int_total_region
  double precision :: elemental_contribution_1,elemental_contribution_2,elemental_contribution_3, &
                      elemental_contribution_4,elemental_contribution_5,elemental_contribution_6, &
                      elemental_contribution_7,elemental_contribution_8,elemental_contribution_9

  ! take into account the fact that the density and the radius of the Earth have previously been non-dimensionalized
  ! for the gravity vector force, a distance is involved in the dimensions
  double precision, parameter :: nondimensionalizing_factor_gi  = RHOAV * R_EARTH
  ! for the second-order gravity tensor, no distance is involved in the dimensions
  double precision, parameter :: nondimensionalizing_factor_Gij = RHOAV

  double precision, parameter :: scaling_factor_gi  = GRAV * nondimensionalizing_factor_gi
  double precision, parameter :: scaling_factor_Gij = GRAV * nondimensionalizing_factor_Gij

  ! initializes
  Roland_Sylvain_int_local(:) = ZERO

  ! calculates volume of all elements in mesh
  do ispec = 1,nspec

    ! suppress fictitious elements in central cube
    if(idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          weight = wxgll(i)*wygll(j)*wzgll(k)

          ! compute the jacobian
          xixl = xixstore(i,j,k,ispec)
          xiyl = xiystore(i,j,k,ispec)
          xizl = xizstore(i,j,k,ispec)
          etaxl = etaxstore(i,j,k,ispec)
          etayl = etaystore(i,j,k,ispec)
          etazl = etazstore(i,j,k,ispec)
          gammaxl = gammaxstore(i,j,k,ispec)
          gammayl = gammaystore(i,j,k,ispec)
          gammazl = gammazstore(i,j,k,ispec)

!! DK DK do this in double precision for accuracy
          jacobianl = 1.d0 / dble(xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))

    x_meshpoint = xstore(i,j,k,ispec)
    y_meshpoint = ystore(i,j,k,ispec)
    z_meshpoint = zstore(i,j,k,ispec)

    rho_meshpoint = rhostore(i,j,k,ispec)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! beginning of loop on all the data to create
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    xval = x_meshpoint - x_observation
    yval = y_meshpoint - y_observation
    zval = z_meshpoint - z_observation

    xval_squared = xval**2
    yval_squared = yval**2
    zval_squared = zval**2

    distance_squared = xval_squared + yval_squared + zval_squared
    distance = sqrt(distance_squared)
    distance_cubed = distance_squared*distance
    distance_fifth_power = distance_squared*distance_cubed

    three_over_distance_squared = 3.d0 / distance_squared
    one_over_distance_cubed = 1.d0 / distance_cubed
    three_over_distance_fifth_power = three_over_distance_squared * one_over_distance_cubed

! g_x
    elemental_contribution_1 = xval * one_over_distance_cubed

! g_y
    elemental_contribution_2 = yval * one_over_distance_cubed

! g_z
    elemental_contribution_3 = zval * one_over_distance_cubed

! G_xx
    elemental_contribution_4 = (xval_squared * three_over_distance_squared - 1.d0) * one_over_distance_cubed

! G_yy
    elemental_contribution_5 = (yval_squared * three_over_distance_squared - 1.d0) * one_over_distance_cubed

! G_zz
    elemental_contribution_6 = (zval_squared * three_over_distance_squared - 1.d0) * one_over_distance_cubed

! G_xy
    elemental_contribution_7 = xval*yval * three_over_distance_fifth_power

! G_xz
    elemental_contribution_8 = xval*zval * three_over_distance_fifth_power

! G_yz
    elemental_contribution_9 = yval*zval * three_over_distance_fifth_power

    common_multiplying_factor = jacobianl * weight * rho_meshpoint

    Roland_Sylvain_int_local(1) = Roland_Sylvain_int_local(1) + common_multiplying_factor*elemental_contribution_1
    Roland_Sylvain_int_local(2) = Roland_Sylvain_int_local(2) + common_multiplying_factor*elemental_contribution_2
    Roland_Sylvain_int_local(3) = Roland_Sylvain_int_local(3) + common_multiplying_factor*elemental_contribution_3
    Roland_Sylvain_int_local(4) = Roland_Sylvain_int_local(4) + common_multiplying_factor*elemental_contribution_4
    Roland_Sylvain_int_local(5) = Roland_Sylvain_int_local(5) + common_multiplying_factor*elemental_contribution_5
    Roland_Sylvain_int_local(6) = Roland_Sylvain_int_local(6) + common_multiplying_factor*elemental_contribution_6
    Roland_Sylvain_int_local(7) = Roland_Sylvain_int_local(7) + common_multiplying_factor*elemental_contribution_7
    Roland_Sylvain_int_local(8) = Roland_Sylvain_int_local(8) + common_multiplying_factor*elemental_contribution_8
    Roland_Sylvain_int_local(9) = Roland_Sylvain_int_local(9) + common_multiplying_factor*elemental_contribution_9

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! end of loop on all the data to create
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        enddo
      enddo
    enddo
  enddo

  ! multiply by the gravitational constant in S.I. units i.e. in m3 kg-1 s-2
  ! and also take into account the fact that the density and the radius of the Earth have previously been non-dimensionalized
  Roland_Sylvain_int_local(1:3) = Roland_Sylvain_int_local(1:3) * scaling_factor_gi
  Roland_Sylvain_int_local(4:9) = Roland_Sylvain_int_local(4:9) * scaling_factor_Gij

  ! use an MPI reduction to compute the total value of the integral
  Roland_Sylvain_int_total_region(:) = ZERO
!! DK DK could use a single MPI call for the nine values
  call sum_all_dp(Roland_Sylvain_int_local(1),Roland_Sylvain_int_total_region(1))
  call sum_all_dp(Roland_Sylvain_int_local(2),Roland_Sylvain_int_total_region(2))
  call sum_all_dp(Roland_Sylvain_int_local(3),Roland_Sylvain_int_total_region(3))
  call sum_all_dp(Roland_Sylvain_int_local(4),Roland_Sylvain_int_total_region(4))
  call sum_all_dp(Roland_Sylvain_int_local(5),Roland_Sylvain_int_total_region(5))
  call sum_all_dp(Roland_Sylvain_int_local(6),Roland_Sylvain_int_total_region(6))
  call sum_all_dp(Roland_Sylvain_int_local(7),Roland_Sylvain_int_total_region(7))
  call sum_all_dp(Roland_Sylvain_int_local(8),Roland_Sylvain_int_total_region(8))
  call sum_all_dp(Roland_Sylvain_int_local(9),Roland_Sylvain_int_total_region(9))

  !   sum volume over all the regions
  if(myrank == 0) Roland_Sylvain_integr_total(:) = Roland_Sylvain_integr_total(:) + Roland_Sylvain_int_total_region(:)

  end subroutine compute_Roland_Sylvain_integr

