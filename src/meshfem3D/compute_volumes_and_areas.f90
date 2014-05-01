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

