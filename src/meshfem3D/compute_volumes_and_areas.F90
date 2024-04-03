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

  subroutine compute_volumes_and_areas(NCHUNKS,iregion_code,nspec,wxgll,wygll,wzgll, &
                                       xixstore,xiystore,xizstore, &
                                       etaxstore,etaystore,etazstore, &
                                       gammaxstore,gammaystore,gammazstore, &
                                       NSPEC2D_BOTTOM,jacobian2D_bottom,NSPEC2D_TOP,jacobian2D_top,idoubling, &
                                       volume_total,RCMB,RICB,R_CENTRAL_CUBE,RINF)

  use constants, only: NGLLX,NGLLY,NGLLZ,myrank, &
    ZERO,CUSTOM_REAL,PI,R_UNIT_SPHERE,IFLAG_IN_FICTITIOUS_CUBE,IMAIN, &
    IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE, &
    IREGION_TRINFINITE,IREGION_INFINITE

  use shared_parameters, only: R_PLANET
  use meshfem_models_par, only: TOPOGRAPHY

  implicit none

  integer,intent(in) :: nspec
  double precision,intent(in) :: wxgll(NGLLX),wygll(NGLLY),wzgll(NGLLZ)

  integer,intent(in) :: NCHUNKS,iregion_code

  double precision,intent(inout) :: volume_total
  double precision,intent(in) :: RCMB,RICB,R_CENTRAL_CUBE,RINF

  integer,dimension(nspec),intent(in) :: idoubling

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

  integer,intent(in) :: NSPEC2D_BOTTOM,NSPEC2D_TOP
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM),intent(in) :: jacobian2D_bottom
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP),intent(in) :: jacobian2D_top

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
    if (idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

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

          jacobianl = (xixl*(etayl*gammazl-etazl*gammayl) &
                     - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                     + xizl*(etaxl*gammayl-etayl*gammaxl))

          if (jacobianl <= 0.0_CUSTOM_REAL) then
            print *,'Error: rank ',myrank,' found negative Jacobian ',jacobianl,'element',ispec,'ijk',i,j,k,'id',idoubling(ispec)
            print *,'Please check if mesh is okay, exiting...'
            call exit_MPI(myrank,'Error: negative Jacobian found in compute_volumes_and_areas() routine')
          endif

          ! inverts jacobian mapping
          jacobianl = 1._CUSTOM_REAL / jacobianl

          ! sums
          volume_local = volume_local + dble(jacobianl)*weight

        enddo
      enddo
    enddo
  enddo

  ! area of bottom surface
  do ispec = 1,NSPEC2D_BOTTOM
    do i = 1,NGLLX
      do j = 1,NGLLY
        weight = wxgll(i)*wygll(j)
        area_local_bottom = area_local_bottom + dble(jacobian2D_bottom(i,j,ispec))*weight
      enddo
    enddo
  enddo

  ! area of top surface
  do ispec = 1,NSPEC2D_TOP
    do i = 1,NGLLX
      do j = 1,NGLLY
        weight = wxgll(i)*wygll(j)
        area_local_top = area_local_top + dble(jacobian2D_top(i,j,ispec))*weight
      enddo
    enddo
  enddo

  ! use an MPI reduction to compute the total area and volume
  volume_total_region = ZERO
  area_total_bottom   = ZERO
  area_total_top      = ZERO

  call sum_all_dp(area_local_bottom,area_total_bottom)
  call sum_all_dp(area_local_top,area_total_top)
  call sum_all_dp(volume_local,volume_total_region)

  if (myrank == 0) then
    ! sum volume over all the regions
    ! (without transition-to-infinite and infinite region)
    if (iregion_code /= IREGION_TRINFINITE .and. iregion_code /= IREGION_INFINITE) &
      volume_total = volume_total + volume_total_region

    ! check volume of chunk, and bottom and top area
    write(IMAIN,*)
    write(IMAIN,*) 'calculated region volume: ',sngl(volume_total_region)
    write(IMAIN,*) '                top area: ',sngl(area_total_top)

    ! compare to exact theoretical value
    if (NCHUNKS == 6 .and. .not. TOPOGRAPHY) then
      select case (iregion_code)
        case (IREGION_CRUST_MANTLE)
          write(IMAIN,*) '              exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*R_UNIT_SPHERE**2
        case (IREGION_OUTER_CORE)
          write(IMAIN,*) '              exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RCMB/R_PLANET)**2
        case (IREGION_INNER_CORE)
          write(IMAIN,*) '              exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RICB/R_PLANET)**2
        ! TODO: need to fix for transition and infinite layers
        case(IREGION_TRINFINITE)
          write(IMAIN,*) '              exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*RINF**2  ! top at RINF
        case(IREGION_INFINITE)
          write(IMAIN,*) '              exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*RINF**2
        case default
          call exit_MPI(myrank,'incorrect region code')
      end select
    endif

    write(IMAIN,*) '             bottom area: ',sngl(area_total_bottom)

    ! compare to exact theoretical value
    if (NCHUNKS == 6 .and. .not. TOPOGRAPHY) then
      select case (iregion_code)
        case (IREGION_CRUST_MANTLE)
          write(IMAIN,*) '              exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RCMB/R_PLANET)**2
        case (IREGION_OUTER_CORE)
          write(IMAIN,*) '              exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RICB/R_PLANET)**2
        case (IREGION_INNER_CORE)
          write(IMAIN,*) '              more or less similar area (central cube): ', &
                                           dble(NCHUNKS)*(2.*(R_CENTRAL_CUBE / R_PLANET)/sqrt(3.))**2
        ! TODO: need to fix for transition and infinite layers
        case(IREGION_TRINFINITE)
          write(IMAIN,*) '              exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*R_UNIT_SPHERE**2 ! bottom at crust surface
        case(IREGION_INFINITE)
          write(IMAIN,*) '              exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*R_UNIT_SPHERE**2
        case default
          call exit_MPI(myrank,'incorrect region code')
      end select
    endif
    call flush_IMAIN()

  endif

  end subroutine compute_volumes_and_areas

!=====================================================================

  ! compute Earth mass of that part of the slice and then total Earth mass

  subroutine compute_Earth_mass(Earth_mass_total, &
                                Earth_center_of_mass_x_total,Earth_center_of_mass_y_total,Earth_center_of_mass_z_total, &
                                nspec,wxgll,wygll,wzgll,xstore,ystore,zstore,xixstore,xiystore,xizstore, &
                                etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,rhostore,idoubling)

  use constants, only: NGLLX,NGLLY,NGLLZ,ZERO,CUSTOM_REAL,IFLAG_IN_FICTITIOUS_CUBE, &
    SHIFT_TO_THIS_CENTER_OF_MASS,x_shift,y_shift,z_shift, &
    myrank

  use shared_parameters, only: R_PLANET,RHOAV

  implicit none

  double precision :: Earth_mass_total
  double precision :: Earth_center_of_mass_x_total,Earth_center_of_mass_y_total,Earth_center_of_mass_z_total

  integer :: nspec
  double precision :: wxgll(NGLLX),wygll(NGLLY),wzgll(NGLLZ)

  integer,dimension(nspec) :: idoubling

  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec) :: rhostore

  ! local parameters
  double precision :: weight
  double precision :: x_meshpoint,y_meshpoint,z_meshpoint
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl,rhol
  integer :: i,j,k,ispec

  double precision :: Earth_mass_local
  double precision :: Earth_center_of_mass_x_local,Earth_center_of_mass_y_local,Earth_center_of_mass_z_local

  double precision :: Earth_mass_total_region
  double precision :: Earth_center_of_mass_x_tot_reg,Earth_center_of_mass_y_tot_reg,Earth_center_of_mass_z_tot_reg

  ! take into account the fact that the density and the radius of the Earth have previously been non-dimensionalized
  double precision :: non_dimensionalizing_factor1
  double precision :: non_dimensionalizing_factor2

  ! initializes
  Earth_mass_local = ZERO
  Earth_center_of_mass_x_local = ZERO
  Earth_center_of_mass_y_local = ZERO
  Earth_center_of_mass_z_local = ZERO

  ! calculates volume of all elements in mesh
  do ispec = 1,nspec

    ! suppress fictitious elements in central cube
    if (idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

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

          rhol = rhostore(i,j,k,ispec)

          Earth_mass_local = Earth_mass_local + dble(jacobianl) * rhol * weight

          x_meshpoint = xstore(i,j,k,ispec)
          y_meshpoint = ystore(i,j,k,ispec)
          z_meshpoint = zstore(i,j,k,ispec)

          !! for gravity integral calculations we may want to shift the reference frame to a pre-computed center of mass
          if (SHIFT_TO_THIS_CENTER_OF_MASS) then
            x_meshpoint = x_meshpoint - x_shift
            y_meshpoint = y_meshpoint - y_shift
            z_meshpoint = z_meshpoint - z_shift
          endif

          Earth_center_of_mass_x_local = Earth_center_of_mass_x_local + dble(jacobianl)*rhol*x_meshpoint*weight
          Earth_center_of_mass_y_local = Earth_center_of_mass_y_local + dble(jacobianl)*rhol*y_meshpoint*weight
          Earth_center_of_mass_z_local = Earth_center_of_mass_z_local + dble(jacobianl)*rhol*z_meshpoint*weight

        enddo
      enddo
    enddo
  enddo

  ! take into account the fact that the density and the radius of the Earth have previously been non-dimensionalized
  non_dimensionalizing_factor1 = RHOAV*R_PLANET**3
  non_dimensionalizing_factor2 = non_dimensionalizing_factor1 * R_PLANET

  Earth_mass_local = Earth_mass_local * non_dimensionalizing_factor1

  Earth_center_of_mass_x_local = Earth_center_of_mass_x_local * non_dimensionalizing_factor2
  Earth_center_of_mass_y_local = Earth_center_of_mass_y_local * non_dimensionalizing_factor2
  Earth_center_of_mass_z_local = Earth_center_of_mass_z_local * non_dimensionalizing_factor2

  ! use an MPI reduction to compute the total Earth mass
  Earth_mass_total_region = ZERO
  Earth_center_of_mass_x_tot_reg = ZERO
  Earth_center_of_mass_y_tot_reg = ZERO
  Earth_center_of_mass_z_tot_reg = ZERO

  call sum_all_dp(Earth_mass_local,Earth_mass_total_region)
  call sum_all_dp(Earth_center_of_mass_x_local,Earth_center_of_mass_x_tot_reg)
  call sum_all_dp(Earth_center_of_mass_y_local,Earth_center_of_mass_y_tot_reg)
  call sum_all_dp(Earth_center_of_mass_z_local,Earth_center_of_mass_z_tot_reg)

  ! sum volume over all the regions
  if (myrank == 0) then
    Earth_mass_total = Earth_mass_total + Earth_mass_total_region
    Earth_center_of_mass_x_total = Earth_center_of_mass_x_total + Earth_center_of_mass_x_tot_reg
    Earth_center_of_mass_y_total = Earth_center_of_mass_y_total + Earth_center_of_mass_y_tot_reg
    Earth_center_of_mass_z_total = Earth_center_of_mass_z_total + Earth_center_of_mass_z_tot_reg
  endif

  end subroutine compute_Earth_mass

