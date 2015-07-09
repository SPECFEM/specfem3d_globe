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

  subroutine compute_volumes_and_areas(myrank,NCHUNKS,iregion_code,nspec,wxgll,wygll,wzgll, &
                            xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
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

          volume_local = volume_local + dble(jacobianl)*weight

        enddo
      enddo
    enddo
  enddo

  ! area of bottom surface
  do ispec = 1,NSPEC2D_BOTTOM
    do i = 1,NGLLX
      do j = 1,NGLLY
        weight=wxgll(i)*wygll(j)
        area_local_bottom = area_local_bottom + dble(jacobian2D_bottom(i,j,ispec))*weight
      enddo
    enddo
  enddo

  ! area of top surface
  do ispec = 1,NSPEC2D_TOP
    do i = 1,NGLLX
      do j = 1,NGLLY
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

  if (myrank == 0) then
    !   sum volume over all the regions
    volume_total = volume_total + volume_total_region

    !   check volume of chunk, and bottom and top area
    write(IMAIN,*)
    write(IMAIN,*) '   calculated top area: ',area_total_top

    ! compare to exact theoretical value
    if (NCHUNKS == 6 .and. .not. TOPOGRAPHY) then
      select case (iregion_code)
        case (IREGION_CRUST_MANTLE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*R_UNIT_SPHERE**2
        case (IREGION_OUTER_CORE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RCMB/R_EARTH)**2
        case (IREGION_INNER_CORE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RICB/R_EARTH)**2
        case default
          call exit_MPI(myrank,'incorrect region code')
      end select
    endif

    write(IMAIN,*) 'calculated bottom area: ',area_total_bottom

    ! compare to exact theoretical value
    if (NCHUNKS == 6 .and. .not. TOPOGRAPHY) then
      select case (iregion_code)
        case (IREGION_CRUST_MANTLE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RCMB/R_EARTH)**2
        case (IREGION_OUTER_CORE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RICB/R_EARTH)**2
        case (IREGION_INNER_CORE)
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
                            Earth_center_of_mass_x_total,Earth_center_of_mass_y_total,Earth_center_of_mass_z_total, &
                            nspec,wxgll,wygll,wzgll,xstore,ystore,zstore,xixstore,xiystore,xizstore, &
                            etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,rhostore,idoubling)

  use constants

  implicit none

  double precision :: Earth_mass_total
  double precision :: Earth_center_of_mass_x_total,Earth_center_of_mass_y_total,Earth_center_of_mass_z_total

  integer :: myrank
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
  double precision, parameter :: non_dimensionalizing_factor1 = RHOAV*R_EARTH**3
  double precision, parameter :: non_dimensionalizing_factor2 = non_dimensionalizing_factor1 * R_EARTH

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

          Earth_mass_local = Earth_mass_local + dble(jacobianl)*rhol*weight

          x_meshpoint = xstore(i,j,k,ispec)
          y_meshpoint = ystore(i,j,k,ispec)
          z_meshpoint = zstore(i,j,k,ispec)

!! DK DK for gravity integral calculations we may want to shift the reference frame to a pre-computed center of mass
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

!=====================================================================

  ! compute gravity integrals of that part of the slice, and then total integrals for the whole Earth

  subroutine compute_gravity_integrals(myrank,iregion_code,nspec,wxgll,wygll,wzgll,xstore,ystore,zstore, &
                xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,rhostore,idoubling)

  use constants

  use meshfem3D_par,only: OUTPUT_FILES, &
#ifdef FORCE_VECTORIZATION
     x_observation1D,y_observation1D,z_observation1D,g_x1D,g_y1D,g_z1D,G_xx1D,G_yy1D,G_zz1D,G_xy1D,G_xz1D,G_yz1D
#else
     x_observation,y_observation,z_observation,g_x,g_y,g_z,G_xx,G_yy,G_zz,G_xy,G_xz,G_yz
#endif

  implicit none

  integer :: myrank,iregion_code,nspec
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
  double precision :: x_meshpoint,y_meshpoint,z_meshpoint
  double precision :: distance_squared,distance_cubed, &
                      three_over_distance_squared,one_over_distance_cubed,three_over_distance_fifth_power
  double precision :: common_multiplying_factor,common_mult_times_one_over,common_mult_times_three_over

  ! name of the timestamp files
  character(len=MAX_STRING_LEN) :: outputname

#ifdef FORCE_VECTORIZATION
  integer :: ix_iy_ichunk
#else
  integer :: ix,iy,ichunk
#endif

  ! if we do not want to compute the gravity integrals, only the center of mass (computed before)
  if (ONLY_COMPUTE_CENTER_OF_MASS) return

  ! calculates volume of all elements in mesh
  do ispec = 1,nspec

    ! print information about number of elements done so far
    if (myrank == 0 .and. (mod(ispec,NSPEC_DISPLAY_INTERVAL) == 0 .or. ispec == 1 .or. ispec == nspec)) then
       write(IMAIN,*) 'for gravity integrals ',ispec,' elements computed out of ',nspec
       ! write time stamp file to give information about progression of simulation
       write(outputname,"('/timestamp_reg',i1.1,'_ispec',i7.7,'_out_of_',i7.7)") iregion_code,ispec,nspec
       ! timestamp file output
       open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown',action='write')
       write(IOUT,*) ispec,' elements done out of ',nspec,' in region ',iregion_code
       close(unit=IOUT)
    endif

    ! suppress fictitious elements in central cube
    if (idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

    ! see if we compute the contribution of the crust only
    if (COMPUTE_CRUST_CONTRIB_ONLY .and. idoubling(ispec) /= IFLAG_CRUST) cycle

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

          ! do this in double precision for accuracy
          jacobianl = 1.d0 / dble(xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))

          if (CHECK_FOR_NEGATIVE_JACOBIANS .and. jacobianl <= ZERO) stop 'Error: negative Jacobian found in integral calculation'

          x_meshpoint = xstore(i,j,k,ispec)
          y_meshpoint = ystore(i,j,k,ispec)
          z_meshpoint = zstore(i,j,k,ispec)

!! DK DK for gravity integral calculations we may want to shift the reference frame to a pre-computed center of mass
          if (SHIFT_TO_THIS_CENTER_OF_MASS) then
            x_meshpoint = x_meshpoint - x_shift
            y_meshpoint = y_meshpoint - y_shift
            z_meshpoint = z_meshpoint - z_shift
          endif

          common_multiplying_factor = jacobianl * weight * rhostore(i,j,k,ispec)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! beginning of loop on all the data to create
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! loop on all the chunks and then on all the observation nodes in each chunk

#ifdef FORCE_VECTORIZATION
! this works only if the arrays are contiguous in memory (which is always the case for static arrays, as used in the code)
! but the code produced is extremely fast because we get a single and fully-vectorized loop on the whole array
  do ix_iy_ichunk = 1,NTOTAL_OBSERVATION

    xval = x_meshpoint - x_observation1D(ix_iy_ichunk)
    yval = y_meshpoint - y_observation1D(ix_iy_ichunk)
    zval = z_meshpoint - z_observation1D(ix_iy_ichunk)

    xval_squared = xval**2
    yval_squared = yval**2
    zval_squared = zval**2

    distance_squared = xval_squared + yval_squared + zval_squared
    distance_cubed = distance_squared * sqrt(distance_squared)

    three_over_distance_squared = 3.d0 / distance_squared
    one_over_distance_cubed = 1.d0 / distance_cubed
    three_over_distance_fifth_power = three_over_distance_squared * one_over_distance_cubed

    common_mult_times_one_over = common_multiplying_factor * one_over_distance_cubed
    common_mult_times_three_over = common_multiplying_factor * three_over_distance_fifth_power

    g_x1D(ix_iy_ichunk) = g_x1D(ix_iy_ichunk) + common_mult_times_one_over * xval
    g_y1D(ix_iy_ichunk) = g_y1D(ix_iy_ichunk) + common_mult_times_one_over * yval
    g_z1D(ix_iy_ichunk) = g_z1D(ix_iy_ichunk) + common_mult_times_one_over * zval

    G_xx1D(ix_iy_ichunk) = G_xx1D(ix_iy_ichunk) + common_mult_times_one_over * (xval_squared * three_over_distance_squared - 1.d0)
    G_yy1D(ix_iy_ichunk) = G_yy1D(ix_iy_ichunk) + common_mult_times_one_over * (yval_squared * three_over_distance_squared - 1.d0)
    G_zz1D(ix_iy_ichunk) = G_zz1D(ix_iy_ichunk) + common_mult_times_one_over * (zval_squared * three_over_distance_squared - 1.d0)

    G_xy1D(ix_iy_ichunk) = G_xy1D(ix_iy_ichunk) + common_mult_times_three_over * xval*yval
    G_xz1D(ix_iy_ichunk) = G_xz1D(ix_iy_ichunk) + common_mult_times_three_over * xval*zval
    G_yz1D(ix_iy_ichunk) = G_yz1D(ix_iy_ichunk) + common_mult_times_three_over * yval*zval

  enddo

#else
  do ichunk = 1,NCHUNKS_MAX
    do iy = 1,NY_OBSERVATION
      do ix = 1,NX_OBSERVATION

    xval = x_meshpoint - x_observation(ix,iy,ichunk)
    yval = y_meshpoint - y_observation(ix,iy,ichunk)
    zval = z_meshpoint - z_observation(ix,iy,ichunk)

    xval_squared = xval**2
    yval_squared = yval**2
    zval_squared = zval**2

    distance_squared = xval_squared + yval_squared + zval_squared
    distance_cubed = distance_squared * sqrt(distance_squared)

    three_over_distance_squared = 3.d0 / distance_squared
    one_over_distance_cubed = 1.d0 / distance_cubed
    three_over_distance_fifth_power = three_over_distance_squared * one_over_distance_cubed

    common_mult_times_one_over = common_multiplying_factor * one_over_distance_cubed
    common_mult_times_three_over = common_multiplying_factor * three_over_distance_fifth_power

    g_x(ix,iy,ichunk) = g_x(ix,iy,ichunk) + common_mult_times_one_over * xval
    g_y(ix,iy,ichunk) = g_y(ix,iy,ichunk) + common_mult_times_one_over * yval
    g_z(ix,iy,ichunk) = g_z(ix,iy,ichunk) + common_mult_times_one_over * zval

    G_xx(ix,iy,ichunk) = G_xx(ix,iy,ichunk) + common_mult_times_one_over * (xval_squared * three_over_distance_squared - 1.d0)
    G_yy(ix,iy,ichunk) = G_yy(ix,iy,ichunk) + common_mult_times_one_over * (yval_squared * three_over_distance_squared - 1.d0)
    G_zz(ix,iy,ichunk) = G_zz(ix,iy,ichunk) + common_mult_times_one_over * (zval_squared * three_over_distance_squared - 1.d0)

    G_xy(ix,iy,ichunk) = G_xy(ix,iy,ichunk) + common_mult_times_three_over * xval*yval
    G_xz(ix,iy,ichunk) = G_xz(ix,iy,ichunk) + common_mult_times_three_over * xval*zval
    G_yz(ix,iy,ichunk) = G_yz(ix,iy,ichunk) + common_mult_times_three_over * yval*zval

      enddo
    enddo
  enddo
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! end of loop on all the data to create
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        enddo
      enddo
    enddo
  enddo

  end subroutine compute_gravity_integrals

