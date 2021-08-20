!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  8 . 0
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

  subroutine moho_stretching_honor_crust(xelm,yelm,zelm, &
                                         elem_in_crust,elem_in_mantle)

! stretching the moho according to the crust 2.0
! input:  myrank, xelm, yelm, zelm
! Dec, 30, 2009

  use constants, only: myrank, &
    NGNOD,R_UNIT_SPHERE, &
    PI_OVER_TWO,RADIANS_TO_DEGREES,TINYVAL,SMALLVAL,ONE,USE_OLD_VERSION_5_1_5_FORMAT, &
    SUPPRESS_MOHO_STRETCHING,ICRUST_CRUST_SH,ICRUST_SGLOBECRUST

  ! Earth
  use constants, only: EARTH_R,EARTH_R_KM
  ! Mars
  use constants, only: MARS_R,MARS_R_KM
  ! Moon
  use constants, only: MOON_R,MOON_R_KM

  use shared_parameters, only: PLANET_TYPE,IPLANET_EARTH,IPLANET_MARS,IPLANET_MOON,R_PLANET

  use meshfem3D_par, only: &
    RMOHO_FICTITIOUS_IN_MESHER,R220,RMIDDLE_CRUST,REFERENCE_CRUSTAL_MODEL

  use meshfem3D_par, only: &
    TOPOGRAPHY

  implicit none

  double precision,dimension(NGNOD),intent(inout) :: xelm,yelm,zelm
  logical,intent(inout) :: elem_in_crust,elem_in_mantle

  ! local parameters
  double precision :: r,theta,phi,lat,lon
  double precision :: vpvc,vphc,vsvc,vshc,etac,rhoc
  double precision :: c11c,c12c,c13c,c14c,c15c,c16c,c22c,c23c,c24c,c25c,c26c, &
                      c33c,c34c,c35c,c36c,c44c,c45c,c46c,c55c,c56c,c66c
  double precision :: moho,sediment
  double precision :: elevation,gamma
  double precision :: x,y,z
  double precision :: R_moho,R_middlecrust
  integer :: ia,count_crust,count_mantle
  logical :: found_crust,moho_only

  double precision :: MOHO_MAXIMUM_DEFAULT,MOHO_MINIMUM_DEFAULT
  double precision :: MOHO_MAXIMUM,MOHO_MINIMUM

  ! default minimum/maximum allowed moho depths
  ! Earth  (5km/90km non-dimensionalized)
  double precision,parameter :: MOHO_MINIMUM_DEFAULT_EARTH = 5.0 / EARTH_R_KM
  double precision,parameter :: MOHO_MAXIMUM_DEFAULT_EARTH = 90.0 / EARTH_R_KM
  ! Mars
  double precision,parameter :: MOHO_MINIMUM_DEFAULT_MARS = 5.0 / MARS_R_KM
  double precision,parameter :: MOHO_MAXIMUM_DEFAULT_MARS = 150.0 / MARS_R_KM
  ! Moon - todo: needs better estimates
  double precision,parameter :: MOHO_MINIMUM_DEFAULT_MOON = 5.0 / MOON_R_KM
  double precision,parameter :: MOHO_MAXIMUM_DEFAULT_MOON = 100.0 / MOON_R_KM

  ! min/max defaults
  select case(PLANET_TYPE)
  case (IPLANET_MARS)
    ! Mars
    MOHO_MINIMUM_DEFAULT = MOHO_MINIMUM_DEFAULT_MARS
    MOHO_MAXIMUM_DEFAULT = MOHO_MAXIMUM_DEFAULT_MARS
  case (IPLANET_EARTH)
    ! Earth
    MOHO_MINIMUM_DEFAULT = MOHO_MINIMUM_DEFAULT_EARTH
    MOHO_MAXIMUM_DEFAULT = MOHO_MAXIMUM_DEFAULT_EARTH
  case (IPLANET_MOON)
    ! Moon
    MOHO_MINIMUM_DEFAULT = MOHO_MINIMUM_DEFAULT_MOON
    MOHO_MAXIMUM_DEFAULT = MOHO_MAXIMUM_DEFAULT_MOON
  case default
    call exit_MPI(myrank,'Invalid planet, moho stretching not implemented yet')
  end select

  ! sets min/max allowed moho depth
  MOHO_MINIMUM = MOHO_MINIMUM_DEFAULT
  MOHO_MAXIMUM = MOHO_MAXIMUM_DEFAULT

  ! modification for different crustal models
  ! full_sphericalharmonics crustal model
  if (REFERENCE_CRUSTAL_MODEL == ICRUST_CRUST_SH) then
    ! minimum moho < 3.9km
    MOHO_MINIMUM = 3.5d0 / (R_PLANET/1000.d0)
  endif
  ! SGLOBE-rani crustal model
  if (REFERENCE_CRUSTAL_MODEL == ICRUST_SGLOBECRUST) then
    ! minimum moho
    MOHO_MINIMUM = 2.d0 / (R_PLANET/1000.d0)
    ! maximum moho
    MOHO_MAXIMUM = 115.d0 / (R_PLANET/1000.d0)
  endif

  ! radii for stretching criteria
  R_moho = RMOHO_FICTITIOUS_IN_MESHER/R_PLANET
  R_middlecrust = RMIDDLE_CRUST/R_PLANET

  ! loops over element's anchor points
  count_crust = 0
  count_mantle = 0
  do ia = 1,NGNOD
    ! gets anchor point location
    x = xelm(ia)
    y = yelm(ia)
    z = zelm(ia)

    ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
    if (USE_OLD_VERSION_5_1_5_FORMAT) then
      ! note: at this point, the mesh is still perfectly spherical, thus no need to
      !         convert the geocentric colatitude to a geographic colatitude
      call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)
      call reduce(theta,phi)
      lat = (PI_OVER_TWO - theta) * RADIANS_TO_DEGREES
      lon = phi * RADIANS_TO_DEGREES
    else
      ! note: the moho information is given in geographic latitude/longitude (with respect to a reference ellipsoid).
      !       we need to convert the geocentric mesh positions (x/y/z) to geographic ones (lat/lon),
      !       thus correcting geocentric latitude for ellipticity
      !
      ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
      call xyz_2_rlatlon_dble(x,y,z,r,lat,lon)
    endif

    ! sets longitude bounds [-180,180]
    if (lon > 180.d0 ) lon = lon - 360.0d0

    ! initializes
    moho = 0.d0
    sediment = 0.d0
    moho_only = .true.  ! only moho value needed

    ! gets smoothed moho depth
    call meshfem3D_model_crust(lat,lon,r,vpvc,vphc,vsvc,vshc,etac,rhoc,moho,sediment,found_crust,elem_in_crust,moho_only, &
                               c11c,c12c,c13c,c14c,c15c,c16c,c22c,c23c,c24c,c25c,c26c, &
                               c33c,c34c,c35c,c36c,c44c,c45c,c46c,c55c,c56c,c66c)

    !debug
    !if (lon > 140.0 .and. lon < 143. .and. lat <-42. .and. lat > -45. .and. ia == 27) then
    !  print '(a, I4, 4F12.4)', 'crust map:', ia, lat, lon, (ONE-r)*R_PLANET_KM,moho*R_PLANET_KM
    !endif
    !if (r > 1.02d0) then
    !  print '(a,4F12.4)', 'DEBUG moho_stretching_honor_crust_reg(): lat, lon, r, moho =', lat, lon, r, moho
    !  stop 'stop for detecting large r in moho_stretching'
    !endif

    ! checks non-dimensionalized moho depth
    !
    ! note: flag found_crust returns .false. for points below moho,
    !          nevertheless its moho depth should be set and will be used in linear stretching
    if (moho < TINYVAL ) call exit_mpi(myrank,'Error moho depth to honor')

    if (.not. USE_OLD_VERSION_5_1_5_FORMAT) then
      ! limits moho depth to a threshold value to avoid stretching problems
      if (moho < MOHO_MINIMUM) then
        print *,'moho value exceeds minimum (in km): ',moho*R_PLANET/1000.d0,MOHO_MINIMUM*R_PLANET/1000.d0,'lat/lon:',lat,lon
        moho = MOHO_MINIMUM
      endif
      if (moho > MOHO_MAXIMUM) then
        print *,'moho value exceeds maximum (in km): ',moho*R_PLANET/1000.d0,MOHO_MAXIMUM*R_PLANET/1000.d0,'lat/lon:',lat,lon
        moho = MOHO_MAXIMUM
      endif
    endif

    ! radius of moho depth (normalized)
    moho = ONE - moho

    ! checks if moho will be honored by elements
    !
    ! note: we will honor the moho only, if the moho depth is below R_moho (~35km)
    !          or above R_middlecrust (~15km). otherwise, the moho will be "interpolated"
    !          within the element
    if (.not. SUPPRESS_MOHO_STRETCHING) then
      ! globe surface must honor topography for elements to be stretched for moho
      !
      ! note:  if no topography is honored, stretching may lead to distorted elements and invalid Jacobian
      if (TOPOGRAPHY) then
        if (moho < R_moho) then
          ! actual moho below fictitious moho
          ! elements in second layer will stretch down to honor moho topography

          elevation = moho - R_moho

          if (r >= R_moho) then
            ! point above fictitious moho
            ! gamma ranges from 0 (point at surface) to 1 (point at fictitious moho depth)
            gamma = (( R_UNIT_SPHERE - r )/( R_UNIT_SPHERE - R_moho ))
          else
            ! point below fictitious moho
            ! gamma ranges from 0 (point at R220) to 1 (point at fictitious moho depth)
            gamma = (( r - R220/R_PLANET)/( R_moho - R220/R_PLANET))

            ! since not all GLL points are exactly at R220, use a small
            ! tolerance for R220 detection, fix R220
            if (abs(gamma) < SMALLVAL) then
              gamma = 0.0d0
            endif
          endif

          if (gamma < -0.0001d0 .or. gamma > 1.0001d0) &
            call exit_MPI(myrank,'incorrect value of gamma for moho from crust 2.0')

          call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

        else if (moho > R_middlecrust) then
          ! moho above middle crust
          ! elements in first layer will squeeze into crust above moho

          elevation = moho - R_middlecrust

          if (r > R_middlecrust) then
            ! point above middle crust
            ! gamma ranges from 0 (point at surface) to 1 (point at middle crust depth)
            gamma = (R_UNIT_SPHERE-r)/(R_UNIT_SPHERE - R_middlecrust )
          else
            ! point below middle crust
            ! gamma ranges from 0 (point at R220) to 1 (point at middle crust depth)
            gamma = (r - R220/R_PLANET)/( R_middlecrust - R220/R_PLANET )

            ! since not all GLL points are exactly at R220, use a small
            ! tolerance for R220 detection, fix R220
            if (abs(gamma) < SMALLVAL) then
              gamma = 0.0d0
            endif
          endif

          if (gamma < -0.0001d0 .or. gamma > 1.0001d0) &
            call exit_MPI(myrank,'incorrect value of gamma for moho from crust 2.0')

          call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

        endif
      endif ! TOPOGRAPHY
    endif

    ! counts corners in above moho
    ! note: uses a small tolerance
    if (r >= 0.9999d0*moho) then
      count_crust = count_crust + 1
    endif
    ! counts corners below moho
    ! again within a small tolerance
    if (r <= 1.0001d0*moho) then
      count_mantle = count_mantle + 1
    endif

  enddo

  ! sets flag when all corners are above moho
  if (count_crust == NGNOD) then
    elem_in_crust = .true.
  endif
  ! sets flag when all corners are below moho
  if (count_mantle == NGNOD) then
    elem_in_mantle = .true.
  endif

  ! small stretch check: stretching should affect only points above R220
  if (r*R_PLANET < R220) then
    print *,'Error moho stretching: ',r*R_PLANET,R220,moho*R_PLANET
    call exit_mpi(myrank,'incorrect moho stretching')
  endif

  end subroutine moho_stretching_honor_crust

!
!------------------------------------------------------------------------------------------------
!

  subroutine moho_stretching_honor_crust_reg(xelm,yelm,zelm, &
                                             elem_in_crust,elem_in_mantle)

! regional routine: for REGIONAL_MOHO_MESH adaptations
!
! uses a 3-layer crust region
!
! stretching the moho according to the crust 2.0
! input:  myrank, xelm, yelm, zelm
! Dec, 30, 2009

  use constants, only: myrank, &
    NGNOD,PI_OVER_TWO,RADIANS_TO_DEGREES,TINYVAL,ONE,USE_OLD_VERSION_5_1_5_FORMAT, &
    SUPPRESS_MOHO_STRETCHING

  use shared_parameters, only: R_PLANET,HONOR_DEEP_MOHO

  use meshfem3D_par, only: R220

  implicit none

  double precision,dimension(NGNOD),intent(inout) :: xelm,yelm,zelm

  logical,intent(inout) :: elem_in_crust,elem_in_mantle

  ! local parameters
  integer :: ia,count_crust,count_mantle
  double precision :: r,theta,phi,lat,lon
  double precision :: vpvc,vphc,vsvc,vshc,etac,rhoc
  double precision :: c11c,c12c,c13c,c14c,c15c,c16c,c22c,c23c,c24c,c25c,c26c, &
                      c33c,c34c,c35c,c36c,c44c,c45c,c46c,c55c,c56c,c66c
  double precision :: moho,sediment
  logical :: found_crust,moho_only
  double precision :: x,y,z

  ! loops over element's anchor points
  count_crust = 0
  count_mantle = 0
  do ia = 1,NGNOD

    ! anchor point location
    x = xelm(ia)
    y = yelm(ia)
    z = zelm(ia)

    ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
    if (USE_OLD_VERSION_5_1_5_FORMAT) then
      ! note: at this point, the mesh is still perfectly spherical, thus no need to
      !         convert the geocentric colatitude to a geographic colatitude
      call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)
      call reduce(theta,phi)
      lat = (PI_OVER_TWO - theta) * RADIANS_TO_DEGREES
      lon = phi * RADIANS_TO_DEGREES
    else
      ! note: the moho information is given in geographic latitude/longitude (with respect to a reference ellipsoid).
      !       we need to convert the geocentric mesh positions (x/y/z) to geographic ones (lat/lon),
      !       thus correcting geocentric latitude for ellipticity
      !
      ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
      call xyz_2_rlatlon_dble(x,y,z,r,lat,lon)
    endif

    ! sets longitude bounds [-180,180]
    if (lon > 180.d0 ) lon = lon - 360.0d0

    ! initializes
    moho = 0.d0
    sediment = 0.d0
    moho_only = .true.  ! only moho value needed

    ! gets smoothed moho depth
    call meshfem3D_model_crust(lat,lon,r,vpvc,vphc,vsvc,vshc,etac,rhoc,moho,sediment,found_crust,elem_in_crust,moho_only, &
                               c11c,c12c,c13c,c14c,c15c,c16c,c22c,c23c,c24c,c25c,c26c, &
                               c33c,c34c,c35c,c36c,c44c,c45c,c46c,c55c,c56c,c66c)

    ! checks moho depth
    if (abs(moho) < TINYVAL ) call exit_mpi(myrank,'Error moho depth in crust_reg to honor')

    moho = ONE - moho

    ! checks if moho will be honored by elements
    !
    ! note: we will honor the moho, if the moho depth is
    !         - above 15km
    !         - between 25km and 45km
    !         - below 60 km (in HONOR_DEEP_MOHO case)
    !         otherwise, the moho will be "interpolated" within the element
    if (.not. SUPPRESS_MOHO_STRETCHING) then
      ! distinguish between regions with very deep moho (e.g. Himalayan)
      if (HONOR_DEEP_MOHO) then
        call stretch_deep_moho(ia,xelm,yelm,zelm,x,y,z,r,moho)
      else
        call stretch_moho(ia,xelm,yelm,zelm,x,y,z,r,moho)
      endif
    endif

    ! counts corners in above moho
    ! note: uses a small tolerance
    if (r >= 0.9999d0*moho) then
      count_crust = count_crust + 1
    endif
    ! counts corners below moho
    ! again within a small tolerance
    if (r <= 1.0001d0*moho) then
      count_mantle = count_mantle + 1
    endif

  enddo

  ! sets flag when all corners are above moho
  if (count_crust == NGNOD) then
    elem_in_crust = .true.
  endif
  ! sets flag when all corners are below moho
  if (count_mantle == NGNOD) then
    elem_in_mantle = .true.
  endif

  ! small stretch check: stretching should affect only points above R220
  if (r*R_PLANET < R220) then
    print *,'Error moho stretching: ',r*R_PLANET,R220,moho*R_PLANET
    call exit_mpi(myrank,'incorrect moho stretching')
  endif

  end subroutine moho_stretching_honor_crust_reg

!
!-------------------------------------------------------------------------------------------------
!

  subroutine stretch_deep_moho(ia,xelm,yelm,zelm,x,y,z,r,moho)

! honors deep moho (below 60 km), otherwise keeps the mesh boundary at r60 fixed

  use constants, only: NGNOD,EARTH_R,R_UNIT_SPHERE,SMALLVAL
  ! Earth
  use constants, only: EARTH_R

  use shared_parameters, only: PLANET_TYPE,IPLANET_EARTH,IPLANET_MARS,IPLANET_MOON, &
    R_PLANET,RMOHO_STRETCH_ADJUSTMENT

  use meshfem3D_par, only: RMOHO_FICTITIOUS_IN_MESHER,R220,RMIDDLE_CRUST

  implicit none

  integer :: ia

  double precision :: xelm(NGNOD)
  double precision :: yelm(NGNOD)
  double precision :: zelm(NGNOD)

  double precision :: x,y,z

  double precision :: r,moho

  ! local parameters
  double precision :: elevation,gamma
  double precision :: R15,R25,R35,R45,R60,RSURFACE

  ! radii for stretching criteria
  double precision,parameter ::  EARTH_R15 = 6356000.d0/EARTH_R
  double precision,parameter ::  EARTH_R25 = 6346000.d0/EARTH_R
  !double precision,parameter ::  EARTH_R30 = 6341000.d0/EARTH_R
  double precision,parameter ::  EARTH_R35 = 6336000.d0/EARTH_R
  !double precision,parameter ::  EARTH_R40 = 6331000.d0/EARTH_R
  double precision,parameter ::  EARTH_R45 = 6326000.d0/EARTH_R
  !double precision,parameter ::  EARTH_R50 = 6321000.d0/EARTH_R
  !double precision,parameter ::  EARTH_R55 = 6316000.d0/EARTH_R
  double precision,parameter ::  EARTH_R60 = 6311000.d0/EARTH_R

  ! radii
  RSURFACE = R_PLANET
  select case(PLANET_TYPE)
  case (IPLANET_EARTH)
    ! Earth
    R15 = EARTH_R15
    R25 = EARTH_R25
    R35 = EARTH_R35
    R45 = EARTH_R45
    R60 = EARTH_R60
    ! checks moho position: supposed to be at 60 km
    if (RMOHO_STRETCH_ADJUSTMENT /= -20000.d0 ) &
      stop 'wrong moho stretch adjustment for stretch_deep_moho'
    if (RMOHO_FICTITIOUS_IN_MESHER/R_PLANET /= R60 ) &
      stop 'wrong moho depth '
    ! checks middle crust position: supposed to be bottom of first layer at 15 km
    if (RMIDDLE_CRUST/R_PLANET /= R15 ) &
      stop 'wrong middle crust depth'
  case (IPLANET_MARS)
    ! Mars
    R15 = (R_PLANET-15000.d0)/R_PLANET ! 15 km
    R25 = (R_PLANET-25000.d0)/R_PLANET ! 25 km
    R35 = (R_PLANET-45000.d0)/R_PLANET ! 45 km
    R45 = (R_PLANET-65000.d0)/R_PLANET ! 65 km
    R60 = (R_PLANET-120000.d0)/R_PLANET ! 120 km
  case (IPLANET_MOON)
    ! Moon
    R15 = (R_PLANET-12000.d0)/R_PLANET ! 12 km
    R25 = (R_PLANET-28000.d0)/R_PLANET ! 28 km
    R35 = (R_PLANET-38000.d0)/R_PLANET ! 38 km
    R45 = (R_PLANET-48000.d0)/R_PLANET ! 48 km
    R60 = (R_PLANET-60000.d0)/R_PLANET ! 60 km
  case default
    print *,'Invalid planet, stretch deep moho not implemented yet'
    stop 'Invalid planet, stretch deep moho not implemented yet'
  end select

  ! stretches mesh by moving point coordinates
  if (moho < R25 .and. moho > R45) then
    ! moho between r25 and r45

    ! stretches mesh at r35 to moho depth
    elevation = moho - R35
    if (r >= R35 .and. r < R15) then
      gamma=((R15-r)/(R15-R35))
    else if (r < R35 .and. r > R60) then
      gamma = (( r - R60)/( R35 - R60)) ! keeps r60 fixed
      if (abs(gamma) < SMALLVAL) then
        gamma=0.0d0
      endif
    else
      gamma=0.0d0
    endif
    if (gamma < -0.0001d0 .or. gamma > 1.0001d0) &
      stop 'incorrect value of gamma for moho from crust 2.0'

    call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

  else if (moho < R45) then
    ! moho below r45

    ! moves mesh at r35 down to r45
    elevation = R45 - R35
    if (r >= R35 .and. r < R15) then
      gamma=((R15-r)/(R15-R35)) ! moves r35 down to r45
    else if (r < R35 .and. r > R60) then
      gamma=((r-R60)/(R35-R60)) ! keeps r60 fixed
      if (abs(gamma) < SMALLVAL) then
        gamma=0.0d0
      endif
    else
      gamma=0.0d0
    endif
    if (gamma < -0.0001d0 .or. gamma > 1.0001d0) &
      stop 'incorrect value of gamma for moho from crust 2.0'

    call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

    ! add deep moho here
    if (moho < R60) then
      ! moho below r60

      ! stretches mesh at r60 to moho
      elevation = moho - R60
      if (r < R45 .and. r >= R60) then
        gamma=(R45-r)/(R45-R60)
      else if (r < R60) then
        gamma=(r-R220/RSURFACE)/(R60-R220/RSURFACE)
        if (abs(gamma) < SMALLVAL) then
          gamma=0.0d0
        endif
      else
        gamma=0.0d0
      endif

      call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)
    endif

  else if (moho > R25) then
    ! moho above r25

    ! moves mesh at r35 up to r25
    elevation = R25-R35
    if (r >= R35 .and. r < R15) then
      gamma=((R15-r)/(R15-R35)) ! stretches r35 up to r25
    else if (r < R35 .and. r > R60) then
      gamma=(r-R60)/(R35-R60) ! keeps r60 fixed
      if (abs(gamma) < SMALLVAL) then
        gamma=0.0d0
      endif
    else
      gamma=0.0d0
    endif
    if (gamma < -0.0001d0 .or. gamma > 1.0001d0) &
      stop 'incorrect value of gamma for moho from crust 2.0'

    call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

    ! add shallow moho here
    if (moho > R15) then
      ! moho above r15

      ! stretches mesh at r15 to moho depth
      elevation = moho-R15
      if (r >= R15) then
        gamma=(R_UNIT_SPHERE-r)/(R_UNIT_SPHERE-R15)
      else if (r < R15 .and. R > R25) then
        gamma=(r-R25)/(R15-R25) ! keeps mesh at r25 fixed
        if (abs(gamma) < SMALLVAL) then
          gamma=0.0d0
        endif
      else
        gamma=0.0d0
      endif

      call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)
    endif
  endif

  end subroutine stretch_deep_moho

!
!-------------------------------------------------------------------------------------------------
!

  subroutine stretch_moho(ia,xelm,yelm,zelm,x,y,z,r,moho)

! honors shallow and middle depth moho, deep moho will be interpolated within elements
! mesh will get stretched down to r220

  use constants, only: NGNOD,EARTH_R,R_UNIT_SPHERE,SMALLVAL
  ! Earth
  use constants, only: EARTH_R

  use shared_parameters, only: PLANET_TYPE,IPLANET_EARTH,IPLANET_MARS,IPLANET_MOON, &
    R_PLANET,RMOHO_STRETCH_ADJUSTMENT

  use meshfem3D_par, only: RMOHO_FICTITIOUS_IN_MESHER,R220,RMIDDLE_CRUST

  implicit none

  integer :: ia

  double precision :: xelm(NGNOD)
  double precision :: yelm(NGNOD)
  double precision :: zelm(NGNOD)

  double precision :: r,moho
  double precision :: x,y,z

  ! local parameters
  double precision :: elevation,gamma
  double precision :: R15,R25,R35,R45,R55,RSURFACE

  ! radii for stretching criteria
  double precision,parameter ::  EARTH_R15 = 6356000.d0/EARTH_R
  double precision,parameter ::  EARTH_R25 = 6346000.d0/EARTH_R
  !double precision,parameter ::  EARTH_R30 = 6341000.d0/EARTH_R
  double precision,parameter ::  EARTH_R35 = 6336000.d0/EARTH_R
  !double precision,parameter ::  EARTH_R40 = 6331000.d0/EARTH_R
  double precision,parameter ::  EARTH_R45 = 6326000.d0/EARTH_R
  !double precision,parameter ::  EARTH_R50 = 6321000.d0/EARTH_R
  double precision,parameter ::  EARTH_R55 = 6316000.d0/EARTH_R
  !double precision,parameter ::  EARTH_R60 = 6311000.d0/EARTH_R

  ! radii
  RSURFACE = R_PLANET
  select case(PLANET_TYPE)
  case (IPLANET_EARTH)
    ! Earth
    R15 = EARTH_R15
    R25 = EARTH_R25
    R35 = EARTH_R35
    R45 = EARTH_R45
    R55 = EARTH_R55
    ! checks moho position: supposed to be at 55 km
    if (RMOHO_STRETCH_ADJUSTMENT /= -15000.d0 ) &
      stop 'wrong moho stretch adjustment for stretch_moho'
    if (RMOHO_FICTITIOUS_IN_MESHER/R_PLANET /= R55 ) &
      stop 'wrong moho depth '
    ! checks middle crust position: supposed to be bottom of first layer at 15 km
    if (RMIDDLE_CRUST/R_PLANET /= R15 ) &
      stop 'wrong middle crust depth'
  case (IPLANET_MARS)
    ! Mars
    R15 = (R_PLANET-15000.d0)/R_PLANET ! 15 km
    R25 = (R_PLANET-25000.d0)/R_PLANET ! 25 km
    R35 = (R_PLANET-45000.d0)/R_PLANET ! 45 km
    R45 = (R_PLANET-65000.d0)/R_PLANET ! 65 km
    R55 = (R_PLANET-95000.d0)/R_PLANET ! 95 km
  case (IPLANET_MOON)
    ! Moon
    R15 = (R_PLANET-12000.d0)/R_PLANET ! 12 km
    R25 = (R_PLANET-28000.d0)/R_PLANET ! 28 km
    R35 = (R_PLANET-38000.d0)/R_PLANET ! 38 km
    R45 = (R_PLANET-48000.d0)/R_PLANET ! 48 km
    R55 = (R_PLANET-55000.d0)/R_PLANET ! 55 km
  case default
    print *,'Invalid planet, stretch moho not implemented yet'
    stop 'Invalid planet, stretch moho not implemented yet'
  end select

  ! moho between 25km and 45 km
  if (moho < R25 .and. moho > R45) then

    elevation = moho - R35
    if (r >= R35 .and. r < R15) then
      gamma=((R15-r)/(R15-R35))
    else if (r < R35 .and. r > R220/RSURFACE) then
      gamma = ((r-R220/RSURFACE)/(R35-R220/RSURFACE))
      if (abs(gamma) < SMALLVAL) then
        gamma=0.0d0
      endif
    else
      gamma=0.0d0
    endif
    if (gamma < -0.0001d0 .or. gamma > 1.0001d0) &
      stop 'incorrect value of gamma for moho from crust 2.0'

    call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

  else if (moho < R45) then
    ! moho below 45 km

    ! moves mesh at r35 down to r45
    elevation = R45 - R35
    if (r >= R35 .and. r < R15) then
      gamma=((R15-r)/(R15-R35))
    else if (r < R35 .and. r > R220/RSURFACE) then
      gamma=((r-R220/RSURFACE)/(R35-R220/RSURFACE))
      if (abs(gamma) < SMALLVAL) then
        gamma=0.0d0
      endif
    else
      gamma=0.0d0
    endif
    if (gamma < -0.0001d0 .or. gamma > 1.0001d0) &
      stop 'incorrect value of gamma for moho from crust 2.0'

    call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

  else if (moho > R25) then
    ! moho above 25km

    ! moves mesh at r35 up to r25
    elevation = R25-R35
    if (r >= R35 .and. r < R15) then
      gamma=((R15-r)/(R15-R35))
    else if (r < R35 .and. r > R220/RSURFACE) then
      gamma=(r-R220/RSURFACE)/(R35-R220/RSURFACE)
      if (abs(gamma) < SMALLVAL) then
        gamma=0.0d0
      endif
    else
      gamma=0.0d0
    endif
    if (gamma < -0.0001d0 .or. gamma > 1.0001d0) &
      stop 'incorrect value of gamma for moho from crust 2.0'

    call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

    ! add shallow moho here
    if (moho > R15) then
      elevation = moho-R15
      if (r >= R15) then
        gamma=(R_UNIT_SPHERE-r)/(R_UNIT_SPHERE-R15)
      else if (r < R15 .and. R > R25) then
        gamma=(r-R25)/(R15-R25)
        if (abs(gamma) < SMALLVAL) then
          gamma=0.0d0
        endif
      else
        gamma=0.0d0
      endif

      call move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)
    endif
  endif

  end subroutine stretch_moho

!
!-------------------------------------------------------------------------------------------------
!

  subroutine move_point(ia,xelm,yelm,zelm,x,y,z,gamma,elevation,r)

! moves a point to a new location defined by gamma,elevation and r

  use constants, only: NGNOD,ONE

  implicit none

  integer,intent(in) :: ia

  double precision,intent(out) :: xelm(NGNOD)
  double precision,intent(out) :: yelm(NGNOD)
  double precision,intent(out) :: zelm(NGNOD)

  double precision,intent(inout) :: x,y,z

  double precision,intent(in) :: elevation,gamma
  double precision,intent(inout) :: r

  ! local parameters
  double precision :: stretch_factor

  !  stretch factor
  ! offset will be gamma * elevation
  ! scaling Cartesian coordinates xyz rather than spherical r/theta/phi involves division of offset by r
  stretch_factor = ONE + gamma * elevation/r

  ! new point location
  x = x * stretch_factor
  y = y * stretch_factor
  z = z * stretch_factor

  ! stores new point location
  xelm(ia) = x
  yelm(ia) = y
  zelm(ia) = z

  ! new radius
  r = dsqrt(xelm(ia)*xelm(ia) + yelm(ia)*yelm(ia) + zelm(ia)*zelm(ia))

  end subroutine move_point
