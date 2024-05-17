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


  subroutine get_timestep_and_layers()

  use constants
  use shared_parameters

  implicit none

  ! local variables
  integer :: NEX_MAX
  integer :: multiplication_factor
  double precision :: min_chunk_width_in_degrees
  double precision :: dt_auto
  double precision :: MIN_GLL_POINT_SPACING,MIN_GLL_POINT_SPACING_NGLL5
  integer :: nex_max_auto_ner_estimate

  ! initializes
  DT = 0.d0
  MIN_ATTENUATION_PERIOD = 0.d0
  MAX_ATTENUATION_PERIOD = 0.d0
  ATT_F_C_SOURCE         = 0.d0

  ! minimum period estimation
  call get_minimum_period_estimate()

  !----
  !----  case prem_onecrust by default
  !----

  NEX_MAX = max(NEX_XI,NEX_ETA)

  ! to suppress the crustal layers
  ! (replaced by an extension of the mantle: R_PLANET is not modified, but no more crustal doubling)
  if (SUPPRESS_CRUSTAL_MESH) then
    multiplication_factor = 2
  else
    multiplication_factor = 1
  endif

  ! sets time step size, attenuation band and layers
  select case (PLANET_TYPE)
  case (IPLANET_MARS)
    ! Mars
    ! radius of central cube
    R_CENTRAL_CUBE = 395000.d0

    ! Mars minimum period:
    !   Mars deg2km ~ 59km, vs_min ~ 4.0 km/s
    !   90./NEX * 59.0 / (NGLL-1) * 4 / 4.0 -> NEX = 80, NGLL = 5: T_min ~ 16.59

    ! time step / number of element layers
    if (NEX_MAX*multiplication_factor <= 80) then
      DT                       = 0.2d0
      ! attenuation period range
      MIN_ATTENUATION_PERIOD   = 20.d0
      MAX_ATTENUATION_PERIOD   = 1000.d0
      ! number of element layers in each mesh region
      NER_CRUST                = 1
      NER_80_MOHO              = 2
      NER_220_80               = 1
      NER_400_220              = 2
      NER_600_400              = 2
      NER_670_600              = 1
      NER_771_670              = 3
      NER_TOPDDOUBLEPRIME_771  = 3
      NER_CMB_TOPDDOUBLEPRIME  = 1
      NER_OUTER_CORE           = 6
      NER_TOP_CENTRAL_CUBE_ICB = 1
    else if (NEX_MAX*multiplication_factor <= 96) then
      DT                       = 0.2d0
      ! attenuation period range
      MIN_ATTENUATION_PERIOD   = 20.d0
      MAX_ATTENUATION_PERIOD   = 1000.d0
      ! number of element layers in each mesh region
      NER_CRUST                = 5  ! 110 km
      NER_80_MOHO              = 2  ! 225 km
      NER_220_80               = 1  ! 147 km
      NER_400_220              = 2  ! 253 km
      NER_600_400              = 2  ! 200 km
      NER_670_600              = 1  !  95 km
      NER_771_670              = 3  ! 313 km
      NER_TOPDDOUBLEPRIME_771  = 4  ! 530 km
      NER_CMB_TOPDDOUBLEPRIME  = 1  !  35 km
      NER_OUTER_CORE           = 10 !1068 km
      NER_TOP_CENTRAL_CUBE_ICB = 2  ! 300 km
    else if (NEX_MAX*multiplication_factor <= 160) then
      DT                       = 0.2d0
      ! attenuation period range
      ! The shortest period at this resoution is about 20 s compare with
      ! NEX=256, so change the attenuation period accordingly.
      ! attenuation period range
      MIN_ATTENUATION_PERIOD   = 20.d0
      MAX_ATTENUATION_PERIOD   = 1000.d0
      ! number of element layers in each mesh region
      NER_CRUST                = 5 !3   change to 5 for moho stretching
      NER_80_MOHO              = 5
      NER_220_80               = 3
      NER_400_220              = 4
      NER_600_400              = 3
      NER_670_600              = 2
      NER_771_670              = 6
      NER_TOPDDOUBLEPRIME_771  = 8
      NER_CMB_TOPDDOUBLEPRIME  = 1
      NER_OUTER_CORE           = 16
      NER_TOP_CENTRAL_CUBE_ICB = 4
    else if (NEX_MAX*multiplication_factor <= 256) then
      DT                       = 0.15d0
      ! attenuation period range
      MIN_ATTENUATION_PERIOD   = 10.d0
      MAX_ATTENUATION_PERIOD   = 500.d0
      ! number of element layers in each mesh region
      NER_CRUST                = 5
      NER_80_MOHO              = 6
      NER_220_80               = 4
      NER_400_220              = 5
      NER_600_400              = 4
      NER_670_600              = 2
      NER_771_670              = 8
      NER_TOPDDOUBLEPRIME_771  = 10
      NER_CMB_TOPDDOUBLEPRIME  = 1
      NER_OUTER_CORE           = 18
      NER_TOP_CENTRAL_CUBE_ICB = 5
    else if (NEX_MAX*multiplication_factor <= 320) then
      DT                       = 0.07d0
      ! attenuation period range
      MIN_ATTENUATION_PERIOD   = 10.d0
      MAX_ATTENUATION_PERIOD   = 500.d0
      ! number of element layers in each mesh region
      NER_CRUST                = 5
      NER_80_MOHO              = 8
      NER_220_80               = 5
      NER_400_220              = 6
      NER_600_400              = 5
      NER_670_600              = 3
      NER_771_670              = 10
      NER_TOPDDOUBLEPRIME_771  = 13
      NER_CMB_TOPDDOUBLEPRIME  = 1
      NER_OUTER_CORE           = 20
      NER_TOP_CENTRAL_CUBE_ICB = 5
    else if (NEX_MAX*multiplication_factor <= 480) then
      DT                       = 0.008d0
      ! attenuation period range
      MIN_ATTENUATION_PERIOD   = 10.d0
      MAX_ATTENUATION_PERIOD   = 500.d0
      ! number of element layers in each mesh region
      NER_CRUST                = 7
      NER_80_MOHO              = 14
      NER_220_80               = 11
      NER_400_220              = 12
      NER_600_400              = 11
      NER_670_600              = 9
      NER_771_670              = 16
      NER_TOPDDOUBLEPRIME_771  = 19
      NER_CMB_TOPDDOUBLEPRIME  = 1
      NER_OUTER_CORE           = 30
      NER_TOP_CENTRAL_CUBE_ICB = 9
    else
      ! for bigger NEX, uses 480 setting and then automatically adjusts in auto_ner()
      DT                       = 0.008d0
      ! attenuation period range
      MIN_ATTENUATION_PERIOD   = 10.d0
      MAX_ATTENUATION_PERIOD   = 500.d0
      ! number of element layers in each mesh region
      NER_CRUST                = 7
      NER_80_MOHO              = 14
      NER_220_80               = 11
      NER_400_220              = 12
      NER_600_400              = 11
      NER_670_600              = 9
      NER_771_670              = 16
      NER_TOPDDOUBLEPRIME_771  = 19
      NER_CMB_TOPDDOUBLEPRIME  = 1
      NER_OUTER_CORE           = 30
      NER_TOP_CENTRAL_CUBE_ICB = 9
    endif
    ! uses automatic estimates for NEX > 480
    nex_max_auto_ner_estimate = 480

    if (HONOR_1D_SPHERICAL_MOHO) then
      ! 1D models honor 1D spherical moho
      if (.not. ONE_CRUST) then
        ! case 1D + two+ crustal layers
        if (NEX_MAX*multiplication_factor <= 160) then
          DT = 0.16d0
        else if (NEX_MAX*multiplication_factor <= 256) then
          DT = 0.09d0
        else if (NEX_MAX*multiplication_factor <= 320) then
          DT = 0.07d0
        endif
      else
        ! makes time set to maximum value without blow up simulation with 3D crust
        if (NEX_MAX*multiplication_factor <= 80) then
          DT = 0.15d0
        else if (NEX_MAX*multiplication_factor <= 160) then
          ! make it smaller to avoid stability issue
          DT = 0.16d0
        else if (NEX_MAX*multiplication_factor <= 256) then
          DT = 0.1d0 !0.12
        else if (NEX_MAX*multiplication_factor <= 320) then
          DT = 0.08d0
        endif
      endif
    endif

  case (IPLANET_MOON)
    ! Moon
    ! radius of central cube
    R_CENTRAL_CUBE = 200000.d0   ! artificial ICB will be at 250 km

    ! Moon minimum period:
    !   Moon deg2km ~ 1737.1 * 2 * PI / 360 = 30.2km, vs_min ~ 4.0 km/s
    !   90./NEX * 30.2 / (NGLL-1) * 4 / 4.0 -> NEX = 80, NGLL = 5: T_min ~ 8.5 s

    ! time step / number of element layers
    if (NEX_MAX*multiplication_factor <= 80) then
      DT                       = 0.16d0
      ! attenuation period range
      MIN_ATTENUATION_PERIOD   = 15.d0
      MAX_ATTENUATION_PERIOD   = 750d0
      ! number of element layers in each mesh region
      NER_CRUST                = 2   ! earth:  1
      NER_80_MOHO              = 2   !         1
      NER_220_80               = 3   !         2
      NER_400_220              = 2   !         2
      NER_600_400              = 2   !         2
      NER_670_600              = 1   !         1
      NER_771_670              = 1   !         1
      NER_TOPDDOUBLEPRIME_771  = 4   !         15
      NER_CMB_TOPDDOUBLEPRIME  = 2   !         1
      NER_OUTER_CORE           = 4   !         16
      NER_TOP_CENTRAL_CUBE_ICB = 2   !         2
    else
      ! for bigger NEX, automatically adjusts in auto_ner()..
      DT                       = 0.1d0
      ! attenuation period range
      MIN_ATTENUATION_PERIOD   = 15.d0
      MAX_ATTENUATION_PERIOD   = 750d0
      ! number of element layers in each mesh region
      NER_CRUST                = 3
      NER_80_MOHO              = 2  ! 2 okay; setting w/ jacobian error: nex96, rmoho = 28km, r80 = 130km, ner > 2
      NER_220_80               = 3
      NER_400_220              = 2
      NER_600_400              = 2
      NER_670_600              = 1
      NER_771_670              = 1
      NER_TOPDDOUBLEPRIME_771  = 5
      NER_CMB_TOPDDOUBLEPRIME  = 2
      NER_OUTER_CORE           = 4
      NER_TOP_CENTRAL_CUBE_ICB = 2
    endif
    ! uses automatic estimates for NEX > 80
    nex_max_auto_ner_estimate = 80

  case (IPLANET_EARTH)
    ! Earth
    ! default radius of central cube
    R_CENTRAL_CUBE = 982000.0d0

    ! default layering
    ! number of element layers in each mesh region
    NER_CRUST                = 1
    NER_80_MOHO              = 1
    NER_220_80               = 2
    NER_400_220              = 2
    NER_600_400              = 2
    NER_670_600              = 1
    NER_771_670              = 1
    NER_TOPDDOUBLEPRIME_771  = 15
    NER_CMB_TOPDDOUBLEPRIME  = 1
    NER_OUTER_CORE           = 16
    NER_TOP_CENTRAL_CUBE_ICB = 2

    ! Earth minimum period example:
    !   deg2km ~ 111km, vs_min ~ 2.25 km/s
    !   90./NEX * 111.0 / (NGLL-1) * 4 / 2.25 -> NEX = 80, NGLL = 5: T_min ~ 55.5

    ! sets empirical values for time step size, attenuation range (for 3 SLS) and number of element layers
    if (NEX_MAX*multiplication_factor <= 80) then
      ! time step
      DT                       = 0.252d0
      ! attenuation period range
      MIN_ATTENUATION_PERIOD   = 30.d0
      MAX_ATTENUATION_PERIOD   = 1500.d0
      ! radius of central cube
      R_CENTRAL_CUBE = 950000.d0

    else if (NEX_MAX*multiplication_factor <= 96) then
      ! time step
      ! to handle a case that Zhinan Xie found to be unstable for NEX = 96 I reduce the time step to 90% of its value here
      DT                       = 0.252d0 * 0.90d0
      ! attenuation period range
      MIN_ATTENUATION_PERIOD   = 30.d0
      MAX_ATTENUATION_PERIOD   = 1500.d0
      ! radius of central cube
      R_CENTRAL_CUBE = 950000.d0

    ! element width =   0.5625000      degrees =    62.54715      km
    else if (NEX_MAX*multiplication_factor <= 160) then
      ! time step
      DT                       = 0.252d0
      ! attenuation period range
      MIN_ATTENUATION_PERIOD   = 30.d0
      MAX_ATTENUATION_PERIOD   = 1500.d0
      ! radius of central cube
      R_CENTRAL_CUBE = 950000.d0

    ! element width =   0.3515625      degrees =    39.09196      km
    else if (NEX_MAX*multiplication_factor <= 256) then
      DT                       = 0.225d0

      MIN_ATTENUATION_PERIOD   = 20.d0
      MAX_ATTENUATION_PERIOD   = 1000.d0

      ! adding more layering
      NER_400_220              = 3
      NER_600_400              = 3
      NER_TOPDDOUBLEPRIME_771  = 22
      NER_CMB_TOPDDOUBLEPRIME  = 2
      NER_OUTER_CORE           = 24
      NER_TOP_CENTRAL_CUBE_ICB = 3

      R_CENTRAL_CUBE = 965000.d0

    ! element width =   0.2812500      degrees =    31.27357      km
    else if (NEX_MAX*multiplication_factor <= 320) then
      DT                       = 0.16d0

      MIN_ATTENUATION_PERIOD   = 15.d0
      MAX_ATTENUATION_PERIOD   = 750.d0

      ! adding more layering
      NER_220_80               = 3
      NER_400_220              = 4
      NER_600_400              = 4
      NER_771_670              = 2
      NER_TOPDDOUBLEPRIME_771  = 29
      NER_CMB_TOPDDOUBLEPRIME  = 2
      NER_OUTER_CORE           = 32
      NER_TOP_CENTRAL_CUBE_ICB = 4

      R_CENTRAL_CUBE = 940000.d0

    ! element width =   0.1875000      degrees =    20.84905      km
    else if (NEX_MAX*multiplication_factor <= 480) then
      DT                       = 0.11d0

      MIN_ATTENUATION_PERIOD   = 10.d0
      MAX_ATTENUATION_PERIOD   = 500.d0

      ! adding more layering
      NER_80_MOHO              = 2
      NER_220_80               = 4
      NER_400_220              = 5
      NER_600_400              = 6
      NER_670_600              = 2
      NER_771_670              = 2
      NER_TOPDDOUBLEPRIME_771  = 44
      NER_CMB_TOPDDOUBLEPRIME  = 3
      NER_OUTER_CORE           = 48
      NER_TOP_CENTRAL_CUBE_ICB = 5

      R_CENTRAL_CUBE = 988000.d0

    ! element width =   0.1757812      degrees =    19.54598      km
    else if (NEX_MAX*multiplication_factor <= 512) then
      DT                       = 0.1125d0

      MIN_ATTENUATION_PERIOD   = 9.d0
      MAX_ATTENUATION_PERIOD   = 500.d0

      ! adding more layering
      NER_80_MOHO              = 2
      NER_220_80               = 4
      NER_400_220              = 6
      NER_600_400              = 6
      NER_670_600              = 2
      NER_771_670              = 3
      NER_TOPDDOUBLEPRIME_771  = 47
      NER_CMB_TOPDDOUBLEPRIME  = 3
      NER_OUTER_CORE           = 51
      NER_TOP_CENTRAL_CUBE_ICB = 5

      R_CENTRAL_CUBE = 1010000.d0

    ! element width =   0.1406250      degrees =    15.63679      km
    else if (NEX_MAX*multiplication_factor <= 640) then
      DT                       = 0.09d0

      MIN_ATTENUATION_PERIOD   = 8.d0
      MAX_ATTENUATION_PERIOD   = 400.d0

      ! adding more layering
      NER_CRUST                = 2
      NER_80_MOHO              = 3
      NER_220_80               = 5
      NER_400_220              = 7
      NER_600_400              = 8
      NER_670_600              = 3
      NER_771_670              = 3
      NER_TOPDDOUBLEPRIME_771  = 59
      NER_CMB_TOPDDOUBLEPRIME  = 4
      NER_OUTER_CORE           = 64
      NER_TOP_CENTRAL_CUBE_ICB = 6

      R_CENTRAL_CUBE = 1020000.d0

    ! element width =   0.1041667      degrees =    11.58280      km
    else if (NEX_MAX*multiplication_factor <= 864) then
      DT                       = 0.0667d0

      MIN_ATTENUATION_PERIOD   = 6.d0
      MAX_ATTENUATION_PERIOD   = 300.d0

      ! adding more layering
      NER_CRUST                = 2
      NER_80_MOHO              = 4
      NER_220_80               = 6
      NER_400_220              = 10
      NER_600_400              = 10
      NER_670_600              = 3
      NER_771_670              = 4
      NER_TOPDDOUBLEPRIME_771  = 79
      NER_CMB_TOPDDOUBLEPRIME  = 5
      NER_OUTER_CORE           = 86
      NER_TOP_CENTRAL_CUBE_ICB = 9

      R_CENTRAL_CUBE = 990000.d0

    ! element width =   7.8125000E-02  degrees =    8.687103      km
    else if (NEX_MAX*multiplication_factor <= 1152) then
      DT                       = 0.05d0

      MIN_ATTENUATION_PERIOD   = 4.d0
      MAX_ATTENUATION_PERIOD   = 200.d0

      ! adding more layering
      NER_CRUST                = 3
      NER_80_MOHO              = 6
      NER_220_80               = 8
      NER_400_220              = 13
      NER_600_400              = 13
      NER_670_600              = 4
      NER_771_670              = 6
      NER_TOPDDOUBLEPRIME_771  = 106
      NER_CMB_TOPDDOUBLEPRIME  = 7
      NER_OUTER_CORE           = 116
      NER_TOP_CENTRAL_CUBE_ICB = 12

      R_CENTRAL_CUBE = 985000.d0

    ! element width =   7.2115384E-02  degrees =    8.018865      km
    else if (NEX_MAX*multiplication_factor <= 1248) then
      DT                       = 0.0462d0

      MIN_ATTENUATION_PERIOD   = 4.d0
      MAX_ATTENUATION_PERIOD   = 200.d0

      ! adding more layering
      NER_CRUST                = 3
      NER_80_MOHO              = 6
      NER_220_80               = 9
      NER_400_220              = 14
      NER_600_400              = 14
      NER_670_600              = 5
      NER_771_670              = 6
      NER_TOPDDOUBLEPRIME_771  = 114
      NER_CMB_TOPDDOUBLEPRIME  = 8
      NER_OUTER_CORE           = 124
      NER_TOP_CENTRAL_CUBE_ICB = 13

      R_CENTRAL_CUBE = 985000.d0

    else
      ! scale with respect to 1248 if above that limit
      DT                       = 0.0462d0 * 1248.d0 / (2.d0*NEX_MAX)

      MIN_ATTENUATION_PERIOD   = 4.d0
      MAX_ATTENUATION_PERIOD   = 200.d0

      ! adding more layering
      NER_CRUST                = nint(3 * 2.d0*NEX_MAX / 1248.d0)
      NER_80_MOHO              = nint(6 * 2.d0*NEX_MAX / 1248.d0)
      NER_220_80               = nint(9 * 2.d0*NEX_MAX / 1248.d0)
      NER_400_220              = nint(14 * 2.d0*NEX_MAX / 1248.d0)
      NER_600_400              = nint(14 * 2.d0*NEX_MAX / 1248.d0)
      NER_670_600              = nint(5 * 2.d0*NEX_MAX / 1248.d0)
      NER_771_670              = nint(6 * 2.d0*NEX_MAX / 1248.d0)
      NER_TOPDDOUBLEPRIME_771  = nint(114 * 2.d0*NEX_MAX / 1248.d0)
      NER_CMB_TOPDDOUBLEPRIME  = nint(8 * 2.d0*NEX_MAX / 1248.d0)
      NER_OUTER_CORE           = nint(124 * 2.d0*NEX_MAX / 1248.d0)
      NER_TOP_CENTRAL_CUBE_ICB = nint(13 * 2.d0*NEX_MAX / 1248.d0)

      R_CENTRAL_CUBE = 985000.d0
    endif
    ! uses automatic estimates auto_ner() for NEX > 1248
    nex_max_auto_ner_estimate = 1248

    !> Hejun
    ! avoids elongated elements below the 670-discontinuity,
    ! since for model REFERENCE_MODEL_1DREF,
    ! the 670-discontinuity is moved up to 650 km depth.
    if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1DREF) then
      NER_771_670 = NER_771_670 + 1
    endif

    ! case of regular PREM with two crustal layers: change the time step for small meshes
    ! because of a different size of elements in the radial direction in the crust
    ! time step
    if (HONOR_1D_SPHERICAL_MOHO) then
      ! 1D models honor 1D spherical moho
      if (.not. ONE_CRUST) then
        ! makes time step smaller
        if (NEX_MAX*multiplication_factor <= 160) then
          DT = 0.20d0
        else if (NEX_MAX*multiplication_factor <= 256) then
          DT = 0.20d0
        endif
      endif
    else
      ! makes time step smaller
      if (NEX_MAX*multiplication_factor <= 80) then
          DT = 0.125d0
      else if (NEX_MAX*multiplication_factor <= 160) then
          DT = 0.15d0
      else if (NEX_MAX*multiplication_factor <= 256) then
          DT = 0.17d0
      else if (NEX_MAX*multiplication_factor <= 320) then
          DT = 0.155d0
      endif
    endif

  case default
    ! avoiding exit_MPI(), since we also call this routine in create_header_file
    ! which can be compiled without MPI - using stop instead
    !call exit_MPI(myrank,'Invalid planet, timestep and layers not implemented yet')
    print *,'Invalid planet, timestep and layers not implemented yet'
    stop 'Invalid planet, timestep and layers not implemented yet'
  end select

  ! in case we stretch the element layers to account for moho, we need at least 2 element layers in the crust
  if (HONOR_1D_SPHERICAL_MOHO) then
    ! 1D models honor 1D spherical moho
    if (.not. ONE_CRUST) then
      ! case 1D + two crustal layers
      !
      ! note: we use here 2 element layers for the crust to honor also the middle crust
      !       for example, PREM distinguishes different constant velocities for upper crust (vp=5.8km/s, depth down to 15km)
      !       and lower crust (vp=6.8km/s, depth down to 24.4km)
      !
      !       be aware that using 2 element layers in the crust will lead to very thin elements for the lower crust
      !       which decreases significantly the time step size for stability
      if (NER_CRUST < 2 ) NER_CRUST = 2
    endif
  else
    ! 3D models: must have two element layers for crust
    if (NER_CRUST < 2 ) NER_CRUST = 2
  endif

  !----
  !----  change some values in the case of regular PREM with two crustal layers or of 3D models
  !----

  ! minimum width of chunk
  min_chunk_width_in_degrees = min(ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES)

  ! gets attenuation min/max range
  if (.not. ATTENUATION_RANGE_PREDEFINED) then
     call auto_attenuation_periods(min_chunk_width_in_degrees, NEX_MAX, MIN_ATTENUATION_PERIOD, MAX_ATTENUATION_PERIOD)
  endif

  ! note: global estimate above for DT is empirical for chunk sizes of 90 degrees.
  ! if regional chunk size is larger, the time step is still limited by the vertical size of elements,
  ! which depends only on the number of vertical element layers.
  ! this however gets estimated above for 90 degrees already, so time stepping remains unchanged.

  ! adapts number of layer elements and time step size
  ! (for regional simulations with chunk sizes < 90 degrees)
  ! (or very, very large meshes)
  if (min_chunk_width_in_degrees < 90.0d0 .or. NEX_MAX > nex_max_auto_ner_estimate) then
    ! adapts number of layer elements and time step size

    ! note: for global simulations, we set
    !         ANGULAR_WIDTH_XI_IN_DEGREES = 90.0d0 and
    !         ANGULAR_WIDTH_ETA_IN_DEGREES = 90.0d0 in read_parameter_file.f90

    ! gets number of element-layers
    call auto_ner(min_chunk_width_in_degrees, NEX_MAX)

    ! re-sets attenuation min/max range
    call auto_attenuation_periods(min_chunk_width_in_degrees, NEX_MAX, MIN_ATTENUATION_PERIOD, MAX_ATTENUATION_PERIOD)

    ! gets time step size
    call auto_time_stepping(min_chunk_width_in_degrees, NEX_MAX, dt_auto)

    !debug
    !if (myrank == 0) print *,'debug: get_timestep: min width',min_chunk_width_in_degrees,NEX_MAX,'dt_auto',dt_auto,'DT',DT

    ! note: automatic time step might overestimate the time step size for chunk sizes larger ~ 40 degrees
    !       thus we only replace the empirical time step size if the estimate gets smaller

    ! sets new time step size
    if (dt_auto < DT) DT = dt_auto

    ! checks minimum number of element-layers in crust
    if (HONOR_1D_SPHERICAL_MOHO) then
      if (.not. ONE_CRUST) then
        ! case 1D + two crustal layers
        if (NER_CRUST < 2 ) NER_CRUST = 2
      endif
    else
      ! case 3D
      if (NER_CRUST < 2 ) NER_CRUST = 2
    endif

    ! Mars & Moon
    ! sets minimum number of element layers in crust
    if (PLANET_TYPE == IPLANET_MARS .or. PLANET_TYPE == IPLANET_MOON) then
      if (HONOR_1D_SPHERICAL_MOHO) then
        if (.not. ONE_CRUST) then
          ! case 1D + two+ crustal layers
          if (NER_CRUST < 3 ) NER_CRUST = 3
        endif
      else
        ! case 3D
        if (NER_CRUST < 5 ) NER_CRUST = 5
      endif
    endif

  else
    ! adapts the empirical time step estimate to different NGLL settings.
    ! the empirical sizes are estimated for NGLL == 5, here we modify them according to the change
    ! of the minimum spacing between different NGLL settings.
    if (NGLLX /= 5) then
      ! relative minimum distance between two GLL points
      ! the roots x_i are given by the first derivative of the Legendre Polynomial: P_n-1'(x_i) = 0
      !
      ! note: the x_i interval is between [-1,1], thus relative to the full length, we divide by 2
      !
      ! formulas:
      ! see: https://en.wikipedia.org/wiki/Gaussian_quadrature  -> section Gauss-Lobatto rules
      !      http://mathworld.wolfram.com/LobattoQuadrature.html
      !
      ! numerical values:
      ! see: http://keisan.casio.com/exec/system/1280801905

      ! spacing for NGLLX == 5
      MIN_GLL_POINT_SPACING_NGLL5 = 0.5d0 * ( 1.d0 - sqrt(3.d0 / 7.d0) ) ! 0.1726

      ! NGLL choosen
      ! (see also in auto_ner.f90)
      select case (NGLLX)
      case (2)
        MIN_GLL_POINT_SPACING = 0.5d0 * ( 1.d0 - 1.0 ) ! 1.0
      case (3)
        MIN_GLL_POINT_SPACING = 0.5d0 * ( 1.d0 - 0.0 ) ! 0.5
      case (4)
        MIN_GLL_POINT_SPACING = 0.5d0 * ( 1.d0 - sqrt(1.d0 / 5.d0) ) ! 0.2764
      case (5)
        MIN_GLL_POINT_SPACING = 0.5d0 * ( 1.d0 - sqrt(3.d0 / 7.d0) ) ! 0.1726
      case (6)
        MIN_GLL_POINT_SPACING = 0.5d0 * ( 1.d0 - sqrt(1.d0/21.d0*(7.d0 + 2.d0 * sqrt(7.d0))) ) !0.117472
      case (7)
        MIN_GLL_POINT_SPACING = 0.5d0 * ( 1.d0 - 0.830223896278566929872 ) ! 0.084888
      case (8)
        MIN_GLL_POINT_SPACING = 0.5d0 * ( 1.d0 - 0.8717401485096066153374 ) ! 0.0641299
      case (9)
        MIN_GLL_POINT_SPACING = 0.5d0 * ( 1.d0 - 0.8997579954114601573123 ) ! 0.050121
      case (10)
        MIN_GLL_POINT_SPACING = 0.5d0 * ( 1.d0 - 0.9195339081664588138289 ) ! 0.040233
      case default
        ! no formula yet, takes average
        MIN_GLL_POINT_SPACING = 0.5d0 * ( 2.d0 / dble(NGLLX-1) )
        !stop 'get_timestep_and_layers: NGLLX > 10 value not supported yet! please consider adding it...'
      end select

      ! adapts DT setting
      DT = DT * MIN_GLL_POINT_SPACING / MIN_GLL_POINT_SPACING_NGLL5
    endif
  endif

  ! cut-off mesh
  if (REGIONAL_MESH_CUTOFF) then
    ! sets number of element layers to zero below the cut-off depth
    if (REGIONAL_MESH_CUTOFF_DEPTH <= 771.d0) then
      NER_TOPDDOUBLEPRIME_771  = 0
      NER_CMB_TOPDDOUBLEPRIME  = 0
      NER_OUTER_CORE           = 0
      NER_TOP_CENTRAL_CUBE_ICB = 0
    endif
    if (REGIONAL_MESH_CUTOFF_DEPTH <= 670.d0) NER_771_670  = 0
    if (REGIONAL_MESH_CUTOFF_DEPTH <= 600.d0) NER_670_600  = 0
    if (REGIONAL_MESH_CUTOFF_DEPTH <= 400.d0) NER_600_400  = 0
    if (REGIONAL_MESH_CUTOFF_DEPTH <= 220.d0) NER_400_220  = 0
    if (REGIONAL_MESH_CUTOFF_DEPTH <= 80.d0) NER_220_80   = 0
    if (REGIONAL_MESH_CUTOFF_DEPTH <= 24.4d0) NER_80_MOHO = 0

    if (REGIONAL_MESH_CUTOFF_DEPTH < 24.0d0) then
      print *,'Regional mesh cutoff depth ',REGIONAL_MESH_CUTOFF_DEPTH,' is too shallow. Please select a depth >= 24.4 km.'
      stop 'Invalid regional mesh cutoff depth'
    endif
  endif

!---
!
! ADD YOUR MODEL HERE
!
!---


  ! time step reductions are based on empirical values (..somehow)
  select case (PLANET_TYPE)
  case (IPLANET_MARS)
    ! Mars
    if (NCHUNKS == 6) then
      if (CRUSTAL .and. CASE_3D) then
        ! limits time step size for thinner crustal elements in mars_1D model
        if (REFERENCE_1D_MODEL == REFERENCE_MODEL_MARS_1D) then
          if (DT > 0.06) DT = 0.06
        endif
        ! crustmaps safety
        if (REFERENCE_CRUSTAL_MODEL == ICRUST_CRUSTMAPS) then
          DT = DT * (1.d0 - 0.1d0)
        endif
      else
        ! Mars 1D model time step reduction
        if (REFERENCE_1D_MODEL == REFERENCE_MODEL_MARS_1D) then
          DT = DT * (1.d0 - 0.3d0)
        endif
      endif
    endif
    if (MARS_REGIONAL_MOHO_MESH .and. MARS_HONOR_DEEP_MOHO) then
      if (HONOR_1D_SPHERICAL_MOHO) then
        ! spherical moho depth, and nothing to deform, default layering will be sufficient
        continue
      else
        ! enforces 5 element layers for regional simulation
        NER_CRUST = 5
      endif
    endif
    ! takes a 5% safety margin on the maximum stable time step
    ! which was obtained by trial and error
    DT = DT * (1.d0 - 0.05d0)

  case (IPLANET_MOON)
    ! takes a 5% safety margin on the maximum stable time step
    ! which was obtained by trial and error
    DT = DT * (1.d0 - 0.05d0)

  case (IPLANET_EARTH)
    ! Earth
    ! following models need special attention, at least for global simulations:
    if (NCHUNKS == 6) then
      ! makes time step smaller for this ref model, otherwise becomes unstable in fluid
      if (REFERENCE_1D_MODEL == REFERENCE_MODEL_IASP91) &
        DT = DT*(1.d0 - 0.3d0)

      ! using inner core anisotropy
      if (ANISOTROPIC_INNER_CORE) then
        ! note: main limiting time step constraint is usually in the crust/mantle elements for 3D models
        !       for details, you can check in the output_mesher.txt
        continue
      endif

      ! makes time step smaller for certain crustal models, otherwise becomes unstable in solid
      ! CRUSTAL: indicates a 3D crustal model, like CRUST2.0 will be used
      ! CASE_3D: indicates element stretching to honor e.g. moho depths and/or upper/lower crusts
      if (CRUSTAL .and. CASE_3D) then
        ! reduces time step size for CRUST1.0 crustal model
        if (REFERENCE_CRUSTAL_MODEL == ICRUST_CRUST1) &
          DT = DT*(1.d0 - 0.1d0)
        ! reduces time step size for crustmaps crustal model
        if (REFERENCE_CRUSTAL_MODEL == ICRUST_CRUSTMAPS) &
          DT = DT*(1.d0 - 0.3d0)
        ! reduces time step size for SGLOBE-rani crustal model
        if (REFERENCE_CRUSTAL_MODEL == ICRUST_SGLOBECRUST) &
          DT = DT*(1.d0 - 0.4d0)
        ! reduces time step size for SPiRaL crustal model
        ! (this is only needed for meshes with NEX < 144 due to critical element shapes in the crust)
        if (REFERENCE_CRUSTAL_MODEL == ICRUST_SPIRAL .and. NEX_MAX < 144) &
          DT = DT*(1.d0 - 0.4d0)
      endif
    endif

    ! following models need special attention, regardless of number of chunks:
    ! it makes the time step smaller for this ref model, otherwise becomes unstable in fluid
    if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1066A) &
      DT = DT*(1.d0 - 0.8d0)  ! *0.20d0

    ! reduces time step size for "no mud" version of AK135F model
    if (REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135F_NO_MUD) &
      DT = DT*(1.d0 - 0.05d0)

    !  decreases time step as otherwise the solution might become unstable for rougher/unsmoothed models
    if (.false.) then
      if (THREE_D_MODEL == THREE_D_MODEL_PPM ) DT = DT * (1.d0 - 0.2d0)
    endif

    ! takes a 5% safety margin on the maximum stable time step
    ! which was obtained by trial and error
    DT = DT * (1.d0 - 0.05d0)

    ! adapts number of element layers in crust and time step for regional simulations
    if (REGIONAL_MOHO_MESH) then
      ! hard coded number of crustal element layers and time step

      ! checks
      if (NCHUNKS > 1 ) stop 'regional moho mesh: NCHUNKS error in rcp_set_timestep_and_layers'

      ! increases number of layers due to element deformations when honoring large moho depth variations (7km - 70km).
      ! this should lead to a better mesh quality.
      if (HONOR_1D_SPHERICAL_MOHO) then
        ! spherical moho depth, and nothing to deform, default layering will be sufficient
        continue
      else
        ! original values
        !print *,'NER:',NER_CRUST
        !print *,'DT:',DT

        ! enforce 3 element layers
        NER_CRUST = 3

        ! increased stability, empirical
        DT = DT*(1.d0 + 0.5d0)

        ! empirical values for different regions
        if (EARTH_REGIONAL_MOHO_MESH_EUROPE ) DT = 0.17 ! Europe
        if (EARTH_REGIONAL_MOHO_MESH_ASIA ) DT = 0.15 ! Asia & Middle East
      endif
    endif

  case default
    ! avoiding exit_MPI(), since we also call this routine in create_header_file
    ! which can be compiled without MPI - using stop instead
    !call exit_MPI(myrank,'Invalid planet, timestep and layers not implemented yet')
    print *,'Invalid planet, timestep and layers not implemented yet'
    stop 'Invalid planet, timestep and layers not implemented yet'
  end select ! planet_type

  ! Further reduce time step for FULL_GRAVITY
  if (FULL_GRAVITY) then
    if (THREE_D_MODEL > 0) DT = 0.7d0 * DT    ! 0.7 is an arbitrary value
  endif

  ! the maximum CFL of LDDRK is significantly higher than that of the Newmark scheme,
  ! in a ratio that is theoretically 1.327 / 0.697 = 1.15 / 0.604 = 1.903 for a solid with Poisson's ratio = 0.25
  ! and for a fluid (see the manual of the 2D code, SPECFEM2D, Tables 4.1 and 4.2, and that ratio does not
  ! depend on whether we are in 2D or in 3D). However in practice a ratio of about 1.5 to 1.7 is often safer
  ! (for instance for models with a large range of Poisson's ratio values).
  ! Since the code computes the time step using the Newmark scheme, for LDDRK we simply
  ! multiply that time step by this ratio when LDDRK is on and when flag INCREASE_CFL_FOR_LDDRK is true.
  if (USE_LDDRK .and. INCREASE_CFL_FOR_LDDRK) DT = DT * RATIO_BY_WHICH_TO_INCREASE_IT

  ! cut at a significant number of digits (e.g., 2 digits)
  ! in steps of 1/2 digits
  ! example: 0.0734815 -> 0.0730
  !          0.07371   -> 0.0735
  call get_timestep_limit_significant_digit(DT)

  ! overrides DT in case specified in Par_file
  if (USER_DT > 0.d0) then
    ! overrides DT
    if (myrank == 0) then
      print *,'simulation time step size:'
      print *,'  DT determined = ',DT
      print *,'  Par_file: user overrides with specified DT = ',USER_DT
      print *
    endif
    DT = USER_DT
  endif

  end subroutine get_timestep_and_layers

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_timestep_limit_significant_digit(time_step)

  ! cut at a significant number of digits (e.g., 2 digits) using 1/2 rounding
  ! example: 0.0734815 -> 0.0730
  !      and 0.0737777 -> 0.0735
  !
  ! also works with different magnitudes of time step sizes (0.118, 0.00523, ..). always cut of after 2 significant digits:
  ! example: 0.118749999 -> 0.115

  implicit none

  double precision,intent(inout) :: time_step

  ! rounding
  integer :: lpow,ival
  double precision :: fac_pow,dt_cut

  ! initializes
  dt_cut = time_step

  ! cut at a significant number of digits (2 digits)
  ! example: 0.0734815 -> lpow = (2 - (-1) = 3
  lpow = int(2.d0 - log10(dt_cut))

  ! example: -> factor 10**3
  fac_pow = 10.d0**(lpow)

  ! example: -> 73
  ival = int(fac_pow * dt_cut)

  ! adds .5-digit (in case): 73.0 -> 0.073
  if ( (fac_pow * dt_cut - ival) >= 0.5 ) then
    dt_cut = (dble(ival) + 0.5d0) / fac_pow
  else
    dt_cut = dble(ival) / fac_pow
  endif

  time_step = dt_cut

  end subroutine get_timestep_limit_significant_digit

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_minimum_period_estimate()

  use constants, only: NGLLX,PI,NPTS_PER_WAVELENGTH,REFERENCE_MODEL_CASE65TAY

  use shared_parameters, only: T_min_period,estimated_min_wavelength, &
    ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES, &
    NEX_XI,NEX_ETA, &
    PLANET_TYPE,IPLANET_EARTH,IPLANET_MARS,IPLANET_MOON,R_PLANET, &
    REFERENCE_1D_MODEL

  implicit none

  ! local parameters
  integer :: NEX_MAX
  double precision :: tmp,deg2km
  double precision :: width,gll_spacing
  double precision :: S_VELOCITY_MIN

  ! we often use an estimate based on NGLL == 5, assuming that the number of points per wavelength
  ! coincides with the number of GLL points and thus the element size is the same length a the minimum wavelength:
  !
  !   Earth: 2 * PI * 6371km / 4 / 256 / 2.3 km/s ~ 17 s
  !
  !   Mars : 2 * PI * 3390km / 4 / 256 / 4.0 km/s ~ 5 s
  !
  !   Moon : 2 * PI * 1737.1km / 4 /256 / 1.8 km/s ~ 6 s
  !
  !   note: for mars, the current mesh leads to a higher resolution in the crust with the above expected estimation,
  !         but lower resolution in the upper mantle region (below the doubling layer) with about twice the estimation size.
  !
  !   daniel todo: we might want to re-evaluate these mesh and minimum period estimations...
  !                the mesh resolution estimations mostly affects the simulation by
  !                setting the attenuation min/max period band in auto_ner.f90
  !
  !                for now, these minimum period estimations are mostly meant to inform users at which periods
  !                the seismograms should be filtered above.
  !
  !select case (PLANET_TYPE)
  !case (IPLANET_EARTH)
  !  ! Earth
  !  T_min_res = 17.d0
  !case (IPLANET_MARS)
  !  ! Mars
  !  T_min_res = 5.d0
  !case (IPLANET_MOON)
  !  ! Moon
  !  T_min_res = 6.d0
  !case default
  !  ! avoiding exit_MPI(), since we also call this routine in create_header_file
  !  ! which can be compiled without MPI - using stop instead
  !  !call exit_MPI(myrank,'Invalid planet, timestep and layers not implemented yet')
  !  print *,'Invalid planet, timestep and layers not implemented yet'
  !  stop 'Invalid planet, timestep and layers not implemented yet'
  !end select
  !T_min = max(ANGULAR_WIDTH_ETA_IN_DEGREES,ANGULAR_WIDTH_XI_IN_DEGREES)/90.0 * 256.0/min(NEX_ETA,NEX_XI) * T_min_res
  !
  !
  ! Here, we use a slighty more general estimation which allows for different NGLL settings.
  ! Also, we estimate the lower bound of the minimum period resolved by the model mesh
  ! and not the upper bound of the minimum period as in the above formula.
  !
  ! This lower bound on the minimum period (T_min_period) will also be used for the attenuation period range estimation.
  !
  ! minimum period estimation
  select case (PLANET_TYPE)
  case (IPLANET_EARTH)
    ! Earth
    ! minimum velocity (Vs)
    S_VELOCITY_MIN = 2.25d0

  case (IPLANET_MARS)
    ! Mars
    ! minimum velocity (Vs)
    if (REFERENCE_1D_MODEL == REFERENCE_MODEL_CASE65TAY) then
      ! based on case65tay, minimum ~ 2.48 km/s
      S_VELOCITY_MIN = 2.48d0
    else
      ! based on Sohl & Spohn, minimum Vs is ~ 4 km/s
      S_VELOCITY_MIN = 4.d0
    endif

  case (IPLANET_MOON)
    ! Moon
    ! VPREMOON has very slow Vs close to surface (0.5 km/s) (regolith layer)
    ! for this estimate, we use a slightly larger value here (vs corresponds to 10km depth)
    S_VELOCITY_MIN = 1.8d0

  case default
    ! avoiding exit_MPI(), since we also call this routine in create_header_file
    ! which can be compiled without MPI - using stop instead
    print *,'Invalid planet in get_minimum_period_estimate() not implemented yet'
    stop 'Invalid planet in get_minimum_period_estimate() not implemented yet'
  end select

  ! number of elements along one side of a chunk
  NEX_MAX = max(NEX_XI,NEX_ETA)

  ! width of chunk
  width = min(ANGULAR_WIDTH_ETA_IN_DEGREES,ANGULAR_WIDTH_XI_IN_DEGREES)

  ! average spacing between GLL points
  gll_spacing = dble(NGLLX - 1)

  ! check
  if (NEX_MAX <= 0) &
    stop 'Invalid NEX_MAX in get_minimum_period_estimate()'
  if (width <= 0.d0) &
    stop 'Invalid WIDTH in get_minimum_period_estimate()'
  if (gll_spacing <= 0.d0) &
    stop 'Invalid gll_spacing in get_minimum_period_estimate()'

  ! converts degree to km
  ! example Earth: radius 6371 km -> 2 * PI * R / 360 ~ 111.19 km
  deg2km = R_PLANET / 1000.d0 * 2.d0 * PI / 360.d0

  ! computes Min Period
  !
  ! width of element in km = (Angular width in degrees / NEX_MAX) * degrees to km
  tmp = width * deg2km / dble(NEX_MAX)

  ! average grid node spacing in km = Width of an element in km / spacing for GLL point
  tmp = tmp / gll_spacing

  ! minimum resolved wavelength (for fixed number of points per wavelength)
  tmp = tmp * dble(NPTS_PER_WAVELENGTH-1)

  ! minimum period resolved = (minimum wavelength) / V_min
  tmp = tmp/S_VELOCITY_MIN

  ! estimated minimum period resolved
  T_min_period = tmp

  ! estimated minimum wavelength
  estimated_min_wavelength = T_min_period * S_VELOCITY_MIN

  end subroutine get_minimum_period_estimate

!
!-------------------------------------------------------------------------------------------------
!

  subroutine band_instrument_code(DT,bic)

! This subroutine is to choose the appropriate band and instrument codes for channel names of seismograms
! based on the IRIS convention (first two letters of channel codes which were LH(Z/E/N) previously).
! For consistency with observed data, we now use the IRIS convention for band codes (first letter in channel codes)of
! SEM seismograms governed by their sampling rate.
! Instrument code (second letter in channel codes) is fixed to "X" which is assigned by IRIS for synthetic seismograms.
! See the manual for further explanations!
! Ebru, November 2010

  implicit none

  double precision,intent(in) :: DT
  character(len=2),intent(out) :: bic

  bic = ''

  if (1.0d0 <= DT)  bic = 'LX'
  if (0.1d0 < DT .and. DT < 1.0d0) bic = 'MX'
  if (0.0125d0 < DT .and. DT <= 0.1d0) bic = 'BX'
  if (0.004d0 < DT .and. DT <= 0.0125d0) bic = 'HX'
  if (0.001d0 < DT .and. DT <= 0.004d0) bic = 'CX'
  if (DT <= 0.001d0) bic = 'FX'

  end subroutine band_instrument_code
