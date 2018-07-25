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


  subroutine get_timestep_and_layers()

  use constants
  use shared_parameters

  implicit none

  ! local variables
  integer :: NEX_MAX
  integer :: multiplication_factor
  double precision :: min_chunk_width_in_degrees
  double precision :: dt_auto

  !----
  !----  case prem_onecrust by default
  !----

  NEX_MAX = max(NEX_XI,NEX_ETA)

  ! to suppress the crustal layers
  ! (replaced by an extension of the mantle: R_EARTH is not modified, but no more crustal doubling)
  if (SUPPRESS_CRUSTAL_MESH) then
    multiplication_factor = 2
  else
    multiplication_factor = 1
  endif

  ! sets empirical values for time step size, attenuation range (for 3 SLS) and number of element layers
  if (NEX_MAX*multiplication_factor <= 80) then
    ! time step
    DT                       = 0.252d0

    ! attenuation period range
    MIN_ATTENUATION_PERIOD   = 30
    MAX_ATTENUATION_PERIOD   = 1500

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

    ! radius of central cube
    R_CENTRAL_CUBE = 950000.d0

  else if (NEX_MAX*multiplication_factor <= 96) then
    ! time step
!! DK DK to handle a case that Zhinan Xie found to be unstable for NEX = 96 I reduce the time step to 90% of its value here
    DT                       = 0.252d0 * 0.90d0

    ! attenuation period range
    MIN_ATTENUATION_PERIOD   = 30
    MAX_ATTENUATION_PERIOD   = 1500

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

    ! radius of central cube
    R_CENTRAL_CUBE = 950000.d0

  ! element width =   0.5625000      degrees =    62.54715      km
  else if (NEX_MAX*multiplication_factor <= 160) then
    ! time step
    DT                       = 0.252d0

    ! attenuation period range
    MIN_ATTENUATION_PERIOD   = 30
    MAX_ATTENUATION_PERIOD   = 1500

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

    ! radius of central cube
    R_CENTRAL_CUBE = 950000.d0

  ! element width =   0.3515625      degrees =    39.09196      km
  else if (NEX_MAX*multiplication_factor <= 256) then
    DT                       = 0.225d0

    MIN_ATTENUATION_PERIOD   = 20
    MAX_ATTENUATION_PERIOD   = 1000

    NER_CRUST                = 1
    NER_80_MOHO              = 1
    NER_220_80               = 2
    NER_400_220              = 3
    NER_600_400              = 3
    NER_670_600              = 1
    NER_771_670              = 1
    NER_TOPDDOUBLEPRIME_771  = 22
    NER_CMB_TOPDDOUBLEPRIME  = 2
    NER_OUTER_CORE           = 24
    NER_TOP_CENTRAL_CUBE_ICB = 3
    R_CENTRAL_CUBE = 965000.d0

  ! element width =   0.2812500      degrees =    31.27357      km
  else if (NEX_MAX*multiplication_factor <= 320) then
    DT                       = 0.16d0

    MIN_ATTENUATION_PERIOD   = 15
    MAX_ATTENUATION_PERIOD   = 750

    NER_CRUST                = 1
    NER_80_MOHO              = 1
    NER_220_80               = 3
    NER_400_220              = 4
    NER_600_400              = 4
    NER_670_600              = 1
    NER_771_670              = 2
    NER_TOPDDOUBLEPRIME_771  = 29
    NER_CMB_TOPDDOUBLEPRIME  = 2
    NER_OUTER_CORE           = 32
    NER_TOP_CENTRAL_CUBE_ICB = 4
    R_CENTRAL_CUBE = 940000.d0

  ! element width =   0.1875000      degrees =    20.84905      km
  else if (NEX_MAX*multiplication_factor <= 480) then
    DT                       = 0.11d0

    MIN_ATTENUATION_PERIOD   = 10
    MAX_ATTENUATION_PERIOD   = 500

    NER_CRUST                = 1
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

    MIN_ATTENUATION_PERIOD   = 9
    MAX_ATTENUATION_PERIOD   = 500

    NER_CRUST                = 1
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

    MIN_ATTENUATION_PERIOD   = 8
    MAX_ATTENUATION_PERIOD   = 400

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

    MIN_ATTENUATION_PERIOD   = 6
    MAX_ATTENUATION_PERIOD   = 300

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

    MIN_ATTENUATION_PERIOD   = 4
    MAX_ATTENUATION_PERIOD   = 200

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

    MIN_ATTENUATION_PERIOD   = 4
    MAX_ATTENUATION_PERIOD   = 200

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

    MIN_ATTENUATION_PERIOD   = 4
    MAX_ATTENUATION_PERIOD   = 200

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

  !! removed this limit           else
  !! removed this limit             stop 'problem with this value of NEX_MAX'
  endif

  !> Hejun
  ! avoids elongated elements below the 670-discontinuity,
  ! since for model REFERENCE_MODEL_1DREF,
  ! the 670-discontinuity is moved up to 650 km depth.
  if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1DREF) then
    NER_771_670 = NER_771_670 + 1
  endif

  !----
  !----  change some values in the case of regular PREM with two crustal layers or of 3D models
  !----

  ! case of regular PREM with two crustal layers: change the time step for small meshes
  ! because of a different size of elements in the radial direction in the crust
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

  ! minimum width of chunk
  min_chunk_width_in_degrees = min(ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES)

  ! gets attenuation min/max range
  if (.not. ATTENUATION_RANGE_PREDEFINED) then
     call auto_attenuation_periods(min_chunk_width_in_degrees, NEX_MAX)
  endif

  ! note: global estimate above for DT is empirical for chunk sizes of 90 degrees.
  ! if regional chunk size is larger, the time step is still limited by the vertical size of elements,
  ! which depends only on the number of vertical element layers.
  ! this however gets estimated above for 90 degrees already, so time stepping remains unchanged.

  ! adapts number of layer elements and time step size
  ! (for regional simulations with chunk sizes < 90 degrees)
  ! (or very, very large meshes)
  if (min_chunk_width_in_degrees < 90.0d0 .or. NEX_MAX > 1248) then
    ! adapts number of layer elements and time step size

    ! note: for global simulations, we set
    !         ANGULAR_WIDTH_XI_IN_DEGREES = 90.0d0 and
    !         ANGULAR_WIDTH_ETA_IN_DEGREES = 90.0d0 in read_parameter_file.f90

    ! gets number of element-layers
    call auto_ner(min_chunk_width_in_degrees, NEX_MAX)

    ! gets attenuation min/max range
    call auto_attenuation_periods(min_chunk_width_in_degrees, NEX_MAX)

    ! gets time step size
    call auto_time_stepping(min_chunk_width_in_degrees, NEX_MAX, dt_auto)

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

  endif

!---
!
! ADD YOUR MODEL HERE
!
!---


  ! time step reductions are based on empirical values (..somehow)

  ! following models need special attention, at least for global simulations:
  if (NCHUNKS == 6) then
    ! makes time step smaller for this ref model, otherwise becomes unstable in fluid
    if (REFERENCE_1D_MODEL == REFERENCE_MODEL_IASP91) &
      DT = DT*(1.d0 - 0.3d0)

    ! using inner core anisotropy, simulations might become unstable in solid
    if (ANISOTROPIC_INNER_CORE) then
      ! DT = DT*(1.d0 - 0.1d0) not working yet...
      stop 'anisotropic inner core - unstable feature, uncomment this line in get_timestep_and_layers.f90'
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
      if (REGIONAL_MOHO_MESH_EUROPE ) DT = 0.17 ! Europe
      if (REGIONAL_MOHO_MESH_ASIA ) DT = 0.15 ! Asia & Middle East
    endif
  endif

! the maximum CFL of LDDRK is significantly higher than that of the Newmark scheme,
! in a ratio that is theoretically 1.327 / 0.697 = 1.15 / 0.604 = 1.903 for a solid with Poisson's ratio = 0.25
! and for a fluid (see the manual of the 2D code, SPECFEM2D, Tables 4.1 and 4.2, and that ratio does not
! depend on whether we are in 2D or in 3D). However in practice a ratio of about 1.5 to 1.7 is often safer
! (for instance for models with a large range of Poisson's ratio values).
! Since the code computes the time step using the Newmark scheme, for LDDRK we simply
! multiply that time step by this ratio when LDDRK is on and when flag INCREASE_CFL_FOR_LDDRK is true.
  if (USE_LDDRK .and. INCREASE_CFL_FOR_LDDRK) DT = DT * RATIO_BY_WHICH_TO_INCREASE_IT


  ! for future usage:
  ! cut at a significant number of digits (2 digits)
  ! in steps of 1/2 digits
  ! example: 0.0734815 -> 0.0730
  !          0.07371   -> 0.0735
  !call get_timestep_limit_significant_digit(DT)

  end subroutine get_timestep_and_layers

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_timestep_limit_significant_digit(time_step)

  ! cut at a significant number of digits (2 digits) using 1/2 rounding
  ! example: 0.0734815 -> 0.0730
  !      and 0.0737777 -> 0.0735

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
