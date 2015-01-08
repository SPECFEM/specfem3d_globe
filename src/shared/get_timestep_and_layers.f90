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


  subroutine get_timestep_and_layers(NEX_MAX)

  use constants
  use shared_parameters

  implicit none

  integer,intent(in) :: NEX_MAX

  ! local variables
  integer :: multiplication_factor
  double precision :: min_chunk_width_in_degrees

  !----
  !----  case prem_onecrust by default
  !----
  if (SUPPRESS_CRUSTAL_MESH) then
    multiplication_factor = 2
  else
    multiplication_factor = 1
  endif

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
      if (NER_CRUST < 2 ) NER_CRUST = 2
      ! makes time step smaller
      if (NEX_MAX*multiplication_factor <= 160) then
        DT = 0.20d0
      else if (NEX_MAX*multiplication_factor <= 256) then
        DT = 0.20d0
      endif
    endif
  else
    ! 3D models: must have two element layers for crust
    if (NER_CRUST < 2 ) NER_CRUST = 2
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
     call auto_attenuation_periods(min_chunk_width_in_degrees, NEX_MAX, &
                                   MIN_ATTENUATION_PERIOD, MAX_ATTENUATION_PERIOD)
  endif

  ! adapts number of layer elements and time step size
  if (ANGULAR_WIDTH_XI_IN_DEGREES  < 90.0d0 .or. &
      ANGULAR_WIDTH_ETA_IN_DEGREES < 90.0d0 .or. &
      NEX_MAX > 1248) then

    ! gets number of element-layers
    call auto_ner(min_chunk_width_in_degrees, NEX_MAX, &
                  NER_CRUST, NER_80_MOHO, NER_220_80, NER_400_220, NER_600_400, &
                  NER_670_600, NER_771_670, NER_TOPDDOUBLEPRIME_771, &
                  NER_CMB_TOPDDOUBLEPRIME, NER_OUTER_CORE, NER_TOP_CENTRAL_CUBE_ICB, &
                  R_CENTRAL_CUBE, CASE_3D, CRUSTAL, &
                  HONOR_1D_SPHERICAL_MOHO, REFERENCE_1D_MODEL)

    ! gets attenuation min/max range
    call auto_attenuation_periods(min_chunk_width_in_degrees, NEX_MAX, &
                                  MIN_ATTENUATION_PERIOD, MAX_ATTENUATION_PERIOD)

    ! gets time step size
    call auto_time_stepping(min_chunk_width_in_degrees, NEX_MAX, DT)

    !! DK DK suppressed because this routine should not write anything to the screen
    !    write(*,*)'##############################################################'
    !    write(*,*)
    !    write(*,*)' Auto Radial Meshing Code '
    !    write(*,*)' Consult read_compute_parameters.f90 and auto_ner.f90 '
    !    write(*,*)' This should only be invoked for chunks less than 90 degrees'
    !    write(*,*)' and for chunks greater than 1248 elements wide'
    !    write(*,*)
    !    write(*,*)'CHUNK WIDTH:              ', min_chunk_width_in_degrees
    !    write(*,*)'NEX:                      ', NEX_MAX
    !    write(*,*)'NER_CRUST:                ', NER_CRUST
    !    write(*,*)'NER_80_MOHO:              ', NER_80_MOHO
    !    write(*,*)'NER_220_80:               ', NER_220_80
    !    write(*,*)'NER_400_220:              ', NER_400_220
    !    write(*,*)'NER_600_400:              ', NER_600_400
    !    write(*,*)'NER_670_600:              ', NER_670_600
    !    write(*,*)'NER_771_670:              ', NER_771_670
    !    write(*,*)'NER_TOPDDOUBLEPRIME_771:  ', NER_TOPDDOUBLEPRIME_771
    !    write(*,*)'NER_CMB_TOPDDOUBLEPRIME:  ', NER_CMB_TOPDDOUBLEPRIME
    !    write(*,*)'NER_OUTER_CORE:           ', NER_OUTER_CORE
    !    write(*,*)'NER_TOP_CENTRAL_CUBE_ICB: ', NER_TOP_CENTRAL_CUBE_ICB
    !    write(*,*)'R_CENTRAL_CUBE:           ', R_CENTRAL_CUBE
    !    write(*,*)'multiplication factor:    ', multiplication_factor
    !    write(*,*)
    !    write(*,*)'DT:                       ',DT
    !    write(*,*)'MIN_ATTENUATION_PERIOD    ',MIN_ATTENUATION_PERIOD
    !    write(*,*)'MAX_ATTENUATION_PERIOD    ',MAX_ATTENUATION_PERIOD
    !    write(*,*)
    !    write(*,*)'##############################################################'

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
  endif

  ! following models need special attention, regardless of number of chunks:
  ! it makes the time step smaller for this ref model, otherwise becomes unstable in fluid
  if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1066A) &
    DT = DT*(1.d0 - 0.8d0)  ! *0.20d0

  ! reduces time step size for "no mud" version of AK135F model
  if (REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135F_NO_MUD) &
    DT = DT*(1.d0 - 0.05d0)

  ! reduces time step size for crustmaps crustal model
  if (ITYPE_CRUSTAL_MODEL == ICRUST_CRUSTMAPS ) &
    DT = DT*(1.d0 - 0.3d0)

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
    if (HONOR_1D_SPHERICAL_MOHO ) return

    ! original values
    !print*,'NER:',NER_CRUST
    !print*,'DT:',DT

    ! enforce 3 element layers
    NER_CRUST = 3

    ! increased stability, empirical
    DT = DT*(1.d0 + 0.5d0)

    if (REGIONAL_MOHO_MESH_EUROPE ) DT = 0.17 ! Europe
    if (REGIONAL_MOHO_MESH_ASIA ) DT = 0.15 ! Asia & Middle East

  endif

! the maximum CFL of LDDRK is significantly higher than that of the Newmark scheme,
! in a ratio that is theoretically 1.327 / 0.697 = 1.15 / 0.604 = 1.903 for a solid with Poisson's ratio = 0.25
! and for a fluid (see the manual of the 2D code, SPECFEM2D, Tables 4.1 and 4.2, and that ratio does not
! depend on whether we are in 2D or in 3D). However in practice a ratio of about 1.5 to 1.7 is often safer
! (for instance for models with a large range of Poisson's ratio values).
! Since the code computes the time step using the Newmark scheme, for LDDRK we simply
! multiply that time step by this ratio when LDDRK is on and when flag INCREASE_CFL_FOR_LDDRK is true.
  if (USE_LDDRK .and. INCREASE_CFL_FOR_LDDRK) DT = DT * RATIO_BY_WHICH_TO_INCREASE_IT


  end subroutine get_timestep_and_layers
