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
!
!  This portion of the SPECFEM3D Code was written by:
!  Brian Savage while at
!     California Institute of Technology
!     Department of Terrestrial Magnetism / Carnegie Institute of Washington
!     University of Rhode Island
!
!  It is based partially upon formulation in:
!
! @ARTICLE{KoTr02a,
! author={D. Komatitsch and J. Tromp},
! year=2002,
! title={Spectral-Element Simulations of Global Seismic Wave Propagation{-I. V}alidation},
! journal={Geophys. J. Int.},
! volume=149,
! number=2,
! pages={390-412},
! doi={10.1046/j.1365-246X.2002.01653.x}}
!
!  and the core determination was developed.
!

  subroutine auto_time_stepping(WIDTH, NEX_MAX, DT)

  use constants, only: DEGREES_TO_RADIANS, NGLLX, &
    REFERENCE_MODEL_PREM2,REFERENCE_MODEL_IASP91,REFERENCE_MODEL_AK135F_NO_MUD, &
    REFERENCE_MODEL_1066A,REFERENCE_MODEL_1DREF,REFERENCE_MODEL_JP1D,REFERENCE_MODEL_SEA1D, &
    REFERENCE_MODEL_CCREM, &
    REFERENCE_MODEL_SOHL,REFERENCE_MODEL_SOHL_B,REFERENCE_MODEL_CASE65TAY, &
    REFERENCE_MODEL_VPREMOON,REFERENCE_MODEL_MOON_MEENA

  use shared_parameters, only: REFERENCE_1D_MODEL,R_PLANET,RICB

  implicit none

  integer,intent(in) :: NEX_MAX
  double precision, intent(in) :: WIDTH
  double precision, intent(inout) :: DT

  ! local parameters
  double precision :: RADIAL_LEN_RATIO_CENTRAL_CUBE
  double precision :: RADIUS_INNER_CORE
  double precision :: DOUBLING_INNER_CORE
  double precision :: P_VELOCITY_MAX ! Located Near the inner Core Boundary
  double precision :: MAXIMUM_STABILITY_CONDITION
  double precision :: MIN_GLL_POINT_SPACING
  double precision :: elem_size,min_grid_dx
  double precision :: dt_suggested
  double precision :: RADIAL_LEN_RATIO_CRUST,RADIUS_SURFACE,P_VELOCITY_MAX_CRUST,dt_suggested_crust
  logical :: check_crust_DT

  ! initializes defaults
  RADIAL_LEN_RATIO_CENTRAL_CUBE   =     0.40d0
  RADIAL_LEN_RATIO_CRUST          =     0.40d0

  ! conservative stability limit
  MAXIMUM_STABILITY_CONDITION     =     0.40d0

  ! note: the following assumes that the time step is mostly limited by the inner core,
  !       where the highest P-velocities appear.
  !       this might however not always be the case, e.g., when there is strong crustal variations
  !       and element sizes vary due to moho stretching.
  !       maybe this could be improved in future...

  DOUBLING_INNER_CORE             =      8.0d0

  ! default for PREM (near inner core boundary)
  RADIUS_INNER_CORE               = RICB / 1000.d0 ! RICB in km
  P_VELOCITY_MAX                  = 11.02827d0 ! km/s vp: 11.26220 - 6.36400 * (1221.49/6371.)**2

  ! default for NGLLX == 5
  MIN_GLL_POINT_SPACING         =   0.1730d0

  ! radius surface/crust
  RADIUS_SURFACE = R_PLANET / 1000.d0    ! in km

  ! also check constraint on DT from crust
  check_crust_DT = .false.

  ! modifies maximum velocity according to reference 1D model
  select case (REFERENCE_1D_MODEL)
  case (REFERENCE_MODEL_PREM2)
    P_VELOCITY_MAX = 11.06003d0   ! vp: 11.3041 - 1.2730 * (1221.5/6371.)

  case (REFERENCE_MODEL_IASP91)
    P_VELOCITY_MAX = 11.09147d0   ! vp: 11.24094 - 4.09689 * (1216.9/6371.)**2

  case (REFERENCE_MODEL_AK135F_NO_MUD)
    P_VELOCITY_MAX = 11.0427d0    ! vp

  case (REFERENCE_MODEL_1066A)
    P_VELOCITY_MAX = 10.9687d0    ! vp

  case (REFERENCE_MODEL_1DREF)
    P_VELOCITY_MAX = 11.02827d0   ! vpv (PREM)

  case (REFERENCE_MODEL_JP1D)
    P_VELOCITY_MAX = 11.09147d0   ! vp: 11.24094 - 4.09689 * x**2 (IASP91)

  case (REFERENCE_MODEL_SEA1D)
    P_VELOCITY_MAX = 11.09142d0   ! vp

  case (REFERENCE_MODEL_CCREM)
    P_VELOCITY_MAX = 11.2636d0    ! vp

  ! Mars models
  case (REFERENCE_MODEL_SOHL, &
        REFERENCE_MODEL_SOHL_B, &
        REFERENCE_MODEL_CASE65TAY)
    ! Mars
    ! note: for mars, the time stepping is mostly affected by crustal elements.
    !       we will use two estimates, one for inner core and another for the crust to determine a minimum time step.
    ! inner core
    P_VELOCITY_MAX = 7.3d0 * 1.1d0        ! vp: 11.26220 - 6.36400 * (1221.49/6371.)**2; daniel: increase by a factor 1.1x
    RADIAL_LEN_RATIO_CENTRAL_CUBE = 0.76  ! for an aspect ratio around 1.3
    ! surface/crust
    check_crust_DT = .true.
    P_VELOCITY_MAX_CRUST = 7.73d0   ! according to crustmap marscrustp7.cmap (lower crust layer) files
    ! empirical factor to account for aspect ratio in crust
    if (NEX_MAX < 480) then
      ! allows for larger time steps
      RADIAL_LEN_RATIO_CRUST = 0.85
    else
      ! takes stretching effect into account which will lead to thinner elements closer to surface
      RADIAL_LEN_RATIO_CRUST = 0.46
    endif

  ! Moon models
  case (REFERENCE_MODEL_VPREMOON, &
        REFERENCE_MODEL_MOON_MEENA)
    ! Moon
    ! uses two estimates, one for inner core and another for the crust to determine a minimum time step.
    ! inner core
    P_VELOCITY_MAX = 7.3d0 * 1.1d0
    RADIAL_LEN_RATIO_CENTRAL_CUBE = 0.76
    ! surface/crust
    check_crust_DT = .true.
    P_VELOCITY_MAX_CRUST = 5.5d0   ! according to VPREMOON (lower crust layer)
    ! empirical factor to account for aspect ratio in crust
    if (NEX_MAX < 480) then
      ! allows for larger time steps
      RADIAL_LEN_RATIO_CRUST = 0.85
    else
      ! takes stretching effect into account which will lead to thinner elements closer to surface
      RADIAL_LEN_RATIO_CRUST = 0.46
    endif

  end select

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
    stop 'auto_time_stepping: NGLLX > 10 value not supported yet! please consider adding it...'

  end select

  ! element at inner core
  ! size in horizontal direction
  elem_size = RADIAL_LEN_RATIO_CENTRAL_CUBE * ((WIDTH * DEGREES_TO_RADIANS) * RADIUS_INNER_CORE) / &
                ( dble(NEX_MAX) / DOUBLING_INNER_CORE )

  ! minimum grid point spacing
  min_grid_dx = elem_size * MIN_GLL_POINT_SPACING

  ! estimated time step
  dt_suggested = MAXIMUM_STABILITY_CONDITION * min_grid_dx / P_VELOCITY_MAX

  ! return as DT
  DT = dt_suggested

  !debug
  !print *,'debug: auto_time_stepping: inner core elem size',elem_size,'width/nex',WIDTH,NEX_MAX,'DT',DT

  ! crust time step
  if (check_crust_DT) then
    ! Mars & Moon models
    ! crustal elements
    elem_size = RADIAL_LEN_RATIO_CRUST * ((WIDTH * DEGREES_TO_RADIANS) * RADIUS_SURFACE) / dble(NEX_MAX)
    ! estimated time step
    dt_suggested_crust = MAXIMUM_STABILITY_CONDITION * MIN_GLL_POINT_SPACING * elem_size / P_VELOCITY_MAX_CRUST
    ! minimum suggested time step
    DT = min(dt_suggested,dt_suggested_crust)
    !debug
    !print *,'debug: auto_time_stepping: mars crust elem size ',elem_size, &
    !        'dt_suggested,dt_suggested_crust',dt_suggested,dt_suggested_crust
  endif

  end subroutine auto_time_stepping

!
!-------------------------------------------------------------------------------------------------
!
  subroutine auto_attenuation_periods(WIDTH, NEX_MAX, MIN_ATTENUATION_PERIOD, MAX_ATTENUATION_PERIOD)

  use constants, only: N_SLS

  use shared_parameters, only: T_min_period

  implicit none

  double precision,intent(in) :: WIDTH
  integer, intent(in) :: NEX_MAX
  double precision, intent(inout) :: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD

  ! local parameters
  double precision :: tmp
  double precision :: THETA(5)
  double precision,parameter :: TOL_ZERO = 1.d-30

  ! safety check
  if (N_SLS < 2 .or. N_SLS > 5) &
    stop 'N_SLS must be greater than 1 or less than 6'
  if (NEX_MAX <= 0) &
    stop 'Invalid NEX_MAX in auto_attenuation_periods()'
  if (WIDTH <= 0.d0) &
    stop 'Invalid WIDTH in auto_attenuation_periods()'

  ! The Minimum attenuation period (integer)
  MIN_ATTENUATION_PERIOD = T_min_period

  ! THETA defines the width of the Attenuation Range in Decades
  !   The number defined here were determined by minimizing
  !   the "flatness" of the absorption spectrum.  Each THETA
  !   is defined for a particular N_SLS (constants.h)
  !   THETA(2) is for N_SLS = 2
  THETA(1)           =   0.00d0
  THETA(2)           =   0.75d0
  THETA(3)           =   1.75d0
  THETA(4)           =   2.25d0
  THETA(5)           =   2.85d0

  ! Compute Max Attenuation Period
  !
  ! The max attenuation period for 3 SLS is optimally
  !   1.75 decades from the min attenuation period, see THETA above
  tmp = MIN_ATTENUATION_PERIOD * 10.0d0**THETA(N_SLS)

  ! maximum period
  MAX_ATTENUATION_PERIOD = tmp

  ! check
  if (MIN_ATTENUATION_PERIOD <= 0.d0) then
    print *,'Error: invalid attenuation minimum: ',MIN_ATTENUATION_PERIOD
    stop 'Invalid attenuation minimum'
  endif
  if (MAX_ATTENUATION_PERIOD <= 0.d0) then
    print *,'Error: invalid attenuation maximum: ',MAX_ATTENUATION_PERIOD
    stop 'Invalid attenuation maximum'
  endif
  if (abs(MAX_ATTENUATION_PERIOD - MIN_ATTENUATION_PERIOD) < TOL_ZERO .or. &
      MAX_ATTENUATION_PERIOD < MIN_ATTENUATION_PERIOD) then
    print *,'Error: invalid attenuation range min/max: ',MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD
    stop 'Invalid attenuation range min/max'
  endif

  end subroutine auto_attenuation_periods

!
!-------------------------------------------------------------------------------------------------
!

  subroutine auto_ner(WIDTH, NEX_MAX)

  use shared_parameters, only: PLANET_TYPE,IPLANET_EARTH,IPLANET_MARS,IPLANET_MOON,R_PLANET

  use shared_parameters, only: &
    NER_CRUST, NER_80_MOHO, NER_220_80, NER_400_220, NER_600_400, &
    NER_670_600, NER_771_670, NER_TOPDDOUBLEPRIME_771, &
    NER_CMB_TOPDDOUBLEPRIME, NER_OUTER_CORE, NER_TOP_CENTRAL_CUBE_ICB, &
    R_CENTRAL_CUBE,CASE_3D

  use shared_parameters, only: &
    R80,R220,R400,R600,R670,R771, &
    RTOPDDOUBLEPRIME,RCMB, &
    RMOHO_FICTITIOUS_IN_MESHER

  implicit none

  double precision,intent(in) :: WIDTH
  integer,intent(in) :: NEX_MAX

  ! local parameters
  integer,          parameter                :: NUM_REGIONS = 14
  integer,          dimension(NUM_REGIONS)   :: scaling
  double precision, dimension(NUM_REGIONS)   :: radius
  double precision, dimension(NUM_REGIONS-1) :: ratio_top
  double precision, dimension(NUM_REGIONS-1) :: ratio_bottom
  integer,          dimension(NUM_REGIONS-1) :: NER

! uses model specific radii to determine number of elements in radial direction
! (set by earlier call to routine get_model_parameters_radii())

  ! radii
  radius(1)  = R_PLANET ! Surface radius
  radius(2)  = RMOHO_FICTITIOUS_IN_MESHER   ! Moho - 1st Mesh Doubling Interface
  radius(3)  = R80
  radius(4)  = R220
  radius(5)  = R400
  radius(6)  = R600
  radius(7)  = R670
  radius(8)  = R771

  ! specific radii for planets
  select case (PLANET_TYPE)
  case (IPLANET_EARTH)
    ! Earth
    ! This is PREM in Kilometers, well ... kinda, not really ....
    !radius(1)  = 6371.00d0 ! Surface
    !radius(2)  = 6346.60d0 !    Moho - 1st Mesh Doubling Interface
    !radius(3)  = 6291.60d0 !      80
    !radius(4)  = 6151.00d0 !     220
    !radius(5)  = 5971.00d0 !     400
    !radius(6)  = 5771.00d0 !     600
    !radius(7)  = 5701.00d0 !     670
    !radius(8)  = 5600.00d0 !     771
    !radius(9)  = 4712.00d0 !    1650 - 2nd Mesh Doubling: Geochemical Layering; Kellogg et al. 1999, Science
    !radius(10) = 3630.00d0 !     D_double_prime
    !radius(11) = 3480.00d0 !     CMB
    !radius(12) = 2511.00d0 !    3860 - 3rd Mesh Doubling Interface
    !radius(13) = 1371.00d0 !    5000 - 4th Mesh Doubling Interface
    !radius(14) =  982.00d0 ! Top Central Cube

    radius(9)  = 4712000.0d0 !    1650 - 2nd Mesh Doubling: Geochemical Layering; Kellogg et al. 1999, Science
    radius(10) = RTOPDDOUBLEPRIME   !     D_double_prime ~ 3630
    radius(11) = RCMB   !     CMB ~ 3480

    radius(12) = 2511000.0d0 !    3860 - 3rd Mesh Doubling Interface
    radius(13) = 1371000.0d0 !    5000 - 4th Mesh Doubling Interface
    radius(14) = R_CENTRAL_CUBE ! Top Central Cube

  case (IPLANET_MARS)
    ! Mars
    ! note: radii R_.. are set according to Mars model geometry
    radius(9)  = 1900000.0d0          ! in between R771 (at 2033km) and RTOPDOUBLEPRIME (at 1503km)
    radius(10) = RTOPDDOUBLEPRIME     ! depth = 1887 km, radius 1503000.0 m
    radius(11) = RCMB                 ! depth = 1922 km

    radius(12) = 1300000.0d0          ! depth = 2090 km - 3rd Mesh Doubling Interface
    radius(13) = 700000.0d0           ! depth = 2690 km - 4th Mesh Doubling Interface
    radius(14) = R_CENTRAL_CUBE       ! Top Central Cube (and ICB at RICB = 515000.0)

  case (IPLANET_MOON)
    ! Moon
    ! note: radii R_.. are set according to Moon model geometry
    radius(9)  = 800000.0d0           ! in between R771 (radius = 966.1 km) and RTOPDOUBLEPRIME
    radius(10) = RTOPDDOUBLEPRIME     ! D" radius = 480 km
    radius(11) = RCMB                 ! RCMB radius = 380 km

    radius(12) = 335000.0d0           ! 3rd Mesh Doubling Interface
    radius(13) = 285000.0d0           ! 4th Mesh Doubling Interface
    radius(14) = R_CENTRAL_CUBE       ! Top Central Cube (and ICB at RICB radius = 240 km)

  case default
    ! avoiding exit_MPI(), since we also call this routine in create_header_file
    ! which can be compiled without MPI - using stop instead
    !call exit_MPI(myrank,'Invalid planet, auto_ner() not implemented yet')
    print *,'Invalid planet, auto_ner() not implemented yet'
    stop 'Invalid planet, auto_ner() not implemented yet'
  end select

  ! radii in km
  radius(:) = radius(:) / 1000.0d0

  call find_r_central_cube(NEX_MAX, radius(14))

  ! Mesh Doubling
  scaling(1)     = 1  ! SURFACE TO MOHO
  scaling(2:8)   = 2  ! MOHO    TO G_double_prime (Geochemical Mantle 1650)
  scaling(9:11)  = 4  ! G_double_prime    TO MIC (Middle Inner Core)
  scaling(12)    = 8  ! MIC     TO MIC-II
  scaling(13:14) = 16 ! MIC-II  TO Central Cube TO Center of the Earth

  ! initializes minimum Number of Elements a Region must have
  NER(:)    = 1
  NER(3:5)  = 2
  if (CASE_3D) then
    NER(1) = 2
  endif

  ! specifies minimum element layers (in vertical direction) for top layer
  if (PLANET_TYPE == IPLANET_MARS .or. PLANET_TYPE == IPLANET_MOON) then
    ! Mars & Moon
    NER(1) = 2
  endif

  ! starts from input arguments of a 90-degree chunk
  ! (where NER values are set empirically for a good mesh design)
  NER(1) = NER_CRUST
  NER(2) = NER_80_MOHO
  NER(3) = NER_220_80
  NER(4) = NER_400_220
  NER(5) = NER_600_400
  NER(6) = NER_670_600
  NER(7) = NER_771_670
  ! distributes NER_TOPDDOUBLEPRIME_771 onto two element layer regions depending on vertical sizes of layers
  NER(8) = max(int( NER_TOPDDOUBLEPRIME_771 * (radius(8) - radius(9)) / (radius(8) - radius(10)) ), 1)
  NER(9) = NER_TOPDDOUBLEPRIME_771 - NER(8)
  NER(10) = NER_CMB_TOPDDOUBLEPRIME
  ! distributes NER_OUTER_CORE onto two element layer regions depending on vertical sizes of layers
  NER(11) = max(int( NER_OUTER_CORE * (radius(11) - radius(12)) / (radius(11) - radius(13)) ), 1)
  NER(12) = NER_OUTER_CORE - NER(11)
  NER(13) = NER_TOP_CENTRAL_CUBE_ICB

  ! debug
  !print *,'debug: input NER:',NER(:)

  ! Find the Number of Radial Elements in a region based upon
  ! the aspect ratio of the elements
  call auto_optimal_ner(NUM_REGIONS, WIDTH, NEX_MAX, radius, scaling, NER, ratio_top, ratio_bottom)

  ! debug
  !print *,'debug: output NER:',NER(:)

  ! Set Output arguments
  NER_CRUST                = NER(1)
  NER_80_MOHO              = NER(2)
  NER_220_80               = NER(3)
  NER_400_220              = NER(4)
  NER_600_400              = NER(5)
  NER_670_600              = NER(6)
  NER_771_670              = NER(7)
  NER_TOPDDOUBLEPRIME_771  = NER(8) + NER(9)
  NER_CMB_TOPDDOUBLEPRIME  = NER(10)
  NER_OUTER_CORE           = NER(11) + NER(12)
  NER_TOP_CENTRAL_CUBE_ICB = NER(13)

  R_CENTRAL_CUBE           = radius(14) * 1000.0d0 ! converts to m

  end subroutine auto_ner

!
!-------------------------------------------------------------------------------------------------
!

  subroutine auto_optimal_ner(NUM_REGIONS, width, NEX, r, scaling, NER, rt, rb)

  use constants, only: DEGREES_TO_RADIANS
  use shared_parameters, only: PLANET_TYPE,IPLANET_EARTH,IPLANET_MARS,IPLANET_MOON

  implicit none

  integer,intent(in) :: NUM_REGIONS
  integer,intent(in) :: NEX
  double precision,intent(in) ::  width   ! Width of the Chunk in Degrees

  integer,          dimension(NUM_REGIONS-1),intent(inout) :: NER      ! Elements per Region    - IN-N-OUT - Yummy !
  integer,          dimension(NUM_REGIONS)  ,intent(in)    :: scaling  ! Element Doubling       - INPUT
  double precision, dimension(NUM_REGIONS)  ,intent(in)    :: r        ! Radius                 - INPUT
  double precision, dimension(NUM_REGIONS-1),intent(out)   :: rt       ! Ratio at Top           - OUTPUT
  double precision, dimension(NUM_REGIONS-1),intent(out)   :: rb       ! Ratio at Bottom        - OUTPUT

  ! local parameters
  double precision :: dr, w, ratio, xi, ximin, wt, wb
  integer :: ner_test
  integer :: i

  ! target ratio
  double precision :: aspect_ratio

  ! Earth: optimal ratio around 1 (almost cubic shapes)
  double precision,parameter :: ASPECT_RATIO_EARTH = 1.d0
  ! Mars: due to a smaller radius, the optimal ratio around 1 is too optimistic and
  !       tends to drastically increase the number of element layers;
  !       we thus allow for more elongated elements with an empirical ratio ~ 1.5
  double precision,parameter :: ASPECT_RATIO_MARS = 1.5d0

  !debug
  logical, parameter :: DEBUG = .false.

  ! sets target aspect ratio
  select case(PLANET_TYPE)
  case (IPLANET_MARS,IPLANET_MOON)
    ! Mars
    aspect_ratio = ASPECT_RATIO_MARS
  case (IPLANET_EARTH)
    ! Earth
    aspect_ratio = ASPECT_RATIO_EARTH
  case default
    print *,'Invalid planet, auto_optimal_ner not implemented yet'
    stop 'Invalid planet, auto_optimal_ner not implemented yet'
  end select

  ! debug
  if (DEBUG) print *,'planet: ',PLANET_TYPE,'(1 == earth,2 == mars, 3 == moon) - target ratio ',aspect_ratio

  ! Find optimal elements per region
  do i = 1,NUM_REGIONS-1
    dr = r(i) - r(i+1)              ! Radial Length of Region
    wt = width * DEGREES_TO_RADIANS * r(i)   / (NEX*1.0d0 / scaling(i)*1.0d0) ! Element Width Top
    wb = width * DEGREES_TO_RADIANS * r(i+1) / (NEX*1.0d0 / scaling(i)*1.0d0) ! Element Width Bottom
    w  = (wt + wb) * 0.5d0          ! Average Width of Region
    ner_test = NER(i)               ! Initial solution

    ! checks
    if (ner_test == 0) then
      print *,'Error auto_ner: region',i,'ner_test = ',ner_test,'width w/wtop/wbottom =',w,wt,wb,'radial length dr = ',dr
      stop 'Error invalid ner value in auto_ner'
    endif

    ratio = (dr / ner_test) / w     ! Aspect Ratio of Element
    xi = dabs(ratio - aspect_ratio) ! Aspect Ratio should be near 1.0
    ximin = 1.e7                    ! Initial Minimum

    !debug
    if (DEBUG) print *,'debug auto_optimal_ner: region ',i, &
                       'element initial ratio: ',sngl(ratio),'xi = ',sngl(xi),'width = ',sngl(w),'ner',ner_test

    ! check
    if (ratio < 0.d0) then
      print *,'Error: auto optimal ner has negative aspect ratio of element:',ratio
      print *,'       region',i,'radius',r(i),r(i+1),dr,'width',w
      stop 'Error auto_optimal_ner() with negative aspect ratio of element'
    endif

    ! increases NER to reach vertical/horizontal element ratio of about 1
    do while(xi <= ximin)
      NER(i) = ner_test                 ! Found a better solution
      ximin = xi
      ner_test = ner_test + 1           ! Increment ner_test
      ratio = (dr / ner_test) / w
      xi = dabs(ratio - aspect_ratio)   ! look for a better solution
    enddo
    rt(i) = dr / NER(i) / wt        ! Find the Ratio of Top
    rb(i) = dr / NER(i) / wb        ! and Bottom for completeness

    !debug
    if (DEBUG) then
      print *,'debug auto_optimal_ner: region ',i, &
              'element final ratio: top = ',sngl(rt(i)),'bottom = ',sngl(rb(i)),'ner',NER(i)
      print *
    endif
  enddo

  end subroutine auto_optimal_ner

!
!-------------------------------------------------------------------------------------------------
!

  subroutine find_r_central_cube(nex_xi_in, rcube)

  use shared_parameters, only: PLANET_TYPE,IPLANET_EARTH,IPLANET_MARS,IPLANET_MOON

  implicit none

  integer, parameter :: NBNODE = 8
  double precision, parameter :: alpha = 0.41d0

  integer,intent(in) :: nex_xi_in
  double precision,intent(inout) :: rcube

  ! local parameters
  integer :: npts
  integer :: nex_xi
  integer :: nex_eta
  double precision :: rcubestep, rcube_test, rcubemax
  double precision :: xi, ximin
  double precision, allocatable, dimension(:,:) :: points
  double precision :: elem(NBNODE+1, 2)
  integer :: nspec_cube, nspec_chunks, ispec, nspec
  double precision :: edgemax, edgemin
  double precision :: max_edgemax, min_edgemin
  double precision :: aspect_ratio, max_aspect_ratio

  ! debug
  !print *,'debug: inner core radius initial = ',sngl(rcube)

  ! inner core radius (in km)
  select case(PLANET_TYPE)
  case (IPLANET_EARTH)
    ! Earth
    rcube_test   =  930.0d0  ! start
    rcubemax     = MAX(1100.0d0,rcube)  ! maximum size

  case (IPLANET_MARS)
    ! Mars
    rcube_test   =  200.0d0  ! start
    rcubemax     =  MAX(300.0d0,rcube)  ! maximum size

  case (IPLANET_MOON)
    ! Moon
    rcube_test   =  150.0d0  ! start
    rcubemax     =  MAX(220.0d0,rcube)  ! maximum size, RICB at 240km

  case default
    ! avoiding exit_MPI(), since we also call this routine in create_header_file
    ! which can be compiled without MPI - using stop instead
    !call exit_MPI(myrank,'Invalid planet, find r central cube not implemented yet')
    print *,'Invalid planet, find r central cube not implemented yet'
    stop 'Invalid planet, find r central cube not implemented yet'
  end select

  nex_xi = nex_xi_in / 16
  ximin        = 1e7

  rcube        = rcube_test
  rcubestep    = 1.0d0

  do while(rcube_test <= rcubemax)
    max_edgemax = -1e7
    min_edgemin = 1e7
    max_aspect_ratio = 0.0d0

    call compute_nex(nex_xi, rcube_test, alpha, nex_eta)

    npts = (4 * nex_xi * nex_eta * NBNODE) + (nex_xi * nex_xi * NBNODE)

    allocate(points(npts, 2))

    call compute_IC_mesh(rcube_test, points, npts, nspec_cube, nspec_chunks, nex_xi, nex_eta)

    nspec = nspec_cube + nspec_chunks
    do ispec = 1,nspec
      call get_element(points, ispec, npts, elem)
      call get_size_min_max(elem, edgemax, edgemin)

      aspect_ratio = edgemax / edgemin
      max_edgemax = MAX(max_edgemax, edgemax)
      min_edgemin = MIN(min_edgemin, edgemin)
      max_aspect_ratio = MAX(max_aspect_ratio, aspect_ratio)
    enddo
    ! maximum aspect ratio
    xi = (max_edgemax / min_edgemin)

    !xi = abs(rcube_test - 981.0d0) / 45.0d0

    ! debug
    !print '(a,6(f14.4,2x))','debug: rcube, xi, ximin, min/max/ratio', &
    !                        rcube_test, xi, ximin, min_edgemin,max_edgemax,max_aspect_ratio

    ! stores better aspect ratio
    if (xi < ximin) then
      ximin      = xi
      rcube      = rcube_test
    endif
    rcube_test = rcube_test + rcubestep

    deallocate(points)
  enddo

  ! debug
  !print *,'debug: inner core radius adjusted to = ',sngl(rcube)

  end subroutine find_r_central_cube

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_nex(nex_xi, rcube, alpha, ner)

  use constants, only: PI,PI_OVER_TWO,PI_OVER_FOUR
  use shared_parameters, only: RICB

  implicit none

  integer,intent(in) :: nex_xi
  double precision,intent(in) :: rcube, alpha
  integer,intent(out) :: ner

  ! local parameters
  integer :: ix
  double precision :: ratio_x, factx, xi
  double precision :: x, y
  double precision :: surfx, surfy
  double precision :: dist_cc_icb, somme, dist_moy
  double precision :: RICB_KM

  ! inner-core boundary in km
  RICB_KM = RICB / 1000.d0

  somme = 0.0d0

  do ix = 0,nex_xi/2,1
     ratio_x = (ix * 1.0d0) / ( nex_xi * 1.0d0)
     factx = 2.0d0 * ratio_x - 1.0d0
     xi = PI_OVER_TWO * factx
     x = (rcube / sqrt(2.0d0)) * factx
     y = (rcube / sqrt(2.0d0)) * (1 + cos(xi) * alpha / PI_OVER_TWO)

     surfx = RICB_KM * cos(3 * PI_OVER_FOUR - ratio_x * PI_OVER_TWO)
     surfy = RICB_KM * sin(3 * PI_OVER_FOUR - ratio_x * PI_OVER_TWO)

     dist_cc_icb = sqrt((surfx -x)**2 + (surfy - y)**2)
     if (ix /= nex_xi/2) then
        dist_cc_icb = dist_cc_icb * 2
     endif
     somme = somme + dist_cc_icb
  enddo
  dist_moy = somme / (nex_xi + 1)
  ner = nint(dist_moy / ((PI * RICB_KM) / (2*nex_xi)))

  end subroutine compute_nex

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_element(points, ispec, npts, pts)

  implicit none

  integer, parameter :: NBNODE = 8

  integer,intent(in) :: npts,ispec
  double precision,intent(in) :: points(npts,2)
  double precision,intent(out) :: pts(NBNODE+1,2)

  ! local parameters
  integer :: istart_left,istart_right,i

  istart_left = 1
  istart_right = (ispec-1)*NBNODE + 1
  do i = 0,NBNODE-1
    pts(istart_left + i,1) = points(istart_right + i,1)
    pts(istart_left + i,2) = points(istart_right + i,2)
  enddo
  pts(NBNODE+1,:) = pts(1,:)  ! use first point as the last point

  end subroutine get_element

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_size_min_max(pts, edgemax, edgemin)

  implicit none

  integer, parameter :: NBNODE = 8

  double precision,intent(in) :: pts(NBNODE+1, 2)
  double precision,intent(out) :: edgemax, edgemin

  ! local parameters
  integer :: ie, ix1,ix2,ix3
  double precision :: edge

  edgemax = -1e7
  edgemin = -edgemax
  do ie = 1,NBNODE/2,1
    ix1 = (ie * 2) - 1
    ix2 = ix1 + 1
    ix3 = ix1 + 2
    edge = sqrt( (pts(ix1,1) - pts(ix2,1))**2 + (pts(ix1,2) - pts(ix2,2))**2 ) + &
           sqrt( (pts(ix2,1) - pts(ix3,1))**2 + (pts(ix2,2) - pts(ix3,2))**2 )
    edgemax = MAX(edgemax, edge)
    edgemin = MIN(edgemin, edge)
  enddo

  end subroutine get_size_min_max

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_IC_mesh(rcube, points, npts, nspec_cube, nspec_chunks, nex_xi, nex_eta)

  implicit none

  integer, parameter :: NBNODE = 8

  integer,intent(in) :: npts
  integer,intent(out) :: nspec_chunks, nspec_cube

  double precision,intent(in) :: rcube
  double precision,intent(out) :: points(npts, 2)

  integer :: nex_eta, nex_xi

  ! local parameters
  double precision :: alpha
  double precision :: x, y
  integer :: ic, ix, iy, in
  ! topology of the elements
  integer, parameter, dimension(NBNODE) :: iaddx(NBNODE) = (/0,1,2,2,2,1,0,0/)
  integer, parameter, dimension(NBNODE) :: iaddy(NBNODE) = (/0,0,0,1,2,2,2,1/)
  integer :: k

  k = 1
  alpha = 0.41d0
  nspec_chunks = 0
  do ic = 0,3
     do ix = 0,(nex_xi-1)*2,2
        do iy = 0,(nex_eta-1)*2,2
           do in = 1,NBNODE
              call compute_coordinate(ix+iaddx(in), iy+iaddy(in), nex_xi*2, nex_eta*2, rcube, ic, alpha, x,y)
              points(k,1) = x
              points(k,2) = y
              k = k + 1
           enddo
           nspec_chunks = nspec_chunks + 1
        enddo
     enddo
  enddo

  nspec_cube = 0
  do ix = 0,(nex_xi-1)*2,2
     do iy = 0,(nex_xi-1)*2,2
        do in = 1,NBNODE
           call compute_coordinate_central_cube(ix+iaddx(in), iy+iaddy(in), nex_xi*2, nex_xi*2, rcube, alpha,x,y)
           points(k,1) = x
           points(k,2) = y
           k = k + 1
        enddo
        nspec_cube = nspec_cube + 1
     enddo
  enddo

  end subroutine compute_IC_mesh

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_coordinate_central_cube(ix,iy,nbx,nby,radius, alpha, x, y)

  use constants, only: PI_OVER_TWO

  implicit none

  integer,intent(in) :: ix, iy, nbx, nby
  double precision,intent(in) :: radius, alpha
  double precision,intent(out) :: x, y

  ! local parameters
  double precision :: ratio_x, ratio_y
  double precision :: factx, facty
  double precision :: xi, eta

  ratio_x = (ix * 1.0d0) / (nbx * 1.0d0)
  ratio_y = (iy * 1.0d0) / (nby * 1.0d0)

  factx = 2.0d0 * ratio_x - 1.0d0
  facty = 2.0d0 * ratio_y - 1.0d0

  xi  = PI_OVER_TWO * factx
  eta = PI_OVER_TWO * facty

  x = (radius / sqrt(2.0d0)) * factx * ( 1 + cos(eta) * alpha / PI_OVER_TWO )
  y = (radius / sqrt(2.0d0)) * facty * ( 1 + cos(xi)  * alpha / PI_OVER_TWO )

  end subroutine compute_coordinate_central_cube

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_coordinate(ix,iy,nbx, nby, rcube, ic, alpha, x, y)

  use constants, only: PI_OVER_TWO,PI_OVER_FOUR
  use shared_parameters, only: RICB

  implicit none

  integer,intent(in) :: ix, iy, nbx, nby, ic
  double precision,intent(in) :: rcube, alpha
  double precision,intent(out) :: x, y

  ! local parameters
  double precision :: ratio_x, ratio_y
  double precision :: factx, xi
  double precision :: xcc, ycc
  double precision :: xsurf, ysurf
  double precision :: deltax, deltay
  double precision :: temp

  double precision :: RICB_KM

  ! inner-core boundary in km
  RICB_KM = RICB / 1000.d0

  ratio_x = (ix * 1.0d0) / (nbx * 1.0d0)
  ratio_y = (iy * 1.0d0) / (nby * 1.0d0)

  factx = 2.0d0 * ratio_x - 1.0d0
  xi = PI_OVER_TWO * factx

  xcc = (rcube / sqrt(2.0d0)) * factx
  ycc = (rcube / sqrt(2.0d0)) * (1 + cos(xi) * alpha / PI_OVER_TWO)

  xsurf = RICB_KM * cos(3.0d0 * PI_OVER_FOUR - ratio_x * PI_OVER_TWO)
  ysurf = RICB_KM * sin(3.0d0 * PI_OVER_FOUR - ratio_x * PI_OVER_TWO)

  deltax = xsurf - xcc
  deltay = ysurf - ycc

  x = xsurf - ratio_y * deltax
  y = ysurf - ratio_y * deltay

  if (ic == 1) then
     temp = x
     x    = y
     y    = temp
  else if (ic == 2) then
     x = -x
     y = -y
  else if (ic == 3) then
     temp = x
     x    = -y
     y    = temp
  endif

  end subroutine compute_coordinate
