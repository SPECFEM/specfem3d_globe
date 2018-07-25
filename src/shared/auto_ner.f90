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

  subroutine auto_time_stepping(WIDTH,  NEX_MAX, DT)

  use constants, only: DEGREES_TO_RADIANS, NGLLX, &
    REFERENCE_MODEL_PREM,REFERENCE_MODEL_IASP91,REFERENCE_MODEL_AK135F_NO_MUD, &
    REFERENCE_MODEL_1066A,REFERENCE_MODEL_1DREF,REFERENCE_MODEL_JP1D,REFERENCE_MODEL_SEA1D

  use shared_parameters, only: REFERENCE_1D_MODEL

  implicit none

  integer,intent(in) :: NEX_MAX
  double precision, intent(in) :: WIDTH
  double precision, intent(out) :: DT

  ! local parameters
  double precision :: RADIAL_LEN_RATIO_CENTRAL_CUBE
  double precision :: RADIUS_INNER_CORE
  double precision :: DOUBLING_INNER_CORE
  double precision :: P_VELOCITY_MAX     ! Located Near the inner Core Boundary
  double precision :: MAXIMUM_STABILITY_CONDITION
  double precision :: MIN_GLL_POINT_SPACING
  double precision :: elem_size,min_grid_dx

  ! initializes defaults
  RADIAL_LEN_RATIO_CENTRAL_CUBE   =     0.40d0

  ! conservative stability limit
  MAXIMUM_STABILITY_CONDITION     =     0.40d0

  DOUBLING_INNER_CORE             =      8.0d0

  ! default for PREM (near inner core boundary)
  RADIUS_INNER_CORE               =   1221.0d0 ! RICB in km
  P_VELOCITY_MAX                  = 11.02827d0 ! in km/s

  ! default for NGLLX == 5
  MIN_GLL_POINT_SPACING         =   0.1730d0

  ! modifies maximum velocity according to reference 1D model
  select case (REFERENCE_1D_MODEL)
  case (REFERENCE_MODEL_PREM)
    RADIUS_INNER_CORE = 1221.0d0
    P_VELOCITY_MAX = 11.02827d0 ! vp: 11.26220 - 6.36400 * (1221.49/6371.)**2

  case (REFERENCE_MODEL_IASP91)
    RADIUS_INNER_CORE = 1217.0d0
    P_VELOCITY_MAX = 11.09147d0 ! vp: 11.24094 - 4.09689 * (1216.9/6371.)**2

  case (REFERENCE_MODEL_AK135F_NO_MUD)
    RADIUS_INNER_CORE = 1217.5d0
    P_VELOCITY_MAX = 11.0427d0 ! vp

  case (REFERENCE_MODEL_1066A)
    RADIUS_INNER_CORE = 1229.48d0
    P_VELOCITY_MAX = 10.9687d0 ! vp

  case (REFERENCE_MODEL_1DREF)
    RADIUS_INNER_CORE = 1221.491d0
    P_VELOCITY_MAX = 11.02827d0  ! vpv (PREM)

  case (REFERENCE_MODEL_JP1D)
    RADIUS_INNER_CORE = 1217.0d0
    P_VELOCITY_MAX = 11.09147d0 ! vp: 11.24094 - 4.09689 * x**2 (IASP91)

  case (REFERENCE_MODEL_SEA1D)
    RADIUS_INNER_CORE = 1217.1d0
    P_VELOCITY_MAX = 11.09142d0 ! vp

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
  DT = MAXIMUM_STABILITY_CONDITION * min_grid_dx / P_VELOCITY_MAX

  end subroutine auto_time_stepping

!
!-------------------------------------------------------------------------------------------------
!
  subroutine auto_attenuation_periods(WIDTH, NEX_MAX)

  use constants, only: N_SLS,NGLLX

  use shared_parameters, only: MIN_ATTENUATION_PERIOD, MAX_ATTENUATION_PERIOD

  implicit none

  double precision,intent(in) :: WIDTH
  integer, intent(in) :: NEX_MAX

  ! local parameters
  double precision :: TMP
  double precision :: GLL_SPACING
  double precision :: S_VELOCITY_MIN
  double precision :: THETA(5)

  ! factor degree to km
  double precision,parameter :: DEG2KM = 111.d0

  ! required points per wavelength
  double precision,parameter :: PTS_PER_WAVELENGTH = 4.d0

  ! safety check
  if (N_SLS < 2 .or. N_SLS > 5) then
     stop 'N_SLS must be greater than 1 or less than 6'
  endif

  ! average spacing between GLL points
  GLL_SPACING = dble(NGLLX - 1)

  ! minimum velocity (Vs)
  S_VELOCITY_MIN = 2.25d0

  ! Compute Min Attenuation Period
  !
  ! width of element in km = (Angular width in degrees / NEX_MAX) * degrees to km
  TMP = WIDTH / dble(NEX_MAX) * DEG2KM

  ! average grid node spacing in km = Width of an element in km / spacing for GLL point
  TMP = TMP / GLL_SPACING

  ! minimum resolved wavelength (for fixed number of points per wavelength)
  TMP = TMP * PTS_PER_WAVELENGTH

  ! The Minimum attenuation period = (minimum wavelength) / V_min
  TMP = TMP/S_VELOCITY_MIN

  ! The Minimum attenuation period (integer)
  MIN_ATTENUATION_PERIOD = int(TMP)

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
  TMP = TMP * 10.0d0**THETA(N_SLS)

  MAX_ATTENUATION_PERIOD = int(TMP)

  ! debug
  !print *,'attenuation range min/max: ',MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD

  end subroutine auto_attenuation_periods

!
!-------------------------------------------------------------------------------------------------
!

  subroutine auto_ner(WIDTH, NEX_MAX)

  use constants, only: R_EARTH

  use shared_parameters, only: &
    CASE_3D !, CRUSTAL, HONOR_1D_SPHERICAL_MOHO, REFERENCE_1D_MODEL

  use shared_parameters, only: &
    NER_CRUST, NER_80_MOHO, NER_220_80, NER_400_220, NER_600_400, &
    NER_670_600, NER_771_670, NER_TOPDDOUBLEPRIME_771, &
    NER_CMB_TOPDDOUBLEPRIME, NER_OUTER_CORE, NER_TOP_CENTRAL_CUBE_ICB, &
    R_CENTRAL_CUBE

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

!  double precision :: ROCEAN,RMIDDLE_CRUST, &
!          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
!          RMOHO_FICTITIOUS_IN_MESHER,R80_FICTITIOUS_IN_MESHER
!  double precision :: RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS

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

  ! uses model specific radii to determine number of elements in radial direction
  ! (set by earlier call to routine get_model_parameters_radii())

  radius(1)  = R_EARTH ! Surface
  radius(2)  = RMOHO_FICTITIOUS_IN_MESHER !    Moho - 1st Mesh Doubling Interface
  radius(3)  = R80    !      80
  radius(4)  = R220   !     220
  radius(5)  = R400   !     400
  radius(6)  = R600   !     600
  radius(7)  = R670   !     670
  radius(8)  = R771   !     771

  radius(9)  = 4712000.0d0 !    1650 - 2nd Mesh Doubling: Geochemical Layering; Kellogg et al. 1999, Science

  radius(10) = RTOPDDOUBLEPRIME   !     D_double_prime ~ 3630
  radius(11) = RCMB   !     CMB ~ 3480

  radius(12) = 2511000.0d0 !    3860 - 3rd Mesh Doubling Interface
  radius(13) = 1371000.0d0 !    5000 - 4th Mesh Doubling Interface
  radius(14) =  982000.0d0 ! Top Central Cube

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
  NER(8) = NER_TOPDDOUBLEPRIME_771 * (radius(8) - radius(9)) / (radius(8) - radius(10))
  NER(9) = NER_TOPDDOUBLEPRIME_771 - NER(8)
  NER(10) = NER_CMB_TOPDDOUBLEPRIME
  ! distributes NER_OUTER_CORE onto two element layer regions depending on vertical sizes of layers
  NER(11) = NER_OUTER_CORE * (radius(11) - radius(12)) / (radius(11) - radius(13))
  NER(12) = NER_OUTER_CORE - NER(11)
  NER(13) = NER_TOP_CENTRAL_CUBE_ICB

  ! debug
  !print *,'input NER:',NER(:)

  ! Find the Number of Radial Elements in a region based upon
  ! the aspect ratio of the elements
  call auto_optimal_ner(NUM_REGIONS, WIDTH, NEX_MAX, radius, scaling, NER, ratio_top, ratio_bottom)

  ! debug
  !print *,'output NER:',NER(:)

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

  R_CENTRAL_CUBE           = radius(14) * 1000.0d0

  end subroutine auto_ner

!
!-------------------------------------------------------------------------------------------------
!

  subroutine auto_optimal_ner(NUM_REGIONS, width, NEX, r, scaling, NER, rt, rb)

  use constants, only: DEGREES_TO_RADIANS

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

  ! Find optimal elements per region
  do i = 1,NUM_REGIONS-1
     dr = r(i) - r(i+1)              ! Radial Length of Region
     wt = width * DEGREES_TO_RADIANS * r(i)   / (NEX*1.0d0 / scaling(i)*1.0d0) ! Element Width Top
     wb = width * DEGREES_TO_RADIANS * r(i+1) / (NEX*1.0d0 / scaling(i)*1.0d0) ! Element Width Bottom
     w  = (wt + wb) * 0.5d0          ! Average Width of Region
     ner_test = NER(i)               ! Initial solution
     ratio = (dr / ner_test) / w     ! Aspect Ratio of Element
     xi = dabs(ratio - 1.0d0)        ! Aspect Ratio should be near 1.0
     ximin = 1.e7                    ! Initial Minimum

     !debug
     !print *,'region ',i,'element ratio: ',ratio,'xi = ',xi,'width = ',w

     ! increases NER to reach vertical/horizontal element ratio of about 1
     do while(xi <= ximin)
        NER(i) = ner_test            ! Found a better solution
        ximin = xi                   !
        ner_test = ner_test + 1      ! Increment ner_test and
        ratio = (dr / ner_test) / w  ! look for a better
        xi = dabs(ratio - 1.0d0)     ! solution
     enddo
     rt(i) = dr / NER(i) / wt        ! Find the Ratio of Top
     rb(i) = dr / NER(i) / wb        ! and Bottom for completeness

     !debug
     !print *,'region ',i,'element ratio: top = ',rt(i),'bottom = ',rb(i)
  enddo

  end subroutine auto_optimal_ner

!
!-------------------------------------------------------------------------------------------------
!

  subroutine find_r_central_cube(nex_xi_in, rcube)

  implicit none

  integer, parameter :: NBNODE = 8
  double precision, parameter :: alpha = 0.41d0

  integer,intent(in) :: nex_xi_in
  double precision,intent(out) :: rcube

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

  nex_xi = nex_xi_in / 16

  rcubestep    = 1.0d0
  rcube_test   =  930.0d0
  rcubemax     = 1100.0d0
  ximin        = 1e7
  rcube        = rcube_test

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
     xi = (max_edgemax / min_edgemin)
!       xi = abs(rcube_test - 981.0d0) / 45.0d0
!       write(*,'(a,5(f14.4,2x))')'rcube, xi, ximin:-',rcube_test, xi, min_edgemin,max_edgemax,max_aspect_ratio
     deallocate(points)

     if (xi < ximin) then
        ximin      = xi
        rcube      = rcube_test
     endif
     rcube_test = rcube_test + rcubestep
  enddo

  end subroutine find_r_central_cube

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_nex(nex_xi, rcube, alpha, ner)

  use constants, only: PI,PI_OVER_TWO,PI_OVER_FOUR

  implicit none

  double precision, parameter :: RICB_KM = 1221.0d0

  integer,intent(in) :: nex_xi
  double precision,intent(in) :: rcube, alpha
  integer,intent(out) :: ner

  ! local parameters
  integer :: ix
  double precision :: ratio_x, factx, xi
  double precision :: x, y
  double precision :: surfx, surfy
  double precision :: dist_cc_icb, somme, dist_moy

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

  implicit none

  double precision, parameter :: RICB_KM = 1221.0d0

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
