!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 5
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine read_parameter_file(MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_220_MOHO,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,NER_DOUBLING_OUTER_CORE, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
          NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS, &
          NUMBER_OF_THIS_RUN,NCHUNKS,DT,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC, &
          ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
          TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
          ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST, &
          ROTATION,ISOTROPIC_3D_MANTLE,TOPOGRAPHY,OCEANS,MOVIE_SURFACE, &
          MOVIE_VOLUME,ATTENUATION_3D,RECEIVERS_CAN_BE_BURIED, &
          PRINT_SOURCE_TIME_FUNCTION,SAVE_MESH_FILES, &
          ATTENUATION,IASPEI,ABSORBING_CONDITIONS, &
          INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,LOCAL_PATH,MODEL,SIMULATION_TYPE,SAVE_FORWARD)

  implicit none

  include "constants.h"

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_220_MOHO,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,NER_DOUBLING_OUTER_CORE, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
          NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS, &
          NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE

  double precision DT,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC, &
          ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE, &
          TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ATTENUATION,IASPEI, &
          ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD

  character(len=150) OUTPUT_FILES,LOCAL_PATH,MODEL,CMTSOLUTION

! local variables
  integer ios,icounter,isource,idummy,NEX_MAX

  double precision RECORD_LENGTH_IN_MINUTES,hdur,minval_hdur

  character(len=150) dummystring

  double precision ELEMENT_WIDTH

  integer, external :: err_occurred

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

  call open_parameter_file

  call read_value_integer(SIMULATION_TYPE, 'solver.SIMULATION_TYPE')
  if(err_occurred() /= 0) return
  call read_value_logical(SAVE_FORWARD, 'solver.SAVE_FORWARD')
  if(err_occurred() /= 0) return

  call read_value_integer(NCHUNKS, 'mesher.NCHUNKS')
  if(err_occurred() /= 0) return
  if(NCHUNKS /= 1 .and. NCHUNKS /= 2 .and. NCHUNKS /= 3 .and. NCHUNKS /= 6) stop 'NCHUNKS must be either 1, 2, 3 or 6'

  call read_value_double_precision(ANGULAR_WIDTH_XI_IN_DEGREES, 'mesher.ANGULAR_WIDTH_XI_IN_DEGREES')
  if(err_occurred() /= 0) return
  call read_value_double_precision(ANGULAR_WIDTH_ETA_IN_DEGREES, 'mesher.ANGULAR_WIDTH_ETA_IN_DEGREES')
  if(err_occurred() /= 0) return
  call read_value_double_precision(CENTER_LATITUDE_IN_DEGREES, 'mesher.CENTER_LATITUDE_IN_DEGREES')
  if(err_occurred() /= 0) return
  call read_value_double_precision(CENTER_LONGITUDE_IN_DEGREES, 'mesher.CENTER_LONGITUDE_IN_DEGREES')
  if(err_occurred() /= 0) return
  call read_value_double_precision(GAMMA_ROTATION_AZIMUTH, 'mesher.GAMMA_ROTATION_AZIMUTH')
  if(err_occurred() /= 0) return

! this MUST be 90 degrees for two chunks or more to match geometrically
  if(NCHUNKS > 1 .and. ANGULAR_WIDTH_XI_IN_DEGREES /= 90.d0) &
    stop 'ANGULAR_WIDTH_XI_IN_DEGREES must be 90 for more than one chunk'

! this can be any value in the case of two chunks
  if(NCHUNKS > 2 .and. ANGULAR_WIDTH_ETA_IN_DEGREES /= 90.d0) &
    stop 'ANGULAR_WIDTH_ETA_IN_DEGREES must be 90 for more than two chunks'

! include central cube or not
! use regular cubed sphere instead of cube for large distances
  if(NCHUNKS == 6) then
    INCLUDE_CENTRAL_CUBE = .true.
    INFLATE_CENTRAL_CUBE = .false.
  else
    INCLUDE_CENTRAL_CUBE = .false.
    INFLATE_CENTRAL_CUBE = .true.
  endif

! number of elements at the surface along the two sides of the first chunk
  call read_value_integer(NEX_XI, 'mesher.NEX_XI')
  if(err_occurred() /= 0) return
  call read_value_integer(NEX_ETA, 'mesher.NEX_ETA')
  if(err_occurred() /= 0) return
  call read_value_integer(NPROC_XI, 'mesher.NPROC_XI')
  if(err_occurred() /= 0) return
  call read_value_integer(NPROC_ETA, 'mesher.NPROC_ETA')
  if(err_occurred() /= 0) return

! set time step, radial distribution of elements, and attenuation period range
! right distribution is determined based upon maximum value of NEX
  NEX_MAX = max(NEX_XI,NEX_ETA)

! standard mesh on Caltech cluster
  if(NEX_MAX <= 160) then   !  Element Width = 0.5625 degrees =~ 62 km

    DT                       = 0.26d0

    MIN_ATTENUATION_PERIOD   = 20
    MAX_ATTENUATION_PERIOD   = 1000

    NER_CRUST                = 1
    NER_220_MOHO             = 3
    NER_400_220              = 2
    NER_600_400              = 2
    NER_670_600              = 1
    NER_771_670              = 1
    NER_TOPDDOUBLEPRIME_771  = 16
    NER_CMB_TOPDDOUBLEPRIME  = 1
    RATIO_TOP_DBL_OC         = 0.40d0
    RATIO_BOTTOM_DBL_OC      = 0.27d0
    NER_TOPDBL_CMB           = 12
    NER_ICB_BOTTOMDBL        = 6
    NER_TOP_CENTRAL_CUBE_ICB = 3

! Par_file_240_20sec_90x90
  else if(NEX_MAX <= 240) then ! Element Width = 0.3750 degrees =~ 41 km

    DT                       = 0.20d0

    MIN_ATTENUATION_PERIOD   = 20
    MAX_ATTENUATION_PERIOD   = 1000

    NER_CRUST                = 1
    NER_220_MOHO             = 3
    NER_400_220              = 2
    NER_600_400              = 2
    NER_670_600              = 1
    NER_771_670              = 1
    NER_TOPDDOUBLEPRIME_771  = 16
    NER_CMB_TOPDDOUBLEPRIME  = 1
    RATIO_TOP_DBL_OC         = 0.40d0
    RATIO_BOTTOM_DBL_OC      = 0.27d0
    NER_TOPDBL_CMB           = 12
    NER_ICB_BOTTOMDBL        = 6
    NER_TOP_CENTRAL_CUBE_ICB = 3

! Par_file_tsuboi_movie_320
  else if(NEX_MAX <= 320) then ! Element Width = 0.2812 degrees =~ 31 km

    DT                       = 0.10d0

    MIN_ATTENUATION_PERIOD   = 20
    MAX_ATTENUATION_PERIOD   = 1000

    NER_CRUST                = 2
    NER_220_MOHO             = 6
    NER_400_220              = 4
    NER_600_400              = 4
    NER_670_600              = 2
    NER_771_670              = 2
    NER_TOPDDOUBLEPRIME_771  = 32
    NER_CMB_TOPDDOUBLEPRIME  = 2
    RATIO_TOP_DBL_OC         = 0.40d0
    RATIO_BOTTOM_DBL_OC      = 0.27d0
    NER_TOPDBL_CMB           = 24
    NER_ICB_BOTTOMDBL        = 12
    NER_TOP_CENTRAL_CUBE_ICB = 6

! Par_file_480_10sec_90x90
  else if(NEX_MAX <= 480) then ! Element Width = 0.1875 degrees =~ 20 km

    DT                       = 0.10d0

    MIN_ATTENUATION_PERIOD   = 20
    MAX_ATTENUATION_PERIOD   = 1000

    NER_CRUST                = 2
    NER_220_MOHO             = 6
    NER_400_220              = 4
    NER_600_400              = 4
    NER_670_600              = 2
    NER_771_670              = 2
    NER_TOPDDOUBLEPRIME_771  = 32
    NER_CMB_TOPDDOUBLEPRIME  = 2
    RATIO_TOP_DBL_OC         = 0.40d0
    RATIO_BOTTOM_DBL_OC      = 0.27d0
    NER_TOPDBL_CMB           = 24
    NER_ICB_BOTTOMDBL        = 12
    NER_TOP_CENTRAL_CUBE_ICB = 6

! Par_file_ES_512_8.5sec
  else if(NEX_MAX <= 512) then ! Element Width = 0.1758 degrees =~ 19 km

    DT                       = 0.125d0

    MIN_ATTENUATION_PERIOD   = 8
    MAX_ATTENUATION_PERIOD   = 1000

    NER_CRUST                = 2
    NER_220_MOHO             = 5
    NER_400_220              = 4
    NER_600_400              = 4
    NER_670_600              = 2
    NER_771_670              = 2
    NER_TOPDDOUBLEPRIME_771  = 32
    NER_CMB_TOPDDOUBLEPRIME  = 3
    RATIO_TOP_DBL_OC         = 0.40d0
    RATIO_BOTTOM_DBL_OC      = 0.30d0
    NER_TOPDBL_CMB           = 22
    NER_ICB_BOTTOMDBL        = 12
    NER_TOP_CENTRAL_CUBE_ICB = 4

! Par_file_ES_640_6.75sec
  else if(NEX_MAX <= 640) then ! Element Width = 0.1460 degrees =~ 16 km

    DT                       = 0.125d0

    MIN_ATTENUATION_PERIOD   = 5
    MAX_ATTENUATION_PERIOD   = 400

    NER_CRUST                = 2
    NER_220_MOHO             = 5
    NER_400_220              = 4
    NER_600_400              = 4
    NER_670_600              = 2
    NER_771_670              = 2
    NER_TOPDDOUBLEPRIME_771  = 32
    NER_CMB_TOPDDOUBLEPRIME  = 3
    RATIO_TOP_DBL_OC         = 0.40d0
    RATIO_BOTTOM_DBL_OC      = 0.30d0
    NER_TOPDBL_CMB           = 22
    NER_ICB_BOTTOMDBL        = 12
    NER_TOP_CENTRAL_CUBE_ICB = 4

! Par_file_ES_1944procs_243nodes_5sec
  else if(NEX_MAX <= 864) then ! Element Width = 0.1042 degrees =~ 11 km

    DT                       = 0.05d0

    MIN_ATTENUATION_PERIOD   = 5
    MAX_ATTENUATION_PERIOD   = 400

    NER_CRUST                = 3
    NER_220_MOHO             = 7
    NER_400_220              = 7
    NER_600_400              = 7
    NER_670_600              = 3
    NER_771_670              = 3
    NER_TOPDDOUBLEPRIME_771  = 46
    NER_CMB_TOPDDOUBLEPRIME  = 4
    RATIO_TOP_DBL_OC         = 0.48d0
    RATIO_BOTTOM_DBL_OC      = 0.43d0
    NER_TOPDBL_CMB           = 38
    NER_ICB_BOTTOMDBL        = 30
    NER_TOP_CENTRAL_CUBE_ICB = 6

! Par_file_ES_1152_4sec
  else if(NEX_MAX <= 1152) then ! Element Width = 0.0781 =~ 8.5 km

    DT                       = 0.0555555555d0

    MIN_ATTENUATION_PERIOD   = 4
    MAX_ATTENUATION_PERIOD   = 400

    NER_CRUST                = 4
    NER_220_MOHO             = 8
    NER_400_220              = 7
    NER_600_400              = 7
    NER_670_600              = 3
    NER_771_670              = 4
    NER_TOPDDOUBLEPRIME_771  = 62
    NER_CMB_TOPDDOUBLEPRIME  = 5
    RATIO_TOP_DBL_OC         = 0.48d0
    RATIO_BOTTOM_DBL_OC      = 0.43d0
    NER_TOPDBL_CMB           = 38
    NER_ICB_BOTTOMDBL        = 30
    NER_TOP_CENTRAL_CUBE_ICB = 6

! Par_file_ES_4056procs_507nodes_3.5sec
  else if(NEX_MAX <= 1248) then ! Element Width = 0.0721 =~ 8 km

    DT                       = 0.05d0

    MIN_ATTENUATION_PERIOD   = 4
    MAX_ATTENUATION_PERIOD   = 400

    NER_CRUST                = 3
    NER_220_MOHO             = 9
    NER_400_220              = 9
    NER_600_400              = 9
    NER_670_600              = 4
    NER_771_670              = 4
    NER_TOPDDOUBLEPRIME_771  = 60
    NER_CMB_TOPDDOUBLEPRIME  = 5
    RATIO_TOP_DBL_OC         = 0.44d0
    RATIO_BOTTOM_DBL_OC      = 0.41d0
    NER_TOPDBL_CMB           = 50
    NER_ICB_BOTTOMDBL        = 40
    NER_TOP_CENTRAL_CUBE_ICB = 8

!  else
!    stop 'this value of NEX_MAX is not in the database, edit read_parameter_file.f90 and recompile'
  endif

  if(ANGULAR_WIDTH_XI_IN_DEGREES < 90.0d0 .OR. NEX_MAX > 1248) then

    call auto_ner(ANGULAR_WIDTH_XI_IN_DEGREES, NEX_MAX, &
          NER_CRUST, NER_220_MOHO, NER_400_220, NER_600_400, &
          NER_670_600, NER_771_670, NER_TOPDDOUBLEPRIME_771, &
          NER_CMB_TOPDDOUBLEPRIME, RATIO_TOP_DBL_OC, RATIO_BOTTOM_DBL_OC, &
          NER_TOPDBL_CMB, NER_ICB_BOTTOMDBL, NER_TOP_CENTRAL_CUBE_ICB)

    ! Determine in appropriate period range for the current mesh
    !   Note: we are using DT as a temporary variable
    ! The Minimum attenuation period = (Grid Spacing in km) / V_min
    !  Grid spacing in km     = Width of an element in km * spacing for GLL point * points per wavelength
    !  Width of element in km = (Angular width in degrees / NEX_MAX) * degrees to km
    !  degrees to km          = 111.11d0
    !  spacing for GLL point  = 4
    !  points per wavelength  = 4
    DT = (ANGULAR_WIDTH_XI_IN_DEGREES / (4.0d0 * dble(NEX_MAX)) * 111.11d0 * 4.0d0 ) / 2.25d0
    MIN_ATTENUATION_PERIOD = DT

    ! The max attenuation period for 3 SLS is optimally
    !   1.75 decades from the min attenuation period
    DT = dble(MIN_ATTENUATION_PERIOD) * 10.0d0**1.75d0
    MAX_ATTENUATION_PERIOD = DT

    ! 0.173 is the minimum spacing between GLL points for NGLL = 5
    ! This should be changed in the future, placed in a header file
    ! 0.40 is the ratio of radial lengths of elements inside the
    ! central cube to those just outside the central cube
    ! 1221.0 is the Radius of the inncer core in km
    ! 0.40 is the maximum stability condition
    ! 11.02827 is Vp near the inner core boundary
    ! See equation 48 in Komatitsch and Tromp (2002, Part I)
    DT = (0.40d0 * &
         ((ANGULAR_WIDTH_XI_IN_DEGREES * (PI/180.0d0)) * 1221.0d0) / &
         (dble(NEX_MAX) / 8.0d0) / 11.02827d0 ) * 0.173d0 * 0.4d0

  endif

! scale radial mesh parameters according to definitions used in mesher
! in order to implement mesh doubling
  NER_220_MOHO = NER_220_MOHO * 2
  NER_400_220 = NER_400_220 * 2
  NER_600_400 = NER_600_400 * 2
  NER_670_600 = NER_670_600 * 2

  NER_771_670 = NER_771_670 * 4
  NER_TOPDDOUBLEPRIME_771 = NER_TOPDDOUBLEPRIME_771 * 4
  NER_CMB_TOPDDOUBLEPRIME = NER_CMB_TOPDDOUBLEPRIME * 4
  NER_TOPDBL_CMB = NER_TOPDBL_CMB * 4
  NER_ICB_BOTTOMDBL = NER_ICB_BOTTOMDBL * 4
  NER_TOP_CENTRAL_CUBE_ICB = NER_TOP_CENTRAL_CUBE_ICB * 4

  NER_ICB_CMB = NER_ICB_BOTTOMDBL + NER_BOTTOMDBL_TOPDBL + NER_TOPDBL_CMB
  NER_DOUBLING_OUTER_CORE = NER_TOP_CENTRAL_CUBE_ICB + NER_ICB_BOTTOMDBL + NER_BOTTOMDBL_TOPDBL

! cannot honor middle crust if only one spectral element in radial direction
  if(NER_CRUST == 1) then
    ONE_CRUST = .true.
  else
    ONE_CRUST = .false.
  endif

! define the velocity model
  call read_value_string(MODEL, 'MODEL')
  if(err_occurred() /= 0) return

  if(MODEL == 'isotropic_prem') then
    IASPEI = .false.
    TRANSVERSE_ISOTROPY = .false.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.

  else if(MODEL == 'transversly_isotropic_prem') then
    IASPEI = .false.
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.

  else if(MODEL == 'iaspei') then
    IASPEI = .true.
    TRANSVERSE_ISOTROPY = .false.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.

  else if(MODEL == 's20rts') then
    IASPEI = .false.
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .true.
    ATTENUATION_3D = .false.

  else if(MODEL == 'Brian_Savage') then
    IASPEI = .false.
    TRANSVERSE_ISOTROPY = .false.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .true.

  else if(MODEL == 'Min_Chen') then
    IASPEI = .false.
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .true.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.

  else
    stop 'model not implemented, edit read_parameter_file.f90 and recompile'
  endif

  call read_value_logical(OCEANS, 'model.OCEANS')
  if(err_occurred() /= 0) return
  call read_value_logical(ELLIPTICITY, 'model.ELLIPTICITY')
  if(err_occurred() /= 0) return
  call read_value_logical(TOPOGRAPHY, 'model.TOPOGRAPHY')
  if(err_occurred() /= 0) return
  call read_value_logical(GRAVITY, 'model.GRAVITY')
  if(err_occurred() /= 0) return
  call read_value_logical(ROTATION, 'model.ROTATION')
  if(err_occurred() /= 0) return
  call read_value_logical(ATTENUATION, 'model.ATTENUATION')
  if(err_occurred() /= 0) return

  call read_value_logical(ABSORBING_CONDITIONS, 'solver.ABSORBING_CONDITIONS')
  if(err_occurred() /= 0) return

  if(ABSORBING_CONDITIONS .and. NCHUNKS == 6) stop 'cannot have absorbing conditions in the full Earth'

  if(ABSORBING_CONDITIONS .and. NCHUNKS == 3) stop 'absorbing conditions not supported for three chunks yet'

  if(ATTENUATION_3D .and. .not. ATTENUATION) stop 'need ATTENUATION to use ATTENUATION_3D'

! radii in PREM or IASPEI
! and normalized density at fluid-solid interface on fluid size for coupling
! ROCEAN: radius of the ocean (m)
! RMIDDLE_CRUST: radius of the middle crust (m)
! RMOHO: radius of the Moho (m)
! R80: radius of 80 km discontinuity (m)
! R220: radius of 220 km discontinuity (m)
! R400: radius of 400 km discontinuity (m)
! R600: radius of 600 km 2nd order discontinuity (m)
! R670: radius of 670 km discontinuity (m)
! R771: radius of 771 km 2nd order discontinuity (m)
! RTOPDDOUBLEPRIME: radius of top of D" 2nd order discontinuity (m)
! RCMB: radius of CMB (m)
! RICB: radius of ICB (m)

! values common to PREM and IASPEI
  ROCEAN = 6368000.d0
  RMIDDLE_CRUST = 6356000.d0
  R80 = 6291000.d0

  RHO_OCEANS = 1020.0 / RHOAV

  if(IASPEI) then

! IASPEI
    RMOHO = 6341000.d0
    R220 = 6161000.d0
    R400 = 5961000.d0
    R600 = 5781000.d0
!    R670 = 5711000.d0
    R670 = 5701000.d0
    R771 = 5611000.d0
    RTOPDDOUBLEPRIME = 3631000.d0
    RCMB = 3482000.d0
    RICB = 1217000.d0

    RHO_TOP_OC = 9900.2379 / RHOAV
    RHO_BOTTOM_OC = 12168.6383 / RHOAV

  else

! PREM
    RMOHO = 6346600.d0
    R220 = 6151000.d0
    R400 = 5971000.d0
    R600 = 5771000.d0
    R670 = 5701000.d0
    R771 = 5600000.d0
    RTOPDDOUBLEPRIME = 3630000.d0
    RCMB = 3480000.d0
    RICB = 1221000.d0

    RHO_TOP_OC = 9903.4384 / RHOAV
    RHO_BOTTOM_OC = 12166.5885 / RHOAV

  endif

! non-dimensionalized size of central cube in the inner core
! This is where the central cube in the inner core and the rest of the mesh
! are matched (150 km below the ICB is optimal)
  R_CENTRAL_CUBE = (RICB - 150000.d0) / R_EARTH

  call read_value_double_precision(RECORD_LENGTH_IN_MINUTES, 'solver.RECORD_LENGTH_IN_MINUTES')
  if(err_occurred() /= 0) return

! compute total number of time steps, rounded to next multiple of 100
  NSTEP = 100 * (int(RECORD_LENGTH_IN_MINUTES * 60.d0 / (100.d0*DT)) + 1)

! compute the total number of sources in the CMTSOLUTION file
! there are NLINES_PER_CMTSOLUTION_SOURCE lines per source in that file
  call get_value_string(CMTSOLUTION, 'solver.CMTSOLUTION', 'DATA/CMTSOLUTION')
  open(unit=1,file=CMTSOLUTION,iostat=ios,status='old')
  if(ios /= 0) stop 'error opening CMTSOLUTION file'
  icounter = 0
  do while(ios == 0)
    read(1,"(a)",iostat=ios) dummystring
    if(ios == 0) icounter = icounter + 1
  enddo
  close(1)
  if(mod(icounter,NLINES_PER_CMTSOLUTION_SOURCE) /= 0) &
    stop 'total number of lines in CMTSOLUTION file should be a multiple of NLINES_PER_CMTSOLUTION_SOURCE'
  NSOURCES = icounter / NLINES_PER_CMTSOLUTION_SOURCE
  if(NSOURCES < 1) stop 'need at least one source in CMTSOLUTION file'

  call read_value_logical(MOVIE_SURFACE, 'solver.MOVIE_SURFACE')
  if(err_occurred() /= 0) return
  call read_value_logical(MOVIE_VOLUME, 'solver.MOVIE_VOLUME')
  if(err_occurred() /= 0) return
  call read_value_integer(NTSTEP_BETWEEN_FRAMES, 'solver.NTSTEP_BETWEEN_FRAMES')
  if(err_occurred() /= 0) return
  call read_value_double_precision(HDUR_MOVIE, 'solver.HDUR_MOVIE')
  if(err_occurred() /= 0) return

! computes a default hdur_movie that creates nice looking movies.
! Sets HDUR_MOVIE as the minimum period the mesh can resolve
  if(HDUR_MOVIE <=TINYVAL) &
    HDUR_MOVIE = 1.1d0*max(240.d0/NEX_XI*18.d0*ANGULAR_WIDTH_XI_IN_DEGREES/90.d0, &
                           240.d0/NEX_ETA*18.d0*ANGULAR_WIDTH_ETA_IN_DEGREES/90.d0)


! compute the minimum value of hdur in CMTSOLUTION file
  open(unit=1,file=CMTSOLUTION,status='old')
  minval_hdur = HUGEVAL
  do isource = 1,NSOURCES

! skip other information
    do idummy = 1,3
      read(1,"(a)") dummystring
    enddo

! read half duration and compute minimum
    read(1,"(a)") dummystring
    read(dummystring(15:len_trim(dummystring)),*) hdur
    minval_hdur = min(minval_hdur,hdur)

! skip other information
    do idummy = 1,9
      read(1,"(a)") dummystring
    enddo

  enddo
  close(1)

! one cannot use a Heaviside source for the movies
!  if((MOVIE_SURFACE .or. MOVIE_VOLUME) .and. minval_hdur < TINYVAL) &
!    stop 'hdur too small for movie creation, movies do not make sense for Heaviside source'

  call read_value_logical(SAVE_MESH_FILES, 'mesher.SAVE_MESH_FILES')
  if(err_occurred() /= 0) return
  call read_value_integer(NUMBER_OF_RUNS, 'solver.NUMBER_OF_RUNS')
  if(err_occurred() /= 0) return
  call read_value_integer(NUMBER_OF_THIS_RUN, 'solver.NUMBER_OF_THIS_RUN')
  if(err_occurred() /= 0) return
  call read_value_string(LOCAL_PATH, 'LOCAL_PATH')
  if(err_occurred() /= 0) return
  call read_value_integer(NTSTEP_BETWEEN_OUTPUT_INFO, 'solver.NTSTEP_BETWEEN_OUTPUT_INFO')
  if(err_occurred() /= 0) return
  call read_value_integer(NTSTEP_BETWEEN_OUTPUT_SEISMOS, 'solver.NTSTEP_BETWEEN_OUTPUT_SEISMOS')
  if(err_occurred() /= 0) return
  call read_value_logical(RECEIVERS_CAN_BE_BURIED, 'solver.RECEIVERS_CAN_BE_BURIED')
  if(err_occurred() /= 0) return
  call read_value_logical(PRINT_SOURCE_TIME_FUNCTION, 'solver.PRINT_SOURCE_TIME_FUNCTION')
  if(err_occurred() /= 0) return

! close parameter file
  call close_parameter_file

!--- check that parameters make sense

! check that reals are either 4 or 8 bytes
  if(CUSTOM_REAL /= SIZE_REAL .and. CUSTOM_REAL /= SIZE_DOUBLE) stop 'wrong size of CUSTOM_REAL for reals'

! check that the parameter file is correct
  if(NGNOD /= 27) stop 'number of control nodes must be 27'
  if(NGNOD == 27 .and. NGNOD2D /= 9) stop 'elements with 27 points should have NGNOD2D = 9'

! for the number of standard linear solids for attenuation
  if(N_SLS /= 3) stop 'number of SLS must be 3'

! check number of slices in each direction
  if(NCHUNKS < 1) stop 'must have at least one chunk'
  if(NPROC_XI < 1) stop 'NPROC_XI must be at least 1'
  if(NPROC_ETA < 1) stop 'NPROC_ETA must be at least 1'

! check number of chunks
  if(NCHUNKS /= 1 .and. NCHUNKS /= 2 .and. NCHUNKS /= 3 .and. NCHUNKS /= 6) &
     stop 'only one, two, three or six chunks can be meshed'

! check that the central cube can be included
  if(INCLUDE_CENTRAL_CUBE .and. NCHUNKS /= 6) stop 'need six chunks to include central cube'

! check that topology is correct if more than two chunks
  if(NCHUNKS > 2 .and. NPROC_XI /= NPROC_ETA) stop 'must have NPROC_XI=NPROC_ETA for more than two chunks'

! check that size can be coarsened in depth twice (block size multiple of 8)
  if(mod(NEX_XI/8,NPROC_XI) /= 0) stop 'NEX_XI must be a multiple of 8*NPROC_XI'

  if(mod(NEX_ETA/8,NPROC_ETA) /= 0) stop 'NEX_ETA must be a multiple of 8*NPROC_ETA'

  if(mod(NEX_XI,8) /= 0) stop 'NEX_XI must be a multiple of 8'

  if(mod(NEX_ETA,8) /= 0) stop 'NEX_ETA must be a multiple of 8'

! check that sphere can be cut into slices without getting negative Jacobian
  if(NEX_XI < 48) stop 'NEX_XI must be greater than 48 to cut the sphere into slices with positive Jacobian'
  if(NEX_ETA < 48) stop 'NEX_ETA must be greater than 48 to cut the sphere into slices with positive Jacobian'

! check that mesh doubling can be implemented
  if(NER_220_MOHO/2<3) stop 'NER_220_MOHO should be at least 3'
  if(NER_400_220/2<2) stop 'NER_400_220 should be at least 2'

! check that IASPEI is isotropic
  if(IASPEI .and. TRANSVERSE_ISOTROPY) stop 'IASPEI is currently isotropic'

  ELEMENT_WIDTH = ANGULAR_WIDTH_XI_IN_DEGREES/dble(NEX_MAX) * DEGREES_TO_RADIANS

  end subroutine read_parameter_file

  subroutine auto_ner(WIDTH, NEX_MAX, &
       NER_CRUST, NER_220_MOHO, NER_400_220, NER_600_400, &
       NER_670_600, NER_771_670, NER_TOPDDOUBLEPRIME_771, &
       NER_CMB_TOPDDOUBLEPRIME, RATIO_TOP_DBL_OC, RATIO_BOTTOM_DBL_OC, &
       NER_TOPDBL_CMB, NER_ICB_BOTTOMDBL, NER_TOP_CENTRAL_CUBE_ICB)

    implicit none

    include 'constants.h'

    double precision WIDTH
    integer NEX_MAX
    integer NER_CRUST, NER_220_MOHO, NER_400_220, NER_600_400, &
         NER_670_600, NER_771_670, NER_TOPDDOUBLEPRIME_771, &
         NER_CMB_TOPDDOUBLEPRIME, NER_TOPDBL_CMB, NER_ICB_BOTTOMDBL, &
         NER_TOP_CENTRAL_CUBE_ICB
    double precision RATIO_TOP_DBL_OC, RATIO_BOTTOM_DBL_OC

    integer, parameter                         :: NUM_REGIONS = 13
    integer, dimension(NUM_REGIONS)            :: scaling
    double precision, dimension(NUM_REGIONS)   :: radius
    double precision, dimension(NUM_REGIONS)   :: element_width
    double precision, dimension(NUM_REGIONS)   :: chunk_width
    double precision, dimension(NUM_REGIONS-1) :: ratio_top
    double precision, dimension(NUM_REGIONS-1) :: ratio_bottom
    integer, dimension(NUM_REGIONS-1)          :: NER
    integer NER_FLUID

    ! This is PREM in Kilometers
    radius(1)  = 6371.00d0 ! Surface
    radius(2)  = 6346.60d0 ! Moho
    radius(3)  = 6151.00d0 ! 220
    radius(4)  = 5971.00d0 ! 400
    radius(5)  = 5771.00d0 ! 600
    radius(6)  = 5701.00d0 ! 670
    radius(7)  = 5600.00d0 ! 771
    radius(8)  = 3630.00d0 ! D''
    radius(9)  = 3480.00d0 ! CMB
    radius(10) =    0.00d0 ! Top Double Fluid
    radius(11) =    0.00d0 ! Bottom Double Fluid
    radius(12) = 1221.00d0 ! ICB
    radius(13) = 1071.00d0 ! Top Central Cube

    ! Mesh Doubling
    scaling(1:1)   = 1
    scaling(2:5)   = 2
    scaling(6:11)  = 4
    scaling(12:13) = 8

    ! Minimum Number of Elements a Region must have
    NER(:)  = 1
    NER(2)  = 3
    NER(3)  = 2
    NER(12) = 2

    NER_FLUID = 6


    ! Find the Optimal Height of the Fluid Region based on the
    ! Aspect ratio of elements within the fluid and the total
    ! number of elements within the fluid
    call auto_fluid_double(WIDTH, NEX_MAX, NUM_REGIONS, radius, scaling, NER_FLUID, RATIO_TOP_DBL_OC,  RATIO_BOTTOM_DBL_OC)

    ! Determine the Radius of Top and Bottom of Fluid Doubling Region
    radius(10) = radius(12) + RATIO_TOP_DBL_OC    * (radius(9) - radius(12))
    radius(11) = radius(12) + RATIO_BOTTOM_DBL_OC * (radius(9) - radius(12))

    ! Horizontal Width of a Chunk
    chunk_width(:) = WIDTH * (PI/180.0d0) * radius(:)

    ! Horizontal Width of the elements within the chunk
    element_width(:) = chunk_width(:) / (NEX_MAX / scaling(:))


    ! Find the Number of Radial Elements in a region based upon
    ! the aspect ratio of the elements
    call auto_optimal_ner(NUM_REGIONS, radius, element_width,NER, ratio_top, ratio_bottom)

    ! Set Output arguments
    NER_CRUST                = NER(1)
    NER_220_MOHO             = NER(2)
    NER_400_220              = NER(3)
    NER_600_400              = NER(4)
    NER_670_600              = NER(5)
    NER_771_670              = NER(6)
    NER_TOPDDOUBLEPRIME_771  = NER(7)
    NER_CMB_TOPDDOUBLEPRIME  = NER(8)
    NER_TOPDBL_CMB           = NER(9)
    NER_FLUID                = NER(10)
    NER_ICB_BOTTOMDBL        = NER(11)
    NER_TOP_CENTRAL_CUBE_ICB = NER(12)

  end subroutine auto_ner

  subroutine auto_optimal_ner(NUM_REGIONS, r, ew, NER, rt, rb)
    implicit none

    integer NUM_REGIONS
    integer,          dimension(NUM_REGIONS-1) :: NER ! Elements per Region
    double precision, dimension(NUM_REGIONS)   :: r   ! Radius
    double precision, dimension(NUM_REGIONS)   :: ew  ! Element Width
    double precision, dimension(NUM_REGIONS-1) :: rt  ! Ratio at Top
    double precision, dimension(NUM_REGIONS-1) :: rb  ! Ratio at Bottom

    double precision dr, w, ratio, xi, ximin
    integer ner_test
    integer i

    ! Find optimal elements per region
    do i = 1,NUM_REGIONS-1
       dr = r(i) - r(i+1)              ! Radial Length of Ragion
       w  = (ew(i) + ew(i+1)) / 2.0d0  ! Average Width of Region
       ner_test = NER(i)               ! Initial solution
       ratio = (dr / ner_test) / w     ! Aspect Ratio of Element
       xi = dabs(ratio - 1.0d0)        ! Aspect Ratio should be near 1.0
       ximin = 1e7                     ! Initial Minimum

       do while(xi <= ximin)
          NER(i) = ner_test            ! Found a better solution
          ximin = xi                   !
          ner_test = ner_test + 1      ! Increment ner_test and
          ratio = (dr / ner_test) / w  ! look for a better
          xi = dabs(ratio - 1.0d0)     ! solution
       end do
       rt(i) = dr / NER(i) / ew(i)     ! Find the Ratio of Top
       rb(i) = dr / NER(i) / ew(i+1)   ! and Bottom for completeness
    end do

  end subroutine auto_optimal_ner


  subroutine auto_fluid_double(WIDTH, NEX_MAX, NUM_REGIONS, r, s, &
       NER_FLUID, RATIO_TOP_DBL_OC, RATIO_BOTTOM_DBL_OC)
    implicit none
    include 'constants.h'

    double precision WIDTH
    integer          NEX_MAX
    integer          NUM_REGIONS
    integer          NER_FLUID
    double precision RATIO_TOP_DBL_OC
    double precision RATIO_BOTTOM_DBL_OC
    double precision, dimension(NUM_REGIONS)  :: r  ! Radius
    integer,          dimension(NUM_REGIONS)  :: s  ! Mesh Scaling

    integer i, j
    double precision r1, r2, fluid_radius
    double precision rtop, rbot, wtop, wbot
    double precision wave, xi, ximin

    ! Find width of Fluid region
    ximin = 1.d7  ! Initial Minimum

    do i = 1,91
       do j = 1,91
          ! 0.05 <= R(1,2) <= 0.96
          r1 = 0.05d0 + i * 0.01d0
          r2 = 0.05d0 + j * 0.01d0

          ! R2 is defined to be less than R1 by definition
          if(r2 < r1) then
             ! Radii ( top, bottom, and element_radius)
             rtop = r(12) + r1 * (r(9) - r(12))                  ! Top
             rbot = r(12) + r2 * (r(9) - r(12))                  ! Bottom
             fluid_radius = (rtop - rbot) / NER_FLUID            ! Radius Element

             ! Element Widths ( top, bottom and average )
             wtop = (WIDTH * (PI/180) * rtop) / (NEX_MAX/s(10))  ! Top
             wbot = (WIDTH * (PI/180) * rbot) / (NEX_MAX/s(11))  ! Bottom
             wave = (wtop + wbot) / 2                            ! Average

             ! Aspect Ration should be near 1.0 and
             ! Centered around the middle of the Fluid Outer Core (Ratio = 0.5)
             xi = dabs(fluid_radius/wave - 1.0d0) + dabs((r1 + r2)/2.0d0 - 0.50d0)

             if(xi < ximin) then
                ! Set our current best solution
                ximin = xi
                RATIO_TOP_DBL_OC    = r1
                RATIO_BOTTOM_DBL_OC = r2
             end if
          end if
       enddo
    enddo
  end subroutine auto_fluid_double

