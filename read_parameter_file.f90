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
          NER_TOP_CENTRAL_CUBE_ICB,NEX_ETA,NEX_XI,NER_DOUBLING_OUTER_CORE, &
          NPROC_ETA,NPROC_XI,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
          NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS, &
          NUMBER_OF_THIS_RUN,NCHUNKS,DT,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC, &
          ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS, &
          TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
          ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST, &
          ROTATION,ISOTROPIC_3D_MANTLE,TOPOGRAPHY,OCEANS,MOVIE_SURFACE, &
          MOVIE_VOLUME,ATTENUATION_3D,RECEIVERS_CAN_BE_BURIED, &
          PRINT_SOURCE_TIME_FUNCTION,SAVE_AVS_DX_MESH_FILES, &
          ATTENUATION,IASPEI,ABSORBING_CONDITIONS, &
          INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,LOCAL_PATH,MODEL)

  implicit none

  include "constants.h"

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_220_MOHO,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_ETA,NEX_XI,NER_DOUBLING_OUTER_CORE, &
          NPROC_ETA,NPROC_XI,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
          NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS, &
          NUMBER_OF_THIS_RUN,NCHUNKS

  double precision DT,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC, &
          ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE, &
          TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_AVS_DX_MESH_FILES,ATTENUATION,IASPEI, &
          ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE

  character(len=150) LOCAL_PATH,MODEL

! first 34 characters of each line in the file are a comment
  character(len=34) junk

! local variables
  integer ios,icounter,isource,idummy,NEX_MAX

  double precision RECORD_LENGTH_IN_MINUTES,hdur,minval_hdur

  character(len=150) dummystring

  open(unit=IIN,file='DATA/Par_file',status='old')

! ignore header
  do idummy=1,11
    read(IIN,*)
  enddo

  read(IIN,1) junk,NCHUNKS
  if(NCHUNKS /= 1 .and. NCHUNKS /= 2 .and. NCHUNKS /= 3 .and. NCHUNKS /= 6) stop 'NCHUNKS must be either 1, 2, 3 or 6'

  read(IIN,*)
  read(IIN,*)
  read(IIN,2) junk,ANGULAR_WIDTH_XI_IN_DEGREES
  read(IIN,2) junk,ANGULAR_WIDTH_ETA_IN_DEGREES
  read(IIN,2) junk,CENTER_LATITUDE_IN_DEGREES
  read(IIN,2) junk,CENTER_LONGITUDE_IN_DEGREES
  read(IIN,2) junk,GAMMA_ROTATION_AZIMUTH

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
  read(IIN,*)
  read(IIN,*)
  read(IIN,*)
  read(IIN,1) junk,NEX_XI
  read(IIN,1) junk,NEX_ETA

  read(IIN,*)
  read(IIN,*)
  read(IIN,1) junk,NPROC_XI
  read(IIN,1) junk,NPROC_ETA

! set time step, radial distribution of elements, and attenuation period range
! right distribution is determined based upon maximum value of NEX
  NEX_MAX = max(NEX_XI,NEX_ETA)

! standard mesh on Caltech cluster
  if(NEX_MAX <= 160) then

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
  else if(NEX_MAX <= 240) then

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
  else if(NEX_MAX <= 320) then

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
  else if(NEX_MAX <= 480) then

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
  else if(NEX_MAX <= 512) then

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
  else if(NEX_MAX <= 640) then

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
  else if(NEX_MAX <= 864) then

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
  else if(NEX_MAX <= 1152) then

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
  else if(NEX_MAX <= 1248) then

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

  else
    stop 'this value of NEX_MAX is not in the database, edit read_parameter_file.f90 and recompile'
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
  read(IIN,*)
  read(IIN,*)
  read(IIN,*)
  read(IIN,4) junk,MODEL

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
    ISOTROPIC_3D_MANTLE = .true.
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

  read(IIN,*)
  read(IIN,*)
  read(IIN,3) junk,OCEANS
  read(IIN,3) junk,ELLIPTICITY
  read(IIN,3) junk,TOPOGRAPHY
  read(IIN,3) junk,GRAVITY
  read(IIN,3) junk,ROTATION
  read(IIN,3) junk,ATTENUATION

  read(IIN,*)
  read(IIN,*)
  read(IIN,3) junk,ABSORBING_CONDITIONS

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
    R670 = 5711000.d0
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

  read(IIN,*)
  read(IIN,*)
  read(IIN,2) junk,RECORD_LENGTH_IN_MINUTES

! compute total number of time steps, rounded to next multiple of 100
  NSTEP = 100 * (int(RECORD_LENGTH_IN_MINUTES * 60.d0 / (100.d0*DT)) + 1)

! compute the total number of sources in the CMTSOLUTION file
! there are NLINES_PER_CMTSOLUTION_SOURCE lines per source in that file
  open(unit=1,file='DATA/CMTSOLUTION',iostat=ios,status='old')
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

  read(IIN,*)
  read(IIN,*)
  read(IIN,3) junk,MOVIE_SURFACE
  read(IIN,3) junk,MOVIE_VOLUME
  read(IIN,1) junk,NTSTEP_BETWEEN_FRAMES

! compute the minimum value of hdur in CMTSOLUTION file
  open(unit=1,file='DATA/CMTSOLUTION',status='old')
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
  if((MOVIE_SURFACE .or. MOVIE_VOLUME) .and. minval_hdur < TINYVAL) &
    stop 'hdur too small for movie creation, movies do not make sense for Heaviside source'

  read(IIN,*)
  read(IIN,*)
  read(IIN,3) junk,SAVE_AVS_DX_MESH_FILES

  read(IIN,*)
  read(IIN,*)
  read(IIN,1) junk,NUMBER_OF_RUNS
  read(IIN,1) junk,NUMBER_OF_THIS_RUN

  read(IIN,*)
  read(IIN,*)
  read(IIN,4) junk,LOCAL_PATH

! ignore name of machine file (used by scripts but not by mesher nor solver)
  read(IIN,*)
  read(IIN,*)
  read(IIN,*)

  read(IIN,*)
  read(IIN,*)
  read(IIN,1) junk,NTSTEP_BETWEEN_OUTPUT_INFO

  read(IIN,*)
  read(IIN,*)
  read(IIN,1) junk,NTSTEP_BETWEEN_OUTPUT_SEISMOS

  read(IIN,*)
  read(IIN,*)
  read(IIN,3) junk,RECEIVERS_CAN_BE_BURIED

  read(IIN,*)
  read(IIN,*)
  read(IIN,3) junk,PRINT_SOURCE_TIME_FUNCTION

! close parameter file
  close(IIN)

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

! formats
 1 format(a,i20)
 2 format(a,f20.8)
 3 format(a,l20)
 4 format(a,a)

  end subroutine read_parameter_file

