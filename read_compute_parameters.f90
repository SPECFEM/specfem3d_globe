!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                    and University of Pau, France
! (c) California Institute of Technology and University of Pau, April 2007
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine read_compute_parameters(MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
         NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
         NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
         NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,RMOHO_FICTITIOUS_IN_MESHER, &
         NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
         NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
         NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,DT, &
         ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
         CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
         RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
         R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
         TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
         ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST, &
         ROTATION,ISOTROPIC_3D_MANTLE,TOPOGRAPHY,OCEANS,MOVIE_SURFACE, &
         MOVIE_VOLUME,ATTENUATION_3D,RECEIVERS_CAN_BE_BURIED, &
         PRINT_SOURCE_TIME_FUNCTION,SAVE_MESH_FILES, &
         ATTENUATION,REFERENCE_1D_MODEL,THREE_D_MODEL,ABSORBING_CONDITIONS, &
         INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,LOCAL_PATH,MODEL,SIMULATION_TYPE,SAVE_FORWARD, &
         NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
         NSPEC, &
         NSPEC2D_XI, &
         NSPEC2D_ETA, &
         NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
         NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
         NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
         NGLOB, &
         ratio_sampling_array, ner, doubling_index,r_bottom,r_top,this_region_has_a_doubling,rmins,rmaxs,CASE_3D, &
         OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
         ROTATE_SEISMOGRAMS_RT)


  implicit none

  include "constants.h"

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,RMOHO_FICTITIOUS_IN_MESHER, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE, &
          REFERENCE_1D_MODEL,THREE_D_MODEL

  double precision DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE, &
          TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ATTENUATION, &
          ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
          ROTATE_SEISMOGRAMS_RT

  character(len=150) OUTPUT_FILES,LOCAL_PATH,MODEL

! local variables
  integer NEX_MAX

  double precision RECORD_LENGTH_IN_MINUTES,ELEMENT_WIDTH

  integer, external :: err_occurred

! parameters to be computed based upon parameters above read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC, &
      NSPEC2D_XI, &
      NSPEC2D_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
      NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
      NGLOB

  integer nblocks_xi,nblocks_eta

  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling

  integer :: ielem,elem_doubling_mantle,elem_doubling_middle_outer_core,elem_doubling_bottom_outer_core
  double precision :: DEPTH_SECOND_DOUBLING_REAL,DEPTH_THIRD_DOUBLING_REAL, &
                          DEPTH_FOURTH_DOUBLING_REAL,distance,distance_min,zval

! honor PREM Moho or not
! doing so drastically reduces the stability condition and therefore the time step
  logical :: HONOR_1D_SPHERICAL_MOHO,CASE_3D

  integer :: ifirst_region, ilast_region, iter_region, iter_layer, doubling, padding, tmp_sum, tmp_sum_xi, tmp_sum_eta
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs
  integer ::  NUMBER_OF_MESH_LAYERS,layer_offset,nspec2D_xi_sb,nspec2D_eta_sb, &
              nb_lay_sb, nspec_sb, nglob_vol, nglob_surf, nglob_edge

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

  call open_parameter_file

  call read_value_integer(SIMULATION_TYPE, 'solver.SIMULATION_TYPE')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(SAVE_FORWARD, 'solver.SAVE_FORWARD')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

  call read_value_integer(NCHUNKS, 'mesher.NCHUNKS')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  if(NCHUNKS /= 1 .and. NCHUNKS /= 2 .and. NCHUNKS /= 3 .and. NCHUNKS /= 6) &
    stop 'NCHUNKS must be either 1, 2, 3 or 6'

  call read_value_double_precision(ANGULAR_WIDTH_XI_IN_DEGREES, 'mesher.ANGULAR_WIDTH_XI_IN_DEGREES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_double_precision(ANGULAR_WIDTH_ETA_IN_DEGREES, 'mesher.ANGULAR_WIDTH_ETA_IN_DEGREES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_double_precision(CENTER_LATITUDE_IN_DEGREES, 'mesher.CENTER_LATITUDE_IN_DEGREES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_double_precision(CENTER_LONGITUDE_IN_DEGREES, 'mesher.CENTER_LONGITUDE_IN_DEGREES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_double_precision(GAMMA_ROTATION_AZIMUTH, 'mesher.GAMMA_ROTATION_AZIMUTH')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

! this MUST be 90 degrees for two chunks or more to match geometrically
  if(NCHUNKS > 1 .and. abs(ANGULAR_WIDTH_XI_IN_DEGREES - 90.d0) > 0.00000001d0) &
    stop 'ANGULAR_WIDTH_XI_IN_DEGREES must be 90 for more than one chunk'

! this can be any value in the case of two chunks
  if(NCHUNKS > 2 .and. abs(ANGULAR_WIDTH_ETA_IN_DEGREES - 90.d0) > 0.00000001d0) &
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
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(NEX_ETA, 'mesher.NEX_ETA')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(NPROC_XI, 'mesher.NPROC_XI')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(NPROC_ETA, 'mesher.NPROC_ETA')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'


! define the velocity model
  call read_value_string(MODEL, 'model.name')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

! use PREM as the 1D reference model by default
  REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM
  HONOR_1D_SPHERICAL_MOHO = .false.
  ONE_CRUST = .false.
  CASE_3D = .false.
! default is no 3D model
  THREE_D_MODEL = 0

  if(MODEL == '1D_isotropic_prem') then
    TRANSVERSE_ISOTROPY = .false.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.
    HONOR_1D_SPHERICAL_MOHO = .true.

  else if(MODEL == '1D_transversely_isotropic_prem') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.
    HONOR_1D_SPHERICAL_MOHO = .true.

  else if(MODEL == '1D_iasp91' .or. MODEL == '1D_1066a' .or. MODEL == '1D_ak135') then
    if(MODEL == '1D_iasp91') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_IASP91
    else if(MODEL == '1D_1066a') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_1066A
    else if(MODEL == '1D_ak135') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_AK135
    else
      stop 'reference 1D Earth model unknown'
    endif
    TRANSVERSE_ISOTROPY = .false.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.
    HONOR_1D_SPHERICAL_MOHO = .true.

  else if(MODEL == '1D_ref') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.
    HONOR_1D_SPHERICAL_MOHO = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_REF

  else if(MODEL == '1D_isotropic_prem_onecrust') then
    TRANSVERSE_ISOTROPY = .false.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.
    HONOR_1D_SPHERICAL_MOHO = .true.
    ONE_CRUST = .true.

  else if(MODEL == '1D_transversely_isotropic_prem_onecrust') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.
    HONOR_1D_SPHERICAL_MOHO = .true.
    ONE_CRUST = .true.

  else if(MODEL == '1D_iasp91_onecrust' .or. MODEL == '1D_1066a_onecrust' .or. MODEL == '1D_ak135_onecrust') then
    if(MODEL == '1D_iasp91_onecrust') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_IASP91
    else if(MODEL == '1D_1066a_onecrust') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_1066A
    else if(MODEL == '1D_ak135_onecrust') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_AK135
    else
      stop 'reference 1D Earth model unknown'
    endif
    TRANSVERSE_ISOTROPY = .false.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.
    HONOR_1D_SPHERICAL_MOHO = .true.
    ONE_CRUST = .true.

  else if(MODEL == 'transversely_isotropic_prem_plus_3D_crust_2.0') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .true.
    ATTENUATION_3D = .false.
    ONE_CRUST = .true.
    CASE_3D = .true.

  else if(MODEL == 's20rts') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .true.
    ATTENUATION_3D = .false.
    ONE_CRUST = .true.
    CASE_3D = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM
    THREE_D_MODEL = THREE_D_MODEL_S20RTS

  else if(MODEL == 's362ani') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .true.
    ATTENUATION_3D = .false.
    ONE_CRUST = .true.
    CASE_3D = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_REF
    THREE_D_MODEL = THREE_D_MODEL_S362ANI

  else if(MODEL == 's362wmani') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .true.
    ATTENUATION_3D = .false.
    ONE_CRUST = .true.
    CASE_3D = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_REF
    THREE_D_MODEL = THREE_D_MODEL_S362WMANI

  else if(MODEL == 's362ani_prem') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .true.
    ATTENUATION_3D = .false.
    ONE_CRUST = .true.
    CASE_3D = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM
    THREE_D_MODEL = THREE_D_MODEL_S362ANI_PREM

  else if(MODEL == 's29ea') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .true.
    ATTENUATION_3D = .false.
    ONE_CRUST = .true.
    CASE_3D = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_REF
    THREE_D_MODEL = THREE_D_MODEL_S29EA

  else if(MODEL == '3D_attenuation') then
    TRANSVERSE_ISOTROPY = .false.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .true.
    ONE_CRUST = .true.
    CASE_3D = .true.

  else if(MODEL == '3D_anisotropic') then
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .false.
    ANISOTROPIC_3D_MANTLE = .true.
    ANISOTROPIC_INNER_CORE = .false.
    CRUSTAL = .false.
    ATTENUATION_3D = .false.
    ONE_CRUST = .true.
    CASE_3D = .true.

  else
    stop 'model not implemented, edit read_compute_parameters.f90 and recompile'
  endif

! set time step, radial distribution of elements, and attenuation period range
! right distribution is determined based upon maximum value of NEX
  NEX_MAX = max(NEX_XI,NEX_ETA)

  if(ANGULAR_WIDTH_XI_IN_DEGREES <= 89.9999d0 .OR. NEX_MAX > 1248) then

!!!!!! DK DK
!!!!!! DK DK  this section written by Brian Savage, commented out by Dimitri Komatitsch
!!!!!! DK DK  because it is based on the old mesher and therefore does not work with the new
!!!!!! DK DK  mesher. Brian will update it and put it back.
!!!!!! DK DK

    stop 'auto_ner commented out by Dimitri Komatitsch for now; Brian Savage will put it back'

!   call auto_ner(ANGULAR_WIDTH_XI_IN_DEGREES, NEX_MAX, &
!         NER_CRUST, NER_220_MOHO, NER_400_220, NER_600_400, &
!         NER_670_600, NER_771_670, NER_TOPDDOUBLEPRIME_771, &
!         NER_CMB_TOPDDOUBLEPRIME, NER_TOP_CENTRAL_CUBE_ICB)

    ! Determine in appropriate period range for the current mesh
    !   Note: we are using DT as a temporary variable
    ! The Minimum attenuation period = (Grid Spacing in km) / V_min
    !  Grid spacing in km     = Width of an element in km * spacing for GLL point * points per wavelength
    !  Width of element in km = (Angular width in degrees / NEX_MAX) * degrees to km
    !  degrees to km          = 111.11d0
    !  spacing for GLL point  = 4
    !  points per wavelength  = 4
!   DT = (ANGULAR_WIDTH_XI_IN_DEGREES / (4.0d0 * dble(NEX_MAX)) * 111.11d0 * 4.0d0 ) / 2.25d0
!   MIN_ATTENUATION_PERIOD = DT

    ! The max attenuation period for 3 SLS is optimally
    !   1.75 decades from the min attenuation period
!   DT = dble(MIN_ATTENUATION_PERIOD) * 10.0d0**1.75d0
!   MAX_ATTENUATION_PERIOD = DT

    ! 0.173 is the minimum spacing between GLL points for NGLL = 5
    ! This should be changed in the future, placed in a header file
    ! 0.40 is the ratio of radial lengths of elements inside the
    ! central cube to those just outside the central cube
    ! 1221.0 is the Radius of the inncer core in km
    ! 0.40 is the maximum stability condition
    ! 11.02827 is Vp near the inner core boundary
    ! See equation 48 in Komatitsch and Tromp (2002, Part I)
!   DT = (0.40d0 * &
!        ((ANGULAR_WIDTH_XI_IN_DEGREES * (PI/180.0d0)) * 1221.0d0) / &
!        (dble(NEX_MAX) / 8.0d0) / 11.02827d0 ) * 0.173d0 * 0.4d0

  else

!----
!----  case prem_onecrust by default
!----

     ! element width =   0.5625000      degrees =    62.54715      km
      if(NEX_MAX <= 160) then
        DT                       = 0.252d0

        MIN_ATTENUATION_PERIOD   = 20
        MAX_ATTENUATION_PERIOD   = 1000

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
        R_CENTRAL_CUBE = 950000.d0

    ! element width =   0.3515625      degrees =    39.09196      km
      else if(NEX_MAX <= 256) then
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
      else if(NEX_MAX <= 320) then
        DT                       = 0.16d0

        MIN_ATTENUATION_PERIOD   = 20
        MAX_ATTENUATION_PERIOD   = 1000

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
      else if(NEX_MAX <= 480) then
        DT                       = 0.12d0

        MIN_ATTENUATION_PERIOD   = 20
        MAX_ATTENUATION_PERIOD   = 1000

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
      else if(NEX_MAX <= 512) then
        DT                       = 0.1125d0

        MIN_ATTENUATION_PERIOD   = 8
        MAX_ATTENUATION_PERIOD   = 1000

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
      else if(NEX_MAX <= 640) then
        DT                       = 0.09d0

        MIN_ATTENUATION_PERIOD   = 5
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
      else if(NEX_MAX <= 864) then
        DT                       = 0.0667d0

        MIN_ATTENUATION_PERIOD   = 5
        MAX_ATTENUATION_PERIOD   = 400

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
      else if(NEX_MAX <= 1152) then
        DT                       = 0.05d0

        MIN_ATTENUATION_PERIOD   = 4
        MAX_ATTENUATION_PERIOD   = 400

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
      else if(NEX_MAX <= 1248) then
        DT                       = 0.0462d0

        MIN_ATTENUATION_PERIOD   = 4
        MAX_ATTENUATION_PERIOD   = 400

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
        stop 'problem with this value of NEX_MAX'
      endif

!----
!----  change some values in the case of regular PREM with two crustal layers or of 3D models
!----
    if (HONOR_1D_SPHERICAL_MOHO .and. .not. ONE_CRUST) then

! case of regular PREM with two crustal layers: change the time step for small meshes
! because of a different size of elements in the radial direction in the crust
      if(NEX_MAX <= 160) then
        DT = 0.20d0
      else if(NEX_MAX <= 256) then
        DT = 0.20d0
      endif

    else

! case 3D: change the upper part of the mesh and the time step because of
! the 3D model with faster of slower velocities in the upper mantle and crust
      if(NEX_MAX <= 160) then
        DT                       = 0.15d0
        NER_CRUST                = 2
      else if(NEX_MAX <= 256) then
        DT                       = 0.17d0
        NER_CRUST                = 2
      else if(NEX_MAX <= 320) then
        DT                       = 0.155d0
        NER_CRUST                = 2
      else if(NEX_MAX <= 480) then
        NER_CRUST                = 2
        NER_80_MOHO              = 1
      else if(NEX_MAX <= 512) then
        NER_CRUST                = 2
        NER_80_MOHO              = 1
      else if(NEX_MAX <= 640) then
        NER_CRUST                = 3
        NER_80_MOHO              = 2
      else if(NEX_MAX <= 864) then
        NER_CRUST                = 4
        NER_80_MOHO              = 3
      else if(NEX_MAX <= 1152) then
        NER_CRUST                = 4
        NER_80_MOHO              = 4
      else if(NEX_MAX <= 1248) then
        NER_CRUST                = 5
        NER_80_MOHO              = 4
      endif

    endif

  endif

! take a 5% safety margin on the maximum stable time step
! which was obtained by trial and error
  DT = DT * (1.d0 - 0.05d0)

  call read_value_logical(OCEANS, 'model.OCEANS')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(ELLIPTICITY, 'model.ELLIPTICITY')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(TOPOGRAPHY, 'model.TOPOGRAPHY')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(GRAVITY, 'model.GRAVITY')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(ROTATION, 'model.ROTATION')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(ATTENUATION, 'model.ATTENUATION')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

  call read_value_logical(ABSORBING_CONDITIONS, 'solver.ABSORBING_CONDITIONS')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

  if(ABSORBING_CONDITIONS .and. NCHUNKS == 6) stop 'cannot have absorbing conditions in the full Earth'

  if(ABSORBING_CONDITIONS .and. NCHUNKS == 3) stop 'absorbing conditions not supported for three chunks yet'

  if(ATTENUATION_3D .and. .not. ATTENUATION) stop 'need ATTENUATION to use ATTENUATION_3D'

! radii in PREM or IASP91
! and normalized density at fluid-solid interface on fluid size for coupling
! ROCEAN: radius of the ocean (m)
! RMIDDLE_CRUST: radius of the middle crust (m)
! RMOHO: radius of the Moho (m)
! R80: radius of 80 km discontinuity (m)
! R120: radius of 120 km discontinuity (m) in IASP91
! R220: radius of 220 km discontinuity (m)
! R400: radius of 400 km discontinuity (m)
! R600: radius of 600 km 2nd order discontinuity (m)
! R670: radius of 670 km discontinuity (m)
! R771: radius of 771 km 2nd order discontinuity (m)
! RTOPDDOUBLEPRIME: radius of top of D" 2nd order discontinuity (m)
! RCMB: radius of CMB (m)
! RICB: radius of ICB (m)

! by default there is no d120 discontinuity, except in IASP91, therefore set to fictitious value
  R120 = -1.d0

! value common to all models
  RHO_OCEANS = 1020.0 / RHOAV

  if(REFERENCE_1D_MODEL == REFERENCE_MODEL_IASP91) then

! IASP91
    ROCEAN = 6371000.d0
    RMIDDLE_CRUST = 6351000.d0
    RMOHO = 6336000.d0
    R80  = 6291000.d0
    R120 = 6251000.d0
    R220 = 6161000.d0
    R400 = 5961000.d0
! there is no d600 discontinuity in IASP91 therefore this value is useless
! but it needs to be there for compatibility with other subroutines
    R600 = R_EARTH - 600000.d0
    R670 = 5711000.d0
    R771 = 5611000.d0
    RTOPDDOUBLEPRIME = 3631000.d0
    RCMB = 3482000.d0
    RICB = 1217000.d0

    RHO_TOP_OC = 9900.2379 / RHOAV
    RHO_BOTTOM_OC = 12168.6383 / RHOAV

  else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135) then

!! DK DK UGLY our implementation of AK135 has not been checked carefully yet
!! DK DK UGLY therefore let us doublecheck it carefully one day before using it

! values below corrected by Ying Zhou <yingz@gps.caltech.edu>

! AK135 without the 300 meters of mud layer
   ROCEAN = 6368000.d0
   RMIDDLE_CRUST = 6361000.d0
   RMOHO  = 6353000.d0
   R80    = 6291000.d0
   R220   = 6161000.d0
   R400   = 5961000.d0
   R670   = 5711000.d0
   RTOPDDOUBLEPRIME = 3631000.d0
   RCMB   = 3479500.d0
   RICB   = 1217500.d0

! values for AK135 that are not discontinuities
   R600 = 5771000.d0
   R771 = 5611000.d0

   RHO_TOP_OC = 9914.5000 / RHOAV
   RHO_BOTTOM_OC = 12139.1000 / RHOAV

  else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_1066A) then

! values below corrected by Ying Zhou <yingz@gps.caltech.edu>

! 1066A
   RMOHO = 6360000.d0
   R400 = 5950000.d0
   R600 = 5781000.d0
   R670 = 5700000.d0
   RCMB = 3484300.d0
   RICB = 1229480.d0

! values for 1066A that are not discontinuities
   RTOPDDOUBLEPRIME = 3631000.d0
   R220 = 6161000.d0
   R771 = 5611000.d0
! RMIDDLE_CRUST used only for high resolution FFSW1C model, with 3 elements crust simulations
! mid_crust = 10 km
   RMIDDLE_CRUST = 6361000.d0
   R80 = 6291000.d0

! model 1066A has no oceans, therefore we use the radius of the Earth instead
   ROCEAN = R_EARTH

   RHO_TOP_OC = 9917.4500 / RHOAV
   RHO_BOTTOM_OC = 12160.6500 / RHOAV

  else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_REF) then

! REF
    ROCEAN = 6368000.d0
    RMIDDLE_CRUST = 6356000.d0
    RMOHO = 6346600.d0
    R80  = 6291000.d0
    R220 = 6151000.d0
    R400 = 5961000.d0
    R600 = 5771000.d0
    R670 = 5721000.d0
    R771 = 5600000.d0
    RTOPDDOUBLEPRIME = 3630000.d0
    RCMB = 3479958.d0
    RICB = 1221491.d0

    RHO_TOP_OC = 9903.48 / RHOAV
    RHO_BOTTOM_OC = 12166.35 / RHOAV

  else

! PREM
    ROCEAN = 6368000.d0
    RMIDDLE_CRUST = 6356000.d0
    RMOHO = 6346600.d0
    R80  = 6291000.d0
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

! honor the PREM Moho or define a fictitious Moho in order to have even radial sampling
! from the d220 to the Earth surface
  if(HONOR_1D_SPHERICAL_MOHO) then
    RMOHO_FICTITIOUS_IN_MESHER = RMOHO
  else
    RMOHO_FICTITIOUS_IN_MESHER = (R80 + R_EARTH) / 2
  endif

! non-dimensionalized size of central cube in the inner core
! This is where the central cube in the inner core and the rest of the mesh
! are matched (150 km below the ICB is optimal)
!  R_CENTRAL_CUBE = RICB - 150000.d0
!  R_CENTRAL_CUBE = 965000.d0

  call read_value_double_precision(RECORD_LENGTH_IN_MINUTES, 'solver.RECORD_LENGTH_IN_MINUTES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

! compute total number of time steps, rounded to next multiple of 100
  NSTEP = 100 * (int(RECORD_LENGTH_IN_MINUTES * 60.d0 / (100.d0*DT)) + 1)

  call read_value_logical(MOVIE_SURFACE, 'solver.MOVIE_SURFACE')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(MOVIE_VOLUME, 'solver.MOVIE_VOLUME')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(NTSTEP_BETWEEN_FRAMES, 'solver.NTSTEP_BETWEEN_FRAMES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_double_precision(HDUR_MOVIE, 'solver.HDUR_MOVIE')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

! computes a default hdur_movie that creates nice looking movies.
! Sets HDUR_MOVIE as the minimum period the mesh can resolve
  if(HDUR_MOVIE <= TINYVAL) &
    HDUR_MOVIE = 1.1d0*max(240.d0/NEX_XI*18.d0*ANGULAR_WIDTH_XI_IN_DEGREES/90.d0, &
                           240.d0/NEX_ETA*18.d0*ANGULAR_WIDTH_ETA_IN_DEGREES/90.d0)


  call read_value_logical(SAVE_MESH_FILES, 'mesher.SAVE_MESH_FILES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(NUMBER_OF_RUNS, 'solver.NUMBER_OF_RUNS')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(NUMBER_OF_THIS_RUN, 'solver.NUMBER_OF_THIS_RUN')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_string(LOCAL_PATH, 'LOCAL_PATH')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(NTSTEP_BETWEEN_OUTPUT_INFO, 'solver.NTSTEP_BETWEEN_OUTPUT_INFO')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(NTSTEP_BETWEEN_OUTPUT_SEISMOS, 'solver.NTSTEP_BETWEEN_OUTPUT_SEISMOS')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(NTSTEP_BETWEEN_READ_ADJSRC, 'solver.NTSTEP_BETWEEN_READ_ADJSRC')
  if(err_occurred() /= 0) return

  call read_value_logical(OUTPUT_SEISMOS_ASCII_TEXT, 'solver.OUTPUT_SEISMOS_ASCII_TEXT')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(OUTPUT_SEISMOS_SAC_ALPHANUM, 'solver.OUTPUT_SEISMOS_SAC_ALPHANUM')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(OUTPUT_SEISMOS_SAC_BINARY, 'solver.OUTPUT_SEISMOS_SAC_BINARY')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(ROTATE_SEISMOGRAMS_RT, 'solver.ROTATE_SEISMOGRAMS_RT')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

  call read_value_logical(RECEIVERS_CAN_BE_BURIED, 'solver.RECEIVERS_CAN_BE_BURIED')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_logical(PRINT_SOURCE_TIME_FUNCTION, 'solver.PRINT_SOURCE_TIME_FUNCTION')

  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

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

! check that sphere can be cut into slices without getting negative Jacobian
  if(NEX_XI < 48) stop 'NEX_XI must be greater than 48 to cut the sphere into slices with positive Jacobian'
  if(NEX_ETA < 48) stop 'NEX_ETA must be greater than 48 to cut the sphere into slices with positive Jacobian'

! check that mesh can be coarsened in depth four times (block size must be a multiple of 32)
  if(mod(NEX_XI,32) /= 0) stop 'NEX_XI must be a multiple of 32'
  if(mod(NEX_ETA,32) /= 0) stop 'NEX_ETA must be a multiple of 32'
  if(mod(NEX_XI/32,NPROC_XI) /= 0) stop 'NEX_XI must be a multiple of 32*NPROC_XI'
  if(mod(NEX_ETA/32,NPROC_ETA) /= 0) stop 'NEX_ETA must be a multiple of 32*NPROC_ETA'

! check that topology is correct if more than two chunks
  if(NCHUNKS > 2 .and. NEX_XI /= NEX_ETA) stop 'must have NEX_XI = NEX_ETA for more than two chunks'
  if(NCHUNKS > 2 .and. NPROC_XI /= NPROC_ETA) stop 'must have NPROC_XI = NPROC_ETA for more than two chunks'

! check that IASP91, AK135, or 1066A is isotropic
  if((REFERENCE_1D_MODEL == REFERENCE_MODEL_IASP91 .or. &
      REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135 .or. &
      REFERENCE_1D_MODEL == REFERENCE_MODEL_1066A) .and. TRANSVERSE_ISOTROPY) &
        stop 'models IASP91, AK135 and 1066A are currently isotropic'

  ELEMENT_WIDTH = ANGULAR_WIDTH_XI_IN_DEGREES/dble(NEX_MAX) * DEGREES_TO_RADIANS

!
!--- compute additional parameters
!

! number of elements horizontally in each slice (i.e. per processor)
! these two values MUST be equal in all cases
  NEX_PER_PROC_XI = NEX_XI / NPROC_XI
  NEX_PER_PROC_ETA = NEX_ETA / NPROC_ETA

! total number of processors in each of the six chunks
  NPROC = NPROC_XI * NPROC_ETA

! total number of processors in the full Earth composed of the six chunks
  NPROCTOT = NCHUNKS * NPROC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  definition of general mesh parameters below
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! find element below top of which we should implement the second doubling in the mantle
! locate element closest to optimal value
  distance_min = HUGEVAL
  do ielem = 2,NER_TOPDDOUBLEPRIME_771
    zval = RTOPDDOUBLEPRIME + ielem * (R771 - RTOPDDOUBLEPRIME) / dble(NER_TOPDDOUBLEPRIME_771)
    distance = abs(zval - (R_EARTH - DEPTH_SECOND_DOUBLING_OPTIMAL))
    if(distance < distance_min) then
      elem_doubling_mantle = ielem
      distance_min = distance
      DEPTH_SECOND_DOUBLING_REAL = R_EARTH - zval
    endif
  enddo

! find element below top of which we should implement the third doubling in the middle of the outer core
! locate element closest to optimal value
  distance_min = HUGEVAL
! start at element number 4 because we need at least two elements below for the fourth doubling
! implemented at the bottom of the outer core
  do ielem = 4,NER_OUTER_CORE
    zval = RICB + ielem * (RCMB - RICB) / dble(NER_OUTER_CORE)
    distance = abs(zval - (R_EARTH - DEPTH_THIRD_DOUBLING_OPTIMAL))
    if(distance < distance_min) then
      elem_doubling_middle_outer_core = ielem
      distance_min = distance
      DEPTH_THIRD_DOUBLING_REAL = R_EARTH - zval
    endif
  enddo

! find element below top of which we should implement the fourth doubling in the middle of the outer core
! locate element closest to optimal value
  distance_min = HUGEVAL
! end two elements before the top because we need at least two elements above for the third doubling
! implemented in the middle of the outer core
  do ielem = 2,NER_OUTER_CORE-2
    zval = RICB + ielem * (RCMB - RICB) / dble(NER_OUTER_CORE)
    distance = abs(zval - (R_EARTH - DEPTH_FOURTH_DOUBLING_OPTIMAL))
    if(distance < distance_min) then
      elem_doubling_bottom_outer_core = ielem
      distance_min = distance
      DEPTH_FOURTH_DOUBLING_REAL = R_EARTH - zval
    endif
  enddo

! make sure that the two doublings in the outer core are found in the right order
  if(elem_doubling_bottom_outer_core >= elem_doubling_middle_outer_core) &
                  stop 'error in location of the two doublings in the outer core'

! define all the layers of the mesh
  if (SUPPRESS_CRUSTAL_MESH) then

    NER_80_MOHO = nint(NER_220_80*((R_EARTH-R80)*1.d0)/((R80-R220)*1.d0))
    NER_CRUST = 0
    RMOHO_FICTITIOUS_IN_MESHER = R_EARTH

    OCEANS= .false.
    TOPOGRAPHY = .false.
    CRUSTAL = .false.

    NUMBER_OF_MESH_LAYERS = 14
    layer_offset = 0

    ner( 1) = NER_CRUST
    ner( 2) = NER_80_MOHO
    ner( 3) = NER_220_80
    ner( 4) = NER_400_220
    ner( 5) = NER_600_400
    ner( 6) = NER_670_600
    ner( 7) = NER_771_670
    ner( 8) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
    ner( 9) = elem_doubling_mantle
    ner(10) = NER_CMB_TOPDDOUBLEPRIME
    ner(11) = NER_OUTER_CORE - elem_doubling_middle_outer_core
    ner(12) = elem_doubling_middle_outer_core - elem_doubling_bottom_outer_core
    ner(13) = elem_doubling_bottom_outer_core
    ner(14) = NER_TOP_CENTRAL_CUBE_ICB

  ! value of the doubling ratio in each radial region of the mesh
    ratio_sampling_array(1:8) = 2
    ratio_sampling_array(9:11) = 4
    ratio_sampling_array(12) = 8
    ratio_sampling_array(13:14) = 16

  ! value of the doubling index flag in each radial region of the mesh
    doubling_index(1:2) = IFLAG_80_MOHO
    doubling_index(3) = IFLAG_220_80
    doubling_index(4:6) = IFLAG_670_220
    doubling_index(7:10) = IFLAG_MANTLE_NORMAL
    doubling_index(11:13) = IFLAG_OUTER_CORE_NORMAL
    doubling_index(14) = IFLAG_INNER_CORE_NORMAL

  ! define the three regions in which we implement a mesh doubling at the top of that region
    this_region_has_a_doubling(:)  = .false.
    this_region_has_a_doubling(9)  = .true.
    this_region_has_a_doubling(12) = .true.
    this_region_has_a_doubling(13) = .true.

    r_top(1) = R_EARTH
    r_bottom(1) = RMOHO_FICTITIOUS_IN_MESHER

    r_top(2) = RMOHO_FICTITIOUS_IN_MESHER
    r_bottom(2) = R80

    r_top(3) = R80
    r_bottom(3) = R220

    r_top(4) = R220
    r_bottom(4) = R400

    r_top(5) = R400
    r_bottom(5) = R600

    r_top(6) = R600
    r_bottom(6) = R670

    r_top(7) = R670
    r_bottom(7) = R771

    r_top(8) = R771
    r_bottom(8) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL

    r_top(9) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL
    r_bottom(9) = RTOPDDOUBLEPRIME

    r_top(10) = RTOPDDOUBLEPRIME
    r_bottom(10) = RCMB

    r_top(11) = RCMB
    r_bottom(11) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL

    r_top(12) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL
    r_bottom(12) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL

    r_top(13) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL
    r_bottom(13) = RICB

    r_top(14) = RICB
    r_bottom(14) = R_CENTRAL_CUBE

  !!! DM new definition of rmins & rmaxs in replacement of mesh_radial
    rmaxs(1) = ONE
    rmins(1) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH

    rmaxs(2) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH
    rmins(2) = R80 / R_EARTH

    rmaxs(3) = R80 / R_EARTH
    rmins(3) = R220 / R_EARTH

    rmaxs(4) = R220 / R_EARTH
    rmins(4) = R400 / R_EARTH

    rmaxs(5) = R400 / R_EARTH
    rmins(5) = R600 / R_EARTH

    rmaxs(6) = R600 / R_EARTH
    rmins(6) = R670 / R_EARTH

    rmaxs(7) = R670 / R_EARTH
    rmins(7) = R771 / R_EARTH

    rmaxs(8:9) = R771 / R_EARTH
    rmins(8:9) = RTOPDDOUBLEPRIME / R_EARTH

    rmaxs(10) = RTOPDDOUBLEPRIME / R_EARTH
    rmins(10) = RCMB / R_EARTH

    rmaxs(11:13) = RCMB / R_EARTH
    rmins(11:13) = RICB / R_EARTH

    rmaxs(14) = RICB / R_EARTH
    rmins(14) = R_CENTRAL_CUBE / R_EARTH

  elseif (ONE_CRUST) then


    NUMBER_OF_MESH_LAYERS = 14
    layer_offset = 0

    ner( 1) = NER_CRUST
    ner( 2) = NER_80_MOHO
    ner( 3) = NER_220_80
    ner( 4) = NER_400_220
    ner( 5) = NER_600_400
    ner( 6) = NER_670_600
    ner( 7) = NER_771_670
    ner( 8) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
    ner( 9) = elem_doubling_mantle
    ner(10) = NER_CMB_TOPDDOUBLEPRIME
    ner(11) = NER_OUTER_CORE - elem_doubling_middle_outer_core
    ner(12) = elem_doubling_middle_outer_core - elem_doubling_bottom_outer_core
    ner(13) = elem_doubling_bottom_outer_core
    ner(14) = NER_TOP_CENTRAL_CUBE_ICB

  ! value of the doubling ratio in each radial region of the mesh
    ratio_sampling_array(1) = 1
    ratio_sampling_array(2:8) = 2
    ratio_sampling_array(9:11) = 4
    ratio_sampling_array(12) = 8
    ratio_sampling_array(13:14) = 16

  ! value of the doubling index flag in each radial region of the mesh
    doubling_index(1) = IFLAG_CRUST
    doubling_index(2) = IFLAG_80_MOHO
    doubling_index(3) = IFLAG_220_80
    doubling_index(4:6) = IFLAG_670_220
    doubling_index(7:10) = IFLAG_MANTLE_NORMAL
    doubling_index(11:13) = IFLAG_OUTER_CORE_NORMAL
    doubling_index(14) = IFLAG_INNER_CORE_NORMAL

  ! define the three regions in which we implement a mesh doubling at the top of that region
    this_region_has_a_doubling(:)  = .false.
    this_region_has_a_doubling(2)  = .true.
    this_region_has_a_doubling(9)  = .true.
    this_region_has_a_doubling(12) = .true.
    this_region_has_a_doubling(13) = .true.

  ! define the top and bottom radii of all the regions of the mesh in the radial direction
  ! the first region is the crust at the surface of the Earth
  ! the last region is in the inner core near the center of the Earth

!!!!!!!!!!! DK DK UGLY: beware, is there a bug when 3D crust crosses anisotropy in the mantle?
!!!!!!!!!!! DK DK UGLY: i.e. if there is no thick crust there, some elements above the Moho
!!!!!!!!!!! DK DK UGLY: should be anisotropic but anisotropy is currently only
!!!!!!!!!!! DK DK UGLY: stored between d220 and MOHO to save memory? Clarify this one day.
!!!!!!!!!!! DK DK UGLY: The Moho stretching and squishing that Jeroen added to V4.0
!!!!!!!!!!! DK DK UGLY: should partly deal with this problem.

    r_top(1) = R_EARTH
    r_bottom(1) = RMOHO_FICTITIOUS_IN_MESHER

    r_top(2) = RMOHO_FICTITIOUS_IN_MESHER
    r_bottom(2) = R80

    r_top(3) = R80
    r_bottom(3) = R220

    r_top(4) = R220
    r_bottom(4) = R400

    r_top(5) = R400
    r_bottom(5) = R600

    r_top(6) = R600
    r_bottom(6) = R670

    r_top(7) = R670
    r_bottom(7) = R771

    r_top(8) = R771
    r_bottom(8) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL

    r_top(9) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL
    r_bottom(9) = RTOPDDOUBLEPRIME

    r_top(10) = RTOPDDOUBLEPRIME
    r_bottom(10) = RCMB

    r_top(11) = RCMB
    r_bottom(11) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL

    r_top(12) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL
    r_bottom(12) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL

    r_top(13) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL
    r_bottom(13) = RICB

    r_top(14) = RICB
    r_bottom(14) = R_CENTRAL_CUBE

  ! new definition of rmins & rmaxs
    rmaxs(1) = ONE
    rmins(1) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH

    rmaxs(2) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH
    rmins(2) = R80 / R_EARTH

    rmaxs(3) = R80 / R_EARTH
    rmins(3) = R220 / R_EARTH

    rmaxs(4) = R220 / R_EARTH
    rmins(4) = R400 / R_EARTH

    rmaxs(5) = R400 / R_EARTH
    rmins(5) = R600 / R_EARTH

    rmaxs(6) = R600 / R_EARTH
    rmins(6) = R670 / R_EARTH

    rmaxs(7) = R670 / R_EARTH
    rmins(7) = R771 / R_EARTH

    rmaxs(8:9) = R771 / R_EARTH
    rmins(8:9) = RTOPDDOUBLEPRIME / R_EARTH

    rmaxs(10) = RTOPDDOUBLEPRIME / R_EARTH
    rmins(10) = RCMB / R_EARTH

    rmaxs(11:13) = RCMB / R_EARTH
    rmins(11:13) = RICB / R_EARTH

    rmaxs(14) = RICB / R_EARTH
    rmins(14) = R_CENTRAL_CUBE / R_EARTH
  else

    NUMBER_OF_MESH_LAYERS = 15
    layer_offset = 1
! DM a revoir
    if (NER_CRUST<2) NER_CRUST=2
    if ((RMIDDLE_CRUST-RMOHO_FICTITIOUS_IN_MESHER)<(R_EARTH-RMIDDLE_CRUST)) then
      ner( 1) = ceiling (NER_CRUST / 2.d0)
      ner( 2) = floor (NER_CRUST / 2.d0)
    else
      ner( 1) = floor (NER_CRUST / 2.d0)
      ner( 2) = ceiling (NER_CRUST / 2.d0)
    endif
    ner( 3) = NER_80_MOHO
    ner( 4) = NER_220_80
    ner( 5) = NER_400_220
    ner( 6) = NER_600_400
    ner( 7) = NER_670_600
    ner( 8) = NER_771_670
    ner( 9) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
    ner(10) = elem_doubling_mantle
    ner(11) = NER_CMB_TOPDDOUBLEPRIME
    ner(12) = NER_OUTER_CORE - elem_doubling_middle_outer_core
    ner(13) = elem_doubling_middle_outer_core - elem_doubling_bottom_outer_core
    ner(14) = elem_doubling_bottom_outer_core
    ner(15) = NER_TOP_CENTRAL_CUBE_ICB

  ! value of the doubling ratio in each radial region of the mesh
    ratio_sampling_array(1:2) = 1
    ratio_sampling_array(3:9) = 2
    ratio_sampling_array(10:12) = 4
    ratio_sampling_array(13) = 8
    ratio_sampling_array(14:15) = 16

  ! value of the doubling index flag in each radial region of the mesh
    doubling_index(1:2) = IFLAG_CRUST
    doubling_index(3) = IFLAG_80_MOHO
    doubling_index(4) = IFLAG_220_80
    doubling_index(5:7) = IFLAG_670_220
    doubling_index(8:11) = IFLAG_MANTLE_NORMAL
    doubling_index(12:14) = IFLAG_OUTER_CORE_NORMAL
    doubling_index(15) = IFLAG_INNER_CORE_NORMAL

  ! define the three regions in which we implement a mesh doubling at the top of that region
    this_region_has_a_doubling(:)  = .false.
    this_region_has_a_doubling(3)  = .true.
    this_region_has_a_doubling(10) = .true.
    this_region_has_a_doubling(13) = .true.
    this_region_has_a_doubling(14) = .true.

  ! define the top and bottom radii of all the regions of the mesh in the radial direction
  ! the first region is the crust at the surface of the Earth
  ! the last region is in the inner core near the center of the Earth

    r_top(1) = R_EARTH
    r_bottom(1) = RMIDDLE_CRUST

    r_top(2) = RMIDDLE_CRUST
    r_bottom(2) = RMOHO_FICTITIOUS_IN_MESHER

    r_top(3) = RMOHO_FICTITIOUS_IN_MESHER
    r_bottom(3) = R80

    r_top(4) = R80
    r_bottom(4) = R220

    r_top(5) = R220
    r_bottom(5) = R400

    r_top(6) = R400
    r_bottom(6) = R600

    r_top(7) = R600
    r_bottom(7) = R670

    r_top(8) = R670
    r_bottom(8) = R771

    r_top(9) = R771
    r_bottom(9) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL

    r_top(10) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL
    r_bottom(10) = RTOPDDOUBLEPRIME

    r_top(11) = RTOPDDOUBLEPRIME
    r_bottom(11) = RCMB

    r_top(12) = RCMB
    r_bottom(12) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL

    r_top(13) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL
    r_bottom(13) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL

    r_top(14) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL
    r_bottom(14) = RICB

    r_top(15) = RICB
    r_bottom(15) = R_CENTRAL_CUBE

  ! new definition of rmins & rmaxs
    rmaxs(1) = ONE
    rmins(1) = RMIDDLE_CRUST / R_EARTH

    rmaxs(2) = RMIDDLE_CRUST / R_EARTH
    rmins(2) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH

    rmaxs(3) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH
    rmins(3) = R80 / R_EARTH

    rmaxs(4) = R80 / R_EARTH
    rmins(4) = R220 / R_EARTH

    rmaxs(5) = R220 / R_EARTH
    rmins(5) = R400 / R_EARTH

    rmaxs(6) = R400 / R_EARTH
    rmins(6) = R600 / R_EARTH

    rmaxs(7) = R600 / R_EARTH
    rmins(7) = R670 / R_EARTH

    rmaxs(8) = R670 / R_EARTH
    rmins(8) = R771 / R_EARTH

    rmaxs(9:10) = R771 / R_EARTH
    rmins(9:10) = RTOPDDOUBLEPRIME / R_EARTH

    rmaxs(11) = RTOPDDOUBLEPRIME / R_EARTH
    rmins(11) = RCMB / R_EARTH

    rmaxs(12:14) = RCMB / R_EARTH
    rmins(12:14) = RICB / R_EARTH

    rmaxs(15) = RICB / R_EARTH
    rmins(15) = R_CENTRAL_CUBE / R_EARTH
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  calculation of number of elements (NSPEC) below
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  1D case
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! theoretical number of spectral elements in radial direction
do iter_region = IREGION_CRUST_MANTLE,IREGION_INNER_CORE
        if(iter_region == IREGION_CRUST_MANTLE) then
                ifirst_region = 1
                ilast_region = 10 + layer_offset
        else if(iter_region == IREGION_OUTER_CORE) then
                ifirst_region = 11 + layer_offset
                ilast_region = NUMBER_OF_MESH_LAYERS - 1
        else if(iter_region == IREGION_INNER_CORE) then
                ifirst_region = NUMBER_OF_MESH_LAYERS
                ilast_region = NUMBER_OF_MESH_LAYERS
        else
                stop 'incorrect region code detected'
        endif
        NSPEC1D_RADIAL(iter_region) = sum(ner(ifirst_region:ilast_region))
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  2D case
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! exact number of surface elements for faces along XI and ETA

do iter_region = IREGION_CRUST_MANTLE,IREGION_INNER_CORE
    if(iter_region == IREGION_CRUST_MANTLE) then
        ifirst_region = 1
        ilast_region = 10 + layer_offset
    else if(iter_region == IREGION_OUTER_CORE) then
        ifirst_region = 11 + layer_offset
        ilast_region = NUMBER_OF_MESH_LAYERS - 1
    else if(iter_region == IREGION_INNER_CORE) then
        ifirst_region = NUMBER_OF_MESH_LAYERS
        ilast_region = NUMBER_OF_MESH_LAYERS
    else
        stop 'incorrect region code detected'
    endif
    tmp_sum_xi = 0
    tmp_sum_eta = 0
    do iter_layer = ifirst_region, ilast_region
        if (this_region_has_a_doubling(iter_layer)) then
            if (ner(iter_layer) == 1) then
              nb_lay_sb = 1
              nspec2D_xi_sb = NSPEC2D_XI_SUPERBRICK_1L
              nspec2D_eta_sb = NSPEC2D_ETA_SUPERBRICK_1L
            else
              nb_lay_sb = 2
              nspec2D_xi_sb = NSPEC2D_XI_SUPERBRICK
              nspec2D_eta_sb = NSPEC2D_ETA_SUPERBRICK
            endif
            doubling = 1
        else
            doubling = 0
            nb_lay_sb = 0
            nspec2D_xi_sb = 0
            nspec2D_eta_sb = 0
        endif

        tmp_sum_xi = tmp_sum_xi + ((NEX_PER_PROC_XI / ratio_sampling_array(iter_layer)) * &
                (ner(iter_layer) - doubling*nb_lay_sb)) + &
                doubling * ((NEX_PER_PROC_XI / ratio_sampling_array(iter_layer)/2) * nspec2D_xi_sb)

        tmp_sum_eta = tmp_sum_eta + ((NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer)) * &
                (ner(iter_layer) - doubling*nb_lay_sb)) + &
                doubling * ((NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer)/2) * nspec2D_eta_sb)
    enddo
    NSPEC2D_XI(iter_region) = tmp_sum_xi
    NSPEC2D_ETA(iter_region) = tmp_sum_eta
    if (iter_region == IREGION_INNER_CORE .and. INCLUDE_CENTRAL_CUBE) then
        NSPEC2D_XI(iter_region) = NSPEC2D_XI(iter_region) + &
        ((NEX_PER_PROC_XI / 16)*(NEX_XI / 16))
        NSPEC2D_ETA(iter_region) = NSPEC2D_ETA(iter_region) + &
        ((NEX_PER_PROC_ETA / 16)*(NEX_XI / 16))
    endif
enddo

! exact number of surface elements on the bottom and top boundaries

! in the crust and mantle
  NSPEC2D_TOP(IREGION_CRUST_MANTLE) = (NEX_XI/ratio_sampling_array(1))*(NEX_ETA/ratio_sampling_array(1))/NPROC
  NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE) = (NEX_XI/ratio_sampling_array(10+layer_offset))*&
                                         (NEX_ETA/ratio_sampling_array(10+layer_offset))/NPROC

! in the outer core with mesh doubling
  NSPEC2D_TOP(IREGION_OUTER_CORE) = (NEX_XI/4)*(NEX_ETA/4)/NPROC
  NSPEC2D_BOTTOM(IREGION_OUTER_CORE) = (NEX_XI/16)*(NEX_ETA/16)/NPROC

! in the top of the inner core
  NSPEC2D_TOP(IREGION_INNER_CORE) = (NEX_XI/16)*(NEX_ETA/16)/NPROC
  NSPEC2D_BOTTOM(IREGION_INNER_CORE) = NSPEC2D_TOP(IREGION_INNER_CORE)

! maximum number of surface elements on vertical boundaries of the slices
  NSPEC2DMAX_XMIN_XMAX(:) = NSPEC2D_ETA(:)
  NSPEC2DMAX_YMIN_YMAX(:) = NSPEC2D_XI(:)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  3D case
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! exact number of spectral elements in each region

do iter_region = IREGION_CRUST_MANTLE,IREGION_INNER_CORE
    if(iter_region == IREGION_CRUST_MANTLE) then
        ifirst_region = 1
        ilast_region = 10 + layer_offset
    else if(iter_region == IREGION_OUTER_CORE) then
        ifirst_region = 11 + layer_offset
        ilast_region = NUMBER_OF_MESH_LAYERS - 1
    else if(iter_region == IREGION_INNER_CORE) then
        ifirst_region = NUMBER_OF_MESH_LAYERS
        ilast_region = NUMBER_OF_MESH_LAYERS
    else
        stop 'incorrect region code detected'
    endif
    tmp_sum = 0;
    do iter_layer = ifirst_region, ilast_region
        if (this_region_has_a_doubling(iter_layer)) then
            if (ner(iter_layer) == 1) then
              nb_lay_sb = 1
              nspec_sb = NSPEC_SUPERBRICK_1L
            else
              nb_lay_sb = 2
              nspec_sb = NSPEC_DOUBLING_SUPERBRICK
            endif
            doubling = 1
        else
            doubling = 0
            nb_lay_sb = 0
            nspec_sb = 0
        endif
        tmp_sum = tmp_sum + ((NEX_XI / ratio_sampling_array(iter_layer)) * (NEX_ETA / ratio_sampling_array(iter_layer)) * &
                (ner(iter_layer) - doubling*nb_lay_sb)) + &
                doubling * ((NEX_XI / ratio_sampling_array(iter_layer)/2) * (NEX_ETA / ratio_sampling_array(iter_layer)/2) * &
                nspec_sb)
    enddo
    NSPEC(iter_region) = tmp_sum / NPROC
enddo

  if(INCLUDE_CENTRAL_CUBE) NSPEC(IREGION_INNER_CORE) = NSPEC(IREGION_INNER_CORE) + &
         (NEX_PER_PROC_XI / 16) * (NEX_PER_PROC_ETA / 16) * (NEX_XI / 16)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  calculation of number of points (NGLOB) below
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  1D case
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! theoretical number of Gauss-Lobatto points in radial direction
  NGLOB1D_RADIAL(:) = NSPEC1D_RADIAL(:)*(NGLLZ-1)+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  2D case
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 2-D addressing and buffers for summation between slices
! we add one to number of points because of the flag after the last point
  NGLOB2DMAX_XMIN_XMAX(:) = NSPEC2DMAX_XMIN_XMAX(:)*NGLLY*NGLLZ + 1
  NGLOB2DMAX_YMIN_YMAX(:) = NSPEC2DMAX_YMIN_YMAX(:)*NGLLX*NGLLZ + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  3D case
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! exact number of global points in each region

! initialize array
  NGLOB(:) = 0

! in the inner core (no doubling region + eventually central cube)
  if(INCLUDE_CENTRAL_CUBE) then
    NGLOB(IREGION_INNER_CORE) = ((NEX_PER_PROC_XI/16) &
      *(NGLLX-1)+1)*((NEX_PER_PROC_ETA/16) &
      *(NGLLY-1)+1)*((NER_TOP_CENTRAL_CUBE_ICB + NEX_XI / 16)*(NGLLZ-1)+1)
  else
    NGLOB(IREGION_INNER_CORE) = ((NEX_PER_PROC_XI/16) &
      *(NGLLX-1)+1)*((NEX_PER_PROC_ETA/16) &
      *(NGLLY-1)+1)*((NER_TOP_CENTRAL_CUBE_ICB)*(NGLLZ-1)+1)
  endif

! in the crust-mantle and outercore
  do iter_region = IREGION_CRUST_MANTLE,IREGION_OUTER_CORE
      if(iter_region == IREGION_CRUST_MANTLE) then
            ifirst_region = 1
            ilast_region = 10 + layer_offset
      else if(iter_region == IREGION_OUTER_CORE) then
            ifirst_region = 11 + layer_offset
            ilast_region = NUMBER_OF_MESH_LAYERS - 1
      else
            stop 'incorrect region code detected'
      endif
      tmp_sum = 0;
      do iter_layer = ifirst_region, ilast_region
        if (this_region_has_a_doubling(iter_layer)) then
            if (ner(iter_layer) == 1) then
              nb_lay_sb = 1
              nglob_vol = 28*NGLLX**3 - 62*NGLLX**2 + 47*NGLLX - 12
              nglob_surf = 6*NGLLX**2-8*NGLLX+3
              nglob_edge = NGLLX
            else
              nb_lay_sb = 2
              nglob_vol = 32*NGLLX**3 - 70*NGLLX**2 + 52*NGLLX - 13
              nglob_surf = 8*NGLLX**2-11*NGLLX+4
              nglob_edge = 2*NGLLX-1
            endif
            doubling = 1
            padding = -1
        else
            doubling = 0
            padding = 0
            nb_lay_sb = 0
            nglob_vol = 0
            nglob_surf = 0
            nglob_edge = 0
        endif
        if (iter_layer == ilast_region) padding = padding +1
        nblocks_xi = NEX_PER_PROC_XI / ratio_sampling_array(iter_layer)
        nblocks_eta = NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer)

        tmp_sum = tmp_sum + &
        ((nblocks_xi)*(NGLLX-1)+1) * ((nblocks_eta)*(NGLLX-1)+1) * ((ner(iter_layer) - doubling*nb_lay_sb)*(NGLLX-1)+padding)+&
        doubling * (((nblocks_xi*nblocks_eta/4)*nglob_vol) - &
        (((nblocks_eta/2-1)*nblocks_xi/2+(nblocks_xi/2-1)*nblocks_eta/2)*nglob_surf) + &
        ((nblocks_eta/2-1)*(nblocks_xi/2-1)*nglob_edge))
      enddo
      NGLOB(iter_region) = tmp_sum
  enddo

!!! example :
!!!                        nblocks_xi/2=5
!!!                  ____________________________________
!!!                  I      I      I      I      I      I
!!!                  I      I      I      I      I      I
!!!                  I      I      I      I      I      I
!!! nblocks_eta/2=3  I______+______+______+______+______I
!!!                  I      I      I      I      I      I
!!!                  I      I      I      I      I      I
!!!                  I      I      I      I      I      I
!!!                  I______+______+______+______+______I
!!!                  I      I      I      I      I      I
!!!                  I      I      I      I      I      I
!!!                  I      I      I      I      I      I
!!!                  I______I______I______I______I______I
!!!
!!! NGLOB for this doubling layer = 3*5*Volume - ((3-1)*5+(5-1)*3)*Surface + (3-1)*(5-1)*Edge
!!!
!!! 32*NGLLX**3 - 70*NGLLX**2 + 52*NGLLX - 13 -> nb GLL points in a superbrick (Volume)
!!! 8*NGLLX**2-11*NGLLX+4 -> nb GLL points on a superbrick side (Surface)
!!! 2*NGLLX-1 -> nb GLL points on a corner edge of a superbrick (Edge)

!!! for the one layer superbrick :
!!! NGLOB = 28.NGLL^3 - 62.NGLL^2 + 47.NGLL - 12 (Volume)
!!! NGLOB = 6.NGLL^2 - 8.NGLL + 3.NGLL (Surface)
!!! NGLOB = NGLL (Edge)
!!!
!!! those results were obtained by using the script UTILS/doubling_brick/count_nglob_analytical.pl
!!! with an opendx file of the superbrick's geometry

  end subroutine read_compute_parameters

!!!!!! DK DK
!!!!!! DK DK  this section written by Brian Savage, commented out by Dimitri Komatitsch
!!!!!! DK DK  because it is based on the old mesher and therefore does not work with the new
!!!!!! DK DK  mesher. Brian should update it and put it back.
!!!!!! DK DK
!
!----
!
!
!  subroutine auto_ner(WIDTH, NEX_MAX, &
!       NER_CRUST, NER_220_MOHO, NER_400_220, NER_600_400, &
!       NER_670_600, NER_771_670, NER_TOPDDOUBLEPRIME_771, &
!       NER_CMB_TOPDDOUBLEPRIME, NER_TOP_CENTRAL_CUBE_ICB)
!
!    implicit none
!
!    include 'constants.h'
!
!    double precision WIDTH
!    integer NEX_MAX
!    integer NER_CRUST, NER_220_MOHO, NER_400_220, NER_600_400, &
!         NER_670_600, NER_771_670, NER_TOPDDOUBLEPRIME_771, &
!         NER_CMB_TOPDDOUBLEPRIME, NER_TOP_CENTRAL_CUBE_ICB
!
!    integer, parameter                         :: NUM_REGIONS = 13
!    integer, dimension(NUM_REGIONS)            :: scaling
!    double precision, dimension(NUM_REGIONS)   :: radius
!    double precision, dimension(NUM_REGIONS)   :: element_width
!    double precision, dimension(NUM_REGIONS)   :: chunk_width
!    double precision, dimension(NUM_REGIONS-1) :: ratio_top
!    double precision, dimension(NUM_REGIONS-1) :: ratio_bottom
!    integer, dimension(NUM_REGIONS-1)          :: NER
!    integer NER_FLUID
!
!    ! This is PREM in Kilometers
!    radius(1)  = 6371.00d0 ! Surface
!    radius(2)  = 6346.60d0 ! Moho
!    radius(3)  = 6151.00d0 ! 220
!    radius(4)  = 5971.00d0 ! 400
!    radius(5)  = 5771.00d0 ! 600
!    radius(6)  = 5701.00d0 ! 670
!    radius(7)  = 5600.00d0 ! 771
!    radius(8)  = 3630.00d0 ! D''
!    radius(9)  = 3480.00d0 ! CMB
!    radius(10) =    0.00d0 ! Top Double Fluid
!    radius(11) =    0.00d0 ! Bottom Double Fluid
!    radius(12) = 1221.00d0 ! ICB
!    radius(13) = 1071.00d0 ! Top Central Cube
!
!    ! Mesh Doubling
!    scaling(1:1)   = 1
!    scaling(2:5)   = 2
!    scaling(6:11)  = 4
!    scaling(12:13) = 8
!
!    ! Minimum Number of Elements a Region must have
!    NER(:)  = 1
!    NER(2)  = 3
!    NER(3)  = 2
!    NER(12) = 2
!
!    NER_FLUID = 6
!
!    ! Determine the Radius of Top and Bottom of Fluid Doubling Region
!!!!!!!!! DK DK suppressed this    radius(10) = radius(12) + RATIO_TOP_DBL_OC    * (radius(9) - radius(12))
!!!!!!!!! DK DK suppressed this    radius(11) = radius(12) + RATIO_BOTTOM_DBL_OC * (radius(9) - radius(12))
!    radius(10) = 0
!    radius(11) = 0
!
!    ! Horizontal Width of a Chunk
!    chunk_width(:) = WIDTH * (PI/180.0d0) * radius(:)
!
!    ! Horizontal Width of the elements within the chunk
!    element_width(:) = chunk_width(:) / (NEX_MAX / scaling(:))
!
!
!    ! Find the Number of Radial Elements in a region based upon
!    ! the aspect ratio of the elements
!    call auto_optimal_ner(NUM_REGIONS, radius, element_width,NER, ratio_top, ratio_bottom)
!
!    ! Set Output arguments
!    NER_CRUST                = NER(1)
!    NER_220_MOHO             = NER(2)
!    NER_400_220              = NER(3)
!    NER_600_400              = NER(4)
!    NER_670_600              = NER(5)
!    NER_771_670              = NER(6)
!    NER_TOPDDOUBLEPRIME_771  = NER(7)
!    NER_CMB_TOPDDOUBLEPRIME  = NER(8)
!    NER_FLUID                = NER(10)
!    NER_TOP_CENTRAL_CUBE_ICB = NER(12)
!
!  end subroutine auto_ner
!
!!
!!----
!!
!
!  subroutine auto_optimal_ner(NUM_REGIONS, r, ew, NER, rt, rb)
!
!    implicit none
!
!    integer NUM_REGIONS
!    integer,          dimension(NUM_REGIONS-1) :: NER ! Elements per Region
!    double precision, dimension(NUM_REGIONS)   :: r   ! Radius
!    double precision, dimension(NUM_REGIONS)   :: ew  ! Element Width
!    double precision, dimension(NUM_REGIONS-1) :: rt  ! Ratio at Top
!    double precision, dimension(NUM_REGIONS-1) :: rb  ! Ratio at Bottom
!
!    double precision dr, w, ratio, xi, ximin
!    integer ner_test
!    integer i
!
!    ! Find optimal elements per region
!    do i = 1,NUM_REGIONS-1
!       dr = r(i) - r(i+1)              ! Radial Length of Ragion
!       w  = (ew(i) + ew(i+1)) / 2.0d0  ! Average Width of Region
!       ner_test = NER(i)               ! Initial solution
!       ratio = (dr / ner_test) / w     ! Aspect Ratio of Element
!       xi = dabs(ratio - 1.0d0)        ! Aspect Ratio should be near 1.0
!       ximin = 1e7                     ! Initial Minimum
!
!       do while(xi <= ximin)
!          NER(i) = ner_test            ! Found a better solution
!          ximin = xi                   !
!          ner_test = ner_test + 1      ! Increment ner_test and
!          ratio = (dr / ner_test) / w  ! look for a better
!          xi = dabs(ratio - 1.0d0)     ! solution
!       end do
!       rt(i) = dr / NER(i) / ew(i)     ! Find the Ratio of Top
!       rb(i) = dr / NER(i) / ew(i+1)   ! and Bottom for completeness
!    end do
!
!  end subroutine auto_optimal_ner
