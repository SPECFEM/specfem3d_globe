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

  subroutine broadcast_computed_parameters()

  use constants, only: myrank
  use shared_parameters

  implicit none

  ! local parameters
  ! broadcast parameter arrays
  integer, parameter :: nparam_i = 47
  integer, dimension(nparam_i) :: bcast_integer

  integer, parameter :: nparam_l = 67
  logical, dimension(nparam_l) :: bcast_logical

  integer, parameter :: nparam_dp = 38
  double precision, dimension(nparam_dp) :: bcast_double_precision

  ! initializes containers
  bcast_integer(:) = 0
  bcast_logical(:) = .false.
  bcast_double_precision(:) = 0.d0

  ! master process prepares broadcasting arrays
  if (myrank == 0) then
    ! simple way to pass parameters in arrays from master to all other processes
    ! rather than single values one by one to reduce MPI communication calls:
    ! sets up broadcasting array
    bcast_integer = (/ &
            NER_CRUST, &
            NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
            NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
            NER_TOP_CENTRAL_CUBE_ICB, &
            NEX_XI,NEX_ETA, &
            NPROC_XI,NPROC_ETA, &
            NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
            NTSTEP_BETWEEN_READ_ADJSRC, &
            NSTEP,NSOURCES, &
            NTSTEP_BETWEEN_FRAMES, &
            NTSTEP_BETWEEN_OUTPUT_INFO, &
            NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS, &
            SIMULATION_TYPE, &
            REFERENCE_1D_MODEL,REFERENCE_CRUSTAL_MODEL, &
            THREE_D_MODEL,THREE_D_MODEL_IC, &
            NPROC,NPROCTOT, &
            NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube, &
            MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP, &
            NOISE_TOMOGRAPHY, &
            ATT1,ATT2,ATT3,ATT4,ATT5, &
            GPU_RUNTIME,NUMBER_OF_SIMULTANEOUS_RUNS, &
            MODEL_GLL_TYPE,USER_NSTEP /)

    bcast_logical = (/ &
            TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
            CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,MODEL_3D_MANTLE_PERTUBATIONS,HETEROGEN_3D_MANTLE, &
            TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,ATTENUATION_3D, ATTENUATION_GLL, &
            RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
            SAVE_MESH_FILES,ATTENUATION, &
            ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD,CASE_3D, &
            CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,SAVE_ALL_SEISMOS_IN_ONE_FILE, &
            HONOR_1D_SPHERICAL_MOHO,MOVIE_COARSE, &
            USE_FORCE_POINT_SOURCE,SAVE_SEISMOGRAMS_STRAIN,SAVE_SEISMOGRAMS_IN_ADJOINT_RUN, &
            OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
            OUTPUT_SEISMOS_ASDF, &
            ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER,USE_BINARY_FOR_LARGE_FILE, &
            READ_ADJSRC_ASDF,SAVE_REGULAR_KL, &
            PARTIAL_PHYS_DISPERSION_ONLY,UNDO_ATTENUATION, &
            USE_LDDRK,INCREASE_CFL_FOR_LDDRK, &
            ANISOTROPIC_KL, &
            SAVE_TRANSVERSE_KL_ONLY, &
            SAVE_AZIMUTHAL_ANISO_KL_ONLY, &
            APPROXIMATE_HESS_KL, &
            USE_FULL_TISO_MANTLE,SAVE_SOURCE_MASK, &
            EXACT_MASS_MATRIX_FOR_ROTATION, &
            GPU_MODE, &
            ADIOS_ENABLED,ADIOS_FOR_FORWARD_ARRAYS, &
            ADIOS_FOR_MPI_ARRAYS,ADIOS_FOR_ARRAYS_SOLVER, &
            ADIOS_FOR_SOLVER_MESHFILES,ADIOS_FOR_AVS_DX, &
            ADIOS_FOR_KERNELS,ADIOS_FOR_MODELS,ADIOS_FOR_UNDO_ATTENUATION, &
            CEM_REQUEST,CEM_ACCEPT,BROADCAST_SAME_MESH_AND_MODEL,MODEL_GLL /)

    bcast_double_precision = (/ &
            DT, &
            MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD, &
            ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
            CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
            RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
            R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
            MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH, &
            RMOHO_FICTITIOUS_IN_MESHER,RATIO_BY_WHICH_TO_INCREASE_IT, &
            MEMORY_INSTALLED_PER_CORE_IN_GB,PERCENT_OF_MEM_TO_USE_PER_CORE, &
            RECORD_LENGTH_IN_MINUTES, USER_DT /)
  endif

  ! broadcasts the information read on the master to the nodes
  call bcast_all_i(bcast_integer,nparam_i)
  call bcast_all_l(bcast_logical,nparam_l)
  call bcast_all_dp(bcast_double_precision,nparam_dp)

  ! broadcasts non-single value parameters
  call bcast_all_ch(LOCAL_PATH,MAX_STRING_LEN)
  call bcast_all_ch(LOCAL_TMP_PATH,MAX_STRING_LEN)

  call bcast_all_ch(MODEL,MAX_STRING_LEN)
  call bcast_all_ch(MODEL_NAME,MAX_STRING_LEN)

  call bcast_all_ch(GPU_PLATFORM,128)
  call bcast_all_ch(GPU_DEVICE,128)

  call bcast_all_i(ner_mesh_layers,MAX_NUMBER_OF_MESH_LAYERS)
  call bcast_all_i(ratio_sampling_array,MAX_NUMBER_OF_MESH_LAYERS)
  call bcast_all_i(doubling_index,MAX_NUMBER_OF_MESH_LAYERS)

  call bcast_all_dp(r_bottom,MAX_NUMBER_OF_MESH_LAYERS)
  call bcast_all_dp(r_top,MAX_NUMBER_OF_MESH_LAYERS)
  call bcast_all_dp(rmins,MAX_NUMBER_OF_MESH_LAYERS)
  call bcast_all_dp(rmaxs,MAX_NUMBER_OF_MESH_LAYERS)

  call bcast_all_l(this_region_has_a_doubling,MAX_NUMBER_OF_MESH_LAYERS)

  call bcast_all_i(NSPEC_REGIONS,MAX_NUM_REGIONS)
  call bcast_all_i(NSPEC2D_XI,MAX_NUM_REGIONS)
  call bcast_all_i(NSPEC2D_ETA,MAX_NUM_REGIONS)
  call bcast_all_i(NSPEC2DMAX_XMIN_XMAX,MAX_NUM_REGIONS)
  call bcast_all_i(NSPEC2DMAX_YMIN_YMAX,MAX_NUM_REGIONS)
  call bcast_all_i(NSPEC2D_BOTTOM,MAX_NUM_REGIONS)
  call bcast_all_i(NSPEC2D_TOP,MAX_NUM_REGIONS)
  call bcast_all_i(NSPEC1D_RADIAL,MAX_NUM_REGIONS)
  call bcast_all_i(NGLOB1D_RADIAL,MAX_NUM_REGIONS)
  call bcast_all_i(NGLOB2DMAX_XMIN_XMAX,MAX_NUM_REGIONS)
  call bcast_all_i(NGLOB2DMAX_YMIN_YMAX,MAX_NUM_REGIONS)
  call bcast_all_i(NGLOB_REGIONS,MAX_NUM_REGIONS)
  call bcast_all_i(DIFF_NSPEC1D_RADIAL,NB_SQUARE_CORNERS*NB_CUT_CASE)
  call bcast_all_i(DIFF_NSPEC2D_ETA,NB_SQUARE_EDGES_ONEDIR*NB_CUT_CASE)
  call bcast_all_i(DIFF_NSPEC2D_XI,NB_SQUARE_EDGES_ONEDIR*NB_CUT_CASE)

  ! broadcasts shared parameters for planet setup
  call bcast_all_singlei(PLANET_TYPE)
  call bcast_all_singledp(R_EARTH)
  call bcast_all_singledp(R_EARTH_KM)
  call bcast_all_singledp(RHOAV)
  call bcast_all_singledp(STANDARD_GRAVITY)
  call bcast_all_singledp(ONE_MINUS_F_SQUARED)
  call bcast_all_ch(PATHNAME_TOPO_FILE,MAX_STRING_LEN)
  call bcast_all_singledp(RESOLUTION_TOPO_FILE)
  call bcast_all_singlei(NX_BATHY)
  call bcast_all_singlei(NY_BATHY)
  call bcast_all_singlei(TOPO_MINIMUM)
  call bcast_all_singlei(TOPO_MAXIMUM)
  call bcast_all_singledp(R_DEEPEST_CRUST)
  call bcast_all_singledp(HOURS_PER_DAY)
  call bcast_all_singledp(SECONDS_PER_HOUR)
  call bcast_all_singledp(MAX_RATIO_CRUST_STRETCHING)
  call bcast_all_singledp(RMOHO_STRETCH_ADJUSTMENT)
  call bcast_all_singledp(R80_STRETCH_ADJUSTMENT)
  call bcast_all_singlel(REGIONAL_MOHO_MESH)
  call bcast_all_singlel(HONOR_DEEP_MOHO)

  ! empirical minimum period resolved estimation
  call bcast_all_singledp(T_min)

  ! non-master processes set their parameters
  if (myrank /= 0) then

    ! please, be careful with ordering and counting here
    ! integers
    NER_CRUST = bcast_integer(1)
    NER_80_MOHO = bcast_integer(2)
    NER_220_80 = bcast_integer(3)
    NER_400_220 = bcast_integer(4)
    NER_600_400 = bcast_integer(5)
    NER_670_600 = bcast_integer(6)
    NER_771_670 = bcast_integer(7)
    NER_TOPDDOUBLEPRIME_771 = bcast_integer(8)
    NER_CMB_TOPDDOUBLEPRIME = bcast_integer(9)
    NER_OUTER_CORE = bcast_integer(10)
    NER_TOP_CENTRAL_CUBE_ICB = bcast_integer(11)
    NEX_XI = bcast_integer(12)
    NEX_ETA = bcast_integer(13)
    NPROC_XI = bcast_integer(14)
    NPROC_ETA = bcast_integer(15)
    NTSTEP_BETWEEN_OUTPUT_SEISMOS = bcast_integer(16)
    NTSTEP_BETWEEN_READ_ADJSRC = bcast_integer(17)
    NSTEP = bcast_integer(18)
    NSOURCES = bcast_integer(19)
    NTSTEP_BETWEEN_FRAMES = bcast_integer(20)
    NTSTEP_BETWEEN_OUTPUT_INFO = bcast_integer(21)
    NUMBER_OF_RUNS = bcast_integer(22)
    NUMBER_OF_THIS_RUN = bcast_integer(23)
    NCHUNKS = bcast_integer(24)
    SIMULATION_TYPE = bcast_integer(25)
    REFERENCE_1D_MODEL = bcast_integer(26)
    REFERENCE_CRUSTAL_MODEL = bcast_integer(27)
    THREE_D_MODEL = bcast_integer(28)
    THREE_D_MODEL_IC = bcast_integer(29)
    NPROC = bcast_integer(30)
    NPROCTOT = bcast_integer(31)
    NEX_PER_PROC_XI = bcast_integer(32)
    NEX_PER_PROC_ETA = bcast_integer(33)
    ratio_divide_central_cube = bcast_integer(34)
    MOVIE_VOLUME_TYPE = bcast_integer(35)
    MOVIE_START = bcast_integer(36)
    MOVIE_STOP = bcast_integer(37)
    NOISE_TOMOGRAPHY = bcast_integer(38)
    ATT1 = bcast_integer(39)
    ATT2 = bcast_integer(40)
    ATT3 = bcast_integer(41)
    ATT4 = bcast_integer(42)
    ATT5 = bcast_integer(43)
    GPU_RUNTIME = bcast_integer(44)
    NUMBER_OF_SIMULTANEOUS_RUNS = bcast_integer(45)
    MODEL_GLL_TYPE  = bcast_integer(46)
    USER_NSTEP = bcast_integer(47)

    ! logicals
    TRANSVERSE_ISOTROPY = bcast_logical(1)
    ANISOTROPIC_3D_MANTLE = bcast_logical(2)
    ANISOTROPIC_INNER_CORE = bcast_logical(3)
    CRUSTAL = bcast_logical(4)
    ELLIPTICITY = bcast_logical(5)
    GRAVITY = bcast_logical(6)
    ONE_CRUST = bcast_logical(7)
    ROTATION = bcast_logical(8)
    MODEL_3D_MANTLE_PERTUBATIONS = bcast_logical(9)
    HETEROGEN_3D_MANTLE = bcast_logical(10)
    TOPOGRAPHY = bcast_logical(11)
    OCEANS = bcast_logical(12)
    MOVIE_SURFACE = bcast_logical(13)
    MOVIE_VOLUME = bcast_logical(14)
    ATTENUATION_3D = bcast_logical(15)
    ATTENUATION_GLL = bcast_logical(16)
    RECEIVERS_CAN_BE_BURIED = bcast_logical(17)
    PRINT_SOURCE_TIME_FUNCTION = bcast_logical(18)
    SAVE_MESH_FILES = bcast_logical(19)
    ATTENUATION = bcast_logical(20)
    ABSORBING_CONDITIONS = bcast_logical(21)
    INCLUDE_CENTRAL_CUBE = bcast_logical(22)
    INFLATE_CENTRAL_CUBE = bcast_logical(23)
    SAVE_FORWARD = bcast_logical(24)
    CASE_3D = bcast_logical(25)
    CUT_SUPERBRICK_XI = bcast_logical(26)
    CUT_SUPERBRICK_ETA = bcast_logical(27)
    SAVE_ALL_SEISMOS_IN_ONE_FILE = bcast_logical(28)
    HONOR_1D_SPHERICAL_MOHO = bcast_logical(29)
    MOVIE_COARSE = bcast_logical(30)
    USE_FORCE_POINT_SOURCE= bcast_logical(31)
    SAVE_SEISMOGRAMS_STRAIN= bcast_logical(32)
    SAVE_SEISMOGRAMS_IN_ADJOINT_RUN= bcast_logical(33)
    OUTPUT_SEISMOS_ASCII_TEXT= bcast_logical(34)
    OUTPUT_SEISMOS_SAC_ALPHANUM= bcast_logical(35)
    OUTPUT_SEISMOS_SAC_BINARY= bcast_logical(36)
    OUTPUT_SEISMOS_ASDF = bcast_logical(37)
    ROTATE_SEISMOGRAMS_RT= bcast_logical(38)
    WRITE_SEISMOGRAMS_BY_MASTER= bcast_logical(39)
    USE_BINARY_FOR_LARGE_FILE= bcast_logical(40)
    READ_ADJSRC_ASDF = bcast_logical(41)
    SAVE_REGULAR_KL = bcast_logical(42)
    PARTIAL_PHYS_DISPERSION_ONLY = bcast_logical(43)
    UNDO_ATTENUATION = bcast_logical(44)
    USE_LDDRK = bcast_logical(45)
    INCREASE_CFL_FOR_LDDRK = bcast_logical(46)
    ANISOTROPIC_KL = bcast_logical(47)
    SAVE_TRANSVERSE_KL_ONLY = bcast_logical(48)
    SAVE_AZIMUTHAL_ANISO_KL_ONLY = bcast_logical(49)
    APPROXIMATE_HESS_KL = bcast_logical(50)
    USE_FULL_TISO_MANTLE = bcast_logical(51)
    SAVE_SOURCE_MASK = bcast_logical(52)
    EXACT_MASS_MATRIX_FOR_ROTATION = bcast_logical(53)
    GPU_MODE = bcast_logical(54)
    ADIOS_ENABLED = bcast_logical(55)
    ADIOS_FOR_FORWARD_ARRAYS = bcast_logical(56)
    ADIOS_FOR_MPI_ARRAYS = bcast_logical(57)
    ADIOS_FOR_ARRAYS_SOLVER = bcast_logical(58)
    ADIOS_FOR_SOLVER_MESHFILES = bcast_logical(59)
    ADIOS_FOR_AVS_DX = bcast_logical(60)
    ADIOS_FOR_KERNELS = bcast_logical(61)
    ADIOS_FOR_MODELS = bcast_logical(62)
    ADIOS_FOR_UNDO_ATTENUATION = bcast_logical(63)
    CEM_REQUEST = bcast_logical(64)
    CEM_ACCEPT = bcast_logical(65)
    BROADCAST_SAME_MESH_AND_MODEL = bcast_logical(66)
    MODEL_GLL = bcast_logical(67)

    ! double precisions
    DT = bcast_double_precision(1)
    MIN_ATTENUATION_PERIOD = bcast_double_precision(2)
    MAX_ATTENUATION_PERIOD = bcast_double_precision(3)
    ANGULAR_WIDTH_XI_IN_DEGREES = bcast_double_precision(4)
    ANGULAR_WIDTH_ETA_IN_DEGREES = bcast_double_precision(5)
    CENTER_LONGITUDE_IN_DEGREES = bcast_double_precision(6)
    CENTER_LATITUDE_IN_DEGREES = bcast_double_precision(7)
    GAMMA_ROTATION_AZIMUTH = bcast_double_precision(8)
    ROCEAN = bcast_double_precision(9)
    RMIDDLE_CRUST = bcast_double_precision(10)
    RMOHO = bcast_double_precision(11)
    R80 = bcast_double_precision(12)
    R120 = bcast_double_precision(13)
    R220 = bcast_double_precision(14)
    R400 = bcast_double_precision(15)
    R600 = bcast_double_precision(16)
    R670 = bcast_double_precision(17)
    R771 = bcast_double_precision(18)
    RTOPDDOUBLEPRIME = bcast_double_precision(19)
    RCMB = bcast_double_precision(20)
    RICB = bcast_double_precision(21)
    R_CENTRAL_CUBE = bcast_double_precision(22)
    RHO_TOP_OC = bcast_double_precision(23)
    RHO_BOTTOM_OC = bcast_double_precision(24)
    RHO_OCEANS = bcast_double_precision(25)
    HDUR_MOVIE = bcast_double_precision(26)
    MOVIE_TOP = bcast_double_precision(27)
    MOVIE_BOTTOM = bcast_double_precision(28)
    MOVIE_WEST = bcast_double_precision(29)
    MOVIE_EAST = bcast_double_precision(30)
    MOVIE_NORTH = bcast_double_precision(31)
    MOVIE_SOUTH = bcast_double_precision(32)
    RMOHO_FICTITIOUS_IN_MESHER = bcast_double_precision(33)
    RATIO_BY_WHICH_TO_INCREASE_IT = bcast_double_precision(34)
    MEMORY_INSTALLED_PER_CORE_IN_GB = bcast_double_precision(35)
    PERCENT_OF_MEM_TO_USE_PER_CORE = bcast_double_precision(36)
    RECORD_LENGTH_IN_MINUTES = bcast_double_precision(37)
    USER_DT = bcast_double_precision(38)

  endif

  end subroutine broadcast_computed_parameters

