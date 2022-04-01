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

  module constants

  include "constants.h"

  ! proc number for MPI process
  integer :: myrank

  ! a negative initial value is a convention that indicates that groups
  ! (i.e. sub-communicators, one per run) are off by default
  integer :: mygroup = -1

  ! if doing simultaneous runs for the same mesh and model, see who
  ! should read the mesh and the model and broadcast it to others
  ! we put a default value here
  logical :: I_should_read_the_database = .true.

  end module constants

!
!-------------------------------------------------------------------------------------------------
!

  module shared_input_parameters

! holds input parameters given in DATA/Par_file

  use constants, only: MAX_STRING_LEN

  implicit none

  ! parameters read from parameter file

  ! seismograms
  integer :: NTSTEP_BETWEEN_OUTPUT_SEISMOS,NTSTEP_BETWEEN_READ_ADJSRC,NTSTEP_BETWEEN_FRAMES, &
             NTSTEP_BETWEEN_OUTPUT_INFO,NTSTEP_BETWEEN_OUTPUT_SAMPLE

  double precision :: RECORD_LENGTH_IN_MINUTES

  logical :: RECEIVERS_CAN_BE_BURIED
  logical :: OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
             OUTPUT_SEISMOS_ASDF,OUTPUT_SEISMOS_3D_ARRAY, &
             ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MAIN, &
             SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,READ_ADJSRC_ASDF

  logical :: SAVE_SEISMOGRAMS_STRAIN ! option under development (see Hom Nath commit 2deb0fa89), no functionality implemented yet
  logical :: SAVE_SEISMOGRAMS_IN_ADJOINT_RUN

  ! sources
  logical :: USE_FORCE_POINT_SOURCE
  logical :: USE_MONOCHROMATIC_CMT_SOURCE,PRINT_SOURCE_TIME_FUNCTION

  ! checkpointing/restart
  integer :: NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN

  ! simulation parameters
  integer :: NCHUNKS,SIMULATION_TYPE
  integer :: NEX_XI_read,NEX_ETA_read,NPROC_XI_read,NPROC_ETA_read
  integer :: NOISE_TOMOGRAPHY

  double precision :: ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES, &
                      CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH

  logical :: SAVE_MESH_FILES,SAVE_FORWARD

  ! movies
  integer :: MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP
  double precision :: HDUR_MOVIE,MOVIE_TOP_KM,MOVIE_BOTTOM_KM, &
          MOVIE_EAST_DEG,MOVIE_WEST_DEG,MOVIE_NORTH_DEG, &
          MOVIE_SOUTH_DEG
  logical :: MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_COARSE

  ! physical parameters
  logical :: ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS, &
             ATTENUATION

  ! regional mesh cut-off
  logical :: REGIONAL_MESH_CUTOFF
  ! regional mesh cut-off depth (in km)
  ! possible selections: 24.4d0, 80.d0, 220.d0, 400.d0, 600.d0, 670.d0, 771.d0
  double precision :: REGIONAL_MESH_CUTOFF_DEPTH = 400.d0
  ! regional mesh cut-off w/ a second doubling layer below 220km interface
  logical :: REGIONAL_MESH_ADD_2ND_DOUBLING = .false.

  ! absorbing boundary conditions
  logical :: ABSORBING_CONDITIONS

  ! absorbing sponge layer
  logical :: ABSORB_USING_GLOBAL_SPONGE
  double precision :: SPONGE_LATITUDE_IN_DEGREES,SPONGE_LONGITUDE_IN_DEGREES,SPONGE_RADIUS_IN_DEGREES

  ! file directories
  character(len=MAX_STRING_LEN) :: OUTPUT_FILES
  character(len=MAX_STRING_LEN) :: MODEL,MODEL_NAME
  character(len=MAX_STRING_LEN) :: LOCAL_PATH,LOCAL_TMP_PATH

  ! attenuation parameters
  logical :: UNDO_ATTENUATION,PARTIAL_PHYS_DISPERSION_ONLY

  ! exact (full) undoing of attenuation
  ! How much memory (in GB) is installed on the machine per CPU core (or per GPU card or per INTEL MIC Phi board)
  double precision :: MEMORY_INSTALLED_PER_CORE_IN_GB

  ! exact (full) undoing of attenuation
  ! What percentage of this total do you allow us to use for arrays to undo attenuation, keeping in mind that you
  ! need to leave some memory available for the GNU/Linux system to run
  double precision :: PERCENT_OF_MEM_TO_USE_PER_CORE

  ! LDD Runge-Kutta time scheme
  logical :: USE_LDDRK,INCREASE_CFL_FOR_LDDRK
  double precision :: RATIO_BY_WHICH_TO_INCREASE_IT

  logical :: EXACT_MASS_MATRIX_FOR_ROTATION

  ! adjoint kernels
  logical :: SAVE_REGULAR_KL, &
             ANISOTROPIC_KL, &
             SAVE_TRANSVERSE_KL_ONLY, &
             SAVE_AZIMUTHAL_ANISO_KL_ONLY, &
             APPROXIMATE_HESS_KL

  logical :: USE_FULL_TISO_MANTLE,SAVE_SOURCE_MASK

  logical :: STEADY_STATE_KERNEL
  double precision :: STEADY_STATE_LENGTH_IN_MINUTES

  ! for simultaneous runs from the same batch job
  integer :: NUMBER_OF_SIMULTANEOUS_RUNS
  logical :: BROADCAST_SAME_MESH_AND_MODEL

  ! GPU simulations
  integer :: GPU_RUNTIME
  character(len=128) :: GPU_PLATFORM
  character(len=128) :: GPU_DEVICE
  logical :: GPU_MODE

  ! adios file output
  logical :: ADIOS_ENABLED
  logical :: ADIOS_FOR_FORWARD_ARRAYS, &
             ADIOS_FOR_MPI_ARRAYS,ADIOS_FOR_ARRAYS_SOLVER,ADIOS_FOR_SOLVER_MESHFILES, &
             ADIOS_FOR_AVS_DX,ADIOS_FOR_KERNELS,ADIOS_FOR_MODELS,ADIOS_FOR_UNDO_ATTENUATION

  ! (optional) parameters
  double precision :: USER_DT = -1.0  ! negative values will be ignored by default
  integer :: USER_NSTEP = -1          ! negative to ignore by default

  end module shared_input_parameters

!
!-------------------------------------------------------------------------------------------------
!

  module shared_compute_parameters

  ! parameters to be computed based upon parameters above read from file

  use constants, only: MAX_NUM_REGIONS,MAX_NUMBER_OF_MESH_LAYERS, &
    NB_SQUARE_CORNERS,NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE,MAX_STRING_LEN

  use constants, only: IPLANET_EARTH,IPLANET_MARS,IPLANET_MOON

  use constants, only: &
    EARTH_R,EARTH_R_KM, &
    EARTH_STANDARD_GRAVITY,EARTH_RHOAV, &
    EARTH_NX_BATHY,EARTH_NY_BATHY,EARTH_TOPO_MAXIMUM,EARTH_TOPO_MINIMUM, &
    EARTH_PATHNAME_TOPO_FILE,EARTH_RESOLUTION_TOPO_FILE, &
    EARTH_HOURS_PER_DAY,EARTH_SECONDS_PER_HOUR,EARTH_ONE_MINUS_F_SQUARED, &
    EARTH_HONOR_DEEP_MOHO,EARTH_R_DEEPEST_CRUST,EARTH_REGIONAL_MOHO_MESH, &
    EARTH_MAX_RATIO_CRUST_STRETCHING,EARTH_RMOHO_STRETCH_ADJUSTMENT,EARTH_R80_STRETCH_ADJUSTMENT

  implicit none

  ! parameters deduced from parameters read from file
  integer :: NPROC,NPROCTOT
  integer :: NPROC_XI,NPROC_ETA

  ! number of time steps
  integer :: NSTEP
  double precision :: DT

  ! number of steady state time steps
  integer :: NSTEP_STEADY_STATE

  ! shortest minimum period resolved by mesh (empirical formula)
  double precision :: T_min_period

  ! number of sources given in CMTSOLUTION file
  integer :: NSOURCES

  ! number of elements
  integer :: NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer :: NER_CRUST, &
             NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
             NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
             NER_TOP_CENTRAL_CUBE_ICB, &
             NEX_XI,NEX_ETA

  ! attenuation
  ! attenuation period band min/max
  double precision :: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD
  ! logarithmic center frequency (center of attenuation band)
  double precision :: ATT_F_C_SOURCE
  ! attenuation array sizes
  integer :: ATT1,ATT2,ATT3,ATT4,ATT5

  ! radii of layers
  double precision :: ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R120,R220,R400, &
                      R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
                      R_CENTRAL_CUBE, &
                      RMOHO_FICTITIOUS_IN_MESHER,R80_FICTITIOUS_IN_MESHER

  ! densities
  double precision :: RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS

  ! movies
  double precision :: MOVIE_TOP,MOVIE_BOTTOM,MOVIE_EAST,MOVIE_WEST, &
                      MOVIE_NORTH,MOVIE_SOUTH
  ! model flags
  integer :: REFERENCE_1D_MODEL,REFERENCE_CRUSTAL_MODEL
  integer :: THREE_D_MODEL,THREE_D_MODEL_IC

  logical :: TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
             CRUSTAL,ONE_CRUST
  logical :: MODEL_3D_MANTLE_PERTUBATIONS,HETEROGEN_3D_MANTLE
  logical :: CEM_REQUEST,CEM_ACCEPT

  logical :: MODEL_GLL
  integer :: MODEL_GLL_TYPE

  logical :: ATTENUATION_3D
  logical :: ATTENUATION_GLL
  logical :: INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE

  logical :: EMULATE_ONLY = .false.

! honor PREM Moho or not
! doing so drastically reduces the stability condition and therefore the time step
  logical :: HONOR_1D_SPHERICAL_MOHO,CASE_3D

  ! this for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC_REGIONS
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC2D_XI,NSPEC2D_ETA
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC2D_BOTTOM,NSPEC2D_TOP
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC1D_RADIAL

  integer, dimension(MAX_NUM_REGIONS) :: NGLOB_REGIONS
  integer, dimension(MAX_NUM_REGIONS) :: NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX
  integer, dimension(MAX_NUM_REGIONS) :: NGLOB1D_RADIAL

  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner_mesh_layers
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ratio_sampling_array
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index

  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling

  integer :: ratio_divide_central_cube

  ! for the cut doublingbrick improvement
  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA
  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA

! default planet
!
! note: for different planets, we will re-set R_EARTH, RHOAV, .. after reading the Par_file
!       to avoid problems when they are used to non-dimensionalize all parameters.
!
!       see routine get_model_planet_constants() in get_model_parameters.F90
!
  ! default planet Earth
  integer :: PLANET_TYPE = IPLANET_EARTH

  ! gravity
  double precision :: STANDARD_GRAVITY = EARTH_STANDARD_GRAVITY

  ! flattening / eccentricity
  double precision :: ONE_MINUS_F_SQUARED = EARTH_ONE_MINUS_F_SQUARED

  ! topo
  character (len=MAX_STRING_LEN) :: PATHNAME_TOPO_FILE = EARTH_PATHNAME_TOPO_FILE
  integer :: NX_BATHY = EARTH_NX_BATHY
  integer :: NY_BATHY = EARTH_NY_BATHY
  double precision :: RESOLUTION_TOPO_FILE = EARTH_RESOLUTION_TOPO_FILE
  integer :: TOPO_MINIMUM = EARTH_TOPO_MINIMUM
  integer :: TOPO_MAXIMUM = EARTH_TOPO_MAXIMUM

  ! planet constants
  ! radius of globe
  ! we still use R_EARTH, but will mostly shift to R_PLANET to allow for different planets
  double precision :: R_EARTH = EARTH_R    ! default: earth
  double precision :: R_EARTH_KM = EARTH_R_KM

  ! physical surface radius
  double precision :: R_PLANET = EARTH_R
  double precision :: R_PLANET_KM = EARTH_R / 1000.d0

  ! average density
  double precision :: RHOAV = EARTH_RHOAV  ! default: earth density

  ! crust
  double precision :: R_DEEPEST_CRUST = EARTH_R_DEEPEST_CRUST

  ! rotation
  double precision :: HOURS_PER_DAY = EARTH_HOURS_PER_DAY
  double precision :: SECONDS_PER_HOUR = EARTH_SECONDS_PER_HOUR

  ! mesh
  double precision :: MAX_RATIO_CRUST_STRETCHING = EARTH_MAX_RATIO_CRUST_STRETCHING
  double precision :: RMOHO_STRETCH_ADJUSTMENT = EARTH_RMOHO_STRETCH_ADJUSTMENT
  double precision :: R80_STRETCH_ADJUSTMENT = EARTH_R80_STRETCH_ADJUSTMENT
  logical :: REGIONAL_MOHO_MESH = EARTH_REGIONAL_MOHO_MESH
  logical :: HONOR_DEEP_MOHO = EARTH_HONOR_DEEP_MOHO

  end module shared_compute_parameters

!
!-------------------------------------------------------------------------------------------------
!

  module shared_parameters

  use shared_input_parameters
  use shared_compute_parameters

  implicit none

  end module shared_parameters

