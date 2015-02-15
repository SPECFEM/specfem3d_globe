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

  module constants

  include "constants.h"

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

  use constants,only: MAX_STRING_LEN

  implicit none

  ! parameters read from parameter file

  ! seismograms
  integer :: NTSTEP_BETWEEN_OUTPUT_SEISMOS,NTSTEP_BETWEEN_READ_ADJSRC,NTSTEP_BETWEEN_FRAMES, &
             NTSTEP_BETWEEN_OUTPUT_INFO
  double precision :: RECORD_LENGTH_IN_MINUTES
  logical :: RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION
  logical :: OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
             OUTPUT_SEISMOS_ASDF,&
             ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER,&
             SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE

  ! checkpointing/restart
  integer :: NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN

  ! simulation parameters
  integer :: NCHUNKS,SIMULATION_TYPE
  integer :: NEX_XI_read,NEX_ETA_read,NPROC_XI_read,NPROC_ETA_read
  integer :: NOISE_TOMOGRAPHY

  double precision :: ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,&
                      CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH
  logical :: SAVE_MESH_FILES,SAVE_FORWARD

  ! movies
  integer :: MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP
  double precision :: HDUR_MOVIE,MOVIE_TOP_KM,MOVIE_BOTTOM_KM, &
          MOVIE_EAST_DEG,MOVIE_WEST_DEG,MOVIE_NORTH_DEG,&
          MOVIE_SOUTH_DEG
  logical :: MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_COARSE

  ! physical parameters
  logical :: ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS, &
             ATTENUATION

  logical :: ABSORBING_CONDITIONS

  ! file directories
  character(len=MAX_STRING_LEN) :: OUTPUT_FILES,LOCAL_PATH,LOCAL_TMP_PATH,MODEL

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
  logical :: SAVE_REGULAR_KL,ANISOTROPIC_KL,SAVE_TRANSVERSE_KL_ONLY, &
             APPROXIMATE_HESS_KL,USE_FULL_TISO_MANTLE,SAVE_SOURCE_MASK

  ! for simultaneous runs from the same batch job
  integer :: NUMBER_OF_SIMULTANEOUS_RUNS
  logical :: BROADCAST_SAME_MESH_AND_MODEL

  ! GPU simulations
  logical :: GPU_MODE
  integer :: GPU_RUNTIME
  character(len=11) :: GPU_PLATFORM
  character(len=11) :: GPU_DEVICE

  ! adios file output
  logical :: ADIOS_ENABLED
  logical :: ADIOS_FOR_FORWARD_ARRAYS, &
             ADIOS_FOR_MPI_ARRAYS,ADIOS_FOR_ARRAYS_SOLVER,ADIOS_FOR_SOLVER_MESHFILES, &
             ADIOS_FOR_AVS_DX,ADIOS_FOR_KERNELS,ADIOS_FOR_MODELS,ADIOS_FOR_UNDO_ATTENUATION

  end module shared_input_parameters

!
!-------------------------------------------------------------------------------------------------
!

  module shared_compute_parameters

  ! parameters to be computed based upon parameters above read from file

  use constants,only: MAX_NUM_REGIONS,MAX_NUMBER_OF_MESH_LAYERS, &
    NB_SQUARE_CORNERS,NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE

  implicit none

  ! parameters deduced from parameters read from file
  integer :: NPROC,NPROCTOT
  integer :: NPROC_XI,NPROC_ETA

  ! number of time steps
  integer :: NSTEP
  double precision :: DT

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
  integer :: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD
  integer :: ATT1,ATT2,ATT3,ATT4,ATT5

  ! radii of layers
  double precision :: ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R120,R220,R400, &
                      R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
                      R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS, &
                      RMOHO_FICTITIOUS_IN_MESHER,R80_FICTITIOUS_IN_MESHER

  ! movies
  double precision :: MOVIE_TOP,MOVIE_BOTTOM,MOVIE_EAST,MOVIE_WEST,&
                      MOVIE_NORTH,MOVIE_SOUTH
  ! model flags
  integer :: REFERENCE_1D_MODEL,THREE_D_MODEL
  logical :: TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
             CRUSTAL,ONE_CRUST,ISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE, &
             CEM_REQUEST,CEM_ACCEPT
  logical :: ATTENUATION_3D
  logical :: INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE
  logical :: EMULATE_ONLY = .false.


! honor PREM Moho or not
! doing so drastically reduces the stability condition and therefore the time step
  logical :: HONOR_1D_SPHERICAL_MOHO,CASE_3D

  ! this for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC2D_XI,NSPEC2D_ETA
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC2D_BOTTOM,NSPEC2D_TOP
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC1D_RADIAL

  integer, dimension(MAX_NUM_REGIONS) :: NGLOB
  integer, dimension(MAX_NUM_REGIONS) :: NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX
  integer, dimension(MAX_NUM_REGIONS) :: NGLOB1D_RADIAL

  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index

  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling

  integer :: ratio_divide_central_cube

  ! for the cut doublingbrick improvement
  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA
  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA

  end module shared_compute_parameters

!
!-------------------------------------------------------------------------------------------------
!

  module shared_parameters

  use shared_input_parameters
  use shared_compute_parameters

  implicit none

  end module shared_parameters

