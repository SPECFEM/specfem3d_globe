!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
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


  subroutine read_compute_parameters()

! call this read_compute_parameters() routine to read the Par_file and get the necessary parameter setup for the computations
!
! note: we split the read_compute_parameters() routine into two separate routine calls
!       to make it easier for testing the reading of the parameter file, i.e., read_parameter_file(), and
!       calling the compute_parameters routine, i.e., rcp_set_compute_parameters(), by unit testing:
!       > make tests
!       for example for test programs in tests/meshfem3D/

  implicit none

  ! reads in Par_file values
  call read_parameter_file()

  ! sets parameters for computation
  call rcp_set_compute_parameters()

  ! sets mesh parameters
  call rcp_set_mesh_parameters()

  end subroutine read_compute_parameters

!
!-------------------------------------------------------------------------------------------------
!

  subroutine rcp_set_compute_parameters()

! sets parameters for computation based on Par_file settings

  use constants
  use shared_parameters

  implicit none

  ! local parameters
  double precision :: T_min_res
  character(len=MAX_STRING_LEN) :: path_to_add

  ! chunk sizes
  if (NCHUNKS == 6) then
    ! global simulations
    ANGULAR_WIDTH_XI_IN_DEGREES = 90.d0
    ANGULAR_WIDTH_ETA_IN_DEGREES = 90.d0
    CENTER_LATITUDE_IN_DEGREES = 0.d0
    CENTER_LONGITUDE_IN_DEGREES = 0.d0
    GAMMA_ROTATION_AZIMUTH = 0.d0
  endif

  ! include central cube or not
  ! use regular cubed sphere instead of cube for large distances
  if (NCHUNKS == 6) then
    INCLUDE_CENTRAL_CUBE = .true.
    INFLATE_CENTRAL_CUBE = .false.
  else
    INCLUDE_CENTRAL_CUBE = .false.
    INFLATE_CENTRAL_CUBE = .true.
  endif

  if (.not. EMULATE_ONLY) then
    NEX_XI = NEX_XI_read
    NEX_ETA = NEX_ETA_read
    NPROC_XI = NPROC_XI_read
    NPROC_ETA = NPROC_ETA_read
  else
    ! this is used in UTILS/estimate_best_values_runs.f90 only, to estimate memory use
    NEX_ETA = NEX_XI
    NPROC_ETA = NPROC_XI
  endif

  ! make sure single run simulation setting valid
  if (NUMBER_OF_RUNS == 1) NUMBER_OF_THIS_RUN = 1

  ! re-sets regional mesh cut-off
  if (NCHUNKS == 6) then
    REGIONAL_MESH_CUTOFF = .false.
    REGIONAL_MESH_CUTOFF_DEPTH = 400.d0
    REGIONAL_MESH_ADD_2ND_DOUBLING = .false.
  endif

  ! sponge layer
  if (.not. ABSORB_USING_GLOBAL_SPONGE) then
    SPONGE_LATITUDE_IN_DEGREES = 0.d0
    SPONGE_LONGITUDE_IN_DEGREES = 0.d0
    SPONGE_RADIUS_IN_DEGREES = 0.d0
  endif

  ! steady state simulations
  if (.not. STEADY_STATE_KERNEL) then
    STEADY_STATE_LENGTH_IN_MINUTES = 0.d0
  endif

  ! ignore EXACT_MASS_MATRIX_FOR_ROTATION if rotation is not included in the simulations
  if (.not. ROTATION) EXACT_MASS_MATRIX_FOR_ROTATION = .false.

  ! re-sets attenuation flags
  if (.not. ATTENUATION) then
    ! turns off PARTIAL_PHYS_DISPERSION_ONLY when ATTENUATION is off in the Par_file
    PARTIAL_PHYS_DISPERSION_ONLY = .false.
  endif

  ! re-sets ADIOS flags
  if (.not. ADIOS_ENABLED) then
    ADIOS_FOR_FORWARD_ARRAYS = .false.
    ADIOS_FOR_MPI_ARRAYS = .false.
    ADIOS_FOR_ARRAYS_SOLVER = .false.
    ADIOS_FOR_SOLVER_MESHFILES = .false.
    ADIOS_FOR_AVS_DX = .false.
    ADIOS_FOR_KERNELS = .false.
    ADIOS_FOR_MODELS = .false.
    ADIOS_FOR_UNDO_ATTENUATION = .false.
  endif

  ! ADIOS is very useful for very large simulations (say using 2000 MPI tasks or more)
  ! but slows down the code if used for simulations that are small or medium size, because of the overhead any library has.
  if (ADIOS_ENABLED .and. NCHUNKS * NPROC_XI_read * NPROC_ETA_read < 2000 .and. myrank == 0) then
    print *
    print *,'**************'
    print *,'**************'
    print *,'ADIOS significantly slows down small or medium-size runs, which is the case here, please consider turning it off'
    print *,'**************'
    print *,'**************'
    print *
  endif

  ! checks flags when perfect sphere is set
  if (ASSUME_PERFECT_SPHERE) then
    if (ELLIPTICITY) then
      stop 'ELLIPTICITY not supported when ASSUME_PERFECT_SPHERE is set .true. in constants.h, please check...'
    endif
    if (TOPOGRAPHY) then
      stop 'TOPOGRAPHY not supported when ASSUME_PERFECT_SPHERE is set .true. in constants.h, please check...'
    endif
  endif

  ! see if we are running several independent runs in parallel
  ! if so, add the right directory for that run (group numbers start at zero, but directory names start at run0001, thus we add one)
  ! a negative value for "mygroup" is a convention that indicates that groups (i.e. sub-communicators, one per run) are off
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    LOCAL_PATH = path_to_add(1:len_trim(path_to_add))//LOCAL_PATH(1:len_trim(LOCAL_PATH))
    LOCAL_TMP_PATH = path_to_add(1:len_trim(path_to_add))//LOCAL_TMP_PATH(1:len_trim(LOCAL_TMP_PATH))
  endif

  ! turns on/off corresponding 1-D/3-D model flags
  ! and sets radius for each discontinuity and ocean density values
  call get_model_parameters()

  ! sets time step size and number of layers
  ! right distribution is determined based upon maximum value of NEX
  call get_timestep_and_layers()

  ! time steps: this is an initial estimate based on the record length.
  !             we will need to add additional time steps for reaching the start time at -t0,
  !             which is only known when reading in the CMT source(s).
  !             (see also routine setup_timesteps() in setup_sources_receivers.f90)
  !
  ! initial guess : compute total number of time steps, rounded to next multiple of 100
  if (RECORD_LENGTH_IN_MINUTES < TINYVAL) then
    ! zero length, uses a minimum of 5 steps for testing
    NSTEP = 5
  else
    NSTEP = 100 * (int(RECORD_LENGTH_IN_MINUTES * 60.d0 / (100.d0*DT)) + 1)
  endif

  ! steady state time step
  if (STEADY_STATE_KERNEL) then
    NSTEP_STEADY_STATE = nint(STEADY_STATE_LENGTH_IN_MINUTES * 60.d0 / DT)

    ! checks length
    if (NSTEP_STEADY_STATE == 0) then
      print *, '*****************************************************************'
      print *, 'Warning: STEADY_STATE_KERNEL disabled because STEADY_STATE_LENGTH_IN_MINUTES is zero'
      print *, '*****************************************************************'
      STEADY_STATE_KERNEL = .false.    ! not used any further, but doesn't hurt to reset flag to .false. ...
    endif
  else
    NSTEP_STEADY_STATE = 0
  endif

  ! noise simulations
  if (NOISE_TOMOGRAPHY /= 0) then
    ! time steps needs to be doubled, due to +/- branches (symmetric around zero)
    NSTEP = 2 * NSTEP - 1
  endif

  ! if doing benchmark runs to measure scaling of the code for a limited number of time steps only
  if (DO_BENCHMARK_RUN_ONLY) NSTEP = NSTEP_FOR_BENCHMARK

  ! overrides NSTEP in case specified in Par_file
  if (USER_NSTEP > 0) then
    ! overrides NSTEP
    if (myrank == 0) then
      print *,'simulation number of time steps:'
      print *,'  NSTEP determined = ',NSTEP
      print *,'  Par_file: user overrides with specified NSTEP = ',USER_NSTEP
      print *
    endif
    NSTEP = USER_NSTEP
  endif

  ! debug
  !print *,'initial time steps = ',NSTEP,' record length = ',RECORD_LENGTH_IN_MINUTES,' DT = ',DT

  ! movies: converts values to radians
  MOVIE_EAST = MOVIE_EAST_DEG * DEGREES_TO_RADIANS
  MOVIE_WEST = MOVIE_WEST_DEG * DEGREES_TO_RADIANS
  MOVIE_NORTH = (90.0d0 - MOVIE_NORTH_DEG) * DEGREES_TO_RADIANS ! converting from latitude to colatitude
  MOVIE_SOUTH = (90.0d0 - MOVIE_SOUTH_DEG) * DEGREES_TO_RADIANS
  ! converts movie top/bottom depths to radii
  MOVIE_TOP = (R_PLANET_KM-MOVIE_TOP_KM)/R_PLANET_KM
  MOVIE_BOTTOM = (R_PLANET_KM-MOVIE_BOTTOM_KM)/R_PLANET_KM

  ! half-time duration
  !
  ! computes a default hdur_movie that creates nice looking movies.
  ! Sets HDUR_MOVIE as the minimum period the mesh can resolve
  if (HDUR_MOVIE <= TINYVAL) then
    ! for an estimate based on NGLL == 5, assuming that the number of points per wavelength
    ! coincides with the number of GLL points and thus the element size is the same length a the minimum wavelength:
    !
    !   Earth: 2 * PI * 6371km / 4 / 256 / 2.3 km/s ~ 17 s
    !
    !   Mars : 2 * PI * 3390km / 4 / 256 / 4.0 km/s ~ 5 s
    !
    !   Moon : 2 * PI * 1737.1km / 4 /256 / 1.8 km/s ~ 6 s
    !
    ! adding a second for smoother wavefields
    select case(PLANET_TYPE)
    case (IPLANET_EARTH)
      T_min_res = 17.0 + 1.0
    case (IPLANET_MARS)
      T_min_res = 5.0 + 1.0
    case (IPLANET_MOON)
      T_min_res = 6.0 + 1.0
    case default
      stop 'Invalid planet, type for HDUR_MOVIE estimation not recognized yet'
    end select
    ! minimum hdur for movies
    HDUR_MOVIE = 1.2d0*max(240.d0/NEX_XI * T_min_res * ANGULAR_WIDTH_XI_IN_DEGREES/90.d0, &
                           240.d0/NEX_ETA * T_min_res * ANGULAR_WIDTH_ETA_IN_DEGREES/90.d0)
  endif
  ! noise simulations require MOVIE_SURFACE flag to output wavefield at Earth's surface;
  ! however they don't need to convolve the source time function with any HDUR_MOVIE
  ! since they employ a separate noise-spectrum source S_squared
  if (NOISE_TOMOGRAPHY /= 0) HDUR_MOVIE = 0.d0

  ! checks simulation setup parameters
  call rcp_check_parameters()

  end subroutine rcp_set_compute_parameters

!
!-------------------------------------------------------------------------------------------------
!

  subroutine rcp_set_mesh_parameters()

  use constants
  use shared_parameters

  implicit none

  ! local parameters
  ! doubling layers
  integer :: ielem,elem_doubling_mantle,elem_doubling_middle_outer_core,elem_doubling_bottom_outer_core
  double precision :: DEPTH_SECOND_DOUBLING_REAL,DEPTH_THIRD_DOUBLING_REAL, &
                      DEPTH_FOURTH_DOUBLING_REAL,distance,distance_min,zval
  ! layers
  integer :: NUMBER_OF_MESH_LAYERS,layer_offset
  ! for the cut doublingbrick improvement
  integer :: last_doubling_layer

  ! check that mesh can be coarsened in depth three or four times
  CUT_SUPERBRICK_XI = .false.
  CUT_SUPERBRICK_ETA = .false.

  if (SUPPRESS_CRUSTAL_MESH .and. .not. ADD_4TH_DOUBLING) then
    if (mod(NEX_XI,8) /= 0) stop 'NEX_XI must be a multiple of 8'
    if (mod(NEX_ETA,8) /= 0) stop 'NEX_ETA must be a multiple of 8'
    if (mod(NEX_XI/4,NPROC_XI) /= 0) stop 'NEX_XI must be a multiple of 4*NPROC_XI'
    if (mod(NEX_ETA/4,NPROC_ETA) /= 0) stop 'NEX_ETA must be a multiple of 4*NPROC_ETA'
    if (mod(NEX_XI/8,NPROC_XI) /= 0) CUT_SUPERBRICK_XI = .true.
    if (mod(NEX_ETA/8,NPROC_ETA) /= 0) CUT_SUPERBRICK_ETA = .true.
  else if (SUPPRESS_CRUSTAL_MESH .or. .not. ADD_4TH_DOUBLING) then
    if (mod(NEX_XI,16) /= 0) stop 'NEX_XI must be a multiple of 16'
    if (mod(NEX_ETA,16) /= 0) stop 'NEX_ETA must be a multiple of 16'
    if (mod(NEX_XI/8,NPROC_XI) /= 0) stop 'NEX_XI must be a multiple of 8*NPROC_XI'
    if (mod(NEX_ETA/8,NPROC_ETA) /= 0) stop 'NEX_ETA must be a multiple of 8*NPROC_ETA'
    if (mod(NEX_XI/16,NPROC_XI) /= 0) CUT_SUPERBRICK_XI = .true.
    if (mod(NEX_ETA/16,NPROC_ETA) /= 0) CUT_SUPERBRICK_ETA = .true.
  else
    if (mod(NEX_XI,32) /= 0) stop 'NEX_XI must be a multiple of 32'
    if (mod(NEX_ETA,32) /= 0) stop 'NEX_ETA must be a multiple of 32'
    if (mod(NEX_XI/16,NPROC_XI) /= 0) stop 'NEX_XI must be a multiple of 16*NPROC_XI'
    if (mod(NEX_ETA/16,NPROC_ETA) /= 0) stop 'NEX_ETA must be a multiple of 16*NPROC_ETA'
    if (mod(NEX_XI/32,NPROC_XI) /= 0) CUT_SUPERBRICK_XI = .true.
    if (mod(NEX_ETA/32,NPROC_ETA) /= 0) CUT_SUPERBRICK_ETA = .true.
  endif

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

  !  definition of general mesh parameters
  call define_all_layers(NUMBER_OF_MESH_LAYERS,layer_offset,last_doubling_layer, &
                         ielem,elem_doubling_mantle,elem_doubling_middle_outer_core, &
                         elem_doubling_bottom_outer_core, &
                         DEPTH_SECOND_DOUBLING_REAL,DEPTH_THIRD_DOUBLING_REAL, &
                         DEPTH_FOURTH_DOUBLING_REAL,distance,distance_min,zval, &
                         rmins,rmaxs)

  ! calculates number of elements (NSPEC_REGIONS)
  call count_elements(NEX_XI,NEX_ETA,NEX_PER_PROC_XI,NPROC, &
                        NEX_PER_PROC_ETA,ratio_divide_central_cube, &
                        NSPEC_REGIONS, &
                        NSPEC2D_XI,NSPEC2D_ETA, &
                        NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                        NSPEC1D_RADIAL, &
                        NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
                        NUMBER_OF_MESH_LAYERS,layer_offset, &
                        CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,INCLUDE_CENTRAL_CUBE, &
                        last_doubling_layer)

  ! calculates number of points (NGLOB_REGIONS)
  call count_points(NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube, &
                        NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
                        NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
                        NGLOB_REGIONS, &
                        INCLUDE_CENTRAL_CUBE,NER_TOP_CENTRAL_CUBE_ICB,NEX_XI, &
                        NUMBER_OF_MESH_LAYERS, layer_offset, &
                        CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, &
                        last_doubling_layer)

  if (ATTENUATION) then
    ! to save a huge amount of memory, when 3D attenuation is off it is sufficient to save a single point
    ! per spectral element because the Q attenuation factor is then constant per layer of the geological model
    if (ATTENUATION_3D .or. ATTENUATION_1D_WITH_3D_STORAGE) then
      ATT1 = NGLLX
      ATT2 = NGLLY
      ATT3 = NGLLZ
    else
      ATT1 = 1
      ATT2 = 1
      ATT3 = 1
    endif
    ATT4 = NSPEC_REGIONS(IREGION_CRUST_MANTLE)  ! only used for header file in save_header_file.F90
    ATT5 = NSPEC_REGIONS(IREGION_INNER_CORE)    ! only used for header file
  else
     ATT1 = 1
     ATT2 = 1
     ATT3 = 1
     ATT4 = 1
     ATT5 = 1
  endif

  ! full gravity support
  if (FULL_GRAVITY) call rcp_SIEM_set_mesh_parameters()

  end subroutine rcp_set_mesh_parameters

!
!-------------------------------------------------------------------------------------------------
!

  subroutine rcp_check_parameters()

  use constants
  use shared_parameters

  implicit none

  ! local parameter
  integer :: nex_minimum

! checks parameters

  if (SIMULATION_TYPE /= 1 .and. SIMULATION_TYPE /= 2 .and. SIMULATION_TYPE /= 3) &
    stop 'SIMULATION_TYPE must be either 1, 2 or 3'

  if (NOISE_TOMOGRAPHY < 0 .or. NOISE_TOMOGRAPHY > 3) &
    stop 'NOISE_TOMOGRAPHY must be either 0, 1, 2 or 3'

  if (NCHUNKS /= 1 .and. NCHUNKS /= 2 .and. NCHUNKS /= 3 .and. NCHUNKS /= 6) &
    stop 'NCHUNKS must be either 1, 2, 3 or 6'

  ! this MUST be 90 degrees for two chunks or more to match geometrically
  if (NCHUNKS > 1 .and. abs(ANGULAR_WIDTH_XI_IN_DEGREES - 90.d0) > 0.00000001d0) &
    stop 'ANGULAR_WIDTH_XI_IN_DEGREES must be 90 for more than one chunk'

  ! this can be any value in the case of two chunks
  if (NCHUNKS > 2 .and. abs(ANGULAR_WIDTH_ETA_IN_DEGREES - 90.d0) > 0.00000001d0) &
    stop 'ANGULAR_WIDTH_ETA_IN_DEGREES must be 90 for more than two chunks'

  if (NUMBER_OF_RUNS < 1) &
    stop 'NUMBER_OF_RUNS must be at least 1'

  if (NUMBER_OF_THIS_RUN > NUMBER_OF_RUNS) &
    stop 'NUMBER_OF_THIS_RUN cannot be larger than NUMBER_OF_RUNS'

  if (SIMULATION_TYPE /= 1 .and. NUMBER_OF_RUNS /= 1) &
    stop 'Only 1 run for SIMULATION_TYPE = 2/3'

  if (ABSORBING_CONDITIONS .and. NCHUNKS == 6) &
    stop 'cannot have absorbing conditions in the full Earth'

  if (ABSORBING_CONDITIONS .and. NCHUNKS == 3) &
    stop 'absorbing conditions not supported for three chunks yet'

  if (ABSORB_USING_GLOBAL_SPONGE .and. NCHUNKS /= 6) &
    stop 'Please set NCHUNKS to 6 in Par_file to use ABSORB_USING_GLOBAL_SPONGE'

  if (SAVE_TRANSVERSE_KL_ONLY .and. .not. ANISOTROPIC_KL) &
    stop 'Please set ANISOTROPIC_KL to .true. in Par_file to use SAVE_TRANSVERSE_KL_ONLY'

  if (SAVE_AZIMUTHAL_ANISO_KL_ONLY .and. .not. ANISOTROPIC_KL) &
    stop 'Please set ANISOTROPIC_KL to .true. in Par_file to use SAVE_AZIMUTHAL_ANISO_KL_ONLY'

  if (SAVE_TRANSVERSE_KL_ONLY .and. SAVE_AZIMUTHAL_ANISO_KL_ONLY) &
    stop 'Please set either SAVE_TRANSVERSE_KL_ONLY or SAVE_AZIMUTHAL_ANISO_KL_ONLY to .true., keep the other one .false.'

  if (PARTIAL_PHYS_DISPERSION_ONLY .and. UNDO_ATTENUATION) &
    stop 'cannot have both PARTIAL_PHYS_DISPERSION_ONLY and UNDO_ATTENUATION, they are mutually exclusive'

  ! simulations with undoing attenuation
  if (UNDO_ATTENUATION .and. MOVIE_VOLUME .and. MOVIE_VOLUME_TYPE == 4 ) &
    stop 'MOVIE_VOLUME_TYPE == 4 is not implemented for UNDO_ATTENUATION in order to save memory'

  ! this should not be difficult to fix and test, but not done yet by lack of time
  if (UNDO_ATTENUATION .and. NUMBER_OF_RUNS /= 1) &
    stop 'NUMBER_OF_RUNS should be == 1 for now when using UNDO_ATTENUATION'

  ! this should not be difficult to fix and test, but not done yet by lack of time
  if (UNDO_ATTENUATION .and. NUMBER_OF_THIS_RUN > 1) &
    stop 'we currently do not support NUMBER_OF_THIS_RUN > 1 in the case of UNDO_ATTENUATION'

  if (STEADY_STATE_KERNEL) then
    if (.not. UNDO_ATTENUATION) &
      stop 'STEADY_STATE_KERNEL currently works only when UNDO_ATTENUATION is enabled'
    if (STEADY_STATE_LENGTH_IN_MINUTES > RECORD_LENGTH_IN_MINUTES) &
      stop 'STEADY_STATE_LENGTH_IN_MINUTES cannot be greater than RECORD_LENGTH_IN_MINUTES'
  endif

  if (USE_LDDRK .and. NUMBER_OF_RUNS > 1) &
    stop 'NUMBER_OF_RUNS should be == 1 for now when using USE_LDDRK'

  ! check that reals are either 4 or 8 bytes
  if (CUSTOM_REAL /= SIZE_REAL .and. CUSTOM_REAL /= SIZE_DOUBLE) &
    stop 'wrong size of CUSTOM_REAL for reals'

  ! check that the parameter file is correct
  if (NGNOD /= 27) &
    stop 'number of control nodes must be 27'
  if (NGNOD == 27 .and. NGNOD2D /= 9) &
    stop 'elements with 27 points should have NGNOD2D = 9'

  ! for the number of standard linear solids for attenuation
  if (N_SLS /= 3) &
    stop 'number of SLS must be 3'

  ! check number of slices in each direction
  if (NCHUNKS < 1) &
    stop 'must have at least one chunk'
  if (NPROC_XI < 1) &
    stop 'NPROC_XI must be at least 1'
  if (NPROC_ETA < 1) &
    stop 'NPROC_ETA must be at least 1'

  ! check number of chunks
  if (NCHUNKS /= 1 .and. NCHUNKS /= 2 .and. NCHUNKS /= 3 .and. NCHUNKS /= 6) &
    stop 'only one, two, three or six chunks can be meshed'

  ! check that the central cube can be included
  if (INCLUDE_CENTRAL_CUBE .and. NCHUNKS /= 6) &
    stop 'need six chunks to include central cube'

  ! check that sphere can be cut into slices without getting negative Jacobian
  if (NCHUNKS == 6) then
    ! sets minimum NEX allowed for simulation
    if (TOPOGRAPHY) then
      ! mesh with topography leads to negative Jacobian for NEX < 48, regardless if 1D or 3D model
      nex_minimum = 48
    else
      ! for flat topography, NEX = 32 setting still okay
      nex_minimum = 32
    endif
    ! checks nex
    if (NEX_XI < nex_minimum) &
      stop 'NEX_XI must be greater to cut the sphere into slices with positive Jacobian'
    if (NEX_ETA < nex_minimum) &
      stop 'NEX_ETA must be greater to cut the sphere into slices with positive Jacobian'
  endif

  ! check that topology is correct if more than two chunks
  if (NCHUNKS > 2 .and. NEX_XI /= NEX_ETA) &
    stop 'must have NEX_XI = NEX_ETA for more than two chunks'

  if (NCHUNKS > 2 .and. NPROC_XI /= NPROC_ETA) &
    stop 'must have NPROC_XI = NPROC_ETA for more than two chunks'

  ! small meshes useful for testing, also for GPU version
  if (NCHUNKS > 1 .and. (NPROC_XI == 1 .or. NPROC_ETA == 1)) then
    if (NUMFACES_SHARED < 4 ) &
      stop 'NPROC_XI,NPROC_ETA == 1: please set in constants.h NUMFACES_SHARED and NUMCORNERS_SHARED equal to 4 and recompile'
    if (NUMCORNERS_SHARED < 4 ) &
      stop 'NPROC_XI,NPROC_ETA == 1: please set in constants.h NUMFACES_SHARED and NUMCORNERS_SHARED equal to 4 and recompile'
  endif

  ! checks movie setup
  if (MOVIE_SURFACE) then
    if (MOVIE_COARSE .and. NGLLX /= NGLLY) &
      stop 'MOVIE_COARSE together with MOVIE_SURFACE requires NGLLX == NGLLY'
  endif
  if (MOVIE_VOLUME) then
    if (MOVIE_VOLUME_TYPE < 1 .or. MOVIE_VOLUME_TYPE > 9) &
      stop 'MOVIE_VOLUME_TYPE has to be in range from 1 to 9'
  endif

  ! noise simulations
  if (NOISE_TOMOGRAPHY /= 0) then
    if ((NOISE_TOMOGRAPHY == 1 .or. NOISE_TOMOGRAPHY == 2) .and. SIMULATION_TYPE /= 1) &
      stop 'Noise simulations with NOISE_TOMOGRAPHY == 1 / 2 must have SIMULATION_TYPE == 1'
    if (NOISE_TOMOGRAPHY == 3 .and. SIMULATION_TYPE /= 3) &
      stop 'Noise simulations with NOISE_TOMOGRAPHY == 3 must have SIMULATION_TYPE == 3'
    if (NUMBER_OF_RUNS /= 1 .or. NUMBER_OF_THIS_RUN /= 1) &
      stop 'NUMBER_OF_RUNS and NUMBER_OF_THIS_RUN must be 1 for NOISE TOMOGRAPHY simulation'
    if (ROTATE_SEISMOGRAMS_RT) &
      stop 'Do NOT rotate seismograms in the code, change ROTATE_SEISMOGRAMS_RT in Par_file for noise simulation'
    if (SAVE_ALL_SEISMOS_IN_ONE_FILE .or. USE_BINARY_FOR_LARGE_FILE) &
      stop 'Please set SAVE_ALL_SEISMOS_IN_ONE_FILE and USE_BINARY_FOR_LARGE_FILE to be .false. for noise simulation'
  endif

  ! gravity
  ! makes sure to turn off full gravity flag if no gravity simulation selected
  if (.not. GRAVITY) FULL_GRAVITY = .false.

  ! gravity integrals
  ! in the case of GRAVITY_INTEGRALS we should always use double precision
  if (GRAVITY_INTEGRALS .and. CUSTOM_REAL /= SIZE_DOUBLE) &
    stop 'For GRAVITY_INTEGRALS use double precision i.e. configure the code with --enable-double-precision'

  ! adjoint simulations: seismogram output only works if each process writes out its local seismos
  if (WRITE_SEISMOGRAMS_BY_MAIN .and. SIMULATION_TYPE == 2) &
    stop 'For SIMULATION_TYPE == 2, please set WRITE_SEISMOGRAMS_BY_MAIN to .false.'

  if (NTSTEP_BETWEEN_OUTPUT_SAMPLE < 1) &
    stop 'Invalid NTSTEP_BETWEEN_OUTPUT_SAMPLE, must be >= 1'

!----------------------------------------------
!
! status of implementation
!
!----------------------------------------------
!
! please remove these security checks only after validating new features

  ! July 2013: temporary, the time for Matthieu Lefebvre to merge his ADIOS implementation
  if (ADIOS_ENABLED .and. SAVE_REGULAR_KL ) &
    stop 'ADIOS_ENABLED support not implemented yet for SAVE_REGULAR_KL'

  ! LDDRK
  if (USE_LDDRK .and. (ABSORBING_CONDITIONS .and. .not. UNDO_ATTENUATION) ) &
    stop 'USE_LDDRK support requires to use UNDO_ATTENUATION when absorbing boundaries are turned on'

  if (UNDO_ATTENUATION .and. MOVIE_VOLUME .and. MOVIE_VOLUME_TYPE == 4 ) &
    stop 'UNDO_ATTENUATION support not implemented yet for MOVIE_VOLUME_TYPE == 4 simulations'
  if (UNDO_ATTENUATION .and. SIMULATION_TYPE == 3 .and. (MOVIE_VOLUME .or. MOVIE_SURFACE) ) &
    stop 'UNDO_ATTENUATION support not implemented yet for SIMULATION_TYPE == 3 and movie simulations'

  end subroutine rcp_check_parameters

!
!-------------------------------------------------------------------------------------------------
!

  subroutine rcp_SIEM_set_mesh_parameters()

! sets mesh parameters for spectral-infinite element mesh regions

  use constants
  use shared_parameters

  implicit none

  ! local parameters
  integer :: iregion0,iter_region,nlayer,nspec1layer,nspec1d,nnode1d

  ! safety check
  if (.not. FULL_GRAVITY) return

  ! start region
  iregion0 = IREGION_CRUST_MANTLE

  ! loops over transition-to-infinite and infinite regions
  do iter_region = IREGION_TRINFINITE,IREGION_INFINITE
    ! checks
    if (iter_region == IREGION_TRINFINITE .and. .not. ADD_TRINF) cycle

    ! sets region layers
    if (iter_region == IREGION_TRINFINITE) then
      ! transition-to-infinite region
      nlayer = NLAYER_TRINF
    else
      ! infinite region
      nlayer = 1      ! single layer
    endif

    nspec1layer = NSPEC2D_TOP(iregion0)
    nspec1d = sqrt(real(nspec1layer))

    ! checks if nspec2d top is squared
    if (nspec1d*nspec1d /= nspec1layer) then
      print *,'Error: full gravity 2d surface elements is invalid for region ',iter_region
      print *,'  nspec1D**2 = ',nspec1d*nspec1d,' should be ',nspec1layer
      print *,'Please make sure to have NEX_XI == NEX_ETA.'
      stop 'Invalid number of full gravity 2d surface elements'
    endif

    ! total number of elements
    NSPEC_REGIONS(iter_region) = nlayer * nspec1layer

    ! boundary flags
    NSPEC2D_BOTTOM(iter_region) = nspec1layer
    NSPEC2D_TOP(iter_region) = nspec1layer

    nnode1d = nspec1d * (NGLLX-1)+1 ! ngllx = nglly

    ! total number of global nodes
    NGLOB_REGIONS(iter_region) = (nlayer*(NGLLZ-1)+1) * (nnode1d*nnode1d)

    ! MPI cut-planes
    NSPEC2D_XI(iter_region) = nlayer * nspec1d
    NSPEC2D_ETA(iter_region) = nlayer * nspec1d

    NSPEC2DMAX_XMIN_XMAX(iter_region) = nlayer * nspec1d
    NSPEC2DMAX_YMIN_YMAX(iter_region) = nlayer * nspec1d

    NGLOB2DMAX_XMIN_XMAX(iter_region) = (nlayer*(NGLLZ-1)+1)*nnode1d
    NGLOB2DMAX_YMIN_YMAX(iter_region) = (nlayer*(NGLLZ-1)+1)*nnode1d

    NGLOB1D_RADIAL(iter_region) = (nlayer*(NGLLZ-1)+1)
    NSPEC1D_RADIAL(iter_region) = nlayer

    ! start region for next region
    iregion0 = iter_region
  enddo

  ! full gravity kernels
  ! TODO: check if this can be put into the solver and/or if CALC_GRAVITY_KERNELS is needed
  !if (SIMULATION_TYPE == 3 .and. FULL_GRAVITY) CALC_GRAVITY_KERNELS = .true.

  end subroutine rcp_SIEM_set_mesh_parameters
