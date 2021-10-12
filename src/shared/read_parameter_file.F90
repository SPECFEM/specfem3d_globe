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

!!!!!! VERY IMPORTANT
!!!!!! VERY IMPORTANT
!!!!!! VERY IMPORTANT if you add new parameters to DATA/Par_file, remember to also
!!!!!! VERY IMPORTANT broadcast them with MPI_BCAST in src/shared/broadcast_computed_parameters.f90
!!!!!! VERY IMPORTANT otherwise the code will *NOT* work
!!!!!! VERY IMPORTANT
!!!!!! VERY IMPORTANT

  subroutine read_parameter_file()

  use constants
  use shared_input_parameters

  implicit none

! local variables
  integer :: ier

  ! opens the parameter file: DATA/Par_file
  call open_parameter_file(ier)
  if (ier /= 0) stop 'an error occurred while opening the parameter file'

  ! reads in values
  call read_value_integer(SIMULATION_TYPE, 'SIMULATION_TYPE', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: SIMULATION_TYPE'
  call read_value_integer(NOISE_TOMOGRAPHY, 'NOISE_TOMOGRAPHY', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: NOISE_TOMOGRAPHY'
  call read_value_logical(SAVE_FORWARD, 'SAVE_FORWARD', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: SAVE_FORWARD'
  call read_value_integer(NCHUNKS, 'NCHUNKS', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: NCHUNKS'

  if (NCHUNKS /= 6) then
    ! 1-chunk or 2-chunk simulations
    call read_value_double_precision(ANGULAR_WIDTH_XI_IN_DEGREES, 'ANGULAR_WIDTH_XI_IN_DEGREES', ier)
    if (ier /= 0) stop 'an error occurred while reading the parameter file: ANGULAR_WIDTH_XI_IN_DEGREES...'
    call read_value_double_precision(ANGULAR_WIDTH_ETA_IN_DEGREES, 'ANGULAR_WIDTH_ETA_IN_DEGREES', ier)
    if (ier /= 0) stop 'an error occurred while reading the parameter file: ANGULAR_WIDTH_ETA_IN_DEGREES...'
    call read_value_double_precision(CENTER_LATITUDE_IN_DEGREES, 'CENTER_LATITUDE_IN_DEGREES', ier)
    if (ier /= 0) stop 'an error occurred while reading the parameter file: CENTER_LATITUDE_IN_DEGREES...'
    call read_value_double_precision(CENTER_LONGITUDE_IN_DEGREES, 'CENTER_LONGITUDE_IN_DEGREES', ier)
    if (ier /= 0) stop 'an error occurred while reading the parameter file: CENTER_LONGITUDE_IN_DEGREES...'
    call read_value_double_precision(GAMMA_ROTATION_AZIMUTH, 'GAMMA_ROTATION_AZIMUTH', ier)
    if (ier /= 0) stop 'an error occurred while reading the parameter file: GAMMA_ROTATION_AZIMUTH...'
  endif

  ! number of elements at the surface along the two sides of the first chunk
  call read_value_integer(NEX_XI_read, 'NEX_XI', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: NEX_XI'
  call read_value_integer(NEX_ETA_read, 'NEX_ETA', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: NEX_ETA'
  call read_value_integer(NPROC_XI_read, 'NPROC_XI', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: NPROC_XI'
  call read_value_integer(NPROC_ETA_read, 'NPROC_ETA', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: NPROC_ETA'

  ! physical parameters
  call read_value_logical(OCEANS, 'OCEANS', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: OCEANS'
  call read_value_logical(ELLIPTICITY, 'ELLIPTICITY', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: ELLIPTICITIY'
  call read_value_logical(TOPOGRAPHY, 'TOPOGRAPHY', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: TOPOGRAPHY'
  call read_value_logical(GRAVITY, 'GRAVITY', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: GRAVITY'
  call read_value_logical(ROTATION, 'ROTATION', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: ROTATION'
  call read_value_logical(ATTENUATION, 'ATTENUATION', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: ATTENUATION'

  call read_value_logical(ABSORBING_CONDITIONS, 'ABSORBING_CONDITIONS', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: ABSORBING_CONDITIONS'

  call read_value_logical(ABSORB_USING_GLOBAL_SPONGE, 'ABSORB_USING_GLOBAL_SPONGE', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: ABSORB_USING_GLOBAL_SPONGE'

  ! sponge abosbing boundary properties
  if (ABSORB_USING_GLOBAL_SPONGE) then
    call read_value_double_precision(SPONGE_LATITUDE_IN_DEGREES, 'SPONGE_LATITUDE_IN_DEGREES', ier)
    if (ier /= 0) stop 'an error occurred while reading the parameter file: SPONGE_LATITUDE_IN_DEGREES...'
    call read_value_double_precision(SPONGE_LONGITUDE_IN_DEGREES, 'SPONGE_LONGITUDE_IN_DEGREES', ier)
    if (ier /= 0) stop 'an error occurred while reading the parameter file: SPONGE_LONGITUDE_IN_DEGREES...'
    call read_value_double_precision(SPONGE_RADIUS_IN_DEGREES, 'SPONGE_RADIUS_IN_DEGREES', ier)
    if (ier /= 0) stop 'an error occurred while reading the parameter file: SPONGE_RADIUS_IN_DEGREES...'
  endif

  ! define the velocity model
  call read_value_string(MODEL, 'MODEL', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: MODEL'
  call read_value_double_precision(RECORD_LENGTH_IN_MINUTES, 'RECORD_LENGTH_IN_MINUTES', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: RECORD_LENGTH_IN_MINUTES'

  ! attenuation parameters
  call read_value_logical(PARTIAL_PHYS_DISPERSION_ONLY, 'PARTIAL_PHYS_DISPERSION_ONLY', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: PARTIAL_PHYS_DISPERSION_ONLY'
  call read_value_logical(UNDO_ATTENUATION, 'UNDO_ATTENUATION', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: UNDO_ATTENUATION'
  call read_value_double_precision(MEMORY_INSTALLED_PER_CORE_IN_GB, 'MEMORY_INSTALLED_PER_CORE_IN_GB', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: MEMORY_INSTALLED_PER_CORE_IN_GB'
  call read_value_double_precision(PERCENT_OF_MEM_TO_USE_PER_CORE, 'PERCENT_OF_MEM_TO_USE_PER_CORE', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: PERCENT_OF_MEM_TO_USE_PER_CORE'

  ! mass matrix corrections
  call read_value_logical(EXACT_MASS_MATRIX_FOR_ROTATION, 'EXACT_MASS_MATRIX_FOR_ROTATION', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: EXACT_MASS_MATRIX_FOR_ROTATION'

  ! low-memory Runge-Kutta time scheme
  call read_value_logical(USE_LDDRK, 'USE_LDDRK', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: USE_LDDRK'
  call read_value_logical(INCREASE_CFL_FOR_LDDRK, 'INCREASE_CFL_FOR_LDDRK', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: INCREASE_CFL_FOR_LDDRK'
  call read_value_double_precision(RATIO_BY_WHICH_TO_INCREASE_IT, 'RATIO_BY_WHICH_TO_INCREASE_IT', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: RATIO_BY_WHICH_TO_INCREASE_IT'

  ! movie options
  call read_value_logical(MOVIE_SURFACE, 'MOVIE_SURFACE', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: MOVIE_SURFACE'
  call read_value_logical(MOVIE_VOLUME, 'MOVIE_VOLUME', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: MOVIE_VOLUME'
  call read_value_logical(MOVIE_COARSE,'MOVIE_COARSE', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: MOVIE_COARSE'
  call read_value_integer(NTSTEP_BETWEEN_FRAMES, 'NTSTEP_BETWEEN_FRAMES', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: NTSTEP_BETWEEN_FRAMES'
  call read_value_double_precision(HDUR_MOVIE, 'HDUR_MOVIE', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: HDUR_MOVIE'
  call read_value_integer(MOVIE_VOLUME_TYPE, 'MOVIE_VOLUME_TYPE', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: MOVIE_VOLUME_TYPE'
  call read_value_double_precision(MOVIE_TOP_KM, 'MOVIE_TOP_KM', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: MOVIE_TOP_KM'
  call read_value_double_precision(MOVIE_BOTTOM_KM, 'MOVIE_BOTTOM_KM', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: MOVIE_BOTTOM_KM'
  call read_value_double_precision(MOVIE_WEST_DEG, 'MOVIE_WEST_DEG', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: MOVIE_WEST_DEG'
  call read_value_double_precision(MOVIE_EAST_DEG, 'MOVIE_EAST_DEG', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: MOVIE_EAST_DEG'
  call read_value_double_precision(MOVIE_NORTH_DEG, 'MOVIE_NORTH_DEG', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: MOVIE_NORTH_DEG'
  call read_value_double_precision(MOVIE_SOUTH_DEG, 'MOVIE_SOUTH_DEG', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: MOVIE_SOUTH_DEG'
  call read_value_integer(MOVIE_START, 'MOVIE_START', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: MOVIE_START'
  call read_value_integer(MOVIE_STOP, 'MOVIE_STOP', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: MOVIE_STOP'

  ! run checkpointing
  call read_value_logical(SAVE_MESH_FILES, 'SAVE_MESH_FILES', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: SAVE_MESH_FILES'
  call read_value_integer(NUMBER_OF_RUNS, 'NUMBER_OF_RUNS', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: NUMBER_OF_RUNS'
  call read_value_integer(NUMBER_OF_THIS_RUN, 'NUMBER_OF_THIS_RUN', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: NUMBER_OF_THIS_RUN'

  ! data file output directories
  call read_value_string(LOCAL_PATH, 'LOCAL_PATH', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: LOCAL_PATH'
  call read_value_string(LOCAL_TMP_PATH, 'LOCAL_TMP_PATH', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: LOCAL_TMP_PATH'

  ! user output
  call read_value_integer(NTSTEP_BETWEEN_OUTPUT_INFO, 'NTSTEP_BETWEEN_OUTPUT_INFO', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: NTSTEP_BETWEEN_OUTPUT_INFO'
  call read_value_integer(NTSTEP_BETWEEN_OUTPUT_SEISMOS, 'NTSTEP_BETWEEN_OUTPUT_SEISMOS', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: NTSTEP_BETWEEN_OUTPUT_SEISMOS'
  call read_value_integer(NTSTEP_BETWEEN_READ_ADJSRC, 'NTSTEP_BETWEEN_READ_ADJSRC', ier)
  if (ier /= 0) return

  ! point force sourse
  call read_value_logical(USE_FORCE_POINT_SOURCE, 'USE_FORCE_POINT_SOURCE', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: USE_FORCE_POINT_SOURCE'

  ! monochromatic double couple source
  call read_value_logical(USE_MONOCHROMATIC_CMT_SOURCE, 'USE_MONOCHROMATIC_CMT_SOURCE', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: USE_MONOCHROMATIC_CMT_SOURCE'

  ! option to save strain seismograms
  call read_value_logical(SAVE_SEISMOGRAMS_STRAIN, 'SAVE_SEISMOGRAMS_STRAIN', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: SAVE_SEISMOGRAMS_STRAIN'

  call read_value_logical(SAVE_SEISMOGRAMS_IN_ADJOINT_RUN, 'SAVE_SEISMOGRAMS_IN_ADJOINT_RUN', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: SAVE_SEISMOGRAMS_IN_ADJOINT_RUN'

  ! seismogram output
  call read_value_logical(OUTPUT_SEISMOS_ASCII_TEXT, 'OUTPUT_SEISMOS_ASCII_TEXT', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: OUTPUT_SIESMOS_ASCII_TEXT'
  call read_value_logical(OUTPUT_SEISMOS_SAC_ALPHANUM, 'OUTPUT_SEISMOS_SAC_ALPHANUM', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: OUTPUT_SEISMOS_SAC_ALPHANUM'
  call read_value_logical(OUTPUT_SEISMOS_SAC_BINARY, 'OUTPUT_SEISMOS_SAC_BINARY', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: OUTPUT_SEISMOS_SAC_BINARY'
  call read_value_logical(OUTPUT_SEISMOS_ASDF, 'OUTPUT_SEISMOS_ASDF', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: OUTPUT_ASDF'
  call read_value_logical(OUTPUT_SEISMOS_3D_ARRAY, 'OUTPUT_SEISMOS_3D_ARRAY', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: OUTPUT_3D_ARRAY'
  call read_value_logical(ROTATE_SEISMOGRAMS_RT, 'ROTATE_SEISMOGRAMS_RT', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: ROTATE_SEISMOGRAMS_RT'
  call read_value_logical(WRITE_SEISMOGRAMS_BY_MAIN, 'WRITE_SEISMOGRAMS_BY_MAIN', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: WRITE_SEISMOGRAMS_BY_MAIN'
  call read_value_logical(SAVE_ALL_SEISMOS_IN_ONE_FILE, 'SAVE_ALL_SEISMOS_IN_ONE_FILE', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: SAVE_ALL_SEISMOS_IN_ONE_FILE'
  call read_value_logical(USE_BINARY_FOR_LARGE_FILE, 'USE_BINARY_FOR_LARGE_FILE', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: USE_BINARY_FOR_LARGE_FILE'
  call read_value_logical(RECEIVERS_CAN_BE_BURIED, 'RECEIVERS_CAN_BE_BURIED', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: RECEIVERS_CAN_BE_BURIED'
  call read_value_logical(PRINT_SOURCE_TIME_FUNCTION, 'PRINT_SOURCE_TIME_FUNCTION', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: PRINT_SOURCE_TIME_FUNCTION'

  ! adjoint kernels
  call read_value_logical(READ_ADJSRC_ASDF, 'READ_ADJSRC_ASDF', ier)
  if (ier /= 0) stop 'an error occured while reading the parameter file: READ_ADJSRC_ASDF'
  call read_value_logical(ANISOTROPIC_KL, 'ANISOTROPIC_KL', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: ANISOTROPIC_KL'
  call read_value_logical(SAVE_TRANSVERSE_KL_ONLY, 'SAVE_TRANSVERSE_KL_ONLY', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: SAVE_TRANSVERSE_KL_ONLY'
  call read_value_logical(SAVE_AZIMUTHAL_ANISO_KL_ONLY,'SAVE_AZIMUTHAL_ANISO_KL_ONLY', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: SAVE_AZIMUTHAL_ANISO_KL_ONLY'
  call read_value_logical(APPROXIMATE_HESS_KL, 'APPROXIMATE_HESS_KL', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: APPROXIMATE_HESS_KL'
  call read_value_logical(USE_FULL_TISO_MANTLE, 'USE_FULL_TISO_MANTLE', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: USE_FULL_TISO_MANTLE'
  call read_value_logical(SAVE_SOURCE_MASK, 'SAVE_SOURCE_MASK', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: SAVE_SOURCE_MASK'
  call read_value_logical(SAVE_REGULAR_KL, 'SAVE_REGULAR_KL', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: SAVE_REGULAR_KL'

  call read_value_logical(STEADY_STATE_KERNEL, 'STEADY_STATE_KERNEL', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: STEADY_STATE_KERNEL'
  if (STEADY_STATE_KERNEL) then
    call read_value_double_precision(STEADY_STATE_LENGTH_IN_MINUTES, 'STEADY_STATE_LENGTH_IN_MINUTES', ier)
    if (ier /= 0) stop 'an error occurred while reading the parameter file: STEADY_STATE_LENGTH_IN_MINUTES'
  endif

  ! for simultaneous runs from the same batch job
  call read_value_integer(NUMBER_OF_SIMULTANEOUS_RUNS, 'NUMBER_OF_SIMULTANEOUS_RUNS', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter NUMBER_OF_SIMULTANEOUS_RUNS'
  call read_value_logical(BROADCAST_SAME_MESH_AND_MODEL, 'BROADCAST_SAME_MESH_AND_MODEL', ier)
  if (ier /= 0) stop 'Error reading Par_file parameter BROADCAST_SAME_MESH_AND_MODEL'

  ! GPU simulations
  call read_value_logical(GPU_MODE, 'GPU_MODE', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: GPU_MODE'
  if (GPU_MODE) then
     call read_value_integer(GPU_RUNTIME, 'GPU_RUNTIME', ier)
     if (ier /= 0) stop 'an error occurred while reading the parameter file: GPU_RUNTIME'
     call read_value_string(GPU_PLATFORM, 'GPU_PLATFORM', ier)
     if (ier /= 0) stop 'an error occurred while reading the parameter file: GPU_PLATFORM'
     call read_value_string(GPU_DEVICE, 'GPU_DEVICE', ier)
     if (ier /= 0) stop 'an error occurred while reading the parameter file: GPU_DEVICE'
  endif

  ! ADIOS file format output
  call read_value_logical(ADIOS_ENABLED, 'ADIOS_ENABLED', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: ADIOS_ENABLED'
  call read_value_logical(ADIOS_FOR_FORWARD_ARRAYS, 'ADIOS_FOR_FORWARD_ARRAYS', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: ADIOS_FOR_FORWARD_ARRAYS'
  call read_value_logical(ADIOS_FOR_MPI_ARRAYS, 'ADIOS_FOR_MPI_ARRAYS', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: ADIOS_FOR_MPI_ARRAYS'
  call read_value_logical(ADIOS_FOR_ARRAYS_SOLVER, 'ADIOS_FOR_ARRAYS_SOLVER', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: ADIOS_FOR_ARRAYS_SOLVER'
  call read_value_logical(ADIOS_FOR_SOLVER_MESHFILES, 'ADIOS_FOR_SOLVER_MESHFILES', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: ADIOS_FOR_SOLVER_MESHFILES'
  call read_value_logical(ADIOS_FOR_AVS_DX, 'ADIOS_FOR_AVS_DX', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: ADIOS_FOR_AVS_DX'
  call read_value_logical(ADIOS_FOR_KERNELS, 'ADIOS_FOR_KERNELS', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: ADIOS_FOR_KERNELS'
  call read_value_logical(ADIOS_FOR_MODELS, 'ADIOS_FOR_MODELS', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: ADIOS_FOR_MODELS'
  call read_value_logical(ADIOS_FOR_UNDO_ATTENUATION, 'ADIOS_FOR_UNDO_ATTENUATION', ier)
  if (ier /= 0) stop 'an error occurred while reading the parameter file: ADIOS_FOR_UNDO_ATTENUATION'

  ! optional (sort of an eastern egg type feature) parameters
  !
  ! these parameters may or may not be present in the Par_file.
  ! in case they are in the Par_file, use a negative value to indicate that they can be ignored
  ! and the code should use the automatically determined value as by default
  !
  ! DT specified to override automatic DT
  call read_value_double_precision(USER_DT, 'DT', ier); ier = 0

  ! NSTEP specified to override automatic NSTEP
  call read_value_integer(USER_NSTEP, 'NSTEP', ier); ier = 0
  ! no error checking, continue if not available

  ! closes parameter file
  call close_parameter_file()

#if !defined(USE_ADIOS) && !defined(USE_ADIOS2)
  if (ADIOS_ENABLED) then
    print *
    print *,'**************'
    print *,'**************'
    print *,'ADIOS is enabled in parameter file but the code was not compiled with ADIOS'
    print *,'See --with-adios2 and --with-adios configure options.'
    print *,'**************'
    print *,'**************'
    print *
    stop 'an error occurred while reading the parameter file: ADIOS is enabled but code not built with ADIOS'
  endif
#endif

  end subroutine read_parameter_file
