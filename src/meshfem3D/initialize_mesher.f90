!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            August 2013
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

  subroutine initialize_mesher()

  use meshfem3D_par
  use meshfem3D_models_par

  implicit none

  ! local parameters
  ! timing
  double precision, external :: wtime

! sizeprocs returns number of processes started (should be equal to NPROCTOT).
! myrank is the rank of each process, between 0 and NPROCTOT-1.
! as usual in MPI, process 0 is in charge of coordinating everything
! and also takes care of the main output
! do not create anything for the inner core here, will be done in solver
  call world_size(sizeprocs)
  call world_rank(myrank)

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

! open main output file, only written to by process 0
  if(myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) &
    open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_mesher.txt',status='unknown')

! get MPI starting time
  time_start = wtime()

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '****************************'
    write(IMAIN,*) '*** Specfem3D MPI Mesher ***'
    write(IMAIN,*) '****************************'
    write(IMAIN,*)
  endif

  if (myrank==0) then
    ! reads the parameter file and computes additional parameters
    call read_compute_parameters()
!          MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
!          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
!          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
!          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,RMOHO_FICTITIOUS_IN_MESHER, &
!          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
!          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
!          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,DT, &
!          ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
!          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
!          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
!          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE,MOVIE_VOLUME_TYPE, &
!          MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH,MOVIE_START,MOVIE_STOP, &
!          TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
!          ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST, &
!          ROTATION,ISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE,TOPOGRAPHY,OCEANS,MOVIE_SURFACE, &
!          MOVIE_VOLUME,MOVIE_COARSE,ATTENUATION_3D,RECEIVERS_CAN_BE_BURIED, &
!          PRINT_SOURCE_TIME_FUNCTION,SAVE_MESH_FILES, &
!          ATTENUATION,REFERENCE_1D_MODEL,THREE_D_MODEL,ABSORBING_CONDITIONS, &
!          INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE, &
!          LOCAL_PATH,LOCAL_TMP_PATH,MODEL, &
!          SIMULATION_TYPE,SAVE_FORWARD, &
!          NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
!          NSPEC,NSPEC2D_XI,NSPEC2D_ETA,NSPEC2DMAX_XMIN_XMAX, &
!          NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
!          NSPEC1D_RADIAL,NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB, &
!          ratio_sampling_array, ner, doubling_index,r_bottom,r_top,&
!          this_region_has_a_doubling,rmins,rmaxs,CASE_3D, &
!          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
!          ROTATE_SEISMOGRAMS_RT,ratio_divide_central_cube, &
!          HONOR_1D_SPHERICAL_MOHO,CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,&
!          DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA,&
!          WRITE_SEISMOGRAMS_BY_MASTER,SAVE_ALL_SEISMOS_IN_ONE_FILE, &
!          USE_BINARY_FOR_LARGE_FILE,.false.,NOISE_TOMOGRAPHY)
!

!    ! ADIOS_ENABLED: parameter is optional, may not be in the Par_file
!    call read_adios_parameters(ADIOS_ENABLED, ADIOS_FOR_FORWARD_ARRAYS, &
!        ADIOS_FOR_MPI_ARRAYS, ADIOS_FOR_ARRAYS_SOLVER, &
!        ADIOS_FOR_SOLVER_MESHFILES, ADIOS_FOR_AVS_DX)

  endif

  ! distributes parameters from master to all processes
  call broadcast_computed_parameters()
!                myrank,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
!                NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
!                NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
!                NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA, &
!                NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
!                NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
!                NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE, &
!                MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP, &
!                DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
!                CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
!                RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
!                R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
!                MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH, &
!                RMOHO_FICTITIOUS_IN_MESHER, &
!                MOVIE_SURFACE,MOVIE_VOLUME,RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
!                SAVE_MESH_FILES,ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD, &
!                SAVE_ALL_SEISMOS_IN_ONE_FILE,MOVIE_COARSE,OUTPUT_SEISMOS_ASCII_TEXT, &
!                OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
!                ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER,USE_BINARY_FOR_LARGE_FILE, &
!                LOCAL_PATH,LOCAL_TMP_PATH,MODEL, &
!                NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
!                NSPEC,NSPEC2D_XI,NSPEC2D_ETA, &
!                NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
!                NSPEC1D_RADIAL,NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB, &
!                ratio_sampling_array, ner, doubling_index,r_bottom,r_top, &
!                this_region_has_a_doubling,rmins,rmaxs, &
!                ratio_divide_central_cube,CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, &
!                DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA, &
!                REFERENCE_1D_MODEL,THREE_D_MODEL,ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS, &
!                HONOR_1D_SPHERICAL_MOHO,CRUSTAL,ONE_CRUST,CASE_3D,TRANSVERSE_ISOTROPY, &
!                ISOTROPIC_3D_MANTLE,ANISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE, &
!                ATTENUATION,ATTENUATION_3D,ANISOTROPIC_INNER_CORE,NOISE_TOMOGRAPHY)
!
  ! broadcasts optional ADIOS_ENABLED
!  call broadcast_adios_parameters(myrank,ADIOS_ENABLED, &
!                                  ADIOS_FOR_FORWARD_ARRAYS, ADIOS_FOR_MPI_ARRAYS, ADIOS_FOR_ARRAYS_SOLVER, &
!                                  ADIOS_FOR_SOLVER_MESHFILES, ADIOS_FOR_AVS_DX)

  ! check that the code is running with the requested number of processes
  if(sizeprocs /= NPROCTOT) call exit_MPI(myrank,'wrong number of MPI processes')

  ! compute rotation matrix from Euler angles
  ANGULAR_WIDTH_XI_RAD = ANGULAR_WIDTH_XI_IN_DEGREES * DEGREES_TO_RADIANS
  ANGULAR_WIDTH_ETA_RAD = ANGULAR_WIDTH_ETA_IN_DEGREES * DEGREES_TO_RADIANS

  if(NCHUNKS /= 6) call euler_angles(rotation_matrix,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)

  if (ADIOS_ENABLED) then
    call adios_setup()
  endif

  end subroutine initialize_mesher
