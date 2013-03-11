!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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


  subroutine initialize_simulation()

  use specfem_par
  use specfem_par_movie

  implicit none

  ! local parameters
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC_computed,NGLOB_computed, &
    NSPEC2D_XI,NSPEC2D_ETA,NSPEC1D_RADIAL,NGLOB1D_RADIAL
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
    NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX
  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA
  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA
  integer :: ratio_divide_central_cube
  integer :: sizeprocs
  integer :: ios
  integer :: NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NCHUNKS,NPROC_XI,NPROC_ETA
  double precision :: RMOHO_FICTITIOUS_IN_MESHER,R120,R_CENTRAL_CUBE,CENTER_LONGITUDE_IN_DEGREES,&
    CENTER_LATITUDE_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,&
    GAMMA_ROTATION_AZIMUTH
  integer :: REFERENCE_1D_MODEL,THREE_D_MODEL
  logical :: TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,OCEANS, &
    ATTENUATION,ATTENUATION_NEW,ATTENUATION_3D,ROTATION,ELLIPTICITY, &
    GRAVITY,CASE_3D,ISOTROPIC_3D_MANTLE, &
    HETEROGEN_3D_MANTLE,CRUSTAL,INFLATE_CENTRAL_CUBE
  character(len=150) :: dummystring
  integer, external :: err_occurred

  ! sizeprocs returns number of processes started (should be equal to NPROCTOT).
  ! myrank is the rank of each process, between 0 and sizeprocs-1.
  ! as usual in MPI, process 0 is in charge of coordinating everything
  ! and also takes care of the main output
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (myrank == 0) then

    ! read the parameter file and compute additional parameters
    call read_compute_parameters(MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
         NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
         NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
         NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,RMOHO_FICTITIOUS_IN_MESHER, &
         NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
         NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
         NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,DT, &
         ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
         CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
         RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
         R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE,MOVIE_VOLUME_TYPE, &
         MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST, &
         MOVIE_NORTH,MOVIE_SOUTH,MOVIE_START,MOVIE_STOP, &
         TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
         ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST, &
         ROTATION,ISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE,TOPOGRAPHY,OCEANS,MOVIE_SURFACE, &
         MOVIE_VOLUME,MOVIE_COARSE,ATTENUATION_3D,RECEIVERS_CAN_BE_BURIED, &
         PRINT_SOURCE_TIME_FUNCTION,SAVE_MESH_FILES, &
         ATTENUATION,ATTENUATION_NEW,REFERENCE_1D_MODEL,THREE_D_MODEL,ABSORBING_CONDITIONS, &
         INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE, &
         LOCAL_PATH,LOCAL_TMP_PATH,MODEL, &
         SIMULATION_TYPE,SAVE_FORWARD, &
         NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
         NSPEC_computed,NSPEC2D_XI,NSPEC2D_ETA,NSPEC2DMAX_XMIN_XMAX, &
         NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
         NSPEC1D_RADIAL,NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB_computed, &
         ratio_sampling_array, ner, doubling_index,r_bottom,r_top, &
         this_region_has_a_doubling,rmins,rmaxs,CASE_3D, &
         OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
         ROTATE_SEISMOGRAMS_RT,ratio_divide_central_cube, &
         HONOR_1D_SPHERICAL_MOHO,CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,&
         DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA,&
         WRITE_SEISMOGRAMS_BY_MASTER,SAVE_ALL_SEISMOS_IN_ONE_FILE, &
         USE_BINARY_FOR_LARGE_FILE,.false.,NOISE_TOMOGRAPHY)

    if(err_occurred() /= 0) then
      call exit_MPI(myrank,'an error occurred while reading the parameter file')
    endif

    ! GPU_MODE: parameter is optional, may not be in the Par_file
    call read_gpu_mode(GPU_MODE)
    ! ADIOS_ENABLED: parameter is optional, may not be in the Par_file
    call read_adios_enabled(ADIOS_ENABLED)
  endif

  ! distributes parameters from master to all processes
  ! note: uses NSPEC_computed,NGLOB_computed as arguments
  call broadcast_compute_parameters(myrank,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
                NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
                NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
                NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA, &
                NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
                NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
                NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE, &
                MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP, &
                DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
                CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
                RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
                R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
                MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH, &
                RMOHO_FICTITIOUS_IN_MESHER, &
                MOVIE_SURFACE,MOVIE_VOLUME,RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
                SAVE_MESH_FILES,ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD, &
                SAVE_ALL_SEISMOS_IN_ONE_FILE,MOVIE_COARSE,OUTPUT_SEISMOS_ASCII_TEXT, &
                OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
                ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER,USE_BINARY_FOR_LARGE_FILE, &
                LOCAL_PATH,LOCAL_TMP_PATH,MODEL, &
                NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                NSPEC_computed,NSPEC2D_XI,NSPEC2D_ETA, &
                NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                NSPEC1D_RADIAL,NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB_computed, &
                ratio_sampling_array, ner, doubling_index,r_bottom,r_top, &
                this_region_has_a_doubling,rmins,rmaxs, &
                ratio_divide_central_cube,CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, &
                DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA, &
                REFERENCE_1D_MODEL,THREE_D_MODEL,ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS, &
                HONOR_1D_SPHERICAL_MOHO,CRUSTAL,ONE_CRUST,CASE_3D,TRANSVERSE_ISOTROPY, &
                ISOTROPIC_3D_MANTLE,ANISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE, &
                ATTENUATION,ATTENUATION_NEW,ATTENUATION_3D,ANISOTROPIC_INNER_CORE,NOISE_TOMOGRAPHY)

  ! broadcasts optional GPU_MODE
  call broadcast_gpu_parameters(myrank,GPU_MODE)
  ! broadcasts optional ADIOS_ENABLED 
  call broadcast_adios_parameters(myrank,ADIOS_ENABLED)

  ! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

  ! open main output file, only written to by process 0
  if(myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) &
    open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_solver.txt',status='unknown',action='write')

  if(myrank == 0) then

    write(IMAIN,*)
    write(IMAIN,*) '******************************'
    write(IMAIN,*) '**** Specfem3D MPI Solver ****'
    write(IMAIN,*) '******************************'
    write(IMAIN,*)
    write(IMAIN,*)

    if(FIX_UNDERFLOW_PROBLEM) write(IMAIN,*) 'Fixing slow underflow trapping problem using small initial field'

    write(IMAIN,*)
    write(IMAIN,*) 'There are ',sizeprocs,' MPI processes'
    write(IMAIN,*) 'Processes are numbered from 0 to ',sizeprocs-1
    write(IMAIN,*)

    write(IMAIN,*) 'There are ',NEX_XI,' elements along xi in each chunk'
    write(IMAIN,*) 'There are ',NEX_ETA,' elements along eta in each chunk'
    write(IMAIN,*)
    write(IMAIN,*) 'There are ',NPROC_XI,' slices along xi in each chunk'
    write(IMAIN,*) 'There are ',NPROC_ETA,' slices along eta in each chunk'
    write(IMAIN,*) 'There is a total of ',NPROC,' slices in each chunk'
    write(IMAIN,*) 'There are ',NCHUNKS,' chunks'
    write(IMAIN,*) 'There is a total of ',NPROCTOT,' slices in all the chunks'

    write(IMAIN,*)
    write(IMAIN,*) 'NDIM = ',NDIM
    write(IMAIN,*)
    write(IMAIN,*) 'NGLLX = ',NGLLX
    write(IMAIN,*) 'NGLLY = ',NGLLY
    write(IMAIN,*) 'NGLLZ = ',NGLLZ
    write(IMAIN,*)

    ! write information about precision used for floating-point operations
    if(CUSTOM_REAL == SIZE_REAL) then
      write(IMAIN,*) 'using single precision for the calculations'
    else
      write(IMAIN,*) 'using double precision for the calculations'
    endif
    write(IMAIN,*)
    write(IMAIN,*) 'smallest and largest possible floating-point numbers are: ', &
      tiny(1._CUSTOM_REAL),huge(1._CUSTOM_REAL)
    write(IMAIN,*)

    ! model user parameters
    write(IMAIN,*) 'model: ',trim(MODEL)
    if(OCEANS) then
      write(IMAIN,*) '  incorporating the oceans using equivalent load'
    else
      write(IMAIN,*) '  no oceans'
    endif
    if(ELLIPTICITY) then
      write(IMAIN,*) '  incorporating ellipticity'
    else
      write(IMAIN,*) '  no ellipticity'
    endif
    if(TOPOGRAPHY) then
      write(IMAIN,*) '  incorporating surface topography'
    else
      write(IMAIN,*) '  no surface topography'
    endif
    if(GRAVITY) then
      write(IMAIN,*) '  incorporating self-gravitation (Cowling approximation)'
    else
      write(IMAIN,*) '  no self-gravitation'
    endif
    if(ROTATION) then
      write(IMAIN,*) '  incorporating rotation'
    else
      write(IMAIN,*) '  no rotation'
    endif
    if(ATTENUATION) then
      write(IMAIN,*) '  incorporating attenuation using ',N_SLS,' standard linear solids'
      if(ATTENUATION_3D) write(IMAIN,*)'  using 3D attenuation model'
    else
      write(IMAIN,*) '  no attenuation'
    endif
    write(IMAIN,*)

    ! model mesh parameters
    if(ISOTROPIC_3D_MANTLE) then
      write(IMAIN,*) '  incorporating 3-D lateral variations'
    else
      write(IMAIN,*) '  no 3-D lateral variations'
    endif
    if(HETEROGEN_3D_MANTLE) then
      write(IMAIN,*) '  incorporating heterogeneities in the mantle'
    else
      write(IMAIN,*) '  no heterogeneities in the mantle'
    endif
    if(CRUSTAL) then
      write(IMAIN,*) '  incorporating crustal variations'
    else
      write(IMAIN,*) '  no crustal variations'
    endif
    if(ONE_CRUST) then
      write(IMAIN,*) '  using one layer only in PREM crust'
    else
      write(IMAIN,*) '  using unmodified 1D crustal model with two layers'
    endif
    if(TRANSVERSE_ISOTROPY) then
      write(IMAIN,*) '  incorporating transverse isotropy'
    else
      write(IMAIN,*) '  no transverse isotropy'
    endif
    if(ANISOTROPIC_INNER_CORE) then
      write(IMAIN,*) '  incorporating anisotropic inner core'
    else
      write(IMAIN,*) '  no inner-core anisotropy'
    endif
    if(ANISOTROPIC_3D_MANTLE) then
      write(IMAIN,*) '  incorporating anisotropic mantle'
    else
      write(IMAIN,*) '  no general mantle anisotropy'
    endif
    write(IMAIN,*)
    write(IMAIN,*)

  endif

  ! checks flags
  call initialize_simulation_check(sizeprocs,NPROCTOT,NSPEC_COMPUTED, &
                                  ATTENUATION,ATTENUATION_NEW,ATTENUATION_3D,NCHUNKS,GRAVITY,ROTATION, &
                                  ELLIPTICITY,OCEANS,NPROC_XI,NPROC_ETA, &
                                  TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
                                  ANISOTROPIC_INNER_CORE)

  ! counts receiver stations
  if (SIMULATION_TYPE == 1) then
    rec_filename = 'DATA/STATIONS'
  else
    rec_filename = 'DATA/STATIONS_ADJOINT'
  endif
  call get_value_string(STATIONS, 'solver.STATIONS', rec_filename)

  ! get total number of receivers
  if(myrank == 0) then
    open(unit=IIN,file=STATIONS,iostat=ios,status='old',action='read')
    nrec = 0
    do while(ios == 0)
      read(IIN,"(a)",iostat=ios) dummystring
      if(ios == 0) nrec = nrec + 1
    enddo
    close(IIN)
  endif

  ! broadcast the information read on the master to the nodes
  call bcast_all_singlei(nrec)

  ! checks number of total receivers
  if(nrec < 1) call exit_MPI(myrank,trim(STATIONS)//': need at least one receiver')

  ! initializes GPU cards
  if( GPU_MODE ) call initialize_GPU()

  ! initializes VTK window
  if( VTK_MODE ) then
    if(myrank == 0 ) call initialize_vtkwindow(GPU_MODE)
  endif

  ! save simulation info to ADIOS header
  if (ADIOS_ENABLED) then
    call adios_setup()
  endif
  if (ADIOS_ENABLED) then
    ! TODO use only one ADIOS group to write simulation parameters
    !      i.e. merge write_solver... write_par_... into
    !      write_specfem3D_globe_adios_header()
    !call write_solver_info_header_ADIOS()
    call write_par_file_header_ADIOS()
  endif

  end subroutine initialize_simulation

!
!-------------------------------------------------------------------------------------------------
!

  subroutine initialize_simulation_check(sizeprocs,NPROCTOT,NSPEC_COMPUTED, &
                                        ATTENUATION,ATTENUATION_NEW,ATTENUATION_3D,NCHUNKS,GRAVITY,ROTATION, &
                                        ELLIPTICITY,OCEANS,NPROC_XI,NPROC_ETA, &
                                        TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
                                        ANISOTROPIC_INNER_CORE)

  use specfem_par
  implicit none

  integer :: sizeprocs
  integer :: NPROCTOT,NCHUNKS,NPROC_XI,NPROC_ETA
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC_computed

  logical :: TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,OCEANS, &
    ATTENUATION,ATTENUATION_NEW,ATTENUATION_3D,ROTATION,ELLIPTICITY,GRAVITY


  ! check that the code is running with the requested nb of processes
  if(sizeprocs /= NPROCTOT) call exit_MPI(myrank,'wrong number of MPI processes(initialization specfem)')

  ! check that the code has been compiled with the right values
  if (NSPEC_computed(IREGION_CRUST_MANTLE) /= NSPEC_CRUST_MANTLE) then
      write(IMAIN,*) 'NSPEC_CRUST_MANTLE:',NSPEC_computed(IREGION_CRUST_MANTLE),NSPEC_CRUST_MANTLE
      call exit_MPI(myrank,'error in compiled parameters, please recompile solver 1')
  endif
  if (NSPEC_computed(IREGION_OUTER_CORE) /= NSPEC_OUTER_CORE) then
      write(IMAIN,*) 'NSPEC_OUTER_CORE:',NSPEC_computed(IREGION_OUTER_CORE),NSPEC_OUTER_CORE
      call exit_MPI(myrank,'error in compiled parameters, please recompile solver 2')
  endif
  if (NSPEC_computed(IREGION_INNER_CORE) /= NSPEC_INNER_CORE) then
      write(IMAIN,*) 'NSPEC_INNER_CORE:',NSPEC_computed(IREGION_INNER_CORE),NSPEC_INNER_CORE
      call exit_MPI(myrank,'error in compiled parameters, please recompile solver 3')
  endif
  if (ATTENUATION_3D .NEQV. ATTENUATION_3D_VAL) then
      write(IMAIN,*) 'ATTENUATION_3D:',ATTENUATION_3D,ATTENUATION_3D_VAL
      call exit_MPI(myrank,'error in compiled parameters ATTENUATION_3D, please recompile solver')
  endif
  if (NCHUNKS /= NCHUNKS_VAL) then
      write(IMAIN,*) 'NCHUNKS:',NCHUNKS,NCHUNKS_VAL
      call exit_MPI(myrank,'error in compiled parameters NCHUNKS, please recompile solver')
  endif
  if (GRAVITY .NEQV. GRAVITY_VAL) then
      write(IMAIN,*) 'GRAVITY:',GRAVITY,GRAVITY_VAL
      call exit_MPI(myrank,'error in compiled parameters GRAVITY, please recompile solver')
  endif
  if (ROTATION .NEQV. ROTATION_VAL) then
      write(IMAIN,*) 'ROTATION:',ROTATION,ROTATION_VAL
      call exit_MPI(myrank,'error in compiled parameters ROTATION, please recompile solver')
  endif
  if (ATTENUATION .NEQV. ATTENUATION_VAL) then
      write(IMAIN,*) 'ATTENUATION:',ATTENUATION,ATTENUATION_VAL
      call exit_MPI(myrank,'error in compiled parameters ATTENUATION, please recompile solver')
  endif
  if (ATTENUATION_NEW .NEQV. ATTENUATION_NEW_VAL) then
      write(IMAIN,*) 'ATTENUATION_NEW:',ATTENUATION_NEW,ATTENUATION_NEW_VAL
      call exit_MPI(myrank,'error in compiled parameters ATTENUATION_NEW, please recompile solver')
  endif
  if (ELLIPTICITY .NEQV. ELLIPTICITY_VAL) then
      write(IMAIN,*) 'ELLIPTICITY:',ELLIPTICITY,ELLIPTICITY_VAL
      call exit_MPI(myrank,'error in compiled parameters ELLIPTICITY, please recompile solver')
  endif
  if (OCEANS .NEQV. OCEANS_VAL) then
      write(IMAIN,*) 'OCEANS:',OCEANS,OCEANS_VAL
      call exit_MPI(myrank,'error in compiled parameters OCEANS, please recompile solver')
  endif
  if (NPROC_XI /= NPROC_XI_VAL) then
      write(IMAIN,*) 'NPROC_XI:',NPROC_XI,NPROC_XI_VAL
      call exit_MPI(myrank,'error in compiled parameters NPROC_XI, please recompile solver')
  endif
  if (NPROC_ETA /= NPROC_ETA_VAL) then
      write(IMAIN,*) 'NPROC_ETA:',NPROC_ETA,NPROC_ETA_VAL
      call exit_MPI(myrank,'error in compiled parameters NPROC_ETA, please recompile solver')
  endif
  if (NPROCTOT /= NPROCTOT_VAL) then
      write(IMAIN,*) 'NPROCTOT:',NPROCTOT,NPROCTOT_VAL
      call exit_MPI(myrank,'error in compiled parameters NPROCTOT, please recompile solver')
  endif
  if (NEX_XI /= NEX_XI_VAL) then
      write(IMAIN,*) 'NEX_XI:',NEX_XI,NEX_XI_VAL
      call exit_MPI(myrank,'error in compiled parameters NEX_XI, please recompile solver')
  endif
  if (NEX_ETA /= NEX_ETA_VAL) then
      write(IMAIN,*) 'NEX_ETA:',NEX_ETA,NEX_ETA_VAL
      call exit_MPI(myrank,'error in compiled parameters NEX_ETA, please recompile solver')
  endif
  if (TRANSVERSE_ISOTROPY .NEQV. TRANSVERSE_ISOTROPY_VAL) then
      write(IMAIN,*) 'TRANSVERSE_ISOTROPY:',TRANSVERSE_ISOTROPY,TRANSVERSE_ISOTROPY_VAL
      call exit_MPI(myrank,'error in compiled parameters, please recompile solver 14')
  endif
  if (ANISOTROPIC_3D_MANTLE .NEQV. ANISOTROPIC_3D_MANTLE_VAL) then
      write(IMAIN,*) 'ANISOTROPIC_3D_MANTLE:',ANISOTROPIC_3D_MANTLE,ANISOTROPIC_3D_MANTLE_VAL
      call exit_MPI(myrank,'error in compiled parameters, please recompile solver 15')
  endif
  if (ANISOTROPIC_INNER_CORE .NEQV. ANISOTROPIC_INNER_CORE_VAL) then
      write(IMAIN,*) 'ANISOTROPIC_INNER_CORE:',ANISOTROPIC_INNER_CORE,ANISOTROPIC_INNER_CORE_VAL
      call exit_MPI(myrank,'error in compiled parameters, please recompile solver 16')
  endif

  ! check simulation pararmeters
  if (SIMULATION_TYPE /= 1 .and.  SIMULATION_TYPE /= 2 .and. SIMULATION_TYPE /= 3) &
    call exit_MPI(myrank, 'SIMULATION_TYPE can only be 1, 2, or 3')

  if (SIMULATION_TYPE /= 1 .and. NSOURCES > 999999)  &
    call exit_MPI(myrank, &
    'for adjoint simulations, NSOURCES <= 999999, if you need more change i6.6 in write_seismograms.f90')

  if((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
    if ( ATTENUATION_VAL) then
      ! checks mimic flag:
      ! attenuation for adjoint simulations must have USE_ATTENUATION_MIMIC set by xcreate_header_file

!daniel: att - debug todo check attenuation mimick
      if( USE_ATTENUATION_MIMIC .eqv. .false. ) &
        call exit_MPI(myrank,'error in compiled attenuation parameters, please recompile solver 17b')

      ! user output
      if( myrank == 0 ) then
        write(IMAIN,*) 'incorporates ATTENUATION for time-reversed simulation'
        write(IMAIN,*)
      endif
    endif

    ! checks adjoint array dimensions
    if(NSPEC_CRUST_MANTLE_ADJOINT /= NSPEC_CRUST_MANTLE &
      .or. NSPEC_OUTER_CORE_ADJOINT /= NSPEC_OUTER_CORE &
      .or. NSPEC_INNER_CORE_ADJOINT /= NSPEC_INNER_CORE &
      .or. NGLOB_CRUST_MANTLE_ADJOINT /= NGLOB_CRUST_MANTLE &
      .or. NGLOB_OUTER_CORE_ADJOINT /= NGLOB_OUTER_CORE &
      .or. NGLOB_INNER_CORE_ADJOINT /= NGLOB_INNER_CORE) &
      call exit_MPI(myrank, 'improper dimensions of adjoint arrays, please recompile solver 18')
  endif

  ! checks attenuation
  if( ATTENUATION_VAL ) then
    if (NSPEC_CRUST_MANTLE_ATTENUAT /= NSPEC_CRUST_MANTLE) &
       call exit_MPI(myrank, 'NSPEC_CRUST_MANTLE_ATTENUAT /= NSPEC_CRUST_MANTLE, exit')
    if (NSPEC_INNER_CORE_ATTENUATION /= NSPEC_INNER_CORE) &
       call exit_MPI(myrank, 'NSPEC_INNER_CORE_ATTENUATION /= NSPEC_INNER_CORE, exit')
  endif

  ! checks strain storage
  if (ATTENUATION_VAL .or. SIMULATION_TYPE /= 1 .or. SAVE_FORWARD &
    .or. (MOVIE_VOLUME .and. SIMULATION_TYPE /= 3)) then
    if( COMPUTE_AND_STORE_STRAIN .neqv. .true. ) &
      call exit_MPI(myrank, 'error in compiled compute_and_store_strain parameter, please recompile solver 19')
  else
    if( COMPUTE_AND_STORE_STRAIN .neqv. .false. ) &
      call exit_MPI(myrank, 'error in compiled compute_and_store_strain parameter, please recompile solver 20')
  endif

  if (SIMULATION_TYPE == 3 .and. (ANISOTROPIC_3D_MANTLE_VAL .or. ANISOTROPIC_INNER_CORE_VAL)) &
     call exit_MPI(myrank, 'anisotropic model is not implemented for kernel simulations yet')

  ! checks model for transverse isotropic kernel computation
  if( SAVE_TRANSVERSE_KL ) then
    if( ANISOTROPIC_3D_MANTLE_VAL ) then
        call exit_mpi(myrank,'error SAVE_TRANSVERSE_KL: Earth model not supported yet')
    endif
    if( SIMULATION_TYPE == 3 ) then
      if( .not. ANISOTROPIC_KL ) then
        call exit_mpi(myrank,'error SAVE_TRANSVERSE_KL: needs anisotropic kernel calculations')
      endif
    endif
  endif

  ! check for GPU runs
  if( GPU_MODE ) then
    if( NGLLX /= 5 .or. NGLLY /= 5 .or. NGLLZ /= 5 ) &
      call exit_mpi(myrank,'GPU mode can only be used if NGLLX == NGLLY == NGLLZ == 5')
    if( CUSTOM_REAL /= 4 ) &
      call exit_mpi(myrank,'GPU mode runs only with CUSTOM_REAL == 4')
    if( ATTENUATION_VAL ) then
      if( N_SLS /= 3 ) &
        call exit_mpi(myrank,'GPU mode does not support N_SLS /= 3 yet')
    endif
  endif

  end subroutine initialize_simulation_check

!
!-------------------------------------------------------------------------------------------------
!

  subroutine initialize_GPU()

! initialization for GPU cards

  use specfem_par
  implicit none
  ! local parameters
  integer :: ncuda_devices,ncuda_devices_min,ncuda_devices_max

  ! GPU_MODE now defined in Par_file
  if(myrank == 0 ) then
    write(IMAIN,*)
    write(IMAIN,*) "GPU_MODE Active."
  endif

  ! check for GPU runs
  if( NGLLX /= 5 .or. NGLLY /= 5 .or. NGLLZ /= 5 ) &
    stop 'GPU mode can only be used if NGLLX == NGLLY == NGLLZ == 5'
  if( CUSTOM_REAL /= 4 ) &
    stop 'GPU mode runs only with CUSTOM_REAL == 4'
  if( ATTENUATION_VAL ) then
    if( N_SLS /= 3 ) &
      stop 'GPU mode does not support N_SLS /= 3 yet'
  endif

  ! initializes GPU and outputs info to files for all processes
  call initialize_cuda_device(myrank,ncuda_devices)

  ! collects min/max of local devices found for statistics
  call sync_all()
  call min_all_i(ncuda_devices,ncuda_devices_min)
  call max_all_i(ncuda_devices,ncuda_devices_max)

  if( myrank == 0 ) then
    write(IMAIN,*) "GPU number of devices per node: min =",ncuda_devices_min
    write(IMAIN,*) "                                max =",ncuda_devices_max
    write(IMAIN,*)
  endif

  end subroutine initialize_GPU
