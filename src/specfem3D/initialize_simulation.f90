!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
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

  subroutine initialize_simulation(myrank,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
                NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
                NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
                NER_TOP_CENTRAL_CUBE_ICB,ANGULAR_WIDTH_XI_IN_DEGREES,NEX_XI,NEX_ETA, &
                NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
                NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
                NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,SIMULATION_TYPE, &
                DT,ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R220,R400,R600,R670,R771,&
                RTOPDDOUBLEPRIME,RCMB,RICB, &
                RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS, &
                MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP, &
                HDUR_MOVIE,MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST, &
                MOVIE_NORTH,MOVIE_SOUTH,MOVIE_SURFACE,MOVIE_VOLUME, &
                RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
                SAVE_MESH_FILES,ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,SAVE_FORWARD, &
                SAVE_ALL_SEISMOS_IN_ONE_FILE,MOVIE_COARSE,OUTPUT_SEISMOS_ASCII_TEXT, &
                OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
                ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER,USE_BINARY_FOR_LARGE_FILE, &
                LOCAL_PATH,MODEL,OUTPUT_FILES, &
                ratio_sampling_array, ner, doubling_index,r_bottom,r_top, &
                this_region_has_a_doubling,rmins,rmaxs, &
                TOPOGRAPHY,HONOR_1D_SPHERICAL_MOHO,ONE_CRUST, &
                nspl,rspl,espl,espl2,ibathy_topo, &
                NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
                xigll,yigll,zigll,wxgll,wygll,wzgll,wgll_cube, &
                hprime_xx,hprime_yy,hprime_zz,hprime_xxT, &
                hprimewgll_xx,hprimewgll_yy,hprimewgll_zz,hprimewgll_xxT, &
                wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                rec_filename,STATIONS,nrec,NOISE_TOMOGRAPHY,SAVE_REGULAR_KL, &
                PARTIAL_PHYS_DISPERSION_ONLY,UNDO_ATTENUATION,NT_DUMP_ATTENUATION, &
          USE_LDDRK,INCREASE_CFL_FOR_LDDRK,ANISOTROPIC_KL,SAVE_TRANSVERSE_KL_ONLY,APPROXIMATE_HESS_KL, &
          USE_FULL_TISO_MANTLE,SAVE_SOURCE_MASK,GPU_MODE,ADIOS_ENABLED,ADIOS_FOR_FORWARD_ARRAYS, &
          ADIOS_FOR_MPI_ARRAYS,ADIOS_FOR_ARRAYS_SOLVER,ADIOS_FOR_AVS_DX,RATIO_BY_WHICH_TO_INCREASE_IT, &
          ATT1,ATT2,ATT3,ATT4,ATT5,EXACT_MASS_MATRIX_FOR_ROTATION)

  use mpi

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank

  ! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA, &
          NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,SIMULATION_TYPE, &
          MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP,NOISE_TOMOGRAPHY,NT_DUMP_ATTENUATION, &
          ATT1,ATT2,ATT3,ATT4,ATT5

  double precision DT,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
          MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH

  logical MOVIE_SURFACE,MOVIE_VOLUME,RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,SAVE_FORWARD, &
          SAVE_ALL_SEISMOS_IN_ONE_FILE,MOVIE_COARSE,OUTPUT_SEISMOS_ASCII_TEXT,&
          OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY,&
          ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER,USE_BINARY_FOR_LARGE_FILE, &
          SAVE_REGULAR_KL,PARTIAL_PHYS_DISPERSION_ONLY,UNDO_ATTENUATION

  character(len=150) LOCAL_PATH,OUTPUT_FILES

  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ratio_sampling_array,ner
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index

  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top

  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs


  ! mesh model parameters
  logical TOPOGRAPHY,HONOR_1D_SPHERICAL_MOHO,ONE_CRUST

  ! for ellipticity
  integer nspl
  double precision rspl(NR),espl(NR),espl2(NR)

  ! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision, dimension(NGLLZ) :: zigll,wzgll
  ! product of weights for gravity term
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube
  ! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xxT,hprimewgll_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  character(len=150) rec_filename,STATIONS
  integer nrec

  ! local parameters
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC_computed,NGLOB_computed, &
    NSPEC2D_XI,NSPEC2D_ETA,NSPEC1D_RADIAL

  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA

  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA

  integer :: ratio_divide_central_cube
  integer :: sizeprocs
  integer :: ier,i,j
  integer :: ios
  integer :: NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NCHUNKS,NPROC_XI,NPROC_ETA

  double precision :: RMOHO_FICTITIOUS_IN_MESHER,R120,R_CENTRAL_CUBE,CENTER_LONGITUDE_IN_DEGREES,&
    CENTER_LATITUDE_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,ANGULAR_WIDTH_XI_IN_DEGREES,&
    GAMMA_ROTATION_AZIMUTH, &
          RATIO_BY_WHICH_TO_INCREASE_IT

  integer :: REFERENCE_1D_MODEL,THREE_D_MODEL

  logical :: TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,OCEANS, &
    ATTENUATION,ATTENUATION_3D,ROTATION,ELLIPTICITY, &
    GRAVITY,CASE_3D,ISOTROPIC_3D_MANTLE, &
    HETEROGEN_3D_MANTLE,CRUSTAL,INFLATE_CENTRAL_CUBE, &
          USE_LDDRK,INCREASE_CFL_FOR_LDDRK,ANISOTROPIC_KL,SAVE_TRANSVERSE_KL_ONLY,APPROXIMATE_HESS_KL, &
          USE_FULL_TISO_MANTLE,SAVE_SOURCE_MASK,GPU_MODE,ADIOS_ENABLED,ADIOS_FOR_FORWARD_ARRAYS, &
          ADIOS_FOR_MPI_ARRAYS,ADIOS_FOR_ARRAYS_SOLVER,ADIOS_FOR_AVS_DX, &
          EXACT_MASS_MATRIX_FOR_ROTATION,ATTENUATION_1D_WITH_3D_STORAGE

  character(len=150) :: MODEL,dummystring

  ! sizeprocs returns number of processes started (should be equal to NPROCTOT).
  ! myrank is the rank of each process, between 0 and sizeprocs-1.
  ! as usual in MPI, process 0 is in charge of coordinating everything
  ! and also takes care of the main output
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

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
         ATTENUATION,REFERENCE_1D_MODEL,THREE_D_MODEL,ABSORBING_CONDITIONS, &
         INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE, &
         LOCAL_PATH,MODEL, &
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
         USE_BINARY_FOR_LARGE_FILE,.false.,NOISE_TOMOGRAPHY,SAVE_REGULAR_KL, &
         PARTIAL_PHYS_DISPERSION_ONLY,UNDO_ATTENUATION,NT_DUMP_ATTENUATION, &
          USE_LDDRK,INCREASE_CFL_FOR_LDDRK,ANISOTROPIC_KL,SAVE_TRANSVERSE_KL_ONLY,APPROXIMATE_HESS_KL, &
          USE_FULL_TISO_MANTLE,SAVE_SOURCE_MASK,GPU_MODE,ADIOS_ENABLED,ADIOS_FOR_FORWARD_ARRAYS, &
          ADIOS_FOR_MPI_ARRAYS,ADIOS_FOR_ARRAYS_SOLVER,ADIOS_FOR_AVS_DX,RATIO_BY_WHICH_TO_INCREASE_IT, &
          ATT1,ATT2,ATT3,ATT4,ATT5,EXACT_MASS_MATRIX_FOR_ROTATION,ATTENUATION_1D_WITH_3D_STORAGE)

  endif

  ! distributes parameters from master to all processes
  ! note: uses NSPEC_computed,NGLOB_computed as arguments
  call broadcast_computed_parameters(myrank,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
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
                LOCAL_PATH,MODEL, &
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
                ATTENUATION,ATTENUATION_3D,ANISOTROPIC_INNER_CORE,NOISE_TOMOGRAPHY,SAVE_REGULAR_KL, &
                PARTIAL_PHYS_DISPERSION_ONLY,UNDO_ATTENUATION,NT_DUMP_ATTENUATION, &
          USE_LDDRK,INCREASE_CFL_FOR_LDDRK,ANISOTROPIC_KL,SAVE_TRANSVERSE_KL_ONLY,APPROXIMATE_HESS_KL, &
          USE_FULL_TISO_MANTLE,SAVE_SOURCE_MASK,GPU_MODE,ADIOS_ENABLED,ADIOS_FOR_FORWARD_ARRAYS, &
          ADIOS_FOR_MPI_ARRAYS,ADIOS_FOR_ARRAYS_SOLVER,ADIOS_FOR_AVS_DX,RATIO_BY_WHICH_TO_INCREASE_IT, &
          ATT1,ATT2,ATT3,ATT4,ATT5,EXACT_MASS_MATRIX_FOR_ROTATION,ATTENUATION_1D_WITH_3D_STORAGE)

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
    if(ANISOTROPIC_INNER_CORE_VAL) then
      write(IMAIN,*) '  incorporating anisotropic inner core'
    else
      write(IMAIN,*) '  no inner-core anisotropy'
    endif
    if(ANISOTROPIC_3D_MANTLE_VAL) then
      write(IMAIN,*) '  incorporating anisotropic mantle'
    else
      write(IMAIN,*) '  no general mantle anisotropy'
    endif

    write(IMAIN,*)
    write(IMAIN,*)

  endif

  ! check that the code is running with the requested nb of processes
  if(sizeprocs /= NPROCTOT) call exit_MPI(myrank,'wrong number of MPI processes(initialization specfem)')

  ! check that the code has been compiled with the right values
  if (NSPEC_computed(IREGION_CRUST_MANTLE) /= NSPEC_CRUST_MANTLE) then
      call exit_MPI(myrank,'error in compiled parameters, please recompile solver 1')
  endif
  if (NSPEC_computed(IREGION_OUTER_CORE) /= NSPEC_OUTER_CORE) then
       call exit_MPI(myrank,'error in compiled parameters, please recompile solver 2')
  endif
  if (NSPEC_computed(IREGION_INNER_CORE) /= NSPEC_INNER_CORE) then
       call exit_MPI(myrank,'error in compiled parameters, please recompile solver 3')
  endif
  if (ATTENUATION_3D .NEQV. ATTENUATION_3D_VAL) then
       call exit_MPI(myrank,'error in compiled parameters, please recompile solver 4')
  endif
  if (NCHUNKS /= NCHUNKS_VAL) then
       call exit_MPI(myrank,'error in compiled parameters, please recompile solver 6')
  endif
  if (GRAVITY .NEQV. GRAVITY_VAL) then
       call exit_MPI(myrank,'error in compiled parameters, please recompile solver 7')
  endif
  if (ROTATION .NEQV. ROTATION_VAL) then
       call exit_MPI(myrank,'error in compiled parameters, please recompile solver 8')
  endif
  if (ATTENUATION .NEQV. ATTENUATION_VAL) then
       call exit_MPI(myrank,'error in compiled parameters, please recompile solver 9')
  endif
  if (ELLIPTICITY .NEQV. ELLIPTICITY_VAL) then
       call exit_MPI(myrank,'error in compiled parameters, please recompile solver 10')
  endif
  if (OCEANS .NEQV. OCEANS_VAL) then
      call exit_MPI(myrank,'error in compiled parameters, please recompile solver 10')
  endif
  if (NPROCTOT /= NPROCTOT_VAL) then
       call exit_MPI(myrank,'error in compiled parameters, please recompile solver 11')
  endif
  if (NPROC_XI /= NPROC_XI_VAL) then
       call exit_MPI(myrank,'error in compiled parameters, please recompile solver 11')
  endif
  if (NPROC_ETA /= NPROC_ETA_VAL) then
       call exit_MPI(myrank,'error in compiled parameters, please recompile solver 11')
  endif
  if (NEX_XI /= NEX_XI_VAL) then
       call exit_MPI(myrank,'error in compiled parameters, please recompile solver 12')
  endif
  if (NEX_ETA /= NEX_ETA_VAL) then
       call exit_MPI(myrank,'error in compiled parameters, please recompile solver 13')
  endif
  if (TRANSVERSE_ISOTROPY .NEQV. TRANSVERSE_ISOTROPY_VAL) then
       call exit_MPI(myrank,'error in compiled parameters, please recompile solver 14')
  endif
  if (ANISOTROPIC_3D_MANTLE .NEQV. ANISOTROPIC_3D_MANTLE_VAL) then
       call exit_MPI(myrank,'error in compiled parameters, please recompile solver 15')
  endif
  if (ANISOTROPIC_INNER_CORE .NEQV. ANISOTROPIC_INNER_CORE_VAL) then
       call exit_MPI(myrank,'error in compiled parameters, please recompile solver 16')
  endif

  ! check simulation pararmeters
  if (SIMULATION_TYPE /= 1 .and.  SIMULATION_TYPE /= 2 .and. SIMULATION_TYPE /= 3) &
    call exit_MPI(myrank, 'SIMULATION_TYPE can only be 1, 2, or 3')

  if (SIMULATION_TYPE /= 1 .and. NSOURCES > 999999)  &
    call exit_MPI(myrank, 'for adjoint simulations, NSOURCES <= 999999, if you need more change i6.6 in write_seismograms.f90')

  if(UNDO_ATTENUATION .and. PARTIAL_PHYS_DISPERSION_ONLY) &
          call exit_MPI(myrank,'cannot have both UNDO_ATTENUATION and PARTIAL_PHYS_DISPERSION_ONLY')

  if((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
    if ( ATTENUATION_VAL) then
      ! checks mimic flag:
      ! attenuation for adjoint simulations must have PARTIAL_PHYS_DISPERSION_ONLY set by xcreate_header_file
      if(.not. UNDO_ATTENUATION)then
        if(.not. PARTIAL_PHYS_DISPERSION_ONLY) &
     call exit_MPI(myrank,'ATTENUATION for adjoint runs or SAVE_FORWARD requires UNDO_ATTENUATION or PARTIAL_PHYS_DISPERSION_ONLY')
      endif

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
    if (NSPEC_CRUST_MANTLE_ATTENUATION /= NSPEC_CRUST_MANTLE) &
       call exit_MPI(myrank, 'NSPEC_CRUST_MANTLE_ATTENUATION /= NSPEC_CRUST_MANTLE, exit')
    if (NSPEC_INNER_CORE_ATTENUATION /= NSPEC_INNER_CORE) &
       call exit_MPI(myrank, 'NSPEC_INNER_CORE_ATTENUATION /= NSPEC_INNER_CORE, exit')
  endif

  ! checks strain storage
  if (ATTENUATION_VAL .or. SIMULATION_TYPE /= 1 .or. SAVE_FORWARD &
    .or. (MOVIE_VOLUME .and. SIMULATION_TYPE /= 3)) then
    if( COMPUTE_AND_STORE_STRAIN_VAL .neqv. .true. ) &
      call exit_MPI(myrank, 'error in compiled compute_and_store_strain parameter, please recompile solver 19')
  else
    if( COMPUTE_AND_STORE_STRAIN_VAL .neqv. .false. ) &
      call exit_MPI(myrank, 'error in compiled compute_and_store_strain parameter, please recompile solver 20')
  endif

  if (SIMULATION_TYPE == 3 .and. (ANISOTROPIC_3D_MANTLE_VAL .or. ANISOTROPIC_INNER_CORE_VAL)) &
     call exit_MPI(myrank, 'anisotropic model is not implemented for kernel simulations yet')

  ! checks model for transverse isotropic kernel computation
  if( SAVE_TRANSVERSE_KL_ONLY ) then
    if( ANISOTROPIC_3D_MANTLE_VAL ) then
        call exit_mpi(myrank,'error SAVE_TRANSVERSE_KL_ONLY: Earth model not supported yet')
    endif
    if( SIMULATION_TYPE == 3 ) then
      if( .not. ANISOTROPIC_KL ) then
        call exit_mpi(myrank,'error SAVE_TRANSVERSE_KL_ONLY: needs anisotropic kernel calculations')
      endif
    endif
  endif

  ! make ellipticity
  if(ELLIPTICITY_VAL) call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)

  ! read topography and bathymetry file
  if(myrank == 0 .and. (TOPOGRAPHY .or. OCEANS_VAL)) call read_topo_bathy_file(ibathy_topo)
  ! broadcast the information read on the master to the nodes
  call MPI_BCAST(ibathy_topo,NX_BATHY_VAL*NY_BATHY_VAL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  ! set up GLL points, weights and derivation matrices
  call define_derivation_matrices(xigll,yigll,zigll,wxgll,wygll,wzgll, &
         hprime_xx,hprime_yy,hprime_zz, &
         hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
         wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube)

  if( USE_DEVILLE_PRODUCTS_VAL ) then

  ! check that optimized routines from Deville et al. (2002) can be used
    if(NGLLX /= 5 .or. NGLLY /= 5 .or. NGLLZ /= 5) &
      stop 'Deville et al. (2002) routines can only be used if NGLLX = NGLLY = NGLLZ = 5'

    ! define transpose of derivation matrix
    do j = 1,NGLLY
      do i = 1,NGLLX
        hprime_xxT(j,i) = hprime_xx(i,j)
        hprimewgll_xxT(j,i) = hprimewgll_xx(i,j)
      enddo
    enddo
  endif

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
  call MPI_BCAST(nrec,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  if(nrec < 1) call exit_MPI(myrank,trim(STATIONS)//': need at least one receiver')

  end subroutine initialize_simulation

