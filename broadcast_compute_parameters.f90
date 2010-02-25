!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology / Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology / Princeton University and University of Pau / CNRS / INRIA
!                            March 2010
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

  subroutine broadcast_compute_parameters(myrank,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
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
                NSPEC,NSPEC2D_XI,NSPEC2D_ETA, &
                NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                NSPEC1D_RADIAL,NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB, &
                ratio_sampling_array, ner, doubling_index,r_bottom,r_top, &
                this_region_has_a_doubling,rmins,rmaxs, &
                ratio_divide_central_cube,CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, &
                DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA, &
                REFERENCE_1D_MODEL,THREE_D_MODEL,ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS, &  
                HONOR_1D_SPHERICAL_MOHO,CRUSTAL,ONE_CRUST,CASE_3D,TRANSVERSE_ISOTROPY, &
                ISOTROPIC_3D_MANTLE,ANISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE, &
                ATTENUATION,ATTENUATION_3D,ANISOTROPIC_INNER_CORE)
  
  implicit none

! standard include of the MPI library
  include 'mpif.h'
  include "constants.h"
  include "precision.h"
  
  integer myrank

  ! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE, &
          MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP

  double precision DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
          MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH, &
          RMOHO_FICTITIOUS_IN_MESHER

  logical   MOVIE_SURFACE,MOVIE_VOLUME,RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD, &
          SAVE_ALL_SEISMOS_IN_ONE_FILE,MOVIE_COARSE,OUTPUT_SEISMOS_ASCII_TEXT,&
          OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY,&
          ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER,USE_BINARY_FOR_LARGE_FILE
          
  character(len=150) LOCAL_PATH,MODEL

  ! parameters to be computed based upon parameters above read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC,NSPEC2D_XI,NSPEC2D_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NSPEC1D_RADIAL,NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
      NGLOB

  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ratio_sampling_array,ner
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top

  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

  integer ratio_divide_central_cube

  ! for the cut doublingbrick improvement
  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA
  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA

  ! mesh model parameters
  integer REFERENCE_1D_MODEL,THREE_D_MODEL  

  logical ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS, &  
    HONOR_1D_SPHERICAL_MOHO,CRUSTAL,ONE_CRUST,CASE_3D,TRANSVERSE_ISOTROPY, &
    ISOTROPIC_3D_MANTLE,ANISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE, &
    ATTENUATION,ATTENUATION_3D,ANISOTROPIC_INNER_CORE

  ! local parameters
  double precision, dimension(31) :: bcast_double_precision
  integer, dimension(38) :: bcast_integer
  logical, dimension(35) :: bcast_logical
  integer ier
  
  ! master process prepares broadcasting arrays
  if (myrank==0) then
    ! count the total number of sources in the CMTSOLUTION file
    call count_number_of_sources(NSOURCES)

    ! funny way to pass parameters in arrays from master to all other processes
    ! rather than single values one by one to reduce MPI communication calls:
    ! sets up broadcasting array
    bcast_integer = (/MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
            NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
            NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
            NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA, &
            NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
            NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
            NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,&
            SIMULATION_TYPE,REFERENCE_1D_MODEL,THREE_D_MODEL,NPROC,NPROCTOT, &
            NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube,&
            MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP,NSOURCES/)

    bcast_logical = (/TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
            CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE, &
            TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,ATTENUATION_3D, &
            RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
            SAVE_MESH_FILES,ATTENUATION, &
            ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD,CASE_3D,&
            CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,SAVE_ALL_SEISMOS_IN_ONE_FILE, &
            HONOR_1D_SPHERICAL_MOHO,MOVIE_COARSE, &
            OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY,&
            ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER,USE_BINARY_FOR_LARGE_FILE/)

    bcast_double_precision = (/DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
            CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
            RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
            R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
            MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH,&
            RMOHO_FICTITIOUS_IN_MESHER /)
  endif

  ! broadcasts the information read on the master to the nodes
  call MPI_BCAST(bcast_integer,38,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(bcast_double_precision,31,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(bcast_logical,35,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)

  ! broadcasts non-single value parameters
  call MPI_BCAST(LOCAL_PATH,150,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(MODEL,150,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(ner,MAX_NUMBER_OF_MESH_LAYERS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(ratio_sampling_array,MAX_NUMBER_OF_MESH_LAYERS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(doubling_index,MAX_NUMBER_OF_MESH_LAYERS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(r_bottom,MAX_NUMBER_OF_MESH_LAYERS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(r_top,MAX_NUMBER_OF_MESH_LAYERS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(rmins,MAX_NUMBER_OF_MESH_LAYERS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(rmaxs,MAX_NUMBER_OF_MESH_LAYERS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(this_region_has_a_doubling,MAX_NUMBER_OF_MESH_LAYERS,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(NSPEC,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(NSPEC2D_XI,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(NSPEC2D_ETA,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(NSPEC2DMAX_XMIN_XMAX,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(NSPEC2DMAX_YMIN_YMAX,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(NSPEC2D_BOTTOM,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(NSPEC2D_TOP,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(NSPEC1D_RADIAL,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(NGLOB1D_RADIAL,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(NGLOB2DMAX_XMIN_XMAX,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(NGLOB2DMAX_YMIN_YMAX,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(NGLOB,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(DIFF_NSPEC1D_RADIAL,NB_SQUARE_CORNERS*NB_CUT_CASE,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(DIFF_NSPEC2D_ETA,NB_SQUARE_EDGES_ONEDIR*NB_CUT_CASE,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(DIFF_NSPEC2D_XI,NB_SQUARE_EDGES_ONEDIR*NB_CUT_CASE,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  ! non-master processes set their parameters
  if (myrank /=0) then
    
    ! please, be careful with ordering and counting here
    ! integers
    MIN_ATTENUATION_PERIOD = bcast_integer(1)
    MAX_ATTENUATION_PERIOD = bcast_integer(2)
    NER_CRUST = bcast_integer(3)
    NER_80_MOHO = bcast_integer(4)
    NER_220_80 = bcast_integer(5)
    NER_400_220 = bcast_integer(6)
    NER_600_400 = bcast_integer(7)
    NER_670_600 = bcast_integer(8)
    NER_771_670 = bcast_integer(9)
    NER_TOPDDOUBLEPRIME_771 = bcast_integer(10)
    NER_CMB_TOPDDOUBLEPRIME = bcast_integer(11)
    NER_OUTER_CORE = bcast_integer(12)
    NER_TOP_CENTRAL_CUBE_ICB = bcast_integer(13)
    NEX_XI = bcast_integer(14)
    NEX_ETA = bcast_integer(15)
    NPROC_XI = bcast_integer(16)
    NPROC_ETA = bcast_integer(17)
    NTSTEP_BETWEEN_OUTPUT_SEISMOS = bcast_integer(18)
    NTSTEP_BETWEEN_READ_ADJSRC = bcast_integer(19)
    NSTEP = bcast_integer(20)
    NSOURCES = bcast_integer(21)
    NTSTEP_BETWEEN_FRAMES = bcast_integer(22)
    NTSTEP_BETWEEN_OUTPUT_INFO = bcast_integer(23)
    NUMBER_OF_RUNS = bcast_integer(24)
    NUMBER_OF_THIS_RUN = bcast_integer(25)
    NCHUNKS = bcast_integer(26)
    SIMULATION_TYPE = bcast_integer(27)
    REFERENCE_1D_MODEL = bcast_integer(28)
    THREE_D_MODEL = bcast_integer(29)
    NPROC = bcast_integer(30)
    NPROCTOT = bcast_integer(31)
    NEX_PER_PROC_XI = bcast_integer(32)
    NEX_PER_PROC_ETA = bcast_integer(33)
    ratio_divide_central_cube = bcast_integer(34)
    MOVIE_VOLUME_TYPE = bcast_integer(35)
    MOVIE_START = bcast_integer(36)
    MOVIE_STOP = bcast_integer(37)
    NSOURCES = bcast_integer(38)
    
    ! logicals
    TRANSVERSE_ISOTROPY = bcast_logical(1)
    ANISOTROPIC_3D_MANTLE = bcast_logical(2)
    ANISOTROPIC_INNER_CORE = bcast_logical(3)
    CRUSTAL = bcast_logical(4)
    ELLIPTICITY = bcast_logical(5)
    GRAVITY = bcast_logical(6)
    ONE_CRUST = bcast_logical(7)
    ROTATION = bcast_logical(8)
    ISOTROPIC_3D_MANTLE = bcast_logical(9)
    HETEROGEN_3D_MANTLE = bcast_logical(10)
    TOPOGRAPHY = bcast_logical(11)
    OCEANS = bcast_logical(12)
    MOVIE_SURFACE = bcast_logical(13)
    MOVIE_VOLUME = bcast_logical(14)
    ATTENUATION_3D = bcast_logical(15)
    RECEIVERS_CAN_BE_BURIED = bcast_logical(16)
    PRINT_SOURCE_TIME_FUNCTION = bcast_logical(17)
    SAVE_MESH_FILES = bcast_logical(18)
    ATTENUATION = bcast_logical(19)
    ABSORBING_CONDITIONS = bcast_logical(20)
    INCLUDE_CENTRAL_CUBE = bcast_logical(21)
    INFLATE_CENTRAL_CUBE = bcast_logical(22)
    SAVE_FORWARD = bcast_logical(23)
    CASE_3D = bcast_logical(24)
    CUT_SUPERBRICK_XI = bcast_logical(25)
    CUT_SUPERBRICK_ETA = bcast_logical(26)
    SAVE_ALL_SEISMOS_IN_ONE_FILE = bcast_logical(27)
    HONOR_1D_SPHERICAL_MOHO = bcast_logical(28)
    MOVIE_COARSE= bcast_logical(29)
    OUTPUT_SEISMOS_ASCII_TEXT= bcast_logical(30)
    OUTPUT_SEISMOS_SAC_ALPHANUM= bcast_logical(31)
    OUTPUT_SEISMOS_SAC_BINARY= bcast_logical(32)
    ROTATE_SEISMOGRAMS_RT= bcast_logical(33)
    WRITE_SEISMOGRAMS_BY_MASTER= bcast_logical(34)
    USE_BINARY_FOR_LARGE_FILE= bcast_logical(35)

    ! double precisions
    DT = bcast_double_precision(1)
    ANGULAR_WIDTH_XI_IN_DEGREES = bcast_double_precision(2)
    ANGULAR_WIDTH_ETA_IN_DEGREES = bcast_double_precision(3)
    CENTER_LONGITUDE_IN_DEGREES = bcast_double_precision(4)
    CENTER_LATITUDE_IN_DEGREES = bcast_double_precision(5)
    GAMMA_ROTATION_AZIMUTH = bcast_double_precision(6)
    ROCEAN = bcast_double_precision(7)
    RMIDDLE_CRUST = bcast_double_precision(8)
    RMOHO = bcast_double_precision(9)
    R80 = bcast_double_precision(10)
    R120 = bcast_double_precision(11)
    R220 = bcast_double_precision(12)
    R400 = bcast_double_precision(13)
    R600 = bcast_double_precision(14)
    R670 = bcast_double_precision(15)
    R771 = bcast_double_precision(16)
    RTOPDDOUBLEPRIME = bcast_double_precision(17)
    RCMB = bcast_double_precision(18)
    RICB = bcast_double_precision(19)
    R_CENTRAL_CUBE = bcast_double_precision(20)
    RHO_TOP_OC = bcast_double_precision(21)
    RHO_BOTTOM_OC = bcast_double_precision(22)
    RHO_OCEANS = bcast_double_precision(23)
    HDUR_MOVIE = bcast_double_precision(24)
    MOVIE_TOP = bcast_double_precision(25)
    MOVIE_BOTTOM = bcast_double_precision(26)
    MOVIE_WEST = bcast_double_precision(27)
    MOVIE_EAST = bcast_double_precision(28)
    MOVIE_NORTH = bcast_double_precision(29)
    MOVIE_SOUTH = bcast_double_precision(30)
    RMOHO_FICTITIOUS_IN_MESHER = bcast_double_precision(31)

  endif
    
  end subroutine broadcast_compute_parameters
