!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

! output a table with the best parameters to use for large runs on a given machine
! in order to use the total amount of memory available in the best possible way

  program estimate_best_values_runs

  implicit none

  include "../../../constants.h"

  integer, parameter :: NB_COLUMNS_TABLE = 10

! maximum total number of processors we want to see in the table
! integer, parameter :: MAX_NUMBER_OF_PROCS = 100000
! integer, parameter :: MAX_NUMBER_OF_PROCS =  4000  ! current maximum on pangu at Caltech
! integer, parameter :: MAX_NUMBER_OF_PROCS = 10240  ! current maximum on MareNostrum in Barcelona
! integer, parameter :: MAX_NUMBER_OF_PROCS = 62976  ! current maximum on Ranger in Texas
! integer, parameter :: MAX_NUMBER_OF_PROCS = 212992 ! current maximum on BlueGene at LLNL
  integer, parameter :: MAX_NUMBER_OF_PROCS = 12150  ! current maximum at CINES in France

! minimum total number of processors we want to see in the table
! integer, parameter :: MIN_NUMBER_OF_PROCS =  100
  integer, parameter :: MIN_NUMBER_OF_PROCS = 1000

! amount of memory installed per processor core on the system (in gigabytes)
! double precision, parameter :: MAX_MEMORY_PER_CORE = 1.5d0 ! pangu at Caltech
! double precision, parameter :: MAX_MEMORY_PER_CORE = 2.0d0 ! Ranger in Texas, MareNostrum in Barcelona
  double precision, parameter :: MAX_MEMORY_PER_CORE = 4.0d0 ! CINES in France

! base value depends if we implement three or four doublings (default is three)
  integer, parameter :: NB_DOUBLING = 3
  integer :: BASE_VALUE

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE, &
          REFERENCE_1D_MODEL,THREE_D_MODEL,MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP

  double precision DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
          MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH,RMOHO_FICTITIOUS_IN_MESHER

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE, &
          TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_COARSE,ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ATTENUATION,CASE_3D, &
          ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
          ROTATE_SEISMOGRAMS_RT,HONOR_1D_SPHERICAL_MOHO,WRITE_SEISMOGRAMS_BY_MASTER,&
          SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,first_time,HETEROGEN_3D_MANTLE

  character(len=150) LOCAL_PATH,MODEL

! parameters deduced from parameters read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

! this for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC,NSPEC2D_XI,NSPEC2D_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
               NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,nglob

  double precision :: static_memory_size
  integer :: NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
         NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUAT, &
         NSPEC_INNER_CORE_ATTENUATION, &
         NSPEC_CRUST_MANTLE_STR_OR_ATT,NSPEC_INNER_CORE_STR_OR_ATT, &
         NSPEC_CRUST_MANTLE_STR_AND_ATT,NSPEC_INNER_CORE_STR_AND_ATT, &
         NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_INNER_CORE_STRAIN_ONLY, &
         NSPEC_CRUST_MANTLE_ADJOINT, &
         NSPEC_OUTER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
         NGLOB_CRUST_MANTLE_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
         NGLOB_INNER_CORE_ADJOINT,NSPEC_OUTER_CORE_ROT_ADJOINT, &
         NSPEC_CRUST_MANTLE_STACEY,NSPEC_OUTER_CORE_STACEY, &
         NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION

  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,EMULATE_ONLY,already_printed
  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA

  integer :: c,counter,offset

  double precision :: mem_per_core,percent

! ************** PROGRAM STARTS HERE **************

! count the total number of sources in the CMTSOLUTION file
  NSOURCES = 1

  BASE_VALUE = 2**NB_DOUBLING

  first_time = .true.

  print *
  print *,'Number of GB of memory per core installed in the machine: ',MAX_MEMORY_PER_CORE
  print *

! make sure we do not start below the lower limit
  if(6 * int(sqrt(MIN_NUMBER_OF_PROCS / 6.d0))**2 == MIN_NUMBER_OF_PROCS) then
    offset = 0
  else
    offset = 1
  endif

! the total number of processors is 6 * NPROC_XI^2
  do NPROC_XI = int(sqrt(MIN_NUMBER_OF_PROCS / 6.d0)) + offset,int(sqrt(MAX_NUMBER_OF_PROCS / 6.d0))

    counter = 1
    c = 0
    NEX_XI = 0
    already_printed = .false.

! loop on the columns of the table
    do while (counter <= NB_COLUMNS_TABLE)

      c = c + 1
      NEX_XI = BASE_VALUE * c * NPROC_XI

! impose an upper limit on current "reasonably large" values to avoid an infinite loop
! these values should be okay for the next few years (i.e., on machines around 1 petaflops)
! (this written in August 2008)
      if(c > 15 .or. NEX_XI > 6000) exit

! also exclude very coarse meshes below NEX_XI = 64
      if(NEX_XI >= 64 .and. mod(NEX_XI,2*BASE_VALUE) == 0) then

        counter = counter + 1

! read the parameter file and compute additional parameters
    EMULATE_ONLY = .true.
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
          MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH,MOVIE_START,MOVIE_STOP, &
          TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
          ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST, &
          ROTATION,ISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE,TOPOGRAPHY,OCEANS,MOVIE_SURFACE, &
          MOVIE_VOLUME,MOVIE_COARSE,ATTENUATION_3D,RECEIVERS_CAN_BE_BURIED, &
          PRINT_SOURCE_TIME_FUNCTION,SAVE_MESH_FILES, &
          ATTENUATION,REFERENCE_1D_MODEL,THREE_D_MODEL,ABSORBING_CONDITIONS, &
          INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,LOCAL_PATH,MODEL,SIMULATION_TYPE,SAVE_FORWARD, &
          NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
          NSPEC,NSPEC2D_XI,NSPEC2D_ETA, &
          NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
          NSPEC1D_RADIAL,NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB, &
          ratio_sampling_array, ner, doubling_index,r_bottom,r_top,this_region_has_a_doubling,rmins,rmaxs,CASE_3D, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
          ROTATE_SEISMOGRAMS_RT,ratio_divide_central_cube,HONOR_1D_SPHERICAL_MOHO,CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,&
          DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA,&
          WRITE_SEISMOGRAMS_BY_MASTER,SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,EMULATE_ONLY)

  if(first_time) then
    first_time = .false.
    if(ATTENUATION) then
      print *,'ATTENUATION = .true.'
    else
      print *,'ATTENUATION = .false.'
    endif
  endif

!! compute memory usage per processor core
! evaluate the amount of static memory needed by the solver
  call memory_eval(OCEANS,ABSORBING_CONDITIONS,ATTENUATION,ANISOTROPIC_3D_MANTLE,&
                   TRANSVERSE_ISOTROPY,ANISOTROPIC_INNER_CORE,ROTATION,&
                   ONE_CRUST,doubling_index,this_region_has_a_doubling,&
                   ner,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_sampling_array,&
                   NSPEC,nglob,SIMULATION_TYPE,MOVIE_VOLUME,SAVE_FORWARD, &
         NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
         NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUAT, &
         NSPEC_INNER_CORE_ATTENUATION, &
         NSPEC_CRUST_MANTLE_STR_OR_ATT,NSPEC_INNER_CORE_STR_OR_ATT, &
         NSPEC_CRUST_MANTLE_STR_AND_ATT,NSPEC_INNER_CORE_STR_AND_ATT, &
         NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_INNER_CORE_STRAIN_ONLY, &
         NSPEC_CRUST_MANTLE_ADJOINT, &
         NSPEC_OUTER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
         NGLOB_CRUST_MANTLE_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
         NGLOB_INNER_CORE_ADJOINT,NSPEC_OUTER_CORE_ROT_ADJOINT, &
         NSPEC_CRUST_MANTLE_STACEY,NSPEC_OUTER_CORE_STACEY, &
         NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION,static_memory_size)

  mem_per_core = static_memory_size/1073741824.d0

  percent = 100.d0 * mem_per_core / MAX_MEMORY_PER_CORE

  if(percent > 100.d0) goto 777

  if(percent < 0.d0) then
    if(.not. already_printed) then
      print *
      print *,'NPROC, % of the machine, NPROC_XI, NEX_XI, GB used/core, % mem used/core:'
    endif
    already_printed = .true.
    write(*,"(' ',i5,'  ',f6.2,'% ',i3,'  ',i4,'  ',f6.2,'GB  ',f6.2,'% **mesher fails/bug**')") &
      6*NPROC_XI**2,100.d0*6*NPROC_XI**2/dble(MAX_NUMBER_OF_PROCS),NPROC_XI,NEX_XI,mem_per_core,percent
    goto 777

  else if(percent >= 93.d0) then
    if(.not. already_printed) then
      print *
      print *,'NPROC, % of the machine, NPROC_XI, NEX_XI, GB used/core, % mem used/core:'
    endif
    already_printed = .true.
    write(*,"(' ',i5,'  ',f6.2,'% ',i3,'  ',i4,'  ',f6.2,'GB  ',f6.2,'% **too high**')") &
      6*NPROC_XI**2,100.d0*6*NPROC_XI**2/dble(MAX_NUMBER_OF_PROCS),NPROC_XI,NEX_XI,mem_per_core,percent
    goto 777

  else if(percent >= 85.d0) then
    if(.not. already_printed) then
      print *
      print *,'NPROC, % of the machine, NPROC_XI, NEX_XI, GB used/core, % mem used/core:'
    endif
    already_printed = .true.
    write(*,"(' ',i5,'  ',f6.2,'% ',i3,'  ',i4,'  ',f6.2,'GB  ',f6.2,'% **excellent**')") &
      6*NPROC_XI**2,100.d0*6*NPROC_XI**2/dble(MAX_NUMBER_OF_PROCS),NPROC_XI,NEX_XI,mem_per_core,percent
    goto 777

  else if(percent >= 80.d0) then
    if(.not. already_printed) then
      print *
      print *,'NPROC, % of the machine, NPROC_XI, NEX_XI, GB used/core, % mem used/core:'
    endif
    already_printed = .true.
    write(*,"(' ',i5,'  ',f6.2,'% ',i3,'  ',i4,'  ',f6.2,'GB  ',f6.2,'% **very good**')") &
      6*NPROC_XI**2,100.d0*6*NPROC_XI**2/dble(MAX_NUMBER_OF_PROCS),NPROC_XI,NEX_XI,mem_per_core,percent
    goto 777

  else if(percent >= 70.d0) then
    if(.not. already_printed) then
      print *
      print *,'NPROC, % of the machine, NPROC_XI, NEX_XI, GB used/core, % mem used/core:'
    endif
    already_printed = .true.
    write(*,"(' ',i5,'  ',f6.2,'% ',i3,'  ',i4,'  ',f6.2,'GB  ',f6.2,'% **not bad**')") &
      6*NPROC_XI**2,100.d0*6*NPROC_XI**2/dble(MAX_NUMBER_OF_PROCS),NPROC_XI,NEX_XI,mem_per_core,percent
    goto 777

  else if(percent >= 60.d0) then
    if(.not. already_printed) then
      print *
      print *,'NPROC, % of the machine, NPROC_XI, NEX_XI, GB used/core, % mem used/core:'
    endif
    already_printed = .true.
    write(*,"(' ',i5,'  ',f6.2,'% ',i3,'  ',i4,'  ',f6.2,'GB  ',f6.2,'% **possible**')") &
      6*NPROC_XI**2,100.d0*6*NPROC_XI**2/dble(MAX_NUMBER_OF_PROCS),NPROC_XI,NEX_XI,mem_per_core,percent

! uncomment line below to completely suppress cases that give an occupancy below 50%
  else ! if(percent >= 50.d0) then
    if(.not. already_printed) then
      print *
      print *,'NPROC, % of the machine, NPROC_XI, NEX_XI, GB used/core, % mem used/core:'
    endif
    already_printed = .true.
    write(*,"(' ',i5,'  ',f6.2,'% ',i3,'  ',i4,'  ',f6.2,'GB  ',f6.2,'%')") &
      6*NPROC_XI**2,100.d0*6*NPROC_XI**2/dble(MAX_NUMBER_OF_PROCS),NPROC_XI,NEX_XI,mem_per_core,percent
  endif

      endif
    enddo

 777 continue

  enddo

  print *

  end program estimate_best_values_runs

!
!----  include subroutines from the main code
!

  include "../../read_compute_parameters.f90"

  include "../../memory_eval.f90"

  include "../../read_value_parameters.f90"

  include "../../get_value_parameters.f90"

  include "../../auto_ner.f90"

