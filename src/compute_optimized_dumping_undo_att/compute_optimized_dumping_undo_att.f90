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

! compute the optimal interval at which to dump restart files to disk to undo attenuation in an exact way

! Dimitri Komatitsch and Zhinan Xie, CNRS Marseille, France, June 2013.

  program xcompute_optimized_dumping

  implicit none

  include "constants.h"

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE, &
          REFERENCE_1D_MODEL,THREE_D_MODEL,MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP,NOISE_TOMOGRAPHY

  double precision DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
          MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH,RMOHO_FICTITIOUS_IN_MESHER

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE, &
          TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_COARSE,ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ATTENUATION,CASE_3D, &
          ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
          ROTATE_SEISMOGRAMS_RT,HONOR_1D_SPHERICAL_MOHO,WRITE_SEISMOGRAMS_BY_MASTER,&
          SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,SAVE_REGULAR_KL

  character(len=150) LOCAL_PATH,MODEL

! parameters deduced from parameters read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

! this for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC, &
               NSPEC2D_XI, &
               NSPEC2D_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
               NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
               nglob

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

  integer :: iregion
  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA
  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA
  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_CORNERS) :: NGLOB1D_RADIAL_CORNER
  integer, dimension(MAX_NUM_REGIONS) :: NGLOB1D_RADIAL_TEMP

  integer :: NT_DUMP_ATTENUATION_optimal_to_use,number_of_dumpings_to_do
  double precision :: gigabytes_avail_per_core,percentage_to_use_per_core,what_we_can_use_in_GB,size_to_store_at_each_time_step, &
                         disk_size_of_each_dumping

  logical :: PARTIAL_PHYS_DISPERSION_ONLY,UNDO_ATTENUATION
  integer :: NT_DUMP_ATTENUATION

! ************** PROGRAM STARTS HERE **************

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
         WRITE_SEISMOGRAMS_BY_MASTER,SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,.false.,NOISE_TOMOGRAPHY,&
         SAVE_REGULAR_KL,PARTIAL_PHYS_DISPERSION_ONLY,UNDO_ATTENUATION,NT_DUMP_ATTENUATION)

! optimal dumping interval calculation can only be done when SIMULATION_TYPE == 3 in the Par_file,
! thus set it to that value here in this serial code even if it has a different value in the Par_file
  SIMULATION_TYPE = 3

! count the total number of sources in the CMTSOLUTION file
  call count_number_of_sources(NSOURCES)

  do iregion=1,MAX_NUM_REGIONS
    NGLOB1D_RADIAL_CORNER(iregion,:) = NGLOB1D_RADIAL(iregion)
  enddo

  if (CUT_SUPERBRICK_XI .or. CUT_SUPERBRICK_ETA) then
    NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + &
                                                maxval(DIFF_NSPEC1D_RADIAL(:,:))*(NGLLZ-1)
  endif

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

  NGLOB1D_RADIAL_TEMP(:) = &
  (/maxval(NGLOB1D_RADIAL_CORNER(1,:)),maxval(NGLOB1D_RADIAL_CORNER(2,:)),maxval(NGLOB1D_RADIAL_CORNER(3,:))/)

  print *
  print *,'number of processors = ',NPROCTOT
  print *
  print *,'maximum number of points per region = ',nglob(IREGION_CRUST_MANTLE)
  print *
  print *,'total elements per slice = ',sum(NSPEC)
  print *,'total points per slice = ',sum(nglob)
  print *
  print *,'number of time steps = ',NSTEP
  print *

  print *,'approximate static memory needed by the solver:'
  print *,'----------------------------------------------'
  print *
  print *,'(lower bound, usually the real amount used is 5% to 10% higher)'
  print *
  print *,'(you can get a more precise estimate of the size used per MPI process'
  print *,' by typing "size -d bin/xspecfem3D"'
  print *,' after compiling the code with the DATA/Par_file you plan to use)'
  print *
  print *,'size of static arrays per slice = ',static_memory_size/1.d6,' MB'
  print *,'                                = ',static_memory_size/1048576.d0,' MiB'
  print *,'                                = ',static_memory_size/1.d9,' GB'
  print *,'                                = ',static_memory_size/1073741824.d0,' GiB'
  print *

  if(static_memory_size*dble(NPROCTOT)/1.d6 < 10000.d0) then
    print *,'size of static arrays for all slices = ',static_memory_size*dble(NPROCTOT)/1.d6,' MB'
    print *,'                                     = ',static_memory_size*dble(NPROCTOT)/1048576.d0,' MiB'
    print *,'                                     = ',static_memory_size*dble(NPROCTOT)/1.d9,' GB'
  else
    print *,'size of static arrays for all slices = ',static_memory_size*dble(NPROCTOT)/1.d9,' GB'
  endif
  print *,'                                     = ',static_memory_size*dble(NPROCTOT)/1073741824.d0,' GiB'
  print *,'                                     = ',static_memory_size*dble(NPROCTOT)/1.d12,' TB'
  print *,'                                     = ',static_memory_size*dble(NPROCTOT)/1099511627776.d0,' TiB'
  print *

  print *,'How much memory (in GB) is installed on your machine per CPU core?'
  print *,'        (or per GPU card or per INTEL MIC Phi board)?'
  print *,'  (beware, this value MUST be given per core, i.e. per MPI thread, i.e. per MPI rank, NOT per node)'
  print *,'  (this value is for instance 4 GB on Tiger at Princeton, 2 GB on the non-GPU part of Titan at ORNL'
  print *,'   i.e. when using CPUs only there, and 1.5 GB on the GPU cluster in Marseille)'
  read(*,*) gigabytes_avail_per_core

  if(gigabytes_avail_per_core < 0.1d0) stop 'less than 100 MB per core does not seem realistic; exiting...'
  if(gigabytes_avail_per_core > 100.d0) stop 'more than 100 GB per core does not seem realistic; exiting...'

  print *
  print *,'What percentage of this total do you allow us to use, keeping in mind that you'
  print *,'need to leave some memory available for the GNU/Linux system to run?'
  print *,'  (a typical value is 90%; 92% to 95% is probably OK too; 85% is very safe)'
  read(*,*) percentage_to_use_per_core

  if(percentage_to_use_per_core < 50.d0) stop 'less than 50% does not seem realistic; exiting...'
  if(percentage_to_use_per_core > 96.d0) stop 'more than 96% is risky; exiting...'

  what_we_can_use_in_GB = gigabytes_avail_per_core * percentage_to_use_per_core / 100.d0

! convert static_memory_size to GB
  static_memory_size = static_memory_size / 1.d9

  print *
  print *,'without undoing of attenuation you are using ',static_memory_size,' GB per core'
  print *,'  i.e. ',sngl(100.d0 * static_memory_size / gigabytes_avail_per_core),'% of the installed memory'

  if(static_memory_size >= gigabytes_avail_per_core) &
    stop 'you are using more memory than what you told us is installed!!! there is an error'

  if(static_memory_size >= what_we_can_use_in_GB) &
    stop 'you are using more memory than what you allowed us to use!!! there is an error'

! compute the size to store in memory at each time step
  size_to_store_at_each_time_step = 0

! displ_crust_mantle
  size_to_store_at_each_time_step = size_to_store_at_each_time_step + dble(NDIM)*NGLOB(IREGION_CRUST_MANTLE)*dble(CUSTOM_REAL)

! displ_inner_core
  size_to_store_at_each_time_step = size_to_store_at_each_time_step + dble(NDIM)*NGLOB(IREGION_INNER_CORE)*dble(CUSTOM_REAL)

! displ_outer_core and accel_outer_core (both being scalar arrays)
  size_to_store_at_each_time_step = size_to_store_at_each_time_step + 2.d0*NGLOB(IREGION_OUTER_CORE)*dble(CUSTOM_REAL)

! convert to GB
  size_to_store_at_each_time_step = size_to_store_at_each_time_step / 1.d9

  print *
  print *,'each time step to store in memory to undo attenuation requires storing ',size_to_store_at_each_time_step,' GB per core'

  print *
  print *,'*******************************************************************************'
  print *,'the optimal value to put in DATA/Par_file is thus:'
  NT_DUMP_ATTENUATION_optimal_to_use = int((what_we_can_use_in_GB - static_memory_size) / size_to_store_at_each_time_step)
  print *
  print *,'NT_DUMP_ATTENUATION = ',NT_DUMP_ATTENUATION_optimal_to_use
  print *
  print *,'(no need to then recompile the code, just edit the file and change the value)'
  print *,'*******************************************************************************'

! compute the size of files to dump to disk
  disk_size_of_each_dumping = 0

! displ_crust_mantle, veloc_crust_mantle, accel_crust_mantle
  disk_size_of_each_dumping = disk_size_of_each_dumping + 3.d0*dble(NDIM)*NGLOB(IREGION_CRUST_MANTLE)*dble(CUSTOM_REAL)

! displ_inner_core, veloc_inner_core, accel_inner_core
  disk_size_of_each_dumping = disk_size_of_each_dumping + 3.d0*dble(NDIM)*NGLOB(IREGION_INNER_CORE)*dble(CUSTOM_REAL)

! displ_outer_core, veloc_outer_core, accel_outer_core (all scalar arrays)
  disk_size_of_each_dumping = disk_size_of_each_dumping + 3.d0*NGLOB(IREGION_OUTER_CORE)*dble(CUSTOM_REAL)

! A_array_rotation,B_array_rotation
  if (ROTATION) disk_size_of_each_dumping = disk_size_of_each_dumping + &
      dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_OUTER_CORE_ROTATION*2.d0*dble(CUSTOM_REAL)

  if (ATTENUATION) then
! R_memory_crust_mantle
    disk_size_of_each_dumping = disk_size_of_each_dumping + 5.d0*dble(N_SLS)*dble(NGLLX)* &
      dble(NGLLY)*dble(NGLLZ)*NSPEC_CRUST_MANTLE_ATTENUAT*dble(CUSTOM_REAL)

! R_memory_inner_core
    disk_size_of_each_dumping = disk_size_of_each_dumping + 5.d0*dble(N_SLS)*dble(NGLLX)* &
      dble(NGLLY)*dble(NGLLZ)*NSPEC_INNER_CORE_ATTENUATION*dble(CUSTOM_REAL)
  endif

! convert to GB
  disk_size_of_each_dumping = disk_size_of_each_dumping / 1.d9

!! DK DK this formula could be made more precise here; currently in some cases it can probably be off by +1 or -1
  number_of_dumpings_to_do = nint(NSTEP / dble(NT_DUMP_ATTENUATION_optimal_to_use))

  print *
  print *,'we will need to save a total of ',number_of_dumpings_to_do,' dumpings (restart files) to disk'

  print *
  print *,'each dumping on the disk to undo attenuation requires storing ',disk_size_of_each_dumping,' GB per core'

  print *
  print *,'each dumping on the disk requires storing ',disk_size_of_each_dumping*NPROCTOT, &
               ' GB for all cores'

  print *
  print *,'ALL dumpings on the disk require storing ',disk_size_of_each_dumping*number_of_dumpings_to_do,' GB per core'

  print *
  print *,'*******************************************************************************'
  print *,'ALL dumpings on the disk require storing ', &
               disk_size_of_each_dumping*number_of_dumpings_to_do*NPROCTOT,' GB for all cores'
  print *,'  i.e. ',disk_size_of_each_dumping*number_of_dumpings_to_do*NPROCTOT/1000.d0,' TB'
  print *,'*******************************************************************************'
  print *

  if(.not. UNDO_ATTENUATION) then
    print *
    print *,'*******************************************************************************'
    print *,'BEWARE, UNDO_ATTENUATION is .false. and thus undoing is currently'
    print *,'turned off, i.e. the above estimates are currently NOT USED.'
    print *,'*******************************************************************************'
    print *
  endif

  end program xcompute_optimized_dumping

