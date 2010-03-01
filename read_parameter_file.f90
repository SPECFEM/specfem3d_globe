!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
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

  subroutine read_parameter_file(OUTPUT_FILES,LOCAL_PATH,MODEL, &
                                NTSTEP_BETWEEN_OUTPUT_SEISMOS,NTSTEP_BETWEEN_READ_ADJSRC,NTSTEP_BETWEEN_FRAMES, &
                                NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS, &
                                NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE, &
                                MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP, &
                                NEX_XI_read,NEX_ETA_read,NPROC_XI_read,NPROC_ETA_read, &
                                ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,&
                                CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,&
                                HDUR_MOVIE,MOVIE_TOP_KM,MOVIE_BOTTOM_KM,RECORD_LENGTH_IN_MINUTES, &
                                MOVIE_EAST_DEG,MOVIE_WEST_DEG,MOVIE_NORTH_DEG,MOVIE_SOUTH_DEG,&
                                ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS,&
                                MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_COARSE, &
                                RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
                                SAVE_MESH_FILES,ATTENUATION,ABSORBING_CONDITIONS,SAVE_FORWARD, &
                                OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
                                ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER,&
                                SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE)

  implicit none

  include "constants.h"

! parameters read from parameter file
  integer NTSTEP_BETWEEN_OUTPUT_SEISMOS,NTSTEP_BETWEEN_READ_ADJSRC,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE, &
          MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP, &
          NEX_XI_read,NEX_ETA_read,NPROC_XI_read,NPROC_ETA_read

  double precision ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,&
          CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,&
          HDUR_MOVIE,MOVIE_TOP_KM,MOVIE_BOTTOM_KM, &
          MOVIE_EAST_DEG,MOVIE_WEST_DEG,MOVIE_NORTH_DEG,&
          MOVIE_SOUTH_DEG,RECORD_LENGTH_IN_MINUTES

  logical ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS,&
         MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_COARSE, &
         RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
         SAVE_MESH_FILES,ATTENUATION, &
         ABSORBING_CONDITIONS,SAVE_FORWARD, &
         OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
         ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER,&
         SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE

  character(len=150) OUTPUT_FILES,LOCAL_PATH,MODEL

! local variables
  integer, external :: err_occurred

  ! gets the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

  ! opens the parameter file: DATA/Par_file
  call open_parameter_file()

  ! reads in values
  call read_value_integer(SIMULATION_TYPE, 'solver.SIMULATION_TYPE')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: SIMULATION_TYPE'
  call read_value_logical(SAVE_FORWARD, 'solver.SAVE_FORWARD')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: SAVE_FORWARD'
  call read_value_integer(NCHUNKS, 'mesher.NCHUNKS')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: NCHUNKS'
  call read_value_double_precision(ANGULAR_WIDTH_XI_IN_DEGREES, 'mesher.ANGULAR_WIDTH_XI_IN_DEGREES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: ANGULAR_WIDTH_XI...'
  call read_value_double_precision(ANGULAR_WIDTH_ETA_IN_DEGREES, 'mesher.ANGULAR_WIDTH_ETA_IN_DEGREES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: ANGULAR_WIDTH_ETA...'
  call read_value_double_precision(CENTER_LATITUDE_IN_DEGREES, 'mesher.CENTER_LATITUDE_IN_DEGREES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: CENTER_LATITUDE...'
  call read_value_double_precision(CENTER_LONGITUDE_IN_DEGREES, 'mesher.CENTER_LONGITUDE_IN_DEGREES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: CENTER_LONGITUDE...'
  call read_value_double_precision(GAMMA_ROTATION_AZIMUTH, 'mesher.GAMMA_ROTATION_AZIMUTH')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: GAMMA_ROTATION...'
  ! number of elements at the surface along the two sides of the first chunk
  call read_value_integer(NEX_XI_read, 'mesher.NEX_XI')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: NEX_XI'
  call read_value_integer(NEX_ETA_read, 'mesher.NEX_ETA')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: NEX_ETA'
  call read_value_integer(NPROC_XI_read, 'mesher.NPROC_XI')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: NPROC_XI'
  call read_value_integer(NPROC_ETA_read, 'mesher.NPROC_ETA')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: NPROC_ETA'
  call read_value_logical(OCEANS, 'model.OCEANS')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: OCEANS'
  call read_value_logical(ELLIPTICITY, 'model.ELLIPTICITY')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: ELLIPTICITIY'
  call read_value_logical(TOPOGRAPHY, 'model.TOPOGRAPHY')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: TOPOGRAPHY'
  call read_value_logical(GRAVITY, 'model.GRAVITY')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: GRAVITY'
  call read_value_logical(ROTATION, 'model.ROTATION')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: ROTATION'
  call read_value_logical(ATTENUATION, 'model.ATTENUATION')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: ATTENUATION'
  call read_value_logical(ABSORBING_CONDITIONS, 'solver.ABSORBING_CONDITIONS')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: ABSORBING_CONDITIONS'
  ! define the velocity model
  call read_value_string(MODEL, 'model.MODEL')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: MODEL'
  call read_value_double_precision(RECORD_LENGTH_IN_MINUTES, 'solver.RECORD_LENGTH_IN_MINUTES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: RECORD_LENGTH..'
  call read_value_logical(MOVIE_SURFACE, 'solver.MOVIE_SURFACE')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: MOVIE_SURFACE'
  call read_value_logical(MOVIE_VOLUME, 'solver.MOVIE_VOLUME')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: MOVIE_VOLUME'
  call read_value_logical(MOVIE_COARSE,'solver.MOVIE_COARSE')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: MOVIE_COARSE'
  call read_value_integer(NTSTEP_BETWEEN_FRAMES, 'solver.NTSTEP_BETWEEN_FRAMES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: NTSTEP_BETWEEN_FRAMES'
  call read_value_double_precision(HDUR_MOVIE, 'solver.HDUR_MOVIE')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: HDUR_MOVIE'
  call read_value_integer(MOVIE_VOLUME_TYPE, 'solver.MOVIE_VOLUME_TYPE')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: MOVIE_VOLUME_TYPE'
  call read_value_double_precision(MOVIE_TOP_KM, 'solver.MOVIE_TOP_KM')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: MOVIE_TOP_KM'
  call read_value_double_precision(MOVIE_BOTTOM_KM, 'solver.MOVIE_BOTTOM_KM')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: MOVIE_BOTTOM_KM'
  call read_value_double_precision(MOVIE_WEST_DEG, 'solver.MOVIE_WEST_DEG')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: MOVIE_WEST_DEG'
  call read_value_double_precision(MOVIE_EAST_DEG, 'solver.MOVIE_EAST_DEG')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: MOVIE_EAST_DEG'
  call read_value_double_precision(MOVIE_NORTH_DEG, 'solver.MOVIE_NORTH_DEG')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: MOVIE_NORTH_DEG'
  call read_value_double_precision(MOVIE_SOUTH_DEG, 'solver.MOVIE_SOUTH_DEG')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: MOVIE_SOUTH_DEG'
  call read_value_integer(MOVIE_START, 'solver.MOVIE_START')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: MOVIE_START'
  call read_value_integer(MOVIE_STOP, 'solver.MOVIE_STOP')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: MOVIE_STOP'
  call read_value_logical(SAVE_MESH_FILES, 'mesher.SAVE_MESH_FILES')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: SAVE_MESH_FILES'
  call read_value_integer(NUMBER_OF_RUNS, 'solver.NUMBER_OF_RUNS')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: NUMBER_OF_RUNS'
  call read_value_integer(NUMBER_OF_THIS_RUN, 'solver.NUMBER_OF_THIS_RUN')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: NUMBER_OF_THIS_RUN'
  call read_value_string(LOCAL_PATH, 'LOCAL_PATH')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: LOCAL_PATH'
  call read_value_integer(NTSTEP_BETWEEN_OUTPUT_INFO, 'solver.NTSTEP_BETWEEN_OUTPUT_INFO')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: NTSTEP_BETWEEN_OUTPUT_INFO'
  call read_value_integer(NTSTEP_BETWEEN_OUTPUT_SEISMOS, 'solver.NTSTEP_BETWEEN_OUTPUT_SEISMOS')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: NTSTEP_BETWEEN_OUTPUT_SEISMOS'
  call read_value_integer(NTSTEP_BETWEEN_READ_ADJSRC, 'solver.NTSTEP_BETWEEN_READ_ADJSRC')
  if(err_occurred() /= 0) return
  call read_value_logical(OUTPUT_SEISMOS_ASCII_TEXT, 'solver.OUTPUT_SEISMOS_ASCII_TEXT')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: OUTPUT_SIESMOS_ASCII_TEXT'
  call read_value_logical(OUTPUT_SEISMOS_SAC_ALPHANUM, 'solver.OUTPUT_SEISMOS_SAC_ALPHANUM')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: OUTPUT_SEISMOS_SAC_ALPHANUM'
  call read_value_logical(OUTPUT_SEISMOS_SAC_BINARY, 'solver.OUTPUT_SEISMOS_SAC_BINARY')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: OUTPUT_SEISMOS_SAC_BINARY'
  call read_value_logical(ROTATE_SEISMOGRAMS_RT, 'solver.ROTATE_SEISMOGRAMS_RT')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: ROTATE_SEISMOGRAMS_RT'
  call read_value_logical(WRITE_SEISMOGRAMS_BY_MASTER, 'solver.WRITE_SEISMOGRAMS_BY_MASTER')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: WRITE_SEISMOGRAMS_BY_MASTER'
  call read_value_logical(SAVE_ALL_SEISMOS_IN_ONE_FILE, 'solver.SAVE_ALL_SEISMOS_IN_ONE_FILE')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: SAVE_ALL_SEISMOS_IN_ONE_FILE'
  call read_value_logical(USE_BINARY_FOR_LARGE_FILE, 'solver.USE_BINARY_FOR_LARGE_FILE')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: USE_BINARY_FOR_LARGE_FILE'
  call read_value_logical(RECEIVERS_CAN_BE_BURIED, 'solver.RECEIVERS_CAN_BE_BURIED')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: RECEIVERS_CAN_BE_BURIED'
  call read_value_logical(PRINT_SOURCE_TIME_FUNCTION, 'solver.PRINT_SOURCE_TIME_FUNCTION')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file: PRINT_SOURCE_TIME_FUNCTION'

  ! closes parameter file
  call close_parameter_file()

  end subroutine read_parameter_file
  
