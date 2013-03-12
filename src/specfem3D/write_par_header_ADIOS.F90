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

!
!-------------------------------------------------------------------------------------------------
!

!> @brief Write data found in Par_file into ADIOS result file header.
!! Write the ADIOS header containing values to ensure reproductibility of
!! the simulation. These values come form the following files :
!! DATA/Par_file, DATA/CMTSOLUTION, DATA/STATIONS
subroutine write_par_file_header_ADIOS ()
  use mpi
  use adios_write_mod
  use specfem_par, only : myrank, NSOURCES

  implicit none
  include "constants.h"

  !-------------------------------------------------------------------
  ! local parameters
  !-------------------------------------------------------------------
  ! parameters read from parameter file (cf. DATA/Par_file)
  integer  :: NTSTEP_BETWEEN_OUTPUT_SEISMOS,NTSTEP_BETWEEN_READ_ADJSRC,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE, &
          MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP, &
          NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,NOISE_TOMOGRAPHY

  double precision :: ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,&
          CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,&
          HDUR_MOVIE,MOVIE_TOP_KM,MOVIE_BOTTOM_KM, &
          MOVIE_EAST_DEG,MOVIE_WEST_DEG,MOVIE_NORTH_DEG,&
          MOVIE_SOUTH_DEG,RECORD_LENGTH_IN_MINUTES

  logical :: ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS,&
         MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_COARSE, &
         RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
         SAVE_MESH_FILES,ATTENUATION,ATTENUATION_NEW, &
         ABSORBING_CONDITIONS,SAVE_FORWARD, &
         OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
         ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER,&
         SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE
  ! values from CMTSOLUTION
  ! integer          :: NSOURCES -> in specfem_par module
  integer,           dimension(NSOURCES) :: yr, mo, da, ho, mi
  double precision,  dimension(NSOURCES) :: sec, t_shift, hdur, lat, long, depth
  double precision,  dimension(NSOURCES) :: mrr, mtt, mpp, mrt, mrp, mtp 
  integer :: datasource_length
  character(len=5) :: datasource_tmp
  character(len=16):: event_name
  character(len=:), allocatable  :: datasource
  !character(len=256) :: datasource

  ! values from STATIONS
  integer :: NSTATIONS
  integer :: station_name_length, network_name_length
  character(len=MAX_LENGTH_STATION_NAME) :: station_name_tmp
  character(len=MAX_LENGTH_NETWORK_NAME) :: network_name_tmp
  character(len=:),   allocatable :: station_name, network_name 
  double precision, allocatable, dimension(:) :: stlat, stlon, stele, stbur

  character(len=150) :: OUTPUT_FILES,LOCAL_PATH,LOCAL_TMP_PATH,MODEL
  character(len=256) :: string, CMTSOLUTION, STATIONS

  integer                 :: isource, irec, ier
  ! Adios variables
  integer                 :: adios_err
  integer*8               :: adios_group, adios_handle, varid
  integer*8               :: adios_groupsize, adios_totalsize ! for the .fh file
  ! TODO : find a better name once the use of ADIOS is more completely
  ! implemented
  character(len=27)       :: filename = "OUTPUT_FILES/header_tmp.bp" 
  integer*8 :: group_size_inc


  ! ensure that only the master open the adios handle inside MPI_COMM_SELF
  if(myrank == 0) then
    call adios_declare_group (adios_group, "SPECFEM3D_GLOBE_HEADER", "", 0, adios_err)
    call adios_select_method (adios_group, "MPI", "", "", adios_err)

    group_size_inc = 0
    !--*** Values read from DATA/Par_file ***
    ! extract all unmodified values from the Par_file
    call read_parameter_file(OUTPUT_FILES,                                                                    & 
                             LOCAL_PATH,LOCAL_TMP_PATH,MODEL,                                                 & 
                             NTSTEP_BETWEEN_OUTPUT_SEISMOS,NTSTEP_BETWEEN_READ_ADJSRC,NTSTEP_BETWEEN_FRAMES,  & 
                             NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,                                       & 
                             NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE,                                      & 
                             MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP,                                        & 
                             NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,                                               & 
                             ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,                        & 
                             CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,   & 
                             HDUR_MOVIE,MOVIE_TOP_KM,MOVIE_BOTTOM_KM,RECORD_LENGTH_IN_MINUTES,                & 
                             MOVIE_EAST_DEG,MOVIE_WEST_DEG,MOVIE_NORTH_DEG,MOVIE_SOUTH_DEG,                   & 
                             ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS,                                  & 
                             MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_COARSE,                                         & 
                             RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION,                              & 
                             SAVE_MESH_FILES,ATTENUATION,ATTENUATION_NEW,ABSORBING_CONDITIONS,SAVE_FORWARD,   & 
                             OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, & 
                             ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER,                               & 
                             SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,NOISE_TOMOGRAPHY)
    ! define adios variables for the Par_file
    !-- double precision variables
    call define_adios_double_scalar (adios_group, "ANGULAR_WIDTH_XI_IN_DEGREES", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_double_scalar (adios_group, "ANGULAR_WIDTH_ETA_IN_DEGREES", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_double_scalar (adios_group, "CENTER_LONGITUDE_IN_DEGREES", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_double_scalar (adios_group, "CENTER_LATITUDE_IN_DEGREES", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_double_scalar (adios_group, "GAMMA_ROTATION_AZIMUTH", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_double_scalar (adios_group, "HDUR_MOVIE", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_double_scalar (adios_group, "MOVIE_TOP_KM", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_double_scalar (adios_group, "MOVIE_BOTTOM_KM", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_double_scalar (adios_group, "MOVIE_EAST_DEG", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_double_scalar (adios_group, "MOVIE_WEST_DEG", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_double_scalar (adios_group, "MOVIE_NORTH_DEG", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_double_scalar (adios_group, "MOVIE_SOUTH_DEG", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_double_scalar (adios_group, "RECORD_LENGTH_IN_MINUTES", "/specfem3D_globe_parameter_file", group_size_inc)
    !-- integer variables 
    call define_adios_integer_scalar (adios_group, "NTSTEP_BETWEEN_OUTPUT_SEISMOS", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_integer_scalar (adios_group, "NTSTEP_BETWEEN_READ_ADJSRC", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_integer_scalar (adios_group, "NTSTEP_BETWEEN_FRAMES", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_integer_scalar (adios_group, "NTSTEP_BETWEEN_OUTPUT_INFO", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_integer_scalar (adios_group, "NUMBER_OF_RUNS", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_integer_scalar (adios_group, "NUMBER_OF_THIS_RUN", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_integer_scalar (adios_group, "NCHUNKS", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_integer_scalar (adios_group, "SIMULATION_TYPE", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_integer_scalar (adios_group, "MOVIE_VOLUME_TYPE", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_integer_scalar (adios_group, "MOVIE_START", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_integer_scalar (adios_group, "MOVIE_STOP", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_integer_scalar (adios_group, "NEX_XI", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_integer_scalar (adios_group, "NEX_ETA", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_integer_scalar (adios_group, "NPROC_XI", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_integer_scalar (adios_group, "NPROC_ETA", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_integer_scalar (adios_group, "NOISE_TOMOGRAPHY", "/specfem3D_globe_parameter_file", group_size_inc)
    !-- logical variables
    call define_adios_byte_scalar (adios_group, "ELLIPTICITY", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "GRAVITY", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "ROTATION", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "TOPOGRAPHY", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "OCEANS", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "MOVIE_SURFACE", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "MOVIE_VOLUME", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "MOVIE_COARSE", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "RECEIVERS_CAN_BE_BURIED", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "PRINT_SOURCE_TIME_FUNCTION", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "SAVE_MESH_FILES", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "ATTENUATION", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "ATTENUATION_NEW", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "ABSORBING_CONDITIONS", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "SAVE_FORWARD", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "OUTPUT_SEISMOS_ASCII_TEXT", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "OUTPUT_SEISMOS_SAC_ALPHANUM", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "OUTPUT_SEISMOS_SAC_BINARY", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "ROTATE_SEISMOGRAMS_RT", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "WRITE_SEISMOGRAMS_BY_MASTER", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "SAVE_ALL_SEISMOS_IN_ONE_FILE", "/specfem3D_globe_parameter_file", group_size_inc)
    call define_adios_byte_scalar (adios_group, "USE_BINARY_FOR_LARGE_FILE", "/specfem3D_globe_parameter_file", group_size_inc)

    !--*** Values read from DATA/CMTSOLUTION ***--
    ! extract all unmodified values from CMTSOLUTION
    ! get_cmt() routine modify the read values
    ! TODO factorize what follows and get_cmt.f90 and probably one or two other
    !      routines
    call get_value_string(CMTSOLUTION, 'solver.CMTSOLUTION', 'DATA/CMTSOLUTION')
    open(unit=1,file=CMTSOLUTION,status='old',action='read')
    datasource_length = 4*NSOURCES
    allocate(character(len=(datasource_length)) :: datasource, stat=ier)
    if (ier /=0) &
        call exit_MPI (myrank, "error allocating datasource string for adios header")
    datasource = ""
    ! ADIOS only  (1) byte for a string. This may cause data overwriting.
    ! => increase the generate by the string size -1
    adios_groupsize = adios_groupsize + 4*NSOURCES - 1
    do isource=1,NSOURCES

      read(1,"(a256)") string
      ! skips empty lines
      do while( len_trim(string) == 0 )
      read(1,"(a256)") string
      enddo

      ! read header with event information -- TODO there is more to add
      !read(string,"(a4,i5,i3,i3,i3,i3,f6.2)") datasource(:,isource),yr(isource), &
      read(string,"(a4,i5,i3,i3,i3,i3,f6.2)") datasource_tmp,yr(isource), &
          mo(isource),da(isource),ho(isource),mi(isource),sec(isource)
      datasource = datasource // datasource_tmp
      ! read event name 
      read(1,"(a)") string
      read(string(12:len_trim(string)),*) event_name
      ! read time shift
      read(1,"(a)") string
      read(string(12:len_trim(string)),*) t_shift(isource)
      ! read half duration
      read(1,"(a)") string
      read(string(15:len_trim(string)),*) hdur(isource)
      ! read latitude
      read(1,"(a)") string
      read(string(10:len_trim(string)),*) lat(isource)
      ! read longitude
      read(1,"(a)") string
      read(string(11:len_trim(string)),*) long(isource)
      ! read depth
      read(1,"(a)") string
      read(string(7:len_trim(string)),*) depth(isource)
      ! read Mrr
      read(1,"(a)") string
      read(string(5:len_trim(string)),*) mrr(isource)
      ! read Mtt
      read(1,"(a)") string
      read(string(5:len_trim(string)),*) mtt(isource)
      ! read Mpp
      read(1,"(a)") string
      read(string(5:len_trim(string)),*) mpp(isource)
      ! read Mrt
      read(1,"(a)") string
      read(string(5:len_trim(string)),*) mrt(isource)
      ! read Mrp
      read(1,"(a)") string
      read(string(5:len_trim(string)),*) mrp(isource)
      ! read Mtp
      read(1,"(a)") string
      read(string(5:len_trim(string)),*) mtp(isource)
    enddo
    close(1)
    ! define adios variables for the CMTSOLUTION 
    !-- Number of SOURCES inside the CMTSOLUTION file
    call define_adios_integer_scalar (adios_group, "NSOURCES", "/CMTSOLUTION", group_size_inc)
    !-- double precision arrays
    call define_adios_double_local_array1D (adios_group, "second", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
    call define_adios_double_local_array1D (adios_group, "time_shift", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
    call define_adios_double_local_array1D (adios_group, "half_duration", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
    call define_adios_double_local_array1D (adios_group, "latitude", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
    call define_adios_double_local_array1D (adios_group, "longitude", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
    call define_adios_double_local_array1D (adios_group, "depth", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
    call define_adios_double_local_array1D (adios_group, "mrr", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
    call define_adios_double_local_array1D (adios_group, "mtt", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
    call define_adios_double_local_array1D (adios_group, "mpp", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
    call define_adios_double_local_array1D (adios_group, "mrt", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
    call define_adios_double_local_array1D (adios_group, "mrp", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
    call define_adios_double_local_array1D (adios_group, "mtp", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
    !-- integer arrays
    call define_adios_integer_local_array1D (adios_group, "year", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
    call define_adios_integer_local_array1D (adios_group, "month", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
    call define_adios_integer_local_array1D (adios_group, "day", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
    call define_adios_integer_local_array1D (adios_group, "hour", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
    call define_adios_integer_local_array1D (adios_group, "minute", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
    !-- string
    call define_adios_integer_scalar (adios_group, "datasource_length", "/CMTSOLUTION", group_size_inc)
    call define_adios_string (adios_group, "datasource", "/CMTSOLUTION", datasource_length, group_size_inc)


    !--*** Values read from DATA/STATIONS
    ! Extract values from STATIONS File
    call get_value_string(STATIONS, 'solver.STATIONS', 'DATA/STATIONS')
    open(unit=1,file=STATIONS,iostat=ier,status='old',action='read')
    NSTATIONS = 0
    do while(ier == 0)
      read(1,"(a)",iostat=ier) string
      if(ier == 0) NSTATIONS = NSTATIONS + 1
    enddo
    allocate (character (len=(MAX_LENGTH_STATION_NAME*NSTATIONS)) :: station_name)
    allocate (character (len=(MAX_LENGTH_NETWORK_NAME*NSTATIONS)) :: network_name)
    allocate (stlat (NSTATIONS)) 
    allocate (stlon (NSTATIONS)) 
    allocate (stele (NSTATIONS)) 
    allocate (stbur (NSTATIONS)) 
    station_name = ""
    network_name = ""
    rewind(1)
    do irec = 1,NSTATIONS
      !read(1,*,iostat=ier) station_name(irec),network_name(irec),stlat(irec),stlon(irec),stele(irec),stbur(irec)
      read(1,*,iostat=ier) station_name_tmp, network_name_tmp, &
                            stlat(irec),      stlon(irec),      &
                            stele(irec),      stbur(irec)
      if( ier /= 0 ) then
        write(IMAIN,*) 'error reading in station ',irec
        call exit_MPI(myrank,'error reading in station in STATIONS file')
      endif
      station_name = station_name // trim(station_name_tmp) // " " 
      network_name = network_name // trim(network_name_tmp) // " " 
    enddo
    close(1)
    station_name = trim(station_name)
    network_name = trim(network_name)
    station_name_length = len(station_name)
    network_name_length = len(network_name)

    !-- Number of STATIONS inside the STATIONS file 
    call define_adios_integer_scalar (adios_group, "NSTATIONS", "/STATIONS", group_size_inc)
   !-- double precision arrays 
    call define_adios_double_local_array1D (adios_group, "station_latitude", "/STATIONS", NSTATIONS, "NSTATIONS", group_size_inc)
    call define_adios_double_local_array1D (adios_group, "station_longitude", "/STATIONS", NSTATIONS, "NSTATIONS", group_size_inc)
    call define_adios_double_local_array1D (adios_group, "station_elevation", "/STATIONS", NSTATIONS, "NSTATIONS", group_size_inc)
    call define_adios_double_local_array1D (adios_group, "station_burial", "/STATIONS", NSTATIONS, "NSTATIONS", group_size_inc)
    !-- string
    call define_adios_integer_scalar (adios_group, "station_name_length", "/STATIONS", group_size_inc)
    call define_adios_integer_scalar (adios_group, "network_name_length", "/STATIONS", group_size_inc)
    call define_adios_string (adios_group, "station_name", "/CMTSOLUTION", station_name_length, group_size_inc)
    call define_adios_string (adios_group, "network_name", "/CMTSOLUTION", network_name_length, group_size_inc)

    call adios_open (adios_handle, "SPECFEM3D_GLOBE_HEADER", filename, "w", &
                     MPI_COMM_SELF, adios_err);

    ! Automatically generated Fortran code by ADIOS gpp.py script:
    ! $> $ADIOS_DIR/bin/gpp.py par_header.xml
    ! Re-run each time the xml file is modified.
    ! #include "gwrite_SPECFEM3D_GLOBE_HEADER.fh"
    adios_groupsize = group_size_inc

    call adios_group_size (adios_handle, adios_groupsize, adios_totalsize, adios_err)
    call adios_write (adios_handle, "ANGULAR_WIDTH_XI_IN_DEGREES", ANGULAR_WIDTH_XI_IN_DEGREES, adios_err)
    call adios_write (adios_handle, "ANGULAR_WIDTH_ETA_IN_DEGREES", ANGULAR_WIDTH_ETA_IN_DEGREES, adios_err)
    call adios_write (adios_handle, "CENTER_LONGITUDE_IN_DEGREES", CENTER_LONGITUDE_IN_DEGREES, adios_err)
    call adios_write (adios_handle, "CENTER_LATITUDE_IN_DEGREES", CENTER_LATITUDE_IN_DEGREES, adios_err)
    call adios_write (adios_handle, "GAMMA_ROTATION_AZIMUTH", GAMMA_ROTATION_AZIMUTH, adios_err)
    call adios_write (adios_handle, "HDUR_MOVIE", HDUR_MOVIE, adios_err)
    call adios_write (adios_handle, "MOVIE_TOP_KM", MOVIE_TOP_KM, adios_err)
    call adios_write (adios_handle, "MOVIE_BOTTOM_KM", MOVIE_BOTTOM_KM, adios_err)
    call adios_write (adios_handle, "MOVIE_EAST_DEG", MOVIE_EAST_DEG, adios_err)
    call adios_write (adios_handle, "MOVIE_WEST_DEG", MOVIE_WEST_DEG, adios_err)
    call adios_write (adios_handle, "MOVIE_NORTH_DEG", MOVIE_NORTH_DEG, adios_err)
    call adios_write (adios_handle, "MOVIE_SOUTH_DEG", MOVIE_SOUTH_DEG, adios_err)
    call adios_write (adios_handle, "RECORD_LENGTH_IN_MINUTES", RECORD_LENGTH_IN_MINUTES, adios_err)
    call adios_write (adios_handle, "NTSTEP_BETWEEN_OUTPUT_SEISMOS", NTSTEP_BETWEEN_OUTPUT_SEISMOS, adios_err)
    call adios_write (adios_handle, "NTSTEP_BETWEEN_READ_ADJSRC", NTSTEP_BETWEEN_READ_ADJSRC, adios_err)
    call adios_write (adios_handle, "NTSTEP_BETWEEN_FRAMES", NTSTEP_BETWEEN_FRAMES, adios_err)
    call adios_write (adios_handle, "NTSTEP_BETWEEN_OUTPUT_INFO", NTSTEP_BETWEEN_OUTPUT_INFO, adios_err)
    call adios_write (adios_handle, "NUMBER_OF_RUNS", NUMBER_OF_RUNS, adios_err)
    call adios_write (adios_handle, "NUMBER_OF_THIS_RUN", NUMBER_OF_THIS_RUN, adios_err)
    call adios_write (adios_handle, "NCHUNKS", NCHUNKS, adios_err)
    call adios_write (adios_handle, "SIMULATION_TYPE", SIMULATION_TYPE, adios_err)
    call adios_write (adios_handle, "MOVIE_VOLUME_TYPE", MOVIE_VOLUME_TYPE, adios_err)
    call adios_write (adios_handle, "MOVIE_START", MOVIE_START, adios_err)
    call adios_write (adios_handle, "MOVIE_STOP", MOVIE_STOP, adios_err)
    call adios_write (adios_handle, "NEX_XI", NEX_XI, adios_err)
    call adios_write (adios_handle, "NEX_ETA", NEX_ETA, adios_err)
    call adios_write (adios_handle, "NPROC_XI", NPROC_XI, adios_err)
    call adios_write (adios_handle, "NPROC_ETA", NPROC_ETA, adios_err)
    call adios_write (adios_handle, "NOISE_TOMOGRAPHY", NOISE_TOMOGRAPHY, adios_err)
    call adios_write (adios_handle, "ELLIPTICITY", ELLIPTICITY, adios_err)
    call adios_write (adios_handle, "GRAVITY", GRAVITY, adios_err)
    call adios_write (adios_handle, "ROTATION", ROTATION, adios_err)
    call adios_write (adios_handle, "TOPOGRAPHY", TOPOGRAPHY, adios_err)
    call adios_write (adios_handle, "OCEANS", OCEANS, adios_err)
    call adios_write (adios_handle, "MOVIE_SURFACE", MOVIE_SURFACE, adios_err)
    call adios_write (adios_handle, "MOVIE_VOLUME", MOVIE_VOLUME, adios_err)
    call adios_write (adios_handle, "MOVIE_COARSE", MOVIE_COARSE, adios_err)
    call adios_write (adios_handle, "RECEIVERS_CAN_BE_BURIED", RECEIVERS_CAN_BE_BURIED, adios_err)
    call adios_write (adios_handle, "PRINT_SOURCE_TIME_FUNCTION", PRINT_SOURCE_TIME_FUNCTION, adios_err)
    call adios_write (adios_handle, "SAVE_MESH_FILES", SAVE_MESH_FILES, adios_err)
    call adios_write (adios_handle, "ATTENUATION", ATTENUATION, adios_err)
    call adios_write (adios_handle, "ATTENUATION_NEW", ATTENUATION_NEW, adios_err)
    call adios_write (adios_handle, "ABSORBING_CONDITIONS", ABSORBING_CONDITIONS, adios_err)
    call adios_write (adios_handle, "SAVE_FORWARD", SAVE_FORWARD, adios_err)
    call adios_write (adios_handle, "OUTPUT_SEISMOS_ASCII_TEXT", OUTPUT_SEISMOS_ASCII_TEXT, adios_err)
    call adios_write (adios_handle, "OUTPUT_SEISMOS_SAC_ALPHANUM", OUTPUT_SEISMOS_SAC_ALPHANUM, adios_err)
    call adios_write (adios_handle, "OUTPUT_SEISMOS_SAC_BINARY", OUTPUT_SEISMOS_SAC_BINARY, adios_err)
    call adios_write (adios_handle, "ROTATE_SEISMOGRAMS_RT", ROTATE_SEISMOGRAMS_RT, adios_err)
    call adios_write (adios_handle, "WRITE_SEISMOGRAMS_BY_MASTER", WRITE_SEISMOGRAMS_BY_MASTER, adios_err)
    call adios_write (adios_handle, "SAVE_ALL_SEISMOS_IN_ONE_FILE", SAVE_ALL_SEISMOS_IN_ONE_FILE, adios_err)
    call adios_write (adios_handle, "USE_BINARY_FOR_LARGE_FILE", USE_BINARY_FOR_LARGE_FILE, adios_err)
    call adios_write (adios_handle, "NSOURCES", NSOURCES, adios_err)
    call adios_write (adios_handle, "year", yr, adios_err)
    call adios_write (adios_handle, "month", mo, adios_err)
    call adios_write (adios_handle, "day", da, adios_err)
    call adios_write (adios_handle, "hour", ho, adios_err)
    call adios_write (adios_handle, "minute", mi, adios_err)
    call adios_write (adios_handle, "second", sec, adios_err)
    call adios_write (adios_handle, "time_shift", t_shift, adios_err)
    call adios_write (adios_handle, "half_duration", hdur, adios_err)
    call adios_write (adios_handle, "latitude", lat, adios_err)
    call adios_write (adios_handle, "longitude", long, adios_err)
    call adios_write (adios_handle, "depth", depth, adios_err)
    call adios_write (adios_handle, "mrr", mrr, adios_err)
    call adios_write (adios_handle, "mtt", mtt, adios_err)
    call adios_write (adios_handle, "mpp", mpp, adios_err)
    call adios_write (adios_handle, "mrt", mrt, adios_err)
    call adios_write (adios_handle, "mrp", mrp, adios_err)
    call adios_write (adios_handle, "mtp", mtp, adios_err)
    call adios_write (adios_handle, "datasource_length", datasource_length, adios_err)
    call adios_write (adios_handle, "datasource", datasource, adios_err)

    call adios_write (adios_handle, "NSTATIONS", NSTATIONS, adios_err)
    call adios_write (adios_handle, "station_latitude", stlat, adios_err)
    call adios_write (adios_handle, "station_longitude", stlon, adios_err)
    call adios_write (adios_handle, "station_elevation", stele, adios_err)
    call adios_write (adios_handle, "station_burial", stbur, adios_err)
    call adios_write (adios_handle, "station_name_length", station_name_length, adios_err)
    call adios_write (adios_handle, "network_name_length", network_name_length, adios_err)
    call adios_write (adios_handle, "station_name", station_name, adios_err)
    call adios_write (adios_handle, "network_name", network_name, adios_err)

    call adios_close (adios_handle, adios_err)

    deallocate(datasource)
    deallocate(station_name)
    deallocate(network_name)
    deallocate(stlat)
    deallocate(stlon)
    deallocate(stele)
    deallocate(stbur)
  endif

end subroutine write_par_file_header_ADIOS


!subroutine define_adios_header_par_file (adios_group, group_size_inc)
!  implicit none
!end subroutine define_adios_header_par_file


subroutine define_adios_double_scalar (adios_group, name, path, group_size_inc)
  use adios_write_mod
  implicit none
  ! Arguments
  integer(kind=8),  intent(in)   :: adios_group
  character(len=*), intent(in)   :: name, path
  integer(kind=8),  intent(out)  :: group_size_inc
  ! Local Variables
  integer(kind=8)                :: varid ! dummy variable, adios use var name

  call adios_define_var (adios_group, trim(name), trim(path), 6,  "", "", "", varid)
  group_size_inc = group_size_inc + 8
end subroutine define_adios_double_scalar

subroutine define_adios_integer_scalar (adios_group, name, path, group_size_inc)
  use adios_write_mod
  implicit none
  ! Arguments
  integer(kind=8),  intent(in)   :: adios_group
  character(len=*), intent(in)   :: name, path
  integer(kind=8),  intent(out)  :: group_size_inc
  ! Local Variables
  integer(kind=8)                :: varid ! dummy variable, adios use var name

  call adios_define_var (adios_group, trim(name), trim(path), 2,  "", "", "", varid)
  group_size_inc = group_size_inc + 4
end subroutine define_adios_integer_scalar

subroutine define_adios_byte_scalar (adios_group, name, path, group_size_inc)
  use adios_write_mod
  implicit none
  ! Arguments
  integer(kind=8),  intent(in)   :: adios_group
  character(len=*), intent(in)   :: name, path
  integer(kind=8),  intent(out)  :: group_size_inc
  ! Local Variables
  integer(kind=8)                :: varid ! dummy variable, adios use var name

  call adios_define_var (adios_group, trim(name), trim(path), 0,  "", "", "", varid)
  group_size_inc = group_size_inc + 1
end subroutine define_adios_byte_scalar

subroutine define_adios_integer_local_array1D (adios_group, name, path, dim, dim_str, group_size_inc)
  use adios_write_mod
  implicit none
  ! Arguments
  integer(kind=8),  intent(in)   :: adios_group
  character(len=*), intent(in)   :: name, path, dim_str
  integer(kind=8),  intent(out)  :: group_size_inc
  integer, intent(in)            :: dim
  ! Local Variables
  integer(kind=8)                :: varid ! dummy variable, adios use var name

  call adios_define_var (adios_group, name, path, 2,  dim_str, "", "", varid)
  group_size_inc = group_size_inc + 4*dim
end subroutine define_adios_integer_local_array1D

subroutine define_adios_double_local_array1D (adios_group, name, path, dim, dim_str, group_size_inc)
  use adios_write_mod
  implicit none
  ! Arguments
  integer(kind=8),  intent(in)   :: adios_group
  character(len=*), intent(in)   :: name, path, dim_str
  integer(kind=8),  intent(out)  :: group_size_inc
  integer, intent(in)            :: dim
  ! Local Variables
  integer(kind=8)                :: varid ! dummy variable, adios use var name

  call adios_define_var (adios_group, name, path, 6, dim_str, "", "", varid)
  group_size_inc = group_size_inc + 8*dim
end subroutine define_adios_double_local_array1D

subroutine define_adios_string (adios_group, name, path, length, group_size_inc)
  use adios_write_mod
  implicit none
  ! Arguments
  integer(kind=8),  intent(in)   :: adios_group
  character(len=*), intent(in)   :: name, path
  integer(kind=8),  intent(out)  :: group_size_inc
  integer                        :: length
  ! Local Variables
  integer(kind=8)                :: varid ! dummy variable, adios use var name

  call adios_define_var (adios_group, name, path, 9,  "", "", "", varid)
  group_size_inc = group_size_inc + 1*length 
end subroutine define_adios_string
