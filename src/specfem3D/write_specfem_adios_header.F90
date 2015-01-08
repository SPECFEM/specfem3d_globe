!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
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


!> \file write_par_header_ADIOS.F90
!! \brief Write in the adios file a group with all the parameters that insure
!!        reproducibility

#include "config.fh"

!
!-------------------------------------------------------------------------------
!

! TODO routine not used yet...

!> @brief Write simulation parameters into ADIOS result file header.
!!
!! Write the ADIOS header containing values to ensure reproducibility of
!! the simulation. These values come form the following files :
!! DATA/Par_file, DATA/CMTSOLUTION, DATA/STATIONS
subroutine write_specfem_header_adios()

  use shared_parameters

  use adios_helpers_mod
  use adios_write_mod

  use specfem_par, only : myrank, NSOURCES,ADIOS_TRANSPORT_METHOD

  implicit none

  !-------------------------------------------------------------------
  ! local parameters
  !-------------------------------------------------------------------

  ! values from CMTSOLUTION -------------------------------
  ! integer          :: NSOURCES -> in specfem_par module
  integer,           dimension(NSOURCES) :: yr, mo, da, ho, mi
  double precision,  dimension(NSOURCES) :: sec, t_shift, hdur, lat, long, depth
  double precision,  dimension(NSOURCES) :: mrr, mtt, mpp, mrt, mrp, mtp
  integer :: event_name_length, datasource_length
  character(len=16):: event_name
  character(len=:), allocatable  :: datasource  ! F03 feature

  ! values from STATIONS ----------------------------------
  integer :: NSTATIONS
  integer :: station_name_length, network_name_length ! for later reading
  character(len=:),   allocatable :: station_name, network_name
  double precision, allocatable, dimension(:) :: stlat, stlon, stele, stbur

  ! Adios variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, adios_handle
  integer(kind=8)         :: adios_groupsize, adios_totalsize
  ! TODO: find a better name once the use of ADIOS is more completely
  !       implemented
  character(len=MAX_STRING_LEN):: filename = "OUTPUT_FILES/header_specfem3d_globe.bp"
  integer(kind=8)   :: group_size_inc
  integer           :: model_length ! for later reading of MODEL
  integer           :: isource, irec, ier
  integer :: comm

  ! only the master needs to read the values to be written
  if (myrank == 0) then
    call adios_declare_group (adios_group, "SPECFEM3D_GLOBE_HEADER","", 0, adios_err)
    ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
    !call check_adios_err(myrank,adios_err)

    call adios_select_method (adios_group, ADIOS_TRANSPORT_METHOD, "", "", adios_err)
    ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
    !call check_adios_err(myrank,adios_err)

    group_size_inc = 0 ! Adios group size. Incremented by adios_helpers

    !-- *** Define variables used to configure specfem
    call define_solver_info_variables (adios_group, group_size_inc)

    !--*** Values read from DATA/Par_file ***
    ! extract all unmodified values from the Par_file
    call read_parameter_file()

    model_length = len(MODEL)
    ! define adios variables for the Par_file
    call define_par_file_variables (adios_group, group_size_inc, model_length)

    !--*** Values read from DATA/CMTSOLUTION ***--
    call read_raw_cmtsolution (yr, mo, da, ho, mi, sec, t_shift, hdur, lat, &
        long, depth, mrr, mtt, mpp, mrt, mrp, mtp, event_name_length, &
        event_name, datasource_length,  datasource)
    call define_cmtsolution_variables (adios_group, group_size_inc, NSOURCES, &
        event_name_length, datasource_length)

    !--*** Values read from DATA/STATIONS
    call read_raw_stations (NSTATIONS, stlat, stlon, stele, stbur, &
        station_name_length, station_name, network_name_length, network_name)
    call define_stations_variables (adios_group, group_size_inc, NSTATIONS, &
        station_name_length, network_name_length)

    call world_get_comm_self(comm)
    ! open the file where the headers have to be written
    call adios_open (adios_handle, "SPECFEM3D_GLOBE_HEADER", filename, "w",comm, adios_err)
    if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed for SPECFEM3D_GLOBE_HEADER'

    ! The group size have been auto-incremented
    adios_groupsize = group_size_inc
    call adios_group_size (adios_handle, adios_groupsize,adios_totalsize, adios_err)
    if (adios_err /= 0 ) stop 'Error calling adios_group_size() routine failed'

    ! Write variables from 'config.h'
    call write_adios_solver_info_variables (adios_handle)

    ! Write variables from 'Par_file'
    call write_adios_par_file_variables (adios_handle, &
        ANGULAR_WIDTH_XI_IN_DEGREES, ANGULAR_WIDTH_ETA_IN_DEGREES, &
        CENTER_LONGITUDE_IN_DEGREES, CENTER_LATITUDE_IN_DEGREES, &
        GAMMA_ROTATION_AZIMUTH, HDUR_MOVIE, MOVIE_TOP_KM, MOVIE_BOTTOM_KM, &
        MOVIE_EAST_DEG, MOVIE_WEST_DEG, MOVIE_NORTH_DEG, MOVIE_SOUTH_DEG, &
        RECORD_LENGTH_IN_MINUTES, NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
        NTSTEP_BETWEEN_READ_ADJSRC, NTSTEP_BETWEEN_FRAMES, &
        NTSTEP_BETWEEN_OUTPUT_INFO, NUMBER_OF_RUNS, NUMBER_OF_THIS_RUN,NCHUNKS,&
        SIMULATION_TYPE, MOVIE_VOLUME_TYPE, MOVIE_START, MOVIE_STOP, NEX_XI, &
        NEX_ETA, NPROC_XI, NPROC_ETA, NOISE_TOMOGRAPHY, ELLIPTICITY, GRAVITY, &
        ROTATION, TOPOGRAPHY, OCEANS, MOVIE_SURFACE, MOVIE_VOLUME,MOVIE_COARSE,&
        RECEIVERS_CAN_BE_BURIED, PRINT_SOURCE_TIME_FUNCTION, SAVE_MESH_FILES, &
        ATTENUATION, ABSORBING_CONDITIONS, SAVE_FORWARD, &
        OUTPUT_SEISMOS_ASCII_TEXT, OUTPUT_SEISMOS_SAC_ALPHANUM, &
        OUTPUT_SEISMOS_SAC_BINARY, ROTATE_SEISMOGRAMS_RT, &
        WRITE_SEISMOGRAMS_BY_MASTER, SAVE_ALL_SEISMOS_IN_ONE_FILE, &
        USE_BINARY_FOR_LARGE_FILE, model_length, MODEL)

    ! Write variables from 'CMTSOLUTION'
    call write_adios_cmtsolution_variables (adios_handle, &
        NSOURCES, yr, mo, da, ho, mi, sec, t_shift, hdur, lat, long, depth, &
        mrr, mtt, mpp, mrt, mrp, mtp, event_name_length, event_name, &
        datasource_length, datasource)

    ! Write variables from 'STATIONS'
    call write_adios_stations_variables (adios_handle, &
       NSTATIONS, stlat, stlon, stele, stbur, station_name_length, &
       station_name, network_name_length, network_name)

    call adios_close (adios_handle, adios_err)

    deallocate(datasource)
    deallocate(station_name)
    deallocate(network_name)
    deallocate(stlat)
    deallocate(stlon)
    deallocate(stele)
    deallocate(stbur)
  endif

! Imbricated/contained subroutines. The initial thought was to do a module with
! public access to the write_specfem_header_adios routine and private access to
! the other routines. The problem then is the files compilation order that
! should be done very carefully. This require modifications of the Makefile
! which is not currently designed to do that.
contains

!> \brief Define ADIOS variable to store values from 'setup/config.h'. Store
!!        configuration parameters to insure reproducibility
!! \param adios_group The ADIOS entity grouping variables for data transfers
!! \param group_size_inc The group size to increment wrt. the variable size
subroutine define_solver_info_variables (adios_group, group_size_inc)

  implicit none
  ! Parameters
  integer(kind=8), intent(in)    :: adios_group
  integer(kind=8), intent(inout) :: group_size_inc
  ! Variables
  integer :: pkg_str_len, conf_flags_len

  pkg_str_len    = len_trim(PACKAGE_STRING)
  conf_flags_len = len_trim(CONFIGURE_FLAGS)

  call define_adios_integer_scalar_s (adios_group, "package_string_length", &
      "/solver_info", group_size_inc)
  call define_adios_string_s (adios_group, "package_name", "/solver_info", &
      pkg_str_len, group_size_inc)

  call define_adios_integer_scalar_s (adios_group, "conf_flags_len", &
      "/solver_info", group_size_inc)
  call define_adios_string_s (adios_group, "conf_flags", "/solver_info", &
      conf_flags_len, group_size_inc)

end subroutine define_solver_info_variables

!> \brief Define ADIOS variable to store values from the Par_file
!! \param adios_group The ADIOS entity grouping variables for data transfers
!! \param group_size_inc The group size to increment wrt. the variable size
!! \param model_length The number of character of the MODEL string.
!!                     Useful for reading back the MODEL
subroutine define_par_file_variables (adios_group, group_size_inc, model_length)

  implicit none

  ! Parameters
  integer(kind=8), intent(in)    :: adios_group
  integer(kind=8), intent(inout) :: group_size_inc
  integer, intent(in)            :: model_length ! for later reading of MODEL

  !-- double precision variables
  call define_adios_double_scalar_s (adios_group, "ANGULAR_WIDTH_XI_IN_DEGREES", &
      "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_double_scalar_s (adios_group, "ANGULAR_WIDTH_ETA_IN_DEGREES", &
      "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_double_scalar_s (adios_group, "CENTER_LONGITUDE_IN_DEGREES", &
      "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_double_scalar_s (adios_group, "CENTER_LATITUDE_IN_DEGREES", &
      "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_double_scalar_s (adios_group, "GAMMA_ROTATION_AZIMUTH", &
      "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_double_scalar_s (adios_group, "HDUR_MOVIE", &
      "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_double_scalar_s (adios_group, "MOVIE_TOP_KM", &
      "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_double_scalar_s (adios_group, "MOVIE_BOTTOM_KM", &
      "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_double_scalar_s (adios_group, "MOVIE_EAST_DEG", &
      "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_double_scalar_s (adios_group, "MOVIE_WEST_DEG", &
      "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_double_scalar_s (adios_group, "MOVIE_NORTH_DEG", &
      "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_double_scalar_s (adios_group, "MOVIE_SOUTH_DEG", &
      "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_double_scalar_s (adios_group, "RECORD_LENGTH_IN_MINUTES", &
      "/specfem3D_globe_parameter_file", group_size_inc)
  !-- integer variables
  call define_adios_integer_scalar_s (adios_group, "NTSTEP_BETWEEN_OUTPUT_SEISMOS", &
      "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_integer_scalar_s(adios_group, "NTSTEP_BETWEEN_READ_ADJSRC", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_integer_scalar_s(adios_group, "NTSTEP_BETWEEN_FRAMES", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_integer_scalar_s(adios_group, "NTSTEP_BETWEEN_OUTPUT_INFO", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_integer_scalar_s(adios_group, "NUMBER_OF_RUNS", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_integer_scalar_s(adios_group, "NUMBER_OF_THIS_RUN", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_integer_scalar_s(adios_group, "NCHUNKS", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_integer_scalar_s(adios_group, "SIMULATION_TYPE", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_integer_scalar_s(adios_group, "MOVIE_VOLUME_TYPE", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_integer_scalar_s(adios_group, "MOVIE_START", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_integer_scalar_s(adios_group, "MOVIE_STOP", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_integer_scalar_s(adios_group, "NEX_XI", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_integer_scalar_s(adios_group, "NEX_ETA", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_integer_scalar_s(adios_group, "NPROC_XI", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_integer_scalar_s(adios_group, "NPROC_ETA", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_integer_scalar_s(adios_group, "NOISE_TOMOGRAPHY", "/specfem3D_globe_parameter_file", group_size_inc)
  !-- logical variables
  call define_adios_byte_scalar_s(adios_group, "ELLIPTICITY", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "GRAVITY", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "ROTATION", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "TOPOGRAPHY", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "OCEANS", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "MOVIE_SURFACE", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "MOVIE_VOLUME", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "MOVIE_COARSE", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "RECEIVERS_CAN_BE_BURIED", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "PRINT_SOURCE_TIME_FUNCTION", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "SAVE_MESH_FILES", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "ATTENUATION", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "ABSORBING_CONDITIONS", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "SAVE_FORWARD", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "OUTPUT_SEISMOS_ASCII_TEXT", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "OUTPUT_SEISMOS_SAC_ALPHANUM", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "OUTPUT_SEISMOS_SAC_BINARY", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "ROTATE_SEISMOGRAMS_RT", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "WRITE_SEISMOGRAMS_BY_MASTER", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "SAVE_ALL_SEISMOS_IN_ONE_FILE", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_byte_scalar_s(adios_group, "USE_BINARY_FOR_LARGE_FILE", "/specfem3D_globe_parameter_file", group_size_inc)
  !-- string variables
  call define_adios_integer_scalar_s (adios_group, "model_length", "/specfem3D_globe_parameter_file", group_size_inc)
  call define_adios_string_s (adios_group, "MODEL", "/specfem3D_globe_parameter_file", model_length, group_size_inc)

end subroutine define_par_file_variables


!> \brief Define ADIOS variable to store values from the CMTSOLUTION file
!! \param adios_group The ADIOS entity grouping variables for data transfers
!! \param group_size_inc The group size to increment wrt. the variable size
!! \param NSOURCES The number of sources. Needed to define array sizes.
!! \param datasource_length The number of character of the datasource string.
!!                          Useful for reading back the datasources.
subroutine define_cmtsolution_variables (adios_group, group_size_inc, NSOURCES,&
                                         event_name_length, datasource_length)

  implicit none
  integer(kind=8), intent(in)    :: adios_group
  integer(kind=8), intent(inout) :: group_size_inc
  integer, intent(in) :: NSOURCES, datasource_length, event_name_length

  !-- Number of SOURCES inside the CMTSOLUTION file
  call define_adios_integer_scalar_s (adios_group, "NSOURCES", "/CMTSOLUTION", group_size_inc)
  !-- double precision arrays
  call define_adios_double_local_array1D_s (adios_group, "second", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
  call define_adios_double_local_array1D_s (adios_group, "time_shift", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
  call define_adios_double_local_array1D_s (adios_group, "half_duration", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
  call define_adios_double_local_array1D_s (adios_group, "latitude", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
  call define_adios_double_local_array1D_s (adios_group, "longitude", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
  call define_adios_double_local_array1D_s (adios_group, "depth", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
  call define_adios_double_local_array1D_s (adios_group, "mrr", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
  call define_adios_double_local_array1D_s (adios_group, "mtt", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
  call define_adios_double_local_array1D_s (adios_group, "mpp", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
  call define_adios_double_local_array1D_s (adios_group, "mrt", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
  call define_adios_double_local_array1D_s (adios_group, "mrp", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
  call define_adios_double_local_array1D_s (adios_group, "mtp", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
  !-- integer arrays
  call define_adios_integer_local_array1D_s (adios_group, "year", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
  call define_adios_integer_local_array1D_s (adios_group, "month", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
  call define_adios_integer_local_array1D_s (adios_group, "day", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
  call define_adios_integer_local_array1D_s (adios_group, "hour", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
  call define_adios_integer_local_array1D_s (adios_group, "minute", "/CMTSOLUTION", NSOURCES, "NSOURCES", group_size_inc)
  !-- string
  call define_adios_integer_scalar_s (adios_group, "event_name_length", "/CMTSOLUTION", group_size_inc)
  call define_adios_string_s (adios_group, "event_name", "/CMTSOLUTION", event_name_length, group_size_inc)
  call define_adios_integer_scalar_s (adios_group, "datasource_length", "/CMTSOLUTION", group_size_inc)
  call define_adios_string_s (adios_group, "datasource", "/CMTSOLUTION", datasource_length, group_size_inc)

end subroutine define_cmtsolution_variables

!> \brief Define ADIOS variable to store values from the STATIONS file
!! \param adios_group The ADIOS entity grouping variables for data transfers
!! \param group_size_inc The group size to increment wrt. the variable size
!! \param NSTATIONS The number of stations. Needed to define array sizes.
!! \param station_name_length The number of character of the station_name
!!                            string.  Useful for reading back the stations.
!! \param network_name_length The number of character of the station_name
!!                            string.  Useful for reading back the networks.
subroutine define_stations_variables (adios_group, group_size_inc, NSTATIONS,&
                                      station_name_length, network_name_length)

  implicit none
  integer(kind=8), intent(in)    :: adios_group
  integer(kind=8), intent(inout) :: group_size_inc
  integer, intent(in) :: NSTATIONS, station_name_length, network_name_length

  !-- Number of STATIONS inside the STATIONS file
  call define_adios_integer_scalar_s (adios_group, "NSTATIONS", "/STATIONS", group_size_inc)
  !-- double precision arrays
  call define_adios_double_local_array1D_s (adios_group, "station_latitude", "/STATIONS", NSTATIONS, "NSTATIONS", group_size_inc)
  call define_adios_double_local_array1D_s (adios_group, "station_longitude", "/STATIONS", NSTATIONS, "NSTATIONS", group_size_inc)
  call define_adios_double_local_array1D_s (adios_group, "station_elevation", "/STATIONS", NSTATIONS, "NSTATIONS", group_size_inc)
  call define_adios_double_local_array1D_s (adios_group, "station_burial", "/STATIONS", NSTATIONS, "NSTATIONS", group_size_inc)
  !-- string
  call define_adios_integer_scalar_s (adios_group, "station_name_length", "/STATIONS", group_size_inc)
  call define_adios_integer_scalar_s (adios_group, "network_name_length", "/STATIONS", group_size_inc)
  call define_adios_string_s (adios_group, "station_name", "/STATIONS", station_name_length, group_size_inc)
  call define_adios_string_s (adios_group, "network_name", "/STATIONS", network_name_length, group_size_inc)

end subroutine define_stations_variables

!> \brief Read the 'CMTSOLUTION file' and do not modify nor transform variables
!! \param yr Array to store the year of the events
!! \param mo Array to store the month of the events
!! \param da Array to store the day of the events
!! \param ho Array to store the hour of the events
!! \param mi Array to store the minute of the events
!! \param sec Array to store the second of the events
!! \param t_shift Array to store the time shift at the beginning of the events
!! \param hdur Array to store the duration of the events
!! \param lat Array to store the latitude of the events
!! \param long Array to store the longitude of the events
!! \param depth Arrays to store the depth of the events
!! \param mrr Arrays to store the mrr component of the events
!! \param mtt Arrays to store the mtt component of the events
!! \param mpp Arrays to store the mpp component of the events
!! \param mrt Arrays to store the mrt component of the events
!! \param mrp Arrays to store the mrp component of the events
!! \param mtp Arrays to store the mtp component of the events
!! \param event_name_length Variable for keeping the size of the event_name
!!                          string
!! \param event_name Strings to store the event name
!! \param  datasource_length Variable for keeping the size of the datasource
!!                          string
!! \param datasource String in which the different datasource names are
!!                   concatenated
!> \note This subroutine and get_cmt.f90 are redundant. Might be factorized in
!!       the future. For now we do not want the value modification from get_cmt
subroutine read_raw_cmtsolution (yr, mo, da, ho, mi, sec, t_shift, hdur, lat, &
                                 long, depth, mrr, mtt, mpp, mrt, mrp, mtp, event_name_length, event_name, &
                                 datasource_length, datasource)

  use constants,only: IIN
  implicit none
  ! Parameters
  integer,           dimension(NSOURCES), intent(out) :: yr, mo, da, ho, mi
  double precision,  dimension(NSOURCES), intent(out) :: sec, t_shift, hdur, lat, long, depth
  double precision,  dimension(NSOURCES), intent(out) :: mrr, mtt, mpp, mrt, mrp, mtp
  integer, intent(inout) :: event_name_length, datasource_length
  character(len=16), intent(out) :: event_name
  character(len=:), allocatable, intent(out)  :: datasource  ! F03 feature
  ! Local variables
  character(len=5)   :: datasource_tmp
  character(len=256) :: string
  ! extract all unmodified values from CMTSOLUTION
  ! get_cmt() routine modify the read values
  ! TODO factorize what follows and get_cmt.f90 and probably one or two other
  !      routines
  open(unit = IIN,file='DATA/CMTSOLUTION',status='old',action='read',iostat=ier)
  if (ier /= 0) stop 'Error opening file DATA/CMTSOLUTION'

  datasource_length = 4*NSOURCES ! a datasource is 4 character, by convention

  allocate(character(len=(datasource_length)) :: datasource, stat=ier)
  if (ier /= 0) call exit_MPI (myrank, "Error allocating datasource string for adios header")

  datasource = ""
  ! ADIOS only  (1) byte for a string. This may cause data overwriting.
  ! => increase the generate by the string size -1
  adios_groupsize = adios_groupsize + 4*NSOURCES - 1
  do isource = 1,NSOURCES

    read(IIN,"(a256)") string
    ! skips empty lines
    do while( len_trim(string) == 0 )
    read(IIN,"(a256)") string
    enddo
    ! read header with event information
    read(string,"(a4,i5,i3,i3,i3,i3,f6.2)") datasource_tmp,yr(isource), &
        mo(isource),da(isource),ho(isource),mi(isource),sec(isource)
    datasource = datasource // datasource_tmp
    ! read event name
    read(IIN,"(a)") string
    read(string(12:len_trim(string)),*) event_name
    ! read time shift
    read(IIN,"(a)") string
    read(string(12:len_trim(string)),*) t_shift(isource)
    ! read half duration
    read(IIN,"(a)") string
    read(string(15:len_trim(string)),*) hdur(isource)
    ! read latitude
    read(IIN,"(a)") string
    read(string(10:len_trim(string)),*) lat(isource)
    ! read longitude
    read(IIN,"(a)") string
    read(string(11:len_trim(string)),*) long(isource)
    ! read depth
    read(IIN,"(a)") string
    read(string(7:len_trim(string)),*) depth(isource)
    ! read Mrr
    read(IIN,"(a)") string
    read(string(5:len_trim(string)),*) mrr(isource)
    ! read Mtt
    read(IIN,"(a)") string
    read(string(5:len_trim(string)),*) mtt(isource)
    ! read Mpp
    read(IIN,"(a)") string
    read(string(5:len_trim(string)),*) mpp(isource)
    ! read Mrt
    read(IIN,"(a)") string
    read(string(5:len_trim(string)),*) mrt(isource)
    ! read Mrp
    read(IIN,"(a)") string
    read(string(5:len_trim(string)),*) mrp(isource)
    ! read Mtp
    read(IIN,"(a)") string
    read(string(5:len_trim(string)),*) mtp(isource)
  enddo
  close(IIN)
  event_name_length = len_trim(event_name)

end subroutine read_raw_cmtsolution

!> \brief Reads information form the 'STATIONS' file without modifying anything
!! \param NSTATIONS How many stations are used
!! \param stlat Array to store the latitude of the stations
!! \param stlon Array to store the longitude of the stations
!! \param stele Array to store the elevation of the stations
!! \param stbur Array to store the burial of the stations
!! \param station_name_length Variable to keep the length of the station_name
!!                            string
!! \param station_name  String in which the different station names are
!!                      concatenated
!! \param network_name_length Variable to keep the length of the network_name
!!                            string
!! \param network_name String in which the different network names are
!!                     concatenated
subroutine read_raw_stations (NSTATIONS, stlat, stlon, stele, stbur, &
                              station_name_length, station_name, network_name_length, network_name)

  use constants,only: MAX_LENGTH_STATION_NAME,MAX_LENGTH_NETWORK_NAME,IMAIN,IIN

  implicit none
  ! Parameters
  integer :: NSTATIONS
  integer, intent(inout) :: station_name_length, network_name_length ! for later reading
  character(len=:), allocatable, intent(out) :: station_name, network_name
  double precision, allocatable, dimension(:), intent(out) :: stlat, stlon, stele, stbur
  ! Local variables
  character(len=MAX_LENGTH_STATION_NAME) :: station_name_tmp
  character(len=MAX_LENGTH_NETWORK_NAME) :: network_name_tmp
  character(len=256) :: string

  ! Extract values from STATIONS File
  open(unit = IIN,file='DATA/STATIONS',status='old',action='read',iostat=ier)
  if (ier /= 0) stop 'Error opening file DATA/STATIONS'

  NSTATIONS = 0
  do while(ier == 0)
    read(IIN,"(a)",iostat=ier) string
    if (ier == 0) NSTATIONS = NSTATIONS + 1
  enddo

  allocate (character (len=(MAX_LENGTH_STATION_NAME*NSTATIONS)) :: station_name)
  allocate (character (len=(MAX_LENGTH_NETWORK_NAME*NSTATIONS)) :: network_name)
  allocate (stlat (NSTATIONS))
  allocate (stlon (NSTATIONS))
  allocate (stele (NSTATIONS))
  allocate (stbur (NSTATIONS))
  station_name = ""
  network_name = ""

  rewind(IIN)
  do irec = 1,NSTATIONS
    read(IIN,*,iostat=ier) station_name_tmp, network_name_tmp, &
    stlat(irec),      stlon(irec), &
    stele(irec),      stbur(irec)
    if (ier /= 0) then
      write(IMAIN,*) 'Error reading in station ',irec
      call exit_MPI(myrank,'Error reading in station in DATA/STATIONS file')
    endif
    station_name = station_name // trim(station_name_tmp) // " "
    network_name = network_name // trim(network_name_tmp) // " "
  enddo
  close(IIN)

  station_name = trim(station_name)
  network_name = trim(network_name)
  station_name_length = len(station_name)
  network_name_length = len(network_name)

end subroutine read_raw_stations

!> \brief Wrapper to write the 'config.h' variables into the adios header
!! \param adios_handle The handle to the file where the variable should be
!!                     written
subroutine write_adios_solver_info_variables (adios_handle)

  implicit none
  ! Parameters
  integer(kind=8), intent(in)    :: adios_handle
  ! Variables
  integer :: pkg_str_len, conf_flags_len, adios_err
  character(len=len_trim(PACKAGE_STRING)) :: pkg_str
  character(len=len_trim(CONFIGURE_FLAGS)) :: conf_flags

  pkg_str    = trim(PACKAGE_STRING)
  conf_flags = trim(CONFIGURE_FLAGS)

  pkg_str_len    = len_trim(PACKAGE_STRING)
  conf_flags_len = len_trim(CONFIGURE_FLAGS)

  call adios_write (adios_handle, "package_string_length", pkg_str_len, adios_err)
  call adios_write (adios_handle, "package_name", pkg_str, adios_err)
  call adios_write (adios_handle, "conf_flags_len", conf_flags_len, adios_err)
  call adios_write (adios_handle, "conf_flags", conf_flags, adios_err)

end subroutine write_adios_solver_info_variables

!> \brief Wrapper to write the 'Par_file' variables into the adios header
!! \param adios_handle The handle to the file where the variable should be
!!                     written
subroutine write_adios_par_file_variables (adios_handle, &
                                           ANGULAR_WIDTH_XI_IN_DEGREES, ANGULAR_WIDTH_ETA_IN_DEGREES, &
                                           CENTER_LONGITUDE_IN_DEGREES, CENTER_LATITUDE_IN_DEGREES, &
                                           GAMMA_ROTATION_AZIMUTH, HDUR_MOVIE, MOVIE_TOP_KM, MOVIE_BOTTOM_KM, &
                                           MOVIE_EAST_DEG, MOVIE_WEST_DEG, MOVIE_NORTH_DEG, MOVIE_SOUTH_DEG, &
                                           RECORD_LENGTH_IN_MINUTES, NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
                                           NTSTEP_BETWEEN_READ_ADJSRC, NTSTEP_BETWEEN_FRAMES, &
                                           NTSTEP_BETWEEN_OUTPUT_INFO, NUMBER_OF_RUNS, NUMBER_OF_THIS_RUN, NCHUNKS, &
                                           SIMULATION_TYPE, MOVIE_VOLUME_TYPE, MOVIE_START, MOVIE_STOP, NEX_XI, &
                                           NEX_ETA, NPROC_XI, NPROC_ETA, NOISE_TOMOGRAPHY, ELLIPTICITY, GRAVITY, &
                                           ROTATION, TOPOGRAPHY, OCEANS, MOVIE_SURFACE, MOVIE_VOLUME, MOVIE_COARSE, &
                                           RECEIVERS_CAN_BE_BURIED, PRINT_SOURCE_TIME_FUNCTION, SAVE_MESH_FILES, &
                                           ATTENUATION, ABSORBING_CONDITIONS, SAVE_FORWARD, &
                                           OUTPUT_SEISMOS_ASCII_TEXT, OUTPUT_SEISMOS_SAC_ALPHANUM, &
                                           OUTPUT_SEISMOS_SAC_BINARY, ROTATE_SEISMOGRAMS_RT, &
                                           WRITE_SEISMOGRAMS_BY_MASTER, SAVE_ALL_SEISMOS_IN_ONE_FILE, &
                                           USE_BINARY_FOR_LARGE_FILE, model_length, MODEL)

 implicit none
 ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in)  :: NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
      NTSTEP_BETWEEN_READ_ADJSRC, NTSTEP_BETWEEN_FRAMES, &
      NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS, NUMBER_OF_THIS_RUN,NCHUNKS, &
      SIMULATION_TYPE, MOVIE_VOLUME_TYPE, MOVIE_START,MOVIE_STOP, NEX_XI, &
      NEX_ETA,NPROC_XI,NPROC_ETA, NOISE_TOMOGRAPHY
  double precision, intent(in) :: ANGULAR_WIDTH_XI_IN_DEGREES, &
      ANGULAR_WIDTH_ETA_IN_DEGREES, CENTER_LONGITUDE_IN_DEGREES, &
      CENTER_LATITUDE_IN_DEGREES, GAMMA_ROTATION_AZIMUTH, HDUR_MOVIE, &
      MOVIE_TOP_KM,MOVIE_BOTTOM_KM, MOVIE_EAST_DEG,MOVIE_WEST_DEG, &
      MOVIE_NORTH_DEG,MOVIE_SOUTH_DEG, RECORD_LENGTH_IN_MINUTES
  logical, intent(in) :: ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS, &
      MOVIE_SURFACE, MOVIE_VOLUME,MOVIE_COARSE, RECEIVERS_CAN_BE_BURIED, &
      PRINT_SOURCE_TIME_FUNCTION, SAVE_MESH_FILES,ATTENUATION, &
      ABSORBING_CONDITIONS,SAVE_FORWARD, OUTPUT_SEISMOS_ASCII_TEXT, &
      OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
      ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER, &
      SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE
  integer, intent(in) :: model_length
  character(len=*), intent(in) :: MODEL
  ! Local variables
  integer  :: adios_err

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
  call adios_write (adios_handle, "ABSORBING_CONDITIONS", ABSORBING_CONDITIONS, adios_err)
  call adios_write (adios_handle, "SAVE_FORWARD", SAVE_FORWARD, adios_err)
  call adios_write (adios_handle, "OUTPUT_SEISMOS_ASCII_TEXT", OUTPUT_SEISMOS_ASCII_TEXT, adios_err)
  call adios_write (adios_handle, "OUTPUT_SEISMOS_SAC_ALPHANUM", OUTPUT_SEISMOS_SAC_ALPHANUM, adios_err)
  call adios_write (adios_handle, "OUTPUT_SEISMOS_SAC_BINARY", OUTPUT_SEISMOS_SAC_BINARY, adios_err)
  call adios_write (adios_handle, "ROTATE_SEISMOGRAMS_RT", ROTATE_SEISMOGRAMS_RT, adios_err)
  call adios_write (adios_handle, "WRITE_SEISMOGRAMS_BY_MASTER", WRITE_SEISMOGRAMS_BY_MASTER, adios_err)
  call adios_write (adios_handle, "SAVE_ALL_SEISMOS_IN_ONE_FILE", SAVE_ALL_SEISMOS_IN_ONE_FILE, adios_err)
  call adios_write (adios_handle, "USE_BINARY_FOR_LARGE_FILE", USE_BINARY_FOR_LARGE_FILE, adios_err)
  call adios_write (adios_handle, "model_length", model_length, adios_err)
  call adios_write (adios_handle, "MODEL", MODEL, adios_err)

end subroutine write_adios_par_file_variables

!> \brief Wrapper to write the 'CMTSOLUTION' variables into the adios header
!! \param adios_handle The handle to the file where the variable should be
!!                     written
subroutine write_adios_cmtsolution_variables (adios_handle, &
                                              NSOURCES, yr, mo, da, ho, mi, sec, t_shift, hdur, lat, long, depth, &
                                              mrr, mtt, mpp, mrt, mrp, mtp, event_name_length, event_name, &
                                              datasource_length, datasource)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: NSOURCES
  integer,          dimension(NSOURCES), intent(in) :: yr, mo, da, ho, mi
  double precision, dimension(NSOURCES), intent(in) :: sec, t_shift, hdur, &
                                                       lat, long, depth
  double precision, dimension(NSOURCES), intent(in) :: mrr, mtt, mpp, &
                                                       mrt, mrp, mtp
  integer, intent(in) :: event_name_length, datasource_length
  character(len=16), intent(in) :: event_name
  character(len=:), allocatable, intent(in)  :: datasource  ! F03 feature
  ! Local variables
  integer :: adios_err

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
  call adios_write (adios_handle, "event_name_length", event_name_length, adios_err)
  call adios_write (adios_handle, "event_name", event_name, adios_err)
  call adios_write (adios_handle, "datasource_length", datasource_length, adios_err)
  call adios_write (adios_handle, "datasource", datasource, adios_err)

end subroutine write_adios_cmtsolution_variables

!> \brief Wrapper to write the 'STATIONS' variables into the adios header
!! \param adios_handle The handle to the file where the variable should be
!!                     written
subroutine write_adios_stations_variables (adios_handle, &
                                           NSTATIONS, stlat, stlon, stele, stbur, station_name_length, station_name, &
                                           network_name_length, network_name)

  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in):: NSTATIONS
  integer, intent(in):: station_name_length, network_name_length ! for later reading
  character(len=:), allocatable, intent(in) :: station_name, network_name
  double precision, allocatable, dimension(:), intent(in) :: stlat, stlon, &
                                                             stele, stbur
  ! Local variables
  integer :: adios_err

  call adios_write (adios_handle, "NSTATIONS", NSTATIONS, adios_err)
  call adios_write (adios_handle, "station_latitude", stlat, adios_err)
  call adios_write (adios_handle, "station_longitude", stlon, adios_err)
  call adios_write (adios_handle, "station_elevation", stele, adios_err)
  call adios_write (adios_handle, "station_burial", stbur, adios_err)
  call adios_write (adios_handle, "station_name_length", station_name_length, adios_err)
  call adios_write (adios_handle, "network_name_length", network_name_length, adios_err)
  call adios_write (adios_handle, "station_name", station_name, adios_err)
  call adios_write (adios_handle, "network_name", network_name, adios_err)

end subroutine write_adios_stations_variables


!===============================================================================
!> Define an ADIOS scalar double precision variable and autoincrement
!! the adios group size by (8).
!! \param adios_group The adios group where the variables belongs
!! \param name The variable to be defined
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
subroutine define_adios_double_scalar_s (adios_group, name, path, group_size_inc)

  use adios_write_mod
  implicit none
  ! Arguments
  integer(kind=8),  intent(in)     :: adios_group
  character(len=*), intent(in)     :: name, path
  integer(kind=8),  intent(inout)  :: group_size_inc
  ! Local Variables
  integer(kind=8)                  :: varid ! dummy variable, adios use var name

  ! adios: 6 == real(kind=8)
  call adios_define_var (adios_group, name, path, 6,  "", "", "", varid)
  group_size_inc = group_size_inc + 8

end subroutine define_adios_double_scalar_s

!===============================================================================
!> Define an ADIOS scalar integer variable and autoincrement the adios
!! group size by (4).
!! \param adios_group The adios group where the variables belongs
!! \param name The variable to be defined
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
subroutine define_adios_integer_scalar_s (adios_group, name, path, group_size_inc)

  use adios_write_mod
  implicit none
  ! Arguments
  integer(kind=8),  intent(in)     :: adios_group
  character(len=*), intent(in)     :: name, path
  integer(kind=8),  intent(inout)  :: group_size_inc
  ! Local Variables
  integer(kind=8)                  :: varid ! dummy variable, adios use var name

  ! adios: 2 == integer(kind=4)
  call adios_define_var (adios_group, name, path, adios_integer,  "", "", "", varid)
  group_size_inc = group_size_inc + 4

end subroutine define_adios_integer_scalar_s

!===============================================================================
!> Define an ADIOS scalar byte variable and autoincrement the adios
!! group size by (1).
!! \param adios_group The adios group where the variables belongs
!! \param name The variable to be defined
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
subroutine define_adios_byte_scalar_s (adios_group, name, path, group_size_inc)

  use adios_write_mod
  implicit none
  ! Arguments
  integer(kind=8),  intent(in)     :: adios_group
  character(len=*), intent(in)     :: name, path
  integer(kind=8),  intent(inout)  :: group_size_inc
  ! Local Variables
  integer(kind=8)                  :: varid ! dummy variable, adios use var name

  ! adios: 0 == byte == any_data_type(kind=1)
  call adios_define_var (adios_group, name, path, 0,  "", "", "", varid)
  group_size_inc = group_size_inc + 1

end subroutine define_adios_byte_scalar_s

!===============================================================================
!> Define a local ADIOS string and autoincrement the adios
!! group size by (1 * string's length).
!! \param adios_group The adios group where the variables belongs
!! \param name The variable to be defined
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param len The length of the string(number of character. in Fortran
!!            it does not include a final '\0' -- null -- character)
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \note Adios string are scalar values counting for (1) byte. It is
!!       mandatory to increase the group size by the length of the
!!       string in order not to overlap 'data regions'.
subroutine define_adios_string_s (adios_group, name, path, length, group_size_inc)

  use adios_write_mod
  implicit none
  ! Arguments
  integer(kind=8),  intent(in)     :: adios_group
  character(len=*), intent(in)     :: name, path
  integer(kind=8),  intent(inout)  :: group_size_inc
  integer                          :: length
  ! Local Variables
  integer(kind=8)                  :: varid ! dummy variable, adios use var name

  ! adios: 9 == string
  call adios_define_var (adios_group, name, path, 9,  "", "", "", varid)
  group_size_inc = group_size_inc + 1*length

end subroutine define_adios_string_s

!===============================================================================
!> Define a local ADIOS array of doubles and autoincrement the adios
!! group size by (8 * number of elements).
!! \param adios_group The adios group where the variables belongs
!! \param name The variable to be defined
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param dim The number of elements in the 1D array. Required to
!!            correctly increment adios group size.
!! \param dim_str The "stringified" version of dim. Needed by adios
!!                to define variables
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
subroutine define_adios_double_local_array1D_s (adios_group, name, path, dim, dim_str, group_size_inc)

  use adios_write_mod
  implicit none
  ! Arguments
  integer(kind=8),  intent(in)     :: adios_group
  character(len=*), intent(in)     :: name, path, dim_str
  integer(kind=8),  intent(inout)  :: group_size_inc
  integer, intent(in)              :: dim
  ! Local Variables
  integer(kind=8)                  :: varid ! dummy variable, adios use var name

  ! adios: 6 == real(kind=8)
  call adios_define_var (adios_group, name, path, 6, dim_str, "", "", varid)
  group_size_inc = group_size_inc + 8*dim

end subroutine define_adios_double_local_array1D_s

!===============================================================================
!> Define a local ADIOS array of integers and autoincrement the adios
!! group size by (4 * number of elements).
!! \param adios_group The adios group where the variables belongs
!! \param name The variable to be defined
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param dim The number of elements in the 1D array. Required to
!!            correctly increment adios group size.
!! \param dim_str The "stringified" version of dim. Needed by adios
!!                to define variables
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
subroutine define_adios_integer_local_array1D_s (adios_group, name, path, dim, dim_str, group_size_inc)

  use adios_write_mod
  implicit none
  ! Arguments
  integer(kind=8),  intent(in)     :: adios_group
  character(len=*), intent(in)     :: name, path, dim_str
  integer(kind=8),  intent(inout)  :: group_size_inc
  integer, intent(in)              :: dim
  ! Local Variables
  integer(kind=8)                  :: varid ! dummy variable, adios use var name

  ! adios: 2 == integer
  call adios_define_var (adios_group, name, path, 2,  dim_str, "", "", varid)
  group_size_inc = group_size_inc + 4*dim

end subroutine define_adios_integer_local_array1D_s


end subroutine write_specfem_header_adios
