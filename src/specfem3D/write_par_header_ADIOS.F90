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

!> @brief Write solver name and version into ADIOS result file header.
!! The solver name, version name, version number are written down inside
!! the xml file.
!subroutine write_solver_info_header_ADIOS()
!  use mpi
!  use adios_write_mod
!  use specfem_par, only : myrank
  !-------------------------------------------------------------------
  ! local parameters
  !-------------------------------------------------------------------
!  integer                 :: comm, ier, NUMBER_MPI_PROCESSES
  ! Adios variables
!  integer                 :: adios_err
!  integer*8               :: adios_handle
!  integer*8               :: adios_groupsize, adios_totalsize ! for the .fh file
  ! TODO : find a better name once the use of ADIOS is more completely
  ! implemented
!  character(len=27)       :: filename = "OUTPUT_FILES/header_tmp.bp" 

!    call MPI_COMM_DUP(MPI_COMM_WORLD, comm, ier)
!    call world_size(NUMBER_MPI_PROCESSES)

!  if(myrank == 0) then
!    call adios_open (adios_handle, "solver_info", filename, "w", &
!                     MPI_COMM_SELF, adios_err);
!                     comm, adios_err);

    ! Automatically generated Fortran code by ADIOS gpp.py script:
    ! $> $ADIOS_DIR/bin/gpp.py par_header.xml
    ! Re-run each time the xml file is modified.
!#include "gwrite_solver_info.fh"

!    call adios_close (adios_handle, adios_err)
!  endif

!end subroutine write_solver_info_header_ADIOS

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
  integer :: datasources_length
  character(len=5) :: datasource_tmp
  character(len=16):: event_name
  character(len=:), allocatable  :: datasource
  !character(len=256) :: datasource

  ! values from STATIONS
  integer :: NSTATIONS
  character(len=MAX_LENGTH_STATION_NAME) :: station_name_tmp
  character(len=MAX_LENGTH_NETWORK_NAME) :: network_name_tmp
  character(len=:),   allocatable :: station_name, network_name 
  double precision, allocatable, dimension(:) :: stlat, stlon, stele, stbur

  character(len=150) :: OUTPUT_FILES,LOCAL_PATH,LOCAL_TMP_PATH,MODEL
  character(len=256) :: string, CMTSOLUTION, STATIONS

  integer                 :: isource, irec, ier
  ! Adios variables
  integer                 :: adios_err
  integer*8               :: adios_handle
  integer*8               :: adios_groupsize, adios_totalsize ! for the .fh file
  ! TODO : find a better name once the use of ADIOS is more completely
  ! implemented
  character(len=27)       :: filename = "OUTPUT_FILES/header_tmp.bp" 

  ! ensure that only the master open the adios handle inside MPI_COMM_SELF
  if(myrank == 0) then
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
    ! extract all unmodified values from CMTSOLUTION
    ! get_cmt() routine modify the read values
    ! TODO factorize what follows and get_cmt.f90 and probably one or two other
    !      routines
    call get_value_string(CMTSOLUTION, 'solver.CMTSOLUTION', 'DATA/CMTSOLUTION')
    open(unit=1,file=CMTSOLUTION,status='old',action='read')
    datasources_length = 4*NSOURCES
    allocate(character(len=(datasources_length)) :: datasource, stat=ier)
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

    call adios_open (adios_handle, "SPECFEM3D_GLOBE_HEADER", filename, "w", &
                     MPI_COMM_SELF, adios_err);

    ! Automatically generated Fortran code by ADIOS gpp.py script:
    ! $> $ADIOS_DIR/bin/gpp.py par_header.xml
    ! Re-run each time the xml file is modified.
#include "gwrite_SPECFEM3D_GLOBE_HEADER.fh"

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
