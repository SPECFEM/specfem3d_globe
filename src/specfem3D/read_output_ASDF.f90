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

!==============================================================================
! Copyright 2015 ASDF developers
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!==============================================================================

!
!-------------------------------------------------------------------------------------------------
!

!> Read the ASDF file
subroutine read_asdf()

  use iso_c_binding
  use specfem_par

  implicit none

  integer, parameter :: MAX_STRING_LENGTH = 256

  character(len=MAX_STRING_LENGTH) :: filename

  integer :: num_stations, num_channels_per_station
  integer :: num_waveforms  ! == num_stations * num_channels_per_station
  ! The number of channels per station is constant, as in SPECFEM
  integer :: start_time, nsamples  ! constant, as in SPECFEM
  ! Network names and station names are two different beast, as in SPECFEM
  ! network_names(i) is related to station_names(i)
  character(len=MAX_STRING_LENGTH), dimension(:), allocatable :: network_names
  character(len=MAX_STRING_LENGTH), dimension(:), allocatable :: station_names
  ! data. dimension = nsamples * num_channels_per_station * num_stations
  real, dimension(:, :, :), allocatable :: waveforms

  character(len=MAX_STRING_LENGTH) :: station_name, waveform_name, path

  !-- ASDF variables 
  integer :: file_id   ! HDF5 file id, also root group "/"
  integer :: station_exists, waveform_exists
  integer :: nsamples_infered

  !--- MPI variables
  integer :: myrank, mysize, comm
  !--- Loop variables
  integer :: i, j, k
  !--- Error variable
  integer :: ier

  ! alias mpi communicator
  call world_duplicate(comm)
  call world_size(mysize)

  !--------------------------------------------------------
  ! Setup data on each process.
  !--------------------------------------------------------
  filename = "synthetic.h5"

  num_stations = nrec_local
  num_channels_per_station = 3
  nsamples = seismo_current
  num_waveforms = num_stations * num_channels_per_station

  allocate(network_names(num_stations), stat=ier)
  allocate(station_names(num_stations), stat=ier)
  allocate(waveforms(nsamples, num_channels_per_station, num_stations), &
           stat=ier)

  do i = 1, num_stations
    write(network_names(i), '(a,i0.2)') "NTWK_", myrank
    write(station_names(i), '(a,i0.2,a,i0.2)') "STAT_", myrank, "_", i
  enddo

  !--------------------------------------------------------
  ! Read the ASDF file.
  !--------------------------------------------------------
print *, "initalizing ASDF"
  call ASDF_initialize_hdf5_f(ier);
  call ASDF_open_read_only_f(trim(filename) // C_NULL_CHAR, comm, file_id)

  do j = 1, num_stations
    station_name = trim(network_names(j)) // '.' // &
                   trim(station_names(j))
    call ASDF_station_exists_f(file_id, &
                               trim(station_name) // C_NULL_CHAR, &
                               station_exists)
    if (station_exists > 0) then
      do  i = 1, num_channels_per_station
        ! Generate dummy waveform names to match write_example.f90 output
        ! Note to fortran users:
        !   be sure to trim your strings and to append '\0' (i.e
        !   C_NULL_CHAR) to them.
        write(waveform_name, '(a,i0.2)') &
            trim(network_names(j)) // "." //      &
            trim(station_names(j)) // ".channel_", i
        call ASDF_waveform_exists_f(file_id, &
                                    trim(station_name) // C_NULL_CHAR, &
                                    trim(waveform_name) // C_NULL_CHAR, &
                                    waveform_exists)
        path = "/Waveforms/" // trim(station_name) // "/" &
            // trim(waveform_name) //  C_NULL_CHAR
        call ASDF_get_num_elements_from_path_f(file_id, path, nsamples_infered)
        if (nsamples_infered == nsamples) then
          call ASDF_read_full_waveform_f(file_id, &
                                         trim(path) // C_NULL_CHAR, &
                                         waveforms(:, i, j), ier)
          if (myrank == mysize-1) then
            print *, "-------------------------------------------------"
            print *, trim(waveform_name)
            print *, "-------------------------------------------------"
            print *, waveforms(:, i, j)
            print *, ""
            !call flush()
          endif
        endif
      enddo
    endif
  enddo

  !--------------------------------------------------------
  ! Clean up
  !--------------------------------------------------------
  call ASDF_finalize_hdf5_f(ier);
  call MPI_Finalize(ier)

end subroutine read_asdf



































  ! Parameters
  type(asdf_event),intent(inout) :: asdf_container

  integer, parameter :: MAX_STRING_LENGTH = 256

  character(len=MAX_STRING_LENGTH) :: filename

  !--- Data to be written to the ASDF file
  character(len=MAX_STRING_LENGTH) :: event_name
  character(len=MAX_STRING_LENGTH) :: quakeml
  character(len=MAX_STRING_LENGTH) :: provenance
  character(len=MAX_STRING_LENGTH) :: station_xml

  integer :: num_stations, num_channels_per_station
  integer :: num_waveforms  ! == num_stations * num_channels_per_station
  ! The number of channels per station is constant, as in SPECFEM

  integer :: start_time, nsamples  ! constant, as in SPECFEM
  double precision :: sampling_rate  ! idem

  ! Network names and station names are two different beast, as in SPECFEM
  ! network_names(i) is related to station_names(i)
  character(len=MAX_STRING_LENGTH), dimension(:), allocatable :: networks_names
  character(len=MAX_STRING_LENGTH), dimension(:), allocatable :: stations_names
  character(len=MAX_STRING_LENGTH), dimension(:), allocatable :: component_names

  ! data. dimension = nsamples * num_channels_per_station * num_stations
  real, dimension(:, :, :), allocatable :: waveforms

  !-- ASDF variables 
  !   These variables are used to know where further writes should be done.
  !   They have to be cleaned as soon as they become useless
  integer :: file_id   ! HDF5 file id, also root group "/"
  integer :: waveforms_grp  ! Group "/Waveforms/" 
  integer, dimension(:, :, :), allocatable :: data_ids

  !--- MPI variables
  integer :: mysize, comm
  !--- Loop variables
  integer :: i, j, k
  !--- Error variable
  integer :: ier

  !--- 'allgather' arrays. Variables that needs to be known by everyone in
  !    order to define ASDF groups and datasets or write them as attributes.
  integer, dimension(:), allocatable :: num_stations_gather
  integer :: max_num_stations_gather
  character(len=MAX_STRING_LENGTH), dimension(:,:), allocatable :: &
      station_names_gather, network_names_gather, component_names_gather
  integer, dimension(:,:), allocatable :: station_grps_gather
  integer, dimension(:), allocatable :: displs, rcounts

  ! temporary name built from network, station and channel names.
  character(len=MAX_STRING_LENGTH) :: waveform_name

  ! alias mpi communicator for ADIOS
  call world_duplicate(comm)
  call world_size(mysize)

  !--------------------------------------------------------
  ! Setup data on each process.
  !--------------------------------------------------------

  filename = "synthetic.h5"
  event_name = trim(event_name_SAC)
  provenance = "<provenance>"
! CMT to quakeml converter
  station_xml = "<station_xml>"

  num_stations = nrec_local
  num_channels_per_station = 3
  sampling_rate = DT
  nsamples = seismo_current
  start_time = seismo_offset*DT-t0+t_cmt_SAC
  num_waveforms = num_stations * num_channels_per_station

  call cmt_to_quakeml(quakeml)
  call generate_provenance(provenance)

  allocate(networks_names(num_stations), stat=ier)
  allocate(stations_names(num_stations), stat=ier)
  allocate(component_names(num_stations*num_channels_per_station), stat=ier)
  allocate(waveforms(nsamples, num_channels_per_station, num_stations), &
           stat=ier)

  do i = 1, num_stations
   write(networks_names(i), '(a)') asdf_container%network_array(i)
   write(stations_names(i), '(a)') asdf_container%receiver_name_array(i)
  enddo

  do i = 1, num_stations*num_channels_per_station
   write(component_names(i), '(a)') asdf_container%component_array(i)
  enddo

  ! -- We do not care about seeding.
  call random_number(waveforms)

  !--------------------------------------------------------
  ! ASDF variables
  !--------------------------------------------------------
  ! Find how many stations are managed by each allgatheress
  allocate(num_stations_gather(mysize))
  call all_gather_all_i(num_stations, num_stations_gather, mysize)
  ! find the largest number of stations per allgatheress
  max_num_stations_gather = maxval(num_stations_gather)

  allocate(displs(mysize))
  allocate(rcounts(mysize))

  ! Everyone should know about each and every station name
  allocate(station_names_gather(max_num_stations_gather, mysize))
  allocate(network_names_gather(max_num_stations_gather, mysize))
  allocate(component_names_gather(max_num_stations_gather*3, mysize))

  ! The number of stations is not constant across processes
  do i = 1, mysize
    displs(i) = (i-1) * max_num_stations_gather * max_string_length
    rcounts(i) = num_stations_gather(i) * max_string_length
  enddo
  call all_gather_all_ch(stations_names, &
                         num_stations * MAX_STRING_LENGTH, &
                         station_names_gather, &
                         rcounts, &
                         displs, &
                         max_num_stations_gather, &
                         MAX_STRING_LENGTH, &
                         mysize)
  call all_gather_all_ch(networks_names, &
                         num_stations * MAX_STRING_LENGTH, &
                         network_names_gather, &
                         rcounts, &
                         displs, &
                         max_num_stations_gather, &
                         MAX_STRING_LENGTH, &
                         mysize)
  call all_gather_all_ch(component_names, &
                         num_stations*3*MAX_STRING_LENGTH, &
                         component_names_gather, &
                         rcounts*3, &
                         displs*3, &
                         max_num_stations_gather*3, &
                         MAX_STRING_LENGTH, &
                         mysize)
  deallocate(displs)
  deallocate(rcounts)

  allocate(station_grps_gather(max_num_stations_gather, mysize))

  allocate(data_ids(num_channels_per_station, &
                    max_num_stations_gather, &
                    mysize))

  !--------------------------------------------------------
  ! write ASDF 
  !--------------------------------------------------------
print *, "initializing ASDF"
  call ASDF_initialize_hdf5_f(ier);
  call ASDF_create_new_file_f(trim(filename), comm, file_id)

  call ASDF_write_string_attribute_f(file_id, "file_format" // C_NULL_CHAR, &
                                     "ASDF" // C_NULL_CHAR, ier)
  call ASDF_write_string_attribute_f(file_id, "file_format_version" // C_NULL_CHAR, &
                                     "0.0.2" // C_NULL_CHAR, ier)

  call ASDF_write_quakeml_f(file_id, trim(quakeml), ier)
  call ASDF_write_provenance_data_f(file_id, trim(provenance), ier)
  call ASDF_write_auxiliary_data_f(file_id, ier)

  call ASDF_create_waveforms_group_f(file_id, waveforms_grp)

  print *, "defining waveforms"
  do k = 1, mysize
    do j = 1, num_stations_gather(k)
      call station_to_stationxml(station_names_gather(j,k), station_xml)
      call ASDF_create_stations_group_f(waveforms_grp,   &
           trim(network_names_gather(j, k)) // "." //      &
           trim(station_names_gather(j,k)) // C_NULL_CHAR, &
           trim(station_xml) // C_NULL_CHAR,             &
           station_grps_gather(j, k))
       do  i = 1, num_channels_per_station
       ! Generate unique dummy waveform names
        write(waveform_name, '(a)') &
           trim(network_names_gather(j,k)) // "." //      &
           trim(station_names_gather(j,k)) // ".." // trim(component_names_gather(i+(3*(j-1)),k)) &
           // "__2014-04-04T01:33:37__2014-04-04T02:15:10__synthetic"
        ! Create the dataset where waveform will be written later on.
        call ASDF_define_waveform_f(station_grps_gather(j,k), &
             nsamples, start_time, sampling_rate, &
             trim(event_name) // C_NULL_CHAR, &
             trim(waveform_name) // C_NULL_CHAR, &
             data_ids(i, j, k))
        if (nrec_local > 0) then
          waveforms(:,i,j) = asdf_container%records(i+(3*(j-1)))%record
        endif
      enddo
    enddo
  enddo

  print *, "writing waveforms"
  do j = 1, num_stations
   do i = 1, num_channels_per_station
      call ASDF_write_full_waveform_f(data_ids(i, j, myrank+1), &
                                      waveforms(:, i, j), ier)
    enddo
  enddo

  !--------------------------------------------------------
  ! Clean up
  !--------------------------------------------------------

  do k = 1, mysize
    do j = 1, num_stations_gather(k)
      call ASDF_close_group_f(station_grps_gather(j, k), ier)
      do i = 1, num_channels_per_station
        call ASDF_close_dataset_f(data_ids(i, j, k), ier)
      enddo
    enddo
  enddo

  call ASDF_close_group_f(waveforms_grp, ier)
  call ASDF_finalize_hdf5_f(ier);

  deallocate(data_ids)
  deallocate(station_grps_gather)
  deallocate(station_names_gather)
  deallocate(network_names_gather)
  deallocate(num_stations_gather)

  deallocate(waveforms)
  deallocate(stations_names)
  deallocate(networks_names)

end subroutine read_asdf
