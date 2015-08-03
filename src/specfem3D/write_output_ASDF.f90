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

!> Initializes the data structure for ASDF
!! \param asdf_container The ASDF data structure
!! \param total_seismos_local The number of records on the local processor
subroutine init_asdf_data(asdf_container,total_seismos_local)

  use asdf_data
  use specfem_par, only : event_name_SAC,myrank

  ! Parameters
  type(asdf_event),intent(inout) :: asdf_container
  integer,intent(in) :: total_seismos_local

  ! Variables
  integer :: ier

  asdf_container%event = trim(event_name_SAC)

  allocate (asdf_container%receiver_name_array(total_seismos_local), &
            STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%network_array(total_seismos_local), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%component_array(total_seismos_local), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%records(total_seismos_local), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')

end subroutine init_asdf_data

!
!-------------------------------------------------------------------------------------------------
!

!> Stores the records into the ASDF structure
!! \param asdf_container The ASDF data structure
!! \param seismogram_tmp The current seismogram to store
!! \param irec_local The local index of the receivers on the local processor
!! \param irec The global index of the receiver
!! \param chn The broadband channel simulated
!! \param iorientation The recorded seismogram's orientation direction
!! \param phi The angle used for calculating azimuth and incident angle
subroutine store_asdf_data(asdf_container, seismogram_tmp, irec_local, &
                           irec, chn, iorientation, phi)

  use asdf_data
  use specfem_par,only: &
    station_name,network_name,stlat,stlon,stele,stbur, &
    DT,t0, seismo_offset,seismo_current, NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
    yr=>yr_SAC,jda=>jda_SAC,ho=>ho_SAC,mi=>mi_SAC,sec=>sec_SAC, &
    tshift_cmt=>t_cmt_SAC, &
    cmt_lat=>cmt_lat_SAC,cmt_lon=>cmt_lon_SAC, &
    cmt_depth=>cmt_depth_SAC

  use specfem_par, only: myrank
  use constants

  implicit none

  ! Parameters
  type(asdf_event),intent(inout) :: asdf_container
  character(len=4),intent(in) :: chn
  integer,intent(in) :: irec_local, irec
  real(kind=CUSTOM_REAL),dimension(5,NTSTEP_BETWEEN_OUTPUT_SEISMOS),intent(in) :: seismogram_tmp
  integer,intent(in) :: iorientation
  double precision,intent(in) :: phi
  ! local Variables
  integer :: length_station_name, length_network_name
  integer :: ier, i

  ! trace index
  i = (irec_local-1)*(3) + (iorientation)

  length_station_name = len_trim(station_name(irec))
  length_network_name = len_trim(network_name(irec))
  asdf_container%receiver_name_array(i) = station_name(irec)(1:length_station_name)
  asdf_container%network_array(i) = network_name(irec)(1:length_network_name)
  asdf_container%component_array(i) = chn

  print *, phi
  allocate (asdf_container%records(i)%record(seismo_current), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocating ASDF container failed.')
  asdf_container%records(i)%record(1:seismo_current) = seismogram_tmp(iorientation, 1:seismo_current)

end subroutine store_asdf_data

!
!-------------------------------------------------------------------------------------------------
!

!> Closes the ASDF data structure by deallocating all arrays
!! \param asdf_container The ASDF data structure
!! \param total_seismos_local The number of seismograms on the local processor
subroutine close_asdf_data(asdf_container, total_seismos_local)

  use asdf_data
  ! Parameters
  type(asdf_event),intent(inout) :: asdf_container
  integer,intent(in) :: total_seismos_local
  !Variables
  integer :: i

  do i = 1, total_seismos_local
    deallocate(asdf_container%records(i)%record)
  enddo
  deallocate (asdf_container%receiver_name_array)
  deallocate (asdf_container%network_array)
  deallocate (asdf_container%component_array)

end subroutine close_asdf_data

!
!-------------------------------------------------------------------------------------------------
!

!> Writes the ASDF data structure to the file
!! \param asdf_container The ASDF data structure
subroutine write_asdf(asdf_container)

  use asdf_data
  use iso_c_binding
  use specfem_par

  implicit none
  ! Parameters
  type(asdf_event),intent(inout) :: asdf_container

  integer, parameter :: MAX_STRING_LENGTH = 7555

  character(len=MAX_STRING_LENGTH) :: filename

  !--- Data to be written to the ASDF file
  character(len=MAX_STRING_LENGTH) :: event_name
  character(len=MAX_STRING_LENGTH) :: quakeml
  character(len=MAX_STRING_LENGTH) :: provenance
  character(len=MAX_STRING_LENGTH) :: station_xml
  character(len=MAX_STRING_LENGTH) :: sf_constants, sf_parfile

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
  character(len=5) :: startTime, endTime

  ! C/Fortran interop for C-allocated strings
  integer :: len

  type(c_ptr) :: cptr
  character, pointer :: fptr(:)

  ! alias mpi communicator for ADIOS
  call world_duplicate(comm)
  call world_size(mysize)

  !--------------------------------------------------------
  ! Setup data on each process.
  !--------------------------------------------------------

  filename = "synthetic.h5"
  event_name = trim(event_name_SAC)
  station_xml = "<station_xml>"
  provenance = '<prov:document xmlns:prov="http://www.w3.org/ns/prov#" xmlns:seis_prov="http://seisprov.org/seis_prov/0.1/#"'//&
               ' xmlns:xsd="http://www.w3.org/2001/XMLSchema" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'//&
               '<prov:activity prov:id="seis_prov:sp001_ws_910a6ce"><prov:label>Waveform Simulation</prov:label>' //&
               '<prov:type xsi:type="xsd:string">seis_prov:waveform_simulation</prov:type></prov:activity></prov:document>'

  num_stations = nrec_local
  num_channels_per_station = 3
  sampling_rate = DT
  nsamples = seismo_current
  start_time = seismo_offset*DT-t0+t_cmt_SAC
  num_waveforms = num_stations * num_channels_per_station

  call cmt_to_quakeml(quakeml)

  !write(startTime, "(F5.2)") starddt_time
  !write(endTime, "(F5.2)") start_time
  call ASDF_clean_provenance_f(cptr)
  call ASDF_generate_sf_provenance_f("2014-01-01T12:15:03", "2014-01-01T02:15:03", cptr, len)
  call c_f_pointer(cptr, fptr, [len])
  print *, fptr


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

  call ASDF_initialize_hdf5_f(ier);
  call ASDF_create_new_file_f(trim(filename), comm, file_id)

  call ASDF_write_string_attribute_f(file_id, "file_format" // C_NULL_CHAR, &
                                     "ASDF" // C_NULL_CHAR, ier)
  call ASDF_write_string_attribute_f(file_id, "file_format_version" // C_NULL_CHAR, &
                                     "0.0.2" // C_NULL_CHAR, ier)


  call ASDF_write_quakeml_f(file_id, quakeml, ier)
  call ASDF_write_provenance_data_f(file_id, trim(provenance), ier)
  call read_file("setup/constants.h", sf_constants)
  call read_file("DATA/Par_file", sf_parfile)
  call ASDF_write_auxiliary_data_f(file_id, trim(sf_constants), trim(sf_parfile), ier)

  call ASDF_create_waveforms_group_f(file_id, waveforms_grp)

  print *, "defining waveforms"
  do k = 1, mysize
    do j = 1, num_stations_gather(k)
      call station_to_stationxml(station_names_gather(j,k), k, station_xml)
      print *, trim(station_xml)
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

end subroutine write_asdf

subroutine cmt_to_quakeml(quakemlstring)

  use specfem_par,only:&
    cmt_lat=>cmt_lat_SAC,cmt_lon=>cmt_lon_SAC,cmt_depth=>cmt_depth_SAC

  implicit none
  character(len=*) :: quakemlstring
  character(len=5) :: lon_str, lat_str, dep_str

  write(lon_str, "(F5.2)") cmt_lat
  write(lat_str, "(F5.2)") cmt_lon
  write(dep_str, "(F5.2)") cmt_depth

  quakemlstring = '<q:quakeml xsi:schemaLocation="http://quakeml.org/schema/xsd/QuakeML-1.2.xsd" '//&
                  'xmlns="http://quakeml.org/xmlns/bed/1.2" xmlns:q="http://quakeml.org/xmlns/quakeml/1.2" xmlns:xsi="http://'//&
                  'www.w3.org/2001/XMLSchema-instance">'//&
                  !'<eventParameters publicID="smi:service.iris.edu/fdsnws/event/1/query">'//&
                  !'<longitude><value>'//trim(lon_str)//'</value></longitude><latitude><value>'//trim(lat_str)//&
                  !'</value></latitude><depth><value>'//trim(dep_str)//'</value></depth>'//&
                  !'</eventParameters>'//&
                  '</q:quakeml>'

end subroutine cmt_to_quakeml

subroutine station_to_stationxml(station_name, irec, stationxmlstring)

  use specfem_par,only:&
    stlat, stlon

  implicit none
  character(len=*) :: station_name, stationxmlstring
  character(len=5) :: station_lat, station_lon
  integer, intent(in) :: irec

  write(station_lat, "(F5.2)") stlat(irec)
  write(station_lon, "(F5.2)") stlon(irec) 

  stationxmlstring = '<FDSNStationXML schemaVersion="1.0" xmlns="http://www.fdsn.org/xml/station/1">'//&
                     '<Source>Erdbebendienst Bayern</Source>'//&
                      '<Module>fdsn-stationxml-converter/1.0.0</Module>'//&
                      '<ModuleURI>http://www.iris.edu/fdsnstationconverter</ModuleURI>'//&
                      '<Created>2014-03-03T11:07:06+00:00</Created>'//&
                      '<Network code="IU"><Station code="ANTO" startDate="2006-12-16T00:00:00+00:00">'//&
                      '<Latitude unit="DEGREES">'//trim(station_lat)//'</Latitude>'//&
                      '<Longitude unit="DEGREES">'//trim(station_lon)//'</Longitude>'//&
                      '<Elevation>565.0</Elevation>'//&
                      '<Site><Name>Fuerstenfeldbruck, Bavaria, GR-Net</Name></Site>'//&
                      '<CreationDate>2006-12-16T00:00:00+00:00</CreationDate>'//&
                      '</Station></Network>'//&
                      '</FDSNStationXML>'

end subroutine station_to_stationxml

subroutine read_file(filename, filestring)

  implicit none
  character(len=*) :: filestring
  character(len=*) :: filename
  integer :: file_size

  open(10, file=filename, status='old')
  inquire(unit=10, size=file_size)
  close(10)

  open(10, file=filename, status='old', &
         recl=file_size, form='unformatted', access='direct')
  read (10, rec=1) filestring
  close(10)

  print *, trim(filestring)

end subroutine read_file
