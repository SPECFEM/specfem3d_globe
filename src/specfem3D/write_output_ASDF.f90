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
  use specfem_par, only : event_name_SAC, myrank

  ! Parameters
  type(asdf_event),intent(inout) :: asdf_container
  integer,intent(in) :: total_seismos_local

  ! Variables
  integer :: ier

  asdf_container%event = trim(event_name_SAC)

  allocate (asdf_container%receiver_name_array(total_seismos_local), STAT=ier)
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
subroutine store_asdf_data(asdf_container, seismogram_tmp, irec_local, &
                           irec, chn, iorientation)

  use asdf_data
  use specfem_par,only: &
    station_name,network_name, &
    seismo_current, NTSTEP_BETWEEN_OUTPUT_SEISMOS

  use specfem_par, only: myrank
  use constants

  implicit none

  ! Parameters
  type(asdf_event),intent(inout) :: asdf_container
  character(len=4),intent(in) :: chn
  integer,intent(in) :: irec_local, irec
  real(kind=CUSTOM_REAL),dimension(5,NTSTEP_BETWEEN_OUTPUT_SEISMOS),intent(in) :: seismogram_tmp
  integer,intent(in) :: iorientation
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
!! \param asdf_container The ASDF data container for the waveforms
subroutine write_asdf(asdf_container)

  use asdf_data
  use iso_c_binding
  use iso_fortran_env
  use specfem_par

  implicit none
  type(asdf_event),intent(inout) :: asdf_container

  ! Parameters
  integer, parameter :: MAX_STRING_LENGTH = 7555

  !--- Character strings to be written to the ASDF file
  character(len=MAX_STRING_LENGTH) :: quakeml
  character(len=MAX_STRING_LENGTH) :: stationxml
  character(len=MAX_STRING_LENGTH) :: sf_constants
  character(len=MAX_StRING_LENGTH) :: sf_parfile

  integer :: num_stations !
  integer :: nsamples  ! constant, as in SPECFEM
  integer(int64) :: start_time
  double precision :: sampling_rate
  double precision :: startTime

  ! Network names and station names are two different beasts, as in SPECFEM
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

  ! C/Fortran interop for C-allocated strings
  integer :: len
  type(c_ptr) :: cptr
  character, pointer :: fptr(:)
  character, dimension(:), allocatable, TARGET :: provenance

  ! Time variables
  character(len=25) :: start_time_string, end_time_string

  ! alias mpi communicator
  call world_duplicate(comm)
  call world_size(mysize)

  !--------------------------------------------------------
  ! Setup data on each process.
  !--------------------------------------------------------

  num_stations = nrec_local
  sampling_rate = DT
  nsamples = seismo_current

  ! Generate minimal QuakeML for SPECFEM3D_GLOBE
  call cmt_to_quakeml(quakeml)

  ! Calculate start_time based on SAC header values
  call get_time(startTime, start_time_string, end_time_string)
  start_time = startTime*1000000000_int64 ! convert to nanoseconds

  call ASDF_clean_provenance_f(cptr)
print *, start_time_string(1:19)
print *, end_time_string(1:19)
  call ASDF_generate_sf_provenance_f(start_time_string(1:19)//C_NULL_CHAR, end_time_string(1:19)//C_NULL_CHAR, cptr, len)
  call c_f_pointer(cptr, fptr, [len])
  allocate(provenance(len+1))
  provenance(1:len) = fptr(1:len)
  provenance(len+1) = C_NULL_CHAR

  allocate(networks_names(num_stations), stat=ier)
  allocate(stations_names(num_stations), stat=ier)
  allocate(component_names(num_stations*3), stat=ier)
  allocate(waveforms(nsamples, 3, num_stations), &
           stat=ier)

  do i = 1, num_stations
   write(networks_names(i), '(a)') asdf_container%network_array(i)
   write(stations_names(i), '(a)') asdf_container%receiver_name_array(i)
  enddo

  do i = 1, num_stations*3
   write(component_names(i), '(a)') asdf_container%component_array(i)
  enddo

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

  allocate(data_ids(3, &
                    max_num_stations_gather, &
                    mysize))

  !--------------------------------------------------------
  ! write ASDF 
  !--------------------------------------------------------

  call ASDF_initialize_hdf5_f(ier);
  call ASDF_create_new_file_f("synthetic.h5", comm, file_id)

  call ASDF_write_string_attribute_f(file_id, "file_format" // C_NULL_CHAR, &
                                     "ASDF" // C_NULL_CHAR, ier)
  call ASDF_write_string_attribute_f(file_id, "file_format_version" // C_NULL_CHAR, &
                                     "0.0.2" // C_NULL_CHAR, ier)

  call ASDF_write_quakeml_f(file_id, trim(quakeml), ier)
  call ASDF_write_provenance_data_f(file_id, provenance(1:len), ier)
  call read_file("setup/constants.h", sf_constants)
  call read_file("DATA/Par_file", sf_parfile)
  call ASDF_write_auxiliary_data_f(file_id, trim(sf_constants), trim(sf_parfile), ier)

  call ASDF_create_waveforms_group_f(file_id, waveforms_grp)


  do k = 1, mysize
    do j = 1, num_stations_gather(k)
      call station_to_stationxml(station_names_gather(j,k), k, stationxml)
      call ASDF_create_stations_group_f(waveforms_grp,   &
           trim(network_names_gather(j, k)) // "." //      &
           trim(station_names_gather(j,k)) // C_NULL_CHAR, &
           trim(stationxml) // C_NULL_CHAR,             &
           station_grps_gather(j, k))
       do  i = 1, 3
       ! Generate unique dummy waveform names
        write(waveform_name, '(a)') &
           trim(network_names_gather(j,k)) // "." //      &
           trim(station_names_gather(j,k)) // ".." // trim(component_names_gather(i+(3*(j-1)),k)) &
           // "__"//trim(start_time_string(1:19))//"__"//trim(end_time_string(1:19))//"__synthetic"
        ! Create the dataset where waveform will be written later on.
        call ASDF_define_waveform_f(station_grps_gather(j,k), &
             nsamples, start_time, sampling_rate, &
             trim(event_name_SAC) // C_NULL_CHAR, &
             trim(waveform_name) // C_NULL_CHAR, &
             data_ids(i, j, k))
        if (nrec_local > 0) then
          waveforms(:,i,j) = asdf_container%records(i+(3*(j-1)))%record
        endif
      enddo
    enddo
  enddo

  do j = 1, num_stations
   do i = 1, 3
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
      do i = 1, 3
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
  character(len=25) :: lon_str, lat_str, dep_str

  write(lon_str, "(g12.5)") cmt_lat
  write(lat_str, "(g12.5)") cmt_lon
  write(dep_str, "(g12.5)") cmt_depth

  quakemlstring = '<q:quakeml xmlns="http://quakeml.org/xmlns/bed/1.2"'//&
                  ' xmlns:q="http://quakeml.org/xmlns/quakeml/1.2">'//&
                  '<eventParameters publicID="smi:local/592edbfc-2482-481e-80fd-03810308104b">'//&
                  '<event publicID="smi:service.iris.edu/fdsnws/event/1/query?eventid=656970">'//&
                  '<type>earthquake</type>'//&
                  '<origin publicID="smi:www.iris.edu/spudservice/momenttensor/gcmtid/B090198B#cmtorigin">'//&
                  '<time><value>2014-09-01T10:29:54.500000Z</value></time>'//&
                  '<latitude>'//&
                  '<value>'//trim(lat_str)//'</value>'//&
                  '</latitude>'//&
                  '<longitude>'//&
                  '<value>'//trim(lon_str)//'</value>'//&
                  '</longitude>'//&
                  '<depth>'//&
                  '<value>'//trim(dep_str)//'</value>'//&
                  '</depth>'//&
                  '</origin>'//&
                  '</event>'//&
                  '</eventParameters>'//&
                  '</q:quakeml>'

end subroutine cmt_to_quakeml

subroutine get_time(startTime, start_time_string, end_time_string)

  use specfem_par,only:&
    yr_SAC, jda_SAC, ho_SAC, mi_SAC, sec_SAC, DT, seismo_current

  implicit none
  character(len=*) :: start_time_string, end_time_string
  double precision, intent(inout) :: startTime
  character(len=4) :: yr
  character(len=2) :: hr, minute
  character(len=15) :: second
  double precision :: end_time_sec, end_sec
  integer :: end_min, end_ho, temp

  write(yr, "(I4.4)") yr_SAC
  write(hr, "(I2.2)") ho_SAC
  write(minute, "(I2.2)") mi_SAC
  write(second, "(F5.2)") sec_SAC

  start_time_string = trim(yr)//"-06-09T"//trim(hr)//':'//trim(minute)//':'//trim(second)
  startTime =(yr_SAC-1970)*365.25*24*60*60.0d0+(jda_SAC-1.0)*24*60*60.0d0+ho_SAC*(60.0d0*60.0d0)+mi_SAC*60.0d0+sec_SAC
    
  end_time_sec = (seismo_current-1)/DT
  end_ho = ho_SAC
  end_min = mi_SAC
  end_sec = sec_SAC

  ! add end_time_sec to the start time
  do
    if (end_time_sec/60.0 > 1.0) then
      temp = end_min + int(end_time_sec/60.0)
      end_min = (temp-60.0)
      end_ho = end_ho + 1
      end_sec = end_sec + temp
      end_sec = 60.0*(end_time_sec/60.0-int(end_time_sec/60.0))+sec_SAC
      exit
    else
      end_sec = end_sec + end_time_sec
      exit
    endif
  enddo 

  !todo: fix end_time so if it rolls over to a new year the date is correct 
  write(yr, "(I4.4)") yr_SAC
  write(hr, "(I2.2)") end_ho
  write(minute, "(I2.2)") end_min
  write(second, "(F5.2)") end_sec

  end_time_string = trim(yr)//"-06-09T"//trim(hr)//':'//trim(minute)//':'//trim(second(1:2))

end subroutine get_time

subroutine station_to_stationxml(station_name, irec, stationxmlstring)

  use specfem_par,only:&
    stlat, stlon, network_name

  implicit none
  character(len=*) :: station_name, stationxmlstring
  character(len=25) :: station_lat, station_lon
  integer, intent(in) :: irec

  write(station_lat, "(g12.5)") stlat(irec)
  write(station_lon, "(g12.5)") stlon(irec) 

  stationxmlstring = '<FDSNStationXML schemaVersion="1.0" xmlns="http://www.fdsn.org/xml/station/1">'//&
                     '<Source>Erdbebendienst Bayern</Source>'//&
                      '<Module>fdsn-stationxml-converter/1.0.0</Module>'//&
                      '<ModuleURI>http://www.iris.edu/fdsnstationconverter</ModuleURI>'//&
                      '<Created>2014-03-03T11:07:06+00:00</Created>'//&
                      '<Network code="IU"><Station code="'//trim(station_name)//'"'//&
                      ' startDate="2006-12-16T00:00:00+00:00">'//&
                      '<Latitude unit="DEGREES">'//trim(station_lat)//'</Latitude>'//&
                      '<Longitude unit="DEGREES">'//trim(station_lon)//'</Longitude>'//&
                      '<Elevation>565.0</Elevation>'//&
                      '<Site><Name>SPECFEM3D_GLOBE</Name></Site>'//&
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

end subroutine read_file
