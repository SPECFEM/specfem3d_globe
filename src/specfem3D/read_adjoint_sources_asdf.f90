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
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
! implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!==============================================================================

subroutine read_adjoint_sources_asdf(file_id)

  use iso_c_binding
  use specfem_par, only : asdf_file_handle

  implicit none

  integer, parameter :: MAX_STRING_LENGTH = 256

  character(len=MAX_STRING_LENGTH) :: filename

  integer :: num_stations, num_channels_per_station
  integer :: num_waveforms  ! == num_stations * num_channels_per_station
  ! The number of channels per station is constant, as in SPECFEM
  integer :: nsamples  ! constant, as in SPECFEM
  ! Network names and station names are two different beast, as in
  ! SPECFEM
  ! network_names(i) is related to station_names(i)
  character(len=MAX_STRING_LENGTH), dimension(:), allocatable :: network_names
  character(len=MAX_STRING_LENGTH), dimension(:), allocatable :: station_names
  ! data. dimension = nsamples * num_channels_per_station * num_stations
  real, dimension(:, :, :), allocatable :: waveforms
  real, dimension(:,:,:), allocatable :: adjoint

  character(len=MAX_STRING_LENGTH) :: station_name, waveform_name, path
  character(len=MAX_STRING_LENGTH) :: file_format

  !-- ASDF variables
  integer,intent(in) :: file_id   ! HDF5 file id, also root group "/"
  integer :: station_exists, waveform_exists
  integer :: nsamples_infered

  ! C/Fortran interop for C-allocated strings

  type(c_ptr) :: cptr
  character, pointer :: fptr(:)

  !--- MPI variables
  integer :: myrank, mysize, comm
  !--- Loop variables
  integer :: i, j
  !--- Error variable
  integer :: ier

  ! alias mpi communicator
  call world_get_comm(comm)
  call world_size(mysize)

  !--------------------------------------------------------
  ! Read the ASDF file.
  !--------------------------------------------------------

print *, "b"

  call ASDF_read_str_attr_f(asdf_file_handle, "/", "file_format", cptr, ier)
print *, "c"
  call c_f_pointer(cptr, fptr, [4])
print *, "d"
  ! print "ASDF" on each processer
  print *, fptr

  !allocate(adjoint(20000, 3, 1), &
   !        stat=ier)

    ! get station name
    !station_name = "BW_ALFO" ! read from adjoint station list in files

    !call ASDF_adjoint_source_exists_f(file_id, &
    !                            trim(station_name)//C_NULL_CHAR, &
  !                              station_exists)
   ! if (station_exists) > 0) then
  !do j = 1, num_channels_per_station

   ! write(adjoint_source_name, '(a)') &
  !      trim(station_name)//"_"//trim(component_name(j))
 !     "BW_ALFO
   ! call ASDF_get_num_elements_from_path_f(file_id, "/AuxiliaryData/&
   !       AdjointSource/"//trim(station_name)//"_EHE"//C_NULL_CHAR,nsamples_infered)
   ! print *, nsamples_infered
   ! allocate(adjoint(nsamples_infered, 1, 1), stat=ier)
   ! call ASDF_read_full_waveform_f(file_id, "/AuxiliaryData/&
   !       AdjointSource/"//trim(station_name)//"_EHE"//C_NULL_CHAR,adjoint(:,j,1), ier)
   ! endif

 ! enddo

end subroutine read_adjoint_sources_asdf
