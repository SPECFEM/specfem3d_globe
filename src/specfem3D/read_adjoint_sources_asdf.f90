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

subroutine read_adjoint_sources_asdf(filename, index_start, index_end, icomp, index_i, adj_src)

  use constants,only: CUSTOM_REAL, MAX_STRING_LEN
  use iso_c_binding
  use specfem_par, only : asdf_file_handle, NDIM

  implicit none

  integer, intent(inout) :: icomp, index_start, index_end,index_i
  real(kind=CUSTOM_REAL), dimension(NDIM,*),intent(inout) :: adj_src

  character(len=MAX_STRING_LEN) :: filename

  !-- ASDF variables
  integer :: station_exists
  integer :: nsamples_infered

  ! C/Fortran interop for C-allocated strings

  type(c_ptr) :: cptr
  character, pointer :: fptr(:)

  integer :: i, j
  !--- Error variable
  integer :: ier

  !--------------------------------------------------------
  ! Read the ASDF file.
  !--------------------------------------------------------

print *, asdf_file_handle
  call ASDF_read_str_attr_f(asdf_file_handle, "/", "file_format", cptr, ier)
  call c_f_pointer(cptr, fptr, [4])
print *, fptr
stop
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
