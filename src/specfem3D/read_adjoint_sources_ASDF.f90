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

subroutine read_adjoint_sources_ASDF(adj_source_name, adj_source, index_start, index_end)

  use iso_c_binding, only: C_NULL_CHAR
  use specfem_par, only: NDIM,CUSTOM_REAL, myrank, current_asdf_handle

  implicit none

  integer :: itime
  integer :: index_start, index_end
  real,dimension(*),intent(out) :: adj_source
  character(len=*) :: adj_source_name
  !--- Error variable
  integer ier

  call ASDF_read_full_waveform_f(current_asdf_handle, "AuxiliaryData/AdjointSource/"//trim(adj_source_name) // C_NULL_CHAR, &
        adj_source, ier)

  if (ier /= 0) then
    print *,'Error reading adjoint source: ',trim(adj_source_name)
    print *,'rank ',myrank,' - time step: ',itime,' index_start:',index_start,' index_end: ',index_end
    print *,'  ',trim(adj_source_name)//'has wrong length, please check with your simulation duration'
    call exit_MPI(myrank,'Adjoint source '//trim(adj_source_name)//' has wrong length, please check with your simulation duration')
  endif

end subroutine read_adjoint_sources_ASDF

!------------------------------------------------------------------------------------------

subroutine check_adjoint_sources_ASDF(irec, nadj_sources_found)

  use iso_c_binding, only: C_NULL_CHAR
  use specfem_par
  use write_seismograms_mod, only: band_instrument_code

  implicit none

  integer,intent(in) :: irec
  integer :: nadj_sources_found
  integer :: nsamples_infered
  integer :: icomp
  integer :: ier

  ! local parameters
  character(len=MAX_STRING_LEN) :: adj_filename,adj_source_file
  character(len=3),dimension(NDIM) :: comp
  character(len=2) :: bic

  adj_source_file = trim(network_name(irec))//'_'//trim(station_name(irec))

  ! bandwidth code
  call band_instrument_code(DT,bic)
  comp(1) = bic(1:2)//'N'
  comp(2) = bic(1:2)//'E'
  comp(3) = bic(1:2)//'Z'

  ! loops over file components E/N/Z
  do icomp = 1,NDIM

    ! opens adjoint source file for this component
    adj_filename = trim(adj_source_file) // '_'// comp(icomp)

    ! checks if adjoint source exists in ASDF file
    call ASDF_adjoint_source_exists_f(current_asdf_handle, trim(adj_filename) // C_NULL_CHAR, ier)
   
    if (ier /= 0) then
      ! adjoint source not found
      ! stops simulation
      call exit_MPI(myrank,'adjoint source '//trim(adj_filename)//' not found, please check STATIONS_ADJOINT file')
    endif

    ! checks length of file
    call ASDF_get_num_elements_from_path_f(current_asdf_handle,&
       "AuxiliaryData/AdjointSource/" // trim(adj_filename) // C_NULL_CHAR, nsamples_infered, ier)

    ! checks length
    if (nsamples_infered /= NSTEP) then
      print *,'adjoint source error: ',trim(adj_filename),' has length',nsamples_infered,' but should be',NSTEP
      call exit_MPI(myrank,&
        'file '//trim(adj_filename)//' length is wrong, please check your adjoint sources and your simulation duration')
    endif

    ! updates counter for found files
    nadj_sources_found = nadj_sources_found + 1

  enddo

end subroutine check_adjoint_sources_ASDF
