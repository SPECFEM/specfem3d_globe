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
!> Open the ASDF file for reading adjoint sources
subroutine asdf_setup(asdf_file_handle)

  use iso_c_binding

  implicit none

  integer, parameter :: MAX_STRING_LENGTH = 256
  integer, intent(inout) :: asdf_file_handle

  character(len=MAX_STRING_LENGTH) :: filename

  ! local parameters
  integer :: comm

  call world_get_comm(comm)

  filename = "synthetic.h5"

  call ASDF_open_read_only_f(trim(filename) // C_NULL_CHAR, comm, asdf_file_handle)

end subroutine asdf_setup

!==============================================================================
!> Close the ASDF file for reading adjoint sources
subroutine asdf_cleanup()

  implicit none

  ! local parameters
  integer :: myrank
  integer :: ier

  call world_rank(myrank)
  call synchronize_all()

  call ASDF_finalize_hdf5_f(ier)
  if (ier /= 0 ) stop 'Error cleaning up ASDF: calling asdf_finalize() routine failed'

end subroutine asdf_cleanup
