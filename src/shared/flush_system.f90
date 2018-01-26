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
! the Free Software Foundation; either version 3 of the License, or
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

!-------------------------------------------------------------------------------------------------
!
! I/O wrapper function
!
!-------------------------------------------------------------------------------------------------

  subroutine flush_IMAIN()

  use constants

  implicit none

  ! only master process writes out to main output file
  ! file I/O in Fortran is buffered by default
  !
  ! note: Fortran2003 includes a FLUSH statement
  !          which is implemented by most compilers by now
  !
  ! otherwise:
  !   a) comment out the line below
  !   b) try to use instead: call flush(IMAIN)

  flush(IMAIN)

  end subroutine flush_IMAIN


!-------------------------------------------------------------------------------------------------

  subroutine system_command(command)

  use constants

  implicit none

  character(len=MAX_STRING_LEN), intent(in):: command

  ! note: system() is a GNU Fortran extension, compilers might complain about this command
  !       (non-standard Fortran2003)
  !
  ! otherwise:
  !   a) comment out the line below
  !   b) try to use other function implementation: e.g. Fortran 2008 supports
  !      call execute_command_line(trim(command))

  call system(trim(command))

  end subroutine system_command

