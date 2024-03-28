!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
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

  ! only main process writes out to main output file
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
  !
  !call system(trim(command))
  !
  ! standard f2008
  call execute_command_line(trim(command))

  end subroutine system_command

!
!-------------------------------------------------------------------------------------------------
!

  subroutine flush_stdout()

! flushes possible left-overs from print-statements

  implicit none

  logical :: is_connected

  ! note: Cray systems don't flush print statements before ending with an MPI abort,
  !       which often omits debugging statements with print before it.
  !
  !       to check which unit is used for standard output, one might also use a Fortran2003 module iso_Fortran_env:
  !         use, intrinsic :: iso_Fortran_env, only: output_unit

  ! checks default stdout unit 6
  inquire(unit=6,opened=is_connected)
  if (is_connected) &
    flush(6)

  ! checks Cray stdout unit 101
  inquire(unit=101,opened=is_connected)
  if (is_connected) &
    flush(101)

  end subroutine flush_stdout
