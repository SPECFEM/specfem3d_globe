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

!-------------------------------------------------------------------------------------------------
!
! end the simulation and exit MPI
!
!-------------------------------------------------------------------------------------------------

! version with rank number printed in the error message
  subroutine exit_MPI(myrank,error_msg)

  use constants
  use shared_input_parameters, only: OUTPUT_FILES

  implicit none

  ! identifier for error message file
  integer, parameter :: IERROR = 30

  integer :: myrank
  character(len=*) :: error_msg
  character(len=80) :: outputname

  ! write error message to screen
  write(*,*) error_msg(1:len(error_msg))
  write(*,*) 'Error detected, aborting MPI... proc ',myrank

  ! write error message to file
  write(outputname,"('/error_message',i6.6,'.txt')") myrank
  open(unit=IERROR,file=trim(OUTPUT_FILES)//'/'//outputname,status='unknown')
  write(IERROR,*) error_msg(1:len(error_msg))
  write(IERROR,*) 'Error detected, aborting MPI... proc ',myrank
  close(IERROR)

  ! close output file
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) close(IMAIN)

  ! stop all the MPI processes, and exit
  call abort_mpi()

  ! otherwise: there is no standard behaviour to exit with an error code in Fortran,
  ! however most compilers do recognize this as an error code stop statement;
  ! to check stop code in terminal: > echo $?
  stop 30

  ! or just exit with message:
  !stop 'Error, program ended in exit_MPI'

  end subroutine exit_MPI

!
!-------------------------------------------------------------------------------------------------
!

! version without rank number printed in the error message
  subroutine exit_MPI_without_rank(error_msg)

  use constants

  implicit none

  character(len=*) :: error_msg

  ! write error message to screen
  write(*,*) error_msg(1:len(error_msg))
  write(*,*) 'Error detected, aborting MPI...'

  ! stop all the MPI processes, and exit
  call abort_mpi()

  stop 'Error, program ended in exit_MPI'

  end subroutine exit_MPI_without_rank

