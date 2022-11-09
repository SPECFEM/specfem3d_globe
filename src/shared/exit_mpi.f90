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
! end the simulation and exit MPI
!
!-------------------------------------------------------------------------------------------------

! version with rank number printed in the error message
  subroutine exit_MPI(myrank,error_msg)

  use constants, only: IMAIN,ISTANDARD_OUTPUT
  use shared_input_parameters, only: OUTPUT_FILES

  implicit none

  ! identifier for error message file
  integer, parameter :: IERROR = 30

  integer,intent(in) :: myrank
  character(len=*),intent(in) :: error_msg

  ! local parameters
  character(len=80) :: outputname

  ! write error message to screen
  write(*,*) error_msg(1:len(error_msg))
  write(*,*) 'Error detected, aborting MPI... proc ',myrank

  ! write error message to file
  write(outputname,"('/error_message',i6.6,'.txt')") myrank
  open(unit=IERROR,file=trim(OUTPUT_FILES)//'/'//trim(outputname),status='unknown')
  write(IERROR,*) error_msg(1:len(error_msg))
  write(IERROR,*) 'Error detected, aborting MPI... proc ',myrank
  close(IERROR)

  ! close output file
  if (myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) close(IMAIN)

  ! flushes possible left-overs from print-statements
  call flush_stdout()

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

  ! flushes possible left-overs from print-statements
  call flush_stdout()

  ! stop all the MPI processes, and exit
  call abort_mpi()

  stop 'Error, program ended in exit_MPI'

  end subroutine exit_MPI_without_rank

!-------------------------------------------------------------------------------------------------
!
! shared helper functions
!
!-------------------------------------------------------------------------------------------------
! put here for lack of better places... move to a better file in future

  subroutine print_gll_min_max_all(nspec,array,name)

! prints out minimum/maximum value of given GLL array

  use constants, only: CUSTOM_REAL,IMAIN,NGLLX,NGLLY,NGLLZ,myrank

  implicit none

  integer, intent(in) :: nspec
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: array
  character(len=*),intent(in) :: name

  ! local parameters
  real(kind=CUSTOM_REAL) :: minvalue,maxvalue,min_all,max_all

  ! gets min/max for all slices
  maxvalue = maxval( array )
  minvalue = minval( array )
  call max_all_cr(maxvalue, max_all)
  call min_all_cr(minvalue, min_all)

  if (myrank == 0) then
    write(IMAIN,*) '  '//trim(name)//' min/max: ',min_all,max_all
  endif

  end subroutine print_gll_min_max_all
