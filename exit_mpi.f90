!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                    and University of Pau, France
! (c) California Institute of Technology and University of Pau, April 2007
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

! end the simulation and exit MPI

! version with rank number printed in the error message
  subroutine exit_MPI(myrank,error_msg)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"

! identifier for error message file
  integer, parameter :: IERROR = 30

  integer myrank
  character(len=*) error_msg

  integer ier
  character(len=80) outputname
  character(len=150) OUTPUT_FILES

! write error message to screen
  write(*,*) error_msg(1:len(error_msg))
  write(*,*) 'Error detected, aborting MPI... proc ',myrank

! write error message to file
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')
  write(outputname,"('/error_message',i6.6,'.txt')") myrank
  open(unit=IERROR,file=trim(OUTPUT_FILES)//outputname,status='unknown')
  write(IERROR,*) error_msg(1:len(error_msg))
  write(IERROR,*) 'Error detected, aborting MPI... proc ',myrank
  close(IERROR)

! close output file
  if(myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) close(IMAIN)

! stop all the MPI processes, and exit
  call MPI_ABORT(ier)
  stop 'error, program ended in exit_MPI'

  end subroutine exit_MPI

!
!----
!

! version without rank number printed in the error message
  subroutine exit_MPI_without_rank(error_msg)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"

  character(len=*) error_msg

  integer ier

! write error message to screen
  write(*,*) error_msg(1:len(error_msg))
  write(*,*) 'Error detected, aborting MPI...'

! stop all the MPI processes, and exit
  call MPI_ABORT(ier)
  stop 'error, program ended in exit_MPI'

  end subroutine exit_MPI_without_rank

