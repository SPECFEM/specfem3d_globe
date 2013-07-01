!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            April 2011
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

! end the simulation and exit MPI

! version with rank number printed in the error message
  subroutine exit_MPI(myrank,error_msg)

  implicit none

  ! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"

  ! identifier for error message file
  integer, parameter :: IERROR = 30

  integer :: myrank
  character(len=*) error_msg

  integer :: ier
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
  ! note: MPI_ABORT does not return, and does exit the
  !          program with an error code of 30
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)

  ! otherwise: there is no standard behaviour to exit with an error code in fortran,
  ! however most compilers do recognize this as an error code stop statement;
  ! to check stop code in terminal: > echo $?
  stop 30

  ! or just exit with message:
  !stop 'error, program ended in exit_MPI'

  end subroutine exit_MPI

!
!-------------------------------------------------------------------------------------------------
!

! version without rank number printed in the error message
  subroutine exit_MPI_without_rank(error_msg)

  implicit none

  ! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"

  character(len=*) error_msg

  integer :: ier

  ! write error message to screen
  write(*,*) error_msg(1:len(error_msg))
  write(*,*) 'Error detected, aborting MPI...'

  ! stop all the MPI processes, and exit
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)
  stop 'error, program ended in exit_MPI'

  end subroutine exit_MPI_without_rank

!-------------------------------------------------------------------------------------------------
!
! I/O wrapper function
!
!-------------------------------------------------------------------------------------------------

  subroutine flush_IMAIN()

  implicit none

  include "constants.h"

  ! only master process writes out to main output file
  ! file I/O in fortran is buffered by default
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
!
! MPI wrapper functions
!
!-------------------------------------------------------------------------------------------------


  subroutine sync_all()

  implicit none

  ! standard include of the MPI library
  include 'mpif.h'

  integer :: ier,rank

  ! gets callers rank
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ier)

  ! synchronizes MPI processes
  call MPI_BARRIER(MPI_COMM_WORLD,ier)
  if( ier /= 0 ) call exit_mpi(rank,'error synchronize MPI processes')

  end subroutine sync_all

!
!-------------------------------------------------------------------------------------------------
!

  subroutine irecv_cr(recvbuf, recvcount, dest, recvtag, req)

  implicit none

  ! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer recvcount, dest, recvtag, req
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

  integer ier

  call MPI_IRECV(recvbuf(1),recvcount,CUSTOM_MPI_TYPE,dest,recvtag, &
                  MPI_COMM_WORLD,req,ier)

  end subroutine irecv_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine isend_cr(sendbuf, sendcount, dest, sendtag, req)

  implicit none

  ! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer sendcount, dest, sendtag, req
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf

  integer ier

  call MPI_ISEND(sendbuf(1),sendcount,CUSTOM_MPI_TYPE,dest,sendtag, &
                  MPI_COMM_WORLD,req,ier)

  end subroutine isend_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine wait_req(req)

  implicit none

  ! standard include of the MPI library
  include 'mpif.h'

  integer :: req

  integer, dimension(MPI_STATUS_SIZE) :: req_mpi_status

  integer :: ier

  call mpi_wait(req,req_mpi_status,ier)

  end subroutine wait_req

!
!-------------------------------------------------------------------------------------------------
!

  double precision function wtime()

  implicit none

  ! standard include of the MPI library
  include 'mpif.h'

  wtime = MPI_WTIME()

  end function wtime

!
!-------------------------------------------------------------------------------------------------
!

  subroutine min_all_i(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer:: sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER, &
                  MPI_MIN,0,MPI_COMM_WORLD,ier)

  end subroutine min_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_i(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER, &
                  MPI_MAX,0,MPI_COMM_WORLD,ier)

  end subroutine max_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_cr(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE, &
                  MPI_MAX,0,MPI_COMM_WORLD,ier)

  end subroutine max_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_i(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER, &
                  MPI_SUM,0,MPI_COMM_WORLD,ier)

  end subroutine sum_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_dp(sendbuf, recvbuf)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  double precision sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION, &
                  MPI_SUM,0,MPI_COMM_WORLD,ier)

  end subroutine sum_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_iproc_i(buffer,iproc)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer :: iproc
  integer :: buffer

  integer ier

  call MPI_BCAST(buffer,1,MPI_INTEGER,iproc,MPI_COMM_WORLD,ier)

  end subroutine bcast_iproc_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singlei(buffer)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer :: buffer

  integer ier

  call MPI_BCAST(buffer,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_singlei

!
!-------------------------------------------------------------------------------------------------
!


  subroutine recv_singlei(recvbuf, dest, recvtag)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer :: dest,recvtag
  integer :: recvbuf

  ! MPI status of messages to be received
  integer :: msg_status(MPI_STATUS_SIZE)
  integer :: ier

  call MPI_RECV(recvbuf,1,MPI_INTEGER,dest,recvtag,MPI_COMM_WORLD,msg_status,ier)

  end subroutine recv_singlei

!
!-------------------------------------------------------------------------------------------------
!


  subroutine recv_i(recvbuf, recvcount, dest, recvtag)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer :: dest,recvtag
  integer :: recvcount
  integer,dimension(recvcount) :: recvbuf

  ! MPI status of messages to be received
  integer :: msg_status(MPI_STATUS_SIZE)
  integer :: ier

  call MPI_RECV(recvbuf,recvcount,MPI_INTEGER,dest,recvtag,MPI_COMM_WORLD,msg_status,ier)

  end subroutine recv_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_i(sendbuf, sendcount, dest, sendtag)

  implicit none

  ! standard include of the MPI library
  include 'mpif.h'

  !integer sendbuf,sendcount,dest,sendtag
  integer dest,sendtag
  integer sendcount
  integer,dimension(sendcount):: sendbuf
  integer ier

  call MPI_SEND(sendbuf,sendcount,MPI_INTEGER,dest,sendtag,MPI_COMM_WORLD,ier)

  end subroutine send_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_singlei(sendbuf, dest, sendtag)

  implicit none

  ! standard include of the MPI library
  include 'mpif.h'

  !integer sendbuf,sendcount,dest,sendtag
  integer :: dest,sendtag
  integer :: sendbuf
  integer :: ier

  call MPI_SEND(sendbuf,1,MPI_INTEGER,dest,sendtag,MPI_COMM_WORLD,ier)

  end subroutine send_singlei

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_size(size)

  implicit none

  ! standard include of the MPI library
  include 'mpif.h'

  integer size
  integer ier

  call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ier)

  end subroutine world_size

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_rank(rank)

  implicit none

  ! standard include of the MPI library
  include 'mpif.h'

  integer rank
  integer ier

  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ier)

  end subroutine world_rank

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_i(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  implicit none

  ! standard include of the MPI library
  include 'mpif.h'

  integer sendcnt, recvcount, NPROC
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcount,0:NPROC-1) :: recvbuf

  integer ier

  call MPI_GATHER(sendbuf,sendcnt,MPI_INTEGER, &
                  recvbuf,recvcount,MPI_INTEGER, &
                  0,MPI_COMM_WORLD,ier)

  end subroutine gather_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gatherv_all_i(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)

  implicit none

  ! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer :: sendcnt,recvcounttot,NPROC
  integer, dimension(NPROC) :: recvcount,recvoffset
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcounttot) :: recvbuf

  integer :: ier

  call MPI_GATHERV(sendbuf,sendcnt,MPI_INTEGER, &
                  recvbuf,recvcount,recvoffset,MPI_INTEGER, &
                  0,MPI_COMM_WORLD,ier)

  end subroutine gatherv_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gatherv_all_cr(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)

  implicit none

  ! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer :: sendcnt,recvcounttot,NPROC
  integer, dimension(NPROC) :: recvcount,recvoffset
  real(kind=CUSTOM_REAL), dimension(sendcnt) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcounttot) :: recvbuf

  integer :: ier

  call MPI_GATHERV(sendbuf,sendcnt,CUSTOM_MPI_TYPE, &
                  recvbuf,recvcount,recvoffset,CUSTOM_MPI_TYPE, &
                  0,MPI_COMM_WORLD,ier)

  end subroutine gatherv_all_cr

