!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            August 2013
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

  use mpi
  use constants

  implicit none

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

  ! otherwise: there is no standard behaviour to exit with an error code in Fortran,
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

  use mpi
  use constants

  implicit none

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
!
! MPI wrapper functions
!
!-------------------------------------------------------------------------------------------------

  subroutine init_mpi()

  use mpi

  implicit none

  integer :: ier

  call MPI_INIT(ier)
  if( ier /= 0 ) stop 'error initializing MPI'

  end subroutine init_mpi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine finalize_mpi()

  use mpi

  implicit none

  integer :: ier

  call MPI_FINALIZE(ier)
  if( ier /= 0 ) stop 'error finalizing MPI'

  end subroutine finalize_mpi



!
!-------------------------------------------------------------------------------------------------
!

  subroutine sync_all()

  use mpi

  implicit none

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

  integer function null_process()

  use mpi

  implicit none

  null_process = MPI_PROC_NULL

  end function null_process

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_request(request,flag_result_test)

  use mpi

  implicit none

  integer :: request
  logical :: flag_result_test

  ! MPI status of messages to be received
  integer :: msg_status(MPI_STATUS_SIZE)

  integer :: ier

  call MPI_TEST(request,flag_result_test,msg_status,ier)

  end subroutine test_request
!
!-------------------------------------------------------------------------------------------------
!

  subroutine irecv_cr(recvbuf, recvcount, dest, recvtag, req)

  use mpi
  use constants

  implicit none

  include "precision.h"

  integer :: recvcount, dest, recvtag, req
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

  integer ier

  call MPI_IRECV(recvbuf(1),recvcount,CUSTOM_MPI_TYPE,dest,recvtag, &
                  MPI_COMM_WORLD,req,ier)

  end subroutine irecv_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine irecv_dp(recvbuf, recvcount, dest, recvtag, req)

  use mpi

  implicit none

  integer :: recvcount, dest, recvtag, req
  double precision, dimension(recvcount) :: recvbuf

  integer :: ier

  call MPI_IRECV(recvbuf(1),recvcount,MPI_DOUBLE_PRECISION,dest,recvtag, &
                  MPI_COMM_WORLD,req,ier)

  end subroutine irecv_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine isend_cr(sendbuf, sendcount, dest, sendtag, req)

  use mpi
  use constants

  implicit none

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

  subroutine isend_dp(sendbuf, sendcount, dest, sendtag, req)

  use mpi

  implicit none

  integer :: sendcount, dest, sendtag, req
  double precision, dimension(sendcount) :: sendbuf

  integer :: ier

  call MPI_ISEND(sendbuf(1),sendcount,MPI_DOUBLE_PRECISION,dest,sendtag, &
                  MPI_COMM_WORLD,req,ier)

  end subroutine isend_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine wait_req(req)

  use mpi

  implicit none

  integer :: req

  integer, dimension(MPI_STATUS_SIZE) :: req_mpi_status

  integer :: ier

  call mpi_wait(req,req_mpi_status,ier)

  end subroutine wait_req

!
!-------------------------------------------------------------------------------------------------
!

  double precision function wtime()

  use mpi

  implicit none

  wtime = MPI_WTIME()

  end function wtime

!
!-------------------------------------------------------------------------------------------------
!

  subroutine min_all_i(sendbuf, recvbuf)

  use mpi

  implicit none

  integer:: sendbuf, recvbuf
  integer ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_MIN,0,MPI_COMM_WORLD,ier)

  end subroutine min_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine min_all_cr(sendbuf, recvbuf)

  use mpi
  use constants

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ier)

  end subroutine min_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_i(sendbuf, recvbuf)

  use mpi

  implicit none

  integer :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ier)

  end subroutine max_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_cr(sendbuf, recvbuf)

  use mpi
  use constants

  implicit none

  include "precision.h"

  real(kind=CUSTOM_REAL) :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,CUSTOM_MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ier)

  end subroutine max_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_i(sendbuf, recvbuf)

  use mpi

  implicit none

  integer :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ier)

  end subroutine sum_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sum_all_dp(sendbuf, recvbuf)

  use mpi

  implicit none

  double precision :: sendbuf, recvbuf
  integer :: ier

  call MPI_REDUCE(sendbuf,recvbuf,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ier)

  end subroutine sum_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_iproc_i(buffer,iproc)

  use mpi

  implicit none

  integer :: iproc
  integer :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,MPI_INTEGER,iproc,MPI_COMM_WORLD,ier)

  end subroutine bcast_iproc_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singlei(buffer)

  use mpi

  implicit none

  integer :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_singlei

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_i(buffer, count)

  use :: mpi

  implicit none

  integer :: count
  integer, dimension(count) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,count,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_cr(buffer, count)

  use :: mpi
  use constants

  implicit none

  include "precision.h"

  integer :: count
  real(kind=CUSTOM_REAL), dimension(count) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,count,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_r(buffer, count)

  use :: mpi
  use constants

  implicit none

  include "precision.h"

  integer :: count
  real, dimension(count) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,count,MPI_REAL,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_r

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_dp(buffer, count)

  use :: mpi

  implicit none

  integer :: count
  double precision, dimension(count) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_dp


!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_ch(buffer, count)

  use :: mpi

  implicit none

  integer :: count
  character(len=count) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,count,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_ch

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_l(buffer, count)

  use :: mpi

  implicit none

  integer :: count
  logical,dimension(count) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,count,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_l


!
!-------------------------------------------------------------------------------------------------
!


  subroutine recv_singlei(recvbuf, dest, recvtag)

  use mpi

  implicit none

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

  use mpi

  implicit none

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


  subroutine recv_cr(recvbuf, recvcount, dest, recvtag)

  use mpi
  use constants

  implicit none

  include "precision.h"

  integer :: dest,recvtag
  integer :: recvcount
  real(kind=CUSTOM_REAL),dimension(recvcount) :: recvbuf

  ! MPI status of messages to be received
  integer :: msg_status(MPI_STATUS_SIZE)
  integer :: ier

  call MPI_RECV(recvbuf,recvcount,CUSTOM_MPI_TYPE,dest,recvtag,MPI_COMM_WORLD,msg_status,ier)

  end subroutine recv_cr

!
!-------------------------------------------------------------------------------------------------
!


  subroutine recv_dp(recvbuf, recvcount, dest, recvtag)

  use mpi

  implicit none

  integer :: dest,recvtag
  integer :: recvcount
  double precision,dimension(recvcount) :: recvbuf

  ! MPI status of messages to be received
  integer :: msg_status(MPI_STATUS_SIZE)
  integer :: ier

  call MPI_RECV(recvbuf,recvcount,MPI_DOUBLE_PRECISION,dest,recvtag,MPI_COMM_WORLD,msg_status,ier)

  end subroutine recv_dp


!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_i(sendbuf, sendcount, dest, sendtag)

  use mpi

  implicit none

  integer :: dest,sendtag
  integer :: sendcount
  integer,dimension(sendcount):: sendbuf
  integer :: ier

  call MPI_SEND(sendbuf,sendcount,MPI_INTEGER,dest,sendtag,MPI_COMM_WORLD,ier)

  end subroutine send_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_singlei(sendbuf, dest, sendtag)

  use mpi

  implicit none

  integer :: dest,sendtag
  integer :: sendbuf
  integer :: ier

  call MPI_SEND(sendbuf,1,MPI_INTEGER,dest,sendtag,MPI_COMM_WORLD,ier)

  end subroutine send_singlei


!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_cr(sendbuf, sendcount, dest, sendtag)

  use mpi
  use constants

  implicit none

  include "precision.h"

  integer :: dest,sendtag
  integer :: sendcount
  real(kind=CUSTOM_REAL),dimension(sendcount):: sendbuf
  integer :: ier

  call MPI_SEND(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag,MPI_COMM_WORLD,ier)

  end subroutine send_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_dp(sendbuf, sendcount, dest, sendtag)

  use mpi

  implicit none

  integer :: dest,sendtag
  integer :: sendcount
  double precision,dimension(sendcount):: sendbuf
  integer :: ier

  call MPI_SEND(sendbuf,sendcount,MPI_DOUBLE_PRECISION,dest,sendtag,MPI_COMM_WORLD,ier)

  end subroutine send_dp


!
!-------------------------------------------------------------------------------------------------
!

  subroutine sendrecv_cr(sendbuf, sendcount, dest, sendtag, &
                         recvbuf, recvcount, source, recvtag)

  use mpi
  use constants

  implicit none

  include "precision.h"

  integer :: sendcount, recvcount, dest, sendtag, source, recvtag
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

! MPI status of messages to be received
  integer :: msg_status(MPI_STATUS_SIZE)

  integer :: ier

  call MPI_SENDRECV(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag, &
                    recvbuf,recvcount,CUSTOM_MPI_TYPE,source,recvtag, &
                    MPI_COMM_WORLD,msg_status,ier)

  end subroutine sendrecv_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sendrecv_dp(sendbuf, sendcount, dest, sendtag, &
                         recvbuf, recvcount, source, recvtag)

  use mpi

  implicit none

  integer :: sendcount, recvcount, dest, sendtag, source, recvtag
  double precision, dimension(sendcount) :: sendbuf
  double precision, dimension(recvcount) :: recvbuf

! MPI status of messages to be received
  integer :: msg_status(MPI_STATUS_SIZE)

  integer :: ier

  call MPI_SENDRECV(sendbuf,sendcount,MPI_DOUBLE_PRECISION,dest,sendtag, &
                    recvbuf,recvcount,MPI_DOUBLE_PRECISION,source,recvtag, &
                    MPI_COMM_WORLD,msg_status,ier)

  end subroutine sendrecv_dp

!
!-------------------------------------------------------------------------------------------------
!



  subroutine world_size(size)

  use mpi

  implicit none

  integer :: size
  integer :: ier

  call MPI_COMM_SIZE(MPI_COMM_WORLD,size,ier)

  end subroutine world_size

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_rank(rank)

  use mpi

  implicit none

  integer :: rank
  integer :: ier

  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ier)

  end subroutine world_rank

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_i(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use mpi

  implicit none

  integer :: sendcnt, recvcount, NPROC
  integer, dimension(sendcnt) :: sendbuf
  integer, dimension(recvcount,0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_GATHER(sendbuf,sendcnt,MPI_INTEGER, &
                  recvbuf,recvcount,MPI_INTEGER, &
                  0,MPI_COMM_WORLD,ier)

  end subroutine gather_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_cr(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use mpi
  use constants

  implicit none

  include "precision.h"

  integer :: sendcnt, recvcount, NPROC
  real(kind=CUSTOM_REAL), dimension(sendcnt) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcount,0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_GATHER(sendbuf,sendcnt,CUSTOM_MPI_TYPE, &
                  recvbuf,recvcount,CUSTOM_MPI_TYPE, &
                  0,MPI_COMM_WORLD,ier)

  end subroutine gather_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_dp(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use mpi

  implicit none

  integer :: sendcnt, recvcount, NPROC
  double precision, dimension(sendcnt) :: sendbuf
  double precision, dimension(recvcount,0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_GATHER(sendbuf,sendcnt,MPI_DOUBLE_PRECISION, &
                  recvbuf,recvcount,MPI_DOUBLE_PRECISION, &
                  0,MPI_COMM_WORLD,ier)

  end subroutine gather_all_dp


!
!-------------------------------------------------------------------------------------------------
!

  subroutine gatherv_all_i(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)

  use mpi
  use constants

  implicit none

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

  use mpi
  use constants

  implicit none

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

