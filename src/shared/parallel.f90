!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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
! MPI wrapper functions
!
!-------------------------------------------------------------------------------------------------

  subroutine init_mpi()

  use mpi

  implicit none

  integer :: ier

  call MPI_INIT(ier)
  if (ier /= 0 ) stop 'Error initializing MPI'

  end subroutine init_mpi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine finalize_mpi()

  use mpi

  implicit none

  integer :: ier

  call MPI_FINALIZE(ier)
  if (ier /= 0 ) stop 'Error finalizing MPI'

  end subroutine finalize_mpi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine abort_mpi()

  use mpi

  implicit none

  integer :: ier

  ! note: MPI_ABORT does not return, and does exit the
  !          program with an error code of 30
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)

  end subroutine abort_mpi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine synchronize_all()

  use mpi

  implicit none

  integer :: ier

  ! synchronizes MPI processes
  call MPI_BARRIER(MPI_COMM_WORLD,ier)
  if (ier /= 0 ) stop 'Error synchronize MPI processes'

  end subroutine synchronize_all

!
!-------------------------------------------------------------------------------------------------
!

  subroutine synchronize_all_comm(comm)

  use mpi

  implicit none

  integer,intent(in) :: comm

  ! local parameters
  integer :: ier

  ! synchronizes MPI processes
  call MPI_BARRIER(comm,ier)
  if (ier /= 0 ) stop 'Error synchronize MPI processes for specified communicator'

  end subroutine synchronize_all_comm

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

  integer :: ier

  call MPI_TEST(request,flag_result_test,MPI_STATUS_IGNORE,ier)

  end subroutine test_request

!
!-------------------------------------------------------------------------------------------------
!

  subroutine irecv_cr(recvbuf, recvcount, dest, recvtag, req)

  use constants
  use mpi

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

  use constants
  use mpi

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

  integer :: ier

  call mpi_wait(req,MPI_STATUS_IGNORE,ier)

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

  use constants
  use mpi

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

  subroutine max_allreduce_i(buffer,countval)

  use mpi

  implicit none

  integer :: countval
  integer,dimension(countval),intent(inout) :: buffer

  ! local parameters
  integer :: ier
  integer,dimension(countval) :: send

  ! seems not to be supported on all kind of MPI implementations...
  !! DK DK: yes, I confirm, using MPI_IN_PLACE is tricky
  !! DK DK (see the answer at http://stackoverflow.com/questions/17741574/in-place-mpi-reduce-crashes-with-openmpi
  !! DK DK      for how to use it right)
  !call MPI_ALLREDUCE(MPI_IN_PLACE, buffer, countval, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ier)

  send(:) = buffer(:)

  call MPI_ALLREDUCE(send, buffer, countval, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ier)
  if (ier /= 0 ) stop 'Allreduce to get max values failed.'

  end subroutine max_allreduce_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine max_all_cr(sendbuf, recvbuf)

  use constants
  use mpi

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

  subroutine sum_all_3Darray_dp(sendbuf, recvbuf, nx,ny,nz)

  use mpi

  implicit none

  integer :: nx,ny,nz
  double precision, dimension(nx,ny,nz) :: sendbuf, recvbuf
  integer :: ier

  ! this works only if the arrays are contiguous in memory (which is always the case for static arrays, as used in the code)
  call MPI_REDUCE(sendbuf,recvbuf,nx*ny*nz,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ier)

  end subroutine sum_all_3Darray_dp

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

  subroutine bcast_all_i(buffer, countval)

  use mpi

  implicit none

  integer :: countval
  integer, dimension(countval) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,countval,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_cr(buffer, countval)

  use constants
  use mpi

  implicit none

  include "precision.h"

  integer :: countval
  real(kind=CUSTOM_REAL), dimension(countval) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,countval,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_cr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_r(buffer, countval)

  use mpi

  implicit none

  integer :: countval
  real, dimension(countval) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,countval,MPI_REAL,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_r

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singler(buffer)

  use mpi

  implicit none

  real :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,MPI_REAL,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_singler

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_dp(buffer, countval)

  use mpi

  implicit none

  integer :: countval
  double precision, dimension(countval) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,countval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_singledp(buffer)

  use mpi

  implicit none

  double precision :: buffer

  integer :: ier

  call MPI_BCAST(buffer,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_singledp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_ch(buffer, countval)

  use mpi

  implicit none

  integer :: countval
  character(len=countval) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,countval,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_ch

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_ch_array(buffer,ndim,countval)

  use mpi

  implicit none

  integer :: countval,ndim
  character(len=countval),dimension(ndim) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,ndim*countval,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_ch_array

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_ch_array2(buffer,ndim1,ndim2,countval)

  use mpi

  implicit none

  integer :: countval,ndim1,ndim2
  character(len=countval),dimension(ndim1,ndim2) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,ndim1*ndim2*countval,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_ch_array2

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_l(buffer, countval)

  use mpi

  implicit none

  integer :: countval
  logical,dimension(countval) :: buffer

  integer :: ier

  call MPI_BCAST(buffer,countval,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)

  end subroutine bcast_all_l

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_singlei(recvbuf, dest, recvtag)

  use mpi

  implicit none

  integer :: dest,recvtag
  integer :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,1,MPI_INTEGER,dest,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

  end subroutine recv_singlei

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_singlel(recvbuf, dest, recvtag)

  use mpi

  implicit none

  integer :: dest,recvtag
  logical :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,1,MPI_LOGICAL,dest,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

  end subroutine recv_singlel

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_i(recvbuf, recvcount, dest, recvtag)

  use mpi

  implicit none

  integer :: dest,recvtag
  integer :: recvcount
  integer,dimension(recvcount) :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,recvcount,MPI_INTEGER,dest,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

  end subroutine recv_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_cr(recvbuf, recvcount, dest, recvtag)

  use constants
  use mpi

  implicit none

  include "precision.h"

  integer :: dest,recvtag
  integer :: recvcount
  real(kind=CUSTOM_REAL),dimension(recvcount) :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,recvcount,CUSTOM_MPI_TYPE,dest,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

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

  integer :: ier

  call MPI_RECV(recvbuf,recvcount,MPI_DOUBLE_PRECISION,dest,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

  end subroutine recv_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine recv_ch(recvbuf, recvcount, dest, recvtag)

  use mpi

  implicit none

  integer :: dest,recvtag
  integer :: recvcount
  character(len=recvcount) :: recvbuf

  integer :: ier

  call MPI_RECV(recvbuf,recvcount,MPI_CHARACTER,dest,recvtag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

  end subroutine recv_ch

!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_ch(sendbuf, sendcount, dest, sendtag)

  use mpi

  implicit none

  integer :: dest,sendtag
  integer :: sendcount
  character(len=sendcount) :: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf,sendcount,MPI_CHARACTER,dest,sendtag,MPI_COMM_WORLD,ier)

  end subroutine send_ch


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

  subroutine send_singlel(sendbuf, dest, sendtag)

  use mpi

  implicit none

  integer :: dest,sendtag
  logical :: sendbuf

  integer :: ier

  call MPI_SEND(sendbuf,1,MPI_LOGICAL,dest,sendtag,MPI_COMM_WORLD,ier)

  end subroutine send_singlel

!
!-------------------------------------------------------------------------------------------------
!

  subroutine send_cr(sendbuf, sendcount, dest, sendtag)

  use constants
  use mpi

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

  use constants
  use mpi

  implicit none

  include "precision.h"

  integer :: sendcount, recvcount, dest, sendtag, source, recvtag
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

  integer :: ier

  call MPI_SENDRECV(sendbuf,sendcount,CUSTOM_MPI_TYPE,dest,sendtag, &
                    recvbuf,recvcount,CUSTOM_MPI_TYPE,source,recvtag, &
                    MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

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

  integer :: ier

  call MPI_SENDRECV(sendbuf,sendcount,MPI_DOUBLE_PRECISION,dest,sendtag, &
                    recvbuf,recvcount,MPI_DOUBLE_PRECISION,source,recvtag, &
                    MPI_COMM_WORLD,MPI_STATUS_IGNORE,ier)

  end subroutine sendrecv_dp

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

  subroutine gather_all_singlei(sendbuf, recvbuf, NPROC)

  use mpi

  implicit none

  integer :: NPROC
  integer :: sendbuf
  integer, dimension(0:NPROC-1) :: recvbuf

  integer :: ier

  call MPI_GATHER(sendbuf,1,MPI_INTEGER, &
                  recvbuf,1,MPI_INTEGER, &
                  0,MPI_COMM_WORLD,ier)

  end subroutine gather_all_singlei

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gather_all_cr(sendbuf, sendcnt, recvbuf, recvcount, NPROC)

  use constants
  use mpi

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

  use constants
  use mpi

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

  use constants
  use mpi

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

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gatherv_all_r(sendbuf, sendcnt, recvbuf, recvcount, recvoffset,recvcounttot, NPROC)

  use constants
  use mpi

  implicit none

  integer :: sendcnt,recvcounttot,NPROC
  integer, dimension(NPROC) :: recvcount,recvoffset
  real, dimension(sendcnt) :: sendbuf
  real, dimension(recvcounttot) :: recvbuf

  integer :: ier

  call MPI_GATHERV(sendbuf,sendcnt,MPI_REAL, &
                  recvbuf,recvcount,recvoffset,MPI_REAL, &
                  0,MPI_COMM_WORLD,ier)

  end subroutine gatherv_all_r

!
!-------------------------------------------------------------------------------------------------
!

  subroutine scatter_all_singlei(sendbuf, recvbuf, NPROC)

  use mpi

  implicit none

  integer :: NPROC
  integer, dimension(0:NPROC-1) :: sendbuf
  integer :: recvbuf

  integer :: ier

  call MPI_Scatter(sendbuf, 1, MPI_INTEGER, &
                   recvbuf, 1, MPI_INTEGER, &
                   0, MPI_COMM_WORLD, ier)

  end subroutine scatter_all_singlei

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_size(sizeval)

  use mpi

  implicit none

  integer,intent(out) :: sizeval

  ! local parameters
  integer :: ier

  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeval,ier)
  if (ier /= 0 ) stop 'Error getting MPI world size'

  end subroutine world_size

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_rank(rank)

  use mpi

  implicit none

  integer,intent(out) :: rank

  ! local parameters
  integer :: ier

  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ier)
  if (ier /= 0 ) stop 'Error getting MPI rank'

  end subroutine world_rank

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_duplicate(comm)

  use mpi

  implicit none

  integer,intent(out) :: comm
  integer :: ier

  call MPI_COMM_DUP(MPI_COMM_WORLD,comm,ier)
  if (ier /= 0 ) stop 'Error duplicating MPI_COMM_WORLD communicator'

  end subroutine world_duplicate

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_get_comm(comm)

  use mpi

  implicit none

  integer,intent(out) :: comm

  comm = MPI_COMM_WORLD

  end subroutine world_get_comm

!
!-------------------------------------------------------------------------------------------------
!

  subroutine world_get_comm_self(comm)

  use mpi

  implicit none

  integer,intent(out) :: comm

  comm = MPI_COMM_SELF

  end subroutine world_get_comm_self

