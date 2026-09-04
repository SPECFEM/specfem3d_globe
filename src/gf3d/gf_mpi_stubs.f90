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

!----
!---- Serial replacements for the handful of MPI wrappers that reused
!---- src/shared/ objects reference.
!----
!---- Why this file exists
!---- --------------------
!---- libgf3d.a is a serial library: it must link without an MPI runtime,
!---- and it must never abort the process of a caller that loaded it as a
!---- shared object.
!----
!---- The routines it reuses from src/shared/ do not come in MPI-free
!---- variants, because a static archive pulls in whole object files:
!---- $O/model_topo_bathy.shared.o holds get_topo_bathy() (which we need,
!---- and which is itself MPI-free) together with
!---- model_topo_bathy_broadcast() (which we never call, but which
!---- references wtime() and bcast_all_i()) and six error paths that call
!---- exit_MPI().
!----
!---- These stubs satisfy the linker for those references. They are
!---- deliberately NOT a reimplementation of the parallel code: the
!---- routines that would use them are unreachable in a serial program.
!----
!---- Note this file intentionally replaces exit_MPI() rather than
!---- editing the src/shared/ sources to say "stop". The precedent in
!---- model_prem.f90 does not transfer: every exit_MPI() call in
!---- model_topo_bathy.f90 sits inside an `if (myrank == 0)` block that is
!---- followed by a collective bcast_all_i(), so turning them into `stop`
!---- would replace a clean MPI_ABORT of the whole job with a hang of
!---- ranks 1..n-1. Stubbing here leaves the solver bit-for-bit unchanged.
!----
!---- No target links both this file and $O/parallel.sharedmpi.o /
!---- $O/exit_mpi.shared.o, and none should.
!----

!
!-------------------------------------------------------------------------------------------------
!

  double precision function wtime()

! serial replacement for MPI_WTIME()

  implicit none

  integer(kind=8) :: count,count_rate

  call system_clock(count,count_rate)

  if (count_rate > 0) then
    wtime = dble(count) / dble(count_rate)
  else
    wtime = 0.d0
  endif

  end function wtime

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_all_i(buffer, countval)

! serial replacement for MPI_BCAST() of an integer array
!
! with a single process the broadcast is the identity, so there is nothing to do

  implicit none

  integer :: countval
  integer, dimension(countval) :: buffer

  ! nothing to broadcast to; touches the arguments so compilers do not warn
  if (countval < 0) buffer(1) = 0

  end subroutine bcast_all_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine exit_MPI(myrank,error_msg)

! serial replacement for src/shared/exit_mpi.f90
!
! reports on stderr and stops. Note this must not be reachable from the
! library's own code paths: every routine in src/gf3d/ returns an error
! status instead, precisely so that a Python interpreter holding
! libgf3d.so is never killed from underneath the user.

  implicit none

  integer,intent(in) :: myrank
  character(len=*),intent(in) :: error_msg

  ! local parameters
  integer, parameter :: ISTDERR = 0

  write(ISTDERR,*) 'Error: ',error_msg(1:len(error_msg))
  write(ISTDERR,*) 'Error detected in serial Green function library, rank ',myrank

  stop 30

  end subroutine exit_MPI

!
!-------------------------------------------------------------------------------------------------
!

  subroutine exit_MPI_without_rank(error_msg)

! serial replacement for src/shared/exit_mpi.f90

  implicit none

  character(len=*),intent(in) :: error_msg

  ! local parameters
  integer, parameter :: ISTDERR = 0

  write(ISTDERR,*) 'Error: ',error_msg(1:len(error_msg))

  stop 30

  end subroutine exit_MPI_without_rank
