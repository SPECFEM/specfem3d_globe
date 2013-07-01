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

!> @brief Initialize ADIOS and setup the xml output file
subroutine adios_setup()
  use adios_write_mod, only: adios_init

  implicit none
  integer :: adios_err, sizeMB

  call adios_init_noxml (adios_err);
  sizeMB = 200 ! TODO 200MB is surely not the right size for the adios buffer
  call adios_allocate_buffer (sizeMB , adios_err)
end subroutine adios_setup

!> @brief Finalize ADIOS. Must be called once everything is written down.
subroutine adios_cleanup()
  use mpi
  use adios_write_mod, only: adios_finalize

  implicit none
  integer :: myrank
  integer :: adios_err, ierr

  call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)
  call adios_finalize (myrank, adios_err)
end subroutine adios_cleanup
