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

subroutine read_parameters_tomo()

! reads in parameters needed (only step length for now...)

  use tomography_par

  implicit none
  character(len=150) :: s_step_fac

  ! subjective step length to multiply to the gradient
  ! e.g. step_fac = 0.03

  call getarg(1,s_step_fac)

  if (trim(s_step_fac) == '') then
    call exit_MPI(myrank,'Usage: add_model step_factor')
  endif

  ! read in parameter information
  read(s_step_fac,*) step_fac

  ! safety check
  if (abs(step_fac) < 1.e-15) then
    print*,'Error: step factor ',step_fac,' is too small and will lead to no update...'
    call exit_MPI(myrank,'Error step factor too small')
  endif

end subroutine read_parameters_tomo

