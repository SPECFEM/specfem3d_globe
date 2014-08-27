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

  double precision function comp_source_time_function(t,hdur)

  use constants

  implicit none

  double precision t,hdur

  double precision, external :: netlib_specfun_erf
  double precision, external :: comp_source_time_function_external

  if ( EXTERNAL_SOURCE_TIME_FUNCTION ) then
    comp_source_time_function = comp_source_time_function_external ()
  else
    ! quasi Heaviside
    comp_source_time_function = 0.5d0*(1.0d0 + netlib_specfun_erf(t/hdur))
  end if

  end function comp_source_time_function

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_rickr(t,f0)

  use constants

  implicit none

  double precision t,f0

  ! Ricker
  comp_source_time_function_rickr = (1.d0 - 2.d0*PI*PI*f0*f0*t*t ) &
                                    * exp( -PI*PI*f0*f0*t*t )

  !!! another source time function they have called 'Ricker' in some old papers,
  !!! e.g., 'Finite-Frequency Kernels Based on Adjoint Methods' by Liu & Tromp, BSSA (2006)
  !!! in order to benchmark those simulations, the following formula is needed.
  ! comp_source_time_function_rickr = -2.d0*PI*PI*f0*f0*f0*t * exp(-PI*PI*f0*f0*t*t)

  end function comp_source_time_function_rickr

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_external ( )

  use specfem_par, only: it, stfArray_external
  implicit none

  ! On the first iteration, go get the ASCII file.
  if ( .not. allocated (stfArray_external) ) then
    call get_external_source_time_function ()
  end if

  comp_source_time_function_external = stfArray_external (it)

  end function comp_source_time_function_external

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_external_source_time_function()

  use specfem_par, only: NSTEP, stfArray_external
  implicit none

  integer :: RetCode, iterator
  character (len=200) :: line

  ! Allocate the source time function array to the number of time steps.
  allocate ( stfArray_external (NSTEP) )
  print *, NSTEP

  ! Read in source time function.
  open (unit=10, file='DATA/stf', status='old', form='formatted')

  read_loop: do iterator=1,NSTEP

    read (10, '(A)', iostat = RetCode) line

    if ( RetCode /= 0 ) then
      print *, "ERROR IN SOURCE TIME FUNCTION."
      stop
    end if

    ! Ignore lines with a hash (comments)
    if ( index (line, "#") /= 0 ) cycle read_loop

    read (line, *) stfArray_external(iterator)

  end do read_loop

  close (10)

  end subroutine get_external_source_time_function
