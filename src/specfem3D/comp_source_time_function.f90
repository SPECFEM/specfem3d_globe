!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
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

  double precision function comp_source_time_function(t,hdur)

  implicit none

  include "constants.h"

  double precision t,hdur

  double precision, external :: netlib_specfun_erf

! quasi Heaviside
  comp_source_time_function = 0.5d0*(1.0d0 + netlib_specfun_erf(t/hdur))

  end function comp_source_time_function

!
!-------------------------------------------------------------------------------------------------
!

  double precision function comp_source_time_function_rickr(t,f0)

  implicit none

  include "constants.h"

  double precision t,f0

  ! ricker
  comp_source_time_function_rickr = (1.d0 - 2.d0*PI*PI*f0*f0*t*t ) &
                                    * exp( -PI*PI*f0*f0*t*t )

  !!! another source time function they have called 'ricker' in some old papers,
  !!! e.g., 'Finite-Frequency Kernels Based on Adjoint Methods' by Liu & Tromp, BSSA (2006)
  !!! in order to benchmark those simulations, the following formula is needed.
  ! comp_source_time_function_rickr = -2.d0*PI*PI*f0*f0*f0*t * exp(-PI*PI*f0*f0*t*t)

  end function comp_source_time_function_rickr
