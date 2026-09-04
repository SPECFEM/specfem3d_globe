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
!---- Source time function operators.
!----
!---- Stage 4 needs only the time integration; the Gaussian half-duration
!---- correction, the output time axis and its left zero-padding are Stage
!---- 5's and will be added here beside it.
!----
!---- Why an integration is needed at all
!---- -----------------------------------
!---- The reciprocal runs that filled the database used a Gaussian source
!---- time function of half duration hdur_gf, so what is stored is
!----
!----   G ⊛ g_(hdur_gf)
!----
!---- for the Green function G. A CMT source is a (quasi-)Heaviside, and
!---- specfem's Heaviside is the integral of exactly that Gaussian
!---- (comp_source_time_function.f90:202 is a unit-area Gaussian), so
!----
!----   G ⊛ H_(hdur_gf) = INTEGRAL of ( G ⊛ g_(hdur_gf) )
!----
!---- One cumulative integration therefore converts the stored Gaussian
!---- response into the Heaviside response *at the database's own half
!---- duration*. Moving from hdur_gf to the CMTSOLUTION's half duration is a
!---- separate Gaussian convolution, and belongs to Stage 5.
!----
!---- Trapezoid, not a left-endpoint cumulative sum
!---- ---------------------------------------------
!---- A left Riemann sum lags the trapezoid by exactly (h/2)(f_n - f_0),
!---- i.e. by half a sample. On the 0.4 s subsampled grid of the shipped
!---- global example that is a 0.2 s phase error -- small enough to look
!---- like a mediocre fit rather than a bug, which is what makes it worth
!---- being explicit about. The Python's own comment at reconstruct_cmt
!---- makes the same point.
!----
!---- No `use hdf5`, no `use specfem_par`: this is a kernel module.
!----

  module gf_stf

  implicit none

  private

  public :: gf_cumtrapz

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_cumtrapz(f,n,dt,g)

! cumulative trapezoidal integration on a uniform grid
!
!   g(1) = 0
!   g(i) = g(i-1) + (dt/2)(f(i-1) + f(i))
!
! `dt` is the spacing of the axis being integrated along, which for a
! database trace is dt_sub = dt*subsample_step and *not* the solver step --
! 0.4 s rather than 0.1 s in the shipped global example. Passing the solver
! step here is wrong by a factor of subsample_step and produces a trace of
! entirely plausible shape.
!
! g may alias f; the running value is held in a scalar so that an in-place
! call is well defined.

  implicit none

  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: f
  double precision, intent(in) :: dt
  double precision, dimension(n), intent(out) :: g

  ! local parameters
  integer :: i
  double precision :: acc,fprev,fcur

  if (n < 1) return

  acc = 0.d0
  g(1) = 0.d0
  if (n == 1) return

  fprev = f(1)
  do i = 2,n
    fcur = f(i)
    acc = acc + 0.5d0*dt*(fprev + fcur)
    g(i) = acc
    fprev = fcur
  enddo

  end subroutine gf_cumtrapz

  end module gf_stf
