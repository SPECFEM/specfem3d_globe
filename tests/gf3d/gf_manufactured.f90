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
!---- Shared support for the tests/gf3d/ suite.
!----
!---- Two rules govern what may live here, and both are load-bearing:
!----
!---- 1. **Nothing in src/gf3d/ may be reimplemented here.** These tests
!----    replace a Fortran-versus-Python comparison that was abolished
!----    precisely because testing one implementation against another proves
!----    nothing. Repeating that mistake inside the test directory would be
!----    no better. What belongs here is closed-form arithmetic: polynomial
!----    evaluation, a random-number source, and assertions.
!----
!---- 2. **The random source is a deterministic LCG, not random_number().**
!----    The intrinsic generator is not reproducible across compilers or
!----    library versions, so a CI failure could not be reproduced locally --
!----    which is most of the value of a failing test. The recurrence below
!----    is the classic 2^31 LCG, evaluated in 64-bit so that the product
!----    never overflows (1103515245 * 2^31 is 2.4e18, against a 9.2e18
!----    limit) and therefore never relies on wraparound, which Fortran does
!----    not define.
!----
!---- 3. **Assertions print the measured error even when they pass.** A
!----    pass/fail bit records nothing; a margin lets results.log show drift
!----    before it becomes a failure, and several tolerances in this suite
!----    are documented by what they actually measure -- the float32 noise
!----    floor that Stage 8's finite-difference budget depends on above all.
!----

  module gf_manufactured

  implicit none

  ! LCG state; seeded through gf_lcg_seed()
  integer(kind=8), private :: lcg_state = 12345_8

  integer(kind=8), parameter, private :: LCG_A = 1103515245_8
  integer(kind=8), parameter, private :: LCG_C = 12345_8
  integer(kind=8), parameter, private :: LCG_M = 2147483648_8

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_lcg_seed(seed)

! resets the deterministic random source

  implicit none

  integer, intent(in) :: seed

  lcg_state = modulo(int(seed,kind=8),LCG_M)

  end subroutine gf_lcg_seed

!
!-------------------------------------------------------------------------------------------------
!

  double precision function gf_rand()

! next pseudo-random double in [0,1)

  implicit none

  lcg_state = modulo(LCG_A*lcg_state + LCG_C, LCG_M)
  gf_rand = dble(lcg_state) / dble(LCG_M)

  end function gf_rand

!
!-------------------------------------------------------------------------------------------------
!

  double precision function gf_rand_range(a,b)

! next pseudo-random double in [a,b)

  implicit none

  double precision, intent(in) :: a,b

  gf_rand_range = a + (b-a)*gf_rand()

  end function gf_rand_range

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_report(name,err,tol,nfail)

! records one assertion, printing the measured error whether it passes or not

  implicit none

  character(len=*), intent(in) :: name
  double precision, intent(in) :: err,tol
  integer, intent(inout) :: nfail

  if (err <= tol .and. .not. gf_is_nan(err)) then
    write(*,'(a,a,a,es12.5,a,es12.5,a)') '  ok   ',name,'   error = ',err,'  (tol ',tol,')'
  else
    write(*,'(a,a,a,es12.5,a,es12.5,a)') '  FAIL ',name,'   error = ',err,'  (tol ',tol,')'
    nfail = nfail + 1
  endif

  end subroutine gf_report

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_report_true(name,cond,nfail)

! records one boolean assertion

  implicit none

  character(len=*), intent(in) :: name
  logical, intent(in) :: cond
  integer, intent(inout) :: nfail

  if (cond) then
    write(*,'(a,a)') '  ok   ',name
  else
    write(*,'(a,a)') '  FAIL ',name
    nfail = nfail + 1
  endif

  end subroutine gf_report_true

!
!-------------------------------------------------------------------------------------------------
!

  logical function gf_is_nan(x)

! NaN test that never performs an invalid comparison
!
! The tree is built with -ffpe-trap=invalid in some configurations, where
! `x /= x` on a signalling NaN traps before it can answer. The bit pattern
! is safe to look at: an IEEE double is NaN when its exponent field is all
! ones and its mantissa is non-zero. Same approach as gf_is_finite() in
! src/gf3d/gf_database.F90, and for the same reason.

  implicit none

  double precision, intent(in) :: x

  ! local parameters
  integer(kind=8) :: bits
  integer(kind=8), parameter :: EXP_MASK  = int(z'7FF0000000000000',kind=8)
  integer(kind=8), parameter :: FRAC_MASK = int(z'000FFFFFFFFFFFFF',kind=8)

  bits = transfer(x,bits)

  gf_is_nan = (iand(bits,EXP_MASK) == EXP_MASK) .and. (iand(bits,FRAC_MASK) /= 0_8)

  end function gf_is_nan

!
!-------------------------------------------------------------------------------------------------
!

  double precision function gf_poly3(coef,na,nb,nc,x,y,z)

! evaluates sum_abc coef(a,b,c) * x^a * y^b * z^c
!
! A plain triple sum over monomials: it shares no code with the tensor-
! product Lagrange evaluation it is used to check, which is the property
! that makes it a legitimate oracle.

  implicit none

  integer, intent(in) :: na,nb,nc
  double precision, dimension(0:na,0:nb,0:nc), intent(in) :: coef
  double precision, intent(in) :: x,y,z

  ! local parameters
  integer :: ia,ib,ic
  double precision :: s,xa,yb,zc

  s = 0.d0
  xa = 1.d0
  do ia = 0,na
    yb = 1.d0
    do ib = 0,nb
      zc = 1.d0
      do ic = 0,nc
        s = s + coef(ia,ib,ic)*xa*yb*zc
        zc = zc*z
      enddo
      yb = yb*y
    enddo
    xa = xa*x
  enddo

  gf_poly3 = s

  end function gf_poly3

  end module gf_manufactured
