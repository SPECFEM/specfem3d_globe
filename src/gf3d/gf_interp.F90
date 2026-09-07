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
!---- GLL interpolation of the stored reciprocal displacement.
!----
!---- The array being interpolated is, in Fortran order,
!----
!----     displacement(3_force, 3_disp, i, j, k, nt_sub)
!----
!---- with the force index N=1, E=2, Z=3 and the displacement index
!---- Cartesian x/y/z, and with xi <-> i, eta <-> j, gamma <-> k. The writer
!---- is the authority for that: src/specfem3D/green_function_io.F90:35-49,
!---- and the force ordering again at :506.
!----
!---- **This index order is the likeliest thing in the library to get
!---- wrong.** GF3DF stored the same quantity as
!---- GF%displacement(ista,iforce,idisp,iglob,it) -- station index leading,
!---- and the three GLL indices collapsed into a single iglob -- so any loop
!---- carried over from it transposes. tests/gf3d/test_gf_interp.f90 drives
!---- these routines with a field carrying three independent factors (a
!---- per-time factor, a per-(force,disp) scale, and a polynomial that is
!---- deliberately *not* symmetric under permuting xi/eta/gamma) precisely so
!---- that each possible transposition fails a different assertion.
!----
!---- The precision seam
!---- ------------------
!---- gf_interp_kernel is double precision and is where all the arithmetic
!---- happens. gf_interp_trace() is the thin wrapper that widens the
!---- real(CUSTOM_REAL) array coming out of HDF5, one time slice at a time.
!----
!---- That seam is a requirement, not a convenience. The unit test's oracle
!---- is P(xi*) for a polynomial P, so its nodal values are irrational
!---- numbers; stored as float32 they carry a relative perturbation of up to
!---- 2^-24, which the 3-D Lebesgue constant (about 4-5 at NGLLX = 5)
!---- amplifies to roughly 3e-7. Without a double-precision entry point the
!---- test could assert nothing tighter than that, and an exact test would
!---- have become a loose one. There is no way around it by choosing nicer
!---- nodal values: making them exactly float32-representable forces P into
!---- the Lagrange basis, which makes the oracle circular.
!----
!---- Widening a slice at a time, rather than the whole array, matters too:
!---- one element-station file is 21 MB in the shipped global example and
!---- must not become 42 MB.
!----
!---- No `use hdf5`, no `use specfem_par`: this is a kernel module.
!----

  module gf_interp

  use gf_par, only: GF_NCOMP

  implicit none

  private

  public :: gf_interp_weights
  public :: gf_interp_weights_deriv
  public :: gf_interp_snapshot
  public :: gf_interp_trace_d
  public :: gf_interp_trace

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_interp_weights_deriv(xi,eta,gamma,hxi,hpxi,heta,hpeta,hgam,hpgam)

! the three per-axis Lagrange factors and their derivatives at a point
!
! Split out from the interpolation itself because they are the same for
! every force component, every displacement component and every time sample:
! computing them once per element rather than 4625 times per element is the
! whole reason this is a separate routine.
!
! The abscissae come from zwgljd() with GAUSSALPHA/GAUSSBETA, matching
! define_derivation_matrices.f90:60-62, and the basis from lagrange_any() --
! the solver's own routines, so a change of quadrature cannot desynchronise
! the extraction from the mesh it is reading.
!
! The derivatives are with respect to the *reference* coordinates; turning
! them into physical-space derivatives is gf_strain's job, and needs the
! inverse Jacobian the locator already produced.

  use constants, only: NGLLX,NGLLY,NGLLZ,GAUSSALPHA,GAUSSBETA

  implicit none

  double precision, intent(in) :: xi,eta,gamma
  double precision, dimension(NGLLX), intent(out) :: hxi,hpxi
  double precision, dimension(NGLLY), intent(out) :: heta,hpeta
  double precision, dimension(NGLLZ), intent(out) :: hgam,hpgam

  ! local parameters
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision, dimension(NGLLZ) :: zigll,wzgll

  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  call lagrange_any(xi,NGLLX,xigll,hxi,hpxi)
  call lagrange_any(eta,NGLLY,yigll,heta,hpeta)
  call lagrange_any(gamma,NGLLZ,zigll,hgam,hpgam)

  end subroutine gf_interp_weights_deriv

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_interp_weights(xi,eta,gamma,hxi,heta,hgam)

! the three per-axis Lagrange factors at a point inside an element
!
! The values only; callers that also need the derivatives (gf_strain) call
! gf_interp_weights_deriv directly.

  use constants, only: NGLLX,NGLLY,NGLLZ

  implicit none

  double precision, intent(in) :: xi,eta,gamma
  double precision, dimension(NGLLX), intent(out) :: hxi
  double precision, dimension(NGLLY), intent(out) :: heta
  double precision, dimension(NGLLZ), intent(out) :: hgam

  ! local parameters
  double precision, dimension(NGLLX) :: hpxi
  double precision, dimension(NGLLY) :: hpeta
  double precision, dimension(NGLLZ) :: hpgam

  call gf_interp_weights_deriv(xi,eta,gamma,hxi,hpxi,heta,hpeta,hgam,hpgam)

  end subroutine gf_interp_weights

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_interp_kernel(u,w,out)

! the tensor-product contraction, on one time slice
!
! Every public routine below funnels through this one, so there is a single
! place where the (i,j,k) <-> (xi,eta,gamma) association is written down.

  use constants, only: NGLLX,NGLLY,NGLLZ

  implicit none

  double precision, dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ), intent(in) :: u
  double precision, dimension(NGLLX,NGLLY,NGLLZ), intent(in) :: w
  double precision, dimension(GF_NCOMP,GF_NCOMP), intent(out) :: out

  ! local parameters
  integer :: i,j,k

  out(:,:) = 0.d0

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        out(:,:) = out(:,:) + w(i,j,k)*u(:,:,i,j,k)
      enddo
    enddo
  enddo

  end subroutine gf_interp_kernel

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_interp_outer(hxi,heta,hgam,w)

! the 125 tensor-product weights, w(i,j,k) = hxi(i)*heta(j)*hgam(k)

  use constants, only: NGLLX,NGLLY,NGLLZ

  implicit none

  double precision, dimension(NGLLX), intent(in) :: hxi
  double precision, dimension(NGLLY), intent(in) :: heta
  double precision, dimension(NGLLZ), intent(in) :: hgam
  double precision, dimension(NGLLX,NGLLY,NGLLZ), intent(out) :: w

  ! local parameters
  integer :: i,j,k

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        w(i,j,k) = hxi(i)*heta(j)*hgam(k)
      enddo
    enddo
  enddo

  end subroutine gf_interp_outer

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_interp_snapshot(u,hxi,heta,hgam,out)

! interpolates one time slice, double precision in and out

  use constants, only: NGLLX,NGLLY,NGLLZ

  implicit none

  double precision, dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ), intent(in) :: u
  double precision, dimension(NGLLX), intent(in) :: hxi
  double precision, dimension(NGLLY), intent(in) :: heta
  double precision, dimension(NGLLZ), intent(in) :: hgam
  double precision, dimension(GF_NCOMP,GF_NCOMP), intent(out) :: out

  ! local parameters
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: w

  call gf_interp_outer(hxi,heta,hgam,w)
  call gf_interp_kernel(u,w,out)

  end subroutine gf_interp_snapshot

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_interp_trace_d(u,hxi,heta,hgam,nt,out)

! interpolates a whole trace, double precision in and out
!
! This is the entry point tests/gf3d/test_gf_interp.f90 asserts at 1e-12; see
! the header for why a double-precision core is a requirement on the library
! rather than a testing convenience.

  use constants, only: NGLLX,NGLLY,NGLLZ

  implicit none

  integer, intent(in) :: nt
  double precision, dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ,nt), intent(in) :: u
  double precision, dimension(NGLLX), intent(in) :: hxi
  double precision, dimension(NGLLY), intent(in) :: heta
  double precision, dimension(NGLLZ), intent(in) :: hgam
  double precision, dimension(GF_NCOMP,GF_NCOMP,nt), intent(out) :: out

  ! local parameters
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: w
  integer :: it

  call gf_interp_outer(hxi,heta,hgam,w)

  do it = 1,nt
    call gf_interp_kernel(u(:,:,:,:,:,it),w,out(:,:,it))
  enddo

  end subroutine gf_interp_trace_d

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_interp_trace(u,hxi,heta,hgam,nt,out)

! interpolates a whole trace read from the database
!
! The only difference from gf_interp_trace_d is the dble() on the way in,
! applied one 225-element time slice at a time so that a 21 MB element-station
! array never has a 42 MB double-precision twin. The arithmetic that follows
! is bit-for-bit the same as the double core's, since widening is exact.

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  integer, intent(in) :: nt
  real(kind=CUSTOM_REAL), dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ,nt), intent(in) :: u
  double precision, dimension(NGLLX), intent(in) :: hxi
  double precision, dimension(NGLLY), intent(in) :: heta
  double precision, dimension(NGLLZ), intent(in) :: hgam
  double precision, dimension(GF_NCOMP,GF_NCOMP,nt), intent(out) :: out

  ! local parameters
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: w
  double precision, dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ) :: ud
  integer :: it

  call gf_interp_outer(hxi,heta,hgam,w)

  do it = 1,nt
    ud(:,:,:,:,:) = dble(u(:,:,:,:,:,it))
    call gf_interp_kernel(ud,w,out(:,:,it))
  enddo

  end subroutine gf_interp_trace

  end module gf_interp
