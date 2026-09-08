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
!---- Strain of the stored reciprocal displacement, at the source point.
!----
!---- For each force component a applied at the station, the database holds a
!---- displacement field u^a_p(x). Its strain at the source is
!----
!----   eps^a_pq = 1/2 ( du^a_p/dx_q + du^a_q/dx_p )
!----
!---- and by reciprocity the seismogram at the station in component a, from a
!---- moment tensor M at the source, is the contraction
!----
!----   seis(a,t) = SUM_pq M_pq eps^a_pq(x_s,t)
!----
!---- which gf_moment.F90 performs. This module produces the strain and
!---- stops there, deliberately: Stage 6's moment-tensor partial derivatives
!---- are the *same* strain re-contracted with six unit tensors, so it must
!---- be possible to re-contract without recomputing.
!----
!---- The chain rule, and why the weights are precomputed
!---- --------------------------------------------------
!---- The basis derivatives that come out of lagrange_any are with respect to
!---- the reference coordinates, so
!----
!----   du_p/dx_q = SUM_ijk u_p(i,j,k) *
!----       [ h'(i)h(j)h(k) dxi/dx_q + h(i)h'(j)h(k) deta/dx_q
!----                                + h(i)h(j)h'(k) dgamma/dx_q ]
!----
!---- The bracket depends only on the point, not on the field, the force
!---- component or the time sample. So it is built once as dw(i,j,k,q) --
!---- 375 numbers -- and every one of the 4625 time slices then costs a plain
!---- contraction. That is the fused structure GF3DF used, kept for the same
!---- reason: it never materialises a per-GLL-point intermediate.
!----
!---- Voigt packing is [xx, yy, zz, xy, xz, yz]. The factor-2 weights on the
!---- off-diagonals belong to the *contraction*, not to the strain, and live
!---- in gf_moment.F90 -- so what this module returns is the strain itself,
!---- and eps(4) really is eps_xy and not 2*eps_xy.
!----
!---- The precision seam, as in gf_interp
!---- -----------------------------------
!---- gf_strain_kernel is double precision; gf_strain_trace widens the
!---- real(CUSTOM_REAL) array from HDF5 one time slice at a time. Without a
!---- double entry point tests/gf3d/test_gf_strain.f90 could assert nothing
!---- tighter than the float32 floor, and the geometric chain it checks --
!---- shape functions, Jacobian, inversion, chain rule, symmetrisation,
!---- Voigt packing -- is exactly the part that must be exact.
!----
!---- No `use hdf5`, no `use specfem_par`: this is a kernel module.
!----

  module gf_strain

  use gf_par, only: GF_NCOMP

  implicit none

  private

  ! Voigt slots, named so that a reader never has to count
  integer, parameter, public :: GF_VOIGT = 6
  integer, parameter, public :: GF_XX = 1, GF_YY = 2, GF_ZZ = 3
  integer, parameter, public :: GF_XY = 4, GF_XZ = 5, GF_YZ = 6

  public :: gf_strain_dweights
  public :: gf_strain_ddweights
  public :: gf_strain_snapshot
  public :: gf_strain_trace_d
  public :: gf_strain_trace

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_strain_dweights(hxi,hpxi,heta,hpeta,hgam,hpgam,jinv,dw)

! physical-space derivatives of the 125 basis functions at a point
!
! dw(i,j,k,q) = d/dx_q of the (i,j,k) basis function, i.e. the reference
! derivatives contracted with the inverse Jacobian that gf_shape3D_map
! returned at this same point. `jinv` has rows xi/eta/gamma and columns
! x/y/z, so jinv(1,q) is dxi/dx_q.
!
! This is where a transposed inverse Jacobian would enter, which is why the
! test's oracle differentiates in physical space and never forms a jinv of
! its own.

  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM

  implicit none

  double precision, dimension(NGLLX), intent(in) :: hxi,hpxi
  double precision, dimension(NGLLY), intent(in) :: heta,hpeta
  double precision, dimension(NGLLZ), intent(in) :: hgam,hpgam
  double precision, dimension(NDIM,NDIM), intent(in) :: jinv
  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM), intent(out) :: dw

  ! local parameters
  integer :: i,j,k,q
  double precision :: dxi,deta,dgam

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        dxi  = hpxi(i)*heta(j)*hgam(k)
        deta = hxi(i)*hpeta(j)*hgam(k)
        dgam = hxi(i)*heta(j)*hpgam(k)
        do q = 1,NDIM
          dw(i,j,k,q) = dxi*jinv(1,q) + deta*jinv(2,q) + dgam*jinv(3,q)
        enddo
      enddo
    enddo
  enddo

  end subroutine gf_strain_dweights

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_strain_ddweights(hxi,hpxi,hppxi,heta,hpeta,hppeta,hgam,hpgam,hppgam,jinv,djinv,ddw)

! the reference-coordinate derivatives of the weight table above
!
! ddw(i,j,k,q,a) = d dw(i,j,k,q) / d xi_a. With dw = SUM_b g_b jinv(b,q),
! where g_b is the reference derivative of the (i,j,k) basis function
! along xi_b,
!
!   ddw(..,q,a) = SUM_b (d g_b / d xi_a) jinv(b,q) + SUM_b g_b djinv(b,q,a)
!
! `djinv(:,:,a)` is d jinv / d xi_a from gf_shape3D_map_2nd, the exact
! derivative of the 27-anchor geometry. This is the whole of Stage 8's
! strain gradient: gf_strain_kernel applied with ddw(:,:,:,:,a) in place
! of dw returns d eps / d xi_a (the symmetrisation is linear), and the
! physical gradient is d eps/dx_m = SUM_a jinv(a,m) d eps/d xi_a. Nothing
! else in this module changes, which is what makes the derivative of the
! strain the derivative of *this* interpolant and not of another one.

  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM

  implicit none

  double precision, dimension(NGLLX), intent(in) :: hxi,hpxi,hppxi
  double precision, dimension(NGLLY), intent(in) :: heta,hpeta,hppeta
  double precision, dimension(NGLLZ), intent(in) :: hgam,hpgam,hppgam
  double precision, dimension(NDIM,NDIM), intent(in) :: jinv
  double precision, dimension(NDIM,NDIM,NDIM), intent(in) :: djinv
  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM,NDIM), intent(out) :: ddw

  ! local parameters
  integer :: i,j,k,q,a,b
  double precision, dimension(NDIM) :: g
  double precision, dimension(NDIM,NDIM) :: dg   ! dg(b,a) = d g_b / d xi_a

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        g(1) = hpxi(i)*heta(j)*hgam(k)
        g(2) = hxi(i)*hpeta(j)*hgam(k)
        g(3) = hxi(i)*heta(j)*hpgam(k)

        dg(1,1) = hppxi(i)*heta(j)*hgam(k)
        dg(2,2) = hxi(i)*hppeta(j)*hgam(k)
        dg(3,3) = hxi(i)*heta(j)*hppgam(k)
        dg(1,2) = hpxi(i)*hpeta(j)*hgam(k)
        dg(1,3) = hpxi(i)*heta(j)*hpgam(k)
        dg(2,3) = hxi(i)*hpeta(j)*hpgam(k)
        dg(2,1) = dg(1,2)
        dg(3,1) = dg(1,3)
        dg(3,2) = dg(2,3)

        do a = 1,NDIM
          do q = 1,NDIM
            ddw(i,j,k,q,a) = 0.d0
            do b = 1,NDIM
              ddw(i,j,k,q,a) = ddw(i,j,k,q,a) + dg(b,a)*jinv(b,q) + g(b)*djinv(b,q,a)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  end subroutine gf_strain_ddweights

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_strain_kernel(u,dw,eps)

! the strain of one time slice, for all three force components
!
! Every public routine funnels through this one, so the gradient, the
! symmetrisation and the Voigt packing are each written down exactly once.

  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM

  implicit none

  double precision, dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ), intent(in) :: u
  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM), intent(in) :: dw
  double precision, dimension(GF_VOIGT,GF_NCOMP), intent(out) :: eps

  ! local parameters
  ! grad(p,q) = du_p/dx_q for one force component
  double precision, dimension(NDIM,NDIM) :: grad
  integer :: i,j,k,a,p,q

  do a = 1,GF_NCOMP

    grad(:,:) = 0.d0

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          do q = 1,NDIM
            do p = 1,NDIM
              ! u(a,p,...) is the p-th Cartesian component of the field
              ! driven by the a-th force component at the station
              grad(p,q) = grad(p,q) + u(a,p,i,j,k)*dw(i,j,k,q)
            enddo
          enddo
        enddo
      enddo
    enddo

    ! symmetrise into Voigt [xx, yy, zz, xy, xz, yz]
    !
    ! The diagonal terms need no factor: eps_xx = du_x/dx exactly. The
    ! off-diagonals carry the 1/2 of the symmetric gradient, so eps(4) is
    ! eps_xy and not 2*eps_xy. The factor 2 that appears in the contraction
    ! SUM_pq M_pq eps_pq belongs to gf_moment, not here.
    eps(GF_XX,a) = grad(1,1)
    eps(GF_YY,a) = grad(2,2)
    eps(GF_ZZ,a) = grad(3,3)
    eps(GF_XY,a) = 0.5d0*(grad(1,2) + grad(2,1))
    eps(GF_XZ,a) = 0.5d0*(grad(1,3) + grad(3,1))
    eps(GF_YZ,a) = 0.5d0*(grad(2,3) + grad(3,2))

  enddo

  end subroutine gf_strain_kernel

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_strain_snapshot(u,dw,eps)

! strain of one time slice, double precision in and out

  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM

  implicit none

  double precision, dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ), intent(in) :: u
  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM), intent(in) :: dw
  double precision, dimension(GF_VOIGT,GF_NCOMP), intent(out) :: eps

  call gf_strain_kernel(u,dw,eps)

  end subroutine gf_strain_snapshot

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_strain_trace_d(u,dw,nt,eps)

! strain of a whole trace, double precision in and out
!
! The entry point tests/gf3d/test_gf_strain.f90 asserts at 1e-12.

  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM

  implicit none

  integer, intent(in) :: nt
  double precision, dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ,nt), intent(in) :: u
  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM), intent(in) :: dw
  double precision, dimension(GF_VOIGT,GF_NCOMP,nt), intent(out) :: eps

  ! local parameters
  integer :: it

  do it = 1,nt
    call gf_strain_kernel(u(:,:,:,:,:,it),dw,eps(:,:,it))
  enddo

  end subroutine gf_strain_trace_d

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_strain_trace(u,dw,nt,eps)

! strain of a whole trace read from the database
!
! Widening is applied one 225-element time slice at a time, so the 21 MB
! element-station array never has a 42 MB double-precision twin, and the
! arithmetic that follows is bit-for-bit the double core's.

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM

  implicit none

  integer, intent(in) :: nt
  real(kind=CUSTOM_REAL), dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ,nt), intent(in) :: u
  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM), intent(in) :: dw
  double precision, dimension(GF_VOIGT,GF_NCOMP,nt), intent(out) :: eps

  ! local parameters
  double precision, dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ) :: ud
  integer :: it

  do it = 1,nt
    ud(:,:,:,:,:) = dble(u(:,:,:,:,:,it))
    call gf_strain_kernel(ud,dw,eps(:,:,it))
  enddo

  end subroutine gf_strain_trace

  end module gf_strain
