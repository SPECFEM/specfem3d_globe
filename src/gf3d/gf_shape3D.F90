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
!---- Tri-quadratic (27-node) element geometry: shape functions, the
!---- Jacobian matrix and its inverse.
!----
!---- This is a *fork* of src/shared/recompute_jacobian.f90, not a caller
!---- of it, for one reason: that routine ends a degenerate element with
!---- `stop '3D Jacobian undefined'` (recompute_jacobian.f90:261), and a
!---- `stop` inside a shared object that Python has loaded (Stage 9) kills
!---- the interpreter with no traceback. The whole library returns `ierr`
!---- instead; this file is where that promise would otherwise break.
!----
!---- Forking is cheap here because it buys a free oracle: the original is
!---- still in the tree and still linked, so tests/gf3d/test_gf_shape3D.f90
!---- compares the two over random anchors and random (xi,eta,gamma). The
!---- arithmetic below is transcribed in the original's exact expression
!---- and accumulation order rather than tidied up, so that under a
!---- value-safe floating-point model the two agree bit for bit. The test
!---- asserts a derived tolerance rather than exact equality, because the
!---- compiler is not bound by that order: ifort/ifx at -O3 -xHost contract
!---- and vectorise the two compilation units differently and disagree in
!---- the last bits (see the test's header). The bit-for-bit count is still
!---- reported there.
!----
!---- Two deliberate additions over the original:
!----
!----   * `jacobian` is returned. recompute_jacobian keeps it as a local,
!----     but Stage 4's strain needs it and it is already computed.
!----   * the nine inverse-Jacobian terms are packed into jinv(NDIM,NDIM)
!----     with rows xi/eta/gamma and columns x/y/z, so a caller cannot
!----     transpose them by mis-ordering nine positional arguments. The
!----     packing itself is covered by the bit-for-bit test.
!----
!---- The second derivatives Stage 8 needs are gf_shape3D_functions_2nd
!---- and gf_shape3D_map_2nd, siblings of the two above the way
!---- lagrange_any_2nd sits beside lagrange_any in
!---- src/shared/lagrange_poly.f90 -- not optional arguments. The
!---- tri-quadratic shape functions have constant second derivatives per
!---- axis (l1'' = 1, l2'' = -2, l3'' = 1), so the second derivatives of
!---- the geometry map are exact and cost nothing in conditioning; the
!---- derivative of the inverse Jacobian follows from
!---- d(J^-1)/dxi_a = -J^-1 (dJ/dxi_a) J^-1.
!----
!---- No `use hdf5`, no `use specfem_par`: this is a kernel module and
!---- tests/gf3d/0.configure.default_make.sh enforces that with `nm`.
!----

  module gf_shape3D

  use gf_par, only: gf_set_error,GF_OK,GF_ERR_MISMATCH,GF_ERR_GEOMETRY

  implicit none

  private

  public :: gf_shape3D_functions
  public :: gf_shape3D_functions_2nd
  public :: gf_shape3D_map
  public :: gf_shape3D_map_2nd

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_shape3D_functions(xi,eta,gamma,shape3D,dershape3D)

! tri-quadratic shape functions and their first derivatives at (xi,eta,gamma)
!
! Transcribed from recompute_jacobian.f90:67-218. The node ordering is
! specfem's 27-anchor convention -- 1-8 corners, 9-20 edge midpoints,
! 21-26 face centres, 27 the centre -- matching hex_nodes() in
! src/shared/hex_nodes.f90, so anchors gathered with hex_nodes_anchor_ijk()
! line up with these functions without any reordering.

  use constants, only: NGNOD,NDIM,HALF,ONE,TWO

  implicit none

  double precision, intent(in) :: xi,eta,gamma
  double precision, dimension(NGNOD), intent(out) :: shape3D
  double precision, dimension(NDIM,NGNOD), intent(out) :: dershape3D

  ! local parameters
  double precision :: l1xi,l2xi,l3xi
  double precision :: l1eta,l2eta,l3eta
  double precision :: l1gamma,l2gamma,l3gamma
  double precision :: l1pxi,l2pxi,l3pxi
  double precision :: l1peta,l2peta,l3peta
  double precision :: l1pgamma,l2pgamma,l3pgamma

  l1xi = HALF*xi*(xi-ONE)
  l2xi = ONE-xi**2
  l3xi = HALF*xi*(xi+ONE)

  l1pxi = xi-HALF
  l2pxi=-TWO*xi
  l3pxi = xi+HALF

  l1eta = HALF*eta*(eta-ONE)
  l2eta = ONE-eta**2
  l3eta = HALF*eta*(eta+ONE)

  l1peta = eta-HALF
  l2peta=-TWO*eta
  l3peta = eta+HALF

  l1gamma = HALF*gamma*(gamma-ONE)
  l2gamma = ONE-gamma**2
  l3gamma = HALF*gamma*(gamma+ONE)

  l1pgamma = gamma-HALF
  l2pgamma=-TWO*gamma
  l3pgamma = gamma+HALF

  ! corner nodes

  shape3D(1)=l1xi*l1eta*l1gamma
  shape3D(2)=l3xi*l1eta*l1gamma
  shape3D(3)=l3xi*l3eta*l1gamma
  shape3D(4)=l1xi*l3eta*l1gamma
  shape3D(5)=l1xi*l1eta*l3gamma
  shape3D(6)=l3xi*l1eta*l3gamma
  shape3D(7)=l3xi*l3eta*l3gamma
  shape3D(8)=l1xi*l3eta*l3gamma

  dershape3D(1,1)=l1pxi*l1eta*l1gamma
  dershape3D(1,2)=l3pxi*l1eta*l1gamma
  dershape3D(1,3)=l3pxi*l3eta*l1gamma
  dershape3D(1,4)=l1pxi*l3eta*l1gamma
  dershape3D(1,5)=l1pxi*l1eta*l3gamma
  dershape3D(1,6)=l3pxi*l1eta*l3gamma
  dershape3D(1,7)=l3pxi*l3eta*l3gamma
  dershape3D(1,8)=l1pxi*l3eta*l3gamma

  dershape3D(2,1)=l1xi*l1peta*l1gamma
  dershape3D(2,2)=l3xi*l1peta*l1gamma
  dershape3D(2,3)=l3xi*l3peta*l1gamma
  dershape3D(2,4)=l1xi*l3peta*l1gamma
  dershape3D(2,5)=l1xi*l1peta*l3gamma
  dershape3D(2,6)=l3xi*l1peta*l3gamma
  dershape3D(2,7)=l3xi*l3peta*l3gamma
  dershape3D(2,8)=l1xi*l3peta*l3gamma

  dershape3D(3,1)=l1xi*l1eta*l1pgamma
  dershape3D(3,2)=l3xi*l1eta*l1pgamma
  dershape3D(3,3)=l3xi*l3eta*l1pgamma
  dershape3D(3,4)=l1xi*l3eta*l1pgamma
  dershape3D(3,5)=l1xi*l1eta*l3pgamma
  dershape3D(3,6)=l3xi*l1eta*l3pgamma
  dershape3D(3,7)=l3xi*l3eta*l3pgamma
  dershape3D(3,8)=l1xi*l3eta*l3pgamma

  ! midside nodes

  shape3D(9)=l2xi*l1eta*l1gamma
  shape3D(10)=l3xi*l2eta*l1gamma
  shape3D(11)=l2xi*l3eta*l1gamma
  shape3D(12)=l1xi*l2eta*l1gamma
  shape3D(13)=l1xi*l1eta*l2gamma
  shape3D(14)=l3xi*l1eta*l2gamma
  shape3D(15)=l3xi*l3eta*l2gamma
  shape3D(16)=l1xi*l3eta*l2gamma
  shape3D(17)=l2xi*l1eta*l3gamma
  shape3D(18)=l3xi*l2eta*l3gamma
  shape3D(19)=l2xi*l3eta*l3gamma
  shape3D(20)=l1xi*l2eta*l3gamma

  dershape3D(1,9)=l2pxi*l1eta*l1gamma
  dershape3D(1,10)=l3pxi*l2eta*l1gamma
  dershape3D(1,11)=l2pxi*l3eta*l1gamma
  dershape3D(1,12)=l1pxi*l2eta*l1gamma
  dershape3D(1,13)=l1pxi*l1eta*l2gamma
  dershape3D(1,14)=l3pxi*l1eta*l2gamma
  dershape3D(1,15)=l3pxi*l3eta*l2gamma
  dershape3D(1,16)=l1pxi*l3eta*l2gamma
  dershape3D(1,17)=l2pxi*l1eta*l3gamma
  dershape3D(1,18)=l3pxi*l2eta*l3gamma
  dershape3D(1,19)=l2pxi*l3eta*l3gamma
  dershape3D(1,20)=l1pxi*l2eta*l3gamma

  dershape3D(2,9)=l2xi*l1peta*l1gamma
  dershape3D(2,10)=l3xi*l2peta*l1gamma
  dershape3D(2,11)=l2xi*l3peta*l1gamma
  dershape3D(2,12)=l1xi*l2peta*l1gamma
  dershape3D(2,13)=l1xi*l1peta*l2gamma
  dershape3D(2,14)=l3xi*l1peta*l2gamma
  dershape3D(2,15)=l3xi*l3peta*l2gamma
  dershape3D(2,16)=l1xi*l3peta*l2gamma
  dershape3D(2,17)=l2xi*l1peta*l3gamma
  dershape3D(2,18)=l3xi*l2peta*l3gamma
  dershape3D(2,19)=l2xi*l3peta*l3gamma
  dershape3D(2,20)=l1xi*l2peta*l3gamma

  dershape3D(3,9)=l2xi*l1eta*l1pgamma
  dershape3D(3,10)=l3xi*l2eta*l1pgamma
  dershape3D(3,11)=l2xi*l3eta*l1pgamma
  dershape3D(3,12)=l1xi*l2eta*l1pgamma
  dershape3D(3,13)=l1xi*l1eta*l2pgamma
  dershape3D(3,14)=l3xi*l1eta*l2pgamma
  dershape3D(3,15)=l3xi*l3eta*l2pgamma
  dershape3D(3,16)=l1xi*l3eta*l2pgamma
  dershape3D(3,17)=l2xi*l1eta*l3pgamma
  dershape3D(3,18)=l3xi*l2eta*l3pgamma
  dershape3D(3,19)=l2xi*l3eta*l3pgamma
  dershape3D(3,20)=l1xi*l2eta*l3pgamma

  ! side center nodes

  shape3D(21)=l2xi*l2eta*l1gamma
  shape3D(22)=l2xi*l1eta*l2gamma
  shape3D(23)=l3xi*l2eta*l2gamma
  shape3D(24)=l2xi*l3eta*l2gamma
  shape3D(25)=l1xi*l2eta*l2gamma
  shape3D(26)=l2xi*l2eta*l3gamma

  dershape3D(1,21)=l2pxi*l2eta*l1gamma
  dershape3D(1,22)=l2pxi*l1eta*l2gamma
  dershape3D(1,23)=l3pxi*l2eta*l2gamma
  dershape3D(1,24)=l2pxi*l3eta*l2gamma
  dershape3D(1,25)=l1pxi*l2eta*l2gamma
  dershape3D(1,26)=l2pxi*l2eta*l3gamma

  dershape3D(2,21)=l2xi*l2peta*l1gamma
  dershape3D(2,22)=l2xi*l1peta*l2gamma
  dershape3D(2,23)=l3xi*l2peta*l2gamma
  dershape3D(2,24)=l2xi*l3peta*l2gamma
  dershape3D(2,25)=l1xi*l2peta*l2gamma
  dershape3D(2,26)=l2xi*l2peta*l3gamma

  dershape3D(3,21)=l2xi*l2eta*l1pgamma
  dershape3D(3,22)=l2xi*l1eta*l2pgamma
  dershape3D(3,23)=l3xi*l2eta*l2pgamma
  dershape3D(3,24)=l2xi*l3eta*l2pgamma
  dershape3D(3,25)=l1xi*l2eta*l2pgamma
  dershape3D(3,26)=l2xi*l2eta*l3pgamma

  ! center node
  shape3D(27) = l2xi*l2eta*l2gamma

  dershape3D(1,27) = l2pxi*l2eta*l2gamma
  dershape3D(2,27) = l2xi*l2peta*l2gamma
  dershape3D(3,27) = l2xi*l2eta*l2pgamma

  end subroutine gf_shape3D_functions

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_shape3D_map(xelm,yelm,zelm,xi,eta,gamma,xyz,jinv,jacobian,ierr)

! maps (xi,eta,gamma) to a physical position and returns the inverse Jacobian
!
! The forward map is x = sum_ia shape3D(ia) * xelm(ia) over the 27 anchors;
! `jinv` is d(xi,eta,gamma)/d(x,y,z), rows xi/eta/gamma, columns x/y/z --
! so jinv(1,2) is what recompute_jacobian calls `xiy`.
!
! Errors, where the original would have stopped:
!   GF_ERR_MISMATCH  the build has NGNOD /= 27
!   GF_ERR_GEOMETRY  jacobian <= 0, i.e. a degenerate or inverted element
! On either, `xyz` still holds the forward map (it is computed before the
! Jacobian is tested and is meaningful even for a degenerate element) while
! `jinv` is zeroed, so a caller that ignores ierr gets zeros rather than
! plausible garbage.

  use constants, only: NGNOD,NDIM,ZERO

  implicit none

  double precision, dimension(NGNOD), intent(in) :: xelm,yelm,zelm
  double precision, intent(in) :: xi,eta,gamma
  double precision, dimension(NDIM), intent(out) :: xyz
  double precision, dimension(NDIM,NDIM), intent(out) :: jinv
  double precision, intent(out) :: jacobian
  integer, intent(out) :: ierr

  ! local parameters
  double precision :: shape3D(NGNOD)
  double precision :: dershape3D(NDIM,NGNOD)
  double precision :: x,y,z
  double precision :: xxi,yxi,zxi
  double precision :: xeta,yeta,zeta
  double precision :: xgamma,ygamma,zgamma
  integer :: ia

  xyz(:) = ZERO
  jinv(:,:) = ZERO
  jacobian = ZERO

  ! the shape functions above are the 27-node forms and nothing else;
  ! recompute_jacobian.f90:65 stops here instead
  if (NGNOD /= 27) then
    call gf_set_error(ierr,GF_ERR_MISMATCH,'elements must have 27 control nodes (NGNOD /= 27)')
    return
  endif

  call gf_shape3D_functions(xi,eta,gamma,shape3D,dershape3D)

  ! compute coordinates and Jacobian matrix
  x = ZERO
  y = ZERO
  z = ZERO

  xxi = ZERO
  xeta = ZERO
  xgamma = ZERO
  yxi = ZERO
  yeta = ZERO
  ygamma = ZERO
  zxi = ZERO
  zeta = ZERO
  zgamma = ZERO

  do ia = 1,NGNOD
    x = x+shape3D(ia)*xelm(ia)
    y = y+shape3D(ia)*yelm(ia)
    z = z+shape3D(ia)*zelm(ia)

    xxi = xxi+dershape3D(1,ia)*xelm(ia)
    xeta = xeta+dershape3D(2,ia)*xelm(ia)
    xgamma = xgamma+dershape3D(3,ia)*xelm(ia)
    yxi = yxi+dershape3D(1,ia)*yelm(ia)
    yeta = yeta+dershape3D(2,ia)*yelm(ia)
    ygamma = ygamma+dershape3D(3,ia)*yelm(ia)
    zxi = zxi+dershape3D(1,ia)*zelm(ia)
    zeta = zeta+dershape3D(2,ia)*zelm(ia)
    zgamma = zgamma+dershape3D(3,ia)*zelm(ia)
  enddo

  xyz(1) = x
  xyz(2) = y
  xyz(3) = z

  jacobian = xxi*(yeta*zgamma-ygamma*zeta) - xeta*(yxi*zgamma-ygamma*zxi) + &
             xgamma*(yxi*zeta-yeta*zxi)

  if (jacobian <= ZERO) then
    call gf_set_error(ierr,GF_ERR_GEOMETRY,'3D Jacobian undefined: element is degenerate or inverted')
    return
  endif

  ! invert the relation (Fletcher p. 50 vol. 2)
  jinv(1,1) = (yeta*zgamma-ygamma*zeta)/jacobian     ! xix
  jinv(1,2) = (xgamma*zeta-xeta*zgamma)/jacobian     ! xiy
  jinv(1,3) = (xeta*ygamma-xgamma*yeta)/jacobian     ! xiz
  jinv(2,1) = (ygamma*zxi-yxi*zgamma)/jacobian       ! etax
  jinv(2,2) = (xxi*zgamma-xgamma*zxi)/jacobian       ! etay
  jinv(2,3) = (xgamma*yxi-xxi*ygamma)/jacobian       ! etaz
  jinv(3,1) = (yxi*zeta-yeta*zxi)/jacobian           ! gammax
  jinv(3,2) = (xeta*zxi-xxi*zeta)/jacobian           ! gammay
  jinv(3,3) = (xxi*yeta-xeta*yxi)/jacobian           ! gammaz

  ierr = GF_OK

  end subroutine gf_shape3D_map

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_shape3D_functions_2nd(xi,eta,gamma,d2shape3D)

! second derivatives of the tri-quadratic shape functions at (xi,eta,gamma)
!
! d2shape3D(a,b,ia) = d^2 h_ia / (dxi_a dxi_b), symmetric in (a,b), with
! the same 27-node ordering as gf_shape3D_functions. Each shape function
! is a product l_i(xi) l_j(eta) l_k(gamma) of one-dimensional quadratics,
! so the pure second derivatives use l'' = (1, -2, 1) for (l1, l2, l3) and
! the mixed ones the products of two first derivatives. The one-dimensional
! factors are evaluated as in gf_shape3D_functions, and the per-node
! (i,j,k) assignment is read off that routine's table rather than
! re-derived, so a node reordering there breaks tests/gf3d/test_gf_partials
! here instead of silently decoupling the two.

  use constants, only: NGNOD,NDIM,HALF,ONE,TWO

  implicit none

  double precision, intent(in) :: xi,eta,gamma
  double precision, dimension(NDIM,NDIM,NGNOD), intent(out) :: d2shape3D

  ! local parameters
  ! the one-dimensional factors: value, first and second derivative, for
  ! l1 (node at -1), l2 (node at 0), l3 (node at +1), per axis
  double precision, dimension(3) :: lx,lpx,lppx,ly,lpy,lppy,lz,lpz,lppz
  ! which one-dimensional factor (1, 2 or 3) each of the 27 nodes uses per
  ! axis, transcribed from the shape3D table of gf_shape3D_functions
  integer, dimension(NGNOD), parameter :: ix = (/ 1,3,3,1,1,3,3,1, 2,3,2,1,1,3,3,1,2,3,2,1, 2,2,3,2,1,2, 2 /)
  integer, dimension(NGNOD), parameter :: iy = (/ 1,1,3,3,1,1,3,3, 1,2,3,2,1,1,3,3,1,2,3,2, 2,1,2,3,2,2, 2 /)
  integer, dimension(NGNOD), parameter :: iz = (/ 1,1,1,1,3,3,3,3, 1,1,1,1,2,2,2,2,3,3,3,3, 1,2,2,2,2,3, 2 /)
  integer :: ia

  lx(1) = HALF*xi*(xi-ONE) ; lx(2) = ONE-xi**2 ; lx(3) = HALF*xi*(xi+ONE)
  lpx(1) = xi-HALF ; lpx(2) = -TWO*xi ; lpx(3) = xi+HALF
  lppx(1) = ONE ; lppx(2) = -TWO ; lppx(3) = ONE

  ly(1) = HALF*eta*(eta-ONE) ; ly(2) = ONE-eta**2 ; ly(3) = HALF*eta*(eta+ONE)
  lpy(1) = eta-HALF ; lpy(2) = -TWO*eta ; lpy(3) = eta+HALF
  lppy(1) = ONE ; lppy(2) = -TWO ; lppy(3) = ONE

  lz(1) = HALF*gamma*(gamma-ONE) ; lz(2) = ONE-gamma**2 ; lz(3) = HALF*gamma*(gamma+ONE)
  lpz(1) = gamma-HALF ; lpz(2) = -TWO*gamma ; lpz(3) = gamma+HALF
  lppz(1) = ONE ; lppz(2) = -TWO ; lppz(3) = ONE

  do ia = 1,NGNOD
    d2shape3D(1,1,ia) = lppx(ix(ia)) * ly(iy(ia))   * lz(iz(ia))
    d2shape3D(2,2,ia) = lx(ix(ia))   * lppy(iy(ia)) * lz(iz(ia))
    d2shape3D(3,3,ia) = lx(ix(ia))   * ly(iy(ia))   * lppz(iz(ia))
    d2shape3D(1,2,ia) = lpx(ix(ia))  * lpy(iy(ia))  * lz(iz(ia))
    d2shape3D(1,3,ia) = lpx(ix(ia))  * ly(iy(ia))   * lpz(iz(ia))
    d2shape3D(2,3,ia) = lx(ix(ia))   * lpy(iy(ia))  * lpz(iz(ia))
    d2shape3D(2,1,ia) = d2shape3D(1,2,ia)
    d2shape3D(3,1,ia) = d2shape3D(1,3,ia)
    d2shape3D(3,2,ia) = d2shape3D(2,3,ia)
  enddo

  end subroutine gf_shape3D_functions_2nd

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_shape3D_map_2nd(xelm,yelm,zelm,xi,eta,gamma,xyz,jinv,jacobian,djinv,ierr)

! gf_shape3D_map, plus the derivative of the inverse Jacobian
!
! djinv(:,:,a) = d jinv / d xi_a, packed like jinv (rows xi/eta/gamma,
! columns x/y/z). With J(d,b) = dx_d/dxi_b = SUM_ia dershape3D(b,ia) x_d(ia)
! and dJ_a(d,b) = SUM_ia d2shape3D(a,b,ia) x_d(ia),
!
!   d(J^-1)/dxi_a = -J^-1 (dJ/dxi_a) J^-1
!
! The first-order part is gf_shape3D_map itself, so the two cannot drift;
! on a degenerate element it returns the error and djinv is zero.

  use constants, only: NGNOD,NDIM,ZERO

  implicit none

  double precision, dimension(NGNOD), intent(in) :: xelm,yelm,zelm
  double precision, intent(in) :: xi,eta,gamma
  double precision, dimension(NDIM), intent(out) :: xyz
  double precision, dimension(NDIM,NDIM), intent(out) :: jinv
  double precision, intent(out) :: jacobian
  double precision, dimension(NDIM,NDIM,NDIM), intent(out) :: djinv
  integer, intent(out) :: ierr

  ! local parameters
  double precision, dimension(NDIM,NDIM,NGNOD) :: d2shape3D
  double precision, dimension(NDIM,NDIM) :: dj,tmp
  integer :: ia,a,b,d,m,n

  djinv(:,:,:) = ZERO

  call gf_shape3D_map(xelm,yelm,zelm,xi,eta,gamma,xyz,jinv,jacobian,ierr)
  if (ierr /= GF_OK) return

  call gf_shape3D_functions_2nd(xi,eta,gamma,d2shape3D)

  do a = 1,NDIM
    ! dJ_a(d,b) = SUM_ia d2shape3D(a,b,ia) x_d(ia)
    dj(:,:) = ZERO
    do ia = 1,NGNOD
      do b = 1,NDIM
        dj(1,b) = dj(1,b) + d2shape3D(a,b,ia)*xelm(ia)
        dj(2,b) = dj(2,b) + d2shape3D(a,b,ia)*yelm(ia)
        dj(3,b) = dj(3,b) + d2shape3D(a,b,ia)*zelm(ia)
      enddo
    enddo
    ! tmp = dJ_a J^-1  (x by x), then djinv_a = -J^-1 tmp  (xi by x)
    do n = 1,NDIM
      do d = 1,NDIM
        tmp(d,n) = ZERO
        do b = 1,NDIM
          tmp(d,n) = tmp(d,n) + dj(d,b)*jinv(b,n)
        enddo
      enddo
    enddo
    do n = 1,NDIM
      do m = 1,NDIM
        do d = 1,NDIM
          djinv(m,n,a) = djinv(m,n,a) - jinv(m,d)*tmp(d,n)
        enddo
      enddo
    enddo
  enddo

  end subroutine gf_shape3D_map_2nd

  end module gf_shape3D
