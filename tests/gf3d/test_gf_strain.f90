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
!---- test_gf_strain -- src/gf3d/gf_strain.F90 and src/gf3d/gf_moment.F90
!----
!---- Covers the whole geometric chain: shape functions -> Jacobian ->
!---- inversion -> the xi-to-x chain rule -> symmetrisation -> Voigt packing
!---- -> contraction with a moment tensor.
!----
!---- **The oracle differentiates in physical space.** For a manufactured
!---- displacement field u(x) written as a polynomial in x, the expected
!---- strain is the analytic symmetric gradient of that polynomial evaluated
!---- at the source point -- so no inverse Jacobian appears on the expected
!---- side at all. The code can only agree if every link in the chain is
!---- right, and a transposed or mis-scaled J^-1 has nowhere to hide.
!----
!---- Three exactness rules, and they are not the same
!---- ------------------------------------------------
!---- Getting these wrong makes a test either vacuous or falsely failing:
!----
!----   reference space (xi)        per-variable degree <= 4  (the basis span)
!----   physical, affine element    *total* degree <= 4
!----   physical, Q2-curved element *total* degree <= 2
!----
!---- because an affine pullback preserves total degree, while a Q2 map makes
!---- each x_j(xi) per-variable quadratic, so a product of two already
!---- reaches the edge of the basis. Note x^4 y^4 z^4 lies in Q4 but its
!---- affine pullback has total degree 12 and is not representable.
!----
!---- The sharpest assertions
!---- -----------------------
!---- * A **rigid rotation** u = W x with W antisymmetric must give eps == 0.
!----   An expected value of exactly zero is something a transposed J^-1
!----   cannot produce.
!---- * The **curved patch test**. On an affine element a wrong J^-1 -- a
!----   missing 1/det, transposed cofactors, dx/dxi where dxi/dx was meant --
!----   can still yield a *constant* wrong strain that a careless test lets
!----   through. With a position-dependent Jacobian a constant-strain field
!----   must still recover a constant eps, and a wrong inverse now varies
!----   from point to point and fails loudly.
!----

  program test_gf_strain

  use gf_par, only: GF_NCOMP,GF_OK

  use gf_shape3D, only: gf_shape3D_map

  use gf_interp, only: gf_interp_weights_deriv

  use gf_strain, only: gf_strain_dweights,gf_strain_snapshot, &
                       GF_VOIGT,GF_XX,GF_YY,GF_ZZ,GF_XY,GF_XZ,GF_YZ

  use gf_moment, only: gf_rotate_moment_tensor,gf_moment_contract, &
                       gf_moment_contract_full,gf_moment_unit_tensor

  use gf_manufactured, only: gf_lcg_seed,gf_rand_range,gf_report,gf_report_true

  use constants, only: NGNOD,NGLLX,NGLLY,NGLLZ,NDIM,GAUSSALPHA,GAUSSBETA,PI

  implicit none

  ! polynomial fields are stored as coefficients c(p, ix,iy,iz) of
  ! u_p(x) = SUM c * x^ix y^iy z^iz, with total degree at most MAXDEG
  integer, parameter :: MAXDEG = 4

  ! local parameters
  integer :: nfail,ierr,i,j,k,ia,p,q,v
  double precision, dimension(NDIM,NGNOD) :: xiref
  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision, dimension(NGLLZ) :: zigll,wzgll
  double precision, dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ) :: u
  double precision, dimension(GF_VOIGT,GF_NCOMP) :: eps
  double precision, dimension(GF_VOIGT) :: expected
  double precision, dimension(NDIM,NDIM) :: amat,bmat,wmat,m_cart,m_ref,eps_full
  double precision, dimension(NDIM) :: x0,xstar
  double precision, dimension(NDIM,0:MAXDEG,0:MAXDEG,0:MAXDEG) :: cf
  double precision :: xi,eta,gam,worst,err,scal,val1,val2
  double precision :: err_c1,err_c2,err_c3,ratio1,ratio2
  double precision, dimension(6) :: m_sph

  nfail = 0

  write(*,'(a)') ''
  write(*,'(a)') 'test_gf_strain'
  write(*,'(a)') ''

  call anchor_reference_coords(xiref)

  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  ! the point inside the element at which everything is evaluated
  xi = 0.3123d0 ; eta = -0.4571d0 ; gam = 0.6219d0

  !--------------------------------------------------------------------
  ! the affine reference element
  !
  ! A is deliberately non-symmetric and not a scalar multiple of an
  ! orthogonal matrix: with a symmetric or conformal A the rigid-rotation
  ! case below cannot distinguish J^-1 from its transpose, which is most of
  ! what that case exists to detect.
  !--------------------------------------------------------------------

  x0(1) = 0.11d0 ; x0(2) = -0.27d0 ; x0(3) = 0.43d0

  amat(1,1) = 1.30d0 ; amat(1,2) = 0.40d0 ; amat(1,3) = -0.20d0
  amat(2,1) = -0.30d0 ; amat(2,2) = 0.90d0 ; amat(2,3) = 0.50d0
  amat(3,1) = 0.20d0 ; amat(3,2) = -0.60d0 ; amat(3,3) = 1.10d0

  call build_affine_element(x0,amat,xiref,xelm,yelm,zelm)
  call map_affine(x0,amat,xi,eta,gam,xstar)

  !--------------------------------------------------------------------
  ! 1. a linear field: u = b + B x, so eps = (B + B^T)/2 exactly
  !
  ! B is non-symmetric and chosen so that all six Voigt slots come out
  ! distinct -- which is what pins the slot ordering, the factor 1/2 and
  ! the symmetrisation all at once. If any two were equal, a permutation of
  ! the Voigt slots would pass.
  !--------------------------------------------------------------------

  bmat(1,1) = 0.7d0  ; bmat(1,2) = 0.3d0  ; bmat(1,3) = -0.5d0
  bmat(2,1) = 0.9d0  ; bmat(2,2) = -1.1d0 ; bmat(2,3) = 0.2d0
  bmat(3,1) = 0.4d0  ; bmat(3,2) = -0.8d0 ; bmat(3,3) = 1.3d0

  cf(:,:,:,:) = 0.d0
  do p = 1,NDIM
    cf(p,1,0,0) = bmat(p,1)
    cf(p,0,1,0) = bmat(p,2)
    cf(p,0,0,1) = bmat(p,3)
  enddo
  cf(1,0,0,0) = 0.05d0
  cf(2,0,0,0) = -0.13d0
  cf(3,0,0,0) = 0.21d0

  call fill_field(cf,x0,amat,.false.,0.d0,xigll,yigll,zigll,u)
  call lib_strain(xelm,yelm,zelm,u,xi,eta,gam,eps,ierr)
  call gf_report_true('linear field: strain evaluated    ',ierr == GF_OK,nfail)

  call analytic_strain(cf,xstar,expected)

  ! all six slots distinct, or the ordering is not actually being tested
  call gf_report_true('  ... six Voigt slots are distinct',all_distinct(expected),nfail)

  worst = 0.d0
  do ia = 1,GF_NCOMP
    do v = 1,GF_VOIGT
      worst = max(worst,abs(eps(v,ia) - expected(v)))
    enddo
  enddo
  call gf_report('affine, u = b + Bx                ',worst,1.d-12,nfail)

  !--------------------------------------------------------------------
  ! 2. rigid rotation: eps must vanish
  !--------------------------------------------------------------------

  wmat(:,:) = 0.d0
  wmat(1,2) = 0.6d0  ; wmat(2,1) = -0.6d0
  wmat(1,3) = -0.4d0 ; wmat(3,1) = 0.4d0
  wmat(2,3) = 0.9d0  ; wmat(3,2) = -0.9d0

  cf(:,:,:,:) = 0.d0
  do p = 1,NDIM
    cf(p,1,0,0) = wmat(p,1)
    cf(p,0,1,0) = wmat(p,2)
    cf(p,0,0,1) = wmat(p,3)
  enddo

  call fill_field(cf,x0,amat,.false.,0.d0,xigll,yigll,zigll,u)
  call lib_strain(xelm,yelm,zelm,u,xi,eta,gam,eps,ierr)

  worst = 0.d0
  do ia = 1,GF_NCOMP
    do v = 1,GF_VOIGT
      worst = max(worst,abs(eps(v,ia)))
    enddo
  enddo
  call gf_report('rigid rotation gives eps = 0      ',worst,1.d-14,nfail)

  !--------------------------------------------------------------------
  ! 3. a random field of total degree <= 4 on the affine element
  !
  ! Total degree, not per-variable: the affine pullback of a total-degree-4
  ! polynomial has total degree 4 in xi and so lies inside the basis, while
  ! the pullback of x^4y^4z^4 would have total degree 12 and would not.
  !--------------------------------------------------------------------

  call gf_lcg_seed(2718281)
  call random_poly(cf,MAXDEG,1.d0)

  call fill_field(cf,x0,amat,.false.,0.d0,xigll,yigll,zigll,u)
  call lib_strain(xelm,yelm,zelm,u,xi,eta,gam,eps,ierr)
  call analytic_strain(cf,xstar,expected)

  worst = 0.d0
  do ia = 1,GF_NCOMP
    do v = 1,GF_VOIGT
      worst = max(worst,abs(eps(v,ia) - expected(v)))
    enddo
  enddo
  call gf_report('affine, random total-degree-4     ',worst,1.d-12,nfail)

  !--------------------------------------------------------------------
  ! 4. the same, at a realistic scale
  !
  ! |x0| = 1 and an element half-width of 1e-2: a unit-radius planet with a
  ! typical element. The looser tolerance is the explainable cost of the
  ! SUM dershape3D = 0 cancellation against an offset origin, and is worth
  ! logging because Stage 8's second derivative amplifies exactly this.
  !--------------------------------------------------------------------

  x0(1) = 0.5773502691896258d0
  x0(2) = 0.5773502691896258d0
  x0(3) = 0.5773502691896258d0
  scal = 1.d-2

  call build_affine_element(x0,scal*amat,xiref,xelm,yelm,zelm)
  call map_affine(x0,scal*amat,xi,eta,gam,xstar)

  call gf_lcg_seed(161803398)
  call random_poly(cf,MAXDEG,1.d0)

  call fill_field(cf,x0,scal*amat,.false.,0.d0,xigll,yigll,zigll,u)
  call lib_strain(xelm,yelm,zelm,u,xi,eta,gam,eps,ierr)
  call analytic_strain(cf,xstar,expected)

  worst = 0.d0
  do ia = 1,GF_NCOMP
    do v = 1,GF_VOIGT
      worst = max(worst,abs(eps(v,ia) - expected(v)))
    enddo
  enddo
  call gf_report('affine at realistic scale (1e-2)  ',worst,1.d-10,nfail)

  !--------------------------------------------------------------------
  ! 5. curved Q2 element: the patch test
  !
  ! A constant-strain field on an element whose Jacobian varies from point
  ! to point. This is the most diagnostic single case in the suite: a wrong
  ! inverse Jacobian that happened to produce a plausible constant on an
  ! affine element now varies with position and cannot.
  !--------------------------------------------------------------------

  x0(1) = 0.11d0 ; x0(2) = -0.27d0 ; x0(3) = 0.43d0

  call build_curved_element(x0,amat,0.15d0,xiref,xelm,yelm,zelm)

  cf(:,:,:,:) = 0.d0
  do p = 1,NDIM
    cf(p,1,0,0) = bmat(p,1)
    cf(p,0,1,0) = bmat(p,2)
    cf(p,0,0,1) = bmat(p,3)
  enddo

  call fill_field(cf,x0,amat,.true.,0.15d0,xigll,yigll,zigll,u)

  ! sampled at several points: a constant strain must come back constant
  worst = 0.d0
  do k = -1,1
    do j = -1,1
      do i = -1,1
        call map_curved(x0,amat,0.15d0,0.4d0*dble(i),0.4d0*dble(j),0.4d0*dble(k),xstar)
        call lib_strain(xelm,yelm,zelm,u,0.4d0*dble(i),0.4d0*dble(j),0.4d0*dble(k),eps,ierr)
        if (ierr /= GF_OK) cycle
        call analytic_strain(cf,xstar,expected)
        do ia = 1,GF_NCOMP
          do v = 1,GF_VOIGT
            worst = max(worst,abs(eps(v,ia) - expected(v)))
          enddo
        enddo
      enddo
    enddo
  enddo
  call gf_report('curved Q2 patch test, 27 points   ',worst,1.d-11,nfail)

  !--------------------------------------------------------------------
  ! 6. the 27 anchors reproduce all 125 GLL coordinates of a Q2 geometry
  !
  ! The manufactured analogue of tests/gf3d/test_gf_anchors.f90, which needs
  ! a database. Here everything is double precision, so this asserts the
  ! geometric premise itself rather than the float32 storage floor -- and it
  ! runs in CI on every commit.
  !--------------------------------------------------------------------

  worst = 0.d0
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        call map_curved(x0,amat,0.15d0,xigll(i),yigll(j),zigll(k),xstar)
        call eval_map_from_anchors(xelm,yelm,zelm,xigll(i),yigll(j),zigll(k),eps_full(:,1),ierr)
        if (ierr /= GF_OK) cycle
        do p = 1,NDIM
          worst = max(worst,abs(eps_full(p,1) - xstar(p)))
        enddo
      enddo
    enddo
  enddo
  call gf_report('27 anchors reproduce 125 GLL (Q2) ',worst,1.d-14,nfail)

  !--------------------------------------------------------------------
  ! 7. the Voigt contraction against the plain 3x3 double sum
  !
  ! This is what pins the factor-2 weights on the off-diagonals, which is
  ! the easiest thing in a moment-tensor contraction to get wrong.
  !--------------------------------------------------------------------

  call gf_lcg_seed(577215)
  do q = 1,NDIM
    do p = 1,NDIM
      m_cart(p,q) = 0.d0
    enddo
  enddo
  do p = 1,NDIM
    do q = p,NDIM
      m_cart(p,q) = gf_rand_range(-1.d0,1.d0)
      m_cart(q,p) = m_cart(p,q)
    enddo
  enddo

  do v = 1,GF_VOIGT
    expected(v) = gf_rand_range(-1.d0,1.d0)
  enddo
  eps_full(1,1) = expected(GF_XX)
  eps_full(2,2) = expected(GF_YY)
  eps_full(3,3) = expected(GF_ZZ)
  eps_full(1,2) = expected(GF_XY) ; eps_full(2,1) = expected(GF_XY)
  eps_full(1,3) = expected(GF_XZ) ; eps_full(3,1) = expected(GF_XZ)
  eps_full(2,3) = expected(GF_YZ) ; eps_full(3,2) = expected(GF_YZ)

  call gf_moment_contract(m_cart,expected,val1)
  call gf_moment_contract_full(m_cart,eps_full,val2)
  call gf_report('Voigt contraction == full 3x3 sum ',abs(val1-val2)/abs(val2),1.d-14,nfail)

  !--------------------------------------------------------------------
  ! 8. linearity in M -- the identity Stage 6's partials rest on
  !
  ! SUM_v M_v * contract(unit_v, eps) must reproduce contract(M, eps). If it
  ! does, the six moment-tensor partial derivatives are the same strain
  ! re-contracted, exactly, and Stage 6 is a loop rather than a derivation.
  !--------------------------------------------------------------------

  val2 = 0.d0
  do v = 1,GF_VOIGT
    call gf_moment_unit_tensor(v,m_ref)
    call gf_moment_contract(m_ref,expected,val1)
    select case (v)
    case (GF_XX) ; val2 = val2 + m_cart(1,1)*val1
    case (GF_YY) ; val2 = val2 + m_cart(2,2)*val1
    case (GF_ZZ) ; val2 = val2 + m_cart(3,3)*val1
    case (GF_XY) ; val2 = val2 + m_cart(1,2)*val1
    case (GF_XZ) ; val2 = val2 + m_cart(1,3)*val1
    case (GF_YZ) ; val2 = val2 + m_cart(2,3)*val1
    end select
  enddo
  call gf_moment_contract(m_cart,expected,val1)
  call gf_report('linearity: SUM M_v dp(v) == full  ',abs(val1-val2)/abs(val1),1.d-14,nfail)

  !--------------------------------------------------------------------
  ! 9. the spherical-to-Cartesian moment tensor rotation
  !
  ! Checked against an independent route: M_cart = SUM_ab M_ab e_a (x) e_b
  ! built from the spherical basis vectors, a matrix construction that
  ! shares no code with the expanded scalar form in gf_moment. A rotation
  ! also preserves the trace and the Frobenius norm, so both are asserted --
  ! a single wrong sign in an off-diagonal breaks the norm.
  !--------------------------------------------------------------------

  call gf_lcg_seed(1414213)
  do v = 1,6
    m_sph(v) = gf_rand_range(-1.d0,1.d0)
  enddo

  worst = 0.d0
  do i = 1,20
    ! avoid the poles, where the spherical basis is degenerate and reduce()
    ! is what the production path relies on
    xi  = gf_rand_range(0.15d0,PI-0.15d0)
    eta = gf_rand_range(0.d0,2.d0*PI)

    call gf_rotate_moment_tensor(xi,eta,m_sph,m_cart)
    call rotate_reference(xi,eta,m_sph,m_ref)

    do q = 1,NDIM
      do p = 1,NDIM
        worst = max(worst,abs(m_cart(p,q) - m_ref(p,q)))
      enddo
    enddo
  enddo
  call gf_report('MT rotation vs basis-vector form  ',worst,1.d-14,nfail)

  ! invariants, on the last drawn orientation
  err = abs((m_cart(1,1)+m_cart(2,2)+m_cart(3,3)) - (m_sph(1)+m_sph(2)+m_sph(3)))
  call gf_report('MT rotation preserves the trace   ',err,1.d-14,nfail)

  val1 = 0.d0
  do q = 1,NDIM
    do p = 1,NDIM
      val1 = val1 + m_cart(p,q)**2
    enddo
  enddo
  val2 = m_sph(1)**2 + m_sph(2)**2 + m_sph(3)**2 &
       + 2.d0*(m_sph(4)**2 + m_sph(5)**2 + m_sph(6)**2)
  call gf_report('MT rotation preserves |M|_F       ',abs(val1-val2)/abs(val2),1.d-14,nfail)

  !--------------------------------------------------------------------
  ! 10. a cubic field on a curved element is NOT exact -- a convergence check
  !
  ! Labelled convergence rather than accuracy so that nobody later reads its
  ! residual as a bug: the error is order 1 here, and correctly so.
  !
  ! The rate is **second** order in the curvature amplitude, not first, and
  ! that is worth writing down because testing.md guessed first order and a
  ! [0.4,0.6] band fails. Expand u(x(xi)) with x = L(xi) + a*Q(xi), L affine
  ! and Q quadratic, u cubic in x:
  !
  !   a^0 : u(L)              per-variable degree 3  -> in the basis
  !   a^1 : grad u . (a Q)    degree 2 + 2 = 4       -> still in the basis
  !   a^2 : hess u : (a Q)^2  degree 1 + 4 = 5       -> NOT representable
  !
  ! So the first term the degree-4 basis cannot reproduce is O(a^2), and
  ! halving the amplitude quarters the error. Three amplitudes are used
  ! rather than two, so the *rate* is measured rather than assumed.
  !--------------------------------------------------------------------

  cf(:,:,:,:) = 0.d0
  call gf_lcg_seed(1732050)
  do p = 1,NDIM
    cf(p,3,0,0) = gf_rand_range(-1.d0,1.d0)
    cf(p,0,3,0) = gf_rand_range(-1.d0,1.d0)
    cf(p,0,0,3) = gf_rand_range(-1.d0,1.d0)
    cf(p,2,1,0) = gf_rand_range(-1.d0,1.d0)
    cf(p,1,0,2) = gf_rand_range(-1.d0,1.d0)
  enddo

  call curved_error(cf,x0,amat,0.20d0,xiref,xigll,yigll,zigll,xi,eta,gam,err_c1)
  call curved_error(cf,x0,amat,0.10d0,xiref,xigll,yigll,zigll,xi,eta,gam,err_c2)
  call curved_error(cf,x0,amat,0.05d0,xiref,xigll,yigll,zigll,xi,eta,gam,err_c3)

  ratio1 = err_c2 / err_c1
  ratio2 = err_c3 / err_c2

  write(*,'(a,3es12.5)') '     curved cubic error, amp 0.20/0.10/0.05 = ', &
                         err_c1,err_c2,err_c3
  write(*,'(a,2f8.4,a)') '     successive ratios = ',ratio1,ratio2, &
                         '   (second order = 0.25)'

  ! the band brackets 0.25; the coarsest amplitude still carries visible
  ! third-order contamination, which is why the first ratio sits above it
  call gf_report_true('cubic on curved: 2nd order in amp ', &
                      ratio1 > 0.18d0 .and. ratio1 < 0.36d0 .and. &
                      ratio2 > 0.18d0 .and. ratio2 < 0.36d0,nfail)

  !--------------------------------------------------------------------
  ! 11. a degenerate element is refused, not evaluated
  !--------------------------------------------------------------------

  do ia = 1,NGNOD
    xelm(ia) = 0.5d0 ; yelm(ia) = 0.5d0 ; zelm(ia) = 0.5d0
  enddo
  call lib_strain(xelm,yelm,zelm,u,0.1d0,0.2d0,0.3d0,eps,ierr)
  call gf_report_true('degenerate element returns ierr   ',ierr /= GF_OK,nfail)

  !--------------------------------------------------------------------

  write(*,'(a)') ''
  if (nfail /= 0) then
    write(*,'(a,i0,a)') 'test_gf_strain: ',nfail,' assertion(s) FAILED'
    stop 1
  endif

  write(*,'(a)') 'test_gf_strain: all assertions passed'

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine anchor_reference_coords(xr)

! reference (xi,eta,gamma) of the 27 anchors, from hex_nodes() rather than a
! hand table, so a reordering of the control nodes breaks this test

  implicit none

  double precision, dimension(NDIM,NGNOD), intent(out) :: xr

  ! local parameters
  integer, dimension(NGNOD) :: iaddx,iaddy,iaddz
  integer :: ial

  call hex_nodes(iaddx,iaddy,iaddz)

  do ial = 1,NGNOD
    xr(1,ial) = dble(iaddx(ial)) - 1.d0
    xr(2,ial) = dble(iaddy(ial)) - 1.d0
    xr(3,ial) = dble(iaddz(ial)) - 1.d0
  enddo

  end subroutine anchor_reference_coords

!
!-------------------------------------------------------------------------------------------------
!

  subroutine map_affine(x0l,al,xl,el,gl,xout)

! x = x0 + A xi

  implicit none

  double precision, dimension(NDIM), intent(in) :: x0l
  double precision, dimension(NDIM,NDIM), intent(in) :: al
  double precision, intent(in) :: xl,el,gl
  double precision, dimension(NDIM), intent(out) :: xout

  ! local parameters
  integer :: pl

  do pl = 1,NDIM
    xout(pl) = x0l(pl) + al(pl,1)*xl + al(pl,2)*el + al(pl,3)*gl
  enddo

  end subroutine map_affine

!
!-------------------------------------------------------------------------------------------------
!

  subroutine map_curved(x0l,al,amp,xl,el,gl,xout)

! a genuinely Q2 map: affine plus quadratic terms
!
! Every term has per-variable degree at most 2, so the tri-quadratic shape
! functions through the 27 anchors reproduce this map *exactly* -- which is
! what makes case 6 an assertion about the geometry rather than about
! interpolation error.

  implicit none

  double precision, dimension(NDIM), intent(in) :: x0l
  double precision, dimension(NDIM,NDIM), intent(in) :: al
  double precision, intent(in) :: amp,xl,el,gl
  double precision, dimension(NDIM), intent(out) :: xout

  call map_affine(x0l,al,xl,el,gl,xout)

  xout(1) = xout(1) + amp*(0.7d0*xl*xl - 0.4d0*el*el + 0.5d0*xl*el)
  xout(2) = xout(2) + amp*(-0.6d0*el*el + 0.3d0*gl*gl + 0.2d0*el*gl)
  xout(3) = xout(3) + amp*(0.9d0*gl*gl - 0.5d0*xl*xl + 0.4d0*xl*gl)

  end subroutine map_curved

!
!-------------------------------------------------------------------------------------------------
!

  subroutine build_affine_element(x0l,al,xr,xe,ye,ze)

! the 27 anchor positions of an affine element

  implicit none

  double precision, dimension(NDIM), intent(in) :: x0l
  double precision, dimension(NDIM,NDIM), intent(in) :: al
  double precision, dimension(NDIM,NGNOD), intent(in) :: xr
  double precision, dimension(NGNOD), intent(out) :: xe,ye,ze

  ! local parameters
  integer :: ial
  double precision, dimension(NDIM) :: xx

  do ial = 1,NGNOD
    call map_affine(x0l,al,xr(1,ial),xr(2,ial),xr(3,ial),xx)
    xe(ial) = xx(1) ; ye(ial) = xx(2) ; ze(ial) = xx(3)
  enddo

  end subroutine build_affine_element

!
!-------------------------------------------------------------------------------------------------
!

  subroutine build_curved_element(x0l,al,amp,xr,xe,ye,ze)

! the 27 anchor positions of a curved Q2 element

  implicit none

  double precision, dimension(NDIM), intent(in) :: x0l
  double precision, dimension(NDIM,NDIM), intent(in) :: al
  double precision, intent(in) :: amp
  double precision, dimension(NDIM,NGNOD), intent(in) :: xr
  double precision, dimension(NGNOD), intent(out) :: xe,ye,ze

  ! local parameters
  integer :: ial
  double precision, dimension(NDIM) :: xx

  do ial = 1,NGNOD
    call map_curved(x0l,al,amp,xr(1,ial),xr(2,ial),xr(3,ial),xx)
    xe(ial) = xx(1) ; ye(ial) = xx(2) ; ze(ial) = xx(3)
  enddo

  end subroutine build_curved_element

!
!-------------------------------------------------------------------------------------------------
!

  subroutine fill_field(c,x0l,al,curved,amp,xg,yg,zg,uu)

! samples u_p(x) at the 125 GLL points of the element
!
! The same polynomial is used for all three force components -- the (a)
! index is exercised by test_gf_interp, and keeping it uniform here means an
! assertion failure points at the geometry rather than at bookkeeping.

  implicit none

  double precision, dimension(NDIM,0:MAXDEG,0:MAXDEG,0:MAXDEG), intent(in) :: c
  double precision, dimension(NDIM), intent(in) :: x0l
  double precision, dimension(NDIM,NDIM), intent(in) :: al
  logical, intent(in) :: curved
  double precision, intent(in) :: amp
  double precision, dimension(NGLLX), intent(in) :: xg
  double precision, dimension(NGLLY), intent(in) :: yg
  double precision, dimension(NGLLZ), intent(in) :: zg
  double precision, dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ), intent(out) :: uu

  ! local parameters
  integer :: il,jl,kl,al2,pl
  double precision, dimension(NDIM) :: xx

  do kl = 1,NGLLZ
    do jl = 1,NGLLY
      do il = 1,NGLLX
        if (curved) then
          call map_curved(x0l,al,amp,xg(il),yg(jl),zg(kl),xx)
        else
          call map_affine(x0l,al,xg(il),yg(jl),zg(kl),xx)
        endif
        do pl = 1,NDIM
          do al2 = 1,GF_NCOMP
            uu(al2,pl,il,jl,kl) = poly_val(c(pl,:,:,:),xx)
          enddo
        enddo
      enddo
    enddo
  enddo

  end subroutine fill_field

!
!-------------------------------------------------------------------------------------------------
!

  double precision function poly_val(c,x)

! SUM c(ix,iy,iz) x^ix y^iy z^iz

  implicit none

  double precision, dimension(0:MAXDEG,0:MAXDEG,0:MAXDEG), intent(in) :: c
  double precision, dimension(NDIM), intent(in) :: x

  ! local parameters
  integer :: ix,iy,iz
  double precision :: s

  s = 0.d0
  do iz = 0,MAXDEG
    do iy = 0,MAXDEG
      do ix = 0,MAXDEG
        if (c(ix,iy,iz) == 0.d0) cycle
        s = s + c(ix,iy,iz) * x(1)**ix * x(2)**iy * x(3)**iz
      enddo
    enddo
  enddo

  poly_val = s

  end function poly_val

!
!-------------------------------------------------------------------------------------------------
!

  double precision function poly_deriv(c,x,idir)

! d/dx_idir of the same polynomial, differentiated analytically in
! *physical* space -- this is the whole point: no Jacobian appears here

  implicit none

  double precision, dimension(0:MAXDEG,0:MAXDEG,0:MAXDEG), intent(in) :: c
  double precision, dimension(NDIM), intent(in) :: x
  integer, intent(in) :: idir

  ! local parameters
  integer :: ix,iy,iz
  double precision :: s

  s = 0.d0
  do iz = 0,MAXDEG
    do iy = 0,MAXDEG
      do ix = 0,MAXDEG
        if (c(ix,iy,iz) == 0.d0) cycle
        select case (idir)
        case (1)
          if (ix > 0) s = s + dble(ix)*c(ix,iy,iz) * x(1)**(ix-1) * x(2)**iy * x(3)**iz
        case (2)
          if (iy > 0) s = s + dble(iy)*c(ix,iy,iz) * x(1)**ix * x(2)**(iy-1) * x(3)**iz
        case (3)
          if (iz > 0) s = s + dble(iz)*c(ix,iy,iz) * x(1)**ix * x(2)**iy * x(3)**(iz-1)
        end select
      enddo
    enddo
  enddo

  poly_deriv = s

  end function poly_deriv

!
!-------------------------------------------------------------------------------------------------
!

  subroutine analytic_strain(c,x,e)

! the expected Voigt strain: the symmetric gradient, in physical space

  implicit none

  double precision, dimension(NDIM,0:MAXDEG,0:MAXDEG,0:MAXDEG), intent(in) :: c
  double precision, dimension(NDIM), intent(in) :: x
  double precision, dimension(GF_VOIGT), intent(out) :: e

  ! local parameters
  double precision, dimension(NDIM,NDIM) :: gr
  integer :: pl,ql

  do ql = 1,NDIM
    do pl = 1,NDIM
      gr(pl,ql) = poly_deriv(c(pl,:,:,:),x,ql)
    enddo
  enddo

  e(GF_XX) = gr(1,1)
  e(GF_YY) = gr(2,2)
  e(GF_ZZ) = gr(3,3)
  e(GF_XY) = 0.5d0*(gr(1,2) + gr(2,1))
  e(GF_XZ) = 0.5d0*(gr(1,3) + gr(3,1))
  e(GF_YZ) = 0.5d0*(gr(2,3) + gr(3,2))

  end subroutine analytic_strain

!
!-------------------------------------------------------------------------------------------------
!

  subroutine random_poly(c,ndeg,amp)

! a random polynomial of *total* degree at most ndeg

  implicit none

  double precision, dimension(NDIM,0:MAXDEG,0:MAXDEG,0:MAXDEG), intent(out) :: c
  integer, intent(in) :: ndeg
  double precision, intent(in) :: amp

  ! local parameters
  integer :: ix,iy,iz,pl

  c(:,:,:,:) = 0.d0
  do pl = 1,NDIM
    do iz = 0,ndeg
      do iy = 0,ndeg
        do ix = 0,ndeg
          if (ix + iy + iz > ndeg) cycle
          c(pl,ix,iy,iz) = amp*gf_rand_range(-1.d0,1.d0)
        enddo
      enddo
    enddo
  enddo

  end subroutine random_poly

!
!-------------------------------------------------------------------------------------------------
!

  subroutine lib_strain(xe,ye,ze,uu,xl,el,gl,e,ier)

! the library's whole geometric chain, at one point
!
! shape functions -> Jacobian -> inversion -> chain rule -> symmetrisation
! -> Voigt packing. Anything wrong anywhere in it shows up in `e`.

  implicit none

  double precision, dimension(NGNOD), intent(in) :: xe,ye,ze
  double precision, dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ), intent(in) :: uu
  double precision, intent(in) :: xl,el,gl
  double precision, dimension(GF_VOIGT,GF_NCOMP), intent(out) :: e
  integer, intent(out) :: ier

  ! local parameters
  double precision, dimension(NDIM) :: xyzl
  double precision, dimension(NDIM,NDIM) :: jinvl
  double precision, dimension(NGLLX) :: hxil,hpxil
  double precision, dimension(NGLLY) :: hetal,hpetal
  double precision, dimension(NGLLZ) :: hgaml,hpgaml
  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM) :: dwl
  double precision :: jacl

  e(:,:) = 0.d0

  call gf_shape3D_map(xe,ye,ze,xl,el,gl,xyzl,jinvl,jacl,ier)
  if (ier /= GF_OK) return

  call gf_interp_weights_deriv(xl,el,gl,hxil,hpxil,hetal,hpetal,hgaml,hpgaml)
  call gf_strain_dweights(hxil,hpxil,hetal,hpetal,hgaml,hpgaml,jinvl,dwl)
  call gf_strain_snapshot(uu,dwl,e)

  end subroutine lib_strain

!
!-------------------------------------------------------------------------------------------------
!

  subroutine eval_map_from_anchors(xe,ye,ze,xl,el,gl,xout,ier)

! the element map evaluated through the 27 anchors

  implicit none

  double precision, dimension(NGNOD), intent(in) :: xe,ye,ze
  double precision, intent(in) :: xl,el,gl
  double precision, dimension(NDIM), intent(out) :: xout
  integer, intent(out) :: ier

  ! local parameters
  double precision, dimension(NDIM,NDIM) :: jinvl
  double precision :: jacl

  call gf_shape3D_map(xe,ye,ze,xl,el,gl,xout,jinvl,jacl,ier)

  end subroutine eval_map_from_anchors

!
!-------------------------------------------------------------------------------------------------
!

  subroutine curved_error(c,x0l,al,amp,xr,xg,yg,zg,xl,el,gl,errout)

! strain error of a cubic-in-x field on a curved element of given curvature

  implicit none

  double precision, dimension(NDIM,0:MAXDEG,0:MAXDEG,0:MAXDEG), intent(in) :: c
  double precision, dimension(NDIM), intent(in) :: x0l
  double precision, dimension(NDIM,NDIM), intent(in) :: al
  double precision, intent(in) :: amp
  double precision, dimension(NDIM,NGNOD), intent(in) :: xr
  double precision, dimension(NGLLX), intent(in) :: xg
  double precision, dimension(NGLLY), intent(in) :: yg
  double precision, dimension(NGLLZ), intent(in) :: zg
  double precision, intent(in) :: xl,el,gl
  double precision, intent(out) :: errout

  ! local parameters
  double precision, dimension(NGNOD) :: xe,ye,ze
  double precision, dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ) :: uu
  double precision, dimension(GF_VOIGT,GF_NCOMP) :: e
  double precision, dimension(GF_VOIGT) :: ex
  double precision, dimension(NDIM) :: xx
  integer :: ier2,vl

  call build_curved_element(x0l,al,amp,xr,xe,ye,ze)
  call fill_field(c,x0l,al,.true.,amp,xg,yg,zg,uu)
  call map_curved(x0l,al,amp,xl,el,gl,xx)

  call lib_strain(xe,ye,ze,uu,xl,el,gl,e,ier2)
  call analytic_strain(c,xx,ex)

  errout = 0.d0
  do vl = 1,GF_VOIGT
    errout = max(errout,abs(e(vl,1) - ex(vl)))
  enddo

  end subroutine curved_error

!
!-------------------------------------------------------------------------------------------------
!

  subroutine rotate_reference(th,ph,ms,mc)

! M_cart = SUM_ab M_sph(a,b) e_a (x) e_b, from the spherical basis vectors
!
! An independent route to the same rotation: this builds the basis vectors
! and forms outer products, where gf_rotate_moment_tensor writes the
! expanded scalar form transcribed from locate_sources.f90. They share no
! code, which is what makes this an oracle.

  implicit none

  double precision, intent(in) :: th,ph
  double precision, dimension(6), intent(in) :: ms
  double precision, dimension(NDIM,NDIM), intent(out) :: mc

  ! local parameters
  double precision, dimension(NDIM,NDIM) :: ev
  double precision, dimension(NDIM,NDIM) :: msm
  integer :: al2,bl,pl,ql

  ! rows: e_r, e_theta, e_phi in Cartesian
  ev(1,1) = sin(th)*cos(ph) ; ev(1,2) = sin(th)*sin(ph) ; ev(1,3) = cos(th)
  ev(2,1) = cos(th)*cos(ph) ; ev(2,2) = cos(th)*sin(ph) ; ev(2,3) = -sin(th)
  ev(3,1) = -sin(ph)        ; ev(3,2) = cos(ph)         ; ev(3,3) = 0.d0

  ! (Mrr,Mtt,Mpp,Mrt,Mrp,Mtp) as a symmetric matrix in (r,theta,phi)
  msm(1,1) = ms(1) ; msm(2,2) = ms(2) ; msm(3,3) = ms(3)
  msm(1,2) = ms(4) ; msm(2,1) = ms(4)
  msm(1,3) = ms(5) ; msm(3,1) = ms(5)
  msm(2,3) = ms(6) ; msm(3,2) = ms(6)

  mc(:,:) = 0.d0
  do bl = 1,NDIM
    do al2 = 1,NDIM
      do ql = 1,NDIM
        do pl = 1,NDIM
          mc(pl,ql) = mc(pl,ql) + msm(al2,bl)*ev(al2,pl)*ev(bl,ql)
        enddo
      enddo
    enddo
  enddo

  end subroutine rotate_reference

!
!-------------------------------------------------------------------------------------------------
!

  logical function all_distinct(e)

! true when the six Voigt slots hold six different values
!
! Guards against a manufactured field that would let a permutation of the
! Voigt slots pass unnoticed.

  implicit none

  double precision, dimension(GF_VOIGT), intent(in) :: e

  ! local parameters
  integer :: il,jl

  all_distinct = .true.
  do jl = 2,GF_VOIGT
    do il = 1,jl-1
      if (abs(e(il) - e(jl)) < 1.d-6) all_distinct = .false.
    enddo
  enddo

  end function all_distinct

  end program test_gf_strain
