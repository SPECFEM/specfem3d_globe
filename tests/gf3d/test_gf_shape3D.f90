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
!---- test_gf_shape3D -- src/gf3d/gf_shape3D.F90
!----
!---- gf_shape3D_map() is a fork of src/shared/recompute_jacobian.f90, made
!---- because the original ends a degenerate element with `stop`, which
!---- inside a shared object loaded by Python kills the interpreter. This
!---- test is the payoff for forking rather than reimplementing: the
!---- original is still in the tree, so it is a free oracle.
!----
!---- The fork transcribes the original's expression and accumulation order
!---- deliberately, so under a value-safe floating-point model the two agree
!---- bit for bit, and this test used to assert exactly that. It no longer
!---- does, because the compiler is not obliged to honour that order: ifort
!---- and ifx at -O3 -xHost (flags.guess sets no -fp-model) contract
!---- multiply-adds into FMAs and vectorise the 27-term reductions, and they
!---- do so *differently* for the two compilation units -- a module
!---- procedure whose shape functions arrive from a call, against an
!---- external with intent(inout) scalars. CI measured 9849 of 10000 trials
!---- differing in the last bits under ifort, with every closed-form
!---- identity intact; a structural error (a transposed jinv packing, a
!---- wrong cofactor) would have failed all 10000.
!----
!---- So the assertion is agreement to a tolerance derived from the
!---- arithmetic (see section 3), and the bit-for-bit mismatch count is
!---- still *reported*: it reads zero under gfortran and non-zero under an
!---- FMA-contracting build, which is useful to know and costs nothing.
!----
!---- Ordering matters in the loop below: gf_shape3D_map() is called FIRST
!---- and the oracle only when it reports success. recompute_jacobian.f90:261
!---- stops on a non-positive Jacobian, so feeding it a degenerate element --
!---- which the degenerate-geometry cases below do on purpose -- would take
!---- the test process down with it rather than failing an assertion.
!----

  program test_gf_shape3D

  use gf_par, only: GF_OK,GF_ERR_GEOMETRY
  use gf_shape3D, only: gf_shape3D_functions,gf_shape3D_map

  use gf_manufactured, only: gf_lcg_seed,gf_rand_range,gf_report,gf_report_true

  use constants, only: NGNOD,NDIM

  implicit none

  ! number of random (anchor set, xi/eta/gamma) pairs
  integer, parameter :: NTRIAL = 10000

  ! local parameters
  integer :: nfail,itrial,ia,i,j,nskipped,nmismatch,ierr
  double precision, dimension(NDIM,NGNOD) :: xiref
  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  double precision, dimension(NGNOD) :: shape3D
  double precision, dimension(NDIM,NGNOD) :: dershape3D
  double precision, dimension(NDIM) :: xyz,x0
  double precision, dimension(NDIM,NDIM) :: jinv,amat,jref
  double precision :: jacobian
  double precision :: xi,eta,gamma
  double precision :: x,y,z,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  double precision :: det,worst_shape,worst_der,s,sd
  double precision :: xscale,worst_xyz,worst_jinv

  nfail = 0

  write(*,'(a)') ''
  write(*,'(a)') 'test_gf_shape3D'
  write(*,'(a)') ''
  write(*,'(a,i0)') '  NGNOD = ',NGNOD
  write(*,'(a,i0)') '  random trials = ',NTRIAL
  write(*,'(a)') ''

  if (NGNOD /= 27) then
    write(*,'(a)') '  this test requires NGNOD = 27'
    stop 1
  endif

  call anchor_reference_coords(xiref)

  !--------------------------------------------------------------------
  ! 1. shape function identities
  !
  ! Partition of unity and the vanishing derivative sums hold for any
  ! interpolatory basis and cost nothing to check, but they localise a
  ! typo in one of the 27 nodes to the basis rather than to the Jacobian.
  !--------------------------------------------------------------------

  call gf_lcg_seed(20260904)

  worst_shape = 0.d0
  worst_der = 0.d0
  do itrial = 1,1000
    xi    = gf_rand_range(-1.1d0,1.1d0)
    eta   = gf_rand_range(-1.1d0,1.1d0)
    gamma = gf_rand_range(-1.1d0,1.1d0)

    call gf_shape3D_functions(xi,eta,gamma,shape3D,dershape3D)

    s = sum(shape3D)
    worst_shape = max(worst_shape,abs(s - 1.d0))
    do i = 1,NDIM
      sd = sum(dershape3D(i,:))
      worst_der = max(worst_der,abs(sd))
    enddo
  enddo

  call gf_report('partition of unity, sum h = 1     ',worst_shape,1.d-14,nfail)
  call gf_report('derivative sums, sum dh = 0       ',worst_der,1.d-14,nfail)

  !--------------------------------------------------------------------
  ! 2. the shape functions are interpolatory at the 27 anchors
  !
  ! h_a(xi_b) = delta_ab. This pins the association between a node's
  ! position in the 1..27 ordering and the basis function that carries it,
  ! which is what lets gf_gather_anchors() use hex_nodes_anchor_ijk()
  ! indices against gf_shape3D_functions() without a translation table.
  !--------------------------------------------------------------------

  worst_shape = 0.d0
  do ia = 1,NGNOD
    call gf_shape3D_functions(xiref(1,ia),xiref(2,ia),xiref(3,ia),shape3D,dershape3D)
    do i = 1,NGNOD
      if (i == ia) then
        worst_shape = max(worst_shape,abs(shape3D(i) - 1.d0))
      else
        worst_shape = max(worst_shape,abs(shape3D(i)))
      endif
    enddo
  enddo

  call gf_report('Kronecker delta at the 27 anchors ',worst_shape,1.d-14,nfail)

  !--------------------------------------------------------------------
  ! 3. against recompute_jacobian
  !
  ! Two relative errors are asserted, each against a bound derived from
  ! the arithmetic the two routines share, with two orders of headroom:
  !
  !   * the mapped position: a 27-term sum of shape functions times
  !     anchor coordinates, no cancellation (the shape functions sum to
  !     one), so the two routines can differ by at most ~27 u times the
  !     summand magnitude, ~1e-14 relative. Asserted at 1e-12.
  !
  !   * the inverse Jacobian: the x_xi-type sums cancel -- sum(dershape3D)
  !     is zero, so an anchor offset of order 1 against an element
  !     half-width of 1e-2 loses two to three digits -- and the cofactor
  !     over determinant division compounds it, ~1e-12 to 1e-11 relative.
  !     Asserted at 1e-9, relative to the largest entry of the reference
  !     matrix rather than entry by entry, so that a cofactor that
  !     happens to be near zero cannot inflate the ratio.
  !
  ! The bit-for-bit count is reported after them, for the record.
  !--------------------------------------------------------------------

  call gf_lcg_seed(987654321)

  nskipped = 0
  nmismatch = 0
  worst_xyz = 0.d0
  worst_jinv = 0.d0

  do itrial = 1,NTRIAL

    ! a random affine element, at a realistic position and size: centred
    ! about a unit-radius point with a half-width of order 1e-2, which is
    ! where the cancellation in sum(dershape3D) = 0 actually bites
    do i = 1,NDIM
      x0(i) = gf_rand_range(-1.d0,1.d0)
      do j = 1,NDIM
        amat(i,j) = gf_rand_range(-1.d-2,1.d-2)
      enddo
    enddo

    ! keeps the map orientation-preserving; a negative determinant is
    ! tested separately below, and must never reach the oracle
    det = amat(1,1)*(amat(2,2)*amat(3,3) - amat(2,3)*amat(3,2)) &
        - amat(1,2)*(amat(2,1)*amat(3,3) - amat(2,3)*amat(3,1)) &
        + amat(1,3)*(amat(2,1)*amat(3,2) - amat(2,2)*amat(3,1))
    if (det < 0.d0) then
      do i = 1,NDIM
        amat(i,1) = -amat(i,1)
      enddo
    endif

    do ia = 1,NGNOD
      xelm(ia) = x0(1) + amat(1,1)*xiref(1,ia) + amat(1,2)*xiref(2,ia) + amat(1,3)*xiref(3,ia)
      yelm(ia) = x0(2) + amat(2,1)*xiref(1,ia) + amat(2,2)*xiref(2,ia) + amat(2,3)*xiref(3,ia)
      zelm(ia) = x0(3) + amat(3,1)*xiref(1,ia) + amat(3,2)*xiref(2,ia) + amat(3,3)*xiref(3,ia)
    enddo

    ! sampled slightly outside [-1,1] as well: the solver accepts local
    ! coordinates up to 1.1 (locate_point.f90:513) and the shipped global
    ! example genuinely uses gamma = 1.053
    xi    = gf_rand_range(-1.1d0,1.1d0)
    eta   = gf_rand_range(-1.1d0,1.1d0)
    gamma = gf_rand_range(-1.1d0,1.1d0)

    ! the routine under test goes first, so that a degenerate draw is
    ! skipped rather than handed to an oracle that would stop on it
    call gf_shape3D_map(xelm,yelm,zelm,xi,eta,gamma,xyz,jinv,jacobian,ierr)
    if (ierr /= GF_OK) then
      nskipped = nskipped + 1
      cycle
    endif

    x = 0.d0 ; y = 0.d0 ; z = 0.d0
    xix = 0.d0 ; xiy = 0.d0 ; xiz = 0.d0
    etax = 0.d0 ; etay = 0.d0 ; etaz = 0.d0
    gammax = 0.d0 ; gammay = 0.d0 ; gammaz = 0.d0

    call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                            xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

    ! the reference packed the way gf_shape3D_map packs it -- rows
    ! xi/eta/gamma, columns x/y/z -- which is the one thing the fork does
    ! that the original does not, and which a transposition would fail at
    ! order one
    jref(1,1) = xix ; jref(1,2) = xiy ; jref(1,3) = xiz
    jref(2,1) = etax ; jref(2,2) = etay ; jref(2,3) = etaz
    jref(3,1) = gammax ; jref(3,2) = gammay ; jref(3,3) = gammaz

    ! the summand magnitude of the position sums
    xscale = max(maxval(abs(xelm)),maxval(abs(yelm)),maxval(abs(zelm)))

    worst_xyz = max(worst_xyz,max(abs(xyz(1) - x),abs(xyz(2) - y),abs(xyz(3) - z))/xscale)
    worst_jinv = max(worst_jinv,maxval(abs(jinv - jref))/maxval(abs(jref)))

    if (xyz(1) /= x .or. xyz(2) /= y .or. xyz(3) /= z .or. &
        jinv(1,1) /= xix .or. jinv(1,2) /= xiy .or. jinv(1,3) /= xiz .or. &
        jinv(2,1) /= etax .or. jinv(2,2) /= etay .or. jinv(2,3) /= etaz .or. &
        jinv(3,1) /= gammax .or. jinv(3,2) /= gammay .or. jinv(3,3) /= gammaz) then
      nmismatch = nmismatch + 1
    endif

  enddo

  write(*,'(a,i0,a,i0,a)') '  compared ',NTRIAL - nskipped,' of ',NTRIAL, &
                           ' trials (the rest were degenerate draws)'

  call gf_report('position vs recompute_jacobian    ',worst_xyz,1.d-12,nfail)
  call gf_report('inverse Jacobian vs recompute_jac.',worst_jinv,1.d-9,nfail)

  ! informational: zero under a value-safe floating-point model, non-zero
  ! under an FMA-contracting or vectorising build (ifort/ifx -O3 -xHost)
  write(*,'(a,i0,a,i0,a)') '  bit-for-bit mismatches = ',nmismatch,' of ',NTRIAL - nskipped, &
                           ' (informational: 0 under a value-safe FP model)'

  ! a run in which everything was skipped would report success while
  ! comparing nothing
  call gf_report_true('the comparison ran at all         ',NTRIAL - nskipped > NTRIAL/2,nfail)

  !--------------------------------------------------------------------
  ! 4. degenerate geometry returns ierr and does not stop
  !
  ! This is the assertion the fork exists for. recompute_jacobian is
  ! deliberately not called on either case: it would stop, and a Fortran
  ! stop exits with status zero, so the runner would see a *passing* test
  ! that had silently stopped executing.
  !
  ! Reaching the report line at all is half the assertion.
  !--------------------------------------------------------------------

  ! mirrored element: det A < 0, so jacobian < 0
  do ia = 1,NGNOD
    xelm(ia) = -1.d-2*xiref(1,ia)
    yelm(ia) =  1.d-2*xiref(2,ia)
    zelm(ia) =  1.d-2*xiref(3,ia)
  enddo
  call gf_shape3D_map(xelm,yelm,zelm,0.1d0,0.2d0,0.3d0,xyz,jinv,jacobian,ierr)
  call gf_report_true('mirrored element returns ierr     ',ierr == GF_ERR_GEOMETRY,nfail)
  call gf_report_true('  ... and zeroes jinv             ',maxval(abs(jinv)) == 0.d0,nfail)

  ! collapsed element: every anchor at the same point, so jacobian = 0
  do ia = 1,NGNOD
    xelm(ia) = 0.5d0
    yelm(ia) = 0.5d0
    zelm(ia) = 0.5d0
  enddo
  call gf_shape3D_map(xelm,yelm,zelm,0.1d0,0.2d0,0.3d0,xyz,jinv,jacobian,ierr)
  call gf_report_true('collapsed element returns ierr    ',ierr == GF_ERR_GEOMETRY,nfail)

  ! a flat element: three anchors' worth of extent in two directions only
  do ia = 1,NGNOD
    xelm(ia) = 1.d-2*xiref(1,ia)
    yelm(ia) = 1.d-2*xiref(2,ia)
    zelm(ia) = 0.d0
  enddo
  call gf_shape3D_map(xelm,yelm,zelm,0.1d0,0.2d0,0.3d0,xyz,jinv,jacobian,ierr)
  call gf_report_true('flat element returns ierr         ',ierr == GF_ERR_GEOMETRY,nfail)

  write(*,'(a)') ''
  write(*,'(a)') '  the process is still running, so no path taken above called stop'
  write(*,'(a)') ''

  !--------------------------------------------------------------------

  if (nfail /= 0) then
    write(*,'(a,i0,a)') 'test_gf_shape3D: ',nfail,' assertion(s) FAILED'
    ! `stop 1`, not a bare stop: a bare stop exits with status zero and the
    ! runner would call that a pass
    stop 1
  endif

  write(*,'(a)') 'test_gf_shape3D: all assertions passed'

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine anchor_reference_coords(xiref_out)

! reference (xi,eta,gamma) of the 27 anchors
!
! Taken from hex_nodes() rather than a hand-written table, so that a
! reordering of the 27 control nodes breaks this test instead of silently
! decoupling it from gf_shape3D_functions(). hex_nodes returns the topology
! as 0/1/2 per direction; the reference coordinate is that minus one.

  implicit none

  double precision, dimension(NDIM,NGNOD), intent(out) :: xiref_out

  ! local parameters
  integer, dimension(NGNOD) :: iaddx,iaddy,iaddz
  integer :: ia_l

  call hex_nodes(iaddx,iaddy,iaddz)

  do ia_l = 1,NGNOD
    xiref_out(1,ia_l) = dble(iaddx(ia_l)) - 1.d0
    xiref_out(2,ia_l) = dble(iaddy(ia_l)) - 1.d0
    xiref_out(3,ia_l) = dble(iaddz(ia_l)) - 1.d0
  enddo

  end subroutine anchor_reference_coords

  end program test_gf_shape3D
