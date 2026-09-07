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
!---- test_gf_interp -- src/gf3d/gf_interp.F90
!----
!---- NGLLX = 5 gives a degree-4 Lagrange basis, which reproduces
!---- polynomials of per-variable degree <= 4 exactly. So a manufactured
!---- polynomial field has an expected value that can be written down in
!---- closed form, needs no database, no HDF5 and no MPI, and runs in CI on
!---- every commit.
!----
!---- The oracle here is gf_poly3() and direct powers -- a plain sum over
!---- monomials. It never evaluates a Lagrange basis, so it shares no code
!---- with the routine under test. That property is what makes this a
!---- legitimate replacement for the Python comparison that was abolished,
!---- rather than a repeat of the same mistake.
!----
!---- The three factors, and what each one catches
!---- -------------------------------------------
!---- The manufactured field is
!----
!----     u(f,d,i,j,k,t) = s(f,d) * q(t) * P_(f,d)(xi_i, eta_j, gamma_k)
!----
!---- and every factor is there to fail a specific transposition:
!----
!----   q(t)              a wrong time stride
!----   s(f,d) = f + 10d  the force and displacement axes swapped -- which is
!----                     invisible if the field is symmetric in them, and
!----                     both extents are 3
!----   P_(f,d) *not* symmetric under permuting its three arguments -- the
!----                     (i,j,k) <-> (xi,eta,gamma) mix-up
!----
!---- That last one is the likeliest bug in the whole port: GF3DF collapsed
!---- the GLL indices into a single iglob, so every loop reindexed from it is
!---- an opportunity to transpose.
!----
!---- Tolerances are derived, not fitted, and are **relative** wherever the
!---- manufactured value carries a magnitude -- the s(f,d) factor reaches 33,
!---- so an absolute bound would be testing the size of the scale rather than
!---- the accuracy of the interpolation. 125 fused multiply-adds bounded by
!---- 125*u*Lambda^3 is about 7e-14 relative, so 1e-12 leaves two orders of
!---- headroom against -O3 reassociation. The CUSTOM_REAL wrapper is at 5e-7
!---- and **the measured value is printed**: it is the float32 noise floor
!---- that Stage 8's finite-difference budget depends on.
!----

  program test_gf_interp

  use gf_par, only: GF_NCOMP

  use gf_interp, only: gf_interp_weights,gf_interp_snapshot, &
                       gf_interp_trace_d,gf_interp_trace

  use gf_manufactured, only: gf_lcg_seed,gf_rand_range,gf_report,gf_report_true,gf_poly3

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,GAUSSALPHA,GAUSSBETA

  implicit none

  ! a short time axis: the point of the time index here is that it exists
  ! and varies, not that it is long
  integer, parameter :: NT = 4

  ! local parameters
  integer :: nfail,i,j,k,f,d,it,ia,ib,ic
  double precision, dimension(NGLLX) :: xigll,wxgll,hxi
  double precision, dimension(NGLLY) :: yigll,wygll,heta
  double precision, dimension(NGLLZ) :: zigll,wzgll,hgam
  double precision, dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ,NT) :: u
  real(kind=CUSTOM_REAL), dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ,NT) :: ur
  double precision, dimension(GF_NCOMP,GF_NCOMP,NT) :: out,out_r
  double precision, dimension(GF_NCOMP,GF_NCOMP) :: snap
  double precision, dimension(GF_NCOMP,GF_NCOMP) :: s
  double precision, dimension(NT) :: q
  double precision, dimension(0:4,0:4,0:4) :: coef
  double precision :: xi,eta,gam,expected,err,worst,worst_rel,denom
  integer :: aworst,bworst,cworst

  nfail = 0

  write(*,'(a)') ''
  write(*,'(a)') 'test_gf_interp'
  write(*,'(a)') ''
  write(*,'(a,i0,a,i0,a,i0)') '  NGLLX,NGLLY,NGLLZ = ',NGLLX,',',NGLLY,',',NGLLZ
  write(*,'(a,i0)') '  CUSTOM_REAL = ',CUSTOM_REAL
  write(*,'(a)') ''

  ! GLL abscissae from the solver's own quadrature, never hard-coded: a
  ! change of quadrature must desynchronise this test loudly rather than
  ! leave it passing while testing a different basis than production uses
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  ! the per-pair scale: distinct for every (force,displacement) pair, and
  ! deliberately not symmetric under swapping them
  do d = 1,GF_NCOMP
    do f = 1,GF_NCOMP
      s(f,d) = dble(f) + 10.d0*dble(d)
    enddo
  enddo

  ! the per-sample time factor
  do it = 1,NT
    q(it) = 1.d0 + 0.5d0*dble(it)
  enddo

  !--------------------------------------------------------------------
  ! 1. partition of unity
  !
  ! P = 1 and q = 1, so the answer is s(f,d) itself. This is the simplest
  ! possible check that the (f,d) pair is carried through untransposed.
  !--------------------------------------------------------------------

  do it = 1,NT
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          do d = 1,GF_NCOMP
            do f = 1,GF_NCOMP
              u(f,d,i,j,k,it) = s(f,d)
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  xi = 0.3123d0 ; eta = -0.4571d0 ; gam = 0.8219d0
  call gf_interp_weights(xi,eta,gam,hxi,heta,hgam)
  call gf_interp_trace_d(u,hxi,heta,hgam,NT,out)

  worst = 0.d0
  do it = 1,NT
    do d = 1,GF_NCOMP
      do f = 1,GF_NCOMP
        worst = max(worst,abs(out(f,d,it) - s(f,d))/abs(s(f,d)))
      enddo
    enddo
  enddo
  call gf_report('partition of unity                ',worst,1.d-14,nfail)

  !--------------------------------------------------------------------
  ! 2. Kronecker property at all 125 nodes
  !
  ! A field that is 1 at one node and 0 elsewhere must interpolate to
  ! exactly 1 at that node's coordinates. Looping all 125 pins the
  ! association between the storage index (i,j,k) and the coordinate
  ! (xi,eta,gamma) node by node -- a transposition shows up as a 1 landing
  ! at the wrong node.
  !--------------------------------------------------------------------

  worst = 0.d0
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        u(:,:,:,:,:,1) = 0.d0
        u(:,:,i,j,k,1) = 1.d0

        call gf_interp_weights(xigll(i),yigll(j),zigll(k),hxi,heta,hgam)
        call gf_interp_snapshot(u(:,:,:,:,:,1),hxi,heta,hgam,snap)

        do d = 1,GF_NCOMP
          do f = 1,GF_NCOMP
            worst = max(worst,abs(snap(f,d) - 1.d0))
          enddo
        enddo

      enddo
    enddo
  enddo
  call gf_report('Kronecker at all 125 GLL nodes    ',worst,1.d-14,nfail)

  !--------------------------------------------------------------------
  ! 3. monomial sweep
  !
  ! P = xi^a eta^b gamma^c for every (a,b,c) with a,b,c <= 4 -- the exact
  ! span of the tensor-product basis. Sweeping rather than testing one
  ! polynomial localises a failure to a direction and a degree, which is
  ! what turns "the interpolator is wrong" into "the gamma axis is wrong".
  !--------------------------------------------------------------------

  xi = 0.3123d0 ; eta = -0.4571d0 ; gam = 0.8219d0
  call gf_interp_weights(xi,eta,gam,hxi,heta,hgam)

  worst = 0.d0
  aworst = -1 ; bworst = -1 ; cworst = -1

  do ic = 0,4
    do ib = 0,4
      do ia = 0,4

        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              u(:,:,i,j,k,1) = xigll(i)**ia * yigll(j)**ib * zigll(k)**ic
            enddo
          enddo
        enddo

        call gf_interp_snapshot(u(:,:,:,:,:,1),hxi,heta,hgam,snap)

        expected = xi**ia * eta**ib * gam**ic
        err = abs(snap(1,1) - expected)
        if (err > worst) then
          worst = err
          aworst = ia ; bworst = ib ; cworst = ic
        endif

      enddo
    enddo
  enddo

  call gf_report('monomial sweep, all (a,b,c) <= 4  ',worst,1.d-12,nfail)
  if (aworst >= 0) then
    write(*,'(a,i0,a,i0,a,i0,a)') '     worst monomial: xi^',aworst,' eta^',bworst,' gamma^',cworst
  endif

  !--------------------------------------------------------------------
  ! 4. a dense random Q4 field, with all three factors at once
  !
  ! The sweep tests one monomial at a time and so cannot see cross-term
  ! mixing; a dense field can. The coefficients are drawn from the
  ! deterministic LCG, so a CI failure reproduces locally.
  !--------------------------------------------------------------------

  call gf_lcg_seed(31415926)

  do ic = 0,4
    do ib = 0,4
      do ia = 0,4
        coef(ia,ib,ic) = gf_rand_range(-1.d0,1.d0)
      enddo
    enddo
  enddo

  ! break the symmetry between the three directions explicitly: without
  ! this a coefficient array that happened to be near-symmetric could let a
  ! transposition through
  coef(4,0,0) = 3.0d0
  coef(0,4,0) = -2.0d0
  coef(0,0,4) = 1.5d0
  coef(3,1,0) = 2.5d0
  coef(0,1,3) = -1.75d0

  call build_field(u,coef,s,q,xigll,yigll,zigll)

  xi = -0.6417d0 ; eta = 0.2288d0 ; gam = -0.9013d0
  call gf_interp_weights(xi,eta,gam,hxi,heta,hgam)
  call gf_interp_trace_d(u,hxi,heta,hgam,NT,out)

  expected = gf_poly3(coef,4,4,4,xi,eta,gam)

  worst = 0.d0
  do it = 1,NT
    do d = 1,GF_NCOMP
      do f = 1,GF_NCOMP
        worst = max(worst,abs(out(f,d,it) - s(f,d)*q(it)*expected) &
                          / abs(s(f,d)*q(it)*expected))
      enddo
    enddo
  enddo
  call gf_report('dense random Q4, all three factors',worst,1.d-12,nfail)

  !--------------------------------------------------------------------
  ! 5. extrapolation beyond the element
  !
  ! The solver accepts local coordinates outside [-1,1] -- the shipped
  ! global example genuinely locates its source at gamma = 1.053, and the
  ! containment rule in gf_par deliberately follows it there. This is the
  ! only place in the suite where the basis is checked in that regime, so
  ! the point used is that example's own solver-reported location.
  !--------------------------------------------------------------------

  xi = 0.755875349d0 ; eta = 0.525376797d0 ; gam = 1.05304706d0
  call gf_interp_weights(xi,eta,gam,hxi,heta,hgam)
  call gf_interp_trace_d(u,hxi,heta,hgam,NT,out)

  expected = gf_poly3(coef,4,4,4,xi,eta,gam)

  worst = 0.d0
  do it = 1,NT
    do d = 1,GF_NCOMP
      do f = 1,GF_NCOMP
        worst = max(worst,abs(out(f,d,it) - s(f,d)*q(it)*expected) &
                          / abs(s(f,d)*q(it)*expected))
      enddo
    enddo
  enddo
  call gf_report('extrapolated to gamma = 1.053     ',worst,1.d-11,nfail)

  !--------------------------------------------------------------------
  ! 6. negative control
  !
  ! xi^5 is one degree above the span, so the interpolant must *not*
  ! reproduce it. Without this, an implementation that somehow returned the
  ! right answer by construction -- or a test that had quietly stopped
  ! exercising the basis at all -- would pass everything above.
  !--------------------------------------------------------------------

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        u(:,:,i,j,k,1) = xigll(i)**5
      enddo
    enddo
  enddo

  xi = 0.3123d0 ; eta = -0.4571d0 ; gam = 0.8219d0
  call gf_interp_weights(xi,eta,gam,hxi,heta,hgam)
  call gf_interp_snapshot(u(:,:,:,:,:,1),hxi,heta,hgam,snap)

  err = abs(snap(1,1) - xi**5)
  call gf_report_true('xi^5 is NOT reproduced (control)  ',err > 1.d-4,nfail)
  write(*,'(a,es12.5)') '     deviation from xi^5 = ',err

  !--------------------------------------------------------------------
  ! 7. the CUSTOM_REAL wrapper against the double core
  !
  ! Same field, same point, the only difference being that the wrapper's
  ! input has been through float32. The measured value is the noise floor
  ! this library imposes on anything driven from the database, and Stage 8's
  ! finite-difference step-size budget depends on knowing it -- so it is
  ! printed, not merely asserted.
  !--------------------------------------------------------------------

  call build_field(u,coef,s,q,xigll,yigll,zigll)
  ur(:,:,:,:,:,:) = real(u(:,:,:,:,:,:),kind=CUSTOM_REAL)

  xi = -0.6417d0 ; eta = 0.2288d0 ; gam = -0.9013d0
  call gf_interp_weights(xi,eta,gam,hxi,heta,hgam)
  call gf_interp_trace_d(u,hxi,heta,hgam,NT,out)
  call gf_interp_trace(ur,hxi,heta,hgam,NT,out_r)

  worst_rel = 0.d0
  do it = 1,NT
    do d = 1,GF_NCOMP
      do f = 1,GF_NCOMP
        denom = abs(out(f,d,it))
        if (denom < 1.d-30) cycle
        worst_rel = max(worst_rel,abs(out_r(f,d,it) - out(f,d,it))/denom)
      enddo
    enddo
  enddo
  call gf_report('CUSTOM_REAL wrapper vs double core',worst_rel,5.d-7,nfail)
  write(*,'(a,es12.5)') '     measured float32 noise floor (relative) = ',worst_rel

  ! and with a double-precision input the wrapper must agree with the core
  ! exactly, since widening is exact -- this is what makes the seam a seam
  ! rather than a second implementation
  if (CUSTOM_REAL == 8) then
    worst = 0.d0
    do it = 1,NT
      do d = 1,GF_NCOMP
        do f = 1,GF_NCOMP
          worst = max(worst,abs(out_r(f,d,it) - out(f,d,it)))
        enddo
      enddo
    enddo
    call gf_report('wrapper == core at CUSTOM_REAL = 8',worst,0.d0,nfail)
  endif

  !--------------------------------------------------------------------
  ! 8. the snapshot and trace entry points must agree
  !--------------------------------------------------------------------

  call gf_interp_snapshot(u(:,:,:,:,:,2),hxi,heta,hgam,snap)
  worst = 0.d0
  do d = 1,GF_NCOMP
    do f = 1,GF_NCOMP
      worst = max(worst,abs(snap(f,d) - out(f,d,2)))
    enddo
  enddo
  call gf_report('snapshot == trace, same time slice',worst,0.d0,nfail)

  !--------------------------------------------------------------------

  write(*,'(a)') ''
  if (nfail /= 0) then
    write(*,'(a,i0,a)') 'test_gf_interp: ',nfail,' assertion(s) FAILED'
    stop 1
  endif

  write(*,'(a)') 'test_gf_interp: all assertions passed'

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine build_field(uu,cf,ss,qq,xg,yg,zg)

! fills the manufactured field in the database's own storage layout
!
! u(f,d,i,j,k,t) = s(f,d) * q(t) * P(xi_i, eta_j, gamma_k)
!
! P comes from gf_poly3 -- a plain monomial sum -- so the field and the
! expected value are produced by the same closed-form route, and neither
! touches a Lagrange basis.

  implicit none

  double precision, dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ,NT), intent(out) :: uu
  double precision, dimension(0:4,0:4,0:4), intent(in) :: cf
  double precision, dimension(GF_NCOMP,GF_NCOMP), intent(in) :: ss
  double precision, dimension(NT), intent(in) :: qq
  double precision, dimension(NGLLX), intent(in) :: xg
  double precision, dimension(NGLLY), intent(in) :: yg
  double precision, dimension(NGLLZ), intent(in) :: zg

  ! local parameters
  integer :: il,jl,kl,fl,dl,itl
  double precision :: p

  do kl = 1,NGLLZ
    do jl = 1,NGLLY
      do il = 1,NGLLX
        p = gf_poly3(cf,4,4,4,xg(il),yg(jl),zg(kl))
        do itl = 1,NT
          do dl = 1,GF_NCOMP
            do fl = 1,GF_NCOMP
              uu(fl,dl,il,jl,kl,itl) = ss(fl,dl)*qq(itl)*p
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  end subroutine build_field

  end program test_gf_interp
