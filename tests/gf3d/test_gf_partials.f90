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
!---- test_gf_partials -- src/gf3d/gf_partials.F90 and the derivative kernel
!---- in src/gf3d/gf_stf.F90
!----
!---- Tier 1: no database, no HDF5, no MPI. The strain traces are
!---- manufactured, the plan and the axis are the shipped regional example's
!---- (as in test_gf_stf), and the oracles are identities:
!----
!----   * linearity -- SUM_v M_v dp(v) reproduces the seismogram assembled
!----     the production way, and each dp(v) *is* the seismogram of the v-th
!----     unit tensor, bitwise, because it runs the same statements;
!----   * the derivative kernel sums to one for every width, with the sum
!----     before normalisation equal to the Poisson-summation closed form
!----     1 + 2 SUM_k exp(-(pi k h/dt)^2), which is what pins the claim that
!----     the normalisation is a no-op when the kernel is resolved;
!----   * the centroid-time partial returns the kernel for a delta, and
!----     against a central difference of the converted trace it converges
!----     at second order under step halving -- the finite difference is the
!----     one with the O(dt^2) error, and the ratio says so.
!----
!---- Stage 8 adds its sections to this program.
!----

  program test_gf_partials

  use constants, only: PI

  use gf_par, only: t_gf_stf,t_gf_taxis,GF_OK,GF_NCOMP, &
                    GF_SRC_CMT,GF_SRC_FORCE,GF_STF_TRUNC,GF_STF_HEAVI,GF_STF_GAUSS

  use gf_strain, only: GF_VOIGT

  use gf_moment, only: gf_rotate_moment_tensor,gf_moment_contract

  use gf_stf, only: gf_stf_plan,gf_taxis_plan,gf_stf_kernel,gf_stf_kernel_gauss, &
                    gf_stf_kernel_gauss_unit,gf_stf_khalf,gf_pad_left,gf_cumsum,gf_stf_apply

  use gf_partials

  use gf_manufactured

  implicit none

  ! the shipped regional example's widths and grid
  double precision, parameter :: HDB  = 6.94968291528492d0        ! station attribute hdur
  double precision, parameter :: HCMT = 60.d0                      ! CMTSOLUTION half duration
  double precision, parameter :: T0DB = 34.74841457642461d0        ! mesh_info.h5 t0
  double precision, parameter :: DTR  = 3.4d0                      ! stored spacing
  integer, parameter :: NTDB = 544, SS = 34

  integer :: nfail

  nfail = 0

  write(*,'(a)') 'test_gf_partials: moment-tensor and centroid-time partials'
  write(*,'(a)') ''

  call test_ndp(nfail)
  call test_unit_kernel(nfail)
  call test_mt_linearity(nfail)
  call test_time_partial(nfail)

  write(*,'(a)') ''
  if (nfail > 0) then
    write(*,'(a,i0,a)') 'test_gf_partials: ',nfail,' assertion(s) FAILED'
    write(0,'(a,i0,a)') 'test_gf_partials: ',nfail,' assertion(s) FAILED'
    stop 1
  endif
  write(*,'(a)') 'test_gf_partials: all assertions passed'

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_ndp(nfail)

! the kernel-type to count mapping

  implicit none
  integer, intent(inout) :: nfail

  integer :: ndp,ierr

  write(*,'(a)') '1. itypsokern'

  call gf_partials_ndp(0,ndp,ierr)
  call gf_report_true('itypsokern 0: no partials',ierr == GF_OK .and. ndp == 0,nfail)
  call gf_partials_ndp(1,ndp,ierr)
  call gf_report_true('itypsokern 1: six',ierr == GF_OK .and. ndp == GF_NDP_MT .and. ndp == 6,nfail)
  call gf_partials_ndp(2,ndp,ierr)
  call gf_report_true('itypsokern 2: ten',ierr == GF_OK .and. ndp == GF_NDP_LOC .and. ndp == 10,nfail)
  call gf_partials_ndp(3,ndp,ierr)
  call gf_report_true('itypsokern 3 (half duration) is refused',ierr /= GF_OK .and. ndp == 0,nfail)
  call gf_report_true('slot names in GF3DF order', &
                      GF_DP_NAME(GF_DP_MRR) == 'Mrr' .and. GF_DP_NAME(GF_DP_MTP) == 'Mtp' .and. &
                      GF_DP_NAME(GF_DP_LAT) == 'lat' .and. GF_DP_NAME(GF_DP_DEP) == 'dep' .and. &
                      GF_DP_NAME(GF_DP_TIM) == 'tim',nfail)

  end subroutine test_ndp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_unit_kernel(nfail)

! the normalised sampled Gaussian: unit sum, symmetry, the Poisson closed
! form for the raw sum, the identity at zero width, and bitwise identity
! with the sampled Gaussian once the kernel is resolved

  implicit none
  integer, intent(inout) :: nfail

  integer, parameter :: NR = 8
  double precision, dimension(NR), parameter :: ratios = &
    (/ 0.d0, 0.1d0, 0.25d0, 0.5d0, 1.d0, 1.5d0, 2.d0, 10.d0 /)
  double precision, dimension(:), allocatable :: w,wg
  double precision :: h,wsum,s,closed,worst_sum,worst_poisson,worst_resolved
  integer :: ir,k,j,khalf,nbad_sym,nbad_id,nbad_res
  character(len=64) :: label

  write(*,'(a)') '2. the derivative kernel'

  worst_sum = 0.d0
  worst_poisson = 0.d0
  worst_resolved = 0.d0
  nbad_sym = 0
  nbad_id = 0
  nbad_res = 0

  do ir = 1,NR
    h = ratios(ir)*DTR
    khalf = gf_stf_khalf(h,DTR,GF_STF_TRUNC)
    allocate(w(-khalf:khalf),wg(-khalf:khalf))

    call gf_stf_kernel_gauss_unit(h,DTR,khalf,w,wsum)

    ! unit sum: positive terms, summed from the tails inwards
    s = 0.d0
    do j = khalf,1,-1
      s = s + w(j) + w(-j)
    enddo
    s = s + w(0)
    worst_sum = max(worst_sum,abs(s - 1.d0))

    ! symmetry, bitwise
    do j = 1,khalf
      if (w(-j) /= w(j)) nbad_sym = nbad_sym + 1
    enddo

    if (ratios(ir) == 0.d0) then
      ! the identity: derivative of the trapezoid limit
      if (khalf /= 0 .or. w(0) /= 1.d0 .or. wsum /= 1.d0) nbad_id = nbad_id + 1
    else
      ! Poisson summation for the raw sum of the sampled Gaussian: the
      ! aliasing is a theta-function tail, and it is what the normalisation
      ! removes
      closed = 1.d0
      do k = 1,60
        closed = closed + 2.d0*exp(-(PI*dble(k)*ratios(ir))**2)
      enddo
      worst_poisson = max(worst_poisson,abs(wsum - closed))
      write(label,'(a,f5.2,a,es10.3)') '     h/dt = ',ratios(ir),'  raw sum - 1 = ',wsum - 1.d0
      write(*,'(a)') trim(label)
    endif

    ! resolved: the normalisation is below one ulp of 1 and changes nothing
    if (ratios(ir) >= 2.d0) then
      call gf_stf_kernel_gauss(h,DTR,khalf,wg)
      do j = -khalf,khalf
        if (w(j) /= wg(j)) nbad_res = nbad_res + 1
      enddo
    endif

    ! barely resolved: within the Poisson bound of the sampled Gaussian
    if (ratios(ir) == 1.5d0) then
      call gf_stf_kernel_gauss(h,DTR,khalf,wg)
      do j = -khalf,khalf
        worst_resolved = max(worst_resolved,abs(w(j) - wg(j))/wg(0))
      enddo
    endif

    deallocate(w,wg)
  enddo

  call gf_report('unit sum for h/dt in {0 .. 10}      ',worst_sum,1.d-13,nfail)
  call gf_report_true('symmetric, bitwise                  ',nbad_sym == 0,nfail)
  call gf_report_true('h = 0 is the identity, w(0) = 1     ',nbad_id == 0,nfail)
  call gf_report('raw sum vs Poisson closed form      ',worst_poisson,1.d-12,nfail)
  call gf_report_true('h >= 2 dt: equals the sampled Gaussian, bitwise',nbad_res == 0,nfail)
  call gf_report('h = 1.5 dt: within 2 exp(-(1.5 pi)^2) = 4.6e-10 of it', &
                 worst_resolved,1.d-9,nfail)

  end subroutine test_unit_kernel

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_mt_linearity(nfail)

! the six moment-tensor partials against the seismogram assembled the
! production way, on manufactured strain traces

  implicit none
  integer, intent(inout) :: nfail

  double precision, dimension(:,:,:), allocatable :: eps,dp
  double precision, dimension(:,:), allocatable :: seis,recon,seis_v
  double precision, dimension(:), allocatable :: trace,xpad,p,y,w
  double precision, dimension(6) :: m_sph,e_sph
  double precision, dimension(3,3) :: m_cart
  type(t_gf_stf) :: stf,stf_g
  type(t_gf_taxis) :: tax
  double precision :: theta,phi,scale_amp,scale_moment,scale_mt,t,worst,ref
  integer :: v,icomp,it,ierr,nt,nbad

  write(*,'(a)') '3. moment-tensor partials: linearity'

  ! the regional CMT plan and axis
  call gf_stf_plan(GF_SRC_CMT,0,HCMT,HDB,DTR,GF_STF_TRUNC,stf,ierr)
  call gf_report_true('regional CMT plan: Heaviside, khalf = 64', &
                      ierr == GF_OK .and. stf%kind_stf == GF_STF_HEAVI .and. stf%khalf == 64,nfail)
  call gf_taxis_plan(NTDB,0.1d0,SS,T0DB,90.d0,tax,ierr)
  call gf_report_true('regional CMT axis: npad = 18, nt = 562', &
                      ierr == GF_OK .and. tax%npad == 18 .and. tax%nt == 562,nfail)
  nt = tax%nt

  allocate(eps(GF_VOIGT,GF_NCOMP,NTDB),dp(GF_NDP_MT,GF_NCOMP,nt), &
           seis(GF_NCOMP,nt),recon(GF_NCOMP,nt),seis_v(GF_NCOMP,nt), &
           trace(NTDB),xpad(nt),p(0:nt),y(nt),w(-stf%khalf:stf%khalf))

  ! a smooth, distinct trace per Voigt slot and component: a wave packet
  ! whose period and centre depend on the slot, so that no two partials
  ! are proportional and a slot mix-up cannot cancel
  call gf_lcg_seed(606)
  do it = 1,NTDB
    t = (dble(it*SS) - 1.d0)*0.1d0 - T0DB
    do icomp = 1,GF_NCOMP
      do v = 1,GF_VOIGT
        eps(v,icomp,it) = exp(-((t - 600.d0 - 40.d0*v)/150.d0)**2) &
                          * cos(2.d0*PI*t/(45.d0 + 7.d0*v + 3.d0*icomp)) &
                          * (1.d0 + 0.1d0*gf_rand_range(-1.d0,1.d0))
      enddo
    enddo
  enddo

  ! a moment tensor in "non-dimensional" units, its scale, and the station
  ! amplitude factor: numbers of the shipped example's order, none of them
  ! round
  do v = 1,6
    m_sph(v) = gf_rand_range(-1.d0,1.d0)
  enddo
  scale_moment = 2.6299730036637251d28 / 0.73d0
  scale_amp = 1.d0 / 1.0537d15
  scale_mt = scale_amp / scale_moment
  theta = 1.672d0
  phi = 4.971d0

  call gf_stf_kernel(stf,tax%dt_sub,w)

  ! the seismogram, the production way (gf_seis_cmt's statements)
  call gf_rotate_moment_tensor(theta,phi,m_sph,m_cart)
  do icomp = 1,GF_NCOMP
    do it = 1,NTDB
      call gf_moment_contract(m_cart,eps(:,icomp,it),trace(it))
      trace(it) = scale_amp * trace(it)
    enddo
    call gf_pad_left(trace,NTDB,tax%npad,xpad)
    call gf_cumsum(xpad,nt,p)
    call gf_stf_apply(stf,tax%dt_sub,w,p,xpad,nt,y)
    seis(icomp,:) = y(:)
  enddo

  ! the partials
  call gf_partials_mt(eps,NTDB,theta,phi,scale_mt,tax,stf,w,dp,ierr)
  call gf_report_true('gf_partials_mt returns GF_OK',ierr == GF_OK,nfail)

  ! linearity: the CMTSOLUTION's own numbers times the partials
  recon(:,:) = 0.d0
  do it = 1,nt
    do icomp = 1,GF_NCOMP
      do v = 1,6
        recon(icomp,it) = recon(icomp,it) + (m_sph(v)*scale_moment)*dp(v,icomp,it)
      enddo
    enddo
  enddo
  ref = maxval(abs(seis))
  worst = maxval(abs(recon - seis))/ref
  call gf_report('SUM_v M_v dp(v) == seismogram (rel)  ',worst,1.d-12,nfail)

  ! each partial is the seismogram of its unit tensor, bitwise: the same
  ! statements on the same numbers
  nbad = 0
  do v = 1,6
    e_sph(:) = 0.d0
    e_sph(v) = 1.d0
    call gf_rotate_moment_tensor(theta,phi,e_sph,m_cart)
    do icomp = 1,GF_NCOMP
      do it = 1,NTDB
        call gf_moment_contract(m_cart,eps(:,icomp,it),trace(it))
        trace(it) = scale_mt * trace(it)
      enddo
      call gf_pad_left(trace,NTDB,tax%npad,xpad)
      call gf_cumsum(xpad,nt,p)
      call gf_stf_apply(stf,tax%dt_sub,w,p,xpad,nt,y)
      do it = 1,nt
        if (y(it) /= dp(v,icomp,it)) nbad = nbad + 1
      enddo
    enddo
  enddo
  call gf_report_true('dp(v) == seismogram of unit tensor v, bitwise',nbad == 0,nfail)

  ! the partials are not degenerate: no two proportional, none zero
  worst = 1.d0
  do v = 1,6
    worst = min(worst,maxval(abs(dp(v,:,:)))/maxval(abs(dp)))
  enddo
  call gf_report_true('every partial carries signal (min/max > 1e-3)',worst > 1.d-3,nfail)

  ! a Gaussian plan (a force source) is refused
  call gf_stf_plan(GF_SRC_FORCE,0,45.d0,HDB,DTR,GF_STF_TRUNC,stf_g,ierr)
  call gf_partials_mt(eps,NTDB,theta,phi,scale_mt,tax,stf_g,w,dp,ierr)
  call gf_report_true('a Gaussian plan is refused          ',ierr /= GF_OK,nfail)

  deallocate(eps,dp,seis,recon,seis_v,trace,xpad,p,y,w)

  end subroutine test_mt_linearity

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_time_partial(nfail)

! the centroid-time partial: the kernel for a delta, second-order agreement
! with a central difference of the converted trace, the guard limit

  implicit none
  integer, intent(inout) :: nfail

  double precision, dimension(:), allocatable :: x,p,y,w,wg,dp10
  type(t_gf_stf) :: stf
  double precision :: wsum,dt,t,err_coarse,err_fine,ref,ratio
  integer :: n,i,j,k0,khalf,ierr,nbad,igrid

  write(*,'(a)') '4. centroid-time partial'

  !--- a delta returns minus the kernel, bitwise ---------------------------

  call gf_stf_plan(GF_SRC_CMT,0,HCMT,HDB,DTR,GF_STF_TRUNC,stf,ierr)
  n = 400
  k0 = 200
  allocate(x(n),dp10(n),wg(-stf%khalf:stf%khalf))
  x(:) = 0.d0
  x(k0) = 1.d0
  call gf_partials_time(x,n,DTR,stf,dp10,wsum,ierr)
  call gf_report_true('gf_partials_time returns GF_OK      ',ierr == GF_OK,nfail)
  call gf_stf_kernel_gauss_unit(stf%hdur_corr,DTR,stf%khalf,wg,ref)
  nbad = 0
  do j = -stf%khalf,stf%khalf
    if (dp10(k0+j) /= -wg(j)) nbad = nbad + 1
  enddo
  do i = 1,n
    if (abs(i-k0) > stf%khalf .and. dp10(i) /= 0.d0) nbad = nbad + 1
  enddo
  call gf_report_true('delta in, minus the kernel out, bitwise',nbad == 0,nfail)
  call gf_report('  raw kernel sum - 1 (aliasing, h = 10.6 dt)',abs(wsum - 1.d0),1.d-15,nfail)
  deallocate(x,dp10,wg)

  !--- against a central difference of the converted trace ----------------
  !
  ! A Gaussian pulse of width 60 s on [0, 1800] s, converted with the
  ! regional CMT plan at dt and at dt/2. The analytic partial has no
  ! O(dt^2) term (the resolved kernel is spectrally exact); the central
  ! difference has (dt^2/6) y''', so their disagreement must fall by four
  ! when the step halves. The last khalf + 1 samples use zero-extended data
  ! and are left out, as is the first sample.

  do igrid = 1,2
    dt = DTR / dble(igrid)
    n = nint(1800.d0/dt)
    call gf_stf_plan(GF_SRC_CMT,0,HCMT,HDB,dt,GF_STF_TRUNC,stf,ierr)
    khalf = stf%khalf
    allocate(x(n),p(0:n),y(n),w(-khalf:khalf),dp10(n))
    do i = 1,n
      t = dble(i-1)*dt
      x(i) = exp(-((t - 800.d0)/60.d0)**2)
    enddo
    call gf_stf_kernel(stf,dt,w)
    call gf_cumsum(x,n,p)
    call gf_stf_apply(stf,dt,w,p,x,n,y)
    call gf_partials_time(x,n,dt,stf,dp10,wsum,ierr)

    ref = maxval(abs(dp10))
    err_fine = 0.d0
    do i = 2,n-khalf-1
      err_fine = max(err_fine,abs((y(i+1) - y(i-1))/(2.d0*dt) + dp10(i)))
    enddo
    err_fine = err_fine/ref
    if (igrid == 1) err_coarse = err_fine
    deallocate(x,p,y,w,dp10)
  enddo

  ratio = err_coarse/err_fine
  write(*,'(a,es10.3,a,es10.3,a,f6.3)') '     FD vs analytic: dt ',err_coarse,'  dt/2 ',err_fine, &
                                        '  ratio ',ratio
  call gf_report('central difference vs analytic at dt (rel)',err_coarse,1.d-2,nfail)
  call gf_report_true('  the difference is the FD''s: ratio in [3.5, 4.5]', &
                      ratio > 3.5d0 .and. ratio < 4.5d0,nfail)

  !--- the guard limit: minus the trace itself, bitwise -------------------

  call gf_stf_plan(GF_SRC_CMT,0,10.d0,HDB,DTR,GF_STF_TRUNC,stf,ierr)
  call gf_report_true('hdur = 10 s: the guard is set, khalf = 0',stf%guard .and. stf%khalf == 0,nfail)
  n = 300
  allocate(x(n),dp10(n))
  call gf_lcg_seed(1010)
  do i = 1,n
    x(i) = gf_rand_range(-1.d0,1.d0)
  enddo
  call gf_partials_time(x,n,DTR,stf,dp10,wsum,ierr)
  nbad = 0
  do i = 1,n
    if (dp10(i) /= -x(i)) nbad = nbad + 1
  enddo
  call gf_report_true('guard: dp(10) == -x, bitwise        ',ierr == GF_OK .and. nbad == 0,nfail)
  deallocate(x,dp10)

  !--- a Gaussian plan is refused ------------------------------------------

  call gf_stf_plan(GF_SRC_FORCE,0,45.d0,HDB,DTR,GF_STF_TRUNC,stf,ierr)
  allocate(x(10),dp10(10))
  x(:) = 1.d0
  call gf_partials_time(x,10,DTR,stf,dp10,wsum,ierr)
  call gf_report_true('a Gaussian plan is refused          ',ierr /= GF_OK,nfail)
  deallocate(x,dp10)

  end subroutine test_time_partial

  end program test_gf_partials
