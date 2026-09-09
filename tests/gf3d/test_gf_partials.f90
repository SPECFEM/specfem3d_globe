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
!----     unit tensor, to round-off, because it runs the same routines;
!----   * the derivative kernel sums to one for every width, with the sum
!----     before normalisation equal to the Poisson-summation closed form
!----     1 + 2 SUM_k exp(-(pi k h/dt)^2), which is what pins the claim that
!----     the normalisation is a no-op when the kernel is resolved;
!----   * the centroid-time partial returns the kernel for a delta, and
!----     against a central difference of the converted trace it converges
!----     at second order under step halving -- the finite difference is the
!----     one with the O(dt^2) error, and the ratio says so.
!----
!---- Stage 8 (sections 5 to 10): the analytic centroid partials, ingredient
!---- by ingredient and then as a whole. The rule of this suite is that the
!---- expected side shares no code with the routine under test; here that
!---- means closed forms where they exist (the second derivatives of a
!---- polynomial field on an affine element, the second derivatives of the
!---- tri-quadratic map, which any central difference reproduces exactly)
!---- and, for the geographic and rotational chains, Richardson-extrapolated
!---- central differences of the *forward* routines -- finite differences
!---- validating the analytic derivative, never the other way round. The
!---- last section runs the production pipeline itself, from a geographic
!---- position through location, strain, contraction and conversion, at
!---- perturbed positions on a manufactured element, and compares with
!---- gf_partials_loc.
!----

  program test_gf_partials

  use constants, only: PI,NGLLX,NGLLY,NGLLZ,NGNOD,NDIM,GAUSSALPHA,GAUSSBETA

  use gf_par, only: t_gf_stf,t_gf_taxis,GF_OK,GF_NCOMP, &
                    GF_SRC_CMT,GF_SRC_FORCE,GF_STF_TRUNC,GF_STF_HEAVI,GF_STF_GAUSS

  use gf_shape3D, only: gf_shape3D_functions,gf_shape3D_functions_2nd,gf_shape3D_map,gf_shape3D_map_2nd

  use gf_geometry, only: gf_geographic_to_cartesian,gf_find_local_coords

  use gf_geo_chain, only: gf_spline_derivative,gf_geographic_jacobian

  use gf_interp, only: gf_interp_weights_deriv,gf_interp_weights_deriv2

  use gf_strain, only: GF_VOIGT,GF_XX,GF_YY,GF_ZZ,GF_XY,GF_XZ,GF_YZ, &
                       gf_strain_dweights,gf_strain_ddweights,gf_strain_snapshot,gf_strain_trace_d

  use gf_moment, only: gf_rotate_moment_tensor,gf_rotate_moment_tensor_deriv,gf_moment_contract

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

  ! Earth, for the geographic chain
  double precision, parameter :: R_EARTH = 6371000.d0

  ! a polynomial in three variables of total degree <= 4: coef(i,j,k) x^i y^j z^k
  integer, parameter :: PDEG = 4

  ! context shared with the program-level helpers geo_forward() and
  ! pipeline() of sections 9 and 10 (an internal procedure cannot have
  ! internal procedures of its own)
  integer, parameter :: NSPL = 80
  double precision, dimension(NSPL) :: ct_rspl,ct_ell,ct_ell2
  double precision, dimension(1) :: ct_nospl
  double precision, dimension(3) :: ct_s0
  double precision :: ct_elev0,ct_glat,ct_glon
  logical :: ct_ell_on
  integer, parameter :: NTC = 60
  double precision, dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ,NTC) :: ct_u
  double precision, dimension(NGNOD) :: ct_xelm,ct_yelm,ct_zelm
  double precision, dimension(NGLLX) :: ct_xigll
  double precision, dimension(NGLLY) :: ct_yigll
  double precision, dimension(NGLLZ) :: ct_zigll
  double precision, dimension(6) :: ct_msph
  type(t_gf_stf) :: ct_stf
  type(t_gf_taxis) :: ct_tax
  double precision, dimension(:), allocatable :: ct_w
  integer :: ct_nt

  integer :: nfail

  nfail = 0

  write(*,'(a)') 'test_gf_partials: moment-tensor, centroid-time and centroid-position partials'
  write(*,'(a)') ''

  call test_ndp(nfail)
  call test_unit_kernel(nfail)
  call test_mt_linearity(nfail)
  call test_time_partial(nfail)
  call test_shape_2nd(nfail)
  call test_strain_gradient(nfail)
  call test_mt_deriv(nfail)
  call test_spline_deriv(nfail)
  call test_geo_jacobian(nfail)
  call test_full_chain(nfail)

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
  integer :: v,icomp,it,ierr,nt

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

  ! each partial is the seismogram of its unit tensor: the same library
  ! routines on the same numbers, with the scaling statement between them
  ! written here and in gf_partials_mt -- two compilation units, so a
  ! derived tolerance rather than equality
  worst = 0.d0
  ref = maxval(abs(dp))
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
        worst = max(worst,abs(y(it) - dp(v,icomp,it))/ref)
      enddo
    enddo
  enddo
  call gf_report('dp(v) == seismogram of unit tensor v (rel)',worst,1.d-14,nfail)

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

!
!-------------------------------------------------------------------------------------------------
!
!  Stage 8 helpers
!
!-------------------------------------------------------------------------------------------------
!

  subroutine anchor_reference_coords(xiref)

! reference (xi,eta,gamma) of the 27 anchors, from hex_nodes() -- as in
! test_gf_shape3D, so a node reordering breaks this test rather than
! decoupling it

  implicit none
  double precision, dimension(NDIM,NGNOD), intent(out) :: xiref

  integer, dimension(NGNOD) :: iaddx,iaddy,iaddz
  integer :: ia

  call hex_nodes(iaddx,iaddy,iaddz)
  do ia = 1,NGNOD
    xiref(1,ia) = dble(iaddx(ia)) - 1.d0
    xiref(2,ia) = dble(iaddy(ia)) - 1.d0
    xiref(3,ia) = dble(iaddz(ia)) - 1.d0
  enddo

  end subroutine anchor_reference_coords

!
!-------------------------------------------------------------------------------------------------
!

  subroutine random_affine(scale,x0,amat)

! a random orientation-preserving affine map x = x0 + A xi, with A
! non-symmetric and of size `scale`

  implicit none
  double precision, intent(in) :: scale
  double precision, dimension(NDIM), intent(out) :: x0
  double precision, dimension(NDIM,NDIM), intent(out) :: amat

  double precision :: det
  integer :: i,j

  do i = 1,NDIM
    x0(i) = gf_rand_range(-1.d0,1.d0)
    do j = 1,NDIM
      amat(i,j) = scale*(gf_rand_range(-0.3d0,0.3d0))
    enddo
    amat(i,i) = amat(i,i) + scale
  enddo
  det = amat(1,1)*(amat(2,2)*amat(3,3) - amat(2,3)*amat(3,2)) &
      - amat(1,2)*(amat(2,1)*amat(3,3) - amat(2,3)*amat(3,1)) &
      + amat(1,3)*(amat(2,1)*amat(3,2) - amat(2,2)*amat(3,1))
  if (det < 0.d0) amat(:,1) = -amat(:,1)

  end subroutine random_affine

!
!-------------------------------------------------------------------------------------------------
!

  subroutine q2_anchors(x0,amat,cmat,xiref,xelm,yelm,zelm)

! the 27 anchors of the quadratic map x = x0 + A xi + C(xi,xi), with
! C(xi,xi)_d = SUM_{a<=b} cmat(d,a,b) xi_a xi_b; cmat = 0 is affine

  implicit none
  double precision, dimension(NDIM), intent(in) :: x0
  double precision, dimension(NDIM,NDIM), intent(in) :: amat
  double precision, dimension(NDIM,NDIM,NDIM), intent(in) :: cmat
  double precision, dimension(NDIM,NGNOD), intent(in) :: xiref
  double precision, dimension(NGNOD), intent(out) :: xelm,yelm,zelm

  double precision, dimension(NDIM) :: x
  integer :: ia

  do ia = 1,NGNOD
    call q2_map(x0,amat,cmat,xiref(:,ia),x)
    xelm(ia) = x(1) ; yelm(ia) = x(2) ; zelm(ia) = x(3)
  enddo

  end subroutine q2_anchors

!
!-------------------------------------------------------------------------------------------------
!

  subroutine q2_map(x0,amat,cmat,xi,x)

  implicit none
  double precision, dimension(NDIM), intent(in) :: x0,xi
  double precision, dimension(NDIM,NDIM), intent(in) :: amat
  double precision, dimension(NDIM,NDIM,NDIM), intent(in) :: cmat
  double precision, dimension(NDIM), intent(out) :: x

  integer :: d,a,b

  do d = 1,NDIM
    x(d) = x0(d) + amat(d,1)*xi(1) + amat(d,2)*xi(2) + amat(d,3)*xi(3)
    do a = 1,NDIM
      do b = a,NDIM
        x(d) = x(d) + cmat(d,a,b)*xi(a)*xi(b)
      enddo
    enddo
  enddo

  end subroutine q2_map

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gll_points(xigll,yigll,zigll)

  implicit none
  double precision, dimension(NGLLX), intent(out) :: xigll
  double precision, dimension(NGLLY), intent(out) :: yigll
  double precision, dimension(NGLLZ), intent(out) :: zigll

  double precision, dimension(NGLLX) :: wx
  double precision, dimension(NGLLY) :: wy
  double precision, dimension(NGLLZ) :: wz

  call zwgljd(xigll,wx,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wy,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wz,NGLLZ,GAUSSALPHA,GAUSSBETA)

  end subroutine gll_points

!
!-------------------------------------------------------------------------------------------------
!

  subroutine random_poly(maxdeg,coef)

! random coefficients for the monomials x^i y^j z^k with i+j+k <= maxdeg,
! zero above

  implicit none
  integer, intent(in) :: maxdeg
  double precision, dimension(0:PDEG,0:PDEG,0:PDEG), intent(out) :: coef

  integer :: i,j,k

  coef(:,:,:) = 0.d0
  do k = 0,PDEG
    do j = 0,PDEG
      do i = 0,PDEG
        if (i+j+k <= maxdeg) coef(i,j,k) = gf_rand_range(-1.d0,1.d0)
      enddo
    enddo
  enddo

  end subroutine random_poly

!
!-------------------------------------------------------------------------------------------------
!

  subroutine poly_eval(coef,x,val,grad,hess)

! value, gradient and Hessian of the polynomial, by the plain monomial sum
! -- no basis, no Jacobian: the oracle side

  implicit none
  double precision, dimension(0:PDEG,0:PDEG,0:PDEG), intent(in) :: coef
  double precision, dimension(NDIM), intent(in) :: x
  double precision, intent(out) :: val
  double precision, dimension(NDIM), intent(out) :: grad
  double precision, dimension(NDIM,NDIM), intent(out) :: hess

  integer :: i,j,k
  double precision :: c

  val = 0.d0
  grad(:) = 0.d0
  hess(:,:) = 0.d0
  do k = 0,PDEG
    do j = 0,PDEG
      do i = 0,PDEG
        c = coef(i,j,k)
        if (c == 0.d0) cycle
        val = val + c*x(1)**i*x(2)**j*x(3)**k
        if (i >= 1) grad(1) = grad(1) + c*i*x(1)**(i-1)*x(2)**j*x(3)**k
        if (j >= 1) grad(2) = grad(2) + c*j*x(1)**i*x(2)**(j-1)*x(3)**k
        if (k >= 1) grad(3) = grad(3) + c*k*x(1)**i*x(2)**j*x(3)**(k-1)
        if (i >= 2) hess(1,1) = hess(1,1) + c*i*(i-1)*x(1)**(i-2)*x(2)**j*x(3)**k
        if (j >= 2) hess(2,2) = hess(2,2) + c*j*(j-1)*x(1)**i*x(2)**(j-2)*x(3)**k
        if (k >= 2) hess(3,3) = hess(3,3) + c*k*(k-1)*x(1)**i*x(2)**j*x(3)**(k-2)
        if (i >= 1 .and. j >= 1) hess(1,2) = hess(1,2) + c*i*j*x(1)**(i-1)*x(2)**(j-1)*x(3)**k
        if (i >= 1 .and. k >= 1) hess(1,3) = hess(1,3) + c*i*k*x(1)**(i-1)*x(2)**j*x(3)**(k-1)
        if (j >= 1 .and. k >= 1) hess(2,3) = hess(2,3) + c*j*k*x(1)**i*x(2)**(j-1)*x(3)**(k-1)
      enddo
    enddo
  enddo
  hess(2,1) = hess(1,2)
  hess(3,1) = hess(1,3)
  hess(3,2) = hess(2,3)

  end subroutine poly_eval

!
!-------------------------------------------------------------------------------------------------
!

  subroutine ell_table(n,rspl,ell,ell2,r_at,ell_at,dell_at)

! a manufactured ellipticity table: a smooth ell(r) on [0.55, 1.02] with
! its exact end tangents, built by the solver's own spline_construction

  implicit none
  integer, intent(in) :: n
  double precision, dimension(n), intent(out) :: rspl,ell,ell2
  double precision, intent(in) :: r_at
  double precision, intent(out) :: ell_at,dell_at

  double precision :: yp1,ypn
  integer :: i

  do i = 1,n
    rspl(i) = 0.55d0 + (1.02d0 - 0.55d0)*dble(i-1)/dble(n-1)
    ell(i) = ell_fun(rspl(i))
  enddo
  yp1 = dell_fun(rspl(1))
  ypn = dell_fun(rspl(n))
  call spline_construction(rspl,ell,n,yp1,ypn,ell2)
  ell_at = ell_fun(r_at)
  dell_at = dell_fun(r_at)

  end subroutine ell_table

  double precision function ell_fun(r)
  implicit none
  double precision, intent(in) :: r
  ell_fun = 3.35d-3*(1.d0 + 0.3d0*r**2 - 0.15d0*r**3 + 0.05d0*sin(7.d0*r))
  end function ell_fun

  double precision function dell_fun(r)
  implicit none
  double precision, intent(in) :: r
  dell_fun = 3.35d-3*(0.6d0*r - 0.45d0*r**2 + 0.35d0*cos(7.d0*r))
  end function dell_fun

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_shape_2nd(nfail)

! the second derivatives of the tri-quadratic map: sums to zero, exact
! against a central difference (any step, the functions are quadratic),
! and the derivative of the inverse Jacobian against a Richardson
! difference of gf_shape3D_map on a curved element, zero on an affine one

  implicit none
  integer, intent(inout) :: nfail

  double precision, parameter :: HSTEP = 0.1d0
  double precision, dimension(NDIM,NDIM,NGNOD) :: d2
  double precision, dimension(NGNOD) :: shp,xelm,yelm,zelm
  double precision, dimension(NDIM,NGNOD) :: dp,dm,xiref
  double precision, dimension(NDIM) :: xi,xp,xm,x0,xyz
  double precision, dimension(NDIM,NDIM) :: amat,jinv,jp,jm,jp2,jm2,d1,d2r
  double precision, dimension(NDIM,NDIM,NDIM) :: cmat,djinv
  double precision :: worst_sum,worst_fd,worst_sym,worst_dj,worst_aff,jacobian,h,s,ref
  integer :: itrial,a,b,ia,ierr,nbad

  write(*,'(a)') '5. second derivatives of the element map'

  call anchor_reference_coords(xiref)
  call gf_lcg_seed(808)

  worst_sum = 0.d0
  worst_fd = 0.d0
  worst_sym = 0.d0
  do itrial = 1,200
    xi(1) = gf_rand_range(-1.1d0,1.1d0)
    xi(2) = gf_rand_range(-1.1d0,1.1d0)
    xi(3) = gf_rand_range(-1.1d0,1.1d0)
    call gf_shape3D_functions_2nd(xi(1),xi(2),xi(3),d2)
    do a = 1,NDIM
      do b = 1,NDIM
        s = 0.d0
        do ia = 1,NGNOD
          s = s + d2(a,b,ia)
          worst_sym = max(worst_sym,abs(d2(a,b,ia) - d2(b,a,ia)))
        enddo
        worst_sum = max(worst_sum,abs(s))
      enddo
      ! central difference of the first derivatives along xi_a
      xp(:) = xi(:) ; xp(a) = xp(a) + HSTEP
      xm(:) = xi(:) ; xm(a) = xm(a) - HSTEP
      call gf_shape3D_functions(xp(1),xp(2),xp(3),shp,dp)
      call gf_shape3D_functions(xm(1),xm(2),xm(3),shp,dm)
      do b = 1,NDIM
        do ia = 1,NGNOD
          worst_fd = max(worst_fd,abs((dp(b,ia) - dm(b,ia))/(2.d0*HSTEP) - d2(a,b,ia)))
        enddo
      enddo
    enddo
  enddo
  call gf_report('sum over the 27 anchors of d2shape = 0',worst_sum,1.d-13,nfail)
  call gf_report('d2shape symmetric in (a,b)             ',worst_sym,0.d0,nfail)
  call gf_report('d2shape == central difference of dshape',worst_fd,1.d-13,nfail)

  ! the inverse Jacobian's derivative on curved Q2 elements, against a
  ! Richardson-extrapolated central difference of gf_shape3D_map's jinv
  worst_dj = 0.d0
  worst_aff = 0.d0
  nbad = 0
  do itrial = 1,20
    call random_affine(1.d0,x0,amat)
    do a = 1,NDIM
      do b = 1,NDIM
        cmat(:,a,b) = 0.d0
        if (b >= a) then
          cmat(1,a,b) = gf_rand_range(-0.08d0,0.08d0)
          cmat(2,a,b) = gf_rand_range(-0.08d0,0.08d0)
          cmat(3,a,b) = gf_rand_range(-0.08d0,0.08d0)
        endif
      enddo
    enddo
    call q2_anchors(x0,amat,cmat,xiref,xelm,yelm,zelm)
    xi(1) = gf_rand_range(-0.9d0,0.9d0)
    xi(2) = gf_rand_range(-0.9d0,0.9d0)
    xi(3) = gf_rand_range(-0.9d0,0.9d0)
    call gf_shape3D_map_2nd(xelm,yelm,zelm,xi(1),xi(2),xi(3),xyz,jinv,jacobian,djinv,ierr)
    if (ierr /= GF_OK) then
      nbad = nbad + 1
      cycle
    endif
    ref = maxval(abs(jinv))
    do a = 1,NDIM
      h = 1.d-2
      xp(:) = xi(:) ; xp(a) = xp(a) + h
      xm(:) = xi(:) ; xm(a) = xm(a) - h
      call gf_shape3D_map(xelm,yelm,zelm,xp(1),xp(2),xp(3),xyz,jp,jacobian,ierr)
      call gf_shape3D_map(xelm,yelm,zelm,xm(1),xm(2),xm(3),xyz,jm,jacobian,ierr)
      d1 = (jp - jm)/(2.d0*h)
      xp(:) = xi(:) ; xp(a) = xp(a) + h/2.d0
      xm(:) = xi(:) ; xm(a) = xm(a) - h/2.d0
      call gf_shape3D_map(xelm,yelm,zelm,xp(1),xp(2),xp(3),xyz,jp2,jacobian,ierr)
      call gf_shape3D_map(xelm,yelm,zelm,xm(1),xm(2),xm(3),xyz,jm2,jacobian,ierr)
      d2r = (4.d0*(jp2 - jm2)/h - d1)/3.d0
      worst_dj = max(worst_dj,maxval(abs(d2r - djinv(:,:,a)))/ref)
    enddo

    ! the same element without curvature: djinv vanishes
    cmat(:,:,:) = 0.d0
    call q2_anchors(x0,amat,cmat,xiref,xelm,yelm,zelm)
    call gf_shape3D_map_2nd(xelm,yelm,zelm,xi(1),xi(2),xi(3),xyz,jinv,jacobian,djinv,ierr)
    worst_aff = max(worst_aff,maxval(abs(djinv))/maxval(abs(jinv)))
  enddo
  call gf_report_true('curved Q2 elements evaluated          ',nbad == 0,nfail)
  call gf_report('djinv vs Richardson FD of jinv, curved ',worst_dj,1.d-9,nfail)
  call gf_report('djinv = 0 on an affine element         ',worst_aff,1.d-12,nfail)

  end subroutine test_shape_2nd

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_strain_gradient(nfail)

! d eps/dx from the differentiated weight table, against the closed-form
! second derivatives of a polynomial field: exact on an affine element for
! total degree <= 4, and on a curved Q2 element for a quadratic field

  implicit none
  integer, intent(inout) :: nfail

  double precision, dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ) :: u
  double precision, dimension(0:PDEG,0:PDEG,0:PDEG,GF_NCOMP,GF_NCOMP) :: coef
  double precision, dimension(NGLLX) :: xigll,hxi,hpxi,hppxi
  double precision, dimension(NGLLY) :: yigll,heta,hpeta,hppeta
  double precision, dimension(NGLLZ) :: zigll,hgam,hpgam,hppgam
  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM) :: dw
  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM,NDIM) :: ddw
  double precision, dimension(GF_VOIGT,GF_NCOMP) :: eps
  double precision, dimension(GF_VOIGT,GF_NCOMP,NDIM) :: deps_xi,deps_x,expect
  double precision, dimension(NDIM,NDIM,NDIM) :: cmat,djinv,hess
  double precision, dimension(NDIM,NDIM) :: amat,jinv,hp
  double precision, dimension(NDIM,NGNOD) :: xiref
  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  double precision, dimension(NDIM) :: x0,xi,xnode,xs,xyz,grad,xloc
  double precision :: worst,jacobian,val,scale,ref
  integer :: icase,a,p,m,b,i,j,k,ierr

  write(*,'(a)') '6. strain gradient vs the polynomial oracle'

  call anchor_reference_coords(xiref)
  call gll_points(xigll,yigll,zigll)
  call gf_lcg_seed(2718)

  do icase = 1,4

    ! 1: affine, unit scale, linear field (gradient of the strain is zero)
    ! 2: affine, unit scale, random total-degree-4 field
    ! 3: affine at realistic scale (|x0| = 1, half-width 1e-2), degree 4
    ! 4: curved Q2 element, quadratic field (constant strain gradient)
    select case (icase)
    case (1,2)
      scale = 1.d0
    case (3)
      scale = 1.d-2
    case (4)
      scale = 1.d0
    end select
    call random_affine(scale,x0,amat)
    if (icase /= 3) x0(:) = 0.d0
    cmat(:,:,:) = 0.d0
    if (icase == 4) then
      do a = 1,NDIM
        do b = a,NDIM
          cmat(1,a,b) = gf_rand_range(-0.08d0,0.08d0)
          cmat(2,a,b) = gf_rand_range(-0.08d0,0.08d0)
          cmat(3,a,b) = gf_rand_range(-0.08d0,0.08d0)
        enddo
      enddo
    endif
    call q2_anchors(x0,amat,cmat,xiref,xelm,yelm,zelm)

    ! the field, per force and displacement component, in the local
    ! coordinate (x - x0)/scale so that the values stay O(1)
    do a = 1,GF_NCOMP
      do p = 1,GF_NCOMP
        select case (icase)
        case (1)
          call random_poly(1,coef(:,:,:,a,p))
        case (4)
          call random_poly(2,coef(:,:,:,a,p))
        case default
          call random_poly(4,coef(:,:,:,a,p))
        end select
      enddo
    enddo

    ! nodal values at the GLL points of the element
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          xi(1) = xigll(i) ; xi(2) = yigll(j) ; xi(3) = zigll(k)
          call q2_map(x0,amat,cmat,xi,xnode)
          xloc(:) = (xnode(:) - x0(:))/scale
          do a = 1,GF_NCOMP
            do p = 1,GF_NCOMP
              call poly_eval(coef(:,:,:,a,p),xloc,val,grad,hp)
              u(a,p,i,j,k) = val
            enddo
          enddo
        enddo
      enddo
    enddo

    ! the evaluation point, its physical position, and the routines under test
    xi(1) = gf_rand_range(-0.9d0,0.9d0)
    xi(2) = gf_rand_range(-0.9d0,0.9d0)
    xi(3) = gf_rand_range(-0.9d0,0.9d0)
    call q2_map(x0,amat,cmat,xi,xs)

    call gf_shape3D_map_2nd(xelm,yelm,zelm,xi(1),xi(2),xi(3),xyz,jinv,jacobian,djinv,ierr)
    call gf_interp_weights_deriv2(xi(1),xi(2),xi(3),hxi,hpxi,hppxi,heta,hpeta,hppeta,hgam,hpgam,hppgam)
    call gf_strain_dweights(hxi,hpxi,heta,hpeta,hgam,hpgam,jinv,dw)
    call gf_strain_ddweights(hxi,hpxi,hppxi,heta,hpeta,hppeta,hgam,hpgam,hppgam,jinv,djinv,ddw)

    call gf_strain_snapshot(u,dw,eps)
    do b = 1,NDIM
      call gf_strain_snapshot(u,ddw(:,:,:,:,b),deps_xi(:,:,b))
    enddo
    do m = 1,NDIM
      deps_x(:,:,m) = jinv(1,m)*deps_xi(:,:,1) + jinv(2,m)*deps_xi(:,:,2) + jinv(3,m)*deps_xi(:,:,3)
    enddo

    ! the oracle: d eps_pq / dx_m = (H_p(q,m) + H_q(p,m))/2 from the
    ! polynomial's Hessian, with the 1/scale^2 of the local coordinate
    xloc(:) = (xs(:) - x0(:))/scale
    do a = 1,GF_NCOMP
      do p = 1,GF_NCOMP
        call poly_eval(coef(:,:,:,a,p),xloc,val,grad,hess(:,:,p))
        hess(:,:,p) = hess(:,:,p)/scale**2
      enddo
      do m = 1,NDIM
        expect(GF_XX,a,m) = hess(1,m,1)
        expect(GF_YY,a,m) = hess(2,m,2)
        expect(GF_ZZ,a,m) = hess(3,m,3)
        expect(GF_XY,a,m) = 0.5d0*(hess(2,m,1) + hess(1,m,2))
        expect(GF_XZ,a,m) = 0.5d0*(hess(3,m,1) + hess(1,m,3))
        expect(GF_YZ,a,m) = 0.5d0*(hess(3,m,2) + hess(2,m,3))
      enddo
    enddo

    ref = maxval(abs(expect))
    if (icase == 1) then
      worst = maxval(abs(deps_x))/maxval(abs(eps))
      call gf_report('linear field: d eps/dx = 0 (rel. to eps)',worst,1.d-13,nfail)
    else
      worst = maxval(abs(deps_x - expect))/ref
      select case (icase)
      case (2)
        call gf_report('affine, random total-degree-4 field    ',worst,1.d-12,nfail)
      case (3)
        call gf_report('affine at realistic scale (1e-2)       ',worst,1.d-10,nfail)
      case (4)
        call gf_report('curved Q2, quadratic field             ',worst,1.d-11,nfail)
      end select
    endif

  enddo

  end subroutine test_strain_gradient

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_mt_deriv(nfail)

! d M_cart / d theta and d phi against a Richardson central difference of
! gf_rotate_moment_tensor, over random orientations and tensors

  implicit none
  integer, intent(inout) :: nfail

  double precision, dimension(6) :: m_sph
  double precision, dimension(NDIM,NDIM) :: dmt,dmp,mp,mm,d1,d2,mc
  double precision :: theta,phi,h,worst_t,worst_p,ref
  integer :: itrial,v

  write(*,'(a)') '7. moment-tensor rotation derivative'

  call gf_lcg_seed(3141)
  worst_t = 0.d0
  worst_p = 0.d0
  h = 1.d-3
  do itrial = 1,30
    theta = gf_rand_range(0.15d0,PI - 0.15d0)
    phi = gf_rand_range(0.d0,2.d0*PI)
    do v = 1,6
      m_sph(v) = gf_rand_range(-1.d0,1.d0)
    enddo
    call gf_rotate_moment_tensor(theta,phi,m_sph,mc)
    ref = maxval(abs(mc))
    call gf_rotate_moment_tensor_deriv(theta,phi,m_sph,dmt,dmp)

    call gf_rotate_moment_tensor(theta + h,phi,m_sph,mp)
    call gf_rotate_moment_tensor(theta - h,phi,m_sph,mm)
    d1 = (mp - mm)/(2.d0*h)
    call gf_rotate_moment_tensor(theta + h/2.d0,phi,m_sph,mp)
    call gf_rotate_moment_tensor(theta - h/2.d0,phi,m_sph,mm)
    d2 = (4.d0*(mp - mm)/h - d1)/3.d0
    worst_t = max(worst_t,maxval(abs(d2 - dmt))/ref)

    call gf_rotate_moment_tensor(theta,phi + h,m_sph,mp)
    call gf_rotate_moment_tensor(theta,phi - h,m_sph,mm)
    d1 = (mp - mm)/(2.d0*h)
    call gf_rotate_moment_tensor(theta,phi + h/2.d0,m_sph,mp)
    call gf_rotate_moment_tensor(theta,phi - h/2.d0,m_sph,mm)
    d2 = (4.d0*(mp - mm)/h - d1)/3.d0
    worst_p = max(worst_p,maxval(abs(d2 - dmp))/ref)
  enddo
  call gf_report('dM/dtheta vs Richardson FD of the rotation',worst_t,1.d-10,nfail)
  call gf_report('dM/dphi   vs Richardson FD of the rotation',worst_p,1.d-10,nfail)

  end subroutine test_mt_deriv

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_spline_deriv(nfail)

! the spline's derivative against a Richardson difference of
! spline_evaluation on a manufactured table

  implicit none
  integer, intent(inout) :: nfail

  integer, parameter :: N = 80
  double precision, dimension(N) :: rspl,ell,ell2
  double precision :: r,dy,yp,ym,d1,d2,h,worst,worst_fun,ell_at,dell_at,ref
  integer :: itrial,ierr

  write(*,'(a)') '8. spline derivative'

  call ell_table(N,rspl,ell,ell2,1.d0,ell_at,dell_at)
  call gf_lcg_seed(577)
  worst = 0.d0
  worst_fun = 0.d0
  ref = abs(dell_fun(0.8d0))
  h = 1.d-4
  do itrial = 1,50
    r = gf_rand_range(0.56d0,1.01d0)
    call gf_spline_derivative(rspl,ell,ell2,N,r,dy,ierr)
    call spline_evaluation(rspl,ell,ell2,N,r + h,yp)
    call spline_evaluation(rspl,ell,ell2,N,r - h,ym)
    d1 = (yp - ym)/(2.d0*h)
    call spline_evaluation(rspl,ell,ell2,N,r + h/2.d0,yp)
    call spline_evaluation(rspl,ell,ell2,N,r - h/2.d0,ym)
    d2 = (4.d0*(yp - ym)/h - d1)/3.d0
    worst = max(worst,abs(dy - d2)/ref)
    worst_fun = max(worst_fun,abs(dy - dell_fun(r))/ref)
  enddo
  call gf_report_true('gf_spline_derivative returns GF_OK        ',ierr == GF_OK,nfail)
  call gf_report('spline derivative vs Richardson FD        ',worst,1.d-9,nfail)
  write(*,'(a,es10.3,a)') '     (vs the function it interpolates: ',worst_fun,', the table''s own error)'

  end subroutine test_spline_deriv

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_geo_jacobian(nfail)

! d(x,y,z)/d(lat,lon,depth) against a Richardson difference of the forward
! map: spherical; elliptical with the manufactured table; and with a
! synthetic elevation field whose gradient is known

  implicit none
  integer, intent(inout) :: nfail

  double precision, dimension(NDIM,3) :: dxds,fd,fd2
  double precision, dimension(NDIM) :: xyz,xp,xm,e_r
  double precision, dimension(3) :: sp,sm,hs
  double precision :: theta,phi,r_surface,dtheta_dlat,dphi_dlon,worst,ref,ell_at,dell_at,worst_dep,worst_exact
  integer :: icase,ia,ierr,itrial

  write(*,'(a)') '9. geographic chain rule'

  call ell_table(NSPL,ct_rspl,ct_ell,ct_ell2,1.d0,ell_at,dell_at)
  ct_nospl(1) = 0.d0
  hs(1) = 1.d-3 ; hs(2) = 1.d-3 ; hs(3) = 5.d-2
  call gf_lcg_seed(1618)

  do icase = 1,3
    ct_ell_on = (icase >= 2)
    worst = 0.d0
    worst_dep = 0.d0
    worst_exact = 0.d0
    do itrial = 1,6
      if (itrial == 1) then
        ct_s0(1) = -5.812d0 ; ct_s0(2) = -75.27d0 ; ct_s0(3) = 122.6d0
      else
        ct_s0(1) = gf_rand_range(-80.d0,80.d0)
        ct_s0(2) = gf_rand_range(-180.d0,180.d0)
        ct_s0(3) = gf_rand_range(1.d0,600.d0)
      endif
      ! a synthetic elevation with a known gradient, in metres per degree
      if (icase == 3) then
        ct_elev0 = 800.d0 ; ct_glat = 1500.d0 ; ct_glon = -900.d0
      else
        ct_elev0 = 0.d0 ; ct_glat = 0.d0 ; ct_glon = 0.d0
      endif

      if (ct_ell_on) then
        call gf_geographic_jacobian(ct_s0(1),ct_s0(2),ct_s0(3),ct_ell_on,ct_elev0,ct_glat,ct_glon, &
                                    NSPL,ct_rspl,ct_ell,ct_ell2,R_EARTH,dxds,dtheta_dlat,dphi_dlon,ierr)
      else
        call gf_geographic_jacobian(ct_s0(1),ct_s0(2),ct_s0(3),ct_ell_on,ct_elev0,ct_glat,ct_glon, &
                                    0,ct_nospl,ct_nospl,ct_nospl,R_EARTH,dxds,dtheta_dlat,dphi_dlon,ierr)
      endif
      if (ierr /= GF_OK) then
        call gf_report_true('gf_geographic_jacobian returns GF_OK',.false.,nfail)
        return
      endif

      ! Richardson central differences of the forward map
      do ia = 1,3
        sp(:) = ct_s0(:) ; sp(ia) = sp(ia) + hs(ia)
        sm(:) = ct_s0(:) ; sm(ia) = sm(ia) - hs(ia)
        call geo_forward(sp,xp,theta,phi,r_surface)
        call geo_forward(sm,xm,theta,phi,r_surface)
        fd(:,ia) = (xp(:) - xm(:))/(2.d0*hs(ia))
        sp(:) = ct_s0(:) ; sp(ia) = sp(ia) + hs(ia)/2.d0
        sm(:) = ct_s0(:) ; sm(ia) = sm(ia) - hs(ia)/2.d0
        call geo_forward(sp,xp,theta,phi,r_surface)
        call geo_forward(sm,xm,theta,phi,r_surface)
        fd2(:,ia) = (xp(:) - xm(:))/hs(ia)
        fd(:,ia) = (4.d0*fd2(:,ia) - fd(:,ia))/3.d0
      enddo
      ref = maxval(abs(dxds(:,1:2)))
      worst = max(worst,maxval(abs(fd(:,1:2) - dxds(:,1:2)))/ref)
      worst_dep = max(worst_dep,maxval(abs(fd(:,3) - dxds(:,3)))/maxval(abs(dxds(:,3))))

      ! depth: exactly -(1000/R) e_r
      call geo_forward(ct_s0,xyz,theta,phi,r_surface)
      e_r(:) = xyz(:)/sqrt(sum(xyz**2))
      worst_exact = max(worst_exact,maxval(abs(dxds(:,3) + (1000.d0/R_EARTH)*e_r(:)))/(1000.d0/R_EARTH))
    enddo
    select case (icase)
    case (1)
      call gf_report('spherical: dx/dlat, dx/dlon vs Richardson FD',worst,1.d-9,nfail)
    case (2)
      call gf_report('elliptical: dx/dlat, dx/dlon vs Richardson FD',worst,1.d-9,nfail)
    case (3)
      call gf_report('with topography gradient: vs Richardson FD  ',worst,1.d-9,nfail)
    end select
    call gf_report('  dx/ddepth vs Richardson FD                ',worst_dep,1.d-9,nfail)
    call gf_report('  dx/ddepth == -(1000/R) e_r, closed form   ',worst_exact,1.d-12,nfail)
  enddo

  end subroutine test_geo_jacobian

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_full_chain(nfail)

! the whole chain: the production pipeline -- geographic map, Newton
! location in a manufactured element, weights, strain, contraction with
! the rotated moment tensor, conversion -- run at perturbed positions and
! Richardson-differenced, against gf_partials_loc at the base position

  implicit none
  integer, intent(inout) :: nfail

  double precision, dimension(0:PDEG,0:PDEG,0:PDEG,GF_NCOMP,GF_NCOMP) :: coef
  double precision, dimension(NGLLX) :: hxi,hpxi,hppxi
  double precision, dimension(NGLLY) :: heta,hpeta,hppeta
  double precision, dimension(NGLLZ) :: hgam,hpgam,hppgam
  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM) :: dw
  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM,NDIM) :: ddw
  double precision, dimension(:,:,:), allocatable :: eps,dp,d1,dr
  double precision, dimension(:,:), allocatable :: seis_p,seis_m,seis_p2,seis_m2
  double precision, dimension(:,:,:,:), allocatable :: deps
  double precision, dimension(NDIM,NDIM,NDIM) :: cmat,djinv
  double precision, dimension(NDIM,NDIM) :: amat,jinv,m_cart,dm_dtheta,dm_dphi,hp
  double precision, dimension(NDIM,3) :: dxds
  double precision, dimension(NDIM,NGNOD) :: xiref
  double precision, dimension(NDIM) :: xc,xi0,xyz0,xyz,xnode,xloc,grad
  double precision, dimension(3) :: sp,sm,hs
  double precision :: theta0,phi0,r_surface,scale,val,jacobian,t,q,worst,ref
  double precision :: dtheta_dlat,dphi_dlon,ell_at,dell_at
  integer :: a,p,i,j,k,it,ia,ierr,b
  character(len=8), dimension(3), parameter :: pname = (/ 'lat     ','lon     ','depth   ' /)

  write(*,'(a)') '10. the whole chain: gf_partials_loc vs Richardson FD of the pipeline'

  call anchor_reference_coords(xiref)
  call gll_points(ct_xigll,ct_yigll,ct_zigll)
  call ell_table(NSPL,ct_rspl,ct_ell,ct_ell2,1.d0,ell_at,dell_at)
  ct_ell_on = .true.
  call gf_lcg_seed(4242)

  ! the source, on an elliptical planet with a synthetic elevation field
  ct_s0(1) = -5.812d0 ; ct_s0(2) = -75.27d0 ; ct_s0(3) = 122.6d0
  ct_elev0 = 800.d0 ; ct_glat = 1500.d0 ; ct_glon = -900.d0
  call geo_forward(ct_s0,xyz0,theta0,phi0,r_surface)

  ! an affine element of half-width 1e-2 (64 km) with the source inside it
  ! but off centre, and a curvature-free geometry so that the interpolant
  ! is exact for the quartic field
  call random_affine(1.d-2,xc,amat)
  xi0(1) = -0.1d0 ; xi0(2) = 0.2d0 ; xi0(3) = -0.15d0
  xc(:) = xyz0(:) - (amat(:,1)*xi0(1) + amat(:,2)*xi0(2) + amat(:,3)*xi0(3))
  cmat(:,:,:) = 0.d0
  call q2_anchors(xc,amat,cmat,xiref,ct_xelm,ct_yelm,ct_zelm)

  ! the field: a quartic per component in the local coordinate, times a
  ! smooth pulse in time
  scale = 1.d-2
  do a = 1,GF_NCOMP
    do p = 1,GF_NCOMP
      call random_poly(4,coef(:,:,:,a,p))
    enddo
  enddo
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        xnode(1) = ct_xigll(i) ; xnode(2) = ct_yigll(j) ; xnode(3) = ct_zigll(k)
        call q2_map(xc,amat,cmat,xnode,xloc)
        xloc(:) = (xloc(:) - xc(:))/scale
        do a = 1,GF_NCOMP
          do p = 1,GF_NCOMP
            call poly_eval(coef(:,:,:,a,p),xloc,val,grad,hp)
            do it = 1,NTC
              t = (dble(it*SS) - 1.d0)*0.1d0 - T0DB
              q = exp(-((t - 90.d0)/35.d0)**2)
              ct_u(a,p,i,j,k,it) = 1.d-6*val*q
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  ! the moment tensor, the plan, the axis
  do a = 1,6
    ct_msph(a) = gf_rand_range(-1.d0,1.d0)
  enddo
  call gf_stf_plan(GF_SRC_CMT,0,HCMT,HDB,DTR,GF_STF_TRUNC,ct_stf,ierr)
  call gf_taxis_plan(NTC,0.1d0,SS,T0DB,90.d0,ct_tax,ierr)
  ct_nt = ct_tax%nt
  allocate(ct_w(-ct_stf%khalf:ct_stf%khalf),eps(GF_VOIGT,GF_NCOMP,NTC),deps(GF_VOIGT,GF_NCOMP,NTC,NDIM), &
           seis_p(GF_NCOMP,ct_nt),seis_m(GF_NCOMP,ct_nt),seis_p2(GF_NCOMP,ct_nt),seis_m2(GF_NCOMP,ct_nt), &
           dp(3,GF_NCOMP,ct_nt),d1(3,GF_NCOMP,ct_nt),dr(3,GF_NCOMP,ct_nt))
  call gf_stf_kernel(ct_stf,ct_tax%dt_sub,ct_w)

  !--- the analytic route at s0 --------------------------------------------

  call gf_find_local_coords(ct_xelm,ct_yelm,ct_zelm,ct_xigll,ct_yigll,ct_zigll,xyz0,3,3,3, &
                            xi0(1),xi0(2),xi0(3),xyz,jinv,jacobian,ierr)
  call gf_shape3D_map_2nd(ct_xelm,ct_yelm,ct_zelm,xi0(1),xi0(2),xi0(3),xyz,jinv,jacobian,djinv,ierr)
  call gf_interp_weights_deriv2(xi0(1),xi0(2),xi0(3),hxi,hpxi,hppxi,heta,hpeta,hppeta,hgam,hpgam,hppgam)
  call gf_strain_dweights(hxi,hpxi,heta,hpeta,hgam,hpgam,jinv,dw)
  call gf_strain_ddweights(hxi,hpxi,hppxi,heta,hpeta,hppeta,hgam,hpgam,hppgam,jinv,djinv,ddw)
  call gf_strain_trace_d(ct_u,dw,NTC,eps)
  do b = 1,NDIM
    call gf_strain_trace_d(ct_u,ddw(:,:,:,:,b),NTC,deps(:,:,:,b))
  enddo
  call gf_rotate_moment_tensor(theta0,phi0,ct_msph,m_cart)
  call gf_rotate_moment_tensor_deriv(theta0,phi0,ct_msph,dm_dtheta,dm_dphi)
  call gf_geographic_jacobian(ct_s0(1),ct_s0(2),ct_s0(3),.true.,ct_elev0,ct_glat,ct_glon, &
                              NSPL,ct_rspl,ct_ell,ct_ell2,R_EARTH,dxds,dtheta_dlat,dphi_dlon,ierr)
  call gf_partials_loc(eps,deps,NTC,m_cart,dm_dtheta,dm_dphi,dtheta_dlat,dphi_dlon,jinv,dxds, &
                       1.d0,ct_tax,ct_stf,ct_w,dp,ierr)
  call gf_report_true('gf_partials_loc returns GF_OK              ',ierr == GF_OK,nfail)

  !--- the finite difference of the production pipeline ----------------------

  hs(1) = 1.d-3 ; hs(2) = 1.d-3 ; hs(3) = 5.d-2
  do ia = 1,3
    sp(:) = ct_s0(:) ; sp(ia) = sp(ia) + hs(ia)
    sm(:) = ct_s0(:) ; sm(ia) = sm(ia) - hs(ia)
    call pipeline(sp,seis_p)
    call pipeline(sm,seis_m)
    sp(:) = ct_s0(:) ; sp(ia) = sp(ia) + hs(ia)/2.d0
    sm(:) = ct_s0(:) ; sm(ia) = sm(ia) - hs(ia)/2.d0
    call pipeline(sp,seis_p2)
    call pipeline(sm,seis_m2)
    d1(ia,:,:) = (seis_p(:,:) - seis_m(:,:))/(2.d0*hs(ia))
    dr(ia,:,:) = (4.d0*(seis_p2(:,:) - seis_m2(:,:))/hs(ia) - d1(ia,:,:))/3.d0
    ref = maxval(abs(dp(ia,:,:)))
    worst = maxval(abs(dr(ia,:,:) - dp(ia,:,:)))/ref
    write(*,'(a,a,a,es10.3,a,es10.3)') '     ',trim(pname(ia)),': FD(h) vs analytic ', &
          maxval(abs(d1(ia,:,:) - dp(ia,:,:)))/ref,'   Richardson vs analytic ',worst
    call gf_report('d seis/d '//trim(pname(ia))//' vs Richardson FD of the pipeline', &
                   worst,1.d-8,nfail)
  enddo

  deallocate(ct_w,eps,deps,seis_p,seis_m,seis_p2,seis_m2,dp,d1,dr)

  end subroutine test_full_chain

!
!-------------------------------------------------------------------------------------------------
!

  subroutine geo_forward(s,xyz,theta,phi,r_surface)

! the forward geographic map with the synthetic elevation field of the
! program-level context: elev = elev0 + g_lat (lat - lat0) + g_lon (lon - lon0)

  implicit none
  double precision, dimension(3), intent(in) :: s
  double precision, dimension(NDIM), intent(out) :: xyz
  double precision, intent(out) :: theta,phi,r_surface

  double precision :: elev
  integer :: ierr

  elev = ct_elev0 + ct_glat*(s(1) - ct_s0(1)) + ct_glon*(s(2) - ct_s0(2))
  if (ct_ell_on) then
    call gf_geographic_to_cartesian(s(1),s(2),s(3),ct_ell_on,elev,NSPL,ct_rspl,ct_ell,ct_ell2,R_EARTH, &
                                    xyz,theta,phi,r_surface,ierr)
  else
    call gf_geographic_to_cartesian(s(1),s(2),s(3),ct_ell_on,elev,0,ct_nospl,ct_nospl,ct_nospl,R_EARTH, &
                                    xyz,theta,phi,r_surface,ierr)
  endif

  end subroutine geo_forward

!
!-------------------------------------------------------------------------------------------------
!

  subroutine pipeline(s,seis)

! the production sequence on the manufactured element of section 10:
! geographic map, Newton location, weights, strain, contraction with the
! rotated moment tensor, conversion

  implicit none
  double precision, dimension(3), intent(in) :: s
  double precision, dimension(GF_NCOMP,ct_nt), intent(out) :: seis

  double precision, dimension(NGLLX) :: hx,hpx
  double precision, dimension(NGLLY) :: hy,hpy
  double precision, dimension(NGLLZ) :: hz,hpz
  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM) :: dwl
  double precision, dimension(NDIM,NDIM) :: jl,mcl
  double precision, dimension(NDIM) :: xt,xm,xil
  double precision, dimension(GF_VOIGT,GF_NCOMP,NTC) :: epsl
  double precision, dimension(NTC) :: trace
  double precision, dimension(ct_nt) :: xpad,y
  double precision, dimension(0:ct_nt) :: pp
  double precision :: th,ph,rs,jac
  integer :: ic,itl,ie

  call geo_forward(s,xt,th,ph,rs)
  call gf_find_local_coords(ct_xelm,ct_yelm,ct_zelm,ct_xigll,ct_yigll,ct_zigll,xt,3,3,3, &
                            xil(1),xil(2),xil(3),xm,jl,jac,ie)
  call gf_interp_weights_deriv(xil(1),xil(2),xil(3),hx,hpx,hy,hpy,hz,hpz)
  call gf_strain_dweights(hx,hpx,hy,hpy,hz,hpz,jl,dwl)
  call gf_strain_trace_d(ct_u,dwl,NTC,epsl)
  call gf_rotate_moment_tensor(th,ph,ct_msph,mcl)
  do ic = 1,GF_NCOMP
    do itl = 1,NTC
      call gf_moment_contract(mcl,epsl(:,ic,itl),trace(itl))
    enddo
    call gf_pad_left(trace,NTC,ct_tax%npad,xpad)
    call gf_cumsum(xpad,ct_nt,pp)
    call gf_stf_apply(ct_stf,ct_tax%dt_sub,ct_w,pp,xpad,ct_nt,y)
    seis(ic,:) = y(:)
  enddo

  end subroutine pipeline

  end program test_gf_partials
