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
!---- test_gf_stf: the source time function operators and the output axis
!----
!---- Everything here is closed form. The oracles are specfem's own erf
!---- (netlib_specfun_erf, the function comp_source_time_function_heavi is
!---- built on) and its analytic source spectrum (comp_source_spectrum),
!---- plus the Gaussian composition identity
!----
!----   g_h1 * g_h2 = g_h3,   h3^2 = h1^2 + h2^2
!----
!---- which *is* the hdur_corr relation Stage 5 rests on, so the test pins
!---- the physical premise and not merely a convolution routine. The two
!---- widths used throughout are the shipped examples' own: hdur_db =
!---- 6.94968291528492 s and the CMTSOLUTION's 60/1.628 s, on both stored
!---- grids (0.4 s global, 3.4 s regional).
!----
!---- The one thing this file cannot show is that the conversion is the
!---- right one for the database -- that integrate-and-smooth is what turns
!---- the stored field into a forward CMT trace, and that t0, dt_sub and
!---- subsample_step were read correctly. Only the forward-run comparison
!---- does that (utils/green_function/gf_compare.py).
!----

  program test_gf_stf

  use constants, only: PI,SOURCE_DECAY_MIMIC_TRIANGLE

  use gf_par, only: t_gf_stf,t_gf_taxis,GF_OK,GF_ERR_ARG, &
                    GF_STF_NONE,GF_STF_GAUSS,GF_STF_HEAVI,GF_STF_TRUNC, &
                    GF_SRC_FORCE,GF_SRC_CMT

  use gf_stf

  use gf_manufactured

  implicit none

  ! the oracles, from src/specfem3D/
  double precision, external :: netlib_specfun_erf,comp_source_spectrum

  ! the shipped examples' widths and grids
  double precision, parameter :: HDB  = 6.94968291528492d0        ! station attribute hdur
  double precision, parameter :: HCMT = 60.d0                      ! CMTSOLUTION half duration
  double precision, parameter :: T0DB = 34.74841457642461d0        ! mesh_info.h5 t0
  double precision, parameter :: DTG  = 0.4d0, DTR = 3.4d0         ! stored spacings

  integer :: nfail

  nfail = 0

  write(*,'(a)') 'test_gf_stf: source time function operators and output axis'
  write(*,'(a)') ''

  call test_hdur_pin(nfail)
  call test_gauss_kernel(nfail)
  call test_heavi_kernel(nfail)
  call test_delta(nfail)
  call test_centroid(nfail)
  call test_composition(nfail)
  call test_trapezoid(nfail)
  call test_heavi_limit(nfail)
  call test_premise_and_routes(nfail)
  call test_guard(nfail)
  call test_force_kinds(nfail)
  call test_taxis(nfail)
  call test_padding(nfail)
  call test_neumaier(nfail)
  call test_onset(nfail)
  call test_errors(nfail)

  write(*,'(a)') ''
  if (nfail > 0) then
    write(*,'(a,i0,a)') 'test_gf_stf: ',nfail,' assertion(s) FAILED'
    write(0,'(a,i0,a)') 'test_gf_stf: ',nfail,' assertion(s) FAILED'
    stop 1
  endif
  write(*,'(a)') 'test_gf_stf: all assertions passed'

  contains

!
!-------------------------------------------------------------------------------------------------
!

  double precision function gauss(t,h)

! specfem's unit-area Gaussian, comp_source_time_function.f90:202

  implicit none
  double precision, intent(in) :: t,h

  gauss = exp(-(t/h)**2) / (h*sqrt(PI))

  end function gauss

!
!-------------------------------------------------------------------------------------------------
!

  double precision function heavi(t,h)

! specfem's quasi-Heaviside, comp_source_time_function.f90:80, on its own erf

  implicit none
  double precision, intent(in) :: t,h

  heavi = 0.5d0*(1.d0 + netlib_specfun_erf(t/h))

  end function heavi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine axis(n,dt,tstart,t)

  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: dt,tstart
  double precision, dimension(n), intent(out) :: t
  integer :: i

  do i = 1,n
    t(i) = tstart + dble(i-1)*dt
  enddo

  end subroutine axis

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_hdur_pin(nfail)

! the 1.628 placement, as arithmetic

  implicit none
  integer, intent(inout) :: nfail

  write(*,'(a)') '1. hdur_Gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE'

  call gf_report_true('gf_hdur_gaussian(60) == 60/1.628, bitwise', &
                      gf_hdur_gaussian(HCMT) == HCMT/SOURCE_DECAY_MIMIC_TRIANGLE,nfail)

  ! the force forward run prints 'Gaussian half duration: 27.641277641277643'
  ! for its f0 = 45 (forward/OUTPUT_FILES/output_solver.txt:142)
  call gf_report('gf_hdur_gaussian(45) vs the solver print 27.641277641277643', &
                 abs(gf_hdur_gaussian(45.d0) - 27.641277641277643d0)/27.641277641277643d0,1.d-14,nfail)

  end subroutine test_hdur_pin

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_gauss_kernel(nfail)

! the sampled Gaussian: unit sum, the width pin, and the analytic spectrum

  implicit none
  integer, intent(inout) :: nfail

  double precision, dimension(:), allocatable :: w
  double precision :: h,s0,s2,om,dtft,ref
  integer :: k,j,i
  double precision, dimension(3), parameter :: OMEGA = (/ 0.01d0, 0.05d0, 0.1d0 /)

  write(*,'(a)') '2. Gaussian kernel (h = 60/1.628, dt = 0.4, trunc 7)'

  h = gf_hdur_gaussian(HCMT)
  k = gf_stf_khalf(h,DTG,7.d0)
  allocate(w(-k:k))
  call gf_stf_kernel_gauss(h,DTG,k,w)

  s0 = 0.d0
  s2 = 0.d0
  do j = -k,k
    s0 = s0 + w(j)
    s2 = s2 + w(j)*(dble(j)*DTG)**2
  enddo

  call gf_report('sum of the sampled kernel = 1',abs(s0 - 1.d0),1.d-13,nfail)

  ! Var of g_h is h^2/2, so sqrt(2 Var) recovers h = hdur/1.628: this is the
  ! placement pin testing.md asks for
  call gf_report('sqrt(2 Var) = hdur/1.628 (the 1.628 pin)',abs(sqrt(2.d0*s2/s0) - h)/h,1.d-12,nfail)

  ! the DTFT of the kernel against comp_source_spectrum(om,hdur) =
  ! exp(-(om hdur/1.628)^2/4), specfem's own statement of the spectrum
  do i = 1,3
    om = OMEGA(i)
    dtft = 0.d0
    do j = -k,k
      dtft = dtft + w(j)*cos(om*dble(j)*DTG)
    enddo
    ref = comp_source_spectrum(om,HCMT)
    call gf_report('DTFT vs comp_source_spectrum at omega = '//trim(fmt(om)),abs(dtft - ref),1.d-12,nfail)
  enddo

  deallocate(w)

  end subroutine test_gauss_kernel

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_heavi_kernel(nfail)

! the sampled quasi-Heaviside against specfem's erf, and its symmetry

  implicit none
  integer, intent(inout) :: nfail

  double precision, dimension(:), allocatable :: w
  double precision :: h,err,ref
  integer :: k,j,nasym,nmono

  write(*,'(a)') '3. Heaviside kernel vs netlib_specfun_erf (h = hdur_corr, dt = 0.4, trunc 6)'

  h = sqrt(gf_hdur_gaussian(HCMT)**2 - HDB**2)
  k = gf_stf_khalf(h,DTG,GF_STF_TRUNC)
  allocate(w(-k:k))
  call gf_stf_kernel_heavi(h,DTG,k,w)

  err = 0.d0
  do j = -k,k
    ref = 0.5d0*(1.d0 + netlib_specfun_erf(dble(j)*DTG/h))
    err = max(err,abs(w(j) - ref))
  enddo
  call gf_report('all taps vs 0.5(1+netlib_specfun_erf), 4 ulp',err,4.d0*epsilon(1.d0),nfail)

  call gf_report_true('w(0) == 1/2 exactly',w(0) == 0.5d0,nfail)

  nasym = 0
  nmono = 0
  do j = 1,k
    if (w(j) + w(-j) /= 1.d0) nasym = nasym + 1
    if (w(j) < w(j-1)) nmono = nmono + 1
    if (w(-j) > w(-j+1)) nmono = nmono + 1
  enddo
  call gf_report_true('w(j) + w(-j) == 1 exactly, every j',nasym == 0,nfail)
  call gf_report_true('kernel is monotone',nmono == 0,nfail)
  call gf_report('w(khalf) saturates to 1 (erfc(6)/2)',abs(w(k) - 1.d0),1.d-15,nfail)

  deallocate(w)

  end subroutine test_heavi_kernel

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_delta(nfail)

! convolving a delta returns the kernel, sample for sample: pins centring

  implicit none
  integer, intent(inout) :: nfail

  integer, parameter :: N = 2000, I0 = 1000
  double precision, dimension(:), allocatable :: w,x,y,p
  double precision :: h,err,expect
  integer :: k,i,j,nbad

  write(*,'(a)') '4. delta response'

  h = sqrt(gf_hdur_gaussian(HCMT)**2 - HDB**2)
  k = gf_stf_khalf(h,DTG,GF_STF_TRUNC)
  allocate(w(-k:k),x(N),y(N),p(0:N))

  x(:) = 0.d0
  x(I0) = 1.d0

  call gf_stf_kernel_gauss(h,DTG,k,w)
  call gf_conv_sym(x,N,k,w,y)
  nbad = 0
  do i = 1,N
    j = i - I0
    if (abs(j) <= k) then
      if (y(i) /= w(j)) nbad = nbad + 1
    else
      if (y(i) /= 0.d0) nbad = nbad + 1
    endif
  enddo
  call gf_report_true('Gaussian: y(i0+j) == w(j) exactly, 0 outside',nbad == 0,nfail)

  ! the Heaviside response is dt*w(j) inside the kernel, dt above it and 0
  ! below; a derived tolerance rather than equality, because the product is
  ! rounded here and in gf_conv_heavi by different compilation units (the CI
  ! ifort rounded them differently) -- a centring error would be a whole
  ! sample, not a rounding
  call gf_stf_kernel_heavi(h,DTG,k,w)
  call gf_cumsum(x,N,p)
  call gf_conv_heavi(x,N,DTG,k,w,p,y)
  err = 0.d0
  do i = 1,N
    j = i - I0
    if (j > k) then
      expect = DTG
    else if (j < -k) then
      expect = 0.d0
    else
      expect = DTG*w(j)
    endif
    err = max(err,abs(y(i) - expect)/DTG)
  enddo
  call gf_report('Heaviside: y == dt*w inside, dt above, 0 below (rel)',err,1.d-14,nfail)

  deallocate(w,x,y,p)

  end subroutine test_delta

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_centroid(nfail)

! a symmetric kernel preserves the first moment; an off-by-one shifts it by dt

  implicit none
  integer, intent(inout) :: nfail

  integer, parameter :: N = 3000
  double precision, dimension(:), allocatable :: w,x,y,t
  double precision :: h,cx,cy,sx,sy
  integer :: k,i

  write(*,'(a)') '5. centroid preservation'

  h = sqrt(gf_hdur_gaussian(HCMT)**2 - HDB**2)
  k = gf_stf_khalf(h,DTG,GF_STF_TRUNC)
  allocate(w(-k:k),x(N),y(N),t(N))

  call axis(N,DTG,-600.d0,t)
  ! an asymmetric pulse, well inside the record
  do i = 1,N
    x(i) = gauss(t(i),HDB) + 0.5d0*gauss(t(i) - 30.d0,HDB)
  enddo

  call gf_stf_kernel_gauss(h,DTG,k,w)
  call gf_conv_sym(x,N,k,w,y)

  sx = 0.d0; cx = 0.d0; sy = 0.d0; cy = 0.d0
  do i = 1,N
    sx = sx + x(i); cx = cx + t(i)*x(i)
    sy = sy + y(i); cy = cy + t(i)*y(i)
  enddo
  cx = cx/sx
  cy = cy/sy

  call gf_report('centroid shift / dt',abs(cy - cx)/DTG,1.d-9,nfail)

  deallocate(w,x,y,t)

  end subroutine test_centroid

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_composition(nfail)

! g_h1 * g_h2 = g_h3 on both stored grids, at three truncations

  implicit none
  integer, intent(inout) :: nfail

  double precision, dimension(3), parameter :: TRUNCS = (/ 7.d0, 6.d0, 4.d0 /)
  double precision, dimension(2), parameter :: DTS = (/ DTG, DTR /)
  double precision, dimension(:), allocatable :: w,x,y,t
  double precision :: h1,h2,h3,dt,tspan,err,ref0,e4
  integer :: n,k,i,ig,it

  write(*,'(a)') '6. composition identity g_h1 * g_h2 = g_h3 (h1 = hdur_db, h2 = hdur_corr)'

  h1 = HDB
  h2 = sqrt(gf_hdur_gaussian(HCMT)**2 - HDB**2)
  h3 = sqrt(h1**2 + h2**2)
  ref0 = gauss(0.d0,h3)

  do ig = 1,2
    dt = DTS(ig)
    tspan = 8.d0*h3
    n = 2*nint(tspan/dt) + 1
    allocate(x(n),y(n),t(n))
    call axis(n,dt,-dble(n/2)*dt,t)
    do i = 1,n
      x(i) = gauss(t(i),h1)
    enddo

    do it = 1,3
      k = gf_stf_khalf(h2,dt,TRUNCS(it))
      allocate(w(-k:k))
      call gf_stf_kernel_gauss(h2,dt,k,w)
      call gf_conv_sym(x,n,k,w,y)
      err = 0.d0
      do i = k+1,n-k
        err = max(err,abs(y(i) - gauss(t(i),h3)))
      enddo
      err = err/ref0
      if (TRUNCS(it) >= 6.d0) then
        call gf_report('dt = '//trim(fmt(dt))//', trunc '//trim(fmt(TRUNCS(it)))//': relative error',err,1.d-12,nfail)
      else
        call gf_report('dt = '//trim(fmt(dt))//', trunc '//trim(fmt(TRUNCS(it)))//': relative error',err,1.d-7,nfail)
        e4 = err
        ! the 4h truncation is a real floor (erfc(4)/2 = 8e-9), not round-off
        call gf_report_true('dt = '//trim(fmt(dt))//', trunc 4: error exceeds 1e-9 (truncation floor is real)', &
                            e4 > 1.d-9,nfail)
      endif
      deallocate(w)
    enddo

    deallocate(x,y,t)
  enddo

  end subroutine test_composition

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_trapezoid(nfail)

! trapezoid exact on a linear integrand; its lead over a left Riemann sum

  implicit none
  integer, intent(inout) :: nfail

  integer, parameter :: N = 101
  double precision, parameter :: DT = 0.5d0
  double precision, dimension(N) :: f,g,r,t
  double precision :: err,errr,exact,rs
  integer :: i

  write(*,'(a)') '7. cumulative trapezoid'

  call axis(N,DT,-10.d0,t)
  do i = 1,N
    f(i) = 2.d0 + 3.d0*t(i)
  enddo
  call gf_cumtrapz(f,N,DT,g)

  err = 0.d0
  do i = 1,N
    exact = 2.d0*(t(i) - t(1)) + 1.5d0*(t(i)**2 - t(1)**2)
    err = max(err,abs(g(i) - exact)/max(1.d0,abs(exact)))
  enddo
  call gf_report('exact on a linear integrand, relative',err,1.d-15,nfail)

  ! left Riemann sum r(i) = dt (f_1 + ... + f_{i-1}); the trapezoid leads it
  ! by exactly (dt/2)(f_i - f_1) -- Stage 4's half-sample warning, as a test
  rs = 0.d0
  r(1) = 0.d0
  do i = 2,N
    rs = rs + DT*f(i-1)
    r(i) = rs
  enddo
  errr = 0.d0
  do i = 1,N
    errr = max(errr,abs((g(i) - r(i)) - 0.5d0*DT*(f(i) - f(1)))/max(1.d0,abs(g(i))))
  enddo
  call gf_report('trapezoid - left Riemann == (dt/2)(f_n - f_1)',errr,1.d-15,nfail)

  end subroutine test_trapezoid

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_heavi_limit(nfail)

! hdur_corr = 0: the Heaviside operator is the trapezoid over the record
! extended by zeros, i.e. gf_cumtrapz + dt x_1/2

  implicit none
  integer, intent(inout) :: nfail

  integer, parameter :: N = 300
  double precision, parameter :: DT = 0.5d0
  double precision, dimension(N) :: x,y,g
  double precision, dimension(0:N) :: p
  double precision, dimension(0:0) :: w
  double precision :: err,ref
  integer :: i,nbad

  write(*,'(a)') '8. hdur_corr = 0 limit'

  call gf_stf_kernel_heavi(0.d0,DT,0,w)
  call gf_report_true('khalf = 0 kernel is w(0) = 1/2',w(0) == 0.5d0,nfail)

  ! dyadic data: integers below 2^20 with dt = 1/2, so every operation in
  ! both routes is exact and the comparison can be bitwise
  call gf_lcg_seed(2024)
  do i = 1,N
    x(i) = dble(int(gf_rand()*1048576.d0)) - 524288.d0
  enddo
  call gf_cumsum(x,N,p)
  call gf_conv_heavi(x,N,DT,0,w,p,y)
  call gf_cumtrapz(x,N,DT,g)
  nbad = 0
  do i = 1,N
    if (y(i) /= g(i) + 0.5d0*DT*x(1)) nbad = nbad + 1
  enddo
  call gf_report_true('dyadic data: y == cumtrapz + dt x_1/2, bitwise',nbad == 0,nfail)

  ! random data: the same value up to the rounding order
  do i = 1,N
    x(i) = gf_rand_range(-1.d0,1.d0)
  enddo
  call gf_cumsum(x,N,p)
  call gf_conv_heavi(x,N,DT,0,w,p,y)
  call gf_cumtrapz(x,N,DT,g)
  err = 0.d0
  ref = 0.d0
  do i = 1,N
    err = max(err,abs(y(i) - (g(i) + 0.5d0*DT*x(1))))
    ref = max(ref,abs(g(i)))
  enddo
  call gf_report('random data: y vs cumtrapz + dt x_1/2, relative',err/ref,1.d-14,nfail)

  ! with a zero first sample the boundary term vanishes
  x(1) = 0.d0
  call gf_cumsum(x,N,p)
  call gf_conv_heavi(x,N,DT,0,w,p,y)
  call gf_cumtrapz(x,N,DT,g)
  err = 0.d0
  do i = 1,N
    err = max(err,abs(y(i) - g(i)))
  enddo
  call gf_report('x_1 = 0: y vs cumtrapz, relative',err/ref,1.d-14,nfail)

  end subroutine test_heavi_limit

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_premise_and_routes(nfail)

! The premise: a Gaussian pulse of width h1 convolved with the Heaviside
! kernel of width h2 is the quasi-Heaviside of width h3 -- the operator
! that turns the database's Gaussian response into a CMT response.
!
! The negative control: Stage 4's route (cumulative trapezoid, then a
! Gaussian) on the same input, against the same closed form. Its error is
! the trapezoid integrator's (w dt)^2/12, which in the time domain is
! (dt^2/12) times the second derivative of the answer, g'_h3: predicted
! and measured here side by side.

  implicit none
  integer, intent(inout) :: nfail

  double precision, dimension(2), parameter :: DTS = (/ DTG, DTR /)
  double precision, dimension(:), allocatable :: wh,wg,x,y,ya,g,t,p
  double precision :: h1,h2,h3,dt,tspan,errb,erra,pred
  integer :: n,k,i,ig

  write(*,'(a)') '9. the premise: g_h1 * H_h2 = H_h3, one operator (h1 = hdur_db, h2 = hdur_corr)'

  h1 = HDB
  h2 = sqrt(gf_hdur_gaussian(HCMT)**2 - HDB**2)
  h3 = sqrt(h1**2 + h2**2)

  do ig = 1,2
    dt = DTS(ig)
    tspan = 8.d0*h3
    n = 2*nint(tspan/dt) + 1
    k = gf_stf_khalf(h2,dt,GF_STF_TRUNC)
    allocate(wh(-k:k),wg(-k:k),x(n),y(n),ya(n),g(n),t(n),p(0:n))
    call axis(n,dt,-dble(n/2)*dt,t)
    do i = 1,n
      x(i) = gauss(t(i),h1)
    enddo

    ! route B: the production operator
    call gf_stf_kernel_heavi(h2,dt,k,wh)
    call gf_cumsum(x,n,p)
    call gf_conv_heavi(x,n,dt,k,wh,p,y)
    errb = 0.d0
    do i = 1,n-k
      errb = max(errb,abs(y(i) - heavi(t(i),h3)))
    enddo
    call gf_report('dt = '//trim(fmt(dt))//': one-shot Heaviside vs H_h3 (netlib erf)',errb,1.d-12,nfail)

    ! route A: Stage 4's integrate-then-smooth
    call gf_stf_kernel_gauss(h2,dt,k,wg)
    call gf_cumtrapz(x,n,dt,g)
    call gf_conv_sym(g,n,k,wg,ya)
    erra = 0.d0
    do i = k+1,n-k
      erra = max(erra,abs(ya(i) - heavi(t(i),h3)))
    enddo
    ! (dt^2/12) max|g'_h3|,  max|g'_h| = sqrt(2) exp(-1/2) / (h^2 sqrt(pi))
    pred = (dt**2/12.d0) * sqrt(2.d0)*exp(-0.5d0)/(h3**2*sqrt(PI))
    write(*,'(a,es12.5,a,es12.5)') '       cumtrapz-then-Gaussian: measured error ',erra,'  predicted (dt^2/12) max|g''_h3| ',pred
    call gf_report('dt = '//trim(fmt(dt))//': two-step route error / prediction, within 20 %', &
                   abs(erra/pred - 1.d0),0.2d0,nfail)
    call gf_report_true('dt = '//trim(fmt(dt))//': one-shot beats two-step by > 1e3',errb*1.d3 < erra,nfail)
    if (dt > 1.d0) then
      call gf_report_true('dt = 3.4: two-step error exceeds 1e-6 (the regional integrator error is real)', &
                          erra > 1.d-6,nfail)
    endif

    deallocate(wh,wg,x,y,ya,g,t,p)
  enddo

  end subroutine test_premise_and_routes

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_guard(nfail)

! hdur_target <= hdur_db, and the Gaussian below-one-sample floor

  implicit none
  integer, intent(inout) :: nfail

  integer, parameter :: N = 200
  type(t_gf_stf) :: stf
  double precision, dimension(N) :: x,y,g
  double precision, dimension(0:N) :: p
  double precision, dimension(0:0) :: w
  integer :: ierr,i,nbad

  write(*,'(a)') '10. the guard'

  ! a 5 s CMT on a database of hdur_db = 6.95: narrower than the database
  call gf_stf_plan(GF_SRC_CMT,0,5.d0,HDB,DTR,GF_STF_TRUNC,stf,ierr)
  call gf_report_true('CMT hdur 5 s: plan returns GF_OK',ierr == GF_OK,nfail)
  call gf_report_true('CMT hdur 5 s: guard set',stf%guard,nfail)
  call gf_report_true('CMT hdur 5 s: khalf = 0, hdur_corr = 0',stf%khalf == 0 .and. stf%hdur_corr == 0.d0,nfail)
  call gf_report_true('CMT hdur 5 s: no NaN in the plan', &
                      .not. (gf_is_nan(stf%hdur_corr) .or. gf_is_nan(stf%hdur_target)),nfail)
  call gf_report_true('CMT hdur 5 s: note mentions the guard',index(stf%note,'GUARD') > 0,nfail)
  call gf_report_true('CMT hdur 5 s: kind stays Heaviside',stf%kind_stf == GF_STF_HEAVI,nfail)

  ! dyadic data: the guarded Heaviside is exactly the trapezoid limit
  call gf_lcg_seed(77)
  do i = 1,N
    x(i) = dble(int(gf_rand()*1048576.d0)) - 524288.d0
  enddo
  call gf_stf_kernel(stf,0.5d0,w)
  call gf_cumsum(x,N,p)
  call gf_stf_apply(stf,0.5d0,w,p,x,N,y)
  call gf_cumtrapz(x,N,0.5d0,g)
  nbad = 0
  do i = 1,N
    if (y(i) /= g(i) + 0.25d0*x(1)) nbad = nbad + 1
  enddo
  call gf_report_true('guarded Heaviside == cumtrapz + dt x_1/2, bitwise',nbad == 0,nfail)

  ! a force with f0 = 5: guarded Gaussian is the identity
  call gf_stf_plan(GF_SRC_FORCE,0,5.d0,HDB,DTR,GF_STF_TRUNC,stf,ierr)
  call gf_report_true('force f0 5: guard set, Gaussian',stf%guard .and. stf%kind_stf == GF_STF_GAUSS,nfail)
  call gf_stf_kernel(stf,0.5d0,w)
  call gf_stf_apply(stf,0.5d0,w,p,x,N,y)
  nbad = 0
  do i = 1,N
    if (y(i) /= x(i)) nbad = nbad + 1
  enddo
  call gf_report_true('guarded Gaussian == identity, bitwise',nbad == 0,nfail)

  ! a correction width below one stored sample: h_t = 7.0 against h_db =
  ! 6.95 gives hdur_corr = 0.84 s, under the 3.4 s regional grid. Guarded
  ! for a Gaussian (unresolved kernel), not for a Heaviside (degrades into
  ! the trapezoid continuously).
  call gf_stf_plan(GF_SRC_FORCE,0,7.d0*SOURCE_DECAY_MIMIC_TRIANGLE,HDB,DTR,GF_STF_TRUNC,stf,ierr)
  call gf_report_true('Gaussian with hdur_corr < dt_sub: guarded',stf%guard .and. stf%khalf == 0,nfail)
  call gf_stf_plan(GF_SRC_CMT,0,7.d0*SOURCE_DECAY_MIMIC_TRIANGLE,HDB,DTR,GF_STF_TRUNC,stf,ierr)
  call gf_report_true('Heaviside with hdur_corr < dt_sub: not guarded, khalf >= 1', &
                      (.not. stf%guard) .and. stf%khalf >= 1,nfail)
  call gf_report('  its hdur_corr = sqrt(7^2 - hdur_db^2)',abs(stf%hdur_corr - sqrt(49.d0 - HDB**2)),1.d-14,nfail)

  end subroutine test_guard

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_force_kinds(nfail)

! which conversion each force_stf gets

  implicit none
  integer, intent(inout) :: nfail

  integer, parameter :: N = 50
  type(t_gf_stf) :: stf
  double precision, dimension(:), allocatable :: w
  double precision, dimension(N) :: x,y
  double precision, dimension(0:N) :: p
  integer :: ierr,i,nbad

  write(*,'(a)') '11. force_stf dispatch'

  call gf_stf_plan(GF_SRC_FORCE,0,45.d0,HDB,DTG,GF_STF_TRUNC,stf,ierr)
  call gf_report_true('force_stf 0: Gaussian at hdur/1.628', &
                      stf%kind_stf == GF_STF_GAUSS .and. stf%hdur_target == 45.d0/SOURCE_DECAY_MIMIC_TRIANGLE,nfail)
  call gf_report('  hdur_corr = sqrt((45/1.628)^2 - hdur_db^2)', &
                 abs(stf%hdur_corr - sqrt((45.d0/SOURCE_DECAY_MIMIC_TRIANGLE)**2 - HDB**2)),1.d-14,nfail)
  call gf_report_true('  khalf = ceiling(6 hdur_corr/dt) = 402',stf%khalf == 402,nfail)

  call gf_stf_plan(GF_SRC_FORCE,2,45.d0,HDB,DTG,GF_STF_TRUNC,stf,ierr)
  call gf_report_true('force_stf 2: Heaviside at hdur/1.628', &
                      stf%kind_stf == GF_STF_HEAVI .and. stf%hdur_target == 45.d0/SOURCE_DECAY_MIMIC_TRIANGLE,nfail)

  call gf_stf_plan(GF_SRC_FORCE,4,45.d0,HDB,DTG,GF_STF_TRUNC,stf,ierr)
  call gf_report_true('force_stf 4: Gaussian',stf%kind_stf == GF_STF_GAUSS,nfail)
  call gf_report('  at hdur/pi',abs(stf%hdur_target - 45.d0/PI)/(45.d0/PI),1.d-15,nfail)

  call gf_lcg_seed(5)
  do i = 1,N
    x(i) = gf_rand_range(-1.d0,1.d0)
  enddo
  call gf_cumsum(x,N,p)

  call gf_stf_plan(GF_SRC_FORCE,1,45.d0,HDB,DTG,GF_STF_TRUNC,stf,ierr)
  call gf_report_true('force_stf 1 (Ricker): no conversion, GF_OK',stf%kind_stf == GF_STF_NONE .and. ierr == GF_OK,nfail)
  allocate(w(-stf%khalf:stf%khalf))
  call gf_stf_kernel(stf,DTG,w)
  call gf_stf_apply(stf,DTG,w,p,x,N,y)
  nbad = 0
  do i = 1,N
    if (y(i) /= x(i)) nbad = nbad + 1
  enddo
  call gf_report_true('  pass-through is bitwise',nbad == 0,nfail)
  deallocate(w)

  call gf_stf_plan(GF_SRC_FORCE,3,45.d0,HDB,DTG,GF_STF_TRUNC,stf,ierr)
  call gf_report_true('force_stf 3 (monochromatic): no conversion, GF_OK', &
                      stf%kind_stf == GF_STF_NONE .and. ierr == GF_OK,nfail)

  call gf_stf_plan(GF_SRC_CMT,0,HCMT,HDB,DTG,GF_STF_TRUNC,stf,ierr)
  call gf_report_true('CMT 60 s on the global grid: Heaviside, khalf = 543', &
                      stf%kind_stf == GF_STF_HEAVI .and. stf%khalf == 543,nfail)
  call gf_report('  hdur_corr = 36.19386203437154',abs(stf%hdur_corr - 36.19386203437154d0),1.d-13,nfail)
  call gf_stf_plan(GF_SRC_CMT,0,HCMT,HDB,DTR,GF_STF_TRUNC,stf,ierr)
  call gf_report_true('CMT 60 s on the regional grid: khalf = 64',stf%khalf == 64,nfail)

  end subroutine test_force_kinds

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_taxis(nfail)

! the output axis on both shipped examples

  implicit none
  integer, intent(inout) :: nfail

  type(t_gf_taxis) :: tax
  double precision, dimension(:), allocatable :: t
  double precision :: err
  integer :: ierr,i

  write(*,'(a)') '12. output time axis'

  ! The database axis is (i ss - 1) dt - t0_db (gf_time_axis); the padded
  ! axis must contain it at t(npad+i). A derived tolerance, not equality: the
  ! formula written out here and gf_taxis_times are two compilation units,
  ! which the CI ifort rounded differently. An off-by-one would be a whole
  ! dt_sub.

  ! global, CMT: t0_req = 1.5 * 60
  call gf_taxis_plan(4625,0.1d0,4,T0DB,90.d0,tax,ierr)
  call gf_report_true('global CMT: GF_OK, npad = 139, nt = 4764', &
                      ierr == GF_OK .and. tax%npad == 139 .and. tax%nt == 4764,nfail)
  allocate(t(tax%nt))
  call gf_taxis_times(tax,tax%nt,t)
  call gf_report_true('  t(1) <= -t0_req',t(1) <= -90.d0,nfail)
  call gf_report_true('  npad minimal: t(1) + dt_sub > -t0_req',t(1) + tax%dt_sub > -90.d0,nfail)
  call gf_report_true('  t_first == t(1)',tax%t_first == t(1),nfail)
  call gf_report('  t0 = t0_db + npad dt_sub',abs(tax%t0 - (T0DB + 139*0.4d0)),1.d-12,nfail)
  err = 0.d0
  do i = 1,4625
    err = max(err,abs(t(tax%npad+i) - ((dble(i*4) - 1.d0)*0.1d0 - T0DB)))
  enddo
  call gf_report('  t(npad+i) == database axis (4i - 1) dt - t0_db, all 4625',err,1.d-12,nfail)
  err = 0.d0
  do i = 2,tax%nt
    err = max(err,abs((t(i) - t(i-1)) - 0.4d0))
  enddo
  call gf_report('  uniform spacing',err,1.d-12,nfail)
  call gf_report('  t(1) = -90.0484145764246',abs(t(1) + 90.0484145764246d0),1.d-12,nfail)
  deallocate(t)

  ! regional, CMT
  call gf_taxis_plan(544,0.1d0,34,T0DB,90.d0,tax,ierr)
  call gf_report_true('regional CMT: npad = 18, nt = 562',tax%npad == 18 .and. tax%nt == 562,nfail)
  call gf_report('  t(1) = -92.6484145764246',abs(tax%t_first + 92.6484145764246d0),1.d-12,nfail)

  ! force: t0_req = 1.5 * 45
  call gf_taxis_plan(4625,0.1d0,4,T0DB,67.5d0,tax,ierr)
  call gf_report_true('global force: npad = 83',tax%npad == 83,nfail)
  call gf_taxis_plan(544,0.1d0,34,T0DB,67.5d0,tax,ierr)
  call gf_report_true('regional force: npad = 11',tax%npad == 11,nfail)

  ! no extension asked for: the database axis itself
  call gf_taxis_plan(4625,0.1d0,4,T0DB,0.d0,tax,ierr)
  call gf_report_true('t0_req = 0: npad = 0, nt = nt_db',tax%npad == 0 .and. tax%nt == 4625,nfail)
  allocate(t(tax%nt))
  call gf_taxis_times(tax,tax%nt,t)
  err = 0.d0
  do i = 1,4625
    err = max(err,abs(t(i) - ((dble(i*4) - 1.d0)*0.1d0 - T0DB)))
  enddo
  call gf_report('  == database axis (4i - 1) dt - t0_db',err,1.d-12,nfail)
  call gf_report('  t(1) = (ss-1) dt - t0_db = -34.4484145764246',abs(t(1) + 34.4484145764246d0),1.d-12,nfail)
  deallocate(t)

  ! a request that starts after the database does is never truncated
  call gf_taxis_plan(4625,0.1d0,4,T0DB,-1000.d0,tax,ierr)
  call gf_report_true('t0_req before the database start: npad = 0, nothing removed', &
                      tax%npad == 0 .and. tax%nt == 4625,nfail)

  end subroutine test_taxis

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_padding(nfail)

! left padding is shift-equivariant, bitwise, for both kernels

  implicit none
  integer, intent(inout) :: nfail

  integer, parameter :: N = 500, NPAD = 18
  double precision, dimension(:), allocatable :: w,x,y,xp,yp,p,pp
  double precision :: h,brute,err,ref
  integer :: k,i,j,nbad

  write(*,'(a)') '13. padding'

  h = sqrt(gf_hdur_gaussian(HCMT)**2 - HDB**2)
  k = gf_stf_khalf(h,DTR,GF_STF_TRUNC)
  allocate(w(-k:k),x(N),y(N),xp(NPAD+N),yp(NPAD+N),p(0:N),pp(0:NPAD+N))

  call gf_lcg_seed(31)
  do i = 1,N
    x(i) = gf_rand_range(-1.d0,1.d0)
  enddo
  call gf_pad_left(x,N,NPAD,xp)

  nbad = 0
  do i = 1,NPAD
    if (xp(i) /= 0.d0) nbad = nbad + 1
  enddo
  do i = 1,N
    if (xp(NPAD+i) /= x(i)) nbad = nbad + 1
  enddo
  call gf_report_true('gf_pad_left: leading zeros, then x, bitwise',nbad == 0,nfail)

  ! Gaussian
  call gf_stf_kernel_gauss(h,DTR,k,w)
  call gf_conv_sym(x,N,k,w,y)
  call gf_conv_sym(xp,NPAD+N,k,w,yp)
  nbad = 0
  do i = 1,N
    if (yp(NPAD+i) /= y(i)) nbad = nbad + 1
  enddo
  call gf_report_true('Gaussian: y_pad(npad+i) == y(i), bitwise',nbad == 0,nfail)

  ! the leading outputs are the kernel reaching forward into the record,
  ! not zeros: check them against a plain loop
  err = 0.d0
  ref = maxval(abs(y))
  do i = 1,NPAD
    brute = 0.d0
    do j = -k,k
      if (i-j >= 1 .and. i-j <= NPAD+N) brute = brute + w(j)*xp(i-j)
    enddo
    err = max(err,abs(yp(i) - brute))
  enddo
  call gf_report('Gaussian: leading npad outputs vs plain loop',err/ref,1.d-14,nfail)

  ! Heaviside
  call gf_stf_kernel_heavi(h,DTR,k,w)
  call gf_cumsum(x,N,p)
  call gf_cumsum(xp,NPAD+N,pp)
  call gf_conv_heavi(x,N,DTR,k,w,p,y)
  call gf_conv_heavi(xp,NPAD+N,DTR,k,w,pp,yp)
  nbad = 0
  do i = 1,N
    if (yp(NPAD+i) /= y(i)) nbad = nbad + 1
  enddo
  call gf_report_true('Heaviside: y_pad(npad+i) == y(i), bitwise',nbad == 0,nfail)

  deallocate(w,x,y,xp,yp,p,pp)

  end subroutine test_padding

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_neumaier(nfail)

! the compensated prefix sum keeps what a running sum drops

  implicit none
  integer, intent(inout) :: nfail

  integer, parameter :: N = 100001
  double precision, dimension(:), allocatable :: x,p
  double precision :: exact
  ! `volatile` so that the negative control really is a sequential running
  ! sum: ifort -O3 -xHost otherwise vectorises the loop into several
  ! partial accumulators, most of which start from zero and keep the
  ! 1e-17 increments, and the control stops controlling anything. The
  ! production routine guards its compensation the same way (gf_cumsum).
  double precision, volatile :: naive
  integer :: i

  write(*,'(a)') '14. compensated summation'

  allocate(x(N),p(0:N))
  x(1) = 1.d0
  do i = 2,N
    x(i) = 1.d-17
  enddo
  exact = 1.d0 + dble(N-1)*1.d-17

  call gf_cumsum(x,N,p)

  naive = 0.d0
  do i = 1,N
    naive = naive + x(i)
  enddo
  write(*,'(a,es22.15,a,es22.15,a,es22.15)') '       exact ',exact,'  compensated ',p(N),'  naive ',naive

  call gf_report('1 + 1e5 x 1e-17: compensated vs exact',abs(p(N) - exact),1.d-15,nfail)
  call gf_report_true('  the naive running sum loses it entirely',naive == 1.d0,nfail)
  call gf_report_true('  p(0) == 0',p(0) == 0.d0,nfail)

  deallocate(x,p)

  end subroutine test_neumaier

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_onset(nfail)

! the silence-before-the-record check

  implicit none
  integer, intent(inout) :: nfail

  integer, parameter :: N = 544
  double precision, dimension(N) :: x,t
  double precision :: ratio
  integer :: i,nbefore,nexp

  write(*,'(a)') '15. onset guard'

  ! the regional axis: two samples lie before -4 hdur_db = -27.8 s
  do i = 1,N
    t(i) = (dble(i*34) - 1.d0)*0.1d0 - T0DB
  enddo
  nexp = 0
  do i = 1,N
    if (t(i) < -4.d0*HDB) nexp = nexp + 1
  enddo

  x(:) = 0.d0
  do i = 1,N
    if (t(i) > 100.d0) x(i) = sin(0.05d0*t(i))
  enddo
  ! pin the maximum at exactly 1 so the ratio below is an exact division
  x(N) = 1.d0
  call gf_stf_onset(x,N,t,HDB,ratio,nbefore)
  call gf_report_true('silent before -4 hdur_db: ratio 0, nbefore = '//trim(fmti(nexp)), &
                      ratio == 0.d0 .and. nbefore == nexp,nfail)

  x(1) = 0.01d0
  call gf_stf_onset(x,N,t,HDB,ratio,nbefore)
  call gf_report('a 1 % sample before the onset: ratio 0.01',abs(ratio - 0.01d0),1.d-15,nfail)

  x(:) = 0.d0
  call gf_stf_onset(x,N,t,HDB,ratio,nbefore)
  call gf_report_true('all-zero trace: ratio 0, no division',ratio == 0.d0,nfail)

  end subroutine test_onset

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_errors(nfail)

! bad arguments come back as GF_ERR_ARG, never a stop or a trap

  implicit none
  integer, intent(inout) :: nfail

  type(t_gf_stf) :: stf
  type(t_gf_taxis) :: tax
  integer :: ierr

  write(*,'(a)') '16. error paths'

  call gf_stf_plan(GF_SRC_CMT,0,HCMT,HDB,0.d0,GF_STF_TRUNC,stf,ierr)
  call gf_report_true('dt_sub = 0 -> GF_ERR_ARG',ierr == GF_ERR_ARG,nfail)
  call gf_stf_plan(GF_SRC_CMT,0,HCMT,0.d0,DTG,GF_STF_TRUNC,stf,ierr)
  call gf_report_true('hdur_db = 0 -> GF_ERR_ARG',ierr == GF_ERR_ARG,nfail)
  call gf_stf_plan(0,0,HCMT,HDB,DTG,GF_STF_TRUNC,stf,ierr)
  call gf_report_true('source type 0 -> GF_ERR_ARG',ierr == GF_ERR_ARG,nfail)
  call gf_stf_plan(GF_SRC_FORCE,7,HCMT,HDB,DTG,GF_STF_TRUNC,stf,ierr)
  call gf_report_true('force_stf 7 -> GF_ERR_ARG',ierr == GF_ERR_ARG,nfail)
  call gf_stf_plan(GF_SRC_CMT,0,HCMT,HDB,DTG,0.d0,stf,ierr)
  call gf_report_true('trunc 0 -> GF_ERR_ARG',ierr == GF_ERR_ARG,nfail)

  call gf_taxis_plan(0,0.1d0,4,T0DB,90.d0,tax,ierr)
  call gf_report_true('nt_db = 0 -> GF_ERR_ARG',ierr == GF_ERR_ARG,nfail)
  call gf_taxis_plan(4625,0.1d0,0,T0DB,90.d0,tax,ierr)
  call gf_report_true('subsample_step = 0 -> GF_ERR_ARG',ierr == GF_ERR_ARG,nfail)
  call gf_taxis_plan(4625,0.d0,4,T0DB,90.d0,tax,ierr)
  call gf_report_true('dt = 0 -> GF_ERR_ARG',ierr == GF_ERR_ARG,nfail)

  end subroutine test_errors

!
!-------------------------------------------------------------------------------------------------
!

  function fmt(v) result(str)

  implicit none
  double precision, intent(in) :: v
  character(len=16) :: str

  write(str,'(f0.2)') v

  end function fmt

!
!-------------------------------------------------------------------------------------------------
!

  function fmti(v) result(str)

  implicit none
  integer, intent(in) :: v
  character(len=16) :: str

  write(str,'(i0)') v

  end function fmti

  end program test_gf_stf
