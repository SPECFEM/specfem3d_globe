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
!---- Source time function operators and the output time axis.
!----
!---- What is stored, and what is wanted
!---- ---------------------------------
!---- The reciprocal runs that filled the database used a Gaussian source
!---- time function of width hdur_db,
!----
!----   g_h(t) = exp(-(t/h)^2) / (h sqrt(pi))              (comp_source_time_function.f90:202)
!----
!---- lowpassed by a zero-phase Butterworth (green_function_stf.F90:82), so
!---- what the database holds is  G * butter(g_hdur_db)  for the Green
!---- function G. hdur_db is *already* a Gaussian width: the writer sets it
!---- to T_min/10 and uses it undivided.
!----
!---- A forward run with a CMTSOLUTION uses specfem's quasi-Heaviside
!----
!----   H_h(t) = 0.5 (1 + erf(t/h)),   h = hdur/SOURCE_DECAY_MIMIC_TRIANGLE
!----                                                     (comp_source_time_function.f90:80,
!----                                                      setup_sources_receivers.f90:777)
!----
!---- and a forward run with a FORCESOLUTION of force_stf = 0 uses g_h with
!---- the same h. So 1.628 divides the *source file's* field, once, and never
!---- touches the database's width. That is the whole of the convention
!---- question; the two older convolvers in the tree
!---- (src/auxiliaries/convolve_source_timefunction.f90:66 and
!---- utils/scripts/lib/convolve_stf.c:210) disagree with each other about
!---- it, and the solver is the authority.
!----
!---- One convolution, not an integration and a convolution
!---- ------------------------------------------------------
!---- Gaussians compose:  g_a * g_b = g_c  with  c^2 = a^2 + b^2.  And H_h is
!---- the integral of g_h. So with  hdur_corr^2 = hdur_target^2 - hdur_db^2,
!----
!----   (G * g_hdur_db) * H_hdur_corr  =  G * INTEGRAL(g_hdur_db * g_hdur_corr)
!----                                  =  G * H_hdur_target
!----
!---- which is the forward CMT seismogram (lowpassed by the database's own
!---- Butterworth). The stored trace is convolved with specfem's own source
!---- time function at the corrected width, and nothing is integrated.
!----
!---- Stage 4 did it in two steps -- cumulative trapezoid, then a Gaussian
!---- -- and the trapezoid integrator is where the error was. Its transfer
!---- function is (1/i w)(x cot x) with x = w dt/2, a broadband amplitude
!---- error of -(w dt)^2/12; in the time domain that is (dt^2/12) times the
!---- second derivative of the answer, which for a seismogram of dominant
!---- period T is (2 pi/T)^2 dt^2/12 of its amplitude: one percent at T = 60 s
!---- on the 3.4 s grid of the shipped regional example. The one-step form
!---- has no integrator; its only error is the aliasing of the sampled
!---- kernel, exp(-(pi hdur_corr/(2 dt))^2), which is 1e-121 on that grid.
!----
!---- The kernel and its cost
!---- -----------------------
!---- H saturates to 0 below -K dt and to 1 above +K dt, so with P the prefix
!---- sum of the trace,
!----
!----   y(i) = dt * ( P(i-K-1) + SUM_{|j| <= K} H(j dt) x(i-j) )
!----
!---- is O(N K) with K = trunc*hdur_corr/dt, a few million flops per trace on
!---- the shipped examples. No FFT. The kernel is built exactly symmetric --
!---- w(0) = 1/2 and w(-j) = 1 - w(j), which is exact in floating point since
!---- w(j) lies in [1/2, 1] -- so the operator is zero-phase by construction,
!---- and no timing error can originate here. What can still be wrong in
!---- time is t0 and the sample indexing, which is what the forward-run
!---- comparison's cross-correlation lag is for.
!----
!---- hdur_corr -> 0 makes K = 0 and the sum dt (x_1 + ... + x_{i-1} + x_i/2):
!---- the cumulative trapezoid over the record extended by zeros to the left.
!---- That is the guard case (hdur_target <= hdur_db), where the requested
!---- source is narrower than the database's and no width correction exists.
!---- gf_cumtrapz is kept, unchanged, for --dump and as the test oracle.
!----
!---- Why the order strain -> contract -> STF is required
!---- --------------------------------------------------
!---- A force pulse has unit impulse, so the reciprocal run leaves the planet
!---- with permanent linear momentum: the stored displacement carries a rigid
!---- translation growing linearly in time. Its strain is identically zero.
!---- The STF operator therefore acts on the contracted strain trace, after
!---- the spatial derivative has removed the drift, and never on raw
!---- displacement. For the same reason a step force cannot be used to build
!---- the database (permanent force => permanent acceleration, u ~ t^2), and
!---- integrating in the solver at write time would store that drift.
!----
!---- The mean of the strain trace is physics, not a nuisance: its integral
!---- over the record is the static strain under a step force, which by
!---- reciprocity is the static offset of the CMT response. Nothing here
!---- detrends.
!----
!---- Extension, and the time axis
!---- ----------------------------
!---- The trace is extended by zeros on the left, which is physical -- the
!---- reciprocal source has not switched on before the record starts, by the
!---- writer's choice of t0 = T_min/2 -- and by zeros on the right, which is
!---- not: the last K samples use data the record does not have, and are
!---- reported as such. gf_stf_onset measures how well the left assumption
!---- holds.
!----
!---- The output axis is the database's own, extended to the left by whole
!---- samples so that a requested start time is covered:
!----
!----   t(i) = ((i - npad)*subsample_step - 1)*dt - t0_db
!----
!---- With npad = 0 this is exactly the Stage 3 axis (verified against the
!---- writer), and t(npad+j) is bitwise the database's t(j). Note t(1) is
!---- (subsample_step-1)*dt - t0, not -t0: the writer snapshots at solver step
!---- i*subsample_step, so the first stored sample sits (ss-1)*dt after the
!---- start of the run.
!----
!---- No `use hdf5`, no `use specfem_par`: this is a kernel module.
!----

  module gf_stf

  use constants, only: PI,SOURCE_DECAY_MIMIC_TRIANGLE

  use gf_par, only: t_gf_stf,t_gf_taxis,gf_set_error, &
                    GF_OK,GF_ERR_ARG, &
                    GF_STF_NONE,GF_STF_GAUSS,GF_STF_HEAVI,GF_STF_TRUNC, &
                    GF_SRC_FORCE,GF_SRC_CMT

  implicit none

  private

  public :: gf_hdur_gaussian
  public :: gf_stf_khalf
  public :: gf_stf_plan
  public :: gf_stf_kernel_gauss
  public :: gf_stf_kernel_heavi
  public :: gf_stf_kernel_gauss_unit
  public :: gf_stf_kernel
  public :: gf_cumsum
  public :: gf_conv_sym
  public :: gf_conv_heavi
  public :: gf_stf_apply
  public :: gf_taxis_plan
  public :: gf_taxis_times
  public :: gf_pad_left
  public :: gf_stf_onset
  public :: gf_stf_kind_name
  public :: gf_print_stf
  public :: gf_cumtrapz

  ! a kernel longer than this is a mistake in the arguments, not a request
  integer, parameter :: GF_STF_KHALF_MAX = 100000000

  contains

!
!-------------------------------------------------------------------------------------------------
!

  double precision function gf_hdur_gaussian(hdur)

! the Gaussian width specfem uses for a source file's half duration:
! hdur_Gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE (setup_sources_receivers.f90:777)

  implicit none

  double precision, intent(in) :: hdur

  gf_hdur_gaussian = hdur / SOURCE_DECAY_MIMIC_TRIANGLE

  end function gf_hdur_gaussian

!
!-------------------------------------------------------------------------------------------------
!

  integer function gf_stf_khalf(hdur,dt,trunc)

! kernel half length in samples for a width hdur truncated at trunc*hdur

  implicit none

  double precision, intent(in) :: hdur,dt,trunc

  if (hdur <= 0.d0 .or. dt <= 0.d0 .or. trunc <= 0.d0) then
    gf_stf_khalf = 0
  else
    gf_stf_khalf = ceiling(trunc*hdur/dt)
  endif

  end function gf_stf_khalf

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_stf_plan(source_type,force_stf,hdur_src,hdur_db,dt_sub,trunc,stf,ierr)

! decides which conversion a source needs, and its width
!
! Pure arithmetic on scalars: no database, so tests/gf3d/test_gf_stf pins
! the 1.628 placement without one.
!
! The guard: when hdur_target <= hdur_db the requested source is narrower
! than the database's, the correction width is imaginary, and the trace is
! left at the database's own width -- integrated (HEAVI) or as is (GAUSS)
! -- with the fact recorded in `note`. It is checked *before* the sqrt,
! because the tree is built with -ffpe-trap and an invalid sqrt is a trap,
! not a NaN.
!
! A Gaussian kernel narrower than one stored sample is not resolved by the
! grid (its sampled sum is not 1), so for GF_STF_GAUSS a correction width
! below dt_sub is also skipped, and said so. The Heaviside kernel has no
! such floor: it degrades continuously into the trapezoid.

  implicit none

  integer, intent(in) :: source_type,force_stf
  double precision, intent(in) :: hdur_src,hdur_db,dt_sub,trunc
  type(t_gf_stf), intent(out) :: stf
  integer, intent(out) :: ierr

  ! local parameters
  double precision :: hsq
  character(len=64) :: tmp

  ierr = GF_OK

  if (dt_sub <= 0.d0) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_stf_plan: dt_sub must be positive')
    return
  endif
  if (hdur_db <= 0.d0) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_stf_plan: the database half duration must be positive')
    return
  endif
  if (trunc <= 0.d0) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_stf_plan: the kernel truncation must be positive')
    return
  endif

  stf%hdur_src = hdur_src
  stf%hdur_db = hdur_db
  stf%trunc = trunc

  select case (source_type)

  case (GF_SRC_CMT)
    ! comp_source_time_function_heavi at hdur_Gaussian (compute_add_sources.f90:556)
    stf%kind_stf = GF_STF_HEAVI
    stf%hdur_target = gf_hdur_gaussian(hdur_src)
    stf%note = 'CMT: quasi-Heaviside 0.5(1+erf(t/h)) at h = hdur/1.628'

  case (GF_SRC_FORCE)
    ! the force_stf branches of compute_add_sources.f90:528-547
    select case (force_stf)
    case (0)
      stf%kind_stf = GF_STF_GAUSS
      stf%hdur_target = gf_hdur_gaussian(hdur_src)
      stf%note = 'force_stf 0: Gaussian at h = hdur/1.628'
    case (2)
      stf%kind_stf = GF_STF_HEAVI
      stf%hdur_target = gf_hdur_gaussian(hdur_src)
      stf%note = 'force_stf 2: quasi-Heaviside 0.5(1+erf(t/h)) at h = hdur/1.628'
    case (4)
      ! comp_source_time_function_gauss_2, sqrt(pi/h^2) exp(-pi^2 t^2/h^2),
      ! is the unit-area Gaussian of width h/pi
      stf%kind_stf = GF_STF_GAUSS
      stf%hdur_target = hdur_src / PI
      stf%note = 'force_stf 4: Meschede Gaussian, a unit Gaussian at h = hdur/pi'
    case (1)
      stf%kind_stf = GF_STF_NONE
      stf%note = 'force_stf 1 (Ricker) is not in the Gaussian family: no conversion, ' // &
                 'the trace is the response to the database Gaussian'
    case (3)
      stf%kind_stf = GF_STF_NONE
      stf%note = 'force_stf 3 (monochromatic) is not in the Gaussian family: no conversion, ' // &
                 'the trace is the response to the database Gaussian'
    case default
      write(tmp,'(i0)') force_stf
      call gf_set_error(ierr,GF_ERR_ARG,'gf_stf_plan: unsupported force_stf = '//trim(tmp))
      return
    end select

  case default
    call gf_set_error(ierr,GF_ERR_ARG,'gf_stf_plan: the source has no type set')
    return

  end select

  if (stf%kind_stf == GF_STF_NONE) return

  if (stf%hdur_target <= stf%hdur_db) then
    stf%guard = .true.
    stf%hdur_corr = 0.d0
    stf%khalf = 0
    stf%note = trim(stf%note)//'; GUARD: requested width <= database width, no width correction ' // &
               '(the trace keeps the database Gaussian; for a Heaviside it is the plain integral, ' // &
               'which on this grid carries the trapezoid error (w dt)^2/12)'
    return
  endif

  hsq = stf%hdur_target**2 - stf%hdur_db**2
  stf%hdur_corr = sqrt(hsq)

  if (stf%kind_stf == GF_STF_GAUSS .and. stf%hdur_corr < dt_sub) then
    stf%guard = .true.
    stf%hdur_corr = 0.d0
    stf%khalf = 0
    stf%note = trim(stf%note)//'; GUARD: correction width below one stored sample, not resolved ' // &
               'by the grid, no width correction'
    return
  endif

  if (trunc*stf%hdur_corr/dt_sub > dble(GF_STF_KHALF_MAX)) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_stf_plan: the kernel would exceed 1e8 samples; check dt and hdur')
    return
  endif

  stf%khalf = gf_stf_khalf(stf%hdur_corr,dt_sub,trunc)

  end subroutine gf_stf_plan

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_stf_kernel_gauss(hdur,dt,khalf,w)

! w(j) = g_hdur(j dt) dt, the sampled unit-area Gaussian, mirrored exactly
!
! Not renormalised: the sum is 1 up to the aliasing of the sampled Gaussian,
! and leaving it so is what lets the test measure that aliasing. A zero half
! length is the identity.

  implicit none

  integer, intent(in) :: khalf
  double precision, intent(in) :: hdur,dt
  double precision, dimension(-khalf:khalf), intent(out) :: w

  ! local parameters
  integer :: j
  double precision :: a,t

  if (khalf == 0 .or. hdur <= 0.d0) then
    w(0) = 1.d0
    return
  endif

  a = dt / (hdur*sqrt(PI))

  do j = 0,khalf
    t = dble(j)*dt
    w(j) = a*exp(-(t/hdur)**2)
  enddo
  do j = 1,khalf
    w(-j) = w(j)
  enddo

  end subroutine gf_stf_kernel_gauss

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_stf_kernel_heavi(hdur,dt,khalf,w)

! w(j) = 0.5 (1 + erf(j dt/hdur)), specfem's comp_source_time_function_heavi
! sampled on the stored grid, with the F2008 erf intrinsic
!
! Built for j > 0 and reflected as w(-j) = 1 - w(j): with w(j) in [1/2, 1]
! that subtraction is exact (Sterbenz), so w(j) + w(-j) == 1 holds bitwise
! and the kernel is symmetric about j = 0 by construction, not by luck. A
! zero width is the unit step with w(0) = 1/2, i.e. the trapezoid rule.

  implicit none

  integer, intent(in) :: khalf
  double precision, intent(in) :: hdur,dt
  double precision, dimension(-khalf:khalf), intent(out) :: w

  ! local parameters
  integer :: j

  w(0) = 0.5d0
  if (khalf == 0) return

  do j = 1,khalf
    if (hdur > 0.d0) then
      w(j) = 0.5d0 + 0.5d0*erf(dble(j)*dt/hdur)
    else
      w(j) = 1.d0
    endif
    w(-j) = 1.d0 - w(j)
  enddo

  end subroutine gf_stf_kernel_heavi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_stf_kernel_gauss_unit(hdur,dt,khalf,w,wsum_raw)

! gf_stf_kernel_gauss divided by its own sum: the sampled unit Gaussian,
! normalised so that the discrete kernel has unit sum as well
!
! It exists for the centroid-time partial of gf_partials.F90, which needs
! the derivative of the quasi-Heaviside kernel -- and d/dt [0.5 (1 +
! erf(t/h))] = exp(-(t/h)^2)/(h sqrt(pi)) is the unit Gaussian, nothing
! else. What this routine adds is only the normalisation, and here is why:
!
! The sampled Gaussian sums to 1 only up to aliasing. By Poisson summation
!
!   SUM_j dt g_h(j dt) = 1 + 2 SUM_k >= 1 exp(-(pi k h / dt)^2)
!
! which is 7e-18 above 1 at h = 2 dt, 1e-4 at h = dt, and 2.3 at h = dt/4.
! Dividing by the sum is therefore a no-op whenever the kernel is resolved
! -- there the sampled Gaussian is spectrally exact, and a box-averaged
! kernel would attenuate by (w dt)^2/24, half a percent at 60 s on the
! 3.4 s grid -- so above h = 2 dt no division is made at all, and below it
! the division is what keeps the static limit right when the kernel
! is not resolved. That regime exists: gf_stf_plan does not guard the
! Heaviside conversion below one sample (it degrades continuously into the
! trapezoid), so hdur_corr in (0, dt) is reached by a CMT half duration
! between 11.3 and 12.6 s on the shipped grids. The normalised kernel is
! continuous in h down to h = 0, where it is the identity w(0) = 1: the
! derivative of the trapezoid limit, i.e. the integrand itself.
!
! `wsum_raw` is the sum before normalisation, the aliasing measure, so a
! caller can report it.

  implicit none

  integer, intent(in) :: khalf
  double precision, intent(in) :: hdur,dt
  double precision, dimension(-khalf:khalf), intent(out) :: w
  double precision, intent(out) :: wsum_raw

  ! local parameters
  integer :: j
  double precision :: s

  call gf_stf_kernel_gauss(hdur,dt,khalf,w)

  wsum_raw = 1.d0
  if (khalf == 0) return
  if (hdur <= 0.d0) then
    ! the identity, spelled out: gf_stf_kernel_gauss only sets w(0) here
    w(:) = 0.d0
    w(0) = 1.d0
    return
  endif

  ! positive terms in increasing order of magnitude, so a plain running sum
  ! is accurate to khalf ulp
  s = 0.d0
  do j = khalf,1,-1
    s = s + w(j)
  enddo
  s = 2.d0*s + w(0)
  wsum_raw = s

  ! Resolved: at h >= 2 dt the aliasing tail 2 exp(-(pi h/dt)^2) is below
  ! 7e-18, less than half an ulp of 1, so the exact sum *is* 1 and the
  ! sampled Gaussian is returned untouched. The decision is made from the
  ! closed form and not from `s`, because `s` itself carries the rounding
  ! of a few hundred additions and can sit an ulp off 1 when the true sum
  ! does not; dividing by it would perturb a spectrally exact kernel for
  ! nothing. tests/gf3d/test_gf_partials.f90 asserts the identity bitwise.
  if (hdur >= 2.d0*dt) return

  ! the same division for w(j) and w(-j), so the symmetry survives bitwise
  do j = -khalf,khalf
    w(j) = w(j)/s
  enddo

  end subroutine gf_stf_kernel_gauss_unit

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_stf_kernel(stf,dt,w)

! the kernel a plan asks for; w(-stf%khalf:stf%khalf)

  implicit none

  type(t_gf_stf), intent(in) :: stf
  double precision, intent(in) :: dt
  double precision, dimension(-stf%khalf:stf%khalf), intent(out) :: w

  select case (stf%kind_stf)
  case (GF_STF_GAUSS)
    call gf_stf_kernel_gauss(stf%hdur_corr,dt,stf%khalf,w)
  case (GF_STF_HEAVI)
    call gf_stf_kernel_heavi(stf%hdur_corr,dt,stf%khalf,w)
  case default
    w(0) = 1.d0
  end select

  end subroutine gf_stf_kernel

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_cumsum(x,n,p)

! prefix sums, p(i) = x(1) + ... + x(i), with p(0) = 0
!
! Neumaier's compensated summation: the running error term is carried
! separately, so the result does not degrade with n. A plain running sum
! over an oscillating trace loses digits in proportion to the sum of the
! magnitudes, and the quantity this feeds -- the saturated tail of the
! Heaviside kernel -- is exactly the static offset that the oscillations
! nearly cancel to.
!
! The two-sum is written as single-operation statements on `volatile`
! temporaries, not as the textbook one-liner c = c + ((s - t) + x). The
! one-liner is only correct if the compiler honours the parentheses, and
! ifort/ifx at their default -fp-model fast do not: they reassociate
! (s - t) + x into (s + x) - t, which is t - t = 0, and the correction
! silently vanishes -- tests/gf3d/test_gf_stf.f90 (test_neumaier) caught
! exactly that under ifort. A volatile temporary must be stored and
! re-read at every reference, so no two of these operations can be fused
! or reordered, whatever the floating-point model. The cost is three
! memory round trips per sample, on arrays of a few thousand samples.

  implicit none

  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: x
  double precision, dimension(0:n), intent(out) :: p

  ! local parameters
  integer :: i
  double precision :: s,c
  double precision, volatile :: t,e1,e2

  p(0) = 0.d0
  s = 0.d0
  c = 0.d0
  do i = 1,n
    t = s + x(i)
    if (abs(s) >= abs(x(i))) then
      e1 = s - t
      e2 = e1 + x(i)
    else
      e1 = x(i) - t
      e2 = e1 + s
    endif
    c = c + e2
    s = t
    p(i) = s + c
  enddo

  end subroutine gf_cumsum

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_conv_sym(x,n,khalf,w,y)

! y(i) = SUM_{j=-khalf..khalf} w(j) x(i-j), with x taken as zero outside 1..n
!
! Taps that would read before the record are the left extension (silence,
! physical); taps that would read past it are the right extension (unknown
! data), which affects the last khalf outputs.

  implicit none

  integer, intent(in) :: n,khalf
  double precision, dimension(n), intent(in) :: x
  double precision, dimension(-khalf:khalf), intent(in) :: w
  double precision, dimension(n), intent(out) :: y

  ! local parameters
  integer :: i,j,jlo,jhi
  double precision :: s,c
  ! the compensated two-sum on volatile temporaries, for the reason given
  ! at gf_cumsum; `v` is volatile as well so that the product cannot be
  ! fused with the following addition into an FMA, which would make t
  ! something other than the rounded sum the correction is derived from
  double precision, volatile :: t,v,e1,e2

  do i = 1,n
    jlo = max(-khalf,i-n)
    jhi = min(khalf,i-1)
    s = 0.d0
    c = 0.d0
    do j = jlo,jhi
      v = w(j)*x(i-j)
      t = s + v
      if (abs(s) >= abs(v)) then
        e1 = s - t
        e2 = e1 + v
      else
        e1 = v - t
        e2 = e1 + s
      endif
      c = c + e2
      s = t
    enddo
    y(i) = s + c
  enddo

  end subroutine gf_conv_sym

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_conv_heavi(x,n,dt,khalf,w,p,y)

! y(i) = dt ( p(i-khalf-1) + SUM_{j=-khalf..khalf} w(j) x(i-j) )
!
! the convolution with a kernel that is w inside |j| <= khalf, 1 above and
! 0 below -- the quasi-Heaviside -- with p the prefix sums of x from
! gf_cumsum. Every term is the trapezoid rule for the convolution integral
! over the record extended by zeros; with khalf = 0 and w(0) = 1/2 the
! result is dt (x_1 + ... + x_{i-1} + x_i/2), which is gf_cumtrapz plus the
! segment dt x_1/2 that a zero-extended record has before its first sample.

  implicit none

  integer, intent(in) :: n,khalf
  double precision, dimension(n), intent(in) :: x
  double precision, intent(in) :: dt
  double precision, dimension(-khalf:khalf), intent(in) :: w
  double precision, dimension(0:n), intent(in) :: p
  double precision, dimension(n), intent(out) :: y

  ! local parameters
  integer :: i,j,jlo,jhi,m
  double precision :: s,c
  ! volatile two-sum temporaries, as in gf_conv_sym
  double precision, volatile :: t,v,e1,e2

  do i = 1,n
    jlo = max(-khalf,i-n)
    jhi = min(khalf,i-1)
    s = 0.d0
    c = 0.d0
    do j = jlo,jhi
      v = w(j)*x(i-j)
      t = s + v
      if (abs(s) >= abs(v)) then
        e1 = s - t
        e2 = e1 + v
      else
        e1 = v - t
        e2 = e1 + s
      endif
      c = c + e2
      s = t
    enddo
    ! the samples the kernel has already saturated over
    m = max(i-khalf-1,0)
    y(i) = dt*(p(m) + (s + c))
  enddo

  end subroutine gf_conv_heavi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_stf_apply(stf,dt,w,p,x,n,y)

! applies a plan to one trace: x and y are on the same (already padded) axis
!
! `p` is only read for GF_STF_HEAVI but must be a valid (0:n) array; the
! caller computes it once per trace with gf_cumsum.

  implicit none

  integer, intent(in) :: n
  type(t_gf_stf), intent(in) :: stf
  double precision, intent(in) :: dt
  double precision, dimension(-stf%khalf:stf%khalf), intent(in) :: w
  double precision, dimension(0:n), intent(in) :: p
  double precision, dimension(n), intent(in) :: x
  double precision, dimension(n), intent(out) :: y

  select case (stf%kind_stf)
  case (GF_STF_GAUSS)
    call gf_conv_sym(x,n,stf%khalf,w,y)
  case (GF_STF_HEAVI)
    call gf_conv_heavi(x,n,dt,stf%khalf,w,p,y)
  case default
    y(:) = x(:)
  end select

  end subroutine gf_stf_apply

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_taxis_plan(nt_db,dt,subsample_step,t0_db,t0_req,tax,ierr)

! the output axis for a requested start time
!
! npad is the smallest number of stored samples that, prepended, puts t(1)
! at or before -t0_req. Nothing is ever removed: a request that starts
! after the database does gives npad = 0 and the database's own axis, and
! tax%t_first says so.

  implicit none

  integer, intent(in) :: nt_db,subsample_step
  double precision, intent(in) :: dt,t0_db,t0_req
  type(t_gf_taxis), intent(out) :: tax
  integer, intent(out) :: ierr

  ! local parameters
  double precision :: t_db1,v
  double precision, dimension(1) :: t1

  ierr = GF_OK

  if (nt_db < 1) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_taxis_plan: the database has no time samples')
    return
  endif
  if (subsample_step < 1) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_taxis_plan: subsample_step must be at least 1')
    return
  endif
  if (dt <= 0.d0) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_taxis_plan: dt must be positive')
    return
  endif

  tax%nt_db = nt_db
  tax%subsample_step = subsample_step
  tax%dt = dt
  tax%dt_sub = dt*dble(subsample_step)
  tax%t0_db = t0_db
  tax%t0_req = t0_req

  ! the database's first sample time, in the writer's convention
  t_db1 = (dble(subsample_step) - 1.d0)*dt - t0_db

  v = (t_db1 + t0_req)/tax%dt_sub
  if (v > 0.d0) then
    if (v > dble(GF_STF_KHALF_MAX)) then
      call gf_set_error(ierr,GF_ERR_ARG,'gf_taxis_plan: the requested t0 would prepend more than 1e8 samples')
      return
    endif
    tax%npad = ceiling(v)
  else
    tax%npad = 0
  endif

  if (huge(1)/subsample_step < nt_db + tax%npad) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_taxis_plan: nt*subsample_step overflows the default integer')
    return
  endif

  tax%nt = nt_db + tax%npad
  tax%t0 = t0_db + dble(tax%npad)*tax%dt_sub

  call gf_taxis_times(tax,1,t1)
  tax%t_first = t1(1)

  end subroutine gf_taxis_plan

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_taxis_times(tax,nt,t)

! t(i) = ((i - npad)*subsample_step - 1)*dt - t0_db
!
! For i > npad this is, expression for expression, the Stage 3 axis
! (gf_seismograms: gf_time_axis), so the database's own sample times come
! out bitwise.

  implicit none

  type(t_gf_taxis), intent(in) :: tax
  integer, intent(in) :: nt
  double precision, dimension(nt), intent(out) :: t

  ! local parameters
  integer :: i

  do i = 1,nt
    t(i) = (dble((i - tax%npad)*tax%subsample_step) - 1.d0)*tax%dt - tax%t0_db
  enddo

  end subroutine gf_taxis_times

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_pad_left(x,n,npad,y)

! y = (npad zeros, x)

  implicit none

  integer, intent(in) :: n,npad
  double precision, dimension(n), intent(in) :: x
  double precision, dimension(npad+n), intent(out) :: y

  ! local parameters
  integer :: i

  do i = 1,npad
    y(i) = 0.d0
  enddo
  do i = 1,n
    y(npad+i) = x(i)
  enddo

  end subroutine gf_pad_left

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_stf_onset(x,n,t,hdur_db,ratio,nbefore)

! how silent the trace is before the reciprocal source switches on
!
! ratio = max |x| over the samples with t < -4 hdur_db, divided by max |x|
! over the trace. The conversion extends the record with zeros to the left,
! which is right only if the record really is zero there; the writer's
! t0 = T_min/2 is meant to guarantee it, and this is the check that it did.
! Zero when there are no such samples, or the trace is identically zero --
! never a division.

  implicit none

  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: x,t
  double precision, intent(in) :: hdur_db
  double precision, intent(out) :: ratio
  integer, intent(out) :: nbefore

  ! local parameters
  integer :: i
  double precision :: xmax,xpre,tcut

  tcut = -4.d0*hdur_db

  xmax = 0.d0
  xpre = 0.d0
  nbefore = 0
  do i = 1,n
    xmax = max(xmax,abs(x(i)))
    if (t(i) < tcut) then
      nbefore = nbefore + 1
      xpre = max(xpre,abs(x(i)))
    endif
  enddo

  if (nbefore == 0 .or. xmax == 0.d0) then
    ratio = 0.d0
  else
    ratio = xpre/xmax
  endif

  end subroutine gf_stf_onset

!
!-------------------------------------------------------------------------------------------------
!

  function gf_stf_kind_name(kind_stf) result(str)

  implicit none

  integer, intent(in) :: kind_stf
  character(len=9) :: str

  select case (kind_stf)
  case (GF_STF_NONE)  ; str = 'none'
  case (GF_STF_GAUSS) ; str = 'gaussian'
  case (GF_STF_HEAVI) ; str = 'heaviside'
  case default        ; str = 'unknown'
  end select

  end function gf_stf_kind_name

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_print_stf(stf,tax,iunit)

! reports a plan, in the 'key = value' layout of print_location

  implicit none

  type(t_gf_stf), intent(in) :: stf
  type(t_gf_taxis), intent(in) :: tax
  integer, intent(in) :: iunit

  write(iunit,'(a)')         'source time function conversion'
  write(iunit,'(a,a)')       '  kind                 = ',trim(gf_stf_kind_name(stf%kind_stf))
  write(iunit,'(a,es22.14)') '  hdur (source file)   = ',stf%hdur_src
  write(iunit,'(a,es22.14)') '  hdur_target          = ',stf%hdur_target
  write(iunit,'(a,es22.14)') '  hdur_db              = ',stf%hdur_db
  write(iunit,'(a,es22.14)') '  hdur_corr            = ',stf%hdur_corr
  write(iunit,'(a,es22.14)') '  truncation, x hdur   = ',stf%trunc
  write(iunit,'(a,i0)')      '  kernel half length   = ',stf%khalf
  write(iunit,'(a,l1)')      '  guard                = ',stf%guard
  write(iunit,'(a,a)')       '  note                 = ',trim(stf%note)
  write(iunit,'(a)')         ''
  write(iunit,'(a)')         'output time axis'
  write(iunit,'(a,i0)')      '  nt_db                = ',tax%nt_db
  write(iunit,'(a,i0)')      '  npad                 = ',tax%npad
  write(iunit,'(a,i0)')      '  nt                   = ',tax%nt
  write(iunit,'(a,i0)')      '  subsample_step       = ',tax%subsample_step
  write(iunit,'(a,es22.14)') '  dt                   = ',tax%dt
  write(iunit,'(a,es22.14)') '  dt_sub               = ',tax%dt_sub
  write(iunit,'(a,es22.14)') '  t0_db                = ',tax%t0_db
  write(iunit,'(a,es22.14)') '  t0 requested         = ',tax%t0_req
  write(iunit,'(a,es22.14)') '  t0                   = ',tax%t0
  write(iunit,'(a,es22.14)') '  t(1)                 = ',tax%t_first
  write(iunit,'(a,i0)')      '  trailing samples using zero-extended data = ',stf%khalf

  end subroutine gf_print_stf

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_cumtrapz(f,n,dt,g)

! cumulative trapezoidal integration on a uniform grid
!
!   g(1) = 0
!   g(i) = g(i-1) + (dt/2)(f(i-1) + f(i))
!
! Stage 4's integrator, no longer on the production path (see the module
! header for why); it serves --dump and is the oracle in test_gf_stf. Left
! byte for byte as Stage 4 wrote it.
!
! `dt` is the spacing of the axis being integrated along, which for a
! database trace is dt_sub = dt*subsample_step and *not* the solver step --
! 0.4 s rather than 0.1 s in the shipped global example. Passing the solver
! step here is wrong by a factor of subsample_step and produces a trace of
! entirely plausible shape.
!
! g may alias f; the running value is held in a scalar so that an in-place
! call is well defined.

  implicit none

  integer, intent(in) :: n
  double precision, dimension(n), intent(in) :: f
  double precision, intent(in) :: dt
  double precision, dimension(n), intent(out) :: g

  ! local parameters
  integer :: i
  double precision :: acc,fprev,fcur

  if (n < 1) return

  acc = 0.d0
  g(1) = 0.d0
  if (n == 1) return

  fprev = f(1)
  do i = 2,n
    fcur = f(i)
    acc = acc + 0.5d0*dt*(fprev + fcur)
    g(i) = acc
    fprev = fcur
  enddo

  end subroutine gf_cumtrapz

  end module gf_stf
