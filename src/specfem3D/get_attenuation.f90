!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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


  subroutine get_attenuation_model_3D(myrank, prname, one_minus_sum_beta, &
                                factor_common, scale_factor, tau_s, vnspec)

  implicit none

  include 'constants.h'

  integer myrank, vnspec
  character(len=150) prname
  double precision, dimension(NGLLX,NGLLY,NGLLZ,vnspec)       :: one_minus_sum_beta, scale_factor
  double precision, dimension(N_SLS,NGLLX,NGLLY,NGLLZ,vnspec) :: factor_common
  double precision, dimension(N_SLS)                          :: tau_s

  integer i,j,k,ispec

  double precision, dimension(N_SLS) :: tau_e, fc
  double precision  omsb, Q_mu, sf, T_c_source, scale_t

  ! All of the following reads use the output parameters as their temporary arrays
  ! use the filename to determine the actual contents of the read
  open(unit=27, file=prname(1:len_trim(prname))//'attenuation.bin', &
        status='old',action='read',form='unformatted')
  read(27) tau_s
  read(27) factor_common
  read(27) scale_factor
  read(27) T_c_source
  close(27)

  scale_t = ONE/dsqrt(PI*GRAV*RHOAV)

  factor_common(:,:,:,:,:) = factor_common(:,:,:,:,:) / scale_t ! This is really tau_e, not factor_common
  tau_s(:)                 = tau_s(:) / scale_t
  T_c_source               = 1000.0d0 / T_c_source
  T_c_source               = T_c_source / scale_t

  do ispec = 1, vnspec
     do k = 1, NGLLZ
        do j = 1, NGLLY
           do i = 1, NGLLX
              tau_e(:) = factor_common(:,i,j,k,ispec)
              Q_mu     = scale_factor(i,j,k,ispec)

              ! Determine the factor_common and one_minus_sum_beta from tau_s and tau_e
              call get_attenuation_property_values(tau_s, tau_e, fc, omsb)

              factor_common(:,i,j,k,ispec)    = fc(:)
              one_minus_sum_beta(i,j,k,ispec) = omsb

              ! Determine the "scale_factor" from tau_s, tau_e, central source frequency, and Q
              call get_attenuation_scale_factor(myrank, T_c_source, tau_e, tau_s, Q_mu, sf)
              scale_factor(i,j,k,ispec) = sf
           enddo
        enddo
     enddo
  enddo

  end subroutine get_attenuation_model_3D

!
!-------------------------------------------------------------------------------------------------
!
  subroutine get_attenuation_property_values(tau_s, tau_e, factor_common, one_minus_sum_beta)

  implicit none

  include 'constants.h'

  double precision, dimension(N_SLS) :: tau_s, tau_e, beta, factor_common
  double precision  one_minus_sum_beta

  double precision, dimension(N_SLS) :: tauinv
  integer i

  tauinv(:) = -1.0d0 / tau_s(:)

  beta(:) = 1.0d0 - tau_e(:) / tau_s(:)
  one_minus_sum_beta = 1.0d0

  do i = 1,N_SLS
     one_minus_sum_beta = one_minus_sum_beta - beta(i)
  enddo

  factor_common(:) = 2.0d0 * beta(:) * tauinv(:)

  end subroutine get_attenuation_property_values

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_attenuation_scale_factor(myrank, T_c_source, tau_mu, tau_sigma, Q_mu, scale_factor)

  implicit none

  include 'constants.h'

  integer myrank
  double precision scale_factor, Q_mu, T_c_source
  double precision, dimension(N_SLS) :: tau_mu, tau_sigma

  double precision scale_t
  double precision f_c_source, w_c_source, f_0_prem
  double precision factor_scale_mu0, factor_scale_mu
  double precision a_val, b_val
  double precision big_omega
  integer i

  scale_t = ONE/dsqrt(PI*GRAV*RHOAV)

  !--- compute central angular frequency of source (non dimensionalized)
  f_c_source = ONE / T_c_source
  w_c_source = TWO_PI * f_c_source

  !--- non dimensionalize PREM reference of 1 second
  f_0_prem = ONE / ( ONE / scale_t)

!--- quantity by which to scale mu_0 to get mu
! this formula can be found for instance in
! Liu, H. P., Anderson, D. L. and Kanamori, H., Velocity dispersion due to
! anelasticity: implications for seismology and mantle composition,
! Geophys. J. R. Astron. Soc., vol. 47, pp. 41-58 (1976)
! and in Aki, K. and Richards, P. G., Quantitative seismology, theory and methods,
! W. H. Freeman, (1980), second edition, sections 5.5 and 5.5.2, eq. (5.81) p. 170
  factor_scale_mu0 = ONE + TWO * log(f_c_source / f_0_prem) / (PI * Q_mu)

  !--- compute a, b and Omega parameters, also compute one minus sum of betas
  a_val = ONE
  b_val = ZERO

  do i = 1,N_SLS
    a_val = a_val - w_c_source * w_c_source * tau_mu(i) * &
      (tau_mu(i) - tau_sigma(i)) / (1.d0 + w_c_source * w_c_source * tau_mu(i) * tau_mu(i))
    b_val = b_val + w_c_source * (tau_mu(i) - tau_sigma(i)) / &
      (1.d0 + w_c_source * w_c_source * tau_mu(i) * tau_mu(i))
  enddo

  big_omega = a_val*(sqrt(1.d0 + b_val*b_val/(a_val*a_val))-1.d0)

  !--- quantity by which to scale mu to get mu_relaxed
  factor_scale_mu = b_val * b_val / (TWO * big_omega)

  !--- total factor by which to scale mu0
  scale_factor = factor_scale_mu * factor_scale_mu0

  !--- check that the correction factor is close to one
  if(scale_factor < 0.8 .or. scale_factor > 1.2) then
     write(*,*)'scale factor: ', scale_factor
     call exit_MPI(myrank,'incorrect correction factor in attenuation model')
  endif

  end subroutine get_attenuation_scale_factor


!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_attenuation_memory_values(tau_s, deltat, alphaval,betaval,gammaval)

  implicit none

  include 'constants.h'

  double precision, dimension(N_SLS) :: tau_s, alphaval, betaval,gammaval
  real(kind=CUSTOM_REAL) deltat

  double precision, dimension(N_SLS) :: tauinv

  tauinv(:) = - 1.0 / tau_s(:)

  alphaval(:)  = 1 + deltat*tauinv(:) + deltat**2*tauinv(:)**2 / 2. + &
                    deltat**3*tauinv(:)**3 / 6. + deltat**4*tauinv(:)**4 / 24.
  betaval(:)   = deltat / 2. + deltat**2*tauinv(:) / 3. &
                + deltat**3*tauinv(:)**2 / 8. + deltat**4*tauinv(:)**3 / 24.
  gammaval(:)  = deltat / 2. + deltat**2*tauinv(:) / 6. &
                + deltat**3*tauinv(:)**2 / 24.0

  end subroutine get_attenuation_memory_values


!
!-------------------------------------------------------------------------------------------------
!
! not used anymore...
!
!  subroutine get_attenuation_model_1D(myrank, prname, iregion_code, tau_s, one_minus_sum_beta, &
!                                    factor_common, scale_factor, vn,vx,vy,vz, AM_V)
!
!  implicit none
!
!  include 'mpif.h'
!  include 'constants.h'
!
!! model_attenuation_variables
!  type model_attenuation_variables
!    sequence
!    double precision min_period, max_period
!    double precision                          :: QT_c_source        ! Source Frequency
!    double precision, dimension(:), pointer   :: Qtau_s             ! tau_sigma
!    double precision, dimension(:), pointer   :: QrDisc             ! Discontinutitues Defined
!    double precision, dimension(:), pointer   :: Qr                 ! Radius
!    integer, dimension(:), pointer            :: interval_Q                 ! Steps
!    double precision, dimension(:), pointer   :: Qmu                ! Shear Attenuation
!    double precision, dimension(:,:), pointer :: Qtau_e             ! tau_epsilon
!    double precision, dimension(:), pointer   :: Qomsb, Qomsb2      ! one_minus_sum_beta
!    double precision, dimension(:,:), pointer :: Qfc, Qfc2          ! factor_common
!    double precision, dimension(:), pointer   :: Qsf, Qsf2          ! scale_factor
!    integer, dimension(:), pointer            :: Qrmin              ! Max and Mins of idoubling
!    integer, dimension(:), pointer            :: Qrmax              ! Max and Mins of idoubling
!    integer                                   :: Qn                 ! Number of points
!    integer dummy_pad ! padding 4 bytes to align the structure
!  end type model_attenuation_variables
!
!  type (model_attenuation_variables) AM_V
!! model_attenuation_variables
!
!  integer myrank, iregion_code
!  character(len=150) prname
!  integer vn, vx,vy,vz
!  double precision, dimension(N_SLS)              :: tau_s
!  double precision, dimension(vx,vy,vz,vn)        :: scale_factor, one_minus_sum_beta
!  double precision, dimension(N_SLS, vx,vy,vz,vn) :: factor_common
!
!  integer i,j,ier,rmax
!  double precision scale_t
!  double precision Qp1, Qpn, radius, fctmp
!  double precision, dimension(:), allocatable :: Qfctmp, Qfc2tmp
!
!  integer, save :: first_time_called = 1
!
!  if(myrank == 0 .AND. iregion_code == IREGION_CRUST_MANTLE .AND. first_time_called == 1) then
!     first_time_called = 0
!     open(unit=27, file=prname(1:len_trim(prname))//'1D_Q.bin', status='unknown', form='unformatted')
!     read(27) AM_V%QT_c_source
!     read(27) tau_s
!     read(27) AM_V%Qn
!
!     allocate(AM_V%Qr(AM_V%Qn))
!     allocate(AM_V%Qmu(AM_V%Qn))
!     allocate(AM_V%Qtau_e(N_SLS,AM_V%Qn))
!
!     read(27) AM_V%Qr
!     read(27) AM_V%Qmu
!     read(27) AM_V%Qtau_e
!     close(27)
!  endif
!
!  ! Synch up after the Read
!  call MPI_BCAST(AM_V%QT_c_source,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
!  call MPI_BCAST(tau_s,N_SLS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
!  call MPI_BCAST(AM_V%Qn,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
!
!  if(myrank /= 0) then
!     allocate(AM_V%Qr(AM_V%Qn))
!     allocate(AM_V%Qmu(AM_V%Qn))
!     allocate(AM_V%Qtau_e(N_SLS,AM_V%Qn))
!  endif
!
!  call MPI_BCAST(AM_V%Qr,AM_V%Qn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
!  call MPI_BCAST(AM_V%Qmu,AM_V%Qn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
!  call MPI_BCAST(AM_V%Qtau_e,AM_V%Qn*N_SLS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
!
!  scale_t = ONE/dsqrt(PI*GRAV*RHOAV)
!
!  ! Scale the Attenuation Values
!  tau_s(:) = tau_s(:) / scale_t
!  AM_V%Qtau_e(:,:) = AM_V%Qtau_e(:,:) / scale_t
!  AM_V%QT_c_source = 1000.0d0 / AM_V%QT_c_source / scale_t
!  AM_V%Qr(:) = AM_V%Qr(:) / R_EARTH
!
!  allocate(AM_V%Qsf(AM_V%Qn))
!  allocate(AM_V%Qomsb(AM_V%Qn))
!  allocate(AM_V%Qfc(N_SLS,AM_V%Qn))
!
!  allocate(AM_V%Qsf2(AM_V%Qn))
!  allocate(AM_V%Qomsb2(AM_V%Qn))
!  allocate(AM_V%Qfc2(N_SLS,AM_V%Qn))
!
!  allocate(AM_V%interval_Q(AM_V%Qn))
!
!  allocate(Qfctmp(AM_V%Qn))
!  allocate(Qfc2tmp(AM_V%Qn))
!
!  do i = 1,AM_V%Qn
!     if(AM_V%Qmu(i) == 0.0d0) then
!        AM_V%Qomsb(i) = 0.0d0
!        AM_V%Qfc(:,i) = 0.0d0
!        AM_V%Qsf(i)   = 0.0d0
!     else
!        call attenuation_property_values(tau_s, AM_V%Qtau_e(:,i), AM_V%Qfc(:,i), AM_V%Qomsb(i))
!        call attenuation_scale_factor(myrank, AM_V%QT_c_source, AM_V%Qtau_e(:,i), tau_s, AM_V%Qmu(i), AM_V%Qsf(i))
!     endif
!  enddo
!
!  ! Determine the Spline Coefficients or Second Derivatives
!  call pspline_construction(AM_V%Qr, AM_V%Qsf,   AM_V%Qn, Qp1, Qpn, AM_V%Qsf2,   AM_V%interval_Q)
!  call pspline_construction(AM_V%Qr, AM_V%Qomsb, AM_V%Qn, Qp1, Qpn, AM_V%Qomsb2, AM_V%interval_Q)
!  do i = 1,N_SLS
!! copy the sub-arrays to temporary arrays to avoid a warning by some compilers
!! about temporary arrays being created automatically when using this expression
!! directly in the call to the subroutine
!     Qfctmp(:) = AM_V%Qfc(i,:)
!     Qfc2tmp(:) = AM_V%Qfc2(i,:)
!     call pspline_construction(AM_V%Qr, Qfctmp, AM_V%Qn, Qp1, Qpn, Qfc2tmp, AM_V%interval_Q)
!! copy the arrays back to the sub-arrays, since these sub-arrays are used
!! as input and output
!     AM_V%Qfc(i,:) = Qfctmp(:)
!     AM_V%Qfc2(i,:) = Qfc2tmp(:)
!  enddo
!
!  radius = 0.0d0
!  rmax = nint(TABLE_ATTENUATION)
!  do i = 1,rmax
!     call attenuation_lookup_value(i, radius)
!     call pspline_evaluation(AM_V%Qr, AM_V%Qsf,   AM_V%Qsf2,   AM_V%Qn, radius, scale_factor(1,1,1,i),       AM_V%interval_Q)
!     call pspline_evaluation(AM_V%Qr, AM_V%Qomsb, AM_V%Qomsb2, AM_V%Qn, radius, one_minus_sum_beta(1,1,1,i), AM_V%interval_Q)
!     do j = 1,N_SLS
!        Qfctmp  = AM_V%Qfc(j,:)
!        Qfc2tmp = AM_V%Qfc2(j,:)
!        call pspline_evaluation(AM_V%Qr, Qfctmp, Qfc2tmp, AM_V%Qn, radius, fctmp, AM_V%interval_Q)
!        factor_common(j,1,1,1,i) = fctmp
!     enddo
!  enddo
!  do i = rmax+1,NRAD_ATTENUATION
!     scale_factor(1,1,1,i)       = scale_factor(1,1,1,rmax)
!     one_minus_sum_beta(1,1,1,i) = one_minus_sum_beta(1,1,1,rmax)
!     factor_common(1,1,1,1,i)    = factor_common(1,1,1,1,rmax)
!     factor_common(2,1,1,1,i)    = factor_common(2,1,1,1,rmax)
!     factor_common(3,1,1,1,i)    = factor_common(3,1,1,1,rmax)
!  enddo
!
!  deallocate(AM_V%Qfc2)
!  deallocate(AM_V%Qsf2)
!  deallocate(AM_V%Qomsb2)
!  deallocate(AM_V%Qfc)
!  deallocate(AM_V%Qsf)
!  deallocate(AM_V%Qomsb)
!  deallocate(AM_V%Qtau_e)
!  deallocate(Qfctmp)
!  deallocate(Qfc2tmp)
!
!  call MPI_BARRIER(MPI_COMM_WORLD, ier)
!
!  end subroutine get_attenuation_model_1D
!
!
!-------------------------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------------------------
!
! not used anymore...
!
! Piecewise Continuous Splines
!   - Added Steps which describes the discontinuities
!   - Steps must be repeats in the dependent variable, X
!   - Derivates at the steps are computed using the point
!     at the derivate and the closest point within that piece
!   - A point lying directly on the discontinuity will recieve the
!     value of the first or smallest piece in terms of X
!   - Beginning and Ending points of the Function become beginning
!     and ending points of the first and last splines
!   - A Step with a value of zero is undefined
!   - Works with functions with steps or no steps
! See the comment below about the ScS bug
!  subroutine pspline_evaluation(xa, ya, y2a, n, x, y, steps)
!
!  implicit none
!
!  integer n
!  double precision xa(n),ya(n),y2a(n)
!  integer steps(n)
!  double precision x, y
!
!  integer i, l, n1, n2
!
!  do i = 1,n-1,1
!     if(steps(i+1) == 0) return
!     if(x >= xa(steps(i)) .and. x <= xa(steps(i+1))) then
!        call pspline_piece(i,n1,n2,l,n,steps)
!        call spline_evaluation(xa(n1), ya(n1), y2a(n1), l, x, y)
!!        return <-- Commented out to fix ScS bug
!     endif
!  enddo
!
!  end subroutine pspline_evaluation
!
!
!-------------------------------------------------------------------------------------------------
!
! not used anymore...
!
!  subroutine pspline_piece(i,n1,n2,l,n,s)
!
!  implicit none
!
!  integer i, n1, n2, l, n, s(n)
!  n1 = s(i)+1
!  if(i == 1) n1 = s(i)
!  n2 = s(i+1)
!  l = n2 - n1 + 1
!
!  end subroutine pspline_piece
!
!
!-------------------------------------------------------------------------------------------------
!
! not used anymore...
!
!  subroutine pspline_construction(x, y, n, yp1, ypn, y2, steps)
!
!  implicit none
!
!  integer n
!  double precision x(n),y(n),y2(n)
!  double precision yp1, ypn
!  integer steps(n)
!
!  integer i,r, l, n1,n2
!
!  steps(:) = 0
!
!  ! Find steps in x, defining pieces
!  steps(1) = 1
!  r = 2
!  do i = 2,n
!     if(x(i) == x(i-1)) then
!        steps(r) = i-1
!        r = r + 1
!     endif
!  end do
!  steps(r) = n
!
!  ! Run spline for each piece
!  do i = 1,r-1
!     call pspline_piece(i,n1,n2,l,n,steps)
!     ! Determine the First Derivates at Begin/End Points
!     yp1 = ( y(n1+1) - y(n1) ) / ( x(n1+1) - x(n1))
!     ypn = ( y(n2) - y(n2-1) ) / ( x(n2) - x(n2-1))
!     call spline_construction(x(n1),y(n1),l,yp1,ypn,y2(n1))
!  enddo
!
!  end subroutine pspline_construction
!
!
!-------------------------------------------------------------------------------------------------
!
!
! not used anymore...
!
!  subroutine attenuation_lookup_value(i, r)
!
!  implicit none
!
!  include 'constants.h'
!
!  integer i
!  double precision r
!
!  r = dble(i) / TABLE_ATTENUATION
!
!  end subroutine attenuation_lookup_value
!
!
!-------------------------------------------------------------------------------------------------
!
! not used anymore...
!
!  subroutine attenuation_save_arrays(prname, iregion_code, AM_V)
!
!  implicit none
!
!  include 'mpif.h'
!  include 'constants.h'
!
!! model_attenuation_variables
!  type model_attenuation_variables
!    sequence
!    double precision min_period, max_period
!    double precision                          :: QT_c_source        ! Source Frequency
!    double precision, dimension(:), pointer   :: Qtau_s             ! tau_sigma
!    double precision, dimension(:), pointer   :: QrDisc             ! Discontinutitues Defined
!    double precision, dimension(:), pointer   :: Qr                 ! Radius
!    integer, dimension(:), pointer            :: interval_Q                 ! Steps
!    double precision, dimension(:), pointer   :: Qmu                ! Shear Attenuation
!    double precision, dimension(:,:), pointer :: Qtau_e             ! tau_epsilon
!    double precision, dimension(:), pointer   :: Qomsb, Qomsb2      ! one_minus_sum_beta
!    double precision, dimension(:,:), pointer :: Qfc, Qfc2          ! factor_common
!    double precision, dimension(:), pointer   :: Qsf, Qsf2          ! scale_factor
!    integer, dimension(:), pointer            :: Qrmin              ! Max and Mins of idoubling
!    integer, dimension(:), pointer            :: Qrmax              ! Max and Mins of idoubling
!    integer                                   :: Qn                 ! Number of points
!  end type model_attenuation_variables
!
!  type (model_attenuation_variables) AM_V
!! model_attenuation_variables
!
!  integer iregion_code
!  character(len=150) prname
!  integer ier
!  integer myrank
!  integer, save :: first_time_called = 1
!
!  call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ier)
!  if(myrank == 0 .AND. iregion_code == IREGION_CRUST_MANTLE .AND. first_time_called == 1) then
!    first_time_called = 0
!    open(unit=27,file=prname(1:len_trim(prname))//'1D_Q.bin',status='unknown',form='unformatted')
!    write(27) AM_V%QT_c_source
!    write(27) AM_V%Qtau_s
!    write(27) AM_V%Qn
!    write(27) AM_V%Qr
!    write(27) AM_V%Qmu
!    write(27) AM_V%Qtau_e
!    close(27)
!  endif
!
!  end subroutine attenuation_save_arrays
!
!
!-------------------------------------------------------------------------------------------------
!
! not used anymore...
!
!  subroutine get_attenuation_index(iflag, radius, index, inner_core, AM_V)
!
!  implicit none
!
!  include 'constants.h'
!
!! model_attenuation_variables
!  type model_attenuation_variables
!    sequence
!    double precision min_period, max_period
!    double precision                          :: QT_c_source        ! Source Frequency
!    double precision, dimension(:), pointer   :: Qtau_s             ! tau_sigma
!    double precision, dimension(:), pointer   :: QrDisc             ! Discontinutitues Defined
!    double precision, dimension(:), pointer   :: Qr                 ! Radius
!    integer, dimension(:), pointer            :: interval_Q                 ! Steps
!    double precision, dimension(:), pointer   :: Qmu                ! Shear Attenuation
!    double precision, dimension(:,:), pointer :: Qtau_e             ! tau_epsilon
!    double precision, dimension(:), pointer   :: Qomsb, Qomsb2      ! one_minus_sum_beta
!    double precision, dimension(:,:), pointer :: Qfc, Qfc2          ! factor_common
!    double precision, dimension(:), pointer   :: Qsf, Qsf2          ! scale_factor
!    integer, dimension(:), pointer            :: Qrmin              ! Max and Mins of idoubling
!    integer, dimension(:), pointer            :: Qrmax              ! Max and Mins of idoubling
!    integer                                   :: Qn                 ! Number of points
!  end type model_attenuation_variables
!
!  type (model_attenuation_variables) AM_V
!! model_attenuation_variables
!
!  integer iflag, iregion, index
!  double precision radius
!
!  ! Inner Core or not
!  logical inner_core
!
!  index = nint(radius * TABLE_ATTENUATION)
!
!!! DK DK this seems incorrect and is difficult to read anyway
!!! DK DK therefore let me rewrite it better
!! if(inner_core) then
!!   if(iflag >= IFLAG_INNER_CORE_NORMAL) then
!!     iregion = IREGION_ATTENUATION_INNER_CORE
!!   else if(iflag >= IFLAG_OUTER_CORE_NORMAL) then
!!     iregion = 6
!!   endif
!! else
!!   if(iflag >= IFLAG_MANTLE_NORMAL) then
!!     iregion = IREGION_ATTENUATION_CMB_670
!!   else if(iflag == IFLAG_670_220) then
!!     iregion = IREGION_ATTENUATION_670_220
!!   else if(iflag <= IFLAG_220_80) then
!!     iregion = IREGION_ATTENUATION_220_80
!!   else
!!     iregion = IREGION_ATTENUATION_80_SURFACE
!!   endif
!! endif
!  if(inner_core) then
!
!    if(iflag == IFLAG_INNER_CORE_NORMAL .or. iflag == IFLAG_MIDDLE_CENTRAL_CUBE .or. &
!       iflag == IFLAG_BOTTOM_CENTRAL_CUBE .or. iflag == IFLAG_TOP_CENTRAL_CUBE .or. &
!       iflag == IFLAG_IN_FICTITIOUS_CUBE) then
!      iregion = IREGION_ATTENUATION_INNER_CORE
!    else
!! this is fictitious for the outer core, which has no Qmu attenuation since it is fluid
!!      iregion = IREGION_ATTENUATION_80_SURFACE + 1
!       iregion = IREGION_ATTENUATION_UNDEFINED
!    endif
!
!  else
!
!    if(iflag == IFLAG_MANTLE_NORMAL) then
!      iregion = IREGION_ATTENUATION_CMB_670
!    else if(iflag == IFLAG_670_220) then
!      iregion = IREGION_ATTENUATION_670_220
!    else if(iflag == IFLAG_220_80) then
!      iregion = IREGION_ATTENUATION_220_80
!    else if(iflag == IFLAG_CRUST .or. iflag == IFLAG_80_MOHO) then
!      iregion = IREGION_ATTENUATION_80_SURFACE
!    else
!! this is fictitious for the outer core, which has no Qmu attenuation since it is fluid
!!      iregion = IREGION_ATTENUATION_80_SURFACE + 1
!       iregion = IREGION_ATTENUATION_UNDEFINED
!    endif
!
!  endif
!
!! Clamp regions
!  if(index < AM_V%Qrmin(iregion)) index = AM_V%Qrmin(iregion)
!  if(index > AM_V%Qrmax(iregion)) index = AM_V%Qrmax(iregion)
!
!  end subroutine get_attenuation_index
!
!
!-------------------------------------------------------------------------------------------------
!
! not used anymore...
!
!  subroutine set_attenuation_regions_1D(RICB, RCMB, R670, R220, R80, AM_V)
!
!  implicit none
!
!  include 'constants.h'
!
!! model_attenuation_variables
!  type model_attenuation_variables
!    sequence
!    double precision min_period, max_period
!    double precision                          :: QT_c_source        ! Source Frequency
!    double precision, dimension(:), pointer   :: Qtau_s             ! tau_sigma
!    double precision, dimension(:), pointer   :: QrDisc             ! Discontinutitues Defined
!    double precision, dimension(:), pointer   :: Qr                 ! Radius
!    integer, dimension(:), pointer            :: interval_Q                 ! Steps
!    double precision, dimension(:), pointer   :: Qmu                ! Shear Attenuation
!    double precision, dimension(:,:), pointer :: Qtau_e             ! tau_epsilon
!    double precision, dimension(:), pointer   :: Qomsb, Qomsb2      ! one_minus_sum_beta
!    double precision, dimension(:,:), pointer :: Qfc, Qfc2          ! factor_common
!    double precision, dimension(:), pointer   :: Qsf, Qsf2          ! scale_factor
!    integer, dimension(:), pointer            :: Qrmin              ! Max and Mins of idoubling
!    integer, dimension(:), pointer            :: Qrmax              ! Max and Mins of idoubling
!    integer                                   :: Qn                 ! Number of points
!  end type model_attenuation_variables
!
!  type (model_attenuation_variables) AM_V
!! model_attenuation_variables
!
!  double precision RICB, RCMB, R670, R220, R80
!  integer i
!
!  allocate(AM_V%Qrmin(6))
!  allocate(AM_V%Qrmax(6))
!  allocate(AM_V%QrDisc(5))
!
!  AM_V%QrDisc(1) = RICB
!  AM_V%QrDisc(2) = RCMB
!  AM_V%QrDisc(3) = R670
!  AM_V%QrDisc(4) = R220
!  AM_V%QrDisc(5) = R80
!
!   ! INNER CORE
!  AM_V%Qrmin(IREGION_ATTENUATION_INNER_CORE) = 1      ! Center of the Earth
!     i = nint(RICB / 100.d0)   ! === BOUNDARY === INNER CORE / OUTER CORE
!  AM_V%Qrmax(IREGION_ATTENUATION_INNER_CORE) = i - 1  ! Inner Core Boundary (Inner)
!
!  ! OUTER_CORE
!  AM_V%Qrmin(6) = i ! Inner Core Boundary (Outer)
!      i = nint(RCMB / 100.d0)  ! === BOUNDARY === INNER CORE / OUTER CORE
!  AM_V%Qrmax(6) = i - 1
!
!  ! LOWER MANTLE
!  AM_V%Qrmin(IREGION_ATTENUATION_CMB_670) = i
!       i = nint(R670 / 100.d0) ! === BOUNDARY === 670 km
!  AM_V%Qrmax(IREGION_ATTENUATION_CMB_670) = i - 1
!
!  ! UPPER MANTLE
!  AM_V%Qrmin(IREGION_ATTENUATION_670_220) = i
!       i = nint(R220 / 100.d0) ! === BOUNDARY === 220 km
!  AM_V%Qrmax(IREGION_ATTENUATION_670_220) = i - 1
!
!  ! MANTLE ISH LITHOSPHERE
!  AM_V%Qrmin(IREGION_ATTENUATION_220_80) = i
!       i = nint(R80 / 100.d0) ! === BOUNDARY === 80 km
!  AM_V%Qrmax(IREGION_ATTENUATION_220_80) = i - 1
!
!  ! CRUST ISH LITHOSPHERE
!  AM_V%Qrmin(IREGION_ATTENUATION_80_SURFACE) = i
!  AM_V%Qrmax(IREGION_ATTENUATION_80_SURFACE) = NRAD_ATTENUATION
!
!  end subroutine set_attenuation_regions_1D
!

