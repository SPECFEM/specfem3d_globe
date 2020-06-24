!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
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


  subroutine get_attenuation_model_3D(iregion_code, &
                                      one_minus_sum_beta, &
                                      factor_common, &
                                      factor_scale, tau_s, vnspec)

  use constants_solver

  use shared_parameters, only: ATT_F_C_SOURCE

  use specfem_par, only: ATTENUATION_VAL,ADIOS_FOR_ARRAYS_SOLVER,LOCAL_PATH, &
    scale_t_inv

  implicit none

  integer,intent(in) :: iregion_code

  ! note: factor_common,one_minus_sum_beta and factor_scale are real custom.
  !       this is better, it works fine and these arrays are really huge
  !       in the crust_mantle region, thus let us not double their size
  integer,intent(in) :: vnspec
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec),intent(out) :: one_minus_sum_beta, factor_scale
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,N_SLS,vnspec),intent(out) :: factor_common

  double precision, dimension(N_SLS),intent(out) :: tau_s

  ! local parameters
  integer :: i,j,k,ispec,ier,i_sls
  double precision, dimension(N_SLS) :: tau_e, fc
  double precision :: omsb, Q_mu, sf
  double precision :: f_c_source  ! frequency
  double precision :: T_c_source  ! period
  character(len=MAX_STRING_LEN) :: prname

  double precision, parameter :: TOL = 1.d-9

  ! checks if attenuation is on and anything to do
  if (.not. ATTENUATION_VAL) return
  if (.not. I_should_read_the_database) return

  ! initializes
  tau_s(:) = 0.d0
  factor_common(:,:,:,:,:) = 0._CUSTOM_REAL
  factor_scale(:,:,:,:) = 0._CUSTOM_REAL
  f_c_source = 0.d0

  ! All of the following reads use the output parameters as their temporary arrays
  ! use the filename to determine the actual contents of the read
  if (ADIOS_FOR_ARRAYS_SOLVER) then
    ! ADIOS format
    call read_attenuation_adios(iregion_code,factor_common, factor_scale, tau_s, vnspec, f_c_source)
  else
    ! binary format
    ! opens corresponding databases file
    call create_name_database(prname,myrank,iregion_code,LOCAL_PATH)

    open(unit=IIN, file=prname(1:len_trim(prname))//'attenuation.bin', &
          status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file attenuation.bin')

    read(IIN) tau_s
    read(IIN) factor_common ! tau_e_store
    read(IIN) factor_scale  ! Qmu_store
    read(IIN) f_c_source    ! center frequency
    close(IIN)
  endif

  ! checks
  if (f_c_source <= 0.d0) call exit_MPI(myrank,'Error: invalid attenuation center frequency read in from mesh file')

  ! since we read in crust/mantle region first, we check if the center frequency is the same for both regions
  if (iregion_code == IREGION_INNER_CORE) then
    if (abs(f_c_source - ATT_F_C_SOURCE) > TOL) then
      print *,'Error: different center frequencies for crust/mantle and inner core regions:'
      print *,'  crust/mantle center frequency: ',ATT_F_C_SOURCE
      print *,'  inner core   center frequency: ',f_c_source
      print *,'Please check if the mesh files are correct.'
      call exit_MPI(myrank,'Error different center frequencies for crust/mantle and inner core regions')
    endif
  endif

  ! stores center frequency as shared parameter
  ATT_F_C_SOURCE = f_c_source

  ! debug
  !print *,'debug: rank ',myrank,' attenuation center frequency = ',ATT_F_C_SOURCE,'region',iregion_code

  ! This is really tau_e, not factor_common
  factor_common(:,:,:,:,:) = real( factor_common(:,:,:,:,:) * scale_t_inv ,kind=CUSTOM_REAL)
  tau_s(:)                 = tau_s(:) * scale_t_inv

  ! converts center frequency to center period
  ! note: the mesher stores the center frequency f_c_source.
  !       here we invert it to have the center period T_c_source
  if (f_c_source > 0.d0) then
    T_c_source = 1.0d0 / f_c_source
  else
    T_c_source = 0.d0
  endif
  ! non-dimensionalizes
  T_c_source = T_c_source * scale_t_inv

  ! loops over elements
  do ispec = 1, vnspec
    ! loops over GLL points
    do k = 1, ATT3_VAL
      do j = 1, ATT2_VAL
        do i = 1, ATT1_VAL
          ! gets relaxation times for each standard linear solid
          do i_sls = 1,N_SLS
            tau_e(i_sls) = factor_common(i,j,k,i_sls,ispec)
          enddo
          Q_mu = factor_scale(i,j,k,ispec)

          ! Determine the factor_common and one_minus_sum_beta from tau_s and tau_e
          call get_attenuation_property_values(tau_s, tau_e, fc, omsb)

          do i_sls = 1,N_SLS
            factor_common(i,j,k,i_sls,ispec) = real(fc(i_sls), kind=CUSTOM_REAL)
          enddo

          ! unrelaxed moduli: factor to scale from relaxed to unrelaxed moduli
          one_minus_sum_beta(i,j,k,ispec) = real(omsb, kind=CUSTOM_REAL)

          ! physical dispersion factor: scales moduli from reference frequency to simulation (source) center frequency
          ! Determine the "factor_scale" from tau_s, tau_e, central source frequency, and Q
          call get_attenuation_scale_factor(T_c_source, tau_e, tau_s, Q_mu, sf)

          factor_scale(i,j,k,ispec) = real(sf, kind=CUSTOM_REAL)
        enddo
      enddo
    enddo
  enddo

  end subroutine get_attenuation_model_3D

!
!-------------------------------------------------------------------------------------------------
!
  subroutine get_attenuation_property_values(tau_s, tau_e, factor_common, one_minus_sum_beta)

  use constants

  implicit none

  double precision, dimension(N_SLS),intent(in) :: tau_s, tau_e
  double precision, dimension(N_SLS),intent(out) :: factor_common
  double precision,intent(out) ::  one_minus_sum_beta

  ! local parameters
  double precision, dimension(N_SLS) :: tauinv,beta
  integer :: i

  tauinv(:) = 0.d0
  where(abs(tau_s(:)) > 0.d0) tauinv(:) = - 1.0d0 / tau_s(:)

  beta(:) = 1.0d0 + tau_e(:) * tauinv(:)     ! 1 - tau_e / tau_s

  ! factor to scale from relaxed to unrelaxed moduli: see e.g. Komatitsch & Tromp 1999, eq. (7)
  one_minus_sum_beta = 1.0d0
  do i = 1,N_SLS
     one_minus_sum_beta = one_minus_sum_beta - beta(i)
  enddo

!ZN beware, here the expression differs from the strain used in memory variable equation (6) in D. Komatitsch and J. Tromp 1999,
!ZN here Brian Savage uses the engineering strain which are epsilon = 1/2*(grad U + (grad U)^T),
!ZN where U is the displacement vector and grad the gradient operator, i.e. there is a 1/2 factor difference between the two.
!ZN Both expressions are fine, but we need to keep in mind that if we have put the 1/2 factor there we need to remove it
!ZN from the expression in which we use the strain here in the code.
!ZN This is why here Brian Savage multiplies beta(:) * tauinv(:) by 2.0 to compensate for the 1/2 factor used before
  factor_common(:) = 2.0d0 * beta(:) * tauinv(:)

  end subroutine get_attenuation_property_values

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_attenuation_scale_factor(T_c_source, tau_mu, tau_sigma, Q_mu, scale_factor)

  use constants, only: ZERO,ONE,TWO,PI,TWO_PI,N_SLS,ATTENUATION_f0_REFERENCE,myrank
  use specfem_par, only: scale_t

  implicit none

  double precision,intent(in) :: T_c_source  ! center period
  double precision,intent(in) :: Q_mu
  double precision, dimension(N_SLS),intent(in) :: tau_mu, tau_sigma

  double precision,intent(out) :: scale_factor

  ! local parameters
  double precision :: f_c_source, w_c_source,f_0_model
  double precision :: factor_scale_mu0, factor_scale_mu
  double precision :: a_val, b_val, denom
  double precision :: big_omega
  integer :: i

  !--- compute central angular frequency of source (non dimensionalized)
  if (T_c_source > 0.d0) then
    f_c_source = ONE / T_c_source
  else
    f_c_source = 0.d0
  endif
  w_c_source = TWO_PI * f_c_source

  !--- non dimensionalize (e.g., PREM reference of 1 second)
  f_0_model = ATTENUATION_f0_REFERENCE * scale_t  ! original f_0_prem = ONE / ( ONE / scale_t)

!--- quantity by which to scale mu_0 to get mu
! this formula can be found for instance in
! Liu, H. P., Anderson, D. L. and Kanamori, H., Velocity dispersion due to
! anelasticity: implications for seismology and mantle composition,
! Geophys. J. R. Astron. Soc., vol. 47, pp. 41-58 (1976)
! and in Aki, K. and Richards, P. G., Quantitative seismology, theory and methods,
! W. H. Freeman, (1980), second edition, sections 5.5 and 5.5.2, eq. (5.81) p. 170
  if (f_0_model > 0.d0) then
    factor_scale_mu0 = ONE + TWO * log(f_c_source / f_0_model) / (PI * Q_mu)
  else
    factor_scale_mu0 = ONE
  endif

  !--- compute a, b and Omega parameters, also compute one minus sum of betas
  a_val = ONE
  b_val = ZERO

  do i = 1,N_SLS
    denom = 1.d0 + w_c_source * w_c_source * tau_mu(i) * tau_mu(i)
    if (denom /= 0.d0) then
      a_val = a_val - w_c_source * w_c_source * tau_mu(i) * (tau_mu(i) - tau_sigma(i)) / denom
      b_val = b_val + w_c_source * (tau_mu(i) - tau_sigma(i)) / denom
    endif
  enddo

  if (a_val /= 0.d0) then
    big_omega = a_val*(sqrt(1.d0 + b_val*b_val/(a_val*a_val))-1.d0)
  else
    big_omega = 0.d0
  endif

  !--- quantity by which to scale mu to get mu_relaxed
  if (big_omega /= 0.d0) then
    factor_scale_mu = b_val * b_val / (TWO * big_omega)
  else
    factor_scale_mu = ONE
  endif

  !--- total factor by which to scale mu0
  scale_factor = factor_scale_mu * factor_scale_mu0

  !--- check that the correction factor is close to one
  if (scale_factor < 0.8d0 .or. scale_factor > 1.2d0) then
    print *,'Error: incorrect scale factor: ', scale_factor
    call exit_MPI(myrank,'incorrect correction factor in attenuation model')
  endif

  end subroutine get_attenuation_scale_factor

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_attenuation_memory_values(tau_s,deltat, alphaval,betaval,gammaval)

  use constants

  implicit none

  double precision, dimension(N_SLS), intent(in) :: tau_s
  double precision, dimension(N_SLS), intent(out) :: alphaval, betaval,gammaval
  real(kind=CUSTOM_REAL), intent(in) :: deltat

  ! local parameters
  double precision, dimension(N_SLS) :: tauinv

  ! inverse of tau_s
  tauinv(:) = - 1.d0 / tau_s(:)

  ! runge-kutta coefficients
  ! see e.g.: Savage et al. (BSSA, 2010): eq. (11)
  alphaval(:) = 1.d0 + deltat * tauinv(:) &
                  + deltat**2 * tauinv(:)**2 / 2.d0 &
                  + deltat**3 * tauinv(:)**3 / 6.d0 &
                  + deltat**4 * tauinv(:)**4 / 24.d0
  betaval(:)  = deltat / 2.d0 &
                  + deltat**2 * tauinv(:) / 3.d0 &
                  + deltat**3 * tauinv(:)**2 / 8.d0 &
                  + deltat**4 * tauinv(:)**3 / 24.d0
  gammaval(:) = deltat / 2.d0 &
                  + deltat**2 * tauinv(:) / 6.d0 &
                  + deltat**3 * tauinv(:)**2 / 24.d0

  end subroutine get_attenuation_memory_values

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_attenuation_model_3D(one_minus_sum_beta,factor_common,factor_scale,tau_s,vnspec)

  use constants_solver
  use shared_parameters, only: ATT_F_C_SOURCE

  implicit none

  integer :: vnspec

  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec) :: one_minus_sum_beta, factor_scale
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,N_SLS,vnspec) :: factor_common
  double precision, dimension(N_SLS) :: tau_s

  call bcast_all_cr_for_database(one_minus_sum_beta(1,1,1,1), size(one_minus_sum_beta))
  call bcast_all_cr_for_database(factor_scale(1,1,1,1), size(factor_scale))
  call bcast_all_cr_for_database(factor_common(1,1,1,1,1), size(factor_common))
  call bcast_all_dp_for_database(tau_s(1), size(tau_s))
  call bcast_all_dp_for_database(ATT_F_C_SOURCE, 1)

  end subroutine bcast_attenuation_model_3D
