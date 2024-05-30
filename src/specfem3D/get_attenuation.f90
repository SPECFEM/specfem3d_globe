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


  subroutine get_attenuation_model_3D(iregion_code, &
                                      one_minus_sum_beta, &
                                      factor_common, &
                                      factor_scale, &
                                      factor_scale_relaxed, &
                                      tau_s, &
                                      vnspec, f0_reference)

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
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec),intent(out) :: factor_scale_relaxed
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,N_SLS,vnspec),intent(out) :: factor_common

  double precision, dimension(N_SLS),intent(out) :: tau_s

  ! model reference frequency
  double precision,intent(in) :: f0_reference

  ! local parameters
  integer :: i,j,k,ispec,ier,i_sls,nspec_region
  double precision, dimension(N_SLS) :: tau_e, fc
  double precision :: omsb, Q_mu, sf, sf_relaxed
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
  factor_scale_relaxed(:,:,:,:) = 0._CUSTOM_REAL
  f_c_source = 0.d0

  ! gets nspec for this region
  nspec_region = 0
  select case(iregion_code)
  case (IREGION_CRUST_MANTLE)
    nspec_region = NSPEC_CRUST_MANTLE
  case (IREGION_OUTER_CORE)
    nspec_region = NSPEC_OUTER_CORE
  case (IREGION_INNER_CORE)
    nspec_region = NSPEC_INNER_CORE
  end select

  ! checks if anything to do
  if (nspec_region == 0) return

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

          ! "factor_scale" consists of two scaling factors:
          !   physical dispersion factor: scales moduli from reference frequency to simulation (source) center frequency
          !   relaxed factor: scales moduli at center frequency to relaxed moduli (zero frequency value)
          ! Determine the "factor_scale" from tau_s, tau_e, central source frequency, and Q
          call get_attenuation_scale_factor(T_c_source, tau_e, tau_s, Q_mu, sf, sf_relaxed, f0_reference)

          factor_scale(i,j,k,ispec) = real(sf, kind=CUSTOM_REAL)
          factor_scale_relaxed(i,j,k,ispec) = real(sf_relaxed, kind=CUSTOM_REAL)
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
  double precision :: omsb1,omsb2
  integer :: i

  tauinv(:) = 0.d0
  where(abs(tau_s(:)) > 0.d0) tauinv(:) = 1.0d0 / tau_s(:)

  beta(:) = 1.0d0 - tau_e(:) * tauinv(:)     ! 1 - tau_e / tau_s

  ! factor to scale from relaxed to unrelaxed moduli: see e.g. Komatitsch & Tromp 1999, eq. (7)
  one_minus_sum_beta = 1.0d0
  do i = 1,N_SLS
     one_minus_sum_beta = one_minus_sum_beta - beta(i)
  enddo

  ! note: Komatitsch & Tromp 1999, eq. (7) (or Komatitsch 2002 eq. 10), defines the unrelaxed modulus as
  !          C_unrelaxed = C_relaxed * [ 1 - sum(1 - tau_eps/tau_sigma) ]
  !       with a scaling factor one_minus_sum_beta = [ 1 - sum(1 - tau_eps/tau_sigma) ]
  !
  !       this is different to Liu (1976) who mentions
  !          M_U = M_R / [ 1 - sum( (tau_eps - tau_sigma)/tau_eps ) ] = M_R / [ 1 - sum(1 - tau_sigma/tau_eps) ]
  !       with a scaling factor = 1 / [ 1 - sum(1 - tau_sigma/tau_eps) ]
  !
  ! to compare these two factors:
  if (.false.) then
    ! factor to scale from relaxed to unrelaxed moduli: see Komatitsch
    omsb1 = 1.0d0
    do i = 1,N_SLS
      omsb1  = omsb1 - (1.d0 - tau_e(i)/tau_s(i))
    enddo
    ! factor to scale from relaxed to unrelaxed moduli: see Liu
    omsb2 = 1.0d0
    do i = 1,N_SLS
      omsb2 = omsb2 - (1.d0 - tau_s(i)/tau_e(i))
    enddo
    print *,'debug: one_minus_sum_beta1,2 = ',omsb1,1.d0/omsb2
  endif
  ! results:
  ! Q ~    65: one_minus_sum_beta1,2 = 1.0648735275165400        1.0677815840121743   -> ratio = 0.997276
  ! Q ~   184: one_minus_sum_beta1,2 = 1.0233224020801992        1.0236877302481358   -> ratio = 0.999643
  ! Q ~ 19577: one_minus_sum_beta1,2 = 1.0002203256626838        1.0002203577661830   -> ratio = 0.999999967
  !
  ! -> for large enough Q, these factors are almost identical

  ! see Savage (2010), eq. (11) where the two terms with previous strain \epsilon' and current strain \epsilon
  ! both have a factor 2 \Delta M_i / \tau_{sigma_i}, where \Delta M_i is the modulus defect (M_unrelaxed - M_relaxed)_i
  !
  ! C_unrelaxed = C_relaxed * [ 1 - sum(1 - tau_eps/tau_sigma) ] = C_relaxed * [ 1 - sum(beta(:))]
  ! -> \Delta M_i =  (C_unrelaxed - C_relaxed)_i = C_relaxed * [1 - beta_i] - C_relaxed = C_relaxed * [1 - beta_i - 1]
  !               = C_relaxed * [ - beta_i] with beta_i = 1 - tau_e_i/tau_s_i
  ! (see also Komatitsch, 1999, eq. 8)

!ZN beware, here the expression differs from the strain used in memory variable equation (6) in D. Komatitsch and J. Tromp 1999,
!ZN here Brian Savage uses the engineering strain which are epsilon = 1/2*(grad U + (grad U)^T),
!ZN where U is the displacement vector and grad the gradient operator, i.e. there is a 1/2 factor difference between the two.
!ZN Both expressions are fine, but we need to keep in mind that if we have put the 1/2 factor there we need to remove it
!ZN from the expression in which we use the strain here in the code.
!ZN This is why here Brian Savage multiplies beta(:) * tauinv(:) by 2.0 to compensate for the 1/2 factor used before
  factor_common(:) = 2.0d0 * (- beta(:)) * tauinv(:)

  end subroutine get_attenuation_property_values

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_attenuation_scale_factor(T_c_source, tau_mu, tau_sigma, Q_mu, scale_factor, scale_factor_relaxed, f0_reference)

  use constants, only: ZERO,ONE,TWO,PI,TWO_PI,N_SLS,myrank
  use specfem_par, only: scale_t

  implicit none

  double precision,intent(in) :: T_c_source  ! center period
  double precision,intent(in) :: Q_mu
  double precision, dimension(N_SLS),intent(in) :: tau_mu, tau_sigma

  double precision,intent(out) :: scale_factor,scale_factor_relaxed

  ! model reference frequency
  double precision,intent(in) :: f0_reference

  ! local parameters
  double precision :: f_c_source, w_c_source,f_0_model
  double precision :: factor_scale_mu0,factor_scale_mu
  double precision :: a_val, b_val, denom
  double precision :: big_omega
  integer :: i

  ! compute central angular frequency of source (non dimensionalized)
  if (T_c_source > 0.d0) then
    f_c_source = ONE / T_c_source
  else
    f_c_source = 0.d0
  endif
  w_c_source = TWO_PI * f_c_source

  ! model reference frequency (e.g., PREM reference of 1 second)
  ! non-dimensionalizes frequency
  f_0_model = f0_reference * scale_t       ! original f_0_prem = ONE / ( ONE / scale_t) or f_reference = ATTENUATION_f0_REFERENCE

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

  !--- compute a, b and Omega parameters
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
  !
  ! note: by assuming that the input velocities are given at a reference frequency (f_0_model),
  !       we first shift them to the center frequency (f_c_source or angular frequency w_c_source) of the absorption-band
  !       and determine the corresponding relaxed moduli based on:
  !         M_relaxed = M(w_c_source) * B**2 / Omega
  !
  !       see Liu 1976 eq. (20), where here we use the expression for modulus, Omega = Omega(w_c_source), B = B(w_c_source).
  !       eq. 20 would be given for phase velocities v_p = sqrt( M(w_c_source)/rho) and v_e = sqrt( M_relaxed/rho)
  !
  !       This shifting to the center frequency of the absorption-band also guarantees that the relaxation times tau_**
  !       are optimal (for a good Q approximation within the band) and thus should give a reliable scaling to
  !       the relaxed modulus as well.
  !
  !       For the stiffness, we will need the moduli to be scaled to unrelaxed ones for the computations of the stress tensor:
  !         T = C_unrelaxed : grad(s) - sum( R_l)  (see e.g. Komatitsch, 1999, eq. 5)
  !
  !       A relaxed modulus can be seen as the modulus at zero frequency, unrelaxed modulus at infinite frequency (instantaneous)

  !       These factors together with typos and the debate about a missing 1/L factor leads to some confusion.
  !       It has been shown that as long as one uses a consistent set of equations, the 1/L factor expressions are equivalent
  !       to the ones used by Liu (1976) without it and implemented here.
  !
  !       Unfortunately, note that Liu (1976) starts with expressions for compliance J = 1/M rather than M, however turns
  !       to complex modulus M_c which again is a stress relaxation function sigma(t) = M_c(t) epsilon(t).
  !
  !       For Liu 1976, the quality factor Q is defined as tan delta = 1/Q = B/A
  !
  !       We further assume the reference velocities are given at a certain reference frequency - and not as relaxed or unrelaxed
  !       velocities, thus velocities need to be shifted to a absorption-band and the relaxed/unrelaxed moduli must be determined.
  !       This is a convention and could also be defined in different ways, see e.g. Carcione (1993) who defines
  !       the elastic model given by relaxed values.
  !
  if (big_omega /= 0.d0) then
    factor_scale_mu = b_val * b_val / (TWO * big_omega)
  else
    factor_scale_mu = ONE
  endif

  !--- total factor by which to scale mu0
  scale_factor = factor_scale_mu * factor_scale_mu0

  ! factor to get relaxed modulus from shifted value
  scale_factor_relaxed = factor_scale_mu

  ! note:
  !  * if reference frequency and center frequency are the same, i.e., f_0_model == f_c_source, then factor_scale_mu0 == 1
  !  * for f_c_source == 0, i.e. at zero frequency, the scaling to relaxed factor_scale_mu should be 1
  !  * Liu 1976 defines tan delta = 1/Q = B/A -> Q = A/B gives the approximated Q value for the given SLS
  !print *,'debug: factor scale mu0,mu = ',factor_scale_mu0,factor_scale_mu,'approximated Q = A/B',a_val/b_val

  !--- check that the correction factor is close to one
  ! note: for very coarse global simulations, the reference frequency (usually 1Hz for global models)
  !       and center frequency of the simulation can be far apart.
  !       the scale_factor limiting range here is somewhat set arbitrarily, and enlarged to [0.5,1.5]
  !       to allow for very coarse global simulation tests with sponge attenuation.
  if (scale_factor < 0.5d0 .or. scale_factor > 1.5d0) then
    print *,'Error: incorrect scale factor: ', scale_factor
    print *,'  scale factor: ', scale_factor,' should be between 0.5 to 1.5'
    print *,'  factor scale_mu = ',factor_scale_mu,' factor scale_mu0 = ',factor_scale_mu0
    print *,'  Q value = ',Q_mu
    print *,'  central period = ',T_c_source * scale_t,' frequency = ',f_c_source / scale_t
    print *,'  model frequency = ',f_0_model / scale_t
    print *,'Please check your reference frequency ATTENUATION_f0_REFERENCE in file constants.h'
    call exit_MPI(myrank,'incorrect correction factor in attenuation model')
  endif

  end subroutine get_attenuation_scale_factor

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_attenuation_memory_values(tau_s, deltat, alphaval, betaval, gammaval)

  use constants

  implicit none

  double precision, dimension(N_SLS), intent(in) :: tau_s
  double precision, dimension(N_SLS), intent(out) :: alphaval, betaval,gammaval
  real(kind=CUSTOM_REAL), intent(in) :: deltat

  ! local parameters
  double precision, dimension(N_SLS) :: tauinv

! daniel todo: not sure why here we take the negative 1/tau_s_i, Brian's eq. 11 uses positive signs.

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

  subroutine bcast_attenuation_model_3D(one_minus_sum_beta,factor_common,factor_scale,factor_scale_relaxed, &
                                        tau_s,vnspec)

  use constants_solver
  use shared_parameters, only: ATT_F_C_SOURCE

  implicit none

  integer :: vnspec

  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec) :: one_minus_sum_beta, factor_scale
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,N_SLS,vnspec) :: factor_common
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec) :: factor_scale_relaxed
  double precision, dimension(N_SLS) :: tau_s

! note: the size(..) function returns either integer(kind=4) or integer(kind=8)
!       depending on compiler flags (-mcmedium), thus adding a kind argument to have integer(kind=4) output

  call bcast_all_cr_for_database(one_minus_sum_beta(1,1,1,1), size(one_minus_sum_beta,kind=4))
  call bcast_all_cr_for_database(factor_scale(1,1,1,1), size(factor_scale,kind=4))
  call bcast_all_cr_for_database(factor_scale_relaxed(1,1,1,1), size(factor_scale_relaxed,kind=4))
  call bcast_all_cr_for_database(factor_common(1,1,1,1,1), size(factor_common,kind=4))

  call bcast_all_dp_for_database(tau_s(1), N_SLS)
  call bcast_all_dp_for_database(ATT_F_C_SOURCE, 1)

  end subroutine bcast_attenuation_model_3D
