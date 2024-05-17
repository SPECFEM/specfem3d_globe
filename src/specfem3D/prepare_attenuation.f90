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



  subroutine prepare_attenuation()

  ! precomputes attenuation factors

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_movie, only: muvstore_crust_mantle_3dmovie

  implicit none

  ! local parameters
  double precision, dimension(N_SLS) :: alphaval_dble, betaval_dble, gammaval_dble
  double precision :: scale_factor,scale_factor_minus_one
  real(kind=CUSTOM_REAL) :: mul,muvl,muhl,eta_aniso
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: factor_scale_relaxed_crust_mantle,factor_scale_relaxed_inner_core
  ! aniso element
  !real(kind=CUSTOM_REAL) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
  !                          c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  double precision :: g11,g12,g13,g14,g15,g16,g22,g23,g24,g25,g26, &
                      g33,g34,g35,g36,g44,g45,g46,g55,g56,g66
  double precision :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66
  double precision :: phi_dble,theta_dble
  double precision :: A_dble,F_dble,L_dble,N_dble!,C_dble
  double precision :: f0_reference

  integer :: ispec,i,j,k,ier,iglob

  ! needs these allocated for subroutine calls
  allocate(factor_common_crust_mantle(ATT1_VAL,ATT2_VAL,ATT3_VAL,N_SLS,ATT4_VAL),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays factor_common_crust_mantle'
  factor_common_crust_mantle(:,:,:,:,:) = 0.0_CUSTOM_REAL

  allocate(factor_common_inner_core(ATT1_VAL,ATT2_VAL,ATT3_VAL,N_SLS,ATT5_VAL),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays factor_common_inner_core'
  factor_common_inner_core(:,:,:,:,:) = 0.0_CUSTOM_REAL

  ! checks if attenuation is on and anything to do further
  if (.not. ATTENUATION_VAL ) return

  ! get and store PREM attenuation model
  if (myrank == 0) then
    write(IMAIN,*) "preparing attenuation"
    write(IMAIN,*) "  The code uses a constant Q quality factor, but approximated"
    write(IMAIN,*) "  based on a series of Zener standard linear solids (SLS)."
    write(IMAIN,*) "  Approximation is performed in the following frequency band:"
    write(IMAIN,*)
    write(IMAIN,*) "  number of SLS bodies: ",N_SLS
    write(IMAIN,*) "  partial attenuation, physical dispersion only: ",PARTIAL_PHYS_DISPERSION_ONLY_VAL
    write(IMAIN,*)
    write(IMAIN,*) "  Reference frequency of anelastic model (Hz): ",sngl(ATTENUATION_f0_REFERENCE)
    write(IMAIN,*) "                                   period (s): ",sngl(1.0/ATTENUATION_f0_REFERENCE)
    call flush_IMAIN()
  endif

  ! allocates arrays
  allocate(one_minus_sum_beta_crust_mantle(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT4_VAL), &
           factor_scale_crust_mantle(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT4_VAL),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays one_minus_sum_beta_crust_mantle,..'
  one_minus_sum_beta_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  factor_scale_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL

  allocate(one_minus_sum_beta_inner_core(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT5_VAL), &
           factor_scale_inner_core(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT5_VAL),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays one_minus_sum_beta_inner_core,..'
  one_minus_sum_beta_inner_core(:,:,:,:) = 0.0_CUSTOM_REAL
  factor_scale_inner_core(:,:,:,:) = 0.0_CUSTOM_REAL

  ! used for reading in and shifted value output
  allocate(factor_scale_relaxed_crust_mantle(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT4_VAL),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays factor_scale_relaxed_crust_mantle'
  factor_scale_relaxed_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL

  allocate(factor_scale_relaxed_inner_core(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT5_VAL),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays factor_scale_relaxed_inner_core'
  factor_scale_relaxed_inner_core(:,:,:,:) = 0.0_CUSTOM_REAL

  ! reference frequency for model velocities (e.g., PREM at 1 Hz)
  f0_reference = ATTENUATION_f0_REFERENCE

  ! gets attenuation values
  ! CRUST_MANTLE ATTENUATION
  if (NSPEC_CRUST_MANTLE > 0) then
    call get_attenuation_model_3D(IREGION_CRUST_MANTLE, &
                                  one_minus_sum_beta_crust_mantle, &
                                  factor_common_crust_mantle, &
                                  factor_scale_crust_mantle, &
                                  factor_scale_relaxed_crust_mantle, &
                                  tau_sigma_dble, &
                                  NSPEC_CRUST_MANTLE,f0_reference)

    call bcast_attenuation_model_3D(one_minus_sum_beta_crust_mantle, &
                                    factor_common_crust_mantle, &
                                    factor_scale_crust_mantle, &
                                    factor_scale_relaxed_crust_mantle, &
                                    tau_sigma_dble, &
                                    NSPEC_CRUST_MANTLE)
  endif

  ! INNER_CORE ATTENUATION
  if (NSPEC_INNER_CORE > 0) then
    call get_attenuation_model_3D(IREGION_INNER_CORE, &
                                  one_minus_sum_beta_inner_core, &
                                  factor_common_inner_core, &
                                  factor_scale_inner_core, &
                                  factor_scale_relaxed_inner_core, &
                                  tau_sigma_dble, &
                                  NSPEC_INNER_CORE,f0_reference)

    call bcast_attenuation_model_3D(one_minus_sum_beta_inner_core, &
                                    factor_common_inner_core, &
                                    factor_scale_inner_core, &
                                    factor_scale_relaxed_inner_core, &
                                    tau_sigma_dble, &
                                    NSPEC_INNER_CORE)
  endif

  ! debug
  !if (myrank == 0) print *,'debug: original moduli muv',muvstore_crust_mantle(1,1,1,1000),muvstore_crust_mantle(3,3,3,3000)
  !if (myrank == 0) print *,'debug: original moduli c11',c11store_crust_mantle(1,1,1,1000),c11store_crust_mantle(3,3,3,3000)
  !if (myrank == 0) print *,'debug: original moduli c44',c44store_crust_mantle(1,1,1,1000),c44store_crust_mantle(3,3,3,3000)

  ! note: the factor "scale_factor" below includes both the scaling to the center frequency
  !       as well as the conversion to the relaxed moduli. the final result when multiplying with this factor
  !       will be the relaxed modulus value (given the SLS approximation).
  !       we will need the relaxed and unrelaxed values to determine the modulus defect (necessary to propagate
  !       the memory variables) and for the stress computation (which uses the unrelaxed values).
  !
  ! physical dispersion: scales moduli from reference frequency to simulation (source) center frequency
  !
  ! if attenuation is on, shift PREM to right frequency
  ! rescale mu in PREM to average frequency for attenuation
  !
  ! the formulas to implement the scaling can be found for instance in
  ! Liu, H. P., Anderson, D. L. and Kanamori, H., Velocity dispersion due to
  ! anelasticity: implications for seismology and mantle composition,
  ! Geophys. J. R. Astron. Soc., vol. 47, pp. 41-58 (1976)
  !
  ! and in Aki, K. and Richards, P. G., Quantitative seismology, theory and methods,
  ! W. H. Freeman, (1980), second edition, sections 5.5 and 5.5.2, eq. (5.81) p. 170.
  !
  ! Beware that in the book of Aki and Richards eq. (5.81) is given for velocities
  ! while we need an equation for "mu" and thus we have an additional factor of 2
  ! in the scaling factor below and in equation (49) of Komatitsch and Tromp, Geophys. J. Int. (2002) 149, 390-412,
  ! because "mu" is related to the square of velocity.
  !
  ! mu(omega_c) = mu(omega_0)[ 1 + 2/(pi Q_mu) ln(omega_c / omega_0) ]
  !
  ! once this shifted modulus value is given, we approximate the relaxed value given the standard linear solids
  ! in routine get_attenuation_scale_factor() by:
  !
  !  M_relaxed = M(w_c_source) * B**2 / Omega
  !
  ! see Liu 1976 eq. (20), where here we use the expression for modulus, Omega = Omega(w_c_source), B = B(w_c_source).
  ! Liu's eq. (20) would be given for phase velocities v_p = sqrt( M(w_c_source)/rho) and v_e = sqrt( M_relaxed/rho).
  !
  ! rescale in crust and mantle
  do ispec = 1,NSPEC_CRUST_MANTLE
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          ! note: only shear attenuation is implemented
          !       thus, only shear moduli are shifted to center frequency
          if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
            scale_factor = factor_scale_crust_mantle(i,j,k,ispec)
          else
            scale_factor = factor_scale_crust_mantle(1,1,1,ispec)
          endif

          if (ANISOTROPIC_3D_MANTLE_VAL) then
            ! anisotropic element
            scale_factor_minus_one = scale_factor - 1.d0

            ! original routine:
            ! the following shifts moduli to the center frequency.
            ! note: the scaling here assumes directly c44 = muv without rotation and becomes inconsistent with tiso simulations
            !  c11 = c11store_crust_mantle(i,j,k,ispec)
            !  c12 = c12store_crust_mantle(i,j,k,ispec)
            !  c13 = c13store_crust_mantle(i,j,k,ispec)
            !  c22 = c22store_crust_mantle(i,j,k,ispec)
            !  c23 = c23store_crust_mantle(i,j,k,ispec)
            !  c33 = c33store_crust_mantle(i,j,k,ispec)
            !  c44 = c44store_crust_mantle(i,j,k,ispec)
            !  c55 = c55store_crust_mantle(i,j,k,ispec)
            !  c66 = c66store_crust_mantle(i,j,k,ispec)
            !
            !  mul = c44store_crust_mantle(i,j,k,ispec)
            !
            !  c11store_crust_mantle(i,j,k,ispec) = c11 + FOUR_THIRDS * scale_factor_minus_one * mul
            !  c12store_crust_mantle(i,j,k,ispec) = c12 - TWO_THIRDS * scale_factor_minus_one * mul
            !  c13store_crust_mantle(i,j,k,ispec) = c13 - TWO_THIRDS * scale_factor_minus_one * mul
            !  c22store_crust_mantle(i,j,k,ispec) = c22 + FOUR_THIRDS * scale_factor_minus_one * mul
            !  c23store_crust_mantle(i,j,k,ispec) = c23 - TWO_THIRDS * scale_factor_minus_one * mul
            !  c33store_crust_mantle(i,j,k,ispec) = c33 + FOUR_THIRDS * scale_factor_minus_one * mul
            !  c44store_crust_mantle(i,j,k,ispec) = c44 + scale_factor_minus_one * mul
            !  c55store_crust_mantle(i,j,k,ispec) = c55 + scale_factor_minus_one * mul
            !  c66store_crust_mantle(i,j,k,ispec) = c66 + scale_factor_minus_one * mul
            !
            !  ! for attenuation
            !  muvstore_crust_mantle(i,j,k,ispec) = mul * scale_factor

            ! new routine to be consistent with tiso attenuation
            !
            ! local position (d_ij given in radial direction)
            ! only in case needed for rotation
            iglob = ibool_crust_mantle(i,j,k,ispec)
            theta_dble = rstore_crust_mantle(2,iglob)
            phi_dble = rstore_crust_mantle(3,iglob)
            call reduce(theta_dble,phi_dble)

            g11 = c11store_crust_mantle(i,j,k,ispec)
            g12 = c12store_crust_mantle(i,j,k,ispec)
            g13 = c13store_crust_mantle(i,j,k,ispec)
            g14 = c14store_crust_mantle(i,j,k,ispec)
            g15 = c15store_crust_mantle(i,j,k,ispec)
            g16 = c16store_crust_mantle(i,j,k,ispec)
            g22 = c22store_crust_mantle(i,j,k,ispec)
            g23 = c23store_crust_mantle(i,j,k,ispec)
            g24 = c24store_crust_mantle(i,j,k,ispec)
            g25 = c25store_crust_mantle(i,j,k,ispec)
            g26 = c26store_crust_mantle(i,j,k,ispec)
            g33 = c33store_crust_mantle(i,j,k,ispec)
            g34 = c34store_crust_mantle(i,j,k,ispec)
            g35 = c35store_crust_mantle(i,j,k,ispec)
            g36 = c36store_crust_mantle(i,j,k,ispec)
            g44 = c44store_crust_mantle(i,j,k,ispec)
            g45 = c45store_crust_mantle(i,j,k,ispec)
            g46 = c46store_crust_mantle(i,j,k,ispec)
            g55 = c55store_crust_mantle(i,j,k,ispec)
            g56 = c56store_crust_mantle(i,j,k,ispec)
            g66 = c66store_crust_mantle(i,j,k,ispec)

            call rotate_tensor_global_to_radial(theta_dble,phi_dble, &
                                                d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                                d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                                                g11,g12,g13,g14,g15,g16,g22,g23,g24,g25,g26, &
                                                g33,g34,g35,g36,g44,g45,g46,g55,g56,g66)

            ! new: shifts moduli by separating muv and muh factors.
            !      still needs rotations to rotate back and forth from SPECFEM global axis to a radial symmetry axis
            !      since this shift assumes a radial symmetry
            A_dble = 0.125d0 * (3.d0 * d11 + 3.d0 * d22 + 2.d0 * d12 + 4.d0 * d66)
            !unused: C_dble = d33
            N_dble = 0.125d0 * (d11 + d22 - 2.d0 * d12 + 4.d0 * d66)
            L_dble = 0.5d0 * (d44 + d55)
            F_dble = 0.5d0 * (d13 + d23)

            eta_aniso = real(F_dble / (A_dble - 2.d0*L_dble),kind=CUSTOM_REAL)   ! eta = F / (A-2L)

            muvl = real(L_dble * scale_factor_minus_one,kind=CUSTOM_REAL)     ! c44 - > L - > muv
            muhl = real(N_dble * scale_factor_minus_one,kind=CUSTOM_REAL)     ! c66 - > N - > muh

            d11 = d11 + FOUR_THIRDS * muhl ! * minus_sum_beta * mul
            d12 = d12 - TWO_THIRDS * muhl
            d13 = d13 + eta_aniso * (FOUR_THIRDS * muhl - 2.d0*muvl)
            d22 = d22 + FOUR_THIRDS * muhl
            d23 = d23 + eta_aniso * (FOUR_THIRDS * muhl - 2.d0*muvl)
            d33 = d33 + FOUR_THIRDS * muvl
            d44 = d44 + muvl
            d55 = d55 + muvl
            d66 = d66 + muhl

            ! debug
            !if (myrank == 0 .and. ispec == 1000 .and. i == 1 .and. j == 1 .and. k == 1) &
            !  print *,'debug: original moduli scaling A,N,L,F,eta',A_dble,N_dble,L_dble,F_dble,eta_aniso,'mu',muvl,muhl

            ! rotates to global reference system
            call rotate_tensor_radial_to_global(theta_dble,phi_dble, &
                                                d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                                d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                                                g11,g12,g13,g14,g15,g16,g22,g23,g24,g25,g26, &
                                                g33,g34,g35,g36,g44,g45,g46,g55,g56,g66)

            ! stores unrelaxed factors
            c11store_crust_mantle(i,j,k,ispec) = real(g11,kind=CUSTOM_REAL)
            c12store_crust_mantle(i,j,k,ispec) = real(g12,kind=CUSTOM_REAL)
            c13store_crust_mantle(i,j,k,ispec) = real(g13,kind=CUSTOM_REAL)
            c14store_crust_mantle(i,j,k,ispec) = real(g14,kind=CUSTOM_REAL)
            c15store_crust_mantle(i,j,k,ispec) = real(g15,kind=CUSTOM_REAL)
            c16store_crust_mantle(i,j,k,ispec) = real(g16,kind=CUSTOM_REAL)
            c22store_crust_mantle(i,j,k,ispec) = real(g22,kind=CUSTOM_REAL)
            c23store_crust_mantle(i,j,k,ispec) = real(g23,kind=CUSTOM_REAL)
            c24store_crust_mantle(i,j,k,ispec) = real(g24,kind=CUSTOM_REAL)
            c25store_crust_mantle(i,j,k,ispec) = real(g25,kind=CUSTOM_REAL)
            c26store_crust_mantle(i,j,k,ispec) = real(g26,kind=CUSTOM_REAL)
            c33store_crust_mantle(i,j,k,ispec) = real(g33,kind=CUSTOM_REAL)
            c34store_crust_mantle(i,j,k,ispec) = real(g34,kind=CUSTOM_REAL)
            c35store_crust_mantle(i,j,k,ispec) = real(g35,kind=CUSTOM_REAL)
            c36store_crust_mantle(i,j,k,ispec) = real(g36,kind=CUSTOM_REAL)
            c44store_crust_mantle(i,j,k,ispec) = real(g44,kind=CUSTOM_REAL)
            c45store_crust_mantle(i,j,k,ispec) = real(g45,kind=CUSTOM_REAL)
            c46store_crust_mantle(i,j,k,ispec) = real(g46,kind=CUSTOM_REAL)
            c55store_crust_mantle(i,j,k,ispec) = real(g55,kind=CUSTOM_REAL)
            c56store_crust_mantle(i,j,k,ispec) = real(g56,kind=CUSTOM_REAL)
            c66store_crust_mantle(i,j,k,ispec) = real(g66,kind=CUSTOM_REAL)

            ! for attenuation
            muvstore_crust_mantle(i,j,k,ispec) = real(L_dble * scale_factor,kind=CUSTOM_REAL)

          else
            ! isotropic or transverse isotropic element
            if (MOVIE_VOLUME .and. SIMULATION_TYPE == 3) then
              ! store the original value of \mu to compute \mu*\eps
              muvstore_crust_mantle_3dmovie(i,j,k,ispec) = muvstore_crust_mantle(i,j,k,ispec)
            endif

            muvstore_crust_mantle(i,j,k,ispec) = real(muvstore_crust_mantle(i,j,k,ispec) * scale_factor,kind=CUSTOM_REAL)

            ! scales transverse isotropic values for mu_h
            if (ispec_is_tiso_crust_mantle(ispec)) then
              muhstore_crust_mantle(i,j,k,ispec) = real(muhstore_crust_mantle(i,j,k,ispec) * scale_factor,kind=CUSTOM_REAL)
            endif
          endif

        enddo
      enddo
    enddo
  enddo ! enddo CRUST MANTLE

  ! note: no scaling for outer core arrays.
  !       we only implement shear attenuation so far, can be neglected in the outer core where vs == 0.

  ! rescale in inner core
  do ispec = 1,NSPEC_INNER_CORE
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
            scale_factor = factor_scale_inner_core(i,j,k,ispec)
          else
            scale_factor = factor_scale_inner_core(1,1,1,ispec)
          endif

          if (ANISOTROPIC_INNER_CORE_VAL) then
            scale_factor_minus_one = scale_factor - 1.d0

            mul = c44store_inner_core(i,j,k,ispec)
            c11store_inner_core(i,j,k,ispec) = c11store_inner_core(i,j,k,ispec) &
                    + real(FOUR_THIRDS * scale_factor_minus_one * mul,kind=CUSTOM_REAL)
            c12store_inner_core(i,j,k,ispec) = c12store_inner_core(i,j,k,ispec) &
                    - real(TWO_THIRDS * scale_factor_minus_one * mul,kind=CUSTOM_REAL)
            c13store_inner_core(i,j,k,ispec) = c13store_inner_core(i,j,k,ispec) &
                    - real(TWO_THIRDS * scale_factor_minus_one * mul,kind=CUSTOM_REAL)
            c33store_inner_core(i,j,k,ispec) = c33store_inner_core(i,j,k,ispec) &
                    + real(FOUR_THIRDS * scale_factor_minus_one * mul,kind=CUSTOM_REAL)
            c44store_inner_core(i,j,k,ispec) = c44store_inner_core(i,j,k,ispec) &
                    + real(scale_factor_minus_one * mul,kind=CUSTOM_REAL)
            ! for attenuation
            muvstore_inner_core(i,j,k,ispec) = real(mul * scale_factor,kind=CUSTOM_REAL)
          else
            muvstore_inner_core(i,j,k,ispec) = real(muvstore_inner_core(i,j,k,ispec) * scale_factor,kind=CUSTOM_REAL)
          endif
        enddo
      enddo
    enddo
  enddo ! enddo INNER CORE

  ! precompute Runge-Kutta coefficients
  call get_attenuation_memory_values(tau_sigma_dble, deltat, alphaval_dble, betaval_dble, gammaval_dble)
  alphaval = real(alphaval_dble, kind=CUSTOM_REAL)
  betaval  = real(betaval_dble, kind=CUSTOM_REAL)
  gammaval = real(gammaval_dble, kind=CUSTOM_REAL)

  if (SIMULATION_TYPE == 3) then
    call get_attenuation_memory_values(tau_sigma_dble, b_deltat, alphaval_dble, betaval_dble, gammaval_dble)
    b_alphaval = real(alphaval_dble, kind=CUSTOM_REAL)
    b_betaval  = real(betaval_dble, kind=CUSTOM_REAL)
    b_gammaval = real(gammaval_dble, kind=CUSTOM_REAL)
  endif

  if (USE_LDDRK) then
    ! inverts tau_sigma as we only need 1/tau_sigma factors
    tau_sigmainv_CUSTOM_REAL(:) = real(1.d0/tau_sigma_dble(:), kind=CUSTOM_REAL)
  endif

  if (UNDO_ATTENUATION) then
    b_alphaval = alphaval
    b_betaval = betaval
    b_gammaval = gammaval
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  Attenuation frequency band min/max (Hz):",sngl(1.0/MAX_ATTENUATION_PERIOD), &
                                                            '/',sngl(1.0/MIN_ATTENUATION_PERIOD)
    write(IMAIN,*) "              period band    min/max (s) :",sngl(MIN_ATTENUATION_PERIOD), &
                                                            '/',sngl(MAX_ATTENUATION_PERIOD)
    write(IMAIN,*) "  Logarithmic center frequency (Hz):",sngl(ATT_F_C_SOURCE)
    write(IMAIN,*) "                     period     (s):",sngl(1.0/ATT_F_C_SOURCE)
    write(IMAIN,*)
    write(IMAIN,*) "  using shear attenuation Q_mu"
    write(IMAIN,*)
    write(IMAIN,*) "  ATTENUATION_1D_WITH_3D_STORAGE  : ",ATTENUATION_1D_WITH_3D_STORAGE_VAL
    write(IMAIN,*) "  ATTENUATION_3D                  : ",ATTENUATION_3D_VAL
    call flush_IMAIN()
  endif

  if (ATTENUATION_SAVE_MODEL_AT_SHIFTED_CENTER_FREQ) then
    ! user output
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  saving velocity model at shifted center frequency...'
      write(IMAIN,*) '  new reference frequency (ATTENUATION_f0_REFERENCE) for anelastic model: ',ATT_F_C_SOURCE,'(Hz)'
      if (TRANSVERSE_ISOTROPY) then
        write(IMAIN,*) '  shear attenuation affects: (vpv,vph,vsv,vsh)'
      else
        write(IMAIN,*) '  shear attenuation affects: (vp,vs)'
      endif
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! outputs model files
    if (ADIOS_FOR_SOLVER_MESHFILES) then
      ! adios file output
      call save_forward_model_at_shifted_frequency_adios(factor_scale_relaxed_crust_mantle,factor_scale_relaxed_inner_core)
    else
      ! outputs model files in binary format
      call save_forward_model_at_shifted_frequency(factor_scale_relaxed_crust_mantle,factor_scale_relaxed_inner_core)
    endif
  endif

  ! frees memory
  deallocate(factor_scale_relaxed_crust_mantle,factor_scale_relaxed_inner_core)

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_attenuation

