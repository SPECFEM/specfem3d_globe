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
  real(kind=CUSTOM_REAL) :: mul
  integer :: ispec,i,j,k,ier

  ! checks if attenuation is on and anything to do
  if (.not. ATTENUATION_VAL ) return

  ! get and store PREM attenuation model
  if (myrank == 0) then
    write(IMAIN,*) "preparing attenuation"
    call flush_IMAIN()
  endif

  ! allocates arrays
  allocate(one_minus_sum_beta_crust_mantle(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT4_VAL), &
           factor_scale_crust_mantle(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT4_VAL), &
           factor_common_crust_mantle(ATT1_VAL,ATT2_VAL,ATT3_VAL,N_SLS,ATT4_VAL),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays one_minus_sum_beta_crust_mantle,..'

  allocate(one_minus_sum_beta_inner_core(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT5_VAL), &
           factor_scale_inner_core(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT5_VAL), &
           factor_common_inner_core(ATT1_VAL,ATT2_VAL,ATT3_VAL,N_SLS,ATT5_VAL),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays one_minus_sum_beta_inner_core,..'

  ! reads in attenuation values
  ! CRUST_MANTLE ATTENUATION
  call get_attenuation_model_3D(IREGION_CRUST_MANTLE, &
                                one_minus_sum_beta_crust_mantle, &
                                factor_common_crust_mantle, &
                                factor_scale_crust_mantle,tau_sigma_dble, &
                                NSPEC_CRUST_MANTLE)
  call bcast_attenuation_model_3D(one_minus_sum_beta_crust_mantle, &
                                  factor_common_crust_mantle, &
                                  factor_scale_crust_mantle, &
                                  tau_sigma_dble, &
                                  NSPEC_CRUST_MANTLE)

  ! INNER_CORE ATTENUATION
  call get_attenuation_model_3D(IREGION_INNER_CORE, &
                                one_minus_sum_beta_inner_core, &
                                factor_common_inner_core, &
                                factor_scale_inner_core,tau_sigma_dble, &
                                NSPEC_INNER_CORE)
  call bcast_attenuation_model_3D(one_minus_sum_beta_inner_core, &
                                  factor_common_inner_core, &
                                  factor_scale_inner_core, &
                                  tau_sigma_dble, &
                                  NSPEC_INNER_CORE)

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

  ! rescale in crust and mantle
  do ispec = 1,NSPEC_CRUST_MANTLE
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
            scale_factor = factor_scale_crust_mantle(i,j,k,ispec)
          else
            scale_factor = factor_scale_crust_mantle(1,1,1,ispec)
          endif

          if (ANISOTROPIC_3D_MANTLE_VAL) then
            scale_factor_minus_one = scale_factor - 1.d0

            mul = c44store_crust_mantle(i,j,k,ispec)

            c11store_crust_mantle(i,j,k,ispec) = c11store_crust_mantle(i,j,k,ispec) &
                    + FOUR_THIRDS * scale_factor_minus_one * mul
            c12store_crust_mantle(i,j,k,ispec) = c12store_crust_mantle(i,j,k,ispec) &
                    - TWO_THIRDS * scale_factor_minus_one * mul
            c13store_crust_mantle(i,j,k,ispec) = c13store_crust_mantle(i,j,k,ispec) &
                    - TWO_THIRDS * scale_factor_minus_one * mul
            c22store_crust_mantle(i,j,k,ispec) = c22store_crust_mantle(i,j,k,ispec) &
                    + FOUR_THIRDS * scale_factor_minus_one * mul
            c23store_crust_mantle(i,j,k,ispec) = c23store_crust_mantle(i,j,k,ispec) &
                    - TWO_THIRDS * scale_factor_minus_one * mul
            c33store_crust_mantle(i,j,k,ispec) = c33store_crust_mantle(i,j,k,ispec) &
                    + FOUR_THIRDS * scale_factor_minus_one * mul
            c44store_crust_mantle(i,j,k,ispec) = c44store_crust_mantle(i,j,k,ispec) &
                    + scale_factor_minus_one * mul
            c55store_crust_mantle(i,j,k,ispec) = c55store_crust_mantle(i,j,k,ispec) &
                    + scale_factor_minus_one * mul
            c66store_crust_mantle(i,j,k,ispec) = c66store_crust_mantle(i,j,k,ispec) &
                    + scale_factor_minus_one * mul
          else
            if (MOVIE_VOLUME .and. SIMULATION_TYPE == 3) then
              ! store the original value of \mu to compute \mu*\eps
              muvstore_crust_mantle_3dmovie(i,j,k,ispec) = muvstore_crust_mantle(i,j,k,ispec)
            endif

            muvstore_crust_mantle(i,j,k,ispec) = muvstore_crust_mantle(i,j,k,ispec) * scale_factor

            ! scales transverse isotropic values for mu_h
            if (ispec_is_tiso_crust_mantle(ispec)) then
              muhstore_crust_mantle(i,j,k,ispec) = muhstore_crust_mantle(i,j,k,ispec) * scale_factor
            endif
          endif

        enddo
      enddo
    enddo
  enddo ! enddo CRUST MANTLE

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
                    + FOUR_THIRDS * scale_factor_minus_one * mul
            c12store_inner_core(i,j,k,ispec) = c12store_inner_core(i,j,k,ispec) &
                    - TWO_THIRDS * scale_factor_minus_one * mul
            c13store_inner_core(i,j,k,ispec) = c13store_inner_core(i,j,k,ispec) &
                    - TWO_THIRDS * scale_factor_minus_one * mul
            c33store_inner_core(i,j,k,ispec) = c33store_inner_core(i,j,k,ispec) &
                    + FOUR_THIRDS * scale_factor_minus_one * mul
            c44store_inner_core(i,j,k,ispec) = c44store_inner_core(i,j,k,ispec) &
                    + scale_factor_minus_one * mul
          endif

          muvstore_inner_core(i,j,k,ispec) = muvstore_inner_core(i,j,k,ispec) * scale_factor

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
    tau_sigma_CUSTOM_REAL(:) = real(tau_sigma_dble(:), kind=CUSTOM_REAL)
  endif

  if (UNDO_ATTENUATION) then
   b_alphaval = alphaval
   b_betaval = betaval
   b_gammaval = gammaval
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  attenuation period range min/max: ",MIN_ATTENUATION_PERIOD,'/',MAX_ATTENUATION_PERIOD,' (s)'
    write(IMAIN,*) "  ATTENUATION_1D_WITH_3D_STORAGE  : ",ATTENUATION_1D_WITH_3D_STORAGE_VAL
    write(IMAIN,*) "  ATTENUATION_3D                  : ",ATTENUATION_3D_VAL
    call flush_IMAIN()
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_attenuation

