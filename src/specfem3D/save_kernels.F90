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

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"



  subroutine save_kernels()

  use constants_solver, only: SAVE_KERNELS_BOUNDARY,SAVE_KERNELS_OC,SAVE_KERNELS_IC

  use specfem_par, only: NOISE_TOMOGRAPHY,SIMULATION_TYPE,nrec_local, &
    APPROXIMATE_HESS_KL,ADIOS_FOR_KERNELS

  use specfem_par_innercore, only: rhostore_inner_core,muvstore_inner_core,kappavstore_inner_core, &
    rho_kl_inner_core,alpha_kl_inner_core,beta_kl_inner_core

  use specfem_par_outercore, only: rhostore_outer_core,kappavstore_outer_core,rho_kl_outer_core,alpha_kl_outer_core

  use manager_adios

  implicit none

  ! Open an handler to the ADIOS file in which kernel variables are written.
  if (ADIOS_FOR_KERNELS) then
    if ((SIMULATION_TYPE == 3) .or. (SIMULATION_TYPE == 2 .and. nrec_local > 0)) then
      call define_kernel_adios_variables()
    endif
  endif

  ! dump kernel arrays
  if (SIMULATION_TYPE == 3) then

    ! restores original reference moduli (before shifting and unrelaxing)
    call restore_original_moduli()

    ! crust mantle
    call save_kernels_crust_mantle()

    ! noise strength kernel
    if (NOISE_TOMOGRAPHY == 3) then
      call save_kernels_strength_noise()
    endif

    ! outer core
    if (SAVE_KERNELS_OC) then
      call save_kernels_outer_core(rhostore_outer_core,kappavstore_outer_core,rho_kl_outer_core,alpha_kl_outer_core)
    endif

    ! inner core
    if (SAVE_KERNELS_IC) then
      call save_kernels_inner_core(rhostore_inner_core,muvstore_inner_core,kappavstore_inner_core, &
                                   rho_kl_inner_core,alpha_kl_inner_core,beta_kl_inner_core)
    endif

    ! boundary kernel
    if (SAVE_KERNELS_BOUNDARY) then
      call save_kernels_boundary_kl()
    endif

    ! approximate Hessian
    if (APPROXIMATE_HESS_KL) then
      call save_kernels_Hessian()
    endif
  endif

  ! save source derivatives for adjoint simulations
  if (SIMULATION_TYPE == 2 .and. nrec_local > 0) then
    call save_kernels_source_derivatives()
  endif

  ! Write ADIOS defined variables to disk.
  if (ADIOS_FOR_KERNELS) then
    if ((SIMULATION_TYPE == 3) .or. (SIMULATION_TYPE == 2 .and. nrec_local > 0)) then
      call close_kernel_adios_file()
    endif
  endif

  end subroutine save_kernels

!
!-------------------------------------------------------------------------------------------------
!

  subroutine restore_original_moduli()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! checks if anything to do
  ! only needs to restore original moduli when attenuation was used in prepare_attenuation.f90 to change moduli
  if (.not. ATTENUATION_VAL) return

  ! debug
  !if (myrank == 0) print *,'debug: unrestored moduli muv',muvstore_crust_mantle(1,1,1,1000),muvstore_crust_mantle(3,3,3,3000)
  !if (myrank == 0) print *,'debug: unrestored moduli c11',c11store_crust_mantle(1,1,1,1000),c11store_crust_mantle(3,3,3,3000)
  !if (myrank == 0) print *,'debug: unrestored moduli c44',c44store_crust_mantle(1,1,1,1000),c44store_crust_mantle(3,3,3,3000)


  ! note: we could save the original moduli in additional **store arrays.
  !       however, this will almost double the memory requirements for the kernel simulations.
  !       we thus tend to revert back the changes done to the elastic moduli,
  !       reusing the same arrays.

  ! beware the order here: in the preparation of the time loop, we first shift moduli, then unrelax them.
  !                        let's do now the opposite...

  ! restores relaxed elastic moduli after setting unrelaxed values before time loop
  call restore_relaxed_moduli()

  ! restore un-shifted elastic moduli (which correspond to the original moduli read in from mesher)
  call restore_unshifted_reference_moduli()

  ! debug
  !if (myrank == 0) print *,'debug: restored moduli muv',muvstore_crust_mantle(1,1,1,1000),muvstore_crust_mantle(3,3,3,3000)
  !if (myrank == 0) print *,'debug: restored moduli c11',c11store_crust_mantle(1,1,1,1000),c11store_crust_mantle(3,3,3,3000)
  !if (myrank == 0) print *,'debug: restored moduli c44',c44store_crust_mantle(1,1,1,1000),c44store_crust_mantle(3,3,3,3000)

  end subroutine restore_original_moduli

!
!-------------------------------------------------------------------------------------------------
!

  subroutine restore_relaxed_moduli()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: one_minus_sum_beta_use,minus_sum_beta,muvl,muhl

  ! for rotations
  double precision :: g11,g12,g13,g14,g15,g16,g22,g23,g24,g25,g26, &
                      g33,g34,g35,g36,g44,g45,g46,g55,g56,g66
  double precision :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66
  double precision :: phi_dble,theta_dble
  double precision :: A_dble,F_dble,L_dble,N_dble,eta_aniso !,C_dble

  integer :: ispec,iglob
#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
#else
  integer :: i,j,k
#endif

  ! check if anything to do
  if (.not. ATTENUATION_VAL) return

  ! only muvstore is used further and needs to be restored

  ! crust/mantle
  do ispec = 1,NSPEC_CRUST_MANTLE
    ! isotropic and tiso elements
    DO_LOOP_IJK
      ! precompute terms for attenuation if needed
      if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
        one_minus_sum_beta_use = one_minus_sum_beta_crust_mantle(INDEX_IJK,ispec)
      else
        one_minus_sum_beta_use = one_minus_sum_beta_crust_mantle(1,1,1,ispec)
      endif

      ! zero-shift case
      if (abs(one_minus_sum_beta_use) < TINYVAL) cycle

      if (ANISOTROPIC_3D_MANTLE_VAL) then
        ! aniso
        minus_sum_beta =  one_minus_sum_beta_use - 1.0_CUSTOM_REAL

        ! zero-shift
        if (abs(one_minus_sum_beta_use) < TINYVAL) cycle

        ! original routine
        !
        ! shifting in prepare_elastic_elements:
        ! c44 = c44store_crust_mantle(INDEX_IJK,ispec)
        ! mul = c44 * minus_sum_beta
        ! c44 = c44 + mul
        ! restores:
        !
        ! mul = c44store_crust_mantle(INDEX_IJK,ispec) / one_minus_sum_beta_use
        ! c44store_crust_mantle(INDEX_IJK,ispec) = mul

        ! shifting in prepare_elastic_elements:
        ! mul = mul * minus_sum_beta
        ! c11 = c11 + FOUR_THIRDS * mul ! * minus_sum_beta * mul
        ! c12 = c12 - TWO_THIRDS * mul
        ! c13 = c13 - TWO_THIRDS * mul
        ! c22 = c22 + FOUR_THIRDS * mul
        ! c23 = c23 - TWO_THIRDS * mul
        ! c33 = c33 + FOUR_THIRDS * mul
        ! c55 = c55 + mul
        ! c66 = c66 + mul
        !
        ! restores:
        ! mul = mul * minus_sum_beta
        ! c11store_crust_mantle(INDEX_IJK,ispec) = c11store_crust_mantle(INDEX_IJK,ispec) - FOUR_THIRDS * mul
        ! c12store_crust_mantle(INDEX_IJK,ispec) = c12store_crust_mantle(INDEX_IJK,ispec) + TWO_THIRDS * mul
        ! c13store_crust_mantle(INDEX_IJK,ispec) = c13store_crust_mantle(INDEX_IJK,ispec) + TWO_THIRDS * mul
        ! c22store_crust_mantle(INDEX_IJK,ispec) = c22store_crust_mantle(INDEX_IJK,ispec) - FOUR_THIRDS * mul
        ! c23store_crust_mantle(INDEX_IJK,ispec) = c23store_crust_mantle(INDEX_IJK,ispec) + TWO_THIRDS * mul
        ! c33store_crust_mantle(INDEX_IJK,ispec) = c33store_crust_mantle(INDEX_IJK,ispec) - FOUR_THIRDS * mul
        ! c55store_crust_mantle(INDEX_IJK,ispec) = c55store_crust_mantle(INDEX_IJK,ispec) - mul
        ! c66store_crust_mantle(INDEX_IJK,ispec) = c66store_crust_mantle(INDEX_IJK,ispec) - mul

        ! new shift version shifts moduli by separating muv and muh factors.
        ! unrelaxed moduli shift: see prepare_elastic_elements to shift back accordingly.

        ! local position (d_ij given in radial direction)
        ! only in case needed for rotation
        iglob = ibool_crust_mantle(INDEX_IJK,ispec)
        theta_dble = rstore_crust_mantle(2,iglob)
        phi_dble = rstore_crust_mantle(3,iglob)
        call reduce(theta_dble,phi_dble)

        g11 = c11store_crust_mantle(INDEX_IJK,ispec)
        g12 = c12store_crust_mantle(INDEX_IJK,ispec)
        g13 = c13store_crust_mantle(INDEX_IJK,ispec)
        g14 = c14store_crust_mantle(INDEX_IJK,ispec)
        g15 = c15store_crust_mantle(INDEX_IJK,ispec)
        g16 = c16store_crust_mantle(INDEX_IJK,ispec)
        g22 = c22store_crust_mantle(INDEX_IJK,ispec)
        g23 = c23store_crust_mantle(INDEX_IJK,ispec)
        g24 = c24store_crust_mantle(INDEX_IJK,ispec)
        g25 = c25store_crust_mantle(INDEX_IJK,ispec)
        g26 = c26store_crust_mantle(INDEX_IJK,ispec)
        g33 = c33store_crust_mantle(INDEX_IJK,ispec)
        g34 = c34store_crust_mantle(INDEX_IJK,ispec)
        g35 = c35store_crust_mantle(INDEX_IJK,ispec)
        g36 = c36store_crust_mantle(INDEX_IJK,ispec)
        g44 = c44store_crust_mantle(INDEX_IJK,ispec)
        g45 = c45store_crust_mantle(INDEX_IJK,ispec)
        g46 = c46store_crust_mantle(INDEX_IJK,ispec)
        g55 = c55store_crust_mantle(INDEX_IJK,ispec)
        g56 = c56store_crust_mantle(INDEX_IJK,ispec)
        g66 = c66store_crust_mantle(INDEX_IJK,ispec)

        call rotate_tensor_global_to_radial(theta_dble,phi_dble, &
                                            d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                            d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                                            g11,g12,g13,g14,g15,g16,g22,g23,g24,g25,g26, &
                                            g33,g34,g35,g36,g44,g45,g46,g55,g56,g66)

        ! new: tries to shift moduli by separating muv and muh factors.
        !      still needs rotations to rotate back and forth from SPECFEM global axis to a radial symmetry axis
        !      since this shift assumes a radial symmetry

        ! orig:
        !         L = 1/2 (d44 + d55)
        !         muvl' = L * minus_sum_beta
        !
        !         d44' = d44 + 1/2 (d44+d55)*minus_sum_beta = d44 + muvl'
        !         d55' = d55 + 1/2 (d44+d55)*minus_sum_beta
        !
        !         1/2 (d44' + d55') = 1/2 (d44+d55 + (d44+d55)*minus_sum_beta )
        !                           = 1/2 (d44+d55) * (1 + minus_sum_beta)
        !         -> L = 1/2 (d44+d55) = 1/2 (d44'+d55')/(1 + minus_sum_beta)
        !              = muv
        L_dble = 0.5d0 * (d44 + d55) / one_minus_sum_beta_use

        ! orig:
        !        N = 0.125d0 * (d11 + d22 - 2.d0 * d12 + 4.d0 * d66)
        !        muhl' = N * minus_sum_beta
        !
        !        d11' = d11 + 4/3 muhl'
        !        d12' = d12 - 2/3 muhl'
        !        d22' = d22 + 4/3 muhl'
        !        d66' = d66 + muhl'
        !
        !        N' = muh' = N * (1 + minus_sum_beta)
        !        -> N = N'/(1 + minus_sum_beta)
        !             = muh
        N_dble = 0.125d0 * (d11 + d22 - 2.d0 * d12 + 4.d0 * d66) / one_minus_sum_beta_use

        ! shifting back
        muvl = real( L_dble * minus_sum_beta ,kind=CUSTOM_REAL)
        muhl = real( N_dble * minus_sum_beta ,kind=CUSTOM_REAL)

        d11 = d11 - FOUR_THIRDS * muhl ! * minus_sum_beta * mul
        d12 = d12 + TWO_THIRDS * muhl
        d22 = d22 - FOUR_THIRDS * muhl
        d33 = d33 - FOUR_THIRDS * muvl
        d44 = d44 - muvl
        d55 = d55 - muvl
        d66 = d66 - muhl

        ! orig:
        !       A = 1/8 (3 d11 + 3 d22 + 2 d12 + 4 d66)
        !
        !       d11' = d11 - 4/3 muhl
        !       d22' = d22 - 4/3 muhl
        !       d12' = d12 + 2/3 muhl
        !       d66' = d66 - muhl
        !
        A_dble = 0.125d0 * (3.d0 * d11 + 3.d0 * d22 + 2.d0 * d12 + 4.d0 * d66)

        !unused: C_dble = d33

        ! orig:
        !       F = 1/2 (d13 + d23)
        !       eta_aniso = F / (A - 2 L)
        !
        !       d13' = d13 + eta_aniso * (4/3 muhl' - 2 muvl')
        !            = d13 + eta_aniso * (4/3 muhl*minus_sum_beta - 2 muvl *minus_sum_beta)
        !
        !       d23' = d23 + eta_aniso * (4/3 muhl' - 2 muvl')
        !
        !       F' = 1/2 (d13' + d23')
        !          = 1/2 (d13 + eta_aniso * (4/3 muhl' - 2 muvl') + d23 + eta_aniso * (4/3 muhl' - 2 muvl')
        !          = 1/2 (d13 + d23) + eta_aniso (4/3 muhl' - 2muvl')
        !          = F + F/(A-2L) * (4/3muhl' - 2muvl')
        !          = F (1 +(4/3muhl' - 2muvl')/(A-2L))
        !       -> F = F' / (1 + (4/3muhl' - 2muvl') / (A-2L))
        F_dble = 0.5d0 * (d13 + d23) / (1.d0 + (FOUR_THIRDS * muhl - 2.d0 * muvl) / (A_dble - 2.d0 * L_dble))
        eta_aniso = F_dble / (A_dble - 2.d0*L_dble)   ! eta = F / (A-2L)

        d13 = d13 - eta_aniso * (FOUR_THIRDS * muhl - 2.d0*muvl)
        d23 = d23 - eta_aniso * (FOUR_THIRDS * muhl - 2.d0*muvl)

        ! debug
        !if (myrank == 0 .and. ispec == 1000 .and. ijk == 1) &
        !  print *,'debug: original moduli unrelaxing A,N,L,F,eta',A_dble,N_dble,L_dble,F_dble,eta_aniso,'mu',muvl,muhl

        ! rotates to global reference system
        call rotate_tensor_radial_to_global(theta_dble,phi_dble, &
                                            d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                            d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                                            g11,g12,g13,g14,g15,g16,g22,g23,g24,g25,g26, &
                                            g33,g34,g35,g36,g44,g45,g46,g55,g56,g66)

        ! stores unrelaxed factors
        c11store_crust_mantle(INDEX_IJK,ispec) = real(g11,kind=CUSTOM_REAL)
        c12store_crust_mantle(INDEX_IJK,ispec) = real(g12,kind=CUSTOM_REAL)
        c13store_crust_mantle(INDEX_IJK,ispec) = real(g13,kind=CUSTOM_REAL)
        c14store_crust_mantle(INDEX_IJK,ispec) = real(g14,kind=CUSTOM_REAL)
        c15store_crust_mantle(INDEX_IJK,ispec) = real(g15,kind=CUSTOM_REAL)
        c16store_crust_mantle(INDEX_IJK,ispec) = real(g16,kind=CUSTOM_REAL)
        c22store_crust_mantle(INDEX_IJK,ispec) = real(g22,kind=CUSTOM_REAL)
        c23store_crust_mantle(INDEX_IJK,ispec) = real(g23,kind=CUSTOM_REAL)
        c24store_crust_mantle(INDEX_IJK,ispec) = real(g24,kind=CUSTOM_REAL)
        c25store_crust_mantle(INDEX_IJK,ispec) = real(g25,kind=CUSTOM_REAL)
        c26store_crust_mantle(INDEX_IJK,ispec) = real(g26,kind=CUSTOM_REAL)
        c33store_crust_mantle(INDEX_IJK,ispec) = real(g33,kind=CUSTOM_REAL)
        c34store_crust_mantle(INDEX_IJK,ispec) = real(g34,kind=CUSTOM_REAL)
        c35store_crust_mantle(INDEX_IJK,ispec) = real(g35,kind=CUSTOM_REAL)
        c36store_crust_mantle(INDEX_IJK,ispec) = real(g36,kind=CUSTOM_REAL)
        c44store_crust_mantle(INDEX_IJK,ispec) = real(g44,kind=CUSTOM_REAL)
        c45store_crust_mantle(INDEX_IJK,ispec) = real(g45,kind=CUSTOM_REAL)
        c46store_crust_mantle(INDEX_IJK,ispec) = real(g46,kind=CUSTOM_REAL)
        c55store_crust_mantle(INDEX_IJK,ispec) = real(g55,kind=CUSTOM_REAL)
        c56store_crust_mantle(INDEX_IJK,ispec) = real(g56,kind=CUSTOM_REAL)
        c66store_crust_mantle(INDEX_IJK,ispec) = real(g66,kind=CUSTOM_REAL)

        muvstore_crust_mantle(INDEX_IJK,ispec) = real(L_dble,kind=CUSTOM_REAL)

      else
        ! layer with both iso and transverse isotropy elements, use kappav and muv
        muvl = muvstore_crust_mantle(INDEX_IJK,ispec)
        ! returns to the relaxed moduli Mu_r = Mu_u / [1 - sum(1 - tau_strain/tau_stress) ]
        muvl = muvl / one_minus_sum_beta_use
        ! stores relaxed shear moduli for kernel computations
        muvstore_crust_mantle(INDEX_IJK,ispec) = muvl

        ! tiso elements also use muh
        if (ispec_is_tiso_crust_mantle(ispec)) then
          muhl = muhstore_crust_mantle(INDEX_IJK,ispec)
          muhl = muhl / one_minus_sum_beta_use
          muhstore_crust_mantle(INDEX_IJK,ispec) = muhl
        endif
      endif
    ENDDO_LOOP_IJK
  enddo
  ! see prepare_elastic_elements:
  ! since we scale muv and c11,.. stores we must divide with this factor to use the relaxed moduli for the modulus defect
  ! calculation in updating the memory variables
  !
  ! note: we won't need to restores factor_common here.
  !       it is only needed for time-stepping of the attenuation memory variables.

  ! inner core
  if (SAVE_KERNELS_IC) then
    do ispec = 1,NSPEC_INNER_CORE
      ! exclude fictitious elements in central cube
      if (idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle
      ! isotropic element
      DO_LOOP_IJK
        ! precompute terms for attenuation if needed
        if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
          one_minus_sum_beta_use = one_minus_sum_beta_inner_core(INDEX_IJK,ispec)
        else
          one_minus_sum_beta_use = one_minus_sum_beta_inner_core(1,1,1,ispec)
        endif

        ! layer with no transverse isotropy, use kappav and muv
        muvl = muvstore_inner_core(INDEX_IJK,ispec)
        ! returns to the relaxed moduli Mu_r = Mu_u / [1 - sum(1 - tau_strain/tau_stress) ]
        muvl = muvl / one_minus_sum_beta_use
        ! stores relaxed shear moduli for kernel computations
        muvstore_inner_core(INDEX_IJK,ispec) = muvl
      ENDDO_LOOP_IJK
    enddo
    ! note: we won't need to restores factor_common here.
    !       it is only needed for time-stepping of the attenuation memory variables.
  endif ! SAVE_KERNELS_IC

  end subroutine restore_relaxed_moduli

!
!-------------------------------------------------------------------------------------------------
!

  subroutine restore_unshifted_reference_moduli()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: scale_factor,scale_factor_minus_one,mul,muvl,muhl

  ! for rotations
  double precision :: g11,g12,g13,g14,g15,g16,g22,g23,g24,g25,g26, &
                      g33,g34,g35,g36,g44,g45,g46,g55,g56,g66
  double precision :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66
  double precision :: phi_dble,theta_dble
  double precision :: A_dble,F_dble,L_dble,N_dble,eta_aniso !,C_dble

  integer :: ispec,iglob
#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
#else
  integer :: i,j,k
#endif

  ! check if anything to do
  if (.not. ATTENUATION_VAL) return

  ! crust/mantle
  do ispec = 1,NSPEC_CRUST_MANTLE
    ! isotropic and tiso elements
    DO_LOOP_IJK
      if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
        scale_factor = factor_scale_crust_mantle(INDEX_IJK,ispec)
      else
        scale_factor = factor_scale_crust_mantle(1,1,1,ispec)
      endif

      ! zero-scale factor
      if (abs(scale_factor) < TINYVAL) cycle

      if (ANISOTROPIC_3D_MANTLE_VAL) then
        ! anisotropic element
        scale_factor_minus_one = scale_factor - 1.0_CUSTOM_REAL

        ! shifting: (in prepare_attenuation.f90)

        ! old routine:
        ! mul = c44store_crust_mantle(i,j,k,ispec)
        ! c44store_crust_mantle(i,j,k,ispec) = c44store_crust_mantle(i,j,k,ispec) + scale_factor_minus_one * mul
        ! restores like:
        ! mul = c44store_crust_mantle(INDEX_IJK,ispec)/scale_factor ! equals mu = c44store_crust_mantle/(1+scale_factor_minus_one)
        ! c44store_crust_mantle(INDEX_IJK,ispec) = mul
        !
        ! see shifting in prepare_attenuation.f90:
        ! c11store_crust_mantle(i,j,k,ispec) = c11store_crust_mantle(i,j,k,ispec) + FOUR_THIRDS * scale_factor_minus_one * mul
        ! c12store_crust_mantle(i,j,k,ispec) = c12store_crust_mantle(i,j,k,ispec) - TWO_THIRDS * scale_factor_minus_one * mul
        ! c13store_crust_mantle(i,j,k,ispec) = c13store_crust_mantle(i,j,k,ispec) - TWO_THIRDS * scale_factor_minus_one * mul
        ! c22store_crust_mantle(i,j,k,ispec) = c22store_crust_mantle(i,j,k,ispec) + FOUR_THIRDS * scale_factor_minus_one * mul
        ! c23store_crust_mantle(i,j,k,ispec) = c23store_crust_mantle(i,j,k,ispec) - TWO_THIRDS * scale_factor_minus_one * mul
        ! c33store_crust_mantle(i,j,k,ispec) = c33store_crust_mantle(i,j,k,ispec) + FOUR_THIRDS * scale_factor_minus_one * mul
        ! c44store_crust_mantle(i,j,k,ispec) = c44store_crust_mantle(i,j,k,ispec) + scale_factor_minus_one * mul
        ! c55store_crust_mantle(i,j,k,ispec) = c55store_crust_mantle(i,j,k,ispec) + scale_factor_minus_one * mul
        ! c66store_crust_mantle(i,j,k,ispec) = c66store_crust_mantle(i,j,k,ispec) + scale_factor_minus_one * mul
        !
        ! becomes (sign change of last term)
        !c11store_crust_mantle(INDEX_IJK,ispec) = c11store_crust_mantle(INDEX_IJK,ispec) &
        !                                         - FOUR_THIRDS * scale_factor_minus_one * mul
        !c12store_crust_mantle(INDEX_IJK,ispec) = c12store_crust_mantle(INDEX_IJK,ispec) &
        !                                         + TWO_THIRDS * scale_factor_minus_one * mul
        !c13store_crust_mantle(INDEX_IJK,ispec) = c13store_crust_mantle(INDEX_IJK,ispec) &
        !                                         + TWO_THIRDS * scale_factor_minus_one * mul
        !c22store_crust_mantle(INDEX_IJK,ispec) = c22store_crust_mantle(INDEX_IJK,ispec) &
        !                                         - FOUR_THIRDS * scale_factor_minus_one * mul
        !c23store_crust_mantle(INDEX_IJK,ispec) = c23store_crust_mantle(INDEX_IJK,ispec) &
        !                                         + TWO_THIRDS * scale_factor_minus_one * mul
        !c33store_crust_mantle(INDEX_IJK,ispec) = c33store_crust_mantle(INDEX_IJK,ispec) &
        !                                         - FOUR_THIRDS * scale_factor_minus_one * mul
        !c55store_crust_mantle(INDEX_IJK,ispec) = c55store_crust_mantle(INDEX_IJK,ispec) &
        !                                         - scale_factor_minus_one * mul
        !c66store_crust_mantle(INDEX_IJK,ispec) = c66store_crust_mantle(INDEX_IJK,ispec) &
        !                                         - scale_factor_minus_one * mul


        ! new shift version shifts moduli by separating muv and muh factors.
        ! local position (d_ij given in radial direction)
        ! only in case needed for rotation
        iglob = ibool_crust_mantle(INDEX_IJK,ispec)
        theta_dble = rstore_crust_mantle(2,iglob)
        phi_dble = rstore_crust_mantle(3,iglob)
        call reduce(theta_dble,phi_dble)

        g11 = c11store_crust_mantle(INDEX_IJK,ispec)
        g12 = c12store_crust_mantle(INDEX_IJK,ispec)
        g13 = c13store_crust_mantle(INDEX_IJK,ispec)
        g14 = c14store_crust_mantle(INDEX_IJK,ispec)
        g15 = c15store_crust_mantle(INDEX_IJK,ispec)
        g16 = c16store_crust_mantle(INDEX_IJK,ispec)
        g22 = c22store_crust_mantle(INDEX_IJK,ispec)
        g23 = c23store_crust_mantle(INDEX_IJK,ispec)
        g24 = c24store_crust_mantle(INDEX_IJK,ispec)
        g25 = c25store_crust_mantle(INDEX_IJK,ispec)
        g26 = c26store_crust_mantle(INDEX_IJK,ispec)
        g33 = c33store_crust_mantle(INDEX_IJK,ispec)
        g34 = c34store_crust_mantle(INDEX_IJK,ispec)
        g35 = c35store_crust_mantle(INDEX_IJK,ispec)
        g36 = c36store_crust_mantle(INDEX_IJK,ispec)
        g44 = c44store_crust_mantle(INDEX_IJK,ispec)
        g45 = c45store_crust_mantle(INDEX_IJK,ispec)
        g46 = c46store_crust_mantle(INDEX_IJK,ispec)
        g55 = c55store_crust_mantle(INDEX_IJK,ispec)
        g56 = c56store_crust_mantle(INDEX_IJK,ispec)
        g66 = c66store_crust_mantle(INDEX_IJK,ispec)

        call rotate_tensor_global_to_radial(theta_dble,phi_dble, &
                                            d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                            d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                                            g11,g12,g13,g14,g15,g16,g22,g23,g24,g25,g26, &
                                            g33,g34,g35,g36,g44,g45,g46,g55,g56,g66)

        ! forward shifts like: (see prepare_attenuation)
        !  A_dble = 0.125d0 * (3.d0 * d11 + 3.d0 * d22 + 2.d0 * d12 + 4.d0 * d66)
        !  N_dble = 0.125d0 * (d11 + d22 - 2.d0 * d12 + 4.d0 * d66)
        !  L_dble = 0.5d0 * (d44 + d55)
        !  F_dble = 0.5d0 * (d13 + d23)
        !  eta_aniso = F_dble / (A_dble - 2.d0*L_dble)   ! eta = F / (A-2L)
        !  muvl = L_dble * scale_factor_minus_one     ! c44 -> L -> muv
        !  muhl = N_dble * scale_factor_minus_one     ! c66 -> N -> muh
        !  d11 = d11 + FOUR_THIRDS * muhl ! * minus_sum_beta * mul
        !  d12 = d12 - TWO_THIRDS * muhl
        !  d13 = d13 + eta_aniso * (FOUR_THIRDS * muhl - 2.d0*muvl)
        !  d22 = d22 + FOUR_THIRDS * muhl
        !  d23 = d23 + eta_aniso * (FOUR_THIRDS * muhl - 2.d0*muvl)
        !  d33 = d33 + FOUR_THIRDS * muvl
        !  d44 = d44 + muvl
        !  d55 = d55 + muvl
        !  d66 = d66 + muhl
        !
        ! scaling similar to unrelaxed one, but with different scaling factor
        !
        ! orig:
        !         L = 1/2 (d44 + d55)
        !         muvl' = L * scale_factor_minus_one
        !
        !         d44' = d44 + 1/2 (d44+d55)*scale_factor_minus_one = d44 + muvl'
        !         d55' = d55 + 1/2 (d44+d55)*scale_factor_minus_one
        !
        !         1/2 (d44' + d55') = 1/2 (d44+d55 + (d44+d55)*scale_factor_minus_one )
        !                           = 1/2 (d44+d55) * (1 + scale_factor_minus_one)
        !         -> L = 1/2 (d44+d55) = 1/2 (d44'+d55')/(1 + scale_factor_minus_one)
        !              = muv
        L_dble = 0.5d0 * (d44 + d55) / scale_factor

        ! orig:
        !        N = 0.125d0 * (d11 + d22 - 2.d0 * d12 + 4.d0 * d66)
        !        muhl' = N * scale_factor_minus_one
        !
        !        d11' = d11 + 4/3 muhl'
        !        d12' = d12 - 2/3 muhl'
        !        d22' = d22 + 4/3 muhl'
        !        d66' = d66 + muhl'
        !
        !        N' = muh' = N * (1 + minus_sum_beta)
        !        -> N = N'/(1 + minus_sum_beta)
        !             = muh
        N_dble = 0.125d0 * (d11 + d22 - 2.d0 * d12 + 4.d0 * d66) / scale_factor

        ! shifting back
        muvl = real( L_dble * scale_factor_minus_one ,kind=CUSTOM_REAL)
        muhl = real( N_dble * scale_factor_minus_one ,kind=CUSTOM_REAL)

        d11 = d11 - FOUR_THIRDS * muhl ! * minus_sum_beta * mul
        d12 = d12 + TWO_THIRDS * muhl
        d22 = d22 - FOUR_THIRDS * muhl
        d33 = d33 - FOUR_THIRDS * muvl
        d44 = d44 - muvl
        d55 = d55 - muvl
        d66 = d66 - muhl

        ! orig:
        !       A = 1/8 (3 d11 + 3 d22 + 2 d12 + 4 d66)
        !
        !       d11' = d11 - 4/3 muhl
        !       d22' = d22 - 4/3 muhl
        !       d12' = d12 + 2/3 muhl
        !       d66' = d66 - muhl
        !
        A_dble = 0.125d0 * (3.d0 * d11 + 3.d0 * d22 + 2.d0 * d12 + 4.d0 * d66)

        !unused: C_dble = d33

        ! orig:
        !       F = 1/2 (d13 + d23)
        !       eta_aniso = F / (A - 2 L)
        !
        !       d13' = d13 + eta_aniso * (4/3 muhl' - 2 muvl')
        !            = d13 + eta_aniso * (4/3 muhl*scale_factor_minus_one) - 2 muvl *scale_factor_minus_one)
        !
        !       d23' = d23 + eta_aniso * (4/3 muhl' - 2 muvl')
        !
        !       F' = 1/2 (d13' + d23')
        !          = 1/2 (d13 + eta_aniso * (4/3 muhl' - 2 muvl') + d23 + eta_aniso * (4/3 muhl' - 2 muvl')
        !          = 1/2 (d13 + d23) + eta_aniso (4/3 muhl' - 2muvl')
        !          = F + F/(A-2L) (4/3muhl' - 2muvl')
        !          = F (1 +(4/3muhl' - 2muvl')/(A-2L))
        !       -> F = F' / (1 + (4/3muhl' - 2muvl') / (A-2L))
        F_dble = 0.5d0 * (d13 + d23) / (1.d0 + (FOUR_THIRDS * muhl - 2.d0 * muvl) / (A_dble - 2.d0 * L_dble))
        eta_aniso = F_dble / (A_dble - 2.d0*L_dble)   ! eta = F / (A-2L)

        d13 = d13 - eta_aniso * (FOUR_THIRDS * muhl - 2.d0*muvl)
        d23 = d23 - eta_aniso * (FOUR_THIRDS * muhl - 2.d0*muvl)

        ! debug
        !if (myrank == 0 .and. ispec == 1000 .and. ijk == 1) &
        !  print *,'debug: original moduli unscaling A,N,L,F,eta',A_dble,N_dble,L_dble,F_dble,eta_aniso,'mu',muvl,muhl

        ! rotates to global reference system
        call rotate_tensor_radial_to_global(theta_dble,phi_dble, &
                                            d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                            d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                                            g11,g12,g13,g14,g15,g16,g22,g23,g24,g25,g26, &
                                            g33,g34,g35,g36,g44,g45,g46,g55,g56,g66)

        ! stores unrelaxed factors
        c11store_crust_mantle(INDEX_IJK,ispec) = real(g11,kind=CUSTOM_REAL)
        c12store_crust_mantle(INDEX_IJK,ispec) = real(g12,kind=CUSTOM_REAL)
        c13store_crust_mantle(INDEX_IJK,ispec) = real(g13,kind=CUSTOM_REAL)
        c14store_crust_mantle(INDEX_IJK,ispec) = real(g14,kind=CUSTOM_REAL)
        c15store_crust_mantle(INDEX_IJK,ispec) = real(g15,kind=CUSTOM_REAL)
        c16store_crust_mantle(INDEX_IJK,ispec) = real(g16,kind=CUSTOM_REAL)
        c22store_crust_mantle(INDEX_IJK,ispec) = real(g22,kind=CUSTOM_REAL)
        c23store_crust_mantle(INDEX_IJK,ispec) = real(g23,kind=CUSTOM_REAL)
        c24store_crust_mantle(INDEX_IJK,ispec) = real(g24,kind=CUSTOM_REAL)
        c25store_crust_mantle(INDEX_IJK,ispec) = real(g25,kind=CUSTOM_REAL)
        c26store_crust_mantle(INDEX_IJK,ispec) = real(g26,kind=CUSTOM_REAL)
        c33store_crust_mantle(INDEX_IJK,ispec) = real(g33,kind=CUSTOM_REAL)
        c34store_crust_mantle(INDEX_IJK,ispec) = real(g34,kind=CUSTOM_REAL)
        c35store_crust_mantle(INDEX_IJK,ispec) = real(g35,kind=CUSTOM_REAL)
        c36store_crust_mantle(INDEX_IJK,ispec) = real(g36,kind=CUSTOM_REAL)
        c44store_crust_mantle(INDEX_IJK,ispec) = real(g44,kind=CUSTOM_REAL)
        c45store_crust_mantle(INDEX_IJK,ispec) = real(g45,kind=CUSTOM_REAL)
        c46store_crust_mantle(INDEX_IJK,ispec) = real(g46,kind=CUSTOM_REAL)
        c55store_crust_mantle(INDEX_IJK,ispec) = real(g55,kind=CUSTOM_REAL)
        c56store_crust_mantle(INDEX_IJK,ispec) = real(g56,kind=CUSTOM_REAL)
        c66store_crust_mantle(INDEX_IJK,ispec) = real(g66,kind=CUSTOM_REAL)

        ! for solving memory-variables, modulus defect \delta \mu_l (Komatitsch, 2002, eq. (11) & (13))
        ! note: for solving the memory variables, we will only use the modulus defect
        !       associated with muv. this is consistent with the implementation for tiso below.
        !
        !       however, to properly account for shear attenuation, one might have to add also
        !       memory-variables for a modulus defect associated with muh.
        muvstore_crust_mantle(INDEX_IJK,ispec) = real(L_dble,kind=CUSTOM_REAL)

      else
        ! isotropic or transverse isotropic element
        muvstore_crust_mantle(INDEX_IJK,ispec) = muvstore_crust_mantle(INDEX_IJK,ispec) / scale_factor
        ! scales transverse isotropic values for mu_h
        if (ispec_is_tiso_crust_mantle(ispec)) then
          muhstore_crust_mantle(INDEX_IJK,ispec) = muhstore_crust_mantle(INDEX_IJK,ispec) / scale_factor
        endif
      endif
    ENDDO_LOOP_IJK
  enddo

  ! inner core
  if (SAVE_KERNELS_IC) then
    do ispec = 1,NSPEC_INNER_CORE
      ! exclude fictitious elements in central cube
      if (idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

      DO_LOOP_IJK
        if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
          scale_factor = factor_scale_inner_core(INDEX_IJK,ispec)
        else
          scale_factor = factor_scale_inner_core(1,1,1,ispec)
        endif

        ! zero-scale factor
        if (abs(scale_factor) < TINYVAL) cycle

        if (ANISOTROPIC_INNER_CORE_VAL) then
          scale_factor_minus_one = scale_factor - 1.0_CUSTOM_REAL

          mul = c44store_inner_core(INDEX_IJK,ispec) / scale_factor
          c44store_inner_core(INDEX_IJK,ispec) = mul

          c11store_inner_core(INDEX_IJK,ispec) = c11store_inner_core(INDEX_IJK,ispec) &
                                                 - FOUR_THIRDS * scale_factor_minus_one * mul
          c12store_inner_core(INDEX_IJK,ispec) = c12store_inner_core(INDEX_IJK,ispec) &
                                                 + TWO_THIRDS * scale_factor_minus_one * mul
          c13store_inner_core(INDEX_IJK,ispec) = c13store_inner_core(INDEX_IJK,ispec) &
                                                 + TWO_THIRDS * scale_factor_minus_one * mul
          c33store_inner_core(INDEX_IJK,ispec) = c33store_inner_core(INDEX_IJK,ispec) &
                                                 - FOUR_THIRDS * scale_factor_minus_one * mul
        endif

        muvstore_inner_core(INDEX_IJK,ispec) = muvstore_inner_core(INDEX_IJK,ispec) / scale_factor

      ENDDO_LOOP_IJK
    enddo
  endif

  end subroutine restore_unshifted_reference_moduli


!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_kernels_crust_mantle()

  use specfem_par, only: SAVE_REGULAR_KL,ANISOTROPIC_KL,FULL_GRAVITY_VAL

  implicit none

  ! outputs sensitivity kernels (in SEM-format) to file
  if (ANISOTROPIC_KL) then
    ! anisotropic kernels
    call save_kernels_crust_mantle_ani()
  else
    ! isotropic kernels
    call save_kernels_crust_mantle_iso()
  endif

  ! stores additional kernels on a regular grid
  if (SAVE_REGULAR_KL) call save_regular_kernels_cm()

  ! additional full gravity kernels
  if (FULL_GRAVITY_VAL) call SIEM_save_crust_mantle_kernels()

  end subroutine save_kernels_crust_mantle

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_kernels_crust_mantle_ani()

! stores kernels for anisotropic/transversely isotropic parameterizations

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(21) ::  cijkl_kl_local
  real(kind=CUSTOM_REAL) :: scale_kl,scale_kl_ani,scale_kl_rho,scaleval,scale_GPa
  real(kind=CUSTOM_REAL) :: rhol,mul,kappal
  real(kind=CUSTOM_REAL) :: theta,phi

  integer :: ispec,i,j,k,iglob
  integer :: ier

  ! primary rho kernel
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: &
    rhonotprime_kl_crust_mantle

  ! transverse isotropic parameters
  real(kind=CUSTOM_REAL), dimension(21) :: an_kl
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: &
    alphav_kl_crust_mantle,alphah_kl_crust_mantle, &
    betav_kl_crust_mantle,betah_kl_crust_mantle, &
    eta_kl_crust_mantle

  ! bulk parameterization
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: &
    bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle, &
    bulk_betav_kl_crust_mantle,bulk_betah_kl_crust_mantle
  real(kind=CUSTOM_REAL) :: A,C,F,L,N,eta,mu0
  real(kind=CUSTOM_REAL) :: muvl,kappavl,muhl,kappahl
  real(kind=CUSTOM_REAL) :: alphav_sq,alphah_sq,betav_sq,betah_sq,bulk_sq

  ! azimuthal
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: &
    Gc_kl_crust_mantle, Gs_kl_crust_mantle, &
    Gc_prime_kl_crust_mantle, Gs_prime_kl_crust_mantle, &
    A_kl_crust_mantle, C_kl_crust_mantle, F_kl_crust_mantle, &
    N_kl_crust_mantle, L_kl_crust_mantle, &
    Jc_kl_crust_mantle, Kc_kl_crust_mantle, Mc_kl_crust_mantle, &
    Bc_kl_crust_mantle, Hc_kl_crust_mantle, Ec_kl_crust_mantle, &
    Dc_kl_crust_mantle
  real(kind=CUSTOM_REAL), dimension(21) :: cij,cij_radial

#ifdef USE_CEM
  character(len=MAX_STRING_LEN) :: filename
#endif

  ! checks if anything to do
  if (.not. ANISOTROPIC_KL) return

  ! scaling factors: note that this scaling has been introduced by Qinya Liu (2006)
  !                  with the intent to dimensionalize kernel values to [ s km^(-3) ]
  !
  ! kernel unit [ s / km^3 ]
  scale_kl = real(scale_t * scale_displ_inv * 1.d9,kind=CUSTOM_REAL)
  ! For anisotropic kernels
  ! final unit : [s km^(-3) GPa^(-1)]
  scale_kl_ani = real(scale_t**3 / (RHOAV*R_PLANET**3) * 1.d18,kind=CUSTOM_REAL)
  ! final unit : [s km^(-3) (kg/m^3)^(-1)]
  scale_kl_rho = real(scale_t * scale_displ_inv / RHOAV * 1.d9,kind=CUSTOM_REAL)
  ! the scale of GPa--[g/cm^3][(km/s)^2]
  scaleval = real(sqrt(PI*GRAV*RHOAV),kind=CUSTOM_REAL)
  scale_GPa = real((RHOAV/1000.d0)*((R_PLANET*scaleval/1000.d0)**2),kind=CUSTOM_REAL)

  ! debug
  !if (myrank == 0) print *,'debug: save kernels: scaling factors',scale_kl,scale_kl_ani,scale_kl_rho

  ! allocates temporary arrays
  if (SAVE_TRANSVERSE_KL_ONLY) then
    ! transverse isotropic kernel arrays for file output
    allocate(alphav_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             alphah_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             betav_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             betah_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             eta_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating transverse kernels alphav_kl_crust_mantle,...'

    ! isotropic kernel arrays for file output
    allocate(bulk_betav_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             bulk_betah_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating transverse kernels bulk_betav_kl_crust_mantle,...'

    ! bulk velocity kernels
    allocate(bulk_c_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             bulk_beta_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating transverse kernels bulk_c_kl_crust_mantle,...'

    ! only dummy
    allocate(Gc_prime_kl_crust_mantle(1,1,1,1), &
             Gs_prime_kl_crust_mantle(1,1,1,1))

  else if (SAVE_AZIMUTHAL_ANISO_KL_ONLY) then
    ! azimuthal anisotropic kernels
    allocate(Gc_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             Gs_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             A_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             C_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             F_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             L_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             N_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             Jc_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             Kc_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             Mc_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             Bc_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             Hc_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             Ec_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             Dc_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating azimuthal kernels alphav_kl_crust_mantle,...'

    ! azimuthally anisotropic kernel arrays for file output
    allocate(alphav_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             alphah_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             betav_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             betah_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             eta_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             Gc_prime_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             Gs_prime_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating azimuthal kernels alphav_kl_crust_mantle,...'

    allocate(bulk_betav_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             bulk_betah_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
    if (ier /= 0 ) stop 'Error allocating transverse kernels bulk_betav_kl_crust_mantle,...'

    ! bulk velocity kernels
    allocate(bulk_c_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             bulk_beta_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
    if (ier /= 0 ) stop 'Error allocating transverse kernels bulk_c_kl_crust_mantle,...'

  else
    ! dummy allocation
    allocate(alphav_kl_crust_mantle(1,1,1,1), &
             alphah_kl_crust_mantle(1,1,1,1), &
             betav_kl_crust_mantle(1,1,1,1), &
             betah_kl_crust_mantle(1,1,1,1), &
             eta_kl_crust_mantle(1,1,1,1))
    ! isotropic kernel arrays for file output
    allocate(bulk_betav_kl_crust_mantle(1,1,1,1), &
             bulk_betah_kl_crust_mantle(1,1,1,1))
    ! only dummy
    allocate(bulk_c_kl_crust_mantle(1,1,1,1), &
             bulk_beta_kl_crust_mantle(1,1,1,1))
    ! only dummy
    allocate(Gc_prime_kl_crust_mantle(1,1,1,1), &
             Gs_prime_kl_crust_mantle(1,1,1,1))
  endif

  allocate(rhonotprime_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
  if (ier /= 0 ) stop 'Error allocating transverse kernels rhonotprime_kl_crust_mantle'

  ! crust_mantle
  do ispec = 1, NSPEC_CRUST_MANTLE_ADJOINT
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

          ! For anisotropic kernels
          iglob = ibool_crust_mantle(i,j,k,ispec)

          ! The Cartesian global cijkl_kl are rotated into the local (radial) cijkl_kl
          ! ystore and zstore are thetaval and phival (line 2252) -- dangerous
          theta = rstore_crust_mantle(2,iglob)
          phi = rstore_crust_mantle(3,iglob)

          call rotate_tensor_global_to_radial_vector(cijkl_kl_crust_mantle(:,i,j,k,ispec),cijkl_kl_local(:),theta,phi)

          cijkl_kl_crust_mantle(:,i,j,k,ispec) = cijkl_kl_local(:) * scale_kl_ani
          rho_kl_crust_mantle(i,j,k,ispec) = rho_kl_crust_mantle(i,j,k,ispec) * scale_kl_rho

          ! transverse isotropic kernel calculations
          if (SAVE_TRANSVERSE_KL_ONLY .or. SAVE_AZIMUTHAL_ANISO_KL_ONLY) then
            ! note: transverse isotropic kernels are calculated for all elements
            !
            !          however, the factors A,C,L,N,F are based only on transverse elements
            !          in between Moho and 220 km layer, otherwise they are taken from isotropic values

            rhol = rhostore_crust_mantle(i,j,k,ispec)

            ! transverse isotropic parameters from compute_force_crust_mantle.f90
            ! C=rhovpvsq A=rhovphsq L=rhovsvsq N=rhovshsq eta=F/(A - 2 L)

            if (.not. ANISOTROPIC_3D_MANTLE) then
              ! Get A,C,F,L,N,eta from kappa,mu
              ! element can have transverse isotropy if between d220 and Moho
              if (.not. ispec_is_tiso_crust_mantle(ispec)) then
                ! isotropic element
                ! layer with no transverse isotropy
                ! A,C,L,N,F from isotropic model
                mul = muvstore_crust_mantle(i,j,k,ispec)
                kappal = kappavstore_crust_mantle(i,j,k,ispec)
                muvl = mul
                muhl = mul

                A = kappal + FOUR_THIRDS * mul
                C = A
                L = mul
                N = mul
                F = kappal - 2._CUSTOM_REAL/3._CUSTOM_REAL * mul
                eta = 1._CUSTOM_REAL
              else
                ! tiso element
                ! A,C,L,N,F from transverse isotropic model
                kappavl = kappavstore_crust_mantle(i,j,k,ispec)
                kappahl = kappahstore_crust_mantle(i,j,k,ispec)
                muvl = muvstore_crust_mantle(i,j,k,ispec)
                muhl = muhstore_crust_mantle(i,j,k,ispec)
                kappal = kappavl

                A = kappahl + FOUR_THIRDS * muhl
                C = kappavl + FOUR_THIRDS * muvl
                L = muvl
                N = muhl
                eta = eta_anisostore_crust_mantle(i,j,k,ispec)  ! that is  F / (A - 2 L)
                F = eta * ( A - 2._CUSTOM_REAL * L )
              endif
            else
              ! anisotropic mantle
              !
              ! Sieminski, 2007: assuming vertical(radial) symmetry axis
              ! A = 1/8 (3 C11 + 3 C22 + 2 C12 + 4 C66)       -> kernel K_A = K_C11 + K_C12 + K_C22
              ! C = C33                                       -> kernel K_C = K_C33
              ! N = 1/8 (C11 + C22 - 2 C12 + 4 C66)           -> kernel K_N = K_C66 - 2 K_C12
              ! L = 1/2 (C44 + C55)                           -> kernel K_L = K_C44 + K_C55
              ! F = 1/2 (C13 + C23)                           -> kernel K_F = K_C13 + K_C23
              ! eta = F / (A - 2 L)
              !
              ! and (also, Anderson & Dziewonski):
              ! A = rho * vph**2
              ! C = rho * vpv**2
              ! N = rho * vsh**2 = muh
              ! L = rho * vsv**2 = muv
              ! F = eta * (A - 2*L)

              ! daniel todo:
              ! Love parameters derived from elastic tensor cij
              !
              ! please check, this assumes a radial symmetry axis, thus we should again first rotate back
              ! the c11store,.. coefficients from the SPECFEM reference coordinate system
              cij(1) = c11store_crust_mantle(i,j,k,ispec); cij(2) = c12store_crust_mantle(i,j,k,ispec)
              cij(3) = c13store_crust_mantle(i,j,k,ispec); cij(4) = c14store_crust_mantle(i,j,k,ispec)
              cij(5) = c15store_crust_mantle(i,j,k,ispec); cij(6) = c16store_crust_mantle(i,j,k,ispec)
              cij(7) = c22store_crust_mantle(i,j,k,ispec); cij(8) = c23store_crust_mantle(i,j,k,ispec)
              cij(9) = c24store_crust_mantle(i,j,k,ispec); cij(10) = c25store_crust_mantle(i,j,k,ispec)
              cij(11) = c26store_crust_mantle(i,j,k,ispec); cij(12) = c33store_crust_mantle(i,j,k,ispec)
              cij(13) = c34store_crust_mantle(i,j,k,ispec); cij(14) = c35store_crust_mantle(i,j,k,ispec)
              cij(15) = c36store_crust_mantle(i,j,k,ispec); cij(16) = c44store_crust_mantle(i,j,k,ispec)
              cij(17) = c45store_crust_mantle(i,j,k,ispec); cij(18) = c46store_crust_mantle(i,j,k,ispec)
              cij(19) = c55store_crust_mantle(i,j,k,ispec); cij(20) = c56store_crust_mantle(i,j,k,ispec)
              cij(21) = c66store_crust_mantle(i,j,k,ispec)

              call rotate_tensor_global_to_radial_vector(cij(:),cij_radial(:),theta,phi)

              ! A = 1/8 ( 3 c11 + 3 c22 + 2 c12 + 4 c66)
              A = 0.125_CUSTOM_REAL * ( 3.0 * cij_radial(1) + 3.0 * cij_radial(7) &
                                       + 2.0 * cij_radial(2) + 4.0 * cij_radial(21) )
              ! C = c33
              C = cij_radial(12)
              ! N = 1/8 ( c11 + c22 - 2 c12 + 4 c66)
              N = 0.125_CUSTOM_REAL * ( cij_radial(1) + cij_radial(7) &
                                       - 2.0 * cij_radial(2) + 4.0 * cij_radial(21))
              ! L = 1/2 (c44 + c55)
              L = 0.5_CUSTOM_REAL * ( cij_radial(16) + cij_radial(19) )
              ! F = 1/2 (c13 + c23)
              F = 0.5_CUSTOM_REAL * ( cij_radial(3) + cij_radial(8) )
              ! eta = F / (A - 2 L)
              eta = F / (A - 2.0_CUSTOM_REAL * L)

              ! mu & kappa for bulk kernels
              ! derived from: kappav = rho*(vpv*vpv - 4.d0/3.0*vsv*vsv) and muv = rho * vsv**2
              !               kappah = rho*(vph*vph - 4.d0/3.0*vsh*vsh) and muh = rho * vsh**2
              muvl = L
              muhl = N
              kappavl = C - 4.0/3.0 * L
              kappahl = A - 4.0/3.0 * N
              kappal = kappavl
            endif
            ! note: cijkl_kl_local() is fully anisotropic C_ij kernel components (non-dimensionalized)
            !          for GLL point at (i,j,k,ispec)

            ! Purpose : compute the kernels for the An coeffs (an_kl)
            ! from the kernels for Cij (cijkl_kl_local)
            ! At r,theta,phi fixed
            ! kernel def : dx = k_cij * dcij + k_rho * drho
            !                 = k_An * dAn  + k_rho * drho

            ! Definition of the input array cij_kl :
            ! cij_kl(1)  = C11 ; cij_kl(2)  = C12 ; cij_kl(3)  = C13
            ! cij_kl(4)  = C14 ; cij_kl(5)  = C15 ; cij_kl(6)  = C16
            ! cij_kl(7)  = C22 ; cij_kl(8)  = C23 ; cij_kl(9)  = C24
            ! cij_kl(10) = C25 ; cij_kl(11) = C26 ; cij_kl(12) = C33
            ! cij_kl(13) = C34 ; cij_kl(14) = C35 ; cij_kl(15) = C36
            ! cij_kl(16) = C44 ; cij_kl(17) = C45 ; cij_kl(18) = C46
            ! cij_kl(19) = C55 ; cij_kl(20) = C56 ; cij_kl(21) = C66
            ! where the Cij (Voigt's notation) are defined as function of
            ! the components of the elastic tensor in spherical coordinates
            ! by eq. (A.1) of Chen & Tromp, GJI 168 (2007)

            ! From the relations giving Cij in function of An
            ! Checked with Min Chen's results (routine build_cij)

            ! note: kernel with respect to anisotropic model parameters
            !       see Sieminski (2007b) and Zhu (2015)
            !
            ! kernels for model parameterization (A,C,N,L,F):
            !
            ! A = 1/8 (3 C11 + 3 C22 + 2 C12 + 4 C66)       -> kernel K_A = K_C11 + K_C12 + K_C22
            ! C = C33                                       -> kernel K_C = K_C33
            ! N = 1/8 (C11 + C22 - 2 C12 + 4 C66)           -> kernel K_N = K_C66 - 2 K_C12
            ! L = 1/2 (C44 + C55)                           -> kernel K_L = K_C44 + K_C55
            ! F = 1/2 (C13 + C23)                           -> kernel K_F = K_C13 + K_C23

            an_kl(1) = cijkl_kl_local(1) + cijkl_kl_local(2) + cijkl_kl_local(7)  ! A
            an_kl(2) = cijkl_kl_local(12)                                         ! C
            an_kl(3) = -2 * cijkl_kl_local(2) + cijkl_kl_local(21)                ! N
            an_kl(4) = cijkl_kl_local(16) + cijkl_kl_local(19)                    ! L
            an_kl(5) = cijkl_kl_local(3) + cijkl_kl_local(8)                      ! F

            ! not used yet
            !
            ! additional primitive kernels for "asymptotic parameters" (Chen & Tromp 2007):
            !
            if (SAVE_AZIMUTHAL_ANISO_KL_ONLY) then
              an_kl(6)  = 2*cijkl_kl_local(5) + 2*cijkl_kl_local(10)+ 2*cijkl_kl_local(14)       !Jc
              an_kl(7)  = 2*cijkl_kl_local(4) + 2*cijkl_kl_local(9) + 2*cijkl_kl_local(13)       !Js
              an_kl(8)  = -2*cijkl_kl_local(14)                                                  !Kc
              an_kl(9)  = -2*cijkl_kl_local(13)                                                  !Ks
              an_kl(10) = -2*cijkl_kl_local(10) + cijkl_kl_local(18)                             !Mc
              an_kl(11) = 2*cijkl_kl_local(4) - cijkl_kl_local(20)                               !Ms
              an_kl(12) = cijkl_kl_local(1) - cijkl_kl_local(7)                                  !Bc
              an_kl(13) = -1./2.*(cijkl_kl_local(6) + cijkl_kl_local(11))                        !Bs
              an_kl(14) = cijkl_kl_local(3) - cijkl_kl_local(8)                                  !Hc
              an_kl(15) = -cijkl_kl_local(15)                                                    !Hs
              an_kl(16) = -cijkl_kl_local(16) + cijkl_kl_local(19)                               !Gc
              an_kl(17) = -cijkl_kl_local(17)                                                    !Gs
              an_kl(18) = cijkl_kl_local(5) - cijkl_kl_local(10) - cijkl_kl_local(18)            !Dc
              an_kl(19) = cijkl_kl_local(4) - cijkl_kl_local(9) + cijkl_kl_local(20)             !Ds
              an_kl(20) = cijkl_kl_local(1) - cijkl_kl_local(2) + cijkl_kl_local(7) - cijkl_kl_local(21)   !Ec
              an_kl(21) = -cijkl_kl_local(6) + cijkl_kl_local(11)                                !Es

              ! more parameterizations:
              !
              ! see Zhu (2015, GJI, appendix A2):
              ! - kernels for model parameterization (L, N, Gc, Gs):
              !
              ! L  = 1/2 (C44 + C55)                          -> kernel K_L  = K_C44 + K_C55
              ! N  = 1/8 (C11 + C22 - 2 C12 + 4 C66)          -> kernel K_N  = K_C66 - 2 K_C12
              ! Gc = 1/2 (C55 - C44)                          -> kernel K_Gc = K_C55 - K_C44
              ! Gs = - C45                                    -> kernel K_Gs = - K_C45
              !
              ! - kernels for model parameterization (dln(beta_v), dln(beta_h), Gc_prime, Gs_prime) (dimension-less):
              !
              ! beta_v = sqrt(L/rho)                          -> kernel K_beta_v = 2 L K_L - 4 L eta K_F
              ! beta_h = sqrt(N/rho)                          -> kernel K_beta_h = 2 N K_N
              ! Gc_prime = Gc / (rho beta_0**2)               -> kernel K_Gc_prime = rho beta_0**2 K_Gc
              ! Gs_prime = Gs / (rho beta_0**2)               -> kernel K_Gs_prime = rho beta_0**2 K_Gs
              !
              ! with beta_0 being the isotropic shear wave speed in the 1-D reference model
              !
              ! note: for azimuthal anisotropy, Gs and Gc will provide the fast axis angle \zeta = 1/2 arctan( Gs / Gc )
              !
              !       Convention here is a Cartesian reference frame (x,y,z) where x points East, y points North and z points up.
              !       And
              !         Gs = - C54    (as compared to Gs = C54 used e.g. by Montagner, 2002, Seismic Anisotropy Tomography)
              !       and
              !         angle \zeta is measured counter-clockwise from South
              !
              Gc_kl_crust_mantle(i,j,k,ispec) = -an_kl(16) * scale_kl_ani
              Gs_kl_crust_mantle(i,j,k,ispec) = -an_kl(17) * scale_kl_ani

              ! daniel todo:
              ! scaling with actual values?
              !beta_v = sqrt(L/rhol)
              !beta_h = sqrt(N/rhol)
              ! isotropic: Voigt' average vs = sqrt( (2.d0*vsv**2 + vsh**2)/3.d0 )
              !beta_0 = sqrt( (2.0 * beta_v**2 + beta_h**2)/3.0 )

              !see model_gll.f90:
              ! test scaling with arbitrary shear moduli
              ! test: choosing PREM crustal values: rho=2.6 g/cm3, vp=5.8 km/s, vs=3.2 km/s -> mu0 = 26.624 GPa
              !Gc = Gc_prime * 26.6/scale_GPa
              !Gs = Gs_prime * 26.6/scale_GPa

              ! scaling with shear moduli from 1D background reference
              !Gc_prime = Gc / (rho beta_0**2) = Gc / mu0
              !Gs_prime = Gs / (rho beta_0**2) = Gs / mu0

              mu0 = mu0store_crust_mantle(i,j,k,ispec) ! original values from 1D background reference model
              if (abs(mu0) > TINYVAL) then
                mu0 = mu0 * scale_GPa  ! scales to GPa
                Gc_prime_kl_crust_mantle(i,j,k,ispec) = Gc_kl_crust_mantle(i,j,k,ispec) * mu0
                Gs_prime_kl_crust_mantle(i,j,k,ispec) = Gs_kl_crust_mantle(i,j,k,ispec) * mu0
              else
                Gc_prime_kl_crust_mantle(i,j,k,ispec) = 0.0_CUSTOM_REAL
                Gs_prime_kl_crust_mantle(i,j,k,ispec) = 0.0_CUSTOM_REAL
              endif
              A_kl_crust_mantle(i,j,k,ispec) = -an_kl(1)*scale_kl_ani
              C_kl_crust_mantle(i,j,k,ispec) = -an_kl(2)*scale_kl_ani
              N_kl_crust_mantle(i,j,k,ispec) = -an_kl(3)*scale_kl_ani
              L_kl_crust_mantle(i,j,k,ispec) = -an_kl(4)*scale_kl_ani
              F_kl_crust_mantle(i,j,k,ispec) = -an_kl(5)*scale_kl_ani
              Jc_kl_crust_mantle(i,j,k,ispec) = -an_kl(6)*scale_kl_ani
              Kc_kl_crust_mantle(i,j,k,ispec) = -an_kl(8)*scale_kl_ani
              Mc_kl_crust_mantle(i,j,k,ispec) = -an_kl(10)*scale_kl_ani
              Bc_kl_crust_mantle(i,j,k,ispec) = -an_kl(12)*scale_kl_ani
              Hc_kl_crust_mantle(i,j,k,ispec) = -an_kl(14)*scale_kl_ani
              Ec_kl_crust_mantle(i,j,k,ispec) = -an_kl(20)*scale_kl_ani
              Dc_kl_crust_mantle(i,j,k,ispec) = -an_kl(18)*scale_kl_ani
            endif

            ! K_rho (primary kernel, for a parameterization (A,C,L,N,F,rho) )
            rhonotprime_kl_crust_mantle(i,j,k,ispec) = rhol * rho_kl_crust_mantle(i,j,k,ispec) / scale_kl_rho

            ! note: transverse isotropic kernels are calculated for ALL elements,
            !          and not just transverse elements
            !
            ! note: the kernels are for relative perturbations (delta ln (m_i) = (m_i - m_0)/m_i )
            !
            ! Gets transverse isotropic kernels
            ! (see Appendix B of Sieminski et al., GJI 171, 2007)

            ! for parameterization: ( alpha_v, alpha_h, beta_v, beta_h, eta, rho )
            ! K_alpha_v
            alphav_kl_crust_mantle(i,j,k,ispec) = 2 * C * an_kl(2)
            ! K_alpha_h
            alphah_kl_crust_mantle(i,j,k,ispec) = 2 * A * an_kl(1) + 2 * A * eta * an_kl(5)
            ! K_beta_v
            betav_kl_crust_mantle(i,j,k,ispec) = 2 * L * an_kl(4) - 4 * L * eta * an_kl(5)
            ! K_beta_h
            betah_kl_crust_mantle(i,j,k,ispec) = 2 * N * an_kl(3)
            ! K_eta
            eta_kl_crust_mantle(i,j,k,ispec) = F * an_kl(5)

            ! K_rhoprime  (for a parameterization (alpha_v, ..., rho) )
            rho_kl_crust_mantle(i,j,k,ispec) = A * an_kl(1) + C * an_kl(2) &
                                             + N * an_kl(3) + L * an_kl(4) + F * an_kl(5) &
                                             + rhonotprime_kl_crust_mantle(i,j,k,ispec)

            ! write the kernel in physical units
            rhonotprime_kl_crust_mantle(i,j,k,ispec) = - rhonotprime_kl_crust_mantle(i,j,k,ispec) * scale_kl

            alphav_kl_crust_mantle(i,j,k,ispec) = - alphav_kl_crust_mantle(i,j,k,ispec) * scale_kl
            alphah_kl_crust_mantle(i,j,k,ispec) = - alphah_kl_crust_mantle(i,j,k,ispec) * scale_kl
            betav_kl_crust_mantle(i,j,k,ispec) = - betav_kl_crust_mantle(i,j,k,ispec) * scale_kl
            betah_kl_crust_mantle(i,j,k,ispec) = - betah_kl_crust_mantle(i,j,k,ispec) * scale_kl
            eta_kl_crust_mantle(i,j,k,ispec) = - eta_kl_crust_mantle(i,j,k,ispec) * scale_kl
            rho_kl_crust_mantle(i,j,k,ispec) = - rho_kl_crust_mantle(i,j,k,ispec) * scale_kl

            ! for parameterization: ( bulk, beta_v, beta_h, eta, rho )
            ! where kappa_v = kappa_h = kappa and bulk c = sqrt( kappa / rho )
            betav_sq = muvl / rhol
            betah_sq = muhl / rhol
            alphav_sq = ( kappal + FOUR_THIRDS * muvl ) / rhol
            alphah_sq = ( kappal + FOUR_THIRDS * muhl ) / rhol
            bulk_sq = kappal / rhol

            bulk_c_kl_crust_mantle(i,j,k,ispec) = &
              bulk_sq / alphav_sq * alphav_kl_crust_mantle(i,j,k,ispec) + &
              bulk_sq / alphah_sq * alphah_kl_crust_mantle(i,j,k,ispec)

            bulk_betah_kl_crust_mantle(i,j,k,ispec) = &
              betah_kl_crust_mantle(i,j,k,ispec) + &
              FOUR_THIRDS * betah_sq / alphah_sq * alphah_kl_crust_mantle(i,j,k,ispec)

            bulk_betav_kl_crust_mantle(i,j,k,ispec) = &
              betav_kl_crust_mantle(i,j,k,ispec) + &
              FOUR_THIRDS * betav_sq / alphav_sq * alphav_kl_crust_mantle(i,j,k,ispec)
            ! the rest, K_eta and K_rho are the same as above

            ! to check: isotropic kernels from transverse isotropic ones
            alpha_kl_crust_mantle(i,j,k,ispec) = alphav_kl_crust_mantle(i,j,k,ispec) &
                                                + alphah_kl_crust_mantle(i,j,k,ispec)
            beta_kl_crust_mantle(i,j,k,ispec) = betav_kl_crust_mantle(i,j,k,ispec) &
                                                + betah_kl_crust_mantle(i,j,k,ispec)
            !rho_kl_crust_mantle(i,j,k,ispec) = rhonotprime_kl_crust_mantle(i,j,k,ispec) &
            !                                    + alpha_kl_crust_mantle(i,j,k,ispec) &
            !                                    + beta_kl_crust_mantle(i,j,k,ispec)
            bulk_beta_kl_crust_mantle(i,j,k,ispec) = bulk_betah_kl_crust_mantle(i,j,k,ispec) &
                                                  + bulk_betav_kl_crust_mantle(i,j,k,ispec)

            ! to check: Sieminski, 2007
            !
            ! isotropic kernels
            ! K_alpha = 2 rho alpha**2 ( K_A + K_C + K_F)
            ! K_beta = 2 rho beta**2 ( K_L + K_N - 2 K_F)

          else

            ! fully anisotropic kernels

            ! note: the C_ij and density kernels are not for relative perturbations (delta ln( m_i) = delta m_i / m_i),
            !          but absolute perturbations (delta m_i = m_i - m_0)
            rho_kl_crust_mantle(i,j,k,ispec) = - rho_kl_crust_mantle(i,j,k,ispec)
            cijkl_kl_crust_mantle(:,i,j,k,ispec) = - cijkl_kl_crust_mantle(:,i,j,k,ispec)

          endif ! SAVE_TRANSVERSE_KL_ONLY .or. SAVE_AZIMUTHAL_ANISO_KL_ONLY

        enddo
      enddo
    enddo
  enddo

  ! writes out kernels to files
  if (ADIOS_FOR_KERNELS) then
    call write_kernels_cm_ani_adios(alphav_kl_crust_mantle,alphah_kl_crust_mantle, &
                                    betav_kl_crust_mantle,betah_kl_crust_mantle, &
                                    eta_kl_crust_mantle, &
                                    bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle, &
                                    bulk_betav_kl_crust_mantle,bulk_betah_kl_crust_mantle, &
                                    Gc_prime_kl_crust_mantle, Gs_prime_kl_crust_mantle)
  else
    ! binary file output
    call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_TMP_PATH)

    ! For anisotropic kernels

    ! outputs transverse isotropic kernels only
    if (SAVE_TRANSVERSE_KL_ONLY) then
      ! transverse isotropic kernels
      ! (alpha_v, alpha_h, beta_v, beta_h, eta, rho ) parameterization
      open(unit=IOUT,file=trim(prname)//'alphav_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) alphav_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'alphah_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) alphah_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'betav_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) betav_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'betah_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) betah_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'eta_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) eta_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'rho_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) rho_kl_crust_mantle
      close(IOUT)

      ! in case one is interested in primary kernel K_rho
      !open(unit=IOUT,file=trim(prname)//'rhonotprime_kernel.bin',status='unknown',form='unformatted',action='write')
      !write(IOUT) rhonotprime_kl_crust_mantle
      !close(IOUT)

      ! (bulk, beta_v, beta_h, eta, rho ) parameterization: K_eta and K_rho same as above
      open(unit=IOUT,file=trim(prname)//'bulk_c_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) bulk_c_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'bulk_betav_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) bulk_betav_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'bulk_betah_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) bulk_betah_kl_crust_mantle
      close(IOUT)

      ! to check: isotropic kernels
      open(unit=IOUT,file=trim(prname)//'alpha_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) alpha_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'beta_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) beta_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'bulk_beta_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) bulk_beta_kl_crust_mantle
      close(IOUT)

    else if (SAVE_AZIMUTHAL_ANISO_KL_ONLY) then
      ! kernels for inversions involving azimuthal anisotropy
      ! (alpha_v, alpha_h, beta_v, beta_h, eta, rho ) parameterization
      open(unit=IOUT,file=trim(prname)//'alphav_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) alphav_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'alphah_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) alphah_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'betav_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) betav_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'betah_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) betah_kl_crust_mantle
      close(IOUT)

      ! (bulk_c, beta_v, beta_h, eta, Gc', Gs', rho ) parameterization
      open(unit=IOUT,file=trim(prname)//'bulk_c_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) bulk_c_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'bulk_betav_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) bulk_betav_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'bulk_betah_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) bulk_betah_kl_crust_mantle
      close(IOUT)

      open(unit=IOUT,file=trim(prname)//'eta_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) eta_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'rho_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) rho_kl_crust_mantle
      close(IOUT)

      ! note: Gc' & Gs' are the normalized Gc & Gs kernels
      open(unit=IOUT,file=trim(prname)//'Gc_prime_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) Gc_prime_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'Gs_prime_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) Gs_prime_kl_crust_mantle
      close(IOUT)

      ! to check: isotropic kernels
      if (.false.) then
        open(unit=IOUT,file=trim(prname)//'alpha_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) alpha_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'beta_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) beta_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'bulk_beta_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) bulk_beta_kl_crust_mantle
        close(IOUT)
      endif

      ! to check: all anisotropic kernels
      if (.false.) then
        open(unit=IOUT,file=trim(prname)//'A_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) A_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'C_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) C_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'L_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) L_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'N_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) N_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'F_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) F_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'Gc_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) Gc_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'Gs_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) Gs_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'Jc_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) Jc_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'Kc_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) Kc_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'Mc_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) Mc_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'Bc_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) Bc_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'Hc_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) Hc_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'Ec_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) Ec_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'Dc_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) Dc_kl_crust_mantle
        close(IOUT)
      endif

    else

      ! fully anisotropic kernels
      open(unit=IOUT,file=trim(prname)//'rho_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) rho_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'cijkl_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) cijkl_kl_crust_mantle
      close(IOUT)

    endif

  endif ! ADIOS_FOR_KERNELS

  ! Output these kernels as netcdf files -- one per processor.
#ifdef USE_CEM
  if (SAVE_TRANSVERSE_KL_ONLY) then
    filename = trim(OUTPUT_FILES)//'/alphavKernelCrustMantle.nc'
    call write_kernel_netcdf(filename, alphav_kl_crust_mantle)

    filename = trim(OUTPUT_FILES)//'/alphahKernelCrustMantle.nc'
    call write_kernel_netcdf(filename, alphah_kl_crust_mantle)

    filename = trim(OUTPUT_FILES)//'/betavKernelCrustMantle.nc'
    call write_kernel_netcdf(filename,betav_kl_crust_mantle)

    filename = trim(OUTPUT_FILES)//'/betahKernelCrustMantle.nc'
    call write_kernel_netcdf(filename,betah_kl_crust_mantle)

    filename = trim(OUTPUT_FILES)//'/etaKernelCrustMantle.nc'
    call write_kernel_netcdf(filename,eta_kl_crust_mantle)

    filename = trim(OUTPUT_FILES)//'/rhoKernelCrustMantle.nc'
    call write_kernel_netcdf(filename,rho_kl_crust_mantle)

    filename = trim(OUTPUT_FILES)//'/xyzCrustMantle.nc'
    call write_coordinates_netcdf(filename)

  else if (SAVE_AZIMUTHAL_ANISO_KL_ONLY) then
    filename = trim(OUTPUT_FILES)//'/alphavKernelCrustMantle.nc'
    call write_kernel_netcdf(filename, alphav_kl_crust_mantle)

    filename = trim(OUTPUT_FILES)//'/alphahKernelCrustMantle.nc'
    call write_kernel_netcdf(filename, alphah_kl_crust_mantle)

    filename = trim(OUTPUT_FILES)//'/betavKernelCrustMantle.nc'
    call write_kernel_netcdf(filename,betav_kl_crust_mantle)

    filename = trim(OUTPUT_FILES)//'/betahKernelCrustMantle.nc'
    call write_kernel_netcdf(filename,betah_kl_crust_mantle)

    filename = trim(OUTPUT_FILES)//'/etaKernelCrustMantle.nc'
    call write_kernel_netcdf(filename,eta_kl_crust_mantle)

    filename = trim(OUTPUT_FILES)//'/rhoKernelCrustMantle.nc'
    call write_kernel_netcdf(filename,rho_kl_crust_mantle)

    filename = trim(OUTPUT_FILES)//'/GcprimeKernelCrustMantle.nc'
    call write_kernel_netcdf(filename,Gc_prime_kl_crust_mantle)

    filename = trim(OUTPUT_FILES)//'/GsprimeKernelCrustMantle.nc'
    call write_kernel_netcdf(filename,Gs_prime_kl_crust_mantle)

    filename = trim(OUTPUT_FILES)//'/xyzCrustMantle.nc'
    call write_coordinates_netcdf(filename)
  endif
#endif

  ! cleans up temporary kernel arrays
  deallocate(alphav_kl_crust_mantle,alphah_kl_crust_mantle, &
             betav_kl_crust_mantle,betah_kl_crust_mantle, &
             eta_kl_crust_mantle)
  deallocate(bulk_betah_kl_crust_mantle, &
             bulk_betav_kl_crust_mantle)
  deallocate(bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle)
  deallocate(rhonotprime_kl_crust_mantle)
  deallocate(Gc_prime_kl_crust_mantle,Gs_prime_kl_crust_mantle)

  if (SAVE_AZIMUTHAL_ANISO_KL_ONLY) then
    deallocate(Gc_kl_crust_mantle,Gs_kl_crust_mantle)
    deallocate(A_kl_crust_mantle,C_kl_crust_mantle,F_kl_crust_mantle, &
               L_kl_crust_mantle,N_kl_crust_mantle,Jc_kl_crust_mantle, &
               Kc_kl_crust_mantle,Mc_kl_crust_mantle,Bc_kl_crust_mantle, &
               Hc_kl_crust_mantle,Ec_kl_crust_mantle,Dc_kl_crust_mantle)
  endif

  end subroutine save_kernels_crust_mantle_ani

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_kernels_crust_mantle_iso()

! stores kernels for isotropic parameterization

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: scale_kl,scale_kl_ani,scale_kl_rho
  real(kind=CUSTOM_REAL) :: rhol,mul,kappal
  real(kind=CUSTOM_REAL) :: rho_kl,alpha_kl,beta_kl
  integer :: ispec,i,j,k
  integer :: ier

  ! primary kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: &
    mu_kl_crust_mantle, kappa_kl_crust_mantle, rhonotprime_kl_crust_mantle

  ! bulk parameterization
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: &
    bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle

  ! checks if anything to do
  if (ANISOTROPIC_KL) return

  ! scaling factors: note that this scaling has been introduced by Qinya Liu (2006)
  !                  with the intent to dimensionalize kernel values to [ s km^(-3) ]
  !
  ! kernel unit [ s / km^3 ]
  scale_kl = real(scale_t * scale_displ_inv * 1.d9,kind=CUSTOM_REAL)
  ! For anisotropic kernels
  ! final unit : [s km^(-3) GPa^(-1)]
  scale_kl_ani = real(scale_t**3 / (RHOAV*R_PLANET**3) * 1.d18,kind=CUSTOM_REAL)
  ! final unit : [s km^(-3) (kg/m^3)^(-1)]
  scale_kl_rho = real(scale_t * scale_displ_inv / RHOAV * 1.d9,kind=CUSTOM_REAL)

  ! isotropic kernels
  !
  ! allocates temporary arrays
  ! primary kernels
  allocate(mu_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
           kappa_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
           rhonotprime_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
  if (ier /= 0 ) stop 'Error allocating transverse kernels bulk_c_kl_crust_mantle,...'
  mu_kl_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  kappa_kl_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  rhonotprime_kl_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL

  ! bulk velocity kernels
  allocate(bulk_c_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
           bulk_beta_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
  if (ier /= 0 ) stop 'Error allocating transverse kernels bulk_c_kl_crust_mantle,...'
  bulk_c_kl_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  bulk_beta_kl_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL

  ! crust_mantle
  do ispec = 1, NSPEC_CRUST_MANTLE_ADJOINT
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          ! gets material properties
          rhol = rhostore_crust_mantle(i,j,k,ispec)
          mul = muvstore_crust_mantle(i,j,k,ispec)
          kappal = kappavstore_crust_mantle(i,j,k,ispec)

          ! kernel values for rho, kappa, mu (primary kernel values)
          rho_kl = - rhol * rho_kl_crust_mantle(i,j,k,ispec)
          alpha_kl = - kappal * alpha_kl_crust_mantle(i,j,k,ispec) ! note: alpha_kl corresponds to K_kappa
          beta_kl =  - 2 * mul * beta_kl_crust_mantle(i,j,k,ispec) ! note: beta_kl corresponds to K_mu

          ! for a parameterization: (rho,mu,kappa) "primary" kernels
          rhonotprime_kl_crust_mantle(i,j,k,ispec) = rho_kl * scale_kl
          mu_kl_crust_mantle(i,j,k,ispec) = beta_kl * scale_kl
          kappa_kl_crust_mantle(i,j,k,ispec) = alpha_kl * scale_kl

          ! for a parameterization: (rho,alpha,beta)
          ! kernels rho^prime, beta, alpha
          rho_kl_crust_mantle(i,j,k,ispec) = (rho_kl + alpha_kl + beta_kl) * scale_kl
          beta_kl_crust_mantle(i,j,k,ispec) = &
            2._CUSTOM_REAL * (beta_kl - FOUR_THIRDS * mul / kappal * alpha_kl) * scale_kl
          alpha_kl_crust_mantle(i,j,k,ispec) = &
            2._CUSTOM_REAL * (1 +  FOUR_THIRDS * mul / kappal) * alpha_kl * scale_kl

          ! for a parameterization: (rho,bulk, beta)
          ! where bulk wave speed is c = sqrt( kappa / rho)
          ! note: rhoprime is the same as for (rho,alpha,beta) parameterization
          bulk_c_kl_crust_mantle(i,j,k,ispec) = 2._CUSTOM_REAL * alpha_kl * scale_kl
          bulk_beta_kl_crust_mantle(i,j,k,ispec ) = 2._CUSTOM_REAL * beta_kl * scale_kl

        enddo
      enddo
    enddo
  enddo

  ! writes out kernels to files
  if (ADIOS_FOR_KERNELS) then
    call write_kernels_cm_iso_adios(mu_kl_crust_mantle, kappa_kl_crust_mantle, rhonotprime_kl_crust_mantle, &
                                    bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle)
  else

    call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_TMP_PATH)

    ! primary kernels: (rho,kappa,mu) parameterization
    open(unit=IOUT,file=trim(prname)//'rhonotprime_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) rhonotprime_kl_crust_mantle
    close(IOUT)
    open(unit=IOUT,file=trim(prname)//'kappa_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) kappa_kl_crust_mantle
    close(IOUT)
    open(unit=IOUT,file=trim(prname)//'mu_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) mu_kl_crust_mantle
    close(IOUT)

    ! (rho, alpha, beta ) parameterization
    open(unit=IOUT,file=trim(prname)//'rho_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) rho_kl_crust_mantle
    close(IOUT)
    open(unit=IOUT,file=trim(prname)//'alpha_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) alpha_kl_crust_mantle
    close(IOUT)
    open(unit=IOUT,file=trim(prname)//'beta_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) beta_kl_crust_mantle
    close(IOUT)

    ! (rho, bulk, beta ) parameterization, K_rho same as above
    open(unit=IOUT,file=trim(prname)//'bulk_c_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) bulk_c_kl_crust_mantle
    close(IOUT)
    open(unit=IOUT,file=trim(prname)//'bulk_beta_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) bulk_beta_kl_crust_mantle
    close(IOUT)

  endif ! ADIOS_FOR_KERNELS

  ! cleans up temporary kernel arrays
  deallocate(bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle)
  deallocate(mu_kl_crust_mantle,kappa_kl_crust_mantle,rhonotprime_kl_crust_mantle)

  end subroutine save_kernels_crust_mantle_iso



!
!-------------------------------------------------------------------------------------------------
!

! put the list of parameters back here to avoid a warning / error from the gfortran compiler
! about undefined behavior when aggressive loop vectorization is used by the compiler

  subroutine save_kernels_outer_core(rhostore_outer_core,kappavstore_outer_core,rho_kl_outer_core,alpha_kl_outer_core)

  use specfem_par

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),intent(in) :: &
    rhostore_outer_core,kappavstore_outer_core

  ! adjoint kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT),intent(inout) :: &
    rho_kl_outer_core,alpha_kl_outer_core

  ! local parameters
  real(kind=CUSTOM_REAL):: scale_kl
  real(kind=CUSTOM_REAL) :: rhol,kappal,rho_kl,alpha_kl
  integer :: ispec,i,j,k

  ! saftey check
  if (.not. SAVE_KERNELS_OC) return

  scale_kl = real(scale_t * scale_displ_inv * 1.d9,kind=CUSTOM_REAL)

  ! outer_core
  do ispec = 1, NSPEC_OUTER_CORE_ADJOINT
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          rhol = rhostore_outer_core(i,j,k,ispec)
          kappal = kappavstore_outer_core(i,j,k,ispec)

          rho_kl = - rhol * rho_kl_outer_core(i,j,k,ispec)
          alpha_kl = - kappal * alpha_kl_outer_core(i,j,k,ispec)

          rho_kl_outer_core(i,j,k,ispec) = (rho_kl + alpha_kl) * scale_kl
          alpha_kl_outer_core(i,j,k,ispec) = 2 * alpha_kl * scale_kl
        enddo
      enddo
    enddo
  enddo

  ! writes out kernels to file
  if (ADIOS_FOR_KERNELS) then
    call write_kernels_oc_adios()
  else
    call create_name_database(prname,myrank,IREGION_OUTER_CORE,LOCAL_TMP_PATH)

    open(unit=IOUT,file=trim(prname)//'rho_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) rho_kl_outer_core
    close(IOUT)
    open(unit=IOUT,file=trim(prname)//'alpha_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) alpha_kl_outer_core
    close(IOUT)

  endif

  end subroutine save_kernels_outer_core

!
!-------------------------------------------------------------------------------------------------
!

! put the list of parameters back here to avoid a warning / error from the gfortran compiler
! about undefined behavior when aggressive loop vectorization is used by the compiler

  subroutine save_kernels_inner_core(rhostore_inner_core,muvstore_inner_core,kappavstore_inner_core, &
                                     rho_kl_inner_core,alpha_kl_inner_core,beta_kl_inner_core)

  use specfem_par

  implicit none

  ! material parameters
  ! (note: muvstore also needed for attenuation in case of anisotropic inner core)
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),intent(in) :: &
    rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core

  ! adjoint kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT),intent(inout) :: &
    rho_kl_inner_core,beta_kl_inner_core, alpha_kl_inner_core

  ! local parameters
  real(kind=CUSTOM_REAL):: scale_kl
  real(kind=CUSTOM_REAL) :: rhol,mul,kappal,rho_kl,alpha_kl,beta_kl
  integer :: ispec,i,j,k

  ! safety check
  if (.not. SAVE_KERNELS_IC) return

  ! scaling to units
  scale_kl = real(scale_t * scale_displ_inv * 1.d9,kind=CUSTOM_REAL)

  ! inner_core
  do ispec = 1, NSPEC_INNER_CORE_ADJOINT
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          rhol = rhostore_inner_core(i,j,k,ispec)
          mul = muvstore_inner_core(i,j,k,ispec)
          kappal = kappavstore_inner_core(i,j,k,ispec)

          rho_kl = -rhol * rho_kl_inner_core(i,j,k,ispec)
          alpha_kl = -kappal * alpha_kl_inner_core(i,j,k,ispec)
          beta_kl =  - 2._CUSTOM_REAL * mul * beta_kl_inner_core(i,j,k,ispec)

          rho_kl_inner_core(i,j,k,ispec) = (rho_kl + alpha_kl + beta_kl) * scale_kl
          beta_kl_inner_core(i,j,k,ispec) = 2._CUSTOM_REAL * (beta_kl - FOUR_THIRDS * mul * alpha_kl / kappal) * scale_kl
          alpha_kl_inner_core(i,j,k,ispec) = 2._CUSTOM_REAL * (1._CUSTOM_REAL +  FOUR_THIRDS * mul / kappal) * alpha_kl * scale_kl
        enddo
      enddo
    enddo
  enddo

  ! writes out kernels to file
  if (ADIOS_FOR_KERNELS) then
    call write_kernels_ic_adios()
  else
    call create_name_database(prname,myrank,IREGION_INNER_CORE,LOCAL_TMP_PATH)

    open(unit=IOUT,file=trim(prname)//'rho_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) rho_kl_inner_core
    close(IOUT)
    open(unit=IOUT,file=trim(prname)//'alpha_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) alpha_kl_inner_core
    close(IOUT)
    open(unit=IOUT,file=trim(prname)//'beta_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) beta_kl_inner_core
    close(IOUT)
  endif

  end subroutine save_kernels_inner_core

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_kernels_boundary_kl()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL):: scale_kl

  ! saftey check
  if (.not. SAVE_KERNELS_BOUNDARY) return

  ! kernel unit [ s / km^3 ]
  scale_kl = real(scale_t * scale_displ_inv * 1.d9,kind=CUSTOM_REAL)

  ! scale the boundary kernels properly: *scale_kl gives s/km^3 and 1.d3 gives
  ! the relative boundary kernels (for every 1 km) in s/km^2
  moho_kl(:,:,:) = moho_kl(:,:,:) * real(scale_kl * 1.d3,kind=CUSTOM_REAL)
  d400_kl(:,:,:) = d400_kl(:,:,:) * real(scale_kl * 1.d3,kind=CUSTOM_REAL)
  d670_kl(:,:,:) = d670_kl(:,:,:) * real(scale_kl * 1.d3,kind=CUSTOM_REAL)
  cmb_kl(:,:,:) = cmb_kl(:,:,:) * real(scale_kl * 1.d3,kind=CUSTOM_REAL)
  icb_kl(:,:,:) = icb_kl(:,:,:) * real(scale_kl * 1.d3,kind=CUSTOM_REAL)

  ! writes out kernels to file
  if (ADIOS_FOR_KERNELS) then
    call write_kernels_boundary_kl_adios()
  else
    call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_TMP_PATH)

    if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
      open(unit=IOUT,file=trim(prname)//'moho_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) moho_kl
      close(IOUT)
    endif

    open(unit=IOUT,file=trim(prname)//'d400_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) d400_kl
    close(IOUT)

    open(unit=IOUT,file=trim(prname)//'d670_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) d670_kl
    close(IOUT)

    open(unit=IOUT,file=trim(prname)//'CMB_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) cmb_kl
    close(IOUT)

    call create_name_database(prname,myrank,IREGION_OUTER_CORE,LOCAL_PATH)

    open(unit=IOUT,file=trim(prname)//'ICB_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) icb_kl
    close(IOUT)
  endif

  end subroutine save_kernels_boundary_kl

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_kernels_source_derivatives()

  use specfem_par

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: scale_mass
  integer :: irec_local
  character(len=MAX_STRING_LEN) :: outputname

  ! scaling factor
  scale_mass = real(RHOAV * (R_EARTH**3),kind=CUSTOM_REAL)

  ! computes derivatives
  do irec_local = 1, nrec_local
    ! rotate and scale the location derivatives to correspond to dn,de,dz
    sloc_der(:,irec_local) = real(matmul(transpose(nu_source(:,:,irec_local)),sloc_der(:,irec_local)) &
                                  * scale_displ * scale_t,kind=CUSTOM_REAL)

    ! rotate scale the moment derivatives to correspond to M[n,e,z][n,e,z]
    moment_der(:,:,irec_local) = real(matmul(matmul(transpose(nu_source(:,:,irec_local)),moment_der(:,:,irec_local)), &
               nu_source(:,:,irec_local)) * scale_t ** 3 / scale_mass,kind=CUSTOM_REAL)

    ! *nu_source* is the rotation matrix from ECEF to local N-E-UP as defined in src/specfem3D/locate_sources.f90

! From Qinya Liu, Toronto University, Canada:
! these derivatives are basically derivatives of the misfit function phi with respect to
! source parameters, which means, if the nu is the rotation matrix that
! transforms coordinates from the global system (x,y,z) to the local
! coordinate system (N,E,V), e.g., the moment tensor is transformed as
!
! M_L = \nu * M_g * \nu^T,
!
! then the derivative should be transformed as
!
! \partial{\phi}{M_L} = \nu^T \partial{\phi}{M_g} \nu
!
! which is in the opposite sense from the transformation of M.

    ! derivatives for time shift and hduration
    stshift_der(irec_local) = stshift_der(irec_local) * real(scale_displ**2,kind=CUSTOM_REAL)
    shdur_der(irec_local) = shdur_der(irec_local) * real(scale_displ**2,kind=CUSTOM_REAL)
  enddo

  ! writes out kernels to file
  if (ADIOS_FOR_KERNELS) then
    call write_kernels_source_derivatives_adios()
  else
    ! kernel file output
    do irec_local = 1, nrec_local
      write(outputname,'(a,i6.6)') trim(OUTPUT_FILES)//'/src_frechet.',number_receiver_global(irec_local)
      open(unit=IOUT,file=trim(outputname),status='unknown',action='write')
      !
      ! r -> z, theta -> -n, phi -> e, plus factor 2 for Mrt,Mrp,Mtp, and 1e-7 to dyne.cm
      !  Mrr =  Mzz
      !  Mtt =  Mnn
      !  Mpp =  Mee
      !  Mrt = -Mzn
      !  Mrp =  Mze
      !  Mtp = -Mne
      ! for consistency, location derivatives are in the order of [Xr,Xt,Xp]
      ! minus sign for sloc_der(3,irec_local) to get derivative for depth instead of radius
      write(IOUT,'(g16.5)') moment_der(3,3,irec_local) * 1e-7
      write(IOUT,'(g16.5)') moment_der(1,1,irec_local) * 1e-7
      write(IOUT,'(g16.5)') moment_der(2,2,irec_local) * 1e-7
      write(IOUT,'(g16.5)') -2*moment_der(1,3,irec_local) * 1e-7
      write(IOUT,'(g16.5)') 2*moment_der(2,3,irec_local) * 1e-7
      write(IOUT,'(g16.5)') -2*moment_der(1,2,irec_local) * 1e-7
      write(IOUT,'(g16.5)') sloc_der(2,irec_local)
      write(IOUT,'(g16.5)') sloc_der(1,irec_local)
      write(IOUT,'(g16.5)') -sloc_der(3,irec_local)

      write(IOUT,'(g16.5)') stshift_der(irec_local)
      write(IOUT,'(g16.5)') shdur_der(irec_local)

      close(IOUT)
    enddo
  endif

  end subroutine save_kernels_source_derivatives

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_kernels_Hessian()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: scale_kl

  ! scaling factors
  scale_kl = real(scale_t * scale_displ_inv * 1.d9,kind=CUSTOM_REAL)

  ! scales approximate Hessian
  hess_kl_crust_mantle(:,:,:,:) = 2._CUSTOM_REAL * hess_kl_crust_mantle(:,:,:,:) * scale_kl

  ! writes out kernels to file
  if (ADIOS_FOR_KERNELS) then
    call write_kernels_Hessian_adios()
  else
    ! stores into file
    call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_TMP_PATH)

    open(unit=IOUT,file=trim(prname)//'hess_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) hess_kl_crust_mantle
    close(IOUT)

    open(unit=IOUT,file=trim(prname)//'hess_rho_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) hess_rho_kl_crust_mantle
    close(IOUT)

    open(unit=IOUT,file=trim(prname)//'hess_kappa_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) hess_kappa_kl_crust_mantle
    close(IOUT)

    open(unit=IOUT,file=trim(prname)//'hess_mu_kernel.bin',status='unknown',form='unformatted',action='write')
    write(IOUT) hess_mu_kl_crust_mantle
    close(IOUT)

  endif

  end subroutine save_kernels_Hessian


