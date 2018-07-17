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

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"



  subroutine save_kernels()

  use constants_solver, only: SAVE_BOUNDARY_MESH

  use specfem_par, only: NOISE_TOMOGRAPHY,SIMULATION_TYPE,nrec_local, &
    APPROXIMATE_HESS_KL,ADIOS_FOR_KERNELS

  use specfem_par_innercore, only: rhostore_inner_core,muvstore_inner_core,kappavstore_inner_core, &
    rho_kl_inner_core,alpha_kl_inner_core,beta_kl_inner_core

  use specfem_par_outercore, only: rhostore_outer_core,kappavstore_outer_core,rho_kl_outer_core,alpha_kl_outer_core

  use manager_adios

  implicit none

  ! Open an handler to the ADIOS file in which kernel variables are written.
  if (ADIOS_FOR_KERNELS) then
    if ((SIMULATION_TYPE == 3) .or. (SIMULATION_TYPE == 2 .and. nrec_local > 0)) &
      call define_kernel_adios_variables()
  endif

  ! dump kernel arrays
  if (SIMULATION_TYPE == 3) then

    ! restores relaxed elastic moduli after being shifted to unrelaxed values before time loop
    call restore_relaxed_moduli()

    ! crust mantle
    call save_kernels_crust_mantle()

    ! noise strength kernel
    if (NOISE_TOMOGRAPHY == 3) then
      call save_kernels_strength_noise()
    endif

    ! outer core
    call save_kernels_outer_core(rhostore_outer_core,kappavstore_outer_core,rho_kl_outer_core,alpha_kl_outer_core)

    ! inner core
    call save_kernels_inner_core(rhostore_inner_core,muvstore_inner_core,kappavstore_inner_core, &
                                     rho_kl_inner_core,alpha_kl_inner_core,beta_kl_inner_core)

    ! boundary kernel
    if (SAVE_BOUNDARY_MESH) then
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
    if ((SIMULATION_TYPE == 3) .or. (SIMULATION_TYPE == 2 .and. nrec_local > 0)) &
      call close_file_adios()
  endif

  end subroutine save_kernels

!
!-------------------------------------------------------------------------------------------------
!

  subroutine restore_relaxed_moduli()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: one_minus_sum_beta_use,mul
  integer :: ispec
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
  if (GPU_MODE) return

  ! only muvstore is used further and needs to be restored

  ! crust/mantle
  do ispec = 1,NSPEC_CRUST_MANTLE
    ! isotropic and tiso elements
    DO_LOOP_IJK

      ! layer with no transverse isotropy, use kappav and muv
      mul = muvstore_crust_mantle(INDEX_IJK,ispec)

      ! precompute terms for attenuation if needed
      if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
        one_minus_sum_beta_use = one_minus_sum_beta_crust_mantle(INDEX_IJK,ispec)
      else
        one_minus_sum_beta_use = one_minus_sum_beta_crust_mantle(1,1,1,ispec)
      endif

      ! returns to the relaxed moduli Mu_r = Mu_u / [1 - sum(1 - tau_strain/tau_stress) ]
      mul = mul / one_minus_sum_beta_use

      ! stores relaxed shear moduli for kernel computations
      muvstore_crust_mantle(INDEX_IJK,ispec) = mul
    ENDDO_LOOP_IJK
  enddo

  ! inner core
  do ispec = 1,NSPEC_INNER_CORE
    ! exclude fictitious elements in central cube
    if (idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle
    ! isotropic element
    DO_LOOP_IJK

      ! layer with no transverse isotropy, use kappav and muv
      mul = muvstore_inner_core(INDEX_IJK,ispec)

      ! precompute terms for attenuation if needed
      if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
        one_minus_sum_beta_use = one_minus_sum_beta_inner_core(INDEX_IJK,ispec)
      else
        one_minus_sum_beta_use = one_minus_sum_beta_inner_core(1,1,1,ispec)
      endif

      ! returns to the relaxed moduli Mu_r = Mu_u / [1 - sum(1 - tau_strain/tau_stress) ]
      mul = mul / one_minus_sum_beta_use

      ! stores relaxed shear moduli for kernel computations
      muvstore_inner_core(INDEX_IJK,ispec) = mul
    ENDDO_LOOP_IJK
  enddo

  end subroutine restore_relaxed_moduli

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_kernels_crust_mantle()

  use specfem_par, only: SAVE_REGULAR_KL,ANISOTROPIC_KL

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
  real(kind=CUSTOM_REAL) :: scale_kl,scale_kl_ani,scale_kl_rho
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
  real(kind=CUSTOM_REAL) :: A,C,F,L,N,eta
  real(kind=CUSTOM_REAL) :: muvl,kappavl,muhl,kappahl
  real(kind=CUSTOM_REAL) :: alphav_sq,alphah_sq,betav_sq,betah_sq,bulk_sq

#ifdef CEM
  character(len=MAX_STRING_LEN) :: filename
#endif

  ! checks if anything to do
  if (.not. ANISOTROPIC_KL) return

  ! scaling factors: note that this scaling has been introduced by Qinya Liu (2006)
  !                  with the intent to dimensionalize kernel values to [ s km^(-3) ]
  !
  ! kernel unit [ s / km^3 ]
  scale_kl = scale_t * scale_displ_inv * 1.d9
  ! For anisotropic kernels
  ! final unit : [s km^(-3) GPa^(-1)]
  scale_kl_ani = scale_t**3 / (RHOAV*R_EARTH**3) * 1.d18
  ! final unit : [s km^(-3) (kg/m^3)^(-1)]
  scale_kl_rho = scale_t * scale_displ_inv / RHOAV * 1.d9

  ! allocates temporary arrays
  if (SAVE_TRANSVERSE_KL_ONLY) then
    ! transverse isotropic kernel arrays for file output
    allocate(alphav_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             alphah_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             betav_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             betah_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             eta_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
    if (ier /= 0 ) stop 'Error allocating transverse kernels alphav_kl_crust_mantle,...'

    ! isotropic kernel arrays for file output
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
  endif

  allocate(rhonotprime_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
  if (ier /= 0 ) stop 'Error allocating transverse kernels rhonotprime_kl_crust_mantle'

  ! crust_mantle
  do ispec = 1, NSPEC_CRUST_MANTLE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

          ! For anisotropic kernels
          iglob = ibool_crust_mantle(i,j,k,ispec)

          ! The Cartesian global cijkl_kl are rotated into the spherical local cijkl_kl
          ! ystore and zstore are thetaval and phival (line 2252) -- dangerous
          theta = rstore_crust_mantle(2,iglob)
          phi = rstore_crust_mantle(3,iglob)

          call rotate_kernels_dble(cijkl_kl_crust_mantle(:,i,j,k,ispec),cijkl_kl_local(:),theta,phi)

          cijkl_kl_crust_mantle(:,i,j,k,ispec) = cijkl_kl_local(:) * scale_kl_ani
          rho_kl_crust_mantle(i,j,k,ispec) = rho_kl_crust_mantle(i,j,k,ispec) * scale_kl_rho

          ! transverse isotropic kernel calculations
          if (SAVE_TRANSVERSE_KL_ONLY) then
            ! note: transverse isotropic kernels are calculated for all elements
            !
            !          however, the factors A,C,L,N,F are based only on transverse elements
            !          in between Moho and 220 km layer, otherwise they are taken from isotropic values

            rhol = rhostore_crust_mantle(i,j,k,ispec)

            ! transverse isotropic parameters from compute_force_crust_mantle.f90
            ! C=rhovpvsq A=rhovphsq L=rhovsvsq N=rhovshsq eta=F/(A - 2 L)

            ! Get A,C,F,L,N,eta from kappa,mu
            ! element can have transverse isotropy if between d220 and Moho
            if (.not. ispec_is_tiso_crust_mantle(ispec)) then

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
            !an_kl(6)  = 2*cijkl_kl_local(5)+2*cijkl_kl_local(10)+2*cijkl_kl_local(14)          !Jc
            !an_kl(7)  = 2*cijkl_kl_local(4)+2*cijkl_kl_local(9)+2*cijkl_kl_local(13)           !Js
            !an_kl(8)  = -2*cijkl_kl_local(14)                                                  !Kc
            !an_kl(9)  = -2*cijkl_kl_local(13)                                                  !Ks
            !an_kl(10) = -2*cijkl_kl_local(10)+cijkl_kl_local(18)                               !Mc
            !an_kl(11) = 2*cijkl_kl_local(4)-cijkl_kl_local(20)                                 !Ms
            !an_kl(12) = cijkl_kl_local(1)-cijkl_kl_local(7)                                    !Bc
            !an_kl(13) = -1./2.*(cijkl_kl_local(6)+cijkl_kl_local(11))                          !Bs
            !an_kl(14) = cijkl_kl_local(3)-cijkl_kl_local(8)                                    !Hc
            !an_kl(15) = -cijkl_kl_local(15)                                                    !Hs
            !an_kl(16) = -cijkl_kl_local(16)+cijkl_kl_local(19)                                 !Gc
            !an_kl(17) = -cijkl_kl_local(17)                                                    !Gs
            !an_kl(18) = cijkl_kl_local(5)-cijkl_kl_local(10)-cijkl_kl_local(18)                !Dc
            !an_kl(19) = cijkl_kl_local(4)-cijkl_kl_local(9)+cijkl_kl_local(20)                 !Ds
            !an_kl(20) = cijkl_kl_local(1)-cijkl_kl_local(2)+cijkl_kl_local(7)-cijkl_kl_local(21)      !Ec
            !an_kl(21) = -cijkl_kl_local(6)+cijkl_kl_local(11)                                  !Es

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

          endif ! SAVE_TRANSVERSE_KL_ONLY

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
                                    bulk_betav_kl_crust_mantle,bulk_betah_kl_crust_mantle)
  else

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

    else

      ! fully anisotropic kernels
      ! note: the C_ij and density kernels are not for relative perturbations (delta ln( m_i) = delta m_i / m_i),
      !          but absolute perturbations (delta m_i = m_i - m_0)
      open(unit=IOUT,file=trim(prname)//'rho_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) - rho_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'cijkl_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) - cijkl_kl_crust_mantle
      close(IOUT)

    endif

  endif ! ADIOS_FOR_KERNELS

  ! Output these kernels as netcdf files -- one per processor.
#ifdef CEM
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
  scale_kl = scale_t * scale_displ_inv * 1.d9
  ! For anisotropic kernels
  ! final unit : [s km^(-3) GPa^(-1)]
  scale_kl_ani = scale_t**3 / (RHOAV*R_EARTH**3) * 1.d18
  ! final unit : [s km^(-3) (kg/m^3)^(-1)]
  scale_kl_rho = scale_t * scale_displ_inv / RHOAV * 1.d9


  ! isotropic kernels
  !
  ! allocates temporary arrays
  ! primary kernels
  allocate(mu_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
           kappa_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
           rhonotprime_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
  if (ier /= 0 ) stop 'Error allocating transverse kernels bulk_c_kl_crust_mantle,...'

  ! bulk velocity kernels
  allocate(bulk_c_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
           bulk_beta_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
  if (ier /= 0 ) stop 'Error allocating transverse kernels bulk_c_kl_crust_mantle,...'

  ! crust_mantle
  do ispec = 1, NSPEC_CRUST_MANTLE
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
            2._CUSTOM_REAL * (beta_kl - FOUR_THIRDS * mul * alpha_kl / kappal) * scale_kl
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

!! DK DK put the list of parameters back here to avoid a warning / error from the gfortran compiler
!! DK DK about undefined behavior when aggressive loop vectorization is used by the compiler
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

  scale_kl = scale_t * scale_displ_inv * 1.d9

  ! outer_core
  do ispec = 1, NSPEC_OUTER_CORE
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

!! DK DK put the list of parameters back here to avoid a warning / error from the gfortran compiler
!! DK DK about undefined behavior when aggressive loop vectorization is used by the compiler
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

  ! scaling to units
  scale_kl = scale_t * scale_displ_inv * 1.d9

  ! inner_core
  do ispec = 1, NSPEC_INNER_CORE
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

  ! kernel unit [ s / km^3 ]
  scale_kl = scale_t * scale_displ_inv * 1.d9

  ! scale the boundary kernels properly: *scale_kl gives s/km^3 and 1.d3 gives
  ! the relative boundary kernels (for every 1 km) in s/km^2
  moho_kl(:,:,:) = moho_kl(:,:,:) * scale_kl * 1.d3
  d400_kl(:,:,:) = d400_kl(:,:,:) * scale_kl * 1.d3
  d670_kl(:,:,:) = d670_kl(:,:,:) * scale_kl * 1.d3
  cmb_kl(:,:,:) = cmb_kl(:,:,:) * scale_kl * 1.d3
  icb_kl(:,:,:) = icb_kl(:,:,:) * scale_kl * 1.d3

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
  real(kind=CUSTOM_REAL),parameter :: scale_mass = RHOAV * (R_EARTH**3)
  integer :: irec_local
  character(len=MAX_STRING_LEN) :: outputname

  !scale_mass = RHOAV * (R_EARTH**3)

  ! computes derivatives
  do irec_local = 1, nrec_local
    ! rotate and scale the location derivatives to correspond to dn,de,dz
    sloc_der(:,irec_local) = matmul(transpose(nu_source(:,:,irec_local)),sloc_der(:,irec_local)) &
                             * scale_displ * scale_t

    ! rotate scale the moment derivatives to correspond to M[n,e,z][n,e,z]
    moment_der(:,:,irec_local) = matmul(matmul(transpose(nu_source(:,:,irec_local)),moment_der(:,:,irec_local)), &
               nu_source(:,:,irec_local)) * scale_t ** 3 / scale_mass

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
    stshift_der(irec_local) = stshift_der(irec_local) * scale_displ**2
    shdur_der(irec_local) = shdur_der(irec_local) * scale_displ**2
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
  scale_kl = scale_t * scale_displ_inv * 1.d9

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
  endif

  end subroutine save_kernels_Hessian

!
!-------------------------------------------------------------------------------------------------
!

  subroutine rotate_kernels_dble(cij_kl,cij_kl_spherical,theta_in,phi_in)

! Purpose : compute the kernels in r,theta,phi (cij_kl_spherical)
! from the kernels in x,y,z (cij_kl) (x,y,z to r,theta,phi)
! At r,theta,phi fixed
! theta and phi are in radians

! Coeff from Min's routine rotate_anisotropic_tensor
! with the help of Collect[Expand[cij],{dij}] in Mathematica

! Definition of the output array cij_kl_spherical :
! cij_kl_spherical(1) = C11 ; cij_kl_spherical(2) = C12 ; cij_kl_spherical(3) = C13
! cij_kl_spherical(4) = C14 ; cij_kl_spherical(5) = C15 ; cij_kl_spherical(6) = C16
! cij_kl_spherical(7) = C22 ; cij_kl_spherical(8) = C23 ; cij_kl_spherical(9) = C24
! cij_kl_spherical(10) = C25 ; cij_kl_spherical(11) = C26 ; cij_kl_spherical(12) = C33
! cij_kl_spherical(13) = C34 ; cij_kl_spherical(14) = C35 ; cij_kl_spherical(15) = C36
! cij_kl_spherical(16) = C44 ; cij_kl_spherical(17) = C45 ; cij_kl_spherical(18) = C46
! cij_kl_spherical(19) = C55 ; cij_kl_spherical(20) = C56 ; cij_kl_spherical(21) = C66
! where the Cij (Voigt's notation) are defined as function of
! the components of the elastic tensor in spherical coordinates
! by eq. (A.1) of Chen & Tromp, GJI 168 (2007)

  use constants

  implicit none

  real(kind=CUSTOM_REAL), dimension(21), intent(in) :: cij_kl
  real(kind=CUSTOM_REAL), dimension(21), intent(out) :: cij_kl_spherical

  real(kind=CUSTOM_REAL), intent(in) :: theta_in,phi_in

  ! local parameters
  double precision :: theta,phi
  double precision :: costheta,sintheta,cosphi,sinphi
  double precision :: costhetasq,sinthetasq,cosphisq,sinphisq
  double precision :: costwotheta,sintwotheta,costwophi,sintwophi
  double precision :: cosfourtheta,sinfourtheta,cosfourphi,sinfourphi
  double precision :: costhetafour,sinthetafour,cosphifour,sinphifour
  double precision :: sintwophisq,sintwothetasq
  double precision :: costhreetheta,sinthreetheta,costhreephi,sinthreephi

  theta = dble(theta_in)
  phi = dble(phi_in)

  costheta = dcos(theta)
  sintheta = dsin(theta)
  cosphi = dcos(phi)
  sinphi = dsin(phi)

  costhetasq = costheta * costheta
  sinthetasq = sintheta * sintheta
  cosphisq = cosphi * cosphi
  sinphisq = sinphi * sinphi

  costhetafour = costhetasq * costhetasq
  sinthetafour = sinthetasq * sinthetasq
  cosphifour = cosphisq * cosphisq
  sinphifour = sinphisq * sinphisq

  costwotheta = dcos(2.d0*theta)
  sintwotheta = dsin(2.d0*theta)
  costwophi = dcos(2.d0*phi)
  sintwophi = dsin(2.d0*phi)

  costhreetheta=dcos(3.d0*theta)
  sinthreetheta=dsin(3.d0*theta)
  costhreephi=dcos(3.d0*phi)
  sinthreephi=dsin(3.d0*phi)

  cosfourtheta = dcos(4.d0*theta)
  sinfourtheta = dsin(4.d0*theta)
  cosfourphi = dcos(4.d0*phi)
  sinfourphi = dsin(4.d0*phi)
  sintwothetasq = sintwotheta * sintwotheta
  sintwophisq = sintwophi * sintwophi

 cij_kl_spherical(1) = ONE_SIXTEENTH * (cij_kl(16) - cij_kl(16)* costwophi + &
     16.d0* cosphi*cosphisq* costhetafour* (cij_kl(1)* cosphi + cij_kl(6)* sinphi) + &
     2.d0* (cij_kl(15) + cij_kl(17))* sintwophi* sintwothetasq - &
     2.d0* (cij_kl(16)* cosfourtheta* sinphisq + &
     2.d0* costhetafour* (-4* cij_kl(7)* sinphifour - &
     (cij_kl(2) + cij_kl(21))* sintwophisq) + &
     8.d0* cij_kl(5)* cosphi*cosphisq* costheta*costhetasq* sintheta - &
     8.d0* cij_kl(8)* costhetasq* sinphisq* sinthetasq - &
     8.d0* cij_kl(12)* sinthetafour + &
     8.d0* cosphisq* costhetasq* sintheta* ((cij_kl(4) + &
     cij_kl(20))* costheta* sinphi - &
     (cij_kl(3) + cij_kl(19))*sintheta) + &
     8.d0* cosphi* costheta* (-cij_kl(11)* costheta*costhetasq* &
     sinphi*sinphisq + (cij_kl(10) + cij_kl(18))* costhetasq* sinphisq* sintheta + &
     cij_kl(14)* sintheta*sinthetasq) + 2.d0* sinphi* (cij_kl(13) + &
     cij_kl(9)* sinphisq)* sintwotheta + &
     sinphi* (-cij_kl(13) + cij_kl(9)* sinphisq)* sinfourtheta))

 cij_kl_spherical(2) = ONE_FOURTH * (costhetasq* (cij_kl(1) + 3.d0* cij_kl(2) + cij_kl(7) - &
      cij_kl(21) + (-cij_kl(1) + cij_kl(2) - cij_kl(7) + &
      cij_kl(21))* cosfourphi + (-cij_kl(6) + cij_kl(11))* sinfourphi) + &
      4.d0* (cij_kl(8)* cosphisq - cij_kl(15)* cosphi* sinphi + &
      cij_kl(3)* sinphisq)* sinthetasq - &
      2.d0* (cij_kl(10)* cosphisq*cosphi + &
      (cij_kl(9) - cij_kl(20))* cosphisq* sinphi + &
      (cij_kl(5) - cij_kl(18))* cosphi* sinphisq + &
      cij_kl(4)* sinphisq*sinphi)* sintwotheta)

 cij_kl_spherical(3) = ONE_EIGHTH * (sintwophi* (3.d0* cij_kl(15) - cij_kl(17) + &
     4.d0* (cij_kl(2) + cij_kl(21))* costhetasq* sintwophi* sinthetasq) + &
     4.d0* cij_kl(12)* sintwothetasq + 4.d0* cij_kl(1)* cosphifour* sintwothetasq + &
     2.d0* cosphi*cosphisq* (8.d0* cij_kl(6)* costhetasq* sinphi* sinthetasq + &
     cij_kl(5)* sinfourtheta) + 2.d0* cosphisq* (3.d0* cij_kl(3) -  cij_kl(19) + &
     (cij_kl(3) + cij_kl(19))* cosfourtheta + &
     (cij_kl(4) + cij_kl(20))* sinphi* sinfourtheta) + &
     2.d0* sinphi* (sinphi* (3.d0* cij_kl(8) - &
     cij_kl(16) + (cij_kl(8) + cij_kl(16))* cosfourtheta + &
     2.d0* cij_kl(7)* sinphisq* sintwothetasq)+ &
     (-cij_kl(13) + cij_kl(9)* sinphisq)* sinfourtheta)+ &
     2.d0* cosphi* ((cij_kl(15) + cij_kl(17))* cosfourtheta* sinphi + &
     8.d0* cij_kl(11)* costhetasq* sinphi*sinphisq* sinthetasq + &
     (-cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq)*sinfourtheta))

 cij_kl_spherical(4) = ONE_EIGHTH * (cosphi* costheta *(5.d0* cij_kl(4) - &
     cij_kl(9) + 4.d0* cij_kl(13) - &
     3.d0* cij_kl(20) + (cij_kl(4) + 3.d0* cij_kl(9) - &
     4.d0* cij_kl(13) + cij_kl(20))* costwotheta) + &
     ONE_HALF* (cij_kl(4) - cij_kl(9) + &
     cij_kl(20))* costhreephi * (costheta + 3.d0* costhreetheta) - &
     costheta* (-cij_kl(5) + 5.d0* cij_kl(10) + &
     4.d0* cij_kl(14) - 3.d0* cij_kl(18) + &
     (3.d0* cij_kl(5) + cij_kl(10) - &
     4.d0* cij_kl(14) + cij_kl(18))* costwotheta)* sinphi - &
     ONE_HALF* (cij_kl(5) - cij_kl(10) - cij_kl(18))* (costheta + &
     3.d0* costhreetheta)* sinthreephi + &
     4.d0* (cij_kl(6) - cij_kl(11))* cosfourphi* costhetasq* sintheta - &
     4.d0* (cij_kl(1) + cij_kl(3) - cij_kl(7) - cij_kl(8) + cij_kl(16) - cij_kl(19) + &
     (cij_kl(1) - cij_kl(3) - cij_kl(7) + cij_kl(8) + &
     cij_kl(16) - cij_kl(19))* costwotheta)* sintwophi* sintheta - &
     4.d0* (cij_kl(1) - cij_kl(2) + cij_kl(7) - &
     cij_kl(21))* costhetasq* sinfourphi* sintheta + &
     costwophi* ((cij_kl(6) + cij_kl(11) + 6.d0* cij_kl(15) - &
     2.d0* cij_kl(17))* sintheta + &
     (cij_kl(6) + cij_kl(11) - 2.d0* (cij_kl(15) + cij_kl(17)))* sinthreetheta))

 cij_kl_spherical(5) = ONE_FOURTH * (2.d0* (cij_kl(4) + &
     cij_kl(20))* cosphisq* (costwotheta + cosfourtheta)* sinphi + &
     2.d0* cij_kl(9)* (costwotheta + cosfourtheta)* sinphi*sinphisq + &
     16.d0* cij_kl(1)* cosphifour* costheta*costhetasq* sintheta + &
     4.d0* costheta*costhetasq* (-2.d0* cij_kl(8)* sinphisq + &
     4.d0* cij_kl(7)* sinphifour + &
     (cij_kl(2) + cij_kl(21))* sintwophisq)* sintheta + &
     4.d0* cij_kl(13)* (1.d0 + 2.d0* costwotheta)* sinphi* sinthetasq + &
     8.d0* costheta* (-2.d0* cij_kl(12) + cij_kl(8)* sinphisq)* sintheta*sinthetasq + &
     2.d0* cosphi*cosphisq* (cij_kl(5)* (costwotheta + cosfourtheta) + &
     8.d0* cij_kl(6)* costheta*costhetasq* sinphi* sintheta) + &
     2.d0* cosphi* (cosfourtheta* (-cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq) + &
     costwotheta* (cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq) + &
     8.d0* cij_kl(11)* costheta*costhetasq* sinphi*sinphisq* sintheta) - &
     (cij_kl(3) + cij_kl(16) + cij_kl(19) + &
     (cij_kl(3) - cij_kl(16) + cij_kl(19))* costwophi + &
     (cij_kl(15) + cij_kl(17))* sintwophi)* sinfourtheta)

 cij_kl_spherical(6) = ONE_HALF * costheta*costhetasq* ((cij_kl(6) + cij_kl(11))* costwophi + &
      (cij_kl(6) - cij_kl(11))* cosfourphi + 2.d0* (-cij_kl(1) + cij_kl(7))* sintwophi + &
      (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* sinfourphi) + &
      ONE_FOURTH* costhetasq* (-(cij_kl(4) + 3* cij_kl(9) + cij_kl(20))* cosphi - &
      3.d0* (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi + &
      (3.d0* cij_kl(5) + cij_kl(10) + cij_kl(18))* sinphi + &
      3.d0* (cij_kl(5) - cij_kl(10) - cij_kl(18))* sinthreephi)* sintheta + &
      costheta* ((cij_kl(15) + cij_kl(17))* costwophi + &
      (-cij_kl(3) + cij_kl(8) + cij_kl(16) - cij_kl(19))* sintwophi)* sinthetasq + &
      (-cij_kl(13)* cosphi + cij_kl(14)* sinphi)* sintheta*sinthetasq

 cij_kl_spherical(7) = cij_kl(7) * cosphifour - cij_kl(11)* cosphi*cosphisq* sinphi + &
      (cij_kl(2) + cij_kl(21))* cosphisq* sinphisq - &
      cij_kl(6)* cosphi* sinphi*sinphisq + &
      cij_kl(1)* sinphifour

 cij_kl_spherical(8) = ONE_HALF * (2.d0* costhetasq* sinphi* (-cij_kl(15)* cosphi + &
      cij_kl(3)* sinphi) + 2.d0* cij_kl(2)* cosphifour* sinthetasq + &
      (2.d0* cij_kl(2)* sinphifour + &
      (cij_kl(1) + cij_kl(7) - cij_kl(21))* sintwophisq)* sinthetasq + &
      cij_kl(4)* sinphi*sinphisq* sintwotheta + &
      cosphi*cosphisq* (2.d0* (-cij_kl(6) + cij_kl(11))* sinphi* sinthetasq + &
      cij_kl(10)* sintwotheta) + cosphi* sinphisq* (2.d0* (cij_kl(6) - &
      cij_kl(11))* sinphi* sinthetasq + &
      (cij_kl(5) - cij_kl(18))* sintwotheta) + &
      cosphisq* (2.d0* cij_kl(8)* costhetasq + &
      (cij_kl(9) - cij_kl(20))* sinphi* sintwotheta))

 cij_kl_spherical(9) = cij_kl(11) * cosphifour* sintheta - sinphi*sinphisq* (cij_kl(5)* costheta + &
      cij_kl(6)* sinphi* sintheta) +  cosphisq* sinphi* (-(cij_kl(10) + &
      cij_kl(18))* costheta + &
      3.d0* (cij_kl(6) - cij_kl(11))* sinphi* sintheta) + &
      cosphi* sinphisq* ((cij_kl(4) + cij_kl(20))* costheta + &
      2.d0* (-2.d0* cij_kl(1) + cij_kl(2) + cij_kl(21))* sinphi* sintheta) + &
      cosphi*cosphisq* (cij_kl(9)* costheta - 2.d0* (cij_kl(2) - 2.d0* cij_kl(7) + &
      cij_kl(21))* sinphi* sintheta)

 cij_kl_spherical(10) = ONE_FOURTH * (4.d0* costwotheta* (cij_kl(10)* cosphi*cosphisq + &
      (cij_kl(9) - cij_kl(20))* cosphisq* sinphi + &
      (cij_kl(5) - cij_kl(18))* cosphi* sinphisq + &
      cij_kl(4)* sinphi*sinphisq) + (cij_kl(1) + 3.d0* cij_kl(2) - &
      2.d0* cij_kl(3) + cij_kl(7) - &
      2.d0* cij_kl(8) - cij_kl(21) + 2.d0* (cij_kl(3) - cij_kl(8))* costwophi + &
      (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* cosfourphi + &
      2.d0* cij_kl(15)* sintwophi + &
      (-cij_kl(6) + cij_kl(11))* sinfourphi)* sintwotheta)

 cij_kl_spherical(11) = ONE_FOURTH * (2.d0* costheta* ((cij_kl(6) + cij_kl(11))* costwophi + &
      (-cij_kl(6) + cij_kl(11))* cosfourphi + &
      2.d0* (-cij_kl(1) + cij_kl(7))* sintwophi + &
      (cij_kl(1) - cij_kl(2) + cij_kl(7) - cij_kl(21))* sinfourphi) + &
      (-(cij_kl(4) + 3.d0* cij_kl(9) + cij_kl(20))* cosphi + &
      (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi + &
      (3.d0* cij_kl(5) + cij_kl(10) + cij_kl(18))* sinphi + &
      (-cij_kl(5) + cij_kl(10) + cij_kl(18))* sinthreephi)* sintheta)

 cij_kl_spherical(12) = ONE_SIXTEENTH * (cij_kl(16) - 2.d0* cij_kl(16)* cosfourtheta* sinphisq + &
      costwophi* (-cij_kl(16) + 8.d0* costheta* sinthetasq* ((cij_kl(3) - &
      cij_kl(8) + cij_kl(19))* costheta + &
      (cij_kl(5) - cij_kl(10) - cij_kl(18))* cosphi* sintheta)) + &
      2.d0* (cij_kl(15) + cij_kl(17))* sintwophi* sintwothetasq + &
      2.d0* (8.d0* cij_kl(12)* costhetafour + &
      8.d0* cij_kl(14)* cosphi* costheta*costhetasq* sintheta + &
      4.d0* cosphi* costheta* (cij_kl(5) + cij_kl(10) + cij_kl(18) + &
      (cij_kl(4) + cij_kl(20))* sintwophi)* &
      sintheta*sinthetasq + 8.d0* cij_kl(1)* cosphifour* sinthetafour + &
      8.d0* cij_kl(6)* cosphi*cosphisq* sinphi* sinthetafour + &
      8.d0* cij_kl(11)* cosphi* sinphi*sinphisq* sinthetafour + &
      8.d0* cij_kl(7)* sinphifour* sinthetafour + &
      2.d0* cij_kl(2)* sintwophisq* sinthetafour + &
      2.d0* cij_kl(21)* sintwophisq* sinthetafour + &
      2.d0* cij_kl(13)* sinphi* sintwotheta + &
      2.d0* cij_kl(9)* sinphi*sinphisq* sintwotheta + &
      cij_kl(3)* sintwothetasq + cij_kl(8)* sintwothetasq + &
      cij_kl(19)* sintwothetasq + cij_kl(13)* sinphi* sinfourtheta - &
      cij_kl(9)* sinphi*sinphisq* sinfourtheta))

 cij_kl_spherical(13) = ONE_EIGHTH * (cosphi* costheta* (cij_kl(4) + 3.d0* cij_kl(9) + &
      4.d0* cij_kl(13) + cij_kl(20) - (cij_kl(4) + 3.d0* cij_kl(9) - &
      4.d0* cij_kl(13) + cij_kl(20))* costwotheta) + 4.d0* (-cij_kl(1) - &
      cij_kl(3) + cij_kl(7) + cij_kl(8) + cij_kl(16) - cij_kl(19) + &
      (cij_kl(1) - cij_kl(3) - cij_kl(7) + cij_kl(8) + cij_kl(16) - &
      cij_kl(19))* costwotheta)* sintwophi* sintheta + &
      4.d0* (cij_kl(6) - cij_kl(11))* cosfourphi* sinthetasq*sintheta - &
      4.d0* (cij_kl(1) - cij_kl(2) + cij_kl(7) - &
      cij_kl(21))* sinfourphi* sinthetasq*sintheta + &
      costheta* ((-3.d0* cij_kl(5) - cij_kl(10) - 4.d0* cij_kl(14) - &
      cij_kl(18) + (3.d0* cij_kl(5) + cij_kl(10) - 4.d0* cij_kl(14) + &
      cij_kl(18))* costwotheta)* sinphi + 6.d0* ((cij_kl(4) - cij_kl(9) + &
      cij_kl(20))* costhreephi + (-cij_kl(5) + cij_kl(10) + &
      cij_kl(18))* sinthreephi)* sinthetasq) + costwophi* ((3* cij_kl(6) + &
      3.d0* cij_kl(11) + 2.d0* (cij_kl(15) + cij_kl(17)))* sintheta - &
      (cij_kl(6) + cij_kl(11) - 2.d0* (cij_kl(15) + &
      cij_kl(17)))* sinthreetheta))

 cij_kl_spherical(14) = ONE_FOURTH * (2.d0* cij_kl(13)* (costwotheta + cosfourtheta)* sinphi + &
      8.d0* costheta*costhetasq* (-2.d0* cij_kl(12) + cij_kl(8)* sinphisq)* sintheta + &
      4.d0* (cij_kl(4) + cij_kl(20))* cosphisq* (1.d0 + &
      2.d0* costwotheta)* sinphi* sinthetasq + &
      4.d0* cij_kl(9)* (1.d0 + 2.d0* costwotheta)* sinphi*sinphisq* sinthetasq + &
      16.d0* cij_kl(1)* cosphifour* costheta* sintheta*sinthetasq + &
      4.d0* costheta* (-2.d0* cij_kl(8)* sinphisq + 4.d0* cij_kl(7)* sinphifour + &
      (cij_kl(2) + cij_kl(21))* sintwophisq)* sintheta*sinthetasq + &
      4.d0* cosphi*cosphisq* sinthetasq* (cij_kl(5) + 2.d0* cij_kl(5)* costwotheta + &
      4.d0* cij_kl(6)* costheta* sinphi* sintheta) + &
      2.d0* cosphi* (cosfourtheta* (cij_kl(14) - (cij_kl(10) + cij_kl(18))* sinphisq) + &
      costwotheta* (cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq) + &
      8.d0* cij_kl(11)* costheta* sinphi*sinphisq* sintheta*sinthetasq) + &
      (cij_kl(3) + cij_kl(16) + cij_kl(19) + (cij_kl(3) - cij_kl(16) + &
      cij_kl(19))* costwophi + (cij_kl(15) + cij_kl(17))* sintwophi)* sinfourtheta)

 cij_kl_spherical(15) = costwophi * costheta* (-cij_kl(17) + (cij_kl(15) + cij_kl(17))* costhetasq) + &
       ONE_SIXTEENTH* (-((11.d0* cij_kl(4) + cij_kl(9) + 4.d0* cij_kl(13) - &
       5.d0* cij_kl(20))* cosphi + (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi - &
       (cij_kl(5) + 11.d0* cij_kl(10) + 4.d0* cij_kl(14) - &
       5.d0* cij_kl(18))* sinphi + (-cij_kl(5) + cij_kl(10) + &
       cij_kl(18))* sinthreephi)* sintheta + &
       8.d0* costheta* ((-cij_kl(1) - cij_kl(3) + cij_kl(7) + cij_kl(8) - cij_kl(16) +&
       cij_kl(19) + (cij_kl(1) - cij_kl(3) - &
       cij_kl(7) + cij_kl(8) + cij_kl(16) - cij_kl(19))* costwotheta)* sintwophi +&
       ((cij_kl(6) + cij_kl(11))* costwophi + &
       (cij_kl(6) - cij_kl(11))* cosfourphi + (-cij_kl(1) + cij_kl(2) - cij_kl(7) +&
       cij_kl(21))* sinfourphi)* sinthetasq) +&
       ((cij_kl(4) + 3.d0* cij_kl(9) - 4.d0* cij_kl(13) + cij_kl(20))* cosphi + &
       3.d0* (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi - &
       (3.d0* cij_kl(5) + cij_kl(10) - 4.d0* cij_kl(14) + cij_kl(18))* sinphi + &
       3.d0* (-cij_kl(5) + cij_kl(10) + cij_kl(18))* sinthreephi)* sinthreetheta)

 cij_kl_spherical(16) = ONE_FOURTH *(cij_kl(1) - cij_kl(2) + cij_kl(7) + cij_kl(16) + &
       cij_kl(19) + cij_kl(21) + 2.d0*(cij_kl(16) - cij_kl(19))*costwophi* costhetasq + &
       (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(16) + &
       cij_kl(19) - cij_kl(21))*costwotheta - 2.d0* cij_kl(17)* costhetasq* sintwophi + &
       2.d0* ((-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* cosfourphi + &
       (-cij_kl(6) + cij_kl(11))* sinfourphi)* sinthetasq + ((cij_kl(5) - cij_kl(10) +&
       cij_kl(18))* cosphi + (-cij_kl(5) + cij_kl(10) + cij_kl(18))* costhreephi +&
       (-cij_kl(4) + cij_kl(9) + cij_kl(20))* sinphi - &
       (cij_kl(4) - cij_kl(9) + cij_kl(20))* sinthreephi)* sintwotheta)

 cij_kl_spherical(17) = ONE_EIGHTH * (4.d0* costwophi* costheta* (cij_kl(6) + cij_kl(11) - &
       2.d0* cij_kl(15) - (cij_kl(6) + cij_kl(11) - 2.d0* (cij_kl(15) + &
       cij_kl(17)))* costwotheta) - (2.d0* cosphi* (-3.d0* cij_kl(4) +&
       cij_kl(9) + 2.d0* cij_kl(13) + cij_kl(20) + (cij_kl(4) - cij_kl(9) + &
       cij_kl(20))* costwophi) - (cij_kl(5) - 5.d0* cij_kl(10) + &
       4.d0* cij_kl(14) + 3.d0* cij_kl(18))* sinphi + (-cij_kl(5) + cij_kl(10) + &
       cij_kl(18))* sinthreephi)* sintheta + &
       8.d0* costheta* ((-cij_kl(1) + cij_kl(3) + cij_kl(7) - cij_kl(8) + &
       (cij_kl(1) - cij_kl(3) - cij_kl(7) + cij_kl(8) + cij_kl(16) - &
       cij_kl(19))* costwotheta)* sintwophi + ((cij_kl(6) - cij_kl(11))* cosfourphi + &
       (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* sinfourphi)* sinthetasq) +&
       ((cij_kl(4) + 3.d0* cij_kl(9) - 4.d0* cij_kl(13) + cij_kl(20))* cosphi + &
       3.d0* (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi - &
       (3.d0* cij_kl(5) + cij_kl(10) - 4.d0* cij_kl(14) + cij_kl(18))* sinphi + &
       3.d0* (-cij_kl(5) + cij_kl(10) + cij_kl(18))* sinthreephi)* sinthreetheta)

 cij_kl_spherical(18) = ONE_HALF * ((cij_kl(5) - cij_kl(10) + cij_kl(18))* cosphi* costwotheta - &
       (cij_kl(5) - cij_kl(10) - cij_kl(18))* costhreephi* costwotheta - &
       2.d0* (cij_kl(4) - cij_kl(9) + &
       (cij_kl(4) - cij_kl(9) + cij_kl(20))* costwophi)* costwotheta* sinphi + &
       (cij_kl(1) - cij_kl(2) + cij_kl(7) - cij_kl(16) - cij_kl(19) + cij_kl(21) + &
       (-cij_kl(16) + cij_kl(19))* costwophi + &
       (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* cosfourphi + &
       cij_kl(17)* sintwophi + &
       (-cij_kl(6) + cij_kl(11))* sinfourphi)* sintwotheta)

 cij_kl_spherical(19) = ONE_FOURTH * (cij_kl(16) - cij_kl(16)* costwophi + &
      (-cij_kl(15) + cij_kl(17))* sintwophi + &
      4.d0* cij_kl(12)* sintwothetasq + &
      2.d0* (2.d0* cij_kl(1)* cosphifour* sintwothetasq + &
      cosphi*cosphisq* (8.d0* cij_kl(6)* costhetasq* sinphi* sinthetasq + &
      cij_kl(5)* sinfourtheta) + cosphisq* (-cij_kl(3) + cij_kl(19) + (cij_kl(3) +&
      cij_kl(19))* cosfourtheta + (cij_kl(4) + cij_kl(20))* sinphi* sinfourtheta) + &
      sinphi* (cosfourtheta* ((cij_kl(15) + cij_kl(17))* cosphi + &
      cij_kl(16)* sinphi) + (cij_kl(2) + cij_kl(7) - 2.d0* cij_kl(8) + cij_kl(21) + &
      (cij_kl(2) - cij_kl(7) + cij_kl(21))* costwophi)* sinphi* sintwothetasq + &
      (-cij_kl(13) + cij_kl(9)* sinphisq)* sinfourtheta) + &
      cosphi* (8.d0* cij_kl(11)* costhetasq* sinphi*sinphisq* sinthetasq + &
      (-cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq)* sinfourtheta)))

 cij_kl_spherical(20) = ONE_EIGHTH * (2.d0* cosphi* costheta* (-3.d0* cij_kl(4) - cij_kl(9) + &
      4.d0* cij_kl(13) + cij_kl(20) + (cij_kl(4) + 3.d0* cij_kl(9) - &
      4.d0* cij_kl(13) + cij_kl(20))* costwotheta) + &
      (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi* (costheta + &
      3.d0* costhreetheta) - &
      2.d0* costheta* (-cij_kl(5) - 3.d0* cij_kl(10) + 4.d0* cij_kl(14) + &
      cij_kl(18) + (3.d0* cij_kl(5) + &
      cij_kl(10) - 4.d0* cij_kl(14) + cij_kl(18))*costwotheta)* sinphi - &
      (cij_kl(5) - cij_kl(10) - cij_kl(18))* &
      (costheta + 3.d0* costhreetheta)* sinthreephi + 8.d0* (cij_kl(6) - &
      cij_kl(11))* cosfourphi* costhetasq* sintheta - 8.d0* (cij_kl(1) - &
      cij_kl(3) - cij_kl(7) + cij_kl(8) + &
      (cij_kl(1) - cij_kl(3) - cij_kl(7) + cij_kl(8) + cij_kl(16) - &
      cij_kl(19))* costwotheta)* sintwophi* sintheta - &
      8.d0* (cij_kl(1) - cij_kl(2) + cij_kl(7) - &
      cij_kl(21))* costhetasq* sinfourphi* sintheta + &
      2.d0* costwophi* ((cij_kl(6) + cij_kl(11) - 2.d0* cij_kl(15) + &
      2.d0* cij_kl(17))* sintheta + &
      (cij_kl(6) + cij_kl(11) - 2.d0* (cij_kl(15) + cij_kl(17)))* sinthreetheta))

 cij_kl_spherical(21) = ONE_FOURTH * (cij_kl(1) - cij_kl(2) + cij_kl(7) + cij_kl(16) + &
      cij_kl(19) + cij_kl(21) - 2.d0* (cij_kl(1) - cij_kl(2) + cij_kl(7) - &
      cij_kl(21))* cosfourphi* costhetasq + &
      (cij_kl(1) - cij_kl(2) + cij_kl(7) - cij_kl(16) - cij_kl(19) + &
      cij_kl(21))* costwotheta + &
      2.d0* (-cij_kl(6) + cij_kl(11))* costhetasq* sinfourphi - &
      2.d0* ((-cij_kl(16) + cij_kl(19))* costwophi + cij_kl(17)* sintwophi)* sinthetasq - &
      ((cij_kl(5) - cij_kl(10) + cij_kl(18))* cosphi + (-cij_kl(5) + cij_kl(10) +&
      cij_kl(18))* costhreephi + &
      (-cij_kl(4) + cij_kl(9) + cij_kl(20))* sinphi - (cij_kl(4) - cij_kl(9) + &
      cij_kl(20))* sinthreephi)* sintwotheta)

  end subroutine rotate_kernels_dble

