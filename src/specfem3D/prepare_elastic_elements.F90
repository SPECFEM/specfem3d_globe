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


  subroutine prepare_elastic_elements()

! precomputes effective elastic tensor for aniso, tiso and iso elements.
!
! this routine also computes the unrelaxed elastic parameters in case of attenuation.
! this can only be done after the attenuation factors have been determined.
!
! for example:
!   see Komatitsch & Tromp, 1999, eq. (7)
!         C_u = C_r [ 1 - sum( 1 - tau_strain / tau_stress ) ]
!       with C_u unrelaxed and C_r relaxed moduli
!
!   or see Komatitsch & Tromp, 2002, eq. (10)
!         Mu_u = Mu_r * [ 1 - sum( 1 - tau_strain / tau_stress) ]
!       with Mu_u unrelaxed and Mu_r relaxed shear moduli

! note: muvstore is still used in sensitivity kernels computation and needed in the relaxed form, i.e.,
!       c11... stores are not used any further after the time loop.
!       we will restore the relaxed elastic moduli when needed in the save_kernels() routine.

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: one_minus_sum_beta_use,minus_sum_beta

  ! tiso elements
  real(kind=CUSTOM_REAL) :: rhovphsq,rhovpvsq,rhovshsq,rhovsvsq,eta_aniso,phi,theta

  !real(kind=CUSTOM_REAL) :: rhovphsq,sinphifour,cosphisq,sinphisq,costhetasq,rhovsvsq,sinthetasq, &
  !      cosphifour,costhetafour,rhovpvsq,sinthetafour,rhovshsq,cosfourphi, &
  !      costwotheta,cosfourtheta,sintwophisq,costheta,sinphi,sintheta,cosphi, &
  !      sintwotheta,costwophi,sintwophi,costwothetasq,costwophisq,phi,theta

  ! aniso element
  real(kind=CUSTOM_REAL) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                            c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

!  real(kind=CUSTOM_REAL) :: two_rhovsvsq,two_rhovshsq ! two_rhovpvsq,two_rhovphsq
!  real(kind=CUSTOM_REAL) :: four_rhovsvsq,four_rhovshsq ! four_rhovpvsq,four_rhovphsq

!  real(kind=CUSTOM_REAL) :: twoetaminone,etaminone,eta_aniso
!  real(kind=CUSTOM_REAL) :: two_eta_aniso,four_eta_aniso,six_eta_aniso

!  real(kind=CUSTOM_REAL) :: templ1,templ1_cos,templ2,templ2_cos,templ3,templ3_two,templ3_cos
  real(kind=CUSTOM_REAL) :: kappavl,kappahl,muvl,muhl,mul,A,F,L,N

  integer :: ispec,iglob,i_SLS,icount,icount_iso

#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
#else
  integer :: i,j,k
#endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing elastic element arrays"
    call flush_IMAIN()
  endif

  ! checks if anything to do
  ! GPU kernels still use original arrays
  if (GPU_MODE) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  GPU mode with original arrays"
      call flush_IMAIN()
    endif
    ! done
    return
  endif

  ! user output
  if (ATTENUATION_VAL) then
    if (myrank == 0) then
      write(IMAIN,*) "  using attenuation: shifting to unrelaxed moduli"
      call flush_IMAIN()
    endif
  endif

  ! crust/mantle
  ! pre-computes c11.. factors
  if (ANISOTROPIC_3D_MANTLE_VAL) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  crust/mantle anisotropic elements"
      call flush_IMAIN()
    endif

    ! only shifting for attenuation case
    if (ATTENUATION_VAL) then
      icount = 0

! openmp solver
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ispec, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk, &
#else
!$OMP i,j,k, &
#endif
!$OMP minus_sum_beta,mul,muvl,muhl,A,L,F,eta_aniso,c11,c12,c13,c22,c23,c33,c44,c55,c66)
!$OMP DO
      do ispec = 1,NSPEC_CRUST_MANTLE
        ! all elements are fully anisotropic
        DO_LOOP_IJK
          ! Sieminski, 2007:
          ! A = 1/8 (3 C11 + 3 C22 + 2 C12 + 4 C66)
          ! C = C33
          ! N = 1/8 (C11 + C22 - 2 C12 + 4 C66)
          ! L = 1/2 (C44 + C55)
          ! F = 1/2 (C13 + C23)
          !
          ! Anderson & Dziewonski, 1982: "Upper mantle anisotropy: evidence from free oscillations", GJR
          ! A = rho * vph**2
          ! C = rho * vpv**2
          ! N = rho * vsh**2
          ! L = rho * vsv**2
          ! F = eta * (A - 2*L)
          !
          ! and therefore (assuming radial axis symmetry)
          ! C11 = A = rho * vph**2
          ! C33 = C = rho * vpv**2
          ! C44 = L = rho * vsv**2
          ! C13 = F = eta * (A - 2*L)
          ! C12 = C11 - 2 C66 = A - 2*N = rho * (vph**2 - 2 * vsh**2)
          ! C22 = C11
          ! C23 = C13
          ! C55 = C44
          ! C66 = N = rho * vsh**2 = (C11-C12)/2
          !
          ! scaling:
          ! muvl = muvl * one_minus_sum_beta_use = muvl*( 1 - sum_beta) = muvl - muvl*sum_beta = muvl + muvl*minus_sum_beta
          ! muhl = muhl * one_minus_sum_beta_use
          !
          ! rho vsv**2 = muv -> L = L * (1 - sum_beta) = L - L * sum_beta = L + muv * minus_sum_beta
          ! rho vsh**2 = muh -> N = N + muh * minus_sum_beta
          !
          ! rho vph**2 = A = kappah + 4/3 muh -> kappah + 4/3 muh * (1 - sum_beta) = kappah + 4/3 muh - 4/3 muh*sum_beta
          ! rho vpv**2 = C = kappav + 4/3 muv -> kappav + 4/3 muv - 4/3 muv*sum_beta
          !
          ! c11' = A' = c11 - 4/3 muh*sum_beta
          ! c33' = C' = c33 - 4/3 muv*sum_beta
          ! c66' = N' = c66 - muh*sum_beta
          ! c44' = L' = c44 - muv*sum_beta
          ! c55' = c55 - muv*sum_beta     ! to be consistent with tiso case c55 == c44
          ! c22' = c22 - 4/3 muh*sum_beta ! to be consistent with tiso case c22 == c11
          !
          ! c13' = F' = eta*(A' - 2*L') = eta*(A - 4/3 muh*sum_beta - 2*(L - muv*sum_beta))
          !                             = c13 + eta*( 2*muv*sum_beta - 4/3 muh*sum_beta)
          ! c23' = c23 + eta*( 2*muv*sum_beta - 4/3 muh*sum_beta) ! be consistent with tiso case c23 == c13
          ! c12' = A' - 2*N' = A - 4/3 muh*sum_beta - 2*(N - muh*sum_beta) = A - 2*N - 4/3 muh*sum_beta + 2 muh*sum_beta
          !                                                                = c12 + 2/3 muh*sum_beta

          c11 = c11store_crust_mantle(INDEX_IJK,ispec)
          c12 = c12store_crust_mantle(INDEX_IJK,ispec)
          c13 = c13store_crust_mantle(INDEX_IJK,ispec)
          c22 = c22store_crust_mantle(INDEX_IJK,ispec)
          c23 = c23store_crust_mantle(INDEX_IJK,ispec)
          c33 = c33store_crust_mantle(INDEX_IJK,ispec)
          c44 = c44store_crust_mantle(INDEX_IJK,ispec)
          c55 = c55store_crust_mantle(INDEX_IJK,ispec)
          c66 = c66store_crust_mantle(INDEX_IJK,ispec)

          ! precompute terms for attenuation if needed
          if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
            minus_sum_beta =  one_minus_sum_beta_crust_mantle(INDEX_IJK,ispec) - 1.0_CUSTOM_REAL
          else
            minus_sum_beta =  one_minus_sum_beta_crust_mantle(1,1,1,ispec) - 1.0_CUSTOM_REAL
          endif

!daniel todo
          if (.true.) then
            ! this is the original routine to shift moduli, using only muv (from c44) to scale
            mul = c44 * minus_sum_beta
            c11 = c11 + FOUR_THIRDS * mul ! * minus_sum_beta * mul
            c12 = c12 - TWO_THIRDS * mul
            c13 = c13 - TWO_THIRDS * mul
            c22 = c22 + FOUR_THIRDS * mul
            c23 = c23 - TWO_THIRDS * mul
            c33 = c33 + FOUR_THIRDS * mul
            c44 = c44 + mul
            c55 = c55 + mul
            c66 = c66 + mul
          else
            ! new: tries to shift moduli by separating muv and muh factors.
            !      still needs rotations to rotate back and forth from SPECFEM global axis to a radial symmetry axis
            !      since this shift assumes a radial symmetry
            A = 0.125d0 * (3.d0 * c11 + 3.d0 * c22 + 2.d0 * c12 + 4.d0 * c66)
            N = 0.125d0 * (c11 + c22 - 2.d0 * c12 + 4.d0 * c66)
            L = 0.5d0 * (c44 + c55)
            F = 0.5d0 * (c13 + c23)
            eta_aniso = F / (A - 2.d0*L)   ! eta = F / (A-2L)

            muvl = L * minus_sum_beta     ! c44 -> L -> muv
            muhl = N * minus_sum_beta     ! c66 -> N -> muh

            c11 = c11 + FOUR_THIRDS * muhl ! * minus_sum_beta * mul
            c12 = c12 - TWO_THIRDS * muhl
            c13 = c13 + eta_aniso * (FOUR_THIRDS * muhl - 2.d0*muvl)
            c22 = c22 + FOUR_THIRDS * muhl
            c23 = c23 + eta_aniso * (FOUR_THIRDS * muhl - 2.d0*muvl)
            c33 = c33 + FOUR_THIRDS * muvl
            c44 = c44 + muvl
            c55 = c55 + muvl
            c66 = c66 + muhl
          endif

          ! stores unrelaxed factors
          c11store_crust_mantle(INDEX_IJK,ispec) = c11
          c12store_crust_mantle(INDEX_IJK,ispec) = c12
          c13store_crust_mantle(INDEX_IJK,ispec) = c13
          c22store_crust_mantle(INDEX_IJK,ispec) = c22
          c23store_crust_mantle(INDEX_IJK,ispec) = c23
          c33store_crust_mantle(INDEX_IJK,ispec) = c33
          c44store_crust_mantle(INDEX_IJK,ispec) = c44
          c55store_crust_mantle(INDEX_IJK,ispec) = c55
          c66store_crust_mantle(INDEX_IJK,ispec) = c66
        ENDDO_LOOP_IJK

        ! counter
!$OMP ATOMIC
        icount = icount + 1
      enddo
!$OMP ENDDO
!$OMP END PARALLEL

      ! user output
      if (myrank == 0) then
        write(IMAIN,*) "  aniso elements = ",icount
        call flush_IMAIN()
      endif
      ! check
      if (icount > NSPECMAX_ANISO_MANTLE) stop 'Error invalid number of aniso elements in prepare_timerun_aniso'
    endif ! ATTENUATION_VAL

  else
    ! iso/tiso mantle elements
    if (myrank == 0) then
      write(IMAIN,*) "  crust/mantle transverse isotropic and isotropic elements"
      call flush_IMAIN()
    endif

    ! crust/mantle
    icount = 0
    icount_iso = 0

! openmp solver
!$OMP PARALLEL if (NSPEC_CRUST_MANTLE > 10) &
!$OMP DEFAULT(PRIVATE) &
!$OMP SHARED(ispec_is_tiso_crust_mantle, kappavstore_crust_mantle,muvstore_crust_mantle, &
!$OMP kappahstore_crust_mantle,muhstore_crust_mantle, &
!$OMP one_minus_sum_beta_crust_mantle,eta_anisostore_crust_mantle,ibool_crust_mantle,rstore_crust_mantle, &
!$OMP c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
!$OMP c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle,c23store_crust_mantle, &
!$OMP c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
!$OMP c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle,c44store_crust_mantle, &
!$OMP c45store_crust_mantle,c46store_crust_mantle,c55store_crust_mantle,c56store_crust_mantle, &
!$OMP c66store_crust_mantle)
!$OMP DO
    do ispec = 1,NSPEC_CRUST_MANTLE

      if (ispec_is_tiso_crust_mantle(ispec)) then
        ! tiso elements
        ! precomputes stress factors for transversely isotropic elements
        DO_LOOP_IJK
          ! note: the mesh is built such that anisotropic elements are created first in anisotropic layers,
          !           thus they are listed first ( see in create_regions_mesh.f90: perm_layer() ordering )
          !           this is therefore still in bounds of 1:NSPECMAX_TISO_MANTLE even if NSPECMAX_TISO is less than NSPEC

          ! use kappa and mu from transversely isotropic model
          kappavl = kappavstore_crust_mantle(INDEX_IJK,ispec)
          muvl = muvstore_crust_mantle(INDEX_IJK,ispec)

          kappahl = kappahstore_crust_mantle(INDEX_IJK,ispec)
          muhl = muhstore_crust_mantle(INDEX_IJK,ispec)

          ! use unrelaxed parameters if attenuation
          ! eta does not need to be shifted since it is a ratio
          if (ATTENUATION_VAL) then
            ! precompute terms for attenuation if needed
            if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
              one_minus_sum_beta_use = one_minus_sum_beta_crust_mantle(INDEX_IJK,ispec)
            else
              one_minus_sum_beta_use = one_minus_sum_beta_crust_mantle(1,1,1,ispec)
            endif
            muvl = muvl * one_minus_sum_beta_use
            muhl = muhl * one_minus_sum_beta_use
          endif

          ! stores unrelaxed shear moduli
          muvstore_crust_mantle(INDEX_IJK,ispec) = muvl

          rhovpvsq = kappavl + FOUR_THIRDS * muvl  !!! that is C
          rhovphsq = kappahl + FOUR_THIRDS * muhl  !!! that is A

          rhovsvsq = muvl  !!! that is L
          rhovshsq = muhl  !!! that is N

          eta_aniso = eta_anisostore_crust_mantle(INDEX_IJK,ispec)  !!! that is  F / (A - 2 L)

          ! use mesh coordinates to get theta and phi
          ! rstore contains theta and phi
          iglob = ibool_crust_mantle(INDEX_IJK,ispec)

          theta = rstore_crust_mantle(2,iglob)
          phi = rstore_crust_mantle(3,iglob)

          ! rotates radial to global reference
          call rotate_tensor_tiso_to_cij(theta,phi, &
                                         rhovphsq,rhovpvsq,rhovsvsq,rhovshsq,eta_aniso, &
                                         c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                         c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

          c11store_crust_mantle(INDEX_IJK,ispec) = c11
          c12store_crust_mantle(INDEX_IJK,ispec) = c12
          c13store_crust_mantle(INDEX_IJK,ispec) = c13
          c14store_crust_mantle(INDEX_IJK,ispec) = c14
          c15store_crust_mantle(INDEX_IJK,ispec) = c15
          c16store_crust_mantle(INDEX_IJK,ispec) = c16
          c22store_crust_mantle(INDEX_IJK,ispec) = c22
          c23store_crust_mantle(INDEX_IJK,ispec) = c23
          c24store_crust_mantle(INDEX_IJK,ispec) = c24
          c25store_crust_mantle(INDEX_IJK,ispec) = c25
          c26store_crust_mantle(INDEX_IJK,ispec) = c26
          c33store_crust_mantle(INDEX_IJK,ispec) = c33
          c34store_crust_mantle(INDEX_IJK,ispec) = c34
          c35store_crust_mantle(INDEX_IJK,ispec) = c35
          c36store_crust_mantle(INDEX_IJK,ispec) = c36
          c44store_crust_mantle(INDEX_IJK,ispec) = c44
          c45store_crust_mantle(INDEX_IJK,ispec) = c45
          c46store_crust_mantle(INDEX_IJK,ispec) = c46
          c55store_crust_mantle(INDEX_IJK,ispec) = c55
          c56store_crust_mantle(INDEX_IJK,ispec) = c56
          c66store_crust_mantle(INDEX_IJK,ispec) = c66
        ENDDO_LOOP_IJK

        ! counter
!$OMP ATOMIC
        icount = icount + 1
      else
        ! isotropic elements
        ! shift to unrelaxed parameters if attenuation for stiffness computations in time loop
        if (ATTENUATION_VAL) then
          DO_LOOP_IJK
            ! layer with no transverse isotropy, use kappav and muv
            mul = muvstore_crust_mantle(INDEX_IJK,ispec)

            ! precompute terms for attenuation if needed
            if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
              one_minus_sum_beta_use = one_minus_sum_beta_crust_mantle(INDEX_IJK,ispec)
            else
              one_minus_sum_beta_use = one_minus_sum_beta_crust_mantle(1,1,1,ispec)
            endif
            mul = mul * one_minus_sum_beta_use

            ! stores unrelaxed shear moduli
            muvstore_crust_mantle(INDEX_IJK,ispec) = mul
          ENDDO_LOOP_IJK
        endif
        ! counter
!$OMP ATOMIC
        icount_iso = icount_iso + 1
      endif ! ispec_is_tiso_crust_mantle(ispec)
    enddo
!$OMP ENDDO
!$OMP END PARALLEL

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  tiso elements = ",icount
      write(IMAIN,*) "  iso elements  = ",icount_iso
      call flush_IMAIN()
    endif

    if (icount > NSPECMAX_TISO_MANTLE) stop 'Error invalid number of tiso elements in prepare_timerun_aniso'
    if (icount_iso > NSPECMAX_ISO_MANTLE) stop 'Error invalid number of iso elements in prepare_timerun_aniso'

    ! since we scale muv and c11,.. stores we must divide with this factor to use the relaxed moduli for the modulus defect
    ! calculation in updating the memory variables
    if (ATTENUATION_VAL) then

! openmp solver
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ispec,i_SLS, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk, &
#else
!$OMP i,j,k, &
#endif
!$OMP one_minus_sum_beta_use)
!$OMP DO
      do ispec = 1,NSPEC_CRUST_MANTLE
        DO_LOOP_IJK
          ! gets factor
          if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
            one_minus_sum_beta_use = one_minus_sum_beta_crust_mantle(INDEX_IJK,ispec)
          else
            one_minus_sum_beta_use = one_minus_sum_beta_crust_mantle(1,1,1,ispec)
          endif

          ! corrects factor_common to obtain relaxed moduli in moduli defect
          do i_SLS = 1,N_SLS
            if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
              factor_common_crust_mantle(INDEX_IJK,i_SLS,ispec) = &
                  factor_common_crust_mantle(INDEX_IJK,i_SLS,ispec) / one_minus_sum_beta_use
            else
              factor_common_crust_mantle(1,1,1,i_SLS,ispec) = &
                  factor_common_crust_mantle(1,1,1,i_SLS,ispec) / one_minus_sum_beta_use
            endif
          enddo
        ENDDO_LOOP_IJK
      enddo
!$OMP ENDDO
!$OMP END PARALLEL
    endif ! ATTENUATION_VAL
  endif

  ! inner core
  if (ANISOTROPIC_INNER_CORE_VAL) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  inner core anisotropic elements"
      call flush_IMAIN()
    endif

    if (ATTENUATION_VAL) then
      ! only scales for attenuation

      ! anisotropic inner core elements
      icount = 0

! openmp solver
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ispec, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk, &
#else
!$OMP i,j,k, &
#endif
!$OMP minus_sum_beta,mul,c11,c12,c13,c33,c44)
!$OMP DO
      do ispec = 1,NSPEC_INNER_CORE

        ! exclude fictitious elements in central cube
        if (idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

        DO_LOOP_IJK
          ! elastic tensor for hexagonal symmetry in reduced notation:
          !
          !      c11 c12 c13  0   0        0
          !      c12 c11 c13  0   0        0
          !      c13 c13 c33  0   0        0
          !       0   0   0  c44  0        0
          !       0   0   0   0  c44       0
          !       0   0   0   0   0  (c11-c12)/2
          !
          !       in terms of the A, C, L, N and F of Love (1927):
          !
          !       c11 = A
          !       c12 = A-2N
          !       c13 = F
          !       c33 = C
          !       c44 = L
          c11 = c11store_inner_core(INDEX_IJK,ispec)
          c12 = c12store_inner_core(INDEX_IJK,ispec)
          c13 = c13store_inner_core(INDEX_IJK,ispec)
          c33 = c33store_inner_core(INDEX_IJK,ispec)
          c44 = c44store_inner_core(INDEX_IJK,ispec)

          ! precompute terms for attenuation if needed
          if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
            minus_sum_beta =  one_minus_sum_beta_inner_core(INDEX_IJK,ispec) - 1.0_CUSTOM_REAL
          else
            minus_sum_beta =  one_minus_sum_beta_inner_core(1,1,1,ispec) - 1.0_CUSTOM_REAL
          endif

          ! note: muvstore has still original values, whereas c11store,.., have been rescaled due to attenuation
          !
          ! please check:
          ! shear along [100] direction: mul = c44
          ! instead of
          !mul = muvstore(INDEX_IJK,ispec) * minus_sum_beta
          ! this would be following the implementation from above for fully anisotropic elements...
          mul = c44 * minus_sum_beta

          c11 = c11 + FOUR_THIRDS * mul ! * minus_sum_beta * mul
          c12 = c12 - TWO_THIRDS * mul
          c13 = c13 - TWO_THIRDS * mul
          c33 = c33 + FOUR_THIRDS * mul
          c44 = c44 + mul

          ! stores scaled values
          c11store_inner_core(INDEX_IJK,ispec) = c11
          c12store_inner_core(INDEX_IJK,ispec) = c12
          c13store_inner_core(INDEX_IJK,ispec) = c13
          c33store_inner_core(INDEX_IJK,ispec) = c33
          c44store_inner_core(INDEX_IJK,ispec) = c44
        ENDDO_LOOP_IJK
        ! counter
!$OMP ATOMIC
        icount = icount + 1
      enddo
!$OMP ENDDO
!$OMP END PARALLEL

      ! user output
      if (myrank == 0) then
        write(IMAIN,*) "  aniso elements = ",icount
        call flush_IMAIN()
      endif
    endif ! ATTENUATION
  else
    ! isotropic inner core elements
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  inner core isotropic elements"
      call flush_IMAIN()
    endif

    if (ATTENUATION_VAL) then
      ! only scales for attenuation

      icount_iso = 0

! openmp solver
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ispec, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk, &
#else
!$OMP i,j,k, &
#endif
!$OMP one_minus_sum_beta_use,mul)
!$OMP DO
      do ispec = 1,NSPEC_INNER_CORE

        ! exclude fictitious elements in central cube
        if (idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

        DO_LOOP_IJK

          ! inner core with no anisotropy, use kappav and muv for instance
          ! layer with no anisotropy, use kappav and muv for instance
          mul = muvstore_inner_core(INDEX_IJK,ispec)

          ! use unrelaxed parameters if attenuation
          ! precompute terms for attenuation if needed
          if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
            one_minus_sum_beta_use = one_minus_sum_beta_inner_core(INDEX_IJK,ispec)
          else
            one_minus_sum_beta_use = one_minus_sum_beta_inner_core(1,1,1,ispec)
          endif
          mul = mul * one_minus_sum_beta_use

          muvstore_inner_core(INDEX_IJK,ispec) = mul
        ENDDO_LOOP_IJK
        ! counter
!$OMP ATOMIC
        icount_iso = icount_iso + 1
      enddo
!$OMP ENDDO
!$OMP END PARALLEL

      ! user output
      if (myrank == 0) then
        write(IMAIN,*) "  iso elements  = ",icount_iso
        call flush_IMAIN()
      endif
    endif ! ATTENUATION
  endif ! anisotropic/isotropic inner core

  if (ATTENUATION_VAL) then
    ! only scales for attenuation

    ! since we scale muv and c11,.. stores we must divide with this factor to use the relaxed moduli for the modulus defect
    ! calculation in updating the memory variables

! openmp solver
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ispec,i_SLS, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk, &
#else
!$OMP i,j,k, &
#endif
!$OMP one_minus_sum_beta_use)
!$OMP DO
    do ispec = 1,NSPEC_INNER_CORE
      ! exclude fictitious elements in central cube
      if (idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

      DO_LOOP_IJK
        ! gets factor
        if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
          one_minus_sum_beta_use = one_minus_sum_beta_inner_core(INDEX_IJK,ispec)
        else
          one_minus_sum_beta_use = one_minus_sum_beta_inner_core(1,1,1,ispec)
        endif

        ! corrects factor_common to obtain relaxed moduli in moduli defect
        do i_SLS = 1,N_SLS
          if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
            factor_common_inner_core(INDEX_IJK,i_SLS,ispec) = &
                factor_common_inner_core(INDEX_IJK,i_SLS,ispec) / one_minus_sum_beta_use
          else
            factor_common_inner_core(1,1,1,i_SLS,ispec) = &
                factor_common_inner_core(1,1,1,i_SLS,ispec) / one_minus_sum_beta_use
          endif
        enddo
      ENDDO_LOOP_IJK
    enddo
!$OMP ENDDO
!$OMP END PARALLEL
  endif ! ATTENUATION

  ! safety check
  if (GPU_MODE) then
    print *,'!!! Please make sure to have GPU routines adapted to new elastic tensor arrays !!!'
    stop 'Safety stop for GPU mode with modified elastic elements'
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_elastic_elements
