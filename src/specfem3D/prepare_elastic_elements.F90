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
  real(kind=CUSTOM_REAL) :: rhovphsq,sinphifour,cosphisq,sinphisq,costhetasq,rhovsvsq,sinthetasq, &
        cosphifour,costhetafour,rhovpvsq,sinthetafour,rhovshsq,cosfourphi, &
        costwotheta,cosfourtheta,sintwophisq,costheta,sinphi,sintheta,cosphi, &
        sintwotheta,costwophi,sintwophi,costwothetasq,costwophisq,phi,theta

  ! aniso element
  real(kind=CUSTOM_REAL) :: c11,c12,c13,c22,c23,c33,c44,c55,c66

  real(kind=CUSTOM_REAL) :: two_rhovsvsq,two_rhovshsq ! two_rhovpvsq,two_rhovphsq
  real(kind=CUSTOM_REAL) :: four_rhovsvsq,four_rhovshsq ! four_rhovpvsq,four_rhovphsq

  real(kind=CUSTOM_REAL) :: twoetaminone,etaminone,eta_aniso
  real(kind=CUSTOM_REAL) :: two_eta_aniso,four_eta_aniso,six_eta_aniso

  real(kind=CUSTOM_REAL) :: templ1,templ1_cos,templ2,templ2_cos,templ3,templ3_two,templ3_cos
  real(kind=CUSTOM_REAL) :: kappavl,kappahl,muvl,muhl,mul

  integer :: ispec,iglob,i_SLS,icount,icount_iso

#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
#else
  integer :: i,j,k
#endif

  ! checks if anything to do
  ! GPU kernels still use original arrays
  if (GPU_MODE) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing elastic element arrays"
    call flush_IMAIN()
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
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ispec, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk, &
#else
!$OMP i,j,k, &
#endif
!$OMP minus_sum_beta,mul,c11,c12,c13,c22,c23,c33,c44,c55,c66)
!$OMP DO
      do ispec = 1,NSPEC_CRUST_MANTLE
        ! all elements are fully anisotropic
        DO_LOOP_IJK
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
    endif
  else
    ! iso/tiso mantle elements
    if (myrank == 0) then
      write(IMAIN,*) "  crust/mantle transverse isotropic and isotropic elements"
      call flush_IMAIN()
    endif

    ! crust/mantle
    icount = 0
    icount_iso = 0
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ispec,iglob, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk, &
#else
!$OMP i,j,k, &
#endif
!$OMP kappavl,muvl,kappahl,muhl, &
!$OMP rhovphsq,sinphifour,cosphisq,sinphisq,costhetasq,rhovsvsq,sinthetasq, &
!$OMP cosphifour,costhetafour,rhovpvsq,sinthetafour,rhovshsq,cosfourphi, &
!$OMP costwotheta,cosfourtheta,sintwophisq,costheta,sinphi,sintheta,cosphi, &
!$OMP sintwotheta,costwophi,sintwophi,costwothetasq,costwophisq,phi,theta, &
!$OMP twoetaminone,etaminone,eta_aniso, &
!$OMP two_eta_aniso,four_eta_aniso,six_eta_aniso, &
!$OMP templ1,templ1_cos,templ2,templ2_cos,templ3,templ3_two,templ3_cos, &
!$OMP one_minus_sum_beta_use,mul)
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

          ! precompute some products to reduce the CPU time

          costheta = cos(theta)
          sintheta = sin(theta)
          cosphi = cos(phi)
          sinphi = sin(phi)

          costhetasq = costheta * costheta
          sinthetasq = sintheta * sintheta
          cosphisq = cosphi * cosphi
          sinphisq = sinphi * sinphi

          costhetafour = costhetasq * costhetasq
          sinthetafour = sinthetasq * sinthetasq
          cosphifour = cosphisq * cosphisq
          sinphifour = sinphisq * sinphisq

          costwotheta = cos(2.0_CUSTOM_REAL*theta)
          sintwotheta = sin(2.0_CUSTOM_REAL*theta)
          costwophi = cos(2.0_CUSTOM_REAL*phi)
          sintwophi = sin(2.0_CUSTOM_REAL*phi)

          cosfourtheta = cos(4.0_CUSTOM_REAL*theta)
          cosfourphi = cos(4.0_CUSTOM_REAL*phi)

          costwothetasq = costwotheta * costwotheta

          costwophisq = costwophi * costwophi
          sintwophisq = sintwophi * sintwophi

          etaminone = eta_aniso - 1.0_CUSTOM_REAL
          twoetaminone = 2.0_CUSTOM_REAL * eta_aniso - 1.0_CUSTOM_REAL

          ! precompute some products to reduce the CPU time
          two_eta_aniso = 2.0_CUSTOM_REAL*eta_aniso
          four_eta_aniso = 4.0_CUSTOM_REAL*eta_aniso
          six_eta_aniso = 6.0_CUSTOM_REAL*eta_aniso

          two_rhovsvsq = 2.0_CUSTOM_REAL*rhovsvsq
          two_rhovshsq = 2.0_CUSTOM_REAL*rhovshsq
          four_rhovsvsq = 4.0_CUSTOM_REAL*rhovsvsq
          four_rhovshsq = 4.0_CUSTOM_REAL*rhovshsq

          ! pre-compute temporary values
          templ1 = four_rhovsvsq - rhovpvsq + twoetaminone*rhovphsq - four_eta_aniso*rhovsvsq
          templ1_cos = rhovphsq - rhovpvsq + costwotheta*templ1
          templ2 = four_rhovsvsq - rhovpvsq - rhovphsq + two_eta_aniso*rhovphsq - four_eta_aniso*rhovsvsq
          templ2_cos = rhovpvsq - rhovphsq + costwotheta*templ2
          templ3 = rhovphsq + rhovpvsq - two_eta_aniso*rhovphsq + four_eta_aniso*rhovsvsq
          templ3_two = templ3 - two_rhovshsq - two_rhovsvsq
          templ3_cos = templ3_two + costwotheta*templ2

          ! reordering operations to facilitate compilation, avoiding divisions, using locality for temporary values
          c11store_crust_mantle(INDEX_IJK,ispec) = rhovphsq*sinphifour &
                + 2.0_CUSTOM_REAL*cosphisq*sinphisq* &
                ( rhovphsq*costhetasq + sinthetasq*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq) ) &
                + cosphifour*(rhovphsq*costhetafour &
                  + 2.0_CUSTOM_REAL*costhetasq*sinthetasq*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq) &
                  + rhovpvsq*sinthetafour)

          c12store_crust_mantle(INDEX_IJK,ispec) = 0.25_CUSTOM_REAL*costhetasq &
                *(rhovphsq - two_rhovshsq)*(3.0_CUSTOM_REAL + cosfourphi) &
                - four_rhovshsq*cosphisq*costhetasq*sinphisq &
                + 0.03125_CUSTOM_REAL*rhovphsq*sintwophisq*(11.0_CUSTOM_REAL + cosfourtheta + 4.0*costwotheta) &
                + eta_aniso*sinthetasq*(rhovphsq - two_rhovsvsq) &
                           *(cosphifour + sinphifour + 2.0_CUSTOM_REAL*cosphisq*costhetasq*sinphisq) &
                + rhovpvsq*cosphisq*sinphisq*sinthetafour &
                - rhovsvsq*sintwophisq*sinthetafour

          c13store_crust_mantle(INDEX_IJK,ispec) = 0.125_CUSTOM_REAL*cosphisq &
                *(rhovphsq + six_eta_aniso*rhovphsq + rhovpvsq - four_rhovsvsq &
                      - 12.0_CUSTOM_REAL*eta_aniso*rhovsvsq + cosfourtheta*templ1) &
                + sinphisq*(eta_aniso*costhetasq*(rhovphsq - two_rhovsvsq) + sinthetasq*(rhovphsq - two_rhovshsq))

          ! uses temporary templ1 from c13
          c15store_crust_mantle(INDEX_IJK,ispec) = cosphi*costheta*sintheta* &
                ( 0.5_CUSTOM_REAL*cosphisq* (rhovpvsq - rhovphsq + costwotheta*templ1) &
                  + etaminone*sinphisq*(rhovphsq - two_rhovsvsq))

          c14store_crust_mantle(INDEX_IJK,ispec) = costheta*sinphi*sintheta* &
                ( 0.5_CUSTOM_REAL*cosphisq*(templ2_cos + four_rhovshsq - four_rhovsvsq) &
                  + sinphisq*(etaminone*rhovphsq + 2.0_CUSTOM_REAL*(rhovshsq - eta_aniso*rhovsvsq)) )

          ! uses temporary templ2_cos from c14
          c16store_crust_mantle(INDEX_IJK,ispec) = 0.5_CUSTOM_REAL*cosphi*sinphi*sinthetasq* &
                ( cosphisq*templ2_cos &
                  + 2.0_CUSTOM_REAL*etaminone*sinphisq*(rhovphsq - two_rhovsvsq) )

          c22store_crust_mantle(INDEX_IJK,ispec) = rhovphsq*cosphifour + 2.0_CUSTOM_REAL*cosphisq*sinphisq* &
                (rhovphsq*costhetasq + sinthetasq*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)) &
                + sinphifour* &
                (rhovphsq*costhetafour + 2.0_CUSTOM_REAL*costhetasq*sinthetasq*(eta_aniso*rhovphsq &
                      + two_rhovsvsq - two_eta_aniso*rhovsvsq) + rhovpvsq*sinthetafour)

          ! uses temporary templ1 from c13
          c23store_crust_mantle(INDEX_IJK,ispec) = 0.125_CUSTOM_REAL*sinphisq*(rhovphsq + six_eta_aniso*rhovphsq &
                  + rhovpvsq - four_rhovsvsq - 12.0_CUSTOM_REAL*eta_aniso*rhovsvsq + cosfourtheta*templ1) &
                + cosphisq*(eta_aniso*costhetasq*(rhovphsq - two_rhovsvsq) + sinthetasq*(rhovphsq - two_rhovshsq))

          ! uses temporary templ1 from c13
          c24store_crust_mantle(INDEX_IJK,ispec) = costheta*sinphi*sintheta* &
                ( etaminone*cosphisq*(rhovphsq - two_rhovsvsq) &
                  + 0.5_CUSTOM_REAL*sinphisq*(rhovpvsq - rhovphsq + costwotheta*templ1) )

          ! uses temporary templ2_cos from c14
          c25store_crust_mantle(INDEX_IJK,ispec) = cosphi*costheta*sintheta* &
                ( cosphisq*(etaminone*rhovphsq + 2.0_CUSTOM_REAL*(rhovshsq - eta_aniso*rhovsvsq)) &
                  + 0.5_CUSTOM_REAL*sinphisq*(templ2_cos + four_rhovshsq - four_rhovsvsq) )

          ! uses temporary templ2_cos from c14
          c26store_crust_mantle(INDEX_IJK,ispec) = 0.5_CUSTOM_REAL*cosphi*sinphi*sinthetasq* &
                ( 2.0_CUSTOM_REAL*etaminone*cosphisq*(rhovphsq - two_rhovsvsq) &
                  + sinphisq*templ2_cos )

          c33store_crust_mantle(INDEX_IJK,ispec) = rhovpvsq*costhetafour &
                + 2.0_CUSTOM_REAL*costhetasq*sinthetasq*(two_rhovsvsq + eta_aniso*(rhovphsq - two_rhovsvsq)) &
                + rhovphsq*sinthetafour

          ! uses temporary templ1_cos from c13
          c34store_crust_mantle(INDEX_IJK,ispec) = - 0.25_CUSTOM_REAL*sinphi*sintwotheta*templ1_cos

          ! uses temporary templ1_cos from c34
          c35store_crust_mantle(INDEX_IJK,ispec) = - 0.25_CUSTOM_REAL*cosphi*sintwotheta*templ1_cos

          ! uses temporary templ1_cos from c34
          c36store_crust_mantle(INDEX_IJK,ispec) = - 0.25_CUSTOM_REAL*sintwophi*sinthetasq &
                *(templ1_cos - four_rhovshsq + four_rhovsvsq)

          c44store_crust_mantle(INDEX_IJK,ispec) = cosphisq*(rhovsvsq*costhetasq + rhovshsq*sinthetasq) &
                + sinphisq*(rhovsvsq*costwothetasq + costhetasq*sinthetasq*templ3)

          ! uses temporary templ3 from c44
          c46store_crust_mantle(INDEX_IJK,ispec) = - cosphi*costheta*sintheta* &
                  ( cosphisq*(rhovshsq - rhovsvsq) - 0.5_CUSTOM_REAL*sinphisq*templ3_cos  )

          ! uses templ3 from c46
          c45store_crust_mantle(INDEX_IJK,ispec) = 0.25_CUSTOM_REAL*sintwophi*sinthetasq* &
                (templ3_two + costwotheta*(rhovphsq + rhovpvsq - two_eta_aniso*rhovphsq + 4.0_CUSTOM_REAL*etaminone*rhovsvsq))

          c55store_crust_mantle(INDEX_IJK,ispec) = sinphisq*(rhovsvsq*costhetasq + rhovshsq*sinthetasq) &
                + cosphisq*(rhovsvsq*costwothetasq &
                    + costhetasq*sinthetasq*(rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + four_eta_aniso*rhovsvsq) )

          ! uses temporary templ3_cos from c46
          c56store_crust_mantle(INDEX_IJK,ispec) = costheta*sinphi*sintheta* &
                ( 0.5_CUSTOM_REAL*cosphisq*templ3_cos + sinphisq*(rhovsvsq - rhovshsq) )

          c66store_crust_mantle(INDEX_IJK,ispec) = rhovshsq*costwophisq*costhetasq &
                - 2.0_CUSTOM_REAL*cosphisq*costhetasq*sinphisq*(rhovphsq - two_rhovshsq) &
                + 0.03125_CUSTOM_REAL*rhovphsq*sintwophisq*(11.0_CUSTOM_REAL + 4.0_CUSTOM_REAL*costwotheta + cosfourtheta) &
                - 0.125_CUSTOM_REAL*rhovsvsq*sinthetasq* &
                ( -6.0_CUSTOM_REAL - 2.0_CUSTOM_REAL*costwotheta - 2.0_CUSTOM_REAL*cosfourphi &
                          + cos(4.0_CUSTOM_REAL*phi - 2.0_CUSTOM_REAL*theta) &
                          + cos(2.0_CUSTOM_REAL*(2.0_CUSTOM_REAL*phi + theta)) ) &
                + rhovpvsq*cosphisq*sinphisq*sinthetafour &
                - 0.5_CUSTOM_REAL*eta_aniso*sintwophisq*sinthetafour*(rhovphsq - two_rhovsvsq)
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
    endif
  endif

  ! inner core
  if (ATTENUATION_VAL) then
    ! only scales for attenuation
    if (ANISOTROPIC_INNER_CORE_VAL) then
      ! user output
      if (myrank == 0) then
        write(IMAIN,*) "  inner core anisotropic elements"
        call flush_IMAIN()
      endif

      ! anisotropic inner core elements
      icount = 0
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

    else
      ! isotropic inner core elements
      ! user output
      if (myrank == 0) then
        write(IMAIN,*) "  inner core isotropic elements"
        call flush_IMAIN()
      endif

      icount_iso = 0
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
    endif

    ! since we scale muv and c11,.. stores we must divide with this factor to use the relaxed moduli for the modulus defect
    ! calculation in updating the memory variables
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
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_elastic_elements
