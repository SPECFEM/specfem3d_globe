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
  real(kind=CUSTOM_REAL) :: kappavl,kappahl,muvl,muhl,mul

  ! aniso element
  real(kind=CUSTOM_REAL) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                            c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! for rotations
  double precision :: g11,g12,g13,g14,g15,g16,g22,g23,g24,g25,g26, &
                      g33,g34,g35,g36,g44,g45,g46,g55,g56,g66
  double precision :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66
  double precision :: phi_dble,theta_dble
  double precision :: A_dble,F_dble,L_dble,N_dble !,C_dble

  integer :: ispec,iglob,i_SLS,icount,icount_iso
  integer :: count_glob,count_iso_glob

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
!$OMP PRIVATE(ispec,iglob, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk, &
#else
!$OMP i,j,k, &
#endif
!$OMP minus_sum_beta,mul,muvl,muhl,eta_aniso, &
!$OMP theta_dble,phi_dble, &
!$OMP d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26,d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
!$OMP g11,g12,g13,g14,g15,g16,g22,g23,g24,g25,g26,g33,g34,g35,g36,g44,g45,g46,g55,g56,g66, &
!$OMP A_dble,F_dble,L_dble,N_dble)
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
          ! for radial axis symmetry and transversely isotropic media this simplifies to:
          ! A = C11 = rho * vph**2 = kappah + 4/3 muh
          ! C = C33 = rho * vpv**2 = kappav + 4/3 muv
          ! N = C66 = rho * vsh**2 = muh
          ! L = C44 = rho * vsv**2 = muv
          ! F = C13 = eta * (A - 2*L)
          ! and
          ! C12 = C11 - 2 C66 = A - 2*N = rho * vph**2 - 2 * rho * vsh**2 = kappah + 4/3 muh - 2 muh = kappah - 2/3 muh
          ! C22 = C11
          ! C23 = C13
          ! C55 = C44
          ! C66 = N = rho * vsh**2 = muh = (C11-C12)/2
          !
          ! shear moduli scaling:
          ! muvl = muvl * one_minus_sum_beta_use = muvl*( 1 - sum_beta) = muvl - muvl*sum_beta = muvl + muvl*minus_sum_beta
          ! muhl = muhl * one_minus_sum_beta_use
          !
          ! rho vsv**2 = muv -> L = L * (1 - sum_beta) = L - L * sum_beta = L + muv * minus_sum_beta
          ! rho vsh**2 = muh -> N = N + muh * minus_sum_beta
          !
          ! rho vph**2 = A = kappah + 4/3 muh -> kappah + 4/3 muh * (1 - sum_beta) = kappah + 4/3 muh - 4/3 muh*sum_beta
          ! rho vpv**2 = C = kappav + 4/3 muv -> kappav + 4/3 muv - 4/3 muv*sum_beta
          !
          ! -> scaling for Love parameters:
          ! A' = A - 4/3 muh*sum_beta
          ! C' = C - 4/3 muv*sum_beta
          ! N' = N * (1 - sum_beta) = N - muh*sum_beta
          ! L' = L * (1 - sum_beta) = L - muv*sum_beta
          !
          ! -> scaling for tensor elements:
          ! c11' = c11 - 4/3 muh*sum_beta  ! A
          ! c33' = c33 - 4/3 muv*sum_beta  ! C
          ! c66' = c66 - muh*sum_beta      ! N
          ! c44' = c44 - muv*sum_beta      ! L
          !
          ! c55' = c55 - muv*sum_beta      ! to be consistent with tiso case c55 == c44
          ! c22' = c22 - 4/3 muh*sum_beta  ! to be consistent with tiso case c22 == c11
          !
          ! c13' = F' = eta*(A' - 2*L') = eta*(A - 4/3 muh*sum_beta - 2*(L - muv*sum_beta))
          !                             = c13 + eta*( 2*muv*sum_beta - 4/3 muh*sum_beta)
          !                             = c13 + eta*( 4/3 muh * minus_sum_beta - 2 muv * minus_sum_beta)
          ! c23' = c23 + eta*( 2*muv*sum_beta - 4/3 muh*sum_beta) ! be consistent with tiso case c23 == c13
          !
          ! c12' = A' - 2*N' = A - 4/3 muh*sum_beta - 2*(N - muh*sum_beta) = A - 2*N - 4/3 muh*sum_beta + 2 muh*sum_beta
          !                                                                = c12 + 2/3 muh*sum_beta


          ! precompute terms for attenuation if needed
          if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
            minus_sum_beta =  one_minus_sum_beta_crust_mantle(INDEX_IJK,ispec) - 1.0_CUSTOM_REAL
          else
            minus_sum_beta =  one_minus_sum_beta_crust_mantle(1,1,1,ispec) - 1.0_CUSTOM_REAL
          endif

          ! original routine...
          !
          ! shifts moduli to unrelaxed ones
          ! this is the original routine to shift moduli, using only muv (from c44) to scale
          !  c11 = c11store_crust_mantle(INDEX_IJK,ispec)
          !  c12 = c12store_crust_mantle(INDEX_IJK,ispec)
          !  c13 = c13store_crust_mantle(INDEX_IJK,ispec)
          !  c22 = c22store_crust_mantle(INDEX_IJK,ispec)
          !  c23 = c23store_crust_mantle(INDEX_IJK,ispec)
          !  c33 = c33store_crust_mantle(INDEX_IJK,ispec)
          !  c44 = c44store_crust_mantle(INDEX_IJK,ispec)
          !  c55 = c55store_crust_mantle(INDEX_IJK,ispec)
          !  c66 = c66store_crust_mantle(INDEX_IJK,ispec)
          !
          !  mul = c44 * minus_sum_beta
          !
          !  c11 = c11 + FOUR_THIRDS * mul ! * minus_sum_beta * mul
          !  c12 = c12 - TWO_THIRDS * mul
          !  c13 = c13 - TWO_THIRDS * mul
          !  c22 = c22 + FOUR_THIRDS * mul
          !  c23 = c23 - TWO_THIRDS * mul
          !  c33 = c33 + FOUR_THIRDS * mul
          !  c44 = c44 + mul
          !  c55 = c55 + mul
          !  c66 = c66 + mul
          !
          !  ! stores unrelaxed factors
          !  c11store_crust_mantle(INDEX_IJK,ispec) = c11
          !  c12store_crust_mantle(INDEX_IJK,ispec) = c12
          !  c13store_crust_mantle(INDEX_IJK,ispec) = c13
          !  c22store_crust_mantle(INDEX_IJK,ispec) = c22
          !  c23store_crust_mantle(INDEX_IJK,ispec) = c23
          !  c33store_crust_mantle(INDEX_IJK,ispec) = c33
          !  c44store_crust_mantle(INDEX_IJK,ispec) = c44
          !  c55store_crust_mantle(INDEX_IJK,ispec) = c55
          !  c66store_crust_mantle(INDEX_IJK,ispec) = c66

          ! new routine
          ! unrelaxed moduli shift using muv and muh scaling
          ! also, needs adaptation in save_kernels.f90 to shift back accordingly

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

          ! new: shifts shear moduli by separating muv and muh factors.
          !      still needs rotations to rotate back and forth from SPECFEM global axis to a radial symmetry axis
          !      since this shift assumes a radial symmetry
          A_dble = 0.125d0 * (3.d0 * d11 + 3.d0 * d22 + 2.d0 * d12 + 4.d0 * d66)
          !C_dble = d33
          N_dble = 0.125d0 * (d11 + d22 - 2.d0 * d12 + 4.d0 * d66)
          L_dble = 0.5d0 * (d44 + d55)
          F_dble = 0.5d0 * (d13 + d23)

          eta_aniso = real(F_dble / (A_dble - 2.d0*L_dble),kind=CUSTOM_REAL)   ! eta = F / (A-2L)

          muvl = real(L_dble * minus_sum_beta,kind=CUSTOM_REAL)     ! c44 - > L - > muv
          muhl = real(N_dble * minus_sum_beta,kind=CUSTOM_REAL)     ! c66 - > N - > muh

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

          ! for solving memory-variables, modulus defect \delta \mu_l (Komatitsch, 2002, eq. (11) & (13))
          ! note: for solving the memory variables, we will only use the modulus defect
          !       associated with muv. this is consistent with the implementation for tiso below.
          !
          !       however, to properly account for shear attenuation, one might have to add also
          !       memory-variables for a modulus defect associated with muh.
          muvstore_crust_mantle(INDEX_IJK,ispec) = real(d44,kind=CUSTOM_REAL)

        ENDDO_LOOP_IJK

        ! counter
!$OMP ATOMIC
        icount = icount + 1
      enddo
!$OMP ENDDO
!$OMP END PARALLEL

      ! stats
      call sum_all_i(icount,count_glob)

      ! user output
      if (myrank == 0) then
        write(IMAIN,*) "  aniso elements = ",count_glob
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
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ispec,iglob, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk, &
#else
!$OMP i,j,k, &
#endif
!$OMP one_minus_sum_beta_use,mul,muvl,muhl,eta_aniso, &
!$OMP kappavl,kappahl,rhovpvsq,rhovphsq,rhovsvsq,rhovshsq, &
!$OMP c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
!$OMP theta,phi)
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
          muhstore_crust_mantle(INDEX_IJK,ispec) = muhl

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

    ! stats
    call sum_all_i(icount,count_glob)
    call sum_all_i(icount_iso,count_iso_glob)

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  tiso elements = ",count_glob
      write(IMAIN,*) "  iso elements  = ",count_iso_glob
      call flush_IMAIN()
    endif
    if (icount > NSPECMAX_TISO_MANTLE) stop 'Error invalid number of tiso elements in prepare_timerun_aniso'
    if (icount_iso > NSPECMAX_ISO_MANTLE) stop 'Error invalid number of iso elements in prepare_timerun_aniso'
  endif ! ANISOTROPIC_3D_MANTLE_VAL

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
          !   > mul = muvstore(INDEX_IJK,ispec) * minus_sum_beta
          ! this would be following the implementation from above for fully anisotropic elements...
          !
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

          ! for solving memory-variables, modulus defect \delta \mu_l (Komatitsch, 2002, eq. (11) & (13))
          muvstore_inner_core(INDEX_IJK,ispec) = c44
        ENDDO_LOOP_IJK
        ! counter
!$OMP ATOMIC
        icount = icount + 1
      enddo
!$OMP ENDDO
!$OMP END PARALLEL

      ! stats
      call sum_all_i(icount,count_glob)

      ! user output
      if (myrank == 0) then
        write(IMAIN,*) "  aniso elements = ",count_glob
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

      ! stats
      call sum_all_i(icount_iso,count_iso_glob)

      ! user output
      if (myrank == 0) then
        write(IMAIN,*) "  iso elements  = ",count_iso_glob
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

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_elastic_elements
