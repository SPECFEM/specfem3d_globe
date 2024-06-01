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

  subroutine compute_forces_inner_core_Dev( NSPEC_STR_OR_ATT,NGLOB,NSPEC_ATT, &
                                            deltat, &
                                            displ_inner_core, &
                                            accel_inner_core, &
                                            iphase, &
                                            R_xx,R_yy,R_xy,R_xz,R_yz, &
                                            R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                            epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                            epsilondev_xz,epsilondev_yz, &
                                            epsilon_trace_over_3, &
                                            alphaval,betaval,gammaval, &
                                            factor_common,vnspec,sum_terms, &
                                            pgrav_inner_core)

! this routine is optimized for NGLLX = NGLLY = NGLLZ = 5 using the Deville et al. (2002) inlined matrix-matrix products

  use constants_solver, only: &
    CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGLLCUBE,NDIM,IFLAG_IN_FICTITIOUS_CUBE, &
    N_SLS,NSPEC_INNER_CORE_STRAIN_ONLY,NSPEC_INNER_CORE, &
    ATT1_VAL,ATT2_VAL,ATT3_VAL, &
    ANISOTROPIC_INNER_CORE_VAL,ATTENUATION_VAL,PARTIAL_PHYS_DISPERSION_ONLY_VAL,GRAVITY_VAL, &
    FULL_GRAVITY_VAL,DISCARD_GCONTRIB, &
    m1,m2

  use specfem_par, only: &
    hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
    wgll_cube, &
    gravity_pre_store => gravity_pre_store_inner_core,gravity_H => gravity_H_inner_core, &
    COMPUTE_AND_STORE_STRAIN,USE_LDDRK

  use specfem_par_innercore, only: &
    deriv => deriv_mapping_inner_core, &
    kappavstore => kappavstore_inner_core, &
    muvstore => muvstore_inner_core, &
    c11store => c11store_inner_core,c12store => c12store_inner_core,c13store => c13store_inner_core, &
    c33store => c33store_inner_core,c44store => c44store_inner_core, &
    ibool => ibool_inner_core, &
    idoubling => idoubling_inner_core, &
    phase_ispec_inner => phase_ispec_inner_inner_core, &
    nspec_outer => nspec_outer_inner_core, &
    nspec_inner => nspec_inner_inner_core

  use specfem_par, only: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D

#ifdef FORCE_VECTORIZATION
  use specfem_par_innercore, only: &
    ibool_inv_tbl => ibool_inv_tbl_inner_core, &
    ibool_inv_st => ibool_inv_st_inner_core, &
    num_globs => num_globs_inner_core, &
    phase_iglob => phase_iglob_inner_core
#endif

  ! full gravity
  use specfem_par_full_gravity, only: &
    gravity_rho => gravity_rho_inner_core

  implicit none

  integer :: NSPEC_STR_OR_ATT,NGLOB,NSPEC_ATT

  ! time step
  real(kind=CUSTOM_REAL) deltat

  ! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: displ_inner_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: accel_inner_core

  ! for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  ! variable lengths for factor_common (and one_minus_sum_beta)
  integer :: vnspec
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,N_SLS,vnspec) :: factor_common
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATT) :: R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATT) :: &
    R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STR_OR_ATT) :: &
    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STRAIN_ONLY) :: epsilon_trace_over_3

  ! work array with contributions
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),intent(out) :: sum_terms

  ! inner/outer element run flag
  integer,intent(in) :: iphase

  ! full gravity
  real(kind=CUSTOM_REAL), dimension(NGLOB),intent(in) :: pgrav_inner_core

  ! local parameters
  ! Deville
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5) :: epsilondev_loc
  real(kind=CUSTOM_REAL) :: fac1,fac2,fac3

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: rho_s_H

  integer :: ispec,iglob
  integer :: num_elements,ispec_p
  integer :: i,j,k
#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
  integer :: ijk_spec, ip, iglob_p
#endif

! ****************************************************
!   big loop over all spectral elements in the solid
! ****************************************************

!  computed_elements = 0

  if (iphase == 1) then
    ! outer elements (halo region)
    num_elements = nspec_outer
  else
    ! inner elements
    num_elements = nspec_inner
  endif

! openmp solver
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED( deriv, &
!$OMP num_elements,iphase,phase_ispec_inner, &
!$OMP ibool,idoubling, &
!$OMP displ_inner_core,accel_inner_core, &
!$OMP c11store,c12store,c13store,c33store,c44store, &
!$OMP muvstore, kappavstore, &
!$OMP factor_common, &
!$OMP alphaval,betaval,gammaval, &
!$OMP R_xx,R_yy,R_xy,R_xz,R_yz, &
!$OMP epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz,epsilon_trace_over_3, &
!$OMP gravity_pre_store,gravity_H, &
!$OMP gravity_rho,pgrav_inner_core, &
!$OMP R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
!$OMP sum_terms, &
#ifdef FORCE_VECTORIZATION
!$OMP ibool_inv_tbl, ibool_inv_st, num_globs, phase_iglob, &
#endif
!$OMP deltat ) &
!$OMP PRIVATE( ispec,ispec_p,i,j,k,iglob, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk_spec,ip,iglob_p, &
!$OMP ijk, &
#endif
!$OMP fac1,fac2,fac3, &
!$OMP tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
!$OMP newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3, &
!$OMP dummyx_loc,dummyy_loc,dummyz_loc, &
!$OMP rho_s_H,epsilondev_loc ) &
!$OMP FIRSTPRIVATE( hprime_xx, hprime_xxT, hprimewgll_xxT, hprimewgll_xx, &
!$OMP wgllwgll_yz_3D, wgllwgll_xz_3D, wgllwgll_xy_3D, wgll_cube, &
!$OMP NSPEC_INNER_CORE,NGLOB, &
!$OMP USE_LDDRK,ANISOTROPIC_INNER_CORE_VAL,GRAVITY_VAL,FULL_GRAVITY_VAL, &
!$OMP ATTENUATION_VAL,PARTIAL_PHYS_DISPERSION_ONLY_VAL, &
!$OMP att1_val,att2_val,att3_val,vnspec, &
!$OMP COMPUTE_AND_STORE_STRAIN )

!$OMP DO SCHEDULE(GUIDED)
  do ispec_p = 1,num_elements

    ! only compute element which belong to current phase (inner or outer elements)
    ispec = phase_ispec_inner(ispec_p,iphase)

    ! exclude fictitious elements in central cube
    if (idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

    ! note: this loop will not fully vectorize because it contains a dependency (through indirect addressing with array ibool())
    !       thus, instead of DO_LOOP_IJK we use do k=..;do j=..;do i=..,
    !       which helps the compiler to unroll the innermost loop
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          dummyx_loc(i,j,k) = displ_inner_core(1,iglob)
          dummyy_loc(i,j,k) = displ_inner_core(2,iglob)
          dummyz_loc(i,j,k) = displ_inner_core(3,iglob)
        enddo
      enddo
    enddo

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for tempx1,..
    call mxm5_3comp_singleA(hprime_xx,m1,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1,m2)
    ! computes 2. matrix multiplication for tempx2,..
    call mxm5_3comp_3dmat_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m1,hprime_xxT,m1,tempx2,tempy2,tempz2,NGLLX)
    ! computes 3. matrix multiplication for tempx1,..
    call mxm5_3comp_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m2,hprime_xxT,tempx3,tempy3,tempz3,m1)

    !
    ! compute either isotropic, transverse isotropic or anisotropic elements
    !
    if (ANISOTROPIC_INNER_CORE_VAL) then
      ! anisotropic elements
      call compute_element_aniso_ic(ispec, &
                                    gravity_pre_store,gravity_H, &
                                    deriv, &
                                    wgll_cube, &
                                    c11store,c12store,c13store,c33store,c44store, &
                                    ibool, &
                                    R_xx,R_yy,R_xy,R_xz,R_yz, &
                                    epsilon_trace_over_3, &
                                    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                    dummyx_loc,dummyy_loc,dummyz_loc, &
                                    epsilondev_loc,rho_s_H)


    else
      ! isotropic elements
      call compute_element_iso_ic(ispec, &
                                  gravity_pre_store,gravity_H, &
                                  deriv, &
                                  wgll_cube, &
                                  kappavstore,muvstore, &
                                  ibool, &
                                  R_xx,R_yy,R_xy,R_xz,R_yz, &
                                  epsilon_trace_over_3, &
                                  tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                  dummyx_loc,dummyy_loc,dummyz_loc, &
                                  epsilondev_loc,rho_s_H)
    endif

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for newtempx1,..
    call mxm5_3comp_singleA(hprimewgll_xxT,m1,tempx1,tempy1,tempz1,newtempx1,newtempy1,newtempz1,m2)
    ! computes 2. matrix multiplication for tempx2,..
    call mxm5_3comp_3dmat_singleB(tempx2,tempy2,tempz2,m1,hprimewgll_xx,m1,newtempx2,newtempy2,newtempz2,NGLLX)
    ! computes 3. matrix multiplication for newtempx3,..
    call mxm5_3comp_singleB(tempx3,tempy3,tempz3,m2,hprimewgll_xx,newtempx3,newtempy3,newtempz3,m1)

    ! sums contributions
    DO_LOOP_IJK
      fac1 = wgllwgll_yz_3D(INDEX_IJK)
      fac2 = wgllwgll_xz_3D(INDEX_IJK)
      fac3 = wgllwgll_xy_3D(INDEX_IJK)
      sum_terms(1,INDEX_IJK,ispec) = - (fac1*newtempx1(INDEX_IJK) + fac2*newtempx2(INDEX_IJK) + fac3*newtempx3(INDEX_IJK))
      sum_terms(2,INDEX_IJK,ispec) = - (fac1*newtempy1(INDEX_IJK) + fac2*newtempy2(INDEX_IJK) + fac3*newtempy3(INDEX_IJK))
      sum_terms(3,INDEX_IJK,ispec) = - (fac1*newtempz1(INDEX_IJK) + fac2*newtempz2(INDEX_IJK) + fac3*newtempz3(INDEX_IJK))
    ENDDO_LOOP_IJK

    ! adds gravity
    if (GRAVITY_VAL) then
      ! full gravity
      if (FULL_GRAVITY_VAL .and. .not. DISCARD_GCONTRIB) then
        call SIEM_solve_element_add_full_gravity(ispec,NSPEC_INNER_CORE,NGLOB,gravity_rho,deriv(:,:,:,:,ispec),ibool, &
                                                 pgrav_inner_core,rho_s_H)
      endif

#ifdef FORCE_VECTORIZATION
      do ijk = 1,NDIM*NGLLCUBE
        sum_terms(ijk,1,1,1,ispec) = sum_terms(ijk,1,1,1,ispec) + rho_s_H(ijk,1,1,1)
      enddo
#else
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            sum_terms(1,INDEX_IJK,ispec) = sum_terms(1,INDEX_IJK,ispec) + rho_s_H(1,INDEX_IJK)
            sum_terms(2,INDEX_IJK,ispec) = sum_terms(2,INDEX_IJK,ispec) + rho_s_H(2,INDEX_IJK)
            sum_terms(3,INDEX_IJK,ispec) = sum_terms(3,INDEX_IJK,ispec) + rho_s_H(3,INDEX_IJK)
          enddo
        enddo
      enddo
#endif
    endif

    ! use Runge-Kutta scheme to march memory variables in time
    ! convention for attenuation
    ! term in xx = 1
    ! term in yy = 2
    ! term in xy = 3
    ! term in xz = 4
    ! term in yz = 5
    ! term in zz not computed since zero trace
    ! This is because we only implement Q_\mu attenuation and not Q_\kappa.
    ! Note that this does *NOT* imply that there is no attenuation for P waves
    ! because for Q_\kappa = infinity one gets (see for instance Dahlen and Tromp (1998)
    ! equation (9.59) page 350): Q_\alpha = Q_\mu * 3 * (V_p/V_s)^2 / 4
    ! therefore Q_\alpha is not zero; for instance for V_p / V_s = sqrt(3)
    ! we get Q_\alpha = (9 / 4) * Q_\mu = 2.25 * Q_\mu
    if (ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL) then
      ! updates R_memory
      if (USE_LDDRK) then
        call compute_element_att_memory_ic_lddrk(ispec,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                                 R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                                 ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec,factor_common, &
                                                 muvstore, &
                                                 epsilondev_loc, &
                                                 deltat)
      else
        call compute_element_att_memory_ic(ispec,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                           ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec,factor_common, &
                                           alphaval,betaval,gammaval, &
                                           muvstore, &
                                           epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                           epsilondev_xz,epsilondev_yz, &
                                           epsilondev_loc)
      endif
    endif

    ! save deviatoric strain for Runge-Kutta scheme
    if (COMPUTE_AND_STORE_STRAIN) then
      epsilondev_xx(:,:,:,ispec) = epsilondev_loc(:,:,:,1)
      epsilondev_yy(:,:,:,ispec) = epsilondev_loc(:,:,:,2)
      epsilondev_xy(:,:,:,ispec) = epsilondev_loc(:,:,:,3)
      epsilondev_xz(:,:,:,ispec) = epsilondev_loc(:,:,:,4)
      epsilondev_yz(:,:,:,ispec) = epsilondev_loc(:,:,:,5)
    endif

  enddo ! of spectral element loop
!$OMP ENDDO

  ! updates acceleration
  ! sum contributions from each element to the global mesh
#ifdef FORCE_VECTORIZATION
!$OMP DO
  do iglob_p = 1,num_globs(iphase)
    iglob = phase_iglob(iglob_p,iphase)
    ! loops over valence points
    do ip = ibool_inv_st(iglob_p,iphase),ibool_inv_st(iglob_p+1,iphase)-1
      ! local 1D index from array ibool
      ijk_spec = ibool_inv_tbl(ip,iphase)
      ! do NOT use array syntax ":" for the three statements below
      ! otherwise most compilers will not be able to vectorize the outer loop
      accel_inner_core(1,iglob) = accel_inner_core(1,iglob) + sum_terms(1,ijk_spec,1,1,1)
      accel_inner_core(2,iglob) = accel_inner_core(2,iglob) + sum_terms(2,ijk_spec,1,1,1)
      accel_inner_core(3,iglob) = accel_inner_core(3,iglob) + sum_terms(3,ijk_spec,1,1,1)
    enddo
  enddo
!$OMP ENDDO
#else
!$OMP DO
  do ispec_p = 1,num_elements
    ! only compute element which belong to current phase (inner or outer elements)
    ispec = phase_ispec_inner(ispec_p,iphase)
    ! exclude fictitious elements in central cube
    if (idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle
    ! updates acceleration (for non-vectorization case)
! note: Critical OpenMP here might degrade performance,
!       especially for a larger number of threads (>8).
!       Using atomic operations can partially help.
#ifndef USE_OPENMP_ATOMIC_INSTEAD_OF_CRITICAL
!$OMP CRITICAL
#endif
! we can force vectorization using a compiler directive here because we know that there is no dependency
! inside a given spectral element, since all the global points of a local elements are different by definition
! (only common points between different elements can be the same)
! IBM, Portland PGI, and Intel and Cray syntax (Intel and Cray are the same)
!IBM* ASSERT (NODEPS)
!pgi$ ivdep
!DIR$ IVDEP
    DO_LOOP_IJK
      iglob = ibool(INDEX_IJK,ispec)
#ifdef USE_OPENMP_ATOMIC_INSTEAD_OF_CRITICAL
!$OMP ATOMIC
#endif
      accel_inner_core(1,iglob) = accel_inner_core(1,iglob) + sum_terms(1,INDEX_IJK,ispec)
#ifdef USE_OPENMP_ATOMIC_INSTEAD_OF_CRITICAL
!$OMP ATOMIC
#endif
      accel_inner_core(2,iglob) = accel_inner_core(2,iglob) + sum_terms(2,INDEX_IJK,ispec)
#ifdef USE_OPENMP_ATOMIC_INSTEAD_OF_CRITICAL
!$OMP ATOMIC
#endif
      accel_inner_core(3,iglob) = accel_inner_core(3,iglob) + sum_terms(3,INDEX_IJK,ispec)
    ENDDO_LOOP_IJK
#ifndef USE_OPENMP_ATOMIC_INSTEAD_OF_CRITICAL
!$OMP END CRITICAL
#endif
  enddo
!$OMP ENDDO
#endif

!$OMP END PARALLEL

  contains

!--------------------------------------------------------------------------------------------
!
! matrix-matrix multiplications
!
! subroutines adapted from Deville, Fischer and Mund, High-order methods
! for incompressible fluid flow, Cambridge University Press (2002),
! pages 386 and 389 and Figure 8.3.1
!
!--------------------------------------------------------------------------------------------
!
! note: the matrix-matrix multiplications are used for very small matrices ( 5 x 5 x 5 elements);
!       thus, calling external optimized libraries for these multiplications are in general slower
!
! please leave the routines here to help compilers inlining the code

  subroutine mxm5_3comp_singleA(A,n1,B1,B2,B3,C1,C2,C3,n3)

! we can force inlining (Intel compiler)
#if defined __INTEL_COMPILER
!DIR$ ATTRIBUTES FORCEINLINE :: mxm5_3comp_singleA
#else
! cray
! note: with Cray Fortran versions >= 14 on Frontier, inlining this routine together with optimization -O3 leads to problems.
!       for now, will avoid inlining by this directive INLINENEVER to allow for default compilation,
!       otherwise the compilation flag -hipa0 would need to be added to suppress all inlining as well.
!!DIR$ INLINEALWAYS mxm5_3comp_singleA
!DIR$ INLINENEVER mxm5_3comp_singleA
#endif

! 3 different arrays for x/y/z-components, 2-dimensional arrays (25,5)/(5,25), same B matrix for all 3 component arrays

  use constants_solver, only: CUSTOM_REAL

#ifdef USE_XSMM
  !use my_libxsmm, only: libxsmm_smm_5_25_5
  use my_libxsmm, only: xmm1,libxsmm_mmcall_abc => libxsmm_smmcall_abc
#endif

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(5,n3),intent(in) :: B1,B2,B3
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C1,C2,C3

  ! local parameters
  integer :: i,j

#ifdef USE_XSMM
  ! matrix-matrix multiplication C = alpha A * B + beta C
  ! with A(n1,n2) 5x5-matrix, B(n2,n3) 5x25-matrix and C(n1,n3) 5x25-matrix
  ! static version using MNK="5 25, 5" ALPHA=1 BETA=0
  !call libxsmm_smm_5_25_5(a=A, b=B1, c=C1, pa=A, pb=B2, pc=C2)
  !call libxsmm_smm_5_25_5(a=A, b=B2, c=C2, pa=A, pb=B3, pc=C3)
  !call libxsmm_smm_5_25_5(a=A, b=B3, c=C3, pa=A, pb=B1, pc=C1) ! with dummy prefetch
  ! dispatch
  call libxsmm_mmcall_abc(xmm1, A, B1, C1)
  call libxsmm_mmcall_abc(xmm1, A, B2, C2)
  call libxsmm_mmcall_abc(xmm1, A, B3, C3)
  return
#endif

  ! matrix-matrix multiplication
  do j = 1,n3
!DIR$ IVDEP
#if defined __INTEL_COMPILER
!DIR$ SIMD
#endif
    do i = 1,n1
      C1(i,j) =  A(i,1) * B1(1,j) &
               + A(i,2) * B1(2,j) &
               + A(i,3) * B1(3,j) &
               + A(i,4) * B1(4,j) &
               + A(i,5) * B1(5,j)

      C2(i,j) =  A(i,1) * B2(1,j) &
               + A(i,2) * B2(2,j) &
               + A(i,3) * B2(3,j) &
               + A(i,4) * B2(4,j) &
               + A(i,5) * B2(5,j)

      C3(i,j) =  A(i,1) * B3(1,j) &
               + A(i,2) * B3(2,j) &
               + A(i,3) * B3(3,j) &
               + A(i,4) * B3(4,j) &
               + A(i,5) * B3(5,j)
    enddo
  enddo

  end subroutine mxm5_3comp_singleA


!--------------------------------------------------------------------------------------------

  subroutine mxm5_3comp_singleB(A1,A2,A3,n1,B,C1,C2,C3,n3)

! we can force inlining (Intel compiler)
#if defined __INTEL_COMPILER
!DIR$ ATTRIBUTES FORCEINLINE :: mxm5_3comp_singleB
#else
! cray
! note: with Cray Fortran versions >= 14 on Frontier, inlining this routine together with optimization -O3 leads to problems.
!       for now, will avoid inlining by this directive INLINENEVER to allow for default compilation,
!       otherwise the compilation flag -hipa0 would need to be added to suppress all inlining as well.
!!DIR$ INLINEALWAYS mxm5_3comp_singleB
!DIR$ INLINENEVER mxm5_3comp_singleB
#endif

! 3 different arrays for x/y/z-components, 2-dimensional arrays (25,5)/(5,25), same B matrix for all 3 component arrays

  use constants_solver, only: CUSTOM_REAL

#ifdef USE_XSMM
  !use my_libxsmm, only: libxsmm_smm_25_5_5
  use my_libxsmm, only: xmm2,libxsmm_mmcall_abc => libxsmm_smmcall_abc
#endif

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5),intent(in) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(5,n3),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C1,C2,C3

  ! local parameters
  integer :: i,j

#ifdef USE_XSMM
  ! matrix-matrix multiplication C = alpha A * B + beta C
  ! with A(n1,n2) 25x5-matrix, B(n2,n3) 5x5-matrix and C(n1,n3) 25x5-matrix
  ! static version
  !call libxsmm_smm_25_5_5(a=A1, b=B, c=C1, pa=A2, pb=B, pc=C2)
  !call libxsmm_smm_25_5_5(a=A2, b=B, c=C2, pa=A3, pb=B, pc=C3)
  !call libxsmm_smm_25_5_5(a=A3, b=B, c=C3, pa=A1, pb=B, pc=C1)
  ! dispatch
  call libxsmm_mmcall_abc(xmm2, A1, B, C1)
  call libxsmm_mmcall_abc(xmm2, A2, B, C2)
  call libxsmm_mmcall_abc(xmm2, A3, B, C3)
  return
#endif

  ! matrix-matrix multiplication
  do j = 1,n3
!DIR$ IVDEP
#if defined __INTEL_COMPILER
!DIR$ SIMD
#endif
    do i = 1,n1
      C1(i,j) =  A1(i,1) * B(1,j) &
               + A1(i,2) * B(2,j) &
               + A1(i,3) * B(3,j) &
               + A1(i,4) * B(4,j) &
               + A1(i,5) * B(5,j)

      C2(i,j) =  A2(i,1) * B(1,j) &
               + A2(i,2) * B(2,j) &
               + A2(i,3) * B(3,j) &
               + A2(i,4) * B(4,j) &
               + A2(i,5) * B(5,j)

      C3(i,j) =  A3(i,1) * B(1,j) &
               + A3(i,2) * B(2,j) &
               + A3(i,3) * B(3,j) &
               + A3(i,4) * B(4,j) &
               + A3(i,5) * B(5,j)
    enddo
  enddo

  end subroutine mxm5_3comp_singleB


!--------------------------------------------------------------------------------------------

  subroutine mxm5_3comp_3dmat_singleB(A1,A2,A3,n1,B,n2,C1,C2,C3,n3)

! we can force inlining (Intel compiler)
#if defined __INTEL_COMPILER
!DIR$ ATTRIBUTES FORCEINLINE :: mxm5_3comp_3dmat_singleB
#else
! cray
! note: with Cray Fortran versions >= 14 on Frontier, inlining this routine together with optimization -O3 leads to problems.
!       for now, will avoid inlining by this directive INLINENEVER to allow for default compilation,
!       otherwise the compilation flag -hipa0 would need to be added to suppress all inlining as well.
!!DIR$ INLINEALWAYS mxm5_3comp_3dmat_singleB
!DIR$ INLINENEVER mxm5_3comp_3dmat_singleB
#endif

! 3 different arrays for x/y/z-components, 3-dimensional arrays (5,5,5), same B matrix for all 3 component arrays

  use constants_solver, only: CUSTOM_REAL

! note: on CPUs like Haswell or Sandy Bridge, the following will slow down computations
!       however, on Intel Phi (KNC) it is still helpful (speedup +3%)
#if defined(XSMM_FORCE_EVEN_IF_SLOWER) || ( defined(XSMM) && defined(__MIC__) )
  !use my_libxsmm, only: libxsmm_smm_5_5_5
  use my_libxsmm, only: xmm3,libxsmm_mmcall_abc => libxsmm_smmcall_abc
#endif

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,5,n3),intent(in) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(5,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C1,C2,C3

  ! local parameters
  integer :: i,j,k

#if defined(XSMM_FORCE_EVEN_IF_SLOWER) || ( defined(XSMM) && defined(__MIC__) )
  ! matrix-matrix multiplication C = alpha A * B + beta C
  ! with A(n1,n2,n4) 5x5x5-matrix, B(n2,n3) 5x5-matrix and C(n1,n3,n4) 5x5x5-matrix
  ! static version
  !do k = 1,5
  !  call libxsmm_call(xmm3, A1(:,:,k), B, C1(:,:,k))
  !  call libxsmm_call(xmm3, A2(:,:,k), B, C2(:,:,k))
  !  call libxsmm_call(xmm3, A3(:,:,k), B, C3(:,:,k))
  !enddo
  ! dispatch
  do k = 1,5
    call libxsmm_mmcall_abc(xmm3, A1(1,1,k), B, C1(1,1,k))
    call libxsmm_mmcall_abc(xmm3, A2(1,1,k), B, C2(1,1,k))
    call libxsmm_mmcall_abc(xmm3, A3(1,1,k), B, C3(1,1,k))
  enddo
  ! unrolled
  !call libxsmm_smm_5_5_5(a=A1(1,1,1), b=B, c=C1(1,1,1),pa=A1(1,1,1+1), pb=B, pc=C1(1,1,1+1))
  !call libxsmm_smm_5_5_5(a=A1(1,1,2), b=B, c=C1(1,1,2),pa=A1(1,1,2+1), pb=B, pc=C1(1,1,2+1))
  !call libxsmm_smm_5_5_5(a=A1(1,1,3), b=B, c=C1(1,1,3),pa=A1(1,1,3+1), pb=B, pc=C1(1,1,3+1))
  !call libxsmm_smm_5_5_5(a=A1(1,1,4), b=B, c=C1(1,1,4),pa=A1(1,1,4+1), pb=B, pc=C1(1,1,4+1))
  !call libxsmm_smm_5_5_5(a=A1(1,1,5), b=B, c=C1(1,1,5),pa=A2(1,1,1), pb=B, pc=C2(1,1,1))

  !call libxsmm_smm_5_5_5(a=A2(1,1,1), b=B, c=C2(1,1,1),pa=A2(1,1,1+1), pb=B, pc=C2(1,1,1+1))
  !call libxsmm_smm_5_5_5(a=A2(1,1,2), b=B, c=C2(1,1,2),pa=A2(1,1,2+1), pb=B, pc=C2(1,1,2+1))
  !call libxsmm_smm_5_5_5(a=A2(1,1,3), b=B, c=C2(1,1,3),pa=A2(1,1,3+1), pb=B, pc=C2(1,1,3+1))
  !call libxsmm_smm_5_5_5(a=A2(1,1,4), b=B, c=C2(1,1,4),pa=A2(1,1,4+1), pb=B, pc=C2(1,1,4+1))
  !call libxsmm_smm_5_5_5(a=A2(1,1,5), b=B, c=C2(1,1,5),pa=A3(1,1,1), pb=B, pc=C3(1,1,1))

  !call libxsmm_smm_5_5_5(a=A3(1,1,1), b=B, c=C3(1,1,1),pa=A3(1,1,1+1), pb=B, pc=C3(1,1,1+1))
  !call libxsmm_smm_5_5_5(a=A3(1,1,2), b=B, c=C3(1,1,2),pa=A3(1,1,2+1), pb=B, pc=C3(1,1,2+1))
  !call libxsmm_smm_5_5_5(a=A3(1,1,3), b=B, c=C3(1,1,3),pa=A3(1,1,3+1), pb=B, pc=C3(1,1,3+1))
  !call libxsmm_smm_5_5_5(a=A3(1,1,4), b=B, c=C3(1,1,4),pa=A3(1,1,4+1), pb=B, pc=C3(1,1,4+1))
  !call libxsmm_smm_5_5_5(a=A3(1,1,5), b=B, c=C3(1,1,5),pa=A3(1,1,5), pb=B, pc=C3(1,1,5))
  return
#endif

  ! matrix-matrix multiplication
  do k = 1,n3
    do j = 1,n2
!DIR$ IVDEP
#if defined __INTEL_COMPILER
!DIR$ SIMD
#endif
      do i = 1,n1
        C1(i,j,k) =  A1(i,1,k) * B(1,j) &
                   + A1(i,2,k) * B(2,j) &
                   + A1(i,3,k) * B(3,j) &
                   + A1(i,4,k) * B(4,j) &
                   + A1(i,5,k) * B(5,j)

        C2(i,j,k) =  A2(i,1,k) * B(1,j) &
                   + A2(i,2,k) * B(2,j) &
                   + A2(i,3,k) * B(3,j) &
                   + A2(i,4,k) * B(4,j) &
                   + A2(i,5,k) * B(5,j)

        C3(i,j,k) =  A3(i,1,k) * B(1,j) &
                   + A3(i,2,k) * B(2,j) &
                   + A3(i,3,k) * B(3,j) &
                   + A3(i,4,k) * B(4,j) &
                   + A3(i,5,k) * B(5,j)
      enddo
    enddo
  enddo

  end subroutine mxm5_3comp_3dmat_singleB

  end subroutine compute_forces_inner_core_Dev

