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

  subroutine compute_forces_crust_mantle_Dev( NSPEC_STR_OR_ATT,NGLOB,NSPEC_ATT, &
                                              deltat, &
                                              displ_crust_mantle, &
                                              accel_crust_mantle, &
                                              iphase, &
                                              R_xx,R_yy,R_xy,R_xz,R_yz, &
                                              R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                              epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                              epsilondev_xz,epsilondev_yz, &
                                              epsilon_trace_over_3, &
                                              alphaval,betaval,gammaval, &
                                              factor_common,vnspec,sum_terms )

! this routine is optimized for NGLLX = NGLLY = NGLLZ = 5 using the Deville et al. (2002) inlined matrix-matrix products

  use constants_solver, only: &
    CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGLLCUBE,NDIM, &
    N_SLS,NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_CRUST_MANTLE, &
    ATT1_VAL,ATT2_VAL,ATT3_VAL, &
    ANISOTROPIC_3D_MANTLE_VAL,ATTENUATION_VAL,PARTIAL_PHYS_DISPERSION_ONLY_VAL,GRAVITY_VAL, &
    m1,m2

  use specfem_par, only: &
    hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
    wgll_cube, &
    gravity_pre_store => gravity_pre_store_crust_mantle,gravity_H => gravity_H_crust_mantle, &
    COMPUTE_AND_STORE_STRAIN,USE_LDDRK

  use specfem_par_crustmantle, only: &
    deriv => deriv_mapping_crust_mantle, &
    kappavstore => kappavstore_crust_mantle, &
    muvstore => muvstore_crust_mantle, &
    c11store => c11store_crust_mantle,c12store => c12store_crust_mantle,c13store => c13store_crust_mantle, &
    c14store => c14store_crust_mantle,c15store => c15store_crust_mantle,c16store => c16store_crust_mantle, &
    c22store => c22store_crust_mantle,c23store => c23store_crust_mantle,c24store => c24store_crust_mantle, &
    c25store => c25store_crust_mantle,c26store => c26store_crust_mantle,c33store => c33store_crust_mantle, &
    c34store => c34store_crust_mantle,c35store => c35store_crust_mantle,c36store => c36store_crust_mantle, &
    c44store => c44store_crust_mantle,c45store => c45store_crust_mantle,c46store => c46store_crust_mantle, &
    c55store => c55store_crust_mantle,c56store => c56store_crust_mantle,c66store => c66store_crust_mantle, &
    ibool => ibool_crust_mantle, &
    ispec_is_tiso => ispec_is_tiso_crust_mantle, &
    phase_ispec_inner => phase_ispec_inner_crust_mantle, &
    nspec_outer => nspec_outer_crust_mantle, &
    nspec_inner => nspec_inner_crust_mantle

  use specfem_par, only: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D

#ifdef FORCE_VECTORIZATION
  ! optimized arrays
  use specfem_par_crustmantle, only: &
    ibool_inv_tbl => ibool_inv_tbl_crust_mantle, &
    ibool_inv_st => ibool_inv_st_crust_mantle, &
    num_globs => num_globs_crust_mantle, &
    phase_iglob => phase_iglob_crust_mantle
#endif

!daniel: att - debug
!  use specfem_par, only: it,NSTEP

  implicit none

  integer,intent(in) :: NSPEC_STR_OR_ATT,NGLOB,NSPEC_ATT

  ! time step
  real(kind=CUSTOM_REAL),intent(in) :: deltat

  ! displacement, velocity and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB),intent(in) :: displ_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB),intent(inout) :: accel_crust_mantle

  ! variable sized array variables
  integer,intent(in) :: vnspec

  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATT),intent(inout) :: R_xx,R_yy,R_xy,R_xz,R_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATT),intent(inout) :: &
    R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STR_OR_ATT),intent(inout) :: &
    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STRAIN_ONLY),intent(inout) :: epsilon_trace_over_3

  ! [alpha,beta,gamma]val reduced to N_SLS and factor_common to N_SLS*NUM_NODES
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,N_SLS,vnspec),intent(in) :: factor_common
  real(kind=CUSTOM_REAL), dimension(N_SLS),intent(in) :: alphaval,betaval,gammaval

  ! work array with contributions
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),intent(out) :: sum_terms

  ! inner/outer element run flag
  integer,intent(in) :: iphase

  ! local parameters

  ! Deville
  ! manually inline the calls to the Deville et al. (2002) routines
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5) :: epsilondev_loc
  real(kind=CUSTOM_REAL) :: fac1,fac2,fac3

  ! for gravity
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: rho_s_H

  integer :: ispec,iglob
  integer :: num_elements,ispec_p
  integer :: i,j,k
#ifdef FORCE_VECTORIZATION
  integer :: ijk_spec,ip,iglob_p,ijk
#endif

  !integer,parameter :: NGLL2 = NGLLY * NGLLZ
  !integer,parameter :: NGLL3 = NGLLX * NGLLY * NGLLZ

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
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP SHARED( deriv, &
!$OMP num_elements,iphase,phase_ispec_inner, &
!$OMP ibool,ispec_is_tiso, &
!$OMP displ_crust_mantle,accel_crust_mantle, &
!$OMP wgll_cube,hprime_xxt,hprime_xx,hprimewgll_xx,hprimewgll_xxT, &
!$OMP wgllwgll_xy_3D, wgllwgll_xz_3D, wgllwgll_yz_3D, &
!$OMP c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
!$OMP c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
!$OMP c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
!$OMP kappavstore,muvstore, &
!$OMP vnspec, &
!$OMP factor_common, &
!$OMP alphaval,betaval,gammaval, &
!$OMP R_xx,R_yy,R_xy,R_xz,R_yz, &
!$OMP epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz,epsilon_trace_over_3, &
!$OMP gravity_pre_store,gravity_H, &
!$OMP USE_LDDRK, &
!$OMP R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
!$OMP sum_terms, &
#ifdef FORCE_VECTORIZATION
!$OMP ibool_inv_tbl, ibool_inv_st, num_globs, phase_iglob, &
#endif
!$OMP deltat,COMPUTE_AND_STORE_STRAIN ) &
!$OMP PRIVATE( ispec,ispec_p,i,j,k,iglob, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk_spec,ip,iglob_p, &
!$OMP ijk, &
#endif
!$OMP fac1,fac2,fac3, &
!$OMP tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
!$OMP newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3, &
!$OMP dummyx_loc,dummyy_loc,dummyz_loc, &
!$OMP rho_s_H,epsilondev_loc )

!$OMP DO SCHEDULE(GUIDED)
  do ispec_p = 1,num_elements

    ! only compute elements which belong to current phase (inner or outer elements)
    ispec = phase_ispec_inner(ispec_p,iphase)

    ! note: this loop will not fully vectorize because it contains a dependency (through indirect addressing with array ibool())
    !       thus, instead of DO_LOOP_IJK we use do k=..;do j=..;do i=..,
    !       which helps the compiler to unroll the innermost loop
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)

          dummyx_loc(i,j,k) = displ_crust_mantle(1,iglob)
          dummyy_loc(i,j,k) = displ_crust_mantle(2,iglob)
          dummyz_loc(i,j,k) = displ_crust_mantle(3,iglob)
        enddo
      enddo
    enddo

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

#ifdef DANIEL_TEST_LOOP
    ! loop over single x/y/z-component, to test if cache utilization is better
    ! x-comp
!DIR$ FORCEINLINE
    call mxm5_3comp_singleA_1(hprime_xx,m1,dummyx_loc,tempx1,m2)
!DIR$ FORCEINLINE
    call mxm5_3comp_3dmat_singleB_1(dummyx_loc,m1,hprime_xxT,m1,tempx2,NGLLX)
!DIR$ FORCEINLINE
    call mxm5_3comp_singleB_1(dummyx_loc,m2,hprime_xxT,tempx3,m1)
    ! y-comp
!DIR$ FORCEINLINE
    call mxm5_3comp_singleA_1(hprime_xx,m1,dummyy_loc,tempy1,m2)
!DIR$ FORCEINLINE
    call mxm5_3comp_3dmat_singleB_1(dummyy_loc,m1,hprime_xxT,m1,tempy2,NGLLX)
!DIR$ FORCEINLINE
    call mxm5_3comp_singleB_1(dummyy_loc,m2,hprime_xxT,tempy3,m1)
    ! z-comp
!DIR$ FORCEINLINE
    call mxm5_3comp_singleA_1(hprime_xx,m1,dummyz_loc,tempz1,m2)
!DIR$ FORCEINLINE
    call mxm5_3comp_3dmat_singleB_1(dummyz_loc,m1,hprime_xxT,m1,tempz2,NGLLX)
!DIR$ FORCEINLINE
    call mxm5_3comp_singleB_1(dummyz_loc,m2,hprime_xxT,tempz3,m1)
#else
    ! computes 1. matrix multiplication for tempx1,..
    call mxm5_3comp_singleA(hprime_xx,m1,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1,m2)
    ! computes 2. matrix multiplication for tempx2,..
    call mxm5_3comp_3dmat_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m1,hprime_xxT,m1,tempx2,tempy2,tempz2,NGLLX)
    ! computes 3. matrix multiplication for tempx3,..
    call mxm5_3comp_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m2,hprime_xxT,tempx3,tempy3,tempz3,m1)
#endif

    !
    ! compute either isotropic, transverse isotropic or anisotropic elements
    !
    ! note: for OpenMP, then non-anisotropic case leads an imbalance depending on which thread is treating what kind of element.
    !       one could avoid it by separating the loop and looping over tiso elements and iso element separately.
    !       using Intel VTune, it estimates a potential gain of ~4% so it's left the way it is for now...
    if (ANISOTROPIC_3D_MANTLE_VAL) then
      ! anisotropic element
      call compute_element_aniso(ispec, &
                                 gravity_pre_store,gravity_H, &
                                 deriv, &
                                 wgll_cube, &
                                 c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                                 c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                                 c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                                 ibool, &
                                 R_xx,R_yy,R_xy,R_xz,R_yz, &
                                 epsilon_trace_over_3, &
                                 tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                 dummyx_loc,dummyy_loc,dummyz_loc, &
                                 epsilondev_loc,rho_s_H)
    else
       if (.not. ispec_is_tiso(ispec)) then
          ! isotropic element
          call compute_element_iso(ispec, &
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
       else
          ! transverse isotropic element
          call compute_element_tiso(ispec, &
                                    gravity_pre_store,gravity_H, &
                                    deriv, &
                                    wgll_cube, &
                                    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                                    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                                    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                                    ibool, &
                                    R_xx,R_yy,R_xy,R_xz,R_yz, &
                                    epsilon_trace_over_3, &
                                    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                    dummyx_loc,dummyy_loc,dummyz_loc, &
                                    epsilondev_loc,rho_s_H)
       endif ! .not. ispec_is_tiso
    endif

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

#ifdef DANIEL_TEST_LOOP
    ! loop over single x/y/z-component, to test if cache utilization is better
    ! x-comp
!DIR$ FORCEINLINE
    call mxm5_3comp_singleA_1(hprimewgll_xxT,m1,tempx1,newtempx1,m2)
!DIR$ FORCEINLINE
    call mxm5_3comp_3dmat_singleB_1(tempx2,m1,hprimewgll_xx,m1,newtempx2,NGLLX)
!DIR$ FORCEINLINE
    call mxm5_3comp_singleB_1(tempx3,m2,hprimewgll_xx,newtempx3,m1)
    ! y-comp
!DIR$ FORCEINLINE
    call mxm5_3comp_singleA_1(hprimewgll_xxT,m1,tempy1,newtempy1,m2)
!DIR$ FORCEINLINE
    call mxm5_3comp_3dmat_singleB_1(tempy2,m1,hprimewgll_xx,m1,newtempy2,NGLLX)
!DIR$ FORCEINLINE
    call mxm5_3comp_singleB_1(tempy3,m2,hprimewgll_xx,newtempy3,m1)
    ! z-comp
!DIR$ FORCEINLINE
    call mxm5_3comp_singleA_1(hprimewgll_xxT,m1,tempz1,newtempz1,m2)
!DIR$ FORCEINLINE
    call mxm5_3comp_3dmat_singleB_1(tempz2,m1,hprimewgll_xx,m1,newtempz2,NGLLX)
!DIR$ FORCEINLINE
    call mxm5_3comp_singleB_1(tempz3,m2,hprimewgll_xx,newtempz3,m1)
#else
    ! computes 1. matrix multiplication for newtempx1,..
    call mxm5_3comp_singleA(hprimewgll_xxT,m1,tempx1,tempy1,tempz1,newtempx1,newtempy1,newtempz1,m2)
    ! computes 2. matrix multiplication for tempx2,..
    call mxm5_3comp_3dmat_singleB(tempx2,tempy2,tempz2,m1,hprimewgll_xx,m1,newtempx2,newtempy2,newtempz2,NGLLX)
    ! computes 3. matrix multiplication for newtempx3,..
    call mxm5_3comp_singleB(tempx3,tempy3,tempz3,m2,hprimewgll_xx,newtempx3,newtempy3,newtempz3,m1)
#endif

    ! sums contributions
    DO_LOOP_IJK
      fac1 = wgllwgll_yz_3D(INDEX_IJK)
      fac2 = wgllwgll_xz_3D(INDEX_IJK)
      fac3 = wgllwgll_xy_3D(INDEX_IJK)
      sum_terms(1,INDEX_IJK,ispec) = - (fac1*newtempx1(INDEX_IJK) + fac2*newtempx2(INDEX_IJK) + fac3*newtempx3(INDEX_IJK))
      sum_terms(2,INDEX_IJK,ispec) = - (fac1*newtempy1(INDEX_IJK) + fac2*newtempy2(INDEX_IJK) + fac3*newtempy3(INDEX_IJK))
      sum_terms(3,INDEX_IJK,ispec) = - (fac1*newtempz1(INDEX_IJK) + fac2*newtempz2(INDEX_IJK) + fac3*newtempz3(INDEX_IJK))
    ENDDO_LOOP_IJK

    ! adds gravity terms
    if (GRAVITY_VAL) then
#ifdef FORCE_VECTORIZATION
      do ijk = 1,NDIM*NGLLCUBE
        sum_terms(ijk,1,1,1,ispec) = sum_terms(ijk,1,1,1,ispec) + rho_s_H(ijk,1,1,1)
      enddo
#else
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            sum_terms(1,i,j,k,ispec) = sum_terms(1,i,j,k,ispec) + rho_s_H(1,i,j,k)
            sum_terms(2,i,j,k,ispec) = sum_terms(2,i,j,k,ispec) + rho_s_H(2,i,j,k)
            sum_terms(3,i,j,k,ispec) = sum_terms(3,i,j,k,ispec) + rho_s_H(3,i,j,k)
          enddo
        enddo
      enddo
#endif
    endif

    ! update memory variables based upon the Runge-Kutta scheme
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
        call compute_element_att_memory_cm_lddrk(ispec,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                                 R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                                 ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec,factor_common, &
                                                 muvstore, &
                                                 epsilondev_loc, &
                                                 deltat)
      else
        call compute_element_att_memory_cm(ispec,R_xx,R_yy,R_xy,R_xz,R_yz, &
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

  enddo ! of spectral element loop NSPEC_CRUST_MANTLE
!$OMP ENDDO

  ! updates acceleration
#ifdef FORCE_VECTORIZATION
  ! updates for vectorized case
  ! loops over all global nodes in this phase (inner/outer)
!$OMP DO
  do iglob_p = 1,num_globs(iphase)
    ! global node index
    iglob = phase_iglob(iglob_p,iphase)
    ! loops over valence points
    do ip = ibool_inv_st(iglob_p,iphase),ibool_inv_st(iglob_p+1,iphase)-1
      ! local 1D index from array ibool
      ijk_spec = ibool_inv_tbl(ip,iphase)

      ! do NOT use array syntax ":" for the three statements below otherwise most compilers
      ! will not be able to vectorize the outer loop
      accel_crust_mantle(1,iglob) = accel_crust_mantle(1,iglob) + sum_terms(1,ijk_spec,1,1,1)
      accel_crust_mantle(2,iglob) = accel_crust_mantle(2,iglob) + sum_terms(2,ijk_spec,1,1,1)
      accel_crust_mantle(3,iglob) = accel_crust_mantle(3,iglob) + sum_terms(3,ijk_spec,1,1,1)
    enddo
  enddo
!$OMP ENDDO

#else
    ! updates for non-vectorization case
!$OMP DO SCHEDULE(GUIDED)
  do ispec_p = 1,num_elements
    ! only compute elements which belong to current phase (inner or outer elements)
    ispec = phase_ispec_inner(ispec_p,iphase)
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
      accel_crust_mantle(1,iglob) = accel_crust_mantle(1,iglob) + sum_terms(1,INDEX_IJK,ispec)
#ifdef USE_OPENMP_ATOMIC_INSTEAD_OF_CRITICAL
!$OMP ATOMIC
#endif
      accel_crust_mantle(2,iglob) = accel_crust_mantle(2,iglob) + sum_terms(2,INDEX_IJK,ispec)
#ifdef USE_OPENMP_ATOMIC_INSTEAD_OF_CRITICAL
!$OMP ATOMIC
#endif
      accel_crust_mantle(3,iglob) = accel_crust_mantle(3,iglob) + sum_terms(3,INDEX_IJK,ispec)
    ENDDO_LOOP_IJK
#ifndef USE_OPENMP_ATOMIC_INSTEAD_OF_CRITICAL
!$OMP END CRITICAL
#endif
  enddo
!$OMP ENDDO
#endif

!$OMP END PARALLEL


! kept here for reference: updating for non-vectorized case
!
! timing example:
!         update here will take: 1m 25s
!         update in ispec-loop : 1m 20s
! thus a 6% increase...
!
! this gets even slower for MIC and a large number of OpenMP threads
! (for best performance use the vectorized version)
!
! ! similar as above but with i/j/k/ispec-indexing
!  do iglob_p = 1,num_globs(iphase)
!    ! global node index
!    iglob = phase_iglob(iglob_p,iphase)
!    ! loops over valence points
!    do ip = ibool_inv_st(iglob_p,iphase),ibool_inv_st(iglob_p+1,iphase)-1
!      ! local 1D index from array ibool
!      ijk_spec = ibool_inv_tbl(ip,iphase)
!
!      ! uses (i,j,k,ispec) indexing!
!      !
!      ! converts to i/j/k/ispec-indexing (starting from 0)
!      ijk_spec = ijk_spec - 1
!
!      ispec = int(ijk_spec / NGLL3)
!      ijk = ijk_spec - ispec * NGLL3
!
!      k = int(ijk / NGLL2)
!      j = int((ijk - k * NGLL2) / NGLLX)
!      i = ijk - k * NGLL2 - j * NGLLX
!
!      ! converts back to indexing starting from 1
!      ispec = ispec + 1
!      i = i + 1
!      j = j + 1
!      k = k + 1
!
!      ! checks
!      !if (i < 1 .or. i > NGLLX .or. j < 1 .or. j > NGLLY .or. k < 1 .or. k > NGLLZ .or. ispec < 1 .or. ispec > NSPEC) then
!      !  print *,'Error i/j/k-index: ',i,j,k,ispec,'from',ijk_spec,ijk
!      !  stop 'Error invalid i-index'
!      !endif
!
!      accel_crust_mantle(1,iglob) = accel_crust_mantle(1,iglob) + sum_terms(1,i,j,k,ispec)
!      accel_crust_mantle(2,iglob) = accel_crust_mantle(2,iglob) + sum_terms(2,i,j,k,ispec)
!      accel_crust_mantle(3,iglob) = accel_crust_mantle(3,iglob) + sum_terms(3,i,j,k,ispec)
!    enddo
!  enddo

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
!DIR$ INLINEALWAYS mxm5_3comp_singleA
#endif

! 3 different arrays for x/y/z-components, 2-dimensional arrays (25,5)/(5,25), same B matrix for all 3 component arrays

  use constants_solver, only: CUSTOM_REAL

#ifdef USE_XSMM
  use my_libxsmm, only: libxsmm_smm_5_25_5
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
  call libxsmm_smm_5_25_5(a=A, b=B1, c=C1, pa=A, pb=B2, pc=C2)
  call libxsmm_smm_5_25_5(a=A, b=B2, c=C2, pa=A, pb=B3, pc=C3)
  call libxsmm_smm_5_25_5(a=A, b=B3, c=C3, pa=A, pb=B1, pc=C1) ! with dummy prefetch
  return
#endif

  ! matrix-matrix multiplication
  do j = 1,n3
!DIR$ IVDEP
!DIR$ SIMD
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
!DIR$ INLINEALWAYS mxm5_3comp_singleB
#endif

! 3 different arrays for x/y/z-components, 2-dimensional arrays (25,5)/(5,25), same B matrix for all 3 component arrays

  use constants_solver, only: CUSTOM_REAL

#ifdef USE_XSMM
  use my_libxsmm, only: libxsmm_smm_25_5_5
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
  call libxsmm_smm_25_5_5(a=A1, b=B, c=C1, pa=A2, pb=B, pc=C2)
  call libxsmm_smm_25_5_5(a=A2, b=B, c=C2, pa=A3, pb=B, pc=C3)
  call libxsmm_smm_25_5_5(a=A3, b=B, c=C3, pa=A1, pb=B, pc=C1)
  return
#endif

  ! matrix-matrix multiplication
  do j = 1,n3
!DIR$ IVDEP
!DIR$ SIMD
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
!DIR$ INLINEALWAYS mxm5_3comp_3dmat_singleB
#endif

! 3 different arrays for x/y/z-components, 3-dimensional arrays (5,5,5), same B matrix for all 3 component arrays

  use constants_solver, only: CUSTOM_REAL

! note: on CPUs like Haswell or Sandy Bridge, the following will slow down computations
!       however, on Intel Phi (KNC) it is still helpful (speedup +3%)
#if defined(XSMM_FORCE_EVEN_IF_SLOWER) || ( defined(XSMM) && defined(__MIC__) )
  use my_libxsmm, only: libxsmm_smm_5_5_5
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

  ! unrolled
  call libxsmm_smm_5_5_5(a=A1(1,1,1), b=B, c=C1(1,1,1),pa=A1(1,1,1+1), pb=B, pc=C1(1,1,1+1))
  call libxsmm_smm_5_5_5(a=A1(1,1,2), b=B, c=C1(1,1,2),pa=A1(1,1,2+1), pb=B, pc=C1(1,1,2+1))
  call libxsmm_smm_5_5_5(a=A1(1,1,3), b=B, c=C1(1,1,3),pa=A1(1,1,3+1), pb=B, pc=C1(1,1,3+1))
  call libxsmm_smm_5_5_5(a=A1(1,1,4), b=B, c=C1(1,1,4),pa=A1(1,1,4+1), pb=B, pc=C1(1,1,4+1))
  call libxsmm_smm_5_5_5(a=A1(1,1,5), b=B, c=C1(1,1,5),pa=A2(1,1,1), pb=B, pc=C2(1,1,1))

  call libxsmm_smm_5_5_5(a=A2(1,1,1), b=B, c=C2(1,1,1),pa=A2(1,1,1+1), pb=B, pc=C2(1,1,1+1))
  call libxsmm_smm_5_5_5(a=A2(1,1,2), b=B, c=C2(1,1,2),pa=A2(1,1,2+1), pb=B, pc=C2(1,1,2+1))
  call libxsmm_smm_5_5_5(a=A2(1,1,3), b=B, c=C2(1,1,3),pa=A2(1,1,3+1), pb=B, pc=C2(1,1,3+1))
  call libxsmm_smm_5_5_5(a=A2(1,1,4), b=B, c=C2(1,1,4),pa=A2(1,1,4+1), pb=B, pc=C2(1,1,4+1))
  call libxsmm_smm_5_5_5(a=A2(1,1,5), b=B, c=C2(1,1,5),pa=A3(1,1,1), pb=B, pc=C3(1,1,1))

  call libxsmm_smm_5_5_5(a=A3(1,1,1), b=B, c=C3(1,1,1),pa=A3(1,1,1+1), pb=B, pc=C3(1,1,1+1))
  call libxsmm_smm_5_5_5(a=A3(1,1,2), b=B, c=C3(1,1,2),pa=A3(1,1,2+1), pb=B, pc=C3(1,1,2+1))
  call libxsmm_smm_5_5_5(a=A3(1,1,3), b=B, c=C3(1,1,3),pa=A3(1,1,3+1), pb=B, pc=C3(1,1,3+1))
  call libxsmm_smm_5_5_5(a=A3(1,1,4), b=B, c=C3(1,1,4),pa=A3(1,1,4+1), pb=B, pc=C3(1,1,4+1))
  call libxsmm_smm_5_5_5(a=A3(1,1,5), b=B, c=C3(1,1,5),pa=A3(1,1,5), pb=B, pc=C3(1,1,5))
  return
#endif

  ! matrix-matrix multiplication
  do k = 1,n3
    do j = 1,n2
!DIR$ IVDEP
!DIR$ SIMD
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


!--------------------------------------------------------------------------------------------

#ifdef DANIEL_TEST_LOOP

! loops over single x/y/z-component
! test if cache utilization is better

  subroutine mxm5_3comp_singleA_1(A,n1,B,C,n3)
  use constants_solver, only: CUSTOM_REAL
#ifdef USE_XSMM
  use my_libxsmm, only: libxsmm_smm_5_25_5
#endif
  implicit none
  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(5,n3),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C
  ! local parameters
  integer :: i,j
#ifdef USE_XSMM
  ! matrix-matrix multiplication C = alpha A * B + beta C
  ! with A(n1,n2) 5x5-matrix, B(n2,n3) 5x25-matrix and C(n1,n3) 5x25-matrix
  ! static version using MNK="5 25, 5" ALPHA=1 BETA=0
  call libxsmm_smm_5_25_5(a=A, b=B, c=C)
  return
#endif
  ! matrix-matrix multiplication
  do j = 1,n3
!dir$ ivdep
    do i = 1,n1
      C(i,j) =  A(i,1) * B(1,j) &
              + A(i,2) * B(2,j) &
              + A(i,3) * B(3,j) &
              + A(i,4) * B(4,j) &
              + A(i,5) * B(5,j)
    enddo
  enddo

  end subroutine mxm5_3comp_singleA_1


  subroutine mxm5_3comp_singleB_1(A,n1,B,C,n3)
  use constants_solver, only: CUSTOM_REAL
#ifdef USE_XSMM
  use my_libxsmm, only: libxsmm_smm_25_5_5
#endif
  implicit none
  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(5,n3),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C
  ! local parameters
  integer :: i,j
#ifdef USE_XSMM
  ! matrix-matrix multiplication C = alpha A * B + beta C
  ! with A(n1,n2) 25x5-matrix, B(n2,n3) 5x5-matrix and C(n1,n3) 25x5-matrix
  ! static version
  call libxsmm_smm_25_5_5(a=A, b=B, c=C)
  return
#endif
  ! matrix-matrix multiplication
  do j = 1,n3
!dir$ ivdep
    do i = 1,n1
      C(i,j) =  A(i,1) * B(1,j) &
              + A(i,2) * B(2,j) &
              + A(i,3) * B(3,j) &
              + A(i,4) * B(4,j) &
              + A(i,5) * B(5,j)
    enddo
  enddo
  end subroutine mxm5_3comp_singleB_1


  subroutine mxm5_3comp_3dmat_singleB_1(A,n1,B,n2,C,n3)
  use constants_solver, only: CUSTOM_REAL
#if defined(XSMM_FORCE_EVEN_IF_SLOWER) || ( defined(XSMM) && defined(__MIC__) )
  use my_libxsmm, only: libxsmm_smm_5_5_5
#endif
  implicit none
  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,5,n3),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(5,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C
  ! local parameters
  integer :: i,j,k
#if defined(XSMM_FORCE_EVEN_IF_SLOWER) || ( defined(XSMM) && defined(__MIC__) )
  ! matrix-matrix multiplication C = alpha A * B + beta C
  ! with A(n1,n2,n4) 5x5x5-matrix, B(n2,n3) 5x5-matrix and C(n1,n3,n4) 5x5x5-matrix
  call libxsmm_smm_5_5_5(a=A(1,1,1), b=B, c=C(1,1,1))
  call libxsmm_smm_5_5_5(a=A(1,1,2), b=B, c=C(1,1,2))
  call libxsmm_smm_5_5_5(a=A(1,1,3), b=B, c=C(1,1,3))
  call libxsmm_smm_5_5_5(a=A(1,1,4), b=B, c=C(1,1,4))
  call libxsmm_smm_5_5_5(a=A(1,1,5), b=B, c=C(1,1,5))
  return
#endif
  ! matrix-matrix multiplication
  do k = 1,n3
    do j = 1,n2
!dir$ ivdep
      do i = 1,n1
        C(i,j,k) =  A(i,1,k) * B(1,j) &
                  + A(i,2,k) * B(2,j) &
                  + A(i,3,k) * B(3,j) &
                  + A(i,4,k) * B(4,j) &
                  + A(i,5,k) * B(5,j)
      enddo
    enddo
  enddo
  end subroutine mxm5_3comp_3dmat_singleB_1
#endif


!--------------------------------------------------------------------------------------------



  end subroutine compute_forces_crust_mantle_Dev



! please leave for reference...
!!--------------------------------------------------------------------------------------------
!!
!! Deville et al. 2002
!! Higher-Order Methods for Incompressible Fluid Flow
!!
!! subroutines adapted from Deville, Fischer and Mund, High-order methods
!! for incompressible fluid flow, Cambridge University Press (2002),
!! pages 386 and 389 and Figure 8.3.1
!!
!!--------------------------------------------------------------------------------------------
!
!! matrix - matrix multiplications
!
!! single component routines
!
!  subroutine mxm(A,n1,B,n2,C,n3)
!
!  use constants_solver, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n2,n3
!  real(kind=CUSTOM_REAL),dimension(n1,n2) :: A
!  real(kind=CUSTOM_REAL),dimension(n2,n3) :: B
!  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C
!
!  ! chooses optimized version
!  select case (n2)
!
!  case (4)
!    call mxm4(A,n1,B,C,n3)
!
!  case (5)
!    call mxm5(A,n1,B,C,n3)
!
!  case (6)
!    call mxm6(A,n1,B,C,n3)
!
!  case (7)
!    call mxm7(A,n1,B,C,n3)
!
!  case (8)
!    call mxm8(A,n1,B,C,n3)
!
!  case default
!    call mxmN(A,n1,B,n2,C,n3)
!
!  end select
!
!  end subroutine mxm
!
!!--------------------------------------------------------------------------------------------
!
!  subroutine mxm4(A,n1,B,C,n3)
!
!  use constants_solver, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n3
!  real(kind=CUSTOM_REAL),dimension(n1,4) :: A
!  real(kind=CUSTOM_REAL),dimension(4,n3) :: B
!  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C
!
!  ! local parameters
!  integer :: i,j
!
!  ! matrix-matrix multiplication
!  do j = 1,n3
!    do i = 1,n1
!      C(i,j) =  A(i,1) * B(1,j) &
!              + A(i,2) * B(2,j) &
!              + A(i,3) * B(3,j) &
!              + A(i,4) * B(4,j)
!    enddo
!  enddo
!
!  end subroutine mxm4
!
!!--------------------------------------------------------------------------------------------
!
!  subroutine mxm5(A,n1,B,C,n3)
!
!  use constants_solver, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n3
!  real(kind=CUSTOM_REAL),dimension(n1,5) :: A
!  real(kind=CUSTOM_REAL),dimension(5,n3) :: B
!  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C
!
!  ! local parameters
!  integer :: i,j
!
!  ! matrix-matrix multiplication
!  do j = 1,n3
!    do i = 1,n1
!      C(i,j) =  A(i,1) * B(1,j) &
!              + A(i,2) * B(2,j) &
!              + A(i,3) * B(3,j) &
!              + A(i,4) * B(4,j) &
!              + A(i,5) * B(5,j)
!    enddo
!  enddo
!
!  end subroutine mxm5
!
!!--------------------------------------------------------------------------------------------
!
!  subroutine mxm6(A,n1,B,C,n3)
!
!  use constants_solver, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n3
!  real(kind=CUSTOM_REAL),dimension(n1,6) :: A
!  real(kind=CUSTOM_REAL),dimension(6,n3) :: B
!  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C
!
!  ! local parameters
!  integer :: i,j
!
!  ! matrix-matrix multiplication
!  do j = 1,n3
!    do i = 1,n1
!      C(i,j) =  A(i,1) * B(1,j) &
!              + A(i,2) * B(2,j) &
!              + A(i,3) * B(3,j) &
!              + A(i,4) * B(4,j) &
!              + A(i,5) * B(5,j) &
!              + A(i,6) * B(6,j)
!    enddo
!  enddo
!
!  end subroutine mxm6
!
!!--------------------------------------------------------------------------------------------
!
!  subroutine mxm7(A,n1,B,C,n3)
!
!  use constants_solver, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n3
!  real(kind=CUSTOM_REAL),dimension(n1,7) :: A
!  real(kind=CUSTOM_REAL),dimension(7,n3) :: B
!  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C
!
!  ! local parameters
!  integer :: i,j
!
!  ! matrix-matrix multiplication
!  do j = 1,n3
!    do i = 1,n1
!      C(i,j) =  A(i,1) * B(1,j) &
!              + A(i,2) * B(2,j) &
!              + A(i,3) * B(3,j) &
!              + A(i,4) * B(4,j) &
!              + A(i,5) * B(5,j) &
!              + A(i,6) * B(6,j) &
!              + A(i,7) * B(7,j)
!    enddo
!  enddo
!
!  end subroutine mxm7
!
!
!!--------------------------------------------------------------------------------------------
!
!  subroutine mxm8(A,n1,B,C,n3)
!
!  use constants_solver, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n3
!  real(kind=CUSTOM_REAL),dimension(n1,8) :: A
!  real(kind=CUSTOM_REAL),dimension(8,n3) :: B
!  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C
!
!  ! local parameters
!  integer :: i,j
!
!  ! matrix-matrix multiplication
!  do j = 1,n3
!    do i = 1,n1
!      C(i,j) =  A(i,1) * B(1,j) &
!              + A(i,2) * B(2,j) &
!              + A(i,3) * B(3,j) &
!              + A(i,4) * B(4,j) &
!              + A(i,5) * B(5,j) &
!              + A(i,6) * B(6,j) &
!              + A(i,7) * B(7,j) &
!              + A(i,8) * B(8,j)
!    enddo
!  enddo
!
!  end subroutine mxm8
!
!
!!--------------------------------------------------------------------------------------------
!
!  subroutine mxmN(A,n1,B,n2,C,n3)
!
!  use constants_solver, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n2,n3
!  real(kind=CUSTOM_REAL),dimension(n1,n2) :: A
!  real(kind=CUSTOM_REAL),dimension(n2,n3) :: B
!  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C
!
!  ! local parameters
!  integer :: i,j,k
!  real(kind=CUSTOM_REAL) :: tmp
!
!  ! general matrix-matrix multiplication
!  do j = 1,n3
!    do k = 1,n2
!      tmp = B(k,j)
!      do i = 1,n1
!        C(i,j) = C(i,j) + A(i,k) * tmp
!      enddo
!    enddo
!  enddo
!
!  end subroutine mxmN
!
!
!!----------------------------------------------------------------------------------------------
!
!
!
!! 3-component routines: combines arrays A1,A2,A3 which correspond to 3 different components x/y/z
!
!  subroutine mxm_3comp(A1,A2,A3,n1,B1,B2,B3,n2,C1,C2,C3,n3)
!
!  use constants_solver, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n2,n3
!  real(kind=CUSTOM_REAL),dimension(n1,n2) :: A1,A2,A3
!  real(kind=CUSTOM_REAL),dimension(n2,n3) :: B1,B2,B3
!  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C1,C2,C3
!
!  ! chooses optimized version
!  select case (n2)
!
!  case (4)
!    call mxm4_3comp(A1,A2,A3,n1,B1,B2,B3,C1,C2,C3,n3)
!
!  case (5)
!    call mxm5_3comp(A1,A2,A3,n1,B1,B2,B3,C1,C2,C3,n3)
!
!  case (6)
!    call mxm6_3comp(A1,A2,A3,n1,B1,B2,B3,C1,C2,C3,n3)
!
!  case (7)
!    call mxm7_3comp(A1,A2,A3,n1,B1,B2,B3,C1,C2,C3,n3)
!
!  case (8)
!    call mxm8_3comp(A1,A2,A3,n1,B1,B2,B3,C1,C2,C3,n3)
!
!  case default
!    call mxmN_3comp(A1,A2,A3,n1,B1,B2,B3,n2,C1,C2,C3,n3)
!
!  end select
!
!  end subroutine mxm_3comp
!
!!--------------------------------------------------------------------------------------------
!
!  subroutine mxm4_3comp(A1,A2,A3,n1,B1,B2,B3,C1,C2,C3,n3)
!
!  use constants_solver, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n3
!  real(kind=CUSTOM_REAL),dimension(n1,4) :: A1,A2,A3
!  real(kind=CUSTOM_REAL),dimension(4,n3) :: B1,B2,B3
!  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C1,C2,C3
!
!  ! local parameters
!  integer :: i,j
!
!  ! matrix-matrix multiplication
!  do j = 1,n3
!    do i = 1,n1
!      C1(i,j) =  A1(i,1) * B1(1,j) &
!               + A1(i,2) * B1(2,j) &
!               + A1(i,3) * B1(3,j) &
!               + A1(i,4) * B1(4,j)
!
!      C2(i,j) =  A2(i,1) * B2(1,j) &
!               + A2(i,2) * B2(2,j) &
!               + A2(i,3) * B2(3,j) &
!               + A2(i,4) * B2(4,j)
!
!      C3(i,j) =  A3(i,1) * B3(1,j) &
!               + A3(i,2) * B3(2,j) &
!               + A3(i,3) * B3(3,j) &
!               + A3(i,4) * B3(4,j)
!    enddo
!  enddo
!
!  end subroutine mxm4_3comp
!
!!--------------------------------------------------------------------------------------------
!
!  subroutine mxm5_3comp(A1,A2,A3,n1,B1,B2,B3,C1,C2,C3,n3)
!
!  use constants_solver, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n3
!  real(kind=CUSTOM_REAL),dimension(n1,5) :: A1,A2,A3
!  real(kind=CUSTOM_REAL),dimension(5,n3) :: B1,B2,B3
!  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C1,C2,C3
!
!  ! local parameters
!  integer :: i,j
!
!  ! matrix-matrix multiplication
!  do j = 1,n3
!    do i = 1,n1
!      C1(i,j) =  A1(i,1) * B1(1,j) &
!               + A1(i,2) * B1(2,j) &
!               + A1(i,3) * B1(3,j) &
!               + A1(i,4) * B1(4,j) &
!               + A1(i,5) * B1(5,j)
!
!      C2(i,j) =  A2(i,1) * B2(1,j) &
!               + A2(i,2) * B2(2,j) &
!               + A2(i,3) * B2(3,j) &
!               + A2(i,4) * B2(4,j) &
!               + A2(i,5) * B2(5,j)
!
!      C3(i,j) =  A3(i,1) * B3(1,j) &
!               + A3(i,2) * B3(2,j) &
!               + A3(i,3) * B3(3,j) &
!               + A3(i,4) * B3(4,j) &
!               + A3(i,5) * B3(5,j)
!    enddo
!  enddo
!
!  end subroutine mxm5_3comp
!
!
!
!!--------------------------------------------------------------------------------------------
!
!  subroutine mxm6_3comp(A1,A2,A3,n1,B1,B2,B3,C1,C2,C3,n3)
!
!  use constants_solver, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n3
!  real(kind=CUSTOM_REAL),dimension(n1,6) :: A1,A2,A3
!  real(kind=CUSTOM_REAL),dimension(6,n3) :: B1,B2,B3
!  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C1,C2,C3
!
!  ! local parameters
!  integer :: i,j
!
!  ! matrix-matrix multiplication
!  do j = 1,n3
!    do i = 1,n1
!      C1(i,j) =  A1(i,1) * B1(1,j) &
!               + A1(i,2) * B1(2,j) &
!               + A1(i,3) * B1(3,j) &
!               + A1(i,4) * B1(4,j) &
!               + A1(i,5) * B1(5,j) &
!               + A1(i,6) * B1(6,j)
!
!      C2(i,j) =  A2(i,1) * B2(1,j) &
!               + A2(i,2) * B2(2,j) &
!               + A2(i,3) * B2(3,j) &
!               + A2(i,4) * B2(4,j) &
!               + A2(i,5) * B2(5,j) &
!               + A2(i,6) * B2(6,j)
!
!      C3(i,j) =  A3(i,1) * B3(1,j) &
!               + A3(i,2) * B3(2,j) &
!               + A3(i,3) * B3(3,j) &
!               + A3(i,4) * B3(4,j) &
!               + A3(i,5) * B3(5,j) &
!               + A3(i,6) * B3(6,j)
!    enddo
!  enddo
!
!  end subroutine mxm6_3comp
!
!!--------------------------------------------------------------------------------------------
!
!  subroutine mxm7_3comp(A1,A2,A3,n1,B1,B2,B3,C1,C2,C3,n3)
!
!  use constants_solver, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n3
!  real(kind=CUSTOM_REAL),dimension(n1,7) :: A1,A2,A3
!  real(kind=CUSTOM_REAL),dimension(7,n3) :: B1,B2,B3
!  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C1,C2,C3
!
!  ! local parameters
!  integer :: i,j
!
!  ! matrix-matrix multiplication
!  do j = 1,n3
!    do i = 1,n1
!      C1(i,j) =  A1(i,1) * B1(1,j) &
!               + A1(i,2) * B1(2,j) &
!               + A1(i,3) * B1(3,j) &
!               + A1(i,4) * B1(4,j) &
!               + A1(i,5) * B1(5,j) &
!               + A1(i,6) * B1(6,j) &
!               + A1(i,7) * B1(7,j)
!
!      C2(i,j) =  A2(i,1) * B2(1,j) &
!               + A2(i,2) * B2(2,j) &
!               + A2(i,3) * B2(3,j) &
!               + A2(i,4) * B2(4,j) &
!               + A2(i,5) * B2(5,j) &
!               + A2(i,6) * B2(6,j) &
!               + A2(i,7) * B2(7,j)
!
!      C3(i,j) =  A3(i,1) * B3(1,j) &
!               + A3(i,2) * B3(2,j) &
!               + A3(i,3) * B3(3,j) &
!               + A3(i,4) * B3(4,j) &
!               + A3(i,5) * B3(5,j) &
!               + A3(i,6) * B3(6,j) &
!               + A3(i,7) * B3(7,j)
!    enddo
!  enddo
!
!  end subroutine mxm7_3comp
!
!
!!--------------------------------------------------------------------------------------------
!
!  subroutine mxm8_3comp(A1,A2,A3,n1,B1,B2,B3,C1,C2,C3,n3)
!
!  use constants_solver, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n3
!  real(kind=CUSTOM_REAL),dimension(n1,8) :: A1,A2,A3
!  real(kind=CUSTOM_REAL),dimension(8,n3) :: B1,B2,B3
!  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C1,C2,C3
!
!  ! local parameters
!  integer :: i,j
!
!  ! matrix-matrix multiplication
!  do j = 1,n3
!    do i = 1,n1
!      C1(i,j) =  A1(i,1) * B1(1,j) &
!               + A1(i,2) * B1(2,j) &
!               + A1(i,3) * B1(3,j) &
!               + A1(i,4) * B1(4,j) &
!               + A1(i,5) * B1(5,j) &
!               + A1(i,6) * B1(6,j) &
!               + A1(i,7) * B1(7,j) &
!               + A1(i,8) * B1(8,j)
!
!      C2(i,j) =  A2(i,1) * B2(1,j) &
!               + A2(i,2) * B2(2,j) &
!               + A2(i,3) * B2(3,j) &
!               + A2(i,4) * B2(4,j) &
!               + A2(i,5) * B2(5,j) &
!               + A2(i,6) * B2(6,j) &
!               + A2(i,7) * B2(7,j) &
!               + A2(i,8) * B2(8,j)
!
!      C3(i,j) =  A3(i,1) * B3(1,j) &
!               + A3(i,2) * B3(2,j) &
!               + A3(i,3) * B3(3,j) &
!               + A3(i,4) * B3(4,j) &
!               + A3(i,5) * B3(5,j) &
!               + A3(i,6) * B3(6,j) &
!               + A3(i,7) * B3(7,j) &
!               + A3(i,8) * B3(8,j)
!    enddo
!  enddo
!
!  end subroutine mxm8_3comp
!
!
!!--------------------------------------------------------------------------------------------
!
!  subroutine mxmN_3comp(A1,A2,A3,n1,B1,B2,B3,n2,C1,C2,C3,n3)
!
!  use constants_solver, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n2,n3
!  real(kind=CUSTOM_REAL),dimension(n1,n2) :: A1,A2,A3
!  real(kind=CUSTOM_REAL),dimension(n2,n3) :: B1,B2,B3
!  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C1,C2,C3
!
!  ! local parameters
!  integer :: i,j,k
!  real(kind=CUSTOM_REAL) :: tmp1,tmp2,tmp3
!
!  ! general matrix-matrix multiplication
!  do j = 1,n3
!    do k = 1,n2
!      tmp1 = B1(k,j)
!      tmp2 = B2(k,j)
!      tmp3 = B3(k,j)
!      do i = 1,n1
!        C1(i,j) = C1(i,j) + A1(i,k) * tmp1
!        C2(i,j) = C2(i,j) + A2(i,k) * tmp2
!        C3(i,j) = C3(i,j) + A3(i,k) * tmp3
!      enddo
!    enddo
!  enddo
!
!  end subroutine mxmN_3comp

