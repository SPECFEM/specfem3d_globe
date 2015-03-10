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
! the Free Software Foundation; either version 2 of the License, or
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

  subroutine compute_forces_crust_mantle_Dev( NSPEC,NGLOB,NSPEC_ATT, &
                                              deltat, &
                                              displ_crust_mantle, &
                                              accel_crust_mantle, &
                                              phase_is_inner, &
                                              R_xx,R_yy,R_xy,R_xz,R_yz, &
                                              R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                              epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                              epsilondev_xz,epsilondev_yz, &
                                              epsilon_trace_over_3, &
                                              alphaval,betaval,gammaval, &
                                              factor_common,vnspec )

! this routine is optimized for NGLLX = NGLLY = NGLLZ = 5 using the Deville et al. (2002) inlined matrix-matrix products

  use constants_solver

  use specfem_par,only: &
    hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
    wgll_cube, &
    minus_gravity_table,density_table,minus_deriv_gravity_table, &
    COMPUTE_AND_STORE_STRAIN,USE_LDDRK

  use specfem_par_crustmantle,only: &
    xstore => xstore_crust_mantle,ystore => ystore_crust_mantle,zstore => zstore_crust_mantle, &
    xix => xix_crust_mantle,xiy => xiy_crust_mantle,xiz => xiz_crust_mantle, &
    etax => etax_crust_mantle,etay => etay_crust_mantle,etaz => etaz_crust_mantle, &
    gammax => gammax_crust_mantle,gammay => gammay_crust_mantle,gammaz => gammaz_crust_mantle, &
    kappavstore => kappavstore_crust_mantle,kappahstore => kappahstore_crust_mantle, &
    muvstore => muvstore_crust_mantle,muhstore => muhstore_crust_mantle, &
    eta_anisostore => eta_anisostore_crust_mantle, &
    c11store => c11store_crust_mantle,c12store => c12store_crust_mantle,c13store => c13store_crust_mantle, &
    c14store => c14store_crust_mantle,c15store => c15store_crust_mantle,c16store => c16store_crust_mantle, &
    c22store => c22store_crust_mantle,c23store => c23store_crust_mantle,c24store => c24store_crust_mantle, &
    c25store => c25store_crust_mantle,c26store => c26store_crust_mantle,c33store => c33store_crust_mantle, &
    c34store => c34store_crust_mantle,c35store => c35store_crust_mantle,c36store => c36store_crust_mantle, &
    c44store => c44store_crust_mantle,c45store => c45store_crust_mantle,c46store => c46store_crust_mantle, &
    c55store => c55store_crust_mantle,c56store => c56store_crust_mantle,c66store => c66store_crust_mantle, &
    ibool => ibool_crust_mantle, &
    ispec_is_tiso => ispec_is_tiso_crust_mantle, &
    one_minus_sum_beta => one_minus_sum_beta_crust_mantle, &
    phase_ispec_inner => phase_ispec_inner_crust_mantle, &
    nspec_outer => nspec_outer_crust_mantle, &
    nspec_inner => nspec_inner_crust_mantle

  use specfem_par,only: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D

!daniel: att - debug
!  use specfem_par,only: it,NSTEP

  implicit none

  integer :: NSPEC,NGLOB,NSPEC_ATT

  ! time step
  real(kind=CUSTOM_REAL) deltat

  ! displacement, velocity and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: displ_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: accel_crust_mantle

  ! variable sized array variables
  integer :: vnspec

  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATT) :: R_xx,R_yy,R_xy,R_xz,R_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATT) :: &
    R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: epsilon_trace_over_3

  ! [alpha,beta,gamma]val reduced to N_SLS and factor_common to N_SLS*NUM_NODES
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,N_SLS,vnspec) :: factor_common
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval

  ! inner/outer element run flag
  logical :: phase_is_inner

  ! local parameters

  ! Deville
  ! manually inline the calls to the Deville et al. (2002) routines
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NDIM) :: sum_terms

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5) :: epsilondev_loc
  real(kind=CUSTOM_REAL) fac1,fac2,fac3

  ! for gravity
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NDIM) :: rho_s_H

  integer :: ispec,iglob
  integer :: num_elements,ispec_p
  integer :: iphase

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

! ****************************************************
!   big loop over all spectral elements in the solid
! ****************************************************

!  computed_elements = 0
  if (.not. phase_is_inner) then
    iphase = 1
    num_elements = nspec_outer
  else
    iphase = 2
    num_elements = nspec_inner
  endif

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
!$OMP one_minus_sum_beta,epsilon_trace_over_3,c11store,c12store,c13store,c14store,c15store, &
!$OMP c16store,c22store,c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
!$OMP c36store,c44store,c45store,c46store,c55store,c56store,c66store,ispec_is_tiso, &
!$OMP kappavstore,muvstore,kappahstore,muhstore,eta_anisostore,ibool,ystore,zstore, &
!$OMP R_xx,R_yy,R_xy,R_xz,R_yz, &
!$OMP xstore,minus_gravity_table,minus_deriv_gravity_table,density_table, &
!$OMP displ_crust_mantle,wgll_cube,hprime_xxt,hprime_xx, &
!$OMP vnspec, &
!$OMP accel_crust_mantle, &
!$OMP hprimewgll_xx,hprimewgll_xxt, &
!$OMP alphaval,betaval, &
!$OMP epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
!$OMP gammaval,factor_common, &
!$OMP iphase, &
!$OMP phase_ispec_inner, &
!$OMP num_elements, USE_LDDRK, &
!$OMP wgllwgll_xy_3D, wgllwgll_xz_3D, wgllwgll_yz_3D, &
!$OMP R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
!$OMP deltat, COMPUTE_AND_STORE_STRAIN ) &
!$OMP PRIVATE(ispec,fac1,fac2,fac3,sum_terms,ispec_p, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk, &
#else
!$OMP i,j,k, &
#endif
!$OMP tempx1,tempx2,tempx3, &
!$OMP newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3, &
!$OMP dummyx_loc,dummyy_loc,dummyz_loc,rho_s_H,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
!$OMP iglob,epsilondev_loc)

!$OMP DO SCHEDULE(GUIDED)
  do ispec_p = 1,num_elements

    ispec = phase_ispec_inner(ispec_p,iphase)

    ! only compute element which belong to current phase (inner or outer elements)

    DO_LOOP_IJK

      iglob = ibool(INDEX_IJK,ispec)
      dummyx_loc(INDEX_IJK) = displ_crust_mantle(1,iglob)
      dummyy_loc(INDEX_IJK) = displ_crust_mantle(2,iglob)
      dummyz_loc(INDEX_IJK) = displ_crust_mantle(3,iglob)

    ENDDO_LOOP_IJK

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
    if (ANISOTROPIC_3D_MANTLE_VAL) then
       ! anisotropic element
       call compute_element_aniso(ispec, &
                                  minus_gravity_table,density_table,minus_deriv_gravity_table, &
                                  xstore,ystore,zstore, &
                                  xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                                  wgll_cube, &
                                  c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                                  c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                                  c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                                  ibool, &
                                  R_xx,R_yy,R_xy,R_xz,R_yz, &
                                  epsilon_trace_over_3, &
                                  one_minus_sum_beta,vnspec, &
                                  tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                  dummyx_loc,dummyy_loc,dummyz_loc, &
                                  epsilondev_loc,rho_s_H)
    else
       if (.not. ispec_is_tiso(ispec)) then
          ! isotropic element
          call compute_element_iso(ispec, &
                                   minus_gravity_table,density_table,minus_deriv_gravity_table, &
                                   xstore,ystore,zstore, &
                                   xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                                   wgll_cube, &
                                   kappavstore,muvstore, &
                                   ibool, &
                                   R_xx,R_yy,R_xy,R_xz,R_yz, &
                                   epsilon_trace_over_3, &
                                   one_minus_sum_beta,vnspec, &
                                   tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                   dummyx_loc,dummyy_loc,dummyz_loc, &
                                   epsilondev_loc,rho_s_H)
       else
          ! transverse isotropic element
          call compute_element_tiso(ispec, &
                                     minus_gravity_table,density_table,minus_deriv_gravity_table, &
                                     xstore,ystore,zstore, &
                                     xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                                     wgll_cube, &
                                     kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                                     ibool, &
                                     R_xx,R_yy,R_xy,R_xz,R_yz, &
                                     epsilon_trace_over_3, &
                                     one_minus_sum_beta,vnspec, &
                                     tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                     dummyx_loc,dummyy_loc,dummyz_loc, &
                                     epsilondev_loc,rho_s_H)
       endif ! .not. ispec_is_tiso
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
      sum_terms(INDEX_IJK,1) = - (fac1*newtempx1(INDEX_IJK) + fac2*newtempx2(INDEX_IJK) + fac3*newtempx3(INDEX_IJK))
      sum_terms(INDEX_IJK,2) = - (fac1*newtempy1(INDEX_IJK) + fac2*newtempy2(INDEX_IJK) + fac3*newtempy3(INDEX_IJK))
      sum_terms(INDEX_IJK,3) = - (fac1*newtempz1(INDEX_IJK) + fac2*newtempz2(INDEX_IJK) + fac3*newtempz3(INDEX_IJK))

    ENDDO_LOOP_IJK

    ! adds gravity terms
    if (GRAVITY_VAL) then

#ifdef FORCE_VECTORIZATION
      do ijk = 1,NDIM*NGLLCUBE
        sum_terms(ijk,1,1,1) = sum_terms(ijk,1,1,1) + rho_s_H(ijk,1,1,1)
      enddo
#else
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            sum_terms(INDEX_IJK,1) = sum_terms(INDEX_IJK,1) + rho_s_H(INDEX_IJK,1)
            sum_terms(INDEX_IJK,2) = sum_terms(INDEX_IJK,2) + rho_s_H(INDEX_IJK,2)
            sum_terms(INDEX_IJK,3) = sum_terms(INDEX_IJK,3) + rho_s_H(INDEX_IJK,3)
          enddo
        enddo
      enddo
#endif

    endif


    ! updates acceleration

#ifdef FORCE_VECTORIZATION
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
    do ijk = 1,NGLLCUBE
#else
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
#endif

          iglob = ibool(INDEX_IJK,ispec)

          ! do NOT use array syntax ":" for the three statements below otherwise most compilers
          ! will not be able to vectorize the outer loop
#ifdef USE_OPENMP_ATOMIC_INSTEAD_OF_CRITICAL
!$OMP ATOMIC
#endif
          accel_crust_mantle(1,iglob) = accel_crust_mantle(1,iglob) + sum_terms(INDEX_IJK,1)
#ifdef USE_OPENMP_ATOMIC_INSTEAD_OF_CRITICAL
!$OMP ATOMIC
#endif
          accel_crust_mantle(2,iglob) = accel_crust_mantle(2,iglob) + sum_terms(INDEX_IJK,2)
#ifdef USE_OPENMP_ATOMIC_INSTEAD_OF_CRITICAL
!$OMP ATOMIC
#endif
          accel_crust_mantle(3,iglob) = accel_crust_mantle(3,iglob) + sum_terms(INDEX_IJK,3)

#ifdef FORCE_VECTORIZATION
    enddo
#ifndef USE_OPENMP_ATOMIC_INSTEAD_OF_CRITICAL
!$OMP END CRITICAL
#endif
#else
        enddo
      enddo
    enddo
#endif
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
                                                 c44store,muvstore, &
                                                 epsilondev_loc, &
                                                 deltat)
      else
        call compute_element_att_memory_cm(ispec,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                           ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec,factor_common, &
                                           alphaval,betaval,gammaval, &
                                           c44store,muvstore, &
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
!$OMP enddo
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

! 3 different arrays for x/y/z-components, 2-dimensional arrays (25,5)/(5,25), same B matrix for all 3 component arrays

  use constants_solver,only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5) :: A
  real(kind=CUSTOM_REAL),dimension(5,n3) :: B1,B2,B3
  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C1,C2,C3

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
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

! 3 different arrays for x/y/z-components, 2-dimensional arrays (25,5)/(5,25), same B matrix for all 3 component arrays

  use constants_solver,only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(5,n3) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C1,C2,C3

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
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

! 3 different arrays for x/y/z-components, 3-dimensional arrays (5,5,5), same B matrix for all 3 component arrays

  use constants_solver,only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,5,n3) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(5,n2) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3) :: C1,C2,C3

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do j = 1,n2
    do i = 1,n1
      ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
      do k = 1,n3
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

  end subroutine compute_forces_crust_mantle_Dev



