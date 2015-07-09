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

  subroutine compute_forces_inner_core_Dev( NSPEC,NGLOB,NSPEC_ATT, &
                                            deltat, &
                                            displ_inner_core, &
                                            accel_inner_core, &
                                            phase_is_inner, &
                                            R_xx,R_yy,R_xy,R_xz,R_yz, &
                                            R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                            epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                            epsilondev_xz,epsilondev_yz, &
                                            epsilon_trace_over_3,&
                                            alphaval,betaval,gammaval, &
                                            factor_common,vnspec)

! this routine is optimized for NGLLX = NGLLY = NGLLZ = 5 using the Deville et al. (2002) inlined matrix-matrix products

  use constants_solver

  use specfem_par,only: &
    hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
    wgll_cube, &
    minus_gravity_table,density_table,minus_deriv_gravity_table, &
    COMPUTE_AND_STORE_STRAIN,USE_LDDRK

  use specfem_par_innercore,only: &
    xstore => xstore_inner_core,ystore => ystore_inner_core,zstore => zstore_inner_core, &
    xix => xix_inner_core,xiy => xiy_inner_core,xiz => xiz_inner_core, &
    etax => etax_inner_core,etay => etay_inner_core,etaz => etaz_inner_core, &
    gammax => gammax_inner_core,gammay => gammay_inner_core,gammaz => gammaz_inner_core, &
    kappavstore => kappavstore_inner_core, &
    muvstore => muvstore_inner_core, &
    c11store => c11store_inner_core,c12store => c12store_inner_core,c13store => c13store_inner_core, &
    c33store => c33store_inner_core,c44store => c44store_inner_core, &
    ibool => ibool_inner_core,idoubling => idoubling_inner_core, &
    one_minus_sum_beta => one_minus_sum_beta_inner_core, &
    phase_ispec_inner => phase_ispec_inner_inner_core, &
    nspec_outer => nspec_outer_inner_core, &
    nspec_inner => nspec_inner_inner_core

  use specfem_par,only: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D

  implicit none

  integer :: NSPEC,NGLOB,NSPEC_ATT

  ! time step
  real(kind=CUSTOM_REAL) deltat

  ! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: displ_inner_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: accel_inner_core

  ! for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  ! variable lengths for factor_common and one_minus_sum_beta
  integer :: vnspec
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,N_SLS,vnspec) :: factor_common
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATT) :: R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATT) :: &
    R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: epsilon_trace_over_3

  ! inner/outer element run flag
  logical :: phase_is_inner

  ! local parameters
  ! Deville
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NDIM) :: sum_terms
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5) :: epsilondev_loc

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz

  real(kind=CUSTOM_REAL) fac1,fac2,fac3
  real(kind=CUSTOM_REAL) lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) kappal
  real(kind=CUSTOM_REAL) templ

  real(kind=CUSTOM_REAL) minus_sum_beta
  real(kind=CUSTOM_REAL) c11l,c33l,c12l,c13l,c44l

  ! for gravity
  double precision radius,rho,minus_g,minus_dg
  double precision minus_g_over_radius,minus_dg_plus_g_over_radius
  double precision cos_theta,sin_theta,cos_phi,sin_phi
  double precision cos_theta_sq,sin_theta_sq,cos_phi_sq,sin_phi_sq
  double precision theta,phi,factor,gxl,gyl,gzl,sx_l,sy_l,sz_l
  double precision Hxxl,Hyyl,Hzzl,Hxyl,Hxzl,Hyzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NDIM) :: rho_s_H
  real(kind=CUSTOM_REAL) sigma_yx,sigma_zx,sigma_zy

  integer :: int_radius
  integer :: ispec,iglob
  integer :: num_elements,ispec_p
  integer :: iphase

#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
  real(kind=CUSTOM_REAL) :: R_xx_val,R_yy_val
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

!$OMP PARALLEL DEFAULT( NONE ) &
!$OMP SHARED( num_elements, phase_ispec_inner, iphase, idoubling, ibool, displ_inner_core, hprime_xx, hprime_xxT, &
!$OMP xix, xiy, xiz, etax, etay, etaz, gammax, gammay, gammaz, COMPUTE_AND_STORE_STRAIN, c11store,  c12store, c13store, &
!$OMP c33store, c44store, one_minus_sum_beta, muvstore, kappavstore, R_xx, R_yy, R_xy, R_xz, R_yz, xstore, ystore, zstore, &
!$OMP minus_gravity_table, minus_deriv_gravity_table, density_table, wgll_cube, hprimewgll_xxT, hprimewgll_xx, wgllwgll_yz_3D, &
!$OMP wgllwgll_xz_3D, wgllwgll_xy_3D, accel_inner_core, USE_LDDRK, R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
!$OMP vnspec, factor_common, deltat, alphaval,betaval,gammaval, epsilondev_xx,epsilondev_yy,epsilondev_xy, &
!$OMP epsilondev_xz,epsilondev_yz, epsilon_trace_over_3 ) &
!$OMP PRIVATE( ispec_p, ispec, iglob, dummyx_loc, dummyy_loc, dummyz_loc, tempx2, tempy2, tempz2, xixl, xiyl, xizl, &
!$OMP etaxl, etayl, etazl, gammaxl, gammayl, gammazl, jacobianl, duxdxl, tempx1, tempx3, duxdyl, duxdzl, duydxl, &
!$OMP tempy1, tempy3, duydyl, duydzl, tempz1, tempz3, duzdxl, duzdyl, duzdzl, duxdxl_plus_duydyl, duxdxl_plus_duzdzl, &
!$OMP duydyl_plus_duzdzl, duxdyl_plus_duydxl, duzdxl_plus_duxdzl, duzdyl_plus_duydzl, templ, epsilondev_loc, &
!$OMP c11l, c12l, c13l, c33l, c44l, minus_sum_beta, mul, sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz, &
!$OMP kappal, lambdalplus2mul, lambdal, sigma_yx, sigma_zx, sigma_zy, radius, theta, phi, cos_theta, sin_theta, cos_phi, &
!$OMP sin_phi, cos_theta_sq, sin_theta_sq, cos_phi_sq, sin_phi_sq, int_radius, minus_g, rho, gxl, gyl, gzl, minus_dg, &
!$OMP minus_g_over_radius, minus_dg_plus_g_over_radius, Hxxl, Hyyl, Hzzl, Hxyl, Hxzl, Hyzl, sx_l, sy_l, sz_l, &
!$OMP factor, rho_s_H, newtempx2, newtempy2, newtempz2, fac1, fac2, fac3, sum_terms, newtempx1, newtempx3 , newtempy1, &
#ifdef FORCE_VECTORIZATION
!$OMP R_xx_val, R_yy_val, &
#endif
!$OMP newtempy3, newtempz1, newtempz3)

!$OMP DO SCHEDULE(GUIDED)
  do ispec_p = 1,num_elements

    ispec = phase_ispec_inner(ispec_p,iphase)

    ! only compute element which belong to current phase (inner or outer elements)

    ! exclude fictitious elements in central cube
    if (idoubling(ispec) /= IFLAG_IN_FICTITIOUS_CUBE) then

      DO_LOOP_IJK

        iglob = ibool(INDEX_IJK,ispec)
        dummyx_loc(INDEX_IJK) = displ_inner_core(1,iglob)
        dummyy_loc(INDEX_IJK) = displ_inner_core(2,iglob)
        dummyz_loc(INDEX_IJK) = displ_inner_core(3,iglob)

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

!      do j = 1,m2
!         do i = 1,m1
!            C1_m1_m2_5points(i,j) = hprime_xx(i,1)*B1_m1_m2_5points(1,j) + &
!                 hprime_xx(i,2)*B1_m1_m2_5points(2,j) + &
!                 hprime_xx(i,3)*B1_m1_m2_5points(3,j) + &
!                 hprime_xx(i,4)*B1_m1_m2_5points(4,j) + &
!                 hprime_xx(i,5)*B1_m1_m2_5points(5,j)
!
!            C2_m1_m2_5points(i,j) = hprime_xx(i,1)*B2_m1_m2_5points(1,j) + &
!                 hprime_xx(i,2)*B2_m1_m2_5points(2,j) + &
!                 hprime_xx(i,3)*B2_m1_m2_5points(3,j) + &
!                 hprime_xx(i,4)*B2_m1_m2_5points(4,j) + &
!                 hprime_xx(i,5)*B2_m1_m2_5points(5,j)
!
!            C3_m1_m2_5points(i,j) = hprime_xx(i,1)*B3_m1_m2_5points(1,j) + &
!                 hprime_xx(i,2)*B3_m1_m2_5points(2,j) + &
!                 hprime_xx(i,3)*B3_m1_m2_5points(3,j) + &
!                 hprime_xx(i,4)*B3_m1_m2_5points(4,j) + &
!                 hprime_xx(i,5)*B3_m1_m2_5points(5,j)
!         enddo
!      enddo
!
!      do j = 1,m1
!         do i = 1,m1
!            ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
!            do k = 1,NGLLX
!               tempx2(i,j,k) = dummyx_loc(i,1,k)*hprime_xxT(1,j) + &
!                    dummyx_loc(i,2,k)*hprime_xxT(2,j) + &
!                    dummyx_loc(i,3,k)*hprime_xxT(3,j) + &
!                    dummyx_loc(i,4,k)*hprime_xxT(4,j) + &
!                    dummyx_loc(i,5,k)*hprime_xxT(5,j)
!
!               tempy2(i,j,k) = dummyy_loc(i,1,k)*hprime_xxT(1,j) + &
!                    dummyy_loc(i,2,k)*hprime_xxT(2,j) + &
!                    dummyy_loc(i,3,k)*hprime_xxT(3,j) + &
!                    dummyy_loc(i,4,k)*hprime_xxT(4,j) + &
!                    dummyy_loc(i,5,k)*hprime_xxT(5,j)
!
!               tempz2(i,j,k) = dummyz_loc(i,1,k)*hprime_xxT(1,j) + &
!                    dummyz_loc(i,2,k)*hprime_xxT(2,j) + &
!                    dummyz_loc(i,3,k)*hprime_xxT(3,j) + &
!                    dummyz_loc(i,4,k)*hprime_xxT(4,j) + &
!                    dummyz_loc(i,5,k)*hprime_xxT(5,j)
!            enddo
!         enddo
!      enddo
!
!      do j = 1,m1
!         do i = 1,m2
!            C1_mxm_m2_m1_5points(i,j) = A1_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
!                 A1_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
!                 A1_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
!                 A1_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
!                 A1_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)
!
!            C2_mxm_m2_m1_5points(i,j) = A2_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
!                 A2_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
!                 A2_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
!                 A2_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
!                 A2_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)
!
!            C3_mxm_m2_m1_5points(i,j) = A3_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
!                 A3_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
!                 A3_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
!                 A3_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
!                 A3_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)
!         enddo
!      enddo
!
      DO_LOOP_IJK

        ! get derivatives of ux, uy and uz with respect to x, y and z
        xixl = xix(INDEX_IJK,ispec)
        xiyl = xiy(INDEX_IJK,ispec)
        xizl = xiz(INDEX_IJK,ispec)
        etaxl = etax(INDEX_IJK,ispec)
        etayl = etay(INDEX_IJK,ispec)
        etazl = etaz(INDEX_IJK,ispec)
        gammaxl = gammax(INDEX_IJK,ispec)
        gammayl = gammay(INDEX_IJK,ispec)
        gammazl = gammaz(INDEX_IJK,ispec)

        ! compute the Jacobian
        jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                      - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                      + xizl*(etaxl*gammayl-etayl*gammaxl))

        duxdxl = xixl*tempx1(INDEX_IJK) + etaxl*tempx2(INDEX_IJK) + gammaxl*tempx3(INDEX_IJK)
        duxdyl = xiyl*tempx1(INDEX_IJK) + etayl*tempx2(INDEX_IJK) + gammayl*tempx3(INDEX_IJK)
        duxdzl = xizl*tempx1(INDEX_IJK) + etazl*tempx2(INDEX_IJK) + gammazl*tempx3(INDEX_IJK)

        duydxl = xixl*tempy1(INDEX_IJK) + etaxl*tempy2(INDEX_IJK) + gammaxl*tempy3(INDEX_IJK)
        duydyl = xiyl*tempy1(INDEX_IJK) + etayl*tempy2(INDEX_IJK) + gammayl*tempy3(INDEX_IJK)
        duydzl = xizl*tempy1(INDEX_IJK) + etazl*tempy2(INDEX_IJK) + gammazl*tempy3(INDEX_IJK)

        duzdxl = xixl*tempz1(INDEX_IJK) + etaxl*tempz2(INDEX_IJK) + gammaxl*tempz3(INDEX_IJK)
        duzdyl = xiyl*tempz1(INDEX_IJK) + etayl*tempz2(INDEX_IJK) + gammayl*tempz3(INDEX_IJK)
        duzdzl = xizl*tempz1(INDEX_IJK) + etazl*tempz2(INDEX_IJK) + gammazl*tempz3(INDEX_IJK)

        ! precompute some sums to save CPU time
        duxdxl_plus_duydyl = duxdxl + duydyl
        duxdxl_plus_duzdzl = duxdxl + duzdzl
        duydyl_plus_duzdzl = duydyl + duzdzl
        duxdyl_plus_duydxl = duxdyl + duydxl
        duzdxl_plus_duxdzl = duzdxl + duxdzl
        duzdyl_plus_duydzl = duzdyl + duydzl

        ! compute deviatoric strain
        if (COMPUTE_AND_STORE_STRAIN) then
          templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
          if (NSPEC_INNER_CORE_STRAIN_ONLY == 1) then
            if (ispec == 1) then
              epsilon_trace_over_3(INDEX_IJK,1) = templ
            endif
          else
            epsilon_trace_over_3(INDEX_IJK,ispec) = templ
          endif
          epsilondev_loc(INDEX_IJK,1) = duxdxl - templ
          epsilondev_loc(INDEX_IJK,2) = duydyl - templ
          epsilondev_loc(INDEX_IJK,3) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
          epsilondev_loc(INDEX_IJK,4) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
          epsilondev_loc(INDEX_IJK,5) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl
        endif

        if (ANISOTROPIC_INNER_CORE_VAL) then
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
          c11l = c11store(INDEX_IJK,ispec)
          c12l = c12store(INDEX_IJK,ispec)
          c13l = c13store(INDEX_IJK,ispec)
          c33l = c33store(INDEX_IJK,ispec)
          c44l = c44store(INDEX_IJK,ispec)

          ! use unrelaxed parameters if attenuation
          if (ATTENUATION_VAL) then
            if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
              minus_sum_beta =  one_minus_sum_beta(INDEX_IJK,ispec) - 1.0_CUSTOM_REAL
            else
              minus_sum_beta =  one_minus_sum_beta(1,1,1,ispec) - 1.0_CUSTOM_REAL
            endif
            mul = muvstore(INDEX_IJK,ispec)
            c11l = c11l + FOUR_THIRDS * minus_sum_beta * mul
            c12l = c12l - TWO_THIRDS * minus_sum_beta * mul
            c13l = c13l - TWO_THIRDS * minus_sum_beta * mul
            c33l = c33l + FOUR_THIRDS * minus_sum_beta * mul
            c44l = c44l + minus_sum_beta * mul
          endif

          sigma_xx = c11l*duxdxl + c12l*duydyl + c13l*duzdzl
          sigma_yy = c12l*duxdxl + c11l*duydyl + c13l*duzdzl
          sigma_zz = c13l*duxdxl + c13l*duydyl + c33l*duzdzl
          sigma_xy = 0.5_CUSTOM_REAL*(c11l-c12l)*duxdyl_plus_duydxl
          sigma_xz = c44l*duzdxl_plus_duxdzl
          sigma_yz = c44l*duzdyl_plus_duydzl
        else

          ! inner core with no anisotropy, use kappav and muv for instance
          ! layer with no anisotropy, use kappav and muv for instance
          kappal = kappavstore(INDEX_IJK,ispec)
          mul = muvstore(INDEX_IJK,ispec)

          ! use unrelaxed parameters if attenuation
          if (ATTENUATION_VAL) then
            if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
              mul = mul * one_minus_sum_beta(INDEX_IJK,ispec)
            else
              mul = mul * one_minus_sum_beta(1,1,1,ispec)
            endif
          endif

          lambdalplus2mul = kappal + FOUR_THIRDS * mul
          lambdal = lambdalplus2mul - 2._CUSTOM_REAL*mul

          ! compute stress sigma
          sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
          sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
          sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl

          sigma_xy = mul*duxdyl_plus_duydxl
          sigma_xz = mul*duzdxl_plus_duxdzl
          sigma_yz = mul*duzdyl_plus_duydzl

        endif

        ! subtract memory variables if attenuation
        if (ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL) then
#ifdef FORCE_VECTORIZATION
          ! do NOT put this is a subroutine, otherwise the call to the subroutine prevents compilers
          !from vectorizing the outer loop

          ! here we assume that N_SLS == 3 in order to be able to unroll and suppress the loop
          ! in order to vectorize the outer loop
          R_xx_val = R_xx(INDEX_IJK,1,ispec)
          R_yy_val = R_yy(INDEX_IJK,1,ispec)
          sigma_xx = sigma_xx - R_xx_val
          sigma_yy = sigma_yy - R_yy_val
          sigma_zz = sigma_zz + R_xx_val + R_yy_val
          sigma_xy = sigma_xy - R_xy(INDEX_IJK,1,ispec)
          sigma_xz = sigma_xz - R_xz(INDEX_IJK,1,ispec)
          sigma_yz = sigma_yz - R_yz(INDEX_IJK,1,ispec)

          R_xx_val = R_xx(INDEX_IJK,2,ispec)
          R_yy_val = R_yy(INDEX_IJK,2,ispec)
          sigma_xx = sigma_xx - R_xx_val
          sigma_yy = sigma_yy - R_yy_val
          sigma_zz = sigma_zz + R_xx_val + R_yy_val
          sigma_xy = sigma_xy - R_xy(INDEX_IJK,2,ispec)
          sigma_xz = sigma_xz - R_xz(INDEX_IJK,2,ispec)
          sigma_yz = sigma_yz - R_yz(INDEX_IJK,2,ispec)

          R_xx_val = R_xx(INDEX_IJK,3,ispec)
          R_yy_val = R_yy(INDEX_IJK,3,ispec)
          sigma_xx = sigma_xx - R_xx_val
          sigma_yy = sigma_yy - R_yy_val
          sigma_zz = sigma_zz + R_xx_val + R_yy_val
          sigma_xy = sigma_xy - R_xy(INDEX_IJK,3,ispec)
          sigma_xz = sigma_xz - R_xz(INDEX_IJK,3,ispec)
          sigma_yz = sigma_yz - R_yz(INDEX_IJK,3,ispec)
#else
          ! note: Fortran passes pointers to array location, thus R_memory(1,1,...) is fine
          call compute_element_att_stress(i,j,k,R_xx(1,1,1,1,ispec),R_yy(1,1,1,1,ispec),R_xy(1,1,1,1,ispec), &
                                          R_xz(1,1,1,1,ispec),R_yz(1,1,1,1,ispec), &
                                          sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz)
#endif
        endif

        ! define symmetric components of sigma for gravity
        sigma_yx = sigma_xy
        sigma_zx = sigma_xz
        sigma_zy = sigma_yz

        ! compute non-symmetric terms for gravity
        if (GRAVITY_VAL) then

          ! use mesh coordinates to get theta and phi
          ! x y and z contain r theta and phi
          iglob = ibool(INDEX_IJK,ispec)
          radius = dble(xstore(iglob))
          theta = dble(ystore(iglob))
          phi = dble(zstore(iglob))

          ! make sure radius is never zero even for points at center of cube
          ! because we later divide by radius
          if (radius < 100.d0 / R_EARTH) radius = 100.d0 / R_EARTH

          cos_theta = dcos(theta)
          sin_theta = dsin(theta)
          cos_phi = dcos(phi)
          sin_phi = dsin(phi)

          cos_theta_sq = cos_theta**2
          sin_theta_sq = sin_theta**2
          cos_phi_sq = cos_phi**2
          sin_phi_sq = sin_phi**2

          ! get g, rho and dg/dr=dg
          ! spherical components of the gravitational acceleration
          ! for efficiency replace with lookup table every 100 m in radial direction
          ! make sure we never use zero for point exactly at the center of the Earth
          int_radius = max(1,nint(radius * R_EARTH_KM * 10.d0))
          minus_g = minus_gravity_table(int_radius)
          minus_dg = minus_deriv_gravity_table(int_radius)
          rho = density_table(int_radius)

          ! Cartesian components of the gravitational acceleration
          gxl = minus_g*sin_theta*cos_phi
          gyl = minus_g*sin_theta*sin_phi
          gzl = minus_g*cos_theta

          ! Cartesian components of gradient of gravitational acceleration
          ! obtained from spherical components
          minus_g_over_radius = minus_g / radius
          minus_dg_plus_g_over_radius = minus_dg - minus_g_over_radius

          Hxxl = minus_g_over_radius*(cos_phi_sq*cos_theta_sq + sin_phi_sq) + cos_phi_sq*minus_dg*sin_theta_sq
          Hyyl = minus_g_over_radius*(cos_phi_sq + cos_theta_sq*sin_phi_sq) + minus_dg*sin_phi_sq*sin_theta_sq
          Hzzl = cos_theta_sq*minus_dg + minus_g_over_radius*sin_theta_sq
          Hxyl = cos_phi*minus_dg_plus_g_over_radius*sin_phi*sin_theta_sq
          Hxzl = cos_phi*cos_theta*minus_dg_plus_g_over_radius*sin_theta
          Hyzl = cos_theta*minus_dg_plus_g_over_radius*sin_phi*sin_theta

          ! for locality principle, we set iglob again, in order to have it in the cache again
          iglob = ibool(INDEX_IJK,ispec)

          ! get displacement and multiply by density to compute G tensor
          sx_l = rho * dble(displ_inner_core(1,iglob))
          sy_l = rho * dble(displ_inner_core(2,iglob))
          sz_l = rho * dble(displ_inner_core(3,iglob))

          ! compute G tensor from s . g and add to sigma (not symmetric)
          sigma_xx = sigma_xx + real(sy_l*gyl + sz_l*gzl, kind=CUSTOM_REAL)
          sigma_yy = sigma_yy + real(sx_l*gxl + sz_l*gzl, kind=CUSTOM_REAL)
          sigma_zz = sigma_zz + real(sx_l*gxl + sy_l*gyl, kind=CUSTOM_REAL)

          sigma_xy = sigma_xy - real(sx_l * gyl, kind=CUSTOM_REAL)
          sigma_yx = sigma_yx - real(sy_l * gxl, kind=CUSTOM_REAL)

          sigma_xz = sigma_xz - real(sx_l * gzl, kind=CUSTOM_REAL)
          sigma_zx = sigma_zx - real(sz_l * gxl, kind=CUSTOM_REAL)

          sigma_yz = sigma_yz - real(sy_l * gzl, kind=CUSTOM_REAL)
          sigma_zy = sigma_zy - real(sz_l * gyl, kind=CUSTOM_REAL)

          ! precompute vector
          factor = dble(jacobianl) * wgll_cube(INDEX_IJK)
          rho_s_H(INDEX_IJK,1) = real(factor * (sx_l * Hxxl + sy_l * Hxyl + sz_l * Hxzl), kind=CUSTOM_REAL)
          rho_s_H(INDEX_IJK,2) = real(factor * (sx_l * Hxyl + sy_l * Hyyl + sz_l * Hyzl), kind=CUSTOM_REAL)
          rho_s_H(INDEX_IJK,3) = real(factor * (sx_l * Hxzl + sy_l * Hyzl + sz_l * Hzzl), kind=CUSTOM_REAL)

        endif  ! end of section with gravity terms

        ! form dot product with test vector, non-symmetric form
        tempx1(INDEX_IJK) = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl) ! this goes to accel_x
        tempy1(INDEX_IJK) = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl) ! this goes to accel_y
        tempz1(INDEX_IJK) = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl) ! this goes to accel_z

        tempx2(INDEX_IJK) = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl) ! this goes to accel_x
        tempy2(INDEX_IJK) = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl) ! this goes to accel_y
        tempz2(INDEX_IJK) = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl) ! this goes to accel_z

        tempx3(INDEX_IJK) = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl) ! this goes to accel_x
        tempy3(INDEX_IJK) = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl) ! this goes to accel_y
        tempz3(INDEX_IJK) = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl) ! this goes to accel_z

      ENDDO_LOOP_IJK

      ! subroutines adapted from Deville, Fischer and Mund, High-order methods
      ! for incompressible fluid flow, Cambridge University Press (2002),
      ! pages 386 and 389 and Figure 8.3.1

      ! computes 1. matrix multiplication for newtempx1,..
      call mxm5_3comp_singleA(hprimewgll_xxT,m1,tempx1,tempy1,tempz1,newtempx1,newtempy1,newtempz1,m2)
      ! computes 2. matrix multiplication for tempx2,..
      call mxm5_3comp_3dmat_singleB(tempx2,tempy2,tempz2,m1,hprimewgll_xx,m1,newtempx2,newtempy2,newtempz2,NGLLX)
      ! computes 3. matrix multiplication for newtempx3,..
      call mxm5_3comp_singleB(tempx3,tempy3,tempz3,m2,hprimewgll_xx,newtempx3,newtempy3,newtempz3,m1)

!      do j = 1,m2
!        do i = 1,m1
!          E1_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C1_m1_m2_5points(1,j) + &
!                                hprimewgll_xxT(i,2)*C1_m1_m2_5points(2,j) + &
!                                hprimewgll_xxT(i,3)*C1_m1_m2_5points(3,j) + &
!                                hprimewgll_xxT(i,4)*C1_m1_m2_5points(4,j) + &
!                                hprimewgll_xxT(i,5)*C1_m1_m2_5points(5,j)
!
!          E2_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C2_m1_m2_5points(1,j) + &
!                                hprimewgll_xxT(i,2)*C2_m1_m2_5points(2,j) + &
!                                hprimewgll_xxT(i,3)*C2_m1_m2_5points(3,j) + &
!                                hprimewgll_xxT(i,4)*C2_m1_m2_5points(4,j) + &
!                                hprimewgll_xxT(i,5)*C2_m1_m2_5points(5,j)
!
!          E3_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C3_m1_m2_5points(1,j) + &
!                                hprimewgll_xxT(i,2)*C3_m1_m2_5points(2,j) + &
!                                hprimewgll_xxT(i,3)*C3_m1_m2_5points(3,j) + &
!                                hprimewgll_xxT(i,4)*C3_m1_m2_5points(4,j) + &
!                                hprimewgll_xxT(i,5)*C3_m1_m2_5points(5,j)
!        enddo
!      enddo
!
!      do i = 1,m1
!        do j = 1,m1
!          ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
!          do k = 1,NGLLX
!            newtempx2(i,j,k) = tempx2(i,1,k)*hprimewgll_xx(1,j) + &
!                             tempx2(i,2,k)*hprimewgll_xx(2,j) + &
!                             tempx2(i,3,k)*hprimewgll_xx(3,j) + &
!                             tempx2(i,4,k)*hprimewgll_xx(4,j) + &
!                             tempx2(i,5,k)*hprimewgll_xx(5,j)
!
!            newtempy2(i,j,k) = tempy2(i,1,k)*hprimewgll_xx(1,j) + &
!                             tempy2(i,2,k)*hprimewgll_xx(2,j) + &
!                             tempy2(i,3,k)*hprimewgll_xx(3,j) + &
!                             tempy2(i,4,k)*hprimewgll_xx(4,j) + &
!                             tempy2(i,5,k)*hprimewgll_xx(5,j)
!
!            newtempz2(i,j,k) = tempz2(i,1,k)*hprimewgll_xx(1,j) + &
!                             tempz2(i,2,k)*hprimewgll_xx(2,j) + &
!                             tempz2(i,3,k)*hprimewgll_xx(3,j) + &
!                             tempz2(i,4,k)*hprimewgll_xx(4,j) + &
!                             tempz2(i,5,k)*hprimewgll_xx(5,j)
!          enddo
!        enddo
!      enddo
!
!      do j = 1,m1
!        do i = 1,m2
!          E1_mxm_m2_m1_5points(i,j) = C1_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
!                                    C1_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
!                                    C1_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
!                                    C1_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
!                                    C1_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)
!
!          E2_mxm_m2_m1_5points(i,j) = C2_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
!                                    C2_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
!                                    C2_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
!                                    C2_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
!                                    C2_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)
!
!          E3_mxm_m2_m1_5points(i,j) = C3_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
!                                    C3_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
!                                    C3_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
!                                    C3_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
!                                    C3_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)
!        enddo
!      enddo
!
      ! sums contributions

      DO_LOOP_IJK

        fac1 = wgllwgll_yz_3D(INDEX_IJK)
        fac2 = wgllwgll_xz_3D(INDEX_IJK)
        fac3 = wgllwgll_xy_3D(INDEX_IJK)
        sum_terms(INDEX_IJK,1) = - (fac1*newtempx1(INDEX_IJK) + fac2*newtempx2(INDEX_IJK) + fac3*newtempx3(INDEX_IJK))
        sum_terms(INDEX_IJK,2) = - (fac1*newtempy1(INDEX_IJK) + fac2*newtempy2(INDEX_IJK) + fac3*newtempy3(INDEX_IJK))
        sum_terms(INDEX_IJK,3) = - (fac1*newtempz1(INDEX_IJK) + fac2*newtempz2(INDEX_IJK) + fac3*newtempz3(INDEX_IJK))

      ENDDO_LOOP_IJK

      ! adds gravity
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


      ! sum contributions from each element to the global mesh and add gravity terms
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
            ! do NOT use array syntax ":" for the three statements below
            ! otherwise most compilers will not be able to vectorize the outer loop
#ifdef USE_OPENMP_ATOMIC_INSTEAD_OF_CRITICAL
!$OMP ATOMIC
#endif
            accel_inner_core(1,iglob) = accel_inner_core(1,iglob) + sum_terms(INDEX_IJK,1)
#ifdef USE_OPENMP_ATOMIC_INSTEAD_OF_CRITICAL
!$OMP ATOMIC
#endif
            accel_inner_core(2,iglob) = accel_inner_core(2,iglob) + sum_terms(INDEX_IJK,2)
#ifdef USE_OPENMP_ATOMIC_INSTEAD_OF_CRITICAL
!$OMP ATOMIC
#endif
            accel_inner_core(3,iglob) = accel_inner_core(3,iglob) + sum_terms(INDEX_IJK,3)

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

    endif ! end of test to exclude fictitious elements in central cube

  enddo ! of spectral element loop
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

  end subroutine compute_forces_inner_core_Dev

