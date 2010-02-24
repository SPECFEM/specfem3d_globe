!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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


  subroutine compute_forces_crust_mantle_Dev(minus_gravity_table,density_table,minus_deriv_gravity_table, &
          displ,accel,xstore,ystore,zstore, &
          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
          hprime_xx,hprime_xxT, &
          hprimewgll_xx,hprimewgll_xxT, &
          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
          kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
          c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
          c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
          c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
          ibool,idoubling,R_memory,epsilondev,epsilon_trace_over_3,one_minus_sum_beta, &
          alphaval,betaval,gammaval,factor_common,vx,vy,vz,vnspec, &
          COMPUTE_AND_STORE_STRAIN)


! this routine is optimized for NGLLX = NGLLY = NGLLZ = 5

  implicit none

  include "constants.h"

  ! include values created by the mesher
  ! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

  ! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: displ,accel
  ! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool

  ! x y and z contain r theta and phi
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  ! array with derivatives of Lagrange polynomials and precalculated products
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xxT,hprimewgll_xxT  
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  ! store anisotropic properties only where needed to save memory
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE) :: &
        kappahstore,muhstore,eta_anisostore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ISO_MANTLE) :: &
        kappavstore,muvstore

  ! arrays for full anisotropy only when needed
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE) :: &
        c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
        c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
        c36store,c44store,c45store,c46store,c55store,c56store,c66store

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: R_memory
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: epsilondev
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: epsilon_trace_over_3

  integer vx, vy, vz, vnspec

  ! [alpha,beta,gamma]val reduced to N_SLS and factor_common to N_SLS*NUM_NODES
  real(kind=CUSTOM_REAL), dimension(N_SLS, vx, vy, vz, vnspec) :: factor_common
  real(kind=CUSTOM_REAL), dimension(vx, vy, vz, vnspec) :: one_minus_sum_beta
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval

  ! array with the local to global mapping per slice
  integer, dimension(NSPEC_CRUST_MANTLE) :: idoubling

  ! gravity
  double precision, dimension(NRAD_GRAVITY) :: minus_gravity_table,density_table,minus_deriv_gravity_table

  ! for forward or backward simulations
  logical COMPUTE_AND_STORE_STRAIN

! local parameters
  ! Deville  
  ! manually inline the calls to the Deville et al. (2002) routines
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3  
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc, &
    newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,m2) :: B1_m1_m2_5points,B2_m1_m2_5points,B3_m1_m2_5points
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: C1_m1_m2_5points,C2_m1_m2_5points,C3_m1_m2_5points
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: E1_m1_m2_5points,E2_m1_m2_5points,E3_m1_m2_5points

  equivalence(dummyx_loc,B1_m1_m2_5points)
  equivalence(dummyy_loc,B2_m1_m2_5points)
  equivalence(dummyz_loc,B3_m1_m2_5points)
  equivalence(tempx1,C1_m1_m2_5points)
  equivalence(tempy1,C2_m1_m2_5points)
  equivalence(tempz1,C3_m1_m2_5points)
  equivalence(newtempx1,E1_m1_m2_5points)
  equivalence(newtempy1,E2_m1_m2_5points)
  equivalence(newtempz1,E3_m1_m2_5points)

  real(kind=CUSTOM_REAL), dimension(m2,NGLLX) :: &
    A1_mxm_m2_m1_5points,A2_mxm_m2_m1_5points,A3_mxm_m2_m1_5points
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: &
    C1_mxm_m2_m1_5points,C2_mxm_m2_m1_5points,C3_mxm_m2_m1_5points
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: &
    E1_mxm_m2_m1_5points,E2_mxm_m2_m1_5points,E3_mxm_m2_m1_5points

  equivalence(dummyx_loc,A1_mxm_m2_m1_5points)
  equivalence(dummyy_loc,A2_mxm_m2_m1_5points)
  equivalence(dummyz_loc,A3_mxm_m2_m1_5points)
  equivalence(tempx3,C1_mxm_m2_m1_5points)
  equivalence(tempy3,C2_mxm_m2_m1_5points)
  equivalence(tempz3,C3_mxm_m2_m1_5points)
  equivalence(newtempx3,E1_mxm_m2_m1_5points)
  equivalence(newtempy3,E2_mxm_m2_m1_5points)
  equivalence(newtempz3,E3_mxm_m2_m1_5points)

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sum_terms

  ! for attenuation
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ) :: &
    factor_common_c44_muv
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc
  real(kind=CUSTOM_REAL) R_xx_val1,R_yy_val1,R_xx_val2,R_yy_val2,R_xx_val3,R_yy_val3
  real(kind=CUSTOM_REAL) one_minus_sum_beta_use,minus_sum_beta

  ! the 21 coefficients for an anisotropic medium in reduced notation
  real(kind=CUSTOM_REAL) c11,c22,c33,c44,c55,c66,c12,c13,c23,c14,c24,c34,c15,c25,c35,c45,c16,c26,c36,c46,c56

  real(kind=CUSTOM_REAL) rhovphsq,sinphifour,cosphisq,sinphisq,costhetasq,rhovsvsq,sinthetasq, &
        cosphifour,costhetafour,rhovpvsq,sinthetafour,rhovshsq,cosfourphi, &
        costwotheta,cosfourtheta,sintwophisq,costheta,sinphi,sintheta,cosphi, &
        sintwotheta,costwophi,sintwophi,costwothetasq,costwophisq,phi,theta

  real(kind=CUSTOM_REAL) two_rhovpvsq,two_rhovphsq,two_rhovsvsq,two_rhovshsq
  real(kind=CUSTOM_REAL) four_rhovpvsq,four_rhovphsq,four_rhovsvsq,four_rhovshsq

  real(kind=CUSTOM_REAL) twoetaminone,etaminone,eta_aniso
  real(kind=CUSTOM_REAL) two_eta_aniso,four_eta_aniso,six_eta_aniso
  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl
  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl
  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz

  real(kind=CUSTOM_REAL) fac1,fac2,fac3
  real(kind=CUSTOM_REAL) lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) kappal,kappavl,kappahl,muvl,muhl

  ! for gravity
  double precision radius,rho,minus_g,minus_dg
  double precision minus_g_over_radius,minus_dg_plus_g_over_radius
  double precision cos_theta,sin_theta,cos_phi,sin_phi
  double precision cos_theta_sq,sin_theta_sq,cos_phi_sq,sin_phi_sq
  double precision factor,sx_l,sy_l,sz_l,gxl,gyl,gzl
  double precision Hxxl,Hyyl,Hzzl,Hxyl,Hxzl,Hyzl
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: rho_s_H
  real(kind=CUSTOM_REAL) sigma_yx,sigma_zx,sigma_zy

  integer ::i_sls,i_memory
  integer :: ispec,ispec_strain
  integer :: i,j,k 
  integer :: int_radius
  integer :: iglob1,iglob2,iglob3,iglob4,iglob5
  
! ****************************************************
!   big loop over all spectral elements in the solid
! ****************************************************

  do ispec = 1,NSPEC_CRUST_MANTLE

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1
    do k=1,NGLLZ
      do j=1,NGLLY

! way 1:      
!        do i=1,NGLLX
!            iglob = ibool(i,j,k,ispec)
!            dummyx_loc(i,j,k) = displ(1,iglob)
!            dummyy_loc(i,j,k) = displ(2,iglob)
!            dummyz_loc(i,j,k) = displ(3,iglob)
!        enddo

! way 2:
        ! since we know that NGLLX = 5, this should help pipelining
        iglob1 = ibool(1,j,k,ispec)
        iglob2 = ibool(2,j,k,ispec)
        iglob3 = ibool(3,j,k,ispec)
        iglob4 = ibool(4,j,k,ispec)
        iglob5 = ibool(5,j,k,ispec)
        
        dummyx_loc(1,j,k) = displ(1,iglob1)
        dummyy_loc(1,j,k) = displ(2,iglob1)
        dummyz_loc(1,j,k) = displ(3,iglob1)

        dummyx_loc(2,j,k) = displ(1,iglob2)
        dummyy_loc(2,j,k) = displ(2,iglob2)
        dummyz_loc(2,j,k) = displ(3,iglob2)

        dummyx_loc(3,j,k) = displ(1,iglob3)
        dummyy_loc(3,j,k) = displ(2,iglob3)
        dummyz_loc(3,j,k) = displ(3,iglob3)

        dummyx_loc(4,j,k) = displ(1,iglob4)
        dummyy_loc(4,j,k) = displ(2,iglob4)
        dummyz_loc(4,j,k) = displ(3,iglob4)

        dummyx_loc(5,j,k) = displ(1,iglob5)
        dummyy_loc(5,j,k) = displ(2,iglob5)
        dummyz_loc(5,j,k) = displ(3,iglob5)
        
      enddo
    enddo  
    do j=1,m2
      do i=1,m1
        C1_m1_m2_5points(i,j) = hprime_xx(i,1)*B1_m1_m2_5points(1,j) + &
                              hprime_xx(i,2)*B1_m1_m2_5points(2,j) + &
                              hprime_xx(i,3)*B1_m1_m2_5points(3,j) + &
                              hprime_xx(i,4)*B1_m1_m2_5points(4,j) + &
                              hprime_xx(i,5)*B1_m1_m2_5points(5,j)

        C2_m1_m2_5points(i,j) = hprime_xx(i,1)*B2_m1_m2_5points(1,j) + &
                              hprime_xx(i,2)*B2_m1_m2_5points(2,j) + &
                              hprime_xx(i,3)*B2_m1_m2_5points(3,j) + &
                              hprime_xx(i,4)*B2_m1_m2_5points(4,j) + &
                              hprime_xx(i,5)*B2_m1_m2_5points(5,j)

        C3_m1_m2_5points(i,j) = hprime_xx(i,1)*B3_m1_m2_5points(1,j) + &
                              hprime_xx(i,2)*B3_m1_m2_5points(2,j) + &
                              hprime_xx(i,3)*B3_m1_m2_5points(3,j) + &
                              hprime_xx(i,4)*B3_m1_m2_5points(4,j) + &
                              hprime_xx(i,5)*B3_m1_m2_5points(5,j)
      enddo
    enddo  
    do j=1,m1
      do i=1,m1
        ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
        do k = 1,NGLLX
          tempx2(i,j,k) = dummyx_loc(i,1,k)*hprime_xxT(1,j) + &
                        dummyx_loc(i,2,k)*hprime_xxT(2,j) + &
                        dummyx_loc(i,3,k)*hprime_xxT(3,j) + &
                        dummyx_loc(i,4,k)*hprime_xxT(4,j) + &
                        dummyx_loc(i,5,k)*hprime_xxT(5,j)

          tempy2(i,j,k) = dummyy_loc(i,1,k)*hprime_xxT(1,j) + &
                        dummyy_loc(i,2,k)*hprime_xxT(2,j) + &
                        dummyy_loc(i,3,k)*hprime_xxT(3,j) + &
                        dummyy_loc(i,4,k)*hprime_xxT(4,j) + &
                        dummyy_loc(i,5,k)*hprime_xxT(5,j)

          tempz2(i,j,k) = dummyz_loc(i,1,k)*hprime_xxT(1,j) + &
                        dummyz_loc(i,2,k)*hprime_xxT(2,j) + &
                        dummyz_loc(i,3,k)*hprime_xxT(3,j) + &
                        dummyz_loc(i,4,k)*hprime_xxT(4,j) + &
                        dummyz_loc(i,5,k)*hprime_xxT(5,j)
        enddo
      enddo
    enddo
    do j=1,m1
      do i=1,m2
        C1_mxm_m2_m1_5points(i,j) = A1_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                  A1_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                  A1_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                  A1_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                  A1_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)

        C2_mxm_m2_m1_5points(i,j) = A2_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                  A2_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                  A2_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                  A2_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                  A2_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)

        C3_mxm_m2_m1_5points(i,j) = A3_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                  A3_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                  A3_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                  A3_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                  A3_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)
      enddo
    enddo
  
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          ! get derivatives of ux, uy and uz with respect to x, y and z
          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)

          ! compute the jacobian
          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))

          duxdxl = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
          duxdyl = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
          duxdzl = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)

          duydxl = xixl*tempy1(i,j,k) + etaxl*tempy2(i,j,k) + gammaxl*tempy3(i,j,k)
          duydyl = xiyl*tempy1(i,j,k) + etayl*tempy2(i,j,k) + gammayl*tempy3(i,j,k)
          duydzl = xizl*tempy1(i,j,k) + etazl*tempy2(i,j,k) + gammazl*tempy3(i,j,k)

          duzdxl = xixl*tempz1(i,j,k) + etaxl*tempz2(i,j,k) + gammaxl*tempz3(i,j,k)
          duzdyl = xiyl*tempz1(i,j,k) + etayl*tempz2(i,j,k) + gammayl*tempz3(i,j,k)
          duzdzl = xizl*tempz1(i,j,k) + etazl*tempz2(i,j,k) + gammazl*tempz3(i,j,k)

          ! precompute some sums to save CPU time
          duxdxl_plus_duydyl = duxdxl + duydyl
          duxdxl_plus_duzdzl = duxdxl + duzdzl
          duydyl_plus_duzdzl = duydyl + duzdzl
          duxdyl_plus_duydxl = duxdyl + duydxl
          duzdxl_plus_duxdzl = duzdxl + duxdzl
          duzdyl_plus_duydzl = duzdyl + duydzl

          ! compute deviatoric strain
          if (COMPUTE_AND_STORE_STRAIN) then
            if(NSPEC_CRUST_MANTLE_STRAIN_ONLY == 1) then
              ispec_strain = 1
            else
              ispec_strain = ispec
            endif
            epsilon_trace_over_3(i,j,k,ispec_strain) = ONE_THIRD * (duxdxl + duydyl + duzdzl)
            epsilondev_loc(1,i,j,k) = duxdxl - epsilon_trace_over_3(i,j,k,ispec_strain)
            epsilondev_loc(2,i,j,k) = duydyl - epsilon_trace_over_3(i,j,k,ispec_strain)
            epsilondev_loc(3,i,j,k) = 0.5 * duxdyl_plus_duydxl
            epsilondev_loc(4,i,j,k) = 0.5 * duzdxl_plus_duxdzl
            epsilondev_loc(5,i,j,k) = 0.5 * duzdyl_plus_duydzl
          endif

          ! precompute terms for attenuation if needed
          if(ATTENUATION_VAL) then
            one_minus_sum_beta_use = one_minus_sum_beta(i,j,k,ispec)
            minus_sum_beta =  one_minus_sum_beta_use - 1.0
          endif

          !
          ! compute either isotropic or anisotropic elements
          !
          if(ANISOTROPIC_3D_MANTLE_VAL) then

            c11 = c11store(i,j,k,ispec)
            c12 = c12store(i,j,k,ispec)
            c13 = c13store(i,j,k,ispec)
            c14 = c14store(i,j,k,ispec)
            c15 = c15store(i,j,k,ispec)
            c16 = c16store(i,j,k,ispec)
            c22 = c22store(i,j,k,ispec)
            c23 = c23store(i,j,k,ispec)
            c24 = c24store(i,j,k,ispec)
            c25 = c25store(i,j,k,ispec)
            c26 = c26store(i,j,k,ispec)
            c33 = c33store(i,j,k,ispec)
            c34 = c34store(i,j,k,ispec)
            c35 = c35store(i,j,k,ispec)
            c36 = c36store(i,j,k,ispec)
            c44 = c44store(i,j,k,ispec)
            c45 = c45store(i,j,k,ispec)
            c46 = c46store(i,j,k,ispec)
            c55 = c55store(i,j,k,ispec)
            c56 = c56store(i,j,k,ispec)
            c66 = c66store(i,j,k,ispec)

            if(ATTENUATION_VAL) then
              mul = c44
              c11 = c11 + FOUR_THIRDS * minus_sum_beta * mul
              c12 = c12 - TWO_THIRDS * minus_sum_beta * mul
              c13 = c13 - TWO_THIRDS * minus_sum_beta * mul
              c22 = c22 + FOUR_THIRDS * minus_sum_beta * mul
              c23 = c23 - TWO_THIRDS * minus_sum_beta * mul
              c33 = c33 + FOUR_THIRDS * minus_sum_beta * mul
              c44 = c44 + minus_sum_beta * mul
              c55 = c55 + minus_sum_beta * mul
              c66 = c66 + minus_sum_beta * mul
            endif

            !mimik: apparent velocity shift
            if( ATTENUATION_MIMIK) then
              mul = c44
              c11 = c11 + FOUR_THIRDS * (1.0-ATTENUATION_MIMIK_FACTOR) * mul
              c12 = c12 - TWO_THIRDS * (1.0-ATTENUATION_MIMIK_FACTOR) * mul
              c13 = c13 - TWO_THIRDS * (1.0-ATTENUATION_MIMIK_FACTOR) * mul
              c22 = c22 + FOUR_THIRDS * (1.0-ATTENUATION_MIMIK_FACTOR) * mul
              c23 = c23 - TWO_THIRDS * (1.0-ATTENUATION_MIMIK_FACTOR) * mul
              c33 = c33 + FOUR_THIRDS * (1.0-ATTENUATION_MIMIK_FACTOR) * mul
              c44 = c44 + (1.0-ATTENUATION_MIMIK_FACTOR) * mul
              c55 = c55 + (1.0-ATTENUATION_MIMIK_FACTOR) * mul
              c66 = c66 + (1.0-ATTENUATION_MIMIK_FACTOR) * mul
            endif

            sigma_xx = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl + &
                     c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl

            sigma_yy = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl + &
                     c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl

            sigma_zz = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl + &
                     c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl

            sigma_xy = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl + &
                     c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl

            sigma_xz = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl + &
                     c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl

            sigma_yz = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl + &
                     c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl

          else

            ! do not use transverse isotropy except if element is between d220 and Moho
            if(.not. (TRANSVERSE_ISOTROPY_VAL .and. (idoubling(ispec)==IFLAG_220_80 &
                  .or. idoubling(ispec)==IFLAG_80_MOHO))) then

              ! layer with no transverse isotropy, use kappav and muv
              kappal = kappavstore(i,j,k,ispec)
              mul = muvstore(i,j,k,ispec)

              ! use unrelaxed parameters if attenuation
              if(ATTENUATION_VAL) mul = mul * one_minus_sum_beta_use

              !mimik: apparent velocity shift
              if( ATTENUATION_MIMIK) mul = mul * ATTENUATION_MIMIK_FACTOR

              lambdalplus2mul = kappal + FOUR_THIRDS * mul
              lambdal = lambdalplus2mul - 2.*mul

              ! compute stress sigma
              sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
              sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
              sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl

              sigma_xy = mul*duxdyl_plus_duydxl
              sigma_xz = mul*duzdxl_plus_duxdzl
              sigma_yz = mul*duzdyl_plus_duydzl

            else

              ! use Kappa and mu from transversely isotropic model
              kappavl = kappavstore(i,j,k,ispec)
              muvl = muvstore(i,j,k,ispec)

              kappahl = kappahstore(i,j,k,ispec)
              muhl = muhstore(i,j,k,ispec)

              ! use unrelaxed parameters if attenuation
              ! eta does not need to be shifted since it is a ratio
              if(ATTENUATION_VAL) then
                muvl = muvl * one_minus_sum_beta_use
                muhl = muhl * one_minus_sum_beta_use
              endif

              !mimik: apparent velocity shift
              if( ATTENUATION_MIMIK) then
                muvl = muvl * ATTENUATION_MIMIK_FACTOR
                muhl = muhl * ATTENUATION_MIMIK_FACTOR  
              endif

              rhovpvsq = kappavl + FOUR_THIRDS * muvl  !!! that is C
              rhovphsq = kappahl + FOUR_THIRDS * muhl  !!! that is A

              rhovsvsq = muvl  !!! that is L
              rhovshsq = muhl  !!! that is N

              eta_aniso = eta_anisostore(i,j,k,ispec)  !!! that is  F / (A - 2 L)

              ! use mesh coordinates to get theta and phi
              ! ystore and zstore contain theta and phi

              iglob1 = ibool(i,j,k,ispec)
              theta = ystore(iglob1)
              phi = zstore(iglob1)

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

              costwotheta = cos(2.*theta)
              sintwotheta = sin(2.*theta)
              costwophi = cos(2.*phi)
              sintwophi = sin(2.*phi)

              cosfourtheta = cos(4.*theta)
              cosfourphi = cos(4.*phi)

              costwothetasq = costwotheta * costwotheta

              costwophisq = costwophi * costwophi
              sintwophisq = sintwophi * sintwophi

              etaminone = eta_aniso - 1.
              twoetaminone = 2. * eta_aniso - 1.

              ! precompute some products to reduce the CPU time
              two_eta_aniso = 2.*eta_aniso
              four_eta_aniso = 4.*eta_aniso
              six_eta_aniso = 6.*eta_aniso

              two_rhovpvsq = 2.*rhovpvsq
              two_rhovphsq = 2.*rhovphsq
              two_rhovsvsq = 2.*rhovsvsq
              two_rhovshsq = 2.*rhovshsq

              four_rhovpvsq = 4.*rhovpvsq
              four_rhovphsq = 4.*rhovphsq
              four_rhovsvsq = 4.*rhovsvsq
              four_rhovshsq = 4.*rhovshsq

              ! the 21 anisotropic coefficients computed using Mathematica

              c11 = rhovphsq*sinphifour + 2.*cosphisq*sinphisq* &
                  (rhovphsq*costhetasq + (eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)* &
                  sinthetasq) + cosphifour* &
                  (rhovphsq*costhetafour + 2.*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)* &
                  costhetasq*sinthetasq + rhovpvsq*sinthetafour)

              c12 = ((rhovphsq - two_rhovshsq)*(3. + cosfourphi)*costhetasq)/4. - &
                  four_rhovshsq*cosphisq*costhetasq*sinphisq + &
                  (rhovphsq*(11. + 4.*costwotheta + cosfourtheta)*sintwophisq)/32. + &
                  eta_aniso*(rhovphsq - two_rhovsvsq)*(cosphifour + &
                  2.*cosphisq*costhetasq*sinphisq + sinphifour)*sinthetasq + &
                  rhovpvsq*cosphisq*sinphisq*sinthetafour - &
                  rhovsvsq*sintwophisq*sinthetafour

              c13 = (cosphisq*(rhovphsq + six_eta_aniso*rhovphsq + rhovpvsq - four_rhovsvsq - &
                  12.*eta_aniso*rhovsvsq + (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - &
                  four_eta_aniso*rhovsvsq)*cosfourtheta))/8. + &
                  sinphisq*(eta_aniso*(rhovphsq - two_rhovsvsq)*costhetasq + &
                  (rhovphsq - two_rhovshsq)*sinthetasq)

              c14 = costheta*sinphi*((cosphisq* &
                   (-rhovphsq + rhovpvsq + four_rhovshsq - four_rhovsvsq + &
                    (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq - &
                    four_eta_aniso*rhovsvsq)*costwotheta))/2. + &
                    (etaminone*rhovphsq + 2.*(rhovshsq - eta_aniso*rhovsvsq))*sinphisq)* sintheta

              c15 = cosphi*costheta*((cosphisq* (-rhovphsq + rhovpvsq + &
                    (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)* &
                    costwotheta))/2. + etaminone*(rhovphsq - two_rhovsvsq)*sinphisq)*sintheta

              c16 = (cosphi*sinphi*(cosphisq* (-rhovphsq + rhovpvsq + &
                    (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq - &
                    four_eta_aniso*rhovsvsq)*costwotheta) + &
                    2.*etaminone*(rhovphsq - two_rhovsvsq)*sinphisq)*sinthetasq)/2.

              c22 = rhovphsq*cosphifour + 2.*cosphisq*sinphisq* &
                  (rhovphsq*costhetasq + (eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)* &
                  sinthetasq) + sinphifour* &
                  (rhovphsq*costhetafour + 2.*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)* &
                  costhetasq*sinthetasq + rhovpvsq*sinthetafour)

              c23 = ((rhovphsq + six_eta_aniso*rhovphsq + rhovpvsq - four_rhovsvsq - 12.*eta_aniso*rhovsvsq + &
                   (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)* &
                    cosfourtheta)*sinphisq)/8. + &
                    cosphisq*(eta_aniso*(rhovphsq - two_rhovsvsq)*costhetasq + &
                    (rhovphsq - two_rhovshsq)*sinthetasq)

              c24 = costheta*sinphi*(etaminone*(rhovphsq - two_rhovsvsq)*cosphisq + &
                    ((-rhovphsq + rhovpvsq + (twoetaminone*rhovphsq - rhovpvsq + &
                    four_rhovsvsq - four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)/2.)*sintheta

              c25 = cosphi*costheta*((etaminone*rhovphsq + 2.*(rhovshsq - eta_aniso*rhovsvsq))* &
                    cosphisq + ((-rhovphsq + rhovpvsq + four_rhovshsq - four_rhovsvsq + &
                     (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq - &
                    four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)/2.)*sintheta

              c26 = (cosphi*sinphi*(2.*etaminone*(rhovphsq - two_rhovsvsq)*cosphisq + &
                      (-rhovphsq + rhovpvsq + (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq - &
                      four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)*sinthetasq)/2.

              c33 = rhovpvsq*costhetafour + 2.*(eta_aniso*(rhovphsq - two_rhovsvsq) + two_rhovsvsq)* &
                    costhetasq*sinthetasq + rhovphsq*sinthetafour

              c34 = -((rhovphsq - rhovpvsq + (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq &
                       - four_eta_aniso*rhovsvsq)*costwotheta)*sinphi*sintwotheta)/4.

              c35 = -(cosphi*(rhovphsq - rhovpvsq + &
                    (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)* &
                    costwotheta)*sintwotheta)/4.

              c36 = -((rhovphsq - rhovpvsq - four_rhovshsq + four_rhovsvsq + &
                    (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)* &
                    costwotheta)*sintwophi*sinthetasq)/4.

              c44 = cosphisq*(rhovsvsq*costhetasq + rhovshsq*sinthetasq) + &
                    sinphisq*(rhovsvsq*costwothetasq + &
                    (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + four_eta_aniso*rhovsvsq)*costhetasq* sinthetasq)

              c45 = ((rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq - two_rhovshsq - two_rhovsvsq + &
                    four_eta_aniso*rhovsvsq + (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + &
                    4.*etaminone*rhovsvsq)*costwotheta)*sintwophi*sinthetasq)/4.

              c46 = -(cosphi*costheta*((rhovshsq - rhovsvsq)*cosphisq - &
                      ((rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq - two_rhovshsq - two_rhovsvsq + &
                      four_eta_aniso*rhovsvsq + (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + &
                      four_rhovsvsq - four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)/2.)* sintheta)

              c55 = sinphisq*(rhovsvsq*costhetasq + rhovshsq*sinthetasq) + &
                  cosphisq*(rhovsvsq*costwothetasq + &
                  (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + four_eta_aniso*rhovsvsq)*costhetasq* sinthetasq)

              c56 = costheta*sinphi*((cosphisq* &
                  (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq - two_rhovshsq - two_rhovsvsq + &
                  four_eta_aniso*rhovsvsq + (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + &
                  four_rhovsvsq - four_eta_aniso*rhovsvsq)*costwotheta))/2. + &
                  (-rhovshsq + rhovsvsq)*sinphisq)*sintheta

              c66 = rhovshsq*costwophisq*costhetasq - &
                  2.*(rhovphsq - two_rhovshsq)*cosphisq*costhetasq*sinphisq + &
                  (rhovphsq*(11. + 4.*costwotheta + cosfourtheta)*sintwophisq)/32. - &
                  (rhovsvsq*(-6. - 2.*cosfourphi + cos(4.*phi - 2.*theta) - 2.*costwotheta + &
                  cos(2.*(2.*phi + theta)))*sinthetasq)/8. + &
                  rhovpvsq*cosphisq*sinphisq*sinthetafour - &
                  (eta_aniso*(rhovphsq - two_rhovsvsq)*sintwophisq*sinthetafour)/2.

              ! general expression of stress tensor for full Cijkl with 21 coefficients
              sigma_xx = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl + &
                       c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl

              sigma_yy = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl + &
                       c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl

              sigma_zz = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl + &
                       c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl

              sigma_xy = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl + &
                       c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl

              sigma_xz = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl + &
                       c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl

              sigma_yz = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl + &
                       c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl

            endif

          endif   ! end of test whether isotropic or anisotropic element

          ! subtract memory variables if attenuation
          if(ATTENUATION_VAL) then
! way 1:
!            do i_sls = 1,N_SLS
!              R_xx_val = R_memory(1,i_sls,i,j,k,ispec)
!              R_yy_val = R_memory(2,i_sls,i,j,k,ispec)
!              sigma_xx = sigma_xx - R_xx_val
!              sigma_yy = sigma_yy - R_yy_val
!              sigma_zz = sigma_zz + R_xx_val + R_yy_val
!              sigma_xy = sigma_xy - R_memory(3,i_sls,i,j,k,ispec)
!              sigma_xz = sigma_xz - R_memory(4,i_sls,i,j,k,ispec)
!              sigma_yz = sigma_yz - R_memory(5,i_sls,i,j,k,ispec)            
!            enddo

! way 2:
! note: this should help compilers to pipeline the code and make better use of the cache;
!          depending on compilers, it can further decrease the computation time by ~ 30%.
!          by default, N_SLS = 3, therefor we take steps of 3
            do i_sls = 1,mod(N_SLS,3)
              R_xx_val1 = R_memory(1,i_sls,i,j,k,ispec)
              R_yy_val1 = R_memory(2,i_sls,i,j,k,ispec)
              sigma_xx = sigma_xx - R_xx_val1
              sigma_yy = sigma_yy - R_yy_val1
              sigma_zz = sigma_zz + R_xx_val1 + R_yy_val1
              sigma_xy = sigma_xy - R_memory(3,i_sls,i,j,k,ispec)
              sigma_xz = sigma_xz - R_memory(4,i_sls,i,j,k,ispec)
              sigma_yz = sigma_yz - R_memory(5,i_sls,i,j,k,ispec)            
            enddo
            
            do i_sls = mod(N_SLS,3)+1,N_SLS,3
              R_xx_val1 = R_memory(1,i_sls,i,j,k,ispec)
              R_yy_val1 = R_memory(2,i_sls,i,j,k,ispec)
              sigma_xx = sigma_xx - R_xx_val1
              sigma_yy = sigma_yy - R_yy_val1
              sigma_zz = sigma_zz + R_xx_val1 + R_yy_val1
              sigma_xy = sigma_xy - R_memory(3,i_sls,i,j,k,ispec)
              sigma_xz = sigma_xz - R_memory(4,i_sls,i,j,k,ispec)
              sigma_yz = sigma_yz - R_memory(5,i_sls,i,j,k,ispec)

              R_xx_val2 = R_memory(1,i_sls+1,i,j,k,ispec)
              R_yy_val2 = R_memory(2,i_sls+1,i,j,k,ispec)
              sigma_xx = sigma_xx - R_xx_val2
              sigma_yy = sigma_yy - R_yy_val2
              sigma_zz = sigma_zz + R_xx_val2 + R_yy_val2
              sigma_xy = sigma_xy - R_memory(3,i_sls+1,i,j,k,ispec)
              sigma_xz = sigma_xz - R_memory(4,i_sls+1,i,j,k,ispec)
              sigma_yz = sigma_yz - R_memory(5,i_sls+1,i,j,k,ispec)

              R_xx_val3 = R_memory(1,i_sls+2,i,j,k,ispec)
              R_yy_val3 = R_memory(2,i_sls+2,i,j,k,ispec)
              sigma_xx = sigma_xx - R_xx_val3
              sigma_yy = sigma_yy - R_yy_val3
              sigma_zz = sigma_zz + R_xx_val3 + R_yy_val3
              sigma_xy = sigma_xy - R_memory(3,i_sls+2,i,j,k,ispec)
              sigma_xz = sigma_xz - R_memory(4,i_sls+2,i,j,k,ispec)
              sigma_yz = sigma_yz - R_memory(5,i_sls+2,i,j,k,ispec)
            enddo

          endif

          ! define symmetric components of sigma for gravity
          sigma_yx = sigma_xy
          sigma_zx = sigma_xz
          sigma_zy = sigma_yz

          ! compute non-symmetric terms for gravity
          if(GRAVITY_VAL) then

            ! use mesh coordinates to get theta and phi
            ! x y and z contain r theta and phi
            iglob1 = ibool(i,j,k,ispec)
            theta = ystore(iglob1)
            phi = zstore(iglob1)

            cos_theta = dcos(dble(theta))
            sin_theta = dsin(dble(theta))
            cos_phi = dcos(dble(phi))
            sin_phi = dsin(dble(phi))

            cos_theta_sq = cos_theta**2
            sin_theta_sq = sin_theta**2
            cos_phi_sq = cos_phi**2
            sin_phi_sq = sin_phi**2

            ! get g, rho and dg/dr=dg
            ! spherical components of the gravitational acceleration
            ! for efficiency replace with lookup table every 100 m in radial direction
            radius = dble(xstore(iglob1))
            int_radius = nint(radius * R_EARTH_KM * 10.d0)
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
            iglob1 = ibool(i,j,k,ispec)

            ! distinguish between single and double precision for reals
            if(CUSTOM_REAL == SIZE_REAL) then

              ! get displacement and multiply by density to compute G tensor
              sx_l = rho * dble(displ(1,iglob1))
              sy_l = rho * dble(displ(2,iglob1))
              sz_l = rho * dble(displ(3,iglob1))

              ! compute G tensor from s . g and add to sigma (not symmetric)
              sigma_xx = sigma_xx + sngl(sy_l*gyl + sz_l*gzl)
              sigma_yy = sigma_yy + sngl(sx_l*gxl + sz_l*gzl)
              sigma_zz = sigma_zz + sngl(sx_l*gxl + sy_l*gyl)

              sigma_xy = sigma_xy - sngl(sx_l * gyl)
              sigma_yx = sigma_yx - sngl(sy_l * gxl)

              sigma_xz = sigma_xz - sngl(sx_l * gzl)
              sigma_zx = sigma_zx - sngl(sz_l * gxl)

              sigma_yz = sigma_yz - sngl(sy_l * gzl)
              sigma_zy = sigma_zy - sngl(sz_l * gyl)

              ! precompute vector
              factor = dble(jacobianl) * wgll_cube(i,j,k)
              rho_s_H(1,i,j,k) = sngl(factor * (sx_l * Hxxl + sy_l * Hxyl + sz_l * Hxzl))
              rho_s_H(2,i,j,k) = sngl(factor * (sx_l * Hxyl + sy_l * Hyyl + sz_l * Hyzl))
              rho_s_H(3,i,j,k) = sngl(factor * (sx_l * Hxzl + sy_l * Hyzl + sz_l * Hzzl))

            else
          
              ! get displacement and multiply by density to compute G tensor
              sx_l = rho * displ(1,iglob1)
              sy_l = rho * displ(2,iglob1)
              sz_l = rho * displ(3,iglob1)

              ! compute G tensor from s . g and add to sigma (not symmetric)
              sigma_xx = sigma_xx + sy_l*gyl + sz_l*gzl
              sigma_yy = sigma_yy + sx_l*gxl + sz_l*gzl
              sigma_zz = sigma_zz + sx_l*gxl + sy_l*gyl

              sigma_xy = sigma_xy - sx_l * gyl
              sigma_yx = sigma_yx - sy_l * gxl

              sigma_xz = sigma_xz - sx_l * gzl
              sigma_zx = sigma_zx - sz_l * gxl

              sigma_yz = sigma_yz - sy_l * gzl
              sigma_zy = sigma_zy - sz_l * gyl

              ! precompute vector
              factor = jacobianl * wgll_cube(i,j,k)
              rho_s_H(1,i,j,k) = factor * (sx_l * Hxxl + sy_l * Hxyl + sz_l * Hxzl)
              rho_s_H(2,i,j,k) = factor * (sx_l * Hxyl + sy_l * Hyyl + sz_l * Hyzl)
              rho_s_H(3,i,j,k) = factor * (sx_l * Hxzl + sy_l * Hyzl + sz_l * Hzzl)

            endif

          endif  ! end of section with gravity terms

          ! form dot product with test vector, non-symmetric form
          tempx1(i,j,k) = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl)
          tempy1(i,j,k) = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl)
          tempz1(i,j,k) = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl)

          tempx2(i,j,k) = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl)
          tempy2(i,j,k) = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl)
          tempz2(i,j,k) = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl)

          tempx3(i,j,k) = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl)
          tempy3(i,j,k) = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl)
          tempz3(i,j,k) = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl)
        enddo ! NGLLX
      enddo ! NGLLY
    enddo ! NGLLZ

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1
    do j=1,m2
      do i=1,m1
        E1_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C1_m1_m2_5points(1,j) + &
                              hprimewgll_xxT(i,2)*C1_m1_m2_5points(2,j) + &
                              hprimewgll_xxT(i,3)*C1_m1_m2_5points(3,j) + &
                              hprimewgll_xxT(i,4)*C1_m1_m2_5points(4,j) + &
                              hprimewgll_xxT(i,5)*C1_m1_m2_5points(5,j)

        E2_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C2_m1_m2_5points(1,j) + &
                              hprimewgll_xxT(i,2)*C2_m1_m2_5points(2,j) + &
                              hprimewgll_xxT(i,3)*C2_m1_m2_5points(3,j) + &
                              hprimewgll_xxT(i,4)*C2_m1_m2_5points(4,j) + &
                              hprimewgll_xxT(i,5)*C2_m1_m2_5points(5,j)

        E3_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C3_m1_m2_5points(1,j) + &
                              hprimewgll_xxT(i,2)*C3_m1_m2_5points(2,j) + &
                              hprimewgll_xxT(i,3)*C3_m1_m2_5points(3,j) + &
                              hprimewgll_xxT(i,4)*C3_m1_m2_5points(4,j) + &
                              hprimewgll_xxT(i,5)*C3_m1_m2_5points(5,j)
      enddo
    enddo
    do i=1,m1
      do j=1,m1
        ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
        do k = 1,NGLLX
          newtempx2(i,j,k) = tempx2(i,1,k)*hprimewgll_xx(1,j) + &
                           tempx2(i,2,k)*hprimewgll_xx(2,j) + &
                           tempx2(i,3,k)*hprimewgll_xx(3,j) + &
                           tempx2(i,4,k)*hprimewgll_xx(4,j) + &
                           tempx2(i,5,k)*hprimewgll_xx(5,j)

          newtempy2(i,j,k) = tempy2(i,1,k)*hprimewgll_xx(1,j) + &
                           tempy2(i,2,k)*hprimewgll_xx(2,j) + &
                           tempy2(i,3,k)*hprimewgll_xx(3,j) + &
                           tempy2(i,4,k)*hprimewgll_xx(4,j) + &
                           tempy2(i,5,k)*hprimewgll_xx(5,j)

          newtempz2(i,j,k) = tempz2(i,1,k)*hprimewgll_xx(1,j) + &
                           tempz2(i,2,k)*hprimewgll_xx(2,j) + &
                           tempz2(i,3,k)*hprimewgll_xx(3,j) + &
                           tempz2(i,4,k)*hprimewgll_xx(4,j) + &
                           tempz2(i,5,k)*hprimewgll_xx(5,j)
        enddo
      enddo
    enddo
    do j=1,m1
      do i=1,m2
        E1_mxm_m2_m1_5points(i,j) = C1_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
                                  C1_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
                                  C1_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
                                  C1_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
                                  C1_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)

        E2_mxm_m2_m1_5points(i,j) = C2_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
                                  C2_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
                                  C2_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
                                  C2_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
                                  C2_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)

        E3_mxm_m2_m1_5points(i,j) = C3_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
                                  C3_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
                                  C3_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
                                  C3_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
                                  C3_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)
      enddo
    enddo

    do k=1,NGLLZ
      do j=1,NGLLY
      
! way 1:
! this seems to be still the fastest way here.      
        fac1 = wgllwgll_yz(j,k)
        do i=1,NGLLX
          fac2 = wgllwgll_xz(i,k)
          fac3 = wgllwgll_xy(i,j)

          ! sum contributions
          sum_terms(1,i,j,k) = - (fac1*newtempx1(i,j,k) + fac2*newtempx2(i,j,k) + fac3*newtempx3(i,j,k))
          sum_terms(2,i,j,k) = - (fac1*newtempy1(i,j,k) + fac2*newtempy2(i,j,k) + fac3*newtempy3(i,j,k))
          sum_terms(3,i,j,k) = - (fac1*newtempz1(i,j,k) + fac2*newtempz2(i,j,k) + fac3*newtempz3(i,j,k))

          if(GRAVITY_VAL) sum_terms(:,i,j,k) = sum_terms(:,i,j,k) + rho_s_H(:,i,j,k)

        enddo ! NGLLX

      enddo ! NGLLY
    enddo ! NGLLZ

    ! sum contributions from each element to the global mesh and add gravity terms
    do k=1,NGLLZ
      do j=1,NGLLY
! way 1:      
!        do i=1,NGLLX
!          iglob = ibool(i,j,k,ispec)
!          accel(:,iglob) = accel(:,iglob) + sum_terms(:,i,j,k)          
!        enddo

! way 2:
        accel(:,ibool(1,j,k,ispec)) = accel(:,ibool(1,j,k,ispec)) + sum_terms(:,1,j,k)          
        accel(:,ibool(2,j,k,ispec)) = accel(:,ibool(2,j,k,ispec)) + sum_terms(:,2,j,k)          
        accel(:,ibool(3,j,k,ispec)) = accel(:,ibool(3,j,k,ispec)) + sum_terms(:,3,j,k)          
        accel(:,ibool(4,j,k,ispec)) = accel(:,ibool(4,j,k,ispec)) + sum_terms(:,4,j,k)          
        accel(:,ibool(5,j,k,ispec)) = accel(:,ibool(5,j,k,ispec)) + sum_terms(:,5,j,k)          

      enddo
    enddo

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

    if(ATTENUATION_VAL) then

      ! use Runge-Kutta scheme to march in time

      ! get coefficients for that standard linear solid
      ! IMPROVE we use mu_v here even if there is some anisotropy
      ! IMPROVE we should probably use an average value instead

! way 1:
! it still seems to be the fastest way here.
      do i_sls = 1,N_SLS      
        ! reformatted R_memory to handle large factor_common and reduced [alpha,beta,gamma]val
        factor_common_c44_muv = factor_common(i_sls,:,:,:,ispec)

        if(ANISOTROPIC_3D_MANTLE_VAL) then
          factor_common_c44_muv = factor_common_c44_muv * c44store(:,:,:,ispec)
        else
          factor_common_c44_muv = factor_common_c44_muv * muvstore(:,:,:,ispec)
        endif
      
        do i_memory = 1,5
          R_memory(i_memory,i_sls,:,:,:,ispec) = alphaval(i_sls) * &
                    R_memory(i_memory,i_sls,:,:,:,ispec) + &
                    factor_common_c44_muv * &
                    (betaval(i_sls) * epsilondev(i_memory,:,:,:,ispec) + &
                    gammaval(i_sls) * epsilondev_loc(i_memory,:,:,:))
        enddo
      enddo

    endif

    ! save deviatoric strain for Runge-Kutta scheme
    if(COMPUTE_AND_STORE_STRAIN) epsilondev(:,:,:,:,ispec) = epsilondev_loc(:,:,:,:)

  enddo   ! spectral element loop NSPEC_CRUST_MANTLE

  end subroutine compute_forces_crust_mantle_Dev

