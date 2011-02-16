!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2011
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


  subroutine compute_forces_inner_core_Dev(minus_gravity_table,density_table,minus_deriv_gravity_table, &
          displ,accel,xstore,ystore,zstore, &
          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
          hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
          kappavstore,muvstore,ibool,idoubling, &
          c11store,c33store,c12store,c13store,c44store,R_memory,epsilondev,epsilon_trace_over_3,&
          one_minus_sum_beta,alphaval,betaval,gammaval,factor_common, &
          vx,vy,vz,vnspec)

  implicit none

  include "constants.h"

! include values created by the mesher
! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

  ! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: displ,accel
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool

  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE) :: xstore,ystore,zstore

  ! arrays with mesh parameters per slice
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: xix,xiy,xiz, &
                      etax,etay,etaz,gammax,gammay,gammaz

  ! for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  ! variable lengths for factor_common and one_minus_sum_beta
  integer vx, vy, vz, vnspec
  real(kind=CUSTOM_REAL), dimension(N_SLS, vx, vy, vz, vnspec) :: factor_common
  real(kind=CUSTOM_REAL), dimension(vx, vy, vz, vnspec) :: one_minus_sum_beta

  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ) :: factor_common_use
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval

  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ATTENUATION) :: R_memory
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: epsilondev
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: epsilon_trace_over_3

  ! array with derivatives of Lagrange polynomials and precalculated products
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xxT,hprimewgll_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: kappavstore,muvstore

!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: &
!    c11store,c33store,c12store,c13store,c44store
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_IC) :: &
    c11store,c33store,c12store,c13store,c44store

  ! array with the local to global mapping per slice
  integer, dimension(NSPEC_INNER_CORE) :: idoubling

  double precision, dimension(NRAD_GRAVITY) :: minus_gravity_table,density_table,minus_deriv_gravity_table

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
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc
  real(kind=CUSTOM_REAL) R_xx_val1,R_yy_val1,R_xx_val2,R_yy_val2,R_xx_val3,R_yy_val3

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz

  real(kind=CUSTOM_REAL) fac1,fac2,fac3,templ
  real(kind=CUSTOM_REAL) lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) kappal

  real(kind=CUSTOM_REAL) minus_sum_beta
  real(kind=CUSTOM_REAL) c11l,c33l,c12l,c13l,c44l

  ! for gravity
  double precision radius,rho,minus_g,minus_dg
  double precision minus_g_over_radius,minus_dg_plus_g_over_radius
  double precision cos_theta,sin_theta,cos_phi,sin_phi
  double precision cos_theta_sq,sin_theta_sq,cos_phi_sq,sin_phi_sq
  double precision theta,phi,factor,gxl,gyl,gzl,sx_l,sy_l,sz_l
  double precision Hxxl,Hyyl,Hzzl,Hxyl,Hxzl,Hyzl
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: rho_s_H
  real(kind=CUSTOM_REAL) sigma_yx,sigma_zx,sigma_zy

  !real(kind=CUSTOM_REAL), dimension(5) :: dummy

  integer :: int_radius
  integer :: ispec,ispec_strain
  integer :: i,j,k !,l
  integer :: i_SLS,i_memory,imodulo_N_SLS
  integer :: iglob1,iglob2,iglob3,iglob4,iglob5


! ****************************************************
!   big loop over all spectral elements in the solid
! ****************************************************

  imodulo_N_SLS = mod(N_SLS,3)

  do ispec = 1,NSPEC_INNER_CORE

    ! exclude fictitious elements in central cube
    if(idoubling(ispec) /= IFLAG_IN_FICTITIOUS_CUBE) then

      ! subroutines adapted from Deville, Fischer and Mund, High-order methods
      ! for incompressible fluid flow, Cambridge University Press (2002),
      ! pages 386 and 389 and Figure 8.3.1
      do k=1,NGLLZ
        do j=1,NGLLY
! way 1:
!          do i=1,NGLLX
!              iglob = ibool(i,j,k,ispec)
!              dummyx_loc(i,j,k) = displ(1,iglob)
!              dummyy_loc(i,j,k) = displ(2,iglob)
!              dummyz_loc(i,j,k) = displ(3,iglob)
!          enddo

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
              if(NSPEC_INNER_CORE_STRAIN_ONLY == 1) then
                ispec_strain = 1
              else
                ispec_strain = ispec
              endif
              templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
              epsilon_trace_over_3(i,j,k,ispec_strain) = templ
              epsilondev_loc(1,i,j,k) = duxdxl - templ
              epsilondev_loc(2,i,j,k) = duydyl - templ
              epsilondev_loc(3,i,j,k) = 0.5 * duxdyl_plus_duydxl
              epsilondev_loc(4,i,j,k) = 0.5 * duzdxl_plus_duxdzl
              epsilondev_loc(5,i,j,k) = 0.5 * duzdyl_plus_duydzl
            endif

            if(ATTENUATION_VAL) then
              minus_sum_beta =  one_minus_sum_beta(i,j,k,ispec) - 1.0
            endif

            if(ANISOTROPIC_INNER_CORE_VAL) then
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
              c11l = c11store(i,j,k,ispec)
              c12l = c12store(i,j,k,ispec)
              c13l = c13store(i,j,k,ispec)
              c33l = c33store(i,j,k,ispec)
              c44l = c44store(i,j,k,ispec)

              ! use unrelaxed parameters if attenuation
              if(ATTENUATION_VAL) then
                mul = muvstore(i,j,k,ispec)
                c11l = c11l + FOUR_THIRDS * minus_sum_beta * mul
                c12l = c12l - TWO_THIRDS * minus_sum_beta * mul
                c13l = c13l - TWO_THIRDS * minus_sum_beta * mul
                c33l = c33l + FOUR_THIRDS * minus_sum_beta * mul
                c44l = c44l + minus_sum_beta * mul
              endif

              sigma_xx = c11l*duxdxl + c12l*duydyl + c13l*duzdzl
              sigma_yy = c12l*duxdxl + c11l*duydyl + c13l*duzdzl
              sigma_zz = c13l*duxdxl + c13l*duydyl + c33l*duzdzl
              sigma_xy = 0.5*(c11l-c12l)*duxdyl_plus_duydxl
              sigma_xz = c44l*duzdxl_plus_duxdzl
              sigma_yz = c44l*duzdyl_plus_duydzl
            else

              ! inner core with no anisotropy, use kappav and muv for instance
              ! layer with no anisotropy, use kappav and muv for instance
              kappal = kappavstore(i,j,k,ispec)
              mul = muvstore(i,j,k,ispec)

              ! use unrelaxed parameters if attenuation
              if(ATTENUATION_VAL) then
                mul = mul * one_minus_sum_beta(i,j,k,ispec)
              endif

              lambdalplus2mul = kappal + FOUR_THIRDS * mul
              lambdal = lambdalplus2mul - 2.*mul

              ! compute stress sigma
              sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
              sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
              sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl

              sigma_xy = mul*duxdyl_plus_duydxl
              sigma_xz = mul*duzdxl_plus_duxdzl
              sigma_yz = mul*duzdyl_plus_duydzl

            endif

            ! subtract memory variables if attenuation
            if(ATTENUATION_VAL .and. ( USE_ATTENUATION_MIMIC .eqv. .false. ) ) then

! way 1:
!              do i_SLS = 1,N_SLS
!                R_xx_val = R_memory(1,i_SLS,i,j,k,ispec)
!                R_yy_val = R_memory(2,i_SLS,i,j,k,ispec)
!                sigma_xx = sigma_xx - R_xx_val
!                sigma_yy = sigma_yy - R_yy_val
!                sigma_zz = sigma_zz + R_xx_val + R_yy_val
!                sigma_xy = sigma_xy - R_memory(3,i_SLS,i,j,k,ispec)
!                sigma_xz = sigma_xz - R_memory(4,i_SLS,i,j,k,ispec)
!                sigma_yz = sigma_yz - R_memory(5,i_SLS,i,j,k,ispec)
!              enddo

! way 2:
! note: this should help compilers to pipeline the code and make better use of the cache;
!          depending on compilers, it can further decrease the computation time by ~ 30%.
!          by default, N_SLS = 3, there for we take steps of 3
            if(imodulo_N_SLS >= 1) then
              do i_SLS = 1,imodulo_N_SLS
                R_xx_val1 = R_memory(1,i_SLS,i,j,k,ispec)
                R_yy_val1 = R_memory(2,i_SLS,i,j,k,ispec)
                sigma_xx = sigma_xx - R_xx_val1
                sigma_yy = sigma_yy - R_yy_val1
                sigma_zz = sigma_zz + R_xx_val1 + R_yy_val1
                sigma_xy = sigma_xy - R_memory(3,i_SLS,i,j,k,ispec)
                sigma_xz = sigma_xz - R_memory(4,i_SLS,i,j,k,ispec)
                sigma_yz = sigma_yz - R_memory(5,i_SLS,i,j,k,ispec)
              enddo
            endif

            if(N_SLS >= imodulo_N_SLS+1) then
              do i_SLS = imodulo_N_SLS+1,N_SLS,3
                R_xx_val1 = R_memory(1,i_SLS,i,j,k,ispec)
                R_yy_val1 = R_memory(2,i_SLS,i,j,k,ispec)
                sigma_xx = sigma_xx - R_xx_val1
                sigma_yy = sigma_yy - R_yy_val1
                sigma_zz = sigma_zz + R_xx_val1 + R_yy_val1
                sigma_xy = sigma_xy - R_memory(3,i_SLS,i,j,k,ispec)
                sigma_xz = sigma_xz - R_memory(4,i_SLS,i,j,k,ispec)
                sigma_yz = sigma_yz - R_memory(5,i_SLS,i,j,k,ispec)

                R_xx_val2 = R_memory(1,i_SLS+1,i,j,k,ispec)
                R_yy_val2 = R_memory(2,i_SLS+1,i,j,k,ispec)
                sigma_xx = sigma_xx - R_xx_val2
                sigma_yy = sigma_yy - R_yy_val2
                sigma_zz = sigma_zz + R_xx_val2 + R_yy_val2
                sigma_xy = sigma_xy - R_memory(3,i_SLS+1,i,j,k,ispec)
                sigma_xz = sigma_xz - R_memory(4,i_SLS+1,i,j,k,ispec)
                sigma_yz = sigma_yz - R_memory(5,i_SLS+1,i,j,k,ispec)

                R_xx_val3 = R_memory(1,i_SLS+2,i,j,k,ispec)
                R_yy_val3 = R_memory(2,i_SLS+2,i,j,k,ispec)
                sigma_xx = sigma_xx - R_xx_val3
                sigma_yy = sigma_yy - R_yy_val3
                sigma_zz = sigma_zz + R_xx_val3 + R_yy_val3
                sigma_xy = sigma_xy - R_memory(3,i_SLS+2,i,j,k,ispec)
                sigma_xz = sigma_xz - R_memory(4,i_SLS+2,i,j,k,ispec)
                sigma_yz = sigma_yz - R_memory(5,i_SLS+2,i,j,k,ispec)
              enddo
            endif

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
              radius = dble(xstore(iglob1))
              theta = dble(ystore(iglob1))
              phi = dble(zstore(iglob1))

              ! make sure radius is never zero even for points at center of cube
              ! because we later divide by radius
              if(radius < 100.d0 / R_EARTH) radius = 100.d0 / R_EARTH

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

          enddo
        enddo
      enddo

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
          do i=1,NGLLX
            fac1 = wgllwgll_yz(j,k)
            fac2 = wgllwgll_xz(i,k)
            fac3 = wgllwgll_xy(i,j)

            ! sum contributions
            sum_terms(1,i,j,k) = - (fac1*newtempx1(i,j,k) + fac2*newtempx2(i,j,k) + fac3*newtempx3(i,j,k))
            sum_terms(2,i,j,k) = - (fac1*newtempy1(i,j,k) + fac2*newtempy2(i,j,k) + fac3*newtempy3(i,j,k))
            sum_terms(3,i,j,k) = - (fac1*newtempz1(i,j,k) + fac2*newtempz2(i,j,k) + fac3*newtempz3(i,j,k))

            if(GRAVITY_VAL) sum_terms(:,i,j,k) = sum_terms(:,i,j,k) + rho_s_H(:,i,j,k)

          enddo
        enddo
      enddo

      ! sum contributions from each element to the global mesh and add gravity terms
      do k=1,NGLLZ
        do j=1,NGLLY
! way 1:
!          do i=1,NGLLX
!            iglob = ibool(i,j,k,ispec)
!            accel(:,iglob) = accel(:,iglob) + sum_terms(:,i,j,k)
!          enddo

! way 2:
          accel(:,ibool(1,j,k,ispec)) = accel(:,ibool(1,j,k,ispec)) + sum_terms(:,1,j,k)
          accel(:,ibool(2,j,k,ispec)) = accel(:,ibool(2,j,k,ispec)) + sum_terms(:,2,j,k)
          accel(:,ibool(3,j,k,ispec)) = accel(:,ibool(3,j,k,ispec)) + sum_terms(:,3,j,k)
          accel(:,ibool(4,j,k,ispec)) = accel(:,ibool(4,j,k,ispec)) + sum_terms(:,4,j,k)
          accel(:,ibool(5,j,k,ispec)) = accel(:,ibool(5,j,k,ispec)) + sum_terms(:,5,j,k)

        enddo
      enddo

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
      if(ATTENUATION_VAL .and. ( USE_ATTENUATION_MIMIC .eqv. .false. ) ) then
        do i_SLS = 1,N_SLS
          factor_common_use = factor_common(i_SLS,:,:,:,ispec)
          do i_memory = 1,5
             R_memory(i_memory,i_SLS,:,:,:,ispec) = &
                  alphaval(i_SLS) * &
                  R_memory(i_memory,i_SLS,:,:,:,ispec) + muvstore(:,:,:,ispec) * &
                  factor_common_use * &
                  (betaval(i_SLS) * &
                  epsilondev(i_memory,:,:,:,ispec) + gammaval(i_SLS) * epsilondev_loc(i_memory,:,:,:))
          enddo
        enddo

      endif

      ! save deviatoric strain for Runge-Kutta scheme
      if(COMPUTE_AND_STORE_STRAIN) then
! way 1:
        !epsilondev(:,:,:,:,ispec) = epsilondev_loc(:,:,:,:)
! way 2:
        do k=1,NGLLZ
          do j=1,NGLLY
              !dummy(:) = epsilondev_loc(:,1,j,k)

              epsilondev(:,1,j,k,ispec) = epsilondev_loc(:,1,j,k)
              epsilondev(:,2,j,k,ispec) = epsilondev_loc(:,2,j,k)
              epsilondev(:,3,j,k,ispec) = epsilondev_loc(:,3,j,k)
              epsilondev(:,4,j,k,ispec) = epsilondev_loc(:,4,j,k)
              epsilondev(:,5,j,k,ispec) = epsilondev_loc(:,5,j,k)
          enddo
        enddo
      endif

    endif   ! end test to exclude fictitious elements in central cube

  enddo ! spectral element loop

  end subroutine compute_forces_inner_core_Dev

