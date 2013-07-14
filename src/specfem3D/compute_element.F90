!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            April 2011
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

  subroutine compute_element_iso(ispec, &
                    minus_gravity_table,density_table,minus_deriv_gravity_table, &
                    xstore,ystore,zstore, &
                    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                    wgll_cube, &
                    kappavstore,muvstore, &
                    ibool, &
                    R_memory, &
                    one_minus_sum_beta,vnspec, &
                    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                    dummyx_loc,dummyy_loc,dummyz_loc,epsilondev_loc,eps_trace_over_3_loc,rho_s_H,PARTIAL_PHYS_DISPERSION_ONLY)

! this routine is optimized for NGLLX = NGLLY = NGLLZ = 5 using the Deville et al. (2002) inlined matrix-matrix products

  implicit none

  include "constants.h"

  ! include values created by the mesher
  ! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

  ! element id
  integer :: ispec

  ! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool

  ! x y and z contain r theta and phi
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  ! array with derivatives of Lagrange polynomials and precalculated products
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

  ! store anisotropic properties only where needed to save memory
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ISO_MANTLE) :: &
        kappavstore,muvstore

  ! variable sized array variables
  integer :: vnspec

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: R_memory

  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec) :: one_minus_sum_beta

  ! gravity
  double precision, dimension(NRAD_GRAVITY) :: minus_gravity_table,density_table,minus_deriv_gravity_table

  ! element info
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: rho_s_H
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: eps_trace_over_3_loc

! local parameters
  real(kind=CUSTOM_REAL) one_minus_sum_beta_use
  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl
  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl
  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz
  real(kind=CUSTOM_REAL) sigma_yx,sigma_zx,sigma_zy

  real(kind=CUSTOM_REAL) templ
  real(kind=CUSTOM_REAL) lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) kappal

  ! for gravity
  double precision dphi,dtheta
  double precision radius,rho,minus_g,minus_dg
  double precision minus_g_over_radius,minus_dg_plus_g_over_radius
  double precision cos_theta,sin_theta,cos_phi,sin_phi
  double precision cos_theta_sq,sin_theta_sq,cos_phi_sq,sin_phi_sq
  double precision factor,sx_l,sy_l,sz_l,gxl,gyl,gzl
  double precision Hxxl,Hyyl,Hzzl,Hxyl,Hxzl,Hyzl

  integer :: i,j,k
  integer :: int_radius
  integer :: iglob1

  logical :: PARTIAL_PHYS_DISPERSION_ONLY

  ! isotropic element

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
        jacobianl = 1.0_CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
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
!ZN beware, here the expression differs from the strain used in memory variable equation (6) in D. Komatitsch and J. Tromp 1999,
!ZN here Brian Savage uses the engineering strain which are epsilon = 1/2*(grad U + (grad U)^T),
!ZN where U is the displacement vector and grad the gradient operator, i.e. there is a 1/2 factor difference between the two.
!ZN Both expressions are fine, but we need to keep in mind that if we put the 1/2 factor here there we need to remove it
!ZN from the expression in which we use the strain later in the code.
          templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
          eps_trace_over_3_loc(i,j,k) = templ
          epsilondev_loc(1,i,j,k) = duxdxl - templ
          epsilondev_loc(2,i,j,k) = duydyl - templ
          epsilondev_loc(3,i,j,k) = 0.5 * duxdyl_plus_duydxl
          epsilondev_loc(4,i,j,k) = 0.5 * duzdxl_plus_duxdzl
          epsilondev_loc(5,i,j,k) = 0.5 * duzdyl_plus_duydzl
        endif

       ! precompute terms for attenuation if needed
        if(ATTENUATION_3D_VAL) then
          one_minus_sum_beta_use = one_minus_sum_beta(i,j,k,ispec)
        else if(ATTENUATION_VAL) then
          one_minus_sum_beta_use = one_minus_sum_beta(1,1,1,ispec)
        endif

        !
        ! compute  isotropic  elements
        !

        ! layer with no transverse isotropy, use kappav and muv
        kappal = kappavstore(i,j,k,ispec)
        mul = muvstore(i,j,k,ispec)

        ! use unrelaxed parameters if attenuation
        if(ATTENUATION_VAL) mul = mul * one_minus_sum_beta_use

        lambdalplus2mul = kappal + FOUR_THIRDS * mul
        lambdal = lambdalplus2mul - 2.0_CUSTOM_REAL*mul

        ! compute stress sigma
        sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
        sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
        sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl

        sigma_xy = mul*duxdyl_plus_duydxl
        sigma_xz = mul*duzdxl_plus_duxdzl
        sigma_yz = mul*duzdyl_plus_duydzl

        ! subtract memory variables if attenuation
        if(ATTENUATION_VAL .and. ( PARTIAL_PHYS_DISPERSION_ONLY .eqv. .false. )  ) then

          ! note: fortran passes pointers to array location, thus R_memory(1,1,...) should be fine
          call compute_element_att_stress( R_memory(1,1,i,j,k,ispec), &
                    sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz)

        endif ! ATTENUATION_VAL

        ! define symmetric components of sigma for gravity
        sigma_yx = sigma_xy
        sigma_zx = sigma_xz
        sigma_zy = sigma_yz

        ! compute non-symmetric terms for gravity
        if(GRAVITY_VAL) then

            ! use mesh coordinates to get theta and phi
            ! x y and z contain r theta and phi
            iglob1 = ibool(i,j,k,ispec)

            dtheta = dble(ystore(iglob1))
            dphi = dble(zstore(iglob1))

            cos_theta = dcos(dtheta)
            sin_theta = dsin(dtheta)
            cos_phi = dcos(dphi)
            sin_phi = dsin(dphi)

            cos_theta_sq = cos_theta*cos_theta
            sin_theta_sq = sin_theta*sin_theta
            cos_phi_sq = cos_phi*cos_phi
            sin_phi_sq = sin_phi*sin_phi

            ! get g, rho and dg/dr=dg
            ! spherical components of the gravitational acceleration
            ! for efficiency replace with lookup table every 100 m in radial direction
            radius = dble(xstore(iglob1))

            int_radius = nint(10.d0 * radius * R_EARTH_KM )
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

            ! distinguish between single and double precision for reals
            if(CUSTOM_REAL == SIZE_REAL) then

              ! get displacement and multiply by density to compute G tensor
              sx_l = rho * dble(dummyx_loc(i,j,k)) ! dble(displ_crust_mantle(1,iglob1))
              sy_l = rho * dble(dummyy_loc(i,j,k)) ! dble(displ_crust_mantle(2,iglob1))
              sz_l = rho * dble(dummyz_loc(i,j,k)) ! dble(displ_crust_mantle(3,iglob1))

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
              sx_l = rho * dummyx_loc(i,j,k) ! displ_crust_mantle(1,iglob1)
              sy_l = rho * dummyy_loc(i,j,k) ! displ_crust_mantle(2,iglob1)
              sz_l = rho * dummyz_loc(i,j,k) ! displ_crust_mantle(3,iglob1)

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
        tempx1(i,j,k) = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl) ! this goes to accel_x
        tempy1(i,j,k) = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl) ! this goes to accel_y
        tempz1(i,j,k) = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl) ! this goes to accel_z

        tempx2(i,j,k) = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl) ! this goes to accel_x
        tempy2(i,j,k) = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl) ! this goes to accel_y
        tempz2(i,j,k) = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl) ! this goes to accel_z

        tempx3(i,j,k) = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl) ! this goes to accel_x
        tempy3(i,j,k) = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl) ! this goes to accel_y
        tempz3(i,j,k) = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl) ! this goes to accel_z

      enddo ! NGLLX
    enddo ! NGLLY
  enddo ! NGLLZ

  end subroutine compute_element_iso

!
!--------------------------------------------------------------------------------------------------
!

  subroutine compute_element_tiso(ispec, &
                    minus_gravity_table,density_table,minus_deriv_gravity_table, &
                    xstore,ystore,zstore, &
                    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                    wgll_cube, &
                    kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                    ibool, &
                    R_memory, &
                    one_minus_sum_beta,vnspec, &
                    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                    dummyx_loc,dummyy_loc,dummyz_loc,epsilondev_loc,eps_trace_over_3_loc,rho_s_H,PARTIAL_PHYS_DISPERSION_ONLY)

! this routine is optimized for NGLLX = NGLLY = NGLLZ = 5 using the Deville et al. (2002) inlined matrix-matrix products

  implicit none

  include "constants.h"

  ! include values created by the mesher
  ! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

  ! element id
  integer :: ispec

  ! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool

  ! x y and z contain r theta and phi
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  ! array with derivatives of Lagrange polynomials and precalculated products
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

  ! store anisotropic properties only where needed to save memory
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE) :: &
        kappahstore,muhstore,eta_anisostore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ISO_MANTLE) :: &
        kappavstore,muvstore

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: R_memory

  ! variable sized array variables
  integer :: vnspec

  ! [alpha,beta,gamma]val reduced to N_SLS  to N_SLS*NUM_NODES
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec) :: one_minus_sum_beta

  ! gravity
  double precision, dimension(NRAD_GRAVITY) :: minus_gravity_table,density_table,minus_deriv_gravity_table

  ! element info
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: rho_s_H
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: eps_trace_over_3_loc

! local parameters
  real(kind=CUSTOM_REAL) one_minus_sum_beta_use
  ! the 21 coefficients for an anisotropic medium in reduced notation
  real(kind=CUSTOM_REAL) c11,c22,c33,c44,c55,c66,c12,c13,c23,c14,c24,c34,c15,c25,c35,c45,c16,c26,c36,c46,c56

  real(kind=CUSTOM_REAL) rhovphsq,sinphifour,cosphisq,sinphisq,costhetasq,rhovsvsq,sinthetasq, &
        cosphifour,costhetafour,rhovpvsq,sinthetafour,rhovshsq,cosfourphi, &
        costwotheta,cosfourtheta,sintwophisq,costheta,sinphi,sintheta,cosphi, &
        sintwotheta,costwophi,sintwophi,costwothetasq,costwophisq,phi,theta

  real(kind=CUSTOM_REAL) two_rhovsvsq,two_rhovshsq ! two_rhovpvsq,two_rhovphsq
  real(kind=CUSTOM_REAL) four_rhovsvsq,four_rhovshsq ! four_rhovpvsq,four_rhovphsq

  real(kind=CUSTOM_REAL) twoetaminone,etaminone,eta_aniso
  real(kind=CUSTOM_REAL) two_eta_aniso,four_eta_aniso,six_eta_aniso
  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl
  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl
  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz

  real(kind=CUSTOM_REAL) templ
  real(kind=CUSTOM_REAL) templ1,templ1_cos,templ2,templ2_cos,templ3,templ3_two,templ3_cos
  real(kind=CUSTOM_REAL) kappavl,kappahl,muvl,muhl

  ! for gravity
  double precision dphi,dtheta
  double precision radius,rho,minus_g,minus_dg
  double precision minus_g_over_radius,minus_dg_plus_g_over_radius
  double precision cos_theta,sin_theta,cos_phi,sin_phi
  double precision cos_theta_sq,sin_theta_sq,cos_phi_sq,sin_phi_sq
  double precision factor,sx_l,sy_l,sz_l,gxl,gyl,gzl
  double precision Hxxl,Hyyl,Hzzl,Hxyl,Hxzl,Hyzl
  real(kind=CUSTOM_REAL) sigma_yx,sigma_zx,sigma_zy

  integer :: i,j,k
  integer :: int_radius
  integer :: iglob1

  logical :: PARTIAL_PHYS_DISPERSION_ONLY

  ! transverse isotropic element

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
        jacobianl = 1.0_CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
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
!ZN beware, here the expression differs from the strain used in memory variable equation (6) in D. Komatitsch and J. Tromp 1999,
!ZN here Brian Savage uses the engineering strain which are epsilon = 1/2*(grad U + (grad U)^T),
!ZN where U is the displacement vector and grad the gradient operator, i.e. there is a 1/2 factor difference between the two.
!ZN Both expressions are fine, but we need to keep in mind that if we put the 1/2 factor here there we need to remove it
!ZN from the expression in which we use the strain later in the code.
          templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
          eps_trace_over_3_loc(i,j,k) = templ
          epsilondev_loc(1,i,j,k) = duxdxl - templ
          epsilondev_loc(2,i,j,k) = duydyl - templ
          epsilondev_loc(3,i,j,k) = 0.5 * duxdyl_plus_duydxl
          epsilondev_loc(4,i,j,k) = 0.5 * duzdxl_plus_duxdzl
          epsilondev_loc(5,i,j,k) = 0.5 * duzdyl_plus_duydzl
        endif

        ! precompute terms for attenuation if needed
        if(ATTENUATION_3D_VAL) then
          one_minus_sum_beta_use = one_minus_sum_beta(i,j,k,ispec)
        else if(ATTENUATION_VAL) then
          one_minus_sum_beta_use = one_minus_sum_beta(1,1,1,ispec)
        endif

        !
        ! compute either isotropic or anisotropic elements
        !

! note : mesh is built such that anisotropic elements are created first in anisotropic layers,
!           thus they are listed first ( see in create_regions_mesh.f90: perm_layer() ordering )
!           this is therefore still in bounds of 1:NSPECMAX_TISO_MANTLE even if NSPECMAX_TISO is less than NSPEC

        ! uncomment to debug
        !if ( ispec > NSPECMAX_TISO_MANTLE ) then
        !  print*,'error tiso: ispec = ',ispec,'max = ',NSPECMAX_TISO_MANTLE
        !  call exit_mpi(0,'error tiso ispec bounds')
        !endif

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


        ! way 2: pre-compute temporary values
        templ1 = four_rhovsvsq - rhovpvsq + twoetaminone*rhovphsq - four_eta_aniso*rhovsvsq
        templ1_cos = rhovphsq - rhovpvsq + costwotheta*templ1
        templ2 = four_rhovsvsq - rhovpvsq - rhovphsq + two_eta_aniso*rhovphsq - four_eta_aniso*rhovsvsq
        templ2_cos = rhovpvsq - rhovphsq + costwotheta*templ2
        templ3 = rhovphsq + rhovpvsq - two_eta_aniso*rhovphsq + four_eta_aniso*rhovsvsq
        templ3_two = templ3 - two_rhovshsq - two_rhovsvsq
        templ3_cos = templ3_two + costwotheta*templ2

        ! way 2: reordering operations to facilitate compilation, avoiding divisions, using locality for temporary values
        c11 = rhovphsq*sinphifour &
              + 2.0_CUSTOM_REAL*cosphisq*sinphisq* &
              ( rhovphsq*costhetasq + sinthetasq*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq) ) &
              + cosphifour*(rhovphsq*costhetafour &
                + 2.0_CUSTOM_REAL*costhetasq*sinthetasq*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq) &
                + rhovpvsq*sinthetafour)

        c12 = 0.25_CUSTOM_REAL*costhetasq*(rhovphsq - two_rhovshsq)*(3.0_CUSTOM_REAL + cosfourphi) &
              - four_rhovshsq*cosphisq*costhetasq*sinphisq &
              + 0.03125_CUSTOM_REAL*rhovphsq*sintwophisq*(11.0_CUSTOM_REAL + cosfourtheta + 4.0*costwotheta) &
              + eta_aniso*sinthetasq*(rhovphsq - two_rhovsvsq) &
                         *(cosphifour + sinphifour + 2.0_CUSTOM_REAL*cosphisq*costhetasq*sinphisq) &
              + rhovpvsq*cosphisq*sinphisq*sinthetafour &
              - rhovsvsq*sintwophisq*sinthetafour

        c13 = 0.125_CUSTOM_REAL*cosphisq*(rhovphsq + six_eta_aniso*rhovphsq + rhovpvsq - four_rhovsvsq &
                    - 12.0_CUSTOM_REAL*eta_aniso*rhovsvsq + cosfourtheta*templ1) &
              + sinphisq*(eta_aniso*costhetasq*(rhovphsq - two_rhovsvsq) + sinthetasq*(rhovphsq - two_rhovshsq))

        ! uses temporary templ1 from c13
        c15 = cosphi*costheta*sintheta* &
              ( 0.5_CUSTOM_REAL*cosphisq* (rhovpvsq - rhovphsq + costwotheta*templ1) &
                + etaminone*sinphisq*(rhovphsq - two_rhovsvsq))

        c14 = costheta*sinphi*sintheta* &
              ( 0.5_CUSTOM_REAL*cosphisq*(templ2_cos + four_rhovshsq - four_rhovsvsq) &
                + sinphisq*(etaminone*rhovphsq + 2.0_CUSTOM_REAL*(rhovshsq - eta_aniso*rhovsvsq)) )

        ! uses temporary templ2_cos from c14
        c16 = 0.5_CUSTOM_REAL*cosphi*sinphi*sinthetasq* &
              ( cosphisq*templ2_cos &
                + 2.0_CUSTOM_REAL*etaminone*sinphisq*(rhovphsq - two_rhovsvsq) )

        c22 = rhovphsq*cosphifour + 2.0_CUSTOM_REAL*cosphisq*sinphisq* &
              (rhovphsq*costhetasq + sinthetasq*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)) &
              + sinphifour* &
              (rhovphsq*costhetafour + 2.0_CUSTOM_REAL*costhetasq*sinthetasq*(eta_aniso*rhovphsq  &
                    + two_rhovsvsq - two_eta_aniso*rhovsvsq) + rhovpvsq*sinthetafour)

        ! uses temporary templ1 from c13
        c23 = 0.125_CUSTOM_REAL*sinphisq*(rhovphsq + six_eta_aniso*rhovphsq &
                + rhovpvsq - four_rhovsvsq - 12.0_CUSTOM_REAL*eta_aniso*rhovsvsq + cosfourtheta*templ1) &
              + cosphisq*(eta_aniso*costhetasq*(rhovphsq - two_rhovsvsq) + sinthetasq*(rhovphsq - two_rhovshsq))

        ! uses temporary templ1 from c13
        c24 = costheta*sinphi*sintheta* &
              ( etaminone*cosphisq*(rhovphsq - two_rhovsvsq) &
                + 0.5_CUSTOM_REAL*sinphisq*(rhovpvsq - rhovphsq + costwotheta*templ1) )

        ! uses temporary templ2_cos from c14
        c25 = cosphi*costheta*sintheta* &
              ( cosphisq*(etaminone*rhovphsq + 2.0_CUSTOM_REAL*(rhovshsq - eta_aniso*rhovsvsq)) &
                + 0.5_CUSTOM_REAL*sinphisq*(templ2_cos + four_rhovshsq - four_rhovsvsq) )

        ! uses temporary templ2_cos from c14
        c26 = 0.5_CUSTOM_REAL*cosphi*sinphi*sinthetasq* &
              ( 2.0_CUSTOM_REAL*etaminone*cosphisq*(rhovphsq - two_rhovsvsq) &
                + sinphisq*templ2_cos )

        c33 = rhovpvsq*costhetafour &
              + 2.0_CUSTOM_REAL*costhetasq*sinthetasq*(two_rhovsvsq + eta_aniso*(rhovphsq - two_rhovsvsq)) &
              + rhovphsq*sinthetafour

        ! uses temporary templ1_cos from c13
        c34 = - 0.25_CUSTOM_REAL*sinphi*sintwotheta*templ1_cos

        ! uses temporary templ1_cos from c34
        c35 = - 0.25_CUSTOM_REAL*cosphi*sintwotheta*templ1_cos

        ! uses temporary templ1_cos from c34
        c36 = - 0.25_CUSTOM_REAL*sintwophi*sinthetasq*(templ1_cos - four_rhovshsq + four_rhovsvsq)

        c44 = cosphisq*(rhovsvsq*costhetasq + rhovshsq*sinthetasq) &
              + sinphisq*(rhovsvsq*costwothetasq + costhetasq*sinthetasq*templ3)

        ! uses temporary templ3 from c44
        c46 = - cosphi*costheta*sintheta* &
                ( cosphisq*(rhovshsq - rhovsvsq) - 0.5_CUSTOM_REAL*sinphisq*templ3_cos  )

        ! uses templ3 from c46
        c45 = 0.25_CUSTOM_REAL*sintwophi*sinthetasq* &
              (templ3_two + costwotheta*(rhovphsq + rhovpvsq - two_eta_aniso*rhovphsq + 4.0_CUSTOM_REAL*etaminone*rhovsvsq))

        c55 = sinphisq*(rhovsvsq*costhetasq + rhovshsq*sinthetasq) &
              + cosphisq*(rhovsvsq*costwothetasq &
                  + costhetasq*sinthetasq*(rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + four_eta_aniso*rhovsvsq) )

        ! uses temporary templ3_cos from c46
        c56 = costheta*sinphi*sintheta* &
              ( 0.5_CUSTOM_REAL*cosphisq*templ3_cos + sinphisq*(rhovsvsq - rhovshsq) )

        c66 = rhovshsq*costwophisq*costhetasq &
              - 2.0_CUSTOM_REAL*cosphisq*costhetasq*sinphisq*(rhovphsq - two_rhovshsq) &
              + 0.03125_CUSTOM_REAL*rhovphsq*sintwophisq*(11.0_CUSTOM_REAL + 4.0_CUSTOM_REAL*costwotheta + cosfourtheta) &
              - 0.125_CUSTOM_REAL*rhovsvsq*sinthetasq* &
              ( -6.0_CUSTOM_REAL - 2.0_CUSTOM_REAL*costwotheta - 2.0_CUSTOM_REAL*cosfourphi &
                        + cos(4.0_CUSTOM_REAL*phi - 2.0_CUSTOM_REAL*theta) &
                        + cos(2.0_CUSTOM_REAL*(2.0_CUSTOM_REAL*phi + theta)) ) &
              + rhovpvsq*cosphisq*sinphisq*sinthetafour &
              - 0.5_CUSTOM_REAL*eta_aniso*sintwophisq*sinthetafour*(rhovphsq - two_rhovsvsq)

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

        ! subtract memory variables if attenuation
        if(ATTENUATION_VAL .and. ( PARTIAL_PHYS_DISPERSION_ONLY .eqv. .false. )  ) then

          ! note: fortran passes pointers to array location, thus R_memory(1,1,...) should be fine
          call compute_element_att_stress( R_memory(1,1,i,j,k,ispec), &
                    sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz)

        endif ! ATTENUATION_VAL

        ! define symmetric components of sigma for gravity
        sigma_yx = sigma_xy
        sigma_zx = sigma_xz
        sigma_zy = sigma_yz

        ! compute non-symmetric terms for gravity
        if(GRAVITY_VAL) then

             ! use mesh coordinates to get theta and phi
            ! x y and z contain r theta and phi
            iglob1 = ibool(i,j,k,ispec)

            dtheta = dble(ystore(iglob1))
            dphi = dble(zstore(iglob1))
            radius = dble(xstore(iglob1))

            cos_theta = dcos(dtheta)
            sin_theta = dsin(dtheta)
            cos_phi = dcos(dphi)
            sin_phi = dsin(dphi)

            ! way 2
            cos_theta_sq = cos_theta*cos_theta
            sin_theta_sq = sin_theta*sin_theta
            cos_phi_sq = cos_phi*cos_phi
            sin_phi_sq = sin_phi*sin_phi

            ! get g, rho and dg/dr=dg
            ! spherical components of the gravitational acceleration
            ! for efficiency replace with lookup table every 100 m in radial direction
            int_radius = nint(10.d0 * radius * R_EARTH_KM )
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

            ! distinguish between single and double precision for reals
            if(CUSTOM_REAL == SIZE_REAL) then

              ! get displacement and multiply by density to compute G tensor
              sx_l = rho * dble(dummyx_loc(i,j,k))
              sy_l = rho * dble(dummyy_loc(i,j,k))
              sz_l = rho * dble(dummyz_loc(i,j,k))


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
              sx_l = rho * dummyx_loc(i,j,k)
              sy_l = rho * dummyy_loc(i,j,k)
              sz_l = rho * dummyz_loc(i,j,k)

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
        tempx1(i,j,k) = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl) ! this goes to accel_x
        tempy1(i,j,k) = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl) ! this goes to accel_y
        tempz1(i,j,k) = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl) ! this goes to accel_z

        tempx2(i,j,k) = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl) ! this goes to accel_x
        tempy2(i,j,k) = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl) ! this goes to accel_y
        tempz2(i,j,k) = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl) ! this goes to accel_z

        tempx3(i,j,k) = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl) ! this goes to accel_x
        tempy3(i,j,k) = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl) ! this goes to accel_y
        tempz3(i,j,k) = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl) ! this goes to accel_z

      enddo ! NGLLX
    enddo ! NGLLY
  enddo ! NGLLZ

  end subroutine compute_element_tiso

!
!--------------------------------------------------------------------------------------------
!

  subroutine compute_element_aniso(ispec, &
                    minus_gravity_table,density_table,minus_deriv_gravity_table, &
                    xstore,ystore,zstore, &
                    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                    wgll_cube, &
                    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                    ibool, &
                    R_memory, &
                    one_minus_sum_beta,vnspec, &
                    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                    dummyx_loc,dummyy_loc,dummyz_loc,epsilondev_loc,eps_trace_over_3_loc,rho_s_H,PARTIAL_PHYS_DISPERSION_ONLY)

! this routine is optimized for NGLLX = NGLLY = NGLLZ = 5 using the Deville et al. (2002) inlined matrix-matrix products

  implicit none

  include "constants.h"

  ! include values created by the mesher
  ! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

  ! element id
  integer :: ispec

  ! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool

  ! x y and z contain r theta and phi
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  ! array with derivatives of Lagrange polynomials and precalculated products
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

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

  ! variable sized array variables
  integer :: vnspec

  ! [alpha,beta,gamma]val reduced to N_SLS  to N_SLS*NUM_NODES
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec) :: one_minus_sum_beta

  ! gravity
  double precision, dimension(NRAD_GRAVITY) :: minus_gravity_table,density_table,minus_deriv_gravity_table

  ! element info
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: rho_s_H
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: eps_trace_over_3_loc

! local parameters
  real(kind=CUSTOM_REAL) one_minus_sum_beta_use
  real(kind=CUSTOM_REAL) minus_sum_beta,mul
  ! the 21 coefficients for an anisotropic medium in reduced notation
  real(kind=CUSTOM_REAL) c11,c22,c33,c44,c55,c66,c12,c13,c23,c14,c24,c34,c15,c25,c35,c45,c16,c26,c36,c46,c56

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl
  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl
  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz

  real(kind=CUSTOM_REAL) templ

  ! for gravity
  double precision dphi,dtheta
  double precision radius,rho,minus_g,minus_dg
  double precision minus_g_over_radius,minus_dg_plus_g_over_radius
  double precision cos_theta,sin_theta,cos_phi,sin_phi
  double precision cos_theta_sq,sin_theta_sq,cos_phi_sq,sin_phi_sq
  double precision factor,sx_l,sy_l,sz_l,gxl,gyl,gzl
  double precision Hxxl,Hyyl,Hzzl,Hxyl,Hxzl,Hyzl
  real(kind=CUSTOM_REAL) sigma_yx,sigma_zx,sigma_zy

  integer :: i,j,k
  integer :: int_radius
  integer :: iglob1

  logical :: PARTIAL_PHYS_DISPERSION_ONLY

  !  anisotropic elements

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
        jacobianl = 1.0_CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
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
!ZN beware, here the expression differs from the strain used in memory variable equation (6) in D. Komatitsch and J. Tromp 1999,
!ZN here Brian Savage uses the engineering strain which are epsilon = 1/2*(grad U + (grad U)^T),
!ZN where U is the displacement vector and grad the gradient operator, i.e. there is a 1/2 factor difference between the two.
!ZN Both expressions are fine, but we need to keep in mind that if we put the 1/2 factor here there we need to remove it
!ZN from the expression in which we use the strain later in the code.
          templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
          eps_trace_over_3_loc(i,j,k) = templ
          epsilondev_loc(1,i,j,k) = duxdxl - templ
          epsilondev_loc(2,i,j,k) = duydyl - templ
          epsilondev_loc(3,i,j,k) = 0.5 * duxdyl_plus_duydxl
          epsilondev_loc(4,i,j,k) = 0.5 * duzdxl_plus_duxdzl
          epsilondev_loc(5,i,j,k) = 0.5 * duzdyl_plus_duydzl
        endif

        ! precompute terms for attenuation if needed
        if(ATTENUATION_3D_VAL) then
          one_minus_sum_beta_use = one_minus_sum_beta(i,j,k,ispec)
          minus_sum_beta =  one_minus_sum_beta_use - 1.0
        else if(ATTENUATION_VAL) then
          one_minus_sum_beta_use = one_minus_sum_beta(1,1,1,ispec)
          minus_sum_beta =  one_minus_sum_beta_use - 1.0
        endif

        !
        ! compute anisotropic elements
        !

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
          !mul = c44
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

        ! subtract memory variables if attenuation
        if(ATTENUATION_VAL .and. ( PARTIAL_PHYS_DISPERSION_ONLY .eqv. .false. )  ) then

          ! note: fortran passes pointers to array location, thus R_memory(1,1,...) should be fine
          call compute_element_att_stress(R_memory(1,1,i,j,k,ispec), &
                                          sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz)

        endif ! ATTENUATION_VAL

        ! define symmetric components of sigma for gravity
        sigma_yx = sigma_xy
        sigma_zx = sigma_xz
        sigma_zy = sigma_yz

        ! compute non-symmetric terms for gravity
        if(GRAVITY_VAL) then

             ! use mesh coordinates to get theta and phi
            ! x y and z contain r theta and phi
            iglob1 = ibool(i,j,k,ispec)

            dtheta = dble(ystore(iglob1))
            dphi = dble(zstore(iglob1))
            radius = dble(xstore(iglob1))

            cos_theta = dcos(dtheta)
            sin_theta = dsin(dtheta)
            cos_phi = dcos(dphi)
            sin_phi = dsin(dphi)

            cos_theta_sq = cos_theta*cos_theta
            sin_theta_sq = sin_theta*sin_theta
            cos_phi_sq = cos_phi*cos_phi
            sin_phi_sq = sin_phi*sin_phi

            ! get g, rho and dg/dr=dg
            ! spherical components of the gravitational acceleration
            ! for efficiency replace with lookup table every 100 m in radial direction
            int_radius = nint(10.d0 * radius * R_EARTH_KM )
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

            ! distinguish between single and double precision for reals
            if(CUSTOM_REAL == SIZE_REAL) then

              ! get displacement and multiply by density to compute G tensor
              sx_l = rho * dble(dummyx_loc(i,j,k)) ! dble(displ_crust_mantle(1,iglob1))
              sy_l = rho * dble(dummyy_loc(i,j,k)) ! dble(displ_crust_mantle(2,iglob1))
              sz_l = rho * dble(dummyz_loc(i,j,k)) !  dble(displ_crust_mantle(3,iglob1))

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
              sx_l = rho * dummyx_loc(i,j,k)  ! displ_crust_mantle(1,iglob1)
              sy_l = rho * dummyy_loc(i,j,k)  !  displ_crust_mantle(2,iglob1)
              sz_l = rho * dummyz_loc(i,j,k)  ! displ_crust_mantle(3,iglob1)

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
        tempx1(i,j,k) = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl) ! this goes to accel_x
        tempy1(i,j,k) = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl) ! this goes to accel_y
        tempz1(i,j,k) = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl) ! this goes to accel_z

        tempx2(i,j,k) = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl) ! this goes to accel_x
        tempy2(i,j,k) = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl) ! this goes to accel_y
        tempz2(i,j,k) = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl) ! this goes to accel_z

        tempx3(i,j,k) = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl) ! this goes to accel_x
        tempy3(i,j,k) = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl) ! this goes to accel_y
        tempz3(i,j,k) = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl) ! this goes to accel_z

      enddo ! NGLLX
    enddo ! NGLLY
  enddo ! NGLLZ

  end subroutine compute_element_aniso

!
!--------------------------------------------------------------------------------------------
!


  subroutine compute_element_att_stress(R_memory_loc, &
                                       sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz)

  implicit none

  include "constants.h"

  ! include values created by the mesher
  ! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(5,N_SLS) :: R_memory_loc
  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz

! local parameters
  real(kind=CUSTOM_REAL) R_xx_val1,R_yy_val1
  integer :: i_SLS

  do i_SLS = 1,N_SLS
    R_xx_val1 = R_memory_loc(1,i_SLS)
    R_yy_val1 = R_memory_loc(2,i_SLS)
    sigma_xx = sigma_xx - R_xx_val1
    sigma_yy = sigma_yy - R_yy_val1
    sigma_zz = sigma_zz + R_xx_val1 + R_yy_val1
    sigma_xy = sigma_xy - R_memory_loc(3,i_SLS)
    sigma_xz = sigma_xz - R_memory_loc(4,i_SLS)
    sigma_yz = sigma_yz - R_memory_loc(5,i_SLS)
  enddo

  end subroutine compute_element_att_stress

!
!--------------------------------------------------------------------------------------------
!

  subroutine compute_element_att_memory_cr(ispec,R_memory, &
                                        vnspec,factor_common, &
                                        alphaval,betaval,gammaval, &
                                        c44store,muvstore, &
                                        epsilondev_loc_nplus1,epsilondev_loc,&
                                        istage,R_memory_lddrk,tau_sigma_CUSTOM_REAL,deltat,USE_LDDRK)
! crust mantle
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

  implicit none

  include "constants.h"

  ! include values created by the mesher
  ! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

  ! element id
  integer :: ispec

  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: R_memory

  ! variable sized array variables
  integer :: vnspec

  real(kind=CUSTOM_REAL), dimension(N_SLS,ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec) :: factor_common
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE) :: c44store
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ISO_MANTLE) :: muvstore

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc

! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: factor_common_c44_muv
  integer :: i_SLS

  integer :: i_memory,i,j,k

! for LDDRK
  integer :: istage
  logical :: USE_LDDRK
  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: R_memory_lddrk
  real(kind=CUSTOM_REAL),dimension(N_SLS) :: tau_sigma_CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: deltat

  ! use Runge-Kutta scheme to march in time

  ! get coefficients for that standard linear solid
  ! IMPROVE we use mu_v here even if there is some anisotropy
  ! IMPROVE we should probably use an average value instead

  do i_SLS = 1,N_SLS

    ! reformatted R_memory to handle large factor_common and reduced [alpha,beta,gamma]val

    if(ATTENUATION_3D_VAL) then
      do k = 1,NGLLZ
        do j = 1,NGLLZ
          do i = 1,NGLLZ
            factor_common_c44_muv(i,j,k) = factor_common(i_SLS,i,j,k,ispec)
          enddo
        enddo
      enddo
    else
      do k = 1,NGLLZ
        do j = 1,NGLLZ
          do i = 1,NGLLZ
            factor_common_c44_muv(i,j,k) = factor_common(i_SLS,1,1,1,ispec)
          enddo
        enddo
      enddo
    endif

    if(ANISOTROPIC_3D_MANTLE_VAL) then
      factor_common_c44_muv(:,:,:) = factor_common_c44_muv(:,:,:) * c44store(:,:,:,ispec)
    else
      factor_common_c44_muv(:,:,:) = factor_common_c44_muv(:,:,:) * muvstore(:,:,:,ispec)
    endif

    if(USE_LDDRK)then
      do i_memory = 1,5
        R_memory_lddrk(i_memory,i_SLS,:,:,:,ispec) = ALPHA_LDDRK(istage) * R_memory_lddrk(i_memory,i_SLS,:,:,:,ispec) + &
                      deltat * (factor_common_c44_muv(:,:,:)*epsilondev_loc(i_memory,:,:,:) - &
                                R_memory(i_memory,i_SLS,:,:,:,ispec)*(1._CUSTOM_REAL/tau_sigma_CUSTOM_REAL(i_SLS)))
        R_memory(i_memory,i_SLS,:,:,:,ispec) = R_memory(i_memory,i_SLS,:,:,:,ispec) + &
                                               BETA_LDDRK(istage) * R_memory_lddrk(i_memory,i_SLS,:,:,:,ispec)
      enddo
    else
      do i_memory = 1,5
        R_memory(i_memory,i_SLS,:,:,:,ispec) = alphaval(i_SLS) * R_memory(i_memory,i_SLS,:,:,:,ispec) &
                  + factor_common_c44_muv(:,:,:) &
                  * (betaval(i_SLS) * epsilondev_loc(i_memory,:,:,:) + gammaval(i_SLS) * epsilondev_loc_nplus1(i_memory,:,:,:))
      enddo
    endif
  enddo ! i_SLS

  end subroutine compute_element_att_memory_cr

!
!--------------------------------------------------------------------------------------------
!

  subroutine compute_element_att_memory_ic(ispec,R_memory, &
                                        vnspec,factor_common, &
                                        alphaval,betaval,gammaval, &
                                        muvstore, &
                                        epsilondev_loc_nplus1,epsilondev_loc,&
                                        istage,R_memory_lddrk,tau_sigma_CUSTOM_REAL,deltat,USE_LDDRK)
! inner core
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

  implicit none

  include "constants.h"

  ! include values created by the mesher
  ! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

  ! element id
  integer :: ispec

  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ATTENUATION) :: R_memory

  ! variable sized array variables
  integer :: vnspec

  real(kind=CUSTOM_REAL), dimension(N_SLS,ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec) :: factor_common
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: muvstore

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc_nplus1

! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: factor_common_use

  integer :: i_SLS

  integer :: i_memory,i,j,k

! for LDDRK
  integer :: istage
  logical :: USE_LDDRK
  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: R_memory_lddrk
  real(kind=CUSTOM_REAL),dimension(N_SLS) :: tau_sigma_CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: deltat

  ! use Runge-Kutta scheme to march in time

  ! get coefficients for that standard linear solid
  ! IMPROVE we use mu_v here even if there is some anisotropy
  ! IMPROVE we should probably use an average value instead

  do i_SLS = 1,N_SLS

    if(ATTENUATION_3D_VAL) then
      do k = 1,NGLLZ
        do j = 1,NGLLZ
          do i = 1,NGLLZ
            factor_common_use(i,j,k) = factor_common(i_SLS,i,j,k,ispec)
          enddo
        enddo
      enddo
    else
      do k = 1,NGLLZ
        do j = 1,NGLLZ
          do i = 1,NGLLZ
            factor_common_use(i,j,k) = factor_common(i_SLS,1,1,1,ispec)
          enddo
        enddo
      enddo
    endif

    if(USE_LDDRK)then
      do i_memory = 1,5
        R_memory_lddrk(i_memory,i_SLS,:,:,:,ispec) = ALPHA_LDDRK(istage) * R_memory_lddrk(i_memory,i_SLS,:,:,:,ispec) + &
            deltat * (muvstore(:,:,:,ispec) * factor_common_use(:,:,:)*epsilondev_loc(i_memory,:,:,:) - &
                      R_memory(i_memory,i_SLS,:,:,:,ispec)*(1._CUSTOM_REAL/tau_sigma_CUSTOM_REAL(i_SLS)))
        R_memory(i_memory,i_SLS,:,:,:,ispec) = R_memory(i_memory,i_SLS,:,:,:,ispec) + &
                                               BETA_LDDRK(istage) * R_memory_lddrk(i_memory,i_SLS,:,:,:,ispec)
      enddo
    else
      do i_memory = 1,5
         R_memory(i_memory,i_SLS,:,:,:,ispec) = alphaval(i_SLS) * R_memory(i_memory,i_SLS,:,:,:,ispec) &
              + muvstore(:,:,:,ispec) * factor_common_use(:,:,:) * &
              (betaval(i_SLS) * epsilondev_loc_nplus1(i_memory,:,:,:) + gammaval(i_SLS) * epsilondev_loc(i_memory,:,:,:))
      enddo
    endif

  enddo

  end subroutine compute_element_att_memory_ic


!
!--------------------------------------------------------------------------------------------
!
 subroutine compute_element_strain_undo_att_Dev(ispec,nglob,nspec,displ,ibool,hprime_xx,hprime_xxT,&
                                       xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,epsilondev_loc,eps_trace_over_3_loc)

  implicit none
  include "constants.h"

  integer :: ispec,nglob,nspec
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob) :: displ
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprime_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: eps_trace_over_3_loc

!  local variable
  integer :: i,j,k,iglob
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,m2) :: B1_m1_m2_5points,B2_m1_m2_5points,B3_m1_m2_5points
  equivalence(dummyx_loc,B1_m1_m2_5points)
  equivalence(dummyy_loc,B2_m1_m2_5points)
  equivalence(dummyz_loc,B3_m1_m2_5points)
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: C1_m1_m2_5points,C2_m1_m2_5points,C3_m1_m2_5points
  equivalence(tempx1,C1_m1_m2_5points)
  equivalence(tempy1,C2_m1_m2_5points)
  equivalence(tempz1,C3_m1_m2_5points)
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: &
    C1_mxm_m2_m1_5points,C2_mxm_m2_m1_5points,C3_mxm_m2_m1_5points
  equivalence(tempx3,C1_mxm_m2_m1_5points)
  equivalence(tempy3,C2_mxm_m2_m1_5points)
  equivalence(tempz3,C3_mxm_m2_m1_5points)
  real(kind=CUSTOM_REAL), dimension(m2,NGLLX) :: &
    A1_mxm_m2_m1_5points,A2_mxm_m2_m1_5points,A3_mxm_m2_m1_5points
  equivalence(dummyx_loc,A1_mxm_m2_m1_5points)
  equivalence(dummyy_loc,A2_mxm_m2_m1_5points)
  equivalence(dummyz_loc,A3_mxm_m2_m1_5points)

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duydyl,duzdzl,duxdyl,duydxl,duzdxl,duxdzl,duzdyl,duydzl,&
                         duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl,&
                         duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)
            dummyx_loc(i,j,k) = displ(1,iglob)
            dummyy_loc(i,j,k) = displ(2,iglob)
            dummyz_loc(i,j,k) = displ(3,iglob)
        enddo
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

        eps_trace_over_3_loc(i,j,k) = ONE_THIRD * (duxdxl + duydyl + duzdzl)
        epsilondev_loc(1,i,j,k) = duxdxl - eps_trace_over_3_loc(i,j,k)
        epsilondev_loc(2,i,j,k) = duydyl - eps_trace_over_3_loc(i,j,k)
        epsilondev_loc(3,i,j,k) = 0.5 * duxdyl_plus_duydxl
        epsilondev_loc(4,i,j,k) = 0.5 * duzdxl_plus_duxdzl
        epsilondev_loc(5,i,j,k) = 0.5 * duzdyl_plus_duydzl

      enddo
    enddo
  enddo

  end subroutine compute_element_strain_undo_att_Dev

!
!--------------------------------------------------------------------------------------------
!
 subroutine compute_element_strain_att_Dev(ispec,nglob,nspec,displ,veloc,deltat,ibool,hprime_xx,hprime_xxT,&
                                       xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,epsilondev_loc_nplus1,&
                                       eps_trace_over_3_loc_nplus1)

  implicit none
  include "constants.h"

  integer :: ispec,nglob,nspec
  real(kind=CUSTOM_REAL) :: deltat
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob) :: displ,veloc
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprime_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: eps_trace_over_3_loc_nplus1

!  local variable
  integer :: i,j,k,iglob
  real(kind=CUSTOM_REAL) :: templ
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,m2) :: B1_m1_m2_5points,B2_m1_m2_5points,B3_m1_m2_5points
  equivalence(dummyx_loc,B1_m1_m2_5points)
  equivalence(dummyy_loc,B2_m1_m2_5points)
  equivalence(dummyz_loc,B3_m1_m2_5points)
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: C1_m1_m2_5points,C2_m1_m2_5points,C3_m1_m2_5points
  equivalence(tempx1,C1_m1_m2_5points)
  equivalence(tempy1,C2_m1_m2_5points)
  equivalence(tempz1,C3_m1_m2_5points)
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: &
    C1_mxm_m2_m1_5points,C2_mxm_m2_m1_5points,C3_mxm_m2_m1_5points
  equivalence(tempx3,C1_mxm_m2_m1_5points)
  equivalence(tempy3,C2_mxm_m2_m1_5points)
  equivalence(tempz3,C3_mxm_m2_m1_5points)
  real(kind=CUSTOM_REAL), dimension(m2,NGLLX) :: &
    A1_mxm_m2_m1_5points,A2_mxm_m2_m1_5points,A3_mxm_m2_m1_5points
  equivalence(dummyx_loc,A1_mxm_m2_m1_5points)
  equivalence(dummyy_loc,A2_mxm_m2_m1_5points)
  equivalence(dummyz_loc,A3_mxm_m2_m1_5points)

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duydyl,duzdzl,duxdyl,duydxl,duzdxl,duxdzl,duzdyl,duydzl,&
                         duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl,&
                         duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)
            dummyx_loc(i,j,k) = displ(1,iglob) + deltat * veloc(1,iglob)
            dummyy_loc(i,j,k) = displ(2,iglob) + deltat * veloc(2,iglob)
            dummyz_loc(i,j,k) = displ(3,iglob) + deltat * veloc(3,iglob)
        enddo
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

        templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
        eps_trace_over_3_loc_nplus1 = templ
        epsilondev_loc_nplus1(1,i,j,k) = duxdxl - templ
        epsilondev_loc_nplus1(2,i,j,k) = duydyl - templ
        epsilondev_loc_nplus1(3,i,j,k) = 0.5 * duxdyl_plus_duydxl
        epsilondev_loc_nplus1(4,i,j,k) = 0.5 * duzdxl_plus_duxdzl
        epsilondev_loc_nplus1(5,i,j,k) = 0.5 * duzdyl_plus_duydzl

      enddo
    enddo
  enddo

 end subroutine compute_element_strain_att_Dev

!=====================================================================

  subroutine compute_element_strain_undo_att_noDev(ispec,nglob,nspec,displ,hprime_xx,hprime_yy,hprime_zz,ibool,&
          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,epsilondev_loc,eps_trace_over_3_loc)

  implicit none

  include "constants.h"

  ! include values created by the mesher
  ! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

  ! element id
  integer :: ispec,i,j,k

  integer NSPEC,NGLOB

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: displ
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: eps_trace_over_3_loc

! local parameters
  integer iglob
  integer l

  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) tempz1l,tempz2l,tempz3l

  real(kind=CUSTOM_REAL) hp1,hp2,hp3

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          tempy1l = 0._CUSTOM_REAL
          tempy2l = 0._CUSTOM_REAL
          tempy3l = 0._CUSTOM_REAL

          tempz1l = 0._CUSTOM_REAL
          tempz2l = 0._CUSTOM_REAL
          tempz3l = 0._CUSTOM_REAL

          do l=1,NGLLX
            hp1 = hprime_xx(i,l)
            iglob = ibool(l,j,k,ispec)
            tempx1l = tempx1l + displ(1,iglob)*hp1
            tempy1l = tempy1l + displ(2,iglob)*hp1
            tempz1l = tempz1l + displ(3,iglob)*hp1
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
            hp2 = hprime_yy(j,l)
            iglob = ibool(i,l,k,ispec)
            tempx2l = tempx2l + displ(1,iglob)*hp2
            tempy2l = tempy2l + displ(2,iglob)*hp2
            tempz2l = tempz2l + displ(3,iglob)*hp2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
            hp3 = hprime_zz(k,l)
            iglob = ibool(i,j,l,ispec)
            tempx3l = tempx3l + displ(1,iglob)*hp3
            tempy3l = tempy3l + displ(2,iglob)*hp3
            tempz3l = tempz3l + displ(3,iglob)*hp3
          enddo

!         get derivatives of ux, uy and uz with respect to x, y and z

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

          duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
          duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
          duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

          duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
          duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
          duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l

! precompute some sums to save CPU time
          duxdxl_plus_duydyl = duxdxl + duydyl
          duxdxl_plus_duzdzl = duxdxl + duzdzl
          duydyl_plus_duzdzl = duydyl + duzdzl
          duxdyl_plus_duydxl = duxdyl + duydxl
          duzdxl_plus_duxdzl = duzdxl + duxdzl
          duzdyl_plus_duydzl = duzdyl + duydzl

! compute deviatoric strain
          eps_trace_over_3_loc(i,j,k) = ONE_THIRD * (duxdxl + duydyl + duzdzl)
          epsilondev_loc(1,i,j,k) = duxdxl - eps_trace_over_3_loc(i,j,k)
          epsilondev_loc(2,i,j,k) = duydyl - eps_trace_over_3_loc(i,j,k)
          epsilondev_loc(3,i,j,k) = 0.5 * duxdyl_plus_duydxl
          epsilondev_loc(4,i,j,k) = 0.5 * duzdxl_plus_duxdzl
          epsilondev_loc(5,i,j,k) = 0.5 * duzdyl_plus_duydzl

        enddo ! NGLLX
      enddo ! NGLLY
    enddo ! NGLLZ

  end subroutine compute_element_strain_undo_att_noDev

!=====================================================================

  subroutine compute_element_strain_att_noDev(ispec,nglob,nspec,displ,veloc,deltat,hprime_xx,hprime_yy,hprime_zz,ibool,&
          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,epsilondev_loc_nplus1)

  implicit none

  include "constants.h"

! include values created by the mesher
! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

  integer ispec,NSPEC,NGLOB
  real(kind=CUSTOM_REAL) deltat

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz


! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: displ,veloc

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc_nplus1

!local parameters
  integer iglob
  integer i,j,k,l
  real(kind=CUSTOM_REAL) templ

  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) tempz1l,tempz2l,tempz3l

  real(kind=CUSTOM_REAL) hp1,hp2,hp3

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          tempy1l = 0._CUSTOM_REAL
          tempy2l = 0._CUSTOM_REAL
          tempy3l = 0._CUSTOM_REAL

          tempz1l = 0._CUSTOM_REAL
          tempz2l = 0._CUSTOM_REAL
          tempz3l = 0._CUSTOM_REAL

          do l=1,NGLLX
            hp1 = hprime_xx(i,l)
            iglob = ibool(l,j,k,ispec)
            tempx1l = tempx1l + (displ(1,iglob) + deltat * veloc(1,iglob))*hp1
            tempy1l = tempy1l + (displ(2,iglob) + deltat * veloc(2,iglob))*hp1
            tempz1l = tempz1l + (displ(3,iglob) + deltat * veloc(3,iglob))*hp1
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
            hp2 = hprime_yy(j,l)
            iglob = ibool(i,l,k,ispec)
            tempx2l = tempx2l + (displ(1,iglob) + deltat * veloc(1,iglob))*hp2
            tempy2l = tempy2l + (displ(2,iglob) + deltat * veloc(2,iglob))*hp2
            tempz2l = tempz2l + (displ(3,iglob) + deltat * veloc(3,iglob))*hp2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
            hp3 = hprime_zz(k,l)
            iglob = ibool(i,j,l,ispec)
            tempx3l = tempx3l + (displ(1,iglob) + deltat * veloc(1,iglob))*hp3
            tempy3l = tempy3l + (displ(2,iglob) + deltat * veloc(2,iglob))*hp3
            tempz3l = tempz3l + (displ(3,iglob) + deltat * veloc(3,iglob))*hp3
          enddo

!         get derivatives of ux, uy and uz with respect to x, y and z

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

          duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
          duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
          duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

          duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
          duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
          duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l

! precompute some sums to save CPU time
          duxdxl_plus_duydyl = duxdxl + duydyl
          duxdxl_plus_duzdzl = duxdxl + duzdzl
          duydyl_plus_duzdzl = duydyl + duzdzl
          duxdyl_plus_duydxl = duxdyl + duydxl
          duzdxl_plus_duxdzl = duzdxl + duxdzl
          duzdyl_plus_duydzl = duzdyl + duydzl

! compute deviatoric strain
          templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
          epsilondev_loc_nplus1(1,i,j,k) = duxdxl - templ
          epsilondev_loc_nplus1(2,i,j,k) = duydyl - templ
          epsilondev_loc_nplus1(3,i,j,k) = 0.5 * duxdyl_plus_duydxl
          epsilondev_loc_nplus1(4,i,j,k) = 0.5 * duzdxl_plus_duxdzl
          epsilondev_loc_nplus1(5,i,j,k) = 0.5 * duzdyl_plus_duydzl

        enddo ! NGLLX
      enddo ! NGLLY
    enddo ! NGLLZ

  end subroutine compute_element_strain_att_noDev

