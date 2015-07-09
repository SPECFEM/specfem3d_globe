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



!--------------------------------------------------------------------------------------------
!
! isotropic element
!
!--------------------------------------------------------------------------------------------

  subroutine compute_element_iso(ispec, &
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

  use constants_solver
  use specfem_par,only: COMPUTE_AND_STORE_STRAIN

  implicit none

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
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION) :: R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: epsilon_trace_over_3

  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec) :: one_minus_sum_beta

  ! gravity
  double precision, dimension(NRAD_GRAVITY) :: minus_gravity_table,density_table,minus_deriv_gravity_table

  ! element info
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NDIM) :: rho_s_H
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5) :: epsilondev_loc

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

  integer :: int_radius
  integer :: iglob

#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
  real(kind=CUSTOM_REAL) :: R_xx_val,R_yy_val
#else
  integer :: i,j,k
#endif

  ! isotropic element

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
    jacobianl = 1.0_CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
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

!ZN beware, here the expression differs from the strain used in memory variable equation (6) in D. Komatitsch and J. Tromp 1999,
!ZN here Brian Savage uses the engineering strain which are epsilon = 1/2*(grad U + (grad U)^T),
!ZN where U is the displacement vector and grad the gradient operator, i.e. there is a 1/2 factor difference between the two.
!ZN Both expressions are fine, but we need to keep in mind that if we put the 1/2 factor here there we need to remove it
!ZN from the expression in which we use the strain later in the code.

    ! compute deviatoric strain
    if (COMPUTE_AND_STORE_STRAIN) then
      templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
      if (NSPEC_CRUST_MANTLE_STRAIN_ONLY == 1) then
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

    !
    ! compute  isotropic  elements
    !

    ! layer with no transverse isotropy, use kappav and muv
    kappal = kappavstore(INDEX_IJK,ispec)
    mul = muvstore(INDEX_IJK,ispec)

    ! use unrelaxed parameters if attenuation
    if (ATTENUATION_VAL) then
      ! precompute terms for attenuation if needed
      if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
        one_minus_sum_beta_use = one_minus_sum_beta(INDEX_IJK,ispec)
      else
        one_minus_sum_beta_use = one_minus_sum_beta(1,1,1,ispec)
      endif
      mul = mul * one_minus_sum_beta_use
    endif

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
!daniel debug: att - debug update
!          call compute_element_att_mem_up_cm(ispec,INDEX_IJK, &
!                                          R_xx(INDEX_IJK,1,ispec), &
!                                          R_yy(INDEX_IJK,1,ispec), &
!                                          R_xy(INDEX_IJK,1,ispec), &
!                                          R_xz(INDEX_IJK,1,ispec), &
!                                          R_yz(INDEX_IJK,1,ispec), &
!                                          epsilondev_loc(INDEX_IJK,:),muvstore(INDEX_IJK,ispec),is_backward_field)

      ! note: function inlining is generally done by fortran compilers;
      !       compilers decide based on performance heuristics
      ! note: fortran passes pointers to array location, thus R_memory(1,1,...) should be fine
      call compute_element_att_stress(i,j,k,R_xx(1,1,1,1,ispec),R_yy(1,1,1,1,ispec),R_xy(1,1,1,1,ispec), &
                                      R_xz(1,1,1,1,ispec),R_yz(1,1,1,1,ispec), &
                                      sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz)
#endif

    endif ! ATTENUATION_VAL

    ! define symmetric components of sigma for gravity
    sigma_yx = sigma_xy
    sigma_zx = sigma_xz
    sigma_zy = sigma_yz

    ! compute non-symmetric terms for gravity
    if (GRAVITY_VAL) then
      ! use mesh coordinates to get theta and phi
      ! x y and z contain r theta and phi
      iglob = ibool(INDEX_IJK,ispec)

      dtheta = dble(ystore(iglob))
      dphi = dble(zstore(iglob))

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
      radius = dble(xstore(iglob))

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

      ! get displacement and multiply by density to compute G tensor
      sx_l = rho * dble(dummyx_loc(INDEX_IJK)) ! dble(displ_crust_mantle(1,iglob))
      sy_l = rho * dble(dummyy_loc(INDEX_IJK)) ! dble(displ_crust_mantle(2,iglob))
      sz_l = rho * dble(dummyz_loc(INDEX_IJK)) ! dble(displ_crust_mantle(3,iglob))

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

  end subroutine compute_element_iso


!--------------------------------------------------------------------------------------------
!
! transversely isotropic element
!
!--------------------------------------------------------------------------------------------

  subroutine compute_element_tiso(ispec, &
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

! this routine is optimized for NGLLX = NGLLY = NGLLZ = 5 using the Deville et al. (2002) inlined matrix-matrix products

  use constants_solver
  use specfem_par,only: COMPUTE_AND_STORE_STRAIN

  implicit none

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
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION) :: R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: epsilon_trace_over_3

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

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NDIM) :: rho_s_H
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5) :: epsilondev_loc

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

  integer :: int_radius
  integer :: iglob

#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
  real(kind=CUSTOM_REAL) :: R_xx_val,R_yy_val
#else
  integer :: i,j,k
#endif

  ! transverse isotropic element

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
    jacobianl = 1.0_CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
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

!ZN beware, here the expression differs from the strain used in memory variable equation (6) in D. Komatitsch and J. Tromp 1999,
!ZN here Brian Savage uses the engineering strain which are epsilon = 1/2*(grad U + (grad U)^T),
!ZN where U is the displacement vector and grad the gradient operator, i.e. there is a 1/2 factor difference between the two.
!ZN Both expressions are fine, but we need to keep in mind that if we put the 1/2 factor here there we need to remove it
!ZN from the expression in which we use the strain later in the code.

    ! compute deviatoric strain
    if (COMPUTE_AND_STORE_STRAIN) then
      templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
      if (NSPEC_CRUST_MANTLE_STRAIN_ONLY == 1) then
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

    !
    ! compute either isotropic or anisotropic elements
    !

! note: the mesh is built such that anisotropic elements are created first in anisotropic layers,
!           thus they are listed first ( see in create_regions_mesh.f90: perm_layer() ordering )
!           this is therefore still in bounds of 1:NSPECMAX_TISO_MANTLE even if NSPECMAX_TISO is less than NSPEC

    ! use kappa and mu from transversely isotropic model
    kappavl = kappavstore(INDEX_IJK,ispec)
    muvl = muvstore(INDEX_IJK,ispec)

    kappahl = kappahstore(INDEX_IJK,ispec)
    muhl = muhstore(INDEX_IJK,ispec)

    ! use unrelaxed parameters if attenuation
    ! eta does not need to be shifted since it is a ratio
    if (ATTENUATION_VAL) then
      ! precompute terms for attenuation if needed
      if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
        one_minus_sum_beta_use = one_minus_sum_beta(INDEX_IJK,ispec)
      else
        one_minus_sum_beta_use = one_minus_sum_beta(1,1,1,ispec)
      endif
      muvl = muvl * one_minus_sum_beta_use
      muhl = muhl * one_minus_sum_beta_use
    endif

    rhovpvsq = kappavl + FOUR_THIRDS * muvl  !!! that is C
    rhovphsq = kappahl + FOUR_THIRDS * muhl  !!! that is A

    rhovsvsq = muvl  !!! that is L
    rhovshsq = muhl  !!! that is N

    eta_aniso = eta_anisostore(INDEX_IJK,ispec)  !!! that is  F / (A - 2 L)

    ! use mesh coordinates to get theta and phi
    ! ystore and zstore contain theta and phi
    iglob = ibool(INDEX_IJK,ispec)

    theta = ystore(iglob)
    phi = zstore(iglob)

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
          (rhovphsq*costhetafour + 2.0_CUSTOM_REAL*costhetasq*sinthetasq*(eta_aniso*rhovphsq &
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
    if (ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL) then
#ifdef FORCE_VECTORIZATION
      ! do NOT put this is a subroutine, otherwise the call to the subroutine prevents compilers
      ! from vectorizing the outer loop

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
      ! note: function inlining is generally done by fortran compilers;
      !       compilers decide based on performance heuristics
      ! note: fortran passes pointers to array location, thus R_memory(1,1,...) is fine
      call compute_element_att_stress(i,j,k,R_xx(1,1,1,1,ispec),R_yy(1,1,1,1,ispec),R_xy(1,1,1,1,ispec), &
                                      R_xz(1,1,1,1,ispec),R_yz(1,1,1,1,ispec), &
                                      sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz)
#endif

    endif ! ATTENUATION_VAL

    ! define symmetric components of sigma for gravity
    sigma_yx = sigma_xy
    sigma_zx = sigma_xz
    sigma_zy = sigma_yz

    ! compute non-symmetric terms for gravity
    if (GRAVITY_VAL) then
      ! use mesh coordinates to get theta and phi
      ! x y and z contain r theta and phi
      iglob = ibool(INDEX_IJK,ispec)

      dtheta = dble(ystore(iglob))
      dphi = dble(zstore(iglob))
      radius = dble(xstore(iglob))

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

      ! get displacement and multiply by density to compute G tensor
      sx_l = rho * dble(dummyx_loc(INDEX_IJK))
      sy_l = rho * dble(dummyy_loc(INDEX_IJK))
      sz_l = rho * dble(dummyz_loc(INDEX_IJK))

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

  end subroutine compute_element_tiso


!--------------------------------------------------------------------------------------------
!
! anisotropic element
!
!--------------------------------------------------------------------------------------------

  subroutine compute_element_aniso(ispec, &
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

  use constants_solver
  use specfem_par,only: COMPUTE_AND_STORE_STRAIN

  implicit none

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
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION) :: R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: epsilon_trace_over_3

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

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NDIM) :: rho_s_H
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5) :: epsilondev_loc

  ! local parameters
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

  integer :: int_radius
  integer :: iglob

#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
  real(kind=CUSTOM_REAL) :: R_xx_val,R_yy_val
#else
  integer :: i,j,k
#endif

  !  anisotropic elements

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
    jacobianl = 1.0_CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
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

!ZN beware, here the expression differs from the strain used in memory variable equation (6) in D. Komatitsch and J. Tromp 1999,
!ZN here Brian Savage uses the engineering strain which are epsilon = 1/2*(grad U + (grad U)^T),
!ZN where U is the displacement vector and grad the gradient operator, i.e. there is a 1/2 factor difference between the two.
!ZN Both expressions are fine, but we need to keep in mind that if we put the 1/2 factor here there we need to remove it
!ZN from the expression in which we use the strain later in the code.

    ! compute deviatoric strain
    if (COMPUTE_AND_STORE_STRAIN) then
      templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
      if (NSPEC_CRUST_MANTLE_STRAIN_ONLY == 1) then
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

    !
    ! compute anisotropic elements
    !

    c11 = c11store(INDEX_IJK,ispec)
    c12 = c12store(INDEX_IJK,ispec)
    c13 = c13store(INDEX_IJK,ispec)
    c14 = c14store(INDEX_IJK,ispec)
    c15 = c15store(INDEX_IJK,ispec)
    c16 = c16store(INDEX_IJK,ispec)
    c22 = c22store(INDEX_IJK,ispec)
    c23 = c23store(INDEX_IJK,ispec)
    c24 = c24store(INDEX_IJK,ispec)
    c25 = c25store(INDEX_IJK,ispec)
    c26 = c26store(INDEX_IJK,ispec)
    c33 = c33store(INDEX_IJK,ispec)
    c34 = c34store(INDEX_IJK,ispec)
    c35 = c35store(INDEX_IJK,ispec)
    c36 = c36store(INDEX_IJK,ispec)
    c44 = c44store(INDEX_IJK,ispec)
    c45 = c45store(INDEX_IJK,ispec)
    c46 = c46store(INDEX_IJK,ispec)
    c55 = c55store(INDEX_IJK,ispec)
    c56 = c56store(INDEX_IJK,ispec)
    c66 = c66store(INDEX_IJK,ispec)

    if (ATTENUATION_VAL) then
      ! precompute terms for attenuation if needed
      if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
        minus_sum_beta =  one_minus_sum_beta(INDEX_IJK,ispec) - 1.0_CUSTOM_REAL
      else
        minus_sum_beta =  one_minus_sum_beta(1,1,1,ispec) - 1.0_CUSTOM_REAL
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
    if (ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL) then
#ifdef FORCE_VECTORIZATION
      ! do NOT put this is a subroutine, otherwise the call to the subroutine prevents compilers
      ! from vectorizing the outer loop

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
      ! note: function inlining is generally done by fortran compilers;
      !       compilers decide based on performance heuristics

      ! note: Fortran passes pointers to array location, thus R_memory(1,1,...) is fine
      call compute_element_att_stress(i,j,k,R_xx(1,1,1,1,ispec),R_yy(1,1,1,1,ispec),R_xy(1,1,1,1,ispec), &
                                      R_xz(1,1,1,1,ispec),R_yz(1,1,1,1,ispec), &
                                      sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz)
#endif
    endif ! ATTENUATION_VAL

    ! define symmetric components of sigma for gravity
    sigma_yx = sigma_xy
    sigma_zx = sigma_xz
    sigma_zy = sigma_yz

    ! compute non-symmetric terms for gravity
    if (GRAVITY_VAL) then
      ! use mesh coordinates to get theta and phi
      ! x y and z contain r theta and phi
      iglob = ibool(INDEX_IJK,ispec)

      dtheta = dble(ystore(iglob))
      dphi = dble(zstore(iglob))
      radius = dble(xstore(iglob))

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

      ! get displacement and multiply by density to compute G tensor
      sx_l = rho * dble(dummyx_loc(INDEX_IJK)) ! dble(displ_crust_mantle(1,iglob))
      sy_l = rho * dble(dummyy_loc(INDEX_IJK)) ! dble(displ_crust_mantle(2,iglob))
      sz_l = rho * dble(dummyz_loc(INDEX_IJK)) !  dble(displ_crust_mantle(3,iglob))

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

  end subroutine compute_element_aniso



!--------------------------------------------------------------------------------------------
!
! helper functions
!
!--------------------------------------------------------------------------------------------
!
! please leave this routine in this file, to help compilers inlining this function...
!


  subroutine compute_element_att_stress(i,j,k,R_xx_loc,R_yy_loc,R_xy_loc,R_xz_loc,R_yz_loc, &
                                       sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz)

! updates stress with attenuation correction


  use constants_solver,only: CUSTOM_REAL,N_SLS,NGLLX,NGLLY,NGLLZ

  implicit none

  integer, intent(in) :: i,j,k
  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS) :: R_xx_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS) :: R_yy_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS) :: R_xy_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS) :: R_xz_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS) :: R_yz_loc

  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz

! local parameters
  real(kind=CUSTOM_REAL) R_xx_val,R_yy_val
  integer :: i_SLS

  do i_SLS = 1,N_SLS
    R_xx_val = R_xx_loc(i,j,k,i_SLS)
    R_yy_val = R_yy_loc(i,j,k,i_SLS)
    sigma_xx = sigma_xx - R_xx_val
    sigma_yy = sigma_yy - R_yy_val
    sigma_zz = sigma_zz + R_xx_val + R_yy_val
    sigma_xy = sigma_xy - R_xy_loc(i,j,k,i_SLS)
    sigma_xz = sigma_xz - R_xz_loc(i,j,k,i_SLS)
    sigma_yz = sigma_yz - R_yz_loc(i,j,k,i_SLS)
  enddo

  end subroutine compute_element_att_stress

! for reference...
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
!  use constants_solver,only: CUSTOM_REAL
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
!  use constants_solver,only: CUSTOM_REAL
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
!  use constants_solver,only: CUSTOM_REAL
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
!  use constants_solver,only: CUSTOM_REAL
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
!  use constants_solver,only: CUSTOM_REAL
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
!  use constants_solver,only: CUSTOM_REAL
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
!  use constants_solver,only: CUSTOM_REAL
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
!  use constants_solver,only: CUSTOM_REAL
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
!  use constants_solver,only: CUSTOM_REAL
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
!  use constants_solver,only: CUSTOM_REAL
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
!  use constants_solver,only: CUSTOM_REAL
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
!  use constants_solver,only: CUSTOM_REAL
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
!  use constants_solver,only: CUSTOM_REAL
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
!  use constants_solver,only: CUSTOM_REAL
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


