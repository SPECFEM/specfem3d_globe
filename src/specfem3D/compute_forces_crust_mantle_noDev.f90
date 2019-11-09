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

  subroutine compute_forces_crust_mantle_noDev(NSPEC_STR_OR_ATT,NGLOB,NSPEC_ATT, &
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
                                               factor_common,vnspec)

  use constants_solver

  use specfem_par, only: &
    hprime_xx,hprime_yy,hprime_zz,hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
    wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
    gravity_pre_store => gravity_pre_store_crust_mantle,gravity_H => gravity_H_crust_mantle, &
    COMPUTE_AND_STORE_STRAIN,USE_LDDRK

  use specfem_par_crustmantle, only: &
    xix => xix_crust_mantle,xiy => xiy_crust_mantle,xiz => xiz_crust_mantle, &
    etax => etax_crust_mantle,etay => etay_crust_mantle,etaz => etaz_crust_mantle, &
    gammax => gammax_crust_mantle,gammay => gammay_crust_mantle,gammaz => gammaz_crust_mantle, &
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

  ! inner/outer element run flag
  integer,intent(in) :: iphase

  ! local parameters

  ! for attenuation
  real(kind=CUSTOM_REAL) R_xx_val,R_yy_val
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5) :: epsilondev_loc

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  integer ispec,iglob
  integer i,j,k,l
  integer i_SLS

  ! the 21 coefficients for an anisotropic medium in reduced notation
  real(kind=CUSTOM_REAL) c11,c22,c33,c44,c55,c66,c12,c13,c23,c14,c24,c34,c15,c25,c35,c45,c16,c26,c36,c46,c56

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz

  real(kind=CUSTOM_REAL) hp1,hp2,hp3
  real(kind=CUSTOM_REAL) fac1,fac2,fac3
  real(kind=CUSTOM_REAL) lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) kappal
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sum_terms

  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) tempz1l,tempz2l,tempz3l
  real(kind=CUSTOM_REAL) templ

  ! for gravity
  real(kind=CUSTOM_REAL) sigma_yx,sigma_zx,sigma_zy
  real(kind=CUSTOM_REAL) factor,sx_l,sy_l,sz_l,gxl,gyl,gzl
  real(kind=CUSTOM_REAL) Hxxl,Hyyl,Hzzl,Hxyl,Hxzl,Hyzl
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: rho_s_H

!  integer :: computed_elements
  integer :: num_elements,ispec_p

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

  do ispec_p = 1,num_elements

    ispec = phase_ispec_inner(ispec_p,iphase)

    ! only compute element which belong to current phase (inner or outer elements)

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          tempy1l = 0._CUSTOM_REAL
          tempy2l = 0._CUSTOM_REAL
          tempy3l = 0._CUSTOM_REAL

          tempz1l = 0._CUSTOM_REAL
          tempz2l = 0._CUSTOM_REAL
          tempz3l = 0._CUSTOM_REAL

          do l = 1,NGLLX
            hp1 = hprime_xx(i,l)
            iglob = ibool(l,j,k,ispec)
            tempx1l = tempx1l + displ_crust_mantle(1,iglob)*hp1
            tempy1l = tempy1l + displ_crust_mantle(2,iglob)*hp1
            tempz1l = tempz1l + displ_crust_mantle(3,iglob)*hp1
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l = 1,NGLLY
            hp2 = hprime_yy(j,l)
            iglob = ibool(i,l,k,ispec)
            tempx2l = tempx2l + displ_crust_mantle(1,iglob)*hp2
            tempy2l = tempy2l + displ_crust_mantle(2,iglob)*hp2
            tempz2l = tempz2l + displ_crust_mantle(3,iglob)*hp2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l = 1,NGLLZ
            hp3 = hprime_zz(k,l)
            iglob = ibool(i,j,l,ispec)
            tempx3l = tempx3l + displ_crust_mantle(1,iglob)*hp3
            tempy3l = tempy3l + displ_crust_mantle(2,iglob)*hp3
            tempz3l = tempz3l + displ_crust_mantle(3,iglob)*hp3
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

! compute the Jacobian
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
          if (COMPUTE_AND_STORE_STRAIN) then
            templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
            if (NSPEC_CRUST_MANTLE_STRAIN_ONLY == 1) then
              if (ispec == 1) then
                epsilon_trace_over_3(i,j,k,1) = templ
              endif
            else
              epsilon_trace_over_3(i,j,k,ispec) = templ
            endif
            epsilondev_loc(i,j,k,1) = duxdxl - templ
            epsilondev_loc(i,j,k,2) = duydyl - templ
            epsilondev_loc(i,j,k,3) = 0.5 * duxdyl_plus_duydxl
            epsilondev_loc(i,j,k,4) = 0.5 * duzdxl_plus_duxdzl
            epsilondev_loc(i,j,k,5) = 0.5 * duzdyl_plus_duydzl
          endif

          ! precompute terms for attenuation if needed
          !if (ATTENUATION_VAL) then
          !  if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
          !    one_minus_sum_beta_use = one_minus_sum_beta(i,j,k,ispec)
          !  else
          !    one_minus_sum_beta_use = one_minus_sum_beta(1,1,1,ispec)
          !  endif
          !  minus_sum_beta =  one_minus_sum_beta_use - 1.0_CUSTOM_REAL
          !endif

          !
          ! compute either isotropic or anisotropic elements
          !

          if (ANISOTROPIC_3D_MANTLE_VAL) then

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

            ! already done in prepare_timerun...
            !if (ATTENUATION_VAL) then
            !  mul = c44
            !  c11 = c11 + FOUR_THIRDS * minus_sum_beta * mul
            !  c12 = c12 - TWO_THIRDS * minus_sum_beta * mul
            !  c13 = c13 - TWO_THIRDS * minus_sum_beta * mul
            !  c22 = c22 + FOUR_THIRDS * minus_sum_beta * mul
            !  c23 = c23 - TWO_THIRDS * minus_sum_beta * mul
            !  c33 = c33 + FOUR_THIRDS * minus_sum_beta * mul
            !  c44 = c44 + minus_sum_beta * mul
            !  c55 = c55 + minus_sum_beta * mul
            !  c66 = c66 + minus_sum_beta * mul
            !endif

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
            if (.not. ispec_is_tiso(ispec)) then

              ! isotropic element

              ! layer with no transverse isotropy, use kappav and muv
              kappal = kappavstore(i,j,k,ispec)
              mul = muvstore(i,j,k,ispec)

              ! use unrelaxed parameters if attenuation
              ! already done in prepare_timerun...
              !if (ATTENUATION_VAL) mul = mul * one_minus_sum_beta_use

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

              ! transverse isotropic element
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
          if (ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL) then
             do i_SLS = 1,N_SLS
                R_xx_val = R_xx(i,j,k,i_SLS,ispec)
                R_yy_val = R_yy(i,j,k,i_SLS,ispec)
                sigma_xx = sigma_xx - R_xx_val
                sigma_yy = sigma_yy - R_yy_val
                sigma_zz = sigma_zz + R_xx_val + R_yy_val
                sigma_xy = sigma_xy - R_xy(i,j,k,i_SLS,ispec)
                sigma_xz = sigma_xz - R_xz(i,j,k,i_SLS,ispec)
                sigma_yz = sigma_yz - R_yz(i,j,k,i_SLS,ispec)
             enddo
          endif

          ! define symmetric components of sigma for gravity
          sigma_yx = sigma_xy
          sigma_zx = sigma_xz
          sigma_zy = sigma_yz

          ! compute non-symmetric terms for gravity
          if (GRAVITY_VAL) then
            ! use mesh coordinates to get theta and phi
            ! x y and z contain r theta and phi
            iglob = ibool(i,j,k,ispec)

            ! Cartesian components of the gravitational acceleration
            gxl = gravity_pre_store(1,iglob) ! minus_g*sin_theta*cos_phi * rho
            gyl = gravity_pre_store(2,iglob) ! minus_g*sin_theta*sin_phi * rho
            gzl = gravity_pre_store(3,iglob) ! minus_g*cos_theta * rho

            ! Cartesian components of gradient of gravitational acceleration
            ! get displacement and multiply by density to compute G tensor
            sx_l = displ_crust_mantle(1,iglob)
            sy_l = displ_crust_mantle(2,iglob)
            sz_l = displ_crust_mantle(3,iglob)

            ! compute G tensor from s . g and add to sigma (not symmetric)
            sigma_xx = sigma_xx + sy_l * gyl + sz_l * gzl
            sigma_yy = sigma_yy + sx_l * gxl + sz_l * gzl
            sigma_zz = sigma_zz + sx_l * gxl + sy_l * gyl

            sigma_xy = sigma_xy - sx_l * gyl
            sigma_yx = sigma_yx - sy_l * gxl

            sigma_xz = sigma_xz - sx_l * gzl
            sigma_zx = sigma_zx - sz_l * gxl

            sigma_yz = sigma_yz - sy_l * gzl
            sigma_zy = sigma_zy - sz_l * gyl

            Hxxl = gravity_H(1,iglob) ! minus_g_over_radius*(cos_phi_sq*cos_theta_sq + sin_phi_sq)
                                      !   + cos_phi_sq*minus_dg*sin_theta_sq * rho
            Hyyl = gravity_H(2,iglob) ! minus_g_over_radius*(cos_phi_sq + cos_theta_sq*sin_phi_sq)
                                      !   + minus_dg*sin_phi_sq*sin_theta_sq * rho
            Hzzl = gravity_H(3,iglob) ! cos_theta_sq*minus_dg + minus_g_over_radius*sin_theta_sq * rho
            Hxyl = gravity_H(4,iglob) ! cos_phi*minus_dg_plus_g_over_radius*sin_phi*sin_theta_sq * rho
            Hxzl = gravity_H(5,iglob) ! cos_phi*cos_theta*minus_dg_plus_g_over_radius*sin_theta * rho
            Hyzl = gravity_H(6,iglob) ! cos_theta*minus_dg_plus_g_over_radius*sin_phi*sin_theta * rho

            ! precompute vector
            factor = jacobianl * wgll_cube(i,j,k)
            rho_s_H(1,i,j,k) = factor * (sx_l * Hxxl + sy_l * Hxyl + sz_l * Hxzl)
            rho_s_H(2,i,j,k) = factor * (sx_l * Hxyl + sy_l * Hyyl + sz_l * Hyzl)
            rho_s_H(3,i,j,k) = factor * (sx_l * Hxzl + sy_l * Hyzl + sz_l * Hzzl)
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

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempy1l = 0._CUSTOM_REAL
          tempz1l = 0._CUSTOM_REAL

          tempx2l = 0._CUSTOM_REAL
          tempy2l = 0._CUSTOM_REAL
          tempz2l = 0._CUSTOM_REAL

          tempx3l = 0._CUSTOM_REAL
          tempy3l = 0._CUSTOM_REAL
          tempz3l = 0._CUSTOM_REAL

          do l = 1,NGLLX
            fac1 = hprimewgll_xx(l,i)
            tempx1l = tempx1l + tempx1(l,j,k)*fac1
            tempy1l = tempy1l + tempy1(l,j,k)*fac1
            tempz1l = tempz1l + tempz1(l,j,k)*fac1
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l = 1,NGLLY
            fac2 = hprimewgll_yy(l,j)
            tempx2l = tempx2l + tempx2(i,l,k)*fac2
            tempy2l = tempy2l + tempy2(i,l,k)*fac2
            tempz2l = tempz2l + tempz2(i,l,k)*fac2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l = 1,NGLLZ
            fac3 = hprimewgll_zz(l,k)
            tempx3l = tempx3l + tempx3(i,j,l)*fac3
            tempy3l = tempy3l + tempy3(i,j,l)*fac3
            tempz3l = tempz3l + tempz3(i,j,l)*fac3
          enddo

          fac1 = wgllwgll_yz(j,k)
          fac2 = wgllwgll_xz(i,k)
          fac3 = wgllwgll_xy(i,j)

          sum_terms(1,i,j,k) = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l)
          sum_terms(2,i,j,k) = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l)
          sum_terms(3,i,j,k) = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l)

          if (GRAVITY_VAL) then
            sum_terms(1,i,j,k) = sum_terms(1,i,j,k) + rho_s_H(1,i,j,k)
            sum_terms(2,i,j,k) = sum_terms(2,i,j,k) + rho_s_H(2,i,j,k)
            sum_terms(3,i,j,k) = sum_terms(3,i,j,k) + rho_s_H(3,i,j,k)
          endif

        enddo ! NGLLX
      enddo ! NGLLY
    enddo ! NGLLZ

    ! sum contributions from each element to the global mesh and add gravity terms
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          accel_crust_mantle(1,iglob) = accel_crust_mantle(1,iglob) + sum_terms(1,i,j,k)
          accel_crust_mantle(2,iglob) = accel_crust_mantle(2,iglob) + sum_terms(2,i,j,k)
          accel_crust_mantle(3,iglob) = accel_crust_mantle(3,iglob) + sum_terms(3,i,j,k)
        enddo
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
      !epsilondev(:,:,:,:,ispec) = epsilondev_loc(:,:,:,:)
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            epsilondev_xx(i,j,k,ispec) = epsilondev_loc(i,j,k,1)
            epsilondev_yy(i,j,k,ispec) = epsilondev_loc(i,j,k,2)
            epsilondev_xy(i,j,k,ispec) = epsilondev_loc(i,j,k,3)
            epsilondev_xz(i,j,k,ispec) = epsilondev_loc(i,j,k,4)
            epsilondev_yz(i,j,k,ispec) = epsilondev_loc(i,j,k,5)
          enddo
        enddo
      enddo
    endif

  enddo   ! spectral element loop NSPEC_CRUST_MANTLE

  end subroutine compute_forces_crust_mantle_noDev

