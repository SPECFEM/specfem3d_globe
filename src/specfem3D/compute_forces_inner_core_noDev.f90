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

  subroutine compute_forces_inner_core_noDev( NSPEC_STR_OR_ATT,NGLOB,NSPEC_ATT, &
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
                                              factor_common,vnspec, &
                                              pgrav_inner_core)

  use constants_solver

  use specfem_par, only: &
    hprime_xx,hprime_yy,hprime_zz,hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
    wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
    gravity_pre_store => gravity_pre_store_inner_core,gravity_H => gravity_H_inner_core, &
    COMPUTE_AND_STORE_STRAIN,USE_LDDRK

  use specfem_par_innercore, only: &
    xix => xix_inner_core,xiy => xiy_inner_core,xiz => xiz_inner_core, &
    etax => etax_inner_core,etay => etay_inner_core,etaz => etaz_inner_core, &
    gammax => gammax_inner_core,gammay => gammay_inner_core,gammaz => gammaz_inner_core, &
    kappavstore => kappavstore_inner_core, &
    muvstore => muvstore_inner_core, &
    c11store => c11store_inner_core,c12store => c12store_inner_core,c13store => c13store_inner_core, &
    c33store => c33store_inner_core,c44store => c44store_inner_core, &
    ibool => ibool_inner_core,idoubling => idoubling_inner_core, &
    phase_ispec_inner => phase_ispec_inner_inner_core, &
    nspec_outer => nspec_outer_inner_core, &
    nspec_inner => nspec_inner_inner_core
    !not needed: one_minus_sum_beta => one_minus_sum_beta_inner_core

  ! full gravity
  use specfem_par_full_gravity, only: &
    gravity_rho => gravity_rho_inner_core

  implicit none

  integer :: NSPEC_STR_OR_ATT,NGLOB,NSPEC_ATT

  ! time step
  real(kind=CUSTOM_REAL) deltat

  ! displacement, velocity and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: displ_inner_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: accel_inner_core

  ! for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  ! variable lengths for factor_common and one_minus_sum_beta

  ! variable sized array variables
  integer :: vnspec
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,N_SLS,vnspec) :: factor_common
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATT) :: R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATT) :: &
    R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STR_OR_ATT) :: &
    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STRAIN_ONLY) :: epsilon_trace_over_3

  ! inner/outer element run flag
  integer,intent(in) :: iphase

  ! full gravity
  real(kind=CUSTOM_REAL), dimension(NGLOB),intent(in) :: pgrav_inner_core

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5) :: epsilondev_loc

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  integer :: ispec,iglob
  integer :: i,j,k,l

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz

  real(kind=CUSTOM_REAL) hp1,hp2,hp3
  real(kind=CUSTOM_REAL) fac1,fac2,fac3
  real(kind=CUSTOM_REAL) lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) kappal
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NDIM) :: sum_terms

!  real(kind=CUSTOM_REAL) minus_sum_beta
  real(kind=CUSTOM_REAL) c11l,c33l,c12l,c13l,c44l

  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) tempz1l,tempz2l,tempz3l
  real(kind=CUSTOM_REAL) templ

  real(kind=CUSTOM_REAL) R_xx_val,R_yy_val
  integer :: i_SLS

  ! for gravity
  real(kind=CUSTOM_REAL) factor,gxl,gyl,gzl,sx_l,sy_l,sz_l
  real(kind=CUSTOM_REAL) Hxxl,Hyyl,Hzzl,Hxyl,Hxzl,Hyzl
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: rho_s_H
  real(kind=CUSTOM_REAL) sigma_yx,sigma_zx,sigma_zy

  ! full gravity
  real(kind=CUSTOM_REAL), dimension(9,NGLLX,NGLLY,NGLLZ) :: deriv_loc

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

    ! exclude fictitious elements in central cube
    if (idoubling(ispec) /= IFLAG_IN_FICTITIOUS_CUBE) then

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
              tempx1l = tempx1l + displ_inner_core(1,iglob)*hp1
              tempy1l = tempy1l + displ_inner_core(2,iglob)*hp1
              tempz1l = tempz1l + displ_inner_core(3,iglob)*hp1
              !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

              !!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l = 1,NGLLY
              hp2 = hprime_yy(j,l)
              iglob = ibool(i,l,k,ispec)
              tempx2l = tempx2l + displ_inner_core(1,iglob)*hp2
              tempy2l = tempy2l + displ_inner_core(2,iglob)*hp2
              tempz2l = tempz2l + displ_inner_core(3,iglob)*hp2
              !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

              !!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l = 1,NGLLZ
              hp3 = hprime_zz(k,l)
              iglob = ibool(i,j,l,ispec)
              tempx3l = tempx3l + displ_inner_core(1,iglob)*hp3
              tempy3l = tempy3l + displ_inner_core(2,iglob)*hp3
              tempz3l = tempz3l + displ_inner_core(3,iglob)*hp3
            enddo

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

            ! full gravity
            if (FULL_GRAVITY_VAL .and. .not. DISCARD_GCONTRIB) then
              ! note: deriv_mapping_inner_core is only available when chosen to use Deville routines
              !       thus, we need to fill a local deriv array here for the function call to pass on later
              deriv_loc(1,i,j,k) = xixl
              deriv_loc(2,i,j,k) = xiyl
              deriv_loc(3,i,j,k) = xizl
              deriv_loc(4,i,j,k) = etaxl
              deriv_loc(5,i,j,k) = etayl
              deriv_loc(6,i,j,k) = etazl
              deriv_loc(7,i,j,k) = gammaxl
              deriv_loc(8,i,j,k) = gammayl
              deriv_loc(9,i,j,k) = gammazl
            endif

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
              if (NSPEC_INNER_CORE_STRAIN_ONLY > 1) epsilon_trace_over_3(i,j,k,ispec) = templ

              epsilondev_loc(i,j,k,1) = duxdxl - templ
              epsilondev_loc(i,j,k,2) = duydyl - templ
              epsilondev_loc(i,j,k,3) = 0.5 * duxdyl_plus_duydxl
              epsilondev_loc(i,j,k,4) = 0.5 * duzdxl_plus_duxdzl
              epsilondev_loc(i,j,k,5) = 0.5 * duzdyl_plus_duydzl
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
              c11l = c11store(i,j,k,ispec)
              c12l = c12store(i,j,k,ispec)
              c13l = c13store(i,j,k,ispec)
              c33l = c33store(i,j,k,ispec)
              c44l = c44store(i,j,k,ispec)

              ! use unrelaxed parameters if attenuation
              ! already done in prepare_timerun_elastic_elements() routine...
              !if (ATTENUATION_VAL) then
              !  if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
              !    minus_sum_beta =  one_minus_sum_beta(i,j,k,ispec) - 1.0_CUSTOM_REAL
              !  else
              !    minus_sum_beta =  one_minus_sum_beta(1,1,1,ispec) - 1.0_CUSTOM_REAL
              !  endif
              !  mul = muvstore(i,j,k,ispec)
              !  c11l = c11l + FOUR_THIRDS * minus_sum_beta * mul
              !  c12l = c12l - TWO_THIRDS * minus_sum_beta * mul
              !  c13l = c13l - TWO_THIRDS * minus_sum_beta * mul
              !  c33l = c33l + FOUR_THIRDS * minus_sum_beta * mul
              !  c44l = c44l + minus_sum_beta * mul
              !endif

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
              ! already done in prepare_timerun_elastic_elements() routine...
              !if (ATTENUATION_VAL) then
              !  if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
              !    mul = mul * one_minus_sum_beta(i,j,k,ispec)
              !  else
              !    mul = mul * one_minus_sum_beta(1,1,1,ispec)
              !  endif
              !endif

              lambdalplus2mul = kappal + FOUR_THIRDS * mul
              lambdal = lambdalplus2mul - 2._CUSTOM_REAL * mul

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
              sx_l = displ_inner_core(1,iglob)
              sy_l = displ_inner_core(2,iglob)
              sz_l = displ_inner_core(3,iglob)

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

          enddo
        enddo
      enddo

      ! full gravity
      if (FULL_GRAVITY_VAL .and. .not. DISCARD_GCONTRIB) then
        call SIEM_solve_element_add_full_gravity(ispec,NSPEC_INNER_CORE,NGLOB,gravity_rho,deriv_loc,ibool, &
                                                 pgrav_inner_core,rho_s_H)
      endif

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

            sum_terms(i,j,k,1) = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l)
            sum_terms(i,j,k,2) = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l)
            sum_terms(i,j,k,3) = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l)

            if (GRAVITY_VAL) then
              sum_terms(i,j,k,1) = sum_terms(i,j,k,1) + rho_s_H(1,i,j,k)
              sum_terms(i,j,k,2) = sum_terms(i,j,k,2) + rho_s_H(2,i,j,k)
              sum_terms(i,j,k,3) = sum_terms(i,j,k,3) + rho_s_H(3,i,j,k)
            endif
          enddo
        enddo
      enddo

      ! sum contributions from each element to the global mesh and add gravity terms
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)
            accel_inner_core(1,iglob) = accel_inner_core(1,iglob) + sum_terms(i,j,k,1)
            accel_inner_core(2,iglob) = accel_inner_core(2,iglob) + sum_terms(i,j,k,2)
            accel_inner_core(3,iglob) = accel_inner_core(3,iglob) + sum_terms(i,j,k,3)
          enddo
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

    endif   ! end test to exclude fictitious elements in central cube

  enddo ! spectral element loop

  end subroutine compute_forces_inner_core_noDev

