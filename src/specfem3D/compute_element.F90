!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            August 2013
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
                    R_xx,R_yy,R_xy,R_xz,R_yz, &
                    epsilon_trace_over_3, &
                    one_minus_sum_beta,vx,vy,vz,vnspec, &
                    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                    dummyx_loc,dummyy_loc,dummyz_loc, &
                    epsilondev_loc,rho_s_H,is_backward_field)

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
  integer :: vx,vy,vz,vnspec

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUATION) :: R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: epsilon_trace_over_3

  real(kind=CUSTOM_REAL), dimension(vx,vy,vz,vnspec) :: one_minus_sum_beta

  ! gravity
  double precision, dimension(NRAD_GRAVITY) :: minus_gravity_table,density_table,minus_deriv_gravity_table

  ! element info
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: rho_s_H
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc

  logical :: is_backward_field

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

  integer :: ispec_strain
  integer :: i,j,k
  integer :: int_radius
  integer :: iglob

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
          templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
          if(NSPEC_CRUST_MANTLE_STRAIN_ONLY == 1) then
             ispec_strain = 1
             epsilon_trace_over_3(i,j,k,ispec_strain) = templ
          else
             ispec_strain = ispec
             epsilon_trace_over_3(i,j,k,ispec_strain) = templ
          endif
          epsilondev_loc(1,i,j,k) = duxdxl - templ
          epsilondev_loc(2,i,j,k) = duydyl - templ
          epsilondev_loc(3,i,j,k) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
          epsilondev_loc(4,i,j,k) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
          epsilondev_loc(5,i,j,k) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl
        endif

        !
        ! compute  isotropic  elements
        !

        ! layer with no transverse isotropy, use kappav and muv
        kappal = kappavstore(i,j,k,ispec)
        mul = muvstore(i,j,k,ispec)

        ! use unrelaxed parameters if attenuation
        if(ATTENUATION_VAL .and. (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL)) then
          ! precompute terms for attenuation if needed
          one_minus_sum_beta_use = one_minus_sum_beta(i,j,k,ispec)
          mul = mul * one_minus_sum_beta_use
        else if( ATTENUATION_VAL ) then
          one_minus_sum_beta_use = one_minus_sum_beta(1,1,1,ispec)
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
        if(ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL ) then

!daniel: att - debug update
!          call compute_element_att_mem_up_cm(ispec,i,j,k, &
!                                          R_xx(1,i,j,k,ispec), &
!                                          R_yy(1,i,j,k,ispec), &
!                                          R_xy(1,i,j,k,ispec), &
!                                          R_xz(1,i,j,k,ispec), &
!                                          R_yz(1,i,j,k,ispec), &
!                                          epsilondev_loc(:,i,j,k),muvstore(i,j,k,ispec),is_backward_field)
! dummy to avoid compiler warning
          if( is_backward_field ) then
          endif

          ! note: fortran passes pointers to array location, thus R_memory(1,1,...) should be fine
          call compute_element_att_stress(R_xx(1,i,j,k,ispec), &
                                          R_yy(1,i,j,k,ispec), &
                                          R_xy(1,i,j,k,ispec), &
                                          R_xz(1,i,j,k,ispec), &
                                          R_yz(1,i,j,k,ispec), &
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
            iglob = ibool(i,j,k,ispec)

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

            ! distinguish between single and double precision for reals
            if(CUSTOM_REAL == SIZE_REAL) then

              ! get displacement and multiply by density to compute G tensor
              sx_l = rho * dble(dummyx_loc(i,j,k)) ! dble(displ_crust_mantle(1,iglob))
              sy_l = rho * dble(dummyy_loc(i,j,k)) ! dble(displ_crust_mantle(2,iglob))
              sz_l = rho * dble(dummyz_loc(i,j,k)) ! dble(displ_crust_mantle(3,iglob))

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
              sx_l = rho * dummyx_loc(i,j,k) ! displ_crust_mantle(1,iglob)
              sy_l = rho * dummyy_loc(i,j,k) ! displ_crust_mantle(2,iglob)
              sz_l = rho * dummyz_loc(i,j,k) ! displ_crust_mantle(3,iglob)

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
                    R_xx,R_yy,R_xy,R_xz,R_yz, &
                    epsilon_trace_over_3, &
                    one_minus_sum_beta,vx,vy,vz,vnspec, &
                    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                    dummyx_loc,dummyy_loc,dummyz_loc, &
                    epsilondev_loc,rho_s_H,is_backward_field)

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
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUATION) :: R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: epsilon_trace_over_3

  ! variable sized array variables
  integer :: vx,vy,vz,vnspec

  ! [alpha,beta,gamma]val reduced to N_SLS  to N_SLS*NUM_NODES
  real(kind=CUSTOM_REAL), dimension(vx,vy,vz,vnspec) :: one_minus_sum_beta

  ! gravity
  double precision, dimension(NRAD_GRAVITY) :: minus_gravity_table,density_table,minus_deriv_gravity_table

  ! element info
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: rho_s_H
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc

  logical :: is_backward_field

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

  integer :: ispec_strain
  integer :: i,j,k
  integer :: int_radius
  integer :: iglob

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
          templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
          if(NSPEC_CRUST_MANTLE_STRAIN_ONLY == 1) then
             ispec_strain = 1
             epsilon_trace_over_3(i,j,k,ispec_strain) = templ
          else
             ispec_strain = ispec
             epsilon_trace_over_3(i,j,k,ispec_strain) = templ
          endif
          epsilondev_loc(1,i,j,k) = duxdxl - templ
          epsilondev_loc(2,i,j,k) = duydyl - templ
          epsilondev_loc(3,i,j,k) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
          epsilondev_loc(4,i,j,k) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
          epsilondev_loc(5,i,j,k) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl
        endif

        !
        ! compute either isotropic or anisotropic elements
        !

! note: the mesh is built such that anisotropic elements are created first in anisotropic layers,
!           thus they are listed first ( see in create_regions_mesh.f90: perm_layer() ordering )
!           this is therefore still in bounds of 1:NSPECMAX_TISO_MANTLE even if NSPECMAX_TISO is less than NSPEC

        ! use kappa and mu from transversely isotropic model
        kappavl = kappavstore(i,j,k,ispec)
        muvl = muvstore(i,j,k,ispec)

        kappahl = kappahstore(i,j,k,ispec)
        muhl = muhstore(i,j,k,ispec)

        ! use unrelaxed parameters if attenuation
        ! eta does not need to be shifted since it is a ratio
        if(ATTENUATION_VAL .and. (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL)) then
          ! precompute terms for attenuation if needed
          one_minus_sum_beta_use = one_minus_sum_beta(i,j,k,ispec)
          muvl = muvl * one_minus_sum_beta_use
          muhl = muhl * one_minus_sum_beta_use
        else if( ATTENUATION_VAL ) then
          one_minus_sum_beta_use = one_minus_sum_beta(1,1,1,ispec)
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
        iglob = ibool(i,j,k,ispec)

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
        if(ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL ) then

          ! note: Fortran passes pointers to array location, thus R_memory(1,1,...) is fine
          call compute_element_att_stress(R_xx(1,i,j,k,ispec), &
                                          R_yy(1,i,j,k,ispec), &
                                          R_xy(1,i,j,k,ispec), &
                                          R_xz(1,i,j,k,ispec), &
                                          R_yz(1,i,j,k,ispec), &
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
            iglob = ibool(i,j,k,ispec)

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
                    R_xx,R_yy,R_xy,R_xz,R_yz, &
                    epsilon_trace_over_3, &
                    one_minus_sum_beta,vx,vy,vz,vnspec, &
                    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                    dummyx_loc,dummyy_loc,dummyz_loc, &
                    epsilondev_loc,rho_s_H,is_backward_field)

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
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUATION) :: R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: epsilon_trace_over_3

  ! variable sized array variables
  integer :: vx,vy,vz,vnspec

  ! [alpha,beta,gamma]val reduced to N_SLS  to N_SLS*NUM_NODES
  real(kind=CUSTOM_REAL), dimension(vx,vy,vz,vnspec) :: one_minus_sum_beta

  ! gravity
  double precision, dimension(NRAD_GRAVITY) :: minus_gravity_table,density_table,minus_deriv_gravity_table

  ! element info
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: rho_s_H
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc

  logical :: is_backward_field

! local parameters
  ! real(kind=CUSTOM_REAL) one_minus_sum_beta_use
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

  integer :: ispec_strain
  integer :: i,j,k
  integer :: int_radius
  integer :: iglob

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
          templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
          if(NSPEC_CRUST_MANTLE_STRAIN_ONLY == 1) then
             ispec_strain = 1
             epsilon_trace_over_3(i,j,k,ispec_strain) = templ
          else
             ispec_strain = ispec
             epsilon_trace_over_3(i,j,k,ispec_strain) = templ
          endif
          epsilondev_loc(1,i,j,k) = duxdxl - templ
          epsilondev_loc(2,i,j,k) = duydyl - templ
          epsilondev_loc(3,i,j,k) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
          epsilondev_loc(4,i,j,k) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
          epsilondev_loc(5,i,j,k) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl
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

        if(ATTENUATION_VAL .and. (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL)) then
          ! precompute terms for attenuation if needed
          minus_sum_beta =  one_minus_sum_beta(i,j,k,ispec) - 1.0_CUSTOM_REAL
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
        else if( ATTENUATION_VAL ) then
          minus_sum_beta =  one_minus_sum_beta(1,1,1,ispec) - 1.0_CUSTOM_REAL
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
        if(ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL ) then

          ! note: Fortran passes pointers to array location, thus R_memory(1,1,...) is fine
          call compute_element_att_stress(R_xx(1,i,j,k,ispec), &
                                          R_yy(1,i,j,k,ispec), &
                                          R_xy(1,i,j,k,ispec), &
                                          R_xz(1,i,j,k,ispec), &
                                          R_yz(1,i,j,k,ispec), &
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
            iglob = ibool(i,j,k,ispec)

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

            ! distinguish between single and double precision for reals
            if(CUSTOM_REAL == SIZE_REAL) then

              ! get displacement and multiply by density to compute G tensor
              sx_l = rho * dble(dummyx_loc(i,j,k)) ! dble(displ_crust_mantle(1,iglob))
              sy_l = rho * dble(dummyy_loc(i,j,k)) ! dble(displ_crust_mantle(2,iglob))
              sz_l = rho * dble(dummyz_loc(i,j,k)) !  dble(displ_crust_mantle(3,iglob))

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
              sx_l = rho * dummyx_loc(i,j,k)  ! displ_crust_mantle(1,iglob)
              sy_l = rho * dummyy_loc(i,j,k)  !  displ_crust_mantle(2,iglob)
              sz_l = rho * dummyz_loc(i,j,k)  ! displ_crust_mantle(3,iglob)

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


  subroutine compute_element_att_stress(R_xx_loc,R_yy_loc,R_xy_loc,R_xz_loc,R_yz_loc, &
                                       sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz)

  use constants_solver

  implicit none

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
!  real(kind=CUSTOM_REAL), dimension(5,N_SLS) :: R_memory_loc
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: R_xx_loc
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: R_yy_loc
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: R_xy_loc
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: R_xz_loc
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: R_yz_loc

  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz

! local parameters
  real(kind=CUSTOM_REAL) R_xx_val1,R_yy_val1
  integer :: i_SLS

  do i_SLS = 1,N_SLS
    R_xx_val1 = R_xx_loc(i_SLS) ! R_memory(1,i_SLS,i,j,k,ispec)
    R_yy_val1 = R_yy_loc(i_SLS) ! R_memory(2,i_SLS,i,j,k,ispec)
    sigma_xx = sigma_xx - R_xx_val1
    sigma_yy = sigma_yy - R_yy_val1
    sigma_zz = sigma_zz + R_xx_val1 + R_yy_val1
    sigma_xy = sigma_xy - R_xy_loc(i_SLS) ! R_memory(3,i_SLS,i,j,k,ispec)
    sigma_xz = sigma_xz - R_xz_loc(i_SLS) ! R_memory(4,i_SLS,i,j,k,ispec)
    sigma_yz = sigma_yz - R_yz_loc(i_SLS) ! R_memory(5,i_SLS,i,j,k,ispec)
  enddo

  end subroutine compute_element_att_stress

!
!--------------------------------------------------------------------------------------------
!

  subroutine compute_element_att_memory_cm(ispec,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                        vx,vy,vz,vnspec,factor_common, &
                                        alphaval,betaval,gammaval, &
                                        c44store,muvstore, &
                                        epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
                                        epsilondev_loc,is_backward_field)
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

  use constants_solver

  implicit none

  ! element id
  integer :: ispec

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUATION) :: R_xx,R_yy,R_xy,R_xz,R_yz

  ! variable sized array variables
  integer :: vx,vy,vz,vnspec

  real(kind=CUSTOM_REAL), dimension(N_SLS,vx,vy,vz,vnspec) :: factor_common
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE) :: c44store
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ISO_MANTLE) :: muvstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc

  logical :: is_backward_field

! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: factor_common_c44_muv
  integer :: i_SLS

  ! use Runge-Kutta scheme to march in time

  ! get coefficients for that standard linear solid
  ! IMPROVE we use mu_v here even if there is some anisotropy
  ! IMPROVE we should probably use an average value instead

  do i_SLS = 1,N_SLS

    ! reformatted R_memory to handle large factor_common and reduced [alpha,beta,gamma]val
    if( ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL ) then
      if(ANISOTROPIC_3D_MANTLE_VAL) then
        factor_common_c44_muv(:,:,:) = factor_common(i_SLS,:,:,:,ispec) * c44store(:,:,:,ispec)
      else
        factor_common_c44_muv(:,:,:) = factor_common(i_SLS,:,:,:,ispec) * muvstore(:,:,:,ispec)
      endif
    else
      if(ANISOTROPIC_3D_MANTLE_VAL) then
        factor_common_c44_muv(:,:,:) = factor_common(i_SLS,1,1,1,ispec) * c44store(:,:,:,ispec)
      else
        factor_common_c44_muv(:,:,:) = factor_common(i_SLS,1,1,1,ispec) * muvstore(:,:,:,ispec)
      endif
    endif

    R_xx(i_SLS,:,:,:,ispec) = alphaval(i_SLS) * R_xx(i_SLS,:,:,:,ispec) + factor_common_c44_muv(:,:,:) * &
          (betaval(i_SLS) * epsilondev_xx(:,:,:,ispec) + gammaval(i_SLS) * epsilondev_loc(1,:,:,:))

    R_yy(i_SLS,:,:,:,ispec) = alphaval(i_SLS) * R_yy(i_SLS,:,:,:,ispec) + factor_common_c44_muv(:,:,:) * &
          (betaval(i_SLS) * epsilondev_yy(:,:,:,ispec) + gammaval(i_SLS) * epsilondev_loc(2,:,:,:))

    R_xy(i_SLS,:,:,:,ispec) = alphaval(i_SLS) * R_xy(i_SLS,:,:,:,ispec) + factor_common_c44_muv(:,:,:) * &
          (betaval(i_SLS) * epsilondev_xy(:,:,:,ispec) + gammaval(i_SLS) * epsilondev_loc(3,:,:,:))

    R_xz(i_SLS,:,:,:,ispec) = alphaval(i_SLS) * R_xz(i_SLS,:,:,:,ispec) + factor_common_c44_muv(:,:,:) * &
          (betaval(i_SLS) * epsilondev_xz(:,:,:,ispec) + gammaval(i_SLS) * epsilondev_loc(4,:,:,:))

    R_yz(i_SLS,:,:,:,ispec) = alphaval(i_SLS) * R_yz(i_SLS,:,:,:,ispec) + factor_common_c44_muv(:,:,:) * &
          (betaval(i_SLS) * epsilondev_yz(:,:,:,ispec) + gammaval(i_SLS) * epsilondev_loc(5,:,:,:))

  enddo ! i_SLS

  end subroutine compute_element_att_memory_cm

!
!--------------------------------------------------------------------------------------------
!

  subroutine compute_element_att_memory_ic(ispec,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                        vx,vy,vz,vnspec,factor_common, &
                                        alphaval,betaval,gammaval, &
                                        muvstore, &
                                        epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz, &
                                        epsilondev_loc,is_backward_field)
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

  use constants_solver

  implicit none

  ! element id
  integer :: ispec

  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ATTENUATION) :: &
    R_xx,R_yy,R_xy,R_xz,R_yz

  ! variable sized array variables
  integer :: vx,vy,vz,vnspec

  real(kind=CUSTOM_REAL), dimension(N_SLS,vx,vy,vz,vnspec) :: factor_common
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: muvstore

!  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: epsilondev
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: &
    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc

  logical :: is_backward_field

! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: factor_common_use

  integer :: i_SLS

  ! use Runge-Kutta scheme to march in time

  ! get coefficients for that standard linear solid
  ! IMPROVE we use mu_v here even if there is some anisotropy
  ! IMPROVE we should probably use an average value instead

  do i_SLS = 1,N_SLS
    ! reformatted R_memory to handle large factor_common and reduced [alpha,beta,gamma]val
    if( ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL ) then
      factor_common_use(:,:,:) = factor_common(i_SLS,:,:,:,ispec) * muvstore(:,:,:,ispec)
    else
      factor_common_use(:,:,:) = factor_common(i_SLS,1,1,1,ispec) * muvstore(:,:,:,ispec)
    endif

!    do i_memory = 1,5
!       R_memory(i_memory,i_SLS,:,:,:,ispec) = alphaval(i_SLS) * R_memory(i_memory,i_SLS,:,:,:,ispec) &
!            + muvstore(:,:,:,ispec) * factor_common_use(:,:,:) * &
!            (betaval(i_SLS) * epsilondev(i_memory,:,:,:,ispec) + gammaval(i_SLS) * epsilondev_loc(i_memory,:,:,:))
!    enddo

    R_xx(i_SLS,:,:,:,ispec) = alphaval(i_SLS) * R_xx(i_SLS,:,:,:,ispec) + factor_common_use(:,:,:) * &
         (betaval(i_SLS) * epsilondev_xx(:,:,:,ispec) + gammaval(i_SLS) * epsilondev_loc(1,:,:,:))

    R_yy(i_SLS,:,:,:,ispec) = alphaval(i_SLS) * R_yy(i_SLS,:,:,:,ispec) + factor_common_use(:,:,:) * &
         (betaval(i_SLS) * epsilondev_yy(:,:,:,ispec) + gammaval(i_SLS) * epsilondev_loc(2,:,:,:))

    R_xy(i_SLS,:,:,:,ispec) = alphaval(i_SLS) * R_xy(i_SLS,:,:,:,ispec) + factor_common_use(:,:,:) * &
         (betaval(i_SLS) * epsilondev_xy(:,:,:,ispec) + gammaval(i_SLS) * epsilondev_loc(3,:,:,:))

    R_xz(i_SLS,:,:,:,ispec) = alphaval(i_SLS) * R_xz(i_SLS,:,:,:,ispec) + factor_common_use(:,:,:) * &
         (betaval(i_SLS) * epsilondev_xz(:,:,:,ispec) + gammaval(i_SLS) * epsilondev_loc(4,:,:,:))

    R_yz(i_SLS,:,:,:,ispec) = alphaval(i_SLS) * R_yz(i_SLS,:,:,:,ispec) + factor_common_use(:,:,:) * &
         (betaval(i_SLS) * epsilondev_yz(:,:,:,ispec) + gammaval(i_SLS) * epsilondev_loc(5,:,:,:))

  enddo

  end subroutine compute_element_att_memory_ic

!
!--------------------------------------------------------------------------------------------
!

  subroutine compute_element_att_mem_up_cm(ispec,i,j,k, &
                                              R_xx_loc,R_yy_loc,R_xy_loc,R_xz_loc,R_yz_loc, &
                                              epsilondev_loc,c44_muv,is_backward_field)
! crust mantle
! update memory variables based upon the Runge-Kutta scheme


!daniel: att - debug update
  use specfem_par,only: tau_sigma_dble,deltat,b_deltat

  use specfem_par_crustmantle,only: factor_common=>factor_common_crust_mantle

  use constants_solver

  implicit none

  ! element id
  integer :: ispec,i,j,k

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
!  real(kind=CUSTOM_REAL), dimension(5,N_SLS) :: R_memory_loc
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: R_xx_loc
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: R_yy_loc
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: R_xy_loc
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: R_xz_loc
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: R_yz_loc

  real(kind=CUSTOM_REAL), dimension(5) :: epsilondev_loc
  real(kind=CUSTOM_REAL) :: c44_muv

  logical :: is_backward_field
  double precision :: dt,kappa

! local parameters
  real(kind=CUSTOM_REAL) :: factor_common_c44_muv
  integer :: i_SLS

  if( .not. is_backward_field ) then
    dt = dble(deltat)
  else
    ! backward/reconstruction: reverse time
    dt = dble(b_deltat)
  endif

  do i_SLS = 1,N_SLS

    ! runge-kutta scheme to update memory variables R(t)
    if( .false. ) then
! classical RK 4:       R'(t) =  - 1/tau * R(t)
!
! Butcher RK4:
! 0     |
! 1/2   | 1/2
! 1/2   | 0    1/2
! 1     | 0          1
! -----------------------------------------------------------------------------
!         1/6  1/3   1/3   1/6
    kappa = - dt/tau_sigma_dble(i_SLS)

    R_xx_loc(i_SLS) = R_xx_loc(i_SLS) * &
      (1.0d0 + kappa*(1.d0 + 0.5d0*kappa*(1.d0 + 1.0d0/6.0d0*kappa*(1.d0 + 1.0d0/24.0d0*kappa))))
    R_yy_loc(i_SLS) = R_yy_loc(i_SLS) * &
      (1.0d0 + kappa*(1.d0 + 0.5d0*kappa*(1.d0 + 1.0d0/6.0d0*kappa*(1.d0 + 1.0d0/24.0d0*kappa))))
    R_xy_loc(i_SLS) = R_xy_loc(i_SLS) * &
      (1.0d0 + kappa*(1.d0 + 0.5d0*kappa*(1.d0 + 1.0d0/6.0d0*kappa*(1.d0 + 1.0d0/24.0d0*kappa))))
    R_xz_loc(i_SLS) = R_xz_loc(i_SLS) * &
      (1.0d0 + kappa*(1.d0 + 0.5d0*kappa*(1.d0 + 1.0d0/6.0d0*kappa*(1.d0 + 1.0d0/24.0d0*kappa))))
    R_yz_loc(i_SLS) = R_yz_loc(i_SLS) * &
      (1.0d0 + kappa*(1.d0 + 0.5d0*kappa*(1.d0 + 1.0d0/6.0d0*kappa*(1.d0 + 1.0d0/24.0d0*kappa))))
    endif

    ! reformatted R_memory to handle large factor_common and reduced [alpha,beta,gamma]val
    if( ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL ) then
      factor_common_c44_muv = factor_common(i_SLS,i,j,k,ispec) * c44_muv
    else
      factor_common_c44_muv = factor_common(i_SLS,1,1,1,ispec) * c44_muv
    endif

    ! adds contributions from current strain
    R_xx_loc(i_SLS) = R_xx_loc(i_SLS) + 0.5d0 * dt * dble(factor_common_c44_muv) * dble(epsilondev_loc(1))
    R_yy_loc(i_SLS) = R_yy_loc(i_SLS) + 0.5d0 * dt * dble(factor_common_c44_muv) * dble(epsilondev_loc(2))
    R_xy_loc(i_SLS) = R_xy_loc(i_SLS) + 0.5d0 * dt * dble(factor_common_c44_muv) * dble(epsilondev_loc(3))
    R_xz_loc(i_SLS) = R_xz_loc(i_SLS) + 0.5d0 * dt * dble(factor_common_c44_muv) * dble(epsilondev_loc(4))
    R_yz_loc(i_SLS) = R_yz_loc(i_SLS) + 0.5d0 * dt * dble(factor_common_c44_muv) * dble(epsilondev_loc(5))

  enddo ! i_SLS

  end subroutine compute_element_att_mem_up_cm
