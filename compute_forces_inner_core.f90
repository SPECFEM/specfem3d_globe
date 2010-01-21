!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
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

  subroutine compute_forces_inner_core(minus_gravity_table,density_table,minus_deriv_gravity_table, &
          displ,accel,xstore,ystore,zstore, &
          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
          hprime_xx,hprime_yy,hprime_zz, &
          hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
          kappavstore,muvstore,ibool,idoubling, &
          c11store,c33store,c12store,c13store,c44store,R_memory,epsilondev,epsilon_trace_over_3,&
          one_minus_sum_beta,alphaval,betaval,gammaval,factor_common, &
          vx,vy,vz,vnspec,COMPUTE_AND_STORE_STRAIN)

  implicit none

  include "constants.h"

! include values created by the mesher
! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

! attenuation_model_variables
  type attenuation_model_variables
    sequence
    double precision min_period, max_period
    double precision                          :: QT_c_source        ! Source Frequency
    double precision, dimension(:), pointer   :: Qtau_s             ! tau_sigma
    double precision, dimension(:), pointer   :: QrDisc             ! Discontinutitues Defined
    double precision, dimension(:), pointer   :: Qr                 ! Radius
    integer, dimension(:), pointer            :: interval_Q                 ! Steps
    double precision, dimension(:), pointer   :: Qmu                ! Shear Attenuation
    double precision, dimension(:,:), pointer :: Qtau_e             ! tau_epsilon
    double precision, dimension(:), pointer   :: Qomsb, Qomsb2      ! one_minus_sum_beta
    double precision, dimension(:,:), pointer :: Qfc, Qfc2          ! factor_common
    double precision, dimension(:), pointer   :: Qsf, Qsf2          ! scale_factor
    integer, dimension(:), pointer            :: Qrmin              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: Qrmax              ! Max and Mins of idoubling
    integer                                   :: Qn                 ! Number of points
  end type attenuation_model_variables

!> Hejun
!  type (attenuation_model_variables) AM_V
! attenuation_model_variables
!< Hejun

! for forward or backward simulations
  logical COMPUTE_AND_STORE_STRAIN

! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: displ,accel

! for attenuation
! memory variables R_ij are stored at the local rather than global level
! to allow for optimization of cache access by compiler
  integer i_sls,i_memory
  real(kind=CUSTOM_REAL) R_xx_val,R_yy_val

! variable lengths for factor_common and one_minus_sum_beta
  integer vx, vy, vz, vnspec

  real(kind=CUSTOM_REAL), dimension(vx, vy, vz, vnspec) :: one_minus_sum_beta

  real(kind=CUSTOM_REAL), dimension(N_SLS, vx, vy, vz, vnspec) :: factor_common
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ) :: factor_common_use

  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ATTENUATION) :: R_memory
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: epsilondev
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: epsilon_trace_over_3

! array with the local to global mapping per slice
  integer, dimension(NSPEC_INNER_CORE) :: idoubling

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: xix,xiy,xiz, &
                      etax,etay,etaz,gammax,gammay,gammaz

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: kappavstore,muvstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: &
    c11store,c33store,c12store,c13store,c44store

  integer ispec,iglob,ispec_strain
  integer i,j,k,l

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

  real(kind=CUSTOM_REAL) minus_sum_beta
  real(kind=CUSTOM_REAL) c11l,c33l,c12l,c13l,c44l

  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) tempz1l,tempz2l,tempz3l

! for gravity
  integer int_radius
  real(kind=CUSTOM_REAL) sigma_yx,sigma_zx,sigma_zy
  double precision radius,rho,minus_g,minus_dg
  double precision minus_g_over_radius,minus_dg_plus_g_over_radius
  double precision cos_theta,sin_theta,cos_phi,sin_phi
  double precision cos_theta_sq,sin_theta_sq,cos_phi_sq,sin_phi_sq
  double precision theta,phi,factor,gxl,gyl,gzl,sx_l,sy_l,sz_l
  double precision Hxxl,Hyyl,Hzzl,Hxyl,Hxzl,Hyzl
  double precision, dimension(NRAD_GRAVITY) :: minus_gravity_table,density_table,minus_deriv_gravity_table
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: rho_s_H
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube
  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE) :: xstore,ystore,zstore
!> Hejun
!  integer iregion_selected
!  real(kind=CUSTOM_REAL) radius_cr
!< Hejun

! ****************************************************
!   big loop over all spectral elements in the solid
! ****************************************************

! set acceleration to zero
  accel(:,:) = 0._CUSTOM_REAL

  do ispec = 1,NSPEC_INNER_CORE

! exclude fictitious elements in central cube
    if(idoubling(ispec) /= IFLAG_IN_FICTITIOUS_CUBE) then

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
  if (COMPUTE_AND_STORE_STRAIN) then
    if(NSPEC_INNER_CORE_STRAIN_ONLY == 1) then
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

!> Hejun
  if(ATTENUATION_VAL) then
!     if(ATTENUATION_3D_VAL) then
        minus_sum_beta =  one_minus_sum_beta(i,j,k,ispec) - 1.0
!     else
!        radius_cr = xstore(ibool(i,j,k,ispec))
!        call get_attenuation_index(idoubling(ispec), dble(radius_cr), iregion_selected, .TRUE., AM_V)
!        minus_sum_beta =  one_minus_sum_beta(1,1,1,iregion_selected) - 1.0
!     endif ! ATTENUATION_3D_VAL
  endif ! ATTENUATION_VAL
!< Hejun

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

!> Hejun
! use unrelaxed parameters if attenuation
  if(ATTENUATION_VAL) then
!    if(ATTENUATION_3D_VAL) then
      mul = mul * one_minus_sum_beta(i,j,k,ispec)
!    else
!      mul = mul * one_minus_sum_beta(1,1,1,iregion_selected)
!    endif
  endif
!< Hejun
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
  if(ATTENUATION_VAL) then
    do i_sls = 1,N_SLS
      R_xx_val = R_memory(1,i_sls,i,j,k,ispec)
      R_yy_val = R_memory(2,i_sls,i,j,k,ispec)
      sigma_xx = sigma_xx - R_xx_val
      sigma_yy = sigma_yy - R_yy_val
      sigma_zz = sigma_zz + R_xx_val + R_yy_val
      sigma_xy = sigma_xy - R_memory(3,i_sls,i,j,k,ispec)
      sigma_xz = sigma_xz - R_memory(4,i_sls,i,j,k,ispec)
      sigma_yz = sigma_yz - R_memory(5,i_sls,i,j,k,ispec)
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

  iglob = ibool(i,j,k,ispec)
  radius = dble(xstore(iglob))
  theta = dble(ystore(iglob))
  phi = dble(zstore(iglob))

! make sure radius is never zero even for points at center of cube
! because we later divide by radius
  if(radius < 100.d0 / R_EARTH) radius = 100.d0 / R_EARTH

  cos_theta = dcos(theta)
  sin_theta = dsin(theta)
  cos_phi = dcos(phi)
  sin_phi = dsin(phi)

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

  cos_theta_sq = cos_theta**2
  sin_theta_sq = sin_theta**2
  cos_phi_sq = cos_phi**2
  sin_phi_sq = sin_phi**2

  Hxxl = minus_g_over_radius*(cos_phi_sq*cos_theta_sq + sin_phi_sq) + cos_phi_sq*minus_dg*sin_theta_sq
  Hyyl = minus_g_over_radius*(cos_phi_sq + cos_theta_sq*sin_phi_sq) + minus_dg*sin_phi_sq*sin_theta_sq
  Hzzl = cos_theta_sq*minus_dg + minus_g_over_radius*sin_theta_sq
  Hxyl = cos_phi*minus_dg_plus_g_over_radius*sin_phi*sin_theta_sq
  Hxzl = cos_phi*cos_theta*minus_dg_plus_g_over_radius*sin_theta
  Hyzl = cos_theta*minus_dg_plus_g_over_radius*sin_phi*sin_theta

  iglob = ibool(i,j,k,ispec)

! distinguish between single and double precision for reals
  if(CUSTOM_REAL == SIZE_REAL) then

! get displacement and multiply by density to compute G tensor
    sx_l = rho * dble(displ(1,iglob))
    sy_l = rho * dble(displ(2,iglob))
    sz_l = rho * dble(displ(3,iglob))

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
    sx_l = rho * displ(1,iglob)
    sy_l = rho * displ(2,iglob)
    sz_l = rho * displ(3,iglob)

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

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempy1l = 0._CUSTOM_REAL
          tempz1l = 0._CUSTOM_REAL

          tempx2l = 0._CUSTOM_REAL
          tempy2l = 0._CUSTOM_REAL
          tempz2l = 0._CUSTOM_REAL

          tempx3l = 0._CUSTOM_REAL
          tempy3l = 0._CUSTOM_REAL
          tempz3l = 0._CUSTOM_REAL

          do l=1,NGLLX
            fac1 = hprimewgll_xx(l,i)
            tempx1l = tempx1l + tempx1(l,j,k)*fac1
            tempy1l = tempy1l + tempy1(l,j,k)*fac1
            tempz1l = tempz1l + tempz1(l,j,k)*fac1
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
            fac2 = hprimewgll_yy(l,j)
            tempx2l = tempx2l + tempx2(i,l,k)*fac2
            tempy2l = tempy2l + tempy2(i,l,k)*fac2
            tempz2l = tempz2l + tempz2(i,l,k)*fac2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
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

     if(GRAVITY_VAL) sum_terms(:,i,j,k) = sum_terms(:,i,j,k) + rho_s_H(:,i,j,k)

        enddo
      enddo
    enddo

! sum contributions from each element to the global mesh and add gravity terms
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          accel(:,iglob) = accel(:,iglob) + sum_terms(:,i,j,k)
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

    if(ATTENUATION_VAL) then

       do i_sls = 1,N_SLS
!> Hejun
!          if(ATTENUATION_3D_VAL) then
             factor_common_use = factor_common(i_sls,:,:,:,ispec)
!          else
!             factor_common_use(:,:,:) = factor_common(i_sls,1,1,1,iregion_selected)
!          endif
!< Hejun
          do i_memory = 1,5
             R_memory(i_memory,i_sls,:,:,:,ispec) = &
                  alphaval(i_sls) * &
                  R_memory(i_memory,i_sls,:,:,:,ispec) + muvstore(:,:,:,ispec) * &
                  factor_common_use * &
                  (betaval(i_sls) * &
                  epsilondev(i_memory,:,:,:,ispec) + gammaval(i_sls) * epsilondev_loc(i_memory,:,:,:))
          enddo
       enddo

     endif

     if (COMPUTE_AND_STORE_STRAIN) then
! save deviatoric strain for Runge-Kutta scheme
       epsilondev(:,:,:,:,ispec) = epsilondev_loc(:,:,:,:)
     endif

   endif   ! end test to exclude fictitious elements in central cube

  enddo ! spectral element loop

  end subroutine compute_forces_inner_core


!
!--------------------------------------------------------------------------------------------------
!

  subroutine compute_forces_inner_core_Dev(minus_gravity_table,density_table,minus_deriv_gravity_table, &
          displ,accel,xstore,ystore,zstore, &
          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
          hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
          kappavstore,muvstore,ibool,idoubling, &
          c11store,c33store,c12store,c13store,c44store,R_memory,epsilondev,epsilon_trace_over_3,&
          one_minus_sum_beta,alphaval,betaval,gammaval,factor_common, &
          vx,vy,vz,vnspec,COMPUTE_AND_STORE_STRAIN)

  implicit none

  include "constants.h"

! include values created by the mesher
! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

! attenuation_model_variables
  type attenuation_model_variables
    sequence
    double precision min_period, max_period
    double precision                          :: QT_c_source        ! Source Frequency
    double precision, dimension(:), pointer   :: Qtau_s             ! tau_sigma
    double precision, dimension(:), pointer   :: QrDisc             ! Discontinutitues Defined
    double precision, dimension(:), pointer   :: Qr                 ! Radius
    integer, dimension(:), pointer            :: interval_Q                 ! Steps
    double precision, dimension(:), pointer   :: Qmu                ! Shear Attenuation
    double precision, dimension(:,:), pointer :: Qtau_e             ! tau_epsilon
    double precision, dimension(:), pointer   :: Qomsb, Qomsb2      ! one_minus_sum_beta
    double precision, dimension(:,:), pointer :: Qfc, Qfc2          ! factor_common
    double precision, dimension(:), pointer   :: Qsf, Qsf2          ! scale_factor
    integer, dimension(:), pointer            :: Qrmin              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: Qrmax              ! Max and Mins of idoubling
    integer                                   :: Qn                 ! Number of points
  end type attenuation_model_variables

!> Hejun
!  type (attenuation_model_variables) AM_V
! attenuation_model_variables
!< Hejun

! for forward or backward simulations
  logical COMPUTE_AND_STORE_STRAIN

! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: displ,accel

! for attenuation
! memory variables R_ij are stored at the local rather than global level
! to allow for optimization of cache access by compiler
  integer i_sls,i_memory
  real(kind=CUSTOM_REAL) R_xx_val,R_yy_val

! variable lengths for factor_common and one_minus_sum_beta
  integer vx, vy, vz, vnspec

  real(kind=CUSTOM_REAL), dimension(vx, vy, vz, vnspec) :: one_minus_sum_beta

  real(kind=CUSTOM_REAL), dimension(N_SLS, vx, vy, vz, vnspec) :: factor_common
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ) :: factor_common_use

  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ATTENUATION) :: R_memory
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: epsilondev
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: epsilon_trace_over_3

! array with the local to global mapping per slice
  integer, dimension(NSPEC_INNER_CORE) :: idoubling

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: xix,xiy,xiz, &
                      etax,etay,etaz,gammax,gammay,gammaz

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: kappavstore,muvstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: &
    c11store,c33store,c12store,c13store,c44store

  integer ispec,iglob,ispec_strain
  integer i,j,k !,l

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz

  !real(kind=CUSTOM_REAL) hp1,hp2,hp3
  real(kind=CUSTOM_REAL) fac1,fac2,fac3
  real(kind=CUSTOM_REAL) lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) kappal
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sum_terms

  real(kind=CUSTOM_REAL) minus_sum_beta
  real(kind=CUSTOM_REAL) c11l,c33l,c12l,c13l,c44l

! for gravity
  integer int_radius
  real(kind=CUSTOM_REAL) sigma_yx,sigma_zx,sigma_zy
  double precision radius,rho,minus_g,minus_dg
  double precision minus_g_over_radius,minus_dg_plus_g_over_radius
  double precision cos_theta,sin_theta,cos_phi,sin_phi
  double precision cos_theta_sq,sin_theta_sq,cos_phi_sq,sin_phi_sq
  double precision theta,phi,factor,gxl,gyl,gzl,sx_l,sy_l,sz_l
  double precision Hxxl,Hyyl,Hzzl,Hxyl,Hxzl,Hyzl
  double precision, dimension(NRAD_GRAVITY) :: minus_gravity_table,density_table,minus_deriv_gravity_table
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: rho_s_H
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube
  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE) :: xstore,ystore,zstore
!> Hejun
!  integer iregion_selected
!  real(kind=CUSTOM_REAL) radius_cr
!< Hejun

! Deville
! manually inline the calls to the Deville et al. (2002) routines
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc, &
    newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xxT,hprimewgll_xxT

  integer, parameter :: m1 = NGLLX, m2 = NGLLX * NGLLY

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

  real(kind=CUSTOM_REAL), dimension(m2,NGLLX) :: A1_mxm_m2_m1_5points,A2_mxm_m2_m1_5points,A3_mxm_m2_m1_5points
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: C1_mxm_m2_m1_5points,C2_mxm_m2_m1_5points,C3_mxm_m2_m1_5points
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: E1_mxm_m2_m1_5points,E2_mxm_m2_m1_5points,E3_mxm_m2_m1_5points

  equivalence(dummyx_loc,A1_mxm_m2_m1_5points)
  equivalence(dummyy_loc,A2_mxm_m2_m1_5points)
  equivalence(dummyz_loc,A3_mxm_m2_m1_5points)
  equivalence(tempx3,C1_mxm_m2_m1_5points)
  equivalence(tempy3,C2_mxm_m2_m1_5points)
  equivalence(tempz3,C3_mxm_m2_m1_5points)
  equivalence(newtempx3,E1_mxm_m2_m1_5points)
  equivalence(newtempy3,E2_mxm_m2_m1_5points)
  equivalence(newtempz3,E3_mxm_m2_m1_5points)

 

! ****************************************************
!   big loop over all spectral elements in the solid
! ****************************************************


! set acceleration to zero
  accel(:,:) = 0._CUSTOM_REAL

  do ispec = 1,NSPEC_INNER_CORE

    ! exclude fictitious elements in central cube
    if(idoubling(ispec) /= IFLAG_IN_FICTITIOUS_CUBE) then

      ! subroutines adapted from Deville, Fischer and Mund, High-order methods
      ! for incompressible fluid flow, Cambridge University Press (2002),
      ! pages 386 and 389 and Figure 8.3.1
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

            ! compute deviatoric strain
            if (COMPUTE_AND_STORE_STRAIN) then
              if(NSPEC_INNER_CORE_STRAIN_ONLY == 1) then
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
!> Hejun
            if(ATTENUATION_VAL) then
!              if(ATTENUATION_3D_VAL) then
                minus_sum_beta =  one_minus_sum_beta(i,j,k,ispec) - 1.0
!              else
!                radius_cr = xstore(ibool(i,j,k,ispec))
!                call get_attenuation_index(idoubling(ispec), dble(radius_cr), iregion_selected, .TRUE., AM_V)
!                minus_sum_beta =  one_minus_sum_beta(1,1,1,iregion_selected) - 1.0
!              endif ! ATTENUATION_3D_VAL
            endif ! ATTENUATION_VAL
!< Hejun
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
!> Hejun
              ! use unrelaxed parameters if attenuation
              if(ATTENUATION_VAL) then
!                if(ATTENUATION_3D_VAL) then
                  mul = mul * one_minus_sum_beta(i,j,k,ispec)
!                else
!                  mul = mul * one_minus_sum_beta(1,1,1,iregion_selected)
!                endif
              endif
!< Hejun
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
            if(ATTENUATION_VAL) then
              do i_sls = 1,N_SLS
                R_xx_val = R_memory(1,i_sls,i,j,k,ispec)
                R_yy_val = R_memory(2,i_sls,i,j,k,ispec)
                sigma_xx = sigma_xx - R_xx_val
                sigma_yy = sigma_yy - R_yy_val
                sigma_zz = sigma_zz + R_xx_val + R_yy_val
                sigma_xy = sigma_xy - R_memory(3,i_sls,i,j,k,ispec)
                sigma_xz = sigma_xz - R_memory(4,i_sls,i,j,k,ispec)
                sigma_yz = sigma_yz - R_memory(5,i_sls,i,j,k,ispec)
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
              iglob = ibool(i,j,k,ispec)
              radius = dble(xstore(iglob))
              theta = dble(ystore(iglob))
              phi = dble(zstore(iglob))

              ! make sure radius is never zero even for points at center of cube
              ! because we later divide by radius
              if(radius < 100.d0 / R_EARTH) radius = 100.d0 / R_EARTH

              cos_theta = dcos(theta)
              sin_theta = dsin(theta)
              cos_phi = dcos(phi)
              sin_phi = dsin(phi)

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

              cos_theta_sq = cos_theta**2
              sin_theta_sq = sin_theta**2
              cos_phi_sq = cos_phi**2
              sin_phi_sq = sin_phi**2

              Hxxl = minus_g_over_radius*(cos_phi_sq*cos_theta_sq + sin_phi_sq) + cos_phi_sq*minus_dg*sin_theta_sq
              Hyyl = minus_g_over_radius*(cos_phi_sq + cos_theta_sq*sin_phi_sq) + minus_dg*sin_phi_sq*sin_theta_sq
              Hzzl = cos_theta_sq*minus_dg + minus_g_over_radius*sin_theta_sq
              Hxyl = cos_phi*minus_dg_plus_g_over_radius*sin_phi*sin_theta_sq
              Hxzl = cos_phi*cos_theta*minus_dg_plus_g_over_radius*sin_theta
              Hyzl = cos_theta*minus_dg_plus_g_over_radius*sin_phi*sin_theta

              iglob = ibool(i,j,k,ispec)

              ! distinguish between single and double precision for reals
              if(CUSTOM_REAL == SIZE_REAL) then
                ! get displacement and multiply by density to compute G tensor
                sx_l = rho * dble(displ(1,iglob))
                sy_l = rho * dble(displ(2,iglob))
                sz_l = rho * dble(displ(3,iglob))

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
                sx_l = rho * displ(1,iglob)
                sy_l = rho * displ(2,iglob)
                sz_l = rho * displ(3,iglob)

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
          do i=1,NGLLX
            iglob = ibool(i,j,k,ispec)
            accel(:,iglob) = accel(:,iglob) + sum_terms(:,i,j,k)
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
      if(ATTENUATION_VAL) then
!> Hejun      
        do i_sls = 1,N_SLS
!          if(ATTENUATION_3D_VAL) then
             factor_common_use = factor_common(i_sls,:,:,:,ispec)
!          else
!             factor_common_use(:,:,:) = factor_common(i_sls,1,1,1,iregion_selected)
!          endif
!< Hejun
          do i_memory = 1,5
             R_memory(i_memory,i_sls,:,:,:,ispec) = &
                  alphaval(i_sls) * &
                  R_memory(i_memory,i_sls,:,:,:,ispec) + muvstore(:,:,:,ispec) * &
                  factor_common_use * &
                  (betaval(i_sls) * &
                  epsilondev(i_memory,:,:,:,ispec) + gammaval(i_sls) * epsilondev_loc(i_memory,:,:,:))
          enddo
        enddo

      endif

      ! save deviatoric strain for Runge-Kutta scheme
      if(COMPUTE_AND_STORE_STRAIN) epsilondev(:,:,:,:,ispec) = epsilondev_loc(:,:,:,:)

    endif   ! end test to exclude fictitious elements in central cube

  enddo ! spectral element loop

  end subroutine compute_forces_inner_core_Dev

