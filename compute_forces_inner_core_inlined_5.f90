!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology August 2003
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine compute_forces_inner_core(minus_gravity_table,density_table,minus_deriv_gravity_table, &
          displ,accel,xstore,ystore,zstore, &
          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian, &
          hprime_xx,hprime_yy,hprime_zz, &
          hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
          wgll_cube,wgllwgll_yz_no_i,wgllwgll_xz_no_j,wgllwgll_xy_no_k, &
          kappavstore,muvstore,ibool,idoubling, &
          c11store,c33store,c12store,c13store,c44store,R_memory,epsilondev, &
! BS
! BS added more variables, with variable lengths
!          one_minus_sum_beta,alphaval,betaval,gammaval,factor_common,index_i,index_k,index_dim)
          one_minus_sum_beta,&
          alphaval,betaval,gammaval, &
          factor_common, &
          vx, vy, vz, vnspec, &
          index_i, index_k, index_dim)
! BS END

  implicit none

  include "constants.h"

! include values created by the mesher
! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: displ,accel

! for attenuation
! memory variables R_ij are stored at the local rather than global level
! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL) R_xx_val,R_yy_val

! BS
! BS variable lengths for factor_common and one_minus_sum_beta
  integer vx, vy, vz, vnspec
! BS END

! BS
!   real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION) :: one_minus_sum_beta
  real(kind=CUSTOM_REAL), dimension(vx, vy, vz, vnspec) :: one_minus_sum_beta
! BS END

! BS
! real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION,N_SLS) :: alphaval,betaval,gammaval,factor_common
  real(kind=CUSTOM_REAL), dimension(vx, vy, vz, vnspec, N_SLS) :: factor_common
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: factor_common_use
! BS END

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ATTENUATION,N_SLS) :: R_memory
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ATTENUATION) :: epsilondev
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc
  real(kind=CUSTOM_REAL) epsilon_trace_over_3

! array with the local to global mapping per slice
  integer, dimension(NSPEC_INNER_CORE) :: idoubling

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: xix,xiy,xiz, &
                      etax,etay,etaz,gammax,gammay,gammaz,jacobian

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: &
    temp1,temp2,temp3,temp1inline,temp2inline,temp3inline
  real(kind=CUSTOM_REAL) temp11,temp12,temp13,temp21,temp22,temp23,temp31,temp32,temp33

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: kappavstore,muvstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: &
    c11store,c33store,c12store,c13store,c44store

  integer ispec,iglob
  integer i,j,k,d,ijd,ijk
  integer, dimension(NGLLSQUARE_NDIM) :: index_i,index_k,index_dim

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz

  real(kind=CUSTOM_REAL) lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) kappal
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sum_terms

  real(kind=CUSTOM_REAL) minus_sum_beta
  real(kind=CUSTOM_REAL) c11l,c33l,c12l,c13l,c44l

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: disploc

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
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: wgllwgll_yz_no_i,wgllwgll_xz_no_j,wgllwgll_xy_no_k
  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE) :: xstore,ystore,zstore

! ****************************************************
!   big loop over all spectral elements in the solid
! ****************************************************

! BS
! BS removed line, handled later
!  minus_sum_beta =  one_minus_sum_beta(IREGION_ATTENUATION_INNER_CORE) - 1.
! BS END

! set acceleration to zero
  accel(:,:) = 0._CUSTOM_REAL

!CDIR NOVECTOR
  do ispec = 1,NSPEC_INNER_CORE

! exclude fictitious elements in central cube
    if(idoubling(ispec) /= IFLAG_IN_FICTITIOUS_CUBE) then

! copy global displacement to local
    do ijk=1,NGLLCUBE
      iglob = ibool(ijk,1,1,ispec)
      disploc(1,ijk,1,1) = displ(1,iglob)
      disploc(2,ijk,1,1) = displ(2,iglob)
      disploc(3,ijk,1,1) = displ(3,iglob)
    enddo

! inlined first matrix product

!CDIR NODEP(temp1,temp2)
  do ijd = 1,NGLLSQUARE_NDIM

    i = index_i(ijd)
    j = i
    k = index_k(ijd)
    d = index_dim(ijd)

!---
  temp1(d,1,j,k) = &
    disploc(d,1,j,k)*hprime_xx(1,1) + disploc(d,2,j,k)*hprime_xx(2,1) + &
    disploc(d,3,j,k)*hprime_xx(3,1) + disploc(d,4,j,k)*hprime_xx(4,1) + &
    disploc(d,5,j,k)*hprime_xx(5,1)

!---
  temp2(d,i,1,k) = &
    disploc(d,i,1,k)*hprime_yy(1,1) + disploc(d,i,2,k)*hprime_yy(2,1) + &
    disploc(d,i,3,k)*hprime_yy(3,1) + disploc(d,i,4,k)*hprime_yy(4,1) + &
    disploc(d,i,5,k)*hprime_yy(5,1)

!---
  temp3(ijd,1,1,1) = &
    disploc(ijd,1,1,1)*hprime_zz(1,1) + disploc(ijd,1,1,2)*hprime_zz(2,1) + &
    disploc(ijd,1,1,3)*hprime_zz(3,1) + disploc(ijd,1,1,4)*hprime_zz(4,1) + &
    disploc(ijd,1,1,5)*hprime_zz(5,1)

!---
  temp1(d,2,j,k) = &
    disploc(d,1,j,k)*hprime_xx(1,2) + disploc(d,2,j,k)*hprime_xx(2,2) + &
    disploc(d,3,j,k)*hprime_xx(3,2) + disploc(d,4,j,k)*hprime_xx(4,2) + &
    disploc(d,5,j,k)*hprime_xx(5,2)

!---
  temp2(d,i,2,k) = &
    disploc(d,i,1,k)*hprime_yy(1,2) + disploc(d,i,2,k)*hprime_yy(2,2) + &
    disploc(d,i,3,k)*hprime_yy(3,2) + disploc(d,i,4,k)*hprime_yy(4,2) + &
    disploc(d,i,5,k)*hprime_yy(5,2)

!---
  temp3(ijd,1,1,2) = &
    disploc(ijd,1,1,1)*hprime_zz(1,2) + disploc(ijd,1,1,2)*hprime_zz(2,2) + &
    disploc(ijd,1,1,3)*hprime_zz(3,2) + disploc(ijd,1,1,4)*hprime_zz(4,2) + &
    disploc(ijd,1,1,5)*hprime_zz(5,2)

!---
  temp1(d,3,j,k) = &
    disploc(d,1,j,k)*hprime_xx(1,3) + disploc(d,2,j,k)*hprime_xx(2,3) + &
    disploc(d,3,j,k)*hprime_xx(3,3) + disploc(d,4,j,k)*hprime_xx(4,3) + &
    disploc(d,5,j,k)*hprime_xx(5,3)

!---
  temp2(d,i,3,k) = &
    disploc(d,i,1,k)*hprime_yy(1,3) + disploc(d,i,2,k)*hprime_yy(2,3) + &
    disploc(d,i,3,k)*hprime_yy(3,3) + disploc(d,i,4,k)*hprime_yy(4,3) + &
    disploc(d,i,5,k)*hprime_yy(5,3)

!---
  temp3(ijd,1,1,3) = &
    disploc(ijd,1,1,1)*hprime_zz(1,3) + disploc(ijd,1,1,2)*hprime_zz(2,3) + &
    disploc(ijd,1,1,3)*hprime_zz(3,3) + disploc(ijd,1,1,4)*hprime_zz(4,3) + &
    disploc(ijd,1,1,5)*hprime_zz(5,3)

!---
  temp1(d,4,j,k) = &
    disploc(d,1,j,k)*hprime_xx(1,4) + disploc(d,2,j,k)*hprime_xx(2,4) + &
    disploc(d,3,j,k)*hprime_xx(3,4) + disploc(d,4,j,k)*hprime_xx(4,4) + &
    disploc(d,5,j,k)*hprime_xx(5,4)

!---
  temp2(d,i,4,k) = &
    disploc(d,i,1,k)*hprime_yy(1,4) + disploc(d,i,2,k)*hprime_yy(2,4) + &
    disploc(d,i,3,k)*hprime_yy(3,4) + disploc(d,i,4,k)*hprime_yy(4,4) + &
    disploc(d,i,5,k)*hprime_yy(5,4)

!---
  temp3(ijd,1,1,4) = &
    disploc(ijd,1,1,1)*hprime_zz(1,4) + disploc(ijd,1,1,2)*hprime_zz(2,4) + &
    disploc(ijd,1,1,3)*hprime_zz(3,4) + disploc(ijd,1,1,4)*hprime_zz(4,4) + &
    disploc(ijd,1,1,5)*hprime_zz(5,4)

!---
  temp1(d,5,j,k) = &
    disploc(d,1,j,k)*hprime_xx(1,5) + disploc(d,2,j,k)*hprime_xx(2,5) + &
    disploc(d,3,j,k)*hprime_xx(3,5) + disploc(d,4,j,k)*hprime_xx(4,5) + &
    disploc(d,5,j,k)*hprime_xx(5,5)

!---
  temp2(d,i,5,k) = &
    disploc(d,i,1,k)*hprime_yy(1,5) + disploc(d,i,2,k)*hprime_yy(2,5) + &
    disploc(d,i,3,k)*hprime_yy(3,5) + disploc(d,i,4,k)*hprime_yy(4,5) + &
    disploc(d,i,5,k)*hprime_yy(5,5)

!---
  temp3(ijd,1,1,5) = &
    disploc(ijd,1,1,1)*hprime_zz(1,5) + disploc(ijd,1,1,2)*hprime_zz(2,5) + &
    disploc(ijd,1,1,3)*hprime_zz(3,5) + disploc(ijd,1,1,4)*hprime_zz(4,5) + &
    disploc(ijd,1,1,5)*hprime_zz(5,5)

  enddo

! fused the three loops for inlined version
    do ijk=1,NGLLCUBE

!         get derivatives of ux, uy and uz with respect to x, y and z

          xixl = xix(ijk,1,1,ispec)
          xiyl = xiy(ijk,1,1,ispec)
          xizl = xiz(ijk,1,1,ispec)
          etaxl = etax(ijk,1,1,ispec)
          etayl = etay(ijk,1,1,ispec)
          etazl = etaz(ijk,1,1,ispec)
          gammaxl = gammax(ijk,1,1,ispec)
          gammayl = gammay(ijk,1,1,ispec)
          gammazl = gammaz(ijk,1,1,ispec)
          jacobianl = jacobian(ijk,1,1,ispec)

          temp11 = temp1(1,ijk,1,1)
          temp12 = temp1(2,ijk,1,1)
          temp13 = temp1(3,ijk,1,1)

          temp21 = temp2(1,ijk,1,1)
          temp22 = temp2(2,ijk,1,1)
          temp23 = temp2(3,ijk,1,1)

          temp31 = temp3(1,ijk,1,1)
          temp32 = temp3(2,ijk,1,1)
          temp33 = temp3(3,ijk,1,1)

          duxdxl = xixl*temp11 + etaxl*temp21 + gammaxl*temp31
          duxdyl = xiyl*temp11 + etayl*temp21 + gammayl*temp31
          duxdzl = xizl*temp11 + etazl*temp21 + gammazl*temp31

          duydxl = xixl*temp12 + etaxl*temp22 + gammaxl*temp32
          duydyl = xiyl*temp12 + etayl*temp22 + gammayl*temp32
          duydzl = xizl*temp12 + etazl*temp22 + gammazl*temp32

          duzdxl = xixl*temp13 + etaxl*temp23 + gammaxl*temp33
          duzdyl = xiyl*temp13 + etayl*temp23 + gammayl*temp33
          duzdzl = xizl*temp13 + etazl*temp23 + gammazl*temp33

! precompute some sums to save CPU time
          duxdxl_plus_duydyl = duxdxl + duydyl
          duxdxl_plus_duzdzl = duxdxl + duzdzl
          duydyl_plus_duzdzl = duydyl + duzdzl
          duxdyl_plus_duydxl = duxdyl + duydxl
          duzdxl_plus_duxdzl = duzdxl + duxdzl
          duzdyl_plus_duydzl = duzdyl + duydzl

! compute deviatoric strain
  if(ATTENUATION_VAL) then
! BS
     if(ATTENUATION_VAL_3D) then
        minus_sum_beta =  one_minus_sum_beta(ijk,1,1,ispec) - 1.0
     else
        minus_sum_beta =  one_minus_sum_beta(1,1,1,IREGION_ATTENUATION_INNER_CORE) - 1.
     endif
! BS END
    epsilon_trace_over_3 = ONE_THIRD * (duxdxl + duydyl + duzdzl)
    epsilondev_loc(1,ijk,1,1) = duxdxl - epsilon_trace_over_3
    epsilondev_loc(2,ijk,1,1) = duydyl - epsilon_trace_over_3
    epsilondev_loc(3,ijk,1,1) = 0.5 * duxdyl_plus_duydxl
    epsilondev_loc(4,ijk,1,1) = 0.5 * duzdxl_plus_duxdzl
    epsilondev_loc(5,ijk,1,1) = 0.5 * duzdyl_plus_duydzl
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

         c11l = c11store(ijk,1,1,ispec)
         c12l = c12store(ijk,1,1,ispec)
         c13l = c13store(ijk,1,1,ispec)
         c33l = c33store(ijk,1,1,ispec)
         c44l = c44store(ijk,1,1,ispec)

! use unrelaxed parameters if attenuation
         if(ATTENUATION_VAL) then
           mul = muvstore(ijk,1,1,ispec)
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
          kappal = kappavstore(ijk,1,1,ispec)
          mul = muvstore(ijk,1,1,ispec)

! use unrelaxed parameters if attenuation
! BS
!  if(ATTENUATION_VAL) mul = mul * one_minus_sum_beta(IREGION_ATTENUATION_INNER_CORE)
  if(ATTENUATION_VAL) then
     if(ATTENUATION_VAL_3D) then
        mul = mul * one_minus_sum_beta(ijk,1,1,ispec)
     else
        mul = mul * one_minus_sum_beta(1,1,1,1)
     endif
  endif
! BS END
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

    R_xx_val = R_memory(1,ijk,1,1,ispec,1)
    R_yy_val = R_memory(2,ijk,1,1,ispec,1)
    sigma_xx = sigma_xx - R_xx_val
    sigma_yy = sigma_yy - R_yy_val
    sigma_zz = sigma_zz + R_xx_val + R_yy_val
    sigma_xy = sigma_xy - R_memory(3,ijk,1,1,ispec,1)
    sigma_xz = sigma_xz - R_memory(4,ijk,1,1,ispec,1)
    sigma_yz = sigma_yz - R_memory(5,ijk,1,1,ispec,1)

    R_xx_val = R_memory(1,ijk,1,1,ispec,2)
    R_yy_val = R_memory(2,ijk,1,1,ispec,2)
    sigma_xx = sigma_xx - R_xx_val
    sigma_yy = sigma_yy - R_yy_val
    sigma_zz = sigma_zz + R_xx_val + R_yy_val
    sigma_xy = sigma_xy - R_memory(3,ijk,1,1,ispec,2)
    sigma_xz = sigma_xz - R_memory(4,ijk,1,1,ispec,2)
    sigma_yz = sigma_yz - R_memory(5,ijk,1,1,ispec,2)

    R_xx_val = R_memory(1,ijk,1,1,ispec,3)
    R_yy_val = R_memory(2,ijk,1,1,ispec,3)
    sigma_xx = sigma_xx - R_xx_val
    sigma_yy = sigma_yy - R_yy_val
    sigma_zz = sigma_zz + R_xx_val + R_yy_val
    sigma_xy = sigma_xy - R_memory(3,ijk,1,1,ispec,3)
    sigma_xz = sigma_xz - R_memory(4,ijk,1,1,ispec,3)
    sigma_yz = sigma_yz - R_memory(5,ijk,1,1,ispec,3)

  endif

! define symmetric components of sigma for gravity
  sigma_yx = sigma_xy
  sigma_zx = sigma_xz
  sigma_zy = sigma_yz

! compute non-symmetric terms for gravity
  if(GRAVITY_VAL) then

! use mesh coordinates to get theta and phi
! x y and z contain r theta and phi

  iglob = ibool(ijk,1,1,ispec)
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

! distinguish whether single or double precision for reals
  if(CUSTOM_REAL == SIZE_REAL) then

! get displacement and multiply by density to compute G tensor
    sx_l = rho * dble(disploc(1,ijk,1,1))
    sy_l = rho * dble(disploc(2,ijk,1,1))
    sz_l = rho * dble(disploc(3,ijk,1,1))

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
    factor = dble(jacobianl) * wgll_cube(ijk,1,1)
    rho_s_H(1,ijk,1,1) = sngl(factor * (sx_l * Hxxl + sy_l * Hxyl + sz_l * Hxzl))
    rho_s_H(2,ijk,1,1) = sngl(factor * (sx_l * Hxyl + sy_l * Hyyl + sz_l * Hyzl))
    rho_s_H(3,ijk,1,1) = sngl(factor * (sx_l * Hxzl + sy_l * Hyzl + sz_l * Hzzl))

  else

! get displacement and multiply by density to compute G tensor
    sx_l = rho * disploc(1,ijk,1,1)
    sy_l = rho * disploc(2,ijk,1,1)
    sz_l = rho * disploc(3,ijk,1,1)

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
    factor = jacobianl * wgll_cube(ijk,1,1)
    rho_s_H(1,ijk,1,1) = factor * (sx_l * Hxxl + sy_l * Hxyl + sz_l * Hxzl)
    rho_s_H(2,ijk,1,1) = factor * (sx_l * Hxyl + sy_l * Hyyl + sz_l * Hyzl)
    rho_s_H(3,ijk,1,1) = factor * (sx_l * Hxzl + sy_l * Hyzl + sz_l * Hzzl)

  endif

  endif  ! end of section with gravity terms

! form dot product with test vector, non-symmetric form
      temp1(1,ijk,1,1) = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl)
      temp1(2,ijk,1,1) = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl)
      temp1(3,ijk,1,1) = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl)

      temp2(1,ijk,1,1) = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl)
      temp2(2,ijk,1,1) = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl)
      temp2(3,ijk,1,1) = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl)

      temp3(1,ijk,1,1) = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl)
      temp3(2,ijk,1,1) = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl)
      temp3(3,ijk,1,1) = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl)

      enddo

! inlined second matrix product

!CDIR NODEP(temp1inline,temp2inline)
  do ijd = 1,NGLLSQUARE_NDIM

    i = index_i(ijd)
    j = i
    k = index_k(ijd)
    d = index_dim(ijd)

!---
  temp1inline(d,1,j,k) = &
    temp1(d,1,j,k)*hprimewgll_xx(1,1) + temp1(d,2,j,k)*hprimewgll_xx(1,2) + &
    temp1(d,3,j,k)*hprimewgll_xx(1,3) + temp1(d,4,j,k)*hprimewgll_xx(1,4) + &
    temp1(d,5,j,k)*hprimewgll_xx(1,5)

!---
  temp2inline(d,i,1,k) = &
    temp2(d,i,1,k)*hprimewgll_yy(1,1) + temp2(d,i,2,k)*hprimewgll_yy(1,2) + &
    temp2(d,i,3,k)*hprimewgll_yy(1,3) + temp2(d,i,4,k)*hprimewgll_yy(1,4) + &
    temp2(d,i,5,k)*hprimewgll_yy(1,5)

!---
  temp3inline(ijd,1,1,1) = &
   temp3(ijd,1,1,1)*hprimewgll_zz(1,1) + temp3(ijd,1,1,2)*hprimewgll_zz(1,2) + &
   temp3(ijd,1,1,3)*hprimewgll_zz(1,3) + temp3(ijd,1,1,4)*hprimewgll_zz(1,4) + &
   temp3(ijd,1,1,5)*hprimewgll_zz(1,5)

!---
  temp1inline(d,2,j,k) = &
    temp1(d,1,j,k)*hprimewgll_xx(2,1) + temp1(d,2,j,k)*hprimewgll_xx(2,2) + &
    temp1(d,3,j,k)*hprimewgll_xx(2,3) + temp1(d,4,j,k)*hprimewgll_xx(2,4) + &
    temp1(d,5,j,k)*hprimewgll_xx(2,5)

!---
  temp2inline(d,i,2,k) = &
    temp2(d,i,1,k)*hprimewgll_yy(2,1) + temp2(d,i,2,k)*hprimewgll_yy(2,2) + &
    temp2(d,i,3,k)*hprimewgll_yy(2,3) + temp2(d,i,4,k)*hprimewgll_yy(2,4) + &
    temp2(d,i,5,k)*hprimewgll_yy(2,5)

!---
  temp3inline(ijd,1,1,2) = &
   temp3(ijd,1,1,1)*hprimewgll_zz(2,1) + temp3(ijd,1,1,2)*hprimewgll_zz(2,2) + &
   temp3(ijd,1,1,3)*hprimewgll_zz(2,3) + temp3(ijd,1,1,4)*hprimewgll_zz(2,4) + &
   temp3(ijd,1,1,5)*hprimewgll_zz(2,5)

!---
  temp1inline(d,3,j,k) = &
    temp1(d,1,j,k)*hprimewgll_xx(3,1) + temp1(d,2,j,k)*hprimewgll_xx(3,2) + &
    temp1(d,3,j,k)*hprimewgll_xx(3,3) + temp1(d,4,j,k)*hprimewgll_xx(3,4) + &
    temp1(d,5,j,k)*hprimewgll_xx(3,5)

!---
  temp2inline(d,i,3,k) = &
    temp2(d,i,1,k)*hprimewgll_yy(3,1) + temp2(d,i,2,k)*hprimewgll_yy(3,2) + &
    temp2(d,i,3,k)*hprimewgll_yy(3,3) + temp2(d,i,4,k)*hprimewgll_yy(3,4) + &
    temp2(d,i,5,k)*hprimewgll_yy(3,5)

!---
  temp3inline(ijd,1,1,3) = &
   temp3(ijd,1,1,1)*hprimewgll_zz(3,1) + temp3(ijd,1,1,2)*hprimewgll_zz(3,2) + &
   temp3(ijd,1,1,3)*hprimewgll_zz(3,3) + temp3(ijd,1,1,4)*hprimewgll_zz(3,4) + &
   temp3(ijd,1,1,5)*hprimewgll_zz(3,5)

!---
  temp1inline(d,4,j,k) = &
    temp1(d,1,j,k)*hprimewgll_xx(4,1) + temp1(d,2,j,k)*hprimewgll_xx(4,2) + &
    temp1(d,3,j,k)*hprimewgll_xx(4,3) + temp1(d,4,j,k)*hprimewgll_xx(4,4) + &
    temp1(d,5,j,k)*hprimewgll_xx(4,5)

!---
  temp2inline(d,i,4,k) = &
    temp2(d,i,1,k)*hprimewgll_yy(4,1) + temp2(d,i,2,k)*hprimewgll_yy(4,2) + &
    temp2(d,i,3,k)*hprimewgll_yy(4,3) + temp2(d,i,4,k)*hprimewgll_yy(4,4) + &
    temp2(d,i,5,k)*hprimewgll_yy(4,5)

!---
  temp3inline(ijd,1,1,4) = &
   temp3(ijd,1,1,1)*hprimewgll_zz(4,1) + temp3(ijd,1,1,2)*hprimewgll_zz(4,2) + &
   temp3(ijd,1,1,3)*hprimewgll_zz(4,3) + temp3(ijd,1,1,4)*hprimewgll_zz(4,4) + &
   temp3(ijd,1,1,5)*hprimewgll_zz(4,5)

!---
  temp1inline(d,5,j,k) = &
    temp1(d,1,j,k)*hprimewgll_xx(5,1) + temp1(d,2,j,k)*hprimewgll_xx(5,2) + &
    temp1(d,3,j,k)*hprimewgll_xx(5,3) + temp1(d,4,j,k)*hprimewgll_xx(5,4) + &
    temp1(d,5,j,k)*hprimewgll_xx(5,5)

!---
  temp2inline(d,i,5,k) = &
    temp2(d,i,1,k)*hprimewgll_yy(5,1) + temp2(d,i,2,k)*hprimewgll_yy(5,2) + &
    temp2(d,i,3,k)*hprimewgll_yy(5,3) + temp2(d,i,4,k)*hprimewgll_yy(5,4) + &
    temp2(d,i,5,k)*hprimewgll_yy(5,5)

!---
  temp3inline(ijd,1,1,5) = &
   temp3(ijd,1,1,1)*hprimewgll_zz(5,1) + temp3(ijd,1,1,2)*hprimewgll_zz(5,2) + &
   temp3(ijd,1,1,3)*hprimewgll_zz(5,3) + temp3(ijd,1,1,4)*hprimewgll_zz(5,4) + &
   temp3(ijd,1,1,5)*hprimewgll_zz(5,5)

  enddo

  if(GRAVITY_VAL) then
    do ijk=1,NGLLCUBE_NDIM
      sum_terms(ijk,1,1,1) = &
           - (wgllwgll_yz_no_i(ijk,1,1,1)*temp1inline(ijk,1,1,1) + &
              wgllwgll_xz_no_j(ijk,1,1,1)*temp2inline(ijk,1,1,1) + &
              wgllwgll_xy_no_k(ijk,1,1,1)*temp3inline(ijk,1,1,1)) + &
              rho_s_H(ijk,1,1,1)
    enddo
  else
    do ijk=1,NGLLCUBE_NDIM
      sum_terms(ijk,1,1,1) = &
           - (wgllwgll_yz_no_i(ijk,1,1,1)*temp1inline(ijk,1,1,1) + &
              wgllwgll_xz_no_j(ijk,1,1,1)*temp2inline(ijk,1,1,1) + &
              wgllwgll_xy_no_k(ijk,1,1,1)*temp3inline(ijk,1,1,1))
    enddo
  endif

! sum contributions from each element to the global mesh and add gravity terms
!CDIR NODEP(accel)
    do ijk=1,NGLLCUBE
      iglob = ibool(ijk,1,1,ispec)
      accel(1,iglob) = accel(1,iglob) + sum_terms(1,ijk,1,1)
      accel(2,iglob) = accel(2,iglob) + sum_terms(2,ijk,1,1)
      accel(3,iglob) = accel(3,iglob) + sum_terms(3,ijk,1,1)
    enddo

! use Runge-Kutta scheme to march memory variables in time
! convention for attenuation
! term in xx = 1
! term in yy = 2
! term in xy = 3
! term in xz = 4
! term in yz = 5
! term in zz not computed since zero trace

  if(ATTENUATION_VAL) then

    do ijk = 1,NGLLCUBE

!---
! BS 
!    R_memory(1,ijk,1,1,ispec,1) = &
!      alphaval(IREGION_ATTENUATION_INNER_CORE,1) * &
!      R_memory(1,ijk,1,1,ispec,1) + muvstore(ijk,1,1,ispec) * &
!      factor_common(IREGION_ATTENUATION_INNER_CORE,1) * &
!      (betaval(IREGION_ATTENUATION_INNER_CORE,1) * &
!      epsilondev(1,ijk,1,1,ispec) + gammaval(IREGION_ATTENUATION_INNER_CORE,1) * epsilondev_loc(1,ijk,1,1))
!
!    R_memory(2,ijk,1,1,ispec,1) = &
!      alphaval(IREGION_ATTENUATION_INNER_CORE,1) * &
!      R_memory(2,ijk,1,1,ispec,1) + muvstore(ijk,1,1,ispec) * &
!      factor_common(IREGION_ATTENUATION_INNER_CORE,1) * &
!      (betaval(IREGION_ATTENUATION_INNER_CORE,1) * &
!      epsilondev(2,ijk,1,1,ispec) + gammaval(IREGION_ATTENUATION_INNER_CORE,1) * epsilondev_loc(2,ijk,1,1))
!
!    R_memory(3,ijk,1,1,ispec,1) = &
!      alphaval(IREGION_ATTENUATION_INNER_CORE,1) * &
!      R_memory(3,ijk,1,1,ispec,1) + muvstore(ijk,1,1,ispec) * &
!      factor_common(IREGION_ATTENUATION_INNER_CORE,1) * &
!      (betaval(IREGION_ATTENUATION_INNER_CORE,1) * &
!      epsilondev(3,ijk,1,1,ispec) + gammaval(IREGION_ATTENUATION_INNER_CORE,1) * epsilondev_loc(3,ijk,1,1))
!
!    R_memory(4,ijk,1,1,ispec,1) = &
!      alphaval(IREGION_ATTENUATION_INNER_CORE,1) * &
!      R_memory(4,ijk,1,1,ispec,1) + muvstore(ijk,1,1,ispec) * &
!      factor_common(IREGION_ATTENUATION_INNER_CORE,1) * &
!      (betaval(IREGION_ATTENUATION_INNER_CORE,1) * &
!      epsilondev(4,ijk,1,1,ispec) + gammaval(IREGION_ATTENUATION_INNER_CORE,1) * epsilondev_loc(4,ijk,1,1))
!
!    R_memory(5,ijk,1,1,ispec,1) = &
!      alphaval(IREGION_ATTENUATION_INNER_CORE,1) * &
!      R_memory(5,ijk,1,1,ispec,1) + muvstore(ijk,1,1,ispec) * &
!      factor_common(IREGION_ATTENUATION_INNER_CORE,1) * &
!      (betaval(IREGION_ATTENUATION_INNER_CORE,1) * &
!      epsilondev(5,ijk,1,1,ispec) + gammaval(IREGION_ATTENUATION_INNER_CORE,1) * epsilondev_loc(5,ijk,1,1))
!
!!---
!    R_memory(1,ijk,1,1,ispec,2) = &
!      alphaval(IREGION_ATTENUATION_INNER_CORE,2) * &
!      R_memory(1,ijk,1,1,ispec,2) + muvstore(ijk,1,1,ispec) * &
!      factor_common(IREGION_ATTENUATION_INNER_CORE,2) * &
!      (betaval(IREGION_ATTENUATION_INNER_CORE,2) * &
!      epsilondev(1,ijk,1,1,ispec) + gammaval(IREGION_ATTENUATION_INNER_CORE,2) * epsilondev_loc(1,ijk,1,1))
!
!    R_memory(2,ijk,1,1,ispec,2) = &
!      alphaval(IREGION_ATTENUATION_INNER_CORE,2) * &
!      R_memory(2,ijk,1,1,ispec,2) + muvstore(ijk,1,1,ispec) * &
!      factor_common(IREGION_ATTENUATION_INNER_CORE,2) * &
!      (betaval(IREGION_ATTENUATION_INNER_CORE,2) * &
!      epsilondev(2,ijk,1,1,ispec) + gammaval(IREGION_ATTENUATION_INNER_CORE,2) * epsilondev_loc(2,ijk,1,1))
!
!    R_memory(3,ijk,1,1,ispec,2) = &
!      alphaval(IREGION_ATTENUATION_INNER_CORE,2) * &
!      R_memory(3,ijk,1,1,ispec,2) + muvstore(ijk,1,1,ispec) * &
!      factor_common(IREGION_ATTENUATION_INNER_CORE,2) * &
!      (betaval(IREGION_ATTENUATION_INNER_CORE,2) * &
!      epsilondev(3,ijk,1,1,ispec) + gammaval(IREGION_ATTENUATION_INNER_CORE,2) * epsilondev_loc(3,ijk,1,1))
!
!    R_memory(4,ijk,1,1,ispec,2) = &
!      alphaval(IREGION_ATTENUATION_INNER_CORE,2) * &
!      R_memory(4,ijk,1,1,ispec,2) + muvstore(ijk,1,1,ispec) * &
!      factor_common(IREGION_ATTENUATION_INNER_CORE,2) * &
!      (betaval(IREGION_ATTENUATION_INNER_CORE,2) * &
!      epsilondev(4,ijk,1,1,ispec) + gammaval(IREGION_ATTENUATION_INNER_CORE,2) * epsilondev_loc(4,ijk,1,1))
!
!    R_memory(5,ijk,1,1,ispec,2) = &
!      alphaval(IREGION_ATTENUATION_INNER_CORE,2) * &
!      R_memory(5,ijk,1,1,ispec,2) + muvstore(ijk,1,1,ispec) * &
!      factor_common(IREGION_ATTENUATION_INNER_CORE,2) * &
!      (betaval(IREGION_ATTENUATION_INNER_CORE,2) * &
!      epsilondev(5,ijk,1,1,ispec) + gammaval(IREGION_ATTENUATION_INNER_CORE,2) * epsilondev_loc(5,ijk,1,1))
!
!!---
!    R_memory(1,ijk,1,1,ispec,3) = &
!      alphaval(IREGION_ATTENUATION_INNER_CORE,3) * &
!      R_memory(1,ijk,1,1,ispec,3) + muvstore(ijk,1,1,ispec) * &
!      factor_common(IREGION_ATTENUATION_INNER_CORE,3) * &
!      (betaval(IREGION_ATTENUATION_INNER_CORE,3) * &
!      epsilondev(1,ijk,1,1,ispec) + gammaval(IREGION_ATTENUATION_INNER_CORE,3) * epsilondev_loc(1,ijk,1,1))
!
!    R_memory(2,ijk,1,1,ispec,3) = &
!      alphaval(IREGION_ATTENUATION_INNER_CORE,3) * &
!      R_memory(2,ijk,1,1,ispec,3) + muvstore(ijk,1,1,ispec) * &
!      factor_common(IREGION_ATTENUATION_INNER_CORE,3) * &
!      (betaval(IREGION_ATTENUATION_INNER_CORE,3) * &
!      epsilondev(2,ijk,1,1,ispec) + gammaval(IREGION_ATTENUATION_INNER_CORE,3) * epsilondev_loc(2,ijk,1,1))
!
!    R_memory(3,ijk,1,1,ispec,3) = &
!      alphaval(IREGION_ATTENUATION_INNER_CORE,3) * &
!      R_memory(3,ijk,1,1,ispec,3) + muvstore(ijk,1,1,ispec) * &
!      factor_common(IREGION_ATTENUATION_INNER_CORE,3) * &
!      (betaval(IREGION_ATTENUATION_INNER_CORE,3) * &
!      epsilondev(3,ijk,1,1,ispec) + gammaval(IREGION_ATTENUATION_INNER_CORE,3) * epsilondev_loc(3,ijk,1,1))
!
!    R_memory(4,ijk,1,1,ispec,3) = &
!      alphaval(IREGION_ATTENUATION_INNER_CORE,3) * &
!      R_memory(4,ijk,1,1,ispec,3) + muvstore(ijk,1,1,ispec) * &
!      factor_common(IREGION_ATTENUATION_INNER_CORE,3) * &
!      (betaval(IREGION_ATTENUATION_INNER_CORE,3) * &
!      epsilondev(4,ijk,1,1,ispec) + gammaval(IREGION_ATTENUATION_INNER_CORE,3) * epsilondev_loc(4,ijk,1,1))
!
!    R_memory(5,ijk,1,1,ispec,3) = &
!      alphaval(IREGION_ATTENUATION_INNER_CORE,3) * &
!      R_memory(5,ijk,1,1,ispec,3) + muvstore(ijk,1,1,ispec) * &
!      factor_common(IREGION_ATTENUATION_INNER_CORE,3) * &
!      (betaval(IREGION_ATTENUATION_INNER_CORE,3) * &
!      epsilondev(5,ijk,1,1,ispec) + gammaval(IREGION_ATTENUATION_INNER_CORE,3) * epsilondev_loc(5,ijk,1,1))

       if(ATTENUATION_VAL_3D) then
          factor_common_use(1) = factor_common(ijk,1,1,ispec,1)
          factor_common_use(2) = factor_common(ijk,1,1,ispec,2)
          factor_common_use(3) = factor_common(ijk,1,1,ispec,3)
       else
          factor_common_use(1) = factor_common(1,1,1,IREGION_ATTENUATION_INNER_CORE,1)
          factor_common_use(2) = factor_common(1,1,1,IREGION_ATTENUATION_INNER_CORE,2)
          factor_common_use(3) = factor_common(1,1,1,IREGION_ATTENUATION_INNER_CORE,3)
       endif
    R_memory(1,ijk,1,1,ispec,1) = &
      alphaval(1) * &
      R_memory(1,ijk,1,1,ispec,1) + muvstore(ijk,1,1,ispec) * &
      factor_common_use(1) * &
      (betaval(1) * &
      epsilondev(1,ijk,1,1,ispec) + gammaval(1) * epsilondev_loc(1,ijk,1,1))

    R_memory(2,ijk,1,1,ispec,1) = &
      alphaval(1) * &
      R_memory(2,ijk,1,1,ispec,1) + muvstore(ijk,1,1,ispec) * &
      factor_common_use(1) * &
      (betaval(1) * &
      epsilondev(2,ijk,1,1,ispec) + gammaval(1) * epsilondev_loc(2,ijk,1,1))

    R_memory(3,ijk,1,1,ispec,1) = &
      alphaval(1) * &
      R_memory(3,ijk,1,1,ispec,1) + muvstore(ijk,1,1,ispec) * &
      factor_common_use(1) * &
      (betaval(1) * &
      epsilondev(3,ijk,1,1,ispec) + gammaval(1) * epsilondev_loc(3,ijk,1,1))

    R_memory(4,ijk,1,1,ispec,1) = &
      alphaval(1) * &
      R_memory(4,ijk,1,1,ispec,1) + muvstore(ijk,1,1,ispec) * &
      factor_common_use(1) * &
      (betaval(1) * &
      epsilondev(4,ijk,1,1,ispec) + gammaval(1) * epsilondev_loc(4,ijk,1,1))

    R_memory(5,ijk,1,1,ispec,1) = &
      alphaval(1) * &
      R_memory(5,ijk,1,1,ispec,1) + muvstore(ijk,1,1,ispec) * &
      factor_common_use(1) * &
      (betaval(1) * &
      epsilondev(5,ijk,1,1,ispec) + gammaval(1) * epsilondev_loc(5,ijk,1,1))

!---
    R_memory(1,ijk,1,1,ispec,2) = &
      alphaval(2) * &
      R_memory(1,ijk,1,1,ispec,2) + muvstore(ijk,1,1,ispec) * &
      factor_common_use(2) * &
      (betaval(2) * &
      epsilondev(1,ijk,1,1,ispec) + gammaval(2) * epsilondev_loc(1,ijk,1,1))

    R_memory(2,ijk,1,1,ispec,2) = &
      alphaval(2) * &
      R_memory(2,ijk,1,1,ispec,2) + muvstore(ijk,1,1,ispec) * &
      factor_common_use(2) * &
      (betaval(2) * &
      epsilondev(2,ijk,1,1,ispec) + gammaval(2) * epsilondev_loc(2,ijk,1,1))

    R_memory(3,ijk,1,1,ispec,2) = &
      alphaval(2) * &
      R_memory(3,ijk,1,1,ispec,2) + muvstore(ijk,1,1,ispec) * &
      factor_common_use(2) * &
      (betaval(2) * &
      epsilondev(3,ijk,1,1,ispec) + gammaval(2) * epsilondev_loc(3,ijk,1,1))

    R_memory(4,ijk,1,1,ispec,2) = &
      alphaval(2) * &
      R_memory(4,ijk,1,1,ispec,2) + muvstore(ijk,1,1,ispec) * &
      factor_common_use(2) * &
      (betaval(2) * &
      epsilondev(4,ijk,1,1,ispec) + gammaval(2) * epsilondev_loc(4,ijk,1,1))

    R_memory(5,ijk,1,1,ispec,2) = &
      alphaval(2) * &
      R_memory(5,ijk,1,1,ispec,2) + muvstore(ijk,1,1,ispec) * &
      factor_common_use(2) * &
      (betaval(2) * &
      epsilondev(5,ijk,1,1,ispec) + gammaval(2) * epsilondev_loc(5,ijk,1,1))

!---
    R_memory(1,ijk,1,1,ispec,3) = &
      alphaval(3) * &
      R_memory(1,ijk,1,1,ispec,3) + muvstore(ijk,1,1,ispec) * &
      factor_common_use(3) * &
      (betaval(3) * &
      epsilondev(1,ijk,1,1,ispec) + gammaval(3) * epsilondev_loc(1,ijk,1,1))

    R_memory(2,ijk,1,1,ispec,3) = &
      alphaval(3) * &
      R_memory(2,ijk,1,1,ispec,3) + muvstore(ijk,1,1,ispec) * &
      factor_common_use(3) * &
      (betaval(3) * &
      epsilondev(2,ijk,1,1,ispec) + gammaval(3) * epsilondev_loc(2,ijk,1,1))

    R_memory(3,ijk,1,1,ispec,3) = &
      alphaval(3) * &
      R_memory(3,ijk,1,1,ispec,3) + muvstore(ijk,1,1,ispec) * &
      factor_common_use(3) * &
      (betaval(3) * &
      epsilondev(3,ijk,1,1,ispec) + gammaval(3) * epsilondev_loc(3,ijk,1,1))

    R_memory(4,ijk,1,1,ispec,3) = &
      alphaval(3) * &
      R_memory(4,ijk,1,1,ispec,3) + muvstore(ijk,1,1,ispec) * &
      factor_common_use(3) * &
      (betaval(3) * &
      epsilondev(4,ijk,1,1,ispec) + gammaval(3) * epsilondev_loc(4,ijk,1,1))

    R_memory(5,ijk,1,1,ispec,3) = &
      alphaval(3) * &
      R_memory(5,ijk,1,1,ispec,3) + muvstore(ijk,1,1,ispec) * &
      factor_common_use(3) * &
      (betaval(3) * &
      epsilondev(5,ijk,1,1,ispec) + gammaval(3) * epsilondev_loc(5,ijk,1,1))

! BS END
    enddo

! save deviatoric strain for Runge-Kutta scheme
    do ijk = 1,5*NGLLCUBE
      epsilondev(ijk,1,1,1,ispec) = epsilondev_loc(ijk,1,1,1)
    enddo

  endif

  endif   ! end test to exclude fictitious elements in central cube

  enddo ! spectral element loop

  end subroutine compute_forces_inner_core

