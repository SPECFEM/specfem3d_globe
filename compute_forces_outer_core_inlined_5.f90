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

  subroutine compute_forces_outer_core(time,deltat,two_omega_earth, &
          A_array_rotation,B_array_rotation, &
          minus_rho_g_over_kappa_fluid,displfluid,accelfluid, &
          xstore,ystore,zstore, &
          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian, &
          hprime_xx,hprime_yy,hprime_zz, &
          hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
          wgll_cube,wgllwgll_yz_no_i,wgllwgll_xz_no_j,wgllwgll_xy_no_k, &
          ibool,nspec_outer_core,nglob_outer_core,index_fluid_i,index_fluid_k)

  implicit none

  include "constants.h"

! include values created by the mesher
! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

! for doubling in the outer core
  integer nspec_outer_core,nglob_outer_core

! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(nglob_outer_core) :: displfluid,accelfluid

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec_outer_core) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_outer_core) :: xix,xiy,xiz, &
                      etax,etay,etaz,gammax,gammay,gammaz,jacobian

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: wgllwgll_yz_no_i,wgllwgll_xz_no_j,wgllwgll_xy_no_k

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempx1loc,tempx2loc,tempx3loc, &
    tempx1inline,tempx2inline,tempx3inline,disploc,sum_terms

! for gravity
  integer int_radius
  double precision radius,theta,phi,gxl,gyl,gzl
  double precision cos_theta,sin_theta,cos_phi,sin_phi
  double precision, dimension(NRAD_GRAVITY) :: minus_rho_g_over_kappa_fluid
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: gravity_term
  real(kind=CUSTOM_REAL), dimension(nglob_outer_core) :: xstore,ystore,zstore

! for the Euler scheme for rotation
  real(kind=CUSTOM_REAL) time,deltat,two_omega_earth
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_OUTER_CORE_ROTATION) :: &
    A_array_rotation,B_array_rotation

  real(kind=CUSTOM_REAL) two_omega_deltat,cos_two_omega_t,sin_two_omega_t,A_rotation,B_rotation, &
       ux_rotation,uy_rotation,dpotentialdx_with_rot,dpotentialdy_with_rot
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: source_euler_A,source_euler_B

  integer ispec,iglob
  integer i,k,ij,ijk
  integer, dimension(NGLLSQUARE) :: index_fluid_i,index_fluid_k

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) dpotentialdxl,dpotentialdyl,dpotentialdzl

! ****************************************************
!   big loop over all spectral elements in the fluid
! ****************************************************

! set acceleration to zero
  accelfluid(:) = 0._CUSTOM_REAL

!CDIR NOVECTOR
  do ispec = 1,nspec_outer_core

! copy global displacement to local
    do ijk=1,NGLLCUBE
      disploc(ijk,1,1) = displfluid(ibool(ijk,1,1,ispec))
    enddo

! inlined first matrix product

!CDIR NODEP(tempx2loc)
  do ij = 1,NGLLSQUARE

    i = index_fluid_i(ij)
    k = index_fluid_k(ij)

!---
!--- ij is actually jk here
    tempx1loc(1,ij,1) = &
        disploc(1,ij,1)*hprime_xx(1,1) + disploc(2,ij,1)*hprime_xx(2,1) + &
        disploc(3,ij,1)*hprime_xx(3,1) + disploc(4,ij,1)*hprime_xx(4,1) + &
        disploc(5,ij,1)*hprime_xx(5,1)

!---
    tempx2loc(i,1,k) = &
        disploc(i,1,k)*hprime_yy(1,1) + disploc(i,2,k)*hprime_yy(2,1) + &
        disploc(i,3,k)*hprime_yy(3,1) + disploc(i,4,k)*hprime_yy(4,1) + &
        disploc(i,5,k)*hprime_yy(5,1)

!---
    tempx3loc(ij,1,1) = &
        disploc(ij,1,1)*hprime_zz(1,1) + disploc(ij,1,2)*hprime_zz(2,1) + &
        disploc(ij,1,3)*hprime_zz(3,1) + disploc(ij,1,4)*hprime_zz(4,1) + &
        disploc(ij,1,5)*hprime_zz(5,1)


!---
!--- ij is actually jk here
    tempx1loc(2,ij,1) = &
        disploc(1,ij,1)*hprime_xx(1,2) + disploc(2,ij,1)*hprime_xx(2,2) + &
        disploc(3,ij,1)*hprime_xx(3,2) + disploc(4,ij,1)*hprime_xx(4,2) + &
        disploc(5,ij,1)*hprime_xx(5,2)

!---
    tempx2loc(i,2,k) = &
        disploc(i,1,k)*hprime_yy(1,2) + disploc(i,2,k)*hprime_yy(2,2) + &
        disploc(i,3,k)*hprime_yy(3,2) + disploc(i,4,k)*hprime_yy(4,2) + &
        disploc(i,5,k)*hprime_yy(5,2)

!---
    tempx3loc(ij,1,2) = &
        disploc(ij,1,1)*hprime_zz(1,2) + disploc(ij,1,2)*hprime_zz(2,2) + &
        disploc(ij,1,3)*hprime_zz(3,2) + disploc(ij,1,4)*hprime_zz(4,2) + &
        disploc(ij,1,5)*hprime_zz(5,2)


!---
!--- ij is actually jk here
    tempx1loc(3,ij,1) = &
        disploc(1,ij,1)*hprime_xx(1,3) + disploc(2,ij,1)*hprime_xx(2,3) + &
        disploc(3,ij,1)*hprime_xx(3,3) + disploc(4,ij,1)*hprime_xx(4,3) + &
        disploc(5,ij,1)*hprime_xx(5,3)

!---
    tempx2loc(i,3,k) = &
        disploc(i,1,k)*hprime_yy(1,3) + disploc(i,2,k)*hprime_yy(2,3) + &
        disploc(i,3,k)*hprime_yy(3,3) + disploc(i,4,k)*hprime_yy(4,3) + &
        disploc(i,5,k)*hprime_yy(5,3)

!---
    tempx3loc(ij,1,3) = &
        disploc(ij,1,1)*hprime_zz(1,3) + disploc(ij,1,2)*hprime_zz(2,3) + &
        disploc(ij,1,3)*hprime_zz(3,3) + disploc(ij,1,4)*hprime_zz(4,3) + &
        disploc(ij,1,5)*hprime_zz(5,3)


!---
!--- ij is actually jk here
    tempx1loc(4,ij,1) = &
        disploc(1,ij,1)*hprime_xx(1,4) + disploc(2,ij,1)*hprime_xx(2,4) + &
        disploc(3,ij,1)*hprime_xx(3,4) + disploc(4,ij,1)*hprime_xx(4,4) + &
        disploc(5,ij,1)*hprime_xx(5,4)

!---
    tempx2loc(i,4,k) = &
        disploc(i,1,k)*hprime_yy(1,4) + disploc(i,2,k)*hprime_yy(2,4) + &
        disploc(i,3,k)*hprime_yy(3,4) + disploc(i,4,k)*hprime_yy(4,4) + &
        disploc(i,5,k)*hprime_yy(5,4)

!---
    tempx3loc(ij,1,4) = &
        disploc(ij,1,1)*hprime_zz(1,4) + disploc(ij,1,2)*hprime_zz(2,4) + &
        disploc(ij,1,3)*hprime_zz(3,4) + disploc(ij,1,4)*hprime_zz(4,4) + &
        disploc(ij,1,5)*hprime_zz(5,4)


!---
!--- ij is actually jk here
    tempx1loc(5,ij,1) = &
        disploc(1,ij,1)*hprime_xx(1,5) + disploc(2,ij,1)*hprime_xx(2,5) + &
        disploc(3,ij,1)*hprime_xx(3,5) + disploc(4,ij,1)*hprime_xx(4,5) + &
        disploc(5,ij,1)*hprime_xx(5,5)

!---
    tempx2loc(i,5,k) = &
        disploc(i,1,k)*hprime_yy(1,5) + disploc(i,2,k)*hprime_yy(2,5) + &
        disploc(i,3,k)*hprime_yy(3,5) + disploc(i,4,k)*hprime_yy(4,5) + &
        disploc(i,5,k)*hprime_yy(5,5)

!---
    tempx3loc(ij,1,5) = &
        disploc(ij,1,1)*hprime_zz(1,5) + disploc(ij,1,2)*hprime_zz(2,5) + &
        disploc(ij,1,3)*hprime_zz(3,5) + disploc(ij,1,4)*hprime_zz(4,5) + &
        disploc(ij,1,5)*hprime_zz(5,5)

  enddo

! fused the three loops for inlined version
    do ijk=1,NGLLCUBE

!         get derivatives of velocity potential with respect to x, y and z

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

          dpotentialdxl = xixl*tempx1loc(ijk,1,1) + etaxl*tempx2loc(ijk,1,1) + gammaxl*tempx3loc(ijk,1,1)
          dpotentialdyl = xiyl*tempx1loc(ijk,1,1) + etayl*tempx2loc(ijk,1,1) + gammayl*tempx3loc(ijk,1,1)
          dpotentialdzl = xizl*tempx1loc(ijk,1,1) + etazl*tempx2loc(ijk,1,1) + gammazl*tempx3loc(ijk,1,1)

! compute contribution of rotation and add to gradient of potential
! this term has no Z component
    if(ROTATION_VAL) then

! store the source for the Euler scheme for A_rotation and B_rotation
      two_omega_deltat = deltat * two_omega_earth

      cos_two_omega_t = cos(two_omega_earth*time)
      sin_two_omega_t = sin(two_omega_earth*time)

! time step deltat of Euler scheme is included in the source
      source_euler_A(ijk,1,1) = two_omega_deltat * (cos_two_omega_t * dpotentialdyl + sin_two_omega_t * dpotentialdxl)
      source_euler_B(ijk,1,1) = two_omega_deltat * (sin_two_omega_t * dpotentialdyl - cos_two_omega_t * dpotentialdxl)

      A_rotation = A_array_rotation(ijk,1,1,ispec)
      B_rotation = B_array_rotation(ijk,1,1,ispec)

      ux_rotation =   A_rotation*cos_two_omega_t + B_rotation*sin_two_omega_t
      uy_rotation = - A_rotation*sin_two_omega_t + B_rotation*cos_two_omega_t

      dpotentialdx_with_rot = dpotentialdxl + ux_rotation
      dpotentialdy_with_rot = dpotentialdyl + uy_rotation
    else

      dpotentialdx_with_rot = dpotentialdxl
      dpotentialdy_with_rot = dpotentialdyl

    endif  ! end of section with rotation

! precompute and store gravity term
          if(GRAVITY_VAL) then

! use mesh coordinates to get theta and phi
! x y z contain r theta phi

            iglob = ibool(ijk,1,1,ispec)
            radius = dble(xstore(iglob))
            theta = dble(ystore(iglob))
            phi = dble(zstore(iglob))

            cos_theta = dcos(theta)
            sin_theta = dsin(theta)
            cos_phi = dcos(phi)
            sin_phi = dsin(phi)

! get g, rho and dg/dr=dg
! spherical components of the gravitational acceleration
! for efficiency replace with lookup table every 100 m in radial direction
            int_radius = nint(radius * R_EARTH_KM * 10.d0)

! Cartesian components of the gravitational acceleration
! integrate and multiply by rho / Kappa
            gxl = sin_theta*cos_phi
            gyl = sin_theta*sin_phi
            gzl = cos_theta

! distinguish whether single or double precision for reals
            if(CUSTOM_REAL == SIZE_REAL) then
              gravity_term(ijk,1,1) = &
                sngl(minus_rho_g_over_kappa_fluid(int_radius) * &
                dble(jacobianl) * wgll_cube(ijk,1,1) * &
               (dble(dpotentialdx_with_rot) * gxl + &
                dble(dpotentialdy_with_rot) * gyl + dble(dpotentialdzl) * gzl))
            else
              gravity_term(ijk,1,1) = minus_rho_g_over_kappa_fluid(int_radius) * &
                 jacobianl * wgll_cube(ijk,1,1) * (dpotentialdx_with_rot * gxl + &
                 dpotentialdy_with_rot * gyl + dpotentialdzl * gzl)
            endif

          endif

          tempx1(ijk,1,1) = jacobianl*(xixl*dpotentialdx_with_rot + xiyl*dpotentialdy_with_rot + xizl*dpotentialdzl)
          tempx2(ijk,1,1) = jacobianl*(etaxl*dpotentialdx_with_rot + etayl*dpotentialdy_with_rot + etazl*dpotentialdzl)
          tempx3(ijk,1,1) = jacobianl*(gammaxl*dpotentialdx_with_rot + gammayl*dpotentialdy_with_rot + gammazl*dpotentialdzl)

      enddo

! inlined second matrix product

!CDIR NODEP(tempx2inline)
  do ij = 1,NGLLSQUARE

    i = index_fluid_i(ij)
    k = index_fluid_k(ij)

!---
!--- ij is actually jk below
  tempx1inline(1,ij,1) = &
    tempx1(1,ij,1)*hprimewgll_xx(1,1) + tempx1(2,ij,1)*hprimewgll_xx(1,2) + &
    tempx1(3,ij,1)*hprimewgll_xx(1,3) + tempx1(4,ij,1)*hprimewgll_xx(1,4) + &
    tempx1(5,ij,1)*hprimewgll_xx(1,5)

!---
  tempx2inline(i,1,k) = &
    tempx2(i,1,k)*hprimewgll_yy(1,1) + tempx2(i,2,k)*hprimewgll_yy(1,2) + &
    tempx2(i,3,k)*hprimewgll_yy(1,3) + tempx2(i,4,k)*hprimewgll_yy(1,4) + &
    tempx2(i,5,k)*hprimewgll_yy(1,5)

!---
  tempx3inline(ij,1,1) = &
    tempx3(ij,1,1)*hprimewgll_zz(1,1) + tempx3(ij,1,2)*hprimewgll_zz(1,2) + &
    tempx3(ij,1,3)*hprimewgll_zz(1,3) + tempx3(ij,1,4)*hprimewgll_zz(1,4) + &
    tempx3(ij,1,5)*hprimewgll_zz(1,5)


!---
!--- ij is actually jk below
  tempx1inline(2,ij,1) = &
    tempx1(1,ij,1)*hprimewgll_xx(2,1) + tempx1(2,ij,1)*hprimewgll_xx(2,2) + &
    tempx1(3,ij,1)*hprimewgll_xx(2,3) + tempx1(4,ij,1)*hprimewgll_xx(2,4) + &
    tempx1(5,ij,1)*hprimewgll_xx(2,5)

!---
  tempx2inline(i,2,k) = &
    tempx2(i,1,k)*hprimewgll_yy(2,1) + tempx2(i,2,k)*hprimewgll_yy(2,2) + &
    tempx2(i,3,k)*hprimewgll_yy(2,3) + tempx2(i,4,k)*hprimewgll_yy(2,4) + &
    tempx2(i,5,k)*hprimewgll_yy(2,5)

!---
  tempx3inline(ij,1,2) = &
    tempx3(ij,1,1)*hprimewgll_zz(2,1) + tempx3(ij,1,2)*hprimewgll_zz(2,2) + &
    tempx3(ij,1,3)*hprimewgll_zz(2,3) + tempx3(ij,1,4)*hprimewgll_zz(2,4) + &
    tempx3(ij,1,5)*hprimewgll_zz(2,5)


!---
!--- ij is actually jk below
  tempx1inline(3,ij,1) = &
    tempx1(1,ij,1)*hprimewgll_xx(3,1) + tempx1(2,ij,1)*hprimewgll_xx(3,2) + &
    tempx1(3,ij,1)*hprimewgll_xx(3,3) + tempx1(4,ij,1)*hprimewgll_xx(3,4) + &
    tempx1(5,ij,1)*hprimewgll_xx(3,5)

!---
  tempx2inline(i,3,k) = &
    tempx2(i,1,k)*hprimewgll_yy(3,1) + tempx2(i,2,k)*hprimewgll_yy(3,2) + &
    tempx2(i,3,k)*hprimewgll_yy(3,3) + tempx2(i,4,k)*hprimewgll_yy(3,4) + &
    tempx2(i,5,k)*hprimewgll_yy(3,5)

!---
  tempx3inline(ij,1,3) = &
    tempx3(ij,1,1)*hprimewgll_zz(3,1) + tempx3(ij,1,2)*hprimewgll_zz(3,2) + &
    tempx3(ij,1,3)*hprimewgll_zz(3,3) + tempx3(ij,1,4)*hprimewgll_zz(3,4) + &
    tempx3(ij,1,5)*hprimewgll_zz(3,5)


!---
!--- ij is actually jk below
  tempx1inline(4,ij,1) = &
    tempx1(1,ij,1)*hprimewgll_xx(4,1) + tempx1(2,ij,1)*hprimewgll_xx(4,2) + &
    tempx1(3,ij,1)*hprimewgll_xx(4,3) + tempx1(4,ij,1)*hprimewgll_xx(4,4) + &
    tempx1(5,ij,1)*hprimewgll_xx(4,5)

!---
  tempx2inline(i,4,k) = &
    tempx2(i,1,k)*hprimewgll_yy(4,1) + tempx2(i,2,k)*hprimewgll_yy(4,2) + &
    tempx2(i,3,k)*hprimewgll_yy(4,3) + tempx2(i,4,k)*hprimewgll_yy(4,4) + &
    tempx2(i,5,k)*hprimewgll_yy(4,5)

!---
  tempx3inline(ij,1,4) = &
    tempx3(ij,1,1)*hprimewgll_zz(4,1) + tempx3(ij,1,2)*hprimewgll_zz(4,2) + &
    tempx3(ij,1,3)*hprimewgll_zz(4,3) + tempx3(ij,1,4)*hprimewgll_zz(4,4) + &
    tempx3(ij,1,5)*hprimewgll_zz(4,5)


!---
!--- ij is actually jk below
  tempx1inline(5,ij,1) = &
    tempx1(1,ij,1)*hprimewgll_xx(5,1) + tempx1(2,ij,1)*hprimewgll_xx(5,2) + &
    tempx1(3,ij,1)*hprimewgll_xx(5,3) + tempx1(4,ij,1)*hprimewgll_xx(5,4) + &
    tempx1(5,ij,1)*hprimewgll_xx(5,5)

!---
  tempx2inline(i,5,k) = &
    tempx2(i,1,k)*hprimewgll_yy(5,1) + tempx2(i,2,k)*hprimewgll_yy(5,2) + &
    tempx2(i,3,k)*hprimewgll_yy(5,3) + tempx2(i,4,k)*hprimewgll_yy(5,4) + &
    tempx2(i,5,k)*hprimewgll_yy(5,5)

!---
  tempx3inline(ij,1,5) = &
    tempx3(ij,1,1)*hprimewgll_zz(5,1) + tempx3(ij,1,2)*hprimewgll_zz(5,2) + &
    tempx3(ij,1,3)*hprimewgll_zz(5,3) + tempx3(ij,1,4)*hprimewgll_zz(5,4) + &
    tempx3(ij,1,5)*hprimewgll_zz(5,5)

  enddo

  if(GRAVITY_VAL) then
    do ijk=1,NGLLCUBE
      sum_terms(ijk,1,1) = &
            - (wgllwgll_yz_no_i(1,ijk,1,1)*tempx1inline(ijk,1,1) + &
               wgllwgll_xz_no_j(1,ijk,1,1)*tempx2inline(ijk,1,1) + &
               wgllwgll_xy_no_k(1,ijk,1,1)*tempx3inline(ijk,1,1)) + &
               gravity_term(ijk,1,1)
    enddo
  else
    do ijk=1,NGLLCUBE
      sum_terms(ijk,1,1) = &
            - (wgllwgll_yz_no_i(1,ijk,1,1)*tempx1inline(ijk,1,1) + &
               wgllwgll_xz_no_j(1,ijk,1,1)*tempx2inline(ijk,1,1) + &
               wgllwgll_xy_no_k(1,ijk,1,1)*tempx3inline(ijk,1,1))
    enddo
  endif

! sum contributions to the global mesh
!CDIR NODEP(accelfluid)
    do ijk=1,NGLLCUBE
      iglob = ibool(ijk,1,1,ispec)
      accelfluid(iglob) = accelfluid(iglob) + sum_terms(ijk,1,1)
    enddo

! update rotation term with Euler scheme
    if(ROTATION_VAL) then

! use the source saved above
    do ijk=1,NGLLCUBE
      A_array_rotation(ijk,1,1,ispec) = A_array_rotation(ijk,1,1,ispec) + source_euler_A(ijk,1,1)
      B_array_rotation(ijk,1,1,ispec) = B_array_rotation(ijk,1,1,ispec) + source_euler_B(ijk,1,1)
    enddo

    endif

  enddo   ! spectral element loop

  end subroutine compute_forces_outer_core

