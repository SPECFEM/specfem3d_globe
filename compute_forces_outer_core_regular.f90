!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology September 2002
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
          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
          ibool,nspec_outer_core,nglob_outer_core)

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
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempx2,tempx3

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
  integer i,j,k,l

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l,sum_terms

! ****************************************************
!   big loop over all spectral elements in the fluid
! ****************************************************

! set acceleration to zero
  accelfluid(:) = 0._CUSTOM_REAL

  do ispec = 1,nspec_outer_core

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          do l=1,NGLLX
            tempx1l = tempx1l + displfluid(ibool(l,j,k,ispec)) * hprime_xx(l,i)
          enddo

          do l=1,NGLLY
            tempx2l = tempx2l + displfluid(ibool(i,l,k,ispec)) * hprime_yy(l,j)
          enddo

          do l=1,NGLLZ
            tempx3l = tempx3l + displfluid(ibool(i,j,l,ispec)) * hprime_zz(l,k)
          enddo

!         get derivatives of velocity potential with respect to x, y and z

          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)
          jacobianl = jacobian(i,j,k,ispec)

          dpotentialdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          dpotentialdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          dpotentialdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

! compute contribution of rotation and add to gradient of potential
! this term has no Z component
    if(ROTATION_VAL) then

! store the source for the Euler scheme for A_rotation and B_rotation
      two_omega_deltat = deltat * two_omega_earth

      cos_two_omega_t = cos(two_omega_earth*time)
      sin_two_omega_t = sin(two_omega_earth*time)

! time step deltat of Euler scheme is included in the source
      source_euler_A(i,j,k) = two_omega_deltat * (cos_two_omega_t * dpotentialdyl + sin_two_omega_t * dpotentialdxl)
      source_euler_B(i,j,k) = two_omega_deltat * (sin_two_omega_t * dpotentialdyl - cos_two_omega_t * dpotentialdxl)

      A_rotation = A_array_rotation(i,j,k,ispec)
      B_rotation = B_array_rotation(i,j,k,ispec)

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

            iglob = ibool(i,j,k,ispec)
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
              gravity_term(i,j,k) = &
                sngl(minus_rho_g_over_kappa_fluid(int_radius) * &
                dble(jacobianl) * wgll_cube(i,j,k) * &
               (dble(dpotentialdx_with_rot) * gxl + &
                dble(dpotentialdy_with_rot) * gyl + dble(dpotentialdzl) * gzl))
            else
              gravity_term(i,j,k) = minus_rho_g_over_kappa_fluid(int_radius) * &
                 jacobianl * wgll_cube(i,j,k) * (dpotentialdx_with_rot * gxl + &
                 dpotentialdy_with_rot * gyl + dpotentialdzl * gzl)
            endif

          endif

          tempx1(i,j,k) = jacobianl*(xixl*dpotentialdx_with_rot + xiyl*dpotentialdy_with_rot + xizl*dpotentialdzl)
          tempx2(i,j,k) = jacobianl*(etaxl*dpotentialdx_with_rot + etayl*dpotentialdy_with_rot + etazl*dpotentialdzl)
          tempx3(i,j,k) = jacobianl*(gammaxl*dpotentialdx_with_rot + gammayl*dpotentialdy_with_rot + gammazl*dpotentialdzl)

          enddo
        enddo
      enddo

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          do l=1,NGLLX
            tempx1l = tempx1l + tempx1(l,j,k) * hprimewgll_xx(i,l)
          enddo

          do l=1,NGLLY
            tempx2l = tempx2l + tempx2(i,l,k) * hprimewgll_yy(j,l)
          enddo

          do l=1,NGLLZ
            tempx3l = tempx3l + tempx3(i,j,l) * hprimewgll_zz(k,l)
          enddo

! sum contributions from each element to the global mesh and add gravity term

          iglob = ibool(i,j,k,ispec)
          sum_terms = - (wgllwgll_yz(j,k)*tempx1l + wgllwgll_xz(i,k)*tempx2l + wgllwgll_xy(i,j)*tempx3l)
          if(GRAVITY_VAL) sum_terms = sum_terms + gravity_term(i,j,k)
          accelfluid(iglob) = accelfluid(iglob) + sum_terms

        enddo
      enddo
    enddo

! update rotation term with Euler scheme
    if(ROTATION_VAL) then

! use the source saved above
      A_array_rotation(:,:,:,ispec) = A_array_rotation(:,:,:,ispec) + source_euler_A(:,:,:)
      B_array_rotation(:,:,:,ispec) = B_array_rotation(:,:,:,ispec) + source_euler_B(:,:,:)

    endif

  enddo   ! spectral element loop

  end subroutine compute_forces_outer_core

