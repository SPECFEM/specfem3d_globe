!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            August 2008
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

  subroutine compute_forces_outer_core(d_ln_density_dr_table, &
          displfluid,accelfluid,xstore,ystore,zstore, &
          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
          hprime_xx,hprime_yy,hprime_zz, &
          hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,ibool)

  implicit none

  include "constants.h"

! include values created by the mesher
! done for performance only using static allocation to allow for loop unrolling
  include "values_from_mesher.h"

! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(nglob_outer_core) :: displfluid,accelfluid

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec_outer_core) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_outer_core) :: xix,xiy,xiz, &
                      etax,etay,etaz,gammax,gammay,gammaz

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempx2,tempx3

! for gravity
  integer int_radius
  double precision radius,theta,phi
  double precision cos_theta,sin_theta,cos_phi,sin_phi
  double precision, dimension(NRAD_GRAVITY) :: d_ln_density_dr_table
  real(kind=CUSTOM_REAL), dimension(nglob_outer_core) :: xstore,ystore,zstore

  integer ispec,iglob
  integer i,j,k,l

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l

  double precision grad_x_ln_rho,grad_y_ln_rho,grad_z_ln_rho
! ****************************************************
!   big loop over all spectral elements in the fluid
! ****************************************************

! set acceleration to zero
  accelfluid(:) = 0._CUSTOM_REAL

  do ispec = 1,NSPEC_OUTER_CORE

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          iglob = ibool(i,j,k,ispec)

          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          do l=1,NGLLX
            tempx1l = tempx1l + displfluid(ibool(l,j,k,ispec)) * hprime_xx(i,l)
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
            tempx2l = tempx2l + displfluid(ibool(i,l,k,ispec)) * hprime_yy(j,l)
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
            tempx3l = tempx3l + displfluid(ibool(i,j,l,ispec)) * hprime_zz(k,l)
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

! compute the jacobian
          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))

          dpotentialdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          dpotentialdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          dpotentialdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

! add (chi/rho)grad(rho) term in no gravity case

! With regards to the non-gravitating case: we cannot set N^2 = 0 *and* let g = 0.
! We can *either* assume N^2 = 0 but keep gravity g, *or* we can assume that gravity
! is negligible to begin with, as in our GJI 2002a, in which case N does not arise.
! We get:
!
! \ddot\chi = \rho^{-1}\kappa\bdel\cdot(\bdel\chi+\chi\bdel\ln\rho)
!
! Then the displacement is
!
! \bu = \bdel\chi+\chi\bdel\ln\rho = \rho^{-1}\bdel(\rho\chi)
!
! and the pressure is
!
! p = -\rho\ddot{\chi}
!
! Thus in our 2002b GJI paper eqn (21) is wrong, and equation (41)
! in our AGU monograph is incorrect; these equations should be replaced by
!
! \ddot\chi = \rho^{-1}\kappa\bdel\cdot(\bdel\chi+\chi\bdel\ln\rho)
!
! Note that the fluid potential we use in GJI 2002a differs from the one used here:
!
! \chi_GJI2002a = \rho\partial\t\chi
!
! such that
!
! \bv = \partial_t\bu=\rho^{-1}\bdel\chi_GJI2002a  (GJI 2002a eqn 20)
!
! p = - \partial_t\chi_GJI2002a (GJI 2002a eqn 19)

! use mesh coordinates to get theta and phi
! x y z contain r theta phi

      radius = dble(xstore(iglob))
      theta = dble(ystore(iglob))
      phi = dble(zstore(iglob))

      cos_theta = dcos(theta)
      sin_theta = dsin(theta)
      cos_phi = dcos(phi)
      sin_phi = dsin(phi)

      int_radius = nint(radius * R_EARTH_KM * 10.d0)

! grad(rho)/rho in Cartesian components
      grad_x_ln_rho = sin_theta * cos_phi * d_ln_density_dr_table(int_radius)
      grad_y_ln_rho = sin_theta * sin_phi * d_ln_density_dr_table(int_radius)
      grad_z_ln_rho = cos_theta * d_ln_density_dr_table(int_radius)

! adding (chi/rho)grad(rho)
      dpotentialdxl = dpotentialdxl + displfluid(iglob) * grad_x_ln_rho
      dpotentialdyl = dpotentialdyl + displfluid(iglob) * grad_y_ln_rho
      dpotentialdzl = dpotentialdzl + displfluid(iglob) * grad_z_ln_rho

          tempx1(i,j,k) = jacobianl*(xixl*dpotentialdxl + xiyl*dpotentialdyl + xizl*dpotentialdzl)
          tempx2(i,j,k) = jacobianl*(etaxl*dpotentialdxl + etayl*dpotentialdyl + etazl*dpotentialdzl)
          tempx3(i,j,k) = jacobianl*(gammaxl*dpotentialdxl + gammayl*dpotentialdyl + gammazl*dpotentialdzl)

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
            tempx1l = tempx1l + tempx1(l,j,k) * hprimewgll_xx(l,i)
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
            tempx2l = tempx2l + tempx2(i,l,k) * hprimewgll_yy(l,j)
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
            tempx3l = tempx3l + tempx3(i,j,l) * hprimewgll_zz(l,k)
          enddo

! sum contributions from each element to the global mesh and add gravity term

          iglob = ibool(i,j,k,ispec)
          accelfluid(iglob) = accelfluid(iglob) - (wgllwgll_yz(j,k)*tempx1l + wgllwgll_xz(i,k)*tempx2l + wgllwgll_xy(i,j)*tempx3l)

        enddo
      enddo
    enddo

  enddo   ! spectral element loop

  end subroutine compute_forces_outer_core

