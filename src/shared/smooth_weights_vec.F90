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


  subroutine smoothing_weights_vec(x0,y0,z0,sigma_h2,sigma_v2,exp_val,&
                                   xx_elem,yy_elem,zz_elem)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGLLCUBE

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: x0,y0,z0,sigma_h2,sigma_v2
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: exp_val
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: xx_elem, yy_elem, zz_elem

  ! local parameters
  real(kind=CUSTOM_REAL) :: dist_h,dist_v
  real(kind=CUSTOM_REAL) :: r0,r1
  real(kind=CUSTOM_REAL) :: theta,ratio
  real(kind=CUSTOM_REAL) :: x1,y1,z1
  real(kind=CUSTOM_REAL) :: sigma_h2_inv,sigma_v2_inv
  real(kind=CUSTOM_REAL) :: val

#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
#else
  integer :: i,j,k
#endif

  ! >>>>>
  ! uniform sigma
  !exp_val(:,:,:) = exp( -((xx(:,:,:,ispec2)-x0)**2+(yy(:,:,:,ispec2)-y0)**2 &
  !          +(zz(:,:,:,ispec2)-z0)**2 )/(2*sigma2) )*factor(:,:,:)

  ! from basin code smoothing:
  ! gaussian function
  !exp_val(:,:,:) = exp( -(xx(:,:,:,ispec2)-x0)**2/(sigma_h2) &
  !                      -(yy(:,:,:,ispec2)-y0)**2/(sigma_h2) &
  !                      -(zz(:,:,:,ispec2)-z0)**2/(sigma_v2) ) * factor(:,:,:)
  ! >>>>>

  ! helper variables
  sigma_h2_inv = 1.0_CUSTOM_REAL / sigma_h2
  sigma_v2_inv = 1.0_CUSTOM_REAL / sigma_v2

  ! length of first position vector
  r0 = sqrt( x0*x0 + y0*y0 + z0*z0 )

  DO_LOOP_IJK

    ! point in second slice
    x1 = xx_elem(INDEX_IJK)
    y1 = yy_elem(INDEX_IJK)
    z1 = zz_elem(INDEX_IJK)

    ! without explicit function calls to help compiler optimize loops

    ! length of position vector
    r1 = sqrt( x1*x1 + y1*y1 + z1*z1 )

    ! vertical distance (squared)
    dist_v = (r1 - r0)*(r1 - r0)

    ! only for flat earth with z in depth: dist_v = sqrt( (cz(ispec2)-cz0(ispec))** 2)

    ! epicentral distance
    ! (accounting for spherical curvature)
    ! calculates distance of circular segment
    ! angle between r0 and r1 in radian
    ! given by dot-product of two vectors
    ratio = (x0*x1 + y0*y1 + z0*z1) / (r0 * r1)

    ! checks boundaries of ratio (due to numerical inaccuracies)
    if (ratio > 1.0_CUSTOM_REAL) then
      ratio = 1.0_CUSTOM_REAL
    else if (ratio < -1.0_CUSTOM_REAL) then
      ratio = -1.0_CUSTOM_REAL
    endif

    theta = acos( ratio )

    ! segment length at heigth of r1 (squared)
    dist_h = (r1 * theta)*(r1 * theta)

    ! Gaussian function
    val = - dist_h*sigma_h2_inv - dist_v*sigma_v2_inv

    !exp_val(INDEX_IJK) = exp(val) ! * factor(INDEX_IJK)

    ! limits to single precision
    if (val < - 86.0) then
      ! smaller than numerical precision: exp(-86) < 1.e-37
      exp_val(INDEX_IJK) = 0.0_CUSTOM_REAL
    else
      exp_val(INDEX_IJK) = exp(val)    ! * factor(INDEX_IJK)
    endif

    ! debug
    !if (debug) then
    !  print*,INDEX_IJK,'smoothing:',dist_v,dist_h,sigma_h2,sigma_v2,ratio,theta,'val',- dist_h/sigma_h2 - dist_v/sigma_v2
    !endif

  ENDDO_LOOP_IJK

  end subroutine smoothing_weights_vec


!
! -----------------------------------------------------------------------------
!

  subroutine get_distance_vec(dist_h,dist_v,x0,y0,z0,x1,y1,z1)

! returns vector lengths as distances in radial and horizontal direction

  use constants,only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL),intent(out) :: dist_h,dist_v
  real(kind=CUSTOM_REAL),intent(in) :: x0,y0,z0,x1,y1,z1

  ! local parameters
  real(kind=CUSTOM_REAL) :: r0,r1
  real(kind=CUSTOM_REAL) :: theta,ratio
  !real(kind=CUSTOM_REAL) :: vx,vy,vz,alpha

  ! vertical distance
  r0 = sqrt( x0*x0 + y0*y0 + z0*z0 ) ! length of first position vector
  r1 = sqrt( x1*x1 + y1*y1 + z1*z1 )
  dist_v = abs(r1 - r0)

  ! only for flat earth with z in depth: dist_v = sqrt( (cz(ispec2)-cz0(ispec))** 2)

  ! epicentral distance
  ! (accounting for spherical curvature)
  ! calculates distance of circular segment
  ! angle between r0 and r1 in radian
  ! given by dot-product of two vectors
  ratio = (x0*x1 + y0*y1 + z0*z1)/(r0 * r1)

  ! checks boundaries of ratio (due to numerical inaccuracies)
  if (ratio > 1.0_CUSTOM_REAL) ratio = 1.0_CUSTOM_REAL
  if (ratio < -1.0_CUSTOM_REAL) ratio = -1.0_CUSTOM_REAL

  theta = acos( ratio )

  ! segment length at heigth of r1
  dist_h = r1 * theta

  ! vector approximation (fast computation): neglects curvature
  ! horizontal distance
  ! length of vector from point 0 to point 1
  ! assuming small earth curvature  (since only for neighboring elements)

  ! scales r0 to have same length as r1
  !alpha = r1 / r0
  !vx = alpha * x0
  !vy = alpha * y0
  !vz = alpha * z0

  ! vector in horizontal between new r0 and r1
  !vx = x1 - vx
  !vy = y1 - vy
  !vz = z1 - vz

  ! distance is vector length
  !dist_h = sqrt( vx*vx + vy*vy + vz*vz )

  end subroutine get_distance_vec

