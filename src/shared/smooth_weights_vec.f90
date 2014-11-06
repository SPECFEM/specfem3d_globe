!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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



  subroutine smoothing_weights_vec(x0,y0,z0,ispec2,sigma_h2,sigma_v2,exp_val,&
                              xx_elem,yy_elem,zz_elem)


  use constants,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: exp_val
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: xx_elem, yy_elem, zz_elem
  real(kind=CUSTOM_REAL),intent(in) :: x0,y0,z0,sigma_h2,sigma_v2
  integer,intent(in) :: ispec2

  ! local parameters
  integer :: ii,jj,kk
  real(kind=CUSTOM_REAL) :: dist_h,dist_v

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
  ! to avoid compiler warning
  ii = ispec2

  do kk = 1, NGLLZ
    do jj = 1, NGLLY
      do ii = 1, NGLLX
        ! point in second slice

        ! vector approximation:
        call get_distance_vec(dist_h,dist_v,x0,y0,z0, &
            xx_elem(ii,jj,kk),yy_elem(ii,jj,kk),zz_elem(ii,jj,kk))

        ! Gaussian function
        exp_val(ii,jj,kk) = exp( - (dist_h*dist_h)/sigma_h2 &
                                  - (dist_v*dist_v)/sigma_v2 )    ! * factor(ii,jj,kk)
      enddo
    enddo
  enddo

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
  dist_v = r1 - r0
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
