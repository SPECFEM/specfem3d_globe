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

  subroutine reduce(theta,phi)

! bring theta between 0 and PI, and phi between 0 and 2*PI

  use constants, only: ZERO,PI,TWO_PI,TINYVAL

  implicit none

  double precision,intent(inout) :: theta,phi

  ! local parameters
  integer :: i
  double precision :: th,ph

  ! slightly move points to avoid roundoff problem when exactly on the polar axis
  if (abs(theta) < TINYVAL ) theta = theta + 0.0000001d0
  if (abs(phi) < TINYVAL ) phi = phi + 0.0000001d0

  ! colatitude
  th = theta
  ! longitude
  ph = phi

  ! brings longitude between 0 and 2*PI
  if (ph < ZERO .or. ph > TWO_PI) then
    i = abs(int(ph/TWO_PI))
    if (ph < ZERO) then
      ph = ph+(i+1)*TWO_PI
    else
      if (ph > TWO_PI) ph=ph-i*TWO_PI
    endif
    phi=ph
  endif

  ! brings colatitude between 0 and PI
  if (th < ZERO .or. th > PI) then
    i=int(th/PI)
    if (th > ZERO) then
      if (mod(i,2) /= 0) then
        th=(i+1)*PI-th
        ! switches hemisphere
        if (ph < PI) then
          ph=ph+PI
        else
          ph=ph-PI
        endif
      else
        th=th-i*PI
      endif
    else
      if (mod(i,2) == 0) then
        th=-th+i*PI
        ! switches hemisphere
        if (ph < PI) then
          ph=ph+PI
        else
          ph=ph-PI
        endif
      else
        th=th-i*PI
      endif
    endif
    theta=th
    phi=ph
  endif

  ! checks ranges
  if (theta < ZERO .or. theta > PI) stop 'theta out of range in reduce'
  if (phi < ZERO .or. phi > TWO_PI) stop 'phi out of range in reduce'

  end subroutine reduce

