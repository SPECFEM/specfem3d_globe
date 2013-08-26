!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
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

  program convert_DX_to_AVS_slices

  implicit none

  integer, parameter :: npoin = 1722, nelem = 1612

  integer i1,i2,i3,i4,ipoin,ielem

  real x,y,z,datavalue

! skip DX header
  read(*,*)

! write AVS header
  write(*,*) npoin,' ',nelem,' 0 1 0'

! read and write the points
  do ipoin = 1,npoin
    read(*,*) x,y,z
    write(*,*) ipoin,x,y,z
  enddo

! skip DX header
  read(*,*)

! read and write the elements
! in the case of OpenDX, node numbers start at zero
! in the case of AVS, node numbers start at one
! point order in OpenDX is 1,4,2,3 *not* 1,2,3,4 as in AVS
  do ielem = 1,nelem
    read(*,*) i1,i4,i2,i3
    write(*,*) ielem,' 1 quad ',i1+1,i2+1,i3+1,i4+1
  enddo

! skip DX header
  read(*,*)
  read(*,*)
  read(*,*)

! write AVS header
  write(*,*) '1 1'
  write(*,*) 'Zcoord, meters'

! read and write the element data
  do ielem = 1,nelem
    read(*,*) datavalue
    write(*,*) ielem,datavalue
  enddo

  end program convert_DX_to_AVS_slices

