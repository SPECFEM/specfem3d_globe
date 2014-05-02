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


  subroutine heap_sort( N, array )

! heap sort algorithm
! sorts integer array (in increasing order, like 1 - 5 - 6 - 9 - 12 - 13 - 14 -...)

  implicit none
  integer,intent(in) :: N
  integer,dimension(N),intent(inout) :: array

  ! local parameters
  integer :: tmp
  integer :: i

  ! checks if anything to do
  if( N < 2 ) return

  ! builds heap
  do i = N/2, 1, -1
    call heap_sort_siftdown(N,array,i,N)
  enddo

  ! sorts array
  do i = N, 2, -1
    ! swaps last and first entry in this section
    tmp = array(1)
    array(1) = array(i)
    array(i) = tmp
    call heap_sort_siftdown(N,array,1,i-1)
  enddo

  end subroutine heap_sort

!
!----
!

  subroutine heap_sort_siftdown(N,array,start,bottom)

  implicit none

  integer,intent(in):: N
  integer,dimension(N),intent(inout) :: array
  integer :: start,bottom

  ! local parameters
  integer :: i,j
  integer :: tmp

  i = start
  tmp = array(i)
  j = 2*i
  do while( j <= bottom )
    ! chooses larger value first in this section
    if( j < bottom ) then
      if( array(j) <= array(j+1) ) j = j + 1
    endif

    ! checks if section already smaller than initial value
    if( array(j) < tmp ) exit

    array(i) = array(j)
    i = j
    j = 2*i
  enddo

  array(i) = tmp
  return

  end subroutine heap_sort_siftdown

