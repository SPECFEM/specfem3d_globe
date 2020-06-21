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
! the Free Software Foundation; either version 3 of the License, or
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

  subroutine hex_nodes(iaddx,iaddy,iaddz)

  use constants, only: NGNOD

  implicit none

! topology of the elements
  integer, dimension(NGNOD), intent(out) :: iaddx,iaddy,iaddz

! define the topology of the hexahedral elements

! the topology of the nodes is described in UTILS/chunk_notes_scanned/numbering_convention_27_nodes.tif

  if (NGNOD /= 27) stop 'elements should have 27 control nodes'

! corner nodes

  iaddx(1) = 0
  iaddy(1) = 0
  iaddz(1) = 0

  iaddx(2) = 2
  iaddy(2) = 0
  iaddz(2) = 0

  iaddx(3) = 2
  iaddy(3) = 2
  iaddz(3) = 0

  iaddx(4) = 0
  iaddy(4) = 2
  iaddz(4) = 0

  iaddx(5) = 0
  iaddy(5) = 0
  iaddz(5) = 2

  iaddx(6) = 2
  iaddy(6) = 0
  iaddz(6) = 2

  iaddx(7) = 2
  iaddy(7) = 2
  iaddz(7) = 2

  iaddx(8) = 0
  iaddy(8) = 2
  iaddz(8) = 2

! midside nodes (nodes located in the middle of an edge)

  iaddx(9) = 1
  iaddy(9) = 0
  iaddz(9) = 0

  iaddx(10) = 2
  iaddy(10) = 1
  iaddz(10) = 0

  iaddx(11) = 1
  iaddy(11) = 2
  iaddz(11) = 0

  iaddx(12) = 0
  iaddy(12) = 1
  iaddz(12) = 0

  iaddx(13) = 0
  iaddy(13) = 0
  iaddz(13) = 1

  iaddx(14) = 2
  iaddy(14) = 0
  iaddz(14) = 1

  iaddx(15) = 2
  iaddy(15) = 2
  iaddz(15) = 1

  iaddx(16) = 0
  iaddy(16) = 2
  iaddz(16) = 1

  iaddx(17) = 1
  iaddy(17) = 0
  iaddz(17) = 2

  iaddx(18) = 2
  iaddy(18) = 1
  iaddz(18) = 2

  iaddx(19) = 1
  iaddy(19) = 2
  iaddz(19) = 2

  iaddx(20) = 0
  iaddy(20) = 1
  iaddz(20) = 2

! side center nodes (nodes located in the middle of a face)

  iaddx(21) = 1
  iaddy(21) = 1
  iaddz(21) = 0

  iaddx(22) = 1
  iaddy(22) = 0
  iaddz(22) = 1

  iaddx(23) = 2
  iaddy(23) = 1
  iaddz(23) = 1

  iaddx(24) = 1
  iaddy(24) = 2
  iaddz(24) = 1

  iaddx(25) = 0
  iaddy(25) = 1
  iaddz(25) = 1

  iaddx(26) = 1
  iaddy(26) = 1
  iaddz(26) = 2

! center node (barycenter of the eight corners)

  iaddx(27) = 1
  iaddy(27) = 1
  iaddz(27) = 1

  end subroutine hex_nodes

!
!-------------------------------------------------------------------------------------------------
!


  subroutine hex_nodes_anchor_ijk(anchor_iax,anchor_iay,anchor_iaz)

! gets control point indices
!
! to get coordinates of control points (corners,midpoints) for an element ispec, they can be use as:
!do ia = 1,NGNOD
!  iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
!  xelm(ia) = dble(xstore(iglob))
!  yelm(ia) = dble(ystore(iglob))
!  zelm(ia) = dble(zstore(iglob))
!enddo

  use constants, only: NGNOD,MIDX,MIDY,MIDZ,NGLLX,NGLLY,NGLLZ

  implicit none

  integer, dimension(NGNOD), intent(out) :: anchor_iax,anchor_iay,anchor_iaz

  ! local parameters
  integer :: ia
  integer :: iax,iay,iaz

  ! topology of the control/anchor points of the surface element
  integer :: iaddx(NGNOD),iaddy(NGNOD),iaddr(NGNOD)

  ! define topology of the control element
  call hex_nodes(iaddx,iaddy,iaddr)

  ! define (i,j,k) indices of the control/anchor points of the elements
  do ia = 1,NGNOD
    ! control point index
    iax = 0
    if (iaddx(ia) == 0) then
      iax = 1
    else if (iaddx(ia) == 1) then
      iax = MIDX
    else if (iaddx(ia) == 2) then
      iax = NGLLX
    else
      stop 'incorrect value of iaddx'
    endif
    anchor_iax(ia) = iax

    iay = 0
    if (iaddy(ia) == 0) then
      iay = 1
    else if (iaddy(ia) == 1) then
      iay = MIDY
    else if (iaddy(ia) == 2) then
      iay = NGLLY
    else
      stop 'incorrect value of iaddy'
    endif
    anchor_iay(ia) = iay

    iaz = 0
    if (iaddr(ia) == 0) then
      iaz = 1
    else if (iaddr(ia) == 1) then
      iaz = MIDZ
    else if (iaddr(ia) == 2) then
      iaz = NGLLZ
    else
      stop 'incorrect value of iaddr'
    endif
    anchor_iaz(ia) = iaz
  enddo

  end subroutine hex_nodes_anchor_ijk

