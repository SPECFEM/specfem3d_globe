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

! subroutines to sort MPI buffers to assemble between chunks

  subroutine sort_array_coordinates(npointot,x,y,z, &
                                    ibool,iglob,locval,ifseg,nglob, &
                                    ind,ninseg,iwork,work)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

  use constants

  implicit none

  integer :: npointot,nglob

  double precision,dimension(npointot) :: x,y,z

  integer,dimension(npointot) :: ibool,iglob,locval
  integer,dimension(npointot) :: ind,ninseg
  logical,dimension(npointot) :: ifseg

  integer,dimension(npointot) :: iwork
  double precision,dimension(npointot) :: work

  ! local parameters
  integer :: ipoin,i,j
  integer :: nseg,ioff,iseg,ig

  ! establish initial pointers
  do ipoin=1,npointot
    locval(ipoin)=ipoin
  enddo

  ifseg(:)=.false.

  nseg=1
  ifseg(1)=.true.
  ninseg(1)=npointot

  do j=1,NDIM

    ! sort within each segment
    ioff=1
    do iseg=1,nseg
      if(j == 1) then
        call rank_buffers(x(ioff),ind,ninseg(iseg))
      else if(j == 2) then
        call rank_buffers(y(ioff),ind,ninseg(iseg))
      else
        call rank_buffers(z(ioff),ind,ninseg(iseg))
      endif

      call swap_all_buffers(ibool(ioff),locval(ioff), &
              x(ioff),y(ioff),z(ioff),iwork,work,ind,ninseg(iseg))

      ioff=ioff+ninseg(iseg)
    enddo

    ! check for jumps in current coordinate
    ! define a tolerance, normalized radius is 1., so let's use a small value
    if(j == 1) then
      do i=2,npointot
        if(dabs(x(i)-x(i-1)) > SMALLVALTOL ) ifseg(i)=.true.
      enddo
    else if(j == 2) then
      do i=2,npointot
        if(dabs(y(i)-y(i-1)) > SMALLVALTOL ) ifseg(i)=.true.
      enddo
    else
      do i=2,npointot
        if(dabs(z(i)-z(i-1)) > SMALLVALTOL ) ifseg(i)=.true.
      enddo
    endif

    ! count up number of different segments
    nseg=0
    do i=1,npointot
      if(ifseg(i)) then
        nseg=nseg+1
        ninseg(nseg)=1
      else
        ninseg(nseg)=ninseg(nseg)+1
      endif
    enddo

  enddo

  ! assign global node numbers (now sorted lexicographically)
  ig=0
  do i=1,npointot
    if(ifseg(i)) ig=ig+1
    iglob(locval(i))=ig
  enddo

  nglob=ig

  end subroutine sort_array_coordinates

! -------------------- library for sorting routine ------------------

! sorting routines put here in same file to allow for inlining

  subroutine rank_buffers(A,IND,N)
!
! Use Heap Sort (Numerical Recipes)
!
  implicit none

  integer :: n
  double precision,dimension(n) :: A
  integer,dimension(n) :: IND

  ! local parameters
  integer :: i,j,l,ir,indx
  double precision :: q

  do j=1,n
    IND(j)=j
  enddo

  if(n == 1) return

  L = n/2 + 1
  ir = n

  do while( .true. )

    IF ( l > 1 ) THEN
      l = l-1
      indx = ind(l)
      q = a(indx)
    ELSE
      indx = ind(ir)
      q = a(indx)
      ind(ir) = ind(1)
      ir = ir-1

      ! checks exit criterion
      if (ir == 1) then
         ind(1) = indx
         return
      endif

    endif

    i = l
    j = l+l

    do while( J <= IR )
      IF ( J < IR ) THEN
        IF ( A(IND(j)) < A(IND(j+1)) ) j=j+1
      endif
      IF ( q < A(IND(j)) ) THEN
        IND(I) = IND(J)
        I = J
        J = J+J
      ELSE
        J = IR+1
      endif

    enddo

    IND(I)=INDX
  enddo

  end subroutine rank_buffers

! -------------------------------------------------------------------

  subroutine swap_all_buffers(IA,IB,A,B,C,IW,W,ind,n)
!
! swap arrays IA, IB, A, B and C according to addressing in array IND
!
  implicit none

  integer :: n
  integer,dimension(n) :: IND
  integer,dimension(n) :: IA,IB,IW
  double precision,dimension(n) :: A,B,C,W

  ! local parameter
  integer :: i

  do i=1,n
    W(i)=A(i)
    IW(i)=IA(i)
  enddo

  do i=1,n
    A(i)=W(ind(i))
    IA(i)=IW(ind(i))
  enddo

  do i=1,n
    W(i)=B(i)
    IW(i)=IB(i)
  enddo

  do i=1,n
    B(i)=W(ind(i))
    IB(i)=IW(ind(i))
  enddo

  do i=1,n
    W(i)=C(i)
  enddo

  do i=1,n
    C(i)=W(ind(i))
  enddo

  end subroutine swap_all_buffers


