!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 3
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

! subroutines to sort MPI buffers to assemble between chunks

  subroutine sort_array_coordinates(npointot,x,y,z,ibool,iglob,loc,ifseg,nglob,ind,ninseg,iwork,work)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

  implicit none

  include "constants.h"

  integer npointot,nglob

  integer ibool(npointot),iglob(npointot),loc(npointot)
  integer ind(npointot),ninseg(npointot)
  logical ifseg(npointot)
  double precision x(npointot),y(npointot),z(npointot)
  integer iwork(npointot)
  double precision work(npointot)

  integer ipoin,i,j
  integer nseg,ioff,iseg,ig
  double precision xtol

! establish initial pointers
  do ipoin=1,npointot
    loc(ipoin)=ipoin
  enddo

! define a tolerance, normalized radius is 1., so let's use a small value
  xtol = SMALLVALTOL

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
    call swap_all_buffers(ibool(ioff),loc(ioff), &
            x(ioff),y(ioff),z(ioff),iwork,work,ind,ninseg(iseg))
    ioff=ioff+ninseg(iseg)
  enddo

! check for jumps in current coordinate
  if(j == 1) then
    do i=2,npointot
      if(dabs(x(i)-x(i-1)) > xtol) ifseg(i)=.true.
    enddo
  else if(j == 2) then
    do i=2,npointot
      if(dabs(y(i)-y(i-1)) > xtol) ifseg(i)=.true.
    enddo
  else
    do i=2,npointot
      if(dabs(z(i)-z(i-1)) > xtol) ifseg(i)=.true.
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
    iglob(loc(i))=ig
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

  integer n
  double precision A(n)
  integer IND(n)

  integer i,j,l,ir,indx
  double precision q

  do j=1,n
    IND(j)=j
  enddo

  if(n == 1) return

  L=n/2+1
  ir=n
  100 CONTINUE
   IF(l>1) THEN
      l=l-1
      indx=ind(l)
      q=a(indx)
   ELSE
      indx=ind(ir)
      q=a(indx)
      ind(ir)=ind(1)
      ir=ir-1
      if (ir == 1) then
         ind(1)=indx
         return
      endif
   ENDIF
   i=l
   j=l+l
  200    CONTINUE
   IF(J <= IR) THEN
      IF(J < IR) THEN
         IF(A(IND(j)) < A(IND(j+1))) j=j+1
      ENDIF
      IF (q < A(IND(j))) THEN
         IND(I)=IND(J)
         I=J
         J=J+J
      ELSE
         J=IR+1
      ENDIF
   goto 200
   ENDIF
   IND(I)=INDX
  goto 100
  end subroutine rank_buffers

! -------------------------------------------------------------------

  subroutine swap_all_buffers(IA,IB,A,B,C,IW,W,ind,n)
!
! swap arrays IA, IB, A, B and C according to addressing in array IND
!
  implicit none

  integer n

  integer IND(n)
  integer IA(n),IB(n),IW(n)
  double precision A(n),B(n),C(n),W(n)

  integer i

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

