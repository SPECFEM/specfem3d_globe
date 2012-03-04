
! extract only the points that are needed for the new doubling brick
! from the modified version of Emanuele's tripling brick

  program extract_unique_points

  implicit none

  integer, parameter :: OLD_NPOINTS = 32, NEW_NPOINTS = 23

! there are seven spectral elements in the basic doubling brick
  integer, parameter :: NSPEC = 7

! bricks are made of eight points
  integer, parameter :: NGNOD = 8

! compute maximum number of points for call to get_global()
  integer, parameter :: npointot = nspec * NGNOD

  integer, dimension(0:OLD_NPOINTS-1) :: new_number_of_unique_point

! point numbers of all the nodes of each spectral element
  integer, dimension(NGNOD) :: old_number

  integer :: ispec,ipoin

  double precision, dimension(0:OLD_NPOINTS-1) :: x,y,z

! variables for call to get_global()
  integer, dimension(npointot) :: iglob,locval
  logical, dimension(npointot) :: ifseg
  double precision, dimension(npointot) :: xp,yp,zp
  integer nglob,ieoff,ilocnum

! define flags for points that must be excluded
  new_number_of_unique_point(:) = -1

  new_number_of_unique_point(0)  = 00
  new_number_of_unique_point(1)  = 01
  new_number_of_unique_point(2)  = 02
  new_number_of_unique_point(3)  = 03
  new_number_of_unique_point(4)  = 04
  new_number_of_unique_point(5)  = 05
  new_number_of_unique_point(6)  = 06
  new_number_of_unique_point(7)  = 07
  new_number_of_unique_point(8)  = 08
  new_number_of_unique_point(9)  = 09
  new_number_of_unique_point(10) = 10
  new_number_of_unique_point(11) = 11
  new_number_of_unique_point(12) = 12
  new_number_of_unique_point(13) = 13
  new_number_of_unique_point(14) = 14
  new_number_of_unique_point(15) = 15
  new_number_of_unique_point(16) = 16
  new_number_of_unique_point(17) = 17
  new_number_of_unique_point(22) = 18
  new_number_of_unique_point(23) = 19
  new_number_of_unique_point(25) = 20
  new_number_of_unique_point(27) = 21
  new_number_of_unique_point(29) = 22

! read the modified brick
  open(unit=27,file='dimitri_brick_with_duplicates.dx',status='old',action='read')

! skip header
  read(27,*)

! read the OLD_NPOINTS old modified points, print them only if they are in the new list
  do ipoin = 0,OLD_NPOINTS-1
    read(27,*) x(ipoin),y(ipoin),z(ipoin)
    if(new_number_of_unique_point(ipoin) >= 0) print *,x(ipoin),y(ipoin),z(ipoin)
  enddo

! skip header
  read(27,*)

! read the NGNOD points of each of the NSPEC old elements and write their new numbers
! also fill the arrays for future call to get_global()
  do ispec = 1, NSPEC
    ieoff = NGNOD * (ispec-1)
    ilocnum = 0
    read(27,*) (old_number(ipoin), ipoin = 1, NGNOD)
    write(*,*) (new_number_of_unique_point(old_number(ipoin)), ipoin = 1, NGNOD)
    do ipoin=1,NGNOD
      ilocnum = ilocnum + 1
      xp(ilocnum+ieoff) = x(old_number(ipoin))
      yp(ilocnum+ieoff) = y(old_number(ipoin))
      zp(ilocnum+ieoff) = z(old_number(ipoin))
    enddo
  enddo

  close(27)

! call get_global() to check that all the points are unique
  call get_global(NSPEC,xp,yp,zp,iglob,locval,ifseg,nglob,npointot,NGNOD)
  print *,'number of unique points detected = ',nglob
  if(nglob /= NEW_NPOINTS) then
    stop 'error in the total number of unique points'
  else
    print *,'this number is OK'
  endif

  end program extract_unique_points

!
!=====================================================================
!

  subroutine get_global(nspec,xp,yp,zp,iglob,loc,ifseg,nglob,npointot,NGNOD)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

! non-structured global numbering software provided by Paul F. Fischer

! leave sorting subroutines in same source file to allow for inlining

  implicit none

! 3-D simulation
  integer, parameter :: NDIM = 3

! geometry tolerance parameter to calculate number of independent grid points
! sensitive to actual size of model, assumes reference size of 1
  double precision, parameter :: SMALLVALTOL = 1.d-10

  integer npointot,NGNOD
  integer iglob(npointot),loc(npointot)
  logical ifseg(npointot)
  double precision xp(npointot),yp(npointot),zp(npointot)
  integer nspec,nglob

  integer ispec,i,j
  integer ieoff,ilocnum,nseg,ioff,iseg,ig

  integer, dimension(:), allocatable :: ind,ninseg,iwork
  double precision, dimension(:), allocatable :: work

! dynamically allocate arrays
  allocate(ind(npointot))
  allocate(ninseg(npointot))
  allocate(iwork(npointot))
  allocate(work(npointot))

! establish initial pointers
  do ispec=1,nspec
    ieoff=NGNOD * (ispec-1)
    do ilocnum=1,NGNOD
      loc(ilocnum+ieoff)=ilocnum+ieoff
    enddo
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
      call rank(xp(ioff),ind,ninseg(iseg))
    else if(j == 2) then
      call rank(yp(ioff),ind,ninseg(iseg))
    else
      call rank(zp(ioff),ind,ninseg(iseg))
    endif
    call swap_all(loc(ioff),xp(ioff),yp(ioff),zp(ioff),iwork,work,ind,ninseg(iseg))
    ioff=ioff+ninseg(iseg)
  enddo

! check for jumps in current coordinate
! compare the coordinates of the points within a small tolerance
  if(j == 1) then
    do i=2,npointot
      if(abs(xp(i)-xp(i-1)) > SMALLVALTOL) ifseg(i)=.true.
    enddo
  else if(j == 2) then
    do i=2,npointot
      if(abs(yp(i)-yp(i-1)) > SMALLVALTOL) ifseg(i)=.true.
    enddo
  else
    do i=2,npointot
      if(abs(zp(i)-zp(i-1)) > SMALLVALTOL) ifseg(i)=.true.
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

! deallocate arrays
  deallocate(ind)
  deallocate(ninseg)
  deallocate(iwork)
  deallocate(work)

  end subroutine get_global

! -----------------------------------

! sorting routines put in same file to allow for inlining

  subroutine rank(A,IND,N)
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

  if (n == 1) return

  L=n/2+1
  ir=n
  100 CONTINUE
   IF (l>1) THEN
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
   IF (J <= IR) THEN
      IF (J<IR) THEN
         IF ( A(IND(j))<A(IND(j+1)) ) j=j+1
      ENDIF
      IF (q<A(IND(j))) THEN
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
  end subroutine rank

! ------------------------------------------------------------------

  subroutine swap_all(IA,A,B,C,IW,W,ind,n)
!
! swap arrays IA, A, B and C according to addressing in array IND
!
  implicit none

  integer n

  integer IND(n)
  integer IA(n),IW(n)
  double precision A(n),B(n),C(n),W(n)

  integer i

  IW(:) = IA(:)
  W(:) = A(:)

  do i=1,n
    IA(i)=IW(ind(i))
    A(i)=W(ind(i))
  enddo

  W(:) = B(:)

  do i=1,n
    B(i)=W(ind(i))
  enddo

  W(:) = C(:)

  do i=1,n
    C(i)=W(ind(i))
  enddo

  end subroutine swap_all

