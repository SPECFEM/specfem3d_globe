
! extract only the points that are needed for the new doubling brick
! from the modified version of Emanuele's tripling brick

  program create_superbrick_duplicates

  implicit none

! there are twenty-seven unique points in the basic doubling brick
  integer, parameter :: NPOINTS = 27

! there are eight spectral elements in the basic doubling brick
  integer, parameter :: NSPEC = 8

! bricks are made of eight points
  integer, parameter :: NGNOD = 8

! number of basic blocks in superbrick
  integer, parameter :: NBLOCKS = 4

! compute maximum number of points for call to get_global()
  integer, parameter :: npointot = NBLOCKS * NSPEC * NGNOD

  integer :: ispec,ipoin,ignod,iblock

! coordinates of all the mesh points
  double precision, dimension(0:NBLOCKS*NPOINTS-1) :: x,y,z

! point numbers of all the nodes of each spectral element
  integer, dimension(NGNOD,NBLOCKS*NSPEC) :: point_number

! variables for call to get_global()
  integer, dimension(npointot) :: iglob,locval
  logical, dimension(npointot) :: ifseg
  double precision, dimension(npointot) :: xp,yp,zp
  integer nglob,ieoff,ilocnum

! read the basic brick
  open(unit=11,file='dimitri_brick_without_duplicates_two_layers.dx',status='old',action='read')

! skip header
  read(11,*)

! read the points
  do ipoin = 0,NPOINTS-1

    read(11,*) x(ipoin),y(ipoin),z(ipoin)

! generate the first symmetric block (z is always unchanged)
    x(ipoin + 1*NPOINTS) = + x(ipoin)
    y(ipoin + 1*NPOINTS) = - y(ipoin)
    z(ipoin + 1*NPOINTS) = + z(ipoin)

! generate the second symmetric block (z is always unchanged)
    x(ipoin + 2*NPOINTS) = - x(ipoin)
    y(ipoin + 2*NPOINTS) = + y(ipoin)
    z(ipoin + 2*NPOINTS) = + z(ipoin)

! generate the third symmetric block (z is always unchanged)
    x(ipoin + 3*NPOINTS) = - x(ipoin)
    y(ipoin + 3*NPOINTS) = - y(ipoin)
    z(ipoin + 3*NPOINTS) = + z(ipoin)

  enddo

! skip header
  read(11,*)

! read the NGNOD points of each of the NSPEC elements
  do ispec = 1, NSPEC

    read(11,*) (point_number(ignod,ispec), ignod = 1, NGNOD)

! create the three symmetric elements
    point_number(:,ispec + 1*NSPEC) = point_number(:,ispec) + 1*NPOINTS
    point_number(:,ispec + 2*NSPEC) = point_number(:,ispec) + 2*NPOINTS
    point_number(:,ispec + 3*NSPEC) = point_number(:,ispec) + 3*NPOINTS

  enddo

  close(11)

! shift the four blocks from [-1,1] to [0,2] in the two horizontal directions
  x(:) = x(:) + 1.d0
  y(:) = y(:) + 1.d0

  do ispec=1,NBLOCKS*NSPEC
    ieoff = NGNOD * (ispec-1)
    ilocnum = 0
    do ignod = 1,NGNOD
      ilocnum = ilocnum + 1
      xp(ilocnum+ieoff) = x(point_number(ignod,ispec))
      yp(ilocnum+ieoff) = y(point_number(ignod,ispec))
      zp(ilocnum+ieoff) = z(point_number(ignod,ispec))
    enddo
  enddo

! call get_global() to check that all the points are unique
  call get_global(NBLOCKS*NSPEC,xp,yp,zp,iglob,locval,ifseg,nglob,npointot,NGNOD)
  print *,'number of unique points detected in superbrick = ',nglob,' instead of ',NBLOCKS*NPOINTS,' with multiples'

! create an OpenDX file with the new superbrick

! write OpenDX header with element data
  open(unit=11,file='dimitri_superbrick_with_duplicates.dx',status='unknown')
  write(11,*) 'object 1 class array type float rank 1 shape 3 items ',NBLOCKS*NPOINTS,' data follows'

! output all the points, with duplicates
  do ipoin = 0,NBLOCKS*NPOINTS-1
    write(11,*) sngl(x(ipoin)),sngl(y(ipoin)),sngl(z(ipoin))
  enddo

! write element header
  write(11,*) 'object 2 class array type int rank 1 shape 8 items ',NBLOCKS*NSPEC,' data follows'

! output all the elements
  do ispec = 1,NBLOCKS*NSPEC
! point order in OpenDX is 4,1,8,5,3,2,7,6, *not* 1,2,3,4,5,6,7,8 as in AVS
    write(11,"(i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6)") (point_number(ignod,ispec), ignod = 1, NGNOD)
  enddo

! output OpenDX header for data
! label for hexahedra in OpenDX is "cubes"
  write(11,*) 'attribute "element type" string "cubes"'
  write(11,*) 'attribute "ref" string "positions"'
  write(11,*) 'object 3 class array type float rank 0 items ',NBLOCKS*NSPEC,' data follows'

! write element data (use the same colors for each of the four blocks)
  do iblock = 1,NBLOCKS
    do ispec=1,nspec
      write(11,*) ispec
    enddo
  enddo

! define OpenDX field
  write(11,*) 'attribute "dep" string "connections"'
  write(11,*) 'object "irregular positions irregular connections" class field'
  write(11,*) 'component "positions" value 1'
  write(11,*) 'component "connections" value 2'
  write(11,*) 'component "data" value 3'
  write(11,*) 'end'

  close(11)

  end program create_superbrick_duplicates

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

