!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                    and University of Pau, France
! (c) California Institute of Technology and University of Pau, April 2007
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine get_MPI_cutplanes_xi(myrank,prname,nspec,iMPIcut_xi,ibool, &
                        xstore,ystore,zstore,mask_ibool,npointot, &
                        NSPEC2D_ETA)

! this routine detects cut planes along xi
! In principle the left cut plane of the first slice
! and the right cut plane of the last slice are not used
! in the solver except if we want to have periodic conditions

  implicit none

  include "constants.h"

  integer nspec,myrank
  integer NSPEC2D_ETA

  logical iMPIcut_xi(2,nspec)

  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

! logical mask used to create arrays iboolleft_xi and iboolright_xi
  integer npointot
  logical mask_ibool(npointot)

! global element numbering
  integer ispec

! MPI cut-plane element numbering
  integer ispecc1,ispecc2,npoin2D_xi,ix,iy,iz
  integer nspec2Dtheor

! processor identification
  character(len=150) prname,errmsg

! theoretical number of surface elements in the buffers
! cut planes along xi=constant correspond to ETA faces
      nspec2Dtheor = NSPEC2D_ETA

! write the MPI buffers for the left and right edges of the slice
! and the position of the points to check that the buffers are fine

!
! determine if the element falls on the left MPI cut plane
!

! global point number and coordinates left MPI cut-plane
  open(unit=10,file=prname(1:len_trim(prname))//'iboolleft_xi.txt',status='unknown')

! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

! nb of global points shared with the other slice
  npoin2D_xi = 0

! nb of elements in this cut-plane
  ispecc1=0

  do ispec=1,nspec
  if(iMPIcut_xi(1,ispec)) then

    ispecc1=ispecc1+1

! loop on all the points in that 2-D element, including edges
  ix = 1
  do iy=1,NGLLY
      do iz=1,NGLLZ

! select point, if not already selected
  if(.not. mask_ibool(ibool(ix,iy,iz,ispec))) then
      mask_ibool(ibool(ix,iy,iz,ispec)) = .true.
      npoin2D_xi = npoin2D_xi + 1

      write(10,*) ibool(ix,iy,iz,ispec),xstore(ix,iy,iz,ispec), &
              ystore(ix,iy,iz,ispec),zstore(ix,iy,iz,ispec)
  endif

      enddo
  enddo

  endif
  enddo

! put flag to indicate end of the list of points
  write(10,*) '0 0  0.  0.  0.'

! write total number of points
  write(10,*) npoin2D_xi

  close(10)

! compare number of surface elements detected to analytical value
  if(ispecc1 /= nspec2Dtheor) then
    write(errmsg,*) 'error MPI cut-planes detection in xi=left T=',nspec2Dtheor,' C=',ispecc1
    call exit_MPI(myrank,errmsg)
  endif
!
! determine if the element falls on the right MPI cut plane
!

! global point number and coordinates right MPI cut-plane
  open(unit=10,file=prname(1:len_trim(prname))//'iboolright_xi.txt',status='unknown')

! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

! nb of global points shared with the other slice
  npoin2D_xi = 0

! nb of elements in this cut-plane
  ispecc2=0

  do ispec=1,nspec
  if(iMPIcut_xi(2,ispec)) then

    ispecc2=ispecc2+1

! loop on all the points in that 2-D element, including edges
  ix = NGLLX
  do iy=1,NGLLY
      do iz=1,NGLLZ

! select point, if not already selected
  if(.not. mask_ibool(ibool(ix,iy,iz,ispec))) then
      mask_ibool(ibool(ix,iy,iz,ispec)) = .true.
      npoin2D_xi = npoin2D_xi + 1

      write(10,*) ibool(ix,iy,iz,ispec),xstore(ix,iy,iz,ispec), &
              ystore(ix,iy,iz,ispec),zstore(ix,iy,iz,ispec)
  endif

      enddo
  enddo

  endif
  enddo

! put flag to indicate end of the list of points
  write(10,*) '0 0  0.  0.  0.'

! write total number of points
  write(10,*) npoin2D_xi

  close(10)

! compare number of surface elements detected to analytical value
  if(ispecc2 /= nspec2Dtheor) then
    write(errmsg,*) 'error MPI cut-planes detection in xi=right T=',nspec2Dtheor,' C=',ispecc2
    call exit_MPI(myrank,errmsg)
  endif

  end subroutine get_MPI_cutplanes_xi

