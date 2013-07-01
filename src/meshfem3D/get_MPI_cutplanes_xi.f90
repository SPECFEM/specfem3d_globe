!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            April 2011
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

  subroutine get_MPI_cutplanes_xi(myrank,prname,nspec,iMPIcut_xi,ibool, &
                                  xstore,ystore,zstore,mask_ibool,npointot, &
                                  NSPEC2D_ETA_FACE,iregion,npoin2D_xi)

! this routine detects cut planes along xi
! In principle the left cut plane of the first slice
! and the right cut plane of the last slice are not used
! in the solver except if we want to have periodic conditions

  implicit none

  include "constants.h"

  integer :: nspec,myrank

  logical,dimension(2,nspec) :: iMPIcut_xi

  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  ! logical mask used to create arrays iboolleft_xi and iboolright_xi
  integer :: npointot
  logical,dimension(npointot) :: mask_ibool

  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_EDGES_ONEDIR) :: NSPEC2D_ETA_FACE

  integer :: iregion
  integer :: npoin2D_xi

  ! processor identification
  character(len=150) :: prname

  ! local parameters
  ! global element numbering
  integer :: ispec
  ! MPI cut-plane element numbering
  integer :: ispecc1,ispecc2,ix,iy,iz
  integer :: nspec2Dtheor
  integer :: ier

  character(len=150) :: errmsg

  ! debug: file output
  logical,parameter :: DEBUG = .false.

  ! theoretical number of surface elements in the buffers
  ! cut planes along xi=constant correspond to ETA faces
  nspec2Dtheor = NSPEC2D_ETA_FACE(iregion,1)

! write the MPI buffers for the left and right edges of the slice
! and the position of the points to check that the buffers are fine

!
! determine if the element falls on the left MPI cut plane
!

! global point number and coordinates left MPI cut-plane
  open(unit=10,file=prname(1:len_trim(prname))//'iboolleft_xi.txt', &
       status='unknown',iostat=ier)
  if( ier /= 0 ) then
    if( myrank == 0 ) then
      write(IMAIN,*)
      write(IMAIN,*) 'error creating file: '
      write(IMAIN,*) prname(1:len_trim(prname))//'iboolleft_xi.txt'
      write(IMAIN,*)
      write(IMAIN,*) 'please make sure that the directory specified in Par_file as LOCAL_PATH exists'
      write(IMAIN,*)
    endif
    call exit_mpi(myrank,'error creating iboolleft_xi.txt, please check your Par_file LOCAL_PATH setting')
  endif

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
                write(10,*) ibool(ix,iy,iz,ispec), xstore(ix,iy,iz,ispec), &
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
  nspec2Dtheor = NSPEC2D_ETA_FACE(iregion,2)

! global point number and coordinates right MPI cut-plane
  open(unit=10,file=prname(1:len_trim(prname))//'iboolright_xi.txt', &
        status='unknown',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error creating iboolright_xi.txt for this process')

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
              write(10,*) ibool(ix,iy,iz,ispec), xstore(ix,iy,iz,ispec), &
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

