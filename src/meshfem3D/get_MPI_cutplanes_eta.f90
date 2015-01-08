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

  subroutine get_MPI_cutplanes_eta(myrank,prname,nspec,iMPIcut_eta,ibool, &
                                  xstore,ystore,zstore,mask_ibool,npointot, &
                                  NSPEC2D_XI_FACE,iregion,npoin2D_eta, &
                                  iboolleft_eta,iboolright_eta, &
                                  npoin2D_eta_all,NGLOB2DMAX_YMIN_YMAX)

! this routine detects cut planes along eta
! In principle the left cut plane of the first slice
! and the right cut plane of the last slice are not used
! in the solver except if we want to have periodic conditions

  use constants

  implicit none

  integer :: nspec,myrank

  logical,dimension(2,nspec) :: iMPIcut_eta

  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  ! logical mask used to create arrays iboolleft_eta and iboolright_eta
  integer :: npointot
  logical,dimension(npointot) :: mask_ibool

  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_EDGES_ONEDIR) :: NSPEC2D_XI_FACE

  integer :: iregion
  integer :: npoin2D_eta

  integer :: NGLOB2DMAX_YMIN_YMAX
  integer, dimension(NGLOB2DMAX_YMIN_YMAX) :: iboolleft_eta,iboolright_eta

  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_eta_all

  ! processor identification
  character(len=MAX_STRING_LEN) :: prname

  ! local parameters
  ! global element numbering
  integer :: ispec

  ! MPI cut-plane element numbering
  integer :: ispecc1,ispecc2,ix,iy,iz
  integer :: nspec2Dtheor

  ! debug: file output
  logical,parameter :: DEBUG = .false.

  ! theoretical number of surface elements in the buffers
  ! cut planes along eta=constant correspond to XI faces
  nspec2Dtheor = NSPEC2D_XI_FACE(iregion,1)

! write the MPI buffers for the left and right edges of the slice
! and the position of the points to check that the buffers are fine

!
! determine if the element falls on the left MPI cut plane
!

  if (DEBUG) then
    ! global point number and coordinates left MPI cut-plane
    open(unit=IOUT,file=prname(1:len_trim(prname))//'iboolleft_eta.txt',status='unknown')
  endif

  ! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  ! nb of global points shared with the other slice
  iboolleft_eta(:) = 0
  npoin2D_eta = 0
  npoin2D_eta_all(1) = 1

  ! nb of elements in this cut-plane
  ispecc1 = 0

  do ispec = 1,nspec
    if (iMPIcut_eta(1,ispec)) then
      ispecc1=ispecc1+1
      ! loop on all the points in that 2-D element, including edges
      iy = 1
      do ix = 1,NGLLX
          do iz = 1,NGLLZ
            ! select point, if not already selected
            if (.not. mask_ibool(ibool(ix,iy,iz,ispec))) then
              mask_ibool(ibool(ix,iy,iz,ispec)) = .true.
              npoin2D_eta = npoin2D_eta + 1

              ! fills buffer arrays
              iboolleft_eta(npoin2D_eta) = ibool(ix,iy,iz,ispec)

              npoin2D_eta_all(1) = npoin2D_eta_all(1) + 1

              ! debug file output
              if (DEBUG) then
                write(IOUT,*) ibool(ix,iy,iz,ispec), xstore(ix,iy,iz,ispec), &
                              ystore(ix,iy,iz,ispec),zstore(ix,iy,iz,ispec)
              endif
            endif
          enddo
      enddo
    endif
  enddo

  if (DEBUG) then
    ! put flag to indicate end of the list of points
    write(IOUT,*) '0 0  0.  0.  0.'
    ! write total number of points
    write(IOUT,*) npoin2D_eta
    close(IOUT)
  endif

! compare number of surface elements detected to analytical value
  if (ispecc1 /= nspec2Dtheor) call exit_MPI(myrank,'Error MPI cut-planes detection in eta=left')

  ! subtract the line that contains the flag after the last point
  npoin2D_eta_all(1) = npoin2D_eta_all(1) - 1
  if (npoin2D_eta_all(1) > NGLOB2DMAX_YMIN_YMAX .or. npoin2D_eta_all(1) /= npoin2D_eta) &
    call exit_MPI(myrank,'incorrect iboolleft_eta read')


!
! determine if the element falls on the right MPI cut plane
!
  nspec2Dtheor = NSPEC2D_XI_FACE(iregion,2)

  if (DEBUG) then
    ! global point number and coordinates right MPI cut-plane
    open(unit=IOUT,file=prname(1:len_trim(prname))//'iboolright_eta.txt',status='unknown')
  endif

  ! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  ! nb of global points shared with the other slice
  iboolright_eta(:) = 0
  npoin2D_eta = 0
  npoin2D_eta_all(2) = 1

  ! nb of elements in this cut-plane
  ispecc2 = 0

  do ispec = 1,nspec
    if (iMPIcut_eta(2,ispec)) then
      ispecc2=ispecc2+1
      ! loop on all the points in that 2-D element, including edges
      iy = NGLLY
      do ix = 1,NGLLX
          do iz = 1,NGLLZ
          ! select point, if not already selected
          if (.not. mask_ibool(ibool(ix,iy,iz,ispec))) then
            mask_ibool(ibool(ix,iy,iz,ispec)) = .true.
            npoin2D_eta = npoin2D_eta + 1

            ! fills buffer arrays
            iboolright_eta(npoin2D_eta) = ibool(ix,iy,iz,ispec)

            npoin2D_eta_all(2) = npoin2D_eta_all(2) + 1

            ! debug file output
            if (DEBUG) then
              write(IOUT,*) ibool(ix,iy,iz,ispec), xstore(ix,iy,iz,ispec), &
                            ystore(ix,iy,iz,ispec),zstore(ix,iy,iz,ispec)
            endif
          endif
        enddo
      enddo
    endif
  enddo

  if (DEBUG) then
    ! put flag to indicate end of the list of points
    write(IOUT,*) '0 0  0.  0.  0.'
    ! write total number of points
    write(IOUT,*) npoin2D_eta
    close(IOUT)
  endif

  ! compare number of surface elements detected to analytical value
  if (ispecc2 /= nspec2Dtheor) call exit_MPI(myrank,'Error MPI cut-planes detection in eta=right')

  ! subtract the line that contains the flag after the last point
  npoin2D_eta_all(2) = npoin2D_eta_all(2) - 1
  if (npoin2D_eta_all(2) > NGLOB2DMAX_YMIN_YMAX .or. npoin2D_eta_all(2) /= npoin2D_eta) &
      call exit_MPI(myrank,'incorrect iboolright_eta read')

  end subroutine get_MPI_cutplanes_eta

