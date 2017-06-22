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

  subroutine get_MPI_cutplanes_xi(prname,nspec,iMPIcut_xi,ibool, &
                                  xstore,ystore,zstore,mask_ibool,npointot, &
                                  NSPEC2D_ETA_FACE,iregion,npoin2D_xi, &
                                  iboolleft_xi,iboolright_xi, &
                                  npoin2D_xi_all,NGLOB2DMAX_XMIN_XMAX)

! this routine detects cut planes along xi
! In principle the left cut plane of the first slice
! and the right cut plane of the last slice are not used
! in the solver except if we want to have periodic conditions

  use constants

  implicit none

  integer :: nspec

  logical,dimension(2,nspec) :: iMPIcut_xi

  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  ! logical mask used to create arrays iboolleft_xi and iboolright_xi
  integer :: npointot
  logical,dimension(npointot) :: mask_ibool

  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_EDGES_ONEDIR) :: NSPEC2D_ETA_FACE

  integer :: iregion
  integer :: npoin2D_xi

  integer :: NGLOB2DMAX_XMIN_XMAX
  integer, dimension(NGLOB2DMAX_XMIN_XMAX) :: iboolleft_xi,iboolright_xi

  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_all

  ! processor identification
  character(len=MAX_STRING_LEN) :: prname

  ! local parameters
  ! global element numbering
  integer :: ispec
  ! MPI cut-plane element numbering
  integer :: ispecc1,ispecc2,ix,iy,iz
  integer :: nspec2Dtheor
  integer :: ier

  character(len=MAX_STRING_LEN) :: errmsg

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

  if (DEBUG) then
    ! global point number and coordinates left MPI cut-plane
    open(unit=IOUT,file=prname(1:len_trim(prname))//'iboolleft_xi.txt',status='unknown',iostat=ier)
    if (ier /= 0) then
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'Error creating file: '
        write(IMAIN,*) prname(1:len_trim(prname))//'iboolleft_xi.txt'
        write(IMAIN,*)
        write(IMAIN,*) 'please make sure that the directory specified in Par_file as LOCAL_PATH exists'
        write(IMAIN,*)
        call flush_IMAIN()
      endif
      call exit_mpi(myrank,'Error creating iboolleft_xi.txt, please check your Par_file LOCAL_PATH setting')
    endif
  endif

  ! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  ! nb of global points shared with the other slice
  iboolleft_xi(:) = 0
  npoin2D_xi = 0
  npoin2D_xi_all(1) = 1

  ! nb of elements in this cut-plane
  ispecc1 = 0
  do ispec = 1,nspec
    if (iMPIcut_xi(1,ispec)) then
      ispecc1=ispecc1+1
      ! loop on all the points in that 2-D element, including edges
      ix = 1
      do iy = 1,NGLLY
          do iz = 1,NGLLZ
            ! select point, if not already selected
            if (.not. mask_ibool(ibool(ix,iy,iz,ispec))) then
              mask_ibool(ibool(ix,iy,iz,ispec)) = .true.
              npoin2D_xi = npoin2D_xi + 1

              ! fills buffer arrays
              iboolleft_xi(npoin2D_xi) = ibool(ix,iy,iz,ispec)

              npoin2D_xi_all(1) = npoin2D_xi_all(1) + 1

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
    write(IOUT,*) npoin2D_xi
    close(IOUT)
  endif

  ! compare number of surface elements detected to analytical value
  if (ispecc1 /= nspec2Dtheor) then
    write(errmsg,*) 'Error MPI cut-planes detection in xi=left T=',nspec2Dtheor,' C=',ispecc1
    call exit_MPI(myrank,errmsg)
  endif

  ! subtract the line that contains the flag after the last point
  npoin2D_xi_all(1) = npoin2D_xi_all(1) - 1
  if (npoin2D_xi_all(1) > NGLOB2DMAX_XMIN_XMAX .or. npoin2D_xi_all(1) /= npoin2D_xi) &
    call exit_MPI(myrank,'incorrect iboolleft_xi read')


!
! determine if the element falls on the right MPI cut plane
!
  nspec2Dtheor = NSPEC2D_ETA_FACE(iregion,2)

  if (DEBUG) then
    ! global point number and coordinates right MPI cut-plane
    open(unit=IOUT,file=prname(1:len_trim(prname))//'iboolright_xi.txt',status='unknown',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error creating iboolright_xi.txt for this process')
  endif

  ! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  ! nb of global points shared with the other slice
  iboolright_xi(:) = 0
  npoin2D_xi = 0
  npoin2D_xi_all(2) = 1

  ! nb of elements in this cut-plane
  ispecc2 = 0
  do ispec = 1,nspec
    if (iMPIcut_xi(2,ispec)) then
      ispecc2=ispecc2+1
      ! loop on all the points in that 2-D element, including edges
      ix = NGLLX
      do iy = 1,NGLLY
        do iz = 1,NGLLZ
          ! select point, if not already selected
          if (.not. mask_ibool(ibool(ix,iy,iz,ispec))) then
            mask_ibool(ibool(ix,iy,iz,ispec)) = .true.
            npoin2D_xi = npoin2D_xi + 1

            ! fills buffer arrays
            iboolright_xi(npoin2D_xi) = ibool(ix,iy,iz,ispec)

            npoin2D_xi_all(2) = npoin2D_xi_all(2) + 1

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
    write(IOUT,*) npoin2D_xi
    close(IOUT)
  endif

  ! compare number of surface elements detected to analytical value
  if (ispecc2 /= nspec2Dtheor) then
    write(errmsg,*) 'Error MPI cut-planes detection in xi=right T=',nspec2Dtheor,' C=',ispecc2
    call exit_MPI(myrank,errmsg)
  endif

  ! subtract the line that contains the flag after the last point
  npoin2D_xi_all(2) = npoin2D_xi_all(2) - 1
  if (npoin2D_xi_all(2) > NGLOB2DMAX_XMIN_XMAX .or. npoin2D_xi_all(2) /= npoin2D_xi) &
      call exit_MPI(myrank,'incorrect iboolright_xi read')

  end subroutine get_MPI_cutplanes_xi

