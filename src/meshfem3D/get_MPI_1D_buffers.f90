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

  subroutine get_MPI_1D_buffers(myrank,prname,nspec,iMPIcut_xi,iMPIcut_eta, &
                                ibool,idoubling,xstore,ystore,zstore,mask_ibool,npointot, &
                                NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER,iregion, &
                                ibool1D_leftxi_lefteta,ibool1D_rightxi_lefteta, &
                                ibool1D_leftxi_righteta,ibool1D_rightxi_righteta, &
                                xyz1D_leftxi_lefteta,xyz1D_rightxi_lefteta, &
                                xyz1D_leftxi_righteta,xyz1D_rightxi_righteta, &
                                NGLOB1D_RADIAL_MAX)

! routine to create the MPI 1D chunk buffers for edges

  use constants

  implicit none

  integer :: nspec,myrank

  logical,dimension(2,nspec) :: iMPIcut_xi,iMPIcut_eta

  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  integer,dimension(nspec) :: idoubling

  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  ! logical mask used to create arrays ibool1D
  integer :: npointot
  logical,dimension(npointot) :: mask_ibool

  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_CORNERS) :: NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER
  integer :: iregion

  integer :: NGLOB1D_RADIAL_MAX
  integer,dimension(NGLOB1D_RADIAL_MAX) :: ibool1D_leftxi_lefteta,ibool1D_rightxi_lefteta, &
                                           ibool1D_leftxi_righteta,ibool1D_rightxi_righteta
  double precision,dimension(NGLOB1D_RADIAL_MAX,NDIM) :: xyz1D_leftxi_lefteta,xyz1D_rightxi_lefteta, &
                                                         xyz1D_leftxi_righteta,xyz1D_rightxi_righteta

  ! processor identification
  character(len=MAX_STRING_LEN) :: prname

  ! local parameters
  ! global element numbering
  integer :: ispec
  ! MPI 1D buffer element numbering
  integer :: ispeccount,npoin1D,ix,iy,iz

  ! debug file output
  logical,parameter :: DEBUG = .false.

! write the MPI buffers for the left and right edges of the slice
! and the position of the points to check that the buffers are fine

! *****************************************************************
! ****************** generate for eta = eta_min *******************
! *****************************************************************

! determine if the element falls on the left MPI cut plane

  if (DEBUG) then
    ! global point number and coordinates left MPI 1D buffer
    open(unit=IOUT,file=prname(1:len_trim(prname))//'ibool1D_leftxi_lefteta.txt',status='unknown')
  endif

  ! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  ! nb of global points shared with the other slice
  npoin1D = 0
  ibool1D_leftxi_lefteta(:) = 0

  ! nb of elements in this 1D buffer
  ispeccount = 0

  do ispec = 1,nspec
    ! remove central cube for chunk buffers
    if (idoubling(ispec) == IFLAG_MIDDLE_CENTRAL_CUBE .or. &
      idoubling(ispec) == IFLAG_BOTTOM_CENTRAL_CUBE .or. &
      idoubling(ispec) == IFLAG_TOP_CENTRAL_CUBE .or. &
      idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

    ! corner detection here
    if (iMPIcut_xi(1,ispec) .and. iMPIcut_eta(1,ispec)) then
      ispeccount=ispeccount+1
      ! loop on all the points
      ix = 1
      iy = 1
      do iz = 1,NGLLZ
        ! select point, if not already selected
        if (.not. mask_ibool(ibool(ix,iy,iz,ispec))) then
          ! adds this point
          mask_ibool(ibool(ix,iy,iz,ispec)) = .true.
          npoin1D = npoin1D + 1

          ! fills buffer array
          ibool1D_leftxi_lefteta(npoin1D) = ibool(ix,iy,iz,ispec)
          xyz1D_leftxi_lefteta(npoin1D,1) = xstore(ix,iy,iz,ispec)
          xyz1D_leftxi_lefteta(npoin1D,2) = ystore(ix,iy,iz,ispec)
          xyz1D_leftxi_lefteta(npoin1D,3) = zstore(ix,iy,iz,ispec)

          ! debug file output
          if (DEBUG) then
            write(IOUT,*) ibool(ix,iy,iz,ispec), xstore(ix,iy,iz,ispec), &
                          ystore(ix,iy,iz,ispec),zstore(ix,iy,iz,ispec)
          endif
        endif
      enddo
    endif
  enddo

  if (DEBUG) then
    ! put flag to indicate end of the list of points
    write(IOUT,*) '0  0  0.  0.  0.'
    ! write total number of points
    write(IOUT,*) npoin1D
    close(IOUT)
  endif

  ! compare number of edge elements detected to analytical value
  if (ispeccount /= NSPEC1D_RADIAL_CORNER(iregion,1) .or. npoin1D /= NGLOB1D_RADIAL_CORNER(iregion,1)) &
    call exit_MPI(myrank,'Error MPI 1D buffer detection in xi=left')

! determine if the element falls on the right MPI cut plane

  if (DEBUG) then
    ! global point number and coordinates right MPI 1D buffer
    open(unit=IOUT,file=prname(1:len_trim(prname))//'ibool1D_rightxi_lefteta.txt',status='unknown')
  endif

  ! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  ! nb of global points shared with the other slice
  npoin1D = 0
  ibool1D_rightxi_lefteta(:) = 0

  ! nb of elements in this 1D buffer
  ispeccount = 0
  do ispec = 1,nspec
    ! remove central cube for chunk buffers
    if (idoubling(ispec) == IFLAG_MIDDLE_CENTRAL_CUBE .or. &
      idoubling(ispec) == IFLAG_BOTTOM_CENTRAL_CUBE .or. &
      idoubling(ispec) == IFLAG_TOP_CENTRAL_CUBE .or. &
      idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

    ! corner detection here
    if (iMPIcut_xi(2,ispec) .and. iMPIcut_eta(1,ispec)) then
      ispeccount=ispeccount+1
      ! loop on all the points
      ix = NGLLX
      iy = 1
      do iz = 1,NGLLZ
        ! select point, if not already selected
        if (.not. mask_ibool(ibool(ix,iy,iz,ispec))) then
          mask_ibool(ibool(ix,iy,iz,ispec)) = .true.
          npoin1D = npoin1D + 1

          ! fills buffer array
          ibool1D_rightxi_lefteta(npoin1D) = ibool(ix,iy,iz,ispec)
          xyz1D_rightxi_lefteta(npoin1D,1) = xstore(ix,iy,iz,ispec)
          xyz1D_rightxi_lefteta(npoin1D,2) = ystore(ix,iy,iz,ispec)
          xyz1D_rightxi_lefteta(npoin1D,3) = zstore(ix,iy,iz,ispec)

          ! debug file output
          if (DEBUG) then
            write(IOUT,*) ibool(ix,iy,iz,ispec), xstore(ix,iy,iz,ispec), &
                          ystore(ix,iy,iz,ispec),zstore(ix,iy,iz,ispec)
          endif
        endif
      enddo
    endif
  enddo

  if (DEBUG) then
    ! put flag to indicate end of the list of points
    write(IOUT,*) '0  0  0.  0.  0.'
    ! write total number of points
    write(IOUT,*) npoin1D
    close(IOUT)
  endif

  ! compare number of edge elements and points detected to analytical value
  if (ispeccount /= NSPEC1D_RADIAL_CORNER(iregion,2) .or. npoin1D /= NGLOB1D_RADIAL_CORNER(iregion,2)) &
    call exit_MPI(myrank,'Error MPI 1D buffer detection in xi=right')

! *****************************************************************
! ****************** generate for eta = eta_max *******************
! *****************************************************************

! determine if the element falls on the left MPI cut plane

  if (DEBUG) then
    ! global point number and coordinates left MPI 1D buffer
    open(unit=IOUT,file=prname(1:len_trim(prname))//'ibool1D_leftxi_righteta.txt',status='unknown')
  endif

  ! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  ! nb of global points shared with the other slice
  npoin1D = 0
  ibool1D_leftxi_righteta(:) = 0

  ! nb of elements in this 1D buffer
  ispeccount = 0

  do ispec = 1,nspec

    ! remove central cube for chunk buffers
    if (idoubling(ispec) == IFLAG_MIDDLE_CENTRAL_CUBE .or. &
      idoubling(ispec) == IFLAG_BOTTOM_CENTRAL_CUBE .or. &
      idoubling(ispec) == IFLAG_TOP_CENTRAL_CUBE .or. &
      idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

    ! corner detection here
    if (iMPIcut_xi(1,ispec) .and. iMPIcut_eta(2,ispec)) then

      ispeccount=ispeccount+1

      ! loop on all the points
      ix = 1
      iy = NGLLY
      do iz = 1,NGLLZ
        ! select point, if not already selected
        if (.not. mask_ibool(ibool(ix,iy,iz,ispec))) then
          mask_ibool(ibool(ix,iy,iz,ispec)) = .true.
          npoin1D = npoin1D + 1

          ! fills buffer array
          ibool1D_leftxi_righteta(npoin1D) = ibool(ix,iy,iz,ispec)
          xyz1D_leftxi_righteta(npoin1D,1) = xstore(ix,iy,iz,ispec)
          xyz1D_leftxi_righteta(npoin1D,2) = ystore(ix,iy,iz,ispec)
          xyz1D_leftxi_righteta(npoin1D,3) = zstore(ix,iy,iz,ispec)

          ! debug file output
          if (DEBUG) then
            write(IOUT,*) ibool(ix,iy,iz,ispec), xstore(ix,iy,iz,ispec), &
                          ystore(ix,iy,iz,ispec),zstore(ix,iy,iz,ispec)
          endif
        endif
      enddo
    endif
  enddo

  if (DEBUG) then
    ! put flag to indicate end of the list of points
    write(IOUT,*) '0  0  0.  0.  0.'
    ! write total number of points
    write(IOUT,*) npoin1D
    close(IOUT)
  endif

  ! compare number of edge elements detected to analytical value
  if (ispeccount /= NSPEC1D_RADIAL_CORNER(iregion,4) .or. npoin1D /= NGLOB1D_RADIAL_CORNER(iregion,4)) &
    call exit_MPI(myrank,'Error MPI 1D buffer detection in xi=left')

! determine if the element falls on the right MPI cut plane

  if (DEBUG) then
    ! global point number and coordinates right MPI 1D buffer
    open(unit=IOUT,file=prname(1:len_trim(prname))//'ibool1D_rightxi_righteta.txt',status='unknown')
  endif

  ! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  ! nb of global points shared with the other slice
  npoin1D = 0
  ibool1D_rightxi_righteta(:) = 0

  ! nb of elements in this 1D buffer
  ispeccount = 0

  do ispec = 1,nspec

    ! remove central cube for chunk buffers
    if (idoubling(ispec) == IFLAG_MIDDLE_CENTRAL_CUBE .or. &
      idoubling(ispec) == IFLAG_BOTTOM_CENTRAL_CUBE .or. &
      idoubling(ispec) == IFLAG_TOP_CENTRAL_CUBE .or. &
      idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

    ! corner detection here
    if (iMPIcut_xi(2,ispec) .and. iMPIcut_eta(2,ispec)) then

      ispeccount=ispeccount+1

      ! loop on all the points
      ix = NGLLX
      iy = NGLLY
      do iz = 1,NGLLZ

        ! select point, if not already selected
        if (.not. mask_ibool(ibool(ix,iy,iz,ispec))) then
          mask_ibool(ibool(ix,iy,iz,ispec)) = .true.
          npoin1D = npoin1D + 1

          ! fills buffer array
          ibool1D_rightxi_righteta(npoin1D) = ibool(ix,iy,iz,ispec)
          xyz1D_rightxi_righteta(npoin1D,1) = xstore(ix,iy,iz,ispec)
          xyz1D_rightxi_righteta(npoin1D,2) = ystore(ix,iy,iz,ispec)
          xyz1D_rightxi_righteta(npoin1D,3) = zstore(ix,iy,iz,ispec)

          ! debug file output
          if (DEBUG) then
            write(IOUT,*) ibool(ix,iy,iz,ispec), xstore(ix,iy,iz,ispec), &
                          ystore(ix,iy,iz,ispec),zstore(ix,iy,iz,ispec)
          endif
        endif
      enddo
    endif
  enddo

  if (DEBUG) then
    ! put flag to indicate end of the list of points
    write(IOUT,*) '0  0  0.  0.  0.'
    ! write total number of points
    write(IOUT,*) npoin1D
    close(IOUT)
  endif

  ! compare number of edge elements and points detected to analytical value
  if (ispeccount /= NSPEC1D_RADIAL_CORNER(iregion,3) .or. npoin1D /= NGLOB1D_RADIAL_CORNER(iregion,3)) &
    call exit_MPI(myrank,'Error MPI 1D buffer detection in xi=right')

  end subroutine get_MPI_1D_buffers

