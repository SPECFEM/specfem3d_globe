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

  subroutine get_global(npointot,xp,yp,zp,iglob,locval,ifseg,nglob)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

! non-structured global numbering software provided by Paul F. Fischer

  use constants

  implicit none

  ! input parameters
  integer, intent(in) :: npointot

  double precision, dimension(npointot), intent(inout) :: xp,yp,zp

  integer, dimension(npointot), intent(out) :: iglob,locval
  logical, dimension(npointot), intent(out) :: ifseg
  integer, intent(out) :: nglob

  ! local variables
  integer, dimension(:), allocatable :: ninseg,idummy
  integer :: ier

  ! dynamically allocate arrays
  allocate(ninseg(npointot),stat=ier)
  if (ier /= 0) stop 'Error while allocating ninseg'

  allocate(idummy(npointot),stat=ier)
  if (ier /= 0) stop 'Error while allocating idummy'

  call sort_array_coordinates(npointot,xp,yp,zp,idummy,iglob,locval,ifseg, &
                              nglob,ninseg)

  ! deallocate arrays
  deallocate(ninseg,idummy)

  end subroutine get_global

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_global_indirect_addressing(nspec,nglob,ibool)

!
!- we can create a new indirect addressing to reduce cache misses
! (put into this subroutine but compiler keeps on complaining that it can't vectorize loops...)

  use constants

  implicit none

  integer,intent(in) :: nspec,nglob
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  ! local parameters
  ! mask to sort ibool
  integer, dimension(:), allocatable :: mask_ibool
  integer, dimension(:,:,:,:), allocatable :: copy_ibool_ori
  integer :: inumber
  integer:: i,j,k,ispec,ier

  ! copies original array
  allocate(copy_ibool_ori(NGLLX,NGLLY,NGLLZ,nspec), &
          mask_ibool(nglob), &
          stat=ier)
  if (ier /= 0) stop 'Error allocating local arrays in get_global_indirect_addressing'

  ! initializes arrays
  mask_ibool(:) = -1
  copy_ibool_ori(:,:,:,:) = ibool(:,:,:,:)

  ! reduces misses
  inumber = 0
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          if (mask_ibool(copy_ibool_ori(i,j,k,ispec)) == -1) then
            ! creates a new point
            inumber = inumber + 1
            ibool(i,j,k,ispec) = inumber
            mask_ibool(copy_ibool_ori(i,j,k,ispec)) = inumber
          else
            ! uses an existing point created previously
            ibool(i,j,k,ispec) = mask_ibool(copy_ibool_ori(i,j,k,ispec))
          endif
        enddo
      enddo
    enddo
  enddo

  ! cleanup
  deallocate(copy_ibool_ori)
  deallocate(mask_ibool)

  end subroutine get_global_indirect_addressing

