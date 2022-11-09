!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
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

  subroutine get_global(npointot,xp,yp,zp,iglob,nglob_new)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

! non-structured global numbering software provided by Paul F. Fischer

#ifdef USE_PARALLEL_STL_SORTING
  use constants, only: SMALLVALTOL
#endif

  implicit none

  ! input parameters
  integer, intent(in) :: npointot

  double precision, dimension(npointot), intent(inout) :: xp,yp,zp

  integer, dimension(npointot), intent(inout) :: iglob
  integer, intent(out) :: nglob_new

#ifdef USE_PARALLEL_STL_SORTING
  ! calls C++ routine w/ parallel STL routines
  call sort_array_coordinates_c(npointot,xp,yp,zp,iglob,SMALLVALTOL,nglob_new)

#else
  ! default sorting/ibool
  ! local variables
  integer, dimension(:), allocatable :: ninseg,idummy
  integer, dimension(:), allocatable :: locval
  logical, dimension(:), allocatable :: ifseg

  integer :: ier

  ! initializes
  nglob_new = 0

  ! checks if anything to do
  if (npointot == 0) return

  ! dynamically allocate arrays
  allocate(ninseg(npointot),stat=ier)
  if (ier /= 0) stop 'Error while allocating ninseg'
  ninseg(:) = 0

  allocate(idummy(npointot),stat=ier)
  if (ier /= 0) stop 'Error while allocating idummy'
  idummy(:) = 0

  allocate(locval(npointot), &
           ifseg(npointot),stat=ier)
  if (ier /= 0) stop 'Error in allocate 20a'
  locval(:) = 0
  ifseg(:) = .false.

  call sort_array_coordinates(npointot,xp,yp,zp,idummy,iglob,locval,ifseg,nglob_new,ninseg)

  ! deallocate arrays
  deallocate(ninseg,idummy)
  deallocate(locval,ifseg)
#endif

  end subroutine get_global

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_global_indirect_addressing(nspec,nglob,ibool)

!
!- we can create a new indirect addressing to reduce cache misses
! (put into this subroutine but compiler keeps on complaining that it can't vectorize loops...)

  use constants, only: NGLLX,NGLLY,NGLLZ

  implicit none

  integer,intent(in) :: nspec,nglob
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout) :: ibool

  ! local parameters
  ! mask to sort ibool
  integer, dimension(:), allocatable :: mask_ibool
  integer, dimension(:,:,:,:), allocatable :: copy_ibool_ori
  integer :: inumber
  integer:: i,j,k,ispec,ier,iglob

  ! checks if anything to do
  if (nspec == 0) return

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
          iglob = copy_ibool_ori(i,j,k,ispec)
          if (mask_ibool(iglob) == -1) then
            ! creates a new point
            inumber = inumber + 1
            ibool(i,j,k,ispec) = inumber
            mask_ibool(iglob) = inumber
          else
            ! uses an existing point created previously
            ibool(i,j,k,ispec) = mask_ibool(iglob)
          endif
        enddo
      enddo
    enddo
  enddo

  ! cleanup
  deallocate(copy_ibool_ori)
  deallocate(mask_ibool)

  end subroutine get_global_indirect_addressing


!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_global_indirect_addressing_phases(nspec,nglob,ibool, &
                                                   phase_ispec_inner,num_phase_ispec,nspec_inner,nspec_outer)

!
! creates indirect addressing to reduce cache misses following inner/outer elements
!

  use constants, only: NGLLX,NGLLY,NGLLZ

  implicit none

  integer,intent(in) :: nspec,nglob
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout) :: ibool

  integer,intent(in) :: num_phase_ispec
  integer, dimension(num_phase_ispec,2),intent(in) :: phase_ispec_inner
  integer,intent(in) :: nspec_inner,nspec_outer

  ! local parameters
  ! mask to sort ibool
  integer, dimension(:), allocatable :: mask_ibool
  integer, dimension(:,:,:,:), allocatable :: copy_ibool_ori
  integer :: inumber
  integer :: i,j,k,ispec,ier,iglob
  integer :: ispec_p,iphase,counter_nspec,num_elements

  ! checks if anything to do
  if (nspec == 0) return

  ! copies original array
  allocate(copy_ibool_ori(NGLLX,NGLLY,NGLLZ,nspec), &
           mask_ibool(nglob), &
           stat=ier)
  if (ier /= 0) stop 'Error allocating local arrays in get_global_indirect_addressing'

  ! initializes arrays
  mask_ibool(:) = -1
  copy_ibool_ori(:,:,:,:) = ibool(:,:,:,:)

  ! reduces misses
  counter_nspec = 0
  inumber = 0
  do iphase = 1,2

    if (iphase == 1) then
      ! outer elements (halo region)
      num_elements = nspec_outer
    else
      ! inner elements
      num_elements = nspec_inner
    endif

    do ispec_p = 1,num_elements
      ! counter
      counter_nspec = counter_nspec + 1

      ! only compute elements which belong to current phase (inner or outer elements)
      ispec = phase_ispec_inner(ispec_p,iphase)

      ! checks
      if (ispec < 1 .or. ispec > nspec ) stop 'Error ispec bounds in get_global_indirect_addressing_phases()'

      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = copy_ibool_ori(i,j,k,ispec)
            if (mask_ibool(iglob) == -1) then
              ! creates a new point
              inumber = inumber + 1
              ibool(i,j,k,ispec) = inumber
              mask_ibool(iglob) = inumber
            else
              ! uses an existing point created previously
              ibool(i,j,k,ispec) = mask_ibool(iglob)
            endif
          enddo
        enddo
      enddo

    enddo ! ispec_p

  enddo ! iphase

  ! checks
  if (counter_nspec /= nspec) stop 'Error invalid nspec counted in get_global_indirect_addressing_phases()'

  end subroutine get_global_indirect_addressing_phases


