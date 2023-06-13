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

! fix the non blocking arrays to assemble the slices inside each chunk: elements
! in contact with the MPI faces by an edge or a corner only but not
! a full face are missing, therefore let us add them
  subroutine fix_non_blocking_slices(is_on_a_slice_edge,iboolright_xi, &
                                     iboolleft_xi,iboolright_eta,iboolleft_eta, &
                                     npoin2D_xi,npoin2D_eta,ibool, &
                                     nspec,nglob,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX)

  use constants

  implicit none

  integer, intent(in) :: nspec,nglob,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX

  integer, dimension(NB_SQUARE_EDGES_ONEDIR), intent(in) :: npoin2D_xi,npoin2D_eta

  logical, dimension(nspec), intent(inout) :: is_on_a_slice_edge

  integer, dimension(NGLOB2DMAX_XMIN_XMAX), intent(in) :: iboolleft_xi,iboolright_xi
  integer, dimension(NGLOB2DMAX_YMIN_YMAX), intent(in) :: iboolleft_eta,iboolright_eta

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: ibool

  ! local parameters
  logical, dimension(:),allocatable :: mask_ibool
  integer :: ipoin,ispec,i,j,k,ier

  ! clean the mask
  allocate(mask_ibool(nglob),stat=ier)
  if (ier /= 0) stop 'Error allocating mask in fix_non_blocking_slices routine'
  mask_ibool(:) = .false.

  ! mark all the points that are in the MPI buffers to assemble inside each chunk
  do ipoin = 1,npoin2D_xi(1)
    mask_ibool(iboolleft_xi(ipoin)) = .true.
  enddo

  do ipoin = 1,npoin2D_eta(1)
    mask_ibool(iboolleft_eta(ipoin)) = .true.
  enddo

  do ipoin = 1,npoin2D_xi(2)
    mask_ibool(iboolright_xi(ipoin)) = .true.
  enddo

  do ipoin = 1,npoin2D_eta(2)
    mask_ibool(iboolright_eta(ipoin)) = .true.
  enddo

  ! now label all the elements that have at least one corner belonging
  ! to any of these buffers as elements that must contribute to the
  ! first step of the calculations (performed on the edges before starting
  ! the non blocking communications); there is no need to examine the inside
  ! of the elements, checking their eight corners is sufficient
  do ispec = 1,nspec
    do k = 1,NGLLZ,NGLLZ-1
      do j = 1,NGLLY,NGLLY-1
        do i = 1,NGLLX,NGLLX-1
          if (mask_ibool(ibool(i,j,k,ispec))) then
            is_on_a_slice_edge(ispec) = .true.
            goto 888
          endif
        enddo
      enddo
    enddo
  888 continue
  enddo

  deallocate(mask_ibool)

  end subroutine fix_non_blocking_slices

!
!-------------------------------------------------------------------------------------------------
!

! fix the non blocking arrays to assemble the central cube: elements
! in contact with the MPI faces by an edge or a corner only but not
! a full face are missing, therefore let us add them
  subroutine fix_non_blocking_central_cube(is_on_a_slice_edge, &
                                           ibool,nspec,nglob,nb_msgs_theor_in_cube,ibelm_bottom_inner_core, &
                                           idoubling_inner_core,npoin2D_cube_from_slices, &
                                           ibool_central_cube,NSPEC2D_BOTTOM_INNER_CORE, &
                                           ichunk,NPROC_XI)

  use constants

  implicit none

  integer, intent(in) :: nspec,nglob,nb_msgs_theor_in_cube,NSPEC2D_BOTTOM_INNER_CORE
  integer, intent(in) :: ichunk,npoin2D_cube_from_slices,NPROC_XI

  logical, dimension(nspec), intent(inout) :: is_on_a_slice_edge

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(in) :: ibool

  integer, dimension(nb_msgs_theor_in_cube,npoin2D_cube_from_slices), intent(in) :: ibool_central_cube

  integer, dimension(NSPEC2D_BOTTOM_INNER_CORE), intent(in) :: ibelm_bottom_inner_core

  integer, dimension(nspec), intent(in) :: idoubling_inner_core

  ! local parameters
  logical, dimension(:),allocatable :: mask_ibool
  integer :: ipoin,ispec,i,j,k,imsg,ispec2D,ier

  if (ichunk /= CHUNK_AB .and. ichunk /= CHUNK_AB_ANTIPODE) then
    do ispec2D = 1,NSPEC2D_BOTTOM_INNER_CORE
      ispec = ibelm_bottom_inner_core(ispec2D)
      is_on_a_slice_edge(ispec) = .true.
    enddo
  endif

  if (ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    do ispec = 1,nspec
      if (idoubling_inner_core(ispec) == IFLAG_BOTTOM_CENTRAL_CUBE .or. &
         idoubling_inner_core(ispec) == IFLAG_TOP_CENTRAL_CUBE) &
        is_on_a_slice_edge(ispec) = .true.
    enddo
  endif

  if (ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then

    ! clean the mask
    allocate(mask_ibool(nglob),stat=ier)
    if (ier /= 0) stop 'Error allocating mask_ibool array'
    mask_ibool(:) = .false.

    do imsg = 1,nb_msgs_theor_in_cube
      do ipoin = 1,npoin2D_cube_from_slices
        if (NPROC_XI == 1) then
          if (ibool_central_cube(imsg,ipoin) > 0) then
            mask_ibool(ibool_central_cube(imsg,ipoin)) = .true.
          endif
        else
          mask_ibool(ibool_central_cube(imsg,ipoin)) = .true.
        endif
      enddo
    enddo

    ! now label all the elements that have at least one corner belonging
    ! to any of these buffers as elements that must contribute to the
    ! first step of the calculations (performed on the edges before starting
    ! the non blocking communications); there is no need to examine the inside
    ! of the elements, checking their eight corners is sufficient
    do ispec = 1,nspec
      do k = 1,NGLLZ,NGLLZ-1
        do j = 1,NGLLY,NGLLY-1
          do i = 1,NGLLX,NGLLX-1
            if (mask_ibool(ibool(i,j,k,ispec))) then
              is_on_a_slice_edge(ispec) = .true.
              goto 888
            endif
          enddo
        enddo
      enddo
    888 continue
    enddo

    ! free array
    deallocate(mask_ibool)

  endif

  end subroutine fix_non_blocking_central_cube

