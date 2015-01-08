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

  subroutine assemble_MPI_central_cube_block(ichunk,nb_msgs_theor_in_cube, sender_from_slices_to_cube, &
                                          npoin2D_cube_from_slices, &
                                          buffer_all_cube_from_slices, buffer_slices, buffer_slices2, &
                                          ibool_central_cube, &
                                          receiver_cube_from_slices, ibool_inner_core, &
                                          idoubling_inner_core, NSPEC_INNER_CORE, &
                                          ibelm_bottom_inner_core, NSPEC2D_BOTTOM_INNER_CORE,NGLOB_INNER_CORE, &
                                          vector_assemble,ndim_assemble, &
                                          iproc_eta,addressing,NCHUNKS,NPROC_XI,NPROC_ETA)

  ! this version of the routine is based on blocking MPI calls

  use constants

  implicit none

  ! for matching with central cube in inner core
  integer ichunk, nb_msgs_theor_in_cube, npoin2D_cube_from_slices
  integer, dimension(nb_msgs_theor_in_cube) :: sender_from_slices_to_cube
  double precision, dimension(npoin2D_cube_from_slices,NDIM) :: &
    buffer_slices,buffer_slices2
  double precision, dimension(nb_msgs_theor_in_cube,npoin2D_cube_from_slices,NDIM) :: &
    buffer_all_cube_from_slices
  integer, dimension(nb_msgs_theor_in_cube,npoin2D_cube_from_slices):: ibool_central_cube
  integer receiver_cube_from_slices

  ! local to global mapping
  integer NSPEC_INNER_CORE,NSPEC2D_BOTTOM_INNER_CORE, NGLOB_INNER_CORE
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core
  integer, dimension(NSPEC_INNER_CORE) :: idoubling_inner_core
  integer, dimension(NSPEC2D_BOTTOM_INNER_CORE) :: ibelm_bottom_inner_core

  ! vector
  integer ndim_assemble
  real(kind=CUSTOM_REAL), dimension(ndim_assemble,NGLOB_INNER_CORE) :: vector_assemble

  !for addressing of the slices
  integer, intent(in) :: NCHUNKS,NPROC_XI,NPROC_ETA
  integer, dimension(NCHUNKS,0:NPROC_XI-1,0:NPROC_ETA-1), intent(in) :: addressing
  integer, intent(in) :: iproc_eta

  integer ipoin,idimension, ispec2D, ispec
  integer i,j,k
  integer sender,receiver,imsg

  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE) :: array_central_cube

  ! mask
  logical, dimension(NGLOB_INNER_CORE) :: mask

  !---
  !---  now use buffers to assemble mass matrix with central cube once and for all
  !---

  ! on chunks AB and AB_ANTIPODE, receive all the messages from slices
  if (ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then

    do imsg = 1,nb_msgs_theor_in_cube-1

  ! receive buffers from slices
    sender = sender_from_slices_to_cube(imsg)
    call recv_dp(buffer_slices,ndim_assemble*npoin2D_cube_from_slices,sender,itag)

  ! copy buffer in 2D array for each slice
    buffer_all_cube_from_slices(imsg,:,1:ndim_assemble) = buffer_slices(:,1:ndim_assemble)

    enddo
  endif

  ! send info to central cube from all the slices except those in CHUNK_AB & CHUNK_AB_ANTIPODE
  if (ichunk /= CHUNK_AB .and. ichunk /= CHUNK_AB_ANTIPODE) then
    ! for bottom elements in contact with central cube from the slices side
    ipoin = 0
    do ispec2D = 1,NSPEC2D_BOTTOM_INNER_CORE

      ispec = ibelm_bottom_inner_core(ispec2D)

      ! only for DOFs exactly on surface of central cube (bottom of these elements)
      k = 1
      do j = 1,NGLLY
        do i = 1,NGLLX
          ipoin = ipoin + 1
          buffer_slices(ipoin,1:ndim_assemble) = dble(vector_assemble(1:ndim_assemble,ibool_inner_core(i,j,k,ispec)))
        enddo
      enddo
    enddo

    ! send buffer to central cube
    receiver = receiver_cube_from_slices
    call send_dp(buffer_slices,ndim_assemble*npoin2D_cube_from_slices,receiver,itag)

    ! in case NPROC_XI == 1, the other chunks exchange all bottom points with
    ! CHUNK_AB **and** CHUNK_AB_ANTIPODE
    if (NPROC_XI == 1) then
      call send_dp(buffer_slices,ndim_assemble*npoin2D_cube_from_slices,addressing(CHUNK_AB_ANTIPODE,0,iproc_eta),itag)
    endif

  endif  ! end sending info to central cube


  ! exchange of their bottom faces between chunks AB and AB_ANTIPODE
  if (ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then

    ipoin = 0
    do ispec = NSPEC_INNER_CORE, 1, -1
      if (idoubling_inner_core(ispec) == IFLAG_BOTTOM_CENTRAL_CUBE) then
        k = 1
        do j = 1,NGLLY
          do i = 1,NGLLX
            ipoin = ipoin + 1
            buffer_slices(ipoin,1:ndim_assemble) = dble(vector_assemble(1:ndim_assemble,ibool_inner_core(i,j,k,ispec)))
          enddo
        enddo
      endif
    enddo

    sender = sender_from_slices_to_cube(nb_msgs_theor_in_cube)

    call sendrecv_dp(buffer_slices,ndim_assemble*npoin2D_cube_from_slices,receiver_cube_from_slices,itag, &
                     buffer_slices2,ndim_assemble*npoin2D_cube_from_slices,sender,itag)

   buffer_all_cube_from_slices(nb_msgs_theor_in_cube,:,1:ndim_assemble) = buffer_slices2(:,1:ndim_assemble)

  endif

  !--- now we need to assemble the contributions

  if (ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then

    do idimension = 1,ndim_assemble
  ! erase contributions to central cube array
      array_central_cube(:) = 0._CUSTOM_REAL

  ! use indirect addressing to store contributions only once
  ! distinguish between single and double precision for reals
      do imsg = 1,nb_msgs_theor_in_cube-1
        do ipoin = 1,npoin2D_cube_from_slices
          if (NPROC_XI == 1) then
            if (ibool_central_cube(imsg,ipoin) > 0) then
              array_central_cube(ibool_central_cube(imsg,ipoin)) = real(buffer_all_cube_from_slices(imsg,ipoin,idimension), &
                                                                        kind=CUSTOM_REAL)
            endif
          else
            array_central_cube(ibool_central_cube(imsg,ipoin)) = real(buffer_all_cube_from_slices(imsg,ipoin,idimension), &
                                                                      kind=CUSTOM_REAL)
          endif
        enddo
      enddo
  ! add the contribution of AB or AB_ANTIPODE to sum with the external slices on the edges
  ! use a mask to avoid taking the same point into account several times.
      mask(:) = .false.
      do ipoin = 1,npoin2D_cube_from_slices
        if (NPROC_XI == 1) then
          if (ibool_central_cube(nb_msgs_theor_in_cube,ipoin) > 0) then
            if (.not. mask(ibool_central_cube(nb_msgs_theor_in_cube,ipoin))) then
              array_central_cube(ibool_central_cube(nb_msgs_theor_in_cube,ipoin)) = &
                array_central_cube(ibool_central_cube(nb_msgs_theor_in_cube,ipoin)) + &
                real(buffer_all_cube_from_slices(nb_msgs_theor_in_cube,ipoin,idimension), kind=CUSTOM_REAL)
            endif
            mask(ibool_central_cube(nb_msgs_theor_in_cube,ipoin)) = .true.
          endif
        else
          if (.not. mask(ibool_central_cube(nb_msgs_theor_in_cube,ipoin))) then
            array_central_cube(ibool_central_cube(nb_msgs_theor_in_cube,ipoin)) = &
              array_central_cube(ibool_central_cube(nb_msgs_theor_in_cube,ipoin)) + &
              real(buffer_all_cube_from_slices(nb_msgs_theor_in_cube,ipoin,idimension), kind=CUSTOM_REAL)
          endif
          mask(ibool_central_cube(nb_msgs_theor_in_cube,ipoin)) = .true.
        endif
      enddo

  ! suppress degrees of freedom already assembled at top of cube on edges
      do ispec = 1,NSPEC_INNER_CORE
        if (idoubling_inner_core(ispec) == IFLAG_TOP_CENTRAL_CUBE) then
          k = NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              array_central_cube(ibool_inner_core(i,j,k,ispec)) = 0._CUSTOM_REAL
            enddo
          enddo
        endif
      enddo

  ! assemble contributions
      vector_assemble(idimension,:) = vector_assemble(idimension,:) + array_central_cube(:)

  ! copy sum back
      do imsg = 1,nb_msgs_theor_in_cube-1
        do ipoin = 1,npoin2D_cube_from_slices
          if (NPROC_XI == 1) then
            if (ibool_central_cube(imsg,ipoin) > 0) then
              buffer_all_cube_from_slices(imsg,ipoin,idimension) = &
                      vector_assemble(idimension,ibool_central_cube(imsg,ipoin))
            else
              buffer_all_cube_from_slices(imsg,ipoin,idimension) = 0._CUSTOM_REAL
            endif
          else
            buffer_all_cube_from_slices(imsg,ipoin,idimension) = &
                    vector_assemble(idimension,ibool_central_cube(imsg,ipoin))
          endif
        enddo
      enddo

    enddo

  endif

  !----------

  ! receive info from central cube on all the slices except those in CHUNK_AB & CHUNK_AB_ANTIPODE
  if (ichunk /= CHUNK_AB .and. ichunk /= CHUNK_AB_ANTIPODE) then
    ! receive buffers from slices
    sender = receiver_cube_from_slices
    call recv_dp(buffer_slices,ndim_assemble*npoin2D_cube_from_slices,sender,itag)

    ! in case NPROC_XI == 1, the other chunks exchange all bottom points with
    ! CHUNK_AB **and** CHUNK_AB_ANTIPODE
    if (NPROC_XI == 1) then
      call recv_dp(buffer_slices2,ndim_assemble*npoin2D_cube_from_slices,addressing(CHUNK_AB_ANTIPODE,0,iproc_eta),itag)

      buffer_slices = buffer_slices + buffer_slices2
    endif

    ! for bottom elements in contact with central cube from the slices side
    ipoin = 0
    do ispec2D = 1,NSPEC2D_BOTTOM_INNER_CORE

      ispec = ibelm_bottom_inner_core(ispec2D)

      ! only for DOFs exactly on surface of central cube (bottom of these elements)
      k = 1
      do j = 1,NGLLY
        do i = 1,NGLLX
          ipoin = ipoin + 1

          ! distinguish between single and double precision for reals
          vector_assemble(1:ndim_assemble,ibool_inner_core(i,j,k,ispec)) = real(buffer_slices(ipoin,1:ndim_assemble), &
                                                                                kind=CUSTOM_REAL)

        enddo
      enddo
    enddo

  endif  ! end receiving info from central cube

  !------- send info back from central cube to slices

  ! on chunk AB & CHUNK_AB_ANTIPODE, send all the messages to slices
  if (ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then

   do imsg = 1,nb_msgs_theor_in_cube-1

  ! copy buffer in 2D array for each slice
   buffer_slices(:,1:ndim_assemble) = buffer_all_cube_from_slices(imsg,:,1:ndim_assemble)

  ! send buffers to slices
    receiver = sender_from_slices_to_cube(imsg)
    call send_dp(buffer_slices,ndim_assemble*npoin2D_cube_from_slices,receiver,itag)

   enddo
   endif

  end subroutine assemble_MPI_central_cube_block

