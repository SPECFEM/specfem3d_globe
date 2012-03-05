!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
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

  subroutine assemble_MPI_central_cube(ichunk,nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
                                      npoin2D_cube_from_slices, &
                                      buffer_all_cube_from_slices,buffer_slices, &
                                      request_send_cc,request_receive_cc, &
                                      request_send_array_cc,request_receive_array_cc, &
                                      ibool_central_cube, &
                                      receiver_cube_from_slices,ibool_inner_core,idoubling_inner_core, &
                                      ibelm_bottom_inner_core,NSPEC2D_BOTTOM_INNER_CORE, &
                                      vector_assemble,ndim_assemble, &
                                      iphase_comm_CC)

  implicit none

  ! standard include of the MPI library
  include 'mpif.h'
  include 'constants.h'

  ! include values created by the mesher
  include "OUTPUT_FILES/values_from_mesher.h"

  ! for matching with central cube in inner core
  integer, intent(in) :: ichunk, nb_msgs_theor_in_cube, npoin2D_cube_from_slices
  integer, intent(in) :: ndim_assemble
  integer, intent(in) :: receiver_cube_from_slices
  integer, intent(inout) :: iphase_comm_CC
  integer, dimension(nb_msgs_theor_in_cube), intent(in) :: sender_from_slices_to_cube

  double precision, dimension(npoin2D_cube_from_slices,ndim_assemble), intent(inout) :: &
    buffer_slices
  double precision, dimension(npoin2D_cube_from_slices,ndim_assemble,nb_msgs_theor_in_cube), intent(inout) :: &
    buffer_all_cube_from_slices

  ! note: these parameters are "saved" now as global parameters
  ! MPI status of messages to be received
  integer, intent(inout) :: request_send_cc,request_receive_cc
  ! maximum value of nb_msgs_theor_in_cube is 5 (when NPROC_XI == 1)
  ! therefore NPROC_XI+4 is always large enough
  integer, dimension(NPROC_XI_VAL+4), intent(inout) :: request_send_array_cc,request_receive_array_cc

  integer, dimension(nb_msgs_theor_in_cube,npoin2D_cube_from_slices), intent(in) :: ibool_central_cube

  ! local to global mapping
  integer, intent(in) :: NSPEC2D_BOTTOM_INNER_CORE
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), intent(in) :: ibool_inner_core
  integer, dimension(NSPEC_INNER_CORE), intent(in) :: idoubling_inner_core
  integer, dimension(NSPEC2D_BOTTOM_INNER_CORE), intent(in) :: ibelm_bottom_inner_core

  ! vector
  real(kind=CUSTOM_REAL), dimension(ndim_assemble,NGLOB_INNER_CORE), intent(inout) :: vector_assemble

! local parameters

  integer ipoin,idimension, ispec2D, ispec
  integer i,j,k
  integer sender,receiver,imsg

  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE) :: array_central_cube

  logical :: flag_result_test
  integer, dimension(MPI_STATUS_SIZE) :: msg_status
  integer :: ier

  ! mask
  logical, dimension(NGLOB_INNER_CORE) :: mask

  !---
  !---  use buffers to assemble mass matrix with central cube once and for all
  !---

  select case( iphase_comm_CC )

  case( 1 )

    ! non-central-cube chunks send values to receiver central cube chunks AB or AB_ANTIPODE

    ! note: only chunks AB and AB_ANTIPODE contain valid central cube elements,
    !          all other have only fictitious central cube elements

    ! on chunks AB and AB_ANTIPODE, receive all the messages from slices
    if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
      do imsg = 1,nb_msgs_theor_in_cube-1
        ! receive buffers from slices
        sender = sender_from_slices_to_cube(imsg)
        call MPI_IRECV(buffer_all_cube_from_slices(1,1,imsg), &
                  ndim_assemble*npoin2D_cube_from_slices,MPI_DOUBLE_PRECISION,sender, &
                  itag,MPI_COMM_WORLD,request_receive_array_cc(imsg),ier)
      enddo
    endif

    ! send info to central cube from all the slices except those in CHUNK_AB & CHUNK_AB_ANTIPODE
    if(ichunk /= CHUNK_AB .and. ichunk /= CHUNK_AB_ANTIPODE) then
      ! for bottom elements in contact with central cube from the slices side
      ipoin = 0
      do ispec2D = 1,NSPEC2D_BOTTOM_INNER_CORE
        ispec = ibelm_bottom_inner_core(ispec2D)
        ! only for DOFs exactly on surface of central cube (bottom of these elements)
        k = 1
        do j = 1,NGLLY
          do i = 1,NGLLX
            ipoin = ipoin + 1
            buffer_slices(ipoin,:) = dble(vector_assemble(1:ndim_assemble,ibool_inner_core(i,j,k,ispec)))
          enddo
        enddo
      enddo
      ! send buffer to central cube
      receiver = receiver_cube_from_slices
      call MPI_ISEND(buffer_slices,ndim_assemble*npoin2D_cube_from_slices, &
                MPI_DOUBLE_PRECISION,receiver,itag,MPI_COMM_WORLD,request_send_cc,ier)
    endif  ! end sending info to central cube

    iphase_comm_CC = iphase_comm_CC + 1
    return ! exit because we have started some communications therefore we need some time

  case( 2 )

    ! central cube chunks AB and AB_ANTIPODE send values to each other

    ! checks that chunks have sent out messages
    if(ichunk /= CHUNK_AB .and. ichunk /= CHUNK_AB_ANTIPODE) then
      call MPI_TEST(request_send_cc,flag_result_test,msg_status,ier)
      if(.not. flag_result_test) return ! exit if message not sent yet
    endif

    ! checks that central cube chunks have received all (requested) messages
    if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
      do imsg = 1,nb_msgs_theor_in_cube-1
        call MPI_TEST(request_receive_array_cc(imsg),flag_result_test,msg_status,ier)
        if(.not. flag_result_test) return ! exit if message not received yet
      enddo
    endif

    ! exchange of their bottom faces between chunks AB and AB_ANTIPODE
    if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
      ipoin = 0
      do ispec = NSPEC_INNER_CORE, 1, -1
        if (idoubling_inner_core(ispec) == IFLAG_BOTTOM_CENTRAL_CUBE) then
          k = 1
          do j = 1,NGLLY
            do i = 1,NGLLX
              ipoin = ipoin + 1
              buffer_slices(ipoin,:) = dble(vector_assemble(1:ndim_assemble,ibool_inner_core(i,j,k,ispec)))
            enddo
          enddo
        endif
      enddo
      sender = sender_from_slices_to_cube(nb_msgs_theor_in_cube)
    !   call MPI_SENDRECV(buffer_slices,ndim_assemble*npoin2D_cube_from_slices,MPI_DOUBLE_PRECISION,receiver_cube_from_slices, &
    !       itag,buffer_slices2,ndim_assemble*npoin2D_cube_from_slices,&
    !       MPI_DOUBLE_PRECISION,sender,itag,MPI_COMM_WORLD,msg_status,ier)

      call MPI_IRECV(buffer_all_cube_from_slices(1,1,nb_msgs_theor_in_cube), &
          ndim_assemble*npoin2D_cube_from_slices,MPI_DOUBLE_PRECISION,sender,itag,MPI_COMM_WORLD,request_receive_cc,ier)
    !! DK DK this merged with previous statement
    !   buffer_all_cube_from_slices(:,:,nb_msgs_theor_in_cube) = buffer_slices2(:,:)

      call MPI_ISEND(buffer_slices,ndim_assemble*npoin2D_cube_from_slices,MPI_DOUBLE_PRECISION,receiver_cube_from_slices, &
          itag,MPI_COMM_WORLD,request_send_cc,ier)
    endif

    iphase_comm_CC = iphase_comm_CC + 1
    return ! exit because we have started some communications therefore we need some time

  case( 3 )

    !--- now we need to assemble the contributions

    ! central cube chunks AB and AB_ANTIPODE assemble values and send them out to others

    if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then

      ! checks that messages between AB and AB_ANTIPODE have been sent and received
      call MPI_TEST(request_send_cc,flag_result_test,msg_status,ier)
      if(.not. flag_result_test) return ! exit if message not sent yet
      call MPI_TEST(request_receive_cc,flag_result_test,msg_status,ier)
      if(.not. flag_result_test) return ! exit if message not received yet

      do idimension = 1,ndim_assemble
        ! erase contributions to central cube array
        array_central_cube(:) = 0._CUSTOM_REAL

        ! use indirect addressing to store contributions only once
        ! distinguish between single and double precision for reals
        do imsg = 1,nb_msgs_theor_in_cube-1
          do ipoin = 1,npoin2D_cube_from_slices
            if(CUSTOM_REAL == SIZE_REAL) then
              array_central_cube(ibool_central_cube(imsg,ipoin)) = sngl(buffer_all_cube_from_slices(ipoin,idimension,imsg))
            else
              array_central_cube(ibool_central_cube(imsg,ipoin)) = buffer_all_cube_from_slices(ipoin,idimension,imsg)
            endif
          enddo
        enddo
        ! add the constribution of AB or AB_ANTIPODE to sum with the external slices on the edges
        ! use a mask to avoid taking the same point into account several times.
        mask(:) = .false.
        do ipoin = 1,npoin2D_cube_from_slices
          if (.not. mask(ibool_central_cube(nb_msgs_theor_in_cube,ipoin))) then
            if(CUSTOM_REAL == SIZE_REAL) then
              array_central_cube(ibool_central_cube(nb_msgs_theor_in_cube,ipoin)) = &
              array_central_cube(ibool_central_cube(nb_msgs_theor_in_cube,ipoin)) + &
              sngl(buffer_all_cube_from_slices(ipoin,idimension,nb_msgs_theor_in_cube))
            else
              array_central_cube(ibool_central_cube(nb_msgs_theor_in_cube,ipoin)) = &
              array_central_cube(ibool_central_cube(nb_msgs_theor_in_cube,ipoin)) + &
              buffer_all_cube_from_slices(ipoin,idimension,nb_msgs_theor_in_cube)
            endif
            mask(ibool_central_cube(nb_msgs_theor_in_cube,ipoin)) = .true.
          endif
        enddo

        ! suppress degrees of freedom already assembled at top of cube on edges
        do ispec = 1,NSPEC_INNER_CORE
          if(idoubling_inner_core(ispec) == IFLAG_TOP_CENTRAL_CUBE) then
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
            buffer_all_cube_from_slices(ipoin,idimension,imsg) = vector_assemble(idimension,ibool_central_cube(imsg,ipoin))
          enddo
        enddo

      enddo

    endif

    !----------

    ! receive info from central cube on all the slices except those in CHUNK_AB & CHUNK_AB_ANTIPODE
    if(ichunk /= CHUNK_AB .and. ichunk /= CHUNK_AB_ANTIPODE) then
      ! receive buffers from slices
      sender = receiver_cube_from_slices
      call MPI_IRECV(buffer_slices, &
                  ndim_assemble*npoin2D_cube_from_slices,MPI_DOUBLE_PRECISION,sender, &
                  itag,MPI_COMM_WORLD,request_receive_cc,ier)
      ! for bottom elements in contact with central cube from the slices side
      !   ipoin = 0
      !   do ispec2D = 1,NSPEC2D_BOTTOM_INNER_CORE
      !     ispec = ibelm_bottom_inner_core(ispec2D)
      ! only for DOFs exactly on surface of central cube (bottom of these elements)
      !     k = 1
      !     do j = 1,NGLLY
      !       do i = 1,NGLLX
      !         ipoin = ipoin + 1
      ! distinguish between single and double precision for reals
      !         if(CUSTOM_REAL == SIZE_REAL) then
      !           vector_assemble(:,ibool_inner_core(i,j,k,ispec)) = sngl(buffer_slices(ipoin,:))
      !         else
      !           vector_assemble(:,ibool_inner_core(i,j,k,ispec)) = buffer_slices(ipoin,:)
      !         endif
      !       enddo
      !     enddo
      !   enddo
    endif  ! end receiving info from central cube

    !------- send info back from central cube to slices

    ! on chunk AB & CHUNK_AB_ANTIPODE, send all the messages to slices
    if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
      do imsg = 1,nb_msgs_theor_in_cube-1
        ! send buffers to slices
        receiver = sender_from_slices_to_cube(imsg)
        call MPI_ISEND(buffer_all_cube_from_slices(1,1,imsg),ndim_assemble*npoin2D_cube_from_slices, &
                MPI_DOUBLE_PRECISION,receiver,itag,MPI_COMM_WORLD,request_send_array_cc(imsg),ier)
      enddo
    endif

    iphase_comm_CC = iphase_comm_CC + 1
    return ! exit because we have started some communications therefore we need some time

  case( 4 )

    ! all non-central cube chunks set the values at the common points with central cube

    ! checks that messages were sent out by central cube chunks AB and AB_ANTIPODE
    if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
      do imsg = 1,nb_msgs_theor_in_cube-1
        call MPI_TEST(request_send_array_cc(imsg),flag_result_test,msg_status,ier)
        if(.not. flag_result_test) return ! exit if message not sent yet
      enddo
    endif

    ! checks that messages have been received
    if(ichunk /= CHUNK_AB .and. ichunk /= CHUNK_AB_ANTIPODE) then
      call MPI_TEST(request_receive_cc,flag_result_test,msg_status,ier)
      if(.not. flag_result_test) return ! exit if message not received yet
    endif

    ! receive info from central cube on all the slices except those in CHUNK_AB & CHUNK_AB_ANTIPODE
    if(ichunk /= CHUNK_AB .and. ichunk /= CHUNK_AB_ANTIPODE) then
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
            if(CUSTOM_REAL == SIZE_REAL) then
              vector_assemble(1:ndim_assemble,ibool_inner_core(i,j,k,ispec)) = sngl(buffer_slices(ipoin,:))
            else
              vector_assemble(1:ndim_assemble,ibool_inner_core(i,j,k,ispec)) = buffer_slices(ipoin,:)
            endif
          enddo
        enddo
      enddo
    endif  ! end receiving info from central cube

    ! this is the exit condition, to go beyond the last phase number
    iphase_comm_CC = iphase_comm_CC + 1

  end select

  end subroutine assemble_MPI_central_cube

