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

!----
!---- assemble the contributions between slices and chunks using MPI
!----

  subroutine assemble_MPI_scalar(myrank,array_val,nglob, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi,iboolright_xi,iboolleft_eta,iboolright_eta, &
            npoin2D_faces,npoin2D_xi,npoin2D_eta, &
            iboolfaces,iboolcorner, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar, &
            npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            request_send,request_receive, &
            request_send_array,request_receive_array, &
            NUMMSGS_FACES,NCORNERSCHUNKS, &
            NPROC_XI,NPROC_ETA,NGLOB1D_RADIAL, &
            NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB2DMAX_XY,NCHUNKS, &
            iphase_comm)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer myrank,nglob,NCHUNKS,iphase_comm

! array to assemble
  real(kind=CUSTOM_REAL), dimension(nglob), intent(inout) :: array_val

  integer, intent(in) :: iproc_xi,iproc_eta,ichunk
  integer, dimension(NB_SQUARE_EDGES_ONEDIR), intent(in) :: npoin2D_xi,npoin2D_eta
  integer, intent(in) :: npoin2D_faces(NUMFACES_SHARED)

  integer, intent(in) :: NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB2DMAX_XY
  integer, intent(in) :: NPROC_XI,NPROC_ETA,NGLOB1D_RADIAL
  integer, intent(in) :: NUMMSGS_FACES,NCORNERSCHUNKS

! for addressing of the slices
  integer, dimension(NCHUNKS,0:NPROC_XI-1,0:NPROC_ETA-1), intent(in) :: addressing

! 2-D addressing and buffers for summation between slices
  integer, dimension(NGLOB2DMAX_XMIN_XMAX), intent(in) :: iboolleft_xi,iboolright_xi
  integer, dimension(NGLOB2DMAX_YMIN_YMAX), intent(in) :: iboolleft_eta,iboolright_eta

! indirect addressing for each corner of the chunks
  integer, dimension(NGLOB1D_RADIAL,NUMCORNERS_SHARED), intent(in) :: iboolcorner

  integer, intent(in) :: npoin2D_max_all_CM_IC
  integer, dimension(NGLOB2DMAX_XY,NUMFACES_SHARED), intent(in) :: iboolfaces
  real(kind=CUSTOM_REAL), dimension(npoin2D_max_all_CM_IC,NUMFACES_SHARED), intent(inout) :: &
    buffer_send_faces_scalar, &
    buffer_received_faces_scalar

  ! buffers for send and receive between corners of the chunks
  real(kind=CUSTOM_REAL), dimension(NGLOB1D_RADIAL), intent(inout) :: buffer_send_chunkcorn_scalar, &
                                                                      buffer_recv_chunkcorn_scalar

  ! ---- arrays to assemble between chunks

  ! communication pattern for faces between chunks
  integer, dimension(NUMMSGS_FACES), intent(in) :: iprocfrom_faces,iprocto_faces

  ! communication pattern for corners between chunks
  integer, dimension(NCORNERSCHUNKS), intent(in) :: iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners

  ! stored as global parameters
  integer, intent(inout) :: request_send,request_receive
  integer, dimension(NUMFACES_SHARED), intent(inout) :: request_send_array,request_receive_array

! local parameters

  ! MPI status of messages to be received
  integer, dimension(MPI_STATUS_SIZE) :: msg_status

  integer :: ipoin,ipoin2D,ipoin1D
  integer :: sender,receiver
  integer :: imsg
  integer :: icount_faces,npoin2D_chunks
  integer :: icount_corners

  integer :: ier
  logical :: flag_result_test


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! check flag to see if we need to assemble (might be turned off when debugging)
  if (.not. ACTUALLY_ASSEMBLE_MPI_SLICES)  then
    iphase_comm = 9999 ! this means everything is finished
    return
  endif

! here we have to assemble all the contributions between slices using MPI

!----
!---- assemble the contributions between slices using MPI
!----

!----
!---- first assemble along xi using the 2-D topology
!----

  if(iphase_comm == 1) then

    ! slices copy the right face into the buffer
    do ipoin=1,npoin2D_xi(2)
      buffer_send_faces_scalar(ipoin,1) = array_val(iboolright_xi(ipoin))
    enddo

    ! send messages forward along each row
    if(iproc_xi == 0) then
      sender = MPI_PROC_NULL
    else
      sender = addressing(ichunk,iproc_xi - 1,iproc_eta)
    endif
    if(iproc_xi == NPROC_XI-1) then
      receiver = MPI_PROC_NULL
    else
      receiver = addressing(ichunk,iproc_xi + 1,iproc_eta)
    endif
    call MPI_IRECV(buffer_received_faces_scalar,npoin2D_xi(1),CUSTOM_MPI_TYPE,sender, &
          itag,MPI_COMM_WORLD,request_receive,ier)

    call MPI_ISEND(buffer_send_faces_scalar,npoin2D_xi(2),CUSTOM_MPI_TYPE,receiver, &
          itag2,MPI_COMM_WORLD,request_send,ier)

    iphase_comm = iphase_comm + 1
    return ! exit because we have started some communications therefore we need some time

  endif !!!!!!!!! end of iphase_comm 1

  if(iphase_comm == 2) then

    ! call MPI_WAIT(request_send,msg_status,ier)
    ! call MPI_WAIT(request_receive,msg_status,ier)
    call MPI_TEST(request_send,flag_result_test,msg_status,ier)
    if(.not. flag_result_test) return ! exit if message not sent yet
    call MPI_TEST(request_receive,flag_result_test,msg_status,ier)
    if(.not. flag_result_test) return ! exit if message not received yet

    ! all slices add the buffer received to the contributions on the left face
    if(iproc_xi > 0) then
      do ipoin=1,npoin2D_xi(1)
        array_val(iboolleft_xi(ipoin)) = array_val(iboolleft_xi(ipoin)) + &
                                buffer_received_faces_scalar(ipoin,1)
      enddo
    endif

    ! the contributions are correctly assembled on the left side of each slice
    ! now we have to send the result back to the sender
    ! all slices copy the left face into the buffer
    do ipoin=1,npoin2D_xi(1)
      buffer_send_faces_scalar(ipoin,1) = array_val(iboolleft_xi(ipoin))
    enddo

    ! send messages backward along each row
    if(iproc_xi == NPROC_XI-1) then
      sender = MPI_PROC_NULL
    else
      sender = addressing(ichunk,iproc_xi + 1,iproc_eta)
    endif
    if(iproc_xi == 0) then
      receiver = MPI_PROC_NULL
    else
      receiver = addressing(ichunk,iproc_xi - 1,iproc_eta)
    endif
    call MPI_IRECV(buffer_received_faces_scalar,npoin2D_xi(2),CUSTOM_MPI_TYPE,sender, &
          itag,MPI_COMM_WORLD,request_receive,ier)

    call MPI_ISEND(buffer_send_faces_scalar,npoin2D_xi(1),CUSTOM_MPI_TYPE,receiver, &
          itag2,MPI_COMM_WORLD,request_send,ier)

    iphase_comm = iphase_comm + 1
    return ! exit because we have started some communications therefore we need some time

  endif !!!!!!!!! end of iphase_comm 2

  if(iphase_comm == 3) then

    ! call MPI_WAIT(request_send,msg_status,ier)
    ! call MPI_WAIT(request_receive,msg_status,ier)
    call MPI_TEST(request_send,flag_result_test,msg_status,ier)
    if(.not. flag_result_test) return ! exit if message not sent yet
    call MPI_TEST(request_receive,flag_result_test,msg_status,ier)
    if(.not. flag_result_test) return ! exit if message not received yet

    ! all slices copy the buffer received to the contributions on the right face
    if(iproc_xi < NPROC_XI-1) then
      do ipoin=1,npoin2D_xi(2)
        array_val(iboolright_xi(ipoin)) = buffer_received_faces_scalar(ipoin,1)
      enddo
    endif

    !----
    !---- then assemble along eta using the 2-D topology
    !----

    ! slices copy the right face into the buffer
    do ipoin=1,npoin2D_eta(2)
      buffer_send_faces_scalar(ipoin,1) = array_val(iboolright_eta(ipoin))
    enddo

    ! send messages forward along each row
    if(iproc_eta == 0) then
      sender = MPI_PROC_NULL
    else
      sender = addressing(ichunk,iproc_xi,iproc_eta - 1)
    endif
    if(iproc_eta == NPROC_ETA-1) then
      receiver = MPI_PROC_NULL
    else
      receiver = addressing(ichunk,iproc_xi,iproc_eta + 1)
    endif
    call MPI_IRECV(buffer_received_faces_scalar,npoin2D_eta(1),CUSTOM_MPI_TYPE,sender, &
      itag,MPI_COMM_WORLD,request_receive,ier)

    call MPI_ISEND(buffer_send_faces_scalar,npoin2D_eta(2),CUSTOM_MPI_TYPE,receiver, &
      itag2,MPI_COMM_WORLD,request_send,ier)

    iphase_comm = iphase_comm + 1
    return ! exit because we have started some communications therefore we need some time

  endif !!!!!!!!! end of iphase_comm 3

  if(iphase_comm == 4) then

    ! call MPI_WAIT(request_send,msg_status,ier)
    ! call MPI_WAIT(request_receive,msg_status,ier)
    call MPI_TEST(request_send,flag_result_test,msg_status,ier)
    if(.not. flag_result_test) return ! exit if message not sent yet
    call MPI_TEST(request_receive,flag_result_test,msg_status,ier)
    if(.not. flag_result_test) return ! exit if message not received yet

    ! all slices add the buffer received to the contributions on the left face
    if(iproc_eta > 0) then
      do ipoin=1,npoin2D_eta(1)
        array_val(iboolleft_eta(ipoin)) = array_val(iboolleft_eta(ipoin)) + &
                                buffer_received_faces_scalar(ipoin,1)
      enddo
    endif

    ! the contributions are correctly assembled on the left side of each slice
    ! now we have to send the result back to the sender
    ! all slices copy the left face into the buffer
    do ipoin=1,npoin2D_eta(1)
      buffer_send_faces_scalar(ipoin,1) = array_val(iboolleft_eta(ipoin))
    enddo

    ! send messages backward along each row
    if(iproc_eta == NPROC_ETA-1) then
      sender = MPI_PROC_NULL
    else
      sender = addressing(ichunk,iproc_xi,iproc_eta + 1)
    endif
    if(iproc_eta == 0) then
      receiver = MPI_PROC_NULL
    else
      receiver = addressing(ichunk,iproc_xi,iproc_eta - 1)
    endif
    call MPI_IRECV(buffer_received_faces_scalar,npoin2D_eta(2),CUSTOM_MPI_TYPE,sender, &
      itag,MPI_COMM_WORLD,request_receive,ier)

    call MPI_ISEND(buffer_send_faces_scalar,npoin2D_eta(1),CUSTOM_MPI_TYPE,receiver, &
      itag2,MPI_COMM_WORLD,request_send,ier)

    iphase_comm = iphase_comm + 1
    return ! exit because we have started some communications therefore we need some time

  endif !!!!!!!!! end of iphase_comm 4

  if(iphase_comm == 5) then

    ! call MPI_WAIT(request_send,msg_status,ier)
    ! call MPI_WAIT(request_receive,msg_status,ier)
    call MPI_TEST(request_send,flag_result_test,msg_status,ier)
    if(.not. flag_result_test) return ! exit if message not sent yet
    call MPI_TEST(request_receive,flag_result_test,msg_status,ier)
    if(.not. flag_result_test) return ! exit if message not received yet

    ! all slices copy the buffer received to the contributions on the right face
    if(iproc_eta < NPROC_ETA-1) then
      do ipoin=1,npoin2D_eta(2)
        array_val(iboolright_eta(ipoin)) = buffer_received_faces_scalar(ipoin,1)
      enddo
    endif

    !----
    !---- start MPI assembling phase between chunks
    !----

    ! check flag to see if we need to assemble (might be turned off when debugging)
    ! and do not assemble if only one chunk
    if (.not. ACTUALLY_ASSEMBLE_MPI_CHUNKS .or. NCHUNKS == 1) then
      iphase_comm = 9999 ! this means everything is finished
      return
    endif

    ! ***************************************************************
    !  transmit messages in forward direction (iprocfrom -> iprocto)
    ! ***************************************************************

    !---- put slices in receive mode
    !---- a given slice can belong to at most two faces

    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if(myrank==iprocfrom_faces(imsg) .or. myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if(myrank==iprocto_faces(imsg)) then
        sender = iprocfrom_faces(imsg)
        npoin2D_chunks = npoin2D_faces(icount_faces)
        call MPI_IRECV(buffer_received_faces_scalar(1,icount_faces), &
                      npoin2D_chunks,CUSTOM_MPI_TYPE,sender, &
                      itag,MPI_COMM_WORLD,request_receive_array(icount_faces),ier)
    !   do ipoin2D=1,npoin2D_chunks
    !     array_val(iboolfaces(ipoin2D,icount_faces)) = &
    !        array_val(iboolfaces(ipoin2D,icount_faces)) + buffer_received_faces_scalar(ipoin2D)
    !   enddo
      endif
    enddo

    !---- put slices in send mode
    !---- a given slice can belong to at most two faces
    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if(myrank==iprocfrom_faces(imsg) .or. myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if(myrank==iprocfrom_faces(imsg)) then
        receiver = iprocto_faces(imsg)
        npoin2D_chunks = npoin2D_faces(icount_faces)
        do ipoin2D=1,npoin2D_chunks
          buffer_send_faces_scalar(ipoin2D,icount_faces) = array_val(iboolfaces(ipoin2D,icount_faces))
        enddo
        call MPI_ISEND(buffer_send_faces_scalar(1,icount_faces),npoin2D_chunks, &
                      CUSTOM_MPI_TYPE,receiver,itag,MPI_COMM_WORLD,request_send_array(icount_faces),ier)
      endif
    enddo

    iphase_comm = iphase_comm + 1
    return ! exit because we have started some communications therefore we need some time

  endif !!!!!!!!! end of iphase_comm 5

  if(iphase_comm == 6) then

    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if(myrank==iprocfrom_faces(imsg) .or. myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if(myrank==iprocto_faces(imsg)) then
        call MPI_TEST(request_receive_array(icount_faces),flag_result_test,msg_status,ier)
        if(.not. flag_result_test) return ! exit if message not received yet
      endif
    enddo

    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if(myrank==iprocfrom_faces(imsg) .or. myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if(myrank==iprocfrom_faces(imsg)) then
        call MPI_TEST(request_send_array(icount_faces),flag_result_test,msg_status,ier)
        if(.not. flag_result_test) return ! exit if message not sent yet
      endif
    enddo

    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if(myrank==iprocfrom_faces(imsg) .or. myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if(myrank==iprocto_faces(imsg)) then
        do ipoin2D=1,npoin2D_faces(icount_faces)
          array_val(iboolfaces(ipoin2D,icount_faces)) = &
            array_val(iboolfaces(ipoin2D,icount_faces)) + buffer_received_faces_scalar(ipoin2D,icount_faces)
        enddo
      endif
    enddo

    ! *********************************************************************
    !  transmit messages back in opposite direction (iprocto -> iprocfrom)
    ! *********************************************************************

    !---- put slices in receive mode
    !---- a given slice can belong to at most two faces

    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if(myrank==iprocfrom_faces(imsg) .or. myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if(myrank==iprocfrom_faces(imsg)) then
        sender = iprocto_faces(imsg)
        npoin2D_chunks = npoin2D_faces(icount_faces)
        call MPI_IRECV(buffer_received_faces_scalar(1,icount_faces), &
                      npoin2D_chunks,CUSTOM_MPI_TYPE,sender, &
                      itag,MPI_COMM_WORLD,request_receive_array(icount_faces),ier)
      !   do ipoin2D=1,npoin2D_chunks
      !     array_val(iboolfaces(ipoin2D,icount_faces)) = buffer_received_faces_scalar(ipoin2D)
      !   enddo
      endif
    enddo

    !---- put slices in send mode
    !---- a given slice can belong to at most two faces
    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if(myrank==iprocfrom_faces(imsg) .or. myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if(myrank==iprocto_faces(imsg)) then
        receiver = iprocfrom_faces(imsg)
        npoin2D_chunks = npoin2D_faces(icount_faces)
        do ipoin2D=1,npoin2D_chunks
          buffer_send_faces_scalar(ipoin2D,icount_faces) = array_val(iboolfaces(ipoin2D,icount_faces))
        enddo
        call MPI_ISEND(buffer_send_faces_scalar(1,icount_faces),npoin2D_chunks, &
                      CUSTOM_MPI_TYPE,receiver,itag,MPI_COMM_WORLD,request_send_array(icount_faces),ier)
      endif
    enddo

    iphase_comm = iphase_comm + 1
    return ! exit because we have started some communications therefore we need some time

  endif !!!!!!!!! end of iphase_comm 6

  if(iphase_comm == 7) then

    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if(myrank==iprocfrom_faces(imsg) .or. myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if(myrank==iprocto_faces(imsg)) then
        call MPI_TEST(request_send_array(icount_faces),flag_result_test,msg_status,ier)
        if(.not. flag_result_test) return ! exit if message not received yet
      endif
    enddo

    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if(myrank==iprocfrom_faces(imsg) .or. myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if(myrank==iprocfrom_faces(imsg)) then
        call MPI_TEST(request_receive_array(icount_faces),flag_result_test,msg_status,ier)
        if(.not. flag_result_test) return ! exit if message not sent yet
      endif
    enddo

    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if(myrank==iprocfrom_faces(imsg) .or. myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if(myrank==iprocfrom_faces(imsg)) then
        do ipoin2D=1,npoin2D_faces(icount_faces)
          array_val(iboolfaces(ipoin2D,icount_faces)) = buffer_received_faces_scalar(ipoin2D,icount_faces)
        enddo
      endif
    enddo

    ! this is the exit condition, to go beyond the last phase number
    iphase_comm = iphase_comm + 1

    !! DK DK do the rest in blocking for now, for simplicity

    !----
    !---- start MPI assembling corners
    !----

    ! scheme for corners cannot deadlock even if NPROC_XI = NPROC_ETA = 1

    ! ***************************************************************
    !  transmit messages in forward direction (two workers -> master)
    ! ***************************************************************

    icount_corners = 0

    do imsg = 1,NCORNERSCHUNKS

      if(myrank == iproc_master_corners(imsg) .or. &
         myrank == iproc_worker1_corners(imsg) .or. &
         (NCHUNKS /= 2 .and. myrank == iproc_worker2_corners(imsg))) icount_corners = icount_corners + 1

      !---- receive messages from the two workers on the master
      if(myrank==iproc_master_corners(imsg)) then

        ! receive from worker #1 and add to local array
        sender = iproc_worker1_corners(imsg)
        call MPI_RECV(buffer_recv_chunkcorn_scalar,NGLOB1D_RADIAL, &
              CUSTOM_MPI_TYPE,sender,itag,MPI_COMM_WORLD,msg_status,ier)
        do ipoin1D=1,NGLOB1D_RADIAL
          array_val(iboolcorner(ipoin1D,icount_corners)) = array_val(iboolcorner(ipoin1D,icount_corners)) + &
                   buffer_recv_chunkcorn_scalar(ipoin1D)
        enddo

        ! receive from worker #2 and add to local array
        if(NCHUNKS /= 2) then
          sender = iproc_worker2_corners(imsg)
          call MPI_RECV(buffer_recv_chunkcorn_scalar,NGLOB1D_RADIAL, &
                      CUSTOM_MPI_TYPE,sender,itag,MPI_COMM_WORLD,msg_status,ier)
          do ipoin1D=1,NGLOB1D_RADIAL
            array_val(iboolcorner(ipoin1D,icount_corners)) = array_val(iboolcorner(ipoin1D,icount_corners)) + &
                   buffer_recv_chunkcorn_scalar(ipoin1D)
          enddo
        endif

      endif

      !---- send messages from the two workers to the master
      if(myrank==iproc_worker1_corners(imsg) .or. &
                  (NCHUNKS /= 2 .and. myrank==iproc_worker2_corners(imsg))) then

        receiver = iproc_master_corners(imsg)
        do ipoin1D=1,NGLOB1D_RADIAL
          buffer_send_chunkcorn_scalar(ipoin1D) = array_val(iboolcorner(ipoin1D,icount_corners))
        enddo
        call MPI_SEND(buffer_send_chunkcorn_scalar,NGLOB1D_RADIAL,CUSTOM_MPI_TYPE, &
                  receiver,itag,MPI_COMM_WORLD,ier)

      endif

      ! *********************************************************************
      !  transmit messages back in opposite direction (master -> two workers)
      ! *********************************************************************

      !---- receive messages from the master on the two workers
      if(myrank==iproc_worker1_corners(imsg) .or. &
                  (NCHUNKS /= 2 .and. myrank==iproc_worker2_corners(imsg))) then

      ! receive from master and copy to local array
        sender = iproc_master_corners(imsg)
        call MPI_RECV(buffer_recv_chunkcorn_scalar,NGLOB1D_RADIAL, &
              CUSTOM_MPI_TYPE,sender,itag,MPI_COMM_WORLD,msg_status,ier)
        do ipoin1D=1,NGLOB1D_RADIAL
          array_val(iboolcorner(ipoin1D,icount_corners)) = buffer_recv_chunkcorn_scalar(ipoin1D)
        enddo

      endif

      !---- send messages from the master to the two workers
      if(myrank==iproc_master_corners(imsg)) then

        do ipoin1D=1,NGLOB1D_RADIAL
          buffer_send_chunkcorn_scalar(ipoin1D) = array_val(iboolcorner(ipoin1D,icount_corners))
        enddo

        ! send to worker #1
        receiver = iproc_worker1_corners(imsg)
        call MPI_SEND(buffer_send_chunkcorn_scalar,NGLOB1D_RADIAL,CUSTOM_MPI_TYPE, &
                  receiver,itag,MPI_COMM_WORLD,ier)

        ! send to worker #2
        if(NCHUNKS /= 2) then
          receiver = iproc_worker2_corners(imsg)
          call MPI_SEND(buffer_send_chunkcorn_scalar,NGLOB1D_RADIAL,CUSTOM_MPI_TYPE, &
                      receiver,itag,MPI_COMM_WORLD,ier)
        endif

      endif

    enddo

  endif !!!!!!!!! end of iphase_comm 7

  end subroutine assemble_MPI_scalar


!
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!

! daniel: TODO - still local versions


  subroutine assemble_MPI_scalar_ext_mesh_s(NPROC,NGLOB_AB,array_val, &
                        buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbours_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)

! non-blocking MPI send

  implicit none

  include "constants.h"

  integer :: NPROC
  integer :: NGLOB_AB
  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

! array to send
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: array_val


  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbours_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh

  integer ipoin,iinterface

! sends only if more than one partition
  if(NPROC > 1) then

    ! partition border copy into the buffer
    do iinterface = 1, num_interfaces_ext_mesh
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        buffer_send_scalar_ext_mesh(ipoin,iinterface) = &
          array_val(ibool_interfaces_ext_mesh(ipoin,iinterface))
      enddo
    enddo

    ! send messages
    do iinterface = 1, num_interfaces_ext_mesh
      ! non-blocking synchronous send request
      call isend_cr(buffer_send_scalar_ext_mesh(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
           nibool_interfaces_ext_mesh(iinterface), &
           my_neighbours_ext_mesh(iinterface), &
           itag, &
           request_send_scalar_ext_mesh(iinterface) &
           )
      ! receive request
      call irecv_cr(buffer_recv_scalar_ext_mesh(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
           nibool_interfaces_ext_mesh(iinterface), &
           my_neighbours_ext_mesh(iinterface), &
           itag, &
           request_recv_scalar_ext_mesh(iinterface) &
           )

    enddo

  endif

  end subroutine assemble_MPI_scalar_ext_mesh_s

!
!----
!

  subroutine assemble_MPI_scalar_ext_mesh_w(NPROC,NGLOB_AB,array_val, &
                        buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh, &
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)

! waits for send/receiver to be completed and assembles contributions

  implicit none

  include "constants.h"

  integer :: NPROC
  integer :: NGLOB_AB
  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh
! array to assemble
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: array_val


  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_scalar_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh

  integer ipoin,iinterface

! assemble only if more than one partition
  if(NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_recv_scalar_ext_mesh(iinterface))
    enddo

    ! adding contributions of neighbours
    do iinterface = 1, num_interfaces_ext_mesh
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        array_val(ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
             array_val(ibool_interfaces_ext_mesh(ipoin,iinterface)) &
             + buffer_recv_scalar_ext_mesh(ipoin,iinterface)
      enddo
    enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_send_scalar_ext_mesh(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_scalar_ext_mesh_w



!
!-------------------------------------------------------------------------------------------------
!

!daniel: TODO - cuda scalar assembly...

  subroutine assemble_MPI_scalar_send_cuda(NPROC, &
                                          buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                                          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                          nibool_interfaces_ext_mesh, &
                                          my_neighbours_ext_mesh, &
                                          request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh, &
                                          FORWARD_OR_ADJOINT)

! non-blocking MPI send

  ! sends data
  ! note: assembling data already filled into buffer_send_scalar_ext_mesh array

  use constants
  use specfem_par,only: Mesh_pointer

  implicit none

  integer :: NPROC
  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbours_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh

  integer :: FORWARD_OR_ADJOINT

  ! local parameters
  integer iinterface

! sends only if more than one partition
  if(NPROC > 1) then

    ! preparation of the contribution between partitions using MPI
    ! transfers mpi buffers to CPU
    call transfer_boun_pot_from_device(Mesh_pointer, &
                                      buffer_send_scalar_ext_mesh, &
                                      FORWARD_OR_ADJOINT)

    ! send messages
    do iinterface = 1, num_interfaces_ext_mesh
      ! non-blocking synchronous send request
      call isend_cr(buffer_send_scalar_ext_mesh(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
           nibool_interfaces_ext_mesh(iinterface), &
           my_neighbours_ext_mesh(iinterface), &
           itag, &
           request_send_scalar_ext_mesh(iinterface) &
           )
      ! receive request
      call irecv_cr(buffer_recv_scalar_ext_mesh(1:nibool_interfaces_ext_mesh(iinterface),iinterface), &
           nibool_interfaces_ext_mesh(iinterface), &
           my_neighbours_ext_mesh(iinterface), &
           itag, &
           request_recv_scalar_ext_mesh(iinterface) &
           )

    enddo

  endif

  end subroutine assemble_MPI_scalar_send_cuda

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_scalar_write_cuda(NPROC, &
                                            buffer_recv_scalar_ext_mesh, &
                                            num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                            request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh, &
                                            FORWARD_OR_ADJOINT)

! waits for send/receiver to be completed and assembles contributions

  use constants
  use specfem_par,only: Mesh_pointer

  implicit none

  integer :: NPROC

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_scalar_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh

  integer :: FORWARD_OR_ADJOINT

  ! local parameters
  integer :: iinterface

! assemble only if more than one partition
  if(NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_recv_scalar_ext_mesh(iinterface))
    enddo

    ! adding contributions of neighbours
    call transfer_asmbl_pot_to_device(Mesh_pointer, &
                                    buffer_recv_scalar_ext_mesh, &
                                    FORWARD_OR_ADJOINT)

    ! note: adding contributions of neighbours has been done just above for cuda
    !do iinterface = 1, num_interfaces_ext_mesh
    !  do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
    !    array_val(ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
    !         array_val(ibool_interfaces_ext_mesh(ipoin,iinterface)) &
    !         + buffer_recv_scalar_ext_mesh(ipoin,iinterface)
    !  enddo
    !enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_send_scalar_ext_mesh(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_scalar_write_cuda
