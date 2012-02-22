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
!---- we handle two regions (crust/mantle and inner core) in the same MPI call
!---- to reduce the total number of MPI calls
!----

  subroutine assemble_MPI_vector(myrank,accel_crust_mantle,accel_inner_core, &
                iproc_xi,iproc_eta,ichunk,addressing, &
                iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,&
                iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
                npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
                iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
                iboolleft_xi_inner_core,iboolright_xi_inner_core,&
                iboolleft_eta_inner_core,iboolright_eta_inner_core, &
                npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
                iboolfaces_inner_core,iboolcorner_inner_core, &
                iprocfrom_faces,iprocto_faces, &
                iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
                buffer_send_faces_vector,buffer_received_faces_vector, &
                npoin2D_max_all_CM_IC, &
                buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector, &
                request_send,request_receive, &
                request_send_array,request_receive_array, &
                NUMMSGS_FACES,NCORNERSCHUNKS, &
                NPROC_XI,NPROC_ETA,NGLOB1D_RADIAL_crust_mantle, &
                NGLOB1D_RADIAL_inner_core,NCHUNKS,iphase_comm)

  implicit none

  ! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  ! include values created by the mesher
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank,NCHUNKS,iphase_comm

  ! the two arrays to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE), intent(inout) :: accel_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE), intent(inout) :: accel_inner_core

  integer, intent(in) :: iproc_xi,iproc_eta,ichunk
  integer, intent(in) :: npoin2D_faces_crust_mantle(NUMFACES_SHARED)
  integer, intent(in) :: npoin2D_faces_inner_core(NUMFACES_SHARED)

  integer, dimension(NB_SQUARE_EDGES_ONEDIR), intent(in) :: npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
        npoin2D_xi_inner_core,npoin2D_eta_inner_core

  integer, intent(in) :: NGLOB1D_RADIAL_crust_mantle,NGLOB1D_RADIAL_inner_core,NPROC_XI,NPROC_ETA
  integer, intent(in) :: NUMMSGS_FACES,NCORNERSCHUNKS

  ! for addressing of the slices
  integer, dimension(NCHUNKS,0:NPROC_XI-1,0:NPROC_ETA-1), intent(in) :: addressing

  ! 2-D addressing and buffers for summation between slices
  integer, dimension(NGLOB2DMAX_XMIN_XMAX_CM), intent(in) :: iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_CM), intent(in) :: iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle
  integer, dimension(NGLOB2DMAX_XMIN_XMAX_IC), intent(in) :: iboolleft_xi_inner_core,iboolright_xi_inner_core
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_IC), intent(in) :: iboolleft_eta_inner_core,iboolright_eta_inner_core

  ! indirect addressing for each corner of the chunks
  integer, dimension(NGLOB1D_RADIAL_crust_mantle,NUMCORNERS_SHARED), intent(in) :: iboolcorner_crust_mantle
  integer, dimension(NGLOB1D_RADIAL_inner_core,NUMCORNERS_SHARED), intent(in) :: iboolcorner_inner_core

  integer, intent(in) :: npoin2D_max_all_CM_IC
  integer, dimension(NGLOB2DMAX_XY_VAL,NUMFACES_SHARED), intent(in) :: iboolfaces_crust_mantle
  integer, dimension(NGLOB2DMAX_XY_VAL,NUMFACES_SHARED), intent(in) :: iboolfaces_inner_core
  real(kind=CUSTOM_REAL), dimension(NDIM,npoin2D_max_all_CM_IC,NUMFACES_SHARED), intent(inout) :: &
      buffer_send_faces_vector,buffer_received_faces_vector

  ! buffers for send and receive between corners of the chunks
  ! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB1D_RADIAL_crust_mantle + NGLOB1D_RADIAL_inner_core), intent(inout) :: &
    buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector

  ! stored as global parameters
  integer, intent(inout) :: request_send,request_receive
  integer, dimension(NUMFACES_SHARED), intent(inout) :: request_send_array,request_receive_array


  ! ---- arrays to assemble between chunks
  ! communication pattern for faces between chunks
  integer, dimension(NUMMSGS_FACES), intent(in) :: iprocfrom_faces,iprocto_faces
  ! communication pattern for corners between chunks
  integer, dimension(NCORNERSCHUNKS), intent(in) :: iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners

! local parameters

  integer :: icount_corners

  ! MPI status of messages to be received
  integer, dimension(MPI_STATUS_SIZE) :: msg_status

  integer :: ipoin,ipoin2D,ipoin1D
  integer :: sender,receiver
  integer :: imsg
  integer :: icount_faces,npoin2D_chunks_all

  integer :: NGLOB1D_RADIAL_all
  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_all,npoin2D_eta_all

  logical :: flag_result_test
  integer :: ier

  ! daniel: TODO - comment below might be obsolete..
  ! do not remove the "save" statement because this routine is non blocking
  ! therefore it needs to find the right value of ioffset when it re-enters
  ! the routine later to perform the next communication step
  integer, save :: ioffset


  ! check flag to see if we need to assemble (might be turned off when debugging)
  if (.not. ACTUALLY_ASSEMBLE_MPI_SLICES) then
    iphase_comm = 9999 ! this means everything is finished
    return
  endif

  ! here we have to assemble all the contributions between slices using MPI

  ! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
  npoin2D_xi_all(:) = npoin2D_xi_crust_mantle(:) + npoin2D_xi_inner_core(:)
  npoin2D_eta_all(:) = npoin2D_eta_crust_mantle(:) + npoin2D_eta_inner_core(:)

  !----
  !---- assemble the contributions between slices using MPI
  !----

  !----
  !---- first assemble along xi using the 2-D topology
  !----
  select case( iphase_comm )

  case( 1 )

    ! slices send out values along xi sides (right face) forward along each row

    ! slices copy the right face into the buffer
    do ipoin = 1,npoin2D_xi_crust_mantle(2)
      buffer_send_faces_vector(:,ipoin,1) = accel_crust_mantle(:,iboolright_xi_crust_mantle(ipoin))
    enddo

    ! the buffer for the inner core starts right after the buffer for the crust and mantle
    ioffset = npoin2D_xi_crust_mantle(2)

    do ipoin = 1,npoin2D_xi_inner_core(2)
      buffer_send_faces_vector(:,ioffset + ipoin,1) = accel_inner_core(:,iboolright_xi_inner_core(ipoin))
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
    ! requests to receive message
    call MPI_IRECV(buffer_received_faces_vector,NDIM*npoin2D_xi_all(1),CUSTOM_MPI_TYPE,sender, &
          itag,MPI_COMM_WORLD,request_receive,ier)

    ! sends out buffer
    call MPI_ISEND(buffer_send_faces_vector,NDIM*npoin2D_xi_all(2),CUSTOM_MPI_TYPE,receiver, &
          itag2,MPI_COMM_WORLD,request_send,ier)

    iphase_comm = iphase_comm + 1
    return ! exit because we have started some communications therefore we need some time


  case( 2 )

    ! slices assemble along (left face) xi sides
    ! and send out values along xi sides (left face) backward along each row

    ! checks if messages have been sent out and requested ones received
    ! call MPI_WAIT(request_send,msg_status,ier)
    ! call MPI_WAIT(request_receive,msg_status,ier)
    call MPI_TEST(request_send,flag_result_test,msg_status,ier)
    if(.not. flag_result_test) return ! exit if message not sent yet
    call MPI_TEST(request_receive,flag_result_test,msg_status,ier)
    if(.not. flag_result_test) return ! exit if message not received yet

    ! all slices add the buffer received to the contributions on the left face
    if(iproc_xi > 0) then

      do ipoin = 1,npoin2D_xi_crust_mantle(1)
        accel_crust_mantle(:,iboolleft_xi_crust_mantle(ipoin)) = &
                accel_crust_mantle(:,iboolleft_xi_crust_mantle(ipoin)) + &
                                  buffer_received_faces_vector(:,ipoin,1)
      enddo

      ! the buffer for the inner core starts right after the buffer for the crust and mantle
      ioffset = npoin2D_xi_crust_mantle(1)

      do ipoin = 1,npoin2D_xi_inner_core(1)
        accel_inner_core(:,iboolleft_xi_inner_core(ipoin)) = &
                  accel_inner_core(:,iboolleft_xi_inner_core(ipoin)) + &
                                  buffer_received_faces_vector(:,ioffset + ipoin,1)
      enddo

    endif

    ! the contributions are correctly assembled on the left side of each slice
    ! now we have to send the result back to the sender
    ! all slices copy the left face into the buffer
    do ipoin = 1,npoin2D_xi_crust_mantle(1)
      buffer_send_faces_vector(:,ipoin,1) = accel_crust_mantle(:,iboolleft_xi_crust_mantle(ipoin))
    enddo

    ! the buffer for the inner core starts right after the buffer for the crust and mantle
    ioffset = npoin2D_xi_crust_mantle(1)

    do ipoin = 1,npoin2D_xi_inner_core(1)
      buffer_send_faces_vector(:,ioffset + ipoin,1) = accel_inner_core(:,iboolleft_xi_inner_core(ipoin))
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
    call MPI_IRECV(buffer_received_faces_vector,NDIM*npoin2D_xi_all(2),CUSTOM_MPI_TYPE,sender, &
          itag,MPI_COMM_WORLD,request_receive,ier)

    call MPI_ISEND(buffer_send_faces_vector,NDIM*npoin2D_xi_all(1),CUSTOM_MPI_TYPE,receiver, &
          itag2,MPI_COMM_WORLD,request_send,ier)

    iphase_comm = iphase_comm + 1
    return ! exit because we have started some communications therefore we need some time

  case( 3 )

    ! slices set (right face) xi sides
    ! and sent out values along eta sides (right face) forward along each column

    ! checks if messages have been sent out and received
    ! call MPI_WAIT(request_send,msg_status,ier)
    ! call MPI_WAIT(request_receive,msg_status,ier)
    call MPI_TEST(request_send,flag_result_test,msg_status,ier)
    if(.not. flag_result_test) return ! exit if message not sent yet
    call MPI_TEST(request_receive,flag_result_test,msg_status,ier)
    if(.not. flag_result_test) return ! exit if message not received yet

    ! all slices copy the buffer received to the contributions on the right face
    if(iproc_xi < NPROC_XI-1) then

      do ipoin = 1,npoin2D_xi_crust_mantle(2)
        accel_crust_mantle(:,iboolright_xi_crust_mantle(ipoin)) = buffer_received_faces_vector(:,ipoin,1)
      enddo

      ! the buffer for the inner core starts right after the buffer for the crust and mantle
      ioffset = npoin2D_xi_crust_mantle(2)

      do ipoin = 1,npoin2D_xi_inner_core(2)
        accel_inner_core(:,iboolright_xi_inner_core(ipoin)) = buffer_received_faces_vector(:,ioffset + ipoin,1)
      enddo

    endif

    !----
    !---- then assemble along eta using the 2-D topology
    !----

    ! slices copy the right face into the buffer
    do ipoin = 1,npoin2D_eta_crust_mantle(2)
      buffer_send_faces_vector(:,ipoin,1) = accel_crust_mantle(:,iboolright_eta_crust_mantle(ipoin))
    enddo

    ! the buffer for the inner core starts right after the buffer for the crust and mantle
    ioffset = npoin2D_eta_crust_mantle(2)

    do ipoin = 1,npoin2D_eta_inner_core(2)
      buffer_send_faces_vector(:,ioffset + ipoin,1) = accel_inner_core(:,iboolright_eta_inner_core(ipoin))
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
    call MPI_IRECV(buffer_received_faces_vector,NDIM*npoin2D_eta_all(1),CUSTOM_MPI_TYPE,sender, &
      itag,MPI_COMM_WORLD,request_receive,ier)

    call MPI_ISEND(buffer_send_faces_vector,NDIM*npoin2D_eta_all(2),CUSTOM_MPI_TYPE,receiver, &
      itag2,MPI_COMM_WORLD,request_send,ier)

    iphase_comm = iphase_comm + 1
    return ! exit because we have started some communications therefore we need some time


 case( 4 )

    ! slices assemble along (left face) eta sides
    ! and sent out values along eta sides (left face) backward along each column

    ! checks if messages have been sent out and received
    ! call MPI_WAIT(request_send,msg_status,ier)
    ! call MPI_WAIT(request_receive,msg_status,ier)
    call MPI_TEST(request_send,flag_result_test,msg_status,ier)
    if(.not. flag_result_test) return ! exit if message not sent yet
    call MPI_TEST(request_receive,flag_result_test,msg_status,ier)
    if(.not. flag_result_test) return ! exit if message not received yet

    ! all slices add the buffer received to the contributions on the left face
    if(iproc_eta > 0) then

      do ipoin = 1,npoin2D_eta_crust_mantle(1)
        accel_crust_mantle(:,iboolleft_eta_crust_mantle(ipoin)) = &
                  accel_crust_mantle(:,iboolleft_eta_crust_mantle(ipoin)) + &
                                  buffer_received_faces_vector(:,ipoin,1)
      enddo

      ! the buffer for the inner core starts right after the buffer for the crust and mantle
      ioffset = npoin2D_eta_crust_mantle(1)

      do ipoin = 1,npoin2D_eta_inner_core(1)
        accel_inner_core(:,iboolleft_eta_inner_core(ipoin)) = &
                  accel_inner_core(:,iboolleft_eta_inner_core(ipoin)) + &
                                  buffer_received_faces_vector(:,ioffset + ipoin,1)
      enddo

    endif

    ! the contributions are correctly assembled on the left side of each slice
    ! now we have to send the result back to the sender
    ! all slices copy the left face into the buffer
    do ipoin = 1,npoin2D_eta_crust_mantle(1)
      buffer_send_faces_vector(:,ipoin,1) = accel_crust_mantle(:,iboolleft_eta_crust_mantle(ipoin))
    enddo

    ! the buffer for the inner core starts right after the buffer for the crust and mantle
    ioffset = npoin2D_eta_crust_mantle(1)

    do ipoin = 1,npoin2D_eta_inner_core(1)
      buffer_send_faces_vector(:,ioffset + ipoin,1) = accel_inner_core(:,iboolleft_eta_inner_core(ipoin))
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
    call MPI_IRECV(buffer_received_faces_vector,NDIM*npoin2D_eta_all(2),CUSTOM_MPI_TYPE,sender, &
      itag,MPI_COMM_WORLD,request_receive,ier)

    call MPI_ISEND(buffer_send_faces_vector,NDIM*npoin2D_eta_all(1),CUSTOM_MPI_TYPE,receiver, &
      itag2,MPI_COMM_WORLD,request_send,ier)

    iphase_comm = iphase_comm + 1
    return ! exit because we have started some communications therefore we need some time


  case( 5 )

    ! slices set (right face) eta sides
    ! and sent out values for neighbor chunks (iboolfaces)

    ! checks if messages have been sent out and received
    ! call MPI_WAIT(request_send,msg_status,ier)
    ! call MPI_WAIT(request_receive,msg_status,ier)
    call MPI_TEST(request_send,flag_result_test,msg_status,ier)
    if(.not. flag_result_test) return ! exit if message not sent yet
    call MPI_TEST(request_receive,flag_result_test,msg_status,ier)
    if(.not. flag_result_test) return ! exit if message not received yet

    ! all slices copy the buffer received to the contributions on the right face
    if(iproc_eta < NPROC_ETA-1) then

      do ipoin = 1,npoin2D_eta_crust_mantle(2)
        accel_crust_mantle(:,iboolright_eta_crust_mantle(ipoin)) = buffer_received_faces_vector(:,ipoin,1)
      enddo

      ! the buffer for the inner core starts right after the buffer for the crust and mantle
      ioffset = npoin2D_eta_crust_mantle(2)

      do ipoin = 1,npoin2D_eta_inner_core(2)
        accel_inner_core(:,iboolright_eta_inner_core(ipoin)) = buffer_received_faces_vector(:,ioffset + ipoin,1)
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

        ! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
        npoin2D_chunks_all = npoin2D_faces_crust_mantle(icount_faces) + npoin2D_faces_inner_core(icount_faces)

        call MPI_IRECV(buffer_received_faces_vector(1,1,icount_faces),NDIM*npoin2D_chunks_all,CUSTOM_MPI_TYPE,sender, &
                itag,MPI_COMM_WORLD,request_receive_array(icount_faces),ier)

        !   do ipoin2D = 1,npoin2D_faces_crust_mantle(icount_faces)
        !     accel_crust_mantle(1,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = &
        !        accel_crust_mantle(1,iboolfaces_crust_mantle(ipoin2D,icount_faces)) &
        !           + buffer_received_faces_vector(1,ipoin2D,icount_faces)
        !     accel_crust_mantle(2,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = &
        !        accel_crust_mantle(2,iboolfaces_crust_mantle(ipoin2D,icount_faces)) &
        !           + buffer_received_faces_vector(2,ipoin2D,icount_faces)
        !     accel_crust_mantle(3,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = &
        !        accel_crust_mantle(3,iboolfaces_crust_mantle(ipoin2D,icount_faces)) &
        !           + buffer_received_faces_vector(3,ipoin2D,icount_faces)
        !   enddo

        ! the buffer for the inner core starts right after the buffer for the crust and mantle
        !   ioffset = npoin2D_faces_crust_mantle(icount_faces)

        !   do ipoin2D = 1,npoin2D_faces_inner_core(icount_faces)
        !     accel_inner_core(1,iboolfaces_inner_core(ipoin2D,icount_faces)) = &
        !        accel_inner_core(1,iboolfaces_inner_core(ipoin2D,icount_faces)) + &
        !          buffer_received_faces_vector(1,ioffset + ipoin2D,icount_faces)
        !     accel_inner_core(2,iboolfaces_inner_core(ipoin2D,icount_faces)) = &
        !        accel_inner_core(2,iboolfaces_inner_core(ipoin2D,icount_faces)) + &
        !          buffer_received_faces_vector(2,ioffset + ipoin2D,icount_faces)
        !     accel_inner_core(3,iboolfaces_inner_core(ipoin2D,icount_faces)) = &
        !        accel_inner_core(3,iboolfaces_inner_core(ipoin2D,icount_faces)) + &
        !          buffer_received_faces_vector(3,ioffset + ipoin2D,icount_faces)
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

        ! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
        npoin2D_chunks_all = npoin2D_faces_crust_mantle(icount_faces) + npoin2D_faces_inner_core(icount_faces)

        do ipoin2D = 1,npoin2D_faces_crust_mantle(icount_faces)
          buffer_send_faces_vector(:,ipoin2D,icount_faces) = &
            accel_crust_mantle(:,iboolfaces_crust_mantle(ipoin2D,icount_faces))
        enddo

        ! the buffer for the inner core starts right after the buffer for the crust and mantle
        ioffset = npoin2D_faces_crust_mantle(icount_faces)

        do ipoin2D = 1,npoin2D_faces_inner_core(icount_faces)
          buffer_send_faces_vector(:,ioffset + ipoin2D,icount_faces) = &
            accel_inner_core(:,iboolfaces_inner_core(ipoin2D,icount_faces))
        enddo

        call MPI_ISEND(buffer_send_faces_vector(1,1,icount_faces),NDIM*npoin2D_chunks_all,CUSTOM_MPI_TYPE,receiver,itag, &
                         MPI_COMM_WORLD,request_send_array(icount_faces),ier)
      endif
    enddo

    iphase_comm = iphase_comm + 1
    return ! exit because we have started some communications therefore we need some time

  case( 6 )

    ! receiver slices on chunk faces assemble values (iboolfaces)
    ! and send values back to senders

    ! checks if messages have been received
    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if(myrank==iprocfrom_faces(imsg) .or. myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if(myrank==iprocto_faces(imsg)) then
        call MPI_TEST(request_receive_array(icount_faces),flag_result_test,msg_status,ier)
        if(.not. flag_result_test) return ! exit if message not received yet
      endif
    enddo

    ! checks if messages have been sent out
    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if(myrank==iprocfrom_faces(imsg) .or. myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if(myrank==iprocfrom_faces(imsg)) then
        call MPI_TEST(request_send_array(icount_faces),flag_result_test,msg_status,ier)
        if(.not. flag_result_test) return ! exit if message not sent yet
      endif
    enddo

    ! assembles values on chunk faces
    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if(myrank==iprocfrom_faces(imsg) .or. myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if(myrank==iprocto_faces(imsg)) then

        do ipoin2D = 1,npoin2D_faces_crust_mantle(icount_faces)
          accel_crust_mantle(:,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = &
             accel_crust_mantle(:,iboolfaces_crust_mantle(ipoin2D,icount_faces)) &
             + buffer_received_faces_vector(:,ipoin2D,icount_faces)
        enddo

        ! the buffer for the inner core starts right after the buffer for the crust and mantle
        ioffset = npoin2D_faces_crust_mantle(icount_faces)

        do ipoin2D = 1,npoin2D_faces_inner_core(icount_faces)
          accel_inner_core(:,iboolfaces_inner_core(ipoin2D,icount_faces)) = &
             accel_inner_core(:,iboolfaces_inner_core(ipoin2D,icount_faces)) + &
               buffer_received_faces_vector(:,ioffset + ipoin2D,icount_faces)
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

        ! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
        npoin2D_chunks_all = npoin2D_faces_crust_mantle(icount_faces) + npoin2D_faces_inner_core(icount_faces)

        call MPI_IRECV(buffer_received_faces_vector(1,1,icount_faces),NDIM*npoin2D_chunks_all,CUSTOM_MPI_TYPE,sender, &
                  itag,MPI_COMM_WORLD,request_receive_array(icount_faces),ier)

        !   do ipoin2D = 1,npoin2D_faces_crust_mantle(icount_faces)
        !     accel_crust_mantle(1,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = &
        !       buffer_received_faces_vector(1,ipoin2D,icount_faces)
        !     accel_crust_mantle(2,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = &
        !       buffer_received_faces_vector(2,ipoin2D,icount_faces)
        !     accel_crust_mantle(3,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = &
        !       buffer_received_faces_vector(3,ipoin2D,icount_faces)
        !   enddo

        ! the buffer for the inner core starts right after the buffer for the crust and mantle
        !   ioffset = npoin2D_faces_crust_mantle(icount_faces)

        !   do ipoin2D = 1,npoin2D_faces_inner_core(icount_faces)
        !     accel_inner_core(1,iboolfaces_inner_core(ipoin2D,icount_faces)) = &
        !       buffer_received_faces_vector(1,ioffset + ipoin2D,icount_faces)
        !     accel_inner_core(2,iboolfaces_inner_core(ipoin2D,icount_faces)) = &
        !       buffer_received_faces_vector(2,ioffset + ipoin2D,icount_faces)
        !     accel_inner_core(3,iboolfaces_inner_core(ipoin2D,icount_faces)) = &
        !       buffer_received_faces_vector(3,ioffset + ipoin2D,icount_faces)
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

        ! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
        npoin2D_chunks_all = npoin2D_faces_crust_mantle(icount_faces) + npoin2D_faces_inner_core(icount_faces)

        do ipoin2D = 1,npoin2D_faces_crust_mantle(icount_faces)
          buffer_send_faces_vector(:,ipoin2D,icount_faces) = &
            accel_crust_mantle(:,iboolfaces_crust_mantle(ipoin2D,icount_faces))
        enddo

        ! the buffer for the inner core starts right after the buffer for the crust and mantle
        ioffset = npoin2D_faces_crust_mantle(icount_faces)

        do ipoin2D = 1,npoin2D_faces_inner_core(icount_faces)
          buffer_send_faces_vector(:,ioffset + ipoin2D,icount_faces) = &
            accel_inner_core(:,iboolfaces_inner_core(ipoin2D,icount_faces))
        enddo

        call MPI_ISEND(buffer_send_faces_vector(1,1,icount_faces),NDIM*npoin2D_chunks_all, &
                        CUSTOM_MPI_TYPE,receiver,itag, &
                         MPI_COMM_WORLD,request_send_array(icount_faces),ier)
      endif
    enddo

    iphase_comm = iphase_comm + 1
    return ! exit because we have started some communications therefore we need some time

  case( 7 )

    ! sender slices on chunk faces set values (iboolfaces)
    ! and slices on corners assemble values on corners with other chunks

    ! checks if messages have been sent
    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if(myrank==iprocfrom_faces(imsg) .or. myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if(myrank==iprocto_faces(imsg)) then
        call MPI_TEST(request_send_array(icount_faces),flag_result_test,msg_status,ier)
        if(.not. flag_result_test) return ! exit if message not received yet
      endif
    enddo

    ! checks if messages have been received
    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if(myrank==iprocfrom_faces(imsg) .or. myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if(myrank==iprocfrom_faces(imsg)) then
        call MPI_TEST(request_receive_array(icount_faces),flag_result_test,msg_status,ier)
        if(.not. flag_result_test) return ! exit if message not sent yet
      endif
    enddo

    ! sets values on faces (iboolfaces)
    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if(myrank==iprocfrom_faces(imsg) .or. myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if(myrank==iprocfrom_faces(imsg)) then
        do ipoin2D = 1,npoin2D_faces_crust_mantle(icount_faces)
          accel_crust_mantle(:,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = &
            buffer_received_faces_vector(:,ipoin2D,icount_faces)
        enddo

        ! the buffer for the inner core starts right after the buffer for the crust and mantle
        ioffset = npoin2D_faces_crust_mantle(icount_faces)

        do ipoin2D = 1,npoin2D_faces_inner_core(icount_faces)
          accel_inner_core(:,iboolfaces_inner_core(ipoin2D,icount_faces)) = &
            buffer_received_faces_vector(:,ioffset + ipoin2D,icount_faces)
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

    ! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
    NGLOB1D_RADIAL_all = NGLOB1D_RADIAL_crust_mantle + NGLOB1D_RADIAL_inner_core

    ! the buffer for the inner core starts right after the buffer for the crust and mantle
    ioffset = NGLOB1D_RADIAL_crust_mantle

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

        call MPI_RECV(buffer_recv_chunkcorn_vector,NDIM*NGLOB1D_RADIAL_all, &
              CUSTOM_MPI_TYPE,sender,itag,MPI_COMM_WORLD,msg_status,ier)

        do ipoin1D = 1,NGLOB1D_RADIAL_crust_mantle
          accel_crust_mantle(:,iboolcorner_crust_mantle(ipoin1D,icount_corners)) = &
                   accel_crust_mantle(:,iboolcorner_crust_mantle(ipoin1D,icount_corners)) + &
                   buffer_recv_chunkcorn_vector(:,ipoin1D)
        enddo

        do ipoin1D = 1,NGLOB1D_RADIAL_inner_core
          accel_inner_core(:,iboolcorner_inner_core(ipoin1D,icount_corners)) = &
                   accel_inner_core(:,iboolcorner_inner_core(ipoin1D,icount_corners)) + &
                   buffer_recv_chunkcorn_vector(:,ioffset + ipoin1D)
        enddo

        ! receive from worker #2 and add to local array
        if(NCHUNKS /= 2) then

          sender = iproc_worker2_corners(imsg)

          call MPI_RECV(buffer_recv_chunkcorn_vector,NDIM*NGLOB1D_RADIAL_all, &
                CUSTOM_MPI_TYPE,sender,itag,MPI_COMM_WORLD,msg_status,ier)

          do ipoin1D = 1,NGLOB1D_RADIAL_crust_mantle
            accel_crust_mantle(:,iboolcorner_crust_mantle(ipoin1D,icount_corners)) = &
                     accel_crust_mantle(:,iboolcorner_crust_mantle(ipoin1D,icount_corners)) + &
                     buffer_recv_chunkcorn_vector(:,ipoin1D)
          enddo

          do ipoin1D = 1,NGLOB1D_RADIAL_inner_core
            accel_inner_core(:,iboolcorner_inner_core(ipoin1D,icount_corners)) = &
                     accel_inner_core(:,iboolcorner_inner_core(ipoin1D,icount_corners)) + &
                     buffer_recv_chunkcorn_vector(:,ioffset + ipoin1D)
          enddo

        endif

      endif

      !---- send messages from the two workers to the master
      if(myrank==iproc_worker1_corners(imsg) .or. &
                  (NCHUNKS /= 2 .and. myrank==iproc_worker2_corners(imsg))) then

        receiver = iproc_master_corners(imsg)

        do ipoin1D = 1,NGLOB1D_RADIAL_crust_mantle
          buffer_send_chunkcorn_vector(:,ipoin1D) = &
            accel_crust_mantle(:,iboolcorner_crust_mantle(ipoin1D,icount_corners))
        enddo

        do ipoin1D = 1,NGLOB1D_RADIAL_inner_core
          buffer_send_chunkcorn_vector(:,ioffset + ipoin1D) = &
            accel_inner_core(:,iboolcorner_inner_core(ipoin1D,icount_corners))
        enddo

        call MPI_SEND(buffer_send_chunkcorn_vector,NDIM*NGLOB1D_RADIAL_all, &
                      CUSTOM_MPI_TYPE,receiver,itag,MPI_COMM_WORLD,ier)

      endif

      ! *********************************************************************
      !  transmit messages back in opposite direction (master -> two workers)
      ! *********************************************************************

      !---- receive messages from the master on the two workers
      if(myrank==iproc_worker1_corners(imsg) .or. &
                  (NCHUNKS /= 2 .and. myrank==iproc_worker2_corners(imsg))) then

        ! receive from master and copy to local array
        sender = iproc_master_corners(imsg)

        call MPI_RECV(buffer_recv_chunkcorn_vector,NDIM*NGLOB1D_RADIAL_all, &
              CUSTOM_MPI_TYPE,sender,itag,MPI_COMM_WORLD,msg_status,ier)

        do ipoin1D = 1,NGLOB1D_RADIAL_crust_mantle
          accel_crust_mantle(:,iboolcorner_crust_mantle(ipoin1D,icount_corners)) = &
            buffer_recv_chunkcorn_vector(:,ipoin1D)
        enddo

        do ipoin1D = 1,NGLOB1D_RADIAL_inner_core
          accel_inner_core(:,iboolcorner_inner_core(ipoin1D,icount_corners)) = &
            buffer_recv_chunkcorn_vector(:,ioffset + ipoin1D)
        enddo

      endif

      !---- send messages from the master to the two workers
      if(myrank==iproc_master_corners(imsg)) then

        do ipoin1D = 1,NGLOB1D_RADIAL_crust_mantle
          buffer_send_chunkcorn_vector(:,ipoin1D) = &
            accel_crust_mantle(:,iboolcorner_crust_mantle(ipoin1D,icount_corners))
        enddo

        do ipoin1D = 1,NGLOB1D_RADIAL_inner_core
          buffer_send_chunkcorn_vector(:,ioffset + ipoin1D) = &
            accel_inner_core(:,iboolcorner_inner_core(ipoin1D,icount_corners))
        enddo

        ! send to worker #1
        receiver = iproc_worker1_corners(imsg)
        call MPI_SEND(buffer_send_chunkcorn_vector,NDIM*NGLOB1D_RADIAL_all,&
                    CUSTOM_MPI_TYPE,receiver,itag,MPI_COMM_WORLD,ier)

        ! send to worker #2
        if(NCHUNKS /= 2) then
          receiver = iproc_worker2_corners(imsg)
          call MPI_SEND(buffer_send_chunkcorn_vector,NDIM*NGLOB1D_RADIAL_all,&
                      CUSTOM_MPI_TYPE,receiver,itag,MPI_COMM_WORLD,ier)
        endif

      endif

    enddo

  end select

  end subroutine assemble_MPI_vector

!
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!
!

!daniel: todo - still local versions

  subroutine assemble_MPI_vector_ext_mesh_s(NPROC,NGLOB_AB, &
                                           array_val, &
                                           buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                                           num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                           nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                           my_neighbours_ext_mesh, &
                                           request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

  ! sends data

  implicit none

  include "constants.h"

  integer :: NPROC
  integer :: NGLOB_AB

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: array_val

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbours_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: request_send_vector_ext_mesh,request_recv_vector_ext_mesh

  integer ipoin,iinterface

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if(NPROC > 1) then

     ! partition border copy into the buffer
     do iinterface = 1, num_interfaces_ext_mesh
        do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
           buffer_send_vector_ext_mesh(:,ipoin,iinterface) = &
                array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface))
        enddo
     enddo

     ! send messages
     do iinterface = 1, num_interfaces_ext_mesh
        call isend_cr(buffer_send_vector_ext_mesh(1,1,iinterface), &
             NDIM*nibool_interfaces_ext_mesh(iinterface), &
             my_neighbours_ext_mesh(iinterface), &
             itag, &
             request_send_vector_ext_mesh(iinterface) &
             )
        call irecv_cr(buffer_recv_vector_ext_mesh(1,1,iinterface), &
             NDIM*nibool_interfaces_ext_mesh(iinterface), &
             my_neighbours_ext_mesh(iinterface), &
             itag, &
             request_recv_vector_ext_mesh(iinterface) &
             )
     enddo

  endif

  end subroutine assemble_MPI_vector_ext_mesh_s

!
!-------------------------------------------------------------------------------------------------
!

!daniel: TODO - cuda scalar assembly...

! interrupt might improve MPI performance
! see: https://computing.llnl.gov/tutorials/mpi_performance/#Sender-ReceiverSync
!
! check: MP_CSS_INTERRUPT environment variable on IBM systems


  subroutine assemble_MPI_vector_send_cuda(NPROC, &
                                          buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh, &
                                          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                          nibool_interfaces_ext_mesh, &
                                          my_neighbours_ext_mesh, &
                                          request_send_vector_ext_mesh,request_recv_vector_ext_mesh,&
                                          IREGION,FORWARD_OR_ADJOINT)

  ! sends data
  ! note: array to assemble already filled into buffer_send_vector_ext_mesh array
  use constants
  use specfem_par,only: Mesh_pointer

  implicit none

  integer :: NPROC

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_send_vector_ext_mesh,buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh,my_neighbours_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: request_send_vector_ext_mesh,request_recv_vector_ext_mesh

  integer :: IREGION
  integer :: FORWARD_OR_ADJOINT

  ! local parameters
  integer iinterface

  ! send only if more than one partition
  if(NPROC > 1) then

    ! preparation of the contribution between partitions using MPI
    ! transfers mpi buffers to CPU
    call transfer_boun_accel_from_device(Mesh_pointer, &
                                        buffer_send_vector_ext_mesh,&
                                        IREGION,FORWARD_OR_ADJOINT)

     ! send messages
     do iinterface = 1, num_interfaces_ext_mesh
        call isend_cr(buffer_send_vector_ext_mesh(1,1,iinterface), &
             NDIM*nibool_interfaces_ext_mesh(iinterface), &
             my_neighbours_ext_mesh(iinterface), &
             itag, &
             request_send_vector_ext_mesh(iinterface) &
             )
        call irecv_cr(buffer_recv_vector_ext_mesh(1,1,iinterface), &
             NDIM*nibool_interfaces_ext_mesh(iinterface), &
             my_neighbours_ext_mesh(iinterface), &
             itag, &
             request_recv_vector_ext_mesh(iinterface) &
             )
     enddo

  endif

  end subroutine assemble_MPI_vector_send_cuda

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_vector_ext_mesh_w(NPROC,NGLOB_AB, &
                                           array_val, &
                                           buffer_recv_vector_ext_mesh, &
                                           num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                           nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                                           request_send_vector_ext_mesh,request_recv_vector_ext_mesh)

! waits for data to receive and assembles

  implicit none

  include "constants.h"

  integer :: NPROC
  integer :: NGLOB_AB

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: array_val

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_vector_ext_mesh

  integer, dimension(num_interfaces_ext_mesh) :: nibool_interfaces_ext_mesh
  integer, dimension(max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: ibool_interfaces_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: request_send_vector_ext_mesh,request_recv_vector_ext_mesh

  integer ipoin,iinterface

! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if(NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_recv_vector_ext_mesh(iinterface))
    enddo

    ! adding contributions of neighbours
    do iinterface = 1, num_interfaces_ext_mesh
      do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
        array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
             array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) &
             + buffer_recv_vector_ext_mesh(:,ipoin,iinterface)
      enddo
    enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_send_vector_ext_mesh(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_vector_ext_mesh_w


!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_vector_write_cuda(NPROC, &
                                            buffer_recv_vector_ext_mesh, &
                                            num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                                            request_send_vector_ext_mesh,request_recv_vector_ext_mesh, &
                                            IREGION,FORWARD_OR_ADJOINT )

! waits for data to receive and assembles
  use constants
  use specfem_par,only: Mesh_pointer

  implicit none

  integer :: NPROC

  integer :: num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces_ext_mesh,num_interfaces_ext_mesh) :: &
       buffer_recv_vector_ext_mesh
  integer, dimension(num_interfaces_ext_mesh) :: request_send_vector_ext_mesh,request_recv_vector_ext_mesh

  integer :: IREGION
  integer :: FORWARD_OR_ADJOINT

  ! local parameters

  integer iinterface

! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if(NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_recv_vector_ext_mesh(iinterface))
    enddo

    ! adding contributions of neighbours
    call transfer_asmbl_accel_to_device(Mesh_pointer, &
                                      buffer_recv_vector_ext_mesh, &
                                      IREGION,FORWARD_OR_ADJOINT)

    ! This step is done via previous function transfer_and_assemble...
    ! do iinterface = 1, num_interfaces_ext_mesh
    !   do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
    !     array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
    !          array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_vector_ext_mesh(:,ipoin,iinterface)
    !   enddo
    ! enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces_ext_mesh
      call wait_req(request_send_vector_ext_mesh(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_vector_write_cuda

!
!----
!

  subroutine irecv_cr(recvbuf, recvcount, dest, recvtag, req)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer recvcount, dest, recvtag, req
  real(kind=CUSTOM_REAL), dimension(recvcount) :: recvbuf

  integer ier

  call MPI_IRECV(recvbuf(1),recvcount,CUSTOM_MPI_TYPE,dest,recvtag, &
                  MPI_COMM_WORLD,req,ier)

  end subroutine irecv_cr

!
!----
!

  subroutine isend_cr(sendbuf, sendcount, dest, sendtag, req)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer sendcount, dest, sendtag, req
  real(kind=CUSTOM_REAL), dimension(sendcount) :: sendbuf

  integer ier

  call MPI_ISEND(sendbuf(1),sendcount,CUSTOM_MPI_TYPE,dest,sendtag, &
                  MPI_COMM_WORLD,req,ier)

  end subroutine isend_cr

!
!----
!

  subroutine wait_req(req)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer :: req

  integer, dimension(MPI_STATUS_SIZE) :: req_mpi_status

  integer :: ier

  call mpi_wait(req,req_mpi_status,ier)

  end subroutine wait_req

