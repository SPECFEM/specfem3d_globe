!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2011
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
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_vector,buffer_received_faces_vector,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorners_vector,buffer_recv_chunkcorners_vector, &
            NUMMSGS_FACES,NCORNERSCHUNKS, &
            NPROC_XI,NPROC_ETA,NGLOB1D_RADIAL_crust_mantle, &
            NGLOB1D_RADIAL_inner_core,NCHUNKS,iphase)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

! include values created by the mesher
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank,NCHUNKS,iphase

! the two arrays to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: accel_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: accel_inner_core

  integer iproc_xi,iproc_eta,ichunk
  integer npoin2D_faces_crust_mantle(NUMFACES_SHARED)
  integer npoin2D_faces_inner_core(NUMFACES_SHARED)

  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
        npoin2D_xi_inner_core,npoin2D_eta_inner_core

  integer NGLOB1D_RADIAL_crust_mantle,NGLOB1D_RADIAL_inner_core,NPROC_XI,NPROC_ETA
  integer NUMMSGS_FACES,NCORNERSCHUNKS

! for addressing of the slices
  integer, dimension(NCHUNKS,0:NPROC_XI-1,0:NPROC_ETA-1) :: addressing

! 2-D addressing and buffers for summation between slices
  integer, dimension(NGLOB2DMAX_XMIN_XMAX_CM) :: iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_CM) :: iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle
  integer, dimension(NGLOB2DMAX_XMIN_XMAX_IC) :: iboolleft_xi_inner_core,iboolright_xi_inner_core
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_IC) :: iboolleft_eta_inner_core,iboolright_eta_inner_core

! indirect addressing for each corner of the chunks
  integer, dimension(NGLOB1D_RADIAL_crust_mantle,NUMCORNERS_SHARED) :: iboolcorner_crust_mantle
  integer, dimension(NGLOB1D_RADIAL_inner_core,NUMCORNERS_SHARED) :: iboolcorner_inner_core
  integer icount_corners

  integer :: npoin2D_max_all_CM_IC
  integer, dimension(NGLOB2DMAX_XY_VAL,NUMFACES_SHARED) :: iboolfaces_crust_mantle
  integer, dimension(NGLOB2DMAX_XY_VAL,NUMFACES_SHARED) :: iboolfaces_inner_core
  real(kind=CUSTOM_REAL), dimension(NDIM,npoin2D_max_all_CM_IC,NUMFACES_SHARED) :: &
      buffer_send_faces_vector,buffer_received_faces_vector

! buffers for send and receive between corners of the chunks
! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB1D_RADIAL_crust_mantle + NGLOB1D_RADIAL_inner_core) :: &
    buffer_send_chunkcorners_vector,buffer_recv_chunkcorners_vector

! ---- arrays to assemble between chunks

! communication pattern for faces between chunks
  integer, dimension(NUMMSGS_FACES) :: iprocfrom_faces,iprocto_faces

! communication pattern for corners between chunks
  integer, dimension(NCORNERSCHUNKS) :: iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners

! MPI status of messages to be received
  integer, dimension(MPI_STATUS_SIZE) :: msg_status

  integer :: ipoin,ipoin2D,ipoin1D
  integer :: sender,receiver
  integer :: imsg
  integer :: icount_faces,npoin2D_chunks_all

  integer :: NGLOB1D_RADIAL_all
  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_all,npoin2D_eta_all

! do not remove the "save" statement because this routine is non blocking
! therefore it needs to find the right value of ioffset when it re-enters
! the routine later to perform the next communication step
  integer, save :: ioffset

  integer :: ier
! do not remove the "save" statement because this routine is non blocking
  integer, save :: request_send,request_receive
  integer, dimension(NUMFACES_SHARED), save :: request_send_array,request_receive_array
  logical :: flag_result_test

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! check flag to see if we need to assemble (might be turned off when debugging)
  if (.not. ACTUALLY_ASSEMBLE_MPI_SLICES) then
    iphase = 9999 ! this means everything is finished
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

  if(iphase == 1) then

! slices copy the right face into the buffer
  do ipoin = 1,npoin2D_xi_crust_mantle(2)
    buffer_send_faces_vector(1,ipoin,1) = accel_crust_mantle(1,iboolright_xi_crust_mantle(ipoin))
    buffer_send_faces_vector(2,ipoin,1) = accel_crust_mantle(2,iboolright_xi_crust_mantle(ipoin))
    buffer_send_faces_vector(3,ipoin,1) = accel_crust_mantle(3,iboolright_xi_crust_mantle(ipoin))
  enddo

! the buffer for the inner core starts right after the buffer for the crust and mantle
  ioffset = npoin2D_xi_crust_mantle(2)

  do ipoin = 1,npoin2D_xi_inner_core(2)
    buffer_send_faces_vector(1,ioffset + ipoin,1) = accel_inner_core(1,iboolright_xi_inner_core(ipoin))
    buffer_send_faces_vector(2,ioffset + ipoin,1) = accel_inner_core(2,iboolright_xi_inner_core(ipoin))
    buffer_send_faces_vector(3,ioffset + ipoin,1) = accel_inner_core(3,iboolright_xi_inner_core(ipoin))
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
  call MPI_IRECV(buffer_received_faces_vector,NDIM*npoin2D_xi_all(1),CUSTOM_MPI_TYPE,sender, &
        itag,MPI_COMM_WORLD,request_receive,ier)

  call MPI_ISSEND(buffer_send_faces_vector,NDIM*npoin2D_xi_all(2),CUSTOM_MPI_TYPE,receiver, &
        itag2,MPI_COMM_WORLD,request_send,ier)

  iphase = iphase + 1
  return ! exit because we have started some communications therefore we need some time

  endif !!!!!!!!! end of iphase 1

  if(iphase == 2) then

! call MPI_WAIT(request_send,msg_status,ier)
! call MPI_WAIT(request_receive,msg_status,ier)
  call MPI_TEST(request_send,flag_result_test,msg_status,ier)
  if(.not. flag_result_test) return ! exit if message not sent yet
  call MPI_TEST(request_receive,flag_result_test,msg_status,ier)
  if(.not. flag_result_test) return ! exit if message not received yet

! all slices add the buffer received to the contributions on the left face
  if(iproc_xi > 0) then

  do ipoin = 1,npoin2D_xi_crust_mantle(1)
    accel_crust_mantle(1,iboolleft_xi_crust_mantle(ipoin)) = accel_crust_mantle(1,iboolleft_xi_crust_mantle(ipoin)) + &
                              buffer_received_faces_vector(1,ipoin,1)
    accel_crust_mantle(2,iboolleft_xi_crust_mantle(ipoin)) = accel_crust_mantle(2,iboolleft_xi_crust_mantle(ipoin)) + &
                              buffer_received_faces_vector(2,ipoin,1)
    accel_crust_mantle(3,iboolleft_xi_crust_mantle(ipoin)) = accel_crust_mantle(3,iboolleft_xi_crust_mantle(ipoin)) + &
                              buffer_received_faces_vector(3,ipoin,1)
  enddo

! the buffer for the inner core starts right after the buffer for the crust and mantle
  ioffset = npoin2D_xi_crust_mantle(1)

  do ipoin = 1,npoin2D_xi_inner_core(1)
    accel_inner_core(1,iboolleft_xi_inner_core(ipoin)) = accel_inner_core(1,iboolleft_xi_inner_core(ipoin)) + &
                              buffer_received_faces_vector(1,ioffset + ipoin,1)
    accel_inner_core(2,iboolleft_xi_inner_core(ipoin)) = accel_inner_core(2,iboolleft_xi_inner_core(ipoin)) + &
                              buffer_received_faces_vector(2,ioffset + ipoin,1)
    accel_inner_core(3,iboolleft_xi_inner_core(ipoin)) = accel_inner_core(3,iboolleft_xi_inner_core(ipoin)) + &
                              buffer_received_faces_vector(3,ioffset + ipoin,1)
  enddo

  endif

! the contributions are correctly assembled on the left side of each slice
! now we have to send the result back to the sender
! all slices copy the left face into the buffer
  do ipoin = 1,npoin2D_xi_crust_mantle(1)
    buffer_send_faces_vector(1,ipoin,1) = accel_crust_mantle(1,iboolleft_xi_crust_mantle(ipoin))
    buffer_send_faces_vector(2,ipoin,1) = accel_crust_mantle(2,iboolleft_xi_crust_mantle(ipoin))
    buffer_send_faces_vector(3,ipoin,1) = accel_crust_mantle(3,iboolleft_xi_crust_mantle(ipoin))
  enddo

! the buffer for the inner core starts right after the buffer for the crust and mantle
  ioffset = npoin2D_xi_crust_mantle(1)

  do ipoin = 1,npoin2D_xi_inner_core(1)
    buffer_send_faces_vector(1,ioffset + ipoin,1) = accel_inner_core(1,iboolleft_xi_inner_core(ipoin))
    buffer_send_faces_vector(2,ioffset + ipoin,1) = accel_inner_core(2,iboolleft_xi_inner_core(ipoin))
    buffer_send_faces_vector(3,ioffset + ipoin,1) = accel_inner_core(3,iboolleft_xi_inner_core(ipoin))
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

  call MPI_ISSEND(buffer_send_faces_vector,NDIM*npoin2D_xi_all(1),CUSTOM_MPI_TYPE,receiver, &
        itag2,MPI_COMM_WORLD,request_send,ier)

  iphase = iphase + 1
  return ! exit because we have started some communications therefore we need some time

  endif !!!!!!!!! end of iphase 2

  if(iphase == 3) then

! call MPI_WAIT(request_send,msg_status,ier)
! call MPI_WAIT(request_receive,msg_status,ier)
  call MPI_TEST(request_send,flag_result_test,msg_status,ier)
  if(.not. flag_result_test) return ! exit if message not sent yet
  call MPI_TEST(request_receive,flag_result_test,msg_status,ier)
  if(.not. flag_result_test) return ! exit if message not received yet

! all slices copy the buffer received to the contributions on the right face
  if(iproc_xi < NPROC_XI-1) then

  do ipoin = 1,npoin2D_xi_crust_mantle(2)
    accel_crust_mantle(1,iboolright_xi_crust_mantle(ipoin)) = buffer_received_faces_vector(1,ipoin,1)
    accel_crust_mantle(2,iboolright_xi_crust_mantle(ipoin)) = buffer_received_faces_vector(2,ipoin,1)
    accel_crust_mantle(3,iboolright_xi_crust_mantle(ipoin)) = buffer_received_faces_vector(3,ipoin,1)
  enddo

! the buffer for the inner core starts right after the buffer for the crust and mantle
  ioffset = npoin2D_xi_crust_mantle(2)

  do ipoin = 1,npoin2D_xi_inner_core(2)
    accel_inner_core(1,iboolright_xi_inner_core(ipoin)) = buffer_received_faces_vector(1,ioffset + ipoin,1)
    accel_inner_core(2,iboolright_xi_inner_core(ipoin)) = buffer_received_faces_vector(2,ioffset + ipoin,1)
    accel_inner_core(3,iboolright_xi_inner_core(ipoin)) = buffer_received_faces_vector(3,ioffset + ipoin,1)
  enddo

  endif

!----
!---- then assemble along eta using the 2-D topology
!----

! slices copy the right face into the buffer
  do ipoin = 1,npoin2D_eta_crust_mantle(2)
    buffer_send_faces_vector(1,ipoin,1) = accel_crust_mantle(1,iboolright_eta_crust_mantle(ipoin))
    buffer_send_faces_vector(2,ipoin,1) = accel_crust_mantle(2,iboolright_eta_crust_mantle(ipoin))
    buffer_send_faces_vector(3,ipoin,1) = accel_crust_mantle(3,iboolright_eta_crust_mantle(ipoin))
  enddo

! the buffer for the inner core starts right after the buffer for the crust and mantle
  ioffset = npoin2D_eta_crust_mantle(2)

  do ipoin = 1,npoin2D_eta_inner_core(2)
    buffer_send_faces_vector(1,ioffset + ipoin,1) = accel_inner_core(1,iboolright_eta_inner_core(ipoin))
    buffer_send_faces_vector(2,ioffset + ipoin,1) = accel_inner_core(2,iboolright_eta_inner_core(ipoin))
    buffer_send_faces_vector(3,ioffset + ipoin,1) = accel_inner_core(3,iboolright_eta_inner_core(ipoin))
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

  call MPI_ISSEND(buffer_send_faces_vector,NDIM*npoin2D_eta_all(2),CUSTOM_MPI_TYPE,receiver, &
    itag2,MPI_COMM_WORLD,request_send,ier)

  iphase = iphase + 1
  return ! exit because we have started some communications therefore we need some time

  endif !!!!!!!!! end of iphase 3

  if(iphase == 4) then

! call MPI_WAIT(request_send,msg_status,ier)
! call MPI_WAIT(request_receive,msg_status,ier)
  call MPI_TEST(request_send,flag_result_test,msg_status,ier)
  if(.not. flag_result_test) return ! exit if message not sent yet
  call MPI_TEST(request_receive,flag_result_test,msg_status,ier)
  if(.not. flag_result_test) return ! exit if message not received yet

! all slices add the buffer received to the contributions on the left face
  if(iproc_eta > 0) then

  do ipoin = 1,npoin2D_eta_crust_mantle(1)
    accel_crust_mantle(1,iboolleft_eta_crust_mantle(ipoin)) = accel_crust_mantle(1,iboolleft_eta_crust_mantle(ipoin)) + &
                              buffer_received_faces_vector(1,ipoin,1)
    accel_crust_mantle(2,iboolleft_eta_crust_mantle(ipoin)) = accel_crust_mantle(2,iboolleft_eta_crust_mantle(ipoin)) + &
                              buffer_received_faces_vector(2,ipoin,1)
    accel_crust_mantle(3,iboolleft_eta_crust_mantle(ipoin)) = accel_crust_mantle(3,iboolleft_eta_crust_mantle(ipoin)) + &
                              buffer_received_faces_vector(3,ipoin,1)
  enddo

! the buffer for the inner core starts right after the buffer for the crust and mantle
  ioffset = npoin2D_eta_crust_mantle(1)

  do ipoin = 1,npoin2D_eta_inner_core(1)
    accel_inner_core(1,iboolleft_eta_inner_core(ipoin)) = accel_inner_core(1,iboolleft_eta_inner_core(ipoin)) + &
                              buffer_received_faces_vector(1,ioffset + ipoin,1)
    accel_inner_core(2,iboolleft_eta_inner_core(ipoin)) = accel_inner_core(2,iboolleft_eta_inner_core(ipoin)) + &
                              buffer_received_faces_vector(2,ioffset + ipoin,1)
    accel_inner_core(3,iboolleft_eta_inner_core(ipoin)) = accel_inner_core(3,iboolleft_eta_inner_core(ipoin)) + &
                              buffer_received_faces_vector(3,ioffset + ipoin,1)
  enddo

  endif

! the contributions are correctly assembled on the left side of each slice
! now we have to send the result back to the sender
! all slices copy the left face into the buffer
  do ipoin = 1,npoin2D_eta_crust_mantle(1)
    buffer_send_faces_vector(1,ipoin,1) = accel_crust_mantle(1,iboolleft_eta_crust_mantle(ipoin))
    buffer_send_faces_vector(2,ipoin,1) = accel_crust_mantle(2,iboolleft_eta_crust_mantle(ipoin))
    buffer_send_faces_vector(3,ipoin,1) = accel_crust_mantle(3,iboolleft_eta_crust_mantle(ipoin))
  enddo

! the buffer for the inner core starts right after the buffer for the crust and mantle
  ioffset = npoin2D_eta_crust_mantle(1)

  do ipoin = 1,npoin2D_eta_inner_core(1)
    buffer_send_faces_vector(1,ioffset + ipoin,1) = accel_inner_core(1,iboolleft_eta_inner_core(ipoin))
    buffer_send_faces_vector(2,ioffset + ipoin,1) = accel_inner_core(2,iboolleft_eta_inner_core(ipoin))
    buffer_send_faces_vector(3,ioffset + ipoin,1) = accel_inner_core(3,iboolleft_eta_inner_core(ipoin))
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

  call MPI_ISSEND(buffer_send_faces_vector,NDIM*npoin2D_eta_all(1),CUSTOM_MPI_TYPE,receiver, &
    itag2,MPI_COMM_WORLD,request_send,ier)

  iphase = iphase + 1
  return ! exit because we have started some communications therefore we need some time

  endif !!!!!!!!! end of iphase 4

  if(iphase == 5) then

! call MPI_WAIT(request_send,msg_status,ier)
! call MPI_WAIT(request_receive,msg_status,ier)
  call MPI_TEST(request_send,flag_result_test,msg_status,ier)
  if(.not. flag_result_test) return ! exit if message not sent yet
  call MPI_TEST(request_receive,flag_result_test,msg_status,ier)
  if(.not. flag_result_test) return ! exit if message not received yet

! all slices copy the buffer received to the contributions on the right face
  if(iproc_eta < NPROC_ETA-1) then

  do ipoin = 1,npoin2D_eta_crust_mantle(2)
    accel_crust_mantle(1,iboolright_eta_crust_mantle(ipoin)) = buffer_received_faces_vector(1,ipoin,1)
    accel_crust_mantle(2,iboolright_eta_crust_mantle(ipoin)) = buffer_received_faces_vector(2,ipoin,1)
    accel_crust_mantle(3,iboolright_eta_crust_mantle(ipoin)) = buffer_received_faces_vector(3,ipoin,1)
  enddo

! the buffer for the inner core starts right after the buffer for the crust and mantle
  ioffset = npoin2D_eta_crust_mantle(2)

  do ipoin = 1,npoin2D_eta_inner_core(2)
    accel_inner_core(1,iboolright_eta_inner_core(ipoin)) = buffer_received_faces_vector(1,ioffset + ipoin,1)
    accel_inner_core(2,iboolright_eta_inner_core(ipoin)) = buffer_received_faces_vector(2,ioffset + ipoin,1)
    accel_inner_core(3,iboolright_eta_inner_core(ipoin)) = buffer_received_faces_vector(3,ioffset + ipoin,1)
  enddo

  endif

!----
!---- start MPI assembling phase between chunks
!----

! check flag to see if we need to assemble (might be turned off when debugging)
! and do not assemble if only one chunk
  if (.not. ACTUALLY_ASSEMBLE_MPI_CHUNKS .or. NCHUNKS == 1) then
    iphase = 9999 ! this means everything is finished
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

    call MPI_IRECV(buffer_received_faces_vector(:,:,icount_faces),NDIM*npoin2D_chunks_all,CUSTOM_MPI_TYPE,sender, &
              itag,MPI_COMM_WORLD,request_receive_array(icount_faces),ier)

!   do ipoin2D = 1,npoin2D_faces_crust_mantle(icount_faces)
!     accel_crust_mantle(1,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = &
!        accel_crust_mantle(1,iboolfaces_crust_mantle(ipoin2D,icount_faces)) + buffer_received_faces_vector(1,ipoin2D,icount_faces)
!     accel_crust_mantle(2,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = &
!        accel_crust_mantle(2,iboolfaces_crust_mantle(ipoin2D,icount_faces)) + buffer_received_faces_vector(2,ipoin2D,icount_faces)
!     accel_crust_mantle(3,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = &
!        accel_crust_mantle(3,iboolfaces_crust_mantle(ipoin2D,icount_faces)) + buffer_received_faces_vector(3,ipoin2D,icount_faces)
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
      buffer_send_faces_vector(1,ipoin2D,icount_faces) = accel_crust_mantle(1,iboolfaces_crust_mantle(ipoin2D,icount_faces))
      buffer_send_faces_vector(2,ipoin2D,icount_faces) = accel_crust_mantle(2,iboolfaces_crust_mantle(ipoin2D,icount_faces))
      buffer_send_faces_vector(3,ipoin2D,icount_faces) = accel_crust_mantle(3,iboolfaces_crust_mantle(ipoin2D,icount_faces))
    enddo

! the buffer for the inner core starts right after the buffer for the crust and mantle
    ioffset = npoin2D_faces_crust_mantle(icount_faces)

    do ipoin2D = 1,npoin2D_faces_inner_core(icount_faces)
      buffer_send_faces_vector(1,ioffset + ipoin2D,icount_faces) = accel_inner_core(1,iboolfaces_inner_core(ipoin2D,icount_faces))
      buffer_send_faces_vector(2,ioffset + ipoin2D,icount_faces) = accel_inner_core(2,iboolfaces_inner_core(ipoin2D,icount_faces))
      buffer_send_faces_vector(3,ioffset + ipoin2D,icount_faces) = accel_inner_core(3,iboolfaces_inner_core(ipoin2D,icount_faces))
    enddo

    call MPI_ISSEND(buffer_send_faces_vector(:,:,icount_faces),NDIM*npoin2D_chunks_all,CUSTOM_MPI_TYPE,receiver,itag, &
                     MPI_COMM_WORLD,request_send_array(icount_faces),ier)
  endif
  enddo

  iphase = iphase + 1
  return ! exit because we have started some communications therefore we need some time

  endif !!!!!!!!! end of iphase 5

  if(iphase == 6) then

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

    do ipoin2D = 1,npoin2D_faces_crust_mantle(icount_faces)
      accel_crust_mantle(1,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = &
         accel_crust_mantle(1,iboolfaces_crust_mantle(ipoin2D,icount_faces)) + buffer_received_faces_vector(1,ipoin2D,icount_faces)
      accel_crust_mantle(2,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = &
         accel_crust_mantle(2,iboolfaces_crust_mantle(ipoin2D,icount_faces)) + buffer_received_faces_vector(2,ipoin2D,icount_faces)
      accel_crust_mantle(3,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = &
         accel_crust_mantle(3,iboolfaces_crust_mantle(ipoin2D,icount_faces)) + buffer_received_faces_vector(3,ipoin2D,icount_faces)
    enddo

! the buffer for the inner core starts right after the buffer for the crust and mantle
    ioffset = npoin2D_faces_crust_mantle(icount_faces)

    do ipoin2D = 1,npoin2D_faces_inner_core(icount_faces)
      accel_inner_core(1,iboolfaces_inner_core(ipoin2D,icount_faces)) = &
         accel_inner_core(1,iboolfaces_inner_core(ipoin2D,icount_faces)) + &
           buffer_received_faces_vector(1,ioffset + ipoin2D,icount_faces)
      accel_inner_core(2,iboolfaces_inner_core(ipoin2D,icount_faces)) = &
         accel_inner_core(2,iboolfaces_inner_core(ipoin2D,icount_faces)) + &
           buffer_received_faces_vector(2,ioffset + ipoin2D,icount_faces)
      accel_inner_core(3,iboolfaces_inner_core(ipoin2D,icount_faces)) = &
         accel_inner_core(3,iboolfaces_inner_core(ipoin2D,icount_faces)) + &
           buffer_received_faces_vector(3,ioffset + ipoin2D,icount_faces)
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

    call MPI_IRECV(buffer_received_faces_vector(:,:,icount_faces),NDIM*npoin2D_chunks_all,CUSTOM_MPI_TYPE,sender, &
              itag,MPI_COMM_WORLD,request_receive_array(icount_faces),ier)

!   do ipoin2D = 1,npoin2D_faces_crust_mantle(icount_faces)
!     accel_crust_mantle(1,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = buffer_received_faces_vector(1,ipoin2D,icount_faces)
!     accel_crust_mantle(2,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = buffer_received_faces_vector(2,ipoin2D,icount_faces)
!     accel_crust_mantle(3,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = buffer_received_faces_vector(3,ipoin2D,icount_faces)
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
      buffer_send_faces_vector(1,ipoin2D,icount_faces) = accel_crust_mantle(1,iboolfaces_crust_mantle(ipoin2D,icount_faces))
      buffer_send_faces_vector(2,ipoin2D,icount_faces) = accel_crust_mantle(2,iboolfaces_crust_mantle(ipoin2D,icount_faces))
      buffer_send_faces_vector(3,ipoin2D,icount_faces) = accel_crust_mantle(3,iboolfaces_crust_mantle(ipoin2D,icount_faces))
    enddo

! the buffer for the inner core starts right after the buffer for the crust and mantle
    ioffset = npoin2D_faces_crust_mantle(icount_faces)

    do ipoin2D = 1,npoin2D_faces_inner_core(icount_faces)
      buffer_send_faces_vector(1,ioffset + ipoin2D,icount_faces) = accel_inner_core(1,iboolfaces_inner_core(ipoin2D,icount_faces))
      buffer_send_faces_vector(2,ioffset + ipoin2D,icount_faces) = accel_inner_core(2,iboolfaces_inner_core(ipoin2D,icount_faces))
      buffer_send_faces_vector(3,ioffset + ipoin2D,icount_faces) = accel_inner_core(3,iboolfaces_inner_core(ipoin2D,icount_faces))
    enddo

    call MPI_ISSEND(buffer_send_faces_vector(:,:,icount_faces),NDIM*npoin2D_chunks_all,CUSTOM_MPI_TYPE,receiver,itag, &
                     MPI_COMM_WORLD,request_send_array(icount_faces),ier)
  endif
  enddo

  iphase = iphase + 1
  return ! exit because we have started some communications therefore we need some time

  endif !!!!!!!!! end of iphase 6

  if(iphase == 7) then

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
    do ipoin2D = 1,npoin2D_faces_crust_mantle(icount_faces)
      accel_crust_mantle(1,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = buffer_received_faces_vector(1,ipoin2D,icount_faces)
      accel_crust_mantle(2,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = buffer_received_faces_vector(2,ipoin2D,icount_faces)
      accel_crust_mantle(3,iboolfaces_crust_mantle(ipoin2D,icount_faces)) = buffer_received_faces_vector(3,ipoin2D,icount_faces)
    enddo

! the buffer for the inner core starts right after the buffer for the crust and mantle
    ioffset = npoin2D_faces_crust_mantle(icount_faces)

    do ipoin2D = 1,npoin2D_faces_inner_core(icount_faces)
      accel_inner_core(1,iboolfaces_inner_core(ipoin2D,icount_faces)) = &
        buffer_received_faces_vector(1,ioffset + ipoin2D,icount_faces)
      accel_inner_core(2,iboolfaces_inner_core(ipoin2D,icount_faces)) = &
        buffer_received_faces_vector(2,ioffset + ipoin2D,icount_faces)
      accel_inner_core(3,iboolfaces_inner_core(ipoin2D,icount_faces)) = &
        buffer_received_faces_vector(3,ioffset + ipoin2D,icount_faces)
    enddo
  endif
  enddo

! this is the exit condition, to go beyond the last phase number
  iphase = iphase + 1

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

    call MPI_RECV(buffer_recv_chunkcorners_vector,NDIM*NGLOB1D_RADIAL_all, &
          CUSTOM_MPI_TYPE,sender,itag,MPI_COMM_WORLD,msg_status,ier)

    do ipoin1D = 1,NGLOB1D_RADIAL_crust_mantle
      accel_crust_mantle(1,iboolcorner_crust_mantle(ipoin1D,icount_corners)) = &
               accel_crust_mantle(1,iboolcorner_crust_mantle(ipoin1D,icount_corners)) + &
               buffer_recv_chunkcorners_vector(1,ipoin1D)
      accel_crust_mantle(2,iboolcorner_crust_mantle(ipoin1D,icount_corners)) = &
               accel_crust_mantle(2,iboolcorner_crust_mantle(ipoin1D,icount_corners)) + &
               buffer_recv_chunkcorners_vector(2,ipoin1D)
      accel_crust_mantle(3,iboolcorner_crust_mantle(ipoin1D,icount_corners)) = &
               accel_crust_mantle(3,iboolcorner_crust_mantle(ipoin1D,icount_corners)) + &
               buffer_recv_chunkcorners_vector(3,ipoin1D)
    enddo

    do ipoin1D = 1,NGLOB1D_RADIAL_inner_core
      accel_inner_core(1,iboolcorner_inner_core(ipoin1D,icount_corners)) = &
               accel_inner_core(1,iboolcorner_inner_core(ipoin1D,icount_corners)) + &
               buffer_recv_chunkcorners_vector(1,ioffset + ipoin1D)
      accel_inner_core(2,iboolcorner_inner_core(ipoin1D,icount_corners)) = &
               accel_inner_core(2,iboolcorner_inner_core(ipoin1D,icount_corners)) + &
               buffer_recv_chunkcorners_vector(2,ioffset + ipoin1D)
      accel_inner_core(3,iboolcorner_inner_core(ipoin1D,icount_corners)) = &
               accel_inner_core(3,iboolcorner_inner_core(ipoin1D,icount_corners)) + &
               buffer_recv_chunkcorners_vector(3,ioffset + ipoin1D)
    enddo

! receive from worker #2 and add to local array
  if(NCHUNKS /= 2) then

    sender = iproc_worker2_corners(imsg)

    call MPI_RECV(buffer_recv_chunkcorners_vector,NDIM*NGLOB1D_RADIAL_all, &
          CUSTOM_MPI_TYPE,sender,itag,MPI_COMM_WORLD,msg_status,ier)

    do ipoin1D = 1,NGLOB1D_RADIAL_crust_mantle
      accel_crust_mantle(1,iboolcorner_crust_mantle(ipoin1D,icount_corners)) = &
               accel_crust_mantle(1,iboolcorner_crust_mantle(ipoin1D,icount_corners)) + &
               buffer_recv_chunkcorners_vector(1,ipoin1D)
      accel_crust_mantle(2,iboolcorner_crust_mantle(ipoin1D,icount_corners)) = &
               accel_crust_mantle(2,iboolcorner_crust_mantle(ipoin1D,icount_corners)) + &
               buffer_recv_chunkcorners_vector(2,ipoin1D)
      accel_crust_mantle(3,iboolcorner_crust_mantle(ipoin1D,icount_corners)) = &
               accel_crust_mantle(3,iboolcorner_crust_mantle(ipoin1D,icount_corners)) + &
               buffer_recv_chunkcorners_vector(3,ipoin1D)
    enddo

    do ipoin1D = 1,NGLOB1D_RADIAL_inner_core
      accel_inner_core(1,iboolcorner_inner_core(ipoin1D,icount_corners)) = &
               accel_inner_core(1,iboolcorner_inner_core(ipoin1D,icount_corners)) + &
               buffer_recv_chunkcorners_vector(1,ioffset + ipoin1D)
      accel_inner_core(2,iboolcorner_inner_core(ipoin1D,icount_corners)) = &
               accel_inner_core(2,iboolcorner_inner_core(ipoin1D,icount_corners)) + &
               buffer_recv_chunkcorners_vector(2,ioffset + ipoin1D)
      accel_inner_core(3,iboolcorner_inner_core(ipoin1D,icount_corners)) = &
               accel_inner_core(3,iboolcorner_inner_core(ipoin1D,icount_corners)) + &
               buffer_recv_chunkcorners_vector(3,ioffset + ipoin1D)
    enddo

  endif

  endif

!---- send messages from the two workers to the master
  if(myrank==iproc_worker1_corners(imsg) .or. &
              (NCHUNKS /= 2 .and. myrank==iproc_worker2_corners(imsg))) then

    receiver = iproc_master_corners(imsg)

    do ipoin1D = 1,NGLOB1D_RADIAL_crust_mantle
      buffer_send_chunkcorners_vector(1,ipoin1D) = accel_crust_mantle(1,iboolcorner_crust_mantle(ipoin1D,icount_corners))
      buffer_send_chunkcorners_vector(2,ipoin1D) = accel_crust_mantle(2,iboolcorner_crust_mantle(ipoin1D,icount_corners))
      buffer_send_chunkcorners_vector(3,ipoin1D) = accel_crust_mantle(3,iboolcorner_crust_mantle(ipoin1D,icount_corners))
    enddo

    do ipoin1D = 1,NGLOB1D_RADIAL_inner_core
      buffer_send_chunkcorners_vector(1,ioffset + ipoin1D) = accel_inner_core(1,iboolcorner_inner_core(ipoin1D,icount_corners))
      buffer_send_chunkcorners_vector(2,ioffset + ipoin1D) = accel_inner_core(2,iboolcorner_inner_core(ipoin1D,icount_corners))
      buffer_send_chunkcorners_vector(3,ioffset + ipoin1D) = accel_inner_core(3,iboolcorner_inner_core(ipoin1D,icount_corners))
    enddo

    call MPI_SEND(buffer_send_chunkcorners_vector,NDIM*NGLOB1D_RADIAL_all,CUSTOM_MPI_TYPE,receiver,itag,MPI_COMM_WORLD,ier)

  endif

! *********************************************************************
!  transmit messages back in opposite direction (master -> two workers)
! *********************************************************************

!---- receive messages from the master on the two workers
  if(myrank==iproc_worker1_corners(imsg) .or. &
              (NCHUNKS /= 2 .and. myrank==iproc_worker2_corners(imsg))) then

! receive from master and copy to local array
    sender = iproc_master_corners(imsg)

    call MPI_RECV(buffer_recv_chunkcorners_vector,NDIM*NGLOB1D_RADIAL_all, &
          CUSTOM_MPI_TYPE,sender,itag,MPI_COMM_WORLD,msg_status,ier)

    do ipoin1D = 1,NGLOB1D_RADIAL_crust_mantle
      accel_crust_mantle(1,iboolcorner_crust_mantle(ipoin1D,icount_corners)) = buffer_recv_chunkcorners_vector(1,ipoin1D)
      accel_crust_mantle(2,iboolcorner_crust_mantle(ipoin1D,icount_corners)) = buffer_recv_chunkcorners_vector(2,ipoin1D)
      accel_crust_mantle(3,iboolcorner_crust_mantle(ipoin1D,icount_corners)) = buffer_recv_chunkcorners_vector(3,ipoin1D)
    enddo

    do ipoin1D = 1,NGLOB1D_RADIAL_inner_core
      accel_inner_core(1,iboolcorner_inner_core(ipoin1D,icount_corners)) = buffer_recv_chunkcorners_vector(1,ioffset + ipoin1D)
      accel_inner_core(2,iboolcorner_inner_core(ipoin1D,icount_corners)) = buffer_recv_chunkcorners_vector(2,ioffset + ipoin1D)
      accel_inner_core(3,iboolcorner_inner_core(ipoin1D,icount_corners)) = buffer_recv_chunkcorners_vector(3,ioffset + ipoin1D)
    enddo

  endif

!---- send messages from the master to the two workers
  if(myrank==iproc_master_corners(imsg)) then

    do ipoin1D = 1,NGLOB1D_RADIAL_crust_mantle
      buffer_send_chunkcorners_vector(1,ipoin1D) = accel_crust_mantle(1,iboolcorner_crust_mantle(ipoin1D,icount_corners))
      buffer_send_chunkcorners_vector(2,ipoin1D) = accel_crust_mantle(2,iboolcorner_crust_mantle(ipoin1D,icount_corners))
      buffer_send_chunkcorners_vector(3,ipoin1D) = accel_crust_mantle(3,iboolcorner_crust_mantle(ipoin1D,icount_corners))
    enddo

    do ipoin1D = 1,NGLOB1D_RADIAL_inner_core
      buffer_send_chunkcorners_vector(1,ioffset + ipoin1D) = accel_inner_core(1,iboolcorner_inner_core(ipoin1D,icount_corners))
      buffer_send_chunkcorners_vector(2,ioffset + ipoin1D) = accel_inner_core(2,iboolcorner_inner_core(ipoin1D,icount_corners))
      buffer_send_chunkcorners_vector(3,ioffset + ipoin1D) = accel_inner_core(3,iboolcorner_inner_core(ipoin1D,icount_corners))
    enddo

! send to worker #1
    receiver = iproc_worker1_corners(imsg)
    call MPI_SEND(buffer_send_chunkcorners_vector,NDIM*NGLOB1D_RADIAL_all,CUSTOM_MPI_TYPE,receiver,itag,MPI_COMM_WORLD,ier)

! send to worker #2
  if(NCHUNKS /= 2) then
    receiver = iproc_worker2_corners(imsg)
    call MPI_SEND(buffer_send_chunkcorners_vector,NDIM*NGLOB1D_RADIAL_all,CUSTOM_MPI_TYPE,receiver,itag,MPI_COMM_WORLD,ier)

  endif

  endif

  enddo

  endif !!!!!!!!! end of iphase 7

  end subroutine assemble_MPI_vector

