!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 5
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
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
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_slave1_corners,iproc_slave2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar, &
            buffer_send_chunkcorners_scalar,buffer_recv_chunkcorners_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI,NPROC_ETA,NPOIN1D_RADIAL, &
            NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NPOIN2DMAX_XY,NCHUNKS)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer myrank,nglob,NCHUNKS

! array to assemble
  real(kind=CUSTOM_REAL), dimension(nglob) :: array_val

  integer iproc_xi,iproc_eta,ichunk
  integer npoin2D_xi,npoin2D_eta
  integer npoin2D_faces(NUMFACES_SHARED)

  integer NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NPOIN2DMAX_XY
  integer NPROC_XI,NPROC_ETA,NPOIN1D_RADIAL
  integer NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS

! for addressing of the slices
  integer, dimension(NCHUNKS,0:NPROC_XI-1,0:NPROC_ETA-1) :: addressing

! 2-D addressing and buffers for summation between slices
  integer, dimension(NPOIN2DMAX_XMIN_XMAX) :: iboolleft_xi,iboolright_xi
  integer, dimension(NPOIN2DMAX_YMIN_YMAX) :: iboolleft_eta,iboolright_eta

! indirect addressing for each corner of the chunks
  integer, dimension(NPOIN1D_RADIAL,NUMCORNERS_SHARED) :: iboolcorner
  integer icount_corners

  integer, dimension(NPOIN2DMAX_XY,NUMFACES_SHARED) :: iboolfaces
  real(kind=CUSTOM_REAL), dimension(NPOIN2DMAX_XY) :: buffer_send_faces_scalar,buffer_received_faces_scalar

! buffers for send and receive between corners of the chunks
  real(kind=CUSTOM_REAL), dimension(NPOIN1D_RADIAL) :: buffer_send_chunkcorners_scalar,buffer_recv_chunkcorners_scalar

! ---- arrays to assemble between chunks

! communication pattern for faces between chunks
  integer, dimension(NUMMSGS_FACES) :: iprocfrom_faces,iprocto_faces,imsg_type

! communication pattern for corners between chunks
  integer, dimension(NCORNERSCHUNKS) :: iproc_master_corners,iproc_slave1_corners,iproc_slave2_corners

! MPI status of messages to be received
  integer msg_status(MPI_STATUS_SIZE)

  integer ipoin,ipoin2D,ipoin1D
  integer isender,ireceiver,ier
  integer imsg,imsg_loop
  integer icount_faces,npoin2D_chunks

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! check flag to see if we need to assemble (might be turned off when debugging)
  if (.not. ACTUALLY_ASSEMBLE_MPI_SLICES) return

! here we have to assemble all the contributions between slices using MPI

!----
!---- assemble the contributions between slices using MPI
!----

!----
!---- first assemble along xi using the 2-D topology
!----

! assemble along xi only if more than one slice
  if(NPROC_XI > 1) then

! slices copy the right face into the buffer
  do ipoin=1,npoin2D_xi
    buffer_send_faces_scalar(ipoin) = array_val(iboolright_xi(ipoin))
  enddo

! send messages forward along each row
  if(iproc_xi == 0) then
    isender = MPI_PROC_NULL
  else
    isender = addressing(ichunk,iproc_xi - 1,iproc_eta)
  endif
  if(iproc_xi == NPROC_XI-1) then
    ireceiver = MPI_PROC_NULL
  else
    ireceiver = addressing(ichunk,iproc_xi + 1,iproc_eta)
  endif
  call MPI_SENDRECV(buffer_send_faces_scalar,npoin2D_xi,CUSTOM_MPI_TYPE,ireceiver, &
        itag2,buffer_received_faces_scalar,npoin2D_xi,CUSTOM_MPI_TYPE,isender, &
        itag,MPI_COMM_WORLD,msg_status,ier)

! all slices add the buffer received to the contributions on the left face
  if(iproc_xi > 0) then
  do ipoin=1,npoin2D_xi
    array_val(iboolleft_xi(ipoin)) = array_val(iboolleft_xi(ipoin)) + &
                              buffer_received_faces_scalar(ipoin)
  enddo
  endif

! the contributions are correctly assembled on the left side of each slice
! now we have to send the result back to the sender
! all slices copy the left face into the buffer
  do ipoin=1,npoin2D_xi
    buffer_send_faces_scalar(ipoin) = array_val(iboolleft_xi(ipoin))
  enddo

! send messages backward along each row
  if(iproc_xi == NPROC_XI-1) then
    isender = MPI_PROC_NULL
  else
    isender = addressing(ichunk,iproc_xi + 1,iproc_eta)
  endif
  if(iproc_xi == 0) then
    ireceiver = MPI_PROC_NULL
  else
    ireceiver = addressing(ichunk,iproc_xi - 1,iproc_eta)
  endif
  call MPI_SENDRECV(buffer_send_faces_scalar,npoin2D_xi,CUSTOM_MPI_TYPE,ireceiver, &
        itag2,buffer_received_faces_scalar,npoin2D_xi,CUSTOM_MPI_TYPE,isender, &
        itag,MPI_COMM_WORLD,msg_status,ier)

! all slices copy the buffer received to the contributions on the right face
  if(iproc_xi < NPROC_XI-1) then
  do ipoin=1,npoin2D_xi
    array_val(iboolright_xi(ipoin)) = buffer_received_faces_scalar(ipoin)
  enddo
  endif

  endif

!----
!---- then assemble along eta using the 2-D topology
!----

! assemble along eta only if more than one slice
  if(NPROC_ETA > 1) then

! slices copy the right face into the buffer
  do ipoin=1,npoin2D_eta
    buffer_send_faces_scalar(ipoin) = array_val(iboolright_eta(ipoin))
  enddo

! send messages forward along each row
  if(iproc_eta == 0) then
    isender = MPI_PROC_NULL
  else
    isender = addressing(ichunk,iproc_xi,iproc_eta - 1)
  endif
  if(iproc_eta == NPROC_ETA-1) then
    ireceiver = MPI_PROC_NULL
  else
    ireceiver = addressing(ichunk,iproc_xi,iproc_eta + 1)
  endif
  call MPI_SENDRECV(buffer_send_faces_scalar,npoin2D_eta,CUSTOM_MPI_TYPE,ireceiver, &
    itag2,buffer_received_faces_scalar,npoin2D_eta,CUSTOM_MPI_TYPE,isender, &
    itag,MPI_COMM_WORLD,msg_status,ier)

! all slices add the buffer received to the contributions on the left face
  if(iproc_eta > 0) then
  do ipoin=1,npoin2D_eta
    array_val(iboolleft_eta(ipoin)) = array_val(iboolleft_eta(ipoin)) + &
                              buffer_received_faces_scalar(ipoin)
  enddo
  endif

! the contributions are correctly assembled on the left side of each slice
! now we have to send the result back to the sender
! all slices copy the left face into the buffer
  do ipoin=1,npoin2D_eta
    buffer_send_faces_scalar(ipoin) = array_val(iboolleft_eta(ipoin))
  enddo

! send messages backward along each row
  if(iproc_eta == NPROC_ETA-1) then
    isender = MPI_PROC_NULL
  else
    isender = addressing(ichunk,iproc_xi,iproc_eta + 1)
  endif
  if(iproc_eta == 0) then
    ireceiver = MPI_PROC_NULL
  else
    ireceiver = addressing(ichunk,iproc_xi,iproc_eta - 1)
  endif
  call MPI_SENDRECV(buffer_send_faces_scalar,npoin2D_eta,CUSTOM_MPI_TYPE,ireceiver, &
    itag2,buffer_received_faces_scalar,npoin2D_eta,CUSTOM_MPI_TYPE,isender, &
    itag,MPI_COMM_WORLD,msg_status,ier)

! all slices copy the buffer received to the contributions on the right face
  if(iproc_eta < NPROC_ETA-1) then
  do ipoin=1,npoin2D_eta
    array_val(iboolright_eta(ipoin)) = buffer_received_faces_scalar(ipoin)
  enddo
  endif

  endif

!----
!---- start MPI assembling phase between chunks
!----

! check flag to see if we need to assemble (might be turned off when debugging)
! and do not assemble if only one chunk
  if (.not. ACTUALLY_ASSEMBLE_MPI_CHUNKS .or. NCHUNKS == 1) return

! ***************************************************************
!  transmit messages in forward direction (iprocfrom -> iprocto)
! ***************************************************************

!---- put slices in receive mode
!---- a given slice can belong to at most two faces

! use three step scheme that can never deadlock
! scheme for faces cannot deadlock even if NPROC_XI = NPROC_ETA = 1
  do imsg_loop = 1,NUM_MSG_TYPES

  icount_faces = 0
  do imsg = 1,NUMMSGS_FACES
  if(myrank==iprocfrom_faces(imsg) .or. &
       myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
  if(myrank==iprocto_faces(imsg) .and. imsg_type(imsg) == imsg_loop) then
    isender = iprocfrom_faces(imsg)
    npoin2D_chunks = npoin2D_faces(icount_faces)
    call MPI_RECV(buffer_received_faces_scalar, &
              npoin2D_chunks,CUSTOM_MPI_TYPE,isender, &
              itag,MPI_COMM_WORLD,msg_status,ier)
    do ipoin2D=1,npoin2D_chunks
      array_val(iboolfaces(ipoin2D,icount_faces)) = &
         array_val(iboolfaces(ipoin2D,icount_faces)) + buffer_received_faces_scalar(ipoin2D)
    enddo
  endif
  enddo

!---- put slices in send mode
!---- a given slice can belong to at most two faces
  icount_faces = 0
  do imsg = 1,NUMMSGS_FACES
  if(myrank==iprocfrom_faces(imsg) .or. &
       myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
  if(myrank==iprocfrom_faces(imsg) .and. imsg_type(imsg) == imsg_loop) then
    ireceiver = iprocto_faces(imsg)
    npoin2D_chunks = npoin2D_faces(icount_faces)
    do ipoin2D=1,npoin2D_chunks
      buffer_send_faces_scalar(ipoin2D) = array_val(iboolfaces(ipoin2D,icount_faces))
    enddo
    call MPI_SEND(buffer_send_faces_scalar,npoin2D_chunks, &
              CUSTOM_MPI_TYPE,ireceiver,itag,MPI_COMM_WORLD,ier)
  endif
  enddo

! *********************************************************************
!  transmit messages back in opposite direction (iprocto -> iprocfrom)
! *********************************************************************

!---- put slices in receive mode
!---- a given slice can belong to at most two faces

  icount_faces = 0
  do imsg = 1,NUMMSGS_FACES
  if(myrank==iprocfrom_faces(imsg) .or. &
       myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
  if(myrank==iprocfrom_faces(imsg) .and. imsg_type(imsg) == imsg_loop) then
    isender = iprocto_faces(imsg)
    npoin2D_chunks = npoin2D_faces(icount_faces)
    call MPI_RECV(buffer_received_faces_scalar, &
              npoin2D_chunks,CUSTOM_MPI_TYPE,isender, &
              itag,MPI_COMM_WORLD,msg_status,ier)
    do ipoin2D=1,npoin2D_chunks
      array_val(iboolfaces(ipoin2D,icount_faces)) = buffer_received_faces_scalar(ipoin2D)
    enddo
  endif
  enddo

!---- put slices in send mode
!---- a given slice can belong to at most two faces
  icount_faces = 0
  do imsg = 1,NUMMSGS_FACES
  if(myrank==iprocfrom_faces(imsg) .or. &
       myrank==iprocto_faces(imsg)) icount_faces = icount_faces + 1
  if(myrank==iprocto_faces(imsg) .and. imsg_type(imsg) == imsg_loop) then
    ireceiver = iprocfrom_faces(imsg)
    npoin2D_chunks = npoin2D_faces(icount_faces)
    do ipoin2D=1,npoin2D_chunks
      buffer_send_faces_scalar(ipoin2D) = array_val(iboolfaces(ipoin2D,icount_faces))
    enddo
    call MPI_SEND(buffer_send_faces_scalar,npoin2D_chunks, &
              CUSTOM_MPI_TYPE,ireceiver,itag,MPI_COMM_WORLD,ier)
  endif
  enddo

! end of anti-deadlocking loop
  enddo

!----
!---- start MPI assembling corners
!----

! scheme for corners cannot deadlock even if NPROC_XI = NPROC_ETA = 1

! ***************************************************************
!  transmit messages in forward direction (two slaves -> master)
! ***************************************************************

  icount_corners = 0

  do imsg = 1,NCORNERSCHUNKS

  if(myrank == iproc_master_corners(imsg) .or. &
     myrank == iproc_slave1_corners(imsg) .or. &
     (NCHUNKS /= 2 .and. myrank == iproc_slave2_corners(imsg))) icount_corners = icount_corners + 1

!---- receive messages from the two slaves on the master
  if(myrank==iproc_master_corners(imsg)) then

! receive from slave #1 and add to local array
    isender = iproc_slave1_corners(imsg)
    call MPI_RECV(buffer_recv_chunkcorners_scalar,NPOIN1D_RADIAL, &
          CUSTOM_MPI_TYPE,isender,itag,MPI_COMM_WORLD,msg_status,ier)
    do ipoin1D=1,NPOIN1D_RADIAL
      array_val(iboolcorner(ipoin1D,icount_corners)) = array_val(iboolcorner(ipoin1D,icount_corners)) + &
               buffer_recv_chunkcorners_scalar(ipoin1D)
    enddo

! receive from slave #2 and add to local array
  if(NCHUNKS /= 2) then
    isender = iproc_slave2_corners(imsg)
    call MPI_RECV(buffer_recv_chunkcorners_scalar,NPOIN1D_RADIAL, &
          CUSTOM_MPI_TYPE,isender,itag,MPI_COMM_WORLD,msg_status,ier)
    do ipoin1D=1,NPOIN1D_RADIAL
      array_val(iboolcorner(ipoin1D,icount_corners)) = array_val(iboolcorner(ipoin1D,icount_corners)) + &
               buffer_recv_chunkcorners_scalar(ipoin1D)
    enddo
  endif

  endif

!---- send messages from the two slaves to the master
  if(myrank==iproc_slave1_corners(imsg) .or. &
              (NCHUNKS /= 2 .and. myrank==iproc_slave2_corners(imsg))) then

    ireceiver = iproc_master_corners(imsg)
    do ipoin1D=1,NPOIN1D_RADIAL
      buffer_send_chunkcorners_scalar(ipoin1D) = array_val(iboolcorner(ipoin1D,icount_corners))
    enddo
    call MPI_SEND(buffer_send_chunkcorners_scalar,NPOIN1D_RADIAL,CUSTOM_MPI_TYPE, &
              ireceiver,itag,MPI_COMM_WORLD,ier)

  endif

! *********************************************************************
!  transmit messages back in opposite direction (master -> two slaves)
! *********************************************************************

!---- receive messages from the master on the two slaves
  if(myrank==iproc_slave1_corners(imsg) .or. &
              (NCHUNKS /= 2 .and. myrank==iproc_slave2_corners(imsg))) then

! receive from master and copy to local array
    isender = iproc_master_corners(imsg)
    call MPI_RECV(buffer_recv_chunkcorners_scalar,NPOIN1D_RADIAL, &
          CUSTOM_MPI_TYPE,isender,itag,MPI_COMM_WORLD,msg_status,ier)
    do ipoin1D=1,NPOIN1D_RADIAL
      array_val(iboolcorner(ipoin1D,icount_corners)) = buffer_recv_chunkcorners_scalar(ipoin1D)
    enddo

  endif

!---- send messages from the master to the two slaves
  if(myrank==iproc_master_corners(imsg)) then

    do ipoin1D=1,NPOIN1D_RADIAL
      buffer_send_chunkcorners_scalar(ipoin1D) = array_val(iboolcorner(ipoin1D,icount_corners))
    enddo

! send to slave #1
    ireceiver = iproc_slave1_corners(imsg)
    call MPI_SEND(buffer_send_chunkcorners_scalar,NPOIN1D_RADIAL,CUSTOM_MPI_TYPE, &
              ireceiver,itag,MPI_COMM_WORLD,ier)

! send to slave #2
  if(NCHUNKS /= 2) then
    ireceiver = iproc_slave2_corners(imsg)
    call MPI_SEND(buffer_send_chunkcorners_scalar,NPOIN1D_RADIAL,CUSTOM_MPI_TYPE, &
              ireceiver,itag,MPI_COMM_WORLD,ier)
  endif

  endif

  enddo

  end subroutine assemble_MPI_scalar

