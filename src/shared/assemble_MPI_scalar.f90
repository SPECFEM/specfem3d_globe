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

!----
!---- assemble the contributions between slices and chunks using MPI
!----

  subroutine assemble_MPI_scalar(NPROC,nglob, &
                                 array_val, &
                                 num_interfaces,max_nibool_interfaces, &
                                 nibool_interfaces,ibool_interfaces, &
                                 my_neighbors)

! blocking send/receive

  use constants, only: CUSTOM_REAL,itag

  implicit none

  integer, intent(in) :: NPROC
  integer, intent(in) :: nglob

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(nglob), intent(inout) :: array_val

  integer, intent(in) :: num_interfaces,max_nibool_interfaces
  integer, dimension(num_interfaces), intent(in) :: nibool_interfaces,my_neighbors
  integer, dimension(max_nibool_interfaces,num_interfaces), intent(in) :: ibool_interfaces

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_send_scalar
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_recv_scalar
  integer, dimension(:), allocatable :: request_send_scalar
  integer, dimension(:), allocatable :: request_recv_scalar

  integer :: ipoin,iinterface,ier

! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC > 1) then

    allocate(buffer_send_scalar(max_nibool_interfaces,num_interfaces),stat=ier)
    if (ier /= 0 ) stop 'Error allocating array buffer_send_scalar'
    allocate(buffer_recv_scalar(max_nibool_interfaces,num_interfaces),stat=ier)
    if (ier /= 0 ) stop 'Error allocating array buffer_recv_scalar'
    allocate(request_send_scalar(num_interfaces),stat=ier)
    if (ier /= 0 ) stop 'Error allocating array request_send_scalar'
    allocate(request_recv_scalar(num_interfaces),stat=ier)
    if (ier /= 0 ) stop 'Error allocating array request_recv_scalar'

    ! partition border copy into the buffer
    do iinterface = 1, num_interfaces
      do ipoin = 1, nibool_interfaces(iinterface)
        buffer_send_scalar(ipoin,iinterface) = array_val(ibool_interfaces(ipoin,iinterface))
      enddo
    enddo

    ! send messages
    do iinterface = 1, num_interfaces
      ! non-blocking synchronous send request
      call isend_cr(buffer_send_scalar(1:nibool_interfaces(iinterface),iinterface), &
                    nibool_interfaces(iinterface), &
                    my_neighbors(iinterface), &
                    itag,request_send_scalar(iinterface) )

      ! receive request
      call irecv_cr(buffer_recv_scalar(1:nibool_interfaces(iinterface),iinterface), &
                    nibool_interfaces(iinterface), &
                    my_neighbors(iinterface), &
                    itag,request_recv_scalar(iinterface) )
    enddo

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces
      call wait_req(request_recv_scalar(iinterface))
    enddo

    ! adding contributions of neighbors
    do iinterface = 1, num_interfaces
      do ipoin = 1, nibool_interfaces(iinterface)
        array_val(ibool_interfaces(ipoin,iinterface)) = &
             array_val(ibool_interfaces(ipoin,iinterface)) + buffer_recv_scalar(ipoin,iinterface)
      enddo
    enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces
      call wait_req(request_send_scalar(iinterface))
    enddo

    deallocate(buffer_send_scalar)
    deallocate(buffer_recv_scalar)
    deallocate(request_send_scalar)
    deallocate(request_recv_scalar)

  endif

  end subroutine assemble_MPI_scalar


!-------------------------------------------------------------------------------------------------
!
! non-blocking routines
!
!-------------------------------------------------------------------------------------------------

  subroutine assemble_MPI_scalar_s(NPROC,nglob, &
                                   array_val, &
                                   buffer_send_scalar,buffer_recv_scalar, &
                                   num_interfaces,max_nibool_interfaces, &
                                   nibool_interfaces,ibool_interfaces, &
                                   my_neighbors, &
                                   request_send_scalar,request_recv_scalar)

! non-blocking MPI send

  use constants, only: CUSTOM_REAL,itag

  implicit none

  integer, intent(in) :: NPROC
  integer, intent(in) :: nglob
  integer, intent(in) :: num_interfaces,max_nibool_interfaces

  ! array to send
  real(kind=CUSTOM_REAL), dimension(nglob) :: array_val

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces,num_interfaces) :: buffer_send_scalar,buffer_recv_scalar

  integer, dimension(num_interfaces),intent(in) :: nibool_interfaces,my_neighbors
  integer, dimension(max_nibool_interfaces,num_interfaces),intent(in) :: ibool_interfaces
  integer, dimension(num_interfaces) :: request_send_scalar,request_recv_scalar

  ! local parameters
  integer :: ipoin,iinterface

! sends only if more than one partition
  if (NPROC > 1) then

    ! partition border copy into the buffer
    do iinterface = 1, num_interfaces
      do ipoin = 1, nibool_interfaces(iinterface)
        buffer_send_scalar(ipoin,iinterface) = array_val(ibool_interfaces(ipoin,iinterface))
      enddo
    enddo

    ! send messages
    do iinterface = 1, num_interfaces
      ! non-blocking synchronous send request
      call isend_cr(buffer_send_scalar(1:nibool_interfaces(iinterface),iinterface), &
                    nibool_interfaces(iinterface), &
                    my_neighbors(iinterface), &
                    itag,request_send_scalar(iinterface))

      ! receive request
      call irecv_cr(buffer_recv_scalar(1:nibool_interfaces(iinterface),iinterface), &
                    nibool_interfaces(iinterface), &
                    my_neighbors(iinterface), &
                    itag,request_recv_scalar(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_scalar_s

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_scalar_w(NPROC,nglob, &
                                   array_val, &
                                   buffer_recv_scalar,num_interfaces, &
                                   max_nibool_interfaces, &
                                   nibool_interfaces,ibool_interfaces, &
                                   my_neighbors, &
                                   request_send_scalar,request_recv_scalar)

! waits for send/receiver to be completed and assembles contributions

  use constants, only: myrank,CUSTOM_REAL,itag,DO_ORDERED_ASSEMBLY

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: nglob
  integer,intent(in) :: num_interfaces,max_nibool_interfaces
! array to assemble
  real(kind=CUSTOM_REAL), dimension(nglob) :: array_val

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces,num_interfaces) :: buffer_recv_scalar

  integer, dimension(num_interfaces),intent(in) :: nibool_interfaces
  integer, dimension(max_nibool_interfaces,num_interfaces),intent(in) :: ibool_interfaces
  integer, dimension(num_interfaces) :: request_send_scalar,request_recv_scalar

  integer, dimension(num_interfaces),intent(in) :: my_neighbors

  ! local parameters
  integer :: ipoin,iinterface
  integer :: iglob
  ! ordered assembly
  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces,num_interfaces) :: mybuffer
  logical :: need_add_my_contrib

! assemble only if more than one partition
  if (NPROC > 1) then

    ! ordered assembly
    if (DO_ORDERED_ASSEMBLY) then
      ! stores own MPI buffer values and set them to zero right away to avoid counting it more than once during assembly:
      ! buffers of higher rank get zeros on nodes shared with current buffer
      !
      ! move interface values of array_val to local buffers
      do iinterface = 1, num_interfaces
        do ipoin = 1, nibool_interfaces(iinterface)
          iglob = ibool_interfaces(ipoin,iinterface)
          mybuffer(ipoin,iinterface) = array_val(iglob)
          array_val(iglob) = 0._CUSTOM_REAL
        enddo
      enddo
    endif

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces
      call wait_req(request_recv_scalar(iinterface))
    enddo

    ! adding contributions of neighbors
    ! ordered assembly
    if (DO_ORDERED_ASSEMBLY) then
      ! The goal of this version is to avoid different round-off errors in different processors.
      ! The contribution of each processor is added following the order of its rank.
      ! This guarantees that the sums are done in the same order on all processors.
      !
      ! adding all contributions in order of processor rank
      need_add_my_contrib = .true.
      do iinterface = 1, num_interfaces
        if (need_add_my_contrib .and. myrank < my_neighbors(iinterface)) call add_my_contrib_scalar()
        do ipoin = 1, nibool_interfaces(iinterface)
          iglob = ibool_interfaces(ipoin,iinterface)
          array_val(iglob) = array_val(iglob) + buffer_recv_scalar(ipoin,iinterface)
        enddo
      enddo
      if (need_add_my_contrib) call add_my_contrib_scalar()

    else
      ! default assembly
      do iinterface = 1, num_interfaces
        do ipoin = 1, nibool_interfaces(iinterface)
          iglob = ibool_interfaces(ipoin,iinterface)
          array_val(iglob) = array_val(iglob) + buffer_recv_scalar(ipoin,iinterface)
        enddo
      enddo
    endif

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces
      call wait_req(request_send_scalar(iinterface))
    enddo

  endif

  contains

    subroutine add_my_contrib_scalar()

    integer :: my_iinterface,my_ipoin

    do my_iinterface = 1, num_interfaces
      do my_ipoin = 1, nibool_interfaces(my_iinterface)
        iglob = ibool_interfaces(my_ipoin,my_iinterface)
        array_val(iglob) = array_val(iglob) + mybuffer(my_ipoin,my_iinterface)
      enddo
    enddo
    need_add_my_contrib = .false.

    end subroutine add_my_contrib_scalar

  end subroutine assemble_MPI_scalar_w

!-------------------------------------------------------------------------------------------------
!
! blocking routines, ordered communication
!
!-------------------------------------------------------------------------------------------------

  subroutine assemble_MPI_scalar_block(array_val,nglob, &
                                       iproc_xi,iproc_eta,ichunk,addressing, &
                                       iboolleft_xi,iboolright_xi,iboolleft_eta,iboolright_eta, &
                                       npoin2D_faces,npoin2D_xi,npoin2D_eta, &
                                       iboolfaces,iboolcorner, &
                                       iprocfrom_faces,iprocto_faces,imsg_type, &
                                       iproc_main_corners,iproc_worker1_corners,iproc_worker2_corners, &
                                       buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
                                       buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
                                       NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
                                       NPROC_XI,NPROC_ETA,NGLOB1D_RADIAL, &
                                       NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB2DMAX_XY,NCHUNKS)

! this version of the routine is based on blocking MPI calls

  use constants, only: CUSTOM_REAL,NB_SQUARE_EDGES_ONEDIR,NUMFACES_SHARED,NUMCORNERS_SHARED, &
    ACTUALLY_ASSEMBLE_MPI_SLICES,ACTUALLY_ASSEMBLE_MPI_CHUNKS,itag,itag2,myrank

  implicit none

  integer,intent(in) :: nglob,NCHUNKS

! array to assemble
  real(kind=CUSTOM_REAL), dimension(nglob) :: array_val

  integer,intent(in) :: iproc_xi,iproc_eta,ichunk
  integer, dimension(NB_SQUARE_EDGES_ONEDIR),intent(in) :: npoin2D_xi,npoin2D_eta
  integer,intent(in) :: npoin2D_faces(NUMFACES_SHARED)

  integer,intent(in) :: NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB2DMAX_XY
  integer,intent(in) :: NPROC_XI,NPROC_ETA,NGLOB1D_RADIAL
  integer,intent(in) :: NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS

! for addressing of the slices
  integer, dimension(NCHUNKS,0:NPROC_XI-1,0:NPROC_ETA-1),intent(in) :: addressing

! 2-D addressing and buffers for summation between slices
  integer, dimension(NGLOB2DMAX_XMIN_XMAX),intent(in) :: iboolleft_xi,iboolright_xi
  integer, dimension(NGLOB2DMAX_YMIN_YMAX),intent(in) :: iboolleft_eta,iboolright_eta

! indirect addressing for each corner of the chunks
  integer, dimension(NGLOB1D_RADIAL,NUMCORNERS_SHARED),intent(in) :: iboolcorner

  integer,intent(in) :: npoin2D_max_all_CM_IC
  integer, dimension(NGLOB2DMAX_XY,NUMFACES_SHARED),intent(in) :: iboolfaces
  real(kind=CUSTOM_REAL), dimension(npoin2D_max_all_CM_IC) :: buffer_send_faces_scalar,buffer_received_faces_scalar

! buffers for send and receive between corners of the chunks
  real(kind=CUSTOM_REAL), dimension(NGLOB1D_RADIAL) :: buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar

! ---- arrays to assemble between chunks

! communication pattern for faces between chunks
  integer, dimension(NUMMSGS_FACES),intent(in) :: iprocfrom_faces,iprocto_faces,imsg_type

! communication pattern for corners between chunks
  integer, dimension(NCORNERSCHUNKS),intent(in) :: iproc_main_corners,iproc_worker1_corners,iproc_worker2_corners

  ! local parameters
  integer :: icount_corners
  integer :: ipoin,ipoin2D,ipoin1D
  integer :: sender,receiver
  integer :: imsg,imsg_loop
  integer :: icount_faces,npoin2D_chunks

  integer, external :: null_process

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
  if (NPROC_XI > 1) then

    ! slices copy the right face into the buffer
    do ipoin = 1,npoin2D_xi(2)
      buffer_send_faces_scalar(ipoin) = array_val(iboolright_xi(ipoin))
    enddo

    ! send messages forward along each row
    if (iproc_xi == 0) then
      sender = null_process()
    else
      sender = addressing(ichunk,iproc_xi - 1,iproc_eta)
    endif
    if (iproc_xi == NPROC_XI-1) then
      receiver = null_process()
    else
      receiver = addressing(ichunk,iproc_xi + 1,iproc_eta)
    endif
    call sendrecv_cr(buffer_send_faces_scalar,npoin2D_xi(2),receiver,itag2, &
                     buffer_received_faces_scalar,npoin2D_xi(1),sender,itag)

    ! all slices add the buffer received to the contributions on the left face
    if (iproc_xi > 0) then
      do ipoin = 1,npoin2D_xi(1)
        array_val(iboolleft_xi(ipoin)) = array_val(iboolleft_xi(ipoin)) + &
                                buffer_received_faces_scalar(ipoin)
      enddo
    endif

    ! the contributions are correctly assembled on the left side of each slice
    ! now we have to send the result back to the sender
    ! all slices copy the left face into the buffer
    do ipoin = 1,npoin2D_xi(1)
      buffer_send_faces_scalar(ipoin) = array_val(iboolleft_xi(ipoin))
    enddo

    ! send messages backward along each row
    if (iproc_xi == NPROC_XI-1) then
      sender = null_process()
    else
      sender = addressing(ichunk,iproc_xi + 1,iproc_eta)
    endif
    if (iproc_xi == 0) then
      receiver = null_process()
    else
      receiver = addressing(ichunk,iproc_xi - 1,iproc_eta)
    endif
    call sendrecv_cr(buffer_send_faces_scalar,npoin2D_xi(1),receiver,itag2, &
                     buffer_received_faces_scalar,npoin2D_xi(2),sender,itag)

    ! all slices copy the buffer received to the contributions on the right face
    if (iproc_xi < NPROC_XI-1) then
      do ipoin = 1,npoin2D_xi(2)
        array_val(iboolright_xi(ipoin)) = buffer_received_faces_scalar(ipoin)
      enddo
    endif

  endif

!----
!---- then assemble along eta using the 2-D topology
!----

! assemble along eta only if more than one slice
  if (NPROC_ETA > 1) then

    ! slices copy the right face into the buffer
    do ipoin = 1,npoin2D_eta(2)
      buffer_send_faces_scalar(ipoin) = array_val(iboolright_eta(ipoin))
    enddo

    ! send messages forward along each row
    if (iproc_eta == 0) then
      sender = null_process()
    else
      sender = addressing(ichunk,iproc_xi,iproc_eta - 1)
    endif
    if (iproc_eta == NPROC_ETA-1) then
      receiver = null_process()
    else
      receiver = addressing(ichunk,iproc_xi,iproc_eta + 1)
    endif
    call sendrecv_cr(buffer_send_faces_scalar,npoin2D_eta(2),receiver,itag2, &
                     buffer_received_faces_scalar,npoin2D_eta(1),sender,itag)

    ! all slices add the buffer received to the contributions on the left face
    if (iproc_eta > 0) then
      do ipoin = 1,npoin2D_eta(1)
        array_val(iboolleft_eta(ipoin)) = array_val(iboolleft_eta(ipoin)) + &
                                buffer_received_faces_scalar(ipoin)
      enddo
    endif

    ! the contributions are correctly assembled on the left side of each slice
    ! now we have to send the result back to the sender
    ! all slices copy the left face into the buffer
    do ipoin = 1,npoin2D_eta(1)
      buffer_send_faces_scalar(ipoin) = array_val(iboolleft_eta(ipoin))
    enddo

    ! send messages backward along each row
    if (iproc_eta == NPROC_ETA-1) then
      sender = null_process()
    else
      sender = addressing(ichunk,iproc_xi,iproc_eta + 1)
    endif
    if (iproc_eta == 0) then
      receiver = null_process()
    else
      receiver = addressing(ichunk,iproc_xi,iproc_eta - 1)
    endif
    call sendrecv_cr(buffer_send_faces_scalar,npoin2D_eta(1),receiver,itag2, &
                     buffer_received_faces_scalar,npoin2D_eta(2),sender,itag)

    ! all slices copy the buffer received to the contributions on the right face
    if (iproc_eta < NPROC_ETA-1) then
      do ipoin = 1,npoin2D_eta(2)
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
      if (myrank == iprocfrom_faces(imsg) .or. &
           myrank == iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if (myrank == iprocto_faces(imsg) .and. imsg_type(imsg) == imsg_loop) then
        sender = iprocfrom_faces(imsg)
        npoin2D_chunks = npoin2D_faces(icount_faces)
        call recv_cr(buffer_received_faces_scalar,npoin2D_chunks,sender,itag)

        do ipoin2D = 1,npoin2D_chunks
          array_val(iboolfaces(ipoin2D,icount_faces)) = &
             array_val(iboolfaces(ipoin2D,icount_faces)) + buffer_received_faces_scalar(ipoin2D)
        enddo
      endif
    enddo

    !---- put slices in send mode
    !---- a given slice can belong to at most two faces
    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if (myrank == iprocfrom_faces(imsg) .or. &
           myrank == iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if (myrank == iprocfrom_faces(imsg) .and. imsg_type(imsg) == imsg_loop) then
        receiver = iprocto_faces(imsg)
        npoin2D_chunks = npoin2D_faces(icount_faces)
        do ipoin2D = 1,npoin2D_chunks
          buffer_send_faces_scalar(ipoin2D) = array_val(iboolfaces(ipoin2D,icount_faces))
        enddo
        call send_cr(buffer_send_faces_scalar,npoin2D_chunks,receiver,itag)
      endif
    enddo

    ! *********************************************************************
    !  transmit messages back in opposite direction (iprocto -> iprocfrom)
    ! *********************************************************************

    !---- put slices in receive mode
    !---- a given slice can belong to at most two faces

    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if (myrank == iprocfrom_faces(imsg) .or. &
           myrank == iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if (myrank == iprocfrom_faces(imsg) .and. imsg_type(imsg) == imsg_loop) then
        sender = iprocto_faces(imsg)
        npoin2D_chunks = npoin2D_faces(icount_faces)
        call recv_cr(buffer_received_faces_scalar,npoin2D_chunks,sender,itag)

        do ipoin2D = 1,npoin2D_chunks
          array_val(iboolfaces(ipoin2D,icount_faces)) = buffer_received_faces_scalar(ipoin2D)
        enddo
      endif
    enddo

    !---- put slices in send mode
    !---- a given slice can belong to at most two faces
    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if (myrank == iprocfrom_faces(imsg) .or. &
           myrank == iprocto_faces(imsg)) icount_faces = icount_faces + 1
      if (myrank == iprocto_faces(imsg) .and. imsg_type(imsg) == imsg_loop) then
        receiver = iprocfrom_faces(imsg)
        npoin2D_chunks = npoin2D_faces(icount_faces)
        do ipoin2D = 1,npoin2D_chunks
          buffer_send_faces_scalar(ipoin2D) = array_val(iboolfaces(ipoin2D,icount_faces))
        enddo
        call send_cr(buffer_send_faces_scalar,npoin2D_chunks,receiver,itag)
      endif
    enddo

! end of anti-deadlocking loop
  enddo

!----
!---- start MPI assembling corners
!----

! scheme for corners cannot deadlock even if NPROC_XI = NPROC_ETA = 1

! ***************************************************************
!  transmit messages in forward direction (two workers -> main process)
! ***************************************************************

  icount_corners = 0

  do imsg = 1,NCORNERSCHUNKS

    if (myrank == iproc_main_corners(imsg) .or. &
       myrank == iproc_worker1_corners(imsg) .or. &
       (NCHUNKS /= 2 .and. myrank == iproc_worker2_corners(imsg))) icount_corners = icount_corners + 1

    !---- receive messages from the two workers on the main
    if (myrank == iproc_main_corners(imsg)) then

      ! receive from worker #1 and add to local array
      sender = iproc_worker1_corners(imsg)
      call recv_cr(buffer_recv_chunkcorn_scalar,NGLOB1D_RADIAL,sender,itag)

      do ipoin1D = 1,NGLOB1D_RADIAL
        array_val(iboolcorner(ipoin1D,icount_corners)) = array_val(iboolcorner(ipoin1D,icount_corners)) + &
                 buffer_recv_chunkcorn_scalar(ipoin1D)
      enddo

      ! receive from worker #2 and add to local array
      if (NCHUNKS /= 2) then
        sender = iproc_worker2_corners(imsg)
        call recv_cr(buffer_recv_chunkcorn_scalar,NGLOB1D_RADIAL,sender,itag)

        do ipoin1D = 1,NGLOB1D_RADIAL
          array_val(iboolcorner(ipoin1D,icount_corners)) = array_val(iboolcorner(ipoin1D,icount_corners)) + &
                   buffer_recv_chunkcorn_scalar(ipoin1D)
        enddo
      endif

    endif

    !---- send messages from the two workers to the main
    if (myrank == iproc_worker1_corners(imsg) .or. &
                (NCHUNKS /= 2 .and. myrank == iproc_worker2_corners(imsg))) then

      receiver = iproc_main_corners(imsg)
      do ipoin1D = 1,NGLOB1D_RADIAL
        buffer_send_chunkcorn_scalar(ipoin1D) = array_val(iboolcorner(ipoin1D,icount_corners))
      enddo
      call send_cr(buffer_send_chunkcorn_scalar,NGLOB1D_RADIAL,receiver,itag)

    endif

    ! *********************************************************************
    !  transmit messages back in opposite direction (main process -> two workers)
    ! *********************************************************************

    !---- receive messages from the main on the two workers
    if (myrank == iproc_worker1_corners(imsg) .or. &
                (NCHUNKS /= 2 .and. myrank == iproc_worker2_corners(imsg))) then

      ! receive from main and copy to local array
      sender = iproc_main_corners(imsg)
      call recv_cr(buffer_recv_chunkcorn_scalar,NGLOB1D_RADIAL,sender,itag)

      do ipoin1D = 1,NGLOB1D_RADIAL
        array_val(iboolcorner(ipoin1D,icount_corners)) = buffer_recv_chunkcorn_scalar(ipoin1D)
      enddo

    endif

    !---- send messages from the main to the two workers
    if (myrank == iproc_main_corners(imsg)) then

      do ipoin1D = 1,NGLOB1D_RADIAL
        buffer_send_chunkcorn_scalar(ipoin1D) = array_val(iboolcorner(ipoin1D,icount_corners))
      enddo

      ! send to worker #1
      receiver = iproc_worker1_corners(imsg)
      call send_cr(buffer_send_chunkcorn_scalar,NGLOB1D_RADIAL,receiver,itag)

      ! send to worker #2
      if (NCHUNKS /= 2) then
        receiver = iproc_worker2_corners(imsg)
        call send_cr(buffer_send_chunkcorn_scalar,NGLOB1D_RADIAL,receiver,itag)
      endif

    endif

  enddo

  end subroutine assemble_MPI_scalar_block
