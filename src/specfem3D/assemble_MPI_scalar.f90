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

!----
!---- assemble the contributions between slices and chunks using MPI
!----

  subroutine assemble_MPI_scalar(NPROC,nglob,array_val, &
                        num_interfaces,max_nibool_interfaces, &
                        nibool_interfaces,ibool_interfaces, &
                        my_neighbours)

! blocking send/receive

  use constants

  implicit none

  integer :: NPROC
  integer :: nglob

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(nglob), intent(inout) :: array_val

  integer :: num_interfaces,max_nibool_interfaces
  integer, dimension(num_interfaces) :: nibool_interfaces,my_neighbours
  integer, dimension(max_nibool_interfaces,num_interfaces) :: ibool_interfaces

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_send_scalar
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_recv_scalar
  integer, dimension(:), allocatable :: request_send_scalar
  integer, dimension(:), allocatable :: request_recv_scalar

  integer ipoin,iinterface,ier

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
                   nibool_interfaces(iinterface),my_neighbours(iinterface), &
                   itag,request_send_scalar(iinterface) )
      ! receive request
      call irecv_cr(buffer_recv_scalar(1:nibool_interfaces(iinterface),iinterface), &
                   nibool_interfaces(iinterface),my_neighbours(iinterface), &
                   itag,request_recv_scalar(iinterface) )
    enddo

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces
      call wait_req(request_recv_scalar(iinterface))
    enddo

    ! adding contributions of neighbours
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

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_scalar_s(NPROC,nglob,array_val, &
                        buffer_send_scalar,buffer_recv_scalar, &
                        num_interfaces,max_nibool_interfaces, &
                        nibool_interfaces,ibool_interfaces, &
                        my_neighbours, &
                        request_send_scalar,request_recv_scalar)

! non-blocking MPI send

  use constants

  implicit none

  integer :: NPROC
  integer :: nglob
  integer :: num_interfaces,max_nibool_interfaces

! array to send
  real(kind=CUSTOM_REAL), dimension(nglob) :: array_val


  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces,num_interfaces) :: &
       buffer_send_scalar,buffer_recv_scalar

  integer, dimension(num_interfaces) :: nibool_interfaces,my_neighbours
  integer, dimension(max_nibool_interfaces,num_interfaces) :: ibool_interfaces
  integer, dimension(num_interfaces) :: request_send_scalar,request_recv_scalar

  integer ipoin,iinterface

! sends only if more than one partition
  if (NPROC > 1) then

    ! partition border copy into the buffer
    do iinterface = 1, num_interfaces
      do ipoin = 1, nibool_interfaces(iinterface)
        buffer_send_scalar(ipoin,iinterface) = &
          array_val(ibool_interfaces(ipoin,iinterface))
      enddo
    enddo

    ! send messages
    do iinterface = 1, num_interfaces
      ! non-blocking synchronous send request
      call isend_cr(buffer_send_scalar(1:nibool_interfaces(iinterface),iinterface), &
                   nibool_interfaces(iinterface),my_neighbours(iinterface), &
                   itag,request_send_scalar(iinterface))
      ! receive request
      call irecv_cr(buffer_recv_scalar(1:nibool_interfaces(iinterface),iinterface), &
                   nibool_interfaces(iinterface),my_neighbours(iinterface), &
                   itag,request_recv_scalar(iinterface))

    enddo

  endif

  end subroutine assemble_MPI_scalar_s

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_scalar_w(NPROC,nglob,array_val, &
                        buffer_recv_scalar,num_interfaces, &
                        max_nibool_interfaces, &
                        nibool_interfaces,ibool_interfaces, &
                        request_send_scalar,request_recv_scalar)

! waits for send/receiver to be completed and assembles contributions

  use constants

  implicit none

  integer :: NPROC
  integer :: nglob
  integer :: num_interfaces,max_nibool_interfaces
! array to assemble
  real(kind=CUSTOM_REAL), dimension(nglob) :: array_val


  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces,num_interfaces) :: &
       buffer_recv_scalar

  integer, dimension(num_interfaces) :: nibool_interfaces
  integer, dimension(max_nibool_interfaces,num_interfaces) :: ibool_interfaces
  integer, dimension(num_interfaces) :: request_send_scalar,request_recv_scalar

  integer ipoin,iinterface

! assemble only if more than one partition
  if (NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces
      call wait_req(request_recv_scalar(iinterface))
    enddo

    ! adding contributions of neighbours
    do iinterface = 1, num_interfaces
      do ipoin = 1, nibool_interfaces(iinterface)
        array_val(ibool_interfaces(ipoin,iinterface)) = &
             array_val(ibool_interfaces(ipoin,iinterface)) &
             + buffer_recv_scalar(ipoin,iinterface)
      enddo
    enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces
      call wait_req(request_send_scalar(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_scalar_w


!-------------------------------------------------------------------------------------------------
!
! CUDA routines
!
!-------------------------------------------------------------------------------------------------


  subroutine assemble_MPI_scalar_send_gpu(NPROC, &
                                          buffer_send_scalar,buffer_recv_scalar, &
                                          num_interfaces,max_nibool_interfaces, &
                                          nibool_interfaces, &
                                          my_neighbours, &
                                          request_send_scalar,request_recv_scalar)

! non-blocking MPI send

  ! sends data
  ! note: assembling data already filled into buffer_send_scalar array

  use constants

  implicit none

  integer :: NPROC
  integer :: num_interfaces,max_nibool_interfaces

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces,num_interfaces) :: &
       buffer_send_scalar,buffer_recv_scalar

  integer, dimension(num_interfaces) :: nibool_interfaces,my_neighbours
  integer, dimension(num_interfaces) :: request_send_scalar,request_recv_scalar

  ! local parameters
  integer iinterface

  ! note: preparation of the contribution between partitions using MPI
  !       transfers MPI buffers to CPU in already done in transfer_boun_pot_from_device() call

  ! sends only if more than one partition
  if (NPROC > 1) then

    ! send messages
    do iinterface = 1, num_interfaces
      ! non-blocking synchronous send request
      call isend_cr(buffer_send_scalar(1:nibool_interfaces(iinterface),iinterface), &
                   nibool_interfaces(iinterface),my_neighbours(iinterface), &
                   itag,request_send_scalar(iinterface))
      ! receive request
      call irecv_cr(buffer_recv_scalar(1:nibool_interfaces(iinterface),iinterface), &
                   nibool_interfaces(iinterface),my_neighbours(iinterface), &
                   itag,request_recv_scalar(iinterface) )

    enddo

  endif

  end subroutine assemble_MPI_scalar_send_gpu

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_scalar_write_gpu(Mesh_pointer,NPROC, &
                                            buffer_recv_scalar, &
                                            num_interfaces,max_nibool_interfaces, &
                                            request_send_scalar,request_recv_scalar, &
                                            FORWARD_OR_ADJOINT)

! waits for send/receiver to be completed and assembles contributions

  use constants

  implicit none

  integer(kind=8) :: Mesh_pointer

  integer :: NPROC

  integer :: num_interfaces,max_nibool_interfaces

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces,num_interfaces) :: buffer_recv_scalar
  integer, dimension(num_interfaces) :: request_send_scalar,request_recv_scalar

  integer :: FORWARD_OR_ADJOINT

  ! local parameters
  integer :: iinterface

! assemble only if more than one partition
  if (NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces
      call wait_req(request_recv_scalar(iinterface))
    enddo

    ! adding contributions of neighbours
    call transfer_asmbl_pot_to_device(Mesh_pointer,buffer_recv_scalar,FORWARD_OR_ADJOINT)

    ! note: adding contributions of neighbours has been done just above for gpu
    !do iinterface = 1, num_interfaces
    !  do ipoin = 1, nibool_interfaces(iinterface)
    !    array_val(ibool_interfaces(ipoin,iinterface)) = &
    !         array_val(ibool_interfaces(ipoin,iinterface)) &
    !         + buffer_recv_scalar(ipoin,iinterface)
    !  enddo
    !enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces
      call wait_req(request_send_scalar(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_scalar_write_gpu

!
!-------------------------------------------------------------------------------------------------
!

! with gpu functions...

  subroutine transfer_boundarypot_to_device(Mesh_pointer, NPROC, &
                                            buffer_recv_scalar, &
                                            num_interfaces,max_nibool_interfaces,&
                                            request_recv_scalar, &
                                            IREGION,FORWARD_OR_ADJOINT)

  use constants

  implicit none

  integer(kind=8) :: Mesh_pointer

  integer :: NPROC

  ! array to assemble
  integer :: num_interfaces,max_nibool_interfaces

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces,num_interfaces) :: buffer_recv_scalar

  integer, dimension(num_interfaces) :: request_recv_scalar

  integer :: IREGION
  integer :: FORWARD_OR_ADJOINT

  ! local parameters
  integer :: iinterface

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! waits for communications completion (recv)
    do iinterface = 1, num_interfaces
      call wait_req(request_recv_scalar(iinterface))
    enddo

    ! sends contributions to GPU
    call transfer_buffer_to_device_async(Mesh_pointer, &
                                         buffer_recv_scalar, &
                                         IREGION,FORWARD_OR_ADJOINT)
  endif

  ! This step is done via previous function transfer_and_assemble...
  ! do iinterface = 1, num_interfaces
  !   do ipoin = 1, nibool_interfaces(iinterface)
  !     array_val(ibool_interfaces(ipoin,iinterface)) = &
  !          array_val(ibool_interfaces(ipoin,iinterface)) + buffer_recv_scalar(ipoin,iinterface)
  !   enddo
  ! enddo

  end subroutine transfer_boundarypot_to_device


