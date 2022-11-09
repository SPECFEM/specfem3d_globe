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


! assemble the contributions between slices and chunks using MPI
! with GPU wrapper routines

!-------------------------------------------------------------------------------------------------
!
! CUDA routines - scalar field assembly
!
!-------------------------------------------------------------------------------------------------

  subroutine assemble_MPI_scalar_send_gpu(NPROC, &
                                          buffer_send_scalar,buffer_recv_scalar, &
                                          num_interfaces,max_nibool_interfaces, &
                                          nibool_interfaces, &
                                          my_neighbors, &
                                          request_send_scalar,request_recv_scalar)

! non-blocking MPI send

  ! sends data
  ! note: assembling data already filled into buffer_send_scalar array

  use constants, only: CUSTOM_REAL,itag

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: num_interfaces,max_nibool_interfaces

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces,num_interfaces),intent(inout) :: &
       buffer_send_scalar,buffer_recv_scalar

  integer, dimension(num_interfaces),intent(in) :: nibool_interfaces,my_neighbors
  integer, dimension(num_interfaces),intent(inout) :: request_send_scalar,request_recv_scalar

  ! local parameters
  integer :: iinterface

  ! note: preparation of the contribution between partitions using MPI
  !       transfers MPI buffers to CPU is already done in transfer_boun_pot_from_device() call

  ! sends only if more than one partition
  if (NPROC > 1) then

    ! send messages
    do iinterface = 1, num_interfaces
      ! non-blocking synchronous send request
      call isend_cr(buffer_send_scalar(1:nibool_interfaces(iinterface),iinterface), &
                    nibool_interfaces(iinterface),my_neighbors(iinterface), &
                    itag,request_send_scalar(iinterface))
      ! receive request
      call irecv_cr(buffer_recv_scalar(1:nibool_interfaces(iinterface),iinterface), &
                    nibool_interfaces(iinterface),my_neighbors(iinterface), &
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

  use constants, only: CUSTOM_REAL

  implicit none

  integer(kind=8),intent(in) :: Mesh_pointer

  integer,intent(in) :: NPROC

  integer,intent(in) :: num_interfaces,max_nibool_interfaces

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces,num_interfaces),intent(inout) :: buffer_recv_scalar
  integer, dimension(num_interfaces),intent(in) :: request_send_scalar,request_recv_scalar

  integer,intent(in) :: FORWARD_OR_ADJOINT

  ! local parameters
  integer :: iinterface

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces
      call wait_req(request_recv_scalar(iinterface))
    enddo

    ! adding contributions of neighbors
    call transfer_asmbl_pot_to_device(Mesh_pointer,buffer_recv_scalar,FORWARD_OR_ADJOINT)

    ! note: adding contributions of neighbors has been done just above for gpu
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
                                            num_interfaces,max_nibool_interfaces, &
                                            request_recv_scalar, &
                                            IREGION,FORWARD_OR_ADJOINT)

  use constants, only: CUSTOM_REAL

  implicit none

  integer(kind=8),intent(in) :: Mesh_pointer

  integer,intent(in) :: NPROC

  ! array to assemble
  integer,intent(in) :: num_interfaces,max_nibool_interfaces

  real(kind=CUSTOM_REAL), dimension(max_nibool_interfaces,num_interfaces),intent(inout) :: buffer_recv_scalar

  integer, dimension(num_interfaces),intent(in) :: request_recv_scalar

  integer,intent(in) :: IREGION
  integer,intent(in) :: FORWARD_OR_ADJOINT

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



!-------------------------------------------------------------------------------------------------
!
! CUDA routines  - vector field assembly
!
!-------------------------------------------------------------------------------------------------

! interrupt might improve MPI performance
! see: https://computing.llnl.gov/tutorials/mpi_performance/#Sender-ReceiverSync
!
! check: MP_CSS_INTERRUPT environment variable on IBM systems

  subroutine assemble_MPI_vector_send_gpu(NPROC, &
                                          buffer_send_vector,buffer_recv_vector, &
                                          num_interfaces,max_nibool_interfaces, &
                                          nibool_interfaces, &
                                          my_neighbors, &
                                          request_send_vector,request_recv_vector)

  ! sends data
  ! note: array to assemble already filled into buffer_send_vector array

  use constants, only: CUSTOM_REAL,NDIM,itag

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: num_interfaces,max_nibool_interfaces

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces,num_interfaces),intent(inout) :: &
       buffer_send_vector,buffer_recv_vector

  integer, dimension(num_interfaces),intent(in) :: nibool_interfaces,my_neighbors
  integer, dimension(num_interfaces),intent(inout) :: request_send_vector,request_recv_vector

  ! local parameters
  integer :: iinterface

  ! note: preparation of the contribution between partitions using MPI
  !          already done in transfer_boun_from_device() routines

  ! send only if more than one partition
  if (NPROC > 1) then

    ! send messages
    do iinterface = 1, num_interfaces
      call isend_cr(buffer_send_vector(1,1,iinterface), &
                    NDIM*nibool_interfaces(iinterface), &
                    my_neighbors(iinterface), &
                    itag, &
                    request_send_vector(iinterface))

      call irecv_cr(buffer_recv_vector(1,1,iinterface), &
                    NDIM*nibool_interfaces(iinterface), &
                    my_neighbors(iinterface), &
                    itag, &
                    request_recv_vector(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_vector_send_gpu

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_vector_write_gpu(Mesh_pointer,NPROC, &
                                            buffer_recv_vector, &
                                            num_interfaces,max_nibool_interfaces, &
                                            request_send_vector,request_recv_vector, &
                                            IREGION,FORWARD_OR_ADJOINT )

! waits for data to receive and assembles

  use constants, only: CUSTOM_REAL,NDIM

  implicit none

  integer(kind=8),intent(in) :: Mesh_pointer

  integer,intent(in) :: NPROC

  integer,intent(in) :: num_interfaces,max_nibool_interfaces

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces,num_interfaces),intent(inout) :: buffer_recv_vector
  integer, dimension(num_interfaces),intent(in) :: request_send_vector,request_recv_vector

  integer,intent(in) :: IREGION
  integer,intent(in) :: FORWARD_OR_ADJOINT

  ! local parameters
  integer :: iinterface

! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces
      call wait_req(request_recv_vector(iinterface))
    enddo

    ! adding contributions of neighbors
    call transfer_asmbl_accel_to_device(Mesh_pointer, &
                                        buffer_recv_vector, &
                                        IREGION,FORWARD_OR_ADJOINT)

    ! This step is done via previous function transfer_and_assemble...
    ! do iinterface = 1, num_interfaces
    !   do ipoin = 1, nibool_interfaces(iinterface)
    !     array_val(:,ibool_interfaces(ipoin,iinterface)) = &
    !          array_val(:,ibool_interfaces(ipoin,iinterface)) + buffer_recv_vector(:,ipoin,iinterface)
    !   enddo
    ! enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces
      call wait_req(request_send_vector(iinterface))
    enddo

  endif

  end subroutine assemble_MPI_vector_write_gpu


!
!-------------------------------------------------------------------------------------------------
!

! with gpu functions...

  subroutine transfer_boundary_to_device(Mesh_pointer, NPROC, &
                                         buffer_recv_vector, &
                                         num_interfaces,max_nibool_interfaces, &
                                         request_recv_vector, &
                                         IREGION,FORWARD_OR_ADJOINT)

  use constants, only: CUSTOM_REAL,NDIM

  implicit none

  integer(kind=8),intent(in) :: Mesh_pointer

  integer,intent(in) :: NPROC

  ! array to assemble
  integer,intent(in) :: num_interfaces,max_nibool_interfaces

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces,num_interfaces),intent(inout) :: buffer_recv_vector

  integer, dimension(num_interfaces),intent(in) :: request_recv_vector

  integer,intent(in) :: IREGION
  integer,intent(in) :: FORWARD_OR_ADJOINT

  ! local parameters
  integer :: iinterface

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC > 1) then

    ! waits for communications completion (recv)
    do iinterface = 1, num_interfaces
      call wait_req(request_recv_vector(iinterface))
    enddo

    ! sends contributions to GPU
    call transfer_buffer_to_device_async(Mesh_pointer, &
                                         buffer_recv_vector, &
                                         IREGION,FORWARD_OR_ADJOINT)
  endif

  ! This step is done via previous function transfer_and_assemble...
  ! do iinterface = 1, num_interfaces
  !   do ipoin = 1, nibool_interfaces(iinterface)
  !     array_val(:,ibool_interfaces(ipoin,iinterface)) = &
  !          array_val(:,ibool_interfaces(ipoin,iinterface)) + buffer_recv_vector(:,ipoin,iinterface)
  !   enddo
  ! enddo

  end subroutine transfer_boundary_to_device



