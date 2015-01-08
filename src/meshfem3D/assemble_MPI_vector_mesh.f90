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


  subroutine assemble_MPI_vector(NPROC,NGLOB_AB,array_val, &
                        num_interfaces,max_nibool_interfaces, &
                        nibool_interfaces,ibool_interfaces, &
                        my_neighbours)

  use constants

  implicit none

  integer :: NPROC
  integer :: NGLOB_AB

! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: array_val

  integer :: num_interfaces,max_nibool_interfaces
  integer, dimension(num_interfaces) :: nibool_interfaces,my_neighbours
  integer, dimension(max_nibool_interfaces,num_interfaces) :: ibool_interfaces

  ! local parameters

  ! send/receive temporary buffers
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_recv_vector

  ! requests
  integer, dimension(:), allocatable :: request_send_vector
  integer, dimension(:), allocatable :: request_recv_vector

  integer ipoin,iinterface,ier


! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC > 1) then

    allocate(buffer_send_vector(NDIM,max_nibool_interfaces,num_interfaces),stat=ier)
    if (ier /= 0 ) stop 'Error allocating array buffer_send_vector'
    allocate(buffer_recv_vector(NDIM,max_nibool_interfaces,num_interfaces),stat=ier)
    if (ier /= 0 ) stop 'Error allocating array buffer_recv_vector'
    allocate(request_send_vector(num_interfaces),stat=ier)
    if (ier /= 0 ) stop 'Error allocating array request_send_vector'
    allocate(request_recv_vector(num_interfaces),stat=ier)
    if (ier /= 0 ) stop 'Error allocating array request_recv_vector'

    ! partition border copy into the buffer
    do iinterface = 1, num_interfaces
      do ipoin = 1, nibool_interfaces(iinterface)
        buffer_send_vector(:,ipoin,iinterface) = &
          array_val(:,ibool_interfaces(ipoin,iinterface))
      enddo
    enddo

    ! send messages
    do iinterface = 1, num_interfaces
      call isend_cr(buffer_send_vector(1,1,iinterface), &
           NDIM*nibool_interfaces(iinterface), &
           my_neighbours(iinterface), &
           itag, &
           request_send_vector(iinterface) &
           )
      call irecv_cr(buffer_recv_vector(1,1,iinterface), &
           NDIM*nibool_interfaces(iinterface), &
           my_neighbours(iinterface), &
           itag, &
           request_recv_vector(iinterface) &
           )
    enddo

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces
      call wait_req(request_recv_vector(iinterface))
    enddo

    ! adding contributions of neighbours
    do iinterface = 1, num_interfaces
      do ipoin = 1, nibool_interfaces(iinterface)
        array_val(:,ibool_interfaces(ipoin,iinterface)) = &
             array_val(:,ibool_interfaces(ipoin,iinterface)) &
             + buffer_recv_vector(:,ipoin,iinterface)
      enddo
    enddo

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces
      call wait_req(request_send_vector(iinterface))
    enddo

    deallocate(buffer_send_vector)
    deallocate(buffer_recv_vector)
    deallocate(request_send_vector)
    deallocate(request_recv_vector)

  endif

  end subroutine assemble_MPI_vector
