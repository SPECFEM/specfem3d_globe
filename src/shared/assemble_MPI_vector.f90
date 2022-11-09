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

!----
!---- assemble the contributions between slices and chunks using MPI
!----


  subroutine assemble_MPI_vector(NPROC,NGLOB_AB,array_val, &
                                 num_interfaces,max_nibool_interfaces, &
                                 nibool_interfaces,ibool_interfaces, &
                                 my_neighbors)

  use constants, only: CUSTOM_REAL,NDIM,itag

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: NGLOB_AB

! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_AB) :: array_val

  integer,intent(in) :: num_interfaces,max_nibool_interfaces
  integer, dimension(num_interfaces),intent(in) :: nibool_interfaces,my_neighbors
  integer, dimension(max_nibool_interfaces,num_interfaces),intent(in) :: ibool_interfaces

  ! local parameters

  ! send/receive temporary buffers
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_recv_vector

  ! requests
  integer, dimension(:), allocatable :: request_send_vector
  integer, dimension(:), allocatable :: request_recv_vector

  integer :: ipoin,iinterface,ier


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
        buffer_send_vector(:,ipoin,iinterface) = array_val(:,ibool_interfaces(ipoin,iinterface))
      enddo
    enddo

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
                    request_recv_vector(iinterface) )
    enddo

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces
      call wait_req(request_recv_vector(iinterface))
    enddo

    ! adding contributions of neighbors
    do iinterface = 1, num_interfaces
      do ipoin = 1, nibool_interfaces(iinterface)
        array_val(:,ibool_interfaces(ipoin,iinterface)) = &
             array_val(:,ibool_interfaces(ipoin,iinterface)) + buffer_recv_vector(:,ipoin,iinterface)
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


!-------------------------------------------------------------------------------------------------
!
! non-blocking routines
!
!-------------------------------------------------------------------------------------------------

  subroutine assemble_MPI_vector_s(NPROC,nglob, &
                                   array_val, &
                                   buffer_send_vector,buffer_recv_vector, &
                                   num_interfaces,max_nibool_interfaces, &
                                   nibool_interfaces,ibool_interfaces, &
                                   my_neighbors, &
                                   request_send_vector,request_recv_vector)

! sends data

  use constants, only: CUSTOM_REAL,NDIM,itag

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: nglob

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob) :: array_val

  integer,intent(in) :: num_interfaces,max_nibool_interfaces

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces,num_interfaces) :: buffer_send_vector,buffer_recv_vector

  integer, dimension(num_interfaces),intent(in) :: nibool_interfaces,my_neighbors
  integer, dimension(max_nibool_interfaces,num_interfaces),intent(in) :: ibool_interfaces
  integer, dimension(num_interfaces) :: request_send_vector,request_recv_vector

  ! local parameters
  integer :: ipoin,iinterface

  ! here we have to assemble all the contributions between partitions using MPI

  ! assemble only if more than one partition
  if (NPROC > 1) then

     ! partition border copy into the buffer
     do iinterface = 1, num_interfaces
        do ipoin = 1, nibool_interfaces(iinterface)
           buffer_send_vector(:,ipoin,iinterface) = array_val(:,ibool_interfaces(ipoin,iinterface))
        enddo
     enddo

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

  end subroutine assemble_MPI_vector_s

!
!-------------------------------------------------------------------------------------------------
!

  subroutine assemble_MPI_vector_w(NPROC,nglob, &
                                   array_val, &
                                   buffer_recv_vector, &
                                   num_interfaces,max_nibool_interfaces, &
                                   nibool_interfaces,ibool_interfaces, &
                                   my_neighbors, &
                                   request_send_vector,request_recv_vector)

! waits for data to receive and assembles

  use constants, only: myrank,CUSTOM_REAL,NDIM,itag,DO_ORDERED_ASSEMBLY

  implicit none

  integer,intent(in) :: NPROC
  integer,intent(in) :: nglob

  ! array to assemble
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob) :: array_val

  integer,intent(in) :: num_interfaces,max_nibool_interfaces

  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces,num_interfaces) :: buffer_recv_vector

  integer, dimension(num_interfaces),intent(in) :: nibool_interfaces
  integer, dimension(max_nibool_interfaces,num_interfaces),intent(in) :: ibool_interfaces
  integer, dimension(num_interfaces) :: request_send_vector,request_recv_vector

  integer, dimension(num_interfaces),intent(in) :: my_neighbors

  ! local parameters
  integer :: ipoin,iinterface
  integer :: iglob
  ! ordered assembly
  real(kind=CUSTOM_REAL), dimension(NDIM,max_nibool_interfaces,num_interfaces) :: mybuffer
  logical :: need_add_my_contrib

! here we have to assemble all the contributions between partitions using MPI

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
          mybuffer(:,ipoin,iinterface) = array_val(:,iglob)
          array_val(:,iglob) = 0._CUSTOM_REAL
        enddo
      enddo
    endif

    ! wait for communications completion (recv)
    do iinterface = 1, num_interfaces
      call wait_req(request_recv_vector(iinterface))
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
        if (need_add_my_contrib .and. myrank < my_neighbors(iinterface)) call add_my_contrib()
        do ipoin = 1, nibool_interfaces(iinterface)
          iglob = ibool_interfaces(ipoin,iinterface)
          array_val(:,iglob) = array_val(:,iglob) + buffer_recv_vector(:,ipoin,iinterface)
        enddo
      enddo
      if (need_add_my_contrib) call add_my_contrib()

    else
      ! default assembly
      do iinterface = 1, num_interfaces
        do ipoin = 1, nibool_interfaces(iinterface)
          iglob = ibool_interfaces(ipoin,iinterface)
          array_val(:,iglob) = array_val(:,iglob) + buffer_recv_vector(:,ipoin,iinterface)
        enddo
      enddo
    endif

    ! wait for communications completion (send)
    do iinterface = 1, num_interfaces
      call wait_req(request_send_vector(iinterface))
    enddo

  endif

  contains

    subroutine add_my_contrib()

    integer :: my_iinterface,my_ipoin

    do my_iinterface = 1, num_interfaces
      do my_ipoin = 1, nibool_interfaces(my_iinterface)
        iglob = ibool_interfaces(my_ipoin,my_iinterface)
        array_val(:,iglob) = array_val(:,iglob) + mybuffer(:,my_ipoin,my_iinterface)
      enddo
    enddo
    need_add_my_contrib = .false.

    end subroutine add_my_contrib

  end subroutine assemble_MPI_vector_w
