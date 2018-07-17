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

  subroutine get_MPI_interfaces(myrank,NGLOB,NSPEC, &
                                    test_flag,my_neighbors,nibool_neighbors,ibool_neighbors, &
                                    num_interfaces,max_nibool_interfaces, &
                                    max_nibool,MAX_NEIGHBORS, &
                                    ibool, &
                                    is_on_a_slice_edge, &
                                    IREGION,add_central_cube,idoubling,INCLUDE_CENTRAL_CUBE, &
                                    xstore,ystore,zstore,NPROCTOT)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,IREGION_INNER_CORE,IFLAG_IN_FICTITIOUS_CUBE
  implicit none

  integer,intent(in) :: myrank,NGLOB,NSPEC

  real(kind=CUSTOM_REAL),dimension(NGLOB),intent(in) :: test_flag

  integer,intent(in) :: max_nibool,MAX_NEIGHBORS
  integer, dimension(MAX_NEIGHBORS),intent(inout) :: my_neighbors,nibool_neighbors
  integer, dimension(max_nibool,MAX_NEIGHBORS),intent(inout) :: ibool_neighbors

  integer,intent(inout) :: num_interfaces,max_nibool_interfaces

  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: ibool

  logical,dimension(NSPEC),intent(inout) :: is_on_a_slice_edge

  integer,intent(in) :: IREGION
  logical,intent(in) :: add_central_cube
  integer,dimension(NSPEC),intent(in) :: idoubling

  logical,intent(in) :: INCLUDE_CENTRAL_CUBE

  real(kind=CUSTOM_REAL),dimension(NGLOB),intent(in) :: xstore,ystore,zstore

  integer :: NPROCTOT

  ! local parameters
  integer :: ispec,iglob,j,k
  integer :: iface,iedge,icorner
  integer :: ii,iinterface,icurrent,rank
  integer :: npoin
  logical :: is_done,ispec_is_outer
  integer,dimension(NGLOB) :: work_test_flag
  logical,dimension(NSPEC) :: work_ispec_is_outer

  integer,parameter :: MID = (NGLLX+1)/2

  ! initializes
  if (add_central_cube) then
    ! adds points to existing inner_core interfaces
    iinterface = num_interfaces
    work_ispec_is_outer(:) = is_on_a_slice_edge(:)
  else
    ! creates new interfaces
    iinterface = 0
    num_interfaces = 0
    max_nibool_interfaces = 0
    my_neighbors(:) = -1
    nibool_neighbors(:) = 0
    ibool_neighbors(:,:) = 0
    work_ispec_is_outer(:) = .false.
  endif

  ! makes working copy (converted to nearest integers)
  work_test_flag(:) = nint( test_flag(:) )

  ! loops over all elements
  do ispec = 1,NSPEC

    ! exclude elements in inner part of slice
    !if (.not. is_on_a_slice_edge(ispec) ) cycle

    ! exclude elements in fictitious core
    if (IREGION == IREGION_INNER_CORE) then
      if (idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE ) cycle
    endif

    ! sets flag if element has global points shared with other processes
    ispec_is_outer = .false.

    ! 1. finds neighbors which share a whole face with this process
    ! (faces are shared only with 1 other neighbor process)

    ! loops over all faces of element
    do iface = 1, 6

      ! chooses a point inside face
      select case (iface)
      case (1)
        ! face I == 1
        iglob = ibool(1,MID,MID,ispec)
      case (2)
        ! face I == NGLLX
        iglob = ibool(NGLLX,MID,MID,ispec)
      case (3)
        ! face J == 1
        iglob = ibool(MID,1,MID,ispec)
      case (4)
        ! face J == NGLLY
        iglob = ibool(MID,NGLLY,MID,ispec)
      case (5)
        ! face K == 1
        iglob = ibool(MID,MID,1,ispec)
      case (6)
        ! face K == NGLLZ
        iglob = ibool(MID,MID,NGLLZ,ispec)
      end select

      ! checks assembled flag on global point
      if (work_test_flag(iglob) > 0) then
        ispec_is_outer = .true.

        ! rank of neighbor process
        rank = work_test_flag(iglob) - 1

        ! checks ranks range
        if (rank < 0 .or. rank >= NPROCTOT) then
          print *,'Error face rank: ',myrank,'ispec=',ispec
          print *,'  neighbor rank = ',rank,'exceeds total nproc:',NPROCTOT
          print *,'  face ',iface
          call exit_mpi(myrank,'Error face neighbor MPI rank')
        endif

        ! checks if already stored
        icurrent = 0
        is_done = .false.
        do ii = 1,iinterface
          if (rank == my_neighbors(ii)) then
            icurrent = ii
            is_done = .true.
            exit
          endif
        enddo

        ! updates interfaces array
        if (.not. is_done) then
          iinterface = iinterface + 1
          if (iinterface > MAX_NEIGHBORS) then
            print *,'Error interfaces rank:',myrank,'iinterface = ',iinterface,MAX_NEIGHBORS
            call exit_mpi(myrank,'interface face exceeds MAX_NEIGHBORS range')
          endif
          ! adds as neighbor new interface
          my_neighbors(iinterface) = rank
          icurrent = iinterface
        endif
        if (icurrent == 0 ) &
          call exit_mpi(myrank,'could not find current interface for this neighbor, please check my_neighbors')

        ! adds interface points and removes neighbor flag from face
        ! assumes NGLLX == NGLLY == NGLLZ
        do k = 1,NGLLX
          do j = 1,NGLLX
            select case (iface)
            case (1)
              ! face I == 1
              iglob = ibool(1,j,k,ispec)
            case (2)
              ! face I == NGLLX
              iglob = ibool(NGLLX,j,k,ispec)
            case (3)
              ! face J == 1
              iglob = ibool(j,1,k,ispec)
            case (4)
              ! face J == NGLLY
              iglob = ibool(j,NGLLY,k,ispec)
            case (5)
              ! face K == 1
              iglob = ibool(j,k,1,ispec)
            case (6)
              ! face K == NGLLZ
              iglob = ibool(j,k,NGLLZ,ispec)
            end select

            ! checks that we take each global point (on edges and corners) only once
            call add_interface_point(iglob,rank,icurrent, &
                                     nibool_neighbors,MAX_NEIGHBORS, &
                                     ibool_neighbors,max_nibool, &
                                     work_test_flag,NGLOB,myrank, &
                                     .true.,add_central_cube)
            ! debug
            if (work_test_flag(iglob) < 0) then
              if (IREGION == IREGION_INNER_CORE .and. INCLUDE_CENTRAL_CUBE) then
                ! we might have missed an interface point on an edge, just re-set to missing value
                print *,'warning face flag:',myrank,'ispec=',ispec,'rank=',rank
                print *,'  flag=',work_test_flag(iglob),'iface jk=',iface,j,k,'missed iglob=',iglob
                !work_test_flag(iglob) = 0
              else
                print *,'Error face flag:',myrank,'ispec=',ispec,'rank=',rank
                print *,'  flag=',work_test_flag(iglob),'iface jk=',iface,j,k,'iglob=',iglob
                call exit_mpi(myrank,'Error face flag')
              endif
            endif

          enddo
        enddo
      endif
    enddo ! iface

    ! 2. finds neighbors which share a single edge with this process
    ! note: by now, faces have subtracted their neighbors, edges can hold only one more process info

    ! loops over all edges of element
    do iedge = 1, 12

      ! chooses a point inside edge but not corner
      select case (iedge)
      case (1)
        ! face I == 1, J == 1
        iglob = ibool(1,1,MID,ispec)
      case (2)
        ! face I == 1, J == NGLLY
        iglob = ibool(1,NGLLY,MID,ispec)
      case (3)
        ! face I == 1, K == 1
        iglob = ibool(1,MID,1,ispec)
      case (4)
        ! face I == 1, K == NGLLZ
        iglob = ibool(1,MID,NGLLZ,ispec)
      case (5)
        ! face I == NGLLX, J == 1
        iglob = ibool(NGLLX,1,MID,ispec)
      case (6)
        ! face I == NGLLX, J == NGLLY
        iglob = ibool(NGLLX,NGLLY,MID,ispec)
      case (7)
        ! face I == NGLLX, K == 1
        iglob = ibool(NGLLX,MID,1,ispec)
      case (8)
        ! face I == NGLLX, K == NGLLZ
        iglob = ibool(NGLLX,MID,NGLLZ,ispec)
      case (9)
        ! face J == 1, K == 1
        iglob = ibool(MID,1,1,ispec)
      case (10)
        ! face J == 1, K == NGLLZ
        iglob = ibool(MID,1,NGLLZ,ispec)
      case (11)
        ! face J == NGLLY, K == 1
        iglob = ibool(MID,NGLLY,1,ispec)
      case (12)
        ! face J == NGLLY, K == NGLLZ
        iglob = ibool(MID,NGLLY,NGLLZ,ispec)
      end select

      ! checks assembled flag on global point
      if (work_test_flag(iglob) > 0) then
        ispec_is_outer = .true.

        ! rank of neighbor process
        rank = work_test_flag(iglob) - 1

        ! checks ranks range
        if (rank < 0 .or. rank >= NPROCTOT) then
          print *,'Error egde rank: ',myrank
          print *,'  neighbor rank = ',rank,'exceeds total nproc:',NPROCTOT
          print *,'  edge ',iedge
          call exit_mpi(myrank,'Error edge neighbor MPI rank')
        endif

        ! checks if already stored
        icurrent = 0
        is_done = .false.
        do ii = 1,iinterface
          if (rank == my_neighbors(ii)) then
            icurrent = ii
            is_done = .true.
            exit
          endif
        enddo

        ! updates interfaces array
        if (.not. is_done) then
          iinterface = iinterface + 1
          if (iinterface > MAX_NEIGHBORS) then
            print *,'Error interfaces rank:',myrank,'iinterface = ',iinterface,MAX_NEIGHBORS
            call exit_mpi(myrank,'interface edge exceeds MAX_NEIGHBORS range')
          endif
          ! adds as neighbor new interface
          my_neighbors(iinterface) = rank
          icurrent = iinterface
        endif
        if (icurrent == 0 ) &
          call exit_mpi(myrank,'could not find current interface for this neighbor, please check my_neighbors')

        ! adds interface points and removes neighbor flag from edge
        ! assumes NGLLX == NGLLY == NGLLZ
        do k = 1,NGLLX
          select case (iedge)
          case (1)
            ! face I == 1, J == 1
            iglob = ibool(1,1,k,ispec)
          case (2)
            ! face I == 1, J == NGLLY
            iglob = ibool(1,NGLLY,k,ispec)
          case (3)
            ! face I == 1, K == 1
            iglob = ibool(1,k,1,ispec)
          case (4)
            ! face I == 1, K == NGLLZ
            iglob = ibool(1,k,NGLLZ,ispec)
          case (5)
            ! face I == NGLLX, J == 1
            iglob = ibool(NGLLX,1,k,ispec)
          case (6)
            ! face I == NGLLX, J == NGLLY
            iglob = ibool(NGLLX,NGLLY,k,ispec)
          case (7)
            ! face I == NGLLX, K == 1
            iglob = ibool(NGLLX,k,1,ispec)
          case (8)
            ! face I == NGLLX, K == NGLLZ
            iglob = ibool(NGLLX,k,NGLLZ,ispec)
          case (9)
            ! face J == 1, K == 1
            iglob = ibool(k,1,1,ispec)
          case (10)
            ! face J == 1, K == NGLLZ
            iglob = ibool(k,1,NGLLZ,ispec)
          case (11)
            ! face J == NGLLY, K == 1
            iglob = ibool(k,NGLLY,1,ispec)
          case (12)
            ! face J == NGLLY, K == NGLLZ
            iglob = ibool(k,NGLLY,NGLLZ,ispec)
          end select

          ! checks that we take each global point (on edges and corners) only once
          call add_interface_point(iglob,rank,icurrent, &
                                   nibool_neighbors,MAX_NEIGHBORS, &
                                   ibool_neighbors,max_nibool, &
                                   work_test_flag,NGLOB,myrank, &
                                   .true.,add_central_cube)

          ! debug
          if (work_test_flag(iglob) < 0) then
            if (IREGION == IREGION_INNER_CORE .and. INCLUDE_CENTRAL_CUBE) then
              ! we might have missed an interface point on an edge, just re-set to missing value
              print *,'warning edge flag:',myrank,'ispec=',ispec,'rank=',rank
              print *,'  flag=',work_test_flag(iglob),'iedge jk=',iedge,k,'missed iglob=',iglob
              !work_test_flag(iglob) = 0
            else
              print *,'Error edge flag:',myrank,'ispec=',ispec,'rank=',rank
              print *,'  flag=',work_test_flag(iglob),'iedge jk=',iedge,k,'iglob=',iglob
              call exit_mpi(myrank,'Error edge flag')
            endif
          endif

        enddo
      endif
    enddo ! iedge


    ! 3. finds neighbors which share a single corner with this process
    ! note: faces and edges have subtracted their neighbors, only one more process left possible

    ! loops over all corners of element
    do icorner = 1, 8

      ! chooses a corner point
      select case (icorner)
      case (1)
        ! face I == 1
        iglob = ibool(1,1,1,ispec)
      case (2)
        ! face I == 1
        iglob = ibool(1,NGLLY,1,ispec)
      case (3)
        ! face I == 1
        iglob = ibool(1,1,NGLLZ,ispec)
      case (4)
        ! face I == 1
        iglob = ibool(1,NGLLY,NGLLZ,ispec)
      case (5)
        ! face I == NGLLX
        iglob = ibool(NGLLX,1,1,ispec)
      case (6)
        ! face I == NGLLX
        iglob = ibool(NGLLX,NGLLY,1,ispec)
      case (7)
        ! face I == NGLLX
        iglob = ibool(NGLLX,1,NGLLZ,ispec)
      case (8)
        ! face I == NGLLX
        iglob = ibool(NGLLX,NGLLY,NGLLZ,ispec)
      end select

      ! makes sure that all elements on MPI interfaces are included
      ! uses original test_flag array, since the working copy reduces values
      ! note: there can be elements which have an edge or corner shared with
      !          other MPI partitions, but have the work_test_flag value already set to zero
      !          since the iglob point was found before.
      !          also, this check here would suffice to determine the outer flag, but we also include the
      !          check everywhere we encounter it too
      if (test_flag(iglob) > 0.5) then
        ispec_is_outer = .true.
      endif

      ! checks assembled flag on global point
      if (work_test_flag(iglob) > 0) then
        ispec_is_outer = .true.

        ! rank of neighbor process
        rank = work_test_flag(iglob) - 1

        ! checks ranks range
        if (rank < 0 .or. rank >= NPROCTOT) then
          print *,'Error corner: ',myrank
          print *,'  neighbor rank = ',rank,'exceeds total nproc:',NPROCTOT
          print *,'  corner ',icorner
          call exit_mpi(myrank,'Error corner neighbor MPI rank')
        endif

        ! checks if already stored
        icurrent = 0
        is_done = .false.
        do ii = 1,iinterface
          if (rank == my_neighbors(ii)) then
            icurrent = ii
            is_done = .true.
            exit
          endif
        enddo

        ! updates interfaces array
        if (.not. is_done) then
          iinterface = iinterface + 1
          if (iinterface > MAX_NEIGHBORS) then
            print *,'Error interfaces rank:',myrank,'iinterface = ',iinterface,MAX_NEIGHBORS
            call exit_mpi(myrank,'interface corner exceed MAX_NEIGHBORS range')
          endif
          ! adds as neighbor new interface
          my_neighbors(iinterface) = rank
          icurrent = iinterface
        endif
        if (icurrent == 0 ) &
          call exit_mpi(myrank,'could not find current interface for this neighbor, please check my_neighbors')

        ! adds this corner as interface point and removes neighbor flag from face,
        ! checks that we take each global point (on edges and corners) only once
        call add_interface_point(iglob,rank,icurrent, &
                                 nibool_neighbors,MAX_NEIGHBORS, &
                                 ibool_neighbors,max_nibool, &
                                 work_test_flag,NGLOB,myrank, &
                                 .false.,add_central_cube)

        ! debug
        if (work_test_flag(iglob) < 0 ) call exit_mpi(myrank,'Error corner flag')

      endif

    enddo ! icorner

    ! stores flags for outer elements when recognized as such
    ! (inner/outer elements separated for non-blocking MPI communications)
    if (ispec_is_outer) then
      work_ispec_is_outer(ispec) = .true.
    endif

  enddo

  ! number of outer elements (on MPI interfaces)
  npoin = count( work_ispec_is_outer )

  ! debug: user output
  if (add_central_cube) then
    print *, 'rank',myrank,'interfaces : ',iinterface
    do j = 1,iinterface
      print *, '  my_neighbors: ',my_neighbors(j),nibool_neighbors(j)
    enddo
    print *, '  test flag min/max: ',minval(work_test_flag),maxval(work_test_flag)
    print *, '  outer elements: ',npoin
    print *
  endif

  ! checks if all points were recognized
  if (minval(work_test_flag) < 0 .or. maxval(work_test_flag) > 0) then
    print *,'Error MPI interface rank: ',myrank
    print *,'  work_test_flag min/max :',minval(work_test_flag),maxval(work_test_flag)
    call exit_mpi(myrank,'Error: MPI points remain unrecognized, please check mesh interfaces')
  endif

  ! sets interfaces info
  num_interfaces = iinterface
  max_nibool_interfaces = maxval( nibool_neighbors(1:num_interfaces) )

  ! checks if unique set of neighbors
  do ii = 1,num_interfaces-1
    rank = my_neighbors(ii)
    do j = ii+1,num_interfaces
      if (rank == my_neighbors(j)) then
        print *,'test MPI: rank ',myrank,'my_neighbors:',rank,my_neighbors(j),'interfaces:',ii,j
        call exit_mpi(myrank,'Error test my_neighbors not unique')
      endif
    enddo
  enddo

  ! sorts buffers obtained to be conforming with neighbors in other slices
  do iinterface = 1,num_interfaces
    ! sorts ibool values in increasing order
    ! used to check if we have duplicates in array
    npoin = nibool_neighbors(iinterface)
    call heap_sort( npoin, ibool_neighbors(1:npoin,iinterface) )

    ! checks if unique set of iglob values
    do j = 1,npoin-1
      if (ibool_neighbors(j,iinterface) == ibool_neighbors(j+1,iinterface)) then
        if (IREGION == IREGION_INNER_CORE .and. INCLUDE_CENTRAL_CUBE) then
          ! missing points might have been counted more than once
          if (ibool_neighbors(j,iinterface) > 0) then
            print *,'warning MPI interface rank:',myrank
            print *,'  interface: ',my_neighbors(iinterface),'point: ',j,'of',npoin,'iglob=',ibool_neighbors(j,iinterface)
            ! decrease number of points
            nibool_neighbors(iinterface) = nibool_neighbors(iinterface) - 1
            if (nibool_neighbors(iinterface) <= 0) then
              print *,'Error zero MPI interface rank:',myrank,'interface=',my_neighbors(iinterface)
              call exit_mpi(myrank,'Error: zero MPI points on interface')
            endif
            ! shift values
            do k = j+1,npoin-1
              ii = ibool_neighbors(k+1,iinterface)
              ibool_neighbors(k,iinterface) = ii
            enddo
            ! re-sets values
            ibool_neighbors(npoin,iinterface) = 0
            npoin = nibool_neighbors(iinterface)
            max_nibool_interfaces = maxval( nibool_neighbors(1:num_interfaces) )
          endif
        else
          print *,'Error MPI interface rank:',myrank
          print *,'  interface: ',my_neighbors(iinterface),'point: ',j,'of',npoin,'iglob=',ibool_neighbors(j,iinterface)
          call exit_mpi(myrank,'Error: MPI points not unique on interface')
        endif
      endif
    enddo

    ! sort buffer obtained to be conforming with neighbor in other chunk
    npoin = nibool_neighbors(iinterface)
    call sort_MPI_interface(myrank,npoin,ibool_neighbors(1:npoin,iinterface), &
                                NGLOB,xstore,ystore,zstore)

  enddo

  ! re-sets flags for outer elements
  is_on_a_slice_edge(:) = work_ispec_is_outer(:)

  end subroutine get_MPI_interfaces

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sort_MPI_interface(myrank,npoin,ibool_n, &
                                    NGLOB,xstore,ystore,zstore)

  use constants, only: CUSTOM_REAL,SIZE_REAL

  implicit none

  integer,intent(in) :: myrank,npoin
  integer,dimension(npoin),intent(inout) :: ibool_n

  integer,intent(in) :: NGLOB
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: xstore,ystore,zstore

  ! local parameters
  ! arrays for sorting routine
  double precision, dimension(:), allocatable :: xstore_selected,ystore_selected,zstore_selected
  integer, dimension(:), allocatable :: ibool_selected
  integer, dimension(:), allocatable :: ninseg,iglob,locval
  logical, dimension(:), allocatable :: ifseg
  integer :: nglob_selected,i,ipoin,ier

  ! allocate arrays for buffers with maximum size
  allocate(ibool_selected(npoin), &
          xstore_selected(npoin), &
          ystore_selected(npoin), &
          zstore_selected(npoin), &
          ninseg(npoin), &
          iglob(npoin), &
          locval(npoin), &
          ifseg(npoin),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error sort MPI interface: allocating temporary sorting arrays')

  ! sets up working arrays
  do i = 1,npoin
    ipoin = ibool_n(i)

    ibool_selected(i) = ipoin

    xstore_selected(i) = dble(xstore(ipoin))
    ystore_selected(i) = dble(ystore(ipoin))
    zstore_selected(i) = dble(zstore(ipoin))
  enddo

  ! sort buffer obtained to be conforming with neighbor in other chunk
  ! sort on x, y and z, the other arrays will be swapped as well
  call sort_array_coordinates(npoin,xstore_selected,ystore_selected,zstore_selected, &
                              ibool_selected,iglob,locval,ifseg,nglob_selected,ninseg)

  ! check that no duplicate has been detected
  if (nglob_selected /= npoin) call exit_MPI(myrank,'Error sort MPI interface: duplicates detected in buffer')

  ! stores new ibool ordering
  ibool_n(1:npoin) = ibool_selected(1:npoin)

  ! frees array memory
  deallocate(ibool_selected,xstore_selected,ystore_selected,zstore_selected, &
             ninseg,iglob,locval,ifseg)


  end subroutine sort_MPI_interface

!
!-------------------------------------------------------------------------------------------------
!

  subroutine add_interface_point(iglob,rank,icurrent, &
                                 nibool_neighbors,MAX_NEIGHBORS, &
                                 ibool_neighbors,max_nibool, &
                                 work_test_flag,NGLOB,myrank, &
                                 is_face_edge,add_central_cube)


  implicit none

  integer,intent(in) :: iglob,rank,icurrent
  integer,intent(in) :: myrank

  integer,intent(in) :: MAX_NEIGHBORS,max_nibool
  integer, dimension(MAX_NEIGHBORS),intent(inout) :: nibool_neighbors
  integer, dimension(max_nibool,MAX_NEIGHBORS),intent(inout) :: ibool_neighbors

  integer,intent(in) :: NGLOB
  integer,dimension(NGLOB) :: work_test_flag

  logical,intent(in) :: is_face_edge,add_central_cube

  ! local parameters
  integer :: i
  logical :: is_done

  ! let's check and be sure for central cube
  !if (work_test_flag(iglob) <= 0 ) cycle ! continues to next point

  ! checks that we take each global point (on edges and corners) only once
  is_done = .false.
  do i = 1,nibool_neighbors(icurrent)
    if (ibool_neighbors(i,icurrent) == iglob) then
      is_done = .true.
      exit
    endif
  enddo

  ! checks if anything to do
  if (is_done) then
    ! special handling for central cube: removes rank if already added in inner core
    if (add_central_cube) then
      if (is_face_edge .and. work_test_flag(iglob) < (rank + 1)) then
        ! re-sets if we missed this rank number
        work_test_flag(iglob) = work_test_flag(iglob) + (rank + 1)
      endif
      ! re-sets flag
      work_test_flag(iglob) = work_test_flag(iglob) - ( rank + 1 )
      if (is_face_edge .and. work_test_flag(iglob) < 0) then
        ! re-sets to zero if we missed this rank number
        if (work_test_flag(iglob) == - (rank + 1 ) ) work_test_flag(iglob) = 0
      endif
    endif
    return
  endif

  ! checks if flag was set correctly
  if (work_test_flag(iglob) <= 0) then
    ! we might have missed an interface point on an edge, just re-set to missing value
    print *,'warning ',myrank,' flag: missed rank=',rank
    print *,'  flag=',work_test_flag(iglob),'missed iglob=',iglob,'interface=',icurrent
    print *
  endif
  ! we might have missed an interface point on an edge, just re-set to missing value
  if (is_face_edge) then
    if (work_test_flag(iglob) < (rank + 1)) then
      ! re-sets if we missed this rank number
      work_test_flag(iglob) = work_test_flag(iglob) + (rank + 1)
    endif
  endif

  ! adds point
  ! increases number of total points on this interface
  nibool_neighbors(icurrent) = nibool_neighbors(icurrent) + 1
  if (nibool_neighbors(icurrent) > max_nibool) &
      call exit_mpi(myrank,'interface face exceeds max_nibool range')

  ! stores interface iglob index
  ibool_neighbors( nibool_neighbors(icurrent),icurrent ) = iglob

  ! re-sets flag
  work_test_flag(iglob) = work_test_flag(iglob) - ( rank + 1 )

  ! checks
  if (is_face_edge .and. work_test_flag(iglob) < 0) then
    ! re-sets to zero if we missed this rank number
    if (work_test_flag(iglob) == - (rank + 1 ) ) work_test_flag(iglob) = 0
  endif

  end subroutine add_interface_point

