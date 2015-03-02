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


  subroutine test_MPI_neighbours(iregion_code, &
                                 num_interfaces,max_nibool_interfaces, &
                                 my_neighbours,nibool_interfaces, &
                                 ibool_interfaces)

  use constants
  use meshfem3D_par,only: NPROCTOT,myrank
  use MPI_crust_mantle_par,only: NGLOB_CRUST_MANTLE
  use MPI_outer_core_par,only: NGLOB_OUTER_CORE
  use MPI_inner_core_par,only: NGLOB_INNER_CORE

  implicit none

  integer,intent(in) :: iregion_code
  integer,intent(in) :: num_interfaces,max_nibool_interfaces
  integer,dimension(num_interfaces),intent(in) :: my_neighbours,nibool_interfaces
  integer,dimension(max_nibool_interfaces,num_interfaces),intent(in):: ibool_interfaces

  ! local parameters
  integer,dimension(:),allocatable :: dummy_i
  integer,dimension(:,:),allocatable :: test_interfaces
  integer,dimension(:,:),allocatable :: test_interfaces_nibool
  integer :: ineighbour,iproc,inum,i,j,ier,ipoints,max_num,iglob
  logical :: is_okay
  logical,dimension(:),allocatable :: mask

  ! debug output
  !do iproc = 0,NPROCTOT-1
  !  if (myrank == iproc) then
  !    print*, 'MPI rank',myrank,'interfaces : ',num_interfaces,'region',iregion_code
  !    do j = 1,num_interfaces
  !      print*, '  my_neighbours: ',my_neighbours(j),nibool_interfaces(j)
  !    enddo
  !    print*
  !  endif
  !  call synchronize_all()
  !enddo

  ! checks maximum number of interface points
  if (max_nibool_interfaces == 0 .and. NPROCTOT > 1) then
    print*,'test MPI: rank ',myrank,'max_nibool_interfaces is zero'
    call exit_mpi(myrank,'Error test max_nibool_interfaces zero')
  endif

  ! allocates global mask
  select case (iregion_code)
  case (IREGION_CRUST_MANTLE)
    allocate(mask(NGLOB_CRUST_MANTLE),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating mask for testing mpi neighbors')
  case (IREGION_OUTER_CORE)
    allocate(mask(NGLOB_OUTER_CORE),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating mask for testing mpi neighbors')
  case (IREGION_INNER_CORE)
    allocate(mask(NGLOB_INNER_CORE),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating mask for testing mpi neighbors')
  case default
    call exit_mpi(myrank,'Error test MPI: iregion_code not recognized')
  end select

  ! test ibool entries
  ! (must be non-zero and unique)
  do i = 1,num_interfaces
    ! number of interface points
    if (nibool_interfaces(i) > max_nibool_interfaces) then
      print*,'Error test MPI: rank',myrank,'nibool values:',nibool_interfaces(i),max_nibool_interfaces
      call exit_mpi(myrank,'Error test MPI: nibool exceeds max_nibool_interfaces')
    endif

    mask(:) = .false.

    ! ibool entries
    do j = 1,nibool_interfaces(i)
      iglob = ibool_interfaces(j,i)

      ! checks zero entry
      if (iglob <= 0) then
        print*,'Error test MPI: rank ',myrank,'ibool value:',iglob,'interface:',i,'point:',j
        call exit_mpi(myrank,'Error test MPI: ibool values invalid')
      endif

      ! checks duplicate
      if (j < nibool_interfaces(i)) then
        if (iglob == ibool_interfaces(j+1,i)) then
          print*,'Error test MPI: rank',myrank,'ibool duplicate:',iglob,'interface:',i,'point:',j
          call exit_mpi(myrank,'Error test MPI: ibool duplicates')
        endif
      endif

      ! checks if unique global value
      if (.not. mask(iglob)) then
        mask(iglob) = .true.
      else
        print*,'Error test MPI: rank',myrank,'ibool masked:',iglob,'interface:',i,'point:',j
        call exit_mpi(myrank,'Error test MPI: ibool masked already')
      endif
    enddo
  enddo
  deallocate(mask)

  ! checks neighbors
  ! gets maximum interfaces from all processes
  call max_all_i(num_interfaces,max_num)

  ! master gathers info
  if (myrank == 0) then
    ! user output
    write(IMAIN,*) '  maximum interfaces:',max_num

    ! array for gathering info
    allocate(test_interfaces(max_num,0:NPROCTOT),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating test_interfaces')
    test_interfaces(:,:) = -1

    allocate(test_interfaces_nibool(max_num,0:NPROCTOT),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating test_interfaces_nibool')
    test_interfaces_nibool(:,:) = 0

    ! used to store number of interfaces per proc
    allocate(dummy_i(0:NPROCTOT),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating dummy_i for test interfaces')
    dummy_i(:) = 0

    ! sets info for master process
    test_interfaces(1:num_interfaces,0) = my_neighbours(1:num_interfaces)
    test_interfaces_nibool(1:num_interfaces,0) = nibool_interfaces(1:num_interfaces)
    dummy_i(0) = num_interfaces

    ! collects from other processes
    do iproc = 1,NPROCTOT-1
      ! gets number of interfaces
      call recv_singlei(inum,iproc,itag)
      dummy_i(iproc) = inum
      if (inum > 0) then
        call recv_i(test_interfaces(1:inum,iproc),inum,iproc,itag)
        call recv_i(test_interfaces_nibool(1:inum,iproc),inum,iproc,itag)
      endif
    enddo
  else
    ! sends info to master process
    call send_singlei(num_interfaces,0,itag)
    if (num_interfaces > 0) then
      call send_i(my_neighbours(1:num_interfaces),num_interfaces,0,itag)
      call send_i(nibool_interfaces(1:num_interfaces),num_interfaces,0,itag)
    endif
  endif
  call synchronize_all()

  ! checks if addressing is okay
  if (myrank == 0) then
    ! for each process
    do iproc = 0,NPROCTOT-1
      ! loops over all neighbors
      do i = 1,dummy_i(iproc)
        ! gets neighbour rank and number of points on interface with it
        ineighbour = test_interfaces(i,iproc)
        ipoints = test_interfaces_nibool(i,iproc)

        ! checks values
        if (ineighbour < 0 .or. ineighbour > NPROCTOT-1) then
          print*,'Error neighbour:',iproc,ineighbour
          call exit_mpi(myrank,'Error ineighbour')
        endif
        if (ipoints <= 0) then
          print*,'Error neighbour points:',iproc,ipoints
          call exit_mpi(myrank,'Error ineighbour points')
        endif

        ! looks up corresponding entry in neighbour array
        is_okay = .false.
        do j = 1,dummy_i(ineighbour)
          if (test_interfaces(j,ineighbour) == iproc) then
            ! checks if same number of interface points with this neighbour
            if (test_interfaces_nibool(j,ineighbour) == ipoints) then
              is_okay = .true.
            else
              print*,'Error ',iproc,'neighbour ',ineighbour,' points =',ipoints
              print*,'  ineighbour has points = ',test_interfaces_nibool(j,ineighbour)
              print*
              call exit_mpi(myrank,'Error ineighbour points differ')
            endif
            exit
          endif
        enddo
        if (.not. is_okay) then
          print*,'Error ',iproc,' neighbour not found: ',ineighbour
          print*,'iproc ',iproc,' interfaces:'
          print*,test_interfaces(1:dummy_i(iproc),iproc)
          print*,'ineighbour ',ineighbour,' interfaces:'
          print*,test_interfaces(1:dummy_i(ineighbour),ineighbour)
          print*
          call exit_mpi(myrank,'Error ineighbour not found')
        endif
      enddo
    enddo

    ! user output
    write(IMAIN,*) '  MPI addressing maximum interfaces:',maxval(dummy_i)
    write(IMAIN,*) '  MPI addressing : all interfaces okay'
    write(IMAIN,*)
    call flush_IMAIN()

    deallocate(dummy_i)
    deallocate(test_interfaces)
    deallocate(test_interfaces_nibool)
  endif
  call synchronize_all()

  end subroutine test_MPI_neighbours

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_MPI_cm()

  use meshfem3D_par,only: NPROCTOT,myrank
  use create_MPI_interfaces_par
  use MPI_crust_mantle_par

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: test_flag_vector
  integer :: i,j,iglob,ier
  integer :: inum,icount
  integer :: num_unique,max_valence
  integer,dimension(:),allocatable :: valence

  ! crust mantle
  allocate(test_flag_vector(NDIM,NGLOB_CRUST_MANTLE),stat=ier)
  if (ier /= 0 ) stop 'Error allocating array test_flag crust/mantle'
  allocate(valence(NGLOB_CRUST_MANTLE),stat=ier)
  if (ier /= 0 ) stop 'Error allocating array valence'

  ! points defined by interfaces
  valence(:) = 0
  test_flag_vector(:,:) = 0.0
  do i = 1,num_interfaces_crust_mantle
    do j = 1,nibool_interfaces_crust_mantle(i)
      iglob = ibool_interfaces_crust_mantle(j,i)
      ! sets flag on
      test_flag_vector(1,iglob) = 1.0_CUSTOM_REAL
      ! counts valence (occurrences)
      valence(iglob) = valence(iglob) + 1
    enddo
  enddo
  ! total number of  interface points
  i = sum(nibool_interfaces_crust_mantle)
  call sum_all_i(i,inum)

  ! total number of unique points (some could be shared between different processes)
  i = nint( sum(test_flag_vector) )
  num_unique= i
  call sum_all_i(i,icount)

  ! maximum valence
  i = maxval( valence(:) )
  call max_all_i(i,max_valence)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  total MPI interface points : ',inum
    write(IMAIN,*) '  unique MPI interface points: ',icount
    write(IMAIN,*) '  maximum valence            : ',max_valence
    call flush_IMAIN()
  endif

  ! initializes for assembly
  test_flag_vector(:,:) = 1.0_CUSTOM_REAL

  ! adds contributions from different partitions to flag arrays
  call assemble_MPI_vector(NPROCTOT,NGLOB_CRUST_MANTLE, &
                      test_flag_vector, &
                      num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                      nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle,&
                      my_neighbours_crust_mantle)

  ! removes initial flag
  test_flag_vector(:,:) = test_flag_vector(:,:) - 1.0_CUSTOM_REAL

  ! checks number of interface points
  i = 0
  do iglob = 1,NGLOB_CRUST_MANTLE
    ! only counts flags with MPI contributions
    if (test_flag_vector(1,iglob) > 0.0 ) i = i + 1

    ! checks valence
    if (valence(iglob) /= nint(test_flag_vector(1,iglob)) .or. &
       valence(iglob) /= nint(test_flag_vector(2,iglob)) .or. &
       valence(iglob) /= nint(test_flag_vector(3,iglob))) then
      print*,'Error test MPI: rank',myrank,'valence:',valence(iglob),'flag:',test_flag_vector(:,:)
      call exit_mpi(myrank,'Error test MPI crust mantle valence')
    endif
  enddo

  ! checks within slice
  if (i /= num_unique) then
    print*,'Error test crust mantle : rank',myrank,'unique MPI points:',i,num_unique
    call exit_mpi(myrank,'Error MPI assembly crust mantle')
  endif

  ! total number of assembly points
  call sum_all_i(i,inum)

  ! points defined by interfaces
  if (myrank == 0) then
    ! checks
    if (inum /= icount) then
      print*,'Error crust mantle : total MPI points:',myrank,'total: ',inum,icount
      call exit_mpi(myrank,'Error MPI assembly crust mantle')
    endif

    ! user output
    write(IMAIN,*) '  total unique MPI interface points:',inum
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  deallocate(test_flag_vector)
  deallocate(valence)

  call synchronize_all()

  end subroutine test_MPI_cm

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_MPI_oc()

  use meshfem3D_par,only: NPROCTOT,myrank
  use create_MPI_interfaces_par
  use MPI_outer_core_par

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: test_flag
  integer :: i,j,iglob,ier
  integer :: inum,icount
  integer :: max_valence,num_unique
  integer,dimension(:),allocatable :: valence

  ! outer core
  allocate(test_flag(NGLOB_OUTER_CORE),stat=ier)
  if (ier /= 0 ) stop 'Error allocating array test_flag outer core'
  allocate(valence(NGLOB_OUTER_CORE),stat=ier)
  if (ier /= 0 ) stop 'Error allocating array valence'

  ! points defined by interfaces
  valence(:) = 0
  test_flag = 0.0
  do i = 1,num_interfaces_outer_core
    do j = 1,nibool_interfaces_outer_core(i)
      iglob = ibool_interfaces_outer_core(j,i)
      test_flag(iglob) = 1.0_CUSTOM_REAL
      ! counts valence (occurrences)
      valence(iglob) = valence(iglob) + 1
    enddo
  enddo
  i = sum(nibool_interfaces_outer_core)
  call sum_all_i(i,inum)

  i = nint( sum(test_flag) )
  num_unique = i
  call sum_all_i(i,icount)

  ! maximum valence
  i = maxval( valence(:) )
  call max_all_i(i,max_valence)

  if (myrank == 0) then
    write(IMAIN,*) '  total MPI interface points : ',inum
    write(IMAIN,*) '  unique MPI interface points: ',icount
    write(IMAIN,*) '  maximum valence            : ',max_valence
  endif

  ! initialized for assembly
  test_flag(:) = 1.0_CUSTOM_REAL

  ! adds contributions from different partitions to flag arrays
  call assemble_MPI_scalar(NPROCTOT,NGLOB_OUTER_CORE, &
                                test_flag, &
                                num_interfaces_outer_core,max_nibool_interfaces_oc, &
                                nibool_interfaces_outer_core,ibool_interfaces_outer_core,&
                                my_neighbours_outer_core)


  ! removes initial flag
  test_flag(:) = test_flag(:) - 1.0_CUSTOM_REAL

  ! checks number of interface points
  i = 0
  do iglob = 1,NGLOB_OUTER_CORE
    ! only counts flags with MPI contributions
    if (test_flag(iglob) > 0.0 ) i = i + 1

    ! checks valence
    if (valence(iglob) /= nint(test_flag(iglob))) then
      print*,'Error test MPI: rank',myrank,'valence:',valence(iglob),'flag:',test_flag(iglob)
      call exit_mpi(myrank,'Error test outer core valence')
    endif
  enddo

  ! checks within slice
  if (i /= num_unique) then
    print*,'Error test outer core : rank',myrank,'unique MPI points:',i,num_unique
    call exit_mpi(myrank,'Error MPI assembly outer core')
  endif
  call sum_all_i(i,inum)

  ! output
  if (myrank == 0) then
    ! checks
    if (inum /= icount) then
      print*,'Error outer core : total MPI points:',myrank,'total: ',inum,icount
      call exit_mpi(myrank,'Error MPI assembly outer_core')
    endif

    ! user output
    write(IMAIN,*) '  total assembled MPI interface points:',inum
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  deallocate(test_flag)
  deallocate(valence)

  call synchronize_all()

  end subroutine test_MPI_oc


!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_MPI_ic()

  use meshfem3D_par,only: NPROCTOT,myrank
  use create_MPI_interfaces_par
  use MPI_inner_core_par

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: test_flag_vector
  integer :: i,j,iglob,ier
  integer :: inum,icount
  integer :: num_unique,max_valence
  integer,dimension(:),allocatable :: valence

  ! inner core
  allocate(test_flag_vector(NDIM,NGLOB_INNER_CORE),stat=ier)
  if (ier /= 0 ) stop 'Error allocating array test_flag inner core'
  allocate(valence(NGLOB_INNER_CORE),stat=ier)
  if (ier /= 0 ) stop 'Error allocating array valence'

  ! points defined by interfaces
  valence(:) = 0
  test_flag_vector(:,:) = 0.0
  do i = 1,num_interfaces_inner_core
    do j = 1,nibool_interfaces_inner_core(i)
      iglob = ibool_interfaces_inner_core(j,i)
      ! sets flag on
      test_flag_vector(1,iglob) = 1.0_CUSTOM_REAL
      ! counts valence (occurrences)
      valence(iglob) = valence(iglob) + 1
    enddo
  enddo
  i = sum(nibool_interfaces_inner_core)
  call sum_all_i(i,inum)

  i = nint( sum(test_flag_vector) )
  num_unique= i
  call sum_all_i(i,icount)

  ! maximum valence
  i = maxval( valence(:) )
  call max_all_i(i,max_valence)

  if (myrank == 0) then
    write(IMAIN,*) '  total MPI interface points : ',inum
    write(IMAIN,*) '  unique MPI interface points: ',icount
    write(IMAIN,*) '  maximum valence            : ',max_valence
  endif

  ! initializes for assembly
  test_flag_vector = 1.0_CUSTOM_REAL

  ! adds contributions from different partitions to flag arrays
  call assemble_MPI_vector(NPROCTOT,NGLOB_INNER_CORE, &
                      test_flag_vector, &
                      num_interfaces_inner_core,max_nibool_interfaces_ic, &
                      nibool_interfaces_inner_core,ibool_interfaces_inner_core,&
                      my_neighbours_inner_core)

  ! removes initial flag
  test_flag_vector(:,:) = test_flag_vector(:,:) - 1.0_CUSTOM_REAL

  ! checks number of interface points
  i = 0
  do iglob = 1,NGLOB_INNER_CORE
    ! only counts flags with MPI contributions
    if (test_flag_vector(1,iglob) > 0.0 ) i = i + 1

    ! checks valence
    if (valence(iglob) /= nint(test_flag_vector(1,iglob)) .or. &
       valence(iglob) /= nint(test_flag_vector(2,iglob)) .or. &
       valence(iglob) /= nint(test_flag_vector(3,iglob))) then
      print*,'Error test MPI: rank',myrank,'valence:',valence(iglob),'flag:',test_flag_vector(:,:)
      call exit_mpi(myrank,'Error test MPI inner core valence')
    endif

  enddo

  ! checks within slice
  if (i /= num_unique) then
    print*,'Error test inner core : rank',myrank,'unique MPI points:',i,num_unique
    call exit_mpi(myrank,'Error MPI assembly inner core')
  endif
  call sum_all_i(i,inum)

  if (myrank == 0) then
    ! checks
    if (inum /= icount) then
      print*,'Error inner core : total MPI points:',myrank,'total: ',inum,icount
      call exit_mpi(myrank,'Error MPI assembly inner core')
    endif

    ! user output
    write(IMAIN,*) '  total assembled MPI interface points:',inum
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  deallocate(test_flag_vector)
  deallocate(valence)

  call synchronize_all()

  end subroutine test_MPI_ic
