
! display mesh topology (how processor numbers are assigned in the cubed sphere)

  program display_mesh_topology

  implicit none

  integer, parameter :: NCHUNKS = 6

! define block type based upon chunk number (between 1 and 6)
! do not change this numbering, chunk AB must be number 1 for central cube
  integer, parameter :: CHUNK_AB = 1
  integer, parameter :: CHUNK_AC = 2
  integer, parameter :: CHUNK_BC = 3
  integer, parameter :: CHUNK_AC_ANTIPODE = 4
  integer, parameter :: CHUNK_BC_ANTIPODE = 5
  integer, parameter :: CHUNK_AB_ANTIPODE = 6

  integer NPROC_XI,NPROC_ETA,NPROC
  integer iproc_xi,iproc_eta,ichunk

! addressing for all the slices
  integer, dimension(:,:,:), allocatable :: addressing

  print *
  print *,'enter the number of processors in each direction of each of the 6 blocks (NPROC_XI):'
  print *
  read(*,*) NPROC_XI

  NPROC_ETA = NPROC_XI

! total number of processors in each of the six chunks
  NPROC = NPROC_XI * NPROC_ETA

  print *
  print *,'the total number of processors used is ',NCHUNKS*NPROC
  print *

  allocate(addressing(NCHUNKS,0:NPROC_XI-1,0:NPROC_ETA-1))

! loop on all the chunks to create global slice addressing for solver
  do ichunk = 1,NCHUNKS
    do iproc_eta=0,NPROC_ETA-1
      do iproc_xi=0,NPROC_XI-1
        addressing(ichunk,iproc_xi,iproc_eta) = (ichunk-1)*NPROC + iproc_eta * NPROC_XI + iproc_xi
      enddo
    enddo
  enddo

! output a topology map of slices
    write(*,*) 'Bloc 1'
    do iproc_xi = NPROC_XI-1, 0, -1
      write(*,'(2x)',advance='no')
      do iproc_eta = NPROC_ETA -1, 0, -1
        write(*,'(i5)',advance='no') addressing(CHUNK_AB,iproc_xi,iproc_eta)
      enddo
      write(*,'(1x)',advance='yes')
    enddo
    write(*, *) ' '

    write(*,*) 'Bloc 2'
    do iproc_xi = NPROC_XI-1, 0, -1
      write(*,'(2x)',advance='no')
      do iproc_eta = NPROC_ETA -1, 0, -1
        write(*,'(i5)',advance='no') addressing(CHUNK_BC,iproc_xi,iproc_eta)
      enddo
      write(*,'(1x)',advance='yes')
    enddo
    write(*, *) ' '

    write(*,*) 'Bloc 3'
    do iproc_xi = NPROC_XI-1, 0, -1
      write(*,'(2x)',advance='no')
      do iproc_eta = NPROC_ETA -1, 0, -1
        write(*,'(i5)',advance='no') addressing(CHUNK_AC,iproc_xi,iproc_eta)
      enddo
      write(*,'(1x)',advance='yes')
    enddo
    write(*, *) ' '

    write(*,*) 'Bloc 4'
    do iproc_xi = NPROC_XI-1, 0, -1
      write(*,'(2x)',advance='no')
      do iproc_eta = NPROC_ETA -1, 0, -1
        write(*,'(i5)',advance='no') addressing(CHUNK_BC_ANTIPODE,iproc_xi,iproc_eta)
      enddo
      write(*,'(1x)',advance='yes')
    enddo
    write(*, *) ' '

    write(*,*) 'Bloc 5'
    do iproc_xi = NPROC_XI-1, 0, -1
      write(*,'(2x)',advance='no')
      do iproc_eta = NPROC_ETA -1, 0, -1
        write(*,'(i5)',advance='no') addressing(CHUNK_AB_ANTIPODE,iproc_xi,iproc_eta)
      enddo
      write(*,'(1x)',advance='yes')
    enddo
    write(*, *) ' '

    write(*,*) 'Bloc 6'
    do iproc_xi = NPROC_XI-1, 0, -1
      write(*,'(2x)',advance='no')
      do iproc_eta = NPROC_ETA -1, 0, -1
        write(*,'(i5)',advance='no') addressing(CHUNK_AC_ANTIPODE,iproc_xi,iproc_eta)
      enddo
      write(*,'(1x)',advance='yes')
    enddo
    write(*, *) ' '

  end program display_mesh_topology

