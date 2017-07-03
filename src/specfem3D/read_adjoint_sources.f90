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

  subroutine read_adjoint_sources()

! reads in adjoint source files
!
! note: depending on the number of adjoint source files and cluster filesystem,
!       this can take a very long time and slow down adjoint simulations considerably.
!
!       when flag IO_ASYNC_COPY is set, a new parallel thread is created to read in adjoint source files,
!       while the main processes can continue computing. we thus overlap file i/o with computations.
!       in case the heavy computations take long enough, this overlap should perfectly hide the i/o latency.
!
!       the costs are then for the additional memory, i.e. buffer array (buffer_source_adjoint), and the
!       copying operation of data from the buffer array to the actual adj_source_adjoint array.

  use specfem_par

  implicit none

  ! local parameters
  integer :: it_sub_adj
  ! debug timing
  !double precision, external :: wtime
  !double precision :: tstart


  ! we already checked that nadj_rec_local > 0 before calling this routine

  ! debug timing
  !tstart = wtime()

  ! determines chunk_number
  it_sub_adj = ceiling( dble(it)/dble(NTSTEP_BETWEEN_READ_ADJSRC) )

  ! debug
  !print *,'read adjoint sources: it_sub_adj = ',it_sub_adj

  ! asynchronously reads in adjoint source files
  if (IO_ASYNC_COPY .and. NSTEP_SUB_ADJ > 1) then
    ! handles file input/output thread
    if (it == it_begin) then
      ! creates new io thread for reading in sources
      call read_adj_io_thread(it_sub_adj)

      ! first chunk of adjoint sources must ready at beginning, so we wait.
      ! waits for previous read to finish and
      ! copy over buffered data into tmp_source_adjoint
      call sync_adj_io_thread(source_adjoint)

    else
      ! waits for previous read to finish and
      ! copy over buffered data into tmp_source_adjoint
      call sync_adj_io_thread(source_adjoint)
    endif

    ! checks if next chunk necessary
    if (it_sub_adj < NSTEP_SUB_ADJ) then
      ! starts thread to read in next junk
      it_sub_adj = it_sub_adj + 1

      ! creates new io thread for reading in sources
      call read_adj_io_thread(it_sub_adj)
    endif

  else
    ! synchronous read routine

    ! reads in local adjoint sources
    call read_adjoint_sources_local(source_adjoint,nadj_rec_local,it_sub_adj)

  endif

  ! debug timing
  !print *,'read adjoint sources: elapsed time = ',wtime() - tstart
  !print *

  end subroutine read_adjoint_sources

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_adjoint_sources_local(source_adjoint,nadj_rec_local,it_sub_adj)

! reads in local adjoint source files

  use specfem_par, only: myrank,NPROCTOT_VAL, &
    nrec,islice_selected_rec,station_name,network_name, &
    nu, &
    iadjsrc_len,iadjsrc,NSTEP_SUB_ADJ, &
    DT,CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ,NTSTEP_BETWEEN_READ_ADJSRC,MAX_STRING_LEN,READ_ADJSRC_ASDF

  implicit none

  integer,intent(in) :: nadj_rec_local
  integer,intent(in) :: it_sub_adj
  real(kind=CUSTOM_REAL), dimension(NDIM,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC) :: source_adjoint

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: tmp_source_adjoint
  integer :: irec,irec_local,itime,ier
  character(len=MAX_STRING_LEN) :: adj_source_file

  ! debug
  !print *,'reading adjoint sources local:',myrank,' - chunk ',it_sub_adj,'out of ',NSTEP_SUB_ADJ, &
  !       ' for local adjoint sources = ',nadj_rec_local

  ! checks chunk number
  if (it_sub_adj < 1 .or. it_sub_adj > NSTEP_SUB_ADJ) then
    print *,'Error reading adjoint sources: chunk number ',it_sub_adj,'is invalid'
    call exit_MPI(myrank,'Error reading adjoint sources with invalid chunk number')
  endif

  ! allocates temporary source array
  allocate(tmp_source_adjoint(NDIM,NTSTEP_BETWEEN_READ_ADJSRC),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating array tmp_source_adjoint')

  ! initializes
  tmp_source_adjoint(:,:) = 0._CUSTOM_REAL

  ! counter
  irec_local = 0

  ! loops over all adjoint sources
  do irec = 1, nrec
    ! checks that the source slice number is okay
    if (islice_selected_rec(irec) < 0 .or. islice_selected_rec(irec) > NPROCTOT_VAL-1) then
      print *,'Error rank ',myrank,': adjoint source slice index ',islice_selected_rec(irec), &
             ' is out of bounds ',NPROCTOT_VAL-1
      call exit_MPI(myrank,'Error adjoint source has wrong source slice number in adjoint simulation')
    endif

    ! compute source arrays for adjoint sources within this rank's slice
    if (myrank == islice_selected_rec(irec)) then
      ! increases counter
      irec_local = irec_local + 1

      ! adjoint source file name **net**.**sta**
      if (READ_ADJSRC_ASDF) then
        adj_source_file = trim(network_name(irec))//'_'//trim(station_name(irec))
      else
        adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))
      endif

      ! reads in *.adj files
      call compute_arrays_source_adjoint(myrank,adj_source_file, &
                                         nu(:,:,irec),tmp_source_adjoint, &
                                         iadjsrc_len(it_sub_adj),iadjsrc,it_sub_adj, &
                                         NSTEP_SUB_ADJ,NTSTEP_BETWEEN_READ_ADJSRC,DT)

      ! stores source array
      ! note: the source_adjoint has a time stepping from 1 to NTSTEP_BETWEEN_READ_ADJSRC
      !          this gets overwritten every time a new block/chunk is read in
      do itime = 1,NTSTEP_BETWEEN_READ_ADJSRC
        source_adjoint(:,irec_local,itime) = tmp_source_adjoint(:,itime)
      enddo

    endif
  enddo

  ! checks that number of read sources is valid
  if (irec_local /= nadj_rec_local) then
    call exit_MPI(myrank,'irec_local /= nadj_rec_local in adjoint simulation')
  endif

  ! frees temporary array
  deallocate(tmp_source_adjoint)

  end subroutine read_adjoint_sources_local

!
!-------------------------------------------------------------------------------------------------
!

  subroutine check_adjoint_sources(irec,nadj_files_found)

  use specfem_par

  implicit none

  integer,intent(in) :: irec
  ! counter
  integer,intent(inout) :: nadj_files_found

  ! local parameters
  double precision :: junk
  integer :: icomp,itime
  integer :: ier
  character(len=MAX_STRING_LEN) :: filename,adj_source_file,path_to_add
  character(len=3),dimension(NDIM) :: comp
  character(len=2) :: bic

  ! checks **net**.**sta**.**MX**.adj files for correct number of time steps
  adj_source_file = trim(network_name(irec))//'.'//trim(station_name(irec))
  ! adjoint source file name **net**.**sta**

  ! bandwidth code
  call band_instrument_code(DT,bic)
  comp(1) = bic(1:2)//'N'
  comp(2) = bic(1:2)//'E'
  comp(3) = bic(1:2)//'Z'

  ! loops over file components E/N/Z
  do icomp = 1,NDIM

    ! opens adjoint source file for this component
    filename = 'SEM/'//trim(adj_source_file) // '.'// comp(icomp) // '.adj'

    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
      write(path_to_add,"('run',i4.4,'/')") mygroup + 1
      filename = path_to_add(1:len_trim(path_to_add))//filename(1:len_trim(filename))
    endif

    open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ier)

    ! checks if file opens/exists
    if (ier /= 0) then
      ! adjoint source file not found
      ! stops simulation
      call exit_MPI(myrank, &
          'file '//trim(filename)//' not found, please check with your STATIONS_ADJOINT file')
    endif

    ! checks length of file
    itime = 0
    do while(ier == 0)
      read(IIN,*,iostat=ier) junk,junk
      if (ier == 0 ) itime = itime + 1
    enddo

    ! checks length
    if (itime /= NSTEP) then
      print *,'adjoint source error: ',trim(filename),' has length',itime,' but should be',NSTEP
      call exit_MPI(myrank, &
        'file '//trim(filename)//' length is wrong, please check your adjoint sources and your simulation duration')
    endif

    ! updates counter for found files
    nadj_files_found = nadj_files_found + 1

    ! closes file
    close(IIN)

  enddo

  end subroutine check_adjoint_sources

