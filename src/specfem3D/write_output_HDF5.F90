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


  subroutine split_string(input_string, delimiter, part1, part2)

  implicit none
  character(len=*), intent(in) :: input_string
  character(len=*), intent(in) :: delimiter
  character(len=*), intent(out) :: part1
  character(len=*), intent(out) :: part2
  integer :: delim_pos

  ! Find the position of the delimiter
  delim_pos = index(input_string, delimiter)

  ! Split the string at the delimiter
  if (delim_pos > 0) then
      part1 = input_string(1:delim_pos-1)
      part2 = input_string(delim_pos+1:)
  else
      part1 = input_string
      part2 = ''
  endif

  end subroutine split_string

!
!-------------------------------------------------------------------------------------------------
!

#ifdef USE_HDF5

  subroutine write_hdf5_seismogram_init(nlength_total_seismogram)

  use specfem_par, only: CUSTOM_REAL, MAX_LENGTH_NETWORK_NAME, MAX_LENGTH_STATION_NAME, &
                         MAX_STRING_LEN, t0, DT, IMAIN, IOUT, myrank, &
                         NSTEP, nrec, hdf5_seismo_fname, &
                         OUTPUT_FILES, SIMULATION_TYPE, WRITE_SEISMOGRAMS_BY_MAIN

  use shared_parameters, only: IO_compute_task, NTSTEP_BETWEEN_OUTPUT_SAMPLE

  use manager_hdf5

  implicit none

  ! local parameters
  integer :: i,irec,ier
  integer, intent(in) :: nlength_total_seismogram

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: time_array
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rec_dists ! store the distances
  character(len=MAX_LENGTH_STATION_NAME), dimension(:), allocatable :: stations
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(:), allocatable :: networks
  integer, parameter :: n_header_lines = 3
  character(len=MAX_STRING_LEN) :: tmpstr, dummy, line

  ! only main process writes out
  if (myrank /= 0) return

  ! create the file and all the datasets for the seismograms
  ! set the name of the file
  hdf5_seismo_fname = trim(OUTPUT_FILES)//'/seismograms.h5'

  ! user output
  if (myrank == 0 .and. IO_compute_task) then
    write(IMAIN,*) 'Creating seismograms in HDF5 file format'
    if (WRITE_SEISMOGRAMS_BY_MAIN) then
      write(IMAIN,*) '  writing waveforms by main...'
    else
      write(IMAIN,*) '  writing waveforms in parallel...'
    endif
    write(IMAIN,*) '  seismogram file: ',trim(hdf5_seismo_fname)
    call flush_IMAIN()
  endif

  ! initialize the HDF5 manager
  call h5_initialize()

  ! create the file
  call h5_create_file(hdf5_seismo_fname)

  ! create time dataset it = 1 ~ NTSTEP
  allocate(time_array(nlength_total_seismogram),stat=ier)
  if (ier /= 0) stop 'Error: write_hdf5_seismogram_init: time_array allocation failed'
  time_array(:) = 0.0_CUSTOM_REAL

  do i = 1,nlength_total_seismogram
    if (SIMULATION_TYPE == 1) then
      ! forward simulation
      time_array(i) = real( dble((i-1)*NTSTEP_BETWEEN_OUTPUT_SAMPLE) * DT - t0, kind=CUSTOM_REAL)
    else
      ! adjoint simulation: backward/reconstructed wavefields
      time_array(i) = real( dble((NSTEP-i)*NTSTEP_BETWEEN_OUTPUT_SAMPLE) * DT - t0, kind=CUSTOM_REAL)
    endif
  enddo

  ! write the time dataset
  call h5_write_dataset_no_group('time',time_array)
  ! free the memory
  deallocate(time_array)

  ! read output_list_stations.txt generated at locate_receivers.f90:613 here to write in the h5 file.
  allocate(stations(nrec), &
           networks(nrec), &
           rec_dists(nrec),stat=ier)
  if (ier /= 0) stop 'Error: write_hdf5_seismogram_init: stations allocation failed'
  rec_dists(:) = 0.0_CUSTOM_REAL

  open(unit=IOUT,file=trim(OUTPUT_FILES)//'/output_list_stations.txt', &
       status='unknown',action='read',iostat=ier)
  if (ier /= 0) then
    stop 'Error: write_hdf5_seismogram_init: error opening output_list_stations.txt'
  endif

  ! skip the header (3 lines)
  do i = 1,n_header_lines
    read(IOUT,*)
  enddo

  ! read the station information
  ! each line includes: network.station dummy dummy dist dummy
  do irec = 1,nrec
    ! read network.station and separate them by '.'
    !read(IOUT) tmpstr, dummy, dummy, rec_dists(irec), dummy
    !read(IOUT, '(A, 1X, A, 1X, A, 3X, F9.7, 6X, A)') tmpstr, dummy, dummy, rec_dists(irec), dummy
    read(IOUT, '(A)', iostat=ier) line
    if (ier /= 0) stop 'Error: write_hdf5_seismogram_init: error reading output_list_stations.txt'
    read(line, *) tmpstr, dummy, dummy, rec_dists(irec), dummy
    call split_string(trim(tmpstr), '.', networks(irec), stations(irec))
  enddo

  close(IOUT)

  ! write the station information
  call h5_write_dataset_no_group('stations',stations)
  call h5_write_dataset_no_group('networks',networks)
  call h5_write_dataset_no_group('dists',rec_dists)

  ! free the memory
  deallocate(stations,networks,rec_dists)

  ! close the file
  call h5_close_file()

  end subroutine write_hdf5_seismogram_init

#endif

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_output_hdf5(seismogram_tmp_in, irec_local, irec, chn, iorientation)

  use specfem_par, only: &
    nlength_seismogram, &
    CUSTOM_REAL

#ifdef USE_HDF5
  use specfem_par, only: &
    myrank, seismo_current, seismo_offset, nrec, &
    NTSTEP_BETWEEN_OUTPUT_SAMPLE, &
    WRITE_SEISMOGRAMS_BY_MAIN, hdf5_seismo_fname
  use shared_parameters, only: &
    NSTEP, OUTPUT_SEISMOS_HDF5
  use manager_hdf5
#endif

  implicit none

  ! input/output variables
  character(len=4),intent(in) :: chn
  integer,intent(in) :: irec_local, irec, iorientation
  real(kind=CUSTOM_REAL),dimension(5,nlength_seismogram),intent(in) :: seismogram_tmp_in

#ifdef USE_HDF5

  real(kind=CUSTOM_REAL), dimension(nlength_seismogram,1) :: seismogram_tmp
  logical, save :: is_initialized = .false.
  integer :: i, nlength_total_seismogram
  logical :: if_dataset_exists
  ! debugging
  logical, parameter :: DEBUG = .false.

  ! check if anything to do
  if (.not. OUTPUT_SEISMOS_HDF5) return

  ! safety check
  if (.not. WRITE_SEISMOGRAMS_BY_MAIN) &
    stop 'Error: WRITE_SEIMSMOGRAMS_BY_MAIN must be true to use HDF5 seismogram output'

  ! total length of the seismogram
  nlength_total_seismogram = NSTEP / NTSTEP_BETWEEN_OUTPUT_SAMPLE

  ! debug
  if (DEBUG .and. myrank == 0) then
    print *, '*** HDF5 seismogram outputs: station irec = ',irec,'***'
    print *, '    nlength_total_seismogram = ', nlength_total_seismogram
    print *, '    NSTEP = ', NSTEP
    print *, '    NTSTEP_BETWEEN_OUTPUT_SAMPLE = ', NTSTEP_BETWEEN_OUTPUT_SAMPLE
    print *, '    nrec = ', nrec
    print *, '    irec = ', irec, 'irec_local = ', irec_local
    print *, '    seismo_offset/seismo_current = ', seismo_offset,'/',seismo_current
    print *, '    shape(seismogram_tmp) = ', shape(seismogram_tmp)
    print *, '    shape(seismogram_tmp(iorientation,1:seismo_current)) = ', shape(seismogram_tmp(iorientation,1:seismo_current))
    print *
  endif

  ! convert  array with shape (seimo_current) to (nlength_total_seismogram, 1)
  do i = 1, seismo_current
    seismogram_tmp(i,1) = seismogram_tmp_in(iorientation,i)
  enddo

  ! initialize
  if (.not. is_initialized) then
    call write_hdf5_seismogram_init(nlength_total_seismogram)
    is_initialized = .true.
  endif

  ! only main process writes out
  if (myrank /= 0) return

  ! write the seismograms
  call h5_open_file(hdf5_seismo_fname)

  ! check if the target dataset components are already created
  call h5_check_dataset_exists(trim(chn),if_dataset_exists)

  ! if the dataset does not exist, create it
  if (.not. if_dataset_exists) then
    call h5_create_dataset_gen(trim(chn),(/nlength_total_seismogram, nrec/), 2, CUSTOM_REAL)
  endif

  ! write the seismogram
  call h5_write_dataset_collect_hyperslab(trim(chn), seismogram_tmp(1:seismo_current, :), (/seismo_offset,irec-1/), .false.)

  ! close the file
  call h5_close_file()

#else
  ! no HDF5 support

  ! to avoid compiler warnings
  integer :: idummy

  idummy = iorientation
  idummy = len_trim(chn)

  idummy = irec
  idummy = irec_local
  idummy = size(seismogram_tmp_in,kind=4)

  write(*,*) 'Error: HDF5 support not enabled in this version of Specfem3D_Globe'
  write(*,*) 'Please recompile with the --with-hdf5 option'
  stop

#endif

  end subroutine write_output_hdf5
