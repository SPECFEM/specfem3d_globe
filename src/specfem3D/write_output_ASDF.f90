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

!> Initializes the data structure for ASDF
!! \param nrec_local The number of receivers on the local processor
  subroutine init_asdf_data(nrec_local)

  use specfem_par, only: myrank

  use asdf_data, only: asdf_container

  implicit none
  ! Parameters
  integer,intent(in) :: nrec_local

  ! Variables
  integer :: total_seismos_local, ier

  total_seismos_local = nrec_local*3 ! 3 components

  asdf_container%nrec_local = nrec_local

  allocate (asdf_container%receiver_name_array(nrec_local), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%network_array(nrec_local), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%component_array(total_seismos_local), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%receiver_lat(nrec_local), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%receiver_lo(nrec_local), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%receiver_el(nrec_local), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%receiver_dpt(nrec_local), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%records(total_seismos_local), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')

  end subroutine init_asdf_data

!
!-------------------------------------------------------------------------------------------------
!

!> Stores the records into the ASDF structure
!! \param seismogram_tmp The current seismogram to store
!! \param irec_local The local index of the receiver
!! \param irec The global index of the receiver
!! \param chn The broadband channel simulated
!! \param iorientation The recorded seismogram's orientation direction
  subroutine store_asdf_data(seismogram_tmp, irec_local, irec, chn, iorientation)

  use constants, only: CUSTOM_REAL

  use specfem_par, only: &
    station_name,network_name, &
    seismo_current, NTSTEP_BETWEEN_OUTPUT_SEISMOS, myrank, &
    stlat, stlon, stele, stbur, &
    ROTATE_SEISMOGRAMS_RT ! BS BS

  use asdf_data, only: asdf_container

  implicit none

  ! Parameters
  character(len=4),intent(in) :: chn
  integer,intent(in) :: irec_local, irec, iorientation
  real(kind=CUSTOM_REAL),dimension(5,NTSTEP_BETWEEN_OUTPUT_SEISMOS),intent(in) :: seismogram_tmp

  ! local Variables
  integer :: length_station_name, length_network_name
  integer :: ier, i, index_increment


  ! BS BS: In case that components are rotated to radial and transverse, we need to
  ! change the way the trace index is incremented below. Otherwise, the code crashes
  ! in close_asdf_data() when trying to deallocate the records in asdf_container
  if (ROTATE_SEISMOGRAMS_RT) then
    index_increment = iorientation - 2
  else
    index_increment = iorientation
  endif

  ! trace index
  !i = (irec_local-1)*(3) + (iorientation)
  i = (irec_local-1)*(3) + (index_increment)

  length_station_name = len_trim(station_name(irec))
  length_network_name = len_trim(network_name(irec))

  asdf_container%receiver_name_array(irec_local) = station_name(irec)(1:length_station_name)
  asdf_container%network_array(irec_local) = network_name(irec)(1:length_network_name)
  asdf_container%component_array(i) = chn(1:3)
  asdf_container%receiver_lat(irec_local) = stlat(irec)
  asdf_container%receiver_lo(irec_local) = stlon(irec)
  asdf_container%receiver_el(irec_local) = stele(irec)
  asdf_container%receiver_dpt(irec_local) = stbur(irec)

  allocate (asdf_container%records(i)%record(seismo_current), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocating ASDF container failed.')

  asdf_container%records(i)%record(1:seismo_current) = seismogram_tmp(iorientation, 1:seismo_current)

  end subroutine store_asdf_data

!
!-------------------------------------------------------------------------------------------------
!

!> Closes the ASDF data structure by deallocating all arrays
  subroutine close_asdf_data()

  use asdf_data, only: asdf_container

  implicit none

  !Variables
  integer :: i

  do i = 1, asdf_container%nrec_local*3 ! 3 components
    deallocate(asdf_container%records(i)%record)
  enddo
  deallocate (asdf_container%receiver_name_array)
  deallocate (asdf_container%network_array)
  deallocate (asdf_container%component_array)

  end subroutine close_asdf_data

!
!-------------------------------------------------------------------------------------------------
!

!> Writes the ASDF data structure to the file
  subroutine write_asdf()

  use constants, only: OUTPUT_PROVENANCE
  use constants_solver, only: NDIM,itag

  use asdf_data, only: asdf_container

  use iso_c_binding
  use iso_Fortran_env

  use specfem_par

  !use my_mpi

  implicit none

  ! Parameters
  integer, parameter :: MAX_STRING_LENGTH = 1024
  integer, parameter :: MAX_QUAKEML_LENGTH = 8096
  integer, parameter :: MAX_STATIONXML_LENGTH = 16182
  integer, parameter :: MAX_PARFILE_LENGTH = 20000
  integer, parameter :: MAX_CONSTANTS_LENGTH = 45000
  integer, parameter :: MAX_TIME_STRING_LENGTH = 22

  !--- Character strings to be written to the ASDF file
  character(len=MAX_QUAKEML_LENGTH) :: quakeml
  character(len=MAX_STATIONXML_LENGTH) :: stationxml

  integer :: num_stations !
  integer :: stationxml_length
  integer :: nsamples  ! constant, as in SPECFEM
  double precision :: sampling_rate
  double precision :: startTime
  integer(kind=8) :: start_time

  ! Network names and station names are two different beast, as in SPECFEM
  ! network_names(i) is related to station_names(i)
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(:), allocatable :: networks_names
  character(len=MAX_LENGTH_STATION_NAME), dimension(:), allocatable :: stations_names
  character(len=3), dimension(:), allocatable :: component_names

  ! data. dimension = nsamples * num_channels_per_station * num_stations

  !-- ASDF variables
  !   These variables are used to know where further writes should be done.
  !   They have to be cleaned as soon as they become useless
  integer :: waveforms_grp  ! Group "/Waveforms/"
  integer, dimension(3) ::  data_ids   ! BS BS

  integer :: station_grp, stationxml_grp, current_proc, sender, receiver ! BS BS
  real (kind=CUSTOM_REAL), dimension(:,:), allocatable :: one_seismogram ! BS BS

  !--- MPI variables
  integer :: mysize, comm
  !--- Loop variables
  integer :: i, j, k, l
  !--- Error variable
  integer :: ier

  !--- 'allgather' arrays. Variables that needs to be known by everyone in
  !    order to define ASDF groups and datasets or write them as attributes.
  integer, dimension(:), allocatable :: num_stations_gather
  integer :: max_num_stations_gather
  character(len=MAX_LENGTH_STATION_NAME), dimension(:,:), allocatable :: &
      station_names_gather
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(:,:), allocatable :: &
      network_names_gather
  character(len=3), dimension(:,:), allocatable :: &
      component_names_gather
  real, dimension(:,:), allocatable :: &
      station_lats_gather, station_longs_gather, station_elevs_gather, &
      station_depths_gather
  integer, dimension(:), allocatable :: displs, rcounts

  ! temporary name built from network, station and channel names.
  character(len=MAX_STRING_LENGTH) :: waveform_name

  ! C/Fortran interop for C-allocated strings
  integer :: len_prov, len_constants, len_Parfile
  type(c_ptr) :: cptr
  character, pointer :: fptr(:)
  character, dimension(:), allocatable, TARGET :: provenance
  character(len=MAX_CONSTANTS_LENGTH) :: sf_constants
  character(len=MAX_PARFILE_LENGTH) :: sf_parfile

  ! Time variables
  character(len=MAX_TIME_STRING_LENGTH) :: start_time_string, end_time_string, &
                       cmt_start_time_string, pde_start_time_string

  ! alias MPI communicator
  call world_duplicate(comm)
  call world_size(mysize)

  num_stations = nrec_local
  sampling_rate = 1.0/DT
  !nsamples = seismo_current * (NSTEP / NTSTEP_BETWEEN_OUTPUT_SEISMOS)
  nsamples = NSTEP ! BS BS: The total number of samples to be written to the
                   !        ASDF file should be NSTEP.
                   !        seismo_current is equivalent to
                   !        NTSTEP_BETWEEN_OUTPUT_SEISMOS (except when
                   !        it==it_end), as only in this case
                   !        the routine is called from write_seismograms()

  ! Calculate start_time
  call get_time(startTime, start_time_string, pde_start_time_string, cmt_start_time_string, end_time_string)
  start_time = startTime*(int(1000000000,kind=8)) ! convert to nanoseconds

  !--------------------------------------------------------
  ! Setup data on each process.
  !--------------------------------------------------------

  ! Generate minimal QuakeML for SPECFEM3D_GLOBE
  call cmt_to_quakeml(quakeml, pde_start_time_string, cmt_start_time_string)

  if (OUTPUT_PROVENANCE) then
    ! Generate specfem provenance string
    call ASDF_generate_sf_provenance_f(trim(start_time_string)//C_NULL_CHAR, &
                                   trim(end_time_string)//C_NULL_CHAR, cptr, len_prov)
    call c_f_pointer(cptr, fptr, [len_prov])
    allocate(provenance(len_prov+1))
    provenance(1:len_prov) = fptr(1:len_prov)
    provenance(len_prov+1) = C_NULL_CHAR
  endif

  allocate(networks_names(num_stations), stat=ier)
  allocate(stations_names(num_stations), stat=ier)
  allocate(component_names(num_stations*3), stat=ier)

  !--------------------------------------------------------
  ! ASDF variables
  !--------------------------------------------------------
  ! Find how many stations are managed by each allgatheress
  allocate(num_stations_gather(mysize))
  call all_gather_all_i(num_stations, num_stations_gather, mysize)

  ! find the largest number of stations per allgatheress
  max_num_stations_gather = maxval(num_stations_gather)

  allocate(displs(mysize))
  allocate(rcounts(mysize))

  ! Everyone should know about each and every station name and its coordinates
  allocate(station_names_gather(max_num_stations_gather, mysize))
  allocate(network_names_gather(max_num_stations_gather, mysize))
  allocate(station_lats_gather(max_num_stations_gather,mysize))
  allocate(station_longs_gather(max_num_stations_gather,mysize))
  allocate(station_elevs_gather(max_num_stations_gather,mysize))
  allocate(station_depths_gather(max_num_stations_gather,mysize))
  allocate(component_names_gather(max_num_stations_gather*3, mysize))

  ! This needs to be done because asdf_data is a pointer
  do i = 1, num_stations
    write(networks_names(i), '(a)') asdf_container%network_array(i)
    write(stations_names(i), '(a)') asdf_container%receiver_name_array(i)
  enddo

  do i = 1, num_stations*3
    write(component_names(i), '(a)') asdf_container%component_array(i)
  enddo

  ! The number of stations is not constant across processes
  do i = 1, mysize
    displs(i) = (i-1) * max_num_stations_gather * MAX_LENGTH_STATION_NAME
    rcounts(i) = num_stations_gather(i) * MAX_LENGTH_STATION_NAME
  enddo

  call all_gather_all_ch(stations_names, &
                       num_stations * MAX_LENGTH_STATION_NAME, &
                       station_names_gather, &
                       rcounts, &
                       displs, &
                       max_num_stations_gather, &
                       MAX_LENGTH_STATION_NAME, &
                       mysize)


  do i = 1, mysize
    displs(i) = (i-1) * max_num_stations_gather * MAX_LENGTH_NETWORK_NAME
    rcounts(i) = num_stations_gather(i) * MAX_LENGTH_NETWORK_NAME
  enddo

  call all_gather_all_ch(networks_names, &
                       num_stations * MAX_LENGTH_NETWORK_NAME, &
                       network_names_gather, &
                       rcounts, &
                       displs, &
                       max_num_stations_gather, &
                       MAX_LENGTH_NETWORK_NAME, &
                       mysize)

  do i = 1, mysize
    displs(i) = (i-1) * max_num_stations_gather * 3
    rcounts(i) = num_stations_gather(i) * 3
  enddo

  call all_gather_all_ch(component_names, &
                       num_stations*3*3, &
                       component_names_gather, &
                       rcounts*3, &
                       displs*3, &
                       max_num_stations_gather*3, &
                       3, &
                       mysize)


  ! Now gather all the coordiante information for these stations
  do i = 1, mysize
    displs(i) = (i-1) * max_num_stations_gather
    rcounts(i) = num_stations_gather(i)
  enddo

  call all_gather_all_r(asdf_container%receiver_lat, &
                       num_stations, &
                       station_lats_gather, &
                       rcounts, &
                       displs, &
                       max_num_stations_gather, &
                       mysize)
  call all_gather_all_r(asdf_container%receiver_lo, &
                       num_stations, &
                       station_longs_gather, &
                       rcounts, &
                       displs, &
                       max_num_stations_gather, &
                       mysize)
  call all_gather_all_r(asdf_container%receiver_el, &
                       num_stations, &
                       station_elevs_gather, &
                       rcounts, &
                       displs, &
                       max_num_stations_gather, &
                       mysize)
  call all_gather_all_r(asdf_container%receiver_dpt, &
                       num_stations, &
                       station_depths_gather, &
                       rcounts, &
                       displs, &
                       max_num_stations_gather, &
                       mysize)

  deallocate(stations_names)
  deallocate(networks_names)
  deallocate(component_names)
  deallocate(displs)
  deallocate(rcounts)

  allocate(one_seismogram(NDIM,seismo_current),stat=ier)

  !--------------------------------------------------------
  ! write ASDF
  !--------------------------------------------------------

  ! we only want to do these steps one time
  if (seismo_offset == 0) then
    if (myrank == 0) then
        call ASDF_initialize_hdf5_f(ier);
        call ASDF_create_new_file_serial_f(trim(OUTPUT_FILES) // "synthetic.h5" // C_NULL_CHAR, current_asdf_handle)

        call ASDF_write_string_attribute_f(current_asdf_handle, "file_format" // C_NULL_CHAR, &
                                           "ASDF" // C_NULL_CHAR, ier)
        call ASDF_write_string_attribute_f(current_asdf_handle, "file_format_version" // C_NULL_CHAR, &
                                           "1.0.0" // C_NULL_CHAR, ier)

        call ASDF_write_quakeml_f(current_asdf_handle, trim(quakeml) // C_NULL_CHAR, ier)

        if (OUTPUT_PROVENANCE) then
          call ASDF_write_provenance_data_f(current_asdf_handle, provenance(1:len_prov+1), ier)

          call read_file("setup/constants.h", sf_constants, len_constants)
          call read_file("DATA/Par_file", sf_parfile, len_Parfile)


          call ASDF_write_auxiliary_data_f(current_asdf_handle, trim(sf_constants) // &
                                         C_NULL_CHAR, trim(sf_parfile(1:len_Parfile)) // &
                                         C_NULL_CHAR, ier)
        endif

        call ASDF_create_waveforms_group_f(current_asdf_handle, waveforms_grp)


        do k = 1, mysize ! Need to set up metadata for all processes

          do j = 1, num_stations_gather(k) ! loop over number of stations on that process
            call ASDF_create_stations_group_f(waveforms_grp, &
               trim(network_names_gather(j, k)) // "." //      &
               trim(station_names_gather(j, k)) // C_NULL_CHAR, &
               station_grp)
            stationxml_length = 1423 + len(trim(station_names_gather(j,k))) + len(trim(network_names_gather(j,k)))

            call ASDF_define_station_xml_f(station_grp, stationxml_length, &
                                         stationxml_grp)

            call station_to_stationxml(station_names_gather(j,k), network_names_gather(j,k), &
                                   station_lats_gather(j,k), station_longs_gather(j,k), &
                                   station_elevs_gather(j,k), station_depths_gather(j,k), &
                                   start_time_string, stationxml)
            call ASDF_write_station_xml_f(stationxml_grp, trim(stationxml)//C_NULL_CHAR, ier)

            do  i = 1, 3 ! loop over each component
              ! Generate unique waveform name
              write(waveform_name, '(a)') &
                trim(network_names_gather(j,k)) // "." // &
                trim(station_names_gather(j,k)) // ".S3." //trim(component_names_gather(i+(3*(j-1)),k)) &
                  //"__"//trim(start_time_string(1:19))//"__"//trim(end_time_string(1:19))//"__synthetic"
                call ASDF_define_waveform_f(station_grp, &
                  nsamples, start_time, sampling_rate, &
                  trim(event_name_SAC) // C_NULL_CHAR, &
                  trim(waveform_name) // C_NULL_CHAR, &
                  data_ids(i))
                call ASDF_close_dataset_f(data_ids(i), ier)
            enddo

            call ASDF_close_dataset_f(stationxml_grp, ier)
            call ASDF_close_group_f(station_grp, ier)

          enddo
        enddo

        call ASDF_close_group_f(waveforms_grp, ier)
        call ASDF_close_file_f(current_asdf_handle, ier)
        call ASDF_finalize_hdf5_f(ier)

    endif ! (myrank == 0)

  endif ! end (seismo_offset == 0) steps

  call synchronize_all()

  ! Now write waveforms
  if (WRITE_SEISMOGRAMS_BY_MASTER) then

    if (myrank == 0) then

      call ASDF_initialize_hdf5_f(ier);
      call ASDF_open_serial_f(trim(OUTPUT_FILES) // "synthetic.h5" // C_NULL_CHAR, current_asdf_handle)
      call ASDF_open_waveforms_group_f(current_asdf_handle, waveforms_grp)

    endif

    do k = 1, mysize ! Need to write data from all processes

        current_proc = k - 1
        sender=current_proc
        receiver=0 ! the master proc does all the writing

      do j = 1, num_stations_gather(k) ! loop over number of stations on that process

        l = (j-1)*(NDIM) ! Index of current receiver in asdf_container%records

        ! First get the information to the master proc
        if (current_proc == 0) then ! current_proc is master proc

          !one_seismogram(:,:) = seismograms(:,j,:)
          if (myrank == 0) then
            do i = 1, NDIM
            !  write(*,*) j, l, l+i, size(asdf_container%records)
              one_seismogram(i,:) = asdf_container%records(l+i)%record(1:seismo_current)
            enddo
          endif

        else ! current_proc is not master proc

          if (myrank == current_proc) then

            !one_seismogram(:,:) = seismograms(:,j,:)
            do i = 1, NDIM
              one_seismogram(i,:) = asdf_container%records(l+i)%record(1:seismo_current)
            enddo

            call send_cr(one_seismogram,NDIM*seismo_current,receiver,itag)

          else if (myrank == 0) then

            call recv_cr(one_seismogram,NDIM*seismo_current,sender,itag)

          endif
        endif

        ! Now do the actual writing
        if (myrank == 0) then

          call ASDF_open_stations_group_f(waveforms_grp, &
            trim(network_names_gather(j, k)) // "." //      &
            trim(station_names_gather(j, k)) // C_NULL_CHAR, &
            station_grp)

          do  i = 1, NDIM ! loop over each component
            ! Generate unique waveform name
            write(waveform_name, '(a)') &
              trim(network_names_gather(j,k)) // "." // &
              trim(station_names_gather(j,k)) // ".S3." //trim(component_names_gather(i+(3*(j-1)),k)) &
                //"__"//trim(start_time_string(1:19))//"__"//trim(end_time_string(1:19))//"__synthetic"

            call ASDF_open_waveform_f(station_grp, &
              trim(waveform_name) // C_NULL_CHAR, &
              data_ids(i))

            call ASDF_write_partial_waveform_f(data_ids(i), &
                                        one_seismogram(i,1:seismo_current), seismo_offset, seismo_current, ier)
            call ASDF_close_dataset_f(data_ids(i), ier)

          enddo
          call ASDF_close_group_f(station_grp, ier)
        endif

      enddo
    enddo

    if (myrank == 0) then

      call ASDF_close_group_f(waveforms_grp, ier)
      call ASDF_close_file_f(current_asdf_handle, ier)
      call ASDF_finalize_hdf5_f(ier)

    endif

  else ! write seismograms in parallel

    call ASDF_initialize_hdf5_f(ier);
    call ASDF_open_f(trim(OUTPUT_FILES) // "synthetic.h5" // C_NULL_CHAR, comm, current_asdf_handle)

    call ASDF_open_waveforms_group_f(current_asdf_handle, waveforms_grp)

    do k = 1, mysize ! Need to open ASDF groups on all processes
      do j = 1, num_stations_gather(k) ! loop over number of stations on that process
        call ASDF_open_stations_group_f(waveforms_grp, &
          trim(network_names_gather(j, k)) // "." //      &
          trim(station_names_gather(j, k)) // C_NULL_CHAR, &
          station_grp)

        l = (j-1)*(NDIM) ! Index of current receiver in asdf_container%records

        do  i = 1, NDIM ! loop over each component

          ! Generate unique waveform name
          write(waveform_name, '(a)') &
            trim(network_names_gather(j,k)) // "." // &
            trim(station_names_gather(j,k)) // ".S3." //trim(component_names_gather(i+(3*(j-1)),k)) &
              //"__"//trim(start_time_string(1:19))//"__"//trim(end_time_string(1:19))//"__synthetic"

          call ASDF_open_waveform_f(station_grp, &
            trim(waveform_name) // C_NULL_CHAR, &
            data_ids(i))

          if (k == myrank+1) then

            one_seismogram(i,:) = asdf_container%records(l+i)%record(1:seismo_current)

            call ASDF_write_partial_waveform_f(data_ids(i), &
                                        one_seismogram(i,1:seismo_current), seismo_offset, seismo_current, ier)
          endif

          call ASDF_close_dataset_f(data_ids(i), ier)
        enddo
        call ASDF_close_group_f(station_grp, ier)
      enddo
    enddo

    call synchronize_all()


    call ASDF_close_group_f(waveforms_grp, ier)
    call ASDF_close_file_f(current_asdf_handle, ier)
    call ASDF_finalize_hdf5_f(ier)

  endif ! WRITE_SEISMOGRAMS_BY_MASTER


  !--------------------------------------------------------
  ! Clean up
  !--------------------------------------------------------


  if (OUTPUT_PROVENANCE) deallocate(provenance)
  deallocate(station_names_gather)
  deallocate(network_names_gather)
  deallocate(component_names_gather)
  deallocate(station_lats_gather)
  deallocate(station_longs_gather)
  deallocate(station_elevs_gather)
  deallocate(station_depths_gather)
  deallocate(num_stations_gather)

  deallocate(one_seismogram)

  end subroutine write_asdf

!
!-------------------------------------------------------------------------------------------------
!

!> Converts the CMT source file read by SPECFEM to a QuakeML file for the ASDF container
!! \param quakemlstring The QuakeML string to store in the ASDF file
!! \start_time_string The start date stored as a character string
  subroutine cmt_to_quakeml(quakemlstring, pde_start_time_string, cmt_start_time_string)

  use specfem_par, only: &
    cmt_lat => cmt_lat_SAC,cmt_lon => cmt_lon_SAC,cmt_depth => cmt_depth_SAC, &
    hdur => cmt_hdur_SAC,M0,Mrr,Mtt,Mpp,Mrt,Mrp,Mtp,event_name_SAC, &
    pde_lat => elat_SAC,pde_lon => elon_SAC,pde_depth => depth_SAC, &
    mb => mb_SAC,ms => ms_SAC,Mw

  implicit none
  character(len=*) :: quakemlstring
  character(len=*) :: pde_start_time_string
  character(len=*) :: cmt_start_time_string
  character(len=13) :: cmt_lon_str, cmt_lat_str, cmt_depth_str, hdur_str
  character(len=13) :: pde_lat_str, pde_lon_str, pde_depth_str
  character(len=25) :: M0_str, mb_str, ms_str, Mw_str
  character(len=25) :: Mrr_str, Mtt_str, Mpp_str, Mrt_str, Mrp_str, Mtp_str

  ! Convert the CMT values to strings for the QuakeML string
  write(pde_lat_str, "(g12.5)") pde_lat
  write(pde_lon_str, "(g12.5)") pde_lon
  write(pde_depth_str, "(g12.5)") pde_depth*1000 ! km to m conversion
  write(cmt_lat_str, "(g12.5)") cmt_lat
  write(cmt_lon_str, "(g12.5)") cmt_lon
  write(cmt_depth_str, "(g12.5)") cmt_depth*1000 ! km to m conversion
  write(hdur_str, "(g12.5)") hdur*2
  write(M0_str, "(g12.5)") M0*1e-7 ! dyn-cm to N-m conversion
  write(mb_str, "(g12.5)") mb
  write(ms_str, "(g12.5)") ms
  write(Mw_str, "(g12.5)") Mw
  write(Mrr_str, "(g12.5)") Mrr*1e-7
  write(Mtt_str, "(g12.5)") Mtt*1e-7
  write(Mpp_str, "(g12.5)") Mpp*1e-7
  write(Mrt_str, "(g12.5)") Mrt*1e-7
  write(Mrp_str, "(g12.5)") Mrp*1e-7
  write(Mtp_str, "(g12.5)") Mtp*1e-7

  quakemlstring = '<q:quakeml xmlns="http://quakeml.org/xmlns/bed/1.2"'//&
                  ' xmlns:q="http://quakeml.org/xmlns/quakeml/1.2">'//&
                  '<eventParameters publicID="smi:local/'//trim(event_name_SAC)//'#eventPrm">'//&
                  '<event publicID="smi:local/'//trim(event_name_SAC)//'#eventID">'//&
                  '<preferredOriginID>smi:local/'//trim(event_name_SAC)//'/origin#cmtorigin</preferredOriginID>'//&
                  '<preferredMagnitudeID>smi:local/'//trim(event_name_SAC)//'/magnitude#moment_mag</preferredMagnitudeID>'//&
                  '<preferredFocalMechanismID>smi:local/'//trim(event_name_SAC)//'/focal_mechanism</preferredFocalMechanismID>'//&
                  '<type>earthquake</type>'//&
                  '<typeCertainty>known</typeCertainty>'//&
                  '<description>'//&
                  '<text>'//trim(event_name_SAC)//'</text>'//&
                  '<type>earthquake name</type>'//&
                  '</description>'//&
                  '<origin publicID="smi:local/'//trim(event_name_SAC)//'/origin#reforigin">'//&
                  '<time>'//&
                  '<value>'//trim(pde_start_time_string)//'</value>'//&
                  '</time>'//&
                  '<latitude>'//&
                  '<value>'//trim(pde_lat_str)//'</value>'//&
                  '</latitude>'//&
                  '<longitude>'//&
                  '<value>'//trim(pde_lon_str)//'</value>'//&
                  '</longitude>'//&
                  '<depth>'//&
                  '<value>'//trim(pde_depth_str)//'</value>'//&
                  '</depth>'//&
                  '<type>hypocenter</type>'//&
                  '<comment id="smi:local/'//trim(event_name_SAC)//'/comment#ref_origin">'//&
                  '<text>Hypocenter catalog: PDE</text>'//&
                  '</comment>'//&
                  '</origin>'//&
                  '<origin publicID="smi:local/'//trim(event_name_SAC)//'/origin#cmtorigin">'//&
                  '<time>'//&
                  '<value>'//trim(cmt_start_time_string)//'</value>'//&
                  '</time>'//&
                  '<latitude>'//&
                  '<value>'//trim(cmt_lat_str)//'</value>'//&
                  '</latitude>'//&
                  '<longitude>'//&
                  '<value>'//trim(cmt_lon_str)//'</value>'//&
                  '</longitude>'//&
                  '<depth>'//&
                  '<value>'//trim(cmt_depth_str)//'</value>'//&
                  '</depth>'//&
                  '</origin>'//&
                  '<focalMechanism publicID="smi:local/'//trim(event_name_SAC)//'/focal_mechanism">'//&
                  '<momentTensor publicID="smi:local/'//trim(event_name_SAC)//'/momenttensor">'//&
                  '<derivedOriginID>smi:local/'//trim(event_name_SAC)//'/origin#cmtorigin'//&
                  '</derivedOriginID>'//&
                  '<momentMagnitudeID>smi:local/'//trim(event_name_SAC)//'/magnitude#moment_mag'//&
                  '</momentMagnitudeID>'//&
                  '<scalarMoment>'//&
                  '<value>'//trim(M0_str)//'</value>'//&
                  '</scalarMoment>'//&
                  '<tensor>'//&
                  '<Mrr>'//&
                  '<value>'//trim(Mrr_str)//'</value>'//&
                  '<uncertainty>0</uncertainty>'//&
                  '</Mrr>'//&
                  '<Mtt>'//&
                  '<value>'//trim(Mtt_str)//'</value>'//&
                  '<uncertainty>0</uncertainty>'//&
                  '</Mtt>'//&
                  '<Mpp>'//&
                  '<value>'//trim(Mpp_str)//'</value>'//&
                  '<uncertainty>0</uncertainty>'//&
                  '</Mpp>'//&
                  '<Mrt>'//&
                  '<value>'//trim(Mrt_str)//'</value>'//&
                  '<uncertainty>0</uncertainty>'//&
                  '</Mrt>'//&
                  '<Mrp>'//&
                  '<value>'//trim(Mrp_str)//'</value>'//&
                  '<uncertainty>0</uncertainty>'//&
                  '</Mrp>'//&
                  '<Mtp>'//&
                  '<value>'//trim(Mtp_str)//'</value>'//&
                  '<uncertainty>0</uncertainty>'//&
                  '</Mtp>'//&
                  '</tensor>'//&
                  '<sourceTimeFunction>'//&
                  '<type>triangle</type>'//&
                  '<duration>'//trim(hdur_str)//'</duration>'//&
                  '</sourceTimeFunction>'//&
                  '</momentTensor>'//&
                  '</focalMechanism>'//&
                  '<magnitude publicID="smi:local/'//trim(event_name_SAC)//'/magnitude#moment_mag">'//&
                  '<mag>'//&
                  '<value>'//trim(Mw_str)//'</value>'//&
                  '</mag>'//&
                  '<type>Mwc</type>'//&
                  '</magnitude>'//&
                  '<magnitude publicID="smi:local/'//trim(event_name_SAC)//'/magnitude#mb">'//&
                  '<mag>'//&
                  '<value>'//trim(mb_str)//'</value>'//&
                  '</mag>'//&
                  '<type>mb</type>'//&
                  '</magnitude>'//&
                  '<magnitude publicID="smi:local/'//trim(event_name_SAC)//'/magnitude#MS">'//&
                  '<mag>'//&
                  '<value>'//trim(ms_str)//'</value>'//&
                  '</mag>'//&
                  '<type>MS</type>'//&
                  '</magnitude>'//&
                  '</event>'//&
                  '</eventParameters>'//&
                  '</q:quakeml>'

  end subroutine cmt_to_quakeml

!
!-------------------------------------------------------------------------------------------------
!

! Convert system time into a string
!! \param time system time
!! \param time_string a string representation of the input time
  subroutine convert_systime_to_string(time, time_string)

  implicit none
  double precision, intent(in) :: time
  character(len=*), intent(out) :: time_string
  real :: fraction_sec
  integer :: iatime(9)
  character(len=4) :: yr
  character(len=2) :: mo, da, hr, minute
  character(len=15) :: second
  real :: real_sec

  ! extract msec
  fraction_sec = time - int(time)

  call gmtime(int(time), iatime)
  write(yr, "(I4.4)") iatime(6) + 1900
  write(mo, "(I2.2)") iatime(5) + 1
  write(da, "(I2.2)") iatime(4)
  write(hr, "(I2.2)") iatime(3)
  write(minute, "(I2.2)") iatime(2)

  real_sec = iatime(1) + fraction_sec
  write(second, "(I2.2, F0.4)") int(real_sec), real_sec-int(real_sec)

  time_string = trim(yr)//"-"//trim(mo)//"-"//trim(da)//"T"//&
                  trim(hr)//':'//trim(minute)//':'//trim(second)

  end subroutine convert_systime_to_string

!
!-------------------------------------------------------------------------------------------------
!

!> Uses the time in the CMTSOLUTION file to calculate the number of seconds since the epoch
!! \param starttime The start time of the simulation from the epoch
!! \param start_time_string A string for defining the waveform name start time
!! \param pde_start_time_string A string for defining the waveform name start time using PDE
!! \param cmt_start_time_string A string for defining the waveform name start time using CMT
!! \param end_time_string A string for defining the waveform name end time
  subroutine get_time(starttime, start_time_string, pde_start_time_string, cmt_start_time_string, end_time_string)

  use specfem_par, only: &
    NSTEP,DT,t_shift_SAC,hdur, &
    yr_SAC,ho_SAC,mi_SAC,sec_SAC,jda_SAC

  implicit none
  character(len=*) :: start_time_string
  character(len=*) :: end_time_string
  character(len=*) :: cmt_start_time_string
  character(len=*) :: pde_start_time_string
  double precision, intent(inout) :: starttime
  double precision :: trace_length_in_sec
  integer :: year
  double precision :: pdetime, cmttime, endtime
  integer,dimension(8) :: values

  ! Calculates the start time since the epoch in seconds
  ! Reference:
  ! http://pubs.opengroup.org/onlinepubs/009695399/basedefs/xbd_chap04.html#tag_04_14
  year = yr_SAC-1900
  pdetime = (year-70)*31536000.0d0+((year-69)/4)*86400.0d0 -((year-1)/100)*86400.0d0+&
              ((year+299)/400)*86400.0d0+(jda_SAC-1)*86400.0d0+ho_SAC*(3600.0d0)+&
              mi_SAC*60.0d0+sec_SAC
  call convert_systime_to_string(pdetime, pde_start_time_string)

  ! cmt centroid time
  cmttime = pdetime + t_shift_SAC
  call convert_systime_to_string(cmttime, cmt_start_time_string)

  ! trace start time
  starttime = cmttime - 1.5*hdur(1)
  call convert_systime_to_string(starttime, start_time_string)

  ! Calculates the number of seconds to add to the start_time
  trace_length_in_sec = DT*NSTEP
  endtime = starttime + trace_length_in_sec
  call convert_systime_to_string(endtime, end_time_string)

  ! Calculates time in seconds since the epoch
  call date_and_time(VALUES=values)
  !time8()

  !write(*,*) fmtdate(values,'The CPU time used by this program is now %c seconds')

  end subroutine get_time

!
!-------------------------------------------------------------------------------------------------
!

!> Converts the Station information to a StationXML string
!! \param station_name The name of the station based on SEED
!! \param network_name The name of the network based on SEED
!! \param stationxmlstring The StationXML string that will be written to the ASDF file
  subroutine station_to_stationxml(station_name, network_name, latitude, longitude, elevation, &
                                   burial_depth, start_time_string, stationxmlstring)

  implicit none
  character(len=*) :: stationxmlstring
  character(len=*) :: station_name
  character(len=*) :: network_name
  character(len=*) :: start_time_string
  character(len=15) :: station_lat, station_lon, station_ele, station_depth
  integer :: len_station_name, len_network_name, len_station_lat
  integer :: len_station_lon, len_station_depth, len_station_ele
  real, intent(in) :: latitude, longitude, elevation, burial_depth

  ! Convert double precision to character strings for the StationXML string

  write(station_lat, "(F11.4)") latitude
  write(station_lon, "(F11.4)") longitude
  write(station_ele, "(F8.1)") elevation
  write(station_depth, "(F8.1)") burial_depth

  ! print *, trim(station_lat), trim(station_lon), trim(station_depth), trim(station_ele)

  len_network_name = len(network_name)
  len_station_name = len(station_name)
  len_station_lat = len(trim(station_lat))
  len_station_lon = len(trim(station_lon))
  len_station_depth = len(trim(station_depth))
  len_station_ele = len(trim(station_ele))

  stationxmlstring = '<FDSNStationXML schemaVersion="1.0" xmlns="http://www.fdsn.org/xml/station/1">'//&
                     '<Source>SPECFEM3D_GLOBE</Source>'//&
                     '<Module>SPECFEM3D_GLOBE/asdf-library</Module>'//&
                     '<ModuleURI>http://seismic-data.org</ModuleURI>'//&
                     '<Created>'//trim(start_time_string)//'</Created>'//&
                     '<Network code="'//trim(network_name(1:len(network_name)))//'"'//&
                     '><Station code="'//trim(station_name(1:len(station_name)))//'">'//&
                     '<Latitude unit="DEGREES">'//trim(station_lat(1:len_station_lat))//'</Latitude>'//&
                     '<Longitude unit="DEGREES">'//trim(station_lon(1:len_station_lon))//'</Longitude>'//&
                     '<Elevation>'//trim(station_ele(1:len_station_ele))//'</Elevation>'//&
                     '<Site>'//&
                     '<Name>N/A</Name>'//&
                     '</Site>'//&
                     '<CreationDate>'//trim(start_time_string)//'</CreationDate>'//&
                     '<TotalNumberChannels>3</TotalNumberChannels>'//&
                     '<SelectedNumberChannels>3</SelectedNumberChannels>'//&
                     '<Channel locationCode="S3" code="MXN"'//&
                     ' startDate="'//trim(start_time_string)//'">'//&
                     '<Latitude unit="DEGREES">'//trim(station_lat(1:len_station_lat))//'</Latitude>'//&
                     '<Longitude unit="DEGREES">'//trim(station_lon(1:len_station_lon))//'</Longitude>'//&
                     '<Elevation>'//trim(station_ele(1:len_station_ele))//'</Elevation>'//&
                     '<Depth>'//trim(station_depth(1:len_station_depth))//'</Depth>'//&
                     '<Azimuth>0.0</Azimuth>'//&
                     '<Dip>0.0</Dip>'//&
                     '</Channel>'//&
                     '<Channel locationCode="S3" code="MXE"'//&
                     ' startDate="'//trim(start_time_string)//'">'//&
                     '<Latitude unit="DEGREES">'//trim(station_lat(1:len_station_lat))//'</Latitude>'//&
                     '<Longitude unit="DEGREES">'//trim(station_lon(1:len_station_lon))//'</Longitude>'//&
                     '<Elevation>'//trim(station_ele(1:len_station_ele))//'</Elevation>'//&
                     '<Depth>'//trim(station_depth(1:len_station_depth))//'</Depth>'//&
                     '<Azimuth>90.0</Azimuth>'//&
                     '<Dip>0.0</Dip>'//&
                     '</Channel>'//&
                     '<Channel locationCode="S3" code="MXZ"'//&
                     ' startDate="'//trim(start_time_string)//'">'//&
                     '<Latitude unit="DEGREES">'//trim(station_lat(1:len_station_lat))//'</Latitude>'//&
                     '<Longitude unit="DEGREES">'//trim(station_lon(1:len_station_lon))//'</Longitude>'//&
                     '<Elevation>'//trim(station_ele(1:len_station_ele))//'</Elevation>'//&
                     '<Depth>'//trim(station_depth(1:len_station_depth))//'</Depth>'//&
                     '<Azimuth>0.0</Azimuth>'//&
                     '<Dip>90.0</Dip>'//&
                     '</Channel>'//&
                     '</Station>'//&
                     '</Network>'//&
                     '</FDSNStationXML>'

  end subroutine station_to_stationxml

!
!-------------------------------------------------------------------------------------------------
!

!> Reads an external file and stores it in filestring
!! \param filename The name of the file to read
!! \param filestring The string that the file is stored
  subroutine read_file(filename, filestring, filesize)

  implicit none
  character(len=*) :: filestring
  character(len=*) :: filename
  integer,intent(out) :: filesize

  ! Get the size of the file using Fortran2003 feature
  open(10, file=filename, status='old')
  inquire(unit=10, size=filesize)
  close(10)

  ! Read in the size of the file using direct access
  open(10, file=filename, status='old', &
         recl=filesize, form='unformatted', access='direct')
  read (10, rec=1) filestring(1:filesize)
  close(10)

  end subroutine read_file
