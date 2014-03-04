!--------------------------------------------------------------------------------------------------------
!> \file write_output_ASDF.F90
!! \brief Write subroutines for writing ASDF seismograms to file using the ADIOS library
!! \author JAS and Wenjie Lei
!--------------------------------------------------------------------------------------------------------

#include "config.fh"

!> Initializes the data structure for asdf
!! \param my_asdf The asdf data structure
!! \param total_seismos_local The number of records on the local processer
subroutine init_asdf_data(my_asdf,total_seismos_local)

  use asdf_data
  use specfem_par, only : event_name_SAC,myrank
  type(asdf_event) :: my_asdf
  integer :: total_seismos_local

  my_asdf%nrecords = total_seismos_local
  my_asdf%event = trim(event_name_SAC)
  allocate (my_asdf%npoints(my_asdf%nrecords))
  allocate (my_asdf%gmt_year(my_asdf%nrecords))
  allocate (my_asdf%gmt_hour(my_asdf%nrecords))
  allocate (my_asdf%gmt_day(my_asdf%nrecords))
  allocate (my_asdf%gmt_min(my_asdf%nrecords))
  allocate (my_asdf%gmt_sec(my_asdf%nrecords))
  allocate (my_asdf%gmt_msec(my_asdf%nrecords))
  allocate (my_asdf%event_lat(my_asdf%nrecords))
  allocate (my_asdf%event_lo(my_asdf%nrecords))
  allocate (my_asdf%event_dpt(my_asdf%nrecords))
  allocate (my_asdf%receiver_lat(my_asdf%nrecords))
  allocate (my_asdf%receiver_lo(my_asdf%nrecords))
  allocate (my_asdf%receiver_el(my_asdf%nrecords))
  allocate (my_asdf%receiver_dpt(my_asdf%nrecords))
  allocate (my_asdf%begin_value(my_asdf%nrecords))
  allocate (my_asdf%end_value(my_asdf%nrecords))
  allocate (my_asdf%cmp_azimuth(my_asdf%nrecords))
  allocate (my_asdf%cmp_incident_ang(my_asdf%nrecords))
  allocate (my_asdf%sample_rate(my_asdf%nrecords))
  allocate (my_asdf%scale_factor(my_asdf%nrecords))
  allocate (my_asdf%ev_to_sta_AZ(my_asdf%nrecords))
  allocate (my_asdf%sta_to_ev_AZ(my_asdf%nrecords))
  allocate (my_asdf%great_circle_arc(my_asdf%nrecords))
  allocate (my_asdf%dist(my_asdf%nrecords))
  allocate (my_asdf%P_pick(my_asdf%nrecords))
  allocate (my_asdf%S_pick(my_asdf%nrecords))
  allocate (my_asdf%records(my_asdf%nrecords))
  allocate (my_asdf%receiver_name_array(my_asdf%nrecords))
  allocate (my_asdf%network_array(my_asdf%nrecords))
  allocate (my_asdf%component_array(my_asdf%nrecords))
  allocate (my_asdf%receiver_id_array(my_asdf%nrecords))

end subroutine init_asdf_data

!> Stores the records into the asdf structure
!! \param my_asdf The asdf data structure
!! \param seismogram_tmp The current seismogram to store
!! \param irec_local The local index of the receivers on the local processor
!! \param irec The global index of the receiver
!! \param chn The broadband channel simulated
!! \param iorientation The recorded seismogram's orientation direction
subroutine store_asdf_data(my_asdf,seismogram_tmp,irec_local,irec,chn,iorientation)

  use asdf_data
  use specfem_par,only: &
          station_name,network_name,stlat,stlon,stele,stbur, &
          DT,t0, &
          seismo_offset,seismo_current,it_end, &
          NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          yr=>yr_SAC,jda=>jda_SAC,ho=>ho_SAC,mi=>mi_SAC,sec=>sec_SAC, &
          tshift_cmt=>t_cmt_SAC,t_shift=>t_shift_SAC, &
          elat=>elat_SAC,elon=>elon_SAC,depth=>depth_SAC, &
          event_name=>event_name_SAC,cmt_lat=>cmt_lat_SAC,cmt_lon=>cmt_lon_SAC,&
          cmt_depth=>cmt_depth_SAC,cmt_hdur=>cmt_hdur_SAC

  implicit none
  include "constants.h"
  character(len=4) :: chn
  integer :: i, irec_local, irec
  integer :: length_station_name, length_network_name
  double precision, allocatable :: trace(:)
  real(kind=CUSTOM_REAL),dimension(5,NTSTEP_BETWEEN_OUTPUT_SEISMOS) :: seismogram_tmp
  real :: B
  integer iorientation
  integer :: adios_err
  type(asdf_event) :: my_asdf

  i = (irec_local-1)*(3) + (iorientation)
  my_asdf%npoints(i) = seismo_current
  my_asdf%gmt_year(i) = yr
  my_asdf%gmt_day(i) = jda
  my_asdf%gmt_hour(i) = ho
  my_asdf%gmt_min(i) = mi
  my_asdf%gmt_sec(i) = sec
  my_asdf%gmt_msec(i) = 0
  my_asdf%event_lat(i) = cmt_lat
  my_asdf%event_lo(i) = cmt_lon
  my_asdf%event_dpt(i) = cmt_depth
  my_asdf%receiver_lat(i) = stlat(irec_local)
  my_asdf%receiver_lo(i) = stlon(irec_local)
  my_asdf%receiver_el(i) = stele(irec_local)
  my_asdf%receiver_dpt(i) = stbur(irec_local)
  my_asdf%begin_value(i) = seismo_offset*DT-t0+tshift_cmt 
  my_asdf%end_value(i) = -12345
  my_asdf%cmp_azimuth(i) = 0.0
  my_asdf%cmp_incident_ang(i) = 0.0
  my_asdf%sample_rate(i) = DT
  my_asdf%scale_factor(i) = -12345
  my_asdf%ev_to_sta_AZ(i) = -12345
  my_asdf%sta_to_ev_AZ(i) = -12345
  my_asdf%great_circle_arc(i) = -12345
  my_asdf%dist(i) = -12345
  my_asdf%P_pick(i) = -12345
  my_asdf%S_pick(i) = -12345
  length_station_name = len_trim(station_name(irec))
  length_network_name = len_trim(network_name(irec))
  my_asdf%receiver_name_array(i) = station_name(irec)(1:length_station_name)
  my_asdf%network_array(i) = network_name(irec)(1:length_network_name)
  my_asdf%component_array(i) = chn
  my_asdf%receiver_id_array(i) = ""
  allocate (my_asdf%records(i)%record(seismo_current))
  my_asdf%records(i)%record(1:seismo_current) = seismogram_tmp(iorientation, 1:seismo_current)

end subroutine store_asdf_data

!> Closes the asdf data structure by deallocating all arrays
!! \param my_asdf The asdf data structure
!! \param total_seismos_local The number of seismograms on the local processor
subroutine close_asdf_data(my_asdf, total_seismos_local)

  use asdf_data
  type(asdf_event) :: my_asdf
  integer :: i, total_seismos_local

  deallocate (my_asdf%npoints)
  deallocate (my_asdf%gmt_year)
  deallocate (my_asdf%gmt_hour)
  deallocate (my_asdf%gmt_day)
  deallocate (my_asdf%gmt_min)
  deallocate (my_asdf%gmt_sec)
  deallocate (my_asdf%gmt_msec)
  deallocate (my_asdf%event_lat)
  deallocate (my_asdf%event_lo)
  deallocate (my_asdf%event_dpt)
  deallocate (my_asdf%receiver_lat)
  deallocate (my_asdf%receiver_lo)
  deallocate (my_asdf%receiver_el)
  deallocate (my_asdf%receiver_dpt)
  deallocate (my_asdf%begin_value)
  deallocate (my_asdf%end_value)
  deallocate (my_asdf%cmp_azimuth)
  deallocate (my_asdf%cmp_incident_ang)
  deallocate (my_asdf%sample_rate)
  deallocate (my_asdf%scale_factor)
  deallocate (my_asdf%ev_to_sta_AZ)
  deallocate (my_asdf%sta_to_ev_AZ)
  deallocate (my_asdf%great_circle_arc)
  deallocate (my_asdf%dist)
  deallocate (my_asdf%P_pick)
  deallocate (my_asdf%S_pick)
  do i = 1, total_seismos_local
    deallocate(my_asdf%records(i)%record)
  enddo 
  deallocate (my_asdf%receiver_name_array)
  deallocate (my_asdf%network_array)
  deallocate (my_asdf%component_array)
  deallocate (my_asdf%receiver_id_array)

end subroutine close_asdf_data

!> Writes the asdf data structure to the file
!! \param my_asdf The asdf data structure
subroutine write_asdf(my_asdf)

  use asdf_data
  use adios_write_mod
  use specfem_par, only : event_name_SAC,myrank

  implicit none
  integer :: adios_err, comm, ierr, sizeprocs
  integer(kind=8) :: adios_group
  character(len=200) :: ASDF_FN
  type(asdf_event) :: my_asdf

  call world_duplicate(comm)
  call world_size(sizeprocs)

  ! declare new group that uses MPI
  call adios_declare_group (adios_group, "EVENTS", "iter", 1, adios_err)
  call adios_select_method (adios_group, "MPI", "", "", adios_err)
  
  ASDF_FN="OUTPUT_FILES/"//trim(event_name_SAC)//"_sem.bp"
  call write_asdf_data (ASDF_FN, my_asdf, adios_group, myrank, sizeprocs, comm, ierr)

end subroutine write_asdf

!> Writes the asdf data structure to asdf_fn using parallel write
!! \param asdf_fn The file name for asdf
!! \param my_asdf The asdf data structure
!! \param adios_group The adios group for the file
!! \param rank The rank of the processor
!! \param nproc The number of processors
!! \param comm The communication group of processors
!! \param ierr The error for adios subroutine calls
subroutine write_asdf_data(asdf_fn, my_asdf, adios_group, rank, nproc, comm, ierr)

  use asdf_data
  use adios_write_mod

  character(len=*) :: asdf_fn
  type(asdf_event),intent(in) :: my_asdf
  integer :: rank, nproc, comm, ierr

  integer         :: adios_err
  integer(kind=8) :: adios_groupsize, adios_totalsize, varid
  integer(kind=8) :: adios_handle, adios_group

  !calculate size
  adios_groupsize = 0
  call define_asdf_data (adios_group, adios_groupsize, my_asdf,&
                                                rank, nproc, comm, ierr)
  call adios_open (adios_handle, "EVENTS", asdf_fn, "w", comm, adios_err)
  call adios_group_size (adios_handle, adios_groupsize, adios_totalsize, adios_err)

  !call the write sub
  call write_asdf_data_sub (my_asdf, adios_handle, adios_group,&
                                adios_groupsize, rank, nproc, comm, ierr)

  !adios close
  call adios_close(adios_handle, adios_err)

end subroutine write_asdf_data

!> Defines the asdf structure using adios
!! \param adios_group The adios group
!! \param my_group_size The adios group size
!! \param my_asdf The asdf data structure
!! \param rank The rank of the processor
!! \param nproc The number of processors
!! \param comm The communication group of processors
!! \param ierr The error for adios subroutine calls
subroutine define_asdf_data (adios_group, my_group_size, my_asdf, &
                                         rank, nproc, comm, ierr)

  use adios_write_mod
  use asdf_helpers_mod
  use asdf_data
  use specfem_par,only: nrec
  implicit none

  integer(kind=8), intent(in) :: adios_group
  integer(kind=8) :: my_group_size
  type(asdf_event) :: my_asdf
  integer, intent(in) :: rank, nproc, comm
  integer :: ierr

  integer :: i, nerr, string_total_length
  integer, parameter :: STRING_COMMON_LENGTH = 20
  integer :: adios_err, stat

  integer(kind=8) :: varid

  integer :: nrecords

  character(len=2)             :: data_type
  character(len=32)            :: header, record
  character(len=6)             :: npts_string
  character(len=10)            :: i_string
  character(len=200)           :: command, dummy, record_path

  integer :: dum_int, int_array(10)
  real    :: dum_real, real_array(10)
  character(len=10) :: dum_string

  integer :: nrecords_total, offset

  !gather info. Here, we only need nrecords_total
  nrecords=my_asdf%nrecords
  call gather_offset_info(nrecords,nrecords_total,offset,rank,nproc,comm,ierr)

  call define_adios_local_string_1d_array (adios_group, my_group_size, &
                                         13,"", "event", dummy)
  !nrecords info
  call define_adios_scalar (adios_group, my_group_size, "", "nreceivers",&
                        dum_int)
  call define_adios_scalar (adios_group, my_group_size, "", "nrecords",&
                        dum_int)
  !frequency(period) info
  call define_adios_scalar (adios_group, my_group_size, "", "min_period", &
                        dum_real)
  call define_adios_scalar (adios_group, my_group_size, "", "max_period", &
                        dum_real)

  !string info
  call define_adios_scalar (adios_group, my_group_size, "", "receiver_name_len", &
                        dum_int)
  call define_adios_scalar (adios_group, my_group_size, "", "network_len", &
                        dum_int)
  call define_adios_scalar (adios_group, my_group_size, "", "receiver_id_len", &
                        dum_int)
  call define_adios_scalar (adios_group, my_group_size, "", "component_len", &
                        dum_int)

  call define_adios_global_integer_1d_array (adios_group, my_group_size,&
                   nrecords, "", "npoints", int_array)
  call define_adios_global_integer_1d_array (adios_group, my_group_size,&
                   nrecords, "", "gmt_year", int_array)
  call define_adios_global_integer_1d_array (adios_group, my_group_size,&
                   nrecords, "", "gmt_day", int_array)
  call define_adios_global_integer_1d_array (adios_group, my_group_size,&
                   nrecords, "", "gmt_hour", int_array)
  call define_adios_global_integer_1d_array (adios_group, my_group_size,&
                   nrecords, "", "gmt_min", int_array)
  call define_adios_global_integer_1d_array (adios_group, my_group_size,&
                   nrecords, "", "gmt_sec", int_array)
  call define_adios_global_integer_1d_array (adios_group, my_group_size,&
                   nrecords, "", "gmt_msec", int_array)

  string_total_length = STRING_COMMON_LENGTH * nrecords_total
  call define_adios_local_string_1d_array (adios_group, my_group_size,&
                    string_total_length, "", "receiver_name", dum_string)
  call define_adios_local_string_1d_array (adios_group, my_group_size,&
                    string_total_length, "", "network", dum_string)
  call define_adios_local_string_1d_array (adios_group, my_group_size,&
                    string_total_length, "", "component", dum_string)
  call define_adios_local_string_1d_array (adios_group, my_group_size,&
                    string_total_length, "", "receiver_id", dum_string)

  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "event_lat", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "event_lo", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "event_dpt", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "receiver_lat", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "receiver_lo", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "receiver_el", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "receiver_dpt", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "begin_value", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "end_value", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "cmp_azimuth", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "cmp_incident_ang", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "sample_rate", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "scale_factor", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "ev_to_sta_AZ", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "sta_to_ev_AZ", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "great_circle_arc", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "dist", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "P_pick", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, &
                   nrecords, "", "S_pick", real_array)

  !DISPLACEMENT
  do i = 1, nrecords
    write(i_string, '(I10)' ) i+offset
    record=trim(my_asdf%receiver_name_array(i))//"."//&
           trim(my_asdf%network_array(i))//"."//&
           trim(my_asdf%component_array(i))//"."//&
           trim(my_asdf%receiver_id_array(i))
    call define_adios_global_real_1d_array (adios_group, my_group_size,&
         my_asdf%npoints(i), "", trim(record), real_array)
  enddo

  !define attribute
  call adios_define_attribute ( adios_group , "nreceivers", "desc", &
        adios_string, "Number of receivers ", "" , adios_err )
  call adios_define_attribute ( adios_group , "nrecords", "desc", &
        adios_string, "Number of records ", "" , adios_err )
  call adios_define_attribute ( adios_group , "min_period", "desc", &
        adios_string, "Low pass filter in Hz (0 if none applied)  ", "" , adios_err )
  call adios_define_attribute ( adios_group , "max_period", "desc", &
        adios_string, "High pass filter in Hz (0 if none applied)  ", "" , adios_err )
  call adios_define_attribute ( adios_group , "event_lat", "desc", adios_string, &
        "Event CMT latitude (degrees, north positive) ", "", adios_err )
  call adios_define_attribute ( adios_group , "event_lo", "desc", adios_string, &
        "Event CMT longitude (degrees, east positive) ", "", adios_err )
  call adios_define_attribute ( adios_group , "event_dpt", "desc", adios_string, &
        "Event CMT depth (km) ", "" , adios_err )
  call adios_define_attribute ( adios_group , "event_dpt", "desc", adios_string, &
        "Event CMT depth (km) ", "" , adios_err )
  call adios_define_attribute ( adios_group , "component", "desc", adios_string, &
        "Record component ", "" , adios_err)
  call adios_define_attribute ( adios_group, "gmt_year", "desc", adios_string, &
        "GMT year corresponding to reference (zero) time in file. ", "" , adios_err)
  call adios_define_attribute ( adios_group, "gmt_day", "desc", adios_string, &
        "GMT julian day corresponding to reference (zero) time in file. ", "" , adios_err)
  call adios_define_attribute ( adios_group, "gmt_hour", "desc", adios_string, &
        "GMT hour corresponding to reference (zero) time in file. ", "" , adios_err)
  call adios_define_attribute ( adios_group, "gmt_min", "desc", adios_string, &
        "GMT minute corresponding to reference (zero) time in file. ", "" , adios_err)
  call adios_define_attribute ( adios_group, "gmt_sec", "desc", adios_string, &
        "GMT second corresponding to reference (zero) time in file. ", "" , adios_err)
  call adios_define_attribute ( adios_group, "gmt_msec", "desc", adios_string, &
        "GMT millisecond corresponding to reference (zero) time in file. ", "" , adios_err)
  call adios_define_attribute ( adios_group , "receiver_lat", "desc", adios_string, &
        "Receiver latitude (degrees, north positive)  ", "" , adios_err )
  call adios_define_attribute ( adios_group , "receiver_lo", "desc", adios_string, &
        "Receiver longitude (degrees, east positive) ", "" , adios_err )
  call adios_define_attribute ( adios_group , "receiver_dpt", "desc", adios_string, &
        "Receiver depth below surface (meters) ", "" , adios_err )
  call adios_define_attribute ( adios_group , "receiver_el", "desc", adios_string, &
        "Receiver elevation (meters) ", "" , adios_err )
  call adios_define_attribute ( adios_group , "begin_value", "desc", adios_string, &
        "Beginning value of time array ", "" , adios_err )
  call adios_define_attribute ( adios_group , "end_value", "desc", adios_string, &
        "End value of time array ", "" , adios_err )
  call adios_define_attribute ( adios_group , "cmp_azimuth", "desc", adios_string, &
        "Component azimuth (degrees clockwise from north) ", "", adios_err )
  call adios_define_attribute ( adios_group , "cmp_incident_ang", "desc", adios_string,&
        "Component incident angle (degrees from vertical) ", "", adios_err )
  call adios_define_attribute ( adios_group , "sample_rate", "desc", adios_string, &
        "Sampling rate (s) ", "" , adios_err )
  call adios_define_attribute ( adios_group , "scale_factor", "desc", adios_string, &
        "Scale factor to convert the unit of synthetics from meters to nanometer ", &
        "" , adios_err )
  call adios_define_attribute ( adios_group , "ev_to_sta_AZ", "desc", adios_string, &
        "Event to station azimuth (degrees) ", "" , adios_err )
  call adios_define_attribute ( adios_group , "sta_to_ev_AZ", "desc", adios_string, &
        "Station to event azimuth (backazimuth, degrees) ", "", adios_err )
  call adios_define_attribute ( adios_group , "great_circle_dist", "desc", adios_string, &
        "Great circle distance between event and station (degrees) ", "", adios_err )
  call adios_define_attribute ( adios_group , "receiver_name", "desc", adios_string, &
        "Receiver name ", "" , adios_err )
  call adios_define_attribute( adios_group , "network", "desc", adios_string, &
        "Receiver network name ", "" , adios_err )
  call adios_define_attribute( adios_group , "receiver_id", "desc", adios_string, &
        "Receiver number ", "" , adios_err )
  call adios_define_attribute ( adios_group , "component", "desc", adios_string,&
        "Receiver component name ", "" , adios_err )
 
end subroutine define_asdf_data

!> Writes the asdf data structure to the adios arrays
!! \param my_asdf The asdf data structure
!! \param adios_handle The asdf file name
!! \param adios_group The adios group
!! \param adios_groupsize The adios group size
!! \param rank The rank of the processor
!! \param nproc The number of processors
!! \param comm The communication group of processors
!! \param ierr The error for adios subroutine calls
subroutine write_asdf_data_sub (my_asdf, adios_handle, my_adios_group,&
                        adios_groupsize, rank, nproc, comm, ierr)

  use adios_write_mod
  use asdf_data
  use asdf_helpers_writers_mod
  use mpi

  implicit none
  integer                       :: adios_err, i
  integer(kind=8),intent(in)    :: my_adios_group, adios_groupsize
  integer(kind=8),intent(in)    :: adios_handle
  integer,intent(in)            :: rank, nproc, comm, ierr
  integer :: nrecords_total, offset, nreceivers
  integer :: receiver_name_len, network_len, component_len, receiver_id_len
  integer :: rn_len_total, nw_len_total, rid_len_total, comp_len_total
  integer :: rn_offset, nw_offset, rid_offset, comp_offset
  character(len=32)              :: loc_string

  character(len=:), allocatable :: receiver_name, network, component, receiver_id
  character(len=:), allocatable :: receiver_name_total, network_total, &
                                  component_total, receiver_id_total

  type(asdf_event) :: my_asdf

  !gather array offset info
  call gather_offset_info(my_asdf%nrecords,nrecords_total,offset,&
                                        rank, nproc, comm, ierr)

  !ensemble the string for receiver_name, network, componen and receiver_id
  allocate(character(len=6*my_asdf%nrecords) :: receiver_name)
  allocate(character(len=6*my_asdf%nrecords) :: network)
  allocate(character(len=6*my_asdf%nrecords) :: component)
  allocate(character(len=6*my_asdf%nrecords) :: receiver_id)
  receiver_name=''
  network=''
  component=''
  receiver_id=''

  do i=1, my_asdf%nrecords
    receiver_name=trim(receiver_name)//trim(my_asdf%receiver_name_array(i))//'.'
    network=trim(network)//trim(my_asdf%network_array(i))//'.'
    component=trim(component)//trim(my_asdf%component_array(i))//'.'
    receiver_id=trim(receiver_id)//trim(my_asdf%receiver_id_array(i))//'.'
  enddo
  receiver_name_len = len_trim(receiver_name)
  network_len = len_trim(network)
  component_len = len_trim(component)
  receiver_id_len = len_trim(receiver_id)

  !get global dimensions for strings
  call gather_string_total_length(receiver_name_len, rn_len_total,&
                                          rank, nproc, comm, ierr)
  call gather_string_total_length(network_len, nw_len_total,&
                                          rank, nproc, comm, ierr)
  call gather_string_total_length(receiver_id_len, rid_len_total,&
                                          rank, nproc, comm, ierr)
  call gather_string_total_length(component_len, comp_len_total,&
                                          rank, nproc, comm, ierr)
  allocate(character(len=rn_len_total) :: receiver_name_total)
  allocate(character(len=nw_len_total) :: network_total)
  allocate(character(len=rid_len_total) :: receiver_id_total)
  allocate(character(len=comp_len_total) :: component_total)

  !write all local strings into global string
  call gather_string_offset_info(receiver_name_len, rn_len_total,rn_offset, &
                                        receiver_name, receiver_name_total,&
                                        rank, nproc, comm, ierr)
  call gather_string_offset_info(network_len, nw_len_total, nw_offset, &
                                        network, network_total,&
                                        rank, nproc, comm, ierr)
  call gather_string_offset_info(component_len, comp_len_total, comp_offset, &
                                        component, component_total,&
                                        rank, nproc, comm, ierr)
  call gather_string_offset_info(receiver_id_len, rid_len_total,rid_offset, &
                                        receiver_id, receiver_id_total,&
                                        rank, nproc, comm, ierr)
  !===========================
  !write out the string info
  if(rank.eq.0)then
    call adios_write(adios_handle, "receiver_name", trim(receiver_name_total),adios_err)
    call adios_write(adios_handle, "network", trim(network_total), adios_err)
    call adios_write(adios_handle, "component", trim(component_total), adios_err)
    call adios_write(adios_handle, "receiver_id", trim(receiver_id_total),adios_err)
  endif
    deallocate(receiver_name_total)
    deallocate(network_total)
    deallocate(receiver_id_total)
    deallocate(component_total)

  !===========================
  ! write seismic records
  do i = 1, my_asdf%nrecords
    write( loc_string, '(I10)' ) i+offset
    loc_string=trim(my_asdf%receiver_name_array(i))//"."//&
               trim(my_asdf%network_array(i))//"."//&
               trim(my_asdf%component_array(i))//"."//&
               trim(my_asdf%receiver_id_array(i))
     call write_adios_global_1d_array(adios_handle, rank, nproc, &
                       my_asdf%npoints(i), my_asdf%npoints(i), 0, &
                       loc_string, my_asdf%records(i)%record)
  enddo

  !===========================
  !scalar
  if(rank.eq.0)then
    call adios_write(adios_handle, "nrecords", nrecords_total, adios_err)
    call adios_write(adios_handle, "receiver_name_len", rn_len_total, adios_err)
    call adios_write(adios_handle, "network_len", nw_len_total, adios_err)
    call adios_write(adios_handle, "component_len", comp_len_total, adios_err)
    call adios_write(adios_handle, "receiver_id_len", rid_len_total, adios_err)
    call adios_write(adios_handle, "nreceivers", nreceivers, adios_err)
    call adios_write(adios_handle, "min_period", 0, adios_err)
    call adios_write(adios_handle, "max_period", 0, adios_err)
    call adios_write(adios_handle, "event", my_asdf%event, adios_err)
  endif

 !===========================
  !write out the array
  call write_adios_global_integer_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "npoints", my_asdf%npoints)
  call write_adios_global_integer_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "gmt_year", my_asdf%gmt_year)
  call write_adios_global_integer_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "gmt_day", my_asdf%gmt_day)
  call write_adios_global_integer_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "gmt_hour", my_asdf%gmt_hour)
  call write_adios_global_integer_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "gmt_min", my_asdf%gmt_min)
  call write_adios_global_integer_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "gmt_sec", my_asdf%gmt_sec)
  call write_adios_global_integer_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "gmt_msec", my_asdf%gmt_msec)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "event_lat", my_asdf%event_lat)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "event_lo", my_asdf%event_lo)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "event_dpt", my_asdf%event_dpt)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "receiver_lat", my_asdf%receiver_lat)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "receiver_lo", my_asdf%receiver_lo)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "receiver_el", my_asdf%receiver_el)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "receiver_dpt", my_asdf%receiver_dpt)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "begin_value", my_asdf%begin_value)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "end_value", my_asdf%end_value)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "cmp_azimuth", my_asdf%cmp_azimuth)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "cmp_incident_ang", my_asdf%cmp_incident_ang)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "sample_rate", my_asdf%sample_rate)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "scale_factor", my_asdf%scale_factor)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "ev_to_sta_AZ", my_asdf%ev_to_sta_AZ)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "sta_to_ev_AZ", my_asdf%sta_to_ev_AZ)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "great_circle_arc", my_asdf%great_circle_arc)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "dist", my_asdf%dist)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "P_pick", my_asdf%P_pick)
  call write_adios_global_real_1d_array(adios_handle, rank, nproc, my_asdf%nrecords,&
        nrecords_total, offset, "S_pick", my_asdf%S_pick)

  deallocate(receiver_name)
  deallocate(network)
  deallocate(receiver_id)
  deallocate(component)
end subroutine write_asdf_data_sub

!> Gets offset values for arrays
!! \param local_dim The local dimension on the processor
!! \param global_dim The global dimension of the array
!! \param The offset for the processor
!! \param rank The rank of the processor
!! \param nproc The number of processors
!! \param comm The communication group of processors
!! \param ierr The error for adios subroutine calls
subroutine gather_offset_info(local_dim, global_dim, offset,&
                                                rank, nproc, comm, ierr)

  use mpi
  implicit none

  integer :: local_dim, global_dim, offset
  integer :: rank, nproc, comm, ierr

  integer, allocatable :: local_dim_all_proc(:)
  integer, allocatable :: offset_all_proc(:)
  integer :: i

  allocate(local_dim_all_proc(nproc))
  allocate(offset_all_proc(nproc))
        
  call synchronize_all()
  call MPI_Gather(local_dim, 1, MPI_INTEGER, local_dim_all_proc, 1, &
                                        MPI_INTEGER, 0, comm, ierr)

  if(rank.eq.0)then
    offset_all_proc(1)=0
    do i=2, nproc
      offset_all_proc(i)=sum(local_dim_all_proc(1:(i-1)))
    enddo
    global_dim=sum(local_dim_all_proc(1:nproc))
  endif

  call MPI_Scatter(offset_all_proc, 1, MPI_INTEGER, offset, &
                              1, MPI_INTEGER, 0, comm, ierr)
  call MPI_Bcast(global_dim, 1, MPI_INTEGER, 0, comm, ierr)

  deallocate(local_dim_all_proc)
  deallocate(offset_all_proc)

end subroutine gather_offset_info

!> Gets total length of strings from each processor
!! \param local_dim The local dimension on the processor
!! \param global_dim The global dimension of the array
!! \param rank The rank of the processor
!! \param nproc The number of processors
!! \param comm The communication group of processors
!! \param ierr The error for adios subroutine calls
subroutine gather_string_total_length(local_dim, global_dim,&
                                          rank, nproc, comm, ierr)

  use mpi
  implicit none

  integer :: local_dim, global_dim
  integer :: rank, nproc, comm, ierr

  integer, allocatable :: local_dim_all_proc(:)
  integer :: i, tag, mpi_status(MPI_STATUS_SIZE)

  if(rank.eq.0)then
    allocate(local_dim_all_proc(nproc))
  endif

  call synchronize_all()
  call MPI_Gather(local_dim, 1, MPI_INTEGER, local_dim_all_proc, 1, &
                                        MPI_INTEGER, 0, comm, ierr)
  call synchronize_all()
  if(rank.eq.0)then
    global_dim=sum(local_dim_all_proc(1:nproc))
    deallocate(local_dim_all_proc)
  endif

end subroutine gather_string_total_length

!> Gets offset values for strings
!! \param local_dim The local dimension on the processor
!! \param global_dim The global dimension of the array
!! \param offset The offset for the string
!! \param string_piece The local string
!! \param string_total The combined string from all processors
!! \param rank The rank of the processor
!! \param nproc The number of processors
!! \param comm The communication group of processors
!! \param ierr The error
subroutine gather_string_offset_info(local_dim, global_dim, offset,&
                                         string_piece, string_total,&
                                                rank, nproc, comm, ierr)

  use mpi
  implicit none

  integer :: local_dim, global_dim, offset
  character(len=*) :: string_piece, string_total
  character(len=10000) :: buffer_string
  integer :: rank, nproc, comm, ierr

  integer, allocatable :: local_dim_all_proc(:)
  integer, allocatable :: offset_all_proc(:)
  integer :: i, tag, mpi_status(MPI_STATUS_SIZE)

  if(rank.eq.0)then
    allocate(local_dim_all_proc(nproc))
    allocate(offset_all_proc(nproc))
  endif

  call synchronize_all()
  call MPI_Gather(local_dim, 1, MPI_INTEGER, local_dim_all_proc, 1, &
                                        MPI_INTEGER, 0, comm, ierr)

  call synchronize_all()
  if(rank.eq.0)then
    offset_all_proc(1)=0
    do i=2, nproc
      offset_all_proc(i)=sum(local_dim_all_proc(1:(i-1)))
    enddo
    !print *, "offset_all_proc:", offset_all_proc(:)
    string_total=''
    buffer_string=''
    string_total=trim(string_total)//trim(string_piece(1:local_dim))
  endif
        
  !print *,"TAG1"
  !if(rank.eq.0) then
  !  print *,"global_dim",global_dim
  !endif

  call synchronize_all()
  if(rank.eq.0)then
    offset_all_proc(1)=0
    do i=2, nproc
      offset_all_proc(i)=sum(local_dim_all_proc(1:(i-1)))
    enddo
    !print *, "offset_all_proc:", offset_all_proc(:)
    string_total=''
    buffer_string=''
    string_total=trim(string_total)//trim(string_piece(1:local_dim))
  endif
        
  !print *,"TAG1"
  !if(rank.eq.0) then
  !print *,"global_dim",global_dim
  !endif

  !call synchronize_all()
  if(rank.eq.0)then
    do i=1,nproc-1
      !print *, "buffer_before:",trim(buffer_string)
      !print *, "local_dim_all_proc:",local_dim_all_proc(i+1)
      call MPI_Recv(buffer_string, local_dim_all_proc(i+1),MPI_CHARACTER,&
                                                        i, 1, comm, mpi_status,ierr)
      !print *,"buffer_string:", trim(buffer_string)
      string_total=trim(string_total)//buffer_string(1:local_dim_all_proc(i+1))
    enddo
  else
   !print *, "local_dim:", local_dim
   !print *,"string_piece:", trim(string_piece)
   call MPI_Send(string_piece, local_dim, MPI_CHARACTER,&
                                                        0, 1, comm, ierr)
  endif
  !print *,"TAG", rank

  call synchronize_all()
  call MPI_Scatter(offset_all_proc, 1, MPI_INTEGER, offset, &
                                       1, MPI_INTEGER, 0, comm, ierr)
  call MPI_Bcast(global_dim, 1, MPI_INTEGER, 0, comm, ierr)

  !print *,"rank, local dim, global_dim,offset:", rank, local_dim, &
  !                                               global_dim, offset
  if (rank.eq.0) then
    deallocate(local_dim_all_proc)
    deallocate(offset_all_proc)
  endif

end subroutine gather_string_offset_info
