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

!-------------------------------------------------------------------------------
!> \file write_output_ASDF.F90
!! \brief Write subroutines for writing ASDF seismograms to file using
!!        the ADIOS library
!! \author JAS and Wenjie Lei
!------------------------------------------------------------------------------

#include "config.fh"

!> Initializes the data structure for ASDF
!! \param asdf_container The ASDF data structure
!! \param total_seismos_local The number of records on the local processor
subroutine init_asdf_data(asdf_container,total_seismos_local)

  use asdf_data
  use specfem_par, only : event_name_SAC,myrank

  ! Parameters
  type(asdf_event),intent(inout) :: asdf_container
  integer,intent(in) :: total_seismos_local

  ! Variables
  integer :: ier

  asdf_container%nrecords = total_seismos_local
  asdf_container%event = trim(event_name_SAC)

  allocate (asdf_container%npoints(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%gmt_year(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%gmt_hour(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%gmt_day(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%gmt_min(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%gmt_sec(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%gmt_msec(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%event_lat(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%event_lo(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%event_dpt(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%receiver_lat(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%receiver_lo(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%receiver_el(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%receiver_dpt(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%begin_value(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%end_value(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%cmp_azimuth(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%cmp_incident_ang(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%sample_rate(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%scale_factor(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%ev_to_sta_AZ(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%sta_to_ev_AZ(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%great_circle_arc(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%dist(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%P_pick(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%S_pick(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%records(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%receiver_name_array(asdf_container%nrecords), &
            STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%network_array(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%component_array(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')
  allocate (asdf_container%receiver_id_array(asdf_container%nrecords), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocate failed.')

end subroutine init_asdf_data

!
!-------------------------------------------------------------------------------------------------
!

!> Stores the records into the ASDF structure
!! \param asdf_container The ASDF data structure
!! \param seismogram_tmp The current seismogram to store
!! \param irec_local The local index of the receivers on the local processor
!! \param irec The global index of the receiver
!! \param chn The broadband channel simulated
!! \param iorientation The recorded seismogram's orientation direction
!! \param phi The angle used for calculating azimuth and incident angle
subroutine store_asdf_data(asdf_container, seismogram_tmp, irec_local, &
                           irec, chn, iorientation, phi)

  use asdf_data
  use specfem_par,only: &
    station_name,network_name,stlat,stlon,stele,stbur, &
    DT,t0, seismo_offset,seismo_current, NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
    yr=>yr_SAC,jda=>jda_SAC,ho=>ho_SAC,mi=>mi_SAC,sec=>sec_SAC, &
    tshift_cmt=>t_cmt_SAC, &
    cmt_lat=>cmt_lat_SAC,cmt_lon=>cmt_lon_SAC, &
    cmt_depth=>cmt_depth_SAC

  use specfem_par, only: myrank
  use constants

  implicit none

  ! Parameters
  type(asdf_event),intent(inout) :: asdf_container
  character(len=4),intent(in) :: chn
  integer,intent(in) :: irec_local, irec
  real(kind=CUSTOM_REAL),dimension(5,NTSTEP_BETWEEN_OUTPUT_SEISMOS),intent(in) :: seismogram_tmp
  integer,intent(in) :: iorientation
  double precision,intent(in) :: phi
  ! local Variables
  integer :: length_station_name, length_network_name
  integer :: ier, i

  ! trace index
  i = (irec_local-1)*(3) + (iorientation)

  asdf_container%npoints(i) = seismo_current
  asdf_container%gmt_year(i) = yr
  asdf_container%gmt_day(i) = jda
  asdf_container%gmt_hour(i) = ho
  asdf_container%gmt_min(i) = mi
  asdf_container%gmt_sec(i) = sec
  asdf_container%gmt_msec(i) = 0
  asdf_container%event_lat(i) = cmt_lat
  asdf_container%event_lo(i) = cmt_lon
  asdf_container%event_dpt(i) = cmt_depth
  asdf_container%receiver_lat(i) = stlat(irec_local)
  asdf_container%receiver_lo(i) = stlon(irec_local)
  asdf_container%receiver_el(i) = stele(irec_local)
  asdf_container%receiver_dpt(i) = stbur(irec_local)
  asdf_container%begin_value(i) = seismo_offset*DT-t0+tshift_cmt
  asdf_container%end_value(i) = -12345

  ! instrument orientation
  if (iorientation == 1) then !N
    asdf_container%cmp_azimuth(i)  = 0.00
    asdf_container%cmp_incident_ang(i) =90.00
  else if (iorientation == 2) then !E
    asdf_container%cmp_azimuth(i)  =90.00
    asdf_container%cmp_incident_ang(i) =90.00
  else if (iorientation == 3) then !Z
    asdf_container%cmp_azimuth(i)  = 0.00
    asdf_container%cmp_incident_ang(i) = 0.00
  else if (iorientation == 4) then !R
    asdf_container%cmp_azimuth(i) = sngl(modulo(phi,360.0d+0))
    asdf_container%cmp_incident_ang(i) =90.00
  else if (iorientation == 5) then !T
    asdf_container%cmp_azimuth(i) = sngl(modulo(phi+90.0,360.0d+0))
    asdf_container%cmp_incident_ang(i) =90.00
  endif

  asdf_container%sample_rate(i) = DT
  asdf_container%scale_factor(i) = 1000000000
  asdf_container%ev_to_sta_AZ(i) = -12345
  asdf_container%sta_to_ev_AZ(i) = -12345
  asdf_container%great_circle_arc(i) = -12345
  asdf_container%dist(i) = -12345
  asdf_container%P_pick(i) = -12345
  asdf_container%S_pick(i) = -12345

  length_station_name = len_trim(station_name(irec))
  length_network_name = len_trim(network_name(irec))
  asdf_container%receiver_name_array(i) = station_name(irec)(1:length_station_name)
  asdf_container%network_array(i) = network_name(irec)(1:length_network_name)
  asdf_container%component_array(i) = chn

  ! note: this array of strings is not set with other new values yet...
  asdf_container%receiver_id_array(i) = ""

  allocate (asdf_container%records(i)%record(seismo_current), STAT=ier)
  if (ier /= 0) call exit_MPI (myrank, 'Allocating ASDF container failed.')

  asdf_container%records(i)%record(1:seismo_current) = seismogram_tmp(iorientation, 1:seismo_current)

end subroutine store_asdf_data

!
!-------------------------------------------------------------------------------------------------
!

!> Closes the ASDF data structure by deallocating all arrays
!! \param asdf_container The ASDF data structure
!! \param total_seismos_local The number of seismograms on the local processor
subroutine close_asdf_data(asdf_container, total_seismos_local)

  use asdf_data
  ! Parameters
  type(asdf_event),intent(inout) :: asdf_container
  integer,intent(in) :: total_seismos_local
  !Variables
  integer :: i

  deallocate (asdf_container%npoints)
  deallocate (asdf_container%gmt_year)
  deallocate (asdf_container%gmt_hour)
  deallocate (asdf_container%gmt_day)
  deallocate (asdf_container%gmt_min)
  deallocate (asdf_container%gmt_sec)
  deallocate (asdf_container%gmt_msec)
  deallocate (asdf_container%event_lat)
  deallocate (asdf_container%event_lo)
  deallocate (asdf_container%event_dpt)
  deallocate (asdf_container%receiver_lat)
  deallocate (asdf_container%receiver_lo)
  deallocate (asdf_container%receiver_el)
  deallocate (asdf_container%receiver_dpt)
  deallocate (asdf_container%begin_value)
  deallocate (asdf_container%end_value)
  deallocate (asdf_container%cmp_azimuth)
  deallocate (asdf_container%cmp_incident_ang)
  deallocate (asdf_container%sample_rate)
  deallocate (asdf_container%scale_factor)
  deallocate (asdf_container%ev_to_sta_AZ)
  deallocate (asdf_container%sta_to_ev_AZ)
  deallocate (asdf_container%great_circle_arc)
  deallocate (asdf_container%dist)
  deallocate (asdf_container%P_pick)
  deallocate (asdf_container%S_pick)
  do i = 1, total_seismos_local
    deallocate(asdf_container%records(i)%record)
  enddo
  deallocate (asdf_container%receiver_name_array)
  deallocate (asdf_container%network_array)
  deallocate (asdf_container%component_array)
  deallocate (asdf_container%receiver_id_array)

end subroutine close_asdf_data

!
!-------------------------------------------------------------------------------------------------
!

!> Writes the ASDF data structure to the file
!! \param asdf_container The ASDF data structure
subroutine write_asdf(asdf_container)

  use adios_write_mod,only: adios_declare_group,adios_select_method

  use asdf_data,only: asdf_event

  use specfem_par, only : event_name_SAC,myrank,ADIOS_TRANSPORT_METHOD, OUTPUT_FILES

  implicit none
  ! Parameters
  type(asdf_event),intent(inout) :: asdf_container
  ! Variables
  integer :: adios_err, comm, sizeprocs
  integer(kind=8) :: adios_group
  character(len=200) :: ASDF_FN

  ! alias mpi communicator for ADIOS
  call world_duplicate(comm)
  call world_size(sizeprocs)

  ! declare new group that uses MPI
  call adios_declare_group (adios_group, "EVENTS", "iter", 1, adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  call adios_select_method (adios_group, ADIOS_TRANSPORT_METHOD, "", "", adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  ! output file name
  ASDF_FN=trim(OUTPUT_FILES)//"/"//trim(event_name_SAC)//"_sem.bp"

  call write_asdf_data (ASDF_FN, asdf_container, adios_group, myrank, sizeprocs, comm)

end subroutine write_asdf

!
!-------------------------------------------------------------------------------------------------
!

!> Writes the ASDF data structure to asdf_fn using parallel write
!! \param asdf_fn The file name for ASDF
!! \param asdf_container The ASDF data structure
!! \param adios_group The adios group for the file
!! \param rank The rank of the processor
!! \param nproc The number of processors
!! \param comm The communication group of processors
subroutine write_asdf_data(asdf_fn, asdf_container, adios_group, rank, nproc, comm)

  use adios_write_mod, only: adios_open,adios_group_size

  use adios_helpers_mod,only: check_adios_err

  use asdf_data,only: asdf_event

  implicit none
  ! Parameters
  character(len=*),intent(inout) :: asdf_fn
  type(asdf_event),intent(inout) :: asdf_container
  integer(kind=8),intent(inout) :: adios_group
  integer,intent(inout) :: rank, nproc, comm
  ! Variables
  integer         :: adios_err
  integer(kind=8) :: adios_groupsize, adios_totalsize
  integer(kind=8) :: adios_handle

  !calculate size
  adios_groupsize = 0
  call define_asdf_data(adios_group, adios_groupsize, asdf_container,rank, nproc)

  ! Open the handle to file containing all the ADIOS variables
  call adios_open(adios_handle, "EVENTS", asdf_fn, "w", comm, adios_err)
  if (adios_err /= 0) then
    print*,'Error: rank ',rank,' could not open adios file ',trim(asdf_fn)
    stop 'Error calling adios_open() routine failed for EVENTS'
  endif

  call adios_group_size(adios_handle, adios_groupsize, adios_totalsize,adios_err)
  if (adios_err /= 0 ) stop 'Error calling adios_group_size() routine failed'

  !call the write sub
  call write_asdf_data_sub (asdf_container, adios_handle, rank, nproc)

  !adios close
  call adios_close(adios_handle, adios_err)
  call check_adios_err(rank,adios_err)

end subroutine write_asdf_data

!
!-------------------------------------------------------------------------------------------------
!

!> Defines the ASDF structure using adios
!! \param adios_group The adios group
!! \param my_group_size The adios group size
!! \param asdf_container The ASDF data structure
!! \param rank The rank of the processor
!! \param nproc The number of processors
subroutine define_asdf_data (adios_group, my_group_size, asdf_container, rank, nproc)

  use adios_write_mod,only: adios_string

  use adios_helpers_mod,only: define_adios_scalar,define_adios_global_integer_1d_array, &
    define_adios_global_real_1d_array,define_adios_local_string_1d_array !,check_adios_err

  use asdf_data,only: asdf_event

  implicit none

  ! Parameters
  integer(kind=8), intent(inout) :: adios_group, my_group_size
  type(asdf_event), intent(inout) :: asdf_container
  integer, intent(in) :: rank, nproc

  ! Variables
  integer :: i,string_total_length
  integer, parameter :: STRING_COMMON_LENGTH = 20
  integer :: adios_err

  integer :: nrecords

  character(len=80)            :: str_record
  character(len=10)            :: i_string
  character(len=200)           :: dummy

  integer :: dum_int, int_array(10)
  real    :: dum_real, real_array(10)
  character(len=10) :: dum_string

  integer :: nrecords_total, offset

  integer,parameter :: nparam_desc = 32
  character(len=200),dimension(2,nparam_desc) :: description

  !gather info. Here, we only need nrecords_total
  nrecords=asdf_container%nrecords
  call gather_offset_info(nrecords,nrecords_total,offset,rank,nproc)

  call define_adios_local_string_1d_array (adios_group, my_group_size, &
                                           13,"", "event", dummy)
  !nrecords info
  call define_adios_scalar (adios_group, my_group_size, "", "nreceivers", dum_int)
  call define_adios_scalar (adios_group, my_group_size, "", "nrecords", dum_int)
  !frequency(period) info
  call define_adios_scalar (adios_group, my_group_size, "", "min_period", dum_real)
  call define_adios_scalar (adios_group, my_group_size, "", "max_period", dum_real)

  !string info
  call define_adios_scalar (adios_group, my_group_size, "", "receiver_name_len", dum_int)
  call define_adios_scalar (adios_group, my_group_size, "", "network_len", dum_int)
  call define_adios_scalar (adios_group, my_group_size, "", "receiver_id_len", dum_int)
  call define_adios_scalar (adios_group, my_group_size, "", "component_len", dum_int)

  call define_adios_global_integer_1d_array (adios_group, my_group_size, nrecords, "", "npoints", int_array)
  call define_adios_global_integer_1d_array (adios_group, my_group_size, nrecords, "", "gmt_year", int_array)
  call define_adios_global_integer_1d_array (adios_group, my_group_size, nrecords, "", "gmt_day", int_array)
  call define_adios_global_integer_1d_array (adios_group, my_group_size, nrecords, "", "gmt_hour", int_array)
  call define_adios_global_integer_1d_array (adios_group, my_group_size, nrecords, "", "gmt_min", int_array)
  call define_adios_global_integer_1d_array (adios_group, my_group_size, nrecords, "", "gmt_sec", int_array)
  call define_adios_global_integer_1d_array (adios_group, my_group_size, nrecords, "", "gmt_msec", int_array)

  string_total_length = STRING_COMMON_LENGTH * nrecords_total

  call define_adios_local_string_1d_array (adios_group, my_group_size, string_total_length, "", "receiver_name", dum_string)
  call define_adios_local_string_1d_array (adios_group, my_group_size, string_total_length, "", "network", dum_string)
  call define_adios_local_string_1d_array (adios_group, my_group_size, string_total_length, "", "component", dum_string)
  call define_adios_local_string_1d_array (adios_group, my_group_size, string_total_length, "", "receiver_id", dum_string)

  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "event_lat", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "event_lo", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "event_dpt", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "receiver_lat", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "receiver_lo", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "receiver_el", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "receiver_dpt", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "begin_value", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "end_value", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "cmp_azimuth", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "cmp_incident_ang", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "sample_rate", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "scale_factor", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "ev_to_sta_AZ", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "sta_to_ev_AZ", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "great_circle_arc", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "dist", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "P_pick", real_array)
  call define_adios_global_real_1d_array (adios_group, my_group_size, nrecords, "", "S_pick", real_array)

  !DISPLACEMENT
  do i = 1, nrecords
    write(i_string, '(I10)' ) i+offset
    str_record = trim(asdf_container%receiver_name_array(i))//"."// &
                 trim(asdf_container%network_array(i))//"."// &
                 trim(asdf_container%component_array(i))//"."// &
                 trim(asdf_container%receiver_id_array(i))

    call define_adios_global_real_1d_array (adios_group, my_group_size,&
                                            asdf_container%npoints(i), "", trim(str_record), real_array)
  enddo

  ! defines attributes descriptions

  DATA description  / &
    "nreceivers", "Number of receivers ", &
    "nrecords"  , "Number of records ", &
    "min_period", "Low pass filter in Hz (0 if none applied)  ", &
    "max_period", "High pass filter in Hz (0 if none applied)  ", &
    "event_lat" , "Event CMT latitude (degrees, north positive) ", &
    "event_lo"  , "Event CMT longitude (degrees, east positive) ", &
    "event_dpt" , "Event CMT depth (km) ", &
    "event_dpt" , "Event CMT depth (km) ", &
    "component" , "Record component ", &
    "gmt_year"  , "GMT year corresponding to reference (zero) time in file. ", &
    "gmt_day"   , "GMT julian day corresponding to reference (zero) time in file. ", &
    "gmt_hour"  , "GMT hour corresponding to reference (zero) time in file. ", &
    "gmt_min"   , "GMT minute corresponding to reference (zero) time in file. ", &
    "gmt_sec"   , "GMT second corresponding to reference (zero) time in file. ", &
    "gmt_msec"  , "GMT millisecond corresponding to reference (zero) time in file. ", &
    "receiver_lat", "Receiver latitude (degrees, north positive)  ", &
    "receiver_lo" , "Receiver longitude (degrees, east positive) ", &
    "receiver_dpt", "Receiver depth below surface (meters) ", &
    "receiver_el" , "Receiver elevation (meters) ", &
    "begin_value" , "Beginning value of time array ", &
    "end_value"   , "End value of time array ", &
    "cmp_azimuth" , "Component azimuth (degrees clockwise from north) ", &
    "cmp_incident_ang", "Component incident angle (degrees from vertical) ", &
    "sample_rate"     , "Sampling rate (s) ", &
    "scale_factor"    , "Scale factor to convert the unit of synthetics from meters to nanometer ", &
    "ev_to_sta_AZ"    , "Event to station azimuth (degrees) ", &
    "sta_to_ev_AZ"    , "Station to event azimuth (backazimuth, degrees) ", &
    "great_circle_dist" , "Great circle distance between event and station (degrees) ", &
    "receiver_name"     , "Receiver name ", &
    "network"           , "Receiver network name ", &
    "receiver_id"       , "Receiver number ", &
    "component"         , "Receiver component name " &
  /

  do i = 1,nparam_desc
    call adios_define_attribute(adios_group,trim(description(1,i)),"desc",adios_string,trim(description(2,i)),"",adios_err)
    ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
    !call check_adios_err(rank,adios_err)
  enddo

end subroutine define_asdf_data

!
!-------------------------------------------------------------------------------------------------
!

!> Writes the ASDF data structure to the adios arrays
!! \param asdf_container The ASDF data structure
!! \param adios_handle The ASDF file name
!! \param adios_group The adios group
!! \param adios_groupsize The adios group size
!! \param rank The rank of the processor
!! \param nproc The number of processors
subroutine write_asdf_data_sub(asdf_container, adios_handle, rank, nproc)

  use adios_write_mod,only:adios_string,adios_write

  use adios_helpers_mod,only: write_adios_global_1d_array_offset, &
    write_adios_global_integer_1d_array_offset,write_adios_global_real_1d_array_offset

  use asdf_data,only: asdf_event

  implicit none

  type(asdf_event),intent(inout) :: asdf_container
  integer(kind=8),intent(in)    :: adios_handle
  integer,intent(in)            :: rank, nproc

  ! local parameters
  integer :: adios_err,i !,ierr
  integer :: nrecords_total, offset, nreceivers
  integer :: receiver_name_len, network_len, component_len, receiver_id_len
  integer :: rn_len_total, nw_len_total, rid_len_total, comp_len_total
  integer :: rn_offset, nw_offset, rid_offset, comp_offset
  character(len=32) :: loc_string

  ! note: this is fortran 2003 standard
  !       and works e.g. by intel ifort compilers and newer gfortran versions;
  !
  !       gfortran 4.7 complains about the allocate statement later on (a gfortran bug which needs a constant at compile time)
  !       which is needed on Titan, thus we provide an alternative way of a character array
  ! way 1
!  character(len=:), allocatable :: receiver_name, network, component, receiver_id
!  character(len=:), allocatable :: receiver_name_total, network_total, &
!                                  component_total, receiver_id_total
  !
  ! way 2: fortran 95 workaround, static buffers
  integer,parameter :: BUFFER_LENGTH = 100000
  integer,parameter :: BUFFER_LENGTH_TOTAL = 600000
  character(len=BUFFER_LENGTH) :: receiver_name, network,component, receiver_id
  character(len=BUFFER_LENGTH_TOTAL) :: receiver_name_total,network_total,component_total,receiver_id_total

  !gather array offset info
  call gather_offset_info(asdf_container%nrecords,nrecords_total,offset,rank, nproc)

  !ensemble the string for receiver_name, network, component and receiver_id
  ! way 1: fortran 2003
!  allocate(character(len=6*asdf_container%nrecords) :: receiver_name, STAT=ierr)
!  if (ierr /= 0) call exit_MPI (rank, 'Allocate failed.')
!  allocate(character(len=6*asdf_container%nrecords) :: network, STAT=ierr)
!  if (ierr /= 0) call exit_MPI (rank, 'Allocate failed.')
!  allocate(character(len=6*asdf_container%nrecords) :: component, STAT=ierr)
!  if (ierr /= 0) call exit_MPI (rank, 'Allocate failed.')
!  allocate(character(len=6*asdf_container%nrecords) :: receiver_id, STAT=ierr)
!  if (ierr /= 0) call exit_MPI (rank, 'Allocate failed.')
  ! way 2: fortran 95
  if (6*asdf_container%nrecords > BUFFER_LENGTH) then
    print*,'Error: buffer length too small - minimum length is ',6*asdf_container%nrecords
    stop 'Error in write_asdf_data_sub() routine, BUFFER_LENGTH too small'
  endif

  ! initializes strings
  receiver_name=''
  network=''
  component=''
  receiver_id=''

  ! appends all strings
  do i = 1, asdf_container%nrecords
    receiver_name = trim(receiver_name) // trim(asdf_container%receiver_name_array(i))  // '.'
    network       = trim(network)       // trim(asdf_container%network_array(i))        // '.'
    component     = trim(component)     // trim(asdf_container%component_array(i))      // '.'
    receiver_id   = trim(receiver_id)   // trim(asdf_container%receiver_id_array(i))    // '.'
  enddo

  ! local string lengths
  receiver_name_len = len_trim(receiver_name)
  network_len = len_trim(network)
  component_len = len_trim(component)
  receiver_id_len = len_trim(receiver_id)

  ! synchronize processes
  call synchronize_all()

  !get global dimensions for strings
  call gather_string_total_length(receiver_name_len, rn_len_total,rank, nproc)
  call gather_string_total_length(network_len, nw_len_total,rank, nproc)
  call gather_string_total_length(receiver_id_len, rid_len_total,rank, nproc)
  call gather_string_total_length(component_len, comp_len_total,rank, nproc)

  if (rank == 0) then
    ! way 1: fortran 2003
!    allocate(character(len=rn_len_total) :: receiver_name_total, STAT=ierr)
!    if (ierr /= 0) call exit_MPI (rank, 'Allocate failed.')
!    allocate(character(len=nw_len_total) :: network_total, STAT=ierr)
!    if (ierr /= 0) call exit_MPI (rank, 'Allocate failed.')
!    allocate(character(len=rid_len_total) :: receiver_id_total, STAT=ierr)
!    if (ierr /= 0) call exit_MPI (rank, 'Allocate failed.')
!    allocate(character(len=comp_len_total) :: component_total, STAT=ierr)
!    if (ierr /= 0) call exit_MPI (rank, 'Allocate failed.')

    ! way 2: fortran 95
    if (rn_len_total > BUFFER_LENGTH_TOTAL .or. nw_len_total > BUFFER_LENGTH_TOTAL .or. &
        rid_len_total > BUFFER_LENGTH_TOTAL .or. comp_len_total > BUFFER_LENGTH_TOTAL) then
      print*,'Error: buffer length total too small - lengths are ',rn_len_total,nw_len_total,rid_len_total,comp_len_total
      stop 'Error in write_asdf_data_sub() routine, BUFFER_LENGTH_TOTAL too small'
    endif
  else
    ! dummy allocation
    ! way 1: fortran 2003
!    allocate( character(len=1) :: receiver_name_total )
!    allocate( character(len=1) :: network_total )
!    allocate( character(len=1) :: receiver_id_total )
!    allocate( character(len=1) :: component_total )
  endif
  call synchronize_all()

  !write all local strings into global string
  call gather_string_offset_info(receiver_name_len, rn_len_total,rn_offset, &
                                 receiver_name, receiver_name_total, &
                                 rank, nproc)
  call gather_string_offset_info(network_len, nw_len_total, nw_offset, &
                                 network, network_total, &
                                 rank, nproc)
  call gather_string_offset_info(component_len, comp_len_total, comp_offset, &
                                 component, component_total, &
                                 rank, nproc)
  call gather_string_offset_info(receiver_id_len, rid_len_total,rid_offset, &
                                 receiver_id, receiver_id_total, &
                                 rank, nproc)
  !==========================
  !write out the string info
  if (rank == 0) then
    call adios_write(adios_handle, "receiver_name", trim(receiver_name_total),adios_err)
    call adios_write(adios_handle, "network", trim(network_total), adios_err)
    call adios_write(adios_handle, "component",trim(component_total), adios_err)
    call adios_write(adios_handle, "receiver_id", trim(receiver_id_total), adios_err)
  endif

  ! way 1: fortran 2003
!  deallocate(receiver_name_total)
!  deallocate(network_total)
!  deallocate(receiver_id_total)
!  deallocate(component_total)

  !===========================
  ! write seismic records
  do i = 1, asdf_container%nrecords
    write( loc_string, '(I10)' ) i+offset
    loc_string=trim(asdf_container%receiver_name_array(i))  // "." // &
               trim(asdf_container%network_array(i))        // "." // &
               trim(asdf_container%component_array(i))      // "." // &
               trim(asdf_container%receiver_id_array(i))

    call write_adios_global_1d_array_offset(adios_handle, rank, nproc, &
                                            asdf_container%npoints(i), asdf_container%npoints(i), 0, &
                                            loc_string, asdf_container%records(i)%record)
  enddo

  !===========================
  !scalar
  if (rank == 0) then
    call adios_write(adios_handle, "nrecords", nrecords_total, adios_err)
    call adios_write(adios_handle, "receiver_name_len", rn_len_total, adios_err)
    call adios_write(adios_handle, "network_len", nw_len_total, adios_err)
    call adios_write(adios_handle, "component_len", comp_len_total, adios_err)
    call adios_write(adios_handle, "receiver_id_len", rid_len_total, adios_err)
    call adios_write(adios_handle, "nreceivers", nreceivers, adios_err)
    call adios_write(adios_handle, "min_period", 0, adios_err)
    call adios_write(adios_handle, "max_period", 0, adios_err)
    call adios_write(adios_handle, "event", asdf_container%event, adios_err)
  endif

  !===========================
  !write out the array values
  ! integer values
  call write_adios_global_integer_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                                  STRINGIFY_VAR_TYPE(asdf_container,npoints) )
  call write_adios_global_integer_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                                  STRINGIFY_VAR_TYPE(asdf_container,gmt_year) )
  call write_adios_global_integer_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                                  STRINGIFY_VAR_TYPE(asdf_container,gmt_day) )
  call write_adios_global_integer_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                                  STRINGIFY_VAR_TYPE(asdf_container,gmt_hour) )
  call write_adios_global_integer_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                                  STRINGIFY_VAR_TYPE(asdf_container,gmt_min) )
  call write_adios_global_integer_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                                  STRINGIFY_VAR_TYPE(asdf_container,gmt_sec) )
  call write_adios_global_integer_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                                  STRINGIFY_VAR_TYPE(asdf_container,gmt_msec) )

  ! real values
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,event_lat) )
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,event_lo) )
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,event_dpt) )
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,receiver_lat) )
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,receiver_lo) )
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,receiver_el) )
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,receiver_dpt) )
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,begin_value) )
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,end_value) )
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,cmp_azimuth) )
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,cmp_incident_ang) )
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,sample_rate) )
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,scale_factor) )
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,ev_to_sta_AZ) )
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,sta_to_ev_AZ) )
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,great_circle_arc) )
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,dist) )
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,P_pick) )
  call write_adios_global_real_1d_array_offset(adios_handle, rank, nproc,asdf_container%nrecords,nrecords_total,offset, &
                                               STRINGIFY_VAR_TYPE(asdf_container,S_pick) )

  ! way 1: fortran 2003
!  deallocate(receiver_name)
!  deallocate(network)
!  deallocate(receiver_id)
!  deallocate(component)

end subroutine write_asdf_data_sub

!
!-------------------------------------------------------------------------------------------------
!

!> Gets offset values for arrays
!! \param local_dim The local dimension on the processor
!! \param global_dim The global dimension of the array
!! \param The offset for the processor
!! \param rank The rank of the processor
!! \param nproc The number of processors
subroutine gather_offset_info(local_dim, global_dim, offset, rank, nproc)

  implicit none

  integer,intent(inout) :: local_dim, global_dim, offset
  integer,intent(in) :: rank, nproc

  ! local parameters
  integer, allocatable :: dim_all_proc(:)
  integer, allocatable :: offset_proc(:)
  integer :: i,ierr

  if (rank == 0) then
    allocate(dim_all_proc(nproc), STAT=ierr)
    if (ierr /= 0) call exit_MPI (rank, 'Allocate failed.')
    allocate(offset_proc(nproc), STAT=ierr)
    if (ierr /= 0) call exit_MPI (rank, 'Allocate failed.')
  else
    ! dummy allocation
    allocate(dim_all_proc(1))
    allocate(offset_proc(1))
  endif

  call gather_all_singlei(local_dim,dim_all_proc,nproc)

  if (rank == 0) then
    offset_proc(1) = 0
    do i=2, nproc
      offset_proc(i)=sum(dim_all_proc(1:(i-1)))
    enddo
    global_dim=sum(dim_all_proc(1:nproc))
  endif

  call scatter_all_singlei(offset_proc,offset,nproc)

  call bcast_all_singlei(global_dim)

  deallocate(dim_all_proc)
  deallocate(offset_proc)

end subroutine gather_offset_info

!
!-------------------------------------------------------------------------------------------------
!

!> Gets total length of strings from each processor
!! \param local_dim The local dimension on the processor
!! \param global_dim The global dimension of the array
!! \param rank The rank of the processor
!! \param nproc The number of processors
subroutine gather_string_total_length(local_dim, global_dim, rank, nproc)

  implicit none

  integer,intent(inout) :: local_dim, global_dim
  integer,intent(in) :: rank, nproc

  ! local parameters
  integer, allocatable :: local_dim_all_proc(:)
  integer :: ierr

  if (rank == 0) then
    allocate(local_dim_all_proc(nproc),STAT=ierr)
    if (ierr /= 0) call exit_MPI (rank, 'Allocate failed.')
  else
    ! dummy allocation
    allocate(local_dim_all_proc(1))
  endif

  call gather_all_singlei(local_dim,local_dim_all_proc,nproc)

  if (rank == 0) then
    global_dim=sum(local_dim_all_proc(1:nproc))
  endif

  deallocate(local_dim_all_proc)

end subroutine gather_string_total_length

!
!-------------------------------------------------------------------------------------------------
!

!> Gets offset values for strings
!! \param local_dim The local dimension on the processor
!! \param global_dim The global dimension of the array
!! \param offset The offset for the string
!! \param string_piece The local string
!! \param string_total The combined string from all processors
!! \param rank The rank of the processor
!! \param nproc The number of processors
subroutine gather_string_offset_info(local_dim, global_dim, offset, &
                                     string_piece, string_total, &
                                     rank, nproc)

  use constants,only: itag

  implicit none

  integer,intent(inout) :: local_dim, global_dim, offset
  character(len=*),intent(in) :: string_piece
  character(len=*),intent(inout) :: string_total

  integer,intent(in) :: rank, nproc

  ! local parameters
  integer,parameter :: BUFFER_LENGTH = 100000
  character(len=BUFFER_LENGTH) :: buffer_string

  integer,dimension(:),allocatable :: local_dim_all_proc,offset_all_proc
  integer :: i,ierr

  ! checks local string
  if (len_trim(string_piece) /= local_dim ) stop 'Error local string and local dim have different lengths'

  ! temporary arrays
  if (rank == 0) then
    allocate(local_dim_all_proc(nproc),STAT=ierr)
    if (ierr /= 0) call exit_MPI (rank, 'Allocate failed.')
    allocate(offset_all_proc(nproc),STAT=ierr)
    if (ierr /= 0) call exit_MPI (rank, 'Allocate failed.')
  else
    ! dummy allocation
    allocate(local_dim_all_proc(1))
    allocate(offset_all_proc(1))
  endif
  call synchronize_all()

  ! master gets all local string lengths
  call gather_all_singlei(local_dim,local_dim_all_proc,nproc)

  if (rank == 0) then
    offset_all_proc(1) = 0
    do i=2, nproc
      offset_all_proc(i)=sum(local_dim_all_proc(1:(i-1)))
    enddo
    ! adds strings from master process 0
    string_total=''
    string_total=trim(string_total)//trim(string_piece(1:local_dim))
  endif

  if (rank == 0) then
    do i = 1,nproc-1
      ! checks if buffer length is sufficient
      if (local_dim_all_proc(i+1) > BUFFER_LENGTH) &
        stop 'Error send/recv buffer length too small in gather_string_offset_info() routine'

      ! receives string
      buffer_string=''
      call recv_ch(buffer_string,local_dim_all_proc(i+1),i,itag)

      ! appends string
      string_total = trim(string_total) // buffer_string(1:local_dim_all_proc(i+1))
    enddo
  else
    ! sends string
    call send_ch(string_piece, local_dim, 0, itag)
  endif

  ! master sends offset values to corresponding processes
  call scatter_all_singlei(offset_all_proc,offset,nproc)

  ! broadcasts global length
  call bcast_all_singlei(global_dim)

  ! frees temporary arrays
  deallocate(local_dim_all_proc)
  deallocate(offset_all_proc)

end subroutine gather_string_offset_info
