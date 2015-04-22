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

!==============================================================================
! Copyright 2015 ASDF developers
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!==============================================================================

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

end subroutine write_asdf
