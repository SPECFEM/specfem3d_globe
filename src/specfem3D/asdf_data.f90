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

!> Module defining the data structure for ASDF
!!     * Waveforms defined per event
!!     * Waveform attributes defined per seismogram per event
!! \author JAS and Wenjie

module asdf_data

  type asdf_record
    real, pointer :: record(:)
  end type asdf_record

  type asdf_event
    !scalars
    character(len=13)     :: event
    real, pointer     :: event_lat(:), event_lo(:), event_dpt(:)

    !processing info
    real              :: min_period, max_period

    !size info
    integer           :: nrecords
    integer           :: nreceivers

    !time info
    integer, pointer    :: gmt_year(:), gmt_day(:), gmt_hour(:)
    integer, pointer    :: gmt_min(:), gmt_sec(:), gmt_msec(:)

    !seismic record info
    integer, pointer    :: npoints(:)
    real, pointer       :: receiver_lat(:), receiver_lo(:)
    real, pointer       :: receiver_el(:),  receiver_dpt(:)
    real, pointer       :: begin_value(:),  end_value(:)
    real, pointer       :: cmp_azimuth(:),  cmp_incident_ang(:)
    real, pointer       :: sample_rate(:),  scale_factor(:)

    real, pointer       :: ev_to_sta_AZ(:), sta_to_ev_AZ(:)
    real, pointer       :: great_circle_arc(:)
    real, pointer       :: dist(:)
    real, pointer       :: P_pick(:), S_pick(:)

    character(len=20),pointer :: receiver_name_array(:), network_array(:)
    character(len=20),pointer :: component_array(:)
    character(len=20),pointer :: receiver_id_array(:)

    !seismograms
    type (asdf_record), pointer :: records(:)

  end type asdf_event

end module asdf_data
