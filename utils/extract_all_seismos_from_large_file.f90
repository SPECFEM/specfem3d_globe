!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            August 2013
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


  program extract_all_seismos_large_file

! extract all the seismograms when they are stored in a unique large file

! Dimitri Komatitsch, University of Pau, France, November 2007

  implicit none

  include "../constants.h"

! binary or ASCII storage format (see its value in DATA/Par_file instead)
  logical, parameter :: USE_BINARY_FOR_LARGE_FILE = .true.

! number of seismogram files stored in the unique large file
  integer, parameter :: N_STATIONS = 647 ! see number of lines in DATA/STATIONS
  integer, parameter :: N_COMPONENTS = 3
  integer, parameter :: NREC = N_STATIONS * N_COMPONENTS

! number of time steps in each seismogram file
  integer, parameter :: NSTEP = 37200 ! see its value in OUTPUT_FILES/values_from_mesher.h

  integer :: irec,istep,irepeat
  real :: time,U_value

  character(len=150) :: station_name

! open the large seismogram file
  if(USE_BINARY_FOR_LARGE_FILE) then
    open(unit=30,file='all_seismograms.bin',status='old',form='unformatted',action='read')
  else
    open(unit=30,file='all_seismograms.ascii',status='old',action='read')
  endif

! loop on all the seismogram files
  do irec = 1,NREC

    if(USE_BINARY_FOR_LARGE_FILE) then
      read(30) station_name
    else
      read(30,*) station_name
    endif

! suppress leading white spaces, if any
    station_name = adjustl(station_name)

! suppress two leading '\' (ASCII code 92), if any
    do irepeat = 1,2
      if(station_name(1:1) == achar(92)) station_name = station_name(2:len_trim(station_name))
      station_name = adjustl(station_name)
    enddo

    print *,'extracting seismogram file ',irec,': ',station_name(1:len_trim(station_name)),' out of ',NREC

    open(unit=27,file=station_name(1:len_trim(station_name)),status='unknown')

! loop on all the time steps in each seismogram
    do istep = 1,NSTEP
      if(USE_BINARY_FOR_LARGE_FILE) then
        read(30) time, U_value
      else
        read(30,*) time, U_value
      endif
      write(27,*) time, U_value
    enddo

    close(27)

  enddo

  close(30)

  end program extract_all_seismos_large_file

