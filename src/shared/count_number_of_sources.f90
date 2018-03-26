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

  subroutine count_number_of_sources(NSOURCES)

! count the total number of sources in the CMTSOLUTION file
! there are NLINES_PER_CMTSOLUTION_SOURCE lines per source in that file

  use constants
  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,USE_FORCE_POINT_SOURCE

  implicit none

  integer, intent(out) :: NSOURCES

  integer :: ios,icounter,nline

  character(len=MAX_STRING_LEN) :: SOURCE_FILE, path_to_add
  character(len=MAX_STRING_LEN) :: dummystring

  if (USE_FORCE_POINT_SOURCE) then
    ! FORCESOLUTION file
    SOURCE_FILE = 'DATA/FORCESOLUTION'
    nline = NLINES_PER_FORCESOLUTION_SOURCE
  else
    ! CMTSOLUTION file
    SOURCE_FILE = 'DATA/CMTSOLUTION'
    nline = NLINES_PER_CMTSOLUTION_SOURCE
  endif

  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    SOURCE_FILE = path_to_add(1:len_trim(path_to_add))//SOURCE_FILE(1:len_trim(SOURCE_FILE))
  endif

  open(unit=IIN,file=trim(SOURCE_FILE),status='old',action='read',iostat=ios)
  if (ios /= 0) then
    print *,'Error opening source SOLUTION file: ',trim(SOURCE_FILE)
    stop 'Error opening source SOLUTION file'
  endif

  icounter = 0
  do while(ios == 0)
    read(IIN,"(a)",iostat=ios) dummystring
    if (ios == 0) icounter = icounter + 1
  enddo
  close(IIN)

  if (mod(icounter,nline) /= 0) then
    if (USE_FORCE_POINT_SOURCE) then
      stop 'total number of lines in FORCESOLUTION file should be a multiple of NLINES_PER_FORCESOLUTION_SOURCE'
    else
      stop 'total number of lines in CMTSOLUTION file should be a multiple of NLINES_PER_CMTSOLUTION_SOURCE'
    endif
  endif

  NSOURCES = icounter / nline

  if (NSOURCES < 1) then
    print *,'Error: ',trim(SOURCE_FILE),' has ',icounter,'lines but need ',nline,'per source... ',NSOURCES
    stop 'need at least one source in CMTSOLUTION or FORCESOLUTION file'
  endif

  end subroutine count_number_of_sources

