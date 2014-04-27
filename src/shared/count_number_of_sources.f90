!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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

  subroutine count_number_of_sources(NSOURCES)

! count the total number of sources in the CMTSOLUTION file
! there are NLINES_PER_CMTSOLUTION_SOURCE lines per source in that file

  use constants

  implicit none

  integer, intent(out) :: NSOURCES

  integer :: ios,icounter

  character(len=150) :: CMTSOLUTION,dummystring

  call get_value_string(CMTSOLUTION, 'solver.CMTSOLUTION', 'DATA/CMTSOLUTION')

  open(unit=IIN,file=trim(CMTSOLUTION),status='old',action='read',iostat=ios)
  if(ios /= 0) stop 'error opening CMTSOLUTION file'

  icounter = 0
  do while(ios == 0)
    read(IIN,"(a)",iostat=ios) dummystring
    if(ios == 0) icounter = icounter + 1
  enddo
  close(IIN)

  if(mod(icounter,NLINES_PER_CMTSOLUTION_SOURCE) /= 0) &
    stop 'total number of lines in CMTSOLUTION file should be a multiple of NLINES_PER_CMTSOLUTION_SOURCE'

  NSOURCES = icounter / NLINES_PER_CMTSOLUTION_SOURCE

  if(NSOURCES < 1) stop 'need at least one source in CMTSOLUTION file'

  end subroutine count_number_of_sources

