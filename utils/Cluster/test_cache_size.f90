!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                    and University of Pau, France
! (c) California Institute of Technology / Princeton University and University of Pau, October 2007
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

  program test_cache_size

! test processor cache size as in the book of Deville et al. (2002), Chapter 8

! Dimitri Komatitsch, University of Pau, France, May 2007

  implicit none

! timer to count elapsed time
  character(len=8) datein
  character(len=10) timein
  character(len=5)  :: zone
  integer, dimension(8) :: time_values_in,time_values_out
  double precision :: time_start,time_end

  integer, parameter :: NELEM_ARRAY_MAX = 10000000

  real, dimension(NELEM_ARRAY_MAX) :: a

  real :: c

  double precision :: megaflops

  integer :: NELEM_ARRAY,l,chiffre,puissance,NLOOP

!***********************************************************************

! count elapsed wall-clock time
  datein = ''
  timein = ''
  zone = ''

  NELEM_ARRAY = 1000

  c = 1.0

  do puissance = 0,6

  do chiffre = 1,9

  NELEM_ARRAY = chiffre * 10**puissance

  a(:) = 0.0

  call c_add(a,c,NELEM_ARRAY)

  if(NELEM_ARRAY <= 1000) then
    NLOOP = 10000000
  else if(NELEM_ARRAY <= 10000) then
    NLOOP = 600000
  else if(NELEM_ARRAY <= 100000) then
    NLOOP = 60000
  else
    NLOOP = 5000
  endif

  call date_and_time(datein,timein,zone,time_values_in)
  do l = 1,NLOOP
    call c_add(a,c,NELEM_ARRAY)
  enddo
  call date_and_time(datein,timein,zone,time_values_out)

! time_values_in(3): day of the month
! time_values_in(5): hour of the day
! time_values_in(6): minutes of the hour
! time_values_in(7): seconds of the minute
! time_values_in(8): milliseconds of the second

! this fails if we cross the end of the hour
  time_start = 60000.d0*time_values_in(6) + 1000.d0*time_values_in(7) + time_values_in(8)
  time_end   = 60000.d0*time_values_out(6) + 1000.d0*time_values_out(7) + time_values_out(8)

  megaflops = dble(NELEM_ARRAY) * dble(NLOOP) / (time_end - time_start) / 1.d3

  print *,NELEM_ARRAY,' elements -> megaflops = ',megaflops

  enddo

  enddo

  end program test_cache_size

!---

  subroutine c_add(a,c,NELEM_ARRAY)

  implicit none

  integer :: NELEM_ARRAY

  integer :: i

  real, dimension(NELEM_ARRAY) :: a

  real :: c

  do i = 1,NELEM_ARRAY
   a(i) = a(i) + c
  enddo

  end subroutine c_add

