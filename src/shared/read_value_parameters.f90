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

! read values from parameter file, ignoring white lines and comments

  subroutine read_value_integer(value_to_read, name, ier)

  implicit none

  integer value_to_read
  character(len=*) name
  character(len=100) string_read
  integer ier

  call param_read(string_read, len(string_read), name, len(name), ier)
  if (ier /= 0) return
  read(string_read,*,iostat=ier) value_to_read

  end subroutine read_value_integer

!--------------------

  subroutine read_value_double_precision(value_to_read, name, ier)

  implicit none

  double precision value_to_read
  character(len=*) name
  character(len=100) string_read
  integer ier

  call param_read(string_read, len(string_read), name, len(name), ier)
  if (ier /= 0) return
  read(string_read,*,iostat=ier) value_to_read

  end subroutine read_value_double_precision

!--------------------

  subroutine read_value_logical(value_to_read, name, ier)

  implicit none

  logical value_to_read
  character(len=*) name
  character(len=100) string_read
  integer ier

  call param_read(string_read, len(string_read), name, len(name), ier)
  if (ier /= 0) return
  read(string_read,*,iostat=ier) value_to_read

  end subroutine read_value_logical

!--------------------

  subroutine read_value_string(value_to_read, name, ier)

  implicit none

  character(len=*) value_to_read
  character(len=*) name
  character(len=256) string_read
  integer ier

  call param_read(string_read, len(string_read), name, len(name), ier)
  if (ier /= 0) return
  value_to_read = string_read

  end subroutine read_value_string

!--------------------

  subroutine open_parameter_file_from_master_only(ier)

  implicit none

  integer :: ier
  character(len=50) :: filename_main,filename_run0001
  logical :: exists_main_Par_file,exists_run0001_Par_file

  filename_main = 'DATA/Par_file'

! also see if we are running several independent runs in parallel
! to do so, add the right directory for that run for the master process only here
  filename_run0001 = 'run0001/'//filename_main(1:len_trim(filename_main))

  call param_open(filename_main, len(filename_main), ier)
  if (ier == 0) then
    exists_main_Par_file = .true.
    call close_parameter_file()
  else
    exists_main_Par_file    = .false.
  endif

  call param_open(filename_run0001, len(filename_run0001), ier)
  if (ier == 0) then
    exists_run0001_Par_file = .true.
    call close_parameter_file()
  else
    exists_run0001_Par_file = .false.
  endif

  if (exists_main_Par_file .and. exists_run0001_Par_file) then
    print *
    print *,'cannot have both DATA/Par_file and run0001/DATA/Par_file present, please remove one of them'
    stop 'error: two different copies of the Par_file'
  endif

  call param_open(filename_main, len(filename_main), ier)
  if (ier /= 0) then
    call param_open(filename_run0001, len(filename_run0001), ier)
    if (ier /= 0) then
      print *
      print *,'opening file failed, please check your file path and run-directory.'
      stop 'error opening Par_file'
    endif
  endif

  end subroutine open_parameter_file_from_master_only

!--------------------

  subroutine open_parameter_file(ier)

  integer ier
  character(len=50) filename
  filename = 'DATA/Par_file'

  call param_open(filename, len(filename), ier)
  if (ier /= 0) return

  end subroutine open_parameter_file

!--------------------

  subroutine close_parameter_file

  call param_close()

  end subroutine close_parameter_file
