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

  subroutine create_name_database(prname,iproc,iregion_code,LOCAL_PATH)

! create the name of the database for the mesher and the solver

  use constants,only: MAX_STRING_LEN

  implicit none

  integer iproc,iregion_code

! name of the database file
  character(len=MAX_STRING_LEN) prname,procname,LOCAL_PATH

! create the name for the database of the current slide and region
  write(procname,"('/proc',i6.6,'_reg',i1,'_')") iproc,iregion_code

! create full name with path
  prname = trim(LOCAL_PATH) // procname

  end subroutine create_name_database

subroutine create_name_database_adios(prname,iregion_code,LOCAL_PATH)

  ! create the name of the database for the mesher and the solver

  use constants,only: MAX_STRING_LEN

  implicit none

  integer iregion_code

! name of the database file
  character(len=MAX_STRING_LEN) prname,procname,LOCAL_PATH

! create the name for the database of the current slide and region
  write(procname,"('/reg',i1,'_')") iregion_code

! create full name with path
  prname = trim(LOCAL_PATH) // procname

end subroutine create_name_database_adios
