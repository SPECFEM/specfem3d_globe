!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                 and University of Pau / CNRS, France
! (c) California Institute of Technology and University of Pau / CNRS, November 2007
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

  program trampoline

! Perform Fortran mojo, and then run the Python script.

! With ifort v9 in particular, this function (i.e., MAIN__) will call
! the undocumented function __intel_new_proc_init or
! __intel_new_proc_init_P.  Without this, SPECFEM runs several
! times slower (!).

  call FC_PY_MAIN()

  end program trampoline
