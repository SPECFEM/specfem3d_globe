!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
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

module gindex_par

  use constants, only: myrank,NGLLX,NGLLY,NGLLZ

  use shared_parameters, only: NPROCTOT

  use specfem_par

  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use specfem_par_trinfinite
  use specfem_par_infinite

  use specfem_par_full_gravity

  !use specfem_par_movie
  implicit none

  ! global nodes for NGLLX = 5
  integer :: ignode_end
  ! global gdof for NGLLX = 5
  integer :: gnf_end
  ! global gdof for NGLLX_INF = 3
  integer :: gnf_end1

end module gindex_par


