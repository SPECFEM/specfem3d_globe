!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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

  subroutine read_topography_bathymetry()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  implicit none

  include 'mpif.h'

  ! local parameters
  integer :: ier
  
  ! get MPI starting time
  time_start = MPI_WTIME()

  ! make ellipticity
  if(ELLIPTICITY_VAL) then
    call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)
  endif

  ! read topography and bathymetry file
  if(myrank == 0 .and. (TOPOGRAPHY .or. OCEANS_VAL)) then
    call read_topo_bathy_file(ibathy_topo)
  endif

  ! broadcast the information read on the master to the nodes
  call MPI_BCAST(ibathy_topo,NX_BATHY*NY_BATHY,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  ! user output
  if( myrank == 0 .and. (TOPOGRAPHY .or. OCEANS_VAL)) then
    ! elapsed time since beginning of mesh generation
    tCPU = MPI_WTIME() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for reading topo/bathy in seconds = ',sngl(tCPU)
    write(IMAIN,*)
  endif

  end subroutine read_topography_bathymetry
