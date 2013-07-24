!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
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

  use mpi
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! local parameters
  integer :: ier
  ! timing
  double precision, external :: wtime

  ! get MPI starting time
  time_start = wtime()

  ! make ellipticity
  if( ELLIPTICITY_VAL ) then
    ! splines used for locating exact source/receivers positions
    ! in locate_sources() and locate_receivers() routines
    call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)
  endif

  ! read topography and bathymetry file
  if( TOPOGRAPHY ) then
    ! allocates topography array
    allocate(ibathy_topo(NX_BATHY,NY_BATHY),stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating ibathy_topo array')

    ! initializes
    ibathy_topo(:,:) = 0

    ! master reads file
    if(myrank == 0 ) then
      ! user output
      write(IMAIN,*) 'topography:'

      ! reads topo file
      call read_topo_bathy_database(ibathy_topo,LOCAL_PATH)
    endif

    ! broadcast the information read on the master to the nodes
    call MPI_BCAST(ibathy_topo,NX_BATHY*NY_BATHY,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  endif

  ! user output
  call sync_all()
  if( myrank == 0 .and. (TOPOGRAPHY .or. OCEANS_VAL .or. ELLIPTICITY_VAL)) then
    ! elapsed time since beginning of mesh generation
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for reading topo/bathy in seconds = ',sngl(tCPU)
    write(IMAIN,*)
  endif

  end subroutine read_topography_bathymetry
