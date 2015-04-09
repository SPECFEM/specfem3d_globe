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

  subroutine read_topography_bathymetry()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! local parameters
  integer :: ier
  ! timing
  double precision :: tCPU
  double precision, external :: wtime
  character(len=MAX_STRING_LEN) :: topo_path

  ! get MPI starting time
  time_start = wtime()

  ! make ellipticity
  if (ELLIPTICITY_VAL) then
    ! splines used for locating exact source/receivers positions
    ! in locate_sources() and locate_receivers() routines
    call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)
  endif

  ! read topography and bathymetry file
  if (TOPOGRAPHY) then
    ! allocates topography array
    allocate(ibathy_topo(NX_BATHY,NY_BATHY),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating ibathy_topo array')

    ! initializes
    ibathy_topo(:,:) = 0

    if (I_should_read_the_database) then
      ! master reads file
      if (myrank == 0) then
        ! user output
        write(IMAIN,*) 'topography:'
        call flush_IMAIN()

        topo_path = LOCAL_PATH

        ! reads topo file
        call read_topo_bathy_database(ibathy_topo, topo_path)
      endif

      ! broadcast the information read on the master to the nodes
      call bcast_all_i(ibathy_topo,NX_BATHY*NY_BATHY)
    endif
    call bcast_ibathy_topo()
  endif

  ! user output
  call synchronize_all()
  if (myrank == 0 .and. (TOPOGRAPHY .or. OCEANS_VAL .or. ELLIPTICITY_VAL)) then
    ! elapsed time since beginning of mesh generation
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for reading topo/bathy in seconds = ',sngl(tCPU)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine read_topography_bathymetry

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_ibathy_topo()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  call bcast_all_i_for_database(ibathy_topo(1,1), size(ibathy_topo))
  end subroutine bcast_ibathy_topo
