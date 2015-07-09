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


  subroutine setup_inner_outer(iregion_code)

  use meshfem3D_par,only: &
    myrank,OUTPUT_FILES,IMAIN, &
    IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE,MAX_STRING_LEN

  use meshfem3D_par,only: ibool,is_on_a_slice_edge

  use MPI_crust_mantle_par
  use MPI_outer_core_par
  use MPI_inner_core_par

  implicit none

  integer,intent(in) :: iregion_code

  ! local parameters
  real :: percentage_edge
  integer :: ier,ispec,iinner,iouter
  ! debug file output
  character(len=MAX_STRING_LEN) :: filename
  logical,parameter :: DEBUG = .false.

  ! stores inner / outer elements
  !
  ! note: arrays is_on_a_slice_edge_.. have flags set for elements which need to
  !         communicate with other MPI processes
  select case (iregion_code)
  case (IREGION_CRUST_MANTLE)
    ! crust_mantle
    nspec_outer_crust_mantle = count( is_on_a_slice_edge )
    nspec_inner_crust_mantle = NSPEC_CRUST_MANTLE - nspec_outer_crust_mantle

    num_phase_ispec_crust_mantle = max(nspec_inner_crust_mantle,nspec_outer_crust_mantle)

    allocate(phase_ispec_inner_crust_mantle(num_phase_ispec_crust_mantle,2),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_crust_mantle')

    phase_ispec_inner_crust_mantle(:,:) = 0
    iinner = 0
    iouter = 0
    do ispec = 1,NSPEC_CRUST_MANTLE
      if (is_on_a_slice_edge(ispec)) then
        ! outer element
        iouter = iouter + 1
        phase_ispec_inner_crust_mantle(iouter,1) = ispec
      else
        ! inner element
        iinner = iinner + 1
        phase_ispec_inner_crust_mantle(iinner,2) = ispec
      endif
    enddo

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'for overlapping of communications with calculations:'
      write(IMAIN,*)
      percentage_edge = 100. * nspec_outer_crust_mantle / real(NSPEC_CRUST_MANTLE)
      write(IMAIN,*) 'percentage of edge elements in crust/mantle ',percentage_edge,'%'
      write(IMAIN,*) 'percentage of volume elements in crust/mantle ',100. - percentage_edge,'%'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! debug: saves element flags
    if (DEBUG) then
      write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_innerouter_crust_mantle_proc',myrank
      call write_VTK_data_elem_l(NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
                                xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                                ibool, &
                                is_on_a_slice_edge,filename)
    endif

  case (IREGION_OUTER_CORE)
    ! outer_core
    nspec_outer_outer_core = count( is_on_a_slice_edge )
    nspec_inner_outer_core = NSPEC_OUTER_CORE - nspec_outer_outer_core

    num_phase_ispec_outer_core = max(nspec_inner_outer_core,nspec_outer_outer_core)

    allocate(phase_ispec_inner_outer_core(num_phase_ispec_outer_core,2),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_outer_core')

    phase_ispec_inner_outer_core(:,:) = 0
    iinner = 0
    iouter = 0
    do ispec = 1,NSPEC_OUTER_CORE
      if (is_on_a_slice_edge(ispec)) then
        ! outer element
        iouter = iouter + 1
        phase_ispec_inner_outer_core(iouter,1) = ispec
      else
        ! inner element
        iinner = iinner + 1
        phase_ispec_inner_outer_core(iinner,2) = ispec
      endif
    enddo

    ! user output
    if (myrank == 0) then
      percentage_edge = 100.* nspec_outer_outer_core / real(NSPEC_OUTER_CORE)
      write(IMAIN,*) 'percentage of edge elements in outer core ',percentage_edge,'%'
      write(IMAIN,*) 'percentage of volume elements in outer core ',100. - percentage_edge,'%'
      write(IMAIN,*)
    endif

    ! debug: saves element flags
    if (DEBUG) then
      write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_innerouter_outer_core_proc',myrank
      call write_VTK_data_elem_l(NSPEC_OUTER_CORE,NGLOB_OUTER_CORE, &
                                xstore_outer_core,ystore_outer_core,zstore_outer_core, &
                                ibool, &
                                is_on_a_slice_edge,filename)
    endif

  case (IREGION_INNER_CORE)
    ! inner_core
    nspec_outer_inner_core = count( is_on_a_slice_edge )
    nspec_inner_inner_core = NSPEC_INNER_CORE - nspec_outer_inner_core

    num_phase_ispec_inner_core = max(nspec_inner_inner_core,nspec_outer_inner_core)

    allocate(phase_ispec_inner_inner_core(num_phase_ispec_inner_core,2),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_inner_core')

    phase_ispec_inner_inner_core(:,:) = 0
    iinner = 0
    iouter = 0
    do ispec = 1,NSPEC_INNER_CORE
      if (is_on_a_slice_edge(ispec)) then
        ! outer element
        iouter = iouter + 1
        phase_ispec_inner_inner_core(iouter,1) = ispec
      else
        ! inner element
        iinner = iinner + 1
        phase_ispec_inner_inner_core(iinner,2) = ispec
      endif
    enddo

    ! user output
    if (myrank == 0) then
      percentage_edge = 100. * nspec_outer_inner_core / real(NSPEC_INNER_CORE)
      write(IMAIN,*) 'percentage of edge elements in inner core ',percentage_edge,'%'
      write(IMAIN,*) 'percentage of volume elements in inner core ',100. - percentage_edge,'%'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! debug: saves element flags
    if (DEBUG) then
      write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_innerouter_inner_core_proc',myrank
      call write_VTK_data_elem_l(NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
                                xstore_inner_core,ystore_inner_core,zstore_inner_core, &
                                ibool, &
                                is_on_a_slice_edge,filename)
    endif

  end select

  end subroutine setup_Inner_Outer
