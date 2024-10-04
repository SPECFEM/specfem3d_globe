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


  subroutine setup_inner_outer(iregion_code)

  use meshfem_par, only: &
    myrank,OUTPUT_FILES,IMAIN,MAX_STRING_LEN, &
    IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE, &
    IREGION_TRINFINITE,IREGION_INFINITE, &
    NPROCTOT, &
    IFLAG_IN_FICTITIOUS_CUBE

  use meshfem_par, only: ibool,is_on_a_slice_edge,xstore_glob,ystore_glob,zstore_glob, &
    idoubling

  use MPI_crust_mantle_par
  use MPI_outer_core_par
  use MPI_inner_core_par

  use MPI_trinfinite_par
  use MPI_infinite_par

  implicit none

  integer,intent(in) :: iregion_code

  ! local parameters
  real :: percentage_edge
  integer :: ier,ispec,iinner,iouter
  ! debug file output
  character(len=MAX_STRING_LEN) :: filename
  logical,parameter :: DEBUG = .false.

  ! explicitly exclude ficitious inner core elements from phase_ispec_* array
  logical,parameter :: DO_EXCLUDE_FICTITIOUS_ELEMENTS = .false.

  ! stores inner / outer elements
  !
  ! note: arrays is_on_a_slice_edge_.. have flags set for elements which need to
  !         communicate with other MPI processes

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'for overlapping of communications with calculations:'
    call flush_IMAIN()
  endif

  select case (iregion_code)
  case (IREGION_CRUST_MANTLE)
    ! crust_mantle
    if (NPROCTOT > 1) then
      nspec_outer_crust_mantle = count( is_on_a_slice_edge )
    else
      nspec_outer_crust_mantle = 0
    endif
    nspec_inner_crust_mantle = NSPEC_CRUST_MANTLE - nspec_outer_crust_mantle

    num_phase_ispec_crust_mantle = max(nspec_inner_crust_mantle,nspec_outer_crust_mantle)

    allocate(phase_ispec_inner_crust_mantle(num_phase_ispec_crust_mantle,2),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_crust_mantle')

    phase_ispec_inner_crust_mantle(:,:) = 0
    iinner = 0
    iouter = 0
    do ispec = 1,NSPEC_CRUST_MANTLE
      if (is_on_a_slice_edge(ispec) .and. NPROCTOT > 1) then
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
      percentage_edge = 100. * nspec_outer_crust_mantle / real(NSPEC_CRUST_MANTLE)
      write(IMAIN,*) '  percentage of edge elements in crust/mantle ',percentage_edge,'%'
      write(IMAIN,*) '  percentage of volume elements in crust/mantle ',100. - percentage_edge,'%'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! debug: saves element flags
    if (DEBUG) then
      write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_innerouter_crust_mantle_proc',myrank
      call write_VTK_data_elem_l(NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
                                xstore_glob,ystore_glob,zstore_glob, &
                                ibool,is_on_a_slice_edge,filename)
    endif

    !debug
    !do ispec = 0,NPROCTOT-1
    !  if (myrank == ispec) then
    !    print *,myrank,'phase_ispec_inner 1: ',phase_ispec_inner_crust_mantle(1:5,1)
    !    print *,myrank,'phase_ispec_inner 2: ',phase_ispec_inner_crust_mantle(1:5,2)
    !  endif
    !  call synchronize_all()
    !enddo

  case (IREGION_OUTER_CORE)
    ! outer_core
    if (NPROCTOT > 1) then
      nspec_outer_outer_core = count( is_on_a_slice_edge )
    else
      nspec_outer_outer_core = 0
    endif
    nspec_inner_outer_core = NSPEC_OUTER_CORE - nspec_outer_outer_core

    num_phase_ispec_outer_core = max(nspec_inner_outer_core,nspec_outer_outer_core)

    allocate(phase_ispec_inner_outer_core(num_phase_ispec_outer_core,2),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_outer_core')

    phase_ispec_inner_outer_core(:,:) = 0
    iinner = 0
    iouter = 0
    do ispec = 1,NSPEC_OUTER_CORE
      if (is_on_a_slice_edge(ispec) .and. NPROCTOT > 1) then
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
      write(IMAIN,*) '  percentage of edge elements in outer core ',percentage_edge,'%'
      write(IMAIN,*) '  percentage of volume elements in outer core ',100. - percentage_edge,'%'
      write(IMAIN,*)
    endif

    ! debug: saves element flags
    if (DEBUG) then
      write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_innerouter_outer_core_proc',myrank
      call write_VTK_data_elem_l(NSPEC_OUTER_CORE,NGLOB_OUTER_CORE, &
                                xstore_glob,ystore_glob,zstore_glob, &
                                ibool,is_on_a_slice_edge,filename)
    endif

  case (IREGION_INNER_CORE)
    ! inner_core
    if (NPROCTOT > 1) then
      nspec_outer_inner_core = count( is_on_a_slice_edge )
    else
      nspec_outer_inner_core = 0
    endif
    nspec_inner_inner_core = NSPEC_INNER_CORE - nspec_outer_inner_core

    ! note: for fictitious elements in the inner core, is_on_a_slice_edge(ispec) is set to .false.
    !       in routine create_regions_elements().
    !       thus, counting the number of elements on a slice edge will only count "active" inner core elements.
    !       and there is no need to explicitly exclude fictitious elements from this nspec_outer_inner_core count again.
    !
    !       however, at the moment fictitious elements are still included in the nspec_inner_inner_core count.
    !       we could further exclude those such that when we loop over phase_ispec_inner_inner_core(*,*) elements
    !       only "active" elements get considered.
    !       this would lead to a total count (nspec_inner_inner_core + nspec_outer_inner_core) < NSPEC_INNER_CORE
    !
    ! excludes fictitious elements from count
    if (DO_EXCLUDE_FICTITIOUS_ELEMENTS) then
      ! user output
      if (myrank == 0) then
        write(IMAIN,*) '  excluding fictitious elements from inner/outer elements'
        write(IMAIN,*)
        call flush_IMAIN()
      endif

      do ispec = 1,NSPEC_INNER_CORE
        if (is_on_a_slice_edge(ispec) .and. NPROCTOT > 1) then
          ! subtract fictitious element
          if (idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) nspec_outer_inner_core = nspec_outer_inner_core - 1
        else
          ! subtract fictitious element
          if (idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) nspec_inner_inner_core = nspec_inner_inner_core - 1
        endif
      enddo
    endif

    num_phase_ispec_inner_core = max(nspec_inner_inner_core,nspec_outer_inner_core)

    allocate(phase_ispec_inner_inner_core(num_phase_ispec_inner_core,2),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_inner_core')

    phase_ispec_inner_inner_core(:,:) = 0
    iinner = 0
    iouter = 0

    if (DO_EXCLUDE_FICTITIOUS_ELEMENTS) then
      ! counts only "active" inner core elements and excludes ficititous elements
      do ispec = 1,NSPEC_INNER_CORE
        if (is_on_a_slice_edge(ispec) .and. NPROCTOT > 1) then
          ! outer element
          if (idoubling(ispec) /= IFLAG_IN_FICTITIOUS_CUBE) then
            iouter = iouter + 1
            phase_ispec_inner_inner_core(iouter,1) = ispec
          endif
        else
          ! inner element
          if (idoubling(ispec) /= IFLAG_IN_FICTITIOUS_CUBE) then
            iinner = iinner + 1
            phase_ispec_inner_inner_core(iinner,2) = ispec
          endif
        endif
      enddo
    else
      ! default
      do ispec = 1,NSPEC_INNER_CORE
        if (is_on_a_slice_edge(ispec) .and. NPROCTOT > 1) then
          ! outer element
          iouter = iouter + 1
          phase_ispec_inner_inner_core(iouter,1) = ispec
        else
          ! inner element
          iinner = iinner + 1
          phase_ispec_inner_inner_core(iinner,2) = ispec
        endif
      enddo
    endif

    ! user output
    if (myrank == 0) then
      percentage_edge = 100. * nspec_outer_inner_core / real(NSPEC_INNER_CORE)
      write(IMAIN,*) '  percentage of edge elements in inner core ',percentage_edge,'%'
      write(IMAIN,*) '  percentage of volume elements in inner core ',100. - percentage_edge,'%'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! debug: saves element flags
    if (DEBUG) then
      write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_innerouter_inner_core_proc',myrank
      call write_VTK_data_elem_l(NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
                                xstore_glob,ystore_glob,zstore_glob, &
                                ibool,is_on_a_slice_edge,filename)
    endif

  case (IREGION_TRINFINITE)
    ! transition infinite region
    if (NPROCTOT > 1) then
      nspec_outer_trinfinite = count( is_on_a_slice_edge )
    else
      nspec_outer_trinfinite = 0
    endif
    nspec_inner_trinfinite = NSPEC_TRINFINITE - nspec_outer_trinfinite

    num_phase_ispec_trinfinite = max(nspec_inner_trinfinite,nspec_outer_trinfinite)

    allocate(phase_ispec_inner_trinfinite(num_phase_ispec_trinfinite,2),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_trinfinite')

    phase_ispec_inner_trinfinite(:,:) = 0
    iinner = 0
    iouter = 0
    do ispec = 1,NSPEC_TRINFINITE
      if (is_on_a_slice_edge(ispec) .and. NPROCTOT > 1) then
        ! outer element
        iouter = iouter + 1
        phase_ispec_inner_trinfinite(iouter,1) = ispec
      else
        ! inner element
        iinner = iinner + 1
        phase_ispec_inner_trinfinite(iinner,2) = ispec
      endif
    enddo

    ! user output
    if (myrank == 0) then
      percentage_edge = 100. * nspec_outer_trinfinite / real(NSPEC_TRINFINITE)
      write(IMAIN,*) '  percentage of edge elements in transition infinite region ',percentage_edge,'%'
      write(IMAIN,*) '  percentage of volume elements in transition infinite region ',100. - percentage_edge,'%'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! debug: saves element flags
    if (DEBUG) then
      write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_innerouter_trinfinite_proc',myrank
      call write_VTK_data_elem_l(NSPEC_TRINFINITE,NGLOB_TRINFINITE, &
                                xstore_glob,ystore_glob,zstore_glob, &
                                ibool,is_on_a_slice_edge,filename)
    endif

  case (IREGION_INFINITE)
    ! infinite region
    if (NPROCTOT > 1) then
      nspec_outer_infinite = count( is_on_a_slice_edge )
    else
      nspec_outer_infinite = 0
    endif
    nspec_inner_infinite = NSPEC_INFINITE - nspec_outer_infinite

    num_phase_ispec_infinite = max(nspec_inner_infinite,nspec_outer_infinite)

    allocate(phase_ispec_inner_infinite(num_phase_ispec_infinite,2),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_infinite')

    phase_ispec_inner_infinite(:,:) = 0
    iinner = 0
    iouter = 0
    do ispec = 1,NSPEC_INFINITE
      if (is_on_a_slice_edge(ispec) .and. NPROCTOT > 1) then
        ! outer element
        iouter = iouter + 1
        phase_ispec_inner_infinite(iouter,1) = ispec
      else
        ! inner element
        iinner = iinner + 1
        phase_ispec_inner_infinite(iinner,2) = ispec
      endif
    enddo

    ! user output
    if (myrank == 0) then
      percentage_edge = 100. * nspec_outer_infinite / real(NSPEC_INFINITE)
      write(IMAIN,*) '  percentage of edge elements in infinite region ',percentage_edge,'%'
      write(IMAIN,*) '  percentage of volume elements in infinite region ',100. - percentage_edge,'%'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! debug: saves element flags
    if (DEBUG) then
      write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_innerouter_infinite_proc',myrank
      call write_VTK_data_elem_l(NSPEC_INFINITE,NGLOB_INFINITE, &
                                xstore_glob,ystore_glob,zstore_glob, &
                                ibool,is_on_a_slice_edge,filename)
    endif

  end select

  end subroutine setup_inner_outer
