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

  subroutine setup_color_perm(iregion_code)

  use constants, only: myrank

  use meshfem3D_par, only: &
    IMAIN,USE_MESH_COLORING_GPU,SAVE_MESH_FILES, &
    IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE

  use meshfem3D_par, only: ibool,is_on_a_slice_edge

  use MPI_crust_mantle_par
  use MPI_outer_core_par
  use MPI_inner_core_par

  implicit none

  integer,intent(in) :: iregion_code

  ! local parameters
  integer, dimension(:), allocatable :: perm
  integer :: ier
  integer :: nspec,nglob
  integer :: idomain

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  mesh coloring: ',USE_MESH_COLORING_GPU
    call flush_IMAIN()
  endif

  select case (iregion_code)
  case (IREGION_CRUST_MANTLE)
    ! crust mantle
    ! initializes
    num_colors_outer_crust_mantle = 0
    num_colors_inner_crust_mantle = 0

    ! mesh coloring
    if (USE_MESH_COLORING_GPU) then

      ! user output
      if (myrank == 0) write(IMAIN,*) '  coloring crust mantle... '

      ! crust/mantle region
      nspec = NSPEC_CRUST_MANTLE
      nglob = NGLOB_CRUST_MANTLE
      idomain = IREGION_CRUST_MANTLE

      ! creates coloring of elements
      allocate(perm(nspec),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating temporary perm crust mantle array')
      perm(:) = 0

      call setup_color(nspec,nglob,ibool,perm, &
                      idomain,is_on_a_slice_edge, &
                      num_phase_ispec_crust_mantle,phase_ispec_inner_crust_mantle, &
                      SAVE_MESH_FILES)

      ! checks
      if (minval(perm) /= 1) &
        call exit_MPI(myrank, 'minval(perm) should be 1')
      if (maxval(perm) /= num_phase_ispec_crust_mantle) &
        call exit_MPI(myrank, 'maxval(perm) should be num_phase_ispec_crust_mantle')

      ! sorts array according to permutation
      call synchronize_all()
      if (myrank == 0) then
        write(IMAIN,*) '     mesh permutation:'
      endif
      call setup_permutation(nspec,nglob,ibool, &
                            idomain,perm, &
                            num_colors_outer_crust_mantle,num_colors_inner_crust_mantle, &
                            num_elem_colors_crust_mantle, &
                            num_phase_ispec_crust_mantle,phase_ispec_inner_crust_mantle, &
                            SAVE_MESH_FILES)

      deallocate(perm)
    else
      ! dummy array
      allocate(num_elem_colors_crust_mantle(num_colors_outer_crust_mantle+num_colors_inner_crust_mantle),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_crust_mantle array')
    endif

  case (IREGION_OUTER_CORE)
    ! outer core
    ! initializes
    num_colors_outer_outer_core = 0
    num_colors_inner_outer_core = 0

    ! mesh coloring
    if (USE_MESH_COLORING_GPU) then

      ! user output
      if (myrank == 0) write(IMAIN,*) '  coloring outer core... '

      ! outer core region
      nspec = NSPEC_OUTER_CORE
      nglob = NGLOB_OUTER_CORE
      idomain = IREGION_OUTER_CORE

      ! creates coloring of elements
      allocate(perm(nspec),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating temporary perm outer_core array')
      perm(:) = 0

      call setup_color(nspec,nglob,ibool,perm, &
                      idomain,is_on_a_slice_edge, &
                      num_phase_ispec_outer_core,phase_ispec_inner_outer_core, &
                      SAVE_MESH_FILES)

      ! checks
      if (minval(perm) /= 1) &
        call exit_MPI(myrank, 'minval(perm) should be 1')
      if (maxval(perm) /= num_phase_ispec_outer_core) &
        call exit_MPI(myrank, 'maxval(perm) should be num_phase_ispec_outer_core')

      ! sorts array according to permutation
      call synchronize_all()
      if (myrank == 0) then
        write(IMAIN,*) '     mesh permutation:'
      endif
      call setup_permutation(nspec,nglob,ibool, &
                            idomain,perm, &
                            num_colors_outer_outer_core,num_colors_inner_outer_core, &
                            num_elem_colors_outer_core, &
                            num_phase_ispec_outer_core,phase_ispec_inner_outer_core, &
                            SAVE_MESH_FILES)

      deallocate(perm)
    else
      ! dummy array
      allocate(num_elem_colors_outer_core(num_colors_outer_outer_core+num_colors_inner_outer_core),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_outer_core array')
    endif

  case (IREGION_INNER_CORE)
    ! inner core
    ! initializes
    num_colors_outer_inner_core = 0
    num_colors_inner_inner_core = 0

    ! mesh coloring
    if (USE_MESH_COLORING_GPU) then

      ! user output
      if (myrank == 0) write(IMAIN,*) '  coloring inner core... '

      ! inner core region
      nspec = NSPEC_INNER_CORE
      nglob = NGLOB_INNER_CORE
      idomain = IREGION_INNER_CORE

      ! creates coloring of elements
      allocate(perm(nspec),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating temporary perm inner_core array')
      perm(:) = 0

      call setup_color(nspec,nglob,ibool,perm, &
                      idomain,is_on_a_slice_edge, &
                      num_phase_ispec_inner_core,phase_ispec_inner_inner_core, &
                      SAVE_MESH_FILES)

      ! checks
      ! inner core contains fictitious elements not counted for
      if (minval(perm) < 0) &
        call exit_MPI(myrank, 'minval(perm) should be at least 0')
      if (maxval(perm) > num_phase_ispec_inner_core) then
        print *,'Error perm inner core:',minval(perm),maxval(perm),num_phase_ispec_inner_core
        call exit_MPI(myrank, 'maxval(perm) should be num_phase_ispec_inner_core')
      endif

      ! sorts array according to permutation
      call synchronize_all()
      if (myrank == 0) then
        write(IMAIN,*) '     mesh permutation:'
      endif
      call setup_permutation(nspec,nglob,ibool, &
                            idomain,perm, &
                            num_colors_outer_inner_core,num_colors_inner_inner_core, &
                            num_elem_colors_inner_core, &
                            num_phase_ispec_inner_core,phase_ispec_inner_inner_core, &
                            SAVE_MESH_FILES)

      deallocate(perm)
    else
      ! dummy array
      allocate(num_elem_colors_inner_core(num_colors_outer_inner_core+num_colors_inner_inner_core),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_inner_core array')
    endif

  end select

  end subroutine setup_color_perm

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_color(nspec,nglob,ibool,perm, &
                            idomain,is_on_a_slice_edge, &
                            num_phase_ispec_d,phase_ispec_inner_d, &
                            SAVE_MESH_FILES)

! sets up mesh coloring

  use constants, only: myrank

  use meshfem3D_par, only: &
    LOCAL_PATH,MAX_NUMBER_OF_COLORS,IMAIN,NGLLX,NGLLY,NGLLZ,IFLAG_IN_FICTITIOUS_CUBE, &
    IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE,MAX_STRING_LEN,IOUT

  use meshfem3D_par, only: &
    idoubling,xstore_glob,ystore_glob,zstore_glob

  use MPI_crust_mantle_par, only: &
    num_colors_outer_crust_mantle,num_colors_inner_crust_mantle,num_elem_colors_crust_mantle

  use MPI_outer_core_par, only: &
    num_colors_outer_outer_core,num_colors_inner_outer_core,num_elem_colors_outer_core

  use MPI_inner_core_par, only: &
    num_colors_outer_inner_core,num_colors_inner_inner_core,num_elem_colors_inner_core

  implicit none

  integer :: nspec,nglob
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  integer, dimension(nspec) :: perm

  ! wrapper array for ispec is in domain:
  ! idomain: 1 == crust/mantle, 2 == outer core, 3 == inner core
  integer :: idomain
  logical, dimension(nspec) :: is_on_a_slice_edge
  integer :: num_phase_ispec_d
  integer, dimension(num_phase_ispec_d,2) :: phase_ispec_inner_d

  logical :: SAVE_MESH_FILES

  ! local parameters
  ! added for color permutation
  integer :: nb_colors_outer_elements,nb_colors_inner_elements
  integer, dimension(:), allocatable :: num_of_elems_in_this_color
  integer, dimension(:), allocatable :: color
  integer, dimension(:), allocatable :: first_elem_number_in_this_color
  logical, dimension(:), allocatable :: ispec_is_d

  integer :: nspec_outer,nspec_inner,nspec_domain
  integer :: nspec_outer_min_global,nspec_outer_max_global
  integer :: nspec_inner_min_global,nspec_inner_max_global
  integer :: min_elem_global,max_elem_global

  integer :: nb_colors
  integer :: nb_colors_min,nb_colors_max

  integer :: icolor,ispec,ispec_counter
  integer :: ispec_inner,ispec_outer
  integer :: ier

  character(len=2),dimension(3) :: str_domain = (/ "cm", "oc", "ic" /)
  character(len=MAX_STRING_LEN) :: filename
  character(len=MAX_STRING_LEN) :: prname

  ! debug file output
  logical, parameter :: DEBUG = .false.
  ! debug coloring : creates dummy mesh coloring, separating only inner/outer elements into colors
  logical, parameter :: DEBUG_COLOR = .false.

  !!!! David Michea: detection of the edges, coloring and permutation separately

  ! implement mesh coloring for GPUs if needed, to create subsets of disconnected elements
  ! to remove dependencies and the need for atomic operations in the sum of
  ! elemental contributions in the solver

  ! allocates temporary array with colors
  allocate(color(nspec),stat=ier)
  if (ier /= 0 ) stop 'Error allocating temporary color array'
  allocate(first_elem_number_in_this_color(MAX_NUMBER_OF_COLORS + 1),stat=ier)
  if (ier /= 0 ) stop 'Error allocating first_elem_number_in_this_color array'

  ! flags for elements in this domain
  ! for compatibility with SPECFEM3D mesh coloring routine
  allocate(ispec_is_d(nspec),stat=ier)
  if (ier /= 0 ) stop 'Error allocating ispec_is_d array'

  ! sets up domain coloring arrays
  select case (idomain)
  case (IREGION_CRUST_MANTLE,IREGION_OUTER_CORE)
    ! crust/mantle and outer core region meshes use all elements
    ispec_is_d(:) = .true.
  case (IREGION_INNER_CORE)
    ! initializes
    ispec_is_d(:) = .true.
    ! excludes fictitious elements from coloring
    where(idoubling == IFLAG_IN_FICTITIOUS_CUBE) ispec_is_d = .false.
    ! checks
    if (count(ispec_is_d) == 0) then
      stop 'Error no inner core elements'
    endif
  case default
    stop 'Error idomain in setup_color'
  end select

  ! fast element coloring scheme
  call get_perm_color_faster(is_on_a_slice_edge,ispec_is_d, &
                            ibool,perm,color, &
                            nspec,nglob, &
                            nb_colors_outer_elements,nb_colors_inner_elements, &
                            nspec_outer,nspec_inner,nspec_domain, &
                            first_elem_number_in_this_color)

  ! debug: file output
  if (SAVE_MESH_FILES .and. DEBUG .and. idomain == IREGION_CRUST_MANTLE) then
    call create_name_database(prname,myrank,idomain,LOCAL_PATH)
    filename = prname(1:len_trim(prname))//'color_'//str_domain(idomain)
    call write_VTK_data_elem_i(nspec,nglob, &
                               xstore_glob,ystore_glob,zstore_glob, &
                               ibool,color,filename)
  endif
  deallocate(color)

  ! for the last color, the next color is fictitious and its first (fictitious) element number is nspec + 1
  first_elem_number_in_this_color(nb_colors_outer_elements + nb_colors_inner_elements + 1) &
    = nspec_domain + 1

  allocate(num_of_elems_in_this_color(nb_colors_outer_elements + nb_colors_inner_elements),stat=ier)
  if (ier /= 0) then
    print *,'Error',myrank,' allocating num_of_elems_in_this_color:',nb_colors_outer_elements,nb_colors_inner_elements, &
          nb_colors_outer_elements + nb_colors_inner_elements
    call exit_MPI(myrank,'Error allocating num_of_elems_in_this_color array')
  endif

  num_of_elems_in_this_color(:) = 0
  do icolor = 1, nb_colors_outer_elements + nb_colors_inner_elements
    num_of_elems_in_this_color(icolor) = first_elem_number_in_this_color(icolor+1) - first_elem_number_in_this_color(icolor)
  enddo
  deallocate(first_elem_number_in_this_color)

  ! check that the sum of all the numbers of elements found in each color is equal
  ! to the total number of elements in the mesh
  if (sum(num_of_elems_in_this_color) /= nspec_domain) then
    print *,'Error number of elements in this color:',idomain
    print *,'rank: ',myrank,' nspec = ',nspec_domain
    print *,'  total number of elements in all the colors of the mesh = ', &
      sum(num_of_elems_in_this_color)
    call exit_MPI(myrank, 'incorrect total number of elements in all the colors of the mesh')
  endif

  ! check that the sum of all the numbers of elements found in each color for the outer elements is equal
  ! to the total number of outer elements found in the mesh
  if (sum(num_of_elems_in_this_color(1:nb_colors_outer_elements)) /= nspec_outer) then
    print *,'Error number of outer elements in this color:',idomain
    print *,'rank: ',myrank,' nspec_outer = ',nspec_outer
    print *,'nb_colors_outer_elements = ',nb_colors_outer_elements
    print *,'total number of elements in all the colors of the mesh for outer elements = ', &
      sum(num_of_elems_in_this_color(1:nb_colors_outer_elements))
    call exit_MPI(myrank, 'incorrect total number of elements in all the colors of the mesh for outer elements')
  endif

  ! debug: no mesh coloring, only creates dummy coloring arrays
  if (DEBUG_COLOR) then
    nb_colors_outer_elements = 0
    nb_colors_inner_elements = 0
    ispec_counter = 0

    ! first generate all the outer elements
    do ispec = 1,nspec
      if (ispec_is_d(ispec)) then
        if (is_on_a_slice_edge(ispec) .eqv. .true.) then
          ispec_counter = ispec_counter + 1
          perm(ispec) = ispec_counter
        endif
      endif
    enddo

    ! store total number of outer elements
    nspec_outer = ispec_counter

    ! only single color
    if (nspec_outer > 0 ) nb_colors_outer_elements = 1

    ! then generate all the inner elements
    do ispec = 1,nspec
      if (ispec_is_d(ispec)) then
        if (is_on_a_slice_edge(ispec) .eqv. .false.) then
          ispec_counter = ispec_counter + 1
          perm(ispec) = ispec_counter - nspec_outer ! starts again at 1
        endif
      endif
    enddo
    nspec_inner = ispec_counter - nspec_outer

    ! only single color
    if (nspec_inner > 0 ) nb_colors_inner_elements = 1

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'debugging mesh coloring:'
      write(IMAIN,*) 'nb_colors inner / outer: ',nb_colors_inner_elements,nb_colors_outer_elements
    endif

    ! re-allocate
    if (allocated(num_of_elems_in_this_color) ) deallocate(num_of_elems_in_this_color)
    allocate(num_of_elems_in_this_color(nb_colors_outer_elements + nb_colors_inner_elements),stat=ier)
    if (ier /= 0) then
      print *,'Error',myrank,' allocating num_of_elems_in_this_color:',nb_colors_outer_elements,nb_colors_inner_elements, &
          nb_colors_outer_elements + nb_colors_inner_elements
      call exit_MPI(myrank,'Error allocating num_of_elems_in_this_color array')
    endif

    if (nspec_outer > 0 ) num_of_elems_in_this_color(1) = nspec_outer
    if (nspec_inner > 0 ) num_of_elems_in_this_color(2) = nspec_inner
  endif ! debug_color

  ! debug: saves mesh coloring numbers into files
  if (DEBUG) then
    ! debug file output
    call create_name_database(prname,myrank,idomain,LOCAL_PATH)
    filename = prname(1:len_trim(prname))//'num_of_elems_in_this_color_'//str_domain(idomain)//'.dat'
    open(unit=IOUT,file=trim(filename),status='unknown',iostat=ier)
    if (ier /= 0 ) stop 'Error opening num_of_elems_in_this_color file'
    ! number of colors for outer elements
    write(IOUT,*) nb_colors_outer_elements
    ! number of colors for inner elements
    write(IOUT,*) nb_colors_inner_elements
    ! number of elements in each color
    ! outer elements
    do icolor = 1, nb_colors_outer_elements + nb_colors_inner_elements
      write(IOUT,*) num_of_elems_in_this_color(icolor)
    enddo
    close(IOUT)
  endif

  ! checks non-zero elements in colors
  do icolor = 1,nb_colors_outer_elements + nb_colors_inner_elements
    ! checks
    if (num_of_elems_in_this_color(icolor) == 0) then
      print *,'rank: ',myrank,'domain:',idomain,' nspec = ',nspec_domain
      print *,'Error zero elements in this color:',icolor
      print *,'total number of elements in all the colors of the mesh = ', &
        sum(num_of_elems_in_this_color)
      call exit_MPI(myrank, 'zero elements in a color of the mesh')
    endif
  enddo



  ! sets up domain coloring arrays
  select case (idomain)
  case (IREGION_CRUST_MANTLE)
    ! crust/mantle domains
    num_colors_outer_crust_mantle = nb_colors_outer_elements
    num_colors_inner_crust_mantle = nb_colors_inner_elements

    allocate(num_elem_colors_crust_mantle(num_colors_outer_crust_mantle + num_colors_inner_crust_mantle),stat=ier)
    if (ier /= 0 ) stop 'Error allocating num_elem_colors_crust_mantle array'

    num_elem_colors_crust_mantle(:) = num_of_elems_in_this_color(:)

  case (IREGION_OUTER_CORE)
    ! outer core domains
    num_colors_outer_outer_core = nb_colors_outer_elements
    num_colors_inner_outer_core = nb_colors_inner_elements

    allocate(num_elem_colors_outer_core(num_colors_outer_outer_core + num_colors_inner_outer_core),stat=ier)
    if (ier /= 0 ) stop 'Error allocating num_elem_colors_outer_core array'

    num_elem_colors_outer_core(:) = num_of_elems_in_this_color(:)

  case (IREGION_INNER_CORE)
    ! inner core domains
    num_colors_outer_inner_core = nb_colors_outer_elements
    num_colors_inner_inner_core = nb_colors_inner_elements

    allocate(num_elem_colors_inner_core(num_colors_outer_inner_core + num_colors_inner_inner_core),stat=ier)
    if (ier /= 0 ) stop 'Error allocating num_elem_colors_inner_core array'

    num_elem_colors_inner_core(:) = num_of_elems_in_this_color(:)

  case default
    stop 'Error idomain not recognized'
  end select

  ! sets up elements for loops in simulations
  do ispec = 1, nspec
    ! only elements in this domain
    if (ispec_is_d(ispec)) then

      ! sets phase_ispec arrays with ordering of elements
      if (is_on_a_slice_edge(ispec) .eqv. .true.) then
        ! outer elements
        ispec_outer = perm(ispec)

        ! checks
        if (ispec_outer < 1 .or. ispec_outer > num_phase_ispec_d) then
          print *,'Error outer permutation:',idomain
          print *,'rank:',myrank,'  ispec_outer = ',ispec_outer
          print *,'num_phase_ispec_d = ',num_phase_ispec_d
          call exit_MPI(myrank,'Error outer permutation')
        endif

        phase_ispec_inner_d(ispec_outer,1) = ispec

      else
        ! inner elements
        ispec_inner = perm(ispec)

        ! checks
        if (ispec_inner < 1 .or. ispec_inner > num_phase_ispec_d) then
          print *,'Error inner permutation:',idomain
          print *,'rank:',myrank,'  ispec_inner = ',ispec_inner
          print *,'num_phase_ispec_d = ',num_phase_ispec_d
          call exit_MPI(myrank,'Error inner permutation')
        endif

        phase_ispec_inner_d(ispec_inner,2) = ispec

      endif
    endif
  enddo

  ! total number of colors
  nb_colors = nb_colors_inner_elements + nb_colors_outer_elements
  call min_all_i(nb_colors,nb_colors_min)
  call max_all_i(nb_colors,nb_colors_max)

  ! min/max of elements per color
  call min_all_i(minval(num_of_elems_in_this_color(:)),min_elem_global)
  call max_all_i(maxval(num_of_elems_in_this_color(:)),max_elem_global)

  ! min/max of inner/outer elements
  call min_all_i(nspec_inner,nspec_inner_min_global)
  call max_all_i(nspec_inner,nspec_inner_max_global)
  call min_all_i(nspec_outer,nspec_outer_min_global)
  call max_all_i(nspec_outer,nspec_outer_max_global)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     total colors:'
    write(IMAIN,*) '       total colors min/max = ',nb_colors_min,nb_colors_max
    write(IMAIN,*) '       elements per color min/max = ',min_elem_global,max_elem_global
    write(IMAIN,*) '       inner elements min/max = ',nspec_inner_min_global,nspec_inner_max_global
    write(IMAIN,*) '       outer elements min/max = ',nspec_outer_min_global,nspec_outer_max_global
    call flush_IMAIN()
  endif

  ! debug: outputs permutation array as VTK file
  if (DEBUG .and. idomain == IREGION_CRUST_MANTLE) then
    call create_name_database(prname,myrank,idomain,LOCAL_PATH)
    filename = prname(1:len_trim(prname))//'perm_'//str_domain(idomain)
    call write_VTK_data_elem_i(nspec,nglob, &
                               xstore_glob,ystore_glob,zstore_glob, &
                               ibool,perm,filename)
  endif

  deallocate(ispec_is_d)
  deallocate(num_of_elems_in_this_color)

  end subroutine setup_color

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_permutation(nspec,nglob,ibool, &
                              idomain,perm, &
                              num_colors_outer,num_colors_inner, &
                              num_elem_colors, &
                              num_phase_ispec_d,phase_ispec_inner_d, &
                              SAVE_MESH_FILES)

  use constants

  use meshfem3D_models_par, only: &
    TRANSVERSE_ISOTROPY,HETEROGEN_3D_MANTLE,ANISOTROPIC_3D_MANTLE, &
    ANISOTROPIC_INNER_CORE,ATTENUATION,SAVE_BOUNDARY_MESH, &
    ATTENUATION_3D,ATTENUATION_1D_WITH_3D_STORAGE

  use meshfem3D_par, only: &
    ABSORBING_CONDITIONS, &
    LOCAL_PATH, &
    NCHUNKS,NSPEC2D_TOP,NSPEC2D_BOTTOM, &
    xstore,ystore,zstore,idoubling,xstore_glob,ystore_glob,zstore_glob

  use regions_mesh_par2, only: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
    gammaxstore,gammaystore,gammazstore, &
    rhostore,dvpstore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
    ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
    rho_vp,rho_vs, &
    nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
    ispec_is_tiso,tau_e_store,Qmu_store, &
    NSPEC2D_MOHO, NSPEC2D_400, NSPEC2D_670, &
    ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot, &
    ibelm_670_top,ibelm_670_bot

  use MPI_crust_mantle_par, only: NSPEC_CRUST_MANTLE
  use MPI_outer_core_par, only: NSPEC_OUTER_CORE
  use MPI_inner_core_par, only: NSPEC_INNER_CORE

  implicit none

  integer,intent(in) :: nspec,nglob
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  integer,intent(in) :: idomain
  integer, dimension(nspec),intent(inout) :: perm

  integer :: num_colors_outer,num_colors_inner
  integer, dimension(num_colors_outer + num_colors_inner) :: num_elem_colors
  integer :: num_phase_ispec_d
  integer, dimension(num_phase_ispec_d,2) :: phase_ispec_inner_d

  logical :: SAVE_MESH_FILES

  ! local parameters
  ! added for sorting
  double precision, dimension(:,:,:,:), allocatable :: temp_array_dble
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: temp_array_real
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: temp_array_real_sls
  integer, dimension(:,:,:,:), allocatable :: temp_array_int
  integer, dimension(:), allocatable :: temp_array_int_1D
  integer, dimension(:), allocatable :: temp_perm_global
  logical, dimension(:), allocatable :: temp_array_logical_1D
  logical, dimension(:), allocatable :: mask_global

  integer :: icolor,icounter,ispec,ielem,ier,i
  integer :: iface,old_ispec,new_ispec

  character(len=MAX_STRING_LEN) :: filename
  character(len=MAX_STRING_LEN) :: prname

  ! debug file output
  logical,parameter :: DEBUG = .false.

  ! sorts array according to permutation
  allocate(temp_perm_global(nspec),stat=ier)
  if (ier /= 0 ) stop 'Error temp_perm_global array'

  ! global ordering
  temp_perm_global(:) = 0
  icounter = 0

  ! fills global permutation array

  ! first outer elements coloring
  ! phase element counter
  ielem = 0
  do icolor = 1,num_colors_outer
    ! loops through elements
    do i = 1,num_elem_colors(icolor)
      ielem = ielem + 1
      ispec = phase_ispec_inner_d(ielem,1) ! 1 -- first phase, outer elements
      ! reorders elements
      icounter = icounter + 1
      temp_perm_global(ispec) = icounter
      ! resets to new order
      phase_ispec_inner_d(ielem,1) = icounter
    enddo
  enddo
  ! inner elements coloring
  ielem = 0
  do icolor = num_colors_outer+1,num_colors_outer+num_colors_inner
    ! loops through elements
    do i = 1,num_elem_colors(icolor)
      ielem = ielem + 1
      ispec = phase_ispec_inner_d(ielem,2) ! 2 -- second phase, inner elements
      ! reorders elements
      icounter = icounter + 1
      temp_perm_global(ispec) = icounter
      ! resets to new order
      phase_ispec_inner_d(ielem,2) = icounter
    enddo
  enddo

  ! handles fictitious cube elements for inner core
  ! which contains fictitious elements not counted for
  if (idomain == IREGION_INNER_CORE) then
    ! fills up permutation with fictitious numbering
    do ispec = 1,nspec
      if (temp_perm_global(ispec) == 0) then
        icounter = icounter + 1
        temp_perm_global(ispec) = icounter
      endif
    enddo
  endif

  ! checks counter
  if (icounter /= nspec) then
    print *,'Error temp perm: ',icounter,nspec
    stop 'Error temporary global permutation incomplete'
  endif
  ! checks values
  if (minval(temp_perm_global) /= 1) call exit_MPI(myrank, 'minval(temp_perm_global) should be 1')
  if (maxval(temp_perm_global) /= nspec) call exit_MPI(myrank, 'maxval(temp_perm_global) should be nspec')

  ! checks if every element was uniquely set
  allocate(mask_global(nspec),stat=ier)
  if (ier /= 0 ) stop 'Error allocating temporary mask_global'
  mask_global(:) = .false.

  icounter = 0 ! counts permutations
  do ispec = 1, nspec
    new_ispec = temp_perm_global(ispec)
    ! checks bounds
    if (new_ispec < 1 .or. new_ispec > nspec ) call exit_MPI(myrank,'Error temp_perm_global ispec bounds')
    ! checks if already set
    if (mask_global(new_ispec)) then
      print *,'Error temp_perm_global:',ispec,new_ispec,'element already set'
      call exit_MPI(myrank,'Error global permutation')
    else
      mask_global(new_ispec) = .true.
    endif
    ! counts permutations
    if (new_ispec /= ispec ) icounter = icounter + 1
  enddo

  ! checks number of set elements
  if (count(mask_global(:)) /= nspec) then
    print *,'Error temp_perm_global:',count(mask_global(:)),nspec,'permutation incomplete'
    call exit_MPI(myrank,'Error global permutation incomplete')
  endif
  deallocate(mask_global)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '       number of permutations = ',icounter
    call flush_IMAIN()
  endif

  ! outputs permutation array as VTK file
  if (SAVE_MESH_FILES .and. DEBUG .and. idomain == IREGION_CRUST_MANTLE) then
    call create_name_database(prname,myrank,idomain,LOCAL_PATH)
    filename = prname(1:len_trim(prname))//'perm_global'
    call write_VTK_data_elem_i(nspec,nglob, &
                               xstore_glob,ystore_glob,zstore_glob, &
                               ibool,temp_perm_global,filename)
  endif

  ! store as new permutation
  perm(:) = temp_perm_global(:)
  deallocate(temp_perm_global)

  ! permutes all required mesh arrays according to new ordering

  ! permutation of ibool
  allocate(temp_array_int(NGLLX,NGLLY,NGLLZ,nspec))
  call permute_elements_integer(ibool,temp_array_int,perm,nspec)
  deallocate(temp_array_int)

  ! element idoubling flags
  allocate(temp_array_int_1D(nspec))
  call permute_elements_integer1D(idoubling,temp_array_int_1D,perm,nspec)
  deallocate(temp_array_int_1D)

  ! element domain flags
  allocate(temp_array_logical_1D(nspec))
  call permute_elements_logical1D(ispec_is_tiso,temp_array_logical_1D,perm,nspec)
  deallocate(temp_array_logical_1D)

  ! mesh arrays
  ! double precision
  allocate(temp_array_dble(NGLLX,NGLLY,NGLLZ,nspec))
  call permute_elements_dble(xstore,temp_array_dble,perm,nspec)
  call permute_elements_dble(ystore,temp_array_dble,perm,nspec)
  call permute_elements_dble(zstore,temp_array_dble,perm,nspec)
  deallocate(temp_array_dble)
  ! custom precision
  allocate(temp_array_real(NGLLX,NGLLY,NGLLZ,nspec))
  call permute_elements_real(xixstore,temp_array_real,perm,nspec)
  call permute_elements_real(xiystore,temp_array_real,perm,nspec)
  call permute_elements_real(xizstore,temp_array_real,perm,nspec)
  call permute_elements_real(etaxstore,temp_array_real,perm,nspec)
  call permute_elements_real(etaystore,temp_array_real,perm,nspec)
  call permute_elements_real(etazstore,temp_array_real,perm,nspec)
  call permute_elements_real(gammaxstore,temp_array_real,perm,nspec)
  call permute_elements_real(gammaystore,temp_array_real,perm,nspec)
  call permute_elements_real(gammazstore,temp_array_real,perm,nspec)

  ! material parameters
  call permute_elements_real(rhostore,temp_array_real,perm,nspec)
  call permute_elements_real(kappavstore,temp_array_real,perm,nspec)
  deallocate(temp_array_real)

  ! boundary surfaces
  ! note: only arrays pointing to ispec will have to be permuted since value of ispec will be different
  !
  ! xmin
  do iface = 1,nspec2D_xmin
      old_ispec = ibelm_xmin(iface)
      new_ispec = perm(old_ispec)
      ibelm_xmin(iface) = new_ispec
  enddo
  ! xmax
  do iface = 1,nspec2D_xmax
      old_ispec = ibelm_xmax(iface)
      new_ispec = perm(old_ispec)
      ibelm_xmax(iface) = new_ispec
  enddo
  ! ymin
  do iface = 1,nspec2D_ymin
      old_ispec = ibelm_ymin(iface)
      new_ispec = perm(old_ispec)
      ibelm_ymin(iface) = new_ispec
  enddo
  ! ymax
  do iface = 1,nspec2D_ymax
      old_ispec = ibelm_ymax(iface)
      new_ispec = perm(old_ispec)
      ibelm_ymax(iface) = new_ispec
  enddo
  ! bottom
  do iface = 1,NSPEC2D_BOTTOM(idomain)
      old_ispec = ibelm_bottom(iface)
      new_ispec = perm(old_ispec)
      ibelm_bottom(iface) = new_ispec
  enddo
  ! top
  do iface = 1,NSPEC2D_TOP(idomain)
      old_ispec = ibelm_top(iface)
      new_ispec = perm(old_ispec)
      ibelm_top(iface) = new_ispec
  enddo

  ! attenuation arrays
  if (ATTENUATION) then
    if (ATTENUATION_3D .or. ATTENUATION_1D_WITH_3D_STORAGE) then
      allocate(temp_array_real(NGLLX,NGLLY,NGLLZ,nspec))
      allocate(temp_array_real_sls(NGLLX,NGLLY,NGLLZ,N_SLS,nspec))
      call permute_elements_real(Qmu_store,temp_array_real,perm,nspec)
      call permute_elements_real_sls(tau_e_store,temp_array_real_sls,perm,nspec)
      deallocate(temp_array_real,temp_array_real_sls)
    else
      allocate(temp_array_real(1,1,1,nspec))
      allocate(temp_array_real_sls(1,1,1,N_SLS,nspec))
      call permute_elements_real1(Qmu_store,temp_array_real,perm,nspec)
      call permute_elements_real_sls1(tau_e_store,temp_array_real_sls,perm,nspec)
      deallocate(temp_array_real,temp_array_real_sls)
    endif
  endif

  select case (idomain)
  case (IREGION_CRUST_MANTLE)
    ! checks number of elements
    if (nspec /= NSPEC_CRUST_MANTLE ) &
      call exit_MPI(myrank,'Error in permutation nspec should be NSPEC_CRUST_MANTLE')

    allocate(temp_array_real(NGLLX,NGLLY,NGLLZ,nspec))

    if (ANISOTROPIC_3D_MANTLE) then
      call permute_elements_real(c11store,temp_array_real,perm,nspec)
      call permute_elements_real(c11store,temp_array_real,perm,nspec)
      call permute_elements_real(c12store,temp_array_real,perm,nspec)
      call permute_elements_real(c13store,temp_array_real,perm,nspec)
      call permute_elements_real(c14store,temp_array_real,perm,nspec)
      call permute_elements_real(c15store,temp_array_real,perm,nspec)
      call permute_elements_real(c16store,temp_array_real,perm,nspec)
      call permute_elements_real(c22store,temp_array_real,perm,nspec)
      call permute_elements_real(c23store,temp_array_real,perm,nspec)
      call permute_elements_real(c24store,temp_array_real,perm,nspec)
      call permute_elements_real(c25store,temp_array_real,perm,nspec)
      call permute_elements_real(c26store,temp_array_real,perm,nspec)
      call permute_elements_real(c33store,temp_array_real,perm,nspec)
      call permute_elements_real(c34store,temp_array_real,perm,nspec)
      call permute_elements_real(c35store,temp_array_real,perm,nspec)
      call permute_elements_real(c36store,temp_array_real,perm,nspec)
      call permute_elements_real(c44store,temp_array_real,perm,nspec)
      call permute_elements_real(c45store,temp_array_real,perm,nspec)
      call permute_elements_real(c46store,temp_array_real,perm,nspec)
      call permute_elements_real(c55store,temp_array_real,perm,nspec)
      call permute_elements_real(c56store,temp_array_real,perm,nspec)
      call permute_elements_real(c66store,temp_array_real,perm,nspec)
    else
      call permute_elements_real(muvstore,temp_array_real,perm,nspec)

      if (TRANSVERSE_ISOTROPY) then
        call permute_elements_real(kappahstore,temp_array_real,perm,nspec)
        call permute_elements_real(muhstore,temp_array_real,perm,nspec)
        call permute_elements_real(eta_anisostore,temp_array_real,perm,nspec)
      endif
    endif

    if (HETEROGEN_3D_MANTLE) then
      call permute_elements_real(dvpstore,temp_array_real,perm,nspec)
    endif

    if (ABSORBING_CONDITIONS .and. NCHUNKS /= 6) then
      call permute_elements_real(rho_vp,temp_array_real,perm,nspec)
      call permute_elements_real(rho_vs,temp_array_real,perm,nspec)
    endif

    deallocate(temp_array_real)

    ! discontinuities boundary surface
    if (SAVE_BOUNDARY_MESH) then
      ! moho
      do iface = 1,nspec2D_MOHO
        ! top
        old_ispec = ibelm_moho_top(iface)
        new_ispec = perm(old_ispec)
        ibelm_moho_top(iface) = new_ispec
        ! bottom
        old_ispec = ibelm_moho_bot(iface)
        new_ispec = perm(old_ispec)
        ibelm_moho_bot(iface) = new_ispec
      enddo
      ! 400
      do iface = 1,nspec2D_400
        ! top
        old_ispec = ibelm_400_top(iface)
        new_ispec = perm(old_ispec)
        ibelm_400_top(iface) = new_ispec
        ! bottom
        old_ispec = ibelm_400_bot(iface)
        new_ispec = perm(old_ispec)
        ibelm_400_bot(iface) = new_ispec
      enddo
      ! 670
      do iface = 1,nspec2D_670
        ! top
        old_ispec = ibelm_670_top(iface)
        new_ispec = perm(old_ispec)
        ibelm_670_top(iface) = new_ispec
        ! bottom
        old_ispec = ibelm_670_bot(iface)
        new_ispec = perm(old_ispec)
        ibelm_670_bot(iface) = new_ispec
      enddo
    endif

  case (IREGION_OUTER_CORE)
    ! checks number of elements
    if (nspec /= NSPEC_OUTER_CORE ) &
      call exit_MPI(myrank,'Error in permutation nspec should be NSPEC_OUTER_CORE')

    if (ABSORBING_CONDITIONS .and. NCHUNKS /= 6) then
      allocate(temp_array_real(NGLLX,NGLLY,NGLLZ,nspec))

      call permute_elements_real(rho_vp,temp_array_real,perm,nspec)

      deallocate(temp_array_real)
    endif

  case (IREGION_INNER_CORE)
    ! checks number of elements
    if (nspec /= NSPEC_INNER_CORE ) &
      call exit_MPI(myrank,'Error in permutation nspec should be NSPEC_INNER_CORE')

    allocate(temp_array_real(NGLLX,NGLLY,NGLLZ,nspec))

    ! note: muvstore needed for attenuation also for anisotropic inner core
    call permute_elements_real(muvstore,temp_array_real,perm,nspec)

    !  anisotropy in the inner core only
    if (ANISOTROPIC_INNER_CORE) then
      call permute_elements_real(c11store,temp_array_real,perm,nspec)
      call permute_elements_real(c33store,temp_array_real,perm,nspec)
      call permute_elements_real(c12store,temp_array_real,perm,nspec)
      call permute_elements_real(c13store,temp_array_real,perm,nspec)
      call permute_elements_real(c44store,temp_array_real,perm,nspec)
    endif

    deallocate(temp_array_real)

  case default
    stop 'Error idomain in setup_permutation'
  end select

  end subroutine setup_permutation

!
!-------------------------------------------------------------------------------------------------
!

! deprecated ...
!
!  subroutine setup_color_perm(iregion_code,nspec,nglob, &
!                              ibool,is_on_a_slice_edge,prname, &
!                              npoin2D_xi,npoin2D_eta)
!
!  use constants
!  use meshfem3D_par, only: NSTEP,DT,NPROC_XI,NPROC_ETA
!  implicit none
!
!  integer :: iregion_code
!
!  integer :: nspec,nglob
!  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
!
!  ! this for non blocking MPI
!  logical, dimension(nspec) :: is_on_a_slice_edge
!
!  ! name of the database file
!  character(len=MAX_STRING_LEN) :: prname
!
!  integer :: npoin2D_xi,npoin2D_eta
!
!  ! local parameters
!  integer :: nb_colors_outer_elements,nb_colors_inner_elements,nspec_outer
!  integer, dimension(:), allocatable :: perm
!  integer, dimension(:), allocatable :: first_elem_number_in_this_color
!  integer, dimension(:), allocatable :: num_of_elems_in_this_color
!
!  integer :: icolor,ispec_counter
!  integer :: nspec_outer_min_global,nspec_outer_max_global
!  integer :: ispec,ier
!
!  !!!! David Michea: detection of the edges, coloring and permutation separately
!  allocate(perm(nspec))
!
!  ! implement mesh coloring for GPUs if needed, to create subsets of disconnected elements
!  ! to remove dependencies and the need for atomic operations in the sum of elemental contributions in the solver
!  if (USE_MESH_COLORING_GPU) then
!
!    ! user output
!    if (myrank == 0 ) write(IMAIN,*) '  creating mesh coloring'
!
!    allocate(first_elem_number_in_this_color(MAX_NUMBER_OF_COLORS + 1))
!
!    call get_perm_color_faster(is_on_a_slice_edge,ibool,perm,nspec,nglob, &
!                              nb_colors_outer_elements,nb_colors_inner_elements,nspec_outer, &
!                              first_elem_number_in_this_color)
!
!    ! for the last color, the next color is fictitious and its first (fictitious) element number is nspec + 1
!    first_elem_number_in_this_color(nb_colors_outer_elements + nb_colors_inner_elements + 1) = nspec + 1
!
!    allocate(num_of_elems_in_this_color(nb_colors_outer_elements + nb_colors_inner_elements))
!
!    ! save mesh coloring
!    open(unit=IOUT,file=prname(1:len_trim(prname))//'num_of_elems_in_this_color.dat', &
!         status='unknown',iostat=ier)
!    if (ier /= 0 ) call exit_mpi(myrank,'Error opening num_of_elems_in_this_color file')
!
!    ! number of colors for outer elements
!    write(IOUT,*) nb_colors_outer_elements
!
!    ! number of colors for inner elements
!    write(IOUT,*) nb_colors_inner_elements
!
!    ! number of elements in each color
!    do icolor = 1, nb_colors_outer_elements + nb_colors_inner_elements
!      num_of_elems_in_this_color(icolor) = first_elem_number_in_this_color(icolor+1) &
!                                          - first_elem_number_in_this_color(icolor)
!      write(IOUT,*) num_of_elems_in_this_color(icolor)
!    enddo
!    close(IOUT)
!
!    ! check that the sum of all the numbers of elements found in each color is equal
!    ! to the total number of elements in the mesh
!    if (sum(num_of_elems_in_this_color) /= nspec) then
!      print *,'nspec = ',nspec
!      print *,'total number of elements in all the colors of the mesh = ',sum(num_of_elems_in_this_color)
!      call exit_mpi(myrank,'incorrect total number of elements in all the colors of the mesh')
!    endif
!
!    ! check that the sum of all the numbers of elements found in each color for the outer elements is equal
!    ! to the total number of outer elements found in the mesh
!    if (sum(num_of_elems_in_this_color(1:nb_colors_outer_elements)) /= nspec_outer) then
!      print *,'nspec_outer = ',nspec_outer
!      print *,'total number of elements in all the colors of the mesh for outer elements = ', &
!        sum(num_of_elems_in_this_color)
!      call exit_mpi(myrank,'incorrect total number of elements in all the colors of the mesh for outer elements')
!    endif
!
!    call MPI_ALLREDUCE(nspec_outer,nspec_outer_min_global,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ier)
!    call MPI_ALLREDUCE(nspec_outer,nspec_outer_max_global,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
!
!    deallocate(first_elem_number_in_this_color)
!    deallocate(num_of_elems_in_this_color)
!
!  else
!
!    !! DK DK for regular C + MPI version for CPUs: do not use colors but nonetheless put all the outer elements
!    !! DK DK first in order to be able to overlap non-blocking MPI communications with calculations
!
!    !! DK DK nov 2010, for Rosa Badia / StarSs:
!    !! no need for mesh coloring, but need to implement inner/outer subsets for non blocking MPI for StarSs
!    ispec_counter = 0
!    perm(:) = 0
!
!    ! first generate all the outer elements
!    do ispec = 1,nspec
!      if (is_on_a_slice_edge(ispec)) then
!        ispec_counter = ispec_counter + 1
!        perm(ispec) = ispec_counter
!      endif
!    enddo
!
!    ! make sure we have detected some outer elements
!    if (ispec_counter <= 0) stop 'fatal error: no outer elements detected!'
!
!    ! store total number of outer elements
!    nspec_outer = ispec_counter
!
!    ! then generate all the inner elements
!    do ispec = 1,nspec
!      if (.not. is_on_a_slice_edge(ispec)) then
!        ispec_counter = ispec_counter + 1
!        perm(ispec) = ispec_counter
!      endif
!    enddo
!
!    ! test that all the elements have been used once and only once
!    if (ispec_counter /= nspec) stop 'fatal error: ispec_counter not equal to nspec'
!
!    ! do basic checks
!    if (minval(perm) /= 1) stop 'minval(perm) should be 1'
!    if (maxval(perm) /= nspec) stop 'maxval(perm) should be nspec'
!
!    call MPI_ALLREDUCE(nspec_outer,nspec_outer_min_global,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ier)
!    call MPI_ALLREDUCE(nspec_outer,nspec_outer_max_global,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
!
!  endif ! USE_MESH_COLORING_GPU
!
!  !! DK DK and Manh Ha, Nov 2011: added this to use the new mesher in the CUDA or C / StarSs test codes
!
!  if (myrank == 0 .and. iregion_code == IREGION_CRUST_MANTLE) then
!    ! write a header file for the Fortran version of the solver
!    open(unit=IOUT,file=prname(1:len_trim(prname))//'values_from_mesher_f90.h', &
!          status='unknown',iostat=ier)
!    if (ier /= 0 ) call exit_mpi(myrank,'Error opening file values_from_mesher_f90.h')
!
!    write(IOUT,*) 'integer, parameter :: NSPEC = ',nspec
!    write(IOUT,*) 'integer, parameter :: NGLOB = ',nglob
!    !!! DK DK use 1000 time steps only for the scaling tests
!    write(IOUT,*) 'integer, parameter :: NSTEP = 1000 !!!!!!!!!!! ',nstep
!    write(IOUT,*) 'real(kind=4), parameter :: deltat = ',DT
!    write(IOUT,*)
!    write(IOUT,*) 'integer, parameter ::  NGLOB2DMAX_XMIN_XMAX = ',npoin2D_xi
!    write(IOUT,*) 'integer, parameter ::  NGLOB2DMAX_YMIN_YMAX = ',npoin2D_eta
!    write(IOUT,*) 'integer, parameter ::  NGLOB2DMAX_ALL = ',max(npoin2D_xi,npoin2D_eta)
!    write(IOUT,*) 'integer, parameter ::  NPROC_XI = ',NPROC_XI
!    write(IOUT,*) 'integer, parameter ::  NPROC_ETA = ',NPROC_ETA
!    write(IOUT,*)
!    write(IOUT,*) '! element number of the source and of the station'
!    write(IOUT,*) '! after permutation of the elements by mesh coloring'
!    write(IOUT,*) '! and inner/outer set splitting in the mesher'
!    write(IOUT,*) 'integer, parameter :: NSPEC_SOURCE = ',perm(NSPEC/3)
!    write(IOUT,*) 'integer, parameter :: RANK_SOURCE = 0'
!    write(IOUT,*)
!    write(IOUT,*) 'integer, parameter :: RANK_STATION = (NPROC_XI*NPROC_ETA - 1)'
!    write(IOUT,*) 'integer, parameter :: NSPEC_STATION = ',perm(2*NSPEC/3)
!
!    ! save coordinates of the seismic source
!    !   write(IOUT,*) xstore(2,2,2,10);
!    !   write(IOUT,*) ystore(2,2,2,10);
!    !   write(IOUT,*) zstore(2,2,2,10);
!
!    ! save coordinates of the seismic station
!    !   write(IOUT,*) xstore(2,2,2,nspec-10);
!    !   write(IOUT,*) ystore(2,2,2,nspec-10);
!    !   write(IOUT,*) zstore(2,2,2,nspec-10);
!    close(IOUT)
!
!    !! write a header file for the C version of the solver
!    open(unit=IOUT,file=prname(1:len_trim(prname))//'values_from_mesher_C.h', &
!          status='unknown',iostat=ier)
!    if (ier /= 0 ) call exit_mpi(myrank,'Error opening file values_from_mesher_C.h')
!
!    write(IOUT,*) '#define NSPEC ',nspec
!    write(IOUT,*) '#define NGLOB ',nglob
!    !!    write(IOUT,*) '#define NSTEP ',nstep
!    !!! DK DK use 1000 time steps only for the scaling tests
!    write(IOUT,*) '// #define NSTEP ',nstep
!    write(IOUT,*) '#define NSTEP 1000'
!    ! put an "f" at the end to force single precision
!    write(IOUT,"('#define deltat ',e18.10,'f')") DT
!    write(IOUT,*) '#define NGLOB2DMAX_XMIN_XMAX ',npoin2D_xi
!    write(IOUT,*) '#define NGLOB2DMAX_YMIN_YMAX ',npoin2D_eta
!    write(IOUT,*) '#define NGLOB2DMAX_ALL ',max(npoin2D_xi,npoin2D_eta)
!    write(IOUT,*) '#define NPROC_XI ',NPROC_XI
!    write(IOUT,*) '#define NPROC_ETA ',NPROC_ETA
!    write(IOUT,*)
!    write(IOUT,*) '// element and MPI slice number of the source and the station'
!    write(IOUT,*) '// after permutation of the elements by mesh coloring'
!    write(IOUT,*) '// and inner/outer set splitting in the mesher'
!    write(IOUT,*) '#define RANK_SOURCE 0'
!    write(IOUT,*) '#define NSPEC_SOURCE ',perm(NSPEC/3)
!    write(IOUT,*)
!    write(IOUT,*) '#define RANK_STATION (NPROC_XI*NPROC_ETA - 1)'
!    write(IOUT,*) '#define NSPEC_STATION ',perm(2*NSPEC/3)
!    close(IOUT)
!
!    open(unit=IOUT,file=prname(1:len_trim(prname))//'values_from_mesher_nspec_outer.h', &
!          status='unknown',iostat=ier)
!    if (ier /= 0 ) call exit_mpi(myrank,'Error opening values_from_mesher_nspec_outer.h file')
!
!    write(IOUT,*) '#define NSPEC_OUTER ',nspec_outer_max_global
!    write(IOUT,*) '// NSPEC_OUTER_min = ',nspec_outer_min_global
!    write(IOUT,*) '// NSPEC_OUTER_max = ',nspec_outer_max_global
!    close(IOUT)
!
!  endif
!
!  !! DK DK and Manh Ha, Nov 2011: added this to use the new mesher in the CUDA or C / StarSs test codes
!
!  deallocate(perm)
!
!
!  end subroutine setup_color_perm
