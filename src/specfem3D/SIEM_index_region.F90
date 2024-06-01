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


  subroutine SIEM_get_index_region()

  use constants, only: myrank,NGLLX,NGLLY,NGLLZ,NGLLCUBE,IMAIN,SIZE_INTEGER,IFLAG_IN_FICTITIOUS_CUBE, &
    NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF,ADD_TRINF

  use constants_solver, only: NSPEC_INNER_CORE, &
    NSPEC_OUTER_CORE,NSPEC_CRUST_MANTLE,NSPEC_TRINFINITE,NSPEC_INFINITE, &
    NGLOB_INNER_CORE,NGLOB_OUTER_CORE,NGLOB_CRUST_MANTLE,NGLOB_TRINFINITE, &
    NGLOB_INFINITE,NSPEC2D_BOTTOM_OC,NSPEC2D_BOTTOM_CM,NSPEC2D_BOTTOM_TRINF, &
    NSPEC2D_BOTTOM_INF

  use specfem_par, only: num_interfaces_inner_core, &
    max_nibool_interfaces_ic,my_neighbors_inner_core, &
    nibool_interfaces_inner_core,ibool_interfaces_inner_core

  use specfem_par, only: num_interfaces_outer_core, &
    max_nibool_interfaces_oc,my_neighbors_outer_core, &
    nibool_interfaces_outer_core,ibool_interfaces_outer_core

  use specfem_par, only: num_interfaces_crust_mantle, &
    max_nibool_interfaces_cm,my_neighbors_crust_mantle, &
    nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle

  use specfem_par_full_gravity, only: num_interfaces_trinfinite, &
    max_nibool_interfaces_trinfinite,my_neighbors_trinfinite, &
    nibool_interfaces_trinfinite,ibool_interfaces_trinfinite

  use specfem_par_full_gravity, only: num_interfaces_infinite, &
    max_nibool_interfaces_infinite,my_neighbors_infinite, &
    nibool_interfaces_infinite,ibool_interfaces_infinite

  use specfem_par_full_gravity, only: num_interfaces_inner_core1, &
    max_nibool_interfaces_inner_core1,my_neighbors_inner_core1, &
    nibool_interfaces_inner_core1,ibool_interfaces_inner_core1

  use specfem_par_full_gravity, only: num_interfaces_outer_core1, &
    max_nibool_interfaces_outer_core1,my_neighbors_outer_core1, &
    nibool_interfaces_outer_core1,ibool_interfaces_outer_core1

  use specfem_par_full_gravity, only: num_interfaces_crust_mantle1, &
    max_nibool_interfaces_crust_mantle1,my_neighbors_crust_mantle1, &
    nibool_interfaces_crust_mantle1,ibool_interfaces_crust_mantle1

  use specfem_par_full_gravity, only: num_interfaces_trinfinite1, &
    max_nibool_interfaces_trinfinite1,my_neighbors_trinfinite1, &
    nibool_interfaces_trinfinite1,ibool_interfaces_trinfinite1

  use specfem_par_full_gravity, only: num_interfaces_infinite1, &
    max_nibool_interfaces_infinite1,my_neighbors_infinite1, &
    nibool_interfaces_infinite1,ibool_interfaces_infinite1

  use specfem_par_crustmantle, only: ibool_crust_mantle,ibelm_bottom_crust_mantle,ibelm_top_crust_mantle
  use specfem_par_outercore, only: ibool_outer_core,ibelm_bottom_outer_core,ibelm_top_outer_core
  use specfem_par_innercore, only: ibool_inner_core,ibelm_top_inner_core,idoubling_inner_core
  use specfem_par_trinfinite, only: ibool_trinfinite,ibelm_bottom_trinfinite,ibelm_top_trinfinite
  use specfem_par_infinite, only: ibool_infinite,ibelm_bottom_infinite

  use specfem_par_full_gravity, only: is_active_gll,igll_active_on, &
    gdof_cm,gdof_cm1,inode_elmt_cm,inode_elmt_cm1,inode_map_cm, &
    gdof_oc,gdof_oc1,inode_elmt_oc,inode_elmt_oc1,inode_map_oc, &
    gdof_ic,gdof_ic1,inode_elmt_ic,inode_elmt_ic1,inode_map_ic, &
    gdof_trinf,gdof_trinf1,inode_elmt_trinf,inode_elmt_trinf1,inode_map_trinf, &
    gdof_inf,gdof_inf1,inode_elmt_inf,inode_elmt_inf1,inode_map_inf

  use specfem_par_full_gravity, only: neq,neq1,nnode,nnode1, &
    nnode_ic1,nnode_oc1,nnode_cm1,nnode_trinf1,nnode_inf1, &
    nmir_ic,nmir_oc,nmir_cm,nmir_trinf,nmir_inf

  implicit none

  integer :: i,j,k,ier
  integer :: i_elmt,i_node
  integer :: ispec_ic,ispec_oc,ispec_cm,ispec_trinf,ispec_inf
  integer :: ibool_ic,ibool_oc,ibool_cm,ibool_trinf,ibool_inf
  integer :: k_ic,k_oc,k_cm,k_trinf,k_inf
  integer :: ibool,inode,ispec,nnode_icb,nnode_cmb,nnode_trinfb,nnode_infb

  integer,dimension(:),allocatable :: inode_ic,inode_oc,inode_cm,inode_trinf,inode_inf
  integer,dimension(:),allocatable :: inode_ic1,inode_oc1,inode_cm1,inode_trinf1,inode_inf1
  logical,dimension(:),allocatable :: isnode_ic,isnode_oc,isnode_cm,isnode_trinf,isnode_inf
  integer,allocatable :: nf(:,:),nf1(:,:),nmir(:)

  integer :: inode1,inum,igll
  logical,allocatable :: isnode(:)

  logical,allocatable :: isibool_interface_ic(:,:),isibool_interface_oc(:,:), &
                         isibool_interface_cm(:,:),isibool_interface_trinf(:,:),isibool_interface_inf(:,:)

  double precision :: sizeval

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  creating solver array indexing"
    call flush_IMAIN()
  endif

  ! estimated memory size required in MB
  ! inode_elmt_cm,..
  sizeval = dble(NGLLCUBE)*dble(NSPEC_CRUST_MANTLE)*dble(SIZE_INTEGER)
  sizeval = sizeval + dble(NGLLCUBE)*dble(NSPEC_INNER_CORE)*dble(SIZE_INTEGER)
  sizeval = sizeval + dble(NGLLCUBE)*dble(NSPEC_OUTER_CORE)*dble(SIZE_INTEGER)
  sizeval = sizeval + dble(NGLLCUBE)*dble(NSPEC_TRINFINITE)*dble(SIZE_INTEGER)
  sizeval = sizeval + dble(NGLLCUBE)*dble(NSPEC_INFINITE)*dble(SIZE_INTEGER)
  ! inode_elmt_cm1,..
  sizeval = sizeval + dble(NGLLCUBE_INF)*dble(NSPEC_CRUST_MANTLE)*dble(SIZE_INTEGER)
  sizeval = sizeval + dble(NGLLCUBE_INF)*dble(NSPEC_INNER_CORE)*dble(SIZE_INTEGER)
  sizeval = sizeval + dble(NGLLCUBE_INF)*dble(NSPEC_OUTER_CORE)*dble(SIZE_INTEGER)
  sizeval = sizeval + dble(NGLLCUBE_INF)*dble(NSPEC_TRINFINITE)*dble(SIZE_INTEGER)
  sizeval = sizeval + dble(NGLLCUBE_INF)*dble(NSPEC_INFINITE)*dble(SIZE_INTEGER)
  ! inode_map_cm,..
  sizeval = sizeval + 2.d0*dble(NGLOB_CRUST_MANTLE)*dble(SIZE_INTEGER)
  sizeval = sizeval + 2.d0*dble(NGLOB_INNER_CORE)*dble(SIZE_INTEGER)
  sizeval = sizeval + 2.d0*dble(NGLOB_OUTER_CORE)*dble(SIZE_INTEGER)
  sizeval = sizeval + 2.d0*dble(NGLOB_TRINFINITE)*dble(SIZE_INTEGER)
  sizeval = sizeval + 2.d0*dble(NGLOB_INFINITE)*dble(SIZE_INTEGER)
  ! nmir_cm,..
  sizeval = sizeval + dble(NGLOB_CRUST_MANTLE)*dble(SIZE_INTEGER)
  sizeval = sizeval + dble(NGLOB_INNER_CORE)*dble(SIZE_INTEGER)
  sizeval = sizeval + dble(NGLOB_OUTER_CORE)*dble(SIZE_INTEGER)
  sizeval = sizeval + dble(NGLOB_TRINFINITE)*dble(SIZE_INTEGER)
  sizeval = sizeval + dble(NGLOB_INFINITE)*dble(SIZE_INTEGER)

  ! in MB
  sizeval = sizeval / 1024.d0 / 1024.d0

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    size of indexing arrays = ',sngl(sizeval),'MB'
    write(IMAIN,*) '                            = ',sngl(sizeval / 1024.d0),'GB'
    call flush_IMAIN()
  endif

  ! allocate inode arrays
  allocate(inode_elmt_cm(NGLLCUBE,NSPEC_CRUST_MANTLE), &
           inode_elmt_ic(NGLLCUBE,NSPEC_INNER_CORE), &
           inode_elmt_oc(NGLLCUBE,NSPEC_OUTER_CORE), &
           inode_elmt_trinf(NGLLCUBE,NSPEC_TRINFINITE), &
           inode_elmt_inf(NGLLCUBE,NSPEC_INFINITE),stat=ier)
  if (ier /= 0) stop 'Error allocating inode_elmt_cm,.. arrays'
  inode_elmt_cm(:,:) = 0; inode_elmt_ic(:,:) = 0; inode_elmt_oc(:,:) = 0
  inode_elmt_trinf(:,:) = 0; inode_elmt_inf(:,:) = 0

  allocate(inode_elmt_cm1(NGLLCUBE_INF,NSPEC_CRUST_MANTLE), &
           inode_elmt_ic1(NGLLCUBE_INF,NSPEC_INNER_CORE), &
           inode_elmt_oc1(NGLLCUBE_INF,NSPEC_OUTER_CORE), &
           inode_elmt_trinf1(NGLLCUBE_INF,NSPEC_TRINFINITE), &
           inode_elmt_inf1(NGLLCUBE_INF,NSPEC_INFINITE),stat=ier)
  if (ier /= 0) stop 'Error allocating inode_elmt_cm1,.. arrays'
  inode_elmt_cm1(:,:) = 0; inode_elmt_ic1(:,:) = 0; inode_elmt_oc1(:,:) = 0
  inode_elmt_trinf1(:,:) = 0; inode_elmt_inf1(:,:) = 0

  allocate(inode_map_ic(2,NGLOB_INNER_CORE), &
           inode_map_oc(2,NGLOB_OUTER_CORE), &
           inode_map_cm(2,NGLOB_CRUST_MANTLE), &
           inode_map_trinf(2,NGLOB_TRINFINITE), &
           inode_map_inf(2,NGLOB_INFINITE),stat=ier)
  if (ier /= 0) stop 'Error allocating inode_map_ic,.. arrays'
  inode_map_ic(:,:) = 0; inode_map_oc(:,:) = 0; inode_map_cm(:,:) = 0
  inode_map_trinf(:,:) = 0; inode_map_inf(:,:) = 0

  allocate(nmir_ic(NGLOB_INNER_CORE), &
           nmir_oc(NGLOB_OUTER_CORE), &
           nmir_cm(NGLOB_CRUST_MANTLE), &
           nmir_trinf(NGLOB_TRINFINITE), &
           nmir_inf(NGLOB_INFINITE),stat=ier)
  if (ier /= 0) stop 'Error allocating nmir_ic,.. arrays'
  nmir_ic(:) = 0; nmir_oc(:) = 0; nmir_cm(:) = 0
  nmir_trinf(:) = 0; nmir_inf(:) = 0

  ! allocate temporary arrays
  allocate(inode_ic(NGLOB_INNER_CORE), &
           inode_oc(NGLOB_OUTER_CORE), &
           inode_cm(NGLOB_CRUST_MANTLE), &
           inode_trinf(NGLOB_TRINFINITE), &
           inode_inf(NGLOB_INFINITE),stat=ier)
  if (ier /= 0) stop 'Error allocating inode_ic,.. arrays'
  inode_ic(:) = 0; inode_oc(:) = 0; inode_cm(:) = 0;
  inode_trinf(:) = 0; inode_inf(:) = 0

  allocate(isnode_ic(NGLOB_INNER_CORE), &
           isnode_oc(NGLOB_OUTER_CORE), &
           isnode_cm(NGLOB_CRUST_MANTLE), &
           isnode_trinf(NGLOB_TRINFINITE), &
           isnode_inf(NGLOB_INFINITE),stat=ier)
  if (ier /= 0) stop 'Error allocating isnode_ic,.. arrays'
  isnode_ic(:) = .false.; isnode_oc(:) = .false.; isnode_cm(:) = .false.
  isnode_trinf(:) = .false.; isnode_inf(:) = .false.

  ! store ibool in 1D linear mapping
  ! inner core
  do i_elmt = 1,NSPEC_INNER_CORE
    inode_elmt_ic(:,i_elmt) = reshape(ibool_inner_core(:,:,:,i_elmt),(/NGLLCUBE/))
  enddo
  ! outer core
  do i_elmt = 1,NSPEC_OUTER_CORE
    inode_elmt_oc(:,i_elmt) = reshape(ibool_outer_core(:,:,:,i_elmt),(/NGLLCUBE/))
  enddo
  ! inner core
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    inode_elmt_cm(:,i_elmt) = reshape(ibool_crust_mantle(:,:,:,i_elmt),(/NGLLCUBE/))
  enddo
  ! transition infinite
  if (ADD_TRINF) then
  do i_elmt = 1,NSPEC_TRINFINITE
    inode_elmt_trinf(:,i_elmt) = reshape(ibool_trinfinite(:,:,:,i_elmt),(/NGLLCUBE/))
  enddo
  endif
  ! infinite
  do i_elmt = 1,NSPEC_INFINITE
    inode_elmt_inf(:,i_elmt) = reshape(ibool_infinite(:,:,:,i_elmt),(/NGLLCUBE/))
  enddo

  ! count global node numbers
  nnode = NGLOB_INNER_CORE + NGLOB_OUTER_CORE + NGLOB_CRUST_MANTLE + NGLOB_TRINFINITE + NGLOB_INFINITE

  ! identify duplicate nodes on the boundary
  ! inner core - outer core boundary (ICB)
  do i_elmt = 1,NSPEC2D_BOTTOM_OC
    ispec = ibelm_bottom_outer_core(i_elmt)
    k = 1 ! bottom face
    do j = 1,NGLLY
      do i = 1,NGLLX
        isnode_oc(ibool_outer_core(i,j,k,ispec)) = .true.
      enddo
    enddo
  enddo
  nnode_icb = count(isnode_oc)

  ! outer core - crust mantle boundary (CMB)
  do i_elmt = 1,NSPEC2D_BOTTOM_CM
    ispec = ibelm_bottom_crust_mantle(i_elmt)
    k = 1 ! bottom face
    do j = 1,NGLLY
      do i = 1,NGLLX
        isnode_cm(ibool_crust_mantle(i,j,k,ispec)) = .true.
      enddo
    enddo
  enddo
  nnode_cmb = count(isnode_cm)

  ! crust mantle - transition infinite boundary (FS: free surface)
  nnode_trinfb = 0
  if (ADD_TRINF) then
    do i_elmt = 1,NSPEC2D_BOTTOM_TRINF
      ispec = ibelm_bottom_trinfinite(i_elmt)
      k = 1 ! bottom face
      do j = 1,NGLLY
        do i = 1,NGLLX
          isnode_trinf(ibool_trinfinite(i,j,k,ispec)) = .true.
        enddo
      enddo
    enddo
    nnode_trinfb = count(isnode_trinf)
  endif

  ! crust mantle - infinite boundary (FS: free surface)
  do i_elmt = 1,NSPEC2D_BOTTOM_INF
    ispec = ibelm_bottom_infinite(i_elmt)
    k = 1 ! bottom face
    do j = 1,NGLLY
      do i = 1,NGLLX
        isnode_inf(ibool_infinite(i,j,k,ispec)) = .true.
      enddo
    enddo
  enddo
  nnode_infb = count(isnode_inf)

  ! number of unique global nodes
  nnode = nnode - nnode_icb - nnode_cmb - nnode_trinfb - nnode_infb

  ! indexify global nodes and store in a region array
  ! inner core
  inode_ic(:) = (/ (inode, inode = 1,NGLOB_INNER_CORE) /)
  inode = NGLOB_INNER_CORE

  ! outer core
  ! ICB
  ! copy common boundary nodes
  isnode_oc(:) = .false.
  do i_elmt = 1,NSPEC2D_BOTTOM_OC
    ispec_oc = ibelm_bottom_outer_core(i_elmt)
    ispec_ic = ibelm_top_inner_core(i_elmt)
    k_oc = 1;    ! bottom face
    k_ic = NGLLZ ! top face
    do j = 1,NGLLY
      do i = 1,NGLLX
        ibool_oc = ibool_outer_core(i,j,k_oc,ispec_oc)
        ibool_ic = ibool_inner_core(i,j,k_ic,ispec_ic)
        inode_oc(ibool_oc) = inode_ic(ibool_ic)
        isnode_oc(ibool_oc) = .true.
      enddo
    enddo
  enddo
  ! now loop through all nodes
  do i_elmt = 1,NSPEC_OUTER_CORE
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool_oc = ibool_outer_core(i,j,k,i_elmt)
          if (.not. isnode_oc(ibool_oc)) then
            isnode_oc(ibool_oc) = .true.
            inode = inode+1
            inode_oc(ibool_oc) = inode
          endif
        enddo
      enddo
    enddo
  enddo

  ! crust-mantle
  ! CMB
  ! copy common boundary nodes
  isnode_cm(:) = .false.
  do i_elmt = 1,NSPEC2D_BOTTOM_CM
    ispec_cm = ibelm_bottom_crust_mantle(i_elmt)
    ispec_oc = ibelm_top_outer_core(i_elmt)
    k_cm = 1; k_oc = NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        ibool_cm = ibool_crust_mantle(i,j,k_cm,ispec_cm)
        ibool_oc = ibool_outer_core(i,j,k_oc,ispec_oc)
        inode_cm(ibool_cm) = inode_oc(ibool_oc)
        isnode_cm(ibool_cm) = .true.
      enddo
    enddo
  enddo
  ! now loop through all nodes
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool_cm = ibool_crust_mantle(i,j,k,i_elmt)
          if (.not. isnode_cm(ibool_cm)) then
            isnode_cm(ibool_cm) = .true.
            inode = inode+1
            inode_cm(ibool_cm) = inode
          endif
        enddo
      enddo
    enddo
  enddo

  if (ADD_TRINF) then
    ! transition infinite
    ! FS
    ! copy common boundary nodes
    isnode_trinf(:) = .false.
    do i_elmt = 1,NSPEC2D_BOTTOM_TRINF
      ispec_trinf = ibelm_bottom_trinfinite(i_elmt)
      ispec_cm = ibelm_top_crust_mantle(i_elmt)
      k_trinf = 1; k_cm = NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool_trinf = ibool_trinfinite(i,j,k_trinf,ispec_trinf)
          ibool_cm = ibool_crust_mantle(i,j,k_cm,ispec_cm)
          inode_trinf(ibool_trinf) = inode_cm(ibool_cm)
          isnode_trinf(ibool_trinf) = .true.
        enddo
      enddo
    enddo
    ! now loop through all nodes
    do i_elmt = 1,NSPEC_TRINFINITE
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            ibool_trinf = ibool_trinfinite(i,j,k,i_elmt)
            if (.not. isnode_trinf(ibool_trinf)) then
              isnode_trinf(ibool_trinf) = .true.
              inode = inode+1
              inode_trinf(ibool_trinf) = inode
            endif
          enddo
        enddo
      enddo
    enddo

    ! infinite
    ! FS
    ! copy common boundary nodes
    isnode_inf(:) = .false.
    do i_elmt = 1,NSPEC2D_BOTTOM_INF
      ispec_inf = ibelm_bottom_infinite(i_elmt)
      ispec_trinf = ibelm_top_trinfinite(i_elmt)
      k_inf = 1; k_trinf = NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool_inf = ibool_infinite(i,j,k_inf,ispec_inf)
          ibool_trinf = ibool_trinfinite(i,j,k_trinf,ispec_trinf)
          inode_inf(ibool_inf) = inode_trinf(ibool_trinf)
          isnode_inf(ibool_inf) = .true.
        enddo
      enddo
    enddo
    ! now loop through all nodes
    do i_elmt = 1,NSPEC_INFINITE
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            ibool_inf = ibool_infinite(i,j,k,i_elmt)
            if (.not. isnode_inf(ibool_inf)) then
              isnode_inf(ibool_inf) = .true.
              inode = inode+1
              inode_inf(ibool_inf) = inode
            endif
          enddo
        enddo
      enddo
    enddo
  else ! if (ADD_TRINF)
    ! infinite
    ! FS
    ! copy common boundary nodes
    isnode_inf(:) = .false.
    do i_elmt = 1,NSPEC2D_BOTTOM_INF
      ispec_inf = ibelm_bottom_infinite(i_elmt)
      ispec_cm = ibelm_top_crust_mantle(i_elmt)
      k_inf = 1; k_cm = NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool_inf = ibool_infinite(i,j,k_inf,ispec_inf)
          ibool_cm = ibool_crust_mantle(i,j,k_cm,ispec_cm)
          inode_inf(ibool_inf) = inode_cm(ibool_cm)
          isnode_inf(ibool_inf) = .true.
        enddo
      enddo
    enddo
    ! now loop through all nodes
    do i_elmt = 1,NSPEC_INFINITE
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            ibool_inf = ibool_infinite(i,j,k,i_elmt)
            if (.not. isnode_inf(ibool_inf)) then
              isnode_inf(ibool_inf) = .true.
              inode = inode+1
              inode_inf(ibool_inf) = inode
            endif
          enddo
        enddo
      enddo
    enddo
  endif ! if (ADD_TRINF)

  if (inode /= nnode) then
    print *,'ERROR: numbers of global nodes mismatch!',inode,nnode
    call synchronize_all()
    call exit_MPI(myrank,'Invalid number of global nodes')
  endif

  ! global degrees of freedoms
  allocate(nf(1,nnode)) ! only gravitational potential
  nf(:,:) = 0

  ! activate freedoms

  ! freedoms of fictitious cube in inner core are deactivated
  do i_elmt = 1,NSPEC_INNER_CORE
    if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) cycle
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool = ibool_inner_core(i,j,k,i_elmt)
          nf(1,inode_ic(ibool)) = 1
        enddo
      enddo
    enddo
  enddo

  ! outer core
  ! all freedoms are active
  nf(1,inode_oc(:)) = 1

  ! crust-mantle
  ! all freedoms are active
  nf(1,inode_cm(:)) = 1

  ! transition infinite
  if (ADD_TRINF) then
    ! all freedoms are active
    nf(1,inode_trinf(:)) = 1
  endif

  ! infinite element
  ! all but surface nodes are activated
  do i_elmt = 1,NSPEC_INFINITE
    do k = 1,NGLLZ-1
      do j = 1,NGLLY
        do i = 1,NGLLX
          ibool = ibool_infinite(i,j,k,i_elmt)
          nf(1,inode_inf(ibool)) = 1
        enddo
      enddo
    enddo
  enddo

  !! infinite boundary conditions
  !k=NGLLZ !(iface=6)
  !do i_elmt=1,NSPEC_INFINITE
  !  do j=1,NGLLY
  !    do i=1,NGLLX
  !      ibool = ibool_infinite(i,j,k,i_elmt)
  !      nf(1,inode_inf(ibool)) = 0
  !    enddo
  !  enddo
  !enddo

  inode = 0 ! gdof
  do i_node = 1,nnode
    if (nf(1,i_node) > 0) then
      inode = inode+1
      nf(1,i_node) = inode
    endif
  enddo

  ! Level-2 solver
  ! number of active global degrees of freedom
  neq = inode

  allocate(gdof_ic(NGLOB_INNER_CORE), &
           gdof_oc(NGLOB_OUTER_CORE), &
           gdof_cm(NGLOB_CRUST_MANTLE), &
           gdof_trinf(NGLOB_TRINFINITE), &
           gdof_inf(NGLOB_INFINITE),stat=ier)
  if (ier /= 0) stop 'Error allocating gdof_ic,.. arrays'
  gdof_ic(:) = 0; gdof_oc(:) = 0; gdof_cm(:) = 0
  gdof_trinf(:) = 0; gdof_inf(:) = 0

  ! store gdof in a region array
  gdof_ic(:) = nf(1,inode_ic(:))
  gdof_oc(:) = nf(1,inode_oc(:))
  gdof_cm(:) = nf(1,inode_cm(:))
  if (ADD_TRINF) gdof_trinf(:) = nf(1,inode_trinf(:))
  gdof_inf(:) = nf(1,inode_inf(:))

  deallocate(nf)

  !-------------------------------------------------------------------------------

  ! WARNING: ONLY APPLICABLE FOR 5 TO 3 SOLVER
  if (NGLLX /= 5 .or. NGLLY /= 5 .or. NGLLZ /= 5 .or. &
      NGLLX_INF /= 3 .or. NGLLY_INF /= 3 .or. NGLLZ_INF /= 3) then
    print *,'Error: invalid NGLL setting for SIEM indexing'
    print *,'       NGLLX/NGLLY/NGLLZ = ',NGLLX,NGLLY,NGLLZ, ' - all must be equal to 5 for SIEM'
    print *,'       NGLLX_INF/NGLLY_INF/NGLLZ_INF = ',NGLLX_INF,NGLLY_INF,NGLLZ_INF,' - all must be equal to 3 for SIEM'
    stop 'SIEM indexing only works for NGLLX == 5 and NGLLX_INF == 3 setting for now!'
  endif

  ! Level-1 solver---------------
  ! count nodes for 1st level solver
  ! active GLL points
  is_active_gll(:) = .false.
  inum = 0
  do k = 1,NGLLZ,2
    do j = 1,NGLLY,2
      do i = 1,NGLLX,2     ! 1,3,5
        inum = inum+1
        igll = NGLLY*NGLLX*(k-1)+NGLLX*(j-1)+i    ! 1,3,5,(5+1),(5+3),(5+5),..
        is_active_gll(igll) = .true.
        igll_active_on(inum) = igll
      enddo
    enddo
  enddo

  ! active nodes
  ! inner core
  isnode_ic(:) = .false.
  do i_elmt = 1,NSPEC_INNER_CORE
    isnode_ic(inode_elmt_ic(igll_active_on(:),i_elmt)) = .true.
  enddo
  nnode_ic1 = count(isnode_ic)

  ! outer core
  isnode_oc(:) = .false.
  do i_elmt = 1,NSPEC_OUTER_CORE
    isnode_oc(inode_elmt_oc(igll_active_on(:),i_elmt)) = .true.
  enddo
  nnode_oc1 = count(isnode_oc)

  ! crust mantle
  isnode_cm(:) = .false.
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    isnode_cm(inode_elmt_cm(igll_active_on(:),i_elmt)) = .true.
  enddo
  nnode_cm1 = count(isnode_cm)

  ! transition infinite
  if (ADD_TRINF) then
    isnode_trinf(:) = .false.
    do i_elmt = 1,NSPEC_TRINFINITE
      isnode_trinf(inode_elmt_trinf(igll_active_on(:),i_elmt)) = .true.
    enddo
    nnode_trinf1 = count(isnode_trinf)
  endif

  ! infinite
  isnode_inf(:) = .false.
  do i_elmt = 1,NSPEC_INFINITE
    isnode_inf(inode_elmt_inf(igll_active_on(:),i_elmt)) = .true.
  enddo
  nnode_inf1 = count(isnode_inf)

  ! find active nodes and mirror to orginal nodes
  ! inner core
  nmir_ic(:) = 0
  inode = 0
  do i_node = 1,NGLOB_INNER_CORE
    if (isnode_ic(i_node)) then
      inode = inode+1
      nmir_ic(i_node) = inode
    endif
  enddo
  if (inode /= nnode_ic1) then
    print *,inode,nnode_ic1,size(isnode_ic),NGLOB_INNER_CORE
    write(*,'(/,a)')'ERROR: counted level-1 active nodes mismatch in inner core!'
    stop
  endif

  ! outer core
  nmir_oc(:) = 0
  inode = 0
  do i_node = 1,NGLOB_OUTER_CORE
    if (isnode_oc(i_node)) then
      inode = inode+1
      nmir_oc(i_node) = inode
    endif
  enddo
  if (inode /= nnode_oc1) then
    write(*,'(/,a)')'ERROR: counted level-1 active nodes mismatch in outer core!'
    stop
  endif

  ! crust mantle
  nmir_cm(:) = 0
  inode = 0
  do i_node = 1,NGLOB_CRUST_MANTLE
    if (isnode_cm(i_node)) then
      inode = inode+1
      nmir_cm(i_node) = inode
    endif
  enddo
  if (inode /= nnode_cm1) then
    write(*,'(/,a)')'ERROR: counted level-1 active nodes mismatch in crust mantle!'
    stop
  endif

  ! transition infinite
  if (ADD_TRINF) then
    nmir_trinf(:) = 0
    inode = 0
    do i_node = 1,NGLOB_TRINFINITE
      if (isnode_trinf(i_node)) then
        inode = inode+1
        nmir_trinf(i_node) = inode
      endif
    enddo
    if (inode /= nnode_trinf1) then
      write(*,'(/,a)')'ERROR: counted level-1 active nodes mismatch in transition infinite!'
      stop
    endif
  endif

  ! infinite
  nmir_inf(:) = 0
  inode = 0
  do i_node = 1,NGLOB_INFINITE
    if (isnode_inf(i_node)) then
      inode = inode+1
      nmir_inf(i_node) = inode
    endif
  enddo
  if (inode /= nnode_inf1) then
    write(*,'(/,a)')'ERROR: counted level-1 active nodes mismatch in infinite!'
    stop
  endif

  ! store ibool1 in 1D linear mapping
  ! inner core
  do i_elmt = 1,NSPEC_INNER_CORE
    inode_elmt_ic1(:,i_elmt) = nmir_ic(inode_elmt_ic(igll_active_on(:),i_elmt))
  enddo
  call synchronize_all()

  ! outer core
  do i_elmt = 1,NSPEC_OUTER_CORE
    inode_elmt_oc1(:,i_elmt) = nmir_oc(inode_elmt_oc(igll_active_on(:),i_elmt))
  enddo
  call synchronize_all()

  ! crust mantle
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    inode_elmt_cm1(:,i_elmt) = nmir_cm(inode_elmt_cm(igll_active_on(:),i_elmt))
  enddo
  call synchronize_all()

  ! transtion infinite
  if (ADD_TRINF) then
    do i_elmt = 1,NSPEC_TRINFINITE
      inode_elmt_trinf1(:,i_elmt) = nmir_trinf(inode_elmt_trinf(igll_active_on(:),i_elmt))
    enddo
  endif
  call synchronize_all()

  ! infinite
  do i_elmt = 1,NSPEC_INFINITE
    inode_elmt_inf1(:,i_elmt) = nmir_inf(inode_elmt_inf(igll_active_on(:),i_elmt))
  enddo
  call synchronize_all()

  ! counts all nodes
  allocate(isnode(nnode))
  isnode(:) = .false.

  do i_elmt = 1,NSPEC_INNER_CORE
    isnode(inode_ic(inode_elmt_ic(igll_active_on(:),i_elmt))) = .true.
  enddo
  do i_elmt = 1,NSPEC_OUTER_CORE
    isnode(inode_oc(inode_elmt_oc(igll_active_on(:),i_elmt))) = .true.
  enddo
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    isnode(inode_cm(inode_elmt_cm(igll_active_on(:),i_elmt))) = .true.
  enddo
  if (ADD_TRINF) then
    do i_elmt = 1,NSPEC_TRINFINITE
      isnode(inode_trinf(inode_elmt_trinf(igll_active_on(:),i_elmt))) = .true.
    enddo
  endif
  do i_elmt = 1,NSPEC_INFINITE
    isnode(inode_inf(inode_elmt_inf(igll_active_on(:),i_elmt))) = .true.
  enddo

  nnode1 = count(isnode)

  ! node mirror of the entire model
  allocate(nmir(nnode))
  nmir(:) = 0
  inode1 = 0
  do i_node = 1,nnode
    if (isnode(i_node)) then
      inode1 = inode1+1
      nmir(i_node) = inode1
    endif
  enddo
  if (inode1 /= nnode1) then
    print *,inode1,nnode1
    write(*,'(/,a)')'ERROR: counted level-1 active nodes mismatch!'
    stop
  endif

  call synchronize_all()

  ! find inode_ic1,inode_oc1,inode_cm1,inode_inf1
  allocate(inode_ic1(nnode_ic1), &
           inode_oc1(nnode_oc1), &
           inode_cm1(nnode_cm1), &
           inode_trinf1(nnode_trinf1), &
           inode_inf1(nnode_inf1),stat=ier)
  if (ier /= 0) stop 'Error allocating inode_ic1,.. arrays'
  inode_ic1(:) = 0; inode_oc1(:) = 0; inode_cm1(:) = 0
  inode_trinf1(:) = 0; inode_inf1(:) = 0

  ! inner core
  do i_node = 1,NGLOB_INNER_CORE
    if (isnode_ic(i_node)) then
      inode_ic1(nmir_ic(i_node)) = nmir(inode_ic(i_node))
    endif
  enddo
  call synchronize_all()

  ! outer core
  do i_node = 1,NGLOB_OUTER_CORE
    if (isnode_oc(i_node)) then
      inode_oc1(nmir_oc(i_node)) = nmir(inode_oc(i_node))
    endif
  enddo
  call synchronize_all()

  ! crust mantle
  do i_node = 1,NGLOB_CRUST_MANTLE
    if (isnode_cm(i_node)) then
      inode_cm1(nmir_cm(i_node)) = nmir(inode_cm(i_node))
    endif
  enddo
  call synchronize_all()

  ! transition infinite
  if (ADD_TRINF) then
    do i_node = 1,NGLOB_TRINFINITE
      if (isnode_trinf(i_node)) then
        inode_trinf1(nmir_trinf(i_node)) = nmir(inode_trinf(i_node))
      endif
    enddo
  endif
  call synchronize_all()

  ! infinite
  do i_node = 1,NGLOB_INFINITE
    if (isnode_inf(i_node)) then
      inode_inf1(nmir_inf(i_node)) = nmir(inode_inf(i_node))
    endif
  enddo
  call synchronize_all()

  deallocate(nmir)

  ! global degrees of freedoms
  allocate(nf1(1,nnode1)) ! only gravitational potential
  nf1(:,:) = 0

  ! activate freedoms

  ! freedoms of fictitious cube in inner core are deactivated
  do i_elmt = 1,NSPEC_INNER_CORE
    if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) cycle
    do i = 1,NGLLCUBE_INF
      ibool = inode_elmt_ic1(i,i_elmt)
      nf1(1,inode_ic1(ibool)) = 1
    enddo
  enddo
  call synchronize_all()

  ! outer core
  ! all freedoms are active
  nf1(1,inode_oc1(:)) = 1

  ! crust-mantle
  ! all freedoms are active
  nf1(1,inode_cm1(:)) = 1

  ! transition infinite
  if (ADD_TRINF) then
    nf1(1,inode_trinf1(:)) = 1
  endif

  call synchronize_all()

  ! infinite element
  ! surface nodes are deactivated
  do i_elmt = 1,NSPEC_INFINITE
    do k = 1,NGLLZ_INF-1
      do j = 1,NGLLY_INF
        do i = 1,NGLLX_INF
          ! simple sum should also work here
          !ibool = NGLLX_INF*NGLLY_INF*(k-1)+NGLLX_INF*(j-1)+i  !ibool_infinite(i,j,k,i_elmt)
          !nf1(1,inode_inf1(ibool)) = 1

          igll = NGLLX_INF*NGLLY_INF*(k-1) + NGLLX_INF*(j-1) + i
          ibool = inode_elmt_inf1(igll,i_elmt)
          nf1(1,inode_inf1(ibool)) = 1
        enddo
      enddo
    enddo
  enddo
  call synchronize_all()

  inode1 = 0 ! gdof
  do i_node = 1,nnode1
    if (nf1(1,i_node) > 0) then
      inode1 = inode1+1
      nf1(1,i_node) = inode1
    endif
  enddo

  ! Level-1 solver
  ! number of active global degrees of freedom
  neq1 = inode1

  allocate(gdof_ic1(nnode_ic1), &
           gdof_oc1(nnode_oc1), &
           gdof_cm1(nnode_cm1), &
           gdof_trinf1(nnode_trinf1), &
           gdof_inf1(nnode_inf1),stat=ier)
  if (ier /= 0) stop 'Error allocating gdof_ic1,.. arrays'
  gdof_ic1(:) = 0; gdof_oc1(:) = 0; gdof_cm1(:) = 0
  gdof_trinf1(:) = 0; gdof_inf1(:) = 0

  ! store gdof in a region array
  gdof_ic1(:) = nf1(1,inode_ic1(:))
  gdof_oc1(:) = nf1(1,inode_oc1(:))
  gdof_cm1(:) = nf1(1,inode_cm1(:))
  if (ADD_TRINF) gdof_trinf1(:) = nf1(1,inode_trinf1(:))
  gdof_inf1(:) = nf1(1,inode_inf1(:))

  deallocate(nf1)

  !-------------------------------------------------------------------------------

  call synchronize_all()

  ! modify MPI interfaces removing deactivated nodes
  ! NOTE: this has to be done before isnode_ic, isnode_oc, etc. change

  ! inner core
  num_interfaces_inner_core1 = num_interfaces_inner_core
  allocate(my_neighbors_inner_core1(num_interfaces_inner_core1))
  if (num_interfaces_inner_core1 > 0) my_neighbors_inner_core1(:) = my_neighbors_inner_core(:)

  allocate(isibool_interface_ic(max_nibool_interfaces_ic,num_interfaces_inner_core1), &
           nibool_interfaces_inner_core1(num_interfaces_inner_core1))
  isibool_interface_ic(:,:) = .false.
  nibool_interfaces_inner_core1(:) = 0
  do i = 1,num_interfaces_inner_core1
    do j = 1,nibool_interfaces_inner_core(i)
      isibool_interface_ic(j,i) = isnode_ic(ibool_interfaces_inner_core(j,i))
    enddo
    nibool_interfaces_inner_core1(i) = count(isibool_interface_ic(:,i))
  enddo
  if (num_interfaces_inner_core1 > 0) then
    max_nibool_interfaces_inner_core1 = maxval(nibool_interfaces_inner_core1)
  else
    max_nibool_interfaces_inner_core1 = 0
  endif

  allocate(ibool_interfaces_inner_core1(max_nibool_interfaces_inner_core1,num_interfaces_inner_core1))
  ibool_interfaces_inner_core1(:,:) = 0
  do i = 1,num_interfaces_inner_core1
    inum = 0
    do j = 1,nibool_interfaces_inner_core(i)
      if (isibool_interface_ic(j,i)) then
        inum = inum+1
        ibool_interfaces_inner_core1(inum,i) = nmir_ic(ibool_interfaces_inner_core(j,i))
      endif
    enddo
  enddo
  deallocate(isibool_interface_ic)

  ! outer core
  num_interfaces_outer_core1 = num_interfaces_outer_core
  allocate(my_neighbors_outer_core1(num_interfaces_outer_core1))
  if (num_interfaces_outer_core1 > 0) my_neighbors_outer_core1(:) = my_neighbors_outer_core(:)

  allocate(isibool_interface_oc(max_nibool_interfaces_oc,num_interfaces_outer_core1), &
           nibool_interfaces_outer_core1(num_interfaces_outer_core1))
  isibool_interface_oc(:,:) = .false.
  nibool_interfaces_outer_core1(:) = 0
  do i = 1,num_interfaces_outer_core1
    do j = 1,nibool_interfaces_outer_core(i)
      isibool_interface_oc(j,i) = isnode_oc(ibool_interfaces_outer_core(j,i))
    enddo
    nibool_interfaces_outer_core1(i) = count(isibool_interface_oc(:,i))
  enddo
  if (num_interfaces_outer_core1 > 0) then
    max_nibool_interfaces_outer_core1 = maxval(nibool_interfaces_outer_core1)
  else
    max_nibool_interfaces_outer_core1 = 0
  endif

  allocate(ibool_interfaces_outer_core1(max_nibool_interfaces_outer_core1,num_interfaces_outer_core1))
  ibool_interfaces_outer_core1(:,:) = 0
  do i = 1,num_interfaces_outer_core1
    inum = 0
    do j = 1,nibool_interfaces_outer_core(i)
      if (isibool_interface_oc(j,i)) then
        inum = inum+1
        ibool_interfaces_outer_core1(inum,i) = nmir_oc(ibool_interfaces_outer_core(j,i))
      endif
    enddo
  enddo
  deallocate(isibool_interface_oc)

  ! crust mantle
  num_interfaces_crust_mantle1 = num_interfaces_crust_mantle
  allocate(my_neighbors_crust_mantle1(num_interfaces_crust_mantle1))
  if (num_interfaces_crust_mantle1 > 0) my_neighbors_crust_mantle1(:) = my_neighbors_crust_mantle(:)

  allocate(isibool_interface_cm(max_nibool_interfaces_cm,num_interfaces_crust_mantle1), &
           nibool_interfaces_crust_mantle1(num_interfaces_crust_mantle1))
  isibool_interface_cm(:,:) = .false.
  nibool_interfaces_crust_mantle1(:) = 0
  do i = 1,num_interfaces_crust_mantle1
    do j = 1,nibool_interfaces_crust_mantle(i)
      isibool_interface_cm(j,i) = isnode_cm(ibool_interfaces_crust_mantle(j,i))
    enddo
    nibool_interfaces_crust_mantle1(i) = count(isibool_interface_cm(:,i))
  enddo
  if (num_interfaces_crust_mantle1 > 0) then
    max_nibool_interfaces_crust_mantle1 = maxval(nibool_interfaces_crust_mantle1)
  else
    max_nibool_interfaces_crust_mantle1 = 0
  endif

  allocate(ibool_interfaces_crust_mantle1(max_nibool_interfaces_crust_mantle1,num_interfaces_crust_mantle1))
  ibool_interfaces_crust_mantle1(:,:) = 0
  do i = 1,num_interfaces_crust_mantle1
    inum = 0
    do j = 1,nibool_interfaces_crust_mantle(i)
      if (isibool_interface_cm(j,i)) then
        inum = inum+1
        ibool_interfaces_crust_mantle1(inum,i) = nmir_cm(ibool_interfaces_crust_mantle(j,i))
      endif
    enddo
  enddo
  deallocate(isibool_interface_cm)

  ! transition infinite
  if (ADD_TRINF) then
    num_interfaces_trinfinite1 = num_interfaces_trinfinite
    allocate(my_neighbors_trinfinite1(num_interfaces_trinfinite1))
    if (num_interfaces_trinfinite1 > 0) my_neighbors_trinfinite1(:) = my_neighbors_trinfinite(:)

    allocate(isibool_interface_trinf(max_nibool_interfaces_trinfinite,num_interfaces_trinfinite1), &
             nibool_interfaces_trinfinite1(num_interfaces_trinfinite1))
    isibool_interface_trinf(:,:) = .false.
    nibool_interfaces_trinfinite1(:) = 0
    do i = 1,num_interfaces_trinfinite1
      do j = 1,nibool_interfaces_trinfinite(i)
        isibool_interface_trinf(j,i) = isnode_trinf(ibool_interfaces_trinfinite(j,i))
      enddo
      nibool_interfaces_trinfinite1(i) = count(isibool_interface_trinf(:,i))
    enddo
    if (num_interfaces_trinfinite1 > 0) then
      max_nibool_interfaces_trinfinite1 = maxval(nibool_interfaces_trinfinite1)
    else
      max_nibool_interfaces_trinfinite1 = 0
    endif

    allocate(ibool_interfaces_trinfinite1(max_nibool_interfaces_trinfinite1,num_interfaces_trinfinite1))
    ibool_interfaces_trinfinite1(:,:) = 0
    do i = 1,num_interfaces_trinfinite1
      inum = 0
      do j = 1,nibool_interfaces_trinfinite(i)
        if (isibool_interface_trinf(j,i)) then
          inum = inum+1
          ibool_interfaces_trinfinite1(inum,i) = nmir_trinf(ibool_interfaces_trinfinite(j,i))
        endif
      enddo
    enddo
    deallocate(isibool_interface_trinf)
  endif

  ! infinite
  num_interfaces_infinite1 = num_interfaces_infinite
  allocate(my_neighbors_infinite1(num_interfaces_infinite1))
  if (num_interfaces_infinite1 > 0) my_neighbors_infinite1(:) = my_neighbors_infinite(:)

  allocate(isibool_interface_inf(max_nibool_interfaces_infinite,num_interfaces_infinite1), &
           nibool_interfaces_infinite1(num_interfaces_infinite1))
  isibool_interface_inf(:,:) = .false.
  nibool_interfaces_infinite1(:) = 0
  do i = 1,num_interfaces_infinite1
    do j = 1,nibool_interfaces_infinite(i)
      isibool_interface_inf(j,i) = isnode_inf(ibool_interfaces_infinite(j,i))
    enddo
    nibool_interfaces_infinite1(i) = count(isibool_interface_inf(:,i))
  enddo
  if (num_interfaces_infinite1 > 0) then
    max_nibool_interfaces_infinite1 = maxval(nibool_interfaces_infinite1)
  else
    max_nibool_interfaces_infinite1 = 0
  endif

  allocate(ibool_interfaces_infinite1(max_nibool_interfaces_infinite1,num_interfaces_infinite1))
  ibool_interfaces_infinite1(:,:) = 0
  do i = 1,num_interfaces_infinite1
    inum = 0
    do j = 1,nibool_interfaces_infinite(i)
      if (isibool_interface_inf(j,i)) then
        inum = inum+1
        ibool_interfaces_infinite1(inum,i) = nmir_inf(ibool_interfaces_infinite(j,i))
      endif
    enddo
  enddo
  deallocate(isibool_interface_inf)

  call synchronize_all()

  !-------------------------------------------------------------------------------

  ! node map with unique elements required for interpolation
  !
  ! WARNING: isnode_ic, isnode_oc etc. will no longer represent the
  !          active/deactive nodes for Level-1 solver after following segment
  inode_map_ic(:,:) = -1
  isnode_ic(:) = .false.
  do i_elmt = 1,NSPEC_INNER_CORE
    if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) cycle
    do i = 1,NGLLCUBE
      inode = inode_elmt_ic(i,i_elmt)
      if (.not. isnode_ic(inode)) then
        inode_map_ic(1,inode) = i_elmt
        inode_map_ic(2,inode) = i
        isnode_ic(inode) = .true.
      endif
    enddo
  enddo

  call synchronize_all()

  inode_map_oc(:,:) = -1
  isnode_oc(:) = .false.
  do i_elmt = 1,NSPEC_OUTER_CORE
    do i = 1,NGLLCUBE
      inode = inode_elmt_oc(i,i_elmt)
      if (.not. isnode_oc(inode)) then
        inode_map_oc(1,inode) = i_elmt
        inode_map_oc(2,inode) = i
        isnode_oc(inode) = .true.
      endif
    enddo
  enddo

  inode_map_cm(:,:) = -1
  isnode_cm(:) = .false.
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    do i = 1,NGLLCUBE
      inode = inode_elmt_cm(i,i_elmt)
      if (.not. isnode_cm(inode)) then
        inode_map_cm(1,inode) = i_elmt
        inode_map_cm(2,inode) = i
        isnode_cm(inode) = .true.
      endif
    enddo
  enddo

  if (ADD_TRINF) then
    inode_map_trinf(:,:) = -1
    isnode_trinf(:) = .false.
    do i_elmt = 1,NSPEC_TRINFINITE
      do i = 1,NGLLCUBE
        inode = inode_elmt_trinf(i,i_elmt)
        if (.not. isnode_trinf(inode)) then
          inode_map_trinf(1,inode) = i_elmt
          inode_map_trinf(2,inode) = i
          isnode_trinf(inode) = .true.
        endif
      enddo
    enddo
  endif

  inode_map_inf(:,:) = -1
  isnode_inf(:) = .false.
  do i_elmt = 1,NSPEC_INFINITE
    do i = 1,NGLLCUBE
      inode = inode_elmt_inf(i,i_elmt)
      if (.not. isnode_inf(inode)) then
        inode_map_inf(1,inode) = i_elmt
        inode_map_inf(2,inode) = i
        isnode_inf(inode) = .true.
      endif
    enddo
  enddo

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) '    Level-1 solver number of active DOFs: ',neq1
    write(IMAIN,*) '    Level-2 solver number of active DOFs: ',neq
    call flush_IMAIN()
  endif

  call synchronize_all()

  ! free temporary arrays
  deallocate(isnode_ic,isnode_oc,isnode_cm,isnode_trinf,isnode_inf)
  deallocate(inode_ic,inode_oc,inode_cm,inode_trinf,inode_inf)
  deallocate(inode_ic1,inode_oc1,inode_cm1,inode_trinf1,inode_inf1)

  end subroutine SIEM_get_index_region

