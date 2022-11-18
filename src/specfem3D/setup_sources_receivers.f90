!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  8 . 0
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

  subroutine setup_sources_receivers()

  use specfem_par, only: myrank,IMAIN,NSOURCES,NSTEP, &
    theta_source,phi_source, &
    TOPOGRAPHY,ibathy_topo, &
    USE_DISTANCE_CRITERION,xyz_midpoints,xadj,adjncy, &
    SAVE_GREEN_FUNCTIONS, &
    ispec_selected_rec, ispec_cm2gf, hxir_store, hetar_store, hgammar_store

  use kdtree_search, only: kdtree_delete,kdtree_nodes_location,kdtree_nodes_index

  implicit none

  ! setup for point search
  call setup_point_search_arrays()

  ! locates sources and determines simulation start time t0
  call setup_sources()

  ! reads in stations file and locates receivers
  call setup_receivers()

  ! setup green function location arrays
  ! This function resets the KDTree that is used for the source location and
  ! should therefore only ever be used after all other kdtree calls are done
  if (SAVE_GREEN_FUNCTIONS) then
    call setup_green_locations()
  endif

  ! write source and receiver VTK files for Paraview
  call setup_sources_receivers_VTKfile()

  ! pre-compute source arrays
  call setup_sources_precompute_arrays()

  ! pre-compute receiver interpolation factors
  call setup_receivers_precompute_intp()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Total number of samples for seismograms = ',NSTEP
    write(IMAIN,*)
    if (NSOURCES > 1) write(IMAIN,*) 'Using ',NSOURCES,' point sources'
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! frees arrays
  deallocate(theta_source,phi_source)

  ! topography array no more needed
  if (TOPOGRAPHY) then
    if (SAVE_GREEN_FUNCTIONS .eqv. .false.) then
      if (allocated(ibathy_topo) ) deallocate(ibathy_topo)
    endif
  endif

  ! frees memory
  if (USE_DISTANCE_CRITERION) deallocate(xyz_midpoints)
  deallocate(xadj,adjncy)
  ! deletes tree arrays
  deallocate(kdtree_nodes_location)
  deallocate(kdtree_nodes_index)
  ! deletes search tree nodes
  call kdtree_delete()

  end subroutine setup_sources_receivers

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_point_search_arrays()

  use constants, only: &
    NDIM,NGLLX,NGLLY,NGLLZ,MIDX,MIDY,MIDZ,IMAIN,TWO_PI,R_UNIT_SPHERE,DEGREES_TO_RADIANS, &
    USE_DISTANCE_CRITERION

  use specfem_par, only: &
    NCHUNKS_VAL,NEX_XI_VAL,NEX_ETA_VAL,ANGULAR_WIDTH_XI_IN_DEGREES_VAL,ANGULAR_WIDTH_ETA_IN_DEGREES_VAL, &
    LAT_LON_MARGIN,myrank

  use specfem_par, only: &
    nspec => NSPEC_CRUST_MANTLE,nglob => NGLOB_CRUST_MANTLE

  use specfem_par_crustmantle, only: &
    ibool => ibool_crust_mantle, &
    xstore => xstore_crust_mantle,ystore => ystore_crust_mantle,zstore => zstore_crust_mantle

  ! for point search
  use specfem_par, only: &
    element_size,typical_size_squared, &
    anchor_iax,anchor_iay,anchor_iaz, &
    lat_min,lat_max,lon_min,lon_max,xyz_midpoints

  use kdtree_search, only: kdtree_setup,kdtree_set_verbose, &
    kdtree_num_nodes,kdtree_nodes_location,kdtree_nodes_index

  implicit none

  ! local parameters
  integer :: ispec,iglob,ier
  integer :: i,j,k,inodes
  double precision ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD
  ! determines tree points
  logical :: use_midpoints_only

  ! compute typical size of elements at the surface
  ! (normalized)
  if (NCHUNKS_VAL == 6) then
    ! estimation for global meshes (assuming 90-degree chunks)
    element_size = TWO_PI * R_UNIT_SPHERE / (4.d0 * NEX_XI_VAL)
  else
    ! estimation for 1-chunk meshes
    ANGULAR_WIDTH_XI_RAD = ANGULAR_WIDTH_XI_IN_DEGREES_VAL * DEGREES_TO_RADIANS
    ANGULAR_WIDTH_ETA_RAD = ANGULAR_WIDTH_ETA_IN_DEGREES_VAL * DEGREES_TO_RADIANS
    element_size = max( ANGULAR_WIDTH_XI_RAD/NEX_XI_VAL,ANGULAR_WIDTH_ETA_RAD/NEX_ETA_VAL ) * R_UNIT_SPHERE
  endif

  ! use 10 times the distance as a criterion for source detection
  typical_size_squared = (10.d0 * element_size)**2

  ! limits receiver search
  if (USE_DISTANCE_CRITERION) then
    ! retrieves latitude/longitude range of this slice
    call xyz_2_latlon_minmax(nspec,nglob,ibool,xstore,ystore,zstore,lat_min,lat_max,lon_min,lon_max)

    ! adds search margin
    lat_min = lat_min - LAT_LON_MARGIN
    lat_max = lat_max + LAT_LON_MARGIN

    lon_min = lon_min - LAT_LON_MARGIN
    lon_max = lon_max + LAT_LON_MARGIN

    ! limits latitude to [-90.0,90.0]
    if (lat_min < -90.d0 ) lat_min = -90.d0
    if (lat_max > 90.d0 ) lat_max = 90.d0

    ! limits longitude to [0.0,360.0]
    if (lon_min < 0.d0 ) lon_min = 0.d0
    if (lon_min > 360.d0 ) lon_min = 360.d0
    if (lon_max < 0.d0 ) lon_max = 0.d0
    if (lon_max > 360.d0 ) lon_max = 360.d0

    ! prepares midpoints coordinates
    allocate(xyz_midpoints(NDIM,nspec),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating array xyz_midpoints')

    ! store x/y/z coordinates of center point
    do ispec = 1,nspec
      iglob = ibool(MIDX,MIDY,MIDZ,ispec)
      xyz_midpoints(1,ispec) =  dble(xstore(iglob))
      xyz_midpoints(2,ispec) =  dble(ystore(iglob))
      xyz_midpoints(3,ispec) =  dble(zstore(iglob))
    enddo
  endif

  ! define (i,j,k) indices of the control/anchor points
  call hex_nodes_anchor_ijk(anchor_iax,anchor_iay,anchor_iaz)

  ! setups adjacency array to search neighbors
  call setup_adjacency_neighbors()

  ! kd-tree setup for point localization
  !
  ! determines tree size
  if (NEX_ETA_VAL > 100 .and. NEX_XI_VAL > 100) then
    ! high-resolution mesh
    ! only midpoints for search, should be sufficient to get accurate location
    use_midpoints_only = .true.
  else
    ! low-resolution mesh
    ! uses element's inner points
    use_midpoints_only = .false.
  endif

  ! sets total number of tree points
  if (use_midpoints_only) then
    ! small tree size
    kdtree_num_nodes = nspec
  else
    ! uses all internal GLL points for search tree
    ! internal GLL points ( 2 to NGLLX-1 )
    kdtree_num_nodes = nspec * (NGLLX-2)*(NGLLY-2)*(NGLLZ-2)
  endif

  ! allocates tree arrays
  allocate(kdtree_nodes_location(NDIM,kdtree_num_nodes),stat=ier)
  if (ier /= 0) stop 'Error allocating kdtree_nodes_location arrays'
  allocate(kdtree_nodes_index(kdtree_num_nodes),stat=ier)
  if (ier /= 0) stop 'Error allocating kdtree_nodes_index arrays'

  ! prepares search arrays, each element takes its internal GLL points for tree search
  kdtree_nodes_index(:) = 0
  kdtree_nodes_location(:,:) = 0.0
  ! adds tree nodes
  inodes = 0
  if (use_midpoints_only) then
    ! sets up tree nodes
    do ispec = 1,nspec
      iglob = ibool(MIDX,MIDY,MIDZ,ispec)

      ! counts nodes
      inodes = inodes + 1
      if (inodes > kdtree_num_nodes) stop 'Error index inodes bigger than kdtree_num_nodes'

      ! adds node index (index points to same ispec for all internal GLL points)
      kdtree_nodes_index(inodes) = ispec

      ! adds node location
      kdtree_nodes_location(1,inodes) = xstore(iglob)
      kdtree_nodes_location(2,inodes) = ystore(iglob)
      kdtree_nodes_location(3,inodes) = zstore(iglob)
    enddo
  else
    ! all internal GLL points
    do ispec = 1,nspec
      do k = 2,NGLLZ-1
        do j = 2,NGLLY-1
          do i = 2,NGLLX-1
            iglob = ibool(i,j,k,ispec)

            ! counts nodes
            inodes = inodes + 1
            if (inodes > kdtree_num_nodes) stop 'Error index inodes bigger than kdtree_num_nodes'

            ! adds node index (index points to same ispec for all internal GLL points)
            kdtree_nodes_index(inodes) = ispec

            ! adds node location
            kdtree_nodes_location(1,inodes) = xstore(iglob)
            kdtree_nodes_location(2,inodes) = ystore(iglob)
            kdtree_nodes_location(3,inodes) = zstore(iglob)
          enddo
        enddo
      enddo
    enddo
  endif
  if (inodes /= kdtree_num_nodes) stop 'Error index inodes does not match kdtree_num_nodes'

  ! tree verbosity
  if (myrank == 0) call kdtree_set_verbose(IMAIN)

  ! creates kd-tree for searching point locations in locate_point() routine
  call kdtree_setup()

  end subroutine setup_point_search_arrays

!
!-------------------------------------------------------------------------------------------------
!


  subroutine setup_adjacency_neighbors()

  use constants, only: &
    NDIM,NGLLX,NGLLY,NGLLZ,MIDX,MIDY,MIDZ,IMAIN, &
    USE_DISTANCE_CRITERION

  use shared_parameters, only: R_PLANET_KM

  use specfem_par, only: &
    myrank, &
    nspec => NSPEC_CRUST_MANTLE

  use specfem_par_crustmantle, only: &
    ibool => ibool_crust_mantle, &
    xstore => xstore_crust_mantle,ystore => ystore_crust_mantle,zstore => zstore_crust_mantle

  ! for point search
  use specfem_par, only: element_size,typical_size_squared,xyz_midpoints, &
    xadj,adjncy,num_neighbors_all

  use kdtree_search, only: kdtree_setup,kdtree_delete, &
    kdtree_nodes_location,kdtree_nodes_index,kdtree_num_nodes, &
    kdtree_count_nearest_n_neighbors,kdtree_get_nearest_n_neighbors, &
    kdtree_search_index,kdtree_search_num_nodes

  implicit none
  integer :: num_neighbors,num_neighbors_max

  integer,dimension(8) :: iglob_corner,iglob_corner2
  integer :: ispec_ref,ispec,iglob,icorner,ier !,jj

  ! temporary
  integer,parameter :: MAX_NEIGHBORS = 50   ! maximum number of neighbors (around 37 should be sufficient for crust/mantle)
  integer,dimension(:),allocatable :: tmp_adjncy ! temporary adjacency
  integer :: inum_neighbor

  ! timer MPI
  double precision :: time1,tCPU
  double precision, external :: wtime

  ! kd-tree search
  integer :: ielem,inodes
  integer :: nsearch_points
  integer :: num_elements,num_elements_max
  integer :: ielem_counter,num_elements_actual_max
  !integer, parameter :: max_search_points = 2000

  double precision :: r_search
  double precision :: xyz_target(NDIM)
  double precision :: dist_squared,dist_squared_max

  logical :: is_neighbor
  logical,parameter :: DO_BRUTE_FORCE_SEARCH = .false.

  if (myrank == 0) then
    write(IMAIN,*) 'adjacency:'
    write(IMAIN,*) '  total number of elements in this slice = ',nspec
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! get MPI starting time
  time1 = wtime()

  ! adjacency arrays
  !
  ! how to use:
  !  num_neighbors = xadj(ispec+1)-xadj(ispec)
  !  do i = 1,num_neighbors
  !    ! get neighbor
  !    ispec_neighbor = adjncy(xadj(ispec) + i)
  !    ..
  !  enddo
  allocate(xadj(nspec + 1),stat=ier)
  if (ier /= 0) stop 'Error allocating xadj'
  allocate(tmp_adjncy(MAX_NEIGHBORS*nspec),stat=ier)
  if (ier /= 0) stop 'Error allocating tmp_adjncy'
  xadj(:) = 0
  tmp_adjncy(:) = 0

  ! kd-tree search
  if (.not. DO_BRUTE_FORCE_SEARCH) then
    ! kd-tree search

    ! search radius around element midpoints
    !
    ! note: typical search size is using 10 times the size of a surface element;
    !       we take here 6 times the surface element size
    !       - since at low resolutions NEX < 64 and large element sizes, this search radius needs to be large enough, and,
    !       - due to doubling layers (elements at depth will become bigger, however radius shrinks)
    !
    !       the search radius r_search given as routine argument must be non-squared
    r_search = 6.0 * element_size

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  using kd-tree search radius = ',r_search * R_PLANET_KM,'(km)'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! kd-tree setup for adjacency search
    !
    ! uses only element midpoint location
    kdtree_num_nodes = nspec

    ! allocates tree arrays
    allocate(kdtree_nodes_location(NDIM,kdtree_num_nodes),stat=ier)
    if (ier /= 0) stop 'Error allocating kdtree_nodes_location arrays'
    allocate(kdtree_nodes_index(kdtree_num_nodes),stat=ier)
    if (ier /= 0) stop 'Error allocating kdtree_nodes_index arrays'

    ! prepares search arrays, each element takes its internal GLL points for tree search
    kdtree_nodes_index(:) = 0
    kdtree_nodes_location(:,:) = 0.0
    ! adds tree nodes
    inodes = 0
    do ispec = 1,nspec
      ! sets up tree nodes
      iglob = ibool(MIDX,MIDY,MIDZ,ispec)

      ! counts nodes
      inodes = inodes + 1
      if (inodes > kdtree_num_nodes ) stop 'Error index inodes bigger than kdtree_num_nodes'

      ! adds node index (index points to same ispec for all internal GLL points)
      kdtree_nodes_index(inodes) = ispec

      ! adds node location
      kdtree_nodes_location(1,inodes) = xstore(iglob)
      kdtree_nodes_location(2,inodes) = ystore(iglob)
      kdtree_nodes_location(3,inodes) = zstore(iglob)
    enddo
    if (inodes /= kdtree_num_nodes ) stop 'Error index inodes does not match nnodes_local'

    ! alternative: to avoid allocating/deallocating search index arrays, though there is hardly a speedup
    !allocate(kdtree_search_index(max_search_points),stat=ier)
    !if (ier /= 0) stop 'Error allocating array kdtree_search_index'

    ! creates kd-tree for searching
    call kdtree_setup()
  endif

  ! gets maximum number of neighbors
  inum_neighbor = 0
  num_neighbors_max = 0
  num_neighbors_all = 0

  num_elements_max = 0
  num_elements_actual_max = 0
  dist_squared_max = 0.d0

  do ispec_ref = 1,nspec
    ! the eight corners of the current element
    iglob_corner(1) = ibool(1,1,1,ispec_ref)
    iglob_corner(2) = ibool(NGLLX,1,1,ispec_ref)
    iglob_corner(3) = ibool(NGLLX,NGLLY,1,ispec_ref)
    iglob_corner(4) = ibool(1,NGLLY,1,ispec_ref)
    iglob_corner(5) = ibool(1,1,NGLLZ,ispec_ref)
    iglob_corner(6) = ibool(NGLLX,1,NGLLZ,ispec_ref)
    iglob_corner(7) = ibool(NGLLX,NGLLY,NGLLZ,ispec_ref)
    iglob_corner(8) = ibool(1,NGLLY,NGLLZ,ispec_ref)

    ! midpoint for search radius
    iglob = ibool(MIDX,MIDY,MIDZ,ispec_ref)
    xyz_target(1) = xstore(iglob)
    xyz_target(2) = ystore(iglob)
    xyz_target(3) = zstore(iglob)

    if (DO_BRUTE_FORCE_SEARCH) then
      ! loops over all other elements to find closest neighbors
      num_elements = nspec
    else
      ! looks only at elements in kd-tree search radius

      ! gets number of tree points within search radius
      ! (within search sphere)
      call kdtree_count_nearest_n_neighbors(xyz_target,r_search,nsearch_points)

      ! debug
      !print *,'  total number of search elements: ',nsearch_points,'ispec',ispec_ref

      ! alternative: limits search results
      !if (nsearch_points > max_search_points) nsearch_points = max_search_points

      ! sets number of search nodes to get
      kdtree_search_num_nodes = nsearch_points

      ! allocates search index
      allocate(kdtree_search_index(kdtree_search_num_nodes),stat=ier)
      if (ier /= 0) stop 'Error allocating array kdtree_search_index'

      ! gets closest n points around target (within search sphere)
      call kdtree_get_nearest_n_neighbors(xyz_target,r_search,nsearch_points)

      ! loops over search radius
      num_elements = nsearch_points
    endif

    ! statistics
    if (num_elements > num_elements_max) num_elements_max = num_elements
    ielem_counter = 0

    ! counts number of neighbors
    num_neighbors = 0

    do ielem = 1,num_elements

      ! gets element index
      if (DO_BRUTE_FORCE_SEARCH) then
        ispec = ielem
      else
        ! kd-tree search radius
        ! gets search point/element index
        ispec = kdtree_search_index(ielem)
        ! checks index
        if (ispec < 1 .or. ispec > nspec) stop 'Error element index is invalid'
      endif

      ! skip reference element
      if (ispec == ispec_ref) cycle

      ! distance to reference element
      dist_squared = (xyz_target(1) - xyz_midpoints(1,ispec))*(xyz_target(1) - xyz_midpoints(1,ispec)) &
                   + (xyz_target(2) - xyz_midpoints(2,ispec))*(xyz_target(2) - xyz_midpoints(2,ispec)) &
                   + (xyz_target(3) - xyz_midpoints(3,ispec))*(xyz_target(3) - xyz_midpoints(3,ispec))

      ! exclude elements that are too far from target
      if (USE_DISTANCE_CRITERION) then
        !  we compare squared distances instead of distances themselves to significantly speed up calculations
        if (dist_squared > typical_size_squared) cycle
      endif

      ielem_counter = ielem_counter + 1

      ! checks if element has a corner iglob from reference element
      is_neighbor = .false.

      iglob_corner2(1) = ibool(1,1,1,ispec)
      iglob_corner2(2) = ibool(NGLLX,1,1,ispec)
      iglob_corner2(3) = ibool(NGLLX,NGLLY,1,ispec)
      iglob_corner2(4) = ibool(1,NGLLY,1,ispec)
      iglob_corner2(5) = ibool(1,1,NGLLZ,ispec)
      iglob_corner2(6) = ibool(NGLLX,1,NGLLZ,ispec)
      iglob_corner2(7) = ibool(NGLLX,NGLLY,NGLLZ,ispec)
      iglob_corner2(8) = ibool(1,NGLLY,NGLLZ,ispec)

      do icorner = 1,8
        ! checks if corner also has reference element
        if (any(iglob_corner(:) == iglob_corner2(icorner))) then
          is_neighbor = .true.
          exit
        endif
        ! alternative: (slightly slower with 12.4s compared to 11.4s with any() intrinsic function)
        !do jj = 1,8
        !  if (iglob == iglob_corner(jj)) then
        !    is_neighbor = .true.
        !    exit
        !  endif
        !enddo
        !if (is_neighbor) exit
      enddo

      ! counts neighbors to reference element
      if (is_neighbor) then
        ! adds to adjacency
        inum_neighbor = inum_neighbor + 1
        ! checks
        if (inum_neighbor > MAX_NEIGHBORS*nspec) stop 'Error maximum neighbors exceeded'
        ! adds element
        tmp_adjncy(inum_neighbor) = ispec

        ! for statistics
        num_neighbors = num_neighbors + 1

        ! maximum distance to reference element
        if (dist_squared > dist_squared_max) dist_squared_max = dist_squared
      endif
    enddo ! ielem

    ! checks if neighbors were found
    if (num_neighbors == 0) then
      print *,'Error: rank ',myrank,' - element ',ispec_ref,'has no neighbors!'
      print *,'  element midpoint location: ',xyz_target(:)*R_PLANET_KM
      print *,'  typical element size     : ',element_size*R_PLANET_KM,'(km)'
      print *,'  brute force search       : ',DO_BRUTE_FORCE_SEARCH
      print *,'  distance criteria        : ',USE_DISTANCE_CRITERION
      print *,'  typical search distance  : ',typical_size_squared*R_PLANET_KM,'(km)'
      print *,'  kd-tree r_search         : ',r_search*R_PLANET_KM,'(km)'
      print *,'  search elements          : ',num_elements
      call exit_MPI(myrank,'Error adjacency invalid')
    endif

    ! statistics
    if (num_neighbors > num_neighbors_max) num_neighbors_max = num_neighbors
    if (ielem_counter > num_elements_actual_max) num_elements_actual_max = ielem_counter

    ! adjacency indexing
    xadj(ispec_ref + 1) = inum_neighbor
    ! how to use:
    !num_neighbors = xadj(ispec+1)-xadj(ispec)
    !do i = 1,num_neighbors
    !  ! get neighbor
    !  ispec_neighbor = adjncy(xadj(ispec) + i)
    !enddo

    ! frees kdtree search array
    if (.not. DO_BRUTE_FORCE_SEARCH) then
      deallocate(kdtree_search_index)
    endif

  enddo ! ispec_ref

  ! total number of neighbors
  num_neighbors_all = inum_neighbor

  ! allocates compacted array
  allocate(adjncy(num_neighbors_all),stat=ier)
  if (ier /= 0) stop 'Error allocating tmp_adjncy'

  adjncy(1:num_neighbors_all) = tmp_adjncy(1:num_neighbors_all)

  ! checks
  if (minval(adjncy(:)) < 1 .or. maxval(adjncy(:)) > nspec) stop 'Invalid adjncy array'

  ! frees temporary array
  deallocate(tmp_adjncy)

  if (.not. DO_BRUTE_FORCE_SEARCH) then
    ! frees current tree memory
    ! deletes tree arrays
    deallocate(kdtree_nodes_location)
    deallocate(kdtree_nodes_index)
    ! deletes search tree nodes
    call kdtree_delete()
  endif

  if (myrank == 0) then
    ! elapsed time since beginning of neighbor detection
    tCPU = wtime() - time1
    write(IMAIN,*) '  maximum search elements                                      = ',num_elements_max
    write(IMAIN,*) '  maximum of actual search elements (after distance criterion) = ',num_elements_actual_max
    write(IMAIN,*)
    write(IMAIN,*) '  estimated typical element size at surface = ',element_size*R_PLANET_KM,'(km)'
    write(IMAIN,*) '  maximum distance between neighbor centers = ',sqrt(dist_squared_max)*R_PLANET_KM,'(km)'
    if (.not. DO_BRUTE_FORCE_SEARCH) then
      if (sqrt(dist_squared_max) > r_search - 0.5*element_size) then
          write(IMAIN,*) '***'
          write(IMAIN,*) '*** Warning: consider increasing the kd-tree search radius to improve this neighbor setup ***'
          write(IMAIN,*) '***'
      endif
    endif
    write(IMAIN,*)
    write(IMAIN,*) '  maximum neighbors found per element = ',num_neighbors_max,'(should be 37 for globe meshes)'
    write(IMAIN,*) '  total number of neighbors           = ',num_neighbors_all
    write(IMAIN,*)
    write(IMAIN,*) '  Elapsed time for detection of neighbors in seconds = ',tCPU
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine setup_adjacency_neighbors

! -------------------------------------------------------------------------

! subroutine setup_adjacency_neighbors_gf()
! ! The routine `setup_adjacency_neighbors_gf()` is used to compute the
! ! neighbours of the subselected elements, so that source location can be done
! ! using the neighbors of the spectral elements as well.

!   use constants, only: &
!     NDIM,NGLLX,NGLLY,NGLLZ,MIDX,MIDY,MIDZ,IMAIN, &
!     USE_DISTANCE_CRITERION, CUSTOM_REAL

!   use shared_parameters, only: R_PLANET_KM

!   use specfem_par, only: &
!     myrank, &
!     iglob_cm2gf, &
!     nspec => ngf_unique_local, &
!     nglob => NGLOB_GF, &
!     ibool => ibool_GF


!   use specfem_par_crustmantle, only: xstore_crust_mantle, &
!     ystore_crust_mantle, zstore_crust_mantle

!   ! for point search
!   use specfem_par, only: element_size,typical_size_squared, &
!     xadj_gf,adjncy_gf,num_neighbors_all_gf

!   use kdtree_search, only: kdtree_setup,kdtree_delete, &
!     kdtree_nodes_location,kdtree_nodes_index,kdtree_num_nodes, &
!     kdtree_count_nearest_n_neighbors,kdtree_get_nearest_n_neighbors, &
!     kdtree_search_index,kdtree_search_num_nodes

!   implicit none
!   integer :: num_neighbors,num_neighbors_max

!   integer,dimension(8) :: iglob_corner,iglob_corner2
!   integer :: ispec_ref,ispec,iglob,icorner,ier !,jj

!   real(kind=CUSTOM_REAL), dimension(nglob) :: xstore
!   real(kind=CUSTOM_REAL), dimension(nglob) :: ystore
!   real(kind=CUSTOM_REAL), dimension(nglob) :: zstore
!   real(kind=CUSTOM_REAL), dimension(nglob,3) :: xyz_midpoints

!   ! temporary
!   integer,parameter :: MAX_NEIGHBORS = 50   ! maximum number of neighbors (around 37 should be sufficient for crust/mantle)
!   integer,dimension(:),allocatable :: tmp_adjncy ! temporary adjacency
!   integer :: inum_neighbor

!   ! timer MPI
!   double precision :: time1,tCPU
!   double precision, external :: wtime

!   ! kd-tree search
!   integer :: ielem,inodes
!   integer :: nsearch_points
!   integer :: num_elements,num_elements_max
!   integer :: ielem_counter,num_elements_actual_max
!   !integer, parameter :: max_search_points = 2000

!   double precision :: r_search
!   double precision :: xyz_target(NDIM)
!   double precision :: dist_squared,dist_squared_max

!   logical :: is_neighbor
!   logical,parameter :: DO_BRUTE_FORCE_SEARCH = .false.

!   if (myrank == 0) then
!     write(IMAIN,*) 'SGT Element adjacency:'
!     write(IMAIN,*) '  total number of tagged elements in this slice = ',nspec
!     write(IMAIN,*)
!     call flush_IMAIN()
!   endif

!   ! get MPI starting time
!   time1 = wtime()

!   ! Get local GF global coordinates
!   write(*,*) 'rank', myrank, 'xstore', xstore_crust_mantle(1)
!   write(*,*) 'rank', myrank, 'iglob', iglob_cm2gf(1)
!   xstore = xstore_crust_mantle(iglob_cm2gf)
!   ystore = ystore_crust_mantle(iglob_cm2gf)
!   zstore = zstore_crust_mantle(iglob_cm2gf)

!   ! Get local GF element midpoints
!   xyz_midpoints(:, 1) = xstore(ibool(MIDX,MIDY,MIDZ,:))
!   xyz_midpoints(:, 2) = ystore(ibool(MIDX,MIDY,MIDZ,:))
!   xyz_midpoints(:, 3) = zstore(ibool(MIDX,MIDY,MIDZ,:))
!   ! adjacency arrays
!   !
!   ! how to use:
!   !  num_neighbors = xadj(ispec+1)-xadj(ispec)
!   !  do i = 1,num_neighbors
!   !    ! get neighbor
!   !    ispec_neighbor = adjncy(xadj(ispec) + i)
!   !    ..
!   !  enddo
!   allocate(xadj_gf(nspec + 1),stat=ier)
!   if (ier /= 0) stop 'Error allocating xadj'
!   allocate(tmp_adjncy(MAX_NEIGHBORS*nspec),stat=ier)
!   if (ier /= 0) stop 'Error allocating tmp_adjncy'
!   xadj_gf(:) = 0
!   tmp_adjncy(:) = 0

!   ! kd-tree search
!   if (.not. DO_BRUTE_FORCE_SEARCH) then
!     ! kd-tree search

!     ! search radius around element midpoints
!     !
!     ! note: typical search size is using 10 times the size of a surface element;
!     !       we take here 6 times the surface element size
!     !       - since at low resolutions NEX < 64 and large element sizes, this search radius needs to be large enough, and,
!     !       - due to doubling layers (elements at depth will become bigger, however radius shrinks)
!     !
!     !       the search radius r_search given as routine argument must be non-squared
!     r_search = 6.0 * element_size

!     ! user output
!     if (myrank == 0) then
!       write(IMAIN,*) '  using kd-tree search radius = ',r_search * R_PLANET_KM,'(km)'
!       write(IMAIN,*)
!       call flush_IMAIN()
!     endif

!     ! kd-tree setup for adjacency search
!     !
!     ! uses only element midpoint location
!     kdtree_num_nodes = nspec

!     write (*,*) 'rank', myrank, 'kdtree_num_nodes', kdtree_num_nodes
!     write (*,*) 'rank', myrank, 'nspec', nspec

!     write (*,*) 'rank', myrank, 'node_locash', kdtree_nodes_location(1,1)
!     ! allocates tree arrays
!     allocate(kdtree_nodes_location(NDIM,kdtree_num_nodes),stat=ier)
!     if (ier /= 0) stop 'Error allocating kdtree_nodes_location arrays'
!     allocate(kdtree_nodes_index(kdtree_num_nodes),stat=ier)
!     if (ier /= 0) stop 'Error allocating kdtree_nodes_index arrays'

!     ! prepares search arrays, each element takes its internal GLL points for tree search
!     kdtree_nodes_index(:) = 0
!     kdtree_nodes_location(:,:) = 0.0
!     ! adds tree nodes
!     inodes = 0
!     do ispec = 1,nspec
!       ! sets up tree nodes
!       iglob = ibool(MIDX,MIDY,MIDZ,ispec)

!       ! counts nodes
!       inodes = inodes + 1
!       if (inodes > kdtree_num_nodes ) stop 'Error index inodes bigger than kdtree_num_nodes'

!       ! adds node index (index points to same ispec for all internal GLL points)
!       kdtree_nodes_index(inodes) = ispec

!       ! adds node location
!       kdtree_nodes_location(1,inodes) = xstore(iglob)
!       kdtree_nodes_location(2,inodes) = ystore(iglob)
!       kdtree_nodes_location(3,inodes) = zstore(iglob)
!     enddo
!     if (inodes /= kdtree_num_nodes ) stop 'Error index inodes does not match nnodes_local'

!     ! alternative: to avoid allocating/deallocating search index arrays, though there is hardly a speedup
!     !allocate(kdtree_search_index(max_search_points),stat=ier)
!     !if (ier /= 0) stop 'Error allocating array kdtree_search_index'

!     ! creates kd-tree for searching
!     call kdtree_setup()
!   endif

!   ! gets maximum number of neighbors
!   inum_neighbor = 0
!   num_neighbors_max = 0
!   num_neighbors_all_gf = 0

!   num_elements_max = 0
!   num_elements_actual_max = 0
!   dist_squared_max = 0.d0

!   do ispec_ref = 1,nspec

!     ! the eight corners of the current element
!     iglob_corner(1) = ibool(1,1,1,ispec_ref)
!     iglob_corner(2) = ibool(NGLLX,1,1,ispec_ref)
!     iglob_corner(3) = ibool(NGLLX,NGLLY,1,ispec_ref)
!     iglob_corner(4) = ibool(1,NGLLY,1,ispec_ref)
!     iglob_corner(5) = ibool(1,1,NGLLZ,ispec_ref)
!     iglob_corner(6) = ibool(NGLLX,1,NGLLZ,ispec_ref)
!     iglob_corner(7) = ibool(NGLLX,NGLLY,NGLLZ,ispec_ref)
!     iglob_corner(8) = ibool(1,NGLLY,NGLLZ,ispec_ref)

!     ! midpoint for search radius
!     iglob = ibool(MIDX,MIDY,MIDZ,ispec_ref)
!     xyz_target(1) = xstore(iglob)
!     xyz_target(2) = ystore(iglob)
!     xyz_target(3) = zstore(iglob)

!     if (DO_BRUTE_FORCE_SEARCH) then
!       ! loops over all other elements to find closest neighbors
!       num_elements = nspec
!     else
!       ! looks only at elements in kd-tree search radius

!       ! gets number of tree points within search radius
!       ! (within search sphere)
!       call kdtree_count_nearest_n_neighbors(xyz_target,r_search,nsearch_points)

!       ! debug
!       !print *,'  total number of search elements: ',nsearch_points,'ispec',ispec_ref

!       ! alternative: limits search results
!       !if (nsearch_points > max_search_points) nsearch_points = max_search_points

!       ! sets number of search nodes to get
!       kdtree_search_num_nodes = nsearch_points

!       ! allocates search index
!       allocate(kdtree_search_index(kdtree_search_num_nodes),stat=ier)
!       if (ier /= 0) stop 'Error allocating array kdtree_search_index'

!       ! gets closest n points around target (within search sphere)
!       call kdtree_get_nearest_n_neighbors(xyz_target,r_search,nsearch_points)

!       ! loops over search radius
!       num_elements = nsearch_points
!     endif

!     ! statistics
!     if (num_elements > num_elements_max) num_elements_max = num_elements
!     ielem_counter = 0

!     ! counts number of neighbors
!     num_neighbors = 0

!     do ielem = 1,num_elements

!       ! gets element index
!       if (DO_BRUTE_FORCE_SEARCH) then
!         ispec = ielem
!       else
!         ! kd-tree search radius
!         ! gets search point/element index
!         ispec = kdtree_search_index(ielem)
!         ! checks index
!         if (ispec < 1 .or. ispec > nspec) stop 'Error element index is invalid'
!       endif

!       ! skip reference element
!       if (ispec == ispec_ref) cycle

!       ! distance to reference element
!       dist_squared = (xyz_target(1) - xyz_midpoints(1,ispec))*(xyz_target(1) - xyz_midpoints(1,ispec)) &
!                    + (xyz_target(2) - xyz_midpoints(2,ispec))*(xyz_target(2) - xyz_midpoints(2,ispec)) &
!                    + (xyz_target(3) - xyz_midpoints(3,ispec))*(xyz_target(3) - xyz_midpoints(3,ispec))

!       ! exclude elements that are too far from target
!       if (USE_DISTANCE_CRITERION) then
!         !  we compare squared distances instead of distances themselves to significantly speed up calculations
!         if (dist_squared > typical_size_squared) cycle
!       endif

!       ielem_counter = ielem_counter + 1

!       ! checks if element has a corner iglob from reference element
!       is_neighbor = .false.

!       iglob_corner2(1) = ibool(1,1,1,ispec)
!       iglob_corner2(2) = ibool(NGLLX,1,1,ispec)
!       iglob_corner2(3) = ibool(NGLLX,NGLLY,1,ispec)
!       iglob_corner2(4) = ibool(1,NGLLY,1,ispec)
!       iglob_corner2(5) = ibool(1,1,NGLLZ,ispec)
!       iglob_corner2(6) = ibool(NGLLX,1,NGLLZ,ispec)
!       iglob_corner2(7) = ibool(NGLLX,NGLLY,NGLLZ,ispec)
!       iglob_corner2(8) = ibool(1,NGLLY,NGLLZ,ispec)

!       do icorner = 1,8
!         ! checks if corner also has reference element
!         if (any(iglob_corner(:) == iglob_corner2(icorner))) then
!           is_neighbor = .true.
!           exit
!         endif
!         ! alternative: (slightly slower with 12.4s compared to 11.4s with any() intrinsic function)
!         !do jj = 1,8
!         !  if (iglob == iglob_corner(jj)) then
!         !    is_neighbor = .true.
!         !    exit
!         !  endif
!         !enddo
!         !if (is_neighbor) exit
!       enddo

!       ! counts neighbors to reference element
!       if (is_neighbor) then
!         ! adds to adjacency
!         inum_neighbor = inum_neighbor + 1
!         ! checks
!         if (inum_neighbor > MAX_NEIGHBORS*nspec) stop 'Error maximum neighbors exceeded'
!         ! adds element
!         tmp_adjncy(inum_neighbor) = ispec

!         ! for statistics
!         num_neighbors = num_neighbors + 1

!         ! maximum distance to reference element
!         if (dist_squared > dist_squared_max) dist_squared_max = dist_squared
!       endif
!     enddo ! ielem

!     ! checks if neighbors were found
!     if (num_neighbors == 0) then
!       print *, 'element', ispec_ref, 'has not neighbors'
!       print *,'This is ok if there is no buffer element, but that has to be check'
!       ! print *,'Error: rank ',myrank,' - element ',ispec_ref,'has no neighbors!'
!       ! print *,'  element midpoint location: ',xyz_target(:)*R_PLANET_KM
!       ! print *,'  typical element size     : ',element_size*R_PLANET_KM,'(km)'
!       ! print *,'  brute force search       : ',DO_BRUTE_FORCE_SEARCH
!       ! print *,'  distance criteria        : ',USE_DISTANCE_CRITERION
!       ! print *,'  typical search distance  : ',typical_size_squared*R_PLANET_KM,'(km)'
!       ! print *,'  kd-tree r_search         : ',r_search*R_PLANET_KM,'(km)'
!       ! print *,'  search elements          : ',num_elements
!       ! call exit_MPI(myrank,'Error adjacency invalid')
!       cycle
!     endif

!     ! statistics
!     if (num_neighbors > num_neighbors_max) num_neighbors_max = num_neighbors
!     if (ielem_counter > num_elements_actual_max) num_elements_actual_max = ielem_counter

!     ! adjacency indexing
!     xadj_gf(ispec_ref + 1) = inum_neighbor
!     ! how to use:
!     !num_neighbors = xadj(ispec+1)-xadj(ispec)
!     !do i = 1,num_neighbors
!     !  ! get neighbor
!     !  ispec_neighbor = adjncy(xadj(ispec) + i)
!     !enddo

!     ! frees kdtree search array
!     if (.not. DO_BRUTE_FORCE_SEARCH) then
!       deallocate(kdtree_search_index)
!     endif

!   enddo ! ispec_ref

!   ! total number of neighbors
!   num_neighbors_all_gf = inum_neighbor

!   ! allocates compacted array
!   if (num_neighbors_all_gf > 0) then
!   !   cycle
!     allocate(adjncy_gf(num_neighbors_all_gf),stat=ier)
!     if (ier /= 0) stop 'Error allocating tmp_adjncy'

!     adjncy_gf(1:num_neighbors_all_gf) = tmp_adjncy(1:num_neighbors_all_gf)

!     ! checks
!     if (minval(adjncy_gf(:)) < 1 .or. maxval(adjncy_gf(:)) > nspec) stop 'Invalid adjncy array'
!   endif

!   ! frees temporary array
!   deallocate(tmp_adjncy)

!   if (.not. DO_BRUTE_FORCE_SEARCH) then
!     ! frees current tree memory
!     ! deletes tree arrays
!     deallocate(kdtree_nodes_location)
!     deallocate(kdtree_nodes_index)
!     ! deletes search tree nodes
!     call kdtree_delete()
!   endif

!   if (myrank == 0) then
!     ! elapsed time since beginning of neighbor detection
!     tCPU = wtime() - time1
!     write(IMAIN,*) '  maximum search elements                                      = ',num_elements_max
!     write(IMAIN,*) '  maximum of actual search elements (after distance criterion) = ',num_elements_actual_max
!     write(IMAIN,*)
!     write(IMAIN,*) '  estimated typical element size at surface = ',element_size*R_PLANET_KM,'(km)'
!     write(IMAIN,*) '  maximum distance between neighbor centers = ',sqrt(dist_squared_max)*R_PLANET_KM,'(km)'
!     if (.not. DO_BRUTE_FORCE_SEARCH) then
!       if (sqrt(dist_squared_max) > r_search - 0.5*element_size) then
!           write(IMAIN,*) '***'
!           write(IMAIN,*) '*** Warning: consider increasing the kd-tree search radius to improve this neighbor setup ***'
!           write(IMAIN,*) '***'
!       endif
!     endif
!     write(IMAIN,*)
!     write(IMAIN,*) '  maximum neighbors found per element = ',num_neighbors_max,'(should be 37 for globe meshes)'
!     write(IMAIN,*) '  total number of neighbors           = ',num_neighbors_all_gf
!     write(IMAIN,*)
!     write(IMAIN,*) '  Elapsed time for detection of neighbors in seconds = ',tCPU
!     write(IMAIN,*)
!     call flush_IMAIN()
!   endif

!   end subroutine setup_adjacency_neighbors_gf


!
!-------------------------------------------------------------------------------------------------
!


  subroutine setup_sources()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_movie
  implicit none

  ! local parameters
  integer :: isource,ier
  character(len=MAX_STRING_LEN) :: filename

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'sources:',NSOURCES
    call flush_IMAIN()
  endif

  ! allocate arrays for source
  allocate(islice_selected_source(NSOURCES), &
           ispec_selected_source(NSOURCES), &
           Mxx(NSOURCES), &
           Myy(NSOURCES), &
           Mzz(NSOURCES), &
           Mxy(NSOURCES), &
           Mxz(NSOURCES), &
           Myz(NSOURCES),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating source arrays')
  ! initializes arrays
  islice_selected_source(:) = -1
  ispec_selected_source(:) = 0
  Mxx(:) = 0.d0; Myy(:) = 0.d0; Mzz(:) = 0.d0
  Mxy(:) = 0.d0; Mxz(:) = 0.d0; Myz(:) = 0.d0

  allocate(xi_source(NSOURCES), &
           eta_source(NSOURCES), &
           gamma_source(NSOURCES),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating source arrays')
  xi_source(:) = 0.d0; eta_source(:) = 0.d0; gamma_source(:) = 0.d0

  allocate(tshift_src(NSOURCES), &
           hdur(NSOURCES), &
           hdur_Gaussian(NSOURCES),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating source arrays')
  tshift_src(:) = 0.d0; hdur(:) = 0.d0; hdur_Gaussian(:) = 0.d0

  allocate(theta_source(NSOURCES), &
           phi_source(NSOURCES),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating source arrays')
  theta_source(:) = 0.d0; phi_source(:) = 0.d0

  allocate(nu_source(NDIM,NDIM,NSOURCES),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating source arrays')
  nu_source(:,:,:) = 0.d0

  if (USE_FORCE_POINT_SOURCE) then
    allocate(force_stf(NSOURCES),factor_force_source(NSOURCES), &
             comp_dir_vect_source_E(NSOURCES), &
             comp_dir_vect_source_N(NSOURCES), &
             comp_dir_vect_source_Z_UP(NSOURCES),stat=ier)
    if (ier /= 0) stop 'error allocating arrays for force point sources'
    force_stf(:) = 0
    factor_force_source(:) = 0.d0
    comp_dir_vect_source_E(:) = 0.d0
    comp_dir_vect_source_N(:) = 0.d0
    comp_dir_vect_source_Z_UP(:) = 0.d0
  endif

  ! sources
  ! BS BS moved open statement and writing of first lines into sr.vtk before the
  ! call to locate_sources, where further write statements to that file follow
  if (myrank == 0) then
  ! write source and receiver VTK files for Paraview
    filename = trim(OUTPUT_FILES)//'/sr_tmp.vtk'
    open(IOUT_VTK,file=trim(filename),status='unknown',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening temporary file sr_temp.vtk')
    write(IOUT_VTK,'(a)') '# vtk DataFile Version 2.0'
    write(IOUT_VTK,'(a)') 'Source and Receiver VTK file'
    write(IOUT_VTK,'(a)') 'ASCII'
    write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
    !  LQY -- won't be able to know NSOURCES+nrec at this point...
    write(IOUT_VTK, '(a,i6,a)') 'POINTS ', NSOURCES, ' float'
    ! closing file, rest of information will be appended later on
    close(IOUT_VTK)
  endif

  ! locate sources in the mesh
  call locate_sources()

  ! determines onset time
  call setup_stf_constants()

  ! count number of sources located in this slice
  nsources_local = 0
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    do isource = 1,NSOURCES
      if (myrank == islice_selected_source(isource)) nsources_local = nsources_local + 1
    enddo
  endif

  ! determines number of times steps for simulation
  call setup_timesteps()

  ! prints source time functions and spectrum to output files
  if (PRINT_SOURCE_TIME_FUNCTION) call print_stf_file()

  ! get information about event name and location
  ! (e.g. needed for SAC seismograms)

  ! The following line is added for get_event_info subroutine.
  ! Because the way NSOURCES_SAC was declared has been changed.
  ! The rest of the changes in this program is just the updates of the subroutines that
  ! I did changes, e.g., adding/removing parameters. by Ebru Bozdag
  call get_event_info_parallel(yr_SAC,jda_SAC,mo_SAC, da_SAC, ho_SAC,mi_SAC,sec_SAC, &
                               event_name_SAC,t_cmt_SAC,t_shift_SAC, &
                               elat_SAC,elon_SAC,depth_SAC,mb_SAC,ms_SAC,cmt_lat_SAC, &
                               cmt_lon_SAC,cmt_depth_SAC,cmt_hdur_SAC,NSOURCES, &
                               Mrr,Mtt,Mpp,Mrt,Mrp,Mtp)

  ! noise simulations ignore the CMTSOLUTIONS sources but employ a noise-spectrum source S_squared instead
  ! checks if anything to do for noise simulations
  if (NOISE_TOMOGRAPHY /= 0) then
    if (myrank == 0) then
      write(IMAIN,*) 'noise simulation will ignore CMT sources'
      call flush_IMAIN()
    endif
  endif

  ! syncs after source setup
  call synchronize_all()

  end subroutine setup_sources

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_stf_constants()

  use specfem_par
  use specfem_par_movie
  implicit none

  ! local parameters
  integer :: isource

  ! makes smaller hdur for movies
  logical,parameter :: USE_SMALLER_HDUR_MOVIE = .false.  ! by default off, to use same HDUR_MOVIE as specified in Par_file

  if (abs(minval(tshift_src)) > TINYVAL) &
    call exit_MPI(myrank,'one tshift_src must be zero, others must be positive')

  ! filter source time function by Gaussian with hdur = HDUR_MOVIE when writing movies or shakemaps
  if (MOVIE_SURFACE .or. MOVIE_VOLUME) then
    ! smaller hdur_movie will do
    if (USE_SMALLER_HDUR_MOVIE) then
      ! hdur_movie gets assigned an automatic value based on the simulation resolution
      ! this will make that a bit smaller to have a higher-frequency movie output
      HDUR_MOVIE = 0.5 * HDUR_MOVIE
    endif

    ! new hdur for simulation
    hdur = sqrt(hdur**2 + HDUR_MOVIE**2)
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Each source is being convolved with HDUR_MOVIE = ',sngl(HDUR_MOVIE)
      write(IMAIN,*) 'Total source hdur = ',sngl(hdur),'(s)'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! convert the half duration for triangle STF to the one for Gaussian STF
  hdur_Gaussian(:) = hdur(:)/SOURCE_DECAY_MIMIC_TRIANGLE

  ! define t0 as the earliest start time
  if (USE_FORCE_POINT_SOURCE) then
    ! point force sources
    ! (might start depending on the frequency given by hdur)
    ! note: point force sources will give the dominant frequency in hdur, thus the main period is 1/hdur.
    !       also, these sources might use a Ricker source time function instead of a Gaussian.
    !       For a Ricker source time function, a start time ~1.2 * main_period is a good choice.
    t0 = 0.d0
    do isource = 1,NSOURCES
      select case(force_stf(isource))
      case (0)
        ! Gaussian source time function
        t0 = min(t0,1.5d0 * (tshift_src(isource) - hdur(isource)))
      case (1)
        ! Ricker source time function
        t0 = min(t0,1.2d0 * (tshift_src(isource) - 1.0d0/hdur(isource)))
      case (2)
        ! Heaviside
        t0 = min(t0,1.5d0 * (tshift_src(isource) - hdur(isource)))
      case (3)
        ! Monochromatic
        t0 = 0.d0
      case (4)
        ! Gaussian source time function by Meschede et al. (2011)
        t0 = min(t0,1.5d0 * (tshift_src(isource) - hdur(isource)))
      case default
        stop 'unsupported force_stf value!'
      end select
    enddo
    ! start time defined as positive value, will be subtracted
    t0 = - t0
  else
    ! moment tensors
    if (USE_MONOCHROMATIC_CMT_SOURCE) then
    ! (based on monochromatic functions)
      t0 = 0.d0
    else
    ! (based on Heaviside functions)
      t0 = - 1.5d0 * minval( tshift_src(:) - hdur(:) )
    endif
  endif

  ! uses an external file for source time function, which starts at time 0.0
  if (EXTERNAL_SOURCE_TIME_FUNCTION) then
    hdur(:) = 0.d0
    t0      = 0.d0
  endif

  ! checks if user set USER_T0 to fix simulation start time
  ! note: USER_T0 has to be positive
  if (USER_T0 > 0.d0) then
    ! user cares about origin time and time shifts of the CMTSOLUTION
    ! and wants to fix simulation start time to a constant start time
    ! time 0 on time axis will correspond to given origin time

    ! notifies user
    if (myrank == 0) then
      write(IMAIN,*) 'USER_T0: ',USER_T0
      write(IMAIN,*) 't0: ',t0,'min_tshift_src_original: ',min_tshift_src_original
      write(IMAIN,*)
    endif

    ! checks if automatically set t0 is too small
    ! note: min_tshift_src_original can be a positive or negative time shift (minimum from all tshift)
    if (t0 <= USER_T0 + min_tshift_src_original) then
      ! by default, tshift_src(:) holds relative time shifts with a minimum time shift set to zero
      ! re-adds (minimum) original time shift such that sources will kick in
      ! according to their absolute time shift
      tshift_src(:) = tshift_src(:) + min_tshift_src_original

      ! sets new simulation start time such that
      ! simulation starts at t = - t0 = - USER_T0
      t0 = USER_T0

      ! notifies user
      if (myrank == 0) then
        write(IMAIN,*) '  set new simulation start time: ', - t0
        write(IMAIN,*)
      endif
    else
      ! start time needs to be at least t0 for numerical stability
      ! notifies user
      if (myrank == 0) then
        write(IMAIN,*) 'Error: USER_T0 is too small'
        write(IMAIN,*) '       must make one of three adjustments:'
        write(IMAIN,*) '       - increase USER_T0 to be at least: ',t0-min_tshift_src_original
        write(IMAIN,*) '       - decrease time shift in CMTSOLUTION file'
        write(IMAIN,*) '       - decrease hdur in CMTSOLUTION file'
        call flush_IMAIN()
      endif
      call exit_mpi(myrank,'Error USER_T0 is set but too small')
    endif
  else if (USER_T0 < 0.d0) then
    if (myrank == 0) then
      write(IMAIN,*) 'Error: USER_T0 is negative, must be set zero or positive!'
    endif
    call exit_mpi(myrank,'Error negative USER_T0 parameter in constants.h')
  endif

  end subroutine setup_stf_constants

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_timesteps()

  use specfem_par
  implicit none

  ! local parameters
  logical :: is_initial_guess

  ! checks if set by initial guess from read_compute_parameters() routine
  if (NTSTEP_BETWEEN_OUTPUT_SEISMOS == NSTEP) then
    is_initial_guess = .true.
  else
    is_initial_guess = .false.
  endif

  ! from initial guess in read_compute_parameters:
  !    compute total number of time steps, rounded to next multiple of 100
  !    NSTEP = 100 * (int(RECORD_LENGTH_IN_MINUTES * 60.d0 / (100.d0*DT)) + 1)
  !
  ! adds initial t0 time to update number of time steps and reach full record length
  if (abs(t0) > 0.d0) then
    ! note: for zero length, nstep has minimal of 5 timesteps for testing
    !       we won't extend this
    !
    ! careful: do not use RECORD_LENGTH_IN_MINUTES here, as it is only read by the main process
    !          when reading the parameter file, but it is not broadcasted to all other processes
    !          NSTEP gets broadcasted, so we work with this values
    if (NSTEP /= 5) then
      ! extend by bulk of 100 steps to account for half-duration rise time
      NSTEP = NSTEP + 100 * (int( abs(t0) / (100.d0*DT)) + 1)
    endif
  endif

  ! if doing benchmark runs to measure scaling of the code for a limited number of time steps only
  if (DO_BENCHMARK_RUN_ONLY) NSTEP = NSTEP_FOR_BENCHMARK

  ! overrides NSTEP in case specified in Par_file
  if (USER_NSTEP > 0) then
    ! overrides NSTEP
    NSTEP = USER_NSTEP
  endif

  ! checks length for symmetry in case of noise simulations
  if (NOISE_TOMOGRAPHY /= 0) then
    if (mod(NSTEP+1,2) /= 0) then
      print *,'Error noise simulation: invalid time steps = ',NSTEP,', NSTEP + 1 must be a multiple of 2 due to branch symmetry'
      call exit_MPI(myrank,'Error noise simulation: number of timesteps must be symmetric, due to +/- branches')
    endif
  endif

  ! time loop increments end
  it_end = NSTEP

  ! subsets used to save seismograms must not be larger than the whole time series,
  ! otherwise we waste memory
  if (NTSTEP_BETWEEN_OUTPUT_SEISMOS > NSTEP .or. is_initial_guess) NTSTEP_BETWEEN_OUTPUT_SEISMOS = NSTEP

  ! subsets used to save adjoint sources must not be larger than the whole time series,
  ! otherwise we waste memory
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! the default value of NTSTEP_BETWEEN_READ_ADJSRC (0) is to read the whole trace at the same time
    if (NTSTEP_BETWEEN_READ_ADJSRC == 0) NTSTEP_BETWEEN_READ_ADJSRC = NSTEP
    ! limits length
    if (NTSTEP_BETWEEN_READ_ADJSRC > NSTEP) NTSTEP_BETWEEN_READ_ADJSRC = NSTEP
    ! total times steps must be dividable by adjoint source chunks/blocks
    if (mod(NSTEP,NTSTEP_BETWEEN_READ_ADJSRC) /= 0) then
      print *,'Error: NSTEP ',NSTEP,' not a multiple of NTSTEP_BETWEEN_READ_ADJSRC ',NTSTEP_BETWEEN_READ_ADJSRC
      print *,'       Please change NTSTEP_BETWEEN_READ_ADJSRC in the Par_file!'
      stop 'Error: mod(NSTEP,NTSTEP_BETWEEN_READ_ADJSRC) must be zero! Please modify Par_file and rerun solver'
    endif
  endif

  ! buffering with undo_attenuation
  NT_DUMP_ATTENUATION = NT_DUMP_ATTENUATION_VAL
  if (UNDO_ATTENUATION) then
    ! makes sure buffer size is not too big for total time length
    !
    ! note: NSTEP must not be a multiple of NT_DUMP_ATTENUATION.
    !       the value from the header file NT_DUMP_ATTENUATION_VAL gives the optimal (maximum) number of time steps for buffering
    if (NSTEP < NT_DUMP_ATTENUATION) NT_DUMP_ATTENUATION = NSTEP
  endif

  ! debug
  !if (myrank == 0 ) print *,'setup time steps = ',NSTEP,' t0 = ',t0,' DT = ',DT

  end subroutine setup_timesteps

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_receivers()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

! check for imbalance of distribution of receivers or of adjoint sources
  logical, parameter :: CHECK_FOR_IMBALANCE = .false.

  ! local parameters
  integer :: irec,isource,nrec_tot_found,i,iproc,ier
  integer :: nrec_simulation
  integer :: nadj_files_found,nadj_files_found_tot
  integer,dimension(0:NPROCTOT_VAL-1) :: tmp_rec_local_all
  integer :: maxrec,maxproc(1)
  double precision :: sizeval

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'receivers:'
    call flush_IMAIN()
  endif

  ! allocate memory for receiver arrays
  allocate(islice_selected_rec(nrec), &
           ispec_selected_rec(nrec), &
           xi_receiver(nrec), &
           eta_receiver(nrec), &
           gamma_receiver(nrec), &
           nu_rec(NDIM,NDIM,nrec), &
           stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating receiver arrays')
  ! initializes arrays
  islice_selected_rec(:) = -1
  ispec_selected_rec(:) = 0
  xi_receiver(:) = 0.d0; eta_receiver(:) = 0.d0; gamma_receiver(:) = 0.d0
  nu_rec(:,:,:) = 0.0d0

  allocate(station_name(nrec), &
           network_name(nrec), &
           stlat(nrec), &
           stlon(nrec), &
           stele(nrec), &
           stbur(nrec),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating receiver arrays')
  stlat(:) = 0.d0; stlon(:) = 0.d0; stele(:) = 0.d0; stbur(:) = 0.d0

  !  receivers
  if (myrank == 0) then
    write(IMAIN,*)
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      write(IMAIN,*) 'Total number of receivers = ', nrec
    else
      write(IMAIN,*) 'Total number of adjoint sources = ', nrec
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! locate receivers in the crust in the mesh
  call locate_receivers(yr_SAC,jda_SAC,ho_SAC,mi_SAC,sec_SAC, &
                        theta_source(1),phi_source(1) )

  ! count number of receivers located in this slice
  nrec_local = 0
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    ! note: for 1-chunk simulations, nrec is now the actual number of receivers found in this chunk
    !       (excludes stations located outside of chunk)
    nrec_simulation = nrec
    do irec = 1,nrec
      if (myrank == islice_selected_rec(irec)) nrec_local = nrec_local + 1
    enddo
  else
    nrec_simulation = NSOURCES
    do isource = 1, NSOURCES
      if (myrank == islice_selected_source(isource)) nrec_local = nrec_local + 1
    enddo
  endif

  ! counter for adjoint receiver stations in local slice, used to allocate adjoint source arrays
  nadj_rec_local = 0

  ! counts receivers for adjoint simulations
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! temporary counter to check if any files are found at all
    nadj_files_found = 0

    do irec = 1,nrec
      ! checks if slice is valid
      if (islice_selected_rec(irec) < 0 .or. islice_selected_rec(irec) > NPROCTOT_VAL-1) &
        call exit_MPI(myrank,'something is wrong with the source slice number in adjoint simulation')

      ! adjoint receiver station in this process slice
      if (myrank == islice_selected_rec(irec)) then
        ! updates counter
        nadj_rec_local = nadj_rec_local + 1

        ! checks **net**.**sta**.**MX**.adj files for correct number of time steps
        if (READ_ADJSRC_ASDF) then
          ! ASDF format
          call check_adjoint_sources_asdf(irec,nadj_files_found)
        else
          ! ASCII format
          call check_adjoint_sources(irec,nadj_files_found)
        endif
      endif
    enddo

    ! checks if any adjoint source files found at all
    call sum_all_i(nadj_files_found,nadj_files_found_tot)
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '    ',nadj_files_found_tot,' adjoint component traces found in all slices'
      call flush_IMAIN()

      ! main process checks if any adjoint files found
      if (nadj_files_found_tot == 0) then
        print *,'Error no adjoint traces found: ',nadj_files_found_tot
        print *,'in directory : ','SEM/'
        print *
        call exit_MPI(myrank,'no adjoint traces found, please check adjoint sources in directory SEM/')
      endif
    endif
  endif

  ! check that the sum of the number of receivers in each slice is nrec (or nsources for adjoint simulations)
  call sum_all_i(nrec_local,nrec_tot_found)
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'found a total of ', nrec_tot_found,' receivers in all slices'
    ! checks for number of receivers
    ! note: for 1-chunk simulations, nrec_simulations is the number of receivers/sources found in this chunk
    if (nrec_tot_found /= nrec_simulation) then
      call exit_MPI(myrank,'problem when dispatching the receivers')
    else
      write(IMAIN,*) 'this total is okay'
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! check that the sum of the number of receivers in each slice is nrec
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    if (myrank == 0 .and. nrec_tot_found /= nrec) &
      call exit_MPI(myrank,'total number of receivers is incorrect')
  endif

  ! synchronizes before info output
  call synchronize_all()

  ! statistics about allocation memory for seismograms & adj_sourcearrays
  ! user output info
  ! sources
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    ! user output
    if (myrank == 0) then
      ! note: all process allocate the full sourcearrays array
      ! sourcearrays(NDIM,NGLLX,NGLLY,NGLLZ,NSOURCES)
      sizeval = dble(NSOURCES) * dble(NDIM * NGLLX * NGLLY * NGLLZ * CUSTOM_REAL / 1024. / 1024. )
      ! outputs info
      write(IMAIN,*) 'source arrays:'
      write(IMAIN,*) '  number of sources is ',NSOURCES
      write(IMAIN,*) '  size of source array                 = ', sngl(sizeval),'MB'
      write(IMAIN,*) '                                       = ', sngl(sizeval/1024.d0),'GB'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! for main process seismogram output
  if (myrank == 0 .and. WRITE_SEISMOGRAMS_BY_MAIN) then
    ! counts number of local receivers for each slice
    allocate(islice_num_rec_local(0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating islice_num_rec_local')

    islice_num_rec_local(:) = 0
    do i = 1,nrec_simulation
      if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
        iproc = islice_selected_rec(i)
      else
        ! adjoint simulations
        iproc = islice_selected_source(i)
      endif

      ! checks iproc value
      if (iproc < 0 .or. iproc >= NPROCTOT_VAL) then
        print *,'Error :',myrank,'iproc = ',iproc,'NPROCTOT = ',NPROCTOT_VAL
        call exit_mpi(myrank,'Error iproc in islice_num_rec_local')
      endif

      ! sums number of receivers for each slice
      islice_num_rec_local(iproc) = islice_num_rec_local(iproc) + 1
    enddo
  endif

  ! seismograms
  ! gather from secondary processes on main
  tmp_rec_local_all(:) = 0
  tmp_rec_local_all(0) = nrec_local
  if (NPROCTOT_VAL > 1) then
    call gather_all_singlei(nrec_local,tmp_rec_local_all,NPROCTOT_VAL)
  endif
  ! user output
  if (myrank == 0) then
    ! determines maximum number of local receivers and corresponding rank
    maxrec = maxval(tmp_rec_local_all(:))
    ! note: MAXLOC will determine the lower bound index as '1'.
    maxproc = maxloc(tmp_rec_local_all(:)) - 1
    ! seismograms array size in MB
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      ! seismograms need seismograms(NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS)
      sizeval = dble(maxrec) * dble(NDIM * NTSTEP_BETWEEN_OUTPUT_SEISMOS * CUSTOM_REAL / 1024. / 1024. )
    else
      ! adjoint seismograms need seismograms(NDIM*NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS)
      sizeval = dble(maxrec) * dble(NDIM * NDIM * NTSTEP_BETWEEN_OUTPUT_SEISMOS * CUSTOM_REAL / 1024. / 1024. )
    endif
    ! outputs info
    write(IMAIN,*) 'seismograms:'
    if (WRITE_SEISMOGRAMS_BY_MAIN) then
      write(IMAIN,*) '  seismograms written by main process only'
    else
      write(IMAIN,*) '  seismograms written by all processes'
    endif
    write(IMAIN,*) '  writing out seismograms at every NTSTEP_BETWEEN_OUTPUT_SEISMOS = ',NTSTEP_BETWEEN_OUTPUT_SEISMOS
    write(IMAIN,*) '  maximum number of local receivers is ',maxrec,' in slice ',maxproc(1)
    write(IMAIN,*) '  size of maximum seismogram array       = ', sngl(sizeval),'MB'
    write(IMAIN,*) '                                         = ', sngl(sizeval/1024.d0),'GB'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! adjoint sources
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! gather from secondary processes on main
    tmp_rec_local_all(:) = 0
    tmp_rec_local_all(0) = nadj_rec_local
    if (NPROCTOT_VAL > 1) then
      call gather_all_singlei(nadj_rec_local,tmp_rec_local_all,NPROCTOT_VAL)
    endif
    ! user output
    if (myrank == 0) then
      ! determines maximum number of local receivers and corresponding rank
      maxrec = maxval(tmp_rec_local_all(:))
      ! note: MAXLOC will determine the lower bound index as '1'.
      maxproc = maxloc(tmp_rec_local_all(:)) - 1
      !do i = 1, NPROCTOT_VAL
      !  if (tmp_rec_local_all(i) > maxrec) then
      !    maxrec = tmp_rec_local_all(i)
      !    maxproc = i-1
      !  endif
      !enddo
      ! source_adjoint size in MB
      ! source_adjoint(NDIM,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC)
      sizeval = dble(maxrec) * dble(NDIM * NTSTEP_BETWEEN_READ_ADJSRC * CUSTOM_REAL / 1024. / 1024. )
      ! note: in case IO_ASYNC_COPY is set, and depending of NSTEP_SUB_ADJ,
      !       this memory requirement might double.
      !       at this point, NSTEP_SUB_ADJ is not set yet...
      if (IO_ASYNC_COPY .and. ceiling( dble(NSTEP)/dble(NTSTEP_BETWEEN_READ_ADJSRC) ) > 1) then
        !buffer_source_adjoint(NDIM,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC)
        sizeval = sizeval + dble(maxrec) * dble(NDIM * NTSTEP_BETWEEN_READ_ADJSRC * CUSTOM_REAL / 1024. / 1024. )
      endif
      ! outputs info
      write(IMAIN,*) 'adjoint source arrays:'
      write(IMAIN,*) '  reading adjoint sources at every NTSTEP_BETWEEN_READ_ADJSRC = ',NTSTEP_BETWEEN_READ_ADJSRC
      if (IO_ASYNC_COPY) then
        write(IMAIN,*) '  using asynchronous buffer for file I/O of adjoint sources'
      endif
      write(IMAIN,*) '  maximum number of local adjoint sources is ',maxrec,' in slice ',maxproc(1)
      write(IMAIN,*) '  size of maximum adjoint source array = ', sngl(sizeval),'MB'
      write(IMAIN,*) '                                       = ', sngl(sizeval/1024.d0),'GB'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

! check for imbalance of distribution of receivers or of adjoint sources
  if (CHECK_FOR_IMBALANCE .and. NPROCTOT_VAL > 1) then
    call gather_all_singlei(nrec_local,tmp_rec_local_all,NPROCTOT_VAL)
    if (myrank == 0) then
      open(unit=9977,file='imbalance_of_nrec_local.dat',status='unknown')
      do i = 0,NPROCTOT_VAL-1
        write(9977,*) i,tmp_rec_local_all(i)
      enddo
      close(9977)
    endif

    call gather_all_singlei(nadj_rec_local,tmp_rec_local_all,NPROCTOT_VAL)
    if (myrank == 0) then
      open(unit=9977,file='imbalance_of_nadj_rec_local.dat',status='unknown')
      do i = 0,NPROCTOT_VAL-1
        write(9977,*) i,tmp_rec_local_all(i)
      enddo
      close(9977)
    endif

    if (myrank == 0) then
      open(unit=9977,file='plot_imbalance_histogram.gnu',status='unknown')
      write(9977,*) '#set terminal x11'
      write(9977,*) 'set terminal wxt'
      write(9977,*) '#set terminal gif'
      write(9977,*) '#set output "imbalance_histogram.gif"'
      write(9977,*)
      write(9977,*) 'set xrange [1:',NPROCTOT_VAL,']'
      write(9977,*) '#set xtics 0,0.1,1'
      write(9977,*) 'set boxwidth 1.'
      write(9977,*) 'set xlabel "Mesh slice number"'
      write(9977,*)
      write(9977,*) 'set ylabel "Number of receivers in that mesh slice"'
      write(9977,*) 'plot "imbalance_of_nrec_local.dat" with boxes'
      write(9977,*) 'pause -1 "hit any key..."'
      write(9977,*)
      write(9977,*) 'set ylabel "Number of adjoint sources in that mesh slice"'
      write(9977,*) 'plot "imbalance_of_nadj_rec_local.dat" with boxes'
      write(9977,*) 'pause -1 "hit any key..."'
      close(9977)
    endif

    call synchronize_all()
  endif

  end subroutine setup_receivers


!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources_receivers_VTKfile()

  use specfem_par, only: myrank,OUTPUT_FILES,NSOURCES,nrec,MAX_STRING_LEN
  implicit none

  ! local parameters
  character(len=MAX_STRING_LEN) :: filename,filename_new
  character(len=MAX_STRING_LEN) :: command
  ! user output
  if (myrank == 0) then

    ! finishes VTK file
    !  we should know NSOURCES+nrec at this point...
    ! creates source/receiver location file
    filename = trim(OUTPUT_FILES)//'/sr_tmp.vtk'
    filename_new = trim(OUTPUT_FILES)//'/sr.vtk'
    write(command, &
  "('sed -e ',a1,'s/POINTS.*/POINTS',i6,' float/',a1,'<',a,'>',a)")&
      "'",NSOURCES + nrec,"'",trim(filename),trim(filename_new)

    ! note: this system() routine is non-standard Fortran
    call system_command(command)

    ! only extract receiver locations and remove temporary file
    filename_new = trim(OUTPUT_FILES)//'/receiver.vtk'
    write(command, &
  "('awk ',a1,'{if (NR < 5) print $0;if (NR == 6)&
   &print ',a1,'POINTS',i6,' float',a1,';if (NR > 5+',i6,')print $0}',a1,'<',a,'>',a)")&
      "'",'"',nrec,'"',NSOURCES,"'",trim(filename),trim(filename_new)

    ! note: this system() routine is non-standard Fortran
    call system_command(command)

    ! only extract source locations and remove temporary file
    filename_new = trim(OUTPUT_FILES)//'/source.vtk'
    write(command, &
  "('awk ',a1,'{if (NR < 6 + ',i6,') print $0}END{print}',a1,'<',a,'>',a,'; rm -f ',a)")&
      "'",NSOURCES,"'",trim(filename),trim(filename_new),trim(filename)

    ! note: this system() routine is non-standard Fortran
    call system_command(command)

  endif

  end subroutine setup_sources_receivers_VTKfile


!-------------------------------------------------------------------------------------------------


  subroutine unique_inverse(A, Nin, unique, unique_idx, inv, Nout)

    ! This subroutine computes the uniqe array of integers contained in an
    ! integer array. The function is not performance optimized. So it may be
    ! slow for large arrays. I use it below in setup_green_locations, to compute
    ! a new, subset ibool array, to get only the tagged elements.
    ! Note that you have to grab the unqiue values likeso:
    ! >> uniq     = unique(1:Nout)
    ! >> uniq_idx = unique_idx(1:Nout)

    implicit none

    ! In
    integer, intent(in) :: Nin
    integer, dimension(Nin), intent(in) :: A

    ! out
    integer, intent(out) :: Nout
    ! integer, dimension(:), intent(out), allocatable :: unique, unique_idx

    ! Local
    integer :: i
    integer, dimension(:), allocatable :: idx
    integer, dimension(Nin) :: unique
    integer, dimension(Nin) :: unique_idx
    integer, dimension(Nin), intent(out) :: inv

    Nout = 0

    unique(:) = 0

    do i=1,Nin

      if (any( A(i) == unique) .eqv. .false.) then
        Nout = Nout + 1
        unique(Nout) = A(i)
        unique_idx(Nout) = i

        inv(i) = Nout


      ! If the value is already in bigunique, check where it is using minloc
      else
        idx = minloc(abs(unique(1:Nout)-A(i)))
        inv(i) = idx(1)
      endif
    enddo

    ! allocate(unique(Nout), unique_idx(Nout))

    ! unique(:) = big_unique(1:Nout)
    ! unique_idx(:) = big_unique_idx(1:Nout)

  end subroutine
!-------------------------------------------------------------------------------------------------

  subroutine setup_green_locations()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  character(len=MAX_STRING_LEN) :: filename
  ! check for imbalance of distribution of receivers or of adjoint sources
  logical, parameter :: CHECK_FOR_IMBALANCE = .false.

  ! local parameters
  integer :: igf, ngf_tot_found,ier,ix
  integer :: i,j,k,l, ispec
  integer, dimension(:), allocatable :: idx
  integer :: inum_neighbor,MAX_NEIGHBORS
  integer :: igllx, iglly, igllz, ispec_sub, ibel, counter,ispec_neighbor, num_neighbors
  integer :: iglob, iglob_counter, igf_counter
  integer :: ngf_simulation
  integer,dimension(0:NPROCTOT_VAL-1) :: tmp_gf_local_all
  double precision :: sizeval
  integer, dimension(NGLLX, NGLLY, NGLLZ) :: checkarray
  integer, dimension(:,:), allocatable :: islicespec_selected_gf_loc
  integer, dimension(:), allocatable :: iglob_tmp
  logical, dimension(:), allocatable :: mask
  integer, dimension(:), allocatable :: rmask
  integer, dimension(:), allocatable :: index_vector
  integer, dimension(:), allocatable :: ispec_unique_gf_loc_local
  integer, dimension(:), allocatable :: islice_unique_gf_loc_local
  integer, dimension(:), allocatable :: ibool_flat, inv, iglob_cm2gf_tmp
  integer, dimension(NSPEC_CRUST_MANTLE) :: ispec_mask, islice_mask
  integer, dimension(NPROCTOT_VAL) :: ngf_unique_array, offset
  integer, dimension(:), allocatable :: tmp_adjacency

  ! timer MPI
  double precision :: time1,tCPU
  double precision, external :: wtime

  MAX_NEIGHBORS=50

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Green function locations:'
    call flush_IMAIN()
  endif
  ! write(*,*) 'myrank', myrank, 'nspec', NSPEC_CRUST_MANTLE


  ! allocate memory for green function arrays
  allocate(xi_gf_loc(ngf), &
           eta_gf_loc(ngf), &
           gamma_gf_loc(ngf),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating green function arrays 0.0')

  allocate(nu_gf_loc(NDIM,NDIM,ngf),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating green function arrays 0.1')

  allocate(islice_selected_gf_loc(ngf), &
           ispec_selected_gf_loc(ngf), &
           islicespec_selected_gf_loc(2,ngf), &
           stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating green function arrays 1')

  ! initializes arrays
  islice_selected_gf_loc(:) = -1
  ispec_selected_gf_loc(:) = 0
  islicespec_selected_gf_loc(:,:) = 0
  xi_gf_loc(:) = 0.d0; eta_gf_loc(:) = 0.d0; gamma_gf_loc(:) = 0.d0
  nu_gf_loc(:,:,:) = 0.0d0

  allocate(gf_loc_lat(ngf), &
           gf_loc_lon(ngf), &
           gf_loc_depth(ngf),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating green function arrays 1 arrays')
  gf_loc_lat(:) = 0.d0; gf_loc_lon(:) = 0.d0; gf_loc_depth(:) = 0.d0

  !  Green Function locations
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Total number of green function locations = ', ngf
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  if (myrank == 0) then
  ! write source and receiver VTK files for Paraview
    filename = trim(OUTPUT_FILES)//'/gf_tmp.vtk'
    open(IOUT_VTK,file=trim(filename),status='unknown',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening temporary file gf_tmp.vtk in setup_sources_reveivers.')
    write(IOUT_VTK,'(a)') '# vtk DataFile Version 2.0'
    write(IOUT_VTK,'(a)') 'Green Function Locations VTK file'
    write(IOUT_VTK,'(a)') 'ASCII'
    write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
    !  LQY -- won't be able to know NSOURCES+nrec at this point...
    write(IOUT_VTK, '(a,i6,a)') 'POINTS ', NSOURCES, ' float'
    ! closing file, rest of information will be appended later on
    close(IOUT_VTK)
  endif

  ! locate receivers in the crust in the mesh
  call locate_green_locations()

! Create column array of slice and spectral element combinations
  islicespec_selected_gf_loc(1,:) = islice_selected_gf_loc
  islicespec_selected_gf_loc(2,:) = ispec_selected_gf_loc

  ! count number of receivers located in this slice
  ngf_local = 0
  ngf_simulation = ngf

  do igf = 1,ngf
    if (myrank == islice_selected_gf_loc(igf)) ngf_local = ngf_local + 1
  enddo

  ! check that the sum of the number of receivers in each slice is nrec (or nsources for adjoint simulations)
  call sum_all_i(ngf_local,ngf_tot_found)

  ! Get unique number of elements and allocate new array of slices and sources
  if (myrank == 0) then

    allocate( &
      mask(ngf), &
      rmask(ngf), stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating mask and rmask')

    mask(:) = .true.

    do ix = ngf,2,-1
      mask(ix) = .not.(any(&
        islicespec_selected_gf_loc(1,:ix-1)==islicespec_selected_gf_loc(1,ix).and.&
        islicespec_selected_gf_loc(2,:ix-1)==islicespec_selected_gf_loc(2,ix)))
    enddo

    ! Get total number of unique elements (IMPORTANT FOR ADIOS FILE ALLOCATION)
    rmask(:) = 0

    where(mask) rmask = 1
    ngf_unique = sum(rmask)

  endif

  ! Broadcast the total number of unique elements
  call bcast_all_singlei(ngf_unique)

  ! Make an index vector
  allocate(index_vector(ngf_unique), stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating mask index array')
  index_vector(:) = 0

  ! Now copy the unique elements of the slice and spec arrays into the unique arrays
  allocate(islice_unique_gf_loc(ngf_unique), stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating unique slice array')

  allocate(ispec_unique_gf_loc(ngf_unique), stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating unique element array')

  if (myrank==0) then
    index_vector = pack([(ix, ix=1,ngf)],mask)
    deallocate(mask, rmask)
  endif

  call synchronize_all()
  call bcast_all_i(index_vector, ngf_unique)

  do igf=1,ngf_unique
    islice_unique_gf_loc(igf) = islice_selected_gf_loc(index_vector(igf))
    ispec_unique_gf_loc(igf) = ispec_selected_gf_loc(index_vector(igf))
  enddo

  if (myrank==0) then
    ! outputs info
    write(IMAIN,*)
    write(IMAIN,*) 'Found a total of ', ngf_unique, 'tagged elements.'
  endif

  call synchronize_all()

  ! This section is adding a neighbouring elements by brute-force checking the
  ! number of overlapping points
  ispec_mask(:) = 0
  islice_mask(:) = 9999
  ibel=0

  ! If Use buffer elements create buffer around elements
  if (USE_BUFFER_ELEMENTS) then

    if (myrank==0) then
      ! outputs info
      write(IMAIN,*)
      write(IMAIN,*) 'Computing element buffer for ...'
      ! get MPI starting time
      time1 = wtime()
    endif


    ! Expand 1 buffer element at a time (this is a very exponential loop)
    do ibel=1,NUMBER_OF_BUFFER_ELEMENTS

      if (myrank==0) then
        write(IMAIN,*) '    ---> ', ibel, 'element.'
      endif

      ! Loop over unique elements
      do igf=1, ngf_unique

        ! Check whether element is in current slice
        if (islice_unique_gf_loc(igf)==myrank) then

          ! Get element
          ispec = ispec_unique_gf_loc(igf)

          ! Get number of neighbors from adjacency offset vector
          num_neighbors = xadj(ispec+1)-xadj(ispec)

          ! Loop over number of neighbors and add to mask
          do i = 1,num_neighbors
            ! Get neighbor element index
            ispec_neighbor = adjncy(xadj(ispec) + i)

            ! Set index of mask to true.
            ispec_mask(ispec_neighbor) = 1
            islice_mask(ispec_neighbor) = myrank

          enddo

        endif

      enddo

      call synchronize_all()

      ! We can deallocate all of the variables since ispec_mask (a local array)
      ! should contain all necessary info
      deallocate(ispec_unique_gf_loc, islice_unique_gf_loc, stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error deallocating ispec_unique_gf_loc, islice_unique_gf_loc.')

      ! Get total number of local elements
      ngf_unique_local = sum(ispec_mask)

      ! Allocate local ispec, and slice
      allocate(&
          ispec_unique_gf_loc_local(ngf_unique_local), &
          islice_unique_gf_loc_local(ngf_unique_local), &
          stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating array when getting buffer elements')

      ispec_unique_gf_loc_local(:) = 0
      islice_unique_gf_loc_local(:) = 0

      ! Get elements
      counter = 1

      ! Loop over elements
      do ispec=1,NSPEC_CRUST_MANTLE
        if (ispec_mask(ispec)==1) then
          ! if in mask is true add element and slice.
          ispec_unique_gf_loc_local(counter) = ispec
          islice_unique_gf_loc_local(counter) = myrank
          counter = counter + 1
        endif
      enddo

      ! Get complete number of unique elements
      ngf_unique = 0
      call sum_all_i(ngf_unique_local, ngf_unique)
      call bcast_all_singlei(ngf_unique)

      ! Get array of unique elements on each proc at root
      ngf_unique_array(:) = 0
      call gather_all_singlei(ngf_unique_local, ngf_unique_array, NPROCTOT_VAL)
      call bcast_all_i(ngf_unique_array, NPROCTOT_VAL)

      ! Get ispecs and slices from all processes to root
      allocate(ispec_unique_gf_loc(ngf_unique), islice_unique_gf_loc(ngf_unique), stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating array when getting buffer elements 2')

      call synchronize_all()

      ! Set Gatherv displacements
      offset(:) = 0
      do i=2,NPROCTOT_VAL
        offset(i) = offset(i-1) + ngf_unique_array(i-1)
      enddo

      ! Bringing back all info to the root process
      call gatherv_all_i(&
          ispec_unique_gf_loc_local, ngf_unique_local, &
          ispec_unique_gf_loc, ngf_unique_array, offset, ngf_unique, NPROCTOT_VAL)

      call gatherv_all_i(&
          islice_unique_gf_loc_local, ngf_unique_local, &
          islice_unique_gf_loc, ngf_unique_array, offset, ngf_unique, NPROCTOT_VAL)

      ! Broadcasting it back to all processes
      call bcast_all_i(ispec_unique_gf_loc, ngf_unique)
      call bcast_all_i(islice_unique_gf_loc, ngf_unique)

      deallocate(ispec_unique_gf_loc_local, islice_unique_gf_loc_local)

    enddo

    if (myrank==0) then
        tCPU = wtime() - time1
        ! outputs info
        write(IMAIN,*)
        write(IMAIN,*) 'After computing buffer elements, we have ', ngf_unique, ' elements across all slices.'
        write(IMAIN,*) 'Buffering took', tCPU, ' seconds.'
    endif

  endif

  ! main process broadcasts the results to all the slices
  call bcast_all_singlei(ngf_unique)
  call bcast_all_i(islice_unique_gf_loc, ngf_unique)
  call bcast_all_i(ispec_unique_gf_loc, ngf_unique)

  ! synchronizes to get right timing
  call synchronize_all()

  ! Get number of unique entries in a slice
  ngf_unique_local = 0
  do igf=1,ngf_unique
    if (islice_unique_gf_loc(igf)==myrank) then
      ngf_unique_local = ngf_unique_local + 1
    endif
  enddo

  call synchronize_all()

  ! ------------
  ! This section is quite important since it sets up the global to local
  ! addressing for the subset of elements.

  allocate(islice_out_gf_loc(ngf), &
           ispec_out_gf_loc(ngf), &
           stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating islice/ispec_out')

  if (ngf_unique_local .gt. 0) then

    ! Get ibool array for output Green function database
    allocate(ibool_GF(NGLLX, NGLLY, NGLLZ, ngf_unique_local), &
             ibool_flat(NGLLX*NGLLY*NGLLZ*ngf_unique_local), &
             inv(NGLLX*NGLLY*NGLLZ*ngf_unique_local), &
             ispec_cm2gf(ngf_unique_local), &
             islice_cm2gf(ngf_unique_local), &
             stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating ibool_GF, or iglob_tmp array')

    ! Conversion arrays from full crust_mantle element array to small
    ! Green function array
    ibool_GF(:,:,:,:) = 0
    ispec_cm2gf(:) = 0
    islice_cm2gf(:) = 0
    islice_out_gf_loc(:) = 0
    ispec_out_gf_loc(:) = 0

    ! Initialize counters to count local elements and local coordinates
    iglob_counter = 0
    igf_counter = 0

    do igf=1,ngf_unique
      if (islice_unique_gf_loc(igf)==myrank) then

        ! Count new elements
        igf_counter = igf_counter + 1

        ! For each element in the new local array get the element of
        ! the original crust_mantle array
        ispec_cm2gf(igf_counter) = ispec_unique_gf_loc(igf)
      endif
    enddo

    ! Allocate temporary placeholder for unique, and unique,index array
    allocate(iglob_cm2gf_tmp(size(ibool_GF)), iglob_tmp(size(ibool_GF)), stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating unique temp arrays')

    ! Flatten subset ibool array for unique-ing
    ibool_flat = pack(ibool_crust_mantle(:,:,:,ispec_cm2gf), .true.)

    if (myrank==0) then
        ! outputs info
        write(IMAIN,*)
        write(IMAIN,*) 'Computing the unique values in the subset ibool array and the inverse'
        write(IMAIN,*) 'of the unique indexing for readdressing ...'
        ! get MPI starting time
        time1 = wtime()
    endif



    ! Get unique values of ibool, indeces of those values, and inverse.
    !                    array      arraysize       unique values, idx,    inverse, Nuniq
    call unique_inverse(ibool_flat, size(ibool_GF), iglob_cm2gf_tmp, iglob_tmp, inv, iglob_counter)

    ! Large array better deallocate quickly.
    deallocate(ibool_flat, iglob_tmp, stat=ier) !!!!! REUSING iglob_tmp in next line because lazy
    if (ier /= 0 ) call exit_MPI(myrank,'Error deallocating bool_flat, iglob_tmp')

    ! Allocate array and placeholder array.
    allocate(iglob_cm2gf(iglob_counter), iglob_tmp(iglob_counter), stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating iglob_cm2gf')

    ! Get subset of unique values
    iglob_cm2gf(:) = iglob_cm2gf_tmp(1:iglob_counter)

    ! Deallocate large array
    deallocate(iglob_cm2gf_tmp, stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error deallocating iglob_cm2gf_tmp')

    ! Create temporary sequence for ibool_GF creation
    do i=1,iglob_counter
      iglob_tmp(i) = i
    enddo

    ! Created ibool GF from inverse unique array
    ibool_GF(:,:,:,:) = reshape(iglob_tmp(inv), shape(ibool_GF))

    deallocate(inv, iglob_tmp, stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error deallocating iglob_tmp')

    ! Total number of Green function coordinates in terms of elements
    NGLOB_GF = iglob_counter

    if (myrank==0) then
        tCPU = wtime() - time1
        ! outputs info
        write(IMAIN,*)
        write(IMAIN,*) '    ... readdressing of the subset elements took ', tCPU, 'seconds.'

    endif

  endif

  ! Get ibool array for  output Green function database
  call synchronize_all()


  ! Convert crust mantle ispec to green function database ispec
  do igf=1,ngf
    if (islice_selected_gf_loc(igf)==myrank) then

        islice_out_gf_loc(igf) = myrank

        do i=1,ngf_unique_local
          if (ispec_cm2gf(i)==ispec_selected_gf_loc(igf)) then
            ispec_out_gf_loc(igf) = i
          endif
        enddo
    endif
  enddo

  ! Collect the results on rank one
  call synchronize_all()

  do igf=1,ngf
    if (islice_selected_gf_loc(igf) /= 0) then
      if (islice_selected_gf_loc(igf) == myrank) then
        call send_singlei(ispec_out_gf_loc(igf), 0, igf)
        call send_singlei(islice_out_gf_loc(igf), 0, igf+ngf)
      endif
      if (myrank==0) then
        call recv_singlei(ispec_out_gf_loc(igf), islice_selected_gf_loc(igf), igf)
        call recv_singlei(islice_out_gf_loc(igf), islice_selected_gf_loc(igf), igf+ngf)
      endif
    endif
  enddo

  ! Get adjacency vectors for to enable source inversion arrays.
  allocate(xadj_gf(ngf_unique_local + 1),stat=ier)
  if (ier /= 0) stop 'Error allocating xadj'
  allocate(tmp_adjacency(MAX_NEIGHBORS*ngf_unique_local),stat=ier)
  if (ier /= 0) stop 'Error allocating temp adj'

  ! Get neighbours.
  xadj_gf(:) = 0
  inum_neighbor = 0
  tmp_adjacency(:) = 0

  if (myrank==0) then
        ! outputs info
        write(IMAIN,*)
        write(IMAIN,*) 'Gathering neighbors of the subset elements ... '
        time1 = wtime()
  endif

  do igf=1,ngf_unique_local

    ! Get element
    ispec = ispec_cm2gf(igf)

    ! Get neighbors
    num_neighbors = xadj(ispec+1)-xadj(ispec)

    ! Loop over neighbors in full mesh
    do i = 1,num_neighbors

      ! get neighbor from global adjacency
      ispec_neighbor = adjncy(xadj(ispec) + i)

      ! Check whether global neighbor is also a local neighbor.
      if (any(ispec_neighbor==ispec_cm2gf)) then

        ! If is neighbor increase total neighbor counter.
        inum_neighbor = inum_neighbor + 1

        ! Get indeces of neighbors
        idx = pack([(ix,ix=1,size(ispec_cm2gf))],ispec_neighbor==ispec_cm2gf)

        ! Add neighbor to adjacency vector
        tmp_adjacency(inum_neighbor) = idx(1)

      endif

    enddo

    ! Add total event counter to adjacency vetor
    xadj_gf(igf+1) = inum_neighbor

  enddo

  if (myrank==0) then
        tCPU = wtime() - time1
        ! outputs info
        write(IMAIN,*) '    ... took ', tCPU, ' seconds.'
  endif

  ! Allocate final adjacency array
  allocate(adjncy_gf(inum_neighbor),stat=ier)
  if (ier /= 0) stop 'Error allocating temp adj'

  ! Define final adjacency array
  adjncy_gf(1:inum_neighbor) = tmp_adjacency(1:inum_neighbor)

  ! Define neighbor counter
  num_neighbors_all_gf = inum_neighbor

  ! Deallocate tmp adjacency
  deallocate(tmp_adjacency)

  call synchronize_all()

  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'found a total of ',ngf_tot_found,' green function locations in all slices'
    write(IMAIN,*) 'found a total of ',ngf_unique,' unique green function elements in all slices'

    ! checks for number of receivers
    ! note: for 1-chunk simulations, nrec_simulations is the number of receivers/sources found in this chunk
    if (ngf_tot_found /= ngf_simulation) then
      call exit_MPI(myrank,'problem when dispatching the receivers')
    else
      write(IMAIN,*) 'this total is okay'
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! check that the sum of the number of receivers in each slice is ngf
  if (myrank == 0 .and. ngf_tot_found /= ngf) then
    call exit_MPI(myrank,'total number of green function slices is incorrect')
  endif
  ! synchronizes before info output
  call synchronize_all()

  ! seismograms
  ! gather from secondary processes on main
  tmp_gf_local_all(:) = 0
  tmp_gf_local_all(0) = ngf_local
  if (NPROCTOT_VAL > 1) then
    call gather_all_singlei(ngf_local,tmp_gf_local_all,NPROCTOT_VAL)
  endif

  ! user output
  if (myrank == 0) then
    ! seismograms need seismograms(NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS)
    sizeval = dble(ngf) * dble(NGLLX * NGLLY * NGLLZ * CUSTOM_REAL / 1024. / 1024. )
    ! outputs info
    write(IMAIN,*) 'Green Function:'
    ! NOTE: Needs to be edited in the future ACCOUNTING for total of SAVED steps
    write(IMAIN,*) '  writing out Green functions at every NTSTEP_BETWEEN_OUTPUT_SEISMOS = ',NTSTEP_BETWEEN_OUTPUT_SEISMOS
    write(IMAIN,*) '  Number of unique elements                 = ', ngf_unique
    write(IMAIN,*) '  Storage per element per step (no overlap) = ', sngl(sizeval),'MB'
    write(IMAIN,*) '                                            = ', sngl(sizeval/1024.d0),'GB'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! synchronizes to get right timing
  call synchronize_all()

  end subroutine setup_green_locations


!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources_precompute_arrays()

  use specfem_par
  use specfem_par_crustmantle
  implicit none

  ! local parameters
  integer :: ier
  integer(kind=8) :: arraysize

  ! allocates source arrays
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    ! source interpolated on all GLL points in source element
    allocate(sourcearrays(NDIM,NGLLX,NGLLY,NGLLZ,NSOURCES),stat=ier)
    if (ier /= 0 ) then
      print *,'Error rank ',myrank,': allocating sourcearrays failed! number of sources = ',NSOURCES
      call exit_MPI(myrank,'Error allocating sourcearrays')
    endif
    ! initializes
    sourcearrays(:,:,:,:,:) = 0._CUSTOM_REAL

    ! stores source arrays
    call setup_sources_receivers_srcarr()
  else
    ! dummy array
    allocate(sourcearrays(1,1,1,1,1))
  endif

  ! adjoint source arrays
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! initializes adjoint source buffer
    ! reverse indexing
    allocate(iadj_vec(NSTEP),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating iadj_vec')
    ! initializes iadj_vec
    do it = 1,NSTEP
       iadj_vec(it) = NSTEP-it+1  ! default is for reversing entire record, e.g. 3000,2999,..,1
    enddo

    ! total number of adjoint source blocks to read in
    NSTEP_SUB_ADJ = ceiling( dble(NSTEP)/dble(NTSTEP_BETWEEN_READ_ADJSRC) )

    if (nadj_rec_local > 0) then
      allocate(source_adjoint(NDIM,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC),stat=ier)
      if (ier /= 0 ) then
        print *,'Error rank ',myrank,': allocating source_adjoint failed! Please check your memory usage...'
        print *,'  failed number of local adjoint sources = ',nadj_rec_local,' steps = ',NTSTEP_BETWEEN_READ_ADJSRC
        call exit_MPI(myrank,'Error allocating adjoint source_adjoint')
      endif
      source_adjoint(:,:,:) = 0.0_CUSTOM_REAL

      ! additional buffer for asynchronous file i/o
      if (IO_ASYNC_COPY .and. NSTEP_SUB_ADJ > 1) then
        ! allocates read buffer
        allocate(buffer_source_adjoint(NDIM,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC),stat=ier)
        if (ier /= 0 ) call exit_MPI(myrank,'Error allocating array buffer_source_adjoint')
        buffer_source_adjoint(:,:,:) = 0.0_CUSTOM_REAL

        ! array size in bytes (note: the multiplication is split into two line to avoid integer-overflow)
        arraysize = NDIM *  CUSTOM_REAL
        arraysize = arraysize * nadj_rec_local * NTSTEP_BETWEEN_READ_ADJSRC

        ! debug
        !print *,'buffer_sourcearrays: size = ',arraysize,' Bytes = ',arraysize/1024./1024.,'MB'

        ! initializes io thread
        call prepare_adj_io_thread(buffer_source_adjoint,arraysize,nadj_rec_local)
      endif

      ! allocate indexing arrays
      allocate(iadjsrc(NSTEP_SUB_ADJ,2), &
               iadjsrc_len(NSTEP_SUB_ADJ),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating adjoint indexing arrays')
      iadjsrc(:,:) = 0; iadjsrc_len(:) = 0

      ! initializes iadjsrc, iadjsrc_len and iadj_vec
      call setup_sources_receivers_adjindx(NSTEP,NSTEP_SUB_ADJ, &
                                           NTSTEP_BETWEEN_READ_ADJSRC, &
                                           iadjsrc,iadjsrc_len,iadj_vec)
    endif
  endif

  end subroutine setup_sources_precompute_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources_receivers_srcarr()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  integer :: isource,i,j,k,ispec !,iglob

  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearray
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd

  double precision :: xi,eta,gamma
  double precision :: hlagrange
  double precision :: norm

  do isource = 1,NSOURCES

    ! initializes
    sourcearray(:,:,:,:) = 0._CUSTOM_REAL

    !   check that the source slice number is okay
    if (islice_selected_source(isource) < 0 .or. islice_selected_source(isource) > NPROCTOT_VAL-1) &
      call exit_MPI(myrank,'Error: source slice number invalid')

    !   compute source arrays in source slice
    if (myrank == islice_selected_source(isource)) then

      ! element id which holds source
      ispec = ispec_selected_source(isource)

      ! checks bounds
      if (ispec < 1 .or. ispec > NSPEC_CRUST_MANTLE ) &
        call exit_MPI(myrank,'Error: source ispec number invalid')

      ! gets source location
      xi = xi_source(isource)
      eta = eta_source(isource)
      gamma = gamma_source(isource)

!      ! pre-computes source contribution on GLL points
!      call compute_arrays_source(sourcearray,xi,eta,gamma, &
!                          Mxx(isource),Myy(isource),Mzz(isource),Mxy(isource),Mxz(isource),Myz(isource), &
!                          xix_crust_mantle(:,:,:,ispec),xiy_crust_mantle(:,:,:,ispec),xiz_crust_mantle(:,:,:,ispec), &
!                          etax_crust_mantle(:,:,:,ispec),etay_crust_mantle(:,:,:,ispec),etaz_crust_mantle(:,:,:,ispec), &
!                          gammax_crust_mantle(:,:,:,ispec),gammay_crust_mantle(:,:,:,ispec),gammaz_crust_mantle(:,:,:,ispec), &
!                          xigll,yigll,zigll)
!
!      ! point forces, initializes sourcearray, used for simplified CUDA routines
!    !-------------POINT FORCE-----------------------------------------------
!      if (USE_FORCE_POINT_SOURCE) then
!        ! note: for use_force_point_source xi/eta/gamma are in the range [1,NGLL*]
!        iglob = ibool_crust_mantle(nint(xi),nint(eta),nint(gamma),ispec)
!
!        ! sets sourcearrays
!        do k = 1,NGLLZ
!          do j = 1,NGLLY
!            do i = 1,NGLLX
!              if (ibool_crust_mantle(i,j,k,ispec) == iglob) then
!                ! elastic source components
!                sourcearray(:,i,j,k) = nu_source(COMPONENT_FORCE_SOURCE,:,isource)
!              endif
!            enddo
!          enddo
!        enddo
!      endif
!    !-------------POINT FORCE-----------------------------------------------
!
!      ! stores source excitations
!      sourcearrays(:,:,:,:,isource) = sourcearray(:,:,:,:)
!    endif

      ! compute Lagrange polynomials at the source location
      call lagrange_any(xi,NGLLX,xigll,hxis,hpxis)
      call lagrange_any(eta,NGLLY,yigll,hetas,hpetas)
      call lagrange_any(gamma,NGLLZ,zigll,hgammas,hpgammas)

      if (USE_FORCE_POINT_SOURCE) then ! use of FORCESOLUTION files

        ! note: for use_force_point_source xi/eta/gamma are also in the range [-1,1], for exact positioning

        ! initializes source array
        sourcearrayd(:,:,:,:) = 0.0d0

        ! calculates source array for interpolated location
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX
              hlagrange = hxis(i) * hetas(j) * hgammas(k)

              ! elastic source
              norm = sqrt( comp_dir_vect_source_E(isource)**2 &
                         + comp_dir_vect_source_N(isource)**2 &
                         + comp_dir_vect_source_Z_UP(isource)**2 )

              ! checks norm of component vector
              if (norm < TINYVAL) then
                call exit_MPI(myrank,'error force point source: component vector has (almost) zero norm')
              endif

              ! normalizes vector
              comp_dir_vect_source_E(isource) = comp_dir_vect_source_E(isource) / norm
              comp_dir_vect_source_N(isource) = comp_dir_vect_source_N(isource) / norm
              comp_dir_vect_source_Z_UP(isource) = comp_dir_vect_source_Z_UP(isource) / norm

              ! we use a tilted force defined by its magnitude and the projections
              ! of an arbitrary (non-unitary) direction vector on the E/N/Z_UP basis
              !
              ! note: nu_source(iorientation,:,isource) is the rotation matrix from ECEF to local N-E-UP
              !       (defined in src/specfem3D/locate_sources.f90)
              sourcearrayd(:,i,j,k) = factor_force_source(isource) * hlagrange * &
                                      ( nu_source(1,:,isource) * comp_dir_vect_source_N(isource) + &
                                        nu_source(2,:,isource) * comp_dir_vect_source_E(isource) + &
                                        nu_source(3,:,isource) * comp_dir_vect_source_Z_UP(isource) )
            enddo
          enddo
        enddo

        ! distinguish between single and double precision for reals
        sourcearray(:,:,:,:) = real(sourcearrayd(:,:,:,:),kind=CUSTOM_REAL)

      else ! use of CMTSOLUTION files

        call compute_arrays_source(sourcearray,xi,eta,gamma, &
                          Mxx(isource),Myy(isource),Mzz(isource),Mxy(isource), &
                          Mxz(isource),Myz(isource), &
                          xix_crust_mantle(:,:,:,ispec), &
                          xiy_crust_mantle(:,:,:,ispec), &
                          xiz_crust_mantle(:,:,:,ispec), &
                          etax_crust_mantle(:,:,:,ispec), &
                          etay_crust_mantle(:,:,:,ispec), &
                          etaz_crust_mantle(:,:,:,ispec), &
                          gammax_crust_mantle(:,:,:,ispec), &
                          gammay_crust_mantle(:,:,:,ispec), &
                          gammaz_crust_mantle(:,:,:,ispec), &
                          xigll,yigll,zigll)

      endif

      ! stores source excitations
      sourcearrays(:,:,:,:,isource) = sourcearray(:,:,:,:)

    endif
  enddo

  end subroutine setup_sources_receivers_srcarr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources_receivers_adjindx(NSTEP,NSTEP_SUB_ADJ, &
                                             NTSTEP_BETWEEN_READ_ADJSRC, &
                                             iadjsrc,iadjsrc_len,iadj_vec)

  use constants

  implicit none

  integer :: NSTEP,NSTEP_SUB_ADJ,NTSTEP_BETWEEN_READ_ADJSRC

  integer, dimension(NSTEP_SUB_ADJ,2) :: iadjsrc ! to read input in chunks
  integer, dimension(NSTEP_SUB_ADJ) :: iadjsrc_len
  integer, dimension(NSTEP) :: iadj_vec

  ! local parameters
  integer :: iadj_block,it,it_sub_adj
  integer :: istart,iend

  ! initializes
  iadjsrc(:,:) = 0
  iadjsrc_len(:) = 0

  ! setting up chunks of NTSTEP_BETWEEN_READ_ADJSRC to read adjoint source traces
  ! i.e. as an example: total length NSTEP = 3000, chunk length NTSTEP_BETWEEN_READ_ADJSRC= 1000
  !                                then it will set first block from 2001 to 3000,
  !                                second block from 1001 to 2000 and so on...
  !
  ! see routine: compute_arrays_source_adjoint()
  !                     how we read in the adjoint source trace in blocks/chunk sizes
  !
  ! see routine: compute_add_sources_adjoint()
  !                     how the adjoint source is added to the (adjoint) acceleration field
  !counts blocks
  ! block number
  ! e.g. increases from 1 (case it=1-1000), 2 (case it=1001-2000) to 3 (case it=2001-3000)
  it_sub_adj = 0
  iadj_block = 1
  do it = 1,NSTEP
    ! we are at the edge of a block
    if (mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC) == 0) then
      ! sets it_sub_adj subset number
      it_sub_adj = iadj_block

      ! block start time ( e.g. 2001)
      istart = NSTEP-it_sub_adj*NTSTEP_BETWEEN_READ_ADJSRC+1
      ! final adj src array
      ! e.g. will be from 1000 to 1, but doesn't go below 1 in cases where NSTEP isn't
      ! a multiple of NTSTEP_BETWEEN_READ_ADJSRC
      if (istart < 1 ) istart = 1

      ! block end time (e.g. 3000)
      iend = NSTEP-(it_sub_adj-1)*NTSTEP_BETWEEN_READ_ADJSRC

      iadjsrc(iadj_block,1) = istart
      iadjsrc(iadj_block,2) = iend

      ! actual block length
      iadjsrc_len(iadj_block) = iend - istart + 1

      ! increases block number
      iadj_block = iadj_block + 1
    endif

    ! time stepping for adjoint sources:
    ! adjoint time step that corresponds to time step in simulation (it).
    ! note, that adjoint source has to be time-reversed with respect to the forward wavefield
    ! e.g.: first block 1 has iadjsrc_len = 1000 with start at 2001 and end at 3000
    !         so iadj_vec(1) = 1000 - 0, iadj_vec(2) = 1000 - 1, ..., to iadj_vec(1000) = 1000 - 999 = 1
    !         then for block 2, iadjsrc_len = 1000 with start at 1001 and end at 2000
    !         so iadj_vec(1001) = 1000 - 0, iadj_vec(1002) = 1000 - 1, .. and so on again down to 1
    !         then block 3 and your guess is right now... iadj_vec(2001) to iadj_vec(3000) is 1000 down to 1. :)
    iadj_vec(it) = iadjsrc_len(it_sub_adj) - mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC)

    ! checks that index is non-negative
    if (iadj_vec(it) < 1 ) iadj_vec(it) = 1
  enddo

  end subroutine setup_sources_receivers_adjindx

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_receivers_precompute_intp()

  use specfem_par
  implicit none

  ! local parameters
  integer :: ier
  integer :: nadj_hprec_local

  double precision, dimension(NGLLX) :: hxir,hpxir
  double precision, dimension(NGLLY) :: hpetar,hetar
  double precision, dimension(NGLLZ) :: hgammar,hpgammar

  integer :: i,j,k,irec,irec_local
  real(kind=CUSTOM_REAL) :: hxi,heta,hgamma

  ! note: for adjoint simulations (SIMULATION_TYPE == 2),
  !         nrec_local     - is set to the number of sources (CMTSOLUTIONs), which act as "receiver" locations
  !                          for storing seismograms or strains
  !
  !         nadj_rec_local - determines the number of adjoint sources, i.e., number of station locations (STATIONS_ADJOINT), which
  !                          act as sources to drive the adjoint wavefield

  ! define local to global receiver numbering mapping
  ! needs to be allocated for subroutine calls (even if nrec_local == 0)
  allocate(number_receiver_global(nrec_local),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating global receiver numbering')
  number_receiver_global(:) = 0

  ! receivers
  ! allocates receiver interpolators
  if (nrec_local > 0) then
    ! allocates Lagrange interpolators for receivers
    allocate(hxir_store(NGLLX,nrec_local), &
             hetar_store(NGLLY,nrec_local), &
             hgammar_store(NGLLZ,nrec_local),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating receiver interpolators')

    ! defines and stores Lagrange interpolators at all the receivers
    if (SIMULATION_TYPE == 2) then
      nadj_hprec_local = nrec_local
    else
      nadj_hprec_local = 1
    endif
    allocate(hpxir_store(NGLLX,nadj_hprec_local), &
             hpetar_store(NGLLY,nadj_hprec_local), &
             hpgammar_store(NGLLZ,nadj_hprec_local),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating derivative interpolators')

    ! stores interpolators for receiver positions
    call setup_sources_receivers_intp(NSOURCES, &
                                      islice_selected_source, &
                                      xi_source,eta_source,gamma_source, &
                                      xigll,yigll,zigll, &
                                      SIMULATION_TYPE,nrec,nrec_local, &
                                      islice_selected_rec,number_receiver_global, &
                                      xi_receiver,eta_receiver,gamma_receiver, &
                                      hxir_store,hetar_store,hgammar_store, &
                                      nadj_hprec_local,hpxir_store,hpetar_store,hpgammar_store)

    ! allocates seismogram array
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      allocate(seismograms(NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS),stat=ier)
      if (ier /= 0) stop 'Error while allocating seismograms'
    else
      ! adjoint seismograms
      allocate(seismograms(NDIM*NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS),stat=ier)
      if (ier /= 0) stop 'Error while allocating adjoint seismograms'

      ! allocates Frechet derivatives array
      allocate(moment_der(NDIM,NDIM,nrec_local), &
               sloc_der(NDIM,nrec_local), &
               stshift_der(nrec_local), &
               shdur_der(nrec_local),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating Frechet derivatives arrays')

      moment_der(:,:,:) = 0._CUSTOM_REAL
      sloc_der(:,:) = 0._CUSTOM_REAL
      stshift_der(:) = 0._CUSTOM_REAL
      shdur_der(:) = 0._CUSTOM_REAL
    endif
    ! initializes seismograms
    seismograms(:,:,:) = 0._CUSTOM_REAL
    ! adjoint seismograms
    it_adj_written = 0
  else
    ! dummy arrays
    ! allocates dummy array since we need it to pass as argument e.g. in write_seismograms() routine
    ! note: nrec_local is zero, Fortran 90/95 should allow zero-sized array allocation...
    allocate(seismograms(NDIM,0,NTSTEP_BETWEEN_OUTPUT_SEISMOS),stat=ier)
    if (ier /= 0) stop 'Error while allocating zero seismograms'
    ! dummy allocation
    allocate(hxir_store(1,1), &
             hetar_store(1,1), &
             hgammar_store(1,1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating dummy receiver interpolators')
  endif

  ! strain seismograms
  if (SAVE_SEISMOGRAMS_STRAIN) then
    if (nrec_local > 0) then
      allocate(seismograms_eps(6,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS),stat=ier)
      if (ier /= 0) stop 'Error while allocating strain seismograms'
      seismograms_eps(:,:,:) = 0._CUSTOM_REAL
    else
      ! dummy
      allocate(seismograms_eps(1,1,1))
    endif
  endif

  ! adjoint sources
  ! optimizing arrays for adjoint sources
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! local adjoint sources arrays
    if (nadj_rec_local > 0) then
      ! determines adjoint sources arrays
      if (SIMULATION_TYPE == 2) then
        ! pure adjoint simulations
        allocate(number_adjsources_global(nadj_rec_local),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating number_adjsources_global array')
        number_adjsources_global(:) = 0

        ! addressing from local to global receiver index
        irec_local = 0
        do irec = 1,nrec
          ! add the source (only if this proc carries the source)
          if (myrank == islice_selected_rec(irec)) then
            irec_local = irec_local + 1
            number_adjsources_global(irec_local) = irec
          endif
        enddo
        if (irec_local /= nadj_rec_local) stop 'Error invalid number of nadj_rec_local found'

        ! allocate Lagrange interpolators for adjoint sources
        !
        ! note: adjoint sources for SIMULATION_TYPE == 2 and 3 are located at the receivers,
        !       however, the interpolator arrays hxir_store are used for "receiver" locations which are different
        !       for pure adjoint or kernel simulations
        !
        !       we will thus allocate interpolator arrays especially for adjoint source locations. for kernel simulations,
        !       these would be the same as hxir_store, but not for pure adjoint simulations.
        !
        ! storing these arrays is cheaper than storing a full (i,j,k) array for each element
        allocate(hxir_adjstore(NGLLX,nadj_rec_local),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating hxir_adjstore array')
        allocate(hetar_adjstore(NGLLY,nadj_rec_local),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating hetar_adjstore array')
        allocate(hgammar_adjstore(NGLLZ,nadj_rec_local),stat=ier)
        if (ier /= 0) call exit_MPI_without_rank('error allocating hgammar_adjstore array')

        ! define and store Lagrange interpolators at all the adjoint source locations
        do irec_local = 1,nadj_rec_local
          irec = number_adjsources_global(irec_local)

          ! receiver positions (become adjoint source locations)
          call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
          call lagrange_any(eta_receiver(irec),NGLLY,yigll,hetar,hpetar)
          call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)

          ! stores interpolators
          hxir_adjstore(:,irec_local) = hxir(:)
          hetar_adjstore(:,irec_local) = hetar(:)
          hgammar_adjstore(:,irec_local) = hgammar(:)
        enddo
      else
        ! kernel simulations (SIMULATION_TYPE == 3)
        ! adjoint source arrays and receiver arrays are the same, no need to allocate new arrays, just point to the existing ones
        number_adjsources_global => number_receiver_global
        hxir_adjstore => hxir_store
        hetar_adjstore => hetar_store
        hgammar_adjstore => hgammar_store
      endif
    else
      ! dummy arrays
      number_adjsources_global => number_receiver_global
      hxir_adjstore => hxir_store
      hetar_adjstore => hetar_store
      hgammar_adjstore => hgammar_store
    endif ! nadj_rec_local
  endif

  ! safety check
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! adjoint source in this partitions
    if (nadj_rec_local > 0) then
      do irec_local = 1, nadj_rec_local
        irec = number_adjsources_global(irec_local)
        if (irec <= 0) stop 'Error invalid irec for local adjoint source'
        ! adds source array
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              hxi = hxir_adjstore(i,irec_local)
              heta = hetar_adjstore(j,irec_local)
              hgamma = hgammar_adjstore(k,irec_local)
              ! checks if array values valid
              ! Lagrange interpolators shoud be about in a range ~ [-0.2,1.2]
              if (abs(hxi) > 2.0 .or. abs(heta) > 2.0 .or. abs(hgamma) > 2.0) then
                print *,'hxi/heta/hgamma = ',hxi,heta,hgamma,irec_local,i,j,k
                print *,'ERROR: trying to use arrays hxir_adjstore/hetar_adjstore/hgammar_adjstore with irec_local = ', &
                        irec_local,' but these array values are invalid!'
                call exit_MPI_without_rank('ERROR: trying to use arrays hxir_adjstore/hetar_adjstore/hgammar_adjstore &
                                           &but these arrays are invalid!')
              endif
            enddo
          enddo
        enddo
      enddo
    endif
  endif

  ! ASDF seismograms
  if (OUTPUT_SEISMOS_ASDF) then
    if (.not. (SIMULATION_TYPE == 3 .and. (.not. SAVE_SEISMOGRAMS_IN_ADJOINT_RUN)) ) then
      ! initializes the ASDF data structure by allocating arrays
      call init_asdf_data(nrec_local)
      call synchronize_all()
    endif
  endif

  end subroutine setup_receivers_precompute_intp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_sources_receivers_intp(NSOURCES, &
                      islice_selected_source, &
                      xi_source,eta_source,gamma_source, &
                      xigll,yigll,zigll, &
                      SIMULATION_TYPE,nrec,nrec_local, &
                      islice_selected_rec,number_receiver_global, &
                      xi_receiver,eta_receiver,gamma_receiver, &
                      hxir_store,hetar_store,hgammar_store, &
                      nadj_hprec_local,hpxir_store,hpetar_store,hpgammar_store)

  use constants

  implicit none

  integer :: NSOURCES

  integer, dimension(NSOURCES) :: islice_selected_source

  double precision, dimension(NSOURCES) :: xi_source,eta_source,gamma_source
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll


  integer :: SIMULATION_TYPE

  integer :: nrec,nrec_local
  integer, dimension(nrec) :: islice_selected_rec
  integer, dimension(nrec_local) :: number_receiver_global
  double precision, dimension(nrec) :: xi_receiver,eta_receiver,gamma_receiver

  double precision, dimension(NGLLX,nrec_local) :: hxir_store
  double precision, dimension(NGLLY,nrec_local) :: hetar_store
  double precision, dimension(NGLLZ,nrec_local) :: hgammar_store

  integer :: nadj_hprec_local
  double precision, dimension(NGLLX,nadj_hprec_local) :: hpxir_store
  double precision, dimension(NGLLY,nadj_hprec_local) :: hpetar_store
  double precision, dimension(NGLLZ,nadj_hprec_local) :: hpgammar_store


  ! local parameters
  integer :: isource,irec,irec_local
  double precision, dimension(NGLLX) :: hxir,hpxir
  double precision, dimension(NGLLY) :: hpetar,hetar
  double precision, dimension(NGLLZ) :: hgammar,hpgammar

  ! select local receivers

  ! define local to global receiver numbering mapping
  irec_local = 0
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    ! forward/kernel simulations
    do irec = 1,nrec
      if (myrank == islice_selected_rec(irec)) then
        irec_local = irec_local + 1
        ! checks counter
        if (irec_local > nrec_local) call exit_MPI(myrank,'Error receiver interpolators: irec_local exceeds bounds')
        ! stores local to global receiver ids
        number_receiver_global(irec_local) = irec
      endif
    enddo
  else
    ! adjoint simulations
    do isource = 1,NSOURCES
      if (myrank == islice_selected_source(isource)) then
        irec_local = irec_local + 1
        ! checks counter
        if (irec_local > nrec_local) call exit_MPI(myrank,'Error adjoint source interpolators: irec_local exceeds bounds')
        ! stores local to global receiver/source ids
        number_receiver_global(irec_local) = isource
      endif
    enddo
  endif
  ! checks if all local receivers have been found
  if (irec_local /= nrec_local) call exit_MPI(myrank,'Error number of local receivers do not match')

  ! define and store Lagrange interpolators at all the receivers
  do irec_local = 1,nrec_local
    irec = number_receiver_global(irec_local)

    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      ! receiver positions
      call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
      call lagrange_any(eta_receiver(irec),NGLLY,yigll,hetar,hpetar)
      call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
    else
      ! source positions
      call lagrange_any(xi_source(irec),NGLLX,xigll,hxir,hpxir)
      call lagrange_any(eta_source(irec),NGLLY,yigll,hetar,hpetar)
      call lagrange_any(gamma_source(irec),NGLLZ,zigll,hgammar,hpgammar)
    endif

    ! stores interpolators
    hxir_store(:,irec_local) = hxir(:)
    hetar_store(:,irec_local) = hetar(:)
    hgammar_store(:,irec_local) = hgammar(:)

    ! stores derivatives
    if (SIMULATION_TYPE == 2) then
      hpxir_store(:,irec_local) = hpxir(:)
      hpetar_store(:,irec_local) = hpetar(:)
      hpgammar_store(:,irec_local) = hpgammar(:)
    endif
  enddo

  end subroutine setup_sources_receivers_intp

