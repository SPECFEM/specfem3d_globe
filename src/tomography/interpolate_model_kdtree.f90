!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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

!------------------------------------------------------------------------------
!
! kdtree - nearest neighbor search
!
!------------------------------------------------------------------------------

module kdtree_search

  ! kd-tree for searching nearest neighbors

  private

  ! single tree node
  type :: kdtree_node
    ! id (used for debugging tree)
    !integer :: id
    ! cut dimension index ( 1/2/3 for 3D point set )
    integer :: idim
    ! cut value
    double precision :: cut_value
    ! bounds value on cut dimension
    !double precision :: cut_min,cut_max
    ! index of associated data point
    integer :: ipoint
    ! index of lower and upper bounds of point range
    integer :: ibound_lower,ibound_upper
    ! child nodes in sub level
    type (kdtree_node),pointer :: left, right
  end type kdtree_node

  ! data associated kd-tree structure (root node)
  ! note: in general, a save attribute should be added to this root pointer
  !       to have the tree variable stored globally
  !       however, putting it here into the module declaration will have the same effect
  type (kdtree_node),pointer :: kdtree

  ! tree search arrays
  ! total number of tree nodes
  integer :: kdtree_nnodes_local

  ! tree node locations
  double precision, dimension(:,:),allocatable,target :: kdtree_nodes_local

  ! associated node index
  integer, dimension(:),allocatable :: kdtree_index_local

  ! info output
  logical :: be_verbose = .false.

  ! public routines
  public :: kdtree_setup,kdtree_find_nearest_neighbor,kdtree_delete
  public :: kdtree_set_verbose

  ! public parameters/arrays
  public :: kdtree_nnodes_local
  public :: kdtree_nodes_local
  public :: kdtree_index_local

contains

! example:
!
! creates kd-tree for searching
!  .. prepare point array kdtree_nodes_local,kdtree_nnodes_local
!  call kdtree_setup()
!
! finds closest point
!  .. do-loop points at xyz
!     call kdtree_find_nearest_neighbor(xyz,iglob_min,dist_min)
!  .. enddo
!
! deletes search tree
!  call kdtree_delete()


  subroutine kdtree_setup()

  ! sets up the kd-tree structure
  !
  ! needs:
  !   kdtree_nnodes_local  - total number of points
  !   kdtree_nodes_local   - 3D array of points
  !
  ! returns:
  !   creates internal tree representation
  !
  !
  ! example tree creation timing = 0.008 s (for data points = 38704 -> tree nodes = 77407 )

  implicit none

  ! local parameters
  integer :: npoints
  integer, dimension(:),allocatable :: points_index
  double precision,dimension(:,:),pointer :: points_data

  ! tree statistics
  integer :: depth
  integer :: numnodes,maxdepth
  integer :: i,ier

  ! timing
  real :: ct_start,ct_end

  ! test search
  double precision, dimension(3) :: xyz_target
  integer :: ipoint_min
  double precision :: dist_min

  !------------------------------------------------------

  ! debugging: performs a test search
  logical,parameter :: DEBUG = .false.

  !------------------------------------------------------

  ! checks
  if (kdtree_nnodes_local <= 0 ) stop 'Error creating kdtree with zero nodes is invalid'
  if (.not. allocated(kdtree_nodes_local) ) stop 'Error array kdtree_nodes_local not allocated yet'

  ! timing
  call cpu_time(ct_start)

  ! number of data points
  npoints = kdtree_nnodes_local

  ! 3D point coordinates
  points_data => kdtree_nodes_local(:,:)

  if (be_verbose) then
    print*,'kd-tree:'
    print*,'  total data points: ',npoints
    !print*,'  box boundaries   : x min/max = ',minval(points_data(1,:)),maxval(points_data(1,:))
    !print*,'                     y min/max = ',minval(points_data(2,:)),maxval(points_data(2,:))
    !print*,'                     z min/max = ',minval(points_data(3,:)),maxval(points_data(3,:))
  endif

  ! theoretical number of node for totally balanced tree
  numnodes = npoints
  i = npoints
  do while ( i >= 1 )
    i = i / 2
    numnodes = numnodes + i
    ! integer 32-bit limitation
    if (numnodes > 2147483646 - i ) stop 'Error number of nodes might exceed integer limit'
  enddo
  if (be_verbose) then
    print*,'  theoretical number of nodes: ',numnodes
  endif

  ! local ordering
  allocate(points_index(kdtree_nnodes_local),stat=ier)
  if (ier /= 0) stop 'Error allocating array points_index'

  ! initial point ordering
  do i = 1,npoints
    points_index(i) = i
  enddo

  ! builds tree structure
  nullify(kdtree)
  depth = 0
  numnodes = 0
  maxdepth = -1

  call create_kdtree(npoints,points_data,points_index,kdtree, &
                     depth,1,npoints,numnodes,maxdepth)

  ! checks root tree node
  if (.not. associated(kdtree) ) stop 'Error creation of kd-tree failed'

  if (be_verbose) then
    print*,'  actual      number of nodes: ',numnodes
    ! tree node size: 4 (idim) + 8 (cut_value) + 4 (ipoint) + 2*4 (ibound_**) + 2*4 (left,right) = 32 bytes
    print*,'  tree memory size: ', ( numnodes * 32 )/1024./1024.,'MB'
    print*,'  maximum depth   : ',maxdepth

    ! timing
    call cpu_time(ct_end)
    print*,'  creation timing : ',ct_end - ct_start, '(s)'
    print*
  endif


  ! debugging
  if (DEBUG) then
    ! outputs tree
    if (be_verbose) then
      numnodes = 0
      call print_kdtree(npoints,points_data,points_index,kdtree,numnodes)
    endif

    ! test search
    print*,'search tree:'
    xyz_target(1) = 0.13261298835277557
    xyz_target(2) = -8.4083788096904755E-002
    xyz_target(3) = 0.97641450166702271

    print*,'search : ',xyz_target(:)

    ipoint_min = -1
    dist_min = 1.d30

    call find_nearest_kdtree_node(npoints,points_data,kdtree, &
                                  xyz_target,ipoint_min,dist_min)

    dist_min = sqrt(dist_min)
    print*,'found : ',ipoint_min,'distance:',dist_min

    if (ipoint_min < 1 ) stop 'Error search kd-tree found no point'

    print*,'target  : ',xyz_target(:)
    print*,'nearest : ',points_data(:,ipoint_min)
    print*
    ! safety stop
    stop 'kdtree_setup safety stop'
  endif

  ! frees temporary arrays
  deallocate(points_index)
  nullify(points_data)

  end subroutine kdtree_setup

!===================================================================================================

  subroutine kdtree_find_nearest_neighbor(xyz_target,iglob_min,dist_min)

  ! kd-tree nearest neighbor search
  !
  ! input:
  !
  ! returns: global index iglob_min and distance dist_min to nearest point

  implicit none

  double precision, dimension(3),intent(in) :: xyz_target

  double precision, intent(out) :: dist_min
  integer, intent(out) :: iglob_min

  ! local parameters
  integer :: ipoint_min

  ! initializes
  ipoint_min = -1
  iglob_min = -1
  dist_min = 1.d30

  ! searches closest node in kd-tree
  call find_nearest_kdtree_node(kdtree_nnodes_local,kdtree_nodes_local,kdtree, &
                                xyz_target,ipoint_min,dist_min)

  if (ipoint_min < 1 .or. ipoint_min > kdtree_nnodes_local ) stop 'Error search kd-tree found no point'

  ! gets global index
  iglob_min = kdtree_index_local(ipoint_min)

  ! checks global index
  if (iglob_min < 1 ) stop 'Error minimum location has wrong global index in kdtree_find_nearest_neighbor'

  ! returns distance (non-squared)
  dist_min = sqrt( dist_min )

  ! debug
  !if (be_verbose) then
  !  print*,'target  : ',xyz_target(:)
  !  print*,'nearest : ',kdtree_nodes_local(:,ipoint_min),'distance:',dist_min*6371.,'(km)',ipoint_min,iglob_min
  !endif

  end subroutine kdtree_find_nearest_neighbor

!===================================================================================================

  recursive subroutine create_kdtree(npoints,points_data,points_index,node, &
                                     depth,ibound_lower,ibound_upper,numnodes,maxdepth)

  ! creates node in kd-tree structure

  implicit none

  integer,intent(in) :: npoints
  double precision,dimension(3,npoints),intent(in) :: points_data
  integer,dimension(npoints),intent(inout) :: points_index

  type (kdtree_node), pointer,intent(inout) :: node

  integer,intent(in) :: depth
  integer,intent(in) :: ibound_lower,ibound_upper

  integer,intent(inout) :: numnodes,maxdepth

  ! local parameters
  double precision :: cut_value
  double precision :: range,range_max,min,max
  integer :: i,ier,idim
  integer :: iloc,ilower,iupper
  integer :: l,u
  integer,dimension(:),allocatable :: workindex

  ! checks if anything to sort
  if (ibound_lower > ibound_upper ) then
    nullify(node)
    return
  endif

  ! creates new node
  allocate(node,stat=ier)
  if (ier /= 0) stop 'Error allocating kd-tree node'

  ! initializes new node
  node%idim = -1
  node%ipoint = -1
  node%cut_value = 0.d0

  nullify(node%left)
  nullify(node%right)

  ! tree statistics
  numnodes = numnodes + 1
  if (maxdepth < depth) maxdepth = depth

  !node%id = numnodes

  ! checks if final node
  if (ibound_lower == ibound_upper ) then
    node%idim = 0
    node%ipoint = points_index(ibound_lower)
    ! done with this node
    return
  endif

  ! sets cut dimension index

  ! version 1: varies between 1 and 3 depending on depth for 3D data set
  ! (leads to some unneccessary unbalanced nodes)
  !idim = mod(depth,3) + 1
  ! determines cut value
  ! range in this dimension
  !min = HUGEVAL
  !max = - HUGEVAL
  !do i = ibound_lower,ibound_upper
  !  iloc = points_index(i)
  !  val = points_data(idim,iloc)
  !  if (val < min ) min = val
  !  if (val > max ) max = val
  !enddo
  !min = minval(points_data(idim,points_index(ibound_lower:ibound_upper)))
  !max = maxval(points_data(idim,points_index(ibound_lower:ibound_upper)))
  !cut_value = 0.5d0 * ( min + max )

  ! version 2: selects cut dimension where biggest range is occurring
  ! (better balances tree)
  cut_value = 0.d0
  idim = -1
  range_max = 0.d0
  do i = 1,3
    min = minval(points_data(i,points_index(ibound_lower:ibound_upper)))
    max = maxval(points_data(i,points_index(ibound_lower:ibound_upper)))
    range = max - min
    ! sets cut dimension where data has maximum range
    if (range > range_max ) then
      range_max = range
      idim = i
      cut_value = 0.5d0 * ( min + max )
      ! stores bounds of cut dimension
      !node%cut_min = min
      !node%cut_max = max
    endif
  enddo
  node%idim = idim
  node%cut_value = cut_value

  !debug
  !print*,'index ',numnodes,'dim:',idim,'range:',ibound_lower,ibound_upper
  !print*,'  data:',points_data(idim,points_index(ibound_lower)),points_data(idim,points_index(ibound_upper))
  !print*,'  min/max:',min,max,'cut value:',cut_value

  ! temporary index array for sorting
  allocate(workindex(ibound_upper - ibound_lower + 1),stat=ier)
  if (ier /= 0) stop 'Error allocating workindex array'

  ! sorts point indices
  ! to have all points with value < cut_value on left side, all others to the right
  ilower = 0
  iupper = 0
  do i = ibound_lower,ibound_upper
    iloc = points_index(i)
    if (points_data(idim,iloc) < cut_value ) then
      ilower = ilower + 1
      workindex(ilower) = iloc
    else
      iupper = iupper + 1
      workindex(ibound_upper - ibound_lower + 2 - iupper) = iloc
    endif
  enddo
  !debug
  !print*,'  ilower/iupper:',ilower,iupper

  ! checks if we catched all
  if (ilower + iupper /= ibound_upper - ibound_lower + 1 ) stop 'Error sorting data points invalid'
  if (ilower == 0 .and. iupper == 0 .and. npoints > 1 ) stop 'Error confusing node counts, please check kdtree...'

  ! replaces index range with new sorting order
  points_index(ibound_lower:ibound_upper) = workindex(:)

  ! frees temporary array
  deallocate(workindex)

  ! lower hemisphere
  if (ilower > 0) then
    ! lower/upper bounds
    l = ibound_lower
    u = ibound_lower + (ilower - 1)
    if (l < 1 .or. u > npoints ) stop 'Error lower hemisphere tree bounds'
    ! adds new node
    call create_kdtree(npoints,points_data,points_index,node%left,depth+1,l,u,numnodes,maxdepth)
  endif

  ! upper hemisphere points
  if (iupper > 0) then
    ! lower/upper bounds
    l = ibound_upper - (iupper - 1)
    u = ibound_upper
    if (l < 1 .or. u > npoints ) stop 'Error upper hemisphere tree bounds'
    ! adds new node
    call create_kdtree(npoints,points_data,points_index,node%right,depth+1,l,u,numnodes,maxdepth)
  endif

  end subroutine create_kdtree

!===================================================================================================


  recursive subroutine print_kdtree(npoints,points_data,points_index,node,numnodes)

  ! prints out all final nodes in kd-tree structure

  implicit none

  integer,intent(in) :: npoints
  double precision,dimension(3,npoints),intent(in) :: points_data
  integer,dimension(npoints),intent(in) :: points_index

  type (kdtree_node), pointer,intent(inout) :: node

  integer,intent(inout) :: numnodes

  ! local parameters
  integer,parameter :: OUTPUT_LENGTH = 50

  ! checks if valid pointer (must have been nullified initially to be able to check with associated())
  if (.not. associated(node) ) return

  ! statistics
  numnodes = numnodes + 1
  if (numnodes == 1) then
    print*,'printing kd-tree: total number of points      = ',npoints
    !print*,'         index array = ',points_index(:)
  endif

  ! outputs infos for a final node
  if (.not. associated(node%left) .and. .not. associated(node%right) ) then
    ! checks info
    if (node%idim /= 0 ) then
      print*,'problem kd-tree node:',node%idim,node%ipoint,numnodes
      print*,'point x/y/z: ',points_data(:,node%ipoint)
      stop 'Error kd-tree node not correct'
    endif

    ! outputs infos
    if (numnodes < OUTPUT_LENGTH) &
      print*,'node:',numnodes,'index:',node%ipoint,' x/y/z = ',points_data(:,node%ipoint)
  else
    ! outputs infos
    if (numnodes < OUTPUT_LENGTH) &
      print*,'node:',numnodes,'dim:',node%idim,'cut = ',node%cut_value
  endif

  ! checks child nodes
  if (associated(node%left) ) then
    call print_kdtree(npoints,points_data,points_index,node%left,numnodes)
  endif
  if (associated(node%right) ) then
    call print_kdtree(npoints,points_data,points_index,node%right,numnodes)
  endif

  end subroutine print_kdtree

!===================================================================================================


  subroutine kdtree_delete()

  ! deletes all (child) nodes in this given tree

  implicit none

  ! deletes tree starting with top node
  call delete_node(kdtree)

  end subroutine kdtree_delete

!===================================================================================================


  recursive subroutine delete_node(node)

  ! deletes all (child) nodes in this given tree

  implicit none

  type (kdtree_node), pointer :: node

  ! delete left hemisphere
  if (associated(node%left)) then
     call delete_node(node%left)
     nullify (node%left)
  end if

  ! deletes right hemisphere
  if (associated(node%right)) then
     call delete_node(node%right)
     nullify (node%right)
  end if

  deallocate(node)

  end subroutine delete_node

!===================================================================================================

  recursive subroutine find_nearest_kdtree_node(npoints,points_data,node, &
                                                xyz_target,ipoint_min,dist_min)

  ! searches for node point closest to given location
  implicit none
  integer,intent(in) :: npoints
  double precision,dimension(3,npoints),intent(in) :: points_data

  type (kdtree_node), pointer,intent(inout) :: node

  double precision,dimension(3),intent(in) :: xyz_target

  integer,intent(inout) :: ipoint_min
  double precision,intent(inout) :: dist_min

  ! local parameters
  double precision :: dist
  double precision,dimension(3) :: xyz

  ! debug
  !if (node%idim == 0) then
  !  print*,'node',node%id,points_data(:,node%ipoint)
  !else
  !  print*,'node',node%id,node%idim,node%cut_value
  !endif
  !if (ipoint_min > 0) &
  !  print*,'node distance',node%id,ipoint_min,dist_min

  ! in case this is a final node
  if ( .not. associated(node%left) .and. .not. associated(node%right) ) then
    ! checks node
    if (node%idim /= 0 ) stop 'Error searched node is not final node'
    if (node%ipoint < 1 ) stop 'Error searched node has wrong point index'

    ! squared distance to associated data point
    xyz(:) = xyz_target(:) - points_data(:,node%ipoint)
    dist = xyz(1) * xyz(1) + xyz(2)*xyz(2) + xyz(3)*xyz(3)
    if (dist < dist_min) then
      ! debug
      !if (ipoint_min < 1 ) then
      !  print*,'new node distance',node%id,node%ipoint,dist
      !else
      !  print*,'     new distance',node%id,node%ipoint,dist
      !endif
      ! stores minimum point
      dist_min = dist
      ipoint_min = node%ipoint
    endif

    ! done
    return
  endif

  ! checks cut dimension
  if (node%idim < 1 .or. node%idim > 3 ) stop 'Error searched node has invalid cut dimension'

  ! compares cut value
  if (xyz_target(node%idim) < node%cut_value ) then
    ! finds closer node in lower hemisphere
    if (associated(node%left) ) then
      call find_nearest_kdtree_node(npoints,points_data,node%left, &
                                    xyz_target,ipoint_min,dist_min)
    endif
  else
    ! finds closer node in upper hemisphere
    if (associated(node%right) ) then
      call find_nearest_kdtree_node(npoints,points_data,node%right, &
                                    xyz_target,ipoint_min,dist_min)
    endif
  endif

  ! at this point, dist_min is the distance to the closest point in the initial hemisphere search
  ! we might need to search in other hemisphere as well if distances are closer

  ! squared distance to cut plane
  dist = ( xyz_target(node%idim) - node%cut_value )**2

  if (xyz_target(node%idim) < node%cut_value ) then
    if (associated(node%right) ) then
      ! checks right node as a final node
      if (node%right%idim == 0 ) then
        dist = sum((xyz_target(:) - points_data(:,node%right%ipoint))**2)
        if (dist <= dist_min) then
          ! stores minimum point
          dist_min = dist
          ipoint_min = node%right%ipoint
          return
        endif
      endif
      ! checks if points beyond cut plane could be closer
      if (dist < dist_min) then
        call find_nearest_kdtree_node(npoints,points_data,node%right, &
                                        xyz_target,ipoint_min,dist_min)
      endif
    endif
  else
    if (associated(node%left) ) then
      ! checks left node as a final node
      if (node%left%idim == 0 ) then
        dist = sum((xyz_target(:) - points_data(:,node%left%ipoint))**2)
        if (dist <= dist_min) then
          ! stores minimum point
          dist_min = dist
          ipoint_min = node%left%ipoint
          return
        endif
      endif
      ! checks if points beyond cut plane could be closer
      if (dist < dist_min) then
        call find_nearest_kdtree_node(npoints,points_data,node%left, &
                                    xyz_target,ipoint_min,dist_min)
      endif
    endif
  endif

  end subroutine find_nearest_kdtree_node


!===================================================================================================

  subroutine kdtree_set_verbose()

  implicit none

  ! sets verbosity on
  be_verbose = .true.

  end subroutine kdtree_set_verbose

end module