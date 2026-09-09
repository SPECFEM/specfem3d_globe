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

!----
!---- Source location against a Green function database.
!----
!---- (lat,lon,depth) -> (element, xi, eta, gamma, inverse Jacobian, nu),
!---- reproducing what the solver itself does. The pieces live one level
!---- down in gf_geometry.F90 and gf_shape3D.F90; this module supplies the
!---- database: which elements to try, in what order, and which answer to
!---- accept.
!----
!---- Element search
!---- --------------
!---- Candidates come from a kd-tree over the element centroids, built with
!---- src/shared/search_kdtree.f90 -- the solver's own tree, not a second
!---- implementation. The alternative the format invites, and that
!---- gf_cross_validate.py:find_containing_element takes, is to open every
!---- every element's own coordinates.h5 and read three attributes; that
!---- is 70,000
!---- file opens per source at production scale, which is exactly what
!---- centroids.bin exists to avoid.
!----
!---- kdtree_search is a **singleton**: one tree per process, held in module
!---- variables that the caller also owns. Stage 9 loads this library into a
!---- Python interpreter where two open databases are ordinary, and there
!---- the singleton would silently hand database B's elements to database
!---- A's queries. So the tree records which handle built it and is rebuilt
!---- when a different one asks -- O(n log n), milliseconds, against a 21 MB
!---- element read. Wrong answers become a rebuild.
!----
!---- On `stop`. Nothing in src/gf3d/ may stop, because a stop inside a
!---- shared object kills the interpreter that loaded it, and kdtree_search
!---- has thirty of them. Auditing the ones this module can reach, every one
!---- guards a precondition it is our job to meet:
!----   * kdtree_num_nodes > 0 and kdtree_nodes_location allocated (:163-164)
!----     -- gf_tree_ensure checks nelem before calling kdtree_setup;
!----   * a non-negative search radius (:356,420,462) -- always positive here;
!----   * kdtree_search_index allocated with a matching
!----     kdtree_search_num_nodes (:360-364,466-470) -- allocated from the
!----     count immediately before the get;
!----   * a valid point index (:305,311,381,496) -- follows from a non-empty
!----     tree over a validated element index.
!---- The verbose block at :245-266 is unreachable (be_verbose is .false.).
!---- What remains is allocation failure (:203,221,678), which is a stop we
!---- cannot take away without forking the module.
!----
!---- Containment
!---- -----------
!---- The first candidate with max(|xi|,|eta|,|gamma|) <= GF_XI_TOL wins.
!---- See the comment on GF_XI_TOL in gf_par.F90 for why that constant is
!---- 1.099 and why replicating the solver is the point. If no candidate
!---- passes, this is a hard error: find_containing_element warns and returns
!---- xi=eta=gamma=0 at the nearest centroid, which yields a plausible-looking
!---- and completely wrong seismogram.
!----

  module gf_locate

  use gf_par, only: t_gfdb,t_gf_location,gf_set_error,MAX_STRING_LEN, &
                    GF_OK,GF_ERR_ARG,GF_ERR_ALLOC,GF_ERR_MISMATCH, &
                    GF_ERR_NO_ELEMENT,GF_XI_TOL,GF_NCAND,GF_ANCHOR_TOL

  use gf_database, only: gf_load_topo,gf_topo_elevation

  use gf_element_io, only: gf_read_element_coords

  use gf_geometry, only: gf_geographic_to_cartesian,gf_source_nu, &
                         gf_gather_anchors,gf_find_local_coords

  use gf_shape3D, only: gf_shape3D_map,gf_shape3D_map_2nd

  implicit none

  private

  public :: gf_locate_source
  public :: gf_locate_release
  public :: gf_check_anchors
  public :: gf_check_anchors_all

  !-----------------------------------------------------------------
  ! kd-tree ownership
  !
  ! kdtree_search keeps its tree in module variables, so these track which
  ! database the live tree describes. Empty `tree_owner` means no tree.
  !-----------------------------------------------------------------

  character(len=MAX_STRING_LEN) :: tree_owner = ''
  logical :: tree_ready = .false.

  ! representative element size, non-dimensional, from the centroid cloud;
  ! only used to seed the search radius, which grows if it comes up empty
  double precision :: tree_typical_size = 0.d0

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_tree_ensure(db,ierr)

! makes sure the live kd-tree describes this database
!
! kdtree_nodes_index is filled with the identity. That is not cosmetic:
! kdtree_find_nearest_neighbor returns kdtree_nodes_index(ipoint) (:308)
! while the radius searches return the raw ipoint (:1114), so the two APIs
! only agree when the mapping is the identity. Keeping it that way lets both
! be used interchangeably below.

  use kdtree_search, only: kdtree_setup,kdtree_num_nodes, &
                           kdtree_nodes_location,kdtree_nodes_index

  implicit none

  type(t_gfdb), intent(in) :: db
  integer, intent(out) :: ierr

  ! local parameters
  integer :: i,ier
  double precision, dimension(3) :: cmin,cmax,ext
  double precision :: diag

  if (tree_ready) then
    if (trim(tree_owner) == trim(db%path)) then
      ierr = GF_OK
      return
    endif
  endif

  call gf_locate_release()

  ! kdtree_setup stops on an empty tree (search_kdtree.f90:163)
  if (db%nelem < 1) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_tree_ensure: database has no elements')
    return
  endif
  if (.not. allocated(db%centroid)) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_tree_ensure: database has no element centroids')
    return
  endif

  allocate(kdtree_nodes_location(3,db%nelem),kdtree_nodes_index(db%nelem),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate the element kd-tree')
    return
  endif

  kdtree_num_nodes = db%nelem
  do i = 1,db%nelem
    kdtree_nodes_location(1,i) = db%centroid(1,i)
    kdtree_nodes_location(2,i) = db%centroid(2,i)
    kdtree_nodes_location(3,i) = db%centroid(3,i)
    kdtree_nodes_index(i) = i
  enddo

  call kdtree_setup()

  ! a representative element size, used only to seed the search radius.
  ! The bounding-box form is deliberately crude but never degenerate: a
  ! database that is one element thick still gives a sane number, which
  ! volume**(1/3) would not.
  do i = 1,3
    cmin(i) = minval(kdtree_nodes_location(i,1:db%nelem))
    cmax(i) = maxval(kdtree_nodes_location(i,1:db%nelem))
    ext(i) = cmax(i) - cmin(i)
  enddo
  diag = sqrt(ext(1)*ext(1) + ext(2)*ext(2) + ext(3)*ext(3))
  tree_typical_size = diag / max(1.d0,dble(db%nelem)**(1.d0/3.d0))
  if (tree_typical_size <= 0.d0) tree_typical_size = 1.d-3

  tree_owner = db%path
  tree_ready = .true.

  ierr = GF_OK

  end subroutine gf_tree_ensure

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_locate_release()

! frees the kd-tree
!
! kdtree_delete() only walks the nodes; the two location/index arrays are
! the caller's, as setup_sources_receivers.f90:94 shows, so they are freed
! here too.
!
! Callers: gf3d_main before exit, and Stage 9's binding from its destructor.
! gf_close does *not* call this -- gf_locate depends on gf_database, so the
! reverse dependency would be circular. Leaving a tree behind is harmless:
! the next gf_tree_ensure for a different path rebuilds it.

  use kdtree_search, only: kdtree_delete,kdtree_num_nodes, &
                           kdtree_nodes_location,kdtree_nodes_index, &
                           kdtree_search_index,kdtree_search_num_nodes

  implicit none

  if (tree_ready) call kdtree_delete()

  if (allocated(kdtree_nodes_location)) deallocate(kdtree_nodes_location)
  if (allocated(kdtree_nodes_index)) deallocate(kdtree_nodes_index)
  if (allocated(kdtree_search_index)) deallocate(kdtree_search_index)

  kdtree_num_nodes = 0
  kdtree_search_num_nodes = 0

  tree_owner = ''
  tree_ready = .false.
  tree_typical_size = 0.d0

  end subroutine gf_locate_release

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_find_candidates(db,xyz,cand,nfound,ierr)

! the GF_NCAND elements whose centroids are nearest a point, nearest first
!
! Follows the count -> allocate -> get idiom of
! src/tomography/postprocess_sensitivity_kernels/create_cross_section.F90:1544.
! The radius searches are radius searches, not k-nearest, so the radius is
! seeded from the actual nearest-centroid distance and doubled while it comes
! back short.

  use kdtree_search, only: kdtree_find_nearest_neighbor, &
                           kdtree_count_nearest_n_neighbors, &
                           kdtree_get_nearest_n_neighbors, &
                           kdtree_nodes_location, &
                           kdtree_search_index,kdtree_search_num_nodes

  implicit none

  type(t_gfdb), intent(in) :: db
  double precision, dimension(3), intent(in) :: xyz
  integer, dimension(GF_NCAND), intent(out) :: cand
  integer, intent(out) :: nfound
  integer, intent(out) :: ierr

  ! how many times the radius may double before we give up on finding
  ! GF_NCAND neighbours and work with what we have
  integer, parameter :: NGROW = 6

  ! local parameters
  integer :: n,ngot,nwant,igrow,i,j,ibest,ier,inear
  double precision :: r,dist_min,dx,dy,dz,dbest
  double precision, dimension(:), allocatable :: dist
  integer, dimension(:), allocatable :: idx
  logical, dimension(:), allocatable :: taken

  cand(:) = 0
  nfound = 0

  call gf_tree_ensure(db,ierr)
  if (ierr /= GF_OK) return

  nwant = min(GF_NCAND,db%nelem)

  ! nearest centroid: gives a length scale without having to guess one
  call kdtree_find_nearest_neighbor(xyz,inear,dist_min)

  r = max(3.d0*dist_min,2.d0*tree_typical_size)
  if (r <= 0.d0) r = tree_typical_size

  n = 0
  do igrow = 1,NGROW
    call kdtree_count_nearest_n_neighbors(xyz,r,n)
    if (n >= nwant) exit
    r = 2.d0*r
  enddo

  if (n < 1) then
    ! the tree is non-empty and the radius has grown 2^NGROW times, so this
    ! means the point is nowhere near the database
    call gf_set_error(ierr,GF_ERR_NO_ELEMENT, &
      'no element centroid found near the requested position')
    return
  endif

  allocate(kdtree_search_index(n),dist(n),idx(n),taken(n),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate the kd-tree search result')
    return
  endif
  kdtree_search_num_nodes = n

  ngot = n
  call kdtree_get_nearest_n_neighbors(xyz,r,ngot)

  ! kdtree_search_index holds point indices, which are element indices here
  ! because gf_tree_ensure set kdtree_nodes_index to the identity. The module
  ! returns them sorted by *index*, not by distance, so order them properly.
  do i = 1,ngot
    idx(i) = kdtree_search_index(i)
    dx = xyz(1) - kdtree_nodes_location(1,idx(i))
    dy = xyz(2) - kdtree_nodes_location(2,idx(i))
    dz = xyz(3) - kdtree_nodes_location(3,idx(i))
    dist(i) = dx*dx + dy*dy + dz*dz
    taken(i) = .false.
  enddo

  ! selection sort of the first nwant only: ngot is a few dozen and nwant is
  ! ten, so this is cheaper than sorting and needs no scratch beyond `taken`
  nfound = 0
  do j = 1,min(nwant,ngot)
    ibest = 0
    dbest = huge(1.d0)
    do i = 1,ngot
      if (taken(i)) cycle
      if (dist(i) < dbest) then
        dbest = dist(i)
        ibest = i
      endif
    enddo
    if (ibest == 0) exit
    taken(ibest) = .true.
    nfound = nfound + 1
    cand(nfound) = idx(ibest)
  enddo

  deallocate(dist,idx,taken)
  if (allocated(kdtree_search_index)) deallocate(kdtree_search_index)
  kdtree_search_num_nodes = 0

  ierr = GF_OK

  end subroutine gf_find_candidates

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_locate_source(db,lat,lon,depth_km,loc,ierr)

! locates a geographic position in the database
!
! `db` is intent(inout) because the topography grid is loaded on demand:
! it is 58 MB and most of gf_open's callers never need it.

  use constants, only: NGNOD,NGLLX,NGLLY,NGLLZ,NDIM,GAUSSALPHA,GAUSSBETA

  implicit none

  type(t_gfdb), intent(inout) :: db
  double precision, intent(in) :: lat,lon,depth_km
  type(t_gf_location), intent(out) :: loc
  integer, intent(out) :: ierr

  ! local parameters
  double precision, dimension(NDIM) :: xyz_target,xyz
  double precision, dimension(NDIM,NDIM) :: jinv,nu
  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: xyz_elem
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision, dimension(NGLLZ) :: zigll,wzgll
  double precision, dimension(1) :: no_spline
  integer, dimension(GF_NCAND) :: cand
  double precision :: theta,phi,r_surface,elevation
  double precision :: xi,eta,gamma,jacobian
  double precision :: xi_max,xi_max_best,dist,dist_best,anchor_err,r_km
  integer :: nfound,icand,ielem,ix0,iy0,iz0,ielem_best

  loc = t_gf_location()

  if (.not. db%is_open) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_locate_source: database is not open')
    return
  endif

  ! distances are reported in km, on the database's own planet rather than
  ! through shared_parameters, so this stays planet-agnostic
  r_km = db%R_PLANET / 1000.d0

  !--- geographic chain -------------------------------------------------

  ! surface elevation; the grid is only read the first time it is wanted
  elevation = 0.d0
  if (db%topography) then
    if (.not. db%topo_loaded) then
      call gf_load_topo(db,ierr)
      if (ierr /= GF_OK) return
    endif
    call gf_topo_elevation(db,lat,lon,elevation)
  endif

  ! the spline arrays only exist when the database has ELLIPTICITY set, so
  ! the two cases are separate calls rather than one call with a conditional
  ! actual argument -- an unallocated allocatable cannot be passed at all
  if (db%ellipticity) then
    call gf_geographic_to_cartesian(lat,lon,depth_km,db%ellipticity,elevation, &
                                    db%nspl,db%rspl,db%ellipicity_spline,db%ellipicity_spline2, &
                                    db%R_PLANET,xyz_target,theta,phi,r_surface,ierr)
  else
    no_spline(1) = 0.d0
    call gf_geographic_to_cartesian(lat,lon,depth_km,db%ellipticity,elevation, &
                                    0,no_spline,no_spline,no_spline, &
                                    db%R_PLANET,xyz_target,theta,phi,r_surface,ierr)
  endif
  if (ierr /= GF_OK) return

  call gf_source_nu(theta,phi,nu)

  !--- element search ---------------------------------------------------

  call gf_find_candidates(db,xyz_target,cand,nfound,ierr)
  if (ierr /= GF_OK) return

  ! GLL abscissae for the initial guess; GAUSSALPHA/GAUSSBETA are 0, i.e.
  ! plain Gauss-Lobatto-Legendre, matching define_derivation_matrices.f90:60
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  xi_max_best = huge(1.d0)
  dist_best = huge(1.d0)
  ielem_best = 0

  do icand = 1,nfound
    ielem = cand(icand)

    call gf_read_element_coords(db,ielem,xyz_elem,ierr)
    if (ierr /= GF_OK) return

    call gf_gather_anchors(xyz_elem,xelm,yelm,zelm)

    call gf_initial_guess(xyz_elem,xyz_target,ix0,iy0,iz0)

    call gf_find_local_coords(xelm,yelm,zelm,xigll,yigll,zigll,xyz_target, &
                              ix0,iy0,iz0,xi,eta,gamma,xyz,jinv,jacobian,ierr)
    if (ierr /= GF_OK) then
      ! a degenerate element is a property of that element, not of the
      ! request: try the next candidate rather than failing the whole locate
      cycle
    endif

    xi_max = max(abs(xi),abs(eta),abs(gamma))
    dist = sqrt((xyz(1)-xyz_target(1))**2 + (xyz(2)-xyz_target(2))**2 + (xyz(3)-xyz_target(3))**2)

    ! remembers the best near-miss, so a failure can say how close it got
    if (xi_max < xi_max_best) then
      xi_max_best = xi_max
      dist_best = dist
      ielem_best = ielem
    endif

    if (xi_max <= GF_XI_TOL) then
      ! accepted: check the element geometry is the tri-quadratic one this
      ! whole approach assumes, then fill the result
      call gf_check_anchors(db,ielem,anchor_err,ierr)
      if (ierr /= GF_OK) return

      if (anchor_err > GF_ANCHOR_TOL) then
        call gf_set_error(ierr,GF_ERR_MISMATCH, &
          'the 27 anchors do not reproduce the stored GLL coordinates of the located element; ' &
          //'this database was probably written from a USE_GLL = .true. mesh')
        return
      endif

      loc%ielem = ielem
      loc%morton_hex = db%morton_hex(ielem)
      loc%xi = xi
      loc%eta = eta
      loc%gamma = gamma
      loc%xyz(:) = xyz(:)
      loc%xyz_target(:) = xyz_target(:)
      loc%jinv(:,:) = jinv(:,:)
      loc%jacobian = jacobian

      ! the derivative of the inverse Jacobian, once, on the accepted
      ! element (Stage 8); its first-order part is the map just made, so
      ! jinv above and the one inside agree bitwise
      call gf_shape3D_map_2nd(xelm,yelm,zelm,xi,eta,gamma,xyz,jinv,jacobian,loc%djinv,ierr)
      if (ierr /= GF_OK) return
      loc%nu(:,:) = nu(:,:)
      loc%theta = theta
      loc%phi = phi
      loc%r_surface = r_surface
      loc%distance_km = dist * r_km
      loc%anchor_err = anchor_err

      ierr = GF_OK
      return
    endif
  enddo

  !--- nothing accepted -------------------------------------------------
  !
  ! Deliberately an error. find_containing_element (gf_cross_validate.py:189)
  ! warns and returns the nearest centroid with xi=eta=gamma=0, which
  ! produces a seismogram that looks entirely reasonable and is wrong.

  if (ielem_best > 0) then
    call gf_set_error(ierr,GF_ERR_NO_ELEMENT, &
      'no database element contains this position: the best of ' &
      //trim(gf_itoa(nfound))//' candidates reached max|xi,eta,gamma| = ' &
      //trim(gf_ftoa(xi_max_best))//' at a distance of '//trim(gf_ftoa(dist_best*r_km)) &
      //' km, against a tolerance of '//trim(gf_ftoa(GF_XI_TOL)))
  else
    call gf_set_error(ierr,GF_ERR_NO_ELEMENT, &
      'no database element could be evaluated at this position')
  endif

  end subroutine gf_locate_source

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_initial_guess(xyz_elem,xyz_target,ix0,iy0,iz0)

! nearest interior GLL point of an element, as the Newton starting point
!
! locate_point.f90:310-330 does the same scan, over ibool; the stored
! element array makes the indirection unnecessary. Interior points only
! (2..NGLL-1), as there: starting on a face would put the first iterate on
! the clamp boundary.

  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM

  implicit none

  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ), intent(in) :: xyz_elem
  double precision, dimension(NDIM), intent(in) :: xyz_target
  integer, intent(out) :: ix0,iy0,iz0

  ! local parameters
  integer :: i,j,k
  double precision :: dx,dy,dz,d_sq,d_min_sq

  ix0 = (NGLLX+1)/2
  iy0 = (NGLLY+1)/2
  iz0 = (NGLLZ+1)/2
  d_min_sq = huge(1.d0)

  do k = 2,NGLLZ-1
    do j = 2,NGLLY-1
      do i = 2,NGLLX-1
        dx = xyz_target(1) - xyz_elem(1,i,j,k)
        dy = xyz_target(2) - xyz_elem(2,i,j,k)
        dz = xyz_target(3) - xyz_elem(3,i,j,k)
        d_sq = dx*dx + dy*dy + dz*dz
        if (d_sq < d_min_sq) then
          d_min_sq = d_sq
          ix0 = i
          iy0 = j
          iz0 = k
        endif
      enddo
    enddo
  enddo

  end subroutine gf_initial_guess

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_check_anchors(db,ielem,max_err,ierr)

! reconstructs every stored GLL coordinate of an element from its 27 anchors
!
! This is the guard on the plan's central geometric premise. With
! USE_GLL = .false. the mesher applies topography and ellipticity to the 27
! anchors and re-interpolates the GLL points with the tri-quadratic shape
! functions, so the stored xyz *is* a tri-quadratic sampled at GLL points and
! the anchors must reproduce it. If they do not, the database came from a
! USE_GLL = .true. mesh and the whole anchor route -- the coordinate map, the
! Jacobian, the strain -- is invalid.
!
! The floor is float32 storage, not arithmetic: see GF_ANCHOR_TOL in
! gf_par.F90. `max_err` is returned rather than compared here so that
! callers can report the measured value; the shipped examples sit at ~6e-8.

  use constants, only: NGNOD,NGLLX,NGLLY,NGLLZ,NDIM,GAUSSALPHA,GAUSSBETA

  implicit none

  type(t_gfdb), intent(in) :: db
  integer, intent(in) :: ielem
  double precision, intent(out) :: max_err
  integer, intent(out) :: ierr

  ! local parameters
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: xyz_elem
  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision, dimension(NGLLZ) :: zigll,wzgll
  double precision, dimension(NDIM) :: xyz
  double precision, dimension(NDIM,NDIM) :: jinv
  double precision :: jacobian,err
  integer :: i,j,k,idim

  max_err = 0.d0

  call gf_read_element_coords(db,ielem,xyz_elem,ierr)
  if (ierr /= GF_OK) return

  call gf_gather_anchors(xyz_elem,xelm,yelm,zelm)

  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        call gf_shape3D_map(xelm,yelm,zelm,xigll(i),yigll(j),zigll(k),xyz,jinv,jacobian,ierr)
        if (ierr /= GF_OK) return

        do idim = 1,NDIM
          err = abs(xyz(idim) - xyz_elem(idim,i,j,k))
          if (err > max_err) max_err = err
        enddo
      enddo
    enddo
  enddo

  ierr = GF_OK

  end subroutine gf_check_anchors

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_check_anchors_all(db,worst_err,ielem_worst,ierr)

! gf_check_anchors over the whole database
!
! Not run at open time: it opens every element's coordinates.h5, which is
! 70,000 file opens at production scale. gf_locate_source checks the element
! it accepted, once; this sweep is the deliberate, explicit form behind
! `xgf3d --check-anchors` and tests/gf3d/test_gf_anchors.f90.

  implicit none

  type(t_gfdb), intent(in) :: db
  double precision, intent(out) :: worst_err
  integer, intent(out) :: ielem_worst
  integer, intent(out) :: ierr

  ! local parameters
  double precision :: err
  integer :: ielem

  worst_err = 0.d0
  ielem_worst = 0

  if (.not. db%is_open) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_check_anchors_all: database is not open')
    return
  endif

  do ielem = 1,db%nelem
    call gf_check_anchors(db,ielem,err,ierr)
    if (ierr /= GF_OK) return
    if (err > worst_err) then
      worst_err = err
      ielem_worst = ielem
    endif
  enddo

  ierr = GF_OK

  end subroutine gf_check_anchors_all

!
!-------------------------------------------------------------------------------------------------
!

  function gf_itoa(i) result(str)

! integer to string, for error messages

  implicit none

  integer, intent(in) :: i
  character(len=16) :: str

  write(str,'(i0)') i

  end function gf_itoa

!
!-------------------------------------------------------------------------------------------------
!

  function gf_ftoa(x) result(str)

! double to a short string, for error messages

  implicit none

  double precision, intent(in) :: x
  character(len=16) :: str

  write(str,'(f0.6)') x

  end function gf_ftoa

  end module gf_locate
