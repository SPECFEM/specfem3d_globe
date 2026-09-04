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
!---- Source geometry: geographic coordinates to the mesh, and back out of
!---- a single element.
!----
!---- Ported from src/specfem3D/locate_sources.f90:220-354 (the geographic
!---- chain and the `nu` orientation matrix) and
!---- src/specfem3D/locate_point.f90:435-529 (the Newton iteration), with
!---- the solver's module arrays turned into arguments.
!----
!---- **specfem3d_globe is the reference implementation, not GF3DF.** GF3DF
!---- carries its own copy of this chain, and that copy is stale in two ways
!---- that would put a systematic, depth-dependent error into every source
!---- position: it implements the deprecated USE_OLD_VERSION_FORMAT branch
!---- (ellipticity applied *after* subtracting depth, locate_sources.f90:329
!---- -- and USE_OLD_VERSION_FORMAT is .false. in setup/constants.h), and it
!---- calls lat_2_geocentric_colat_dble with two arguments where it now takes
!---- three. Neither is ported here.
!----
!---- No `use hdf5`, no `use specfem_par`: this is a kernel module. The
!---- database handle is passed in as plain arrays for the same reason --
!---- everything HDF5-shaped (loading the topography grid, reading an
!---- element) lives in gf_locate.F90 one level up.
!----

  module gf_geometry

  use gf_par, only: gf_set_error,GF_OK,GF_ERR_ARG,GF_ERR_MISMATCH

  use gf_shape3D, only: gf_shape3D_map

  implicit none

  private

  public :: gf_geographic_to_cartesian
  public :: gf_source_nu
  public :: gf_find_local_coords
  public :: gf_gather_anchors

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_geographic_to_cartesian(lat,lon,depth_km,is_elliptical,elevation, &
                                        nspl,rspl,ell1,ell2,R_PLANET, &
                                        xyz_target,theta,phi,r_surface,ierr)

! converts geographic (lat,lon,depth) to a non-dimensional Cartesian position
!
! Follows locate_sources.f90:220-354 step for step: geocentric colatitude,
! topography on the surface radius, ellipticity on the *surface* radius, and
! only then the depth subtraction.
!
! `is_elliptical` is the database's own ELLIPTICITY attribute, standing in
! for the solver's ELLIPTICITY_VAL (which comes from values_from_mesher.h).
! It is passed to lat_2_geocentric_colat_dble as well as gating the
! ellipticity correction, exactly as in the solver.
!
! `elevation` is passed in already evaluated, in metres, rather than the
! topography grid, so that this module needs no get_topo_bathy() and
! therefore no model_topo_bathy.shared.o and no MPI stubs behind it. The
! caller supplies gf_topo_elevation(db,...) or zero. The solver's own
! TOPOGRAPHY switch shows up here as elevation == 0.

  use constants, only: DEGREES_TO_RADIANS,R_UNIT_SPHERE,NDIM,NR_DENSITY

  implicit none

  double precision, intent(in) :: lat,lon,depth_km
  logical, intent(in) :: is_elliptical
  double precision, intent(in) :: elevation
  integer, intent(in) :: nspl
  double precision, dimension(:), intent(in) :: rspl,ell1,ell2
  double precision, intent(in) :: R_PLANET
  double precision, dimension(NDIM), intent(out) :: xyz_target
  double precision, intent(out) :: theta,phi,r_surface
  integer, intent(out) :: ierr

  ! local parameters
  double precision :: xlat,xlon,depth,r0,r_target
  double precision :: sint,cost,sinp,cosp

  ! add_ellipticity_rtheta() declares its spline dummies with the explicit
  ! shape rspl(NR_DENSITY) (make_ellipticity.f90:495), while the database
  ! allocates them to nspl -- 628 against NR_DENSITY = 640 for the shipped
  ! examples. Passing the short arrays straight through is not standard
  ! conforming and trips -fcheck=bounds, so they are copied into full-length
  ! buffers here. 640 doubles per call is nothing next to a locate.
  double precision, dimension(NR_DENSITY) :: rspl_p,ell1_p,ell2_p

  xyz_target(:) = 0.d0
  theta = 0.d0
  phi = 0.d0
  r_surface = 0.d0

  if (R_PLANET <= 0.d0) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_geographic_to_cartesian: R_PLANET must be positive')
    return
  endif

  xlat = lat
  xlon = lon

  ! limits longitude to [0.0,360.0]
  if (xlon < 0.d0 ) xlon = xlon + 360.d0
  if (xlon > 360.d0 ) xlon = xlon - 360.d0

  ! convert geographic latitude lat (degrees) to geocentric colatitude theta (radians)
  call lat_2_geocentric_colat_dble(xlat,theta,is_elliptical)

  ! longitude
  phi = xlon*DEGREES_TO_RADIANS

  ! theta to [0,PI] and phi to [0,2PI]
  call reduce(theta,phi)

  sint = sin(theta)
  cost = cos(theta)
  sinp = sin(phi)
  cosp = cos(phi)

  ! point depth (in m)
  depth = depth_km*1000.0d0

  ! normalized source radius
  r0 = R_UNIT_SPHERE

  ! elevation of position (zero when the database has no topography)
  r0 = r0 + elevation/R_PLANET

  ! ellipticity, applied to the surface radius before the depth is removed
  if (is_elliptical) then
    if (nspl <= 0 .or. nspl > NR_DENSITY) then
      call gf_set_error(ierr,GF_ERR_MISMATCH, &
        'ellipticity requested but the spline table has an unusable length')
      return
    endif
    if (size(rspl) < nspl .or. size(ell1) < nspl .or. size(ell2) < nspl) then
      call gf_set_error(ierr,GF_ERR_ARG,'ellipticity spline arrays are shorter than nspl')
      return
    endif

    rspl_p(:) = 0.d0
    ell1_p(:) = 0.d0
    ell2_p(:) = 0.d0
    rspl_p(1:nspl) = rspl(1:nspl)
    ell1_p(1:nspl) = ell1(1:nspl)
    ell2_p(1:nspl) = ell2(1:nspl)

    call add_ellipticity_rtheta(r0,theta,nspl,rspl_p,ell1_p,ell2_p)
  endif

  ! stores surface radius above the point
  r_surface = r0

  ! subtracts depth (given in m)
  r0 = r0 - depth/R_PLANET
  r_target = r0

  ! compute the Cartesian position
  xyz_target(1) = r_target*sint*cosp
  xyz_target(2) = r_target*sint*sinp
  xyz_target(3) = r_target*cost

  ierr = GF_OK

  end subroutine gf_geographic_to_cartesian

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_source_nu(theta,phi,nu)

! orientation matrix at a point: rows are North, East and vertical (up),
! each expressed in global Cartesian x/y/z
!
! Transcribed from locate_sources.f90:275-312, including the Harvard
! normal-mode convention for `n`. The solver builds the same matrix from a
! (stazi,stdip) pair per component; those pairs are compile-time constants
! there -- (0,0), (90,0), (0,-90) -- so they are written out here.
!
! This is the matrix compute_arrays_source_forcesolution() contracts a
! FORCESOLUTION direction vector with (compute_arrays_source.f90:172-177),
! and it is why row order N,E,Z matters: the reciprocal runs that filled the
! database applied their forces along exactly these three directions at the
! station, so the same ordering indexes the stored force component.

  use constants, only: DEGREES_TO_RADIANS,NDIM

  implicit none

  double precision, intent(in) :: theta,phi
  double precision, dimension(NDIM,NDIM), intent(out) :: nu

  ! local parameters
  double precision :: sint,cost,sinp,cosp
  double precision :: stazi,stdip,thetan,phin
  double precision, dimension(NDIM) :: n
  integer :: iorientation

  sint = sin(theta)
  cost = cos(theta)
  sinp = sin(phi)
  cosp = cos(phi)

  do iorientation = 1,3
    ! initializes azimuth/dip
    stazi = 0.d0
    stdip = 0.d0
    ! North
    if (iorientation == 1) then
      stazi = 0.d0
      stdip = 0.d0
    ! East
    else if (iorientation == 2) then
      stazi = 90.d0
      stdip = 0.d0
    ! Vertical
    else
      stazi = 0.d0
      stdip = - 90.d0
    endif

    ! get the orientation of the position
    thetan = (90.0d0+stdip)*DEGREES_TO_RADIANS
    phin = stazi*DEGREES_TO_RADIANS

    ! we use the same convention as in Harvard normal modes for the orientation
    ! vertical component
    n(1) = cos(thetan)
    ! N-S component
    n(2) = - sin(thetan)*cos(phin)
    ! E-W component
    n(3) = sin(thetan)*sin(phin)

    ! get the Cartesian components of n in the model: nu
    nu(iorientation,1) = n(1)*sint*cosp + n(2)*cost*cosp - n(3)*sinp
    nu(iorientation,2) = n(1)*sint*sinp + n(2)*cost*sinp + n(3)*cosp
    nu(iorientation,3) = n(1)*cost - n(2)*sint
  enddo

  end subroutine gf_source_nu

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_gather_anchors(xyz_elem,xelm,yelm,zelm)

! extracts the 27 anchor coordinates from a stored element
!
! The solver reaches its anchors through ibool
! (locate_point.f90:436-441); here the element is already a dense
! xyz(3,NGLLX,NGLLY,NGLLZ) array, so the same anchor_ia* indices apply
! directly. hex_nodes_anchor_ijk() is called rather than a hand-written
! table so that a future reordering of the 27 nodes changes both this and
! gf_shape3D_functions together.

  use constants, only: NGNOD,NGLLX,NGLLY,NGLLZ

  implicit none

  double precision, dimension(3,NGLLX,NGLLY,NGLLZ), intent(in) :: xyz_elem
  double precision, dimension(NGNOD), intent(out) :: xelm,yelm,zelm

  ! local parameters
  integer, dimension(NGNOD) :: anchor_iax,anchor_iay,anchor_iaz
  integer :: ia

  call hex_nodes_anchor_ijk(anchor_iax,anchor_iay,anchor_iaz)

  do ia = 1,NGNOD
    xelm(ia) = xyz_elem(1,anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia))
    yelm(ia) = xyz_elem(2,anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia))
    zelm(ia) = xyz_elem(3,anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia))
  enddo

  end subroutine gf_gather_anchors

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_find_local_coords(xelm,yelm,zelm,xigll,yigll,zigll, &
                                  xyz_target,ix0,iy0,iz0, &
                                  xi,eta,gamma,xyz,jinv,jacobian,ierr)

! Newton iteration for the local coordinates of a point inside one element
!
! This is find_local_coordinates() (locate_point.f90:395-531) with the
! solver's module arrays -- ibool, xstore/ystore/zstore, anchor_ia*, xigll --
! turned into arguments, and recompute_jacobian replaced by gf_shape3D_map.
! Everything numerical is preserved deliberately, because matching the
! solver's element choice is the whole point of the containment rule in
! gf_par (GF_XI_TOL):
!
!   * NUM_ITER iterations, no residual-based convergence test;
!   * an iterate that is *worse* than the previous one ends the loop, leaving
!     xi/eta/gamma at the values that produced it;
!   * a step longer than 1 in (dxi,deta,dgamma) is damped by 0.33333333333;
!   * each coordinate is clamped to +-1.10 after every update;
!   * a final map so the returned xyz is consistent with the returned
!     xi/eta/gamma rather than with the last iterate.
!
! The solver's POINT_CAN_BE_BURIED branch is not carried over: it forces
! gamma = 1 for receivers pinned to the free surface, and every point this
! library locates is a source at depth.

  use constants, only: NGNOD,NGLLX,NGLLY,NGLLZ,NDIM,NUM_ITER,HUGEVAL

  implicit none

  double precision, dimension(NGNOD), intent(in) :: xelm,yelm,zelm
  double precision, dimension(NGLLX), intent(in) :: xigll
  double precision, dimension(NGLLY), intent(in) :: yigll
  double precision, dimension(NGLLZ), intent(in) :: zigll
  double precision, dimension(NDIM), intent(in) :: xyz_target
  integer, intent(in) :: ix0,iy0,iz0
  double precision, intent(out) :: xi,eta,gamma
  double precision, dimension(NDIM), intent(out) :: xyz
  double precision, dimension(NDIM,NDIM), intent(out) :: jinv
  double precision, intent(out) :: jacobian
  integer, intent(out) :: ierr

  ! local parameters
  integer :: iter_loop
  double precision :: dx,dy,dz
  double precision :: d_sq,d_min_sq
  double precision :: dxi,deta,dgamma

  ! use initial guess in xi and eta
  xi = xigll(ix0)
  eta = yigll(iy0)
  gamma = zigll(iz0)

  d_min_sq = HUGEVAL

  ! iterate to solve the non linear system
  do iter_loop = 1,NUM_ITER

    ! recompute Jacobian for the new point
    call gf_shape3D_map(xelm,yelm,zelm,xi,eta,gamma,xyz,jinv,jacobian,ierr)
    if (ierr /= GF_OK) return

    ! compute distance to target location
    dx = - (xyz(1) - xyz_target(1))
    dy = - (xyz(2) - xyz_target(2))
    dz = - (xyz(3) - xyz_target(3))

    ! distance squared
    d_sq = dx*dx + dy*dy + dz*dz

    ! compute increments
    if (d_sq < d_min_sq) then
      d_min_sq = d_sq

      dxi = jinv(1,1)*dx + jinv(1,2)*dy + jinv(1,3)*dz
      deta = jinv(2,1)*dx + jinv(2,2)*dy + jinv(2,3)*dz
      dgamma = jinv(3,1)*dx + jinv(3,2)*dy + jinv(3,3)*dz
    else
      ! new position is worse than old one, no change necessary
      ! stop, no further improvements
      exit
    endif

    ! decreases step length if step is large
    if ((dxi*dxi + deta*deta + dgamma*dgamma) > 1.0d0) then
      dxi = dxi * 0.33333333333d0
      deta = deta * 0.33333333333d0
      dgamma = dgamma * 0.33333333333d0
    endif

    ! update values
    xi = xi + dxi
    eta = eta + deta
    gamma = gamma + dgamma

    ! impose that we stay in that element
    ! we can go slightly outside the [1,1] segment since with finite elements
    ! the polynomial solution is defined everywhere
    if (abs(xi) > 1.10d0) xi = sign(1.10d0,xi)
    if (abs(eta) > 1.10d0) eta = sign(1.10d0,eta)
    if (abs(gamma) > 1.10d0) gamma = sign(1.10d0,gamma)

  ! end of non linear iterations
  enddo

  ! compute final coordinates of point found
  call gf_shape3D_map(xelm,yelm,zelm,xi,eta,gamma,xyz,jinv,jacobian,ierr)

  end subroutine gf_find_local_coords

  end module gf_geometry
