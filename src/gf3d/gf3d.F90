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
!---- The public face of libgf3d.
!----
!---- One module, so that a program outside this tree needs
!----
!----     use gf3d
!----
!---- and nothing else. Everything below is re-exported from the gf_*
!---- modules that implement it; the only code here is gf_release() and
!---- get_seismograms(), both of which exist to spare a caller a sequence it
!---- would otherwise have to get right itself.
!----
!---- Build and link
!---- --------------
!----     make gf3d
!---- installs include/gf3d.mod beside include/gf3d.h, so
!----
!----     <FC> myprog.f90 -I<repo>/include -L<repo>/lib -lgf3d <HDF5 link line>
!----
!---- links against the static library. lib/libgf3d.so is the same library as
!---- a shared object, for ctypes and for anything that dlopen()s it.
!----
!---- Relation to GF3DF
!---- -----------------
!---- The name `gf3d` is GF3DF's, deliberately: a downstream GCMT-style
!---- caller writes `use gf3d` there too. The type aliases t_GF and t_source
!---- and the generic get_seismograms() below carry that as far as it goes
!---- honestly. It does not go all the way, and the differences are not
!---- cosmetic:
!----
!----   * t_source there carries Mrr..Mtp in dyne-cm and a `time_shift`;
!----     t_gf_source here carries a non-dimensional moment_tensor(6) and
!----     min_tshift_src_original, because get_cmt() has already divided by
!----     scaleM. Multiply by src%scale_moment to get dyne-cm back.
!----   * read_GF(file) + GF%get_kdtree() becomes gf_open(dir,db,ierr): the
!----     database is a directory tree, not one file, and the kd-tree is
!----     built on the first locate.
!----   * every routine returns ierr instead of stopping.
!----
!---- So the port of a downstream caller is the open call, the source fields
!---- and an error argument -- not a drop-in.
!----
!---- Errors
!---- ------
!---- No routine reachable from here stops. Each returns one of the GF_*
!---- codes; gf_error_string(code) is a short name and gf_errmsg the message
!---- from the last failure.
!----

  module gf3d

  !--- types, constants and error handling
  use gf_par, only: &
    t_gfdb, t_gf_station, t_gf_source, t_gf_location, t_gf_taxis, t_gf_stf, &
    GF3D_VERSION, GF_NCOMP, GF_MORTON_HEXLEN, GF_STATION_ID_LEN, &
    GF_XI_TOL, GF_NCAND, GF_ANCHOR_TOL, &
    GF_OK, GF_ERR_NO_HDF5, GF_ERR_NO_PATH, GF_ERR_NO_FILE, GF_ERR_HDF5, &
    GF_ERR_IO, GF_ERR_FORMAT, GF_ERR_MISMATCH, GF_ERR_INCOMPLETE, &
    GF_ERR_ALLOC, GF_ERR_ARG, GF_ERR_NO_ELEMENT, GF_ERR_GEOMETRY, &
    GF_SRC_FORCE, GF_SRC_CMT, &
    GF_STF_NONE, GF_STF_GAUSS, GF_STF_HEAVI, GF_STF_TRUNC, &
    gf_set_error, gf_error_string, gf_errmsg, gf_is_finite

  !--- GF3DF-compatible aliases, so that `use gf3d, only: t_GF, t_source`
  !--- keeps naming something
  use gf_par, only: t_GF => t_gfdb, t_source => t_gf_source

  !--- opening and interrogating a database
  use gf_database, only: &
    gf_open, gf_close, gf_load_topo, gf_check_completion, gf_print_info, &
    gf_topo_elevation, gf_topo_gradient

  !--- sources
  use gf_source, only: &
    gf_read_source, gf_read_cmt_source, gf_read_force_source, &
    gf_source_set_cmt, gf_source_set_force, &
    gf_detect_source_type, gf_force_direction, gf_print_source

  !--- locating a source in the mesh
  use gf_locate, only: &
    gf_locate_source, gf_locate_release, gf_locate_tree_owner, &
    gf_check_anchors, gf_check_anchors_all, gf_print_location

  !--- the source time function and the output time axis
  use gf_stf, only: &
    gf_stf_plan, gf_taxis_plan, gf_taxis_times, gf_stf_kind_name, gf_print_stf, &
    gf_hdur_gaussian, gf_default_t0

  !--- extraction
  use gf_seismograms, only: &
    gf_time_axis, gf_seis_plan, gf_seis, &
    gf_write_seis, gf_write_partials, gf_write_dump

  !--- partial derivatives: the count, the slot order, the names and units
  use gf_partials, only: &
    gf_partials_ndp, GF_NDP_MT, GF_NDP_LOC, &
    GF_DP_MRR, GF_DP_MTT, GF_DP_MPP, GF_DP_MRT, GF_DP_MRP, GF_DP_MTP, &
    GF_DP_LAT, GF_DP_LON, GF_DP_DEP, GF_DP_TIM, &
    GF_DP_NAME, GF_DP_UNIT

  !--- SAC output. Only the whole-database writer: the header record and the
  !--- single-trace writer are for xgf3d and for the tests, and exporting
  !--- them would mean exporting t_gf_sac_header with them.
  use gf_sac, only: gf_write_sac, GF_SAC_CHANNEL

  implicit none

  public

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_release(db)

! closes a database and frees the process-wide search tree
!
! gf_close() alone is not enough. The kd-tree lives in module variables of
! src/shared/search_kdtree.f90, which gf_locate owns on the database's
! behalf; gf_close() cannot release it, because gf_locate depends on
! gf_database and the reverse call would be circular. So a caller that opens
! and closes databases over a long-lived process -- an inversion loop, a
! Python interpreter -- must call both, in this order, and this routine is
! that pair.

  implicit none

  type(t_gfdb), intent(inout) :: db

  call gf_locate_release()
  call gf_close(db)

  end subroutine gf_release

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_seismograms(db,src,t0,synt,ierr,t)

! seismograms at every station for one source
!
! The whole of xgf3d's --seis path in one call: locate, plan, extract.
!
!   synt(nsta, 3, nt)   metres, components N/E/Z, allocated here
!   t(nt)               seconds relative to the centroid time; optional
!
! `t0` is where the output axis starts, in seconds before the centroid time,
! and must be non-negative. gf_default_t0() returns what specfem's own
! forward run would use (1.5*hdur for a CMT); pass that unless there is a
! reason to choose another.
!
! For the partial derivatives as well, call get_partials. `src` is
! intent(in): nothing is written back into it.

  implicit none

  type(t_gfdb), intent(inout) :: db          ! inout: the topography grid loads lazily
  type(t_gf_source), intent(in) :: src
  double precision, intent(in) :: t0
  double precision, dimension(:,:,:), allocatable, intent(out) :: synt
  integer, intent(out) :: ierr
  double precision, dimension(:), allocatable, intent(out), optional :: t

  ! local parameters
  double precision, dimension(:,:,:,:), allocatable :: dp

  call extract(db,src,t0,0,synt,dp,ierr,t)

  if (allocated(dp)) deallocate(dp)

  end subroutine get_seismograms

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_partials(db,src,t0,itypsokern,synt,dp,ierr,t)

! seismograms and their partial derivatives for one source
!
!   synt(nsta, 3, nt)      metres, components N/E/Z
!   dp(ndp, nsta, 3, nt)   ndp = 6 for itypsokern = 1, 10 for 2
!   t(nt)                  seconds relative to the centroid time; optional
!
! all allocated here. itypsokern = 1 gives the six moment-tensor partials;
! 2 gives those plus latitude, longitude, depth and centroid time. The
! parameter names and units are GF_DP_NAME and GF_DP_UNIT.
!
! The seismogram comes back too because it is computed on the way -- the
! moment-tensor partials contracted with the moment tensor *are* the
! seismogram -- so there is nothing to be saved by asking for one without
! the other.
!
! Partials are defined for a moment-tensor source only. itypsokern = 3
! (GF3DF's half-duration partial) is not supported; see gf_partials.

  implicit none

  type(t_gfdb), intent(inout) :: db
  type(t_gf_source), intent(in) :: src
  double precision, intent(in) :: t0
  integer, intent(in) :: itypsokern
  double precision, dimension(:,:,:), allocatable, intent(out) :: synt
  double precision, dimension(:,:,:,:), allocatable, intent(out) :: dp
  integer, intent(out) :: ierr
  double precision, dimension(:), allocatable, intent(out), optional :: t

  if (itypsokern < 1) then
    call gf_set_error(ierr,GF_ERR_ARG, &
      'get_partials: itypsokern must be 1 or 2; use get_seismograms for none')
    return
  endif

  call extract(db,src,t0,itypsokern,synt,dp,ierr,t)

  end subroutine get_partials

!
!-------------------------------------------------------------------------------------------------
!

  subroutine extract(db,src,t0,itypsokern,synt,dp,ierr,t)

! locate, plan, allocate, extract -- the body both entry points share
!
! `dp` is always allocated, with a first extent of zero when no partials
! were asked for, because gf_seis takes it by explicit shape. The caller
! that does not want it throws it away.

  implicit none

  type(t_gfdb), intent(inout) :: db
  type(t_gf_source), intent(in) :: src
  double precision, intent(in) :: t0
  integer, intent(in) :: itypsokern
  double precision, dimension(:,:,:), allocatable, intent(out) :: synt
  double precision, dimension(:,:,:,:), allocatable, intent(out) :: dp
  integer, intent(out) :: ierr
  double precision, dimension(:), allocatable, intent(out), optional :: t

  ! local parameters
  type(t_gf_location) :: loc
  type(t_gf_taxis) :: tax
  type(t_gf_stf) :: stf
  double precision, dimension(:), allocatable :: onset,tsec
  integer :: ndp,ier

  if (.not. db%is_open) then
    call gf_set_error(ierr,GF_ERR_ARG,'database is not open')
    return
  endif

  call gf_partials_ndp(itypsokern,ndp,ierr)
  if (ierr /= GF_OK) return

  call gf_locate_source(db,src%latitude,src%longitude,src%depth,loc,ierr)
  if (ierr /= GF_OK) return

  call gf_seis_plan(db,src,t0,tax,stf,ierr)
  if (ierr /= GF_OK) return

  allocate(synt(db%nstations,GF_NCOMP,tax%nt), &
           dp(ndp,db%nstations,GF_NCOMP,tax%nt), &
           tsec(tax%nt),onset(db%nstations),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate the output arrays')
    return
  endif

  ! gf_seis refuses a source type that has no partials, so the check the
  ! facade used to make itself is not repeated here
  call gf_seis(db,src,loc,tax,stf,itypsokern,ndp,synt,dp,tsec,onset,ierr)

  if (present(t)) then
    allocate(t(tax%nt),stat=ier)
    if (ier /= 0) then
      call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate the time axis')
    else
      t(:) = tsec(:)
    endif
  endif

  deallocate(onset,tsec)

  end subroutine extract

  end module gf3d
