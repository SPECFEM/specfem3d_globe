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
!---- Seismograms from a force or a moment-tensor source, and the ASCII output.
!----
!---- Reciprocity, written out once because every index below follows from it
!---- --------------------------------------------------------------------
!---- The database holds, per element and station,
!----
!----   displacement(a, d, i,j,k, t)
!----     = the displacement at the GLL point, in Cartesian direction d,
!----       produced by a force at the *station* along its local direction a
!----       (a: N=1, E=2, Z=3).
!----
!---- Reciprocity exchanges source and receiver:
!----
!----   u_a(x_r ; force at x_s along d) = u_d(x_s ; force at x_r along a)
!----
!---- so the seismogram at the station in component a, from a source force
!---- along the Cartesian unit vector fhat, is
!----
!----   seis(a,t) = SUM_d fhat(d) * displacement(a, d, x_s, t)
!----
!---- The first index is therefore the *output component* and the second is
!---- contracted with the source direction. Getting that pair the wrong way
!---- round is not a crash; it is a plausible seismogram with the wrong
!---- radiation pattern.
!----
!---- Amplitude -- the factor gf_cross_validate.py omits
!---- ------------------------------------------------
!---- The stored field was produced by a force of magnitude
!---- `factor_force_source` at the station, so per unit force it is that
!---- divided out, and the source's own magnitude multiplied back in:
!----
!----   scale = src%factor_force_source / db%stations(ista)%factor_force_source
!----
!---- gf_cross_validate.py:reconstruct_force applies no such factor. That is
!---- correct only by accident on the shipped examples: the Snakefile's
!---- reciprocal runs use FACTOR_FORCE = 1.0d15 and validation_data's
!---- FORCESOLUTION also says 1.0d15, so the ratio is one. Both values are
!---- non-dimensional (get_force divides by scaleF, and the writer stores the
!---- attribute after the same division), so the ratio is unit-consistent.
!----
!---- For a moment-tensor source the same reciprocal normalisation applies,
!---- but there is no source *force* to multiply back in -- the moment tensor
!---- carries the magnitude itself:
!----
!----   scale = 1 / db%stations(ista)%factor_force_source
!----
!---- with M already non-dimensionalised by scaleM inside get_cmt. This is
!---- the factor gf_cross_validate.py writes as
!---- mt_scale = 1/(SCALE_M * factor_force_nondim): our M supplies the first
!---- half, this supplies the second. Unlike the force case the ratio is
!---- *not* one on the shipped examples -- it is 1.05e10 -- so omitting it is
!---- not a subtle error but a trace ten orders of magnitude too small.
!----
!---- The source time function, and the order of operations
!---- ------------------------------------------------------
!---- The stored field is the response to the reciprocal run's Gaussian.
!---- The seismogram a caller wants is the response to the source time
!---- function specfem would use for the source file in hand, and getting
!---- from one to the other is a single convolution, planned by
!---- gf_seis_plan and applied by gf_stf (see the header of gf_stf.F90 for
!---- the derivation and for why there is no time integration).
!----
!---- That convolution acts on the *contracted* trace, after the spatial
!---- derivative for a moment tensor and after the direction contraction for
!---- a force, never on the raw displacement. For a moment tensor that order
!---- is required, not preferred: the reciprocal force pulse leaves the
!---- planet with permanent momentum, so the raw displacement carries a
!---- rigid translation growing linearly in time, and it is the strain that
!---- removes it.
!----
!---- Time axis
!---- ---------
!---- The stored axis is
!----
!----   t(i) = (i*subsample_step - 1)*dt - t0,   i = 1..nt_subsampled
!----
!---- verified against the writer rather than assumed: green_function_io.F90
!---- :394 snapshots when mod(it,GF_SUBSAMPLE_STEP) == 0, so snapshot i is
!---- solver step i*subsample_step, whose simulation time is (it-1)*dt - t0.
!----
!---- The output axis is that axis extended to the left by whole samples,
!---- with zeros, so that a requested start time -- by default 1.5*hdur, the
!---- forward run's own -- is covered (gf_taxis_plan). Nothing is resampled
!---- and nothing removed: every stored sample and its time survive bitwise.
!----
!---- Use dt_sub = dt*subsample_step for anything that steps along this axis:
!---- 0.4 s, not 0.1 s, in the shipped global example. GF3DF's GF%dt was
!---- already the stored spacing, so a literal port silently uses the solver
!---- step and is wrong by a factor of subsample_step.
!----

  module gf_seismograms

  ! at module scope because the derived types below are sized by them; the
  ! procedures re-import the same names, which is legal and keeps each one's
  ! dependencies readable
  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM

  use gf_par, only: t_gfdb,t_gf_location,t_gf_source,t_gf_stf,t_gf_taxis,gf_set_error, &
                    gf_is_finite, &
                    GF_OK,GF_ERR_ARG,GF_ERR_ALLOC,GF_ERR_IO,GF_ERR_MISMATCH, &
                    GF_NCOMP,GF3D_VERSION,GF_STF_TRUNC,GF_STF_HEAVI, &
                    GF_SRC_FORCE,GF_SRC_CMT

  use gf_database, only: gf_dir_exists,gf_topo_elevation,gf_topo_gradient

  use gf_shared_params, only: gf_init_shared_params

  use gf_element_io, only: gf_read_element_displ

  use gf_interp, only: gf_interp_weights,gf_interp_weights_deriv,gf_interp_weights_deriv2, &
                       gf_interp_trace

  use gf_strain, only: gf_strain_dweights,gf_strain_ddweights,gf_strain_trace,GF_VOIGT

  use gf_moment, only: gf_rotate_moment_tensor,gf_rotate_moment_tensor_deriv,gf_moment_contract

  use gf_stf, only: gf_stf_plan,gf_taxis_plan,gf_taxis_times,gf_stf_onset,gf_stf_kind_name, &
                    t_gf_stf_work,gf_stf_work_init,gf_stf_convert,gf_stf_work_free

  use gf_source, only: gf_force_direction

  use gf_partials, only: gf_partials_mt,gf_partials_time,gf_partials_loc,gf_partials_ndp, &
                         GF_NDP_MT,GF_NDP_LOC,GF_DP_NAME,GF_DP_UNIT,GF_DP_LAT,GF_DP_TIM

  use gf_geo_chain, only: gf_geographic_jacobian

  implicit none

  private

  public :: gf_time_axis
  public :: gf_seis_plan
  public :: gf_seis
  public :: gf_write_seis
  public :: gf_write_partials
  public :: gf_write_dump

  ! component labels, in the stored force order (green_function_io.F90:38)
  character(len=1), dimension(GF_NCOMP), parameter :: GF_COMP_NAME = (/ 'N','E','Z' /)

  !-----------------------------------------------------------------
  ! everything about a source's position that does not depend on the
  ! station, built once by gf_seis_geometry
  !
  ! Which parts are filled depends on what is being extracted: the force
  ! direction for a force source, the strain weights and the rotated moment
  ! tensor for a moment tensor, and the second-derivative tables and the
  ! geographic Jacobian only when the centroid partials are wanted.
  ! Fixed-size rather than allocatable -- ddw is the big one at ~9 kB, which
  ! gf_seis_cmt_partials already carried as a local.
  !-----------------------------------------------------------------

  type :: t_gf_seis_geom
    !--- every kind ---
    ! from gf_interp_weights_deriv, i.e. lagrange_any. See gf_seis_geometry
    ! for why these must NOT come from gf_interp_weights_deriv2.
    double precision, dimension(NGLLX) :: hxi  = 0.d0
    double precision, dimension(NGLLY) :: heta = 0.d0
    double precision, dimension(NGLLZ) :: hgam = 0.d0

    !--- GF_SRC_FORCE ---
    double precision, dimension(NDIM) :: fhat = 0.d0

    !--- GF_SRC_CMT ---
    double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM) :: dw = 0.d0
    double precision, dimension(NDIM,NDIM) :: m_cart = 0.d0

    !--- itypsokern = 2 only ---
    logical :: want_loc = .false.
    double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM,NDIM) :: ddw = 0.d0
    double precision, dimension(NDIM,NDIM) :: dm_dtheta = 0.d0, dm_dphi = 0.d0
    double precision, dimension(NDIM,3) :: dxds = 0.d0
    double precision :: dtheta_dlat = 0.d0, dphi_dlon = 0.d0
  end type t_gf_seis_geom

  !-----------------------------------------------------------------
  ! the per-extraction buffers, allocated once and reused for every station
  !
  ! `displ` is the large one -- 21 MB in the shipped global example -- which
  ! is why the element is read once and the stations looped inside it. It is
  ! a component rather than a dummy argument so that there is a single place
  ! for a later element cache to fill instead of reading.
  !
  ! Which of the rest exist depends on the extraction: `g` for a force
  ! source, `eps` for a moment tensor, and `deps`/`dpl`/`dp10` only for the
  ! centroid partials.
  !-----------------------------------------------------------------

  type :: t_gf_seis_work
    real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: displ
    double precision, dimension(:,:,:), allocatable :: g       ! (3,3,nt_db)      force
    double precision, dimension(:,:,:), allocatable :: eps     ! (6,3,nt_db)      cmt
    double precision, dimension(:,:,:,:), allocatable :: deps  ! (6,3,nt_db,NDIM) itypsokern 2
    double precision, dimension(:,:,:), allocatable :: dpm     ! (6,3,nt)         itypsokern >= 1
    double precision, dimension(:,:,:), allocatable :: dpl     ! (3,3,nt)         itypsokern 2
    double precision, dimension(:), allocatable :: dp10        ! (nt)             itypsokern 2
    type(t_gf_stf_work) :: stfw
  end type t_gf_seis_work

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_time_axis(db,nt,t)

! sample times of the stored traces, in seconds relative to the origin
!
! The output axis with no extension: one expression in the tree for the
! axis, in gf_taxis_times, so the two cannot drift apart.

  implicit none

  type(t_gfdb), intent(in) :: db
  integer, intent(in) :: nt
  double precision, dimension(nt), intent(out) :: t

  ! local parameters
  type(t_gf_taxis) :: tax

  tax%npad = 0
  tax%subsample_step = db%subsample_step
  tax%dt = db%dt
  tax%t0_db = db%t0

  call gf_taxis_times(tax,nt,t)

  end subroutine gf_time_axis

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_seis_plan(db,src,t0_req,tax,stf,ierr)

! decides the source time function conversion and the output axis for a
! source, before any element is read
!
! `t0_req` is the start time in seconds before the origin, and must be a
! real one: a caller that wants specfem's own rule calls gf_default_t0
! first and passes the result. Only the two boundaries with a human user --
! `xgf3d --t0` left off, and a negative t0_req through the C ABI -- still
! carry a sentinel, and each resolves it before calling here.
!
! The database's half duration is a per-station attribute; it is T_min/10
! by construction and therefore the same for every station of one mesh,
! and a database where it is not was assembled from two runs.

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(in) :: src
  double precision, intent(in) :: t0_req
  type(t_gf_taxis), intent(out) :: tax
  type(t_gf_stf), intent(out) :: stf
  integer, intent(out) :: ierr

  ! local parameters
  double precision :: hdur_db,t0,dt_sub
  integer :: ista

  if (.not. db%is_open) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_plan: database is not open')
    return
  endif
  if (db%nstations < 1) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_plan: the database holds no stations')
    return
  endif

  hdur_db = db%stations(1)%hdur
  do ista = 2,db%nstations
    if (db%stations(ista)%hdur /= hdur_db) then
      call gf_set_error(ierr,GF_ERR_MISMATCH, &
        'stations '//trim(db%stations(1)%id)//' and '//trim(db%stations(ista)%id)// &
        ' were built with different half durations')
      return
    endif
  enddo

  ! see gf_locate_source: the per-process globals must be this handle's
  call gf_init_shared_params(db,ierr)
  if (ierr /= GF_OK) return

  ! t0_req and hdur both reach the kernel widths and the padding count; a
  ! NaN there produces an array bound, not a NaN, so it is screened here
  if (.not. (gf_is_finite(t0_req) .and. gf_is_finite(src%hdur))) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_plan: t0 or the source half duration is not finite')
    return
  endif

  if (t0_req < 0.d0) then
    call gf_set_error(ierr,GF_ERR_ARG, &
      'gf_seis_plan: t0_req must be a resolved start time; call gf_default_t0 for specfem''s rule')
    return
  endif
  t0 = t0_req

  dt_sub = db%dt*dble(db%subsample_step)

  call gf_stf_plan(src%source_type,src%force_stf,src%hdur,hdur_db,dt_sub,GF_STF_TRUNC,stf,ierr)
  if (ierr /= GF_OK) return

  call gf_taxis_plan(db%nt_subsampled,db%dt,db%subsample_step,db%t0,t0,tax,ierr)

  end subroutine gf_seis_plan

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_seis_geometry(db,src,loc,itypsokern,geom,ierr)

! everything about the source's position that does not depend on the station
!
! The three extraction routines each built their own copy of this, in the
! same order, from the same inputs. `itypsokern` selects how much is needed:
! 0 gives the interpolation weights and, for a force source, the direction;
! 1 adds the strain weights and the rotated moment tensor; 2 adds the
! second-derivative tables and the geographic Jacobian.
!
! `dw` comes from gf_interp_weights_deriv (lagrange_any, hand-unrolled for
! NGLL = 5 with a fixed association order, lagrange_poly.f90:67-81) and `ddw`
! from gf_interp_weights_deriv2 (lagrange_any_2nd, a generic loop
! accumulating in a different order, :245). The two are kept in separate
! variables rather than sharing one table.
!
! Measured, not assumed: building dw from deriv2's tables instead leaves
! every reference output bitwise unchanged on both examples, so the two
! routines do agree to the last bit there. The separation is kept anyway
! because nothing guarantees that for every xi, and because the code this
! replaces got the same effect by overwriting hxi/hpxi in place *after* dw
! had been built -- correct, but resting on the order of two statements
! twenty lines apart.

  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(in) :: src
  type(t_gf_location), intent(in) :: loc
  integer, intent(in) :: itypsokern
  type(t_gf_seis_geom), intent(out) :: geom
  integer, intent(out) :: ierr

  ! local parameters
  double precision, dimension(NGLLX) :: hpxi
  double precision, dimension(NGLLY) :: hpeta
  double precision, dimension(NGLLZ) :: hpgam
  ! the deriv2 tables, deliberately separate from geom%hxi -- see above
  double precision, dimension(NGLLX) :: hxi2,hpxi2,hppxi
  double precision, dimension(NGLLY) :: heta2,hpeta2,hppeta
  double precision, dimension(NGLLZ) :: hgam2,hpgam2,hppgam
  double precision, dimension(1) :: no_spline
  double precision :: elevation,delev_dlat,delev_dlon
  integer :: nspl_use

  geom%want_loc = (itypsokern == 2)

  call gf_interp_weights_deriv(loc%xi,loc%eta,loc%gamma, &
                               geom%hxi,hpxi,geom%heta,hpeta,geom%hgam,hpgam)

  if (src%source_type == GF_SRC_FORCE) then
    call gf_force_direction(src,loc%nu,geom%fhat,ierr)
    if (ierr /= GF_OK) return
  else
    call gf_strain_dweights(geom%hxi,hpxi,geom%heta,hpeta,geom%hgam,hpgam,loc%jinv,geom%dw)
    call gf_rotate_moment_tensor(loc%theta,loc%phi,src%moment_tensor,geom%m_cart)
  endif

  if (geom%want_loc) then
    ! the differentiated weight table, the rotation's derivative and the
    ! geographic map's derivative, once: none depends on the station
    call gf_interp_weights_deriv2(loc%xi,loc%eta,loc%gamma,hxi2,hpxi2,hppxi,heta2,hpeta2,hppeta, &
                                  hgam2,hpgam2,hppgam)
    call gf_strain_ddweights(hxi2,hpxi2,hppxi,heta2,hpeta2,hppeta,hgam2,hpgam2,hppgam, &
                             loc%jinv,loc%djinv,geom%ddw)

    call gf_rotate_moment_tensor_deriv(loc%theta,loc%phi,src%moment_tensor, &
                                       geom%dm_dtheta,geom%dm_dphi)

    ! the elevation and its gradient, as gf_locate evaluated the elevation
    elevation = 0.d0
    delev_dlat = 0.d0
    delev_dlon = 0.d0
    if (db%topography) then
      call gf_topo_elevation(db,src%latitude,src%longitude,elevation)
      call gf_topo_gradient(db,src%latitude,src%longitude,delev_dlat,delev_dlon)
    endif

    ! the spline arrays only exist when the database has ELLIPTICITY set;
    ! an unallocated allocatable cannot be passed, so the no-spline case
    ! passes a length-one dummy instead
    if (db%ellipticity) then
      call gf_geographic_jacobian(src%latitude,src%longitude,src%depth,db%ellipticity, &
                                  elevation,delev_dlat,delev_dlon, &
                                  db%nspl,db%rspl,db%ellipicity_spline,db%ellipicity_spline2, &
                                  db%R_PLANET,geom%dxds,geom%dtheta_dlat,geom%dphi_dlon,ierr)
    else
      no_spline(1) = 0.d0
      nspl_use = 0
      call gf_geographic_jacobian(src%latitude,src%longitude,src%depth,db%ellipticity, &
                                  elevation,delev_dlat,delev_dlon, &
                                  nspl_use,no_spline,no_spline,no_spline, &
                                  db%R_PLANET,geom%dxds,geom%dtheta_dlat,geom%dphi_dlon,ierr)
    endif
    if (ierr /= GF_OK) return
  endif

  ierr = GF_OK

  end subroutine gf_seis_geometry

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_seis_work_init(db,src,tax,stf,itypsokern,swork,ierr)

! the buffers one extraction needs, sized from the plan

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(in) :: src
  type(t_gf_taxis), intent(in) :: tax
  type(t_gf_stf), intent(in) :: stf
  integer, intent(in) :: itypsokern
  type(t_gf_seis_work), intent(inout) :: swork
  integer, intent(out) :: ierr

  ! local parameters
  integer :: ier,nt_db,nt

  call gf_seis_work_free(swork)

  nt_db = db%nt_subsampled
  nt = tax%nt

  allocate(swork%displ(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ,nt_db),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate the element displacement buffer')
    return
  endif

  if (src%source_type == GF_SRC_FORCE) then
    allocate(swork%g(GF_NCOMP,GF_NCOMP,nt_db),stat=ier)
  else
    allocate(swork%eps(GF_VOIGT,GF_NCOMP,nt_db),stat=ier)
  endif
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate the interpolated field buffer')
    return
  endif

  if (itypsokern >= 1) then
    allocate(swork%dpm(GF_NDP_MT,GF_NCOMP,nt),stat=ier)
    if (ier /= 0) then
      call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate the moment-tensor partials buffer')
      return
    endif
  endif

  if (itypsokern == 2) then
    allocate(swork%deps(GF_VOIGT,GF_NCOMP,nt_db,NDIM),swork%dpl(3,GF_NCOMP,nt), &
             swork%dp10(nt),stat=ier)
    if (ier /= 0) then
      call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate the strain gradient buffer')
      return
    endif
  endif

  call gf_stf_work_init(stf,tax,swork%stfw,ierr)

  end subroutine gf_seis_work_init

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_seis_work_free(swork)

! releases the buffers; safe on an object that was never initialised

  implicit none

  type(t_gf_seis_work), intent(inout) :: swork

  if (allocated(swork%displ)) deallocate(swork%displ)
  if (allocated(swork%g)) deallocate(swork%g)
  if (allocated(swork%eps)) deallocate(swork%eps)
  if (allocated(swork%deps)) deallocate(swork%deps)
  if (allocated(swork%dpm)) deallocate(swork%dpm)
  if (allocated(swork%dpl)) deallocate(swork%dpl)
  if (allocated(swork%dp10)) deallocate(swork%dp10)
  call gf_stf_work_free(swork%stfw)

  end subroutine gf_seis_work_free

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_seis_station(db,src,loc,tax,stf,geom,itypsokern,ista,swork, &
                             seis,ndp,dp,onset,ierr)

! one station's traces, from the element block already in swork%displ
!
! Interpolate or differentiate, contract, convert, and -- when partials are
! wanted -- the three partial families. The caller has read the block and
! built `geom`; everything here depends on the station.
!
! `seis`, `dp` and `onset` are intent(inout): onset(ista) accumulates the
! worst component with max(), and the zeroing belongs to gf_seis, which owns
! the whole array.

  use constants, only: NDIM

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(in) :: src
  type(t_gf_location), intent(in) :: loc
  type(t_gf_taxis), intent(in) :: tax
  type(t_gf_stf), intent(in) :: stf
  type(t_gf_seis_geom), intent(in) :: geom
  integer, intent(in) :: itypsokern,ista,ndp
  type(t_gf_seis_work), intent(inout) :: swork
  double precision, dimension(db%nstations,GF_NCOMP,tax%nt), intent(inout) :: seis
  double precision, dimension(ndp,db%nstations,GF_NCOMP,tax%nt), intent(inout) :: dp
  double precision, dimension(db%nstations), intent(inout) :: onset
  integer, intent(out) :: ierr

  ! local parameters
  double precision :: scale_amp,ratio,wsum_raw
  integer :: it,icomp,idisp,ip,b,nt_db,nt,nbefore

  nt_db = db%nt_subsampled
  nt = tax%nt

  !--- the field at the source -------------------------------------------

  if (src%source_type == GF_SRC_FORCE) then
    ! g(a,d,t): the interpolated field at the source, still carrying both
    ! the station-force index a and the displacement index d
    call gf_interp_trace(swork%displ,geom%hxi,geom%heta,geom%hgam,nt_db,swork%g)
  else
    call gf_strain_trace(swork%displ,geom%dw,nt_db,swork%eps)
    if (itypsokern == 2) then
      ! d eps / d xi_b: the same kernel with the differentiated table
      do b = 1,NDIM
        call gf_strain_trace(swork%displ,geom%ddw(:,:,:,:,b),nt_db,swork%deps(:,:,:,b))
      enddo
    endif
  endif

  if (db%stations(ista)%factor_force_source == 0.d0) then
    call gf_set_error(ierr,GF_ERR_ARG, &
      'station '//trim(db%stations(ista)%id)//' has factor_force_source = 0')
    return
  endif

  ! two derivations, two statements: see the amplitude section of the module
  ! header. They are not one expression because they are not one idea.
  if (src%source_type == GF_SRC_FORCE) then
    scale_amp = src%factor_force_source / db%stations(ista)%factor_force_source
  else
    scale_amp = 1.d0 / db%stations(ista)%factor_force_source
  endif

  !--- the seismogram, component by component ----------------------------

  do icomp = 1,GF_NCOMP

    if (src%source_type == GF_SRC_FORCE) then
      do it = 1,nt_db
        swork%stfw%trace(it) = 0.d0
        do idisp = 1,GF_NCOMP
          swork%stfw%trace(it) = swork%stfw%trace(it) + geom%fhat(idisp)*swork%g(icomp,idisp,it)
        enddo
        swork%stfw%trace(it) = scale_amp * swork%stfw%trace(it)
      enddo
    else
      do it = 1,nt_db
        call gf_moment_contract(geom%m_cart,swork%eps(:,icomp,it),swork%stfw%trace(it))
        swork%stfw%trace(it) = scale_amp * swork%stfw%trace(it)
      enddo
    endif

    call gf_stf_onset(swork%stfw%trace,nt_db,swork%stfw%tdb,stf%hdur_db,ratio,nbefore)
    onset(ista) = max(onset(ista),ratio)

    ! extend, then convert: the kernel reaches into the padded region. For a
    ! moment tensor this is also the Gaussian -> Heaviside step.
    call gf_stf_convert(stf,tax,swork%stfw,ierr)
    if (ierr /= GF_OK) return

    do it = 1,nt
      seis(ista,icomp,it) = swork%stfw%y(it)
    enddo

    !--- the centroid-time partial, from the same padded trace -----------

    if (itypsokern == 2) then
      call gf_partials_time(swork%stfw%xpad,nt,tax%dt_sub,stf,swork%dp10,wsum_raw,ierr)
      if (ierr /= GF_OK) return
      do it = 1,nt
        dp(GF_DP_TIM,ista,icomp,it) = swork%dp10(it)
      enddo
    endif

  enddo

  !--- the moment-tensor partials, from the same strain -------------------

  if (itypsokern >= 1) then
    call gf_partials_mt(swork%eps,nt_db,loc%theta,loc%phi,scale_amp/src%scale_moment, &
                        tax,stf,swork%stfw,swork%dpm,ierr)
    if (ierr /= GF_OK) return

    do it = 1,nt
      do icomp = 1,GF_NCOMP
        do ip = 1,GF_NDP_MT
          dp(ip,ista,icomp,it) = swork%dpm(ip,icomp,it)
        enddo
      enddo
    enddo
  endif

  !--- the centroid-position partials --------------------------------------

  if (itypsokern == 2) then
    call gf_partials_loc(swork%eps,swork%deps,nt_db,geom%m_cart,geom%dm_dtheta,geom%dm_dphi, &
                         geom%dtheta_dlat,geom%dphi_dlon, &
                         loc%jinv,geom%dxds,scale_amp,tax,stf,swork%stfw,swork%dpl,ierr)
    if (ierr /= GF_OK) return
    do it = 1,nt
      do icomp = 1,GF_NCOMP
        do ip = 1,3
          dp(GF_DP_LAT+ip-1,ista,icomp,it) = swork%dpl(ip,icomp,it)
        enddo
      enddo
    enddo
  endif

  ierr = GF_OK

  end subroutine gf_seis_station

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_seis(db,src,loc,tax,stf,itypsokern,ndp,seis,dp,t,onset,ierr)

! seismograms, and optionally their partial derivatives, for either kind of
! source on a planned axis
!
!   itypsokern = 0   seismograms only; ndp = 0 and dp is a zero-sized array
!              = 1   plus the six moment-tensor partials
!              = 2   plus latitude, longitude, depth and centroid time
!
! `ndp` must be what gf_partials_ndp returns for `itypsokern`; it is the
! first extent of `dp` and is checked rather than inferred, so that a caller
! who sized its array from the wrong kind is told rather than silently given
! the wrong slots.
!
! `seis` is (nstations, 3 components N/E/Z, tax%nt) in metres, on the output
! axis `t`. `onset(ista)` is the worst silence-before-the-record ratio over
! the station's three components; see gf_stf_onset.
!
! The element is read once and the stations looped inside it, because the
! per-(element,station) array is the large object here -- 21 MB in the
! shipped global example. Nothing here holds two of them, and raw
! displacement is never handed back to a caller.

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(in) :: src
  type(t_gf_location), intent(in) :: loc
  type(t_gf_taxis), intent(in) :: tax
  type(t_gf_stf), intent(in) :: stf
  integer, intent(in) :: itypsokern,ndp
  double precision, dimension(db%nstations,GF_NCOMP,tax%nt), intent(out) :: seis
  double precision, dimension(ndp,db%nstations,GF_NCOMP,tax%nt), intent(out) :: dp
  double precision, dimension(tax%nt), intent(out) :: t
  double precision, dimension(db%nstations), intent(out) :: onset
  integer, intent(out) :: ierr

  ! local parameters
  type(t_gf_seis_geom) :: geom
  type(t_gf_seis_work) :: swork
  integer :: ista,nt_db,nt,ndp_want,ier

  seis(:,:,:) = 0.d0
  dp(:,:,:,:) = 0.d0
  onset(:) = 0.d0

  !--- what this request must satisfy -------------------------------------

  if (.not. db%is_open) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis: database is not open')
    return
  endif
  if (loc%ielem < 1) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis: the source has not been located')
    return
  endif
  if (db%nstations < 1) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis: the database holds no stations')
    return
  endif
  if (src%source_type /= GF_SRC_FORCE .and. src%source_type /= GF_SRC_CMT) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis: the source has no type set')
    return
  endif
  if (tax%nt_db /= db%nt_subsampled .or. tax%nt < tax%nt_db) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis: the time axis was not planned for this database')
    return
  endif

  ! the kind, and the array size it implies
  call gf_partials_ndp(itypsokern,ndp_want,ierr)
  if (ierr /= GF_OK) return
  if (ndp /= ndp_want) then
    call gf_set_error(ierr,GF_ERR_ARG, &
      'gf_seis: dp is not sized for this kind of partial; ask gf_partials_ndp for the count')
    return
  endif

  if (itypsokern > 0) then
    if (src%source_type /= GF_SRC_CMT) then
      call gf_set_error(ierr,GF_ERR_ARG, &
        'gf_seis: partial derivatives are defined for a moment-tensor source only')
      return
    endif
    if (src%scale_moment <= 0.d0) then
      call gf_set_error(ierr,GF_ERR_ARG,'gf_seis: the source carries no moment scale')
      return
    endif
  endif

  ! a force source plans a Gaussian; only the moment-tensor route integrates
  if (src%source_type == GF_SRC_CMT .and. stf%kind_stf /= GF_STF_HEAVI) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis: the plan is not a Heaviside conversion')
    return
  endif

  if (itypsokern == 2 .and. db%topography .and. .not. db%topo_loaded) then
    call gf_set_error(ierr,GF_ERR_ARG, &
      'gf_seis: the topography grid is not loaded; locate the source first')
    return
  endif

  !--- everything that does not depend on the station ---------------------

  nt_db = db%nt_subsampled
  nt = tax%nt

  call gf_taxis_times(tax,nt,t)

  call gf_seis_geometry(db,src,loc,itypsokern,geom,ierr)
  if (ierr /= GF_OK) return

  call gf_seis_work_init(db,src,tax,stf,itypsokern,swork,ier)
  if (ier /= GF_OK) then
    ierr = ier
    goto 99
  endif

  ! the database's own axis, for the onset check
  call gf_time_axis(db,nt_db,swork%stfw%tdb)

  !--- and the stations ----------------------------------------------------

  do ista = 1,db%nstations

    call gf_read_element_displ(db,loc%ielem,ista,swork%displ,ierr)
    if (ierr /= GF_OK) goto 99

    call gf_seis_station(db,src,loc,tax,stf,geom,itypsokern,ista,swork, &
                         seis,ndp,dp,onset,ierr)
    if (ierr /= GF_OK) goto 99

  enddo

  ierr = GF_OK

99 continue
  call gf_seis_work_free(swork)

  end subroutine gf_seis

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_write_seis(db,src,loc,tax,stf,seis,t,onset,outdir,ierr)

! writes one ASCII file per station
!
! The format is a maintained interface, not a debugging convenience: the
! comparison harness (utils/green_function/gf_compare.py) reads it, and
! Stage 9 diffs against it. Keep the '#' header block, the fixed-order
! `# key : values` lines that describe the conversion and the axis, and the
! four columns.

  use constants, only: MAX_STRING_LEN

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(in) :: src
  type(t_gf_location), intent(in) :: loc
  type(t_gf_taxis), intent(in) :: tax
  type(t_gf_stf), intent(in) :: stf
  double precision, dimension(db%nstations,GF_NCOMP,tax%nt), intent(in) :: seis
  double precision, dimension(tax%nt), intent(in) :: t
  double precision, dimension(db%nstations), intent(in) :: onset
  character(len=*), intent(in) :: outdir
  integer, intent(out) :: ierr

  ! local parameters
  character(len=MAX_STRING_LEN) :: filename
  integer :: ista,it,iout,ios
  logical :: exists

  call gf_dir_exists(outdir,exists)
  if (.not. exists) then
    call gf_set_error(ierr,GF_ERR_ARG,'output directory does not exist: '//trim(outdir))
    return
  endif

  do ista = 1,db%nstations

    filename = trim(outdir)//'/'//trim(db%stations(ista)%id)//'.gf3d.txt'

    open(newunit=iout,file=trim(filename),status='replace',action='write',iostat=ios)
    if (ios /= 0) then
      call gf_set_error(ierr,GF_ERR_IO,'could not open for writing: '//trim(filename))
      return
    endif

    call gf_write_header(db,src,loc,ista,iout,tax,stf,onset(ista))
    write(iout,'(a)') '# columns: t[s]  '//GF_COMP_NAME(1)//'[m]  ' &
                      //GF_COMP_NAME(2)//'[m]  '//GF_COMP_NAME(3)//'[m]'

    do it = 1,tax%nt
      write(iout,'(4es24.16)') t(it),seis(ista,1,it),seis(ista,2,it),seis(ista,3,it)
    enddo

    close(iout)

  enddo

  ierr = GF_OK

  end subroutine gf_write_seis

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_write_partials(db,src,loc,tax,stf,ndp,dp,t,outdir,ierr)

! writes one ASCII file of partial derivatives per station
!
! <NET>.<STA>.partials.txt, beside the seismogram: the same provenance and
! conversion header as gf_write_seis, a `# partials` line naming the slots
! and a `# units` line, then `t` followed by dp(ip,comp) for comp = N,E,Z
! and ip = 1..ndp, ip varying fastest -- 1 + 3 ndp columns. Like the
! seismogram format this is an interface, read by the tests and by Stage
! 9's binding checks.

  use constants, only: MAX_STRING_LEN

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(in) :: src
  type(t_gf_location), intent(in) :: loc
  type(t_gf_taxis), intent(in) :: tax
  type(t_gf_stf), intent(in) :: stf
  integer, intent(in) :: ndp
  double precision, dimension(ndp,db%nstations,GF_NCOMP,tax%nt), intent(in) :: dp
  double precision, dimension(tax%nt), intent(in) :: t
  character(len=*), intent(in) :: outdir
  integer, intent(out) :: ierr

  ! local parameters
  character(len=MAX_STRING_LEN) :: filename,line
  integer :: ista,it,ip,icomp,iout,ios
  logical :: exists

  if (ndp < 1 .or. ndp > GF_NDP_LOC) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_write_partials: ndp must be between 1 and 10')
    return
  endif

  call gf_dir_exists(outdir,exists)
  if (.not. exists) then
    call gf_set_error(ierr,GF_ERR_ARG,'output directory does not exist: '//trim(outdir))
    return
  endif

  do ista = 1,db%nstations

    filename = trim(outdir)//'/'//trim(db%stations(ista)%id)//'.partials.txt'

    open(newunit=iout,file=trim(filename),status='replace',action='write',iostat=ios)
    if (ios /= 0) then
      call gf_set_error(ierr,GF_ERR_IO,'could not open for writing: '//trim(filename))
      return
    endif

    call gf_write_header(db,src,loc,ista,iout,tax,stf)

    line = ''
    do ip = 1,ndp
      line = trim(line)//' '//GF_DP_NAME(ip)
    enddo
    write(iout,'(a,i0,a)') '# partials   : ',ndp,trim(line)
    line = ''
    do ip = 1,ndp
      line = trim(line)//' '//trim(GF_DP_UNIT(ip))
    enddo
    write(iout,'(a)') '# units      :'//trim(line)
    write(iout,'(a)') '# columns: t[s]  then dp(ip,comp) for comp = '//GF_COMP_NAME(1)//','// &
                      GF_COMP_NAME(2)//','//GF_COMP_NAME(3)//' and ip = 1..ndp, ip varying fastest'

    do it = 1,tax%nt
      write(iout,'(*(es24.16))') t(it),((dp(ip,ista,icomp,it),ip = 1,ndp),icomp = 1,GF_NCOMP)
    enddo

    close(iout)

  enddo

  ierr = GF_OK

  end subroutine gf_write_partials

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_write_dump(db,src,loc,outdir,station_filter,ierr)

! writes the interpolated displacement, before any contraction
!
! This is the honest replacement for what a Fortran-versus-Python differ used
! to provide. The manufactured-solution unit tests pin the *operators* on
! synthetic input; they cannot show that gf_element_io read the real HDF5
! layout in the right index order, nor that the nu convention is right. Both
! are forward-run-only and both fail non-locally at Stage 5. This mode
! answers "which operator is wrong", not "is the answer right" -- a debugging
! instrument, not a second implementation.
!
! Writes <NET>.<STA>.dump.txt with nine columns after the time: g(a,d), with
! a the station force component (N,E,Z) and d the Cartesian displacement
! component (x,y,z), d varying fastest.
!
! For a moment-tensor source it also writes <NET>.<STA>.strain.txt with the
! Voigt strain (6 slots x 3 force components) and the contracted trace
! *before* the source time function conversion -- the two intermediates
! between the raw field and the seismogram, so a disagreement can be
! attributed to the geometric chain, the contraction, or the conversion
! separately.
!
! Everything here is on the database's own axis: these are the stored
! quantities, and the output axis is a property of the seismogram.
!
! `station_filter` of '' dumps every station; 21 MB is read per station, so
! naming one is usually what is wanted.

  use gf_par, only: GF_SRC_CMT

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,NGLLX,NGLLY,NGLLZ,NDIM

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(in) :: src
  type(t_gf_location), intent(in) :: loc
  character(len=*), intent(in) :: outdir
  character(len=*), intent(in) :: station_filter
  integer, intent(out) :: ierr

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: displ
  double precision, dimension(:,:,:), allocatable :: g,eps
  double precision, dimension(:,:), allocatable :: pre_stf
  double precision, dimension(:), allocatable :: t
  type(t_gf_seis_geom) :: geom
  character(len=MAX_STRING_LEN) :: filename
  integer :: ista,it,ia,id,iv,iout,ios,ier,nt,nwritten
  logical :: exists,do_strain

  if (.not. db%is_open) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_write_dump: database is not open')
    return
  endif
  if (loc%ielem < 1) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_write_dump: the source has not been located')
    return
  endif

  call gf_dir_exists(outdir,exists)
  if (.not. exists) then
    call gf_set_error(ierr,GF_ERR_ARG,'output directory does not exist: '//trim(outdir))
    return
  endif

  nt = db%nt_subsampled

  do_strain = (src%source_type == GF_SRC_CMT)

  ! the same weights, strain table and rotated moment tensor the seismogram
  ! path uses, from the same routine, so that a disagreement between --dump
  ! and --seis cannot come from the geometry
  call gf_seis_geometry(db,src,loc,0,geom,ierr)
  if (ierr /= GF_OK) return

  allocate(displ(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ,nt),g(GF_NCOMP,GF_NCOMP,nt), &
           eps(GF_VOIGT,GF_NCOMP,nt),pre_stf(GF_NCOMP,nt),t(nt),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate the element displacement buffer')
    return
  endif

  call gf_time_axis(db,nt,t)

  nwritten = 0

  do ista = 1,db%nstations

    if (len_trim(station_filter) > 0) then
      if (trim(db%stations(ista)%id) /= trim(station_filter)) cycle
    endif

    call gf_read_element_displ(db,loc%ielem,ista,displ,ierr)
    if (ierr /= GF_OK) goto 99

    call gf_interp_trace(displ,geom%hxi,geom%heta,geom%hgam,nt,g)

    filename = trim(outdir)//'/'//trim(db%stations(ista)%id)//'.dump.txt'

    open(newunit=iout,file=trim(filename),status='replace',action='write',iostat=ios)
    if (ios /= 0) then
      call gf_set_error(ierr,GF_ERR_IO,'could not open for writing: '//trim(filename))
      goto 99
    endif

    call gf_write_header(db,src,loc,ista,iout)
    write(iout,'(a)') '# interpolated displacement, before the source contraction'
    write(iout,'(a)') '# columns: t[s]  then g(a,d) for a = N,E,Z (force at the station)'
    write(iout,'(a)') '#          and d = x,y,z (Cartesian displacement), d varying fastest'

    do it = 1,nt
      write(iout,'(10es24.16)') t(it),((g(ia,id,it),id = 1,GF_NCOMP),ia = 1,GF_NCOMP)
    enddo

    close(iout)

    !--- the moment-tensor intermediates -------------------------------

    if (do_strain) then

      call gf_strain_trace(displ,geom%dw,nt,eps)

      ! The seismogram path refuses this rather than dividing by zero
      ! (gf_seis_station); --dump used to divide anyway.
      if (db%stations(ista)%factor_force_source == 0.d0) then
        call gf_set_error(ierr,GF_ERR_ARG, &
          'station '//trim(db%stations(ista)%id)//' has factor_force_source = 0')
        goto 99
      endif

      ! In the seismogram's units, so the dumped trace differs from it only
      ! by the source time function conversion.
      !
      ! Note this divides where gf_seis_station multiplies by a reciprocal
      ! it formed once: x/f and x*(1/f) differ in the last bit whenever 1/f
      ! is inexact, so the two are not bitwise equal and are not meant to be.
      ! The .strain.txt files are in the reference set, so making them agree
      ! is an output change with its own entry in
      ! allowed_output_changes.md -- not something to tidy in passing.
      do ia = 1,GF_NCOMP
        do it = 1,nt
          call gf_moment_contract(geom%m_cart,eps(:,ia,it),pre_stf(ia,it))
          pre_stf(ia,it) = pre_stf(ia,it) / db%stations(ista)%factor_force_source
        enddo
      enddo

      filename = trim(outdir)//'/'//trim(db%stations(ista)%id)//'.strain.txt'

      open(newunit=iout,file=trim(filename),status='replace',action='write',iostat=ios)
      if (ios /= 0) then
        call gf_set_error(ierr,GF_ERR_IO,'could not open for writing: '//trim(filename))
        goto 99
      endif

      call gf_write_header(db,src,loc,ista,iout)
      write(iout,'(a)') '# Voigt strain of the reciprocal field, and the moment-tensor'
      write(iout,'(a)') '# contraction before the time integration'
      write(iout,'(a)') '# columns: t[s]  then eps(v,a) for a = N,E,Z (force at the station)'
      write(iout,'(a)') '#          and v = xx,yy,zz,xy,xz,yz, v varying fastest (18 columns)'
      write(iout,'(a)') '#          then the pre-integration trace for a = N,E,Z (3 columns)'

      do it = 1,nt
        write(iout,'(22es24.16)') t(it), &
          ((eps(iv,ia,it),iv = 1,GF_VOIGT),ia = 1,GF_NCOMP), &
          (pre_stf(ia,it),ia = 1,GF_NCOMP)
      enddo

      close(iout)

    endif

    nwritten = nwritten + 1

  enddo

  if (nwritten == 0) then
    call gf_set_error(ierr,GF_ERR_ARG,'no station matched: '//trim(station_filter))
    goto 99
  endif

  ierr = GF_OK

99 continue
  if (allocated(displ)) deallocate(displ)
  if (allocated(g)) deallocate(g)
  if (allocated(eps)) deallocate(eps)
  if (allocated(pre_stf)) deallocate(pre_stf)
  if (allocated(t)) deallocate(t)

  end subroutine gf_write_dump

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_write_header(db,src,loc,ista,iout,tax,stf,onset)

! the '#' provenance block shared by every output format
!
! Everything needed to reproduce the trace: which database, which source,
! which element and where in it, and the time axis. A trace whose provenance
! has to be reconstructed from the filename is a trace nobody can check.
!
! `tax`, `stf` and `onset` are optional because the three writers do not all
! have them. --dump runs before a plan exists and passes none; the partials
! file has no per-station onset. The order of the lines below is the order
! both formats have always emitted and is compared byte for byte by
! check_reference.sh, so an `if (present(...))` that moves a line is a
! defect even though every line is still written.

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(in) :: src
  type(t_gf_location), intent(in) :: loc
  integer, intent(in) :: ista
  integer, intent(in) :: iout
  type(t_gf_taxis), intent(in), optional :: tax
  type(t_gf_stf), intent(in), optional :: stf
  double precision, intent(in), optional :: onset

  write(iout,'(a)')  '# xgf3d '//GF3D_VERSION
  write(iout,'(a)')  '# database    : '//trim(db%path)
  write(iout,'(a)')  '# source      : '//trim(src%filename)
  write(iout,'(a,3es24.16)') '# source lat/lon/depth_km : ', &
                             src%latitude,src%longitude,src%depth
  write(iout,'(a)')  '# element     : '//loc%morton_hex
  write(iout,'(a,3es24.16)') '# xi,eta,gamma            : ',loc%xi,loc%eta,loc%gamma
  write(iout,'(a)')  '# station     : '//trim(db%stations(ista)%id)
  write(iout,'(a,3es24.16)') '# station lat/lon/depth_m : ', &
                             db%stations(ista)%latitude,db%stations(ista)%longitude, &
                             db%stations(ista)%depth
  write(iout,'(a,2es24.16,i12)') '# t0, dt_sub, nt          : ', &
                             db%t0,db%dt*dble(db%subsample_step),db%nt_subsampled

  ! the conversion and the axis, machine-parsable
  if (present(stf)) then
    write(iout,'(a,a)')               '# stf kind   : ',trim(gf_stf_kind_name(stf%kind_stf))
    write(iout,'(a,4es24.16)')        '# stf hdur   : ',stf%hdur_src,stf%hdur_target,stf%hdur_db,stf%hdur_corr
    write(iout,'(a,es24.16,i12,l4)')  '# stf kernel : ',stf%trunc,stf%khalf,stf%guard
  endif

  if (present(onset)) then
    write(iout,'(a,es24.16)')         '# stf onset  : ',onset
    if (present(stf)) write(iout,'(a,a)') '# stf note   : ',trim(stf%note)
  endif

  if (present(tax)) then
    write(iout,'(a,es24.16,4i12)')    '# axis       : ',tax%dt,tax%subsample_step,tax%nt_db,tax%npad,tax%nt
    write(iout,'(a,4es24.16)')        '# axis t0    : ',tax%t0_db,tax%t0_req,tax%t0,tax%t_first
  endif

  if (present(onset) .and. present(stf)) then
    write(iout,'(a,i0,a)')            '# edge       : trailing ',stf%khalf, &
                                      ' samples use zero-extended data; leading samples assume silence before the database start'
  endif

  end subroutine gf_write_header

  end module gf_seismograms
