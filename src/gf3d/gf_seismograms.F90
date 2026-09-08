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

  use gf_par, only: t_gfdb,t_gf_location,t_gf_source,t_gf_stf,t_gf_taxis,gf_set_error, &
                    GF_OK,GF_ERR_ARG,GF_ERR_ALLOC,GF_ERR_IO,GF_ERR_MISMATCH, &
                    GF_NCOMP,GF3D_VERSION,GF_STF_TRUNC,GF_STF_HEAVI, &
                    GF_SRC_FORCE,GF_SRC_CMT

  use gf_database, only: gf_dir_exists

  use gf_element_io, only: gf_read_element_displ

  use gf_interp, only: gf_interp_weights,gf_interp_weights_deriv,gf_interp_trace

  use gf_strain, only: gf_strain_dweights,gf_strain_trace,GF_VOIGT

  use gf_moment, only: gf_rotate_moment_tensor,gf_moment_contract

  use gf_stf, only: gf_stf_plan,gf_taxis_plan,gf_taxis_times,gf_pad_left,gf_cumsum, &
                    gf_stf_kernel,gf_stf_apply,gf_stf_onset,gf_stf_kind_name

  use gf_source, only: gf_force_direction

  use gf_partials, only: gf_partials_mt,GF_NDP_MT,GF_NDP_LOC,GF_DP_NAME,GF_DP_UNIT

  implicit none

  private

  public :: gf_time_axis
  public :: gf_seis_plan
  public :: gf_seis_force
  public :: gf_seis_cmt
  public :: gf_seis_cmt_partials
  public :: gf_seis
  public :: gf_write_seis
  public :: gf_write_partials
  public :: gf_write_dump

  ! component labels, in the stored force order (green_function_io.F90:38)
  character(len=1), dimension(GF_NCOMP), parameter :: GF_COMP_NAME = (/ 'N','E','Z' /)

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
! `t0_req` is the requested start time in seconds before the origin. A
! negative value asks for specfem's own rule for the forward run
! (setup_sources_receivers.f90:784-809): 1.5*hdur for a CMT and for the
! Gaussian and Heaviside force types, 1.2/f0 for a Ricker, 0 for a
! monochromatic force. tshift_src is zero for a single source, so that is
! the forward run's t0 exactly (SAC header b = -90 for the shipped
! CMTSOLUTION, -67.5 for the FORCESOLUTION).
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

  t0 = t0_req
  if (t0 < 0.d0) then
    select case (src%source_type)
    case (GF_SRC_CMT)
      t0 = 1.5d0*src%hdur
    case (GF_SRC_FORCE)
      select case (src%force_stf)
      case (1)
        ! Ricker: hdur holds the dominant frequency
        if (src%hdur <= 0.d0) then
          call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_plan: a Ricker force needs a positive f0')
          return
        endif
        t0 = 1.2d0/src%hdur
      case (3)
        t0 = 0.d0
      case default
        t0 = 1.5d0*src%hdur
      end select
    case default
      call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_plan: the source has no type set')
      return
    end select
  endif

  dt_sub = db%dt*dble(db%subsample_step)

  call gf_stf_plan(src%source_type,src%force_stf,src%hdur,hdur_db,dt_sub,GF_STF_TRUNC,stf,ierr)
  if (ierr /= GF_OK) return

  call gf_taxis_plan(db%nt_subsampled,db%dt,db%subsample_step,db%t0,t0,tax,ierr)

  end subroutine gf_seis_plan

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_seis_force(db,src,loc,tax,stf,seis,t,onset,ierr)

! seismograms at every station for a force source
!
! `seis` is (nstations, 3 components N/E/Z, tax%nt) in metres, on the
! output axis `t`, as the response to the source time function specfem
! would use for this FORCESOLUTION. `onset(ista)` is the worst
! silence-before-the-record ratio over the station's three components; see
! gf_stf_onset.
!
! The element is read once and the stations looped inside it, because the
! per-(element,station) array is the large object here -- 21 MB in the
! shipped global example. Nothing in this routine holds two of them, and raw
! displacement is never handed back to a caller.

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(in) :: src
  type(t_gf_location), intent(in) :: loc
  type(t_gf_taxis), intent(in) :: tax
  type(t_gf_stf), intent(in) :: stf
  double precision, dimension(db%nstations,GF_NCOMP,tax%nt), intent(out) :: seis
  double precision, dimension(tax%nt), intent(out) :: t
  double precision, dimension(db%nstations), intent(out) :: onset
  integer, intent(out) :: ierr

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: displ
  double precision, dimension(:,:,:), allocatable :: g
  double precision, dimension(:), allocatable :: trace,tdb,xpad,p,y,w
  double precision, dimension(NGLLX) :: hxi
  double precision, dimension(NGLLY) :: heta
  double precision, dimension(NGLLZ) :: hgam
  double precision, dimension(NDIM) :: fhat
  double precision :: scale_amp,ratio
  integer :: ista,it,icomp,idisp,ier,nt_db,nt,nbefore

  seis(:,:,:) = 0.d0
  onset(:) = 0.d0

  if (.not. db%is_open) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_force: database is not open')
    return
  endif
  if (loc%ielem < 1) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_force: the source has not been located')
    return
  endif
  if (db%nstations < 1) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_force: the database holds no stations')
    return
  endif
  if (tax%nt_db /= db%nt_subsampled .or. tax%nt < tax%nt_db) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_force: the time axis was not planned for this database')
    return
  endif

  nt_db = db%nt_subsampled
  nt = tax%nt

  call gf_taxis_times(tax,nt,t)

  call gf_interp_weights(loc%xi,loc%eta,loc%gamma,hxi,heta,hgam)

  call gf_force_direction(src,loc%nu,fhat,ierr)
  if (ierr /= GF_OK) return

  allocate(displ(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ,nt_db),g(GF_NCOMP,GF_NCOMP,nt_db), &
           trace(nt_db),tdb(nt_db),xpad(nt),p(0:nt),y(nt),w(-stf%khalf:stf%khalf),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate the element displacement buffer')
    return
  endif

  ! the database's own axis, for the onset check
  call gf_time_axis(db,nt_db,tdb)

  ! the conversion kernel, once: it does not depend on the station
  call gf_stf_kernel(stf,tax%dt_sub,w)

  do ista = 1,db%nstations

    call gf_read_element_displ(db,loc%ielem,ista,displ,ierr)
    if (ierr /= GF_OK) goto 99

    ! g(a,d,t): the interpolated field at the source, still carrying both
    ! the station-force index a and the displacement index d
    call gf_interp_trace(displ,hxi,heta,hgam,nt_db,g)

    if (db%stations(ista)%factor_force_source == 0.d0) then
      call gf_set_error(ierr,GF_ERR_ARG, &
        'station '//trim(db%stations(ista)%id)//' has factor_force_source = 0')
      goto 99
    endif
    scale_amp = src%factor_force_source / db%stations(ista)%factor_force_source

    do icomp = 1,GF_NCOMP

      do it = 1,nt_db
        trace(it) = 0.d0
        do idisp = 1,GF_NCOMP
          trace(it) = trace(it) + fhat(idisp)*g(icomp,idisp,it)
        enddo
        trace(it) = scale_amp * trace(it)
      enddo

      call gf_stf_onset(trace,nt_db,tdb,stf%hdur_db,ratio,nbefore)
      onset(ista) = max(onset(ista),ratio)

      ! extend, then convert: the kernel reaches into the padded region
      call gf_pad_left(trace,nt_db,tax%npad,xpad)
      call gf_cumsum(xpad,nt,p)
      call gf_stf_apply(stf,tax%dt_sub,w,p,xpad,nt,y)

      do it = 1,nt
        seis(ista,icomp,it) = y(it)
      enddo

    enddo

  enddo

  ierr = GF_OK

99 continue
  if (allocated(displ)) deallocate(displ)
  if (allocated(g)) deallocate(g)
  if (allocated(trace)) deallocate(trace)
  if (allocated(tdb)) deallocate(tdb)
  if (allocated(xpad)) deallocate(xpad)
  if (allocated(p)) deallocate(p)
  if (allocated(y)) deallocate(y)
  if (allocated(w)) deallocate(w)

  end subroutine gf_seis_force

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_seis_cmt(db,src,loc,tax,stf,seis,t,onset,ierr)

! seismograms at every station for a moment-tensor source
!
! `seis` is (nstations, 3 components N/E/Z, tax%nt) in metres, on the
! output axis `t`, as the response to specfem's quasi-Heaviside at the
! CMTSOLUTION's half duration -- what the forward run computes, lowpassed
! by the database's own Butterworth.
!
! Three steps, in this order:
!
!   1. the strain of the reciprocal field at the source, eps(6,3,nt);
!   2. contraction with the Cartesian moment tensor, giving the response to
!      the reciprocal run's Gaussian source time function;
!   3. one convolution with the quasi-Heaviside of the corrected width
!      (gf_stf.F90), which is both the Gaussian -> Heaviside conversion and
!      the change of half duration.
!
! Step 3 replaces Stage 4's cumulative trapezoid, whose integrator error on
! the 3.4 s regional grid is one percent at a 60 s period. Nothing here
! integrates.

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(in) :: src
  type(t_gf_location), intent(in) :: loc
  type(t_gf_taxis), intent(in) :: tax
  type(t_gf_stf), intent(in) :: stf
  double precision, dimension(db%nstations,GF_NCOMP,tax%nt), intent(out) :: seis
  double precision, dimension(tax%nt), intent(out) :: t
  double precision, dimension(db%nstations), intent(out) :: onset
  integer, intent(out) :: ierr

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: displ
  double precision, dimension(:,:,:), allocatable :: eps
  double precision, dimension(:), allocatable :: trace,tdb,xpad,p,y,w
  double precision, dimension(NGLLX) :: hxi,hpxi
  double precision, dimension(NGLLY) :: heta,hpeta
  double precision, dimension(NGLLZ) :: hgam,hpgam
  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM) :: dw
  double precision, dimension(NDIM,NDIM) :: m_cart
  double precision :: scale_amp,ratio
  integer :: ista,it,icomp,ier,nt_db,nt,nbefore

  seis(:,:,:) = 0.d0
  onset(:) = 0.d0

  if (.not. db%is_open) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_cmt: database is not open')
    return
  endif
  if (loc%ielem < 1) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_cmt: the source has not been located')
    return
  endif
  if (db%nstations < 1) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_cmt: the database holds no stations')
    return
  endif
  if (tax%nt_db /= db%nt_subsampled .or. tax%nt < tax%nt_db) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_cmt: the time axis was not planned for this database')
    return
  endif
  if (stf%kind_stf /= GF_STF_HEAVI) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_cmt: the plan is not a Heaviside conversion')
    return
  endif

  nt_db = db%nt_subsampled
  nt = tax%nt

  call gf_taxis_times(tax,nt,t)

  ! basis values and reference derivatives at the source, then the
  ! physical-space derivatives through the locator's inverse Jacobian
  call gf_interp_weights_deriv(loc%xi,loc%eta,loc%gamma,hxi,hpxi,heta,hpeta,hgam,hpgam)
  call gf_strain_dweights(hxi,hpxi,heta,hpeta,hgam,hpgam,loc%jinv,dw)

  ! the moment tensor, rotated once: it does not depend on the station
  call gf_rotate_moment_tensor(loc%theta,loc%phi,src%moment_tensor,m_cart)

  allocate(displ(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ,nt_db),eps(GF_VOIGT,GF_NCOMP,nt_db), &
           trace(nt_db),tdb(nt_db),xpad(nt),p(0:nt),y(nt),w(-stf%khalf:stf%khalf),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate the element displacement buffer')
    return
  endif

  call gf_time_axis(db,nt_db,tdb)

  call gf_stf_kernel(stf,tax%dt_sub,w)

  do ista = 1,db%nstations

    call gf_read_element_displ(db,loc%ielem,ista,displ,ierr)
    if (ierr /= GF_OK) goto 99

    call gf_strain_trace(displ,dw,nt_db,eps)

    ! per unit reciprocal force at the station; see the module header
    if (db%stations(ista)%factor_force_source == 0.d0) then
      call gf_set_error(ierr,GF_ERR_ARG, &
        'station '//trim(db%stations(ista)%id)//' has factor_force_source = 0')
      goto 99
    endif
    scale_amp = 1.d0 / db%stations(ista)%factor_force_source

    do icomp = 1,GF_NCOMP

      do it = 1,nt_db
        call gf_moment_contract(m_cart,eps(:,icomp,it),trace(it))
        trace(it) = scale_amp * trace(it)
      enddo

      call gf_stf_onset(trace,nt_db,tdb,stf%hdur_db,ratio,nbefore)
      onset(ista) = max(onset(ista),ratio)

      ! Gaussian response -> Heaviside response at the CMT's half duration,
      ! in one convolution on the extended axis
      call gf_pad_left(trace,nt_db,tax%npad,xpad)
      call gf_cumsum(xpad,nt,p)
      call gf_stf_apply(stf,tax%dt_sub,w,p,xpad,nt,y)

      do it = 1,nt
        seis(ista,icomp,it) = y(it)
      enddo

    enddo

  enddo

  ierr = GF_OK

99 continue
  if (allocated(displ)) deallocate(displ)
  if (allocated(eps)) deallocate(eps)
  if (allocated(trace)) deallocate(trace)
  if (allocated(tdb)) deallocate(tdb)
  if (allocated(xpad)) deallocate(xpad)
  if (allocated(p)) deallocate(p)
  if (allocated(y)) deallocate(y)
  if (allocated(w)) deallocate(w)

  end subroutine gf_seis_cmt

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_seis_cmt_partials(db,src,loc,tax,stf,itypsokern,ndp,seis,dp,t,onset,ierr)

! seismograms and their partial derivatives at every station, for a
! moment-tensor source
!
! `seis` is what gf_seis_cmt returns, produced by the same statement
! sequence on the same arrays, and `dp(ndp,nstations,3,tax%nt)` holds the
! partials in the order of GF_DP_NAME (gf_partials.F90), `ndp` being what
! gf_partials_ndp returns for `itypsokern`: 6 for the moment tensor, 10
! with the centroid position and time.
!
! gf_seis_cmt is left untouched on purpose. The comparison gate in
! EXAMPLES/ rests on its output, and identical statements in different
! contexts can round differently under a fast floating-point model (the
! lesson of tests/gf3d/test_gf_shape3D.f90), so the two are pinned equal by
! tests/gf3d/test_gf_partials_db.f90 rather than assumed equal by
! construction. Per station the partials come from the *same* strain trace
! the seismogram was contracted from: one element read, as before.
!
! Stage 6 fills slots 1..6. itypsokern = 2 -- slots 7..9 (Stage 8) and 10
! -- is refused until then rather than returned half empty.

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM

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
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: displ
  double precision, dimension(:,:,:), allocatable :: eps,dpm
  double precision, dimension(:), allocatable :: trace,tdb,xpad,p,y,w
  double precision, dimension(NGLLX) :: hxi,hpxi
  double precision, dimension(NGLLY) :: heta,hpeta
  double precision, dimension(NGLLZ) :: hgam,hpgam
  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM) :: dw
  double precision, dimension(NDIM,NDIM) :: m_cart
  double precision :: scale_amp,ratio
  integer :: ista,it,icomp,ier,nt_db,nt,nbefore,ip

  seis(:,:,:) = 0.d0
  dp(:,:,:,:) = 0.d0
  onset(:) = 0.d0

  if (.not. db%is_open) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_cmt_partials: database is not open')
    return
  endif
  if (loc%ielem < 1) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_cmt_partials: the source has not been located')
    return
  endif
  if (db%nstations < 1) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_cmt_partials: the database holds no stations')
    return
  endif
  if (src%source_type /= GF_SRC_CMT) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_cmt_partials: partials are defined for moment-tensor sources')
    return
  endif
  if (tax%nt_db /= db%nt_subsampled .or. tax%nt < tax%nt_db) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_cmt_partials: the time axis was not planned for this database')
    return
  endif
  if (stf%kind_stf /= GF_STF_HEAVI) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_cmt_partials: the plan is not a Heaviside conversion')
    return
  endif
  select case (itypsokern)
  case (1)
    if (ndp /= GF_NDP_MT) then
      call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_cmt_partials: itypsokern = 1 returns 6 partials')
      return
    endif
  case (2)
    call gf_set_error(ierr,GF_ERR_ARG, &
      'gf_seis_cmt_partials: the centroid partials (itypsokern = 2) arrive with Stage 8')
    return
  case default
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_cmt_partials: itypsokern must be 1 or 2')
    return
  end select
  if (src%scale_moment <= 0.d0) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis_cmt_partials: the source carries no moment scale')
    return
  endif

  nt_db = db%nt_subsampled
  nt = tax%nt

  call gf_taxis_times(tax,nt,t)

  call gf_interp_weights_deriv(loc%xi,loc%eta,loc%gamma,hxi,hpxi,heta,hpeta,hgam,hpgam)
  call gf_strain_dweights(hxi,hpxi,heta,hpeta,hgam,hpgam,loc%jinv,dw)

  call gf_rotate_moment_tensor(loc%theta,loc%phi,src%moment_tensor,m_cart)

  allocate(displ(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ,nt_db),eps(GF_VOIGT,GF_NCOMP,nt_db), &
           trace(nt_db),tdb(nt_db),xpad(nt),p(0:nt),y(nt),w(-stf%khalf:stf%khalf), &
           dpm(GF_NDP_MT,GF_NCOMP,nt),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate the element displacement buffer')
    return
  endif

  call gf_time_axis(db,nt_db,tdb)

  call gf_stf_kernel(stf,tax%dt_sub,w)

  do ista = 1,db%nstations

    call gf_read_element_displ(db,loc%ielem,ista,displ,ierr)
    if (ierr /= GF_OK) goto 99

    call gf_strain_trace(displ,dw,nt_db,eps)

    if (db%stations(ista)%factor_force_source == 0.d0) then
      call gf_set_error(ierr,GF_ERR_ARG, &
        'station '//trim(db%stations(ista)%id)//' has factor_force_source = 0')
      goto 99
    endif
    scale_amp = 1.d0 / db%stations(ista)%factor_force_source

    !--- the seismogram: gf_seis_cmt's statements, verbatim ----------------

    do icomp = 1,GF_NCOMP

      do it = 1,nt_db
        call gf_moment_contract(m_cart,eps(:,icomp,it),trace(it))
        trace(it) = scale_amp * trace(it)
      enddo

      call gf_stf_onset(trace,nt_db,tdb,stf%hdur_db,ratio,nbefore)
      onset(ista) = max(onset(ista),ratio)

      call gf_pad_left(trace,nt_db,tax%npad,xpad)
      call gf_cumsum(xpad,nt,p)
      call gf_stf_apply(stf,tax%dt_sub,w,p,xpad,nt,y)

      do it = 1,nt
        seis(ista,icomp,it) = y(it)
      enddo

    enddo

    !--- the moment-tensor partials, from the same strain -------------------

    call gf_partials_mt(eps,nt_db,loc%theta,loc%phi,scale_amp/src%scale_moment, &
                        tax,stf,w,dpm,ierr)
    if (ierr /= GF_OK) goto 99

    do it = 1,nt
      do icomp = 1,GF_NCOMP
        do ip = 1,GF_NDP_MT
          dp(ip,ista,icomp,it) = dpm(ip,icomp,it)
        enddo
      enddo
    enddo

  enddo

  ierr = GF_OK

99 continue
  if (allocated(displ)) deallocate(displ)
  if (allocated(eps)) deallocate(eps)
  if (allocated(dpm)) deallocate(dpm)
  if (allocated(trace)) deallocate(trace)
  if (allocated(tdb)) deallocate(tdb)
  if (allocated(xpad)) deallocate(xpad)
  if (allocated(p)) deallocate(p)
  if (allocated(y)) deallocate(y)
  if (allocated(w)) deallocate(w)

  end subroutine gf_seis_cmt_partials

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_seis(db,src,loc,tax,stf,seis,t,onset,ierr)

! seismograms for either kind of source, on a planned axis

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(in) :: src
  type(t_gf_location), intent(in) :: loc
  type(t_gf_taxis), intent(in) :: tax
  type(t_gf_stf), intent(in) :: stf
  double precision, dimension(db%nstations,GF_NCOMP,tax%nt), intent(out) :: seis
  double precision, dimension(tax%nt), intent(out) :: t
  double precision, dimension(db%nstations), intent(out) :: onset
  integer, intent(out) :: ierr

  select case (src%source_type)
  case (GF_SRC_FORCE)
    call gf_seis_force(db,src,loc,tax,stf,seis,t,onset,ierr)
  case (GF_SRC_CMT)
    call gf_seis_cmt(db,src,loc,tax,stf,seis,t,onset,ierr)
  case default
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis: the source has no type set')
  end select

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

    call gf_write_header(db,src,loc,ista,iout)

    ! the conversion and the axis, machine-parsable
    write(iout,'(a,a)')               '# stf kind   : ',trim(gf_stf_kind_name(stf%kind_stf))
    write(iout,'(a,4es24.16)')        '# stf hdur   : ',stf%hdur_src,stf%hdur_target,stf%hdur_db,stf%hdur_corr
    write(iout,'(a,es24.16,i12,l4)')  '# stf kernel : ',stf%trunc,stf%khalf,stf%guard
    write(iout,'(a,es24.16)')         '# stf onset  : ',onset(ista)
    write(iout,'(a,a)')               '# stf note   : ',trim(stf%note)
    write(iout,'(a,es24.16,4i12)')    '# axis       : ',tax%dt,tax%subsample_step,tax%nt_db,tax%npad,tax%nt
    write(iout,'(a,4es24.16)')        '# axis t0    : ',tax%t0_db,tax%t0_req,tax%t0,tax%t_first
    write(iout,'(a,i0,a)')            '# edge       : trailing ',stf%khalf, &
                                      ' samples use zero-extended data; leading samples assume silence before the database start'
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

    call gf_write_header(db,src,loc,ista,iout)

    write(iout,'(a,a)')               '# stf kind   : ',trim(gf_stf_kind_name(stf%kind_stf))
    write(iout,'(a,4es24.16)')        '# stf hdur   : ',stf%hdur_src,stf%hdur_target,stf%hdur_db,stf%hdur_corr
    write(iout,'(a,es24.16,i12,l4)')  '# stf kernel : ',stf%trunc,stf%khalf,stf%guard
    write(iout,'(a,es24.16,4i12)')    '# axis       : ',tax%dt,tax%subsample_step,tax%nt_db,tax%npad,tax%nt
    write(iout,'(a,4es24.16)')        '# axis t0    : ',tax%t0_db,tax%t0_req,tax%t0,tax%t_first

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
  double precision, dimension(NGLLX) :: hxi,hpxi
  double precision, dimension(NGLLY) :: heta,hpeta
  double precision, dimension(NGLLZ) :: hgam,hpgam
  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM) :: dw
  double precision, dimension(NDIM,NDIM) :: m_cart
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

  call gf_interp_weights_deriv(loc%xi,loc%eta,loc%gamma,hxi,hpxi,heta,hpeta,hgam,hpgam)

  if (do_strain) then
    call gf_strain_dweights(hxi,hpxi,heta,hpeta,hgam,hpgam,loc%jinv,dw)
    call gf_rotate_moment_tensor(loc%theta,loc%phi,src%moment_tensor,m_cart)
  endif

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

    call gf_interp_trace(displ,hxi,heta,hgam,nt,g)

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

      call gf_strain_trace(displ,dw,nt,eps)

      ! scaled exactly as gf_seis_cmt scales it, so the dumped trace is in
      ! the seismogram's units and differs from it only by the conversion
      do ia = 1,GF_NCOMP
        do it = 1,nt
          call gf_moment_contract(m_cart,eps(:,ia,it),pre_stf(ia,it))
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

  subroutine gf_write_header(db,src,loc,ista,iout)

! the '#' provenance block shared by both output formats
!
! Everything needed to reproduce the trace: which database, which source,
! which element and where in it, and the time axis. A trace whose provenance
! has to be reconstructed from the filename is a trace nobody can check.

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(in) :: src
  type(t_gf_location), intent(in) :: loc
  integer, intent(in) :: ista
  integer, intent(in) :: iout

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

  end subroutine gf_write_header

  end module gf_seismograms
