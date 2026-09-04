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
!---- Seismograms from a force source, and the ASCII output.
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
!---- Time axis
!---- ---------
!----   t(i) = (i*subsample_step - 1)*dt - t0,   i = 1..nt_subsampled
!----
!---- verified against the writer rather than assumed: green_function_io.F90
!---- :394 snapshots when mod(it,GF_SUBSAMPLE_STEP) == 0, so snapshot i is
!---- solver step i*subsample_step, whose simulation time is (it-1)*dt - t0.
!----
!---- Use dt_sub = dt*subsample_step for anything that steps along this axis:
!---- 0.4 s, not 0.1 s, in the shipped global example. GF3DF's GF%dt was
!---- already the stored spacing, so a literal port silently uses the solver
!---- step and is wrong by a factor of subsample_step.
!----
!---- The output-time-axis question -- a requested t0_out, and the left
!---- zero-padding it needs before the STF convolution -- belongs to Stage 5.
!---- For a force source with no STF correction the database's own axis is
!---- the right one.
!----

  module gf_seismograms

  use gf_par, only: t_gfdb,t_gf_location,t_gf_source,gf_set_error, &
                    GF_OK,GF_ERR_ARG,GF_ERR_ALLOC,GF_ERR_IO, &
                    GF_NCOMP,GF3D_VERSION

  use gf_database, only: gf_dir_exists

  use gf_element_io, only: gf_read_element_displ

  use gf_interp, only: gf_interp_weights,gf_interp_weights_deriv,gf_interp_trace

  use gf_strain, only: gf_strain_dweights,gf_strain_trace,GF_VOIGT

  use gf_moment, only: gf_rotate_moment_tensor,gf_moment_contract

  use gf_stf, only: gf_cumtrapz

  use gf_source, only: gf_force_direction

  implicit none

  private

  public :: gf_time_axis
  public :: gf_seis_force
  public :: gf_seis_cmt
  public :: gf_seis
  public :: gf_write_seis
  public :: gf_write_dump

  ! component labels, in the stored force order (green_function_io.F90:38)
  character(len=1), dimension(GF_NCOMP), parameter :: GF_COMP_NAME = (/ 'N','E','Z' /)

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_time_axis(db,nt,t)

! sample times of the stored traces, in seconds relative to the origin

  implicit none

  type(t_gfdb), intent(in) :: db
  integer, intent(in) :: nt
  double precision, dimension(nt), intent(out) :: t

  ! local parameters
  integer :: i

  do i = 1,nt
    t(i) = (dble(i*db%subsample_step) - 1.d0)*db%dt - db%t0
  enddo

  end subroutine gf_time_axis

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_seis_force(db,src,loc,seis,t,ierr)

! seismograms at every station for a force source
!
! `seis` is (nstations, 3 components N/E/Z, nt_subsampled) in metres.
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
  double precision, dimension(db%nstations,GF_NCOMP,db%nt_subsampled), intent(out) :: seis
  double precision, dimension(db%nt_subsampled), intent(out) :: t
  integer, intent(out) :: ierr

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: displ
  double precision, dimension(:,:,:), allocatable :: g
  double precision, dimension(NGLLX) :: hxi
  double precision, dimension(NGLLY) :: heta
  double precision, dimension(NGLLZ) :: hgam
  double precision, dimension(NDIM) :: fhat
  double precision :: scale_amp
  integer :: ista,it,icomp,idisp,ier,nt

  seis(:,:,:) = 0.d0

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

  nt = db%nt_subsampled

  call gf_time_axis(db,nt,t)

  call gf_interp_weights(loc%xi,loc%eta,loc%gamma,hxi,heta,hgam)

  call gf_force_direction(src,loc%nu,fhat,ierr)
  if (ierr /= GF_OK) return

  allocate(displ(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ,nt),g(GF_NCOMP,GF_NCOMP,nt),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate the element displacement buffer')
    return
  endif

  do ista = 1,db%nstations

    call gf_read_element_displ(db,loc%ielem,ista,displ,ierr)
    if (ierr /= GF_OK) goto 99

    ! g(a,d,t): the interpolated field at the source, still carrying both
    ! the station-force index a and the displacement index d
    call gf_interp_trace(displ,hxi,heta,hgam,nt,g)

    if (db%stations(ista)%factor_force_source == 0.d0) then
      call gf_set_error(ierr,GF_ERR_ARG, &
        'station '//trim(db%stations(ista)%id)//' has factor_force_source = 0')
      goto 99
    endif
    scale_amp = src%factor_force_source / db%stations(ista)%factor_force_source

    do it = 1,nt
      do icomp = 1,GF_NCOMP
        do idisp = 1,GF_NCOMP
          seis(ista,icomp,it) = seis(ista,icomp,it) + fhat(idisp)*g(icomp,idisp,it)
        enddo
        seis(ista,icomp,it) = scale_amp * seis(ista,icomp,it)
      enddo
    enddo

  enddo

  ierr = GF_OK

99 continue
  if (allocated(displ)) deallocate(displ)
  if (allocated(g)) deallocate(g)

  end subroutine gf_seis_force

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_seis_cmt(db,src,loc,seis,t,ierr)

! seismograms at every station for a moment-tensor source
!
! `seis` is (nstations, 3 components N/E/Z, nt_subsampled) in metres, as the
! *Heaviside* response at the database's own half duration. Correcting that
! half duration to the CMTSOLUTION's -- the hdur/1.628 Gaussian convolution
! -- is Stage 5's, and until it exists this trace cannot be compared to a
! forward run.
!
! Three steps, in this order:
!
!   1. the strain of the reciprocal field at the source, eps(6,3,nt);
!   2. contraction with the Cartesian moment tensor, giving the response to
!      a Gaussian source time function;
!   3. one cumulative trapezoidal integration, turning that into the
!      response to the Heaviside that a CMT source is.
!
! Step 3 is where a left-endpoint cumulative sum would lag by half a sample
! -- 1.7 s on the regional example's 3.4 s grid. See gf_stf.F90.

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(in) :: src
  type(t_gf_location), intent(in) :: loc
  double precision, dimension(db%nstations,GF_NCOMP,db%nt_subsampled), intent(out) :: seis
  double precision, dimension(db%nt_subsampled), intent(out) :: t
  integer, intent(out) :: ierr

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: displ
  double precision, dimension(:,:,:), allocatable :: eps
  double precision, dimension(:), allocatable :: trace,trace_int
  double precision, dimension(NGLLX) :: hxi,hpxi
  double precision, dimension(NGLLY) :: heta,hpeta
  double precision, dimension(NGLLZ) :: hgam,hpgam
  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM) :: dw
  double precision, dimension(NDIM,NDIM) :: m_cart
  double precision :: dt_sub,scale_amp
  integer :: ista,it,icomp,ier,nt

  seis(:,:,:) = 0.d0

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

  nt = db%nt_subsampled

  ! the spacing of the *stored* axis, not the solver step
  dt_sub = db%dt*dble(db%subsample_step)

  call gf_time_axis(db,nt,t)

  ! basis values and reference derivatives at the source, then the
  ! physical-space derivatives through the locator's inverse Jacobian
  call gf_interp_weights_deriv(loc%xi,loc%eta,loc%gamma,hxi,hpxi,heta,hpeta,hgam,hpgam)
  call gf_strain_dweights(hxi,hpxi,heta,hpeta,hgam,hpgam,loc%jinv,dw)

  ! the moment tensor, rotated once: it does not depend on the station
  call gf_rotate_moment_tensor(loc%theta,loc%phi,src%moment_tensor,m_cart)

  allocate(displ(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ,nt),eps(GF_VOIGT,GF_NCOMP,nt), &
           trace(nt),trace_int(nt),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate the element displacement buffer')
    return
  endif

  do ista = 1,db%nstations

    call gf_read_element_displ(db,loc%ielem,ista,displ,ierr)
    if (ierr /= GF_OK) goto 99

    call gf_strain_trace(displ,dw,nt,eps)

    ! per unit reciprocal force at the station; see the module header
    if (db%stations(ista)%factor_force_source == 0.d0) then
      call gf_set_error(ierr,GF_ERR_ARG, &
        'station '//trim(db%stations(ista)%id)//' has factor_force_source = 0')
      goto 99
    endif
    scale_amp = 1.d0 / db%stations(ista)%factor_force_source

    do icomp = 1,GF_NCOMP
      do it = 1,nt
        call gf_moment_contract(m_cart,eps(:,icomp,it),trace(it))
        trace(it) = scale_amp * trace(it)
      enddo

      ! Gaussian response -> Heaviside response
      call gf_cumtrapz(trace,nt,dt_sub,trace_int)

      do it = 1,nt
        seis(ista,icomp,it) = trace_int(it)
      enddo
    enddo

  enddo

  ierr = GF_OK

99 continue
  if (allocated(displ)) deallocate(displ)
  if (allocated(eps)) deallocate(eps)
  if (allocated(trace)) deallocate(trace)
  if (allocated(trace_int)) deallocate(trace_int)

  end subroutine gf_seis_cmt

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_seis(db,src,loc,seis,t,ierr)

! seismograms for either kind of source

  use gf_par, only: GF_SRC_FORCE,GF_SRC_CMT

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(in) :: src
  type(t_gf_location), intent(in) :: loc
  double precision, dimension(db%nstations,GF_NCOMP,db%nt_subsampled), intent(out) :: seis
  double precision, dimension(db%nt_subsampled), intent(out) :: t
  integer, intent(out) :: ierr

  select case (src%source_type)
  case (GF_SRC_FORCE)
    call gf_seis_force(db,src,loc,seis,t,ierr)
  case (GF_SRC_CMT)
    call gf_seis_cmt(db,src,loc,seis,t,ierr)
  case default
    call gf_set_error(ierr,GF_ERR_ARG,'gf_seis: the source has no type set')
  end select

  end subroutine gf_seis

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_write_seis(db,src,loc,seis,t,outdir,ierr)

! writes one ASCII file per station
!
! The format is a maintained interface, not a debugging convenience: Stages
! 4, 5 and 9 diff against it, and the Stage 5 comparison harness reads it
! with numpy.loadtxt. Keep the '#' header block and the four columns.

  use constants, only: MAX_STRING_LEN

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(in) :: src
  type(t_gf_location), intent(in) :: loc
  double precision, dimension(db%nstations,GF_NCOMP,db%nt_subsampled), intent(in) :: seis
  double precision, dimension(db%nt_subsampled), intent(in) :: t
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
    write(iout,'(a)') '# columns: t[s]  '//GF_COMP_NAME(1)//'[m]  ' &
                      //GF_COMP_NAME(2)//'[m]  '//GF_COMP_NAME(3)//'[m]'

    do it = 1,db%nt_subsampled
      write(iout,'(4es24.16)') t(it),seis(ista,1,it),seis(ista,2,it),seis(ista,3,it)
    enddo

    close(iout)

  enddo

  ierr = GF_OK

  end subroutine gf_write_seis

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
! *before* the time integration -- the two intermediates between the raw
! field and the seismogram, so a disagreement can be attributed to the
! geometric chain, the contraction, or the integration separately.
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
      ! the seismogram's units and differs from it only by the integration
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
