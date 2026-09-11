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
!---- The C ABI of libgf3d. Its contract is include/gf3d.h, which is the
!---- document to read first and to keep in step with this file.
!----
!---- This is the first bind(C) code in specfem3d_globe -- the tree's other
!---- C interoperability is all FC_FUNC_ in the other direction, C called
!---- from Fortran -- so the conventions are set here rather than inherited:
!----
!----   * every entry is a *function* returning integer(c_int), one of the
!----     GF_* codes, and writes its outputs through arguments;
!----   * the caller allocates every output array; nothing allocatable
!----     crosses the boundary;
!----   * strings come in null-terminated as character(kind=c_char),
!----     dimension(*) and go out as a (buffer, length) pair, truncated and
!----     always terminated;
!----   * a nullable output is type(c_ptr), value, tested with
!----     c_associated. Not `optional`: an optional bind(C) dummy is
!----     Fortran 2018, and the tests build with -std=f2008;
!----   * every real that arrives from outside is screened with
!----     gf_is_finite before anything arithmetic touches it.
!----
!---- Why an integer handle rather than an opaque pointer
!---- -------------------------------------------------
!---- The stage plan called for type(c_ptr) from c_loc on a saved target.
!---- An index into a table is used instead, because the whole reason this
!---- facade exists is that nothing reachable from a Python interpreter may
!---- crash it: a stale or forged pointer is undefined behaviour, whereas a
!---- stale or forged index is a bounds test and a GF_ERR_ARG. The cost is a
!---- fixed ceiling on simultaneously open databases (GF3D_MAX_HANDLES),
!---- which no caller is near.
!----
!---- Why values, not file paths
!---- --------------------------
!---- There is deliberately no "open this CMTSOLUTION" entry point. Reading
!---- one means get_cmt()/get_force(), which stop on malformed input; the
!---- caller passes numbers instead and gf_source_set_cmt/force applies the
!---- readers' own arithmetic to them. Parsing the file is the Python side's
!---- job. See the header of gf_source.F90.
!----
!---- Process-wide state, and what it means for a long-lived caller
!---- ------------------------------------------------------------
!---- gf_errmsg, the kd-tree in src/shared/search_kdtree.f90, and specfem's
!---- shared_parameters are all one per process, not one per handle. So:
!----     - the facade is not thread-safe, and says so in the header;
!----     - every computing entry re-installs shared_parameters from its own
!----       handle first, so that opening a second database cannot leave the
!----       first one indexing a topography grid with the wrong dimensions.
!----

  module gf3d_capi

  use, intrinsic :: iso_c_binding, only: &
    c_int, c_double, c_char, c_ptr, c_null_char, c_associated, c_f_pointer, c_sizeof

  use constants, only: MAX_STRING_LEN

  ! GF3D_VERSION is renamed on import: Fortran is case-insensitive, so the
  ! constant and the gf3d_version() entry point below would otherwise be the
  ! same name. Renaming the constant rather than the function keeps every
  ! procedure here spelled exactly like the C symbol it binds.
  use gf_par, only: t_gfdb, t_gf_source, t_gf_location, t_gf_taxis, t_gf_stf, &
                    GF_VERSION_STRING => GF3D_VERSION, GF_NCOMP, &
                    GF_OK, GF_ERR_ARG, GF_ERR_ALLOC, &
                    GF_SRC_FORCE, GF_SRC_CMT, &
                    gf_set_error, gf_error_string, gf_errmsg, gf_is_finite

  use gf_shared_params, only: gf_init_shared_params

  use gf_database, only: gf_open, gf_close

  use gf_source, only: gf_source_set_cmt, gf_source_set_force

  use gf_locate, only: gf_locate_source, gf_locate_release

  use gf_seismograms, only: gf_seis_plan, gf_seis, gf_seis_cmt_partials

  use gf_partials, only: gf_partials_ndp, GF_NDP_LOC, GF_DP_NAME, GF_DP_UNIT

  use shared_parameters, only: R_PLANET, RHOAV

  implicit none

  private

  !--- must match include/gf3d.h ---
  integer, parameter :: GF3D_MAX_HANDLES = 32
  integer, parameter :: GF3D_STRLEN = 64
  integer, parameter :: GF3D_MORTON_STRLEN = 24

  !-----------------------------------------------------------------
  ! the interoperable mirrors of the library's derived types
  !
  ! Field for field, and in the same order as the structs in gf3d.h. Note
  ! the integers are grouped ahead of the doubles, with an explicit pad
  ! where the count is odd, so that the layout a C compiler chooses and the
  ! layout gfortran chooses cannot differ by an alignment hole. gf3d_sizeof
  ! below lets a binding check that at run time rather than trust it.
  !
  ! No default initialisation: these are plain data.
  !-----------------------------------------------------------------

  type, bind(C) :: gf3d_source_t
    integer(c_int) :: source_type
    integer(c_int) :: force_stf
    real(c_double) :: latitude
    real(c_double) :: longitude
    real(c_double) :: depth_km
    real(c_double) :: hdur
    real(c_double) :: time_shift
    real(c_double), dimension(6) :: moment
    real(c_double) :: force_factor
    real(c_double), dimension(3) :: force_dir
  end type gf3d_source_t

  type, bind(C) :: gf3d_info_t
    integer(c_int) :: nelem
    integer(c_int) :: nstations
    integer(c_int) :: nstep
    integer(c_int) :: nt_subsampled
    integer(c_int) :: subsample_step
    integer(c_int) :: ngllx, nglly, ngllz
    integer(c_int) :: topography
    integer(c_int) :: ellipticity
    integer(c_int) :: rotation
    integer(c_int) :: attenuation
    integer(c_int) :: gravity
    integer(c_int) :: pad_
    real(c_double) :: dt
    real(c_double) :: t0
    real(c_double) :: r_planet
    real(c_double) :: rhoav
    real(c_double) :: scale_displ
  end type gf3d_info_t

  type, bind(C) :: gf3d_station_t
    character(kind=c_char), dimension(GF3D_STRLEN) :: id
    character(kind=c_char), dimension(GF3D_STRLEN) :: network
    character(kind=c_char), dimension(GF3D_STRLEN) :: station
    real(c_double) :: latitude
    real(c_double) :: longitude
    real(c_double) :: depth_m
    real(c_double) :: hdur
    real(c_double) :: f_cutoff
    real(c_double) :: factor_force_source
    real(c_double) :: time_shift
  end type gf3d_station_t

  type, bind(C) :: gf3d_location_t
    integer(c_int) :: ielem
    character(kind=c_char), dimension(GF3D_MORTON_STRLEN) :: morton_hex
    real(c_double) :: xi, eta, gamma
    real(c_double), dimension(3) :: xyz
    real(c_double), dimension(3) :: xyz_target
    real(c_double) :: distance_km
    real(c_double) :: anchor_err
    real(c_double) :: theta
    real(c_double) :: phi
    real(c_double) :: r_surface
  end type gf3d_location_t

  type, bind(C) :: gf3d_plan_t
    integer(c_int) :: nt
    integer(c_int) :: nt_db
    integer(c_int) :: npad
    integer(c_int) :: subsample_step
    integer(c_int) :: khalf
    integer(c_int) :: guard
    integer(c_int) :: kind_stf
    integer(c_int) :: pad_
    real(c_double) :: dt
    real(c_double) :: dt_sub
    real(c_double) :: t0_db
    real(c_double) :: t0_req
    real(c_double) :: t0
    real(c_double) :: t_first
    real(c_double) :: hdur_src
    real(c_double) :: hdur_target
    real(c_double) :: hdur_db
    real(c_double) :: hdur_corr
    real(c_double) :: trunc
  end type gf3d_plan_t

  !-----------------------------------------------------------------
  ! the handle table
  !-----------------------------------------------------------------

  type(t_gfdb), dimension(GF3D_MAX_HANDLES), save :: handles
  logical, dimension(GF3D_MAX_HANDLES), save :: in_use = .false.

  public :: gf3d_version, gf3d_sizeof, gf3d_last_error, gf3d_error_string
  public :: gf3d_open, gf3d_close, gf3d_get_info, gf3d_get_station
  public :: gf3d_locate, gf3d_get_plan, gf3d_ndp, gf3d_partial_name
  public :: gf3d_seismograms, gf3d_partials

  contains

!
!===================================================================
! helpers (not exported)
!===================================================================
!

  subroutine put_string(fstr,buf,buflen)

! copies a Fortran string into a C buffer, truncating, always terminating

  implicit none

  character(len=*), intent(in) :: fstr
  character(kind=c_char), dimension(*), intent(out) :: buf
  integer, intent(in) :: buflen

  ! local parameters
  integer :: i,n

  if (buflen < 1) return

  n = min(len_trim(fstr),buflen-1)
  do i = 1,n
    buf(i) = fstr(i:i)
  enddo
  buf(n+1) = c_null_char

  end subroutine put_string

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_string(buf,fstr)

! copies a null-terminated C string into a Fortran string
!
! A string longer than fstr is truncated rather than refused; every caller
! below then fails on the truncated path, which is a clearer error than a
! length complaint.

  implicit none

  character(kind=c_char), dimension(*), intent(in) :: buf
  character(len=*), intent(out) :: fstr

  ! local parameters
  integer :: i

  fstr = ''
  do i = 1,len(fstr)
    if (buf(i) == c_null_char) return
    fstr(i:i) = buf(i)
  enddo

  end subroutine get_string

!
!-------------------------------------------------------------------------------------------------
!

  subroutine use_handle(h,ierr)

! validates a handle and re-installs specfem's globals from it
!
! The re-installation is the part that is easy to miss: R_PLANET, NX_BATHY,
! NY_BATHY and RESOLUTION_TOPO_FILE live in shared_parameters, one set per
! process, and gf_open() writes them. With two databases open, the second
! open would otherwise leave the first handle's topography grid being
! indexed with the second's dimensions -- and get_topo_bathy() does not
! range-check, it returns plausible-looking garbage.

  implicit none

  integer(c_int), intent(in) :: h
  integer, intent(out) :: ierr

  if (h < 1 .or. h > GF3D_MAX_HANDLES) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf3d: handle out of range')
    return
  endif

  if (.not. in_use(h)) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf3d: handle is not open')
    return
  endif

  call gf_init_shared_params(handles(h),ierr)

  end subroutine use_handle

!
!-------------------------------------------------------------------------------------------------
!

  subroutine build_source(csrc,db,src,ierr)

! turns the C struct into the library's own source type

  implicit none

  type(gf3d_source_t), intent(in) :: csrc
  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(out) :: src
  integer, intent(out) :: ierr

  select case (csrc%source_type)
  case (GF_SRC_CMT)
    call gf_source_set_cmt(src,csrc%latitude,csrc%longitude,csrc%depth_km, &
                           csrc%hdur,csrc%time_shift,csrc%moment,db%dt,ierr)
  case (GF_SRC_FORCE)
    call gf_source_set_force(src,csrc%latitude,csrc%longitude,csrc%depth_km, &
                             csrc%hdur,csrc%time_shift,int(csrc%force_stf), &
                             csrc%force_factor, &
                             csrc%force_dir(1),csrc%force_dir(2),csrc%force_dir(3), &
                             db%dt,ierr)
  case default
    call gf_set_error(ierr,GF_ERR_ARG, &
      'gf3d: source_type must be GF_SRC_CMT (2) or GF_SRC_FORCE (1)')
  end select

  end subroutine build_source

!
!-------------------------------------------------------------------------------------------------
!

  subroutine fill_location(loc,cloc)

  implicit none

  type(t_gf_location), intent(in) :: loc
  type(gf3d_location_t), intent(out) :: cloc

  cloc%ielem = int(loc%ielem,kind=c_int)
  call put_string(loc%morton_hex,cloc%morton_hex,GF3D_MORTON_STRLEN)
  cloc%xi = loc%xi
  cloc%eta = loc%eta
  cloc%gamma = loc%gamma
  cloc%xyz(1:3) = loc%xyz(1:3)
  cloc%xyz_target(1:3) = loc%xyz_target(1:3)
  cloc%distance_km = loc%distance_km
  cloc%anchor_err = loc%anchor_err
  cloc%theta = loc%theta
  cloc%phi = loc%phi
  cloc%r_surface = loc%r_surface

  end subroutine fill_location

!
!-------------------------------------------------------------------------------------------------
!

  subroutine fill_plan(tax,stf,t0_req,cplan)

  implicit none

  type(t_gf_taxis), intent(in) :: tax
  type(t_gf_stf), intent(in) :: stf
  double precision, intent(in) :: t0_req
  type(gf3d_plan_t), intent(out) :: cplan

  cplan%nt = int(tax%nt,kind=c_int)
  cplan%nt_db = int(tax%nt_db,kind=c_int)
  cplan%npad = int(tax%npad,kind=c_int)
  cplan%subsample_step = int(tax%subsample_step,kind=c_int)
  cplan%khalf = int(stf%khalf,kind=c_int)
  cplan%guard = 0
  if (stf%guard) cplan%guard = 1
  cplan%kind_stf = int(stf%kind_stf,kind=c_int)
  cplan%pad_ = 0
  cplan%dt = tax%dt
  cplan%dt_sub = tax%dt_sub
  cplan%t0_db = tax%t0_db
  cplan%t0_req = t0_req
  cplan%t0 = tax%t0
  cplan%t_first = tax%t_first
  cplan%hdur_src = stf%hdur_src
  cplan%hdur_target = stf%hdur_target
  cplan%hdur_db = stf%hdur_db
  cplan%hdur_corr = stf%hdur_corr
  cplan%trunc = stf%trunc

  end subroutine fill_plan

!
!-------------------------------------------------------------------------------------------------
!

  subroutine pack_seis(nsta,nt,seis,seis_c)

! (nsta,3,nt) Fortran-ordered -> [ista][icomp][it] C-ordered
!
! One copy, made here rather than left to the caller, so that a C or numpy
! user indexes seis[ista][icomp][it] with no transposition of their own.

  implicit none

  integer, intent(in) :: nsta,nt
  double precision, dimension(nsta,GF_NCOMP,nt), intent(in) :: seis
  real(c_double), dimension(nt,GF_NCOMP,nsta), intent(out) :: seis_c

  ! local parameters
  integer :: ista,icomp,it

  do ista = 1,nsta
    do icomp = 1,GF_NCOMP
      do it = 1,nt
        seis_c(it,icomp,ista) = seis(ista,icomp,it)
      enddo
    enddo
  enddo

  end subroutine pack_seis

!
!-------------------------------------------------------------------------------------------------
!

  subroutine pack_dp(ndp,nsta,nt,dp,dp_c)

! (ndp,nsta,3,nt) Fortran-ordered -> [ista][ip][icomp][it] C-ordered

  implicit none

  integer, intent(in) :: ndp,nsta,nt
  double precision, dimension(ndp,nsta,GF_NCOMP,nt), intent(in) :: dp
  real(c_double), dimension(nt,GF_NCOMP,ndp,nsta), intent(out) :: dp_c

  ! local parameters
  integer :: ista,icomp,it,ip

  do ista = 1,nsta
    do ip = 1,ndp
      do icomp = 1,GF_NCOMP
        do it = 1,nt
          dp_c(it,icomp,ip,ista) = dp(ip,ista,icomp,it)
        enddo
      enddo
    enddo
  enddo

  end subroutine pack_dp

!
!===================================================================
! version, sizes and errors
!===================================================================
!

  integer(c_int) function gf3d_version(buf,buflen) bind(C,name='gf3d_version')

  implicit none

  character(kind=c_char), dimension(*), intent(out) :: buf
  integer(c_int), value :: buflen

  call put_string(GF_VERSION_STRING,buf,int(buflen))

  gf3d_version = GF_OK

  end function gf3d_version

!
!-------------------------------------------------------------------------------------------------
!

  integer(c_int) function gf3d_sizeof(source,info,station,location,plan) bind(C,name='gf3d_sizeof')

! the size of each interoperable struct, so that a binding written against
! gf3d.h can check at load time that it agrees with the library it found
!
! Any argument may be NULL.

  implicit none

  type(c_ptr), value :: source,info,station,location,plan

  ! local parameters
  integer(c_int), pointer :: p
  type(gf3d_source_t) :: a
  type(gf3d_info_t) :: b
  type(gf3d_station_t) :: c
  type(gf3d_location_t) :: d
  type(gf3d_plan_t) :: e

  if (c_associated(source)) then
    call c_f_pointer(source,p) ; p = int(c_sizeof(a),kind=c_int)
  endif
  if (c_associated(info)) then
    call c_f_pointer(info,p) ; p = int(c_sizeof(b),kind=c_int)
  endif
  if (c_associated(station)) then
    call c_f_pointer(station,p) ; p = int(c_sizeof(c),kind=c_int)
  endif
  if (c_associated(location)) then
    call c_f_pointer(location,p) ; p = int(c_sizeof(d),kind=c_int)
  endif
  if (c_associated(plan)) then
    call c_f_pointer(plan,p) ; p = int(c_sizeof(e),kind=c_int)
  endif

  gf3d_sizeof = GF_OK

  end function gf3d_sizeof

!
!-------------------------------------------------------------------------------------------------
!

  integer(c_int) function gf3d_last_error(buf,buflen) bind(C,name='gf3d_last_error')

  implicit none

  character(kind=c_char), dimension(*), intent(out) :: buf
  integer(c_int), value :: buflen

  call put_string(gf_errmsg,buf,int(buflen))

  gf3d_last_error = GF_OK

  end function gf3d_last_error

!
!-------------------------------------------------------------------------------------------------
!

  integer(c_int) function gf3d_error_string(code,buf,buflen) bind(C,name='gf3d_error_string')

  implicit none

  integer(c_int), value :: code
  character(kind=c_char), dimension(*), intent(out) :: buf
  integer(c_int), value :: buflen

  call put_string(gf_error_string(int(code)),buf,int(buflen))

  gf3d_error_string = GF_OK

  end function gf3d_error_string

!
!===================================================================
! opening and interrogating a database
!===================================================================
!

  integer(c_int) function gf3d_open(path,check_completion,h) bind(C,name='gf3d_open')

  implicit none

  character(kind=c_char), dimension(*), intent(in) :: path
  integer(c_int), value :: check_completion
  integer(c_int), intent(out) :: h

  ! local parameters
  character(len=MAX_STRING_LEN) :: fpath
  integer :: islot,i,ierr

  h = 0

  call get_string(path,fpath)

  if (len_trim(fpath) == 0) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf3d_open: no database path given')
    gf3d_open = int(ierr,kind=c_int)
    return
  endif

  islot = 0
  do i = 1,GF3D_MAX_HANDLES
    if (.not. in_use(i)) then
      islot = i
      exit
    endif
  enddo

  if (islot == 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'gf3d_open: too many databases open at once')
    gf3d_open = int(ierr,kind=c_int)
    return
  endif

  call gf_open(trim(fpath),handles(islot),ierr,check_completion = (check_completion /= 0))

  if (ierr /= GF_OK) then
    ! gf_open closes what it opened; the slot stays free
    gf3d_open = int(ierr,kind=c_int)
    return
  endif

  in_use(islot) = .true.
  h = int(islot,kind=c_int)

  gf3d_open = GF_OK

  end function gf3d_open

!
!-------------------------------------------------------------------------------------------------
!

  integer(c_int) function gf3d_close(h) bind(C,name='gf3d_close')

! closes a database and releases the process-wide search tree
!
! The tree is released unconditionally, even when another handle is still
! open: it belongs to whichever database located last, and the next locate
! rebuilds it for whoever asks. Leaking it instead would leave
! src/shared/search_kdtree.f90's module arrays allocated for the life of a
! Python interpreter.

  implicit none

  integer(c_int), value :: h

  ! local parameters
  integer :: ierr

  if (h < 1 .or. h > GF3D_MAX_HANDLES) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf3d_close: handle out of range')
    gf3d_close = int(ierr,kind=c_int)
    return
  endif

  if (.not. in_use(h)) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf3d_close: handle is not open')
    gf3d_close = int(ierr,kind=c_int)
    return
  endif

  call gf_locate_release()
  call gf_close(handles(h))

  in_use(h) = .false.

  gf3d_close = GF_OK

  end function gf3d_close

!
!-------------------------------------------------------------------------------------------------
!

  integer(c_int) function gf3d_get_info(h,info) bind(C,name='gf3d_get_info')

  implicit none

  integer(c_int), value :: h
  type(gf3d_info_t), intent(out) :: info

  ! local parameters
  integer :: ierr

  call use_handle(h,ierr)
  if (ierr /= GF_OK) then
    gf3d_get_info = int(ierr,kind=c_int)
    return
  endif

  info%nelem          = int(handles(h)%nelem,kind=c_int)
  info%nstations      = int(handles(h)%nstations,kind=c_int)
  info%nstep          = int(handles(h)%nstep,kind=c_int)
  info%nt_subsampled  = int(handles(h)%nt_subsampled,kind=c_int)
  info%subsample_step = int(handles(h)%subsample_step,kind=c_int)
  info%ngllx          = int(handles(h)%ngllx,kind=c_int)
  info%nglly          = int(handles(h)%nglly,kind=c_int)
  info%ngllz          = int(handles(h)%ngllz,kind=c_int)

  info%topography  = 0 ; if (handles(h)%topography)  info%topography  = 1
  info%ellipticity = 0 ; if (handles(h)%ellipticity) info%ellipticity = 1
  info%rotation    = 0 ; if (handles(h)%rotation)    info%rotation    = 1
  info%attenuation = 0 ; if (handles(h)%attenuation) info%attenuation = 1
  info%gravity     = 0 ; if (handles(h)%gravity)     info%gravity     = 1
  info%pad_        = 0

  info%dt          = handles(h)%dt
  info%t0          = handles(h)%t0
  info%scale_displ = handles(h)%scale_displ

  ! the values the library is actually using, not the (possibly zero) ones
  ! the database recorded: the current writer stores neither RHOAV nor the
  ! flattening, so these fall back to specfem's Earth defaults
  info%r_planet = R_PLANET
  info%rhoav    = RHOAV

  gf3d_get_info = GF_OK

  end function gf3d_get_info

!
!-------------------------------------------------------------------------------------------------
!

  integer(c_int) function gf3d_get_station(h,ista,sta) bind(C,name='gf3d_get_station')

! ista is 0-based, as C expects; the Fortran array below is 1-based

  implicit none

  integer(c_int), value :: h
  integer(c_int), value :: ista
  type(gf3d_station_t), intent(out) :: sta

  ! local parameters
  integer :: ierr,i

  call use_handle(h,ierr)
  if (ierr /= GF_OK) then
    gf3d_get_station = int(ierr,kind=c_int)
    return
  endif

  i = int(ista) + 1

  if (i < 1 .or. i > handles(h)%nstations) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf3d_get_station: station index out of range')
    gf3d_get_station = int(ierr,kind=c_int)
    return
  endif

  call put_string(handles(h)%stations(i)%id,sta%id,GF3D_STRLEN)
  call put_string(handles(h)%stations(i)%network,sta%network,GF3D_STRLEN)
  call put_string(handles(h)%stations(i)%station,sta%station,GF3D_STRLEN)

  sta%latitude            = handles(h)%stations(i)%latitude
  sta%longitude           = handles(h)%stations(i)%longitude
  sta%depth_m             = handles(h)%stations(i)%depth
  sta%hdur                = handles(h)%stations(i)%hdur
  sta%f_cutoff            = handles(h)%stations(i)%f_cutoff
  sta%factor_force_source = handles(h)%stations(i)%factor_force_source
  sta%time_shift          = handles(h)%stations(i)%time_shift

  gf3d_get_station = GF_OK

  end function gf3d_get_station

!
!===================================================================
! locating and planning
!===================================================================
!

  integer(c_int) function gf3d_locate(h,lat,lon,depth_km,loc) bind(C,name='gf3d_locate')

  implicit none

  integer(c_int), value :: h
  real(c_double), value :: lat,lon,depth_km
  type(gf3d_location_t), intent(out) :: loc

  ! local parameters
  type(t_gf_location) :: floc
  integer :: ierr

  call use_handle(h,ierr)
  if (ierr /= GF_OK) then
    gf3d_locate = int(ierr,kind=c_int)
    return
  endif

  ! a NaN would pass reduce()'s range test -- both comparisons are false --
  ! and reach the kd-tree, which stops when it finds no point
  if (.not. (gf_is_finite(lat) .and. gf_is_finite(lon) .and. gf_is_finite(depth_km))) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf3d_locate: latitude, longitude or depth is not finite')
    gf3d_locate = int(ierr,kind=c_int)
    return
  endif

  call gf_locate_source(handles(h),lat,lon,depth_km,floc,ierr)
  if (ierr /= GF_OK) then
    gf3d_locate = int(ierr,kind=c_int)
    return
  endif

  call fill_location(floc,loc)

  gf3d_locate = GF_OK

  end function gf3d_locate

!
!-------------------------------------------------------------------------------------------------
!

  integer(c_int) function gf3d_get_plan(h,src,t0_req,plan) bind(C,name='gf3d_get_plan')

  implicit none

  integer(c_int), value :: h
  type(gf3d_source_t), intent(in) :: src
  real(c_double), value :: t0_req
  type(gf3d_plan_t), intent(out) :: plan

  ! local parameters
  type(t_gf_source) :: fsrc
  type(t_gf_taxis) :: tax
  type(t_gf_stf) :: stf
  integer :: ierr

  call use_handle(h,ierr)
  if (ierr /= GF_OK) then
    gf3d_get_plan = int(ierr,kind=c_int)
    return
  endif

  if (.not. gf_is_finite(t0_req)) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf3d_get_plan: t0 is not finite')
    gf3d_get_plan = int(ierr,kind=c_int)
    return
  endif

  call build_source(src,handles(h),fsrc,ierr)
  if (ierr /= GF_OK) then
    gf3d_get_plan = int(ierr,kind=c_int)
    return
  endif

  call gf_seis_plan(handles(h),fsrc,t0_req,tax,stf,ierr)
  if (ierr /= GF_OK) then
    gf3d_get_plan = int(ierr,kind=c_int)
    return
  endif

  call fill_plan(tax,stf,t0_req,plan)

  gf3d_get_plan = GF_OK

  end function gf3d_get_plan

!
!-------------------------------------------------------------------------------------------------
!

  integer(c_int) function gf3d_ndp(itypsokern,ndp) bind(C,name='gf3d_ndp')

  implicit none

  integer(c_int), value :: itypsokern
  integer(c_int), intent(out) :: ndp

  ! local parameters
  integer :: n,ierr

  call gf_partials_ndp(int(itypsokern),n,ierr)

  ndp = int(n,kind=c_int)

  gf3d_ndp = int(ierr,kind=c_int)

  end function gf3d_ndp

!
!-------------------------------------------------------------------------------------------------
!

  integer(c_int) function gf3d_partial_name(ip,name,namelen,unit,unitlen) &
    bind(C,name='gf3d_partial_name')

! ip is 0-based; either buffer may be NULL

  implicit none

  integer(c_int), value :: ip
  type(c_ptr), value :: name
  integer(c_int), value :: namelen
  type(c_ptr), value :: unit
  integer(c_int), value :: unitlen

  ! local parameters
  character(kind=c_char), dimension(:), pointer :: buf
  integer :: i,ierr

  i = int(ip) + 1

  if (i < 1 .or. i > GF_NDP_LOC) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf3d_partial_name: partial index out of range')
    gf3d_partial_name = int(ierr,kind=c_int)
    return
  endif

  if (c_associated(name) .and. namelen > 0) then
    call c_f_pointer(name,buf,(/ int(namelen) /))
    call put_string(GF_DP_NAME(i),buf,int(namelen))
  endif

  if (c_associated(unit) .and. unitlen > 0) then
    call c_f_pointer(unit,buf,(/ int(unitlen) /))
    call put_string(GF_DP_UNIT(i),buf,int(unitlen))
  endif

  gf3d_partial_name = GF_OK

  end function gf3d_partial_name

!
!===================================================================
! extraction
!===================================================================
!

  integer(c_int) function gf3d_seismograms(h,src,t0_req,nt,seis,t,onset,loc) &
    bind(C,name='gf3d_seismograms')

! seismograms at every station, in C order [ista][icomp][it]
!
! Mirrors gf3d_main.F90's --seis path exactly: locate (which also loads the
! topography grid, lazily), plan, extract. `loc` may be NULL.

  implicit none

  integer(c_int), value :: h
  type(gf3d_source_t), intent(in) :: src
  real(c_double), value :: t0_req
  integer(c_int), value :: nt
  real(c_double), dimension(*), intent(out) :: seis
  real(c_double), dimension(*), intent(out) :: t
  real(c_double), dimension(*), intent(out) :: onset
  type(c_ptr), value :: loc

  ! local parameters
  type(t_gf_source) :: fsrc
  type(t_gf_location) :: floc
  type(gf3d_location_t), pointer :: ploc
  type(t_gf_taxis) :: tax
  type(t_gf_stf) :: stf
  double precision, dimension(:,:,:), allocatable :: fseis
  double precision, dimension(:), allocatable :: ft,fonset
  integer :: ierr,ier,nsta,it,ista

  call use_handle(h,ierr)
  if (ierr /= GF_OK) then
    gf3d_seismograms = int(ierr,kind=c_int)
    return
  endif

  if (.not. gf_is_finite(t0_req)) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf3d_seismograms: t0 is not finite')
    gf3d_seismograms = int(ierr,kind=c_int)
    return
  endif

  call build_source(src,handles(h),fsrc,ierr)
  if (ierr /= GF_OK) then
    gf3d_seismograms = int(ierr,kind=c_int)
    return
  endif

  call gf_locate_source(handles(h),fsrc%latitude,fsrc%longitude,fsrc%depth,floc,ierr)
  if (ierr /= GF_OK) then
    gf3d_seismograms = int(ierr,kind=c_int)
    return
  endif

  call gf_seis_plan(handles(h),fsrc,t0_req,tax,stf,ierr)
  if (ierr /= GF_OK) then
    gf3d_seismograms = int(ierr,kind=c_int)
    return
  endif

  if (int(nt) /= tax%nt) then
    call gf_set_error(ierr,GF_ERR_ARG, &
      'gf3d_seismograms: nt does not match the plan; call gf3d_get_plan first')
    gf3d_seismograms = int(ierr,kind=c_int)
    return
  endif

  nsta = handles(h)%nstations

  allocate(fseis(nsta,GF_NCOMP,tax%nt),ft(tax%nt),fonset(nsta),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'gf3d_seismograms: could not allocate the work arrays')
    gf3d_seismograms = int(ierr,kind=c_int)
    return
  endif

  call gf_seis(handles(h),fsrc,floc,tax,stf,fseis,ft,fonset,ierr)

  if (ierr == GF_OK) then
    call pack_seis(nsta,tax%nt,fseis,seis)
    do it = 1,tax%nt
      t(it) = ft(it)
    enddo
    do ista = 1,nsta
      onset(ista) = fonset(ista)
    enddo
    if (c_associated(loc)) then
      call c_f_pointer(loc,ploc)
      call fill_location(floc,ploc)
    endif
  endif

  deallocate(fseis,ft,fonset)

  gf3d_seismograms = int(ierr,kind=c_int)

  end function gf3d_seismograms

!
!-------------------------------------------------------------------------------------------------
!

  integer(c_int) function gf3d_partials(h,src,t0_req,itypsokern,nt,ndp, &
                                        seis,dp,t,onset,loc) bind(C,name='gf3d_partials')

! seismograms and their partial derivatives, in C order
! [ista][icomp][it] and [ista][ip][icomp][it]

  implicit none

  integer(c_int), value :: h
  type(gf3d_source_t), intent(in) :: src
  real(c_double), value :: t0_req
  integer(c_int), value :: itypsokern
  integer(c_int), value :: nt
  integer(c_int), value :: ndp
  real(c_double), dimension(*), intent(out) :: seis
  real(c_double), dimension(*), intent(out) :: dp
  real(c_double), dimension(*), intent(out) :: t
  real(c_double), dimension(*), intent(out) :: onset
  type(c_ptr), value :: loc

  ! local parameters
  type(t_gf_source) :: fsrc
  type(t_gf_location) :: floc
  type(gf3d_location_t), pointer :: ploc
  type(t_gf_taxis) :: tax
  type(t_gf_stf) :: stf
  double precision, dimension(:,:,:), allocatable :: fseis
  double precision, dimension(:,:,:,:), allocatable :: fdp
  double precision, dimension(:), allocatable :: ft,fonset
  integer :: ierr,ier,nsta,ndp_want,it,ista

  call use_handle(h,ierr)
  if (ierr /= GF_OK) then
    gf3d_partials = int(ierr,kind=c_int)
    return
  endif

  if (.not. gf_is_finite(t0_req)) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf3d_partials: t0 is not finite')
    gf3d_partials = int(ierr,kind=c_int)
    return
  endif

  if (int(itypsokern) < 1) then
    call gf_set_error(ierr,GF_ERR_ARG, &
      'gf3d_partials: itypsokern must be 1 or 2; use gf3d_seismograms for 0')
    gf3d_partials = int(ierr,kind=c_int)
    return
  endif

  call gf_partials_ndp(int(itypsokern),ndp_want,ierr)
  if (ierr /= GF_OK) then
    gf3d_partials = int(ierr,kind=c_int)
    return
  endif

  if (int(ndp) /= ndp_want) then
    call gf_set_error(ierr,GF_ERR_ARG, &
      'gf3d_partials: ndp does not match itypsokern; ask gf3d_ndp for it')
    gf3d_partials = int(ierr,kind=c_int)
    return
  endif

  call build_source(src,handles(h),fsrc,ierr)
  if (ierr /= GF_OK) then
    gf3d_partials = int(ierr,kind=c_int)
    return
  endif

  if (fsrc%source_type /= GF_SRC_CMT) then
    call gf_set_error(ierr,GF_ERR_ARG, &
      'gf3d_partials: partial derivatives are defined for a moment-tensor source only')
    gf3d_partials = int(ierr,kind=c_int)
    return
  endif

  call gf_locate_source(handles(h),fsrc%latitude,fsrc%longitude,fsrc%depth,floc,ierr)
  if (ierr /= GF_OK) then
    gf3d_partials = int(ierr,kind=c_int)
    return
  endif

  call gf_seis_plan(handles(h),fsrc,t0_req,tax,stf,ierr)
  if (ierr /= GF_OK) then
    gf3d_partials = int(ierr,kind=c_int)
    return
  endif

  if (int(nt) /= tax%nt) then
    call gf_set_error(ierr,GF_ERR_ARG, &
      'gf3d_partials: nt does not match the plan; call gf3d_get_plan first')
    gf3d_partials = int(ierr,kind=c_int)
    return
  endif

  nsta = handles(h)%nstations

  allocate(fseis(nsta,GF_NCOMP,tax%nt),fdp(ndp_want,nsta,GF_NCOMP,tax%nt), &
           ft(tax%nt),fonset(nsta),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'gf3d_partials: could not allocate the work arrays')
    gf3d_partials = int(ierr,kind=c_int)
    return
  endif

  call gf_seis_cmt_partials(handles(h),fsrc,floc,tax,stf,int(itypsokern),ndp_want, &
                            fseis,fdp,ft,fonset,ierr)

  if (ierr == GF_OK) then
    call pack_seis(nsta,tax%nt,fseis,seis)
    call pack_dp(ndp_want,nsta,tax%nt,fdp,dp)
    do it = 1,tax%nt
      t(it) = ft(it)
    enddo
    do ista = 1,nsta
      onset(ista) = fonset(ista)
    enddo
    if (c_associated(loc)) then
      call c_f_pointer(loc,ploc)
      call fill_location(floc,ploc)
    endif
  endif

  deallocate(fseis,fdp,ft,fonset)

  gf3d_partials = int(ierr,kind=c_int)

  end function gf3d_partials

  end module gf3d_capi
