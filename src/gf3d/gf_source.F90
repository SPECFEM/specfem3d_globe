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
!---- Reading a seismic source.
!----
!---- A thin wrapper over src/specfem3D/get_force.f90 (and, from Stage 4,
!---- get_cmt.f90) rather than a second parser. That is the whole point of
!---- going in-repo: those readers also apply the non-dimensionalisation --
!----
!----     scaleF = RHOAV * R_PLANET**4 * PI*GRAV*RHOAV      (get_force.f90:289)
!----     scaleM = 1.d7 * RHOAV * R_PLANET**5 * PI*GRAV*RHOAV  (get_cmt.f90:426)
!----
!---- so the amplitude scaling comes from the solver's own source and cannot
!---- drift away from it. A pure scale error is invisible in waveform shape
!---- and fatal for a moment inversion, which is exactly the kind of mistake
!---- a re-implementation would be free to make.
!----
!---- PR 1 gave both readers a required `filename` argument for this; pass
!---- '' to get the solver's usual DATA/FORCESOLUTION behaviour.
!----
!---- What lives outside the readers, and is a silent bug if missed
!---- ------------------------------------------------------------
!---- * get_force stores the dominant frequency f0 in `hdur`, not a half
!----   duration. get_cmt stores the raw *triangle* half duration; the
!----   conversion to a Gaussian width is hdur/SOURCE_DECAY_MIMIC_TRIANGLE
!----   and happens at setup_sources_receivers.f90:777, in the caller.
!---- * Both zero `tshift_src` when NSOURCES == 1 and return the original in
!----   min_tshift_src_original (29.0 s for the shipped example), so t = 0 is
!----   the centroid time. Those seconds are origin-time metadata for Stage
!----   10's SAC headers, not part of the trace.
!----
!---- The `stop`s in the readers, and the way round them (Stage 9)
!---- -----------------------------------------------------------
!---- get_force `stop`s on malformed input (get_force.f90:240,248,272) and
!---- get_cmt on a missing file (:103), which breaks this library's no-stop
!---- invariant. A `stop` inside lib/libgf3d.so kills the calling Python
!---- interpreter with no traceback, so the C facade must not be able to
!---- reach one.
!----
!---- The common case -- a missing file -- is pre-empted below with an
!---- inquire(). For the rest, Stage 9 does not fork the readers (that would
!---- give up the scaleF/scaleM reuse that motivates using them) and does not
!---- add a validating pre-pass (that would amount to a second parser).
!---- Instead the *file* route stays exactly as it is, for xgf3d and for a
!---- Fortran caller, and the C/Python route never enters it: the caller
!---- passes values, and gf_source_set_cmt/gf_source_set_force apply the
!---- readers' own post-parse arithmetic to them, returning an error code
!---- where the reader would stop. Parsing a CMTSOLUTION is then the Python
!---- side's job, which is where a malformed file should be diagnosed anyway.
!----
!---- The two routes are pinned against each other by tests/gf3d/
!---- 4d.test_gf_source, at the suite's derived tolerance: they are separate
!---- compilation units evaluating the same expressions.
!----
!---- No `use hdf5`: this is a kernel module.
!----

  module gf_source

  use gf_par, only: t_gf_source,gf_set_error,gf_is_finite,gf_all_finite, &
                    GF_OK,GF_ERR_ARG,GF_ERR_NO_FILE,GF_ERR_FORMAT, &
                    GF_SRC_FORCE,GF_SRC_CMT

  implicit none

  private

  public :: gf_read_force_source
  public :: gf_read_cmt_source
  public :: gf_read_source
  public :: gf_source_set_cmt
  public :: gf_source_set_force
  public :: gf_detect_source_type
  public :: gf_force_direction
  public :: gf_print_source

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_read_force_source(filename,dt,src,ierr)

! reads a FORCESOLUTION through the solver's own get_force()
!
! `dt` is the solver time step the database was written with, from
! mesh_info.h5. get_force needs it because several of its source time
! function branches clamp hdur to 5*DT (get_force.f90:223,234,246), and using
! a different DT here than the reciprocal runs used would put a different
! floor on the source it reports.

  use constants, only: MAX_STRING_LEN

  implicit none

  character(len=*), intent(in) :: filename
  double precision, intent(in) :: dt
  type(t_gf_source), intent(out) :: src
  integer, intent(out) :: ierr

  ! one source per file: this library extracts for one source at a time, and
  ! get_force's NSOURCES > 1 machinery exists for the solver's multi-source
  ! runs, where it also rebases tshift_src onto the earliest source
  integer, parameter :: NSOURCES = 1

  ! local parameters
  double precision, dimension(NSOURCES) :: tshift_src,hdur,lat,lon,depth
  double precision, dimension(NSOURCES) :: factor_force_source
  double precision, dimension(NSOURCES) :: comp_E,comp_N,comp_Z_UP
  integer, dimension(NSOURCES) :: force_stf
  double precision :: min_tshift_src_original
  logical :: exists

  src = t_gf_source()

  if (len_trim(filename) == 0) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_read_force_source: no FORCESOLUTION path given')
    return
  endif

  if (dt <= 0.d0) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_read_force_source: the database time step must be positive')
    return
  endif

  ! get_force.f90:103 stops on a missing file, so catch that here: it is by
  ! far the most likely way a caller gets this wrong
  inquire(file=trim(filename),exist=exists)
  if (.not. exists) then
    call gf_set_error(ierr,GF_ERR_NO_FILE,'no such FORCESOLUTION file: '//trim(filename))
    return
  endif

  call get_force(tshift_src,hdur,lat,lon,depth,dt,NSOURCES, &
                 min_tshift_src_original,force_stf,factor_force_source, &
                 comp_E,comp_N,comp_Z_UP,trim(filename))

  src%source_type = GF_SRC_FORCE
  src%filename = filename

  src%latitude  = lat(1)
  src%longitude = lon(1)
  src%depth     = depth(1)

  ! for a force source this is the dominant frequency f0, not a half
  ! duration; the name is get_force's
  src%hdur = hdur(1)

  src%tshift_src = tshift_src(1)
  src%min_tshift_src_original = min_tshift_src_original

  src%force_stf = force_stf(1)
  src%factor_force_source = factor_force_source(1)
  src%comp_dir_vect_source_E    = comp_E(1)
  src%comp_dir_vect_source_N    = comp_N(1)
  src%comp_dir_vect_source_Z_UP = comp_Z_UP(1)

  ierr = GF_OK

  end subroutine gf_read_force_source

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_read_cmt_source(filename,dt,src,ierr)

! reads a CMTSOLUTION through the solver's own get_cmt()
!
! Three things live *outside* get_cmt and are silent bugs if a caller
! forgets them; all three are recorded here rather than in a stage document
! that a future reader may not have:
!
!  1. `hdur` comes back as the raw *triangle* half duration. The conversion
!     to a Gaussian width is hdur/SOURCE_DECAY_MIMIC_TRIANGLE and happens at
!     setup_sources_receivers.f90:777, in the caller. Miss it and every
!     trace is wrong by a factor of 1.628 in duration. gf_stf applies it in
!     Stage 5; nothing here does.
!  2. `tshift_src` is zeroed when NSOURCES == 1 and the original returned in
!     min_tshift_src_original (29.0 s for the shipped example), so t = 0 is
!     the centroid time. Those seconds are origin-time metadata for Stage
!     10's SAC headers, not part of the trace.
!  3. EXTERNAL_SOURCE_TIME_FUNCTION (constants.h:337, currently .false.)
!     zeroes hdur entirely if it is ever flipped.

  use constants, only: MAX_STRING_LEN,PI,GRAV

  use shared_parameters, only: RHOAV,R_PLANET

  implicit none

  character(len=*), intent(in) :: filename
  double precision, intent(in) :: dt
  type(t_gf_source), intent(out) :: src
  integer, intent(out) :: ierr

  ! one source per file, as for the force case
  integer, parameter :: NSOURCES = 1

  ! local parameters
  double precision, dimension(NSOURCES) :: tshift_src,hdur,lat,lon,depth
  double precision, dimension(6,NSOURCES) :: moment_tensor
  double precision :: min_tshift_src_original,sec
  integer :: yr,jda,mo,da,ho,mi
  logical :: exists

  src = t_gf_source()

  if (len_trim(filename) == 0) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_read_cmt_source: no CMTSOLUTION path given')
    return
  endif

  if (dt <= 0.d0) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_read_cmt_source: the database time step must be positive')
    return
  endif

  ! get_cmt.f90:103 stops on a missing file
  inquire(file=trim(filename),exist=exists)
  if (.not. exists) then
    call gf_set_error(ierr,GF_ERR_NO_FILE,'no such CMTSOLUTION file: '//trim(filename))
    return
  endif

  call get_cmt(yr,jda,mo,da,ho,mi,sec, &
               tshift_src,hdur,lat,lon,depth,moment_tensor, &
               dt,NSOURCES,min_tshift_src_original,trim(filename))

  src%source_type = GF_SRC_CMT
  src%filename = filename

  src%latitude  = lat(1)
  src%longitude = lon(1)
  src%depth     = depth(1)

  ! the raw triangle half duration -- see note 1 above
  src%hdur = hdur(1)

  src%tshift_src = tshift_src(1)
  src%min_tshift_src_original = min_tshift_src_original

  ! (Mrr,Mtt,Mpp,Mrt,Mrp,Mtp), non-dimensional: get_cmt.f90:426 has already
  ! divided by scaleM = 1.d7 * RHOAV * R_PLANET**5 * PI*GRAV*RHOAV
  src%moment_tensor(1:6) = moment_tensor(1:6,1)

  ! the same scaleM, in get_cmt's own expression (get_cmt.f90:426) from the
  ! same module variables, which gf_shared_params set from the database
  ! before the source was read. get_cmt does not return it, and the
  ! partials of Stage 6 are per dyne-cm, i.e. per unit of what the file
  ! says, so the factor has to be known.
  src%scale_moment = 1.d7 * RHOAV * (R_PLANET**5) * PI*GRAV*RHOAV

  src%yr = yr ; src%jda = jda ; src%mo = mo
  src%da = da ; src%ho = ho ; src%mi = mi
  src%sec = sec

  ! the one header field get_cmt reads past
  call gf_cmt_event_name(filename,src%event_name)

  ierr = GF_OK

  end subroutine gf_read_cmt_source

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_source_set_cmt(src,lat,lon,depth_km,hdur,time_shift,moment_dynecm,dt,ierr)

! builds a CMT source from values rather than from a file
!
! This is what the C facade calls, so that no path reachable from Python
! enters get_cmt() and its stop statements (see the module header).
!
! Everything get_cmt does to the numbers *after* it has parsed them is
! applied here, line by line, so that the two routes cannot drift:
!
!   get_cmt.f90:397     hdur = max(hdur, 5*DT). Unconditional in the current
!                       source: the USE_FORCE_POINT_SOURCE branch above it
!                       (:387-393) is commented out. A null half duration
!                       means a Heaviside, replaced by a very short error
!                       function.
!   get_cmt.f90:405     NOISE_TOMOGRAPHY /= 0 zeroes hdur entirely
!   get_cmt.f90:408     EXTERNAL_SOURCE_TIME_FUNCTION does the same
!                       (constants.h:337, .false.); both are no-ops under
!                       the settings gf_shared_params installs, and both are
!                       mirrored so that flipping either cannot make the two
!                       routes disagree
!   get_cmt.f90:411-413 NSOURCES == 1 zeroes tshift_src and returns the file's
!                       value in min_tshift_src_original, so t = 0 is the
!                       centroid time
!   get_cmt.f90:427-428 scaleM, and the moment tensor divided by it
!
! `moment_dynecm` is (Mrr,Mtt,Mpp,Mrt,Mrp,Mtp) in dyne-cm, i.e. the
! CMTSOLUTION's own numbers and units. `depth_km` is km, `hdur` the raw
! triangle half duration in seconds, `time_shift` the file's `time shift:`
! in seconds. `dt` is the database's solver step, needed for the clamp.
!
! The PDE header fields and the event name are left blank: they exist for
! Stage 10's SAC headers, which the in-memory API does not write.

  use constants, only: PI,GRAV,EXTERNAL_SOURCE_TIME_FUNCTION

  use shared_parameters, only: RHOAV,R_PLANET,NOISE_TOMOGRAPHY

  implicit none

  type(t_gf_source), intent(out) :: src
  double precision, intent(in) :: lat,lon,depth_km,hdur,time_shift,dt
  double precision, dimension(6), intent(in) :: moment_dynecm
  integer, intent(out) :: ierr

  ! local parameters
  double precision :: scale_moment,hdur_use

  src = t_gf_source()

  ! screened before anything arithmetic touches them: a NaN latitude would
  ! pass reduce()'s range test in gf_locate (both comparisons are false)
  ! and reach the kd-tree, which stops when it finds no point
  if (.not. (gf_is_finite(lat) .and. gf_is_finite(lon) .and. gf_is_finite(depth_km) &
       .and. gf_is_finite(hdur) .and. gf_is_finite(time_shift) .and. gf_is_finite(dt))) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_source_set_cmt: a source parameter is not finite')
    return
  endif
  if (.not. gf_all_finite(moment_dynecm)) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_source_set_cmt: a moment tensor component is not finite')
    return
  endif

  if (dt <= 0.d0) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_source_set_cmt: the database time step must be positive')
    return
  endif

  if (RHOAV <= 0.d0 .or. R_PLANET <= 0.d0) then
    call gf_set_error(ierr,GF_ERR_ARG, &
      'gf_source_set_cmt: planet constants are unset; open a database first')
    return
  endif

  src%source_type = GF_SRC_CMT
  src%filename = ''

  src%latitude  = lat
  src%longitude = lon
  src%depth     = depth_km

  ! get_cmt.f90:397, and then :405/:408
  hdur_use = hdur
  if (hdur_use < 5.d0 * dt) hdur_use = 5.d0 * dt
  if (NOISE_TOMOGRAPHY /= 0) hdur_use = 0.d0
  if (EXTERNAL_SOURCE_TIME_FUNCTION) hdur_use = 0.d0

  src%hdur = hdur_use

  ! get_cmt.f90:411-413
  src%tshift_src = 0.d0
  src%min_tshift_src_original = time_shift

  ! get_cmt.f90:427-428, in get_cmt's own expression from the same module
  ! variables gf_shared_params set from the database
  scale_moment = 1.d7 * RHOAV * (R_PLANET**5) * PI*GRAV*RHOAV

  src%moment_tensor(1:6) = moment_dynecm(1:6) / scale_moment
  src%scale_moment = scale_moment

  ierr = GF_OK

  end subroutine gf_source_set_cmt

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_source_set_force(src,lat,lon,depth_km,f0,time_shift,force_stf, &
                                 factor_newton,dir_E,dir_N,dir_Z_UP,dt,ierr)

! builds a force source from values rather than from a file
!
! The force half of gf_source_set_cmt, and the one that matters more for the
! no-stop invariant: get_force stops on three malformed inputs, and each of
! them is an error code here.
!
!   get_force.f90:218-249  per force_stf, the half-duration check:
!                            0, 2, 4  hdur = max(f0, 5*DT)
!                            1        Ricker, hdur = max(f0, TINYVAL)
!                            3        monochromatic; f0 < TINYVAL *stops*
!                            other    *stops*
!                          (the reader writes `5. * DT`, a default-real
!                           literal that is exactly 5, so 5.d0*dt agrees
!                           bit for bit)
!   get_force.f90:260      a second, unconditional hdur = max(hdur, TINYVAL)
!   get_force.f90:265-273  a (near) zero-length direction vector *stops*
!   get_force.f90:277-279  NSOURCES == 1 zeroes tshift_src
!   get_force.f90:289-290  scaleF, and the factor divided by it
!
! `f0` is the FORCESOLUTION's `f0:` field -- a dominant frequency for a
! Ricker, a Gaussian width otherwise -- which the reader stores in `hdur`.
! `factor_newton` is `factor force source:` in Newtons. The direction vector
! has arbitrary length; gf_force_direction normalises it.

  use constants, only: PI,GRAV,TINYVAL

  use shared_parameters, only: RHOAV,R_PLANET

  implicit none

  type(t_gf_source), intent(out) :: src
  double precision, intent(in) :: lat,lon,depth_km,f0,time_shift
  integer, intent(in) :: force_stf
  double precision, intent(in) :: factor_newton,dir_E,dir_N,dir_Z_UP,dt
  integer, intent(out) :: ierr

  ! local parameters
  double precision :: scaleF,hdur_use,norm

  src = t_gf_source()

  if (.not. (gf_is_finite(lat) .and. gf_is_finite(lon) .and. gf_is_finite(depth_km) &
       .and. gf_is_finite(f0) .and. gf_is_finite(time_shift) .and. gf_is_finite(dt) &
       .and. gf_is_finite(factor_newton) &
       .and. gf_is_finite(dir_E) .and. gf_is_finite(dir_N) .and. gf_is_finite(dir_Z_UP))) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_source_set_force: a source parameter is not finite')
    return
  endif

  if (dt <= 0.d0) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_source_set_force: the database time step must be positive')
    return
  endif

  if (RHOAV <= 0.d0 .or. R_PLANET <= 0.d0) then
    call gf_set_error(ierr,GF_ERR_ARG, &
      'gf_source_set_force: planet constants are unset; open a database first')
    return
  endif

  ! get_force.f90:218-249, with the two stops turned into error codes
  hdur_use = f0
  select case (force_stf)
  case (0,2,4)
    ! Gaussian, step, or Meschede Gaussian: a null width means a Dirac,
    ! replaced by a very short function
    if (hdur_use < 5.d0 * dt) hdur_use = 5.d0 * dt
  case (1)
    ! Ricker: f0 is the dominant frequency
    if (hdur_use < TINYVAL) hdur_use = TINYVAL
  case (3)
    ! monochromatic: f0 is the period, and there is no sensible default
    if (hdur_use < TINYVAL) then
      call gf_set_error(ierr,GF_ERR_ARG, &
        'gf_source_set_force: a monochromatic force (source time function 3) needs a non-zero period')
      return
    endif
  case default
    call gf_set_error(ierr,GF_ERR_ARG, &
      'gf_source_set_force: unsupported source time function type (force_stf), expected 0 to 4')
    return
  end select

  ! get_force.f90:260
  if (hdur_use < TINYVAL) hdur_use = TINYVAL

  ! get_force.f90:265-273
  norm = sqrt(dir_E**2 + dir_N**2 + dir_Z_UP**2)
  if (norm < TINYVAL) then
    call gf_set_error(ierr,GF_ERR_ARG, &
      'gf_source_set_force: the force direction vector has (almost) zero length')
    return
  endif

  src%source_type = GF_SRC_FORCE
  src%filename = ''

  src%latitude  = lat
  src%longitude = lon
  src%depth     = depth_km

  src%hdur = hdur_use

  ! get_force.f90:277-279
  src%tshift_src = 0.d0
  src%min_tshift_src_original = time_shift

  src%force_stf = force_stf

  ! get_force.f90:289-290
  scaleF = RHOAV * (R_PLANET**4) * PI*GRAV*RHOAV
  src%factor_force_source = factor_newton / scaleF

  src%comp_dir_vect_source_E    = dir_E
  src%comp_dir_vect_source_N    = dir_N
  src%comp_dir_vect_source_Z_UP = dir_Z_UP

  ierr = GF_OK

  end subroutine gf_source_set_force

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_cmt_event_name(filename,event_name)

! the `event name:` line of a CMTSOLUTION, or '' if there is none
!
! get_cmt parses everything else in the file but skips this line, and the
! solver takes the name from get_event_info_serial(), which reads a
! hard-coded DATA/CMTSOLUTION. The SAC KEVNM header (Stage 10) wants it, so
! it is scanned here: one line, one field, no numbers.

  use constants, only: MAX_STRING_LEN

  implicit none

  character(len=*), intent(in) :: filename
  character(len=*), intent(out) :: event_name

  ! local parameters
  character(len=MAX_STRING_LEN) :: line
  integer :: iunit,ios,icolon

  event_name = ''

  open(newunit=iunit,file=trim(filename),status='old',action='read',iostat=ios)
  if (ios /= 0) return

  do
    read(iunit,'(a)',iostat=ios) line
    if (ios /= 0) exit
    line = adjustl(line)
    if (line(1:11) == 'event name:') then
      icolon = index(line,':')
      event_name = adjustl(line(icolon+1:))
      exit
    endif
  enddo

  close(iunit)

  end subroutine gf_cmt_event_name

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_detect_source_type(filename,source_type,ierr)

! decides whether a file is a FORCESOLUTION or a CMTSOLUTION
!
! The two formats are distinguishable on their first non-blank line: a
! FORCESOLUTION opens with the literal token FORCE (get_force.f90 skips it as
! a header label), while a CMTSOLUTION opens with a PDE line whose first
! token is the data source, e.g. " PDE 1994 6 9 ...". Sniffing for FORCE is
! therefore both sufficient and stable -- and it beats asking the user to
! say which they meant, since getting it wrong would otherwise surface as a
! parse error deep inside a solver reader that answers with `stop`.

  use constants, only: MAX_STRING_LEN

  implicit none

  character(len=*), intent(in) :: filename
  integer, intent(out) :: source_type
  integer, intent(out) :: ierr

  ! local parameters
  character(len=MAX_STRING_LEN) :: line
  integer :: iin,ios
  logical :: exists

  source_type = 0

  inquire(file=trim(filename),exist=exists)
  if (.not. exists) then
    call gf_set_error(ierr,GF_ERR_NO_FILE,'no such source file: '//trim(filename))
    return
  endif

  open(newunit=iin,file=trim(filename),status='old',action='read',iostat=ios)
  if (ios /= 0) then
    call gf_set_error(ierr,GF_ERR_NO_FILE,'could not open the source file: '//trim(filename))
    return
  endif

  do
    read(iin,'(a)',iostat=ios) line
    if (ios /= 0) then
      close(iin)
      call gf_set_error(ierr,GF_ERR_FORMAT,'source file is empty: '//trim(filename))
      return
    endif
    if (len_trim(line) > 0) exit
  enddo

  close(iin)

  line = adjustl(line)
  if (line(1:5) == 'FORCE' .or. line(1:5) == 'force') then
    source_type = GF_SRC_FORCE
  else
    source_type = GF_SRC_CMT
  endif

  ierr = GF_OK

  end subroutine gf_detect_source_type

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_read_source(filename,dt,src,ierr)

! reads either kind of source file, deciding which from its contents

  implicit none

  character(len=*), intent(in) :: filename
  double precision, intent(in) :: dt
  type(t_gf_source), intent(out) :: src
  integer, intent(out) :: ierr

  ! local parameters
  integer :: source_type

  call gf_detect_source_type(filename,source_type,ierr)
  if (ierr /= GF_OK) return

  if (source_type == GF_SRC_FORCE) then
    call gf_read_force_source(filename,dt,src,ierr)
  else
    call gf_read_cmt_source(filename,dt,src,ierr)
  endif

  end subroutine gf_read_source

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_force_direction(src,nu,fhat,ierr)

! the Cartesian unit direction of a force source
!
! Transcribed from the solver, which is the authority on this and where the
! two halves live apart:
!
!   setup_sources_receivers.f90:1585-1596   normalise (E,N,Z_UP), then
!     comp_x = N/norm, comp_y = E/norm, comp_z = Z_UP/norm
!   compute_arrays_source.f90:172-177       contract with nu:
!     f = nu(1,:)*comp_x + nu(2,:)*comp_y + nu(3,:)*comp_z
!
! Note the normalisation: FORCESOLUTION gives a direction vector of arbitrary
! length and the magnitude lives entirely in factor_force_source, so the
! length of the given vector must not leak into the amplitude. The solver
! divides it out; so does this.
!
! `nu` rows are N, E, Z-up in Cartesian, as gf_source_nu builds them -- the
! comp_x/comp_y/comp_z ordering above pairs N with row 1 and E with row 2,
! which is why the solver's local names look transposed against the file's
! E/N/Z_UP field order. Getting this pairing backwards mirrors the source
! about the local meridian and is invisible until a forward comparison.

  use constants, only: NDIM,TINYVAL

  implicit none

  type(t_gf_source), intent(in) :: src
  double precision, dimension(NDIM,NDIM), intent(in) :: nu
  double precision, dimension(NDIM), intent(out) :: fhat
  integer, intent(out) :: ierr

  ! local parameters
  double precision :: norm,comp_x,comp_y,comp_z
  integer :: i

  fhat(:) = 0.d0

  norm = sqrt(src%comp_dir_vect_source_E**2 &
            + src%comp_dir_vect_source_N**2 &
            + src%comp_dir_vect_source_Z_UP**2)

  ! get_force.f90:272 stops on this; reached only if a caller built a
  ! t_gf_source by hand, since the reader would have stopped first
  if (norm < TINYVAL) then
    call gf_set_error(ierr,GF_ERR_ARG, &
      'force source has a (almost) zero-length direction vector')
    return
  endif

  comp_x = src%comp_dir_vect_source_N / norm
  comp_y = src%comp_dir_vect_source_E / norm
  comp_z = src%comp_dir_vect_source_Z_UP / norm

  do i = 1,NDIM
    fhat(i) = nu(1,i)*comp_x + nu(2,i)*comp_y + nu(3,i)*comp_z
  enddo

  ierr = GF_OK

  end subroutine gf_force_direction

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_print_source(src,iunit)

! reports a source, in the 'key = value' layout gf_print_info uses

  implicit none

  type(t_gf_source), intent(in) :: src
  integer, intent(in) :: iunit

  write(iunit,'(a)')         'source'
  write(iunit,'(a,a)')       '  file                 = ',trim(src%filename)
  if (src%source_type == GF_SRC_FORCE) then
    write(iunit,'(a)')       '  type                 = FORCESOLUTION'
  else
    write(iunit,'(a)')       '  type                 = CMTSOLUTION'
  endif
  write(iunit,'(a,es22.14)') '  latitude             = ',src%latitude
  write(iunit,'(a,es22.14)') '  longitude            = ',src%longitude
  write(iunit,'(a,es22.14)') '  depth, km            = ',src%depth
  write(iunit,'(a,es22.14)') '  tshift_src           = ',src%tshift_src
  write(iunit,'(a,es22.14)') '  origin time shift    = ',src%min_tshift_src_original

  if (src%source_type == GF_SRC_FORCE) then
    write(iunit,'(a,i0)')      '  force_stf            = ',src%force_stf
    write(iunit,'(a,es22.14)') '  f0                   = ',src%hdur
    write(iunit,'(a,es22.14)') '  factor (non-dim)     = ',src%factor_force_source
    write(iunit,'(a,es22.14)') '  comp dir vect E      = ',src%comp_dir_vect_source_E
    write(iunit,'(a,es22.14)') '  comp dir vect N      = ',src%comp_dir_vect_source_N
    write(iunit,'(a,es22.14)') '  comp dir vect Z_UP   = ',src%comp_dir_vect_source_Z_UP
  else
    ! the raw triangle half duration; the /1.628 conversion to a Gaussian
    ! width is Stage 5's, so both are printed to keep them distinguishable
    write(iunit,'(a,es22.14)') '  hdur (triangle)      = ',src%hdur
    write(iunit,'(a,es22.14)') '  Mrr (non-dim)        = ',src%moment_tensor(1)
    write(iunit,'(a,es22.14)') '  Mtt (non-dim)        = ',src%moment_tensor(2)
    write(iunit,'(a,es22.14)') '  Mpp (non-dim)        = ',src%moment_tensor(3)
    write(iunit,'(a,es22.14)') '  Mrt (non-dim)        = ',src%moment_tensor(4)
    write(iunit,'(a,es22.14)') '  Mrp (non-dim)        = ',src%moment_tensor(5)
    write(iunit,'(a,es22.14)') '  Mtp (non-dim)        = ',src%moment_tensor(6)
  endif

  end subroutine gf_print_source

  end module gf_source
