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
!---- Known wart, deferred to Stage 9
!---- -------------------------------
!---- get_force `stop`s on malformed input (get_force.f90:240,248,272), which
!---- breaks this library's no-stop invariant. The common case -- a missing
!---- file, :103 -- is pre-empted below with an inquire(); the rest need a
!---- validating pre-pass that would amount to a second parser, and forking
!---- the reader would give up the scaleF reuse that motivates using it. The
!---- trade is recorded rather than hidden.
!----
!---- No `use hdf5`: this is a kernel module.
!----

  module gf_source

  use gf_par, only: t_gf_source,gf_set_error, &
                    GF_OK,GF_ERR_ARG,GF_ERR_NO_FILE,GF_ERR_FORMAT, &
                    GF_SRC_FORCE,GF_SRC_CMT

  implicit none

  private

  public :: gf_read_force_source
  public :: gf_read_cmt_source
  public :: gf_read_source
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

  use constants, only: MAX_STRING_LEN

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

  src%yr = yr ; src%jda = jda ; src%mo = mo
  src%da = da ; src%ho = ho ; src%mi = mi
  src%sec = sec

  ierr = GF_OK

  end subroutine gf_read_cmt_source

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
