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
!---- xgf3d — the Green function extraction tool.
!----
!---- One executable with several modes rather than a family of small
!---- binaries, so that every mode survives as a debugging entry point:
!----
!----   xgf3d --version
!----   xgf3d --info          <GFDB> [--topo] [--no-check]
!----   xgf3d --locate        <GFDB> <lat> <lon> <depth_km>
!----   xgf3d --check-anchors <GFDB>
!----   xgf3d --seis          <GFDB> <SOURCE> <outdir> [--format ascii|sac|sacan|all]
!----                                                  [--partials 1|2] [--t0 <seconds>]
!----   xgf3d --dump          <GFDB> <SOURCE> <outdir> [--station NET.STA]
!----
!---- <SOURCE> is a FORCESOLUTION or a CMTSOLUTION; which one is decided
!---- from the file's own first non-blank line, so the caller never has to
!---- say. --seis is *the* extraction mode: --format picks the writer (SAC
!---- by default, the deliverable; ASCII is what the comparison harness
!---- reads), and --partials adds the derivative products in the same
!---- format, so there is no separate SAC mode and no --seis-partials.
!----
!---- This is the only file under src/gf3d/ that may stop, and it does so in
!---- exactly two places: die() for a failure and usage_error() for a
!---- malformed command line. The library itself returns codes, because a
!---- stop inside lib/libgf3d.so would kill the Python interpreter that
!---- loaded it.
!----

  program xgf3d

! the same module a program outside this tree uses; the executable is a
! caller of the library, not a privileged part of it
  use gf3d, only: t_gfdb,t_gf_location,t_gf_source,t_gf_stf,t_gf_taxis, &
                  gf_errmsg,gf_error_string, &
                  GF_OK,GF3D_VERSION,GF_ANCHOR_TOL,GF_NCOMP,GF_SRC_CMT, &
                  gf_open,gf_close,gf_release,gf_print_info, &
                  gf_locate_source,gf_locate_release,gf_check_anchors_all, &
                  gf_print_location, &
                  gf_read_source,gf_print_source, &
                  gf_extract, &
                  gf_write_seis,gf_write_partials,gf_write_dump, &
                  gf_write_sac,gf_print_stf,gf_default_t0

  use constants, only: MAX_STRING_LEN,NGLLX,NGNOD

  implicit none

  ! local parameters
  character(len=MAX_STRING_LEN) :: arg,mode,dbpath,srcfile,outdir,station,format
  type(t_gfdb) :: db
  type(t_gf_location) :: loc
  type(t_gf_source) :: src
  type(t_gf_taxis) :: tax
  type(t_gf_stf) :: stf
  integer :: nargs,iarg,ierr,ios
  logical :: with_topo,do_check,want_ascii,want_sac,want_sacan
  logical :: have_station,have_t0,have_partials,have_format
  double precision :: lat,lon,depth_km,worst_err,t0_req
  integer :: ielem_worst,ista,ista_worst,itypsokern,ndp
  double precision, dimension(:,:,:), allocatable :: seis
  double precision, dimension(:,:,:,:), allocatable :: dp
  double precision, dimension(:), allocatable :: tsec,onset

  ! the mode/option cross-check table, filled below
  integer, parameter :: NXOPT = 4
  logical, dimension(NXOPT) :: opt_given,opt_for_seis
  character(len=48), dimension(NXOPT) :: opt_msg
  integer :: iopt

  ! standard error, used for anything that is not the requested output
  integer, parameter :: ISTDERR = 0

  nargs = command_argument_count()

  if (nargs < 1) call usage_error()

  call get_command_argument(1,mode)

  select case (trim(mode))

  case ('--version','-v')
    write(*,'(a)') 'xgf3d version '//GF3D_VERSION//' (specfem3d_globe)'
#ifdef USE_HDF5
    write(*,'(a)') '  HDF5 support : yes'
#else
    write(*,'(a)') '  HDF5 support : no  (re-run configure with --with-hdf5)'
#endif
    write(*,'(a,i0,a,i0)') '  built with     NGLLX = ',NGLLX,'  NGNOD = ',NGNOD
    stop

  case ('--help','-h')
    call print_usage(6)
    stop

  case ('--info')

    if (nargs < 2) call usage_error('Error: --info needs a database directory')

    call get_command_argument(2,dbpath)

    with_topo = .false.
    do_check = .true.
    do iarg = 3,nargs
      call get_command_argument(iarg,arg)
      select case (trim(arg))
      case ('--topo')
        with_topo = .true.
      case ('--no-check')
        do_check = .false.
      case default
        call die('Error: unknown option for --info: '//trim(arg))
      end select
    enddo

    call gf_open(dbpath,db,ierr,check_completion=do_check)
    if (ierr /= GF_OK) call die('Error opening the Green function database',ierr)

    call gf_print_info(db,6,with_topo)

    call gf_close(db)

  case ('--locate')

    if (nargs < 5) call usage_error( &
      'Error: --locate needs a database directory, a latitude, a longitude and a depth in km')

    call get_command_argument(2,dbpath)

    call read_double_arg(3,'latitude',lat)
    call read_double_arg(4,'longitude',lon)
    call read_double_arg(5,'depth',depth_km)

    call gf_open(dbpath,db,ierr,check_completion=.false.)
    if (ierr /= GF_OK) call die('Error opening the Green function database',ierr)

    call gf_locate_source(db,lat,lon,depth_km,loc,ierr)
    if (ierr /= GF_OK) call die('Error locating the source',ierr)

    call gf_print_location(db,loc,lat,lon,depth_km,6)

    call gf_locate_release()
    call gf_close(db)

  case ('--seis','--dump')

    ! xgf3d --seis <GFDB> <FORCESOLUTION> <outdir>
    ! xgf3d --dump <GFDB> <FORCESOLUTION> <outdir> [--station NET.STA]
    !
    ! Both modes open the database and read the source the same way; --dump
    ! then locates it directly and dumps the stored quantities on the
    ! database's own axis, while --seis calls gf_extract, which locates,
    ! plans and extracts in one sequence -- the same one gf3d_capi and a
    ! Fortran caller run. --dump exists because the manufactured-solution
    ! tests pin the operators on synthetic input, but nothing except a
    ! forward comparison can show that the HDF5 layout was read in the right
    ! index order or that the nu convention is right, and those fail
    ! non-locally.

    if (nargs < 4) call usage_error( &
      'Error: '//trim(mode)//' needs a database, a source file and an output directory')

    call get_command_argument(2,dbpath)
    call get_command_argument(3,srcfile)
    call get_command_argument(4,outdir)

    station = ''
    ! resolved below, once the source is read: absent --t0 means specfem's
    ! own rule for this source
    t0_req = -1.d0
    itypsokern = 0
    ! SAC is the deliverable and the default; the comparison harness asks
    ! for ascii explicitly (the Snakefiles under EXAMPLES/green_function_database)
    format = 'sac'
    have_station = .false.
    have_t0 = .false.
    have_partials = .false.
    have_format = .false.
    iarg = 5
    do while (iarg <= nargs)
      call get_command_argument(iarg,arg)
      select case (trim(arg))
      case ('--format')
        if (iarg == nargs) call die('Error: --format needs one of ascii, sac, sacan, all')
        iarg = iarg + 1
        call get_command_argument(iarg,format)
        select case (trim(format))
        case ('ascii','sac','sacan','all')
          continue
        case default
          call die('Error: --format takes ascii, sac, sacan or all, not '//trim(format))
        end select
        have_format = .true.
      case ('--station')
        if (iarg == nargs) call die('Error: --station needs a NET.STA identifier')
        iarg = iarg + 1
        call get_command_argument(iarg,station)
        have_station = .true.
      case ('--t0')
        if (iarg == nargs) call die('Error: --t0 needs a start time in seconds before the origin')
        iarg = iarg + 1
        call read_double_arg(iarg,'t0',t0_req)
        have_t0 = .true.
      case ('--partials')
        if (iarg == nargs) call die( &
          'Error: --partials needs a kernel type: 1 (moment tensor) or 2 (and centroid)')
        iarg = iarg + 1
        call get_command_argument(iarg,arg)
        read(arg,*,iostat=ios) itypsokern
        if (ios /= 0 .or. itypsokern < 1 .or. itypsokern > 2) call die( &
          'Error: --partials takes 1 (moment tensor) or 2 (and centroid), not '//trim(arg))
        have_partials = .true.
      case default
        call die('Error: unknown option for '//trim(mode)//': '//trim(arg))
      end select
      iarg = iarg + 1
    enddo

    ! which mode each option belongs to, as a table over the four flags the
    ! parser set, in the order the four tests used to be made.
    !
    ! --dump writes the stored quantities on the database's own axis by
    ! design, so the output axis, the partials and the writer choice are all
    ! --seis's. The last of these used to be written `format /= 'sac'`, which
    ! could not tell an explicit --format sac from the default and therefore
    ! let --dump --format sac through; have_format can.
    opt_given(1:NXOPT)    = (/ have_station, have_t0, have_partials, have_format /)
    opt_for_seis(1:NXOPT) = (/ .false.,      .true.,  .true.,        .true.      /)
    opt_msg(1) = 'Error: --station applies to --dump only'
    opt_msg(2) = 'Error: --t0 applies to --seis only'
    opt_msg(3) = 'Error: --partials applies to --seis only'
    opt_msg(4) = 'Error: --format applies to --seis only'

    do iopt = 1,NXOPT
      if (opt_given(iopt) .and. (opt_for_seis(iopt) .neqv. (trim(mode) == '--seis'))) &
        call die(trim(opt_msg(iopt)))
    enddo

    want_ascii = (trim(format) == 'ascii' .or. trim(format) == 'all')
    want_sac   = (trim(format) == 'sac'   .or. trim(format) == 'all')
    want_sacan = (trim(format) == 'sacan' .or. trim(format) == 'all')

    call gf_open(dbpath,db,ierr,check_completion=.false.)
    if (ierr /= GF_OK) call die('Error opening the Green function database',ierr)

    ! FORCESOLUTION or CMTSOLUTION, decided from the file's own first line
    call gf_read_source(db,srcfile,src,ierr)
    if (ierr /= GF_OK) call die('Error reading the source',ierr)

    call gf_print_source(src,6)
    write(*,'(a)') ''

    ! cross-check #5: needs src, so it sits after the read
    if (itypsokern > 0 .and. src%source_type /= GF_SRC_CMT) call die( &
      'Error: --partials is defined for a CMTSOLUTION, not a FORCESOLUTION')

    if (trim(mode) == '--dump') then

      call gf_locate_source(db,src%latitude,src%longitude,src%depth,loc,ierr)
      if (ierr /= GF_OK) call die('Error locating the source',ierr)

      call gf_print_location(db,loc,src%latitude,src%longitude,src%depth,6)
      write(*,'(a)') ''

      call gf_write_dump(db,src,loc,outdir,station,ierr)
      if (ierr /= GF_OK) call die('Error writing the dump',ierr)
      write(*,'(a,a)') 'wrote interpolated displacement to ',trim(outdir)

    else

      ! --t0 left off asks for specfem's own start time for this source;
      ! gf_seis_plan itself takes a resolved value only
      if (.not. have_t0) then
        call gf_default_t0(src,t0_req,ierr)
        if (ierr /= GF_OK) call die('Error choosing the default start time',ierr)
      endif

      ! locate, plan and extract in one call -- the same sequence gf3d_capi
      ! and a Fortran caller run. loc and stf are valid for every stage that
      ! completed even when ierr /= GF_OK, so a failure inside the
      ! extraction itself still reports where the source was and what
      ! conversion was planned.
      call gf_extract(db,src,t0_req,itypsokern,seis,dp,ierr,t=tsec,onset=onset, &
                       loc=loc,tax=tax,stf=stf)

      if (loc%ielem > 0) then
        call gf_print_location(db,loc,src%latitude,src%longitude,src%depth,6)
        write(*,'(a)') ''
      endif
      if (tax%nt > 0) then
        call gf_print_stf(stf,tax,6)
        write(*,'(a)') ''
      endif
      if (ierr /= GF_OK) call die('Error computing the seismograms',ierr)

      ! the partials array exists in every case (zero-size without
      ! --partials), so that the SAC writer has one interface
      ndp = size(dp,1)

      !--- the writers: every requested format, every product -------------

      if (want_ascii) then
        call gf_write_seis(db,src,loc,tax,stf,seis,tsec,onset,outdir,ierr)
        if (ierr /= GF_OK) call die('Error writing the seismograms',ierr)
        write(*,'(a,i0,a,a)') 'wrote ',db%nstations,' ASCII seismogram files to ',trim(outdir)

        if (itypsokern > 0) then
          call gf_write_partials(db,src,loc,tax,stf,ndp,dp,tsec,outdir,ierr)
          if (ierr /= GF_OK) call die('Error writing the partials',ierr)
          write(*,'(a,i0,a,i0,a,a)') 'wrote ',ndp,' ASCII partials per component for ',db%nstations, &
                                     ' stations to ',trim(outdir)
        endif
      endif

      if (want_sac .or. want_sacan) then
        call gf_write_sac(db,src,tax,stf,seis,ndp,dp,outdir,want_sac,want_sacan,ierr)
        if (ierr /= GF_OK) call die('Error writing the SAC files',ierr)
        write(*,'(a,i0,a,i0,a,a)') 'wrote ',db%nstations*GF_NCOMP*(1 + ndp),' SAC traces (', &
                                   ndp,' partials per component) to ',trim(outdir)
      endif

      deallocate(dp)

      ! the silence-before-the-record ratio the conversion rests on
      ista_worst = 1
      do ista = 2,db%nstations
        if (onset(ista) > onset(ista_worst)) ista_worst = ista
      enddo
      write(*,'(a,es22.14,a,a)') 'onset ratio (worst)  = ',onset(ista_worst), &
                                 '  at ',trim(db%stations(ista_worst)%id)
      if (onset(ista_worst) > 1.d-3) then
        write(*,'(a)') '  WARNING: the trace is not silent before the reciprocal source switches on;'
        write(*,'(a)') '           the conversion extends it with zeros there and will be in error'
      endif

      deallocate(seis,tsec,onset)

    endif

    call gf_locate_release()
    call gf_close(db)

  case ('--check-anchors')

    if (nargs < 2) call usage_error('Error: --check-anchors needs a database directory')

    call get_command_argument(2,dbpath)

    call gf_open(dbpath,db,ierr,check_completion=.false.)
    if (ierr /= GF_OK) call die('Error opening the Green function database',ierr)

    call gf_check_anchors_all(db,worst_err,ielem_worst,ierr)
    if (ierr /= GF_OK) call die('Error checking the element anchors',ierr)

    write(*,'(a)')             'anchor consistency check'
    write(*,'(a,i0)')          '  elements checked     = ',db%nelem
    write(*,'(a,es22.14)')     '  worst residual       = ',worst_err
    if (ielem_worst > 0) then
      write(*,'(a,a)')         '  worst element        = ',db%morton_hex(ielem_worst)
    endif
    write(*,'(a,es22.14)')     '  tolerance            = ',GF_ANCHOR_TOL
    write(*,'(a,es22.14)')     '  worst residual, m    = ',worst_err*db%R_PLANET
    ! the verdict is two lines: the first is this mode's own finding, the
    ! second is what to do about it, so die() writes only the second
    if (worst_err > GF_ANCHOR_TOL) then
      write(ISTDERR,'(a)') 'the 27 anchors do not reproduce the stored GLL coordinates'
      call die('  this database was probably written from a USE_GLL = .true. mesh')
    endif
    write(*,'(a)')             '  result               = ok'

    call gf_close(db)

  case default
    call usage_error('Error: unknown mode: '//trim(mode))

  end select

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine die(context,ierr)

! reports a failure and leaves, with the database released
!
! The one way this program fails. `context` is what was being attempted; the
! optional `ierr` adds the library's own code and message on a second,
! indented line, which is the format every error path here has always
! printed and which the test runners read.
!
! gf_release() rather than gf_close(): the kd-tree is process-wide and
! gf_close() cannot reach it. Both are no-ops on a handle that was never
! opened (db%is_open defaults to .false.) and the release is a no-op without
! a tree, so this is safe on every path, including the ones that die before
! gf_open.

  implicit none

  character(len=*), intent(in) :: context
  integer, intent(in), optional :: ierr

  write(ISTDERR,'(a)') context
  if (present(ierr)) then
    write(ISTDERR,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
  endif

  if (db%is_open) call gf_release(db)

  stop 1

  end subroutine die

!
!-------------------------------------------------------------------------------------------------
!

  subroutine usage_error(message)

! reports a malformed command line, with the usage summary
!
! `message` is optional because the no-argument invocation is not a mistake
! to explain, only a request for the usage -- and it prints none.

  implicit none

  character(len=*), intent(in), optional :: message

  if (present(message)) write(ISTDERR,'(a)') message

  call print_usage(ISTDERR)

  stop 1

  end subroutine usage_error

!
!-------------------------------------------------------------------------------------------------
!

  subroutine print_usage(iunit)

  implicit none

  integer, intent(in) :: iunit

  write(iunit,'(a)') 'usage: xgf3d <mode> [arguments]'
  write(iunit,'(a)') ''
  write(iunit,'(a)') '  --version                       print the version and build settings'
  write(iunit,'(a)') '  --help                          print this message'
  write(iunit,'(a)') '  --info <GFDB> [options]         describe a Green function database'
  write(iunit,'(a)') '  --locate <GFDB> <lat> <lon> <depth_km>'
  write(iunit,'(a)') '                                  locate a position in the database'
  write(iunit,'(a)') '  --check-anchors <GFDB>          verify the 27-anchor geometry of every element'
  write(iunit,'(a)') '  --seis <GFDB> <SOURCE> <outdir> [--format ascii|sac|sacan|all] [--partials 1|2]'
  write(iunit,'(a)') '                                  [--t0 <seconds>]'
  write(iunit,'(a)') '                                  seismograms at every station, converted to the'
  write(iunit,'(a)') '                                  source time function specfem would use'
  write(iunit,'(a)') '  --dump <GFDB> <SOURCE> <outdir> [--station NET.STA]'
  write(iunit,'(a)') '                                  interpolated displacement, and for a CMT source'
  write(iunit,'(a)') '                                  the strain and the trace before the conversion'
  write(iunit,'(a)') ''
  write(iunit,'(a)') '  <SOURCE> is a FORCESOLUTION or a CMTSOLUTION, detected from its contents'
  write(iunit,'(a)') ''
  write(iunit,'(a)') '  options for --info:'
  write(iunit,'(a)') '    --topo       also load the topography grid and probe it at each station'
  write(iunit,'(a)') '    --no-check   skip the per-element completion scan'
  write(iunit,'(a)') ''
  write(iunit,'(a)') '  options for --seis:'
  write(iunit,'(a)') '    --format F   sac (default): NET.STA.BXN.sem.sac binary SAC, the solver''s header'
  write(iunit,'(a)') '                 rules; sacan: alphanumeric SAC; ascii: NET.STA.gf3d.txt columns,'
  write(iunit,'(a)') '                 what the comparison harness reads; all: every format'
  write(iunit,'(a)') '    --t0 <s>     start the output axis at or before <s> seconds before the origin,'
  write(iunit,'(a)') '                 extending the stored axis with zeros (default: 1.5*hdur, the'
  write(iunit,'(a)') '                 forward run''s own); the stored samples are never resampled'
  write(iunit,'(a)') '    --partials N also write the partial derivatives of a CMTSOLUTION''s seismograms,'
  write(iunit,'(a)') '                 in the same format (NET.STA.BXN.Mrr.sem.sac, NET.STA.partials.txt):'
  write(iunit,'(a)') '                 N = 1 the six moment-tensor components (m per dyne-cm); N = 2 also'
  write(iunit,'(a)') '                 latitude, longitude, depth and centroid time (m per degree,'
  write(iunit,'(a)') '                 degree, km, second)'
  write(iunit,'(a)') ''
  write(iunit,'(a)') '  the output directory must exist; xgf3d does not create it'
  write(iunit,'(a)') ''

  end subroutine print_usage

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_double_arg(iarg_in,name,val)

! reads a command-line argument as a double, or dies with a message

  implicit none

  integer, intent(in) :: iarg_in
  character(len=*), intent(in) :: name
  double precision, intent(out) :: val

  ! local parameters
  character(len=MAX_STRING_LEN) :: str
  integer :: ios

  call get_command_argument(iarg_in,str)
  read(str,*,iostat=ios) val
  if (ios /= 0) call die('Error: could not read the '//trim(name)//' from "'//trim(str)//'"')

  end subroutine read_double_arg

  end program xgf3d
