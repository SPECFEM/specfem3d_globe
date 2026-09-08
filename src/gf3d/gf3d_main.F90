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

  program xgf3d

  use gf_par, only: t_gfdb,t_gf_location,t_gf_source,t_gf_stf,t_gf_taxis, &
                    gf_errmsg,gf_error_string, &
                    GF_OK,GF3D_VERSION,GF_XI_TOL,GF_ANCHOR_TOL,GF_NCOMP,GF_SRC_CMT
  use gf_database, only: gf_open,gf_close,gf_print_info
  use gf_locate, only: gf_locate_source,gf_locate_release,gf_check_anchors_all
  use gf_source, only: gf_read_source,gf_print_source
  use gf_seismograms, only: gf_seis_plan,gf_seis,gf_seis_cmt_partials, &
                            gf_write_seis,gf_write_partials,gf_write_dump
  use gf_partials, only: gf_partials_ndp
  use gf_sac, only: gf_write_sac
  use gf_stf, only: gf_print_stf

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
  logical :: with_topo,do_check,have_t0,want_ascii,want_sac,want_sacan
  double precision :: lat,lon,depth_km,worst_err,t0_req
  integer :: ielem_worst,ista,ista_worst,itypsokern,ndp
  double precision, dimension(:,:,:), allocatable :: seis
  double precision, dimension(:,:,:,:), allocatable :: dp
  double precision, dimension(:), allocatable :: tsec,onset

  ! standard error, used for anything that is not the requested output
  integer, parameter :: ISTDERR = 0

  nargs = command_argument_count()

  if (nargs < 1) then
    call print_usage(ISTDERR)
    stop 1
  endif

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

    if (nargs < 2) then
      write(ISTDERR,'(a)') 'Error: --info needs a database directory'
      call print_usage(ISTDERR)
      stop 1
    endif

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
        write(ISTDERR,'(a)') 'Error: unknown option for --info: '//trim(arg)
        stop 1
      end select
    enddo

    call gf_open(dbpath,db,ierr,check_completion=do_check)
    if (ierr /= GF_OK) then
      write(ISTDERR,'(a)') 'Error opening the Green function database'
      write(ISTDERR,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
      stop 1
    endif

    call gf_print_info(db,6,with_topo)

    call gf_close(db)

  case ('--locate')

    if (nargs < 5) then
      write(ISTDERR,'(a)') 'Error: --locate needs a database directory, a latitude, a longitude and a depth in km'
      call print_usage(ISTDERR)
      stop 1
    endif

    call get_command_argument(2,dbpath)

    call read_double_arg(3,'latitude',lat)
    call read_double_arg(4,'longitude',lon)
    call read_double_arg(5,'depth',depth_km)

    call gf_open(dbpath,db,ierr,check_completion=.false.)
    if (ierr /= GF_OK) then
      write(ISTDERR,'(a)') 'Error opening the Green function database'
      write(ISTDERR,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
      stop 1
    endif

    call gf_locate_source(db,lat,lon,depth_km,loc,ierr)
    if (ierr /= GF_OK) then
      write(ISTDERR,'(a)') 'Error locating the source'
      write(ISTDERR,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
      call gf_locate_release()
      call gf_close(db)
      stop 1
    endif

    call print_location(db,loc,lat,lon,depth_km,6)

    call gf_locate_release()
    call gf_close(db)

  case ('--seis','--dump')

    ! xgf3d --seis <GFDB> <FORCESOLUTION> <outdir>
    ! xgf3d --dump <GFDB> <FORCESOLUTION> <outdir> [--station NET.STA]
    !
    ! Both modes share everything up to the located source; they differ only
    ! in what they write. --dump exists to localise a Stage 5 disagreement:
    ! the manufactured-solution tests pin the operators on synthetic input,
    ! but nothing except a forward comparison can show that the HDF5 layout
    ! was read in the right index order or that the nu convention is right,
    ! and those fail non-locally.

    if (nargs < 4) then
      write(ISTDERR,'(a)') 'Error: '//trim(mode)//' needs a database, a source file and an output directory'
      call print_usage(ISTDERR)
      stop 1
    endif

    call get_command_argument(2,dbpath)
    call get_command_argument(3,srcfile)
    call get_command_argument(4,outdir)

    station = ''
    ! a negative request asks gf_seis_plan for specfem's own rule, 1.5*hdur
    t0_req = -1.d0
    have_t0 = .false.
    itypsokern = 0
    ! SAC is the deliverable and the default; the comparison harness asks
    ! for ascii explicitly (the Snakefiles under EXAMPLES/green_function_database)
    format = 'sac'
    iarg = 5
    do while (iarg <= nargs)
      call get_command_argument(iarg,arg)
      select case (trim(arg))
      case ('--format')
        if (iarg == nargs) then
          write(ISTDERR,'(a)') 'Error: --format needs one of ascii, sac, sacan, all'
          stop 1
        endif
        iarg = iarg + 1
        call get_command_argument(iarg,format)
        select case (trim(format))
        case ('ascii','sac','sacan','all')
          continue
        case default
          write(ISTDERR,'(a)') 'Error: --format takes ascii, sac, sacan or all, not '//trim(format)
          stop 1
        end select
      case ('--station')
        if (iarg == nargs) then
          write(ISTDERR,'(a)') 'Error: --station needs a NET.STA identifier'
          stop 1
        endif
        iarg = iarg + 1
        call get_command_argument(iarg,station)
      case ('--t0')
        if (iarg == nargs) then
          write(ISTDERR,'(a)') 'Error: --t0 needs a start time in seconds before the origin'
          stop 1
        endif
        iarg = iarg + 1
        call read_double_arg(iarg,'t0',t0_req)
        have_t0 = .true.
      case ('--partials')
        if (iarg == nargs) then
          write(ISTDERR,'(a)') 'Error: --partials needs a kernel type: 1 (moment tensor) or 2 (and centroid)'
          stop 1
        endif
        iarg = iarg + 1
        call get_command_argument(iarg,arg)
        read(arg,*,iostat=ios) itypsokern
        if (ios /= 0 .or. itypsokern < 1 .or. itypsokern > 2) then
          write(ISTDERR,'(a)') 'Error: --partials takes 1 (moment tensor) or 2 (and centroid), not '//trim(arg)
          stop 1
        endif
      case default
        write(ISTDERR,'(a)') 'Error: unknown option for '//trim(mode)//': '//trim(arg)
        stop 1
      end select
      iarg = iarg + 1
    enddo

    if (trim(mode) == '--seis' .and. len_trim(station) > 0) then
      write(ISTDERR,'(a)') 'Error: --station applies to --dump only'
      stop 1
    endif
    ! --dump writes the stored quantities on the database's own axis by
    ! design; the output axis is a property of the seismogram
    if (trim(mode) == '--dump' .and. have_t0) then
      write(ISTDERR,'(a)') 'Error: --t0 applies to --seis only'
      stop 1
    endif
    if (trim(mode) == '--dump' .and. itypsokern > 0) then
      write(ISTDERR,'(a)') 'Error: --partials applies to --seis only'
      stop 1
    endif
    if (trim(mode) == '--dump' .and. trim(format) /= 'sac') then
      write(ISTDERR,'(a)') 'Error: --format applies to --seis only'
      stop 1
    endif
    want_ascii = (trim(format) == 'ascii' .or. trim(format) == 'all')
    want_sac   = (trim(format) == 'sac'   .or. trim(format) == 'all')
    want_sacan = (trim(format) == 'sacan' .or. trim(format) == 'all')

    call gf_open(dbpath,db,ierr,check_completion=.false.)
    if (ierr /= GF_OK) then
      write(ISTDERR,'(a)') 'Error opening the Green function database'
      write(ISTDERR,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
      stop 1
    endif

    ! FORCESOLUTION or CMTSOLUTION, decided from the file's own first line
    call gf_read_source(srcfile,db%dt,src,ierr)
    if (ierr /= GF_OK) then
      write(ISTDERR,'(a)') 'Error reading the source'
      write(ISTDERR,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
      call gf_close(db)
      stop 1
    endif

    call gf_locate_source(db,src%latitude,src%longitude,src%depth,loc,ierr)
    if (ierr /= GF_OK) then
      write(ISTDERR,'(a)') 'Error locating the source'
      write(ISTDERR,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
      call gf_locate_release()
      call gf_close(db)
      stop 1
    endif

    call gf_print_source(src,6)
    write(*,'(a)') ''
    call print_location(db,loc,src%latitude,src%longitude,src%depth,6)
    write(*,'(a)') ''

    if (trim(mode) == '--dump') then

      call gf_write_dump(db,src,loc,outdir,station,ierr)
      if (ierr /= GF_OK) then
        write(ISTDERR,'(a)') 'Error writing the dump'
        write(ISTDERR,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
        call gf_locate_release()
        call gf_close(db)
        stop 1
      endif
      write(*,'(a,a)') 'wrote interpolated displacement to ',trim(outdir)

    else

      ! the conversion and the axis are decided before any element is read,
      ! so a bad request fails in milliseconds and the plan is on record
      call gf_seis_plan(db,src,t0_req,tax,stf,ierr)
      if (ierr /= GF_OK) then
        write(ISTDERR,'(a)') 'Error planning the seismograms'
        write(ISTDERR,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
        call gf_locate_release()
        call gf_close(db)
        stop 1
      endif

      call gf_print_stf(stf,tax,6)
      write(*,'(a)') ''

      allocate(seis(db%nstations,GF_NCOMP,tax%nt),tsec(tax%nt),onset(db%nstations),stat=ierr)
      if (ierr /= 0) then
        write(ISTDERR,'(a)') 'Error: could not allocate the seismogram array'
        call gf_locate_release()
        call gf_close(db)
        stop 1
      endif

      ! the partials array exists in every case (zero-size without
      ! --partials), so that the SAC writer has one interface
      call gf_partials_ndp(itypsokern,ndp,ierr)
      allocate(dp(ndp,db%nstations,GF_NCOMP,tax%nt),stat=ierr)
      if (ierr /= 0) then
        write(ISTDERR,'(a)') 'Error: could not allocate the partials array'
        call gf_locate_release()
        call gf_close(db)
        stop 1
      endif

      if (itypsokern > 0) then

        ! seismograms and partials from one pass over the elements
        if (src%source_type /= GF_SRC_CMT) then
          write(ISTDERR,'(a)') 'Error: --partials is defined for a CMTSOLUTION, not a FORCESOLUTION'
          call gf_locate_release()
          call gf_close(db)
          stop 1
        endif

        call gf_seis_cmt_partials(db,src,loc,tax,stf,itypsokern,ndp,seis,dp,tsec,onset,ierr)
        if (ierr /= GF_OK) then
          write(ISTDERR,'(a)') 'Error computing the seismograms and partials'
          write(ISTDERR,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
          call gf_locate_release()
          call gf_close(db)
          stop 1
        endif

      else

        call gf_seis(db,src,loc,tax,stf,seis,tsec,onset,ierr)
        if (ierr /= GF_OK) then
          write(ISTDERR,'(a)') 'Error computing the seismograms'
          write(ISTDERR,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
          call gf_locate_release()
          call gf_close(db)
          stop 1
        endif

      endif

      !--- the writers: every requested format, every product -------------

      if (want_ascii) then
        call gf_write_seis(db,src,loc,tax,stf,seis,tsec,onset,outdir,ierr)
        if (ierr /= GF_OK) then
          write(ISTDERR,'(a)') 'Error writing the seismograms'
          write(ISTDERR,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
          call gf_locate_release()
          call gf_close(db)
          stop 1
        endif
        write(*,'(a,i0,a,a)') 'wrote ',db%nstations,' ASCII seismogram files to ',trim(outdir)

        if (itypsokern > 0) then
          call gf_write_partials(db,src,loc,tax,stf,ndp,dp,tsec,outdir,ierr)
          if (ierr /= GF_OK) then
            write(ISTDERR,'(a)') 'Error writing the partials'
            write(ISTDERR,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
            call gf_locate_release()
            call gf_close(db)
            stop 1
          endif
          write(*,'(a,i0,a,i0,a,a)') 'wrote ',ndp,' ASCII partials per component for ',db%nstations, &
                                     ' stations to ',trim(outdir)
        endif
      endif

      if (want_sac .or. want_sacan) then
        call gf_write_sac(db,src,tax,stf,seis,ndp,dp,outdir,want_sac,want_sacan,ierr)
        if (ierr /= GF_OK) then
          write(ISTDERR,'(a)') 'Error writing the SAC files'
          write(ISTDERR,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
          call gf_locate_release()
          call gf_close(db)
          stop 1
        endif
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

    if (nargs < 2) then
      write(ISTDERR,'(a)') 'Error: --check-anchors needs a database directory'
      call print_usage(ISTDERR)
      stop 1
    endif

    call get_command_argument(2,dbpath)

    call gf_open(dbpath,db,ierr,check_completion=.false.)
    if (ierr /= GF_OK) then
      write(ISTDERR,'(a)') 'Error opening the Green function database'
      write(ISTDERR,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
      stop 1
    endif

    call gf_check_anchors_all(db,worst_err,ielem_worst,ierr)
    if (ierr /= GF_OK) then
      write(ISTDERR,'(a)') 'Error checking the element anchors'
      write(ISTDERR,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
      call gf_close(db)
      stop 1
    endif

    write(*,'(a)')             'anchor consistency check'
    write(*,'(a,i0)')          '  elements checked     = ',db%nelem
    write(*,'(a,es22.14)')     '  worst residual       = ',worst_err
    if (ielem_worst > 0) then
      write(*,'(a,a)')         '  worst element        = ',db%morton_hex(ielem_worst)
    endif
    write(*,'(a,es22.14)')     '  tolerance            = ',GF_ANCHOR_TOL
    write(*,'(a,es22.14)')     '  worst residual, m    = ',worst_err*db%R_PLANET
    if (worst_err > GF_ANCHOR_TOL) then
      write(ISTDERR,'(a)') 'the 27 anchors do not reproduce the stored GLL coordinates'
      write(ISTDERR,'(a)') '  this database was probably written from a USE_GLL = .true. mesh'
      call gf_close(db)
      stop 1
    endif
    write(*,'(a)')             '  result               = ok'

    call gf_close(db)

  case default
    write(ISTDERR,'(a)') 'Error: unknown mode: '//trim(mode)
    call print_usage(ISTDERR)
    stop 1

  end select

  contains

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

! reads a command-line argument as a double, or exits with a message

  implicit none

  integer, intent(in) :: iarg_in
  character(len=*), intent(in) :: name
  double precision, intent(out) :: val

  ! local parameters
  character(len=MAX_STRING_LEN) :: str
  integer :: ios

  call get_command_argument(iarg_in,str)
  read(str,*,iostat=ios) val
  if (ios /= 0) then
    write(ISTDERR,'(a)') 'Error: could not read the '//trim(name)//' from "'//trim(str)//'"'
    stop 1
  endif

  end subroutine read_double_arg

!
!-------------------------------------------------------------------------------------------------
!

  subroutine print_location(db,loc,lat,lon,depth_km,iunit)

! reports a located source
!
! The layout follows gf_print_info: 'key = value', one per line, es22.14 for
! doubles, so the output can be diffed against the solver's own
! OUTPUT_FILES/output_solver.txt -- which is the only oracle this step has.
! Note the solver prints its position through sngl(), so a comparison there
! is bounded by float32 (~3e-8 at these magnitudes), not by our precision.

  use gf_par, only: t_gfdb,t_gf_location

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_location), intent(in) :: loc
  double precision, intent(in) :: lat,lon,depth_km
  integer, intent(in) :: iunit

  ! local parameters
  integer :: i

  write(iunit,'(a)')         'source location'
  write(iunit,'(a,a)')       '  database             = ',trim(db%path)
  write(iunit,'(a,es22.14)') '  latitude             = ',lat
  write(iunit,'(a,es22.14)') '  longitude            = ',lon
  write(iunit,'(a,es22.14)') '  depth, km            = ',depth_km
  write(iunit,'(a)')         ''
  write(iunit,'(a,a)')       '  morton_hex           = ',loc%morton_hex
  write(iunit,'(a,i0,a,i0)') '  element              = ',loc%ielem,' of ',db%nelem
  write(iunit,'(a)')         ''
  write(iunit,'(a,es22.14)') '  xi                   = ',loc%xi
  write(iunit,'(a,es22.14)') '  eta                  = ',loc%eta
  write(iunit,'(a,es22.14)') '  gamma                = ',loc%gamma
  write(iunit,'(a,es22.14)') '  max|xi,eta,gamma|    = ',max(abs(loc%xi),abs(loc%eta),abs(loc%gamma))
  write(iunit,'(a,es22.14)') '  containment tolerance= ',GF_XI_TOL
  write(iunit,'(a)')         ''
  write(iunit,'(a,es22.14)') '  x                    = ',loc%xyz(1)
  write(iunit,'(a,es22.14)') '  y                    = ',loc%xyz(2)
  write(iunit,'(a,es22.14)') '  z                    = ',loc%xyz(3)
  write(iunit,'(a,es22.14)') '  x_target             = ',loc%xyz_target(1)
  write(iunit,'(a,es22.14)') '  y_target             = ',loc%xyz_target(2)
  write(iunit,'(a,es22.14)') '  z_target             = ',loc%xyz_target(3)
  write(iunit,'(a,es22.14)') '  location error, km   = ',loc%distance_km
  write(iunit,'(a)')         ''
  write(iunit,'(a,es22.14)') '  theta                = ',loc%theta
  write(iunit,'(a,es22.14)') '  phi                  = ',loc%phi
  write(iunit,'(a,es22.14)') '  surface radius       = ',loc%r_surface
  write(iunit,'(a)')         ''
  do i = 1,3
    select case (i)
    case (1) ; write(iunit,'(a,3es22.14)') '  nu(N,:)              = ',loc%nu(i,1),loc%nu(i,2),loc%nu(i,3)
    case (2) ; write(iunit,'(a,3es22.14)') '  nu(E,:)              = ',loc%nu(i,1),loc%nu(i,2),loc%nu(i,3)
    case (3) ; write(iunit,'(a,3es22.14)') '  nu(Z,:)              = ',loc%nu(i,1),loc%nu(i,2),loc%nu(i,3)
    end select
  enddo
  write(iunit,'(a)')         ''
  write(iunit,'(a,es22.14)') '  jacobian             = ',loc%jacobian
  write(iunit,'(a,es22.14)') '  anchor residual      = ',loc%anchor_err

  end subroutine print_location

  end program xgf3d
