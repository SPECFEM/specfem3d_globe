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
!---- SAC output: the header rules of the solver's own writer, applied to
!---- the library's traces.
!----
!---- src/specfem3D/write_output_SAC.f90 is welded to specfem_par and cannot
!---- be called from a serial library, so the writing is transcribed --
!---- from GF3DF's copy of it (src/sac/sac.F90), which had already done the
!---- de-welding -- but the header *values* follow the solver, not GF3DF,
!---- wherever the two differ. They differ in two places that matter, and
!---- the forward SAC files in EXAMPLES/ are the oracle for both:
!----
!----   * B: the solver writes seismo_offset*DT - t0 + tshift_src, i.e. -t0
!----     for a single source (tshift_src is zeroed and its original value
!----     returned separately); GF3DF wrote -t_shift, the 29 s that belong
!----     in the origin time, not on the axis. Here B is the first sample of
!----     the planned output axis, tax%t_first, which sits within one stored
!----     sample below the requested -t0 (gf_stf.F90).
!----   * NZSEC/NZMSEC: the solver adds the CMT time shift to the PDE
!----     seconds, with a rollover into the next minute/hour/day/year;
!----     GF3DF's copy has that addition commented out. Those seconds are
!----     the origin-time metadata Stage 4 kept out of the trace
!----     (min_tshift_src_original), and the header is where they go. The
!----     solver's exact expressions are used, float arithmetic and all --
!----     16.40 + 29 gives nzsec 45, nzmsec 399, not 400, and the forward
!----     file says 399.
!----
!---- DELTA is the stored spacing dt*subsample_step and NPTS the output
!---- axis length; neither is comparable with the forward run, which is on
!---- the solver grid. STEL is undefined: the station file stores burial
!---- only (STDP). The event name is the CMTSOLUTION's, read by gf_source.
!----
!---- One writer for every product: a seismogram and a partial derivative
!---- (Stage 6) share the header record and differ in the file name and in
!---- KUSER1/KUSER2, which name the parameter and its unit.
!----
!---- Binary files go through src/shared/binary_c_io.c (the solver's own
!---- writer, one file open at a time), alphanumeric ones through Fortran
!---- formatted I/O. No `use hdf5`, no `use specfem_par`: this is a kernel
!---- module, and tests/gf3d/test_gf_sac.f90 re-reads what it writes without
!---- a database.
!----

  module gf_sac

  use gf_par, only: t_gfdb,t_gf_source,t_gf_taxis,t_gf_stf,gf_set_error, &
                    GF_OK,GF_ERR_ARG,GF_ERR_IO,GF_NCOMP,GF3D_VERSION

  use gf_partials, only: GF_DP_NAME,GF_NDP_LOC

  implicit none

  private

  ! the header fields this library sets; everything else is SAC's
  ! undefined value -12345
  type, public :: t_gf_sac_header
    ! floating point (single precision on disk)
    double precision :: delta = 0.d0, b = 0.d0, o = 0.d0
    double precision :: stla = 0.d0, stlo = 0.d0, stdp = 0.d0
    double precision :: evla = 0.d0, evlo = 0.d0, evdp = 0.d0
    double precision :: user0 = 0.d0, user1 = 0.d0, user2 = 0.d0
    double precision :: cmpaz = 0.d0, cmpinc = 0.d0
    ! integers
    integer :: nzyear = 0, nzjday = 0, nzhour = 0, nzmin = 0, nzsec = 0, nzmsec = 0
    integer :: npts = 0
    ! strings
    character(len=8)  :: kstnm = '', knetwk = '', kcmpnm = '', khole = ''
    character(len=8)  :: kuser0 = '', kuser1 = '', kuser2 = ''
    character(len=16) :: kevnm = ''
  end type t_gf_sac_header

  ! SAC's undefined value, for every field not set above
  real, parameter :: SAC_UNDEF = -12345.0
  character(len=8), parameter :: SAC_UNDEF_STR = '-12345  '

  ! channel names in the stored component order N,E,Z (gf_seismograms)
  character(len=3), dimension(GF_NCOMP), parameter, public :: GF_SAC_CHANNEL = (/ 'BXN','BXE','BXZ' /)

  ! short units for KUSER2 (8 characters)
  character(len=8), dimension(GF_NDP_LOC), parameter :: GF_SAC_DP_UNIT = &
    (/ 'm/dynecm','m/dynecm','m/dynecm','m/dynecm','m/dynecm','m/dynecm', &
       'm/deg   ','m/deg   ','m/km    ','m/s     ' /)

  public :: gf_sac_origin_time
  public :: gf_sac_header
  public :: gf_write_sac_trace
  public :: gf_write_sac

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_sac_origin_time(yr,jda,ho,mi,sec,t_shift,nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec)

! the SAC reference time: the PDE time plus the source's time shift
!
! Transcribed from write_output_SAC.f90:265-300, expression for expression,
! because the forward files were written by it: NZSEC and NZMSEC are formed
! from sec + t_shift in single-statement float arithmetic (16.40 + 29 gives
! 45 and 399), and a sum past 60 s rolls into the minute, hour, day and
! year the way the solver rolls it -- including the day-of-year bookkeeping
! it does through is_leap_year (src/shared/calendar.f90).

  implicit none

  integer, intent(in) :: yr,jda,ho,mi
  double precision, intent(in) :: sec,t_shift
  integer, intent(out) :: nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec

  ! local parameters
  integer :: time_sec
  logical, external :: is_leap_year

  nzyear = yr
  nzjday = jda
  nzhour = ho
  nzmin  = mi

  ! adds time-shift to get the CMT time in the headers as origin time of events
  nzsec  = int(sec+t_shift)
  nzmsec = int((sec+t_shift-int(sec+t_shift))*1000)

  ! Adjust event time and date after t_shift is added
  if (nzsec >= 60) then
    time_sec = jda*24*3600 + ho*3600 + mi*60 + int(sec+t_shift)
    nzjday   = int(time_sec/(24*3600))
    nzhour   = int(mod(time_sec,24*3600)/3600)
    nzmin    = int(mod(time_sec,3600)/60)
    nzsec    = mod(time_sec,60)
    if (nzjday > 365 .and. .not. is_leap_year(nzyear)) then
      nzjday = mod(nzjday,365)
      nzyear = yr + 1
    else if (nzjday > 366 .and. is_leap_year(nzyear)) then
      nzjday = mod(nzjday,366)
      nzyear = yr + 1
    else if (nzjday == 366 .and. is_leap_year(nzyear)) then
      nzjday = 366
    endif
  endif

  end subroutine gf_sac_origin_time

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_sac_header(db,src,tax,stf,ista,icomp,hdr)

! the header record for one station and component, with the solver's rules
!
! KUSER1/KUSER2 are set for a seismogram ('gf3d', the version); a caller
! writing a partial overrides them with the parameter name and its unit.

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(in) :: src
  type(t_gf_taxis), intent(in) :: tax
  type(t_gf_stf), intent(in) :: stf
  integer, intent(in) :: ista,icomp
  type(t_gf_sac_header), intent(out) :: hdr

  hdr = t_gf_sac_header()

  ! the axis: the planned output axis, whose first sample is at or before
  ! the requested -t0 (write_output_SAC.f90:171 has -t0 exactly, for a
  ! single source)
  hdr%delta = tax%dt_sub
  hdr%b     = tax%t_first
  hdr%o     = 0.d0
  hdr%npts  = tax%nt

  ! station values (write_output_SAC.f90:179-182); the station file stores
  ! the burial only, so STEL stays undefined
  hdr%stla = db%stations(ista)%latitude
  hdr%stlo = db%stations(ista)%longitude
  hdr%stdp = db%stations(ista)%depth

  ! event values: the CMT location, as the solver (:187-190)
  hdr%evla = src%latitude
  hdr%evlo = src%longitude
  hdr%evdp = src%depth

  ! USER0 the half duration (:193); USER1/USER2 the shortest and longest
  ! periods the simulation is accurate at (:198-206): the writer set the
  ! database's Gaussian width to T_min/10, so T_min is ten of them, and
  ! the solver's 500 s stays
  hdr%user0 = src%hdur
  hdr%user1 = 10.d0*stf%hdur_db
  hdr%user2 = 500.d0

  ! instrument orientation (:251-258)
  select case (icomp)
  case (1)
    hdr%cmpaz  = 0.d0
    hdr%cmpinc = 90.d0
  case (2)
    hdr%cmpaz  = 90.d0
    hdr%cmpinc = 90.d0
  case default
    hdr%cmpaz  = 0.d0
    hdr%cmpinc = 0.d0
  end select

  ! reference time: the PDE time plus the original time shift (:265-300)
  call gf_sac_origin_time(src%yr,src%jda,src%ho,src%mi,src%sec,src%min_tshift_src_original, &
                          hdr%nzyear,hdr%nzjday,hdr%nzhour,hdr%nzmin,hdr%nzsec,hdr%nzmsec)

  ! names (:325-355). KHOLE marks 3-D synthetics in the solver from the
  ! model name, which the database does not record; 'S3' is its usual
  ! value. KUSER0 is IRIS's network code for synthetics.
  hdr%kstnm  = db%stations(ista)%station(1:min(8,len(db%stations(ista)%station)))
  hdr%knetwk = db%stations(ista)%network(1:min(8,len(db%stations(ista)%network)))
  hdr%kcmpnm = GF_SAC_CHANNEL(icomp)
  hdr%kevnm  = src%event_name
  hdr%khole  = 'S3'
  hdr%kuser0 = 'SY'
  hdr%kuser1 = 'gf3d'
  hdr%kuser2 = GF3D_VERSION

  end subroutine gf_sac_header

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_write_sac_trace(hdr,trace,nt,basename,binary,alphanum,ierr)

! writes one trace as <basename>.sac (binary) and/or <basename>.sacan
! (alphanumeric), in single precision as the solver does
!
! The field order is the SAC file format's (632-byte header: 70 reals, 40
! integers, 24 strings), transcribed from GF3DF's sac.F90:497-640, itself a
! transcription of write_output_SAC.f90.

  implicit none

  type(t_gf_sac_header), intent(in) :: hdr
  integer, intent(in) :: nt
  double precision, dimension(nt), intent(in) :: trace
  character(len=*), intent(in) :: basename
  logical, intent(in) :: binary,alphanum
  integer, intent(out) :: ierr

  ! local parameters
  real, dimension(:), allocatable :: tmp
  real :: undef,bysac,internal,unused,scale_f
  real :: delta,b,e,o,a,stla,stlo,stel,stdp,evla,evlo,evel,evdp,mag
  real :: user0,user1,user2,dist,az,baz,gcarc,depmin,depmax,depmen,cmpaz,cmpinc,odelta
  integer :: nvhdr,norid,nevid,iftype,idep,iztype,ievtyp,iqual,isynth,imagtyp
  integer :: leven,lpspol,lovrok,lcalda
  integer :: iout,ios,isample,imodulo_5,ier

  ierr = GF_OK

  if (nt < 1) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_write_sac_trace: an empty trace')
    return
  endif
  if (hdr%npts /= nt) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_write_sac_trace: header npts does not match the trace')
    return
  endif
  if (.not. binary .and. .not. alphanum) return

  ! the solver's defaults (write_output_SAC.f90:150-172, 300-321)
  undef    = SAC_UNDEF
  bysac    = SAC_UNDEF
  internal = SAC_UNDEF
  unused   = SAC_UNDEF
  scale_f  = 1.0e9          ! factor to nm, as the solver writes it

  delta  = real(hdr%delta)
  b      = real(hdr%b)
  e      = bysac
  o      = real(hdr%o)
  a      = undef
  odelta = undef
  depmin = bysac
  depmax = bysac
  depmen = bysac
  stla = real(hdr%stla) ; stlo = real(hdr%stlo) ; stel = undef ; stdp = real(hdr%stdp)
  evla = real(hdr%evla) ; evlo = real(hdr%evlo) ; evel = undef ; evdp = real(hdr%evdp)
  mag  = undef
  user0 = real(hdr%user0) ; user1 = real(hdr%user1) ; user2 = real(hdr%user2)
  dist = bysac ; az = bysac ; baz = bysac ; gcarc = bysac
  cmpaz = real(hdr%cmpaz) ; cmpinc = real(hdr%cmpinc)

  nvhdr   = 6
  norid   = int(undef)
  nevid   = int(undef)
  iftype  = 1        ! ITIME: a time series
  idep    = 6        ! displacement in nm
  iztype  = 11       ! IO: reference time is the origin time
  ievtyp  = 40       ! earthquake
  iqual   = int(undef)
  isynth  = int(undef)
  imagtyp = int(undef)
  leven   = 1
  lpspol  = 1
  lovrok  = 1
  lcalda  = 1

  !--- alphanumeric ----------------------------------------------------------

  if (alphanum) then

    open(newunit=iout,file=trim(basename)//'.sacan',status='replace',action='write',iostat=ios)
    if (ios /= 0) then
      call gf_set_error(ierr,GF_ERR_IO,'could not open for writing: '//trim(basename)//'.sacan')
      return
    endif

    ! the 14 real cards, 5 per card, in file order
    write(iout,510) delta,    depmin,  depmax,  scale_f, odelta
    write(iout,510) b,        e,       o,       a,       internal
    write(iout,510) undef,    undef,   undef,   undef,   undef
    write(iout,510) undef,    undef,   undef,   undef,   undef
    write(iout,510) undef,    undef,   undef,   undef,   undef
    write(iout,510) undef,    undef,   undef,   undef,   undef
    write(iout,510) undef,    stla,    stlo,    stel,    stdp
    write(iout,510) evla,     evlo,    evel,    evdp,    mag
    write(iout,510) user0,    user1,   user2,   undef,   undef
    write(iout,510) undef,    undef,   undef,   undef,   undef
    write(iout,510) dist,     az,      baz,     gcarc,   internal
    write(iout,510) internal, depmen,  cmpaz,   cmpinc,  undef
    write(iout,510) undef,    undef,   undef,   undef,   undef
    write(iout,510) unused,   unused,  unused,  unused,  unused
    ! the 8 integer cards
    write(iout,520) hdr%nzyear, hdr%nzjday, hdr%nzhour, hdr%nzmin, hdr%nzsec
    write(iout,520) hdr%nzmsec, nvhdr, norid, nevid, hdr%npts
    write(iout,520) int(undef), int(undef), int(undef), int(undef), int(undef)
    write(iout,520) iftype, idep, iztype, int(unused), int(undef)
    write(iout,520) int(undef), int(undef), ievtyp, int(undef), isynth
    write(iout,520) imagtyp, int(undef), int(undef), int(undef), int(undef)
    write(iout,520) int(unused), int(unused), int(unused), int(unused), int(unused)
    write(iout,520) leven, lpspol, lovrok, lcalda, int(unused)
    ! the 8 string cards
    write(iout,530) hdr%kstnm, hdr%kevnm
    write(iout,540) hdr%khole, SAC_UNDEF_STR, SAC_UNDEF_STR
    write(iout,540) SAC_UNDEF_STR, SAC_UNDEF_STR, SAC_UNDEF_STR
    write(iout,540) SAC_UNDEF_STR, SAC_UNDEF_STR, SAC_UNDEF_STR
    write(iout,540) SAC_UNDEF_STR, SAC_UNDEF_STR, SAC_UNDEF_STR
    write(iout,540) SAC_UNDEF_STR, SAC_UNDEF_STR, hdr%kuser0
    write(iout,540) hdr%kuser1, hdr%kuser2, hdr%kcmpnm
    write(iout,540) hdr%knetwk, SAC_UNDEF_STR, SAC_UNDEF_STR

    ! the data, five per line
    imodulo_5 = mod(nt,5)
    do isample = 1,nt-imodulo_5,5
      write(iout,510) real(trace(isample)),real(trace(isample+1)),real(trace(isample+2)), &
                      real(trace(isample+3)),real(trace(isample+4))
    enddo
    if (imodulo_5 > 0) then
      write(iout,510) (real(trace(isample)),isample = nt-imodulo_5+1,nt)
    endif

    close(iout)

  endif

510 format(5G15.7)
520 format(5I10)
530 format(A8,A16)
540 format(A8,A8,A8)

  !--- binary ----------------------------------------------------------------

  if (binary) then

    allocate(tmp(nt),stat=ier)
    if (ier /= 0) then
      call gf_set_error(ierr,GF_ERR_IO,'gf_write_sac_trace: could not allocate the output buffer')
      return
    endif

    call open_file_create(trim(basename)//'.sac'//char(0))

    ! reals 1:70
    call write_real(delta)         !(1)
    call write_real(depmin)        !(2)
    call write_real(depmax)        !(3)
    call write_real(scale_f)       !(4)
    call write_real(odelta)        !(5)
    call write_real(b)             !(6)
    call write_real(e)             !(7)
    call write_real(o)             !(8)
    call write_real(a)             !(9)
    call write_real(internal)      !(10)
    do isample = 11,31             !(11:31) T0..T9, F, RESP0..RESP9
      call write_real(undef)
    enddo
    call write_real(stla)          !(32)
    call write_real(stlo)          !(33)
    call write_real(stel)          !(34)
    call write_real(stdp)          !(35)
    call write_real(evla)          !(36)
    call write_real(evlo)          !(37)
    call write_real(evel)          !(38)
    call write_real(evdp)          !(39)
    call write_real(mag)           !(40)
    call write_real(user0)         !(41)
    call write_real(user1)         !(42)
    call write_real(user2)         !(43)
    do isample = 44,50             !(44:50) USER3..USER9
      call write_real(undef)
    enddo
    call write_real(dist)          !(51)
    call write_real(az)            !(52)
    call write_real(baz)           !(53)
    call write_real(gcarc)         !(54)
    call write_real(internal)      !(55)
    call write_real(internal)      !(56)
    call write_real(depmen)        !(57)
    call write_real(cmpaz)         !(58)
    call write_real(cmpinc)        !(59)
    do isample = 60,70             !(60:70) XMINIMUM..UNUSED
      call write_real(undef)
    enddo

    ! integers 71:110
    call write_integer(hdr%nzyear)    !(71)
    call write_integer(hdr%nzjday)    !(72)
    call write_integer(hdr%nzhour)    !(73)
    call write_integer(hdr%nzmin)     !(74)
    call write_integer(hdr%nzsec)     !(75)
    call write_integer(hdr%nzmsec)    !(76)
    call write_integer(nvhdr)         !(77)
    call write_integer(norid)         !(78)
    call write_integer(nevid)         !(79)
    call write_integer(hdr%npts)      !(80)
    do isample = 81,85                !(81:85) UNUSED, NWFID, NXSIZE, NYSIZE, UNUSED
      call write_integer(int(undef))
    enddo
    call write_integer(iftype)        !(86)
    call write_integer(idep)          !(87)
    call write_integer(iztype)        !(88)
    call write_integer(int(undef))    !(89)
    call write_integer(int(undef))    !(90) IINST
    call write_integer(int(undef))    !(91) ISTREG
    call write_integer(int(undef))    !(92) IEVREG
    call write_integer(ievtyp)        !(93)
    call write_integer(iqual)         !(94)
    call write_integer(isynth)        !(95)
    call write_integer(imagtyp)       !(96)
    call write_integer(int(undef))    !(97) IMAGSRC
    do isample = 98,105               !(98:105) UNUSED
      call write_integer(int(unused))
    enddo
    call write_integer(leven)         !(106)
    call write_integer(lpspol)        !(107)
    call write_integer(lovrok)        !(108)
    call write_integer(lcalda)        !(109)
    call write_integer(int(unused))   !(110)

    ! strings 111:302
    call write_character(hdr%kstnm,8)        !(111:118)
    call write_character(hdr%kevnm,16)       !(119:134)
    call write_character(hdr%khole,8)        !(135:142)
    do isample = 1,13                        !(143:246) KO, KA, KT0..KT9, KF
      call write_character(SAC_UNDEF_STR,8)
    enddo
    call write_character(hdr%kuser0,8)       !(247:254)
    call write_character(hdr%kuser1,8)       !(255:262)
    call write_character(hdr%kuser2,8)       !(263:270)
    call write_character(hdr%kcmpnm,8)       !(271:278)
    call write_character(hdr%knetwk,8)       !(279:286)
    call write_character(SAC_UNDEF_STR,8)    !(287:294) KDATRD
    call write_character(SAC_UNDEF_STR,8)    !(295:302) KINST

    ! the data
    tmp(1:nt) = real(trace(1:nt))
    call write_n_real(tmp,nt)

    call close_file()

    deallocate(tmp)

  endif

  end subroutine gf_write_sac_trace

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_write_sac(db,src,tax,stf,seis,ndp,dp,outdir,binary,alphanum,ierr)

! every station and component as SAC: the seismograms, and the partials
! when ndp > 0
!
! File names follow the solver: <outdir>/NET.STA.BXN.sem.sac (and .sacan);
! a partial is <outdir>/NET.STA.BXN.<name>.sem.sac with <name> from
! GF_DP_NAME, and its KUSER1/KUSER2 carry the name and the unit.

  implicit none

  type(t_gfdb), intent(in) :: db
  type(t_gf_source), intent(in) :: src
  type(t_gf_taxis), intent(in) :: tax
  type(t_gf_stf), intent(in) :: stf
  double precision, dimension(db%nstations,GF_NCOMP,tax%nt), intent(in) :: seis
  integer, intent(in) :: ndp
  double precision, dimension(ndp,db%nstations,GF_NCOMP,tax%nt), intent(in) :: dp
  character(len=*), intent(in) :: outdir
  logical, intent(in) :: binary,alphanum
  integer, intent(out) :: ierr

  ! local parameters
  type(t_gf_sac_header) :: hdr
  double precision, dimension(:), allocatable :: trace
  character(len=len(outdir)+64) :: basename
  integer :: ista,icomp,ip,it,ier

  ierr = GF_OK

  if (ndp < 0 .or. ndp > GF_NDP_LOC) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_write_sac: ndp must be between 0 and 10')
    return
  endif

  allocate(trace(tax%nt),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_IO,'gf_write_sac: could not allocate the trace buffer')
    return
  endif

  do ista = 1,db%nstations
    do icomp = 1,GF_NCOMP

      call gf_sac_header(db,src,tax,stf,ista,icomp,hdr)

      basename = trim(outdir)//'/'//trim(db%stations(ista)%id)//'.'//GF_SAC_CHANNEL(icomp)//'.sem'
      do it = 1,tax%nt
        trace(it) = seis(ista,icomp,it)
      enddo
      call gf_write_sac_trace(hdr,trace,tax%nt,basename,binary,alphanum,ierr)
      if (ierr /= GF_OK) goto 99

      do ip = 1,ndp
        hdr%kuser1 = GF_DP_NAME(ip)
        hdr%kuser2 = GF_SAC_DP_UNIT(ip)
        basename = trim(outdir)//'/'//trim(db%stations(ista)%id)//'.'//GF_SAC_CHANNEL(icomp)// &
                   '.'//GF_DP_NAME(ip)//'.sem'
        do it = 1,tax%nt
          trace(it) = dp(ip,ista,icomp,it)
        enddo
        call gf_write_sac_trace(hdr,trace,tax%nt,basename,binary,alphanum,ierr)
        if (ierr /= GF_OK) goto 99
      enddo

    enddo
  enddo

99 continue
  if (allocated(trace)) deallocate(trace)

  end subroutine gf_write_sac

  end module gf_sac
