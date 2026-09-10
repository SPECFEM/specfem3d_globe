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
!---- test_gf_sac -- src/gf3d/gf_sac.F90
!----
!---- Tier 1: no database, no HDF5, no MPI. A station record, a source and
!---- the regional example's plan are built by hand, one seismogram and one
!---- partial are written as SAC, and the files are read back with stream
!---- I/O against the SAC format's own field layout (632-byte header: 70
!---- reals, 40 integers, 24 strings). What is pinned:
!----
!----   * the header rules taken from the solver's writer: B is the first
!----     sample of the planned axis, DELTA the stored spacing, O = 0, the
!----     reference time is the PDE time plus the source's time shift with
!----     the solver's own rollover -- 16.40 + 29 s giving nzsec 45 and
!----     nzmsec 399, as the forward file in EXAMPLES/ has it -- and the
!----     minute, day and year boundaries the same way;
!----   * the header reals and the data round-trip to single precision
!----     (1e-6 relative; the trace peak as the data's scale); integers
!----     and strings exactly;
!----   * a partial's file name and KUSER1/KUSER2 name the parameter.
!----
!---- The oracle is the file format and the solver's expressions, not
!---- another writer. Every member of a grouped assertion names itself in
!---- the log when it fails: on the CI the log is all there is of a run.
!----

  program test_gf_sac

  use gf_par, only: t_gfdb,t_gf_source,t_gf_stf,t_gf_taxis,t_gf_station, &
                    GF_OK,GF_NCOMP,GF_SRC_CMT,GF_STF_TRUNC,GF3D_VERSION

  use gf_stf, only: gf_stf_plan,gf_taxis_plan

  use gf_partials, only: GF_NDP_MT

  use gf_sac

  use gf_manufactured

  implicit none

  ! the shipped regional example's widths and grid
  double precision, parameter :: HDB  = 6.94968291528492d0
  double precision, parameter :: HCMT = 60.d0
  double precision, parameter :: T0DB = 34.74841457642461d0
  double precision, parameter :: DTR  = 3.4d0
  integer, parameter :: NTDB = 544, SS = 34

  character(len=*), parameter :: OUTDIR = './OUTPUT_FILES'

  ! single-precision storage: unit roundoff 2^-24 = 6e-8 relative, and the
  ! alphanumeric cards' G15.7 carry 7 digits; 1e-6 covers both with headroom.
  ! Not equality: the test's real(x) and the writer's are two compilation
  ! units (the CI ifort failed the bitwise form), and a real defect -- a
  ! shifted sample, a swapped component, a wrong slot -- is a whole sample.
  double precision, parameter :: TOL_SINGLE = 1.d-6

  integer :: nfail

  nfail = 0

  write(*,'(a)') 'test_gf_sac: SAC output'
  write(*,'(a)') ''

  call test_origin_time(nfail)
  call test_files(nfail)

  write(*,'(a)') ''
  if (nfail > 0) then
    write(*,'(a,i0,a)') 'test_gf_sac: ',nfail,' assertion(s) FAILED'
    write(0,'(a,i0,a)') 'test_gf_sac: ',nfail,' assertion(s) FAILED'
    stop 1
  endif
  write(*,'(a)') 'test_gf_sac: all assertions passed'

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_origin_time(nfail)

! the solver's reference-time arithmetic, including the rollovers

  implicit none
  integer, intent(inout) :: nfail

  integer :: yr,jda,ho,mi,sec,msec

  write(*,'(a)') '1. reference time'

  ! the shipped CMTSOLUTION: PDE 1994 160 0:33:16.40, time shift 29 s; the
  ! forward SAC file has nzsec 45 and nzmsec 399 -- the 399 is the float
  ! arithmetic of the solver's expression, and it must be reproduced
  call gf_sac_origin_time(1994,160,0,33,16.40d0,29.d0,yr,jda,ho,mi,sec,msec)
  call gf_report_true('060994A + 29 s: 1994 160 00:33:45.399', &
                      yr == 1994 .and. jda == 160 .and. ho == 0 .and. mi == 33 .and. &
                      sec == 45 .and. msec == 399,nfail)

  ! into the next minute
  call gf_sac_origin_time(1994,160,0,33,50.d0,20.d0,yr,jda,ho,mi,sec,msec)
  call gf_report_true('rollover into the next minute',ho == 0 .and. mi == 34 .and. sec == 10 .and. msec == 0,nfail)

  ! into the next day
  call gf_sac_origin_time(1994,160,23,59,50.d0,20.d0,yr,jda,ho,mi,sec,msec)
  call gf_report_true('rollover into the next day', &
                      yr == 1994 .and. jda == 161 .and. ho == 0 .and. mi == 0 .and. sec == 10,nfail)

  ! into the next year, common and leap
  call gf_sac_origin_time(1994,365,23,59,50.d0,20.d0,yr,jda,ho,mi,sec,msec)
  call gf_report_true('rollover into the next year (1994 -> 1995 day 1)',yr == 1995 .and. jda == 1,nfail)
  call gf_sac_origin_time(1996,366,23,59,50.d0,20.d0,yr,jda,ho,mi,sec,msec)
  call gf_report_true('rollover into the next year (leap 1996 -> 1997 day 1)',yr == 1997 .and. jda == 1,nfail)
  call gf_sac_origin_time(1996,365,23,59,50.d0,20.d0,yr,jda,ho,mi,sec,msec)
  call gf_report_true('day 365 -> 366 of a leap year stays in it',yr == 1996 .and. jda == 366,nfail)

  ! no shift, no change
  call gf_sac_origin_time(2001,42,7,8,9.25d0,0.d0,yr,jda,ho,mi,sec,msec)
  call gf_report_true('zero shift: the PDE time itself', &
                      yr == 2001 .and. jda == 42 .and. ho == 7 .and. mi == 8 .and. sec == 9 .and. msec == 250,nfail)

  end subroutine test_origin_time

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_files(nfail)

! write, then read back with stream I/O against the format

  implicit none
  integer, intent(inout) :: nfail

  type(t_gfdb) :: db
  type(t_gf_source) :: src
  type(t_gf_stf) :: stf
  type(t_gf_taxis) :: tax
  type(t_gf_sac_header) :: hdr
  double precision, dimension(:,:,:), allocatable :: seis
  double precision, dimension(:,:,:,:), allocatable :: dp
  real, dimension(70) :: fh
  integer, dimension(40) :: ih
  character(len=192) :: ch
  real, dimension(:), allocatable :: data
  real, dimension(5) :: card
  integer, dimension(5) :: icard
  character(len=8) :: c8a,c8b,c8c
  character(len=16) :: c16
  double precision :: t,err_real,err_data,peak,e,emax
  integer :: ierr,nt,it,icomp,iunit,ios,nbad,nlines
  integer :: nbad_file,nbad_int,nbad_str,itmax
  integer(kind=8) :: fsize
  logical :: exists,readok
  character(len=256) :: fname
  character(len=3), dimension(GF_NCOMP) :: chn
  character(len=5) :: pre

  write(*,'(a)') '2. files'

  !--- a station, a source, the regional plan -----------------------------

  db%nstations = 1
  allocate(db%stations(1))
  db%stations(1)%network = 'IU'
  db%stations(1)%station = 'SJG'
  db%stations(1)%id = 'IU.SJG'
  db%stations(1)%latitude = 18.1091d0
  db%stations(1)%longitude = -66.15d0
  db%stations(1)%depth = 0.d0
  db%stations(1)%hdur = HDB

  src%source_type = GF_SRC_CMT
  src%latitude = -5.812d0
  src%longitude = -75.27d0
  src%depth = 122.6d0
  src%hdur = HCMT
  src%min_tshift_src_original = 29.d0
  src%yr = 1994 ; src%jda = 160 ; src%ho = 0 ; src%mi = 33 ; src%sec = 16.40d0
  src%event_name = '060994A'

  call gf_stf_plan(GF_SRC_CMT,0,HCMT,HDB,DTR,GF_STF_TRUNC,stf,ierr)
  call gf_taxis_plan(NTDB,0.1d0,SS,T0DB,90.d0,tax,ierr)
  nt = tax%nt
  call gf_report_true('regional CMT plan and axis (nt = 562)',ierr == GF_OK .and. nt == 562,nfail)

  ! a distinct smooth trace per component, and one partial
  allocate(seis(1,GF_NCOMP,nt),dp(1,1,GF_NCOMP,nt))
  do it = 1,nt
    t = tax%t_first + dble(it-1)*tax%dt_sub
    do icomp = 1,GF_NCOMP
      seis(1,icomp,it) = 1.d-5*dble(icomp)*sin(t/(40.d0*icomp))*exp(-((t - 700.d0)/300.d0)**2)
      dp(1,1,icomp,it) = 1.d-33*cos(t/(25.d0 + icomp))
    enddo
  enddo

  !--- the header record ----------------------------------------------------

  call gf_sac_header(db,src,tax,stf,1,2,hdr)
  call gf_report_true('header: B = t_first, DELTA = dt_sub, NPTS = nt, O = 0', &
                      hdr%b == tax%t_first .and. hdr%delta == tax%dt_sub .and. &
                      hdr%npts == nt .and. hdr%o == 0.d0,nfail)
  call gf_report_true('header: E component is cmpaz 90, cmpinc 90, BXE', &
                      hdr%cmpaz == 90.d0 .and. hdr%cmpinc == 90.d0 .and. hdr%kcmpnm == 'BXE',nfail)
  call gf_report_true('header: reference time 1994 160 00:33:45.399', &
                      hdr%nzyear == 1994 .and. hdr%nzjday == 160 .and. hdr%nzsec == 45 .and. hdr%nzmsec == 399,nfail)
  call gf_report_true('header: names', &
                      hdr%kstnm == 'SJG' .and. hdr%knetwk == 'IU' .and. hdr%kevnm == '060994A' .and. &
                      hdr%khole == 'S3' .and. hdr%kuser0 == 'SY' .and. hdr%kuser1 == 'gf3d' .and. &
                      hdr%kuser2 == GF3D_VERSION,nfail)
  call gf_report('header: USER1 = T_min = 10 hdur_db',abs(hdr%user1 - 10.d0*HDB),1.d-12,nfail)

  !--- write everything ------------------------------------------------------

  call gf_write_sac(db,src,tax,stf,seis,1,dp,OUTDIR,.true.,.true.,ierr)
  call gf_report_true('gf_write_sac returns GF_OK',ierr == GF_OK,nfail)

  chn(1) = 'BXN' ; chn(2) = 'BXE' ; chn(3) = 'BXZ'

  !--- the binary seismograms, read back against the format ---------------
  !
  ! Five groups over the three files: the bytes, the header reals, the
  ! header integers, the header strings, the data. The reals and the data
  ! assert TOL_SINGLE against the doubles the writer was given -- not
  ! equality with a real(x) re-evaluated here, which is a claim across two
  ! compilation units and the one the CI ifort failed; integers and strings
  ! are bytes and stay exact. A file that cannot be read fails the first
  ! group and drops out of the others.

  allocate(data(nt))
  nbad_file = 0 ; nbad_int = 0 ; nbad_str = 0
  err_real = 0.d0 ; err_data = 0.d0
  do icomp = 1,GF_NCOMP
    pre = chn(icomp)//': '
    fname = OUTDIR//'/IU.SJG.'//chn(icomp)//'.sem.sac'
    inquire(file=trim(fname),exist=exists)
    call note(exists,pre//'missing '//trim(fname),nbad_file)
    if (.not. exists) cycle
    ! the size first, so that a short file and a failed read are told apart
    inquire(file=trim(fname),size=fsize)
    call note(fsize == 632 + 4*nt,pre//'file size',nbad_file,int(fsize))
    open(newunit=iunit,file=trim(fname),access='stream',form='unformatted',status='old',action='read',iostat=ios)
    call note(ios == 0,pre//'open, ios',nbad_file,ios)
    if (ios /= 0) cycle
    read(iunit,iostat=ios) fh,ih,ch,data
    close(iunit)
    call note(ios == 0,pre//'stream read, ios',nbad_file,ios)
    if (ios /= 0) cycle

    ! reals: DELTA(1) B(6) O(8) STLA(32) STLO(33) STEL(34) STDP(35) EVLA(36)
    ! EVLO(37) EVDP(39) USER0(41) USER2(43) CMPAZ(58) CMPINC(59)
    call note_real(pre//'DELTA',fh(1),tax%dt_sub,err_real)
    call note_real(pre//'B',fh(6),tax%t_first,err_real)
    call note_real(pre//'O',fh(8),0.d0,err_real)
    call note_real(pre//'STLA',fh(32),18.1091d0,err_real)
    call note_real(pre//'STLO',fh(33),-66.15d0,err_real)
    call note_real(pre//'STEL',fh(34),-12345.d0,err_real)
    call note_real(pre//'STDP',fh(35),0.d0,err_real)
    call note_real(pre//'EVLA',fh(36),-5.812d0,err_real)
    call note_real(pre//'EVLO',fh(37),-75.27d0,err_real)
    call note_real(pre//'EVDP',fh(39),122.6d0,err_real)
    call note_real(pre//'USER0',fh(41),60.d0,err_real)
    call note_real(pre//'USER2',fh(43),500.d0,err_real)
    select case (icomp)
    case (1)
      call note_real(pre//'CMPAZ',fh(58),0.d0,err_real)
      call note_real(pre//'CMPINC',fh(59),90.d0,err_real)
    case (2)
      call note_real(pre//'CMPAZ',fh(58),90.d0,err_real)
      call note_real(pre//'CMPINC',fh(59),90.d0,err_real)
    case (3)
      call note_real(pre//'CMPAZ',fh(58),0.d0,err_real)
      call note_real(pre//'CMPINC',fh(59),0.d0,err_real)
    end select
    ! integers: NZYEAR..NZMSEC (1..6), NVHDR(7), NPTS(10), IFTYPE(16),
    ! IDEP(17), IZTYPE(18), LEVEN(36), LCALDA(39)
    call note(ih(1) == 1994,pre//'NZYEAR',nbad_int,ih(1))
    call note(ih(2) == 160,pre//'NZJDAY',nbad_int,ih(2))
    call note(ih(3) == 0,pre//'NZHOUR',nbad_int,ih(3))
    call note(ih(4) == 33,pre//'NZMIN',nbad_int,ih(4))
    call note(ih(5) == 45,pre//'NZSEC',nbad_int,ih(5))
    call note(ih(6) == 399,pre//'NZMSEC',nbad_int,ih(6))
    call note(ih(7) == 6,pre//'NVHDR',nbad_int,ih(7))
    call note(ih(10) == nt,pre//'NPTS',nbad_int,ih(10))
    call note(ih(16) == 1,pre//'IFTYPE',nbad_int,ih(16))
    call note(ih(17) == 6,pre//'IDEP',nbad_int,ih(17))
    call note(ih(18) == 11,pre//'IZTYPE',nbad_int,ih(18))
    call note(ih(36) == 1,pre//'LEVEN',nbad_int,ih(36))
    call note(ih(39) == 1,pre//'LCALDA',nbad_int,ih(39))
    ! strings: KSTNM(1:8) KEVNM(9:24) KHOLE(25:32) KUSER0(137:144) KUSER1(145:152)
    ! KUSER2(153:160) KCMPNM(161:168) KNETWK(169:176); blank-padded, as the
    ! writer pads them
    call note(ch(1:8) == 'SJG',pre//'KSTNM '//ch(1:8),nbad_str)
    call note(ch(9:24) == '060994A',pre//'KEVNM '//ch(9:24),nbad_str)
    call note(ch(25:32) == 'S3',pre//'KHOLE '//ch(25:32),nbad_str)
    call note(ch(137:144) == 'SY',pre//'KUSER0 '//ch(137:144),nbad_str)
    call note(ch(145:152) == 'gf3d',pre//'KUSER1 '//ch(145:152),nbad_str)
    call note(ch(153:160) == GF3D_VERSION,pre//'KUSER2 '//ch(153:160),nbad_str)
    call note(ch(161:168) == chn(icomp),pre//'KCMPNM '//ch(161:168),nbad_str)
    call note(ch(169:176) == 'IU',pre//'KNETWK '//ch(169:176),nbad_str)
    ! the data, against the doubles the writer was given, relative to the
    ! trace's peak: per sample, a zero crossing would divide by ~0 and a
    ! tiny sample flushed under -ftz would count as a defect; the rounding
    ! bound u|x| <= u*peak still holds, and a shifted sample or a swapped
    ! component is a whole sample of the peak, not a rounding
    peak = maxval(abs(seis(1,icomp,:)))
    emax = 0.d0 ; itmax = 0
    do it = 1,nt
      e = abs(dble(data(it)) - seis(1,icomp,it))/peak
      if (e > emax) then
        emax = e ; itmax = it
      endif
    enddo
    if (emax > TOL_SINGLE) write(*,'(a,es10.3,a,i0)') '       mismatch: '//pre//'data, rel. err ',emax,' at sample ',itmax
    err_data = max(err_data,emax)
  enddo
  call gf_report_true('binary seismograms: 3 files of 632 + 4 nt bytes, read back',nbad_file == 0,nfail)
  call gf_report('binary seismograms: header reals to single precision (rel)',err_real,TOL_SINGLE,nfail)
  call gf_report_true('binary seismograms: header integers NZ*, NVHDR, NPTS, IFTYPE, IDEP, IZTYPE, LEVEN, LCALDA', &
                      nbad_int == 0,nfail)
  call gf_report_true('binary seismograms: header strings KSTNM, KEVNM, KHOLE, KUSER0..2, KCMPNM, KNETWK', &
                      nbad_str == 0,nfail)
  call gf_report('binary seismograms: data to single precision (rel. to the trace peak)',err_data,TOL_SINGLE,nfail)

  !--- the binary partial ----------------------------------------------------

  nbad_file = 0 ; nbad_str = 0
  err_real = 0.d0 ; err_data = 0.d0
  readok = .false.
  fname = OUTDIR//'/IU.SJG.BXZ.Mrr.sem.sac'
  inquire(file=trim(fname),exist=exists)
  call note(exists,'missing '//trim(fname),nbad_file)
  if (exists) then
    inquire(file=trim(fname),size=fsize)
    call note(fsize == 632 + 4*nt,'file size',nbad_file,int(fsize))
    open(newunit=iunit,file=trim(fname),access='stream',form='unformatted',status='old',action='read',iostat=ios)
    call note(ios == 0,'open, ios',nbad_file,ios)
    if (ios == 0) then
      read(iunit,iostat=ios) fh,ih,ch,data
      close(iunit)
      call note(ios == 0,'stream read, ios',nbad_file,ios)
      readok = (ios == 0)
    endif
  endif
  call gf_report_true('binary partial: file of 632 + 4 nt bytes, read back',nbad_file == 0,nfail)
  if (readok) then
    call note(ch(145:152) == 'Mrr','KUSER1 '//ch(145:152),nbad_str)
    call note(ch(153:160) == 'm/dynecm','KUSER2 '//ch(153:160),nbad_str)
    call note(ch(161:168) == 'BXZ','KCMPNM '//ch(161:168),nbad_str)
    call note(ih(10) == nt,'NPTS',nbad_str,ih(10))
    call gf_report_true('binary partial: KUSER1 = Mrr, KUSER2 = m/dynecm, KCMPNM = BXZ, NPTS = nt',nbad_str == 0,nfail)
    call note_real('B',fh(6),tax%t_first,err_real)
    call gf_report('binary partial: B = t_first to single precision (rel)',err_real,TOL_SINGLE,nfail)
    peak = maxval(abs(dp(1,1,3,:)))
    emax = 0.d0 ; itmax = 0
    do it = 1,nt
      e = abs(dble(data(it)) - dp(1,1,3,it))/peak
      if (e > emax) then
        emax = e ; itmax = it
      endif
    enddo
    if (emax > TOL_SINGLE) write(*,'(a,es10.3,a,i0)') '       mismatch: data, rel. err ',emax,' at sample ',itmax
    err_data = emax
    call gf_report('binary partial: data to single precision (rel. to the trace peak)',err_data,TOL_SINGLE,nfail)
  else
    write(*,'(a)') '       skipped: the partial''s header and data comparisons (the file was not read)'
  endif

  !--- the alphanumeric file -------------------------------------------------

  nbad = 0
  fname = OUTDIR//'/IU.SJG.BXN.sem.sacan'
  inquire(file=trim(fname),exist=exists)
  if (exists) then
    open(newunit=iunit,file=trim(fname),status='old',action='read',iostat=ios)
    read(iunit,'(5G15.7)',iostat=ios) card            ! DELTA DEPMIN DEPMAX SCALE ODELTA
    if (abs(card(1) - real(tax%dt_sub)) > 1.e-6*real(tax%dt_sub)) nbad = nbad + 1
    read(iunit,'(5G15.7)',iostat=ios) card            ! B E O A INTERNAL
    if (abs(card(1) - real(tax%t_first)) > 1.e-6*abs(real(tax%t_first)) .or. card(3) /= 0.0) nbad = nbad + 1
    do it = 3,14
      read(iunit,'(5G15.7)',iostat=ios) card
    enddo
    read(iunit,'(5I10)',iostat=ios) icard             ! NZYEAR NZJDAY NZHOUR NZMIN NZSEC
    if (icard(1) /= 1994 .or. icard(2) /= 160 .or. icard(5) /= 45) nbad = nbad + 1
    read(iunit,'(5I10)',iostat=ios) icard             ! NZMSEC NVHDR NORID NEVID NPTS
    if (icard(1) /= 399 .or. icard(2) /= 6 .or. icard(5) /= nt) nbad = nbad + 1
    do it = 3,8
      read(iunit,'(5I10)',iostat=ios) icard
    enddo
    read(iunit,'(A8,A16)',iostat=ios) c8a,c16          ! KSTNM KEVNM
    if (c8a /= 'SJG' .or. c16 /= '060994A') nbad = nbad + 1
    read(iunit,'(A8,A8,A8)',iostat=ios) c8a,c8b,c8c    ! card 2: KHOLE KO KA
    if (c8a /= 'S3') nbad = nbad + 1
    do it = 3,6                                        ! cards 3-5: KT0..KT8; card 6: KT9 KF KUSER0
      read(iunit,'(A8,A8,A8)',iostat=ios) c8a,c8b,c8c
    enddo
    if (c8c /= 'SY') nbad = nbad + 1
    read(iunit,'(A8,A8,A8)',iostat=ios) c8a,c8b,c8c    ! card 7: KUSER1 KUSER2 KCMPNM
    if (c8a /= 'gf3d' .or. c8c /= 'BXN') nbad = nbad + 1
    read(iunit,'(A8,A8,A8)',iostat=ios) c8a,c8b,c8c    ! card 8: KNETWK KDATRD KINST
    if (c8a /= 'IU') nbad = nbad + 1
    ! the data: five per line, nt/5 lines plus the remainder
    nlines = 0
    do
      read(iunit,'(a)',iostat=ios) fname
      if (ios /= 0) exit
      nlines = nlines + 1
    enddo
    if (nlines /= (nt + 4)/5) nbad = nbad + 1
    close(iunit)
  else
    nbad = nbad + 1
  endif
  call gf_report_true('alphanumeric file: cards and data lines as written',nbad == 0,nfail)

  !--- errors ----------------------------------------------------------------

  hdr%npts = nt + 1
  call gf_write_sac_trace(hdr,seis(1,1,:),nt,OUTDIR//'/IU.SJG.BXN.bad',.true.,.false.,ierr)
  call gf_report_true('a header/trace length mismatch is refused',ierr /= GF_OK,nfail)

  deallocate(seis,dp,data,db%stations)

  end subroutine test_files

!
!-------------------------------------------------------------------------------------------------
!

  subroutine note(cond,what,nbad,ival)

! one member of a boolean group: names itself in the log when it fails, so
! that the group's verdict line is not the only trace of which field was
! wrong (the CI log is all there is of a run there); ival, when given, is
! the value found

  implicit none
  logical, intent(in) :: cond
  character(len=*), intent(in) :: what
  integer, intent(inout) :: nbad
  integer, intent(in), optional :: ival

  if (cond) return
  if (present(ival)) then
    write(*,'(a,i0)') '       mismatch: '//what//' = ',ival
  else
    write(*,'(a)') '       mismatch: '//what
  endif
  nbad = nbad + 1

  end subroutine note

!
!-------------------------------------------------------------------------------------------------
!

  subroutine note_real(what,f,x,err)

! one single-precision header field against the double the writer was
! given: |f - x| / max(1,|x|), accumulated into the group's error and named
! when above TOL_SINGLE

  implicit none
  character(len=*), intent(in) :: what
  real, intent(in) :: f
  double precision, intent(in) :: x
  double precision, intent(inout) :: err

  double precision :: e

  e = abs(dble(f) - x)/max(1.d0,abs(x))
  if (e > TOL_SINGLE) write(*,'(a,es12.5,a,es12.5)') '       mismatch: '//what//' = ',f,'  expected ',x
  err = max(err,e)

  end subroutine note_real

  end program test_gf_sac
