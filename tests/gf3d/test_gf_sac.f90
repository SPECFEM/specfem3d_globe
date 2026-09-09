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
!----   * the data round-trips to single precision bitwise;
!----   * a partial's file name and KUSER1/KUSER2 name the parameter.
!----
!---- The oracle is the file format and the solver's expressions, not
!---- another writer.
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
  double precision :: t
  integer :: ierr,nt,it,icomp,iunit,ios,nbad,nlines
  logical :: exists
  character(len=256) :: fname
  character(len=3), dimension(GF_NCOMP) :: chn

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

  allocate(data(nt))
  nbad = 0
  do icomp = 1,GF_NCOMP
    fname = OUTDIR//'/IU.SJG.'//chn(icomp)//'.sem.sac'
    inquire(file=trim(fname),exist=exists)
    if (.not. exists) then
      nbad = nbad + 1
      cycle
    endif
    open(newunit=iunit,file=trim(fname),access='stream',form='unformatted',status='old',action='read',iostat=ios)
    read(iunit,iostat=ios) fh,ih,ch,data
    close(iunit)
    if (ios /= 0) then
      nbad = nbad + 1
      cycle
    endif

    ! reals: DELTA(1) B(6) O(8) STLA(32) STLO(33) STEL(34) STDP(35) EVLA(36)
    ! EVLO(37) EVDP(39) USER0(41) CMPAZ(58) CMPINC(59)
    if (fh(1) /= real(tax%dt_sub) .or. fh(6) /= real(tax%t_first) .or. fh(8) /= 0.0) nbad = nbad + 1
    if (fh(32) /= real(18.1091d0) .or. fh(33) /= real(-66.15d0) .or. fh(34) /= -12345.0 .or. fh(35) /= 0.0) nbad = nbad + 1
    if (fh(36) /= real(-5.812d0) .or. fh(37) /= real(-75.27d0) .or. fh(39) /= real(122.6d0)) nbad = nbad + 1
    if (fh(41) /= 60.0 .or. fh(43) /= 500.0) nbad = nbad + 1
    select case (icomp)
    case (1)
      if (fh(58) /= 0.0 .or. fh(59) /= 90.0) nbad = nbad + 1
    case (2)
      if (fh(58) /= 90.0 .or. fh(59) /= 90.0) nbad = nbad + 1
    case (3)
      if (fh(58) /= 0.0 .or. fh(59) /= 0.0) nbad = nbad + 1
    end select
    ! integers: NZYEAR..NZMSEC (1..6), NVHDR(7), NPTS(10), IFTYPE(16),
    ! IDEP(17), IZTYPE(18), LEVEN(36), LCALDA(39)
    if (ih(1) /= 1994 .or. ih(2) /= 160 .or. ih(3) /= 0 .or. ih(4) /= 33 .or. ih(5) /= 45 .or. ih(6) /= 399) nbad = nbad + 1
    if (ih(7) /= 6 .or. ih(10) /= nt .or. ih(16) /= 1 .or. ih(17) /= 6 .or. ih(18) /= 11) nbad = nbad + 1
    if (ih(36) /= 1 .or. ih(39) /= 1) nbad = nbad + 1
    ! strings: KSTNM(1:8) KEVNM(9:24) KHOLE(25:32) KUSER0(137:144) KUSER1(145:152)
    ! KUSER2(153:160) KCMPNM(161:168) KNETWK(169:176)
    if (ch(1:8) /= 'SJG     ' .or. ch(9:24) /= '060994A         ' .or. ch(25:32) /= 'S3      ') nbad = nbad + 1
    if (ch(137:144) /= 'SY      ' .or. ch(145:152) /= 'gf3d    ' .or. ch(161:168) /= chn(icomp)//'     ' .or. &
        ch(169:176) /= 'IU      ') nbad = nbad + 1
    ! the data, single precision, bitwise
    do it = 1,nt
      if (data(it) /= real(seis(1,icomp,it))) nbad = nbad + 1
    enddo
  enddo
  call gf_report_true('binary seismograms: header fields and data as written',nbad == 0,nfail)

  !--- the binary partial ----------------------------------------------------

  nbad = 0
  fname = OUTDIR//'/IU.SJG.BXZ.Mrr.sem.sac'
  inquire(file=trim(fname),exist=exists)
  if (exists) then
    open(newunit=iunit,file=trim(fname),access='stream',form='unformatted',status='old',action='read',iostat=ios)
    read(iunit,iostat=ios) fh,ih,ch,data
    close(iunit)
    if (ios /= 0) nbad = nbad + 1
    if (ch(145:152) /= 'Mrr     ' .or. ch(153:160) /= 'm/dynecm') nbad = nbad + 1
    if (ch(161:168) /= 'BXZ     ' .or. ih(10) /= nt .or. fh(6) /= real(tax%t_first)) nbad = nbad + 1
    do it = 1,nt
      if (data(it) /= real(dp(1,1,3,it))) nbad = nbad + 1
    enddo
  else
    nbad = nbad + 1
  endif
  call gf_report_true('binary partial: named in the file and in KUSER1/KUSER2, data as written',nbad == 0,nfail)

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

  end program test_gf_sac
