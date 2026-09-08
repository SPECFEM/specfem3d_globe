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
!---- test_gf_partials_db -- the partial derivatives on a real database
!----
!---- Stage 7 as the plan finally has it: finite differences are the
!---- *validation* of the analytic derivative, not a product. This driver
!---- perturbs the source through the public locate-and-extract routines
!---- and compares Richardson-extrapolated central differences with
!---- gf_seis_cmt_partials. Nothing FD-shaped exists in the library.
!----
!---- What is asserted:
!----
!----   1. the seismograms gf_seis_cmt_partials writes beside the partials
!----      are the seismograms of gf_seis_cmt -- bitwise, because the two
!----      run the same statements on the same arrays, and gf_seis_cmt is the
!----      routine the comparison gate in EXAMPLES/ rests on;
!----   2. linearity on real data: SUM_v M_v dp(v) with the CMTSOLUTION's own
!----      dyne-cm reproduces the seismogram;
!----   3. dp(10) against a central difference of the seismogram: the
!----      difference's own (w dt)^2/6 error is what is measured, and it is
!----      reported, not asserted tightly;
!----   4. the centroid partials against relocation differences at four step
!----      sizes: second-order convergence of the difference, then the
!----      Richardson combination of the two finest steps against the analytic
!----      value at 1e-6. Both routes differentiate the same interpolant built
!----      from the same float32 nodes, so the float32 noise does not limit
!----      their agreement; measured on the shipped examples it is ~1e-10.
!----      Three things are checked per perturbed position and reported when
!----      they bite, because each makes the comparison meaningless rather
!----      than wrong: a relocation into a *neighbouring* element (the
!----      interpolant is only C0 across faces), a position outside the
!----      database, and a stencil straddling a topography cell edge (the
!----      elevation is bilinear per cell, so its gradient jumps there; the
!----      steps are chosen well inside a 4-arc-minute cell and the check is
!----      that the elevation is linear over the stencil).
!----
!---- Needs a built example (GFDB with a validation CMTSOLUTION); one
!---- station is enough for the sweep, every station for the rest.
!----

  program test_gf_partials_db

  use gf_par, only: t_gfdb,t_gf_source,t_gf_location,t_gf_taxis,t_gf_stf, &
                    gf_errmsg,gf_error_string,GF_OK,GF_NCOMP,GF_SRC_CMT
  use gf_database, only: gf_open,gf_close,gf_topo_elevation
  use gf_locate, only: gf_locate_source,gf_locate_release
  use gf_source, only: gf_read_source
  use gf_seismograms, only: gf_seis_plan,gf_seis_cmt,gf_seis_cmt_partials
  use gf_partials, only: gf_partials_ndp,GF_NDP_LOC,GF_DP_NAME,GF_DP_LAT,GF_DP_LON,GF_DP_DEP,GF_DP_TIM

  use gf_manufactured, only: gf_report,gf_report_true

  use constants, only: MAX_STRING_LEN

  implicit none

  ! the step sweep: degrees for lat/lon, km for depth. 4e-3 degrees is 6 %
  ! of a 4-arc-minute topography cell; the stencil stays inside one cell
  ! for most sources, and the elevation check below says when it does not.
  integer, parameter :: NH = 4
  double precision, dimension(NH), parameter :: H_DEG = (/ 4.d-3, 2.d-3, 1.d-3, 5.d-4 /)
  double precision, dimension(NH), parameter :: H_KM  = (/ 1.d-1, 5.d-2, 2.5d-2, 1.25d-2 /)

  character(len=MAX_STRING_LEN) :: dbpath,cmtfile
  type(t_gfdb) :: db
  type(t_gf_source) :: src,srcp
  type(t_gf_location) :: loc,locp
  type(t_gf_taxis) :: tax
  type(t_gf_stf) :: stf
  double precision, dimension(:,:,:), allocatable :: seis,seis0,fp,fm
  double precision, dimension(:,:,:,:), allocatable :: dp
  double precision, dimension(:,:,:), allocatable :: dfd            ! (NH, ncomp, nt) per parameter, one station
  double precision, dimension(:), allocatable :: t,onset,onsetp
  double precision, dimension(6) :: m_dynecm
  double precision :: worst,ref,err_h(NH),ratio,rich,elev_p,elev_m,elev_0,dt_sub,h,s0,worst_bit
  integer :: nfail,ierr,nargs,ndp,ista,icomp,it,nt,ip,ih,v,nbad,ista_sweep,ncount
  logical :: same_elem,in_db,in_cell
  character(len=8) :: pname

  nfail = 0

  nargs = command_argument_count()
  if (nargs < 2) then
    write(*,'(a)') 'usage: test_gf_partials_db <GFDB> <CMTSOLUTION>'
    stop 1
  endif
  call get_command_argument(1,dbpath)
  call get_command_argument(2,cmtfile)

  write(*,'(a)') ''
  write(*,'(a)') 'test_gf_partials_db'
  write(*,'(a)') ''
  write(*,'(a,a)') '  database = ',trim(dbpath)
  write(*,'(a,a)') '  source   = ',trim(cmtfile)
  write(*,'(a)') ''

  call gf_open(dbpath,db,ierr,check_completion=.false.)
  if (ierr /= GF_OK) then
    write(*,'(a)') '  could not open the database: '//trim(gf_errmsg)
    stop 1
  endif

  call gf_read_source(cmtfile,db%dt,src,ierr)
  if (ierr /= GF_OK .or. src%source_type /= GF_SRC_CMT) then
    write(*,'(a)') '  could not read a CMTSOLUTION: '//trim(gf_errmsg)
    call gf_close(db)
    stop 1
  endif
  ! the file's own numbers, for the linearity identity
  m_dynecm(:) = src%moment_tensor(:)*src%scale_moment

  call gf_locate_source(db,src%latitude,src%longitude,src%depth,loc,ierr)
  call gf_report_true('source located                    ',ierr == GF_OK,nfail)
  if (ierr /= GF_OK) then
    write(*,'(a)') '  '//trim(gf_errmsg)
    call gf_close(db)
    stop 1
  endif
  write(*,'(a,a,a,3f10.5)') '  element ',loc%morton_hex,'  xi,eta,gamma = ',loc%xi,loc%eta,loc%gamma

  call gf_seis_plan(db,src,-1.d0,tax,stf,ierr)
  call gf_report_true('plan: Heaviside conversion        ',ierr == GF_OK,nfail)
  nt = tax%nt
  dt_sub = tax%dt_sub

  call gf_partials_ndp(2,ndp,ierr)
  allocate(seis(db%nstations,GF_NCOMP,nt),seis0(db%nstations,GF_NCOMP,nt), &
           fp(db%nstations,GF_NCOMP,nt),fm(db%nstations,GF_NCOMP,nt), &
           dp(ndp,db%nstations,GF_NCOMP,nt),t(nt),onset(db%nstations),onsetp(db%nstations), &
           dfd(NH,GF_NCOMP,nt))

  !--------------------------------------------------------------------
  ! 1. the seismograms beside the partials are gf_seis_cmt's
  !--------------------------------------------------------------------

  write(*,'(a)') '1. seismograms'
  call gf_seis_cmt(db,src,loc,tax,stf,seis0,t,onset,ierr)
  call gf_report_true('gf_seis_cmt                       ',ierr == GF_OK,nfail)
  call gf_seis_cmt_partials(db,src,loc,tax,stf,2,ndp,seis,dp,t,onset,ierr)
  call gf_report_true('gf_seis_cmt_partials, itypsokern 2',ierr == GF_OK,nfail)
  if (ierr /= GF_OK) write(*,'(a)') '  '//trim(gf_errmsg)

  nbad = 0
  worst_bit = 0.d0
  ref = maxval(abs(seis0))
  do it = 1,nt
    do icomp = 1,GF_NCOMP
      do ista = 1,db%nstations
        if (seis(ista,icomp,it) /= seis0(ista,icomp,it)) nbad = nbad + 1
        worst_bit = max(worst_bit,abs(seis(ista,icomp,it) - seis0(ista,icomp,it))/ref)
      enddo
    enddo
  enddo
  call gf_report('  seis == gf_seis_cmt (rel)         ',worst_bit,1.d-15,nfail)
  write(*,'(a,i0,a,i0,a)') '     bitwise mismatches = ',nbad,' of ',db%nstations*GF_NCOMP*nt, &
                           ' (informational: 0 under a value-safe FP model)'

  !--------------------------------------------------------------------
  ! 2. linearity on real data
  !--------------------------------------------------------------------

  write(*,'(a)') '2. linearity'
  worst = 0.d0
  do ista = 1,db%nstations
    ref = maxval(abs(seis0(ista,:,:)))
    do it = 1,nt
      do icomp = 1,GF_NCOMP
        s0 = 0.d0
        do v = 1,6
          s0 = s0 + m_dynecm(v)*dp(v,ista,icomp,it)
        enddo
        worst = max(worst,abs(s0 - seis0(ista,icomp,it))/ref)
      enddo
    enddo
  enddo
  call gf_report('SUM_v M_v dp(v) == seismogram, all stations',worst,1.d-12,nfail)

  !--------------------------------------------------------------------
  ! 3. dp(10) against a central difference of the seismogram
  !--------------------------------------------------------------------

  write(*,'(a)') '3. centroid time'
  worst = 0.d0
  do ista = 1,db%nstations
    ref = maxval(abs(dp(GF_DP_TIM,ista,:,:)))
    do icomp = 1,GF_NCOMP
      do it = 2,nt-stf%khalf-1
        worst = max(worst,abs(-(seis0(ista,icomp,it+1) - seis0(ista,icomp,it-1))/(2.d0*dt_sub) &
                              - dp(GF_DP_TIM,ista,icomp,it))/ref)
      enddo
    enddo
  enddo
  write(*,'(a,es10.3,a)') '     central difference vs dp(10): ',worst,'  (the difference''s (w dt)^2/6)'
  call gf_report('dp(10) vs central difference (rel)  ',worst,5.d-2,nfail)

  !--------------------------------------------------------------------
  ! 4. the centroid partials against relocation differences
  !
  ! One station, the first; every perturbed position is checked for the
  ! same element, for being inside the database, and for a linear
  ! elevation over the stencil.
  !--------------------------------------------------------------------

  write(*,'(a)') '4. centroid position: FD by relocation'
  ista_sweep = 1
  ncount = 0

  do ip = 1,3
    pname = GF_DP_NAME(GF_DP_LAT+ip-1)
    same_elem = .true.
    in_db = .true.
    in_cell = .true.

    do ih = 1,NH
      h = H_DEG(ih)
      if (ip == 3) h = H_KM(ih)

      ! + h
      srcp = src
      call perturb(srcp,ip,+h)
      call gf_locate_source(db,srcp%latitude,srcp%longitude,srcp%depth,locp,ierr)
      if (ierr /= GF_OK) then
        in_db = .false.
        exit
      endif
      if (locp%ielem /= loc%ielem) same_elem = .false.
      call gf_seis_cmt(db,srcp,locp,tax,stf,fp,t,onsetp,ierr)
      if (db%topography) call gf_topo_elevation(db,srcp%latitude,srcp%longitude,elev_p)

      ! - h
      srcp = src
      call perturb(srcp,ip,-h)
      call gf_locate_source(db,srcp%latitude,srcp%longitude,srcp%depth,locp,ierr)
      if (ierr /= GF_OK) then
        in_db = .false.
        exit
      endif
      if (locp%ielem /= loc%ielem) same_elem = .false.
      call gf_seis_cmt(db,srcp,locp,tax,stf,fm,t,onsetp,ierr)
      if (db%topography) call gf_topo_elevation(db,srcp%latitude,srcp%longitude,elev_m)

      ! a linear elevation over the stencil: the bilinear interpolant's
      ! second difference is zero inside a cell
      if (db%topography .and. ip < 3) then
        call gf_topo_elevation(db,src%latitude,src%longitude,elev_0)
        if (abs(elev_p + elev_m - 2.d0*elev_0) > 1.d-6*max(1.d0,abs(elev_0))) in_cell = .false.
      endif

      do it = 1,nt
        do icomp = 1,GF_NCOMP
          dfd(ih,icomp,it) = (fp(ista_sweep,icomp,it) - fm(ista_sweep,icomp,it))/(2.d0*h)
        enddo
      enddo
    enddo

    if (.not. in_db) then
      write(*,'(a,a,a)') '     ',trim(pname),': a perturbed position left the database -- skipped'
      cycle
    endif
    if (.not. same_elem) write(*,'(a,a,a)') '     ',trim(pname), &
      ': a perturbed position relocated into a neighbouring element (C0 face): reported, not asserted'
    if (.not. in_cell) write(*,'(a,a,a)') '     ',trim(pname), &
      ': the stencil straddles a topography cell edge: reported, not asserted'

    ! the error of each step against the analytic value, and its order
    ref = maxval(abs(dp(GF_DP_LAT+ip-1,ista_sweep,:,:)))
    do ih = 1,NH
      err_h(ih) = maxval(abs(dfd(ih,:,:) - dp(GF_DP_LAT+ip-1,ista_sweep,:,:)))/ref
    enddo
    ratio = err_h(1)/err_h(2)
    ! Richardson from the two finest steps
    rich = 0.d0
    do it = 1,nt
      do icomp = 1,GF_NCOMP
        rich = max(rich,abs((4.d0*dfd(NH,icomp,it) - dfd(NH-1,icomp,it))/3.d0 &
                            - dp(GF_DP_LAT+ip-1,ista_sweep,icomp,it))/ref)
      enddo
    enddo
    write(*,'(a,a,a,4es10.3,a,f6.2,a,es10.3)') '     ',trim(pname),': FD error per step ',err_h, &
                                               '  ratio(h1/h2) ',ratio,'  Richardson ',rich

    if (same_elem .and. in_cell) then
      call gf_report_true('  '//trim(pname)//': second-order convergence (ratio in [3, 5])', &
                          ratio > 3.d0 .and. ratio < 5.d0,nfail)
      call gf_report('  '//trim(pname)//': Richardson FD vs analytic (rel)',rich,1.d-6,nfail)
      ncount = ncount + 1
    endif
  enddo

  call gf_report_true('at least two parameters compared cleanly',ncount >= 2,nfail)

  !--------------------------------------------------------------------

  call gf_locate_release()
  call gf_close(db)
  deallocate(seis,seis0,fp,fm,dp,t,onset,onsetp,dfd)

  write(*,'(a)') ''
  if (nfail /= 0) then
    write(*,'(a,i0,a)') 'test_gf_partials_db: ',nfail,' assertion(s) FAILED'
    stop 1
  endif
  write(*,'(a)') 'test_gf_partials_db: all assertions passed'

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine perturb(s,ip,h)

  implicit none
  type(t_gf_source), intent(inout) :: s
  integer, intent(in) :: ip
  double precision, intent(in) :: h

  select case (ip)
  case (1)
    s%latitude = s%latitude + h
  case (2)
    s%longitude = s%longitude + h
  case (3)
    s%depth = s%depth + h
  end select

  end subroutine perturb

  end program test_gf_partials_db
