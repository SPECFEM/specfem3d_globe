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
!---- and compares Richardson-extrapolated central differences with gf_seis's
!---- analytic partials. Nothing FD-shaped exists in the library.
!----
!---- What is asserted:
!----
!----   1. the seismograms gf_seis writes beside the partials (kind 2) are
!----      its seismograms alone (kind 0) -- bitwise under a value-safe FP
!----      model, because both contract the same block with the same weight
!----      vector and convert it the same way, and the kind 0 seismogram is
!----      what the comparison gate in EXAMPLES/ rests on;
!----   2. linearity on real data: SUM_v M_v dp(v) with the CMTSOLUTION's own
!----      dyne-cm reproduces the seismogram (asserted in section 5, against
!----      the size of the terms summed);
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
!----   5. the weights (gf_seis_weights): every station's element block,
!----      contracted with the seismogram's weight vector and with each of
!----      the nine partial columns, scaled and converted as gf_seis does,
!----      reproduces gf_seis's seismogram and partials to rounding. The
!----      derivative weights sum to zero over the element, so a trace is a
!----      small difference of large terms, and what a different summation
!----      order moves is measured against the size of those terms, the
!----      converted SUM|u w|; the error against each trace's own peak is
!----      printed. 1e-13 of SUM|u w| is a margin over the weights route's
!----      own worst case, 375 eps; a wrong index or term is ~1e-6 or more.
!----      kind 0 and 1 must give kind 2's weights, and wrong requests are
!----      refused.
!----
!---- Needs a built example (GFDB with a validation CMTSOLUTION); one
!---- station is enough for the sweep, every station for the rest.
!----

  program test_gf_partials_db

  use gf_par, only: t_gfdb,t_gf_source,t_gf_location,t_gf_taxis,t_gf_stf, &
                    gf_errmsg,gf_error_string,GF_OK,GF_ERR_ARG,GF_NCOMP,GF_SRC_CMT
  use gf_database, only: gf_open,gf_close,gf_topo_elevation
  use gf_locate, only: gf_locate_source,gf_locate_release
  use gf_source, only: gf_read_source
  use gf_seismograms, only: gf_seis_plan,gf_seis,gf_seis_weights,gf_seis_station_scale
  use gf_element_io, only: gf_read_element_displ
  use gf_weights, only: GF_NW_MT,GF_NW_LOC
  use gf_stf, only: t_gf_stf_work,gf_stf_work_init,gf_stf_convert,gf_stf_work_free

  use gf_stf, only: gf_default_t0
  use gf_partials, only: gf_partials_ndp,GF_NDP_LOC,GF_DP_NAME,GF_DP_LAT,GF_DP_LON,GF_DP_DEP,GF_DP_TIM

  use gf_manufactured, only: gf_report,gf_report_true,gf_quiet_nan

  use constants, only: MAX_STRING_LEN,CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  ! the step sweep: degrees for lat/lon, km for depth. 4e-3 degrees is 6 %
  ! of a 4-arc-minute topography cell; the stencil stays inside one cell
  ! for most sources, and the elevation check below says when it does not.
  integer, parameter :: NH = 4
  double precision, dimension(NH), parameter :: H_DEG = (/ 4.d-3, 2.d-3, 1.d-3, 5.d-4 /)
  double precision, dimension(NH), parameter :: H_KM  = (/ 1.d-1, 5.d-2, 2.5d-2, 1.25d-2 /)

  character(len=MAX_STRING_LEN) :: dbpath,cmtfile
  type(t_gfdb) :: db
  type(t_gf_source) :: src,srcp,src_bad
  type(t_gf_location) :: loc,locp
  type(t_gf_taxis) :: tax,tax_bad
  type(t_gf_stf) :: stf,stf_bad
  double precision, dimension(:,:,:), allocatable :: seis,seis0,fp,fm
  double precision, dimension(:,:,:,:), allocatable :: dp
  double precision, dimension(:,:,:), allocatable :: dfd            ! (NH, ncomp, nt) per parameter, one station
  double precision, dimension(:), allocatable :: t,onset,onsetp
  ! gf_seis takes the partials array by explicit shape; the seismogram-only
  ! calls below pass a zero-sized one
  double precision, dimension(:,:,:,:), allocatable :: dp_none
  double precision, dimension(6) :: m_dynecm
  double precision :: worst,ref,err_h(NH),ratio,rich,elev_p,elev_m,elev_0,dt_sub,h,s0,worst_bit,t0
  integer :: nfail,ierr,nargs,ndp,ista,icomp,it,nt,ip,ih,v,nbad,ista_sweep,ncount
  logical :: same_elem,in_db,in_cell
  character(len=8) :: pname
  ! section 5
  integer, parameter :: NW = GF_NW_MT + GF_NW_LOC
  double precision, dimension(GF_NCOMP,NGLLX,NGLLY,NGLLZ) :: w
  double precision, dimension(GF_NCOMP,NGLLX,NGLLY,NGLLZ,NW) :: wcols
  double precision, dimension(GF_NCOMP,NGLLX,NGLLY,NGLLZ) :: wv
  double precision, dimension(NW) :: colscale
  double precision, dimension(0:NW) :: worst_col,worst_cond,cond_col
  double precision, dimension(0:NW,GF_NCOMP) :: bound_col
  double precision :: worst_lin
  double precision, dimension(:), allocatable :: y_w
  double precision, dimension(:,:,:,:,:), allocatable :: wk1
  double precision :: bound
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: displ
  type(t_gf_stf_work) :: work
  double precision :: sta,sc
  integer :: icol

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

  call gf_read_source(db,cmtfile,src,ierr)
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

  call gf_default_t0(src,t0,ierr)
  call gf_report_true('gf_default_t0                     ',ierr == GF_OK,nfail)

  call gf_seis_plan(db,src,t0,tax,stf,ierr)
  call gf_report_true('plan: Heaviside conversion        ',ierr == GF_OK,nfail)
  if (ierr /= GF_OK) then
    write(*,'(a)') '  '//trim(gf_errmsg)
    call gf_locate_release()
    call gf_close(db)
    stop 1
  endif

  ! t0 and hdur set the padding count and the kernel half length, so a NaN
  ! there becomes an array bound rather than a NaN in the output. The screen
  ! is in gf_seis_plan, which is what both routes call.
  call gf_seis_plan(db,src,gf_quiet_nan(),tax_bad,stf_bad,ierr)
  call gf_report_true('plan: a NaN t0 is refused         ',ierr == GF_ERR_ARG,nfail)

  src_bad = src
  src_bad%hdur = gf_quiet_nan()
  call gf_seis_plan(db,src_bad,t0,tax_bad,stf_bad,ierr)
  call gf_report_true('plan: a NaN hdur is refused       ',ierr == GF_ERR_ARG,nfail)

  call gf_seis_plan(db,src,-1.d0,tax_bad,stf_bad,ierr)
  call gf_report_true('plan: an unresolved t0 is refused ',ierr == GF_ERR_ARG,nfail)
  nt = tax%nt
  dt_sub = tax%dt_sub

  call gf_partials_ndp(2,ndp,ierr)
  allocate(seis(db%nstations,GF_NCOMP,nt),seis0(db%nstations,GF_NCOMP,nt), &
           fp(db%nstations,GF_NCOMP,nt),fm(db%nstations,GF_NCOMP,nt), &
           dp(ndp,db%nstations,GF_NCOMP,nt),t(nt),onset(db%nstations),onsetp(db%nstations), &
           dfd(NH,GF_NCOMP,nt))

  !--------------------------------------------------------------------
  ! 1. the seismograms beside the partials are the kind 0 ones
  !--------------------------------------------------------------------

  write(*,'(a)') '1. seismograms'
  allocate(dp_none(0,db%nstations,GF_NCOMP,nt))
  call gf_seis(db,src,loc,tax,stf,0,0,seis0,dp_none,t,onset,ierr)
  call gf_report_true('gf_seis, kind 0             ',ierr == GF_OK,nfail)
  call gf_seis(db,src,loc,tax,stf,2,ndp,seis,dp,t,onset,ierr)
  call gf_report_true('gf_seis, kind 2             ',ierr == GF_OK,nfail)
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
  call gf_report('  seis == kind 0 (rel)        ',worst_bit,1.d-15,nfail)
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
  ! The seismogram and each partial are contractions with their own weight
  ! vectors, which are linear in the tensor only to rounding; the traces
  ! are much smaller than what is summed into them, so the identity is
  ! asserted in section 5 against that size. Against the peak: printed.
  write(*,'(a,es10.3,a)') '     SUM_v M_v dp(v) vs seismogram: ',worst, &
                          ' of the trace peak (asserted in section 5)'

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
      call gf_seis(db,srcp,locp,tax,stf,0,0,fp,dp_none,t,onsetp,ierr)
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
      call gf_seis(db,srcp,locp,tax,stf,0,0,fm,dp_none,t,onsetp,ierr)
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

      ! The convergence ratio is only a statement about truncation while
      ! truncation is what is being measured. Once the finite difference has
      ! bottomed out on round-off the errors stop falling -- they wander, and
      ! may even rise as h shrinks -- and a ratio taken there is a ratio of
      ! two noise samples. Asserting it would be asserting noise, so the
      ! ratio is claimed only when the coarsest step is well clear of the
      ! floor, estimated as the smallest error in the sweep.
      !
      ! Both cases are real. On the shipped example all three parameters are
      ! truncation-dominated. On a field smooth enough that the third
      ! derivative along one direction is small -- the synthetic fixture's
      ! depth direction, for instance -- that parameter is noise-limited from
      ! the coarsest step onwards.
      !
      ! The Richardson value is asserted either way: it compares the
      ! extrapolated derivative against the analytic one, which is the
      ! property this section exists for, and a noise-limited sweep simply
      ! makes it a tighter check rather than a meaningless one.
      if (err_h(1) >= 20.d0*minval(err_h)) then
        call gf_report_true('  '//trim(pname)//': second-order convergence (ratio in [3, 5])', &
                            ratio > 3.d0 .and. ratio < 5.d0,nfail)
      else
        write(*,'(a,a,a,es10.3,a)') '     ',trim(pname), &
          ': finite difference is round-off limited (floor ',minval(err_h), &
          '); convergence ratio reported, not asserted'
      endif

      call gf_report('  '//trim(pname)//': Richardson FD vs analytic (rel)',rich,1.d-6,nfail)
      ncount = ncount + 1
    endif
  enddo

  ! How many of the three qualify is a property of this database's element
  ! size and topography grid, not of the library: on a coarser mesh a 0.01
  ! degree step stays inside one element for all three, on a finer one for
  ! none. Each parameter that did qualify was asserted above; the count is
  ! reported so that a run comparing nothing is visible in the log.
  write(*,'(a,i0,a)') '     ',ncount,' of 3 position parameters stayed in-element and in-cell'

  !--------------------------------------------------------------------
  ! 5. the weights: the block contracted, then converted
  !--------------------------------------------------------------------

  write(*,'(a)') '5. weights'

  ! the source's own position again, and its seismogram and partials
  call gf_locate_source(db,src%latitude,src%longitude,src%depth,loc,ierr)
  if (ierr /= GF_OK) call die5('relocating the source')
  call gf_seis(db,src,loc,tax,stf,2,ndp,seis,dp,t,onset,ierr)
  if (ierr /= GF_OK) call die5('gf_seis, kind 2')

  call gf_seis_weights(db,src,loc,2,NW,w,wcols,colscale,ierr)
  call gf_report_true('gf_seis_weights, kind 2           ',ierr == GF_OK,nfail)
  if (ierr /= GF_OK) call die5('gf_seis_weights')

  allocate(displ(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ,db%nt_subsampled),y_w(nt))
  call gf_stf_work_init(stf,tax,work,ierr)
  if (ierr /= GF_OK) call die5('gf_stf_work_init')

  worst_col(:) = 0.d0
  worst_cond(:) = 0.d0
  cond_col(:) = 0.d0
  worst_lin = 0.d0
  do ista = 1,db%nstations
    call gf_read_element_displ(db,loc%ielem,ista,displ,ierr)
    if (ierr /= GF_OK) call die5('gf_read_element_displ')
    sta = gf_seis_station_scale(db,src,ista)
    do icol = 0,NW
      if (icol == 0) then
        wv(:,:,:,:) = w(:,:,:,:)
        sc = sta
      else
        wv(:,:,:,:) = wcols(:,:,:,:,icol)
        sc = sta*colscale(icol)
      endif
      do icomp = 1,GF_NCOMP
        ! the trace, through the weights
        do it = 1,db%nt_subsampled
          work%trace(it) = sc * sum(dble(displ(icomp,:,:,:,:,it))*wv(:,:,:,:))
        enddo
        call gf_stf_convert(stf,tax,work,ierr)
        if (ierr /= GF_OK) call die5('gf_stf_convert')
        y_w(:) = work%y(1:nt)
        ! the same with every term made positive: the size of what is summed,
        ! which bounds what a different summation order can change. The
        ! conversion's kernel is positive, so this bounds it after conversion
        ! too.
        do it = 1,db%nt_subsampled
          work%trace(it) = abs(sc) * sum(abs(dble(displ(icomp,:,:,:,:,it)))*abs(wv(:,:,:,:)))
        enddo
        call gf_stf_convert(stf,tax,work,ierr)
        if (ierr /= GF_OK) call die5('gf_stf_convert')
        bound = maxval(abs(work%y(1:nt)))
        bound_col(icol,icomp) = bound
        if (icol == 0) then
          ref = maxval(abs(seis(ista,icomp,:)))
          err_h(1) = maxval(abs(y_w(:) - seis(ista,icomp,:)))
        else
          ref = maxval(abs(dp(icol,ista,icomp,:)))
          err_h(1) = maxval(abs(y_w(:) - dp(icol,ista,icomp,:)))
        endif
        ! gated on the terms, not on the library's trace: a column that came
        ! back all zeros must fail, not be skipped
        if (bound > 0.d0) worst_cond(icol) = max(worst_cond(icol),err_h(1)/bound)
        if (ref > 0.d0) then
          worst_col(icol) = max(worst_col(icol),err_h(1)/ref)
          cond_col(icol) = max(cond_col(icol),bound/ref)
        endif
      enddo
    enddo
    ! linearity, against the size of every term it sums
    do icomp = 1,GF_NCOMP
      bound = bound_col(0,icomp)
      do v = 1,6
        bound = bound + abs(m_dynecm(v))*bound_col(v,icomp)
      enddo
      do it = 1,nt
        s0 = 0.d0
        do v = 1,6
          s0 = s0 + m_dynecm(v)*dp(v,ista,icomp,it)
        enddo
        if (bound > 0.d0) worst_lin = max(worst_lin,abs(s0 - seis(ista,icomp,it))/bound)
      enddo
    enddo
  enddo
  ! The derivative weights sum to zero over the element, so a trace is a
  ! difference of terms that can be far larger than it: the fixture's field
  ! has a large smooth part. A different summation order then moves the
  ! result by rounding of the terms, eps * SUM |u w|, not of the result. The
  ! assertion is against that size; the error against the trace's own peak
  ! and the ratio of the two sizes are printed beside it.
  do icol = 0,NW
    write(*,'(a,a4,a,es10.3,a,es10.3)') '     ',merge('seis',GF_DP_NAME(max(icol,1))//' ',icol == 0), &
      ': error / trace peak ',worst_col(icol),'   SUM|u w| / trace peak ',cond_col(icol)
  enddo
  call gf_report('weights -> seismogram, / SUM|u w|      ',worst_cond(0),1.d-13,nfail)
  call gf_report('SUM_v M_v dp(v) == seismogram, / SUM|u w|',worst_lin,1.d-13,nfail)
  do icol = 1,NW
    call gf_report('weights -> partial '//GF_DP_NAME(icol)//', / SUM|u w|',worst_cond(icol),1.d-13,nfail)
  enddo

  ! the other kinds give the same weights, and wrong requests are refused
  call gf_seis_weights(db,src,loc,0,0,wv,wcols(:,:,:,:,1:0),colscale(1:0),ierr)
  call gf_report_true('kind 0: the seismogram weights, bitwise ', &
                      ierr == GF_OK .and. all(wv == w),nfail)
  allocate(wk1(GF_NCOMP,NGLLX,NGLLY,NGLLZ,GF_NW_MT))
  call gf_seis_weights(db,src,loc,1,GF_NW_MT,wv,wk1,colscale(1:GF_NW_MT),ierr)
  call gf_report_true('kind 1: and the six MT columns, bitwise  ', &
                      ierr == GF_OK .and. all(wv == w) .and. all(wk1 == wcols(:,:,:,:,1:GF_NW_MT)),nfail)
  call gf_seis_weights(db,src,loc,1,NW,wv,wcols,colscale,ierr)
  call gf_report_true('kind 1 with 9 columns is refused        ',ierr == GF_ERR_ARG,nfail)
  call gf_seis_weights(db,src,loc,3,0,wv,wcols(:,:,:,:,1:0),colscale(1:0),ierr)
  call gf_report_true('kind 3 is refused                       ',ierr == GF_ERR_ARG,nfail)
  deallocate(wk1)

  call gf_stf_work_free(work)
  deallocate(displ,y_w)

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

  subroutine die5(what)

  implicit none
  character(len=*), intent(in) :: what

  write(*,'(a,a,a,a)') '  section 5: ',what,' failed: ',trim(gf_errmsg)
  stop 1

  end subroutine die5

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
