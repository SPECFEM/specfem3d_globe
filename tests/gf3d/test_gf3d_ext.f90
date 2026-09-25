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
!---- test_gf3d_ext -- the external program's view of the library
!----
!---- This is the "downstream caller" test: it is compiled with -I./include
!---- and linked against ./lib/libgf3d.a and nothing else -- no ./obj, no
!---- gf_* module files -- so it can only see what an installed library
!---- offers. If `use gf3d` does not by itself name everything a caller
!---- needs, this test does not compile, which is the point.
!----
!---- It also reaches the C entry points, and it does so the way a foreign
!---- caller would: by transcribing include/gf3d.h into its own interface
!---- blocks and bind(C) types. gf3d_capi.mod is deliberately not installed,
!---- so this is not a shortcut but the real exercise -- and it means the
!---- header and the Fortran facade are checked against a third, independent
!---- transcription of the same layout.
!----
!---- What it asserts: the two routes to the same numbers agree. The file
!---- route reads the CMTSOLUTION through the solver's own reader and calls
!---- the public Fortran API; the C route passes the same values through the
!---- ABI. Seismograms and all ten partials must match to 1e-12 relative --
!---- a derived tolerance, not equality, because the two paths are separate
!---- compilation units evaluating the same expressions (testing.md).
!----
!---- Both routes are one call to gf_extract, so section 4 asserts the part
!---- of its contract neither of them exposes: that its optional loc, tax and
!---- stf outputs describe every stage that completed, including when the
!---- call failed at a later one.
!----
!---- Usage: test_gf3d_ext <GFDB directory> <CMTSOLUTION> [FORCESOLUTION]
!----

  program test_gf3d_ext

! everything below comes from this one module: that is the claim being tested
  use gf3d

  use, intrinsic :: iso_c_binding, only: c_int,c_double,c_char,c_ptr,c_null_char, &
                                         c_null_ptr,c_loc

  implicit none

  !--- the C ABI, transcribed from include/gf3d.h ---

  integer, parameter :: NCOMP_C = 3
  integer, parameter :: MORTON_STRLEN_C = 24

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

  type, bind(C) :: gf3d_plan_t
    integer(c_int) :: nt,nt_db,npad,subsample_step,khalf,guard,kind_stf,pad_
    real(c_double) :: dt,dt_sub,t0_db,t0_req,t0,t_first
    real(c_double) :: hdur_src,hdur_target,hdur_db,hdur_corr,trunc
  end type gf3d_plan_t

  type, bind(C) :: gf3d_info_t
    integer(c_int) :: nelem,nstations,nstep,nt_subsampled,subsample_step
    integer(c_int) :: ngllx,nglly,ngllz
    integer(c_int) :: topography,ellipticity,rotation,attenuation,gravity,pad_
    real(c_double) :: dt,t0,r_planet,rhoav,scale_displ
  end type gf3d_info_t

  type, bind(C) :: gf3d_location_t
    integer(c_int) :: ielem
    character(kind=c_char), dimension(MORTON_STRLEN_C) :: morton_hex
    real(c_double) :: xi,eta,gamma
    real(c_double), dimension(3) :: xyz
    real(c_double), dimension(3) :: xyz_target
    real(c_double) :: distance_km,anchor_err,theta,phi,r_surface
  end type gf3d_location_t

  interface

    integer(c_int) function c_gf3d_open(path,check_completion,h) bind(C,name='gf3d_open')
      import :: c_int,c_char
      character(kind=c_char), dimension(*), intent(in) :: path
      integer(c_int), value :: check_completion
      integer(c_int), intent(out) :: h
    end function c_gf3d_open

    integer(c_int) function c_gf3d_close(h) bind(C,name='gf3d_close')
      import :: c_int
      integer(c_int), value :: h
    end function c_gf3d_close

    integer(c_int) function c_gf3d_get_info(h,info) bind(C,name='gf3d_get_info')
      import :: c_int,gf3d_info_t
      integer(c_int), value :: h
      type(gf3d_info_t), intent(out) :: info
    end function c_gf3d_get_info

    integer(c_int) function c_gf3d_get_plan(h,src,t0_req,plan) bind(C,name='gf3d_get_plan')
      import :: c_int,c_double,gf3d_source_t,gf3d_plan_t
      integer(c_int), value :: h
      type(gf3d_source_t), intent(in) :: src
      real(c_double), value :: t0_req
      type(gf3d_plan_t), intent(out) :: plan
    end function c_gf3d_get_plan

    integer(c_int) function c_gf3d_seismograms(h,src,t0_req,nt,seis,t,onset,loc) &
      bind(C,name='gf3d_seismograms')
      import :: c_int,c_double,c_ptr,gf3d_source_t
      integer(c_int), value :: h
      type(gf3d_source_t), intent(in) :: src
      real(c_double), value :: t0_req
      integer(c_int), value :: nt
      real(c_double), dimension(*), intent(out) :: seis,t,onset
      type(c_ptr), value :: loc
    end function c_gf3d_seismograms

    integer(c_int) function c_gf3d_partials(h,src,t0_req,kind,nt,ndp, &
                                            seis,dp,t,onset,loc) bind(C,name='gf3d_partials')
      import :: c_int,c_double,c_ptr,gf3d_source_t
      integer(c_int), value :: h
      type(gf3d_source_t), intent(in) :: src
      real(c_double), value :: t0_req
      integer(c_int), value :: kind,nt,ndp
      real(c_double), dimension(*), intent(out) :: seis,dp,t,onset
      type(c_ptr), value :: loc
    end function c_gf3d_partials

  end interface

  !--- the test ---

  double precision, parameter :: TOL = 1.d-12

  type(t_gfdb) :: db
  type(t_gf_source) :: src
  double precision, dimension(:,:,:), allocatable :: synt,synt_ref,synt_again
  double precision, dimension(:,:,:,:), allocatable :: dp

  ! section 7: walk north in these steps until the element changes
  double precision, parameter :: WALK_STEP = 0.25d0
  integer, parameter :: WALK_MAX = 40
  type(t_gf_source) :: far
  type(t_gf_location) :: loc_a,loc_b,loc_back
  double precision, dimension(:,:,:), allocatable :: synt_far
  integer :: iw
  logical :: crossed
  double precision, dimension(:), allocatable :: t,t_again
  type(t_gfdb) :: db2
  type(t_gf_location) :: loc2

  type(gf3d_source_t) :: csrc
  type(gf3d_plan_t) :: cplan
  type(gf3d_info_t) :: cinfo
  ! target: section 4 hands its address across the C boundary
  type(gf3d_location_t), target :: cloc
  real(c_double), dimension(:), allocatable :: cseis,cdp,ct,consetd
  integer(c_int) :: ch,cerr

  ! section 4: gf_extract's own outputs
  type(t_gf_location) :: locx
  type(t_gf_taxis) :: taxx
  type(t_gf_stf) :: stfx
  double precision, dimension(:), allocatable :: onset

  character(len=512) :: dbpath,cmtpath,forcepath
  double precision :: t0,worst
  integer :: ierr,nfail,narg,ista,ic,it

  nfail = 0

  write(*,*)
  write(*,*) '******************************'
  write(*,*) 'test_gf3d_ext'
  write(*,*) '******************************'
  write(*,*)

  narg = command_argument_count()
  if (narg < 2) then
    write(*,*) 'usage: test_gf3d_ext <GFDB> <CMTSOLUTION> [FORCESOLUTION]'
    stop 1
  endif
  call get_command_argument(1,dbpath)
  call get_command_argument(2,cmtpath)
  forcepath = ''
  if (narg >= 3) call get_command_argument(3,forcepath)

  !--- the Fortran route ---

  write(*,*) '1. the public Fortran module'

  call gf_open(trim(dbpath),db,ierr)
  call report_true('   gf_open',ierr == GF_OK,nfail)
  if (ierr /= GF_OK) then
    write(*,*) '   ',trim(gf_errmsg)
    stop 1
  endif

  call gf_read_cmt_source(db,trim(cmtpath),src,ierr)
  call report_true('   gf_read_cmt_source',ierr == GF_OK,nfail)
  if (ierr /= GF_OK) stop 1

  call gf_default_t0(src,t0,ierr)
  call report_true('   gf_default_t0',ierr == GF_OK,nfail)
  if (ierr /= GF_OK) stop 1

  call get_partials(db,src,t0,2,synt,dp,ierr,t)
  call report_true('   get_partials, kind 2   ',ierr == GF_OK,nfail)
  if (ierr /= GF_OK) then
    write(*,*) '   ',trim(gf_errmsg)
    stop 1
  endif

  call report_true('   synt is (nsta,3,nt)', &
    size(synt,1) == db%nstations .and. size(synt,2) == GF_NCOMP .and. &
    size(synt,3) == size(t),nfail)
  call report_true('   dp is (10,nsta,3,nt)', &
    size(dp,1) == GF_NDP_LOC .and. size(dp,2) == db%nstations,nfail)
  call report_true('   the partial names are the GF3DF order', &
    GF_DP_NAME(1) == 'Mrr' .and. GF_DP_NAME(GF_DP_TIM) == 'tim',nfail)

  write(*,'(a,i0,a,i0,a)') '       ',db%nstations,' stations, ',size(t),' samples'

  !--- the same thing through the C ABI ---

  write(*,*) '2. the same source through the C entry points'

  cerr = c_gf3d_open(trim(dbpath)//c_null_char,0_c_int,ch)
  call report_true('   gf3d_open',cerr == GF_OK,nfail)
  if (cerr /= GF_OK) stop 1

  cerr = c_gf3d_get_info(ch,cinfo)
  call report_true('   gf3d_get_info agrees with the handle', &
    cerr == GF_OK .and. int(cinfo%nstations) == db%nstations .and. &
    int(cinfo%nelem) == db%nelem,nfail)

  ! the values a caller would have parsed out of the CMTSOLUTION, recovered
  ! from what the reader produced: the moment tensor is stored
  ! non-dimensional, so scale_moment puts it back into dyne-cm, and the
  ! original time shift is in min_tshift_src_original
  csrc%source_type = GF_SRC_CMT
  csrc%force_stf = 0
  csrc%latitude = src%latitude
  csrc%longitude = src%longitude
  csrc%depth_km = src%depth
  csrc%hdur = src%hdur
  csrc%time_shift = src%min_tshift_src_original
  csrc%moment(1:6) = src%moment_tensor(1:6)*src%scale_moment
  csrc%force_factor = 0.d0
  csrc%force_dir(1:3) = 0.d0

  cerr = c_gf3d_get_plan(ch,csrc,-1.0_c_double,cplan)
  call report_true('   gf3d_get_plan',cerr == GF_OK,nfail)
  call report_true('   the plan lengths agree',int(cplan%nt) == size(t),nfail)
  if (cerr /= GF_OK) stop 1

  allocate(cseis(db%nstations*GF_NCOMP*int(cplan%nt)), &
           cdp(db%nstations*GF_NDP_LOC*GF_NCOMP*int(cplan%nt)), &
           ct(int(cplan%nt)),consetd(db%nstations))

  cerr = c_gf3d_partials(ch,csrc,-1.0_c_double,2_c_int,cplan%nt, &
                         int(GF_NDP_LOC,kind=c_int),cseis,cdp,ct,consetd,c_loc(cloc))
  call report_true('   gf3d_partials',cerr == GF_OK,nfail)
  if (cerr /= GF_OK) stop 1

  !--- and they must be the same numbers ---

  write(*,*) '3. the two routes against each other'

  call compare_seis(db%nstations,int(cplan%nt),synt,cseis,nfail)
  call compare_time(int(cplan%nt),t,ct,nfail)
  call compare_dp(db%nstations,int(cplan%nt),dp,cdp,nfail)

  allocate(synt_ref,source=synt)

  deallocate(synt,dp,t)

  !--- 4. gf_extract, the sequence both routes above run ---
  !
  ! get_partials and gf3d_partials are each one call to gf_extract, so the
  ! numbers already agree by construction. What is asserted here is the part
  ! of gf_extract's contract that neither entry point exposes: the optional
  ! loc, tax and onset outputs, and their validity after a failure. That
  ! contract is what lets xgf3d report a location for an extraction that
  ! never finished, without a printer inside the compute path.

  write(*,*) '4. gf_extract and its optional outputs'

  call gf_extract(db,src,t0,2,synt,dp,ierr,t=t,onset=onset, &
                  loc=locx,tax=taxx,stf=stfx)
  call report_true('   gf_extract, kind 2',ierr == GF_OK,nfail)

  if (ierr == GF_OK) then
    call report_true('   it reports the element it used',locx%ielem > 0,nfail)
    call report_true('   and the axis it returned',taxx%nt == size(t),nfail)
    call report_true('   the conversion plan came back too', &
                     stfx%kind_stf == GF_STF_HEAVI,nfail)
    call report_true('   the element is the C route''s', &
                     locx%ielem == int(cloc%ielem),nfail)
    call report('   xi,eta,gamma vs the C route   ', &
      max(abs(locx%xi - cloc%xi),abs(locx%eta - cloc%eta), &
          abs(locx%gamma - cloc%gamma)),TOL,nfail)
    worst = 0.d0
    do ista = 1,db%nstations
      worst = max(worst,abs(onset(ista) - consetd(ista)))
    enddo
    call report('   onset vs the C route          ',worst,TOL,nfail)
    deallocate(synt,dp,t,onset)
  endif

  ! a failure after the locate: an unresolved t0 is refused by the plan, so
  ! the location is complete and the axis is not
  call gf_extract(db,src,-1.d0,0,synt,dp,ierr,loc=locx,tax=taxx)
  call report_true('   an unresolved t0 is refused',ierr /= GF_OK,nfail)
  call report_true('   the location survives the failure',locx%ielem > 0,nfail)
  call report_true('   and the axis reports none',taxx%nt == 0,nfail)

  ! a failure at the locate: nothing located, so nothing to report
  far = src
  far%depth = 1.d6
  call gf_extract(db,far,t0,0,synt,dp,ierr,loc=locx)
  call report_true('   a source outside the mesh is refused',ierr /= GF_OK,nfail)
  call report_true('   with no element to report',locx%ielem < 1,nfail)

  deallocate(cseis,cdp,ct,consetd)

  !--- a force source, if the example ships one ---

  if (len_trim(forcepath) > 0) then
    write(*,*) '5. a force source, both routes'
    call test_force(db,ch,forcepath,nfail)
  endif

  cerr = c_gf3d_close(ch)
  call report_true('   gf3d_close',cerr == GF_OK,nfail)

  !--- a second handle on the same database, through the Fortran route ---
  !
  ! The second open is of the same directory on purpose. Ownership of the one
  ! process-wide kd-tree is by open id, not by path, so two handles on one
  ! database exercise it exactly as two databases would -- and this runs
  ! wherever a single database exists, which a second example does not: both
  ! shipped examples are gitignored, so a check needing two of them never ran
  ! anywhere but a machine that had built both.
  !
  ! The rebuild is asserted through gf_locate_tree_owner() and not inferred
  ! from the result. Two handles on one directory hold identical centroids,
  ! so a tree that was never rebuilt would still give the right answer: only
  ! the owner says whether the mechanism worked.

  write(*,*) '6. a second handle on the same database'

  call gf_open(trim(dbpath),db2,ierr,check_completion=.false.)
  call report_true('   a second handle opens',ierr == GF_OK,nfail)
  call report_true('   with a different open id',db2%open_id /= db%open_id,nfail)

  if (ierr == GF_OK) then
    call gf_locate_source(db2,src%latitude,src%longitude,src%depth,loc2,ierr)
    call report_true('   the second handle locates',ierr == GF_OK,nfail)
    call report_true('   and the tree is now its own', &
                     gf_locate_tree_owner() == db2%open_id,nfail)

    call get_seismograms(db,src,t0,synt_again,ierr,t_again)
    call report_true('   the first still extracts',ierr == GF_OK,nfail)
    call report_true('   taking the tree back', &
                     gf_locate_tree_owner() == db%open_id,nfail)

    if (ierr == GF_OK) then
      worst = 0.d0
      do it = 1,size(t_again)
        do ic = 1,GF_NCOMP
          do ista = 1,db%nstations
            worst = max(worst,abs(synt_again(ista,ic,it) - synt_ref(ista,ic,it)))
          enddo
        enddo
      enddo
      call report('   and returns what it did before',worst, &
                  1.d-15*maxval(abs(synt_ref)),nfail)
      deallocate(synt_again,t_again)
    endif

    call gf_close(db2)
  endif

  deallocate(synt_ref)

  !--- 7. a source that leaves its element ---
  !
  ! Extraction from a neighbouring element must work exactly as from the
  ! first: the element index, the kd-tree and the read path are all shared
  ! state, and a source that moves is the whole point of the library.
  !
  ! The crossing distance is found by walking rather than hard-coded, so
  ! this runs on the fixture and on either shipped example. A step that
  ! lands between elements is skipped -- a database covers the region its
  ! source needed, not the whole planet.

  write(*,*) '7. a source in a second element'

  call gf_locate_source(db,src%latitude,src%longitude,src%depth,loc_a,ierr)
  call report_true('   the original position locates',ierr == GF_OK,nfail)

  far = src
  crossed = .false.
  do iw = 1,WALK_MAX
    far%latitude = src%latitude + dble(iw)*WALK_STEP
    call gf_locate_source(db,far%latitude,far%longitude,far%depth,loc_b,ierr)
    if (ierr /= GF_OK) cycle
    if (loc_b%morton_hex /= loc_a%morton_hex) then
      crossed = .true.
      exit
    endif
  enddo

  call report_true('   a second element is reachable    ',crossed,nfail)

  if (crossed) then
    call report_true('   its Morton code differs          ', &
                     loc_b%morton_hex /= loc_a%morton_hex,nfail)
    call report_true('   its element index differs        ', &
                     loc_b%ielem /= loc_a%ielem,nfail)

    call get_seismograms(db,far,t0,synt_far,ierr)
    call report_true('   extraction from it returns GF_OK ',ierr == GF_OK,nfail)

    if (ierr == GF_OK) then
      call report_true('   its traces are finite            ', &
                       all(synt_far == synt_far) .and. &
                       maxval(abs(synt_far)) < huge(1.d0),nfail)
      call report_true('   and not identically zero         ', &
                       maxval(abs(synt_far)) > 0.d0,nfail)
      deallocate(synt_far)
    endif

    ! moving back must recover the first element: nothing about the crossing
    ! may leave the index or the tree wedged
    call gf_locate_source(db,src%latitude,src%longitude,src%depth,loc_back,ierr)
    call report_true('   moving back recovers element A   ', &
                     ierr == GF_OK .and. loc_back%ielem == loc_a%ielem,nfail)
  endif

  ! the pair a long-lived caller must use: gf_close alone leaves the
  ! kd-tree allocated
  call gf_release(db)
  call report_true('   gf_release',.not. db%is_open,nfail)

  write(*,*)
  if (nfail == 0) then
    write(*,*) 'test_gf3d_ext: all assertions passed'
  else
    write(*,*) 'test_gf3d_ext: ',nfail,' assertion(s) FAILED'
    stop 1
  endif
  write(*,*)

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine report_true(name,cond,nfail)

  implicit none

  character(len=*), intent(in) :: name
  logical, intent(in) :: cond
  integer, intent(inout) :: nfail

  if (cond) then
    write(*,'(a,a)') '  ok   ',name
  else
    write(*,'(a,a)') '  FAIL ',name
    nfail = nfail + 1
  endif

  end subroutine report_true

!
!-------------------------------------------------------------------------------------------------
!

  subroutine report(name,err,tol,nfail)

  implicit none

  character(len=*), intent(in) :: name
  double precision, intent(in) :: err,tol
  integer, intent(inout) :: nfail

  if (err <= tol) then
    write(*,'(a,a,a,es12.5,a,es12.5,a)') '  ok   ',name,'   error = ',err,'  (tol ',tol,')'
  else
    write(*,'(a,a,a,es12.5,a,es12.5,a)') '  FAIL ',name,'   error = ',err,'  (tol ',tol,')'
    nfail = nfail + 1
  endif

  end subroutine report

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compare_seis(nsta,nt,synt,cseis,nfail)

! the Fortran (nsta,3,nt) against the C [ista][icomp][it]

  implicit none

  integer, intent(in) :: nsta,nt
  double precision, dimension(nsta,GF_NCOMP,nt), intent(in) :: synt
  real(c_double), dimension(nt,GF_NCOMP,nsta), intent(in) :: cseis
  integer, intent(inout) :: nfail

  ! local parameters
  integer :: ista,icomp,it
  double precision :: worst,peak

  worst = 0.d0
  peak = 0.d0
  do ista = 1,nsta
    do icomp = 1,GF_NCOMP
      do it = 1,nt
        peak = max(peak,abs(synt(ista,icomp,it)))
        worst = max(worst,abs(synt(ista,icomp,it) - cseis(it,icomp,ista)))
      enddo
    enddo
  enddo

  call report('   seismograms   ',worst/max(peak,1.d-300),TOL,nfail)

  end subroutine compare_seis

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compare_time(nt,t,ct,nfail)

  implicit none

  integer, intent(in) :: nt
  double precision, dimension(nt), intent(in) :: t
  real(c_double), dimension(nt), intent(in) :: ct
  integer, intent(inout) :: nfail

  ! local parameters
  integer :: it
  double precision :: worst

  worst = 0.d0
  do it = 1,nt
    worst = max(worst,abs(t(it) - ct(it)))
  enddo

  call report('   time axis     ',worst,1.d-12,nfail)

  end subroutine compare_time

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compare_dp(nsta,nt,dp,cdp,nfail)

! the Fortran (ndp,nsta,3,nt) against the C [ista][ip][icomp][it]
!
! Reported per partial, because the ten of them carry four different units
! and a mistake in the packing would show in one slot, not all ten.

  implicit none

  integer, intent(in) :: nsta,nt
  double precision, dimension(GF_NDP_LOC,nsta,GF_NCOMP,nt), intent(in) :: dp
  real(c_double), dimension(nt,GF_NCOMP,GF_NDP_LOC,nsta), intent(in) :: cdp
  integer, intent(inout) :: nfail

  ! local parameters
  integer :: ista,icomp,it,ip
  double precision :: worst,peak
  character(len=18) :: label

  do ip = 1,GF_NDP_LOC
    worst = 0.d0
    peak = 0.d0
    do ista = 1,nsta
      do icomp = 1,GF_NCOMP
        do it = 1,nt
          peak = max(peak,abs(dp(ip,ista,icomp,it)))
          worst = max(worst,abs(dp(ip,ista,icomp,it) - cdp(it,icomp,ip,ista)))
        enddo
      enddo
    enddo
    write(label,'(a,a,a)') '   dp ',GF_DP_NAME(ip),'        '
    call report(label,worst/max(peak,1.d-300),TOL,nfail)
  enddo

  end subroutine compare_dp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_force(db,ch,forcepath,nfail)

! the force path, both routes, seismograms only

  implicit none

  type(t_gfdb), intent(inout) :: db
  integer(c_int), intent(in) :: ch
  character(len=*), intent(in) :: forcepath
  integer, intent(inout) :: nfail

  ! local parameters
  type(t_gf_source) :: fsrc
  double precision, dimension(:,:,:), allocatable :: fsynt
  double precision, dimension(:), allocatable :: ft
  double precision :: ft0
  type(gf3d_source_t) :: cf
  type(gf3d_plan_t) :: fplan
  real(c_double), dimension(:), allocatable :: fseis,fct,fonset
  integer(c_int) :: cerr
  integer :: ierr

  call gf_read_force_source(db,trim(forcepath),fsrc,ierr)
  call report_true('   gf_read_force_source',ierr == GF_OK,nfail)
  if (ierr /= GF_OK) return

  call gf_default_t0(fsrc,ft0,ierr)
  call report_true('   gf_default_t0, force',ierr == GF_OK,nfail)
  if (ierr /= GF_OK) return

  call get_seismograms(db,fsrc,ft0,fsynt,ierr,ft)
  call report_true('   get_seismograms, force',ierr == GF_OK,nfail)
  if (ierr /= GF_OK) return

  ! The C side takes the file's own numbers, in the file's own units, so
  ! they are read from the file rather than recovered from the reader's
  ! output: factor_force_source comes back already divided by scaleF, and
  ! an external caller has no way to name that constant.
  cf%source_type = GF_SRC_FORCE
  cf%moment(1:6) = 0.d0
  call read_forcesolution(forcepath,cf,ierr)
  call report_true('   FORCESOLUTION parsed',ierr == 0,nfail)
  if (ierr /= 0) return

  cerr = c_gf3d_get_plan(ch,cf,-1.0_c_double,fplan)
  call report_true('   gf3d_get_plan, force',cerr == GF_OK,nfail)
  if (cerr /= GF_OK) return

  call report_true('   the force plan lengths agree',int(fplan%nt) == size(ft),nfail)

  allocate(fseis(db%nstations*GF_NCOMP*int(fplan%nt)),fct(int(fplan%nt)), &
           fonset(db%nstations))

  cerr = c_gf3d_seismograms(ch,cf,-1.0_c_double,fplan%nt,fseis,fct,fonset,c_null_ptr)
  call report_true('   gf3d_seismograms, force',cerr == GF_OK,nfail)

  if (cerr == GF_OK) call compare_seis(db%nstations,int(fplan%nt),fsynt,fseis,nfail)

  deallocate(fseis,fct,fonset)
  deallocate(fsynt,ft)

  end subroutine test_force

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_forcesolution(filename,cf,ierr)

! the eleven lines of a FORCESOLUTION, in the file's own units
!
! Every line is "name: value", so the value is what follows the last colon;
! list-directed input accepts the Fortran `d` exponent the factor is
! usually written with. This is what a caller outside the tree has to do,
! and it is the reason the C API takes values rather than a path: reading
! the file with get_force() would put its stop statements back in reach.

  implicit none

  character(len=*), intent(in) :: filename
  type(gf3d_source_t), intent(inout) :: cf
  integer, intent(out) :: ierr

  ! local parameters
  integer, parameter :: IIN_F = 93
  character(len=512) :: line
  integer :: iline,icolon,ios,istf
  double precision :: v

  ierr = 1

  open(unit=IIN_F,file=trim(filename),status='old',action='read',iostat=ios)
  if (ios /= 0) return

  do iline = 1,11
    read(IIN_F,'(a)',iostat=ios) line
    if (ios /= 0) then
      close(IIN_F)
      return
    endif

    ! the first line is 'FORCE <label>', which carries no value
    if (iline == 1) cycle

    icolon = index(line,':',back=.true.)
    if (icolon < 1) then
      close(IIN_F)
      return
    endif

    if (iline == 7) then
      read(line(icolon+1:),*,iostat=ios) istf
    else
      read(line(icolon+1:),*,iostat=ios) v
    endif
    if (ios /= 0) then
      close(IIN_F)
      return
    endif

    select case (iline)
    case (2) ; cf%time_shift = v
    case (3) ; cf%hdur = v
    case (4) ; cf%latitude = v
    case (5) ; cf%longitude = v
    case (6) ; cf%depth_km = v
    case (7) ; cf%force_stf = int(istf,kind=c_int)
    case (8) ; cf%force_factor = v
    case (9) ; cf%force_dir(1) = v
    case (10) ; cf%force_dir(2) = v
    case (11) ; cf%force_dir(3) = v
    end select
  enddo

  close(IIN_F)

  ierr = 0

  end subroutine read_forcesolution

  end program test_gf3d_ext
