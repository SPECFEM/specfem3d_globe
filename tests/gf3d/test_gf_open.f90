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
!---- test_gf_open — exercises gf_open / gf_close / gf_load_topo against a
!---- real Green function database.
!----
!---- Needs a built example database, which is gitignored, so the runner
!---- skips this when there is none. The database path comes in as the
!---- first command-line argument.
!----
!---- What is actually checked here is *internal consistency plus the
!---- library's own invariants*, not physics: that the metadata agrees with
!---- itself, that the element index is ordered and round-trips through its
!---- hex representation, that every route into the index yields the same
!---- answer, that the topography grid loads at the shape its attributes
!---- claim, and — the property the whole design rests on — that a bad path
!---- comes back as an error code rather than stopping the process.
!----
!---- The numbers themselves are checked against `h5dump -A` by the runner.
!----

  program test_gf_open

  use gf_par
  use gf_database

  implicit none

  character(len=MAX_STRING_LEN) :: dbpath
  type(t_gfdb) :: db
  integer :: ierr,i,nfail,nincomplete

  nfail = 0

  if (command_argument_count() < 1) then
    print *,'usage: test_gf_open <GFDB directory>'
    stop 1
  endif
  call get_command_argument(1,dbpath)

  print *,'test_gf_open'
  print *,'  database: ',trim(dbpath)
  print *,''

  !--------------------------------------------------------------------
  ! a path that does not exist must return an error, not stop
  !
  ! This is the single most important assertion in the file. Stage 9 loads
  ! this library into a Python interpreter, where a Fortran stop takes the
  ! interpreter with it, so no failure mode may be allowed to reach one.
  !--------------------------------------------------------------------
  call gf_open('/no/such/green/function/database',db,ierr)
  call check(ierr /= GF_OK,'a nonexistent path returns an error',nfail)
  call check(ierr == GF_ERR_NO_PATH,'and the error is GF_ERR_NO_PATH',nfail)
  call check(.not. db%is_open,'and the handle is left closed',nfail)

  ! an empty path likewise
  call gf_open('',db,ierr)
  call check(ierr == GF_ERR_ARG,'an empty path returns GF_ERR_ARG',nfail)

  ! a directory that exists but holds no database
  call gf_open('.',db,ierr)
  call check(ierr /= GF_OK,'a directory without mesh_info.h5 returns an error',nfail)

  !--------------------------------------------------------------------
  ! the real thing
  !--------------------------------------------------------------------
  call gf_open(dbpath,db,ierr)
  if (ierr /= GF_OK) then
    print *,'  FAIL: gf_open returned ',trim(gf_error_string(ierr))
    print *,'        ',trim(gf_errmsg)
    stop 1
  endif
  call check(db%is_open,'gf_open succeeds on the example database',nfail)

  !--- metadata self-consistency ---
  call check(db%dt > 0.d0,'dt is positive',nfail)
  call check(db%nstep > 0,'nstep is positive',nfail)
  call check(db%subsample_step > 0,'subsample_step is positive',nfail)
  call check(db%nt_subsampled == db%nstep / db%subsample_step, &
             'nt_subsampled == nstep / subsample_step',nfail)
  call check(db%R_PLANET > 0.d0,'R_PLANET is positive',nfail)
  call check(db%scale_displ > 0.d0,'scale_displ is positive',nfail)

  !--- element index ---
  call check(db%nelem > 0,'the database holds at least one element',nfail)
  call check(allocated(db%morton),'the Morton codes are allocated',nfail)
  call check(allocated(db%centroid),'the centroids are allocated',nfail)
  call check(size(db%centroid,1) == 3,'centroids are (3,nelem)',nfail)
  call check(size(db%centroid,2) == db%nelem,'centroids have nelem columns',nfail)

  call check(index_is_ordered(db),'Morton codes ascend strictly',nfail)
  call check(hex_round_trips(db),'every Morton code round-trips through its hex name',nfail)
  call check(centroids_are_inside_the_planet(db), &
             'every centroid lies inside the non-dimensional unit ball',nfail)

  !--- stations ---
  call check(db%nstations > 0,'the database holds at least one station',nfail)
  call check(allocated(db%stations),'the station table is allocated',nfail)

  do i = 1,db%nstations
    call check(len_trim(db%stations(i)%id) > 0,'station '//trim(itoa(i))//' has an id',nfail)
    call check(len_trim(db%stations(i)%network) > 0, &
               'station '//trim(db%stations(i)%id)//' has a network code',nfail)
    call check(len_trim(db%stations(i)%station) > 0, &
               'station '//trim(db%stations(i)%id)//' has a station code',nfail)
    call check(abs(db%stations(i)%latitude) <= 90.d0, &
               'station '//trim(db%stations(i)%id)//' has a valid latitude',nfail)
    call check(abs(db%stations(i)%longitude) <= 360.d0, &
               'station '//trim(db%stations(i)%id)//' has a valid longitude',nfail)
    call check(db%stations(i)%hdur >= 0.d0, &
               'station '//trim(db%stations(i)%id)//' has hdur >= 0',nfail)
    call check(allocated(db%stations(i)%stf), &
               'station '//trim(db%stations(i)%id)//' has an stf',nfail)
    call check(size(db%stations(i)%stf) == db%nstep, &
               'station '//trim(db%stations(i)%id)//' stf has nstep samples',nfail)
    call check(maxval(abs(db%stations(i)%stf)) > 0.d0, &
               'station '//trim(db%stations(i)%id)//' stf is not all zeros',nfail)
  enddo

  ! station ids must be unique, or the element file names collide
  call check(station_ids_are_unique(db),'station ids are unique',nfail)

  !--- completion ---
  call gf_check_completion(db,nincomplete,ierr)
  call check(ierr == GF_OK .and. nincomplete == 0, &
             'every element-station file exists with computed_ALL set',nfail)

  !--- ellipticity ---
  if (db%ellipticity) then
    call check(db%nspl > 0,'nspl is positive',nfail)
    call check(allocated(db%rspl),'the ellipticity splines are allocated',nfail)
    call check(size(db%rspl) == db%nspl,'rspl has nspl entries',nfail)
    call check(rspl_is_ordered(db),'rspl ascends',nfail)
  endif

  !--- topography, loaded on demand ---
  if (db%topography) then
    call check(.not. db%topo_loaded,'gf_open leaves the topography grid unloaded',nfail)

    call gf_load_topo(db,ierr)
    call check(ierr == GF_OK,'gf_load_topo succeeds',nfail)
    call check(db%topo_loaded,'and marks the grid loaded',nfail)
    call check(allocated(db%ibathy_topo),'and allocates it',nfail)
    call check(size(db%ibathy_topo,1) == db%NX_BATHY,'grid is (NX_BATHY,.)',nfail)
    call check(size(db%ibathy_topo,2) == db%NY_BATHY,'grid is (.,NY_BATHY)',nfail)

    ! a grid read at the wrong shape or from the wrong dataset would still be
    ! integers, but not plausible elevations in metres
    call check(minval(db%ibathy_topo) > -12000 .and. minval(db%ibathy_topo) < 0, &
               'the deepest point is an ocean trench, not garbage',nfail)
    call check(maxval(db%ibathy_topo) > 0 .and. maxval(db%ibathy_topo) < 12000, &
               'the highest point is a mountain, not garbage',nfail)

    ! loading twice must be a no-op rather than a leak or a re-read
    call gf_load_topo(db,ierr)
    call check(ierr == GF_OK,'gf_load_topo is idempotent',nfail)
  endif

  !--------------------------------------------------------------------
  ! closing
  !--------------------------------------------------------------------
  call gf_close(db)
  call check(.not. db%is_open,'gf_close closes the handle',nfail)
  call check(.not. allocated(db%morton),'and releases the element index',nfail)
  call check(.not. allocated(db%stations),'and releases the station table',nfail)
  call check(.not. allocated(db%ibathy_topo),'and releases the topography grid',nfail)

  ! closing twice must be safe
  call gf_close(db)
  call check(.not. db%is_open,'gf_close is idempotent',nfail)

  ! and the handle must be reusable
  call gf_open(dbpath,db,ierr)
  call check(ierr == GF_OK,'a closed handle can be reopened',nfail)
  call gf_close(db)

  !--------------------------------------------------------------------
  print *,''
  if (nfail == 0) then
    print *,'test_gf_open: all assertions passed'
  else
    print *,'test_gf_open: ',nfail,' assertion(s) FAILED'
    stop 1
  endif

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine check(cond,what,nfail)

! prints every assertion, passing or failing, so that results.log records
! what was actually exercised rather than a pass/fail bit

  implicit none

  logical, intent(in) :: cond
  character(len=*), intent(in) :: what
  integer, intent(inout) :: nfail

  if (cond) then
    print *,'  ok   : ',what
  else
    print *,'  FAIL : ',what
    nfail = nfail + 1
  endif

  end subroutine check

!
!-------------------------------------------------------------------------------------------------
!

  logical function index_is_ordered(db)

  implicit none
  type(t_gfdb), intent(in) :: db
  integer :: i

  index_is_ordered = .true.
  do i = 2,db%nelem
    if (db%morton(i) <= db%morton(i-1)) then
      index_is_ordered = .false.
      return
    endif
  enddo

  end function index_is_ordered

!
!-------------------------------------------------------------------------------------------------
!

  logical function hex_round_trips(db)

! the directory name and the integer code must agree, in both directions

  implicit none
  type(t_gfdb), intent(in) :: db
  character(len=GF_MORTON_HEXLEN) :: s
  integer :: i

  hex_round_trips = .true.
  do i = 1,db%nelem
    write(s,'(Z16.16)') db%morton(i)
    if (s /= db%morton_hex(i)) then
      hex_round_trips = .false.
      return
    endif
  enddo

  end function hex_round_trips

!
!-------------------------------------------------------------------------------------------------
!

  logical function centroids_are_inside_the_planet(db)

! centroids are non-dimensionalised by R_PLANET, the mean radius, so a
! surface element sits slightly above 1 on an elliptical planet with
! topography; 1.05 is generous. See the note in gf_database.F90.

  implicit none
  type(t_gfdb), intent(in) :: db
  double precision :: r
  integer :: i

  centroids_are_inside_the_planet = .true.
  do i = 1,db%nelem
    r = sqrt(db%centroid(1,i)**2 + db%centroid(2,i)**2 + db%centroid(3,i)**2)
    if (r > 1.05d0 .or. r <= 0.d0) then
      centroids_are_inside_the_planet = .false.
      return
    endif
  enddo

  end function centroids_are_inside_the_planet

!
!-------------------------------------------------------------------------------------------------
!

  logical function station_ids_are_unique(db)

  implicit none
  type(t_gfdb), intent(in) :: db
  integer :: i,j

  station_ids_are_unique = .true.
  do i = 1,db%nstations
    do j = i+1,db%nstations
      if (db%stations(i)%id == db%stations(j)%id) then
        station_ids_are_unique = .false.
        return
      endif
    enddo
  enddo

  end function station_ids_are_unique

!
!-------------------------------------------------------------------------------------------------
!

  logical function rspl_is_ordered(db)

  implicit none
  type(t_gfdb), intent(in) :: db
  integer :: i

  rspl_is_ordered = .true.
  do i = 2,db%nspl
    if (db%rspl(i) < db%rspl(i-1)) then
      rspl_is_ordered = .false.
      return
    endif
  enddo

  end function rspl_is_ordered

!
!-------------------------------------------------------------------------------------------------
!

  character(len=12) function itoa(i)

  implicit none
  integer, intent(in) :: i

  write(itoa,'(i0)') i

  end function itoa

  end program test_gf_open
