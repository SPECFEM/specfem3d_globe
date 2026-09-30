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
!---- test_gf_cache -- what a handle keeps between extractions
!----
!---- The handle owns a cache object from gf_open to gf_close, and counts
!---- what its extractions load: misses (element blocks read from disk),
!---- hits, evictions, the elements held, and files_read (element files
!---- whose data it read). This file pins those counters -- they are what
!---- makes the cache's behaviour testable at all -- and the object's
!---- lifetime. With max_elements > 0 it also pins that the cache changes
!---- no number, that it evicts the least recently used element, and that
!---- a fill that fails part-way leaves nothing behind.
!----
!---- Elements A, B and C are found by walking north from the CMTSOLUTION,
!---- as test_gf3d_ext does, so this runs on the fixture (three elements, 8
!---- degrees apart) and on a real $GF3D_TEST_GFDB alike.
!----

  program test_gf_cache

  use gf3d, only: t_gfdb,t_gf_source,t_gf_location, &
                  gf_open,gf_close,gf_cache_stats,gf_read_source, &
                  gf_locate_source,gf_locate_release,get_seismograms,gf_extract, &
                  gf_default_t0,gf_element_bytes, &
                  gf_errmsg,GF_OK,GF_ERR_ARG

  use gf_manufactured, only: gf_report_true

  implicit none

  ! the northward walk that finds the next element
  double precision, parameter :: WALK_STEP = 0.25d0
  integer, parameter :: WALK_MAX = 80

  character(len=512) :: dbpath,cmtpath,forcepath,broken
  type(t_gfdb) :: db,dbc,never_opened
  type(t_gf_source) :: src,force,s
  type(t_gf_location) :: loc,loc_c
  double precision, dimension(:,:,:), allocatable :: synt,synt0
  double precision, dimension(:,:,:,:), allocatable :: dp,dp0
  double precision :: t0,t0_force
  ! latitudes of elements A, B, C (1..3), and whether each was found
  double precision, dimension(3) :: lat
  integer, dimension(3) :: ielem
  character(len=16) :: hex_b
  logical :: have_b,have_c,staged
  integer(kind=8) :: hits,misses,evictions,files_read
  integer(kind=8) :: misses0,files0,nread_locate,nread_extract
  integer :: n_cached,ierr,nfail,k,nbad
  integer, dimension(4) :: seq

  nfail = 0

  if (command_argument_count() < 3) then
    write(*,'(a)') 'usage: test_gf_cache <GFDB> <CMTSOLUTION> <FORCESOLUTION>'
    stop 1
  endif
  call get_command_argument(1,dbpath)
  call get_command_argument(2,cmtpath)
  call get_command_argument(3,forcepath)

  write(*,'(a)') ''
  write(*,'(a)') 'test_gf_cache'
  write(*,'(a)') ''
  write(*,'(a,a)') '  database = ',trim(dbpath)
  write(*,'(a,a)') '  source   = ',trim(cmtpath)
  write(*,'(a)') ''

  !--------------------------------------------------------------------
  ! 1. the cache belongs to an open handle
  !--------------------------------------------------------------------

  write(*,'(a)') '1. lifetime'

  call gf_cache_stats(never_opened,hits,misses,evictions,n_cached,files_read,ierr)
  call gf_report_true('a handle never opened has no statistics',ierr == GF_ERR_ARG,nfail)
  call gf_report_true('and no cache                           ', &
                      .not. associated(never_opened%cache),nfail)

  call open_or_die(db)
  call gf_report_true('an open handle has a cache             ',associated(db%cache),nfail)

  call gf_cache_stats(db,hits,misses,evictions,n_cached,files_read,ierr)
  call gf_report_true('its statistics can be read             ',ierr == GF_OK,nfail)
  call expect('hits after open      ',hits,0_8,nfail)
  call expect('misses after open    ',misses,0_8,nfail)
  call expect('evictions after open ',evictions,0_8,nfail)
  call expect('n_cached after open  ',int(n_cached,8),0_8,nfail)
  call expect('files_read after open',files_read,0_8,nfail)

  call gf_read_source(db,cmtpath,src,ierr)
  if (ierr /= GF_OK) call die('could not read the CMTSOLUTION')
  call gf_default_t0(src,t0,ierr)
  if (ierr /= GF_OK) call die('no default start time for this source')

  !--------------------------------------------------------------------
  ! 2. what the counters count, with nothing kept
  !
  ! A locate reads coordinates and loads no element block. An extraction
  ! is a locate plus one element block, which is one displacement file per
  ! station. With no cache, a second identical extraction costs the same
  ! again. The locate's own count depends on how many candidates it tried,
  ! so it is measured, then required to be the same inside the extraction.
  !--------------------------------------------------------------------

  write(*,'(a)') ''
  write(*,'(a)') '2. counters with no cache'

  call snapshot(db,misses0,files0)
  call gf_locate_source(db,src%latitude,src%longitude,src%depth,loc,ierr)
  if (ierr /= GF_OK) call die('could not locate the source')
  call gf_cache_stats(db,hits,misses,evictions,n_cached,files_read,ierr)
  nread_locate = files_read - files0
  write(*,'(a,i0)') '     coordinate files a locate reads = ',nread_locate
  call gf_report_true('a locate reads coordinates             ',nread_locate >= 1,nfail)
  call expect('misses from a locate ',misses - misses0,0_8,nfail)

  nread_extract = nread_locate + int(db%nstations,8)

  call snapshot(db,misses0,files0)
  call get_seismograms(db,src,t0,synt,ierr)
  if (ierr /= GF_OK) call die('the first extraction failed')
  call gf_cache_stats(db,hits,misses,evictions,n_cached,files_read,ierr)
  call expect('misses, extraction 1 ',misses - misses0,1_8,nfail)
  call expect('files, extraction 1  ',files_read - files0,nread_extract,nfail)

  call snapshot(db,misses0,files0)
  call get_seismograms(db,src,t0,synt,ierr)
  if (ierr /= GF_OK) call die('the second extraction failed')
  call gf_cache_stats(db,hits,misses,evictions,n_cached,files_read,ierr)
  call expect('misses, extraction 2 ',misses - misses0,1_8,nfail)
  call expect('files, extraction 2  ',files_read - files0,nread_extract,nfail)
  call expect('hits with no cache   ',hits,0_8,nfail)
  call expect('n_cached, no cache   ',int(n_cached,8),0_8,nfail)
  call expect('evictions, no cache  ',evictions,0_8,nfail)

  !--------------------------------------------------------------------
  ! 3. three elements to move between
  !--------------------------------------------------------------------

  write(*,'(a)') ''
  write(*,'(a)') '3. elements A, B, C'

  lat(:) = src%latitude
  ielem(:) = 0
  ielem(1) = loc%ielem
  call find_next(lat(1),ielem(1),lat(2),ielem(2),have_b)
  have_c = .false.
  if (have_b) call find_next(lat(2),ielem(2),lat(3),ielem(3),have_c)
  ! a walk north from B could in principle come back to A's column
  if (have_c) have_c = (ielem(3) /= ielem(1))

  write(*,'(a,i0,a,f9.3)') '     A: element ',ielem(1),' at latitude ',lat(1)
  if (have_b) write(*,'(a,i0,a,f9.3)') '     B: element ',ielem(2),' at latitude ',lat(2)
  if (have_c) write(*,'(a,i0,a,f9.3)') '     C: element ',ielem(3),' at latitude ',lat(3)
  write(*,'(a,i0,a)') '     one element takes ',gf_element_bytes(db),' bytes'

  call gf_report_true('a second element B is reachable        ',have_b,nfail)
  if (.not. have_c) then
    write(*,'(a)') '     no third element within reach: the sequences that need C'
    write(*,'(a)') '     are not exercised on this database'
  endif

  if (.not. have_b) call die('everything below moves between elements')

  !--------------------------------------------------------------------
  ! 4. the cache changes no number
  !
  ! A handle keeping two elements against the one keeping none, for
  ! seismograms and all ten partials, at A, B, A, B: two misses, then a hit
  ! on each slot. Bitwise, not to a tolerance: both handles run the same
  ! gf_seis_station call on the same bytes, the cache's copy being exactly
  ! what gf_read_element_displ wrote into it. And B's traces must differ
  ! from A's, or a cache that served the wrong element would pass too.
  !--------------------------------------------------------------------

  write(*,'(a)') ''
  write(*,'(a)') '4. max_elements = 2 against 0, bitwise'

  call gf_open(dbpath,dbc,ierr,check_completion=.false.,max_elements=2)
  if (ierr /= GF_OK) call die('could not open the database with max_elements = 2')

  seq = (/ 1, 2, 1, 2 /)
  do k = 1,4
    s = src
    s%latitude = lat(seq(k))
    call gf_extract(db,s,t0,2,synt0,dp0,ierr)
    if (ierr /= GF_OK) call die('extraction without a cache failed')
    call gf_extract(dbc,s,t0,2,synt,dp,ierr)
    if (ierr /= GF_OK) call die('extraction with a cache failed')
    nbad = count(synt /= synt0) + count(dp /= dp0)
    call expect('CMT + 10 partials, '//label(seq(k))//' (call '//trim(itoa8(int(k,8)))// &
                '), differing samples',int(nbad,8),0_8,nfail)
    if (k == 2) then
      ! synt0 holds B's traces: A's, extracted again, must differ
      s%latitude = lat(1)
      call get_seismograms(db,s,t0,synt,ierr)
      if (ierr /= GF_OK) call die('extraction at A failed')
      call gf_report_true('B''s traces differ from A''s              ',any(synt /= synt0),nfail)
    endif
  enddo

  call gf_cache_stats(dbc,hits,misses,evictions,n_cached,files_read,ierr)
  call expect('misses over A B A B  ',misses,2_8,nfail)
  call expect('hits over A B A B    ',hits,2_8,nfail)

  ! the force route goes through the interpolation rather than the strain
  call gf_read_source(db,forcepath,force,ierr)
  if (ierr /= GF_OK) call die('could not read the FORCESOLUTION')
  call gf_default_t0(force,t0_force,ierr)
  if (ierr /= GF_OK) call die('no default start time for the force')
  do k = 1,2
    s = force
    s%latitude = lat(k)
    call get_seismograms(db,s,t0_force,synt0,ierr)
    if (ierr /= GF_OK) call die('force extraction without a cache failed')
    call get_seismograms(dbc,s,t0_force,synt,ierr)
    if (ierr /= GF_OK) call die('force extraction with a cache failed')
    call expect('force, '//label(k)//', differing samples', &
                int(count(synt /= synt0),8),0_8,nfail)
  enddo

  call gf_close(dbc)

  !--------------------------------------------------------------------
  ! 5. the least recently used element is the one that goes
  !
  ! Each sequence on a fresh handle. The last is the one that tells LRU
  ! from FIFO: after A B A, B is the older in use and A the older in
  ! insertion, so C evicts B and the final A is a hit (FIFO: 4 misses,
  ! 1 hit, 2 evictions).
  !--------------------------------------------------------------------

  write(*,'(a)') ''
  write(*,'(a)') '5. eviction order'

  call sequence(1,(/ 1, 2, 1 /),'N=1 A B A    ',3_8,0_8,2_8,1)
  call sequence(2,(/ 1, 2, 1 /),'N=2 A B A    ',2_8,1_8,0_8,2)
  if (have_c) then
    call sequence(2,(/ 1, 2, 3, 1 /),'N=2 A B C A  ',4_8,0_8,2_8,2)
    call sequence(2,(/ 1, 2, 1, 3, 1 /),'N=2 A B A C A',3_8,2_8,1_8,2)
  endif

  ! What a caching handle reads. The first extraction at A reads what an
  ! uncached one does; returning to A reads nothing at all -- not the
  ! displacement, and not the coordinates the locate needs either, which
  ! the handle kept from the first time. Nor does a bare locate there.
  ! Through C when there is one: the coordinate store starts with room for
  ! two elements, so the third makes it grow, and A's coordinates are then
  ! read back from the grown copy.
  call gf_open(dbpath,dbc,ierr,check_completion=.false.,max_elements=3)
  if (ierr /= GF_OK) call die('could not open the database with max_elements = 3')
  call snapshot(dbc,misses0,files0)
  call extract_at(dbc,1)
  call gf_cache_stats(dbc,hits,misses,evictions,n_cached,files_read,ierr)
  call expect('files, first at A    ',files_read - files0,nread_extract,nfail)
  call extract_at(dbc,2)
  if (have_c) call extract_at(dbc,3)
  call snapshot(dbc,misses0,files0)
  call extract_at(dbc,1)
  call gf_cache_stats(dbc,hits,misses,evictions,n_cached,files_read,ierr)
  call expect('files read by a hit  ',files_read - files0,0_8,nfail)
  call snapshot(dbc,misses0,files0)
  call gf_locate_source(dbc,src%latitude,src%longitude,src%depth,loc_c,ierr)
  call gf_cache_stats(dbc,hits,misses,evictions,n_cached,files_read,ierr)
  call expect('files, locate again  ',files_read - files0,0_8,nfail)
  ! and the stored coordinates put the source where a fresh read does
  call gf_report_true('same element, xi, eta, gamma as uncached', &
                      loc_c%ielem == loc%ielem .and. loc_c%xi == loc%xi .and. &
                      loc_c%eta == loc%eta .and. loc_c%gamma == loc%gamma,nfail)
  call gf_close(dbc)

  !--------------------------------------------------------------------
  ! 6. what max_elements accepts
  !--------------------------------------------------------------------

  write(*,'(a)') ''
  write(*,'(a)') '6. max_elements'

  call gf_open(dbpath,dbc,ierr,check_completion=.false.,max_elements=-1)
  call gf_report_true('a negative max_elements is refused      ',ierr == GF_ERR_ARG,nfail)
  call gf_report_true('and leaves the handle closed            ',.not. dbc%is_open,nfail)

  call gf_open(dbpath,dbc,ierr,check_completion=.false.,max_elements=huge(1))
  call gf_report_true('more than nelem opens                   ',ierr == GF_OK,nfail)
  if (ierr == GF_OK) then
    call expect('and keeps at most nelem',int(dbc%cache%capacity,8),int(dbc%nelem,8),nfail)
  endif
  call gf_close(dbc)

  !--------------------------------------------------------------------
  ! 7. a fill that fails part-way leaves nothing behind
  !
  ! A copy of the database, made of symbolic links, whose element B lacks
  ! its last station's displacement: the fill reads every other station
  ! and then fails. B must then fail again rather than come back as a hit
  ! on a half-filled slot, and A must be untouched.
  !--------------------------------------------------------------------

  write(*,'(a)') ''
  write(*,'(a)') '7. a failed fill'

  broken = trim(dbpath)//'.gf3d_test_cache_broken'
  hex_b = db%morton_hex(ielem(2))
  call stage_broken_db(trim(dbpath),trim(broken),hex_b, &
                       trim(db%stations(db%nstations)%id),staged)

  if (.not. staged) then
    write(*,'(a)') '     note: could not stage the broken copy; not exercised'
  else
    call gf_open(trim(broken),dbc,ierr,check_completion=.false.,max_elements=2)
    if (ierr /= GF_OK) call die('could not open the broken copy')

    call extract_at(dbc,1)
    s = src
    s%latitude = lat(2)
    call get_seismograms(dbc,s,t0,synt,ierr)
    call gf_report_true('B''s fill fails                         ',ierr /= GF_OK,nfail)
    call get_seismograms(dbc,s,t0,synt,ierr)
    call gf_report_true('and fails again, not a hit              ',ierr /= GF_OK,nfail)
    call extract_at(dbc,1)

    call gf_cache_stats(dbc,hits,misses,evictions,n_cached,files_read,ierr)
    call expect('misses: A, B, B      ',misses,3_8,nfail)
    call expect('hits: A again        ',hits,1_8,nfail)
    call expect('elements held: A     ',int(n_cached,8),1_8,nfail)

    call gf_close(dbc)
    call execute_command_line('rm -rf '''//trim(broken)//'''')
  endif

  !--------------------------------------------------------------------
  ! 8. closing frees it, and a reopened handle starts again
  !--------------------------------------------------------------------

  write(*,'(a)') ''
  write(*,'(a)') '8. close'

  call gf_locate_release()
  call gf_close(db)
  call gf_report_true('a closed handle has no cache           ',.not. associated(db%cache),nfail)
  call gf_cache_stats(db,hits,misses,evictions,n_cached,files_read,ierr)
  call gf_report_true('and no statistics                      ',ierr == GF_ERR_ARG,nfail)

  call open_or_die(db)
  call gf_cache_stats(db,hits,misses,evictions,n_cached,files_read,ierr)
  call expect('misses after reopen  ',misses,0_8,nfail)
  call expect('files after reopen   ',files_read,0_8,nfail)

  call gf_locate_release()
  call gf_close(db)
  if (allocated(synt)) deallocate(synt)
  if (allocated(synt0)) deallocate(synt0)
  if (allocated(dp)) deallocate(dp)
  if (allocated(dp0)) deallocate(dp0)

  write(*,'(a)') ''
  if (nfail /= 0) then
    write(*,'(a,i0,a)') 'test_gf_cache: ',nfail,' assertion(s) FAILED'
    stop 1
  endif
  write(*,'(a)') 'test_gf_cache: all assertions passed'

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine find_next(lat_from,ielem_from,lat_to,ielem_to,found)

! walks north from lat_from until the source locates in another element;
! a step that lands between elements is skipped

  implicit none
  double precision, intent(in) :: lat_from
  integer, intent(in) :: ielem_from
  double precision, intent(out) :: lat_to
  integer, intent(out) :: ielem_to
  logical, intent(out) :: found

  type(t_gf_location) :: l
  integer :: iw,ier

  found = .false.
  lat_to = lat_from
  ielem_to = 0
  do iw = 1,WALK_MAX
    lat_to = lat_from + dble(iw)*WALK_STEP
    if (lat_to > 90.d0) exit
    call gf_locate_source(db,lat_to,src%longitude,src%depth,l,ier)
    if (ier /= GF_OK) cycle
    if (l%ielem /= ielem_from) then
      ielem_to = l%ielem
      found = .true.
      return
    endif
  enddo

  end subroutine find_next

!
!-------------------------------------------------------------------------------------------------
!

  subroutine extract_at(h,iabc)

! a seismogram-only extraction with the source moved to element A, B or C

  implicit none
  type(t_gfdb), intent(inout) :: h
  integer, intent(in) :: iabc

  type(t_gf_source) :: sm
  double precision, dimension(:,:,:), allocatable :: out
  integer :: ier

  sm = src
  sm%latitude = lat(iabc)
  call get_seismograms(h,sm,t0,out,ier)
  if (ier /= GF_OK) call die('extraction at '//label(iabc)//' failed')

  end subroutine extract_at

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sequence(nkeep,order,name,want_misses,want_hits,want_evictions,want_held)

! one sequence of extractions on a fresh handle, and the counters it ends on

  implicit none
  integer, intent(in) :: nkeep
  integer, dimension(:), intent(in) :: order
  character(len=*), intent(in) :: name
  integer(kind=8), intent(in) :: want_misses,want_hits,want_evictions
  integer, intent(in) :: want_held

  type(t_gfdb) :: h
  integer(kind=8) :: nh,nm,ne,nf
  integer :: nc,i,ier

  call gf_open(dbpath,h,ier,check_completion=.false.,max_elements=nkeep)
  if (ier /= GF_OK) call die('could not open the database for '//name)

  do i = 1,size(order)
    call extract_at(h,order(i))
  enddo

  call gf_cache_stats(h,nh,nm,ne,nc,nf,ier)
  call expect(name//' misses   ',nm,want_misses,nfail)
  call expect(name//' hits     ',nh,want_hits,nfail)
  call expect(name//' evictions',ne,want_evictions,nfail)
  call expect(name//' held     ',int(nc,8),int(want_held,8),nfail)

  call gf_close(h)

  end subroutine sequence

!
!-------------------------------------------------------------------------------------------------
!

  function label(i) result(name)

  implicit none
  integer, intent(in) :: i
  character(len=1) :: name

  character(len=3), parameter :: NAMES = 'ABC'

  name = NAMES(i:i)

  end function label

!
!-------------------------------------------------------------------------------------------------
!

  subroutine stage_broken_db(from,to,hex,last_station,ok)

! a copy of the database, all symbolic links, in which element `hex` lacks
! the last station's displacement file
!
! execute_command_line, as test_gf_open stages its partial database: this
! test has no business writing HDF5. `ok` false means the case is not
! exercised, not that anything failed.

  implicit none
  character(len=*), intent(in) :: from,to,hex,last_station
  logical, intent(out) :: ok

  integer :: stat

  ok = .false.
  ! set before the call: libgfortran reads exitstat's incoming value
  stat = -1
  call execute_command_line( &
    'set -e; src=$(cd '''//from//''' && pwd); dst='''//to//'''; ' // &
    'rm -rf "$dst"; mkdir -p "$dst/elements"; ' // &
    'cp "$src/mesh_info.h5" "$dst/"; ' // &
    'for f in manifest.csv centroids.bin; do [ -e "$src/$f" ] && cp "$src/$f" "$dst/" || true; done; ' // &
    'ln -s "$src/stations" "$dst/stations"; ' // &
    'for d in "$src"/elements/*; do ln -s "$d" "$dst/elements/"; done; ' // &
    'rm "$dst/elements/'//hex//'"; mkdir "$dst/elements/'//hex//'"; ' // &
    'for f in "$src/elements/'//hex//'"/*; do ln -s "$f" "$dst/elements/'//hex//'/"; done; ' // &
    'rm "$dst/elements/'//hex//'/'//last_station//'.h5"', &
    exitstat=stat)
  ok = (stat == 0)

  end subroutine stage_broken_db

!
!-------------------------------------------------------------------------------------------------
!

  subroutine open_or_die(h)

  implicit none
  type(t_gfdb), intent(inout) :: h

  integer :: ier

  call gf_open(dbpath,h,ier,check_completion=.false.)
  if (ier /= GF_OK) call die('could not open the database')

  end subroutine open_or_die

!
!-------------------------------------------------------------------------------------------------
!

  subroutine snapshot(h,m,f)

  implicit none
  type(t_gfdb), intent(in) :: h
  integer(kind=8), intent(out) :: m,f

  integer(kind=8) :: hh,ee
  integer :: nc,ier

  call gf_cache_stats(h,hh,m,ee,nc,f,ier)
  if (ier /= GF_OK) call die('could not read the cache statistics')

  end subroutine snapshot

!
!-------------------------------------------------------------------------------------------------
!

  subroutine expect(name,got,want,nf)

! an exact count, with both numbers printed when it is wrong

  implicit none
  character(len=*), intent(in) :: name
  integer(kind=8), intent(in) :: got,want
  integer, intent(inout) :: nf

  call gf_report_true(name//' = '//trim(itoa8(want)),got == want,nf)
  if (got /= want) write(*,'(a,a,a,i0)') '     mismatch: ',name,' = ',got

  end subroutine expect

!
!-------------------------------------------------------------------------------------------------
!

  function itoa8(i) result(s)

  implicit none
  integer(kind=8), intent(in) :: i
  character(len=24) :: s

  write(s,'(i0)') i

  end function itoa8

!
!-------------------------------------------------------------------------------------------------
!

  subroutine die(msg)

  implicit none
  character(len=*), intent(in) :: msg

  write(*,'(a)') '  '//msg
  write(*,'(a)') '  '//trim(gf_errmsg)
  stop 1

  end subroutine die

  end program test_gf_cache
