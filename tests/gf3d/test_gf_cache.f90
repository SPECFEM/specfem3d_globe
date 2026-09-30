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
!---- lifetime.
!----
!---- Runs on the fixture database, or on $GF3D_TEST_GFDB.
!----

  program test_gf_cache

  use gf3d, only: t_gfdb,t_gf_source,t_gf_location, &
                  gf_open,gf_close,gf_cache_stats,gf_read_source, &
                  gf_locate_source,gf_locate_release,get_seismograms,gf_default_t0, &
                  gf_errmsg,GF_OK,GF_ERR_ARG

  use gf_manufactured, only: gf_report_true

  implicit none

  character(len=512) :: dbpath,cmtpath
  type(t_gfdb) :: db,never_opened
  type(t_gf_source) :: src
  type(t_gf_location) :: loc
  double precision, dimension(:,:,:), allocatable :: synt
  double precision :: t0
  integer(kind=8) :: hits,misses,evictions,files_read
  integer(kind=8) :: misses0,files0,nread_locate,nread_extract
  integer :: n_cached,ierr,nfail

  nfail = 0

  if (command_argument_count() < 2) then
    write(*,'(a)') 'usage: test_gf_cache <GFDB> <CMTSOLUTION>'
    stop 1
  endif
  call get_command_argument(1,dbpath)
  call get_command_argument(2,cmtpath)

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
  ! 3. closing frees it, and a reopened handle starts again
  !--------------------------------------------------------------------

  write(*,'(a)') ''
  write(*,'(a)') '3. close'

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
