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
!---- test_gf_anchors -- the 27-anchor premise, against a real database
!----
!---- The whole extraction rests on one geometric claim: that the element
!---- coordinate map is the tri-quadratic one through 27 anchor nodes. That
!---- is true because setup/constants.h has USE_GLL = .false., so
!---- compute_element_properties.f90 applies topography and ellipticity to
!---- the 27 anchors and re-interpolates the GLL points with shape3D --
!---- meaning the stored xyz(3,5,5,5) *is* a tri-quadratic sampled at GLL
!---- points, and the anchors must reproduce all 125 of them.
!----
!---- This test checks that against every element of a real database, which
!---- makes it 125 x nelem independently stored values versus one map. If it
!---- fails, the database came from a USE_GLL = .true. mesh and the anchor
!---- route -- coordinates, Jacobian, strain, everything downstream -- is
!---- invalid.
!----
!---- On the tolerance. The stage plan originally asked for 1e-12. That is
!---- unreachable and the reason is worth recording: the solver holds
!---- xstore_crust_mantle as real(CUSTOM_REAL) and the writer copies it
!---- straight through (green_function_metadata.F90:82,119-121), so with
!---- CUSTOM_REAL = 4 both the anchors and the values they must reproduce
!---- carry a float32 rounding. Measured on the shipped examples the residual
!---- is ~6.2e-8, which is float32 epsilon times the O(0.5) magnitude of a
!---- non-dimensional coordinate -- about 0.4 m. That is the mesh's own
!---- precision, shared with the solver, not an error the database adds; and
!---- since the acceptance criterion for this work is reproducing a forward
!---- run that used exactly those float32 anchors, matching them is correct.
!----
!---- The check keeps its power regardless: a USE_GLL = .true. database
!---- misses by 1e-3 to 1e-2, four orders above the tolerance.
!----

  program test_gf_anchors

  use gf_par, only: t_gfdb,gf_errmsg,gf_error_string,GF_OK,GF_ANCHOR_TOL
  use gf_database, only: gf_open,gf_close
  use gf_locate, only: gf_check_anchors,gf_check_anchors_all,gf_locate_release

  use gf_manufactured, only: gf_report,gf_report_true

  use constants, only: MAX_STRING_LEN

  implicit none

  ! local parameters
  character(len=MAX_STRING_LEN) :: dbpath
  type(t_gfdb) :: db
  integer :: nfail,ierr,nargs,ielem_worst
  double precision :: worst_err,err_one

  nfail = 0

  nargs = command_argument_count()
  if (nargs < 1) then
    write(*,'(a)') 'usage: test_gf_anchors <GFDB>'
    stop 1
  endif
  call get_command_argument(1,dbpath)

  write(*,'(a)') ''
  write(*,'(a)') 'test_gf_anchors'
  write(*,'(a)') ''
  write(*,'(a,a)') '  database = ',trim(dbpath)
  write(*,'(a)') ''

  ! the per-element completion scan is not what is being tested here, and it
  ! is the slow part of gf_open on a large database
  call gf_open(dbpath,db,ierr,check_completion=.false.)
  if (ierr /= GF_OK) then
    write(*,'(a)') '  could not open the database'
    write(*,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
    stop 1
  endif

  write(*,'(a,i0)') '  elements = ',db%nelem
  write(*,'(a)') ''

  call gf_report_true('database has elements             ',db%nelem > 0,nfail)

  !--------------------------------------------------------------------
  ! every element, every GLL point
  !--------------------------------------------------------------------

  call gf_check_anchors_all(db,worst_err,ielem_worst,ierr)
  call gf_report_true('anchor sweep completed            ',ierr == GF_OK,nfail)

  if (ierr /= GF_OK) then
    write(*,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
    call gf_close(db)
    stop 1
  endif

  call gf_report('27 anchors reproduce 125 GLL coords',worst_err,GF_ANCHOR_TOL,nfail)

  write(*,'(a,es22.14)') '     worst residual in metres = ',worst_err*db%R_PLANET
  if (ielem_worst > 0) then
    write(*,'(a,a)')     '     worst element            = ',db%morton_hex(ielem_worst)
  endif

  ! Records the measured value against the float32 floor rather than only
  ! against the tolerance. float32 epsilon is 1.19e-7 and the coordinates are
  ! O(1), so a residual far *below* ~1e-8 would mean the database was written
  ! in double precision (a CUSTOM_REAL = 8 solver), and one far above 1e-7
  ! would mean the map is not tri-quadratic at all. Both are worth seeing in
  ! results.log rather than inferring later.
  call gf_report_true('residual sits at the float32 floor',worst_err < 1.d-6,nfail)
  if (worst_err < 1.d-10) then
    write(*,'(a)') '     note: residual well below the float32 floor;'
    write(*,'(a)') '           this database was probably written by a CUSTOM_REAL = 8 solver'
  endif

  !--------------------------------------------------------------------
  ! the single-element entry point must agree with the sweep
  !
  ! gf_locate_source() calls gf_check_anchors() on the element it accepts,
  ! once per locate, rather than sweeping the database at open time. The two
  ! paths must not drift apart.
  !--------------------------------------------------------------------

  if (ielem_worst > 0) then
    call gf_check_anchors(db,ielem_worst,err_one,ierr)
    call gf_report_true('single-element check completed    ',ierr == GF_OK,nfail)
    call gf_report_true('single-element check == sweep     ',err_one == worst_err,nfail)
  endif

  !--------------------------------------------------------------------

  call gf_locate_release()
  call gf_close(db)

  write(*,'(a)') ''
  if (nfail /= 0) then
    write(*,'(a,i0,a)') 'test_gf_anchors: ',nfail,' assertion(s) FAILED'
    stop 1
  endif

  write(*,'(a)') 'test_gf_anchors: all assertions passed'

  end program test_gf_anchors
