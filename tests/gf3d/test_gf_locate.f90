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
!---- test_gf_locate -- source location against the solver's own output
!----
!---- This is the one step in the whole project with no second
!---- implementation to lean on. gf_cross_validate.py never converts
!---- lat/lon/depth to Cartesian at all: parse_solver_output() reads x, y, z
!---- and nu straight out of OUTPUT_FILES/output_solver.txt. So the only
!---- oracle is the solver's own text output, and this test is where the
!---- geographic chain -- geocentric colatitude, topography, ellipticity on
!---- the surface radius, depth -- is pinned.
!----
!---- The reference values are extracted by the runner and passed as
!---- arguments; grep and awk are the right tools for that and keep the
!---- parsing visible.
!----
!---- On the tolerance. locate_sources.f90:444 prints its position through
!---- sngl(), so the reference carries float32 precision: about 3e-8 at the
!---- magnitude of these coordinates. Asking for better than the oracle can
!---- express would be asking for noise, so the assertion is 1e-7 and the
!---- measured value is printed. On the shipped global example the agreement
!---- is 2.1e-8, i.e. below the last digit the solver prints.
!----
!---- usage:
!----   test_gf_locate <GFDB> <lat> <lon> <depth_km> <x_ref> <y_ref> <z_ref> [<GFDB2>]
!----
!---- The optional second database exercises the kd-tree ownership guard;
!---- see the header of src/gf3d/gf_locate.F90.
!----

  program test_gf_locate

  use gf_par, only: t_gfdb,t_gf_location,gf_errmsg,gf_error_string, &
                    GF_OK,GF_ERR_ARG,GF_ERR_NO_ELEMENT,GF_XI_TOL,GF_ANCHOR_TOL
  use gf_database, only: gf_open,gf_close
  use gf_locate, only: gf_locate_source,gf_locate_release,gf_locate_tree_owner

  use gf_manufactured, only: gf_report,gf_report_true,gf_quiet_nan

  use constants, only: MAX_STRING_LEN,NDIM

  implicit none

  ! agreement with a reference the solver printed through sngl()
  double precision, parameter :: TOL_POSITION = 1.d-7
  ! output_solver.txt prints xi,eta,gamma through sngl()
  double precision, parameter :: TOL_LOCAL_COORD = 1.d-7

  ! |mapped - target| after the Newton iteration, in km. The solver reports
  ! 7.1e-13 km for the shipped global example; 1e-9 km is a micrometre and
  ! still four orders of slack.
  double precision, parameter :: TOL_CLOSURE_KM = 1.d-9

  ! local parameters
  character(len=MAX_STRING_LEN) :: dbpath
  type(t_gfdb) :: db,db2
  type(t_gf_location) :: loc,loc2,loc_again
  integer :: nfail,ierr,nargs,i,j
  double precision :: lat,lon,depth_km
  double precision, dimension(NDIM) :: xyz_ref
  double precision :: err,worst,dot
  integer :: owner_first,owner_second
  logical :: have_reference,have_local_ref
  double precision, dimension(NDIM) :: xi_ref

  nfail = 0

  nargs = command_argument_count()
  if (nargs < 4) then
    write(*,'(a)') 'usage: test_gf_locate <GFDB> <lat> <lon> <depth_km> [<x_ref> <y_ref> <z_ref>]'
    stop 1
  endif

  call get_command_argument(1,dbpath)
  call read_double_arg(2,lat)
  call read_double_arg(3,lon)
  call read_double_arg(4,depth_km)

  ! The solver reference is optional. It describes one specific database --
  ! the forward run that produced it -- so it is supplied only when the
  ! runner was given one to go with the database it is testing. Without it
  ! everything here still runs except the two comparisons against the
  ! solver's own Cartesian position, which are reported instead.
  have_reference = (nargs >= 7)
  xyz_ref(:) = 0.d0
  if (have_reference) then
    call read_double_arg(5,xyz_ref(1))
    call read_double_arg(6,xyz_ref(2))
    call read_double_arg(7,xyz_ref(3))
  endif

  ! The solver's own xi,eta,gamma, when the reference carries them. This is
  ! the only assertion in the suite that says we chose the element the
  ! forward run chose: the Cartesian position above is the same point
  ! whichever element it is expressed in, so it cannot tell them apart.
  have_local_ref = (nargs >= 10)
  xi_ref(:) = 0.d0
  if (have_local_ref) then
    call read_double_arg(8,xi_ref(1))
    call read_double_arg(9,xi_ref(2))
    call read_double_arg(10,xi_ref(3))
  endif

  write(*,'(a)') ''
  write(*,'(a)') 'test_gf_locate'
  write(*,'(a)') ''
  write(*,'(a,a)')       '  database   = ',trim(dbpath)
  write(*,'(a,3es22.14)') '  request    = ',lat,lon,depth_km
  if (have_reference) then
    write(*,'(a,3es22.14)') '  reference  = ',xyz_ref(1),xyz_ref(2),xyz_ref(3)
  else
    write(*,'(a)')          '  reference  = (none supplied; solver comparison not asserted)'
  endif
  write(*,'(a)') ''

  call gf_open(dbpath,db,ierr,check_completion=.false.)
  if (ierr /= GF_OK) then
    write(*,'(a)') '  could not open the database'
    write(*,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
    stop 1
  endif

  !--------------------------------------------------------------------
  ! 1. the locate itself
  !--------------------------------------------------------------------

  call gf_locate_source(db,lat,lon,depth_km,loc,ierr)
  call gf_report_true('locate succeeded                  ',ierr == GF_OK,nfail)
  if (ierr /= GF_OK) then
    write(*,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
    call gf_locate_release()
    call gf_close(db)
    stop 1
  endif

  write(*,'(a,a,a,i0,a,i0,a)') '  element    = ',loc%morton_hex,'  (',loc%ielem,' of ',db%nelem,')'
  write(*,'(a,3es22.14)')      '  xi,eta,gam = ',loc%xi,loc%eta,loc%gamma
  write(*,'(a,3es22.14)')      '  x,y,z      = ',loc%xyz(1),loc%xyz(2),loc%xyz(3)
  write(*,'(a)') ''

  call gf_report_true('element index is in range         ', &
                      loc%ielem >= 1 .and. loc%ielem <= db%nelem,nfail)
  call gf_report_true('morton_hex is set                 ',len_trim(loc%morton_hex) > 0,nfail)

  !--------------------------------------------------------------------
  ! 2. physical position against the solver
  !
  ! Both the target -- our own geographic chain, which never touches the
  ! stored element coordinates -- and the mapped point, which comes out of
  ! the element map at the located xi/eta/gamma. They test different things:
  ! the first the geographic conversion, the second that the Newton iteration
  ! actually converged onto it.
  !--------------------------------------------------------------------

  if (have_reference) then
    worst = 0.d0
    do i = 1,NDIM
      worst = max(worst,abs(loc%xyz_target(i) - xyz_ref(i)))
    enddo
    call gf_report('geographic chain vs output_solver ',worst,TOL_POSITION,nfail)

    worst = 0.d0
    do i = 1,NDIM
      worst = max(worst,abs(loc%xyz(i) - xyz_ref(i)))
    enddo
    call gf_report('mapped position vs output_solver  ',worst,TOL_POSITION,nfail)
  else
    write(*,'(a)') '     no solver reference for this database: the geographic'
    write(*,'(a)') '     chain is not compared (it has no other oracle)'
  endif

  ! The element the solver used, and where in it. A source may legitimately
  ! sit outside its element -- the shipped global example is at
  ! gamma = 1.053 because the forward run put it there -- so the test is
  ! agreement with the solver, never |xi| <= 1.
  if (have_local_ref) then
    worst = max(abs(loc%xi - xi_ref(1)), &
                abs(loc%eta - xi_ref(2)), &
                abs(loc%gamma - xi_ref(3)))
    call gf_report('xi,eta,gamma vs output_solver     ',worst,TOL_LOCAL_COORD,nfail)
  else
    write(*,'(a)') '     no solver xi,eta,gamma for this database: the element'
    write(*,'(a)') '     choice is not compared'
  endif

  call gf_report('Newton closure |mapped-target| km ',loc%distance_km,TOL_CLOSURE_KM,nfail)

  !--------------------------------------------------------------------
  ! 3. containment and geometry
  !--------------------------------------------------------------------

  call gf_report_true('containment rule satisfied        ', &
                      max(abs(loc%xi),abs(loc%eta),abs(loc%gamma)) <= GF_XI_TOL,nfail)

  call gf_report('accepted element anchor residual  ',loc%anchor_err,GF_ANCHOR_TOL,nfail)

  call gf_report_true('jacobian is positive              ',loc%jacobian > 0.d0,nfail)

  !--------------------------------------------------------------------
  ! 4. nu is a rotation
  !
  ! The rows are the local North, East and vertical directions in Cartesian,
  ! so nu nu^T must be the identity with determinant +1. A sign or ordering
  ! slip in gf_source_nu shows up here rather than as a mirrored seismogram
  ! component three stages later.
  !--------------------------------------------------------------------

  worst = 0.d0
  do i = 1,NDIM
    do j = 1,NDIM
      dot = loc%nu(i,1)*loc%nu(j,1) + loc%nu(i,2)*loc%nu(j,2) + loc%nu(i,3)*loc%nu(j,3)
      if (i == j) then
        worst = max(worst,abs(dot - 1.d0))
      else
        worst = max(worst,abs(dot))
      endif
    enddo
  enddo
  call gf_report('nu is orthonormal                 ',worst,1.d-14,nfail)

  ! right-handed: det(nu) = +1 for (N,E,Z-up)? The Harvard convention in
  ! locate_sources.f90:300-311 makes (N,E,Z) left-handed, so the determinant
  ! is -1. Asserting the value rather than the magnitude is the point: it
  ! pins the handedness that the reciprocal force components inherit.
  dot = loc%nu(1,1)*(loc%nu(2,2)*loc%nu(3,3) - loc%nu(2,3)*loc%nu(3,2)) &
      - loc%nu(1,2)*(loc%nu(2,1)*loc%nu(3,3) - loc%nu(2,3)*loc%nu(3,1)) &
      + loc%nu(1,3)*(loc%nu(2,1)*loc%nu(3,2) - loc%nu(2,2)*loc%nu(3,1))
  call gf_report('det(nu) = -1 (N,E,Z-up handedness)',abs(dot + 1.d0),1.d-14,nfail)

  !--------------------------------------------------------------------
  ! 5. a position outside the database is an error, not a guess
  !
  ! find_containing_element (gf_cross_validate.py:189) warns and returns
  ! xi=eta=gamma=0 at the nearest centroid, which yields a seismogram that
  ! looks entirely plausible and is wrong. This asserts the library refuses.
  !--------------------------------------------------------------------

  call gf_locate_source(db,-lat,lon+180.d0,depth_km,loc2,ierr)
  call gf_report_true('antipode is rejected, not guessed ',ierr == GF_ERR_NO_ELEMENT,nfail)
  if (ierr == GF_ERR_NO_ELEMENT) then
    write(*,'(a,a)') '     reported: ',trim(gf_errmsg)
  endif

  ! A NaN passes reduce()'s range test -- both comparisons are false -- and
  ! reaches the kd-tree, whose "found no point" branch is a stop. The screen
  ! is in gf_locate_source, so the Fortran route is covered by it too and
  ! the process is still here to assert that.
  call gf_locate_source(db,gf_quiet_nan(),lon,depth_km,loc2,ierr)
  call gf_report_true('a NaN latitude is refused         ',ierr == GF_ERR_ARG,nfail)
  call gf_locate_source(db,lat,lon,gf_quiet_nan(),loc2,ierr)
  call gf_report_true('a NaN depth is refused            ',ierr == GF_ERR_ARG,nfail)
  call gf_report_true('and the process is still running  ',.true.,nfail)

  !--------------------------------------------------------------------
  ! 6. the kd-tree ownership guard
  !
  ! src/shared/search_kdtree.f90 keeps one tree per process in module
  ! variables, so which open it describes is tracked by t_gfdb%open_id.
  ! Opening the same directory a second time is a different open, and is
  ! therefore the same exercise as a second database would be -- without
  ! needing one. Both shipped examples are gitignored, so a check that needed
  ! two of them ran only on a machine that had built both.
  !
  ! The owner is what is asserted. Two handles on one directory hold
  ! identical centroids, so a tree that was never rebuilt would still return
  ! the right element: only gf_locate_tree_owner() distinguishes the two.
  !--------------------------------------------------------------------

  write(*,'(a)') ''

  call gf_open(dbpath,db2,ierr,check_completion=.false.)
  call gf_report_true('the same database opens twice     ',ierr == GF_OK,nfail)

  if (ierr == GF_OK) then
    call gf_report_true('the second open has its own id    ', &
                        db2%open_id /= db%open_id .and. db2%open_id /= 0,nfail)

    call gf_locate_source(db2,lat,lon,depth_km,loc2,ierr)
    call gf_report_true('locate through the second handle  ',ierr == GF_OK,nfail)
    call gf_report_true('the tree belongs to the second    ', &
                        gf_locate_tree_owner() == db2%open_id,nfail)

    ! back to the first
    call gf_locate_source(db,lat,lon,depth_km,loc_again,ierr)
    call gf_report_true('locate again through the first    ',ierr == GF_OK,nfail)
    call gf_report_true('and the tree is the first''s again ', &
                        gf_locate_tree_owner() == db%open_id,nfail)

    if (ierr == GF_OK) then
      call gf_report_true('same element after the switch     ', &
                          loc_again%ielem == loc%ielem .and. &
                          loc_again%morton_hex == loc%morton_hex,nfail)
      err = max(abs(loc_again%xi - loc%xi), &
                abs(loc_again%eta - loc%eta), &
                abs(loc_again%gamma - loc%gamma))
      call gf_report('identical xi,eta,gamma after switch',err,0.d0,nfail)
    endif

    call gf_close(db2)
  endif

  !--------------------------------------------------------------------
  ! 7. the same path, opened twice, is two different opens
  !
  ! The tree used to be owned by the database's path, which made a handle
  ! closed and reopened indistinguishable from the one that built the tree:
  ! gf_tree_ensure would keep a tree built for arrays that gf_close had
  ! already deallocated. Ownership is t_gfdb%open_id, so the reopen is a
  ! different owner and the tree is rebuilt -- onto the same answer.
  !--------------------------------------------------------------------

  owner_first = gf_locate_tree_owner()
  call gf_report_true('the tree belongs to this open     ', &
                      owner_first == db%open_id .and. db%open_id /= 0,nfail)

  call gf_close(db)
  call gf_open(dbpath,db,ierr,check_completion=.false.)
  call gf_report_true('the same path reopens             ',ierr == GF_OK,nfail)
  call gf_report_true('with a different open id          ',db%open_id /= owner_first,nfail)

  call gf_locate_source(db,lat,lon,depth_km,loc_again,ierr)
  call gf_report_true('and locates again                 ',ierr == GF_OK,nfail)

  owner_second = gf_locate_tree_owner()
  call gf_report_true('the tree was rebuilt for it       ', &
                      owner_second == db%open_id .and. owner_second /= owner_first,nfail)

  if (ierr == GF_OK) then
    call gf_report_true('onto the same element             ', &
                        loc_again%ielem == loc%ielem .and. &
                        loc_again%morton_hex == loc%morton_hex,nfail)
    err = max(abs(loc_again%xi - loc%xi), &
              abs(loc_again%eta - loc%eta), &
              abs(loc_again%gamma - loc%gamma))
    call gf_report('identical xi,eta,gamma after reopen',err,0.d0,nfail)
  endif

  !--------------------------------------------------------------------

  call gf_locate_release()
  call gf_close(db)

  write(*,'(a)') ''
  if (nfail /= 0) then
    write(*,'(a,i0,a)') 'test_gf_locate: ',nfail,' assertion(s) FAILED'
    stop 1
  endif

  write(*,'(a)') 'test_gf_locate: all assertions passed'

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_double_arg(iarg,val)

! reads a command-line argument as a double

  implicit none

  integer, intent(in) :: iarg
  double precision, intent(out) :: val

  ! local parameters
  character(len=MAX_STRING_LEN) :: str
  integer :: ios

  call get_command_argument(iarg,str)
  read(str,*,iostat=ios) val
  if (ios /= 0) then
    write(*,'(a,i0,a,a,a)') 'could not read argument ',iarg,' ("',trim(str),'") as a number'
    stop 1
  endif

  end subroutine read_double_arg

  end program test_gf_locate
