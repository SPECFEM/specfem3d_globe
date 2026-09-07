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
!---- test_gf_cmt -- the moment-tensor amplitude, against the solver
!----
!---- This is the amplitude pin. A pure scale error is invisible in waveform
!---- shape and fatal for a moment inversion, so it must not be able to hide
!---- behind a shape-only comparison -- and the whole chain from CMTSOLUTION
!---- dyne-cm through scaleM to a non-dimensional moment tensor is exactly
!---- the kind of thing that goes wrong by a constant.
!----
!---- **The oracle is specfem itself**, in two independent senses:
!----
!---- 1. OUTPUT_FILES/output_solver.txt records the scalar moment and moment
!----    magnitude the solver computed for this very CMTSOLUTION, at full
!----    double precision. Our source, read through get_cmt and
!----    non-dimensionalised by scaleM, must reproduce them.
!---- 2. get_cmt_scalar_moment() is the solver's own routine and multiplies
!----    scaleM back in. Feeding it our moment tensor closes the loop on the
!----    non-dimensionalisation without this test hard-coding a single
!----    constant of its own.
!----
!---- The Euclidean norm of a symmetric tensor is invariant under rotation
!---- (get_cmt.f90 says so where it notes the routine accepts either the
!---- spherical or the Cartesian components), which gives a third assertion
!---- for free: gf_rotate_moment_tensor must not change M0. That pins the
!---- rotation against a specfem routine rather than against a re-derivation.
!----
!---- usage:
!----   test_gf_cmt <GFDB> <CMTSOLUTION> <M0_ref> <Mw_ref> <hdur_ref> <tshift_ref>
!----

  program test_gf_cmt

  use gf_par, only: t_gfdb,t_gf_source,t_gf_location,gf_errmsg,gf_error_string, &
                    GF_OK,GF_SRC_CMT
  use gf_database, only: gf_open,gf_close
  use gf_locate, only: gf_locate_source,gf_locate_release
  use gf_source, only: gf_read_source
  use gf_moment, only: gf_rotate_moment_tensor

  use gf_manufactured, only: gf_report,gf_report_true

  use constants, only: MAX_STRING_LEN,NDIM

  implicit none

  ! local parameters
  character(len=MAX_STRING_LEN) :: dbpath,cmtfile
  type(t_gfdb) :: db
  type(t_gf_source) :: src
  type(t_gf_location) :: loc
  double precision, dimension(NDIM,NDIM) :: m_cart
  double precision :: m0_ref,mw_ref,hdur_ref,tshift_ref
  double precision :: m0,mw,m0_cart
  integer :: nfail,ierr,nargs

  double precision, external :: get_cmt_scalar_moment
  double precision, external :: get_cmt_moment_magnitude

  nfail = 0

  nargs = command_argument_count()
  if (nargs < 6) then
    write(*,'(a)') 'usage: test_gf_cmt <GFDB> <CMTSOLUTION> <M0> <Mw> <hdur> <tshift>'
    stop 1
  endif

  call get_command_argument(1,dbpath)
  call get_command_argument(2,cmtfile)
  call read_double_arg(3,m0_ref)
  call read_double_arg(4,mw_ref)
  call read_double_arg(5,hdur_ref)
  call read_double_arg(6,tshift_ref)

  write(*,'(a)') ''
  write(*,'(a)') 'test_gf_cmt'
  write(*,'(a)') ''
  write(*,'(a,a)')       '  database  = ',trim(dbpath)
  write(*,'(a,a)')       '  source    = ',trim(cmtfile)
  write(*,'(a,es24.16)') '  solver M0 = ',m0_ref
  write(*,'(a,es24.16)') '  solver Mw = ',mw_ref
  write(*,'(a)') ''

  call gf_open(dbpath,db,ierr,check_completion=.false.)
  if (ierr /= GF_OK) then
    write(*,'(a)') '  could not open the database'
    write(*,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
    stop 1
  endif

  !--------------------------------------------------------------------
  ! 1. reading the source
  !
  ! get_cmt is reused rather than re-implemented precisely so that scaleM
  ! comes from the solver; these assertions check that the reuse is wired
  ! up, including the two semantics that live outside the reader.
  !--------------------------------------------------------------------

  call gf_read_source(cmtfile,db%dt,src,ierr)
  call gf_report_true('CMTSOLUTION read                  ',ierr == GF_OK,nfail)
  if (ierr /= GF_OK) then
    write(*,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
    call gf_close(db)
    stop 1
  endif

  call gf_report_true('detected as a moment tensor       ',src%source_type == GF_SRC_CMT,nfail)

  ! get_cmt returns the raw *triangle* half duration; the /1.628 conversion
  ! to a Gaussian width is the caller's, and is Stage 5's
  call gf_report('half duration vs output_solver    ', &
                 abs(src%hdur - hdur_ref),1.d-12,nfail)

  ! get_cmt zeroes tshift_src when NSOURCES == 1 and returns the original in
  ! min_tshift_src_original, so t = 0 is the centroid time. The solver
  ! prints the zeroed value, which is what this compares against.
  call gf_report('time shift vs output_solver       ', &
                 abs(src%tshift_src - tshift_ref),1.d-12,nfail)

  write(*,'(a,es24.16)') '     min_tshift_src_original = ',src%min_tshift_src_original
  call gf_report_true('origin time shift is non-zero     ', &
                      abs(src%min_tshift_src_original) > 0.d0,nfail)

  !--------------------------------------------------------------------
  ! 2. the amplitude, against the solver's own numbers
  !
  ! get_cmt_scalar_moment multiplies scaleM back in, so this is a closed
  ! loop over the non-dimensionalisation: a wrong scaleM anywhere shows up
  ! immediately, and no constant is hard-coded here.
  !--------------------------------------------------------------------

  m0 = get_cmt_scalar_moment(src%moment_tensor(1),src%moment_tensor(2), &
                             src%moment_tensor(3),src%moment_tensor(4), &
                             src%moment_tensor(5),src%moment_tensor(6))

  write(*,'(a,es24.16)') '     our M0                  = ',m0

  call gf_report('scalar moment M0 vs the solver    ', &
                 abs(m0 - m0_ref)/abs(m0_ref),1.d-12,nfail)

  mw = get_cmt_moment_magnitude(src%moment_tensor(1),src%moment_tensor(2), &
                                src%moment_tensor(3),src%moment_tensor(4), &
                                src%moment_tensor(5),src%moment_tensor(6))

  call gf_report('moment magnitude Mw vs the solver ', &
                 abs(mw - mw_ref)/abs(mw_ref),1.d-12,nfail)

  !--------------------------------------------------------------------
  ! 3. the spherical-to-Cartesian rotation preserves M0
  !
  ! The Euclidean norm of a symmetric tensor is rotation invariant, so
  ! feeding the rotated components to the same solver routine must give the
  ! same M0. A sign slip in any one of the six expanded expressions in
  ! gf_rotate_moment_tensor changes the norm and fails here.
  !
  ! The source has to be located first, because the rotation is about the
  ! source's own geocentric colatitude and longitude -- which is itself the
  ! quantity test_gf_locate pinned against output_solver.txt.
  !--------------------------------------------------------------------

  call gf_locate_source(db,src%latitude,src%longitude,src%depth,loc,ierr)
  call gf_report_true('source located                    ',ierr == GF_OK,nfail)

  if (ierr == GF_OK) then
    call gf_rotate_moment_tensor(loc%theta,loc%phi,src%moment_tensor,m_cart)

    m0_cart = get_cmt_scalar_moment(m_cart(1,1),m_cart(2,2),m_cart(3,3), &
                                    m_cart(1,2),m_cart(1,3),m_cart(2,3))

    call gf_report('rotation preserves M0             ', &
                   abs(m0_cart - m0)/abs(m0),1.d-12,nfail)

    ! and the rotated tensor is still symmetric and still traceless-preserving
    call gf_report('rotated tensor is symmetric       ', &
                   max(abs(m_cart(1,2)-m_cart(2,1)), &
                       abs(m_cart(1,3)-m_cart(3,1)), &
                       abs(m_cart(2,3)-m_cart(3,2))),0.d0,nfail)

    call gf_report('rotation preserves the trace      ', &
                   abs((m_cart(1,1)+m_cart(2,2)+m_cart(3,3)) &
                       - (src%moment_tensor(1)+src%moment_tensor(2)+src%moment_tensor(3))) &
                   / abs(m0/1.d0),1.d-12,nfail)
  endif

  !--------------------------------------------------------------------

  call gf_locate_release()
  call gf_close(db)

  write(*,'(a)') ''
  if (nfail /= 0) then
    write(*,'(a,i0,a)') 'test_gf_cmt: ',nfail,' assertion(s) FAILED'
    stop 1
  endif

  write(*,'(a)') 'test_gf_cmt: all assertions passed'

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

  end program test_gf_cmt
