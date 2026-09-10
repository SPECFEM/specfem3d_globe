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
!---- test_gf_source -- the values route against the file route
!----
!---- Stage 9 gives the C/Python side a way to build a source from numbers
!---- rather than from a file, because reading a file means get_cmt() and
!---- get_force(), which `stop` on malformed input -- and a stop inside
!---- lib/libgf3d.so kills the calling interpreter with no traceback.
!----
!---- That leaves two code paths where there was one, and they must agree.
!---- This test writes a source file, reads it through the solver's own
!---- reader, builds the same source from the same numbers through
!---- gf_source_set_*, and compares. The oracle is the reader; what is being
!---- checked is the transcription of everything the reader does *after*
!---- parsing -- the half-duration clamps, the time-shift convention and the
!---- non-dimensionalisation -- since that is what was copied out.
!----
!---- Tolerance: 1e-12 relative, not equality. The two routes are separate
!---- compilation units evaluating the same expressions, and testing.md
!---- forbids asserting equality of a value a test recomputes: the Intel CI
!---- has twice failed such checks that no local target reproduces. An
!---- error here of the size that matters -- a missing clamp, a forgotten
!---- 1e7, a dropped division -- is a factor, not an ulp.
!----
!---- The second half checks the paths that used to be a `stop`: every one
!---- returns an error code, sets a message, and leaves the program running.
!---- That is the property the whole facade rests on, so it is asserted
!---- rather than assumed.
!----
!---- Tier 1: no database, no HDF5, no MPI. Runs on every commit.
!----

  program test_gf_source

  use gf_par, only: t_gf_source,gf_errmsg,GF_OK,GF_ERR_ARG,GF_SRC_CMT,GF_SRC_FORCE

  use gf_shared_params, only: gf_init_shared_params

  use gf_source, only: gf_read_cmt_source,gf_read_force_source, &
                       gf_source_set_cmt,gf_source_set_force

  use gf_manufactured, only: gf_report,gf_report_true

  implicit none

  ! the shipped regional example's own CMTSOLUTION values, so that the
  ! numbers under test are the ones the project is validated with
  double precision, parameter :: LAT = -5.8120d0
  double precision, parameter :: LON = -75.2700d0
  double precision, parameter :: DEP = 122.6000d0
  double precision, parameter :: HDUR = 60.0000d0
  double precision, parameter :: TSHIFT = 29.0000d0

  double precision, parameter :: MRR = -7.590000d27
  double precision, parameter :: MTT =  7.750000d27
  double precision, parameter :: MPP = -1.600000d26
  double precision, parameter :: MRT = -2.503000d28
  double precision, parameter :: MRP =  4.200000d26
  double precision, parameter :: MTP = -2.480000d27

  ! the regional database's solver time step
  double precision, parameter :: DT = 0.1d0

  double precision, parameter :: TOL = 1.d-12

  character(len=*), parameter :: CMTFILE = './OUTPUT_FILES/test_gf_source.CMTSOLUTION'
  character(len=*), parameter :: FORCEFILE = './OUTPUT_FILES/test_gf_source.FORCESOLUTION'

  integer :: nfail

  nfail = 0

  write(*,*)
  write(*,*) '******************************'
  write(*,*) 'test_gf_source'
  write(*,*) '******************************'
  write(*,*)

  ! get_cmt reads NUMBER_OF_SIMULTANEOUS_RUNS and NOISE_TOMOGRAPHY, neither
  ! of which has an initialiser in shared_par.f90, and both routes read
  ! R_PLANET and RHOAV. A default t_gfdb leaves the Earth defaults standing,
  ! which is what a tier-1 test wants: no database exists here.
  call init_globals(nfail)

  call test_cmt_agrees(nfail)
  call test_cmt_clamp(nfail)
  call test_force_agrees(nfail)
  call test_refusals(nfail)

  write(*,*)
  if (nfail == 0) then
    write(*,*) 'test_gf_source: all assertions passed'
  else
    write(*,*) 'test_gf_source: ',nfail,' assertion(s) FAILED'
    stop 1
  endif
  write(*,*)

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine init_globals(nfail)

  use gf_par, only: t_gfdb

  implicit none

  integer, intent(inout) :: nfail

  ! local parameters
  type(t_gfdb) :: db
  integer :: ier

  call gf_init_shared_params(db,ier)

  call gf_report_true('shared parameters installed from a blank handle',ier == GF_OK,nfail)

  end subroutine init_globals

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_cmt(filename,hdur_in,tshift_in)

! a 13-line CMTSOLUTION, in the format get_cmt expects

  implicit none

  character(len=*), intent(in) :: filename
  double precision, intent(in) :: hdur_in,tshift_in

  ! local parameters
  integer, parameter :: IOUT = 91

  open(unit=IOUT,file=trim(filename),status='unknown',action='write')
  write(IOUT,'(a)') 'PDE 1994  6  9  0 33 16.40 -13.8300  -67.5600 637.0 6.9 6.8 NORTHERN BOLIVIA'
  write(IOUT,'(a)') 'event name:     test0001'
  write(IOUT,'(a,f12.4)') 'time shift:  ',tshift_in
  write(IOUT,'(a,f12.4)') 'half duration:',hdur_in
  write(IOUT,'(a,f12.4)') 'latitude:    ',LAT
  write(IOUT,'(a,f12.4)') 'longitude:   ',LON
  write(IOUT,'(a,f12.4)') 'depth:       ',DEP
  write(IOUT,'(a,es16.6)') 'Mrr:       ',MRR
  write(IOUT,'(a,es16.6)') 'Mtt:       ',MTT
  write(IOUT,'(a,es16.6)') 'Mpp:       ',MPP
  write(IOUT,'(a,es16.6)') 'Mrt:       ',MRT
  write(IOUT,'(a,es16.6)') 'Mrp:       ',MRP
  write(IOUT,'(a,es16.6)') 'Mtp:       ',MTP
  close(IOUT)

  end subroutine write_cmt

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_force(filename,f0_in,stf_in,factor_in,dE,dN,dZ)

! an 11-line FORCESOLUTION, in the globe's format (f0, not a half duration)

  implicit none

  character(len=*), intent(in) :: filename
  double precision, intent(in) :: f0_in,factor_in,dE,dN,dZ
  integer, intent(in) :: stf_in

  ! local parameters
  integer, parameter :: IOUT = 92

  open(unit=IOUT,file=trim(filename),status='unknown',action='write')
  write(IOUT,'(a)') 'FORCE  001'
  write(IOUT,'(a,f12.4)') 'time shift:  ',0.d0
  write(IOUT,'(a,f12.4)') 'f0:          ',f0_in
  write(IOUT,'(a,f12.4)') 'latitude:    ',LAT
  write(IOUT,'(a,f12.4)') 'longitude:   ',LON
  write(IOUT,'(a,f12.4)') 'depth:       ',DEP
  write(IOUT,'(a,i6)')    'source time function:',stf_in
  write(IOUT,'(a,es16.6)') 'factor force source:',factor_in
  write(IOUT,'(a,f12.4)') 'comp dir vect source E:',dE
  write(IOUT,'(a,f12.4)') 'comp dir vect source N:',dN
  write(IOUT,'(a,f12.4)') 'comp dir vect source Z:',dZ
  close(IOUT)

  end subroutine write_force

!
!-------------------------------------------------------------------------------------------------
!

  double precision function relerr(a,b)

! relative difference, falling back to absolute when both are ~zero

  implicit none

  double precision, intent(in) :: a,b

  ! local parameters
  double precision :: scale

  scale = max(abs(a),abs(b))
  if (scale < 1.d-30) then
    relerr = abs(a-b)
  else
    relerr = abs(a-b)/scale
  endif

  end function relerr

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_cmt_agrees(nfail)

! the two routes on the shipped example's own moment tensor

  implicit none

  integer, intent(inout) :: nfail

  ! local parameters
  type(t_gf_source) :: sfile,sval
  double precision, dimension(6) :: m
  integer :: ier,iv
  character(len=16) :: label

  write(*,*) '1. CMT: values route against get_cmt'

  m = (/ MRR,MTT,MPP,MRT,MRP,MTP /)

  call write_cmt(CMTFILE,HDUR,TSHIFT)

  call gf_read_cmt_source(CMTFILE,DT,sfile,ier)
  call gf_report_true('   get_cmt read the file',ier == GF_OK,nfail)
  if (ier /= GF_OK) return

  call gf_source_set_cmt(sval,LAT,LON,DEP,HDUR,TSHIFT,m,DT,ier)
  call gf_report_true('   gf_source_set_cmt succeeded',ier == GF_OK,nfail)
  if (ier /= GF_OK) return

  call gf_report_true('   source_type is CMT',sval%source_type == GF_SRC_CMT,nfail)

  call gf_report('   latitude   ',relerr(sval%latitude,sfile%latitude),TOL,nfail)
  call gf_report('   longitude  ',relerr(sval%longitude,sfile%longitude),TOL,nfail)
  call gf_report('   depth      ',relerr(sval%depth,sfile%depth),TOL,nfail)
  call gf_report('   hdur       ',relerr(sval%hdur,sfile%hdur),TOL,nfail)
  call gf_report('   tshift_src ',abs(sval%tshift_src - sfile%tshift_src),TOL,nfail)
  call gf_report('   time shift ',relerr(sval%min_tshift_src_original, &
                                         sfile%min_tshift_src_original),TOL,nfail)
  call gf_report('   scale_moment',relerr(sval%scale_moment,sfile%scale_moment),TOL,nfail)

  ! the six non-dimensional components, named one by one: the CI log is the
  ! only trace of a run there, and "the moment tensor differs" would not say
  ! which convention slipped
  do iv = 1,6
    write(label,'(a,i1,a)') '   moment(',iv,')  '
    call gf_report(label,relerr(sval%moment_tensor(iv),sfile%moment_tensor(iv)),TOL,nfail)
  enddo

  ! and the round trip back to dyne-cm, which is what the partials are per
  do iv = 1,6
    write(label,'(a,i1,a)') '   dyne-cm(',iv,') '
    call gf_report(label,relerr(sval%moment_tensor(iv)*sval%scale_moment,m(iv)),TOL,nfail)
  enddo

  end subroutine test_cmt_agrees

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_cmt_clamp(nfail)

! the half-duration clamp: get_cmt.f90:397 replaces a null half duration by
! 5*DT, and a values route that forgot it would give a Heaviside where the
! solver gives a short error function

  implicit none

  integer, intent(inout) :: nfail

  ! local parameters
  type(t_gf_source) :: sfile,sval
  double precision, dimension(6) :: m
  integer :: ier

  write(*,*) '2. CMT: the half-duration clamp'

  m = (/ MRR,MTT,MPP,MRT,MRP,MTP /)

  call write_cmt(CMTFILE,0.d0,TSHIFT)

  call gf_read_cmt_source(CMTFILE,DT,sfile,ier)
  if (ier /= GF_OK) then
    call gf_report_true('   get_cmt read the zero-hdur file',.false.,nfail)
    return
  endif

  call gf_source_set_cmt(sval,LAT,LON,DEP,0.d0,TSHIFT,m,DT,ier)
  if (ier /= GF_OK) then
    call gf_report_true('   gf_source_set_cmt on zero hdur',.false.,nfail)
    return
  endif

  call gf_report('   hdur clamped equally',relerr(sval%hdur,sfile%hdur),TOL,nfail)
  call gf_report('   hdur is 5*dt        ',relerr(sval%hdur,5.d0*DT),TOL,nfail)

  end subroutine test_cmt_clamp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_force_agrees(nfail)

! every source time function type the reader accepts

  implicit none

  integer, intent(inout) :: nfail

  ! local parameters
  type(t_gf_source) :: sfile,sval
  integer :: ier,k,istf
  integer, dimension(4), parameter :: STF_KINDS = (/ 0,1,2,4 /)
  double precision, parameter :: FACTOR = 1.d15
  double precision, parameter :: DE = 1.d0, DN = 0.5d0, DZ = -2.d0
  double precision :: f0
  character(len=20) :: label

  write(*,*) '3. force: values route against get_force, per source time function'

  do k = 1,4
    istf = STF_KINDS(k)

    ! a Ricker reads f0 as a frequency, the others as a width; either way
    ! the number is small enough that the 5*DT clamp is exercised on the
    ! Gaussian branches
    f0 = 0.05d0

    call write_force(FORCEFILE,f0,istf,FACTOR,DE,DN,DZ)

    call gf_read_force_source(FORCEFILE,DT,sfile,ier)
    write(label,'(a,i1,a)') '   stf ',istf,' read      '
    call gf_report_true(label,ier == GF_OK,nfail)
    if (ier /= GF_OK) cycle

    call gf_source_set_force(sval,LAT,LON,DEP,f0,0.d0,istf,FACTOR,DE,DN,DZ,DT,ier)
    write(label,'(a,i1,a)') '   stf ',istf,' built     '
    call gf_report_true(label,ier == GF_OK,nfail)
    if (ier /= GF_OK) cycle

    write(label,'(a,i1,a)') '   stf ',istf,' hdur      '
    call gf_report(label,relerr(sval%hdur,sfile%hdur),TOL,nfail)

    write(label,'(a,i1,a)') '   stf ',istf,' factor    '
    call gf_report(label,relerr(sval%factor_force_source,sfile%factor_force_source),TOL,nfail)

    write(label,'(a,i1,a)') '   stf ',istf,' dir E     '
    call gf_report(label,relerr(sval%comp_dir_vect_source_E, &
                                sfile%comp_dir_vect_source_E),TOL,nfail)
    write(label,'(a,i1,a)') '   stf ',istf,' dir N     '
    call gf_report(label,relerr(sval%comp_dir_vect_source_N, &
                                sfile%comp_dir_vect_source_N),TOL,nfail)
    write(label,'(a,i1,a)') '   stf ',istf,' dir Z     '
    call gf_report(label,relerr(sval%comp_dir_vect_source_Z_UP, &
                                sfile%comp_dir_vect_source_Z_UP),TOL,nfail)

    write(label,'(a,i1,a)') '   stf ',istf,' type      '
    call gf_report_true(label,sval%source_type == GF_SRC_FORCE .and. &
                              sval%force_stf == sfile%force_stf,nfail)
  enddo

  end subroutine test_force_agrees

!
!-------------------------------------------------------------------------------------------------
!

  subroutine test_refusals(nfail)

! the paths that are a `stop` in the readers
!
! Each of these would end the process inside get_force(); here each must
! return GF_ERR_ARG, leave a message behind, and let the program continue --
! which it demonstrably does, because the ones after it still run.

  implicit none

  integer, intent(inout) :: nfail

  ! local parameters
  type(t_gf_source) :: s
  double precision, dimension(6) :: m
  double precision :: nan,inf
  integer :: ier

  write(*,*) '4. refusals: what used to be a stop'

  m = (/ MRR,MTT,MPP,MRT,MRP,MTP /)

  ! built by bit pattern rather than by 0/0, which would trap under
  ! -ffpe-trap=invalid before it could produce anything
  nan = transfer(int(z'7FF8000000000000',kind=8),nan)
  inf = transfer(int(z'7FF0000000000000',kind=8),inf)

  ! get_force.f90:248 -- an unsupported source time function type
  call gf_source_set_force(s,LAT,LON,DEP,0.05d0,0.d0,5,1.d15,1.d0,0.d0,0.d0,DT,ier)
  call gf_report_true('   force_stf = 5 refused         ',ier == GF_ERR_ARG,nfail)
  call gf_report_true('   ... with a message            ',len_trim(gf_errmsg) > 0,nfail)

  ! get_force.f90:240 -- a monochromatic force with no period
  call gf_source_set_force(s,LAT,LON,DEP,0.d0,0.d0,3,1.d15,1.d0,0.d0,0.d0,DT,ier)
  call gf_report_true('   monochromatic f0 = 0 refused  ',ier == GF_ERR_ARG,nfail)

  ! get_force.f90:272 -- a direction vector of zero length
  call gf_source_set_force(s,LAT,LON,DEP,0.05d0,0.d0,0,1.d15,0.d0,0.d0,0.d0,DT,ier)
  call gf_report_true('   zero direction vector refused ',ier == GF_ERR_ARG,nfail)

  ! not a reader case: the library's own precondition
  call gf_source_set_cmt(s,LAT,LON,DEP,HDUR,TSHIFT,m,0.d0,ier)
  call gf_report_true('   dt = 0 refused (CMT)          ',ier == GF_ERR_ARG,nfail)

  ! the ones that would reach the kd-tree's own stop through gf_locate
  call gf_source_set_cmt(s,nan,LON,DEP,HDUR,TSHIFT,m,DT,ier)
  call gf_report_true('   NaN latitude refused          ',ier == GF_ERR_ARG,nfail)

  call gf_source_set_cmt(s,LAT,LON,inf,HDUR,TSHIFT,m,DT,ier)
  call gf_report_true('   infinite depth refused        ',ier == GF_ERR_ARG,nfail)

  m(3) = nan
  call gf_source_set_cmt(s,LAT,LON,DEP,HDUR,TSHIFT,m,DT,ier)
  call gf_report_true('   NaN moment component refused  ',ier == GF_ERR_ARG,nfail)
  m(3) = MPP

  call gf_source_set_force(s,LAT,nan,DEP,0.05d0,0.d0,0,1.d15,1.d0,0.d0,0.d0,DT,ier)
  call gf_report_true('   NaN longitude refused (force) ',ier == GF_ERR_ARG,nfail)

  ! and the program is still here to say so
  call gf_source_set_cmt(s,LAT,LON,DEP,HDUR,TSHIFT,m,DT,ier)
  call gf_report_true('   a good source still builds    ',ier == GF_OK,nfail)

  end subroutine test_refusals

  end program test_gf_source
