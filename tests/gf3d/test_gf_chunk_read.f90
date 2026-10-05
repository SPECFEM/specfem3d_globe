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
!---- test_gf_chunk_read -- the direct chunk reader against h5dread_f
!----
!---- gf_h5_read_displ_chunks reads a displacement dataset chunk by chunk
!---- through H5Dread_chunk when the file is laid out as the solver writes
!---- it, and through h5dread_f otherwise. Either way it must hand back
!---- exactly the stored floats, in the order the caller asked for. The
!---- oracle is this program's own h5dread_f of the whole dataset, the way
!---- the library read every displacement file before this reader existed.
!---- For every station of the first NELEM_MAX elements of a database (all of
!---- the fixture's):
!----
!----   0. gf_read_element_displ, which every extraction reads through and
!----      which now goes through gf_h5_read_displ_chunks;
!----   1. both orders, at time prefixes 1, 7, 8 and nt -- 7 is the fixture's
!----      chunk length, so these are one row of a chunk, one whole chunk,
!----      one more than a chunk, and everything including the partial last
!----      chunk -- with a sentinel after the buffer, so that copying a
!----      partial chunk at full length is caught;
!----      and the route taken, which the caller states: `raw` for a chunked
!----      database (the solver's own, or the fixture's chunked variant),
!----      `fallback` for the contiguous fixture, `any` when unknown. A silent
!----      fallback must not pass for the chunk path;
!----   2. a file this program writes with one chunk never written, as an
!----      incomplete database has: H5Dread_chunk refuses that chunk, so the
!----      reader must fall back and return the fill value there. Its twin,
!----      written completely, must take the raw route, which is what makes
!----      the first assertion mean anything.
!----
!---- Every comparison is bitwise: nothing here is computed, only read.
!----

  program test_gf_chunk_read

  use hdf5

  use constants, only: CUSTOM_REAL,SIZE_REAL,NGLLX,NGLLY,NGLLZ,MAX_STRING_LEN

  use gf3d, only: t_gfdb,gf_open,gf_close,gf_errmsg,GF_OK

  use gf_element_io, only: gf_element_path,gf_read_element_displ

  use gf_hdf5_read, only: GF_HID,gf_h5_file_open,gf_h5_file_close,gf_h5_read_displ_chunks, &
                          GF_H5_ORDER_DISPL,GF_H5_ORDER_ATM

  use gf_manufactured, only: gf_report_true

  implicit none

  integer, parameter :: NC = 3
  integer, parameter :: NM = NC*NGLLX*NGLLY*NGLLZ
  integer, parameter :: NPREFIX = 4
  real(kind=CUSTOM_REAL), parameter :: SENTINEL = -12345._CUSTOM_REAL

  ! a real $GF3D_TEST_GFDB has thousands of elements; the fixture has 3
  integer, parameter :: NELEM_MAX = 3

  ! section 2's file: a short record, chunk length 7, so 20 = 2*7 + 6
  integer, parameter :: NT3 = 20, NCHUNK3 = 7

  character(len=MAX_STRING_LEN) :: dbpath,expect,filename
  type(t_gfdb) :: db
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: ref,displ
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: buf
  integer, dimension(NPREFIX) :: prefixes
  integer :: nt,ip,nt_out,ie,ista,iord,ierr,nfail,nbad,nroute,ncheck,nelem
  integer(kind=GF_HID) :: fid
  logical :: raw,want_raw

  nfail = 0

  if (command_argument_count() < 2) then
    write(*,'(a)') 'usage: test_gf_chunk_read <GFDB> <raw|fallback|any>'
    stop 1
  endif
  call get_command_argument(1,dbpath)
  call get_command_argument(2,expect)
  if (trim(expect) /= 'raw' .and. trim(expect) /= 'fallback' .and. trim(expect) /= 'any') then
    write(*,'(a)') 'usage: test_gf_chunk_read <GFDB> <raw|fallback|any>'
    stop 1
  endif
  ! a CUSTOM_REAL = 8 build never reads raw, whatever the file
  want_raw = (trim(expect) == 'raw' .and. CUSTOM_REAL == SIZE_REAL)

  write(*,'(a)') ''
  write(*,'(a)') 'test_gf_chunk_read'
  write(*,'(a)') ''
  write(*,'(a,a)') '  database = ',trim(dbpath)
  write(*,'(a,a)') '  route    = ',trim(expect)
  write(*,'(a)') ''

  call gf_open(dbpath,db,ierr,check_completion=.false.)
  if (ierr /= GF_OK) call die('could not open the database')

  nt = db%nt_subsampled
  prefixes = (/ 1, 7, 8, nt /)

  nelem = min(db%nelem,NELEM_MAX)

  allocate(ref(NC,NC,NGLLX,NGLLY,NGLLZ,nt),displ(NC,NC,NGLLX,NGLLY,NGLLZ,nt), &
           buf(NC*NM*nt+1),stat=ierr)
  if (ierr /= 0) call die('could not allocate')

  !--------------------------------------------------------------------
  ! 0. the reader every extraction uses
  !--------------------------------------------------------------------

  write(*,'(a)') '0. gf_read_element_displ against h5dread_f'

  nbad = 0
  ncheck = 0
  do ie = 1,nelem
    do ista = 1,db%nstations
      call reference(ie,ista,ref)
      call gf_read_element_displ(db,ie,ista,displ,ierr)
      if (ierr /= GF_OK) call die('gf_read_element_displ failed')
      ncheck = ncheck + 1
      if (any(displ /= ref)) then
        nbad = nbad + 1
        if (nbad == 1) write(*,'(a,i0,a,i0)') '     mismatch: first at element ',ie,', station ',ista
      endif
    enddo
  enddo
  call gf_report_true('gf_read_element_displ: bitwise, every station',nbad == 0 .and. ncheck > 0,nfail)

  !--------------------------------------------------------------------
  ! 1. both orders, four prefixes, the route
  !--------------------------------------------------------------------

  write(*,'(a)') '1. gf_h5_read_displ_chunks against h5dread_f'

  do iord = 1,2
    do ip = 1,NPREFIX
      nt_out = prefixes(ip)
      if (nt_out > nt) cycle
      nbad = 0
      nroute = 0
      ncheck = 0
      do ie = 1,nelem
        do ista = 1,db%nstations
          call reference(ie,ista,ref)

          call gf_element_path(db,ie,ista,filename,ierr)
          if (ierr /= GF_OK) call die('gf_element_path failed')
          call gf_h5_file_open(filename,fid,ierr)
          if (ierr /= GF_OK) call die('could not open '//trim(filename))

          buf(:) = SENTINEL
          call gf_h5_read_displ_chunks(fid,'displacement',nt,nt_out,iord,buf,raw,ierr)
          if (ierr /= GF_OK) call die('gf_h5_read_displ_chunks failed on '//trim(filename))
          call gf_h5_file_close(fid,ierr)

          ncheck = ncheck + 1
          if (.not. same(iord,ref,nt,nt_out,buf)) then
            nbad = nbad + 1
            if (nbad == 1) write(*,'(a,i0,a,i0)') '     mismatch: first at element ',ie,', station ',ista
          endif
          if (trim(expect) /= 'any' .and. (raw .neqv. want_raw)) then
            nroute = nroute + 1
            if (nroute == 1) write(*,'(a,l1,a,i0,a,i0)') '     mismatch: raw = ',raw, &
                                    ' at element ',ie,', station ',ista
          endif
        enddo
      enddo
      call gf_report_true(label(iord,nt_out,'bitwise, buffer end untouched'), &
                          nbad == 0 .and. ncheck > 0,nfail)
      ! says which route is asserted: on a CUSTOM_REAL = 8 build a chunked
      ! database is read through h5dread_f too
      if (trim(expect) /= 'any') &
        call gf_report_true(label(iord,nt_out,'route is '//route_name(want_raw)),nroute == 0,nfail)
    enddo
  enddo

  call gf_close(db)

  !--------------------------------------------------------------------
  ! 2. a chunk never written
  !--------------------------------------------------------------------

  write(*,'(a)') '2. a file with a chunk never written'

  call written_file_check(.true.)
  call written_file_check(.false.)

  write(*,'(a)') ''
  if (nfail /= 0) then
    write(*,'(a,i0,a)') '  ',nfail,' check(s) FAILED'
    stop 1
  endif
  write(*,'(a)') '  all checks passed'

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine reference(ie,is,d)

! the oracle: h5dread_f of the whole dataset into the dataset's own order,
! asking for the buffer's type as the library always has

  implicit none
  integer, intent(in) :: ie,is
  real(kind=CUSTOM_REAL), dimension(NC,NC,NGLLX,NGLLY,NGLLZ,nt), intent(out) :: d

  character(len=MAX_STRING_LEN) :: fn
  integer(kind=HID_T) :: f,ds,mem_type
  integer(kind=HSIZE_T), dimension(6) :: dims
  integer :: ier

  call gf_element_path(db,ie,is,fn,ier)
  if (ier /= GF_OK) call die('gf_element_path failed')
  if (CUSTOM_REAL == SIZE_REAL) then
    mem_type = H5T_NATIVE_REAL
  else
    mem_type = H5T_NATIVE_DOUBLE
  endif
  dims = (/ int(NC,HSIZE_T),int(NC,HSIZE_T),int(NGLLX,HSIZE_T),int(NGLLY,HSIZE_T), &
            int(NGLLZ,HSIZE_T),int(nt,HSIZE_T) /)
  call h5fopen_f(trim(fn),H5F_ACC_RDONLY_F,f,ier)
  if (ier /= 0) call die('could not open '//trim(fn))
  call h5dopen_f(f,'displacement',ds,ier)
  if (ier /= 0) call die('no displacement in '//trim(fn))
  call h5dread_f(ds,mem_type,d,dims,ier)
  if (ier /= 0) call die('h5dread_f failed on '//trim(fn))
  call h5dclose_f(ds,ier)
  call h5fclose_f(f,ier)

  end subroutine reference

!
!-------------------------------------------------------------------------------------------------
!

  logical function same(iorder,d,ntd,nto,b)

! b holds what the order says d(:,:,:,:,:,1:nto) is, and nothing after it

  implicit none
  integer, intent(in) :: iorder,ntd,nto
  real(kind=CUSTOM_REAL), dimension(NC,NM,ntd), intent(in) :: d
  real(kind=CUSTOM_REAL), dimension(:), intent(in) :: b

  integer :: ia,it,base

  same = (b(NC*NM*nto+1) == SENTINEL)
  if (iorder == GF_H5_ORDER_DISPL) then
    ! displ(a,m,t): the dataset's own order
    same = same .and. all(b(1:NC*NM*nto) == reshape(d(:,:,1:nto),(/ NC*NM*nto /)))
  else
    ! blk(m,t,a)
    do ia = 1,NC
      do it = 1,nto
        base = ((ia-1)*nto + (it-1))*NM
        same = same .and. all(b(base+1:base+NM) == d(ia,:,it))
      enddo
    enddo
  endif

  end function same

!
!-------------------------------------------------------------------------------------------------
!

  function route_name(r) result(s)

  implicit none
  logical, intent(in) :: r
  character(len=8) :: s

  if (r) then
    s = 'raw'
  else
    s = 'fallback'
  endif

  end function route_name

!
!-------------------------------------------------------------------------------------------------
!

  function label(iorder,nto,what) result(s)

  implicit none
  integer, intent(in) :: iorder,nto
  character(len=*), intent(in) :: what
  character(len=96) :: s

  if (iorder == GF_H5_ORDER_DISPL) then
    write(s,'(a,i4,a,a)') 'displ(a,p,i,j,k,t), nt_out =',nto,': ',what
  else
    write(s,'(a,i4,a,a)') 'blk(m,t,a),         nt_out =',nto,': ',what
  endif

  end function label

!
!-------------------------------------------------------------------------------------------------
!

  subroutine written_file_check(complete)

! writes displacement(3,3,5,5,5,NT3) chunked as the solver does, with the
! last chunk of force component 3 left unwritten unless `complete`, reads it
! back through gf_h5_read_displ_chunks in both orders, and compares with what
! was written -- the fill value 0 where nothing was

  implicit none
  logical, intent(in) :: complete

  real, dimension(NC,NC,NGLLX,NGLLY,NGLLZ,NT3) :: d,want
  real(kind=CUSTOM_REAL), dimension(NC*NM*NT3+1) :: b
  integer(kind=HID_T) :: f,s,ds,p,ms,fs
  integer(kind=HSIZE_T), dimension(6) :: dims,cdims,cnt,off
  integer :: ia,i,ier,iorder,nwritten
  character(len=*), parameter :: fname = './bin/test_gf_chunk_read.h5'
  logical :: r

  ! distinct, exactly representable, never the fill value
  d = reshape((/ (real(i), i = 1,size(d)) /),shape(d))

  ! (NT3/NCHUNK3)*NCHUNK3 = 14 samples, two whole chunks: the third, partial
  ! chunk of component 3 is never written
  nwritten = (NT3/NCHUNK3)*NCHUNK3
  want = d
  if (.not. complete) want(3,:,:,:,:,nwritten+1:NT3) = 0.

  dims  = (/ int(NC,HSIZE_T),int(NC,HSIZE_T),int(NGLLX,HSIZE_T),int(NGLLY,HSIZE_T), &
             int(NGLLZ,HSIZE_T),int(NT3,HSIZE_T) /)
  cdims = (/ 1_HSIZE_T,int(NC,HSIZE_T),int(NGLLX,HSIZE_T),int(NGLLY,HSIZE_T), &
             int(NGLLZ,HSIZE_T),int(NCHUNK3,HSIZE_T) /)

  call h5fcreate_f(fname,H5F_ACC_TRUNC_F,f,ier)
  if (ier /= 0) call die('could not create '//fname)
  call h5screate_simple_f(6,dims,s,ier)
  call h5pcreate_f(H5P_DATASET_CREATE_F,p,ier)
  call h5pset_chunk_f(p,6,cdims,ier)
  call h5pset_fill_value_f(p,H5T_NATIVE_REAL,0.,ier)
  call h5dcreate_f(f,'displacement',H5T_NATIVE_REAL,s,ds,ier,dcpl_id=p)
  call h5pclose_f(p,ier)
  ! one force component at a time, as the solver writes
  do ia = 1,NC
    cnt = dims
    cnt(1) = 1
    if (ia == NC .and. .not. complete) cnt(6) = int(nwritten,HSIZE_T)
    off(:) = 0
    off(1) = int(ia-1,HSIZE_T)
    call h5dget_space_f(ds,fs,ier)
    call h5sselect_hyperslab_f(fs,H5S_SELECT_SET_F,off,cnt,ier)
    call h5screate_simple_f(6,cnt,ms,ier)
    call h5dwrite_f(ds,H5T_NATIVE_REAL,d(ia:ia,:,:,:,:,1:cnt(6)),cnt,ier, &
                    mem_space_id=ms,file_space_id=fs)
    if (ier /= 0) call die('could not write '//fname)
    call h5sclose_f(ms,ier)
    call h5sclose_f(fs,ier)
  enddo
  call h5dclose_f(ds,ier)
  call h5sclose_f(s,ier)
  call h5fclose_f(f,ier)

  call gf_h5_file_open(fname,fid,ier)
  if (ier /= GF_OK) call die('could not reopen '//fname)
  do iorder = 1,2
    b(:) = SENTINEL
    call gf_h5_read_displ_chunks(fid,'displacement',NT3,NT3,iorder,b,r,ier)
    if (ier /= GF_OK) call die('gf_h5_read_displ_chunks failed on '//fname)
    if (complete) then
      call gf_report_true(label(iorder,NT3,'complete file: bitwise'), &
                          same(iorder,real(want,CUSTOM_REAL),NT3,NT3,b),nfail)
      call gf_report_true(label(iorder,NT3,'complete file: raw route'), &
                          r .eqv. (CUSTOM_REAL == SIZE_REAL),nfail)
    else
      call gf_report_true(label(iorder,NT3,'unwritten chunk: bitwise, fill value'), &
                          same(iorder,real(want,CUSTOM_REAL),NT3,NT3,b),nfail)
      call gf_report_true(label(iorder,NT3,'unwritten chunk: fallback route'), &
                          .not. r,nfail)
    endif
  enddo
  call gf_h5_file_close(fid,ier)

  open(unit=11,file=fname,status='old',iostat=ier)
  if (ier == 0) close(11,status='delete')

  end subroutine written_file_check

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

  end program test_gf_chunk_read
