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
!---- Opening, validating and inspecting a Green function database.
!----
!---- Layout written by src/specfem3D/green_function_*.F90:
!----
!----   {GFDB}/mesh_info.h5                       simulation-wide parameters
!----   {GFDB}/elements/{morton_hex}/coordinates.h5     xyz(3,5,5,5) + centroid
!----   {GFDB}/elements/{morton_hex}/{net}.{sta}.h5     displacement(3,3,5,5,5,nt_sub)
!----   {GFDB}/stations/{net}.{sta}.h5            station metadata + stf
!----
!---- and, only after utils/green_function/gf_build_manifest.py has been
!---- run over it:
!----
!----   {GFDB}/centroids.bin                      nelem * (int64 morton + 3 * float64)
!----   {GFDB}/manifest.csv                       morton_hex,cx,cy,cz
!----
!---- The element index therefore has to be recoverable without those two
!---- files: a database straight out of the solver has neither. (The
!---- solver contains gf_write_manifest() in green_function_morton.F90,
!---- which would write its own copies under {GFDB}/elements/, but nothing
!---- calls it.) gf_open() tries, in order:
!----
!----   1. {GFDB}/centroids.bin        fast path, one read
!----   2. {GFDB}/manifest.csv         text, one read
!----   3. {GFDB}/elements/centroids.bin, {GFDB}/elements/manifest.csv
!----                                  the paths gf_write_manifest() uses
!----   4. a scan of the {GFDB}/elements/ subdirectories
!----                                  always works, costs nelem file opens
!----
!---- Every route yields the same index in the same order (ascending
!---- Morton code), so nothing downstream depends on which one ran.
!----
!---- Precision policy: read float32, compute and store double.
!----

  module gf_database

  use gf_par

  use gf_hdf5_read, only: GF_HID,gf_h5_init,gf_h5_file_open,gf_h5_file_close, &
                          gf_h5_has_attr,gf_h5_has_dset, &
                          gf_h5_read_attr_i,gf_h5_read_attr_d, &
                          gf_h5_dset_dims,gf_h5_read_1d_d,gf_h5_read_2d_i

  use gf_shared_params, only: gf_init_shared_params

  implicit none

  private

  public :: gf_open
  public :: gf_close
  public :: gf_load_topo
  public :: gf_check_completion
  public :: gf_print_info

  ! unit numbers used for the plain-text / stream reads
  integer, parameter :: IIN_GF = 71

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_open(path,db,ierr,check_completion)

! opens a Green function database and reads all of its metadata
!
! The contract is "open once, extract many": this reads the element index,
! every station's metadata and source time function, and validates the
! database against the current build — but none of the bulk arrays. The
! 58 MB topography grid is left to gf_load_topo().
!
! check_completion (default .true.) additionally opens every
! element-station file to verify its computed_ALL flag. That is
! nelem * nstations file opens; later stages that have already validated a
! handle may pass .false.

  implicit none

  character(len=*), intent(in) :: path
  type(t_gfdb), intent(inout) :: db
  integer, intent(out) :: ierr
  logical, intent(in), optional :: check_completion

  ! local parameters
  logical :: do_check,exists
  integer :: nincomplete

  ! starts from a clean handle
  call gf_close(db)

  if (len_trim(path) == 0) then
    call gf_set_error(ierr,GF_ERR_ARG,'empty database path')
    return
  endif
  if (len_trim(path) > len(db%path)) then
    call gf_set_error(ierr,GF_ERR_ARG,'database path is too long')
    return
  endif

  db%path = trim(path)

  ! strips a trailing slash so that path joining below is unambiguous
  do while (len_trim(db%path) > 1)
    if (db%path(len_trim(db%path):len_trim(db%path)) /= '/') exit
    db%path = db%path(1:len_trim(db%path)-1)
  enddo

  call gf_dir_exists(db%path,exists)
  if (.not. exists) then
    call gf_set_error(ierr,GF_ERR_NO_PATH,'no such database directory: '//trim(db%path))
    return
  endif

  call gf_h5_init(ierr)
  if (ierr /= GF_OK) return

  call gf_read_mesh_info(db,ierr)
  if (ierr /= GF_OK) return

  call gf_read_element_index(db,ierr)
  if (ierr /= GF_OK) return

  call gf_read_stations(db,ierr)
  if (ierr /= GF_OK) return

  call gf_validate(db,ierr)
  if (ierr /= GF_OK) return

  ! the reused src/shared/ routines read module state that nothing else
  ! initialises in a standalone program
  call gf_init_shared_params(db,ierr)
  if (ierr /= GF_OK) return

  do_check = .true.
  if (present(check_completion)) do_check = check_completion

  if (do_check) then
    call gf_check_completion(db,nincomplete,ierr)
    if (ierr /= GF_OK) return
  endif

  db%is_open = .true.
  ierr = GF_OK

  end subroutine gf_open

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_close(db)

! releases everything held by a handle
!
! Safe to call on a handle that was never opened, and safe to call twice.

  implicit none

  type(t_gfdb), intent(inout) :: db

  ! local parameters
  integer :: i

  if (allocated(db%rspl)) deallocate(db%rspl)
  if (allocated(db%ellipicity_spline)) deallocate(db%ellipicity_spline)
  if (allocated(db%ellipicity_spline2)) deallocate(db%ellipicity_spline2)
  if (allocated(db%ibathy_topo)) deallocate(db%ibathy_topo)
  if (allocated(db%morton)) deallocate(db%morton)
  if (allocated(db%morton_hex)) deallocate(db%morton_hex)
  if (allocated(db%centroid)) deallocate(db%centroid)

  if (allocated(db%stations)) then
    do i = 1,size(db%stations)
      if (allocated(db%stations(i)%stf)) deallocate(db%stations(i)%stf)
    enddo
    deallocate(db%stations)
  endif

  db%is_open = .false.
  db%path = ''
  db%topo_loaded = .false.
  db%nelem = 0
  db%nstations = 0
  db%nspl = 0
  db%index_source = ''

  end subroutine gf_close

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_read_mesh_info(db,ierr)

! reads {GFDB}/mesh_info.h5

  implicit none

  type(t_gfdb), intent(inout) :: db
  integer, intent(out) :: ierr

  ! local parameters
  character(len=MAX_STRING_LEN) :: filename
  integer(kind=GF_HID) :: fid
  integer :: ival,ierr2

  filename = trim(db%path)//'/mesh_info.h5'

  call gf_h5_file_open(filename,fid,ierr)
  if (ierr /= GF_OK) then
    call gf_set_error(ierr,GF_ERR_NO_FILE,'cannot open '//trim(filename))
    return
  endif

  !--- time axis ---
  call gf_h5_read_attr_d(fid,'dt',db%dt,ierr)
  if (ierr /= GF_OK) goto 99
  call gf_h5_read_attr_d(fid,'t0',db%t0,ierr)
  if (ierr /= GF_OK) goto 99
  call gf_h5_read_attr_i(fid,'nstep',db%nstep,ierr)
  if (ierr /= GF_OK) goto 99
  call gf_h5_read_attr_i(fid,'nt_subsampled',db%nt_subsampled,ierr)
  if (ierr /= GF_OK) goto 99
  call gf_h5_read_attr_i(fid,'subsample_step',db%subsample_step,ierr)
  if (ierr /= GF_OK) goto 99

  !--- storage parameters ---
  call gf_h5_read_attr_i(fid,'buffer_size',db%buffer_size,ierr)
  if (ierr /= GF_OK) goto 99
  call gf_h5_read_attr_i(fid,'neighbor_shells',db%neighbor_shells,ierr)
  if (ierr /= GF_OK) goto 99

  !--- element discretisation ---
  ! note: nglly/ngllz were added to the writer after the shipped example
  !       databases were generated, so fall back to ngllx when absent
  call gf_h5_read_attr_i(fid,'ngllx',db%ngllx,ierr)
  if (ierr /= GF_OK) goto 99

  db%nglly = db%ngllx
  db%ngllz = db%ngllx
  if (gf_h5_has_attr(fid,'nglly')) then
    call gf_h5_read_attr_i(fid,'nglly',db%nglly,ierr)
    if (ierr /= GF_OK) goto 99
  endif
  if (gf_h5_has_attr(fid,'ngllz')) then
    call gf_h5_read_attr_i(fid,'ngllz',db%ngllz,ierr)
    if (ierr /= GF_OK) goto 99
  endif

  !--- planet constants ---
  call gf_h5_read_attr_d(fid,'scale_displ',db%scale_displ,ierr)
  if (ierr /= GF_OK) goto 99
  call gf_h5_read_attr_d(fid,'R_PLANET',db%R_PLANET,ierr)
  if (ierr /= GF_OK) goto 99

  ! not written by the current writer; read it if a later one adds it, so
  ! that the library is not silently Earth-only (see gf_shared_params.F90)
  db%RHOAV = 0.d0
  if (gf_h5_has_attr(fid,'RHOAV')) then
    call gf_h5_read_attr_d(fid,'RHOAV',db%RHOAV,ierr)
    if (ierr /= GF_OK) goto 99
  endif

  !--- simulation flags ---
  call gf_read_flag(fid,'TOPOGRAPHY',db%topography,ierr)
  if (ierr /= GF_OK) goto 99
  call gf_read_flag(fid,'ELLIPTICITY',db%ellipticity,ierr)
  if (ierr /= GF_OK) goto 99
  call gf_read_flag(fid,'rotation',db%rotation,ierr)
  if (ierr /= GF_OK) goto 99
  call gf_read_flag(fid,'attenuation',db%attenuation,ierr)
  if (ierr /= GF_OK) goto 99
  call gf_read_flag(fid,'gravity',db%gravity,ierr)
  if (ierr /= GF_OK) goto 99

  !--- topography grid geometry (the grid itself is loaded lazily) ---
  if (db%topography) then
    call gf_h5_read_attr_i(fid,'NX_BATHY',db%NX_BATHY,ierr)
    if (ierr /= GF_OK) goto 99
    call gf_h5_read_attr_i(fid,'NY_BATHY',db%NY_BATHY,ierr)
    if (ierr /= GF_OK) goto 99
    call gf_h5_read_attr_d(fid,'RESOLUTION_TOPO_FILE',db%RESOLUTION_TOPO_FILE,ierr)
    if (ierr /= GF_OK) goto 99

    if (.not. gf_h5_has_dset(fid,'ibathy_topo')) then
      call gf_set_error(ierr,GF_ERR_FORMAT, &
        'TOPOGRAPHY is set but mesh_info.h5 has no ibathy_topo dataset')
      goto 99
    endif
  endif

  !--- ellipticity splines (small: read them now) ---
  if (db%ellipticity) then
    call gf_h5_read_attr_i(fid,'nspl',db%nspl,ierr)
    if (ierr /= GF_OK) goto 99

    if (db%nspl <= 0) then
      call gf_set_error(ierr,GF_ERR_FORMAT,'ELLIPTICITY is set but nspl <= 0')
      goto 99
    endif

    allocate(db%rspl(db%nspl),db%ellipicity_spline(db%nspl), &
             db%ellipicity_spline2(db%nspl),stat=ival)
    if (ival /= 0) then
      call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate the ellipticity splines')
      goto 99
    endif

    call gf_h5_read_1d_d(fid,'rspl',db%nspl,db%rspl,ierr)
    if (ierr /= GF_OK) goto 99
    call gf_h5_read_1d_d(fid,'ellipicity_spline',db%nspl,db%ellipicity_spline,ierr)
    if (ierr /= GF_OK) goto 99
    call gf_h5_read_1d_d(fid,'ellipicity_spline2',db%nspl,db%ellipicity_spline2,ierr)
    if (ierr /= GF_OK) goto 99
  endif

  call gf_h5_file_close(fid,ierr)
  return

99 continue
  call gf_h5_file_close(fid,ierr2)

  end subroutine gf_read_mesh_info

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_read_flag(fid,name,flag,ierr)

! reads an integer attribute written as merge(1,0,logical)

  implicit none

  integer(kind=GF_HID), intent(in) :: fid
  character(len=*), intent(in) :: name
  logical, intent(out) :: flag
  integer, intent(out) :: ierr

  ! local parameters
  integer :: ival

  flag = .false.

  call gf_h5_read_attr_i(fid,name,ival,ierr)
  if (ierr /= GF_OK) return

  flag = (ival /= 0)

  end subroutine gf_read_flag

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_read_element_index(db,ierr)

! builds the element index: Morton codes and centroids, ascending
!
! See the header of this file for the order in which the four possible
! sources are tried.

  implicit none

  type(t_gfdb), intent(inout) :: db
  integer, intent(out) :: ierr

  ! local parameters
  character(len=MAX_STRING_LEN) :: elemdir
  logical :: exists,ok
  integer :: nelem_csv

  elemdir = trim(db%path)//'/elements'

  call gf_dir_exists(elemdir,exists)
  if (.not. exists) then
    call gf_set_error(ierr,GF_ERR_NO_PATH,'no such directory: '//trim(elemdir))
    return
  endif

  ! the reference count, when a manifest is available at all: the plan's
  ! "nelem must equal wc -l manifest.csv minus 1"
  nelem_csv = -1
  call gf_count_manifest(trim(db%path)//'/manifest.csv',nelem_csv)
  if (nelem_csv < 0) call gf_count_manifest(trim(elemdir)//'/manifest.csv',nelem_csv)

  ! 1. {GFDB}/centroids.bin
  call gf_try_centroids_bin(trim(db%path)//'/centroids.bin',nelem_csv,db,ok)
  if (ok) then
    db%index_source = 'centroids.bin'
    goto 10
  endif

  ! 2. {GFDB}/manifest.csv
  call gf_try_manifest_csv(trim(db%path)//'/manifest.csv',db,ok)
  if (ok) then
    db%index_source = 'manifest.csv'
    goto 10
  endif

  ! 3. the paths gf_write_manifest() in the solver would use
  call gf_try_centroids_bin(trim(elemdir)//'/centroids.bin',nelem_csv,db,ok)
  if (ok) then
    db%index_source = 'centroids.bin'
    goto 10
  endif

  call gf_try_manifest_csv(trim(elemdir)//'/manifest.csv',db,ok)
  if (ok) then
    db%index_source = 'manifest.csv'
    goto 10
  endif

  ! 4. scan the element directories themselves
  call gf_scan_element_dirs(db,ierr)
  if (ierr /= GF_OK) return
  db%index_source = 'element scan'

10 continue

  if (db%nelem <= 0) then
    call gf_set_error(ierr,GF_ERR_FORMAT,'the database contains no elements')
    return
  endif

  ! cross-check against the manifest line count when there was one
  if (nelem_csv >= 0 .and. nelem_csv /= db%nelem) then
    call gf_set_error(ierr,GF_ERR_FORMAT, &
      'element count disagrees with manifest.csv; the database index is stale, ' // &
      'please re-run utils/green_function/gf_build_manifest.py')
    return
  endif

  ierr = GF_OK

  end subroutine gf_read_element_index

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_try_centroids_bin(filename,nelem_expected,db,ok)

! reads the packed element index: nelem records of int64 Morton + 3 float64
!
! Written by Python as little-endian '<Qddd' (gf_build_manifest.py) and by
! Fortran stream I/O (green_function_morton.F90); both are native on every
! platform this runs on. Rather than assume that, the records are sanity
! checked after the read: non-dimensional centroids lie inside the unit
! ball and Morton codes ascend. A byte-swapped or otherwise unreadable file
! fails those and we fall through to the text manifest, which cannot be
! misread.

  implicit none

  character(len=*), intent(in) :: filename
  integer, intent(in) :: nelem_expected
  type(t_gfdb), intent(inout) :: db
  logical, intent(out) :: ok

  ! local parameters
  integer(kind=8) :: fsize,mval
  integer :: nrec,i,ier
  double precision :: cx,cy,cz
  logical :: exists

  ok = .false.

  call gf_file_exists(filename,exists)
  if (.not. exists) return

  inquire(file=trim(filename),size=fsize)
  if (fsize <= 0) return

  ! 8 bytes of Morton code + 3 * 8 bytes of centroid
  if (mod(fsize,32_8) /= 0) return

  nrec = int(fsize / 32_8)
  if (nrec <= 0) return
  if (nelem_expected >= 0 .and. nrec /= nelem_expected) return

  open(unit=IIN_GF,file=trim(filename),form='unformatted',access='stream', &
       status='old',action='read',iostat=ier)
  if (ier /= 0) return

  allocate(db%morton(nrec),db%morton_hex(nrec),db%centroid(3,nrec),stat=ier)
  if (ier /= 0) then
    close(IIN_GF)
    return
  endif

  do i = 1,nrec
    read(IIN_GF,iostat=ier) mval,cx,cy,cz
    if (ier /= 0) goto 99

    db%morton(i) = mval
    db%centroid(1,i) = cx
    db%centroid(2,i) = cy
    db%centroid(3,i) = cz
    write(db%morton_hex(i),'(Z16.16)') mval
  enddo

  close(IIN_GF)

  ! raw bytes: apply the range test too, so that an unreadable file makes us
  ! fall through to the text manifest rather than trusting nonsense
  if (.not. gf_index_is_sane(db,nrec,check_range=.true.)) goto 98

  db%nelem = nrec
  ok = .true.
  return

99 continue
  close(IIN_GF)
98 continue
  deallocate(db%morton,db%morton_hex,db%centroid)

  end subroutine gf_try_centroids_bin

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_try_manifest_csv(filename,db,ok)

! reads the text element index: a header line then 'morton_hex,cx,cy,cz'

  implicit none

  character(len=*), intent(in) :: filename
  type(t_gfdb), intent(inout) :: db
  logical, intent(out) :: ok

  ! local parameters
  character(len=MAX_STRING_LEN) :: line
  integer :: nrec,i,ier,ic1,ic2,ic3
  logical :: exists

  ok = .false.

  call gf_file_exists(filename,exists)
  if (.not. exists) return

  call gf_count_manifest(filename,nrec)
  if (nrec <= 0) return

  open(unit=IIN_GF,file=trim(filename),form='formatted',status='old',action='read',iostat=ier)
  if (ier /= 0) return

  ! header
  read(IIN_GF,'(a)',iostat=ier) line
  if (ier /= 0) then
    close(IIN_GF)
    return
  endif

  allocate(db%morton(nrec),db%morton_hex(nrec),db%centroid(3,nrec),stat=ier)
  if (ier /= 0) then
    close(IIN_GF)
    return
  endif

  do i = 1,nrec
    read(IIN_GF,'(a)',iostat=ier) line
    if (ier /= 0) goto 99

    ic1 = index(line,',')
    if (ic1 <= 1) goto 99
    ic2 = index(line(ic1+1:),',')
    if (ic2 <= 0) goto 99
    ic2 = ic1 + ic2
    ic3 = index(line(ic2+1:),',')
    if (ic3 <= 0) goto 99
    ic3 = ic2 + ic3

    if (ic1 - 1 /= GF_MORTON_HEXLEN) goto 99
    db%morton_hex(i) = line(1:ic1-1)

    read(db%morton_hex(i),'(Z16)',iostat=ier) db%morton(i)
    if (ier /= 0) goto 99

    read(line(ic1+1:ic2-1),*,iostat=ier) db%centroid(1,i)
    if (ier /= 0) goto 99
    read(line(ic2+1:ic3-1),*,iostat=ier) db%centroid(2,i)
    if (ier /= 0) goto 99
    read(line(ic3+1:),*,iostat=ier) db%centroid(3,i)
    if (ier /= 0) goto 99
  enddo

  close(IIN_GF)

  ! as for centroids.bin: a failure here falls through to the element scan
  if (.not. gf_index_is_sane(db,nrec,check_range=.true.)) then
    deallocate(db%morton,db%morton_hex,db%centroid)
    return
  endif

  db%nelem = nrec
  ok = .true.
  return

99 continue
  close(IIN_GF)
  deallocate(db%morton,db%morton_hex,db%centroid)

  end subroutine gf_try_manifest_csv

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_scan_element_dirs(db,ierr)

! builds the element index straight from {GFDB}/elements/, by reading the
! cx/cy/cz attributes of each coordinates.h5
!
! This is the only route that works on a database that has never been
! post-processed, which is what the solver actually produces today.
! Directory names are the Morton codes in Z16.16 format, and sorting them
! as strings gives the same ascending order as sorting the codes: they are
! fixed width, upper case, and never exceed 2^63 (21 bits per axis,
! interleaved into bits 0..62).

  implicit none

  type(t_gfdb), intent(inout) :: db
  integer, intent(out) :: ierr

  ! local parameters
  character(len=MAX_STRING_LEN) :: elemdir,filename
  character(len=GF_MORTON_HEXLEN), dimension(:), allocatable :: names
  integer(kind=GF_HID) :: fid
  integer :: nfound,i,ier,ierr2

  elemdir = trim(db%path)//'/elements'

  call gf_list_dir(elemdir,'',.true.,GF_MORTON_HEXLEN,names,nfound,ierr)
  if (ierr /= GF_OK) return

  if (nfound <= 0) then
    call gf_set_error(ierr,GF_ERR_FORMAT,'no element directories under '//trim(elemdir))
    return
  endif

  allocate(db%morton(nfound),db%morton_hex(nfound),db%centroid(3,nfound),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate the element index')
    return
  endif

  do i = 1,nfound
    db%morton_hex(i) = names(i)

    read(db%morton_hex(i),'(Z16)',iostat=ier) db%morton(i)
    if (ier /= 0) then
      call gf_set_error(ierr,GF_ERR_FORMAT, &
        'element directory name is not a 16-digit Morton code: '//trim(names(i)))
      return
    endif

    filename = trim(elemdir)//'/'//db%morton_hex(i)//'/coordinates.h5'

    call gf_h5_file_open(filename,fid,ierr)
    if (ierr /= GF_OK) return

    call gf_h5_read_attr_d(fid,'cx',db%centroid(1,i),ierr)
    if (ierr == GF_OK) call gf_h5_read_attr_d(fid,'cy',db%centroid(2,i),ierr)
    if (ierr == GF_OK) call gf_h5_read_attr_d(fid,'cz',db%centroid(3,i),ierr)

    call gf_h5_file_close(fid,ierr2)
    if (ierr /= GF_OK) return
  enddo

  ! last route, nothing to fall through to: structural invariants only, so an
  ! unusual but legitimate database is never rejected on its centroid range
  if (.not. gf_index_is_sane(db,nfound,check_range=.false.)) then
    call gf_set_error(ierr,GF_ERR_FORMAT,'implausible element centroids or Morton ordering')
    return
  endif

  db%nelem = nfound
  ierr = GF_OK

  end subroutine gf_scan_element_dirs

!
!-------------------------------------------------------------------------------------------------
!

  logical function gf_index_is_sane(db,n,check_range)

! sanity check on a freshly parsed element index
!
! Two tiers, because a failure means different things on different routes.
!
! Structural, always applied: centroids must be finite, and Morton codes
! must be non-negative and strictly ascending. Codes are 21 bits per axis
! interleaved into bits 0..62, so the sign bit is never set; both writers
! sort, and the directory scan sorts lexicographically over fixed-width
! upper-case hex, which is the same order. A duplicate would mean two
! element directories collided. None of this can misfire on a database that
! was written correctly.
!
! Range, applied only where falling through to another route is possible:
! each centroid component must satisfy |c| <= RADIUS_MAX. This is what
! catches a byte-swapped or otherwise unreadable centroids.bin, or a
! manifest.csv whose columns are not what we assume — in both cases we then
! read the index from somewhere safer. It is deliberately *not* applied to
! the element-directory scan, which has nowhere to fall through to: values
! read from HDF5 attributes cannot be misinterpreted the way raw bytes can,
! so an unusual-but-legitimate database must not be rejected there.
!
! On RADIUS_MAX: centroids are non-dimensionalised by R_PLANET, which is
! the *mean* radius, so a surface element on an elliptical, topographic
! planet does exceed 1. Earth's worst case is ellipticity (+1/(3*299.8) =
! +0.0011 at the equator) plus Everest (+9000/6371000 = +0.0014), so
! 1.0025 for a surface point and less for an element centroid, which sits
! half an element below it. Measured over both shipped example databases
! the largest centroid radius is 0.99963. RADIUS_MAX = 1.05 therefore has
! roughly twenty times the headroom that is physically reachable, and is
! already looser than the format's own domain: gf_morton_encode()
! quantises (c+1)/2 over [0,2^21-1] and *clamps*, so a database with
! |c| > 1 would have colliding Morton codes regardless of this reader.
!
! Note the order of the tests. This runs on bytes that came straight off
! disk and have not been validated by anything, so it must reach no
! floating-point operation that could trap: the project compiles with
! -ffpe-trap=invalid,zero,overflow, under which squaring a NaN or a 1e300
! read out of a corrupt file kills the process. Finiteness is therefore
! established first, by the non-signalling bit test in gf_is_finite(), and
! each component is bounded before anything is squared.

  implicit none

  type(t_gfdb), intent(in) :: db
  integer, intent(in) :: n
  logical, intent(in) :: check_range

  ! local parameters
  double precision, parameter :: RADIUS_MAX = 1.05d0
  double precision :: r2
  integer :: i,j

  gf_index_is_sane = .false.

  do i = 1,n
    do j = 1,3
      if (.not. gf_is_finite(db%centroid(j,i))) return
      if (check_range) then
        if (abs(db%centroid(j,i)) > RADIUS_MAX) return
      endif
    enddo

    if (check_range) then
      ! safe now: every component is finite and bounded by RADIUS_MAX
      r2 = db%centroid(1,i)**2 + db%centroid(2,i)**2 + db%centroid(3,i)**2
      if (r2 > RADIUS_MAX*RADIUS_MAX) return
    endif

    if (db%morton(i) < 0) return
    if (i > 1) then
      if (db%morton(i) <= db%morton(i-1)) return
    endif
  enddo

  gf_index_is_sane = .true.

  end function gf_index_is_sane

!
!-------------------------------------------------------------------------------------------------
!

  logical function gf_is_finite(x)

! true when x is neither NaN nor infinite
!
! Written as a bit test rather than with ieee_is_finite() because the
! intrinsic IEEE modules are not universally available — the Red Hat
! gcc-toolset gfortran packages, among others, ship without them — and
! because a comparison against a NaN is itself an invalid operation under
! -ffpe-trap=invalid. transfer() and iand() touch no floating-point unit.
!
! An IEEE-754 binary64 is NaN or infinite exactly when its 11 exponent
! bits (52..62) are all set. ISHFT is defined as a *logical* shift, so the
! sign bit does not smear into the result.

  implicit none

  double precision, intent(in) :: x

  ! local parameters
  integer(kind=8) :: bits
  integer(kind=8), parameter :: EXPONENT_MASK = 2047_8

  bits = transfer(x,bits)

  gf_is_finite = (iand(ishft(bits,-52),EXPONENT_MASK) /= EXPONENT_MASK)

  end function gf_is_finite

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_count_manifest(filename,nrec)

! counts the data lines of a manifest.csv, i.e. its line count minus the header
!
! returns -1 when the file does not exist or cannot be read

  implicit none

  character(len=*), intent(in) :: filename
  integer, intent(out) :: nrec

  ! local parameters
  character(len=MAX_STRING_LEN) :: line
  integer :: ier,n
  logical :: exists

  nrec = -1

  call gf_file_exists(filename,exists)
  if (.not. exists) return

  open(unit=IIN_GF,file=trim(filename),form='formatted',status='old',action='read',iostat=ier)
  if (ier /= 0) return

  ! header
  read(IIN_GF,'(a)',iostat=ier) line
  if (ier /= 0) then
    close(IIN_GF)
    return
  endif

  n = 0
  do
    read(IIN_GF,'(a)',iostat=ier) line
    if (ier /= 0) exit
    if (len_trim(line) == 0) cycle
    n = n + 1
  enddo

  close(IIN_GF)

  nrec = n

  end subroutine gf_count_manifest

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_read_stations(db,ierr)

! enumerates the .h5 files under {GFDB}/stations/ and reads each station's
! metadata and stf
!
! Nothing in the database records the station list, so this is a directory
! scan; see src/gf3d/gf_dirlist.c.

  implicit none

  type(t_gfdb), intent(inout) :: db
  integer, intent(out) :: ierr

  ! local parameters
  character(len=MAX_STRING_LEN) :: stadir,filename
  character(len=GF_STATION_ID_LEN+3), dimension(:), allocatable :: names
  integer(kind=GF_HID) :: fid
  integer :: nfound,i,ier,ierr2,idot,nlen
  integer :: ndims
  integer(kind=8), dimension(1) :: dims
  logical :: exists

  stadir = trim(db%path)//'/stations'

  call gf_dir_exists(stadir,exists)
  if (.not. exists) then
    call gf_set_error(ierr,GF_ERR_NO_PATH,'no such directory: '//trim(stadir))
    return
  endif

  call gf_list_dir(stadir,'.h5',.false.,GF_STATION_ID_LEN+3,names,nfound,ierr)
  if (ierr /= GF_OK) return

  if (nfound <= 0) then
    call gf_set_error(ierr,GF_ERR_FORMAT,'no station files under '//trim(stadir))
    return
  endif

  allocate(db%stations(nfound),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate the station table')
    return
  endif
  db%nstations = nfound

  do i = 1,nfound
    ! '{net}.{sta}.h5' -> network, station
    nlen = len_trim(names(i)) - 3
    if (nlen <= 0) then
      call gf_set_error(ierr,GF_ERR_FORMAT,'unexpected station file name: '//trim(names(i)))
      return
    endif

    db%stations(i)%id = names(i)(1:nlen)

    idot = index(db%stations(i)%id,'.')
    if (idot <= 1 .or. idot >= nlen) then
      call gf_set_error(ierr,GF_ERR_FORMAT, &
        'station file name is not NET.STA.h5: '//trim(names(i)))
      return
    endif
    db%stations(i)%network = db%stations(i)%id(1:idot-1)
    db%stations(i)%station = db%stations(i)%id(idot+1:nlen)

    filename = trim(stadir)//'/'//trim(names(i))

    call gf_h5_file_open(filename,fid,ierr)
    if (ierr /= GF_OK) return

    call gf_h5_read_attr_d(fid,'latitude',db%stations(i)%latitude,ierr)
    if (ierr == GF_OK) call gf_h5_read_attr_d(fid,'longitude',db%stations(i)%longitude,ierr)
    if (ierr == GF_OK) call gf_h5_read_attr_d(fid,'depth',db%stations(i)%depth,ierr)
    if (ierr == GF_OK) call gf_h5_read_attr_d(fid,'hdur',db%stations(i)%hdur,ierr)
    if (ierr == GF_OK) call gf_h5_read_attr_d(fid,'f_cutoff',db%stations(i)%f_cutoff,ierr)
    if (ierr == GF_OK) call gf_h5_read_attr_d(fid,'factor_force_source', &
                                              db%stations(i)%factor_force_source,ierr)
    if (ierr == GF_OK) call gf_h5_read_attr_d(fid,'time_shift',db%stations(i)%time_shift,ierr)
    if (ierr /= GF_OK) then
      call gf_h5_file_close(fid,ierr2)
      return
    endif

    ! source time function, on the unsubsampled solver time axis
    ndims = 1
    call gf_h5_dset_dims(fid,'stf',ndims,dims,ierr)
    if (ierr /= GF_OK) then
      call gf_h5_file_close(fid,ierr2)
      return
    endif

    if (ndims /= 1 .or. dims(1) /= int(db%nstep,kind=8)) then
      call gf_h5_file_close(fid,ierr2)
      call gf_set_error(ierr,GF_ERR_FORMAT, &
        'stf of station '//trim(db%stations(i)%id)//' does not have nstep samples')
      return
    endif

    allocate(db%stations(i)%stf(db%nstep),stat=ier)
    if (ier /= 0) then
      call gf_h5_file_close(fid,ierr2)
      call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate an stf array')
      return
    endif

    call gf_h5_read_1d_d(fid,'stf',db%nstep,db%stations(i)%stf,ierr)

    call gf_h5_file_close(fid,ierr2)
    if (ierr /= GF_OK) return
  enddo

  ierr = GF_OK

  end subroutine gf_read_stations

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_validate(db,ierr)

! hard failures, not warnings
!
! Anything here that were allowed through as a warning would corrupt memory
! or produce a plausible wrong answer further down.

  use constants, only: NGLLX,NGLLY,NGLLZ,NGNOD

  implicit none

  type(t_gfdb), intent(in) :: db
  integer, intent(out) :: ierr

  ! local parameters
  character(len=MAX_STRING_LEN) :: msg

  ! The GLL count is a compile-time parameter here (setup/constants.h:57) and
  ! the reused shape functions and anchor tables read it from `constants`. A
  ! database written with a different one would silently corrupt memory.
  if (db%ngllx /= NGLLX .or. db%nglly /= NGLLY .or. db%ngllz /= NGLLZ) then
    write(msg,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)') &
      'database was written with NGLL = ',db%ngllx,' x ',db%nglly,' x ',db%ngllz, &
      ' but this build has ',NGLLX,' x ',NGLLY,' x ',NGLLZ, &
      '; rebuild with a matching NGLLX in setup/constants.h.in'
    call gf_set_error(ierr,GF_ERR_MISMATCH,trim(msg))
    return
  endif

  ! the element geometry uses the 27 anchor nodes, matching the solver
  if (NGNOD /= 27) then
    call gf_set_error(ierr,GF_ERR_MISMATCH, &
      'this build has NGNOD /= 27; the Green function extraction needs the 27-anchor geometry')
    return
  endif

  if (db%nstep <= 0 .or. db%nt_subsampled <= 0 .or. db%subsample_step <= 0) then
    call gf_set_error(ierr,GF_ERR_FORMAT,'non-positive nstep / nt_subsampled / subsample_step')
    return
  endif

  if (db%nt_subsampled /= db%nstep / db%subsample_step) then
    call gf_set_error(ierr,GF_ERR_FORMAT, &
      'nt_subsampled does not equal nstep / subsample_step')
    return
  endif

  if (db%dt <= 0.d0) then
    call gf_set_error(ierr,GF_ERR_FORMAT,'non-positive time step dt')
    return
  endif

  if (db%R_PLANET <= 0.d0) then
    call gf_set_error(ierr,GF_ERR_FORMAT,'non-positive R_PLANET')
    return
  endif

  ierr = GF_OK

  end subroutine gf_validate

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_check_completion(db,nincomplete,ierr)

! verifies that every element-station file exists and has computed_ALL set
!
! This matters because the HDF5 fill value is 0: a force component that was
! never written reads back as a plausible-looking array of zeros rather than
! as an error. The writer sets computed_N / computed_E / computed_Z as each
! reciprocal run finishes and computed_ALL once all three have.

  implicit none

  type(t_gfdb), intent(in) :: db
  integer, intent(out) :: nincomplete
  integer, intent(out) :: ierr

  ! local parameters
  character(len=MAX_STRING_LEN) :: filename,first_bad
  integer(kind=GF_HID) :: fid
  integer :: ielem,ista,ival,ierr2
  logical :: exists

  nincomplete = 0
  first_bad = ''

  do ista = 1,db%nstations
    do ielem = 1,db%nelem

      filename = trim(db%path)//'/elements/'//db%morton_hex(ielem)//'/' &
                 //trim(db%stations(ista)%id)//'.h5'

      call gf_file_exists(filename,exists)
      if (.not. exists) then
        nincomplete = nincomplete + 1
        if (len_trim(first_bad) == 0) first_bad = trim(filename)//' (missing)'
        cycle
      endif

      call gf_h5_file_open(filename,fid,ierr)
      if (ierr /= GF_OK) then
        nincomplete = nincomplete + 1
        if (len_trim(first_bad) == 0) first_bad = trim(filename)//' (unreadable)'
        cycle
      endif

      ival = 0
      if (gf_h5_has_attr(fid,'computed_ALL')) then
        call gf_h5_read_attr_i(fid,'computed_ALL',ival,ierr)
        if (ierr /= GF_OK) ival = 0
      endif

      call gf_h5_file_close(fid,ierr2)

      if (ival == 0) then
        nincomplete = nincomplete + 1
        if (len_trim(first_bad) == 0) first_bad = trim(filename)//' (computed_ALL not set)'
      endif
    enddo
  enddo

  if (nincomplete > 0) then
    call gf_set_error(ierr,GF_ERR_INCOMPLETE, &
      'incomplete database: first offender '//trim(first_bad))
    return
  endif

  ierr = GF_OK

  end subroutine gf_check_completion

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_load_topo(db,ierr)

! loads the topography grid from mesh_info.h5
!
! Deliberately not done by gf_open(): ibathy_topo dominates the 58 MB of
! mesh_info.h5 and is only needed once a source position actually has to be
! placed relative to the surface. Calling this twice is a no-op.

  implicit none

  type(t_gfdb), intent(inout) :: db
  integer, intent(out) :: ierr

  ! local parameters
  character(len=MAX_STRING_LEN) :: filename
  integer(kind=GF_HID) :: fid
  integer :: ndims,ier,ierr2
  integer(kind=8), dimension(2) :: dims

  ierr = GF_OK

  if (db%topo_loaded) return

  if (.not. db%topography) then
    call gf_set_error(ierr,GF_ERR_ARG,'this database was written without topography')
    return
  endif

  filename = trim(db%path)//'/mesh_info.h5'

  call gf_h5_file_open(filename,fid,ierr)
  if (ierr /= GF_OK) return

  ndims = 2
  call gf_h5_dset_dims(fid,'ibathy_topo',ndims,dims,ierr)
  if (ierr /= GF_OK) then
    call gf_h5_file_close(fid,ierr2)
    return
  endif

  if (ndims /= 2 .or. dims(1) /= int(db%NX_BATHY,kind=8) &
                 .or. dims(2) /= int(db%NY_BATHY,kind=8)) then
    call gf_h5_file_close(fid,ierr2)
    call gf_set_error(ierr,GF_ERR_FORMAT, &
      'ibathy_topo does not have shape (NX_BATHY,NY_BATHY)')
    return
  endif

  allocate(db%ibathy_topo(db%NX_BATHY,db%NY_BATHY),stat=ier)
  if (ier /= 0) then
    call gf_h5_file_close(fid,ierr2)
    call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate the topography grid')
    return
  endif

  call gf_h5_read_2d_i(fid,'ibathy_topo',db%NX_BATHY,db%NY_BATHY,db%ibathy_topo,ierr)

  call gf_h5_file_close(fid,ierr2)
  if (ierr /= GF_OK) then
    deallocate(db%ibathy_topo)
    return
  endif

  db%topo_loaded = .true.

  end subroutine gf_load_topo

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_print_info(db,iunit,with_topo)

! prints every number in the handle, one 'key = value' per line
!
! Kept diffable on purpose: the Stage 1 acceptance test is a comparison of
! this output against `h5dump -A`, and this stays the first thing to run
! when a database will not load.

  implicit none

  type(t_gfdb), intent(inout) :: db
  integer, intent(in) :: iunit
  logical, intent(in) :: with_topo

  ! local parameters
  integer :: i,ierr
  double precision :: elevation

  write(iunit,'(a)') 'database'
  write(iunit,'(a,a)')      '  path                 = ',trim(db%path)
  write(iunit,'(a,a)')      '  element index from   = ',trim(db%index_source)
  write(iunit,'(a)') ''

  write(iunit,'(a)') 'time axis'
  write(iunit,'(a,es22.14)') '  dt                   = ',db%dt
  write(iunit,'(a,es22.14)') '  t0                   = ',db%t0
  write(iunit,'(a,i0)')      '  nstep                = ',db%nstep
  write(iunit,'(a,i0)')      '  nt_subsampled        = ',db%nt_subsampled
  write(iunit,'(a,i0)')      '  subsample_step       = ',db%subsample_step
  write(iunit,'(a,es22.14)') '  dt_subsampled        = ',db%dt * dble(db%subsample_step)
  write(iunit,'(a)') ''

  write(iunit,'(a)') 'mesh'
  write(iunit,'(a,i0,a,i0,a,i0)') '  ngll                 = ',db%ngllx,' x ',db%nglly,' x ',db%ngllz
  write(iunit,'(a,i0)')     '  neighbor_shells      = ',db%neighbor_shells
  write(iunit,'(a,i0)')     '  buffer_size          = ',db%buffer_size
  write(iunit,'(a,es22.14)') '  scale_displ          = ',db%scale_displ
  write(iunit,'(a,es22.14)') '  R_PLANET             = ',db%R_PLANET
  if (db%RHOAV > 0.d0) then
    write(iunit,'(a,es22.14)') '  RHOAV                = ',db%RHOAV
  else
    write(iunit,'(a)')      '  RHOAV                = (not stored; using the build default)'
  endif
  write(iunit,'(a)') ''

  write(iunit,'(a)') 'simulation flags'
  write(iunit,'(a,l1)')     '  TOPOGRAPHY           = ',db%topography
  write(iunit,'(a,l1)')     '  ELLIPTICITY          = ',db%ellipticity
  write(iunit,'(a,l1)')     '  ROTATION             = ',db%rotation
  write(iunit,'(a,l1)')     '  ATTENUATION          = ',db%attenuation
  write(iunit,'(a,l1)')     '  GRAVITY              = ',db%gravity
  write(iunit,'(a)') ''

  if (db%topography) then
    write(iunit,'(a)') 'topography'
    write(iunit,'(a,i0)')      '  NX_BATHY             = ',db%NX_BATHY
    write(iunit,'(a,i0)')      '  NY_BATHY             = ',db%NY_BATHY
    write(iunit,'(a,es22.14)') '  RESOLUTION_TOPO_FILE = ',db%RESOLUTION_TOPO_FILE
    if (with_topo) then
      call gf_load_topo(db,ierr)
      if (ierr /= GF_OK) then
        write(iunit,'(a,a)') '  grid                 = NOT LOADED: ',trim(gf_errmsg)
      else
        write(iunit,'(a,i0,a,i0)') '  ibathy_topo min/max  = ', &
          minval(db%ibathy_topo),' / ',maxval(db%ibathy_topo)
      endif
    endif
    write(iunit,'(a)') ''
  endif

  if (db%ellipticity) then
    write(iunit,'(a)') 'ellipticity'
    write(iunit,'(a,i0)')      '  nspl                 = ',db%nspl
    write(iunit,'(a,es22.14)') '  rspl(1)              = ',db%rspl(1)
    write(iunit,'(a,es22.14)') '  rspl(nspl)           = ',db%rspl(db%nspl)
    write(iunit,'(a)') ''
  endif

  write(iunit,'(a)') 'elements'
  write(iunit,'(a,i0)')  '  nelem                = ',db%nelem
  if (db%nelem > 0) then
    write(iunit,'(a,a,a,3(es22.14,1x))') '  first                = ',db%morton_hex(1), &
      '  ',db%centroid(1,1),db%centroid(2,1),db%centroid(3,1)
    write(iunit,'(a,a,a,3(es22.14,1x))') '  last                 = ',db%morton_hex(db%nelem), &
      '  ',db%centroid(1,db%nelem),db%centroid(2,db%nelem),db%centroid(3,db%nelem)
  endif
  write(iunit,'(a)') ''

  write(iunit,'(a)') 'stations'
  write(iunit,'(a,i0)')  '  nstations            = ',db%nstations
  do i = 1,db%nstations
    write(iunit,'(a)') ''
    write(iunit,'(a,a)')       '  station              = ',trim(db%stations(i)%id)
    write(iunit,'(a,es22.14)') '    latitude           = ',db%stations(i)%latitude
    write(iunit,'(a,es22.14)') '    longitude          = ',db%stations(i)%longitude
    write(iunit,'(a,es22.14)') '    depth              = ',db%stations(i)%depth
    write(iunit,'(a,es22.14)') '    hdur               = ',db%stations(i)%hdur
    write(iunit,'(a,es22.14)') '    f_cutoff           = ',db%stations(i)%f_cutoff
    write(iunit,'(a,es22.14)') '    factor_force_source= ',db%stations(i)%factor_force_source
    write(iunit,'(a,es22.14)') '    time_shift         = ',db%stations(i)%time_shift
    write(iunit,'(a,i0)')      '    stf samples        = ',size(db%stations(i)%stf)
    write(iunit,'(a,es22.14)') '    stf max            = ',maxval(abs(db%stations(i)%stf))

    ! probing the topography here is not decoration: get_topo_bathy() returns
    ! plausible garbage rather than an error when the shared_parameters
    ! bootstrap in gf_shared_params.F90 has the grid geometry wrong, and this
    ! is the cheapest place to notice
    if (db%topography .and. with_topo .and. db%topo_loaded) then
      call gf_topo_elevation(db,db%stations(i)%latitude,db%stations(i)%longitude,elevation)
      write(iunit,'(a,es22.14,a)') '    surface elevation  = ',elevation,' m'
    endif
  enddo

  end subroutine gf_print_info

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_topo_elevation(db,lat,lon,elevation)

! surface elevation at a geographic position, through the solver's own
! get_topo_bathy() (src/shared/model_topo_bathy.f90)

  implicit none

  type(t_gfdb), intent(in) :: db
  double precision, intent(in) :: lat,lon
  double precision, intent(out) :: elevation

  elevation = 0.d0
  if (.not. db%topo_loaded) return

  call get_topo_bathy(lat,lon,elevation,db%ibathy_topo)

  end subroutine gf_topo_elevation

!
!-------------------------------------------------------------------------------------------------
!
!  filesystem helpers
!
!  Fortran's inquire(file=...) is not required to report on directories and
!  compilers disagree, and there is no portable way to read a directory at
!  all — hence src/gf3d/gf_dirlist.c.
!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_file_exists(filename,exists)

  implicit none

  character(len=*), intent(in) :: filename
  logical, intent(out) :: exists

  ! local parameters
  integer :: n,want_dirs,iexists

  n = len_trim(filename)
  want_dirs = 0
  iexists = 0

  call gf_path_exists(trim(filename),n,want_dirs,iexists)

  exists = (iexists /= 0)

  end subroutine gf_file_exists

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_dir_exists(dirname,exists)

  implicit none

  character(len=*), intent(in) :: dirname
  logical, intent(out) :: exists

  ! local parameters
  integer :: n,want_dirs,iexists

  n = len_trim(dirname)
  want_dirs = 1
  iexists = 0

  call gf_path_exists(trim(dirname),n,want_dirs,iexists)

  exists = (iexists /= 0)

  end subroutine gf_dir_exists

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_list_dir(dirname,suffix,want_dirs,namelen,names,nfound,ierr)

! sorted listing of a directory
!
! Scans twice: once to count, once to fill, so the caller never has to guess
! a capacity.

  implicit none

  character(len=*), intent(in) :: dirname
  character(len=*), intent(in) :: suffix
  logical, intent(in) :: want_dirs
  integer, intent(in) :: namelen
  character(len=namelen), dimension(:), allocatable, intent(out) :: names
  integer, intent(out) :: nfound
  integer, intent(out) :: ierr

  ! local parameters
  ! note: `sfx` is a fixed-size buffer rather than trim(suffix) so that C is
  !       never handed the address of a zero-length string; `slen` alone says
  !       how much of it to look at
  character(len=32) :: sfx
  integer :: ndir,slen,nmax,ier,iwant

  nfound = 0
  ndir = len_trim(dirname)
  slen = len_trim(suffix)
  if (slen > len(sfx)) then
    call gf_set_error(ierr,GF_ERR_ARG,'suffix filter is too long')
    return
  endif
  sfx = suffix

  iwant = 0
  if (want_dirs) iwant = 1

  ! first pass: count only. Uses the same declared type for `names` as the
  ! second pass, so that there is one consistent implicit interface here.
  allocate(names(1),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate a directory listing')
    return
  endif
  names(:) = ''

  nmax = 0
  call gf_dir_scan(trim(dirname),ndir,sfx,slen,iwant,nmax,namelen,names,nfound,ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_IO,'could not list directory: '//trim(dirname))
    return
  endif

  if (nfound <= 0) then
    deallocate(names)
    allocate(names(0),stat=ier)
    if (ier /= 0) then
      call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate a directory listing')
      return
    endif
    ierr = GF_OK
    return
  endif

  deallocate(names)
  allocate(names(nfound),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'could not allocate a directory listing')
    return
  endif
  names(:) = ''

  nmax = nfound
  call gf_dir_scan(trim(dirname),ndir,sfx,slen,iwant,nmax,namelen,names,nfound,ier)
  if (ier == 2) then
    call gf_set_error(ierr,GF_ERR_FORMAT, &
      'an entry name in '//trim(dirname)//' is longer than this reader expects')
    return
  else if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_IO,'could not list directory: '//trim(dirname))
    return
  endif

  ierr = GF_OK

  end subroutine gf_list_dir

  end module gf_database
