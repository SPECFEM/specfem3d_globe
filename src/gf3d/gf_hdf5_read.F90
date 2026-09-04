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
!---- Thin serial HDF5 readers for the Green function database.
!----
!---- Deliberately *not* built on src/shared/hdf5_manager.F90: that module
!---- is collective and MPI-coupled (world_get_comm, h5_set_mpi_info), and
!---- this library must link without an MPI runtime. The raw-HDF5 template
!---- followed here is src/specfem3D/green_function_io.F90, which is also
!---- the authoritative description of the on-disk layout.
!----
!---- Precision: every floating-point read lands in a `double precision`
!---- buffer and asks HDF5 for H5T_NATIVE_DOUBLE, so a CUSTOM_REAL = 4
!---- database is widened by the HDF5 conversion layer on the way in and
!---- nothing downstream ever sees a float32. Integer reads use
!---- H5T_NATIVE_INTEGER the same way.
!----

  module gf_hdf5_read

#ifdef USE_HDF5
  use hdf5
#endif

  use gf_par, only: gf_set_error,GF_OK,GF_ERR_NO_HDF5,GF_ERR_NO_FILE,GF_ERR_HDF5,GF_ERR_FORMAT

  implicit none

  private

#ifdef USE_HDF5
  ! HDF5 object identifier kind, re-exported so that callers of this module
  ! need not `use hdf5` themselves
  integer, parameter, public :: GF_HID = HID_T
#else
  integer, parameter, public :: GF_HID = 8
#endif

  public :: gf_h5_init
  public :: gf_h5_file_open
  public :: gf_h5_file_close
  public :: gf_h5_has_attr
  public :: gf_h5_has_dset
  public :: gf_h5_read_attr_i
  public :: gf_h5_read_attr_d
  public :: gf_h5_dset_dims
  public :: gf_h5_dset_type_size
  public :: gf_h5_read_1d_d
  public :: gf_h5_read_2d_i
  public :: gf_h5_read_4d_d
  public :: gf_h5_read_6d_r

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_h5_init(ierr)

! initializes the HDF5 Fortran interface
!
! h5open_f() is reference counted, so calling this more than once is harmless

  implicit none

  integer, intent(out) :: ierr

#ifdef USE_HDF5
  integer :: hdferr

  call h5open_f(hdferr)
  if (hdferr /= 0) then
    call gf_set_error(ierr,GF_ERR_HDF5,'could not initialize the HDF5 Fortran interface')
    return
  endif

  ! Silences HDF5's automatic error reporting.
  !
  ! By default a failed HDF5 call prints a thirty-line diagnostic stack to
  ! stderr before returning its status. That is reasonable for an
  ! application and wrong for a library: this one reports failures through
  ! `ierr` and gf_errmsg so that the caller decides what the user sees, and
  ! several of its code paths fail *on purpose* — probing for an optional
  ! file, or trying one route into the element index before another.
  ! Leaving this on would make every such probe look like a disaster.
  call h5eset_auto_f(0,hdferr)
  if (hdferr /= 0) then
    call gf_set_error(ierr,GF_ERR_HDF5,'could not disable HDF5 automatic error reporting')
    return
  endif

  ierr = GF_OK
#else
  call gf_set_error(ierr,GF_ERR_NO_HDF5, &
                    'this build has no HDF5 support; re-run configure with --with-hdf5')
#endif

  end subroutine gf_h5_init

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_h5_file_open(filename,fid,ierr)

! opens an existing file read-only

  implicit none

  character(len=*), intent(in) :: filename
  integer(kind=GF_HID), intent(out) :: fid
  integer, intent(out) :: ierr

#ifdef USE_HDF5
  integer :: hdferr
  logical :: exists

  fid = -1

  ! distinguishes "no such file", which several callers probe for on
  ! purpose, from "this file is not readable HDF5", which is a real fault
  inquire(file=trim(filename),exist=exists)
  if (.not. exists) then
    call gf_set_error(ierr,GF_ERR_NO_FILE,'no such file: '//trim(filename))
    return
  endif

  call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, fid, hdferr)
  if (hdferr /= 0) then
    call gf_set_error(ierr,GF_ERR_HDF5,'could not open HDF5 file: '//trim(filename))
    return
  endif
  ierr = GF_OK
#else
  fid = -1
  call gf_set_error(ierr,GF_ERR_NO_HDF5, &
                    'this build has no HDF5 support; re-run configure with --with-hdf5')
#endif

  end subroutine gf_h5_file_open

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_h5_file_close(fid,ierr)

  implicit none

  integer(kind=GF_HID), intent(inout) :: fid
  integer, intent(out) :: ierr

#ifdef USE_HDF5
  integer :: hdferr

  ierr = GF_OK
  if (fid < 0) return

  call h5fclose_f(fid, hdferr)
  if (hdferr /= 0) call gf_set_error(ierr,GF_ERR_HDF5,'could not close HDF5 file')
  fid = -1
#else
  fid = -1
  call gf_set_error(ierr,GF_ERR_NO_HDF5,'this build has no HDF5 support')
#endif

  end subroutine gf_h5_file_close

!
!-------------------------------------------------------------------------------------------------
!

  logical function gf_h5_has_attr(loc_id,name)

! tests for the presence of an attribute
!
! used to keep databases written before an attribute was added readable

  implicit none

  integer(kind=GF_HID), intent(in) :: loc_id
  character(len=*), intent(in) :: name

#ifdef USE_HDF5
  integer :: hdferr
  logical :: exists

  exists = .false.
  call h5aexists_f(loc_id, trim(name), exists, hdferr)
  if (hdferr /= 0) exists = .false.

  gf_h5_has_attr = exists
#else
  gf_h5_has_attr = .false.
#endif

  end function gf_h5_has_attr

!
!-------------------------------------------------------------------------------------------------
!

  logical function gf_h5_has_dset(loc_id,name)

! tests for the presence of a dataset

  implicit none

  integer(kind=GF_HID), intent(in) :: loc_id
  character(len=*), intent(in) :: name

#ifdef USE_HDF5
  integer :: hdferr
  logical :: exists

  exists = .false.
  call h5lexists_f(loc_id, trim(name), exists, hdferr)
  if (hdferr /= 0) exists = .false.

  gf_h5_has_dset = exists
#else
  gf_h5_has_dset = .false.
#endif

  end function gf_h5_has_dset

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_h5_read_attr_i(loc_id,name,val,ierr)

! reads a scalar integer attribute

  implicit none

  integer(kind=GF_HID), intent(in) :: loc_id
  character(len=*), intent(in) :: name
  integer, intent(out) :: val
  integer, intent(out) :: ierr

#ifdef USE_HDF5
  integer(kind=GF_HID) :: attr_id
  integer(HSIZE_T), dimension(1) :: adim
  integer, dimension(1) :: buf
  integer :: hdferr

  val = 0
  adim(1) = 1

  call h5aopen_f(loc_id, trim(name), attr_id, hdferr)
  if (hdferr /= 0) then
    call gf_set_error(ierr,GF_ERR_FORMAT,'missing integer attribute: '//trim(name))
    return
  endif

  call h5aread_f(attr_id, H5T_NATIVE_INTEGER, buf, adim, hdferr)
  if (hdferr /= 0) then
    call h5aclose_f(attr_id, hdferr)
    call gf_set_error(ierr,GF_ERR_HDF5,'could not read integer attribute: '//trim(name))
    return
  endif

  call h5aclose_f(attr_id, hdferr)

  val = buf(1)
  ierr = GF_OK
#else
  val = 0
  call gf_set_error(ierr,GF_ERR_NO_HDF5,'this build has no HDF5 support')
#endif

  end subroutine gf_h5_read_attr_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_h5_read_attr_d(loc_id,name,val,ierr)

! reads a scalar attribute as double precision
!
! asking for H5T_NATIVE_DOUBLE lets HDF5 widen a float32 attribute for us

  implicit none

  integer(kind=GF_HID), intent(in) :: loc_id
  character(len=*), intent(in) :: name
  double precision, intent(out) :: val
  integer, intent(out) :: ierr

#ifdef USE_HDF5
  integer(kind=GF_HID) :: attr_id
  integer(HSIZE_T), dimension(1) :: adim
  double precision, dimension(1) :: buf
  integer :: hdferr

  val = 0.d0
  adim(1) = 1

  call h5aopen_f(loc_id, trim(name), attr_id, hdferr)
  if (hdferr /= 0) then
    call gf_set_error(ierr,GF_ERR_FORMAT,'missing attribute: '//trim(name))
    return
  endif

  call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, buf, adim, hdferr)
  if (hdferr /= 0) then
    call h5aclose_f(attr_id, hdferr)
    call gf_set_error(ierr,GF_ERR_HDF5,'could not read attribute: '//trim(name))
    return
  endif

  call h5aclose_f(attr_id, hdferr)

  val = buf(1)
  ierr = GF_OK
#else
  val = 0.d0
  call gf_set_error(ierr,GF_ERR_NO_HDF5,'this build has no HDF5 support')
#endif

  end subroutine gf_h5_read_attr_d

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_h5_dset_dims(loc_id,name,ndims,dims,ierr)

! returns the shape of a dataset
!
! `ndims` is the maximum rank the caller can accept on input and the actual
! rank on output; `dims` is in Fortran order

  implicit none

  integer(kind=GF_HID), intent(in) :: loc_id
  character(len=*), intent(in) :: name
  integer, intent(inout) :: ndims
  integer(kind=8), dimension(:), intent(out) :: dims
  integer, intent(out) :: ierr

#ifdef USE_HDF5
  integer(kind=GF_HID) :: dset_id,dspace_id
  integer(HSIZE_T), dimension(7) :: d,dmax
  integer :: hdferr,rank,i,nmax

  nmax = ndims
  ndims = 0
  dims(:) = 0

  call h5dopen_f(loc_id, trim(name), dset_id, hdferr)
  if (hdferr /= 0) then
    call gf_set_error(ierr,GF_ERR_FORMAT,'missing dataset: '//trim(name))
    return
  endif

  call h5dget_space_f(dset_id, dspace_id, hdferr)
  if (hdferr /= 0) then
    call h5dclose_f(dset_id, hdferr)
    call gf_set_error(ierr,GF_ERR_HDF5,'could not get dataspace of: '//trim(name))
    return
  endif

  call h5sget_simple_extent_ndims_f(dspace_id, rank, hdferr)
  if (hdferr /= 0 .or. rank < 0 .or. rank > 7) then
    call h5sclose_f(dspace_id, hdferr)
    call h5dclose_f(dset_id, hdferr)
    call gf_set_error(ierr,GF_ERR_HDF5,'could not get rank of: '//trim(name))
    return
  endif

  call h5sget_simple_extent_dims_f(dspace_id, d, dmax, hdferr)
  ! note: h5sget_simple_extent_dims_f returns the rank, not 0, on success
  if (hdferr /= rank) then
    call h5sclose_f(dspace_id, hdferr)
    call h5dclose_f(dset_id, hdferr)
    call gf_set_error(ierr,GF_ERR_HDF5,'could not get dimensions of: '//trim(name))
    return
  endif

  call h5sclose_f(dspace_id, hdferr)
  call h5dclose_f(dset_id, hdferr)

  if (rank > nmax .or. rank > size(dims)) then
    call gf_set_error(ierr,GF_ERR_FORMAT,'unexpected rank for dataset: '//trim(name))
    return
  endif

  ! note: the HDF5 Fortran wrappers transpose for us. h5screate_simple_f()
  !       reverses the Fortran dimensions on the way out (which is why h5py
  !       reports mesh_info.h5's ibathy_topo as (NY_BATHY,NX_BATHY) while the
  !       writer created it as (NX_BATHY,NY_BATHY)), and
  !       h5sget_simple_extent_dims_f() reverses them back on the way in. So
  !       `d` is already in Fortran order and must not be reversed again.
  do i = 1,rank
    dims(i) = int(d(i),kind=8)
  enddo
  ndims = rank

  ierr = GF_OK
#else
  ndims = 0
  dims(:) = 0
  call gf_set_error(ierr,GF_ERR_NO_HDF5,'this build has no HDF5 support')
#endif

  end subroutine gf_h5_dset_dims

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_h5_read_1d_d(loc_id,name,n,arr,ierr)

! reads a rank-1 dataset of n elements as double precision
!
! a float32 dataset (CUSTOM_REAL = 4, the usual case) is widened by HDF5

  implicit none

  integer(kind=GF_HID), intent(in) :: loc_id
  character(len=*), intent(in) :: name
  integer, intent(in) :: n
  double precision, dimension(n), intent(out) :: arr
  integer, intent(out) :: ierr

#ifdef USE_HDF5
  integer(kind=GF_HID) :: dset_id
  integer(HSIZE_T), dimension(1) :: dims
  integer :: hdferr

  arr(:) = 0.d0
  dims(1) = n

  call h5dopen_f(loc_id, trim(name), dset_id, hdferr)
  if (hdferr /= 0) then
    call gf_set_error(ierr,GF_ERR_FORMAT,'missing dataset: '//trim(name))
    return
  endif

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, arr, dims, hdferr)
  if (hdferr /= 0) then
    call h5dclose_f(dset_id, hdferr)
    call gf_set_error(ierr,GF_ERR_HDF5,'could not read dataset: '//trim(name))
    return
  endif

  call h5dclose_f(dset_id, hdferr)

  ierr = GF_OK
#else
  arr(:) = 0.d0
  call gf_set_error(ierr,GF_ERR_NO_HDF5,'this build has no HDF5 support')
#endif

  end subroutine gf_h5_read_1d_d

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_h5_read_2d_i(loc_id,name,n1,n2,arr,ierr)

! reads a rank-2 integer dataset, in Fortran order

  implicit none

  integer(kind=GF_HID), intent(in) :: loc_id
  character(len=*), intent(in) :: name
  integer, intent(in) :: n1,n2
  integer, dimension(n1,n2), intent(out) :: arr
  integer, intent(out) :: ierr

#ifdef USE_HDF5
  integer(kind=GF_HID) :: dset_id
  integer(HSIZE_T), dimension(2) :: dims
  integer :: hdferr

  arr(:,:) = 0
  dims(1) = n1
  dims(2) = n2

  call h5dopen_f(loc_id, trim(name), dset_id, hdferr)
  if (hdferr /= 0) then
    call gf_set_error(ierr,GF_ERR_FORMAT,'missing dataset: '//trim(name))
    return
  endif

  call h5dread_f(dset_id, H5T_NATIVE_INTEGER, arr, dims, hdferr)
  if (hdferr /= 0) then
    call h5dclose_f(dset_id, hdferr)
    call gf_set_error(ierr,GF_ERR_HDF5,'could not read dataset: '//trim(name))
    return
  endif

  call h5dclose_f(dset_id, hdferr)

  ierr = GF_OK
#else
  arr(:,:) = 0
  call gf_set_error(ierr,GF_ERR_NO_HDF5,'this build has no HDF5 support')
#endif

  end subroutine gf_h5_read_2d_i

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_h5_dset_type_size(loc_id,name,tsize,ierr)

! returns the storage size in bytes of a dataset's element type
!
! Used only for reporting. The writer emits H5T_NATIVE_REAL or
! H5T_NATIVE_DOUBLE according to the CUSTOM_REAL the *solver* was built with
! (green_function_metadata.F90:151-154, green_function_io.F90:201-204), which
! need not be the CUSTOM_REAL this reader was built with. Every read below
! goes through the HDF5 conversion layer and is therefore correct either
! way, but a double database read by a single-precision build silently
! narrows -- so `xgf3d --info` says what is actually on disk.

  implicit none

  integer(kind=GF_HID), intent(in) :: loc_id
  character(len=*), intent(in) :: name
  integer, intent(out) :: tsize
  integer, intent(out) :: ierr

#ifdef USE_HDF5
  integer(kind=GF_HID) :: dset_id,type_id
  integer(SIZE_T) :: sz
  integer :: hdferr

  tsize = 0

  call h5dopen_f(loc_id, trim(name), dset_id, hdferr)
  if (hdferr /= 0) then
    call gf_set_error(ierr,GF_ERR_FORMAT,'missing dataset: '//trim(name))
    return
  endif

  call h5dget_type_f(dset_id, type_id, hdferr)
  if (hdferr /= 0) then
    call h5dclose_f(dset_id, hdferr)
    call gf_set_error(ierr,GF_ERR_HDF5,'could not get datatype of: '//trim(name))
    return
  endif

  call h5tget_size_f(type_id, sz, hdferr)
  if (hdferr /= 0) then
    call h5tclose_f(type_id, hdferr)
    call h5dclose_f(dset_id, hdferr)
    call gf_set_error(ierr,GF_ERR_HDF5,'could not get datatype size of: '//trim(name))
    return
  endif

  call h5tclose_f(type_id, hdferr)
  call h5dclose_f(dset_id, hdferr)

  tsize = int(sz)

  ierr = GF_OK
#else
  tsize = 0
  call gf_set_error(ierr,GF_ERR_NO_HDF5,'this build has no HDF5 support')
#endif

  end subroutine gf_h5_dset_type_size

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_h5_read_4d_d(loc_id,name,n1,n2,n3,n4,arr,ierr)

! reads a rank-4 dataset as double precision, in Fortran order
!
! This is the element coordinate array, xyz(3,NGLLX,NGLLY,NGLLZ). It is
! always widened to double on the way in: it is 3 kB per element, the
! precision policy asks for it, and it means one code path covers both a
! CUSTOM_REAL = 4 and a CUSTOM_REAL = 8 database.
!
! Note that widening does not *recover* anything. With CUSTOM_REAL = 4 the
! solver had already rounded these coordinates to float32 before the writer
! saw them, which is why the 27-anchor consistency guard in gf_locate tests
! against GF_ANCHOR_TOL rather than against round-off. See the comment on
! GF_ANCHOR_TOL in gf_par.F90.

  implicit none

  integer(kind=GF_HID), intent(in) :: loc_id
  character(len=*), intent(in) :: name
  integer, intent(in) :: n1,n2,n3,n4
  double precision, dimension(n1,n2,n3,n4), intent(out) :: arr
  integer, intent(out) :: ierr

#ifdef USE_HDF5
  integer(kind=GF_HID) :: dset_id
  integer(HSIZE_T), dimension(4) :: dims
  integer :: hdferr

  arr(:,:,:,:) = 0.d0
  dims(1) = n1
  dims(2) = n2
  dims(3) = n3
  dims(4) = n4

  call h5dopen_f(loc_id, trim(name), dset_id, hdferr)
  if (hdferr /= 0) then
    call gf_set_error(ierr,GF_ERR_FORMAT,'missing dataset: '//trim(name))
    return
  endif

  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, arr, dims, hdferr)
  if (hdferr /= 0) then
    call h5dclose_f(dset_id, hdferr)
    call gf_set_error(ierr,GF_ERR_HDF5,'could not read dataset: '//trim(name))
    return
  endif

  call h5dclose_f(dset_id, hdferr)

  ierr = GF_OK
#else
  arr(:,:,:,:) = 0.d0
  call gf_set_error(ierr,GF_ERR_NO_HDF5,'this build has no HDF5 support')
#endif

  end subroutine gf_h5_read_4d_d

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_h5_read_6d_r(loc_id,name,n1,n2,n3,n4,n5,n6,arr,ierr)

! reads a rank-6 dataset into a real(CUSTOM_REAL) buffer, in Fortran order
!
! This is the bulk array: displacement(3_force,3_disp,NGLLX,NGLLY,NGLLZ,nt),
! 21 MB per element-station file in the shipped global example.
!
! It is the one place the library does *not* widen on read. Asking HDF5 for
! H5T_NATIVE_DOUBLE here would double the resident footprint of the single
! largest allocation in the whole extraction, for no gain: the interpolator
! in gf_interp.F90 widens one 225-element time slice at a time, which is
! where the double-precision core actually begins.
!
! The native type requested matches the *buffer*, not the file, so HDF5
! converts if the writer used the other CUSTOM_REAL. gf_h5_dset_type_size()
! lets a caller report when that conversion narrows.

  use constants, only: CUSTOM_REAL,SIZE_REAL

  implicit none

  integer(kind=GF_HID), intent(in) :: loc_id
  character(len=*), intent(in) :: name
  integer, intent(in) :: n1,n2,n3,n4,n5,n6
  real(kind=CUSTOM_REAL), dimension(n1,n2,n3,n4,n5,n6), intent(out) :: arr
  integer, intent(out) :: ierr

#ifdef USE_HDF5
  integer(kind=GF_HID) :: dset_id
  integer(HSIZE_T), dimension(6) :: dims
  integer :: hdferr

  arr(:,:,:,:,:,:) = 0._CUSTOM_REAL
  dims(1) = n1
  dims(2) = n2
  dims(3) = n3
  dims(4) = n4
  dims(5) = n5
  dims(6) = n6

  call h5dopen_f(loc_id, trim(name), dset_id, hdferr)
  if (hdferr /= 0) then
    call gf_set_error(ierr,GF_ERR_FORMAT,'missing dataset: '//trim(name))
    return
  endif

  if (CUSTOM_REAL == SIZE_REAL) then
    call h5dread_f(dset_id, H5T_NATIVE_REAL, arr, dims, hdferr)
  else
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, arr, dims, hdferr)
  endif
  if (hdferr /= 0) then
    call h5dclose_f(dset_id, hdferr)
    call gf_set_error(ierr,GF_ERR_HDF5,'could not read dataset: '//trim(name))
    return
  endif

  call h5dclose_f(dset_id, hdferr)

  ierr = GF_OK
#else
  arr(:,:,:,:,:,:) = 0._CUSTOM_REAL
  call gf_set_error(ierr,GF_ERR_NO_HDF5,'this build has no HDF5 support')
#endif

  end subroutine gf_h5_read_6d_r

  end module gf_hdf5_read
