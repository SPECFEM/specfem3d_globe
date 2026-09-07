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
!---- Per-element reads: the geometry of one element, and the reciprocal
!---- displacement field of one (element, station) pair.
!----
!---- Layout, from the writer, which is the authoritative description
!---- (src/specfem3D/green_function_metadata.F90:150-162 and
!---- src/specfem3D/green_function_io.F90:35-49):
!----
!----   elements/{morton_hex}/coordinates.h5
!----       xyz(3,NGLLX,NGLLY,NGLLZ)                     Fortran order
!----       attributes cx, cy, cz, morton_hex
!----
!----   elements/{morton_hex}/{NET}.{STA}.h5
!----       displacement(3_force,3_disp,NGLLX,NGLLY,NGLLZ,nt_sub)
!----       dim 1: force component at the station, N=1 E=2 Z=3
!----              (green_function_io.F90:506, hs_offset(1) = gf_force_component - 1)
!----       dim 2: displacement component, Cartesian x/y/z
!----       dims 3-5: GLL indices, xi <-> i, eta <-> j, gamma <-> k
!----       dim 6: subsampled time
!----       attributes computed_N/E/Z/ALL on the file,
!----       dt, subsample_step, nsnap, hdur, f_cutoff, ngll, cx/cy/cz on the dataset
!----
!---- That index order is worth stating twice, because it is the single
!---- likeliest thing to get wrong in this library. Any loop structure
!---- carried over from GF3DF is suspect: it stored the same quantity as
!---- GF%displacement(ista,iforce,idisp,iglob,it), with the station index
!---- leading and the three GLL indices collapsed into one. Here the station
!---- is a separate *file* and (i,j,k) are explicit, so a literal
!---- transcription transposes. tests/gf3d/test_gf_interp.f90 exists to
!---- catch exactly that.
!----
!---- Memory: one element-station file is 21 MB in the shipped global
!---- example (3 x 3 x 125 x 4625 float32). Callers read one element at a
!---- time and loop stations inside; nothing here holds two.
!----

  module gf_element_io

  use gf_par, only: t_gfdb,gf_set_error, &
                    GF_OK,GF_ERR_ARG,GF_ERR_FORMAT,GF_NCOMP

  use gf_hdf5_read, only: GF_HID,gf_h5_file_open,gf_h5_file_close, &
                          gf_h5_dset_dims,gf_h5_dset_type_size, &
                          gf_h5_read_4d_d,gf_h5_read_6d_r

  implicit none

  private

  public :: gf_element_path
  public :: gf_read_element_coords
  public :: gf_read_element_displ
  public :: gf_element_type_sizes

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_element_path(db,ielem,ista,filename,ierr)

! builds the path of an element file
!
! `ista` of 0 asks for the element's coordinates.h5, otherwise for that
! station's displacement file. The Morton directory name is used untrimmed:
! the writer formats it as exactly 16 hexadecimal digits
! (green_function_morton.F90, Z16.16), so a trim() would silently corrupt a
! code whose leading digits happen to be blanks -- which cannot occur with
! Z16.16, but the invariant is worth not relying on twice.

  implicit none

  type(t_gfdb), intent(in) :: db
  integer, intent(in) :: ielem,ista
  character(len=*), intent(out) :: filename
  integer, intent(out) :: ierr

  filename = ''

  if (.not. db%is_open) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_element_path: database is not open')
    return
  endif

  if (ielem < 1 .or. ielem > db%nelem) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_element_path: element index out of range')
    return
  endif

  if (ista == 0) then
    filename = trim(db%path)//'/elements/'//db%morton_hex(ielem)//'/coordinates.h5'
  else
    if (ista < 1 .or. ista > db%nstations) then
      call gf_set_error(ierr,GF_ERR_ARG,'gf_element_path: station index out of range')
      return
    endif
    filename = trim(db%path)//'/elements/'//db%morton_hex(ielem)//'/' &
               //trim(db%stations(ista)%id)//'.h5'
  endif

  ierr = GF_OK

  end subroutine gf_element_path

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_read_element_coords(db,ielem,xyz_elem,ierr)

! reads the GLL coordinates of one element, widened to double precision

  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM,MAX_STRING_LEN

  implicit none

  type(t_gfdb), intent(in) :: db
  integer, intent(in) :: ielem
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ), intent(out) :: xyz_elem
  integer, intent(out) :: ierr

  ! local parameters
  character(len=MAX_STRING_LEN) :: filename
  integer(kind=GF_HID) :: fid
  integer(kind=8), dimension(4) :: dims
  integer :: ndims,ierr2

  xyz_elem(:,:,:,:) = 0.d0

  call gf_element_path(db,ielem,0,filename,ierr)
  if (ierr /= GF_OK) return

  call gf_h5_file_open(filename,fid,ierr)
  if (ierr /= GF_OK) return

  ndims = 4
  call gf_h5_dset_dims(fid,'xyz',ndims,dims,ierr)
  if (ierr /= GF_OK) goto 99

  if (ndims /= 4 .or. dims(1) /= NDIM .or. &
      dims(2) /= NGLLX .or. dims(3) /= NGLLY .or. dims(4) /= NGLLZ) then
    call gf_set_error(ierr,GF_ERR_FORMAT, &
      'unexpected shape for xyz, expected (3,NGLLX,NGLLY,NGLLZ), in '//trim(filename))
    goto 99
  endif

  call gf_h5_read_4d_d(fid,'xyz',NDIM,NGLLX,NGLLY,NGLLZ,xyz_elem,ierr)

99 continue
  call gf_h5_file_close(fid,ierr2)

  end subroutine gf_read_element_coords

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_read_element_displ(db,ielem,ista,displ,ierr)

! reads the full reciprocal displacement field of one (element, station)
!
! 21 MB in the shipped global example, so the caller owns the buffer and is
! expected to reuse it across stations rather than allocate per station.

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,MAX_STRING_LEN

  implicit none

  type(t_gfdb), intent(in) :: db
  integer, intent(in) :: ielem,ista
  real(kind=CUSTOM_REAL), dimension(GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ,db%nt_subsampled), &
    intent(out) :: displ
  integer, intent(out) :: ierr

  ! local parameters
  character(len=MAX_STRING_LEN) :: filename
  integer(kind=GF_HID) :: fid
  integer(kind=8), dimension(6) :: dims
  integer :: ndims,ierr2

  call gf_element_path(db,ielem,ista,filename,ierr)
  if (ierr /= GF_OK) return

  call gf_h5_file_open(filename,fid,ierr)
  if (ierr /= GF_OK) return

  ndims = 6
  call gf_h5_dset_dims(fid,'displacement',ndims,dims,ierr)
  if (ierr /= GF_OK) goto 99

  if (ndims /= 6 .or. dims(1) /= GF_NCOMP .or. dims(2) /= GF_NCOMP .or. &
      dims(3) /= NGLLX .or. dims(4) /= NGLLY .or. dims(5) /= NGLLZ .or. &
      dims(6) /= db%nt_subsampled) then
    call gf_set_error(ierr,GF_ERR_FORMAT, &
      'unexpected shape for displacement, expected (3,3,NGLLX,NGLLY,NGLLZ,nt_subsampled), in ' &
      //trim(filename))
    goto 99
  endif

  call gf_h5_read_6d_r(fid,'displacement',GF_NCOMP,GF_NCOMP,NGLLX,NGLLY,NGLLZ, &
                       db%nt_subsampled,displ,ierr)

99 continue
  call gf_h5_file_close(fid,ierr2)

  end subroutine gf_read_element_displ

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_element_type_sizes(db,tsize_xyz,tsize_displ,ierr)

! probes the on-disk element type sizes, in bytes, from the first element
!
! Reporting only, for `xgf3d --info`. The writer picks its HDF5 type from the
! CUSTOM_REAL the *solver* was built with, which need not match this build's
! -- see the note on gf_h5_dset_type_size. `tsize_displ` comes back as 0 when
! the database holds no stations.

  use constants, only: MAX_STRING_LEN

  implicit none

  type(t_gfdb), intent(in) :: db
  integer, intent(out) :: tsize_xyz,tsize_displ
  integer, intent(out) :: ierr

  ! local parameters
  character(len=MAX_STRING_LEN) :: filename
  integer(kind=GF_HID) :: fid
  integer :: ierr2

  tsize_xyz = 0
  tsize_displ = 0

  call gf_element_path(db,1,0,filename,ierr)
  if (ierr /= GF_OK) return

  call gf_h5_file_open(filename,fid,ierr)
  if (ierr /= GF_OK) return
  call gf_h5_dset_type_size(fid,'xyz',tsize_xyz,ierr)
  call gf_h5_file_close(fid,ierr2)
  if (ierr /= GF_OK) return

  if (db%nstations < 1) then
    ierr = GF_OK
    return
  endif

  call gf_element_path(db,1,1,filename,ierr)
  if (ierr /= GF_OK) return

  call gf_h5_file_open(filename,fid,ierr)
  if (ierr /= GF_OK) return
  call gf_h5_dset_type_size(fid,'displacement',tsize_displ,ierr)
  call gf_h5_file_close(fid,ierr2)

  end subroutine gf_element_type_sizes

  end module gf_element_io
