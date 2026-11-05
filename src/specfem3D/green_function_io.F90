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
!---- Green function database: HDF5 file creation, buffered writing,
!---- and time loop integration.
!----
!---- Each MPI rank independently creates/writes one HDF5 file per
!---- tagged element: {GF_DATABASE_PATH}/elements/{morton_hex}/{net}.{sta}.h5
!----
!---- Dataset "displacement": shape (3, 3, 5, 5, 5, nt_sub) in Fortran order
!----   where nt_sub = NSTEP / GF_SUBSAMPLE_STEP
!----   dim 1: force component (N, E, Z) — GF_NCOMP_FORCE
!----   dim 2: displacement component (x, y, z) — GF_NCOMP_DISP
!----   dims 3-5: GLL indices (NGLLX, NGLLY, NGLLZ)
!----   dim 6: subsampled time snapshots
!----
!---- Each simulation run writes into one force component slice.
!---- Files are created on the first run and opened (not truncated)
!---- on subsequent runs, preserving previously written components.
!----
!---- Chunking: (1, 3, 5, 5, 5, min(GF_BUFFER_SIZE, nt_sub))
!----   one force component per chunk, aligned to the write pattern
!----

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_init_hdf5()

! Initializes HDF5 files and allocates the write buffer.
! Called from prepare_green_function_storage() during solver setup.
!
! On the first run (N component), creates new HDF5 files with fill value 0.
! On subsequent runs (E, Z), opens existing files without truncating.

#ifdef USE_HDF5
  use hdf5

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,IMAIN,MAX_STRING_LEN, &
                        GF_NCOMP_DISP,GF_NCOMP_FORCE

  use specfem_par, only: NSTEP,DT,myrank

  use shared_parameters, only: GF_DATABASE_PATH,GF_SUBSAMPLE_STEP,GF_BUFFER_SIZE,GF_OVERWRITE

  use green_function_par, only: &
    gf_nelem_local, gf_morton_hex, gf_center_xyz, &
    gf_buffer, gf_isnap, gf_ibuf, gf_nt_sub, &
    gf_network_name, gf_station_name, &
    gf_force_component, gf_f_cutoff, gf_hdur, &
    gf_elem_active

  implicit none

  ! local variables
  integer :: i, ier
  integer :: hdferr
  character(len=MAX_STRING_LEN) :: filepath
  logical :: file_exists
  integer :: nskipped

  ! HDF5 identifiers (local — each file is opened and closed immediately)
  integer(HID_T) :: fid, dspace_id, dset_id, dcpl_id
  integer(HID_T) :: aspace_id, attr_id

  ! dataset dimensions (6D)
  integer(HSIZE_T), dimension(6) :: dims, chunk_dims
  integer(HSIZE_T), dimension(1) :: adim

  ! attribute values
  double precision :: attr_dp(1)
  integer :: attr_int(1)

  ! fill value
  real :: fill_val_r
  double precision :: fill_val_d

  ! nothing to do if no local elements
  if (gf_nelem_local == 0) return

  ! compute number of subsampled timesteps
  gf_nt_sub = NSTEP / GF_SUBSAMPLE_STEP

  ! initialize counters
  gf_isnap = 0
  gf_ibuf = 0

  ! allocate write buffer: (3, 5, 5, 5, bufsize, nelem)
  ! element index is LAST so each element's time slices are contiguous
  allocate(gf_buffer(GF_NCOMP_DISP, NGLLX, NGLLY, NGLLZ, GF_BUFFER_SIZE, gf_nelem_local), &
           stat=ier)
  if (ier /= 0) call exit_MPI(myrank, 'Error allocating gf_buffer')
  gf_buffer(:,:,:,:,:,:) = 0.0_CUSTOM_REAL

  ! allocate element active mask (completion tracking)
  allocate(gf_elem_active(gf_nelem_local), stat=ier)
  if (ier /= 0) call exit_MPI(myrank, 'Error allocating gf_elem_active')
  gf_elem_active(:) = .true.

  ! initialize HDF5 Fortran interface
  call h5open_f(hdferr)
  if (hdferr /= 0) call exit_MPI(myrank, 'Error initializing HDF5')

  ! dataset dimensions: (GF_NCOMP_FORCE, GF_NCOMP_DISP, NGLLX, NGLLY, NGLLZ, nt_sub)
  dims(1) = GF_NCOMP_FORCE
  dims(2) = GF_NCOMP_DISP
  dims(3) = NGLLX
  dims(4) = NGLLY
  dims(5) = NGLLZ
  dims(6) = gf_nt_sub

  ! chunk dimensions: (1, 3, 5, 5, 5, min(buffer, nt_sub))
  ! one force component per chunk
  chunk_dims(1) = 1
  chunk_dims(2) = GF_NCOMP_DISP
  chunk_dims(3) = NGLLX
  chunk_dims(4) = NGLLY
  chunk_dims(5) = NGLLZ
  chunk_dims(6) = min(GF_BUFFER_SIZE, gf_nt_sub)

  ! create or open one HDF5 file per local element
  do i = 1, gf_nelem_local
    filepath = trim(GF_DATABASE_PATH) // '/elements/' // &
               gf_morton_hex(i) // '/' // &
               trim(gf_network_name) // '.' // trim(gf_station_name) // '.h5'

    ! check if file already exists (from a previous force component run)
    inquire(file=trim(filepath), exist=file_exists)

    if (file_exists) then
      ! open existing file — preserves data from other force components
      call h5fopen_f(trim(filepath), H5F_ACC_RDWR_F, fid, hdferr)
      if (hdferr /= 0) call exit_MPI(myrank, 'Error opening existing GF HDF5 file')

      ! check completion: if computed_ALL == 1 and not overwriting, skip this element
      if (.not. GF_OVERWRITE) then
        call h5aopen_f(fid, 'computed_ALL', attr_id, hdferr)
        if (hdferr == 0) then
          call h5aread_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
          call h5aclose_f(attr_id, hdferr)
          if (attr_int(1) == 1) then
            gf_elem_active(i) = .false.
          endif
        endif
      endif

      call h5fclose_f(fid, hdferr)
    else
      ! create new file with dataset
      call h5fcreate_f(trim(filepath), H5F_ACC_TRUNC_F, fid, hdferr)
      if (hdferr /= 0) call exit_MPI(myrank, 'Error creating GF HDF5 file')

      ! create dataspace (6D)
      call h5screate_simple_f(6, dims, dspace_id, hdferr)
      if (hdferr /= 0) call exit_MPI(myrank, 'Error creating GF HDF5 dataspace')

      ! create dataset creation property list with chunking and fill value
      call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, hdferr)
      if (hdferr /= 0) call exit_MPI(myrank, 'Error creating GF HDF5 property list')

      call h5pset_chunk_f(dcpl_id, 6, chunk_dims, hdferr)
      if (hdferr /= 0) call exit_MPI(myrank, 'Error setting GF HDF5 chunk size')

      ! set fill value to 0 so unwritten force components are zero
      if (CUSTOM_REAL == 4) then
        fill_val_r = 0.0
        call h5pset_fill_value_f(dcpl_id, H5T_NATIVE_REAL, fill_val_r, hdferr)
      else
        fill_val_d = 0.0d0
        call h5pset_fill_value_f(dcpl_id, H5T_NATIVE_DOUBLE, fill_val_d, hdferr)
      endif

      ! create dataset
      if (CUSTOM_REAL == 4) then
        call h5dcreate_f(fid, 'displacement', H5T_NATIVE_REAL, dspace_id, dset_id, hdferr, &
                          dcpl_id=dcpl_id)
      else
        call h5dcreate_f(fid, 'displacement', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr, &
                          dcpl_id=dcpl_id)
      endif
      if (hdferr /= 0) call exit_MPI(myrank, 'Error creating GF HDF5 dataset')

      ! write attributes on the dataset
      adim(1) = 1

      ! dt (double)
      call h5screate_simple_f(1, adim, aspace_id, hdferr)
      call h5acreate_f(dset_id, 'dt', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
      attr_dp(1) = DT
      call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
      call h5aclose_f(attr_id, hdferr)
      call h5sclose_f(aspace_id, hdferr)

      ! subsample_step (integer)
      call h5screate_simple_f(1, adim, aspace_id, hdferr)
      call h5acreate_f(dset_id, 'subsample_step', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
      attr_int(1) = GF_SUBSAMPLE_STEP
      call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
      call h5aclose_f(attr_id, hdferr)
      call h5sclose_f(aspace_id, hdferr)

      ! nsnap (integer)
      call h5screate_simple_f(1, adim, aspace_id, hdferr)
      call h5acreate_f(dset_id, 'nsnap', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
      attr_int(1) = gf_nt_sub
      call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
      call h5aclose_f(attr_id, hdferr)
      call h5sclose_f(aspace_id, hdferr)

      ! centroid coordinates cx, cy, cz (double)
      call h5screate_simple_f(1, adim, aspace_id, hdferr)
      call h5acreate_f(dset_id, 'cx', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
      attr_dp(1) = dble(gf_center_xyz(1, i))
      call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
      call h5aclose_f(attr_id, hdferr)
      call h5sclose_f(aspace_id, hdferr)

      call h5screate_simple_f(1, adim, aspace_id, hdferr)
      call h5acreate_f(dset_id, 'cy', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
      attr_dp(1) = dble(gf_center_xyz(2, i))
      call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
      call h5aclose_f(attr_id, hdferr)
      call h5sclose_f(aspace_id, hdferr)

      call h5screate_simple_f(1, adim, aspace_id, hdferr)
      call h5acreate_f(dset_id, 'cz', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
      attr_dp(1) = dble(gf_center_xyz(3, i))
      call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
      call h5aclose_f(attr_id, hdferr)
      call h5sclose_f(aspace_id, hdferr)

      ! hdur (double)
      call h5screate_simple_f(1, adim, aspace_id, hdferr)
      call h5acreate_f(dset_id, 'hdur', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
      attr_dp(1) = gf_hdur
      call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
      call h5aclose_f(attr_id, hdferr)
      call h5sclose_f(aspace_id, hdferr)

      ! f_cutoff (double)
      call h5screate_simple_f(1, adim, aspace_id, hdferr)
      call h5acreate_f(dset_id, 'f_cutoff', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
      attr_dp(1) = gf_f_cutoff
      call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
      call h5aclose_f(attr_id, hdferr)
      call h5sclose_f(aspace_id, hdferr)

      ! ngll (integer)
      call h5screate_simple_f(1, adim, aspace_id, hdferr)
      call h5acreate_f(dset_id, 'ngll', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
      attr_int(1) = NGLLX
      call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
      call h5aclose_f(attr_id, hdferr)
      call h5sclose_f(aspace_id, hdferr)

      ! close dataset, property list, dataspace
      call h5dclose_f(dset_id, hdferr)
      call h5pclose_f(dcpl_id, hdferr)
      call h5sclose_f(dspace_id, hdferr)

      ! completion tracking attributes (on the file root, all zero initially)
      attr_int(1) = 0

      call h5screate_simple_f(1, adim, aspace_id, hdferr)
      call h5acreate_f(fid, 'computed_N', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
      call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
      call h5aclose_f(attr_id, hdferr)
      call h5sclose_f(aspace_id, hdferr)

      call h5screate_simple_f(1, adim, aspace_id, hdferr)
      call h5acreate_f(fid, 'computed_E', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
      call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
      call h5aclose_f(attr_id, hdferr)
      call h5sclose_f(aspace_id, hdferr)

      call h5screate_simple_f(1, adim, aspace_id, hdferr)
      call h5acreate_f(fid, 'computed_Z', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
      call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
      call h5aclose_f(attr_id, hdferr)
      call h5sclose_f(aspace_id, hdferr)

      call h5screate_simple_f(1, adim, aspace_id, hdferr)
      call h5acreate_f(fid, 'computed_ALL', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
      call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
      call h5aclose_f(attr_id, hdferr)
      call h5sclose_f(aspace_id, hdferr)

      call h5fclose_f(fid, hdferr)
    endif
  enddo

  ! count skipped elements
  nskipped = 0
  do i = 1, gf_nelem_local
    if (.not. gf_elem_active(i)) nskipped = nskipped + 1
  enddo

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Green function database: HDF5 files initialized'
    write(IMAIN,*) '  force component:       ', gf_force_component
    write(IMAIN,*) '  subsampled timesteps:  ', gf_nt_sub
    write(IMAIN,*) '  buffer size:           ', GF_BUFFER_SIZE
    write(IMAIN,*) '  chunk time dimension:  ', min(GF_BUFFER_SIZE, gf_nt_sub)
    write(IMAIN,*) '  dataset shape:         (', GF_NCOMP_FORCE, ',', GF_NCOMP_DISP, ',', &
                   NGLLX, ',', NGLLY, ',', NGLLZ, ',', gf_nt_sub, ')'
    if (nskipped > 0) then
      write(IMAIN,*) '  elements skipped (already completed): ', nskipped
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

#else
  ! no HDF5 support compiled
  use shared_parameters, only: GF_DATABASE_ENABLED
  use specfem_par, only: myrank
  implicit none
  if (GF_DATABASE_ENABLED) then
    call exit_MPI(myrank, &
      'GF_DATABASE_ENABLED requires HDF5. Recompile with --with-hdf5')
  endif
#endif

  end subroutine gf_init_hdf5

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_write_timestep(it)

! Called from iterate_time() at every timestep.
! Checks if this timestep should be subsampled, extracts displacement
! into the buffer, and flushes when the buffer is full.

#ifdef USE_HDF5

  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ,GF_NCOMP_DISP

  use shared_parameters, only: GF_DATABASE_ENABLED,GF_SUBSAMPLE_STEP,GF_BUFFER_SIZE

  use specfem_par, only: scale_displ,GPU_MODE,Mesh_pointer,NGLOB_CRUST_MANTLE

  use specfem_par_crustmantle, only: &
    displ_crust_mantle, &
    ibool => ibool_crust_mantle, &
    displ => displ_crust_mantle

  use green_function_par, only: &
    gf_nelem_local, gf_local_elements, &
    gf_buffer, gf_isnap, gf_ibuf, gf_elem_active

  implicit none

  integer, intent(in) :: it

  ! local variables
  integer :: ielem, ispec, i, j, k, iglob
  real(kind=CUSTOM_REAL) :: scale

  ! skip if no local elements or GF not enabled
  if (.not. GF_DATABASE_ENABLED) return
  if (gf_nelem_local == 0) return

  ! check if this timestep is a subsampled output step
  if (mod(it, GF_SUBSAMPLE_STEP) /= 0) return

  ! In GPU mode the wavefield lives on the device during the time loop; copy
  ! the crust/mantle displacement back to the host before reading it.
  if (GPU_MODE) then
    call transfer_displ_cm_from_device(NDIM*NGLOB_CRUST_MANTLE, displ_crust_mantle, Mesh_pointer)
  endif

  ! increment buffer position
  gf_ibuf = gf_ibuf + 1
  gf_isnap = gf_isnap + 1

  ! scale_displ converts non-dimensional displacement to meters
  scale = real(scale_displ, kind=CUSTOM_REAL)

  ! extract displacement for all tagged elements into buffer
  do ielem = 1, gf_nelem_local
    if (.not. gf_elem_active(ielem)) cycle
    ispec = gf_local_elements(ielem)
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i, j, k, ispec)
          gf_buffer(1, i, j, k, gf_ibuf, ielem) = displ(1, iglob) * scale
          gf_buffer(2, i, j, k, gf_ibuf, ielem) = displ(2, iglob) * scale
          gf_buffer(3, i, j, k, gf_ibuf, ielem) = displ(3, iglob) * scale
        enddo
      enddo
    enddo
  enddo

  ! flush buffer if full
  if (gf_ibuf == GF_BUFFER_SIZE) then
    call gf_flush_buffer()
  endif

#else
  implicit none
  integer, intent(in) :: it
  ! unused
  integer :: idummy
  idummy = it
#endif

  end subroutine gf_write_timestep

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_flush_buffer()

! Writes the current buffer contents to HDF5 files via hyperslab writes.
! Each element's data is written to its own file.
!
! The buffer contains displacement for ONE force component (the current run).
! The hyperslab offset uses gf_force_component to write into the correct
! slice of the 6D dataset.

#ifdef USE_HDF5
  use hdf5
  use iso_c_binding, only: c_loc

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,MAX_STRING_LEN, &
                        GF_NCOMP_DISP

  use specfem_par, only: myrank

  use shared_parameters, only: GF_DATABASE_PATH

  use green_function_par, only: &
    gf_nelem_local, gf_morton_hex, &
    gf_buffer, gf_isnap, gf_ibuf, &
    gf_network_name, gf_station_name, &
    gf_force_component, gf_elem_active

  implicit none

  ! local variables
  integer :: ielem, hdferr
  character(len=MAX_STRING_LEN) :: filepath

  ! HDF5 identifiers
  integer(HID_T) :: fid, dset_id, fspace_id, mspace_id

  ! hyperslab selection (6D file dataset)
  integer(HSIZE_T), dimension(6) :: hs_offset, hs_count

  ! memory dataspace dimensions
  integer(HSIZE_T), dimension(6) :: mem_dims

  ! nothing to flush
  if (gf_ibuf == 0) return
  if (gf_nelem_local == 0) return

  ! memory dimensions matching the hyperslab count
  mem_dims(1) = 1              ! single force component
  mem_dims(2) = GF_NCOMP_DISP
  mem_dims(3) = NGLLX
  mem_dims(4) = NGLLY
  mem_dims(5) = NGLLZ
  mem_dims(6) = gf_ibuf

  ! hyperslab count in file dataset (6D)
  hs_count(1) = 1              ! one force component
  hs_count(2) = GF_NCOMP_DISP
  hs_count(3) = NGLLX
  hs_count(4) = NGLLY
  hs_count(5) = NGLLZ
  hs_count(6) = gf_ibuf

  ! hyperslab offset (0-based)
  hs_offset(1) = gf_force_component - 1   ! which force component (N = 0, E = 1, Z = 2)
  hs_offset(2) = 0
  hs_offset(3) = 0
  hs_offset(4) = 0
  hs_offset(5) = 0
  hs_offset(6) = gf_isnap - gf_ibuf       ! time offset

  do ielem = 1, gf_nelem_local
    if (.not. gf_elem_active(ielem)) cycle

    filepath = trim(GF_DATABASE_PATH) // '/elements/' // &
               gf_morton_hex(ielem) // '/' // &
               trim(gf_network_name) // '.' // trim(gf_station_name) // '.h5'

    ! open existing file
    call h5fopen_f(trim(filepath), H5F_ACC_RDWR_F, fid, hdferr)
    if (hdferr /= 0) call exit_MPI(myrank, 'Error opening GF HDF5 file for writing')

    ! open dataset
    call h5dopen_f(fid, 'displacement', dset_id, hdferr)
    if (hdferr /= 0) call exit_MPI(myrank, 'Error opening GF displacement dataset')

    ! get file dataspace (6D) and select hyperslab
    call h5dget_space_f(dset_id, fspace_id, hdferr)
    if (hdferr /= 0) call exit_MPI(myrank, 'Error getting GF file dataspace')

    call h5sselect_hyperslab_f(fspace_id, H5S_SELECT_SET_F, hs_offset, hs_count, hdferr)
    if (hdferr /= 0) call exit_MPI(myrank, 'Error selecting GF HDF5 hyperslab')

    ! create memory dataspace matching hyperslab count
    call h5screate_simple_f(6, mem_dims, mspace_id, hdferr)
    if (hdferr /= 0) call exit_MPI(myrank, 'Error creating GF memory dataspace')

    ! write data using c_loc to bypass Fortran rank checking
    ! gf_buffer(:,:,:,:, 1:ibuf, ielem) is contiguous in memory
    if (CUSTOM_REAL == 4) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, &
                       c_loc(gf_buffer(1, 1, 1, 1, 1, ielem)), hdferr, &
                       mem_space_id=mspace_id, file_space_id=fspace_id)
    else
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                       c_loc(gf_buffer(1, 1, 1, 1, 1, ielem)), hdferr, &
                       mem_space_id=mspace_id, file_space_id=fspace_id)
    endif
    if (hdferr /= 0) call exit_MPI(myrank, 'Error writing GF HDF5 data')

    ! close all
    call h5sclose_f(mspace_id, hdferr)
    call h5sclose_f(fspace_id, hdferr)
    call h5dclose_f(dset_id, hdferr)
    call h5fclose_f(fid, hdferr)
  enddo

  ! reset buffer position
  gf_ibuf = 0

#else
  implicit none
  ! stub — nothing to do without HDF5
#endif

  end subroutine gf_flush_buffer

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_finalize_hdf5()

! Flushes any remaining buffer contents and cleans up.
! Called from finalize_simulation().

#ifdef USE_HDF5
  use constants, only: IMAIN
  use specfem_par, only: myrank

  use shared_parameters, only: GF_DATABASE_ENABLED

  use green_function_par, only: &
    gf_nelem_local, gf_buffer, gf_ibuf, gf_isnap, gf_nt_sub, gf_elem_active

  implicit none

  if (.not. GF_DATABASE_ENABLED) return
  if (gf_nelem_local == 0) return

  ! flush remaining data in buffer
  if (gf_ibuf > 0) then
    call gf_flush_buffer()
  endif

  ! update completion flags on element and station files
  call gf_update_completion_flags()

  ! verify all snapshots were written
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Green function database: HDF5 finalized'
    write(IMAIN,*) '  total snapshots written: ', gf_isnap, ' / ', gf_nt_sub
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! deallocate buffer and active mask
  if (allocated(gf_buffer)) deallocate(gf_buffer)
  if (allocated(gf_elem_active)) deallocate(gf_elem_active)

#else
  implicit none
  ! stub — nothing to do without HDF5
#endif

  end subroutine gf_finalize_hdf5
