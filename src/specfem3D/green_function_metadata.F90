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
!---- Green function database: metadata, coordinate, and mesh info writing.
!----
!---- Per-element coordinates.h5: full GLL node coordinates (3,5,5,5)
!---- Per-database mesh_info.h5: topography + ellipticity arrays (write-once)
!---- Per-station stations/{net}.{sta}.h5: simulation params, STF, station info
!---- Completion tracking: computed_N/E/Z/ALL attributes on station displacement files
!----

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_write_coordinates()

! Writes per-element coordinate files:
!   {GF_DATABASE_PATH}/elements/{morton_hex}/coordinates.h5
!
! Contains the full GLL node coordinates (3, NGLLX, NGLLY, NGLLZ)
! needed for point-in-element location and Lagrange interpolation at extraction time.
!
! Idempotent: skips writing if the file already exists (another station run
! may have created it — coordinates are deterministic per element).

#ifdef USE_HDF5
  use hdf5

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,IMAIN,MAX_STRING_LEN

  use specfem_par, only: myrank

  use shared_parameters, only: GF_DATABASE_PATH

  use specfem_par_crustmantle, only: &
    ibool => ibool_crust_mantle, &
    xstore => xstore_crust_mantle, &
    ystore => ystore_crust_mantle, &
    zstore => zstore_crust_mantle

  use green_function_par, only: &
    gf_nelem_local, gf_local_elements, gf_morton_hex, gf_center_xyz

  implicit none

  integer :: ielem, ispec, i, j, k, iglob, hdferr
  character(len=MAX_STRING_LEN) :: filepath
  logical :: file_exists

  ! HDF5 identifiers
  integer(HID_T) :: fid, dspace_id, dset_id, aspace_id, attr_id
  integer(HSIZE_T), dimension(4) :: dims
  integer(HSIZE_T), dimension(1) :: adim

  ! coordinate buffer for one element
  real(kind=CUSTOM_REAL), dimension(3, NGLLX, NGLLY, NGLLZ), target :: xyz_elem

  ! attribute values
  double precision :: attr_dp(1)

  integer :: nwritten

  if (gf_nelem_local == 0) return

  ! initialize HDF5 Fortran interface (reference-counted, safe to call multiple times)
  call h5open_f(hdferr)
  if (hdferr /= 0) call exit_MPI(myrank, 'Error initializing HDF5 in gf_write_coordinates')

  dims(1) = 3
  dims(2) = NGLLX
  dims(3) = NGLLY
  dims(4) = NGLLZ

  adim(1) = 1

  nwritten = 0

  do ielem = 1, gf_nelem_local
    filepath = trim(GF_DATABASE_PATH) // '/elements/' // &
               gf_morton_hex(ielem) // '/coordinates.h5'

    ! skip if file already exists (idempotent)
    inquire(file=trim(filepath), exist=file_exists)
    if (file_exists) cycle

    ispec = gf_local_elements(ielem)

    ! extract GLL node coordinates
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i, j, k, ispec)
          xyz_elem(1, i, j, k) = xstore(iglob)
          xyz_elem(2, i, j, k) = ystore(iglob)
          xyz_elem(3, i, j, k) = zstore(iglob)
        enddo
      enddo
    enddo

    ! create HDF5 file
    !
    ! H5F_ACC_EXCL_F, not H5F_ACC_TRUNC_F, and a lost race is not an error.
    !
    ! coordinates.h5 depends only on the element, not on the station, so
    ! every station's reciprocal run writes the same set of files. When two
    ! of those runs are in flight at once (snakemake -j N over the stations)
    ! the inquire() above is a time-of-check/time-of-use window: both see the
    ! file missing, both create it, and the loser used to die inside HDF5
    ! with "unable to lock file, errno = 11" and abort the whole run.
    !
    ! H5F_ACC_EXCL_F maps to open(O_CREAT|O_EXCL), which is atomic, so
    ! exactly one process creates the file. The loser simply skips: the
    ! contents are a pure function of the element, so whoever won wrote the
    ! same bytes we would have.
    call h5fcreate_f(trim(filepath), H5F_ACC_EXCL_F, fid, hdferr)
    if (hdferr /= 0) then
      ! another process got there first -> nothing to do for this element
      inquire(file=trim(filepath), exist=file_exists)
      if (file_exists) cycle
      call exit_MPI(myrank, 'Error creating coordinates.h5')
    endif

    ! create dataset xyz (3, 5, 5, 5)
    call h5screate_simple_f(4, dims, dspace_id, hdferr)
    if (CUSTOM_REAL == 4) then
      call h5dcreate_f(fid, 'xyz', H5T_NATIVE_REAL, dspace_id, dset_id, hdferr)
    else
      call h5dcreate_f(fid, 'xyz', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
    endif
    if (hdferr /= 0) call exit_MPI(myrank, 'Error creating xyz dataset')

    ! write data
    if (CUSTOM_REAL == 4) then
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, xyz_elem, dims, hdferr)
    else
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, xyz_elem, dims, hdferr)
    endif
    if (hdferr /= 0) call exit_MPI(myrank, 'Error writing xyz dataset')

    call h5dclose_f(dset_id, hdferr)
    call h5sclose_f(dspace_id, hdferr)

    ! write centroid attributes for quick lookup
    call h5screate_simple_f(1, adim, aspace_id, hdferr)
    call h5acreate_f(fid, 'cx', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
    attr_dp(1) = dble(gf_center_xyz(1, ielem))
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
    call h5aclose_f(attr_id, hdferr)
    call h5sclose_f(aspace_id, hdferr)

    call h5screate_simple_f(1, adim, aspace_id, hdferr)
    call h5acreate_f(fid, 'cy', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
    attr_dp(1) = dble(gf_center_xyz(2, ielem))
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
    call h5aclose_f(attr_id, hdferr)
    call h5sclose_f(aspace_id, hdferr)

    call h5screate_simple_f(1, adim, aspace_id, hdferr)
    call h5acreate_f(fid, 'cz', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
    attr_dp(1) = dble(gf_center_xyz(3, ielem))
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
    call h5aclose_f(attr_id, hdferr)
    call h5sclose_f(aspace_id, hdferr)

    ! write morton_hex as attribute
    call h5screate_f(H5S_SCALAR_F, aspace_id, hdferr)
    call h5acreate_f(fid, 'morton_hex', H5T_NATIVE_CHARACTER, aspace_id, attr_id, hdferr)
    call h5aclose_f(attr_id, hdferr)
    call h5sclose_f(aspace_id, hdferr)

    call h5fclose_f(fid, hdferr)

    nwritten = nwritten + 1
  enddo

  if (myrank == 0) then
    write(IMAIN,*) 'Green function database: coordinate files written'
    write(IMAIN,*) '  new files on rank 0: ', nwritten
    write(IMAIN,*)
    call flush_IMAIN()
  endif

#else
  implicit none
  ! stub — nothing to do without HDF5
#endif

  end subroutine gf_write_coordinates

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_write_mesh_info()

! Writes the shared mesh info file:
!   {GF_DATABASE_PATH}/mesh_info.h5
!
! Contains topography (ibathy_topo) and ellipticity spline arrays needed
! for geographic-to-Cartesian coordinate conversion at extraction time.
!
! Written by rank 0 only. Idempotent: skips if file already exists.

#ifdef USE_HDF5
  use hdf5

  use constants, only: CUSTOM_REAL,NGLLX,IMAIN,MAX_STRING_LEN,NR_DENSITY

  use specfem_par, only: myrank, NSTEP, DT, t0, &
    ibathy_topo, NX_BATHY_VAL, NY_BATHY_VAL, &
    nspl_ellip, rspl_ellip, ellipicity_spline, ellipicity_spline2, &
    scale_displ

  use shared_parameters, only: GF_DATABASE_PATH, GF_SUBSAMPLE_STEP, GF_BUFFER_SIZE, &
    GF_NEIGHBOR_SHELLS, &
    TOPOGRAPHY, ELLIPTICITY, GRAVITY, ROTATION, ATTENUATION, &
    RESOLUTION_TOPO_FILE, R_PLANET

  use green_function_par, only: gf_nt_sub

  implicit none

  character(len=MAX_STRING_LEN) :: filepath
  logical :: file_exists
  integer :: hdferr

  ! HDF5 identifiers
  integer(HID_T) :: fid, dspace_id, dset_id, aspace_id, attr_id
  integer(HSIZE_T), dimension(2) :: dims2
  integer(HSIZE_T), dimension(1) :: dims1, adim

  ! attribute values
  double precision :: attr_dp(1)
  integer :: attr_int(1)

  ! validation
  double precision :: existing_dt(1)
  integer :: existing_nstep(1), existing_subsample(1)

  if (myrank /= 0) return

  ! initialize HDF5 Fortran interface (reference-counted, safe to call multiple times)
  call h5open_f(hdferr)
  if (hdferr /= 0) call exit_MPI(myrank, 'Error initializing HDF5 in gf_write_mesh_info')

  filepath = trim(GF_DATABASE_PATH) // '/mesh_info.h5'

  ! idempotent: skip if file already exists and is valid
  inquire(file=trim(filepath), exist=file_exists)
  if (file_exists) then
    ! validate by trying to open and read a known attribute
    call h5fopen_f(trim(filepath), H5F_ACC_RDONLY_F, fid, hdferr)

    ! The file exists but will not open. When several stations' reciprocal
    ! runs are in flight at once this is the normal case rather than an
    ! error: another run created mesh_info.h5 moments ago and still holds
    ! the HDF5 write lock on it, so our read fails with
    ! "unable to lock file, errno = 11". Falling through to the recreate
    ! path below would open it H5F_ACC_TRUNC_F, hit the same lock and abort
    ! the run.
    !
    ! Skipping is safe because mesh_info.h5 describes the mesh and the
    ! simulation parameters, not the station, so the run holding the lock is
    ! writing exactly what we would have written. A file left corrupt by a
    ! crashed run is not silently repaired here any more, but it is caught
    ! when the database is opened for extraction.
    if (hdferr /= 0) then
      write(IMAIN,*) 'Green function database: mesh_info.h5 is being written by a concurrent run, skipping'
      write(IMAIN,*)
      call flush_IMAIN()
      return
    endif

    if (hdferr == 0) then
      call h5aopen_f(fid, 'scale_displ', attr_id, hdferr)
      if (hdferr == 0) then
        call h5aclose_f(attr_id, hdferr)

        ! validate that simulation parameters match the current run
        call h5aopen_f(fid, 'dt', attr_id, hdferr)
        if (hdferr == 0) then
          call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, existing_dt, adim, hdferr)
          call h5aclose_f(attr_id, hdferr)
          if (abs(existing_dt(1) - DT) > 1.0d-12) then
            call h5fclose_f(fid, hdferr)
            write(IMAIN,*) 'ERROR: mesh_info.h5 has dt =', existing_dt(1), ' but current run has DT =', DT
            write(IMAIN,*) 'Delete the existing mesh_info.h5 and restart.'
            call flush_IMAIN()
            call exit_MPI(myrank, 'mesh_info.h5 dt mismatch — delete and restart')
          endif
        endif

        call h5aopen_f(fid, 'nstep', attr_id, hdferr)
        if (hdferr == 0) then
          call h5aread_f(attr_id, H5T_NATIVE_INTEGER, existing_nstep, adim, hdferr)
          call h5aclose_f(attr_id, hdferr)
          if (existing_nstep(1) /= NSTEP) then
            call h5fclose_f(fid, hdferr)
            write(IMAIN,*) 'ERROR: mesh_info.h5 has nstep =', existing_nstep(1), ' but current run has NSTEP =', NSTEP
            write(IMAIN,*) 'Delete the existing mesh_info.h5 and restart.'
            call flush_IMAIN()
            call exit_MPI(myrank, 'mesh_info.h5 nstep mismatch — delete and restart')
          endif
        endif

        call h5aopen_f(fid, 'subsample_step', attr_id, hdferr)
        if (hdferr == 0) then
          call h5aread_f(attr_id, H5T_NATIVE_INTEGER, existing_subsample, adim, hdferr)
          call h5aclose_f(attr_id, hdferr)
          if (existing_subsample(1) /= GF_SUBSAMPLE_STEP) then
            call h5fclose_f(fid, hdferr)
            write(IMAIN,*) 'ERROR: mesh_info.h5 has subsample_step =', existing_subsample(1), &
                           ' but current run has GF_SUBSAMPLE_STEP =', GF_SUBSAMPLE_STEP
            write(IMAIN,*) 'Delete the existing mesh_info.h5 and restart.'
            call flush_IMAIN()
            call exit_MPI(myrank, 'mesh_info.h5 subsample_step mismatch — delete and restart')
          endif
        endif

        call h5fclose_f(fid, hdferr)
        write(IMAIN,*) 'Green function database: mesh_info.h5 already exists, skipping'
        write(IMAIN,*)
        call flush_IMAIN()
        return
      endif
      call h5fclose_f(fid, hdferr)
    endif
    ! file exists but is corrupt/incomplete — remove and recreate
    write(IMAIN,*) 'Green function database: mesh_info.h5 is corrupt, recreating'
    call flush_IMAIN()
  endif

  ! create file
  if (file_exists) then
    ! recovery path: the file above was found to be corrupt or incomplete, so
    ! deliberately replace it
    call h5fcreate_f(trim(filepath), H5F_ACC_TRUNC_F, fid, hdferr)
    if (hdferr /= 0) call exit_MPI(myrank, 'Error creating mesh_info.h5')
  else
    ! mesh_info.h5 is shared by every station's reciprocal run, so when
    ! several of those run at once the inquire() above is a
    ! time-of-check/time-of-use window. Create exclusively and let the first
    ! writer win; the contents do not depend on the station. See the longer
    ! note in gf_write_coordinates().
    call h5fcreate_f(trim(filepath), H5F_ACC_EXCL_F, fid, hdferr)
    if (hdferr /= 0) then
      inquire(file=trim(filepath), exist=file_exists)
      if (file_exists) then
        write(IMAIN,*) 'Green function database: mesh_info.h5 written by a concurrent run, skipping'
        write(IMAIN,*)
        call flush_IMAIN()
        return
      endif
      call exit_MPI(myrank, 'Error creating mesh_info.h5')
    endif
  endif

  adim(1) = 1

  ! write scale_displ
  call h5screate_simple_f(1, adim, aspace_id, hdferr)
  call h5acreate_f(fid, 'scale_displ', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
  attr_dp(1) = scale_displ
  call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
  call h5aclose_f(attr_id, hdferr)
  call h5sclose_f(aspace_id, hdferr)

  ! write R_PLANET
  call h5screate_simple_f(1, adim, aspace_id, hdferr)
  call h5acreate_f(fid, 'R_PLANET', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
  attr_dp(1) = R_PLANET
  call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
  call h5aclose_f(attr_id, hdferr)
  call h5sclose_f(aspace_id, hdferr)

  ! topography flag
  call h5screate_simple_f(1, adim, aspace_id, hdferr)
  call h5acreate_f(fid, 'TOPOGRAPHY', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
  if (TOPOGRAPHY) then
    attr_int(1) = 1
  else
    attr_int(1) = 0
  endif
  call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
  call h5aclose_f(attr_id, hdferr)
  call h5sclose_f(aspace_id, hdferr)

  ! ellipticity flag
  call h5screate_simple_f(1, adim, aspace_id, hdferr)
  call h5acreate_f(fid, 'ELLIPTICITY', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
  if (ELLIPTICITY) then
    attr_int(1) = 1
  else
    attr_int(1) = 0
  endif
  call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
  call h5aclose_f(attr_id, hdferr)
  call h5sclose_f(aspace_id, hdferr)

  ! simulation parameters
  call h5screate_simple_f(1, adim, aspace_id, hdferr)
  call h5acreate_f(fid, 'dt', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
  attr_dp(1) = DT
  call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
  call h5aclose_f(attr_id, hdferr)
  call h5sclose_f(aspace_id, hdferr)

  ! t0: time of the first time step relative to the source origin.
  ! The simulation time axis is timeval(it) = (it-1)*DT - t0, so the first
  ! stored (subsampled) sample sits at -t0 + (first_snap-1)*DT. The GF
  ! reconstruction needs t0 to place its trace on an absolute time axis that
  ! is consistent with a forward simulation; storing it here makes the
  ! database self-contained (no need to parse output_solver.txt).
  call h5screate_simple_f(1, adim, aspace_id, hdferr)
  call h5acreate_f(fid, 't0', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
  attr_dp(1) = t0
  call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
  call h5aclose_f(attr_id, hdferr)
  call h5sclose_f(aspace_id, hdferr)

  call h5screate_simple_f(1, adim, aspace_id, hdferr)
  call h5acreate_f(fid, 'nstep', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
  attr_int(1) = NSTEP
  call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
  call h5aclose_f(attr_id, hdferr)
  call h5sclose_f(aspace_id, hdferr)

  call h5screate_simple_f(1, adim, aspace_id, hdferr)
  call h5acreate_f(fid, 'nt_subsampled', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
  attr_int(1) = NSTEP / GF_SUBSAMPLE_STEP
  call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
  call h5aclose_f(attr_id, hdferr)
  call h5sclose_f(aspace_id, hdferr)

  call h5screate_simple_f(1, adim, aspace_id, hdferr)
  call h5acreate_f(fid, 'subsample_step', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
  attr_int(1) = GF_SUBSAMPLE_STEP
  call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
  call h5aclose_f(attr_id, hdferr)
  call h5sclose_f(aspace_id, hdferr)

  call h5screate_simple_f(1, adim, aspace_id, hdferr)
  call h5acreate_f(fid, 'buffer_size', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
  attr_int(1) = GF_BUFFER_SIZE
  call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
  call h5aclose_f(attr_id, hdferr)
  call h5sclose_f(aspace_id, hdferr)

  call h5screate_simple_f(1, adim, aspace_id, hdferr)
  call h5acreate_f(fid, 'neighbor_shells', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
  attr_int(1) = GF_NEIGHBOR_SHELLS
  call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
  call h5aclose_f(attr_id, hdferr)
  call h5sclose_f(aspace_id, hdferr)

  call h5screate_simple_f(1, adim, aspace_id, hdferr)
  call h5acreate_f(fid, 'ngllx', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
  attr_int(1) = NGLLX
  call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
  call h5aclose_f(attr_id, hdferr)
  call h5sclose_f(aspace_id, hdferr)

  ! physics flags
  call h5screate_simple_f(1, adim, aspace_id, hdferr)
  call h5acreate_f(fid, 'rotation', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
  attr_int(1) = merge(1, 0, ROTATION)
  call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
  call h5aclose_f(attr_id, hdferr)
  call h5sclose_f(aspace_id, hdferr)

  call h5screate_simple_f(1, adim, aspace_id, hdferr)
  call h5acreate_f(fid, 'attenuation', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
  attr_int(1) = merge(1, 0, ATTENUATION)
  call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
  call h5aclose_f(attr_id, hdferr)
  call h5sclose_f(aspace_id, hdferr)

  call h5screate_simple_f(1, adim, aspace_id, hdferr)
  call h5acreate_f(fid, 'gravity', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
  attr_int(1) = merge(1, 0, GRAVITY)
  call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
  call h5aclose_f(attr_id, hdferr)
  call h5sclose_f(aspace_id, hdferr)

  ! topography data
  if (TOPOGRAPHY) then
    ! NX_BATHY, NY_BATHY attributes
    call h5screate_simple_f(1, adim, aspace_id, hdferr)
    call h5acreate_f(fid, 'NX_BATHY', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
    attr_int(1) = NX_BATHY_VAL
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
    call h5aclose_f(attr_id, hdferr)
    call h5sclose_f(aspace_id, hdferr)

    call h5screate_simple_f(1, adim, aspace_id, hdferr)
    call h5acreate_f(fid, 'NY_BATHY', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
    attr_int(1) = NY_BATHY_VAL
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
    call h5aclose_f(attr_id, hdferr)
    call h5sclose_f(aspace_id, hdferr)

    call h5screate_simple_f(1, adim, aspace_id, hdferr)
    call h5acreate_f(fid, 'RESOLUTION_TOPO_FILE', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
    attr_dp(1) = RESOLUTION_TOPO_FILE
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
    call h5aclose_f(attr_id, hdferr)
    call h5sclose_f(aspace_id, hdferr)

    ! ibathy_topo dataset (NX_BATHY, NY_BATHY)
    dims2(1) = NX_BATHY_VAL
    dims2(2) = NY_BATHY_VAL
    call h5screate_simple_f(2, dims2, dspace_id, hdferr)
    call h5dcreate_f(fid, 'ibathy_topo', H5T_NATIVE_INTEGER, dspace_id, dset_id, hdferr)
    if (hdferr /= 0) call exit_MPI(myrank, 'Error creating ibathy_topo dataset')
    call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, ibathy_topo, dims2, hdferr)
    if (hdferr /= 0) call exit_MPI(myrank, 'Error writing ibathy_topo dataset')
    call h5dclose_f(dset_id, hdferr)
    call h5sclose_f(dspace_id, hdferr)
  endif

  ! ellipticity data
  if (ELLIPTICITY) then
    ! nspl attribute
    call h5screate_simple_f(1, adim, aspace_id, hdferr)
    call h5acreate_f(fid, 'nspl', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
    attr_int(1) = nspl_ellip
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
    call h5aclose_f(attr_id, hdferr)
    call h5sclose_f(aspace_id, hdferr)

    ! rspl dataset
    dims1(1) = nspl_ellip
    call h5screate_simple_f(1, dims1, dspace_id, hdferr)
    call h5dcreate_f(fid, 'rspl', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, rspl_ellip(1:nspl_ellip), dims1, hdferr)
    call h5dclose_f(dset_id, hdferr)
    call h5sclose_f(dspace_id, hdferr)

    ! ellipicity_spline dataset
    dims1(1) = nspl_ellip
    call h5screate_simple_f(1, dims1, dspace_id, hdferr)
    call h5dcreate_f(fid, 'ellipicity_spline', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ellipicity_spline(1:nspl_ellip), dims1, hdferr)
    call h5dclose_f(dset_id, hdferr)
    call h5sclose_f(dspace_id, hdferr)

    ! ellipicity_spline2 dataset
    dims1(1) = nspl_ellip
    call h5screate_simple_f(1, dims1, dspace_id, hdferr)
    call h5dcreate_f(fid, 'ellipicity_spline2', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, ellipicity_spline2(1:nspl_ellip), dims1, hdferr)
    call h5dclose_f(dset_id, hdferr)
    call h5sclose_f(dspace_id, hdferr)
  endif

  call h5fclose_f(fid, hdferr)

  write(IMAIN,*) 'Green function database: mesh_info.h5 written'
  if (TOPOGRAPHY) write(IMAIN,*) '  ibathy_topo: ', NX_BATHY_VAL, ' x ', NY_BATHY_VAL
  if (ELLIPTICITY) write(IMAIN,*) '  ellipticity splines: nspl = ', nspl_ellip
  write(IMAIN,*) '  scale_displ = ', scale_displ
  write(IMAIN,*)
  call flush_IMAIN()

#else
  implicit none
  ! stub — nothing to do without HDF5
#endif

  end subroutine gf_write_mesh_info

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_write_station_metadata()

! Writes per-station metadata file:
!   {GF_DATABASE_PATH}/stations/{net}.{sta}.h5
!
! Contains station-specific info, STF array, and source parameters.
! Simulation-wide parameters are in mesh_info.h5.
! Written by rank 0. Creates the file on first component run,
! skips if it already exists (idempotent — station info is the same across N/E/Z runs).

#ifdef USE_HDF5
  use hdf5

  use constants, only: CUSTOM_REAL,IMAIN,MAX_STRING_LEN

  use specfem_par, only: myrank, NSTEP, &
    tshift_src, factor_force_source

  use shared_parameters, only: GF_DATABASE_PATH

  use green_function_par, only: &
    gf_network_name, gf_station_name, &
    gf_hdur, gf_f_cutoff, gf_stf, &
    gf_station_lat, gf_station_lon, gf_station_depth

  implicit none

  character(len=MAX_STRING_LEN) :: filepath, dirpath
  character(len=MAX_STRING_LEN) :: command
  logical :: file_exists
  integer :: hdferr

  ! HDF5 identifiers
  integer(HID_T) :: fid, dspace_id, dset_id, aspace_id, attr_id
  integer(HSIZE_T), dimension(1) :: dims1, adim

  ! attribute values
  double precision :: attr_dp(1)

  if (myrank /= 0) return

  ! initialize HDF5 Fortran interface (reference-counted, safe to call multiple times)
  call h5open_f(hdferr)
  if (hdferr /= 0) call exit_MPI(myrank, 'Error initializing HDF5 in gf_write_station_metadata')

  ! create stations directory
  dirpath = trim(GF_DATABASE_PATH) // '/stations'
  ! note: uses the system_command() wrapper, system() is a non-standard GNU extension
  command = 'mkdir -p ' // trim(dirpath)
  call system_command(command)

  filepath = trim(dirpath) // '/' // &
             trim(gf_network_name) // '.' // trim(gf_station_name) // '.h5'

  adim(1) = 1

  ! idempotent: skip if file already exists (station info is identical across N/E/Z runs)
  inquire(file=trim(filepath), exist=file_exists)

  if (.not. file_exists) then
    ! create new file with station-specific metadata
    call h5fcreate_f(trim(filepath), H5F_ACC_TRUNC_F, fid, hdferr)
    if (hdferr /= 0) call exit_MPI(myrank, 'Error creating station metadata file')

    ! station info
    call h5screate_simple_f(1, adim, aspace_id, hdferr)
    call h5acreate_f(fid, 'latitude', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
    attr_dp(1) = gf_station_lat
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
    call h5aclose_f(attr_id, hdferr)
    call h5sclose_f(aspace_id, hdferr)

    call h5screate_simple_f(1, adim, aspace_id, hdferr)
    call h5acreate_f(fid, 'longitude', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
    attr_dp(1) = gf_station_lon
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
    call h5aclose_f(attr_id, hdferr)
    call h5sclose_f(aspace_id, hdferr)

    call h5screate_simple_f(1, adim, aspace_id, hdferr)
    call h5acreate_f(fid, 'depth', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
    attr_dp(1) = gf_station_depth
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
    call h5aclose_f(attr_id, hdferr)
    call h5sclose_f(aspace_id, hdferr)

    ! STF parameters
    call h5screate_simple_f(1, adim, aspace_id, hdferr)
    call h5acreate_f(fid, 'hdur', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
    attr_dp(1) = gf_hdur
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
    call h5aclose_f(attr_id, hdferr)
    call h5sclose_f(aspace_id, hdferr)

    call h5screate_simple_f(1, adim, aspace_id, hdferr)
    call h5acreate_f(fid, 'f_cutoff', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
    attr_dp(1) = gf_f_cutoff
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
    call h5aclose_f(attr_id, hdferr)
    call h5sclose_f(aspace_id, hdferr)

    call h5screate_simple_f(1, adim, aspace_id, hdferr)
    call h5acreate_f(fid, 'factor_force_source', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
    attr_dp(1) = factor_force_source(1)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
    call h5aclose_f(attr_id, hdferr)
    call h5sclose_f(aspace_id, hdferr)

    call h5screate_simple_f(1, adim, aspace_id, hdferr)
    call h5acreate_f(fid, 'time_shift', H5T_NATIVE_DOUBLE, aspace_id, attr_id, hdferr)
    attr_dp(1) = tshift_src(1)
    call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, attr_dp, adim, hdferr)
    call h5aclose_f(attr_id, hdferr)
    call h5sclose_f(aspace_id, hdferr)

    ! STF array
    dims1(1) = NSTEP
    call h5screate_simple_f(1, dims1, dspace_id, hdferr)
    if (CUSTOM_REAL == 4) then
      call h5dcreate_f(fid, 'stf', H5T_NATIVE_REAL, dspace_id, dset_id, hdferr)
      call h5dwrite_f(dset_id, H5T_NATIVE_REAL, gf_stf, dims1, hdferr)
    else
      call h5dcreate_f(fid, 'stf', H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, gf_stf, dims1, hdferr)
    endif
    call h5dclose_f(dset_id, hdferr)
    call h5sclose_f(dspace_id, hdferr)

    call h5fclose_f(fid, hdferr)

    write(IMAIN,*) 'Green function database: station metadata created'
    write(IMAIN,*) '  file: ', trim(filepath)
    write(IMAIN,*)
    call flush_IMAIN()
  else
    write(IMAIN,*) 'Green function database: station metadata file exists, will update at finalization'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

#else
  implicit none
  ! stub — nothing to do without HDF5
#endif

  end subroutine gf_write_station_metadata

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_update_completion_flags()

! Updates completion flags on both the per-station metadata file and
! the per-element displacement files after successful data writing.
!
! Sets computed_N, computed_E, or computed_Z = 1 for the current force component.
! Called from gf_finalize_hdf5() after all data has been flushed.
! Rank 0 updates the station metadata file; each rank updates its own element files.

#ifdef USE_HDF5
  use hdf5

  use constants, only: IMAIN,MAX_STRING_LEN

  use specfem_par, only: myrank

  use shared_parameters, only: GF_DATABASE_PATH

  use green_function_par, only: &
    gf_nelem_local, gf_morton_hex, &
    gf_network_name, gf_station_name, gf_force_component

  implicit none

  integer :: ielem, hdferr
  character(len=MAX_STRING_LEN) :: filepath
  character(len=16) :: comp_attr_name
  integer :: attr_int(1), comp_N, comp_E, comp_Z
  integer(HSIZE_T), dimension(1) :: adim

  ! HDF5 identifiers
  integer(HID_T) :: fid, attr_id, aspace_id

  adim(1) = 1

  ! determine attribute name for current component
  select case (gf_force_component)
    case (1)
      comp_attr_name = 'computed_N'
    case (2)
      comp_attr_name = 'computed_E'
    case (3)
      comp_attr_name = 'computed_Z'
  end select

  ! update per-element station files
  do ielem = 1, gf_nelem_local
    filepath = trim(GF_DATABASE_PATH) // '/elements/' // &
               gf_morton_hex(ielem) // '/' // &
               trim(gf_network_name) // '.' // trim(gf_station_name) // '.h5'

    call h5fopen_f(trim(filepath), H5F_ACC_RDWR_F, fid, hdferr)
    if (hdferr /= 0) cycle

    ! delete old attribute and write new value = 1
    call h5adelete_f(fid, trim(comp_attr_name), hdferr)
    call h5screate_simple_f(1, adim, aspace_id, hdferr)
    call h5acreate_f(fid, trim(comp_attr_name), H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
    attr_int(1) = 1
    call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
    call h5aclose_f(attr_id, hdferr)
    call h5sclose_f(aspace_id, hdferr)

    ! check if all components are done
    comp_N = 0; comp_E = 0; comp_Z = 0

    call h5aopen_f(fid, 'computed_N', attr_id, hdferr)
    if (hdferr == 0) then
      call h5aread_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
      comp_N = attr_int(1)
      call h5aclose_f(attr_id, hdferr)
    endif

    call h5aopen_f(fid, 'computed_E', attr_id, hdferr)
    if (hdferr == 0) then
      call h5aread_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
      comp_E = attr_int(1)
      call h5aclose_f(attr_id, hdferr)
    endif

    call h5aopen_f(fid, 'computed_Z', attr_id, hdferr)
    if (hdferr == 0) then
      call h5aread_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
      comp_Z = attr_int(1)
      call h5aclose_f(attr_id, hdferr)
    endif

    ! set computed_ALL if all three are done
    if (comp_N == 1 .and. comp_E == 1 .and. comp_Z == 1) then
      ! delete existing computed_ALL if present
      call h5adelete_f(fid, 'computed_ALL', hdferr)
      call h5screate_simple_f(1, adim, aspace_id, hdferr)
      call h5acreate_f(fid, 'computed_ALL', H5T_NATIVE_INTEGER, aspace_id, attr_id, hdferr)
      attr_int(1) = 1
      call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, attr_int, adim, hdferr)
      call h5aclose_f(attr_id, hdferr)
      call h5sclose_f(aspace_id, hdferr)
    endif

    call h5fclose_f(fid, hdferr)
  enddo

  if (myrank == 0) then
    write(IMAIN,*) 'Green function database: completion flags updated'
    write(IMAIN,*) '  component: ', trim(comp_attr_name), ' = 1'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

#else
  implicit none
  ! stub — nothing to do without HDF5
#endif

  end subroutine gf_update_completion_flags
