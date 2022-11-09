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


  subroutine prepare_noise()

  use specfem_par
  use specfem_par_crustmantle, only: NSPEC_TOP
  use specfem_par_noise

  implicit none
  ! local parameters
  integer :: ier
  double precision :: sizeval

  ! NOISE TOMOGRAPHY
  ! checks if anything to do
  if (NOISE_TOMOGRAPHY == 0) return

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) "preparing noise arrays"
    write(IMAIN,*) "  NOISE_TOMOGRAPHY = ",NOISE_TOMOGRAPHY
    call flush_IMAIN()
  endif

  ! checks noise setup
  call check_parameters_noise()

  ! determines file i/o buffer size for surface movies
  ! (needed for better performance on clusters, otherwise i/o will become a serious bottleneck)
  ! size of a single noise movie snapshot at surface (in MB)
  sizeval = dble(CUSTOM_REAL) * dble(NDIM) * dble(NGLLX) * dble(NGLLY) * dble(NSPEC_TOP) / 1024.d0 / 1024.d0
  ! sets file i/o buffer size
  if (NOISE_TOMOGRAPHY == 3 .and. UNDO_ATTENUATION) then
    ! needs to align with attenuation buffer size, otherwise things will get very complicated
    NT_DUMP_NOISE_BUFFER = NT_DUMP_ATTENUATION
  else
    ! sets a user specified maximum size (given in MB)
    NT_DUMP_NOISE_BUFFER = int(MAXIMUM_NOISE_BUFFER_SIZE_IN_MB / sizeval)
    ! limits size
    if (NT_DUMP_NOISE_BUFFER > NSTEP) NT_DUMP_NOISE_BUFFER = NSTEP
  endif

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) "  timing:"
    write(IMAIN,*) "    start time           = ",sngl(-t0)," seconds"
    write(IMAIN,*) "    time step            = ",sngl(DT)," s"
    write(IMAIN,*) "    number of time steps = ",NSTEP
    ! noise surface movie array size
    ! (holds displacement at surface for a single time step)
    write(IMAIN,*) "  arrays:"
    write(IMAIN,*) "    size of noise surface movie array = ",sngl(sizeval),"MB"
    write(IMAIN,*) "                                      = ",sngl(sizeval / 1024.d0),"GB"
    ! buffer size for file i/o
    write(IMAIN,*) "  noise buffer: "
    write(IMAIN,*) "    number of buffered time steps = ",NT_DUMP_NOISE_BUFFER
    write(IMAIN,*) "    size of noise buffer array for each slice = ",sngl(sizeval * dble(NT_DUMP_NOISE_BUFFER)),"MB"
    write(IMAIN,*) "                                              = ",sngl(sizeval * dble(NT_DUMP_NOISE_BUFFER) / 1024.d0),"GB"
    call flush_IMAIN()
  endif

  ! synchronizes processes
  call synchronize_all()

  ! for noise tomography, number of surface (movie) points needed for 'surface movie';
  ! surface output must NOT be coarse (have to be saved on all GLL points)
  ! number of points
  num_noise_surface_points = NGLLX * NGLLY * NSPEC_TOP

  ! in case arrays become too big, the following can lead to segmentation faults and crash.
  ! let's just be more verbose.

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  allocating noise arrays:"
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! allocates noise arrays
  if (myrank == 0) then
    write(IMAIN,*) "    source array"
    call flush_IMAIN()
  endif
  if (NOISE_TOMOGRAPHY == 1) then
    ! main noise source (only needed for 1. step)
    allocate(noise_sourcearray(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating noise source array')
  else
    ! dummy
    allocate(noise_sourcearray(1,1,1,1,1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating noise source array')
  endif
  ! initializes
  noise_sourcearray(:,:,:,:,:) = 0._CUSTOM_REAL

  ! ensemble surface noise
  if (myrank == 0) then
    write(IMAIN,*) "    ensemble surface arrays"
    call flush_IMAIN()
  endif
  call synchronize_all()
  allocate(normal_x_noise(num_noise_surface_points), &
           normal_y_noise(num_noise_surface_points), &
           normal_z_noise(num_noise_surface_points), &
           mask_noise(num_noise_surface_points), &
           noise_surface_movie(NDIM,NGLLX,NGLLY,NSPEC_TOP),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating noise arrays')

  ! initializes
  normal_x_noise(:)            = 0._CUSTOM_REAL
  normal_y_noise(:)            = 0._CUSTOM_REAL
  normal_z_noise(:)            = 0._CUSTOM_REAL
  mask_noise(:)                = 0._CUSTOM_REAL
  noise_surface_movie(:,:,:,:) = 0._CUSTOM_REAL

  ! file i/o buffer
  ! checks integer size limit: size of buffer must fit onto an 4-byte integer ( < 2 GB)
  if (NSPEC_TOP > 2147483646 / (CUSTOM_REAL * NGLLX * NGLLY * NDIM * NT_DUMP_NOISE_BUFFER)) then
    print *,'buffer of noise surface_movie needed exceeds integer 4-byte limit: ',dble(reclen_noise) * dble(NT_DUMP_NOISE_BUFFER)
    print *,'  ',CUSTOM_REAL, NDIM, NGLLX * NGLLY, NSPEC_TOP,NT_DUMP_NOISE_BUFFER
    print *,'bit size Fortran: ',bit_size(NSPEC_TOP)
    print *,'NT_DUMP_NOISE_BUFFER: ',NT_DUMP_NOISE_BUFFER
    print *,'Please reduce size of noise buffer for file i/o ...'
    call exit_MPI(myrank,"Error NT_DUMP_NOISE_BUFFER leads to buffer length exceeding integer limit (2 GB)")
  endif

  ! allocates buffer memory
  if (myrank == 0) then
    write(IMAIN,*) "    noise buffer array"
    call flush_IMAIN()
  endif
  call synchronize_all()
  allocate(noise_buffer(NDIM,NGLLX,NGLLY,NSPEC_TOP,NT_DUMP_NOISE_BUFFER),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating noise buffer array')

  ! initializes buffer and counters
  noise_buffer(:,:,:,:,:) = 0._CUSTOM_REAL
  icounter_noise_buffer = 0
  nstep_subset_noise_buffer = 0

  ! gets noise parameters and sets up arrays
  call read_parameters_noise()

  ! read noise distribution and direction
  call read_noise_distribution_direction()

  ! user output
  if (myrank == 0) then
    ! noise simulations ignore the CMTSOLUTIONS sources but employ a noise-spectrum source S_squared instead
    write(IMAIN,*)
    write(IMAIN,*) "  ignoring CMT sources"
    write(IMAIN,*)
    select case (NOISE_TOMOGRAPHY)
    case (1)
      write(IMAIN,*) "  noise source uses main record id = ",irec_main_noise
      write(IMAIN,*) "  noise main station: ",trim(network_name(irec_main_noise))//'.'//trim(station_name(irec_main_noise))
    case (2)
      write(IMAIN,*) "  noise source uses ensemble forward source"
    case (3)
      write(IMAIN,*) "  reconstructing ensemble forward wavefield"
      write(IMAIN,*) "  noise source uses ensemble adjoint sources"
    end select
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! user output of distances to main station
  if (NOISE_TOMOGRAPHY == 1) call print_main_distances_noise()

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_noise
