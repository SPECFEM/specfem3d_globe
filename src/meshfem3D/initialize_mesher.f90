!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
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

  subroutine initialize_mesher()

  use meshfem3D_par
  use meshfem3D_models_par

  use manager_adios

  implicit none

  include 'version.fh'

  ! local parameters
  integer :: sizeprocs
  ! timing
  double precision, external :: wtime

! sizeprocs returns number of processes started (should be equal to NPROCTOT).
! myrank is the rank of each process, between 0 and NPROCTOT-1.
! as usual in MPI, process 0 is in charge of coordinating everything
! and also takes care of the main output
! do not create anything for the inner core here, will be done in solver
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! set the base pathname for output files
  OUTPUT_FILES = OUTPUT_FILES_BASE

  ! open main output file, only written to by process 0
  if (myrank == 0) then
    if (IMAIN /= ISTANDARD_OUTPUT) &
      open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_mesher.txt',status='unknown')

    write(IMAIN,*)
    write(IMAIN,*) '****************************'
    write(IMAIN,*) '*** Specfem3D MPI Mesher ***'
    write(IMAIN,*) '****************************'
    write(IMAIN,*)
    write(IMAIN,*) 'Version: ', git_version
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! get MPI starting time
  time_start = wtime()

  if (myrank == 0) then
    ! reads the parameter file and computes additional parameters
    call read_compute_parameters()
  endif

  ! broadcast parameters read from main to all processes
  call broadcast_computed_parameters()

  ! check that the code is running with the requested number of processes
  if (sizeprocs /= NPROCTOT) then
    if (myrank == 0) print *,'Error wrong number of MPI processes ',sizeprocs,' should be ',NPROCTOT,', please check...'
    call exit_MPI(myrank,'wrong number of MPI processes')
  endif

  ! synchronizes processes
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    select case(PLANET_TYPE)
    case (IPLANET_EARTH)
      write(IMAIN,*) 'Planet: Earth'
    case (IPLANET_MARS)
      write(IMAIN,*) 'Planet: Mars'
    case (IPLANET_MOON)
      write(IMAIN,*) 'Natural satellite: Moon'
    case default
      call exit_MPI(myrank,'Invalid planet, type not recognized yet')
    end select
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! additional initialization on this system (ADIOS,OpenMP,..)
  call im_initialize_system()

  end subroutine initialize_mesher

!
!-------------------------------------------------------------------------------------------------
!

  subroutine im_initialize_system()

  use meshfem3D_par
  use meshfem3D_models_par

  use manager_adios

  implicit none

  ! compute rotation matrix from Euler angles
  ANGULAR_WIDTH_XI_RAD = ANGULAR_WIDTH_XI_IN_DEGREES * DEGREES_TO_RADIANS
  ANGULAR_WIDTH_ETA_RAD = ANGULAR_WIDTH_ETA_IN_DEGREES * DEGREES_TO_RADIANS

  if (NCHUNKS /= 6) call euler_angles(rotation_matrix,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)

  ! ADIOS
  if (ADIOS_ENABLED) then
    call initialize_adios()
  endif

  ! gravity integrals
  if (GRAVITY_INTEGRALS) then
    call gravity_initialize_integrals()
  endif

  ! OpenMP
  call init_openmp()

  ! synchronizes processes
  call synchronize_all()

  end subroutine im_initialize_system



