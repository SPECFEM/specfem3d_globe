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

  subroutine initialize_gindex()

  use gindex_par

  implicit none

  include 'version.fh'

  ! local parameters
  character(len = 20) :: snproc
  ! mpi
  integer :: sizeprocs
  integer :: ier

  ! runs in single process mode
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! checks number of processes
  ! note: must run as a single process with: mpirun -np 1 ..
  if (sizeprocs /= 1) then
    ! usage info
    if (myrank == 0) then
      print *, 'xgindex3D requires MPI functionality. However, this program executes as sequential program.'
      print *, 'Invalid number of processes used: ', sizeprocs, ' procs'
      print *
      print *, 'Please run: mpirun -np 1 ./bin/xgindex3D <number_of_solver_processes>'
      print *
      print *, '            for example: mpirun -np 1 ./bin/xgindex3D 96'
      print *
    endif
    call abort_mpi()
  endif

  ! reads input parameters
  if (command_argument_count() /= 1) then
    ! usage info
    print *, 'Usage: mpirun -np 1 ./bin/xgindex3D <number_of_solver_processes>'
    print *
    print *, '       for example: mpirun -np 1 ./bin/xgindex3D 96'
    print *
    stop 'Wrong number of arguments'
  endif

  call get_command_argument(1,snproc)
  read(snproc,*) nproc

  ! open main output file, only written to by process 0
  if (myrank == 0) then
    if (IMAIN /= ISTANDARD_OUTPUT) then
      open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_gindex3D.txt',status='unknown',action='write',iostat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening file output_gindex3D.txt for writing output info')
    endif

    write(IMAIN,*)
    write(IMAIN,*) '******************************'
    write(IMAIN,*) '**** Specfem3D gindex3D   ****'
    write(IMAIN,*) '******************************'
    write(IMAIN,*)
    write(IMAIN,*) 'Version: ', git_package_version
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  if (myrank == 0) write(IMAIN,'(a)') '<< xgindex3D...'

  !allocate(ndof_p2p(nproc,nproc))
  !ndof_p2p=0

  ! initializes simulation parameters
  if (myrank == 0) write(IMAIN,'(a)',advance='no') '   initialising...'

  ! reads in Par_file and sets compute parameters
  call read_compute_parameters()

  ! read the mesh parameters for all array setup
  call read_mesh_parameters()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'mesh parameters (from input directory):'
    write(IMAIN,*) '  NSPEC_CRUST_MANTLE = ',NSPEC_CRUST_MANTLE
    write(IMAIN,*) '  NSPEC_OUTER_CORE   = ',NSPEC_OUTER_CORE
    write(IMAIN,*) '  NSPEC_INNER_CORE   = ',NSPEC_INNER_CORE
    write(IMAIN,*)
    write(IMAIN,*) '  NSPEC_TRINFINITE   = ',NSPEC_TRINFINITE
    write(IMAIN,*) '  NSPEC_INFINITE     = ',NSPEC_INFINITE
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine initialize_gindex


