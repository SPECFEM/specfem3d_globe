!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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
! the Free Software Foundation; either version 2 of the License, or
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

! create file OUTPUT_FILES/values_from_mesher.h based upon DATA/Par_file
! in order to compile the solver with the right array sizes

  program xcreate_header_file

  use shared_parameters
  use constants

  implicit none

  double precision :: static_memory_size

  integer :: NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
         NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUATION, &
         NSPEC_INNER_CORE_ATTENUATION, &
         NSPEC_CRUST_MANTLE_STR_OR_ATT,NSPEC_INNER_CORE_STR_OR_ATT, &
         NSPEC_CRUST_MANTLE_STR_AND_ATT,NSPEC_INNER_CORE_STR_AND_ATT, &
         NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_INNER_CORE_STRAIN_ONLY, &
         NSPEC_CRUST_MANTLE_ADJOINT, &
         NSPEC_OUTER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
         NGLOB_CRUST_MANTLE_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
         NGLOB_INNER_CORE_ADJOINT,NSPEC_OUTER_CORE_ROT_ADJOINT, &
         NSPEC_CRUST_MANTLE_STACEY,NSPEC_OUTER_CORE_STACEY, &
         NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION

!! DK DK for UNDO_ATTENUATION
  integer :: saved_SIMULATION_TYPE
  integer :: NT_DUMP_ATTENUATION_optimal,number_of_dumpings_to_do
  double precision :: static_memory_size_GB,size_to_store_at_each_time_step,disk_size_of_each_dumping

! ************** PROGRAM STARTS HERE **************

  print *
  print *,'creating file OUTPUT_FILES/values_from_mesher.h to compile solver with correct values'

! read the parameter file and compute additional parameters
  call read_compute_parameters()

! evaluate the amount of static memory needed by the solver
  call memory_eval(doubling_index,this_region_has_a_doubling, &
                   ner,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                   ratio_sampling_array,NPROCTOT, &
                   NSPEC,NGLOB, &
                   NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
                   NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUATION, &
                   NSPEC_INNER_CORE_ATTENUATION, &
                   NSPEC_CRUST_MANTLE_STR_OR_ATT,NSPEC_INNER_CORE_STR_OR_ATT, &
                   NSPEC_CRUST_MANTLE_STR_AND_ATT,NSPEC_INNER_CORE_STR_AND_ATT, &
                   NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_INNER_CORE_STRAIN_ONLY, &
                   NSPEC_CRUST_MANTLE_ADJOINT, &
                   NSPEC_OUTER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
                   NGLOB_CRUST_MANTLE_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
                   NGLOB_INNER_CORE_ADJOINT,NSPEC_OUTER_CORE_ROT_ADJOINT, &
                   NSPEC_CRUST_MANTLE_STACEY,NSPEC_OUTER_CORE_STACEY, &
                   NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION, &
                   NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                   static_memory_size)

  print *
  print *,'edit file OUTPUT_FILES/values_from_mesher.h to see'
  print *,'some statistics about the mesh'
  print *

  print *,'number of processors = ',NPROCTOT
  print *
  print *,'maximum number of points per region = ',nglob(IREGION_CRUST_MANTLE)
  print *
  print *,'total elements per slice = ',sum(NSPEC)
  print *,'total points per slice = ',sum(nglob)
  print *
  print *,'time-stepping of the solver will be: ',DT
  print *
  if(MOVIE_SURFACE .or. MOVIE_VOLUME) then
    print *,'MOVIE_VOLUME :',MOVIE_VOLUME
    print *,'MOVIE_SURFACE:',MOVIE_SURFACE
    print *,'Saving movie frames every',NTSTEP_BETWEEN_FRAMES
    print *
  endif
  print *,'on NEC SX, make sure "loopcnt=" parameter'
! use fused loops on NEC SX
  print *,'in Makefile is greater than max vector length = ',nglob(IREGION_CRUST_MANTLE)*NDIM
  print *

  print *,'approximate static memory needed by the solver:'
  print *,'----------------------------------------------'
  print *
  print *,'(lower bound, usually the real amount used is 5% to 10% higher)'
  print *
  print *,'(you can get a more precise estimate of the size used per MPI process'
  print *,' by typing "size -d bin/xspecfem3D"'
  print *,' after compiling the code with the DATA/Par_file you plan to use)'
  print *
  print *,'size of static arrays per slice = ',static_memory_size/1.d6,' MB'
  print *,'                                = ',static_memory_size/1048576.d0,' MiB'
  print *,'                                = ',static_memory_size/1.d9,' GB'
  print *,'                                = ',static_memory_size/1073741824.d0,' GiB'
  print *

! note: using less memory becomes an issue only if the strong scaling of the code is poor.
!          Some users will run simulations with an executable using far less than 80% RAM per core
!          if they prefer having a faster computational time (and use a higher number of cores).

  print *,'   (should be below 80% or 90% of the memory installed per core)'
  print *,'   (if significantly more, the job will not run by lack of memory)'
  print *,'   (note that if significantly less, you waste a significant amount'
  print *,'    of memory per processor core)'
  print *,'   (but that can be perfectly acceptable if you can afford it and'
  print *,'    want faster results by using more cores)'
  print *
  if(static_memory_size*dble(NPROCTOT)/1.d6 < 10000.d0) then
    print *,'size of static arrays for all slices = ',static_memory_size*dble(NPROCTOT)/1.d6,' MB'
    print *,'                                     = ',static_memory_size*dble(NPROCTOT)/1048576.d0,' MiB'
    print *,'                                     = ',static_memory_size*dble(NPROCTOT)/1.d9,' GB'
  else
    print *,'size of static arrays for all slices = ',static_memory_size*dble(NPROCTOT)/1.d9,' GB'
  endif
  print *,'                                     = ',static_memory_size*dble(NPROCTOT)/1073741824.d0,' GiB'
  print *,'                                     = ',static_memory_size*dble(NPROCTOT)/1.d12,' TB'
  print *,'                                     = ',static_memory_size*dble(NPROCTOT)/1099511627776.d0,' TiB'
  print *

!! DK DK for UNDO_ATTENUATION

  print *,'*******************************************************************************'
  print *,'Estimating optimal disk dumping interval for UNDO_ATTENUATION:'
  print *,'*******************************************************************************'
  print *

! optimal dumping interval calculation can only be done when SIMULATION_TYPE == 3 in the Par_file,
! thus set it to that value here in this serial code even if it has a different value in the Par_file
  saved_SIMULATION_TYPE = SIMULATION_TYPE
  SIMULATION_TYPE = 3

! evaluate the amount of static memory needed by the solver, but imposing that SIMULATION_TYPE = 3
! because that is by far the most expensive setup for runs in terms of memory usage, thus that is
! the type of run for which we need to make sure that everything fits in memory
  call memory_eval(doubling_index,this_region_has_a_doubling, &
                   ner,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                   ratio_sampling_array,NPROCTOT, &
                   NSPEC,NGLOB, &
                   NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
                   NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUATION, &
                   NSPEC_INNER_CORE_ATTENUATION, &
                   NSPEC_CRUST_MANTLE_STR_OR_ATT,NSPEC_INNER_CORE_STR_OR_ATT, &
                   NSPEC_CRUST_MANTLE_STR_AND_ATT,NSPEC_INNER_CORE_STR_AND_ATT, &
                   NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_INNER_CORE_STRAIN_ONLY, &
                   NSPEC_CRUST_MANTLE_ADJOINT, &
                   NSPEC_OUTER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
                   NGLOB_CRUST_MANTLE_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
                   NGLOB_INNER_CORE_ADJOINT,NSPEC_OUTER_CORE_ROT_ADJOINT, &
                   NSPEC_CRUST_MANTLE_STACEY,NSPEC_OUTER_CORE_STACEY, &
                   NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION, &
                   NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                   static_memory_size)

  call compute_optimized_dumping(static_memory_size,NT_DUMP_ATTENUATION_optimal,number_of_dumpings_to_do, &
                 static_memory_size_GB,size_to_store_at_each_time_step,disk_size_of_each_dumping)

! restore the simulation type that we have temporarily erased
  SIMULATION_TYPE = saved_SIMULATION_TYPE

  print *,'without undoing of attenuation you are using ',static_memory_size_GB,' GB per core'
  print *,'  i.e. ',sngl(100.d0 * static_memory_size_GB / MEMORY_INSTALLED_PER_CORE_IN_GB),'% of the installed memory'

  print *
  print *,'each time step to store in memory to undo attenuation will require storing ', &
                size_to_store_at_each_time_step,' GB per core'
  print *
  print *,'*******************************************************************************'
  print *,'the optimal value is thus NT_DUMP_ATTENUATION = ',NT_DUMP_ATTENUATION_optimal
  print *,'*******************************************************************************'

  print *
  print *,'we will need to save a total of ',number_of_dumpings_to_do,' dumpings (restart files) to disk'

  print *
  print *,'each dumping on the disk to undo attenuation will require storing ',disk_size_of_each_dumping,' GB per core'

  print *
  print *,'each dumping on the disk will require storing ',disk_size_of_each_dumping*NPROCTOT,' GB for all cores'

  print *
  print *,'ALL dumpings on the disk will require storing ',disk_size_of_each_dumping*number_of_dumpings_to_do,' GB per core'

  print *
  print *,'*******************************************************************************'
  print *,'ALL dumpings on the disk will require storing ', &
                 disk_size_of_each_dumping*number_of_dumpings_to_do*NPROCTOT,' GB for all cores'
  print *,'  i.e. ',disk_size_of_each_dumping*number_of_dumpings_to_do*NPROCTOT/1000.d0,' TB'
  print *,'*******************************************************************************'
  print *

  if(.not. UNDO_ATTENUATION) then
    print *
    print *,'*******************************************************************************'
    print *,'BEWARE, UNDO_ATTENUATION is .false. and thus exact undoing of attenuation is currently'
    print *,'turned off, i.e. the above related estimates are currently UNUSED.'
    print *,'*******************************************************************************'
    print *

    ! set to a dummy very large value
    NT_DUMP_ATTENUATION_optimal = 100000000
  endif

! create include file for the solver
  call save_header_file(NSPEC,NGLOB,NPROC,NPROCTOT, &
                        static_memory_size, &
                        NSPEC2D_TOP,NSPEC2D_BOTTOM, &
                        NSPEC2DMAX_YMIN_YMAX,NSPEC2DMAX_XMIN_XMAX, &
                        NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
                        NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUATION, &
                        NSPEC_INNER_CORE_ATTENUATION, &
                        NSPEC_CRUST_MANTLE_STR_OR_ATT,NSPEC_INNER_CORE_STR_OR_ATT, &
                        NSPEC_CRUST_MANTLE_STR_AND_ATT,NSPEC_INNER_CORE_STR_AND_ATT, &
                        NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_INNER_CORE_STRAIN_ONLY, &
                        NSPEC_CRUST_MANTLE_ADJOINT, &
                        NSPEC_OUTER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
                        NGLOB_CRUST_MANTLE_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
                        NGLOB_INNER_CORE_ADJOINT,NSPEC_OUTER_CORE_ROT_ADJOINT, &
                        NSPEC_CRUST_MANTLE_STACEY,NSPEC_OUTER_CORE_STACEY, &
                        NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION,NT_DUMP_ATTENUATION_optimal )

  end program xcreate_header_file

