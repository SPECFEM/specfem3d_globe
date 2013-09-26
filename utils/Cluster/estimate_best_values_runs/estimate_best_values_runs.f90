!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            August 2013
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

! output a table with the best parameters to use for large runs on a given machine
! in order to use the total amount of memory available in the best possible way

  program estimate_best_values_runs

  use constants
  use shared_parameters

  implicit none

  integer, parameter :: NB_COLUMNS_TABLE = 10

! maximum total number of processors we want to see in the table
! integer, parameter :: MAX_NUMBER_OF_PROCS = 100000
! integer, parameter :: MAX_NUMBER_OF_PROCS =  4000  ! current maximum on pangu at Caltech
! integer, parameter :: MAX_NUMBER_OF_PROCS = 10240  ! current maximum on MareNostrum in Barcelona
! integer, parameter :: MAX_NUMBER_OF_PROCS = 62976  ! current maximum on Ranger in Texas
! integer, parameter :: MAX_NUMBER_OF_PROCS = 212992 ! current maximum on BlueGene at LLNL
  integer, parameter :: MAX_NUMBER_OF_PROCS = 12150  ! current maximum at CINES in France

! minimum total number of processors we want to see in the table
! integer, parameter :: MIN_NUMBER_OF_PROCS =  100
  integer, parameter :: MIN_NUMBER_OF_PROCS = 1000

! amount of memory installed per processor core on the system (in gigabytes)
! double precision, parameter :: MAX_MEMORY_PER_CORE = 1.5d0 ! pangu at Caltech
! double precision, parameter :: MAX_MEMORY_PER_CORE = 2.0d0 ! Ranger in Texas, MareNostrum in Barcelona
  double precision, parameter :: MAX_MEMORY_PER_CORE = 4.0d0 ! CINES in France

! base value depends if we implement three or four doublings (default is three)
  integer, parameter :: NB_DOUBLING = 3
  integer :: BASE_VALUE

  double precision :: static_memory_size

  integer :: NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
         NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUAT, &
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

  logical :: already_printed
  integer :: c,counter,offset

  double precision :: mem_per_core,percent

! ************** PROGRAM STARTS HERE **************

! count the total number of sources in the CMTSOLUTION file
  NSOURCES = 1

  BASE_VALUE = 2**NB_DOUBLING

  first_time = .true.

  print *
  print *,'Number of GB of memory per core installed in the machine: ',MAX_MEMORY_PER_CORE
  print *

! make sure we do not start below the lower limit
  if(6 * int(sqrt(MIN_NUMBER_OF_PROCS / 6.d0))**2 == MIN_NUMBER_OF_PROCS) then
    offset = 0
  else
    offset = 1
  endif

! the total number of processors is 6 * NPROC_XI^2
  do NPROC_XI = int(sqrt(MIN_NUMBER_OF_PROCS / 6.d0)) + offset,int(sqrt(MAX_NUMBER_OF_PROCS / 6.d0))

    counter = 1
    c = 0
    NEX_XI = 0
    already_printed = .false.

! loop on the columns of the table
    do while (counter <= NB_COLUMNS_TABLE)

      c = c + 1
      NEX_XI = BASE_VALUE * c * NPROC_XI

! impose an upper limit on current "reasonably large" values to avoid an infinite loop
! these values should be okay for the next few years (i.e., on machines around 1 petaflops)
! (this written in August 2008)
      if(c > 15 .or. NEX_XI > 6000) exit

! also exclude very coarse meshes below NEX_XI = 64
      if(NEX_XI >= 64 .and. mod(NEX_XI,2*BASE_VALUE) == 0) then

        counter = counter + 1

! read the parameter file and compute additional parameters
        EMULATE_ONLY = .true.
        call read_compute_parameters()

        if(first_time) then
          first_time = .false.
          if(ATTENUATION) then
            print *,'ATTENUATION = .true.'
          else
            print *,'ATTENUATION = .false.'
          endif
        endif

!! compute memory usage per processor core
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

          mem_per_core = static_memory_size/1073741824.d0

          percent = 100.d0 * mem_per_core / MAX_MEMORY_PER_CORE

          if(percent > 100.d0) goto 777

  if(percent < 0.d0) then
    if(.not. already_printed) then
      print *
      print *,'NPROC, % of the machine, NPROC_XI, NEX_XI, GB used/core, % mem used/core:'
    endif
    already_printed = .true.
    write(*,"(' ',i5,'  ',f6.2,'% ',i3,'  ',i4,'  ',f6.2,'GB  ',f6.2,'% **mesher fails/bug**')") &
      6*NPROC_XI**2,100.d0*6*NPROC_XI**2/dble(MAX_NUMBER_OF_PROCS),NPROC_XI,NEX_XI,mem_per_core,percent
    goto 777

  else if(percent >= 93.d0) then
    if(.not. already_printed) then
      print *
      print *,'NPROC, % of the machine, NPROC_XI, NEX_XI, GB used/core, % mem used/core:'
    endif
    already_printed = .true.
    write(*,"(' ',i5,'  ',f6.2,'% ',i3,'  ',i4,'  ',f6.2,'GB  ',f6.2,'% **too high**')") &
      6*NPROC_XI**2,100.d0*6*NPROC_XI**2/dble(MAX_NUMBER_OF_PROCS),NPROC_XI,NEX_XI,mem_per_core,percent
    goto 777

  else if(percent >= 85.d0) then
    if(.not. already_printed) then
      print *
      print *,'NPROC, % of the machine, NPROC_XI, NEX_XI, GB used/core, % mem used/core:'
    endif
    already_printed = .true.
    write(*,"(' ',i5,'  ',f6.2,'% ',i3,'  ',i4,'  ',f6.2,'GB  ',f6.2,'% **excellent**')") &
      6*NPROC_XI**2,100.d0*6*NPROC_XI**2/dble(MAX_NUMBER_OF_PROCS),NPROC_XI,NEX_XI,mem_per_core,percent
    goto 777

  else if(percent >= 80.d0) then
    if(.not. already_printed) then
      print *
      print *,'NPROC, % of the machine, NPROC_XI, NEX_XI, GB used/core, % mem used/core:'
    endif
    already_printed = .true.
    write(*,"(' ',i5,'  ',f6.2,'% ',i3,'  ',i4,'  ',f6.2,'GB  ',f6.2,'% **very good**')") &
      6*NPROC_XI**2,100.d0*6*NPROC_XI**2/dble(MAX_NUMBER_OF_PROCS),NPROC_XI,NEX_XI,mem_per_core,percent
    goto 777

  else if(percent >= 70.d0) then
    if(.not. already_printed) then
      print *
      print *,'NPROC, % of the machine, NPROC_XI, NEX_XI, GB used/core, % mem used/core:'
    endif
    already_printed = .true.
    write(*,"(' ',i5,'  ',f6.2,'% ',i3,'  ',i4,'  ',f6.2,'GB  ',f6.2,'% **not bad**')") &
      6*NPROC_XI**2,100.d0*6*NPROC_XI**2/dble(MAX_NUMBER_OF_PROCS),NPROC_XI,NEX_XI,mem_per_core,percent
    goto 777

  else if(percent >= 60.d0) then
    if(.not. already_printed) then
      print *
      print *,'NPROC, % of the machine, NPROC_XI, NEX_XI, GB used/core, % mem used/core:'
    endif
    already_printed = .true.
    write(*,"(' ',i5,'  ',f6.2,'% ',i3,'  ',i4,'  ',f6.2,'GB  ',f6.2,'% **possible**')") &
      6*NPROC_XI**2,100.d0*6*NPROC_XI**2/dble(MAX_NUMBER_OF_PROCS),NPROC_XI,NEX_XI,mem_per_core,percent

! uncomment line below to completely suppress cases that give an occupancy below 50%
  else ! if(percent >= 50.d0) then
    if(.not. already_printed) then
      print *
      print *,'NPROC, % of the machine, NPROC_XI, NEX_XI, GB used/core, % mem used/core:'
    endif
    already_printed = .true.
    write(*,"(' ',i5,'  ',f6.2,'% ',i3,'  ',i4,'  ',f6.2,'GB  ',f6.2,'%')") &
      6*NPROC_XI**2,100.d0*6*NPROC_XI**2/dble(MAX_NUMBER_OF_PROCS),NPROC_XI,NEX_XI,mem_per_core,percent
  endif

      endif
    enddo

 777 continue

  enddo

  print *

  end program estimate_best_values_runs

!
!----  include subroutines from the main code
!

  include "../../read_compute_parameters.f90"

  include "../../memory_eval.f90"

  include "../../read_value_parameters.f90"

  include "../../get_value_parameters.f90"

  include "../../auto_ner.f90"

