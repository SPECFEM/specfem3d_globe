!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            April 2011
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

  subroutine read_forward_arrays_startrun(myrank,NSTEP, &
                    SIMULATION_TYPE,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN, &
                    it_begin,it_end, &
                    displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle, &
                    displ_inner_core,veloc_inner_core,accel_inner_core, &
                    displ_outer_core,veloc_outer_core,accel_outer_core, &
                    R_memory_crust_mantle,R_memory_inner_core, &
                    A_array_rotation,B_array_rotation, &
                    b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle, &
                    b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core, &
                    b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core, &
                    b_R_memory_crust_mantle,b_R_memory_inner_core, &
                    b_A_array_rotation,b_B_array_rotation,LOCAL_PATH)

! reads in saved wavefields

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank,NSTEP

  integer SIMULATION_TYPE

  integer NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,it_begin,it_end

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: &
    displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: &
    displ_inner_core,veloc_inner_core,accel_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: &
    displ_outer_core,veloc_outer_core,accel_outer_core

  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: &
    R_memory_crust_mantle
  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ATTENUATION) :: &
    R_memory_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROTATION) :: &
    A_array_rotation,B_array_rotation

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE_ADJOINT) :: &
    b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE_ADJOINT) :: &
    b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE_ADJOINT) :: &
    b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core

  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_AND_ATT) :: &
    b_R_memory_crust_mantle

  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STR_AND_ATT) :: &
    b_R_memory_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROT_ADJOINT) :: &
    b_A_array_rotation,b_B_array_rotation

  character(len=150) LOCAL_PATH

  !local parameters
  character(len=150) outputname

  ! define correct time steps if restart files
  if(NUMBER_OF_RUNS < 1 .or. NUMBER_OF_RUNS > NSTEP) &
    stop 'number of restart runs can not be less than 1 or greater than NSTEP'
  if(NUMBER_OF_THIS_RUN < 1 .or. NUMBER_OF_THIS_RUN > NUMBER_OF_RUNS) stop 'incorrect run number'
  if (SIMULATION_TYPE /= 1 .and. NUMBER_OF_RUNS /= 1) stop 'Only 1 run for SIMULATION_TYPE = 2/3'

  it_begin = (NUMBER_OF_THIS_RUN - 1) * (NSTEP / NUMBER_OF_RUNS) + 1
  if (NUMBER_OF_THIS_RUN < NUMBER_OF_RUNS) then
    it_end = NUMBER_OF_THIS_RUN * (NSTEP / NUMBER_OF_RUNS)
  else
    ! Last run may be a bit larger
    it_end = NSTEP
  endif

  ! read files back from local disk or MT tape system if restart file
  if(NUMBER_OF_THIS_RUN > 1) then
    write(outputname,"('dump_all_arrays',i6.6)") myrank
    open(unit=55,file=trim(LOCAL_PATH)//'/'//outputname,status='old',action='read',form='unformatted')
    read(55) displ_crust_mantle
    read(55) veloc_crust_mantle
    read(55) accel_crust_mantle
    read(55) displ_inner_core
    read(55) veloc_inner_core
    read(55) accel_inner_core
    read(55) displ_outer_core
    read(55) veloc_outer_core
    read(55) accel_outer_core
    read(55) A_array_rotation
    read(55) B_array_rotation
    read(55) R_memory_crust_mantle
    read(55) R_memory_inner_core
    close(55)
  endif

  if (SIMULATION_TYPE == 3) then
    ! initializes
    b_displ_crust_mantle = 0._CUSTOM_REAL
    b_veloc_crust_mantle = 0._CUSTOM_REAL
    b_accel_crust_mantle = 0._CUSTOM_REAL
    b_displ_inner_core = 0._CUSTOM_REAL
    b_veloc_inner_core = 0._CUSTOM_REAL
    b_accel_inner_core = 0._CUSTOM_REAL
    b_displ_outer_core = 0._CUSTOM_REAL
    b_veloc_outer_core = 0._CUSTOM_REAL
    b_accel_outer_core = 0._CUSTOM_REAL
    if (ROTATION_VAL) then
      b_A_array_rotation = 0._CUSTOM_REAL
      b_B_array_rotation = 0._CUSTOM_REAL
    endif
    if (ATTENUATION_VAL) then
      b_R_memory_crust_mantle = 0._CUSTOM_REAL
      b_R_memory_inner_core = 0._CUSTOM_REAL
    endif
  endif

  end subroutine read_forward_arrays_startrun

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_forward_arrays(myrank, &
                    b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle, &
                    b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core, &
                    b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core, &
                    b_R_memory_crust_mantle,b_R_memory_inner_core, &
                    b_A_array_rotation,b_B_array_rotation,LOCAL_PATH)

! reads in saved wavefields

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank

  ! backward/reconstructed wavefields
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE_ADJOINT) :: &
    b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE_ADJOINT) :: &
    b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE_ADJOINT) :: &
    b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core

  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_AND_ATT) :: &
    b_R_memory_crust_mantle

  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STR_AND_ATT) :: &
    b_R_memory_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROT_ADJOINT) :: &
    b_A_array_rotation,b_B_array_rotation

  character(len=150) LOCAL_PATH

  !local parameters
  character(len=150) outputname

  write(outputname,'(a,i6.6,a)') 'proc',myrank,'_save_forward_arrays.bin'
  open(unit=55,file=trim(LOCAL_PATH)//'/'//outputname,status='old',action='read',form='unformatted')
  read(55) b_displ_crust_mantle
  read(55) b_veloc_crust_mantle
  read(55) b_accel_crust_mantle
  read(55) b_displ_inner_core
  read(55) b_veloc_inner_core
  read(55) b_accel_inner_core
  read(55) b_displ_outer_core
  read(55) b_veloc_outer_core
  read(55) b_accel_outer_core
  if (ROTATION_VAL) then
    read(55) b_A_array_rotation
    read(55) b_B_array_rotation
  endif
  if (ATTENUATION_VAL) then
    read(55) b_R_memory_crust_mantle
    read(55) b_R_memory_inner_core
  endif
  close(55)

  end subroutine read_forward_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_forward_arrays_undoatt(myrank, &
                    b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle, &
                    b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core, &
                    b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core, &
                    b_R_memory_crust_mantle,b_R_memory_inner_core, &
                    b_A_array_rotation,b_B_array_rotation,LOCAL_PATH,iteration_on_subset)

! reads in saved wavefields

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank

  ! backward/reconstructed wavefields
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE_ADJOINT) :: &
    b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE_ADJOINT) :: &
    b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE_ADJOINT) :: &
    b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core

  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_AND_ATT) :: &
    b_R_memory_crust_mantle

  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STR_AND_ATT) :: &
    b_R_memory_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROT_ADJOINT) :: &
    b_A_array_rotation,b_B_array_rotation

  character(len=150) LOCAL_PATH

  integer iteration_on_subset

  !local parameters
  character(len=150) outputname

  write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_save_frame_at',iteration_on_subset,'.bin'
  open(unit=55,file=trim(LOCAL_PATH)//'/'//outputname,status='old',action='read',form='unformatted')
  read(55) b_displ_crust_mantle
  read(55) b_veloc_crust_mantle
  read(55) b_accel_crust_mantle
  read(55) b_displ_inner_core
  read(55) b_veloc_inner_core
  read(55) b_accel_inner_core
  read(55) b_displ_outer_core
  read(55) b_veloc_outer_core
  read(55) b_accel_outer_core
  if (ROTATION_VAL) then
    read(55) b_A_array_rotation
    read(55) b_B_array_rotation
  endif
  if (ATTENUATION_VAL) then
    read(55) b_R_memory_crust_mantle
    read(55) b_R_memory_inner_core
  endif
  close(55)

  end subroutine read_forward_arrays_undoatt
