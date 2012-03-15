!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
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

  subroutine save_forward_arrays()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! local parameters
  character(len=150) outputname

  ! save files to local disk or tape system if restart file
  if(NUMBER_OF_RUNS > 1 .and. NUMBER_OF_THIS_RUN < NUMBER_OF_RUNS) then
    write(outputname,"('dump_all_arrays',i6.6)") myrank
    open(unit=55,file=trim(LOCAL_TMP_PATH)//'/'//outputname,status='unknown',form='unformatted',action='write')

    write(55) displ_crust_mantle
    write(55) veloc_crust_mantle
    write(55) accel_crust_mantle
    write(55) displ_inner_core
    write(55) veloc_inner_core
    write(55) accel_inner_core
    write(55) displ_outer_core
    write(55) veloc_outer_core
    write(55) accel_outer_core

    write(55) epsilondev_xx_crust_mantle
    write(55) epsilondev_yy_crust_mantle
    write(55) epsilondev_xy_crust_mantle
    write(55) epsilondev_xz_crust_mantle
    write(55) epsilondev_yz_crust_mantle

    write(55) epsilondev_xx_inner_core
    write(55) epsilondev_yy_inner_core
    write(55) epsilondev_xy_inner_core
    write(55) epsilondev_xz_inner_core
    write(55) epsilondev_yz_inner_core

    write(55) A_array_rotation
    write(55) B_array_rotation

    write(55) R_xx_crust_mantle
    write(55) R_yy_crust_mantle
    write(55) R_xy_crust_mantle
    write(55) R_xz_crust_mantle
    write(55) R_yz_crust_mantle

    write(55) R_xx_inner_core
    write(55) R_yy_inner_core
    write(55) R_xy_inner_core
    write(55) R_xz_inner_core
    write(55) R_yz_inner_core

    close(55)
  endif

  ! save last frame of the forward simulation
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_save_forward_arrays.bin'
    open(unit=55,file=trim(LOCAL_TMP_PATH)//'/'//outputname,status='unknown',form='unformatted',action='write')

    write(55) displ_crust_mantle
    write(55) veloc_crust_mantle
    write(55) accel_crust_mantle
    write(55) displ_inner_core
    write(55) veloc_inner_core
    write(55) accel_inner_core
    write(55) displ_outer_core
    write(55) veloc_outer_core
    write(55) accel_outer_core

    write(55) epsilondev_xx_crust_mantle
    write(55) epsilondev_yy_crust_mantle
    write(55) epsilondev_xy_crust_mantle
    write(55) epsilondev_xz_crust_mantle
    write(55) epsilondev_yz_crust_mantle

    write(55) epsilondev_xx_inner_core
    write(55) epsilondev_yy_inner_core
    write(55) epsilondev_xy_inner_core
    write(55) epsilondev_xz_inner_core
    write(55) epsilondev_yz_inner_core

    if (ROTATION_VAL) then
      write(55) A_array_rotation
      write(55) B_array_rotation
    endif
    if (ATTENUATION_VAL) then
      write(55) R_xx_crust_mantle
      write(55) R_yy_crust_mantle
      write(55) R_xy_crust_mantle
      write(55) R_xz_crust_mantle
      write(55) R_yz_crust_mantle

      write(55) R_xx_inner_core
      write(55) R_yy_inner_core
      write(55) R_xy_inner_core
      write(55) R_xz_inner_core
      write(55) R_yz_inner_core

    endif
    close(55)
  endif

  end subroutine save_forward_arrays
