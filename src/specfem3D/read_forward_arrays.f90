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

  subroutine read_forward_arrays_startrun()

! reads in saved wavefields

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! local parameters
  character(len=150) outputname

  ! define correct time steps if restart files
  if(NUMBER_OF_RUNS < 1 .or. NUMBER_OF_RUNS > 3) stop 'number of restart runs can be 1, 2 or 3'
  if(NUMBER_OF_THIS_RUN < 1 .or. NUMBER_OF_THIS_RUN > NUMBER_OF_RUNS) stop 'incorrect run number'
  if (SIMULATION_TYPE /= 1 .and. NUMBER_OF_RUNS /= 1) stop 'Only 1 run for SIMULATION_TYPE = 2/3'

  ! set start/end steps for time iteration loop
  select case( NUMBER_OF_RUNS )
  case( 3 )
    if(NUMBER_OF_THIS_RUN == 1) then
      it_begin = 1
      it_end = NSTEP/3
    else if(NUMBER_OF_THIS_RUN == 2) then
      it_begin = NSTEP/3 + 1
      it_end = 2*(NSTEP/3)
    else
      it_begin = 2*(NSTEP/3) + 1
      it_end = NSTEP
    endif
  case( 2 )
    if(NUMBER_OF_THIS_RUN == 1) then
      it_begin = 1
      it_end = NSTEP/2
    else
      it_begin = NSTEP/2 + 1
      it_end = NSTEP
    endif
  case default
    it_begin = 1
    it_end = NSTEP
  end select

  ! read files back from local disk or MT tape system if restart file
  if(NUMBER_OF_THIS_RUN > 1) then
    if (ADIOS_FOR_FORWARD_ARRAYS) then
      call read_intermediate_forward_arrays_adios()
    else
      write(outputname,"('dump_all_arrays',i6.6)") myrank
      open(unit=55,file=trim(LOCAL_TMP_PATH)//'/'//outputname,status='old',action='read',form='unformatted')

      read(55) displ_crust_mantle
      read(55) veloc_crust_mantle
      read(55) accel_crust_mantle
      read(55) displ_inner_core
      read(55) veloc_inner_core
      read(55) accel_inner_core
      read(55) displ_outer_core
      read(55) veloc_outer_core
      read(55) accel_outer_core

      read(55) epsilondev_xx_crust_mantle
      read(55) epsilondev_yy_crust_mantle
      read(55) epsilondev_xy_crust_mantle
      read(55) epsilondev_xz_crust_mantle
      read(55) epsilondev_yz_crust_mantle

      read(55) epsilondev_xx_inner_core
      read(55) epsilondev_yy_inner_core
      read(55) epsilondev_xy_inner_core
      read(55) epsilondev_xz_inner_core
      read(55) epsilondev_yz_inner_core

      read(55) A_array_rotation
      read(55) B_array_rotation

      read(55) R_xx_crust_mantle
      read(55) R_yy_crust_mantle
      read(55) R_xy_crust_mantle
      read(55) R_xz_crust_mantle
      read(55) R_yz_crust_mantle

      read(55) R_xx_inner_core
      read(55) R_yy_inner_core
      read(55) R_xy_inner_core
      read(55) R_xz_inner_core
      read(55) R_yz_inner_core

      close(55)
    endif
  endif

  ! initializes backward/reconstructed arrays
  if (SIMULATION_TYPE == 3) then
    ! initializes wavefields
    b_displ_crust_mantle = 0._CUSTOM_REAL
    b_veloc_crust_mantle = 0._CUSTOM_REAL
    b_accel_crust_mantle = 0._CUSTOM_REAL

    b_displ_inner_core = 0._CUSTOM_REAL
    b_veloc_inner_core = 0._CUSTOM_REAL
    b_accel_inner_core = 0._CUSTOM_REAL

    b_displ_outer_core = 0._CUSTOM_REAL
    b_veloc_outer_core = 0._CUSTOM_REAL
    b_accel_outer_core = 0._CUSTOM_REAL

    b_epsilondev_xx_crust_mantle = 0._CUSTOM_REAL
    b_epsilondev_yy_crust_mantle = 0._CUSTOM_REAL
    b_epsilondev_xy_crust_mantle = 0._CUSTOM_REAL
    b_epsilondev_xz_crust_mantle = 0._CUSTOM_REAL
    b_epsilondev_yz_crust_mantle = 0._CUSTOM_REAL

    b_epsilondev_xx_inner_core = 0._CUSTOM_REAL
    b_epsilondev_yy_inner_core = 0._CUSTOM_REAL
    b_epsilondev_xy_inner_core = 0._CUSTOM_REAL
    b_epsilondev_xz_inner_core = 0._CUSTOM_REAL
    b_epsilondev_yz_inner_core = 0._CUSTOM_REAL

    if (ROTATION_VAL) then
      b_A_array_rotation = 0._CUSTOM_REAL
      b_B_array_rotation = 0._CUSTOM_REAL
    endif

    if (ATTENUATION_VAL) then
      b_R_xx_crust_mantle = 0._CUSTOM_REAL
      b_R_yy_crust_mantle = 0._CUSTOM_REAL
      b_R_xy_crust_mantle = 0._CUSTOM_REAL
      b_R_xz_crust_mantle = 0._CUSTOM_REAL
      b_R_yz_crust_mantle = 0._CUSTOM_REAL

      b_R_xx_inner_core = 0._CUSTOM_REAL
      b_R_yy_inner_core = 0._CUSTOM_REAL
      b_R_xy_inner_core = 0._CUSTOM_REAL
      b_R_xz_inner_core = 0._CUSTOM_REAL
      b_R_yz_inner_core = 0._CUSTOM_REAL
    endif
  endif

  end subroutine read_forward_arrays_startrun

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_forward_arrays()

! reads in saved wavefields to reconstruct/backward wavefield

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! local parameters
  integer :: ier
  character(len=150) outputname

  if (ADIOS_FOR_FORWARD_ARRAYS) then
    call read_forward_arrays_adios()
  else
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_save_forward_arrays.bin'
    open(unit=55,file=trim(LOCAL_TMP_PATH)//'/'//outputname, &
          status='old',action='read',form='unformatted',iostat=ier)
    if( ier /= 0 ) then
      print*,'error: opening proc_****_save_forward_arrays.bin'
      print*,'path: ',trim(LOCAL_TMP_PATH)//'/'//outputname
      call exit_mpi(myrank,'error open file save_forward_arrays.bin')
    endif

    read(55) b_displ_crust_mantle
    read(55) b_veloc_crust_mantle
    read(55) b_accel_crust_mantle

    read(55) b_displ_inner_core
    read(55) b_veloc_inner_core
    read(55) b_accel_inner_core

    read(55) b_displ_outer_core
    read(55) b_veloc_outer_core
    read(55) b_accel_outer_core

    read(55) b_epsilondev_xx_crust_mantle
    read(55) b_epsilondev_yy_crust_mantle
    read(55) b_epsilondev_xy_crust_mantle
    read(55) b_epsilondev_xz_crust_mantle
    read(55) b_epsilondev_yz_crust_mantle

    read(55) b_epsilondev_xx_inner_core
    read(55) b_epsilondev_yy_inner_core
    read(55) b_epsilondev_xy_inner_core
    read(55) b_epsilondev_xz_inner_core
    read(55) b_epsilondev_yz_inner_core

    if (ROTATION_VAL) then
      read(55) b_A_array_rotation
      read(55) b_B_array_rotation
    endif

    if (ATTENUATION_VAL) then
       read(55) b_R_xx_crust_mantle
       read(55) b_R_yy_crust_mantle
       read(55) b_R_xy_crust_mantle
       read(55) b_R_xz_crust_mantle
       read(55) b_R_yz_crust_mantle

       read(55) b_R_xx_inner_core
       read(55) b_R_yy_inner_core
       read(55) b_R_xy_inner_core
       read(55) b_R_xz_inner_core
       read(55) b_R_yz_inner_core
    endif
    close(55)
  endif ! ADIOS_FOR_FORWARD_ARRAYS

  ! transfers fields onto GPU
  if(GPU_MODE) then
    call transfer_b_fields_cm_to_device(NDIM*NGLOB_CRUST_MANTLE, &
                                    b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle, &
                                    Mesh_pointer)

    call transfer_b_fields_ic_to_device(NDIM*NGLOB_INNER_CORE, &
                                    b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core, &
                                    Mesh_pointer)

    call transfer_b_fields_oc_to_device(NGLOB_OUTER_CORE, &
                                    b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core, &
                                    Mesh_pointer)

    call transfer_b_strain_cm_to_device(Mesh_pointer, &
                                    b_epsilondev_xx_crust_mantle,b_epsilondev_yy_crust_mantle, &
                                    b_epsilondev_xy_crust_mantle,b_epsilondev_xz_crust_mantle, &
                                    b_epsilondev_yz_crust_mantle)

    call transfer_b_strain_ic_to_device(Mesh_pointer, &
                                    b_epsilondev_xx_inner_core,b_epsilondev_yy_inner_core, &
                                    b_epsilondev_xy_inner_core,b_epsilondev_xz_inner_core, &
                                    b_epsilondev_yz_inner_core)
    if (ROTATION_VAL) then
      call transfer_b_rotation_to_device(Mesh_pointer,b_A_array_rotation,b_B_array_rotation)
    endif

    if (ATTENUATION_VAL) then
      call transfer_b_rmemory_cm_to_device(Mesh_pointer, &
                                    b_R_xx_crust_mantle,b_R_yy_crust_mantle, &
                                    b_R_xy_crust_mantle,b_R_xz_crust_mantle, &
                                    b_R_yz_crust_mantle)
      call transfer_b_rmemory_ic_to_device(Mesh_pointer, &
                                    b_R_xx_inner_core,b_R_yy_inner_core, &
                                    b_R_xy_inner_core,b_R_xz_inner_core, &
                                    b_R_yz_inner_core)
    endif
  endif

  end subroutine read_forward_arrays
