!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
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


subroutine crm_save_mesh_files_adios(nspec,npointot,iregion_code, &
    num_ibool_AVS_DX, mask_ibool)
  use mpi
  use adios_write_mod

  use meshfem3d_par,only: &
    ibool,idoubling, &
    xstore,ystore,zstore, &
    myrank,NGLLX,NGLLY,NGLLZ, &
    RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
    RMIDDLE_CRUST,ROCEAN, &
    ADIOS_FOR_AVS_DX, LOCAL_PATH


  use meshfem3D_models_par,only: &
    ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
    nspl,rspl,espl,espl2

  use create_regions_mesh_par2

  ! Modules for temporary AVS/DX data
  use AVS_DX_global_mod
  use AVS_DX_global_faces_mod
  use AVS_DX_global_chunks_mod
  use AVS_DX_surface_mod

  implicit none

  ! number of spectral elements in each block
  integer,intent(in) :: nspec,npointot,iregion_code

  ! local parameters
  ! arrays used for AVS or DX files
  integer, dimension(npointot), intent(inout) :: num_ibool_AVS_DX
  logical, dimension(npointot), intent(inout) :: mask_ibool
  ! structures used for ADIOS AVS/DX files
  type(avs_dx_global_t) :: avs_dx_global_vars
  type(avs_dx_global_faces_t) :: avs_dx_global_faces_vars
  type(avs_dx_global_chunks_t) :: avs_dx_global_chunks_vars
  type(avs_dx_surface_t) :: avs_dx_surface_vars

  character(len=150) :: reg_name, outputname, group_name
  integer :: comm, sizeprocs, ier
  integer(kind=8) :: adios_group, group_size_inc, adios_totalsize, adios_handle

  ! create a prefix for the file name such as LOCAL_PATH/regX_
  call create_name_database_adios(reg_name,iregion_code,LOCAL_PATH)
  outputname = trim(reg_name) // "AVS_DX.bp"
  write(group_name,"('SPECFEM3D_GLOBE_AVS_DX_reg',i1)") iregion_code
  call world_size(sizeprocs) ! TODO keep it in parameters
  ! Alias COMM_WORLD to use ADIOS
  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ier)
  group_size_inc = 0
  call adios_declare_group(adios_group, group_name, &
      "", 0, ier)
  ! We set the transport method to 'MPI'. This seems to be the correct choice
  ! for now. We might want to move this to the constant.h file later on.
  call adios_select_method(adios_group, "MPI", "", "", ier)

  !--- Define ADIOS variables -----------------------------
  call define_AVS_DX_global_data_adios(adios_group, myrank, nspec, ibool, &
      npointot, mask_ibool, group_size_inc, avs_dx_global_vars)

  call define_AVS_DX_global_faces_data_adios (adios_group, &
      myrank, prname, nspec, iMPIcut_xi,iMPIcut_eta, &
      ibool,idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool, &
      npointot,rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
      ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
      RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
      RMIDDLE_CRUST,ROCEAN,iregion_code, &
      group_size_inc, avs_dx_global_faces_vars)

  call define_AVS_DX_global_chunks_data(adios_group, &
      myrank,prname,nspec,iboun,ibool, &
      idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool, &
      npointot,rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
      ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
      RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
      RMIDDLE_CRUST,ROCEAN,iregion_code, &
      group_size_inc, avs_dx_global_chunks_vars)

  call define_AVS_DX_surfaces_data_adios(adios_group, &
      myrank,prname,nspec,iboun, &
      ibool,idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool,npointot,&
      rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
      ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
      RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
      RMIDDLE_CRUST,ROCEAN,iregion_code, &
      group_size_inc, avs_dx_surface_vars)

  !--- Open an ADIOS handler to the AVS_DX file. ---------
  call adios_open (adios_handle, group_name, &
      outputname, "w", comm, ier);
  call adios_group_size (adios_handle, group_size_inc, &
                         adios_totalsize, ier)

  !--- Schedule writes for the previously defined ADIOS variables
  call prepare_AVS_DX_global_data_adios(adios_handle, myrank, &
      nspec, ibool, idoubling, xstore, ystore, zstore, num_ibool_AVS_DX, &
      mask_ibool, npointot, avs_dx_global_vars)
  call write_AVS_DX_global_data_adios(adios_handle, myrank, &
      sizeprocs, avs_dx_global_vars)

  call prepare_AVS_DX_global_faces_data_adios (myrank, prname, nspec, &
      iMPIcut_xi,iMPIcut_eta, &
      ibool,idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool, &
      npointot,rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
      ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
      RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
      RMIDDLE_CRUST,ROCEAN,iregion_code, &
      avs_dx_global_faces_vars)
  call write_AVS_DX_global_faces_data_adios(adios_handle, myrank, &
      sizeprocs, avs_dx_global_faces_vars, ISOTROPIC_3D_MANTLE)

  call prepare_AVS_DX_global_chunks_data_adios(myrank,prname,nspec, &
      iboun,ibool, idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool,&
      npointot,rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
      ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
      RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
      RMIDDLE_CRUST,ROCEAN,iregion_code, &
      avs_dx_global_chunks_vars)
  call write_AVS_DX_global_chunks_data_adios(adios_handle, myrank, &
      sizeprocs, avs_dx_global_chunks_vars, ISOTROPIC_3D_MANTLE)

  call prepare_AVS_DX_surfaces_data_adios(myrank,prname,nspec,iboun, &
      ibool,idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool,npointot,&
      rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
      ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
      RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
      RMIDDLE_CRUST,ROCEAN,iregion_code, &
      avs_dx_surface_vars)
  call write_AVS_DX_surfaces_data_adios(adios_handle, myrank, &
      sizeprocs, avs_dx_surface_vars, ISOTROPIC_3D_MANTLE)

  !--- Reset the path to zero and perform the actual write to disk
  call adios_set_path (adios_handle, "", ier)
  call adios_close(adios_handle, ier)

  !--- Clean up temporary arrays -------------------------
  call free_AVS_DX_global_data_adios(myrank, avs_dx_global_vars)
  call free_AVS_DX_global_faces_data_adios(myrank, avs_dx_global_faces_vars, &
      ISOTROPIC_3D_MANTLE)
  call free_AVS_DX_global_chunks_data_adios(myrank, avs_dx_global_chunks_vars, &
      ISOTROPIC_3D_MANTLE)
  call free_AVS_DX_surfaces_data_adios(myrank, avs_dx_surface_vars, &
      ISOTROPIC_3D_MANTLE)
end subroutine crm_save_mesh_files_adios
