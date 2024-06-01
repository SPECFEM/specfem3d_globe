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


  subroutine write_AVS_DX_output_adios(npointot,iregion_code, &
                                       num_ibool_AVS_DX, mask_ibool)

  use meshfem_par, only: &
    nspec,ibool,idoubling, &
    xstore,ystore,zstore, &
    myrank,NGLLX,NGLLY,NGLLZ, &
    RICB,RCMB,RTOPDDOUBLEPRIME,R670,R220,R771,R400,R120,R80,RMOHO, &
    RMIDDLE_CRUST, &
    LOCAL_PATH,IMAIN

  use meshfem_models_par, only: &
    ELLIPTICITY,MODEL_3D_MANTLE_PERTUBATIONS, &
    nspl_ellip,rspl_ellip,ellipicity_spline,ellipicity_spline2

  use regions_mesh_par2

  ! Modules for temporary AVS/DX data
  use AVS_DX_global_mod
  use AVS_DX_global_faces_mod
  use AVS_DX_global_chunks_mod
  use AVS_DX_surface_mod

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! number of spectral elements in each block
  integer,intent(in) :: npointot,iregion_code

  ! arrays used for AVS or DX files
  integer, dimension(npointot), intent(inout) :: num_ibool_AVS_DX
  logical, dimension(npointot), intent(inout) :: mask_ibool

  ! local parameters
  ! structures used for ADIOS AVS/DX files
  type(avs_dx_global_t) :: avs_dx_global_vars
  type(avs_dx_global_faces_t) :: avs_dx_global_faces_vars
  type(avs_dx_global_chunks_t) :: avs_dx_global_chunks_vars
  type(avs_dx_surface_t) :: avs_dx_surface_vars

  character(len=MAX_STRING_LEN) :: reg_name, outputname, group_name
  integer(kind=8) :: group_size_inc

  integer :: sizeprocs

  ! number of MPI processes
  call world_size(sizeprocs)

  ! create a prefix for the file name such as LOCAL_PATH/regX_
  call create_name_database_adios(reg_name,iregion_code,LOCAL_PATH)

  write(group_name,"('SPECFEM3D_GLOBE_AVS_DX_reg',i1)") iregion_code

  ! set the adios group size to 0 before incremented by calls to helpers functions.
  group_size_inc = 0
  call init_adios_group(myadios_group,group_name)

  !--- Define ADIOS variables -----------------------------
  call define_AVS_DX_global_data_adios(nspec, ibool, &
                                       npointot, mask_ibool, group_size_inc, avs_dx_global_vars)

  call define_AVS_DX_global_faces_data_adios(nspec, iMPIcut_xi,iMPIcut_eta, &
                                             ibool,mask_ibool,npointot, &
                                             MODEL_3D_MANTLE_PERTUBATIONS, &
                                             group_size_inc, avs_dx_global_faces_vars)

  call define_AVS_DX_global_chunks_data(nspec,iboun,ibool, &
                                        mask_ibool,npointot, &
                                        MODEL_3D_MANTLE_PERTUBATIONS, &
                                        group_size_inc, avs_dx_global_chunks_vars)

  call define_AVS_DX_surfaces_data_adios(nspec,iboun,ibool, &
                                         mask_ibool,npointot, &
                                         MODEL_3D_MANTLE_PERTUBATIONS, &
                                         group_size_inc, avs_dx_surface_vars)


  !--- Open an ADIOS handler to the AVS_DX file. ---------
  outputname = get_adios_filename(trim(reg_name) // "AVS_DX")

  ! user output
  if (myrank == 0) write(IMAIN,*) '    saving arrays in ADIOS file: ',trim(outputname)

  ! opens file for writing
  call open_file_adios_write(myadios_file,myadios_group,outputname,group_name)
  call set_adios_group_size(myadios_file,group_size_inc)

  !--- Schedule writes for the previously defined ADIOS variables
  call prepare_AVS_DX_global_data_adios(nspec, ibool, idoubling, xstore, ystore, zstore, num_ibool_AVS_DX, &
                                        mask_ibool, npointot, avs_dx_global_vars)

  call write_AVS_DX_global_data_adios(myrank,sizeprocs, avs_dx_global_vars)

  call prepare_AVS_DX_global_faces_data_adios(nspec, &
                                              iMPIcut_xi,iMPIcut_eta, &
                                              ibool,idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool, &
                                              npointot,rhostore,kappavstore,muvstore, &
                                              nspl_ellip,rspl_ellip,ellipicity_spline,ellipicity_spline2,ELLIPTICITY, &
                                              MODEL_3D_MANTLE_PERTUBATIONS, &
                                              RICB,RCMB,RTOPDDOUBLEPRIME,R670,R220,R771,R400,R120,R80,RMOHO, &
                                              RMIDDLE_CRUST,iregion_code, &
                                              avs_dx_global_faces_vars)

  call write_AVS_DX_global_faces_data_adios(myrank, sizeprocs, avs_dx_global_faces_vars, MODEL_3D_MANTLE_PERTUBATIONS)

  call prepare_AVS_DX_global_chunks_data_adios(prname,nspec, &
                                               iboun,ibool, idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool, &
                                               npointot,rhostore,kappavstore,muvstore, &
                                               nspl_ellip,rspl_ellip,ellipicity_spline,ellipicity_spline2,ELLIPTICITY, &
                                               MODEL_3D_MANTLE_PERTUBATIONS, &
                                               RICB,RCMB,RTOPDDOUBLEPRIME,R670,R220,R771,R400,R120,R80,RMOHO, &
                                               RMIDDLE_CRUST,iregion_code, &
                                               avs_dx_global_chunks_vars)

  call write_AVS_DX_global_chunks_data_adios(myrank, sizeprocs, avs_dx_global_chunks_vars, MODEL_3D_MANTLE_PERTUBATIONS)

  call prepare_AVS_DX_surfaces_data_adios(nspec,iboun, &
                                          ibool,idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool,npointot, &
                                          rhostore,kappavstore,muvstore, &
                                          nspl_ellip,rspl_ellip,ellipicity_spline,ellipicity_spline2,ELLIPTICITY, &
                                          MODEL_3D_MANTLE_PERTUBATIONS, &
                                          RICB,RCMB,RTOPDDOUBLEPRIME,R670,R220,R771,R400,R120,R80,RMOHO, &
                                          RMIDDLE_CRUST,iregion_code, &
                                          avs_dx_surface_vars)

  call write_AVS_DX_surfaces_data_adios(myrank, sizeprocs, avs_dx_surface_vars, MODEL_3D_MANTLE_PERTUBATIONS)



  !--- Reset the path to zero and perform the actual write to disk
  call write_adios_perform(myadios_file)
  ! closes file
  call close_file_adios(myadios_file)

  !--- Clean up temporary arrays -------------------------
  call free_AVS_DX_global_data_adios(avs_dx_global_vars)
  call free_AVS_DX_global_faces_data_adios(avs_dx_global_faces_vars, MODEL_3D_MANTLE_PERTUBATIONS)
  call free_AVS_DX_global_chunks_data_adios(avs_dx_global_chunks_vars, MODEL_3D_MANTLE_PERTUBATIONS)
  call free_AVS_DX_surfaces_data_adios(avs_dx_surface_vars, MODEL_3D_MANTLE_PERTUBATIONS)

  end subroutine write_AVS_DX_output_adios
