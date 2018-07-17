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


  subroutine write_AVS_DX_output(npointot,iregion_code)

  use meshfem3d_par, only: &
    nspec,ibool,idoubling, &
    xstore,ystore,zstore, &
    NGLLX,NGLLY,NGLLZ, &
    RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
    RMIDDLE_CRUST,ROCEAN, &
    ADIOS_FOR_AVS_DX


  use meshfem3D_models_par, only: &
    ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
    nspl,rspl,espl,espl2

  use regions_mesh_par2

  ! Modules for temporary AVS/DX data
!  use AVS_DX_global_mod

  implicit none

  ! number of spectral elements in each block
  integer,intent(in) :: npointot,iregion_code

  ! local parameters
  ! arrays used for AVS or DX files
  integer, dimension(:), allocatable :: num_ibool_AVS_DX
  logical, dimension(:), allocatable :: mask_ibool
  integer :: ier

  ! arrays num_ibool_AVS_DX and mask_ibool used to save memory
  ! allocate memory for arrays
  allocate(num_ibool_AVS_DX(npointot), &
           mask_ibool(npointot), &
           stat=ier)
  if (ier /= 0) stop 'Error in allocate 21'

  if (ADIOS_FOR_AVS_DX) then
    call write_AVS_DX_output_adios(npointot,iregion_code, &
                                   num_ibool_AVS_DX, mask_ibool)
  else
    call write_AVS_DX_global_data(prname,nspec,ibool,idoubling, &
        xstore,ystore,zstore, num_ibool_AVS_DX,mask_ibool,npointot)

    call write_AVS_DX_global_faces_data(prname,nspec,iMPIcut_xi, &
        iMPIcut_eta,ibool, idoubling,xstore,ystore,zstore,num_ibool_AVS_DX, &
        mask_ibool,npointot, rhostore,kappavstore,muvstore,nspl,rspl, &
        espl,espl2, ELLIPTICITY,ISOTROPIC_3D_MANTLE, RICB,RCMB, &
        RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
        RMIDDLE_CRUST,ROCEAN,iregion_code)

    call write_AVS_DX_global_chunks_data(prname,nspec,iboun,ibool, &
            idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool, &
            npointot,rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
            ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
            RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
            RMIDDLE_CRUST,ROCEAN,iregion_code)

    call write_AVS_DX_surface_data(prname,nspec,iboun,ibool, &
            idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool,npointot, &
            rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
            ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
            RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
            RMIDDLE_CRUST,ROCEAN,iregion_code)
  endif

  ! Output material information for all GLL points
  ! Can be use to check the mesh
  !    call write_AVS_DX_global_data_gll(prname,nspec,xstore,ystore,zstore, &
  !                rhostore,kappavstore,muvstore,Qmu_store,ATTENUATION)
  deallocate(num_ibool_AVS_DX,mask_ibool)

  end subroutine write_AVS_DX_output

