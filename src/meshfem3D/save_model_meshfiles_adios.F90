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

#include "config.fh"


!===============================================================================
!> \brief Save the meshfiles for visualization that will be used by the solver in an ADIOS format.
!!
!! \param myrank The MPI rank of the current process.
!! \param reg_name Output file prefix with the name of the region included
!! \param nspec Number of GLL points per spectral elements
  subroutine save_model_meshfiles_adios()

  use constants

  use meshfem3D_par, only: &
    LOCAL_PATH,nspec,nglob,iregion_code

  use meshfem3D_models_par, only: &
    TRANSVERSE_ISOTROPY,ATTENUATION, &
    ATTENUATION_3D,ATTENUATION_1D_WITH_3D_STORAGE, &
    ANISOTROPIC_3D_MANTLE

  use regions_mesh_par2, only: &
    rhostore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
    Qmu_store,Gc_prime_store,Gs_prime_store,mu0_store

  use adios_write_mod, only: adios_declare_group,adios_select_method
  use adios_helpers_mod, only: define_adios_global_array1D,define_adios_scalar, &
    write_adios_global_1d_array,check_adios_err
  use manager_adios

  implicit none

  ! local parameters
  integer :: i,j,k,ispec
  real(kind=CUSTOM_REAL) :: scaleval1,scaleval2,scaleval,scale_GPa
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: temp_store
  real(kind=CUSTOM_REAL), dimension(1,1,1,1) :: dummy_ijke

  ! local parameters
  character(len=MAX_STRING_LEN) :: outputname, group_name
  integer :: local_dim
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group
  character(len=128)      :: region_name, region_name_scalar
  !--- Save the number of region written. Open the file in "w" mode if 0, else
  !    in "a"  mode
  integer, save :: num_regions_written = 0

  ! scaling factors to re-dimensionalize units
  scaleval1 = sngl( sqrt(PI*GRAV*RHOAV)*(R_EARTH/1000.0d0) )
  scaleval2 = sngl( RHOAV/1000.0d0 )

  ! isotropic model
  write(region_name,"('reg',i1, '/')") iregion_code
  write(region_name_scalar,"('reg',i1)") iregion_code
  write(group_name,"('SPECFEM3D_GLOBE_solver_meshfiles_reg',i1)") iregion_code

  group_size_inc = 0

  call adios_declare_group(adios_group, group_name, '', 1, adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  call adios_select_method(adios_group, ADIOS_TRANSPORT_METHOD, '', '', adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  ! save nspec and nglob, to be used in combine_paraview_data
  call define_adios_scalar (adios_group, group_size_inc, &
                            region_name_scalar, STRINGIFY_VAR(nspec))
  call define_adios_scalar (adios_group, group_size_inc, &
                            region_name_scalar, STRINGIFY_VAR(nglob))

  !--- Define ADIOS variables -----------------------------
  !--- vp arrays -------------------------------------------
  local_dim = size (kappavstore)
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "vp", dummy_ijke)

  !--- vs arrays -------------------------------------------
  local_dim = size (rhostore)
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "vs", dummy_ijke)

  !--- rho arrays ------------------------------------------
  local_dim = size (rhostore)
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "rho", dummy_ijke)

  ! transverse isotropic model
  if (TRANSVERSE_ISOTROPY) then
    !--- vpv arrays ----------------------------------------
    local_dim = size (kappavstore)
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "vpv", dummy_ijke)

    !--- vph arrays ----------------------------------------
    local_dim = size (kappavstore)
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "vph", dummy_ijke)

    !--- vsv arrays ----------------------------------------
    local_dim = size (rhostore)
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "vsv", dummy_ijke)

    !--- vsh arrays ----------------------------------------
    local_dim = size (rhostore)
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "vsh", dummy_ijke)

    !--- eta arrays ----------------------------------------
    local_dim = size (eta_anisostore)
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "eta", eta_anisostore)
  endif

  ! anisotropic values
  if (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
    ! Gc_prime
    local_dim = size (Gc_prime_store)
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "Gc_prime", Gc_prime_store)
    ! Gs_prime
    local_dim = size (Gs_prime_store)
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "Gs_prime", Gs_prime_store)
    ! mu0
    local_dim = size (mu0_store)
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "mu0", mu0_store)
  endif

  if (ATTENUATION) then
    !--- Qmu arrays ----------------------------------------
    local_dim = NGLLX * NGLLY * NGLLZ * nspec
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "qmu", dummy_ijke)
  endif

  !--- Open an ADIOS handler to the restart file. ---------
  outputname = trim(LOCAL_PATH) // "/solver_meshfiles.bp"

  ! user output
  if (myrank == 0) write(IMAIN,*) '    saving arrays in ADIOS file: ',trim(outputname)

  if (num_regions_written == 0) then
    ! opens file for writing
    call open_file_adios_write(outputname,group_name)
  else
    ! opens file for writing in append mode
    call open_file_adios_write_append(outputname,group_name)
  endif
  call set_adios_group_size(group_size_inc)

  ! save nspec and nglob, to be used in combine_paraview_data
  call adios_write(file_handle_adios, trim(region_name) // "nspec", nspec, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_write(file_handle_adios, trim(region_name) // "nglob", nglob, adios_err)
  call check_adios_err(myrank,adios_err)

  !--- Schedule writes for the previously defined ADIOS variables
  ! TODO Try the new write helpers routines
  !--- vp arrays -------------------------------------------
  local_dim = size (kappavstore)
  call write_adios_global_1d_array(file_handle_adios, myrank, sizeprocs_adios, local_dim, trim(region_name) // "vp", &
                                   sqrt( (kappavstore+4.*muvstore/3.)/rhostore )*scaleval1)

  !--- vs arrays -------------------------------------------
  local_dim = size (rhostore)
  call write_adios_global_1d_array(file_handle_adios, myrank, sizeprocs_adios, local_dim, trim(region_name) // "vs", &
                                   sqrt( muvstore/rhostore )*scaleval1 )

  !--- rho arrays ------------------------------------------
  local_dim = size (rhostore)
  call write_adios_global_1d_array(file_handle_adios, myrank, sizeprocs_adios, local_dim, trim(region_name) // "rho", &
                                   rhostore *scaleval2)

  ! transverse isotropic model
  if (TRANSVERSE_ISOTROPY) then
    !--- vps arrays ----------------------------------------
    local_dim = size (kappavstore)
    call write_adios_global_1d_array(file_handle_adios, myrank, sizeprocs_adios, local_dim, trim(region_name) // "vpv", &
                                     sqrt( (kappavstore+4.*muvstore/3.)/rhostore )*scaleval1)

    !--- vph arrays ----------------------------------------
    local_dim = size (kappavstore)
    call write_adios_global_1d_array(file_handle_adios, myrank, sizeprocs_adios, local_dim, trim(region_name) // "vph", &
                                     sqrt( (kappahstore+4.*muhstore/3.)/rhostore )*scaleval1)

    !--- vsv arrays ----------------------------------------
    local_dim = size (rhostore)
    call write_adios_global_1d_array(file_handle_adios, myrank, sizeprocs_adios, local_dim, trim(region_name) // "vsv", &
                                     sqrt( muvstore/rhostore )*scaleval1)

    !--- vsh arrays ----------------------------------------
    local_dim = size (rhostore)
    call write_adios_global_1d_array(file_handle_adios, myrank, sizeprocs_adios, local_dim, trim(region_name) // "vsh", &
                                     sqrt( muhstore/rhostore )*scaleval1)

    !--- eta arrays ----------------------------------------
    local_dim = size (eta_anisostore)
    call write_adios_global_1d_array(file_handle_adios, myrank, sizeprocs_adios, local_dim, trim(region_name) // "eta", &
                                     eta_anisostore)
  endif ! TRANSVERSE_ISOTROPY

  ! anisotropic values
  if (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
    ! the scale of GPa--[g/cm^3][(km/s)^2]
    scaleval = dsqrt(PI*GRAV*RHOAV)
    scale_GPa = (RHOAV/1000.d0)*((R_EARTH*scaleval/1000.d0)**2)

    ! Gc_prime
    local_dim = size (Gc_prime_store)
    call write_adios_global_1d_array(file_handle_adios, myrank, sizeprocs_adios, local_dim, trim(region_name) // "Gc_prime", &
                                     Gc_prime_store)
    ! Gs_prime
    local_dim = size (Gs_prime_store)
    call write_adios_global_1d_array(file_handle_adios, myrank, sizeprocs_adios, local_dim, trim(region_name) // "Gs_prime", &
                                     Gs_prime_store)
    ! mu0
    local_dim = size (mu0_store)
    call write_adios_global_1d_array(file_handle_adios, myrank, sizeprocs_adios, local_dim, trim(region_name) // "mu0", &
                                     mu0_store * scale_GPa)
  endif

  ! shear attenuation
  if (ATTENUATION) then
    !-------------------------------------------------------
    !--- Qmu arrays ----------------------------------------
    !-------------------------------------------------------
    ! saves Qmu_store to full CUSTOM_REAL array
    ! uses temporary array
    allocate(temp_store(NGLLX,NGLLY,NGLLZ,nspec))
    if (ATTENUATION_3D .or. ATTENUATION_1D_WITH_3D_STORAGE) then
      ! attenuation arrays are fully 3D
      temp_store(:,:,:,:) = Qmu_store(:,:,:,:)
    else
      ! attenuation array dimensions: Q_mustore(1,1,1,nspec)
      do ispec = 1,nspec
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              temp_store(i,j,k,ispec) = Qmu_store(1,1,1,ispec)
            enddo
          enddo
        enddo
      enddo
    endif

    local_dim = NGLLX * NGLLY * NGLLZ * nspec
    call write_adios_global_1d_array(file_handle_adios, myrank, sizeprocs_adios, local_dim, trim(region_name) // "qmu", &
                                     temp_store)

    ! frees temporary memory
    deallocate(temp_store)
  endif ! ATTENUATION

  ! closes file
  call close_file_adios()

  !---------------------------------------------------------
  !--- dvpstore arrays ------------------------------------------
  !---------------------------------------------------------
  !obsolete
  !if (HETEROGEN_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
  !  call model_heterogen_mantle_output_dvp_adios(prname)
  !endif

  num_regions_written = num_regions_written + 1

  end subroutine save_model_meshfiles_adios
