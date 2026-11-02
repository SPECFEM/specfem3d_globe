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


  subroutine save_model_meshfiles_hdf5()

  use constants
#ifdef USE_HDF5
! outputs model files in binary format

  use shared_parameters, only: R_PLANET,RHOAV,LOCAL_PATH,H5_COL

  use meshfem_par, only: nspec,iregion_code,NPROCTOT

  use meshfem_models_par, only: &
    TRANSVERSE_ISOTROPY,ATTENUATION, &
    ATTENUATION_3D,ATTENUATION_3D_BERKELEY,ATTENUATION_1D_WITH_3D_STORAGE, &
    HETEROGEN_3D_MANTLE,ANISOTROPIC_3D_MANTLE

  use regions_mesh_par2, only: &
    rhostore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
    Qmu_store,Gc_prime_store,Gs_prime_store,mu0store

  use model_heterogen_mantle_par

  use manager_hdf5
#endif

  implicit none

#ifdef USE_HDF5
  ! local parameters
  integer :: i,j,k,ispec,ier
  real(kind=CUSTOM_REAL) :: scaleval1,scaleval2,scaleval,scale_GPa
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: temp_store

  ! dset_name and group_name
  character(len=64) :: dset_name, gname_region

  ! offset for nspec
  integer, dimension(0:NPROCTOT-1) :: offset_nelems

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     model in HDF5 file format'
    call flush_IMAIN()
  endif

  ! gather nspec
  call gather_all_all_singlei(nspec, offset_nelems, NPROCTOT)

  ! file name
  name_database_hdf5 = LOCAL_PATH(1:len_trim(LOCAL_PATH))//'/'//'meshfile.h5'
  ! group name
  write(gname_region, "('reg',i1)") iregion_code

  if (myrank == 0) then
    ! create the file
    if (iregion_code == 1) then
      call h5_create_file(name_database_hdf5)
    endif

    ! open the file
    call h5_open_file(name_database_hdf5)
    ! create the group
    call h5_create_group(gname_region)
    ! open the group
    call h5_open_group(gname_region)

    if (TRANSVERSE_ISOTROPY) then
      ! vpv
      dset_name = 'vpv'
      call h5_create_dataset_gen_in_group(dset_name, (/NGLLX,NGLLY,NGLLZ,sum(offset_nelems(:))/), 4, CUSTOM_REAL)
      ! vph
      dset_name = 'vph'
      call h5_create_dataset_gen_in_group(dset_name, (/NGLLX,NGLLY,NGLLZ,sum(offset_nelems(:))/), 4, CUSTOM_REAL)
      ! vsv
      dset_name = 'vsv'
      call h5_create_dataset_gen_in_group(dset_name, (/NGLLX,NGLLY,NGLLZ,sum(offset_nelems(:))/), 4, CUSTOM_REAL)
      ! vsh
      dset_name = 'vsh'
      call h5_create_dataset_gen_in_group(dset_name, (/NGLLX,NGLLY,NGLLZ,sum(offset_nelems(:))/), 4, CUSTOM_REAL)
      ! rho
      dset_name = 'rho'
      call h5_create_dataset_gen_in_group(dset_name, (/NGLLX,NGLLY,NGLLZ,sum(offset_nelems(:))/), 4, CUSTOM_REAL)
      ! eta
      dset_name = 'eta'
      call h5_create_dataset_gen_in_group(dset_name, (/NGLLX,NGLLY,NGLLZ,sum(offset_nelems(:))/), 4, CUSTOM_REAL)

    else

      ! isotropic model
      ! vp
      dset_name = 'vp'
      call h5_create_dataset_gen_in_group(dset_name, (/NGLLX,NGLLY,NGLLZ,sum(offset_nelems(:))/), 4, CUSTOM_REAL)
      ! vs
      dset_name = 'vs'
      call h5_create_dataset_gen_in_group(dset_name, (/NGLLX,NGLLY,NGLLZ,sum(offset_nelems(:))/), 4, CUSTOM_REAL)
      ! rho
      dset_name = 'rho'
      call h5_create_dataset_gen_in_group(dset_name, (/NGLLX,NGLLY,NGLLZ,sum(offset_nelems(:))/), 4, CUSTOM_REAL)

    endif ! TRANSVERSE_ISOTROPY

    ! anisotropic values
    if (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
      ! Gc_prime
      dset_name = 'Gc_prime'
      call h5_create_dataset_gen_in_group(dset_name, (/NGLLX,NGLLY,NGLLZ,sum(offset_nelems(:))/), 4, CUSTOM_REAL)
      ! Gs_prime
      dset_name = 'Gs_prime'
      call h5_create_dataset_gen_in_group(dset_name, (/NGLLX,NGLLY,NGLLZ,sum(offset_nelems(:))/), 4, CUSTOM_REAL)
      ! mu0
      dset_name = 'mu0'
      call h5_create_dataset_gen_in_group(dset_name, (/NGLLX,NGLLY,NGLLZ,sum(offset_nelems(:))/), 4, CUSTOM_REAL)
    endif

    if (ATTENUATION) then
      ! Qmu
      dset_name = 'Qmu'
      call h5_create_dataset_gen_in_group(dset_name, (/NGLLX,NGLLY,NGLLZ,sum(offset_nelems(:))/), 4, CUSTOM_REAL)
    endif

    if (HETEROGEN_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
      ! dvp
      dset_name = 'dvp'
      call h5_create_dataset_gen_in_group(dset_name, (/NGLLX,NGLLY,NGLLZ,sum(offset_nelems(:))/), 4, CUSTOM_REAL)
    endif

    ! close the group
    call h5_close_group()
    ! close the file
    call h5_close_file()

  endif ! myrank == 0

  ! synchronizes MPI processes
  call synchronize_all()

  !
  ! write the data
  !
  ! open the file
  call h5_open_file_p(name_database_hdf5)
  call h5_open_group(gname_region)

  ! scaling factors to re-dimensionalize units
  scaleval1 = real(sqrt(PI*GRAV*RHOAV)*(R_PLANET/1000.0d0),kind=CUSTOM_REAL)
  scaleval2 = real(RHOAV/1000.0d0,kind=CUSTOM_REAL)

  ! uses temporary array
  allocate(temp_store(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) stop 'Error allocating temp_store array'
  temp_store(:,:,:,:) = 0._CUSTOM_REAL

  if (TRANSVERSE_ISOTROPY) then

    ! vpv
    temp_store(:,:,:,:) = sqrt((kappavstore(:,:,:,:) + 4.0_CUSTOM_REAL * muvstore(:,:,:,:)/3.0_CUSTOM_REAL)/rhostore(:,:,:,:)) &
                          * scaleval1
    dset_name = 'vpv'
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, temp_store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
    ! vph
    temp_store(:,:,:,:) = sqrt((kappahstore(:,:,:,:) + 4.0_CUSTOM_REAL * muhstore(:,:,:,:)/3.0_CUSTOM_REAL)/rhostore(:,:,:,:)) &
                          * scaleval1
    dset_name = 'vph'
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, temp_store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
    ! vsv
    temp_store(:,:,:,:) = sqrt( muvstore(:,:,:,:)/rhostore(:,:,:,:) )*scaleval1
    dset_name = 'vsv'
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, temp_store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
    ! vsh
    temp_store(:,:,:,:) = sqrt( muhstore(:,:,:,:)/rhostore(:,:,:,:) )*scaleval1
    dset_name = 'vsh'
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, temp_store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
    ! rho
    temp_store(:,:,:,:) = rhostore(:,:,:,:) * scaleval2
    dset_name = 'rho'
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, temp_store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
    ! eta
    dset_name = 'eta'
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, eta_anisostore, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)

  else
    ! isotropic model
    ! vp
    temp_store(:,:,:,:) = sqrt((kappavstore(:,:,:,:) + 4.0_CUSTOM_REAL * muvstore(:,:,:,:)/3.0_CUSTOM_REAL)/rhostore(:,:,:,:)) &
                          * scaleval1
    dset_name = 'vp'
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, temp_store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
    ! vs
    temp_store(:,:,:,:) = sqrt( muvstore(:,:,:,:)/rhostore(:,:,:,:) )*scaleval1
    dset_name = 'vs'
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, temp_store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
    ! rho
    temp_store(:,:,:,:) = rhostore(:,:,:,:) * scaleval2
    dset_name = 'rho'
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, temp_store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)

  endif ! TRANSVERSE_ISOTROPY

  ! anisotropic values
  if (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
    ! the scale of GPa--[g/cm^3][(km/s)^2]
    scaleval = real(sqrt(PI*GRAV*RHOAV),kind=CUSTOM_REAL)
    scale_GPa = real((RHOAV/1000.d0)*((R_PLANET*scaleval/1000.d0)**2),kind=CUSTOM_REAL)

    ! Gc_prime
    dset_name = 'Gc_prime'
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, Gc_prime_store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
    ! Gs_prime
    dset_name = 'Gs_prime'
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, Gs_prime_store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
    ! mu0
    temp_store(:,:,:,:) = mu0store(:,:,:,:) * scale_GPa
    dset_name = 'mu0'
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, temp_store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)

  endif

  if (ATTENUATION) then
    if (ATTENUATION_3D .or. ATTENUATION_3D_BERKELEY .or. ATTENUATION_1D_WITH_3D_STORAGE) then
      temp_store(:,:,:,:) = Qmu_store(:,:,:,:)
    else
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

    ! Qmu
    dset_name = 'Qmu'
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, temp_store, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)

  endif ! ATTENUATION

  deallocate(temp_store)

  if (HETEROGEN_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
    ! dvp
    dset_name = 'dvp'
    call h5_write_dataset_collect_hyperslab_in_group(dset_name, dvpstore, (/0,0,0,sum(offset_nelems(0:myrank-1))/), H5_COL)
  endif

  ! close the group
  call h5_close_group()
  ! close the file
  call h5_close_file()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     model region ',iregion_code,'done'
    call flush_IMAIN()
  endif

#else

  print *, 'Error: HDF5 support not enabled in this build'
  stop

#endif

  end subroutine save_model_meshfiles_hdf5
