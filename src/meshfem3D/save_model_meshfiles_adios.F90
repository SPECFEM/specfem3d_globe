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

#include "config.fh"


!===============================================================================
!> \brief Save the meshfiles for visualization that will be used by the solver in an ADIOS format.
!!
!! \param myrank The MPI rank of the current process.
!! \param reg_name Output file prefix with the name of the region included
!! \param nspec Number of GLL points per spectral elements
  subroutine save_model_meshfiles_adios()

  use constants
  use shared_parameters, only: R_PLANET,RHOAV

  use meshfem_par, only: &
    LOCAL_PATH,nspec,nglob,iregion_code

  use meshfem_models_par, only: &
    TRANSVERSE_ISOTROPY,ATTENUATION, &
    ATTENUATION_3D,ATTENUATION_1D_WITH_3D_STORAGE, &
    ANISOTROPIC_3D_MANTLE

  use regions_mesh_par2, only: &
    rhostore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
    Qmu_store,Gc_prime_store,Gs_prime_store,mu0store

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! local parameters
  integer :: i,j,k,ispec,ier
  real(kind=CUSTOM_REAL) :: scaleval1,scaleval2,scaleval,scale_GPa
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: temp_store_vp,temp_store_vs,temp_store_rho,temp_store_rho_inv
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: temp_store_vpv,temp_store_vph,temp_store_vsv,temp_store_vsh
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: temp_store_mu0
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: temp_store_Qmu

  ! local parameters
  character(len=MAX_STRING_LEN) :: outputname, group_name
  integer(kind=8) :: local_dim
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  character(len=128) :: region_name, region_name_scalar
  !--- Save the number of region written. Open the file in "w" mode if 0, else
  !    in "a"  mode
  integer, save :: num_regions_written = 0

  ! user output
  if (myrank == 0) then
#if defined(USE_ADIOS)
    write(IMAIN,*) '    model    in ADIOS 1 file format'
#elif defined(USE_ADIOS2)
    write(IMAIN,*) '    model    in ADIOS 2 file format'
#endif
    call flush_IMAIN()
  endif

  ! scaling factors to re-dimensionalize units
  scaleval1 = sngl( sqrt(PI*GRAV*RHOAV)*(R_PLANET/1000.0d0) )
  scaleval2 = sngl( RHOAV/1000.0d0 )

! note: the following uses temporary arrays for array expressions like sqrt( (kappavstore+..)).
!       since the write_adios_** calls might be in deferred mode, these temporary arrays should be valid
!       until a perform/close/end_step call is done.
!
!       as a work-around, we will explicitly allocate temporary arrays and deallocate them after the file close.
  allocate(temp_store_vp(NGLLX,NGLLY,NGLLZ,nspec), &
           temp_store_vs(NGLLX,NGLLY,NGLLZ,nspec), &
           temp_store_rho(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) stop 'Error allocating temp vp,.. arrays'
  temp_store_vp(:,:,:,:) = 0.0_CUSTOM_REAL
  temp_store_vs(:,:,:,:) = 0.0_CUSTOM_REAL
  temp_store_rho(:,:,:,:) = 0.0_CUSTOM_REAL

  allocate(temp_store_rho_inv(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) stop 'Error allocating temp rho_inv array'

  ! this might have issues when rho is zero in fictitious inner core elements:
  !temp_store_rho(:,:,:,:) = rhostore(:,:,:,:) * scaleval2
  !temp_store_vp(:,:,:,:) = sqrt( (kappavstore + 4.0_CUSTOM_REAL * muvstore/3.00_CUSTOM_REAL)/rhostore ) * scaleval1
  !temp_store_vs(:,:,:,:) = sqrt( muvstore/rhostore ) * scaleval1
  !
  ! takes inverse of rho, avoiding zero values in fictitious elements:
  temp_store_rho_inv(:,:,:,:) = rhostore(:,:,:,:)
  where(temp_store_rho_inv(:,:,:,:) <= 0.0_CUSTOM_REAL) temp_store_rho_inv = 1.0_CUSTOM_REAL
  temp_store_rho_inv = 1.0_CUSTOM_REAL / temp_store_rho_inv

  ! rho: for storing, we take the original rho and dimensionalize it
  temp_store_rho(:,:,:,:) = rhostore(:,:,:,:) * scaleval2
  ! vp
  temp_store_vp(:,:,:,:) = sqrt( (kappavstore + 4.0_CUSTOM_REAL * muvstore/3.0_CUSTOM_REAL) * temp_store_rho_inv ) * scaleval1
  ! vs
  temp_store_vs(:,:,:,:) = sqrt( muvstore * temp_store_rho_inv ) * scaleval1

  ! isotropic model
  write(region_name,"('reg',i1, '/')") iregion_code
  write(region_name_scalar,"('reg',i1)") iregion_code
  write(group_name,"('SPECFEM3D_GLOBE_MODEL_reg',i1)") iregion_code

!daniel
! note: on Mac OsX with OpenMPI 3.1, an error can occur at the end when finalizing MPI:
!
!--------------------------------------------------------------------------
!A system call failed during shared memory initialization that should
!not have.  It is likely that your MPI job will now either abort or
!experience performance degradation.
!
!  Local host:  mac.home
!  System call: unlink(2) /var/folders/0b/8r9lhyv48v58lf006s7hkfmr0000gn/T//ompi.mac.501/pid.39081/1/vader_segment.mac.e8f40001.0
!  Error:       No such file or directory (errno 2)
!--------------------------------------------------------------------------
!
! this seems to be related to this routine somehow and the crust/mantle region.
!
! a work-around could be to call mpirun with:
! > mpirun --mca btl_vader_backing_directory /tmp -np 4 ./bin/xmeshfem3D
!
! or either in your shell, use:
! > export OMPI_MCA_btl=self,tcp
! or when calling mpirun:
! > OMPI_MCA_btl=self,tcp mpirun -np 4 ./bin/xmeshfem3D
!
!
! as software work-around, we re-initialize adios.
!
! todo: this will need to be re-evaluated in future as it might be fixed in future versions.
  if (is_adios_version2) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '    re-initializes adios2 for meshfile model output'
      call flush_IMAIN()
    endif
    call synchronize_all()
    ! frees adios main object
    call finalize_adios()
    ! re-initializes
    call initialize_adios()
  endif

  ! initializes i/o group
  call init_adios_group(myadios_val_group,group_name)

  ! save nspec and nglob, to be used in combine_paraview_data
  group_size_inc = 0
  call define_adios_scalar (myadios_val_group, group_size_inc, &
                            region_name_scalar, STRINGIFY_VAR(nspec))
  call define_adios_scalar (myadios_val_group, group_size_inc, &
                            region_name_scalar, STRINGIFY_VAR(nglob))

  ! array sizes
  local_dim = NGLLX * NGLLY * NGLLZ * nspec

  ! checks size
  if (size(kappavstore) /= local_dim) then
    print *,'Error: size kappavstore ',size(kappavstore), ' should be ',local_dim
    call exit_mpi(myrank,'Error size kappavstore for storing meshfiles')
  endif

  !--- Define ADIOS variables -----------------------------
  ! rho
  call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "rho", temp_store_rho)
  ! vp
  call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "vp", temp_store_vp)
  ! vs (will store it even for the outer core, although it should just be zero there)
  call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "vs", temp_store_vs)

  ! transverse isotropic model
  if (TRANSVERSE_ISOTROPY) then
    ! vpv
    ! (using temp_store_vp as dummy array here for definition, but with correct array size and type)
    call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "vpv", temp_store_vp) !dummy
    ! vph
    call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "vph", temp_store_vp) !dummy
    ! vsv
    call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "vsv", temp_store_vp) !dummy
    ! vsh
    call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "vsh", temp_store_vp) !dummy
    ! eta
    call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "eta", temp_store_vp) !dummy
  endif

  ! anisotropic values
  if (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
    ! Gc_prime
    call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "Gc_prime", temp_store_vp) !dummy
    ! Gs_prime
    call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "Gs_prime", temp_store_vp) !dummy
    ! mu0
    call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "mu0", temp_store_vp) !dummy
  endif

  if (ATTENUATION) then
    ! Qmu
    call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "qmu", temp_store_vp) !dummy
  endif

  !--- Open an ADIOS handler to the restart file. ---------
  outputname = get_adios_filename(trim(LOCAL_PATH) // "/model_gll")

  ! user output
  if (myrank == 0) write(IMAIN,*) '    saving model meshfile arrays in ADIOS file: ',trim(outputname)

  if (num_regions_written == 0) then
    ! opens file for writing
    call open_file_adios_write(myadios_val_file,myadios_val_group,outputname,group_name)
  else
    ! opens file for writing in append mode
    call open_file_adios_write_append(myadios_val_file,myadios_val_group,outputname,group_name)
  endif

  call set_adios_group_size(myadios_val_file,group_size_inc)

  ! save nspec and nglob, to be used in combine_paraview_data
  call write_adios_scalar(myadios_val_file,myadios_val_group,trim(region_name) // "nspec",nspec)
  call write_adios_scalar(myadios_val_file,myadios_val_group,trim(region_name) // "nglob",nglob)

  !--- Schedule writes for the previously defined ADIOS variables
  ! rho
  call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // "rho", temp_store_rho)
  ! vp
  call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // "vp", temp_store_vp)
  ! vs
  call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // "vs", temp_store_vs)

  ! transverse isotropic model
  if (TRANSVERSE_ISOTROPY) then
    allocate(temp_store_vpv(NGLLX,NGLLY,NGLLZ,nspec), &
             temp_store_vph(NGLLX,NGLLY,NGLLZ,nspec), &
             temp_store_vsv(NGLLX,NGLLY,NGLLZ,nspec), &
             temp_store_vsh(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
    if (ier /= 0) stop 'Error allocating temp vpv,.. arrays'

    ! original:
    !temp_store_vpv(:,:,:,:) = sqrt( (kappavstore+4.0_CUSTOM_REAL*muvstore/3.0_CUSTOM_REAL)/rhostore ) * scaleval1
    !temp_store_vph(:,:,:,:) = sqrt( (kappahstore+4.0_CUSTOM_REAL*muhstore/3.0_CUSTOM_REAL)/rhostore ) * scaleval1
    !temp_store_vsv(:,:,:,:) = sqrt( muvstore/rhostore ) * scaleval1
    !temp_store_vsh(:,:,:,:) = sqrt( muhstore/rhostore ) * scaleval1
    ! using safer inverse rho array:
    temp_store_vpv(:,:,:,:) = sqrt( (kappavstore + 4.0_CUSTOM_REAL*muvstore/3.0_CUSTOM_REAL) * temp_store_rho_inv ) * scaleval1
    temp_store_vph(:,:,:,:) = sqrt( (kappahstore + 4.0_CUSTOM_REAL*muhstore/3.0_CUSTOM_REAL) * temp_store_rho_inv ) * scaleval1
    temp_store_vsv(:,:,:,:) = sqrt( muvstore * temp_store_rho_inv ) * scaleval1
    temp_store_vsh(:,:,:,:) = sqrt( muhstore * temp_store_rho_inv ) * scaleval1

    ! vpv
    call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // "vpv", temp_store_vpv)
    ! vph
    call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // "vph", temp_store_vph)
    ! vsv
    call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // "vsv", temp_store_vsv)
    ! vsh
    call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // "vsh", temp_store_vsh)
    ! eta
    call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // "eta", eta_anisostore)
  else
    ! dummy
    allocate(temp_store_vpv(1,1,1,1), &
             temp_store_vph(1,1,1,1), &
             temp_store_vsv(1,1,1,1), &
             temp_store_vsh(1,1,1,1))
  endif ! TRANSVERSE_ISOTROPY

  ! anisotropic values
  if (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
    ! the scale of GPa--[g/cm^3][(km/s)^2]
    scaleval = real(sqrt(PI*GRAV*RHOAV),kind=CUSTOM_REAL)
    scale_GPa = real((RHOAV/1000.d0)*((R_PLANET*scaleval/1000.d0)**2),kind=CUSTOM_REAL)

    allocate(temp_store_mu0(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
    if (ier /= 0) stop 'Error allocating temp mu0 array'

    temp_store_mu0(:,:,:,:) = mu0store(:,:,:,:) * scale_GPa

    ! Gc_prime
    call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // "Gc_prime", Gc_prime_store)
    ! Gs_prime
    call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // "Gs_prime", Gs_prime_store)
    ! mu0
    call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // "mu0", temp_store_mu0)
  else
    ! dummy
    allocate(temp_store_mu0(1,1,1,1))
  endif

  ! shear attenuation
  if (ATTENUATION) then
    ! Qmu arrays
    ! saves Qmu_store to full CUSTOM_REAL array
    ! uses temporary array

    ! note: write_adios_* calls could be in deferred mode and not synchronized.
    !       thus, temporary arrays should be valid and unmodified until a perform/close/end_step call is done.
    !
    !       we will thus keep the allocated array until after the file close.
    allocate(temp_store_Qmu(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
    if (ier /= 0) stop 'Error allocating temp qmu array'

    if (ATTENUATION_3D .or. ATTENUATION_1D_WITH_3D_STORAGE) then
      ! attenuation arrays are fully 3D
      temp_store_Qmu(:,:,:,:) = Qmu_store(:,:,:,:)
    else
      ! fills full attenuation array dimensions: Q_mustore(1,1,1,nspec)
      do ispec = 1,nspec
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              temp_store_Qmu(i,j,k,ispec) = Qmu_store(1,1,1,ispec)
            enddo
          enddo
        enddo
      enddo
    endif
    ! Qmu
    call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // "qmu", temp_store_Qmu)
  else
    ! dummy
    allocate(temp_store_Qmu(1,1,1,1))
  endif ! ATTENUATION

  !--- Reset the path to zero and perform the actual write to disk
  call write_adios_perform(myadios_val_file)
  ! flushes all engines related to this io group (not really required, but used to make sure i/o has all written out)
  call flush_adios_group_all(myadios_val_group)
  ! closes file
  call close_file_adios(myadios_val_file)

  !---------------------------------------------------------
  !--- dvpstore arrays ------------------------------------------
  !---------------------------------------------------------
  !obsolete
  !if (HETEROGEN_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
  !  call model_heterogen_mantle_output_dvp_adios(prname)
  !endif

  num_regions_written = num_regions_written + 1

  ! adios should be done with writing memory out.
  call synchronize_all()

  ! frees temporary memory
  deallocate(temp_store_vp,temp_store_vs,temp_store_rho)
  deallocate(temp_store_rho_inv)
  deallocate(temp_store_vpv,temp_store_vph,temp_store_vsv,temp_store_vsh)
  deallocate(temp_store_mu0)
  deallocate(temp_store_Qmu)

  end subroutine save_model_meshfiles_adios
