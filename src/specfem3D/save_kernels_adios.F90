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


!-------------------------------------------------------------------------------
!> \file save_kernels_adios.f90
!! \brief Save kernels arrays to file with the help of the ADIOS library.
!! \author MPBL
!-------------------------------------------------------------------------------

#include "config.fh"

! note: to write out the kernel in ADIOS format, we start defining the group, all variables and open the adios file.
!       the write routines are then called if the corresponding flags have been set.
!       at the very end, we close the ADIOS kernel file.
!
!       for adios2, this implies that the memory pointers to the kernel arrays should not change if we use write_adios_** routines
!       with a deferred mode. unfortunately, that's not the case as we create a lot of temporary arrays in subroutines to
!       construct different kernels for different parameterizations based on a few "primary" kernels.
!
!       to circumvent, we would need either to allocate all necessary kernels and free them only after the closing the file.
!       or, we open/close/append to the file in each subroutine. however, ADIOS will increase a step counter each time
!       when we append to a file which might lead to issues when reading back the kernels.
!
!       using adios2, we do have an explicity perform_write routine to call to solve this. unfortunately, in adios1 there
!       is no such routine, unless we close the path of the group. this will erase the path info for all variables in the group
!       which we unfortunately use to assign variables like
!         my_kernel_kl/local_dim, my_kernel_kl/offset, mykernel_kl/global_dim and mykernel_kl/array
!
!       unless we drop all these additional infos (offset,..) per array, it will corrupt the writing.
!       however, it looks like adios1 is copying the arrays in the adios_write call and we're fine when closing the file.
!
!       so for now, we only call the write_adios_perform() routine for adios2 cases and only at the very end for adios1.
!       todo: check this behavior in future versions...



!==============================================================================
!> Define all the kernels that will be written to the ADIOS file.
!!
!! \note Everything is define in this single function, even the group size.
!!       It is the reason why this function require only an handle on an ADIOS
!!       file as an argument.
  subroutine define_kernel_adios_variables()

  use specfem_par ! Just for dimensions. No need of arrays for now.
  use specfem_par_crustmantle
  use specfem_par_outercore
  use specfem_par_innercore
  use specfem_par_noise

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! local Variables
  character(len=MAX_STRING_LEN) :: outputname, group_name
  integer(kind=8) :: group_size_inc
  integer(kind=8) :: local_dim

  ! Type inference for define_adios_global_array1D. Avoid additional args.
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: dummy_real4d
  integer :: ier

  ! user output
  if (myrank == 0) then
#if defined(USE_ADIOS)
    write(IMAIN,*) '  saving kernels in ADIOS 1 file format'
#elif defined(USE_ADIOS2)
    write(IMAIN,*) '  saving kernels in ADIOS 2 file format'
#endif
    call flush_IMAIN()
  endif

  outputname = get_adios_filename(trim(OUTPUT_FILES)//"/kernels")

  group_name = "SPECFEM3D_GLOBE_KERNELS"

  ! set the adios group size to 0 before incremented by calls to helpers functions.
  group_size_inc = 0
  call init_adios_group(myadios_group,group_name)

  ! scalar for checking
  call define_adios_scalar(myadios_group, group_size_inc, '', "NSPEC", NSPEC_CRUST_MANTLE_ADJOINT)
  call define_adios_scalar(myadios_group, group_size_inc, '', "reg1/nspec", NSPEC_CRUST_MANTLE_ADJOINT)

  if (SIMULATION_TYPE == 3) then
    ! dummy for definitions
    allocate(dummy_real4d(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating dummy array'
    dummy_real4d(:,:,:,:) = 0.0

    ! crust mantle
    if (ANISOTROPIC_KL) then

      ! anisotropic kernels
      local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT

      ! outputs transverse isotropic kernels only
      if (SAVE_TRANSVERSE_KL_ONLY) then
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                         "alphav_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                         "alphah_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                         "betav_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                         "betah_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                         "eta_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                         "rho_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                         "bulk_c_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                         "bulk_betav_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                         "bulk_betah_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                         "alpha_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                         STRINGIFY_VAR(beta_kl_crust_mantle))
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                         "bulk_beta_kl_crust_mantle", dummy_real4d)

      else if (SAVE_AZIMUTHAL_ANISO_KL_ONLY) then
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim,'', &
                                         "betav_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim,'', &
                                         "betah_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim,'', &
                                         "eta_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim,'', &
                                         "rho_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim,'', &
                                         "bulk_c_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim,'', &
                                         "bulk_betav_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim,'', &
                                         "bulk_betah_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim,'', &
                                         "alpha_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim,'', &
                                         STRINGIFY_VAR(beta_kl_crust_mantle))
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim,'', &
                                         "bulk_beta_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim,'', &
                                         "Gc_prime_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim,'', &
                                         "Gs_prime_kl_crust_mantle", dummy_real4d)

      else
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                         STRINGIFY_VAR(rho_kl_crust_mantle))
        local_dim = 21 * NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                         STRINGIFY_VAR(cijkl_kl_crust_mantle))
      endif

    else

      ! isotropic kernels
      local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                       "rhonotprime_kl_crust_mantle", dummy_real4d)
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                       "kappa_kl_crust_mantle", dummy_real4d)
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                       "mu_kl_crust_mantle", dummy_real4d)
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                       STRINGIFY_VAR(rho_kl_crust_mantle))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                       STRINGIFY_VAR(alpha_kl_crust_mantle))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                       STRINGIFY_VAR(beta_kl_crust_mantle))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                       "bulk_c_kl_crust_mantle", dummy_real4d)
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                       "bulk_beta_kl_crust_mantle", dummy_real4d)
    endif

    ! noise strength kernel
    if (NOISE_TOMOGRAPHY == 3) then
      local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', STRINGIFY_VAR(sigma_kl_crust_mantle))
    endif

    ! outer core
    if (SAVE_KERNELS_OC) then
      local_dim = NSPEC_OUTER_CORE * NGLLX * NGLLY * NGLLZ
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', STRINGIFY_VAR(rho_kl_outer_core))
      call define_adios_global_array1D(myadios_group, group_size_inc,local_dim, '', STRINGIFY_VAR(alpha_kl_outer_core))
      if (deviatoric_outercore) then
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', STRINGIFY_VAR(beta_kl_outer_core))
      endif
    endif

    ! inner core
    if (SAVE_KERNELS_IC) then
      local_dim = NSPEC_INNER_CORE * NGLLX * NGLLY * NGLLZ
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', STRINGIFY_VAR(rho_kl_inner_core))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', STRINGIFY_VAR(alpha_kl_inner_core))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', STRINGIFY_VAR(beta_kl_inner_core))
    endif

    ! boundary kernel
    if (SAVE_KERNELS_BOUNDARY) then
      !call save_kernels_boundary_kl()
      if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
        local_dim = NSPEC2D_MOHO * NGLLX * NGLLY * NDIM
        call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', STRINGIFY_VAR(moho_kl))
      endif
      local_dim = NSPEC2D_400 * NGLLX * NGLLY
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', STRINGIFY_VAR(d400_kl))

      local_dim = NSPEC2D_670 * NGLLX * NGLLY
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', STRINGIFY_VAR(d670_kl))

      local_dim = NSPEC2D_CMB * NGLLX * NGLLY
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', STRINGIFY_VAR(cmb_kl))

      local_dim = NSPEC2D_ICB * NGLLX * NGLLY
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', STRINGIFY_VAR(icb_kl))
    endif

    ! approximate Hessian
    if (APPROXIMATE_HESS_KL) then
      !call save_kernels_Hessian()
      local_dim = NSPEC_CRUST_MANTLE_ADJOINT* NGLLX * NGLLY * NGLLZ
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                       STRINGIFY_VAR(hess_kl_crust_mantle))

      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                       STRINGIFY_VAR(hess_rho_kl_crust_mantle))

      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                       STRINGIFY_VAR(hess_kappa_kl_crust_mantle))

      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', &
                                       STRINGIFY_VAR(hess_mu_kl_crust_mantle))
    endif

    deallocate(dummy_real4d)

  endif

  ! save source derivatives for adjoint simulations
  if (SIMULATION_TYPE == 2 .and. nrec_local > 0) then
    local_dim = 3 * 3 * nrec_local
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', STRINGIFY_VAR(moment_der))

    local_dim = 3 * nrec_local
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', STRINGIFY_VAR(sloc_der))

    local_dim = nrec_local
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', STRINGIFY_VAR(stshift_der))

    local_dim = nrec_local
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, '', STRINGIFY_VAR(shdur_der))
  endif

  ! Open the handle to file containing all the ADIOS variables previously defined
  ! opens file for writing
  call open_file_adios_write(myadios_file,myadios_group,outputname,group_name)
  call set_adios_group_size(myadios_file,group_size_inc)

  ! scalar for checking
  ! writes nspec
  call write_adios_scalar(myadios_file,myadios_group,"NSPEC",NSPEC_CRUST_MANTLE_ADJOINT)
  call write_adios_scalar(myadios_file,myadios_group,"reg1/nspec",NSPEC_CRUST_MANTLE_ADJOINT)

  end subroutine define_kernel_adios_variables


!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the crust mantle.
  subroutine write_kernels_cm_ani_adios(alphav_kl_crust_mantle,alphah_kl_crust_mantle, &
                                        betav_kl_crust_mantle,betah_kl_crust_mantle, &
                                        eta_kl_crust_mantle, &
                                        bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle, &
                                        bulk_betav_kl_crust_mantle,bulk_betah_kl_crust_mantle, &
                                        Gc_prime_kl_crust_mantle,Gs_prime_kl_crust_mantle)

  use specfem_par
  use specfem_par_crustmantle

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! input Parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
      alphav_kl_crust_mantle,alphah_kl_crust_mantle, &
      betav_kl_crust_mantle,betah_kl_crust_mantle, &
      eta_kl_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
      bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle, &
      bulk_betav_kl_crust_mantle,bulk_betah_kl_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
      Gc_prime_kl_crust_mantle, Gs_prime_kl_crust_mantle

  ! Variables
  integer(kind=8) :: local_dim

  ! checks if anything to do
  if (.not. ANISOTROPIC_KL) return

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT

  ! For anisotropic kernels
  ! outputs transverse isotropic kernels only
  if (SAVE_TRANSVERSE_KL_ONLY) then
    ! transverse isotropic kernels
    ! (alpha_v, alpha_h, beta_v, beta_h, eta, rho ) parameterization
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     STRINGIFY_VAR(alphav_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     STRINGIFY_VAR(alphah_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     STRINGIFY_VAR(betav_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     STRINGIFY_VAR(betah_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     STRINGIFY_VAR(eta_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     STRINGIFY_VAR(rho_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     STRINGIFY_VAR(bulk_c_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     STRINGIFY_VAR(bulk_betav_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     STRINGIFY_VAR(bulk_betah_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     STRINGIFY_VAR(alpha_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     STRINGIFY_VAR(beta_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     STRINGIFY_VAR(bulk_beta_kl_crust_mantle))
  else if (SAVE_AZIMUTHAL_ANISO_KL_ONLY) then
    ! kernels for inversions involving azimuthal anisotropy
    ! (bulk_c, beta_v, beta_h, eta, Gc', Gs', rho ) parameterization
    ! note: Gc' & Gs' are the normalized Gc & Gs kernels
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios,local_dim, &
                                     STRINGIFY_VAR(betav_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios,local_dim, &
                                     STRINGIFY_VAR(betah_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios,local_dim, &
                                     STRINGIFY_VAR(eta_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios,local_dim, &
                                     STRINGIFY_VAR(rho_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios,local_dim, &
                                     STRINGIFY_VAR(bulk_c_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios,local_dim, &
                                     STRINGIFY_VAR(bulk_betav_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios,local_dim, &
                                     STRINGIFY_VAR(bulk_betah_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios,local_dim, &
                                     STRINGIFY_VAR(alpha_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios,local_dim, &
                                     STRINGIFY_VAR(beta_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios,local_dim, &
                                     STRINGIFY_VAR(bulk_beta_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios,local_dim, &
                                     STRINGIFY_VAR(Gc_prime_kl_crust_mantle))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios,local_dim, &
                                     STRINGIFY_VAR(Gs_prime_kl_crust_mantle))
  else
    ! note: the C_ij and density kernels are not for relative perturbations
    !       (delta ln( m_i) = delta m_i / m_i),
    !       but absolute perturbations (delta m_i = m_i - m_0)
    !
    ! note: the write_adios_** calls might use deferred mode and only synchronize the array when perform/close/end_step is done.
    !       using a temporary kernel here might corrupt the memory pointer when exiting this subroutine and loosing the scope.
    !       we thus will put a minus sign like - rho_kl_.. in the top routine and only pass the array pointer
    !       to this subroutine.
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     "rho_kl_crust_mantle", rho_kl_crust_mantle)

    local_dim = 21 * NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     "cijkl_kl_crust_mantle", cijkl_kl_crust_mantle)
  endif

  ! sync adios2 writes before loosing the temporary scope of the arrays
  if (is_adios_version2) call write_adios_perform(myadios_file)

  end subroutine write_kernels_cm_ani_adios

!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the crust mantle.
  subroutine write_kernels_cm_iso_adios(mu_kl_crust_mantle, kappa_kl_crust_mantle, rhonotprime_kl_crust_mantle, &
                                        bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle)

  use specfem_par
  use specfem_par_crustmantle

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! Parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
      mu_kl_crust_mantle, kappa_kl_crust_mantle, rhonotprime_kl_crust_mantle, &
      bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle

  ! Variables
  integer(kind=8) :: local_dim

  ! checks if anything to do
  if (ANISOTROPIC_KL) return

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT

  ! isotropic kernels
  ! primary kernels: (rho,kappa,mu) parameterization
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(rhonotprime_kl_crust_mantle))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(kappa_kl_crust_mantle))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(mu_kl_crust_mantle))

  ! (rho, alpha, beta ) parameterization
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(rho_kl_crust_mantle))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(alpha_kl_crust_mantle))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(beta_kl_crust_mantle))

  ! (rho, bulk, beta ) parameterization, K_rho same as above
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(bulk_c_kl_crust_mantle))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(bulk_beta_kl_crust_mantle))

  ! sync adios2 writes before loosing the temporary scope of the arrays
  if (is_adios_version2) call write_adios_perform(myadios_file)

  end subroutine write_kernels_cm_iso_adios


!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the outer core.
  subroutine write_kernels_oc_adios()

  use specfem_par
  use specfem_par_outercore

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! Variables
  integer(kind=8) :: local_dim

  ! checks if anything to do
  if (.not. SAVE_KERNELS_OC) return

  local_dim = NSPEC_OUTER_CORE * NGLLX* NGLLY * NGLLZ

  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(rho_kl_outer_core))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(alpha_kl_outer_core))

  !deviatoric kernel check
  if (deviatoric_outercore) then
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     STRINGIFY_VAR(beta_kl_outer_core))
  endif

  ! sync adios2 writes before loosing the temporary scope of the arrays
  if (is_adios_version2) call write_adios_perform(myadios_file)

  end subroutine write_kernels_oc_adios


!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the inner core.
  subroutine write_kernels_ic_adios()

  use specfem_par
  use specfem_par_innercore

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! Variables
  integer(kind=8) :: local_dim

  ! checks if anything to do
  if (.not. SAVE_KERNELS_IC) return

  local_dim = NSPEC_INNER_CORE * NGLLX * NGLLY * NGLLZ

  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(rho_kl_inner_core))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(alpha_kl_inner_core))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(beta_kl_inner_core))

  ! sync adios2 writes before loosing the temporary scope of the arrays
  if (is_adios_version2) call write_adios_perform(myadios_file)

  end subroutine write_kernels_ic_adios


!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the boundaries.
  subroutine write_kernels_boundary_kl_adios()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! Variables
  integer(kind=8) :: local_dim

  ! checks if anything to do
  if (.not. SAVE_KERNELS_BOUNDARY) return

  if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
    local_dim = NSPEC2D_MOHO * NGLLX * NGLLY * NDIM

    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, STRINGIFY_VAR(moho_kl))
  endif

  local_dim = NSPEC2D_400 * NGLLX * NGLLY
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, STRINGIFY_VAR(d400_kl))

  local_dim = NSPEC2D_670 * NGLLX * NGLLY
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, STRINGIFY_VAR(d670_kl))

  local_dim = NSPEC2D_CMB * NGLLX * NGLLY
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, STRINGIFY_VAR(cmb_kl))

  local_dim = NSPEC2D_ICB * NGLLX * NGLLY
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, STRINGIFY_VAR(icb_kl))

  ! sync adios2 writes before loosing the temporary scope of the arrays
  if (is_adios_version2) call write_adios_perform(myadios_file)

  end subroutine write_kernels_boundary_kl_adios


!==============================================================================
!> Schedule writes for the source derivatives (moment tensors and source
!! locations.
!!
!! \note Not to dump one value at a time (as in the non ADIOS version) data are
!!       scaled and oriented but are not written in the same order than in the
!!       non ADIOS version.
!!       (see save_kernels_source_derivatives in save_kernels.f90)
  subroutine write_kernels_source_derivatives_adios()

  use specfem_par

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! Variables
  ! We do not want to change moment_der and sloc as it might introduce future
  ! concerns if we want to use them after.
  ! No use of modified moment_der values since it implies to allocate those
  ! arrays in save_kernels() and to carry them along the way. It might be better
  ! to transform these arrays in the post processing phase.
  !real(kind=CUSTOM_REAL), dimension(3,3,nrec_local) :: moment_der_tmp
  !real(kind=CUSTOM_REAL), dimension(3,nrec_local) :: sloc_der_tmp
  integer(kind=8) :: local_dim

  !moment_der_tmp(:, :, :) = moment_der(:, :, :) * 1e-7
  !moment_der_tmp(1, 3, :) = -2 * moment_der(1, 3, :)
  !moment_der_tmp(2, 3, :) =  2 * moment_der(2, 3, :)
  !moment_der_tmp(1, 2, :) = -2 * moment_der(1, 2, :)
  !moment_der_tmp(3, 1, :) = moment_der(1, 3, :)
  !moment_der_tmp(3, 2, :) = moment_der(2, 3, :)
  !moment_der_tmp(2, 1, :) = moment_der(1, 2, :)
  !sloc_der_tmp(:, :) = - sloc_der(:, :)
  !sloc_der_tmp(3, :) = - sloc_der(3, :)


  local_dim = 3 * 3 * nrec_local
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, STRINGIFY_VAR(moment_der))

  local_dim = 3 * nrec_local
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, STRINGIFY_VAR(sloc_der))

  local_dim = nrec_local
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, STRINGIFY_VAR(stshift_der))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, STRINGIFY_VAR(shdur_der))

  ! sync adios2 writes before loosing the temporary scope of the arrays
  if (is_adios_version2) call write_adios_perform(myadios_file)

  end subroutine write_kernels_source_derivatives_adios


!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the Hessian.
  subroutine write_kernels_Hessian_adios()

  use specfem_par
  use specfem_par_crustmantle

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! Variables
  integer(kind=8) :: local_dim

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT

  ! stores into file
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(hess_kl_crust_mantle))

  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(hess_rho_kl_crust_mantle))

  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(hess_kappa_kl_crust_mantle))

  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(hess_mu_kl_crust_mantle))

  ! sync adios2 writes before loosing the temporary scope of the arrays
  if (is_adios_version2) call write_adios_perform(myadios_file)

  end subroutine write_kernels_Hessian_adios

!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the noise strength kernel.
  subroutine write_kernels_strength_noise_adios()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_noise

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! Variables
  integer(kind=8) :: local_dim

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT

  ! stores into file
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(sigma_kl_crust_mantle))

  ! sync adios2 writes before loosing the temporary scope of the arrays
  if (is_adios_version2) call write_adios_perform(myadios_file)

  end subroutine write_kernels_strength_noise_adios


!==============================================================================
!> closes kernel ADIOS file.
  subroutine close_kernel_adios_file()

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! finalize method and close file
  call close_file_adios(myadios_file)

  end subroutine close_kernel_adios_file
