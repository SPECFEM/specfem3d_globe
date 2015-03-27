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


!-------------------------------------------------------------------------------
!> \file save_kernels_adios.f90
!! \brief Save kernels arrays to file with the help of the ADIOS library.
!! \author MPBL
!-------------------------------------------------------------------------------

#include "config.fh"

!==============================================================================
!> Perform the actual write of all the kernels variables to file.
!! \param[IN] adios_handle The handle pointing on the open ADIOS file intended
!!                         to store kernels.
!!
!! \note Obviously this is a general routine that should be extracted and used
!!       everywhere as the 'adios_handle' argument can be used for any kind of
!!       ADIOS file.
!!       The only reason such a routine is defined is to avoid using
!!       ADIOS modules in non ADIOS file, in case the ADIOS library is not
!!       available on the system.
subroutine perform_write_adios_kernels(adios_handle)

  use adios_write_mod

  implicit none

  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  ! Variables
  integer :: adios_err

  call adios_close(adios_handle, adios_err)

end subroutine perform_write_adios_kernels

!==============================================================================
!> Define all the kernels that will be written to the ADIOS file.
!! \param[INOUT] adios_handle The handle pointing on the open ADIOS file
!!                            intended to store kernels data.
!! \note Everything is define in this single function, even the group size.
!!       It is the reason why this function require only an handle on an ADIOS
!!       file as an argument.
subroutine define_kernel_adios_variables(adios_handle)

  use adios_write_mod

  use specfem_par ! Just for dimensions. No need of arrays for now.
  use specfem_par_crustmantle
  use specfem_par_outercore
  use specfem_par_innercore

  use adios_helpers_mod

  implicit none

  ! Parameters
  integer(kind=8), intent(INOUT) :: adios_handle
  ! Variables
  character(len=MAX_STRING_LEN) :: outputname, group_name
  integer(kind=8) :: adios_group, group_size_inc, adios_totalsize
  integer :: local_dim, comm, adios_err

  ! Type inference for define_adios_global_array1D. Avoid additional args.
  real(kind=CUSTOM_REAL), dimension(1,1,1,1) :: dummy_real4d

  outputname = trim(OUTPUT_FILES)//"/kernels.bp"
  group_name = "SPECFEM3D_GLOBE_KERNELS"

  call world_duplicate(comm)

  group_size_inc = 0
  call adios_declare_group(adios_group, group_name, "", 0, adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  call adios_select_method(adios_group, ADIOS_TRANSPORT_METHOD, "", "", adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  if (SIMULATION_TYPE == 3) then
    ! crust mantle
    if (ANISOTROPIC_KL) then

      ! anisotropic kernels
      local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT

      ! outputs transverse isotropic kernels only
      if (SAVE_TRANSVERSE_KL_ONLY) then
        call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", "alphav_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", "alphah_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", "betav_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", "betah_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", "eta_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", "rho_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", "bulk_c_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", "bulk_betav_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", "bulk_betah_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", "alpha_kl_crust_mantle", dummy_real4d)
        call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(beta_kl_crust_mantle))
        call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", "bulk_beta_kl_crust_mantle", dummy_real4d)

      else
        call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(rho_kl_crust_mantle))
        local_dim = 21 * NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT
        call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(cijkl_kl_crust_mantle))
      endif

    else

      ! isotropic kernels
      local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", "rhonotprime_kl_crust_mantle", dummy_real4d)
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", "kappa_kl_crust_mantle", dummy_real4d)
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", "mu_kl_crust_mantle", dummy_real4d)
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(rho_kl_crust_mantle))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(alpha_kl_crust_mantle))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(beta_kl_crust_mantle))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", "bulk_c_kl_crust_mantle", dummy_real4d)
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", "bulk_beta_kl_crust_mantle", dummy_real4d)
    endif

    ! noise strength kernel
    if (NOISE_TOMOGRAPHY == 3) then
      local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(sigma_kl_crust_mantle))
    endif

    ! outer core
    local_dim = NSPEC_OUTER_CORE * NGLLX * NGLLY * NGLLZ
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(rho_kl_outer_core))
      call define_adios_global_array1D(adios_group, group_size_inc,local_dim, "", STRINGIFY_VAR(alpha_kl_outer_core))
    if (deviatoric_outercore) then
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(beta_kl_outer_core))
    endif

    ! inner core
    local_dim = NSPEC_INNER_CORE * NGLLX * NGLLY * NGLLZ
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(rho_kl_inner_core))
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(alpha_kl_inner_core))
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(beta_kl_inner_core))

    ! boundary kernel
    if (SAVE_BOUNDARY_MESH) then
      !call save_kernels_boundary_kl()
      if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
        local_dim = NSPEC2D_MOHO * NGLLX * NGLLY * NDIM
        call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(moho_kl))
      endif
      local_dim = NSPEC2D_400 * NGLLX * NGLLY
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(d400_kl))

      local_dim = NSPEC2D_670 * NGLLX * NGLLY
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(d670_kl))

      local_dim = NSPEC2D_CMB * NGLLX * NGLLY
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(cmb_kl))

      local_dim = NSPEC2D_ICB * NGLLX * NGLLY
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(icb_kl))
    endif

    ! approximate hessian
    if (APPROXIMATE_HESS_KL) then
      !call save_kernels_hessian()
      local_dim = NSPEC_CRUST_MANTLE_ADJOINT* NGLLX * NGLLY * NGLLZ
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(hess_kl_crust_mantle))
    endif
  endif

  ! save source derivatives for adjoint simulations
  if (SIMULATION_TYPE == 2 .and. nrec_local > 0) then
    local_dim = 3 * 3 * nrec_local
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(moment_der))

    local_dim = 3 * nrec_local
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(sloc_der))

    local_dim = nrec_local
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(stshift_der))

    local_dim = nrec_local
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(shdur_der))
  endif

  ! Open the handle to file containing all the ADIOS variables
  ! previously defined
  call adios_open (adios_handle, group_name, outputname, "w", comm, adios_err);
  if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed'

  call adios_group_size (adios_handle, group_size_inc, adios_totalsize, adios_err)
  if (adios_err /= 0 ) stop 'Error calling adios_group_size() routine failed'

end subroutine define_kernel_adios_variables


!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the crust mantle.
!! \param[INOUT] adios_handle The handle pointing on the open ADIOS file
!!                            intended to store kernels data.
subroutine write_kernels_cm_adios(adios_handle, &
                                            mu_kl_crust_mantle, kappa_kl_crust_mantle, rhonotprime_kl_crust_mantle, &
                                            alphav_kl_crust_mantle,alphah_kl_crust_mantle, &
                                            betav_kl_crust_mantle,betah_kl_crust_mantle, &
                                            eta_kl_crust_mantle, &
                                            bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle, &
                                            bulk_betav_kl_crust_mantle,bulk_betah_kl_crust_mantle)

  use adios_write_mod

  use specfem_par
  use specfem_par_crustmantle

  use adios_helpers_mod

  implicit none

  ! Parameters
  integer(kind=8), intent(INOUT) :: adios_handle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
      mu_kl_crust_mantle, kappa_kl_crust_mantle, rhonotprime_kl_crust_mantle, &
      alphav_kl_crust_mantle,alphah_kl_crust_mantle, &
      betav_kl_crust_mantle,betah_kl_crust_mantle, &
      eta_kl_crust_mantle, &
      bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle, &
      bulk_betav_kl_crust_mantle,bulk_betah_kl_crust_mantle

  ! Variables
  integer :: local_dim !, adios_err
  integer :: sizeprocs

  ! number of MPI processes
  call world_size(sizeprocs)

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT

  ! For anisotropic kernels
  if (ANISOTROPIC_KL) then
    ! outputs transverse isotropic kernels only
    if (SAVE_TRANSVERSE_KL_ONLY) then
      ! transverse isotropic kernels
      ! (alpha_v, alpha_h, beta_v, beta_h, eta, rho ) parameterization
      call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(alphav_kl_crust_mantle))
      call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(alphah_kl_crust_mantle))
      call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(betav_kl_crust_mantle))
      call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(betah_kl_crust_mantle))
      call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(eta_kl_crust_mantle))
      call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(rho_kl_crust_mantle))
      call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(bulk_c_kl_crust_mantle))
      call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(bulk_betav_kl_crust_mantle))
      call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(bulk_betah_kl_crust_mantle))
      call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(alpha_kl_crust_mantle))
      call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(beta_kl_crust_mantle))
      call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                       STRINGIFY_VAR(bulk_beta_kl_crust_mantle))
    else
      ! note: the C_ij and density kernels are not for relative perturbations
      !       (delta ln( m_i) = delta m_i / m_i),
      !       but absolute perturbations (delta m_i = m_i - m_0)
      call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                       "rho_kl_crust_mantle", -rho_kl_crust_mantle)

      local_dim = 21 * NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT
      call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                       "cijkl_kl_crust_mantle",-cijkl_kl_crust_mantle)
    endif
  else
    ! primary kernels: (rho,kappa,mu) parameterization
    call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(rhonotprime_kl_crust_mantle))
    call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(kappa_kl_crust_mantle))
    call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(mu_kl_crust_mantle))

    ! (rho, alpha, beta ) parameterization
    call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(rho_kl_crust_mantle))
    call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(alpha_kl_crust_mantle))
    call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(beta_kl_crust_mantle))

    ! (rho, bulk, beta ) parameterization, K_rho same as above
    call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(bulk_c_kl_crust_mantle))
    call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                     STRINGIFY_VAR(bulk_beta_kl_crust_mantle))
  endif

end subroutine write_kernels_cm_adios


!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the outer core.
!! \param[INOUT] adios_handle The handle pointing on the open ADIOS file
!!                            intended to store kernels data.
subroutine write_kernels_oc_adios(adios_handle)

  use adios_write_mod

  use specfem_par
  use specfem_par_outercore

  use adios_helpers_mod

  implicit none

  ! Parameters
  integer(kind=8), intent(INOUT) :: adios_handle
  ! Variables
  integer :: local_dim !, adios_err
  integer :: sizeprocs

  ! number of MPI processes
  call world_size(sizeprocs)

  local_dim = NSPEC_OUTER_CORE * NGLLX* NGLLY * NGLLZ

  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(rho_kl_outer_core))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(alpha_kl_outer_core))

  !deviatoric kernel check
  if (deviatoric_outercore) then
    call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(beta_kl_outer_core))
  endif

end subroutine write_kernels_oc_adios


!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the inner core.
!! \param[INOUT] adios_handle The handle pointing on the open ADIOS file
!!                            intended to store kernels data.
subroutine write_kernels_ic_adios(adios_handle)

  use adios_write_mod

  use specfem_par
  use specfem_par_innercore

  use adios_helpers_mod

  implicit none

  ! Parameters
  integer(kind=8), intent(INOUT) :: adios_handle
  ! Variables
  integer :: local_dim !, adios_err
  integer :: sizeprocs

  ! number of MPI processes
  call world_size(sizeprocs)

  local_dim = NSPEC_INNER_CORE * NGLLX * NGLLY * NGLLZ

  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(rho_kl_inner_core))

  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(alpha_kl_inner_core))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(beta_kl_inner_core))

end subroutine write_kernels_ic_adios


!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the boundaries.
!! \param[INOUT] adios_handle The handle pointing on the open ADIOS file
!!                            intended to store kernels data.
subroutine write_kernels_boundary_kl_adios(adios_handle)

  use adios_write_mod

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore

  use adios_helpers_mod

  implicit none

  ! Parameters
  integer(kind=8), intent(INOUT) :: adios_handle
  ! Variables
  integer :: local_dim !, adios_err
  integer :: sizeprocs

  ! number of MPI processes
  call world_size(sizeprocs)

  if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
    local_dim = NSPEC2D_MOHO * NGLLX * NGLLY * NDIM

    call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(moho_kl))
  endif

  local_dim = NSPEC2D_400 * NGLLX * NGLLY
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(d400_kl))

  local_dim = NSPEC2D_670 * NGLLX * NGLLY
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(d670_kl))

  local_dim = NSPEC2D_CMB * NGLLX * NGLLY
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(cmb_kl))

  local_dim = NSPEC2D_ICB * NGLLX * NGLLY
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(icb_kl))

end subroutine write_kernels_boundary_kl_adios


!==============================================================================
!> Schedule writes for the source derivatives (moment tensors and source
!! locations.
!! \param[INOUT] adios_handle The handle pointing on the open ADIOS file
!!                            intended to store kernels data.
!!
!! \note Not to dump one value at a time (as in the non ADIOS version) data are
!!       scaled and oriented but are not written in the same order than in the
!!       non ADIOS version.
!!       (see save_kernels_source_derivatives in save_kernels.f90)
subroutine write_kernels_source_derivatives_adios(adios_handle)

  use adios_write_mod

  use specfem_par

  use adios_helpers_mod

  implicit none

  ! Parameters
  integer(kind=8), intent(INOUT) :: adios_handle
  ! Variables
  ! We do not want to change moment_der and sloc as it might introduce future
  ! concerns if we want to use them after.
  ! No use of modified moment_der values since it implies to allocate those
  ! arrays in save_kernels() and to carry them along the way. It might be better
  ! to transform these arrays in the post processing phase.
  !real(kind=CUSTOM_REAL), dimension(3,3,nrec_local) :: moment_der_tmp
  !real(kind=CUSTOM_REAL), dimension(3,nrec_local) :: sloc_der_tmp
  integer :: local_dim !, adios_err
  integer :: sizeprocs

  ! number of MPI processes
  call world_size(sizeprocs)

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
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(moment_der))

  local_dim = 3 * nrec_local
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(sloc_der))

  local_dim = nrec_local
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(stshift_der))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(shdur_der))

end subroutine write_kernels_source_derivatives_adios


!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the Hessian.
!! \param[INOUT] adios_handle The handle pointing on the open ADIOS file
!!                            intended to store kernels data.
subroutine write_kernels_hessian_adios(adios_handle)

  use adios_write_mod

  use specfem_par
  use specfem_par_crustmantle

  use adios_helpers_mod

  implicit none

  ! Parameters
  integer(kind=8), intent(INOUT) :: adios_handle
  ! Variables
  integer :: local_dim !, adios_err
  integer :: sizeprocs

  ! number of MPI processes
  call world_size(sizeprocs)

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT

  ! stores into file
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(hess_kl_crust_mantle))

end subroutine write_kernels_hessian_adios

!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the noise strength kernel.
!! \param[INOUT] adios_handle The handle pointing on the open ADIOS file
!!                            intended to store kernels data.
subroutine write_kernels_strength_noise_adios(adios_handle)

  use adios_write_mod

  use specfem_par
  use specfem_par_crustmantle

  use adios_helpers_mod

  implicit none

  ! Parameters
  integer(kind=8), intent(INOUT) :: adios_handle
  ! Variables
  integer :: local_dim
  integer :: sizeprocs

  ! number of MPI processes
  call world_size(sizeprocs)

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT

  ! stores into file
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(sigma_kl_crust_mantle))

end subroutine write_kernels_strength_noise_adios
