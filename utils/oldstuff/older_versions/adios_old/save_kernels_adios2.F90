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

module kernel_adios2

  use adios2, only: adios2_engine, adios2_io, adios2_variable, adios2_attribute

  use manager_adios2, only: adios2_file_handle

  implicit none

  ! adios2 variables
  type(adios2_variable), public :: v_alphav_kl_crust_mantle, v_alphah_kl_crust_mantle, v_betav_kl_crust_mantle, &
                                   v_betah_kl_crust_mantle, v_eta_kl_crust_mantle, v_rho_kl_crust_mantle, &
                                   v_bulk_c_kl_crust_mantle, v_bulk_betav_kl_crust_mantle, v_bulk_betah_kl_crust_mantle, &
                                   v_Gc_prime_kl_crust_mantle, v_Gs_prime_kl_crust_mantle, &
                                   v_alpha_kl_crust_mantle, v_cijkl_kl_crust_mantle, v_rhonotprime_kl_crust_mantle, &
                                   v_kappa_kl_crust_mantle,  v_mu_kl_crust_mantle, &
                                   v_beta_kl_crust_mantle,  v_sigma_kl_crust_mantle, v_bulk_beta_kl_crust_mantle, &
                                   v_hess_kl_crust_mantle, v_moment_der, v_sloc_der, v_stshift_der, v_shdur_der

  type(adios2_variable), public :: v_rho_kl_outer_core, v_alpha_kl_outer_core, v_beta_kl_outer_core
  type(adios2_variable), public :: v_rho_kl_inner_core, v_alpha_kl_inner_core, v_beta_kl_inner_core
  type(adios2_variable), public :: v_moho_kl, v_d400_kl, v_d670_kl, v_cmb_kl, v_icb_kl,

end module kernel_adios2


!-------------------------------------------------------------------------------
!> \file save_kernels_adios2.f90
!! \brief Save kernels arrays to file with the help of the ADIOS2 library.
!! \author MPBL
!-------------------------------------------------------------------------------

#include "config.fh"

!==============================================================================
!> Define all the kernels that will be written to the ADIOS file.
!!
!! \note Everything is define in this single function, even the group size.
!!       It is the reason why this function require only an handle on an ADIOS
!!       file as an argument.
  subroutine define_kernel_adios2_variables()

  use specfem_par ! Just for dimensions. No need of arrays for now.
  use specfem_par_crustmantle
  use specfem_par_outercore
  use specfem_par_innercore
  use specfem_par_noise

  use adios2
  use manager_adios2
  use kernel_adios2

  implicit none

  ! local Variables
  integer :: ier
  type(adios2_io) :: adios2_group
  integer(kind=8), dimension(1) :: gdim  ! Global shape of array
  integer(kind=8), dimension(1) :: ldim  ! Local size of array
  integer(kind=8), dimension(1) :: offs  ! Starting offset in global array
  character(len=MAX_STRING_LEN) :: outputname, group_name
  integer :: ier

  ! opens file
  outputname = trim(OUTPUT_FILES)//"/kernels.bp"

  group_name = "SPECFEM3D_GLOBE_KERNELS"

  call init_adios_group(adios2_group,group_name)
  call open_file_adios_write(myadios_file,myadios_group,outputname,adios2_group)

  ! defines adios arrays
  if (SIMULATION_TYPE == 3) then
    ! crust mantle
    if (ANISOTROPIC_KL) then

      ! anisotropic kernels
      ldim(1) = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT
      gdim = sizeprocs_adios2 * ldim;
      offs = myrank_adios2 * ldim;

      ! outputs transverse isotropic kernels only
      if (SAVE_TRANSVERSE_KL_ONLY) then
        call adios2_define_variable(v_alphav_kl_crust_mantle, adios2_group, "alphav_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_alphah_kl_crust_mantle, adios2_group, "alphah_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_betav_kl_crust_mantle, adios2_group, "betav_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_betah_kl_crust_mantle, adios2_group, "betah_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_eta_kl_crust_mantle, adios2_group, "eta_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_rho_kl_crust_mantle, adios2_group, "rho_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_bulk_c_kl_crust_mantle, adios2_group, "bulk_c_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_bulk_betav_kl_crust_mantle, adios2_group, "bulk_betav_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_bulk_betah_kl_crust_mantle, adios2_group, "bulk_betah_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_alpha_kl_crust_mantle, adios2_group, "alpha_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_beta_kl_crust_mantle, adios2_group, "beta_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_bulk_beta_kl_crust_mantle, adios2_group, "bulk_beta_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

      else if (SAVE_AZIMUTHAL_ANISO_KL_ONLY) then
        call adios2_define_variable(v_betav_kl_crust_mantle, adios2_group, "betav_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_betah_kl_crust_mantle, adios2_group, "betah_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_eta_kl_crust_mantle, adios2_group, "eta_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_rho_kl_crust_mantle, adios2_group, "rho_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_bulk_c_kl_crust_mantle, adios2_group, "bulk_c_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_bulk_betav_kl_crust_mantle, adios2_group, "bulk_betav_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_bulk_betah_kl_crust_mantle, adios2_group, "bulk_betah_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_alpha_kl_crust_mantle, adios2_group, "alpha_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_beta_kl_crust_mantle, adios2_group, "beta_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_bulk_beta_kl_crust_mantle, adios2_group, "bulk_beta_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_Gc_prime_kl_crust_mantle, adios2_group, "Gc_prime_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
        call adios2_define_variable(v_Gs_prime_kl_crust_mantle, adios2_group, "Gs_prime_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

      else

        call adios2_define_variable(v_rho_kl_crust_mantle, adios2_group, "rho_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

        ldim(1) = 21 * NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT
        gdim = sizeprocs_adios2 * ldim;
        offs = myrank_adios2 * ldim;

        call adios2_define_variable(v_cijkl_kl_crust_mantle, adios2_group, "cijkl_kl_crust_mantle", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)


      endif

    else

      ! isotropic kernels
      ldim(1) = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT
      gdim = sizeprocs_adios2 * ldim;
      offs = myrank_adios2 * ldim;

      call adios2_define_variable(v_rhonotprime_kl_crust_mantle, adios2_group, "rhonotprime_kl_crust_mantle", &
                                  adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
      call adios2_define_variable(v_kappa_kl_crust_mantle, adios2_group, "kappa_kl_crust_mantle", &
                                  adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
      call adios2_define_variable(v_mu_kl_crust_mantle, adios2_group, "mu_kl_crust_mantle", &
                                  adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
      call adios2_define_variable(v_rho_kl_crust_mantle, adios2_group, "rho_kl_crust_mantle", &
                                  adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
      call adios2_define_variable(v_alpha_kl_crust_mantle, adios2_group, "alpha_kl_crust_mantle", &
                                  adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
      call adios2_define_variable(v_beta_kl_crust_mantle, adios2_group, "beta_kl_crust_mantle", &
                                  adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
      call adios2_define_variable(v_bulk_c_kl_crust_mantle, adios2_group, "v_bulk_c_kl_crust_mantle", &
                                  adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
      call adios2_define_variable(v_bulk_beta_kl_crust_mantle, adios2_group, "bulk_beta_kl_crust_mantle", &
                                  adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

    endif

    ! noise strength kernel
    if (NOISE_TOMOGRAPHY == 3) then
      ldim(1) = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ADJOINT
      gdim = sizeprocs_adios2 * ldim;
      offs = myrank_adios2 * ldim;

      call adios2_define_variable(v_sigma_kl_crust_mantle, adios2_group, "sigma_kl_crust_mantle", &
                                  adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

    endif

    ! outer core
    if (SAVE_KERNELS_OC) then
      ldim(1) = NSPEC_OUTER_CORE * NGLLX * NGLLY * NGLLZ
      gdim = sizeprocs_adios2 * ldim;
      offs = myrank_adios2 * ldim;

      call adios2_define_variable(v_rho_kl_outer_core, adios2_group, "rho_kl_outer_core", &
                                  adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
      call adios2_define_variable(v_alpha_kl_outer_core, adios2_group, "alpha_kl_outer_core", &
                                  adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
      if (deviatoric_outercore) then
        call adios2_define_variable(v_beta_kl_outer_core, adios2_group, "beta_kl_outer_core", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
      endif
    endif

    ! inner core
    if (SAVE_KERNELS_IC) then
      ldim(1) = NSPEC_INNER_CORE * NGLLX * NGLLY * NGLLZ
      gdim = sizeprocs_adios2 * ldim;
      offs = myrank_adios2 * ldim;
      call adios2_define_variable(v_rho_kl_inner_core, adios2_group, "rho_kl_inner_core", &
                                  adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
      call adios2_define_variable(v_alpha_kl_inner_core, adios2_group, "alpha_kl_inner_core", &
                                  adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
      call adios2_define_variable(v_beta_kl_inner_core, adios2_group, "beta_kl_inner_core", &
                                  adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
    endif

    ! boundary kernel
    if (SAVE_KERNELS_BOUNDARY) then
      !call save_kernels_boundary_kl()
      if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
        ldim(1) = NSPEC2D_MOHO * NGLLX * NGLLY * NDIM
        gdim = sizeprocs_adios2 * ldim;
        offs = myrank_adios2 * ldim;
        call adios2_define_variable(v_moho_kl, adios2_group, "moho_kl", &
                                    adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
      endif
      ldim(1) = NSPEC2D_400 * NGLLX * NGLLY
      gdim = sizeprocs_adios2 * ldim;
      offs = myrank_adios2 * ldim;
      call adios2_define_variable(v_d400_kl, adios2_group, "d400_kl", &
                                adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

      ldim(1) = NSPEC2D_670 * NGLLX * NGLLY
      gdim = sizeprocs_adios2 * ldim;
      offs = myrank_adios2 * ldim;
      call adios2_define_variable(v_d670_kl, adios2_group, "d670_kl", &
                                  adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

      ldim(1) = NSPEC2D_CMB * NGLLX * NGLLY
      gdim = sizeprocs_adios2 * ldim;
      offs = myrank_adios2 * ldim;
      call adios2_define_variable(v_cmb_kl, adios2_group, "cmb_kl", &
                                  adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

      ldim(1) = NSPEC2D_ICB * NGLLX * NGLLY
      gdim = sizeprocs_adios2 * ldim;
      offs = myrank_adios2 * ldim;
      call adios2_define_variable(v_icb_kl, adios2_group, "icb_kl", &
                                  adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
    endif

    ! approximate Hessian
    if (APPROXIMATE_HESS_KL) then
      !call save_kernels_Hessian()
      ldim(1) = NSPEC_CRUST_MANTLE_ADJOINT* NGLLX * NGLLY * NGLLZ
      gdim = sizeprocs_adios2 * ldim;
      offs = myrank_adios2 * ldim;
      call adios2_define_variable(v_hess_kl_crust_mantle, adios2_group, "hess_kl_crust_mantle", &
                                  adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
    endif
  endif

  ! save source derivatives for adjoint simulations
  if (SIMULATION_TYPE == 2 .and. nrec_local > 0) then
    ldim(1) = 3 * 3 * nrec_local
    gdim = sizeprocs_adios2 * ldim;
    offs = myrank_adios2 * ldim;
    call adios2_define_variable(v_moment_der, adios2_group, "moment_der", &
                                adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

    ldim(1) = 3 * nrec_local
    gdim = sizeprocs_adios2 * ldim;
    offs = myrank_adios2 * ldim;
    call adios2_define_variable(v_sloc_der, adios2_group, "sloc_der", &
                                adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

    ldim(1) = nrec_local
    gdim = sizeprocs_adios2 * ldim;
    offs = myrank_adios2 * ldim;
    call adios2_define_variable(v_stshift_der, adios2_group, "stshift_der", &
                                adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

    ldim(1) = nrec_local
    gdim = sizeprocs_adios2 * ldim;
    offs = myrank_adios2 * ldim;
    call adios2_define_variable(v_shdur_der, adios2_group, "shdur_der", &
                                adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  endif

  end subroutine define_kernel_adios2_variables


!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the crust mantle.
  subroutine write_kernels_cm_ani_adios2(alphav_kl_crust_mantle,alphah_kl_crust_mantle, &
                                        betav_kl_crust_mantle,betah_kl_crust_mantle, &
                                        eta_kl_crust_mantle, &
                                        bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle, &
                                        bulk_betav_kl_crust_mantle,bulk_betah_kl_crust_mantle, &
                                        Gc_prime_kl_crust_mantle, Gs_prime_kl_crust_mantle)

  use specfem_par
  use specfem_par_crustmantle

  use adios2
  use kernel_adios2

  implicit none

  ! Parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
      alphav_kl_crust_mantle,alphah_kl_crust_mantle, &
      betav_kl_crust_mantle,betah_kl_crust_mantle, &
      eta_kl_crust_mantle, &
      bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle, &
      bulk_betav_kl_crust_mantle,bulk_betah_kl_crust_mantle, &
      Gc_prime_kl_crust_mantle, Gs_prime_kl_crust_mantle

  ! Variables
  integer :: ier

  ! checks if anything to do
  if (.not. ANISOTROPIC_KL) return

  ! For anisotropic kernels
  ! outputs transverse isotropic kernels only
  if (SAVE_TRANSVERSE_KL_ONLY) then
    ! transverse isotropic kernels
    ! (alpha_v, alpha_h, beta_v, beta_h, eta, rho ) parameterization
    call adios2_put(adios2_file_handle, v_alphav_kl_crust_mantle, alphav_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_alphah_kl_crust_mantle, alphah_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_betav_kl_crust_mantle, betav_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_betah_kl_crust_mantle, betah_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_eta_kl_crust_mantle, eta_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_rho_kl_crust_mantle, rho_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_bulk_c_kl_crust_mantle, bulk_c_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_bulk_betav_kl_crust_mantle, bulk_betav_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_bulk_betah_kl_crust_mantle, bulk_betah_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_alpha_kl_crust_mantle, alpha_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_beta_kl_crust_mantle, beta_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_bulk_beta_kl_crust_mantle, bulk_beta_kl_crust_mantle, adios2_mode_sync, ier)

  else if (SAVE_AZIMUTHAL_ANISO_KL_ONLY) then
    ! kernels for inversions involving azimuthal anisotropy
    ! (bulk_c, beta_v, beta_h, eta, Gc', Gs', rho ) parameterization
    ! note: Gc' & Gs' are the normalized Gc & Gs kernels
    call adios2_put(adios2_file_handle, v_betav_kl_crust_mantle, betav_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_betah_kl_crust_mantle, betah_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_eta_kl_crust_mantle, eta_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_rho_kl_crust_mantle, rho_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_bulk_c_kl_crust_mantle, bulk_c_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_bulk_betav_kl_crust_mantle, bulk_betav_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_bulk_betah_kl_crust_mantle, bulk_betah_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_alpha_kl_crust_mantle, alpha_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_beta_kl_crust_mantle, beta_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_bulk_beta_kl_crust_mantle, bulk_beta_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_Gc_prime_kl_crust_mantle, Gc_prime_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_Gs_prime_kl_crust_mantle, Gs_prime_kl_crust_mantle, adios2_mode_sync, ier)

  else
    ! note: the C_ij and density kernels are not for relative perturbations
    !       (delta ln( m_i) = delta m_i / m_i),
    !       but absolute perturbations (delta m_i = m_i - m_0)
    ! Note: adios2_put() by default is in deferred mode which does not work if the array
    ! disappears before adios2_end_step/adios2_close
    ! For temporary arrays, we must use adios2_mode_sync as an extra argument
    call adios2_put(adios2_file_handle, v_rho_kl_crust_mantle, -rho_kl_crust_mantle, adios2_mode_sync, ier)
    call adios2_put(adios2_file_handle, v_cijkl_kl_crust_mantle, -cijkl_kl_crust_mantle, adios2_mode_sync, ier)
  endif

  end subroutine write_kernels_cm_ani_adios2

!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the crust mantle.
  subroutine write_kernels_cm_iso_adios2(mu_kl_crust_mantle, kappa_kl_crust_mantle, rhonotprime_kl_crust_mantle, &
                                        bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle)

  use specfem_par
  use specfem_par_crustmantle

  use adios2
  use kernel_adios2

  implicit none

  ! Parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
      mu_kl_crust_mantle, kappa_kl_crust_mantle, rhonotprime_kl_crust_mantle, &
      bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle

  ! Variables
  integer :: ier

  ! checks if anything to do
  if (ANISOTROPIC_KL) return

  ! isotropic kernels
  ! primary kernels: (rho,kappa,mu) parameterization
  call adios2_put(adios2_file_handle, v_rhonotprime_kl_crust_mantle, rhonotprime_kl_crust_mantle, adios2_mode_sync, ier)
  call adios2_put(adios2_file_handle, v_kappa_kl_crust_mantle, kappa_kl_crust_mantle, adios2_mode_sync, ier)
  call adios2_put(adios2_file_handle, v_mu_kl_crust_mantle, mu_kl_crust_mantle, adios2_mode_sync, ier)

  ! (rho, alpha, beta ) parameterization
  call adios2_put(adios2_file_handle, v_rho_kl_crust_mantle, rho_kl_crust_mantle, adios2_mode_sync, ier)
  call adios2_put(adios2_file_handle, v_alpha_kl_crust_mantle, alpha_kl_crust_mantle, adios2_mode_sync, ier)
  call adios2_put(adios2_file_handle, v_beta_kl_crust_mantle, beta_kl_crust_mantle, adios2_mode_sync, ier)

  ! (rho, bulk, beta ) parameterization, K_rho same as above
  call adios2_put(adios2_file_handle, v_bulk_c_kl_crust_mantle, bulk_c_kl_crust_mantle, adios2_mode_sync, ier)
  call adios2_put(adios2_file_handle, v_bulk_beta_kl_crust_mantle, bulk_beta_kl_crust_mantle, adios2_mode_sync, ier)

  end subroutine write_kernels_cm_iso_adios2


!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the outer core.
  subroutine write_kernels_oc_adios2()

  use specfem_par
  use specfem_par_outercore

  use adios2
  use kernel_adios2

  implicit none

  ! Variables
  integer :: ier

  ! checks if anything to do
  if (.not. SAVE_KERNELS_OC) return

  call adios2_put(adios2_file_handle, v_rho_kl_outer_core, rho_kl_outer_core, adios2_mode_sync, ier)
  call adios2_put(adios2_file_handle, v_alpha_kl_outer_core, alpha_kl_outer_core, adios2_mode_sync, ier)

  !deviatoric kernel check
  if (deviatoric_outercore) then
    call adios2_put(adios2_file_handle, v_beta_kl_outer_core, beta_kl_outer_core, adios2_mode_sync, ier)
  endif

  end subroutine write_kernels_oc_adios2


!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the inner core.
  subroutine write_kernels_ic_adios2()

  use specfem_par
  use specfem_par_innercore

  use adios2
  use kernel_adios2

  implicit none

  ! Variables
  integer :: ier

  ! checks if anything to do
  if (.not. SAVE_KERNELS_IC) return

  call adios2_put(adios2_file_handle, v_rho_kl_inner_core, rho_kl_inner_core, adios2_mode_sync, ier)
  call adios2_put(adios2_file_handle, v_alpha_kl_inner_core, alpha_kl_inner_core, adios2_mode_sync, ier)
  call adios2_put(adios2_file_handle, v_beta_kl_inner_core, beta_kl_inner_core, adios2_mode_sync, ier)

  end subroutine write_kernels_ic_adios2


!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the boundaries.
  subroutine write_kernels_boundary_kl_adios2()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore

  use adios2
  use kernel_adios2

  implicit none

  ! Variables
  integer :: ier

  ! checks if anything to do
  if (.not. SAVE_KERNELS_BOUNDARY) return

  if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
    call adios2_put(adios2_file_handle, v_moho_kl, moho_kl, adios2_mode_sync, ier)
  endif

  call adios2_put(adios2_file_handle, v_d400_kl, d400_kl, adios2_mode_sync, ier)
  call adios2_put(adios2_file_handle, v_d670_kl, d670_kl, adios2_mode_sync, ier)
  call adios2_put(adios2_file_handle, v_cmb_kl, cmb_kl, adios2_mode_sync, ier)
  call adios2_put(adios2_file_handle, v_icb_kl, icb_kl, adios2_mode_sync, ier)

  end subroutine write_kernels_boundary_kl_adios2


!==============================================================================
!> Schedule writes for the source derivatives (moment tensors and source
!! locations.
!!
!! \note Not to dump one value at a time (as in the non ADIOS version) data are
!!       scaled and oriented but are not written in the same order than in the
!!       non ADIOS version.
!!       (see save_kernels_source_derivatives in save_kernels.f90)
  subroutine write_kernels_source_derivatives_adios2()

  use specfem_par

  use adios2
  use kernel_adios2

  implicit none

  ! Variables
  ! We do not want to change moment_der and sloc as it might introduce future
  ! concerns if we want to use them after.
  ! No use of modified moment_der values since it implies to allocate those
  ! arrays in save_kernels() and to carry them along the way. It might be better
  ! to transform these arrays in the post processing phase.
  !real(kind=CUSTOM_REAL), dimension(3,3,nrec_local) :: moment_der_tmp
  !real(kind=CUSTOM_REAL), dimension(3,nrec_local) :: sloc_der_tmp
  integer :: ier

  !moment_der_tmp(:, :, :) = moment_der(:, :, :) * 1e-7
  !moment_der_tmp(1, 3, :) = -2 * moment_der(1, 3, :)
  !moment_der_tmp(2, 3, :) =  2 * moment_der(2, 3, :)
  !moment_der_tmp(1, 2, :) = -2 * moment_der(1, 2, :)
  !moment_der_tmp(3, 1, :) = moment_der(1, 3, :)
  !moment_der_tmp(3, 2, :) = moment_der(2, 3, :)
  !moment_der_tmp(2, 1, :) = moment_der(1, 2, :)
  !sloc_der_tmp(:, :) = - sloc_der(:, :)
  !sloc_der_tmp(3, :) = - sloc_der(3, :)


  ! Note: adios2_put() by default is in deferred mode which does not work if the array
  ! disappears before adios2_end_step/adios2_close
  ! For temporary arrays, use adios2_mode_sync as an extra argument

  call adios2_put(adios2_file_handle, v_moment_der, moment_der, adios2_mode_sync, ier)
  call adios2_put(adios2_file_handle, v_sloc_der, sloc_der, adios2_mode_sync, ier)
  call adios2_put(adios2_file_handle, v_stshift_der, stshift_der, adios2_mode_sync, ier)
  call adios2_put(adios2_file_handle, v_shdur_der, shdur_der, adios2_mode_sync, ier)

  end subroutine write_kernels_source_derivatives_adios2


!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the Hessian.
  subroutine write_kernels_Hessian_adios2()

  use specfem_par
  use specfem_par_crustmantle

  use adios2
  use kernel_adios2

  implicit none

  ! Variables
  integer :: ier

  ! stores into file
  call adios2_put(adios2_file_handle, v_hess_kl_crust_mantle, hess_kl_crust_mantle, adios2_mode_sync, ier)

  end subroutine write_kernels_Hessian_adios2

!==============================================================================
!> Schedule ADIOS writes for kernel variables related to the noise strength kernel.
  subroutine write_kernels_strength_noise_adios2()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_noise

  use adios2
  use kernel_adios2

  implicit none

  ! Variables
  integer :: ier

  ! stores into file
  call adios2_put(adios2_file_handle, v_sigma_kl_crust_mantle, sigma_kl_crust_mantle, adios2_mode_sync, ier)

  end subroutine write_kernels_strength_noise_adios2
