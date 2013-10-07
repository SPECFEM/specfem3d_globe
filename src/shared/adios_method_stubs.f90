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


! placeholders for non-adios compilation

  subroutine warn_no_adios()
    stop 'Trying to use ADIOS while not configured for it.'
  end subroutine warn_no_adios

! modules

!  module AVS_DX_global_mod
!    type avs_dx_global_t
!    end type avs_dx_global_t
!  end module AVS_DX_global_mod

!  module AVS_DX_global_faces_mod
!    type avs_dx_global_faces_t
!    end type avs_dx_global_faces_t
!  end module AVS_DX_global_faces_mod

!  module AVS_DX_global_chunks_mod
!    type avs_dx_global_chunks_t
!    end type avs_dx_global_chunks_t
!  end module AVS_DX_global_chunks_mod

!  module AVS_DX_surface_mod
!    type avs_dx_surface_t
!    end type avs_dx_surface_t
!  end module AVS_DX_surface_mod

! for both xmeshfem3D/xspecfem3D compilation

  subroutine adios_setup()
  call warn_no_adios()
  end subroutine

  subroutine adios_cleanup()
  end subroutine

! for xmeshfem3D compilation

  subroutine crm_save_mesh_files_adios()
  end subroutine

  subroutine get_absorb_adios()
  end subroutine

  subroutine save_arrays_solver_adios()
  end subroutine

  subroutine save_arrays_solver_meshfiles_adios()
  end subroutine

  subroutine save_arrays_solver_boundary_adios()
  end subroutine

  subroutine save_MPI_arrays_adios()
  end subroutine

  subroutine read_gll_model_adios
  end subroutine
  
!  subroutine prepare_AVS_DX_global_chunks_data_adios()
!  end subroutine
!
!  subroutine write_AVS_DX_global_chunks_data_adios()
!  end subroutine
!
!  subroutine free_AVS_DX_global_chunks_data_adios()
!  end subroutine
!
!  subroutine define_AVS_DX_global_data_adios()
!  end subroutine
!
!  subroutine prepare_AVS_DX_global_data_adios()
!  end subroutine
!
!  subroutine write_AVS_DX_global_data_adios()
!  end subroutine
!
!  subroutine free_AVS_DX_global_data_adios()
!  end subroutine
!
!  subroutine define_AVS_DX_global_faces_data_adios()
!  end subroutine
!
!  subroutine prepare_AVS_DX_global_faces_data_adios()
!  end subroutine
!
!  subroutine write_AVS_DX_global_faces_data_adios()
!  end subroutine
!
!  subroutine free_AVS_DX_global_faces_data_adios()
!  end subroutine
!
!  subroutine define_AVS_DX_surfaces_data_adios()
!  end subroutine
!
!  subroutine prepare_AVS_DX_surfaces_data_adios()
!  end subroutine
!
!  subroutine write_AVS_DX_surfaces_data_adios()
!  end subroutine
!
!  subroutine free_AVS_DX_surfaces_data_adios()
!  end subroutine
!

! for xspecfem3D compilation

  subroutine read_arrays_solver_adios()
  end subroutine

  subroutine read_attenuation_adios()
  end subroutine

  subroutine read_forward_arrays_adios()
  end subroutine

  subroutine read_intermediate_forward_arrays_adios()
  end subroutine

  subroutine read_mesh_databases_coupling_adios()
  end subroutine

  subroutine read_mesh_databases_mpi_cm_adios()
  end subroutine

  subroutine read_mesh_databases_mpi_ic_adios()
  end subroutine

  subroutine read_mesh_databases_mpi_oc_adios()
  end subroutine

  subroutine read_mesh_databases_stacey_adios()
  end subroutine

  subroutine save_forward_arrays_adios()
  end subroutine

  subroutine save_intermediate_forward_arrays_adios()
  end subroutine

!  subroutine define_common_forward_arrays_adios()
!  end subroutine
!
!  subroutine define_rotation_forward_arrays_adios()
!  end subroutine
!
!  subroutine define_attenuation_forward_arrays_adios()
!  end subroutine
!
!  subroutine write_common_forward_arrays_adios()
!  end subroutine
!
!  subroutine write_rotation_forward_arrays_adios()
!  end subroutine
!
!  subroutine write_attenuation_forward_arrays_adios()
!  end subroutine
!
!  subroutine write_1D_global_array_adios_dims()
!  end subroutine
!

  subroutine perform_write_adios_kernels()
  end subroutine

  subroutine define_kernel_adios_variables()
  end subroutine

  subroutine write_kernels_crust_mantle_adios()
  end subroutine

  subroutine write_kernels_outer_core_adios()
  end subroutine

  subroutine write_kernels_inner_core_adios()
  end subroutine

  subroutine write_kernels_boundary_kl_adios()
  end subroutine

  subroutine write_kernels_source_derivatives_adios()
  end subroutine

  subroutine write_kernels_hessian_adios()
  end subroutine

!
!  subroutine write_specfem_header_adios()
!  end subroutine

