/*
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
*/

#include <stdio.h>
#include <stdlib.h>

#include "config.h"

typedef float realw;


// for xmeshfem3D compilation

void FC_FUNC_(get_absorb_adios,GET_ABSORB_ADIOS)(int* iregion,
                                                 int* nimin, int* nimax, int* njmin, int* njmax, int* nkmin_xi, int* nkmin_eta,
                                                 int* NSPEC2DMAX_XMIN_XMAX, int* NSPEC2DMAX_YMIN_YMAX) {}

void FC_FUNC_(read_gll_model_adios,READ_GLL_MODEL_ADIOS)(void) {}

void FC_FUNC_(save_arrays_solver_adios,SAVE_ARRAYS_SOLVER_ADIOS)(int* idoubling, int* ibool,
                                                                 int* iregion_code,
                                                                 realw* xstore, realw* ystore, realw* zstore,
                                                                 int* NSPEC2DMAX_XMIN_XMAX, int* NSPEC2DMAX_YMIN_YMAX,
                                                                 int* NSPEC2D_TOP, int* NSPEC2D_BOTTOM) {}

void FC_FUNC_(save_arrays_solver_meshfiles_adios,SAVE_ARRAYS_SOLVER_MESHFILES_ADIOS)(void) {}

void FC_FUNC_(save_arrays_boundary_adios,SAVE_ARRAYS_BOUNDARY_ADIOS)(void) {}

void FC_FUNC_(save_mpi_arrays_adios,SAVE_MPI_ARRAYS_ADIOS)(void) {}

void FC_FUNC_(write_avs_dx_output_adios,WRITE_AVS_DX_OUTPUT_ADIOS)(int* npointot, int* iregion_code,
                                                                   int* num_ibool_AVS_DX, int* mask_ibool) {}

// for xspecfem3D compilation

void FC_FUNC_(read_arrays_solver_adios,READ_ARRAYS_SOLVER_ADIOS)(void) {}

void FC_FUNC_(read_attenuation_adios,READ_ATTENUATION_ADIOS)(void) {}

void FC_FUNC_(read_forward_arrays_adios,READ_FORWARD_ARRAYS_ADIOS)(void) {}

void FC_FUNC_(read_forward_arrays_undoatt_adios,READ_FORWARD_ARRAYS_UNDOATT_ADIOS)(void) {}

void FC_FUNC_(read_intermediate_forward_arrays_adios,READ_INTERMEDIATE_FORWARD_ARRAYS_ADIOS)(void) {}

void FC_FUNC_(read_mesh_databases_coupling_adios,READ_MESH_DATABASES_COUPLING_ADIOS)(void) {}

void FC_FUNC_(read_mesh_databases_mpi_cm_adios,READ_MESH_DATABASES_MPI_CM_ADIOS)(void) {}

void FC_FUNC_(read_mesh_databases_mpi_ic_adios,READ_MESH_DATABASES_MPI_IC_ADIOS)(void) {}

void FC_FUNC_(read_mesh_databases_mpi_oc_adios,READ_MESH_DATABASES_MPI_OC_ADIOS)(void) {}

void FC_FUNC_(read_mesh_databases_stacey_adios,READ_MESH_DATABASES_STACEY_ADIOS)(void) {}

void FC_FUNC_(save_forward_arrays_adios,SAVE_FORWARD_ARRAYS_ADIOS)(void) {}

void FC_FUNC_(save_forward_arrays_undoatt_adios,SAVE_FORWARD_ARRAYS_UNDOATT_ADIOS)(void) {}

void FC_FUNC_(save_intermediate_forward_arrays_adios,SAVE_INTERMEDIATE_FORWARD_ARRAYS_ADIOS)(void) {}

void FC_FUNC_(define_kernel_adios_variables,DEFINE_KERNEL_ADIOS_VARIABLES)(void) {}

void FC_FUNC_(write_kernels_cm_ani_adios,WRITE_KERNELS_CM_ADIOS)(void) {}

void FC_FUNC_(write_kernels_cm_iso_adios,WRITE_KERNELS_CM_ADIOS)(void) {}

void FC_FUNC_(write_kernels_oc_adios,WRITE_KERNELS_OC_ADIOS)(void) {}

void FC_FUNC_(write_kernels_ic_adios,WRITE_KERNELS_IC_ADIOS)(void) {}

void FC_FUNC_(write_kernels_boundary_kl_adios,WRITE_KERNELS_BOUNDARY_KL_ADIOS)(void) {}

void FC_FUNC_(write_kernels_source_derivatives_adios,WRITE_KERNELS_SOURCE_DERIVATIVES_ADIOS)(void) {}

void FC_FUNC_(write_kernels_hessian_adios,WRITE_KERNELS_HESSIAN_ADIOS)(void) {}

void FC_FUNC_(write_kernels_strength_noise_adios,WRITE_KERNELS_STRENGTH_NOISE_ADIOS)(void) {}
