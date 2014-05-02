/*
!=====================================================================
!
 !          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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
*/

#include <stdio.h>
#include <stdlib.h>

#include "config.h"

typedef float realw;

// placeholders for non-adios compilation

void warn_no_adios() {
 fprintf(stderr,"ERROR: ADIOS enabled without ADIOS Support. To enable ADIOS support, reconfigure with --with-adios flag.\n");
 exit(1);
}

void FC_FUNC_(adios_setup,ADIOS_SETUP)(){ warn_no_adios(); }

void FC_FUNC_(adios_cleanup,ADIOS_CLEANUP)(){}



// for xmeshfem3D compilation

void FC_FUNC_(crm_save_mesh_files_adios,CRM_SAVE_MESH_FILES_ADIOS)(int* nspec, int* npointot, int* iregion_code,
                                                                   int* num_ibool_AVS_DX, int* mask_ibool){}

void FC_FUNC_(get_absorb_adios,GET_ABSORB_ADIOS)(int* myrank, int* iregion,
                                                 int* nimin, int* nimax, int* njmin, int* njmax, int* nkmin_xi, int* nkmin_eta,
                                                 int* NSPEC2DMAX_XMIN_XMAX, int* NSPEC2DMAX_YMIN_YMAX){}

void FC_FUNC_(save_arrays_solver_adios,SAVE_ARRAYS_SOLVER_ADIOS)(int* myrank, int* nspec, int* nglob,
                                                                 int* idoubling, int* ibool,
                                                                 int* iregion_code,
                                                                 realw* xstore, realw* ystore, realw* zstore,
                                                                 int* NSPEC2DMAX_XMIN_XMAX, int* NSPEC2DMAX_YMIN_YMAX,
                                                                 int* NSPEC2D_TOP, int* NSPEC2D_BOTTOM){}

void FC_FUNC_(save_arrays_solver_meshfiles_adios,SAVE_ARRAYS_SOLVER_MESHFILES_ADIOS)(){}

void FC_FUNC_(save_arrays_solver_boundary_adios,SAVE_ARRAYS_SOLVER_BOUNDARY_ADIOS)(){}

void FC_FUNC_(save_mpi_arrays_adios,SAVE_MPI_ARRAYS_ADIOS)(){}

void FC_FUNC_(read_gll_model_adios,READ_GLL_MODEL_ADIOS)(){}

// for xspecfem3D compilation

void FC_FUNC_(read_arrays_solver_adios,READ_ARRAYS_SOLVER_ADIOS)(){}

void FC_FUNC_(read_attenuation_adios,READ_ATTENUATION_ADIOS)(){}

void FC_FUNC_(read_forward_arrays_adios,READ_FORWARD_ARRAYS_ADIOS)(){}

void FC_FUNC_(read_intermediate_forward_arrays_adios,READ_INTERMEDIATE_FORWARD_ARRAYS_ADIOS)(){}

void FC_FUNC_(read_mesh_databases_coupling_adios,READ_MESH_DATABASES_COUPLING_ADIOS)(){}

void FC_FUNC_(read_mesh_databases_mpi_cm_adios,READ_MESH_DATABASES_MPI_CM_ADIOS)(){}

void FC_FUNC_(read_mesh_databases_mpi_ic_adios,READ_MESH_DATABASES_MPI_IC_ADIOS)(){}

void FC_FUNC_(read_mesh_databases_mpi_oc_adios,READ_MESH_DATABASES_MPI_OC_ADIOS)(){}

void FC_FUNC_(read_mesh_databases_stacey_adios,READ_MESH_DATABASES_STACEY_ADIOS)(){}

void FC_FUNC_(save_forward_arrays_adios,SAVE_FORWARD_ARRAYS_ADIOS)(){}

void FC_FUNC_(save_intermediate_forward_arrays_adios,SAVE_INTERMEDIATE_FORWARD_ARRAYS_ADIOS)(){}

void FC_FUNC_(perform_write_adios_kernels,PERFORM_WRITE_ADIOS_KERNELS)(){}

void FC_FUNC_(define_kernel_adios_variables,DEFINE_KERNEL_ADIOS_VARIABLES)(){}

void FC_FUNC_(write_kernels_crust_mantle_adios,WRITE_KERNELS_CRUST_MANTLE_ADIOS)(){}

void FC_FUNC_(write_kernels_outer_core_adios,WRITE_KERNELS_OUTER_CORE_ADIOS)(){}

void FC_FUNC_(write_kernels_inner_core_adios,WRITE_KERNELS_INNER_CORE_ADIOS)(){}

void FC_FUNC_(write_kernels_boundary_kl_adios,WRITE_KERNELS_BOUNDARY_KL_ADIOS)(){}

void FC_FUNC_(write_kernels_source_derivatives_adios,WRITE_KERNELS_SOURCE_DERIVATIVES_ADIOS)(){}

void FC_FUNC_(write_kernels_hessian_adios,WRITE_KERNELS_HESSIAN_ADIOS)(){}

// For xspecfem3d -- "ASDF" -- seismograms in ADIOS

void FC_FUNC_(init_asdf_data, INIT_ASDF_DATA)(void* asdf_event,
                                             int* total_seismos_local){}

void FC_FUNC_(store_asdf_data, STORE_ASDF_DATA)
    (void* my_asdf, realw* seismogram_tmp, int* irec_local, int *irec,
     char* chn, int* iorientation){}

void FC_FUNC_(close_asdf_data, CLOSE_ASDF_DATA)(void *my_asdf,
                                                int *total_seismos_local){}

void FC_FUNC_(write_asdf, WRITE_ASDF)(void* my_asdf){}


