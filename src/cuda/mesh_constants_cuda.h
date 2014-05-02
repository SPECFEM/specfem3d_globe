/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  2 . 0
 !               ---------------------------------------
 !
 !          Main authors: Dimitri Komatitsch and Jeroen Tromp
 !    Princeton University, USA and University of Pau / CNRS / INRIA
 ! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
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
 */

/* trivia

- for most working arrays we use now "realw" instead of "float" type declarations to make it easier to switch
  between a real or double precision simulation
  (matching CUSTOM_REAL == 4 or 8 in fortran routines).

- instead of boolean "logical" declared in fortran routines, in C (or Cuda-C) we have to use "int" variables.
  ifort / gfortran caveat:
    to check whether it is true or false, do not check for == 1 to test for true values since ifort just uses
    non-zero values for true (e.g. can be -1 for true). however, false will be always == 0.
  thus, rather use: if( var ) {...}  for testing if true instead of if( var == 1){...} (alternative: one could use if( var != 0 ){...}

*/

#ifndef CUDA_MESH_H
#define CUDA_MESH_H

#include <sys/types.h>
#include <unistd.h>

/* ----------------------------------------------------------------------------------------------- */

// for debugging and benchmarking

/* ----------------------------------------------------------------------------------------------- */

// debug: outputs traces
#define DEBUG 0
#if DEBUG == 1
#define TRACE(x) printf("%s\n",x);
#else
#define TRACE(x) // printf("%s\n",x);
#endif

// more outputs
#define MAXDEBUG 0
#if MAXDEBUG == 1
#define LOG(x) printf("%s\n",x)
#define PRINT5(var,offset) for(;print_count<5;print_count++) printf("var(%d)=%2.20f\n",print_count,var[offset+print_count]);
#define PRINT10(var) if(print_count<10) { printf("var=%1.20e\n",var); print_count++; }
#define PRINT10i(var) if(print_count<10) { printf("var=%d\n",var); print_count++; }
#else
#define LOG(x) // printf("%s\n",x);
#define PRINT5(var,offset) // for(i=0;i<10;i++) printf("var(%d)=%f\n",i,var[offset+i]);
#endif

// daniel debug: run backward simulations with/without GPU routines and empty arrays for debugging
#define DEBUG_BACKWARD_SIMULATIONS 0
#if DEBUG_BACKWARD_SIMULATIONS == 1
#define DEBUG_BACKWARD_ASSEMBLY() return;
#define DEBUG_BACKWARD_COUPLING() return;
#define DEBUG_BACKWARD_FORCES() return;
#define DEBUG_BACKWARD_KERNEL() return;
#define DEBUG_BACKWARD_SOURCES() return;
#define DEBUG_BACKWARD_TRANSFER() return;
#define DEBUG_BACKWARD_UPDATE() return;
#else
#define DEBUG_BACKWARD_ASSEMBLY()
#define DEBUG_BACKWARD_COUPLING()
#define DEBUG_BACKWARD_FORCES()
#define DEBUG_BACKWARD_KERNEL()
#define DEBUG_BACKWARD_SOURCES()
#define DEBUG_BACKWARD_TRANSFER()
#define DEBUG_BACKWARD_UPDATE()
#endif

// error checking after cuda function calls
// (note: this synchronizes many calls, thus e.g. no asynchronous memcpy possible)
//#define ENABLE_VERY_SLOW_ERROR_CHECKING

// maximum function
#define MAX(x,y)                    (((x) < (y)) ? (y) : (x))


/* ----------------------------------------------------------------------------------------------- */

// cuda constant arrays

/* ----------------------------------------------------------------------------------------------- */
// (must match constants.h definitions)

// dimensions
#define NDIM 3

// Gauss-Lobatto-Legendre
#define NGLLX 5
#define NGLL2 25
#define NGLL3 125 // no padding: requires same size as in fortran for NGLLX * NGLLY * NGLLZ

// padding: 128 == 2**7 might improve on older graphics cards w/ coalescent memory accesses:
#define NGLL3_PADDED 128
// no padding: 125 == 5*5*5 to avoid allocation of extra memory
//#define NGLL3_PADDED 125

// number of standard linear solids
#define N_SLS 3

// region ids
#define IREGION_CRUST_MANTLE  1
#define IREGION_OUTER_CORE    2
#define IREGION_INNER_CORE    3

// inner core : fictitious elements id (from constants.h)
#define IFLAG_IN_FICTITIOUS_CUBE  11

// R_EARTH_KM is the radius of the bottom of the oceans (radius of Earth in km)
#define R_EARTH_KM 6371.0f
// uncomment line below for PREM with oceans
//#define R_EARTH_KM 6368.0f

// Asynchronous memory copies between GPU and CPU
#define GPU_ASYNC_COPY 1

/* ----------------------------------------------------------------------------------------------- */

// (optional) pre-processing directive used in kernels: if defined check that it is also set in src/shared/constants.h:
// leads up to ~ 5% performance increase
//#define USE_MESH_COLORING_GPU
#ifdef USE_MESH_COLORING_GPU
#pragma message ("\nCompiling with: USE_MESH_COLORING_GPU enabled\n")
#endif

// note: mesh coloring has a tradeoff between extra costs for looping over colors
//          and slowliness of atomic updates;
//          in general, the more elements per color the more efficient
//
// thresholds for coloring :
//   we assume that the crust/mantle region has enough elements for coloring
//
// minimum number of elements required in inner core before mesh coloring becomes attractive
#define COLORING_MIN_NSPEC_INNER_CORE 1000
// minimum number of elements required in outer core before mesh coloring becomes attractive
#define COLORING_MIN_NSPEC_OUTER_CORE 1000

/* ----------------------------------------------------------------------------------------------- */

// Texture memory usage:
// requires CUDA version >= 4.0, see check below
// Use textures for d_displ and d_accel ~ 1% performance boost
//#define USE_TEXTURES_FIELDS

// Using texture memory for the hprime-style constants is slower on
// Fermi generation hardware, but *may* be faster on Kepler
// generation.
// Use textures for hprime_xx
//#define USE_TEXTURES_CONSTANTS

// CUDA version >= 4.0 needed for cudaTextureType1D and cudaDeviceSynchronize()
#if CUDA_VERSION < 4000
#undef USE_TEXTURES_FIELDS
#undef USE_TEXTURES_CONSTANTS
#endif

// compiling info
#ifdef USE_TEXTURES_FIELDS
#pragma message ("\nCompiling with: USE_TEXTURES_FIELDS enabled\n")
#endif
#ifdef USE_TEXTURES_CONSTANTS
#pragma message ("\nCompiling with: USE_TEXTURES_CONSTANTS enabled\n")
#endif

// (optional) unrolling loops
// leads up to ~1% performance increase
//#define MANUALLY_UNROLLED_LOOPS

// compiler specifications
// (optional) use launch_bounds specification to increase compiler optimization
//
// note: main kernel is Kernel_2_crust_mantle_impl() which is limited by register usage to only 5 active blocks
//       while shared memory usage would allow up to 7 blocks (see profiling with nvcc...)
//       here we specifiy to launch 7 blocks to increase occupancy and let the compiler reduce registers
//       (depending on GPU type, register spilling might slow down the performance)
//
// performance statistics: main kernel Kernel_2_crust_mantle_impl():
//       shared memory per block = 6200    for Kepler: total = 49152 -> limits active blocks to 7
//       registers per thread    = 72                                   (limited by LAUNCH_MIN_BLOCKS 7)
//       registers per block     = 9216                total = 65536    (limited by LAUNCH_MIN_BLOCKS 7)
//
// using launch_bounds leads to ~ 20% performance increase on Kepler GPUs
// (uncomment if not desired)
#define USE_LAUNCH_BOUNDS
#define LAUNCH_MIN_BLOCKS 7

/* ----------------------------------------------------------------------------------------------- */

// cuda kernel block size for updating displacements/potential (newmark time scheme)
// current hardware: 128 is slightly faster than 256 ( ~ 4%)
#define BLOCKSIZE_KERNEL1 128
#define BLOCKSIZE_KERNEL3 128
#define BLOCKSIZE_TRANSFER 256

// maximum grid dimension in one direction of GPU
#define MAXIMUM_GRID_DIM 65535

/* ----------------------------------------------------------------------------------------------- */

// indexing

#define INDEX2(isize,i,j) i + isize*j
#define INDEX3(isize,jsize,i,j,k) i + isize*(j + jsize*k)
#define INDEX4(isize,jsize,ksize,i,j,k,x) i + isize*(j + jsize*(k + ksize*x))
#define INDEX5(isize,jsize,ksize,xsize,i,j,k,x,y) i + isize*(j + jsize*(k + ksize*(x + xsize*y)))

//#define INDEX6(isize,jsize,ksize,xsize,ysize,i,j,k,x,y,z) i + isize*(j + jsize*(k + ksize*(x + xsize*(y + ysize*z))))
//#define INDEX4_PADDED(isize,jsize,ksize,i,j,k,x) i + isize*(j + jsize*k) + x*NGLL3_PADDED

/* ----------------------------------------------------------------------------------------------- */

// type of "working" variables: see also CUSTOM_REAL
// double precision temporary variables leads to 10% performance decrease
// in Kernel_2_impl (not very much..)
typedef float realw;

// textures
typedef texture<float, cudaTextureType1D, cudaReadModeElementType> realw_texture;

// restricted pointers: improves performance on Kepler ~ 10%
// see: http://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#restrict
typedef const float* __restrict__ realw_const_p; // otherwise use: //typedef const float* realw_const_p;
typedef float* __restrict__ realw_p; // otherwise use: //typedef float* realw_p;

/* ----------------------------------------------------------------------------------------------- */

// utility functions: defined in check_fields_cuda.cu

/* ----------------------------------------------------------------------------------------------- */

double get_time();
void get_free_memory(double* free_db, double* used_db, double* total_db);
void print_CUDA_error_if_any(cudaError_t err, int num);
void pause_for_debugger(int pause);
void exit_on_cuda_error(char* kernel_name);
void exit_on_error(char* info);
void synchronize_cuda();
void synchronize_mpi();
void start_timing_cuda(cudaEvent_t* start,cudaEvent_t* stop);
void stop_timing_cuda(cudaEvent_t* start,cudaEvent_t* stop, char* info_str);
void get_blocks_xy(int num_blocks,int* num_blocks_x,int* num_blocks_y);
realw get_device_array_maximum_value(realw* array,int size);

/* ----------------------------------------------------------------------------------------------- */

// mesh pointer wrapper structure

/* ----------------------------------------------------------------------------------------------- */

typedef struct mesh_ {

  // mesh resolution
  // ------------------------------------------------------------------ //
  // crust_mantle
  // ------------------------------------------------------------------ //
  int NSPEC_CRUST_MANTLE;
  int NGLOB_CRUST_MANTLE;
  int NSPEC_CRUST_MANTLE_STRAIN_ONLY;

  // interpolators
  realw* d_xix_crust_mantle; realw* d_xiy_crust_mantle; realw* d_xiz_crust_mantle;
  realw* d_etax_crust_mantle; realw* d_etay_crust_mantle; realw* d_etaz_crust_mantle;
  realw* d_gammax_crust_mantle; realw* d_gammay_crust_mantle; realw* d_gammaz_crust_mantle;

  // model parameters
  realw* d_rhostore_crust_mantle;
  realw* d_kappavstore_crust_mantle; realw* d_muvstore_crust_mantle;
  realw* d_kappahstore_crust_mantle; realw* d_muhstore_crust_mantle;
  realw* d_eta_anisostore_crust_mantle;

  // mass matrices
  realw* d_rmassx_crust_mantle;
  realw* d_rmassy_crust_mantle;
  realw* d_rmassz_crust_mantle;
  realw* d_b_rmassx_crust_mantle;
  realw* d_b_rmassy_crust_mantle;
  realw* d_b_rmassz_crust_mantle;

  // global indexing
  int* d_ibool_crust_mantle;
  int* d_ispec_is_tiso_crust_mantle;

  // mesh locations
  realw* d_xstore_crust_mantle; realw* d_ystore_crust_mantle; realw* d_zstore_crust_mantle;

  // anisotropic 3D mantle
  realw* d_c11store_crust_mantle;
  realw* d_c12store_crust_mantle;
  realw* d_c13store_crust_mantle;
  realw* d_c14store_crust_mantle;
  realw* d_c15store_crust_mantle;
  realw* d_c16store_crust_mantle;
  realw* d_c22store_crust_mantle;
  realw* d_c23store_crust_mantle;
  realw* d_c24store_crust_mantle;
  realw* d_c25store_crust_mantle;
  realw* d_c26store_crust_mantle;
  realw* d_c33store_crust_mantle;
  realw* d_c34store_crust_mantle;
  realw* d_c35store_crust_mantle;
  realw* d_c36store_crust_mantle;
  realw* d_c44store_crust_mantle;
  realw* d_c45store_crust_mantle;
  realw* d_c46store_crust_mantle;
  realw* d_c55store_crust_mantle;
  realw* d_c56store_crust_mantle;
  realw* d_c66store_crust_mantle;

  // wavefields
  // displacement, velocity, acceleration
  realw* d_displ_crust_mantle; realw* d_veloc_crust_mantle; realw* d_accel_crust_mantle;
  // backward/reconstructed elastic wavefield
  realw* d_b_displ_crust_mantle; realw* d_b_veloc_crust_mantle; realw* d_b_accel_crust_mantle;

  // attenuation
  realw* d_R_xx_crust_mantle;
  realw* d_R_yy_crust_mantle;
  realw* d_R_xy_crust_mantle;
  realw* d_R_xz_crust_mantle;
  realw* d_R_yz_crust_mantle;

  realw* d_b_R_xx_crust_mantle;
  realw* d_b_R_yy_crust_mantle;
  realw* d_b_R_xy_crust_mantle;
  realw* d_b_R_xz_crust_mantle;
  realw* d_b_R_yz_crust_mantle;

  realw* d_factor_common_crust_mantle;
  realw* d_one_minus_sum_beta_crust_mantle;

  realw* d_epsilondev_xx_crust_mantle;
  realw* d_epsilondev_yy_crust_mantle;
  realw* d_epsilondev_xy_crust_mantle;
  realw* d_epsilondev_xz_crust_mantle;
  realw* d_epsilondev_yz_crust_mantle;

  realw* d_b_epsilondev_xx_crust_mantle;
  realw* d_b_epsilondev_yy_crust_mantle;
  realw* d_b_epsilondev_xy_crust_mantle;
  realw* d_b_epsilondev_xz_crust_mantle;
  realw* d_b_epsilondev_yz_crust_mantle;

  realw* d_eps_trace_over_3_crust_mantle;
  realw* d_b_eps_trace_over_3_crust_mantle;

  // kernels
  realw* d_rho_kl_crust_mantle;
  realw* d_alpha_kl_crust_mantle;
  realw* d_beta_kl_crust_mantle;
  realw* d_cijkl_kl_crust_mantle;
  realw* d_hess_kl_crust_mantle;

  // inner / outer elements
  int* d_phase_ispec_inner_crust_mantle;
  int num_phase_ispec_crust_mantle;

  int nspec_outer_crust_mantle;
  int nspec_inner_crust_mantle;
  int nspec2D_bottom_crust_mantle;

  int num_colors_inner_crust_mantle;
  int num_colors_outer_crust_mantle;
  int* h_num_elem_colors_crust_mantle;

  int* d_ibelm_bottom_crust_mantle;

  // ------------------------------------------------------------------ //
  // outer_core
  // ------------------------------------------------------------------ //
  int NSPEC_OUTER_CORE;
  int NGLOB_OUTER_CORE;

  // interpolators
  realw* d_xix_outer_core; realw* d_xiy_outer_core; realw* d_xiz_outer_core;
  realw* d_etax_outer_core; realw* d_etay_outer_core; realw* d_etaz_outer_core;
  realw* d_gammax_outer_core; realw* d_gammay_outer_core; realw* d_gammaz_outer_core;

  // model parameters
  realw* d_rhostore_outer_core; realw* d_kappavstore_outer_core;

  realw* d_rmass_outer_core;
  realw* d_b_rmass_outer_core;

  // global indexing
  int* d_ibool_outer_core;

  // mesh locations
  realw* d_xstore_outer_core; realw* d_ystore_outer_core; realw* d_zstore_outer_core;

  // wavefields
  // displacement, velocity, acceleration
  realw* d_displ_outer_core; realw* d_veloc_outer_core; realw* d_accel_outer_core;
  // backward/reconstructed elastic wavefield
  realw* d_b_displ_outer_core; realw* d_b_veloc_outer_core; realw* d_b_accel_outer_core;

  // kernels
  realw* d_rho_kl_outer_core;
  realw* d_alpha_kl_outer_core;

  // inner / outer elements
  int* d_phase_ispec_inner_outer_core;
  int num_phase_ispec_outer_core;

  int nspec_outer_outer_core;
  int nspec_inner_outer_core;
  int nspec2D_top_outer_core;
  int nspec2D_bottom_outer_core;

  int num_colors_inner_outer_core;
  int num_colors_outer_outer_core;
  int* h_num_elem_colors_outer_core;

  int* d_ibelm_top_outer_core;
  int* d_ibelm_bottom_outer_core;

  // normals definitions for coupling regions
  realw* d_normal_top_outer_core;
  realw* d_normal_bottom_outer_core;

  // Jacobian definitions
  realw* d_jacobian2D_top_outer_core;
  realw* d_jacobian2D_bottom_outer_core;

  // ------------------------------------------------------------------ //
  // inner_core
  // ------------------------------------------------------------------ //
  int NSPEC_INNER_CORE;
  int NGLOB_INNER_CORE;
  int NSPEC_INNER_CORE_STRAIN_ONLY;

  // interpolators
  realw* d_xix_inner_core; realw* d_xiy_inner_core; realw* d_xiz_inner_core;
  realw* d_etax_inner_core; realw* d_etay_inner_core; realw* d_etaz_inner_core;
  realw* d_gammax_inner_core; realw* d_gammay_inner_core; realw* d_gammaz_inner_core;

  // model parameters
  realw* d_rhostore_inner_core;
  realw* d_kappavstore_inner_core; realw* d_muvstore_inner_core;

  realw* d_rmassx_inner_core;
  realw* d_rmassy_inner_core;
  realw* d_rmassz_inner_core;
  realw* d_b_rmassx_inner_core;
  realw* d_b_rmassy_inner_core;
  realw* d_b_rmassz_inner_core;

  // global indexing
  int* d_ibool_inner_core;
  int* d_idoubling_inner_core;

  // mesh locations
  realw* d_xstore_inner_core; realw* d_ystore_inner_core; realw* d_zstore_inner_core;

  // anisotropic 3D mantle
  realw* d_c11store_inner_core;
  realw* d_c12store_inner_core;
  realw* d_c13store_inner_core;
  realw* d_c33store_inner_core;
  realw* d_c44store_inner_core;

  // wavefields
  // displacement, velocity, acceleration
  realw* d_displ_inner_core; realw* d_veloc_inner_core; realw* d_accel_inner_core;
  // backward/reconstructed elastic wavefield
  realw* d_b_displ_inner_core; realw* d_b_veloc_inner_core; realw* d_b_accel_inner_core;

  // attenuation
  realw* d_R_xx_inner_core;
  realw* d_R_yy_inner_core;
  realw* d_R_xy_inner_core;
  realw* d_R_xz_inner_core;
  realw* d_R_yz_inner_core;

  realw* d_b_R_xx_inner_core;
  realw* d_b_R_yy_inner_core;
  realw* d_b_R_xy_inner_core;
  realw* d_b_R_xz_inner_core;
  realw* d_b_R_yz_inner_core;

  realw* d_factor_common_inner_core;
  realw* d_one_minus_sum_beta_inner_core;

  realw* d_epsilondev_xx_inner_core;
  realw* d_epsilondev_yy_inner_core;
  realw* d_epsilondev_xy_inner_core;
  realw* d_epsilondev_xz_inner_core;
  realw* d_epsilondev_yz_inner_core;

  realw* d_b_epsilondev_xx_inner_core;
  realw* d_b_epsilondev_yy_inner_core;
  realw* d_b_epsilondev_xy_inner_core;
  realw* d_b_epsilondev_xz_inner_core;
  realw* d_b_epsilondev_yz_inner_core;

  realw* d_eps_trace_over_3_inner_core;
  realw* d_b_eps_trace_over_3_inner_core;

  // kernels
  realw* d_rho_kl_inner_core;
  realw* d_alpha_kl_inner_core;
  realw* d_beta_kl_inner_core;

  // inner / outer elements
  int* d_phase_ispec_inner_inner_core;
  int num_phase_ispec_inner_core;

  int nspec_outer_inner_core;
  int nspec_inner_inner_core;
  int nspec2D_top_inner_core;

  int num_colors_inner_inner_core;
  int num_colors_outer_inner_core;
  int* h_num_elem_colors_inner_core;

  int* d_ibelm_top_inner_core;

  // ------------------------------------------------------------------ //
  // oceans
  // ------------------------------------------------------------------ //
  int npoin_oceans;

  // model parameter
  int* d_ibool_ocean_load;
  realw* d_rmass_ocean_load;
  realw* d_normal_ocean_load;

  // ------------------------------------------------------------------ //
  // attenuation
  // ------------------------------------------------------------------ //
  realw* d_alphaval;
  realw* d_betaval;
  realw* d_gammaval;

  realw* d_b_alphaval;
  realw* d_b_betaval;
  realw* d_b_gammaval;

  // ------------------------------------------------------------------ //
  // GLL points & weights
  // ------------------------------------------------------------------ //

  // pointers to constant memory arrays
  realw* d_hprime_xx;
  //realw* d_hprime_yy; // only needed if NGLLX != NGLLY != NGLLZ
  //realw* d_hprime_zz; // only needed if NGLLX != NGLLY != NGLLZ

  realw* d_hprimewgll_xx;
  //realw* d_hprimewgll_yy; // only needed if NGLLX != NGLLY != NGLLZ
  //realw* d_hprimewgll_zz; // only needed if NGLLX != NGLLY != NGLLZ

  realw* d_wgllwgll_xy; realw* d_wgllwgll_xz; realw* d_wgllwgll_yz;
  realw* d_wgll_cube;

  // simulation type: 1 = forward, 2 = adjoint, 3 = kernel
  int simulation_type;

  // mesh coloring flag
  int use_mesh_coloring_gpu;

  // simulation flags
  int save_forward;
  int absorbing_conditions;

  int attenuation;
  int undo_attenuation;
  int partial_phys_dispersion_only;
  int use_3d_attenuation_arrays;

  int compute_and_store_strain;
  int anisotropic_3D_mantle;
  int gravity;
  int rotation;
  int exact_mass_matrix_for_rotation;
  int oceans;
  int anisotropic_inner_core;
  int save_boundary_mesh;

  int anisotropic_kl;
  int approximate_hess_kl;

  realw deltat;
  realw b_deltat;

  // ------------------------------------------------------------------ //
  // gravity
  // ------------------------------------------------------------------ //
  realw* d_d_ln_density_dr_table; // needed also for no gravity case
  realw* d_minus_rho_g_over_kappa_fluid;
  realw* d_minus_gravity_table;
  realw* d_minus_deriv_gravity_table;
  realw* d_density_table;

  realw minus_g_icb;
  realw minus_g_cmb;

  realw RHO_BOTTOM_OC;
  realw RHO_TOP_OC;

  // ------------------------------------------------------------------ //
  // rotation
  // ------------------------------------------------------------------ //
  realw two_omega_earth;
  realw* d_A_array_rotation; realw* d_B_array_rotation;

  // needed for backward/reconstructed fields (kernel runs)
  realw b_two_omega_earth;
  realw* d_b_A_array_rotation; realw* d_b_B_array_rotation;

  // ------------------------------------------------------------------ //
  // sources
  // ------------------------------------------------------------------ //
  int nsources_local;
  realw* d_sourcearrays;
  double* d_stf_pre_compute;
  int* d_islice_selected_source;
  int* d_ispec_selected_source;

  // ------------------------------------------------------------------ //
  // receivers
  // ------------------------------------------------------------------ //
  int* d_number_receiver_global;
  int* d_ispec_selected_rec;
  int* d_islice_selected_rec;

  int nrec_local;

  realw* d_station_seismo_field;
  realw* h_station_seismo_field;

  realw* d_station_strain_field;
  realw* h_station_strain_field;

  // adjoint receivers/sources
  int nadj_rec_local;
  realw* d_adj_sourcearrays;
  realw* h_adj_sourcearrays_slice;
  int* d_pre_computed_irec;

  // ------------------------------------------------------------------ //
  // assembly
  // ------------------------------------------------------------------ //
  int myrank;

  int num_interfaces_crust_mantle;
  int max_nibool_interfaces_cm;
  int* d_nibool_interfaces_crust_mantle;
  int* d_ibool_interfaces_crust_mantle;

  realw* d_send_accel_buffer_crust_mantle;
  realw* d_b_send_accel_buffer_crust_mantle;

  int num_interfaces_inner_core;
  int max_nibool_interfaces_ic;
  int* d_nibool_interfaces_inner_core;
  int* d_ibool_interfaces_inner_core;

  realw* d_send_accel_buffer_inner_core;
  realw* d_b_send_accel_buffer_inner_core;

  int num_interfaces_outer_core;
  int max_nibool_interfaces_oc;
  int* d_nibool_interfaces_outer_core;
  int* d_ibool_interfaces_outer_core;

  realw* d_send_accel_buffer_outer_core;
  realw* d_b_send_accel_buffer_outer_core;

  // ------------------------------------------------------------------ //
  // absorbing boundaries
  // ------------------------------------------------------------------ //

  int nspec2D_xmin_crust_mantle,nspec2D_xmax_crust_mantle;
  int nspec2D_ymin_crust_mantle,nspec2D_ymax_crust_mantle;

  int* d_nimin_crust_mantle, *d_nimax_crust_mantle;
  int* d_njmin_crust_mantle, *d_njmax_crust_mantle;
  int* d_nkmin_xi_crust_mantle, *d_nkmin_eta_crust_mantle;

  int* d_ibelm_xmin_crust_mantle, *d_ibelm_xmax_crust_mantle;
  int* d_ibelm_ymin_crust_mantle, *d_ibelm_ymax_crust_mantle;

  realw* d_normal_xmin_crust_mantle, *d_normal_xmax_crust_mantle;
  realw* d_normal_ymin_crust_mantle, *d_normal_ymax_crust_mantle;

  realw* d_jacobian2D_xmin_crust_mantle, *d_jacobian2D_xmax_crust_mantle;
  realw* d_jacobian2D_ymin_crust_mantle, *d_jacobian2D_ymax_crust_mantle;

  realw* d_absorb_xmin_crust_mantle, *d_absorb_xmax_crust_mantle;
  realw* d_absorb_ymin_crust_mantle, *d_absorb_ymax_crust_mantle;

  realw* d_rho_vp_crust_mantle;
  realw* d_rho_vs_crust_mantle;

  int nspec2D_xmin_outer_core,nspec2D_xmax_outer_core;
  int nspec2D_ymin_outer_core,nspec2D_ymax_outer_core;
  int nspec2D_zmin_outer_core;

  int* d_nimin_outer_core, *d_nimax_outer_core;
  int* d_njmin_outer_core, *d_njmax_outer_core;
  int* d_nkmin_xi_outer_core, *d_nkmin_eta_outer_core;

  int* d_ibelm_xmin_outer_core, *d_ibelm_xmax_outer_core;
  int* d_ibelm_ymin_outer_core, *d_ibelm_ymax_outer_core;

  realw* d_jacobian2D_xmin_outer_core, *d_jacobian2D_xmax_outer_core;
  realw* d_jacobian2D_ymin_outer_core, *d_jacobian2D_ymax_outer_core;

  realw* d_absorb_xmin_outer_core, *d_absorb_xmax_outer_core;
  realw* d_absorb_ymin_outer_core, *d_absorb_ymax_outer_core;
  realw* d_absorb_zmin_outer_core;

  realw* d_vp_outer_core;

  // ------------------------------------------------------------------ //
  // noise tomography
  // ------------------------------------------------------------------ //
  int noise_tomography;

  realw* d_noise_surface_movie;
  realw* d_noise_sourcearray;

  realw* d_normal_x_noise;
  realw* d_normal_y_noise;
  realw* d_normal_z_noise;
  realw* d_mask_noise;

  // free surface
  int nspec2D_top_crust_mantle;
  int* d_ibelm_top_crust_mantle;
  realw* d_jacobian2D_top_crust_mantle;

  // noise sensitivity kernel
  realw* d_Sigma_kl;

  // ------------------------------------------------------------------ //
  // optimizations
  // ------------------------------------------------------------------ //

  // overlapped memcpy streams
  cudaStream_t compute_stream;
  cudaStream_t copy_stream;

  // A buffer for MPI send/recv, which is duplicated in Fortran but is
  // allocated with pinned memory to facilitate asynchronous device <->
  // host memory transfers
  // crust/mantle
  float* h_send_accel_buffer_cm;
  float* h_recv_accel_buffer_cm;
  float* h_b_send_accel_buffer_cm;
  float* h_b_recv_accel_buffer_cm;
  // inner core
  float* h_send_accel_buffer_ic;
  float* h_recv_accel_buffer_ic;
  float* h_b_send_accel_buffer_ic;
  float* h_b_recv_accel_buffer_ic;
  // outer core
  float* h_send_accel_buffer_oc;
  float* h_recv_accel_buffer_oc;
  float* h_b_send_accel_buffer_oc;
  float* h_b_recv_accel_buffer_oc;

} Mesh;


#endif
// CUDA_MESH_H
