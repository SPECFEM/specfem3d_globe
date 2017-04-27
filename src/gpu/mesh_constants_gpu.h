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
   (matching CUSTOM_REAL == 4 or 8 in Fortran routines).

   - instead of boolean "logical" declared in Fortran routines, in C (or OpenCL-C) we have to use "int" variables.
   ifort / gfortran caveat:
   to check whether it is true or false, do not check for == 1 to test for true values since ifort just uses
   non-zero values for true (e.g. can be -1 for true). however, false will be always == 0.
   thus, rather use: if (var) {...}  for testing if true instead of if (var == 1) {...} (alternative: one could use if (var != 0) {...}
*/

#ifndef MESH_CONSTANTS_GPU_H
#define MESH_CONSTANTS_GPU_H

#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"

#ifdef USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif

// type of "working" variables: see also CUSTOM_REAL
// double precision temporary variables leads to 10% performance decrease
// in Kernel_2_impl (not very much..)
typedef float realw;

/*----------------------------------------------------------------------------------------------- */
// for debugging and benchmarking
/*----------------------------------------------------------------------------------------------- */

// debug: outputs traces
#define DEBUG 0
#if DEBUG == 1
#pragma message ("\nCompiling with: DEBUG enabled\n")
#define TRACE(x) printf ("%s\n", x); fflush(stdout);
#define TRACE_EXTENDED(x) printf ("%s --- function %s file %s line %d \n",x,__func__,__FILE__,__LINE__); fflush(stdout);
#else
#define TRACE(x)
#endif

// debug: outputs maximum values of wavefields
#define DEBUG_FIELDS 0

// more outputs
#define MAXDEBUG 0
#if MAXDEBUG == 1
#pragma message ("\nCompiling with: MAXDEBUG enabled\n")
#define LOG(x) printf ("%s\n", x)
#define PRINT5(var, offset) for (;print_count<5;print_count++) printf ("var (%d)=%2.20f\n", print_count, var[offset+print_count]);
#define PRINT10(var) if (print_count<10) { printf ("var=%1.20e\n", var); print_count++; }
#define PRINT10i(var) if (print_count<10) { printf ("var=%d\n", var); print_count++; }
#else
#define LOG(x)   // printf ("%s\n", x);
#define PRINT5(var, offset)   // for (i = 0;i<10;i++) printf ("var (%d)=%f\n", i, var[offset+i]);
#endif

// debug: run backward simulations with/without GPU routines and empty arrays for debugging
#define DEBUG_BACKWARD_SIMULATIONS 0
#if DEBUG_BACKWARD_SIMULATIONS == 1
#pragma message ("\nCompiling with: DEBUG_BACKWARD_SIMULATIONS enabled\n")
#define DEBUG_BACKWARD_ASSEMBLY_OC()  return;
#define DEBUG_BACKWARD_ASSEMBLY_IC()  return;
#define DEBUG_BACKWARD_ASSEMBLY_CM()  return;
#define DEBUG_BACKWARD_COUPLING()     return;
#define DEBUG_BACKWARD_FORCES()       return;
#define DEBUG_BACKWARD_KERNEL()       return;
#define DEBUG_BACKWARD_SOURCES()      return;
#define DEBUG_BACKWARD_TRANSFER()     return;
#define DEBUG_BACKWARD_UPDATE()       return;
#else
#define DEBUG_BACKWARD_ASSEMBLY_OC()
#define DEBUG_BACKWARD_ASSEMBLY_IC()
#define DEBUG_BACKWARD_ASSEMBLY_CM()
#define DEBUG_BACKWARD_COUPLING()
#define DEBUG_BACKWARD_FORCES()
#define DEBUG_BACKWARD_KERNEL()
#define DEBUG_BACKWARD_SOURCES()
#define DEBUG_BACKWARD_TRANSFER()
#define DEBUG_BACKWARD_UPDATE()
#endif

// error checking after cuda function calls
// (note: this synchronizes many calls, thus e.g. no asynchronous memcpy possible)
#define ENABLE_VERY_SLOW_ERROR_CHECKING 0
#if ENABLE_VERY_SLOW_ERROR_CHECKING == 1
#pragma message ("\nCompiling with: ENABLE_VERY_SLOW_ERROR_CHECKING enabled\n")
#define GPU_ERROR_CHECKING(x) exit_on_gpu_error(x);
#else
#define GPU_ERROR_CHECKING(x)
#endif

// maximum function
#define MAX(x, y)                  (((x) < (y)) ? (y) : (x))

/*----------------------------------------------------------------------------------------------- */
// GPU constant arrays
/*----------------------------------------------------------------------------------------------- */
// (must match constants.h definitions)

// dimensions
#define NDIM 3

// Gauss-Lobatto-Legendre
#define NGLLX 5
#define NGLL2 25
#define NGLL3 125 // no padding: requires same size as in Fortran for NGLLX * NGLLY * NGLLZ

// padding: 128 == 2**7 might improve on older graphics cards w/ coalescent memory accesses:
#define NGLL3_PADDED 128
// no padding: 125 == 5*5*5 to avoid allocation of extra memory
//#define NGLL3_PADDED 125

// number of standard linear solids
#define N_SLS 3

// region ids
#define IREGION_CRUST_MANTLE  1
#define IREGION_OUTER_CORE    2
#define IREGION_INNER_CORE  3

// inner core : fictitious elements id (from constants.h)
#define IFLAG_IN_FICTITIOUS_CUBE  11

// R_EARTH_KM is the radius of the bottom of the oceans (radius of Earth in km)
#define R_EARTH_KM 6371.0f

// uncomment line below for PREM with oceans
//#define R_EARTH_KM 6368.0f

// Asynchronous memory copies between GPU and CPU
// (set to 0 for synchronuous/blocking copies, set to 1 for asynchronuous copies)
#define GPU_ASYNC_COPY 1

// Reduce GPU-register pressure by limited the number of thread spread
// (GPU for embedded devices are not powerful enough for big kernels)
// Must match BOAST compiled value (--elem flag)
#ifndef GPU_ELEM_PER_THREAD
#define GPU_ELEM_PER_THREAD 1
#endif

/*----------------------------------------------------------------------------------------------- */

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

/*----------------------------------------------------------------------------------------------- */

// Texture memory usage:
// requires CUDA version >= 4.0, see check below
// Use textures for d_displ and d_accel ~ 1% performance boost
//#define USE_TEXTURES_FIELDS

// Using texture memory for the hprime-style constants is slower on
// Fermi generation hardware, but *may* be faster on Kepler
// generation.
// Use textures for hprime_xx
//#define USE_TEXTURES_CONSTANTS

#ifdef USE_CUDA
  // CUDA version >= 4.0 needed for cudaTextureType1D and cudaDeviceSynchronize()
  #if CUDA_VERSION < 4000
    #undef USE_TEXTURES_FIELDS
    #undef USE_TEXTURES_CONSTANTS
    #pragma message ("\nCompiling for CUDA version < 4.0\n")
  #endif
#endif

// compiling info
#ifdef USE_TEXTURES_FIELDS
#pragma message ("\nCompiling with: USE_TEXTURES_FIELDS enabled\n")
#endif
#ifdef USE_TEXTURES_CONSTANTS
#pragma message ("\nCompiling with: USE_TEXTURES_CONSTANTS enabled\n")
#endif

// (optional) unrolling loops
// leads up to ~10% performance increase in OpenCL and ~1% in Cuda
#define MANUALLY_UNROLLED_LOOPS

// CUDA compiler specifications
// (optional) use launch_bounds specification to increase compiler optimization
//
#ifdef GPU_DEVICE_K20
// note: main kernel is Kernel_2_crust_mantle_impl() which is limited by register usage to only 5 active blocks
//       while shared memory usage would allow up to 7 blocks (see profiling with nvcc...)
//       here we specifiy to launch 7 blocks to increase occupancy and let the compiler reduce registers
//       (depending on GPU type, register spilling might slow down the performance)
//       (single block uses 128 threads -> ptxas info: 72 registers (per thread) -> 72 * 128 = 9216 registers per block)
//
// performance statistics: main kernel Kernel_2_crust_mantle_impl():
//       shared memory per block = 6200    for Kepler: total = 49152 -> limits active blocks to 7
//       registers per thread    = 72                                   (limited by LAUNCH_MIN_BLOCKS 7)
//       registers per block     = 9216                total = 65536    (limited by LAUNCH_MIN_BLOCKS 7)
//
// using launch_bounds leads to ~ 20% performance increase on Kepler GPUs
// (uncomment if not desired)
//#pragma message ("\nCompiling with: USE_LAUNCH_BOUNDS enabled for K20\n")
#define USE_LAUNCH_BOUNDS
#define LAUNCH_MIN_BLOCKS 7
#endif
#ifdef GPU_DEVICE_P100
// Pascal P100: by default, the crust_mantle_impl_kernel_forward kernel uses 80 registers.
//              80 * 128 threads -> 10240 registers    for Pascal: total of 65536 -> limits active blocks to 6
//              using launch bounds to increase the number of blocks will lead to register spilling.
//              for Pascal, the spilling slows down the kernels by ~6%
#undef USE_LAUNCH_BOUNDS
#define LAUNCH_MIN_BLOCKS 6
#endif



/*----------------------------------------------------------------------------------------------- */

// kernel block size for updating displacements/potential (Newmark time scheme)
// current hardware: 128 is slightly faster than 256 (~ 4%)
#define BLOCKSIZE_KERNEL1 128
#define BLOCKSIZE_KERNEL3 128
#define BLOCKSIZE_TRANSFER 256

// maximum grid dimension in one direction of GPU
#define MAXIMUM_GRID_DIM 65535

/*----------------------------------------------------------------------------------------------- */

//MIC (Knights Corner)
#define MIC_KNIGHTS_CORNER 0
#if MIC_KNIGHTS_CORNER == 1
// redefines block sizes
// local work group size: SIMD size 512-bit/64-byte -> multiple of 16 floats (each 4-byte)
#undef BLOCKSIZE_KERNEL1
#undef BLOCKSIZE_KERNEL3
#undef BLOCKSIZE_TRANSFER
#define BLOCKSIZE_KERNEL1 64    // for updating displacement (Newmark)
#define BLOCKSIZE_KERNEL3 64
#define BLOCKSIZE_TRANSFER 128

// max_compute_unit (number of work-groups should be equal to this for best performance)
// Knights Corner has 236 threads (1percore), optimal number of work groups is a multiple
#define MAXIMUM_COMPUTE_UNITS 236

// preferred work group size multiple for kernels (taken from device/kernel info)
#define PREFERRED_WORK_GROUP_SIZE_MULTIPLE 128

// max_work_group_size
// device info: CL_DEVICE_MAX_WORK_ITEM_SIZES == 8192
// kernel info: CL_KERNEL_WORK_GROUP_SIZE == 2048 for crust_mantle kernel
#undef MAXIMUM_GRID_DIM
#define MAXIMUM_GRID_DIM 32 * MAXIMUM_COMPUTE_UNITS // 7552 = 32 * 236 = 7552
#endif // MIC_KNIGHTS_CORNER

/*----------------------------------------------------------------------------------------------- */

// balancing work group x/y-size
#undef BALANCE_WORK_GROUP

// maximum number of work group units in one dimension
#define BALANCE_WORK_GROUP_UNITS 7552 // == 32 * 236 for Knights Corner test

#ifdef BALANCE_WORK_GROUP
#pragma message ("\nCompiling with: BALANCE_WORK_GROUP enabled\n")
#endif

/*----------------------------------------------------------------------------------------------- */

// indexing
#define INDEX2(isize,i,j) i + isize*j
#define INDEX3(isize,jsize,i,j,k) i + isize*(j + jsize*k)
#define INDEX4(isize,jsize,ksize,i,j,k,x) i + isize*(j + jsize*(k + ksize*x))
#define INDEX5(isize,jsize,ksize,xsize,i,j,k,x,y) i + isize*(j + jsize*(k + ksize*(x + xsize*y)))

/*----------------------------------------------------------------------------------------------- */

#ifdef USE_OPENCL
#include "mesh_constants_ocl.h"
#endif
#ifdef USE_CUDA
#include "mesh_constants_cuda.h"
#endif

typedef union {
#ifdef USE_OPENCL
  cl_mem ocl;
#endif
#ifdef USE_CUDA
  int *cuda;
#endif
} gpu_int_mem;

typedef union {
#ifdef USE_OPENCL
  cl_mem ocl;
#endif
#ifdef USE_CUDA
  realw *cuda;
#endif
} gpu_realw_mem;

typedef union {
#ifdef USE_OPENCL
  cl_mem ocl;
#endif
#ifdef USE_CUDA
  double *cuda;
#endif
} gpu_double_mem;

typedef union {
#ifdef USE_OPENCL
  cl_mem ocl;
#endif
#ifdef USE_CUDA
  void *cuda;
#endif
} gpu_mem;


#ifdef __cplusplus
#define EXTERN_LANG "C"
#else
#define EXTERN_LANG
#endif

extern int run_cuda;
extern int run_opencl;

/*----------------------------------------------------------------------------------------------- */
// mesh pointer wrapper structure
/*----------------------------------------------------------------------------------------------- */

typedef struct mesh_ {

  // mesh resolution
  // ------------------------------------------------------------------   //
  // crust_mantle
  // ------------------------------------------------------------------   //
  int NSPEC_CRUST_MANTLE;
  int NGLOB_CRUST_MANTLE;
  int NSPEC_CRUST_MANTLE_STRAIN_ONLY;
  int NSPECMAX_TISO_MANTLE;
  int NSPECMAX_ISO_MANTLE;

  // interpolators
  gpu_realw_mem d_xix_crust_mantle;
  gpu_realw_mem d_xiy_crust_mantle;
  gpu_realw_mem d_xiz_crust_mantle;
  gpu_realw_mem d_etax_crust_mantle;
  gpu_realw_mem d_etay_crust_mantle;
  gpu_realw_mem d_etaz_crust_mantle;
  gpu_realw_mem d_gammax_crust_mantle;
  gpu_realw_mem d_gammay_crust_mantle;
  gpu_realw_mem d_gammaz_crust_mantle;

  // model parameters
  gpu_realw_mem d_rhostore_crust_mantle;
  gpu_realw_mem d_kappavstore_crust_mantle;
  gpu_realw_mem d_muvstore_crust_mantle;
  gpu_realw_mem d_kappahstore_crust_mantle;
  gpu_realw_mem d_muhstore_crust_mantle;
  gpu_realw_mem d_eta_anisostore_crust_mantle;

  // mass matrices
  gpu_realw_mem d_rmassx_crust_mantle;
  gpu_realw_mem d_rmassy_crust_mantle;
  gpu_realw_mem d_rmassz_crust_mantle;
  gpu_realw_mem d_b_rmassx_crust_mantle;
  gpu_realw_mem d_b_rmassy_crust_mantle;
  gpu_realw_mem d_b_rmassz_crust_mantle;

  // global indexing
  gpu_int_mem d_ibool_crust_mantle;
  gpu_int_mem d_ispec_is_tiso_crust_mantle;

  // mesh locations
  gpu_realw_mem d_rstore_crust_mantle;

  // anisotropic 3D mantle
  gpu_realw_mem d_c11store_crust_mantle;
  gpu_realw_mem d_c12store_crust_mantle;
  gpu_realw_mem d_c13store_crust_mantle;
  gpu_realw_mem d_c14store_crust_mantle;
  gpu_realw_mem d_c15store_crust_mantle;
  gpu_realw_mem d_c16store_crust_mantle;
  gpu_realw_mem d_c22store_crust_mantle;
  gpu_realw_mem d_c23store_crust_mantle;
  gpu_realw_mem d_c24store_crust_mantle;
  gpu_realw_mem d_c25store_crust_mantle;
  gpu_realw_mem d_c26store_crust_mantle;
  gpu_realw_mem d_c33store_crust_mantle;
  gpu_realw_mem d_c34store_crust_mantle;
  gpu_realw_mem d_c35store_crust_mantle;
  gpu_realw_mem d_c36store_crust_mantle;
  gpu_realw_mem d_c44store_crust_mantle;
  gpu_realw_mem d_c45store_crust_mantle;
  gpu_realw_mem d_c46store_crust_mantle;
  gpu_realw_mem d_c55store_crust_mantle;
  gpu_realw_mem d_c56store_crust_mantle;
  gpu_realw_mem d_c66store_crust_mantle;

  // wavefields
  // displacement, velocity, acceleration
  gpu_realw_mem d_displ_crust_mantle;
  gpu_realw_mem d_veloc_crust_mantle;
  gpu_realw_mem d_accel_crust_mantle;
  // backward/reconstructed elastic wavefield
  gpu_realw_mem d_b_displ_crust_mantle;
  gpu_realw_mem d_b_veloc_crust_mantle;
  gpu_realw_mem d_b_accel_crust_mantle;

  // attenuation
  gpu_realw_mem d_R_xx_crust_mantle;
  gpu_realw_mem d_R_yy_crust_mantle;
  gpu_realw_mem d_R_xy_crust_mantle;
  gpu_realw_mem d_R_xz_crust_mantle;
  gpu_realw_mem d_R_yz_crust_mantle;

  gpu_realw_mem d_b_R_xx_crust_mantle;
  gpu_realw_mem d_b_R_yy_crust_mantle;
  gpu_realw_mem d_b_R_xy_crust_mantle;
  gpu_realw_mem d_b_R_xz_crust_mantle;
  gpu_realw_mem d_b_R_yz_crust_mantle;

  gpu_realw_mem d_factor_common_crust_mantle;
  gpu_realw_mem d_one_minus_sum_beta_crust_mantle;

  gpu_realw_mem d_epsilondev_xx_crust_mantle;
  gpu_realw_mem d_epsilondev_yy_crust_mantle;
  gpu_realw_mem d_epsilondev_xy_crust_mantle;
  gpu_realw_mem d_epsilondev_xz_crust_mantle;
  gpu_realw_mem d_epsilondev_yz_crust_mantle;

  gpu_realw_mem d_b_epsilondev_xx_crust_mantle;
  gpu_realw_mem d_b_epsilondev_yy_crust_mantle;
  gpu_realw_mem d_b_epsilondev_xy_crust_mantle;
  gpu_realw_mem d_b_epsilondev_xz_crust_mantle;
  gpu_realw_mem d_b_epsilondev_yz_crust_mantle;

  gpu_realw_mem d_eps_trace_over_3_crust_mantle;
  gpu_realw_mem d_b_eps_trace_over_3_crust_mantle;

  // kernels
  gpu_realw_mem d_rho_kl_crust_mantle;
  gpu_realw_mem d_alpha_kl_crust_mantle;
  gpu_realw_mem d_beta_kl_crust_mantle;
  gpu_realw_mem d_cijkl_kl_crust_mantle;
  gpu_realw_mem d_hess_kl_crust_mantle;

  // inner / outer elements
  gpu_int_mem d_phase_ispec_inner_crust_mantle;
  int num_phase_ispec_crust_mantle;

  int nspec_outer_crust_mantle;
  int nspec_inner_crust_mantle;
  int nspec2D_bottom_crust_mantle;

  int num_colors_inner_crust_mantle;
  int num_colors_outer_crust_mantle;
  int *h_num_elem_colors_crust_mantle;

  gpu_int_mem d_ibelm_bottom_crust_mantle;

  // ------------------------------------------------------------------   //
  // outer_core
  // ------------------------------------------------------------------   //
  int NSPEC_OUTER_CORE;
  int NGLOB_OUTER_CORE;

  // interpolators
  gpu_realw_mem d_xix_outer_core;
  gpu_realw_mem d_xiy_outer_core;
  gpu_realw_mem d_xiz_outer_core;
  gpu_realw_mem d_etax_outer_core;
  gpu_realw_mem d_etay_outer_core;
  gpu_realw_mem d_etaz_outer_core;
  gpu_realw_mem d_gammax_outer_core;
  gpu_realw_mem d_gammay_outer_core;
  gpu_realw_mem d_gammaz_outer_core;

  // model parameters
  gpu_realw_mem d_rhostore_outer_core;
  gpu_realw_mem d_kappavstore_outer_core;

  gpu_realw_mem d_rmass_outer_core;
  gpu_realw_mem d_b_rmass_outer_core;

  // global indexing
  gpu_int_mem d_ibool_outer_core;

  // mesh locations
  gpu_realw_mem d_rstore_outer_core;

  // wavefields
  // displacement, velocity, acceleration
  gpu_realw_mem d_displ_outer_core;
  gpu_realw_mem d_veloc_outer_core;
  gpu_realw_mem d_accel_outer_core;
  // backward/reconstructed elastic wavefield
  gpu_realw_mem d_b_displ_outer_core;
  gpu_realw_mem d_b_veloc_outer_core;
  gpu_realw_mem d_b_accel_outer_core;

  // kernels
  gpu_realw_mem d_rho_kl_outer_core;
  gpu_realw_mem d_alpha_kl_outer_core;

  // inner / outer elements
  gpu_int_mem d_phase_ispec_inner_outer_core;
  int num_phase_ispec_outer_core;

  int nspec_outer_outer_core;
  int nspec_inner_outer_core;
  int nspec2D_top_outer_core;
  int nspec2D_bottom_outer_core;

  int num_colors_inner_outer_core;
  int num_colors_outer_outer_core;
  int *h_num_elem_colors_outer_core;

  gpu_int_mem d_ibelm_top_outer_core;
  gpu_int_mem d_ibelm_bottom_outer_core;

  // normals definitions for coupling regions
  gpu_realw_mem d_normal_top_outer_core;
  gpu_realw_mem d_normal_bottom_outer_core;

  // Jacobian definitions
  gpu_realw_mem d_jacobian2D_top_outer_core;
  gpu_realw_mem d_jacobian2D_bottom_outer_core;

  // ------------------------------------------------------------------   //
  // inner_core
  // ------------------------------------------------------------------   //
  int NSPEC_INNER_CORE;
  int NGLOB_INNER_CORE;
  int NSPEC_INNER_CORE_STRAIN_ONLY;

  // interpolators
  gpu_realw_mem d_xix_inner_core;
  gpu_realw_mem d_xiy_inner_core;
  gpu_realw_mem d_xiz_inner_core;
  gpu_realw_mem d_etax_inner_core;
  gpu_realw_mem d_etay_inner_core;
  gpu_realw_mem d_etaz_inner_core;
  gpu_realw_mem d_gammax_inner_core;
  gpu_realw_mem d_gammay_inner_core;
  gpu_realw_mem d_gammaz_inner_core;

  // model parameters
  gpu_realw_mem d_rhostore_inner_core;
  gpu_realw_mem d_kappavstore_inner_core;
  gpu_realw_mem d_muvstore_inner_core;

  gpu_realw_mem d_rmassx_inner_core;
  gpu_realw_mem d_rmassy_inner_core;
  gpu_realw_mem d_rmassz_inner_core;
  gpu_realw_mem d_b_rmassx_inner_core;
  gpu_realw_mem d_b_rmassy_inner_core;
  gpu_realw_mem d_b_rmassz_inner_core;

  // global indexing
  gpu_int_mem d_ibool_inner_core;
  gpu_int_mem d_idoubling_inner_core;

  // mesh locations
  gpu_realw_mem d_rstore_inner_core;

  // anisotropic 3D mantle
  gpu_realw_mem d_c11store_inner_core;
  gpu_realw_mem d_c12store_inner_core;
  gpu_realw_mem d_c13store_inner_core;
  gpu_realw_mem d_c33store_inner_core;
  gpu_realw_mem d_c44store_inner_core;

  // wavefields
  // displacement, velocity, acceleration
  gpu_realw_mem d_displ_inner_core;
  gpu_realw_mem d_veloc_inner_core;
  gpu_realw_mem d_accel_inner_core;
  // backward/reconstructed elastic wavefield
  gpu_realw_mem d_b_displ_inner_core;
  gpu_realw_mem d_b_veloc_inner_core;
  gpu_realw_mem d_b_accel_inner_core;

  // attenuation
  gpu_realw_mem d_R_xx_inner_core;
  gpu_realw_mem d_R_yy_inner_core;
  gpu_realw_mem d_R_xy_inner_core;
  gpu_realw_mem d_R_xz_inner_core;
  gpu_realw_mem d_R_yz_inner_core;

  gpu_realw_mem d_b_R_xx_inner_core;
  gpu_realw_mem d_b_R_yy_inner_core;
  gpu_realw_mem d_b_R_xy_inner_core;
  gpu_realw_mem d_b_R_xz_inner_core;
  gpu_realw_mem d_b_R_yz_inner_core;

  gpu_realw_mem d_factor_common_inner_core;
  gpu_realw_mem d_one_minus_sum_beta_inner_core;

  gpu_realw_mem d_epsilondev_xx_inner_core;
  gpu_realw_mem d_epsilondev_yy_inner_core;
  gpu_realw_mem d_epsilondev_xy_inner_core;
  gpu_realw_mem d_epsilondev_xz_inner_core;
  gpu_realw_mem d_epsilondev_yz_inner_core;

  gpu_realw_mem d_b_epsilondev_xx_inner_core;
  gpu_realw_mem d_b_epsilondev_yy_inner_core;
  gpu_realw_mem d_b_epsilondev_xy_inner_core;
  gpu_realw_mem d_b_epsilondev_xz_inner_core;
  gpu_realw_mem d_b_epsilondev_yz_inner_core;

  gpu_realw_mem d_eps_trace_over_3_inner_core;
  gpu_realw_mem d_b_eps_trace_over_3_inner_core;

  // kernels
  gpu_realw_mem d_rho_kl_inner_core;
  gpu_realw_mem d_alpha_kl_inner_core;
  gpu_realw_mem d_beta_kl_inner_core;

  // inner / outer elements
  gpu_int_mem d_phase_ispec_inner_inner_core;
  int num_phase_ispec_inner_core;

  int nspec_outer_inner_core;
  int nspec_inner_inner_core;
  int nspec2D_top_inner_core;

  int num_colors_inner_inner_core;
  int num_colors_outer_inner_core;
  int *h_num_elem_colors_inner_core;

  gpu_int_mem d_ibelm_top_inner_core;

  // ------------------------------------------------------------------   //
  // oceans
  // ------------------------------------------------------------------   //
  int npoin_oceans;

  // model parameter
  gpu_int_mem d_ibool_ocean_load;
  gpu_realw_mem d_rmass_ocean_load;
  gpu_realw_mem d_normal_ocean_load;

  // ------------------------------------------------------------------   //
  // attenuation
  // ------------------------------------------------------------------   //
  gpu_realw_mem d_alphaval;
  gpu_realw_mem d_betaval;
  gpu_realw_mem d_gammaval;

  gpu_realw_mem d_b_alphaval;
  gpu_realw_mem d_b_betaval;
  gpu_realw_mem d_b_gammaval;

  // ------------------------------------------------------------------   //
  // GLL points & weights
  // ------------------------------------------------------------------   //

  // pointers to constant memory arrays
  gpu_realw_mem d_hprime_xx;
  //gpu_realw_mem d_hprime_yy;   // only needed if NGLLX != NGLLY != NGLLZ
  //gpu_realw_mem d_hprime_zz;   // only needed if NGLLX != NGLLY != NGLLZ

  gpu_realw_mem d_hprimewgll_xx;
  //gpu_realw_mem d_hprimewgll_yy;   // only needed if NGLLX != NGLLY != NGLLZ
  //gpu_realw_mem d_hprimewgll_zz;   // only needed if NGLLX != NGLLY != NGLLZ

  gpu_realw_mem d_wgllwgll_xy;
  gpu_realw_mem d_wgllwgll_xz;
  gpu_realw_mem d_wgllwgll_yz;
  gpu_realw_mem d_wgll_cube;

  // simulation type: 1 = forward, 2 = adjoint, 3 = kernel
  int simulation_type;

  // mesh coloring flag
  int use_mesh_coloring_gpu;

  // simulation flags
  int save_forward;
  int absorbing_conditions;
  int save_stacey;

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

  // ------------------------------------------------------------------   //
  // gravity
  // ------------------------------------------------------------------   //
  gpu_realw_mem d_d_ln_density_dr_table;   // needed also for no gravity case
  gpu_realw_mem d_minus_rho_g_over_kappa_fluid;
  gpu_realw_mem d_minus_gravity_table;
  gpu_realw_mem d_minus_deriv_gravity_table;
  gpu_realw_mem d_density_table;

  realw minus_g_icb;
  realw minus_g_cmb;

  realw RHO_BOTTOM_OC;
  realw RHO_TOP_OC;

  // ------------------------------------------------------------------   //
  // rotation
  // ------------------------------------------------------------------   //
  realw two_omega_earth;
  gpu_realw_mem d_A_array_rotation;
  gpu_realw_mem d_B_array_rotation;

  // needed for backward/reconstructed fields (kernel runs)
  realw b_two_omega_earth;
  gpu_realw_mem d_b_A_array_rotation;
  gpu_realw_mem d_b_B_array_rotation;

  // ------------------------------------------------------------------   //
  // sources
  // ------------------------------------------------------------------   //
  int nsources_local;
  gpu_realw_mem d_sourcearrays;
  gpu_double_mem d_stf_pre_compute;
  gpu_int_mem d_islice_selected_source;
  gpu_int_mem d_ispec_selected_source;

  // ------------------------------------------------------------------   //
  // receivers
  // ------------------------------------------------------------------   //
  gpu_int_mem d_number_receiver_global;
  gpu_int_mem d_ispec_selected_rec;
  gpu_int_mem d_islice_selected_rec;

  int nrec_local;

  gpu_realw_mem d_station_seismo_field;
  realw *h_station_seismo_field;

  gpu_realw_mem d_station_strain_field;
  realw* h_station_strain_field;

  // adjoint receivers/sources
  int nadj_rec_local;
  gpu_realw_mem d_adj_sourcearrays;
  realw *h_adj_sourcearrays_slice;
  gpu_int_mem d_pre_computed_irec;

  // norm checking
  gpu_realw_mem d_norm_max;

  // ------------------------------------------------------------------   //
  // assembly
  // ------------------------------------------------------------------   //
  int myrank;

  int num_interfaces_crust_mantle;
  int max_nibool_interfaces_cm;
  gpu_int_mem d_nibool_interfaces_crust_mantle;
  gpu_int_mem d_ibool_interfaces_crust_mantle;

  gpu_realw_mem d_send_accel_buffer_crust_mantle;
  gpu_realw_mem d_b_send_accel_buffer_crust_mantle;

  int num_interfaces_inner_core;
  int max_nibool_interfaces_ic;
  gpu_int_mem d_nibool_interfaces_inner_core;
  gpu_int_mem d_ibool_interfaces_inner_core;

  gpu_realw_mem d_send_accel_buffer_inner_core;
  gpu_realw_mem d_b_send_accel_buffer_inner_core;

  int num_interfaces_outer_core;
  int max_nibool_interfaces_oc;
  gpu_int_mem d_nibool_interfaces_outer_core;
  gpu_int_mem d_ibool_interfaces_outer_core;

  gpu_realw_mem d_send_accel_buffer_outer_core;
  gpu_realw_mem d_b_send_accel_buffer_outer_core;

  // ------------------------------------------------------------------   //
  // absorbing boundaries
  // ------------------------------------------------------------------   //

  int nspec2D_xmin_crust_mantle, nspec2D_xmax_crust_mantle;
  int nspec2D_ymin_crust_mantle, nspec2D_ymax_crust_mantle;

  gpu_int_mem d_nimin_crust_mantle;
  gpu_int_mem d_nimax_crust_mantle;
  gpu_int_mem d_njmin_crust_mantle;
  gpu_int_mem d_njmax_crust_mantle;
  gpu_int_mem d_nkmin_xi_crust_mantle;
  gpu_int_mem d_nkmin_eta_crust_mantle;

  gpu_int_mem d_ibelm_xmin_crust_mantle;
  gpu_int_mem d_ibelm_xmax_crust_mantle;
  gpu_int_mem d_ibelm_ymin_crust_mantle;
  gpu_int_mem d_ibelm_ymax_crust_mantle;

  gpu_realw_mem d_normal_xmin_crust_mantle;
  gpu_realw_mem d_normal_xmax_crust_mantle;
  gpu_realw_mem d_normal_ymin_crust_mantle;
  gpu_realw_mem d_normal_ymax_crust_mantle;

  gpu_realw_mem d_jacobian2D_xmin_crust_mantle;
  gpu_realw_mem d_jacobian2D_xmax_crust_mantle;
  gpu_realw_mem d_jacobian2D_ymin_crust_mantle;
  gpu_realw_mem d_jacobian2D_ymax_crust_mantle;

  gpu_realw_mem d_absorb_xmin_crust_mantle;
  gpu_realw_mem d_absorb_xmax_crust_mantle;
  gpu_realw_mem d_absorb_ymin_crust_mantle;
  gpu_realw_mem d_absorb_ymax_crust_mantle;

  gpu_realw_mem d_rho_vp_crust_mantle;
  gpu_realw_mem d_rho_vs_crust_mantle;

  int nspec2D_xmin_outer_core, nspec2D_xmax_outer_core;
  int nspec2D_ymin_outer_core, nspec2D_ymax_outer_core;
  int nspec2D_zmin_outer_core;

  gpu_int_mem d_nimin_outer_core;
  gpu_int_mem d_nimax_outer_core;
  gpu_int_mem d_njmin_outer_core;
  gpu_int_mem d_njmax_outer_core;
  gpu_int_mem d_nkmin_xi_outer_core;
  gpu_int_mem d_nkmin_eta_outer_core;

  gpu_int_mem d_ibelm_xmin_outer_core;
  gpu_int_mem d_ibelm_xmax_outer_core;
  gpu_int_mem d_ibelm_ymin_outer_core;
  gpu_int_mem d_ibelm_ymax_outer_core;

  gpu_realw_mem d_jacobian2D_xmin_outer_core;
  gpu_realw_mem d_jacobian2D_xmax_outer_core;
  gpu_realw_mem d_jacobian2D_ymin_outer_core;
  gpu_realw_mem d_jacobian2D_ymax_outer_core;

  gpu_realw_mem d_absorb_xmin_outer_core;
  gpu_realw_mem d_absorb_xmax_outer_core;
  gpu_realw_mem d_absorb_ymin_outer_core;
  gpu_realw_mem d_absorb_ymax_outer_core;
  gpu_realw_mem d_absorb_zmin_outer_core;

  gpu_realw_mem d_vp_outer_core;

  // ------------------------------------------------------------------   //
  // noise tomography
  // ------------------------------------------------------------------   //
  int noise_tomography;

  gpu_realw_mem d_noise_surface_movie;
  gpu_realw_mem d_noise_sourcearray;

  gpu_realw_mem d_normal_x_noise;
  gpu_realw_mem d_normal_y_noise;
  gpu_realw_mem d_normal_z_noise;
  gpu_realw_mem d_mask_noise;

  // free surface
  int nspec2D_top_crust_mantle;
  gpu_int_mem d_ibelm_top_crust_mantle;
  gpu_realw_mem d_jacobian2D_top_crust_mantle;

  // noise sensitivity kernel
  gpu_realw_mem d_Sigma_kl;

  // ------------------------------------------------------------------ //
  // optimizations
  // ------------------------------------------------------------------ //

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

#ifdef USE_OPENCL
  // pinned memory allocated by ALLOC_PINNED_BUFFER_OCL
  cl_mem h_pinned_station_seismo_field;
  cl_mem h_pinned_adj_sourcearrays_slice;

  // crust mantle
  cl_mem h_pinned_send_accel_buffer_cm;
  cl_mem h_pinned_recv_accel_buffer_cm;
  cl_mem h_pinned_b_send_accel_buffer_cm;
  cl_mem h_pinned_b_recv_accel_buffer_cm;

  // inner core
  cl_mem h_pinned_send_accel_buffer_ic;
  cl_mem h_pinned_recv_accel_buffer_ic;
  cl_mem h_pinned_b_send_accel_buffer_ic;
  cl_mem h_pinned_b_recv_accel_buffer_ic;

  // outer core
  cl_mem h_pinned_send_accel_buffer_oc;
  cl_mem h_pinned_recv_accel_buffer_oc;
  cl_mem h_pinned_b_send_accel_buffer_oc;
  cl_mem h_pinned_b_recv_accel_buffer_oc;
#endif

  // streams
#ifdef USE_CUDA
  // overlapped memcpy streams
  cudaStream_t compute_stream;
  cudaStream_t copy_stream;
#endif

#ifdef USE_OPENCL
  cl_event last_copy_evt;
  int has_last_copy_evt;
#endif

  // specific OpenCL texture arrays
#ifdef USE_OPENCL
// note: need to be defined as they are passed as function arguments
#ifdef USE_TEXTURES_FIELDS
  // forward
  cl_mem d_displ_cm_tex;
  cl_mem d_accel_cm_tex;

  cl_mem d_displ_oc_tex;
  cl_mem d_accel_oc_tex;

  cl_mem d_displ_ic_tex;
  cl_mem d_accel_ic_tex;

  // backward/reconstructed
  cl_mem d_b_displ_cm_tex;
  cl_mem d_b_accel_cm_tex;

  cl_mem d_b_displ_oc_tex;
  cl_mem d_b_accel_oc_tex;

  cl_mem d_b_displ_ic_tex;
  cl_mem d_b_accel_ic_tex;
#endif
#ifdef USE_TEXTURES_CONSTANTS
  // hprime
  cl_mem d_hprime_xx_cm_tex;
  // weighted hprime
  cl_mem d_hprimewgll_xx_cm_tex;
#endif
#endif

  // ------------------------------------------------------------------ //
  // LDDRK
  // ------------------------------------------------------------------ //
  int use_lddrk;
  // daniel debug: todo - add lddrk arrays here and in gpu_buffer_list.c for initialization

} Mesh;


/*----------------------------------------------------------------------------------------------- */
// utility functions
/*----------------------------------------------------------------------------------------------- */

// defined in helper_functions_gpu.c
void gpuCreateCopy_todevice_int (gpu_int_mem *d_array_addr_ptr, int *h_array, size_t size);
void gpuCreateCopy_todevice_realw (gpu_realw_mem *d_array_addr_ptr, realw *h_array, size_t size);

void gpuCopy_todevice_realw (gpu_realw_mem *d_array_addr_ptr, realw *h_array, size_t size);
void gpuCopy_todevice_double (gpu_double_mem *d_array_addr_ptr, double *h_array, size_t size);
void gpuCopy_todevice_int (gpu_int_mem *d_array_addr_ptr, int *h_array, size_t size);

void gpuCopy_from_device_realw (gpu_realw_mem *d_array_addr_ptr, realw *h_array, size_t size);

void gpuCopy_todevice_realw_offset (gpu_realw_mem *d_array_addr_ptr, realw *h_array, size_t size, size_t offset);
void gpuCopy_from_device_realw_offset (gpu_realw_mem *d_array_addr_ptr, realw *h_array, size_t size, size_t offset);

void gpuRegisterHost_realw ( realw *h_array, const size_t size);
void gpuUnregisterHost_realw ( realw *h_array);

void gpuMalloc_int (gpu_int_mem *buffer, size_t size);
void gpuMalloc_realw (gpu_realw_mem *buffer, size_t size);
void gpuMalloc_double (gpu_double_mem *buffer, size_t size);

void gpuMemset_realw (gpu_realw_mem *buffer, size_t size, int value);

void gpuSetConst (gpu_realw_mem *buffer, size_t size, realw *array);
void gpuFree (void *d_array_addr_ptr);
void gpuInitialize_buffers (Mesh *mp);
void gpuSynchronize ();
void gpuReset ();

void exit_on_gpu_error (const char *kernel_name);
void exit_on_error (const char *info);

void synchronize_mpi ();
double get_time_val ();

// defined in check_fields_gpu.c
void get_free_memory (double *free_db, double *used_db, double *total_db);
realw get_device_array_maximum_value (gpu_realw_mem d_array, int size);

/* ----------------------------------------------------------------------------------------------- */

// OpenCL / CUDA macro definitions

#ifndef TAKE_REF_OCL
#define TAKE_REF_OCL(_buffer_)
#endif
#ifndef TAKE_REF_CUDA
#define TAKE_REF_CUDA(_buffer_)
#endif

#define gpuTakeRef(_buffer_) _buffer_ ;         \
  TAKE_REF_OCL(_buffer_);                       \
  TAKE_REF_CUDA(_buffer_);

#ifndef INITIALIZE_OFFSET_OCL
#define INITIALIZE_OFFSET_OCL()
#endif
#ifndef INITIALIZE_OFFSET_CUDA
#define INITIALIZE_OFFSET_CUDA()
#endif

#define INITIALIZE_OFFSET()                     \
  INITIALIZE_OFFSET_OCL();                      \
  INITIALIZE_OFFSET_CUDA();

#ifndef INIT_OFFSET_OCL
#define INIT_OFFSET_OCL(_buffer_, _offset_)
#endif
#ifndef INIT_OFFSET_CUDA
#define INIT_OFFSET_CUDA(_buffer_, _offset_)
#endif

#define INIT_OFFSET(_buffer_, _offset_)         \
  __typeof__(mp->_buffer_) _buffer_##_##_offset_;   \
  INIT_OFFSET_OCL(_buffer_, _offset_);           \
  INIT_OFFSET_CUDA(_buffer_, _offset_);

#define PASS_OFFSET(_buffer_, _offset_) _buffer_ ##_##  _offset_

#ifndef RELEASE_OFFSET_OCL
#define RELEASE_OFFSET_OCL(_buffer_, _offset_) {}
#endif
#ifndef RELEASE_OFFSET_CUDA
#define RELEASE_OFFSET_CUDA(_buffer_, _offset_) {}
#endif

#define RELEASE_OFFSET(_buffer_, _offset_)      \
  RELEASE_OFFSET_OCL(_buffer_, _offset_);        \
  RELEASE_OFFSET_CUDA(_buffer_, _offset_);

/* ----------------------------------------------------------------------------------------------- */
// kernel setup function
/* ----------------------------------------------------------------------------------------------- */
// moved here into header to inline function calls if possible

static inline void get_blocks_xy (int num_blocks, int *num_blocks_x, int *num_blocks_y) {
  // Initially sets the blocks_x to be the num_blocks, and adds rows as needed (block size limit of 65535).
  // If an additional row is added, the row length is cut in
  // half. If the block count is odd, there will be 1 too many blocks,
  // which must be managed at runtime with an if statement.

  *num_blocks_x = num_blocks;
  *num_blocks_y = 1;

  while (*num_blocks_x > MAXIMUM_GRID_DIM) {
    *num_blocks_x = (int) ceil (*num_blocks_x * 0.5f);
    *num_blocks_y = *num_blocks_y * 2;
  }

#if DEBUG == 1
  printf("work group - total %d has group size x = %d / y = %d\n",
         num_blocks,*num_blocks_x,*num_blocks_y);
#endif

  // tries to balance x- and y-group
#ifdef BALANCE_WORK_GROUP
  if (*num_blocks_x > BALANCE_WORK_GROUP_UNITS && *num_blocks_y < BALANCE_WORK_GROUP_UNITS){
    while (*num_blocks_x > BALANCE_WORK_GROUP_UNITS && *num_blocks_y < BALANCE_WORK_GROUP_UNITS) {
      *num_blocks_x = (int) ceil (*num_blocks_x * 0.5f);
      *num_blocks_y = *num_blocks_y * 2;
    }
  }

#if DEBUG == 1
  printf("balancing work group with limit size %d - total %d has group size x = %d / y = %d\n",
         BALANCE_WORK_GROUP_UNITS,num_blocks,*num_blocks_x,*num_blocks_y);
#endif

#endif
}


#endif   // MESH_CONSTANTS_GPU_H
