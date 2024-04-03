/*
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
*/

#ifndef MESH_CONSTANTS_CUDA_H
#define MESH_CONSTANTS_CUDA_H

// CUDA specifics

#ifdef USE_CUDA

// (optional) unrolling loops
// leads up to ~10% performance increase in OpenCL and ~1% in Cuda
#define MANUALLY_UNROLLED_LOOPS   // uncomment to use loops

/*----------------------------------------------------------------------------------------------- */

// definitions
typedef cudaEvent_t gpu_event;
typedef cudaStream_t gpu_stream;

// cuda header files
#include "kernel_proto.cu.h"

static inline void print_CUDA_error_if_any(cudaError_t err, int num) {
  if (cudaSuccess != err)
  {
    printf("\nCUDA error !!!!! <%s> !!!!! \nat CUDA call error code: # %d\n",cudaGetErrorString(err),num);
    fflush(stdout);

    // outputs error file
    FILE* fp;
    int myrank;
    char filename[BUFSIZ];
#ifdef WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#else
    myrank = 0;
#endif
    sprintf(filename,"OUTPUT_FILES/error_message_%06d.txt",myrank);
    fp = fopen(filename,"a+");
    if (fp != NULL){
      fprintf(fp,"\nCUDA error !!!!! <%s> !!!!! \nat CUDA call error code: # %d\n",cudaGetErrorString(err),num);
      fclose(fp);
    }

    // stops program
#ifdef WITH_MPI
    MPI_Abort(MPI_COMM_WORLD,1);
#endif
    exit(EXIT_FAILURE);
  }
}


/* ----------------------------------------------------------------------------------------------- */

#ifndef CUSTOM_REAL
#pragma message ("\nmesh_constants_cuda.h: CUSTOM_REAL not defined for textures, using CUSTOM_REAL == 4\n")
#define CUSTOM_REAL 4
#endif

#if CUSTOM_REAL == 4
// textures
// textures
// note: texture templates are supported only for CUDA versions <= 11.x
//       since CUDA 12.x, these are deprecated and texture objects should be used instead
//       see: https://developer.nvidia.com/blog/cuda-pro-tip-kepler-texture-objects-improve-performance-and-flexibility/
#if defined(USE_TEXTURES_FIELDS) || defined(USE_TEXTURES_CONSTANTS)
typedef texture<float, cudaTextureType1D, cudaReadModeElementType> realw_texture;
#endif
// restricted pointers: improves performance on Kepler ~ 10%
// see: http://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#restrict
typedef const float* __restrict__ realw_const_p; // otherwise use: //typedef const float* realw_const_p;
typedef float* __restrict__ realw_p; // otherwise use: //typedef float* realw_p;

#elif CUSTOM_REAL == 8
// textures
// textures
// note: texture templates are supported only for CUDA versions <= 11.x
//       since CUDA 12.x, these are deprecated and texture objects should be used instead
//       see: https://developer.nvidia.com/blog/cuda-pro-tip-kepler-texture-objects-improve-performance-and-flexibility/
#if defined(USE_TEXTURES_FIELDS) || defined(USE_TEXTURES_CONSTANTS)
typedef texture<double, cudaTextureType1D, cudaReadModeElementType> realw_texture;
#endif
// restricted pointers: improves performance on Kepler ~ 10%
// see: http://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#restrict
typedef const double* __restrict__ realw_const_p; // otherwise use: //typedef const double* realw_const_p;
typedef double* __restrict__ realw_p; // otherwise use: //typedef double* realw_p;
#endif

/* ----------------------------------------------------------------------------------------------- */

#define INITIALIZE_OFFSET_CUDA()

#define INIT_OFFSET_CUDA(_buffer_, _offset_)                        \
do {                                                                \
  if (run_cuda) {                                                   \
    if (mp->_buffer_.cuda != NULL){                                 \
      _buffer_##_##_offset_.cuda = mp->_buffer_.cuda + _offset_;    \
    } else {                                                        \
      _buffer_##_##_offset_.cuda = NULL;                            \
    }                                                               \
  }                                                                 \
} while (0)

#define RELEASE_OFFSET_CUDA(_buffer_, _offset_)

#define TAKE_REF_CUDA(_buffer_)

/* ----------------------------------------------------------------------------------------------- */

#ifndef USE_OLDER_CUDA4_GPU
  #ifdef USE_TEXTURES_FIELDS
    // forward
    extern realw_texture d_displ_cm_tex;
    extern realw_texture d_accel_cm_tex;

    extern realw_texture d_displ_oc_tex;
    extern realw_texture d_accel_oc_tex;

    extern realw_texture d_displ_ic_tex;
    extern realw_texture d_accel_ic_tex;

    // backward/reconstructed
    extern realw_texture d_b_displ_cm_tex;
    extern realw_texture d_b_accel_cm_tex;

    extern realw_texture d_b_displ_oc_tex;
    extern realw_texture d_b_accel_oc_tex;

    extern realw_texture d_b_displ_ic_tex;
    extern realw_texture d_b_accel_ic_tex;
  #endif

  #ifdef USE_TEXTURES_CONSTANTS
    // hprime
    extern realw_texture d_hprime_xx_tex;
    extern __constant__ size_t d_hprime_xx_tex_offset;
    // weighted hprime
    extern realw_texture d_hprimewgll_xx_tex;
    extern __constant__ size_t d_hprimewgll_xx_tex_offset;
  #endif
#endif

#endif  // USE_CUDA

#endif

