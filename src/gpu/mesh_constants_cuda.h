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

#ifndef MESH_CONSTANTS_CUDA_H
#define MESH_CONSTANTS_CUDA_H

#include "kernel_proto.cu.h"


void print_CUDA_error_if_any(cudaError_t err, int num);

/* ----------------------------------------------------------------------------------------------- */

// textures
typedef texture<float, cudaTextureType1D, cudaReadModeElementType> realw_texture;

// restricted pointers: improves performance on Kepler ~ 10%
// see: http://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#restrict
typedef const float* __restrict__ realw_const_p; // otherwise use: //typedef const float* realw_const_p;
typedef float* __restrict__ realw_p; // otherwise use: //typedef float* realw_p;

/* ----------------------------------------------------------------------------------------------- */

#define INITIALIZE_OFFSET_CUDA()

#define INIT_OFFSET_CUDA(_buffer_, _offset_)                        \
do {                                                                \
  if (run_opencl) {                                                 \
    _buffer_##_##_offset_.cuda = mp->_buffer_.cuda + _offset_;      \
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

#endif

