/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  3 . 0
 !               ---------------------------------------
 !
 !     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
 !                        Princeton University, USA
 !                and CNRS / University of Marseille, France
 !                 (there are currently many more authors!)
 ! (c) Princeton University and CNRS / University of Marseille, July 2012
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

#ifndef MESH_CONSTANTS_HIP_H
#define MESH_CONSTANTS_HIP_H

// HIP specifics

#ifdef USE_HIP

// (optional) unrolling loops
// leads up to ~3% performance increase with HIP
#define MANUALLY_UNROLLED_LOOPS   // uncomment to use loops

// AMD MI100
#ifdef GPU_DEVICE_MI100
#undef USE_LAUNCH_BOUNDS
#endif

// AMD MI250X
#ifdef GPU_DEVICE_MI250
#undef USE_LAUNCH_BOUNDS
//#define USE_LAUNCH_BOUNDS     // will slow down kernels...
//#define LAUNCH_MIN_BLOCKS 7
#endif

/*----------------------------------------------------------------------------------------------- */

// definitions
typedef hipEvent_t gpu_event;
typedef hipStream_t gpu_stream;

// hip header files
#include "kernel_proto.cu.h"

static inline void print_HIP_error_if_any(hipError_t err, int num) {
  if (hipSuccess != err)
  {
    printf("\nHIP error !!!!! <%s> !!!!! \nat HIP call error code: # %d\n",hipGetErrorString(err),num);
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
      fprintf(fp,"\nHIP error !!!!! <%s> !!!!! \nat HIP call error code: # %d\n",hipGetErrorString(err),num);
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

#define INITIALIZE_OFFSET_HIP()

#define INIT_OFFSET_HIP(_buffer_, _offset_)                        \
do {                                                                \
  if (run_hip) {                                                   \
    if (mp->_buffer_.hip != NULL){                                 \
      _buffer_##_##_offset_.hip = mp->_buffer_.hip + _offset_;    \
    } else {                                                        \
      _buffer_##_##_offset_.hip = NULL;                            \
    }                                                               \
  }                                                                 \
} while (0)

#define RELEASE_OFFSET_HIP(_buffer_, _offset_)

#define TAKE_REF_HIP(_buffer_)

/* ----------------------------------------------------------------------------------------------- */

#endif  // USE_HIP

#endif  // MESH_CONSTANTS_HIP_H
