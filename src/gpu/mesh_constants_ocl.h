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

#ifndef MESH_CONSTANTS_OCL_H
#define MESH_CONSTANTS_OCL_H

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

const char* clewErrorString (cl_int error);

#define INITIALIZE_OFFSET_OCL()                     \
  cl_uint buffer_create_type;                       \
  size_t size;                                      \
  cl_buffer_region region_type;

#define INIT_OFFSET_OCL(_buffer_, _offset_)                             \
do {                                                                    \
  if (run_opencl) {                                                     \
    clCheck (clGetMemObjectInfo (mp->_buffer_.ocl, CL_MEM_FLAGS, sizeof(cl_uint), &buffer_create_type, NULL)); \
    clCheck (clGetMemObjectInfo (mp->_buffer_.ocl, CL_MEM_SIZE , sizeof(size_t), &size, NULL)); \
                                                                        \
    region_type.origin = _offset_ * sizeof(CL_FLOAT);                   \
    region_type.size = size;                                            \
                                                                        \
    _buffer_##_##_offset_.ocl = clCreateSubBuffer (mp->_buffer_.ocl, buffer_create_type, CL_BUFFER_CREATE_TYPE_REGION, \
                                                   (void *) &region_type, clck_(&mocl_errcode)); \
  }                                                                     \
} while (0)

#define RELEASE_OFFSET_OCL(_buffer_, _offset_)                  \
do {                                                            \
  if (run_opencl) {                                             \
    clCheck(clReleaseMemObject(_buffer_##_##_offset_.ocl));     \
  }                                                             \
} while (0)

#define ALLOC_PINNED_BUFFER_OCL(_buffer_, _size_)                               \
do {                                                                            \
  mp->h_pinned_##_buffer_ = clCreateBuffer(mocl.context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, \
                                           _size_, NULL, clck_(&mocl_errcode)); \
  mp->h_##_buffer_ = (realw *) clEnqueueMapBuffer(mocl.command_queue, mp->h_pinned_##_buffer_, CL_TRUE, \
                                                  CL_MAP_READ | CL_MAP_WRITE, 0, _size_, 0, \
                                                  NULL, NULL, clck_(&mocl_errcode)); \
} while (0)

#define RELEASE_PINNED_BUFFER_OCL(_buffer_)                                       \
do {                                                                              \
    clCheck(clEnqueueUnmapMemObject(mocl.command_queue, mp->h_pinned_##_buffer_,  \
                                    mp->h_##_buffer_, 0, NULL, NULL));            \
    clCheck(clReleaseMemObject (mp->h_pinned_##_buffer_));                        \
} while (0)


/* ----------------------------------------------------------------------------------------------- */

extern int mocl_errcode;

static inline cl_int _clCheck(cl_int errcode, const char *file, int line, const char *func) {
  mocl_errcode = errcode;
  if (mocl_errcode != CL_SUCCESS) {
    fprintf (stderr, "OpenCL Error %d/%s at %s:%d %s\n", mocl_errcode,
             clewErrorString(mocl_errcode),
             file, line, func);
    fflush(NULL);
    exit(1);
  }
  return errcode;
}

#define clCheck(to_check) _clCheck(to_check,__FILE__, __LINE__,  __func__)

#define clck_(var) var); clCheck(*var

/* ----------------------------------------------------------------------------------------------- */

#define TAKE_REF_OCL(_buffer_)                                  \
do {                                                            \
  if (run_opencl) {                                             \
    clCheck(clRetainMemObject(_buffer_.ocl));                   \
  }                                                             \
} while (0)


/* ----------------------------------------------------------------------------------------------- */

struct mesh_programs_s {
#undef BOAST_KERNEL
#define BOAST_KERNEL(__kern_name__) cl_program __kern_name__##_program

  #include "kernel_list.h"
};

/* ----------------------------------------------------------------------------------------------- */

struct mesh_kernels_s {
#undef BOAST_KERNEL
#define BOAST_KERNEL(__kern_name__) cl_kernel __kern_name__

  #include "kernel_list.h"
};


/* ----------------------------------------------------------------------------------------------- */

extern struct _mesh_opencl {
  struct mesh_programs_s programs;
  struct mesh_kernels_s kernels;
  cl_command_queue command_queue;
  cl_command_queue copy_queue;
  cl_context context;
  cl_device_id device;
  cl_int nb_devices;
} mocl;

#endif
