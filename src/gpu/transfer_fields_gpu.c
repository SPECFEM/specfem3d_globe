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

#include "mesh_constants_gpu.h"

// crust_mantle
extern EXTERN_LANG
void FC_FUNC_(transfer_fields_cm_to_device,
              TRANSFER_FIELDS_CM_TO_DEVICE)(int *size, realw *displ, realw *veloc, realw *accel, long *Mesh_pointer_f) {

  TRACE("transfer_fields_cm_to_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_displ_crust_mantle.ocl, CL_FALSE, 0,
                                  sizeof (realw) * (*size),
                                  displ, 0, NULL, NULL));
    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_veloc_crust_mantle.ocl, CL_FALSE, 0,
                                  sizeof (realw) * (*size),
                                  veloc, 0, NULL, NULL));
    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_accel_crust_mantle.ocl, CL_FALSE, 0,
                                  sizeof (realw) * (*size),
                                  accel, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(mp->d_displ_crust_mantle.cuda,displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40003);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_veloc_crust_mantle.cuda,veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40004);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_accel_crust_mantle.cuda,accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40005);
  }
#endif
}

// inner_core
extern EXTERN_LANG
void FC_FUNC_(transfer_fields_ic_to_device,
              TRANSFER_FIELDS_IC_TO_DEVICE)(int *size, realw *displ, realw *veloc, realw *accel, long *Mesh_pointer_f) {

  TRACE("transfer_fields_ic_to_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_displ_inner_core.ocl, CL_FALSE, 0,
                                  sizeof (realw) * (*size),
                                  displ, 0, NULL, NULL));

    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_veloc_inner_core.ocl, CL_FALSE, 0,
                                  sizeof (realw) * (*size),
                                  veloc, 0, NULL, NULL));

    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_accel_inner_core.ocl, CL_FALSE, 0,
                                  sizeof (realw) * (*size),
                                  accel, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(mp->d_displ_inner_core.cuda,displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40003);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_veloc_inner_core.cuda,veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40004);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_accel_inner_core.cuda,accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40005);
  }
#endif
}

// outer_core
extern EXTERN_LANG
void FC_FUNC_(transfer_fields_oc_to_device,
              TRANSFER_FIELDS_OC_TO_DEVICE)(int *size, realw *displ, realw *veloc, realw *accel, long *Mesh_pointer_f) {

  TRACE("transfer_fields_oc_to_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_displ_outer_core.ocl, CL_FALSE, 0,
                                  sizeof (realw) * (*size),
                                  displ, 0, NULL, NULL));

    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_veloc_outer_core.ocl, CL_FALSE, 0,
                                  sizeof (realw) * (*size),
                                  veloc, 0, NULL, NULL));

    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_accel_outer_core.ocl, CL_FALSE, 0,
                                  sizeof (realw) * (*size),
                                  accel, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(mp->d_displ_outer_core.cuda,displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40003);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_veloc_outer_core.cuda,veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40004);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_accel_outer_core.cuda,accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40005);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// backward/reconstructed fields

// crust_mantle
extern EXTERN_LANG
void FC_FUNC_(transfer_b_fields_cm_to_device,
              TRANSFER_FIELDS_B_CM_TO_DEVICE)(int *size, realw *b_displ, realw *b_veloc, realw *b_accel,
                                              long *Mesh_pointer_f) {

  TRACE("transfer_fields_b_cm_to_device");
  // debug
  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_displ_crust_mantle.ocl, CL_FALSE, 0,
                                  sizeof (realw) * (*size),
                                  b_displ, 0, NULL, NULL));

    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_veloc_crust_mantle.ocl, CL_FALSE, 0,
                                  sizeof (realw) * (*size),
                                  b_veloc, 0, NULL, NULL));

    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_accel_crust_mantle.ocl, CL_FALSE, 0,
                                  sizeof (realw) * (*size),
                                  b_accel, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_displ_crust_mantle.cuda,b_displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40003);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_veloc_crust_mantle.cuda,b_veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40004);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_accel_crust_mantle.cuda,b_accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40005);
  }
#endif
}

// inner_core
extern EXTERN_LANG
void FC_FUNC_(transfer_b_fields_ic_to_device,
              TRANSFER_FIELDS_B_IC_TO_DEVICE)(int *size, realw *b_displ, realw *b_veloc, realw *b_accel,
                                              long *Mesh_pointer_f) {

  TRACE("transfer_fields_b_ic_to_device");
  // debug
  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_displ_inner_core.ocl, CL_FALSE, 0,
                                  sizeof (realw) * (*size),
                                  b_displ, 0, NULL, NULL));

    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_veloc_inner_core.ocl, CL_FALSE, 0,
                                  sizeof (realw) * (*size),
                                  b_veloc, 0, NULL, NULL));

    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_accel_inner_core.ocl, CL_FALSE, 0,
                                  sizeof (realw) * (*size),
                                  b_accel, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_displ_inner_core.cuda,b_displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40003);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_veloc_inner_core.cuda,b_veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40004);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_accel_inner_core.cuda,b_accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40005);
  }
#endif
}

// outer_core
extern EXTERN_LANG
void FC_FUNC_(transfer_b_fields_oc_to_device,
              TRANSFER_FIELDS_B_OC_TO_DEVICE)(int *size, realw *b_displ, realw *b_veloc, realw *b_accel,
                                              long *Mesh_pointer_f) {

  TRACE("transfer_fields_b_oc_to_device");

  // debug
  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_displ_outer_core.ocl, CL_FALSE, 0,
                                  sizeof (realw) * (*size),
                                  b_displ, 0, NULL, NULL));

    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_veloc_outer_core.ocl, CL_FALSE, 0,
                                  sizeof (realw) * (*size),
                                  b_veloc, 0, NULL, NULL));

    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_accel_outer_core.ocl, CL_TRUE, 0,
                                  sizeof (realw) * (*size),
                                  b_accel, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_displ_outer_core.cuda,b_displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40003);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_veloc_outer_core.cuda,b_veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40004);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_accel_outer_core.cuda,b_accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40005);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// transfer memory from GPU device to CPU host
/* ----------------------------------------------------------------------------------------------- */

// crust_mantle
extern EXTERN_LANG
void FC_FUNC_(transfer_fields_cm_from_device,
              TRANSFER_FIELDS_CM_FROM_DEVICE)(int *size, realw *displ, realw *veloc, realw *accel, long *Mesh_pointer_f) {

  TRACE("transfer_fields_cm_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_displ_crust_mantle.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 displ, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_veloc_crust_mantle.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 veloc, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_accel_crust_mantle.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 accel, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_displ_crust_mantle.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
    print_CUDA_error_if_any(cudaMemcpy(veloc,mp->d_veloc_crust_mantle.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
    print_CUDA_error_if_any(cudaMemcpy(accel,mp->d_accel_crust_mantle.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40008);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// inner_core
extern EXTERN_LANG
void FC_FUNC_(transfer_fields_ic_from_device,
              TRANSFER_FIELDS_IC_FROM_DEVICE)(int *size, realw *displ, realw *veloc, realw *accel, long *Mesh_pointer_f) {

  TRACE("transfer_fields_ic_from_device");
  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_displ_inner_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 displ, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_veloc_inner_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 veloc, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_accel_inner_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 accel, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_displ_inner_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
    print_CUDA_error_if_any(cudaMemcpy(veloc,mp->d_veloc_inner_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
    print_CUDA_error_if_any(cudaMemcpy(accel,mp->d_accel_inner_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40008);
  }
#endif
}
/* ----------------------------------------------------------------------------------------------- */


// outer_core
extern EXTERN_LANG
void FC_FUNC_(transfer_fields_oc_from_device,
              TRANSFER_FIELDS_OC_FROM_DEVICE)(int *size, realw *displ, realw *veloc, realw *accel, long *Mesh_pointer_f) {

  TRACE("transfer_fields_oc_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_displ_outer_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 displ, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_veloc_outer_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 veloc, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_accel_outer_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 accel, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_displ_outer_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
    print_CUDA_error_if_any(cudaMemcpy(veloc,mp->d_veloc_outer_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
    print_CUDA_error_if_any(cudaMemcpy(accel,mp->d_accel_outer_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40008);
  }
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// backward/reconstructed fields

// crust_mantle
extern EXTERN_LANG
void FC_FUNC_(transfer_b_fields_cm_from_device,
              TRANSFER_B_FIELDS_CM_FROM_DEVICE)(int *size, realw *b_displ, realw *b_veloc, realw *b_accel,
                                                long *Mesh_pointer_f) {

  TRACE("transfer_b_fields_cm_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_b_displ_crust_mantle.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 b_displ, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_b_veloc_crust_mantle.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 b_veloc, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_b_accel_crust_mantle.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 b_accel, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(b_displ,mp->d_b_displ_crust_mantle.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
    print_CUDA_error_if_any(cudaMemcpy(b_veloc,mp->d_b_veloc_crust_mantle.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
    print_CUDA_error_if_any(cudaMemcpy(b_accel,mp->d_b_accel_crust_mantle.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40008);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */


// inner_core
extern EXTERN_LANG
void FC_FUNC_(transfer_b_fields_ic_from_device,
              TRANSFER_B_FIELDS_IC_FROM_DEVICE)(int *size, realw *b_displ, realw *b_veloc, realw *b_accel,
                                                long *Mesh_pointer_f) {
  TRACE("transfer_fields_b_ic_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_b_displ_inner_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 b_displ, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_b_veloc_inner_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 b_veloc, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_b_accel_inner_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 b_accel, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(b_displ,mp->d_b_displ_inner_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
    print_CUDA_error_if_any(cudaMemcpy(b_veloc,mp->d_b_veloc_inner_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
    print_CUDA_error_if_any(cudaMemcpy(b_accel,mp->d_b_accel_inner_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40008);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */


// outer_core
extern EXTERN_LANG
void FC_FUNC_(transfer_b_fields_oc_from_device,
              TRANSFER_B_FIELDS_OC_FROM_DEVICE)(int *size, realw *b_displ, realw *b_veloc, realw *b_accel,
                                                long *Mesh_pointer_f) {

  TRACE("transfer_b_fields_oc_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_b_displ_outer_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 b_displ, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_b_veloc_outer_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 b_veloc, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_b_accel_outer_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 b_accel, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(b_displ,mp->d_b_displ_outer_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
    print_CUDA_error_if_any(cudaMemcpy(b_veloc,mp->d_b_veloc_outer_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
    print_CUDA_error_if_any(cudaMemcpy(b_accel,mp->d_b_accel_outer_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40008);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// single wavefield transfers
/* ----------------------------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------------------------- */
// displacements
/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_displ_cm_from_device,
              TRANSFER_DISPL_CM_FROM_DEVICE)(int *size, realw *displ, long *Mesh_pointer_f) {

  TRACE("transfer_displ_cm_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_displ_crust_mantle.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 displ, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_displ_crust_mantle.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_displ_cm_from_device,
              TRANSFER_B_DISPL_CM_FROM_DEVICE)(int *size, realw *displ, long *Mesh_pointer_f) {

  TRACE("transfer_b_displ_cm_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_b_displ_crust_mantle.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 displ, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_b_displ_crust_mantle.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(transfer_displ_ic_from_device,
              TRANSFER_DISPL_IC_FROM_DEVICE)(int *size, realw *displ, long *Mesh_pointer_f) {

  TRACE("transfer_displ_ic_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_displ_inner_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 displ, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_displ_inner_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(transfer_b_displ_ic_from_device,
              TRANSFER_B_DISPL_IC_FROM_DEVICE)(int *size, realw *displ, long *Mesh_pointer_f) {

  TRACE("transfer_b_displ_ic_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_b_displ_inner_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 displ, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_b_displ_inner_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(transfer_displ_oc_from_device,
              TRANSFER_DISPL_OC_FROM_DEVICE)(int *size, realw *displ, long *Mesh_pointer_f) {

  TRACE("transfer_displ_oc_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

#ifdef USE_OPENCL
  if (run_opencl) {

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_displ_outer_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 displ, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_displ_outer_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
  }
#endif
  /* ----------------------------------------------------------------------------------------------- */
}


extern EXTERN_LANG
void FC_FUNC_(transfer_b_displ_oc_from_device,
              TRANSFER_B_DISPL_OC_FROM_DEVICE)(int *size, realw *displ, long *Mesh_pointer_f) {

  TRACE("transfer_b_displ_oc_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_b_displ_outer_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 displ, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_b_displ_outer_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// velocities
/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_veloc_cm_from_device,
              TRANSFER_VELOC_CM_FROM_DEVICE)(int *size, realw *veloc, long *Mesh_pointer_f) {

  TRACE("transfer_veloc_cm_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_veloc_crust_mantle.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 veloc, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(veloc,mp->d_veloc_crust_mantle.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(transfer_veloc_ic_from_device,
              TRANSFER_VELOC_IC_FROM_DEVICE)(int *size, realw *veloc, long *Mesh_pointer_f) {

  TRACE("transfer_veloc_ic_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_veloc_inner_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 veloc, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(veloc,mp->d_veloc_inner_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
  }
#endif
}
/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(transfer_veloc_oc_from_device,
              TRANSFER_VELOC_OC_FROM_DEVICE)(int *size, realw *veloc, long *Mesh_pointer_f) {

  TRACE("transfer_veloc_oc_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_veloc_outer_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 veloc, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(veloc,mp->d_veloc_outer_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
  }
#endif
}


/* ----------------------------------------------------------------------------------------------- */
// accelerations
/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(transfer_accel_cm_to_device,
              TRANSFER_ACCEL_CM_TO_DEVICE)(int *size, realw *accel, long *Mesh_pointer_f) {

  TRACE("transfer_accel_cm_to_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_accel_crust_mantle.ocl, CL_FALSE, 0,
                                  sizeof (realw) * (*size),
                                  accel, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(mp->d_accel_crust_mantle.cuda,accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40016);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_accel_cm_from_device,
              TRANSFER_ACCEL_CM_FROM_DEVICE)(int *size, realw *accel, long *Mesh_pointer_f) {

  TRACE("transfer_accel_cm_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_accel_crust_mantle.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 accel, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(accel,mp->d_accel_crust_mantle.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40026);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */
extern EXTERN_LANG
void FC_FUNC_(transfer_b_accel_cm_from_device,
              TRANSFER_B_ACCEL_CM_FROM_DEVICE)(int *size, realw *b_accel, long *Mesh_pointer_f) {

  TRACE("transfer_b_accel_cm_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_b_accel_crust_mantle.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 b_accel, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(b_accel,mp->d_b_accel_crust_mantle.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40036);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(transfer_accel_ic_from_device,
              TRANSFER_ACCEL_IC_FROM_DEVICE)(int *size, realw *accel, long *Mesh_pointer_f) {

  TRACE("transfer_accel_ic_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_accel_inner_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 accel, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(accel,mp->d_accel_inner_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40026);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_accel_oc_from_device,
              TRANSFER_ACCEL_OC_FROM_DEVICE)(int *size, realw *accel, long *Mesh_pointer_f) {

  TRACE("transfer_accel_oc_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_accel_outer_core.ocl, CL_TRUE, 0,
                                 sizeof (realw) * (*size),
                                 accel, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(accel,mp->d_accel_outer_core.cuda,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40026);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// strain fields
/* ----------------------------------------------------------------------------------------------- */


// crust/mantle
extern EXTERN_LANG
void FC_FUNC_(transfer_strain_cm_from_device,
              TRANSFER_STRAIN_CM_FROM_DEVICE)(long *Mesh_pointer,
                                              realw *eps_trace_over_3,
                                              realw *epsilondev_xx,
                                              realw *epsilondev_yy,
                                              realw *epsilondev_xy,
                                              realw *epsilondev_xz,
                                              realw *epsilondev_yz) {
  TRACE("transfer_strain_cm_from_device");
  //get mesh pointer out of Fortran integer container

  Mesh *mp = (Mesh *) *Mesh_pointer;

  int size = NGLL3*mp->NSPEC_CRUST_MANTLE;

#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_eps_trace_over_3_crust_mantle.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 eps_trace_over_3, 0, NULL, NULL));


    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_epsilondev_xx_crust_mantle.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 epsilondev_xx, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_epsilondev_yy_crust_mantle.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 epsilondev_yy, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_epsilondev_xy_crust_mantle.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 epsilondev_xy, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_epsilondev_xz_crust_mantle.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 epsilondev_xz, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_epsilondev_yz_crust_mantle.ocl, CL_TRUE, 0,
                                 size  *sizeof (realw),
                                 epsilondev_yz, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(eps_trace_over_3,mp->d_eps_trace_over_3_crust_mantle.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),320001);

    print_CUDA_error_if_any(cudaMemcpy(epsilondev_xx,mp->d_epsilondev_xx_crust_mantle.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),320002);
    print_CUDA_error_if_any(cudaMemcpy(epsilondev_yy,mp->d_epsilondev_yy_crust_mantle.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),320003);
    print_CUDA_error_if_any(cudaMemcpy(epsilondev_xy,mp->d_epsilondev_xy_crust_mantle.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),320004);
    print_CUDA_error_if_any(cudaMemcpy(epsilondev_xz,mp->d_epsilondev_xz_crust_mantle.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),320005);
    print_CUDA_error_if_any(cudaMemcpy(epsilondev_yz,mp->d_epsilondev_yz_crust_mantle.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),320006);
  }
#endif

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_gpu_error("after transfer_strain_cm_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// crust/mantle

extern EXTERN_LANG
void FC_FUNC_(transfer_b_strain_cm_to_device,
              TRANSFER_B_STRAIN_CM_TO_DEVICE)(long *Mesh_pointer,
                                              realw *epsilondev_xx,
                                              realw *epsilondev_yy,
                                              realw *epsilondev_xy,
                                              realw *epsilondev_xz,
                                              realw *epsilondev_yz) {
  TRACE("transfer_b_strain_cm_to_device");

  // debug
  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer;

  int size = NGLL3*mp->NSPEC_CRUST_MANTLE;

  if (! mp->undo_attenuation) {
#ifdef USE_OPENCL
    if (run_opencl) {
      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_epsilondev_xx_crust_mantle.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    epsilondev_xx, 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_epsilondev_yy_crust_mantle.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    epsilondev_yy, 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_epsilondev_xy_crust_mantle.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    epsilondev_xy, 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_epsilondev_xz_crust_mantle.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    epsilondev_xz, 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_epsilondev_yz_crust_mantle.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    epsilondev_yz, 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xx_crust_mantle.cuda,epsilondev_xx,size*sizeof(realw),cudaMemcpyHostToDevice),330001);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_yy_crust_mantle.cuda,epsilondev_yy,size*sizeof(realw),cudaMemcpyHostToDevice),330002);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xy_crust_mantle.cuda,epsilondev_xy,size*sizeof(realw),cudaMemcpyHostToDevice),330003);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xz_crust_mantle.cuda,epsilondev_xz,size*sizeof(realw),cudaMemcpyHostToDevice),330004);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_yz_crust_mantle.cuda,epsilondev_yz,size*sizeof(realw),cudaMemcpyHostToDevice),330005);
    }
#endif
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_gpu_error("after transfer_b_strain_cm_to_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// inner core

extern EXTERN_LANG
void FC_FUNC_(transfer_strain_ic_from_device,
              TRANSFER_STRAIN_IC_FROM_DEVICE)(long *Mesh_pointer,
                                              realw *eps_trace_over_3,
                                              realw *epsilondev_xx,
                                              realw *epsilondev_yy,
                                              realw *epsilondev_xy,
                                              realw *epsilondev_xz,
                                              realw *epsilondev_yz) {
  TRACE("transfer_strain_ic_from_device");
  //get mesh pointer out of Fortran integer container

  Mesh *mp = (Mesh *) *Mesh_pointer;

  int size = NGLL3 * mp->NSPEC_INNER_CORE;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_eps_trace_over_3_inner_core.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 eps_trace_over_3, 0, NULL, NULL));


    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_epsilondev_xx_inner_core.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 epsilondev_xx, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_epsilondev_yy_inner_core.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 epsilondev_yy, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_epsilondev_xy_inner_core.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 epsilondev_xy, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_epsilondev_xz_inner_core.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 epsilondev_xz, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_epsilondev_yz_inner_core.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 epsilondev_yz, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(eps_trace_over_3,mp->d_eps_trace_over_3_inner_core.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),340001);

    print_CUDA_error_if_any(cudaMemcpy(epsilondev_xx,mp->d_epsilondev_xx_inner_core.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),340002);
    print_CUDA_error_if_any(cudaMemcpy(epsilondev_yy,mp->d_epsilondev_yy_inner_core.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),340003);
    print_CUDA_error_if_any(cudaMemcpy(epsilondev_xy,mp->d_epsilondev_xy_inner_core.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),340004);
    print_CUDA_error_if_any(cudaMemcpy(epsilondev_xz,mp->d_epsilondev_xz_inner_core.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),340005);
    print_CUDA_error_if_any(cudaMemcpy(epsilondev_yz,mp->d_epsilondev_yz_inner_core.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),340006);
  }
#endif

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_gpu_error("after transfer_strain_ic_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// inner_core

extern EXTERN_LANG
void FC_FUNC_(transfer_b_strain_ic_to_device,
              TRANSFER_B_STRAIN_IC_TO_DEVICE)(long *Mesh_pointer,
                                              realw *epsilondev_xx,
                                              realw *epsilondev_yy,
                                              realw *epsilondev_xy,
                                              realw *epsilondev_xz,
                                              realw *epsilondev_yz) {
  TRACE("transfer_b_strain_cm_to_device");
  // debug
  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer;

  int size = NGLL3*mp->NSPEC_INNER_CORE;

  if (! mp->undo_attenuation) {

#ifdef USE_OPENCL
    if (run_opencl) {
      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_epsilondev_xx_inner_core.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    epsilondev_xx, 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_epsilondev_yy_inner_core.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    epsilondev_yy, 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_epsilondev_xy_inner_core.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    epsilondev_xy, 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_epsilondev_xz_inner_core.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    epsilondev_xz, 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_epsilondev_yz_inner_core.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    epsilondev_yz, 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xx_inner_core.cuda,epsilondev_xx,size*sizeof(realw),cudaMemcpyHostToDevice),350001);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_yy_inner_core.cuda,epsilondev_yy,size*sizeof(realw),cudaMemcpyHostToDevice),350002);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xy_inner_core.cuda,epsilondev_xy,size*sizeof(realw),cudaMemcpyHostToDevice),350003);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xz_inner_core.cuda,epsilondev_xz,size*sizeof(realw),cudaMemcpyHostToDevice),350004);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_yz_inner_core.cuda,epsilondev_yz,size*sizeof(realw),cudaMemcpyHostToDevice),350005);
    }
#endif
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_gpu_error("after transfer_b_strain_ic_to_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// R memory variables
/* ----------------------------------------------------------------------------------------------- */

// crust/mantle

extern EXTERN_LANG
void FC_FUNC_(transfer_rmemory_cm_from_device,
              TRANSFER_RMEMORY_CM_FROM_DEVICE)(long* Mesh_pointer,
                                               realw* R_xx,
                                               realw* R_yy,
                                               realw* R_xy,
                                               realw* R_xz,
                                               realw* R_yz) {
  TRACE("transfer_rmemory_cm_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *)(*Mesh_pointer);

  int size = N_SLS * NGLL3 * mp->NSPEC_CRUST_MANTLE;

#if USE_CUDA
  if(run_cuda) {
  print_CUDA_error_if_any(cudaMemcpy(R_xx,mp->d_R_xx_crust_mantle.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),360011);
  print_CUDA_error_if_any(cudaMemcpy(R_yy,mp->d_R_yy_crust_mantle.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),360012);
  print_CUDA_error_if_any(cudaMemcpy(R_xy,mp->d_R_xy_crust_mantle.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),360013);
  print_CUDA_error_if_any(cudaMemcpy(R_xz,mp->d_R_xz_crust_mantle.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),360014);
  print_CUDA_error_if_any(cudaMemcpy(R_yz,mp->d_R_yz_crust_mantle.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),360015);
}
#endif

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_gpu_error("after transfer_rmemory_cm_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// crust/mantle

extern EXTERN_LANG
void FC_FUNC_(transfer_b_rmemory_cm_to_device,
              TRANSFER_B_RMEMORY_CM_TO_DEVICE)(long *Mesh_pointer,
                                               realw *b_R_xx,
                                               realw *b_R_yy,
                                               realw *b_R_xy,
                                               realw *b_R_xz,
                                               realw *b_R_yz) {
  TRACE("transfer_b_Rmemory_cm_to_device");

  // debug
  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer;

  int size = N_SLS*NGLL3*mp->NSPEC_CRUST_MANTLE;

  if (! mp->partial_phys_dispersion_only) {
#ifdef USE_OPENCL
    if (run_opencl) {
      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_R_xx_crust_mantle.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    b_R_xx, 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_R_yy_crust_mantle.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    b_R_yy, 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_R_xy_crust_mantle.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    b_R_xy, 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_R_xz_crust_mantle.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    b_R_xz, 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_R_yz_crust_mantle.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    b_R_yz, 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_xx_crust_mantle.cuda,b_R_xx,size*sizeof(realw),cudaMemcpyHostToDevice),360001);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_yy_crust_mantle.cuda,b_R_yy,size*sizeof(realw),cudaMemcpyHostToDevice),360002);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_xy_crust_mantle.cuda,b_R_xy,size*sizeof(realw),cudaMemcpyHostToDevice),360003);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_xz_crust_mantle.cuda,b_R_xz,size*sizeof(realw),cudaMemcpyHostToDevice),360004);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_yz_crust_mantle.cuda,b_R_yz,size*sizeof(realw),cudaMemcpyHostToDevice),360005);
    }
#endif
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_gpu_error("after transfer_b_rmemory_cm_to_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */


// inner core

extern EXTERN_LANG
void FC_FUNC_(transfer_rmemory_ic_from_device,
              TRANSFER_RMEMORY_IC_FROM_DEVICE)(long* Mesh_pointer,
                                               realw* R_xx,
                                               realw* R_yy,
                                               realw* R_xy,
                                               realw* R_xz,
                                               realw* R_yz) {
  TRACE("transfer_rmemory_cm_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = N_SLS*NGLL3*mp->NSPEC_INNER_CORE;
#ifdef USE_CUDA
    if (run_cuda) {
  print_CUDA_error_if_any(cudaMemcpy(R_xx,mp->d_R_xx_inner_core.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),370011);
  print_CUDA_error_if_any(cudaMemcpy(R_yy,mp->d_R_yy_inner_core.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),370012);
  print_CUDA_error_if_any(cudaMemcpy(R_xy,mp->d_R_xy_inner_core.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),370013);
  print_CUDA_error_if_any(cudaMemcpy(R_xz,mp->d_R_xz_inner_core.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),370014);
  print_CUDA_error_if_any(cudaMemcpy(R_yz,mp->d_R_yz_inner_core.cuda,size*sizeof(realw),cudaMemcpyDeviceToHost),370015);

}
#endif
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_gpu_error("after transfer_rmemory_ic_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// inner core

extern EXTERN_LANG
void FC_FUNC_(transfer_b_rmemory_ic_to_device,
              TRANSFER_B_RMEMORY_IC_TO_DEVICE)(long *Mesh_pointer,
                                               realw *b_R_xx,
                                               realw *b_R_yy,
                                               realw *b_R_xy,
                                               realw *b_R_xz,
                                               realw *b_R_yz) {
  TRACE("transfer_b_rmemory_ic_to_device");
  // debug

  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of Fortran integer container

  Mesh *mp = (Mesh *) *Mesh_pointer;

  int size = N_SLS*NGLL3*mp->NSPEC_INNER_CORE;

  if (! mp->partial_phys_dispersion_only) {

#ifdef USE_OPENCL
    if (run_opencl) {
      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_R_xx_inner_core.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    b_R_xx, 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_R_yy_inner_core.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    b_R_yy, 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_R_xy_inner_core.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    b_R_xy, 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_R_xz_inner_core.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    b_R_xz, 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_R_yz_inner_core.ocl, CL_FALSE, 0,
                                    size * sizeof (realw),
                                    b_R_yz, 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_xx_inner_core.cuda,b_R_xx,size*sizeof(realw),cudaMemcpyHostToDevice),370001);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_yy_inner_core.cuda,b_R_yy,size*sizeof(realw),cudaMemcpyHostToDevice),370002);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_xy_inner_core.cuda,b_R_xy,size*sizeof(realw),cudaMemcpyHostToDevice),370003);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_xz_inner_core.cuda,b_R_xz,size*sizeof(realw),cudaMemcpyHostToDevice),370004);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_yz_inner_core.cuda,b_R_yz,size*sizeof(realw),cudaMemcpyHostToDevice),370005);
    }
#endif
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_gpu_error("after transfer_b_rmemory_ic_to_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// rotation arrays
/* ----------------------------------------------------------------------------------------------- */

// for outer core

extern EXTERN_LANG
void FC_FUNC_(transfer_rotation_from_device,
              TRANSFER_ROTATION_FROM_DEVICE)(long *Mesh_pointer,
                                             realw *A_array_rotation,
                                             realw *B_array_rotation) {
  TRACE("transfer_rotation_from_device");

  //get mesh pointer out of Fortran integer container

  Mesh *mp = (Mesh *) *Mesh_pointer;

  int size = NGLL3*mp->NSPEC_OUTER_CORE;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_A_array_rotation.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 A_array_rotation, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_B_array_rotation.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 B_array_rotation, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(A_array_rotation,mp->d_A_array_rotation.cuda,
                                       size*sizeof(realw),cudaMemcpyDeviceToHost),380001);
    print_CUDA_error_if_any(cudaMemcpy(B_array_rotation,mp->d_B_array_rotation.cuda,
                                       size*sizeof(realw),cudaMemcpyDeviceToHost),380002);

  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// for outer core

extern EXTERN_LANG
void FC_FUNC_(transfer_b_rotation_to_device,
              TRANSFER_B_ROTATION_TO_DEVICE)(long *Mesh_pointer,
                                             realw *A_array_rotation,
                                             realw *B_array_rotation) {
  TRACE("transfer_b_rotation_to_device");
  // debug

  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of Fortran integer container

  Mesh *mp = (Mesh *) *Mesh_pointer;

  int size = NGLL3*mp->NSPEC_OUTER_CORE;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_A_array_rotation.ocl, CL_FALSE, 0,
                                  size * sizeof (realw),
                                  A_array_rotation, 0, NULL, NULL));

    clCheck (clEnqueueWriteBuffer(mocl.command_queue, mp->d_b_B_array_rotation.ocl, CL_FALSE, 0,
                                  size * sizeof (realw),
                                  B_array_rotation, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_A_array_rotation.cuda,A_array_rotation,
                                       size*sizeof(realw),cudaMemcpyHostToDevice),390001);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_B_array_rotation.cuda,B_array_rotation,
                                       size*sizeof(realw),cudaMemcpyHostToDevice),390002);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// KERNEL transfers
/* ----------------------------------------------------------------------------------------------- */

// crust/mantle

extern EXTERN_LANG
void FC_FUNC_(transfer_kernels_cm_to_host,
              TRANSFER_KERNELS_CM_TO_HOST)(long *Mesh_pointer,
                                           realw *h_rho_kl,
                                           realw *h_alpha_kl,
                                           realw *h_beta_kl,
                                           realw *h_cijkl_kl,
                                           int *NSPEC) {
  TRACE("transfer_kernels_cm_to_host");

  //get mesh pointer out of Fortran integer container

  Mesh *mp = (Mesh *) *Mesh_pointer;

  int size = (*NSPEC)*NGLL3;

  // density kernel


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_rho_kl_crust_mantle.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 h_rho_kl, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(h_rho_kl,mp->d_rho_kl_crust_mantle.cuda,
                                       size*sizeof(realw),cudaMemcpyDeviceToHost),40101);
  }
#endif

  if (! mp->anisotropic_kl) {
    // isotropic kernels
#ifdef USE_OPENCL
    if (run_opencl) {
      clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_alpha_kl_crust_mantle.ocl, CL_TRUE, 0,
                                   size * sizeof (realw),
                                   h_alpha_kl, 0, NULL, NULL));

      clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_beta_kl_crust_mantle.ocl, CL_TRUE, 0,
                                   size * sizeof (realw),
                                   h_beta_kl, 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      print_CUDA_error_if_any(cudaMemcpy(h_alpha_kl,mp->d_alpha_kl_crust_mantle.cuda,
                                         size*sizeof(realw),cudaMemcpyDeviceToHost),40102);
      print_CUDA_error_if_any(cudaMemcpy(h_beta_kl,mp->d_beta_kl_crust_mantle.cuda,
                                         size*sizeof(realw),cudaMemcpyDeviceToHost),40103);
    }
#endif
  } else {
    // anisotropic kernels

#ifdef USE_OPENCL
    if (run_opencl) {
      clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_cijkl_kl_crust_mantle.ocl, CL_TRUE, 0,
                                   21 * size * sizeof (realw),
                                   h_cijkl_kl, 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      print_CUDA_error_if_any(cudaMemcpy(h_cijkl_kl,mp->d_cijkl_kl_crust_mantle.cuda,
                                         21*size*sizeof(realw),cudaMemcpyDeviceToHost),40102);

    }
#endif
  }
}

/* ----------------------------------------------------------------------------------------------- */
// inner core

extern EXTERN_LANG
void FC_FUNC_(transfer_kernels_ic_to_host,
              TRANSFER_KERNELS_IC_TO_HOST)(long *Mesh_pointer,
                                           realw *h_rho_kl,
                                           realw *h_alpha_kl,
                                           realw *h_beta_kl,
                                           int *NSPEC) {
  TRACE("transfer_kernels_ic_to_host");

  //get mesh pointer out of Fortran integer container

  Mesh *mp = (Mesh *) *Mesh_pointer;

  int size = (*NSPEC) * NGLL3;


#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_rho_kl_inner_core.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 h_rho_kl, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_alpha_kl_inner_core.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 h_alpha_kl, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_beta_kl_inner_core.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 h_beta_kl, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(h_rho_kl,mp->d_rho_kl_inner_core.cuda,
                                       size*sizeof(realw),cudaMemcpyDeviceToHost),40101);
    print_CUDA_error_if_any(cudaMemcpy(h_alpha_kl,mp->d_alpha_kl_inner_core.cuda,
                                       size*sizeof(realw),cudaMemcpyDeviceToHost),40102);
    print_CUDA_error_if_any(cudaMemcpy(h_beta_kl,mp->d_beta_kl_inner_core.cuda,
                                       size*sizeof(realw),cudaMemcpyDeviceToHost),40103);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// outer core

extern EXTERN_LANG
void FC_FUNC_(transfer_kernels_oc_to_host,
              TRANSFER_KERNELS_OC_TO_HOST)(long *Mesh_pointer,
                                           realw *h_rho_kl,
                                           realw *h_alpha_kl,
                                           int *NSPEC) {

  TRACE("transfer_kernels_oc_to_host");

  //get mesh pointer out of Fortran integer container

  Mesh *mp = (Mesh *) *Mesh_pointer;

  int size = (*NSPEC)*NGLL3;

  // copies kernel values over to CPU host

#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_rho_kl_outer_core.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 h_rho_kl, 0, NULL, NULL));

    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_alpha_kl_outer_core.ocl, CL_TRUE, 0,
                                 size * sizeof (realw),
                                 h_alpha_kl, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(h_rho_kl,mp->d_rho_kl_outer_core.cuda,
                                       size*sizeof(realw),cudaMemcpyDeviceToHost),54101);
    print_CUDA_error_if_any(cudaMemcpy(h_alpha_kl,mp->d_alpha_kl_outer_core.cuda,
                                       size*sizeof(realw),cudaMemcpyDeviceToHost),54102);
  }
#endif
}

/* ----------------------------------------------------------------------------------------------- */
// for NOISE simulations

extern EXTERN_LANG
void FC_FUNC_(transfer_kernels_noise_to_host,
              TRANSFER_KERNELS_NOISE_TO_HOST)(long *Mesh_pointer,
                                              realw *h_Sigma_kl,
                                              int *NSPEC) {
  TRACE("transfer_kernels_noise_to_host");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer;

#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_Sigma_kl.ocl, CL_TRUE, 0,
                                 NGLL3 * (*NSPEC) * sizeof (realw),
                                 h_Sigma_kl, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(h_Sigma_kl,mp->d_Sigma_kl.cuda,NGLL3*(*NSPEC)*sizeof(realw),
                                       cudaMemcpyDeviceToHost),40201);
  }
#endif
}
/* ----------------------------------------------------------------------------------------------- */


// for Hess kernel calculations

extern EXTERN_LANG
void FC_FUNC_(transfer_kernels_hess_cm_tohost,
              TRANSFER_KERNELS_HESS_CM_TOHOST)(long *Mesh_pointer,
                                               realw *h_hess_kl,
                                               int *NSPEC) {
  TRACE("transfer_kernels_hess_cm_tohost");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer;

#ifdef USE_OPENCL
  if (run_opencl) {
    clCheck (clEnqueueReadBuffer(mocl.command_queue, mp->d_hess_kl_crust_mantle.ocl, CL_TRUE, 0,
                                 NGLL3 * (*NSPEC) * sizeof (realw),
                                 h_hess_kl, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    print_CUDA_error_if_any(cudaMemcpy(h_hess_kl,mp->d_hess_kl_crust_mantle.cuda,NGLL3*(*NSPEC)*sizeof(realw),
                                       cudaMemcpyDeviceToHost),70201);
  }
#endif
}
