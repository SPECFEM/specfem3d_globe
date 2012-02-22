/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  2 . 0
 !               ---------------------------------------
 !
 !          Main authors: Dimitri Komatitsch and Jeroen Tromp
 !    Princeton University, USA and University of Pau / CNRS / INRIA
 ! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
 !                            April 2011
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
#include <cuda.h>
#include <cublas.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"
#include "prepare_constants_cuda.h"

/* ----------------------------------------------------------------------------------------------- */

// Transfer functions

/* ----------------------------------------------------------------------------------------------- */



/* ----------------------------------------------------------------------------------------------- */

// transfer memory from CPU host to GPU device

/* ----------------------------------------------------------------------------------------------- */

// crust_mantle
extern "C"
void FC_FUNC_(transfer_fields_cm_to_device,
              TRANSFER_FIELDS_CM_TO_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer_f) {

TRACE("transfer_fields_cm_to_device_");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(mp->d_displ_crust_mantle,displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40003);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_veloc_crust_mantle,veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40004);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_accel_crust_mantle,accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40005);

}

// inner_core
extern "C"
void FC_FUNC_(transfer_fields_ic_to_device,
              TRANSFER_FIELDS_IC_TO_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer_f) {

  TRACE("transfer_fields_ic_to_device_");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(mp->d_displ_inner_core,displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40003);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_veloc_inner_core,veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40004);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_accel_inner_core,accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40005);

}

// outer_core
extern "C"
void FC_FUNC_(transfer_fields_oc_to_device,
              TRANSFER_FIELDS_OC_TO_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer_f) {

  TRACE("transfer_fields_oc_to_device_");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(mp->d_displ_outer_core,displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40003);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_veloc_outer_core,veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40004);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_accel_outer_core,accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40005);

}

/* ----------------------------------------------------------------------------------------------- */

// backward/reconstructed fields

// crust_mantle
extern "C"
void FC_FUNC_(transfer_b_fields_cm_to_device,
              TRANSFER_FIELDS_B_CM_TO_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                              long* Mesh_pointer_f) {

  TRACE("transfer_fields_b_cm_to_device_");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container
  cudaMemcpy(mp->d_b_displ_crust_mantle,b_displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_veloc_crust_mantle,b_veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_accel_crust_mantle,b_accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice);

}

// inner_core
extern "C"
void FC_FUNC_(transfer_b_fields_ic_to_device,
              TRANSFER_FIELDS_B_IC_TO_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                              long* Mesh_pointer_f) {

  TRACE("transfer_fields_b_ic_to_device_");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container
  cudaMemcpy(mp->d_b_displ_inner_core,b_displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_veloc_inner_core,b_veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_accel_inner_core,b_accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice);

}

// outer_core
extern "C"
void FC_FUNC_(transfer_b_fields_oc_to_device,
              TRANSFER_FIELDS_B_OC_TO_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                              long* Mesh_pointer_f) {

  TRACE("transfer_fields_b_oc_to_device_");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container
  cudaMemcpy(mp->d_b_displ_outer_core,b_displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_veloc_outer_core,b_veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_accel_outer_core,b_accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice);

}



/* ----------------------------------------------------------------------------------------------- */

// transfer memory from GPU device to CPU host

/* ----------------------------------------------------------------------------------------------- */

// crust_mantle
extern "C"
void FC_FUNC_(transfer_fields_cm_from_device,
              TRANSFER_FIELDS_CM_FROM_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer_f) {

  TRACE("transfer_fields_cm_from_device_");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_displ_crust_mantle,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
  print_CUDA_error_if_any(cudaMemcpy(veloc,mp->d_veloc_crust_mantle,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
  print_CUDA_error_if_any(cudaMemcpy(accel,mp->d_accel_crust_mantle,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40008);

}

/* ----------------------------------------------------------------------------------------------- */

// inner_core
extern "C"
void FC_FUNC_(transfer_fields_ic_from_device,
              TRANSFER_FIELDS_IC_FROM_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer_f) {

  TRACE("transfer_fields_ic_from_device_");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_displ_inner_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
  print_CUDA_error_if_any(cudaMemcpy(veloc,mp->d_veloc_inner_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
  print_CUDA_error_if_any(cudaMemcpy(accel,mp->d_accel_inner_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40008);

}

/* ----------------------------------------------------------------------------------------------- */

// outer_core
extern "C"
void FC_FUNC_(transfer_fields_oc_from_device,
              TRANSFER_FIELDS_OC_FROM_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer_f) {

  TRACE("transfer_fields_oc_from_device_");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_displ_outer_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
  print_CUDA_error_if_any(cudaMemcpy(veloc,mp->d_veloc_outer_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
  print_CUDA_error_if_any(cudaMemcpy(accel,mp->d_accel_outer_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40008);

}



/* ----------------------------------------------------------------------------------------------- */

// crust_mantle
extern "C"
void FC_FUNC_(transfer_b_fields_cm_from_device,
              TRANSFER_B_FIELDS_CM_FROM_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                                long* Mesh_pointer_f) {

TRACE("transfer_b_fields_cm_from_device_");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  cudaMemcpy(b_displ,mp->d_b_displ_crust_mantle,sizeof(realw)*(*size),cudaMemcpyDeviceToHost);
  cudaMemcpy(b_veloc,mp->d_b_veloc_crust_mantle,sizeof(realw)*(*size),cudaMemcpyDeviceToHost);
  cudaMemcpy(b_accel,mp->d_b_accel_crust_mantle,sizeof(realw)*(*size),cudaMemcpyDeviceToHost);

}

/* ----------------------------------------------------------------------------------------------- */

// inner_core
extern "C"
void FC_FUNC_(transfer_b_fields_ic_from_device,
              TRANSFER_B_FIELDS_IC_FROM_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                                long* Mesh_pointer_f) {

  TRACE("transfer_fields_b_ic_from_device_");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(b_displ,mp->d_b_displ_inner_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
  print_CUDA_error_if_any(cudaMemcpy(b_veloc,mp->d_b_veloc_inner_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
  print_CUDA_error_if_any(cudaMemcpy(b_accel,mp->d_b_accel_inner_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40008);

}

/* ----------------------------------------------------------------------------------------------- */

// outer_core
extern "C"
void FC_FUNC_(transfer_b_fields_oc_from_device,
              TRANSFER_B_FIELDS_OC_FROM_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                                long* Mesh_pointer_f) {

  TRACE("transfer_b_fields_oc_from_device_");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(b_displ,mp->d_b_displ_outer_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
  print_CUDA_error_if_any(cudaMemcpy(b_veloc,mp->d_b_veloc_outer_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
  print_CUDA_error_if_any(cudaMemcpy(b_accel,mp->d_b_accel_outer_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40008);

}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_accel_cm_to_device,
              TRNASFER_ACCEL_CM_TO_DEVICE)(int* size, realw* accel,long* Mesh_pointer_f) {

TRACE("transfer_accel_cm_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(mp->d_accel_crust_mantle,accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40016);

}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_displ_cm_from_device,
              TRANSFER_DISPL_CM_FROM_DEVICE)(int* size, realw* displ, long* Mesh_pointer_f) {

  TRACE("transfer_displ_cm_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_displ_crust_mantle,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_b_displ_cm_from_device,
              TRANSFER_B_DISPL_CM_FROM_DEVICE)(int* size, realw* displ, long* Mesh_pointer_f) {

  TRACE("transfer_b_displ_cm_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_b_displ_crust_mantle,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_displ_ic_from_device,
              TRANSFER_DISPL_IC_FROM_DEVICE)(int* size, realw* displ, long* Mesh_pointer_f) {

  TRACE("transfer_displ_ic_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_displ_inner_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_b_displ_ic_from_device,
              TRANSFER_B_DISPL_IC_FROM_DEVICE)(int* size, realw* displ, long* Mesh_pointer_f) {

  TRACE("transfer_b_displ_ic_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_b_displ_inner_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_veloc_cm_from_device,
              TRANSFER_DISPL_CM_FROM_DEVICE)(int* size, realw* veloc, long* Mesh_pointer_f) {

  TRACE("transfer_veloc_cm_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(veloc,mp->d_veloc_crust_mantle,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_accel_cm_from_device,
              TRANSFER_ACCEL_CM_FROM_DEVICE)(int* size, realw* accel,long* Mesh_pointer_f) {

  TRACE("transfer_accel_cm_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(accel,mp->d_accel_crust_mantle,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40026);

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_b_accel_cm_from_device,
              TRNASFER_B_ACCEL_CM_FROM_DEVICE)(int* size, realw* b_accel,long* Mesh_pointer_f) {

TRACE("transfer_b_accel_cm_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(b_accel,mp->d_b_accel_crust_mantle,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40036);

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_accel_ic_from_device,
              TRANSFER_ACCEL_IC_FROM_DEVICE)(int* size, realw* accel,long* Mesh_pointer_f) {

  TRACE("transfer_accel_ic_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(accel,mp->d_accel_inner_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40026);

}

/* ----------------------------------------------------------------------------------------------- */

// strain fields

/* ----------------------------------------------------------------------------------------------- */

// crust/mantle
extern "C"
void FC_FUNC_(transfer_strain_cm_from_device,
              TRANSFER_STRAIN_CM_FROM_DEVICE)(long* Mesh_pointer,
                                                  realw* eps_trace_over_3,
                                                  realw* epsilondev_xx,
                                                  realw* epsilondev_yy,
                                                  realw* epsilondev_xy,
                                                  realw* epsilondev_xz,
                                                  realw* epsilondev_yz) {
  TRACE("transfer_strain_cm_from_device");
  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = NGLL3*mp->NSPEC_CRUST_MANTLE;

  cudaMemcpy(eps_trace_over_3,mp->d_eps_trace_over_3_crust_mantle,size*sizeof(realw),cudaMemcpyDeviceToHost);

  cudaMemcpy(epsilondev_xx,mp->d_epsilondev_xx_crust_mantle,size*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_yy,mp->d_epsilondev_yy_crust_mantle,size*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_xy,mp->d_epsilondev_xy_crust_mantle,size*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_xz,mp->d_epsilondev_xz_crust_mantle,size*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_yz,mp->d_epsilondev_yz_crust_mantle,size*sizeof(realw),cudaMemcpyDeviceToHost);


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_strain_cm_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// crust/mantle

extern "C"
void FC_FUNC_(transfer_b_strain_cm_to_device,
              TRANSFER_B_STRAIN_CM_TO_DEVICE)(long* Mesh_pointer,
                                              realw* epsilondev_xx,
                                              realw* epsilondev_yy,
                                              realw* epsilondev_xy,
                                              realw* epsilondev_xz,
                                              realw* epsilondev_yz) {
  TRACE("transfer_b_strain_cm_to_device");
  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = NGLL3*mp->NSPEC_CRUST_MANTLE;

  cudaMemcpy(mp->d_b_epsilondev_xx_crust_mantle,epsilondev_xx,size*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_epsilondev_yy_crust_mantle,epsilondev_yy,size*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_epsilondev_xy_crust_mantle,epsilondev_xy,size*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_epsilondev_xz_crust_mantle,epsilondev_xz,size*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_epsilondev_yz_crust_mantle,epsilondev_yz,size*sizeof(realw),cudaMemcpyHostToDevice);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_b_strain_cm_to_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// inner core

extern "C"
void FC_FUNC_(transfer_strain_ic_from_device,
              TRANSFER_STRAIN_IC_FROM_DEVICE)(long* Mesh_pointer,
                                              realw* eps_trace_over_3,
                                              realw* epsilondev_xx,
                                              realw* epsilondev_yy,
                                              realw* epsilondev_xy,
                                              realw* epsilondev_xz,
                                              realw* epsilondev_yz) {
  TRACE("transfer_strain_ic_from_device");
  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = NGLL3*mp->NSPEC_INNER_CORE;

  cudaMemcpy(eps_trace_over_3,mp->d_eps_trace_over_3_inner_core,size*sizeof(realw),cudaMemcpyDeviceToHost);

  cudaMemcpy(epsilondev_xx,mp->d_epsilondev_xx_inner_core,size*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_yy,mp->d_epsilondev_yy_inner_core,size*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_xy,mp->d_epsilondev_xy_inner_core,size*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_xz,mp->d_epsilondev_xz_inner_core,size*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_yz,mp->d_epsilondev_yz_inner_core,size*sizeof(realw),cudaMemcpyDeviceToHost);


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_strain_ic_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// inner_core

extern "C"
void FC_FUNC_(transfer_b_strain_ic_to_device,
              TRANSFER_B_STRAIN_IC_TO_DEVICE)(long* Mesh_pointer,
                                              realw* epsilondev_xx,
                                              realw* epsilondev_yy,
                                              realw* epsilondev_xy,
                                              realw* epsilondev_xz,
                                              realw* epsilondev_yz) {
  TRACE("transfer_b_strain_cm_to_device");
  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = NGLL3*mp->NSPEC_INNER_CORE;

  cudaMemcpy(mp->d_b_epsilondev_xx_inner_core,epsilondev_xx,size*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_epsilondev_yy_inner_core,epsilondev_yy,size*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_epsilondev_xy_inner_core,epsilondev_xy,size*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_epsilondev_xz_inner_core,epsilondev_xz,size*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_epsilondev_yz_inner_core,epsilondev_yz,size*sizeof(realw),cudaMemcpyHostToDevice);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_b_strain_ic_to_device");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// rotation arrays

/* ----------------------------------------------------------------------------------------------- */

// for outer core

extern "C"
void FC_FUNC_(transfer_rotation_from_device,
              TRANSFER_ROTATION_FROM_DEVICE)(long* Mesh_pointer,
                                             realw* A_array_rotation,
                                             realw* B_array_rotation) {
  TRACE("transfer_rotation_from_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = NGLL3*mp->NSPEC_OUTER_CORE;

  cudaMemcpy(A_array_rotation,mp->d_A_array_rotation,size*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(B_array_rotation,mp->d_B_array_rotation,size*sizeof(realw),cudaMemcpyDeviceToHost);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_rotation_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// for outer core

extern "C"
void FC_FUNC_(transfer_b_rotation_to_device,
              TRANSFER_B_ROTATION_TO_DEVICE)(long* Mesh_pointer,
                                              realw* A_array_rotation,
                                              realw* B_array_rotation) {
  TRACE("transfer_b_rotation_to_device");
  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = NGLL3*mp->NSPEC_OUTER_CORE;

  cudaMemcpy(mp->d_b_A_array_rotation,A_array_rotation,size*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_B_array_rotation,B_array_rotation,size*sizeof(realw),cudaMemcpyHostToDevice);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_b_rotation_to_device");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// attenuation fields

/* ----------------------------------------------------------------------------------------------- */

/*
// feature not used so far ...
// crust_mantle
extern "C"
void FC_FUNC_(transfer_b_att_cm_to_device,
              TRANSFER_B_ATT_CM_TO_DEVICE)(long* Mesh_pointer,
                                           realw* R_xx,
                                           realw* R_yy,
                                           realw* R_xy,
                                           realw* R_xz,
                                           realw* R_yz) {
  TRACE("transfer_b_att_cm_to_device");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  if( ! mp->use_attenuation_mimic){ exit_on_cuda_error("not supported attenuation feature yet");}

  // not used so far...
  // see notes about USE_ATTENUATION_MIMIC
  int size = N_SLS*NGLL3*mp->NSPEC_CRUST_MANTLE;

  cudaMemcpy(mp->d_b_R_xx_crust_mantle,R_xx,size*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_R_yy_crust_mantle,R_yy,size*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_R_xy_crust_mantle,R_xy,size*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_R_xz_crust_mantle,R_xz,size*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_R_yz_crust_mantle,R_yz,size*sizeof(realw),cudaMemcpyHostToDevice);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_b_att_cm_to_device");
#endif
}
*/









//daniel: TODO old code routines...



/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_sigma_from_device,
              TRANSFER_SIGMA_FROM_DEVICE)(int* size, realw* sigma_kl,long* Mesh_pointer_f) {

TRACE("transfer_sigma_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(sigma_kl,mp->d_Sigma_kl,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40046);

}

/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(transfer_b_displ_from_device,
              TRANSFER_B_DISPL_FROM_DEVICE)(int* size, realw* displ,long* Mesh_pointer_f) {

TRACE("transfer_b_displ_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_displ,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40056);

}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(transfer_displ_from_device,
              TRANSFER_DISPL_FROM_DEVICE)(int* size, realw* displ,long* Mesh_pointer_f) {

TRACE("transfer_displ_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_displ,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40066);

}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(transfer_compute_kernel_answers_from_device,
              TRANSFER_COMPUTE_KERNEL_ANSWERS_FROM_DEVICE)(long* Mesh_pointer,
                                                           realw* rho_kl,int* size_rho,
                                                           realw* mu_kl, int* size_mu,
                                                           realw* kappa_kl, int* size_kappa) {
TRACE("transfer_compute_kernel_answers_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  cudaMemcpy(rho_kl,mp->d_rho_kl,*size_rho*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(mu_kl,mp->d_mu_kl,*size_mu*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(kappa_kl,mp->d_kappa_kl,*size_kappa*sizeof(realw),cudaMemcpyDeviceToHost);

}
*/

/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(transfer_compute_kernel_fields_from_device,
              TRANSFER_COMPUTE_KERNEL_FIELDS_FROM_DEVICE)(long* Mesh_pointer,
                                                          realw* accel, int* size_accel,
                                                          realw* b_displ, int* size_b_displ,
                                                          realw* epsilondev_xx,
                                                          realw* epsilondev_yy,
                                                          realw* epsilondev_xy,
                                                          realw* epsilondev_xz,
                                                          realw* epsilondev_yz,
                                                          int* size_epsilondev,
                                                          realw* b_epsilondev_xx,
                                                          realw* b_epsilondev_yy,
                                                          realw* b_epsilondev_xy,
                                                          realw* b_epsilondev_xz,
                                                          realw* b_epsilondev_yz,
                                                          int* size_b_epsilondev,
                                                          realw* rho_kl,int* size_rho,
                                                          realw* mu_kl, int* size_mu,
                                                          realw* kappa_kl, int* size_kappa,
                                                          realw* epsilon_trace_over_3,
                                                          realw* b_epsilon_trace_over_3,
                                                          int* size_epsilon_trace_over_3) {
TRACE("transfer_compute_kernel_fields_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  cudaMemcpy(accel,mp->d_accel,*size_accel*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(b_displ,mp->d_b_displ,*size_b_displ*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_xx,mp->d_epsilondev_xx,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_yy,mp->d_epsilondev_yy,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_xy,mp->d_epsilondev_xy,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_xz,mp->d_epsilondev_xz,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_yz,mp->d_epsilondev_yz,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(b_epsilondev_xx,mp->d_b_epsilondev_xx,*size_b_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(b_epsilondev_yy,mp->d_b_epsilondev_yy,*size_b_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(b_epsilondev_xy,mp->d_b_epsilondev_xy,*size_b_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(b_epsilondev_xz,mp->d_b_epsilondev_xz,*size_b_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(b_epsilondev_yz,mp->d_b_epsilondev_yz,*size_b_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(rho_kl,mp->d_rho_kl,*size_rho*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(mu_kl,mp->d_mu_kl,*size_mu*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(kappa_kl,mp->d_kappa_kl,*size_kappa*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilon_trace_over_3,mp->d_epsilon_trace_over_3,*size_epsilon_trace_over_3*sizeof(realw),
       cudaMemcpyDeviceToHost);
  cudaMemcpy(b_epsilon_trace_over_3,mp->d_b_epsilon_trace_over_3,*size_epsilon_trace_over_3*sizeof(realw),
       cudaMemcpyDeviceToHost);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_compute_kernel_fields_from_device");
#endif
}
*/

/* ----------------------------------------------------------------------------------------------- */
/*
// attenuation fields

extern "C"
void FC_FUNC_(transfer_b_fields_att_to_device,
              TRANSFER_B_FIELDS_ATT_TO_DEVICE)(long* Mesh_pointer,
                                             realw* b_R_xx,realw* b_R_yy,realw* b_R_xy,realw* b_R_xz,realw* b_R_yz,
                                             int* size_R,
                                             realw* b_epsilondev_xx,
                                             realw* b_epsilondev_yy,
                                             realw* b_epsilondev_xy,
                                             realw* b_epsilondev_xz,
                                             realw* b_epsilondev_yz,
                                             int* size_epsilondev) {
  TRACE("transfer_b_fields_att_to_device");
  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  cudaMemcpy(mp->d_b_R_xx,b_R_xx,*size_R*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_R_yy,b_R_yy,*size_R*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_R_xy,b_R_xy,*size_R*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_R_xz,b_R_xz,*size_R*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_R_yz,b_R_yz,*size_R*sizeof(realw),cudaMemcpyHostToDevice);

  cudaMemcpy(mp->d_b_epsilondev_xx,b_epsilondev_xx,*size_epsilondev*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_epsilondev_yy,b_epsilondev_yy,*size_epsilondev*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_epsilondev_xy,b_epsilondev_xy,*size_epsilondev*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_epsilondev_xz,b_epsilondev_xz,*size_epsilondev*sizeof(realw),cudaMemcpyHostToDevice);
  cudaMemcpy(mp->d_b_epsilondev_yz,b_epsilondev_yz,*size_epsilondev*sizeof(realw),cudaMemcpyHostToDevice);


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_b_fields_att_to_device");
#endif
}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
// attenuation fields

extern "C"
void FC_FUNC_(transfer_fields_att_from_device,
              TRANSFER_FIELDS_ATT_FROM_DEVICE)(long* Mesh_pointer,
                                               realw* R_xx,realw* R_yy,realw* R_xy,realw* R_xz,realw* R_yz,
                                               int* size_R,
                                               realw* epsilondev_xx,
                                               realw* epsilondev_yy,
                                               realw* epsilondev_xy,
                                               realw* epsilondev_xz,
                                               realw* epsilondev_yz,
                                               int* size_epsilondev) {
  TRACE("transfer_fields_att_from_device");
  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  cudaMemcpy(R_xx,mp->d_R_xx,*size_R*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(R_yy,mp->d_R_yy,*size_R*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(R_xy,mp->d_R_xy,*size_R*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(R_xz,mp->d_R_xz,*size_R*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(R_yz,mp->d_R_yz,*size_R*sizeof(realw),cudaMemcpyDeviceToHost);

  cudaMemcpy(epsilondev_xx,mp->d_epsilondev_xx,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_yy,mp->d_epsilondev_yy,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_xy,mp->d_epsilondev_xy,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_xz,mp->d_epsilondev_xz,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(epsilondev_yz,mp->d_epsilondev_yz,*size_epsilondev*sizeof(realw),cudaMemcpyDeviceToHost);


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_fields_att_from_device");
#endif
}
*/

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_kernels_el_to_host,
              TRANSFER_KERNELS_EL_TO_HOST)(long* Mesh_pointer,
                                                    realw* h_rho_kl,
                                                    realw* h_mu_kl,
                                                    realw* h_kappa_kl,
                                                    int* NSPEC_AB) {
TRACE("transfer_kernels_el_to_host");
  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  print_CUDA_error_if_any(cudaMemcpy(h_rho_kl,mp->d_rho_kl,*NSPEC_AB*NGLL3*sizeof(realw),
                                     cudaMemcpyDeviceToHost),40101);
  print_CUDA_error_if_any(cudaMemcpy(h_mu_kl,mp->d_mu_kl,*NSPEC_AB*NGLL3*sizeof(realw),
                                     cudaMemcpyDeviceToHost),40102);
  print_CUDA_error_if_any(cudaMemcpy(h_kappa_kl,mp->d_kappa_kl,*NSPEC_AB*NGLL3*sizeof(realw),
                                     cudaMemcpyDeviceToHost),40103);

}

/* ----------------------------------------------------------------------------------------------- */

// for NOISE simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_kernels_noise_to_host,
              TRANSFER_KERNELS_NOISE_TO_HOST)(long* Mesh_pointer,
                                                          realw* h_Sigma_kl,
                                                          int* NSPEC_AB) {
TRACE("transfer_kernels_noise_to_host");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(h_Sigma_kl,mp->d_Sigma_kl,NGLL3*(*NSPEC_AB)*sizeof(realw),
                                     cudaMemcpyDeviceToHost),40201);

}


/* ----------------------------------------------------------------------------------------------- */

// for ACOUSTIC simulations

/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(transfer_fields_ac_to_device,
              TRANSFER_FIELDS_AC_TO_DEVICE)(
                                                  int* size,
                                                  realw* potential_acoustic,
                                                  realw* potential_dot_acoustic,
                                                  realw* potential_dot_dot_acoustic,
                                                  long* Mesh_pointer_f) {
TRACE("transfer_fields_ac_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(mp->d_potential_acoustic,potential_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyHostToDevice),50110);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_potential_dot_acoustic,potential_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyHostToDevice),50120);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_potential_dot_dot_acoustic,potential_dot_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyHostToDevice),50130);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_fields_ac_to_device");
#endif
}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(transfer_b_fields_ac_to_device,
              TRANSFER_B_FIELDS_AC_TO_DEVICE)(
                                                    int* size,
                                                    realw* b_potential_acoustic,
                                                    realw* b_potential_dot_acoustic,
                                                    realw* b_potential_dot_dot_acoustic,
                                                    long* Mesh_pointer_f) {
TRACE("transfer_b_fields_ac_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_potential_acoustic,b_potential_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyHostToDevice),51110);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_potential_dot_acoustic,b_potential_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyHostToDevice),51120);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_potential_dot_dot_acoustic,b_potential_dot_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyHostToDevice),51130);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_b_fields_ac_to_device");
#endif
}
*/

/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(transfer_fields_ac_from_device,TRANSFER_FIELDS_AC_FROM_DEVICE)(
                                                                                         int* size,
                                                                                         realw* potential_acoustic,
                                                                                         realw* potential_dot_acoustic,
                                                                                         realw* potential_dot_dot_acoustic,
                                                                                         long* Mesh_pointer_f) {
TRACE("transfer_fields_ac_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(potential_acoustic,mp->d_potential_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyDeviceToHost),52111);
  print_CUDA_error_if_any(cudaMemcpy(potential_dot_acoustic,mp->d_potential_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyDeviceToHost),52121);
  print_CUDA_error_if_any(cudaMemcpy(potential_dot_dot_acoustic,mp->d_potential_dot_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyDeviceToHost),52131);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_fields_ac_from_device");
#endif
}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(transfer_b_fields_ac_from_device,
              TRANSFER_B_FIELDS_AC_FROM_DEVICE)(
                                                      int* size,
                                                      realw* b_potential_acoustic,
                                                      realw* b_potential_dot_acoustic,
                                                      realw* b_potential_dot_dot_acoustic,
                                                      long* Mesh_pointer_f) {
TRACE("transfer_b_fields_ac_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(b_potential_acoustic,mp->d_b_potential_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyDeviceToHost),53111);
  print_CUDA_error_if_any(cudaMemcpy(b_potential_dot_acoustic,mp->d_b_potential_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyDeviceToHost),53121);
  print_CUDA_error_if_any(cudaMemcpy(b_potential_dot_dot_acoustic,mp->d_b_potential_dot_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyDeviceToHost),53131);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_b_fields_ac_from_device");
#endif
}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(transfer_dot_dot_from_device,
              TRNASFER_DOT_DOT_FROM_DEVICE)(int* size, realw* potential_dot_dot_acoustic,long* Mesh_pointer_f) {

  TRACE("transfer_dot_dot_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(potential_dot_dot_acoustic,mp->d_potential_dot_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyDeviceToHost),50041);

}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(transfer_b_dot_dot_from_device,
              TRNASFER_B_DOT_DOT_FROM_DEVICE)(int* size, realw* b_potential_dot_dot_acoustic,long* Mesh_pointer_f) {

  TRACE("transfer_b_dot_dot_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(b_potential_dot_dot_acoustic,mp->d_b_potential_dot_dot_acoustic,
                                     sizeof(realw)*(*size),cudaMemcpyDeviceToHost),50042);

}
*/

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_kernels_ac_to_host,
              TRANSFER_KERNELS_AC_TO_HOST)(long* Mesh_pointer,
                                                             realw* h_rho_ac_kl,
                                                             realw* h_kappa_ac_kl,
                                                             int* NSPEC_AB) {

  TRACE("transfer_kernels_ac_to_host");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);
  int size = *NSPEC_AB*NGLL3;

  // copies kernel values over to CPU host
  print_CUDA_error_if_any(cudaMemcpy(h_rho_ac_kl,mp->d_rho_ac_kl,size*sizeof(realw),
                                     cudaMemcpyDeviceToHost),54101);
  print_CUDA_error_if_any(cudaMemcpy(h_kappa_ac_kl,mp->d_kappa_ac_kl,size*sizeof(realw),
                                     cudaMemcpyDeviceToHost),54102);
}

/* ----------------------------------------------------------------------------------------------- */

// for Hess kernel calculations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_kernels_hess_cm_tohost,
              TRANSFER_KERNELS_HESS_CM_TOHOST)(long* Mesh_pointer,
                                              realw* h_hess_kl,
                                              int* NSPEC) {
TRACE("transfer_kernels_hess_cm_tohost");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(h_hess_kl,mp->d_hess_kl_crust_mantle,NGLL3*(*NSPEC)*sizeof(realw),
                                     cudaMemcpyDeviceToHost),70201);
}

/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(transfer_kernels_hess_ac_tohost,
              TRANSFER_KERNELS_HESS_AC_TOHOST)(long* Mesh_pointer,
                                             realw* h_hess_ac_kl,
                                             int* NSPEC_AB) {
  TRACE("transfer_kernels_hess_ac_tohost");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(h_hess_ac_kl,mp->d_hess_ac_kl,NGLL3*(*NSPEC_AB)*sizeof(realw),
                                     cudaMemcpyDeviceToHost),70202);
}
*/

