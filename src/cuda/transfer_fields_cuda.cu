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

#include <stdio.h>
#include <cuda.h>
#include <cublas.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"

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

TRACE("transfer_fields_cm_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(mp->d_displ_crust_mantle,displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40003);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_veloc_crust_mantle,veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40004);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_accel_crust_mantle,accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40005);

}

// inner_core
extern "C"
void FC_FUNC_(transfer_fields_ic_to_device,
              TRANSFER_FIELDS_IC_TO_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer_f) {

  TRACE("transfer_fields_ic_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(mp->d_displ_inner_core,displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40003);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_veloc_inner_core,veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40004);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_accel_inner_core,accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40005);

}

// outer_core
extern "C"
void FC_FUNC_(transfer_fields_oc_to_device,
              TRANSFER_FIELDS_OC_TO_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer_f) {

  TRACE("transfer_fields_oc_to_device");

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

  TRACE("transfer_fields_b_cm_to_device");
  // debug
  DEBUG_BACKWARD_TRANSFER();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_displ_crust_mantle,b_displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40003);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_veloc_crust_mantle,b_veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40004);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_accel_crust_mantle,b_accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40005);

}

// inner_core
extern "C"
void FC_FUNC_(transfer_b_fields_ic_to_device,
              TRANSFER_FIELDS_B_IC_TO_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                              long* Mesh_pointer_f) {

  TRACE("transfer_fields_b_ic_to_device");
  // debug
  DEBUG_BACKWARD_TRANSFER();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_displ_inner_core,b_displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40003);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_veloc_inner_core,b_veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40004);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_accel_inner_core,b_accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40005);

}

// outer_core
extern "C"
void FC_FUNC_(transfer_b_fields_oc_to_device,
              TRANSFER_FIELDS_B_OC_TO_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                              long* Mesh_pointer_f) {

  TRACE("transfer_fields_b_oc_to_device");
  // debug
  DEBUG_BACKWARD_TRANSFER();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_displ_outer_core,b_displ,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40003);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_veloc_outer_core,b_veloc,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40004);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_accel_outer_core,b_accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40005);

}

/* ----------------------------------------------------------------------------------------------- */

// transfer memory from GPU device to CPU host

/* ----------------------------------------------------------------------------------------------- */

// crust_mantle
extern "C"
void FC_FUNC_(transfer_fields_cm_from_device,
              TRANSFER_FIELDS_CM_FROM_DEVICE)(int* size, realw* displ, realw* veloc, realw* accel,long* Mesh_pointer_f) {

  TRACE("transfer_fields_cm_from_device");

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

  TRACE("transfer_fields_ic_from_device");

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

  TRACE("transfer_fields_oc_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_displ_outer_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
  print_CUDA_error_if_any(cudaMemcpy(veloc,mp->d_veloc_outer_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
  print_CUDA_error_if_any(cudaMemcpy(accel,mp->d_accel_outer_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40008);

}

/* ----------------------------------------------------------------------------------------------- */

// backward/reconstructed fields

// crust_mantle
extern "C"
void FC_FUNC_(transfer_b_fields_cm_from_device,
              TRANSFER_B_FIELDS_CM_FROM_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                                long* Mesh_pointer_f) {

TRACE("transfer_b_fields_cm_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(b_displ,mp->d_b_displ_crust_mantle,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
  print_CUDA_error_if_any(cudaMemcpy(b_veloc,mp->d_b_veloc_crust_mantle,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
  print_CUDA_error_if_any(cudaMemcpy(b_accel,mp->d_b_accel_crust_mantle,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40008);

}

/* ----------------------------------------------------------------------------------------------- */

// inner_core
extern "C"
void FC_FUNC_(transfer_b_fields_ic_from_device,
              TRANSFER_B_FIELDS_IC_FROM_DEVICE)(int* size, realw* b_displ, realw* b_veloc, realw* b_accel,
                                                long* Mesh_pointer_f) {
  TRACE("transfer_fields_b_ic_from_device");

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

  TRACE("transfer_b_fields_oc_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(b_displ,mp->d_b_displ_outer_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
  print_CUDA_error_if_any(cudaMemcpy(b_veloc,mp->d_b_veloc_outer_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
  print_CUDA_error_if_any(cudaMemcpy(b_accel,mp->d_b_accel_outer_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40008);

}

/* ----------------------------------------------------------------------------------------------- */

// single wavefield transfers

/* ----------------------------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------------------------- */

// displacements

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
void FC_FUNC_(transfer_displ_oc_from_device,
              TRANSFER_DISPL_OC_FROM_DEVICE)(int* size, realw* displ, long* Mesh_pointer_f) {

  TRACE("transfer_displ_oc_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_displ_outer_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_b_displ_oc_from_device,
              TRANSFER_B_DISPL_OC_FROM_DEVICE)(int* size, realw* displ, long* Mesh_pointer_f) {

  TRACE("transfer_b_displ_oc_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(displ,mp->d_b_displ_outer_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40006);
}


/* ----------------------------------------------------------------------------------------------- */

// velocities

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_veloc_cm_from_device,
              TRANSFER_VELOC_CM_FROM_DEVICE)(int* size, realw* veloc, long* Mesh_pointer_f) {

  TRACE("transfer_veloc_cm_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(veloc,mp->d_veloc_crust_mantle,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_veloc_ic_from_device,
              TRANSFER_VELOC_IC_FROM_DEVICE)(int* size, realw* veloc, long* Mesh_pointer_f) {

  TRACE("transfer_veloc_ic_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(veloc,mp->d_veloc_inner_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_veloc_oc_from_device,
              TRANSFER_VELOC_OC_FROM_DEVICE)(int* size, realw* veloc, long* Mesh_pointer_f) {

  TRACE("transfer_veloc_oc_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(veloc,mp->d_veloc_outer_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40007);
}


/* ----------------------------------------------------------------------------------------------- */

// accelerations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_accel_cm_to_device,
              TRANSFER_ACCEL_CM_TO_DEVICE)(int* size, realw* accel,long* Mesh_pointer_f) {

TRACE("transfer_accel_cm_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(mp->d_accel_crust_mantle,accel,sizeof(realw)*(*size),cudaMemcpyHostToDevice),40016);

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
              TRANSFER_B_ACCEL_CM_FROM_DEVICE)(int* size, realw* b_accel,long* Mesh_pointer_f) {

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

extern "C"
void FC_FUNC_(transfer_accel_oc_from_device,
              TRANSFER_ACCEL_OC_FROM_DEVICE)(int* size, realw* accel,long* Mesh_pointer_f) {

  TRACE("transfer_accel_oc_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(accel,mp->d_accel_outer_core,sizeof(realw)*(*size),cudaMemcpyDeviceToHost),40026);

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

  print_CUDA_error_if_any(cudaMemcpy(eps_trace_over_3,mp->d_eps_trace_over_3_crust_mantle,size*sizeof(realw),cudaMemcpyDeviceToHost),320001);

  print_CUDA_error_if_any(cudaMemcpy(epsilondev_xx,mp->d_epsilondev_xx_crust_mantle,size*sizeof(realw),cudaMemcpyDeviceToHost),320002);
  print_CUDA_error_if_any(cudaMemcpy(epsilondev_yy,mp->d_epsilondev_yy_crust_mantle,size*sizeof(realw),cudaMemcpyDeviceToHost),320003);
  print_CUDA_error_if_any(cudaMemcpy(epsilondev_xy,mp->d_epsilondev_xy_crust_mantle,size*sizeof(realw),cudaMemcpyDeviceToHost),320004);
  print_CUDA_error_if_any(cudaMemcpy(epsilondev_xz,mp->d_epsilondev_xz_crust_mantle,size*sizeof(realw),cudaMemcpyDeviceToHost),320005);
  print_CUDA_error_if_any(cudaMemcpy(epsilondev_yz,mp->d_epsilondev_yz_crust_mantle,size*sizeof(realw),cudaMemcpyDeviceToHost),320006);


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
  // debug
  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = NGLL3*mp->NSPEC_CRUST_MANTLE;

  if( ! mp->undo_attenuation ){
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xx_crust_mantle,epsilondev_xx,size*sizeof(realw),cudaMemcpyHostToDevice),330001);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_yy_crust_mantle,epsilondev_yy,size*sizeof(realw),cudaMemcpyHostToDevice),330002);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xy_crust_mantle,epsilondev_xy,size*sizeof(realw),cudaMemcpyHostToDevice),330003);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xz_crust_mantle,epsilondev_xz,size*sizeof(realw),cudaMemcpyHostToDevice),330004);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_yz_crust_mantle,epsilondev_yz,size*sizeof(realw),cudaMemcpyHostToDevice),330005);
  }

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

  print_CUDA_error_if_any(cudaMemcpy(eps_trace_over_3,mp->d_eps_trace_over_3_inner_core,size*sizeof(realw),cudaMemcpyDeviceToHost),340001);

  print_CUDA_error_if_any(cudaMemcpy(epsilondev_xx,mp->d_epsilondev_xx_inner_core,size*sizeof(realw),cudaMemcpyDeviceToHost),340002);
  print_CUDA_error_if_any(cudaMemcpy(epsilondev_yy,mp->d_epsilondev_yy_inner_core,size*sizeof(realw),cudaMemcpyDeviceToHost),340003);
  print_CUDA_error_if_any(cudaMemcpy(epsilondev_xy,mp->d_epsilondev_xy_inner_core,size*sizeof(realw),cudaMemcpyDeviceToHost),340004);
  print_CUDA_error_if_any(cudaMemcpy(epsilondev_xz,mp->d_epsilondev_xz_inner_core,size*sizeof(realw),cudaMemcpyDeviceToHost),340005);
  print_CUDA_error_if_any(cudaMemcpy(epsilondev_yz,mp->d_epsilondev_yz_inner_core,size*sizeof(realw),cudaMemcpyDeviceToHost),340006);


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
  // debug
  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = NGLL3*mp->NSPEC_INNER_CORE;

  if( ! mp->undo_attenuation ){
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xx_inner_core,epsilondev_xx,size*sizeof(realw),cudaMemcpyHostToDevice),350001);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_yy_inner_core,epsilondev_yy,size*sizeof(realw),cudaMemcpyHostToDevice),350002);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xy_inner_core,epsilondev_xy,size*sizeof(realw),cudaMemcpyHostToDevice),350003);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xz_inner_core,epsilondev_xz,size*sizeof(realw),cudaMemcpyHostToDevice),350004);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_yz_inner_core,epsilondev_yz,size*sizeof(realw),cudaMemcpyHostToDevice),350005);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_b_strain_ic_to_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// R memory variables

/* ----------------------------------------------------------------------------------------------- */

// crust/mantle

extern "C"
void FC_FUNC_(transfer_b_rmemory_cm_to_device,
              TRANSFER_B_RMEMORY_CM_TO_DEVICE)(long* Mesh_pointer,
                                               realw* b_R_xx,
                                               realw* b_R_yy,
                                               realw* b_R_xy,
                                               realw* b_R_xz,
                                               realw* b_R_yz) {
  TRACE("transfer_b_Rmemory_cm_to_device");
  // debug
  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = N_SLS*NGLL3*mp->NSPEC_CRUST_MANTLE;

  if( ! mp->partial_phys_dispersion_only ){
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_xx_crust_mantle,b_R_xx,size*sizeof(realw),cudaMemcpyHostToDevice),360001);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_yy_crust_mantle,b_R_yy,size*sizeof(realw),cudaMemcpyHostToDevice),360002);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_xy_crust_mantle,b_R_xy,size*sizeof(realw),cudaMemcpyHostToDevice),360003);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_xz_crust_mantle,b_R_xz,size*sizeof(realw),cudaMemcpyHostToDevice),360004);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_yz_crust_mantle,b_R_yz,size*sizeof(realw),cudaMemcpyHostToDevice),360005);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_b_rmemory_cm_to_device");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// inner core

extern "C"
void FC_FUNC_(transfer_b_rmemory_ic_to_device,
              TRANSFER_B_RMEMORY_IC_TO_DEVICE)(long* Mesh_pointer,
                                               realw* b_R_xx,
                                               realw* b_R_yy,
                                               realw* b_R_xy,
                                               realw* b_R_xz,
                                               realw* b_R_yz) {
  TRACE("transfer_b_rmemory_ic_to_device");
  // debug
  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = N_SLS*NGLL3*mp->NSPEC_INNER_CORE;

  if( ! mp->partial_phys_dispersion_only ){
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_xx_inner_core,b_R_xx,size*sizeof(realw),cudaMemcpyHostToDevice),370001);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_yy_inner_core,b_R_yy,size*sizeof(realw),cudaMemcpyHostToDevice),370002);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_xy_inner_core,b_R_xy,size*sizeof(realw),cudaMemcpyHostToDevice),370003);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_xz_inner_core,b_R_xz,size*sizeof(realw),cudaMemcpyHostToDevice),370004);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_R_yz_inner_core,b_R_yz,size*sizeof(realw),cudaMemcpyHostToDevice),370005);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after transfer_b_rmemory_ic_to_device");
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

  print_CUDA_error_if_any(cudaMemcpy(A_array_rotation,mp->d_A_array_rotation,
                                     size*sizeof(realw),cudaMemcpyDeviceToHost),380001);
  print_CUDA_error_if_any(cudaMemcpy(B_array_rotation,mp->d_B_array_rotation,
                                     size*sizeof(realw),cudaMemcpyDeviceToHost),380002);
}

/* ----------------------------------------------------------------------------------------------- */

// for outer core

extern "C"
void FC_FUNC_(transfer_b_rotation_to_device,
              TRANSFER_B_ROTATION_TO_DEVICE)(long* Mesh_pointer,
                                              realw* A_array_rotation,
                                              realw* B_array_rotation) {
  TRACE("transfer_b_rotation_to_device");
  // debug
  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = NGLL3*mp->NSPEC_OUTER_CORE;

  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_A_array_rotation,A_array_rotation,
                                     size*sizeof(realw),cudaMemcpyHostToDevice),390001);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_b_B_array_rotation,B_array_rotation,
                                     size*sizeof(realw),cudaMemcpyHostToDevice),390002);
}


/* ----------------------------------------------------------------------------------------------- */

// KERNEL transfers

/* ----------------------------------------------------------------------------------------------- */

// crust/mantle

extern "C"
void FC_FUNC_(transfer_kernels_cm_to_host,
              TRANSFER_KERNELS_CM_TO_HOST)(long* Mesh_pointer,
                                           realw* h_rho_kl,
                                           realw* h_alpha_kl,
                                           realw* h_beta_kl,
                                           realw* h_cijkl_kl,
                                           int* NSPEC) {
  TRACE("transfer_kernels_cm_to_host");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = (*NSPEC)*NGLL3;

  // density kernel
  print_CUDA_error_if_any(cudaMemcpy(h_rho_kl,mp->d_rho_kl_crust_mantle,
                                     size*sizeof(realw),cudaMemcpyDeviceToHost),40101);

  if( ! mp->anisotropic_kl){
    // isotropic kernels
    print_CUDA_error_if_any(cudaMemcpy(h_alpha_kl,mp->d_alpha_kl_crust_mantle,
                                       size*sizeof(realw),cudaMemcpyDeviceToHost),40102);
    print_CUDA_error_if_any(cudaMemcpy(h_beta_kl,mp->d_beta_kl_crust_mantle,
                                       size*sizeof(realw),cudaMemcpyDeviceToHost),40103);
  }else{
    // anisotropic kernels
    print_CUDA_error_if_any(cudaMemcpy(h_cijkl_kl,mp->d_cijkl_kl_crust_mantle,
                                       21*size*sizeof(realw),cudaMemcpyDeviceToHost),40102);
  }
}

/* ----------------------------------------------------------------------------------------------- */

// inner core

extern "C"
void FC_FUNC_(transfer_kernels_ic_to_host,
              TRANSFER_KERNELS_IC_TO_HOST)(long* Mesh_pointer,
                                           realw* h_rho_kl,
                                           realw* h_alpha_kl,
                                           realw* h_beta_kl,
                                           int* NSPEC) {
TRACE("transfer_kernels_ic_to_host");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = (*NSPEC)*NGLL3;

  print_CUDA_error_if_any(cudaMemcpy(h_rho_kl,mp->d_rho_kl_inner_core,
                                     size*sizeof(realw),cudaMemcpyDeviceToHost),40101);
  print_CUDA_error_if_any(cudaMemcpy(h_alpha_kl,mp->d_alpha_kl_inner_core,
                                     size*sizeof(realw),cudaMemcpyDeviceToHost),40102);
  print_CUDA_error_if_any(cudaMemcpy(h_beta_kl,mp->d_beta_kl_inner_core,
                                     size*sizeof(realw),cudaMemcpyDeviceToHost),40103);
}


/* ----------------------------------------------------------------------------------------------- */

// outer core

extern "C"
void FC_FUNC_(transfer_kernels_oc_to_host,
              TRANSFER_KERNELS_OC_TO_HOST)(long* Mesh_pointer,
                                           realw* h_rho_kl,
                                           realw* h_alpha_kl,
                                           int* NSPEC) {

  TRACE("transfer_kernels_oc_to_host");

  //get mesh pointer out of fortran integer container
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = (*NSPEC)*NGLL3;

  // copies kernel values over to CPU host
  print_CUDA_error_if_any(cudaMemcpy(h_rho_kl,mp->d_rho_kl_outer_core,
                                     size*sizeof(realw),cudaMemcpyDeviceToHost),54101);
  print_CUDA_error_if_any(cudaMemcpy(h_alpha_kl,mp->d_alpha_kl_outer_core,
                                     size*sizeof(realw),cudaMemcpyDeviceToHost),54102);
}

/* ----------------------------------------------------------------------------------------------- */

// for NOISE simulations

extern "C"
void FC_FUNC_(transfer_kernels_noise_to_host,
              TRANSFER_KERNELS_NOISE_TO_HOST)(long* Mesh_pointer,
                                              realw* h_Sigma_kl,
                                              int* NSPEC) {
  TRACE("transfer_kernels_noise_to_host");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  print_CUDA_error_if_any(cudaMemcpy(h_Sigma_kl,mp->d_Sigma_kl,NGLL3*(*NSPEC)*sizeof(realw),
                                     cudaMemcpyDeviceToHost),40201);
}

/* ----------------------------------------------------------------------------------------------- */

// for Hess kernel calculations

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

