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

#include "mesh_constants_gpu.h"

/* ----------------------------------------------------------------------------------------------- */
// wavefield transfers
/* ----------------------------------------------------------------------------------------------- */

// crust_mantle
extern EXTERN_LANG
void FC_FUNC_(transfer_fields_cm_to_device,
              TRANSFER_FIELDS_CM_TO_DEVICE)(int *size, realw *displ, realw *veloc, realw *accel, long *Mesh_pointer_f) {

  TRACE("transfer_fields_cm_to_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  gpuCopy_todevice_realw (&mp->d_displ_crust_mantle, displ, *size);
  gpuCopy_todevice_realw (&mp->d_veloc_crust_mantle, veloc, *size);
  gpuCopy_todevice_realw (&mp->d_accel_crust_mantle, accel, *size);

  GPU_ERROR_CHECKING ("after transfer_fields_cm_to_device");
}

/* ----------------------------------------------------------------------------------------------- */

// inner_core
extern EXTERN_LANG
void FC_FUNC_(transfer_fields_ic_to_device,
              TRANSFER_FIELDS_IC_TO_DEVICE)(int *size, realw *displ, realw *veloc, realw *accel, long *Mesh_pointer_f) {

  TRACE("transfer_fields_ic_to_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  gpuCopy_todevice_realw (&mp->d_displ_inner_core, displ, *size);
  gpuCopy_todevice_realw (&mp->d_veloc_inner_core, veloc, *size);
  gpuCopy_todevice_realw (&mp->d_accel_inner_core, accel, *size);

  GPU_ERROR_CHECKING ("after transfer_fields_ic_to_device");
}

/* ----------------------------------------------------------------------------------------------- */

// outer_core
extern EXTERN_LANG
void FC_FUNC_(transfer_fields_oc_to_device,
              TRANSFER_FIELDS_OC_TO_DEVICE)(int *size, realw *displ, realw *veloc, realw *accel, long *Mesh_pointer_f) {

  TRACE("transfer_fields_oc_to_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  gpuCopy_todevice_realw (&mp->d_displ_outer_core, displ, *size);
  gpuCopy_todevice_realw (&mp->d_veloc_outer_core, veloc, *size);
  gpuCopy_todevice_realw (&mp->d_accel_outer_core, accel, *size);

  GPU_ERROR_CHECKING ("after transfer_fields_oc_to_device");
}

/* ----------------------------------------------------------------------------------------------- */
// backward/reconstructed fields

// crust_mantle
extern EXTERN_LANG
void FC_FUNC_(transfer_b_fields_cm_to_device,
              TRANSFER_FIELDS_B_CM_TO_DEVICE)(int *size, realw *b_displ, realw *b_veloc, realw *b_accel, long *Mesh_pointer_f) {

  TRACE("transfer_fields_b_cm_to_device");
  // debug
  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies last wavefield snapshots to GPU
  gpuCopy_todevice_realw (&mp->d_b_displ_crust_mantle, b_displ, *size);
  gpuCopy_todevice_realw (&mp->d_b_veloc_crust_mantle, b_veloc, *size);
  gpuCopy_todevice_realw (&mp->d_b_accel_crust_mantle, b_accel, *size);

  GPU_ERROR_CHECKING ("after transfer_b_fields_cm_to_device");
}

/* ----------------------------------------------------------------------------------------------- */

// inner_core
extern EXTERN_LANG
void FC_FUNC_(transfer_b_fields_ic_to_device,
              TRANSFER_FIELDS_B_IC_TO_DEVICE)(int *size, realw *b_displ, realw *b_veloc, realw *b_accel, long *Mesh_pointer_f) {

  TRACE("transfer_fields_b_ic_to_device");
  // debug
  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies last wavefield snapshots to GPU
  gpuCopy_todevice_realw (&mp->d_b_displ_inner_core, b_displ, *size);
  gpuCopy_todevice_realw (&mp->d_b_veloc_inner_core, b_veloc, *size);
  gpuCopy_todevice_realw (&mp->d_b_accel_inner_core, b_accel, *size);

  GPU_ERROR_CHECKING ("after transfer_b_fields_ic_to_device");
}

/* ----------------------------------------------------------------------------------------------- */

// outer_core
extern EXTERN_LANG
void FC_FUNC_(transfer_b_fields_oc_to_device,
              TRANSFER_FIELDS_B_OC_TO_DEVICE)(int *size, realw *b_displ, realw *b_veloc, realw *b_accel, long *Mesh_pointer_f) {

  TRACE("transfer_fields_b_oc_to_device");
  // debug
  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies last wavefield snapshots to GPU
  gpuCopy_todevice_realw (&mp->d_b_displ_outer_core, b_displ, *size);
  gpuCopy_todevice_realw (&mp->d_b_veloc_outer_core, b_veloc, *size);
  gpuCopy_todevice_realw (&mp->d_b_accel_outer_core, b_accel, *size);

  GPU_ERROR_CHECKING ("after transfer_b_fields_oc_to_device");
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

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_displ_crust_mantle, displ, *size);
  gpuCopy_from_device_realw (&mp->d_veloc_crust_mantle, veloc, *size);
  gpuCopy_from_device_realw (&mp->d_accel_crust_mantle, accel, *size);

  GPU_ERROR_CHECKING ("after transfer_fields_cm_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

// inner_core
extern EXTERN_LANG
void FC_FUNC_(transfer_fields_ic_from_device,
              TRANSFER_FIELDS_IC_FROM_DEVICE)(int *size, realw *displ, realw *veloc, realw *accel, long *Mesh_pointer_f) {

  TRACE("transfer_fields_ic_from_device");
  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_displ_inner_core, displ, *size);
  gpuCopy_from_device_realw (&mp->d_veloc_inner_core, veloc, *size);
  gpuCopy_from_device_realw (&mp->d_accel_inner_core, accel, *size);

  GPU_ERROR_CHECKING ("after transfer_fields_ic_from_device");
}
/* ----------------------------------------------------------------------------------------------- */


// outer_core
extern EXTERN_LANG
void FC_FUNC_(transfer_fields_oc_from_device,
              TRANSFER_FIELDS_OC_FROM_DEVICE)(int *size, realw *displ, realw *veloc, realw *accel, long *Mesh_pointer_f) {

  TRACE("transfer_fields_oc_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_displ_outer_core, displ, *size);
  gpuCopy_from_device_realw (&mp->d_veloc_outer_core, veloc, *size);
  gpuCopy_from_device_realw (&mp->d_accel_outer_core, accel, *size);

  GPU_ERROR_CHECKING ("after transfer_fields_oc_from_device");
}


/* ----------------------------------------------------------------------------------------------- */

// backward/reconstructed fields

// crust_mantle
extern EXTERN_LANG
void FC_FUNC_(transfer_b_fields_cm_from_device,
              TRANSFER_B_FIELDS_CM_FROM_DEVICE)(int *size, realw *b_displ, realw *b_veloc, realw *b_accel, long *Mesh_pointer_f) {

  TRACE("transfer_b_fields_cm_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_b_displ_crust_mantle, b_displ, *size);
  gpuCopy_from_device_realw (&mp->d_b_veloc_crust_mantle, b_veloc, *size);
  gpuCopy_from_device_realw (&mp->d_b_accel_crust_mantle, b_accel, *size);

  GPU_ERROR_CHECKING ("after transfer_b_fields_cm_from_device");
}

/* ----------------------------------------------------------------------------------------------- */


// inner_core
extern EXTERN_LANG
void FC_FUNC_(transfer_b_fields_ic_from_device,
              TRANSFER_B_FIELDS_IC_FROM_DEVICE)(int *size, realw *b_displ, realw *b_veloc, realw *b_accel, long *Mesh_pointer_f) {
  TRACE("transfer_fields_b_ic_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_b_displ_inner_core, b_displ, *size);
  gpuCopy_from_device_realw (&mp->d_b_veloc_inner_core, b_veloc, *size);
  gpuCopy_from_device_realw (&mp->d_b_accel_inner_core, b_accel, *size);

  GPU_ERROR_CHECKING ("after transfer_b_fields_ic_from_device");
}

/* ----------------------------------------------------------------------------------------------- */


// outer_core
extern EXTERN_LANG
void FC_FUNC_(transfer_b_fields_oc_from_device,
              TRANSFER_B_FIELDS_OC_FROM_DEVICE)(int *size, realw *b_displ, realw *b_veloc, realw *b_accel, long *Mesh_pointer_f) {

  TRACE("transfer_b_fields_oc_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_b_displ_outer_core, b_displ, *size);
  gpuCopy_from_device_realw (&mp->d_b_veloc_outer_core, b_veloc, *size);
  gpuCopy_from_device_realw (&mp->d_b_accel_outer_core, b_accel, *size);

  GPU_ERROR_CHECKING ("after transfer_b_fields_oc_from_device");
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

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_displ_crust_mantle, displ, *size);

  GPU_ERROR_CHECKING ("after transfer_displ_cm_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_displ_cm_from_device,
              TRANSFER_B_DISPL_CM_FROM_DEVICE)(int *size, realw *displ, long *Mesh_pointer_f) {

  TRACE("transfer_b_displ_cm_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_b_displ_crust_mantle, displ, *size);

  GPU_ERROR_CHECKING ("after transfer_b_displ_cm_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_ofs_b_displ_cm_from_device,
              TRANSFER_OFS_B_DISPL_CM_FROM_DEVICE)(int *size, int *offset, realw *displ, long *Mesh_pointer_f) {

  TRACE("transfer_ofs_b_displ_cm_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw_offset (&mp->d_b_displ_crust_mantle, displ, *size, *offset);

  GPU_ERROR_CHECKING ("after transfer_ofs_b_displ_cm_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_displ_cm_to_device,
              TRANSFER_B_DISPL_CM_TO_DEVICE)(int *size, realw *displ, long *Mesh_pointer_f) {

  TRACE("transfer_b_displ_cm_to_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_todevice_realw (&mp->d_b_displ_crust_mantle, displ, *size);

  GPU_ERROR_CHECKING ("after transfer_b_displ_cm_to_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_ofs_b_displ_cm_to_device,
              TRANSFER_OFS_B_DISPL_CM_TO_DEVICE)(int *size, int *offset, realw *displ, long *Mesh_pointer_f) {

  TRACE("transfer_ofs_b_displ_cm_to_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;
  // copies array to CPU
  gpuCopy_todevice_realw_offset (&mp->d_b_displ_crust_mantle, displ, *size, *offset);

  GPU_ERROR_CHECKING ("after transfer_ofs_b_displ_cm_to_device");
}

/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(transfer_displ_ic_from_device,
              TRANSFER_DISPL_IC_FROM_DEVICE)(int *size, realw *displ, long *Mesh_pointer_f) {

  TRACE("transfer_displ_ic_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_displ_inner_core, displ, *size);

  GPU_ERROR_CHECKING ("after transfer_displ_ic_from_device");
}

/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(transfer_b_displ_ic_from_device,
              TRANSFER_B_DISPL_IC_FROM_DEVICE)(int *size, realw *displ, long *Mesh_pointer_f) {

  TRACE("transfer_b_displ_ic_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_b_displ_inner_core, displ, *size);

  GPU_ERROR_CHECKING ("after transfer_b_displ_ic_from_device");
}

/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(transfer_ofs_b_displ_ic_from_device,
              TRANSFER_OFS_B_DISPL_IC_FROM_DEVICE)(int *size, int *offset, realw *displ, long *Mesh_pointer_f) {

  TRACE("transfer_ofs_b_displ_ic_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw_offset (&mp->d_b_displ_inner_core, displ, *size, *offset);

  GPU_ERROR_CHECKING ("after transfer_ofs_b_displ_ic_from_device");
}

/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(transfer_b_displ_ic_to_device,
              TRANSFER_B_DISPL_IC_TO_DEVICE)(int *size, realw *displ, long *Mesh_pointer_f) {

  TRACE("transfer_b_displ_ic_to_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_todevice_realw (&mp->d_b_displ_inner_core, displ, *size);

  GPU_ERROR_CHECKING ("after transfer_b_displ_ic_to_device");
}

/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(transfer_ofs_b_displ_ic_to_device,
              TRANSFER_OFS_B_DISPL_IC_TO_DEVICE)(int *size, int *offset, realw *displ, long *Mesh_pointer_f) {

  TRACE("transfer_ofs_b_displ_ic_to_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_todevice_realw_offset (&mp->d_b_displ_inner_core, displ, *size, *offset);

  GPU_ERROR_CHECKING ("after transfer_ofs_b_displ_ic_to_device");
}
/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(transfer_displ_oc_from_device,
              TRANSFER_DISPL_OC_FROM_DEVICE)(int *size, realw *displ, long *Mesh_pointer_f) {

  TRACE("transfer_displ_oc_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_displ_outer_core, displ, *size);

  GPU_ERROR_CHECKING ("after transfer_displ_oc_from_device");
}


/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_displ_oc_from_device,
              TRANSFER_B_DISPL_OC_FROM_DEVICE)(int *size, realw *b_displ, long *Mesh_pointer_f) {

  TRACE("transfer_b_displ_oc_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_b_displ_outer_core, b_displ, *size);

  GPU_ERROR_CHECKING ("after transfer_b_displ_oc_from_device");
}

/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(transfer_ofs_b_displ_oc_from_device,
              TRANSFER_OFS_B_DISPL_OC_FROM_DEVICE)(int *size, int *offset, realw *b_displ, long *Mesh_pointer_f) {

  TRACE("transfer_ofs_b_displ_oc_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw_offset (&mp->d_b_displ_outer_core, b_displ, *size, *offset);

  GPU_ERROR_CHECKING ("after transfer_ofs_b_displ_oc_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_displ_oc_to_device,
              TRANSFER_B_DISPL_OC_TO_DEVICE)(int *size, realw *b_displ, long *Mesh_pointer_f) {

  TRACE("transfer_b_displ_oc_to_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_todevice_realw (&mp->d_b_displ_outer_core, b_displ, *size);

  GPU_ERROR_CHECKING ("after transfer_b_displ_oc_to_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_ofs_b_displ_oc_to_device,
              TRANSFER_OFS_B_DISPL_OC_TO_DEVICE)(int *size, int *offset, realw *b_displ, long *Mesh_pointer_f) {

  TRACE("transfer_ofs_b_displ_oc_to_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_todevice_realw_offset (&mp->d_b_displ_outer_core, b_displ, *size, *offset);

  GPU_ERROR_CHECKING ("after transfer_ofs_b_displ_oc_to_device");
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

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_veloc_crust_mantle, veloc, *size);

  GPU_ERROR_CHECKING ("after transfer_veloc_cm_from_device");
}

/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(transfer_veloc_ic_from_device,
              TRANSFER_VELOC_IC_FROM_DEVICE)(int *size, realw *veloc, long *Mesh_pointer_f) {

  TRACE("transfer_veloc_ic_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_veloc_inner_core, veloc, *size);

  GPU_ERROR_CHECKING ("after transfer_veloc_ic_from_device");
}
/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_(transfer_veloc_oc_from_device,
              TRANSFER_VELOC_OC_FROM_DEVICE)(int *size, realw *veloc, long *Mesh_pointer_f) {

  TRACE("transfer_veloc_oc_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_veloc_outer_core, veloc, *size);

  GPU_ERROR_CHECKING ("after transfer_veloc_oc_from_device");
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

  gpuCopy_todevice_realw (&mp->d_accel_crust_mantle, accel, *size);

  GPU_ERROR_CHECKING ("after transfer_accel_cm_to_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_accel_cm_from_device,
              TRANSFER_ACCEL_CM_FROM_DEVICE)(int *size, realw *accel, long *Mesh_pointer_f) {

  TRACE("transfer_accel_cm_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_accel_crust_mantle, accel, *size);

  GPU_ERROR_CHECKING ("after transfer_accel_cm_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_accel_cm_from_device,
              TRANSFER_B_ACCEL_CM_FROM_DEVICE)(int *size, realw *b_accel, long *Mesh_pointer_f) {

  TRACE("transfer_b_accel_cm_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_b_accel_crust_mantle, b_accel, *size);

  GPU_ERROR_CHECKING ("after transfer_b_accel_cm_from_device");
}


/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_accel_ic_from_device,
              TRANSFER_ACCEL_IC_FROM_DEVICE)(int *size, realw *accel, long *Mesh_pointer_f) {

  TRACE("transfer_accel_ic_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_accel_inner_core, accel, *size);

  GPU_ERROR_CHECKING ("after transfer_accel_ic_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_accel_oc_from_device,
              TRANSFER_ACCEL_OC_FROM_DEVICE)(int *size, realw *accel, long *Mesh_pointer_f) {

  TRACE("transfer_accel_oc_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_accel_outer_core, accel, *size);

  GPU_ERROR_CHECKING ("after transfer_accel_oc_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_accel_oc_from_device,
              TRANSFER_B_ACCEL_OC_FROM_DEVICE)(int *size, realw *b_accel, long *Mesh_pointer_f) {

  TRACE("transfer_b_accel_oc_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_b_accel_outer_core, b_accel, *size);

  GPU_ERROR_CHECKING ("after transfer_b_accel_oc_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_ofs_b_accel_oc_from_device,
              TRANSFER_OFS_B_ACCEL_OC_FROM_DEVICE)(int *size, int *offset, realw *b_accel, long *Mesh_pointer_f) {

  TRACE("transfer_ofs_b_accel_oc_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw_offset (&mp->d_b_accel_outer_core, b_accel, *size, *offset);

  GPU_ERROR_CHECKING ("after transfer_ofs_b_accel_oc_from_device");
}
/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_b_accel_oc_to_device,
              TRANSFER_B_ACCEL_OC_TO_DEVICE)(int *size, realw *b_accel, long *Mesh_pointer_f) {

  TRACE("transfer_b_accel_oc_to_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_todevice_realw (&mp->d_b_accel_outer_core, b_accel, *size);

  GPU_ERROR_CHECKING ("after transfer_b_accel_oc_to_device");
}

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_ofs_b_accel_oc_to_device,
              TRANSFER_OFS_B_ACCEL_OC_TO_DEVICE)(int *size, int *offset, realw *b_accel, long *Mesh_pointer_f) {

  TRACE("transfer_ofs_b_accel_oc_to_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_todevice_realw_offset (&mp->d_b_accel_outer_core, b_accel, *size, *offset);

  GPU_ERROR_CHECKING ("after transfer_ofs_b_accel_oc_to_device");
}

/* ----------------------------------------------------------------------------------------------- */
// strain fields
/* ----------------------------------------------------------------------------------------------- */

// crust/mantle
extern EXTERN_LANG
void FC_FUNC_(transfer_strain_cm_from_device,
              TRANSFER_STRAIN_CM_FROM_DEVICE)(long *Mesh_pointer_f,
                                              realw *eps_trace_over_3,
                                              realw *epsilondev_xx,
                                              realw *epsilondev_yy,
                                              realw *epsilondev_xy,
                                              realw *epsilondev_xz,
                                              realw *epsilondev_yz) {
  TRACE("transfer_strain_cm_from_device");
  //get mesh pointer out of Fortran integer container

  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int size = NGLL3*mp->NSPEC_CRUST_MANTLE;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_eps_trace_over_3_crust_mantle, eps_trace_over_3, size);
  gpuCopy_from_device_realw (&mp->d_epsilondev_xx_crust_mantle, epsilondev_xx, size);
  gpuCopy_from_device_realw (&mp->d_epsilondev_yy_crust_mantle, epsilondev_yy, size);
  gpuCopy_from_device_realw (&mp->d_epsilondev_xy_crust_mantle, epsilondev_xy, size);
  gpuCopy_from_device_realw (&mp->d_epsilondev_xz_crust_mantle, epsilondev_xz, size);
  gpuCopy_from_device_realw (&mp->d_epsilondev_yz_crust_mantle, epsilondev_yz, size);

  GPU_ERROR_CHECKING ("after transfer_strain_cm_from_device");
}

/* ----------------------------------------------------------------------------------------------- */
// crust/mantle

extern EXTERN_LANG
void FC_FUNC_(transfer_b_strain_cm_to_device,
              TRANSFER_B_STRAIN_CM_TO_DEVICE)(long *Mesh_pointer_f,
                                              realw *epsilondev_xx,
                                              realw *epsilondev_yy,
                                              realw *epsilondev_xy,
                                              realw *epsilondev_xz,
                                              realw *epsilondev_yz) {
  TRACE("transfer_b_strain_cm_to_device");

  // debug
  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int size = NGLL3*mp->NSPEC_CRUST_MANTLE;

  // checks flag
  if (mp->undo_attenuation) exit_on_error("Error invalid undo_attenuation flag in transfer_b_strain_cm_to_device()");

  // transfers strains
  gpuCopy_todevice_realw (&mp->d_b_epsilondev_xx_crust_mantle, epsilondev_xx, size);
  gpuCopy_todevice_realw (&mp->d_b_epsilondev_yy_crust_mantle, epsilondev_yy, size);
  gpuCopy_todevice_realw (&mp->d_b_epsilondev_xy_crust_mantle, epsilondev_xy, size);
  gpuCopy_todevice_realw (&mp->d_b_epsilondev_xz_crust_mantle, epsilondev_xz, size);
  gpuCopy_todevice_realw (&mp->d_b_epsilondev_yz_crust_mantle, epsilondev_yz, size);

  GPU_ERROR_CHECKING ("after transfer_b_strain_cm_to_device");
}

/* ----------------------------------------------------------------------------------------------- */
// inner core

extern EXTERN_LANG
void FC_FUNC_(transfer_strain_ic_from_device,
              TRANSFER_STRAIN_IC_FROM_DEVICE)(long *Mesh_pointer_f,
                                              realw *eps_trace_over_3,
                                              realw *epsilondev_xx,
                                              realw *epsilondev_yy,
                                              realw *epsilondev_xy,
                                              realw *epsilondev_xz,
                                              realw *epsilondev_yz) {
  TRACE("transfer_strain_ic_from_device");
  //get mesh pointer out of Fortran integer container

  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int size = NGLL3 * mp->NSPEC_INNER_CORE;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_eps_trace_over_3_inner_core, eps_trace_over_3, size);
  gpuCopy_from_device_realw (&mp->d_epsilondev_xx_inner_core, epsilondev_xx, size);
  gpuCopy_from_device_realw (&mp->d_epsilondev_yy_inner_core, epsilondev_yy, size);
  gpuCopy_from_device_realw (&mp->d_epsilondev_xy_inner_core, epsilondev_xy, size);
  gpuCopy_from_device_realw (&mp->d_epsilondev_xz_inner_core, epsilondev_xz, size);
  gpuCopy_from_device_realw (&mp->d_epsilondev_yz_inner_core, epsilondev_yz, size);

  GPU_ERROR_CHECKING ("after transfer_strain_ic_from_device");
}

/* ----------------------------------------------------------------------------------------------- */
// inner_core

extern EXTERN_LANG
void FC_FUNC_(transfer_b_strain_ic_to_device,
              TRANSFER_B_STRAIN_IC_TO_DEVICE)(long *Mesh_pointer_f,
                                              realw *epsilondev_xx,
                                              realw *epsilondev_yy,
                                              realw *epsilondev_xy,
                                              realw *epsilondev_xz,
                                              realw *epsilondev_yz) {
  TRACE("transfer_b_strain_cm_to_device");
  // debug
  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int size = NGLL3*mp->NSPEC_INNER_CORE;

  // checks flag
  if (mp->undo_attenuation) exit_on_error("Error invalid undo_attenuation flag in transfer_b_strain_ic_to_device()");

  // transfers strains to GPU
  gpuCopy_todevice_realw (&mp->d_b_epsilondev_xx_inner_core, epsilondev_xx, size);
  gpuCopy_todevice_realw (&mp->d_b_epsilondev_yy_inner_core, epsilondev_yy, size);
  gpuCopy_todevice_realw (&mp->d_b_epsilondev_xy_inner_core, epsilondev_xy, size);
  gpuCopy_todevice_realw (&mp->d_b_epsilondev_xz_inner_core, epsilondev_xz, size);
  gpuCopy_todevice_realw (&mp->d_b_epsilondev_yz_inner_core, epsilondev_yz, size);

  GPU_ERROR_CHECKING ("after transfer_b_strain_ic_to_device");
}

/* ----------------------------------------------------------------------------------------------- */
// R memory variables
/* ----------------------------------------------------------------------------------------------- */

// crust/mantle

extern EXTERN_LANG
void FC_FUNC_(transfer_rmemory_cm_from_device,
              TRANSFER_RMEMORY_CM_FROM_DEVICE)(long* Mesh_pointer_f,
                                               realw* R_xx,
                                               realw* R_yy,
                                               realw* R_xy,
                                               realw* R_xz,
                                               realw* R_yz) {
  TRACE("transfer_rmemory_cm_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int size = N_SLS * NGLL3 * mp->NSPEC_CRUST_MANTLE;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_R_xx_crust_mantle, R_xx, size);
  gpuCopy_from_device_realw (&mp->d_R_yy_crust_mantle, R_yy, size);
  gpuCopy_from_device_realw (&mp->d_R_xy_crust_mantle, R_xy, size);
  gpuCopy_from_device_realw (&mp->d_R_xz_crust_mantle, R_xz, size);
  gpuCopy_from_device_realw (&mp->d_R_yz_crust_mantle, R_yz, size);

  GPU_ERROR_CHECKING ("after transfer_rmemory_cm_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

// crust/mantle

extern EXTERN_LANG
void FC_FUNC_(transfer_b_rmemory_cm_to_device,
              TRANSFER_B_RMEMORY_CM_TO_DEVICE)(long *Mesh_pointer_f,
                                               realw *b_R_xx,
                                               realw *b_R_yy,
                                               realw *b_R_xy,
                                               realw *b_R_xz,
                                               realw *b_R_yz) {
  TRACE("transfer_b_Rmemory_cm_to_device");

  // debug
  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int size = N_SLS*NGLL3*mp->NSPEC_CRUST_MANTLE;

  // checks if anything to do
  if (mp->partial_phys_dispersion_only) return;

  gpuCopy_todevice_realw (&mp->d_b_R_xx_crust_mantle, b_R_xx, size);
  gpuCopy_todevice_realw (&mp->d_b_R_yy_crust_mantle, b_R_yy, size);
  gpuCopy_todevice_realw (&mp->d_b_R_xy_crust_mantle, b_R_xy, size);
  gpuCopy_todevice_realw (&mp->d_b_R_xz_crust_mantle, b_R_xz, size);
  gpuCopy_todevice_realw (&mp->d_b_R_yz_crust_mantle, b_R_yz, size);

  GPU_ERROR_CHECKING ("after transfer_b_rmemory_cm_to_device");
}

/* ----------------------------------------------------------------------------------------------- */
// inner core

extern EXTERN_LANG
void FC_FUNC_(transfer_rmemory_ic_from_device,
              TRANSFER_RMEMORY_IC_FROM_DEVICE)(long* Mesh_pointer_f,
                                               realw* R_xx,
                                               realw* R_yy,
                                               realw* R_xy,
                                               realw* R_xz,
                                               realw* R_yz) {
  TRACE("transfer_rmemory_cm_from_device");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int size = N_SLS*NGLL3*mp->NSPEC_INNER_CORE;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_R_xx_inner_core, R_xx, size);
  gpuCopy_from_device_realw (&mp->d_R_yy_inner_core, R_yy, size);
  gpuCopy_from_device_realw (&mp->d_R_xy_inner_core, R_xy, size);
  gpuCopy_from_device_realw (&mp->d_R_xz_inner_core, R_xz, size);
  gpuCopy_from_device_realw (&mp->d_R_yz_inner_core, R_yz, size);

  GPU_ERROR_CHECKING ("after transfer_rmemory_ic_from_device");
}

/* ----------------------------------------------------------------------------------------------- */
// inner core

extern EXTERN_LANG
void FC_FUNC_(transfer_b_rmemory_ic_to_device,
              TRANSFER_B_RMEMORY_IC_TO_DEVICE)(long *Mesh_pointer_f,
                                               realw *b_R_xx,
                                               realw *b_R_yy,
                                               realw *b_R_xy,
                                               realw *b_R_xz,
                                               realw *b_R_yz) {
  TRACE("transfer_b_rmemory_ic_to_device");
  // debug

  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of Fortran integer container

  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int size = N_SLS*NGLL3*mp->NSPEC_INNER_CORE;

  // checks if anything to do
  if (mp->partial_phys_dispersion_only) return;

  gpuCopy_todevice_realw (&mp->d_b_R_xx_inner_core, b_R_xx, size);
  gpuCopy_todevice_realw (&mp->d_b_R_yy_inner_core, b_R_yy, size);
  gpuCopy_todevice_realw (&mp->d_b_R_xy_inner_core, b_R_xy, size);
  gpuCopy_todevice_realw (&mp->d_b_R_xz_inner_core, b_R_xz, size);
  gpuCopy_todevice_realw (&mp->d_b_R_yz_inner_core, b_R_yz, size);

  GPU_ERROR_CHECKING ("after transfer_b_rmemory_ic_to_device");
}

/* ----------------------------------------------------------------------------------------------- */
// rotation arrays
/* ----------------------------------------------------------------------------------------------- */

// for outer core

extern EXTERN_LANG
void FC_FUNC_(transfer_rotation_from_device,
              TRANSFER_ROTATION_FROM_DEVICE)(long *Mesh_pointer_f,
                                             realw *A_array_rotation,
                                             realw *B_array_rotation) {
  TRACE("transfer_rotation_from_device");

  //get mesh pointer out of Fortran integer container

  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int size = NGLL3*mp->NSPEC_OUTER_CORE;

  // copies arrays to CPU
  gpuCopy_from_device_realw (&mp->d_A_array_rotation, A_array_rotation, size);
  gpuCopy_from_device_realw (&mp->d_B_array_rotation, B_array_rotation, size);

  GPU_ERROR_CHECKING ("after transfer_rotation_from_device");
}

/* ----------------------------------------------------------------------------------------------- */
// for outer core

extern EXTERN_LANG
void FC_FUNC_(transfer_b_rotation_to_device,
              TRANSFER_B_ROTATION_TO_DEVICE)(long *Mesh_pointer_f,
                                             realw *A_array_rotation,
                                             realw *B_array_rotation) {
  TRACE("transfer_b_rotation_to_device");
  // debug

  DEBUG_BACKWARD_TRANSFER();

  //get mesh pointer out of Fortran integer container

  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int size = NGLL3 * mp->NSPEC_OUTER_CORE;

  // copies arrays to CPU
  gpuCopy_from_device_realw (&mp->d_b_A_array_rotation, A_array_rotation, size);
  gpuCopy_from_device_realw (&mp->d_b_B_array_rotation, B_array_rotation, size);

  GPU_ERROR_CHECKING ("after transfer_b_rotation_to_device");
}

/* ----------------------------------------------------------------------------------------------- */
// KERNEL transfers
/* ----------------------------------------------------------------------------------------------- */

// crust/mantle

extern EXTERN_LANG
void FC_FUNC_(transfer_kernels_cm_to_host,
              TRANSFER_KERNELS_CM_TO_HOST)(long *Mesh_pointer_f,
                                           realw *h_rho_kl,
                                           realw *h_alpha_kl,
                                           realw *h_beta_kl,
                                           int *NSPEC) {
  TRACE("transfer_kernels_cm_to_host");

  //get mesh pointer out of Fortran integer container

  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int size = (*NSPEC)*NGLL3;


  // density kernel
  // copies arrays to CPU
  gpuCopy_from_device_realw (&mp->d_rho_kl_crust_mantle, h_rho_kl, size);

  // isotropic kernels
  if (! mp->anisotropic_kl) {
    gpuCopy_from_device_realw (&mp->d_alpha_kl_crust_mantle, h_alpha_kl, size);
    gpuCopy_from_device_realw (&mp->d_beta_kl_crust_mantle, h_beta_kl, size);
  }

  GPU_ERROR_CHECKING ("after transfer_kernels_cm_to_host");
}


/* ----------------------------------------------------------------------------------------------- */

// crust/mantle

extern EXTERN_LANG
void FC_FUNC_(transfer_kernels_ani_cm_to_host,
              TRANSFER_KERNELS_ANI_CM_TO_HOST)(long *Mesh_pointer_f,
                                               realw *h_cijkl_kl,
                                               int *NSPEC) {
  TRACE("transfer_kernels_ani_cm_to_host");

  //get mesh pointer out of Fortran integer container

  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int size = 21 * (*NSPEC) * NGLL3;

  // anisotropic kernel
  if (! mp->anisotropic_kl) { exit_on_error("Error invalid ANISOTROPIC_KL in transfer_kernels_ani_cm_to_host() routine");}

  // anisotropic kernels
  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_cijkl_kl_crust_mantle, h_cijkl_kl, size);

  GPU_ERROR_CHECKING ("after transfer_kernels_ani_cm_to_host");
}

/* ----------------------------------------------------------------------------------------------- */
// inner core

extern EXTERN_LANG
void FC_FUNC_(transfer_kernels_ic_to_host,
              TRANSFER_KERNELS_IC_TO_HOST)(long *Mesh_pointer_f,
                                           realw *h_rho_kl,
                                           realw *h_alpha_kl,
                                           realw *h_beta_kl,
                                           int *NSPEC) {
  TRACE("transfer_kernels_ic_to_host");

  //get mesh pointer out of Fortran integer container

  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int size = (*NSPEC) * NGLL3;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_rho_kl_inner_core, h_rho_kl, size);
  gpuCopy_from_device_realw (&mp->d_alpha_kl_inner_core, h_alpha_kl, size);
  gpuCopy_from_device_realw (&mp->d_beta_kl_inner_core, h_beta_kl, size);

  GPU_ERROR_CHECKING ("after transfer_kernels_ic_to_host");
}

/* ----------------------------------------------------------------------------------------------- */
// outer core

extern EXTERN_LANG
void FC_FUNC_(transfer_kernels_oc_to_host,
              TRANSFER_KERNELS_OC_TO_HOST)(long *Mesh_pointer_f,
                                           realw *h_rho_kl,
                                           realw *h_alpha_kl,
                                           int *NSPEC) {

  TRACE("transfer_kernels_oc_to_host");

  //get mesh pointer out of Fortran integer container

  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  int size = (*NSPEC)*NGLL3;

  // copies kernel values over to CPU host
  gpuCopy_from_device_realw (&mp->d_rho_kl_outer_core, h_rho_kl, size);
  gpuCopy_from_device_realw (&mp->d_alpha_kl_outer_core, h_alpha_kl, size);

  GPU_ERROR_CHECKING ("after transfer_kernels_oc_to_host");
}

/* ----------------------------------------------------------------------------------------------- */
// for NOISE simulations

extern EXTERN_LANG
void FC_FUNC_(transfer_kernels_noise_to_host,
              TRANSFER_KERNELS_NOISE_TO_HOST)(long *Mesh_pointer_f,
                                              realw *h_Sigma_kl,
                                              int *NSPEC) {
  TRACE("transfer_kernels_noise_to_host");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_Sigma_kl, h_Sigma_kl, NGLL3 * (*NSPEC));

  GPU_ERROR_CHECKING ("after transfer_kernels_noise_to_host");
}

/* ----------------------------------------------------------------------------------------------- */
// for Hess kernel calculations

extern EXTERN_LANG
void FC_FUNC_(transfer_kernels_hess_cm_tohost,
              TRANSFER_KERNELS_HESS_CM_TOHOST)(long *Mesh_pointer_f,
                                               realw *h_hess_kl,
                                               int *NSPEC) {
  TRACE("transfer_kernels_hess_cm_tohost");

  //get mesh pointer out of Fortran integer container
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // copies array to CPU
  gpuCopy_from_device_realw (&mp->d_hess_kl_crust_mantle, h_hess_kl, NGLL3 * (*NSPEC));

  GPU_ERROR_CHECKING ("after transfer_kernels_hess_cm_tohost");
}

/* ----------------------------------------------------------------------------------------------- */
// register host array for pinned memory

extern EXTERN_LANG
void FC_FUNC_(register_host_array,
              REGISTER_HOST_ARRAY)(int *size, realw *h_array) {

  TRACE("register_host_array");

  gpuRegisterHost_realw ( h_array, *size);

  GPU_ERROR_CHECKING ("after register_host_array");
}


extern EXTERN_LANG
void FC_FUNC_(unregister_host_array,
              UNREGISTER_HOST_ARRAY)(realw *h_array) {

  TRACE("unregister_host_array");

  gpuUnregisterHost_realw (h_array);

  GPU_ERROR_CHECKING ("after unregister_host_array");
}
