/*
!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  8 . 0
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

#ifndef SMOOTH_CUDA_H
#define SMOOTH_CUDA_H

#include "mesh_constants_gpu.h"

// prototype definitions
// from file process_smooth.cu
__global__ void process_smooth( realw_const_p xstore_me,
                                realw_const_p ystore_me,
                                realw_const_p zstore_me,
                                realw_const_p xstore_other,
                                realw_const_p ystore_other,
                                realw_const_p zstore_other,
                                realw_const_p data_other,
                                const realw sigma_h2_inv,
                                const realw sigma_v2_inv,
                                const int iker,
                                const int nspec_me,
                                const int nspec_other,
                                const realw v_criterion,
                                const realw h_criterion,
                                realw_const_p integ_factor,
                                realw_p data_smooth,
                                realw_p normalisation);

__global__ void normalize_data( realw_p data_smooth,
                                realw_const_p normalisation,
                                int nker,
                                int nspec_me);



// smoothing data structure
typedef struct Smooth_data_ {

  gpu_realw_mem x_me;
  gpu_realw_mem y_me;
  gpu_realw_mem z_me;

  gpu_realw_mem x_other;
  gpu_realw_mem y_other;
  gpu_realw_mem z_other;

  gpu_realw_mem integ_factor;
  gpu_realw_mem kernel;

  gpu_realw_mem data_smooth;
  gpu_realw_mem normalisation;

  realw sigma_h2_inv;
  realw sigma_v2_inv;

  realw h_criterion;
  realw v_criterion;

  int nspec_me;
  int nker;
} Smooth_data;


#endif  // SMOOTH_CUDA_H
