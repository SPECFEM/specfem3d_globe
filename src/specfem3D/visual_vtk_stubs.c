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

#include <stdio.h>
#include <stdlib.h>

#include "config.h"

/* ----------------------------------------------------------------------------------------------- */

// empty stubs

/* ----------------------------------------------------------------------------------------------- */

void FC_FUNC_(initialize_vtkwindow,INITIALIZE_VTKWINDOW)(int* GPU_MODE) {
  //printf("VTK_MODE not enabled at compile time: use -DHAVE_VTK\n");
}

void FC_FUNC_(prepare_vtksource,PREPARE_VTKSOURCE)(float* xs_x,float* xs_y, float* xs_z) {}

void FC_FUNC_(prepare_vtkfreesurface,PREPARE_VTKFREESURFACE)(int* free_np,
                                                             float* free_x,
                                                             float* free_y,
                                                             float* free_z,
                                                             int* free_nspec,
                                                             int* free_conn) {}

void FC_FUNC_(prepare_vtkfield,PREPARE_VTKFIELD)(int* vol_np,
                                                 float* vol_x,
                                                 float* vol_y,
                                                 float* vol_z,
                                                 int* vol_nspec,
                                                 int* vol_conn) {}

void FC_FUNC_(visualize_vtkdata,VISUALIZE_VTKDATA)(int* it_h,float* time_h, float* data) {}

void FC_FUNC_(finish_vtkwindow,FINISH_VTKWINDOW)() {}
