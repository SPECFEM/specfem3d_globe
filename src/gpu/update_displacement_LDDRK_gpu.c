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

#include "mesh_constants_gpu.h"


/*----------------------------------------------------------------------------------------------- */
// update LDDRK
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (update_displ_lddrk_gpu,
               UPDATE_DISPL_LDDRK_GPU) (long *Mesh_pointer_f,
                                        int *FORWARD_OR_ADJOINT) {

  TRACE ("update_displ_lddrk_gpu");

  //get Mesh from Fortran integer wrapper
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in update_displ_lddrk_gpu() routine");
  }

  // crust/mantle
  int size_cm = NDIM * mp->NGLOB_CRUST_MANTLE;
  int size_oc = mp->NGLOB_OUTER_CORE;
  int size_ic = NDIM * mp->NGLOB_INNER_CORE;

  //debug
#if (DEBUG_BACKWARD_SIMULATIONS == 1 && DEBUG == 1 && DEBUG_FIELDS == 1)
  {
    realw max_d, max_v, max_a;
    max_d = get_device_array_maximum_value(mp->d_b_displ_crust_mantle, size_cm);
    max_v = get_device_array_maximum_value(mp->d_b_veloc_crust_mantle, size_cm);
    max_a = get_device_array_maximum_value(mp->d_b_accel_crust_mantle, size_cm);
    printf ("rank %d - forward/adjoint: %i, max crust_mantle displ: %f veloc: %f accel: %f\n",
            mp->myrank, *FORWARD_OR_ADJOINT, max_d, max_v, max_a);
    fflush (stdout);
    synchronize_mpi ();
  }
#endif

  // sets gpu arrays
  gpu_realw_mem accel_cm,accel_oc,accel_ic;
  if (*FORWARD_OR_ADJOINT == 1) {
    accel_cm = mp->d_accel_crust_mantle;
    accel_oc = mp->d_accel_outer_core;
    accel_ic = mp->d_accel_inner_core;
  } else {
    //debug
    DEBUG_BACKWARD_UPDATE ();
    // for backward/reconstructed fields
    accel_cm = mp->d_b_accel_crust_mantle;
    accel_oc = mp->d_b_accel_outer_core;
    accel_ic = mp->d_b_accel_inner_core;
  }

  // kernel timing
  gpu_event start,stop;
  if (GPU_TIMING){ start_timing_gpu(&start,&stop); }

  gpuMemset_realw (&accel_cm, size_cm, 0);
  gpuMemset_realw (&accel_oc, size_oc, 0);
  gpuMemset_realw (&accel_ic, size_ic, 0);

  /*  or by kernel call - daniel todo: check if kernel call would be faster...

  int size = size_cm; // maximum size
  int blocksize = BLOCKSIZE_KERNEL1;
  int size_padded = ((int)ceil ( ( (double)size) / ( (double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    //launch kernel
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (cl_mem), (void *) &accel_cm.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (cl_mem), (void *) &accel_oc.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (cl_mem), (void *) &accel_ic.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (int), (void *) &size_cm));
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (int), (void *) &size_oc));
    clCheck (clSetKernelArg (mocl.kernels.update_disp_veloc_kernel, idx++, sizeof (int), (void *) &size_ic));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.update_displ_lddrk_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    //launch kernel
    update_displ_lddrk_kernel<<<grid,threads,0,mp->compute_stream>>>(accel_cm.cuda,
                                                                     accel_oc.cuda,
                                                                     accel_ic.cuda,
                                                                     size_cm,size_oc,size_ic);
  }
#endif
#ifdef USE_HIP
  // todo..
#endif
  */

  // kernel timing
  if (GPU_TIMING){ stop_timing_gpu(&start,&stop,"update_displ_lddrk_gpu"); }

  GPU_ERROR_CHECKING ("update_displ_lddrk_gpu");
}



/*----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_ (update_elastic_lddrk_gpu,
               UPDATE_ELASTIC_LDDRK_GPU) (long *Mesh_pointer_f,
                                          realw *alpha_lddrk_f, realw *beta_lddrk_f,
                                          int *FORWARD_OR_ADJOINT) {

  TRACE ("update_elastic_lddrk_gpu");

  //get Mesh from Fortran integer wrapper
  Mesh *mp = (Mesh *) *Mesh_pointer_f;
  realw alpha_lddrk = *alpha_lddrk_f;
  realw beta_lddrk = *beta_lddrk_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in update_elastic_lddrk_gpu() routine");
  }

  // update kernel
  int blocksize = BLOCKSIZE_KERNEL3;

  //updates displacement & velocity
  int size_padded, num_blocks_x, num_blocks_y;

  //crust/mantle region
  size_padded = ((int)ceil ( ( (double)mp->NGLOB_CRUST_MANTLE) / ( (double)blocksize)))*blocksize;
  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

  // sets gpu arrays
  realw deltat;
  gpu_realw_mem displ,veloc,accel;
  gpu_realw_mem displ_lddrk,veloc_lddrk;

  if (*FORWARD_OR_ADJOINT == 1) {
    displ = mp->d_displ_crust_mantle;
    veloc = mp->d_veloc_crust_mantle;
    accel = mp->d_accel_crust_mantle;
    displ_lddrk = mp->d_displ_crust_mantle_lddrk;
    veloc_lddrk = mp->d_veloc_crust_mantle_lddrk;
    deltat = mp->deltat;
  } else {
    //debug
    DEBUG_BACKWARD_UPDATE ();
    // for backward/reconstructed fields
    displ = mp->d_b_displ_crust_mantle;
    veloc = mp->d_b_veloc_crust_mantle;
    accel = mp->d_b_accel_crust_mantle;
    displ_lddrk = mp->d_b_displ_crust_mantle_lddrk;
    veloc_lddrk = mp->d_b_veloc_crust_mantle_lddrk;
    deltat = mp->b_deltat;
  }

  // kernel timing
  gpu_event start,stop;
  if (GPU_TIMING){ start_timing_gpu(&start,&stop); }

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clSetKernelArg (mocl.kernels.update_elastic_lddrk_kernel, idx++, sizeof (cl_mem), (void *) &displ.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_elastic_lddrk_kernel, idx++, sizeof (cl_mem), (void *) &veloc.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_elastic_lddrk_kernel, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_elastic_lddrk_kernel, idx++, sizeof (cl_mem), (void *) &displ_lddrk.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_elastic_lddrk_kernel, idx++, sizeof (cl_mem), (void *) &veloc_lddrk.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_elastic_lddrk_kernel, idx++, sizeof (realw), (void *) &alpha_lddrk));
    clCheck (clSetKernelArg (mocl.kernels.update_elastic_lddrk_kernel, idx++, sizeof (realw), (void *) &beta_lddrk));
    clCheck (clSetKernelArg (mocl.kernels.update_elastic_lddrk_kernel, idx++, sizeof (realw), (void *) &deltat));
    clCheck (clSetKernelArg (mocl.kernels.update_elastic_lddrk_kernel, idx++, sizeof (int), (void *) &mp->NGLOB_CRUST_MANTLE));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.update_elastic_lddrk_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    // graph
#ifdef USE_CUDA_GRAPHS
    if (! mp->use_graph_call_elastic){
#endif
    dim3 grid = dim3(num_blocks_x,num_blocks_y);
    dim3 threads = dim3(blocksize,1,1);

    // launches kernel
    update_elastic_lddrk_kernel<<< grid, threads,0,mp->compute_stream>>>(displ.cuda,
                                                                         veloc.cuda,
                                                                         accel.cuda,
                                                                         displ_lddrk.cuda,
                                                                         veloc_lddrk.cuda,
                                                                         alpha_lddrk,
                                                                         beta_lddrk,
                                                                         deltat,
                                                                         mp->NGLOB_CRUST_MANTLE);
#ifdef USE_CUDA_GRAPHS
    } // graph
#endif
  }
#endif
#ifdef USE_HIP
  if (run_hip) {
    dim3 grid = dim3(num_blocks_x,num_blocks_y);
    dim3 threads = dim3(blocksize,1,1);

    // launches kernel
    hipLaunchKernelGGL(HIP_KERNEL_NAME(update_elastic_lddrk_kernel), grid, threads, 0, mp->compute_stream,
                                                                     displ.hip,
                                                                     veloc.hip,
                                                                     accel.hip,
                                                                     displ_lddrk.hip,
                                                                     veloc_lddrk.hip,
                                                                     alpha_lddrk,
                                                                     beta_lddrk,
                                                                     deltat,
                                                                     mp->NGLOB_CRUST_MANTLE);
  }
#endif

  //inner core region
  size_padded = ((int)ceil ( ( (double)mp->NGLOB_INNER_CORE) / ((double)blocksize))) * blocksize;
  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

  // sets gpu arrays
  if (*FORWARD_OR_ADJOINT == 1) {
    displ = mp->d_displ_inner_core;
    veloc = mp->d_veloc_inner_core;
    accel = mp->d_accel_inner_core;
    displ_lddrk = mp->d_displ_inner_core_lddrk;
    veloc_lddrk = mp->d_veloc_inner_core_lddrk;
    deltat = mp->deltat;
  } else {
    //debug
    DEBUG_BACKWARD_UPDATE ();
    // for backward/reconstructed fields
    displ = mp->d_b_displ_inner_core;
    veloc = mp->d_b_veloc_inner_core;
    accel = mp->d_b_accel_inner_core;
    displ_lddrk = mp->d_b_displ_inner_core_lddrk;
    veloc_lddrk = mp->d_b_veloc_inner_core_lddrk;
    deltat = mp->b_deltat;
  }

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clSetKernelArg (mocl.kernels.update_elastic_lddrk_kernel, idx++, sizeof (cl_mem), (void *) &displ.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_elastic_lddrk_kernel, idx++, sizeof (cl_mem), (void *) &veloc.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_elastic_lddrk_kernel, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_elastic_lddrk_kernel, idx++, sizeof (cl_mem), (void *) &displ_lddrk.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_elastic_lddrk_kernel, idx++, sizeof (cl_mem), (void *) &veloc_lddrk.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_elastic_lddrk_kernel, idx++, sizeof (realw), (void *) &alpha_lddrk));
    clCheck (clSetKernelArg (mocl.kernels.update_elastic_lddrk_kernel, idx++, sizeof (realw), (void *) &beta_lddrk));
    clCheck (clSetKernelArg (mocl.kernels.update_elastic_lddrk_kernel, idx++, sizeof (realw), (void *) &deltat));
    clCheck (clSetKernelArg (mocl.kernels.update_elastic_lddrk_kernel, idx++, sizeof (int), (void *) &mp->NGLOB_INNER_CORE));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.update_elastic_lddrk_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    // graph
#ifdef USE_CUDA_GRAPHS
    if (! mp->use_graph_call_elastic){
#endif
    dim3 grid = dim3(num_blocks_x,num_blocks_y);
    dim3 threads = dim3(blocksize,1,1);

    // launches kernel
    update_elastic_lddrk_kernel<<< grid, threads,0,mp->compute_stream>>>(displ.cuda,
                                                                         veloc.cuda,
                                                                         accel.cuda,
                                                                         displ_lddrk.cuda,
                                                                         veloc_lddrk.cuda,
                                                                         alpha_lddrk,
                                                                         beta_lddrk,
                                                                         deltat,
                                                                         mp->NGLOB_INNER_CORE);
    // graph
#ifdef USE_CUDA_GRAPHS
    } // graph

    // finish creating graph
    if (mp->init_graph_elastic){
      // stop capturing
      print_CUDA_error_if_any(cudaStreamEndCapture(mp->compute_stream, &mp->graph_elastic),930);

      // get graph info
      size_t numNodes = 0;
      print_CUDA_error_if_any(cudaGraphGetNodes(mp->graph_elastic, NULL, &numNodes),931);
      //if (mp->myrank == 0) printf("\nGraph: elastic number of nodes = %zu\n",numNodes);

      print_CUDA_error_if_any(cudaGraphInstantiate(&mp->graphExec_elastic, mp->graph_elastic, NULL, NULL, 0),932);
      //if (mp->myrank == 0) printf("\nGraph: elastic instantiated\n");

      // graph is initialized, ready to be called by graph from now on
      mp->init_graph_elastic = 0;
      mp->use_graph_call_elastic = 1;
    }

    // launches graph instead of separate kernels
    if (mp->use_graph_call_elastic){
      // graph
      print_CUDA_error_if_any(cudaGraphLaunch(mp->graphExec_elastic, mp->compute_stream),935);
      //if (mp->myrank == 0) printf("\nGraph: elastic launch \n");
    }
#endif
  }
#endif
#ifdef USE_HIP
  if (run_hip) {
    dim3 grid = dim3(num_blocks_x,num_blocks_y);
    dim3 threads = dim3(blocksize,1,1);

    // launches kernel
    hipLaunchKernelGGL(HIP_KERNEL_NAME(update_elastic_lddrk_kernel), grid, threads, 0, mp->compute_stream,
                                                                     displ.hip,
                                                                     veloc.hip,
                                                                     accel.hip,
                                                                     displ_lddrk.hip,
                                                                     veloc_lddrk.hip,
                                                                     alpha_lddrk,
                                                                     beta_lddrk,
                                                                     deltat,
                                                                     mp->NGLOB_INNER_CORE);

  }
#endif

  // kernel timing
  if (GPU_TIMING){ stop_timing_gpu(&start,&stop,"update_elastic_lddrk_gpu"); }

  GPU_ERROR_CHECKING ("after update_elastic_lddrk_gpu");
}

/* ----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_ (update_acoustic_lddrk_gpu,
               UPDATE_ACOUSTIC_LDDRK_GPU) (long *Mesh_pointer_f,
                                           realw *alpha_lddrk_f, realw *beta_lddrk_f,
                                           int *FORWARD_OR_ADJOINT) {

  TRACE ("update_acoustic_lddrk_gpu");

  //get Mesh from Fortran integer wrapper
  Mesh *mp = (Mesh *) *Mesh_pointer_f;
  realw alpha_lddrk = *alpha_lddrk_f;
  realw beta_lddrk = *beta_lddrk_f;

  // safety check
  if (*FORWARD_OR_ADJOINT != 1 && *FORWARD_OR_ADJOINT != 3) {
    exit_on_error("Error invalid FORWARD_OR_ADJOINT in update_acoustic_lddrk_gpu() routine");
  }

  // update kernel
  int blocksize = BLOCKSIZE_KERNEL3;

  // outer core
  int size_padded = ((int)ceil ( ( (double)mp->NGLOB_OUTER_CORE) / ( (double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (size_padded/blocksize, &num_blocks_x, &num_blocks_y);

  // sets gpu arrays
  realw deltat;
  gpu_realw_mem displ,veloc,accel;
  gpu_realw_mem displ_lddrk,veloc_lddrk;
  if (*FORWARD_OR_ADJOINT == 1) {
    displ = mp->d_displ_outer_core;
    veloc = mp->d_veloc_outer_core;
    accel = mp->d_accel_outer_core;
    displ_lddrk = mp->d_displ_outer_core_lddrk;
    veloc_lddrk = mp->d_veloc_outer_core_lddrk;
    deltat = mp->deltat;
  } else {
    //debug
    DEBUG_BACKWARD_UPDATE ();
    // for backward/reconstructed fields
    displ = mp->d_b_displ_outer_core;
    veloc = mp->d_b_veloc_outer_core;
    accel = mp->d_b_accel_outer_core;
    displ_lddrk = mp->d_b_displ_outer_core_lddrk;
    veloc_lddrk = mp->d_b_veloc_outer_core_lddrk;
    deltat = mp->b_deltat;
  }

  // kernel timing
  gpu_event start,stop;
  if (GPU_TIMING){ start_timing_gpu(&start,&stop); }

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    //updates velocity
    clCheck (clSetKernelArg (mocl.kernels.update_acoustic_lddrk_kernel, idx++, sizeof (cl_mem), (void *) &displ.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_acoustic_lddrk_kernel, idx++, sizeof (cl_mem), (void *) &veloc.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_acoustic_lddrk_kernel, idx++, sizeof (cl_mem), (void *) &accel.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_acoustic_lddrk_kernel, idx++, sizeof (cl_mem), (void *) &displ_lddrk.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_acoustic_lddrk_kernel, idx++, sizeof (cl_mem), (void *) &veloc_lddrk.ocl));
    clCheck (clSetKernelArg (mocl.kernels.update_acoustic_lddrk_kernel, idx++, sizeof (realw), (void *) &alpha_lddrk));
    clCheck (clSetKernelArg (mocl.kernels.update_acoustic_lddrk_kernel, idx++, sizeof (realw), (void *) &beta_lddrk));
    clCheck (clSetKernelArg (mocl.kernels.update_acoustic_lddrk_kernel, idx++, sizeof (realw), (void *) &deltat));
    clCheck (clSetKernelArg (mocl.kernels.update_acoustic_lddrk_kernel, idx++, sizeof (int), (void *) &mp->NGLOB_OUTER_CORE));

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.update_acoustic_lddrk_kernel, 2, NULL,
                                     global_work_size, local_work_size, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    // graph
#ifdef USE_CUDA_GRAPHS
    if (! mp->use_graph_call_acoustic){
#endif
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // updates velocity
    update_acoustic_lddrk_kernel<<< grid, threads,0,mp->compute_stream>>>(displ.cuda,
                                                                          veloc.cuda,
                                                                          accel.cuda,
                                                                          displ_lddrk.cuda,
                                                                          veloc_lddrk.cuda,
                                                                          alpha_lddrk,
                                                                          beta_lddrk,
                                                                          deltat,
                                                                          mp->NGLOB_OUTER_CORE);
    // graph
#ifdef USE_CUDA_GRAPHS
    } // graph

    // finish creating graph
    if (mp->init_graph_acoustic){
      // stop capturing
      print_CUDA_error_if_any(cudaStreamEndCapture(mp->compute_stream, &mp->graph_acoustic),920);

      // get graph info
      size_t numNodes = 0;
      print_CUDA_error_if_any(cudaGraphGetNodes(mp->graph_acoustic, NULL, &numNodes),921);
      //if (mp->myrank == 0) printf("\nGraph: acoustic number of nodes = %zu\n",numNodes);

      print_CUDA_error_if_any(cudaGraphInstantiate(&mp->graphExec_acoustic, mp->graph_acoustic, NULL, NULL, 0),922);
      //if (mp->myrank == 0) printf("\nGraph: acoustic instantiated\n");

      // graph is initialized, ready to be called by graph from now on
      mp->init_graph_acoustic = 0;
      mp->use_graph_call_acoustic = 1;
    }

    // launches graph instead of separate kernels
    if (mp->use_graph_call_acoustic){
      // graph
      print_CUDA_error_if_any(cudaGraphLaunch(mp->graphExec_acoustic, mp->compute_stream),925);
      //if (mp->myrank == 0) printf("\nGraph: acoustic launch \n");
    }
#endif
  }
#endif
#ifdef USE_HIP
  if (run_hip) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // updates velocity
    hipLaunchKernelGGL(HIP_KERNEL_NAME(update_acoustic_lddrk_kernel), grid, threads, 0, mp->compute_stream,
                                                                      displ.hip,
                                                                      veloc.hip,
                                                                      accel.hip,
                                                                      displ_lddrk.hip,
                                                                      veloc_lddrk.hip,
                                                                      alpha_lddrk,
                                                                      beta_lddrk,
                                                                      deltat,
                                                                      mp->NGLOB_OUTER_CORE);
  }
#endif

  // kernel timing
  if (GPU_TIMING){ stop_timing_gpu(&start,&stop,"update_acoustic_lddrk_gpu"); }

  GPU_ERROR_CHECKING ("after update_acoustic_lddrk_gpu");
}
