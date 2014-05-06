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

void write_seismograms_transfer_from_device (Mesh *mp,
                                             gpu_realw_mem *d_field,
                                             realw *h_field,
                                             int *number_receiver_global,
                                             gpu_int_mem *d_ispec_selected,
                                             int *h_ispec_selected,
                                             int *ibool) {

  TRACE ("write_seismograms_transfer_from_device");

  int irec_local, irec;
  int ispec, iglob, i;

  //checks if anything to do

  if (mp->nrec_local == 0)
    return;

  int blocksize = NGLL3;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->nrec_local, &num_blocks_x, &num_blocks_y);

  //prepare field transfer array on device

#ifdef USE_OPENCL
  if (run_opencl) {
    if (GPU_ASYNC_COPY) {
      clCheck(clFlush(mocl.copy_queue));
    }
    
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;
    
    clCheck (clSetKernelArg (mocl.kernels.write_seismograms_transfer_from_device_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_number_receiver_global.ocl));
    clCheck (clSetKernelArg (mocl.kernels.write_seismograms_transfer_from_device_kernel, idx++, sizeof (cl_mem), (void *) &d_ispec_selected->ocl));
    clCheck (clSetKernelArg (mocl.kernels.write_seismograms_transfer_from_device_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.write_seismograms_transfer_from_device_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_station_seismo_field.ocl));
    clCheck (clSetKernelArg (mocl.kernels.write_seismograms_transfer_from_device_kernel, idx++, sizeof (cl_mem), (void *) &d_field->ocl));
    clCheck (clSetKernelArg (mocl.kernels.write_seismograms_transfer_from_device_kernel, idx++, sizeof (int), (void *) &mp->nrec_local));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;


    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.write_seismograms_transfer_from_device_kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL));

    //copies array to CPU
    if (GPU_ASYNC_COPY) {
      clCheck(clFlush(mocl.command_queue));
      clCheck (clEnqueueReadBuffer (mocl.copy_queue, mp->d_station_seismo_field.ocl, CL_FALSE, 0,
                                    3 * NGLL3 * mp->nrec_local * sizeof (realw),
                                    mp->h_station_seismo_field, 0, NULL, NULL));
    } else {
      clCheck (clEnqueueReadBuffer (mocl.command_queue, mp->d_station_seismo_field.ocl, CL_TRUE, 0,
                                    3 * NGLL3 * mp->nrec_local * sizeof (realw),
                                    mp->h_station_seismo_field, 0, NULL, NULL));
    }
  }
#endif
#if USE_CUDA
  if (run_cuda) {
    // waits for previous copy call to be finished
    if( GPU_ASYNC_COPY ){
      cudaStreamSynchronize(mp->copy_stream);
    }
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // prepare field transfer array on device
    write_seismograms_transfer_from_device_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_number_receiver_global.cuda,
                                                                                         d_ispec_selected->cuda,
                                                                                         mp->d_ibool_crust_mantle.cuda,
                                                                                         mp->d_station_seismo_field.cuda,
                                                                                         d_field->cuda,
                                                                                         mp->nrec_local);

 // copies array to CPU
  if( GPU_ASYNC_COPY ){
    // waits for previous compute call to be finished
    cudaStreamSynchronize(mp->compute_stream);

    // asynchronous copy
    // note: we need to update the host array in a subsequent call to transfer_seismo_from_device_async() routine
    cudaMemcpyAsync(mp->h_station_seismo_field,mp->d_station_seismo_field.cuda,
                    3*NGLL3*(mp->nrec_local)*sizeof(realw),
                    cudaMemcpyDeviceToHost,mp->copy_stream);
  }else{
    // synchronous copy
    print_CUDA_error_if_any(cudaMemcpy(mp->h_station_seismo_field,mp->d_station_seismo_field.cuda,
                                       3*NGLL3*(mp->nrec_local)*sizeof(realw),cudaMemcpyDeviceToHost),77000);

	}
  }
#endif
  if (!GPU_ASYNC_COPY) {
    for (irec_local = 0; irec_local < mp->nrec_local; irec_local++) {
      irec = number_receiver_global[irec_local] - 1;
      ispec = h_ispec_selected[irec] - 1;
      
      for (i = 0; i < NGLL3; i++) {
        iglob = ibool[i+NGLL3*ispec] - 1;
        h_field[0+3*iglob] = mp->h_station_seismo_field[0+3*i+irec_local*NGLL3*3];
        h_field[1+3*iglob] = mp->h_station_seismo_field[1+3*i+irec_local*NGLL3*3];
        h_field[2+3*iglob] = mp->h_station_seismo_field[2+3*i+irec_local*NGLL3*3];
      }
    }
  }
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_gpu_error ("write_seismograms_transfer_from_device");
#endif
}

/*----------------------------------------------------------------------------------------------- */

void write_seismograms_transfer_strain_from_device (Mesh *mp,
                                                    gpu_realw_mem *d_field,
                                                    realw *h_field,
                                                    int *number_receiver_global,
                                                    gpu_int_mem *d_ispec_selected,
                                                    int *h_ispec_selected,
                                                    int *ibool) {

  TRACE ("write_seismograms_transfer_strain_from_device");

  int irec_local, irec;
  int ispec, iglob, i;

  //checks if anything to do
  if (mp->nrec_local == 0)
    return;

  int blocksize = NGLL3;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy (mp->nrec_local, &num_blocks_x, &num_blocks_y);

#ifdef USE_OPENCL
  if (run_opencl) {
    size_t global_work_size[2];
    size_t local_work_size[2];
    cl_uint idx = 0;
    
    //prepare field transfer array on device
    clCheck (clSetKernelArg (mocl.kernels.write_seismograms_transfer_strain_from_device_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_number_receiver_global.ocl));
    clCheck (clSetKernelArg (mocl.kernels.write_seismograms_transfer_strain_from_device_kernel, idx++, sizeof (cl_mem), (void *) &d_ispec_selected->ocl));
    clCheck (clSetKernelArg (mocl.kernels.write_seismograms_transfer_strain_from_device_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_ibool_crust_mantle.ocl));
    clCheck (clSetKernelArg (mocl.kernels.write_seismograms_transfer_strain_from_device_kernel, idx++, sizeof (cl_mem), (void *) &mp->d_station_seismo_field.ocl));
    clCheck (clSetKernelArg (mocl.kernels.write_seismograms_transfer_strain_from_device_kernel, idx++, sizeof (cl_mem), (void *) &d_field->ocl));
    clCheck (clSetKernelArg (mocl.kernels.write_seismograms_transfer_strain_from_device_kernel, idx++, sizeof (int), (void *) &mp->nrec_local));

    local_work_size[0] = blocksize;
    local_work_size[1] = 1;
    global_work_size[0] = num_blocks_x * blocksize;
    global_work_size[1] = num_blocks_y;

    clCheck (clEnqueueNDRangeKernel (mocl.command_queue, mocl.kernels.write_seismograms_transfer_strain_from_device_kernel, 2, NULL, global_work_size, local_work_size, 0, NULL, NULL));


    //copies array to CPU
    clCheck (clEnqueueReadBuffer (mocl.command_queue, mp->d_station_seismo_field.ocl, CL_TRUE, 0,
                                  NGLL3 * mp->nrec_local * sizeof (realw),
                                  mp->h_station_seismo_field, 0, NULL, NULL));
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    // prepare field transfer array on device
    write_seismograms_transfer_strain_from_device_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_number_receiver_global.cuda,
                                                                                                d_ispec_selected->cuda,
                                                                                                mp->d_ibool_crust_mantle.cuda,
                                                                                                mp->d_station_strain_field.cuda,
                                                                                                d_field->cuda,
                                                                                                mp->nrec_local);

    // copies array to CPU
    // synchronous copy
    print_CUDA_error_if_any(cudaMemcpy(mp->h_station_strain_field,mp->d_station_strain_field.cuda,
                                       NGLL3*(mp->nrec_local)*sizeof(realw),cudaMemcpyDeviceToHost),77001);
  }
#endif
  // updates host array
  for (irec_local = 0; irec_local < mp->nrec_local; irec_local++) {
    irec = number_receiver_global[irec_local] - 1;
    ispec = h_ispec_selected[irec] - 1;
    for (i = 0; i < NGLL3; i++) {
      iglob = ibool[i+NGLL3*ispec] - 1;
      h_field[iglob] = mp->h_station_strain_field[i+irec_local*NGLL3];
    }
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_gpu_error ("write_seismograms_transfer_strain_from_device");
#endif
}

/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (write_seismograms_transfer_gpu,
               WRITE_SEISMOGRAMS_TRANSFER_GPU) (long *Mesh_pointer_f,
                                                realw *displ,
                                                realw *b_displ,
                                                realw *eps_trace_over_3,
                                                realw *epsilondev_xx,
                                                realw *epsilondev_yy,
                                                realw *epsilondev_xy,
                                                realw *epsilondev_xz,
                                                realw *epsilondev_yz,
                                                int *number_receiver_global,
                                                int *ispec_selected_rec,
                                                int *ispec_selected_source,
                                                int *ibool) {
  TRACE ("write_seismograms_transfer_gpu");

  //get Mesh from Fortran integer wrapper
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  //checks if anything to do
  if (mp->nrec_local == 0) return;

  // transfers displacement values in receiver elements from GPU to CPU
  switch (mp->simulation_type) {
  case 1:
    write_seismograms_transfer_from_device (mp,
                                            &mp->d_displ_crust_mantle,
                                            displ,
                                            number_receiver_global,
                                            &mp->d_ispec_selected_rec,
                                            ispec_selected_rec,
                                            ibool);
    break;

  case 2:
    write_seismograms_transfer_from_device (mp,
                                            &mp->d_displ_crust_mantle,
                                            displ,
                                            number_receiver_global,
                                            &mp->d_ispec_selected_source,
                                            ispec_selected_source,
                                            ibool);

    write_seismograms_transfer_strain_from_device (mp,
                                                   &mp->d_eps_trace_over_3_crust_mantle,
                                                   eps_trace_over_3,
                                                   number_receiver_global,
                                                   &mp->d_ispec_selected_source,
                                                   ispec_selected_source,
                                                   ibool);
    write_seismograms_transfer_strain_from_device (mp,
                                                   &mp->d_epsilondev_xx_crust_mantle,
                                                   epsilondev_xx,
                                                   number_receiver_global,
                                                   &mp->d_ispec_selected_source,
                                                   ispec_selected_source,
                                                   ibool);
    write_seismograms_transfer_strain_from_device (mp,
                                                   &mp->d_epsilondev_yy_crust_mantle,
                                                   epsilondev_yy,
                                                   number_receiver_global,
                                                   &mp->d_ispec_selected_source,
                                                   ispec_selected_source,
                                                   ibool);
    write_seismograms_transfer_strain_from_device (mp,
                                                   &mp->d_epsilondev_xy_crust_mantle,
                                                   epsilondev_xy,
                                                   number_receiver_global,
                                                   &mp->d_ispec_selected_source,
                                                   ispec_selected_source,
                                                   ibool);
    write_seismograms_transfer_strain_from_device (mp,
                                                   &mp->d_epsilondev_xz_crust_mantle,
                                                   epsilondev_xz,
                                                   number_receiver_global,
                                                   &mp->d_ispec_selected_source,
                                                   ispec_selected_source,
                                                   ibool);
    write_seismograms_transfer_strain_from_device (mp,
                                                   &mp->d_epsilondev_yz_crust_mantle,
                                                   epsilondev_yz,
                                                   number_receiver_global,
                                                   &mp->d_ispec_selected_source,
                                                   ispec_selected_source,
                                                   ibool);
    break;

  case 3:
    write_seismograms_transfer_from_device (mp,
                                            &mp->d_b_displ_crust_mantle,
                                            b_displ,
                                            number_receiver_global,
                                            &mp->d_ispec_selected_rec,
                                            ispec_selected_rec,
                                            ibool);
    break;
  }
}

/* ----------------------------------------------------------------------------------------------- */

// data transfer to CPU host

/* ----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_(transfer_seismo_from_device_async,
              TRANSFER_SEISMO_FROM_DEVICE_ASYNC)(long* Mesh_pointer_f,
                                                 realw* displ,
                                                 realw* b_displ,
                                                 int* number_receiver_global,
                                                 int* ispec_selected_rec,
                                                 int* ispec_selected_source,
                                                 int* ibool) {

  TRACE("transfer_seismo_from_device_async");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from Fortran integer wrapper

  int irec,ispec,iglob,i;
  realw* h_field;
  int* h_ispec_selected;

  // checks if anything to do
  if (mp->nrec_local == 0) {
    return;
  }

  // checks async-memcpy
  if (GPU_ASYNC_COPY ==  0){
    exit_on_error("transfer_seismo_from_device_async must be called with GPU_ASYNC_COPY == 1, please check mesh_constants_cuda.h");
  }

  // waits for previous copy call to be finished
#ifdef USE_CUDA
  if (run_cuda) {
    cudaStreamSynchronize(mp->copy_stream);
  }
#endif

  // transfers displacements
  // select target array on host
  switch (mp->simulation_type) {
    case 1:
      // forward simulation
      h_field = displ;
      h_ispec_selected = ispec_selected_rec;
      break;

    case 2:
      // adjoint simulation
      h_field = displ;
      h_ispec_selected = ispec_selected_source;
      break;

    case 3:
      // kernel simulation
      h_field = b_displ;
      h_ispec_selected = ispec_selected_rec;
      break;
  }

  // updates corresponding array on CPU
  int irec_local;
  for (irec_local = 0 ; irec_local < mp->nrec_local; irec_local++) {
    irec = number_receiver_global[irec_local] - 1;
    ispec = h_ispec_selected[irec] - 1;

    for (i = 0; i < NGLL3; i++) {
      iglob = ibool[i+NGLL3*ispec] - 1;
      h_field[0+3*iglob] = mp->h_station_seismo_field[0+3*i+irec_local*NGLL3*3];
      h_field[1+3*iglob] = mp->h_station_seismo_field[1+3*i+irec_local*NGLL3*3];
      h_field[2+3*iglob] = mp->h_station_seismo_field[2+3*i+irec_local*NGLL3*3];
    }
  }

}
