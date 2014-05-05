// from assemble_MPI_scalar_cuda.cu
typedef float realw;

__global__ void prepare_boundary_potential_on_device(realw* d_potential_dot_dot_acoustic,
                                                     realw* d_send_potential_dot_dot_buffer,
                                                     int num_interfaces,
                                                     int max_nibool_interfaces,
                                                     int* d_nibool_interfaces,
                                                     int* d_ibool_interfaces) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  int iglob,iloc;

  for( int iinterface=0; iinterface < num_interfaces; iinterface++) {
    if(id<d_nibool_interfaces[iinterface]) {

      iloc = id + max_nibool_interfaces*iinterface;
      iglob = d_ibool_interfaces[iloc] - 1;

      // fills buffer
      d_send_potential_dot_dot_buffer[iloc] = d_potential_dot_dot_acoustic[iglob];
    }
  }

}
