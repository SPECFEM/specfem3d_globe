// from update_displacement_cuda.cu
typedef float realw;

__global__ void update_disp_veloc_kernel(realw* displ,
                                         realw* veloc,
                                         realw* accel,
                                         int size,
                                         realw deltat,
                                         realw deltatsqover2,
                                         realw deltatover2) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if(id < size) {
    displ[id] = displ[id] + deltat*veloc[id] + deltatsqover2*accel[id];
    veloc[id] = veloc[id] + deltatover2*accel[id];
    accel[id] = 0.0f; // can do this using memset...not sure if faster
  }
}
