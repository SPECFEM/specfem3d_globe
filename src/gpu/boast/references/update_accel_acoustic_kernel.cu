// from update_displacement_cuda.cu
typedef float realw;

__global__ void update_accel_acoustic_kernel(realw* accel,
                                             int size,
                                             realw* rmass) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if(id < size) {
    // multiplies pressure with the inverse of the mass matrix
    accel[id] = accel[id]*rmass[id];
  }
}
