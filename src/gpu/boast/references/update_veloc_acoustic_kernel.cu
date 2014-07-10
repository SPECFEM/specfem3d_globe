// from update_displacement_cuda.cu
typedef float realw;

__global__ void update_veloc_acoustic_kernel(realw* veloc,
                                             realw* accel,
                                             int size,
                                             realw deltatover2) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if(id < size) {
    // Newmark time scheme: corrector term
    veloc[id] = veloc[id] + deltatover2*accel[id];
  }
}
