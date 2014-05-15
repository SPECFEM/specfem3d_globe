// from update_displacement_cuda.cu
typedef float realw;

__global__ void update_veloc_elastic_kernel(realw* veloc,
                                            realw* accel,
                                            int size,
                                            realw deltatover2) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if(id < size) {
    veloc[3*id] = veloc[3*id] + deltatover2*accel[3*id];
    veloc[3*id + 1] = veloc[3*id + 1] + deltatover2*accel[3*id + 1];
    veloc[3*id + 2] = veloc[3*id + 2] + deltatover2*accel[3*id + 2];
  }
}
