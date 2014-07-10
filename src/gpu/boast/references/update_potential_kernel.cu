// from update_displacement_cuda.cu
typedef float realw;

__global__ void update_potential_kernel(realw* potential_acoustic,
                                        realw* potential_dot_acoustic,
                                        realw* potential_dot_dot_acoustic,
                                        int size,
                                        realw deltat,
                                        realw deltatsqover2,
                                        realw deltatover2) {
  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  /* because of block and grid sizing problems, there is a small */
  /* amount of buffer at the end of the calculation */
  if(id < size) {
    potential_acoustic[id] = potential_acoustic[id]
                            + deltat*potential_dot_acoustic[id]
                            + deltatsqover2*potential_dot_dot_acoustic[id];

    potential_dot_acoustic[id] = potential_dot_acoustic[id]
                                + deltatover2*potential_dot_dot_acoustic[id];

    potential_dot_dot_acoustic[id] = 0.0f;
  }
}
