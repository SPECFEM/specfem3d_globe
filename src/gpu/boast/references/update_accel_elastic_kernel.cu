// from update_displacement_cuda.cu
typedef float realw;

__global__ void update_accel_elastic_kernel(realw* accel,
                                            realw* veloc,
                                            int size,
                                            realw two_omega_earth,
                                            realw* rmassx,
                                            realw* rmassy,
                                            realw* rmassz) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // because of block and grid sizing problems, there is a small
  // amount of buffer at the end of the calculation
  if(id < size) {
    // note: update adds rotational acceleration in case two_omega_earth is non-zero
    accel[3*id] = accel[3*id]*rmassx[id] + two_omega_earth*veloc[3*id + 1]; // (2,i);
    accel[3*id + 1] = accel[3*id + 1]*rmassy[id] - two_omega_earth*veloc[3*id]; //(1,i);
    accel[3*id + 2] = accel[3*id + 2]*rmassz[id];
  }
}
