// from noise_tomography_cuda.cu
#define NGLL3 125

typedef float realw;

__global__ void noise_add_source_master_rec_kernel(int* ibool,
                                                   int* ispec_selected_rec,
                                                   int irec_master_noise,
                                                   realw* accel,
                                                   realw* noise_sourcearray,
                                                   int it) {
  int tx = threadIdx.x;
  int ispec = ispec_selected_rec[irec_master_noise]-1;
  int iglob = ibool[tx + NGLL3*ispec]-1;

  atomicAdd(&accel[iglob*3  ],noise_sourcearray[  3*tx + 3*NGLL3*it]);
  atomicAdd(&accel[iglob*3+1],noise_sourcearray[1+3*tx + 3*NGLL3*it]);
  atomicAdd(&accel[iglob*3+2],noise_sourcearray[2+3*tx + 3*NGLL3*it]);
}
