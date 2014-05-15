// from write_seismograms_cuda.cu
#define NGLL3 125

typedef float realw;

__global__ void write_seismograms_transfer_strain_from_device_kernel(int* number_receiver_global,
                                                                     int* ispec_selected_rec,
                                                                     int* ibool,
                                                                     realw* station_strain_field,
                                                                     realw* d_field,
                                                                     int nrec_local) {

// scalar fields

  int blockID = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;

  int irec,ispec,iglob;

  if(blockID < nrec_local) {
    irec = number_receiver_global[blockID]-1;
    ispec = ispec_selected_rec[irec]-1;
    iglob = ibool[tx + NGLL3*ispec]-1;

    station_strain_field[NGLL3*blockID + tx] = d_field[iglob];
  }
}
