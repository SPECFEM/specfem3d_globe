// from write_seismograms_cuda.cu
#define NGLL3 125

typedef float realw;

__global__ void write_seismograms_transfer_from_device_kernel(int* number_receiver_global,
                                                              int* ispec_selected_rec,
                                                              int* ibool,
                                                              realw* station_seismo_field,
                                                              realw* d_field,
                                                              int nrec_local) {

// vector fields

  int blockID = blockIdx.x + blockIdx.y*gridDim.x;
  int tx = threadIdx.x;

  int irec,ispec,iglob;

  if(blockID < nrec_local) {
    irec = number_receiver_global[blockID]-1;
    ispec = ispec_selected_rec[irec]-1;
    iglob = ibool[tx + NGLL3*ispec]-1;

    station_seismo_field[3*NGLL3*blockID + 3*tx+0] = d_field[3*iglob];
    station_seismo_field[3*NGLL3*blockID + 3*tx+1] = d_field[3*iglob+1];
    station_seismo_field[3*NGLL3*blockID + 3*tx+2] = d_field[3*iglob+2];
  }
}
