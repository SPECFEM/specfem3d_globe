#ifndef INDEX2
#define INDEX2(xsize,x,y) x + (y)*xsize
#endif
#ifndef INDEX3
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)
#endif
#ifndef INDEX4
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))
#endif
#ifndef INDEX5
#define INDEX5(xsize,ysize,zsize,isize,x,y,z,i,j) x + xsize*(y + ysize*(z + zsize*(i + isize*(j))))
#endif
#ifndef NDIM
#define NDIM 3
#endif
#ifndef NGLLX
#define NGLLX 5
#endif
#ifndef NGLL2
#define NGLL2 25
#endif
#ifndef NGLL3
#define NGLL3 125
#endif
#ifndef NGLL3_PADDED
#define NGLL3_PADDED 128
#endif
#ifndef N_SLS
#define N_SLS 3
#endif
#ifndef IREGION_CRUST_MANTLE
#define IREGION_CRUST_MANTLE 1
#endif
#ifndef IREGION_INNER_CORE
#define IREGION_INNER_CORE 3
#endif
#ifndef IFLAG_IN_FICTITIOUS_CUBE
#define IFLAG_IN_FICTITIOUS_CUBE 11
#endif
#ifndef R_EARTH_KM
#define R_EARTH_KM 6371.0f
#endif
#ifndef COLORING_MIN_NSPEC_INNER_CORE
#define COLORING_MIN_NSPEC_INNER_CORE 1000
#endif
#ifndef COLORING_MIN_NSPEC_OUTER_CORE
#define COLORING_MIN_NSPEC_OUTER_CORE 1000
#endif
#ifndef BLOCKSIZE_TRANSFER
#define BLOCKSIZE_TRANSFER 256
#endif
__global__ void write_seismograms_transfer_from_device_kernel(const int * number_receiver_global, const int * ispec_selected_rec, const int * ibool, float * station_seismo_field, const float * d_field, const int nrec_local){
  int tx;
  int irec;
  int ispec;
  int iglob;
  int blockID;
  blockID = blockIdx.x + (blockIdx.y) * (gridDim.x);
  tx = threadIdx.x;
  if(blockID < nrec_local){
    irec = number_receiver_global[blockID - 0] - (1);
    ispec = ispec_selected_rec[irec - 0] - (1);
    iglob = ibool[tx + (NGLL3) * (ispec) - 0] - (1);
    station_seismo_field[((NGLL3) * (3)) * (blockID) + (tx) * (3) + 0 - 0] = d_field[(iglob) * (3) + 0 - 0];
    station_seismo_field[((NGLL3) * (3)) * (blockID) + (tx) * (3) + 1 - 0] = d_field[(iglob) * (3) + 1 - 0];
    station_seismo_field[((NGLL3) * (3)) * (blockID) + (tx) * (3) + 2 - 0] = d_field[(iglob) * (3) + 2 - 0];
  }
}
