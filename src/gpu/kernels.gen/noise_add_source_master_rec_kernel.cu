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
__global__ void noise_add_source_master_rec_kernel(const int * ibool, const int * ispec_selected_rec, const int irec_master_noise, float * accel, const float * noise_sourcearray, const int it){
  int tx;
  int ispec;
  int iglob;
  tx = threadIdx.x;
  ispec = ispec_selected_rec[irec_master_noise - 0] - (1);
  iglob = ibool[tx + (NGLL3) * (ispec) - 0] - (1);
  atomicAdd(accel + (iglob) * (3) + 0, noise_sourcearray[(tx) * (3) + ((NGLL3) * (3)) * (it) + 0 - 0]);
  atomicAdd(accel + (iglob) * (3) + 1, noise_sourcearray[(tx) * (3) + ((NGLL3) * (3)) * (it) + 1 - 0]);
  atomicAdd(accel + (iglob) * (3) + 2, noise_sourcearray[(tx) * (3) + ((NGLL3) * (3)) * (it) + 2 - 0]);
}
