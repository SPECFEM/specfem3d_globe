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
__global__ void compute_coupling_ocean_kernel(float * accel_crust_mantle, const float * rmassx_crust_mantle, const float * rmassy_crust_mantle, const float * rmassz_crust_mantle, const float * rmass_ocean_load, const int npoin_ocean_load, const int * ibool_ocean_load, const float * normal_ocean_load){
  int ipoin;
  int iglob;
  float nx;
  float ny;
  float nz;
  float rmass;
  float force_normal_comp;
  float additional_term_x;
  float additional_term_y;
  float additional_term_z;
  ipoin = threadIdx.x + (blockIdx.x) * (blockDim.x) + ((gridDim.x) * (blockDim.x)) * (threadIdx.y + (blockIdx.y) * (blockDim.y));
  if(ipoin < npoin_ocean_load){
    iglob = ibool_ocean_load[ipoin - 0] - (1);
    nx = normal_ocean_load[INDEX2(NDIM, 0, ipoin) - 0];
    ny = normal_ocean_load[INDEX2(NDIM, 1, ipoin) - 0];
    nz = normal_ocean_load[INDEX2(NDIM, 2, ipoin) - 0];
    force_normal_comp = ((accel_crust_mantle[0 - 0 + (iglob - (0)) * (3)]) * (nx)) / (rmassx_crust_mantle[iglob - 0]) + ((accel_crust_mantle[1 - 0 + (iglob - (0)) * (3)]) * (ny)) / (rmassy_crust_mantle[iglob - 0]) + ((accel_crust_mantle[2 - 0 + (iglob - (0)) * (3)]) * (nz)) / (rmassz_crust_mantle[iglob - 0]);
    rmass = rmass_ocean_load[ipoin - 0];
    additional_term_x = (rmass - (rmassx_crust_mantle[iglob - 0])) * (force_normal_comp);
    additional_term_y = (rmass - (rmassy_crust_mantle[iglob - 0])) * (force_normal_comp);
    additional_term_z = (rmass - (rmassz_crust_mantle[iglob - 0])) * (force_normal_comp);
    accel_crust_mantle[0 - 0 + (iglob - (0)) * (3)] = accel_crust_mantle[0 - 0 + (iglob - (0)) * (3)] + (additional_term_x) * (nx);
    accel_crust_mantle[1 - 0 + (iglob - (0)) * (3)] = accel_crust_mantle[1 - 0 + (iglob - (0)) * (3)] + (additional_term_y) * (ny);
    accel_crust_mantle[2 - 0 + (iglob - (0)) * (3)] = accel_crust_mantle[2 - 0 + (iglob - (0)) * (3)] + (additional_term_z) * (nz);
  }
}
