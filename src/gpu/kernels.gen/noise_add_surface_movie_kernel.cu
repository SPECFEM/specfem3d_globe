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
__global__ void noise_add_surface_movie_kernel(float * accel, const int * ibool, const int * ibelm_top, const int nspec_top, const float * noise_surface_movie, const float * normal_x_noise, const float * normal_y_noise, const float * normal_z_noise, const float * mask_noise, const float * jacobian2D, const float * wgllwgll){
  int igll;
  int iface;
  igll = threadIdx.x;
  iface = blockIdx.x + (blockIdx.y) * (gridDim.x);
  if(iface < nspec_top){
    int i;
    int j;
    int k;
    int ispec;
    int iglob;
    int ipoin;
    float eta;
    float jacobianw;
    float normal_x;
    float normal_y;
    float normal_z;
    ispec = ibelm_top[iface - 0] - (1);
    k = NGLLX - (1);
    j = (igll) / (NGLLX);
    i = igll - ((j) * (NGLLX));
    iglob = ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] - (1);
    ipoin = (NGLL2) * (iface) + igll;
    normal_x = normal_x_noise[ipoin - 0];
    normal_y = normal_y_noise[ipoin - 0];
    normal_z = normal_z_noise[ipoin - 0];
    eta = 0.0f;
    eta = eta + (noise_surface_movie[INDEX3(NDIM, NGLL2, 0, igll, iface) - 0]) * (normal_x);
    eta = eta + (noise_surface_movie[INDEX3(NDIM, NGLL2, 1, igll, iface) - 0]) * (normal_y);
    eta = eta + (noise_surface_movie[INDEX3(NDIM, NGLL2, 2, igll, iface) - 0]) * (normal_z);
    jacobianw = (wgllwgll[(k) * (NGLLX) + i - 0]) * (jacobian2D[igll + (NGLL2) * (iface) - 0]);
    atomicAdd(accel + (iglob) * (3) + 0, (((eta) * (mask_noise[ipoin - 0])) * (normal_x)) * (jacobianw));
    atomicAdd(accel + (iglob) * (3) + 1, (((eta) * (mask_noise[ipoin - 0])) * (normal_y)) * (jacobianw));
    atomicAdd(accel + (iglob) * (3) + 2, (((eta) * (mask_noise[ipoin - 0])) * (normal_z)) * (jacobianw));
  }
}
