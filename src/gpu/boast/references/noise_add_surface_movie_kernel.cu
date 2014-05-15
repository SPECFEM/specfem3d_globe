// from noise_tomography_cuda.cu
#define NDIM 3
#define NGLLX 5
#define NGLL2 25
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))

typedef float realw;

__global__ void noise_add_surface_movie_kernel(realw* accel,
                                               int* ibool,
                                               int* ibelm_top,
                                               int nspec_top,
                                               realw* noise_surface_movie,
                                               realw* normal_x_noise,
                                               realw* normal_y_noise,
                                               realw* normal_z_noise,
                                               realw* mask_noise,
                                               realw* jacobian2D,
                                               realw* wgllwgll) {

  int igll = threadIdx.x;
  int iface = blockIdx.x + gridDim.x*blockIdx.y; // surface element id

  // when nspec_top > MAXIMUM_GRID_DIM, but mod(nspec_top,2) > 0, we end up with an extra block.
  if(iface < nspec_top) {

    int ispec = ibelm_top[iface]-1;

    int k = NGLLX - 1;
    int j = (igll/NGLLX);
    int i = (igll-j*NGLLX);

    int iglob = ibool[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

    int ipoin = NGLL2*iface + igll;
    realw normal_x = normal_x_noise[ipoin];
    realw normal_y = normal_y_noise[ipoin];
    realw normal_z = normal_z_noise[ipoin];

    realw eta = (noise_surface_movie[INDEX3(NDIM,NGLL2,0,igll,iface)]*normal_x +
                 noise_surface_movie[INDEX3(NDIM,NGLL2,1,igll,iface)]*normal_y +
                 noise_surface_movie[INDEX3(NDIM,NGLL2,2,igll,iface)]*normal_z);

    // weighted jacobian
    realw jacobianw = wgllwgll[k*NGLLX+i]*jacobian2D[igll+NGLL2*iface];

    // note: check error from cuda-memcheck and ddt seems "incorrect", because we
    //          are passing a __constant__ variable pointer around like it was
    //          made using cudaMalloc, which *may* be "incorrect", but produces
    //          correct results.

    // note: global version uses jacobian2D arrays which do not include gll weights wgllwgll,
    //          thus we have to explicitly add: wgllwgll(..) * jacobian2D(..)

    atomicAdd(&accel[iglob*3]  ,eta*mask_noise[ipoin]*normal_x*jacobianw);
    atomicAdd(&accel[iglob*3+1],eta*mask_noise[ipoin]*normal_y*jacobianw);
    atomicAdd(&accel[iglob*3+2],eta*mask_noise[ipoin]*normal_z*jacobianw);

  }
}
