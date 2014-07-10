// from compute_kernels_cuda.cu
#define NDIM 3
#define NGLLX 5
#define NGLL2 25
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))

typedef float realw;

__global__ void compute_strength_noise_kernel(realw* displ,
                                              int* ibelm_top,
                                              int* ibool,
                                              realw* noise_surface_movie,
                                              realw* normal_x_noise,
                                              realw* normal_y_noise,
                                              realw* normal_z_noise,
                                              realw* Sigma_kl,
                                              realw deltat,
                                              int nspec_top) {
  int iface = blockIdx.x + blockIdx.y*gridDim.x;

  if(iface < nspec_top) {

    int ispec = ibelm_top[iface]-1;
    int igll = threadIdx.x;
    int ipoin = igll + NGLL2*iface;

    int k = NGLLX-1;
    int j = (igll/NGLLX);
    int i = (igll-j*NGLLX);

    int iglob = ibool[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] - 1 ;

    realw eta = ( noise_surface_movie[INDEX3(NDIM,NGLL2,0,igll,iface)]*normal_x_noise[ipoin]+
                 noise_surface_movie[INDEX3(NDIM,NGLL2,1,igll,iface)]*normal_y_noise[ipoin]+
                 noise_surface_movie[INDEX3(NDIM,NGLL2,2,igll,iface)]*normal_z_noise[ipoin]);

    Sigma_kl[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] += deltat*eta*
                                                      (normal_x_noise[ipoin]*displ[3*iglob]+
                                                       normal_y_noise[ipoin]*displ[1+3*iglob]+
                                                       normal_z_noise[ipoin]*displ[2+3*iglob]);
  }
}
