// from noise_tomography_cuda.cu
#define NDIM 3
#define NGLLX 5
#define NGLL2 25
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))

typedef float realw;

__global__ void noise_transfer_surface_to_host_kernel(int* ibelm_top,
                                                      int nspec_top,
                                                      int* ibool,
                                                      realw* displ,
                                                      realw* noise_surface_movie) {
  int igll = threadIdx.x;
  int iface = blockIdx.x + blockIdx.y*gridDim.x;

  // int id = tx + blockIdx.x*blockDim.x + blockIdx.y*blockDim.x*gridDim.x;

  if(iface < nspec_top) {
    int ispec = ibelm_top[iface]-1; //-1 for C-based indexing

    int k = NGLLX-1;
    int j = (igll/NGLLX);
    int i = (igll-j*NGLLX);

    int iglob = ibool[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

    noise_surface_movie[INDEX3(NDIM,NGLL2,0,igll,iface)] = displ[iglob*3];
    noise_surface_movie[INDEX3(NDIM,NGLL2,1,igll,iface)] = displ[iglob*3+1];
    noise_surface_movie[INDEX3(NDIM,NGLL2,2,igll,iface)] = displ[iglob*3+2];
  }
}
