// from compute_add_sources_elastic_cuda.cu
#define NDIM 3
#define NGLLX 5
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))
#define INDEX5(xsize,ysize,zsize,isize,x,y,z,i,j) x + xsize*(y + ysize*(z + zsize*(i + isize*(j))))

typedef float realw;

__global__ void compute_add_sources_kernel(realw* accel,
                                           int* ibool,
                                           realw* sourcearrays,
                                           double* stf_pre_compute,
                                           int myrank,
                                           int* islice_selected_source,
                                           int* ispec_selected_source,
                                           int NSOURCES) {
  int ispec,iglob;
  realw stf;

  int i = threadIdx.x;
  int j = threadIdx.y;
  int k = threadIdx.z;
  int isource  = blockIdx.x + gridDim.x*blockIdx.y; // bx

  // when NSOURCES > MAXIMUM_GRID_DIM, but mod(nspec_top,2) > 0, we end up with an extra block.
  if(isource < NSOURCES) {
    if(myrank == islice_selected_source[isource]) {

      ispec = ispec_selected_source[isource]-1;

      stf = (realw) stf_pre_compute[isource];
      iglob = ibool[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

      // note: for global version, sourcearrays has dimensions
      //            sourcearrays(NDIM,NGLLX,NGLLY,NGLLZ,NSOURCES)
      atomicAdd(&accel[3*iglob], sourcearrays[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,0,i,j,k,isource)]*stf);
      atomicAdd(&accel[3*iglob+1], sourcearrays[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,1,i,j,k,isource)]*stf);
      atomicAdd(&accel[3*iglob+2], sourcearrays[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,2,i,j,k,isource)]*stf);
    }
  }
}
