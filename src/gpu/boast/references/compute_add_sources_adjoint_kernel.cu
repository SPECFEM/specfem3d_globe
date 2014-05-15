// from compute_add_sources_elastic_cuda.cu
#define NDIM 3
#define NGLLX 5
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))
#define INDEX5(xsize,ysize,zsize,isize,x,y,z,i,j) x + xsize*(y + ysize*(z + zsize*(i + isize*(j))))

typedef float realw;

__global__ void compute_add_sources_adjoint_kernel(realw* accel,
                                                   int nrec,
                                                   realw* adj_sourcearrays,
                                                   int* ibool,
                                                   int* ispec_selected_rec,
                                                   int* pre_computed_irec,
                                                   int nadj_rec_local) {

  int ispec,iglob;
  int irec,i,j,k;

  int irec_local = blockIdx.x + gridDim.x*blockIdx.y;

  // when nrec > MAXIMUM_GRID_DIM, but mod(nspec_top,2) > 0, we end up with an extra block.
  if(irec_local < nadj_rec_local) {
    irec = pre_computed_irec[irec_local];
    ispec = ispec_selected_rec[irec]-1;

    i = threadIdx.x;
    j = threadIdx.y;
    k = threadIdx.z;
    iglob = ibool[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

    // atomic operations are absolutely necessary for correctness!
    atomicAdd(&accel[3*iglob], adj_sourcearrays[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,0,i,j,k,irec_local)]);
    atomicAdd(&accel[3*iglob+1], adj_sourcearrays[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,1,i,j,k,irec_local)]);
    atomicAdd(&accel[3*iglob+2], adj_sourcearrays[INDEX5(NDIM,NGLLX,NGLLX,NGLLX,2,i,j,k,irec_local)]);
  }
}
