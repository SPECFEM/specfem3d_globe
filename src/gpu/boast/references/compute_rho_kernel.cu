// from compute_kernels_cuda.cu
#define NGLL3 125

typedef float realw;

__global__ void compute_rho_kernel(int* ibool,
                                   realw* accel,
                                   realw* b_displ,
                                   realw* rho_kl,
                                   int NSPEC,
                                   realw deltat) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;

  // handles case when there is 1 extra block (due to rectangular grid)
  if(ispec < NSPEC) {

    int ijk_ispec = threadIdx.x + NGLL3*ispec;
    int iglob = ibool[ijk_ispec] - 1 ;

    // density kernel
    rho_kl[ijk_ispec] += deltat * (accel[3*iglob]*b_displ[3*iglob]+
                                   accel[3*iglob+1]*b_displ[3*iglob+1]+
                                   accel[3*iglob+2]*b_displ[3*iglob+2]);
  }
}
