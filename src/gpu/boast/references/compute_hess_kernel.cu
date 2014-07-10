// from compute_kernels_cuda.cu
#define NGLL3 125

typedef float realw;

__global__ void compute_hess_kernel(int* ibool,
                                    realw* accel,
                                    realw* b_accel,
                                    realw* hess_kl,
                                    realw deltat,
                                    int NSPEC_AB) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;

  // handles case when there is 1 extra block (due to rectangular grid)
  if(ispec < NSPEC_AB) {

    int ijk = threadIdx.x;
    int ijk_ispec = ijk + NGLL3*ispec;
    int iglob = ibool[ijk_ispec] - 1 ;

    // approximate hessian
    hess_kl[ijk_ispec] += deltat * (accel[3*iglob]*b_accel[3*iglob] +
                                    accel[3*iglob+1]*b_accel[3*iglob+1] +
                                    accel[3*iglob+2]*b_accel[3*iglob+2]);
  }
}
