//resort_array.cu
//by boast version

__global__ void resort_array(float * old_array, const int NSPEC){
  int ispec;
  int i,id,idx,t_idx,tx,offset;

  __shared__ float sh_tmp[(2625)];

  ispec = blockIdx.x + (blockIdx.y) * (gridDim.x);
  if (ispec < NSPEC) {
    tx = threadIdx.x;
    offset = ((ispec) * (NGLL3)) * (21) + tx;
    for (i = 0; i <= 20; i += 1) {
      sh_tmp[(i) * (NGLL3) + tx] = old_array[(i) * (NGLL3) + offset];
    }
  }

  __syncthreads();

  if (ispec < NSPEC) {
    for (i = 0; i <= 20; i += 1) {
      id = (i) * (NGLL3) + tx;
      idx = (id) / (21);
      t_idx = id % 21;
      old_array[(i) * (NGLL3) + offset] = sh_tmp[idx + (t_idx) * (NGLL3)];
    }
  }
}
