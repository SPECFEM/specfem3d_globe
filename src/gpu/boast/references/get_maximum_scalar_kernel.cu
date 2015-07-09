// from check_fields_cuda.cu
#define BLOCKSIZE_TRANSFER 256

typedef float realw;

__global__ void get_maximum_scalar_kernel(realw* array, int size, realw* d_max){

  /* simplest version: uses only 1 thread
   realw max;
   max = 0;
   // finds maximum value in array
   if (size > 0){
   max = fabs(array[0]);
   for( int i=1; i < size; i++){
   if (fabs(array[i]) > max) max = fabs(array[i]);
   }
   }
   *d_max = max;
   */

  // reduction example:
  __shared__ realw sdata[BLOCKSIZE_TRANSFER] ;

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int bx = blockIdx.y*gridDim.x+blockIdx.x;
  unsigned int i = tid + bx*blockDim.x;

  // loads absolute values into shared memory
  sdata[tid] = (i < size) ? fabs(array[i]) : 0.0f ;

  __syncthreads();

  // do reduction in shared mem
  for(unsigned int s=blockDim.x/2; s>0; s>>=1)
  {
    if (tid < s){
      // summation:
      //sdata[tid] += sdata[tid + s];
      // maximum:
      if (sdata[tid] < sdata[tid + s]) sdata[tid] = sdata[tid + s];
    }
    __syncthreads();
  }

  // write result for this block to global mem
  if (tid == 0) d_max[bx] = sdata[0];

}
