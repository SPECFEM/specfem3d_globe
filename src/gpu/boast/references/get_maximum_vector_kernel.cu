// from check_fields_cuda.cu
#define BLOCKSIZE_TRANSFER 256

typedef float realw;

__global__ void get_maximum_vector_kernel(realw* array, int size, realw* d_max){

  // reduction example:
  __shared__ realw sdata[BLOCKSIZE_TRANSFER] ;

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int bx = blockIdx.y*gridDim.x+blockIdx.x;
  unsigned int i = tid + bx*blockDim.x;

  // loads values into shared memory: assume array is a vector array
  sdata[tid] = (i < size) ? sqrt(array[i*3]*array[i*3]
                                 + array[i*3+1]*array[i*3+1]
                                 + array[i*3+2]*array[i*3+2]) : 0.0f ;

  __syncthreads();

  // do reduction in shared mem
  for(unsigned int s=blockDim.x/2; s>0; s>>=1)
  {
    if (tid < s){
      // summation:
      //sdata[tid] += sdata[tid + s];
      // maximum:
      if( sdata[tid] < sdata[tid + s] ) sdata[tid] = sdata[tid + s];
    }
    __syncthreads();
  }

  // write result for this block to global mem
  if (tid == 0) d_max[bx] = sdata[0];

}
