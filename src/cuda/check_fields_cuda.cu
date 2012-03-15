/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  2 . 0
 !               ---------------------------------------
 !
 !          Main authors: Dimitri Komatitsch and Jeroen Tromp
 !    Princeton University, USA and University of Pau / CNRS / INRIA
 ! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
 !                            April 2011
 !
 ! This program is free software; you can redistribute it and/or modify
 ! it under the terms of the GNU General Public License as published by
 ! the Free Software Foundation; either version 2 of the License, or
 ! (at your option) any later version.
 !
 ! This program is distributed in the hope that it will be useful,
 ! but WITHOUT ANY WARRANTY; without even the implied warranty of
 ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ! GNU General Public License for more details.
 !
 ! You should have received a copy of the GNU General Public License along
 ! with this program; if not, write to the Free Software Foundation, Inc.,
 ! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 !
 !=====================================================================
 */

#include <stdio.h>
#include <cuda.h>
#include <cublas.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"
#include "prepare_constants_cuda.h"


/* ----------------------------------------------------------------------------------------------- */

// Helper functions

/* ----------------------------------------------------------------------------------------------- */

double get_time()
{
  struct timeval t;
  struct timezone tzp;
  gettimeofday(&t, &tzp);
  return t.tv_sec + t.tv_usec*1e-6;
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(pause_for_debug,
              PAUSE_FOR_DEBUG)() {
  TRACE("pause_for_debug");

  pause_for_debugger(1);
}

/* ----------------------------------------------------------------------------------------------- */

void pause_for_debugger(int pause) {
  if(pause) {
    int myrank;
#ifdef WITH_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#else
    myrank = 0;
#endif
    printf("I'm rank %d\n",myrank);
    int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s:%d ready for attach\n", getpid(), hostname,myrank);
    FILE *file = fopen("/scratch/eiger/rietmann/attach_gdb.txt","w+");
    if( file != NULL ){
      fprintf(file,"PID %d on %s:%d ready for attach\n", getpid(), hostname,myrank);
      fclose(file);
    }
    fflush(stdout);
    while (0 == i)
      sleep(5);
  }
}

/* ----------------------------------------------------------------------------------------------- */

void exit_on_cuda_error(char* kernel_name) {
  // sync and check to catch errors from previous async operations
  cudaThreadSynchronize();
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    fprintf(stderr,"Error after: %s --- %s\n", kernel_name, cudaGetErrorString(err));
    pause_for_debugger(0);
    //free(kernel_name);
    exit(EXIT_FAILURE);
  }
}

/* ----------------------------------------------------------------------------------------------- */

void exit_on_error(char* info)
{
  printf("\nERROR: %s\n",info);
  fflush(stdout);
#ifdef WITH_MPI
  MPI_Abort(MPI_COMM_WORLD,1);
#endif
  //free(info);
  exit(EXIT_FAILURE);
  return;
}

/* ----------------------------------------------------------------------------------------------- */

void print_CUDA_error_if_any(cudaError_t err, int num)
{
  if (cudaSuccess != err)
  {
    printf("\nCUDA error !!!!! <%s> !!!!! \nat CUDA call error code: # %d\n",cudaGetErrorString(err),num);
    fflush(stdout);
#ifdef WITH_MPI
    MPI_Abort(MPI_COMM_WORLD,1);
#endif
    exit(EXIT_FAILURE);
  }
  return;
}

/* ----------------------------------------------------------------------------------------------- */

void get_free_memory(double* free_db, double* used_db, double* total_db) {

  // gets memory usage in byte
  size_t free_byte ;
  size_t total_byte ;
  cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;
  if ( cudaSuccess != cuda_status ){
    printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
    exit(EXIT_FAILURE);
  }

  *free_db = (double)free_byte ;
  *total_db = (double)total_byte ;
  *used_db = *total_db - *free_db ;
  return;
}

/* ----------------------------------------------------------------------------------------------- */

// Saves GPU memory usage to file
void output_free_memory(int myrank,char* info_str) {

  FILE* fp;
  char filename[BUFSIZ];
  double free_db,used_db,total_db;

  get_free_memory(&free_db,&used_db,&total_db);

  sprintf(filename,"OUTPUT_FILES/gpu_device_mem_usage_proc_%06d.txt",myrank);
  fp = fopen(filename,"a+");
  if( fp != NULL ){
    fprintf(fp,"%d: @%s GPU memory usage: used = %f MB, free = %f MB, total = %f MB\n", myrank, info_str,
            used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);
    fclose(fp);
  }
}

/* ----------------------------------------------------------------------------------------------- */

// Fortran-callable version of above method
extern "C"
void FC_FUNC_(output_free_device_memory,
              OUTPUT_FREE_DEVICE_MEMORY)(int* myrank) {

  TRACE("output_free_device_memory");

  char info[6];
  sprintf(info,"f %d:",*myrank);
  output_free_memory(*myrank,info);
}

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(get_free_device_memory,
              get_FREE_DEVICE_MEMORY)(realw* free, realw* used, realw* total ) {
  TRACE("get_free_device_memory");

  double free_db,used_db,total_db;

  get_free_memory(&free_db,&used_db,&total_db);

  // converts to MB
  *free = (realw) free_db/1024.0/1024.0;
  *used = (realw) used_db/1024.0/1024.0;
  *total = (realw) total_db/1024.0/1024.0;
  return;
}


/* ----------------------------------------------------------------------------------------------- */
//daniel: helper function
/*
 __global__ void check_phase_ispec_kernel(int num_phase_ispec,
 int* phase_ispec,
 int NSPEC_AB,
 int* ier) {

 int i,ispec,iphase,count0,count1;
 *ier = 0;

 for(iphase=0; iphase < 2; iphase++){
 count0 = 0;
 count1 = 0;

 for(i=0; i < num_phase_ispec; i++){
 ispec = phase_ispec[iphase*num_phase_ispec + i] - 1;
 if( ispec < -1 || ispec >= NSPEC_AB ){
 printf("Error in d_phase_ispec_inner_elastic %d %d\n",i,ispec);
 *ier = 1;
 return;
 }
 if( ispec >= 0 ){ count0++;}
 if( ispec < 0 ){ count1++;}
 }

 printf("check_phase_ispec done: phase %d, count = %d %d \n",iphase,count0,count1);

 }
 }

 void check_phase_ispec(long* Mesh_pointer_f,int type){

 Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

 printf("check phase_ispec for type=%d\n",type);

 dim3 grid(1,1);
 dim3 threads(1,1,1);

 int* h_debug = (int*) calloc(1,sizeof(int));
 int* d_debug;
 cudaMalloc((void**)&d_debug,sizeof(int));

 if( type == 1 ){
 check_phase_ispec_kernel<<<grid,threads>>>(mp->num_phase_ispec_elastic,
 mp->d_phase_ispec_inner_elastic,
 mp->NSPEC_AB,
 d_debug);
 }else if( type == 2 ){
 check_phase_ispec_kernel<<<grid,threads>>>(mp->num_phase_ispec_acoustic,
 mp->d_phase_ispec_inner_acoustic,
 mp->NSPEC_AB,
 d_debug);
 }

 cudaMemcpy(h_debug,d_debug,1*sizeof(int),cudaMemcpyDeviceToHost);
 cudaFree(d_debug);
 if( *h_debug != 0 ){printf("error for type=%d\n",type); exit(1);}
 free(h_debug);
 fflush(stdout);

 #ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
 exit_on_cuda_error("check_phase_ispec");
 #endif

 }
 */

/* ----------------------------------------------------------------------------------------------- */
//daniel: helper function
/*
 __global__ void check_ispec_is_kernel(int NSPEC_AB,
 int* ispec_is,
 int* ier) {

 int ispec,count0,count1;

 *ier = 0;
 count0 = 0;
 count1 = 0;
 for(ispec=0; ispec < NSPEC_AB; ispec++){
 if( ispec_is[ispec] < -1 || ispec_is[ispec] > 1 ){
 printf("Error in ispec_is %d %d\n",ispec,ispec_is[ispec]);
 *ier = 1;
 return;
 //exit(1);
 }
 if( ispec_is[ispec] == 0 ){count0++;}
 if( ispec_is[ispec] != 0 ){count1++;}
 }
 printf("check_ispec_is done: count = %d %d\n",count0,count1);
 }

 void check_ispec_is(long* Mesh_pointer_f,int type){

 Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

 printf("check ispec_is for type=%d\n",type);

 dim3 grid(1,1);
 dim3 threads(1,1,1);

 int* h_debug = (int*) calloc(1,sizeof(int));
 int* d_debug;
 cudaMalloc((void**)&d_debug,sizeof(int));

 if( type == 0 ){
 check_ispec_is_kernel<<<grid,threads>>>(mp->NSPEC_AB,
 mp->d_ispec_is_inner,
 d_debug);
 }else if( type == 1 ){
 check_ispec_is_kernel<<<grid,threads>>>(mp->NSPEC_AB,
 mp->d_ispec_is_elastic,
 d_debug);
 }else if( type == 2 ){
 check_ispec_is_kernel<<<grid,threads>>>(mp->NSPEC_AB,
 mp->d_ispec_is_acoustic,
 d_debug);
 }

 cudaMemcpy(h_debug,d_debug,1*sizeof(int),cudaMemcpyDeviceToHost);
 cudaFree(d_debug);
 if( *h_debug != 0 ){printf("error for type=%d\n",type); exit(1);}
 free(h_debug);
 fflush(stdout);

 #ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
 exit_on_cuda_error("check_ispec_is");
 #endif
 }
 */
/* ----------------------------------------------------------------------------------------------- */
//daniel: helper function
/*
 __global__ void check_array_ispec_kernel(int num_array_ispec,
 int* array_ispec,
 int NSPEC_AB,
 int* ier) {

 int i,ispec,count0,count1;

 *ier = 0;
 count0 = 0;
 count1 = 0;

 for(i=0; i < num_array_ispec; i++){
 ispec = array_ispec[i] - 1;
 if( ispec < -1 || ispec >= NSPEC_AB ){
 printf("Error in d_array_ispec %d %d\n",i,ispec);
 *ier = 1;
 return;
 }
 if( ispec >= 0 ){ count0++;}
 if( ispec < 0 ){ count1++;}
 }

 printf("check_array_ispec done: count = %d %d \n",count0,count1);
 }

 void check_array_ispec(long* Mesh_pointer_f,int type){

 Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

 printf("check array_ispec for type=%d\n",type);

 dim3 grid(1,1);
 dim3 threads(1,1,1);

 int* h_debug = (int*) calloc(1,sizeof(int));
 int* d_debug;
 cudaMalloc((void**)&d_debug,sizeof(int));

 if( type == 1 ){
 check_array_ispec_kernel<<<grid,threads>>>(mp->d_num_abs_boundary_faces,
 mp->d_abs_boundary_ispec,
 mp->NSPEC_AB,
 d_debug);
 }

 cudaMemcpy(h_debug,d_debug,1*sizeof(int),cudaMemcpyDeviceToHost);
 cudaFree(d_debug);
 if( *h_debug != 0 ){printf("error for type=%d\n",type); exit(1);}
 free(h_debug);
 fflush(stdout);

 #ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
 exit_on_cuda_error("check_array_ispec");
 #endif

 }
 */


/* ----------------------------------------------------------------------------------------------- */

// Check functions

/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(check_max_norm_displ_gpu,
              CHECK_MAX_NORM_DISPL_GPU)(int* size, realw* displ,long* Mesh_pointer_f,int* announceID) {

TRACE("check_max_norm_displ_gpu");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  cudaMemcpy(displ, mp->d_displ,*size*sizeof(realw),cudaMemcpyDeviceToHost);
  realw maxnorm=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(displ[i]));
  }
  printf("%d: maxnorm of forward displ = %e\n",*announceID,maxnorm);
}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(check_max_norm_vector,
              CHECK_MAX_NORM_VECTOR)(int* size, realw* vector1, int* announceID) {

TRACE("check_max_norm_vector");

  int procid;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
#else
  procid = 0;
#endif
  realw maxnorm=0;
  int maxloc;
  for(int i=0;i<*size;i++) {
    if(maxnorm<fabsf(vector1[i])) {
      maxnorm = vector1[i];
      maxloc = i;
    }
  }
  printf("%d:maxnorm of vector %d [%d] = %e\n",procid,*announceID,maxloc,maxnorm);
}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(check_max_norm_displ,
              CHECK_MAX_NORM_DISPL)(int* size, realw* displ, int* announceID) {

TRACE("check_max_norm_displ");

  realw maxnorm=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(displ[i]));
  }
  printf("%d: maxnorm of forward displ = %e\n",*announceID,maxnorm);
}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(check_max_norm_b_displ_gpu,
              CHECK_MAX_NORM_B_DISPL_GPU)(int* size, realw* b_displ,long* Mesh_pointer_f,int* announceID) {

TRACE("check_max_norm_b_displ_gpu");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  realw* b_accel = (realw*)malloc(*size*sizeof(realw));

  cudaMemcpy(b_displ, mp->d_b_displ,*size*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(b_accel, mp->d_b_accel,*size*sizeof(realw),cudaMemcpyDeviceToHost);

  realw maxnorm=0;
  realw maxnorm_accel=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(b_displ[i]));
    maxnorm_accel = MAX(maxnorm,fabsf(b_accel[i]));
  }
  free(b_accel);
  printf("%d: maxnorm of backward displ = %e\n",*announceID,maxnorm);
  printf("%d: maxnorm of backward accel = %e\n",*announceID,maxnorm_accel);
}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(check_max_norm_b_accel_gpu,
              CHECK_MAX_NORM_B_ACCEL_GPU)(int* size, realw* b_accel,long* Mesh_pointer_f,int* announceID) {

TRACE("check_max_norm_b_accel_gpu");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  cudaMemcpy(b_accel, mp->d_b_accel,*size*sizeof(realw),cudaMemcpyDeviceToHost);

  realw maxnorm=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(b_accel[i]));
  }
  printf("%d: maxnorm of backward accel = %e\n",*announceID,maxnorm);
}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(check_max_norm_b_veloc_gpu,
              CHECK_MAX_NORM_B_VELOC_GPU)(int* size, realw* b_veloc,long* Mesh_pointer_f,int* announceID) {

TRACE("check_max_norm_b_veloc_gpu");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  cudaMemcpy(b_veloc, mp->d_b_veloc,*size*sizeof(realw),cudaMemcpyDeviceToHost);

  realw maxnorm=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(b_veloc[i]));
  }
  printf("%d: maxnorm of backward veloc = %e\n",*announceID,maxnorm);
}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(check_max_norm_b_displ,
              CHECK_MAX_NORM_B_DISPL)(int* size, realw* b_displ,int* announceID) {

TRACE("check_max_norm_b_displ");

  realw maxnorm=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(b_displ[i]));
  }
  printf("%d:maxnorm of backward displ = %e\n",*announceID,maxnorm);
}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(check_max_norm_b_accel,
              CHECK_MAX_NORM_B_ACCEL)(int* size, realw* b_accel,int* announceID) {

TRACE("check_max_norm_b_accel");

  realw maxnorm=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(b_accel[i]));
  }
  printf("%d:maxnorm of backward accel = %e\n",*announceID,maxnorm);
}
*/
/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(check_error_vectors,
              CHECK_ERROR_VECTORS)(int* sizef, realw* vector1,realw* vector2) {

TRACE("check_error_vectors");

  int size = *sizef;

  double diff2 = 0;
  double sum = 0;
  double temp;
  double maxerr=0;
  int maxerrorloc;

  for(int i=0;i<size;++i) {
    temp = vector1[i]-vector2[i];
    diff2 += temp*temp;
    sum += vector1[i]*vector1[i];
    if(maxerr < fabsf(temp)) {
      maxerr = abs(temp);
      maxerrorloc = i;
    }
  }

  printf("rel error = %f, maxerr = %e @ %d\n",diff2/sum,maxerr,maxerrorloc);
  int myrank;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#else
  myrank = 0;
#endif
  if(myrank == 0) {
    for(int i=maxerrorloc;i>maxerrorloc-5;i--) {
      printf("[%d]: %e vs. %e\n",i,vector1[i],vector2[i]);
    }
  }

}
*/

/* ----------------------------------------------------------------------------------------------- */

// Auxiliary functions

/* ----------------------------------------------------------------------------------------------- */


/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(get_max_accel,
              GET_MAX_ACCEL)(int* itf,int* sizef,long* Mesh_pointer) {

TRACE("get_max_accel");

  Mesh* mp = (Mesh*)(*Mesh_pointer);
  int procid;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
#else
  procid = 0;
#endif
  int size = *sizef;
  int it = *itf;
  realw* accel_cpy = (realw*)malloc(size*sizeof(realw));
  cudaMemcpy(accel_cpy,mp->d_accel,size*sizeof(realw),cudaMemcpyDeviceToHost);
  realw maxval=0;
  for(int i=0;i<size;++i) {
    maxval = MAX(maxval,accel_cpy[i]);
  }
  printf("%d/%d: max=%e\n",it,procid,maxval);
  free(accel_cpy);
}
*/
/* ----------------------------------------------------------------------------------------------- */

// ACOUSTIC simulations

/* ----------------------------------------------------------------------------------------------- */

__global__ void get_maximum_scalar_kernel(realw* array, int size, realw* d_max){

  /* simplest version: uses only 1 thread
   realw max;
   max = 0;
   // finds maximum value in array
   if( size > 0 ){
   max = abs(array[0]);
   for( int i=1; i < size; i++){
   if( abs(array[i]) > max ) max = abs(array[i]);
   }
   }
   *d_max = max;
   */

  // reduction example:
  __shared__ realw sdata[256] ;

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // loads absolute values into shared memory
  sdata[tid] = (i < size) ? fabs(array[i]) : 0.0 ;

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
  if (tid == 0) d_max[blockIdx.x] = sdata[0];

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(check_norm_acoustic_from_device,
              CHECK_NORM_ACOUSTIC_FROM_DEVICE)(realw* norm,
                                                  long* Mesh_pointer_f,
                                                  int* SIMULATION_TYPE) {

TRACE("check_norm_acoustic_from_device");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container
  realw max;
  realw *d_max;

  max = 0;

  /* way 1 : timing Elapsed time: 8.464813e-03
   realw* h_array;
   h_array = (realw*)calloc(mp->NGLOB_AB,sizeof(realw));

   print_CUDA_error_if_any(cudaMemcpy(h_array,mp->d_potential_dot_dot_acoustic,
   sizeof(realw)*(mp->NGLOB_AB),cudaMemcpyDeviceToHost),131);

   // finds maximum value in array
   max = h_array[0];
   for( int i=1; i < mp->NGLOB_AB; i++){
   if( abs(h_array[i]) > max ) max = abs(h_array[i]);
   }
   free(h_array);
   */

  /* way 2: timing Elapsed time: 8.818102e-02
   // launch simple kernel
   cudaMalloc((void**)&d_max,sizeof(realw));

   dim3 grid(1,1);
   dim3 threads(1,1,1);

   get_maximum_kernel<<<grid,threads>>>(mp->d_potential_dot_dot_acoustic,
   mp->NGLOB_AB,
   d_max);
   print_CUDA_error_if_any(cudaMemcpy(&max,d_max, sizeof(realw), cudaMemcpyDeviceToHost),222);

   cudaFree(d_max);
   */

  // way 2 b: timing Elapsed time: 1.236916e-03
  // launch simple reduction kernel
  realw* h_max;
  int blocksize = 256;

  // outer core
  int size = mp->NGLOB_OUTER_CORE;

  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;
  int num_blocks_x = size_padded/blocksize;
  int num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  h_max = (realw*) calloc(num_blocks_x*num_blocks_y,sizeof(realw));
  cudaMalloc((void**)&d_max,num_blocks_x*num_blocks_y*sizeof(realw));

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  if(*SIMULATION_TYPE == 1 ){
    get_maximum_scalar_kernel<<<grid,threads>>>(mp->d_displ_outer_core,size,d_max);
  }else if(*SIMULATION_TYPE == 3 ){
    get_maximum_scalar_kernel<<<grid,threads>>>(mp->d_b_displ_outer_core,size,d_max);
  }

  print_CUDA_error_if_any(cudaMemcpy(h_max,d_max,num_blocks_x*num_blocks_y*sizeof(realw),
                                     cudaMemcpyDeviceToHost),222);

  // determines max for all blocks
  max = h_max[0];
  for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
    if( max < h_max[i]) max = h_max[i];
  }

  cudaFree(d_max);
  free(h_max);

  /* way 3: doesn't work properly...
   cublasStatus status;

   // Initialize CUBLAS
   status = cublasInit();
   if (status != CUBLAS_STATUS_SUCCESS) {
   fprintf (stderr, "!!!! CUBLAS initialization error\n");
   exit(1);
   }

   // cublas function: cublasIsamax
   //       finds the smallest index of the maximum magnitude element of single 
   //      precision vector x
   int incr = 1;
   int imax = 0;
   imax = cublasIsamax(mp->NGLOB_AB,(realw*)mp->d_potential_dot_dot_acoustic, incr);
   status= cublasGetError();
   if (status != CUBLAS_STATUS_SUCCESS) {
   fprintf (stderr, "!!!! CUBLAS error in cublasIsamax\n");
   exit(1);
   }

   print_CUDA_error_if_any(cudaMemcpy(&max,&(mp->d_potential_dot_dot_acoustic[imax]),
                      sizeof(realw), cudaMemcpyDeviceToHost),222);

   printf("maximum %i %i %f \n",mp->NGLOB_AB,imax,max);

   // Shutdown
   status = cublasShutdown();
   if (status != CUBLAS_STATUS_SUCCESS) {
   fprintf (stderr, "!!!! shutdown error (A)\n");
   exit(1);
   }

   */

  // return result
  *norm = max;

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("after check_norm_acoustic_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// ELASTIC simulations

/* ----------------------------------------------------------------------------------------------- */

__global__ void get_maximum_vector_kernel(realw* array, int size, realw* d_max){

  // reduction example:
  __shared__ realw sdata[256] ;

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  // loads values into shared memory: assume array is a vector array
  sdata[tid] = (i < size) ? sqrt(array[i*3]*array[i*3]
                                 + array[i*3+1]*array[i*3+1]
                                 + array[i*3+2]*array[i*3+2]) : 0.0 ;

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
  if (tid == 0) d_max[blockIdx.x] = sdata[0];

}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(check_norm_elastic_from_device,
              CHECK_NORM_ELASTIC_FROM_DEVICE)(realw* norm,
                                              long* Mesh_pointer_f,
                                              int* SIMULATION_TYPE) {

  TRACE("check_norm_elastic_from_device");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  realw max,max_crust_mantle,max_inner_core;
  realw *d_max;
  int num_blocks_x,num_blocks_y;
  int size,size_padded;
  dim3 grid,threads;

  // launch simple reduction kernel
  realw* h_max;
  int blocksize = 256;

  // crust_mantle
  max = 0;
  size = mp->NGLOB_CRUST_MANTLE;

  size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;
  num_blocks_x = size_padded/blocksize;
  num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  h_max = (realw*) calloc(num_blocks_x*num_blocks_y,sizeof(realw));
  cudaMalloc((void**)&d_max,num_blocks_x*num_blocks_y*sizeof(realw));

  grid = dim3(num_blocks_x,num_blocks_y);
  threads = dim3(blocksize,1,1);

  if(*SIMULATION_TYPE == 1 ){
    get_maximum_vector_kernel<<<grid,threads>>>(mp->d_displ_crust_mantle,size,d_max);
  }else if(*SIMULATION_TYPE == 3 ){
    get_maximum_vector_kernel<<<grid,threads>>>(mp->d_b_displ_crust_mantle,size,d_max);
  }

  // copies to CPU
  print_CUDA_error_if_any(cudaMemcpy(h_max,d_max,num_blocks_x*num_blocks_y*sizeof(realw),
                                     cudaMemcpyDeviceToHost),222);

  // determines max for all blocks
  max = h_max[0];
  for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
    if( max < h_max[i]) max = h_max[i];
  }
  max_crust_mantle = max;

  cudaFree(d_max);
  free(h_max);

  // inner_core
  max = 0;
  size = mp->NGLOB_INNER_CORE;

  size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;
  num_blocks_x = size_padded/blocksize;
  num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  h_max = (realw*) calloc(num_blocks_x*num_blocks_y,sizeof(realw));
  cudaMalloc((void**)&d_max,num_blocks_x*num_blocks_y*sizeof(realw));

  grid = dim3(num_blocks_x,num_blocks_y);
  threads = dim3(blocksize,1,1);

  if(*SIMULATION_TYPE == 1 ){
    get_maximum_vector_kernel<<<grid,threads>>>(mp->d_displ_inner_core,size,d_max);
  }else if(*SIMULATION_TYPE == 3 ){
    get_maximum_vector_kernel<<<grid,threads>>>(mp->d_b_displ_inner_core,size,d_max);
  }

  // copies to CPU
  print_CUDA_error_if_any(cudaMemcpy(h_max,d_max,num_blocks_x*num_blocks_y*sizeof(realw),
                                     cudaMemcpyDeviceToHost),222);

  // determines max for all blocks
  max = h_max[0];
  for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
    if( max < h_max[i]) max = h_max[i];
  }
  max_inner_core = max;

  cudaFree(d_max);
  free(h_max);

  // return result
  max = MAX(max_inner_core,max_crust_mantle);
  *norm = max;

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("after check_norm_elastic_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(check_norm_strain_from_device,
              CHECK_NORM_STRAIN_FROM_DEVICE)(realw* strain_norm,
                                             realw* strain_norm2,
                                             long* Mesh_pointer_f) {

  TRACE("check_norm_strain_from_device");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  realw max,max_eps;
  realw *d_max;
  int num_blocks_x,num_blocks_y;
  int size,size_padded;
  dim3 grid;
  dim3 threads;

  // launch simple reduction kernel
  realw* h_max;
  int blocksize = 256;

  // crust_mantle strain arrays
  size = NGLL3*(mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY);

  size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;
  num_blocks_x = size_padded/blocksize;
  num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  h_max = (realw*) calloc(num_blocks_x*num_blocks_y,sizeof(realw));
  cudaMalloc((void**)&d_max,num_blocks_x*num_blocks_y*sizeof(realw));

  grid = dim3(num_blocks_x,num_blocks_y);
  threads = dim3(blocksize,1,1);

  max = 0.0f;

  // determines max for: eps_trace_over_3_crust_mantle
  get_maximum_scalar_kernel<<<grid,threads>>>(mp->d_eps_trace_over_3_crust_mantle,size,d_max);

  print_CUDA_error_if_any(cudaMemcpy(h_max,d_max,num_blocks_x*num_blocks_y*sizeof(realw),
                                     cudaMemcpyDeviceToHost),221);
  max = h_max[0];
  for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
    if( max < h_max[i]) max = h_max[i];
  }
  // strain trace maximum
  *strain_norm = max;

  // frees arrays
  cudaFree(d_max);
  free(h_max);

  // initializes
  // crust_mantle arrays
  size = NGLL3*(mp->NSPEC_CRUST_MANTLE);

  size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;
  num_blocks_x = size_padded/blocksize;
  num_blocks_y = 1;
  while(num_blocks_x > 65535) {
    num_blocks_x = (int) ceil(num_blocks_x*0.5f);
    num_blocks_y = num_blocks_y*2;
  }

  h_max = (realw*) calloc(num_blocks_x*num_blocks_y,sizeof(realw));
  cudaMalloc((void**)&d_max,num_blocks_x*num_blocks_y*sizeof(realw));

  grid = dim3(num_blocks_x,num_blocks_y);
  threads = dim3(blocksize,1,1);

  max_eps = 0.0f;

  // determines max for: epsilondev_xx_crust_mantle
  get_maximum_scalar_kernel<<<grid,threads>>>(mp->d_epsilondev_xx_crust_mantle,size,d_max);

  print_CUDA_error_if_any(cudaMemcpy(h_max,d_max,num_blocks_x*num_blocks_y*sizeof(realw),
                                     cudaMemcpyDeviceToHost),222);
  max = h_max[0];
  for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
    if( max < h_max[i]) max = h_max[i];
  }
  max_eps = MAX(max_eps,max);

  // determines max for: epsilondev_yy_crust_mantle
  get_maximum_scalar_kernel<<<grid,threads>>>(mp->d_epsilondev_yy_crust_mantle,size,d_max);

  print_CUDA_error_if_any(cudaMemcpy(h_max,d_max,num_blocks_x*num_blocks_y*sizeof(realw),
                                     cudaMemcpyDeviceToHost),223);
  max = h_max[0];
  for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
    if( max < h_max[i]) max = h_max[i];
  }
  max_eps = MAX(max_eps,max);

  // determines max for: epsilondev_xy_crust_mantle
  get_maximum_scalar_kernel<<<grid,threads>>>(mp->d_epsilondev_xy_crust_mantle,size,d_max);

  print_CUDA_error_if_any(cudaMemcpy(h_max,d_max,num_blocks_x*num_blocks_y*sizeof(realw),
                                     cudaMemcpyDeviceToHost),224);
  max = h_max[0];
  for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
    if( max < h_max[i]) max = h_max[i];
  }
  max_eps = MAX(max_eps,max);

  // determines max for: epsilondev_xz_crust_mantle
  get_maximum_scalar_kernel<<<grid,threads>>>(mp->d_epsilondev_xz_crust_mantle,size,d_max);

  print_CUDA_error_if_any(cudaMemcpy(h_max,d_max,num_blocks_x*num_blocks_y*sizeof(realw),
                                     cudaMemcpyDeviceToHost),225);
  max = h_max[0];
  for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
    if( max < h_max[i]) max = h_max[i];
  }
  max_eps = MAX(max_eps,max);

  // determines max for: epsilondev_yz_crust_mantle
  get_maximum_scalar_kernel<<<grid,threads>>>(mp->d_epsilondev_yz_crust_mantle,size,d_max);

  print_CUDA_error_if_any(cudaMemcpy(h_max,d_max,num_blocks_x*num_blocks_y*sizeof(realw),
                                     cudaMemcpyDeviceToHost),226);
  max = h_max[0];
  for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
    if( max < h_max[i]) max = h_max[i];
  }
  max_eps = MAX(max_eps,max);

  // strain maximum
  *strain_norm2 = max_eps;

  // frees arrays
  cudaFree(d_max);
  free(h_max);

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("after check_norm_strain_from_device");
#endif
}


