/*
!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
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

#include "mesh_constants_gpu.h"
#include <string.h>

// GPU initialization

#ifdef USE_OPENCL
// macro definitions used in GPU kernels

#define STR(x) #x
#define PASS(x) {#x, STR(x)}

static struct {
  const char *name;
  const char *value;
} _macro_to_kernel[] = {
  // macro values
  PASS(NDIM),
  PASS(NGLLX), PASS(NGLL2), PASS(NGLL3), PASS(NGLL3_PADDED),
  PASS(N_SLS),
  PASS(IREGION_INNER_CORE),
  PASS(IFLAG_IN_FICTITIOUS_CUBE),
  PASS(COLORING_MIN_NSPEC_OUTER_CORE), PASS(COLORING_MIN_NSPEC_INNER_CORE),
  PASS(R_EARTH_KM),

  // macro functions: not working yet, spaces not allowed in OCL compiler

/* PASS(INDEX2(xsize, x, y)),
   PASS(INDEX3(xsize, ysize, x, y, z)),
   PASS(INDEX4(xsize, ysize, zsize, x, y, z, i)),
   PASS(INDEX4_PADDED(xsize, ysize, zsize, x, y, z, i)),
   PASS(INDEX5(xsize, ysize, zsize, isize, x, y, z, i, j)),
   PASS(INDEX6(xsize, ysize, zsize, isize, jsize, x, y, z, i, j, k)), */

  // macro flags, passed only ifdefed
  PASS(MANUALLY_UNROLLED_LOOPS), PASS(USE_TEXTURES_CONSTANTS), PASS(USE_TEXTURES_FIELDS),

  PASS(USE_LAUNCH_BOUNDS),
  PASS(LAUNCH_MIN_BLOCKS),

  {NULL, NULL}
};
#endif

/* ----------------------------------------------------------------------------------------------- */

// gpu runtime flags
int run_cuda = 0;
int run_opencl = 0;

/* ----------------------------------------------------------------------------------------------- */
// CUDA initialization
/* ----------------------------------------------------------------------------------------------- */

#ifdef USE_CUDA

// initializes CUDA devices

static void initialize_cuda_device(const char *platform_filter, const char *device_filter, int myrank, int *nb_devices) {
  int device_count = 0;
  cudaError_t err;

  // Gets number of GPU devices
  cudaGetDeviceCount(&device_count);
  // Do not check if command failed with `exit_on_cuda_error` since it calls cudaDevice()/ThreadSynchronize():
  // If multiple MPI tasks access multiple GPUs per node, they will try to synchronize
  // GPU 0 and depending on the order of the calls, an error will be raised
  // when setting the device number. If MPS is enabled, some GPUs will silently not be used.
  //
  // being verbose and catches error from first call to CUDA runtime function, without synchronize call
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf(stderr,"Error after cudaGetDeviceCount: %s\n", cudaGetErrorString(err));
    exit_on_error("\
CUDA runtime error: cudaGetDeviceCount failed\n\n\
please check if driver and runtime libraries work together\n\n");
  }

  // returns device count to fortran
  if (device_count == 0) {
    exit_on_error("CUDA runtime error: there is no device supporting CUDA\n");
  }

  *nb_devices = device_count;

  // releases previous contexts
  gpuReset();

  // determines device id for this process
  int *matchingDevices = (int *) malloc (sizeof(int) * device_count);
  int nbMatchingDevices = 0;
  struct cudaDeviceProp deviceProp;
  int i;

  for (i = 0; i < device_count; i++) {
    // get device properties
    cudaGetDeviceProperties(&deviceProp, i);
    if (!strcasestr(deviceProp.name, device_filter)) {
      continue;
    }
    // debug
    //printf("device match: %d match %d out of %d - filter platform = %s device = %s\n",
    //        i,nbMatchingDevices,device_count,platform_filter, device_filter);

    // adds match
    matchingDevices[nbMatchingDevices] = i;
    nbMatchingDevices++;
  }

  if (nbMatchingDevices == 0) {
    printf("Error: no matching devices for criteria %s/%s\n", platform_filter, device_filter);
    exit_on_error("Error CUDA found no matching devices (for device filter set in Par_file)\n");
  }

  int myDevice = matchingDevices[myrank % nbMatchingDevices];

  free(matchingDevices);

  // user error info
  const char* err_info = "\
Please check GPU settings on your node \n\
and/or check CUDA MPS setup to use a single GPU with multiple MPI processes,\n\
e.g., on titan enable environment CRAY_CUDA_MPS=1 to use a single GPU with multiple MPI processes\n\n";

  // sets CUDA device for this process
  // note: setting/getting device ids seems to return success also for multiple processes setting the same GPU id
  //       and even if the GPU mode is thread exclusive (only a single process would be allowed to use a single GPU).
  //       we will have to catch the error later on...
  err = cudaSetDevice(myDevice);
  if (err != cudaSuccess) {
    fprintf(stderr,"Error cudaSetDevice: %s\n", cudaGetErrorString(err));
    if (err == cudaErrorDevicesUnavailable){ fprintf(stderr,"\n%s\n", err_info); }
    exit_on_error("CUDA runtime error: cudaSetDevice failed\n\n");
  }

  // checks if setting device was successful
  int device;
  cudaGetDevice(&device);

  err = cudaGetLastError();
  // debug
  //printf("device set/get: rank %d set %d get %d\n - return %s",myrank,myDevice,device,cudaGetErrorString(err));
  if (err != cudaSuccess) {
    fprintf(stderr,"Error cudaGetDevice: %s\n", cudaGetErrorString(err));
    if (err == cudaErrorDevicesUnavailable){ fprintf(stderr,"\n%s\n", err_info); }
    exit_on_error("CUDA runtime error: cudaGetDevice failed\n\n");
  }

  // checks device id
  if ( device != myDevice){
    fprintf(stderr,"Error cudaGetDevice: setting myDevice = %d is differnt to actual device = %d\n", myDevice,device);
    exit_on_error("Error CUDA setting device failed\n");
  }

  // checks device properties
  cudaGetDeviceProperties(&deviceProp, myDevice);

  // exit if the machine has no CUDA-enabled device
  if (deviceProp.major == 9999 && deviceProp.minor == 9999) {
    fprintf(stderr,"No CUDA-enabled device found, exiting...\n\n");
    exit_on_error("CUDA runtime error: there is no CUDA-enabled device found\n");
  }

  // outputs device info to file
  char filename[BUFSIZ];
  FILE* fp;
  int do_output_info = 0;

  // by default, only master process outputs device info to avoid file cluttering
  if (myrank == 0) {
    do_output_info = 1;
    sprintf(filename, "OUTPUT_FILES/gpu_device_info.txt");
  }
  // debugging
  if (DEBUG) {
    do_output_info = 1;
    sprintf(filename,"OUTPUT_FILES/gpu_device_info_proc_%06d.txt",myrank);
  }

  // output to file
  if (do_output_info) {
    fp = fopen(filename,"w");
    if (fp != NULL) {
      // display device properties
      fprintf(fp,"Device Name = %s\n",deviceProp.name);
      fprintf(fp,"memory:\n");
      fprintf(fp,"  totalGlobalMem (in MB): %f\n",(unsigned long) deviceProp.totalGlobalMem / (1024.f * 1024.f));
      fprintf(fp,"  totalGlobalMem (in GB): %f\n",(unsigned long) deviceProp.totalGlobalMem / (1024.f * 1024.f * 1024.f));
      fprintf(fp,"  totalConstMem (in bytes): %lu\n",(unsigned long) deviceProp.totalConstMem);
      fprintf(fp,"  Maximum 1D texture size (in bytes): %lu\n",(unsigned long) deviceProp.maxTexture1D);
      fprintf(fp,"  sharedMemPerBlock (in bytes): %lu\n",(unsigned long) deviceProp.sharedMemPerBlock);
      fprintf(fp,"  regsPerBlock (in bytes): %lu\n",(unsigned long) deviceProp.regsPerBlock);
      fprintf(fp,"blocks:\n");
      fprintf(fp,"  Maximum number of threads per block: %d\n",deviceProp.maxThreadsPerBlock);
      fprintf(fp,"  Maximum size of each dimension of a block: %d x %d x %d\n",
              deviceProp.maxThreadsDim[0],deviceProp.maxThreadsDim[1],deviceProp.maxThreadsDim[2]);
      fprintf(fp,"  Maximum sizes of each dimension of a grid: %d x %d x %d\n",
              deviceProp.maxGridSize[0],deviceProp.maxGridSize[1],deviceProp.maxGridSize[2]);
      fprintf(fp,"features:\n");
      fprintf(fp,"  Compute capability of the device = %d.%d\n", deviceProp.major, deviceProp.minor);
      fprintf(fp,"  multiProcessorCount: %d\n",deviceProp.multiProcessorCount);
      if (deviceProp.canMapHostMemory) {
        fprintf(fp,"  canMapHostMemory: TRUE\n");
      }else{
        fprintf(fp,"  canMapHostMemory: FALSE\n");
      }
      if (deviceProp.deviceOverlap) {
        fprintf(fp,"  deviceOverlap: TRUE\n");
      }else{
        fprintf(fp,"  deviceOverlap: FALSE\n");
      }
      if (deviceProp.concurrentKernels) {
        fprintf(fp,"  concurrentKernels: TRUE\n");
      }else{
        fprintf(fp,"  concurrentKernels: FALSE\n");
      }
      // outputs initial memory info via cudaMemGetInfo()
      double free_db,used_db,total_db;
      get_free_memory(&free_db,&used_db,&total_db);
      fprintf(fp,"memory usage:\n");
      fprintf(fp,"  rank %d: GPU memory usage: used = %f MB, free = %f MB, total = %f MB\n",myrank,
              used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);

      // closes output file
      fclose(fp);
    }
  }

  // make sure that the device has compute capability >= 1.3
  if (deviceProp.major < 1) {
    fprintf(stderr,"Compute capability major number should be at least 1, exiting...\n\n");
    exit_on_error("CUDA Compute capability major number should be at least 1\n");
  }
  if (deviceProp.major == 1 && deviceProp.minor < 3) {
    fprintf(stderr,"Compute capability should be at least 1.3, exiting...\n");
    exit_on_error("CUDA Compute capability major number should be at least 1.3\n");
  }
  // we use pinned memory for asynchronous copy
  if (GPU_ASYNC_COPY) {
    if (! deviceProp.canMapHostMemory) {
      fprintf(stderr,"Device capability should allow to map host memory, exiting...\n");
      exit_on_error("CUDA Device capability canMapHostMemory should be TRUE\n");
    }
  }

  // checks kernel optimization setting
#ifdef USE_LAUNCH_BOUNDS
  // see: mesh_constants_cuda.h
  // performance statistics: main kernel Kernel_2_crust_mantle_impl():
  //       shared memory per block = 6200    for Kepler: total = 49152 -> limits active blocks to 7
  //       registers per thread    = 72                                   (limited by LAUNCH_MIN_BLOCKS 7)
  //       registers per block     = 9216                total = 65536    (limited by LAUNCH_MIN_BLOCKS 7)

  // shared memory
  if (deviceProp.sharedMemPerBlock > 49152 && LAUNCH_MIN_BLOCKS <= 7) {
    if (myrank == 0) {
      printf("GPU non-optimal settings: your setting of using LAUNCH_MIN_BLOCK %i is too low and limits the register usage\n",
             LAUNCH_MIN_BLOCKS);
    }
  }

  // registers
  if (deviceProp.regsPerBlock > 65536 && LAUNCH_MIN_BLOCKS <= 7) {
    if (myrank == 0) {
      printf("GPU non-optimal settings: your setting of using LAUNCH_MIN_BLOCK %i is too low and limits the register usage\n",
             LAUNCH_MIN_BLOCKS);
    }
  }
#endif

  // tests the device with a small memory allocation
  int size = 128;
  int* d_array;
  err = cudaMalloc((void**)&d_array,size*sizeof(int));
  if (err != cudaSuccess) {
    fprintf(stderr,"Error testing memory allocation on device failed\n");
    fprintf(stderr,"Error rank %d: cudaMalloc failed: %s\n", myrank,cudaGetErrorString(err));
    if (err == cudaErrorDevicesUnavailable){ fprintf(stderr,"\n%s\n", err_info); }
    exit_on_error("CUDA runtime error: cudaMalloc failed\n\n");
  }
  err = cudaFree(d_array);
  if (err != cudaSuccess) {
    fprintf(stderr,"Error cudaFree failed: %s\n", cudaGetErrorString(err));
    if (err == cudaErrorDevicesUnavailable){ fprintf(stderr,"\n%s\n", err_info); }
    exit_on_error("CUDA runtime error: cudaFree failed\n\n");
  }

  // synchronizes
  gpuSynchronize();

}
#endif

/* ----------------------------------------------------------------------------------------------- */
// OpenCL initialization
/* ----------------------------------------------------------------------------------------------- */

#ifdef USE_OPENCL

// OpenCL mesh
struct _mesh_opencl mocl;

// function definitions
cl_device_id oclGetMyDevice(int rank);
void ocl_select_device(const char *platform_filter, const char *device_filter, int myrank, int *nb_devices);
void build_kernels (void);

// initializes OpenCL devices

static void initialize_ocl_device(const char *platform_filter, const char *device_filter, int myrank, int *nb_devices) {

  // selects device
  ocl_select_device(platform_filter, device_filter, myrank, nb_devices);

  // outputs device info to file
  char filename[BUFSIZ];
  FILE *fp;
  int do_output_info = 0;

  // by default, only master process outputs device info to avoid file cluttering
  if (myrank == 0) {
    do_output_info = 1;
    sprintf(filename, "OUTPUT_FILES/gpu_device_info.txt");
  }
  // debugging
  if (DEBUG) {
    do_output_info = 1;
    sprintf(filename,"OUTPUT_FILES/gpu_device_info_proc_%06d.txt",myrank);
  }

  // output to file
  if (do_output_info) {
    fp = fopen(filename,"w");
    if (fp != NULL) {
      cl_device_type device_type;
      size_t max_work_group_size;
      cl_ulong mem_size;
      cl_uint units;
      char name[1024];
      size_t image2d_max_size[2];

      // display device properties
      clGetDeviceInfo(mocl.device, CL_DEVICE_NAME, sizeof(name), name, NULL);
      fprintf (fp, "Device Name = %s\n", name);
      clGetDeviceInfo(mocl.device, CL_DEVICE_VENDOR, sizeof(name), name, NULL);
      fprintf (fp, "Device Vendor = %s\n", name);
      fprintf (fp, "Memory:\n");
      clGetDeviceInfo(mocl.device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(mem_size), &mem_size, NULL);
      fprintf (fp, "  local_mem_size (in KB) : %f\n", mem_size / 1024.f);
      clGetDeviceInfo(mocl.device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(mem_size), &mem_size, NULL);
      fprintf (fp, "  global_mem_size (in MB): %f\n", mem_size / (1024.f * 1024.f));
      clGetDeviceInfo(mocl.device, CL_DEVICE_IMAGE2D_MAX_WIDTH, sizeof(size_t), &image2d_max_size[0], NULL);
      clGetDeviceInfo(mocl.device, CL_DEVICE_IMAGE2D_MAX_HEIGHT, sizeof(size_t), &image2d_max_size[1], NULL);
      fprintf (fp, "  image2d_max_size: %zu x %zu\n", image2d_max_size[0], image2d_max_size[1]);
      fprintf(fp,"blocks:\n");
      clGetDeviceInfo(mocl.device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(units), &units, NULL);
      fprintf (fp, "  max_compute_units: %u\n", units);
      clGetDeviceInfo(mocl.device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(max_work_group_size), &max_work_group_size, NULL);
      fprintf (fp, "  max_work_group_size: %lu\n", max_work_group_size);
      fprintf(fp,"features:\n");
      clGetDeviceInfo(mocl.device, CL_DEVICE_VERSION, sizeof(name), name, NULL);
      fprintf (fp, "  device version : %s\n", name);
      clGetDeviceInfo(mocl.device, CL_DEVICE_TYPE, sizeof(device_type), &device_type, NULL);
      fprintf (fp, "  device type: %d\n", (int) device_type);
      clGetDeviceInfo(mocl.device, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(units), &units, NULL);
      fprintf (fp, "  device max_clock_frequency: %u\n", units);
      clGetDeviceInfo(mocl.device, CL_DRIVER_VERSION, sizeof(name), name, NULL);
      fprintf (fp, "  driver version : %s\n", name);

      fclose (fp);
    }
  }

  // builds OpenCL kernels
  build_kernels();

}

#define xQUOTE(str) #str
#define QUOTE(str)  xQUOTE(str)

#ifdef OCL_GPU_CFLAGS
#define _OCL_GPU_CFLAGS QUOTE(OCL_GPU_CFLAGS)
#else
#define _OCL_GPU_CFLAGS ""
#endif

/* ----------------------------------------------------------------------------------------------- */

#define PARAMETER_STR_SIZE 1024

void build_kernels (void) {

  static char parameters[PARAMETER_STR_SIZE] = _OCL_GPU_CFLAGS " ";
  cl_int errcode;
  char *pos = parameters + strlen(_OCL_GPU_CFLAGS) + 1;
  int len = PARAMETER_STR_SIZE;
  int i;

  // adds preprocessor definitions
  // e.g. -DNDIM=3 -DNGLLX=5 ..
  for(i = 0; _macro_to_kernel[i].name != NULL; i++) {
    if (!strcmp(_macro_to_kernel[i].name, _macro_to_kernel[i].value)) {
      continue;
    }
    if (!len) {
      printf("Error: OpenCL buffer for macro parameters is not large enough, please review its size (%s:%d)\n", __FILE__, __LINE__);
      exit(1);
    }
    int written = snprintf(pos, len, "-D%s=%s ", _macro_to_kernel[i].name, _macro_to_kernel[i].value);
    pos += written;
    len -= written;
  }

  // debug
  //printf("building OpenCL kernels: parameters = %s \n",parameters);

  // adds kernels as const char definitions
  #include "kernel_inc_cl.c"

  // defines OpenCL build program macro
#undef BOAST_KERNEL
#define BOAST_KERNEL(__kern_name__)                                     \
  mocl.programs.__kern_name__##_program = clCreateProgramWithSource( mocl.context, 1, \
                                                                     &__kern_name__##_program, NULL, clck_(&errcode));\
  mocl_errcode = clBuildProgram(mocl.programs.__kern_name__##_program, 0, NULL, parameters, NULL, NULL);\
  if (mocl_errcode != CL_SUCCESS) {                                     \
    fprintf(stderr,"OpenCL Error: Failed to build program "#__kern_name__": %s\n", clewErrorString(mocl_errcode)); \
    char cBuildLog[10240];                                              \
    clGetProgramBuildInfo(mocl.programs.__kern_name__##_program, mocl.device, CL_PROGRAM_BUILD_LOG, \
                          sizeof(cBuildLog), cBuildLog, NULL ); \
    fprintf(stderr,"OpenCL Log: %s\n",cBuildLog);                                   \
    exit(1);                                                            \
  }                                                                     \
  mocl.kernels.__kern_name__ = clCreateKernel (mocl.programs.__kern_name__ ## _program, #__kern_name__ , clck_(&errcode));

  // builds each OpenCL kernel
  #include "kernel_list.h"

}


/* ----------------------------------------------------------------------------------------------- */

struct _opencl_version {
  cl_uint minor;
  cl_uint major;
};
struct _opencl_version opencl_version_1_0 = {1,0};
struct _opencl_version opencl_version_1_1 = {1,1};
struct _opencl_version opencl_version_1_2 = {1,2};

/* ----------------------------------------------------------------------------------------------- */

cl_int compare_opencl_version(struct _opencl_version v1, struct _opencl_version v2) {
  if (v1.major > v2.major)
    return 1;
  if (v1.major < v2.major)
    return -1;
  if (v1.minor > v2.minor)
    return 1;
  if (v1.minor < v2.minor)
    return -1;
  return 0;
}

/* ----------------------------------------------------------------------------------------------- */

static void get_platform_version(cl_platform_id platform_id, struct _opencl_version *version) {

  size_t cl_platform_version_size;
  clCheck(clGetPlatformInfo(platform_id, CL_PLATFORM_VERSION, 0, NULL, &cl_platform_version_size));

  char *cl_platform_version;
  cl_platform_version = (char *) malloc(cl_platform_version_size);

  if (cl_platform_version == NULL) {
    fprintf(stderr,"Error: Failed to create string (out of memory)!\n");
    exit(1);
  }

  clCheck(clGetPlatformInfo(platform_id, CL_PLATFORM_VERSION, cl_platform_version_size, cl_platform_version, NULL));

  //OpenCL<space><major_version.minor_version><space><platform-specific information>
  char minor[2], major[2];
  major[0] = cl_platform_version[7];
  major[1] = 0;
  minor[0] = cl_platform_version[9];
  minor[1] = 0;
  version->major = atoi(major);
  version->major = atoi(minor);
  free(cl_platform_version);
}

/* ----------------------------------------------------------------------------------------------- */

#define OCL_DEV_TYPE CL_DEVICE_TYPE_ALL

void ocl_select_device(const char *platform_filter, const char *device_filter, int myrank, int *nb_devices) {

  cl_int errcode = CL_SUCCESS;
  cl_platform_id *platform_ids;
  cl_uint num_platforms;

  // first OpenCL call
  // only gets number of platforms
  clCheck( clGetPlatformIDs(0, NULL, &num_platforms) );

  // checks if OpenCL platforms available
  if (num_platforms == 0) {
    fprintf(stderr,"OpenCL error: No OpenCL platform available!\n");
    exit(1);
  }

  platform_ids = (cl_platform_id *) malloc(num_platforms * sizeof(cl_platform_id));

  // gets platform infos
  clCheck( clGetPlatformIDs(num_platforms, platform_ids, NULL));

  cl_context_properties properties[] = {CL_CONTEXT_PLATFORM, 0, 0 };

  // temporary array to store infos
  int i,j;
  char *info_all[num_platforms][2];
  // initializes pointers
  for (i = 0; i < num_platforms; i++) {
    info_all[i][0] = NULL;
    info_all[i][1] = NULL;
  }

  // looks for platform matching GPU_PLATFORM string given in Par_file
  if (strlen(platform_filter)) {
    cl_uint found = 0;

    for (i = 0; i < num_platforms && !found; i++) {
      size_t info_length;
      char *info;

      // checks vendor and platform names for matching with GPU_PLATFORM
      int props_to_check[] = {CL_PLATFORM_VENDOR, CL_PLATFORM_NAME};
      for (j = 0; j < 2 && !found; j++) {
        // gets property info length
        clCheck( clGetPlatformInfo(platform_ids[i], props_to_check[j], 0, NULL, &info_length));

        // checks info
        if (info_length == 0) {
          fprintf(stderr,"OpenCL error: No OpenCL platform info available!\n");
          exit(1);
        }

        // allocates info buffer and gets info string
        info = (char *) malloc(info_length * sizeof(char));
        clCheck( clGetPlatformInfo(platform_ids[i], props_to_check[j], info_length, info, NULL));

        // stores info
        info_all[i][j] = (char *) malloc( strlen(info) + 1);
        strcpy(info_all[i][j],info);

        // sets matching platform id
        if (strcasestr(info, platform_filter)) {
          properties[1] = (cl_context_properties) platform_ids[i];
          found = 1;
        }
        // frees temporary array
        free(info);
      }
    }

    // checks if platform found
    if (!found) {
      if (myrank == 0) {
        fprintf(stderr, "\nAvailable platforms are:\n");
        for (i = 0; i < num_platforms; i++) {
          if (info_all[i][0]) { fprintf(stderr, "  platform %i: vendor = %s , name = %s\n",i,info_all[i][0],info_all[i][1]);}
        }
        fprintf(stderr, "Please check your parameter GPU_PLATFORM in Par_file\n\n");
      }
      // frees info array
      for (i = 0; i < num_platforms; i++) {
        if (info_all[i][0]) { free(info_all[i][0]); }
        if (info_all[i][1]) { free(info_all[i][1]); }
      }
      // exits
      fprintf(stderr, "No matching OpenCL platform available : %s\n", platform_filter);
      exit(1);
    }

    // frees info array
    for (i = 0; i < num_platforms; i++) {
      if (info_all[i][0]) { free(info_all[i][0]); }
      if (info_all[i][1]) { free(info_all[i][1]); }
    }

  } else {
    // wild-card platform filter given (GPU_PLATFORM set to '*'), takes first platform
    properties[1] = (cl_context_properties) platform_ids[0];
  }

  // searches for device
  if (strlen(device_filter)) {
    cl_uint found = 0;
    cl_uint i;
    cl_uint num_devices;
    cl_device_id *device_ids;
    cl_device_id *matching_device_ids;

    // only gets number of devices for this platform
    clCheck( clGetDeviceIDs((cl_platform_id) properties[1], OCL_DEV_TYPE, 0, NULL, &num_devices));

    // checks
    if (num_devices == 0) {
      fprintf(stderr,"No OpenCL device of type %d!\n", (int) OCL_DEV_TYPE);
      exit(1);
    }

    device_ids = (cl_device_id *) malloc(num_devices * sizeof(cl_device_id));

    matching_device_ids = (cl_device_id *) malloc(num_devices * sizeof(cl_device_id));

    // gets device infos
    clCheck( clGetDeviceIDs((cl_platform_id) properties[1], OCL_DEV_TYPE, num_devices, device_ids, NULL));

    // temporary array to store device infos
    char *info_device_all[num_devices];
    // initializes pointers
    for (i = 0; i < num_devices; i++) {
      info_device_all[i] = NULL;
    }

    // searches device matching GPU_DEVICE string
    for (i = 0; i < num_devices; i++) {
      size_t info_length;
      char *info;

      clCheck( clGetDeviceInfo(device_ids[i], CL_DEVICE_NAME, 0, NULL, &info_length));

      info = (char *) malloc(info_length * sizeof(char));

      clCheck( clGetDeviceInfo(device_ids[i], CL_DEVICE_NAME, info_length, info, NULL));

      // stores info
      info_device_all[i] = (char *) malloc( strlen(info) + 1);
      strcpy(info_device_all[i],info);

      // sets matching device id
      if (strcasestr(info, device_filter)) {
        matching_device_ids[found] = device_ids[i];
        found++;
      }

      free(info);
    }

    if (!found) {
      // user output
      if (myrank == 0) {
        fprintf(stderr, "\nAvailable devices are:\n");
        for (i = 0; i < num_devices; i++) {
          if (info_device_all[i]) { fprintf(stderr, "  device %i: name = %s\n",i,info_device_all[i]);}
        }
        fprintf(stderr, "Please check your parameter GPU_DEVICE in Par_file\n\n");
      }
      // frees info array
      for (i = 0; i < num_devices; i++) {
        if (info_device_all[i]) { free(info_device_all[i]); }
      }
      // exits
      fprintf(stderr, "No matching OpenCL device available : %s\n", device_filter);
      exit(1);
    }

    // creates an OpenCL context
    mocl.context = clCreateContext(properties, found, matching_device_ids, NULL, NULL, clck_(&errcode));

    // frees temporary arrays
    free (matching_device_ids);
    free (device_ids);
    // frees info array
    for (i = 0; i < num_devices; i++) {
      if (info_device_all[i]) { free(info_device_all[i]); }
    }

  } else {
    // wild-card GPU_DEVICE set to '*'
    mocl.context = clCreateContextFromType(properties, OCL_DEV_TYPE, NULL, NULL, clck_(&errcode));
  }

  //get the number of devices available in the context (devices which are of DEVICE_TYPE_GPU of platform platform_ids[0])
  struct _opencl_version  platform_version;
  get_platform_version((cl_platform_id) properties[1], &platform_version);

#ifdef CL_VERSION_1_1
  if (compare_opencl_version(platform_version, opencl_version_1_1) >= 0 ) {
    clGetContextInfo(mocl.context, CL_CONTEXT_NUM_DEVICES, sizeof(*nb_devices), nb_devices, NULL);
  } else
#endif
  {
    size_t nContextDescriptorSize;
    clGetContextInfo(mocl.context, CL_CONTEXT_DEVICES, 0, 0, &nContextDescriptorSize);
    *nb_devices = nContextDescriptorSize / sizeof(cl_device_id);
  }

  // stores info in mesh opencl structure
  mocl.nb_devices = *nb_devices;
  free(platform_ids);

  size_t szParmDataBytes;
  cl_device_id* cdDevices;

  // get the list of GPU devices associated with this context
  clGetContextInfo(mocl.context, CL_CONTEXT_DEVICES, 0, NULL, &szParmDataBytes);
  cdDevices = (cl_device_id *) malloc(szParmDataBytes);

  clGetContextInfo(mocl.context, CL_CONTEXT_DEVICES, szParmDataBytes, cdDevices, NULL);

  mocl.device = cdDevices[myrank % mocl.nb_devices];
  free(cdDevices);

  // command kernel queues
  mocl.command_queue = clCreateCommandQueue(mocl.context, mocl.device, 0, clck_(&errcode));
  if (GPU_ASYNC_COPY) {
    mocl.copy_queue = clCreateCommandQueue(mocl.context, mocl.device, 0, clck_(&errcode));
  }
}
#endif

/* ----------------------------------------------------------------------------------------------- */
// GPU initialization
/* ----------------------------------------------------------------------------------------------- */

#define isspace(c) ((c) == ' ')

static char *trim_and_default(char *s)
{
  // trim before
  while (*s != '\0' && isspace(*s)) { s++; }

  if (*s == '\0') {
    return s;
  }

  // note: the platform_filter argument acts weird on apple platforms, giving a string "NVIDIA   Geforce", instead of just "NVIDIA" and "Geforce"
  //       here we assume that maximum length of GPU_PLATFORM is 11 characters
  //       todo - find better way to avoid this?
  // debug
  //printf("string: %s has length %i \n",s,strlen(s));
  int len = strlen(s);
  if (len > 11 ) len = 11;

  // trim after
  char *back = s + len;
  while (isspace(*--back));
  *(back + 1) = '\0';

  // replace * by empty string
  if (strlen(s) == 1 && *s == '*') {
    *s = '\0';
  }

  return s;
}

/* ----------------------------------------------------------------------------------------------- */

enum gpu_runtime_e {COMPILE, CUDA, OPENCL};

extern EXTERN_LANG
void FC_FUNC_ (initialize_gpu_device,
               INITIALIZE_GPU_DEVICE) (int *runtime_f, char *platform_filter, char *device_filter, int *myrank_f, int *nb_devices) {

  TRACE ("initialize_device");

  enum gpu_runtime_e runtime_type = (enum gpu_runtime_e) *runtime_f;

  // trims GPU_PLATFORM and GPU_DEVICE strings
  platform_filter = trim_and_default(platform_filter);
  device_filter = trim_and_default(device_filter);

  // sets and checks gpu runtime flags
#if defined(USE_OPENCL) && defined(USE_CUDA)
  run_cuda = runtime_type == CUDA;
  run_opencl = runtime_type == OPENCL;
  if (runtime_type == COMPILE) {
    if (*myrank_f == 0) {
      printf("\
Error: GPU_RUNTIME set to compile time decision (%d), but both OpenCL (%d) and CUDA (%d) are compiled.\n\
Please set Par_file accordingly...\n\n", COMPILE, OPENCL, CUDA);
    }
    exit(1);
  }
#elif defined(USE_OPENCL)
  run_opencl = 1;
  if (runtime_type != COMPILE && runtime_type != OPENCL) {
    if (*myrank_f == 0) {
      printf("\
Warning: GPU_RUNTIME parameter in Par_file set to (%d) is incompatible with OpenCL-only compilation (OPENCL=%d, COMPILE=%d).\n\
This simulation will continue using the OpenCL runtime...\n\n", runtime_type, OPENCL, COMPILE);
    }
  }
#elif defined(USE_CUDA)
  run_cuda = 1;
  if (runtime_type != COMPILE && runtime_type != CUDA) {
    if (*myrank_f == 0) {
      printf("\
Warning: GPU_RUNTIME parameter in Par_file set to (%d) is incompatible with Cuda-only compilation (CUDA=%d, COMPILE=%d).\n\
This simulation will continue using the Cuda runtime...\n", runtime_type, CUDA, COMPILE);
    }
  }
#else
  #error "GPU code compiled but neither CUDA nor OpenCL are enabled"
#endif

  // initializes gpu cards
#ifdef USE_OPENCL
  if (run_opencl) {
    initialize_ocl_device(platform_filter, device_filter, *myrank_f, nb_devices);
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    initialize_cuda_device(platform_filter, device_filter, *myrank_f, nb_devices);
  }
#endif
}
