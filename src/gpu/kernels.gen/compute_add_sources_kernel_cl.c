const char * compute_add_sources_kernel_program = "\
#pragma OPENCL EXTENSION cl_khr_fp64: enable\n\
#pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable\n\
inline void atomicAdd(volatile __global float *source, const float val) {\n\
  union {\n\
    unsigned int iVal;\n\
    float fVal;\n\
  } res, orig;\n\
  do {\n\
    orig.fVal = *source;\n\
    res.fVal = orig.fVal + val;\n\
  } while (atomic_cmpxchg((volatile __global unsigned int *)source, orig.iVal, res.iVal) != orig.iVal);\n\
}\n\
#ifndef INDEX2\n\
#define INDEX2(isize,i,j) i + isize*j\n\
#endif\n\
#ifndef INDEX3\n\
#define INDEX3(isize,jsize,i,j,k) i + isize*(j + jsize*k)\n\
#endif\n\
#ifndef INDEX4\n\
#define INDEX4(isize,jsize,ksize,i,j,k,x) i + isize*(j + jsize*(k + ksize*x))\n\
#endif\n\
#ifndef INDEX5\n\
#define INDEX5(isize,jsize,ksize,xsize,i,j,k,x,y) i + isize*(j + jsize*(k + ksize*(x + xsize*y)))\n\
#endif\n\
#ifndef NDIM\n\
#define NDIM 3\n\
#endif\n\
#ifndef NGLLX\n\
#define NGLLX 5\n\
#endif\n\
#ifndef NGLL2\n\
#define NGLL2 25\n\
#endif\n\
#ifndef NGLL3\n\
#define NGLL3 125\n\
#endif\n\
#ifndef NGLL3_PADDED\n\
#define NGLL3_PADDED 128\n\
#endif\n\
#ifndef N_SLS\n\
#define N_SLS 3\n\
#endif\n\
#ifndef IREGION_CRUST_MANTLE\n\
#define IREGION_CRUST_MANTLE 1\n\
#endif\n\
#ifndef IREGION_INNER_CORE\n\
#define IREGION_INNER_CORE 3\n\
#endif\n\
#ifndef IFLAG_IN_FICTITIOUS_CUBE\n\
#define IFLAG_IN_FICTITIOUS_CUBE 11\n\
#endif\n\
#ifndef R_EARTH_KM\n\
#define R_EARTH_KM 6371.0f\n\
#endif\n\
#ifndef COLORING_MIN_NSPEC_INNER_CORE\n\
#define COLORING_MIN_NSPEC_INNER_CORE 1000\n\
#endif\n\
#ifndef COLORING_MIN_NSPEC_OUTER_CORE\n\
#define COLORING_MIN_NSPEC_OUTER_CORE 1000\n\
#endif\n\
#ifndef BLOCKSIZE_TRANSFER\n\
#define BLOCKSIZE_TRANSFER 256\n\
#endif\n\
__kernel void compute_add_sources_kernel(__global float * accel, const __global int * ibool, const __global float * sourcearrays, const __global double * stf_pre_compute, const int myrank, const __global int * islice_selected_source, const __global int * ispec_selected_source, const int NSOURCES){\n\
  int ispec;\n\
  int iglob;\n\
  float stf;\n\
  int isource;\n\
  int i;\n\
  int j;\n\
  int k;\n\
  i = get_local_id(0);\n\
  j = get_local_id(1);\n\
  k = get_local_id(2);\n\
  isource = get_group_id(0) + (get_num_groups(0)) * (get_group_id(1));\n\
  if(isource < NSOURCES){\n\
    if(myrank == islice_selected_source[isource - 0]){\n\
      ispec = ispec_selected_source[isource - 0] - (1);\n\
      stf = stf_pre_compute[isource - 0];\n\
      iglob = ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] - (1);\n\
      atomicAdd(accel + (iglob) * (3) + 0, (sourcearrays[INDEX5(NDIM, NGLLX, NGLLX, NGLLX, 0, i, j, k, isource) - 0]) * (stf));\n\
      atomicAdd(accel + (iglob) * (3) + 1, (sourcearrays[INDEX5(NDIM, NGLLX, NGLLX, NGLLX, 1, i, j, k, isource) - 0]) * (stf));\n\
      atomicAdd(accel + (iglob) * (3) + 2, (sourcearrays[INDEX5(NDIM, NGLLX, NGLLX, NGLLX, 2, i, j, k, isource) - 0]) * (stf));\n\
    }\n\
  }\n\
}\n\
";
