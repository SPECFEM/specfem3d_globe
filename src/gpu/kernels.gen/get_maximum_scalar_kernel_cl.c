const char * get_maximum_scalar_kernel_program = "\
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
#define INDEX2(xsize,x,y) x + (y)*xsize\n\
#endif\n\
#ifndef INDEX3\n\
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)\n\
#endif\n\
#ifndef INDEX4\n\
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))\n\
#endif\n\
#ifndef INDEX5\n\
#define INDEX5(xsize,ysize,zsize,isize,x,y,z,i,j) x + xsize*(y + ysize*(z + zsize*(i + isize*(j))))\n\
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
__kernel void get_maximum_scalar_kernel(const __global float * array, const int size, __global float * d_max){\n\
  __local float sdata[BLOCKSIZE_TRANSFER + 0 - (1) - (0) + 1];\n\
  int tid;\n\
  int bx;\n\
  int i;\n\
  int s;\n\
  tid = get_local_id(0);\n\
  bx = (get_group_id(1)) * (get_num_groups(0)) + get_group_id(0);\n\
  i = tid + (bx) * (get_local_size(0));\n\
  sdata[tid - 0] = (i < size ? fabs(array[i - 0]) : 0.0f);\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  s = (get_local_size(0)) / (2);\n\
  while(s > 0){\n\
    if(tid < s){\n\
      if(sdata[tid - 0] < sdata[tid + s - 0]){\n\
        sdata[tid - 0] = sdata[tid + s - 0];\n\
      }\n\
    }\n\
    s = s >> 1;\n\
    barrier(CLK_LOCAL_MEM_FENCE);\n\
  }\n\
  if(tid == 0){\n\
    d_max[bx - 0] = sdata[0 - 0];\n\
  }\n\
}\n\
";
