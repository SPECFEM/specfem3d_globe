const char * noise_add_surface_movie_kernel_program = "\
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
__kernel void noise_add_surface_movie_kernel(__global float * accel, const __global int * ibool, const __global int * ibelm_top, const int nspec_top, const __global float * noise_surface_movie, const __global float * normal_x_noise, const __global float * normal_y_noise, const __global float * normal_z_noise, const __global float * mask_noise, const __global float * jacobian2D, const __global float * wgllwgll){\n\
  int igll;\n\
  int iface;\n\
  igll = get_local_id(0);\n\
  iface = get_group_id(0) + (get_group_id(1)) * (get_num_groups(0));\n\
  if(iface < nspec_top){\n\
    int i;\n\
    int j;\n\
    int k;\n\
    int ispec;\n\
    int iglob;\n\
    int ipoin;\n\
    float eta;\n\
    float jacobianw;\n\
    float normal_x;\n\
    float normal_y;\n\
    float normal_z;\n\
    ispec = ibelm_top[iface - 0] - (1);\n\
    k = NGLLX - (1);\n\
    j = (igll) / (NGLLX);\n\
    i = igll - ((j) * (NGLLX));\n\
    iglob = ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] - (1);\n\
    ipoin = (NGLL2) * (iface) + igll;\n\
    normal_x = normal_x_noise[ipoin - 0];\n\
    normal_y = normal_y_noise[ipoin - 0];\n\
    normal_z = normal_z_noise[ipoin - 0];\n\
    eta = 0.0f;\n\
    eta = eta + (noise_surface_movie[INDEX3(NDIM, NGLL2, 0, igll, iface) - 0]) * (normal_x);\n\
    eta = eta + (noise_surface_movie[INDEX3(NDIM, NGLL2, 1, igll, iface) - 0]) * (normal_y);\n\
    eta = eta + (noise_surface_movie[INDEX3(NDIM, NGLL2, 2, igll, iface) - 0]) * (normal_z);\n\
    jacobianw = (wgllwgll[(k) * (NGLLX) + i - 0]) * (jacobian2D[igll + (NGLL2) * (iface) - 0]);\n\
    atomicAdd(accel + (iglob) * (3) + 0, (((eta) * (mask_noise[ipoin - 0])) * (normal_x)) * (jacobianw));\n\
    atomicAdd(accel + (iglob) * (3) + 1, (((eta) * (mask_noise[ipoin - 0])) * (normal_y)) * (jacobianw));\n\
    atomicAdd(accel + (iglob) * (3) + 2, (((eta) * (mask_noise[ipoin - 0])) * (normal_z)) * (jacobianw));\n\
  }\n\
}\n\
";
