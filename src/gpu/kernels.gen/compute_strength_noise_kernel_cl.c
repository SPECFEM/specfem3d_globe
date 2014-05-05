const char * compute_strength_noise_kernel_program = "\
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
__kernel void compute_strength_noise_kernel(const __global float * displ, const __global int * ibelm_top, const __global int * ibool, const __global float * noise_surface_movie, const __global float * normal_x_noise, const __global float * normal_y_noise, const __global float * normal_z_noise, __global float * Sigma_kl, const float deltat, const int nspec_top){\n\
  int iface;\n\
  int ispec;\n\
  int igll;\n\
  int ipoin;\n\
  int i;\n\
  int j;\n\
  int k;\n\
  int iglob;\n\
  float eta;\n\
  iface = get_group_id(0) + (get_group_id(1)) * (get_num_groups(0));\n\
  if(iface < nspec_top){\n\
    ispec = ibelm_top[iface - 0] - (1);\n\
    igll = get_local_id(0);\n\
    ipoin = igll + (NGLL2) * (iface);\n\
    k = NGLLX - (1);\n\
    j = (igll) / (NGLLX);\n\
    i = igll - ((j) * (NGLLX));\n\
    iglob = ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] - (1);\n\
    eta = (noise_surface_movie[INDEX3(NDIM, NGLL2, 0, igll, iface) - 0]) * (normal_x_noise[ipoin - 0]) + (noise_surface_movie[INDEX3(NDIM, NGLL2, 1, igll, iface) - 0]) * (normal_y_noise[ipoin - 0]) + (noise_surface_movie[INDEX3(NDIM, NGLL2, 2, igll, iface) - 0]) * (normal_z_noise[ipoin - 0]);\n\
    Sigma_kl[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] = Sigma_kl[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] + ((deltat) * (eta)) * ((normal_x_noise[ipoin - 0]) * (displ[0 - 0 + (iglob - (0)) * (3)]) + (normal_y_noise[ipoin - 0]) * (displ[1 - 0 + (iglob - (0)) * (3)]) + (normal_z_noise[ipoin - 0]) * (displ[2 - 0 + (iglob - (0)) * (3)]));\n\
  }\n\
}\n\
";
