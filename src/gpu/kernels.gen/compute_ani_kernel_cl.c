const char * compute_ani_kernel_program = "\
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
void compute_strain_product(float * prod, const float eps_trace_over_3, const float * epsdev, const float b_eps_trace_over_3, const float * b_epsdev){\n\
  float eps[6];\n\
  float b_eps[6];\n\
  int p;\n\
  int i;\n\
  int j;\n\
  eps[0 - 0] = epsdev[0 - 0] + eps_trace_over_3;\n\
  eps[1 - 0] = epsdev[1 - 0] + eps_trace_over_3;\n\
  eps[2 - 0] =  -(eps[0 - 0] + eps[1 - 0]) + (eps_trace_over_3) * (3.0f);\n\
  eps[3 - 0] = epsdev[4 - 0];\n\
  eps[4 - 0] = epsdev[3 - 0];\n\
  eps[5 - 0] = epsdev[2 - 0];\n\
  b_eps[0 - 0] = b_epsdev[0 - 0] + b_eps_trace_over_3;\n\
  b_eps[1 - 0] = b_epsdev[1 - 0] + b_eps_trace_over_3;\n\
  b_eps[2 - 0] =  -(b_eps[0 - 0] + b_eps[1 - 0]) + (b_eps_trace_over_3) * (3.0f);\n\
  b_eps[3 - 0] = b_epsdev[4 - 0];\n\
  b_eps[4 - 0] = b_epsdev[3 - 0];\n\
  b_eps[5 - 0] = b_epsdev[2 - 0];\n\
  p = 0;\n\
  for(i=0; i<=5; i+=1){\n\
    for(j=0; j<=5; j+=1){\n\
      prod[p - 0] = (eps[i - 0]) * (b_eps[j - 0]);\n\
      if(j > i){\n\
        prod[p - 0] = prod[p - 0] + (eps[j - 0]) * (b_eps[i - 0]);\n\
        if(j > 2 && i < 3){\n\
          prod[p - 0] = (prod[p - 0]) * (2.0f);\n\
        }\n\
      }\n\
    }\n\
  }\n\
}\n\
__kernel void compute_ani_kernel(const __global float * epsilondev_xx, const __global float * epsilondev_yy, const __global float * epsilondev_xy, const __global float * epsilondev_xz, const __global float * epsilondev_yz, const __global float * epsilon_trace_over_3, const __global float * b_epsilondev_xx, const __global float * b_epsilondev_yy, const __global float * b_epsilondev_xy, const __global float * b_epsilondev_xz, const __global float * b_epsilondev_yz, const __global float * b_epsilon_trace_over_3, __global float * cijkl_kl, const int NSPEC, const float deltat){\n\
  int i;\n\
  int ispec;\n\
  int ijk_ispec;\n\
  float eps_trace_over_3;\n\
  float b_eps_trace_over_3;\n\
  float prod[21];\n\
  float epsdev[5];\n\
  float b_epsdev[5];\n\
  ispec = get_group_id(0) + (get_group_id(1)) * (get_num_groups(0));\n\
  if(ispec < NSPEC){\n\
    ijk_ispec = get_local_id(0) + (NGLL3) * (ispec);\n\
    epsdev[0 - 0] = epsilondev_xx[ijk_ispec - 0];\n\
    epsdev[1 - 0] = epsilondev_yy[ijk_ispec - 0];\n\
    epsdev[2 - 0] = epsilondev_xy[ijk_ispec - 0];\n\
    epsdev[3 - 0] = epsilondev_xz[ijk_ispec - 0];\n\
    epsdev[4 - 0] = epsilondev_yz[ijk_ispec - 0];\n\
    epsdev[0 - 0] = b_epsilondev_xx[ijk_ispec - 0];\n\
    epsdev[1 - 0] = b_epsilondev_yy[ijk_ispec - 0];\n\
    epsdev[2 - 0] = b_epsilondev_xy[ijk_ispec - 0];\n\
    epsdev[3 - 0] = b_epsilondev_xz[ijk_ispec - 0];\n\
    epsdev[4 - 0] = b_epsilondev_yz[ijk_ispec - 0];\n\
    eps_trace_over_3 = epsilon_trace_over_3[ijk_ispec - 0];\n\
    b_eps_trace_over_3 = b_epsilon_trace_over_3[ijk_ispec - 0];\n\
    compute_strain_product(prod, eps_trace_over_3, epsdev, b_eps_trace_over_3, b_epsdev);\n\
    for(i=0; i<=20; i+=1){\n\
      cijkl_kl[i - 0 + (ijk_ispec - (0)) * (21)] = cijkl_kl[i - 0 + (ijk_ispec - (0)) * (21)] + (deltat) * (prod[i - 0]);\n\
    }\n\
  }\n\
}\n\
";
