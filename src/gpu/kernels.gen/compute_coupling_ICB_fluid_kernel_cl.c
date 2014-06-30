const char * compute_coupling_ICB_fluid_kernel_program = "\
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
__kernel void compute_coupling_ICB_fluid_kernel(const __global float * displ_inner_core, __global float * accel_inner_core, const __global float * accel_outer_core, const __global int * ibool_inner_core, const __global int * ibelm_top_inner_core, const __global float * normal_bottom_outer_core, const __global float * jacobian2D_bottom_outer_core, const __global float * wgllwgll_xy, const __global int * ibool_outer_core, const __global int * ibelm_bottom_outer_core, const float RHO_BOTTOM_OC, const float minus_g_icb, int GRAVITY, const int NSPEC2D_TOP_IC){\n\
  int i;\n\
  int j;\n\
  int k;\n\
  int iface;\n\
  int k_corresp;\n\
  int iglob_oc;\n\
  int iglob_ic;\n\
  float pressure;\n\
  int ispec;\n\
  int ispec_selected;\n\
  float nx;\n\
  float ny;\n\
  float nz;\n\
  float weight;\n\
  i = get_local_id(0);\n\
  j = get_local_id(1);\n\
  iface = get_group_id(0) + (get_num_groups(0)) * (get_group_id(1));\n\
  if(iface < NSPEC2D_TOP_IC){\n\
    ispec = ibelm_top_inner_core[iface - (0)] - (1);\n\
    ispec_selected = ibelm_bottom_outer_core[iface - (0)] - (1);\n\
    k = NGLLX - (1);\n\
    k_corresp = 0;\n\
    iglob_oc = ibool_outer_core[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k_corresp, ispec_selected) - (0)] - (1);\n\
    nx = normal_bottom_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 0, i, j, iface) - (0)];\n\
    ny = normal_bottom_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 1, i, j, iface) - (0)];\n\
    nz = normal_bottom_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 2, i, j, iface) - (0)];\n\
    weight = (jacobian2D_bottom_outer_core[INDEX3(NGLLX, NGLLX, i, j, iface) - (0)]) * (wgllwgll_xy[INDEX2(NGLLX, i, j) - (0)]);\n\
    iglob_ic = ibool_inner_core[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - (0)] - (1);\n\
    if(GRAVITY){\n\
      pressure = (RHO_BOTTOM_OC) * ((minus_g_icb) * ((displ_inner_core[(iglob_ic) * (3) - (0)]) * (nx) + (displ_inner_core[(iglob_ic) * (3) + 1 - (0)]) * (ny) + (displ_inner_core[(iglob_ic) * (3) + 2 - (0)]) * (nz)) - (accel_outer_core[iglob_oc - (0)]));\n\
    } else {\n\
      pressure = ( -(RHO_BOTTOM_OC)) * (accel_outer_core[iglob_oc - (0)]);\n\
    }\n\
    atomicAdd(accel_inner_core + (iglob_ic) * (3) + 0, (( -(weight)) * (nx)) * (pressure));\n\
    atomicAdd(accel_inner_core + (iglob_ic) * (3) + 1, (( -(weight)) * (ny)) * (pressure));\n\
    atomicAdd(accel_inner_core + (iglob_ic) * (3) + 2, (( -(weight)) * (nz)) * (pressure));\n\
  }\n\
}\n\
";
