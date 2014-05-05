const char * compute_coupling_fluid_CMB_kernel_program = "\
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
__kernel void compute_coupling_fluid_CMB_kernel(const __global float * displ_crust_mantle, __global float * accel_outer_core, const __global int * ibool_crust_mantle, const __global int * ibelm_bottom_crust_mantle, const __global float * normal_top_outer_core, const __global float * jacobian2D_top_outer_core, const __global float * wgllwgll_xy, const __global int * ibool_outer_core, const __global int * ibelm_top_outer_core, const int NSPEC2D_TOP_OC){\n\
  int i;\n\
  int j;\n\
  int k;\n\
  int iface;\n\
  int k_corresp;\n\
  float displ_n;\n\
  int iglob_cm;\n\
  int iglob_oc;\n\
  int ispec;\n\
  int ispec_selected;\n\
  float displ_x;\n\
  float displ_y;\n\
  float displ_z;\n\
  float nx;\n\
  float ny;\n\
  float nz;\n\
  float weight;\n\
  i = get_local_id(0);\n\
  j = get_local_id(1);\n\
  iface = get_group_id(0) + (get_num_groups(0)) * (get_group_id(1));\n\
  if(iface < NSPEC2D_TOP_OC){\n\
    ispec = ibelm_top_outer_core[iface - 0] - (1);\n\
    ispec_selected = ibelm_bottom_crust_mantle[iface - 0] - (1);\n\
    k = NGLLX - (1);\n\
    k_corresp = 0;\n\
    iglob_cm = ibool_crust_mantle[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k_corresp, ispec_selected) - 0] - (1);\n\
    displ_x = displ_crust_mantle[(iglob_cm) * (3) + 0 - 0];\n\
    displ_y = displ_crust_mantle[(iglob_cm) * (3) + 1 - 0];\n\
    displ_z = displ_crust_mantle[(iglob_cm) * (3) + 2 - 0];\n\
    nx = normal_top_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 0, i, j, iface) - 0];\n\
    ny = normal_top_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 1, i, j, iface) - 0];\n\
    nz = normal_top_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 2, i, j, iface) - 0];\n\
    displ_n = (displ_x) * (nx) + (displ_y) * (ny) + (displ_z) * (nz);\n\
    weight = (jacobian2D_top_outer_core[INDEX3(NGLLX, NGLLX, i, j, iface) - 0]) * (wgllwgll_xy[INDEX2(NGLLX, i, j) - 0]);\n\
    iglob_oc = ibool_outer_core[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] - (1);\n\
    atomicAdd(accel_outer_core + iglob_oc, (weight) * (displ_n));\n\
  }\n\
}\n\
";
