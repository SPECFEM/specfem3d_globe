const char * compute_stacey_elastic_backward_kernel_program = "\
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
__kernel void compute_stacey_elastic_backward_kernel(__global float * b_accel, const __global float * b_absorb_field, const int interface_type, const int num_abs_boundary_faces, const __global int * abs_boundary_ispec, const __global int * nkmin_xi, const __global int * nkmin_eta, const __global int * njmin, const __global int * njmax, const __global int * nimin, const __global int * nimax, const __global int * ibool){\n\
  int igll;\n\
  int iface;\n\
  int i;\n\
  int j;\n\
  int k;\n\
  int iglob;\n\
  int ispec;\n\
  igll = get_local_id(0);\n\
  iface = get_group_id(0) + (get_group_id(1)) * (get_num_groups(0));\n\
  if(iface < num_abs_boundary_faces){\n\
    ispec = abs_boundary_ispec[iface - 0] - (1);\n\
    switch(interface_type){\n\
      case 0 :\n\
        if(nkmin_xi[INDEX2(2, 0, iface) - 0] == 0 || njmin[INDEX2(2, 0, iface) - 0] == 0){\n\
           return ;\n\
        }\n\
        i = 0;\n\
        k = (igll) / (NGLLX);\n\
        j = igll - ((k) * (NGLLX));\n\
        if(k < nkmin_xi[INDEX2(2, 0, iface) - 0] - (1) || k > NGLLX - (1)){\n\
           return ;\n\
        }\n\
        if(j < njmin[INDEX2(2, 0, iface) - 0] - (1) || j > NGLLX - (1)){\n\
           return ;\n\
        }\n\
        break;\n\
      case 1 :\n\
        if(nkmin_xi[INDEX2(2, 1, iface) - 0] == 0 || njmin[INDEX2(2, 1, iface) - 0] == 0){\n\
           return ;\n\
        }\n\
        i = NGLLX - (1);\n\
        k = (igll) / (NGLLX);\n\
        j = igll - ((k) * (NGLLX));\n\
        if(k < nkmin_xi[INDEX2(2, 1, iface) - 0] - (1) || k > NGLLX - (1)){\n\
           return ;\n\
        }\n\
        if(j < njmin[INDEX2(2, 1, iface) - 0] - (1) || j > njmax[INDEX2(2, 1, iface) - 0] - (1)){\n\
           return ;\n\
        }\n\
        break;\n\
      case 2 :\n\
        if(nkmin_eta[INDEX2(2, 0, iface) - 0] == 0 || nimin[INDEX2(2, 0, iface) - 0] == 0){\n\
           return ;\n\
        }\n\
        j = 0;\n\
        k = (igll) / (NGLLX);\n\
        i = igll - ((k) * (NGLLX));\n\
        if(k < nkmin_eta[INDEX2(2, 0, iface) - 0] - (1) || k > NGLLX - (1)){\n\
           return ;\n\
        }\n\
        if(i < nimin[INDEX2(2, 0, iface) - 0] - (1) || i > nimax[INDEX2(2, 0, iface) - 0] - (1)){\n\
           return ;\n\
        }\n\
        break;\n\
      case 3 :\n\
        if(nkmin_eta[INDEX2(2, 1, iface) - 0] == 0 || nimin[INDEX2(2, 1, iface) - 0] == 0){\n\
           return ;\n\
        }\n\
        j = NGLLX - (1);\n\
        k = (igll) / (NGLLX);\n\
        i = igll - ((k) * (NGLLX));\n\
        if(k < nkmin_eta[INDEX2(2, 1, iface) - 0] - (1) || k > NGLLX - (1)){\n\
           return ;\n\
        }\n\
        if(i < nimin[INDEX2(2, 1, iface) - 0] - (1) || i > nimax[INDEX2(2, 1, iface) - 0] - (1)){\n\
           return ;\n\
        }\n\
        break;\n\
      }\n\
  }\n\
  iglob = ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] - (1);\n\
  atomicAdd(b_accel + (iglob) * (3) + 0,  -(b_absorb_field[INDEX3(NDIM, NGLL2, 0, igll, iface) - 0]));\n\
  atomicAdd(b_accel + (iglob) * (3) + 1,  -(b_absorb_field[INDEX3(NDIM, NGLL2, 1, igll, iface) - 0]));\n\
  atomicAdd(b_accel + (iglob) * (3) + 2,  -(b_absorb_field[INDEX3(NDIM, NGLL2, 2, igll, iface) - 0]));\n\
}\n\
";
