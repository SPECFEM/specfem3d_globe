const char * compute_stacey_elastic_kernel_program = "\
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
__kernel void compute_stacey_elastic_kernel(const __global float * veloc, __global float * accel, const int interface_type, const int num_abs_boundary_faces, const __global int * abs_boundary_ispec, const __global int * nkmin_xi, const __global int * nkmin_eta, const __global int * njmin, const __global int * njmax, const __global int * nimin, const __global int * nimax, const __global float * abs_boundary_normal, const __global float * abs_boundary_jacobian2D, const __global float * wgllwgll, const __global int * ibool, const __global float * rho_vp, const __global float * rho_vs, const int SAVE_FORWARD, __global float * b_absorb_field){\n\
  int igll;\n\
  int iface;\n\
  int i;\n\
  int j;\n\
  int k;\n\
  int iglob;\n\
  int ispec;\n\
  float vx;\n\
  float vy;\n\
  float vz;\n\
  float vn;\n\
  float nx;\n\
  float ny;\n\
  float nz;\n\
  float rho_vp_temp;\n\
  float rho_vs_temp;\n\
  float tx;\n\
  float ty;\n\
  float tz;\n\
  float jacobianw;\n\
  float fac1;\n\
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
        fac1 = wgllwgll[(k) * (NGLLX) + j - 0];\n\
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
        fac1 = wgllwgll[(k) * (NGLLX) + j - 0];\n\
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
        fac1 = wgllwgll[(k) * (NGLLX) + i - 0];\n\
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
        fac1 = wgllwgll[(k) * (NGLLX) + i - 0];\n\
        break;\n\
      }\n\
  }\n\
  iglob = ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] - (1);\n\
  vx = veloc[(iglob) * (3) + 0 - 0];\n\
  vy = veloc[(iglob) * (3) + 1 - 0];\n\
  vz = veloc[(iglob) * (3) + 2 - 0];\n\
  nx = abs_boundary_normal[INDEX3(NDIM, NGLL2, 0, igll, iface) - 0];\n\
  ny = abs_boundary_normal[INDEX3(NDIM, NGLL2, 1, igll, iface) - 0];\n\
  nz = abs_boundary_normal[INDEX3(NDIM, NGLL2, 2, igll, iface) - 0];\n\
  vn = (vx) * (nx) + (vy) * (ny) + (vz) * (nz);\n\
  rho_vp_temp = rho_vp[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0];\n\
  rho_vs_temp = rho_vs[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0];\n\
  tx = ((rho_vp_temp) * (vn)) * (nx) + (rho_vs_temp) * (vx - ((vn) * (nx)));\n\
  ty = ((rho_vp_temp) * (vn)) * (ny) + (rho_vs_temp) * (vy - ((vn) * (ny)));\n\
  tz = ((rho_vp_temp) * (vn)) * (nz) + (rho_vs_temp) * (vz - ((vn) * (nz)));\n\
  jacobianw = (abs_boundary_jacobian2D[INDEX2(NGLL2, igll, iface) - 0]) * (fac1);\n\
  atomicAdd(accel + (iglob) * (3) + 0, ( -(tx)) * (jacobianw));\n\
  atomicAdd(accel + (iglob) * (3) + 1, ( -(ty)) * (jacobianw));\n\
  atomicAdd(accel + (iglob) * (3) + 2, ( -(tz)) * (jacobianw));\n\
  if(SAVE_FORWARD){\n\
    b_absorb_field[INDEX3(NDIM, NGLL2, 0, igll, iface) - 0] = (tx) * (jacobianw);\n\
    b_absorb_field[INDEX3(NDIM, NGLL2, 1, igll, iface) - 0] = (ty) * (jacobianw);\n\
    b_absorb_field[INDEX3(NDIM, NGLL2, 2, igll, iface) - 0] = (tz) * (jacobianw);\n\
  }\n\
}\n\
";
