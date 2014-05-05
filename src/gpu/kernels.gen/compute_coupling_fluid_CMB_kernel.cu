#ifndef INDEX2
#define INDEX2(xsize,x,y) x + (y)*xsize
#endif
#ifndef INDEX3
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)
#endif
#ifndef INDEX4
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))
#endif
#ifndef INDEX5
#define INDEX5(xsize,ysize,zsize,isize,x,y,z,i,j) x + xsize*(y + ysize*(z + zsize*(i + isize*(j))))
#endif
#ifndef NDIM
#define NDIM 3
#endif
#ifndef NGLLX
#define NGLLX 5
#endif
#ifndef NGLL2
#define NGLL2 25
#endif
#ifndef NGLL3
#define NGLL3 125
#endif
#ifndef NGLL3_PADDED
#define NGLL3_PADDED 128
#endif
#ifndef N_SLS
#define N_SLS 3
#endif
#ifndef IREGION_CRUST_MANTLE
#define IREGION_CRUST_MANTLE 1
#endif
#ifndef IREGION_INNER_CORE
#define IREGION_INNER_CORE 3
#endif
#ifndef IFLAG_IN_FICTITIOUS_CUBE
#define IFLAG_IN_FICTITIOUS_CUBE 11
#endif
#ifndef R_EARTH_KM
#define R_EARTH_KM 6371.0f
#endif
#ifndef COLORING_MIN_NSPEC_INNER_CORE
#define COLORING_MIN_NSPEC_INNER_CORE 1000
#endif
#ifndef COLORING_MIN_NSPEC_OUTER_CORE
#define COLORING_MIN_NSPEC_OUTER_CORE 1000
#endif
#ifndef BLOCKSIZE_TRANSFER
#define BLOCKSIZE_TRANSFER 256
#endif
__global__ void compute_coupling_fluid_CMB_kernel(const float * displ_crust_mantle, float * accel_outer_core, const int * ibool_crust_mantle, const int * ibelm_bottom_crust_mantle, const float * normal_top_outer_core, const float * jacobian2D_top_outer_core, const float * wgllwgll_xy, const int * ibool_outer_core, const int * ibelm_top_outer_core, const int NSPEC2D_TOP_OC){
  int i;
  int j;
  int k;
  int iface;
  int k_corresp;
  float displ_n;
  int iglob_cm;
  int iglob_oc;
  int ispec;
  int ispec_selected;
  float displ_x;
  float displ_y;
  float displ_z;
  float nx;
  float ny;
  float nz;
  float weight;
  i = threadIdx.x;
  j = threadIdx.y;
  iface = blockIdx.x + (gridDim.x) * (blockIdx.y);
  if(iface < NSPEC2D_TOP_OC){
    ispec = ibelm_top_outer_core[iface - 0] - (1);
    ispec_selected = ibelm_bottom_crust_mantle[iface - 0] - (1);
    k = NGLLX - (1);
    k_corresp = 0;
    iglob_cm = ibool_crust_mantle[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k_corresp, ispec_selected) - 0] - (1);
    displ_x = displ_crust_mantle[(iglob_cm) * (3) + 0 - 0];
    displ_y = displ_crust_mantle[(iglob_cm) * (3) + 1 - 0];
    displ_z = displ_crust_mantle[(iglob_cm) * (3) + 2 - 0];
    nx = normal_top_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 0, i, j, iface) - 0];
    ny = normal_top_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 1, i, j, iface) - 0];
    nz = normal_top_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 2, i, j, iface) - 0];
    displ_n = (displ_x) * (nx) + (displ_y) * (ny) + (displ_z) * (nz);
    weight = (jacobian2D_top_outer_core[INDEX3(NGLLX, NGLLX, i, j, iface) - 0]) * (wgllwgll_xy[INDEX2(NGLLX, i, j) - 0]);
    iglob_oc = ibool_outer_core[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] - (1);
    atomicAdd(accel_outer_core + iglob_oc, (weight) * (displ_n));
  }
}
