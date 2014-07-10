#ifndef INDEX2
#define INDEX2(isize,i,j) i + isize*j
#endif
#ifndef INDEX3
#define INDEX3(isize,jsize,i,j,k) i + isize*(j + jsize*k)
#endif
#ifndef INDEX4
#define INDEX4(isize,jsize,ksize,i,j,k,x) i + isize*(j + jsize*(k + ksize*x))
#endif
#ifndef INDEX5
#define INDEX5(isize,jsize,ksize,xsize,i,j,k,x,y) i + isize*(j + jsize*(k + ksize*(x + xsize*y)))
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
__global__ void compute_coupling_fluid_ICB_kernel(const float * displ_inner_core, float * accel_outer_core, const int * ibool_inner_core, const int * ibelm_top_inner_core, const float * normal_bottom_outer_core, const float * jacobian2D_bottom_outer_core, const float * wgllwgll_xy, const int * ibool_outer_core, const int * ibelm_bottom_outer_core, const int NSPEC2D_BOTTOM_OC){
  int i;
  int j;
  int k;
  int iface;
  int k_corresp;
  float displ_n;
  int iglob_ic;
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
  if(iface < NSPEC2D_BOTTOM_OC){
    ispec = ibelm_bottom_outer_core[iface - 0] - (1);
    ispec_selected = ibelm_top_inner_core[iface - 0] - (1);
    k = 0;
    k_corresp = NGLLX - (1);
    iglob_ic = ibool_inner_core[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k_corresp, ispec_selected) - 0] - (1);
    displ_x = displ_inner_core[(iglob_ic) * (3) + 0 - 0];
    displ_y = displ_inner_core[(iglob_ic) * (3) + 1 - 0];
    displ_z = displ_inner_core[(iglob_ic) * (3) + 2 - 0];
    nx = normal_bottom_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 0, i, j, iface) - 0];
    ny = normal_bottom_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 1, i, j, iface) - 0];
    nz = normal_bottom_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 2, i, j, iface) - 0];
    displ_n = (displ_x) * (nx) + (displ_y) * (ny) + (displ_z) * (nz);
    weight = (jacobian2D_bottom_outer_core[INDEX3(NGLLX, NGLLX, i, j, iface) - 0]) * (wgllwgll_xy[INDEX2(NGLLX, i, j) - 0]);
    iglob_oc = ibool_outer_core[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] - (1);
    atomicAdd(accel_outer_core + iglob_oc, ( -(weight)) * (displ_n));
  }
}
