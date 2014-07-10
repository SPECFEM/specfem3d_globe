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
__global__ void compute_coupling_ICB_fluid_kernel(const float * displ_inner_core, float * accel_inner_core, const float * accel_outer_core, const int * ibool_inner_core, const int * ibelm_top_inner_core, const float * normal_bottom_outer_core, const float * jacobian2D_bottom_outer_core, const float * wgllwgll_xy, const int * ibool_outer_core, const int * ibelm_bottom_outer_core, const float RHO_BOTTOM_OC, const float minus_g_icb, int GRAVITY, const int NSPEC2D_TOP_IC){
  int i;
  int j;
  int k;
  int iface;
  int k_corresp;
  int iglob_oc;
  int iglob_ic;
  float pressure;
  int ispec;
  int ispec_selected;
  float nx;
  float ny;
  float nz;
  float weight;
  i = threadIdx.x;
  j = threadIdx.y;
  iface = blockIdx.x + (gridDim.x) * (blockIdx.y);
  if(iface < NSPEC2D_TOP_IC){
    ispec = ibelm_top_inner_core[iface - 0] - (1);
    ispec_selected = ibelm_bottom_outer_core[iface - 0] - (1);
    k = NGLLX - (1);
    k_corresp = 0;
    iglob_oc = ibool_outer_core[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k_corresp, ispec_selected) - 0] - (1);
    nx = normal_bottom_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 0, i, j, iface) - 0];
    ny = normal_bottom_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 1, i, j, iface) - 0];
    nz = normal_bottom_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 2, i, j, iface) - 0];
    weight = (jacobian2D_bottom_outer_core[INDEX3(NGLLX, NGLLX, i, j, iface) - 0]) * (wgllwgll_xy[INDEX2(NGLLX, i, j) - 0]);
    iglob_ic = ibool_inner_core[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] - (1);
    if(GRAVITY){
      pressure = (RHO_BOTTOM_OC) * ((minus_g_icb) * ((displ_inner_core[(iglob_ic) * (3) - 0]) * (nx) + (displ_inner_core[(iglob_ic) * (3) + 1 - 0]) * (ny) + (displ_inner_core[(iglob_ic) * (3) + 2 - 0]) * (nz)) - (accel_outer_core[iglob_oc - 0]));
    } else {
      pressure = ( -(RHO_BOTTOM_OC)) * (accel_outer_core[iglob_oc - 0]);
    }
    atomicAdd(accel_inner_core + (iglob_ic) * (3) + 0, (( -(weight)) * (nx)) * (pressure));
    atomicAdd(accel_inner_core + (iglob_ic) * (3) + 1, (( -(weight)) * (ny)) * (pressure));
    atomicAdd(accel_inner_core + (iglob_ic) * (3) + 2, (( -(weight)) * (nz)) * (pressure));
  }
}
