// from compute_coupling_cuda.cu
#define NDIM 3
#define NGLLX 5
#define INDEX2(xsize,x,y) x + (y)*xsize
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))
#define INDEX5(xsize,ysize,zsize,isize,x,y,z,i,j) x + xsize*(y + ysize*(z + zsize*(i + isize*(j))))

typedef float realw;

__global__ void compute_coupling_fluid_ICB_kernel(realw* displ_inner_core,
                                                  realw* accel_outer_core,
                                                  int* ibool_inner_core,
                                                  int* ibelm_top_inner_core,
                                                  realw* normal_bottom_outer_core,
                                                  realw* jacobian2D_bottom_outer_core,
                                                  realw* wgllwgll_xy,
                                                  int* ibool_outer_core,
                                                  int* ibelm_bottom_outer_core,
                                                  int NSPEC2D_BOTTOM_OC) {

  int i = threadIdx.x;
  int j = threadIdx.y;

  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int k,k_corresp,iglob_ic,iglob_oc,ispec,ispec_selected;
  realw displ_x,displ_y,displ_z,displ_n;
  realw nx,ny,nz;
  realw weight;

  // for surfaces elements exactly at the bottom of the outer core (inner core top)
  if( iface < NSPEC2D_BOTTOM_OC ){

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = ibelm_bottom_outer_core[iface] - 1;
    ispec_selected = ibelm_top_inner_core[iface] - 1;

    // only for DOFs exactly on the ICB (bottom of these elements)
    k = 0;
    // get displacement on the solid side using pointwise matching
    k_corresp = NGLLX - 1;

    // get global point number
    // corresponding points are located at the bottom of the mantle
    iglob_ic = ibool_inner_core[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k_corresp,ispec_selected)] - 1;

    // elastic displacement on global point
    displ_x = displ_inner_core[iglob_ic*3]; // (1,iglob)
    displ_y = displ_inner_core[iglob_ic*3+1]; // (2,iglob)
    displ_z = displ_inner_core[iglob_ic*3+2]; // (3,iglob)

    // get normal on the ICB
    nx = normal_bottom_outer_core[INDEX4(NDIM,NGLLX,NGLLX,0,i,j,iface)]; // (1,i,j,iface)
    ny = normal_bottom_outer_core[INDEX4(NDIM,NGLLX,NGLLX,1,i,j,iface)]; // (2,i,j,iface)
    nz = normal_bottom_outer_core[INDEX4(NDIM,NGLLX,NGLLX,2,i,j,iface)]; // (3,i,j,iface)

    // calculates displacement component along normal
    // (normal points outwards of acoustic element)
    displ_n = displ_x*nx + displ_y*ny + displ_z*nz;

    // formulation with generalized potential: gets associated, weighted jacobian
    weight = jacobian2D_bottom_outer_core[INDEX3(NGLLX,NGLLX,i,j,iface)]*wgllwgll_xy[INDEX2(NGLLX,i,j)];

    // get global point number
    iglob_oc = ibool_outer_core[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] - 1;

    // update fluid acceleration/pressure
    // note: sign changes to minus because of normal pointing down into inner core
    atomicAdd(&accel_outer_core[iglob_oc], - weight*displ_n);
  }
}
