// from compute_coupling_cuda.cu
#define NDIM 3
#define NGLLX 5
#define INDEX2(xsize,x,y) x + (y)*xsize
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))
#define INDEX5(xsize,ysize,zsize,isize,x,y,z,i,j) x + xsize*(y + ysize*(z + zsize*(i + isize*(j))))

typedef float realw;

__global__ void compute_coupling_fluid_CMB_kernel(realw* displ_crust_mantle,
                                                  realw* accel_outer_core,
                                                  int* ibool_crust_mantle,
                                                  int* ibelm_bottom_crust_mantle,
                                                  realw* normal_top_outer_core,
                                                  realw* jacobian2D_top_outer_core,
                                                  realw* wgllwgll_xy,
                                                  int* ibool_outer_core,
                                                  int* ibelm_top_outer_core,
                                                  int NSPEC2D_TOP_OC) {

  int i = threadIdx.x;
  int j = threadIdx.y;

  int iface = blockIdx.x + gridDim.x*blockIdx.y;

  int k,k_corresp,iglob_cm,iglob_oc,ispec,ispec_selected;
  realw displ_x,displ_y,displ_z,displ_n;
  realw nx,ny,nz;
  realw weight;

  // for surfaces elements exactly at the top of the outer core (crust mantle bottom)
  if( iface < NSPEC2D_TOP_OC ){

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = ibelm_top_outer_core[iface] - 1;
    ispec_selected = ibelm_bottom_crust_mantle[iface] - 1;

    // only for DOFs exactly on the CMB (top of these elements)
    k = NGLLX - 1;
    // get displacement on the solid side using pointwise matching
    k_corresp = 0;

    // get global point number
    // corresponding points are located at the bottom of the mantle
    iglob_cm = ibool_crust_mantle[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k_corresp,ispec_selected)] - 1;

    // elastic displacement on global point
    displ_x = displ_crust_mantle[iglob_cm*3]; // (1,iglob)
    displ_y = displ_crust_mantle[iglob_cm*3+1]; // (2,iglob)
    displ_z = displ_crust_mantle[iglob_cm*3+2]; // (3,iglob)

    // get normal on the CMB
    nx = normal_top_outer_core[INDEX4(NDIM,NGLLX,NGLLX,0,i,j,iface)]; // (1,i,j,iface)
    ny = normal_top_outer_core[INDEX4(NDIM,NGLLX,NGLLX,1,i,j,iface)]; // (2,i,j,iface)
    nz = normal_top_outer_core[INDEX4(NDIM,NGLLX,NGLLX,2,i,j,iface)]; // (3,i,j,iface)

    // calculates displacement component along normal
    // (normal points outwards of acoustic element)
    displ_n = displ_x*nx + displ_y*ny + displ_z*nz;

    // formulation with generalized potential: gets associated, weighted jacobian
    weight = jacobian2D_top_outer_core[INDEX3(NGLLX,NGLLX,i,j,iface)]*wgllwgll_xy[INDEX2(NGLLX,i,j)];

    // get global point number
    iglob_oc = ibool_outer_core[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)] - 1;

    // update fluid acceleration/pressure
    atomicAdd(&accel_outer_core[iglob_oc], + weight*displ_n);
  }
}
