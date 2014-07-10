// from compute_coupling_cuda.cu
#define NDIM 3
#define INDEX2(xsize,x,y) x + (y)*xsize

typedef float realw;

__global__ void compute_coupling_ocean_kernel(realw* accel_crust_mantle,
                                              realw* rmassx_crust_mantle,
                                              realw* rmassy_crust_mantle,
                                              realw* rmassz_crust_mantle,
                                              realw* rmass_ocean_load,
                                              int npoin_ocean_load,
                                              int* ibool_ocean_load,
                                              realw* normal_ocean_load) {

  int ipoin = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;

  int iglob;
  realw nx,ny,nz,rmass;
  realw force_normal_comp;
  realw additional_term_x,additional_term_y,additional_term_z;

  // for global points exactly at the top of the crust mantle (ocean bottom)
  if(ipoin < npoin_ocean_load) {

    // get global point number
    // "-1" from index values to convert from Fortran-> C indexing
    iglob = ibool_ocean_load[ipoin] - 1;

    // get normal
    nx = normal_ocean_load[INDEX2(NDIM,0,ipoin)]; // (1,ipoin)
    ny = normal_ocean_load[INDEX2(NDIM,1,ipoin)]; // (1,ipoin)
    nz = normal_ocean_load[INDEX2(NDIM,2,ipoin)]; // (1,ipoin)

    // make updated component of right-hand side
    // we divide by rmass() which is 1 / M
    // we use the total force which includes the Coriolis term above
    force_normal_comp = accel_crust_mantle[iglob*3]*nx / rmassx_crust_mantle[iglob]
                      + accel_crust_mantle[iglob*3+1]*ny / rmassy_crust_mantle[iglob]
                      + accel_crust_mantle[iglob*3+2]*nz / rmassz_crust_mantle[iglob];

    rmass = rmass_ocean_load[ipoin];

    additional_term_x = (rmass - rmassx_crust_mantle[iglob]) * force_normal_comp;
    additional_term_y = (rmass - rmassy_crust_mantle[iglob]) * force_normal_comp;
    additional_term_z = (rmass - rmassz_crust_mantle[iglob]) * force_normal_comp;

    // since we access this global point only once, no need to use atomics ...
    accel_crust_mantle[iglob*3] += additional_term_x * nx;
    accel_crust_mantle[iglob*3+1] += additional_term_y * ny;
    accel_crust_mantle[iglob*3+2] += additional_term_z * nz;
  }
}
