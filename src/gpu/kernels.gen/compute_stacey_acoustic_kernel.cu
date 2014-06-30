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
__global__ void compute_stacey_acoustic_kernel(const float * potential_dot_acoustic, float * potential_dot_dot_acoustic, const int interface_type, const int num_abs_boundary_faces, const int * abs_boundary_ispec, const int * nkmin_xi, const int * nkmin_eta, const int * njmin, const int * njmax, const int * nimin, const int * nimax, const float * abs_boundary_jacobian2D, const float * wgllwgll, const int * ibool, const float * vpstore, const int SAVE_FORWARD, float * b_absorb_potential){
  int igll;
  int iface;
  int i;
  int j;
  int k;
  int iglob;
  int ispec;
  float sn;
  float jacobianw;
  float fac1;
  igll = threadIdx.x;
  iface = blockIdx.x + (blockIdx.y) * (gridDim.x);
  if(iface < num_abs_boundary_faces){
    ispec = abs_boundary_ispec[iface - (0)] - (1);
    switch(interface_type){
      case 4 :
        if(nkmin_xi[INDEX2(2, 0, iface) - (0)] == 0 || njmin[INDEX2(2, 0, iface) - (0)] == 0){
           return ;
        }
        i = 0;
        k = (igll) / (NGLLX);
        j = igll - ((k) * (NGLLX));
        if(k < nkmin_xi[INDEX2(2, 0, iface) - (0)] - (1) || k > NGLLX - (1)){
           return ;
        }
        if(j < njmin[INDEX2(2, 0, iface) - (0)] - (1) || j > njmax[INDEX2(2, 0, iface) - (0)] - (1)){
           return ;
        }
        fac1 = wgllwgll[(k) * (NGLLX) + j - (0)];
        break;
      case 5 :
        if(nkmin_xi[INDEX2(2, 1, iface) - (0)] == 0 || njmin[INDEX2(2, 1, iface) - (0)] == 0){
           return ;
        }
        i = NGLLX - (1);
        k = (igll) / (NGLLX);
        j = igll - ((k) * (NGLLX));
        if(k < nkmin_xi[INDEX2(2, 1, iface) - (0)] - (1) || k > NGLLX - (1)){
           return ;
        }
        if(j < njmin[INDEX2(2, 1, iface) - (0)] - (1) || j > njmax[INDEX2(2, 1, iface) - (0)] - (1)){
           return ;
        }
        fac1 = wgllwgll[(k) * (NGLLX) + j - (0)];
        break;
      case 6 :
        if(nkmin_eta[INDEX2(2, 0, iface) - (0)] == 0 || nimin[INDEX2(2, 0, iface) - (0)] == 0){
           return ;
        }
        j = 0;
        k = (igll) / (NGLLX);
        i = igll - ((k) * (NGLLX));
        if(k < nkmin_eta[INDEX2(2, 0, iface) - (0)] - (1) || k > NGLLX - (1)){
           return ;
        }
        if(i < nimin[INDEX2(2, 0, iface) - (0)] - (1) || i > nimax[INDEX2(2, 0, iface) - (0)] - (1)){
           return ;
        }
        fac1 = wgllwgll[(k) * (NGLLX) + i - (0)];
        break;
      case 7 :
        if(nkmin_eta[INDEX2(2, 1, iface) - (0)] == 0 || nimin[INDEX2(2, 1, iface) - (0)] == 0){
           return ;
        }
        j = NGLLX - (1);
        k = (igll) / (NGLLX);
        i = igll - ((k) * (NGLLX));
        if(k < nkmin_eta[INDEX2(2, 1, iface) - (0)] - (1) || k > NGLLX - (1)){
           return ;
        }
        if(i < nimin[INDEX2(2, 1, iface) - (0)] - (1) || i > nimax[INDEX2(2, 1, iface) - (0)] - (1)){
           return ;
        }
        fac1 = wgllwgll[(k) * (NGLLX) + i - (0)];
        break;
      case 8 :
        k = 0;
        j = (igll) / (NGLLX);
        i = igll - ((j) * (NGLLX));
        if(j < 0 || j > NGLLX - (1)){
           return ;
        }
        if(i < 0 || i > NGLLX - (1)){
           return ;
        }
        fac1 = wgllwgll[(j) * (NGLLX) + i - (0)];
        break;
      }
  }
  iglob = ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - (0)] - (1);
  sn = (potential_dot_acoustic[iglob - (0)]) / (vpstore[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - (0)]);
  jacobianw = (abs_boundary_jacobian2D[INDEX2(NGLL2, igll, iface) - (0)]) * (fac1);
  atomicAdd(potential_dot_dot_acoustic + iglob, ( -(sn)) * (jacobianw));
  if(SAVE_FORWARD){
    b_absorb_potential[INDEX2(NGLL2, igll, iface) - (0)] = (sn) * (jacobianw);
  }
}
