// from compute_stacey_elastic_cuda.cu
#define NDIM 3
#define NGLLX 5
#define NGLL2 25
#define INDEX2(xsize,x,y) x + (y)*xsize
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))

typedef float realw;

__global__ void compute_stacey_elastic_backward_kernel(realw* b_accel,
                                                       realw* b_absorb_field,
                                                       int interface_type,
                                                       int num_abs_boundary_faces,
                                                       int* abs_boundary_ispec,
                                                       int* nkmin_xi, int* nkmin_eta,
                                                       int* njmin, int* njmax,
                                                       int* nimin, int* nimax,
                                                       int* ibool) {

  int igll = threadIdx.x; // tx
  int iface = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int i,j,k,iglob,ispec;

  // don't compute surface faces outside of range
  // and don't compute points outside NGLLSQUARE==NGLL2==25
  //if(igll < NGLL2 && iface < num_abs_boundary_faces) {

  // way 2: only check face, no further check needed since blocksize = 25
  if( iface < num_abs_boundary_faces){

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    // determines indices i,j,k depending on absorbing boundary type
    switch( interface_type ){
      case 0:
        // xmin
        if( nkmin_xi[INDEX2(2,0,iface)] == 0 || njmin[INDEX2(2,0,iface)] == 0 ) return;

        i = 0; // index -1
        k = (igll/NGLLX);
        j = (igll-k*NGLLX);

        if( k < nkmin_xi[INDEX2(2,0,iface)]-1 || k > NGLLX-1 ) return;
        if( j < njmin[INDEX2(2,0,iface)]-1 || j > NGLLX-1 ) return;

        break;

      case 1:
        // xmax
        if( nkmin_xi[INDEX2(2,1,iface)] == 0 || njmin[INDEX2(2,1,iface)] == 0 ) return;

        i = NGLLX-1;
        k = (igll/NGLLX);
        j = (igll-k*NGLLX);

        if( k < nkmin_xi[INDEX2(2,1,iface)]-1 || k > NGLLX-1 ) return;
        if( j < njmin[INDEX2(2,1,iface)]-1 || j > njmax[INDEX2(2,1,iface)]-1 ) return;

        break;

      case 2:
        // ymin
        if( nkmin_eta[INDEX2(2,0,iface)] == 0 || nimin[INDEX2(2,0,iface)] == 0 ) return;

        j = 0;
        k = (igll/NGLLX);
        i = (igll-k*NGLLX);

        if( k < nkmin_eta[INDEX2(2,0,iface)]-1 || k > NGLLX-1 ) return;
        if( i < nimin[INDEX2(2,0,iface)]-1 || i > nimax[INDEX2(2,0,iface)]-1 ) return;

        break;

      case 3:
        // ymax
        if( nkmin_eta[INDEX2(2,1,iface)] == 0 || nimin[INDEX2(2,1,iface)] == 0 ) return;

        j = NGLLX-1;
        k = (igll/NGLLX);
        i = (igll-k*NGLLX);

        if( k < nkmin_eta[INDEX2(2,1,iface)]-1 || k > NGLLX-1 ) return;
        if( i < nimin[INDEX2(2,1,iface)]-1 || i > nimax[INDEX2(2,1,iface)]-1 ) return;

        break;
    }

    iglob = ibool[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

    atomicAdd(&b_accel[iglob*3  ],-b_absorb_field[INDEX3(NDIM,NGLL2,0,igll,iface)]);
    atomicAdd(&b_accel[iglob*3+1],-b_absorb_field[INDEX3(NDIM,NGLL2,1,igll,iface)]);
    atomicAdd(&b_accel[iglob*3+2],-b_absorb_field[INDEX3(NDIM,NGLL2,2,igll,iface)]);

  } // num_abs_boundary_faces
}

