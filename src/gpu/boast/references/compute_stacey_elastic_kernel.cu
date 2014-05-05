// from compute_stacey_elastic_cuda.cu
#define NDIM 3
#define NGLLX 5
#define NGLL2 25
#define INDEX2(xsize,x,y) x + (y)*xsize
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))

typedef float realw;

__global__ void compute_stacey_elastic_kernel(realw* veloc,
                                              realw* accel,
                                              int interface_type,
                                              int num_abs_boundary_faces,
                                              int* abs_boundary_ispec,
                                              int* nkmin_xi, int* nkmin_eta,
                                              int* njmin, int* njmax,
                                              int* nimin, int* nimax,
                                              realw* abs_boundary_normal,
                                              realw* abs_boundary_jacobian2D,
                                              realw* wgllwgll,
                                              int* ibool,
                                              realw* rho_vp,
                                              realw* rho_vs,
                                              int SAVE_FORWARD,
                                              realw* b_absorb_field) {

  int igll = threadIdx.x; // tx
  int iface = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int i,j,k,iglob,ispec;
  realw vx,vy,vz,vn;
  realw nx,ny,nz;
  realw rho_vp_temp,rho_vs_temp;
  realw tx,ty,tz;
  realw jacobianw;
  realw fac1;

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

        fac1 = wgllwgll[k*NGLLX+j];
        break;

      case 1:
        // xmax
        if( nkmin_xi[INDEX2(2,1,iface)] == 0 || njmin[INDEX2(2,1,iface)] == 0 ) return;

        i = NGLLX-1;
        k = (igll/NGLLX);
        j = (igll-k*NGLLX);

        if( k < nkmin_xi[INDEX2(2,1,iface)]-1 || k > NGLLX-1 ) return;
        if( j < njmin[INDEX2(2,1,iface)]-1 || j > njmax[INDEX2(2,1,iface)]-1 ) return;

        fac1 = wgllwgll[k*NGLLX+j];
        break;

      case 2:
        // ymin
        if( nkmin_eta[INDEX2(2,0,iface)] == 0 || nimin[INDEX2(2,0,iface)] == 0 ) return;

        j = 0;
        k = (igll/NGLLX);
        i = (igll-k*NGLLX);

        if( k < nkmin_eta[INDEX2(2,0,iface)]-1 || k > NGLLX-1 ) return;
        if( i < nimin[INDEX2(2,0,iface)]-1 || i > nimax[INDEX2(2,0,iface)]-1 ) return;

        fac1 = wgllwgll[k*NGLLX+i];
        break;

      case 3:
        // ymax
        if( nkmin_eta[INDEX2(2,1,iface)] == 0 || nimin[INDEX2(2,1,iface)] == 0 ) return;

        j = NGLLX-1;
        k = (igll/NGLLX);
        i = (igll-k*NGLLX);

        if( k < nkmin_eta[INDEX2(2,1,iface)]-1 || k > NGLLX-1 ) return;
        if( i < nimin[INDEX2(2,1,iface)]-1 || i > nimax[INDEX2(2,1,iface)]-1 ) return;

        fac1 = wgllwgll[k*NGLLX+i];
        break;
    }

    iglob = ibool[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

    // gets associated velocity
    vx = veloc[iglob*3+0];
    vy = veloc[iglob*3+1];
    vz = veloc[iglob*3+2];

    // gets associated normal
    nx = abs_boundary_normal[INDEX3(NDIM,NGLL2,0,igll,iface)];
    ny = abs_boundary_normal[INDEX3(NDIM,NGLL2,1,igll,iface)];
    nz = abs_boundary_normal[INDEX3(NDIM,NGLL2,2,igll,iface)];

    // // velocity component in normal direction (normal points out of element)
    vn = vx*nx + vy*ny + vz*nz;

    rho_vp_temp = rho_vp[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];
    rho_vs_temp = rho_vs[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];

    tx = rho_vp_temp*vn*nx + rho_vs_temp*(vx-vn*nx);
    ty = rho_vp_temp*vn*ny + rho_vs_temp*(vy-vn*ny);
    tz = rho_vp_temp*vn*nz + rho_vs_temp*(vz-vn*nz);

    jacobianw = abs_boundary_jacobian2D[INDEX2(NGLL2,igll,iface)]*fac1;

    atomicAdd(&accel[iglob*3],-tx*jacobianw);
    atomicAdd(&accel[iglob*3+1],-ty*jacobianw);
    atomicAdd(&accel[iglob*3+2],-tz*jacobianw);

    if( SAVE_FORWARD ){
      b_absorb_field[INDEX3(NDIM,NGLL2,0,igll,iface)] = tx*jacobianw;
      b_absorb_field[INDEX3(NDIM,NGLL2,1,igll,iface)] = ty*jacobianw;
      b_absorb_field[INDEX3(NDIM,NGLL2,2,igll,iface)] = tz*jacobianw;
    } // SIMULATION_TYPE

  } // num_abs_boundary_faces
}

