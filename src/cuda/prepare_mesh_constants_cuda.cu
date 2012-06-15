/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  2 . 0
 !               ---------------------------------------
 !
 !          Main authors: Dimitri Komatitsch and Jeroen Tromp
 !    Princeton University, USA and University of Pau / CNRS / INRIA
 ! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
 !                            April 2011
 !
 ! This program is free software; you can redistribute it and/or modify
 ! it under the terms of the GNU General Public License as published by
 ! the Free Software Foundation; either version 2 of the License, or
 ! (at your option) any later version.
 !
 ! This program is distributed in the hope that it will be useful,
 ! but WITHOUT ANY WARRANTY; without even the implied warranty of
 ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ! GNU General Public License for more details.
 !
 ! You should have received a copy of the GNU General Public License along
 ! with this program; if not, write to the Free Software Foundation, Inc.,
 ! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 !
 !=====================================================================
 */

#include <stdio.h>
#include <cuda.h>
#include <cublas.h>

#ifdef WITH_MPI
#include <mpi.h>
#endif

#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"
#include "prepare_constants_cuda.h"



/* ----------------------------------------------------------------------------------------------- */

// GPU preparation

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_cuda_device,
              PREPARE_CUDA_DEVICE)(int* myrank_f,int* ncuda_devices) {
  TRACE("prepare_cuda_device");

  // Gets rank number of MPI process
  int myrank = *myrank_f;

  // cuda initialization (needs -lcuda library)
  // note:   cuInit initializes the driver API.
  //             it is needed for any following CUDA driver API function call (format cuFUNCTION(..) )
  //             however, for the CUDA runtime API functions (format cudaFUNCTION(..) )
  //             the initialization is implicit, thus cuInit() here would not be needed...
  CUresult status = cuInit(0);
  if ( CUDA_SUCCESS != status ) exit_on_error("CUDA driver API device initialization failed\n");

  // returns a handle to the first cuda compute device
  CUdevice dev;
  status = cuDeviceGet(&dev, 0);
  if ( CUDA_SUCCESS != status ) exit_on_error("CUDA device not found\n");

  // gets device properties
  int major,minor;
  status = cuDeviceComputeCapability(&major,&minor,dev);
  if ( CUDA_SUCCESS != status ) exit_on_error("CUDA device information not found\n");

  // make sure that the device has compute capability >= 1.3
  if (major < 1){
    fprintf(stderr,"Compute capability major number should be at least 1, got: %d \nexiting...\n",major);
    exit_on_error("CUDA Compute capability major number should be at least 1\n");
  }
  if (major == 1 && minor < 3){
    fprintf(stderr,"Compute capability should be at least 1.3, got: %d.%d \nexiting...\n",major,minor);
    exit_on_error("CUDA Compute capability major number should be at least 1.3\n");
  }

  // note: from here on we use the runtime API  ...
  // Gets number of GPU devices
  int device_count = 0;
  cudaGetDeviceCount(&device_count);
  exit_on_cuda_error("CUDA runtime cudaGetDeviceCount: check if driver and runtime libraries work together\nexiting...\n");

  // returns device count to fortran
  if (device_count == 0) exit_on_error("CUDA runtime error: there is no device supporting CUDA\n");
  *ncuda_devices = device_count;


  // Sets the active device
  if(device_count > 1) {
    // generalized for more GPUs per node
    // note: without previous context release, cudaSetDevice will complain with the cuda error
    //         "setting the device when a process is active is not allowed"
    // releases previous contexts
    cudaThreadExit();

    //printf("rank %d: cuda device count = %d sets device = %d \n",myrank,device_count,myrank % device_count);
    //MPI_Barrier(MPI_COMM_WORLD);

    // sets active device
    cudaSetDevice( myrank % device_count );
    exit_on_cuda_error("cudaSetDevice");
  }

  // returns a handle to the active device
  int device;
  cudaGetDevice(&device);

  // get device properties
  struct cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp,device);

  // exit if the machine has no CUDA-enabled device
  if (deviceProp.major == 9999 && deviceProp.minor == 9999){
    fprintf(stderr,"No CUDA-enabled device found, exiting...\n\n");
    exit_on_error("CUDA runtime error: there is no CUDA-enabled device found\n");
  }

  // outputs device infos to file
  char filename[BUFSIZ];
  FILE* fp;
  sprintf(filename,"OUTPUT_FILES/gpu_device_info_proc_%06d.txt",myrank);
  fp = fopen(filename,"a+");
  if (fp != NULL){
    // display device properties
    fprintf(fp,"Device Name = %s\n",deviceProp.name);
    fprintf(fp,"multiProcessorCount: %d\n",deviceProp.multiProcessorCount);
    fprintf(fp,"totalGlobalMem (in MB): %f\n",(unsigned long) deviceProp.totalGlobalMem / (1024.f * 1024.f));
    fprintf(fp,"totalGlobalMem (in GB): %f\n",(unsigned long) deviceProp.totalGlobalMem / (1024.f * 1024.f * 1024.f));
    fprintf(fp,"sharedMemPerBlock (in bytes): %lu\n",(unsigned long) deviceProp.sharedMemPerBlock);
    fprintf(fp,"Maximum number of threads per block: %d\n",deviceProp.maxThreadsPerBlock);
    fprintf(fp,"Maximum size of each dimension of a block: %d x %d x %d\n",
            deviceProp.maxThreadsDim[0],deviceProp.maxThreadsDim[1],deviceProp.maxThreadsDim[2]);
    fprintf(fp,"Maximum sizes of each dimension of a grid: %d x %d x %d\n",
            deviceProp.maxGridSize[0],deviceProp.maxGridSize[1],deviceProp.maxGridSize[2]);
    fprintf(fp,"Compute capability of the device = %d.%d\n", deviceProp.major, deviceProp.minor);
    if(deviceProp.canMapHostMemory){
      fprintf(fp,"canMapHostMemory: TRUE\n");
    }else{
      fprintf(fp,"canMapHostMemory: FALSE\n");
    }
    if(deviceProp.deviceOverlap){
      fprintf(fp,"deviceOverlap: TRUE\n");
    }else{
      fprintf(fp,"deviceOverlap: FALSE\n");
    }

    // make sure that the device has compute capability >= 1.3
    //if (deviceProp.major < 1){
    //  fprintf(stderr,"Compute capability major number should be at least 1, exiting...\n\n");
    //  exit_on_error("CUDA Compute capability major number should be at least 1");
    //}
    //if (deviceProp.major == 1 && deviceProp.minor < 3){
    //  fprintf(stderr,"Compute capability should be at least 1.3, exiting...\n");
    //  exit_on_error("CUDA Compute capability major number should be at least 1.3");
    //}

    // outputs initial memory infos via cudaMemGetInfo()
    double free_db,used_db,total_db;
    get_free_memory(&free_db,&used_db,&total_db);
    fprintf(fp,"%d: GPU memory usage: used = %f MB, free = %f MB, total = %f MB\n",myrank,
            used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);

    fclose(fp);
  }
}

/* ----------------------------------------------------------------------------------------------- */

// SIMULATION constants

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_constants_device,
              PREPARE_CONSTANTS_DEVICE)(long* Mesh_pointer,
                                        int* myrank_f,
                                        int* h_NGLLX,
                                        realw* h_hprime_xx,realw* h_hprime_yy,realw* h_hprime_zz,
                                        realw* h_hprimewgll_xx,realw* h_hprimewgll_yy,realw* h_hprimewgll_zz,
                                        realw* h_wgllwgll_xy,realw* h_wgllwgll_xz,realw* h_wgllwgll_yz,
                                        int* NSOURCES,int* nsources_local,
                                        realw* h_sourcearrays,
                                        int* h_islice_selected_source,int* h_ispec_selected_source,
                                        int* h_number_receiver_global,
                                        int* h_islice_selected_rec,int* h_ispec_selected_rec,
                                        int* nrec,int* nrec_local, int* nadj_rec_local,
                                        int* NSPEC_CRUST_MANTLE, int* NGLOB_CRUST_MANTLE,
                                        int* NSPEC_CRUST_MANTLE_STRAIN_ONLY,
                                        int* NGLOB_CRUST_MANTLE_OCEANS,
                                        int* NSPEC_OUTER_CORE, int* NGLOB_OUTER_CORE,
                                        int* NSPEC_INNER_CORE, int* NGLOB_INNER_CORE,
                                        int* SIMULATION_TYPE,
                                        int* NOISE_TOMOGRAPHY,
                                        int* SAVE_FORWARD_f,
                                        int* ABSORBING_CONDITIONS_f,
                                        int* OCEANS_f,
                                        int* GRAVITY_f,
                                        int* ROTATION_f,
                                        int* ATTENUATION_f,
                                        int* ATTENUATION_NEW_f,
                                        int* USE_ATTENUATION_MIMIC_f,
                                        int* COMPUTE_AND_STORE_STRAIN_f,
                                        int* ANISOTROPIC_3D_MANTLE_f,
                                        int* ANISOTROPIC_INNER_CORE_f,
                                        int* SAVE_BOUNDARY_MESH_f,
                                        int* USE_MESH_COLORING_GPU_f,
                                        int* ANISOTROPIC_KL_f,
                                        int* APPROXIMATE_HESS_KL_f) {

TRACE("prepare_constants_device");

  // allocates mesh parameter structure
  Mesh* mp = (Mesh*) malloc( sizeof(Mesh) );
  if (mp == NULL) exit_on_error("error allocating mesh pointer");
  *Mesh_pointer = (long)mp;

  // checks if NGLLX == 5
  if( *h_NGLLX != NGLLX ){
    exit_on_error("NGLLX must be 5 for CUDA devices");
  }

  // sets constant arrays
  setConst_hprime_xx(h_hprime_xx,mp);
  setConst_hprime_yy(h_hprime_yy,mp);
  setConst_hprime_zz(h_hprime_zz,mp);
  setConst_hprimewgll_xx(h_hprimewgll_xx,mp);
  setConst_hprimewgll_yy(h_hprimewgll_yy,mp);
  setConst_hprimewgll_zz(h_hprimewgll_zz,mp);
  setConst_wgllwgll_xy(h_wgllwgll_xy,mp);
  setConst_wgllwgll_xz(h_wgllwgll_xz,mp);
  setConst_wgllwgll_yz(h_wgllwgll_yz,mp);

  // sets global parameters
  mp->NSPEC_CRUST_MANTLE = *NSPEC_CRUST_MANTLE;
  mp->NGLOB_CRUST_MANTLE = *NGLOB_CRUST_MANTLE;
  mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY = *NSPEC_CRUST_MANTLE_STRAIN_ONLY;
  mp->NGLOB_CRUST_MANTLE_OCEANS = *NGLOB_CRUST_MANTLE_OCEANS;
  mp->NSPEC_OUTER_CORE = *NSPEC_OUTER_CORE;
  mp->NGLOB_OUTER_CORE = *NGLOB_OUTER_CORE;
  mp->NSPEC_INNER_CORE = *NSPEC_INNER_CORE;
  mp->NGLOB_INNER_CORE = *NGLOB_INNER_CORE;

  // simulation type
  mp->simulation_type = *SIMULATION_TYPE;
  mp->noise_tomography = *NOISE_TOMOGRAPHY;

  // simulation flags initialization
  mp->save_forward = *SAVE_FORWARD_f;
  mp->absorbing_conditions = *ABSORBING_CONDITIONS_f;
  mp->oceans = *OCEANS_f;
  mp->gravity = *GRAVITY_f;
  mp->rotation = *ROTATION_f;
  mp->attenuation = *ATTENUATION_f;
  mp->attenuation_new = *ATTENUATION_NEW_f;
  mp->use_attenuation_mimic = *USE_ATTENUATION_MIMIC_f;
  mp->compute_and_store_strain = *COMPUTE_AND_STORE_STRAIN_f;
  mp->anisotropic_3D_mantle = *ANISOTROPIC_3D_MANTLE_f;
  mp->anisotropic_inner_core = *ANISOTROPIC_INNER_CORE_f;
  mp->save_boundary_mesh = *SAVE_BOUNDARY_MESH_f;

  mp->anisotropic_kl = *ANISOTROPIC_KL_f;
  mp->approximate_hess_kl = *APPROXIMATE_HESS_KL_f;

  // mpi process rank
  mp->myrank = *myrank_f;

  // mesh coloring flag
#ifdef USE_MESH_COLORING_GPU
  mp->use_mesh_coloring_gpu = 1;
  if( ! *USE_MESH_COLORING_GPU_f ){exit_on_error("error with USE_MESH_COLORING_GPU constant; please re-compile\n");}
#else
  // mesh coloring
  // note: this here passes the coloring as an option to the kernel routines
  //          the performance seems to be the same if one uses the pre-processing directives above or not
  mp->use_mesh_coloring_gpu = *USE_MESH_COLORING_GPU_f;
#endif

  // sources
  mp->nsources_local = *nsources_local;
  if( mp->simulation_type == 1  || mp->simulation_type == 3 ){
    // not needed in case of pure adjoint simulations (SIMULATION_TYPE == 2)
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_sourcearrays,
                                       sizeof(realw)* *NSOURCES*3*NGLL3),1301);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_sourcearrays, h_sourcearrays,
                                       sizeof(realw)* *NSOURCES*3*NGLL3,cudaMemcpyHostToDevice),1302);
    // buffer for source time function values
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_stf_pre_compute,
                                       *NSOURCES*sizeof(double)),1303);
  }

  print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_islice_selected_source,
                                     sizeof(int) * *NSOURCES),1401);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_islice_selected_source, h_islice_selected_source,
                                     sizeof(int)* *NSOURCES,cudaMemcpyHostToDevice),1402);

  print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_ispec_selected_source,
                                     sizeof(int)* *NSOURCES),1403);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_ispec_selected_source, h_ispec_selected_source,
                                     sizeof(int)* *NSOURCES,cudaMemcpyHostToDevice),1404);


  // receiver stations
  // note that:   size(number_receiver_global) = nrec_local
  //                   size(ispec_selected_rec) = nrec
  // number of receiver located in this partition
  mp->nrec_local = *nrec_local;
  if( mp->nrec_local > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_number_receiver_global),mp->nrec_local*sizeof(int)),1);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_number_receiver_global,h_number_receiver_global,
                                     mp->nrec_local*sizeof(int),cudaMemcpyHostToDevice),1512);

    // for seismograms
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_station_seismo_field),
                                       3*NGLL3*(mp->nrec_local)*sizeof(realw)),4015);

    mp->h_station_seismo_field = (realw*) malloc( 3*NGLL3*(mp->nrec_local)*sizeof(realw) );
    if( mp->h_station_seismo_field == NULL) exit_on_error("h_station_seismo_field not allocated \n");

  }
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_ispec_selected_rec),(*nrec)*sizeof(int)),1513);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_ispec_selected_rec,h_ispec_selected_rec,
                                     (*nrec)*sizeof(int),cudaMemcpyHostToDevice),1514);

  // receiver adjoint source arrays only used for noise and adjoint simulations
  // adjoint source arrays
  mp->nadj_rec_local = *nadj_rec_local;
  if( mp->nadj_rec_local > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_adj_sourcearrays,
                                       (mp->nadj_rec_local)*3*NGLL3*sizeof(realw)),6003);

    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_pre_computed_irec,
                                       (mp->nadj_rec_local)*sizeof(int)),6004);

    // prepares local irec array:
    // the irec_local variable needs to be precomputed (as
    // h_pre_comp..), because normally it is in the loop updating accel,
    // and due to how it's incremented, it cannot be parallelized
    int* h_pre_computed_irec = (int*) malloc( (mp->nadj_rec_local)*sizeof(int) );
    if( h_pre_computed_irec == NULL ) exit_on_error("h_pre_computed_irec not allocated\n");

    int irec_local = 0;
    for(int irec = 0; irec < *nrec; irec++) {
      if(mp->myrank == h_islice_selected_rec[irec]) {
        irec_local++;
        h_pre_computed_irec[irec_local-1] = irec;
      }
    }
    if( irec_local != mp->nadj_rec_local ) exit_on_error("prepare_sim2_or_3_const_device: irec_local not equal\n");
    // copies values onto GPU
    print_CUDA_error_if_any(cudaMemcpy(mp->d_pre_computed_irec,h_pre_computed_irec,
                                       (mp->nadj_rec_local)*sizeof(int),cudaMemcpyHostToDevice),6010);
    free(h_pre_computed_irec);

    // temporary array to prepare extracted source array values
    mp->h_adj_sourcearrays_slice = (realw*) malloc( (mp->nadj_rec_local)*3*NGLL3*sizeof(realw) );
    if( mp->h_adj_sourcearrays_slice == NULL ) exit_on_error("h_adj_sourcearrays_slice not allocated\n");
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_constants_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// ROTATION simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_rotation_device,
              PREPARE_FIELDS_ROTATION_DEVICE)(long* Mesh_pointer_f,
                                              realw* two_omega_earth,
                                              realw* deltat,
                                              realw* A_array_rotation,
                                              realw* B_array_rotation,
                                              realw* b_two_omega_earth,
                                              realw* b_deltat,
                                              realw* b_A_array_rotation,
                                              realw* b_B_array_rotation,
                                              int* NSPEC_OUTER_CORE_ROTATION
                                              ) {

  TRACE("prepare_fields_rotation_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // arrays only needed when rotation is required
  if( ! mp->rotation ){ exit_on_cuda_error("prepare_fields_rotation_device rotation not properly initialized"); }

  // rotation arrays (needed only for outer core region)
  mp->d_two_omega_earth = *two_omega_earth;
  mp->d_deltat = *deltat;

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_A_array_rotation,
                                     NGLL3*(*NSPEC_OUTER_CORE_ROTATION)*sizeof(realw)),9000);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_A_array_rotation, A_array_rotation,
                                     NGLL3*(*NSPEC_OUTER_CORE_ROTATION)*sizeof(realw),cudaMemcpyHostToDevice),9001);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_B_array_rotation,
                                     NGLL3*(*NSPEC_OUTER_CORE_ROTATION)*sizeof(realw)),9002);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_B_array_rotation, B_array_rotation,
                                     NGLL3*(*NSPEC_OUTER_CORE_ROTATION)*sizeof(realw),cudaMemcpyHostToDevice),9003);

  // backward/reconstructed fields
  if( mp->simulation_type == 3 ){
    mp->d_b_two_omega_earth = *b_two_omega_earth;
    mp->d_b_deltat = *b_deltat;

    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_b_A_array_rotation,
                                       NGLL3*(*NSPEC_OUTER_CORE_ROTATION)*sizeof(realw)),9000);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_A_array_rotation, b_A_array_rotation,
                                       NGLL3*(*NSPEC_OUTER_CORE_ROTATION)*sizeof(realw),cudaMemcpyHostToDevice),9001);

    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_b_B_array_rotation,
                                       NGLL3*(*NSPEC_OUTER_CORE_ROTATION)*sizeof(realw)),9002);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_B_array_rotation, b_B_array_rotation,
                                       NGLL3*(*NSPEC_OUTER_CORE_ROTATION)*sizeof(realw),cudaMemcpyHostToDevice),9003);
  }
}


/* ----------------------------------------------------------------------------------------------- */

// GRAVITY simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_gravity_device,
              PREPARE_FIELDS_gravity_DEVICE)(long* Mesh_pointer_f,
                                             realw* d_ln_density_dr_table,
                                             realw* minus_rho_g_over_kappa_fluid,
                                             realw* minus_gravity_table,
                                             realw* minus_deriv_gravity_table,
                                             realw* density_table,
                                             realw* h_wgll_cube,
                                             int* NRAD_GRAVITY
                                             ) {

  TRACE("prepare_fields_gravity_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  if( ! mp->gravity ){
    // no gravity case

    // d ln(rho)/dr needed for the no gravity fluid potential
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_d_ln_density_dr_table,
                                       (*NRAD_GRAVITY)*sizeof(realw)),8000);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_d_ln_density_dr_table, d_ln_density_dr_table,
                                       (*NRAD_GRAVITY)*sizeof(realw),cudaMemcpyHostToDevice),8001);

  }else{
    // gravity case

    // sets up gll weights cubed
    setConst_wgll_cube(h_wgll_cube,mp);

    // prepares gravity arrays
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_minus_rho_g_over_kappa_fluid,
                                       (*NRAD_GRAVITY)*sizeof(realw)),8000);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_minus_rho_g_over_kappa_fluid, minus_rho_g_over_kappa_fluid,
                                       (*NRAD_GRAVITY)*sizeof(realw),cudaMemcpyHostToDevice),8001);

    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_minus_gravity_table,
                                       (*NRAD_GRAVITY)*sizeof(realw)),8000);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_minus_gravity_table, minus_gravity_table,
                                       (*NRAD_GRAVITY)*sizeof(realw),cudaMemcpyHostToDevice),8001);

    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_minus_deriv_gravity_table,
                                       (*NRAD_GRAVITY)*sizeof(realw)),8000);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_minus_deriv_gravity_table, minus_deriv_gravity_table,
                                       (*NRAD_GRAVITY)*sizeof(realw),cudaMemcpyHostToDevice),8001);

    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_density_table,
                                       (*NRAD_GRAVITY)*sizeof(realw)),8000);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_density_table, density_table,
                                       (*NRAD_GRAVITY)*sizeof(realw),cudaMemcpyHostToDevice),8001);
  }
}



/* ----------------------------------------------------------------------------------------------- */

// ATTENUATION simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_attenuat_device,
              PREPARE_FIELDS_ATTENUAT_DEVICE)(long* Mesh_pointer_f,
                                                 realw* R_xx_crust_mantle,
                                                 realw* R_yy_crust_mantle,
                                                 realw* R_xy_crust_mantle,
                                                 realw* R_xz_crust_mantle,
                                                 realw* R_yz_crust_mantle,
                                                 realw* factor_common_crust_mantle,
                                                 realw* one_minus_sum_beta_crust_mantle,
                                                 realw* R_xx_inner_core,
                                                 realw* R_yy_inner_core,
                                                 realw* R_xy_inner_core,
                                                 realw* R_xz_inner_core,
                                                 realw* R_yz_inner_core,
                                                 realw* factor_common_inner_core,
                                                 realw* one_minus_sum_beta_inner_core,
                                                 realw* alphaval,realw* betaval,realw* gammaval,
                                                 realw* b_alphaval,realw* b_betaval,realw* b_gammaval
                                                 ) {

  TRACE("prepare_fields_attenuat_device");
  int R_size1,R_size2,R_size3;

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // checks flag
  if( ! mp->attenuation ){ exit_on_cuda_error("prepare_fields_attenuat_device attenuation not properly initialized"); }

  // crust_mantle
  R_size1 = N_SLS*NGLL3*mp->NSPEC_CRUST_MANTLE;
  R_size2 = NGLL3*mp->NSPEC_CRUST_MANTLE;
  R_size3 = N_SLS*NGLL3*mp->NSPEC_CRUST_MANTLE;

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_one_minus_sum_beta_crust_mantle,
                                     R_size2*sizeof(realw)),4430);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_one_minus_sum_beta_crust_mantle,one_minus_sum_beta_crust_mantle,
                                     R_size2*sizeof(realw),cudaMemcpyHostToDevice),4431);

  if( ! mp->use_attenuation_mimic ){
    // common factor
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_factor_common_crust_mantle,
                                         R_size3*sizeof(realw)),4432);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_factor_common_crust_mantle,factor_common_crust_mantle,
                                         R_size3*sizeof(realw),cudaMemcpyHostToDevice),4433);

    // memory variables
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_R_xx_crust_mantle,
                                       R_size1*sizeof(realw)),4401);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_R_yy_crust_mantle,
                                       R_size1*sizeof(realw)),4401);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_R_xy_crust_mantle,
                                       R_size1*sizeof(realw)),4401);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_R_xz_crust_mantle,
                                       R_size1*sizeof(realw)),4401);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_R_yz_crust_mantle,
                                       R_size1*sizeof(realw)),4401);

    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_xx_crust_mantle,R_xx_crust_mantle,
                                         R_size1*sizeof(realw),cudaMemcpyHostToDevice),4800);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_yy_crust_mantle,R_yy_crust_mantle,
                                         R_size1*sizeof(realw),cudaMemcpyHostToDevice),4800);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_xy_crust_mantle,R_xy_crust_mantle,
                                         R_size1*sizeof(realw),cudaMemcpyHostToDevice),4800);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_xz_crust_mantle,R_xz_crust_mantle,
                                         R_size1*sizeof(realw),cudaMemcpyHostToDevice),4800);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_yz_crust_mantle,R_yz_crust_mantle,
                                         R_size1*sizeof(realw),cudaMemcpyHostToDevice),4800);
  }

  // inner_core
  R_size1 = 5*N_SLS*NGLL3*mp->NSPEC_INNER_CORE;
  R_size2 = NGLL3*mp->NSPEC_INNER_CORE;
  R_size3 = N_SLS*NGLL3*mp->NSPEC_INNER_CORE;

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_one_minus_sum_beta_inner_core,
                                     R_size2*sizeof(realw)),4430);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_one_minus_sum_beta_inner_core,one_minus_sum_beta_inner_core,
                                     R_size2*sizeof(realw),cudaMemcpyHostToDevice),4431);

  if( ! mp->use_attenuation_mimic ){
    // common factor
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_factor_common_inner_core,
                                       R_size3*sizeof(realw)),4432);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_factor_common_inner_core,factor_common_inner_core,
                                       R_size3*sizeof(realw),cudaMemcpyHostToDevice),4433);

    // memory variables
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_R_xx_inner_core,
                                     R_size1*sizeof(realw)),4401);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_R_yy_inner_core,
                                       R_size1*sizeof(realw)),4401);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_R_xy_inner_core,
                                       R_size1*sizeof(realw)),4401);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_R_xz_inner_core,
                                       R_size1*sizeof(realw)),4401);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_R_yz_inner_core,
                                       R_size1*sizeof(realw)),4401);

    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_xx_inner_core,R_xx_inner_core,
                                     R_size1*sizeof(realw),cudaMemcpyHostToDevice),4402);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_yy_inner_core,R_yy_inner_core,
                                       R_size1*sizeof(realw),cudaMemcpyHostToDevice),4402);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_xy_inner_core,R_xy_inner_core,
                                       R_size1*sizeof(realw),cudaMemcpyHostToDevice),4402);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_xz_inner_core,R_xz_inner_core,
                                       R_size1*sizeof(realw),cudaMemcpyHostToDevice),4402);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_yz_inner_core,R_yz_inner_core,
                                       R_size1*sizeof(realw),cudaMemcpyHostToDevice),4402);
  }

  // alpha,beta,gamma factors
  print_CUDA_error_if_any(cudaMalloc((void**) &(mp->d_alphaval),
                                     N_SLS*sizeof(realw)),4434);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_alphaval ,alphaval,
                                     N_SLS*sizeof(realw),cudaMemcpyHostToDevice),4435);

  print_CUDA_error_if_any(cudaMalloc((void**) &(mp->d_betaval),
                                     N_SLS*sizeof(realw)),4436);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_betaval ,betaval,
                                     N_SLS*sizeof(realw),cudaMemcpyHostToDevice),4437);

  print_CUDA_error_if_any(cudaMalloc((void**) &(mp->d_gammaval),
                                     N_SLS*sizeof(realw)),4438);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_gammaval ,gammaval,
                                     N_SLS*sizeof(realw),cudaMemcpyHostToDevice),4439);

  if( mp->simulation_type == 3 ){
    // alpha,beta,gamma factors for backward fields
    print_CUDA_error_if_any(cudaMalloc((void**) &(mp->d_b_alphaval),
                                       N_SLS*sizeof(realw)),5434);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_alphaval ,b_alphaval,
                                       N_SLS*sizeof(realw),cudaMemcpyHostToDevice),5435);

    print_CUDA_error_if_any(cudaMalloc((void**) &(mp->d_b_betaval),
                                       N_SLS*sizeof(realw)),5436);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_betaval ,b_betaval,
                                       N_SLS*sizeof(realw),cudaMemcpyHostToDevice),5437);

    print_CUDA_error_if_any(cudaMalloc((void**) &(mp->d_b_gammaval),
                                       N_SLS*sizeof(realw)),5438);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_gammaval ,b_gammaval,
                                       N_SLS*sizeof(realw),cudaMemcpyHostToDevice),5439);
  }
}

/* ----------------------------------------------------------------------------------------------- */

// STRAIN simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_strain_device,
              PREPARE_FIELDS_STRAIN_DEVICE)(long* Mesh_pointer_f,
                                            realw* epsilondev_xx_crust_mantle,
                                            realw* epsilondev_yy_crust_mantle,
                                            realw* epsilondev_xy_crust_mantle,
                                            realw* epsilondev_xz_crust_mantle,
                                            realw* epsilondev_yz_crust_mantle,
                                            realw* b_epsilondev_xx_crust_mantle,
                                            realw* b_epsilondev_yy_crust_mantle,
                                            realw* b_epsilondev_xy_crust_mantle,
                                            realw* b_epsilondev_xz_crust_mantle,
                                            realw* b_epsilondev_yz_crust_mantle,
                                            realw* eps_trace_over_3_crust_mantle,
                                            realw* b_eps_trace_over_3_crust_mantle,
                                            realw* epsilondev_xx_inner_core,
                                            realw* epsilondev_yy_inner_core,
                                            realw* epsilondev_xy_inner_core,
                                            realw* epsilondev_xz_inner_core,
                                            realw* epsilondev_yz_inner_core,
                                            realw* b_epsilondev_xx_inner_core,
                                            realw* b_epsilondev_yy_inner_core,
                                            realw* b_epsilondev_xy_inner_core,
                                            realw* b_epsilondev_xz_inner_core,
                                            realw* b_epsilondev_yz_inner_core,
                                            realw* eps_trace_over_3_inner_core,
                                            realw* b_eps_trace_over_3_inner_core
                                            ) {

  TRACE("prepare_fields_strain_device");
  int R_size,size_strain_only;

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // checks flag
  if( ! mp->compute_and_store_strain ){ exit_on_cuda_error("prepare_fields_strain_device strain not properly initialized"); }

  // crust_mantle
  R_size = NGLL3*mp->NSPEC_CRUST_MANTLE;
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_epsilondev_xx_crust_mantle,
                                     R_size*sizeof(realw)),4432);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_epsilondev_yy_crust_mantle,
                                     R_size*sizeof(realw)),4432);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_epsilondev_xy_crust_mantle,
                                     R_size*sizeof(realw)),4432);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_epsilondev_xz_crust_mantle,
                                     R_size*sizeof(realw)),4432);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_epsilondev_yz_crust_mantle,
                                     R_size*sizeof(realw)),4432);

  print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_xx_crust_mantle,epsilondev_xx_crust_mantle,
                                     R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_yy_crust_mantle,epsilondev_yy_crust_mantle,
                                     R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_xy_crust_mantle,epsilondev_xy_crust_mantle,
                                     R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_xz_crust_mantle,epsilondev_xz_crust_mantle,
                                     R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_yz_crust_mantle,epsilondev_yz_crust_mantle,
                                     R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);

  // strain
  size_strain_only = NGLL3*(mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_eps_trace_over_3_crust_mantle,
                                      size_strain_only*sizeof(realw)),4401);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_eps_trace_over_3_crust_mantle,eps_trace_over_3_crust_mantle,
                                      size_strain_only*sizeof(realw),cudaMemcpyHostToDevice),4402);

  // backward/reconstructed fields
  if( mp->simulation_type == 3 ){
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_b_epsilondev_xx_crust_mantle,
                                       R_size*sizeof(realw)),4432);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_b_epsilondev_yy_crust_mantle,
                                       R_size*sizeof(realw)),4432);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_b_epsilondev_xy_crust_mantle,
                                       R_size*sizeof(realw)),4432);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_b_epsilondev_xz_crust_mantle,
                                       R_size*sizeof(realw)),4432);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_b_epsilondev_yz_crust_mantle,
                                       R_size*sizeof(realw)),4432);

    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xx_crust_mantle,b_epsilondev_xx_crust_mantle,
                                       R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_yy_crust_mantle,b_epsilondev_yy_crust_mantle,
                                       R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xy_crust_mantle,b_epsilondev_xy_crust_mantle,
                                       R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xz_crust_mantle,b_epsilondev_xz_crust_mantle,
                                       R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_yz_crust_mantle,b_epsilondev_yz_crust_mantle,
                                       R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);


    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_b_eps_trace_over_3_crust_mantle,
                                         R_size*sizeof(realw)),4401);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_eps_trace_over_3_crust_mantle,b_eps_trace_over_3_crust_mantle,
                                         R_size*sizeof(realw),cudaMemcpyHostToDevice),4402);
  }

  // inner_core
  R_size = NGLL3*mp->NSPEC_INNER_CORE;
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_epsilondev_xx_inner_core,
                                     R_size*sizeof(realw)),4432);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_epsilondev_yy_inner_core,
                                     R_size*sizeof(realw)),4432);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_epsilondev_xy_inner_core,
                                     R_size*sizeof(realw)),4432);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_epsilondev_xz_inner_core,
                                     R_size*sizeof(realw)),4432);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_epsilondev_yz_inner_core,
                                     R_size*sizeof(realw)),4432);

  print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_xx_inner_core,epsilondev_xx_inner_core,
                                     R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_yy_inner_core,epsilondev_yy_inner_core,
                                     R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_xy_inner_core,epsilondev_xy_inner_core,
                                     R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_xz_inner_core,epsilondev_xz_inner_core,
                                     R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_yz_inner_core,epsilondev_yz_inner_core,
                                     R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);


  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_eps_trace_over_3_inner_core,
                                     R_size*sizeof(realw)),4401);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_eps_trace_over_3_inner_core,eps_trace_over_3_inner_core,
                                     R_size*sizeof(realw),cudaMemcpyHostToDevice),4402);
  // backward/reconstructed fields
  if( mp->simulation_type == 3 ){
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_b_epsilondev_xx_inner_core,
                                       R_size*sizeof(realw)),4432);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_b_epsilondev_yy_inner_core,
                                       R_size*sizeof(realw)),4432);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_b_epsilondev_xy_inner_core,
                                       R_size*sizeof(realw)),4432);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_b_epsilondev_xz_inner_core,
                                       R_size*sizeof(realw)),4432);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_b_epsilondev_yz_inner_core,
                                       R_size*sizeof(realw)),4432);

    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xx_inner_core,b_epsilondev_xx_inner_core,
                                       R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_yy_inner_core,b_epsilondev_yy_inner_core,
                                       R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xy_inner_core,b_epsilondev_xy_inner_core,
                                       R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_xz_inner_core,b_epsilondev_xz_inner_core,
                                       R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_epsilondev_yz_inner_core,b_epsilondev_yz_inner_core,
                                       R_size*sizeof(realw),cudaMemcpyHostToDevice),4433);


    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_b_eps_trace_over_3_inner_core,
                                       R_size*sizeof(realw)),4401);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_b_eps_trace_over_3_inner_core,b_eps_trace_over_3_inner_core,
                                       R_size*sizeof(realw),cudaMemcpyHostToDevice),4402);
  }
}

/* ----------------------------------------------------------------------------------------------- */

// STRAIN simulations

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_fields_absorb_device,
              PREPARE_FIELDS_ABSORB_DEVICE)(long* Mesh_pointer_f,
                                            int* nspec2D_xmin_crust_mantle,int* nspec2D_xmax_crust_mantle,
                                            int* nspec2D_ymin_crust_mantle,int* nspec2D_ymax_crust_mantle,
                                            int* NSPEC2DMAX_XMIN_XMAX_CM,int* NSPEC2DMAX_YMIN_YMAX_CM,
                                            int* nimin_crust_mantle,int* nimax_crust_mantle,
                                            int* njmin_crust_mantle,int* njmax_crust_mantle,
                                            int* nkmin_xi_crust_mantle,int* nkmin_eta_crust_mantle,
                                            int* ibelm_xmin_crust_mantle,int* ibelm_xmax_crust_mantle,
                                            int* ibelm_ymin_crust_mantle,int* ibelm_ymax_crust_mantle,
                                            realw* normal_xmin_crust_mantle,realw* normal_xmax_crust_mantle,
                                            realw* normal_ymin_crust_mantle,realw* normal_ymax_crust_mantle,
                                            realw* jacobian2D_xmin_crust_mantle, realw* jacobian2D_xmax_crust_mantle,
                                            realw* jacobian2D_ymin_crust_mantle, realw* jacobian2D_ymax_crust_mantle,
                                            realw* rho_vp_crust_mantle,
                                            realw* rho_vs_crust_mantle,
                                            int* nspec2D_xmin_outer_core,int* nspec2D_xmax_outer_core,
                                            int* nspec2D_ymin_outer_core,int* nspec2D_ymax_outer_core,
                                            int* nspec2D_zmin_outer_core,
                                            int* NSPEC2DMAX_XMIN_XMAX_OC,int* NSPEC2DMAX_YMIN_YMAX_OC,
                                            int* nimin_outer_core,int* nimax_outer_core,
                                            int* njmin_outer_core,int* njmax_outer_core,
                                            int* nkmin_xi_outer_core,int* nkmin_eta_outer_core,
                                            int* ibelm_xmin_outer_core,int* ibelm_xmax_outer_core,
                                            int* ibelm_ymin_outer_core,int* ibelm_ymax_outer_core,
                                            int* ibelm_bottom_outer_core,
                                            realw* jacobian2D_xmin_outer_core, realw* jacobian2D_xmax_outer_core,
                                            realw* jacobian2D_ymin_outer_core, realw* jacobian2D_ymax_outer_core,
                                            realw* jacobian2D_bottom_outer_core,
                                            realw* vp_outer_core
                                            ) {

  TRACE("prepare_fields_absorb_device");
  int size;

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // checks flag
  if( ! mp->absorbing_conditions ){ exit_on_cuda_error("prepare_fields_absorb_device absorbing_conditions not properly initialized"); }

  // crust_mantle
  mp->nspec2D_xmin_crust_mantle = *nspec2D_xmin_crust_mantle;
  mp->nspec2D_xmax_crust_mantle = *nspec2D_xmax_crust_mantle;
  mp->nspec2D_ymin_crust_mantle = *nspec2D_ymin_crust_mantle;
  mp->nspec2D_ymax_crust_mantle = *nspec2D_ymax_crust_mantle;

  // vp & vs
  size = NGLL3*(mp->NSPEC_CRUST_MANTLE);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_rho_vp_crust_mantle,
                                     size*sizeof(realw)),2201);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_rho_vp_crust_mantle,rho_vp_crust_mantle,
                                     size*sizeof(realw),cudaMemcpyHostToDevice),2202);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_rho_vs_crust_mantle,
                                     size*sizeof(realw)),2201);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_rho_vs_crust_mantle,rho_vs_crust_mantle,
                                     size*sizeof(realw),cudaMemcpyHostToDevice),2202);

  // ijk index arrays
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_nkmin_xi_crust_mantle,
                                     2*(*NSPEC2DMAX_XMIN_XMAX_CM)*sizeof(int)),1201);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_nkmin_xi_crust_mantle,nkmin_xi_crust_mantle,
                                     2*(*NSPEC2DMAX_XMIN_XMAX_CM)*sizeof(int),cudaMemcpyHostToDevice),1202);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_nkmin_eta_crust_mantle,
                                     2*(*NSPEC2DMAX_YMIN_YMAX_CM)*sizeof(int)),1201);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_nkmin_eta_crust_mantle,nkmin_eta_crust_mantle,
                                     2*(*NSPEC2DMAX_YMIN_YMAX_CM)*sizeof(int),cudaMemcpyHostToDevice),1202);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_njmin_crust_mantle,
                                     2*(*NSPEC2DMAX_XMIN_XMAX_CM)*sizeof(int)),1201);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_njmin_crust_mantle,njmin_crust_mantle,
                                     2*(*NSPEC2DMAX_XMIN_XMAX_CM)*sizeof(int),cudaMemcpyHostToDevice),1202);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_njmax_crust_mantle,
                                     2*(*NSPEC2DMAX_XMIN_XMAX_CM)*sizeof(int)),1201);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_njmax_crust_mantle,njmax_crust_mantle,
                                     2*(*NSPEC2DMAX_XMIN_XMAX_CM)*sizeof(int),cudaMemcpyHostToDevice),1202);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_nimin_crust_mantle,
                                     2*(*NSPEC2DMAX_YMIN_YMAX_CM)*sizeof(int)),1201);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_nimin_crust_mantle,nimin_crust_mantle,
                                     2*(*NSPEC2DMAX_YMIN_YMAX_CM)*sizeof(int),cudaMemcpyHostToDevice),1202);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_nimax_crust_mantle,
                                     2*(*NSPEC2DMAX_YMIN_YMAX_CM)*sizeof(int)),1201);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_nimax_crust_mantle,nimax_crust_mantle,
                                     2*(*NSPEC2DMAX_YMIN_YMAX_CM)*sizeof(int),cudaMemcpyHostToDevice),1202);


  // xmin
  if( mp->nspec2D_xmin_crust_mantle > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibelm_xmin_crust_mantle,
                                       (mp->nspec2D_xmin_crust_mantle)*sizeof(int)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_ibelm_xmin_crust_mantle,ibelm_xmin_crust_mantle,
                                       (mp->nspec2D_xmin_crust_mantle)*sizeof(int),cudaMemcpyHostToDevice),1202);

    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_normal_xmin_crust_mantle,
                            NDIM*NGLL2*(mp->nspec2D_xmin_crust_mantle)*sizeof(realw)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_normal_xmin_crust_mantle,normal_xmin_crust_mantle,
                            NDIM*NGLL2*(mp->nspec2D_xmin_crust_mantle)*sizeof(realw),cudaMemcpyHostToDevice),1202);

    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_jacobian2D_xmin_crust_mantle,
                            NGLL2*(mp->nspec2D_xmin_crust_mantle)*sizeof(realw)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_jacobian2D_xmin_crust_mantle,jacobian2D_xmin_crust_mantle,
                            NGLL2*(mp->nspec2D_xmin_crust_mantle)*sizeof(realw),cudaMemcpyHostToDevice),1202);

    // boundary buffer
    if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
      print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_absorb_xmin_crust_mantle,
                              NDIM*NGLL2*(mp->nspec2D_xmin_crust_mantle)*sizeof(realw)),1202);
    }
  }

  // xmax
  if( mp->nspec2D_xmax_crust_mantle > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibelm_xmax_crust_mantle,
                                       (mp->nspec2D_xmax_crust_mantle)*sizeof(int)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_ibelm_xmax_crust_mantle,ibelm_xmax_crust_mantle,
                                       (mp->nspec2D_xmax_crust_mantle)*sizeof(int),cudaMemcpyHostToDevice),1202);

    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_normal_xmax_crust_mantle,
                                       NDIM*NGLL2*(mp->nspec2D_xmax_crust_mantle)*sizeof(realw)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_normal_xmax_crust_mantle,normal_xmax_crust_mantle,
                                       NDIM*NGLL2*(mp->nspec2D_xmax_crust_mantle)*sizeof(realw),cudaMemcpyHostToDevice),1202);

    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_jacobian2D_xmax_crust_mantle,
                                       NGLL2*(mp->nspec2D_xmax_crust_mantle)*sizeof(realw)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_jacobian2D_xmax_crust_mantle,jacobian2D_xmax_crust_mantle,
                                       NGLL2*(mp->nspec2D_xmax_crust_mantle)*sizeof(realw),cudaMemcpyHostToDevice),1202);

    // boundary buffer
    if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
      print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_absorb_xmax_crust_mantle,
                                         NDIM*NGLL2*(mp->nspec2D_xmax_crust_mantle)*sizeof(realw)),1202);
    }
  }

  // ymin
  if( mp->nspec2D_ymin_crust_mantle > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibelm_ymin_crust_mantle,
                                       (mp->nspec2D_ymin_crust_mantle)*sizeof(int)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_ibelm_ymin_crust_mantle,ibelm_ymin_crust_mantle,
                                       (mp->nspec2D_ymin_crust_mantle)*sizeof(int),cudaMemcpyHostToDevice),1202);

    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_normal_ymin_crust_mantle,
                                       NDIM*NGLL2*(mp->nspec2D_ymin_crust_mantle)*sizeof(realw)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_normal_ymin_crust_mantle,normal_ymin_crust_mantle,
                                       NDIM*NGLL2*(mp->nspec2D_ymin_crust_mantle)*sizeof(realw),cudaMemcpyHostToDevice),1202);

    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_jacobian2D_ymin_crust_mantle,
                                       NGLL2*(mp->nspec2D_ymin_crust_mantle)*sizeof(realw)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_jacobian2D_ymin_crust_mantle,jacobian2D_ymin_crust_mantle,
                                       NGLL2*(mp->nspec2D_ymin_crust_mantle)*sizeof(realw),cudaMemcpyHostToDevice),1202);

    // boundary buffer
    if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
      print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_absorb_ymin_crust_mantle,
                                         NDIM*NGLL2*(mp->nspec2D_ymin_crust_mantle)*sizeof(realw)),1202);
    }
  }

  // ymax
  if( mp->nspec2D_ymax_crust_mantle > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibelm_ymax_crust_mantle,
                                       (mp->nspec2D_ymax_crust_mantle)*sizeof(int)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_ibelm_ymax_crust_mantle,ibelm_ymax_crust_mantle,
                                       (mp->nspec2D_ymax_crust_mantle)*sizeof(int),cudaMemcpyHostToDevice),1202);

    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_normal_ymax_crust_mantle,
                                       NDIM*NGLL2*(mp->nspec2D_ymax_crust_mantle)*sizeof(realw)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_normal_ymax_crust_mantle,normal_ymax_crust_mantle,
                                       NDIM*NGLL2*(mp->nspec2D_ymax_crust_mantle)*sizeof(realw),cudaMemcpyHostToDevice),1202);

    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_jacobian2D_ymax_crust_mantle,
                                       NGLL2*(mp->nspec2D_ymax_crust_mantle)*sizeof(realw)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_jacobian2D_ymax_crust_mantle,jacobian2D_ymax_crust_mantle,
                                       NGLL2*(mp->nspec2D_ymax_crust_mantle)*sizeof(realw),cudaMemcpyHostToDevice),1202);

    // boundary buffer
    if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
      print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_absorb_ymax_crust_mantle,
                                         NDIM*NGLL2*(mp->nspec2D_ymax_crust_mantle)*sizeof(realw)),1202);
    }
  }


  // outer_core
  mp->nspec2D_xmin_outer_core = *nspec2D_xmin_outer_core;
  mp->nspec2D_xmax_outer_core = *nspec2D_xmax_outer_core;
  mp->nspec2D_ymin_outer_core = *nspec2D_ymin_outer_core;
  mp->nspec2D_ymax_outer_core = *nspec2D_ymax_outer_core;
  mp->nspec2D_zmin_outer_core = *nspec2D_zmin_outer_core;

  // vp
  size = NGLL3*(mp->NSPEC_OUTER_CORE);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_vp_outer_core,
                                     size*sizeof(realw)),2201);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_vp_outer_core,vp_outer_core,
                                     size*sizeof(realw),cudaMemcpyHostToDevice),2202);

  // ijk index arrays
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_nkmin_xi_outer_core,
                                     2*(*NSPEC2DMAX_XMIN_XMAX_OC)*sizeof(int)),1201);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_nkmin_xi_outer_core,nkmin_xi_outer_core,
                                     2*(*NSPEC2DMAX_XMIN_XMAX_OC)*sizeof(int),cudaMemcpyHostToDevice),1202);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_nkmin_eta_outer_core,
                                     2*(*NSPEC2DMAX_YMIN_YMAX_OC)*sizeof(int)),1201);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_nkmin_eta_outer_core,nkmin_eta_outer_core,
                                     2*(*NSPEC2DMAX_YMIN_YMAX_OC)*sizeof(int),cudaMemcpyHostToDevice),1202);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_njmin_outer_core,
                                     2*(*NSPEC2DMAX_XMIN_XMAX_OC)*sizeof(int)),1201);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_njmin_outer_core,njmin_outer_core,
                                     2*(*NSPEC2DMAX_XMIN_XMAX_OC)*sizeof(int),cudaMemcpyHostToDevice),1202);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_njmax_outer_core,
                                     2*(*NSPEC2DMAX_XMIN_XMAX_OC)*sizeof(int)),1201);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_njmax_outer_core,njmax_outer_core,
                                     2*(*NSPEC2DMAX_XMIN_XMAX_OC)*sizeof(int),cudaMemcpyHostToDevice),1202);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_nimin_outer_core,
                                     2*(*NSPEC2DMAX_YMIN_YMAX_OC)*sizeof(int)),1201);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_nimin_outer_core,nimin_outer_core,
                                     2*(*NSPEC2DMAX_YMIN_YMAX_OC)*sizeof(int),cudaMemcpyHostToDevice),1202);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_nimax_outer_core,
                                     2*(*NSPEC2DMAX_YMIN_YMAX_OC)*sizeof(int)),1201);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_nimax_outer_core,nimax_outer_core,
                                     2*(*NSPEC2DMAX_YMIN_YMAX_OC)*sizeof(int),cudaMemcpyHostToDevice),1202);

  // xmin
  if( mp->nspec2D_xmin_outer_core > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibelm_xmin_outer_core,
                                       (mp->nspec2D_xmin_outer_core)*sizeof(int)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_ibelm_xmin_outer_core,ibelm_xmin_outer_core,
                                       (mp->nspec2D_xmin_outer_core)*sizeof(int),cudaMemcpyHostToDevice),1202);

    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_jacobian2D_xmin_outer_core,
                                       NGLL2*(mp->nspec2D_xmin_outer_core)*sizeof(realw)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_jacobian2D_xmin_outer_core,jacobian2D_xmin_outer_core,
                                       NGLL2*(mp->nspec2D_xmin_outer_core)*sizeof(realw),cudaMemcpyHostToDevice),1202);

    // boundary buffer
    if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
      print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_absorb_xmin_outer_core,
                                         NGLL2*(mp->nspec2D_xmin_outer_core)*sizeof(realw)),1202);
    }
  }

  // xmax
  if( mp->nspec2D_xmax_outer_core > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibelm_xmax_outer_core,
                                       (mp->nspec2D_xmax_outer_core)*sizeof(int)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_ibelm_xmax_outer_core,ibelm_xmax_outer_core,
                                       (mp->nspec2D_xmax_outer_core)*sizeof(int),cudaMemcpyHostToDevice),1202);

    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_jacobian2D_xmax_outer_core,
                                       NGLL2*(mp->nspec2D_xmax_outer_core)*sizeof(realw)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_jacobian2D_xmax_outer_core,jacobian2D_xmax_outer_core,
                                       NGLL2*(mp->nspec2D_xmax_outer_core)*sizeof(realw),cudaMemcpyHostToDevice),1202);

    // boundary buffer
    if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
      print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_absorb_xmax_outer_core,
                                         NGLL2*(mp->nspec2D_xmax_outer_core)*sizeof(realw)),1202);
    }
  }

  // ymin
  if( mp->nspec2D_ymin_outer_core > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibelm_ymin_outer_core,
                                       (mp->nspec2D_ymin_outer_core)*sizeof(int)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_ibelm_ymin_outer_core,ibelm_ymin_outer_core,
                                       (mp->nspec2D_ymin_outer_core)*sizeof(int),cudaMemcpyHostToDevice),1202);

    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_jacobian2D_ymin_outer_core,
                                       NGLL2*(mp->nspec2D_ymin_outer_core)*sizeof(realw)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_jacobian2D_ymin_outer_core,jacobian2D_ymin_outer_core,
                                       NGLL2*(mp->nspec2D_ymin_outer_core)*sizeof(realw),cudaMemcpyHostToDevice),1202);

    // boundary buffer
    if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
      print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_absorb_ymin_outer_core,
                                         NGLL2*(mp->nspec2D_ymin_outer_core)*sizeof(realw)),1202);
    }
  }

  // ymax
  if( mp->nspec2D_ymax_outer_core > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibelm_ymax_outer_core,
                                       (mp->nspec2D_ymax_outer_core)*sizeof(int)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_ibelm_ymax_outer_core,ibelm_ymax_outer_core,
                                       (mp->nspec2D_ymax_outer_core)*sizeof(int),cudaMemcpyHostToDevice),1202);

    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_jacobian2D_ymax_outer_core,
                                       NGLL2*(mp->nspec2D_ymax_outer_core)*sizeof(realw)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_jacobian2D_ymax_outer_core,jacobian2D_ymax_outer_core,
                                       NGLL2*(mp->nspec2D_ymax_outer_core)*sizeof(realw),cudaMemcpyHostToDevice),1202);

    // boundary buffer
    if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
      print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_absorb_ymax_outer_core,
                                         NGLL2*(mp->nspec2D_ymax_outer_core)*sizeof(realw)),1202);
    }
  }

  // zmin
  if( mp->nspec2D_zmin_outer_core > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibelm_zmin_outer_core,
                                       (mp->nspec2D_zmin_outer_core)*sizeof(int)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_ibelm_zmin_outer_core,ibelm_bottom_outer_core,
                                       (mp->nspec2D_zmin_outer_core)*sizeof(int),cudaMemcpyHostToDevice),1202);

    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_jacobian2D_zmin_outer_core,
                                       NGLL2*(mp->nspec2D_zmin_outer_core)*sizeof(realw)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_jacobian2D_zmin_outer_core,jacobian2D_bottom_outer_core,
                                       NGLL2*(mp->nspec2D_zmin_outer_core)*sizeof(realw),cudaMemcpyHostToDevice),1202);

    // boundary buffer
    if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
      print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_absorb_zmin_outer_core,
                                         NGLL2*(mp->nspec2D_zmin_outer_core)*sizeof(realw)),1202);
    }
  }

}

/* ----------------------------------------------------------------------------------------------- */

// MPI interfaces

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_mpi_buffers_device,
              PREPARE_MPI_BUFFERS_DEVICE)(long* Mesh_pointer_f,
                                          int* num_interfaces_crust_mantle,
                                          int* max_nibool_interfaces_crust_mantle,
                                          int* nibool_interfaces_crust_mantle,
                                          int* ibool_interfaces_crust_mantle,
                                          int* num_interfaces_inner_core,
                                          int* max_nibool_interfaces_inner_core,
                                          int* nibool_interfaces_inner_core,
                                          int* ibool_interfaces_inner_core,
                                          int* num_interfaces_outer_core,
                                          int* max_nibool_interfaces_outer_core,
                                          int* nibool_interfaces_outer_core,
                                          int* ibool_interfaces_outer_core
                                          ){

  TRACE("prepare_mpi_buffers_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // prepares interprocess-edge exchange information

  // crust/mantle mesh
  mp->num_interfaces_crust_mantle = *num_interfaces_crust_mantle;
  mp->max_nibool_interfaces_crust_mantle = *max_nibool_interfaces_crust_mantle;
  if( mp->num_interfaces_crust_mantle > 0 ){
    // number of ibool entries array
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_nibool_interfaces_crust_mantle,
                                       (mp->num_interfaces_crust_mantle)*sizeof(int)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_nibool_interfaces_crust_mantle,nibool_interfaces_crust_mantle,
                                       (mp->num_interfaces_crust_mantle)*sizeof(int),cudaMemcpyHostToDevice),1202);
    // ibool entries (iglob indices) values on interface
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibool_interfaces_crust_mantle,
                                       (mp->num_interfaces_crust_mantle)*(mp->max_nibool_interfaces_crust_mantle)*sizeof(int)),1203);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_ibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle,
                                       (mp->num_interfaces_crust_mantle)*(mp->max_nibool_interfaces_crust_mantle)*sizeof(int),
                                       cudaMemcpyHostToDevice),1204);
    // allocates mpi buffer for exchange with cpu
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_send_accel_buffer_crust_mantle),
                                       3*(mp->max_nibool_interfaces_crust_mantle)*(mp->num_interfaces_crust_mantle)*sizeof(realw)),4004);
  }

  // inner core mesh
  mp->num_interfaces_inner_core = *num_interfaces_inner_core;
  mp->max_nibool_interfaces_inner_core = *max_nibool_interfaces_inner_core;
  if( mp->num_interfaces_inner_core > 0 ){
    // number of ibool entries array
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_nibool_interfaces_inner_core,
                                       (mp->num_interfaces_inner_core)*sizeof(int)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_nibool_interfaces_inner_core,nibool_interfaces_inner_core,
                                       (mp->num_interfaces_inner_core)*sizeof(int),cudaMemcpyHostToDevice),1202);
    // ibool entries (iglob indices) values on interface
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibool_interfaces_inner_core,
                                       (mp->num_interfaces_inner_core)*(mp->max_nibool_interfaces_inner_core)*sizeof(int)),1203);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_ibool_interfaces_inner_core,ibool_interfaces_inner_core,
                                       (mp->num_interfaces_inner_core)*(mp->max_nibool_interfaces_inner_core)*sizeof(int),
                                       cudaMemcpyHostToDevice),1204);
    // allocates mpi buffer for exchange with cpu
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_send_accel_buffer_inner_core),
                                       3*(mp->max_nibool_interfaces_inner_core)*(mp->num_interfaces_inner_core)*sizeof(realw)),4004);
  }

  // outer core mesh
  // note: uses only scalar wavefield arrays
  mp->num_interfaces_outer_core = *num_interfaces_outer_core;
  mp->max_nibool_interfaces_outer_core = *max_nibool_interfaces_outer_core;
  if( mp->num_interfaces_outer_core > 0 ){
    // number of ibool entries array
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_nibool_interfaces_outer_core,
                                       (mp->num_interfaces_outer_core)*sizeof(int)),1201);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_nibool_interfaces_outer_core,nibool_interfaces_outer_core,
                                       (mp->num_interfaces_outer_core)*sizeof(int),cudaMemcpyHostToDevice),1202);
    // ibool entries (iglob indices) values on interface
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibool_interfaces_outer_core,
                                       (mp->num_interfaces_outer_core)*(mp->max_nibool_interfaces_outer_core)*sizeof(int)),1203);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_ibool_interfaces_outer_core,ibool_interfaces_outer_core,
                                       (mp->num_interfaces_outer_core)*(mp->max_nibool_interfaces_outer_core)*sizeof(int),
                                       cudaMemcpyHostToDevice),1204);
    // allocates mpi buffer for exchange with cpu
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_send_accel_buffer_outer_core),
                                       (mp->max_nibool_interfaces_outer_core)*(mp->num_interfaces_outer_core)*sizeof(realw)),4004);
  }


}

/* ----------------------------------------------------------------------------------------------- */

// for NOISE simulations

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(prepare_fields_noise_device,
              PREPARE_FIELDS_NOISE_DEVICE)(long* Mesh_pointer_f,
                                           int* nspec_top,
                                           int* ibelm_top_crust_mantle,
                                           int* NSTEP,
                                           realw* noise_sourcearray,
                                           realw* normal_x_noise,
                                           realw* normal_y_noise,
                                           realw* normal_z_noise,
                                           realw* mask_noise,
                                           realw* jacobian2D_top_crust_mantle) {

  TRACE("prepare_fields_noise_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // free surface
  mp->nspec_top = *nspec_top;
  if( mp->nspec_top > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibelm_top_crust_mantle,
                                       mp->nspec_top*sizeof(int)),7001);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_ibelm_top_crust_mantle,ibelm_top_crust_mantle,
                                       mp->nspec_top*sizeof(int),cudaMemcpyHostToDevice),7002);

    // alloc storage for the surface buffer to be copied
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_noise_surface_movie,
                                       NDIM*NGLL2*(mp->nspec_top)*sizeof(realw)),7005);
  }else{
    // for global mesh: each crust/mantle slice should have at top a free surface
    exit_on_cuda_error("prepare_fields_noise_device nspec_top not properly initialized");
  }


  // prepares noise source array
  if( mp->noise_tomography == 1 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_noise_sourcearray,
                                       NDIM*NGLL3*(*NSTEP)*sizeof(realw)),7101);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_noise_sourcearray,noise_sourcearray,
                                       NDIM*NGLL3*(*NSTEP)*sizeof(realw),cudaMemcpyHostToDevice),7102);
  }

  // prepares noise directions
  if( mp->noise_tomography > 1 ){
    int nface_size = NGLL2*(mp->nspec_top);
    // allocates memory on GPU
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_normal_x_noise,
                                       nface_size*sizeof(realw)),7301);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_normal_x_noise, normal_x_noise,
                                       nface_size*sizeof(realw),cudaMemcpyHostToDevice),7306);

    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_normal_y_noise,
                                       nface_size*sizeof(realw)),7302);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_normal_y_noise, normal_y_noise,
                                       nface_size*sizeof(realw),cudaMemcpyHostToDevice),7307);

    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_normal_z_noise,
                                       nface_size*sizeof(realw)),7303);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_normal_z_noise, normal_z_noise,
                                       nface_size*sizeof(realw),cudaMemcpyHostToDevice),7308);

    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_mask_noise,
                                       nface_size*sizeof(realw)),7304);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_mask_noise, mask_noise,
                                       nface_size*sizeof(realw),cudaMemcpyHostToDevice),7309);

    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_jacobian2D_top_crust_mantle,
                                       nface_size*sizeof(realw)),7305);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_jacobian2D_top_crust_mantle, jacobian2D_top_crust_mantle,
                                       nface_size*sizeof(realw),cudaMemcpyHostToDevice),7310);
  }

  // prepares noise strength kernel
  if( mp->noise_tomography == 3 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_Sigma_kl),
                                       NGLL3*(mp->NSPEC_CRUST_MANTLE)*sizeof(realw)),7401);
    // initializes kernel values to zero
    print_CUDA_error_if_any(cudaMemset(mp->d_Sigma_kl,0,
                                       NGLL3*mp->NSPEC_CRUST_MANTLE*sizeof(realw)),7403);

  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_fields_noise_device");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// Earth regions

// CRUST / MANTLE

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_crust_mantle_device,
             PREPARE_CRUST_MANTLE_DEVICE)(long* Mesh_pointer_f,
             realw* h_xix, realw* h_xiy, realw* h_xiz,
             realw* h_etax, realw* h_etay, realw* h_etaz,
             realw* h_gammax, realw* h_gammay, realw* h_gammaz,
             realw* h_rho,
             realw* h_kappav, realw* h_muv,
             realw* h_kappah, realw* h_muh,
             realw* h_eta_aniso,
             realw* h_rmassx,
             realw* h_rmassy,
             realw* h_rmassz,
             realw* h_normal_top_crust_mantle,
             int* h_ibelm_top_crust_mantle,
             int* h_ibelm_bottom_crust_mantle,
             int* h_ibool,
             realw* h_xstore, realw* h_ystore, realw* h_zstore,
             int* h_ispec_is_tiso,
             realw *c11store,realw *c12store,realw *c13store,
             realw *c14store,realw *c15store,realw *c16store,
             realw *c22store,realw *c23store,realw *c24store,
             realw *c25store,realw *c26store,realw *c33store,
             realw *c34store,realw *c35store,realw *c36store,
             realw *c44store,realw *c45store,realw *c46store,
             realw *c55store,realw *c56store,realw *c66store,
             int* num_phase_ispec,
             int* phase_ispec_inner,
             int* nspec_outer,
             int* nspec_inner,
             int* NSPEC2D_TOP_CM,
	     int* NSPEC2D_BOTTOM_CM,
	     int* NCHUNKS_VAL
             ) {

  TRACE("prepare_crust_mantle_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  /* Assuming NGLLX=5. Padded is then 128 (5^3+3) */
  int size_padded = NGLL3_PADDED * (mp->NSPEC_CRUST_MANTLE);
  int size_glob = mp->NGLOB_CRUST_MANTLE;

  // mesh
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xix_crust_mantle, size_padded*sizeof(realw)),1001);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xiy_crust_mantle, size_padded*sizeof(realw)),1002);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xiz_crust_mantle, size_padded*sizeof(realw)),1003);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etax_crust_mantle, size_padded*sizeof(realw)),1004);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etay_crust_mantle, size_padded*sizeof(realw)),1005);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etaz_crust_mantle, size_padded*sizeof(realw)),1006);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammax_crust_mantle, size_padded*sizeof(realw)),1007);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammay_crust_mantle, size_padded*sizeof(realw)),1008);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammaz_crust_mantle, size_padded*sizeof(realw)),1009);
  // muv needed in case for attenuation (only Q_mu shear attenuation)
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_muvstore_crust_mantle, size_padded*sizeof(realw)),1011);

  // transfer constant element data with padding
  for(int i=0;i < mp->NSPEC_CRUST_MANTLE;i++) {
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xix_crust_mantle + i*NGLL3_PADDED, &h_xix[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1501);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xiy_crust_mantle+i*NGLL3_PADDED,   &h_xiy[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1502);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xiz_crust_mantle+i*NGLL3_PADDED,   &h_xiz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1503);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etax_crust_mantle+i*NGLL3_PADDED,  &h_etax[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1504);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etay_crust_mantle+i*NGLL3_PADDED,  &h_etay[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1505);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etaz_crust_mantle+i*NGLL3_PADDED,  &h_etaz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1506);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammax_crust_mantle+i*NGLL3_PADDED,&h_gammax[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1507);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammay_crust_mantle+i*NGLL3_PADDED,&h_gammay[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1508);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammaz_crust_mantle+i*NGLL3_PADDED,&h_gammaz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1509);
    // muvstore
    print_CUDA_error_if_any(cudaMemcpy(mp->d_muvstore_crust_mantle+i*NGLL3_PADDED,   &h_muv[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1511);
  }

  // global indexing
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibool_crust_mantle, size_padded*sizeof(int)),1021);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_ibool_crust_mantle, h_ibool,
                                     NGLL3*(mp->NSPEC_CRUST_MANTLE)*sizeof(int),cudaMemcpyHostToDevice),1022);

  // transverse isotropic elements
  // only needed if not anisotropic 3D mantle
  if( ! mp->anisotropic_3D_mantle ){
    // no anisotropy

    // transverse isotropy flag
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ispec_is_tiso_crust_mantle, (mp->NSPEC_CRUST_MANTLE)*sizeof(int)),1025);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_ispec_is_tiso_crust_mantle, h_ispec_is_tiso,
                                       (mp->NSPEC_CRUST_MANTLE)*sizeof(int),cudaMemcpyHostToDevice),1025);

    // kappavstore, kappahstore
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_kappavstore_crust_mantle, size_padded*sizeof(realw)),1010);
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_kappahstore_crust_mantle, size_padded*sizeof(realw)),1010);
    // muhstore
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_muhstore_crust_mantle, size_padded*sizeof(realw)),1010);
    // eta_anisostore
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_eta_anisostore_crust_mantle, size_padded*sizeof(realw)),1010);

    // transfer with padding
    for(int i=0;i < mp->NSPEC_CRUST_MANTLE;i++) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_kappavstore_crust_mantle+i*NGLL3_PADDED,&h_kappav[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1510);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_kappahstore_crust_mantle+i*NGLL3_PADDED,&h_kappah[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1510);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_muhstore_crust_mantle+i*NGLL3_PADDED,&h_muh[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1510);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_eta_anisostore_crust_mantle+i*NGLL3_PADDED,&h_eta_aniso[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1510);
    }
  }else{
    // anisotropic 3D mantle

    // allocates memory on GPU
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c11store_crust_mantle),
                                       size_padded*sizeof(realw)),4700);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c12store_crust_mantle),
                                       size_padded*sizeof(realw)),4701);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c13store_crust_mantle),
                                       size_padded*sizeof(realw)),4702);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c14store_crust_mantle),
                                       size_padded*sizeof(realw)),4703);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c15store_crust_mantle),
                                       size_padded*sizeof(realw)),4704);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c16store_crust_mantle),
                                       size_padded*sizeof(realw)),4705);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c22store_crust_mantle),
                                       size_padded*sizeof(realw)),4706);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c23store_crust_mantle),
                                       size_padded*sizeof(realw)),4707);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c24store_crust_mantle),
                                       size_padded*sizeof(realw)),4708);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c25store_crust_mantle),
                                       size_padded*sizeof(realw)),4709);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c26store_crust_mantle),
                                       size_padded*sizeof(realw)),4710);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c33store_crust_mantle),
                                       size_padded*sizeof(realw)),4711);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c34store_crust_mantle),
                                       size_padded*sizeof(realw)),4712);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c35store_crust_mantle),
                                       size_padded*sizeof(realw)),4713);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c36store_crust_mantle),
                                       size_padded*sizeof(realw)),4714);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c44store_crust_mantle),
                                       size_padded*sizeof(realw)),4715);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c45store_crust_mantle),
                                       size_padded*sizeof(realw)),4716);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c46store_crust_mantle),
                                       size_padded*sizeof(realw)),4717);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c55store_crust_mantle),
                                       size_padded*sizeof(realw)),4718);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c56store_crust_mantle),
                                       size_padded*sizeof(realw)),4719);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c66store_crust_mantle),
                                       size_padded*sizeof(realw)),4720);

    // transfer constant element data with padding
    for(int i=0;i < mp->NSPEC_CRUST_MANTLE;i++) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c11store_crust_mantle + i*NGLL3_PADDED, &c11store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4800);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c12store_crust_mantle + i*NGLL3_PADDED, &c12store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4801);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c13store_crust_mantle + i*NGLL3_PADDED, &c13store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4802);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c14store_crust_mantle + i*NGLL3_PADDED, &c14store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4803);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c15store_crust_mantle + i*NGLL3_PADDED, &c15store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4804);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c16store_crust_mantle + i*NGLL3_PADDED, &c16store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4805);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c22store_crust_mantle + i*NGLL3_PADDED, &c22store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4806);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c23store_crust_mantle + i*NGLL3_PADDED, &c23store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4807);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c24store_crust_mantle + i*NGLL3_PADDED, &c24store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4808);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c25store_crust_mantle + i*NGLL3_PADDED, &c25store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4809);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c26store_crust_mantle + i*NGLL3_PADDED, &c26store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4810);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c33store_crust_mantle + i*NGLL3_PADDED, &c33store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4811);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c34store_crust_mantle + i*NGLL3_PADDED, &c34store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4812);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c35store_crust_mantle + i*NGLL3_PADDED, &c35store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4813);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c36store_crust_mantle + i*NGLL3_PADDED, &c36store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4814);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c44store_crust_mantle + i*NGLL3_PADDED, &c44store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4815);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c45store_crust_mantle + i*NGLL3_PADDED, &c45store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4816);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c46store_crust_mantle + i*NGLL3_PADDED, &c46store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4817);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c55store_crust_mantle + i*NGLL3_PADDED, &c55store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4818);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c56store_crust_mantle + i*NGLL3_PADDED, &c56store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4819);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c66store_crust_mantle + i*NGLL3_PADDED, &c66store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4820);
    }
  }

  // needed for boundary kernel calculations
  if( mp->simulation_type == 3 && mp->save_boundary_mesh ){
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_rhostore_crust_mantle, size_padded*sizeof(realw)),1010);
    for(int i=0;i < mp->NSPEC_CRUST_MANTLE;i++) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_rhostore_crust_mantle+i*NGLL3_PADDED, &h_rho[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),2106);
    }
  }

  // mesh locations
  // ystore & zstore needed for tiso elements
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_ystore_crust_mantle),sizeof(realw)*size_glob),2005);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_ystore_crust_mantle,h_ystore,
                                     sizeof(realw)*size_glob,cudaMemcpyHostToDevice),2100);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_zstore_crust_mantle),sizeof(realw)*size_glob),2005);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_zstore_crust_mantle,h_zstore,
                                     sizeof(realw)*size_glob,cudaMemcpyHostToDevice),2100);

  // xstore only needed when gravity is on
  if( mp->gravity ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_xstore_crust_mantle),sizeof(realw)*size_glob),2005);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xstore_crust_mantle,h_xstore,
                                     sizeof(realw)*size_glob,cudaMemcpyHostToDevice),2100);
  }

  // inner/outer elements
  mp->num_phase_ispec_crust_mantle = *num_phase_ispec;
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_phase_ispec_inner_crust_mantle),
                                     mp->num_phase_ispec_crust_mantle*2*sizeof(int)),2008);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_phase_ispec_inner_crust_mantle,phase_ispec_inner,
                                     mp->num_phase_ispec_crust_mantle*2*sizeof(int),cudaMemcpyHostToDevice),2101);

  mp->nspec_outer_crust_mantle = *nspec_outer;
  mp->nspec_inner_crust_mantle = *nspec_inner;

  // CMB/ocean coupling
  mp->nspec2D_top_crust_mantle = *NSPEC2D_TOP_CM;
  mp->nspec2D_bottom_crust_mantle = *NSPEC2D_BOTTOM_CM;
  int size_tcm = NGLL2*(mp->nspec2D_top_crust_mantle);

  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_normal_top_crust_mantle),sizeof(realw)*NDIM*size_tcm),40020);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_normal_top_crust_mantle,h_normal_top_crust_mantle,sizeof(realw)*NDIM*size_tcm,cudaMemcpyHostToDevice),40030);

  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_ibelm_top_crust_mantle),sizeof(int)*(mp->nspec2D_top_crust_mantle)),40021);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_ibelm_bottom_crust_mantle),sizeof(int)*(mp->nspec2D_bottom_crust_mantle)),40021);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_ibelm_top_crust_mantle,h_ibelm_top_crust_mantle,sizeof(int)*(mp->nspec2D_top_crust_mantle),cudaMemcpyHostToDevice),40031);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_ibelm_bottom_crust_mantle,h_ibelm_bottom_crust_mantle,sizeof(int)*(mp->nspec2D_bottom_crust_mantle),cudaMemcpyHostToDevice),40031);

  // wavefield
  int size = NDIM * mp->NGLOB_CRUST_MANTLE;

  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_displ_crust_mantle),sizeof(realw)*size),4001);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_veloc_crust_mantle),sizeof(realw)*size),4002);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_accel_crust_mantle),sizeof(realw)*size),4003);
  // backward/reconstructed wavefield
  if( mp->simulation_type == 3 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_displ_crust_mantle),sizeof(realw)*size),4001);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_veloc_crust_mantle),sizeof(realw)*size),4002);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_accel_crust_mantle),sizeof(realw)*size),4003);
  }

  // mass matrices
  if( *NCHUNKS_VAL != 6 && mp->absorbing_conditions){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rmassx_crust_mantle),sizeof(realw)*size_glob),2005);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_rmassx_crust_mantle,h_rmassx,
				       sizeof(realw)*size_glob,cudaMemcpyHostToDevice),2100);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rmassy_crust_mantle),sizeof(realw)*size_glob),2005);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_rmassy_crust_mantle,h_rmassy,
				       sizeof(realw)*size_glob,cudaMemcpyHostToDevice),2100);
  }

  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rmassz_crust_mantle),sizeof(realw)*size_glob),2005);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_rmassz_crust_mantle,h_rmassz,
                                     sizeof(realw)*size_glob,cudaMemcpyHostToDevice),2100);

  // kernels
  if( mp->simulation_type == 3 ){

    size = NGLL3*(mp->NSPEC_CRUST_MANTLE);

    // density kernel
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rho_kl_crust_mantle),
                                       size*sizeof(realw)),5204);
    // initializes kernel values to zero
    print_CUDA_error_if_any(cudaMemset(mp->d_rho_kl_crust_mantle,0,size*sizeof(realw)),5207);

    if( ! mp->anisotropic_kl){
      // isotropic kernels
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_alpha_kl_crust_mantle),
                                         size*sizeof(realw)),5205);
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_beta_kl_crust_mantle),
                                         size*sizeof(realw)),5206);
      print_CUDA_error_if_any(cudaMemset(mp->d_alpha_kl_crust_mantle,0,size*sizeof(realw)),5208);
      print_CUDA_error_if_any(cudaMemset(mp->d_beta_kl_crust_mantle,0,size*sizeof(realw)),5209);
    }else{
      // anisotropic kernels
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_cijkl_kl_crust_mantle),
                                         21*size*sizeof(realw)),5206);
      print_CUDA_error_if_any(cudaMemset(mp->d_cijkl_kl_crust_mantle,0,size*sizeof(realw)),5209);
    }

    // preconditioner
    if( mp->approximate_hess_kl ){
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_hess_kl_crust_mantle),
                                         size*sizeof(realw)),3030);
      print_CUDA_error_if_any(cudaMemset(mp->d_hess_kl_crust_mantle,0,size*sizeof(realw)),3031);
    }
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_crust_mantle_device");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// OUTER CORE

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_outer_core_device,
              PREPARE_OUTER_CORE_DEVICE)(long* Mesh_pointer_f,
                                         realw* h_xix, realw* h_xiy, realw* h_xiz,
                                         realw* h_etax, realw* h_etay, realw* h_etaz,
                                         realw* h_gammax, realw* h_gammay, realw* h_gammaz,
                                         realw* h_rho, realw* h_kappav,
                                         realw* h_rmass,
					 realw* h_normal_top_outer_core,
					 realw* h_normal_bottom_outer_core,
					 realw* h_jacobian2D_top_outer_core,
					 realw* h_jacobian2D_bottom_outer_core,
					 int* h_ibelm_top_outer_core,
					 int* h_ibelm_bottom_outer_core,
                                         int* h_ibool,
                                         realw* h_xstore, realw* h_ystore, realw* h_zstore,
                                         int* num_phase_ispec,
                                         int* phase_ispec_inner,
                                         int* nspec_outer,
                                         int* nspec_inner,
					 int* NSPEC2D_TOP_OC,
					 int* NSPEC2D_BOTTOM_OC
                                         ) {

  TRACE("prepare_outer_core_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  /* Assuming NGLLX=5. Padded is then 128 (5^3+3) */
  int size_padded = NGLL3_PADDED * (mp->NSPEC_OUTER_CORE);
  int size_glob = mp->NGLOB_OUTER_CORE;

  // mesh
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xix_outer_core, size_padded*sizeof(realw)),1001);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xiy_outer_core, size_padded*sizeof(realw)),1002);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xiz_outer_core, size_padded*sizeof(realw)),1003);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etax_outer_core, size_padded*sizeof(realw)),1004);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etay_outer_core, size_padded*sizeof(realw)),1005);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etaz_outer_core, size_padded*sizeof(realw)),1006);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammax_outer_core, size_padded*sizeof(realw)),1007);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammay_outer_core, size_padded*sizeof(realw)),1008);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammaz_outer_core, size_padded*sizeof(realw)),1009);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_kappavstore_outer_core, size_padded*sizeof(realw)),1010);

  // transfer constant element data with padding
  for(int i=0;i < mp->NSPEC_OUTER_CORE;i++) {
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xix_outer_core + i*NGLL3_PADDED, &h_xix[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1501);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xiy_outer_core+i*NGLL3_PADDED,   &h_xiy[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1502);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xiz_outer_core+i*NGLL3_PADDED,   &h_xiz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1503);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etax_outer_core+i*NGLL3_PADDED,  &h_etax[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1504);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etay_outer_core+i*NGLL3_PADDED,  &h_etay[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1505);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etaz_outer_core+i*NGLL3_PADDED,  &h_etaz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1506);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammax_outer_core+i*NGLL3_PADDED,&h_gammax[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1507);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammay_outer_core+i*NGLL3_PADDED,&h_gammay[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1508);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammaz_outer_core+i*NGLL3_PADDED,&h_gammaz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1509);

    print_CUDA_error_if_any(cudaMemcpy(mp->d_kappavstore_outer_core+i*NGLL3_PADDED,&h_kappav[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1510);
  }

  // needed for kernel calculations
  if( mp->simulation_type == 3 ){
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_rhostore_outer_core, size_padded*sizeof(realw)),1010);
    for(int i=0;i < mp->NSPEC_OUTER_CORE;i++) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_rhostore_outer_core+i*NGLL3_PADDED, &h_rho[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),2106);
    }
  }

  // global indexing
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibool_outer_core, size_padded*sizeof(int)),1021);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_ibool_outer_core, h_ibool,
                                     NGLL3*(mp->NSPEC_OUTER_CORE)*sizeof(int),cudaMemcpyHostToDevice),1022);

  // mesh locations
  // always needed
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_xstore_outer_core),sizeof(realw)*size_glob),2005);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_xstore_outer_core,h_xstore,
                                     sizeof(realw)*size_glob,cudaMemcpyHostToDevice),2100);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_ystore_outer_core),sizeof(realw)*size_glob),2005);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_ystore_outer_core,h_ystore,
                                     sizeof(realw)*size_glob,cudaMemcpyHostToDevice),2100);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_zstore_outer_core),sizeof(realw)*size_glob),2005);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_zstore_outer_core,h_zstore,
                                     sizeof(realw)*size_glob,cudaMemcpyHostToDevice),2100);

  // inner/outer elements
  mp->num_phase_ispec_outer_core = *num_phase_ispec;
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_phase_ispec_inner_outer_core),
                                     mp->num_phase_ispec_outer_core*2*sizeof(int)),2008);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_phase_ispec_inner_outer_core,phase_ispec_inner,
                                     mp->num_phase_ispec_outer_core*2*sizeof(int),cudaMemcpyHostToDevice),2101);

  mp->nspec_outer_outer_core = *nspec_outer;
  mp->nspec_inner_outer_core = *nspec_inner;

  // CMB/ICB coupling
  mp->nspec2D_top_outer_core = *NSPEC2D_TOP_OC;
  mp->nspec2D_bottom_outer_core = *NSPEC2D_BOTTOM_OC;
  int size_toc = NGLL2*(mp->nspec2D_top_outer_core);
  int size_boc = NGLL2*(mp->nspec2D_bottom_outer_core);

  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_normal_top_outer_core),sizeof(realw)*NDIM*size_toc),40020);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_normal_bottom_outer_core),sizeof(realw)*NDIM*size_boc),40021);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_normal_top_outer_core,h_normal_top_outer_core,sizeof(realw)*NDIM*size_toc,cudaMemcpyHostToDevice),40030);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_normal_bottom_outer_core,h_normal_bottom_outer_core,sizeof(realw)*NDIM*size_boc,cudaMemcpyHostToDevice),40031);

  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_jacobian2D_top_outer_core),sizeof(realw)*size_toc),40022);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_jacobian2D_bottom_outer_core),sizeof(realw)*size_boc),40023);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_jacobian2D_top_outer_core,h_jacobian2D_top_outer_core,sizeof(realw)*size_toc,cudaMemcpyHostToDevice),40032);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_jacobian2D_bottom_outer_core,h_jacobian2D_bottom_outer_core,sizeof(realw)*size_boc,cudaMemcpyHostToDevice),40033);

  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_ibelm_top_outer_core),sizeof(int)*(mp->nspec2D_top_outer_core)),40024);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_ibelm_bottom_outer_core),sizeof(int)*(mp->nspec2D_bottom_outer_core)),40025);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_ibelm_top_outer_core,h_ibelm_top_outer_core,sizeof(int)*(mp->nspec2D_top_outer_core),cudaMemcpyHostToDevice),40034);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_ibelm_bottom_outer_core,h_ibelm_bottom_outer_core,sizeof(int)*(mp->nspec2D_bottom_outer_core),cudaMemcpyHostToDevice),40035);

  // wavefield
  int size = mp->NGLOB_OUTER_CORE;
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_displ_outer_core),sizeof(realw)*size),4001);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_veloc_outer_core),sizeof(realw)*size),4002);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_accel_outer_core),sizeof(realw)*size),4003);
  // backward/reconstructed wavefield
  if( mp->simulation_type == 3 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_displ_outer_core),sizeof(realw)*size),4001);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_veloc_outer_core),sizeof(realw)*size),4002);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_accel_outer_core),sizeof(realw)*size),4003);
  }

  // mass matrix
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rmass_outer_core),sizeof(realw)*size_glob),2005);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_rmass_outer_core,h_rmass,
                                     sizeof(realw)*size_glob,cudaMemcpyHostToDevice),2100);

  // kernels
  if( mp->simulation_type == 3 ){

    size = NGLL3*(mp->NSPEC_OUTER_CORE);

    // density kernel
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rho_kl_outer_core),
                                       size*sizeof(realw)),5204);
    print_CUDA_error_if_any(cudaMemset(mp->d_rho_kl_outer_core,0,size*sizeof(realw)),5207);

    // isotropic kernel
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_alpha_kl_outer_core),
                                        size*sizeof(realw)),5205);
    print_CUDA_error_if_any(cudaMemset(mp->d_alpha_kl_outer_core,0,size*sizeof(realw)),5208);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_outer_core_device");
#endif
}


/* ----------------------------------------------------------------------------------------------- */

// INNER CORE

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_inner_core_device,
              PREPARE_INNER_CORE_DEVICE)(long* Mesh_pointer_f,
					 realw* h_xix, realw* h_xiy, realw* h_xiz,
					 realw* h_etax, realw* h_etay, realw* h_etaz,
					 realw* h_gammax, realw* h_gammay, realw* h_gammaz,
					 realw* h_rho, realw* h_kappav, realw* h_muv,
					 realw* h_rmass,
					 int* h_ibelm_top_inner_core,
					 int* h_ibool,
					 realw* h_xstore, realw* h_ystore, realw* h_zstore,
					 realw *c11store,realw *c12store,realw *c13store,
					 realw *c33store,realw *c44store,
					 int* h_idoubling_inner_core,
					 int* num_phase_ispec,
					 int* phase_ispec_inner,
					 int* nspec_outer,
					 int* nspec_inner,
					 int* NSPEC2D_TOP_IC) {
  
  TRACE("prepare_inner_core_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  /* Assuming NGLLX=5. Padded is then 128 (5^3+3) */
  int size_padded = NGLL3_PADDED * (mp->NSPEC_INNER_CORE);
  int size_glob = mp->NGLOB_INNER_CORE;

  // mesh
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xix_inner_core, size_padded*sizeof(realw)),1001);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xiy_inner_core, size_padded*sizeof(realw)),1002);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_xiz_inner_core, size_padded*sizeof(realw)),1003);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etax_inner_core, size_padded*sizeof(realw)),1004);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etay_inner_core, size_padded*sizeof(realw)),1005);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_etaz_inner_core, size_padded*sizeof(realw)),1006);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammax_inner_core, size_padded*sizeof(realw)),1007);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammay_inner_core, size_padded*sizeof(realw)),1008);
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_gammaz_inner_core, size_padded*sizeof(realw)),1009);

  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_muvstore_inner_core, size_padded*sizeof(realw)),1011);

  // transfer constant element data with padding
  for(int i=0;i < mp->NSPEC_INNER_CORE;i++) {
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xix_inner_core + i*NGLL3_PADDED, &h_xix[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1501);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xiy_inner_core+i*NGLL3_PADDED,   &h_xiy[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1502);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xiz_inner_core+i*NGLL3_PADDED,   &h_xiz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1503);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etax_inner_core+i*NGLL3_PADDED,  &h_etax[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1504);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etay_inner_core+i*NGLL3_PADDED,  &h_etay[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1505);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_etaz_inner_core+i*NGLL3_PADDED,  &h_etaz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1506);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammax_inner_core+i*NGLL3_PADDED,&h_gammax[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1507);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammay_inner_core+i*NGLL3_PADDED,&h_gammay[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1508);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammaz_inner_core+i*NGLL3_PADDED,&h_gammaz[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1509);

    print_CUDA_error_if_any(cudaMemcpy(mp->d_muvstore_inner_core+i*NGLL3_PADDED,   &h_muv[i*NGLL3],
                                       NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1511);
  }

  // anisotropy
  if( ! mp->anisotropic_inner_core ){
    // no anisotropy (uses kappav and muv in inner core)
    // kappavstore needed
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_kappavstore_inner_core, size_padded*sizeof(realw)),1010);
    for(int i=0;i < mp->NSPEC_INNER_CORE;i++) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_kappavstore_inner_core+i*NGLL3_PADDED,&h_kappav[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1510);
    }
  }else{
    // anisotropic inner core
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c11store_inner_core),
                                     size_padded*sizeof(realw)),4700);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c12store_inner_core),
                                     size_padded*sizeof(realw)),4701);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c13store_inner_core),
                                     size_padded*sizeof(realw)),4702);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c33store_inner_core),
                                     size_padded*sizeof(realw)),4703);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c44store_inner_core),
                                     size_padded*sizeof(realw)),4704);

    // transfer constant element data with padding
    for(int i=0;i < mp->NSPEC_INNER_CORE;i++) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c11store_inner_core + i*NGLL3_PADDED, &c11store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4800);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c12store_inner_core + i*NGLL3_PADDED, &c12store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4801);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c13store_inner_core + i*NGLL3_PADDED, &c13store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4802);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c33store_inner_core + i*NGLL3_PADDED, &c33store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4803);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c44store_inner_core + i*NGLL3_PADDED, &c44store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4804);
    }
  }

  // needed for boundary kernel calculations
  if( mp->simulation_type == 3 && mp->save_boundary_mesh ){
    print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_rhostore_inner_core, size_padded*sizeof(realw)),1010);
    for(int i=0;i < mp->NSPEC_INNER_CORE;i++) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_rhostore_inner_core+i*NGLL3_PADDED, &h_rho[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),2106);
    }
  }

  // global indexing
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_ibool_inner_core, size_padded*sizeof(int)),1021);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_ibool_inner_core, h_ibool,
                                     NGLL3*(mp->NSPEC_INNER_CORE)*sizeof(int),cudaMemcpyHostToDevice),1022);

  // fictious element flags
  print_CUDA_error_if_any(cudaMalloc((void**) &mp->d_idoubling_inner_core,
                                     mp->NSPEC_INNER_CORE*sizeof(int)),2010);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_idoubling_inner_core, h_idoubling_inner_core,
                                     mp->NSPEC_INNER_CORE*sizeof(int),cudaMemcpyHostToDevice),2011);

  // mesh locations
  // only needed when gravity is on
  if( mp->gravity ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_xstore_inner_core),sizeof(realw)*size_glob),2005);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_xstore_inner_core,h_xstore,
                                       sizeof(realw)*size_glob,cudaMemcpyHostToDevice),2100);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_ystore_inner_core),sizeof(realw)*size_glob),2005);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_ystore_inner_core,h_ystore,
                                       sizeof(realw)*size_glob,cudaMemcpyHostToDevice),2100);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_zstore_inner_core),sizeof(realw)*size_glob),2005);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_zstore_inner_core,h_zstore,
                                       sizeof(realw)*size_glob,cudaMemcpyHostToDevice),2100);
  }

  // inner/outer elements
  mp->num_phase_ispec_inner_core = *num_phase_ispec;
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_phase_ispec_inner_inner_core),
                                     mp->num_phase_ispec_inner_core*2*sizeof(int)),2008);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_phase_ispec_inner_inner_core,phase_ispec_inner,
                                     mp->num_phase_ispec_inner_core*2*sizeof(int),cudaMemcpyHostToDevice),2101);

  mp->nspec_outer_inner_core = *nspec_outer;
  mp->nspec_inner_inner_core = *nspec_inner;
  mp->nspec2D_top_inner_core = *NSPEC2D_TOP_IC;

  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_ibelm_top_inner_core),sizeof(int)*(mp->nspec2D_top_inner_core)),40021);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_ibelm_top_inner_core,h_ibelm_top_inner_core,sizeof(int)*(mp->nspec2D_top_inner_core),cudaMemcpyHostToDevice),40031);

  // wavefield
  int size = NDIM * mp->NGLOB_INNER_CORE;
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_displ_inner_core),sizeof(realw)*size),4001);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_veloc_inner_core),sizeof(realw)*size),4002);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_accel_inner_core),sizeof(realw)*size),4003);
  // backward/reconstructed wavefield
  if( mp->simulation_type == 3 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_displ_inner_core),sizeof(realw)*size),4001);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_veloc_inner_core),sizeof(realw)*size),4002);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_accel_inner_core),sizeof(realw)*size),4003);
  }

  // mass matrix
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rmass_inner_core),sizeof(realw)*size_glob),2005);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_rmass_inner_core,h_rmass,
                                     sizeof(realw)*size_glob,cudaMemcpyHostToDevice),2100);

  // kernels
  if( mp->simulation_type == 3 ){

    size = NGLL3*(mp->NSPEC_INNER_CORE);

    // density kernel
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rho_kl_inner_core),
                                       size*sizeof(realw)),5204);
    print_CUDA_error_if_any(cudaMemset(mp->d_rho_kl_inner_core,0,size*sizeof(realw)),5207);

    // isotropic kernel
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_alpha_kl_inner_core),
                                       size*sizeof(realw)),5205);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_beta_kl_inner_core),
                                       size*sizeof(realw)),5205);
    print_CUDA_error_if_any(cudaMemset(mp->d_alpha_kl_inner_core,0,size*sizeof(realw)),5208);
    print_CUDA_error_if_any(cudaMemset(mp->d_beta_kl_inner_core,0,size*sizeof(realw)),5208);
  }


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_inner_core_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// OCEANS

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_oceans_device,
              PREPARE_OCEANS_DEVICE)(long* Mesh_pointer_f,
             realw* h_rmass_ocean_load) {

  TRACE("prepare_oceans_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // mass matrix
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rmass_ocean_load),
                                     sizeof(realw)*mp->NGLOB_CRUST_MANTLE_OCEANS),4501);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_rmass_ocean_load,h_rmass_ocean_load,
                                     sizeof(realw)*mp->NGLOB_CRUST_MANTLE_OCEANS,cudaMemcpyHostToDevice),4502);

  // temporary global array: used to synchronize updates on global accel array
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_updated_dof_ocean_load),
                                     sizeof(int)*mp->NGLOB_CRUST_MANTLE_OCEANS),4502);


#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_oceans_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// for ELASTIC simulations

/* ----------------------------------------------------------------------------------------------- */

/*
extern "C"
void FC_FUNC_(prepare_fields_elastic_device,
              PREPARE_FIELDS_ELASTIC_DEVICE)(long* Mesh_pointer_f,
                                             int* size,
                                             realw* rmass,
                                             realw* rho_vp,
                                             realw* rho_vs,
                                             int* num_phase_ispec_elastic,
                                             int* phase_ispec_inner_elastic,
                                             int* ispec_is_elastic,
                                             int* ABSORBING_CONDITIONS,
                                             realw* h_b_absorb_field,
                                             int* h_b_reclen_field,
                                             int* SIMULATION_TYPE,int* SAVE_FORWARD,
                                             int* COMPUTE_AND_STORE_STRAIN,
                                             realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                             realw* epsilondev_xz,realw* epsilondev_yz,
                                             int* ATTENUATION,
                                             int* R_size,
                                             realw* R_xx,realw* R_yy,realw* R_xy,realw* R_xz,realw* R_yz,
                                             realw* one_minus_sum_beta,realw* factor_common,
                                             realw* alphaval,realw* betaval,realw* gammaval,
                                             int* OCEANS,
                                             realw* rmass_ocean_load,
                                             int* NOISE_TOMOGRAPHY,
                                             realw* free_surface_normal,
                                             int* free_surface_ispec,
                                             int* free_surface_ijk,
                                             int* num_free_surface_faces,
                                             int* ACOUSTIC_SIMULATION,
                                             int* num_colors_outer_elastic,
                                             int* num_colors_inner_elastic,
                                             int* num_elem_colors_elastic,
                                             int* ANISOTROPY,
                                             realw *c11store,
                                             realw *c12store,
                                             realw *c13store,
                                             realw *c14store,
                                             realw *c15store,
                                             realw *c16store,
                                             realw *c22store,
                                             realw *c23store,
                                             realw *c24store,
                                             realw *c25store,
                                             realw *c26store,
                                             realw *c33store,
                                             realw *c34store,
                                             realw *c35store,
                                             realw *c36store,
                                             realw *c44store,
                                             realw *c45store,
                                             realw *c46store,
                                             realw *c55store,
                                             realw *c56store,
                                             realw *c66store){

TRACE("prepare_fields_elastic_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);
  // Assuming NGLLX==5. Padded is then 128 (5^3+3)
  int size_padded = NGLL3_PADDED * (mp->NSPEC_AB);
  int size_nonpadded = NGLL3 * (mp->NSPEC_AB);

  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_displ),sizeof(realw)*(*size)),4001);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_veloc),sizeof(realw)*(*size)),4002);
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_accel),sizeof(realw)*(*size)),4003);

  // mpi buffer
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_send_accel_buffer),
                        3*(mp->max_nibool_interfaces_ext_mesh)*(mp->num_interfaces_ext_mesh)*sizeof(realw)),4004);

  // mass matrix
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rmass),sizeof(realw)*mp->NGLOB_AB),4005);
  // transfer element data
  print_CUDA_error_if_any(cudaMemcpy(mp->d_rmass,rmass,
                                     sizeof(realw)*mp->NGLOB_AB,cudaMemcpyHostToDevice),4010);


  // element indices
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_ispec_is_elastic),mp->NSPEC_AB*sizeof(int)),4009);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_ispec_is_elastic,ispec_is_elastic,
                                     mp->NSPEC_AB*sizeof(int),cudaMemcpyHostToDevice),4012);

  // phase elements
  mp->num_phase_ispec_elastic = *num_phase_ispec_elastic;
  print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_phase_ispec_inner_elastic),
                                     mp->num_phase_ispec_elastic*2*sizeof(int)),4008);
  print_CUDA_error_if_any(cudaMemcpy(mp->d_phase_ispec_inner_elastic,phase_ispec_inner_elastic,
                                     mp->num_phase_ispec_elastic*2*sizeof(int),cudaMemcpyHostToDevice),4011);

  // for seismograms
  if( mp->nrec_local > 0 ){
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_station_seismo_field),
                                     3*NGLL3*(mp->nrec_local)*sizeof(realw)),4015);

    mp->h_station_seismo_field = (realw*) malloc( 3*NGLL3*(mp->nrec_local)*sizeof(realw) );
    if( mp->h_station_seismo_field == NULL) exit_on_error("h_station_seismo_field not allocated \n");
  }

  // absorbing conditions
  if( *ABSORBING_CONDITIONS && mp->d_num_abs_boundary_faces > 0){
    // non-padded arrays
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rho_vp),size_nonpadded*sizeof(realw)),4006);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rho_vs),size_nonpadded*sizeof(realw)),4007);

    // rho_vp, rho_vs non-padded; they are needed for stacey boundary condition
    print_CUDA_error_if_any(cudaMemcpy(mp->d_rho_vp, rho_vp,
                                       NGLL3*mp->NSPEC_AB*sizeof(realw),cudaMemcpyHostToDevice),4013);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_rho_vs, rho_vs,
                                       NGLL3*mp->NSPEC_AB*sizeof(realw),cudaMemcpyHostToDevice),4014);

    // absorb_field array used for file i/o
    if(*SIMULATION_TYPE == 3 || ( *SIMULATION_TYPE == 1 && *SAVE_FORWARD )){
      mp->d_b_reclen_field = *h_b_reclen_field;
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_b_absorb_field),
                                       mp->d_b_reclen_field),4016);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_absorb_field, h_b_absorb_field,
                                       mp->d_b_reclen_field,cudaMemcpyHostToDevice),4017);
    }
  }

  // strains used for attenuation and kernel simulations
  if( *COMPUTE_AND_STORE_STRAIN ){
    // strains
    int epsilondev_size = NGLL3*mp->NSPEC_AB; // note: non-aligned; if align, check memcpy below and indexing

    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_epsilondev_xx,
                                       epsilondev_size*sizeof(realw)),4301);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_xx,epsilondev_xx,epsilondev_size*sizeof(realw),
                                       cudaMemcpyHostToDevice),4302);
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_epsilondev_yy,
                                       epsilondev_size*sizeof(realw)),4302);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_yy,epsilondev_yy,epsilondev_size*sizeof(realw),
                                       cudaMemcpyHostToDevice),4303);
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_epsilondev_xy,
                                       epsilondev_size*sizeof(realw)),4304);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_xy,epsilondev_xy,epsilondev_size*sizeof(realw),
                                       cudaMemcpyHostToDevice),4305);
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_epsilondev_xz,
                                       epsilondev_size*sizeof(realw)),4306);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_xz,epsilondev_xz,epsilondev_size*sizeof(realw),
                                       cudaMemcpyHostToDevice),4307);
    print_CUDA_error_if_any(cudaMalloc((void**)&mp->d_epsilondev_yz,
                                       epsilondev_size*sizeof(realw)),4308);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_epsilondev_yz,epsilondev_yz,epsilondev_size*sizeof(realw),
                                       cudaMemcpyHostToDevice),4309);

  }

  // attenuation memory variables
  if( *ATTENUATION ){
    // memory arrays
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_R_xx),
                                       (*R_size)*sizeof(realw)),4401);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_xx,R_xx,(*R_size)*sizeof(realw),
                                       cudaMemcpyHostToDevice),4402);

    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_R_yy),
                                       (*R_size)*sizeof(realw)),4403);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_yy,R_yy,(*R_size)*sizeof(realw),
                                       cudaMemcpyHostToDevice),4404);

    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_R_xy),
                                       (*R_size)*sizeof(realw)),4405);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_xy,R_xy,(*R_size)*sizeof(realw),
                                       cudaMemcpyHostToDevice),4406);

    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_R_xz),
                                       (*R_size)*sizeof(realw)),4407);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_xz,R_xz,(*R_size)*sizeof(realw),
                                       cudaMemcpyHostToDevice),4408);

    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_R_yz),
                                       (*R_size)*sizeof(realw)),4409);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_R_yz,R_yz,(*R_size)*sizeof(realw),
                                       cudaMemcpyHostToDevice),4410);

    // attenuation factors
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_one_minus_sum_beta),
                                       NGLL3*mp->NSPEC_AB*sizeof(realw)),4430);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_one_minus_sum_beta ,one_minus_sum_beta,
                                       NGLL3*mp->NSPEC_AB*sizeof(realw),cudaMemcpyHostToDevice),4431);

    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_factor_common),
                                       N_SLS*NGLL3*mp->NSPEC_AB*sizeof(realw)),4432);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_factor_common ,factor_common,
                                       N_SLS*NGLL3*mp->NSPEC_AB*sizeof(realw),cudaMemcpyHostToDevice),4433);

    // alpha,beta,gamma factors
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_alphaval),
                                       N_SLS*sizeof(realw)),4434);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_alphaval ,alphaval,
                                       N_SLS*sizeof(realw),cudaMemcpyHostToDevice),4435);

    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_betaval),
                                       N_SLS*sizeof(realw)),4436);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_betaval ,betaval,
                                       N_SLS*sizeof(realw),cudaMemcpyHostToDevice),4437);

    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_gammaval),
                                       N_SLS*sizeof(realw)),4438);
    print_CUDA_error_if_any(cudaMemcpy(mp->d_gammaval ,gammaval,
                                       N_SLS*sizeof(realw),cudaMemcpyHostToDevice),4439);

  }

  // anisotropy
  if( *ANISOTROPY ){
    // allocates memory on GPU
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c11store),
                                       size_padded*sizeof(realw)),4700);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c12store),
                                       size_padded*sizeof(realw)),4701);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c13store),
                                       size_padded*sizeof(realw)),4702);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c14store),
                                       size_padded*sizeof(realw)),4703);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c15store),
                                       size_padded*sizeof(realw)),4704);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c16store),
                                       size_padded*sizeof(realw)),4705);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c22store),
                                       size_padded*sizeof(realw)),4706);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c23store),
                                       size_padded*sizeof(realw)),4707);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c24store),
                                       size_padded*sizeof(realw)),4708);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c25store),
                                       size_padded*sizeof(realw)),4709);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c26store),
                                       size_padded*sizeof(realw)),4710);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c33store),
                                       size_padded*sizeof(realw)),4711);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c34store),
                                       size_padded*sizeof(realw)),4712);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c35store),
                                       size_padded*sizeof(realw)),4713);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c36store),
                                       size_padded*sizeof(realw)),4714);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c44store),
                                       size_padded*sizeof(realw)),4715);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c45store),
                                       size_padded*sizeof(realw)),4716);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c46store),
                                       size_padded*sizeof(realw)),4717);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c55store),
                                       size_padded*sizeof(realw)),4718);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c56store),
                                       size_padded*sizeof(realw)),4719);
    print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_c66store),
                                       size_padded*sizeof(realw)),4720);

    // transfer constant element data with padding
    for(int i=0;i < mp->NSPEC_AB;i++) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c11store + i*NGLL3_PADDED, &c11store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4800);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c12store + i*NGLL3_PADDED, &c12store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4801);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c13store + i*NGLL3_PADDED, &c13store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4802);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c14store + i*NGLL3_PADDED, &c14store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4803);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c15store + i*NGLL3_PADDED, &c15store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4804);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c16store + i*NGLL3_PADDED, &c16store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4805);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c22store + i*NGLL3_PADDED, &c22store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4806);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c23store + i*NGLL3_PADDED, &c23store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4807);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c24store + i*NGLL3_PADDED, &c24store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4808);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c25store + i*NGLL3_PADDED, &c25store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4809);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c26store + i*NGLL3_PADDED, &c26store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4810);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c33store + i*NGLL3_PADDED, &c33store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4811);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c34store + i*NGLL3_PADDED, &c34store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4812);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c35store + i*NGLL3_PADDED, &c35store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4813);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c36store + i*NGLL3_PADDED, &c36store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4814);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c44store + i*NGLL3_PADDED, &c44store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4815);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c45store + i*NGLL3_PADDED, &c45store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4816);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c46store + i*NGLL3_PADDED, &c46store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4817);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c55store + i*NGLL3_PADDED, &c55store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4818);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c56store + i*NGLL3_PADDED, &c56store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4819);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_c66store + i*NGLL3_PADDED, &c66store[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4820);
    }
  }

  // ocean load approximation
  if( *OCEANS ){
    // oceans needs a free surface
    mp->num_free_surface_faces = *num_free_surface_faces;
    if( mp->num_free_surface_faces > 0 ){
      // mass matrix
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_rmass_ocean_load),
                                         sizeof(realw)*mp->NGLOB_AB),4501);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_rmass_ocean_load,rmass_ocean_load,
                                         sizeof(realw)*mp->NGLOB_AB,cudaMemcpyHostToDevice),4502);
      // surface normal
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_free_surface_normal),
                                         3*NGLL2*(mp->num_free_surface_faces)*sizeof(realw)),4503);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_free_surface_normal,free_surface_normal,
                                         3*NGLL2*(mp->num_free_surface_faces)*sizeof(realw),cudaMemcpyHostToDevice),4504);

      // temporary global array: used to synchronize updates on global accel array
      print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_updated_dof_ocean_load),
                                         sizeof(int)*mp->NGLOB_AB),4505);

      if( *NOISE_TOMOGRAPHY == 0 && *ACOUSTIC_SIMULATION == 0 ){
        print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_free_surface_ispec),
                                          mp->num_free_surface_faces*sizeof(int)),4601);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_free_surface_ispec,free_surface_ispec,
                                          mp->num_free_surface_faces*sizeof(int),cudaMemcpyHostToDevice),4603);
        print_CUDA_error_if_any(cudaMalloc((void**)&(mp->d_free_surface_ijk),
                                          3*NGLL2*mp->num_free_surface_faces*sizeof(int)),4602);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_free_surface_ijk,free_surface_ijk,
                                          3*NGLL2*mp->num_free_surface_faces*sizeof(int),cudaMemcpyHostToDevice),4604);
      }
    }
  }

  // mesh coloring
  if( mp->use_mesh_coloring_gpu ){
    mp->num_colors_outer_elastic = *num_colors_outer_elastic;
    mp->num_colors_inner_elastic = *num_colors_inner_elastic;
    mp->h_num_elem_colors_elastic = (int*) num_elem_colors_elastic;
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("prepare_fields_elastic_device");
#endif
}
*/



/* ----------------------------------------------------------------------------------------------- */

// cleanup

/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(prepare_cleanup_device,
              PREPARE_CLEANUP_DEVICE)(long* Mesh_pointer_f,
				      int* NCHUNKS_VAL) {

TRACE("prepare_cleanup_device");

  // frees allocated memory arrays
  Mesh* mp = (Mesh*)(*Mesh_pointer_f);

  // frees memory on GPU

  //------------------------------------------
  // sources
  //------------------------------------------
  if( mp->simulation_type == 1  || mp->simulation_type == 3 ){
    cudaFree(mp->d_sourcearrays);
    cudaFree(mp->d_stf_pre_compute);
  }

  cudaFree(mp->d_islice_selected_source);
  cudaFree(mp->d_ispec_selected_source);

  //------------------------------------------
  // receivers
  //------------------------------------------
  if( mp->nrec_local > 0 ) {
    cudaFree(mp->d_number_receiver_global);
    cudaFree(mp->d_station_seismo_field);
    free(mp->h_station_seismo_field);
  }
  cudaFree(mp->d_ispec_selected_rec);

  if( mp->nadj_rec_local > 0 ){
    cudaFree(mp->d_adj_sourcearrays);
    cudaFree(mp->d_pre_computed_irec);
    free(mp->h_adj_sourcearrays_slice);
  }

  //------------------------------------------
  // rotation arrays
  //------------------------------------------
  if( mp->rotation ){
    cudaFree(mp->d_A_array_rotation);
    cudaFree(mp->d_B_array_rotation);
    if( mp->simulation_type == 3 ){
      cudaFree(mp->d_b_A_array_rotation);
      cudaFree(mp->d_b_B_array_rotation);
    }
  }

  //------------------------------------------
  // gravity arrays
  //------------------------------------------
  if( ! mp->gravity ){
    cudaFree(mp->d_d_ln_density_dr_table);
  }else{
    cudaFree(mp->d_minus_rho_g_over_kappa_fluid);
    cudaFree(mp->d_minus_gravity_table);
    cudaFree(mp->d_minus_deriv_gravity_table);
    cudaFree(mp->d_density_table);
  }

  //------------------------------------------
  // attenuation arrays
  //------------------------------------------
  if( mp->attenuation ){
    cudaFree(mp->d_one_minus_sum_beta_crust_mantle);
    cudaFree(mp->d_one_minus_sum_beta_inner_core);
    if( ! mp->use_attenuation_mimic ){
      cudaFree(mp->d_factor_common_crust_mantle);
      cudaFree(mp->d_R_xx_crust_mantle);
      cudaFree(mp->d_R_yy_crust_mantle);
      cudaFree(mp->d_R_xy_crust_mantle);
      cudaFree(mp->d_R_xz_crust_mantle);
      cudaFree(mp->d_R_yz_crust_mantle);
      cudaFree(mp->d_factor_common_inner_core);
      cudaFree(mp->d_R_xx_inner_core);
      cudaFree(mp->d_R_yy_inner_core);
      cudaFree(mp->d_R_xy_inner_core);
      cudaFree(mp->d_R_xz_inner_core);
      cudaFree(mp->d_R_yz_inner_core);
    }
    cudaFree(mp->d_alphaval);
    cudaFree(mp->d_betaval);
    cudaFree(mp->d_gammaval);
    if( mp->simulation_type == 3 ){
      cudaFree(mp->d_b_alphaval);
      cudaFree(mp->d_b_betaval);
      cudaFree(mp->d_b_gammaval);
    }
  }

  //------------------------------------------
  // strain
  //------------------------------------------
  if( mp->compute_and_store_strain ){
    cudaFree(mp->d_epsilondev_xx_crust_mantle);
    cudaFree(mp->d_epsilondev_yy_crust_mantle);
    cudaFree(mp->d_epsilondev_xy_crust_mantle);
    cudaFree(mp->d_epsilondev_xz_crust_mantle);
    cudaFree(mp->d_epsilondev_yz_crust_mantle);

    cudaFree(mp->d_epsilondev_xx_inner_core);
    cudaFree(mp->d_epsilondev_yy_inner_core);
    cudaFree(mp->d_epsilondev_xy_inner_core);
    cudaFree(mp->d_epsilondev_xz_inner_core);
    cudaFree(mp->d_epsilondev_yz_inner_core);

    cudaFree(mp->d_eps_trace_over_3_crust_mantle);
    cudaFree(mp->d_eps_trace_over_3_inner_core);
    if( mp->simulation_type == 3 ){
      cudaFree(mp->d_b_epsilondev_xx_crust_mantle);
      cudaFree(mp->d_b_epsilondev_yy_crust_mantle);
      cudaFree(mp->d_b_epsilondev_xy_crust_mantle);
      cudaFree(mp->d_b_epsilondev_xz_crust_mantle);
      cudaFree(mp->d_b_epsilondev_yz_crust_mantle);

      cudaFree(mp->d_b_epsilondev_xx_inner_core);
      cudaFree(mp->d_b_epsilondev_yy_inner_core);
      cudaFree(mp->d_b_epsilondev_xy_inner_core);
      cudaFree(mp->d_b_epsilondev_xz_inner_core);
      cudaFree(mp->d_b_epsilondev_yz_inner_core);

      cudaFree(mp->d_b_eps_trace_over_3_crust_mantle);
      cudaFree(mp->d_b_eps_trace_over_3_inner_core);
    }
  }

  //------------------------------------------
  // absorbing boundaries arrays
  //------------------------------------------
  if( mp->absorbing_conditions){
    cudaFree(mp->d_rho_vp_crust_mantle);
    cudaFree(mp->d_rho_vs_crust_mantle);
    cudaFree(mp->d_nkmin_xi_crust_mantle);
    cudaFree(mp->d_nkmin_eta_crust_mantle);
    cudaFree(mp->d_njmin_crust_mantle);
    cudaFree(mp->d_njmax_crust_mantle);
    cudaFree(mp->d_nimin_crust_mantle);
    cudaFree(mp->d_nimax_crust_mantle);
    if( mp->nspec2D_xmin_crust_mantle > 0 ){
      cudaFree(mp->d_ibelm_xmin_crust_mantle);
      cudaFree(mp->d_normal_xmin_crust_mantle);
      cudaFree(mp->d_jacobian2D_xmin_crust_mantle);
      if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
        cudaFree(mp->d_absorb_xmin_crust_mantle);
      }
    }
    if( mp->nspec2D_xmax_crust_mantle > 0 ){
      cudaFree(mp->d_ibelm_xmax_crust_mantle);
      cudaFree(mp->d_normal_xmax_crust_mantle);
      cudaFree(mp->d_jacobian2D_xmax_crust_mantle);
      if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
        cudaFree(mp->d_absorb_xmax_crust_mantle);
      }
    }
    if( mp->nspec2D_ymin_crust_mantle > 0 ){
      cudaFree(mp->d_ibelm_ymin_crust_mantle);
      cudaFree(mp->d_normal_ymin_crust_mantle);
      cudaFree(mp->d_jacobian2D_ymin_crust_mantle);
      if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
        cudaFree(mp->d_absorb_ymin_crust_mantle);
      }
    }
    if( mp->nspec2D_ymax_crust_mantle > 0 ){
      cudaFree(mp->d_ibelm_ymax_crust_mantle);
      cudaFree(mp->d_normal_ymax_crust_mantle);
      cudaFree(mp->d_jacobian2D_ymax_crust_mantle);
      if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
        cudaFree(mp->d_absorb_ymax_crust_mantle);
      }
    }

    cudaFree(mp->d_vp_outer_core);
    cudaFree(mp->d_nkmin_xi_outer_core);
    cudaFree(mp->d_nkmin_eta_outer_core);
    cudaFree(mp->d_njmin_outer_core);
    cudaFree(mp->d_njmax_outer_core);
    cudaFree(mp->d_nimin_outer_core);
    cudaFree(mp->d_nimax_outer_core);
    if( mp->nspec2D_xmin_outer_core > 0 ){
      cudaFree(mp->d_ibelm_xmin_outer_core);
      cudaFree(mp->d_jacobian2D_xmin_outer_core);
      if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
        cudaFree(mp->d_absorb_xmin_outer_core);
      }
    }
    if( mp->nspec2D_xmax_outer_core > 0 ){
      cudaFree(mp->d_ibelm_xmax_outer_core);
      cudaFree(mp->d_jacobian2D_xmax_outer_core);
      if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
        cudaFree(mp->d_absorb_xmax_outer_core);
      }
    }
    if( mp->nspec2D_ymin_outer_core > 0 ){
      cudaFree(mp->d_ibelm_ymin_outer_core);
      cudaFree(mp->d_jacobian2D_ymin_outer_core);
      if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
        cudaFree(mp->d_absorb_ymin_outer_core);
      }
    }
    if( mp->nspec2D_ymax_outer_core > 0 ){
      cudaFree(mp->d_ibelm_ymax_outer_core);
      cudaFree(mp->d_jacobian2D_ymax_outer_core);
      if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
        cudaFree(mp->d_absorb_ymax_outer_core);
      }
    }
    if( mp->nspec2D_zmin_outer_core > 0 ){
      cudaFree(mp->d_ibelm_zmin_outer_core);
      cudaFree(mp->d_jacobian2D_zmin_outer_core);
      if( (mp->simulation_type == 1 && mp->save_forward ) || (mp->simulation_type == 3) ){
        cudaFree(mp->d_absorb_zmin_outer_core);
      }
    }

  }

  //------------------------------------------
  // mpi buffers
  //------------------------------------------
  if( mp->num_interfaces_crust_mantle > 0 ){
    cudaFree(mp->d_nibool_interfaces_crust_mantle);
    cudaFree(mp->d_ibool_interfaces_crust_mantle);
    cudaFree(mp->d_send_accel_buffer_crust_mantle);
  }
  if( mp->num_interfaces_inner_core > 0 ){
    cudaFree(mp->d_nibool_interfaces_inner_core);
    cudaFree(mp->d_ibool_interfaces_inner_core);
    cudaFree(mp->d_send_accel_buffer_inner_core);
  }
  if( mp->num_interfaces_outer_core > 0 ){
    cudaFree(mp->d_nibool_interfaces_outer_core);
    cudaFree(mp->d_ibool_interfaces_outer_core);
    cudaFree(mp->d_send_accel_buffer_outer_core);
  }

  //------------------------------------------
  // NOISE arrays
  //------------------------------------------
  if( mp->noise_tomography > 0 ){
    cudaFree(mp->d_ibelm_top_crust_mantle);
    cudaFree(mp->d_noise_surface_movie);
    if( mp->noise_tomography == 1 ) cudaFree(mp->d_noise_sourcearray);
    if( mp->noise_tomography > 1 ){
      cudaFree(mp->d_normal_x_noise);
      cudaFree(mp->d_normal_y_noise);
      cudaFree(mp->d_normal_z_noise);
      cudaFree(mp->d_mask_noise);
      cudaFree(mp->d_jacobian2D_top_crust_mantle);
    }
    if( mp->noise_tomography == 3 ) cudaFree(mp->d_Sigma_kl);
  }

  //------------------------------------------
  // crust_mantle
  //------------------------------------------
  cudaFree(mp->d_xix_crust_mantle);
  cudaFree(mp->d_xiy_crust_mantle);
  cudaFree(mp->d_xiz_crust_mantle);
  cudaFree(mp->d_etax_crust_mantle);
  cudaFree(mp->d_etay_crust_mantle);
  cudaFree(mp->d_etaz_crust_mantle);
  cudaFree(mp->d_gammax_crust_mantle);
  cudaFree(mp->d_gammay_crust_mantle);
  cudaFree(mp->d_gammaz_crust_mantle);

  cudaFree(mp->d_muvstore_crust_mantle);
  cudaFree(mp->d_ibool_crust_mantle);

  if( ! mp->anisotropic_3D_mantle ){
    cudaFree(mp->d_kappavstore_crust_mantle);
    cudaFree(mp->d_kappahstore_crust_mantle);
    cudaFree(mp->d_muhstore_crust_mantle);
    cudaFree(mp->d_eta_anisostore_crust_mantle);
    cudaFree(mp->d_ispec_is_tiso_crust_mantle);
  }else{
    cudaFree(mp->d_c11store_crust_mantle);
    cudaFree(mp->d_c12store_crust_mantle);
    cudaFree(mp->d_c13store_crust_mantle);
    cudaFree(mp->d_c14store_crust_mantle);
    cudaFree(mp->d_c15store_crust_mantle);
    cudaFree(mp->d_c16store_crust_mantle);
    cudaFree(mp->d_c22store_crust_mantle);
    cudaFree(mp->d_c23store_crust_mantle);
    cudaFree(mp->d_c24store_crust_mantle);
    cudaFree(mp->d_c25store_crust_mantle);
    cudaFree(mp->d_c26store_crust_mantle);
    cudaFree(mp->d_c33store_crust_mantle);
    cudaFree(mp->d_c34store_crust_mantle);
    cudaFree(mp->d_c35store_crust_mantle);
    cudaFree(mp->d_c36store_crust_mantle);
    cudaFree(mp->d_c44store_crust_mantle);
    cudaFree(mp->d_c45store_crust_mantle);
    cudaFree(mp->d_c46store_crust_mantle);
    cudaFree(mp->d_c55store_crust_mantle);
    cudaFree(mp->d_c56store_crust_mantle);
    cudaFree(mp->d_c66store_crust_mantle);
  }

  if( mp->simulation_type == 3 && mp->save_boundary_mesh ){
    cudaFree(mp->d_rhostore_crust_mantle);
  }

  cudaFree(mp->d_ystore_crust_mantle);
  cudaFree(mp->d_zstore_crust_mantle);
  if( mp->gravity ){
    cudaFree(mp->d_xstore_crust_mantle);
  }
  cudaFree(mp->d_phase_ispec_inner_crust_mantle);

  cudaFree(mp->d_normal_top_crust_mantle);
  cudaFree(mp->d_ibelm_top_crust_mantle);
  cudaFree(mp->d_ibelm_bottom_crust_mantle);

  cudaFree(mp->d_displ_crust_mantle);
  cudaFree(mp->d_veloc_crust_mantle);
  cudaFree(mp->d_accel_crust_mantle);
  if( mp->simulation_type == 3 ){
    cudaFree(mp->d_b_displ_crust_mantle);
    cudaFree(mp->d_b_veloc_crust_mantle);
    cudaFree(mp->d_b_accel_crust_mantle);
    cudaFree(mp->d_rho_kl_crust_mantle);
    if(mp->anisotropic_kl){
      cudaFree(mp->d_cijkl_kl_crust_mantle);
    }else{
      cudaFree(mp->d_alpha_kl_crust_mantle);
      cudaFree(mp->d_beta_kl_crust_mantle);
    }
    if(mp->approximate_hess_kl){ cudaFree(mp->d_hess_kl_crust_mantle);}
  }
  if( *NCHUNKS_VAL != 6 && mp->absorbing_conditions){
    cudaFree(mp->d_rmassx_crust_mantle);
    cudaFree(mp->d_rmassy_crust_mantle);
  }
  cudaFree(mp->d_rmassz_crust_mantle);

  //------------------------------------------
  // outer_core
  //------------------------------------------
  cudaFree(mp->d_xix_outer_core);
  cudaFree(mp->d_xiy_outer_core);
  cudaFree(mp->d_xiz_outer_core);
  cudaFree(mp->d_etax_outer_core);
  cudaFree(mp->d_etay_outer_core);
  cudaFree(mp->d_etaz_outer_core);
  cudaFree(mp->d_gammax_outer_core);
  cudaFree(mp->d_gammay_outer_core);
  cudaFree(mp->d_gammaz_outer_core);

  cudaFree(mp->d_kappavstore_outer_core);
  if( mp->simulation_type == 3 ){
    cudaFree(mp->d_rhostore_outer_core);
  }

  cudaFree(mp->d_xstore_outer_core);
  cudaFree(mp->d_ystore_outer_core);
  cudaFree(mp->d_zstore_outer_core);

  cudaFree(mp->d_ibool_outer_core);
  cudaFree(mp->d_phase_ispec_inner_outer_core);

  cudaFree(mp->d_ibelm_top_outer_core);
  cudaFree(mp->d_ibelm_bottom_outer_core);

  cudaFree(mp->d_normal_top_outer_core);
  cudaFree(mp->d_normal_bottom_outer_core);

  cudaFree(mp->d_jacobian2D_top_outer_core);
  cudaFree(mp->d_jacobian2D_bottom_outer_core);

  cudaFree(mp->d_displ_outer_core);
  cudaFree(mp->d_veloc_outer_core);
  cudaFree(mp->d_accel_outer_core);
  if( mp->simulation_type == 3 ){
    cudaFree(mp->d_b_displ_outer_core);
    cudaFree(mp->d_b_veloc_outer_core);
    cudaFree(mp->d_b_accel_outer_core);
    cudaFree(mp->d_rho_kl_outer_core);
    cudaFree(mp->d_alpha_kl_outer_core);
  }
  cudaFree(mp->d_rmass_outer_core);

  //------------------------------------------
  // inner_core
  //------------------------------------------
  cudaFree(mp->d_xix_inner_core);
  cudaFree(mp->d_xiy_inner_core);
  cudaFree(mp->d_xiz_inner_core);
  cudaFree(mp->d_etax_inner_core);
  cudaFree(mp->d_etay_inner_core);
  cudaFree(mp->d_etaz_inner_core);
  cudaFree(mp->d_gammax_inner_core);
  cudaFree(mp->d_gammay_inner_core);
  cudaFree(mp->d_gammaz_inner_core);


  cudaFree(mp->d_muvstore_inner_core);
  cudaFree(mp->d_ibool_inner_core);

  if( mp->oceans ){
    cudaFree(mp->d_rmass_ocean_load);
    cudaFree(mp->d_updated_dof_ocean_load);
  }

  if( mp->gravity ){
    cudaFree(mp->d_xstore_inner_core);
    cudaFree(mp->d_ystore_inner_core);
    cudaFree(mp->d_zstore_inner_core);
  }

  cudaFree(mp->d_ibelm_top_inner_core);

  if( ! mp->anisotropic_inner_core ){
    cudaFree(mp->d_kappavstore_inner_core);
  }else{
    cudaFree(mp->d_c11store_inner_core);
    cudaFree(mp->d_c12store_inner_core);
    cudaFree(mp->d_c13store_inner_core);
    cudaFree(mp->d_c33store_inner_core);
    cudaFree(mp->d_c44store_inner_core);
  }

  if( mp->simulation_type == 3 && mp->save_boundary_mesh ){
    cudaFree(mp->d_rhostore_inner_core);
  }
  cudaFree(mp->d_idoubling_inner_core);
  if( mp->gravity ){
    cudaFree(mp->d_xstore_inner_core);
    cudaFree(mp->d_ystore_inner_core);
    cudaFree(mp->d_zstore_inner_core);
  }
  cudaFree(mp->d_phase_ispec_inner_inner_core);

  cudaFree(mp->d_displ_inner_core);
  cudaFree(mp->d_veloc_inner_core);
  cudaFree(mp->d_accel_inner_core);
  if( mp->simulation_type == 3 ) {
    cudaFree(mp->d_b_displ_inner_core);
    cudaFree(mp->d_b_veloc_inner_core);
    cudaFree(mp->d_b_accel_inner_core);

    cudaFree(mp->d_rho_kl_inner_core);
    cudaFree(mp->d_alpha_kl_inner_core);
    cudaFree(mp->d_beta_kl_inner_core);
  }
  cudaFree(mp->d_rmass_inner_core);

/*

    if( *OCEANS ){
      if( mp->num_free_surface_faces > 0 ){
        cudaFree(mp->d_rmass_ocean_load);
        cudaFree(mp->d_free_surface_normal);
        cudaFree(mp->d_updated_dof_ocean_load);
        if( *NOISE_TOMOGRAPHY == 0){
          cudaFree(mp->d_free_surface_ispec);
          cudaFree(mp->d_free_surface_ijk);
        }
      }
    }
  } // ELASTIC_SIMULATION

*/

  // releases previous contexts
  cudaThreadExit();

  // mesh pointer - not needed anymore
  free(mp);
}

