/*
!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
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

#include "mesh_constants_gpu.h"

/* ----------------------------------------------------------------------------------------------- */
// OpenCL setup
/* ----------------------------------------------------------------------------------------------- */

#ifdef USE_OPENCL
#ifdef USE_TEXTURES_FIELDS
static cl_mem moclGetDummyImage2D (Mesh *mp) {
  static int inited = 0;
  static cl_mem image2d;
  cl_int errcode;

  if (!inited) {
    cl_image_format format = {CL_RGBA, CL_UNSIGNED_INT32};
    image2d = clCreateImage2D (mocl.context, CL_MEM_READ_ONLY, &format, 10, 10, 0, NULL, clck_(&errcode));
    inited = 1;
  } else {
    clRetainMemObject (image2d);
  }

  return image2d;
}
#endif
/* ----------------------------------------------------------------------------------------------- */

void release_kernels (void) {
#undef BOAST_KERNEL
#define BOAST_KERNEL(__kern_name__)                                     \
  clCheck (clReleaseKernel (mocl.kernels.__kern_name__));               \
  clCheck (clReleaseProgram (mocl.programs.__kern_name__ ## _program));

  #include "kernel_list.h"
}

#endif // USE_OPENCL


/*----------------------------------------------------------------------------------------------- */
// GPU preparation routines
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (prepare_constants_device,
               PREPARE_CONSTANTS_DEVICE) (long *Mesh_pointer,
                                          int *myrank_f,
                                          int *h_NGLLX,
                                          realw *h_hprime_xx, realw *h_hprimewgll_xx,
                                          realw *h_wgllwgll_xy, realw *h_wgllwgll_xz, realw *h_wgllwgll_yz,
                                          int *NSOURCES, int *nsources_local,
                                          realw *h_sourcearrays,
                                          int *h_islice_selected_source, int *h_ispec_selected_source,
                                          int *nrec, int *nrec_local, int *nadj_rec_local,
                                          int *h_number_receiver_global,
                                          int *h_islice_selected_rec, int *h_ispec_selected_rec,
                                          int *NSPEC_CRUST_MANTLE, int *NGLOB_CRUST_MANTLE,
                                          int *NSPEC_CRUST_MANTLE_STRAIN_ONLY,
                                          int *NSPEC_OUTER_CORE, int *NGLOB_OUTER_CORE,
                                          int *NSPEC_INNER_CORE, int *NGLOB_INNER_CORE,
                                          int *NSPEC_INNER_CORE_STRAIN_ONLY,
                                          int *SIMULATION_TYPE, int *NOISE_TOMOGRAPHY,
                                          int *SAVE_FORWARD_f, int *ABSORBING_CONDITIONS_f,
                                          int *OCEANS_f, int *GRAVITY_f,
                                          int *ROTATION_f, int *EXACT_MASS_MATRIX_FOR_ROTATION_f,
                                          int *ATTENUATION_f, int *UNDO_ATTENUATION_f,
                                          int *PARTIAL_PHYS_DISPERSION_ONLY_f, int *USE_3D_ATTENUATION_ARRAYS_f,
                                          int *COMPUTE_AND_STORE_STRAIN_f,
                                          int *ANISOTROPIC_3D_MANTLE_f, int *ANISOTROPIC_INNER_CORE_f,
                                          int *SAVE_BOUNDARY_MESH_f,
                                          int *USE_MESH_COLORING_GPU_f,
                                          int *ANISOTROPIC_KL_f, int *APPROXIMATE_HESS_KL_f,
                                          realw *deltat_f, realw *b_deltat_f,
                                          int *GPU_ASYNC_COPY_f) {

  TRACE ("prepare_constants_device");

  // allocates mesh parameter structure
  Mesh *mp = (Mesh *) malloc (sizeof (Mesh));

  // checks pointer
  if (! mp) exit_on_error ("Error allocating mesh pointer");

  // sets fortran pointer
  *Mesh_pointer = (long) mp;

  // safety checks
  // constants defined in constants.h and mesh_constants_gpu.h have to match
  if (*h_NGLLX != NGLLX) { exit_on_error ("NGLLX must be set to 5 for GPU devices in constants.h; please re-compile"); }
  if (*GPU_ASYNC_COPY_f) {
    if (! GPU_ASYNC_COPY) { exit_on_error ("GPU_ASYNC_COPY must be set to 1 for GPU devices in mesh_constants_gpu.h; please re-compile"); }
  }else{
    if (GPU_ASYNC_COPY) { exit_on_error ("GPU_ASYNC_COPY must be set to 0 for GPU devices in mesh_constants_gpu.h; please re-compile"); }
  }
#ifdef USE_MESH_COLORING_GPU
  if (! *USE_MESH_COLORING_GPU_f) { exit_on_error("Error with USE_MESH_COLORING_GPU constant; please re-compile\n"); }
#endif

  // initializes gpu array pointers
  gpuInitialize_buffers(mp);

  // MPI process rank
  mp->myrank = *myrank_f;

#ifdef USE_CUDA
  if (run_cuda) {
    // setup two streams, one for compute and one for host<->device memory copies
    // uses pinned memory for asynchronous data transfers
    // note: creating streams may fail if multiple processes use a single GPU, and when CUDA devices are set
    //       for exclusive usage. in this case, you will have to setup a CUDA multiple process service (MPS), aka CUDA proxy.
    //       see: https://docs.nvidia.com/deploy/pdf/CUDA_Multi_Process_Service_Overview.pdf
    //            http://cudamusing.blogspot.ch/2013/07/enabling-cuda-multi-process-service-mps.html
    //
    //       for some reason, calling cudaSetDevice() may return success even when devices are exclusive,
    //       and only when calling the cudaStreamCreate() below a cudaError occurs...
    //
    // compute stream
    print_CUDA_error_if_any( cudaStreamCreate(&mp->compute_stream), 101);
    // copy stream (needed to transfer MPI buffers)
    if (GPU_ASYNC_COPY) print_CUDA_error_if_any( cudaStreamCreate(&mp->copy_stream), 102);
  }
#endif

  // sets constant arrays
  gpuSetConst (&mp->d_hprime_xx, NGLL2, h_hprime_xx);
  gpuSetConst (&mp->d_hprimewgll_xx, NGLL2, h_hprimewgll_xx);

  gpuSetConst (&mp->d_wgllwgll_xy, NGLL2, h_wgllwgll_xy);
  gpuSetConst (&mp->d_wgllwgll_xz, NGLL2, h_wgllwgll_xz);
  gpuSetConst (&mp->d_wgllwgll_yz, NGLL2, h_wgllwgll_yz);

  // textures
#ifdef USE_OPENCL
  if (run_opencl) {
#ifdef USE_TEXTURES_CONSTANTS
    cl_int errcode;
    cl_image_format format = {CL_R, CL_UNSIGNED_INT32};

    mp->d_hprime_xx_cm_tex = clCreateImage2D (mocl.context, CL_MEM_READ_ONLY, &format, NGLL2, 1, 0, mp->d_hprime_xx.ocl, clck_(&errcode));
    mp->d_hprimewgll_xx_cm_tex = clCreateImage2D (mocl.context, CL_MEM_READ_ONLY, &format, NGLL2, 1, 0, mp->d_hprimewgll_xx.ocl, clck_(&errcode));
#endif //USE_TEXTURES_CONSTANTS
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    // Using texture memory for the hprime-style constants is slower on
    // Fermi generation hardware, but *may* be faster on Kepler
    // generation. We will reevaluate this again, so might as well leave
    // in the code with #USE_TEXTURES_FIELDS not-defined.
#ifdef USE_TEXTURES_CONSTANTS
    // checks that realw is a float
    if (sizeof(realw) != sizeof(float)) exit_on_error("TEXTURES only work with realw selected as float");

    // note: device memory returned by cudaMalloc guarantee that the offset is 0,
    //       however here we use the global memory array d_hprime_xx and need to provide an offset variable
    size_t offset;

    // binds textures
#ifdef USE_OLDER_CUDA4_GPU
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<realw>();
    const textureReference* d_hprime_xx_tex_ptr;
    print_CUDA_error_if_any(cudaGetTextureReference(&d_hprime_xx_tex_ptr, "d_hprime_xx_tex"), 1101);
    print_CUDA_error_if_any(cudaBindTexture(&offset, d_hprime_xx_tex_ptr, mp->d_hprime_xx,
                                            &channelDesc, sizeof(realw)*(NGLL2)), 1102);
    // d_hprime_xx was cudaMalloced, so offset is 0
    //print_CUDA_error_if_any(cudaMemcpyToSymbol(d_hprime_xx_tex_offset,&offset,sizeof(offset)),11202);

    // weighted
    const textureReference* d_hprimewgll_xx_tex_ptr;
    print_CUDA_error_if_any(cudaGetTextureReference(&d_hprimewgll_xx_tex_ptr, "d_hprimewgll_xx_tex"), 1107);
    print_CUDA_error_if_any(cudaBindTexture(&offset, d_hprimewgll_xx_tex_ptr, mp->d_hprimewgll_xx,
                                            &channelDesc, sizeof(realw)*(NGLL2)), 1108);
    //print_CUDA_error_if_any(cudaMemcpyToSymbol(d_hprimewgll_xx_tex_offset,&offset,sizeof(offset)),11205);
#else
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<realw>();
    print_CUDA_error_if_any(cudaBindTexture(&offset, &d_hprime_xx_tex, mp->d_hprime_xx.cuda,
                                            &channelDesc, sizeof(realw)*(NGLL2)), 11201);
    // d_hprime_xx was cudaMalloced, so offset is 0
    //print_CUDA_error_if_any(cudaMemcpyToSymbol(d_hprime_xx_tex_offset,&offset,sizeof(offset)),11202);
    // debug
    //if (mp->myrank == 0 ) printf("texture constants hprime_xx: offset = %lu \n",offset);

    // weighted
    print_CUDA_error_if_any(cudaBindTexture(&offset, &d_hprimewgll_xx_tex, mp->d_hprimewgll_xx.cuda,
                                            &channelDesc, sizeof(realw)*(NGLL2)), 11204);
    //print_CUDA_error_if_any(cudaMemcpyToSymbol(d_hprimewgll_xx_tex_offset,&offset,sizeof(offset)),11205);
    // debug
    //if (mp->myrank == 0 ) printf("texture constants hprimewgll_xx: offset = %lu \n",offset);
#endif
#endif
  }
#endif

  // sets global parameters
  mp->NSPEC_CRUST_MANTLE = *NSPEC_CRUST_MANTLE;
  mp->NGLOB_CRUST_MANTLE = *NGLOB_CRUST_MANTLE;
  mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY = *NSPEC_CRUST_MANTLE_STRAIN_ONLY;

  mp->NSPEC_OUTER_CORE = *NSPEC_OUTER_CORE;
  mp->NGLOB_OUTER_CORE = *NGLOB_OUTER_CORE;

  mp->NSPEC_INNER_CORE = *NSPEC_INNER_CORE;
  mp->NGLOB_INNER_CORE = *NGLOB_INNER_CORE;
  mp->NSPEC_INNER_CORE_STRAIN_ONLY = *NSPEC_INNER_CORE_STRAIN_ONLY;

  // simulation type
  mp->simulation_type = *SIMULATION_TYPE;
  mp->noise_tomography = *NOISE_TOMOGRAPHY;

  // simulation flags initialization
  mp->save_forward = *SAVE_FORWARD_f;
  mp->absorbing_conditions = *ABSORBING_CONDITIONS_f;
  mp->oceans = *OCEANS_f;
  mp->gravity = *GRAVITY_f;
  mp->rotation = *ROTATION_f;
  mp->exact_mass_matrix_for_rotation = *EXACT_MASS_MATRIX_FOR_ROTATION_f;

  mp->attenuation = *ATTENUATION_f;
  mp->undo_attenuation = *UNDO_ATTENUATION_f;
  mp->partial_phys_dispersion_only = *PARTIAL_PHYS_DISPERSION_ONLY_f;
  mp->use_3d_attenuation_arrays = *USE_3D_ATTENUATION_ARRAYS_f;

  mp->compute_and_store_strain = *COMPUTE_AND_STORE_STRAIN_f;
  mp->anisotropic_3D_mantle = *ANISOTROPIC_3D_MANTLE_f;
  mp->anisotropic_inner_core = *ANISOTROPIC_INNER_CORE_f;
  mp->save_boundary_mesh = *SAVE_BOUNDARY_MESH_f;

  mp->anisotropic_kl = *ANISOTROPIC_KL_f;
  mp->approximate_hess_kl = *APPROXIMATE_HESS_KL_f;

  // mesh coloring flag
#ifdef USE_MESH_COLORING_GPU
  mp->use_mesh_coloring_gpu = 1;
#else
  // mesh coloring
  // note: this here passes the coloring as an option to the kernel routines
  //          the performance seems to be the same if one uses the pre-processing directives above or not
  mp->use_mesh_coloring_gpu = *USE_MESH_COLORING_GPU_f;
#endif

  // sources
  mp->nsources_local = *nsources_local;
  if (mp->simulation_type == 1 || mp->simulation_type == 3) {
    // not needed in case of pure adjoint simulations (SIMULATION_TYPE == 2)
    gpuCreateCopy_todevice_realw (&mp->d_sourcearrays, h_sourcearrays, (*NSOURCES) * NDIM * NGLL3);

    // allocates buffer on GPU for source time function values
    gpuMalloc_double (&mp->d_stf_pre_compute, *NSOURCES);
  }
  gpuCreateCopy_todevice_int (&mp->d_islice_selected_source, h_islice_selected_source, *NSOURCES);
  gpuCreateCopy_todevice_int (&mp->d_ispec_selected_source, h_ispec_selected_source, *NSOURCES);

  // receiver stations
  // note that:   size (number_receiver_global) = nrec_local
  //                   size (ispec_selected_rec) = nrec
  // number of receiver located in this partition

  mp->nrec_local = *nrec_local;
  if (mp->nrec_local > 0) {
    gpuCreateCopy_todevice_int (&mp->d_number_receiver_global, h_number_receiver_global, mp->nrec_local);
    // for seismograms
    gpuMalloc_realw (&mp->d_station_seismo_field, NDIM * NGLL3 * mp->nrec_local);

    // for transferring values from GPU to CPU
    if (GPU_ASYNC_COPY) {
#ifdef USE_OPENCL
      if (run_opencl) {
        ALLOC_PINNED_BUFFER_OCL (station_seismo_field, NDIM * NGLL3 * mp->nrec_local * sizeof(realw));
      }
#endif
#ifdef USE_CUDA
      if (run_cuda) {
        // TODO
        // only pinned memory can handle memcpy calls asynchronously
        print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_station_seismo_field), NDIM*NGLL3*(mp->nrec_local)*sizeof(realw)),4015);
      }
#endif
    } else {
      mp->h_station_seismo_field = (realw *) malloc (NDIM * NGLL3 * mp->nrec_local * sizeof(realw));
      if (mp->h_station_seismo_field == NULL) { exit_on_error("h_station_seismo_field not allocated \n"); }
    }

    // for adjoint strain
    if (mp->simulation_type == 2) {
      gpuMalloc_realw (&mp->d_station_strain_field, NGLL3 * mp->nrec_local);

      mp->h_station_strain_field = (realw *) malloc (NGLL3 * mp->nrec_local * sizeof(realw));
      if (mp->h_station_strain_field == NULL) { exit_on_error("h_station_strain_field not allocated \n"); }
    }
  }

  gpuCreateCopy_todevice_int (&mp->d_ispec_selected_rec, h_ispec_selected_rec, *nrec);

  // receiver adjoint source arrays only used for noise and adjoint simulations
  // adjoint source arrays

  mp->nadj_rec_local = *nadj_rec_local;
  if (mp->nadj_rec_local > 0) {
    // adjoint simulations
    // checks simulation flag
    if (mp->simulation_type != 2 && mp->simulation_type != 3 )
      exit_on_error ("prepare_constants_device: nadj_rec_local invalid with simulation type\n");

    // prepares local irec array:
    gpuMalloc_int (&mp->d_pre_computed_irec, mp->nadj_rec_local);

    // the irec_local variable needs to be precomputed (as
    // h_pre_comp..), because normally it is in the loop updating accel,
    // and due to how it's incremented, it cannot be parallelized
    int *h_pre_computed_irec = (int *) malloc (mp->nadj_rec_local * sizeof (int));
    if (h_pre_computed_irec == NULL) exit_on_error ("h_pre_computed_irec not allocated\n");

    // sets up local irec array
    int irec_local = 0;
    int irec;
    for (irec = 0; irec < *nrec; irec++) {
      if (mp->myrank == h_islice_selected_rec[irec]) {
        irec_local++;
        h_pre_computed_irec[irec_local-1] = irec;
      }
    }
    if (irec_local != mp->nadj_rec_local) exit_on_error ("prepare_constants_device: irec_local not equal to nadj_rec_local\n");

    // copies values onto GPU
    gpuCopy_todevice_int (&mp->d_pre_computed_irec, h_pre_computed_irec, mp->nadj_rec_local);

    free (h_pre_computed_irec);

    // temporary array to prepare extracted source array values
    if (GPU_ASYNC_COPY) {
#ifdef USE_OPENCL
      if (run_opencl) {
        ALLOC_PINNED_BUFFER_OCL(adj_sourcearrays_slice, mp->nadj_rec_local * NDIM * NGLL3 * sizeof(realw));
      }
#endif
#ifdef USE_CUDA
      if (run_cuda) {
        // note: Allocate pinned buffers otherwise cudaMemcpyAsync() will behave like cudaMemcpy(), i.e. synchronously.
        print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_adj_sourcearrays_slice),
                                               (mp->nadj_rec_local)*NDIM*NGLL3*sizeof(realw)),6011);
      }
#endif
    } else {
      mp->h_adj_sourcearrays_slice = (realw *) malloc (mp->nadj_rec_local * NDIM * NGLL3 * sizeof (realw));
      if (mp->h_adj_sourcearrays_slice == NULL) exit_on_error ("h_adj_sourcearrays_slice not allocated\n");
    }

    gpuMalloc_realw (&mp->d_adj_sourcearrays, mp->nadj_rec_local * NDIM * NGLL3);
  }

  // for rotation and new attenuation
  mp->deltat = *deltat_f;
  if (mp->simulation_type == 3) {
    mp->b_deltat = *b_deltat_f;
  }

  // initializes for rotational effects
  mp->two_omega_earth = 0.f;
  mp->b_two_omega_earth = 0.f;

#ifdef USE_OPENCL
  if (run_opencl) {
    mp->has_last_copy_evt = 0;
  }
#endif

  // buffer for norm checking at every timestamp
  int blocksize = BLOCKSIZE_TRANSFER;

  int size,size_padded;
  int num_blocks_x,num_blocks_y;

  // buffer for crust_mantle arrays has maximum size
  if (mp->compute_and_store_strain){
    size = MAX(mp->NGLOB_CRUST_MANTLE, NGLL3 * (mp->NSPEC_CRUST_MANTLE));
  } else{
    size = mp->NGLOB_CRUST_MANTLE;
  }
  size_padded = ((int) ceil (((double) size) / ((double) blocksize))) * blocksize;
  get_blocks_xy (size_padded / blocksize, &num_blocks_x, &num_blocks_y);

  // creates buffer on GPU for maximum array values
  gpuMalloc_realw (&mp->d_norm_max, num_blocks_x * num_blocks_y);

  // synchronizes gpu calls
  gpuSynchronize();

  GPU_ERROR_CHECKING ("prepare_constants_device");
}

/*----------------------------------------------------------------------------------------------- */
// ROTATION simulations
/*----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_ (prepare_fields_rotation_device,
               PREPARE_FIELDS_ROTATION_DEVICE) (long *Mesh_pointer_f,
                                                realw *two_omega_earth_f,
                                                realw *A_array_rotation,
                                                realw *B_array_rotation,
                                                realw *b_two_omega_earth_f,
                                                realw *b_A_array_rotation,
                                                realw *b_B_array_rotation,
                                                int *NSPEC_OUTER_CORE_ROTATION) {

  TRACE ("prepare_fields_rotation_device");
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // arrays only needed when rotation is required

  if (! mp->rotation) {exit_on_error("prepare_fields_rotation_device: rotation flag not properly initialized");}

  // checks array size
  if (*NSPEC_OUTER_CORE_ROTATION != mp->NSPEC_OUTER_CORE) {
    printf ("Error prepare_fields_rotation_device: rotation array has wrong size: %d instead of %d\n",
            *NSPEC_OUTER_CORE_ROTATION, mp->NSPEC_OUTER_CORE);
    exit_on_error ("prepare_fields_rotation_device: rotation array has wrong size");
  }

  // rotation arrays (needed only for outer core region)
  mp->two_omega_earth = *two_omega_earth_f;
  gpuCreateCopy_todevice_realw (&mp->d_A_array_rotation, A_array_rotation, NGLL3 * mp->NSPEC_OUTER_CORE);
  gpuCreateCopy_todevice_realw (&mp->d_B_array_rotation, B_array_rotation, NGLL3 * mp->NSPEC_OUTER_CORE);

  // backward/reconstructed fields
  if (mp->simulation_type == 3) {
    mp->b_two_omega_earth = *b_two_omega_earth_f;
    gpuCreateCopy_todevice_realw (&mp->d_b_A_array_rotation, b_A_array_rotation, NGLL3 * mp->NSPEC_OUTER_CORE);
    gpuCreateCopy_todevice_realw (&mp->d_b_B_array_rotation, b_B_array_rotation, NGLL3 * mp->NSPEC_OUTER_CORE);
  }

  GPU_ERROR_CHECKING ("prepare_fields_rotation_device");
}

/*----------------------------------------------------------------------------------------------- */
// GRAVITY simulations
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (prepare_fields_gravity_device,
               PREPARE_FIELDS_gravity_DEVICE) (long *Mesh_pointer_f,
                                               realw *d_ln_density_dr_table,
                                               realw *minus_rho_g_over_kappa_fluid,
                                               realw *minus_gravity_table,
                                               realw *minus_deriv_gravity_table,
                                               realw *density_table,
                                               realw *h_wgll_cube,
                                               int *NRAD_GRAVITY,
                                               realw *minus_g_icb,
                                               realw *minus_g_cmb,
                                               double *RHO_BOTTOM_OC,
                                               double *RHO_TOP_OC) {

  TRACE ("prepare_fields_gravity_device");
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  if (!mp->gravity) {
    // no gravity case
    // d ln (rho)/dr needed for the no gravity fluid potential
    gpuCreateCopy_todevice_realw (&mp->d_d_ln_density_dr_table, d_ln_density_dr_table, *NRAD_GRAVITY);

  } else {
    // gravity case
    mp->minus_g_icb = *minus_g_icb;
    mp->minus_g_cmb = *minus_g_cmb;

    // sets up GLL weights cubed
    gpuSetConst (&mp->d_wgll_cube, NGLL3, h_wgll_cube);

    // prepares gravity arrays
    gpuCreateCopy_todevice_realw (&mp->d_minus_rho_g_over_kappa_fluid, minus_rho_g_over_kappa_fluid, *NRAD_GRAVITY);
    gpuCreateCopy_todevice_realw (&mp->d_minus_gravity_table, minus_gravity_table, *NRAD_GRAVITY);
    gpuCreateCopy_todevice_realw (&mp->d_minus_deriv_gravity_table, minus_deriv_gravity_table, *NRAD_GRAVITY);
    gpuCreateCopy_todevice_realw (&mp->d_density_table, density_table, *NRAD_GRAVITY);
  }

  // constants
  mp->RHO_BOTTOM_OC = (realw) *RHO_BOTTOM_OC;
  mp->RHO_TOP_OC = (realw) *RHO_TOP_OC;

  GPU_ERROR_CHECKING ("prepare_fields_gravity_device");
}


/*----------------------------------------------------------------------------------------------- */
// ATTENUATION simulations
/*----------------------------------------------------------------------------------------------- */


extern EXTERN_LANG
void FC_FUNC_ (prepare_fields_attenuat_device,
               PREPARE_FIELDS_ATTENUAT_DEVICE) (long *Mesh_pointer_f,
                                                realw *R_xx_crust_mantle,
                                                realw *R_yy_crust_mantle,
                                                realw *R_xy_crust_mantle,
                                                realw *R_xz_crust_mantle,
                                                realw *R_yz_crust_mantle,
                                                realw *b_R_xx_crust_mantle,
                                                realw *b_R_yy_crust_mantle,
                                                realw *b_R_xy_crust_mantle,
                                                realw *b_R_xz_crust_mantle,
                                                realw *b_R_yz_crust_mantle,
                                                realw *factor_common_crust_mantle,
                                                realw *one_minus_sum_beta_crust_mantle,
                                                realw *R_xx_inner_core,
                                                realw *R_yy_inner_core,
                                                realw *R_xy_inner_core,
                                                realw *R_xz_inner_core,
                                                realw *R_yz_inner_core,
                                                realw *b_R_xx_inner_core,
                                                realw *b_R_yy_inner_core,
                                                realw *b_R_xy_inner_core,
                                                realw *b_R_xz_inner_core,
                                                realw *b_R_yz_inner_core,
                                                realw *factor_common_inner_core,
                                                realw *one_minus_sum_beta_inner_core,
                                                realw *alphaval, realw *betaval, realw *gammaval,
                                                realw *b_alphaval, realw *b_betaval, realw *b_gammaval) {

  TRACE ("prepare_fields_attenuat_device");
  int R_size1, R_size2, R_size3;
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // checks flag
  if (! mp->attenuation) { exit_on_error("prepare_fields_attenuat_device attenuation not properly initialized"); }

  // crust_mantle
  R_size1 = N_SLS*NGLL3*mp->NSPEC_CRUST_MANTLE;
  if (mp->use_3d_attenuation_arrays) {
    R_size2 = NGLL3*mp->NSPEC_CRUST_MANTLE;
    R_size3 = N_SLS*NGLL3*mp->NSPEC_CRUST_MANTLE;
  } else {
    R_size2 = 1*mp->NSPEC_CRUST_MANTLE;
    R_size3 = N_SLS*1*mp->NSPEC_CRUST_MANTLE;
  }

  gpuCreateCopy_todevice_realw (&mp->d_one_minus_sum_beta_crust_mantle, one_minus_sum_beta_crust_mantle, R_size2);

  if (! mp->partial_phys_dispersion_only) {
    // common factor
    gpuCreateCopy_todevice_realw (&mp->d_factor_common_crust_mantle, factor_common_crust_mantle, R_size3);

    // memory variables
    gpuCreateCopy_todevice_realw (&mp->d_R_xx_crust_mantle, R_xx_crust_mantle, R_size1);
    gpuCreateCopy_todevice_realw (&mp->d_R_yy_crust_mantle, R_yy_crust_mantle, R_size1);
    gpuCreateCopy_todevice_realw (&mp->d_R_xy_crust_mantle, R_xy_crust_mantle, R_size1);
    gpuCreateCopy_todevice_realw (&mp->d_R_xz_crust_mantle, R_xz_crust_mantle, R_size1);
    gpuCreateCopy_todevice_realw (&mp->d_R_yz_crust_mantle, R_yz_crust_mantle, R_size1);
  }

  if (mp->simulation_type == 3) {
    if (! mp->partial_phys_dispersion_only) {
      // memory variables
      gpuCreateCopy_todevice_realw (&mp->d_b_R_xx_crust_mantle, b_R_xx_crust_mantle, R_size1);
      gpuCreateCopy_todevice_realw (&mp->d_b_R_yy_crust_mantle, b_R_yy_crust_mantle, R_size1);
      gpuCreateCopy_todevice_realw (&mp->d_b_R_xy_crust_mantle, b_R_xy_crust_mantle, R_size1);
      gpuCreateCopy_todevice_realw (&mp->d_b_R_xz_crust_mantle, b_R_xz_crust_mantle, R_size1);
      gpuCreateCopy_todevice_realw (&mp->d_b_R_yz_crust_mantle, b_R_yz_crust_mantle, R_size1);
    }
  }

  // inner_core
  R_size1 = N_SLS*NGLL3*mp->NSPEC_INNER_CORE;
  if (mp->use_3d_attenuation_arrays) {
    R_size2 = NGLL3*mp->NSPEC_INNER_CORE;
    R_size3 = N_SLS*NGLL3*mp->NSPEC_INNER_CORE;
  } else {
    R_size2 = 1*mp->NSPEC_INNER_CORE;
    R_size3 = N_SLS*1*mp->NSPEC_INNER_CORE;
  }

  gpuCreateCopy_todevice_realw (&mp->d_one_minus_sum_beta_inner_core, one_minus_sum_beta_inner_core, R_size2);

  if (! mp->partial_phys_dispersion_only) {
    // common factor
    gpuCreateCopy_todevice_realw (&mp->d_factor_common_inner_core, factor_common_inner_core, R_size3);

    // memory variables
    gpuCreateCopy_todevice_realw (&mp->d_R_xx_inner_core, R_xx_inner_core, R_size1);
    gpuCreateCopy_todevice_realw (&mp->d_R_yy_inner_core, R_yy_inner_core, R_size1);
    gpuCreateCopy_todevice_realw (&mp->d_R_xy_inner_core, R_xy_inner_core, R_size1);
    gpuCreateCopy_todevice_realw (&mp->d_R_xz_inner_core, R_xz_inner_core, R_size1);
    gpuCreateCopy_todevice_realw (&mp->d_R_yz_inner_core, R_yz_inner_core, R_size1);
  }

  if (mp->simulation_type == 3) {
    if (! mp->partial_phys_dispersion_only) {
      // memory variables
      gpuCreateCopy_todevice_realw (&mp->d_b_R_xx_inner_core, b_R_xx_inner_core, R_size1);
      gpuCreateCopy_todevice_realw (&mp->d_b_R_yy_inner_core, b_R_yy_inner_core, R_size1);
      gpuCreateCopy_todevice_realw (&mp->d_b_R_xy_inner_core, b_R_xy_inner_core, R_size1);
      gpuCreateCopy_todevice_realw (&mp->d_b_R_xz_inner_core, b_R_xz_inner_core, R_size1);
      gpuCreateCopy_todevice_realw (&mp->d_b_R_yz_inner_core, b_R_yz_inner_core, R_size1);
    }
  }

  // alpha, beta, gamma factors
  gpuCreateCopy_todevice_realw (&mp->d_alphaval, alphaval, N_SLS);
  gpuCreateCopy_todevice_realw (&mp->d_betaval, betaval, N_SLS);
  gpuCreateCopy_todevice_realw (&mp->d_gammaval, gammaval, N_SLS);

  if (mp->simulation_type == 3) {
    // alpha, beta, gamma factors for backward fields
    gpuCreateCopy_todevice_realw (&mp->d_b_alphaval, b_alphaval, N_SLS);
    gpuCreateCopy_todevice_realw (&mp->d_b_betaval, b_betaval, N_SLS);
    gpuCreateCopy_todevice_realw (&mp->d_b_gammaval, b_gammaval, N_SLS);
  }

  GPU_ERROR_CHECKING ("prepare_fields_attenuat_device");
}

/*----------------------------------------------------------------------------------------------- */
// STRAIN simulations
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (prepare_fields_strain_device,
               PREPARE_FIELDS_STRAIN_DEVICE) (long *Mesh_pointer_f,
                                              realw *epsilondev_xx_crust_mantle,
                                              realw *epsilondev_yy_crust_mantle,
                                              realw *epsilondev_xy_crust_mantle,
                                              realw *epsilondev_xz_crust_mantle,
                                              realw *epsilondev_yz_crust_mantle,
                                              realw *b_epsilondev_xx_crust_mantle,
                                              realw *b_epsilondev_yy_crust_mantle,
                                              realw *b_epsilondev_xy_crust_mantle,
                                              realw *b_epsilondev_xz_crust_mantle,
                                              realw *b_epsilondev_yz_crust_mantle,
                                              realw *eps_trace_over_3_crust_mantle,
                                              realw *b_eps_trace_over_3_crust_mantle,
                                              realw *epsilondev_xx_inner_core,
                                              realw *epsilondev_yy_inner_core,
                                              realw *epsilondev_xy_inner_core,
                                              realw *epsilondev_xz_inner_core,
                                              realw *epsilondev_yz_inner_core,
                                              realw *b_epsilondev_xx_inner_core,
                                              realw *b_epsilondev_yy_inner_core,
                                              realw *b_epsilondev_xy_inner_core,
                                              realw *b_epsilondev_xz_inner_core,
                                              realw *b_epsilondev_yz_inner_core,
                                              realw *eps_trace_over_3_inner_core,
                                              realw *b_eps_trace_over_3_inner_core) {

  TRACE ("prepare_fields_strain_device");
  int R_size, size_strain_only;
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // checks flag
  if (!mp->compute_and_store_strain) {
    exit_on_error ("prepare_fields_strain_device strain not properly initialized");
  }

  // crust_mantle
  R_size = NGLL3*mp->NSPEC_CRUST_MANTLE;
  gpuCreateCopy_todevice_realw (&mp->d_epsilondev_xx_crust_mantle, epsilondev_xx_crust_mantle, R_size);
  gpuCreateCopy_todevice_realw (&mp->d_epsilondev_yy_crust_mantle, epsilondev_yy_crust_mantle, R_size);
  gpuCreateCopy_todevice_realw (&mp->d_epsilondev_xy_crust_mantle, epsilondev_xy_crust_mantle, R_size);
  gpuCreateCopy_todevice_realw (&mp->d_epsilondev_xz_crust_mantle, epsilondev_xz_crust_mantle, R_size);
  gpuCreateCopy_todevice_realw (&mp->d_epsilondev_yz_crust_mantle, epsilondev_yz_crust_mantle, R_size);

  // strain
  size_strain_only = NGLL3 * mp->NSPEC_CRUST_MANTLE_STRAIN_ONLY;
  gpuCreateCopy_todevice_realw (&mp->d_eps_trace_over_3_crust_mantle, eps_trace_over_3_crust_mantle, size_strain_only);

  // backward/reconstructed fields
  if (mp->simulation_type == 3) {
    if (mp->undo_attenuation) {
      // strain will be computed locally based on displacement wavefield
      // only uses pointers to already allocated arrays
      mp->d_b_epsilondev_xx_crust_mantle = gpuTakeRef(mp->d_epsilondev_xx_crust_mantle);
      mp->d_b_epsilondev_yy_crust_mantle = gpuTakeRef(mp->d_epsilondev_yy_crust_mantle);
      mp->d_b_epsilondev_xy_crust_mantle = gpuTakeRef(mp->d_epsilondev_xy_crust_mantle);
      mp->d_b_epsilondev_xz_crust_mantle = gpuTakeRef(mp->d_epsilondev_xz_crust_mantle);
      mp->d_b_epsilondev_yz_crust_mantle = gpuTakeRef(mp->d_epsilondev_yz_crust_mantle);
      mp->d_b_eps_trace_over_3_crust_mantle = gpuTakeRef(mp->d_eps_trace_over_3_crust_mantle);
    } else {
      gpuCreateCopy_todevice_realw (&mp->d_b_epsilondev_xx_crust_mantle, b_epsilondev_xx_crust_mantle, R_size);
      gpuCreateCopy_todevice_realw (&mp->d_b_epsilondev_yy_crust_mantle, b_epsilondev_yy_crust_mantle, R_size);
      gpuCreateCopy_todevice_realw (&mp->d_b_epsilondev_xy_crust_mantle, b_epsilondev_xy_crust_mantle, R_size);
      gpuCreateCopy_todevice_realw (&mp->d_b_epsilondev_xz_crust_mantle, b_epsilondev_xz_crust_mantle, R_size);
      gpuCreateCopy_todevice_realw (&mp->d_b_epsilondev_yz_crust_mantle, b_epsilondev_yz_crust_mantle, R_size);
      //strain
      gpuCreateCopy_todevice_realw (&mp->d_b_eps_trace_over_3_crust_mantle, b_eps_trace_over_3_crust_mantle, R_size);
    }
  }

  // inner_core
  R_size = NGLL3*mp->NSPEC_INNER_CORE;
  gpuCreateCopy_todevice_realw (&mp->d_epsilondev_xx_inner_core, epsilondev_xx_inner_core, R_size);
  gpuCreateCopy_todevice_realw (&mp->d_epsilondev_yy_inner_core, epsilondev_yy_inner_core, R_size);
  gpuCreateCopy_todevice_realw (&mp->d_epsilondev_xy_inner_core, epsilondev_xy_inner_core, R_size);
  gpuCreateCopy_todevice_realw (&mp->d_epsilondev_xz_inner_core, epsilondev_xz_inner_core, R_size);
  gpuCreateCopy_todevice_realw (&mp->d_epsilondev_yz_inner_core, epsilondev_yz_inner_core, R_size);

  // strain
  size_strain_only = NGLL3 * (mp->NSPEC_INNER_CORE_STRAIN_ONLY);
  gpuCreateCopy_todevice_realw (&mp->d_eps_trace_over_3_inner_core, eps_trace_over_3_inner_core, size_strain_only);

  // backward/reconstructed fields
  if (mp->simulation_type == 3) {
    if (mp->undo_attenuation) {
      // strain will be computed locally based on displacement wavefield
      // only uses pointers to already allocated arrays
      mp->d_b_epsilondev_xx_inner_core = gpuTakeRef(mp->d_epsilondev_xx_inner_core);
      mp->d_b_epsilondev_yy_inner_core = gpuTakeRef(mp->d_epsilondev_yy_inner_core);
      mp->d_b_epsilondev_xy_inner_core = gpuTakeRef(mp->d_epsilondev_xy_inner_core);
      mp->d_b_epsilondev_xz_inner_core = gpuTakeRef(mp->d_epsilondev_xz_inner_core);
      mp->d_b_epsilondev_yz_inner_core = gpuTakeRef(mp->d_epsilondev_yz_inner_core);
      mp->d_b_eps_trace_over_3_inner_core = gpuTakeRef(mp->d_eps_trace_over_3_inner_core);
    } else {
      gpuCreateCopy_todevice_realw (&mp->d_b_epsilondev_xx_inner_core, b_epsilondev_xx_inner_core, R_size);
      gpuCreateCopy_todevice_realw (&mp->d_b_epsilondev_yy_inner_core, b_epsilondev_yy_inner_core, R_size);
      gpuCreateCopy_todevice_realw (&mp->d_b_epsilondev_xy_inner_core, b_epsilondev_xy_inner_core, R_size);
      gpuCreateCopy_todevice_realw (&mp->d_b_epsilondev_xz_inner_core, b_epsilondev_xz_inner_core, R_size);
      gpuCreateCopy_todevice_realw (&mp->d_b_epsilondev_yz_inner_core, b_epsilondev_yz_inner_core, R_size);
      // strain
      gpuCreateCopy_todevice_realw (&mp->d_b_eps_trace_over_3_inner_core, b_eps_trace_over_3_inner_core, R_size);
    }
  }

  GPU_ERROR_CHECKING ("prepare_fields_strain_device");
}

/*----------------------------------------------------------------------------------------------- */
// STRAIN simulations
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (prepare_fields_absorb_device,
               PREPARE_FIELDS_ABSORB_DEVICE) (long *Mesh_pointer_f,
                                              int *nspec2D_xmin_crust_mantle, int *nspec2D_xmax_crust_mantle,
                                              int *nspec2D_ymin_crust_mantle, int *nspec2D_ymax_crust_mantle,
                                              int *NSPEC2DMAX_XMIN_XMAX_CM, int *NSPEC2DMAX_YMIN_YMAX_CM,
                                              int *nimin_crust_mantle, int *nimax_crust_mantle,
                                              int *njmin_crust_mantle, int *njmax_crust_mantle,
                                              int *nkmin_xi_crust_mantle, int *nkmin_eta_crust_mantle,
                                              int *ibelm_xmin_crust_mantle, int *ibelm_xmax_crust_mantle,
                                              int *ibelm_ymin_crust_mantle, int *ibelm_ymax_crust_mantle,
                                              realw *normal_xmin_crust_mantle, realw *normal_xmax_crust_mantle,
                                              realw *normal_ymin_crust_mantle, realw *normal_ymax_crust_mantle,
                                              realw *jacobian2D_xmin_crust_mantle, realw *jacobian2D_xmax_crust_mantle,
                                              realw *jacobian2D_ymin_crust_mantle, realw *jacobian2D_ymax_crust_mantle,
                                              realw *rho_vp_crust_mantle,
                                              realw *rho_vs_crust_mantle,
                                              int *nspec2D_xmin_outer_core, int *nspec2D_xmax_outer_core,
                                              int *nspec2D_ymin_outer_core, int *nspec2D_ymax_outer_core,
                                              int *nspec2D_zmin_outer_core,
                                              int *NSPEC2DMAX_XMIN_XMAX_OC, int *NSPEC2DMAX_YMIN_YMAX_OC,
                                              int *nimin_outer_core, int *nimax_outer_core,
                                              int *njmin_outer_core, int *njmax_outer_core,
                                              int *nkmin_xi_outer_core, int *nkmin_eta_outer_core,
                                              int *ibelm_xmin_outer_core, int *ibelm_xmax_outer_core,
                                              int *ibelm_ymin_outer_core, int *ibelm_ymax_outer_core,
                                              realw *jacobian2D_xmin_outer_core, realw *jacobian2D_xmax_outer_core,
                                              realw *jacobian2D_ymin_outer_core, realw *jacobian2D_ymax_outer_core,
                                              realw *vp_outer_core) {

  TRACE ("prepare_fields_absorb_device");
  int size;
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // checks flag
  if (! mp->absorbing_conditions) {
    exit_on_error ("prepare_fields_absorb_device absorbing_conditions not properly initialized");
  }

  // crust_mantle
  mp->nspec2D_xmin_crust_mantle = *nspec2D_xmin_crust_mantle;
  mp->nspec2D_xmax_crust_mantle = *nspec2D_xmax_crust_mantle;
  mp->nspec2D_ymin_crust_mantle = *nspec2D_ymin_crust_mantle;
  mp->nspec2D_ymax_crust_mantle = *nspec2D_ymax_crust_mantle;

  // vp & vs
  size = NGLL3 * (mp->NSPEC_CRUST_MANTLE);
  gpuCreateCopy_todevice_realw (&mp->d_rho_vp_crust_mantle, rho_vp_crust_mantle, size);
  gpuCreateCopy_todevice_realw (&mp->d_rho_vs_crust_mantle, rho_vs_crust_mantle, size);

  // ijk index arrays
  gpuCreateCopy_todevice_int (&mp->d_nkmin_xi_crust_mantle, nkmin_xi_crust_mantle, 2 * (*NSPEC2DMAX_XMIN_XMAX_CM));
  gpuCreateCopy_todevice_int (&mp->d_nkmin_eta_crust_mantle, nkmin_eta_crust_mantle, 2 * (*NSPEC2DMAX_YMIN_YMAX_CM));
  gpuCreateCopy_todevice_int (&mp->d_njmin_crust_mantle, njmin_crust_mantle, 2 * (*NSPEC2DMAX_XMIN_XMAX_CM));
  gpuCreateCopy_todevice_int (&mp->d_njmax_crust_mantle, njmax_crust_mantle, 2 * (*NSPEC2DMAX_XMIN_XMAX_CM));
  gpuCreateCopy_todevice_int (&mp->d_nimin_crust_mantle, nimin_crust_mantle, 2 * (*NSPEC2DMAX_YMIN_YMAX_CM));
  gpuCreateCopy_todevice_int (&mp->d_nimax_crust_mantle, nimax_crust_mantle, 2 * (*NSPEC2DMAX_YMIN_YMAX_CM));

  // xmin
  if (mp->nspec2D_xmin_crust_mantle > 0) {
    gpuCreateCopy_todevice_int (&mp->d_ibelm_xmin_crust_mantle, ibelm_xmin_crust_mantle, mp->nspec2D_xmin_crust_mantle);
    gpuCreateCopy_todevice_realw (&mp->d_normal_xmin_crust_mantle, normal_xmin_crust_mantle, NDIM*NGLL2 * mp->nspec2D_xmin_crust_mantle);
    gpuCreateCopy_todevice_realw (&mp->d_jacobian2D_xmin_crust_mantle, jacobian2D_xmin_crust_mantle, NGLL2 * mp->nspec2D_xmin_crust_mantle);

    // boundary buffer
    if ((mp->simulation_type == 1 && mp->save_forward) || (mp->simulation_type == 3)) {
      gpuMalloc_realw (&mp->d_absorb_xmin_crust_mantle, NDIM * NGLL2 * mp->nspec2D_xmin_crust_mantle);
    }
  }

  // xmax
  if (mp->nspec2D_xmax_crust_mantle > 0) {
    gpuCreateCopy_todevice_int (&mp->d_ibelm_xmax_crust_mantle, ibelm_xmax_crust_mantle, mp->nspec2D_xmax_crust_mantle);
    gpuCreateCopy_todevice_realw (&mp->d_normal_xmax_crust_mantle, normal_xmax_crust_mantle, NDIM*NGLL2 * mp->nspec2D_xmax_crust_mantle);
    gpuCreateCopy_todevice_realw (&mp->d_jacobian2D_xmax_crust_mantle, jacobian2D_xmax_crust_mantle, NGLL2 * mp->nspec2D_xmax_crust_mantle);
    // boundary buffer
    if ((mp->simulation_type == 1 && mp->save_forward) || (mp->simulation_type == 3)) {
      gpuMalloc_realw (&mp->d_absorb_xmax_crust_mantle, NDIM * NGLL2 * mp->nspec2D_xmax_crust_mantle);
    }
  }

  // ymin
  if (mp->nspec2D_ymin_crust_mantle > 0) {
    gpuCreateCopy_todevice_int (&mp->d_ibelm_ymin_crust_mantle, ibelm_ymin_crust_mantle, mp->nspec2D_ymin_crust_mantle);
    gpuCreateCopy_todevice_realw (&mp->d_normal_ymin_crust_mantle, normal_ymin_crust_mantle, NDIM*NGLL2 * mp->nspec2D_ymin_crust_mantle);
    gpuCreateCopy_todevice_realw (&mp->d_jacobian2D_ymin_crust_mantle, jacobian2D_ymin_crust_mantle, NGLL2 * mp->nspec2D_ymin_crust_mantle);
    // boundary buffer
    if ((mp->simulation_type == 1 && mp->save_forward) || (mp->simulation_type == 3)) {
      gpuMalloc_realw (&mp->d_absorb_ymin_crust_mantle, NDIM * NGLL2 * mp->nspec2D_ymin_crust_mantle);
    }
  }

  // ymax
  if (mp->nspec2D_ymax_crust_mantle > 0) {
    gpuCreateCopy_todevice_int (&mp->d_ibelm_ymax_crust_mantle,
                       ibelm_ymax_crust_mantle, mp->nspec2D_ymax_crust_mantle);
    gpuCreateCopy_todevice_realw (&mp->d_normal_ymax_crust_mantle, normal_ymax_crust_mantle, NDIM * NGLL2 * mp->nspec2D_ymax_crust_mantle);
    gpuCreateCopy_todevice_realw (&mp->d_jacobian2D_ymax_crust_mantle, jacobian2D_ymax_crust_mantle, NGLL2 * mp->nspec2D_ymax_crust_mantle);

    // boundary buffer
    if ((mp->simulation_type == 1 && mp->save_forward) || (mp->simulation_type == 3)) {
      gpuMalloc_realw (&mp->d_absorb_ymax_crust_mantle, NDIM * NGLL2 * mp->nspec2D_ymax_crust_mantle);
    }
  }

  // outer_core
  mp->nspec2D_xmin_outer_core = *nspec2D_xmin_outer_core;
  mp->nspec2D_xmax_outer_core = *nspec2D_xmax_outer_core;
  mp->nspec2D_ymin_outer_core = *nspec2D_ymin_outer_core;
  mp->nspec2D_ymax_outer_core = *nspec2D_ymax_outer_core;
  mp->nspec2D_zmin_outer_core = *nspec2D_zmin_outer_core;

  // vp
  size = NGLL3 * (mp->NSPEC_OUTER_CORE);
  gpuCreateCopy_todevice_realw (&mp->d_vp_outer_core, vp_outer_core, size);

  // ijk index arrays
  gpuCreateCopy_todevice_int (&mp->d_nkmin_xi_outer_core, nkmin_xi_outer_core, 2 * (*NSPEC2DMAX_XMIN_XMAX_OC));
  gpuCreateCopy_todevice_int (&mp->d_nkmin_eta_outer_core, nkmin_eta_outer_core, 2 * (*NSPEC2DMAX_YMIN_YMAX_OC));
  gpuCreateCopy_todevice_int (&mp->d_njmin_outer_core, njmin_outer_core, 2 * (*NSPEC2DMAX_XMIN_XMAX_OC));
  gpuCreateCopy_todevice_int (&mp->d_njmax_outer_core, njmax_outer_core, 2 * (*NSPEC2DMAX_XMIN_XMAX_OC));
  gpuCreateCopy_todevice_int (&mp->d_nimin_outer_core, nimin_outer_core, 2 * (*NSPEC2DMAX_YMIN_YMAX_OC));
  gpuCreateCopy_todevice_int (&mp->d_nimax_outer_core, nimax_outer_core, 2 * (*NSPEC2DMAX_YMIN_YMAX_OC));

  // xmin
  if (mp->nspec2D_xmin_outer_core > 0) {
    gpuCreateCopy_todevice_int (&mp->d_ibelm_xmin_outer_core, ibelm_xmin_outer_core, mp->nspec2D_xmin_outer_core);
    gpuCreateCopy_todevice_realw (&mp->d_jacobian2D_xmin_outer_core, jacobian2D_xmin_outer_core, NGLL2 * mp->nspec2D_xmin_outer_core);
    // boundary buffer
    if ((mp->simulation_type == 1 && mp->save_forward) || (mp->simulation_type == 3)) {
      gpuMalloc_realw (&mp->d_absorb_xmin_outer_core, NGLL2 * mp->nspec2D_xmin_outer_core);
    }
  }

  // xmax
  if (mp->nspec2D_xmax_outer_core > 0) {
    gpuCreateCopy_todevice_int (&mp->d_ibelm_xmax_outer_core, ibelm_xmax_outer_core, mp->nspec2D_xmax_outer_core);
    gpuCreateCopy_todevice_realw (&mp->d_jacobian2D_xmax_outer_core, jacobian2D_xmax_outer_core, NGLL2 * mp->nspec2D_xmax_outer_core);

    // boundary buffer
    if ((mp->simulation_type == 1 && mp->save_forward) || (mp->simulation_type == 3)) {
      gpuMalloc_realw (&mp->d_absorb_xmax_outer_core, NGLL2 * mp->nspec2D_xmax_outer_core);
    }
  }

  // ymin
  if (mp->nspec2D_ymin_outer_core > 0) {
    gpuCreateCopy_todevice_int (&mp->d_ibelm_ymin_outer_core, ibelm_ymin_outer_core, mp->nspec2D_ymin_outer_core);
    gpuCreateCopy_todevice_realw (&mp->d_jacobian2D_ymin_outer_core, jacobian2D_ymin_outer_core, NGLL2 * mp->nspec2D_ymin_outer_core);
    // boundary buffer
    if ((mp->simulation_type == 1 && mp->save_forward) || (mp->simulation_type == 3)) {
      gpuMalloc_realw (&mp->d_absorb_ymin_outer_core, NGLL2 * mp->nspec2D_ymin_outer_core);
    }
  }

  // ymax
  if (mp->nspec2D_ymax_outer_core > 0) {
    gpuCreateCopy_todevice_int (&mp->d_ibelm_ymax_outer_core, ibelm_ymax_outer_core, mp->nspec2D_ymax_outer_core);
    gpuCreateCopy_todevice_realw (&mp->d_jacobian2D_ymax_outer_core, jacobian2D_ymax_outer_core, NGLL2 * mp->nspec2D_ymax_outer_core);
    // boundary buffer
    if ((mp->simulation_type == 1 && mp->save_forward) || (mp->simulation_type == 3)) {
      gpuMalloc_realw (&mp->d_absorb_ymax_outer_core, NGLL2 * mp->nspec2D_ymax_outer_core);
    }
  }

  // zmin
  if (mp->nspec2D_zmin_outer_core > 0) {
    // note: ibelm_bottom_outer_core and jacobian2D_bottom_outer_core will be allocated
    //          when preparing the outer core
    // boundary buffer
    if ((mp->simulation_type == 1 && mp->save_forward) || (mp->simulation_type == 3)) {
      gpuMalloc_realw (&mp->d_absorb_zmin_outer_core, NGLL2 * mp->nspec2D_zmin_outer_core);
    }
  }

  GPU_ERROR_CHECKING ("prepare_fields_absorb_device");
}

/*----------------------------------------------------------------------------------------------- */
// MPI interfaces
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (prepare_mpi_buffers_device,
               PREPARE_MPI_BUFFERS_DEVICE) (long *Mesh_pointer_f,
                                            int *num_interfaces_crust_mantle,
                                            int *max_nibool_interfaces_cm,
                                            int *nibool_interfaces_crust_mantle,
                                            int *ibool_interfaces_crust_mantle,
                                            int *num_interfaces_inner_core,
                                            int *max_nibool_interfaces_ic,
                                            int *nibool_interfaces_inner_core,
                                            int *ibool_interfaces_inner_core,
                                            int *num_interfaces_outer_core,
                                            int *max_nibool_interfaces_oc,
                                            int *nibool_interfaces_outer_core,
                                            int *ibool_interfaces_outer_core) {

  TRACE ("prepare_mpi_buffers_device");

  Mesh *mp = (Mesh *) *Mesh_pointer_f;
  int size_mpi_buffer;

  // prepares interprocess-edge exchange information

  // crust/mantle mesh
  mp->num_interfaces_crust_mantle = *num_interfaces_crust_mantle;
  mp->max_nibool_interfaces_cm = *max_nibool_interfaces_cm;

  if (mp->num_interfaces_crust_mantle > 0) {
    // number of ibool entries array
    gpuCreateCopy_todevice_int (&mp->d_nibool_interfaces_crust_mantle, nibool_interfaces_crust_mantle,
                                mp->num_interfaces_crust_mantle);

    // ibool entries (iglob indices) values on interface
    gpuCreateCopy_todevice_int (&mp->d_ibool_interfaces_crust_mantle, ibool_interfaces_crust_mantle,
                                mp->num_interfaces_crust_mantle * mp->max_nibool_interfaces_cm);

    size_mpi_buffer = NDIM*(mp->max_nibool_interfaces_cm)*(mp->num_interfaces_crust_mantle);

    // allocates MPI buffer for exchange with CPU
    gpuMalloc_realw (&mp->d_send_accel_buffer_crust_mantle, size_mpi_buffer);
    if (mp->simulation_type == 3) {
      gpuMalloc_realw (&mp->d_b_send_accel_buffer_crust_mantle, size_mpi_buffer);
    }

    // asynchronous MPI buffer
    if (GPU_ASYNC_COPY) {
#ifdef USE_OPENCL
      if (run_opencl) {
        ALLOC_PINNED_BUFFER_OCL(send_accel_buffer_cm, sizeof(realw)* size_mpi_buffer);
        ALLOC_PINNED_BUFFER_OCL(recv_accel_buffer_cm, sizeof(realw)* size_mpi_buffer);

        if (mp->simulation_type == 3) {
          ALLOC_PINNED_BUFFER_OCL(b_send_accel_buffer_cm, sizeof(realw)* size_mpi_buffer);
          ALLOC_PINNED_BUFFER_OCL(b_recv_accel_buffer_cm, sizeof(realw)* size_mpi_buffer);
        }
      }
#endif
#ifdef USE_CUDA
      if (run_cuda) {
        // note: Allocate pinned MPI buffers.
        //       MPI buffers use pinned memory allocated by cudaMallocHost, which
        //       enables the use of asynchronous memory copies from host <-> device
        // send buffer
        print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_send_accel_buffer_cm),sizeof(realw)* size_mpi_buffer ),8004);
        // receive buffer
        print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_recv_accel_buffer_cm),sizeof(realw)* size_mpi_buffer ),8004);
        if (mp->simulation_type == 3) {
          print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_b_send_accel_buffer_cm),sizeof(realw)* size_mpi_buffer ),8004);
          print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_b_recv_accel_buffer_cm),sizeof(realw)* size_mpi_buffer ),8004);
        }
      }
#endif
    }
  }

  // inner core mesh
  mp->num_interfaces_inner_core = *num_interfaces_inner_core;
  mp->max_nibool_interfaces_ic = *max_nibool_interfaces_ic;
  if (mp->num_interfaces_inner_core > 0) {
    // number of ibool entries array
    gpuCreateCopy_todevice_int (&mp->d_nibool_interfaces_inner_core, nibool_interfaces_inner_core,
                                mp->num_interfaces_inner_core);

    // ibool entries (iglob indices) values on interface
    gpuCreateCopy_todevice_int (&mp->d_ibool_interfaces_inner_core, ibool_interfaces_inner_core,
                                mp->num_interfaces_inner_core * mp->max_nibool_interfaces_ic);

    size_mpi_buffer = NDIM * (mp->max_nibool_interfaces_ic) * (mp->num_interfaces_inner_core);

    gpuMalloc_realw (&mp->d_send_accel_buffer_inner_core, size_mpi_buffer);
    if (mp->simulation_type == 3) {
      gpuMalloc_realw (&mp->d_b_send_accel_buffer_inner_core, size_mpi_buffer);
    }

    // asynchronous MPI buffer
    if (GPU_ASYNC_COPY) {
#ifdef USE_OPENCL
      if (run_opencl) {
        ALLOC_PINNED_BUFFER_OCL(send_accel_buffer_ic, sizeof(realw)* size_mpi_buffer);
        ALLOC_PINNED_BUFFER_OCL(recv_accel_buffer_ic, sizeof(realw)* size_mpi_buffer);

        if (mp->simulation_type == 3) {
          ALLOC_PINNED_BUFFER_OCL(b_send_accel_buffer_ic, sizeof(realw)* size_mpi_buffer);
          ALLOC_PINNED_BUFFER_OCL(b_recv_accel_buffer_ic, sizeof(realw)* size_mpi_buffer);
        }
      }
#endif
#ifdef USE_CUDA
      if (run_cuda) {
        // note: Allocate pinned MPI buffers.
        //       MPI buffers use pinned memory allocated by cudaMallocHost, which
        //       enables the use of asynchronous memory copies from host <-> device
        // send buffer
        print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_send_accel_buffer_ic),sizeof(realw)*size_mpi_buffer ),8004);
        // receive buffer
        print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_recv_accel_buffer_ic),sizeof(realw)*size_mpi_buffer ),8004);
        // adjoint
        if (mp->simulation_type == 3) {
          print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_b_send_accel_buffer_ic),sizeof(realw)*size_mpi_buffer ),8004);
          print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_b_recv_accel_buffer_ic),sizeof(realw)*size_mpi_buffer ),8004);
        }
      }
#endif
    }
  }
  // outer core mesh
  // note: uses only scalar wavefield arrays
  mp->num_interfaces_outer_core = *num_interfaces_outer_core;
  mp->max_nibool_interfaces_oc = *max_nibool_interfaces_oc;
  if (mp->num_interfaces_outer_core > 0) {
    // number of ibool entries array
    gpuCreateCopy_todevice_int (&mp->d_nibool_interfaces_outer_core, nibool_interfaces_outer_core,
                                mp->num_interfaces_outer_core);

    // ibool entries (iglob indices) values on interface
    gpuCreateCopy_todevice_int (&mp->d_ibool_interfaces_outer_core, ibool_interfaces_outer_core,
                                (mp->num_interfaces_outer_core) * (mp->max_nibool_interfaces_oc));

    size_mpi_buffer = (mp->max_nibool_interfaces_oc)*(mp->num_interfaces_outer_core);

    // allocates MPI buffer for exchange with CPU
    gpuMalloc_realw (&mp->d_send_accel_buffer_outer_core, size_mpi_buffer);
    if (mp->simulation_type == 3) {
      gpuMalloc_realw (&mp->d_b_send_accel_buffer_outer_core, size_mpi_buffer);
    }

    // asynchronous MPI buffer
    if (GPU_ASYNC_COPY) {
#ifdef USE_OPENCL
      if (run_opencl) {
        ALLOC_PINNED_BUFFER_OCL(send_accel_buffer_oc, sizeof(realw)* size_mpi_buffer);
        ALLOC_PINNED_BUFFER_OCL(recv_accel_buffer_oc, sizeof(realw)* size_mpi_buffer);

        if (mp->simulation_type == 3) {
          ALLOC_PINNED_BUFFER_OCL(b_send_accel_buffer_oc, sizeof(realw)* size_mpi_buffer);
          ALLOC_PINNED_BUFFER_OCL(b_recv_accel_buffer_oc, sizeof(realw)* size_mpi_buffer);
        }
      }
#endif
#ifdef USE_CUDA
      if (run_cuda) {
        // note: Allocate pinned MPI buffers.
        //       MPI buffers use pinned memory allocated by cudaMallocHost, which
        //       enables the use of asynchronous memory copies from host <-> device
        // send buffer
        print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_send_accel_buffer_oc),sizeof(realw)*size_mpi_buffer ),8004);
        // receive buffer
        print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_recv_accel_buffer_oc),sizeof(realw)*size_mpi_buffer ),8004);
        if (mp->simulation_type == 3) {
          print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_b_send_accel_buffer_oc),sizeof(realw)*size_mpi_buffer ),8004);
          print_CUDA_error_if_any(cudaMallocHost((void**)&(mp->h_b_recv_accel_buffer_oc),sizeof(realw)*size_mpi_buffer ),8004);
        }
      }
#endif
    }
  }

  GPU_ERROR_CHECKING ("prepare_mpi_buffers_device");
}


/*----------------------------------------------------------------------------------------------- */
// for NOISE simulations
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (prepare_fields_noise_device,
               PREPARE_FIELDS_NOISE_DEVICE) (long *Mesh_pointer_f,
                                             int *NSPEC_TOP,
                                             int *NSTEP,
                                             int *h_ibelm_top_crust_mantle,
                                             realw *noise_sourcearray,
                                             realw *normal_x_noise,
                                             realw *normal_y_noise,
                                             realw *normal_z_noise,
                                             realw *mask_noise,
                                             realw *jacobian2D_top_crust_mantle) {

  TRACE ("prepare_fields_noise_device");

  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // free surface
  mp->nspec2D_top_crust_mantle = *NSPEC_TOP;
  if (mp->nspec2D_top_crust_mantle > 0) {
    // note: d_ibelm_top_crust_mantle will only be needed for noise computations
    gpuCreateCopy_todevice_int (&mp->d_ibelm_top_crust_mantle, h_ibelm_top_crust_mantle, mp->nspec2D_top_crust_mantle);

    // alloc storage for the surface buffer to be copied
    gpuMalloc_realw (&mp->d_noise_surface_movie, NDIM * NGLL2 * mp->nspec2D_top_crust_mantle);

  } else {
    // for global mesh: each crust/mantle slice should have at top a free surface
    exit_on_error ("prepare_fields_noise_device NSPEC_TOP not properly initialized");
  }

  // prepares noise source array
  if (mp->noise_tomography == 1) {
    gpuCreateCopy_todevice_realw (&mp->d_noise_sourcearray, noise_sourcearray, NDIM*NGLL3 * (*NSTEP));
  }

  // prepares noise directions
  if (mp->noise_tomography > 1) {
    int nface_size = NGLL2 * mp->nspec2D_top_crust_mantle;

    // allocates memory on GPU
    gpuCreateCopy_todevice_realw (&mp->d_normal_x_noise, normal_x_noise, nface_size);
    gpuCreateCopy_todevice_realw (&mp->d_normal_y_noise, normal_y_noise, nface_size);
    gpuCreateCopy_todevice_realw (&mp->d_normal_z_noise, normal_z_noise, nface_size);

    gpuCreateCopy_todevice_realw (&mp->d_mask_noise, mask_noise, nface_size);
    gpuCreateCopy_todevice_realw (&mp->d_jacobian2D_top_crust_mantle, jacobian2D_top_crust_mantle, nface_size);
  }

  // prepares noise strength kernel
  if (mp->noise_tomography == 3) {
    gpuMalloc_realw (&mp->d_Sigma_kl, NGLL3 * mp->NSPEC_CRUST_MANTLE);

    // initializes kernel values to zero
    gpuMemset_realw (&mp->d_Sigma_kl, NGLL3 * mp->NSPEC_CRUST_MANTLE, 0);
  }

  GPU_ERROR_CHECKING ("prepare_fields_noise_device");
}


/*----------------------------------------------------------------------------------------------- */
// OCEANS
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (prepare_oceans_device,
               PREPARE_OCEANS_DEVICE) (long *Mesh_pointer_f,
                                       int *npoin_oceans,
                                       int *h_iglob_ocean_load,
                                       realw *h_rmass_ocean_load_selected,
                                       realw *h_normal_ocean_load) {

  TRACE ("prepare_oceans_device");
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // arrays with global points on ocean surface
  mp->npoin_oceans = *npoin_oceans;

  // checks for global partitions, each slice must have a top surface with points on i

  if (mp->npoin_oceans == 0) {
    exit_on_error ("prepare_oceans_device has zero npoin_oceans");
  }

  // global point indices
  gpuCreateCopy_todevice_int (&mp->d_ibool_ocean_load, h_iglob_ocean_load, mp->npoin_oceans);

  // mass matrix
  gpuCreateCopy_todevice_realw (&mp->d_rmass_ocean_load, h_rmass_ocean_load_selected, mp->npoin_oceans);

  // normals
  gpuCreateCopy_todevice_realw (&mp->d_normal_ocean_load, h_normal_ocean_load, NDIM * mp->npoin_oceans);

  GPU_ERROR_CHECKING ("prepare_oceans_device");
}

/*----------------------------------------------------------------------------------------------- */
// LDDRK
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (prepare_lddrk_device,
               PREPARE_LDDRK_DEVICE) (long *Mesh_pointer_f) {

  // prepares LDDRK time scheme arrays on GPU
  TRACE ("prepare_lddrk_device");
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // sets flag
  mp->use_lddrk = 1;

  exit_on_error ("prepare_lddrk_device not implemented yet");

  GPU_ERROR_CHECKING ("prepare_lddrk_device");
}



/*----------------------------------------------------------------------------------------------- */
// Earth regions
// CRUST / MANTLE
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (prepare_crust_mantle_device,
               PREPARE_CRUST_MANTLE_DEVICE) (long *Mesh_pointer_f,
                                             realw *h_xix, realw *h_xiy, realw *h_xiz,
                                             realw *h_etax, realw *h_etay, realw *h_etaz,
                                             realw *h_gammax, realw *h_gammay, realw *h_gammaz,
                                             realw *h_rho,
                                             realw *h_kappav, realw *h_muv,
                                             realw *h_kappah, realw *h_muh,
                                             realw *h_eta_aniso,
                                             realw *h_rmassx, realw *h_rmassy, realw *h_rmassz,
                                             realw *h_b_rmassx, realw *h_b_rmassy,
                                             int *h_ibool,
                                             realw *h_xstore, realw *h_ystore, realw *h_zstore,
                                             int *h_ispec_is_tiso,
                                             realw *c11store, realw *c12store, realw *c13store,
                                             realw *c14store, realw *c15store, realw *c16store,
                                             realw *c22store, realw *c23store, realw *c24store,
                                             realw *c25store, realw *c26store, realw *c33store,
                                             realw *c34store, realw *c35store, realw *c36store,
                                             realw *c44store, realw *c45store, realw *c46store,
                                             realw *c55store, realw *c56store, realw *c66store,
                                             int *num_phase_ispec,
                                             int *phase_ispec_inner,
                                             int *nspec_outer,
                                             int *nspec_inner,
                                             int *NSPEC2D_BOTTOM_CM,
                                             int *h_ibelm_bottom_crust_mantle,
                                             int *NCHUNKS_VAL,
                                             int *num_colors_outer,
                                             int *num_colors_inner,
                                             int *num_elem_colors) {

  TRACE ("prepare_crust_mantle_device");

  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  /*Assuming NGLLX=5. Padded is then 128 (5^3+3) */
  int size_padded = NGLL3_PADDED * (mp->NSPEC_CRUST_MANTLE);
  int size_glob = mp->NGLOB_CRUST_MANTLE;

  // mesh
  TRACE ("prepare_crust_mantle gll mesh");
  gpuMalloc_realw (&mp->d_xix_crust_mantle, size_padded);
  gpuMalloc_realw (&mp->d_xiy_crust_mantle, size_padded);
  gpuMalloc_realw (&mp->d_xiz_crust_mantle, size_padded);

  gpuMalloc_realw (&mp->d_etax_crust_mantle, size_padded);
  gpuMalloc_realw (&mp->d_etay_crust_mantle, size_padded);
  gpuMalloc_realw (&mp->d_etaz_crust_mantle, size_padded);

  gpuMalloc_realw (&mp->d_gammax_crust_mantle, size_padded);
  gpuMalloc_realw (&mp->d_gammay_crust_mantle, size_padded);
  gpuMalloc_realw (&mp->d_gammaz_crust_mantle, size_padded);

  // transfer constant element data with padding
  int i;
  for (i = 0; i < mp->NSPEC_CRUST_MANTLE; i++) {
#ifdef USE_OPENCL
    if (run_opencl) {
      int offset = i * NGLL3_PADDED * sizeof (realw);
      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_xix_crust_mantle.ocl, CL_FALSE, offset,
                                     NGLL3 * sizeof (realw), &h_xix[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_xiy_crust_mantle.ocl, CL_FALSE, offset,
                                     NGLL3 * sizeof (realw), &h_xiy[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_xiz_crust_mantle.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_xiz[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_etax_crust_mantle.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_etax[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_etay_crust_mantle.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_etay[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_etaz_crust_mantle.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_etaz[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_gammax_crust_mantle.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_gammax[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_gammay_crust_mantle.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_gammay[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_gammaz_crust_mantle.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw),&h_gammaz[i*NGLL3], 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_xix_crust_mantle.cuda+i*NGLL3_PADDED, &h_xix[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1501);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_xiy_crust_mantle.cuda+i*NGLL3_PADDED, &h_xiy[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1502);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_xiz_crust_mantle.cuda+i*NGLL3_PADDED, &h_xiz[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1503);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_etax_crust_mantle.cuda+i*NGLL3_PADDED, &h_etax[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1504);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_etay_crust_mantle.cuda+i*NGLL3_PADDED, &h_etay[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1505);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_etaz_crust_mantle.cuda+i*NGLL3_PADDED, &h_etaz[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1506);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_gammax_crust_mantle.cuda+i*NGLL3_PADDED, &h_gammax[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1507);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_gammay_crust_mantle.cuda+i*NGLL3_PADDED, &h_gammay[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1508);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_gammaz_crust_mantle.cuda+i*NGLL3_PADDED, &h_gammaz[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1509);
    }
#endif
  }

  // global indexing
  TRACE ("prepare_crust_mantle global indexing");
  gpuCreateCopy_todevice_int (&mp->d_ibool_crust_mantle, h_ibool, NGLL3 * mp->NSPEC_CRUST_MANTLE);

  // Earth model arrays
  TRACE ("prepare_crust_mantle model arrays");
  if (! mp->anisotropic_3D_mantle) {
    // no anisotropy
    // only needed if not anisotropic 3D mantle

    // transverse isotropic elements

    // transverse isotropy flag
    gpuCreateCopy_todevice_int (&mp->d_ispec_is_tiso_crust_mantle, h_ispec_is_tiso, mp->NSPEC_CRUST_MANTLE);

    // kappavstore, kappahstore
    gpuMalloc_realw (&mp->d_kappavstore_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_kappahstore_crust_mantle, size_padded);

    // muvstore, muhstore
    gpuMalloc_realw (&mp->d_muvstore_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_muhstore_crust_mantle, size_padded);

    // eta_anisostore
    gpuMalloc_realw (&mp->d_eta_anisostore_crust_mantle, size_padded);

    // transfer with padding
    int i;
    for (i = 0; i < mp->NSPEC_CRUST_MANTLE; i++) {
#ifdef USE_OPENCL
      if (run_opencl) {
        int offset = i * NGLL3_PADDED * sizeof (realw);
        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_kappavstore_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &h_kappav[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_kappahstore_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &h_kappah[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_muvstore_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &h_muv[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_muhstore_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &h_muh[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_eta_anisostore_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &h_eta_aniso[i*NGLL3], 0, NULL, NULL));
      }
#endif
#ifdef USE_CUDA
      if (run_cuda) {
        print_CUDA_error_if_any(cudaMemcpy(mp->d_kappavstore_crust_mantle.cuda+i*NGLL3_PADDED,&h_kappav[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1510);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_kappahstore_crust_mantle.cuda+i*NGLL3_PADDED,&h_kappah[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1511);

        print_CUDA_error_if_any(cudaMemcpy(mp->d_muvstore_crust_mantle.cuda+i*NGLL3_PADDED,  &h_muv[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1512);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_muhstore_crust_mantle.cuda+i*NGLL3_PADDED,&h_muh[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1513);

        print_CUDA_error_if_any(cudaMemcpy(mp->d_eta_anisostore_crust_mantle.cuda+i*NGLL3_PADDED,&h_eta_aniso[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1514);

      }
#endif
    }
  } else {
    // anisotropic 3D mantle

    // allocates memory on GPU
    gpuMalloc_realw (&mp->d_c11store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c12store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c13store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c14store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c15store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c16store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c22store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c23store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c24store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c25store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c26store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c33store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c34store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c35store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c36store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c44store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c45store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c46store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c55store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c56store_crust_mantle, size_padded);
    gpuMalloc_realw (&mp->d_c66store_crust_mantle, size_padded);

    // transfer constant element data with padding
    int i;
    for (i = 0; i < mp->NSPEC_CRUST_MANTLE; i++) {
#ifdef USE_OPENCL
      if (run_opencl) {
        int offset = i * NGLL3_PADDED * sizeof (realw);

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c11store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw),&c11store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c12store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c12store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c13store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c13store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c14store_crust_mantle.ocl, CL_FALSE,  offset,
                                       NGLL3*sizeof (realw), &c14store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c15store_crust_mantle.ocl, CL_FALSE,  offset,
                                       NGLL3*sizeof (realw), &c15store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c16store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c16store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c22store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c22store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c23store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c23store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c24store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c24store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c25store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c25store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c26store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c26store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c33store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c33store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c34store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c34store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c35store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c35store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c36store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c36store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c44store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c44store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c45store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c45store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c46store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c46store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c55store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c55store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c56store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c56store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c66store_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c66store[i*NGLL3], 0, NULL, NULL));
      }
#endif
#ifdef USE_CUDA
      if (run_cuda) {
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c11store_crust_mantle.cuda + i*NGLL3_PADDED, &c11store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4800);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c12store_crust_mantle.cuda + i*NGLL3_PADDED, &c12store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4801);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c13store_crust_mantle.cuda + i*NGLL3_PADDED, &c13store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4802);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c14store_crust_mantle.cuda + i*NGLL3_PADDED, &c14store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4803);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c15store_crust_mantle.cuda + i*NGLL3_PADDED, &c15store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4804);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c16store_crust_mantle.cuda + i*NGLL3_PADDED, &c16store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4805);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c22store_crust_mantle.cuda + i*NGLL3_PADDED, &c22store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4806);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c23store_crust_mantle.cuda + i*NGLL3_PADDED, &c23store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4807);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c24store_crust_mantle.cuda + i*NGLL3_PADDED, &c24store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4808);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c25store_crust_mantle.cuda + i*NGLL3_PADDED, &c25store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4809);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c26store_crust_mantle.cuda + i*NGLL3_PADDED, &c26store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4810);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c33store_crust_mantle.cuda + i*NGLL3_PADDED, &c33store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4811);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c34store_crust_mantle.cuda + i*NGLL3_PADDED, &c34store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4812);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c35store_crust_mantle.cuda + i*NGLL3_PADDED, &c35store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4813);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c36store_crust_mantle.cuda + i*NGLL3_PADDED, &c36store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4814);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c44store_crust_mantle.cuda + i*NGLL3_PADDED, &c44store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4815);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c45store_crust_mantle.cuda + i*NGLL3_PADDED, &c45store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4816);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c46store_crust_mantle.cuda + i*NGLL3_PADDED, &c46store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4817);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c55store_crust_mantle.cuda + i*NGLL3_PADDED, &c55store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4818);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c56store_crust_mantle.cuda + i*NGLL3_PADDED, &c56store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4819);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c66store_crust_mantle.cuda + i*NGLL3_PADDED, &c66store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4820);
      }
#endif
    }
  }

  // needed for boundary kernel calculations
  if (mp->simulation_type == 3 && mp->save_boundary_mesh) {
    gpuMalloc_realw (&mp->d_rhostore_crust_mantle, size_padded);

    int i;
    for (i = 0; i < mp->NSPEC_CRUST_MANTLE; i++) {
#ifdef USE_OPENCL
      if (run_opencl) {
        int offset = i * NGLL3_PADDED * sizeof(realw);
        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_rhostore_crust_mantle.ocl, CL_FALSE, offset,
                                       NGLL3 * sizeof (realw), &h_rho[i * NGLL3], 0, NULL, NULL));
      }
#endif
#ifdef USE_CUDA
      if (run_cuda) {
        print_CUDA_error_if_any(cudaMemcpy(mp->d_rhostore_crust_mantle.cuda+i*NGLL3_PADDED, &h_rho[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),2106);
      }
#endif
    }
  }

  // mesh locations
  TRACE ("prepare_crust_mantle mesh locations");

  // ystore & zstore needed for tiso elements
  gpuCreateCopy_todevice_realw (&mp->d_ystore_crust_mantle, h_ystore, size_glob);
  gpuCreateCopy_todevice_realw (&mp->d_zstore_crust_mantle, h_zstore, size_glob);

  // xstore only needed when gravity is on
  if (mp->gravity) {
    gpuCreateCopy_todevice_realw (&mp->d_xstore_crust_mantle, h_xstore, size_glob);
  }

  // inner/outer elements
  mp->num_phase_ispec_crust_mantle = *num_phase_ispec;
  gpuCreateCopy_todevice_int (&mp->d_phase_ispec_inner_crust_mantle, phase_ispec_inner, mp->num_phase_ispec_crust_mantle*2);

  mp->nspec_outer_crust_mantle = *nspec_outer;
  mp->nspec_inner_crust_mantle = *nspec_inner;

  // CMB/fluid outer core coupling
  mp->nspec2D_bottom_crust_mantle = *NSPEC2D_BOTTOM_CM;
  gpuCreateCopy_todevice_int (&mp->d_ibelm_bottom_crust_mantle, h_ibelm_bottom_crust_mantle, mp->nspec2D_bottom_crust_mantle);

  // wavefield
  TRACE ("prepare_crust_mantle wavefields");
  int size = NDIM * mp->NGLOB_CRUST_MANTLE;

  gpuMalloc_realw (&mp->d_displ_crust_mantle, size);
  gpuMalloc_realw (&mp->d_veloc_crust_mantle, size);
  gpuMalloc_realw (&mp->d_accel_crust_mantle, size);
  // backward/reconstructed wavefield
  if (mp->simulation_type == 3) {
    gpuMalloc_realw (&mp->d_b_displ_crust_mantle, size);
    gpuMalloc_realw (&mp->d_b_veloc_crust_mantle, size);
    gpuMalloc_realw (&mp->d_b_accel_crust_mantle, size);

    //debugging with empty arrays
    if (DEBUG_BACKWARD_SIMULATIONS == 1) {
      gpuMemset_realw (&mp->d_b_displ_crust_mantle, size, 0);
      gpuMemset_realw (&mp->d_b_veloc_crust_mantle, size, 0);
      gpuMemset_realw (&mp->d_b_accel_crust_mantle, size, 0);
    }
  }

#ifdef USE_OPENCL
  if (run_opencl) {
#ifdef USE_TEXTURES_FIELDS
    cl_int errcode;
    cl_image_format format = {CL_R, CL_UNSIGNED_INT32};

    mp->d_displ_cm_tex = clCreateImage2D (mocl.context, CL_MEM_READ_ONLY, &format, size, 1, 0, mp->d_displ_crust_mantle.ocl, clck_(&errcode));
    mp->d_accel_cm_tex = clCreateImage2D (mocl.context, CL_MEM_READ_ONLY, &format, size, 1, 0, mp->d_accel_crust_mantle.ocl, clck_(&errcode));
    // backward/reconstructed fields
    if (mp->simulation_type == 3) {
      mp->d_b_displ_cm_tex = clCreateImage2D (mocl.context, CL_MEM_READ_ONLY, &format, size, 1, 0, mp->d_b_displ_crust_mantle.ocl, clck_(&errcode));
      mp->d_b_accel_cm_tex = clCreateImage2D (mocl.context, CL_MEM_READ_ONLY, &format, size, 1, 0, mp->d_b_accel_crust_mantle.ocl, clck_(&errcode));
    } else {
      mp->d_b_displ_cm_tex = moclGetDummyImage2D(mp);
      mp->d_b_accel_cm_tex = moclGetDummyImage2D(mp);
    }
#endif
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
#ifdef USE_TEXTURES_FIELDS
    {
    // checks single precision
    if (sizeof(realw) != sizeof(float)) exit_on_error("TEXTURES only work with realw selected as float");
    // binds textures
#ifdef USE_OLDER_CUDA4_GPU
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<realw>();
      const textureReference* d_displ_cm_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_displ_cm_tex_ref_ptr, "d_displ_cm_tex"), 4021);
      print_CUDA_error_if_any(cudaBindTexture(0, d_displ_cm_tex_ref_ptr, mp->d_displ_crust_mantle.cuda,
                                              &channelDesc, sizeof(realw)*size), 4021);

      const textureReference* d_accel_cm_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_accel_cm_tex_ref_ptr, "d_accel_cm_tex"), 4023);
      print_CUDA_error_if_any(cudaBindTexture(0, d_accel_cm_tex_ref_ptr, mp->d_accel_crust_mantle.cuda,
                                              &channelDesc, sizeof(realw)*size), 4023);

      // backward/reconstructed wavefields
      if (mp->simulation_type == 3) {
        const textureReference* d_b_displ_cm_tex_ref_ptr;
        print_CUDA_error_if_any(cudaGetTextureReference(&d_b_displ_cm_tex_ref_ptr, "d_b_displ_cm_tex"), 4021);
        print_CUDA_error_if_any(cudaBindTexture(0, d_b_displ_cm_tex_ref_ptr, mp->d_b_displ_crust_mantle.cuda,
                                                &channelDesc, sizeof(realw)*size), 4021);

        const textureReference* d_b_accel_cm_tex_ref_ptr;
        print_CUDA_error_if_any(cudaGetTextureReference(&d_b_accel_cm_tex_ref_ptr, "d_b_accel_cm_tex"), 4023);
        print_CUDA_error_if_any(cudaBindTexture(0, d_b_accel_cm_tex_ref_ptr, mp->d_b_accel_crust_mantle.cuda,
                                                &channelDesc, sizeof(realw)*size), 4023);
      }
#else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<realw>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_displ_cm_tex, mp->d_displ_crust_mantle.cuda,
                                              &channelDesc, sizeof(realw)*size), 4021);
      print_CUDA_error_if_any(cudaBindTexture(0, &d_accel_cm_tex, mp->d_accel_crust_mantle.cuda,
                                              &channelDesc, sizeof(realw)*size), 4023);
      // backward/reconstructed wavefields
      if (mp->simulation_type == 3) {
        print_CUDA_error_if_any(cudaBindTexture(0, &d_b_displ_cm_tex, mp->d_b_displ_crust_mantle.cuda,
                                                &channelDesc, sizeof(realw)*size), 4021);
        print_CUDA_error_if_any(cudaBindTexture(0, &d_b_accel_cm_tex, mp->d_b_accel_crust_mantle.cuda,
                                                &channelDesc, sizeof(realw)*size), 4023);
      }
#endif
    }
#endif
  }
#endif

  // mass matrices
  TRACE ("prepare_crust_mantle mass matrices");
  gpuCreateCopy_todevice_realw (&mp->d_rmassz_crust_mantle, h_rmassz, size_glob);
  if ((*NCHUNKS_VAL != 6 && mp->absorbing_conditions) || (mp->rotation && mp->exact_mass_matrix_for_rotation)) {
    gpuCreateCopy_todevice_realw (&mp->d_rmassx_crust_mantle, h_rmassx, size_glob);
    gpuCreateCopy_todevice_realw (&mp->d_rmassy_crust_mantle, h_rmassy, size_glob);
  } else {
    mp->d_rmassx_crust_mantle = gpuTakeRef(mp->d_rmassz_crust_mantle);
    mp->d_rmassy_crust_mantle = gpuTakeRef(mp->d_rmassz_crust_mantle);
  }

  // kernel simulations
  if (mp->simulation_type == 3) {
    mp->d_b_rmassz_crust_mantle = gpuTakeRef(mp->d_rmassz_crust_mantle);
    if (mp->rotation && mp->exact_mass_matrix_for_rotation) {
      gpuCreateCopy_todevice_realw (&mp->d_b_rmassx_crust_mantle, h_b_rmassx, size_glob);
      gpuCreateCopy_todevice_realw (&mp->d_b_rmassy_crust_mantle, h_b_rmassy, size_glob);
    } else {
      mp->d_b_rmassx_crust_mantle = gpuTakeRef(mp->d_rmassx_crust_mantle);
      mp->d_b_rmassy_crust_mantle = gpuTakeRef(mp->d_rmassy_crust_mantle);
    }

    // kernels
    size = NGLL3 * (mp->NSPEC_CRUST_MANTLE);

    // density kernel
    gpuMalloc_realw (&mp->d_rho_kl_crust_mantle, size);
    // initializes kernel values to zero
    gpuMemset_realw (&mp->d_rho_kl_crust_mantle, size, 0);

    // wavespeed kernels
    if (! mp->anisotropic_kl) {
      // isotropic kernels
      gpuMalloc_realw (&mp->d_alpha_kl_crust_mantle, size);
      gpuMalloc_realw (&mp->d_beta_kl_crust_mantle, size);
      // sets array values to zero
      gpuMemset_realw (&mp->d_alpha_kl_crust_mantle, size, 0);
      gpuMemset_realw (&mp->d_beta_kl_crust_mantle, size, 0);

    } else {
      // anisotropic kernels
      gpuMalloc_realw (&mp->d_cijkl_kl_crust_mantle, 21 * size);
      gpuMemset_realw (&mp->d_cijkl_kl_crust_mantle, 21 * size, 0);
    }

    // preconditioner
    if (mp->approximate_hess_kl) {
      gpuMalloc_realw (&mp->d_hess_kl_crust_mantle, size);
      gpuMemset_realw (&mp->d_hess_kl_crust_mantle, size, 0);
    }
  }

  // mesh coloring
  mp->num_colors_outer_crust_mantle = *num_colors_outer;
  mp->num_colors_inner_crust_mantle = *num_colors_inner;
  mp->h_num_elem_colors_crust_mantle = (int *) num_elem_colors;

  // executes gpu calls
  // (especially needed for OpenCL kernels to finish executing)
  gpuSynchronize();

  GPU_ERROR_CHECKING ("prepare_crust_mantle_device");
  // debug
  //printf("%d rank - prepare_crust_mantle done",mp->myrank);
}


/*----------------------------------------------------------------------------------------------- */
// OUTER CORE
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (prepare_outer_core_device,
               PREPARE_OUTER_CORE_DEVICE) (long *Mesh_pointer_f,
                                           realw *h_xix, realw *h_xiy, realw *h_xiz,
                                           realw *h_etax, realw *h_etay, realw *h_etaz,
                                           realw *h_gammax, realw *h_gammay, realw *h_gammaz,
                                           realw *h_rho, realw *h_kappav,
                                           realw *h_rmass,
                                           int *h_ibool,
                                           realw *h_xstore, realw *h_ystore, realw *h_zstore,
                                           int *num_phase_ispec,
                                           int *phase_ispec_inner,
                                           int *nspec_outer,
                                           int *nspec_inner,
                                           int *NSPEC2D_TOP_OC,
                                           int *NSPEC2D_BOTTOM_OC,
                                           realw *h_normal_top_outer_core,
                                           realw *h_normal_bottom_outer_core,
                                           realw *h_jacobian2D_top_outer_core,
                                           realw *h_jacobian2D_bottom_outer_core,
                                           int *h_ibelm_top_outer_core,
                                           int *h_ibelm_bottom_outer_core,
                                           int *num_colors_outer,
                                           int *num_colors_inner,
                                           int *num_elem_colors) {

  TRACE ("prepare_outer_core_device");

  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  /*Assuming NGLLX=5. Padded is then 128 (5^3+3) */

  int size_padded = NGLL3_PADDED * (mp->NSPEC_OUTER_CORE);
  int size_glob = mp->NGLOB_OUTER_CORE;

  // mesh
  gpuMalloc_realw (&mp->d_xix_outer_core, size_padded);
  gpuMalloc_realw (&mp->d_xiy_outer_core, size_padded);
  gpuMalloc_realw (&mp->d_xiz_outer_core, size_padded);
  gpuMalloc_realw (&mp->d_etax_outer_core, size_padded);
  gpuMalloc_realw (&mp->d_etay_outer_core, size_padded);
  gpuMalloc_realw (&mp->d_etaz_outer_core, size_padded);
  gpuMalloc_realw (&mp->d_gammax_outer_core, size_padded);
  gpuMalloc_realw (&mp->d_gammay_outer_core, size_padded);
  gpuMalloc_realw (&mp->d_gammaz_outer_core, size_padded);

  gpuMalloc_realw (&mp->d_kappavstore_outer_core, size_padded);

  // transfer constant element data with padding
  int i;
  for (i = 0; i < mp->NSPEC_OUTER_CORE; i++) {
#ifdef USE_OPENCL
    if (run_opencl) {
      int offset = i * NGLL3_PADDED * sizeof (realw);
      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_xix_outer_core.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_xix[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_xiy_outer_core.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_xiy[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_xiz_outer_core.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_xiz[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_etax_outer_core.ocl, CL_FALSE, offset,
                                     NGLL3 * sizeof (realw), &h_etax[i * NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_etay_outer_core.ocl, CL_FALSE, offset,
                                     NGLL3 * sizeof (realw), &h_etay[i * NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_etaz_outer_core.ocl, CL_FALSE, offset,
                                     NGLL3 * sizeof (realw), &h_etaz[i * NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_gammax_outer_core.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_gammax[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_gammay_outer_core.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_gammay[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_gammaz_outer_core.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_gammaz[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_kappavstore_outer_core.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_kappav[i*NGLL3], 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_xix_outer_core.cuda+i*NGLL3_PADDED, &h_xix[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1501);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_xiy_outer_core.cuda+i*NGLL3_PADDED, &h_xiy[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1502);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_xiz_outer_core.cuda+i*NGLL3_PADDED, &h_xiz[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1503);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_etax_outer_core.cuda+i*NGLL3_PADDED, &h_etax[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1504);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_etay_outer_core.cuda+i*NGLL3_PADDED, &h_etay[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1505);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_etaz_outer_core.cuda+i*NGLL3_PADDED, &h_etaz[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1506);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_gammax_outer_core.cuda+i*NGLL3_PADDED, &h_gammax[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1507);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_gammay_outer_core.cuda+i*NGLL3_PADDED, &h_gammay[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1508);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_gammaz_outer_core.cuda+i*NGLL3_PADDED, &h_gammaz[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1509);

      print_CUDA_error_if_any(cudaMemcpy(mp->d_kappavstore_outer_core.cuda+i*NGLL3_PADDED,&h_kappav[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1510);
    }
#endif
  }

  // needed for kernel calculations
  if (mp->simulation_type == 3) {
    gpuMalloc_realw (&mp->d_rhostore_outer_core, size_padded);

    int i;
    for (i = 0; i < mp->NSPEC_OUTER_CORE; i++) {
#ifdef USE_OPENCL
      if (run_opencl) {
        int offset = i * NGLL3_PADDED * sizeof (realw);
        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_rhostore_outer_core.ocl, CL_FALSE, offset,
                                       NGLL3 * sizeof (realw), &h_rho[i*NGLL3], 0, NULL, NULL));
      }
#endif
#ifdef USE_CUDA
      if (run_cuda) {
        print_CUDA_error_if_any(cudaMemcpy(mp->d_rhostore_outer_core.cuda+i*NGLL3_PADDED, &h_rho[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),2106);
      }
#endif
    }
  }

  // global indexing
  gpuCreateCopy_todevice_int (&mp->d_ibool_outer_core, h_ibool, NGLL3 * (mp->NSPEC_OUTER_CORE));

  // mesh locations

  // always needed
  gpuCreateCopy_todevice_realw (&mp->d_xstore_outer_core, h_xstore, size_glob);
  gpuCreateCopy_todevice_realw (&mp->d_ystore_outer_core, h_ystore, size_glob);
  gpuCreateCopy_todevice_realw (&mp->d_zstore_outer_core, h_zstore, size_glob);

  // inner/outer elements

  mp->num_phase_ispec_outer_core = *num_phase_ispec;
  gpuCreateCopy_todevice_int (&mp->d_phase_ispec_inner_outer_core, phase_ispec_inner, mp->num_phase_ispec_outer_core*2);

  mp->nspec_outer_outer_core = *nspec_outer;
  mp->nspec_inner_outer_core = *nspec_inner;

  // CMB/ICB coupling

  mp->nspec2D_top_outer_core = *NSPEC2D_TOP_OC;
  mp->nspec2D_bottom_outer_core = *NSPEC2D_BOTTOM_OC;

  int size_toc = NGLL2 * mp->nspec2D_top_outer_core;
  gpuCreateCopy_todevice_int (&mp->d_ibelm_top_outer_core, h_ibelm_top_outer_core, mp->nspec2D_top_outer_core);
  gpuCreateCopy_todevice_realw (&mp->d_jacobian2D_top_outer_core, h_jacobian2D_top_outer_core, size_toc);
  gpuCreateCopy_todevice_realw (&mp->d_normal_top_outer_core, h_normal_top_outer_core, NDIM*size_toc);

  int size_boc = NGLL2 * mp->nspec2D_bottom_outer_core;
  gpuCreateCopy_todevice_int (&mp->d_ibelm_bottom_outer_core, h_ibelm_bottom_outer_core, mp->nspec2D_bottom_outer_core);
  gpuCreateCopy_todevice_realw (&mp->d_jacobian2D_bottom_outer_core, h_jacobian2D_bottom_outer_core, size_boc);
  gpuCreateCopy_todevice_realw (&mp->d_normal_bottom_outer_core, h_normal_bottom_outer_core, NDIM*size_boc);

  // wavefield
  gpuMalloc_realw (&mp->d_displ_outer_core, size_glob);
  gpuMalloc_realw (&mp->d_veloc_outer_core, size_glob);
  gpuMalloc_realw (&mp->d_accel_outer_core, size_glob);
  // backward/reconstructed wavefield
  if (mp->simulation_type == 3) {
    gpuMalloc_realw (&mp->d_b_displ_outer_core, size_glob);
    gpuMalloc_realw (&mp->d_b_veloc_outer_core, size_glob);
    gpuMalloc_realw (&mp->d_b_accel_outer_core, size_glob);

    //debugging with empty arrays
    if (DEBUG_BACKWARD_SIMULATIONS == 1) {
      gpuMemset_realw (&mp->d_b_displ_outer_core, size_glob, 0);
      gpuMemset_realw (&mp->d_b_veloc_outer_core, size_glob, 0);
      gpuMemset_realw (&mp->d_b_accel_outer_core, size_glob, 0);
    }
  }

#ifdef USE_OPENCL
  if (run_opencl) {
#ifdef USE_TEXTURES_FIELDS
    cl_int errcode;
    cl_image_format format = {CL_R, CL_UNSIGNED_INT32};

    mp->d_displ_oc_tex = clCreateImage2D (mocl.context, CL_MEM_READ_ONLY, &format, size_glob, 1, 0, mp->d_displ_outer_core.ocl, clck_(&errcode));
    mp->d_accel_oc_tex = clCreateImage2D (mocl.context, CL_MEM_READ_ONLY, &format, size_glob, 1, 0, mp->d_accel_outer_core.ocl, clck_(&errcode));
    // backward/reconstructed fields
    if (mp->simulation_type == 3) {
      mp->d_b_displ_oc_tex = clCreateImage2D (mocl.context, CL_MEM_READ_ONLY, &format, size_glob, 1, 0, mp->d_b_displ_outer_core.ocl, clck_(&errcode));
      mp->d_b_accel_oc_tex = clCreateImage2D (mocl.context, CL_MEM_READ_ONLY, &format, size_glob, 1, 0, mp->d_b_accel_outer_core.ocl, clck_(&errcode));
    } else {
      mp->d_b_displ_oc_tex = moclGetDummyImage2D(mp);
      mp->d_b_accel_oc_tex = moclGetDummyImage2D(mp);
    }
#endif
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
#ifdef USE_TEXTURES_FIELDS
#ifdef USE_OLDER_CUDA4_GPU
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<realw>();
    const textureReference* d_displ_oc_tex_ref_ptr;
    print_CUDA_error_if_any(cudaGetTextureReference(&d_displ_oc_tex_ref_ptr, "d_displ_oc_tex"), 5021);
    print_CUDA_error_if_any(cudaBindTexture(0, d_displ_oc_tex_ref_ptr, mp->d_displ_outer_core.cuda,
                                            &channelDesc, sizeof(realw)*size_glob), 5021);

    const textureReference* d_accel_oc_tex_ref_ptr;
    print_CUDA_error_if_any(cudaGetTextureReference(&d_accel_oc_tex_ref_ptr, "d_accel_oc_tex"), 5023);
    print_CUDA_error_if_any(cudaBindTexture(0, d_accel_oc_tex_ref_ptr, mp->d_accel_outer_core.cuda,
                                            &channelDesc, sizeof(realw)*size_glob), 5023);
    // backward/reconstructed wavefields
    if (mp->simulation_type == 3) {
      const textureReference* d_b_displ_oc_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_b_displ_oc_tex_ref_ptr, "d_b_displ_oc_tex"), 5021);
      print_CUDA_error_if_any(cudaBindTexture(0, d_b_displ_oc_tex_ref_ptr, mp->d_b_displ_outer_core.cuda,
                                              &channelDesc, sizeof(realw)*size_glob), 5021);

      const textureReference* d_b_accel_oc_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_b_accel_oc_tex_ref_ptr, "d_b_accel_oc_tex"), 5023);
        print_CUDA_error_if_any(cudaBindTexture(0, d_b_accel_oc_tex_ref_ptr, mp->d_b_accel_outer_core.cuda,
                                                &channelDesc, sizeof(realw)*size_glob), 5023);
    }
#else
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<realw>();
    print_CUDA_error_if_any(cudaBindTexture(0, &d_displ_oc_tex, mp->d_displ_outer_core.cuda,
                                            &channelDesc, sizeof(realw)*size_glob), 5021);
    print_CUDA_error_if_any(cudaBindTexture(0, &d_accel_oc_tex, mp->d_accel_outer_core.cuda,
                                            &channelDesc, sizeof(realw)*size_glob), 5023);
    // backward/reconstructed wavefields
    if (mp->simulation_type == 3) {
      print_CUDA_error_if_any(cudaBindTexture(0, &d_b_displ_oc_tex, mp->d_b_displ_outer_core.cuda,
                                              &channelDesc, sizeof(realw)*size_glob), 5021);
      print_CUDA_error_if_any(cudaBindTexture(0, &d_b_accel_oc_tex, mp->d_b_accel_outer_core.cuda,
                                              &channelDesc, sizeof(realw)*size_glob), 5023);
    }
#endif
#endif
  }
#endif

  // mass matrix
  gpuCreateCopy_todevice_realw (&mp->d_rmass_outer_core, h_rmass, size_glob);

  // kernel simulations
  if (mp->simulation_type == 3) {
    // mass matrix
    mp->d_b_rmass_outer_core = gpuTakeRef(mp->d_rmass_outer_core);

    //kernels
    int size = NGLL3 * (mp->NSPEC_OUTER_CORE);

    // density kernel
    gpuMalloc_realw (&mp->d_rho_kl_outer_core, size);
    gpuMemset_realw (&mp->d_rho_kl_outer_core, size, 0);

    // isotropic kernel
    gpuMalloc_realw (&mp->d_alpha_kl_outer_core, size);
    gpuMemset_realw (&mp->d_alpha_kl_outer_core, size, 0);
  }

  // mesh coloring
  mp->num_colors_outer_outer_core = *num_colors_outer;
  mp->num_colors_inner_outer_core = *num_colors_inner;
  mp->h_num_elem_colors_outer_core = (int *) num_elem_colors;

  // executes gpu calls
  // (especially needed for OpenCL kernels to finish executing)
  gpuSynchronize();

  GPU_ERROR_CHECKING ("prepare_outer_core_device");
}

/*----------------------------------------------------------------------------------------------- */
// INNER CORE
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (prepare_inner_core_device,
               PREPARE_INNER_CORE_DEVICE) (long *Mesh_pointer_f,
                                           realw *h_xix, realw *h_xiy, realw *h_xiz,
                                           realw *h_etax, realw *h_etay, realw *h_etaz,
                                           realw *h_gammax, realw *h_gammay, realw *h_gammaz,
                                           realw *h_rho, realw *h_kappav, realw *h_muv,
                                           realw *h_rmassx, realw *h_rmassy, realw *h_rmassz,
                                           realw *h_b_rmassx, realw *h_b_rmassy,
                                           int *h_ibool,
                                           realw *h_xstore, realw *h_ystore, realw *h_zstore,
                                           realw *c11store, realw *c12store, realw *c13store,
                                           realw *c33store, realw *c44store,
                                           int *h_idoubling_inner_core,
                                           int *num_phase_ispec,
                                           int *phase_ispec_inner,
                                           int *nspec_outer,
                                           int *nspec_inner,
                                           int *NSPEC2D_TOP_IC,
                                           int *h_ibelm_top_inner_core,
                                           int *num_colors_outer,
                                           int *num_colors_inner,
                                           int *num_elem_colors) {

  TRACE ("prepare_inner_core_device");

  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  /* Assuming NGLLX=5. Padded is then 128 (5^3+3) */
  int size_padded = NGLL3_PADDED * (mp->NSPEC_INNER_CORE);
  int size_glob = mp->NGLOB_INNER_CORE;

  // mesh
  gpuMalloc_realw (&mp->d_xix_inner_core, size_padded);
  gpuMalloc_realw (&mp->d_xiy_inner_core, size_padded);
  gpuMalloc_realw (&mp->d_xiz_inner_core, size_padded);
  gpuMalloc_realw (&mp->d_etax_inner_core, size_padded);
  gpuMalloc_realw (&mp->d_etay_inner_core, size_padded);
  gpuMalloc_realw (&mp->d_etaz_inner_core, size_padded);
  gpuMalloc_realw (&mp->d_gammax_inner_core, size_padded);
  gpuMalloc_realw (&mp->d_gammay_inner_core, size_padded);
  gpuMalloc_realw (&mp->d_gammaz_inner_core, size_padded);

  // muvstore needed for attenuation also for anisotropic inner core
  gpuMalloc_realw (&mp->d_muvstore_inner_core, size_padded);

  // transfer constant element data with padding
  int i;
  for (i = 0;i < mp->NSPEC_INNER_CORE;i++) {
#ifdef USE_OPENCL
    if (run_opencl) {
      int offset = i * NGLL3_PADDED * sizeof (realw);

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_xix_inner_core.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_xix[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_xiy_inner_core.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_xiy[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_xiz_inner_core.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_xiz[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_etax_inner_core.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_etax[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_etay_inner_core.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_etay[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_etaz_inner_core.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_etaz[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_gammax_inner_core.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_gammax[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_gammay_inner_core.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw),&h_gammay[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_gammaz_inner_core.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_gammaz[i*NGLL3], 0, NULL, NULL));

      clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_muvstore_inner_core.ocl, CL_FALSE, offset,
                                     NGLL3*sizeof (realw), &h_muv[i*NGLL3], 0, NULL, NULL));
    }
#endif
#ifdef USE_CUDA
    if (run_cuda) {
      print_CUDA_error_if_any(cudaMemcpy(mp->d_xix_inner_core.cuda+i*NGLL3_PADDED, &h_xix[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1501);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_xiy_inner_core.cuda+i*NGLL3_PADDED, &h_xiy[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1502);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_xiz_inner_core.cuda+i*NGLL3_PADDED, &h_xiz[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1503);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_etax_inner_core.cuda+i*NGLL3_PADDED, &h_etax[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1504);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_etay_inner_core.cuda+i*NGLL3_PADDED, &h_etay[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1505);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_etaz_inner_core.cuda+i*NGLL3_PADDED, &h_etaz[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1506);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_gammax_inner_core.cuda+i*NGLL3_PADDED, &h_gammax[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1507);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_gammay_inner_core.cuda+i*NGLL3_PADDED, &h_gammay[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1508);
      print_CUDA_error_if_any(cudaMemcpy(mp->d_gammaz_inner_core.cuda+i*NGLL3_PADDED, &h_gammaz[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1509);

      print_CUDA_error_if_any(cudaMemcpy(mp->d_muvstore_inner_core.cuda+i*NGLL3_PADDED, &h_muv[i*NGLL3],
                                         NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1511);
    }
#endif
  }

  // anisotropy
  if (! mp->anisotropic_inner_core) {
    // no anisotropy (uses kappav and muv in inner core)

    // kappavstore needed
    gpuMalloc_realw (&mp->d_kappavstore_inner_core, size_padded);

    int i;
    for (i = 0;i < mp->NSPEC_INNER_CORE;i++) {
#ifdef USE_OPENCL
      if (run_opencl) {
        int offset = i * NGLL3_PADDED * sizeof (realw);
        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_kappavstore_inner_core.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &h_kappav[i*NGLL3], 0, NULL, NULL));
      }
#endif
#ifdef USE_CUDA
      if (run_cuda) {
        print_CUDA_error_if_any(cudaMemcpy(mp->d_kappavstore_inner_core.cuda+i*NGLL3_PADDED,&h_kappav[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),1510);
      }
#endif
    }
  } else {
    // anisotropic inner core
    gpuMalloc_realw (&mp->d_c11store_inner_core, size_padded);
    gpuMalloc_realw (&mp->d_c12store_inner_core, size_padded);
    gpuMalloc_realw (&mp->d_c13store_inner_core, size_padded);
    gpuMalloc_realw (&mp->d_c33store_inner_core, size_padded);
    gpuMalloc_realw (&mp->d_c44store_inner_core, size_padded);

    // transfer constant element data with padding
    int i;
    for (i = 0;i < mp->NSPEC_INNER_CORE; i++) {
#ifdef USE_OPENCL
      if (run_opencl) {
        int offset = i * NGLL3_PADDED * sizeof (realw);
        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c11store_inner_core.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c11store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c12store_inner_core.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c12store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c13store_inner_core.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c13store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c33store_inner_core.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c33store[i*NGLL3], 0, NULL, NULL));

        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_c44store_inner_core.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &c44store[i*NGLL3], 0, NULL, NULL));
      }
#endif
#ifdef USE_CUDA
      if (run_cuda) {
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c11store_inner_core.cuda + i*NGLL3_PADDED, &c11store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4800);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c12store_inner_core.cuda + i*NGLL3_PADDED, &c12store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4801);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c13store_inner_core.cuda + i*NGLL3_PADDED, &c13store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4802);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c33store_inner_core.cuda + i*NGLL3_PADDED, &c33store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4803);
        print_CUDA_error_if_any(cudaMemcpy(mp->d_c44store_inner_core.cuda + i*NGLL3_PADDED, &c44store[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),4804);
      }
#endif
    }
  }

  // needed for boundary kernel calculations
  if (mp->simulation_type == 3 && mp->save_boundary_mesh) {
    gpuMalloc_realw (&mp->d_rhostore_inner_core, size_padded);

    int i;
    for (i = 0; i < mp->NSPEC_INNER_CORE; i++) {
#ifdef USE_OPENCL
      if (run_opencl) {
        int offset = i * NGLL3_PADDED * sizeof (realw);
        clCheck (clEnqueueWriteBuffer (mocl.command_queue, mp->d_rhostore_inner_core.ocl, CL_FALSE, offset,
                                       NGLL3*sizeof (realw), &h_rho[i*NGLL3], 0, NULL, NULL));
      }
#endif
#ifdef USE_CUDA
      if (run_cuda) {
        print_CUDA_error_if_any(cudaMemcpy(mp->d_rhostore_inner_core.cuda+i*NGLL3_PADDED, &h_rho[i*NGLL3],
                                           NGLL3*sizeof(realw),cudaMemcpyHostToDevice),2106);
      }
#endif
    }
  }

  // global indexing
  gpuCreateCopy_todevice_int (&mp->d_ibool_inner_core, h_ibool, NGLL3 * (mp->NSPEC_INNER_CORE));

  // fictious element flags

  gpuCreateCopy_todevice_int (&mp->d_idoubling_inner_core, h_idoubling_inner_core, mp->NSPEC_INNER_CORE);

  // mesh locations

  // only needed when gravity is on

  if (mp->gravity) {
    gpuCreateCopy_todevice_realw (&mp->d_xstore_inner_core, h_xstore, size_glob);
    gpuCreateCopy_todevice_realw (&mp->d_ystore_inner_core, h_ystore, size_glob);
    gpuCreateCopy_todevice_realw (&mp->d_zstore_inner_core, h_zstore, size_glob);
  }

  // inner/outer elements
  mp->num_phase_ispec_inner_core = *num_phase_ispec;
  gpuCreateCopy_todevice_int (&mp->d_phase_ispec_inner_inner_core, phase_ispec_inner, mp->num_phase_ispec_inner_core*2);

  mp->nspec_outer_inner_core = *nspec_outer;
  mp->nspec_inner_inner_core = *nspec_inner;

  // boundary elements on top

  mp->nspec2D_top_inner_core = *NSPEC2D_TOP_IC;
  gpuCreateCopy_todevice_int (&mp->d_ibelm_top_inner_core, h_ibelm_top_inner_core, mp->nspec2D_top_inner_core);

  // wavefield
  int size = NDIM * mp->NGLOB_INNER_CORE;

  gpuMalloc_realw (&mp->d_displ_inner_core, size);
  gpuMalloc_realw (&mp->d_veloc_inner_core, size);
  gpuMalloc_realw (&mp->d_accel_inner_core, size);
  // backward/reconstructed wavefield
  if (mp->simulation_type == 3) {
    gpuMalloc_realw (&mp->d_b_displ_inner_core, size);
    gpuMalloc_realw (&mp->d_b_veloc_inner_core, size);
    gpuMalloc_realw (&mp->d_b_accel_inner_core, size);

    //debugging with empty arrays
    if (DEBUG_BACKWARD_SIMULATIONS == 1) {
      gpuMemset_realw (&mp->d_displ_inner_core, size, 0);
      gpuMemset_realw (&mp->d_veloc_inner_core, size, 0);
      gpuMemset_realw (&mp->d_accel_inner_core, size, 0);
    }
  }

#ifdef USE_OPENCL
  if (run_opencl) {
#ifdef USE_TEXTURES_FIELDS
    cl_int errcode;
    cl_image_format format = {CL_R, CL_UNSIGNED_INT32};
    mp->d_displ_ic_tex = clCreateImage2D (mocl.context, CL_MEM_READ_ONLY, &format, size, 1, 0, mp->d_displ_inner_core.ocl, clck_(&errcode));
    mp->d_accel_ic_tex = clCreateImage2D (mocl.context, CL_MEM_READ_ONLY, &format, size, 1, 0, mp->d_accel_inner_core.ocl, clck_(&errcode));
    // backward/reconstructed fields
    if (mp->simulation_type == 3) {
      mp->d_b_displ_ic_tex = clCreateImage2D (mocl.context, CL_MEM_READ_ONLY, &format, size, 1, 0, mp->d_b_displ_inner_core.ocl, clck_(&errcode));
      mp->d_b_accel_ic_tex = clCreateImage2D (mocl.context, CL_MEM_READ_ONLY, &format, size, 1, 0, mp->d_b_accel_inner_core.ocl, clck_(&errcode));
    } else {
      mp->d_b_displ_ic_tex = moclGetDummyImage2D(mp);
      mp->d_b_accel_ic_tex = moclGetDummyImage2D(mp);
    }
#endif
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
#ifdef USE_TEXTURES_FIELDS
    {
#ifdef USE_OLDER_CUDA4_GPU
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<realw>();
      const textureReference* d_displ_ic_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_displ_ic_tex_ref_ptr, "d_displ_ic_tex"), 6021);
      print_CUDA_error_if_any(cudaBindTexture(0, d_displ_ic_tex_ref_ptr, mp->d_displ_inner_core,
                                              &channelDesc, sizeof(realw)*size), 6021);

      const textureReference* d_accel_ic_tex_ref_ptr;
      print_CUDA_error_if_any(cudaGetTextureReference(&d_accel_ic_tex_ref_ptr, "d_accel_ic_tex"), 6023);
      print_CUDA_error_if_any(cudaBindTexture(0, d_accel_ic_tex_ref_ptr, mp->d_accel_inner_core,
                                              &channelDesc, sizeof(realw)*size), 6023);
      // backward/reconstructed wavefields
      if (mp->simulation_type == 3) {
        const textureReference* d_b_displ_ic_tex_ref_ptr;
        print_CUDA_error_if_any(cudaGetTextureReference(&d_b_displ_ic_tex_ref_ptr, "d_b_displ_ic_tex"), 6021);
        print_CUDA_error_if_any(cudaBindTexture(0, d_b_displ_ic_tex_ref_ptr, mp->d_b_displ_inner_core,
                                                &channelDesc, sizeof(realw)*size), 6021);

        const textureReference* d_b_accel_ic_tex_ref_ptr;
        print_CUDA_error_if_any(cudaGetTextureReference(&d_b_accel_ic_tex_ref_ptr, "d_b_accel_ic_tex"), 6023);
        print_CUDA_error_if_any(cudaBindTexture(0, d_b_accel_ic_tex_ref_ptr, mp->d_b_accel_inner_core,
                                                &channelDesc, sizeof(realw)*size), 6023);

      }
#else
      cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<realw>();
      print_CUDA_error_if_any(cudaBindTexture(0, &d_displ_ic_tex, mp->d_displ_inner_core.cuda,
                                              &channelDesc, sizeof(realw)*size), 6021);
      print_CUDA_error_if_any(cudaBindTexture(0, &d_accel_ic_tex, mp->d_accel_inner_core.cuda,
                                              &channelDesc, sizeof(realw)*size), 6023);
      // backward/reconstructed wavefields
      if (mp->simulation_type == 3) {
        print_CUDA_error_if_any(cudaBindTexture(0, &d_b_displ_ic_tex, mp->d_b_displ_inner_core.cuda,
                                                &channelDesc, sizeof(realw)*size), 6021);
        print_CUDA_error_if_any(cudaBindTexture(0, &d_b_accel_ic_tex, mp->d_b_accel_inner_core.cuda,
                                                &channelDesc, sizeof(realw)*size), 6023);
      }
#endif
    }
#endif
  }
#endif
  // mass matrix
  gpuCreateCopy_todevice_realw (&mp->d_rmassz_inner_core, h_rmassz, size_glob);
  if (mp->rotation && mp->exact_mass_matrix_for_rotation) {
    gpuCreateCopy_todevice_realw (&mp->d_rmassx_inner_core, h_rmassx, size_glob);
    gpuCreateCopy_todevice_realw (&mp->d_rmassy_inner_core, h_rmassy, size_glob);
  } else {
    mp->d_rmassx_inner_core = gpuTakeRef(mp->d_rmassz_inner_core);
    mp->d_rmassy_inner_core = gpuTakeRef(mp->d_rmassz_inner_core);
  }

  // kernel simulations
  if (mp->simulation_type == 3) {
    // mass matrices
    mp->d_b_rmassz_inner_core = gpuTakeRef(mp->d_rmassz_inner_core);
    if (mp->rotation && mp->exact_mass_matrix_for_rotation) {
      gpuCreateCopy_todevice_realw (&mp->d_b_rmassx_inner_core, h_b_rmassx, size_glob);
      gpuCreateCopy_todevice_realw (&mp->d_b_rmassy_inner_core, h_b_rmassy, size_glob);
    } else {
      mp->d_b_rmassx_inner_core = gpuTakeRef(mp->d_rmassx_inner_core);
      mp->d_b_rmassy_inner_core = gpuTakeRef(mp->d_rmassy_inner_core);
    }


    // kernels
    size = NGLL3 * (mp->NSPEC_INNER_CORE);

    // density kernel
    gpuMalloc_realw (&mp->d_rho_kl_inner_core, size);
    gpuMemset_realw (&mp->d_rho_kl_inner_core, size, 0);

    // isotropic kernel
    gpuMalloc_realw (&mp->d_alpha_kl_inner_core, size);
    gpuMemset_realw (&mp->d_alpha_kl_inner_core, size, 0);

    gpuMalloc_realw (&mp->d_beta_kl_inner_core, size);
    gpuMemset_realw (&mp->d_beta_kl_inner_core, size, 0);
  }

  // mesh coloring
  mp->num_colors_outer_inner_core = *num_colors_outer;
  mp->num_colors_inner_inner_core = *num_colors_inner;
  mp->h_num_elem_colors_inner_core = (int *) num_elem_colors;

  // executes gpu calls
  // (especially needed for OpenCL kernels to finish executing)
  gpuSynchronize();

  GPU_ERROR_CHECKING ("prepare_inner_core_device");
}

/*----------------------------------------------------------------------------------------------- */
// cleanup
/*----------------------------------------------------------------------------------------------- */

extern EXTERN_LANG
void FC_FUNC_ (prepare_cleanup_device,
               PREPARE_CLEANUP_DEVICE) (long *Mesh_pointer_f,
                                        int *NCHUNKS_VAL) {

  TRACE ("prepare_cleanup_device");

  // frees allocated memory arrays
  Mesh *mp = (Mesh *) *Mesh_pointer_f;

  // synchronizes device
  gpuSynchronize();

  // frees memory on GPU
  //------------------------------------------
  // textures
  //------------------------------------------
#ifdef USE_CUDA
  if (run_cuda) {
#ifdef USE_TEXTURES_CONSTANTS
    cudaUnbindTexture(d_hprime_xx_tex);
    cudaUnbindTexture(d_hprimewgll_xx_tex);
#endif
#ifdef USE_TEXTURES_FIELDS
    cudaUnbindTexture(d_displ_cm_tex);
    cudaUnbindTexture(d_accel_cm_tex);
    if (mp->simulation_type == 3) {
      cudaUnbindTexture(d_b_displ_cm_tex);
      cudaUnbindTexture(d_b_accel_cm_tex);
    }
    cudaUnbindTexture(d_displ_ic_tex);
    cudaUnbindTexture(d_accel_ic_tex);
    if (mp->simulation_type == 3) {
      cudaUnbindTexture(d_b_displ_ic_tex);
      cudaUnbindTexture(d_b_accel_ic_tex);
    }
    cudaUnbindTexture(d_displ_oc_tex);
    cudaUnbindTexture(d_accel_oc_tex);
    if (mp->simulation_type == 3) {
      cudaUnbindTexture(d_b_displ_oc_tex);
      cudaUnbindTexture(d_b_accel_oc_tex);
    }
#endif
  }
#endif

  //------------------------------------------
  // pinned memory
  //------------------------------------------
  // receiver seismograms
  if (mp->nrec_local > 0) {
    if (GPU_ASYNC_COPY) {
#ifdef USE_OPENCL
      if (run_opencl) RELEASE_PINNED_BUFFER_OCL (station_seismo_field);
#endif
#ifdef USE_CUDA
      if (run_cuda) cudaFreeHost(mp->h_station_seismo_field);
#endif
    } else {
      free (mp->h_station_seismo_field);
    }
  }

  // adjoint source arrays
  if (mp->nadj_rec_local > 0) {
    if (GPU_ASYNC_COPY) {
#ifdef USE_OPENCL
      if (run_opencl) RELEASE_PINNED_BUFFER_OCL (adj_sourcearrays_slice);
#endif
#ifdef USE_CUDA
      if (run_cuda) cudaFreeHost(mp->h_adj_sourcearrays_slice);
#endif
    } else {
      free (mp->h_adj_sourcearrays_slice);
    }
  }

  // mpi buffers
#ifdef USE_OPENCL
  if (run_opencl) {
    if (mp->num_interfaces_crust_mantle > 0) {
      if (GPU_ASYNC_COPY) {
        RELEASE_PINNED_BUFFER_OCL (send_accel_buffer_cm);
        RELEASE_PINNED_BUFFER_OCL (recv_accel_buffer_cm);

        if (mp->simulation_type == 3) {
          RELEASE_PINNED_BUFFER_OCL (b_send_accel_buffer_cm);
          RELEASE_PINNED_BUFFER_OCL (b_recv_accel_buffer_cm);
        }
      }
    }
    if (mp->num_interfaces_inner_core > 0) {
      if (GPU_ASYNC_COPY) {
        RELEASE_PINNED_BUFFER_OCL (send_accel_buffer_ic);
        RELEASE_PINNED_BUFFER_OCL (recv_accel_buffer_ic);

        if (mp->simulation_type == 3) {
          RELEASE_PINNED_BUFFER_OCL (b_send_accel_buffer_ic);
          RELEASE_PINNED_BUFFER_OCL (b_recv_accel_buffer_ic);
        }
      }
    }
    if (mp->num_interfaces_outer_core > 0) {
      if (GPU_ASYNC_COPY) {
        RELEASE_PINNED_BUFFER_OCL (send_accel_buffer_oc);
        RELEASE_PINNED_BUFFER_OCL (recv_accel_buffer_oc);
        if (mp->simulation_type == 3) {
          RELEASE_PINNED_BUFFER_OCL (b_send_accel_buffer_oc);
          RELEASE_PINNED_BUFFER_OCL (b_recv_accel_buffer_oc);
        }
      }
    }
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    if (mp->num_interfaces_crust_mantle > 0) {
      if (GPU_ASYNC_COPY) {
        cudaFreeHost(mp->h_send_accel_buffer_cm);
        cudaFreeHost(mp->h_recv_accel_buffer_cm);
        if (mp->simulation_type == 3) {
          cudaFreeHost(mp->h_b_send_accel_buffer_cm);
          cudaFreeHost(mp->h_b_recv_accel_buffer_cm);
        }
      }
    }
    if (mp->num_interfaces_inner_core > 0) {
      if (GPU_ASYNC_COPY) {
        cudaFreeHost(mp->h_send_accel_buffer_ic);
        cudaFreeHost(mp->h_recv_accel_buffer_ic);
        if (mp->simulation_type == 3) {
          cudaFreeHost(mp->h_b_send_accel_buffer_ic);
          cudaFreeHost(mp->h_b_recv_accel_buffer_ic);
        }
      }
    }
    if (mp->num_interfaces_outer_core > 0) {
      if (GPU_ASYNC_COPY) {
        cudaFreeHost(mp->h_send_accel_buffer_oc);
        cudaFreeHost(mp->h_recv_accel_buffer_oc);
        if (mp->simulation_type == 3) {
          cudaFreeHost(mp->h_b_send_accel_buffer_oc);
          cudaFreeHost(mp->h_b_recv_accel_buffer_oc);
        }
      }
    }
  }
#endif

  //------------------------------------------
  // constants
  //------------------------------------------
#ifdef USE_OPENCL
  if (run_opencl) {
#ifdef USE_TEXTURES_CONSTANTS
    clReleaseMemObject (mp->d_hprime_xx.ocl);
    clReleaseMemObject (mp->d_hprimewgll_xx.ocl);
#endif

    clReleaseMemObject (mp->d_wgllwgll_xy.ocl);
    clReleaseMemObject (mp->d_wgllwgll_xz.ocl);
    clReleaseMemObject (mp->d_wgllwgll_yz.ocl);
    clReleaseMemObject (mp->d_wgll_cube.ocl);
  }
#endif

  //------------------------------------------
  // sources
  //------------------------------------------
  if (mp->simulation_type == 1 || mp->simulation_type == 3) {
    gpuFree (&mp->d_sourcearrays);
    gpuFree (&mp->d_stf_pre_compute);
  }
  gpuFree (&mp->d_islice_selected_source);
  gpuFree (&mp->d_ispec_selected_source);

  //------------------------------------------
  // receivers
  //------------------------------------------
  if (mp->nrec_local > 0) {
    gpuFree (&mp->d_number_receiver_global);
    gpuFree (&mp->d_station_seismo_field);
    if (mp->simulation_type == 2) {
      gpuFree (&mp->d_station_strain_field);
    }
  }
  gpuFree (&mp->d_ispec_selected_rec);
  if (mp->nadj_rec_local > 0) {
    gpuFree (&mp->d_adj_sourcearrays);
    gpuFree (&mp->d_pre_computed_irec);
  }

  gpuFree (&mp->d_norm_max);

  //------------------------------------------
  // rotation arrays
  //------------------------------------------
  if (mp->rotation) {
    gpuFree (&mp->d_A_array_rotation);
    gpuFree (&mp->d_B_array_rotation);
    if (mp->simulation_type == 3) {
      gpuFree (&mp->d_b_A_array_rotation);
      gpuFree (&mp->d_b_B_array_rotation);
    }
  }

  //------------------------------------------
  // gravity arrays
  //------------------------------------------
  if (! mp->gravity) {
    gpuFree (&mp->d_d_ln_density_dr_table);
  } else {
    gpuFree (&mp->d_minus_rho_g_over_kappa_fluid);
    gpuFree (&mp->d_minus_gravity_table);
    gpuFree (&mp->d_minus_deriv_gravity_table);
    gpuFree (&mp->d_density_table);
  }

  //------------------------------------------
  // attenuation arrays
  //------------------------------------------
  if (mp->attenuation) {
    gpuFree (&mp->d_one_minus_sum_beta_crust_mantle);
    gpuFree (&mp->d_one_minus_sum_beta_inner_core);
    if (! mp->partial_phys_dispersion_only) {
      gpuFree (&mp->d_factor_common_crust_mantle);
      gpuFree (&mp->d_R_xx_crust_mantle);
      gpuFree (&mp->d_R_yy_crust_mantle);
      gpuFree (&mp->d_R_xy_crust_mantle);
      gpuFree (&mp->d_R_xz_crust_mantle);
      gpuFree (&mp->d_R_yz_crust_mantle);
      gpuFree (&mp->d_factor_common_inner_core);
      gpuFree (&mp->d_R_xx_inner_core);
      gpuFree (&mp->d_R_yy_inner_core);
      gpuFree (&mp->d_R_xy_inner_core);
      gpuFree (&mp->d_R_xz_inner_core);
      gpuFree (&mp->d_R_yz_inner_core);
    }
    gpuFree (&mp->d_alphaval);
    gpuFree (&mp->d_betaval);
    gpuFree (&mp->d_gammaval);
    if (mp->simulation_type == 3) {
      gpuFree (&mp->d_b_alphaval);
      gpuFree (&mp->d_b_betaval);
      gpuFree (&mp->d_b_gammaval);
    }
  }

  //------------------------------------------
  // strain
  //------------------------------------------
  if (mp->compute_and_store_strain) {
    gpuFree (&mp->d_epsilondev_xx_crust_mantle);
    gpuFree (&mp->d_epsilondev_yy_crust_mantle);
    gpuFree (&mp->d_epsilondev_xy_crust_mantle);
    gpuFree (&mp->d_epsilondev_xz_crust_mantle);
    gpuFree (&mp->d_epsilondev_yz_crust_mantle);

    gpuFree (&mp->d_epsilondev_xx_inner_core);
    gpuFree (&mp->d_epsilondev_yy_inner_core);
    gpuFree (&mp->d_epsilondev_xy_inner_core);
    gpuFree (&mp->d_epsilondev_xz_inner_core);
    gpuFree (&mp->d_epsilondev_yz_inner_core);

    gpuFree (&mp->d_eps_trace_over_3_crust_mantle);
    gpuFree (&mp->d_eps_trace_over_3_inner_core);

    if (mp->simulation_type == 3 && ! mp->undo_attenuation) {
      gpuFree (&mp->d_b_epsilondev_xx_crust_mantle);
      gpuFree (&mp->d_b_epsilondev_yy_crust_mantle);
      gpuFree (&mp->d_b_epsilondev_xy_crust_mantle);
      gpuFree (&mp->d_b_epsilondev_xz_crust_mantle);
      gpuFree (&mp->d_b_epsilondev_yz_crust_mantle);

      gpuFree (&mp->d_b_epsilondev_xx_inner_core);
      gpuFree (&mp->d_b_epsilondev_yy_inner_core);
      gpuFree (&mp->d_b_epsilondev_xy_inner_core);
      gpuFree (&mp->d_b_epsilondev_xz_inner_core);
      gpuFree (&mp->d_b_epsilondev_yz_inner_core);

      gpuFree (&mp->d_b_eps_trace_over_3_crust_mantle);
      gpuFree (&mp->d_b_eps_trace_over_3_inner_core);
    }
  }

  //------------------------------------------
  // absorbing boundaries arrays
  //------------------------------------------
  if (mp->absorbing_conditions) {
    gpuFree (&mp->d_rho_vp_crust_mantle);
    gpuFree (&mp->d_rho_vs_crust_mantle);
    gpuFree (&mp->d_nkmin_xi_crust_mantle);
    gpuFree (&mp->d_nkmin_eta_crust_mantle);
    gpuFree (&mp->d_njmin_crust_mantle);
    gpuFree (&mp->d_njmax_crust_mantle);
    gpuFree (&mp->d_nimin_crust_mantle);
    gpuFree (&mp->d_nimax_crust_mantle);
    if (mp->nspec2D_xmin_crust_mantle > 0) {
      gpuFree (&mp->d_ibelm_xmin_crust_mantle);
      gpuFree (&mp->d_normal_xmin_crust_mantle);
      gpuFree (&mp->d_jacobian2D_xmin_crust_mantle);
      if ((mp->simulation_type == 1 && mp->save_forward) || (mp->simulation_type == 3)) {
        gpuFree (&mp->d_absorb_xmin_crust_mantle);
      }
    }
    if (mp->nspec2D_xmax_crust_mantle > 0) {
      gpuFree (&mp->d_ibelm_xmax_crust_mantle);
      gpuFree (&mp->d_normal_xmax_crust_mantle);
      gpuFree (&mp->d_jacobian2D_xmax_crust_mantle);
      if ((mp->simulation_type == 1 && mp->save_forward) || (mp->simulation_type == 3)) {
        gpuFree (&mp->d_absorb_xmax_crust_mantle);
      }
    }
    if (mp->nspec2D_ymin_crust_mantle > 0) {
      gpuFree (&mp->d_ibelm_ymin_crust_mantle);
      gpuFree (&mp->d_normal_ymin_crust_mantle);
      gpuFree (&mp->d_jacobian2D_ymin_crust_mantle);
      if ((mp->simulation_type == 1 && mp->save_forward) || (mp->simulation_type == 3)) {
        gpuFree (&mp->d_absorb_ymin_crust_mantle);
      }
    }
    if (mp->nspec2D_ymax_crust_mantle > 0) {
      gpuFree (&mp->d_ibelm_ymax_crust_mantle);
      gpuFree (&mp->d_normal_ymax_crust_mantle);
      gpuFree (&mp->d_jacobian2D_ymax_crust_mantle);
      if ((mp->simulation_type == 1 && mp->save_forward) || (mp->simulation_type == 3)) {
        gpuFree (&mp->d_absorb_ymax_crust_mantle);
      }
    }
    gpuFree (&mp->d_vp_outer_core);
    gpuFree (&mp->d_nkmin_xi_outer_core);
    gpuFree (&mp->d_nkmin_eta_outer_core);
    gpuFree (&mp->d_njmin_outer_core);
    gpuFree (&mp->d_njmax_outer_core);
    gpuFree (&mp->d_nimin_outer_core);
    gpuFree (&mp->d_nimax_outer_core);
    if (mp->nspec2D_xmin_outer_core > 0) {
      gpuFree (&mp->d_ibelm_xmin_outer_core);
      gpuFree (&mp->d_jacobian2D_xmin_outer_core);
      if ((mp->simulation_type == 1 && mp->save_forward) || (mp->simulation_type == 3)) {
        gpuFree (&mp->d_absorb_xmin_outer_core);
      }
    }
    if (mp->nspec2D_xmax_outer_core > 0) {
      gpuFree (&mp->d_ibelm_xmax_outer_core);
      gpuFree (&mp->d_jacobian2D_xmax_outer_core);
      if ((mp->simulation_type == 1 && mp->save_forward) || (mp->simulation_type == 3)) {
        gpuFree (&mp->d_absorb_xmax_outer_core);
      }
    }
    if (mp->nspec2D_ymin_outer_core > 0) {
      gpuFree (&mp->d_ibelm_ymin_outer_core);
      gpuFree (&mp->d_jacobian2D_ymin_outer_core);
      if ((mp->simulation_type == 1 && mp->save_forward) || (mp->simulation_type == 3)) {
        gpuFree (&mp->d_absorb_ymin_outer_core);
      }
    }
    if (mp->nspec2D_ymax_outer_core > 0) {
      gpuFree (&mp->d_ibelm_ymax_outer_core);
      gpuFree (&mp->d_jacobian2D_ymax_outer_core);
      if ((mp->simulation_type == 1 && mp->save_forward) || (mp->simulation_type == 3)) {
        gpuFree (&mp->d_absorb_ymax_outer_core);
      }
    }
    if (mp->nspec2D_zmin_outer_core > 0) {
      if ((mp->simulation_type == 1 && mp->save_forward) || (mp->simulation_type == 3)) {
        gpuFree (&mp->d_absorb_zmin_outer_core);
      }
    }
  }

  //------------------------------------------
  // MPI buffers
  //------------------------------------------
  if (mp->num_interfaces_crust_mantle > 0) {
    gpuFree (&mp->d_nibool_interfaces_crust_mantle);
    gpuFree (&mp->d_ibool_interfaces_crust_mantle);
    gpuFree (&mp->d_send_accel_buffer_crust_mantle);
    if (mp->simulation_type == 3) gpuFree (&mp->d_b_send_accel_buffer_crust_mantle);
  }
  if (mp->num_interfaces_inner_core > 0) {
    gpuFree (&mp->d_nibool_interfaces_inner_core);
    gpuFree (&mp->d_ibool_interfaces_inner_core);
    gpuFree (&mp->d_send_accel_buffer_inner_core);
    if (mp->simulation_type == 3) gpuFree (&mp->d_b_send_accel_buffer_inner_core);
  }
  if (mp->num_interfaces_outer_core > 0) {
    gpuFree (&mp->d_nibool_interfaces_outer_core);
    gpuFree (&mp->d_ibool_interfaces_outer_core);
    gpuFree (&mp->d_send_accel_buffer_outer_core);
    if (mp->simulation_type == 3) gpuFree (&mp->d_b_send_accel_buffer_outer_core);
  }

  //------------------------------------------
  // NOISE arrays
  //------------------------------------------
  if (mp->noise_tomography > 0) {
    gpuFree (&mp->d_ibelm_top_crust_mantle);
    gpuFree (&mp->d_noise_surface_movie);
    if (mp->noise_tomography == 1) gpuFree (&mp->d_noise_sourcearray);
    if (mp->noise_tomography > 1) {
      gpuFree (&mp->d_normal_x_noise);
      gpuFree (&mp->d_normal_y_noise);
      gpuFree (&mp->d_normal_z_noise);
      gpuFree (&mp->d_mask_noise);
      gpuFree (&mp->d_jacobian2D_top_crust_mantle);
    }
    if (mp->noise_tomography == 3) gpuFree (&mp->d_Sigma_kl);
  }

  //------------------------------------------
  // crust_mantle
  //------------------------------------------
  gpuFree (&mp->d_xix_crust_mantle);
  gpuFree (&mp->d_xiy_crust_mantle);
  gpuFree (&mp->d_xiz_crust_mantle);
  gpuFree (&mp->d_etax_crust_mantle);
  gpuFree (&mp->d_etay_crust_mantle);
  gpuFree (&mp->d_etaz_crust_mantle);
  gpuFree (&mp->d_gammax_crust_mantle);
  gpuFree (&mp->d_gammay_crust_mantle);
  gpuFree (&mp->d_gammaz_crust_mantle);

  gpuFree (&mp->d_ibool_crust_mantle);

  if (! mp->anisotropic_3D_mantle) {
    gpuFree (&mp->d_ispec_is_tiso_crust_mantle);
    gpuFree (&mp->d_kappavstore_crust_mantle);
    gpuFree (&mp->d_kappahstore_crust_mantle);
    gpuFree (&mp->d_muvstore_crust_mantle);
    gpuFree (&mp->d_muhstore_crust_mantle);
    gpuFree (&mp->d_eta_anisostore_crust_mantle);
  } else {
    gpuFree (&mp->d_c11store_crust_mantle);
    gpuFree (&mp->d_c12store_crust_mantle);
    gpuFree (&mp->d_c13store_crust_mantle);
    gpuFree (&mp->d_c14store_crust_mantle);
    gpuFree (&mp->d_c15store_crust_mantle);
    gpuFree (&mp->d_c16store_crust_mantle);
    gpuFree (&mp->d_c22store_crust_mantle);
    gpuFree (&mp->d_c23store_crust_mantle);
    gpuFree (&mp->d_c24store_crust_mantle);
    gpuFree (&mp->d_c25store_crust_mantle);
    gpuFree (&mp->d_c26store_crust_mantle);
    gpuFree (&mp->d_c33store_crust_mantle);
    gpuFree (&mp->d_c34store_crust_mantle);
    gpuFree (&mp->d_c35store_crust_mantle);
    gpuFree (&mp->d_c36store_crust_mantle);
    gpuFree (&mp->d_c44store_crust_mantle);
    gpuFree (&mp->d_c45store_crust_mantle);
    gpuFree (&mp->d_c46store_crust_mantle);
    gpuFree (&mp->d_c55store_crust_mantle);
    gpuFree (&mp->d_c56store_crust_mantle);
    gpuFree (&mp->d_c66store_crust_mantle);
  }

  if (mp->simulation_type == 3 && mp->save_boundary_mesh) {
    gpuFree (&mp->d_rhostore_crust_mantle);
  }

  gpuFree (&mp->d_ystore_crust_mantle);
  gpuFree (&mp->d_zstore_crust_mantle);
  if (mp->gravity) gpuFree (&mp->d_xstore_crust_mantle);

  gpuFree (&mp->d_phase_ispec_inner_crust_mantle);
  gpuFree (&mp->d_ibelm_bottom_crust_mantle);

  // wavefield
  gpuFree (&mp->d_displ_crust_mantle);
  gpuFree (&mp->d_veloc_crust_mantle);
  gpuFree (&mp->d_accel_crust_mantle);
  if (mp->simulation_type == 3) {
    gpuFree (&mp->d_b_displ_crust_mantle);
    gpuFree (&mp->d_b_veloc_crust_mantle);
    gpuFree (&mp->d_b_accel_crust_mantle);
  }

  // mass matrix
  gpuFree (&mp->d_rmassz_crust_mantle);
  if (*NCHUNKS_VAL != 6 && mp->absorbing_conditions) {
    gpuFree (&mp->d_rmassx_crust_mantle);
    gpuFree (&mp->d_rmassy_crust_mantle);
  }

  // kernel simulations
  if (mp->simulation_type == 3) {
    if (mp->rotation && mp->exact_mass_matrix_for_rotation) {
      gpuFree (&mp->d_b_rmassx_crust_mantle);
      gpuFree (&mp->d_b_rmassy_crust_mantle);
    }
    // kernels
    gpuFree (&mp->d_rho_kl_crust_mantle);
    if (!mp->anisotropic_kl) {
      gpuFree (&mp->d_alpha_kl_crust_mantle);
      gpuFree (&mp->d_beta_kl_crust_mantle);
    } else {
      gpuFree (&mp->d_cijkl_kl_crust_mantle);
    }
    if (mp->approximate_hess_kl) {
      gpuFree (&mp->d_hess_kl_crust_mantle);
    }
  }

  //------------------------------------------
  // outer_core
  //------------------------------------------
  gpuFree (&mp->d_xix_outer_core);
  gpuFree (&mp->d_xiy_outer_core);
  gpuFree (&mp->d_xiz_outer_core);
  gpuFree (&mp->d_etax_outer_core);
  gpuFree (&mp->d_etay_outer_core);
  gpuFree (&mp->d_etaz_outer_core);
  gpuFree (&mp->d_gammax_outer_core);
  gpuFree (&mp->d_gammay_outer_core);
  gpuFree (&mp->d_gammaz_outer_core);

  gpuFree (&mp->d_kappavstore_outer_core);
  if (mp->simulation_type == 3) {
    gpuFree (&mp->d_rhostore_outer_core);
  }

  gpuFree (&mp->d_ibool_outer_core);

  gpuFree (&mp->d_xstore_outer_core);
  gpuFree (&mp->d_ystore_outer_core);
  gpuFree (&mp->d_zstore_outer_core);

  gpuFree (&mp->d_phase_ispec_inner_outer_core);

  gpuFree (&mp->d_ibelm_top_outer_core);
  gpuFree (&mp->d_jacobian2D_top_outer_core);
  gpuFree (&mp->d_normal_top_outer_core);

  gpuFree (&mp->d_ibelm_bottom_outer_core);
  gpuFree (&mp->d_normal_bottom_outer_core);
  gpuFree (&mp->d_jacobian2D_bottom_outer_core);

  // wavefield
  gpuFree (&mp->d_displ_outer_core);
  gpuFree (&mp->d_veloc_outer_core);
  gpuFree (&mp->d_accel_outer_core);
  if (mp->simulation_type == 3) {
    gpuFree (&mp->d_b_displ_outer_core);
    gpuFree (&mp->d_b_veloc_outer_core);
    gpuFree (&mp->d_b_accel_outer_core);
  }

  // mass matrix
  gpuFree (&mp->d_rmass_outer_core);

  if (mp->simulation_type == 3) {
    gpuFree (&mp->d_rho_kl_outer_core);
    gpuFree (&mp->d_alpha_kl_outer_core);
  }

  //------------------------------------------
  // inner_core
  //------------------------------------------
  gpuFree (&mp->d_xix_inner_core);
  gpuFree (&mp->d_xiy_inner_core);
  gpuFree (&mp->d_xiz_inner_core);
  gpuFree (&mp->d_etax_inner_core);
  gpuFree (&mp->d_etay_inner_core);
  gpuFree (&mp->d_etaz_inner_core);
  gpuFree (&mp->d_gammax_inner_core);
  gpuFree (&mp->d_gammay_inner_core);
  gpuFree (&mp->d_gammaz_inner_core);

  gpuFree (&mp->d_muvstore_inner_core);
  if (! mp->anisotropic_inner_core) {
    gpuFree (&mp->d_kappavstore_inner_core);
  } else {
    gpuFree (&mp->d_c11store_inner_core);
    gpuFree (&mp->d_c12store_inner_core);
    gpuFree (&mp->d_c13store_inner_core);
    gpuFree (&mp->d_c33store_inner_core);
    gpuFree (&mp->d_c44store_inner_core);
  }

  if (mp->simulation_type == 3 && mp->save_boundary_mesh) {
    gpuFree (&mp->d_rhostore_inner_core);
  }

  gpuFree (&mp->d_ibool_inner_core);
  gpuFree (&mp->d_idoubling_inner_core);

  // gravity
  if (mp->gravity) {
    gpuFree (&mp->d_xstore_inner_core);
    gpuFree (&mp->d_ystore_inner_core);
    gpuFree (&mp->d_zstore_inner_core);
  }

  gpuFree (&mp->d_phase_ispec_inner_inner_core);
  gpuFree (&mp->d_ibelm_top_inner_core);

  // wavefield
  gpuFree (&mp->d_displ_inner_core);
  gpuFree (&mp->d_veloc_inner_core);
  gpuFree (&mp->d_accel_inner_core);
  if (mp->simulation_type == 3) {
    gpuFree (&mp->d_b_displ_inner_core);
    gpuFree (&mp->d_b_veloc_inner_core);
    gpuFree (&mp->d_b_accel_inner_core);
  }

  // mass matrix
  gpuFree (&mp->d_rmassz_inner_core);
  if (mp->rotation && mp->exact_mass_matrix_for_rotation) {
    gpuFree (&mp->d_rmassx_inner_core);
    gpuFree (&mp->d_rmassy_inner_core);
  }

  // kernel simulations
  if (mp->simulation_type == 3) {
    if (mp->rotation && mp->exact_mass_matrix_for_rotation) {
      gpuFree (&mp->d_b_rmassx_inner_core);
      gpuFree (&mp->d_b_rmassy_inner_core);
    }
    // kernels
    gpuFree (&mp->d_rho_kl_inner_core);
    gpuFree (&mp->d_alpha_kl_inner_core);
    gpuFree (&mp->d_beta_kl_inner_core);
  }

  //------------------------------------------
  // oceans
  //------------------------------------------
  if (mp->oceans) {
    gpuFree (&mp->d_ibool_ocean_load);
    gpuFree (&mp->d_rmass_ocean_load);
    gpuFree (&mp->d_normal_ocean_load);
  }

#ifdef USE_OPENCL
  if (run_opencl) {
#ifdef USE_TEXTURES_FIELDS
    clReleaseMemObject (mp->d_displ_cm_tex);
    clReleaseMemObject (mp->d_accel_cm_tex);
    clReleaseMemObject (mp->d_b_displ_cm_tex);
    clReleaseMemObject (mp->d_b_accel_cm_tex);

    clReleaseMemObject (mp->d_displ_oc_tex);
    clReleaseMemObject (mp->d_accel_oc_tex);
    clReleaseMemObject (mp->d_b_displ_oc_tex);
    clReleaseMemObject (mp->d_b_accel_oc_tex);

    clReleaseMemObject (mp->d_displ_ic_tex);
    clReleaseMemObject (mp->d_accel_ic_tex);
    clReleaseMemObject (mp->d_b_displ_ic_tex);
    clReleaseMemObject (mp->d_b_accel_ic_tex);
#endif
#ifdef USE_TEXTURES_CONSTANTS
    clReleaseMemObject (mp->d_hprime_xx_cm_tex);
    clReleaseMemObject (mp->d_hprimewgll_xx_cm_tex);
#endif
  }
#endif

  // synchronizes device
  gpuSynchronize();

  // cleans up asynchronous streams/queues
#ifdef USE_OPENCL
  if (run_opencl) {
    // cleans up queues
    clReleaseCommandQueue (mocl.command_queue);
    if (GPU_ASYNC_COPY) clReleaseCommandQueue (mocl.copy_queue);
  }
#endif
#ifdef USE_CUDA
  if (run_cuda) {
    cudaStreamDestroy(mp->compute_stream);
    if (GPU_ASYNC_COPY) cudaStreamDestroy(mp->copy_stream);
  }
#endif

  // specific OpenCL: frees kernels and programs
#ifdef USE_OPENCL
  if (run_opencl) release_kernels();
#endif

  // releases previous contexts
  gpuReset();

  // mesh pointer - not needed anymore
  free (mp);
}
