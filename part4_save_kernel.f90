
  ! synchronize all processes
  call MPI_BARRIER(MPI_COMM_WORLD,ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error synchronize saving forward')

  ! dump kernel arrays
  if (SIMULATION_TYPE == 3) then

    ! crust mantle
    if (SAVE_REGULAR_KL) then
    call save_regular_kernels_crust_mantle(myrank, &
                  npoints_slice, hxir_reg, hetar_reg, hgammar_reg, &
                  scale_t,scale_displ, &
                  cijkl_kl_crust_mantle,rho_kl_crust_mantle, &
                  alpha_kl_crust_mantle,beta_kl_crust_mantle, &
                  ystore_crust_mantle,zstore_crust_mantle, &
                  rhostore_crust_mantle,muvstore_crust_mantle, &
                  kappavstore_crust_mantle,ibool_crust_mantle, &
                  kappahstore_crust_mantle,muhstore_crust_mantle, &
                  eta_anisostore_crust_mantle,ispec_is_tiso_crust_mantle,LOCAL_PATH)
    else
    call save_kernels_crust_mantle(myrank,scale_t,scale_displ, &
                  cijkl_kl_crust_mantle,rho_kl_crust_mantle, &
                  alpha_kl_crust_mantle,beta_kl_crust_mantle, &
                  ystore_crust_mantle,zstore_crust_mantle, &
                  rhostore_crust_mantle,muvstore_crust_mantle, &
                  kappavstore_crust_mantle,ibool_crust_mantle, &
                  kappahstore_crust_mantle,muhstore_crust_mantle, &
                  eta_anisostore_crust_mantle,ispec_is_tiso_crust_mantle,LOCAL_PATH)
    endif

    ! noise strength kernel
    if (NOISE_TOMOGRAPHY == 3) then
       call save_kernels_strength_noise(myrank,LOCAL_PATH,Sigma_kl_crust_mantle)
    endif

    ! outer core
    call save_kernels_outer_core(myrank,scale_t,scale_displ, &
                        rho_kl_outer_core,alpha_kl_outer_core, &
                        rhostore_outer_core,kappavstore_outer_core, &
                        deviatoric_outercore,nspec_beta_kl_outer_core,beta_kl_outer_core,LOCAL_PATH)

    ! inner core
    call save_kernels_inner_core(myrank,scale_t,scale_displ, &
                          rho_kl_inner_core,beta_kl_inner_core,alpha_kl_inner_core, &
                          rhostore_inner_core,muvstore_inner_core,kappavstore_inner_core,LOCAL_PATH)

    ! boundary kernel
    if (SAVE_BOUNDARY_MESH) then
      call save_kernels_boundary_kl(myrank,scale_t,scale_displ, &
                                  moho_kl,d400_kl,d670_kl,cmb_kl,icb_kl,LOCAL_PATH,HONOR_1D_SPHERICAL_MOHO)
    endif

    ! approximate hessian
    if( APPROXIMATE_HESS_KL) then
      stop "APPROXIMATE_HESS_KL is not support for undoing att"  !ZN
      call save_kernels_hessian(myrank,scale_t,scale_displ,hess_kl_crust_mantle,LOCAL_PATH)
    endif
  endif

