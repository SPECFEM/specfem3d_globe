
!! DK DK created this for merged version

  subroutine specfem3D( &
!! DK DK to do later, for attenuation only; not done yet by lack of time
  omsb_crust_mantle_dble,factor_scale_crust_mantle_dble, omsb_inner_core_dble,factor_scale_inner_core_dble, &
  one_minus_sum_beta_crust_mantle,factor_scale_crust_mantle, one_minus_sum_beta_inner_core,factor_scale_inner_core, &
  factor_common_crust_mantle,factor_common_inner_core,factor_common_crust_mantle_dble, factor_common_inner_core_dble, &
!! DK DK to do later, for oceans only
  rmass_ocean_load, &
!! DK DK already computed
  myrank,sizeprocs,addressing,ichunk_slice,iproc_xi_slice,iproc_eta_slice,ibathy_topo, &
  ibool_crust_mantle, ibool_outer_core, ibool_inner_core, idoubling_crust_mantle,idoubling_inner_core, &
ibelm_bottom_crust_mantle, ibelm_bottom_outer_core, ibelm_top_outer_core, &
ibelm_xmin_inner_core,ibelm_xmax_inner_core,ibelm_ymin_inner_core,ibelm_ymax_inner_core,ibelm_bottom_inner_core, &
ibelm_top_inner_core,iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle, iboolleft_eta_crust_mantle, &
iboolright_eta_crust_mantle,iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
  iboolleft_xi_inner_core,iboolright_xi_inner_core, iboolleft_eta_inner_core,iboolright_eta_inner_core,&
  jacobian2D_bottom_outer_core,jacobian2D_top_outer_core, &
  normal_bottom_outer_core, normal_top_outer_core,kappavstore_crust_mantle,muvstore_crust_mantle, &
  kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle,kappavstore_inner_core,muvstore_inner_core, &
  rmass_crust_mantle, rmass_outer_core, rmass_inner_core, &
  nspec2D_xmin_inner_core,nspec2D_xmax_inner_core,nspec2D_ymin_inner_core,nspec2D_ymax_inner_core, &
iprocfrom_faces,iprocto_faces,imsg_type,iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
  iboolfaces_crust_mantle,iboolfaces_outer_core,iboolfaces_inner_core, &
  iboolcorner_crust_mantle,iboolcorner_outer_core,iboolcorner_inner_core, &
  npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
  npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
  npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
!! DK DK recomputed after the end of the mesher and before the beginning of the solver
  xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle,xstore_outer_core,ystore_outer_core,zstore_outer_core, &
xstore_inner_core,ystore_inner_core,zstore_inner_core, &
xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle,etax_crust_mantle,etay_crust_mantle, &
  etaz_crust_mantle,gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
xix_outer_core,xiy_outer_core,xiz_outer_core,etax_outer_core,etay_outer_core, &
etaz_outer_core,gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
  xix_inner_core,xiy_inner_core,xiz_inner_core,etax_inner_core,etay_inner_core, &
  etaz_inner_core,gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
!! DK DK do not need to be initialized
  displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle, &
  displ_outer_core,veloc_outer_core,accel_outer_core,displ_inner_core,veloc_inner_core,accel_inner_core, &
  buffer_send_chunkcorners_scalar,buffer_recv_chunkcorners_scalar, &
  buffer_send_chunkcorners_vector,buffer_recv_chunkcorners_vector, &
R_memory_crust_mantle, epsilondev_crust_mantle, R_memory_inner_core, epsilondev_inner_core)

