
!! DK DK created this for merged version

  call meshfem3D( &
  myrank,sizeprocs,addressing,ichunk_slice,iproc_xi_slice,iproc_eta_slice,ibathy_topo,NSOURCES, &
  NTSTEP_BETWEEN_OUTPUT_SEISMOS,ibool_crust_mantle,ibool_outer_core,ibool_inner_core, &
  idoubling_crust_mantle,idoubling_inner_core,ibelm_bottom_crust_mantle, ibelm_bottom_outer_core, ibelm_top_outer_core, &
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
  npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core)

