
  allocate(xstore_crust_mantle(NGLOB_CRUST_MANTLE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(ystore_crust_mantle(NGLOB_CRUST_MANTLE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(zstore_crust_mantle(NGLOB_CRUST_MANTLE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif


  allocate(xstore_outer_core(NGLOB_OUTER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     
    
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(ystore_outer_core(NGLOB_OUTER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(zstore_outer_core(NGLOB_OUTER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif


  allocate(xstore_inner_core(NGLOB_INNER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(ystore_inner_core(NGLOB_INNER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(zstore_inner_core(NGLOB_INNER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif


!---

  allocate(xix_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(xiy_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(xiz_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(etax_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(etay_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(etaz_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(gammax_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(gammay_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(gammaz_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif


  allocate(xix_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(xiy_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(xiz_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(etax_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(etay_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(etaz_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(gammax_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(gammay_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(gammaz_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif


  allocate(xix_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(xiy_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(xiz_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(etax_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(etay_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(etaz_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(gammax_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(gammay_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(gammaz_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
  if(ier /= 0) then 
     print *,"ABORTING can not allocate in allocate_after_1 ier=",ier
     call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif


