
  allocate(rmass_ocean_load(NGLOB_CRUST_MANTLE_OCEANS),stat=ier)
  if(ier /= 0) then
    print *,"ABORTING can not allocate in allocate_after_2 ier=",ier
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

!---

  allocate(displ_crust_mantle(NDIM,NGLOB_CRUST_MANTLE),stat=ier)
  if(ier /= 0) then
    print *,"ABORTING can not allocate in allocate_after_2 ier=",ier
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(veloc_crust_mantle(NDIM,NGLOB_CRUST_MANTLE),stat=ier)
  if(ier /= 0) then
    print *,"ABORTING can not allocate in allocate_after_2 ier=",ier
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(accel_crust_mantle(NDIM,NGLOB_CRUST_MANTLE),stat=ier)
  if(ier /= 0) then
    print *,"ABORTING can not allocate in allocate_after_2 ier=",ier
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif


  allocate(displ_outer_core(NGLOB_OUTER_CORE),stat=ier)
  if(ier /= 0) then
    print *,"ABORTING can not allocate in allocate_after_2 ier=",ier
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(veloc_outer_core(NGLOB_OUTER_CORE),stat=ier)
  if(ier /= 0) then
    print *,"ABORTING can not allocate in allocate_after_2 ier=",ier
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(accel_outer_core(NGLOB_OUTER_CORE),stat=ier)
  if(ier /= 0) then
    print *,"ABORTING can not allocate in allocate_after_2 ier=",ier
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif


  allocate(displ_inner_core(NDIM,NGLOB_INNER_CORE),stat=ier)
  if(ier /= 0) then
    print *,"ABORTING can not allocate in allocate_after_2 ier=",ier
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(veloc_inner_core(NDIM,NGLOB_INNER_CORE),stat=ier)
  if(ier /= 0) then
    print *,"ABORTING can not allocate in allocate_after_2 ier=",ier
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif

  allocate(accel_inner_core(NDIM,NGLOB_INNER_CORE),stat=ier)
  if(ier /= 0) then
    print *,"ABORTING can not allocate in allocate_after_2 ier=",ier
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif


!---

  allocate(R_memory_crust_mantle(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT),stat=ier)
  if(ier /= 0) then
    print *,"ABORTING can not allocate in allocate_after_2 ier=",ier
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif


  allocate(R_memory_inner_core(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ATTENUATION),stat=ier)
  if(ier /= 0) then
    print *,"ABORTING can not allocate in allocate_after_2 ier=",ier
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif


!---

  allocate(epsilondev_crust_mantle(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_OR_ATT),stat=ier)
  if(ier /= 0) then
    print *,"ABORTING can not allocate in allocate_after_2 ier=",ier
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif


  allocate(epsilondev_inner_core(5,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STR_OR_ATT),stat=ier)
  if(ier /= 0) then
    print *,"ABORTING can not allocate in allocate_after_2 ier=",ier
    call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
  endif


