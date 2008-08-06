
!! DK DK this for the new merged version

  deallocate(xstore_crust_mantle,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(ystore_crust_mantle,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(zstore_crust_mantle,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif


  deallocate(xstore_outer_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(ystore_outer_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(zstore_outer_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif


  deallocate(xstore_inner_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(ystore_inner_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(zstore_inner_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif


!---

  deallocate(xix_crust_mantle,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(xiy_crust_mantle,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(xiz_crust_mantle,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(etax_crust_mantle,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(etay_crust_mantle,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(etaz_crust_mantle,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(gammax_crust_mantle,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(gammay_crust_mantle,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(gammaz_crust_mantle,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif


  deallocate(xix_outer_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(xiy_outer_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(xiz_outer_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(etax_outer_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(etay_outer_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(etaz_outer_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(gammax_outer_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(gammay_outer_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(gammaz_outer_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif


  deallocate(xix_inner_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(xiy_inner_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(xiz_inner_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(etax_inner_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(etay_inner_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(etaz_inner_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(gammax_inner_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(gammay_inner_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

  deallocate(gammaz_inner_core,STAT=ier )
  if (ier /= 0) then
    print *,"ERROR can not deallocate in deallocate.f90 ier=",ier
  endif

