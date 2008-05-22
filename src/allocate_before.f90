
!! DK DK added this for merged version
  allocate(xelm_store_crust_mantle(NGNOD,NSPEC_CRUST_MANTLE),stat=ier)
  if(ier /= 0) stop 'error memory allocation merged version'
  allocate(yelm_store_crust_mantle(NGNOD,NSPEC_CRUST_MANTLE),stat=ier)
  if(ier /= 0) stop 'error memory allocation merged version'
  allocate(zelm_store_crust_mantle(NGNOD,NSPEC_CRUST_MANTLE),stat=ier)
  if(ier /= 0) stop 'error memory allocation merged version'

  allocate(xelm_store_outer_core(NGNOD,NSPEC_OUTER_CORE),stat=ier)
  if(ier /= 0) stop 'error memory allocation merged version'
  allocate(yelm_store_outer_core(NGNOD,NSPEC_OUTER_CORE),stat=ier)
  if(ier /= 0) stop 'error memory allocation merged version'
  allocate(zelm_store_outer_core(NGNOD,NSPEC_OUTER_CORE),stat=ier)
  if(ier /= 0) stop 'error memory allocation merged version'

  allocate(xelm_store_inner_core(NGNOD,NSPEC_INNER_CORE),stat=ier)
  if(ier /= 0) stop 'error memory allocation merged version'
  allocate(yelm_store_inner_core(NGNOD,NSPEC_INNER_CORE),stat=ier)
  if(ier /= 0) stop 'error memory allocation merged version'
  allocate(zelm_store_inner_core(NGNOD,NSPEC_INNER_CORE),stat=ier)
  if(ier /= 0) stop 'error memory allocation merged version'

