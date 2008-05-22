
!! DK DK added for merged version

! distinguish between single and double precision for reals
  if(CUSTOM_REAL == SIZE_REAL) then
    xelm_store(:,ispec) = sngl(xelm(:))
    yelm_store(:,ispec) = sngl(yelm(:))
    zelm_store(:,ispec) = sngl(zelm(:))
  else
    xelm_store(:,ispec) = xelm(:)
    yelm_store(:,ispec) = yelm(:)
    zelm_store(:,ispec) = zelm(:)
  endif

