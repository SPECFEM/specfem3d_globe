
!! DK DK added this for merged version

  if(ipass == 2) then

! suppress fictitious elements in central cube
! also take into account the fact that array idoubling is not allocated for the outer core
  add_contrib_this_element = .true.
  if(iregion_code == IREGION_INNER_CORE) then
    if(idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) add_contrib_this_element = .false.
  endif

  if(add_contrib_this_element) then

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        weight = wxgll(i)*wygll(j)*wzgll(k)

        iglobnum = ibool(i,j,k,ispec)

! compute the jacobian
        xixl = xixstore(i,j,k)
        xiyl = xiystore(i,j,k)
        xizl = xizstore(i,j,k)
        etaxl = etaxstore(i,j,k)
        etayl = etaystore(i,j,k)
        etazl = etazstore(i,j,k)
        gammaxl = gammaxstore(i,j,k)
        gammayl = gammaystore(i,j,k)
        gammazl = gammazstore(i,j,k)

        jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))

! definition depends if region is fluid or solid
  if(iregion_code == IREGION_CRUST_MANTLE .or. iregion_code == IREGION_INNER_CORE) then

! distinguish between single and double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
      rmass(iglobnum) = rmass(iglobnum) + &
             sngl(dble(rhostore_local(i,j,k)) * dble(jacobianl) * weight)
    else
      rmass(iglobnum) = rmass(iglobnum) + rhostore_local(i,j,k) * jacobianl * weight
    endif

! fluid in outer core
  else if(iregion_code == IREGION_OUTER_CORE) then

! no anisotropy in the fluid, use kappav

! distinguish between single and double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
      rmass(iglobnum) = rmass(iglobnum) + &
             sngl(dble(jacobianl) * weight * dble(rhostore_local(i,j,k)) / dble(kappavstore_local(i,j,k)))
    else
      rmass(iglobnum) = rmass(iglobnum) + &
             jacobianl * weight * rhostore_local(i,j,k) / kappavstore_local(i,j,k)
    endif

  else
    call exit_MPI(myrank,'wrong region code')
  endif

      enddo
    enddo
  enddo

  endif ! of exclusion of fictitious inner core elements

  endif ! of ipass == 2

