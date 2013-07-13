
  program expand_compute_strain_product

!! DK DK July 2013: expand the calculations in compute_strain_product() in order to speed up the code

  implicit none

  integer :: i,j,p

  ! Computing the 21 strain products without assuming eps(i)*b_eps(j) = eps(j)*b_eps(i)
  p=1
  do i=1,6
    do j=i,6
      print *,'prod(',p,')=eps',i,'*b_eps',j
      if(j>i) then
        print *,'prod(',p,')=prod(',p,')+eps',j,'*b_eps',i
        if(j>3 .and. i<4) print *,'prod(',p,') = prod(',p,') * 2._CUSTOM_REAL'
      endif
      if(i>3) print *,'prod(',p,') = prod(',p,') * 4._CUSTOM_REAL'
      p=p+1
    enddo
  enddo

  end program expand_compute_strain_product

