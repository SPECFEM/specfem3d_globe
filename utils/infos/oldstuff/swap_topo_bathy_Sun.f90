
!---- swap longitude blocks in topography and bathymetry file once and for all
!---- to switch from the convention used in the original file on the Sun
!---- to the convention used in SPECFEM3D

  program swap_topo_bathy_Sun

  implicit none

  include "../../constants.h"

! use integer array to store values
  integer ibathy_topo(NX_BATHY,NY_BATHY)

  integer itopo_x,itopo_y

  open(unit=13,file='topo_bathy_etopo4_from_etopo2_subsampled.dat',status='old')
  do itopo_y = 1,NY_BATHY
    do itopo_x = 1,NX_BATHY
      read(13,*) ibathy_topo(itopo_x,itopo_y)
    enddo
  enddo
  close(13)

! blocks of longitude in [0,180[ and [180,360[ must be swapped
! in the final file, itopo_x = 1 should correspond to longitude = 0
! therefore one should see Africa on the left of the JPEG image of topography
  open(unit=13,file='topo_bathy_etopo4_from_etopo2_subsampled_2.dat',status='unknown')
  do itopo_y = 1,NY_BATHY
    do itopo_x = NX_BATHY/2+1,NX_BATHY
      write(13,*) ibathy_topo(itopo_x,itopo_y)
    enddo
    do itopo_x = 1,NX_BATHY/2
      write(13,*) ibathy_topo(itopo_x,itopo_y)
    enddo
  enddo
  close(13)

  end program swap_topo_bathy_Sun

