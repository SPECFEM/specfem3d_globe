
  program merge_filter_ori_bathy_topo

!! DK DK smooth original bathy and topo file

  implicit none

  include "../../constants.h"

! filter final surface using box filter
  integer, parameter :: SIZE_FILTER_ONE_SIDE = 2

  integer ix,iy,ixconv,iyconv,i,ic,ixcount,iycount
  integer ixbox,iybox,ixval,iyval

  double precision rlon,rlat,rx,ry,a,b
  double precision sumval,value,dist,sigma

! use integer array to store values
  integer itopo(NX_BATHY,NY_BATHY)
  integer itopo_filtered(NX_BATHY,NY_BATHY)

  integer itopo_x,itopo_y

  open(unit=13,file='topo_bathy_etopo5_harvard.dat_ori',status='old')
  do itopo_y=1,NY_BATHY
    do itopo_x=1,NX_BATHY
      read(13,*) itopo(itopo_x,itopo_y)
    enddo
  enddo
  close(13)

! filter final surface using box filter
    print *,'filtering final surface using box filter'
    do iy = 1,NY_BATHY
      print *,'doing iy = ',iy,' out of ',NY_BATHY
      do ix = 1,NX_BATHY
        sumval = 0.d0
        do iybox = iy-SIZE_FILTER_ONE_SIDE,iy+SIZE_FILTER_ONE_SIDE
          do ixbox = ix-SIZE_FILTER_ONE_SIDE,ix+SIZE_FILTER_ONE_SIDE
            ixval = ixbox
            iyval = iybox
            if(ixval < 1) ixval = 1
            if(iyval < 1) iyval = 1
            if(ixval > NX_BATHY) ixval = NX_BATHY
            if(iyval > NY_BATHY) iyval = NY_BATHY
            sumval = sumval + dble(itopo(ixval,iyval))
          enddo
        enddo
        itopo_filtered(ix,iy) = nint(sumval/dble((2*SIZE_FILTER_ONE_SIDE+1)**2))
      enddo
    enddo
    itopo(:,:) = itopo_filtered(:,:)

  open(unit=13,file='topo_bathy_etopo5_harvard.dat',status='unknown')
  do itopo_y=1,NY_BATHY
    do itopo_x=1,NX_BATHY
      write(13,*) itopo(itopo_x,itopo_y)
    enddo
  enddo
  close(13)

  end program merge_filter_ori_bathy_topo

