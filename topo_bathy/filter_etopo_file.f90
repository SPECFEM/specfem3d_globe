
  program filter_topo_file

! smooth topo file using a window filter

  implicit none

  include "../../constants.h"

! filter final surface using box filter
  integer, parameter :: SIZE_FILTER_ONE_SIDE = 2

! impose min and max of topography and bathymetry
  integer, parameter :: MIN_TOPO = -5000
  integer, parameter :: MAX_TOPO = +5000

  integer ix,iy,ixconv,iyconv,i,ic,ixcount,iycount
  integer ixbox,iybox,ixval,iyval

  double precision rlon,rlat,rx,ry,a,b
  double precision sumval,value,dist,sigma

! use integer array to store values
  integer itopo(NX_BATHY,NY_BATHY)
  integer itopo_filtered(NX_BATHY,NY_BATHY)

  integer itopo_x,itopo_y

  print *
  print *,'size of window filter is ',2*SIZE_FILTER_ONE_SIDE + 1
  print *

  print *,'imposing min and max of topography and bathymetry = ',MIN_TOPO,MAX_TOPO
  print *

  open(unit=13,file='topo_bathy_etopo5.dat',status='old')
  do itopo_y=1,NY_BATHY
    do itopo_x=1,NX_BATHY
      read(13,*) itopo(itopo_x,itopo_y)
    enddo
  enddo
  close(13)

  print *,'min and max of original topo file is ',minval(itopo),maxval(itopo)

! impose min and max values
  do itopo_y=1,NY_BATHY
    do itopo_x=1,NX_BATHY
      if(itopo(itopo_x,itopo_y) < MIN_TOPO) itopo(itopo_x,itopo_y) = MIN_TOPO
      if(itopo(itopo_x,itopo_y) > MAX_TOPO) itopo(itopo_x,itopo_y) = MAX_TOPO
    enddo
  enddo

! filter final surface using box filter
    print *,'filtering topography data file using box filter'
    do iy = 1,NY_BATHY
      print *,'doing iy = ',iy,' out of ',NY_BATHY
      do ix = 1,NX_BATHY
        sumval = 0.d0
        do iybox = iy-SIZE_FILTER_ONE_SIDE,iy+SIZE_FILTER_ONE_SIDE
          do ixbox = ix-SIZE_FILTER_ONE_SIDE,ix+SIZE_FILTER_ONE_SIDE
            ixval = ixbox
            iyval = iybox
            if(ixval < 1) ixval = NX_BATHY - abs(ixval)
            if(iyval < 1) iyval = NY_BATHY - abs(iyval)
            if(ixval > NX_BATHY) ixval = ixval - NX_BATHY
            if(iyval > NY_BATHY) iyval = iyval - NY_BATHY
            sumval = sumval + dble(itopo(ixval,iyval))
          enddo
        enddo
        itopo_filtered(ix,iy) = nint(sumval/dble((2*SIZE_FILTER_ONE_SIDE+1)**2))
      enddo
    enddo

  print *,'min and max of filtered topo file is ',minval(itopo_filtered),maxval(itopo_filtered)

  open(unit=13,file='topo_bathy_etopo5_filtered.dat',status='unknown')
  do itopo_y=1,NY_BATHY
    do itopo_x=1,NX_BATHY
      write(13,*) itopo_filtered(itopo_x,itopo_y)
    enddo
  enddo
  close(13)

  end program filter_topo_file

