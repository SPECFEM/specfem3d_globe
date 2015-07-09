program normal_plane

! this program computes the normal to the source-receiver plane
! Qinya Liu, Caltech, May 2007

  use sub_slice_number

  implicit none

  character(len=100) :: ch_lats, ch_lons, ch_latr, ch_lonr
  real(kind=CUSTOM_REAL) :: lats, lons, latr, lonr, xs, ys, zs, nx, ny, nz, xc, yc, zc

  call getarg(1,ch_lats)
  call getarg(2,ch_lons)
  call getarg(3,ch_latr)
  call getarg(4,ch_lonr)

  if (trim(ch_lats) == '' .or. trim(ch_lons) == '' .or. trim(ch_latr) == '' .or. trim(ch_lonr) == '') &
             stop 'Usage: normal_plane lats lons latr lonr'

  read(ch_lats,*) lats
  read(ch_lons,*) lons
  read(ch_latr,*) latr
  read(ch_lonr,*) lonr

  call latlon2norm(lats,lons,latr,lonr,xs,ys,zs,nx,ny,nz,xc,yc,zc)

  write(*,'(3f10.3,a)') xs, ys, zs, '   --- source location'

  write(*,'(3f10.3,a)') nx, ny, nz, '   --- normal of s-r plane'

  write(*,'(3f10.3,a)') xc, yc, zc, '   --- center of s-r gcarc'


end program normal_plane
