program make_az_stations

! this program computes a station list from a source location at given azimuth
! Qinya Liu, Caltech, May 2007


  use sub_slice_number

  implicit none

  character(len=200) :: ch_elon, ch_elat, ch_az, ch_dist, filename
  real(kind=CUSTOM_REAL) :: elon, elat, az, ths, phs, thn, phn, azr, thn_new, phn_new, phr, dist
  integer np, i
  real(kind=CUSTOM_REAL), dimension(1000) :: thr_array, thr_new, phr_new, latr_new, lonr_new

  !============================

  ! set up the input arguments

  call getarg(1,ch_elon)
  call getarg(2,ch_elat)
  call getarg(3,ch_az)

  if (trim(ch_elon) == '' .or. trim(ch_elat) == '' .or. trim(ch_az) == '')  then
    print *, 'Usage: xmake_az_stations elon elat az [dist(degrees)]'
    stop
  endif

  read(ch_elon,*) elon
  read(ch_elat,*) elat
  read(ch_az,*) az

  ! convert to degrees to radians
  ths = (90-elat)/180. * pi
  phs = elon/180. * pi

  ! another point is the north pole
  thn = (90-90)/180. * pi
  phn = 0/180. * pi

  ! receiver azimuth
  azr = az /180. * pi

  ! rotate (ths,phs) to North Pole, and figure out the new location of old N. Pole
  call norm_rot(ths,phs,thn,phn,thn_new,phn_new)
  print *, 'After rotating the source to the N. Pole, old N. Pole become: '
  print *, thn_new*180/pi, phn_new *180/pi

  phr = phn_new - azr ! here could be + azr as well
  write(filename,'(a,i3.3)') 'STATIONS_',floor(az)
  open(11,file=trim(filename),status='unknown')

  ! read in extra distance argument
  call getarg(4,ch_dist)
  if (trim(ch_dist) /= '') then
    read(ch_dist,*) dist
    np = 1
    thr_array(1) = dist / 180. * pi
  else
  ! discretize the receiver garc, and figure out the slice number for individual point
    np = 180
  endif

  do i = 1, np
    if (trim(ch_dist) == '') thr_array(i) = (i*1.) / np * pi
    call norm_rot_back(ths,phs,thr_array(i),phr,thr_new(i),phr_new(i))
    latr_new(i) = (pi/2 - thr_new(i)) * 180/pi
    lonr_new(i) = phr_new(i)  * 180/pi

    write(11,'(a,i3.3,4x,a,4x,f10.3,4x,f10.3,4x,f5.2,4x,f5.2)') 'T',i,'NT',latr_new(i),lonr_new(i),0.,0.

  enddo

  close(11)


end program make_az_stations
