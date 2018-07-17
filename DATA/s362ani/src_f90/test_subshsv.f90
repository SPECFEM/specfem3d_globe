
  program test_subshsv

! --- subshsv returns velocity and density perturbations in percent
! --- subtopo returns topography of the 410- and 650-km discontinuities in km (depression)

  implicit none

  real(kind=4) :: xlat,xcolat,xlon,xdep,xrad
  real(kind=4) :: vshout,vsvout,vphout,vpvout,etaout,rhoout
  real(kind=4) :: topo410out,topo650out
  integer ifknowmodel

! ---

  print *,'xlat='
  read(5,*)xlat
  print *,'xlon='
  read(5,*)xlon
  print *,'xdep='
  read(5,*)xdep

  xcolat=90.0-xlat
  xrad=6371.0-xdep
  ifknowmodel=0

  call subshsv(xcolat,xlon,xrad, &
                   vshout,vsvout,vphout,vpvout,etaout,rhoout, &
                     ifknowmodel)
  write(*,"('    vsh       vsv       vph       vpv       eta    rho    ')")
  write(*,"(6f10.5)") vshout,vsvout,vphout,vpvout,etaout,rhoout

  call subtopo(xcolat,xlon,topo410out,topo650out,ifknowmodel)
  write(*,"('   topo410    topo650 ')")
  write(*,"(2f11.5)") topo410out,topo650out

  end program test_subshsv

