!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology August 2003
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

! create AVS or DX 2D data for the faces of the global chunks,
! to be recombined in postprocessing
  subroutine write_AVS_DX_global_chunks_data(myrank,prname,nspec,iboun, &
        ibool,idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool, &
        npointot,rhostore,kappavstore,muvstore, &
        nspl,rspl,espl,espl2, &
        ELLIPTICITY,THREE_D,CRUSTAL,ONE_CRUST)

  implicit none

  include "constants.h"

  integer nspec,myrank
  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

  integer idoubling(nspec)

  logical iboun(6,nspec),ELLIPTICITY,THREE_D,CRUSTAL,ONE_CRUST

  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

  real(kind=CUSTOM_REAL) kappavstore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) muvstore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) rhostore(NGLLX,NGLLY,NGLLZ,nspec)

! logical mask used to output global points only once
  integer npointot
  logical mask_ibool(npointot)

! numbering of global AVS or DX points
  integer num_ibool_AVS_DX(npointot)

  integer ispec
  integer i,j,k,np
  integer, dimension(8) :: iglobval
  integer npoin,numpoin,nspecface,ispecface

  real(kind=CUSTOM_REAL) vmin,vmax

  double precision r,rho,vp,vs,Qkappa,Qmu
  double precision x,y,z,theta,phi_dummy,cost,p20,ell,factor
  real(kind=CUSTOM_REAL) dvp,dvs

! for ellipticity
  integer nspl
  double precision rspl(NR),espl(NR),espl2(NR)

! processor identification
  character(len=150) prname

! writing points
  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXpointschunks.txt',status='unknown')
  open(unit=11,file=prname(1:len_trim(prname))//'AVS_DXpointschunks_stability.txt',status='unknown')

! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  nspecface = 0

! mark global AVS or DX points
  do ispec=1,nspec
! only if on face
  if(iboun(1,ispec) .or. iboun(2,ispec) .or. &
              iboun(3,ispec) .or. iboun(4,ispec)) then
    iglobval(1)=ibool(1,1,1,ispec)
    iglobval(2)=ibool(NGLLX,1,1,ispec)
    iglobval(3)=ibool(NGLLX,NGLLY,1,ispec)
    iglobval(4)=ibool(1,NGLLY,1,ispec)
    iglobval(5)=ibool(1,1,NGLLZ,ispec)
    iglobval(6)=ibool(NGLLX,1,NGLLZ,ispec)
    iglobval(7)=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglobval(8)=ibool(1,NGLLY,NGLLZ,ispec)

! face xi = xi_min
  if(iboun(1,ispec)) then
    nspecface = nspecface + 1
    mask_ibool(iglobval(1)) = .true.
    mask_ibool(iglobval(4)) = .true.
    mask_ibool(iglobval(8)) = .true.
    mask_ibool(iglobval(5)) = .true.
  endif

! face xi = xi_max
  if(iboun(2,ispec)) then
    nspecface = nspecface + 1
    mask_ibool(iglobval(2)) = .true.
    mask_ibool(iglobval(3)) = .true.
    mask_ibool(iglobval(7)) = .true.
    mask_ibool(iglobval(6)) = .true.
  endif

! face eta = eta_min
  if(iboun(3,ispec)) then
    nspecface = nspecface + 1
    mask_ibool(iglobval(1)) = .true.
    mask_ibool(iglobval(2)) = .true.
    mask_ibool(iglobval(6)) = .true.
    mask_ibool(iglobval(5)) = .true.
  endif

! face eta = eta_max
  if(iboun(4,ispec)) then
    nspecface = nspecface + 1
    mask_ibool(iglobval(4)) = .true.
    mask_ibool(iglobval(3)) = .true.
    mask_ibool(iglobval(7)) = .true.
    mask_ibool(iglobval(8)) = .true.
  endif

  endif
  enddo

! count global number of AVS or DX points
  npoin = count(mask_ibool(:))

! number of points in AVS or DX file
  write(10,*) npoin

! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

! output global AVS or DX points
  numpoin = 0
  do ispec=1,nspec
! only if on face
  if(iboun(1,ispec) .or. iboun(2,ispec) .or. &
              iboun(3,ispec) .or. iboun(4,ispec)) then
    iglobval(1)=ibool(1,1,1,ispec)
    iglobval(2)=ibool(NGLLX,1,1,ispec)
    iglobval(3)=ibool(NGLLX,NGLLY,1,ispec)
    iglobval(4)=ibool(1,NGLLY,1,ispec)
    iglobval(5)=ibool(1,1,NGLLZ,ispec)
    iglobval(6)=ibool(NGLLX,1,NGLLZ,ispec)
    iglobval(7)=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglobval(8)=ibool(1,NGLLY,NGLLZ,ispec)

! face xi = xi_min
  if(iboun(1,ispec)) then

    if(.not. mask_ibool(iglobval(1))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(1)) = numpoin
      write(10,*) numpoin,sngl(xstore(1,1,1,ispec)), &
              sngl(ystore(1,1,1,ispec)),sngl(zstore(1,1,1,ispec))
      vmax = sqrt((kappavstore(1,1,1,ispec)+4.*muvstore(1,1,1,ispec)/3.)/rhostore(1,1,1,ispec))
      vmin = sqrt(muvstore(1,1,1,ispec)/rhostore(1,1,1,ispec))
! particular case of the outer core (muvstore contains 1/rho)
  if(idoubling(ispec) == IFLAG_TOP_OUTER_CORE .or. &
     idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL .or. &
     idoubling(ispec) == IFLAG_BOTTOM_OUTER_CORE) then
    r = dsqrt(xstore(1,1,1,ispec)**2 + ystore(1,1,1,ispec)**2 + zstore(1,1,1,ispec)**2)
    call prem_display_outer_core(myrank,r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
    vmax = vp
    vmin = vp
  endif
      if(vmin == 0.0) vmin=vmax
      write(11,*) numpoin,vmin,vmax
    endif

    if(.not. mask_ibool(iglobval(4))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(4)) = numpoin
      write(10,*) numpoin,sngl(xstore(1,NGLLY,1,ispec)), &
              sngl(ystore(1,NGLLY,1,ispec)),sngl(zstore(1,NGLLY,1,ispec))
      vmax = sqrt((kappavstore(1,NGLLY,1,ispec)+4.*muvstore(1,NGLLY,1,ispec)/3.)/rhostore(1,NGLLY,1,ispec))
      vmin = sqrt(muvstore(1,NGLLY,1,ispec)/rhostore(1,NGLLY,1,ispec))
! particular case of the outer core (muvstore contains 1/rho)
  if(idoubling(ispec) == IFLAG_TOP_OUTER_CORE .or. &
     idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL .or. &
     idoubling(ispec) == IFLAG_BOTTOM_OUTER_CORE) then
    r = dsqrt(xstore(1,NGLLY,1,ispec)**2 + ystore(1,NGLLY,1,ispec)**2 + zstore(1,NGLLY,1,ispec)**2)
    call prem_display_outer_core(myrank,r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
    vmax = vp
    vmin = vp
  endif
      if(vmin == 0.0) vmin=vmax
      write(11,*) numpoin,vmin,vmax
    endif

    if(.not. mask_ibool(iglobval(8))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(8)) = numpoin
      write(10,*) numpoin,sngl(xstore(1,NGLLY,NGLLZ,ispec)), &
              sngl(ystore(1,NGLLY,NGLLZ,ispec)),sngl(zstore(1,NGLLY,NGLLZ,ispec))
      vmax = sqrt((kappavstore(1,NGLLY,NGLLZ,ispec)+4.*muvstore(1,NGLLY,NGLLZ,ispec)/3.)/rhostore(1,NGLLY,NGLLZ,ispec))
      vmin = sqrt(muvstore(1,NGLLY,NGLLZ,ispec)/rhostore(1,NGLLY,NGLLZ,ispec))
! particular case of the outer core (muvstore contains 1/rho)
  if(idoubling(ispec) == IFLAG_TOP_OUTER_CORE .or. &
     idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL .or. &
     idoubling(ispec) == IFLAG_BOTTOM_OUTER_CORE) then
    r = dsqrt(xstore(1,NGLLY,NGLLZ,ispec)**2 + ystore(1,NGLLY,NGLLZ,ispec)**2 + zstore(1,NGLLY,NGLLZ,ispec)**2)
    call prem_display_outer_core(myrank,r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
    vmax = vp
    vmin = vp
  endif
      if(vmin == 0.0) vmin=vmax
      write(11,*) numpoin,vmin,vmax
    endif

    if(.not. mask_ibool(iglobval(5))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(5)) = numpoin
      write(10,*) numpoin,sngl(xstore(1,1,NGLLZ,ispec)), &
              sngl(ystore(1,1,NGLLZ,ispec)),sngl(zstore(1,1,NGLLZ,ispec))
      vmax = sqrt((kappavstore(1,1,NGLLZ,ispec)+4.*muvstore(1,1,NGLLZ,ispec)/3.)/rhostore(1,1,NGLLZ,ispec))
      vmin = sqrt(muvstore(1,1,NGLLZ,ispec)/rhostore(1,1,NGLLZ,ispec))
! particular case of the outer core (muvstore contains 1/rho)
  if(idoubling(ispec) == IFLAG_TOP_OUTER_CORE .or. &
     idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL .or. &
     idoubling(ispec) == IFLAG_BOTTOM_OUTER_CORE) then
    r = dsqrt(xstore(1,1,NGLLZ,ispec)**2 + ystore(1,1,NGLLZ,ispec)**2 + zstore(1,1,NGLLZ,ispec)**2)
    call prem_display_outer_core(myrank,r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
    vmax = vp
    vmin = vp
  endif
      if(vmin == 0.0) vmin=vmax
      write(11,*) numpoin,vmin,vmax
    endif

    mask_ibool(iglobval(1)) = .true.
    mask_ibool(iglobval(4)) = .true.
    mask_ibool(iglobval(8)) = .true.
    mask_ibool(iglobval(5)) = .true.
  endif

! face xi = xi_max
  if(iboun(2,ispec)) then

    if(.not. mask_ibool(iglobval(2))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(2)) = numpoin
      write(10,*) numpoin,sngl(xstore(NGLLX,1,1,ispec)), &
              sngl(ystore(NGLLX,1,1,ispec)),sngl(zstore(NGLLX,1,1,ispec))
      vmax = sqrt((kappavstore(NGLLX,1,1,ispec)+4.*muvstore(NGLLX,1,1,ispec)/3.)/rhostore(NGLLX,1,1,ispec))
      vmin = sqrt(muvstore(NGLLX,1,1,ispec)/rhostore(NGLLX,1,1,ispec))
! particular case of the outer core (muvstore contains 1/rho)
  if(idoubling(ispec) == IFLAG_TOP_OUTER_CORE .or. &
     idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL .or. &
     idoubling(ispec) == IFLAG_BOTTOM_OUTER_CORE) then
    r = dsqrt(xstore(NGLLX,1,1,ispec)**2 + ystore(NGLLX,1,1,ispec)**2 + zstore(NGLLX,1,1,ispec)**2)
    call prem_display_outer_core(myrank,r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
    vmax = vp
    vmin = vp
  endif
      if(vmin == 0.0) vmin=vmax
      write(11,*) numpoin,vmin,vmax
    endif

    if(.not. mask_ibool(iglobval(3))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(3)) = numpoin
      write(10,*) numpoin,sngl(xstore(NGLLX,NGLLY,1,ispec)), &
              sngl(ystore(NGLLX,NGLLY,1,ispec)),sngl(zstore(NGLLX,NGLLY,1,ispec))
      vmax = sqrt((kappavstore(NGLLX,NGLLY,1,ispec)+4.*muvstore(NGLLX,NGLLY,1,ispec)/3.)/rhostore(NGLLX,NGLLY,1,ispec))
      vmin = sqrt(muvstore(NGLLX,NGLLY,1,ispec)/rhostore(NGLLX,NGLLY,1,ispec))
! particular case of the outer core (muvstore contains 1/rho)
  if(idoubling(ispec) == IFLAG_TOP_OUTER_CORE .or. &
     idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL .or. &
     idoubling(ispec) == IFLAG_BOTTOM_OUTER_CORE) then
    r = dsqrt(xstore(NGLLX,NGLLY,1,ispec)**2 + ystore(NGLLX,NGLLY,1,ispec)**2 + zstore(NGLLX,NGLLY,1,ispec)**2)
    call prem_display_outer_core(myrank,r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
    vmax = vp
    vmin = vp
  endif
      if(vmin == 0.0) vmin=vmax
      write(11,*) numpoin,vmin,vmax
    endif

    if(.not. mask_ibool(iglobval(7))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(7)) = numpoin
      write(10,*) numpoin,sngl(xstore(NGLLX,NGLLY,NGLLZ,ispec)), &
              sngl(ystore(NGLLX,NGLLY,NGLLZ,ispec)),sngl(zstore(NGLLX,NGLLY,NGLLZ,ispec))
      vmax = sqrt((kappavstore(NGLLX,NGLLY,NGLLZ,ispec)+4.*muvstore(NGLLX,NGLLY,NGLLZ,ispec)/3.)/rhostore(NGLLX,NGLLY,NGLLZ,ispec))
      vmin = sqrt(muvstore(NGLLX,NGLLY,NGLLZ,ispec)/rhostore(NGLLX,NGLLY,NGLLZ,ispec))
! particular case of the outer core (muvstore contains 1/rho)
  if(idoubling(ispec) == IFLAG_TOP_OUTER_CORE .or. &
     idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL .or. &
     idoubling(ispec) == IFLAG_BOTTOM_OUTER_CORE) then
    r = dsqrt(xstore(NGLLX,NGLLY,NGLLZ,ispec)**2 + ystore(NGLLX,NGLLY,NGLLZ,ispec)**2 + zstore(NGLLX,NGLLY,NGLLZ,ispec)**2)
    call prem_display_outer_core(myrank,r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
    vmax = vp
    vmin = vp
  endif
      if(vmin == 0.0) vmin=vmax
      write(11,*) numpoin,vmin,vmax
    endif

    if(.not. mask_ibool(iglobval(6))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(6)) = numpoin
      write(10,*) numpoin,sngl(xstore(NGLLX,1,NGLLZ,ispec)), &
              sngl(ystore(NGLLX,1,NGLLZ,ispec)),sngl(zstore(NGLLX,1,NGLLZ,ispec))
      vmax = sqrt((kappavstore(NGLLX,1,NGLLZ,ispec)+4.*muvstore(NGLLX,1,NGLLZ,ispec)/3.)/rhostore(NGLLX,1,NGLLZ,ispec))
      vmin = sqrt(muvstore(NGLLX,1,NGLLZ,ispec)/rhostore(NGLLX,1,NGLLZ,ispec))
! particular case of the outer core (muvstore contains 1/rho)
  if(idoubling(ispec) == IFLAG_TOP_OUTER_CORE .or. &
     idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL .or. &
     idoubling(ispec) == IFLAG_BOTTOM_OUTER_CORE) then
    r = dsqrt(xstore(NGLLX,1,NGLLZ,ispec)**2 + ystore(NGLLX,1,NGLLZ,ispec)**2 + zstore(NGLLX,1,NGLLZ,ispec)**2)
    call prem_display_outer_core(myrank,r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
    vmax = vp
    vmin = vp
  endif
      if(vmin == 0.0) vmin=vmax
      write(11,*) numpoin,vmin,vmax
    endif

    mask_ibool(iglobval(2)) = .true.
    mask_ibool(iglobval(3)) = .true.
    mask_ibool(iglobval(7)) = .true.
    mask_ibool(iglobval(6)) = .true.
  endif

! face eta = eta_min
  if(iboun(3,ispec)) then

    if(.not. mask_ibool(iglobval(1))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(1)) = numpoin
      write(10,*) numpoin,sngl(xstore(1,1,1,ispec)), &
              sngl(ystore(1,1,1,ispec)),sngl(zstore(1,1,1,ispec))
      vmax = sqrt((kappavstore(1,1,1,ispec)+4.*muvstore(1,1,1,ispec)/3.)/rhostore(1,1,1,ispec))
      vmin = sqrt(muvstore(1,1,1,ispec)/rhostore(1,1,1,ispec))
! particular case of the outer core (muvstore contains 1/rho)
  if(idoubling(ispec) == IFLAG_TOP_OUTER_CORE .or. &
     idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL .or. &
     idoubling(ispec) == IFLAG_BOTTOM_OUTER_CORE) then
    r = dsqrt(xstore(1,1,1,ispec)**2 + ystore(1,1,1,ispec)**2 + zstore(1,1,1,ispec)**2)
    call prem_display_outer_core(myrank,r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
    vmax = vp
    vmin = vp
  endif
      if(vmin == 0.0) vmin=vmax
      write(11,*) numpoin,vmin,vmax
    endif

    if(.not. mask_ibool(iglobval(2))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(2)) = numpoin
      write(10,*) numpoin,sngl(xstore(NGLLX,1,1,ispec)), &
              sngl(ystore(NGLLX,1,1,ispec)),sngl(zstore(NGLLX,1,1,ispec))
      vmax = sqrt((kappavstore(NGLLX,1,1,ispec)+4.*muvstore(NGLLX,1,1,ispec)/3.)/rhostore(NGLLX,1,1,ispec))
      vmin = sqrt(muvstore(NGLLX,1,1,ispec)/rhostore(NGLLX,1,1,ispec))
! particular case of the outer core (muvstore contains 1/rho)
  if(idoubling(ispec) == IFLAG_TOP_OUTER_CORE .or. &
     idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL .or. &
     idoubling(ispec) == IFLAG_BOTTOM_OUTER_CORE) then
    r = dsqrt(xstore(NGLLX,1,1,ispec)**2 + ystore(NGLLX,1,1,ispec)**2 + zstore(NGLLX,1,1,ispec)**2)
    call prem_display_outer_core(myrank,r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
    vmax = vp
    vmin = vp
  endif
      if(vmin == 0.0) vmin=vmax
      write(11,*) numpoin,vmin,vmax
    endif

    if(.not. mask_ibool(iglobval(6))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(6)) = numpoin
      write(10,*) numpoin,sngl(xstore(NGLLX,1,NGLLZ,ispec)), &
              sngl(ystore(NGLLX,1,NGLLZ,ispec)),sngl(zstore(NGLLX,1,NGLLZ,ispec))
      vmax = sqrt((kappavstore(NGLLX,1,NGLLZ,ispec)+4.*muvstore(NGLLX,1,NGLLZ,ispec)/3.)/rhostore(NGLLX,1,NGLLZ,ispec))
      vmin = sqrt(muvstore(NGLLX,1,NGLLZ,ispec)/rhostore(NGLLX,1,NGLLZ,ispec))
! particular case of the outer core (muvstore contains 1/rho)
  if(idoubling(ispec) == IFLAG_TOP_OUTER_CORE .or. &
     idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL .or. &
     idoubling(ispec) == IFLAG_BOTTOM_OUTER_CORE) then
    r = dsqrt(xstore(NGLLX,1,NGLLZ,ispec)**2 + ystore(NGLLX,1,NGLLZ,ispec)**2 + zstore(NGLLX,1,NGLLZ,ispec)**2)
    call prem_display_outer_core(myrank,r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
    vmax = vp
    vmin = vp
  endif
      if(vmin == 0.0) vmin=vmax
      write(11,*) numpoin,vmin,vmax
    endif

    if(.not. mask_ibool(iglobval(5))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(5)) = numpoin
      write(10,*) numpoin,sngl(xstore(1,1,NGLLZ,ispec)), &
              sngl(ystore(1,1,NGLLZ,ispec)),sngl(zstore(1,1,NGLLZ,ispec))
      vmax = sqrt((kappavstore(1,1,NGLLZ,ispec)+4.*muvstore(1,1,NGLLZ,ispec)/3.)/rhostore(1,1,NGLLZ,ispec))
      vmin = sqrt(muvstore(1,1,NGLLZ,ispec)/rhostore(1,1,NGLLZ,ispec))
! particular case of the outer core (muvstore contains 1/rho)
  if(idoubling(ispec) == IFLAG_TOP_OUTER_CORE .or. &
     idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL .or. &
     idoubling(ispec) == IFLAG_BOTTOM_OUTER_CORE) then
    r = dsqrt(xstore(1,1,NGLLZ,ispec)**2 + ystore(1,1,NGLLZ,ispec)**2 + zstore(1,1,NGLLZ,ispec)**2)
    call prem_display_outer_core(myrank,r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
    vmax = vp
    vmin = vp
  endif
      if(vmin == 0.0) vmin=vmax
      write(11,*) numpoin,vmin,vmax
    endif

    mask_ibool(iglobval(1)) = .true.
    mask_ibool(iglobval(2)) = .true.
    mask_ibool(iglobval(6)) = .true.
    mask_ibool(iglobval(5)) = .true.
  endif

! face eta = eta_max
  if(iboun(4,ispec)) then

    if(.not. mask_ibool(iglobval(4))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(4)) = numpoin
      write(10,*) numpoin,sngl(xstore(1,NGLLY,1,ispec)), &
              sngl(ystore(1,NGLLY,1,ispec)),sngl(zstore(1,NGLLY,1,ispec))
      vmax = sqrt((kappavstore(1,NGLLY,1,ispec)+4.*muvstore(1,NGLLY,1,ispec)/3.)/rhostore(1,NGLLY,1,ispec))
      vmin = sqrt(muvstore(1,NGLLY,1,ispec)/rhostore(1,NGLLY,1,ispec))
! particular case of the outer core (muvstore contains 1/rho)
  if(idoubling(ispec) == IFLAG_TOP_OUTER_CORE .or. &
     idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL .or. &
     idoubling(ispec) == IFLAG_BOTTOM_OUTER_CORE) then
    r = dsqrt(xstore(1,NGLLY,1,ispec)**2 + ystore(1,NGLLY,1,ispec)**2 + zstore(1,NGLLY,1,ispec)**2)
    call prem_display_outer_core(myrank,r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
    vmax = vp
    vmin = vp
  endif
      if(vmin == 0.0) vmin=vmax
      write(11,*) numpoin,vmin,vmax
    endif

    if(.not. mask_ibool(iglobval(3))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(3)) = numpoin
      write(10,*) numpoin,sngl(xstore(NGLLX,NGLLY,1,ispec)), &
              sngl(ystore(NGLLX,NGLLY,1,ispec)),sngl(zstore(NGLLX,NGLLY,1,ispec))
      vmax = sqrt((kappavstore(NGLLX,NGLLY,1,ispec)+4.*muvstore(NGLLX,NGLLY,1,ispec)/3.)/rhostore(NGLLX,NGLLY,1,ispec))
      vmin = sqrt(muvstore(NGLLX,NGLLY,1,ispec)/rhostore(NGLLX,NGLLY,1,ispec))
! particular case of the outer core (muvstore contains 1/rho)
  if(idoubling(ispec) == IFLAG_TOP_OUTER_CORE .or. &
     idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL .or. &
     idoubling(ispec) == IFLAG_BOTTOM_OUTER_CORE) then
    r = dsqrt(xstore(NGLLX,NGLLY,1,ispec)**2 + ystore(NGLLX,NGLLY,1,ispec)**2 + zstore(NGLLX,NGLLY,1,ispec)**2)
    call prem_display_outer_core(myrank,r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
    vmax = vp
    vmin = vp
  endif
      if(vmin == 0.0) vmin=vmax
      write(11,*) numpoin,vmin,vmax
    endif

    if(.not. mask_ibool(iglobval(7))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(7)) = numpoin
      write(10,*) numpoin,sngl(xstore(NGLLX,NGLLY,NGLLZ,ispec)), &
              sngl(ystore(NGLLX,NGLLY,NGLLZ,ispec)),sngl(zstore(NGLLX,NGLLY,NGLLZ,ispec))
      vmax = sqrt((kappavstore(NGLLX,NGLLY,NGLLZ,ispec)+4.*muvstore(NGLLX,NGLLY,NGLLZ,ispec)/3.)/rhostore(NGLLX,NGLLY,NGLLZ,ispec))
      vmin = sqrt(muvstore(NGLLX,NGLLY,NGLLZ,ispec)/rhostore(NGLLX,NGLLY,NGLLZ,ispec))
! particular case of the outer core (muvstore contains 1/rho)
  if(idoubling(ispec) == IFLAG_TOP_OUTER_CORE .or. &
     idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL .or. &
     idoubling(ispec) == IFLAG_BOTTOM_OUTER_CORE) then
    r = dsqrt(xstore(NGLLX,NGLLY,NGLLZ,ispec)**2 + ystore(NGLLX,NGLLY,NGLLZ,ispec)**2 + zstore(NGLLX,NGLLY,NGLLZ,ispec)**2)
    call prem_display_outer_core(myrank,r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
    vmax = vp
    vmin = vp
  endif
      if(vmin == 0.0) vmin=vmax
      write(11,*) numpoin,vmin,vmax
    endif

    if(.not. mask_ibool(iglobval(8))) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglobval(8)) = numpoin
      write(10,*) numpoin,sngl(xstore(1,NGLLY,NGLLZ,ispec)), &
              sngl(ystore(1,NGLLY,NGLLZ,ispec)),sngl(zstore(1,NGLLY,NGLLZ,ispec))
      vmax = sqrt((kappavstore(1,NGLLY,NGLLZ,ispec)+4.*muvstore(1,NGLLY,NGLLZ,ispec)/3.)/rhostore(1,NGLLY,NGLLZ,ispec))
      vmin = sqrt(muvstore(1,NGLLY,NGLLZ,ispec)/rhostore(1,NGLLY,NGLLZ,ispec))
! particular case of the outer core (muvstore contains 1/rho)
  if(idoubling(ispec) == IFLAG_TOP_OUTER_CORE .or. &
     idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL .or. &
     idoubling(ispec) == IFLAG_BOTTOM_OUTER_CORE) then
    r = dsqrt(xstore(1,NGLLY,NGLLZ,ispec)**2 + ystore(1,NGLLY,NGLLZ,ispec)**2 + zstore(1,NGLLY,NGLLZ,ispec)**2)
    call prem_display_outer_core(myrank,r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec))
    vmax = vp
    vmin = vp
  endif
      if(vmin == 0.0) vmin=vmax
      write(11,*) numpoin,vmin,vmax
    endif

    mask_ibool(iglobval(4)) = .true.
    mask_ibool(iglobval(3)) = .true.
    mask_ibool(iglobval(7)) = .true.
    mask_ibool(iglobval(8)) = .true.
  endif

  endif
  enddo

! check that number of global points output is okay
  if(numpoin /= npoin) &
    call exit_MPI(myrank,'incorrect number of global points in AVS or DX file creation')

  close(10)
  close(11)

! output global AVS or DX elements

! writing elements
  open(unit=10,file=prname(1:len_trim(prname))//'AVS_DXelementschunks.txt',status='unknown')
  if(THREE_D) &
    open(unit=11,file=prname(1:len_trim(prname))//'AVS_DXelementschunks_dvp_dvs.txt',status='unknown')

! number of elements in AVS or DX file
  write(10,*) nspecface

  ispecface = 0
  do ispec=1,nspec
! only if on face
  if(iboun(1,ispec) .or. iboun(2,ispec) .or. &
              iboun(3,ispec) .or. iboun(4,ispec)) then
    iglobval(1)=ibool(1,1,1,ispec)
    iglobval(2)=ibool(NGLLX,1,1,ispec)
    iglobval(3)=ibool(NGLLX,NGLLY,1,ispec)
    iglobval(4)=ibool(1,NGLLY,1,ispec)
    iglobval(5)=ibool(1,1,NGLLZ,ispec)
    iglobval(6)=ibool(NGLLX,1,NGLLZ,ispec)
    iglobval(7)=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglobval(8)=ibool(1,NGLLY,NGLLZ,ispec)

! include lateral variations if needed

  if(THREE_D) then
!   pick a point within the element and get its radius
    r=dsqrt(xstore(2,2,2,ispec)**2+ystore(2,2,2,ispec)**2+zstore(2,2,2,ispec)**2)

    if(r > RCMB/R_EARTH .and. r < R_UNIT_SPHERE) then
!     average over the element
      dvp = 0.0
      dvs = 0.0
      np =0
      do k=2,NGLLZ-1
        do j=2,NGLLY-1
          do i=2,NGLLX-1
            np=np+1
            x=xstore(i,j,k,ispec)
            y=ystore(i,j,k,ispec)
            z=zstore(i,j,k,ispec)
            r=dsqrt(x*x+y*y+z*z)
!           take out ellipticity
            if(ELLIPTICITY) then
              call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi_dummy)
              cost=dcos(theta)
              p20=0.5d0*(3.0d0*cost*cost-1.0d0)
              call splint(rspl,espl,espl2,nspl,r,ell)
              factor=ONE-(TWO/3.0d0)*ell*p20
              r=r/factor
            endif
            call prem_iso(myrank,r,rho,vp,vs,Qkappa,Qmu,idoubling(ispec),CRUSTAL,ONE_CRUST,.true.)
            dvp = dvp + (sqrt((kappavstore(i,j,k,ispec)+4.*muvstore(i,j,k,ispec)/3.)/rhostore(i,j,k,ispec)) - sngl(vp))/sngl(vp)
            dvs = dvs + (sqrt(muvstore(i,j,k,ispec)/rhostore(i,j,k,ispec)) - sngl(vs))/sngl(vs)
          enddo
        enddo
      enddo
      dvp = dvp / np
      dvs = dvs / np
    else
      dvp = 0.0
      dvs = 0.0
    endif
  endif

! face xi = xi_min
  if(iboun(1,ispec)) then
    ispecface = ispecface + 1
    write(10,*) ispecface,idoubling(ispec),num_ibool_AVS_DX(iglobval(1)), &
                  num_ibool_AVS_DX(iglobval(4)),num_ibool_AVS_DX(iglobval(8)), &
                  num_ibool_AVS_DX(iglobval(5))
    if(THREE_D) write(11,*) ispecface,dvp,dvs
  endif

! face xi = xi_max
  if(iboun(2,ispec)) then
    ispecface = ispecface + 1
    write(10,*) ispecface,idoubling(ispec),num_ibool_AVS_DX(iglobval(2)), &
                  num_ibool_AVS_DX(iglobval(3)),num_ibool_AVS_DX(iglobval(7)), &
                  num_ibool_AVS_DX(iglobval(6))
    if(THREE_D) write(11,*) ispecface,dvp,dvs
  endif

! face eta = eta_min
  if(iboun(3,ispec)) then
    ispecface = ispecface + 1
    write(10,*) ispecface,idoubling(ispec),num_ibool_AVS_DX(iglobval(1)), &
                  num_ibool_AVS_DX(iglobval(2)),num_ibool_AVS_DX(iglobval(6)), &
                  num_ibool_AVS_DX(iglobval(5))
    if(THREE_D) write(11,*) ispecface,dvp,dvs
  endif

! face eta = eta_max
  if(iboun(4,ispec)) then
    ispecface = ispecface + 1
    write(10,*) ispecface,idoubling(ispec),num_ibool_AVS_DX(iglobval(4)), &
                  num_ibool_AVS_DX(iglobval(3)),num_ibool_AVS_DX(iglobval(7)), &
                  num_ibool_AVS_DX(iglobval(8))
    if(THREE_D) write(11,*) ispecface,dvp,dvs
  endif

  endif
  enddo

! check that number of surface elements output is okay
  if(ispecface /= nspecface) &
    call exit_MPI(myrank,'incorrect number of surface elements in AVS or DX file creation')

  close(10)
  if(THREE_D) close(11)

  end subroutine write_AVS_DX_global_chunks_data

