!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!--------------------------------------------------------------------------------------------------
! S362ANI
!
! A global shear-wave speed model developed by Kustowski et al. [2006].
!
! In this model, radial anisotropy is confined to the uppermost mantle.
! The model (and the corresponding mesh) incorporate
! tomography on the 650-km and 410-km discontinuities in the 1D reference model REF.
!
! s362wmani: A version of S362ANI with anisotropy allowed throughout the mantle.
!
! s362ani_prem: A version of S362ANI calculated using PREM as the 1D reference model
!
! s29ea: A global model with higher resolution in the upper mantle beneath Eurasia
! calculated using REF as the 1D reference model.
!
! note:
! - statistics on topography perturbations: (depending on mesh resolution, here nex=256)
!   410-km: minimum / maximum = -13.48 km / + 13.24 km
!   650-km: minimum / maximum = -14.34 km / + 19.19 km
!--------------------------------------------------------------------------------------------------


  module model_s362ani_par

  ! used for 3D Harvard models s362ani, s362wmani, s362ani_prem and s2.9ea
  integer, parameter :: maxker = 200
  integer, parameter :: maxl   = 72
  integer, parameter :: maxcoe = 2000
  integer, parameter :: maxver = 1000
  integer, parameter :: maxhpa = 2

  real(kind=4),dimension(:,:),allocatable :: xlaspl,xlospl,radspl,coe
  real(kind=4),dimension(:),allocatable :: vercof,vercofd

  integer, dimension(:),allocatable :: lmxhpa,itypehpa,numcoe

  integer,dimension(:,:),allocatable :: itpspl
  integer,dimension(:),allocatable :: ihpakern,ivarkern

  integer :: numker,numhpa

  character(len=80) :: kerstr
  character(len=80) :: refmdl
  character(len=40) :: varstr(maxker)
  character(len=80) :: hsplfl(maxhpa)
  character(len=40) :: dskker(maxker)

  end module model_s362ani_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_s362ani_broadcast(THREE_D_MODEL)

! standard routine to setup model

  use constants, only: myrank
  use model_s362ani_par

  implicit none

  integer,intent(in) :: THREE_D_MODEL

  ! local parameters
  integer :: ier

  ! allocates model arrays
  allocate(xlaspl(maxcoe,maxhpa), &
           xlospl(maxcoe,maxhpa), &
           radspl(maxcoe,maxhpa), &
           coe(maxcoe,maxker), &
           vercof(maxker), &
           vercofd(maxker), &
           itpspl(maxcoe,maxhpa), &
           ihpakern(maxker), &
           ivarkern(maxker), &
           stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating s362ani arrays')

  ! initializes
  xlaspl(:,:) = 0.0
  xlospl(:,:) = 0.0
  radspl(:,:) = 0.0
  coe(:,:) = 0.0

  vercof(:) = 0.0
  vercofd(:) = 0.0

  itpspl(:,:) = 0
  ihpakern(:) = 0
  ivarkern(:) = 0

  ! allocates
  allocate(lmxhpa(maxhpa),itypehpa(maxhpa),numcoe(maxhpa),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating s362ani lmxhpa, .. arrays')

  ! initializes
  lmxhpa(:) = 0
  itypehpa(:) = 0
  numcoe(:) = 0

  ! master process
  if (myrank == 0) call read_model_s362ani(THREE_D_MODEL)

  call bcast_all_singlei(numker)
  call bcast_all_singlei(numhpa)

  call bcast_all_i(lmxhpa,maxhpa)
  call bcast_all_i(itypehpa,maxhpa)
  call bcast_all_i(ihpakern,maxker)
  call bcast_all_i(numcoe,maxhpa)
  call bcast_all_i(ivarkern,maxker)
  call bcast_all_i(itpspl,maxcoe*maxhpa)

  call bcast_all_r(xlaspl,maxcoe*maxhpa)
  call bcast_all_r(xlospl,maxcoe*maxhpa)
  call bcast_all_r(radspl,maxcoe*maxhpa)
  call bcast_all_r(coe,maxcoe*maxker)

  call bcast_all_ch_array(hsplfl,maxhpa,80)
  call bcast_all_ch_array(dskker,maxker,40)

  call bcast_all_ch(kerstr,80)
  call bcast_all_ch(refmdl,80)

  call bcast_all_ch_array(varstr,maxker,40)

  end subroutine model_s362ani_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_model_s362ani(THREE_D_MODEL)

  use model_s362ani_par
  use constants, only: THREE_D_MODEL_S362ANI,THREE_D_MODEL_S362WMANI, &
    THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA

  implicit none

  integer,intent(in) :: THREE_D_MODEL

  ! local parameters
  character(len=128) :: modeldef
  logical :: exists
  integer :: numvar,ihpa

  ! sets model file name
  if (THREE_D_MODEL == THREE_D_MODEL_S362ANI) then
    modeldef='DATA/s362ani/S362ANI'
  else if (THREE_D_MODEL == THREE_D_MODEL_S362WMANI) then
    modeldef='DATA/s362ani/S362WMANI'
  else if (THREE_D_MODEL == THREE_D_MODEL_S362ANI_PREM) then
    modeldef='DATA/s362ani/S362ANI_PREM'
  else if (THREE_D_MODEL == THREE_D_MODEL_S29EA) then
    modeldef='DATA/s362ani/S2.9EA'
  else
    stop 'Error: unknown 3D model in read_model_s362ani() routine'
  endif

  ! reads in model parameters
  inquire(file=modeldef,exist=exists)
  if (exists) then
    call gt3dmodl(modeldef, &
                  numhpa,numker,numcoe,lmxhpa, &
                  ihpakern,itypehpa,coe, &
                  itpspl,xlaspl,xlospl,radspl, &
                  numvar,ivarkern,varstr, &
                  refmdl,kerstr,hsplfl,dskker)
  else
    write(*,"('model ',a,' does not exist')") modeldef(1:len_trim(modeldef))
    stop 'model does not exist in s362_ani'
  endif

  !  checks arrays maximum values
  if (numker > maxker) stop 'Error: numker > maxker in read_model_s362ani() routine'
  do ihpa = 1,numhpa
    if (itypehpa(ihpa) == 1) then
      if (lmxhpa(ihpa) > maxl) stop 'lmxhpa(ihpa) > maxl in read_model_s362ani() routine'

    else if (itypehpa(ihpa) == 2) then
      if (numcoe(ihpa) > maxcoe) stop 'numcoe(ihpa) > maxcoe in read_model_s362ani() routine'

    else
      stop 'problem with itypehpa'
    endif
  enddo

  end subroutine read_model_s362ani

!
!-------------------------------------------------------------------------------------------------
!

  subroutine evradker(depth,string,nker,vercof,dvercof,ierror)

  use constants, only: R_EARTH_KM

  implicit none

  integer :: nker,ierror

  real(kind=4) :: chebyshev(100)
  real(kind=4) :: chebyshev2(100)
  real(kind=4) :: vercof(nker)
  real(kind=4) :: dvercof(nker)
  real(kind=4) :: splpts(100)

  character(len=80) string

  logical :: upper,upper_650
  logical :: lower,lower_650

  real(kind=4), parameter :: r0 = R_EARTH_KM ! 6371.0
  real(kind=4), parameter :: rmoho = r0 - 24.4  ! subtracting the thickness here
  real(kind=4), parameter :: r670 = r0 - 670.0    ! subtracting the thickness here
  real(kind=4), parameter :: r650 = r0 - 650.0    ! subtracting the thickness here
  real(kind=4), parameter :: rcmb = 3480.0

  integer :: i,nspl,nskip,nlower,nupper,iker,lstr

  real(kind=4) :: u,u2,ddep,radius2,radius,depth

  ierror = 0
  lstr = len_trim(string)

  radius = r0 - depth
  ddep = 0.1
  radius2 = r0 - depth + ddep
  upper = .false.
  lower = .false.
  if (radius > rcmb .and. radius < r670) then
    lower = .true.
  else if (radius >= r670 .and. radius < rmoho) then
    upper = .true.
  endif
  upper_650 = .false.
  lower_650 = .false.
  if (radius > rcmb .and. radius < r650) then
    lower_650 = .true.
  else if (radius >= r650 .and. radius < rmoho) then
    upper_650 = .true.
  endif
  do iker = 1,nker
    vercof(iker) = 0.0
    dvercof(iker) = 0.0
  enddo

  if (string(1:16) == 'WDC+SPC_U4L8CHEB') then
    nupper = 5
    nlower = 9
    nskip = 2
    if (upper) then
      u = (radius+radius-rmoho-r670)/(rmoho-r670)
      u2 = (radius2+radius2-rmoho-r670)/(rmoho-r670)
  !   write(*,"('upper mantle:',2f10.3)") u,u2
      call chebyfun(u,13,chebyshev)
      do i = 1+nskip,nskip+nupper
        vercof(i) = chebyshev(i-nskip)
      enddo
      call chebyfun(u2,13,chebyshev2)
      do i = 1+nskip,nskip+nupper
        dvercof(i) = (chebyshev2(i-nskip)-chebyshev(i-nskip))/ddep
      enddo
    else if (lower) then
      u = (radius+radius-r670-rcmb)/(r670-rcmb)
      u2 = (radius2+radius2-r670-rcmb)/(r670-rcmb)
  !   write(*,"('lower mantle:',2f10.3)") u,u2
      call chebyfun(u,13,chebyshev)
      do i = 1+nskip+nupper,nskip+nupper+nlower
        vercof(i) = chebyshev(i-nskip-nupper)
      enddo
      call chebyfun(u2,13,chebyshev2)
      do i = 1+nskip+nupper,nskip+nupper+nlower
        dvercof(i) = (chebyshev2(i-nskip-nupper) - chebyshev(i-nskip-nupper))/ddep
      enddo
    endif
  else if (string(1:13) == 'WDC+SHSVWM20A') then
    nspl=20
    splpts(1)=0.0
    splpts(2)=50.0
    splpts(3)=100.0
    splpts(4)=150.0
    splpts(5)=200.0
    splpts(6)=250.0
    splpts(7)=300.0
    splpts(8)=400.0
    splpts(9)=500.0
    splpts(10)=600.0
    splpts(11)=700.0
    splpts(12)=850.0
    splpts(13)=1050.0
    splpts(14)=1300.0
    splpts(15)=1600.0
    splpts(16)=1900.0
    splpts(17)=2200.0
    splpts(18)=2500.0
    splpts(19)=2700.0
    splpts(20)=2891.0
    call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
    do i=22,27
      vercof(i)=vercof(i-20)
      dvercof(i)=dvercof(i-20)
    enddo
    vercof(1)=1.0
  else if (string(1:16) == 'WDC+XBS_362_U6L8') then
    if (upper) then
      nspl=6
      splpts(1)=24.4
      splpts(2)=100.0
      splpts(3)=225.0
      splpts(4)=350.0
      splpts(5)=500.0
      splpts(6)=670.0
      call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
    else if (lower) then
      nspl=8
      splpts(1)=670.0
      splpts(2)=820.0
      splpts(3)=1320.0
      splpts(4)=1820.0
      splpts(5)=2320.0
      splpts(6)=2550.0
      splpts(7)=2791.0
      splpts(8)=2891.0
      call vbspl(depth,nspl,splpts,vercof(8),dvercof(8))
    endif
    vercof(1)=1.0
!        vercof(16)=1.
!        vercof(17)=1.
!      else if (string(1:21) == 'WDC+ANI_362_U6L8_TOPO') then
!        if (upper) then
!         nspl=6
!         splpts(1)=24.4
!         splpts(2)=100.
!         splpts(3)=225.
!         splpts(4)=350.
!         splpts(5)=500.
!         splpts(6)=670.
!         call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
!         do i=16,21
!          vercof(i)=vercof(i-14)
!          dvercof(i)=dvercof(i-14)
!         enddo
!     else if (lower) then
!      nspl=8
!         splpts(1)=670.
!         splpts(2)=820.
!         splpts(3)=1320.
!         splpts(4)=1820.
!         splpts(5)=2320.
!         splpts(6)=2550.
!         splpts(7)=2791.
!         splpts(8)=2891.
!         call vbspl(depth,nspl,splpts,vercof(8),dvercof(8))
!     endif
!        vercof(1)=1.
!        vercof(22)=1.
!        vercof(23)=1.
!        vercof(24)=1.
!        vercof(25)=1.
  else if ((string(1:lstr) == 'WDC+ANI_362_U6L8' .and. lstr == 16) &
      .or. (string(1:lstr) == 'WDC+ANI_362_U6L8_TOPO'.and.lstr == 21)) then
    if (upper) then
      nspl=6
      splpts(1)=24.4
      splpts(2)=100.0
      splpts(3)=225.0
      splpts(4)=350.0
      splpts(5)=500.0
      splpts(6)=670.0
      call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
      do i=16,21
        vercof(i)=vercof(i-14)
        dvercof(i)=dvercof(i-14)
      enddo
    else if (lower) then
      nspl=8
      splpts(1)=670.0
      splpts(2)=820.0
      splpts(3)=1320.0
      splpts(4)=1820.0
      splpts(5)=2320.0
      splpts(6)=2550.0
      splpts(7)=2791.0
      splpts(8)=2891.0
      call vbspl(depth,nspl,splpts,vercof(8),dvercof(8))
    endif
    vercof(1)=1.0
    vercof(22)=1.0
    vercof(23)=1.0
  else if (string(1:lstr) == 'WDC+WM_362_U6L8' .and. lstr == 15) then
    if (upper) then
      nspl=6
      splpts(1)=24.4
      splpts(2)=100.0
      splpts(3)=225.0
      splpts(4)=350.0
      splpts(5)=500.0
      splpts(6)=670.0
      call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
      do i=16,21
        vercof(i)=vercof(i-14)
        dvercof(i)=dvercof(i-14)
      enddo
    else if (lower) then
      nspl=8
      splpts(1)=670.0
      splpts(2)=820.0
      splpts(3)=1320.0
      splpts(4)=1820.0
      splpts(5)=2320.0
      splpts(6)=2550.0
      splpts(7)=2791.0
      splpts(8)=2891.0
      call vbspl(depth,nspl,splpts,vercof(8),dvercof(8))
      do i=22,29
        vercof(i)=vercof(i-14)
        dvercof(i)=dvercof(i-14)
      enddo
    endif
    vercof(1)=1.0
    vercof(30)=1.0
    vercof(31)=1.0
    vercof(32)=1.0
  else if ((string(1:lstr) == 'WDC+ANI_362_U6L8_650' .and. lstr == 20) &
      .or. (string(1:lstr) == 'WDC+ANI_362_U6L8_TOPO_650'.and.lstr == 25)) then
    if (upper_650) then
      nspl=6
      splpts(1)=24.4
      splpts(2)=100.0
      splpts(3)=225.0
      splpts(4)=350.0
      splpts(5)=500.0
      splpts(6)=650.0
      call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
      do i=16,21
        vercof(i)=vercof(i-14)
        dvercof(i)=dvercof(i-14)
      enddo
    else if (lower_650) then
      nspl=8
      splpts(1)=650.0
      splpts(2)=820.0
      splpts(3)=1320.0
      splpts(4)=1820.0
      splpts(5)=2320.0
      splpts(6)=2550.0
      splpts(7)=2791.0
      splpts(8)=2891.0
      call vbspl(depth,nspl,splpts,vercof(8),dvercof(8))
    endif
    vercof(1)=1.0
    vercof(22)=1.0
    vercof(23)=1.0
  else if (string(1:lstr) == 'WDC+WM_362_U6L8_650' .and. lstr == 19) then
    if (upper_650) then
      nspl=6
      splpts(1)=24.4
      splpts(2)=100.0
      splpts(3)=225.0
      splpts(4)=350.0
      splpts(5)=500.0
      splpts(6)=650.0
      call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
      do i=16,21
        vercof(i)=vercof(i-14)
        dvercof(i)=dvercof(i-14)
      enddo
    else if (lower_650) then
      nspl=8
      splpts(1)=650.0
      splpts(2)=820.0
      splpts(3)=1320.0
      splpts(4)=1820.0
      splpts(5)=2320.0
      splpts(6)=2550.0
      splpts(7)=2791.0
      splpts(8)=2891.0
      call vbspl(depth,nspl,splpts,vercof(8),dvercof(8))
      do i=22,29
        vercof(i)=vercof(i-14)
        dvercof(i)=dvercof(i-14)
      enddo
    endif
    vercof(1)=1.0
    vercof(30)=1.0
    vercof(31)=1.0
    vercof(32)=1.0
  else if (string(1:lstr) == 'WDC+U8L8_650' .and. lstr == 12) then
    if (upper_650) then
      nspl=8
      splpts(1)=24.4
      splpts(2)=75.0
      splpts(3)=150.0
      splpts(4)=225.0
      splpts(5)=300.0
      splpts(6)=410.0
      splpts(7)=530.0
      splpts(8)=650.0
      call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
      do i=18,25
        vercof(i)=vercof(i-16)
        dvercof(i)=dvercof(i-16)
      enddo
    else if (lower_650) then
      nspl=8
      splpts(1)=650.0
      splpts(2)=820.0
      splpts(3)=1320.0
      splpts(4)=1820.0
      splpts(5)=2320.0
      splpts(6)=2550.0
      splpts(7)=2791.0
      splpts(8)=2891.0
      call vbspl(depth,nspl,splpts,vercof(10),dvercof(10))
      do i=26,33
        vercof(i)=vercof(i-16)
        dvercof(i)=dvercof(i-16)
      enddo
    endif
    vercof(1)=1.0
    vercof(34)=1.0
    vercof(35)=1.0
    vercof(36)=1.0
  else if (string(1:lstr) == 'WDC+U8L8_670' .and. lstr == 12) then
    if (upper) then
      nspl=8
      splpts(1)=24.4
      splpts(2)=75.0
      splpts(3)=150.0
      splpts(4)=225.0
      splpts(5)=300.0
      splpts(6)=410.0
      splpts(7)=530.0
      splpts(8)=670.0
      call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
      do i=18,25
        vercof(i)=vercof(i-16)
        dvercof(i)=dvercof(i-16)
      enddo
    else if (lower) then
      nspl=8
      splpts(1)=670.0
      splpts(2)=820.0
      splpts(3)=1320.0
      splpts(4)=1820.0
      splpts(5)=2320.0
      splpts(6)=2550.0
      splpts(7)=2791.0
      splpts(8)=2891.0
      call vbspl(depth,nspl,splpts,vercof(10),dvercof(10))
      do i=26,33
        vercof(i)=vercof(i-16)
        dvercof(i)=dvercof(i-16)
      enddo
    endif
    vercof(1)=1.0
    vercof(34)=1.0
    vercof(35)=1.0
    vercof(36)=1.0
  else if ((string(1:lstr) == 'WDC+U8L8_I1D_650' .and. lstr == 16) &
      .or. (string(1:lstr) == 'WDC+U8L8_I3D_650'.and.lstr == 16)) then
    if (upper_650) then
      nspl=8
      splpts(1)=24.4
      splpts(2)=75.0
      splpts(3)=150.0
      splpts(4)=225.0
      splpts(5)=300.0
      splpts(6)=410.0
      splpts(7)=530.0
      splpts(8)=650.0
      call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
      do i=18,25
        vercof(i)=vercof(i-16)
        dvercof(i)=dvercof(i-16)
      enddo
      do i=37,40
        vercof(i)=vercof(i-35)
        dvercof(i)=dvercof(i-35)
      enddo
      do i=41,44
        vercof(i)=vercof(i-39)
        dvercof(i)=dvercof(i-39)
      enddo
      do i=45,48
        vercof(i)=vercof(i-43)
        dvercof(i)=dvercof(i-43)
      enddo
      do i=49,52
        vercof(i)=vercof(i-47)
        dvercof(i)=dvercof(i-47)
      enddo
    else if (lower_650) then
      nspl=8
      splpts(1)=650.0
      splpts(2)=820.0
      splpts(3)=1320.0
      splpts(4)=1820.0
      splpts(5)=2320.0
      splpts(6)=2550.0
      splpts(7)=2791.0
      splpts(8)=2891.0
      call vbspl(depth,nspl,splpts,vercof(10),dvercof(10))
      do i=26,33
        vercof(i)=vercof(i-16)
        dvercof(i)=dvercof(i-16)
      enddo
    endif
    vercof(1)=1.0
    vercof(34)=1.0
    vercof(35)=1.0
    vercof(36)=1.0
  else if ((string(1:lstr) == 'WDC+I1D_650' .and. lstr == 11) .or. &
          (string(1:lstr) == 'WDC+I3D_650' .and. lstr == 11)) then
    if (upper_650) then
      nspl=8
      splpts(1)=24.4
      splpts(2)=75.0
      splpts(3)=150.0
      splpts(4)=225.0
      splpts(5)=300.0
      splpts(6)=410.0
      splpts(7)=530.0
      splpts(8)=650.0
      call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
      do i=18,25
        vercof(i)=vercof(i-16)
        dvercof(i)=dvercof(i-16)
      enddo
      do i=37,44
        vercof(i)=vercof(i-35)
        dvercof(i)=dvercof(i-35)
      enddo
      do i=53,60
        vercof(i)=vercof(i-51)
        dvercof(i)=dvercof(i-51)
      enddo
      do i=69,76
        vercof(i)=vercof(i-67)
        dvercof(i)=dvercof(i-67)
      enddo
      do i=85,92
        vercof(i)=vercof(i-83)
        dvercof(i)=dvercof(i-83)
      enddo
    else if (lower_650) then
      nspl=8
      splpts(1)=650.0
      splpts(2)=820.0
      splpts(3)=1320.0
      splpts(4)=1820.0
      splpts(5)=2320.0
      splpts(6)=2550.0
      splpts(7)=2791.0
      splpts(8)=2891.0
      call vbspl(depth,nspl,splpts,vercof(10),dvercof(10))
      do i=26,33
        vercof(i)=vercof(i-16)
        dvercof(i)=dvercof(i-16)
      enddo
      do i=45,52
        vercof(i)=vercof(i-35)
        dvercof(i)=dvercof(i-35)
      enddo
      do i=61,68
        vercof(i)=vercof(i-51)
        dvercof(i)=dvercof(i-51)
      enddo
      do i=77,84
        vercof(i)=vercof(i-67)
        dvercof(i)=dvercof(i-67)
      enddo
      do i=93,100
        vercof(i)=vercof(i-83)
        dvercof(i)=dvercof(i-83)
      enddo
    endif
    vercof(1)=1.0
    vercof(34)=1.0
    vercof(35)=1.0
    vercof(36)=1.0
  else if (string(1:lstr) == 'V16A4_V7A4' .and. lstr == 10) then
    if (upper_650) then
      nspl=8
      splpts(1)=24.4
      splpts(2)=75.0
      splpts(3)=150.0
      splpts(4)=225.0
      splpts(5)=300.0
      splpts(6)=410.0
      splpts(7)=530.0
      splpts(8)=650.0
      call vbspl(depth,nspl,splpts,vercof(1),dvercof(1))
      do i=17,20
        vercof(i)=vercof(i-16)
        dvercof(i)=dvercof(i-16)
      enddo
      do i=23,29
        vercof(i)=vercof(i-22)
        dvercof(i)=dvercof(i-22)
      enddo
      do i=30,33
        vercof(i)=vercof(i-29)
        dvercof(i)=dvercof(i-29)
      enddo
    else if (lower_650) then
      nspl=8
      splpts(1)=650.0
      splpts(2)=820.0
      splpts(3)=1320.0
      splpts(4)=1820.0
      splpts(5)=2320.0
      splpts(6)=2550.0
      splpts(7)=2791.0
      splpts(8)=2891.0
      call vbspl(depth,nspl,splpts,vercof(9),dvercof(9))
    endif
    vercof(21)=1.0
    vercof(22)=1.0
  else
    write(*,"('problem 4')")
    write(*,"(a)")string(1:len_trim(string))
    stop
  endif

  end subroutine evradker

!
!-------------------------------------------------------------------------------------------------
!

  subroutine chebyfun(u,kmax,f)

  implicit none

  integer :: kmax

  real(kind=4),intent(in) :: u
  real(kind=4),intent(out) :: f(0:kmax)

  ! local parameters
  integer :: k
  real(kind=4) :: twou
  real(kind=4), dimension(0:13), parameter :: chebycoeff = &
    (/ 0.70710678118655,1.2247448713916,1.0350983390135,1.0145993123918, &
       1.00803225754840,1.0050890913907,1.0035149493262,1.0025740068320, &
       1.00196657023780,1.0015515913133,1.0012554932754,1.0010368069141, &
       1.00087070107920,1.0007415648034 /)

  if (kmax > 13) stop 'Error kmax exceeds the limit in chebyfun'

  f(0) = 1.0
  f(1) = u
  twou = 2.0*u

  do k = 2,kmax
   f(k) = twou*f(k-1)-f(k-2)
  enddo

  do k = 0,kmax
   f(k) = f(k)*chebycoeff(k)
  enddo

  end subroutine chebyfun

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gt3dmodl(targetfile, &
                      numhpa,numker,numcoe,lmxhpa, &
                      ihpakern,itypehpa,coe, &
                      itpspl,xlatspl,xlonspl,radispl, &
                      numvar,ivarkern,varstr, &
                      refmdl,kerstr,hsplfl,dskker)

  use model_s362ani_par, only: maxcoe,maxhpa,maxker

  implicit none

  character(len=128) :: targetfile

  integer :: numhpa,numker

  integer :: numcoe(maxhpa)
  integer :: lmxhpa(maxhpa)
  integer :: ihpakern(maxker)
  integer :: itypehpa(maxhpa)

  real(kind=4) :: coe(maxcoe,maxker)

  ! splines
  integer :: itpspl(maxcoe,maxhpa)
  real(kind=4) :: xlatspl(maxcoe,maxhpa)
  real(kind=4) :: xlonspl(maxcoe,maxhpa)
  real(kind=4) :: radispl(maxcoe,maxhpa)

  integer :: numvar
  integer :: ivarkern(maxker)

  character(len=40) :: varstr(maxker)
  character(len=80) :: refmdl
  character(len=80) :: kerstr
  character(len=80) :: hsplfl(maxhpa)
  character(len=40) :: dskker(maxker)

  ! local parameters
  character(len=80) :: refmodel
  character(len=80) :: kernstri
  character(len=40) :: desckern(maxker)
  character(len=80) :: hsplfile(maxhpa)

  integer :: ihorpar(maxker)
  integer :: ityphpar(maxhpa)
  integer :: ixlspl(maxcoe,maxhpa)
  integer :: lmaxhor(maxhpa)
  integer :: ncoefhor(maxhpa)

  real(kind=4) :: coef(maxcoe,maxker)
  real(kind=4) :: xlaspl(maxcoe,maxhpa)
  real(kind=4) :: xlospl(maxcoe,maxhpa)
  real(kind=4) :: xraspl(maxcoe,maxhpa)

  character(len=40) :: string

  integer :: nhorpar,nmodkern,i,j,lstr,k
  integer :: ierror

  ierror = 0

  call rd3dmodl(targetfile,ierror, &
                nmodkern,nhorpar,ityphpar, &
                ihorpar,lmaxhor,ncoefhor, &
                xlaspl,xlospl,xraspl,ixlspl,coef, &
                hsplfile,refmodel,kernstri,desckern)

  if (nhorpar <= maxhpa) then
    numhpa = nhorpar
  else
    ierror = ierror+1
  endif

  if (nmodkern <= maxker) then
    numker = nmodkern
  else
    ierror = ierror+1
  endif

  do i = 1,nmodkern
    ihpakern(i) = ihorpar(i)
    dskker(i) = desckern(i)
    do j = 1,ncoefhor(ihpakern(i))
      coe(j,i) = coef(j,i)
      !   if (j == 1) write(*,"(e12.4)") coe(j,i)
    enddo
  enddo

  do i = 1,nhorpar
    numcoe(i) = ncoefhor(i)
    lmxhpa(i) = lmaxhor(i)
    itypehpa(i) = ityphpar(i)
    if (itypehpa(i) == 2) then
      do j = 1,ncoefhor(i)
        itpspl(j,i) = ixlspl(j,i)
        xlatspl(j,i) = xlaspl(j,i)
        xlonspl(j,i) = xlospl(j,i)
        radispl(j,i) = xraspl(j,i)
      enddo
    endif
    hsplfl(i) = hsplfile(i)
  enddo

  numvar = 0
  do i = 1,nmodkern
    string = dskker(i)
    lstr = len_trim(string)
    j = 1
    do while(string(j:j) /= ',' .and. j < lstr)
      j = j+1
    enddo
    ivarkern(i) = 0
    do k = 1,numvar
      if (string(1:j) == varstr(k)(1:j)) then
        ivarkern(i) = k
      endif
    enddo
    if (ivarkern(i) == 0) then
      numvar = numvar+1
      varstr(numvar) = string(1:j)
      ivarkern(i) = numvar
    endif
  enddo

  refmdl = refmodel
  kerstr = kernstri

  ! checks error
  if (ierror /= 0) then
    print *,'Error: reading model in get3dmodel() routine failed with error code ',ierror
    stop 'Error in model s362ani in get3dmodl() routine'
  endif

  end subroutine gt3dmodl

!
!-------------------------------------------------------------------------------------------------
!

  subroutine rd3dmodl(filename,ierror, &
                      nmodkern,nhorpar,ityphpar, &
                      ihorpar,lmaxhor,ncoefhor, &
                      xlaspl,xlospl,xraspl,ixlspl,coef, &
                      hsplfile,refmodel,kernstri,desckern)

  use constants, only: IMAIN,IIN

  use model_s362ani_par, only: maxcoe,maxhpa,maxker

  implicit none

  character(len=128) :: filename
  integer,intent(inout) :: ierror

  integer :: nmodkern,nhorpar

  integer :: ityphpar(maxhpa)
  integer :: ihorpar(maxker)
  integer :: lmaxhor(maxhpa)
  integer :: ncoefhor(maxhpa)

  ! splines
  real(kind=4) :: xlaspl(maxcoe,maxhpa)
  real(kind=4) :: xlospl(maxcoe,maxhpa)
  real(kind=4) :: xraspl(maxcoe,maxhpa)
  integer :: ixlspl(maxcoe,maxhpa)

  real(kind=4) :: coef(maxcoe,maxker)

  character(len=80) :: hsplfile(maxhpa)
  character(len=80) :: refmodel
  character(len=80) :: kernstri
  character(len=40) :: desckern(maxker)

  ! local parameters
  integer :: ncoef,lmax
  integer :: i,ihor,ifst,ilst,ifst1,ios,lstr,idummy
  character(len=128) :: string
  character(len=128) :: substr

  ! opens model file
  open(IIN,file=trim(filename),status='old',action='read',iostat=ios)
  if (ios /= 0) then
    write(IMAIN,*) 'Error opening "', trim(filename), '": ', ios
    call flush_IMAIN()
    call exit_MPI(0, 'Error in model s362ani in rd3dmodl() routine')
  endif

  do while (ios == 0)
  read(IIN,"(a)",iostat=ios) string
  lstr=len_trim(string)
  if (ios == 0) then
    if (string(1:16) == 'REFERENCE MODEL:') then
      substr=string(17:lstr)
      ifst = 1
      ilst=len_trim(substr)
      do while (substr(ifst:ifst) == ' ' .and. ifst < ilst)
        ifst = ifst+1
      enddo
      if (ilst-ifst <= 0) then
        stop 'Error reading model 1'
      else
        refmodel=substr(ifst:ilst)
      endif
    else if (string(1:11) == 'KERNEL SET:') then
      substr=string(12:len_trim(string))
      ifst = 1
      ilst = len_trim(substr)
      do while (substr(ifst:ifst) == ' ' .and. ifst < ilst)
        ifst = ifst+1
      enddo
      if (ilst-ifst <= 0) then
        stop 'Error reading model 2'
      else
        kernstri=substr(ifst:ilst)
      endif
    else if (string(1:25) == 'RADIAL STRUCTURE KERNELS:') then
      substr=string(26:len_trim(string))
      read(substr,*,iostat=ierror) nmodkern
      if (ierror /= 0) then
        stop 'Error reading model 3'
      endif
    else if (string(1:4) == 'DESC' .and. string(9:9) == ':') then
      read(string(5:8),"(i4)") idummy
      substr=string(10:len_trim(string))
      ifst = 1
      ilst = len_trim(substr)
      do while (substr(ifst:ifst) == ' ' .and. ifst < ilst)
        ifst = ifst+1
      enddo
      if (ilst-ifst <= 0) then
        stop 'Error reading model 4'
      else
        desckern(idummy)=substr(ifst:ilst)
      endif
    else if (string(1:29) == 'HORIZONTAL PARAMETERIZATIONS:') then
      substr=string(30:len_trim(string))
      read(substr,*,iostat=ierror) nhorpar
      if (ierror /= 0) then
        stop 'Error reading model 5'
      endif
    else if (string(1:4) == 'HPAR' .and. string(9:9) == ':') then
      read(string(5:8),"(i4)") idummy
      ifst = 10
      ilst = len_trim(string)
      do while (string(ifst:ifst) == ' ' .and. ifst < ilst)
        ifst = ifst+1
      enddo
      if (ilst-ifst <= 0) then
        stop 'Error reading model 6'
      else if (string(ifst:ifst+19) == 'SPHERICAL HARMONICS,') then
        substr=string(20+ifst:len_trim(string))
        read(substr,*) lmax
        ityphpar(idummy) = 1
        lmaxhor(idummy)=lmax
        ncoefhor(idummy)=(lmax+1)**2
      else if (string(ifst:ifst+17) == 'SPHERICAL SPLINES,') then
        ifst1=ifst+18
        ifst=len_trim(string)
        ilst=len_trim(string)
        do while(string(ifst:ifst) /= ',')
          ifst=ifst-1
        enddo
        read(string(ifst+1:ilst),*) ncoef
        substr=string(ifst1:ifst-1)
        do while (string(ifst1:ifst1) == ' ' .and. ifst1 < ifst)
          ifst1=ifst1+1
        enddo
        hsplfile(idummy)=string(ifst1:ifst-1)
        ityphpar(idummy)=2
        lmaxhor(idummy) = 0
        ncoefhor(idummy)=ncoef
        do i = 1,ncoef
          read(IIN,*) ixlspl(i,idummy),xlaspl(i,idummy),xlospl(i,idummy),xraspl(i,idummy)
        enddo
      endif
    else if (string(1:4) == 'STRU' .and. string(9:9) == ':') then
      read(string(5:8),"(i4)") idummy
      substr=string(10:len_trim(string))
      read(substr,*) ihor
      ihorpar(idummy)=ihor
      ncoef=ncoefhor(ihor)
      read(IIN,"(6e12.4)") (coef(i,idummy),i = 1,ncoef)
    endif
  endif
  enddo
  close(IIN)

  end subroutine rd3dmodl

!
!-------------------------------------------------------------------------------------------------
!

  subroutine splcon(xlat,xlon,numcoe,verlat,verlon,verrad,ncon,icon,con)

  use constants, only: DEGREES_TO_RADIANS,RADIANS_TO_DEGREES

  use model_s362ani_par, only: maxver,maxcoe

  implicit none

  real(kind=4),intent(in) :: xlat,xlon

  integer, intent(in) :: numcoe

  real(kind=4), intent(in) :: verlat(numcoe)
  real(kind=4), intent(in) :: verlon(numcoe)
  real(kind=4), intent(in) :: verrad(numcoe)

  integer, intent(out) :: ncon
  integer, intent(out) :: icon(maxver)
  real(kind=4), intent(out) :: con(maxver)

  ! local parameters
  double precision :: rn
  double precision :: dr
  double precision :: ver8
  double precision :: xla8
  double precision, dimension(maxcoe) :: dd

  integer :: iver

  ! safety check
  if (numcoe > maxcoe ) stop 'Error: numcoe > maxver in splcon() routine'

  do iver = 1,numcoe
    if (abs(xlat - verlat(iver)) < 2.*verrad(iver)) then
      ver8 = DEGREES_TO_RADIANS*(verlat(iver))
      xla8 = DEGREES_TO_RADIANS*(xlat)
      dd(iver) = sin(ver8)*sin(xla8) + cos(ver8)*cos(xla8)* cos(DEGREES_TO_RADIANS*(xlon-verlon(iver)))
      dd(iver) = acos(dd(iver)) * RADIANS_TO_DEGREES
    else
      ! acos can never be negative, thus use -1 to mark "invalid"
      dd(iver) = -1.0
    endif
  enddo

  ncon = 0
  do iver = 1,numcoe
    if (dd(iver) >= 0.0 .and. .not. (dd(iver) > (verrad(iver))*2.d0)) then
      ncon = ncon + 1

!! DK DK added this safety test
      if (ncon > maxver) stop 'Error: ncon > maxver in splcon() routine'

      icon(ncon) = iver
      rn = dd(iver)/verrad(iver)
      dr = rn - 1.d0
      if (rn <= 1.d0) then
        con(ncon) = (0.75d0*rn-1.5d0)*(rn**2)+1.d0
      else if (rn > 1.d0) then
        con(ncon) = ((-0.25d0*dr+0.75d0)*dr-0.75d0)*dr+0.25d0
      else
        con(ncon) = 0.0
      endif
    endif
  enddo

  end subroutine splcon

!
!-------------------------------------------------------------------------------------------------
!

! --- evaluate perturbations

  subroutine model_s362ani_subshsv(xcolat,xlon,xrad,dvsh,dvsv,dvph,dvpv)

  use model_s362ani_par
  use constants, only: R_EARTH_KM

  implicit none

  real(kind=4) :: xcolat,xlon,xrad
  real(kind=4) :: dvsh,dvsv,dvph,dvpv

  ! local parameters
  ! --- model evaluation
  integer :: ish ! --- 0 if SV, 1 if SH
  integer :: ieval     ! --- 1 for velocity, 2 for anisotropy
  real(kind=4) :: valu(2)    ! --- valu(1) if S; valu(1)=velo, valu(2)=aniso
  real(kind=4) :: valueval   ! --- used in single evaluation of perturbation
  integer :: isel      ! --- if variable should be included
  real(kind=4) :: depth      ! --- depth
  real(kind=4) :: x,y  ! --- lat lon
  real(kind=4) :: vsh3drel   ! --- relative perturbation
  real(kind=4) :: vsv3drel   ! --- relative perturbation
  ! ---
  integer :: iker,i,ihpa,iver
  integer :: lmax,nylm,numcof
  character(len=40) :: vstr
  integer :: lstr
  integer :: ierror
  ! spherical harmonics
  real(kind=4),dimension((maxl+1)**2,maxhpa) :: ylmcof
  ! splines
  real(kind=4),dimension(maxver,maxhpa) :: conpt
  integer,dimension(maxhpa) :: nconpt
  integer,dimension(maxver,maxhpa) :: iconpt

  real(kind=4), parameter :: r0 = R_EARTH_KM ! 6371.0

  ! initializes
  vsv3drel = 0.0
  vsh3drel = 0.0

  depth = r0 - xrad
  call evradker (depth,kerstr,numker,vercof,vercofd,ierror)
  if (ierror /= 0) stop 'ierror evradker'

  ! loop over sv and sh (sv = 0,sh=1)
  do ish = 0,1

    ! contributing horizontal basis functions at xlat,xlon
    y = 90.0 - xcolat
    x = xlon

    ! sets up coefficients
    do ihpa = 1,numhpa
      if (itypehpa(ihpa) == 1) then
        ! spherical harmonic expansion
        lmax = lmxhpa(ihpa)
        call ylm(y,x,lmax,ylmcof(1,ihpa))

      else if (itypehpa(ihpa) == 2) then
        ! spline setup
        numcof = numcoe(ihpa)
        call splcon(y,x,numcof,xlaspl(1,ihpa), &
                    xlospl(1,ihpa),radspl(1,ihpa), &
                    nconpt(ihpa),iconpt(1,ihpa),conpt(1,ihpa))
      else
        write(*,"('problem 1')")
      endif
    enddo

    ! evaluate 3-D perturbations in velocity and anisotropy

    valu(1) = 0.0 ! --- velocity
    valu(2) = 0.0 ! --- anisotropy

    do ieval = 1,2
      valueval = 0.0

      do iker = 1,numker
        isel = 0
        lstr = len_trim(varstr(ivarkern(iker)))
        vstr = (varstr(ivarkern(iker)))
        if (ieval == 1) then
          if (vstr(1:lstr) == 'UM (SH+SV)*0.5,' .or. vstr(1:lstr) == 'LM (SH+SV)*0.5,' .or. &
             vstr(1:lstr) == 'EA (SH+SV)*0.5,') then
            isel = 1
          endif
        else if (ieval == 2) then
          if (vstr(1:lstr) == 'UM SH-SV,' .or. vstr(1:lstr) == 'LM SH-SV,' .or. &
             vstr(1:lstr) == 'EA SH-SV,') then
            isel = 1
          endif
        endif

        if (isel == 1) then
          if (vercof(iker) /= 0.0) then
            if (itypehpa(ihpakern(iker)) == 1) then
              ! spherical harmonics
              ihpa = ihpakern(iker)
              nylm = (lmxhpa(ihpakern(iker))+1)**2
              do i = 1,nylm
                valueval = valueval + vercof(iker)*ylmcof(i,ihpa)*coe(i,iker)
              enddo

            else if (itypehpa(ihpakern(iker)) == 2) then
              ! splines
              ihpa = ihpakern(iker)
              do i = 1,nconpt(ihpa)
                iver = iconpt(i,ihpa)
                valueval = valueval + vercof(iker)*conpt(i,ihpa)*coe(iver,iker)
              enddo

            else
              stop 'problem 2'
            endif ! --- itypehpa
          endif ! --- vercof(iker) /= 0.
        endif ! --- isel == 1
      enddo ! --- end of do iker = 1,numker

      valu(ieval) = valueval

    enddo ! --- ieval

    ! evaluate perturbations in vsh and vsv
    if (ish == 1) then
      vsh3drel = valu(1) + 0.5*valu(2)
    else if (ish == 0) then
      vsv3drel = valu(1) - 0.5*valu(2)
    else
      stop 'something is wrong in model_s362ani_subshsv'
    endif

  enddo ! --- by ish

  ! evaluate perturbations
  dvsh = vsh3drel
  dvsv = vsv3drel
  dvph = 0.55*dvsh    ! scaling used in the inversion
  dvpv = 0.55*dvsv    ! scaling used in the inversion

  end subroutine model_s362ani_subshsv

!
!-------------------------------------------------------------------------------------------------
!

! --- evaluate depressions of the 410- and 650-km discontinuities in km

  subroutine model_s362ani_subtopo(xcolat,xlon,topo410,topo650)

  use model_s362ani_par, only: numhpa,itypehpa,lmxhpa,maxcoe,numcoe,coe, &
    xlaspl,xlospl,radspl, &
    numker,ihpakern,ivarkern,varstr, &
    maxl,maxhpa,maxver,maxhpa

  implicit none

  real(kind=4),intent(in) :: xcolat,xlon
  real(kind=4),intent(out) :: topo410,topo650

  ! --- model evaluation
  integer :: ieval     ! --- 1 for velocity, 2 for anisotropy
  real(kind=4) :: valu(2)    ! --- valu(1) if S; valu(1)=velo, valu(2)=aniso
  real(kind=4) :: valueval   ! --- used in single evaluation of perturbation
  integer :: isel      ! --- if variable should be included
  real(kind=4) :: x,y  ! --- lat lon
  ! ---
  integer :: iker,i,ihpa,iver
  integer :: lmax,nylm,numcof
  character(len=40) :: vstr
  integer :: lstr
  ! spherical harmonics
  real(kind=4),dimension((maxl+1)**2,maxhpa) :: ylmcof
  ! splines
  real(kind=4),dimension(maxver,maxhpa) :: conpt
  integer,dimension(maxhpa) :: nconpt
  integer,dimension(maxver,maxhpa) :: iconpt

  ! contributing horizontal basis functions at xlat,xlon
  y = 90.0 - xcolat
  x = xlon

  ! sets up coefficients
  do ihpa = 1,numhpa
    if (itypehpa(ihpa) == 1) then
      ! spherical harmonics
      lmax = lmxhpa(ihpa)
      call ylm(y,x,lmax,ylmcof(1,ihpa))

    else if (itypehpa(ihpa) == 2) then
      ! splines
      numcof = numcoe(ihpa)
      call splcon(y,x,numcof,xlaspl(1,ihpa), &
                  xlospl(1,ihpa),radspl(1,ihpa), &
                  nconpt(ihpa),iconpt(1,ihpa),conpt(1,ihpa))

    else
      write(*,"('problem 1')")
    endif
  enddo

  ! evaluate topography (depression) in km

  valu(1) = 0.0 ! --- 410
  valu(2) = 0.0 ! --- 650

  ! evaluates topography perturbation for 410/670 discontinuities
  do ieval = 1,2
    valueval = 0.0

    do iker = 1,numker
      isel = 0
      lstr = len_trim(varstr(ivarkern(iker)))
      vstr = (varstr(ivarkern(iker)))

      ! selects discontinuity
      if (ieval == 1) then
        ! treats discontinuity 410
        if (vstr(1:lstr) == 'Topo 400,') then
          isel = 1
        endif
      else if (ieval == 2) then
        ! treats discontinuity 670
        if (vstr(1:lstr) == 'Topo 670,') then
          isel = 1
        endif
      endif

      ! evaluates value
      if (isel == 1) then
        if (itypehpa(ihpakern(iker)) == 1) then
          ! spherical harmonics
          ihpa = ihpakern(iker)
          nylm = (lmxhpa(ihpakern(iker))+1)**2
          do i = 1,nylm
            valueval = valueval + ylmcof(i,ihpa)*coe(i,iker)
          enddo

        else if (itypehpa(ihpakern(iker)) == 2) then
          ! splines
          ihpa = ihpakern(iker)
          do i = 1,nconpt(ihpa)
            iver = iconpt(i,ihpa)
            valueval = valueval + conpt(i,ihpa)*coe(iver,iker)
          enddo

        else
          stop 'problem 2'
        endif ! --- itypehpa
      endif ! --- isel == 1
    enddo ! --- end of do iker = 1,numker

    valu(ieval) = valueval

  enddo ! --- ieval

  topo410 = valu(1)
  topo650 = valu(2)

  end subroutine model_s362ani_subtopo

!
!-------------------------------------------------------------------------------------------------
!

  subroutine vbspl(x,np,xarr,splcon,splcond)
!
!---- this subroutine returns the spline contributions at a particular value of x
!
  implicit none

  integer :: np

  real(kind=4) :: xarr(np),x
  real(kind=4) :: splcon(np)
  real(kind=4) :: splcond(np)

  real(kind=4) :: r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13
  real(kind=4) :: r1d,r2d,r3d,r4d,r5d,r6d,r7d,r8d,r9d,r10d,r11d,r12d,r13d,val,vald

  real(kind=4) :: rr1,rr2,rr3,rr4,rr5,rr6,rr7,rr8,rr9,rr10,rr11,rr12
  real(kind=4) :: rr1d,rr2d,rr3d,rr4d,rr5d,rr6d,rr7d,rr8d,rr9d,rr10d,rr11d,rr12d

  integer :: iflag,interval,ik,ib

!
!---- iflag = 1, second derivative is 0 at end points
!---- iflag = 0, first derivative is 0 at end points
!
  iflag = 1
!
!---- first, find out within which interval x falls
!
  interval = 0
  ik = 1
  do while(interval == 0 .and. ik < np)
    ik = ik+1
    if (x >= xarr(ik-1) .and. x <= xarr(ik)) interval=ik-1
  enddo
  if (x > xarr(np)) then
    interval = np
  endif

  do ib = 1,np

  val = 0.0
  vald = 0.0

  if (ib == 1) then

    r1 = (x-xarr(1))/(xarr(2)-xarr(1))
    r2 = (xarr(3)-x)/(xarr(3)-xarr(1))
    r4 = (xarr(2)-x)/(xarr(2)-xarr(1))
    r5 = (x-xarr(1))/(xarr(2)-xarr(1))
    r6 = (xarr(3)-x)/(xarr(3)-xarr(1))
    r10 = (xarr(2)-x)/(xarr(2)-xarr(1))
    r11 = (x-xarr(1))  /(xarr(2)-xarr(1))
    r12 = (xarr(3)-x)/(xarr(3)-xarr(2))
    r13 = (xarr(2)-x)/(xarr(2)-xarr(1))

    r1d = 1.0/(xarr(2)-xarr(1))
    r2d = -1.0/(xarr(3)-xarr(1))
    r4d = -1.0/(xarr(2)-xarr(1))
    r5d = 1.0/(xarr(2)-xarr(1))
    r6d = -1.0/(xarr(3)-xarr(1))
    r10d = -1.0/(xarr(2)-xarr(1))
    r11d = 1.0/(xarr(2)-xarr(1))
    r12d = -1.0/(xarr(3)-xarr(2))
    r13d = -1.0/(xarr(2)-xarr(1))

    if (interval == ib .or. interval == 0) then

      if (iflag == 0) then
        val = r1*r4*r10 + r2*r5*r10 + r2*r6*r11 +r13**3
        vald = r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
        vald = vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
        vald = vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
        vald = vald+3.*r13d*r13**2
      else if (iflag == 1) then
        val = 0.6667*(r1*r4*r10 + r2*r5*r10 + r2*r6*r11 + 1.5*r13**3)
        vald = r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
        vald = vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
        vald = vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
        vald = vald+4.5*r13d*r13**2
        vald = 0.6667*vald
      endif
    else if (interval == ib+1) then
      if (iflag == 0) then
        val = r2*r6*r12
        vald = r2d*r6*r12+r2*r6d*r12+r2*r6*r12d
      else if (iflag == 1) then
        val = 0.6667*r2*r6*r12
        vald = 0.6667*(r2d*r6*r12+r2*r6d*r12+r2*r6*r12d)
      endif
    else
      val = 0.0
    endif

  else if (ib == 2) then

    rr1 = (x-xarr(1))/(xarr(2)-xarr(1))
    rr2 = (xarr(3)-x)/(xarr(3)-xarr(1))
    rr4 = (xarr(2)-x)/(xarr(2)-xarr(1))
    rr5 = (x-xarr(1))/(xarr(2)-xarr(1))
    rr6 = (xarr(3)-x)/(xarr(3)-xarr(1))
    rr10 = (xarr(2)-x)/(xarr(2)-xarr(1))
    rr11 = (x-xarr(1))  /(xarr(2)-xarr(1))
    rr12 = (xarr(3)-x)/(xarr(3)-xarr(2))

    rr1d = 1.0/(xarr(2)-xarr(1))
    rr2d = -1.0/(xarr(3)-xarr(1))
    rr4d = -1.0/(xarr(2)-xarr(1))
    rr5d = 1.0/(xarr(2)-xarr(1))
    rr6d = -1.0/(xarr(3)-xarr(1))
    rr10d = -1.0/(xarr(2)-xarr(1))
    rr11d = 1.0/(xarr(2)-xarr(1))
    rr12d = -1.0/(xarr(3)-xarr(2))

    r1 = (x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
    r2 = (xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib-1))
    r3 = (x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
    r4 = (xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
    r5 = (x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
    r6 = (xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib))
    r8 = (xarr(ib)-x)/  (xarr(ib)-xarr(ib-1))
    r9 = (x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
    r10 = (xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
    r11 = (x-xarr(ib))  /(xarr(ib+1)-xarr(ib))
    r12 = (xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib+1))

    r1d = 1.0/(xarr(ib+1)-xarr(ib-1))
    r2d = -1.0/(xarr(ib+2)-xarr(ib-1))
    r3d = 1.0/(xarr(ib)-xarr(ib-1))
    r4d = -1.0/(xarr(ib+1)-xarr(ib-1))
    r5d = 1.0/(xarr(ib+1)-xarr(ib-1))
    r6d = -1.0/(xarr(ib+2)-xarr(ib))
    r8d = -1.0/  (xarr(ib)-xarr(ib-1))
    r9d = 1.0/(xarr(ib)-xarr(ib-1))
    r10d = -1.0/(xarr(ib+1)-xarr(ib))
    r11d = 1.0/(xarr(ib+1)-xarr(ib))
    r12d = -1.0/(xarr(ib+2)-xarr(ib+1))

    if (interval == ib-1 .or. interval == 0) then
      val = r1*r3*r8 + r1*r4*r9 + r2*r5*r9
      vald = r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
      vald = vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
      vald = vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
      if (iflag == 1) then
        val = val+0.3333*(rr1*rr4*rr10 + rr2*rr5*rr10 + rr2*rr6*rr11)
        vald = vald+0.3333*(rr1d*rr4*rr10+rr1*rr4d*rr10+ rr1*rr4*rr10d)
        vald = vald+0.3333*(rr2d*rr5*rr10+rr2*rr5d*rr10+ rr2*rr5*rr10d)
        vald = vald+0.3333*(rr2d*rr6*rr11+rr2*rr6d*rr11+ rr2*rr6*rr11d)
      endif
    else if (interval == ib) then
      val = r1*r4*r10 + r2*r5*r10 + r2*r6*r11
      vald = r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
      vald = vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
      vald = vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
      if (iflag == 1) then
        val = val+0.3333*rr2*rr6*rr12
        vald = vald+0.3333*(rr2d*rr6*rr12+rr2*rr6d*rr12+ rr2*rr6*rr12d)
      endif
    else if (interval == ib+1) then
      val = r2*r6*r12
      vald = r2d*r6*r12+r2*r6d*r12+r2*r6*r12d
    else
      val = 0.0
    endif

  else if (ib == np-1) then

    rr1 = (x-xarr(np-2))/(xarr(np)-xarr(np-2))
    rr2 = (xarr(np)-x)/(xarr(np)-xarr(np-1))
    rr3 = (x-xarr(np-2))/(xarr(np)-xarr(np-2))
    rr4 = (xarr(np)-x)/(xarr(np)-xarr(np-1))
    rr5 = (x-xarr(np-1))/(xarr(np)-xarr(np-1))
    rr7 = (x-xarr(np-2))/(xarr(np-1)-xarr(np-2))
    rr8 = (xarr(np)-x)/  (xarr(np)-xarr(np-1))
    rr9 = (x-xarr(np-1))/(xarr(np)-xarr(np-1))

    rr1d = 1.0/(xarr(np)-xarr(np-2))
    rr2d = -1.0/(xarr(np)-xarr(np-1))
    rr3d = 1.0/(xarr(np)-xarr(np-2))
    rr4d = -1.0/(xarr(np)-xarr(np-1))
    rr5d = 1.0/(xarr(np)-xarr(np-1))
    rr7d = 1.0/(xarr(np-1)-xarr(np-2))
    rr8d = -1.0/  (xarr(np)-xarr(np-1))
    rr9d = 1.0/(xarr(np)-xarr(np-1))

    r1 = (x-xarr(ib-2))/(xarr(ib+1)-xarr(ib-2))
    r2 = (xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
    r3 = (x-xarr(ib-2))/(xarr(ib)-xarr(ib-2))
    r4 = (xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
    r5 = (x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
    r6 = (xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
    r7 = (x-xarr(ib-2))/(xarr(ib-1)-xarr(ib-2))
    r8 = (xarr(ib)-x)/  (xarr(ib)-xarr(ib-1))
    r9 = (x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
    r10 = (xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
    r11 = (x-xarr(ib))  /(xarr(ib+1)-xarr(ib))

    r1d = 1.0/(xarr(ib+1)-xarr(ib-2))
    r2d = -1.0/(xarr(ib+1)-xarr(ib-1))
    r3d = 1.0/(xarr(ib)-xarr(ib-2))
    r4d = -1.0/(xarr(ib+1)-xarr(ib-1))
    r5d = 1.0/(xarr(ib+1)-xarr(ib-1))
    r6d = -1.0/(xarr(ib+1)-xarr(ib))
    r7d = 1.0/(xarr(ib-1)-xarr(ib-2))
    r8d = -1.0/(xarr(ib)-xarr(ib-1))
    r9d = 1.0/(xarr(ib)-xarr(ib-1))
    r10d = -1.0/(xarr(ib+1)-xarr(ib))
    r11d = 1.0/(xarr(ib+1)-xarr(ib))

    if (interval == ib-2) then
      val = r1*r3*r7
      vald = r1d*r3*r7+r1*r3d*r7+r1*r3*r7d
    else if (interval == ib-1) then
      val = r1*r3*r8 + r1*r4*r9 + r2*r5*r9
      vald = r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
      vald = vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
      vald = vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
      if (iflag == 1) then
        val = val+0.3333*rr1*rr3*rr7
        vald = vald+0.3333*(rr1d*rr3*rr7+rr1*rr3d*rr7+ rr1*rr3*rr7d)
      endif
    else if (interval == ib .or. interval == np) then
      val = r1*r4*r10 + r2*r5*r10 + r2*r6*r11
      vald = r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
      vald = vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
      vald = vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
      if (iflag == 1) then
        val = val+0.3333*(rr1*rr3*rr8 + rr1*rr4*rr9 + rr2*rr5*rr9)
        vald = vald+0.3333*(rr1d*rr3*rr8+rr1*rr3d*rr8+ rr1*rr3*rr8d)
        vald = vald+0.3333*(rr1d*rr4*rr9+rr1*rr4d*rr9+ rr1*rr4*rr9d)
        vald = vald+0.3333*(rr2d*rr5*rr9+rr2*rr5d*rr9+ rr2*rr5*rr9d)
      endif
    else
      val = 0.0
    endif

  else if (ib == np) then

    r1 = (x-xarr(np-2))/(xarr(np)-xarr(np-2))
    r2 = (xarr(np)-x)/(xarr(np)-xarr(np-1))
    r3 = (x-xarr(np-2))/(xarr(np)-xarr(np-2))
    r4 = (xarr(np)-x)/(xarr(np)-xarr(np-1))
    r5 = (x-xarr(np-1))/(xarr(np)-xarr(np-1))
    r7 = (x-xarr(np-2))/(xarr(np-1)-xarr(np-2))
    r8 = (xarr(np)-x)/  (xarr(np)-xarr(np-1))
    r9 = (x-xarr(np-1))/(xarr(np)-xarr(np-1))
    r13 = (x-xarr(np-1))/(xarr(np)-xarr(np-1))

    r1d = 1.0/(xarr(np)-xarr(np-2))
    r2d = -1.0/(xarr(np)-xarr(np-1))
    r3d = 1.0/(xarr(np)-xarr(np-2))
    r4d = -1.0/(xarr(np)-xarr(np-1))
    r5d = 1.0/(xarr(np)-xarr(np-1))
    r7d = 1.0/(xarr(np-1)-xarr(np-2))
    r8d = -1.0/  (xarr(np)-xarr(np-1))
    r9d = 1.0/(xarr(np)-xarr(np-1))
    r13d = 1.0/(xarr(np)-xarr(np-1))

    if (interval == np-2) then
      if (iflag == 0) then
        val = r1*r3*r7
        vald = r1d*r3*r7+r1*r3d*r7+r1*r3*r7d
      else if (iflag == 1) then
        val = 0.6667*r1*r3*r7
        vald = 0.6667*(r1d*r3*r7+r1*r3d*r7+r1*r3*r7d)
      endif
    else if (interval == np-1 .or. interval == np) then
      if (iflag == 0) then
        val = r1*r3*r8 + r1*r4*r9 + r2*r5*r9 + r13**3
        vald = r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
        vald = vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
        vald = vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
        vald = vald+3.*r13d*r13**2
      else if (iflag == 1) then
        val = 0.6667*(r1*r3*r8 + r1*r4*r9 + r2*r5*r9 + 1.5*r13**3)
        vald = r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
        vald = vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
        vald = vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
        vald = vald+4.5*r13d*r13**2
        vald = 0.6667*vald
      endif
    else
      val = 0.0
    endif
  else

    r1 = (x-xarr(ib-2))/(xarr(ib+1)-xarr(ib-2))
    r2 = (xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib-1))
    r3 = (x-xarr(ib-2))/(xarr(ib)-xarr(ib-2))
    r4 = (xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
    r5 = (x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
    r6 = (xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib))
    r7 = (x-xarr(ib-2))/(xarr(ib-1)-xarr(ib-2))
    r8 = (xarr(ib)-x)/  (xarr(ib)-xarr(ib-1))
    r9 = (x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
    r10 = (xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
    r11 = (x-xarr(ib))  /(xarr(ib+1)-xarr(ib))
    r12 = (xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib+1))

    r1d = 1.0/(xarr(ib+1)-xarr(ib-2))
    r2d = -1.0/(xarr(ib+2)-xarr(ib-1))
    r3d = 1.0/(xarr(ib)-xarr(ib-2))
    r4d = -1.0/(xarr(ib+1)-xarr(ib-1))
    r5d = 1.0/(xarr(ib+1)-xarr(ib-1))
    r6d = -1.0/(xarr(ib+2)-xarr(ib))
    r7d = 1.0/(xarr(ib-1)-xarr(ib-2))
    r8d = -1.0/  (xarr(ib)-xarr(ib-1))
    r9d = 1.0/(xarr(ib)-xarr(ib-1))
    r10d = -1.0/(xarr(ib+1)-xarr(ib))
    r11d = 1.0/(xarr(ib+1)-xarr(ib))
    r12d = -1.0/(xarr(ib+2)-xarr(ib+1))

    if (interval == ib-2) then
      val = r1*r3*r7
      vald = r1d*r3*r7+r1*r3d*r7+r1*r3*r7d
    else if (interval == ib-1) then
      val = r1*r3*r8 + r1*r4*r9 + r2*r5*r9
      vald = r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
      vald = vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
      vald = vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
    else if (interval == ib) then
      val = r1*r4*r10 + r2*r5*r10 + r2*r6*r11
      vald = r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
      vald = vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
      vald = vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
    else if (interval == ib+1) then
      val = r2*r6*r12
      vald = r2d*r6*r12+r2*r6d*r12+r2*r6*r12d
    else
      val = 0.0
    endif
  endif

  splcon(ib) = val
  splcond(ib) = vald

  enddo

  end subroutine vbspl

