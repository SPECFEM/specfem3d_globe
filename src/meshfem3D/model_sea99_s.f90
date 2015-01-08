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
! the Free Software Foundation; either version 2 of the License, or
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
! SEA 99 model
!
! contains relative Vs anomalies  dVs/Vs from
! SV-velocity model for SE Asia - W Pacific.
!
! defined for:
! -20.00   45.00 -- min, max latitude
!  95.00  160.00 -- min, max longitude
! and depths between 6 km to 860 km
!
! computed by Lebedev and Nolet in 1999, to come out in JGR in 2003.
! reference period: 50 s.
!--------------------------------------------------------------------------------------------------

  module model_sea99_s_par

  double precision,dimension(:,:,:), allocatable :: sea99_vs
  double precision,dimension(:), allocatable :: sea99_depth
  double precision :: sea99_ddeg
  double precision :: alatmin
  double precision :: alatmax
  double precision :: alonmin
  double precision :: alonmax
  integer :: sea99_ndep
  integer :: sea99_nlat
  integer :: sea99_nlon

  end module model_sea99_s_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_sea99_s_broadcast(myrank)

! standard routine to setup model

  use constants
  use model_sea99_s_par

  implicit none

  integer :: myrank

  integer :: ier

  ! allocates model arrays
  allocate(sea99_vs(100,100,100), &
           sea99_depth(100), &
           stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating sea99 arrays')

  ! master proc reads in values
  if (myrank == 0) call read_sea99_s_model()

  ! broadcast the information read on the master to the nodes
  call bcast_all_singlei(sea99_ndep)
  call bcast_all_singlei(sea99_nlat)
  call bcast_all_singlei(sea99_nlon)

  call bcast_all_singledp(sea99_ddeg)
  call bcast_all_singledp(alatmin)
  call bcast_all_singledp(alatmax)
  call bcast_all_singledp(alonmin)
  call bcast_all_singledp(alonmax)

  call bcast_all_dp(sea99_vs,100*100*100)
  call bcast_all_dp(sea99_depth,100)

  end subroutine model_sea99_s_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_sea99_s_model()

  use constants
  use model_sea99_s_par

  implicit none

  integer :: i,ia,io,j,ier

!----------------------- choose input file:  ------------------
! relative anomaly


  open(IIN,file='DATA/Lebedev_sea99/sea99_dvsvs',status='old',action='read',iostat=ier)
  if (ier /= 0 ) call exit_MPI(0,'Error opening file sea99_dvsvs')

!----------------------- read input file:  ------------------

  do i = 1, 6
     read(IIN,*)
  enddo
  read(IIN,*) sea99_ndep
  read(IIN,*) (sea99_depth(i), i = 1, sea99_ndep)
  read(IIN,*)
  read(IIN,*) alatmin, alatmax
  read(IIN,*) alonmin, alonmax
  read(IIN,*) sea99_ddeg,sea99_nlat,sea99_nlon
  if (sea99_nlat /= nint((alatmax-alatmin)/sea99_ddeg)+1) then
     stop 'alatmin,alatmax,sea99_nlat'
  endif
  if (sea99_nlon /= nint((alonmax-alonmin)/sea99_ddeg)+1) then
     stop 'alonmin,alonmax,sea99_nlon'
  endif
  read(IIN,*)
  do j = 1, sea99_ndep
     do ia = 1, sea99_nlat
        read(IIN,*) (sea99_vs(ia,io,j), io = 1, sea99_nlon)
     enddo
  enddo
  close(IIN)

  end subroutine read_sea99_s_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_sea99_s(radius,theta,phi,dvs)

! returns Vs perturbation (dvs) for given position r/theta/phi

  use constants
  use model_sea99_s_par

  implicit none

  integer :: id1,i,ilat,ilon
  double precision :: alat1,alon1,radius,theta,phi,dvs
  double precision :: xxx,yyy,dep,pla,plo,xd1,dd1,dd2,ddd(2)

  ! initializes
  dvs = 0.d0

  id1 = 0
  xd1 = 0

  !----------------------- depth in the model ------------------
  dep=R_EARTH_KM*(R_UNIT_SPHERE - radius)
  if (dep <= sea99_depth(1)) then
     id1 = 1
     xd1 = 0
  else if (dep >= sea99_depth(sea99_ndep)) then
     id1 = sea99_ndep
     xd1 = 0
  else
     do i = 2, sea99_ndep
        if (dep <= sea99_depth(i)) then
           id1 = i-1
           xd1 = (dep-sea99_depth(i-1)) / (sea99_depth(i) - sea99_depth(i-1))
           exit
        endif
     enddo
  endif

  !----------------------- value at a point ---------------------
  !----- approximate interpolation, OK for the (dense) 1-degree sampling ------

  ! latitude / longitude in degree
  pla = 90.0d0 - theta/DEGREES_TO_RADIANS
  plo = phi/DEGREES_TO_RADIANS

  ! model defined for:
  ! -20.00   45.00 -- min, max latitude
  !  95.00  160.00 -- min, max longitude
  ! checks range
  if (pla < alatmin .or. pla > alatmax &
    .or. plo < alonmin .or. plo > alonmax ) return

  ! array indices
  ilat = int((pla - alatmin)/sea99_ddeg) + 1
  ilon = int((plo - alonmin)/sea99_ddeg) + 1
  alat1 = alatmin + (ilat-1)*sea99_ddeg
  alon1 = alonmin + (ilon-1)*sea99_ddeg

  do i = 1, 2
     xxx = (pla-alat1)/sea99_ddeg
     yyy = sea99_vs(ilat+1,ilon,id1+i-1)-sea99_vs(ilat,ilon,id1+i-1)
     dd1 = sea99_vs(ilat,ilon,id1+i-1) + yyy*xxx
     yyy = sea99_vs(ilat+1,ilon+1,id1+i-1)-sea99_vs(ilat,ilon+1,id1+i-1)
     dd2 = sea99_vs(ilat,ilon+1,id1+i-1) + yyy*xxx
     xxx = (plo-alon1)/sea99_ddeg
     yyy = dd2 - dd1
     ddd(i) = dd1 + yyy*xxx
  enddo
  dvs = ddd(1) + (ddd(2)-ddd(1)) * xd1

  ! checks perturbation
  if (dvs > 1.d0) dvs = 0.0d0

  end subroutine model_sea99_s


