!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            December 2010
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


  subroutine model_sea99_s_broadcast(myrank,SEA99M_V)

! standard routine to setup model

  implicit none

  include "constants.h"
  ! standard include of the MPI library
  include 'mpif.h'

  ! model_sea99_s_variables
  type model_sea99_s_variables
    sequence
    double precision :: sea99_vs(100,100,100)
    double precision :: sea99_depth(100)
    double precision :: sea99_ddeg
    double precision :: alatmin
    double precision :: alatmax
    double precision :: alonmin
    double precision :: alonmax
    integer :: sea99_ndep
    integer :: sea99_nlat
    integer :: sea99_nlon
    integer :: dummy_pad ! padding 4 bytes to align the structure
 end type model_sea99_s_variables

  type (model_sea99_s_variables) SEA99M_V
  ! model_sea99_s_variables

  integer :: myrank
  integer :: ier

  if(myrank == 0) call read_sea99_s_model(SEA99M_V)

  ! broadcast the information read on the master to the nodes
  ! SEA99M_V
  call MPI_BCAST(SEA99M_V%sea99_ndep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(SEA99M_V%sea99_nlat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(SEA99M_V%sea99_nlon,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(SEA99M_V%sea99_ddeg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(SEA99M_V%alatmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(SEA99M_V%alatmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(SEA99M_V%alonmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(SEA99M_V%alonmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(SEA99M_V%sea99_vs,100*100*100,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(SEA99M_V%sea99_depth,100,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

  end subroutine model_sea99_s_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_sea99_s_model(SEA99M_V)

  implicit none

  include "constants.h"

  ! model_sea99_s_variables
  type model_sea99_s_variables
    sequence
    double precision :: sea99_vs(100,100,100)
    double precision :: sea99_depth(100)
    double precision :: sea99_ddeg
    double precision :: alatmin
    double precision :: alatmax
    double precision :: alonmin
    double precision :: alonmax
    integer :: sea99_ndep
    integer :: sea99_nlat
    integer :: sea99_nlon
    integer :: dummy_pad ! padding 4 bytes to align the structure
 end type model_sea99_s_variables

  type (model_sea99_s_variables) SEA99M_V
  ! model_sea99_s_variables

  integer :: i,ia,io,j

!----------------------- choose input file:  ------------------
! relative anomaly


  open(1,file='DATA/Lebedev_sea99/sea99_dvsvs')

!----------------------- read input file:  ------------------

  do i = 1, 6
     read(1,*)
  enddo
  read(1,*) SEA99M_V%sea99_ndep
  read(1,*) (SEA99M_V%sea99_depth(i), i = 1, SEA99M_V%sea99_ndep)
  read(1,*)
  read(1,*) SEA99M_V%alatmin, SEA99M_V%alatmax
  read(1,*) SEA99M_V%alonmin, SEA99M_V%alonmax
  read(1,*) SEA99M_V%sea99_ddeg,SEA99M_V%sea99_nlat,SEA99M_V%sea99_nlon
  if (SEA99M_V%sea99_nlat /= nint((SEA99M_V%alatmax-SEA99M_V%alatmin)/SEA99M_V%sea99_ddeg)+1) then
     stop 'alatmin,alatmax,sea99_nlat'
  endif
  if (SEA99M_V%sea99_nlon /= nint((SEA99M_V%alonmax-SEA99M_V%alonmin)/SEA99M_V%sea99_ddeg)+1) then
     stop 'alonmin,alonmax,sea99_nlon'
  endif
  read(1,*)
  do j = 1, SEA99M_V%sea99_ndep
     do ia = 1, SEA99M_V%sea99_nlat
        read (1,*) (SEA99M_V%sea99_vs(ia,io,j), io = 1, SEA99M_V%sea99_nlon)
     enddo
  enddo

  end subroutine read_sea99_s_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_sea99_s(radius,theta,phi,dvs,SEA99M_V)

! returns Vs perturbation (dvs) for given position r/theta/phi

  implicit none

  include "constants.h"

  ! model_sea99_s_variables
  type model_sea99_s_variables
    sequence
    double precision :: sea99_vs(100,100,100)
    double precision :: sea99_depth(100)
    double precision :: sea99_ddeg
    double precision :: alatmin
    double precision :: alatmax
    double precision :: alonmin
    double precision :: alonmax
    integer :: sea99_ndep
    integer :: sea99_nlat
    integer :: sea99_nlon
    integer :: dummy_pad ! padding 4 bytes to align the structure
 end type model_sea99_s_variables

  type (model_sea99_s_variables) SEA99M_V
  ! model_sea99_s_variables

  integer :: id1,i,ilat,ilon
  double precision :: alat1,alon1,radius,theta,phi,dvs
  double precision :: xxx,yyy,dep,pla,plo,xd1,dd1,dd2,ddd(2)

  ! initializes
  dvs = 0.d0

  id1 = 0
  xd1 = 0

  !----------------------- depth in the model ------------------
  dep=R_EARTH_KM*(R_UNIT_SPHERE - radius)
  if (dep .le. SEA99M_V%sea99_depth(1)) then
     id1 = 1
     xd1 = 0
  else if (dep .ge. SEA99M_V%sea99_depth(SEA99M_V%sea99_ndep)) then
     id1 = SEA99M_V%sea99_ndep
     xd1 = 0
  else
     do i = 2, SEA99M_V%sea99_ndep
        if (dep .le. SEA99M_V%sea99_depth(i)) then
           id1 = i-1
           xd1 = (dep-SEA99M_V%sea99_depth(i-1)) / (SEA99M_V%sea99_depth(i) - SEA99M_V%sea99_depth(i-1))
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
  if( pla < SEA99M_V%alatmin .or. pla > SEA99M_V%alatmax &
    .or. plo < SEA99M_V%alonmin .or. plo > SEA99M_V%alonmax ) return

  ! array indices
  ilat = int((pla - SEA99M_V%alatmin)/SEA99M_V%sea99_ddeg) + 1
  ilon = int((plo - SEA99M_V%alonmin)/SEA99M_V%sea99_ddeg) + 1
  alat1 = SEA99M_V%alatmin + (ilat-1)*SEA99M_V%sea99_ddeg
  alon1 = SEA99M_V%alonmin + (ilon-1)*SEA99M_V%sea99_ddeg

  do i = 1, 2
     xxx = (pla-alat1)/SEA99M_V%sea99_ddeg
     yyy = SEA99M_V%sea99_vs(ilat+1,ilon,id1+i-1)-SEA99M_V%sea99_vs(ilat,ilon,id1+i-1)
     dd1 = SEA99M_V%sea99_vs(ilat,ilon,id1+i-1) + yyy*xxx
     yyy = SEA99M_V%sea99_vs(ilat+1,ilon+1,id1+i-1)-SEA99M_V%sea99_vs(ilat,ilon+1,id1+i-1)
     dd2 = SEA99M_V%sea99_vs(ilat,ilon+1,id1+i-1) + yyy*xxx
     xxx = (plo-alon1)/SEA99M_V%sea99_ddeg
     yyy = dd2 - dd1
     ddd(i) = dd1 + yyy*xxx
  enddo
  dvs = ddd(1) + (ddd(2)-ddd(1)) * xd1

  ! checks perturbation
  if(dvs > 1.d0) dvs = 0.0d0

  end subroutine model_sea99_s


