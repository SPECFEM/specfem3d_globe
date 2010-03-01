!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            March 2010
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

  subroutine get_ellipticity(xelm,yelm,zelm,nspl,rspl,espl,espl2)

  implicit none

  include "constants.h"

  integer nspl
  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)
  double precision rspl(NR),espl(NR),espl2(NR)

  integer ia

  double precision ell
  double precision r,theta,phi,factor
  double precision cost,p20

  do ia=1,NGNOD

  call xyz_2_rthetaphi_dble(xelm(ia),yelm(ia),zelm(ia),r,theta,phi)

  cost=dcos(theta)
  p20=0.5d0*(3.0d0*cost*cost-1.0d0)

! get ellipticity using spline evaluation
  call spline_evaluation(rspl,espl,espl2,nspl,r,ell)

  factor=ONE-(TWO/3.0d0)*ell*p20

  xelm(ia)=xelm(ia)*factor
  yelm(ia)=yelm(ia)*factor
  zelm(ia)=zelm(ia)*factor

  enddo

!
!-------------------------------------------------------------------------------------------------
!

  end subroutine get_ellipticity

  !> Hejun 
  ! get ellipticity according to GLL points
  ! JAN08, 2010
  subroutine get_ellipticity_gll(xstore,ystore,zstore,ispec,nspec,nspl,rspl,espl,espl2)

  implicit none

  include "constants.h"

  integer nspl
  integer::ispec,nspec
  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec):: xstore,ystore,zstore
  double precision rspl(NR),espl(NR),espl2(NR)

  integer i,j,k

  double precision ell
  double precision r,theta,phi,factor
  double precision cost,p20

  do k = 1,NGLLZ
     do j = 1,NGLLY
        do i = 1,NGLLX

           call xyz_2_rthetaphi_dble(xstore(i,j,k,ispec),ystore(i,j,k,ispec),zstore(i,j,k,ispec),r,theta,phi)

           cost=dcos(theta)
           p20=0.5d0*(3.0d0*cost*cost-1.0d0)

           ! get ellipticity using spline evaluation
           call spline_evaluation(rspl,espl,espl2,nspl,r,ell)

           factor=ONE-(TWO/3.0d0)*ell*p20

           xstore(i,j,k,ispec)=xstore(i,j,k,ispec)*factor
           ystore(i,j,k,ispec)=ystore(i,j,k,ispec)*factor
           zstore(i,j,k,ispec)=zstore(i,j,k,ispec)*factor

        end do
      end do 
  end do
  end subroutine get_ellipticity_gll

