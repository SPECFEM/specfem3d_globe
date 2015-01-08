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

  subroutine get_ellipticity(xelm,yelm,zelm,nspl,rspl,espl,espl2)

  use constants

  implicit none

  integer :: nspl
  double precision :: xelm(NGNOD)
  double precision :: yelm(NGNOD)
  double precision :: zelm(NGNOD)
  double precision :: rspl(NR),espl(NR),espl2(NR)

  ! local parameters
  integer :: ia

  double precision :: ell
  double precision :: r,theta,phi,factor
  double precision :: x,y,z
  double precision :: cost,p20

  do ia = 1,NGNOD

    x = xelm(ia)
    y = yelm(ia)
    z = zelm(ia)

    call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)

    cost = dcos(theta)
! this is the Legendre polynomial of degree two, P2(cos(theta)), see the discussion above eq (14.4) in Dahlen and Tromp (1998)
    p20 = 0.5d0*(3.0d0*cost*cost-1.0d0)

    ! get ellipticity using spline evaluation
    call spline_evaluation(rspl,espl,espl2,nspl,r,ell)

! this is eq (14.4) in Dahlen and Tromp (1998)
    factor = ONE-(TWO/3.0d0)*ell*p20

    xelm(ia) = x*factor
    yelm(ia) = y*factor
    zelm(ia) = z*factor

  enddo

  end subroutine get_ellipticity

!
!-------------------------------------------------------------------------------------------------
!

  !> Hejun
  ! ellipticity at the GLL points
  ! JAN08, 2010
  subroutine get_ellipticity_gll(xstore,ystore,zstore,ispec,nspec,nspl,rspl,espl,espl2)

  use constants

  implicit none

  integer :: nspl
  integer :: ispec,nspec
  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore
  double precision :: rspl(NR),espl(NR),espl2(NR)

  ! local parameters
  integer :: i,j,k

  double precision :: ell
  double precision :: r,theta,phi,factor
  double precision :: x,y,z
  double precision :: cost,p20

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        x = xstore(i,j,k,ispec)
        y = ystore(i,j,k,ispec)
        z = zstore(i,j,k,ispec)

        call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)

        cost = dcos(theta)
! this is the Legendre polynomial of degree two, P2(cos(theta)), see the discussion above eq (14.4) in Dahlen and Tromp (1998)
        p20 = 0.5d0*(3.0d0*cost*cost-1.0d0)

        ! get ellipticity using spline evaluation
        call spline_evaluation(rspl,espl,espl2,nspl,r,ell)

! this is eq (14.4) in Dahlen and Tromp (1998)
        factor = ONE-(TWO/3.0d0)*ell*p20

        xstore(i,j,k,ispec) = x*factor
        ystore(i,j,k,ispec) = y*factor
        zstore(i,j,k,ispec) = z*factor

      enddo
    enddo
  enddo

  end subroutine get_ellipticity_gll


!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_ellipticity_single_point(x,y,z,nspl,rspl,espl,espl2)

  use constants

  implicit none

  integer :: nspl
  double precision :: x,y,z
  double precision :: rspl(NR),espl(NR),espl2(NR)

  ! local parameters
  double precision :: ell
  double precision :: r,theta,phi,factor
  double precision :: cost,p20

  call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)

  cost = dcos(theta)
! this is the Legendre polynomial of degree two, P2(cos(theta)), see the discussion above eq (14.4) in Dahlen and Tromp (1998)
  p20 = 0.5d0*(3.0d0*cost*cost-1.0d0)

  ! get ellipticity using spline evaluation
  call spline_evaluation(rspl,espl,espl2,nspl,r,ell)

! this is eq (14.4) in Dahlen and Tromp (1998)
  factor = ONE-(TWO/3.0d0)*ell*p20

  x = x*factor
  y = y*factor
  z = z*factor

  end subroutine get_ellipticity_single_point

