!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

!! DK DK added this for merged version
  subroutine recompute_missing_arrays(myrank, &
     xixstore,xiystore,xizstore, &
     etaxstore,etaystore,etazstore, &
     gammaxstore,gammaystore,gammazstore, &
     xstore,ystore,zstore, &
     xelm_store,yelm_store,zelm_store,ibool,nspec,nglob, &
     xigll,yigll,zigll)

  implicit none

  include "constants.h"

  integer nspec,nglob,myrank

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  real(kind=CUSTOM_REAL), dimension(NGNOD,nspec) :: xelm_store,yelm_store,zelm_store

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xixstore,xiystore,xizstore, &
        etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore,ystore,zstore

  integer i,j,k,ia,ispec

  double precision xxi,xeta,xgamma,yxi,yeta,ygamma,zxi,zeta,zgamma
  double precision xmesh,ymesh,zmesh
  double precision xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  double precision jacobian

! Gauss-Lobatto-Legendre points of integration
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll

! allocate these automatic arrays in the memory stack to avoid memory fragmentation with "allocate()"
! 3D shape functions and their derivatives
  double precision, dimension(NGNOD,NGLLX,NGLLY,NGLLZ) :: shape3D
  double precision, dimension(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ) :: dershape3D

! get the 3-D shape functions
  call get_shape3D(myrank,shape3D,dershape3D,xigll,yigll,zigll)

  do ispec = 1,nspec

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX

      xxi = ZERO
      xeta = ZERO
      xgamma = ZERO
      yxi = ZERO
      yeta = ZERO
      ygamma = ZERO
      zxi = ZERO
      zeta = ZERO
      zgamma = ZERO
      xmesh = ZERO
      ymesh = ZERO
      zmesh = ZERO

      do ia=1,NGNOD
        if(CUSTOM_REAL == SIZE_REAL) then
          xxi = xxi + dershape3D(1,ia,i,j,k)*dble(xelm_store(ia,ispec))
          xeta = xeta + dershape3D(2,ia,i,j,k)*dble(xelm_store(ia,ispec))
          xgamma = xgamma + dershape3D(3,ia,i,j,k)*dble(xelm_store(ia,ispec))
          yxi = yxi + dershape3D(1,ia,i,j,k)*dble(yelm_store(ia,ispec))
          yeta = yeta + dershape3D(2,ia,i,j,k)*dble(yelm_store(ia,ispec))
          ygamma = ygamma + dershape3D(3,ia,i,j,k)*dble(yelm_store(ia,ispec))
          zxi = zxi + dershape3D(1,ia,i,j,k)*dble(zelm_store(ia,ispec))
          zeta = zeta + dershape3D(2,ia,i,j,k)*dble(zelm_store(ia,ispec))
          zgamma = zgamma + dershape3D(3,ia,i,j,k)*dble(zelm_store(ia,ispec))
          xmesh = xmesh + shape3D(ia,i,j,k)*dble(xelm_store(ia,ispec))
          ymesh = ymesh + shape3D(ia,i,j,k)*dble(yelm_store(ia,ispec))
          zmesh = zmesh + shape3D(ia,i,j,k)*dble(zelm_store(ia,ispec))
        else
          xxi = xxi + dershape3D(1,ia,i,j,k)*xelm_store(ia,ispec)
          xeta = xeta + dershape3D(2,ia,i,j,k)*xelm_store(ia,ispec)
          xgamma = xgamma + dershape3D(3,ia,i,j,k)*xelm_store(ia,ispec)
          yxi = yxi + dershape3D(1,ia,i,j,k)*yelm_store(ia,ispec)
          yeta = yeta + dershape3D(2,ia,i,j,k)*yelm_store(ia,ispec)
          ygamma = ygamma + dershape3D(3,ia,i,j,k)*yelm_store(ia,ispec)
          zxi = zxi + dershape3D(1,ia,i,j,k)*zelm_store(ia,ispec)
          zeta = zeta + dershape3D(2,ia,i,j,k)*zelm_store(ia,ispec)
          zgamma = zgamma + dershape3D(3,ia,i,j,k)*zelm_store(ia,ispec)
          xmesh = xmesh + shape3D(ia,i,j,k)*xelm_store(ia,ispec)
          ymesh = ymesh + shape3D(ia,i,j,k)*yelm_store(ia,ispec)
          zmesh = zmesh + shape3D(ia,i,j,k)*zelm_store(ia,ispec)
        endif
      enddo

      jacobian = xxi*(yeta*zgamma-ygamma*zeta) - &
             xeta*(yxi*zgamma-ygamma*zxi) + &
             xgamma*(yxi*zeta-yeta*zxi)

      if(jacobian <= ZERO) call exit_MPI(myrank,'3D Jacobian undefined when recomputing missing arrays')

! invert the relation (Fletcher p. 50 vol. 2)
      xix = (yeta*zgamma-ygamma*zeta) / jacobian
      xiy = (xgamma*zeta-xeta*zgamma) / jacobian
      xiz = (xeta*ygamma-xgamma*yeta) / jacobian
      etax = (ygamma*zxi-yxi*zgamma) / jacobian
      etay = (xxi*zgamma-xgamma*zxi) / jacobian
      etaz = (xgamma*yxi-xxi*ygamma) / jacobian
      gammax = (yxi*zeta-yeta*zxi) / jacobian
      gammay = (xeta*zxi-xxi*zeta) / jacobian
      gammaz = (xxi*yeta-xeta*yxi) / jacobian

! save the derivatives and the jacobian
! store mesh coordinates
! distinguish between single and double precision for reals
      if(CUSTOM_REAL == SIZE_REAL) then
        xixstore(i,j,k,ispec) = sngl(xix)
        xiystore(i,j,k,ispec) = sngl(xiy)
        xizstore(i,j,k,ispec) = sngl(xiz)
        etaxstore(i,j,k,ispec) = sngl(etax)
        etaystore(i,j,k,ispec) = sngl(etay)
        etazstore(i,j,k,ispec) = sngl(etaz)
        gammaxstore(i,j,k,ispec) = sngl(gammax)
        gammaystore(i,j,k,ispec) = sngl(gammay)
        gammazstore(i,j,k,ispec) = sngl(gammaz)

        xstore(ibool(i,j,k,ispec)) = sngl(xmesh)
        ystore(ibool(i,j,k,ispec)) = sngl(ymesh)
        zstore(ibool(i,j,k,ispec)) = sngl(zmesh)
      else
        xixstore(i,j,k,ispec) = xix
        xiystore(i,j,k,ispec) = xiy
        xizstore(i,j,k,ispec) = xiz
        etaxstore(i,j,k,ispec) = etax
        etaystore(i,j,k,ispec) = etay
        etazstore(i,j,k,ispec) = etaz
        gammaxstore(i,j,k,ispec) = gammax
        gammaystore(i,j,k,ispec) = gammay
        gammazstore(i,j,k,ispec) = gammaz

        xstore(ibool(i,j,k,ispec)) = xmesh
        ystore(ibool(i,j,k,ispec)) = ymesh
        zstore(ibool(i,j,k,ispec)) = zmesh
      endif

      enddo
    enddo
  enddo

  enddo

  end subroutine recompute_missing_arrays

