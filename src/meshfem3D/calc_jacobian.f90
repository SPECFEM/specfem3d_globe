!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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
!
!> Hejun
! This subroutine recomputes the 3D Jacobian for one element
! based upon 125 GLL points
! Hejun Zhu OCT16,2009

! input: myrank,
!        xstore,ystore,zstore ----- input GLL point coordinate
!        xigll,yigll,zigll ----- GLL points position
!        ispec,nspec       ----- element number

! output: xixstore,xiystore,xizstore,
!         etaxstore,etaystore,etazstore,
!         gammaxstore,gammaystore,gammazstore ------ parameters used to calculate Jacobian


  subroutine recalc_jacobian_gll3D(myrank,xstore,ystore,zstore,xigll,yigll,zigll,&
                                   ispec,nspec,&
                                   xixstore,xiystore,xizstore, &
                                   etaxstore,etaystore,etazstore, &
                                   gammaxstore,gammaystore,gammazstore)

  use constants,only: NGLLX,NGLLY,NGLLZ,CUSTOM_REAL,SIZE_REAL, &
    ZERO,ONE,TINYVAL,VERYSMALLVAL,R_EARTH_KM,RADIANS_TO_DEGREES

  implicit none

  ! input parameter
  integer::myrank,ispec,nspec

  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  double precision, dimension(NGLLX):: xigll
  double precision, dimension(NGLLY):: yigll
  double precision, dimension(NGLLZ):: zigll

  ! output results
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
    gammaxstore,gammaystore,gammazstore

  ! local parameters for this subroutine
  double precision,dimension(NGLLX):: hxir,hpxir
  double precision,dimension(NGLLY):: hetar,hpetar
  double precision,dimension(NGLLZ):: hgammar,hpgammar

  double precision:: xxi,xeta,xgamma,yxi,yeta,ygamma,zxi,zeta,zgamma
  double precision:: xi,eta,gamma

  double precision:: hlagrange,hlagrange_xi,hlagrange_eta,hlagrange_gamma
  double precision:: jacobian,jacobian_inv
  double precision:: xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  double precision:: r,theta,phi
  double precision:: x,y,z

  integer:: i,j,k,i1,j1,k1

  ! test parameters which can be deleted
  double precision:: xmesh,ymesh,zmesh
  double precision:: sumshape,sumdershapexi,sumdershapeeta,sumdershapegamma

  ! first go over all 125 GLL points
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

        xi = xigll(i)
        eta = yigll(j)
        gamma = zigll(k)

        ! calculate Lagrange polynomial and its derivative
        call lagrange_any(xi,NGLLX,xigll,hxir,hpxir)
        call lagrange_any(eta,NGLLY,yigll,hetar,hpetar)
        call lagrange_any(gamma,NGLLZ,zigll,hgammar,hpgammar)

        ! test parameters
        sumshape = ZERO
        sumdershapexi = ZERO
        sumdershapeeta = ZERO
        sumdershapegamma = ZERO

        xmesh = ZERO
        ymesh = ZERO
        zmesh = ZERO

        do k1 = 1,NGLLZ
          do j1 = 1,NGLLY
            do i1 = 1,NGLLX

              hlagrange = hxir(i1) * hetar(j1) * hgammar(k1)
              hlagrange_xi = hpxir(i1) * hetar(j1) * hgammar(k1)
              hlagrange_eta = hxir(i1) * hpetar(j1) * hgammar(k1)
              hlagrange_gamma = hxir(i1) * hetar(j1) * hpgammar(k1)

              x = xstore(i1,j1,k1,ispec)
              y = ystore(i1,j1,k1,ispec)
              z = zstore(i1,j1,k1,ispec)

              xxi = xxi + x * hlagrange_xi
              xeta = xeta + x * hlagrange_eta
              xgamma = xgamma + x * hlagrange_gamma

              yxi = yxi + y * hlagrange_xi
              yeta = yeta + y * hlagrange_eta
              ygamma = ygamma + y * hlagrange_gamma

              zxi = zxi + z * hlagrange_xi
              zeta = zeta + z * hlagrange_eta
              zgamma = zgamma + z * hlagrange_gamma

              ! test the Lagrange polynomial and its derivative
              xmesh = xmesh + x * hlagrange
              ymesh = ymesh + y * hlagrange
              zmesh = zmesh + z * hlagrange

              sumshape = sumshape + hlagrange
              sumdershapexi = sumdershapexi + hlagrange_xi
              sumdershapeeta = sumdershapeeta + hlagrange_eta
              sumdershapegamma = sumdershapegamma + hlagrange_gamma

            enddo
          enddo
        enddo

        ! Check the Lagrange polynomial and its derivative
        if (abs(xmesh - xstore(i,j,k,ispec)) > TINYVAL &
          .or. abs(ymesh - ystore(i,j,k,ispec)) > TINYVAL &
          .or. abs(zmesh - zstore(i,j,k,ispec)) > TINYVAL ) then
          call exit_MPI(myrank,'new mesh is wrong in recalc_jacobian_gll3D.f90')
        endif
        if(abs(sumshape-one) >  TINYVAL) then
          call exit_MPI(myrank,'error shape functions in recalc_jacobian_gll3D.f90')
        endif
        if(abs(sumdershapexi) >  TINYVAL) then
          call exit_MPI(myrank,'error derivative xi in recalc_jacobian_gll3D.f90')
        endif
        if(abs(sumdershapeeta) >  TINYVAL) then
          call exit_MPI(myrank,'error derivative eta in recalc_jacobian_gll3D.f90')
        endif
        if(abs(sumdershapegamma) >  TINYVAL) then
          call exit_MPI(myrank,'error derivative gamma in recalc_jacobian_gll3D.f90')
        endif

        ! Jacobian calculation
        jacobian = xxi*(yeta*zgamma-ygamma*zeta) - &
                   xeta*(yxi*zgamma-ygamma*zxi) + &
                   xgamma*(yxi*zeta-yeta*zxi)

        ! Check the Jacobian
        ! note: when honoring the moho, we squeeze and stretch elements
        !          thus, it can happen that with a coarse mesh resolution, the Jacobian encounters problems
        if(jacobian <= VERYSMALLVAL) then
          ! note: the mesh can have ellipticity, thus the geocentric colatitude might differ from the geographic one
          !
          ! converts position to geocentric coordinates
          call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r,theta,phi)
          print*,'error Jacobian rank:',myrank
          print*,'  location r/lat/lon: ',r*R_EARTH_KM, &
            90.0-(theta*RADIANS_TO_DEGREES),phi*RADIANS_TO_DEGREES
          print*,'  Jacobian: ',jacobian
          call exit_MPI(myrank,'3D Jacobian undefined in recalc_jacobian_gll3D.f90')
        endif

        !     invert the relation (Fletcher p. 50 vol. 2)
        jacobian_inv = ONE / jacobian

        xix = (yeta*zgamma-ygamma*zeta) * jacobian_inv
        xiy = (xgamma*zeta-xeta*zgamma) * jacobian_inv
        xiz = (xeta*ygamma-xgamma*yeta) * jacobian_inv
        etax = (ygamma*zxi-yxi*zgamma) * jacobian_inv
        etay = (xxi*zgamma-xgamma*zxi) * jacobian_inv
        etaz = (xgamma*yxi-xxi*ygamma) * jacobian_inv
        gammax = (yxi*zeta-yeta*zxi) * jacobian_inv
        gammay = (xeta*zxi-xxi*zeta) * jacobian_inv
        gammaz = (xxi*yeta-xeta*yxi) * jacobian_inv

        ! resave the derivatives and the Jacobian
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
        endif

      enddo
    enddo
  enddo

  end subroutine recalc_jacobian_gll3D


!
!-------------------------------------------------------------------------------------------------
!

  ! Hejun Zhu used to recalculate 2D Jacobian according to GLL points rather
  ! than control nodes
  ! Hejun Zhu JAN08, 2010

  ! input parameters:   myrank,ispecb,
  !                     xelm2D,yelm2D,zelm2D,
  !                     xigll,yigll,NSPEC2DMAX_AB,NGLLA,NGLLB

  ! output results:     jacobian2D,normal
  subroutine recalc_jacobian_gll2D(myrank,ispecb, &
                                xelm2D,yelm2D,zelm2D,xigll,yigll,&
                                jacobian2D,normal,NGLLA,NGLLB,NSPEC2DMAX_AB)

  use constants

  implicit none

  ! input parameters
  integer::myrank,ispecb,NSPEC2DMAX_AB,NGLLA,NGLLB

  double precision,dimension(NGLLA,NGLLB)::xelm2D,yelm2D,zelm2D

  double precision,dimension(NGLLA)::xigll
  double precision,dimension(NGLLB)::yigll

  ! output results
  real(kind=CUSTOM_REAL),dimension(NGLLA,NGLLB,NSPEC2DMAX_AB)::jacobian2D
  real(kind=CUSTOM_REAL),dimension(3,NGLLA,NGLLB,NSPEC2DMAX_AB)::normal

  ! local parameters in this subroutine
  integer::i,j,i1,j1
  double precision::xxi,xeta,yxi,yeta,zxi,zeta,&
    xi,eta,xmesh,ymesh,zmesh,hlagrange,hlagrange_xi,hlagrange_eta,&
    sumshape,sumdershapexi,sumdershapeeta,unx,uny,unz,jacobian,jacobian_inv
  double precision,dimension(NGLLA)::hxir,hpxir
  double precision,dimension(NGLLB)::hetar,hpetar

  do j = 1,NGLLB
     do i = 1,NGLLA
        xxi = ZERO
        xeta = ZERO
        yxi = ZERO
        yeta = ZERO
        zxi = ZERO
        zeta = ZERO

        xi=xigll(i)
        eta = yigll(j)

        call lagrange_any(xi,NGLLA,xigll,hxir,hpxir)
        call lagrange_any(eta,NGLLB,yigll,hetar,hpetar)

        xmesh = ZERO
        ymesh = ZERO
        zmesh = ZERO
        sumshape = ZERO
        sumdershapexi = ZERO
        sumdershapeeta = ZERO

        do j1 = 1,NGLLB
           do i1 = 1,NGLLA
              hlagrange = hxir(i1)*hetar(j1)
              hlagrange_xi = hpxir(i1)*hetar(j1)
              hlagrange_eta = hxir(i1)*hpetar(j1)

              xxi = xxi + xelm2D(i1,j1)*hlagrange_xi
              xeta = xeta + xelm2D(i1,j1)*hlagrange_eta

              yxi = yxi + yelm2D(i1,j1)*hlagrange_xi
              yeta = yeta + yelm2D(i1,j1)*hlagrange_eta

              zxi = zxi + zelm2D(i1,j1)*hlagrange_xi
              zeta = zeta + zelm2D(i1,j1)*hlagrange_eta

              xmesh = xmesh + xelm2D(i1,j1)*hlagrange
              ymesh = ymesh + yelm2D(i1,j1)*hlagrange
              zmesh = zmesh + zelm2D(i1,j1)*hlagrange
              sumshape = sumshape + hlagrange
              sumdershapexi = sumdershapexi + hlagrange_xi
              sumdershapeeta = sumdershapeeta + hlagrange_eta
           enddo
        enddo


        ! Check the Lagrange polynomial
        if ( abs(xmesh - xelm2D(i,j)) > TINYVAL &
            .or. abs(ymesh - yelm2D(i,j)) > TINYVAL &
            .or. abs(zmesh - zelm2D(i,j)) > TINYVAL ) then
           call exit_MPI(myrank,'new boundary mesh is wrong in recalc_jacobian_gll2D')
        endif

        if (abs(sumshape-one) >  TINYVAL) then
           call exit_MPI(myrank,'error shape functions in recalc_jacobian_gll2D')
        endif
        if (abs(sumdershapexi) >  TINYVAL) then
           call exit_MPI(myrank,'error derivative xi in recalc_jacobian_gll2D')
        endif
        if (abs(sumdershapeeta) >  TINYVAL) then
           call exit_MPI(myrank,'error derivative eta in recalc_jacobian_gll2D')
        endif

        ! calculates 2D Jacobian
        unx = yxi*zeta - yeta*zxi
        uny = zxi*xeta - zeta*xxi
        unz = xxi*yeta - xeta*yxi
        jacobian = dsqrt(unx*unx + uny*uny + unz*unz)

        ! checks
        if (abs(jacobian) < TINYVAL ) &
          call exit_MPI(myrank,'2D Jacobian undefined in recalc_jacobian_gll2D')

        ! inverts Jacobian
        jacobian_inv = ONE / jacobian

        if (CUSTOM_REAL == SIZE_REAL) then
           jacobian2D(i,j,ispecb) = sngl(jacobian)
           normal(1,i,j,ispecb) = sngl(unx * jacobian_inv)
           normal(2,i,j,ispecb) = sngl(uny * jacobian_inv)
           normal(3,i,j,ispecb) = sngl(unz * jacobian_inv)
        else
           jacobian2D(i,j,ispecb) = jacobian
           normal(1,i,j,ispecb) = unx * jacobian_inv
           normal(2,i,j,ispecb) = uny * jacobian_inv
           normal(3,i,j,ispecb) = unz * jacobian_inv
        endif
     enddo
  enddo

  end subroutine recalc_jacobian_gll2D
