!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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
! This subroutine recomputes the 3D jacobian for one element
! based upon 125 GLL points
! Hejun Zhu OCT16,2009

! input: myrank,
!        xstore,ystore,zstore ----- input GLL point coordinate
!        xigll,yigll,zigll ----- gll points position
!        ispec,nspec       ----- element number
!        ACTUALLY_STORE_ARRAYS   ------ save array or not

! output: xixstore,xiystore,xizstore,
!         etaxstore,etaystore,etazstore,
!         gammaxstore,gammaystore,gammazstore ------ parameters used to calculate jacobian


  subroutine recalc_jacobian_gll3D(myrank,xstore,ystore,zstore,xigll,yigll,zigll,&
                                ispec,nspec,ACTUALLY_STORE_ARRAYS,&
                                xixstore,xiystore,xizstore, &
                                etaxstore,etaystore,etazstore, &
                                gammaxstore,gammaystore,gammazstore)

  implicit none

  include "constants.h"

  ! input parameter
  integer::myrank,ispec,nspec
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore
  double precision, dimension(NGLLX):: xigll
  double precision, dimension(NGLLY):: yigll
  double precision, dimension(NGLLZ):: zigll
  logical::ACTUALLY_STORE_ARRAYS


  ! output results
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
                        xixstore,xiystore,xizstore,&
                        etaxstore,etaystore,etazstore,&
                        gammaxstore,gammaystore,gammazstore


  ! local parameters for this subroutine
  integer:: i,j,k,i1,j1,k1
  double precision:: xxi,xeta,xgamma,yxi,yeta,ygamma,zxi,zeta,zgamma
  double precision:: xi,eta,gamma
  double precision,dimension(NGLLX):: hxir,hpxir
  double precision,dimension(NGLLY):: hetar,hpetar
  double precision,dimension(NGLLZ):: hgammar,hpgammar
  double precision:: hlagrange,hlagrange_xi,hlagrange_eta,hlagrange_gamma
  double precision:: jacobian
  double precision:: xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  double precision:: r,theta,phi


  ! test parameters which can be deleted
  double precision:: xmesh,ymesh,zmesh
  double precision:: sumshape,sumdershapexi,sumdershapeeta,sumdershapegamma

  ! first go over all 125 GLL points
  do k=1,NGLLZ
     do j=1,NGLLY
        do i=1,NGLLX

            xxi = 0.0
            xeta = 0.0
            xgamma = 0.0
            yxi = 0.0
            yeta = 0.0
            ygamma = 0.0
            zxi = 0.0
            zeta = 0.0
            zgamma = 0.0

            xi = xigll(i)
            eta = yigll(j)
            gamma = zigll(k)

            ! calculate lagrange polynomial and its derivative
            call lagrange_any(xi,NGLLX,xigll,hxir,hpxir)
            call lagrange_any(eta,NGLLY,yigll,hetar,hpetar)
            call lagrange_any(gamma,NGLLZ,zigll,hgammar,hpgammar)

            ! test parameters
            sumshape = 0.0
            sumdershapexi = 0.0
            sumdershapeeta = 0.0
            sumdershapegamma = 0.0
            xmesh = 0.0
            ymesh = 0.0
            zmesh = 0.0


            do k1 = 1,NGLLZ
               do j1 = 1,NGLLY
                  do i1 = 1,NGLLX
                     hlagrange = hxir(i1)*hetar(j1)*hgammar(k1)
                     hlagrange_xi = hpxir(i1)*hetar(j1)*hgammar(k1)
                     hlagrange_eta = hxir(i1)*hpetar(j1)*hgammar(k1)
                     hlagrange_gamma = hxir(i1)*hetar(j1)*hpgammar(k1)


                     xxi = xxi + xstore(i1,j1,k1,ispec)*hlagrange_xi
                     xeta = xeta + xstore(i1,j1,k1,ispec)*hlagrange_eta
                     xgamma = xgamma + xstore(i1,j1,k1,ispec)*hlagrange_gamma

                     yxi = yxi + ystore(i1,j1,k1,ispec)*hlagrange_xi
                     yeta = yeta + ystore(i1,j1,k1,ispec)*hlagrange_eta
                     ygamma = ygamma + ystore(i1,j1,k1,ispec)*hlagrange_gamma

                     zxi = zxi + zstore(i1,j1,k1,ispec)*hlagrange_xi
                     zeta = zeta + zstore(i1,j1,k1,ispec)*hlagrange_eta
                     zgamma = zgamma + zstore(i1,j1,k1,ispec)*hlagrange_gamma

                     ! test the lagrange polynomial and its derivate
                     xmesh = xmesh + xstore(i1,j1,k1,ispec)*hlagrange
                     ymesh = ymesh + ystore(i1,j1,k1,ispec)*hlagrange
                     zmesh = zmesh + zstore(i1,j1,k1,ispec)*hlagrange
                     sumshape = sumshape + hlagrange
                     sumdershapexi = sumdershapexi + hlagrange_xi
                     sumdershapeeta = sumdershapeeta + hlagrange_eta
                     sumdershapegamma = sumdershapegamma + hlagrange_gamma

                  end do
               end do
            end do

            ! Check the lagrange polynomial and its derivative
            if (abs(xmesh - xstore(i,j,k,ispec)) > TINYVAL &
              .or. abs(ymesh - ystore(i,j,k,ispec)) > TINYVAL &
              .or. abs(zmesh - zstore(i,j,k,ispec)) > TINYVAL ) then
                    call exit_MPI(myrank,'new mesh are wrong in recalc_jacobian_gall3D.f90')
            end if
            if(abs(sumshape-one) >  TINYVAL) then
                    call exit_MPI(myrank,'error shape functions in recalc_jacobian_gll3D.f90')
            end if
            if(abs(sumdershapexi) >  TINYVAL) then
                    call exit_MPI(myrank,'error derivative xi in recalc_jacobian_gll3D.f90')
            end if
            if(abs(sumdershapeeta) >  TINYVAL) then
                    call exit_MPI(myrank,'error derivative eta in recalc_jacobian_gll3D.f90')
            end if
            if(abs(sumdershapegamma) >  TINYVAL) then
                    call exit_MPI(myrank,'error derivative gamma in recalc_jacobian_gll3D.f90')
            end if


            jacobian = xxi*(yeta*zgamma-ygamma*zeta) - &
                 xeta*(yxi*zgamma-ygamma*zxi) + &
                 xgamma*(yxi*zeta-yeta*zxi)

            ! Check the jacobian
            ! note: when honoring the moho, we squeeze and stretch elements
            !          thus, it can happen that with a coarse mesh resolution, the jacobian encounters problems
            if(jacobian <= ZERO) then
              call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r,theta,phi)
              print*,'error jacobian rank:',myrank
              print*,'  location r/lat/lon: ',r*R_EARTH_KM,90.0-theta*180./PI,phi*180./PI
              print*,'  jacobian: ',jacobian
              call exit_MPI(myrank,'3D Jacobian undefined in recalc_jacobian_gll3D.f90')
            end if

            !     invert the relation (Fletcher p. 50 vol. 2)
            xix = (yeta*zgamma-ygamma*zeta) / jacobian
            xiy = (xgamma*zeta-xeta*zgamma) / jacobian
            xiz = (xeta*ygamma-xgamma*yeta) / jacobian
            etax = (ygamma*zxi-yxi*zgamma) / jacobian
            etay = (xxi*zgamma-xgamma*zxi) / jacobian
            etaz = (xgamma*yxi-xxi*ygamma) / jacobian
            gammax = (yxi*zeta-yeta*zxi) / jacobian
            gammay = (xeta*zxi-xxi*zeta) / jacobian
            gammaz = (xxi*yeta-xeta*yxi) / jacobian


            ! resave the derivatives and the jacobian
            ! distinguish between single and double precision for reals
            if (ACTUALLY_STORE_ARRAYS) then
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
             end if
        enddo
    enddo
  enddo

  end subroutine recalc_jacobian_gll3D


!
!-------------------------------------------------------------------------------------------------
!

  ! Hejun Zhu used to recalculate 2D jacobian according to gll points rather
  ! than control nodes
  ! Hejun Zhu JAN08, 2010

  ! input parameters:   myrank,ispecb,
  !                     xelm2D,yelm2D,zelm2D,
  !                     xigll,yigll,NSPEC2DMAX_AB,NGLLA,NGLLB

  ! output results:     jacobian2D,normal
  subroutine recalc_jacobian_gll2D(myrank,ispecb, &
                                xelm2D,yelm2D,zelm2D,xigll,yigll,&
                                jacobian2D,normal,NGLLA,NGLLB,NSPEC2DMAX_AB)

  implicit none
  include "constants.h"
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
                sumshape,sumdershapexi,sumdershapeeta,unx,uny,unz,jacobian
  double precision,dimension(NGLLA)::hxir,hpxir
  double precision,dimension(NGLLB)::hetar,hpetar

  do j = 1,NGLLB
     do i = 1,NGLLA
        xxi = 0.0
        xeta = 0.0
        yxi = 0.0
        yeta = 0.0
        zxi = 0.0
        zeta = 0.0

        xi=xigll(i)
        eta = yigll(j)

        call lagrange_any(xi,NGLLA,xigll,hxir,hpxir)
        call lagrange_any(eta,NGLLB,yigll,hetar,hpetar)


        xmesh = 0.0
        ymesh = 0.0
        zmesh = 0.0
        sumshape = 0.0
        sumdershapexi = 0.0
        sumdershapeeta = 0.0
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
           end do
        end do


        ! Check the lagrange polynomial
        if ( abs(xmesh - xelm2D(i,j)) > TINYVAL &
            .or. abs(ymesh - yelm2D(i,j)) > TINYVAL &
            .or. abs(zmesh - zelm2D(i,j)) > TINYVAL ) then
           call exit_MPI(myrank,'new boundary mesh is wrong in recalc_jacobian_gll2D')
        end if

        if (abs(sumshape-one) >  TINYVAL) then
           call exit_MPI(myrank,'error shape functions in recalc_jacobian_gll2D')
        end if
        if (abs(sumdershapexi) >  TINYVAL) then
           call exit_MPI(myrank,'error derivative xi in recalc_jacobian_gll2D')
        end if
        if (abs(sumdershapeeta) >  TINYVAL) then
           call exit_MPI(myrank,'error derivative eta in recalc_jacobian_gll2D')
        end if

        unx = yxi*zeta - yeta*zxi
        uny = zxi*xeta - zeta*xxi
        unz = xxi*yeta - xeta*yxi
        jacobian = dsqrt(unx**2+uny**2+unz**2)
        if (abs(jacobian) < TINYVAL ) call exit_MPI(myrank,'2D Jacobian undefined in recalc_jacobian_gll2D')

        if (CUSTOM_REAL == SIZE_REAL) then
           jacobian2D(i,j,ispecb)=sngl(jacobian)
           normal(1,i,j,ispecb)=sngl(unx/jacobian)
           normal(2,i,j,ispecb)=sngl(uny/jacobian)
           normal(3,i,j,ispecb)=sngl(unz/jacobian)
        else
           jacobian2D(i,j,ispecb)=jacobian
           normal(1,i,j,ispecb)=unx/jacobian
           normal(2,i,j,ispecb)=uny/jacobian
           normal(3,i,j,ispecb)=unz/jacobian
        endif
     end do
  end do

  end subroutine recalc_jacobian_gll2D

!
!-------------------------------------------------------------------------------------------------
!
! deprecated...
!
!  subroutine calc_jacobian(myrank,xixstore,xiystore,xizstore, &
!     etaxstore,etaystore,etazstore, &
!     gammaxstore,gammaystore,gammazstore, &
!     xstore,ystore,zstore, &
!     xelm,yelm,zelm,shape3D,dershape3D,ispec,nspec,ACTUALLY_STORE_ARRAYS)
!
!  implicit none
!
!  include "constants.h"
!
!  integer ispec,nspec,myrank
!
!  logical ACTUALLY_STORE_ARRAYS
!
!  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)
!  double precision dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)
!
!  double precision xelm(NGNOD)
!  double precision yelm(NGNOD)
!  double precision zelm(NGNOD)
!
!  real(kind=CUSTOM_REAL) xixstore(NGLLX,NGLLY,NGLLZ,nspec), &
!                         xiystore(NGLLX,NGLLY,NGLLZ,nspec), &
!                         xizstore(NGLLX,NGLLY,NGLLZ,nspec), &
!                         etaxstore(NGLLX,NGLLY,NGLLZ,nspec), &
!                         etaystore(NGLLX,NGLLY,NGLLZ,nspec), &
!                         etazstore(NGLLX,NGLLY,NGLLZ,nspec), &
!                         gammaxstore(NGLLX,NGLLY,NGLLZ,nspec), &
!                         gammaystore(NGLLX,NGLLY,NGLLZ,nspec), &
!                         gammazstore(NGLLX,NGLLY,NGLLZ,nspec)
!
!  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
!  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
!  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)
!
!  integer i,j,k,ia
!
!  double precision xxi,xeta,xgamma,yxi,yeta,ygamma,zxi,zeta,zgamma
!  double precision xmesh,ymesh,zmesh
!  double precision xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
!  double precision jacobian
!
!  do k=1,NGLLZ
!    do j=1,NGLLY
!      do i=1,NGLLX
!
!      xxi = ZERO
!      xeta = ZERO
!      xgamma = ZERO
!      yxi = ZERO
!      yeta = ZERO
!      ygamma = ZERO
!      zxi = ZERO
!      zeta = ZERO
!      zgamma = ZERO
!      xmesh = ZERO
!      ymesh = ZERO
!      zmesh = ZERO
!
!      do ia=1,NGNOD
!        xxi = xxi + dershape3D(1,ia,i,j,k)*xelm(ia)
!        xeta = xeta + dershape3D(2,ia,i,j,k)*xelm(ia)
!        xgamma = xgamma + dershape3D(3,ia,i,j,k)*xelm(ia)
!        yxi = yxi + dershape3D(1,ia,i,j,k)*yelm(ia)
!        yeta = yeta + dershape3D(2,ia,i,j,k)*yelm(ia)
!        ygamma = ygamma + dershape3D(3,ia,i,j,k)*yelm(ia)
!        zxi = zxi + dershape3D(1,ia,i,j,k)*zelm(ia)
!        zeta = zeta + dershape3D(2,ia,i,j,k)*zelm(ia)
!        zgamma = zgamma + dershape3D(3,ia,i,j,k)*zelm(ia)
!        xmesh = xmesh + shape3D(ia,i,j,k)*xelm(ia)
!        ymesh = ymesh + shape3D(ia,i,j,k)*yelm(ia)
!        zmesh = zmesh + shape3D(ia,i,j,k)*zelm(ia)
!      enddo
!
!      jacobian = xxi*(yeta*zgamma-ygamma*zeta) - &
!             xeta*(yxi*zgamma-ygamma*zxi) + &
!             xgamma*(yxi*zeta-yeta*zxi)
!
!      if(jacobian <= ZERO) then
!        print*,'jacobian error:',myrank
!        print*,'  point ijk:',i,j,k,ispec
!        print*,'  xyz:',xmesh,ymesh,zmesh
!        call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,xxi,xeta,xgamma)
!        print*,'  r/lat/lon:',xxi*R_EARTH_KM,90.0-xeta*180./PI,xgamma*180./PI
!        print*,'  nodes:'
!        do ia=1,NGNOD
!          print*,xelm(ia),yelm(ia),zelm(ia)
!        enddo
!        print*
!        print*,'maybe check with CAP smoothing'
!        call exit_MPI(myrank,'3D Jacobian undefined')
!      endif
!
!! invert the relation (Fletcher p. 50 vol. 2)
!      xix = (yeta*zgamma-ygamma*zeta) / jacobian
!      xiy = (xgamma*zeta-xeta*zgamma) / jacobian
!      xiz = (xeta*ygamma-xgamma*yeta) / jacobian
!      etax = (ygamma*zxi-yxi*zgamma) / jacobian
!      etay = (xxi*zgamma-xgamma*zxi) / jacobian
!      etaz = (xgamma*yxi-xxi*ygamma) / jacobian
!      gammax = (yxi*zeta-yeta*zxi) / jacobian
!      gammay = (xeta*zxi-xxi*zeta) / jacobian
!      gammaz = (xxi*yeta-xeta*yxi) / jacobian
!
!! save the derivatives and the jacobian
!! distinguish between single and double precision for reals
!      if(ACTUALLY_STORE_ARRAYS) then
!        if(CUSTOM_REAL == SIZE_REAL) then
!          xixstore(i,j,k,ispec) = sngl(xix)
!          xiystore(i,j,k,ispec) = sngl(xiy)
!          xizstore(i,j,k,ispec) = sngl(xiz)
!          etaxstore(i,j,k,ispec) = sngl(etax)
!          etaystore(i,j,k,ispec) = sngl(etay)
!          etazstore(i,j,k,ispec) = sngl(etaz)
!          gammaxstore(i,j,k,ispec) = sngl(gammax)
!          gammaystore(i,j,k,ispec) = sngl(gammay)
!          gammazstore(i,j,k,ispec) = sngl(gammaz)
!        else
!          xixstore(i,j,k,ispec) = xix
!          xiystore(i,j,k,ispec) = xiy
!          xizstore(i,j,k,ispec) = xiz
!          etaxstore(i,j,k,ispec) = etax
!          etaystore(i,j,k,ispec) = etay
!          etazstore(i,j,k,ispec) = etaz
!          gammaxstore(i,j,k,ispec) = gammax
!          gammaystore(i,j,k,ispec) = gammay
!          gammazstore(i,j,k,ispec) = gammaz
!        endif
!      endif
!
!! store mesh coordinates
!!      xstore(i,j,k,ispec) = xmesh
!!      ystore(i,j,k,ispec) = ymesh
!!      zstore(i,j,k,ispec) = zmesh
!
!      enddo
!    enddo
!  enddo
!
!  end subroutine calc_jacobian
!


