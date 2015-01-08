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

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"



  subroutine compute_element_strain_undoatt_Dev(ispec,nglob,nspec, &
                                                displ,ibool, &
                                                hprime_xx,hprime_xxT,&
                                                xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                                                epsilondev_loc,eps_trace_over_3_loc)

! computes strain for single element based on Deville routine setup (NGLLX == NGLLY == NGLLZ == 5)

  use constants

  implicit none

  integer,intent(in) :: ispec,nglob,nspec

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: displ
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX),intent(in) :: hprime_xx,hprime_xxT

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5),intent(out) :: epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: eps_trace_over_3_loc

  !  local variable
  integer :: iglob

  real(kind=CUSTOM_REAL) :: templ
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) :: duxdxl,duydyl,duzdzl,duxdyl,duydxl,duzdxl,duxdzl,duzdyl,duydzl, &
                            duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl, &
                            duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  DO_LOOP_IJK

    iglob = ibool(INDEX_IJK,ispec)
    dummyx_loc(INDEX_IJK) = displ(1,iglob)
    dummyy_loc(INDEX_IJK) = displ(2,iglob)
    dummyz_loc(INDEX_IJK) = displ(3,iglob)

  ENDDO_LOOP_IJK

  ! Deville optimizations
  ! computes 1. matrix multiplication for tempx1,..
  call mxm5_3comp_singleA(hprime_xx,m1,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1,m2)
  ! computes 2. matrix multiplication for tempx2,..
  call mxm5_3comp_3dmat_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m1,hprime_xxT,m1,tempx2,tempy2,tempz2,NGLLX)
  ! computes 3. matrix multiplication for tempx1,..
  call mxm5_3comp_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m2,hprime_xxT,tempx3,tempy3,tempz3,m1)

  DO_LOOP_IJK

    ! get derivatives of ux, uy and uz with respect to x, y and z
    xixl = xix(INDEX_IJK,ispec)
    xiyl = xiy(INDEX_IJK,ispec)
    xizl = xiz(INDEX_IJK,ispec)
    etaxl = etax(INDEX_IJK,ispec)
    etayl = etay(INDEX_IJK,ispec)
    etazl = etaz(INDEX_IJK,ispec)
    gammaxl = gammax(INDEX_IJK,ispec)
    gammayl = gammay(INDEX_IJK,ispec)
    gammazl = gammaz(INDEX_IJK,ispec)

    duxdxl = xixl*tempx1(INDEX_IJK) + etaxl*tempx2(INDEX_IJK) + gammaxl*tempx3(INDEX_IJK)
    duxdyl = xiyl*tempx1(INDEX_IJK) + etayl*tempx2(INDEX_IJK) + gammayl*tempx3(INDEX_IJK)
    duxdzl = xizl*tempx1(INDEX_IJK) + etazl*tempx2(INDEX_IJK) + gammazl*tempx3(INDEX_IJK)

    duydxl = xixl*tempy1(INDEX_IJK) + etaxl*tempy2(INDEX_IJK) + gammaxl*tempy3(INDEX_IJK)
    duydyl = xiyl*tempy1(INDEX_IJK) + etayl*tempy2(INDEX_IJK) + gammayl*tempy3(INDEX_IJK)
    duydzl = xizl*tempy1(INDEX_IJK) + etazl*tempy2(INDEX_IJK) + gammazl*tempy3(INDEX_IJK)

    duzdxl = xixl*tempz1(INDEX_IJK) + etaxl*tempz2(INDEX_IJK) + gammaxl*tempz3(INDEX_IJK)
    duzdyl = xiyl*tempz1(INDEX_IJK) + etayl*tempz2(INDEX_IJK) + gammayl*tempz3(INDEX_IJK)
    duzdzl = xizl*tempz1(INDEX_IJK) + etazl*tempz2(INDEX_IJK) + gammazl*tempz3(INDEX_IJK)

    ! precompute some sums to save CPU time
    duxdxl_plus_duydyl = duxdxl + duydyl
    duxdxl_plus_duzdzl = duxdxl + duzdzl
    duydyl_plus_duzdzl = duydyl + duzdzl
    duxdyl_plus_duydxl = duxdyl + duydxl
    duzdxl_plus_duxdzl = duzdxl + duxdzl
    duzdyl_plus_duydzl = duzdyl + duydzl

    ! strains
    templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
    eps_trace_over_3_loc(INDEX_IJK) = templ
    epsilondev_loc(INDEX_IJK,1) = duxdxl - templ
    epsilondev_loc(INDEX_IJK,2) = duydyl - templ
    epsilondev_loc(INDEX_IJK,3) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
    epsilondev_loc(INDEX_IJK,4) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
    epsilondev_loc(INDEX_IJK,5) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl

  ENDDO_LOOP_IJK

  contains

!--------------------------------------------------------------------------------------------
!
! matrix-matrix multiplications
!
! subroutines adapted from Deville, Fischer and Mund, High-order methods
! for incompressible fluid flow, Cambridge University Press (2002),
! pages 386 and 389 and Figure 8.3.1
!
!--------------------------------------------------------------------------------------------
!
! note: the matrix-matrix multiplications are used for very small matrices ( 5 x 5 x 5 elements);
!       thus, calling external optimized libraries for these multiplications are in general slower
!
! please leave the routines here to help compilers inlining the code

  subroutine mxm5_3comp_singleA(A,n1,B1,B2,B3,C1,C2,C3,n3)

! 3 different arrays for x/y/z-components, 2-dimensional arrays (25,5)/(5,25), same B matrix for all 3 component arrays

  use constants_solver,only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5) :: A
  real(kind=CUSTOM_REAL),dimension(5,n3) :: B1,B2,B3
  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C1,C2,C3

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
    do i = 1,n1
      C1(i,j) =  A(i,1) * B1(1,j) &
               + A(i,2) * B1(2,j) &
               + A(i,3) * B1(3,j) &
               + A(i,4) * B1(4,j) &
               + A(i,5) * B1(5,j)

      C2(i,j) =  A(i,1) * B2(1,j) &
               + A(i,2) * B2(2,j) &
               + A(i,3) * B2(3,j) &
               + A(i,4) * B2(4,j) &
               + A(i,5) * B2(5,j)

      C3(i,j) =  A(i,1) * B3(1,j) &
               + A(i,2) * B3(2,j) &
               + A(i,3) * B3(3,j) &
               + A(i,4) * B3(4,j) &
               + A(i,5) * B3(5,j)
    enddo
  enddo

  end subroutine mxm5_3comp_singleA


!--------------------------------------------------------------------------------------------

  subroutine mxm5_3comp_singleB(A1,A2,A3,n1,B,C1,C2,C3,n3)

! 3 different arrays for x/y/z-components, 2-dimensional arrays (25,5)/(5,25), same B matrix for all 3 component arrays

  use constants_solver,only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(5,n3) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C1,C2,C3

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
    do i = 1,n1
      C1(i,j) =  A1(i,1) * B(1,j) &
               + A1(i,2) * B(2,j) &
               + A1(i,3) * B(3,j) &
               + A1(i,4) * B(4,j) &
               + A1(i,5) * B(5,j)

      C2(i,j) =  A2(i,1) * B(1,j) &
               + A2(i,2) * B(2,j) &
               + A2(i,3) * B(3,j) &
               + A2(i,4) * B(4,j) &
               + A2(i,5) * B(5,j)

      C3(i,j) =  A3(i,1) * B(1,j) &
               + A3(i,2) * B(2,j) &
               + A3(i,3) * B(3,j) &
               + A3(i,4) * B(4,j) &
               + A3(i,5) * B(5,j)
    enddo
  enddo

  end subroutine mxm5_3comp_singleB


!--------------------------------------------------------------------------------------------

  subroutine mxm5_3comp_3dmat_singleB(A1,A2,A3,n1,B,n2,C1,C2,C3,n3)

! 3 different arrays for x/y/z-components, 3-dimensional arrays (5,5,5), same B matrix for all 3 component arrays

  use constants_solver,only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,5,n3) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(5,n2) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3) :: C1,C2,C3

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do j = 1,n2
    do i = 1,n1
      ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
      do k = 1,n3
        C1(i,j,k) =  A1(i,1,k) * B(1,j) &
                   + A1(i,2,k) * B(2,j) &
                   + A1(i,3,k) * B(3,j) &
                   + A1(i,4,k) * B(4,j) &
                   + A1(i,5,k) * B(5,j)

        C2(i,j,k) =  A2(i,1,k) * B(1,j) &
                   + A2(i,2,k) * B(2,j) &
                   + A2(i,3,k) * B(3,j) &
                   + A2(i,4,k) * B(4,j) &
                   + A2(i,5,k) * B(5,j)

        C3(i,j,k) =  A3(i,1,k) * B(1,j) &
                   + A3(i,2,k) * B(2,j) &
                   + A3(i,3,k) * B(3,j) &
                   + A3(i,4,k) * B(4,j) &
                   + A3(i,5,k) * B(5,j)
      enddo
    enddo
  enddo

  end subroutine mxm5_3comp_3dmat_singleB


  end subroutine compute_element_strain_undoatt_Dev


!
!--------------------------------------------------------------------------------------------
!

  subroutine compute_element_strain_undoatt_noDev(ispec,nglob,nspec, &
                                                  displ, &
                                                  hprime_xx,hprime_yy,hprime_zz, &
                                                  ibool,&
                                                  xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                                                  epsilondev_loc,eps_trace_over_3_loc)

! computes strain for single element

  use constants

  implicit none

  ! element id
  integer,intent(in) :: ispec

  integer,intent(in) :: NSPEC,NGLOB

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB),intent(in) :: displ

  ! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX),intent(in) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY),intent(in) :: hprime_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ),intent(in) :: hprime_zz

  ! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: &
    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5),intent(out) :: epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: eps_trace_over_3_loc

  ! local parameters
  integer :: iglob
  integer :: i,j,k,l

  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) tempz1l,tempz2l,tempz3l

  real(kind=CUSTOM_REAL) hp1,hp2,hp3
  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl
  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        tempx1l = 0._CUSTOM_REAL
        tempx2l = 0._CUSTOM_REAL
        tempx3l = 0._CUSTOM_REAL

        tempy1l = 0._CUSTOM_REAL
        tempy2l = 0._CUSTOM_REAL
        tempy3l = 0._CUSTOM_REAL

        tempz1l = 0._CUSTOM_REAL
        tempz2l = 0._CUSTOM_REAL
        tempz3l = 0._CUSTOM_REAL

        do l = 1,NGLLX
          hp1 = hprime_xx(i,l)
          iglob = ibool(l,j,k,ispec)
          tempx1l = tempx1l + displ(1,iglob)*hp1
          tempy1l = tempy1l + displ(2,iglob)*hp1
          tempz1l = tempz1l + displ(3,iglob)*hp1
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l = 1,NGLLY
          hp2 = hprime_yy(j,l)
          iglob = ibool(i,l,k,ispec)
          tempx2l = tempx2l + displ(1,iglob)*hp2
          tempy2l = tempy2l + displ(2,iglob)*hp2
          tempz2l = tempz2l + displ(3,iglob)*hp2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l = 1,NGLLZ
          hp3 = hprime_zz(k,l)
          iglob = ibool(i,j,l,ispec)
          tempx3l = tempx3l + displ(1,iglob)*hp3
          tempy3l = tempy3l + displ(2,iglob)*hp3
          tempz3l = tempz3l + displ(3,iglob)*hp3
        enddo

        ! get derivatives of ux, uy and uz with respect to x, y and z
        xixl = xix(i,j,k,ispec)
        xiyl = xiy(i,j,k,ispec)
        xizl = xiz(i,j,k,ispec)
        etaxl = etax(i,j,k,ispec)
        etayl = etay(i,j,k,ispec)
        etazl = etaz(i,j,k,ispec)
        gammaxl = gammax(i,j,k,ispec)
        gammayl = gammay(i,j,k,ispec)
        gammazl = gammaz(i,j,k,ispec)

        ! compute the Jacobian
        !jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
        !              - xiyl*(etaxl*gammazl-etazl*gammaxl) &
        !              + xizl*(etaxl*gammayl-etayl*gammaxl))

        duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
        duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
        duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

        duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
        duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
        duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

        duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
        duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
        duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l

        ! precompute some sums to save CPU time
        duxdxl_plus_duydyl = duxdxl + duydyl
        duxdxl_plus_duzdzl = duxdxl + duzdzl
        duydyl_plus_duzdzl = duydyl + duzdzl
        duxdyl_plus_duydxl = duxdyl + duydxl
        duzdxl_plus_duxdzl = duzdxl + duxdzl
        duzdyl_plus_duydzl = duzdyl + duydzl

        ! compute deviatoric strain
        eps_trace_over_3_loc(i,j,k) = ONE_THIRD * (duxdxl + duydyl + duzdzl)
        epsilondev_loc(i,j,k,1) = duxdxl - eps_trace_over_3_loc(i,j,k)
        epsilondev_loc(i,j,k,2) = duydyl - eps_trace_over_3_loc(i,j,k)
        epsilondev_loc(i,j,k,3) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
        epsilondev_loc(i,j,k,4) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
        epsilondev_loc(i,j,k,5) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl

      enddo ! NGLLX
    enddo ! NGLLY
  enddo ! NGLLZ

  end subroutine compute_element_strain_undoatt_noDev


!--------------------------------------------------------------------------------------------
!
! strain separated into single xx,yy,xy,xz,yz-component arrays
!
!--------------------------------------------------------------------------------------------


  subroutine compute_element_strain_att_Dev(ispec,nglob,nspec, &
                                           displ,veloc,deltat, &
                                           ibool, &
                                           hprime_xx,hprime_xxT,&
                                           xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                                           epsilondev_xx_loc_nplus1, &
                                           epsilondev_yy_loc_nplus1, &
                                           epsilondev_xy_loc_nplus1, &
                                           epsilondev_xz_loc_nplus1, &
                                           epsilondev_yz_loc_nplus1, &
                                           nspec_strain_only,eps_trace_over_3_loc_nplus1)

  use constants

  implicit none

  integer :: ispec,nglob,nspec
  real(kind=CUSTOM_REAL) :: deltat

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob) :: displ,veloc

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprime_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: epsilondev_xx_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: epsilondev_yy_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: epsilondev_xy_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: epsilondev_xz_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: epsilondev_yz_loc_nplus1

  integer :: nspec_strain_only
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_strain_only) :: eps_trace_over_3_loc_nplus1

  ! local variable
  integer :: iglob

  real(kind=CUSTOM_REAL) :: templ
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) :: duxdxl,duydyl,duzdzl,duxdyl,duydxl,duzdxl,duxdzl,duzdyl,duydzl, &
                            duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl, &
                            duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  DO_LOOP_IJK

    iglob = ibool(INDEX_IJK,ispec)
    dummyx_loc(INDEX_IJK) = displ(1,iglob) + deltat * veloc(1,iglob)
    dummyy_loc(INDEX_IJK) = displ(2,iglob) + deltat * veloc(2,iglob)
    dummyz_loc(INDEX_IJK) = displ(3,iglob) + deltat * veloc(3,iglob)

  ENDDO_LOOP_IJK

  ! Deville optimizations
  ! computes 1. matrix multiplication for tempx1,..
  call mxm5_3comp_singleA(hprime_xx,m1,dummyx_loc,dummyy_loc,dummyz_loc,tempx1,tempy1,tempz1,m2)
  ! computes 2. matrix multiplication for tempx2,..
  call mxm5_3comp_3dmat_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m1,hprime_xxT,m1,tempx2,tempy2,tempz2,NGLLX)
  ! computes 3. matrix multiplication for tempx1,..
  call mxm5_3comp_singleB(dummyx_loc,dummyy_loc,dummyz_loc,m2,hprime_xxT,tempx3,tempy3,tempz3,m1)

  DO_LOOP_IJK

    ! get derivatives of ux, uy and uz with respect to x, y and z
    xixl = xix(INDEX_IJK,ispec)
    xiyl = xiy(INDEX_IJK,ispec)
    xizl = xiz(INDEX_IJK,ispec)
    etaxl = etax(INDEX_IJK,ispec)
    etayl = etay(INDEX_IJK,ispec)
    etazl = etaz(INDEX_IJK,ispec)
    gammaxl = gammax(INDEX_IJK,ispec)
    gammayl = gammay(INDEX_IJK,ispec)
    gammazl = gammaz(INDEX_IJK,ispec)

    duxdxl = xixl*tempx1(INDEX_IJK) + etaxl*tempx2(INDEX_IJK) + gammaxl*tempx3(INDEX_IJK)
    duxdyl = xiyl*tempx1(INDEX_IJK) + etayl*tempx2(INDEX_IJK) + gammayl*tempx3(INDEX_IJK)
    duxdzl = xizl*tempx1(INDEX_IJK) + etazl*tempx2(INDEX_IJK) + gammazl*tempx3(INDEX_IJK)

    duydxl = xixl*tempy1(INDEX_IJK) + etaxl*tempy2(INDEX_IJK) + gammaxl*tempy3(INDEX_IJK)
    duydyl = xiyl*tempy1(INDEX_IJK) + etayl*tempy2(INDEX_IJK) + gammayl*tempy3(INDEX_IJK)
    duydzl = xizl*tempy1(INDEX_IJK) + etazl*tempy2(INDEX_IJK) + gammazl*tempy3(INDEX_IJK)

    duzdxl = xixl*tempz1(INDEX_IJK) + etaxl*tempz2(INDEX_IJK) + gammaxl*tempz3(INDEX_IJK)
    duzdyl = xiyl*tempz1(INDEX_IJK) + etayl*tempz2(INDEX_IJK) + gammayl*tempz3(INDEX_IJK)
    duzdzl = xizl*tempz1(INDEX_IJK) + etazl*tempz2(INDEX_IJK) + gammazl*tempz3(INDEX_IJK)

    ! precompute some sums to save CPU time
    duxdxl_plus_duydyl = duxdxl + duydyl
    duxdxl_plus_duzdzl = duxdxl + duzdzl
    duydyl_plus_duzdzl = duydyl + duzdzl
    duxdyl_plus_duydxl = duxdyl + duydxl
    duzdxl_plus_duxdzl = duzdxl + duxdzl
    duzdyl_plus_duydzl = duzdyl + duydzl

    templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
    if (nspec_strain_only == 1) then
      if (ispec == 1) then
        eps_trace_over_3_loc_nplus1(INDEX_IJK,1) = templ
      endif
    else
      eps_trace_over_3_loc_nplus1(INDEX_IJK,ispec) = templ
    endif
    epsilondev_xx_loc_nplus1(INDEX_IJK,ispec) = duxdxl - templ
    epsilondev_yy_loc_nplus1(INDEX_IJK,ispec) = duydyl - templ
    epsilondev_xy_loc_nplus1(INDEX_IJK,ispec) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
    epsilondev_xz_loc_nplus1(INDEX_IJK,ispec) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
    epsilondev_yz_loc_nplus1(INDEX_IJK,ispec) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl

  ENDDO_LOOP_IJK

  contains

!--------------------------------------------------------------------------------------------
!
! matrix-matrix multiplications
!
! subroutines adapted from Deville, Fischer and Mund, High-order methods
! for incompressible fluid flow, Cambridge University Press (2002),
! pages 386 and 389 and Figure 8.3.1
!
!--------------------------------------------------------------------------------------------
!
! note: the matrix-matrix multiplications are used for very small matrices ( 5 x 5 x 5 elements);
!       thus, calling external optimized libraries for these multiplications are in general slower
!
! please leave the routines here to help compilers inlining the code

  subroutine mxm5_3comp_singleA(A,n1,B1,B2,B3,C1,C2,C3,n3)

! 3 different arrays for x/y/z-components, 2-dimensional arrays (25,5)/(5,25), same B matrix for all 3 component arrays

  use constants_solver,only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5) :: A
  real(kind=CUSTOM_REAL),dimension(5,n3) :: B1,B2,B3
  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C1,C2,C3

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
    do i = 1,n1
      C1(i,j) =  A(i,1) * B1(1,j) &
               + A(i,2) * B1(2,j) &
               + A(i,3) * B1(3,j) &
               + A(i,4) * B1(4,j) &
               + A(i,5) * B1(5,j)

      C2(i,j) =  A(i,1) * B2(1,j) &
               + A(i,2) * B2(2,j) &
               + A(i,3) * B2(3,j) &
               + A(i,4) * B2(4,j) &
               + A(i,5) * B2(5,j)

      C3(i,j) =  A(i,1) * B3(1,j) &
               + A(i,2) * B3(2,j) &
               + A(i,3) * B3(3,j) &
               + A(i,4) * B3(4,j) &
               + A(i,5) * B3(5,j)
    enddo
  enddo

  end subroutine mxm5_3comp_singleA


!--------------------------------------------------------------------------------------------

  subroutine mxm5_3comp_singleB(A1,A2,A3,n1,B,C1,C2,C3,n3)

! 3 different arrays for x/y/z-components, 2-dimensional arrays (25,5)/(5,25), same B matrix for all 3 component arrays

  use constants_solver,only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(5,n3) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C1,C2,C3

  ! local parameters
  integer :: i,j

  ! matrix-matrix multiplication
  do j = 1,n3
    do i = 1,n1
      C1(i,j) =  A1(i,1) * B(1,j) &
               + A1(i,2) * B(2,j) &
               + A1(i,3) * B(3,j) &
               + A1(i,4) * B(4,j) &
               + A1(i,5) * B(5,j)

      C2(i,j) =  A2(i,1) * B(1,j) &
               + A2(i,2) * B(2,j) &
               + A2(i,3) * B(3,j) &
               + A2(i,4) * B(4,j) &
               + A2(i,5) * B(5,j)

      C3(i,j) =  A3(i,1) * B(1,j) &
               + A3(i,2) * B(2,j) &
               + A3(i,3) * B(3,j) &
               + A3(i,4) * B(4,j) &
               + A3(i,5) * B(5,j)
    enddo
  enddo

  end subroutine mxm5_3comp_singleB


!--------------------------------------------------------------------------------------------

  subroutine mxm5_3comp_3dmat_singleB(A1,A2,A3,n1,B,n2,C1,C2,C3,n3)

! 3 different arrays for x/y/z-components, 3-dimensional arrays (5,5,5), same B matrix for all 3 component arrays

  use constants_solver,only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,5,n3) :: A1,A2,A3
  real(kind=CUSTOM_REAL),dimension(5,n2) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3) :: C1,C2,C3

  ! local parameters
  integer :: i,j,k

  ! matrix-matrix multiplication
  do j = 1,n2
    do i = 1,n1
      ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
      do k = 1,n3
        C1(i,j,k) =  A1(i,1,k) * B(1,j) &
                   + A1(i,2,k) * B(2,j) &
                   + A1(i,3,k) * B(3,j) &
                   + A1(i,4,k) * B(4,j) &
                   + A1(i,5,k) * B(5,j)

        C2(i,j,k) =  A2(i,1,k) * B(1,j) &
                   + A2(i,2,k) * B(2,j) &
                   + A2(i,3,k) * B(3,j) &
                   + A2(i,4,k) * B(4,j) &
                   + A2(i,5,k) * B(5,j)

        C3(i,j,k) =  A3(i,1,k) * B(1,j) &
                   + A3(i,2,k) * B(2,j) &
                   + A3(i,3,k) * B(3,j) &
                   + A3(i,4,k) * B(4,j) &
                   + A3(i,5,k) * B(5,j)
      enddo
    enddo
  enddo

  end subroutine mxm5_3comp_3dmat_singleB

  end subroutine compute_element_strain_att_Dev


!
!--------------------------------------------------------------------------------------------
!

  subroutine compute_element_strain_att_noDev(ispec,nglob,nspec, &
                                              displ,veloc,deltat, &
                                              ibool, &
                                              hprime_xx,hprime_yy,hprime_zz, &
                                              xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                                              epsilondev_xx_loc_nplus1, &
                                              epsilondev_yy_loc_nplus1, &
                                              epsilondev_xy_loc_nplus1, &
                                              epsilondev_xz_loc_nplus1, &
                                              epsilondev_yz_loc_nplus1, &
                                              nspec_strain_only,eps_trace_over_3_loc_nplus1)

  use constants

  implicit none

  integer :: ispec,nglob,nspec
  real(kind=CUSTOM_REAL) :: deltat

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob) :: displ,veloc

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  ! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX),intent(in) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY),intent(in) :: hprime_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ),intent(in) :: hprime_zz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: epsilondev_xx_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: epsilondev_yy_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: epsilondev_xy_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: epsilondev_xz_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: epsilondev_yz_loc_nplus1

  integer :: nspec_strain_only
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_strain_only) :: eps_trace_over_3_loc_nplus1

  ! local variable
  integer :: i,j,k,l,iglob
  real(kind=CUSTOM_REAL) :: templ

  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) tempz1l,tempz2l,tempz3l

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL) hp1,hp2,hp3
  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl
  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          dummyx_loc(i,j,k) = displ(1,iglob) + deltat * veloc(1,iglob)
          dummyy_loc(i,j,k) = displ(2,iglob) + deltat * veloc(2,iglob)
          dummyz_loc(i,j,k) = displ(3,iglob) + deltat * veloc(3,iglob)
      enddo
    enddo
  enddo

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        tempx1l = 0._CUSTOM_REAL
        tempx2l = 0._CUSTOM_REAL
        tempx3l = 0._CUSTOM_REAL

        tempy1l = 0._CUSTOM_REAL
        tempy2l = 0._CUSTOM_REAL
        tempy3l = 0._CUSTOM_REAL

        tempz1l = 0._CUSTOM_REAL
        tempz2l = 0._CUSTOM_REAL
        tempz3l = 0._CUSTOM_REAL

        do l = 1,NGLLX
          hp1 = hprime_xx(i,l)
          tempx1l = tempx1l + dummyx_loc(l,j,k)*hp1
          tempy1l = tempy1l + dummyy_loc(l,j,k)*hp1
          tempz1l = tempz1l + dummyz_loc(l,j,k)*hp1
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l = 1,NGLLY
          hp2 = hprime_yy(j,l)
          tempx2l = tempx2l + dummyx_loc(i,l,k)*hp2
          tempy2l = tempy2l + dummyy_loc(i,l,k)*hp2
          tempz2l = tempz2l + dummyz_loc(i,l,k)*hp2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l = 1,NGLLZ
          hp3 = hprime_zz(k,l)
          tempx3l = tempx3l + dummyx_loc(i,j,l)*hp3
          tempy3l = tempy3l + dummyy_loc(i,j,l)*hp3
          tempz3l = tempz3l + dummyz_loc(i,j,l)*hp3
        enddo

        ! get derivatives of ux, uy and uz with respect to x, y and z
        xixl = xix(i,j,k,ispec)
        xiyl = xiy(i,j,k,ispec)
        xizl = xiz(i,j,k,ispec)
        etaxl = etax(i,j,k,ispec)
        etayl = etay(i,j,k,ispec)
        etazl = etaz(i,j,k,ispec)
        gammaxl = gammax(i,j,k,ispec)
        gammayl = gammay(i,j,k,ispec)
        gammazl = gammaz(i,j,k,ispec)

        duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
        duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
        duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

        duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
        duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
        duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

        duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
        duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
        duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l

        ! precompute some sums to save CPU time
        duxdxl_plus_duydyl = duxdxl + duydyl
        duxdxl_plus_duzdzl = duxdxl + duzdzl
        duydyl_plus_duzdzl = duydyl + duzdzl
        duxdyl_plus_duydxl = duxdyl + duydxl
        duzdxl_plus_duxdzl = duzdxl + duxdzl
        duzdyl_plus_duydzl = duzdyl + duydzl

        templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
        if (nspec_strain_only == 1) then
          if (ispec == 1) then
            eps_trace_over_3_loc_nplus1(i,j,k,1) = templ
          endif
        else
          eps_trace_over_3_loc_nplus1(i,j,k,ispec) = templ
        endif
        epsilondev_xx_loc_nplus1(i,j,k,ispec) = duxdxl - templ
        epsilondev_yy_loc_nplus1(i,j,k,ispec) = duydyl - templ
        epsilondev_xy_loc_nplus1(i,j,k,ispec) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
        epsilondev_xz_loc_nplus1(i,j,k,ispec) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
        epsilondev_yz_loc_nplus1(i,j,k,ispec) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl
      enddo
    enddo
  enddo

  end subroutine compute_element_strain_att_noDev

