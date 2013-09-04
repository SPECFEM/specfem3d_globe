!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            August 2013
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



  subroutine compute_element_strain_undo_att_Dev(ispec,nglob,nspec, &
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

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ),intent(out) :: epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: eps_trace_over_3_loc

!  local variable
  integer :: i,j,k,iglob
  real(kind=CUSTOM_REAL) :: templ
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,m2) :: B1_m1_m2_5points,B2_m1_m2_5points,B3_m1_m2_5points

  equivalence(dummyx_loc,B1_m1_m2_5points)
  equivalence(dummyy_loc,B2_m1_m2_5points)
  equivalence(dummyz_loc,B3_m1_m2_5points)

  real(kind=CUSTOM_REAL), dimension(m1,m2) :: C1_m1_m2_5points,C2_m1_m2_5points,C3_m1_m2_5points

  equivalence(tempx1,C1_m1_m2_5points)
  equivalence(tempy1,C2_m1_m2_5points)
  equivalence(tempz1,C3_m1_m2_5points)

  real(kind=CUSTOM_REAL), dimension(m2,m1) :: &
    C1_mxm_m2_m1_5points,C2_mxm_m2_m1_5points,C3_mxm_m2_m1_5points

  equivalence(tempx3,C1_mxm_m2_m1_5points)
  equivalence(tempy3,C2_mxm_m2_m1_5points)
  equivalence(tempz3,C3_mxm_m2_m1_5points)

  real(kind=CUSTOM_REAL), dimension(m2,NGLLX) :: &
    A1_mxm_m2_m1_5points,A2_mxm_m2_m1_5points,A3_mxm_m2_m1_5points

  equivalence(dummyx_loc,A1_mxm_m2_m1_5points)
  equivalence(dummyy_loc,A2_mxm_m2_m1_5points)
  equivalence(dummyz_loc,A3_mxm_m2_m1_5points)

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duydyl,duzdzl,duxdyl,duydxl,duzdxl,duxdzl,duzdyl,duydzl,&
                         duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl,&
                         duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#endif

#ifdef FORCE_VECTORIZATION
  do ijk=1,NGLLCUBE
    iglob = ibool(ijk,1,1,ispec)
    dummyx_loc(ijk,1,1) = displ(1,iglob)
    dummyy_loc(ijk,1,1) = displ(2,iglob)
    dummyz_loc(ijk,1,1) = displ(3,iglob)
  enddo
#else
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        iglob = ibool(i,j,k,ispec)
        dummyx_loc(i,j,k) = displ(1,iglob)
        dummyy_loc(i,j,k) = displ(2,iglob)
        dummyz_loc(i,j,k) = displ(3,iglob)
      enddo
    enddo
  enddo
#endif

  do j=1,m2
    do i=1,m1
      C1_m1_m2_5points(i,j) = hprime_xx(i,1)*B1_m1_m2_5points(1,j) + &
                              hprime_xx(i,2)*B1_m1_m2_5points(2,j) + &
                              hprime_xx(i,3)*B1_m1_m2_5points(3,j) + &
                              hprime_xx(i,4)*B1_m1_m2_5points(4,j) + &
                              hprime_xx(i,5)*B1_m1_m2_5points(5,j)

      C2_m1_m2_5points(i,j) = hprime_xx(i,1)*B2_m1_m2_5points(1,j) + &
                              hprime_xx(i,2)*B2_m1_m2_5points(2,j) + &
                              hprime_xx(i,3)*B2_m1_m2_5points(3,j) + &
                              hprime_xx(i,4)*B2_m1_m2_5points(4,j) + &
                              hprime_xx(i,5)*B2_m1_m2_5points(5,j)

      C3_m1_m2_5points(i,j) = hprime_xx(i,1)*B3_m1_m2_5points(1,j) + &
                              hprime_xx(i,2)*B3_m1_m2_5points(2,j) + &
                              hprime_xx(i,3)*B3_m1_m2_5points(3,j) + &
                              hprime_xx(i,4)*B3_m1_m2_5points(4,j) + &
                              hprime_xx(i,5)*B3_m1_m2_5points(5,j)
    enddo
  enddo

  do j=1,m1
    do i=1,m1
      ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
      do k = 1,NGLLX
        tempx2(i,j,k) = dummyx_loc(i,1,k)*hprime_xxT(1,j) + &
                        dummyx_loc(i,2,k)*hprime_xxT(2,j) + &
                        dummyx_loc(i,3,k)*hprime_xxT(3,j) + &
                        dummyx_loc(i,4,k)*hprime_xxT(4,j) + &
                        dummyx_loc(i,5,k)*hprime_xxT(5,j)

        tempy2(i,j,k) = dummyy_loc(i,1,k)*hprime_xxT(1,j) + &
                        dummyy_loc(i,2,k)*hprime_xxT(2,j) + &
                        dummyy_loc(i,3,k)*hprime_xxT(3,j) + &
                        dummyy_loc(i,4,k)*hprime_xxT(4,j) + &
                        dummyy_loc(i,5,k)*hprime_xxT(5,j)

        tempz2(i,j,k) = dummyz_loc(i,1,k)*hprime_xxT(1,j) + &
                        dummyz_loc(i,2,k)*hprime_xxT(2,j) + &
                        dummyz_loc(i,3,k)*hprime_xxT(3,j) + &
                        dummyz_loc(i,4,k)*hprime_xxT(4,j) + &
                        dummyz_loc(i,5,k)*hprime_xxT(5,j)
      enddo
    enddo
  enddo

  do j=1,m1
    do i=1,m2
      C1_mxm_m2_m1_5points(i,j) = A1_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                  A1_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                  A1_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                  A1_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                  A1_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)

      C2_mxm_m2_m1_5points(i,j) = A2_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                  A2_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                  A2_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                  A2_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                  A2_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)

      C3_mxm_m2_m1_5points(i,j) = A3_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                  A3_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                  A3_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                  A3_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                  A3_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)
    enddo
  enddo

#ifdef FORCE_VECTORIZATION
  do ijk=1,NGLLCUBE
    ! get derivatives of ux, uy and uz with respect to x, y and z
    xixl = xix(ijk,1,1,ispec)
    xiyl = xiy(ijk,1,1,ispec)
    xizl = xiz(ijk,1,1,ispec)
    etaxl = etax(ijk,1,1,ispec)
    etayl = etay(ijk,1,1,ispec)
    etazl = etaz(ijk,1,1,ispec)
    gammaxl = gammax(ijk,1,1,ispec)
    gammayl = gammay(ijk,1,1,ispec)
    gammazl = gammaz(ijk,1,1,ispec)

    ! compute the jacobian
    jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                  - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                  + xizl*(etaxl*gammayl-etayl*gammaxl))

    duxdxl = xixl*tempx1(ijk,1,1) + etaxl*tempx2(ijk,1,1) + gammaxl*tempx3(ijk,1,1)
    duxdyl = xiyl*tempx1(ijk,1,1) + etayl*tempx2(ijk,1,1) + gammayl*tempx3(ijk,1,1)
    duxdzl = xizl*tempx1(ijk,1,1) + etazl*tempx2(ijk,1,1) + gammazl*tempx3(ijk,1,1)

    duydxl = xixl*tempy1(ijk,1,1) + etaxl*tempy2(ijk,1,1) + gammaxl*tempy3(ijk,1,1)
    duydyl = xiyl*tempy1(ijk,1,1) + etayl*tempy2(ijk,1,1) + gammayl*tempy3(ijk,1,1)
    duydzl = xizl*tempy1(ijk,1,1) + etazl*tempy2(ijk,1,1) + gammazl*tempy3(ijk,1,1)

    duzdxl = xixl*tempz1(ijk,1,1) + etaxl*tempz2(ijk,1,1) + gammaxl*tempz3(ijk,1,1)
    duzdyl = xiyl*tempz1(ijk,1,1) + etayl*tempz2(ijk,1,1) + gammayl*tempz3(ijk,1,1)
    duzdzl = xizl*tempz1(ijk,1,1) + etazl*tempz2(ijk,1,1) + gammazl*tempz3(ijk,1,1)

    ! precompute some sums to save CPU time
    duxdxl_plus_duydyl = duxdxl + duydyl
    duxdxl_plus_duzdzl = duxdxl + duzdzl
    duydyl_plus_duzdzl = duydyl + duzdzl
    duxdyl_plus_duydxl = duxdyl + duydxl
    duzdxl_plus_duxdzl = duzdxl + duxdzl
    duzdyl_plus_duydzl = duzdyl + duydzl

    ! strains
    templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
    eps_trace_over_3_loc(ijk,1,1) = templ
    epsilondev_loc(1,ijk,1,1) = duxdxl - templ
    epsilondev_loc(2,ijk,1,1) = duydyl - templ
    epsilondev_loc(3,ijk,1,1) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
    epsilondev_loc(4,ijk,1,1) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
    epsilondev_loc(5,ijk,1,1) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl
  enddo
#else
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
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

        ! compute the jacobian
        jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                      - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                      + xizl*(etaxl*gammayl-etayl*gammaxl))

        duxdxl = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
        duxdyl = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
        duxdzl = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)

        duydxl = xixl*tempy1(i,j,k) + etaxl*tempy2(i,j,k) + gammaxl*tempy3(i,j,k)
        duydyl = xiyl*tempy1(i,j,k) + etayl*tempy2(i,j,k) + gammayl*tempy3(i,j,k)
        duydzl = xizl*tempy1(i,j,k) + etazl*tempy2(i,j,k) + gammazl*tempy3(i,j,k)

        duzdxl = xixl*tempz1(i,j,k) + etaxl*tempz2(i,j,k) + gammaxl*tempz3(i,j,k)
        duzdyl = xiyl*tempz1(i,j,k) + etayl*tempz2(i,j,k) + gammayl*tempz3(i,j,k)
        duzdzl = xizl*tempz1(i,j,k) + etazl*tempz2(i,j,k) + gammazl*tempz3(i,j,k)

        ! precompute some sums to save CPU time
        duxdxl_plus_duydyl = duxdxl + duydyl
        duxdxl_plus_duzdzl = duxdxl + duzdzl
        duydyl_plus_duzdzl = duydyl + duzdzl
        duxdyl_plus_duydxl = duxdyl + duydxl
        duzdxl_plus_duxdzl = duzdxl + duxdzl
        duzdyl_plus_duydzl = duzdyl + duydzl

        ! strains
        templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
        eps_trace_over_3_loc(i,j,k) = templ
        epsilondev_loc(1,i,j,k) = duxdxl - templ
        epsilondev_loc(2,i,j,k) = duydyl - templ
        epsilondev_loc(3,i,j,k) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
        epsilondev_loc(4,i,j,k) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
        epsilondev_loc(5,i,j,k) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl
      enddo
    enddo
  enddo
#endif

  end subroutine compute_element_strain_undo_att_Dev


!
!--------------------------------------------------------------------------------------------
!

  subroutine compute_element_strain_undo_att_noDev(ispec,nglob,nspec, &
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

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ),intent(out) :: epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: eps_trace_over_3_loc

  ! local parameters
  integer :: iglob
  integer :: i,j,k,l

  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) tempz1l,tempz2l,tempz3l

  real(kind=CUSTOM_REAL) hp1,hp2,hp3
  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl
  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX

        tempx1l = 0._CUSTOM_REAL
        tempx2l = 0._CUSTOM_REAL
        tempx3l = 0._CUSTOM_REAL

        tempy1l = 0._CUSTOM_REAL
        tempy2l = 0._CUSTOM_REAL
        tempy3l = 0._CUSTOM_REAL

        tempz1l = 0._CUSTOM_REAL
        tempz2l = 0._CUSTOM_REAL
        tempz3l = 0._CUSTOM_REAL

        do l=1,NGLLX
          hp1 = hprime_xx(i,l)
          iglob = ibool(l,j,k,ispec)
          tempx1l = tempx1l + displ(1,iglob)*hp1
          tempy1l = tempy1l + displ(2,iglob)*hp1
          tempz1l = tempz1l + displ(3,iglob)*hp1
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
          hp2 = hprime_yy(j,l)
          iglob = ibool(i,l,k,ispec)
          tempx2l = tempx2l + displ(1,iglob)*hp2
          tempy2l = tempy2l + displ(2,iglob)*hp2
          tempz2l = tempz2l + displ(3,iglob)*hp2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
          hp3 = hprime_zz(k,l)
          iglob = ibool(i,j,l,ispec)
          tempx3l = tempx3l + displ(1,iglob)*hp3
          tempy3l = tempy3l + displ(2,iglob)*hp3
          tempz3l = tempz3l + displ(3,iglob)*hp3
        enddo

!         get derivatives of ux, uy and uz with respect to x, y and z
        xixl = xix(i,j,k,ispec)
        xiyl = xiy(i,j,k,ispec)
        xizl = xiz(i,j,k,ispec)
        etaxl = etax(i,j,k,ispec)
        etayl = etay(i,j,k,ispec)
        etazl = etaz(i,j,k,ispec)
        gammaxl = gammax(i,j,k,ispec)
        gammayl = gammay(i,j,k,ispec)
        gammazl = gammaz(i,j,k,ispec)

! compute the jacobian
        jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                      - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                      + xizl*(etaxl*gammayl-etayl*gammaxl))

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
        epsilondev_loc(1,i,j,k) = duxdxl - eps_trace_over_3_loc(i,j,k)
        epsilondev_loc(2,i,j,k) = duydyl - eps_trace_over_3_loc(i,j,k)
        epsilondev_loc(3,i,j,k) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
        epsilondev_loc(4,i,j,k) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
        epsilondev_loc(5,i,j,k) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl

      enddo ! NGLLX
    enddo ! NGLLY
  enddo ! NGLLZ

  end subroutine compute_element_strain_undo_att_noDev


!--------------------------------------------------------------------------------------------
!
! strain separated into single xx,yy,xy,xz,yz-component arrays
!
!--------------------------------------------------------------------------------------------
!

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
                                           eps_trace_over_3_loc_nplus1)

  use constants

  implicit none

  integer :: ispec,nglob,nspec
  real(kind=CUSTOM_REAL) :: deltat

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob) :: displ,veloc

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprime_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  !real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_xx_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_yy_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_xy_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_xz_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_yz_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: eps_trace_over_3_loc_nplus1

  ! local variable
  integer :: i,j,k,iglob
  real(kind=CUSTOM_REAL) :: templ
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,m2) :: B1_m1_m2_5points,B2_m1_m2_5points,B3_m1_m2_5points
  equivalence(dummyx_loc,B1_m1_m2_5points)
  equivalence(dummyy_loc,B2_m1_m2_5points)
  equivalence(dummyz_loc,B3_m1_m2_5points)
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: C1_m1_m2_5points,C2_m1_m2_5points,C3_m1_m2_5points
  equivalence(tempx1,C1_m1_m2_5points)
  equivalence(tempy1,C2_m1_m2_5points)
  equivalence(tempz1,C3_m1_m2_5points)
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: &
    C1_mxm_m2_m1_5points,C2_mxm_m2_m1_5points,C3_mxm_m2_m1_5points
  equivalence(tempx3,C1_mxm_m2_m1_5points)
  equivalence(tempy3,C2_mxm_m2_m1_5points)
  equivalence(tempz3,C3_mxm_m2_m1_5points)
  real(kind=CUSTOM_REAL), dimension(m2,NGLLX) :: &
    A1_mxm_m2_m1_5points,A2_mxm_m2_m1_5points,A3_mxm_m2_m1_5points
  equivalence(dummyx_loc,A1_mxm_m2_m1_5points)
  equivalence(dummyy_loc,A2_mxm_m2_m1_5points)
  equivalence(dummyz_loc,A3_mxm_m2_m1_5points)

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duydyl,duzdzl,duxdyl,duydxl,duzdxl,duxdzl,duzdyl,duydzl,&
                         duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl,&
                         duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#endif

#ifdef FORCE_VECTORIZATION
  do ijk=1,NGLLCUBE
      iglob = ibool(ijk,1,1,ispec)
      dummyx_loc(ijk,1,1) = displ(1,iglob) + deltat * veloc(1,iglob)
      dummyy_loc(ijk,1,1) = displ(2,iglob) + deltat * veloc(2,iglob)
      dummyz_loc(ijk,1,1) = displ(3,iglob) + deltat * veloc(3,iglob)
  enddo
#else
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          dummyx_loc(i,j,k) = displ(1,iglob) + deltat * veloc(1,iglob)
          dummyy_loc(i,j,k) = displ(2,iglob) + deltat * veloc(2,iglob)
          dummyz_loc(i,j,k) = displ(3,iglob) + deltat * veloc(3,iglob)
      enddo
    enddo
  enddo
#endif

  do j=1,m2
    do i=1,m1
      C1_m1_m2_5points(i,j) = hprime_xx(i,1)*B1_m1_m2_5points(1,j) + &
                              hprime_xx(i,2)*B1_m1_m2_5points(2,j) + &
                              hprime_xx(i,3)*B1_m1_m2_5points(3,j) + &
                              hprime_xx(i,4)*B1_m1_m2_5points(4,j) + &
                              hprime_xx(i,5)*B1_m1_m2_5points(5,j)

      C2_m1_m2_5points(i,j) = hprime_xx(i,1)*B2_m1_m2_5points(1,j) + &
                              hprime_xx(i,2)*B2_m1_m2_5points(2,j) + &
                              hprime_xx(i,3)*B2_m1_m2_5points(3,j) + &
                              hprime_xx(i,4)*B2_m1_m2_5points(4,j) + &
                              hprime_xx(i,5)*B2_m1_m2_5points(5,j)

      C3_m1_m2_5points(i,j) = hprime_xx(i,1)*B3_m1_m2_5points(1,j) + &
                              hprime_xx(i,2)*B3_m1_m2_5points(2,j) + &
                              hprime_xx(i,3)*B3_m1_m2_5points(3,j) + &
                              hprime_xx(i,4)*B3_m1_m2_5points(4,j) + &
                              hprime_xx(i,5)*B3_m1_m2_5points(5,j)
    enddo
  enddo

  do j=1,m1
    do i=1,m1
      ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
      do k = 1,NGLLX
        tempx2(i,j,k) = dummyx_loc(i,1,k)*hprime_xxT(1,j) + &
                        dummyx_loc(i,2,k)*hprime_xxT(2,j) + &
                        dummyx_loc(i,3,k)*hprime_xxT(3,j) + &
                        dummyx_loc(i,4,k)*hprime_xxT(4,j) + &
                        dummyx_loc(i,5,k)*hprime_xxT(5,j)

        tempy2(i,j,k) = dummyy_loc(i,1,k)*hprime_xxT(1,j) + &
                        dummyy_loc(i,2,k)*hprime_xxT(2,j) + &
                        dummyy_loc(i,3,k)*hprime_xxT(3,j) + &
                        dummyy_loc(i,4,k)*hprime_xxT(4,j) + &
                        dummyy_loc(i,5,k)*hprime_xxT(5,j)

        tempz2(i,j,k) = dummyz_loc(i,1,k)*hprime_xxT(1,j) + &
                        dummyz_loc(i,2,k)*hprime_xxT(2,j) + &
                        dummyz_loc(i,3,k)*hprime_xxT(3,j) + &
                        dummyz_loc(i,4,k)*hprime_xxT(4,j) + &
                        dummyz_loc(i,5,k)*hprime_xxT(5,j)
      enddo
    enddo
  enddo

  do j=1,m1
    do i=1,m2
      C1_mxm_m2_m1_5points(i,j) = A1_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                  A1_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                  A1_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                  A1_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                  A1_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)

      C2_mxm_m2_m1_5points(i,j) = A2_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                  A2_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                  A2_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                  A2_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                  A2_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)

      C3_mxm_m2_m1_5points(i,j) = A3_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                  A3_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                  A3_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                  A3_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                  A3_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)
    enddo
  enddo

#ifdef FORCE_VECTORIZATION
  do ijk=1,NGLLCUBE
    ! get derivatives of ux, uy and uz with respect to x, y and z
    xixl = xix(ijk,1,1,ispec)
    xiyl = xiy(ijk,1,1,ispec)
    xizl = xiz(ijk,1,1,ispec)
    etaxl = etax(ijk,1,1,ispec)
    etayl = etay(ijk,1,1,ispec)
    etazl = etaz(ijk,1,1,ispec)
    gammaxl = gammax(ijk,1,1,ispec)
    gammayl = gammay(ijk,1,1,ispec)
    gammazl = gammaz(ijk,1,1,ispec)

    ! compute the jacobian
    jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                  - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                  + xizl*(etaxl*gammayl-etayl*gammaxl))

    duxdxl = xixl*tempx1(ijk,1,1) + etaxl*tempx2(ijk,1,1) + gammaxl*tempx3(ijk,1,1)
    duxdyl = xiyl*tempx1(ijk,1,1) + etayl*tempx2(ijk,1,1) + gammayl*tempx3(ijk,1,1)
    duxdzl = xizl*tempx1(ijk,1,1) + etazl*tempx2(ijk,1,1) + gammazl*tempx3(ijk,1,1)

    duydxl = xixl*tempy1(ijk,1,1) + etaxl*tempy2(ijk,1,1) + gammaxl*tempy3(ijk,1,1)
    duydyl = xiyl*tempy1(ijk,1,1) + etayl*tempy2(ijk,1,1) + gammayl*tempy3(ijk,1,1)
    duydzl = xizl*tempy1(ijk,1,1) + etazl*tempy2(ijk,1,1) + gammazl*tempy3(ijk,1,1)

    duzdxl = xixl*tempz1(ijk,1,1) + etaxl*tempz2(ijk,1,1) + gammaxl*tempz3(ijk,1,1)
    duzdyl = xiyl*tempz1(ijk,1,1) + etayl*tempz2(ijk,1,1) + gammayl*tempz3(ijk,1,1)
    duzdzl = xizl*tempz1(ijk,1,1) + etazl*tempz2(ijk,1,1) + gammazl*tempz3(ijk,1,1)

    ! precompute some sums to save CPU time
    duxdxl_plus_duydyl = duxdxl + duydyl
    duxdxl_plus_duzdzl = duxdxl + duzdzl
    duydyl_plus_duzdzl = duydyl + duzdzl
    duxdyl_plus_duydxl = duxdyl + duydxl
    duzdxl_plus_duxdzl = duzdxl + duxdzl
    duzdyl_plus_duydzl = duzdyl + duydzl

    templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
    eps_trace_over_3_loc_nplus1 = templ
    epsilondev_xx_loc_nplus1(ijk,1,1) = duxdxl - templ
    epsilondev_yy_loc_nplus1(ijk,1,1) = duydyl - templ
    epsilondev_xy_loc_nplus1(ijk,1,1) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
    epsilondev_xz_loc_nplus1(ijk,1,1) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
    epsilondev_yz_loc_nplus1(ijk,1,1) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl
  enddo
#else
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
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

        ! compute the jacobian
        jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                      - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                      + xizl*(etaxl*gammayl-etayl*gammaxl))

        duxdxl = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
        duxdyl = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
        duxdzl = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)

        duydxl = xixl*tempy1(i,j,k) + etaxl*tempy2(i,j,k) + gammaxl*tempy3(i,j,k)
        duydyl = xiyl*tempy1(i,j,k) + etayl*tempy2(i,j,k) + gammayl*tempy3(i,j,k)
        duydzl = xizl*tempy1(i,j,k) + etazl*tempy2(i,j,k) + gammazl*tempy3(i,j,k)

        duzdxl = xixl*tempz1(i,j,k) + etaxl*tempz2(i,j,k) + gammaxl*tempz3(i,j,k)
        duzdyl = xiyl*tempz1(i,j,k) + etayl*tempz2(i,j,k) + gammayl*tempz3(i,j,k)
        duzdzl = xizl*tempz1(i,j,k) + etazl*tempz2(i,j,k) + gammazl*tempz3(i,j,k)

        ! precompute some sums to save CPU time
        duxdxl_plus_duydyl = duxdxl + duydyl
        duxdxl_plus_duzdzl = duxdxl + duzdzl
        duydyl_plus_duzdzl = duydyl + duzdzl
        duxdyl_plus_duydxl = duxdyl + duydxl
        duzdxl_plus_duxdzl = duzdxl + duxdzl
        duzdyl_plus_duydzl = duzdyl + duydzl

        templ = ONE_THIRD * (duxdxl + duydyl + duzdzl)
        eps_trace_over_3_loc_nplus1 = templ
        epsilondev_xx_loc_nplus1(i,j,k) = duxdxl - templ
        epsilondev_yy_loc_nplus1(i,j,k) = duydyl - templ
        epsilondev_xy_loc_nplus1(i,j,k) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
        epsilondev_xz_loc_nplus1(i,j,k) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
        epsilondev_yz_loc_nplus1(i,j,k) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl
      enddo
    enddo
  enddo
#endif

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
                                             eps_trace_over_3_loc_nplus1)

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

  !real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_xx_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_yy_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_xy_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_xz_loc_nplus1
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: epsilondev_yz_loc_nplus1

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: eps_trace_over_3_loc_nplus1

  ! local variable
  integer :: i,j,k,l,iglob
  real(kind=CUSTOM_REAL) :: templ

  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) tempz1l,tempz2l,tempz3l

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL) hp1,hp2,hp3
  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl
  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          dummyx_loc(i,j,k) = displ(1,iglob) + deltat * veloc(1,iglob)
          dummyy_loc(i,j,k) = displ(2,iglob) + deltat * veloc(2,iglob)
          dummyz_loc(i,j,k) = displ(3,iglob) + deltat * veloc(3,iglob)
      enddo
    enddo
  enddo

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX

        tempx1l = 0._CUSTOM_REAL
        tempx2l = 0._CUSTOM_REAL
        tempx3l = 0._CUSTOM_REAL

        tempy1l = 0._CUSTOM_REAL
        tempy2l = 0._CUSTOM_REAL
        tempy3l = 0._CUSTOM_REAL

        tempz1l = 0._CUSTOM_REAL
        tempz2l = 0._CUSTOM_REAL
        tempz3l = 0._CUSTOM_REAL

        do l=1,NGLLX
          hp1 = hprime_xx(i,l)
          tempx1l = tempx1l + dummyx_loc(l,j,k)*hp1
          tempy1l = tempy1l + dummyy_loc(l,j,k)*hp1
          tempz1l = tempz1l + dummyz_loc(l,j,k)*hp1
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
          hp2 = hprime_yy(j,l)
          tempx2l = tempx2l + dummyx_loc(i,l,k)*hp2
          tempy2l = tempy2l + dummyy_loc(i,l,k)*hp2
          tempz2l = tempz2l + dummyz_loc(i,l,k)*hp2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
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

        ! compute the jacobian
        jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                      - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                      + xizl*(etaxl*gammayl-etayl*gammaxl))

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
        eps_trace_over_3_loc_nplus1 = templ
        epsilondev_xx_loc_nplus1(i,j,k) = duxdxl - templ
        epsilondev_yy_loc_nplus1(i,j,k) = duydyl - templ
        epsilondev_xy_loc_nplus1(i,j,k) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl
        epsilondev_xz_loc_nplus1(i,j,k) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl
        epsilondev_yz_loc_nplus1(i,j,k) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl
      enddo
    enddo
  enddo

 end subroutine compute_element_strain_att_noDev

