
! generate sed command to inline matrix products
! use on file below
! then cut and paste resulting code in compute_forces*_inline.f90

  program generate_sed

  implicit none

  integer, parameter :: NGLL = 5

  integer i

  do i=1,NGLL
    write(*,130) i,i
    write(*,131) i,i
    write(*,140) i,i
    write(*,141) i,i
  enddo

  130  format('sed -e ''1,$s/aa/',i1.1,'/g'' < inline_scalar_first > file_scalar_first_',i1.1)
  131  format('sed -e ''1,$s/aa/',i1.1,'/g'' < inline_scalar_second > file_scalar_second_',i1.1)

  140  format('sed -e ''1,$s/aa/',i1.1,'/g'' < inline_vector_first > file_vector_first_',i1.1)
  141  format('sed -e ''1,$s/aa/',i1.1,'/g'' < inline_vector_second > file_vector_second_',i1.1)

  end program generate_sed

!----------------
!---------------- source program below for matrix products
!----------------

! scalar is for outer_core, vector is for crust_mantle and inner_core

!!--------------------- first scalar below ----------------------------
!
!!---
!!--- ij is actually jk here
!    tempx1loc(aa,ij,1) = &
!        disploc(1,ij,1)*hprime_xx(1,aa) + disploc(2,ij,1)*hprime_xx(2,aa) + &
!        disploc(3,ij,1)*hprime_xx(3,aa) + disploc(4,ij,1)*hprime_xx(4,aa) + &
!        disploc(5,ij,1)*hprime_xx(5,aa)
!
!!---
!    tempx2loc(i,aa,k) = &
!        disploc(i,1,k)*hprime_yy(1,aa) + disploc(i,2,k)*hprime_yy(2,aa) + &
!        disploc(i,3,k)*hprime_yy(3,aa) + disploc(i,4,k)*hprime_yy(4,aa) + &
!        disploc(i,5,k)*hprime_yy(5,aa)
!
!!---
!    tempx3loc(ij,1,aa) = &
!        disploc(ij,1,1)*hprime_zz(1,aa) + disploc(ij,1,2)*hprime_zz(2,aa) + &
!        disploc(ij,1,3)*hprime_zz(3,aa) + disploc(ij,1,4)*hprime_zz(4,aa) + &
!        disploc(ij,1,5)*hprime_zz(5,aa)
!
!
!!--------------------- second scalar below ----------------------------
!
!!---
!!--- ij is actually jk below
!  tempx1linline(aa,ij,1) = &
!    tempx1(1,ij,1)*hprimewgll_xx(aa,1) + tempx1(2,ij,1)*hprimewgll_xx(aa,2) + &
!    tempx1(3,ij,1)*hprimewgll_xx(aa,3) + tempx1(4,ij,1)*hprimewgll_xx(aa,4) + &
!    tempx1(5,ij,1)*hprimewgll_xx(aa,5)
!
!!---
!  tempx2linline(i,aa,k) = &
!    tempx2(i,1,k)*hprimewgll_yy(aa,1) + tempx2(i,2,k)*hprimewgll_yy(aa,2) + &
!    tempx2(i,3,k)*hprimewgll_yy(aa,3) + tempx2(i,4,k)*hprimewgll_yy(aa,4) + &
!    tempx2(i,5,k)*hprimewgll_yy(aa,5)
!
!!---
!  tempx3linline(ij,1,aa) = &
!    tempx3(ij,1,1)*hprimewgll_zz(aa,1) + tempx3(ij,1,2)*hprimewgll_zz(aa,2) + &
!    tempx3(ij,1,3)*hprimewgll_zz(aa,3) + tempx3(ij,1,4)*hprimewgll_zz(aa,4) + &
!    tempx3(ij,1,5)*hprimewgll_zz(aa,5)
!
!
!!--------------------- first vector below ----------------------------
!
!!---
!  temp1(d,aa,j,k) = &
!    disploc(d,1,j,k)*hprime_xx(1,aa) + disploc(d,2,j,k)*hprime_xx(2,aa) + &
!    disploc(d,3,j,k)*hprime_xx(3,aa) + disploc(d,4,j,k)*hprime_xx(4,aa) + &
!    disploc(d,5,j,k)*hprime_xx(5,aa)
!
!!---
!  temp2(d,i,aa,k) = &
!    disploc(d,i,1,k)*hprime_yy(1,aa) + disploc(d,i,2,k)*hprime_yy(2,aa) + &
!    disploc(d,i,3,k)*hprime_yy(3,aa) + disploc(d,i,4,k)*hprime_yy(4,aa) + &
!    disploc(d,i,5,k)*hprime_yy(5,aa)
!
!!---
!  temp3(ijd,1,1,aa) = &
!    disploc(ijd,1,1,1)*hprime_zz(1,aa) + disploc(ijd,1,1,2)*hprime_zz(2,aa) + &
!    disploc(ijd,1,1,3)*hprime_zz(3,aa) + disploc(ijd,1,1,4)*hprime_zz(4,aa) + &
!    disploc(ijd,1,1,5)*hprime_zz(5,aa)
!
!
!!--------------------- second vector below ----------------------------
!
!!---
!  temp1inline(d,aa,j,k) = &
!    temp1(d,1,j,k)*hprimewgll_xx(aa,1) + temp1(d,2,j,k)*hprimewgll_xx(aa,2) + &
!    temp1(d,3,j,k)*hprimewgll_xx(aa,3) + temp1(d,4,j,k)*hprimewgll_xx(aa,4) + &
!    temp1(d,5,j,k)*hprimewgll_xx(aa,5)
!
!!---
!  temp2inline(d,i,aa,k) = &
!    temp2(d,i,1,k)*hprimewgll_yy(aa,1) + temp2(d,i,2,k)*hprimewgll_yy(aa,2) + &
!    temp2(d,i,3,k)*hprimewgll_yy(aa,3) + temp2(d,i,4,k)*hprimewgll_yy(aa,4) + &
!    temp2(d,i,5,k)*hprimewgll_yy(aa,5)
!
!!---
!  temp3inline(ijd,1,1,aa) = &
!   temp3(ijd,1,1,1)*hprimewgll_zz(aa,1) + temp3(ijd,1,1,2)*hprimewgll_zz(aa,2) + &
!   temp3(ijd,1,1,3)*hprimewgll_zz(aa,3) + temp3(ijd,1,1,4)*hprimewgll_zz(aa,4) + &
!   temp3(ijd,1,1,5)*hprimewgll_zz(aa,5)
!
!
!!--------------------- end of file ----------------------------
!
