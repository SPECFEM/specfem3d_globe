
  subroutine ylm(XLAT,XLON,LMAX,Y,WK1,WK2,WK3)

  implicit none

  complex TEMP,FAC,DFAC

  real(kind=4) WK1(1),WK2(1),WK3(1),Y(1),XLAT,XLON

  integer :: LMAX

!
!     WK1,WK2,WK3 SHOULD BE DIMENSIONED AT LEAST (LMAX+1)*4
!
  real(kind=4), parameter :: RADIAN = 57.2957795

  integer :: IM,IL1,IND,LM1,L

  real(kind=4) :: THETA,PHI

  THETA=(90.-XLAT)/RADIAN
  PHI=XLON/RADIAN

  IND=0
  LM1=LMAX+1

  DO IL1=1,LM1

  L=IL1-1
  CALL legndr(THETA,L,L,WK1,WK2,WK3)

  FAC=(1.,0.)
  DFAC=CEXP(CMPLX(0.,PHI))

  do IM=1,IL1
    TEMP=FAC*CMPLX(WK1(IM),0.)
    IND=IND+1
    Y(IND)=REAL(TEMP)
    if (IM == 1) goto 20
    IND=IND+1
    Y(IND)=AIMAG(TEMP)
 20 FAC=FAC*DFAC
  enddo

  enddo

  end subroutine ylm

!------------------------------------

  subroutine legndr(THETA,L,M,X,XP,XCOSEC)

  implicit none

  real(kind=4) :: X(2),XP(2),XCOSEC(2)

  double precision :: SMALL,SUM,COMPAR,CT,ST,FCT,COT,X1,X2,X3,F1,F2,XM,TH

  double precision, parameter :: FPI = 12.56637062D0

  integer :: i,M,MP1,k,l,LP1

  real(kind=4) :: THETA,DSFL3,COSEC,SFL3

!!!!!! illegal statement, removed by Dimitri Komatitsch   DFLOAT(I)=FLOAT(I)

  SUM=0.D0
  LP1=L+1
  TH=THETA
  CT=DCOS(TH)
  ST=DSIN(TH)
  MP1=M+1
  FCT=DSQRT(dble(2*L+1)/FPI)
  SFL3=SQRT(FLOAT(L*(L+1)))
  COMPAR=dble(2*L+1)/FPI
  DSFL3=SFL3
  SMALL=1.D-16*COMPAR

  do I=1,MP1
    X(I)=0.
    XCOSEC(I)=0.
    XP(I)=0.
  enddo

  if (L > 1 .and. ABS(THETA) > 1.E-5) goto 3
  X(1)=FCT
  if (L == 0) return
  X(1)=CT*FCT
  X(2)=-ST*FCT/DSFL3
  XP(1)=-ST*FCT
  XP(2)=-.5D0*CT*FCT*DSFL3
  if (ABS(THETA) < 1.E-5) XCOSEC(2)=XP(2)
  if (ABS(THETA) >= 1.E-5) XCOSEC(2)=X(2)/ST
  return

 3 X1=1.D0
  X2=CT

  do I=2,L
    X3=(dble(2*I-1)*CT*X2-dble(I-1)*X1)/dble(I)
    X1=X2
    X2=X3
  enddo

  COT=CT/ST
  COSEC=1./ST
  X3=X2*FCT
  X2=dble(L)*(X1-CT*X2)*FCT/ST
  X(1)=X3
  X(2)=X2
  SUM=X3*X3
  XP(1)=-X2
  XP(2)=dble(L*(L+1))*X3-COT*X2
  X(2)=-X(2)/SFL3
  XCOSEC(2)=X(2)*COSEC
  XP(2)=-XP(2)/SFL3
  SUM=SUM+2.D0*X(2)*X(2)
  if (SUM-COMPAR > SMALL) return
  X1=X3
  X2=-X2/DSQRT(dble(L*(L+1)))

  do I=3,MP1
    K=I-1
    F1=DSQRT(dble((L+I-1)*(L-I+2)))
    F2=DSQRT(dble((L+I-2)*(L-I+3)))
    XM=K
    X3=-(2.D0*COT*(XM-1.D0)*X2+F2*X1)/F1
    SUM=SUM+2.D0*X3*X3
    if (SUM-COMPAR > SMALL .and. I /= LP1) return
    X(I)=X3
    XCOSEC(I)=X(I)*COSEC
    X1=X2
    XP(I)=-(F1*X2+XM*COT*X3)
    X2=X3
  enddo

  end subroutine legndr

