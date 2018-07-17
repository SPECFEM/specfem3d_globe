subroutine diff_lagrange(N,DL,rho,xi)
!
!calculates the derivatives of the Lagrange Polynomials at the GLL quadrature
!points and stores them in Matrix DL (coloumn = d (L)/d (xi) = dL, row = dL at
!point xi

implicit none

double precision xi(0:N),LL(0:N),VN(0:N),rho(0:N)
double precision DL(0:N,0:N)
integer N,i,j,k



!Calculate the Gauss-Lobatto-Legendre Quadrature points and the values of the
!Legendre polynomials VN at the GLL nodes

call main_pol(xi,N,VN)

!calculate the weights rho(i) at every GLL node for the GLL integration rule

rho = 2 / ( float((N+1)*N) * VN**2)

!calculating the derivatives of the Polynomials at the  Gauss-Lobatto-Legendre Quadrature points

do i=0,N

    LL=0
    LL(i)=1
    call DELEGL(N,xi,VN,LL,DL(:,i))

enddo

end subroutine diff_lagrange

!*****************************************************************************
subroutine main_pol(xi,N,VN)

implicit none

DOUBLE PRECISION X,xi(0:N),Y,DY,D2Y,ET(0:1000),VN(0:1000)
INTEGER j,k,N

DO j=0,1000
X=-1+REAL(j)/1000.*2
!CALL VALEPO(N,X,Y,DY,D2Y)
CALL ZELEGL(N,ET,VN)
!WRITE(5,*)X,Y,DY,D2Y
enddo
DO j=0,N

  xi(j)=ET(j)

!write(*,*)ET(j),VN(j)
enddo
!CLOSE(5)
!CLOSE(6)

end subroutine main_pol


!****************************************************************************
      subroutine VALEPO(N,X,Y,DY,D2Y)

!*   COMPUTES THE VALUE OF THE LEGENDRE POLYNOMIAL OF DEGREE N
!*   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT
!*   N  = DEGREE OF THE POLYNOMIAL
!*   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED
!*   Y  = VALUE OF THE POLYNOMIAL IN X
!*   DY = VALUE OF THE FIRST DERIVATIVE IN X
!*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X

!      implicit DOUBLE PRECISION (A-H,O-Z)

      implicit none
      DOUBLE PRECISION Y,DY,D2Y,X,YP,DYP,D2YP,C1,C2,C3,C4,YM,DYM,D2YM
      INTEGER N,I

         Y   = 1.
         DY  = 0.
         D2Y = 0.

      if (N == 0) return

         Y   = X
         DY  = 1.
         D2Y = 0.

      if (N == 1) return

         YP   = 1.
         DYP  = 0.
         D2YP = 0.

      DO I=2,N
         C1 = DFLOAT(I)
         C2 = 2.*C1-1.
         C4 = C1-1.
         YM = Y
         Y  = (C2*X*Y-C4*YP)/C1
         YP = YM
         DYM  = DY
         DY   = (C2*X*DY-C4*DYP+C2*YP)/C1
         DYP  = DYM
         D2YM = D2Y
         D2Y  = (C2*X*D2Y-C4*D2YP+2.D0*C2*DYP)/C1
         D2YP = D2YM
       enddo

      return

end subroutine VALEPO


!****************************************************************************

      subroutine ZELEGL(N,ET,VN)

!*   COMPUTES THE NODES RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
!*   N  = ORDER OF THE FORMULA
!*   ET = VECTOR OF THE NODES, ET(I), I=0,N
!*   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N

!      implicit DOUBLE PRECISION (A-H,O-Z)
      implicit none

      DOUBLE PRECISION SN,X,Y,DY,D2Y,PI,C,ETX
      DOUBLE PRECISION ET(0:1000),VN(0:1000)
      INTEGER N,N2,I,IT

!     DIMENSION ET(0:*), VN(0:*)

      if (N == 0) return

         N2 = (N-1)/2
         SN = DFLOAT(2*N-4*N2-3)
         ET(0) = -1.
         ET(N) = 1.
         VN(0) = SN
         VN(N) = 1.
      if (N == 1) return

         ET(N2+1) = 0.
         X = 0.
      CALL VALEPO(N,X,Y,DY,D2Y)
         VN(N2+1) = Y
      if (N == 2) return

         PI = 3.14159265358979323846
         C  = PI/DFLOAT(N)
      DO I=1,N2
        ETX = DCOS(C*DFLOAT(I))
        DO IT=1,8
         CALL VALEPO(N,ETX,Y,DY,D2Y)
         ETX = ETX-DY/D2Y
        enddo
        ET(I) = -ETX
        ET(N-I) = ETX
        VN(I) = Y*SN
        VN(N-I) = Y
      enddo

      return
      end subroutine ZELEGL


!***********************************************************************

      subroutine DELEGL(N,ET,VN,QN,DQN)

!************************************************************************
!*  COMPUTES THE DERIVATIVE OF A POLYNOMIAL AT THE LEGENDRE GAUSS-LOBATTO
!*  NODES FROM THE VALUES OF THE POLYNOMIAL ATTAINED AT THE SAME POINTS
!*   N   = THE DEGREE OF THE POLYNOMIAL
!*   ET  = VECTOR OF THE NODES, ET(I), I=0,N
!*   VN  = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
!*   QN  = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N
!*   DQN = DERIVATIVES OF THE POLYNOMIAL AT THE NODES, DQZ(I), I=0,N
!************************************************************************
      implicit DOUBLE PRECISION (A-H,O-Z)
!      DIMENSION ET(0:*), VN(0:*), QN(0:*), DQN(0:*)
double precision ET(0:*), VN(0:*), QN(0:*), DQN(0:*)
          DQN(0) = 0.D0
      if (N == 0) return

      DO 1 I=0,N
          SU = 0.D0
          VI = VN(I)
          EI = ET(I)
      DO 2 J=0,N
      if (I == J) goto 2
          VJ = VN(J)
          EJ = ET(J)
          SU = SU+QN(J)/(VJ*(EI-EJ))
2     continue
          DQN(I) = VI*SU
1     continue

          DN = DFLOAT(N)
          C  = .25D0*DN*(DN+1.D0)
          DQN(0) = DQN(0)-C*QN(0)
          DQN(N) = DQN(N)+C*QN(N)

      return
      end subroutine DELEGL


