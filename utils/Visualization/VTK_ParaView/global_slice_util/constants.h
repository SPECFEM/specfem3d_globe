  integer, parameter :: SIZE_REAL = 4
  integer, parameter :: SIZE_DOUBLE = 8

! uncomment this to run in single precision
 integer, parameter :: CUSTOM_REAL = SIZE_REAL
! uncomment this to run in double precision (increases memory size by 2)
!  integer, parameter :: CUSTOM_REAL = SIZE_DOUBLE

 real(kind=CUSTOM_REAL), parameter :: pi = 3.1415926
 real(kind=CUSTOM_REAL), parameter :: EPS = 0.01_CUSTOM_REAL

 integer, parameter :: NPTS = 100000
