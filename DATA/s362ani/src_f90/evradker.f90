
  subroutine evradker(depth,string,nker,vercof,dvercof,ierror)

  implicit none

  integer :: nker,ierror

  real(kind=4) :: chebyshev(100)
  real(kind=4) :: chebyshev2(100)
  real(kind=4) :: vercof(nker)
  real(kind=4) :: dvercof(nker)
  real(kind=4) :: splpts(100)

  character(len=80) string

  logical upper,upper_650
  logical lower,lower_650

  real(kind=4), parameter :: r0=6371.
  real(kind=4), parameter :: rmoho=6371.0-24.4
  real(kind=4), parameter :: r670=6371.-670.
  real(kind=4), parameter :: r650=6371.-650.
  real(kind=4), parameter :: rcmb=3480.0

  integer :: i,nspl,nskip,nlower,nupper,iker,lstr

  real(kind=4) :: u,u2,ddep,radius2,radius,depth

  ierror=0
  lstr=len_trim(string)

  radius=r0-depth
  ddep=0.1
  radius2=r0-depth+ddep
  upper=.false.
  lower=.false.
  if (radius > rcmb .and. radius < r670) then
  lower=.true.
  else if (radius >= r670 .and. radius < rmoho) then
  upper=.true.
  endif
  upper_650=.false.
  lower_650=.false.
  if (radius > rcmb .and. radius < r650) then
  lower_650=.true.
  else if (radius >= r650 .and. radius < rmoho) then
  upper_650=.true.
  endif
  do iker=1,nker
  vercof(iker)=0.
  dvercof(iker)=0.
  enddo

  if (string(1:16) == 'WDC+SPC_U4L8CHEB') then
  nupper=5
  nlower=9
  nskip=2
  if (upper) then
    u=(radius+radius-rmoho-r670)/(rmoho-r670)
    u2=(radius2+radius2-rmoho-r670)/(rmoho-r670)
!          write(6,"('upper mantle:',2f10.3)") u,u2
    call chebyfun(u,13,chebyshev)
    do i=1+nskip,nskip+nupper
      vercof(i)=chebyshev(i-nskip)
    enddo
    call chebyfun(u2,13,chebyshev2)
    do i=1+nskip,nskip+nupper
      dvercof(i)=(chebyshev2(i-nskip)-chebyshev(i-nskip))/ddep
    enddo
  else if (lower) then
    u=(radius+radius-r670-rcmb)/(r670-rcmb)
    u2=(radius2+radius2-r670-rcmb)/(r670-rcmb)
!          write(6,"('lower mantle:',2f10.3)") u,u2
    call chebyfun(u,13,chebyshev)
    do i=1+nskip+nupper,nskip+nupper+nlower
      vercof(i)=chebyshev(i-nskip-nupper)
    enddo
    call chebyfun(u2,13,chebyshev2)
    do i=1+nskip+nupper,nskip+nupper+nlower
      dvercof(i)=(chebyshev2(i-nskip-nupper)- &
                    chebyshev(i-nskip-nupper))/ddep
    enddo
  endif
  else if (string(1:13) == 'WDC+SHSVWM20A') then
  nspl=20
  splpts(1)=0.
  splpts(2)=50.
  splpts(3)=100.
  splpts(4)=150.
  splpts(5)=200.
  splpts(6)=250.
  splpts(7)=300.
  splpts(8)=400.
  splpts(9)=500.
  splpts(10)=600.
  splpts(11)=700.
  splpts(12)=850.
  splpts(13)=1050.
  splpts(14)=1300.
  splpts(15)=1600.
  splpts(16)=1900.
  splpts(17)=2200.
  splpts(18)=2500.
  splpts(19)=2700.
  splpts(20)=2891.
  call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
  do i=22,27
    vercof(i)=vercof(i-20)
    dvercof(i)=dvercof(i-20)
  enddo
  vercof(1)=1.
  else if (string(1:16) == 'WDC+XBS_362_U6L8') then
  if (upper) then
   nspl=6
   splpts(1)=24.4
   splpts(2)=100.
   splpts(3)=225.
   splpts(4)=350.
   splpts(5)=500.
   splpts(6)=670.
   call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
  else if (lower) then
 nspl=8
   splpts(1)=670.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(8),dvercof(8))
  endif
  vercof(1)=1.
!        vercof(16)=1.
!        vercof(17)=1.
!      else if (string(1:21)=='WDC+ANI_362_U6L8_TOPO') then
!        if (upper) then
!         nspl=6
!         splpts(1)=24.4
!         splpts(2)=100.
!         splpts(3)=225.
!         splpts(4)=350.
!         splpts(5)=500.
!         splpts(6)=670.
!         call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
!         do i=16,21
!          vercof(i)=vercof(i-14)
!          dvercof(i)=dvercof(i-14)
!         enddo
!     else if (lower) then
!      nspl=8
!         splpts(1)=670.
!         splpts(2)=820.
!         splpts(3)=1320.
!         splpts(4)=1820.
!         splpts(5)=2320.
!         splpts(6)=2550.
!         splpts(7)=2791.
!         splpts(8)=2891.
!         call vbspl(depth,nspl,splpts,vercof(8),dvercof(8))
!     endif
!        vercof(1)=1.
!        vercof(22)=1.
!        vercof(23)=1.
!        vercof(24)=1.
!        vercof(25)=1.
  else if ( &
       (string(1:lstr) == 'WDC+ANI_362_U6L8' .and. lstr == 16) &
       .or. &
           (string(1:lstr) == 'WDC+ANI_362_U6L8_TOPO' .and. lstr == 21) &
       ) then
  if (upper) then
   nspl=6
   splpts(1)=24.4
   splpts(2)=100.
   splpts(3)=225.
   splpts(4)=350.
   splpts(5)=500.
   splpts(6)=670.
   call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
   do i=16,21
    vercof(i)=vercof(i-14)
    dvercof(i)=dvercof(i-14)
   enddo
  else if (lower) then
 nspl=8
   splpts(1)=670.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(8),dvercof(8))
  endif
  vercof(1)=1.
  vercof(22)=1.
  vercof(23)=1.
  else if (string(1:lstr) == 'WDC+WM_362_U6L8' .and. lstr == 15) then
  if (upper) then
   nspl=6
   splpts(1)=24.4
   splpts(2)=100.
   splpts(3)=225.
   splpts(4)=350.
   splpts(5)=500.
   splpts(6)=670.
   call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
   do i=16,21
    vercof(i)=vercof(i-14)
    dvercof(i)=dvercof(i-14)
   enddo
  else if (lower) then
 nspl=8
   splpts(1)=670.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(8),dvercof(8))
   do i=22,29
    vercof(i)=vercof(i-14)
    dvercof(i)=dvercof(i-14)
   enddo
  endif
  vercof(1)=1.
  vercof(30)=1.
  vercof(31)=1.
  vercof(32)=1.
  else if ( &
     (string(1:lstr) == 'WDC+ANI_362_U6L8_650' .and. lstr == 20) &
     .or. &
         (string(1:lstr) == 'WDC+ANI_362_U6L8_TOPO_650' .and. lstr == 25) &
     ) then
  if (upper_650) then
   nspl=6
   splpts(1)=24.4
   splpts(2)=100.
   splpts(3)=225.
   splpts(4)=350.
   splpts(5)=500.
   splpts(6)=650.
   call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
   do i=16,21
    vercof(i)=vercof(i-14)
    dvercof(i)=dvercof(i-14)
   enddo
  else if (lower_650) then
 nspl=8
   splpts(1)=650.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(8),dvercof(8))
  endif
  vercof(1)=1.
  vercof(22)=1.
  vercof(23)=1.
  else if (string(1:lstr) == 'WDC+WM_362_U6L8_650' &
       .and.lstr == 19) then
  if (upper_650) then
   nspl=6
   splpts(1)=24.4
   splpts(2)=100.
   splpts(3)=225.
   splpts(4)=350.
   splpts(5)=500.
   splpts(6)=650.
   call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
   do i=16,21
    vercof(i)=vercof(i-14)
    dvercof(i)=dvercof(i-14)
   enddo
  else if (lower_650) then
 nspl=8
   splpts(1)=650.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(8),dvercof(8))
   do i=22,29
    vercof(i)=vercof(i-14)
    dvercof(i)=dvercof(i-14)
   enddo
  endif
  vercof(1)=1.
  vercof(30)=1.
  vercof(31)=1.
  vercof(32)=1.
  else if (string(1:lstr) == 'WDC+U8L8_650' .and. lstr == 12) then
  if (upper_650) then
   nspl=8
   splpts(1)=24.4
   splpts(2)=75.
   splpts(3)=150.
   splpts(4)=225.
   splpts(5)=300.
   splpts(6)=410.
   splpts(7)=530.
   splpts(8)=650.
   call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
   do i=18,25
    vercof(i)=vercof(i-16)
    dvercof(i)=dvercof(i-16)
   enddo
  else if (lower_650) then
 nspl=8
   splpts(1)=650.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(10),dvercof(10))
   do i=26,33
    vercof(i)=vercof(i-16)
    dvercof(i)=dvercof(i-16)
   enddo
  endif
  vercof(1)=1.
  vercof(34)=1.
  vercof(35)=1.
  vercof(36)=1.
  else if (string(1:lstr) == 'WDC+U8L8_670' .and. lstr == 12) then
  if (upper) then
   nspl=8
   splpts(1)=24.4
   splpts(2)=75.
   splpts(3)=150.
   splpts(4)=225.
   splpts(5)=300.
   splpts(6)=410.
   splpts(7)=530.
   splpts(8)=670.
   call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
   do i=18,25
    vercof(i)=vercof(i-16)
    dvercof(i)=dvercof(i-16)
   enddo
  else if (lower) then
 nspl=8
   splpts(1)=670.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(10),dvercof(10))
   do i=26,33
    vercof(i)=vercof(i-16)
    dvercof(i)=dvercof(i-16)
   enddo
  endif
  vercof(1)=1.
  vercof(34)=1.
  vercof(35)=1.
  vercof(36)=1.
  else if ( &
      (string(1:lstr) == 'WDC+U8L8_I1D_650' .and. lstr == 16) &
      .or. &
      (string(1:lstr) == 'WDC+U8L8_I3D_650' .and. lstr == 16) &
      ) then
  if (upper_650) then
   nspl=8
   splpts(1)=24.4
   splpts(2)=75.
   splpts(3)=150.
   splpts(4)=225.
   splpts(5)=300.
   splpts(6)=410.
   splpts(7)=530.
   splpts(8)=650.
   call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
   do i=18,25
    vercof(i)=vercof(i-16)
    dvercof(i)=dvercof(i-16)
   enddo
   do i=37,40
    vercof(i)=vercof(i-35)
    dvercof(i)=dvercof(i-35)
   enddo
   do i=41,44
    vercof(i)=vercof(i-39)
    dvercof(i)=dvercof(i-39)
   enddo
   do i=45,48
    vercof(i)=vercof(i-43)
    dvercof(i)=dvercof(i-43)
   enddo
   do i=49,52
    vercof(i)=vercof(i-47)
    dvercof(i)=dvercof(i-47)
   enddo
  else if (lower_650) then
 nspl=8
   splpts(1)=650.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(10),dvercof(10))
   do i=26,33
    vercof(i)=vercof(i-16)
    dvercof(i)=dvercof(i-16)
   enddo
  endif
  vercof(1)=1.
  vercof(34)=1.
  vercof(35)=1.
  vercof(36)=1.
  else if ((string(1:lstr) == 'WDC+I1D_650' .and. lstr == 11) .or. &
          (string(1:lstr) == 'WDC+I3D_650' .and. lstr == 11)) then
  if (upper_650) then
   nspl=8
   splpts(1)=24.4
   splpts(2)=75.
   splpts(3)=150.
   splpts(4)=225.
   splpts(5)=300.
   splpts(6)=410.
   splpts(7)=530.
   splpts(8)=650.
   call vbspl(depth,nspl,splpts,vercof(2),dvercof(2))
   do i=18,25
    vercof(i)=vercof(i-16)
    dvercof(i)=dvercof(i-16)
   enddo
   do i=37,44
    vercof(i)=vercof(i-35)
    dvercof(i)=dvercof(i-35)
   enddo
   do i=53,60
    vercof(i)=vercof(i-51)
    dvercof(i)=dvercof(i-51)
   enddo
   do i=69,76
    vercof(i)=vercof(i-67)
    dvercof(i)=dvercof(i-67)
   enddo
   do i=85,92
    vercof(i)=vercof(i-83)
    dvercof(i)=dvercof(i-83)
   enddo
  else if (lower_650) then
 nspl=8
   splpts(1)=650.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(10),dvercof(10))
   do i=26,33
    vercof(i)=vercof(i-16)
    dvercof(i)=dvercof(i-16)
   enddo
   do i=45,52
    vercof(i)=vercof(i-35)
    dvercof(i)=dvercof(i-35)
   enddo
   do i=61,68
    vercof(i)=vercof(i-51)
    dvercof(i)=dvercof(i-51)
   enddo
   do i=77,84
    vercof(i)=vercof(i-67)
    dvercof(i)=dvercof(i-67)
   enddo
   do i=93,100
    vercof(i)=vercof(i-83)
    dvercof(i)=dvercof(i-83)
   enddo
  endif
  vercof(1)=1.
  vercof(34)=1.
  vercof(35)=1.
  vercof(36)=1.
  else if (string(1:lstr) == 'V16A4_V7A4' .and. lstr == 10) then
  if (upper_650) then
   nspl=8
   splpts(1)=24.4
   splpts(2)=75.
   splpts(3)=150.
   splpts(4)=225.
   splpts(5)=300.
   splpts(6)=410.
   splpts(7)=530.
   splpts(8)=650.
   call vbspl(depth,nspl,splpts,vercof(1),dvercof(1))
   do i=17,20
    vercof(i)=vercof(i-16)
    dvercof(i)=dvercof(i-16)
   enddo
   do i=23,29
    vercof(i)=vercof(i-22)
    dvercof(i)=dvercof(i-22)
   enddo
   do i=30,33
    vercof(i)=vercof(i-29)
    dvercof(i)=dvercof(i-29)
   enddo
  else if (lower_650) then
 nspl=8
   splpts(1)=650.
   splpts(2)=820.
   splpts(3)=1320.
   splpts(4)=1820.
   splpts(5)=2320.
   splpts(6)=2550.
   splpts(7)=2791.
   splpts(8)=2891.
   call vbspl(depth,nspl,splpts,vercof(9),dvercof(9))
  endif
  vercof(21)=1.
  vercof(22)=1.
  else
  write(6,"('problem 4')")
  write(6,"(a)")string(1:len_trim(string))
  stop
  endif

  end subroutine evradker

! ---

  subroutine chebyfun(u,kmax,f)

  implicit none

  integer :: kmax

  real(kind=4) :: chebycoeff(0:13),f(0:kmax),u

  integer :: k

  real(kind=4) :: twou

  data chebycoeff / &
   0.70710678118655,1.2247448713916,1.0350983390135,1.0145993123918, &
   1.00803225754840,1.0050890913907,1.0035149493262,1.0025740068320, &
   1.00196657023780,1.0015515913133,1.0012554932754,1.0010368069141, &
   1.00087070107920,1.0007415648034 /

  if (kmax > 13) then
   write(*,"(' kmax exceeds the limit in chebyfun')")
   stop
  endif

  f(0)=1.0
  f(1)=u
  twou=2.0*u

  do k=2,kmax
   f(k) = twou*f(k-1)-f(k-2)
  enddo

  do k=0,kmax
   f(k)=f(k)*chebycoeff(k)
  enddo

  end subroutine chebyfun

