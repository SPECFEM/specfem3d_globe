!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

!--------------------------------------------------------------------------------------------------
!
! Berkeley crust implemented as in csem
!
!--------------------------------------------------------------------------------------------------

module model_crust_berkeley_par

  implicit none

  integer :: NBT,NBP
  real :: drin
  real, parameter :: drfiltre = 2.e0
  !double precision, parameter :: dr_ = 2.d0

  real, dimension(:,:), allocatable :: crust_array
  real, dimension(:,:), allocatable :: moho_start

  ! moho depth from 1D reference model (in km)
  double precision :: moho1D_depth = -1.0d0

end module model_crust_berkeley_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_berkeley_crust_broadcast()

! standard routine to setup model

  use model_crust_berkeley_par
  use constants, only: A3d_folder,myrank,IMAIN

  implicit none

  character(len=100), parameter :: file_crust = trim(A3d_folder) // 'crust2cru2av_2x2.dat'
  character(len=100), parameter :: file_moho = trim(A3d_folder) // 'crust2moho_2x2.dat'
  integer :: ier

  ! user output
  if (myrank == 0) then
    ! user output
    write(IMAIN,*)
    write(IMAIN,*) 'incorporating crustal model: Berkeley'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! moho data
  open(52,file=trim(file_moho),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(file_moho)
    stop 'Error opening file crust2moho_2x2.dat'
  endif
  call read_crustmoho_filtre(52)
  close(52)

  !drfiltre = dr_

  ! crustal data
  open(52,file=trim(file_crust),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(file_crust)
    stop 'Error opening file crust2cru2av_2x2.dat'
  endif
  call read_crust_smooth_variable(52)
  close(52)

  ! gets moho depth from 1D model
  if (myrank == 0) then
    call determine_1dberkeley_moho_depth(moho1D_depth)
    !debug
    !print *,'debug: [model_berkeley_crust_broadcast] moho1d_depth = ',moho1D_depth
  endif
  ! broadcasts to all other processes
  call bcast_all_singledp(moho1D_depth)

  end subroutine model_berkeley_crust_broadcast

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_berkeley_crust(x,theta,phi,vp,vs,rho,moho,found_crust)

! returns isotropic crustal velocities

  use model_crust_berkeley_par
  use constants, only: PI,GRAV,EARTH_RHOAV,EARTH_R,EARTH_R_KM

  implicit none

  double precision,intent(in) :: x,theta,phi
  double precision,intent(out) :: vp,vs,rho,moho
  logical,intent(out) :: found_crust

  ! local parameters
  double precision :: vsv,vsh
  double precision :: depth
  double precision :: moho_depth, scaleval
  !double precision, parameter :: deg2rad = PI/180.d0, rad2deg = 180.d0/PI

  ! initializes
  found_crust = .true.

  !theta = (90.d0 - lat) * deg2rad ! assumed lat range: [-90,90]
  !phi = lon * deg2rad ! assumed lon range: [-180,180]

  depth = (1.d0 - x) * EARTH_R_KM

  call get_crust_val_csem(theta,phi,depth,rho,vp,vsv,vsh,moho_depth)

  !TODO: daniel - check flag usage
  ! sets flag if point within crust
  !if (depth > moho_depth) found_crust = .false.

  !
  ! get equivalent isotropic vs
  !
  vs = dsqrt((2.d0 * vsv*vsv + vsh*vsh)/3.d0)

  !
  ! scale values for specfem
  !
  scaleval = dsqrt(PI*GRAV*EARTH_RHOAV)

  vp = vp / (EARTH_R * scaleval)
  vs = vs / (EARTH_R * scaleval)
  rho = rho / EARTH_RHOAV

  moho = moho_depth / EARTH_R_KM

  end subroutine model_berkeley_crust

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_berkeley_crust_aniso(x,theta,phi,vpv,vph,vsv,vsh,eta_aniso,rho,moho,found_crust)

! returns anisotropic crustal velocities

  use model_crust_berkeley_par
  use constants, only: PI,GRAV,EARTH_RHOAV,EARTH_R,EARTH_R_KM

  implicit none

  double precision,intent(in) :: x,theta,phi
  double precision,intent(out) :: vpv,vph,vsv,vsh,eta_aniso,rho,moho
  logical,intent(out) :: found_crust

  ! local parameters
  double precision :: vp
  double precision :: depth
  double precision :: moho_depth, scaleval
  !double precision, parameter :: deg2rad = PI/180.d0, rad2deg = 180.d0/PI

  ! initializes
  eta_aniso = 1.d0
  found_crust = .true.

  !theta = (90.d0 - lat) * deg2rad ! assumed lat range: [-90,90]
  !phi = lon * deg2rad ! assumed lon range: [-180,180]

  depth = (1.d0 - x) * EARTH_R_KM

  call get_crust_val_csem(theta,phi,depth,rho,vp,vsv,vsh,moho_depth)

  !TODO: daniel - check flag usage
  ! sets flag if point within crust
  !if (depth > moho_depth) found_crust = .false.

  !
  ! scale values for specfem
  !
  scaleval = dsqrt(PI*GRAV*EARTH_RHOAV)

  vph = vp / (EARTH_R * scaleval)
  vpv = vp / (EARTH_R * scaleval)
  vsh = vsh / (EARTH_R * scaleval)
  vsv = vsv / (EARTH_R * scaleval)
  rho = rho / EARTH_RHOAV

  moho = moho_depth / EARTH_R_KM

  end subroutine model_berkeley_crust_aniso

!
!--------------------------------------------------------------------------------------------------
!

  subroutine get_crust_val_csem(theta,phi,z,rho,vp,vsv,vsh,moho_depth)

  use model_crust_berkeley_par

  implicit none

  double precision,intent(in) :: theta,phi,z
  double precision,intent(out) :: rho,vp,vsv,vsh,moho_depth

  ! local parameters
  double precision :: xi(5),wi(5)
  double precision :: rho_cr(5),vp_cr(5),vsv_cr(5),vsh_cr(5)
  double precision :: x
  integer :: i

  real,external :: moho_filtre
  double precision,external :: lagrange

  ! get moho depth
  moho_depth = moho1D_depth - dble(moho_filtre(theta,phi))

  !debug
  !print *,"debug: [get_crust_val_csem] Moho depth:",moho_depth, &
  !        "moho1D_depth:",moho1D_depth,"moho_filtre:",moho_filtre(theta,phi)

  !
  ! horizontal interpolation for all registered depths
  !
  call crust_bilinear_variable(theta,phi,1,rho_cr(1),vp_cr(1),vsv_cr(1),vsh_cr(1))
  call crust_bilinear_variable(theta,phi,2,rho_cr(2),vp_cr(2),vsv_cr(2),vsh_cr(2))
  call crust_bilinear_variable(theta,phi,3,rho_cr(3),vp_cr(3),vsv_cr(3),vsh_cr(3))
  call crust_bilinear_variable(theta,phi,4,rho_cr(4),vp_cr(4),vsv_cr(4),vsh_cr(4))
  call crust_bilinear_variable(theta,phi,5,rho_cr(5),vp_cr(5),vsv_cr(5),vsh_cr(5))

  !
  ! get GLL nodes position
  !
  call LGL_NODES(4,1.d-6,xi,wi)
  xi(:) = moho_depth - (xi(:) + 1.d0)/2.d0 * moho_depth
  x = z
  if (x > maxval(xi)) x = maxval(xi)
  if (x < minval(xi)) x = minval(xi)

  !
  ! init values
  !
  vp  = 0.d0
  vsv = 0.d0
  vsh = 0.d0
  rho = 0.d0

  !
  ! depth interpolation
  !
  do i = 1,5
    rho = rho + lagrange(i,xi,5,x) * rho_cr(i)
    vp  = vp  + lagrange(i,xi,5,x) * vp_cr(i)
    vsv = vsv + lagrange(i,xi,5,x) * vsv_cr(i)
    vsh = vsh + lagrange(i,xi,5,x) * vsh_cr(i)
  enddo

  rho = rho * 1000.d0
  vp  = vp  * 1000.d0
  vsv = vsv * 1000.d0
  vsh = vsh * 1000.d0

  end subroutine get_crust_val_csem

!
!--------------------------------------------------------------------------------------------------
!

  double precision function lagrange(ind,xi,nxi,x)

  implicit none

  integer,intent(in) :: ind,nxi
  double precision,intent(in) :: xi(nxi),x

  integer :: i

  lagrange = 1.d0

  do i = 1,nxi
    if (i /= ind) lagrange = lagrange * ( (x-xi(i))/(xi(ind)-xi(i)) )
  enddo

  end function lagrange

!
!--------------------------------------------------------------------------------------------------
!

  subroutine LGL_NODES(N,EPS,XI,WI)

  implicit none

  INTEGER :: N ! POLYNOMIAL ORDER
  DOUBLE PRECISION :: EPS ! DESIRED ERROR

  DOUBLE PRECISION :: XI(N+1) ! GAUSS-LOBATTO-LEGENDRE INTERPOLATION POINTS
  DOUBLE PRECISION :: WI(N+1) ! GAUSS-LOBATTO-LEGENDRE QUADRATURE WEIGHT

  DOUBLE PRECISION :: P(N+1,N+1),XOLD(N+1)
  DOUBLE PRECISION :: ERMAX
  INTEGER :: I,J,N1
  DOUBLE PRECISION, PARAMETER :: PI = 3.141592653589793d0  ! PI = 4.D0 * DATAN(1.D0)

  N1 = N+1
  do I = 1,N1
    XI(I) = -COS(PI * real(I-1)/real(N1-1))
  enddo

  ERMAX = 2.d0 * EPS

  DO WHILE(ERMAX > EPS)
    do I = 1,N1
      XOLD(I) = XI(I)
      P(I,1) = 1.d0
      P(I,2) = XI(I)
    enddo
    do J = 2,N
      do I = 1,N1
        P(I,J+1) = ((2*J-1) * XI(I) * P(I,J) - (J-1) * P(I,J-1)) / REAL(J)
      enddo
    enddo
    do I = 1,N1
      XI(I) = XOLD(I) - (XI(I) * P(I,N1) - P(I,N)) / (N1 * P(I,N1))
    enddo
    ERMAX = 0.d0
    do I = 1,N1
      ERMAX = MAX(ERMAX, abs(XI(I)-XOLD(I)))
    enddo
  enddo

  do I = 1,N1
    WI(I) = 2.d0 / (N * N1 * P(I,N1)**2)
  enddo

  end subroutine LGL_NODES

!
!--------------------------------------------------------------------------------------------------
!

  real function moho_filtre(theta,phi)

!theta phi en radians
!dr en degre
!reponse en metre

  use model_crust_berkeley_par

  implicit none

  double precision,intent(in) :: theta,phi
  ! local parameter
  real :: t,p
  real, parameter :: pi = 3.141592653589793, deg2rad = pi/180.0, rad2deg = 180.0/pi

  t = sngl(theta) / deg2rad
  p = sngl(phi) / deg2rad
  moho_filtre = gauss_filtre1(moho_start,t,p,drfiltre)

  contains

    !-----------------------------------------------------------------
    real function gauss_filtre1(tin,theta,phi,dr)
    !-----------------------------------------------------------------
    implicit none
    real,intent(in) :: theta,phi,dr
    real, dimension(:,:),intent(in) :: tin
    ! local parameters
    real :: thetar,phir,tmpnorm,inte
    real :: tmp
    integer :: i,ii,j,jj,LARG

    tmp = 0.0
    tmpnorm = 0.0
    LARG = 10
    do i = 1,LARG+1
       do j = 1,LARG+1
          call get_indexloc(phi,theta,i,j,dr,LARG,ii,jj,phir,thetar)
          inte = cos_cylindre(theta,phi,dr,thetar,phir) * (dr/real(LARG/2) * deg2rad)**2 * sin(thetar * deg2rad)
          tmp = tmp + tin(ii,jj) * inte
          tmpnorm = tmpnorm + inte
       enddo
    enddo
    gauss_filtre1 = tmp / tmpnorm

    end function gauss_filtre1

    !----------------------------------------------------------------------
    real function cos_cylindre(t0_,p0_,d0_,theta_,phi_)
    !----------------------------------------------------------------------
    implicit none
    real,intent(in) :: t0_,p0_
    real,intent(in) :: d0_,theta_,phi_
    ! local parameters
    real :: t0,p0,d0,theta,phi,d_ang
    real :: cosa

    t0 = t0_ * deg2rad
    p0 = p0_ * deg2rad
    theta = theta_ * deg2rad
    phi = phi_ * deg2rad
    d0 = d0_ * deg2rad
    !distance angulaire au centre du cylindre:
    cosa = cos(theta) * cos(t0) + sin(theta) * sin(t0) * cos(phi - p0)
    if (cosa >= 1.0) then
      d_ang = 0.0
    else if (cosa <= -1.0) then
      d_ang = 4.0 * atan(1.0)
    else
      d_ang = acos(cos(theta) * cos(t0) + sin(theta) * sin(t0) * cos(phi - p0))
    endif
    if (d_ang > d0) then
      cos_cylindre = 0.0
    else
      cos_cylindre = 0.5 * (1.0 + cos(PI * d_ang/d0))
    endif

    end function cos_cylindre

  end function moho_filtre

!
!--------------------------------------------------------------------------------------------------
!

  subroutine get_indexloc(phi,theta,i,j,dr,LARG,ii,jj,phir,thetar)

  use model_crust_berkeley_par

  implicit none
  real, intent(in) :: theta,phi,dr
  integer, intent(in) :: i,j,LARG
  real, intent(out) :: thetar,phir
  integer, intent(out) :: ii,jj

  ! local parameters
  double precision :: t,p
  double precision, parameter :: eps = 1.d-8

  p  = dble(phi) + (i-1-LARG/2) * dble(dr)/dble(LARG/2)
  t  = dble(theta) + (j-1-LARG/2) * dble(dr)/dble(LARG/2)

  if (p < 0.d0-eps) p = p + 360.d0
  if (p >= 360.d0-eps) p = p - 360.d0

  if (t > 180.d0 - eps) then
    t = t - 180.d0
    p = 360.d0 - p
  else if (t < 0.d0 - eps) then
    t = 180.d0 + t
    p = 360.d0 - p
  endif

  if (p < 0.d0-eps) p = p + 360.d0
  if (p >= 360.d0-eps) p = p - 360.d0

  ii = nint(p/drin) + 1
  if (ii > NBP) ii = NBP

  jj = nint(t/drin) + 1
  if (jj > NBT) jj = NBT

  thetar = sngl(t)
  phir  = sngl(p)

  end subroutine get_indexloc

!
!--------------------------------------------------------------------------------------------------
!

  subroutine crust_bilinear_variable(theta,phi,gll_ind,rho_cr,vp_cr,vsv_cr,vsh_cr)

!theta phi en radians
!dr en degre
!reponse en metre

  use model_crust_berkeley_par
  use constants, only: PI

  implicit none
!    real, dimension(2) :: crust_array
  integer :: gll_ind
  double precision :: theta,phi,dx,dy,t,p,D_DEG
  double precision, intent(out) :: rho_cr,vp_cr,vsv_cr,vsh_cr
  integer :: broj_lokacija,j,m
  double precision, parameter :: deg2rad = PI/180.d0, rad2deg = 180.d0/PI

  ! sfrench 20110103 model grid spacing (degrees)
  D_DEG = 2.d0

  ! latitude in degrees
  t = 90.d0 - theta / deg2rad

  ! sfrench 20110103 new bounds on latitude in accord with new model grid spacing
  !if (t>89.d0)  t = 89.d0
  !if (t <-89.d0) t = -89.d0
  if (t > 90.d0 - D_DEG) t = 90.d0 - D_DEG
  if (t < D_DEG - 90.d0) t = D_DEG - 90.d0

  ! longitude in degrees
  p = phi / deg2rad
  if (p > 180.d0) p = p - 360.d0

  broj_lokacija = size(crust_array, DIM = 1)
  m = 0

  rho_cr = 0.d0
  vp_cr = 0.d0
  vsv_cr = 0.d0
  vsh_cr = 0.d0

  do j = 1,broj_lokacija
    if (m >= 4) then
      exit
    else
      if ( gll_ind == nint(crust_array(j,3) + 1) ) then
        dy = dabs(crust_array(j,1) - t)
        dy = dy / D_DEG ! sfrench 20110103 : normalized to model grid spacing

        if (dy < 1.d0) then
          dx = dabs(crust_array(j,2) - p)
          if (dx > 180.d0) dx = 360.d0 - dx

          dx = dx / D_DEG ! sfrench 20110103 : normalized to model grid spacing

          if (dabs(dx) < 1.d0) then
            ! Increment number of found locations. Must be <= 4
            m = m + 1
            rho_cr = rho_cr + (1.d0-dx) * (1.d0-dy) * crust_array(j,4)
            vp_cr  = vp_cr  + (1.d0-dx) * (1.d0-dy) * crust_array(j,5)
            vsv_cr = vsv_cr + (1.d0-dx) * (1.d0-dy) * crust_array(j,6)
            vsh_cr = vsh_cr + (1.d0-dx) * (1.d0-dy) * crust_array(j,7)

            !  print *,'m= ',m,'dx= ',dx,' dy= ',dy, 'Vsv= ',vsv_cr
          endif
        endif
      endif
    endif
  enddo

  end subroutine crust_bilinear_variable

!
!--------------------------------------------------------------------------------------------------
!

  subroutine read_crust_smooth_variable(unit)

  use model_crust_berkeley_par

  implicit none
  integer, intent(in) :: unit
  integer :: j,l,broj_lokacija

  read(unit,*) broj_lokacija

  ! sfrench 20110103 note: crust_array is < t> < p> < gll_ind> < rho> < vp> < vsv> < vsh>
  allocate(crust_array(broj_lokacija,7))
  crust_array(:,:) = -1000.0

  do j = 1,broj_lokacija
    read(unit,*) (crust_array(j,l),l=1,7)
  enddo

  !debug
  !print *,'read_crust_smooth_variable: I have read ',broj_lokacija,' crustal inputs!'

  end subroutine read_crust_smooth_variable

!
!--------------------------------------------------------------------------------------------------
!

  subroutine read_crustmoho_filtre(unit)

  use model_crust_berkeley_par

  implicit none
  integer, intent(in) :: unit
  integer :: i,j

  read(unit,*) NBP,NBT,drin
  if (drin /= 2.) STOP 'read_crust_filtre: dr muste be == 2'

  allocate(moho_start(NBP,NBT))
  moho_start(:,:) = -1000.0

  do j = 1,NBP
    do i = 1,NBT
      read(unit,*) moho_start(j,i)
    enddo
  enddo

  moho_start(:,:) = moho_start(:,:) / 1000.0
  NBT = NBT - 1

  end subroutine read_crustmoho_filtre

