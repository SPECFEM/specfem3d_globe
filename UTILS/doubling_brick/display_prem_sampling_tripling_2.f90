!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 2
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!          (c) California Institute of Technology July 2002
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

  program prem_iso

  implicit none

  include "constants.h"

! PREM
  double precision, parameter :: RMOHO = 6346600.d0
  double precision, parameter :: R220 = 6151000.d0
  double precision, parameter :: R400 = 5971000.d0
  double precision, parameter :: R600 = 5771000.d0
  double precision, parameter :: R670 = 5701000.d0
  double precision, parameter :: R771 = 5600000.d0
  double precision, parameter :: RTOPDDOUBLEPRIME = 3630000.d0
  double precision, parameter :: RCMB = 3480000.d0
  double precision, parameter :: RICB = 1221000.d0

! values common to PREM and IASP91
  double precision, parameter :: ROCEAN = 6368000.d0
  double precision, parameter :: RMIDDLE_CRUST = 6356000.d0
  double precision, parameter :: R80 = 6291000.d0

  integer i

  double precision r
  double precision x,rho,vp,vs,ratio_sampling,f0,lambda_P,lambda_S,delta_theta, &
      mesh_size,normalization_factor_P,normalization_factor_S
  double precision, dimension(1000) :: num_elements_per_lambda_P,num_elements_per_lambda_P_modif, &
      num_elements_per_lambda_S,num_elements_per_lambda_S_modif,radius,stability

!!!!!! DK DK exclude center of inner core to avoid division by zero
!!!!!! DK DK we have: 1000 * 1200/6371 = 188.35347
!!!!!! DK DK therefore the central cube ends around i = 188
!!!!!! DK DK  do i=1,1000
  integer, parameter :: NINIT = 140

  do i=NINIT,1000

! compute real physical radius in meters
  x = dble(i)/1000.d0

  r = x * R_EARTH

!
! PREM
!
!--- inner core
!
  if(r >= 0.d0 .and. r <= RICB) then
    rho=13.0885d0-8.8381d0*x*x
    vp=11.2622d0-6.3640d0*x*x
    vs=3.6678d0-4.4475d0*x*x
    ratio_sampling = 9.d0
!
!--- outer core
!
  else if(r > RICB .and. r <= RCMB) then
    rho=12.5815d0-1.2638d0*x-3.6426d0*x*x-5.5281d0*x*x*x
    vp=11.0487d0-4.0362d0*x+4.8023d0*x*x-13.5732d0*x*x*x
    vs=0.0d0

!! DK DK it seems optimal to implement the third doubling at a depth of 4300 km in the mantle
    if(r > R_EARTH - 3800000.d0) then
      ratio_sampling = 3.d0
    else
      ratio_sampling = 9.d0
    endif

!
!--- D" at the base of the mantle
!
  else if(r > RCMB .and. r <= RTOPDDOUBLEPRIME) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vp=15.3891d0-5.3181d0*x+5.5242d0*x*x-2.5514d0*x*x*x
    vs=6.9254d0+1.4672d0*x-2.0834d0*x*x+0.9783d0*x*x*x
    ratio_sampling = 3.d0
!
!--- mantle: from top of D" to d670
!
  else if(r > RTOPDDOUBLEPRIME .and. r <= R771) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vp=24.9520d0-40.4673d0*x+51.4832d0*x*x-26.6419d0*x*x*x
    vs=11.1671d0-13.7818d0*x+17.4575d0*x*x-9.2777d0*x*x*x

!! DK DK it seems optimal to implement the second doubling at a depth of 1500 km in the mantle
!   if(r > R_EARTH - 2500000.d0) then
      ratio_sampling = 3.d0
!   else
!     ratio_sampling = 3.d0
!   endif

  else if(r > R771 .and. r <= R670) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vp=29.2766d0-23.6027d0*x+5.5242d0*x*x-2.5514d0*x*x*x
    vs=22.3459d0-17.2473d0*x-2.0834d0*x*x+0.9783d0*x*x*x
    ratio_sampling = 3.d0
!
!--- mantle: above d670
!
  else if(r > R670 .and. r <= R600) then
    rho=5.3197d0-1.4836d0*x
    vp=19.0957d0-9.8672d0*x
    vs=9.9839d0-4.9324d0*x
    ratio_sampling = 3.d0
  else if(r > R600 .and. r <= R400) then
    rho=11.2494d0-8.0298d0*x
    vp=39.7027d0-32.6166d0*x
    vs=22.3512d0-18.5856d0*x
    ratio_sampling = 3.d0
  else if(r > R400 .and. r <= R220) then
    rho=7.1089d0-3.8045d0*x
    vp=20.3926d0-12.2569d0*x
    vs=8.9496d0-4.4597d0*x
    ratio_sampling = 3.d0
  else if(r > R220 .and. r <= R80) then
    rho=2.6910d0+0.6924d0*x
    vp=4.1875d0+3.9382d0*x
    vs=2.1519d0+2.3481d0*x
    ratio_sampling = 3.d0
  else
! use PREM crust
    if(r > R80 .and. r <= RMOHO) then
      rho=2.6910d0+0.6924d0*x
      vp=4.1875d0+3.9382d0*x
      vs=2.1519d0+2.3481d0*x
!!!!!!!!!!!!!!!      ratio_sampling = 3.d0
      ratio_sampling = 3.d0
    else if(r > RMOHO .and. r <= RMIDDLE_CRUST) then
      rho=2.9d0
      vp=6.8d0
      vs=3.9d0
      ratio_sampling = 1.d0
    else if(r > RMIDDLE_CRUST .and. r <= ROCEAN) then
      rho=2.6d0
      vp=5.8d0
      vs=3.2d0
      ratio_sampling = 1.d0
! for density profile for gravity, we do not check that r <= R_EARTH
    else if(r > ROCEAN) then
      rho=2.6d0
      vp=5.8d0
      vs=3.2d0
      ratio_sampling = 1.d0
    endif
  endif

  rho=rho*1000.0d0
  vp=vp*1000.0d0
  vs=vs*1000.0d0

! compute spatial wavelength at a given reference frequency f0
  f0 = 1.d0
  lambda_P = vp / f0
  lambda_S = vs / f0

! compute "horizontal" mesh size at that depth for a given angular size delta_theta
  delta_theta = 1.d0
  mesh_size = (r/1000.d0) * delta_theta * ratio_sampling

  radius(i) = (R_EARTH-r)/1000.d0

! compute number of elements per wavelength at that depth
  num_elements_per_lambda_P(i) = lambda_P / mesh_size
  num_elements_per_lambda_S(i) = lambda_S / mesh_size

! compute stability condition (Courant number)
  stability(i) = vp / mesh_size

!   print *,sngl((R_EARTH-r)/1000.d0),sngl(ratio_sampling),sngl(ratio_sampling)

  enddo

! store normalization factors
  num_elements_per_lambda_P_modif(:) = num_elements_per_lambda_P(:)
  where(num_elements_per_lambda_P < 0.001d0) num_elements_per_lambda_P_modif = 10000.d0
  normalization_factor_P = minval(num_elements_per_lambda_P_modif)!* 1.35d0

  num_elements_per_lambda_S_modif(:) = num_elements_per_lambda_S(:)
  where(num_elements_per_lambda_S < 0.001d0) num_elements_per_lambda_S_modif = 10000.d0
  normalization_factor_S = minval(num_elements_per_lambda_S_modif)!* 1.35d0

  do i=NINIT,1000
    print *,sngl(radius(i)),sngl(num_elements_per_lambda_P(i)/normalization_factor_P), &
      sngl(num_elements_per_lambda_S(i)/normalization_factor_S),sngl(stability(i))
  enddo

! print *
! print *,'use gnuplot to display radius,vp,vs,rho file'
! print *

  end program prem_iso

