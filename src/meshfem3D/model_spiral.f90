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

! SPiRaL model - by Nathan Simmons et al., LLNL
! Model inverted for the 5 Thomsen parameters
! This routine handles crust, mantle and 410-660 topography model files
!
! reference:
!   N. A. Simmons,  S. C. Myers,  C. Morency,  A. Chiang,  D. R. Knapp, 2021.
!   SPiRaL: a multiresolution global tomography model of seismic wave speeds
!   and radial anisotropy variations in the crust and mantle,
!   Geophysical Journal International, Volume 227, Issue 2, November 2021, Pages 1366-1391
!   https://doi.org/10.1093/gji/ggab277
!
!
! How to setup model in DATA/Par_file:
! You can use SPiRaL mantle and/or crust and other available models for mantle and/or crust.
! Examples:
! ---------
! Use SPiRaL for mantle and crust1.0 for crust:
!   MODEL                           = spiral_crust1.0
!
! Use SPiRaL for mantle and for crust:
!   MODEL                           = spiral_crustspiral     # or just: spiral
!
! Use S362ani for mantle and SPiRaL for crust:
!   MODEL                           = s362ani_crustspiral
!
!
! Implementation notes:
! - Version SPiRaL.v1.2: Christina Morency, LLNL, 07/14/2021
!   Updated to latest version of SPiRaL.1.4
!
! - Version SPiRaL.v1.1: Christina Morency, LLNL, 02/25/2020
!   Debugged mantle loop.
!   Added flag_smooth_spiral_crust to smooth the crust if needed
!
! - Version SPiRaL.v1.0: Christina Morency, LLNL, 01/21/2020
!   Using model SPiRaL.1.3.
!   Routine returns non-dimensionalized cij, although at the moment the
!   model is defined as radial isotropic, but we kept the general expression of cij
!   for future development.
!
!  LLNL-CODE-802993
!
!This work was produced under the auspices of the U.S. Department of
!Energy by Lawrence Livermore National Laboratory under Contract
!DE-AC52-07NA27344.
!
!This work was prepared as an account of work sponsored by an agency of
!the United States Government. Neither the United States Government nor
!Lawrence Livermore National Security, LLC, nor any of their employees
!makes any warranty, expressed or implied, or assumes any legal liability
!or responsibility for the accuracy, completeness, or usefulness of any
!information, apparatus, product, or process disclosed, or represents that
!its use would not infringe privately owned rights.
!
!Reference herein to any specific commercial product, process, or service
!by trade name, trademark, manufacturer, or otherwise does not necessarily
!constitute or imply its endorsement, recommendation, or favoring by the
!United States Government or Lawrence Livermore National Security, LLC.
!
!The views and opinions of authors expressed herein do not necessarily
!state or reflect those of the United States Government or Lawrence
!Livermore National Security, LLC, and shall not be used for advertising
!or product endorsement purposes.
!
!BSD 3-Clause License
!
!Copyright (c) 2018, Lawrence Livermore National Security, LLC
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
!
!* Redistributions of source code must retain the above copyright notice, this
!  list of conditions and the following disclaimer.
!
!* Redistributions in binary form must reproduce the above copyright notice,
!  this list of conditions and the following disclaimer in the documentation
!  and/or other materials provided with the distribution.
!
!* Neither the name of the copyright holder nor the names of its
!  contributors may be used to endorse or promote products derived from
!  this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



!-----------------------------------------------------------------------------------------
!
! Crustal model
!
!-----------------------------------------------------------------------------------------

! The 8 crustal layers:
! ====================
! 1) water
! 2) ice
! 3) upper sediments
! 4) middle sediments
! 5) lower sediments
! 6) upper crystalline crust
! 7) middle crystalline crust
! 8) lower crystalline crust
! + a ninth layer gives coefficents and rho below the Moho. The values
!   are associated with LLNL model G3Cv3 on continents and a thermal
!   model in the oceans.
!
! returns non-dimensionalized cij
!--------------------------------------------------------------------------------------------------

  module model_spiral_crust_par

  ! crustal_model_constants

  ! crustal model parameters
  ! spiral1.3_orig
  !integer, parameter :: Nbnd     = 8      ! number of bands in crust_bands_info.txt
  !integer, parameter :: CRUST_NP = 9      ! number of crustal layers
  !integer, parameter :: CRUST_NB = 241086 ! sum_bnd (Nlat_bnd*Nlon_bnd)

  ! spiral1.4
  integer, parameter :: Nbnd     = 13     ! number of bands in crust_bands_info.txt
  integer, parameter :: CRUST_NP = 9      ! number of crustal layers
  integer, parameter :: CRUST_NB = 248291 ! sum_bnd (Nlat_bnd*Nlon_bnd)

  ! model_crust_variables
  ! CIJ: 1=c11, 2=c13, 3=c33, 4=c44, 5=c66
  double precision, dimension(:,:,:), allocatable :: crust_coef

  ! layer thickness and density
  double precision, dimension(:,:), allocatable :: crust_thickness, crust_rho

  ! variables of the bands
  integer, dimension(:), allocatable :: bnd_nlon, bnd_nlat
  double precision, dimension(:), allocatable :: bnd_lat1, bnd_lat2

  ! crustal smoothing
  logical :: flag_smooth_spiral_crust
  ! sampling rate for CAP points
  integer, parameter :: NTHETA_spiral = 4, NPHI_spiral = 20

  ! model value interpolation (between different grid points)
  logical, parameter :: interpolate_crust = .true.

  end module model_spiral_crust_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_crust_spiral_broadcast()

! standard routine to setup model

  use constants
  use shared_parameters, only: NEX_XI,NEX_ETA
  use model_spiral_crust_par

  implicit none

  ! local parameters
  integer :: ier
  integer :: NEX_MAX

  ! allocate crustal arrays
  allocate(crust_thickness(CRUST_NP,CRUST_NB), &
           crust_coef(5,CRUST_NP,CRUST_NB), &
           crust_rho(CRUST_NP,CRUST_NB), &
           bnd_nlat(Nbnd), &
           bnd_nlon(Nbnd), &
           bnd_lat1(Nbnd), &
           bnd_lat2(Nbnd), &
           stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating crustal arrays')

  ! initializes
  crust_coef(:,:,:) = ZERO
  crust_rho(:,:) = ZERO
  crust_thickness(:,:) = ZERO
  bnd_nlat(:) = 0
  bnd_nlon(:) = 0
  bnd_lat1(:) = ZERO
  bnd_lat2(:) = ZERO

  ! by default, crustal smoothing is on
  ! (required for lower NEX models, otherwise negative Jacobian mesh gets produced)
  flag_smooth_spiral_crust = .true.

  ! turns smoothing off when NEX resolution is high enough to produce good quality mesh
  ! note from Christina:
  !   Crustal smoothing does affect waveforms amplitudes. The absence of smoothing better fits the data,
  !   and you can appreciate the differences between with/without crustal smoothing on amplitudes.
  !   So, we strongly recommend removing crustal smoothing when NEX >= 144.
  NEX_MAX = max(NEX_XI,NEX_ETA)
  if (NEX_MAX >= 144) flag_smooth_spiral_crust = .false.

  ! the variables read are declared and stored in structure model_spiral_par
  if (myrank == 0) call read_crust_spiral_model()

  ! broadcast the information read on the master to the nodes
  call bcast_all_dp(crust_thickness,CRUST_NP*CRUST_NB)
  call bcast_all_dp(crust_coef,5*CRUST_NP*CRUST_NB)
  call bcast_all_dp(crust_rho,CRUST_NP*CRUST_NB)

  call bcast_all_i(bnd_nlat,Nbnd)
  call bcast_all_i(bnd_nlon,Nbnd)
  call bcast_all_dp(bnd_lat1,Nbnd)
  call bcast_all_dp(bnd_lat2,Nbnd)

  end subroutine model_crust_spiral_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_crust_spiral_model()

  use constants
  !use shared_parameters, only: R_PLANET_KM,RHOAV

  use model_spiral_crust_par

  implicit none

  ! local variables
  integer :: ier
  integer :: i,k,b,Nbnd_read,NlatNlon,Nbnd_numpoints
  character(len=7), dimension(:),allocatable :: dlat,dlon,lat1,lat2
  ! layer depths
  double precision, dimension(:,:),allocatable :: bnd
  double precision :: moho
  double precision :: h_moho_min,h_moho_max
  double precision, dimension(:,:),allocatable :: single_par

  ! debug
  !double precision :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33, &
  !                   c34,c35,c36,c44,c45,c46,c55,c56,c66
  !double precision :: lat,lon,x,vpvc,vphc,vsvc,vshc,etac,rhoc
  !logical :: found_crust
  !integer, parameter :: CRUST_NLO = 360
  !integer, parameter :: CRUST_NLA = 180

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) 'incorporating crustal model: SPiRaL1.4'
  write(IMAIN,*)
  call flush_IMAIN()

  ! allocates temporary array
  allocate(dlat(Nbnd), &
           dlon(Nbnd), &
           lat1(Nbnd), &
           lat2(Nbnd), &
           bnd(CRUST_NP,CRUST_NB),stat=ier)
  if (ier /= 0 ) call exit_MPI(0,'Error allocating dlat, dlon, depths crustal arrays in read routine')
  dlat(:) = ""; dlon(:) = ""; lat1(:) = ""; lat2(:) = ""

  ! initializes
  bnd(:,:) = ZERO

  ! opens crust info file
  open(IIN,file='DATA/spiral1.4/crust/crust_bands_info.txt',action='read',status='old',iostat=ier)
  if (ier /= 0) then
    write(IMAIN,*) 'Error opening "DATA/spiral1.4/crust/crust_bands_info.txt": ', ier
    call exit_MPI(0,'Error model spiral1.4: file not found DATA/spiral1.4/crust/crust_bands_info.txt')
  endif

  ! checking the number of bands
  read(IIN,*) Nbnd_read
  if (Nbnd_read /= Nbnd) then
    write(IMAIN,*) 'Error opening "DATA/spiral1.4/crust/crust_bands_info.txt", wrong number of bands: ', Nbnd_read
    call exit_MPI(0,'Error model spiral1.4: wrong number of bands in DATA/spiral1.4/crust/crust_bands_info.txt')
  endif

  ! read the variables for each band
  do b = 1,Nbnd_read
      read(IIN,*) lat1(b),lat2(b),dlat(b),dlon(b),bnd_nlat(b),bnd_nlon(b)
  enddo
  ! closes files
  close(IIN)

  read(lat1(:),*) bnd_lat1(:)
  read(lat2(:),*) bnd_lat2(:)

  ! temporary array for reading single parameter
  allocate(single_par(CRUST_NP,CRUST_NB),stat=ier)
  if (ier /= 0) call exit_MPI(0,'Error allocating temporary array single_par')
  single_par(:,:) = 0.d0

  NlatNlon = 0

  ! read CIJ, density and depths from crust files
  do b = 1,Nbnd_read
    Nbnd_numpoints = bnd_nlat(b) * bnd_nlon(b)

    ! checks if anything to do in this band
    if (Nbnd_numpoints == 0) cycle

    call read_general_crust_model(single_par,'C11    ', &
                                  lat1(b),lat2(b),dlat(b),dlon(b),NlatNlon)
    crust_coef(1,:,NlatNlon+1:NlatNlon+Nbnd_numpoints) = single_par(:,NlatNlon+1:NlatNlon+Nbnd_numpoints)

    call read_general_crust_model(single_par,'C13    ', &
                                  lat1(b),lat2(b),dlat(b),dlon(b),NlatNlon)
    crust_coef(2,:,NlatNlon+1:NlatNlon+Nbnd_numpoints) = single_par(:,NlatNlon+1:NlatNlon+Nbnd_numpoints)

    call read_general_crust_model(single_par,'C33    ', &
                                  lat1(b),lat2(b),dlat(b),dlon(b),NlatNlon)
    crust_coef(3,:,NlatNlon+1:NlatNlon+Nbnd_numpoints) = single_par(:,NlatNlon+1:NlatNlon+Nbnd_numpoints)

    call read_general_crust_model(single_par,'C44    ', &
                                  lat1(b),lat2(b),dlat(b),dlon(b),NlatNlon)
    crust_coef(4,:,NlatNlon+1:NlatNlon+Nbnd_numpoints) = single_par(:,NlatNlon+1:NlatNlon+Nbnd_numpoints)

    call read_general_crust_model(single_par,'C66    ', &
                                  lat1(b),lat2(b),dlat(b),dlon(b),NlatNlon)
    crust_coef(5,:,NlatNlon+1:NlatNlon+Nbnd_numpoints) = single_par(:,NlatNlon+1:NlatNlon+Nbnd_numpoints)

    call read_general_crust_model(bnd(:,:),'depths ', &
                                  lat1(b),lat2(b),dlat(b),dlon(b),NlatNlon)

    call read_general_crust_model(crust_rho(:,:),'density', &
                                  lat1(b),lat2(b),dlat(b),dlon(b),NlatNlon)

    NlatNlon = NlatNlon + Nbnd_numpoints
  enddo ! reading files

  ! frees memory
  deallocate(lat1)
  deallocate(lat2)
  deallocate(dlat)
  deallocate(dlon)
  deallocate(single_par)

  ! checks if all read
  if (NlatNlon /= CRUST_NB) then
    print *,'Error: spiral crust total number of parameters invalid',NlatNlon,' should be ',CRUST_NB
    call exit_MPI(0,'Error spiral invalid total number of parameters')
  endif

  h_moho_min = HUGEVAL
  h_moho_max = -HUGEVAL

  ! determines layer thickness
  do i = 1,CRUST_NB
    do k = 1,CRUST_NP - 1
      crust_thickness(k,i) = - (bnd(k+1,i) - bnd(k,i))
    enddo

    ! thickness = ice (layer index 2) + sediment (index 3+4+5) + crystalline crust (index 6+7+8)
    ! crustal thickness without ice
    ! note: etopo1 has topography including ice ("ice surface" version) and at base of ice sheets ("bedrock" version)
    !       see: http://www.ngdc.noaa.gov/mgg/global/global.html
    moho = crust_thickness(3,i) + crust_thickness(4,i) + crust_thickness(5,i) &
           + crust_thickness(6,i) + crust_thickness(7,i) + crust_thickness(8,i)

    ! limit moho thickness
    if (moho > h_moho_max) h_moho_max = moho
    if (moho < h_moho_min) h_moho_min = moho
  enddo

  ! frees memory
  deallocate(bnd)

  ! user output
  write(IMAIN,*) '  Moho crustal thickness (without ice) min/max = ',sngl(h_moho_min),sngl(h_moho_max),' km'
  write(IMAIN,*)
  call flush_IMAIN()

  ! checks min/max
  if (h_moho_min == HUGEVAL .or. h_moho_max == -HUGEVAL) stop 'incorrect moho depths in read_crust_spiral_model'

  ! debug: file output for smoothed data
  !open(77,file='tmp-crust.dat',status='unknown')
  !write(77,*)'#crustal thickness: #lat (degree) #lon (degree) #moho (km) (w/out ice) #vp (at surface) #vs (at surface)'
  ! smoothed version
  !!    do k = 1,CRUST_NLA
  !!      lat = 90.d0 - k + 0.5d0
  !      do i = 1,CRUST_NLO
  !        lon = -180.d0 + i - 0.5d0
  !    do k = 1,CRUST_NLA
  !      lat = 90.d0 - k + 0.5d0
  !        x = 1.0d0
  !!      do i=1,180
  !!      x=0.90 + (1.0-0.90)*(i-1)/(180-1)
  !!      lat=40.00
  !!      lon=-99.00
  !!    write(IMAIN,*) 'k,i,lat,lon',k,i,lat,lon
  !!        call model_crust_spiral(lat,lon,x,vpvc,vphc,vsvc,vshc,etac,rhoc,moho, &
  !!                              c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
  !!                              c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,found_crust,.false.)
  !!        write(77,*) x, vphc*(R_PLANET_KM*dsqrt(PI*GRAV*RHOAV))
  !        write(77,'(6d10.2)')lat,lon,moho*R_PLANET_KM, &
  !         vphc*(R_PLANET_KM*dsqrt(PI*GRAV*RHOAV)),vshc*(R_PLANET_KM*dsqrt(PI*GRAV*RHOAV)),rhoc*(RHOAV/1000.0d0)
  !!      enddo
  !    enddo
  !close(77)

  ! user output
  if (flag_smooth_spiral_crust) then
    write(IMAIN,*) '  using crustal smoothing:'
    write(IMAIN,*) '    CAP smoothing degree default                  = ',sngl(CAP_SMOOTHING_DEGREE_DEFAULT)
    write(IMAIN,*) '    smooth crust even more (for critical regions) = ',SMOOTH_CRUST_EVEN_MORE
  else
    write(IMAIN,*) '  no crustal smoothing'
  endif
  write(IMAIN,*)
  call flush_IMAIN()

  end subroutine read_crust_spiral_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_general_crust_model(var,var_letter,lat1,lat2,dlat,dlon,NlatNlon)

  use constants
  use model_spiral_crust_par, only: CRUST_NP,CRUST_NB

  implicit none

  double precision, intent(inout), dimension(CRUST_NP,CRUST_NB) :: var
  character(len=7), intent(in) :: var_letter
  character(len=7), intent(in) :: lat1,lat2,dlat,dlon
  integer, intent(in) :: NlatNlon

  ! local variables
  character(len=MAX_STRING_LEN) :: filecrust, line
  integer :: ier, ila, ilo, Nlat, Nlon, k, nlines
  double precision :: dlat1,dlat2,ddlat,ddlon

  write(filecrust,'(a37,a,a6,a,a6,a,a6,a,a1,a)') &
    'DATA/spiral1.4/crust/crust_band_lat1_', &
    trim(lat1),'_lat2_',trim(lat2),'_dlat_',trim(dlat),'_dlon_',trim(dlon),'.',var_letter

  read(lat1,*) dlat1
  read(lat2,*) dlat2
  read(dlat,*) ddlat
  read(dlon,*) ddlon

  Nlat = int((dlat2 - dlat1)/ddlat) + 1
  Nlon = int((180 + 180)/ddlon) + 1

  ! debug
  !print *,'debug: reading', trim(filecrust), ' Nlat, Nlon, NlatNlon=',Nlat, Nlon,NlatNlon

  open(unit = IIN,file=trim(filecrust),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    write(IMAIN,*) 'Error opening "', trim(filecrust), '": ', ier
    call flush_IMAIN()
    ! stop
    call exit_MPI(0, 'Error crust model spiral1.4')
  endif

  ! count number of lines
  nlines = 0
  do while (ier == 0)
    read(IIN,*,iostat=ier) line
    if (ier == 0) nlines = nlines + 1
  enddo
  rewind(IIN)

  ! checks line count
  if (nlines /= Nlon * Nlat) then
    print *,'Error: file ',trim(filecrust),' has wrong number of lines ',nlines,' should be ',Nlon * Nlat
    call exit_MPI(0,'Error spiral mantle file w/ wrong number of lines')
  endif

  ! checks crustal layers
  if (CRUST_NP /= 9) then
    print *,'Error: routine read_general_crust_model() assumes a 9-layered crust, but was compiled with ',CRUST_NP
    call exit_MPI(0,'Error spiral crust has invalid number of crustal layers in routine read_general_crust_model()')
  endif

  k = NlatNlon
  do ila = 1,Nlat
    do ilo = 1,Nlon
      k = k+1
      read(IIN,*) var(1,k),var(2,k),var(3,k),var(4,k),var(5,k),var(6,k),var(7,k),var(8,k), var(9,k)
    enddo
  enddo
  close(IIN)

  ! debug
  !print *,'debug: spiral crust reading', trim(filecrust), ' NlatNlon+Nlat*Nlon=',k

  end subroutine read_general_crust_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_crust_spiral(lat_in,lon_in,x_in,vpvc,vphc,vsvc,vshc,etac,rhoc,moho,sediment, &
                                c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                                found_crust,elem_in_crust,moho_only)

  use constants
  use shared_parameters, only: R_PLANET,R_PLANET_KM,RHOAV

  use model_spiral_crust_par

  implicit none

  double precision,intent(in) :: lat_in,lon_in,x_in ! lat = [-90;90], lon = [-180;180]
  double precision,intent(out) :: vpvc,vphc,vsvc,vshc,etac,rhoc
  double precision,intent(out) :: moho,sediment
  double precision,intent(out) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  logical,intent(out) :: found_crust
  logical,intent(in) :: elem_in_crust,moho_only

  ! local parameters
  double precision :: lat,lon,x
  double precision :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66
  double precision :: scale_GPa,scaleval,theta,phi
  double precision :: thicks_2
  double precision :: h_sed,h_uc
  double precision :: x2,x3,x4,x5,x6,x7,x8
  double precision,dimension(CRUST_NP) :: rhosl,thicksl,rhos,thicks
  double precision,dimension(5,CRUST_NP) :: coefsl,coefs
  double precision,dimension(NTHETA_spiral*NPHI_spiral) :: x1,y1,weight
  double precision :: weightl,cap_degree_spiral
  double precision :: dist
  integer :: k,ii

  !-------------------------------
  ! work-around to avoid Jacobian problems when stretching mesh elements;
  ! one could also try to slightly change the shape of the doubling element bricks (which cause the problem)...
  !
  ! defines a "critical" region around the andes to have at least a 2-degree smoothing;
  ! critical region can lead to negative Jacobians for mesh stretching when CAP smoothing is too small
  double precision,parameter :: LAT_CRITICAL_ANDES = -20.0d0
  double precision,parameter :: LON_CRITICAL_ANDES = -70.0d0
  double precision,parameter :: CRITICAL_RANGE = 70.0d0
  !-------------------------------

  ! gets position
  lat = lat_in
  lon = lon_in
  x = x_in

  ! Should be [-90,90], [-180,180] convention
  if (lat > 90.0d0 .or. lat < -90.0d0 .or. lon > 180.0d0 .or. lon < -180.0d0) then
    print *,'Error in lat/lon:',lat,lon,'in routine model_crust_spiral()'
    stop 'Error in latitude/longitude range in model_crust_spiral()'
  endif

  ! makes sure lat/lon are within range
  if (lat == 90.0d0) lat = 89.9999d0
  if (lat == -90.0d0) lat = -89.9999d0
  if (lon == 180.0d0) lon = 179.9999d0
  if (lon == -180.0d0) lon = -179.9999d0

  ! initialization
  c11 = ZERO
  c12 = ZERO
  c13 = ZERO
  c14 = ZERO
  c15 = ZERO
  c16 = ZERO
  c22 = ZERO
  c23 = ZERO
  c24 = ZERO
  c25 = ZERO
  c26 = ZERO
  c33 = ZERO
  c34 = ZERO
  c35 = ZERO
  c36 = ZERO
  c44 = ZERO
  c45 = ZERO
  c46 = ZERO
  c55 = ZERO
  c56 = ZERO
  c66 = ZERO

  vsvc = ZERO
  vshc = ZERO
  vpvc = ZERO
  vphc = ZERO
  rhoc = ZERO
  etac = ZERO

  moho = ZERO
  sediment = ZERO
  found_crust = .true.

  ! gets crustal structure at position lat/lon
  if (.not. flag_smooth_spiral_crust) then
    ! no smoothing
    call read_crust_spiral(lat,lon,coefs,rhos,thicks)
  else
    ! smoothing
    ! sets up smoothing points based on cap smoothing
    cap_degree_spiral = CAP_SMOOTHING_DEGREE_DEFAULT

    ! checks if inside/outside of critical region for mesh stretching
    if (SMOOTH_CRUST_EVEN_MORE) then
      dist = dsqrt( (lon-LON_CRITICAL_ANDES)**2 + (lat-LAT_CRITICAL_ANDES )**2 )
      if (dist < CRITICAL_RANGE) then
        ! increases cap smoothing degree
        ! scales between -1 at center and 0 at border
        dist = dist / CRITICAL_RANGE - ONE
        ! shifts value to 1 at center and 0 to the border with exponential decay
        dist = ONE - exp( - dist*dist*10.0d0 )
        ! increases smoothing degree inside of critical region to 2 degree
        cap_degree_spiral = cap_degree_spiral + dist
      endif
    endif

    ! gets smoothing points and weights
    call spiral_crust_smooth(lon,lat,x1,y1,weight,cap_degree_spiral)

    coefs(1,:) = ZERO
    coefs(2,:) = ZERO
    coefs(3,:) = ZERO
    coefs(4,:) = ZERO
    coefs(5,:) = ZERO
    thicks(:) = ZERO
    rhos(:) = ZERO

    do k = 1,NTHETA_spiral * NPHI_spiral
      call read_crust_spiral(y1(k),x1(k),coefsl,rhosl,thicksl)

      ! sediment thickness
      h_sed = thicksl(3) + thicksl(4) + thicksl(5)

      ! takes upper crust value if sediment too thin
      if (h_sed < MINIMUM_SEDIMENT_THICKNESS) then
        ! layers: 3 - sediment upper
        !         4 - sediment middle
        !         5 - sediment lower
        !         6 - upper crust
        do ii = 3,5
          coefsl(1,ii) = coefsl(1,6)
          coefsl(2,ii) = coefsl(2,6)
          coefsl(3,ii) = coefsl(3,6)
          coefsl(4,ii) = coefsl(4,6)
          coefsl(5,ii) = coefsl(5,6)
          rhosl(ii) = rhosl(6)
        enddo
      endif

      ! weighting value
      weightl = weight(k)

      ! total, smoothed values
      coefs(1,:) = coefs(1,:) + weightl*coefsl(1,:)
      coefs(2,:) = coefs(2,:) + weightl*coefsl(2,:)
      coefs(3,:) = coefs(3,:) + weightl*coefsl(3,:)
      coefs(4,:) = coefs(4,:) + weightl*coefsl(4,:)
      coefs(5,:) = coefs(5,:) + weightl*coefsl(5,:)
      thicks(:) = thicks(:) + weightl*thicksl(:)
      rhos(:) = rhos(:) + weightl*rhosl(:)
    enddo
  endif

  ! note: for seismic wave propagation in general we ignore the water and ice
  ! sheets (oceans are re-added later as an ocean load)
  ! note: but for gravity integral calculations we include the ice
  if (INCLUDE_ICE_IN_CRUST) then
    thicks_2 = thicks(2)
  else
    thicks_2 = ZERO
  endif

  ! whole sediment thickness (with ice if included)
  h_sed = thicks_2 + thicks(3) + thicks(4) + thicks(5)

  ! upper crust thickness (including sediments above, and also ice if included)
  h_uc = h_sed + thicks(6)

  ! non-dimensionalization factor
  scaleval = ONE / R_PLANET_KM

  ! non-dimensionalize thicknesses (given in km)
  ! ice
  x2 = ONE - thicks_2 * scaleval
  ! upper sediment
  x3 = ONE - (thicks_2 + thicks(3)) * scaleval
  ! middle sediment
  x4 = ONE - (thicks_2 + thicks(3) + thicks(4)) * scaleval
  ! all sediments
  x5 = ONE - h_sed * scaleval
  ! upper crust
  x6 = ONE - h_uc * scaleval
  ! middle crust
  x7 = ONE - (h_uc + thicks(7)) * scaleval
  ! lower crust
  x8 = ONE - (h_uc + thicks(7) + thicks(8)) * scaleval

  ! no matter if found_crust is true or false, compute moho thickness
  moho = (h_uc + thicks(7) + thicks(8)) * scaleval

  ! checks if anything further to do
  if (moho_only) return

  ! sediment thickness
  if (INCLUDE_SEDIMENTS_IN_CRUST) then
    sediment = h_sed * scaleval
  endif

  ! initializes local aniso values
  d11 = ZERO
  d12 = ZERO
  d13 = ZERO
  d14 = ZERO
  d15 = ZERO
  d16 = ZERO
  d22 = ZERO
  d23 = ZERO
  d24 = ZERO
  d25 = ZERO
  d26 = ZERO
  d33 = ZERO
  d34 = ZERO
  d35 = ZERO
  d36 = ZERO
  d44 = ZERO
  d45 = ZERO
  d46 = ZERO
  d55 = ZERO
  d56 = ZERO
  d66 = ZERO

  ! gets corresponding crustal velocities and density
  if (x > x2 .and. INCLUDE_ICE_IN_CRUST) then
    d11 = coefs(1,2)
    d13 = coefs(2,2)
    d33 = coefs(3,2)
    d44 = coefs(4,2)
    d66 = coefs(5,2)
    d12 = d11 - 2.d0*d66
    d22 = d11
    d23 = d13
    d55 = d44
    rhoc = rhos(2)
  else if (x > x3 .and. INCLUDE_SEDIMENTS_IN_CRUST) then
    d11 = coefs(1,3)
    d13 = coefs(2,3)
    d33 = coefs(3,3)
    d44 = coefs(4,3)
    d66 = coefs(5,3)
    d12 = d11 - 2.d0*d66
    d22 = d11
    d23 = d13
    d55 = d44
    rhoc = rhos(3)
  else if (x > x4 .and. INCLUDE_SEDIMENTS_IN_CRUST) then
    d11 = coefs(1,4)
    d13 = coefs(2,4)
    d33 = coefs(3,4)
    d44 = coefs(4,4)
    d66 = coefs(5,4)
    d12 = d11 - 2.d0*d66
    d22 = d11
    d23 = d13
    d55 = d44
    rhoc = rhos(4)
  else if (x > x5 .and. INCLUDE_SEDIMENTS_IN_CRUST) then
    d11 = coefs(1,5)
    d13 = coefs(2,5)
    d33 = coefs(3,5)
    d44 = coefs(4,5)
    d66 = coefs(5,5)
    d12 = d11 - 2.d0*d66
    d22 = d11
    d23 = d13
    d55 = d44
    rhoc = rhos(5)
  else if (x > x6) then
    d11 = coefs(1,6)
    d13 = coefs(2,6)
    d33 = coefs(3,6)
    d44 = coefs(4,6)
    d66 = coefs(5,6)
    d12 = d11 - 2.d0*d66
    d22 = d11
    d23 = d13
    d55 = d44
    rhoc = rhos(6)
  else if (x > x7) then
    d11 = coefs(1,7)
    d13 = coefs(2,7)
    d33 = coefs(3,7)
    d44 = coefs(4,7)
    d66 = coefs(5,7)
    d12 = d11 - 2.d0*d66
    d22 = d11
    d23 = d13
    d55 = d44
    rhoc = rhos(7)
  else if (x > x8 .or. elem_in_crust) then
    ! takes lower crustal values only if x is slightly above moho depth or
    ! if elem_in_crust is set
    !
    ! note: it looks like this does distinguish between GLL points at the exact
    ! moho boundary,
    !          where the point is on the interface between both,
    !          oceanic elements and mantle elements below
    d11 = coefs(1,8)
    d13 = coefs(2,8)
    d33 = coefs(3,8)
    d44 = coefs(4,8)
    d66 = coefs(5,8)
    d12 = d11 - 2.d0*d66
    d22 = d11
    d23 = d13
    d55 = d44
    rhoc = rhos(8)
  else
    ! note: if x is exactly the moho depth this will return false
    found_crust = .false.
  endif

  ! non-dimensionalize
  if (found_crust) then
    ! TI velocities
    vphc = sqrt(d11/rhoc)
    vpvc = sqrt(d33/rhoc)
    vshc = sqrt(d66/rhoc)
    vsvc = sqrt(d44/rhoc)
    etac = d13/(d11 - 2.d0*d44)

    scaleval = ONE / ( R_PLANET_KM * dsqrt(PI*GRAV*RHOAV) )
    vsvc = vsvc * scaleval
    vshc = vshc * scaleval
    vpvc = vpvc * scaleval
    vphc = vphc * scaleval
    rhoc = rhoc * 1000.0d0 / RHOAV

    ! debug
    !print *,'debug: vsvc,vshc,vpvc,vphc,rhoc,etac,moho'
    !print *,vsvc,vshc,vpvc,vphc,rhoc,etac,moho

    ! non-dimensionalize the elastic coefficients using
    ! the scale of GPa--[g/cm^3][(km/s)^2]
    scaleval = dsqrt(PI*GRAV*RHOAV)
    scale_GPa = (RHOAV/1000.d0)*((R_PLANET*scaleval/1000.d0)**2)

    d11 = d11/scale_GPa
    d12 = d12/scale_GPa
    d13 = d13/scale_GPa
    d14 = d14/scale_GPa
    d15 = d15/scale_GPa
    d16 = d16/scale_GPa
    d22 = d22/scale_GPa
    d23 = d23/scale_GPa
    d24 = d24/scale_GPa
    d25 = d25/scale_GPa
    d26 = d26/scale_GPa
    d33 = d33/scale_GPa
    d34 = d34/scale_GPa
    d35 = d35/scale_GPa
    d36 = d36/scale_GPa
    d44 = d44/scale_GPa
    d45 = d45/scale_GPa
    d46 = d46/scale_GPa
    d55 = d55/scale_GPa
    d56 = d56/scale_GPa
    d66 = d66/scale_GPa

    ! assign the local (d_ij) or global (c_ij) anisotropic parameters.
    ! The c_ij are the coefficients in the global
    ! reference frame used in SPECFEM3D.
    theta = (90.d0 - lat) * DEGREES_TO_RADIANS
    phi = lon * DEGREES_TO_RADIANS
    call rotate_tensor_radial_to_global(theta,phi, &
                                        d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                        d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                                        c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                        c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
 endif

 end subroutine model_crust_spiral

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_crust_spiral(lat,lon,coefs,rhos,thicks)

! crustal CIJ, layer thickness in km

  use constants
  use model_spiral_crust_par

  implicit none

  ! argument variables
  double precision,intent(in) :: lat,lon
  double precision,intent(out) :: rhos(CRUST_NP),thicks(CRUST_NP)
  double precision,intent(out) :: coefs(5,CRUST_NP)

  ! local variables
  integer :: NlatNlon,l,index,i
  double precision :: dlon,dlat

  ! interpolation variables
  double precision :: a,b  ! weights
  integer :: rec_read  ! position of the record to read in model.dat (direct access file)
  integer :: i_min, i_max, j_min, j_max ! upper and lower bound indices
  double precision, dimension(CRUST_NP) :: thick1,thick2,thick3,thick4   ! thickness corner values in model.dat
  double precision, dimension(CRUST_NP) :: rho1,rho2,rho3,rho4   ! rho corner values in model.dat
  double precision, dimension(5,CRUST_NP) :: coef1,coef2,coef3,coef4   ! CIJ corner values in model.dat

  ! checks latitude/longitude
  if (lat > 90.0d0 .or. lat < -90.0d0 .or. lon > 180.0d0 .or. lon < -180.0d0) then
    print *,'Error in lat/lon:',lat,lon
    stop 'Error in latitude/longitude range in crust spiral1.4'
  endif

  ! initializes
  coefs(:,:) = ZERO
  rhos(:) = ZERO
  thicks(:) = ZERO

  if (.not. interpolate_crust) then
    ! debug
    !print *,'debug: checking we are in the right loop', lon, lat

    ! with no interpolation
    NlatNlon = 0
    do l = 1,Nbnd
      !write(IMAIN,*) 'checking', bnd_lat1(l),bnd_lat2(l),lat, NlatNlon
      if (lat >= bnd_lat1(l) .and. lat <= bnd_lat2(l)) then
        dlat = (bnd_lat2(l) - bnd_lat1(l))/(bnd_nlat(l)-1)
        dlon = 360.d0/(bnd_nlon(l)-1)
        !write(IMAIN,*) 'checking', bnd_lat1(l),bnd_lat2(l),bnd_nlat(l),bnd_nlon(l),dlat,dlon
        index = NlatNlon + &
                  int( (lat - bnd_lat1(l))/dlat ) * bnd_nlon(l) + &
                  int( (lon + 180.d0)/dlon ) + 1
        !write(IMAIN,*) 'checking', bnd_lat1(l),bnd_lat2(l),lat,dlat,dlon,index

        do i = 1,CRUST_NP
          rhos(i) = crust_rho(i,index)
          thicks(i) = crust_thickness(i,index)
          coefs(:,i) = crust_coef(:,i,index)
        enddo

        ! debug
        !print *,'debug: checking', NlatNlon, lat, lon, index, rhos(2)

        return
      endif

      NlatNlon = NlatNlon + bnd_nlat(l)*bnd_nlon(l)
    enddo ! do l = 1,Nbnd, with no interpolation

  else
    ! with interpolation
    NlatNlon = 0
    do l = 1,Nbnd
      if (lat >= bnd_lat1(l) .and. lat <= bnd_lat2(l)) then
        dlat = (bnd_lat2(l) - bnd_lat1(l))/(bnd_nlat(l)-1)
        dlon = 360.d0/(bnd_nlon(l)-1)

        ! locate the corners
        i_min = int( (lon + 180.d0)/dlon ) + 1
        i_max = i_min + 1
        j_min = int( (lat -bnd_lat1(l))/dlat ) + 1
        j_max = j_min + 1

        ! (1) i_min, j_min
        rec_read = NlatNlon + (j_min-1)*bnd_nlon(l) + (i_min-1) + 1
        ! index bounds
        if (rec_read < 1) rec_read = 1
        if (rec_read > CRUST_NB) rec_read = CRUST_NB

        do i = 1,CRUST_NP
          rho1(i) = crust_rho(i,rec_read)
          thick1(i) = crust_thickness(i,rec_read)
          coef1(:,i) = crust_coef(:,i,rec_read)
        enddo

        ! (2) i_min, j_max
        rec_read = NlatNlon + (j_max-1)*bnd_nlon(l) + (i_min-1) + 1
        ! index bounds
        if (rec_read < 1) rec_read = 1
        if (rec_read > CRUST_NB) rec_read = CRUST_NB

        do i = 1,CRUST_NP
          rho2(i) = crust_rho(i,rec_read)
          thick2(i) = crust_thickness(i,rec_read)
          coef2(:,i) = crust_coef(:,i,rec_read)
        enddo

        ! (3) i_max, j_min
        rec_read = NlatNlon + (j_min-1)*bnd_nlon(l) + (i_max-1) + 1
        ! index bounds
        if (rec_read < 1) rec_read = 1
        if (rec_read > CRUST_NB) rec_read = CRUST_NB

        do i = 1,CRUST_NP
          rho3(i) = crust_rho(i,rec_read)
          thick3(i) = crust_thickness(i,rec_read)
          coef3(:,i) = crust_coef(:,i,rec_read)
        enddo

        ! (4) i_max, j_max
        rec_read = NlatNlon + (j_max-1)*bnd_nlon(l) + (i_max-1) + 1
        ! index bounds
        if (rec_read < 1) rec_read = 1
        if (rec_read > CRUST_NB) rec_read = CRUST_NB

        do i = 1,CRUST_NP
          rho4(i) = crust_rho(i,rec_read)
          thick4(i) = crust_thickness(i,rec_read)
          coef4(:,i) = crust_coef(:,i,rec_read)
        enddo

        ! weights
        a = ( lon + 180.d0 - (i_min-1)*dlon ) / dlon
        b = ( lat - bnd_lat1(l) - (j_min-1)*dlat ) / dlat

        do i = 1,CRUST_NP
          rhos(i) = rho1(i)*(1.-a)*(1.-b) + rho2(i)*(1.-a)*b + rho3(i)*a*(1.-b) + rho4(i)*a*b
          thicks(i) = thick1(i)*(1.-a)*(1.-b) + thick2(i)*(1.-a)*b + thick3(i)*a*(1.-b) + thick4(i)*a*b
          coefs(:,i) = coef1(:,i)*(1.-a)*(1.-b) + coef2(:,i)*(1.-a)*b + coef3(:,i)*a*(1.-b) + coef4(:,i)*a*b
        enddo

        return
      endif

      NlatNlon = NlatNlon + bnd_nlat(l)*bnd_nlon(l)
   enddo ! do l = 1,Nbnd, with interpolation

  endif ! if (.not. interpolation)

  end subroutine read_crust_spiral
!
!-----------------------------------------------------------------------------------------
!

  subroutine spiral_crust_smooth(x,y,x1,y1,weight,cap_degree_spiral)

! smooths with a cap of size CAP (in degrees)
! using NTHETA points in the theta direction (latitudinal)
! and NPHI in the phi direction (longitudinal).
! The cap is rotated to the North Pole.

  use constants
  use model_spiral_crust_par, only: NTHETA_spiral,NPHI_spiral

  implicit none

  ! INPUT & OUTPUT
  double precision,intent(in) :: x, y
  double precision,dimension(NTHETA_spiral*NPHI_spiral),intent(out) :: x1,y1,weight
  double precision,intent(in) :: cap_degree_spiral
  ! Local
  double precision :: CAP,dtheta,dphi,cap_area,dweight,pi_over_nphi,total,wght
  double precision :: theta,phi,sint,cost,sinp,cosp
  double precision :: r_rot,theta_rot,phi_rot
  double precision,dimension(3,3) :: rotation_matrix
  double precision,dimension(3) :: xx,xc
  integer :: i,j,k,itheta,iphi

  x1(:) = ZERO
  y1(:) = ZERO
  weight(:) = ZERO

  if (cap_degree_spiral < TINYVAL) then
    print *, 'Error cap:', cap_degree_spiral
    print *, 'lat/lon:', x,y
    stop 'Error cap_degree too small'
  endif

  CAP = cap_degree_spiral * DEGREES_TO_RADIANS
  dtheta = 0.5d0 * CAP / dble(NTHETA_spiral)
  dphi = TWO_PI / dble(NPHI_spiral)

  cap_area = TWO_PI * (1.0d0 - dcos(CAP))
  dweight = CAP / dble(NTHETA_spiral) * dphi / cap_area
  pi_over_nphi = PI / dble(NPHI_spiral)

  phi = x * DEGREES_TO_RADIANS
  theta = (90.0d0 - y) * DEGREES_TO_RADIANS

  sint = dsin(theta)
  cost = dcos(theta)
  sinp = dsin(phi)
  cosp = dcos(phi)

  rotation_matrix(1,1) = cosp*cost
  rotation_matrix(1,2) = -sinp
  rotation_matrix(1,3) = cosp*sint
  rotation_matrix(2,1) = sinp*cost
  rotation_matrix(2,2) = cosp
  rotation_matrix(2,3) = sinp*sint
  rotation_matrix(3,1) = -sint
  rotation_matrix(3,2) = ZERO
  rotation_matrix(3,3) = cost

  i = 0
  total = 0.0d0
  do itheta = 1,NTHETA_spiral
    theta = dble(2*itheta-1)*dtheta
    cost = dcos(theta)
    sint = dsin(theta)
    wght = sint*dweight

    do iphi = 1,NPHI_spiral
      i = i+1
      weight(i) = wght

      total = total+weight(i)
      phi = dble(2*iphi-1)*pi_over_nphi
      cosp = dcos(phi)
      sinp = dsin(phi)

      xc(1) = sint*cosp
      xc(2) = sint*sinp
      xc(3) = cost

      do j = 1,3
        xx(j) = 0.0d0
        do k = 1,3
          xx(j) = xx(j)+rotation_matrix(j,k)*xc(k)
        enddo
      enddo

      call xyz_2_rthetaphi_dble(xx(1),xx(2),xx(3),r_rot,theta_rot,phi_rot)
      call reduce(theta_rot,phi_rot)
      x1(i) = phi_rot * RADIANS_TO_DEGREES                   ! lon
      y1(i) = (PI_OVER_TWO - theta_rot) * RADIANS_TO_DEGREES ! lat
      if (x1(i) > 180.d0) x1(i) = x1(i) - 360.d0
    enddo
  enddo

  if (abs(total - 1.0d0) > 0.001d0) then
    print *,'Error cap:',total,cap_degree_spiral
    stop
  endif

  end subroutine spiral_crust_smooth


!-----------------------------------------------------------------------------------------
!
! Mantle model
!
!-----------------------------------------------------------------------------------------

  module model_spiral_mantle_par

  ! mantle_model_constants
  ! mantle model parameters
  integer, parameter :: Nbnd       = 13     ! number of bands in mantle_bands_info.txt
  integer, parameter :: MANTLE_NB  = 248291 ! sum_bnd (Nlat_bnd*Nlon_bnd)
  integer, parameter :: Nbndz      = 3      ! number of bands in mantle_dzones_info.txt
  integer, parameter :: MANTLE_NBZ = 70     ! sum_bnddz (Ndep_bnddz)

  integer, parameter :: TOPO_NLO   = 721    ! Lon dimension for topo d410 and d660
  integer, parameter :: TOPO_NLA   = 361    ! Lat dimension for topo d410 and d660
  integer, parameter :: TOPO_RES   = 2      ! Spacing for topo d410 and d660 is 0.5 degress

  ! model_mantle_variables
  ! CIJ: 1=c11, 2=c13, 3=c33, 4=c44, 5=c66
  double precision, dimension(:,:), allocatable :: mantle_coef
  double precision, dimension(:), allocatable :: mantle_rho
  double precision, dimension(:), allocatable :: mantle_d410, mantle_d660

  ! variables of the bands
  integer, dimension(:), allocatable :: mtle_bnd_nlon, mtle_bnd_nlat
  integer, dimension(:), allocatable :: mtle_bnd_ndep
  double precision, dimension(:), allocatable :: mtle_bnd_lat1, mtle_bnd_lat2
  double precision, dimension(:), allocatable :: mtle_bnd_dep1, mtle_bnd_dep2

  ! model interpolation (between different grid points)
  logical, parameter :: interpolate_mantle = .true.

  ! 410/660 topography interpolation
  ! needs to be .true. otherwise mesh points will collapse, and global indexing fails
  logical, parameter :: interpolate_topo = .true.

  end module model_spiral_mantle_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_mantle_spiral_broadcast()

! standard routine to setup model

  use constants
  use model_spiral_mantle_par

  implicit none

  integer :: ier

  ! allocate crustal arrays
  allocate(mantle_coef(5,MANTLE_NB*MANTLE_NBZ), &
           mantle_rho(MANTLE_NB*MANTLE_NBZ), &
           mtle_bnd_nlat(Nbnd), &
           mtle_bnd_nlon(Nbnd), &
           mtle_bnd_ndep(Nbndz), &
           mtle_bnd_lat1(Nbnd), &
           mtle_bnd_lat2(Nbnd), &
           mtle_bnd_dep1(Nbndz), &
           mtle_bnd_dep2(Nbndz), &
           mantle_d410(TOPO_NLO*TOPO_NLA), &
           mantle_d660(TOPO_NLO*TOPO_NLA),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating mantle arrays')

  ! initializes
  mantle_coef(:,:) = ZERO
  mantle_rho(:) = ZERO
  mtle_bnd_nlat(:) = 0
  mtle_bnd_nlon(:) = 0
  mtle_bnd_ndep(:) = 0
  mtle_bnd_lat1(:) = ZERO
  mtle_bnd_lat2(:) = ZERO
  mtle_bnd_dep1(:) = ZERO
  mtle_bnd_dep2(:) = ZERO
  mantle_d410(:) = ZERO
  mantle_d660(:) = ZERO

  ! the variables read are declared and stored in structure model_spiral_par
  if (myrank == 0) call read_mantle_spiral_model()

  ! broadcast the information read on the master to the nodes
  call bcast_all_dp(mantle_coef,5*MANTLE_NB*MANTLE_NBZ)
  call bcast_all_dp(mantle_rho,MANTLE_NB*MANTLE_NBZ)

  call bcast_all_i(mtle_bnd_nlat,Nbnd)
  call bcast_all_i(mtle_bnd_nlon,Nbnd)
  call bcast_all_i(mtle_bnd_ndep,Nbndz)
  call bcast_all_dp(mtle_bnd_lat1,Nbnd)
  call bcast_all_dp(mtle_bnd_lat2,Nbnd)
  call bcast_all_dp(mtle_bnd_dep1,Nbndz)
  call bcast_all_dp(mtle_bnd_dep2,Nbndz)

  call bcast_all_dp(mantle_d410,TOPO_NLO*TOPO_NLA)
  call bcast_all_dp(mantle_d660,TOPO_NLO*TOPO_NLA)

  end subroutine model_mantle_spiral_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mantle_spiral_model()

  use constants
  use model_spiral_mantle_par

  implicit none

  ! local variables
  integer :: ier, ila, ilo
  integer :: k,b,m,Nbnd_read,Nbndz_read,NlatNlonNdep,Ndep,NLA,NLO,Nbnd_numpoints
  character(len=7), dimension(:),allocatable :: dlat,dlon,ddep,lat1,lat2,dep1,dep2
  double precision, dimension(:),allocatable :: dp_ddep
  double precision :: dummy
  double precision :: d410_min,d410_max
  double precision :: d660_min,d660_max
  double precision, dimension(:),allocatable :: single_par

  ! debug
  !double precision :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33, &
  !                   c34,c35,c36,c44,c45,c46,c55,c56,c66
  !double precision :: lat,lon,r,vpv,vph,vsv,vsh,eta,rho
  !integer :: i
  !integer, parameter :: NLO = 360
  !integer, parameter :: NLA = 180

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) 'broadcast mantle model: SPiRaL1.4'
  write(IMAIN,*)
  call flush_IMAIN()

  ! allocates temporary array
  allocate(dlat(Nbnd), &
           dlon(Nbnd), &
           ddep(Nbndz), &
           lat1(Nbnd), &
           lat2(Nbnd), &
           dep1(Nbndz), &
           dep2(Nbndz), &
           dp_ddep(Nbndz),stat=ier)
  if (ier /= 0 ) call exit_MPI(0,'Error allocating dlat, dlon, ddepth mantle arrays in read routine')
  dlat(:) = ""; dlon(:) = ""; ddep(:) = ""; lat1(:) = ""; lat2(:) = ""; dep1(:) = ""; dep2(:) = ""
  dp_ddep(:) = 0.d0

  ! opens mantle info files
  open(IIN,file='DATA/spiral1.4/mantle/mantle_bands_info.txt',action='read',status='old',iostat=ier)
  if (ier /= 0) then
    write(IMAIN,*) 'Error opening "DATA/spiral1.4/mantle/mantle_bands_info.txt": ', ier
    call exit_MPI(0,'Error model spiral1.4: file not found DATA/spiral1.4/mantle/mantle_bands_info.txt')
  endif

  ! checking the number of bands
  read(IIN,*) Nbnd_read
  if (Nbnd_read /= Nbnd) then
    write(IMAIN,*) 'Error opening "DATA/spiral1.4/mantle/mantle_bands_info.txt", wrong number of bands: ', Nbnd_read
    call exit_MPI(0,'Error model spiral1.4: wrong number of bands in DATA/spiral1.4/mantle/mantle_bands_info.txt')
  endif

  ! read the variables for each band
  do b = 1,Nbnd_read
    read(IIN,*) lat1(b),lat2(b),dlat(b),dlon(b),mtle_bnd_nlat(b),mtle_bnd_nlon(b)
    ! debug
    !print *,'debug: bands: ', lat1(b),lat2(b),dlat(b),dlon(b),mtle_bnd_nlat(b),mtle_bnd_nlon(b)
  enddo
  ! closes files
  close(IIN)

  read(lat1(:),*) mtle_bnd_lat1(:)
  read(lat2(:),*) mtle_bnd_lat2(:)
  ! debug
  !print *,'debug: bands: mtle_bnd_lat1 ',mtle_bnd_lat1(:),'mtle_bnd_lat2',mtle_bnd_lat2(:)

  ! opens mantle zone info files
  open(IIN,file='DATA/spiral1.4/mantle/mantle_dzones_info.txt',action='read',status='old',iostat=ier)
  if (ier /= 0) then
    write(IMAIN,*) 'Error opening "DATA/spiral1.4/mantle/mantle_dzones_info.txt": ', ier
    call exit_MPI(0,'Error model spiral1.4: file not found DATA/spiral1.4/mantle/mantle_dzones_info.txt')
  endif

  ! checking the number of bands
  read(IIN,*) Nbndz_read
  if (Nbndz_read /= Nbndz) then
    write(IMAIN,*) 'Error opening "DATA/spiral1.4/mantle/mantle_dzones_info.txt", wrong number of bands: ', Nbndz_read
    call exit_MPI(0,'Error model spiral1.4: wrong number of bands in DATA/spiral1.4/mantle/mantle_dzones_info.txt')
  endif

  ! read the variables for each band
  do b = 1,Nbndz_read
    read(IIN,*) dep1(b),dep2(b),ddep(b)
    ! debug
    !print *,'debug: dzones: ',dep1(b),dep2(b),ddep(b)
  enddo
  ! closes files
  close(IIN)

  read(dep1(:),*) mtle_bnd_dep1(:)
  read(dep2(:),*) mtle_bnd_dep2(:)
  read(ddep(:),*) dp_ddep(:)

  mtle_bnd_ndep(:) = int( (mtle_bnd_dep2(:) - mtle_bnd_dep1(:))/dp_ddep(:) ) + 1
  ! debug
  !print *,'debug: dzones: mtle_bnd_ndep ',mtle_bnd_ndep(:), &
  !        'mtle_bnd_dep1',mtle_bnd_dep1(:),'mtle_bnd_dep2',mtle_bnd_dep2(:),'ddep',dp_ddep(:)

  ! temporary array for reading single parameter
  allocate(single_par(MANTLE_NB*MANTLE_NBZ),stat=ier)
  if (ier /= 0) call exit_MPI(0,'Error allocating temporary array single_par')
  single_par(:) = 0.d0

  NlatNlonNdep = 0

  ! read CIJ, density from mantle files
  do m = 1,Nbndz_read
    Ndep = mtle_bnd_ndep(m)
    do b = 1,Nbnd_read
      Nbnd_numpoints = mtle_bnd_nlat(b) * mtle_bnd_nlon(b) * Ndep

      ! checks if anything to do in this band
      if (Nbnd_numpoints == 0) cycle

      call read_general_mantle_model(single_par,'C11    ', &
                                     lat1(b),lat2(b),dlat(b),dlon(b),dep1(m),dep2(m),ddep(m),NlatNlonNdep)
      mantle_coef(1,NlatNlonNdep+1:NlatNlonNdep+Nbnd_numpoints) = single_par(NlatNlonNdep+1:NlatNlonNdep+Nbnd_numpoints)

      call read_general_mantle_model(single_par,'C13    ', &
                                     lat1(b),lat2(b),dlat(b),dlon(b),dep1(m),dep2(m),ddep(m),NlatNlonNdep)
      mantle_coef(2,NlatNlonNdep+1:NlatNlonNdep+Nbnd_numpoints) = single_par(NlatNlonNdep+1:NlatNlonNdep+Nbnd_numpoints)

      call read_general_mantle_model(single_par,'C33    ', &
                                     lat1(b),lat2(b),dlat(b),dlon(b),dep1(m),dep2(m),ddep(m),NlatNlonNdep)
      mantle_coef(3,NlatNlonNdep+1:NlatNlonNdep+Nbnd_numpoints) = single_par(NlatNlonNdep+1:NlatNlonNdep+Nbnd_numpoints)

      call read_general_mantle_model(single_par,'C44    ', &
                                     lat1(b),lat2(b),dlat(b),dlon(b),dep1(m),dep2(m),ddep(m),NlatNlonNdep)
      mantle_coef(4,NlatNlonNdep+1:NlatNlonNdep+Nbnd_numpoints) = single_par(NlatNlonNdep+1:NlatNlonNdep+Nbnd_numpoints)

      call read_general_mantle_model(single_par,'C66    ', &
                                     lat1(b),lat2(b),dlat(b),dlon(b),dep1(m),dep2(m),ddep(m),NlatNlonNdep)
      mantle_coef(5,NlatNlonNdep+1:NlatNlonNdep+Nbnd_numpoints) = single_par(NlatNlonNdep+1:NlatNlonNdep+Nbnd_numpoints)

      call read_general_mantle_model(mantle_rho(:),'density', &
                                     lat1(b),lat2(b),dlat(b),dlon(b),dep1(m),dep2(m),ddep(m),NlatNlonNdep)

      NlatNlonNdep = NlatNlonNdep + Nbnd_numpoints
    enddo
  enddo ! reading files

  ! frees memory
  deallocate(lat1)
  deallocate(lat2)
  deallocate(dlat)
  deallocate(dlon)
  deallocate(dep1)
  deallocate(dep2)
  deallocate(ddep)
  deallocate(dp_ddep)
  deallocate(single_par)

  ! checks if all read
  if (NlatNlonNdep /= MANTLE_NB * MANTLE_NBZ) then
    print *,'Error: spiral mantle total number of parameters invalid',NlatNlonNdep,' should be ',MANTLE_NB * MANTLE_NBZ
    call exit_MPI(0,'Error spiral invalid total number of parameters')
  endif

  ! debug: file output for smoothed data
  !    open(77,file='tmp-mantle.dat',status='unknown')
  !    write(77,*)'#crustal thickness: #lat (degree) #lon (degree) (w/out ice) #vp (at surface) #vs (at surface)'
  !
  !    ! smoothed version
  !    do k = 1,180
  !      lat = 90.d0 - k + 0.5d0
  !      do i = 1,360
  !        lon = -180.d0 + i - 0.5d0
  !        r = 0.97d0
  !!      do i=1,360
  !!         r=0.30 + (0.99-0.3)*(i-1)/(360-1)
  !!     lat=40.00
  !!      lon=-99.00
  !    write(IMAIN,*) 'k,i,lat,lon',k,i,lat,lon
  !!        call model_mantle_spiral(r,lat,lon,vpv,vph,vsv,vsh,eta,rho, &
  !!                             c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
  !!                              c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
  !!
  !!        write(77,*) r, vph*(R_PLANET_KM*dsqrt(PI*GRAV*RHOAV))
  !        write(77,*)lat,lon, eta, &
  !         vph*(R_PLANET_KM*dsqrt(PI*GRAV*RHOAV)),vsh*(R_PLANET_KM*dsqrt(PI*GRAV*RHOAV)),rho*(RHOAV/1000.0d0)
  !!      enddo
  !    enddo
  !    close(77)

  ! opens internal topography info file
  open(IIN,file='DATA/spiral1.4/mantle/transitionzone_topo.txt',action='read',status='old',iostat=ier)
  if (ier /= 0) then
    write(IMAIN,*) 'Error opening "DATA/spiral1.4/mantle/transitionzone_topo.txt": ', ier
    call exit_MPI(0,'Error model spiral1.4: file not found DATA/spiral1.4/mantle/transitionzone_topo.txt')
  endif

  ! checking dimension
  read(IIN,*) NLA, NLO
  if (NLA /= TOPO_NLA .or. NLO /= TOPO_NLO) then
    write(IMAIN,*) 'Error opening "DATA/spiral1.4/mantle/transitionzone_topo.txt", wrong dimension: ', NLA, NLO
    call exit_MPI(0,'Error model spiral1.4: wrong dimension in DATA/spiral1.4/mantle/transitionzone_topo.txt')
  endif

  k = 0
  do ila = 1,TOPO_NLA
    do ilo = 1,TOPO_NLO
      k = k+1
      read(IIN,*) dummy, dummy, mantle_d410(k),mantle_d660(k)
    enddo
  enddo
  close(IIN)

  ! topography statistics
  ! gets min/max values
  d410_min = HUGEVAL
  d410_max = -HUGEVAL
  d660_min = HUGEVAL
  d660_max = -HUGEVAL

  do k = 1,TOPO_NLA*TOPO_NLO
    if (mantle_d410(k) > d410_max) d410_max = mantle_d410(k)
    if (mantle_d410(k) < d410_min) d410_min = mantle_d410(k)
    if (mantle_d660(k) > d660_max) d660_max = mantle_d660(k)
    if (mantle_d660(k) < d660_min) d660_min = mantle_d660(k)
  enddo

  ! user output
  write(IMAIN,*) '  410 topography min/max = ',sngl(d410_min),sngl(d410_max),' km'
  write(IMAIN,*) '  660 topography min/max = ',sngl(d660_min),sngl(d660_max),' km'
  call flush_IMAIN()

  ! checks min/max
  if (d410_min == HUGEVAL .or. d410_max == -HUGEVAL) stop 'incorrect d410 topography in read_mantle_spiral_model'
  if (d660_min == HUGEVAL .or. d660_max == -HUGEVAL) stop 'incorrect d660 topography in read_mantle_spiral_model'

  end subroutine read_mantle_spiral_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_general_mantle_model(var,var_letter,lat1,lat2,dlat,dlon,dep1,dep2,dep,NlatNlonNdep)

  use constants
  use model_spiral_mantle_par, only: MANTLE_NB, MANTLE_NBZ

  implicit none

  double precision, intent(inout), dimension(MANTLE_NB*MANTLE_NBZ) :: var
  character(len=7), intent(in) :: var_letter
  character(len=7), intent(in) :: lat1,lat2,dlat,dlon,dep1,dep2,dep
  integer, intent(in) :: NlatNlonNdep

  ! local variables
  character(len=MAX_STRING_LEN) :: filemantle, line
  integer :: ier, ila, ilo, ide, Nlat, Nlon, Ndep, k, nlines
  double precision :: dlat1,dlat2,ddlat,ddlon,ddep1,ddep2,ddep

  write(filemantle,'(a39,a,a6,a,a6,a,a6,a,a4,a,a4,a,a4,a,a1,a)') &
    'DATA/spiral1.4/mantle/mantle_band_lat1_', &
    trim(lat1),'_lat2_',trim(lat2),'_dlat_',trim(dlat),'_dlon_',trim(dlon), &
    '_d1_',trim(dep1),'_d2_',trim(dep2),'_dZ_',trim(dep),'.',var_letter

  read(lat1,*) dlat1
  read(lat2,*) dlat2
  read(dlat,*) ddlat
  read(dlon,*) ddlon
  read(dep1,*) ddep1
  read(dep2,*) ddep2
  read(dep,*) ddep

  Nlat = int((dlat2 - dlat1)/ddlat) + 1
  Nlon = int((180.d0 + 180.d0)/ddlon) + 1
  Ndep = int((ddep2 - ddep1)/ddep) + 1

  ! debug
  !print *,'debug: mantle reading', trim(filemantle), ' Nlat, Nlon, Ndep, NlatNlonNdep =',Nlat, Nlon , Ndep, NlatNlonNdep

  open(unit = IIN,file=trim(filemantle),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    write(IMAIN,*) 'Error opening "', trim(filemantle), '": ', ier
    call flush_IMAIN()
    ! stop
    call exit_MPI(0, 'Error mantle model spiral1.4')
  endif

  ! count number of lines
  nlines = 0
  do while(ier == 0)
    read(IIN,*,iostat=ier) line
    if (ier == 0) nlines = nlines + 1
  enddo
  rewind(IIN)

  ! checks line count
  if (nlines /= Nlon * Nlat * Ndep) then
    print *,'Error: file ',trim(filemantle),' has wrong number of lines ',nlines,' should be ',Nlon * Nlat * Ndep
    call exit_MPI(0,'Error spiral mantle file w/ wrong number of lines')
  endif

  k = NlatNlonNdep
  do ide = 1,Ndep
    do ila = 1,Nlat
      do ilo = 1,Nlon
        k = k+1
        read(IIN,*) var(k)
      enddo
    enddo
  enddo
  close(IIN)

  ! debug
  !print *,'debug: reading', trim(filemantle), ' NlatNlonNdep+Nlat*Nlon*Ndep=',k

  end subroutine read_general_mantle_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_mantle_spiral(r,lat,lon,vpv,vph,vsv,vsh,eta,rho, &
                                 c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                 c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  use constants
  use shared_parameters, only: R_PLANET,R_PLANET_KM,RHOAV

  use model_spiral_mantle_par

  implicit none

  double precision,intent(in) :: lat,lon,r ! lat = [-90;90], lon = [-180;180]
  double precision,intent(out) :: vpv,vph,vsv,vsh,eta,rho
  double precision,intent(out) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! local parameters
  double precision :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66
  double precision :: scale_GPa,scaleval,theta,phi
  double precision :: mtle_rho
  double precision,dimension(5):: mtle_coefs

  ! Should be [-90,90], [-180,180] convention
  if (lat > 90.0d0 .or. lat < -90.0d0 .or. lon > 180.0d0 .or. lon < -180.0d0) then
    print *,'Error: in lat/lon:',lat,lon,'in routine model_mantle_spiral()'
    stop 'Error in latitude/longitude range in model_mantle_spiral()'
  endif

  ! initializes
  c11 = ZERO
  c12 = ZERO
  c13 = ZERO
  c14 = ZERO
  c15 = ZERO
  c16 = ZERO
  c22 = ZERO
  c23 = ZERO
  c24 = ZERO
  c25 = ZERO
  c26 = ZERO
  c33 = ZERO
  c34 = ZERO
  c35 = ZERO
  c36 = ZERO
  c44 = ZERO
  c45 = ZERO
  c46 = ZERO
  c55 = ZERO
  c56 = ZERO
  c66 = ZERO

  vsv = ZERO
  vsh = ZERO
  vpv = ZERO
  vph = ZERO
  rho = ZERO
  eta = ZERO

  ! gets mantle coefficients
  call read_mantle_spiral(r,lat,lon,mtle_coefs,mtle_rho)

  ! local parameters initialization
  d11 = ZERO
  d12 = ZERO
  d13 = ZERO
  d14 = ZERO
  d15 = ZERO
  d16 = ZERO
  d22 = ZERO
  d23 = ZERO
  d24 = ZERO
  d25 = ZERO
  d26 = ZERO
  d33 = ZERO
  d34 = ZERO
  d35 = ZERO
  d36 = ZERO
  d44 = ZERO
  d45 = ZERO
  d46 = ZERO
  d55 = ZERO
  d56 = ZERO
  d66 = ZERO

  ! get mantle properties for a given lat/lon/r
  ! in local reference:
  ! d11 = A = rho * vph**2
  ! d13 = F = eta * (A - 2*L)
  ! d33 = C = rho * vpv**2
  ! d44 = L = rho * vsv**2
  ! d66 = N = rho * vsh**2 = (C11-C12)/2

  ! d12 = d11 - 2 d66 = A - 2*N = rho * (vph**2 - 2 * vsh**2)
  ! d22 = d11 = A
  ! d23 = d13 = F
  ! d55 = d44 = L

  d11 = mtle_coefs(1)  ! A
  d13 = mtle_coefs(2)  ! F
  d33 = mtle_coefs(3)  ! C
  d44 = mtle_coefs(4)  ! L
  d66 = mtle_coefs(5)  ! N

  d12 = d11 - 2.d0*d66
  d22 = d11
  d23 = d13
  d55 = d44

  rho = mtle_rho

  ! TI velocities
  vph = sqrt(d11/rho)
  vpv = sqrt(d33/rho)
  vsh = sqrt(d66/rho)
  vsv = sqrt(d44/rho)
  eta = d13/(d11 - 2.d0*d44)

  scaleval = ONE / ( R_PLANET_KM * dsqrt(PI*GRAV*RHOAV) )
  vsv = vsv * scaleval
  vsh = vsh * scaleval
  vpv = vpv * scaleval
  vph = vph * scaleval
  rho = rho * 1000.0d0 / RHOAV

  ! debug
  !print *,'debug: vsv,vsh,vpv,vph,rho,eta'
  !print *,vsv,vsh,vpv,vph,rho,eta

  ! non-dimensionalize the elastic coefficients using
  ! the scale of GPa--[g/cm^3][(km/s)^2]
  scaleval = dsqrt(PI*GRAV*RHOAV)
  scale_GPa = (RHOAV/1000.d0)*((R_PLANET*scaleval/1000.d0)**2)

  d11 = d11/scale_GPa
  d12 = d12/scale_GPa
  d13 = d13/scale_GPa
  d14 = d14/scale_GPa
  d15 = d15/scale_GPa
  d16 = d16/scale_GPa
  d22 = d22/scale_GPa
  d23 = d23/scale_GPa
  d24 = d24/scale_GPa
  d25 = d25/scale_GPa
  d26 = d26/scale_GPa
  d33 = d33/scale_GPa
  d34 = d34/scale_GPa
  d35 = d35/scale_GPa
  d36 = d36/scale_GPa
  d44 = d44/scale_GPa
  d45 = d45/scale_GPa
  d46 = d46/scale_GPa
  d55 = d55/scale_GPa
  d56 = d56/scale_GPa
  d66 = d66/scale_GPa

  ! assign the local (d_ij) or global (c_ij) anisotropic parameters.
  ! The c_ij are the coefficients in the global
  ! reference frame used in SPECFEM3D.
  theta = (90.d0 - lat) * DEGREES_TO_RADIANS
  phi = lon * DEGREES_TO_RADIANS
  call rotate_tensor_radial_to_global(theta,phi, &
                                      d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                                      c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                      c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  end subroutine model_mantle_spiral

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mantle_spiral(r_in,lat_in,lon_in,mtle_coefs,mtle_rho)

! get mantle properties at a given lat/lon/r

  use constants
  use shared_parameters, only: R_PLANET_KM

  use model_spiral_mantle_par

  implicit none

  ! argument variables
  double precision,intent(in) :: lat_in,lon_in,r_in
  double precision,intent(out) :: mtle_rho
  double precision,intent(out) :: mtle_coefs(5)

  ! local variables
  integer :: NlatNlonNdep,Ndep,m,l,index
  double precision :: dlon,dlat,ddep,depth
  double precision :: lat,lon,r

  ! interpolation variables
  double precision :: a,b,c  ! weights
  integer :: rec_read  ! position of the record to read in model.dat (direct access file)
  integer :: i_min, i_max, j_min, j_max, k_min, k_max ! upper and lower bound indices
  double precision :: rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8    ! rho corner values in model.dat
  double precision, dimension(5) :: coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8   ! CIJ corner values in model.dat

  ! gets position (in degrees and non-dimensional r)
  lat = lat_in
  lon = lon_in
  r = r_in

  ! checks latitude/longitude
  if (lat > 90.0d0 .or. lat < -90.0d0 .or. lon > 180.0d0 .or. lon < -180.0d0) then
    print *,'Error in lat/lon:',lat,lon
    stop 'Error in latitude/longitude range in mantle spiral1.4'
  endif

  ! makes sure lat/lon are within range
  if (lat == 90.0d0) lat = 89.9999d0
  if (lat == -90.0d0) lat = -89.9999d0
  if (lon == 180.0d0) lon = 179.9999d0
  if (lon == -180.0d0) lon = -179.9999d0

  ! dimensionalize
  depth = R_PLANET_KM * (R_UNIT_SPHERE - r)

  if (depth > 2891.d0) depth = 2891.d0  ! to prevent error due to numerical accuracy.
  if (depth < 27.d0) depth = 27.d0      ! extend mantle up to surface from 27km up.
  if (depth > 2891.0d0 .or. depth < 27.0d0) then
    print *,'Error in depth:',depth
    stop 'Error in depth range in mantle spiral1.4'
  endif

  ! initializes
  mtle_coefs(:) = ZERO
  mtle_rho = ZERO

  NlatNlonNdep = 0

  if (.not. interpolate_mantle) then
    ! debug
    !print *,'debug: checking we are in the right loop', lon, lat, depth

    ! no interpolation
    do m = 1,Nbndz
      Ndep = mtle_bnd_ndep(m)
      do l = 1,Nbnd
        if (lat >= mtle_bnd_lat1(l) .and. lat <= mtle_bnd_lat2(l) .and. &
            depth >= mtle_bnd_dep1(m) .and. depth <= mtle_bnd_dep2(m)) then

          dlat = (mtle_bnd_lat2(l) - mtle_bnd_lat1(l))/(mtle_bnd_nlat(l)-1)
          dlon = 360.d0/(mtle_bnd_nlon(l)-1)
          ddep = (mtle_bnd_dep2(m) - mtle_bnd_dep1(m))/(mtle_bnd_ndep(m)-1)
          ! debug
          !print *,'debug: checking 1', mtle_bnd_lat1(l),mtle_bnd_dep2(m),lat,lon,depth

          index = NlatNlonNdep + &
                  int( (depth - mtle_bnd_dep1(m))/ddep ) * mtle_bnd_nlat(l)*mtle_bnd_nlon(l) + &
                  int( (lat - mtle_bnd_lat1(l))/dlat ) * mtle_bnd_nlon(l) + &
                  int( (lon + 180.d0)/dlon ) + 1

          ! debug
          !print *,'debug: checking', NlatNlon, lat, lon, index

          mtle_rho = mantle_rho(index)
          mtle_coefs(:) = mantle_coef(:,index)

          ! all done
          return
        endif !if (lat >= mtle_bnd_lat1(l) .and. lat <= mtle_bnd_lat2(l)....

        NlatNlonNdep = NlatNlonNdep + mtle_bnd_nlat(l)*mtle_bnd_nlon(l)*Ndep
      enddo ! do l = 1,Nbnd, with no interpolation
    enddo ! do m = 1,Nbndz, with no interpolation

  else
    ! with interpolation
    do m = 1,Nbndz
      Ndep = mtle_bnd_ndep(m)
      do l = 1,Nbnd
        if (lat >= mtle_bnd_lat1(l) .and. lat <= mtle_bnd_lat2(l) .and. &
            depth >= mtle_bnd_dep1(m) .and. depth <= mtle_bnd_dep2(m)) then

          dlat = (mtle_bnd_lat2(l) - mtle_bnd_lat1(l))/(mtle_bnd_nlat(l)-1)
          dlon = 360.d0/(mtle_bnd_nlon(l)-1)
          ddep = (mtle_bnd_dep2(m) - mtle_bnd_dep1(m))/(mtle_bnd_ndep(m)-1)

          ! locate the corners
          i_min = int( (lon + 180.d0)/dlon ) + 1
          i_max = i_min + 1
          j_min = int( (lat - mtle_bnd_lat1(l))/dlat ) + 1
          j_max = j_min + 1
          k_min = int( (depth - mtle_bnd_dep1(m))/ddep ) + 1
          k_max = k_min + 1

          ! (1) i_min, j_min, k_min
          rec_read = NlatNlonNdep + (k_min-1)*mtle_bnd_nlat(l)*mtle_bnd_nlon(l) + (j_min-1)*mtle_bnd_nlon(l) + (i_min-1) + 1
          ! index bounds
          if (rec_read < 1) rec_read = 1
          if (rec_read > MANTLE_NB*MANTLE_NBZ) rec_read = MANTLE_NB*MANTLE_NBZ

          rho1 = mantle_rho(rec_read)
          coef1(:) = mantle_coef(:,rec_read)

          ! (2) i_min, j_max, k_min
          rec_read = NlatNlonNdep + (k_min-1)*mtle_bnd_nlat(l)*mtle_bnd_nlon(l) + (j_max-1)*mtle_bnd_nlon(l) + (i_min-1) + 1
          ! index bounds
          if (rec_read < 1) rec_read = 1
          if (rec_read > MANTLE_NB*MANTLE_NBZ) rec_read = MANTLE_NB*MANTLE_NBZ

          rho2 = mantle_rho(rec_read)
          coef2(:) = mantle_coef(:,rec_read)

          ! (3) i_max, j_min, k_min
          rec_read = NlatNlonNdep + (k_min-1)*mtle_bnd_nlat(l)*mtle_bnd_nlon(l) + (j_min-1)*mtle_bnd_nlon(l) + (i_max-1) + 1
          ! index bounds
          if (rec_read < 1) rec_read = 1
          if (rec_read > MANTLE_NB*MANTLE_NBZ) rec_read = MANTLE_NB*MANTLE_NBZ

          rho3 = mantle_rho(rec_read)
          coef3(:) = mantle_coef(:,rec_read)

          ! (4) i_max, j_max, k_min
          rec_read = NlatNlonNdep + (k_min-1)*mtle_bnd_nlat(l)*mtle_bnd_nlon(l) + (j_max-1)*mtle_bnd_nlon(l) + (i_max-1) + 1
          ! index bounds
          if (rec_read < 1) rec_read = 1
          if (rec_read > MANTLE_NB*MANTLE_NBZ) rec_read = MANTLE_NB*MANTLE_NBZ

          rho4 = mantle_rho(rec_read)
          coef4(:) = mantle_coef(:,rec_read)

          ! (5) i_min, j_min, k_max
          rec_read = NlatNlonNdep + (k_max-1)*mtle_bnd_nlat(l)*mtle_bnd_nlon(l) + (j_min-1)*mtle_bnd_nlon(l) + (i_min-1) + 1
          ! index bounds
          if (rec_read < 1) rec_read = 1
          if (rec_read > MANTLE_NB*MANTLE_NBZ) rec_read = MANTLE_NB*MANTLE_NBZ

          rho5 = mantle_rho(rec_read)
          coef5(:) = mantle_coef(:,rec_read)

          ! (6) i_min, j_max, k_max
          rec_read = NlatNlonNdep + (k_max-1)*mtle_bnd_nlat(l)*mtle_bnd_nlon(l) + (j_max-1)*mtle_bnd_nlon(l) + (i_min-1) + 1
          ! index bounds
          if (rec_read < 1) rec_read = 1
          if (rec_read > MANTLE_NB*MANTLE_NBZ) rec_read = MANTLE_NB*MANTLE_NBZ

          rho6 = mantle_rho(rec_read)
          coef6(:) = mantle_coef(:,rec_read)

          ! (7) i_max, j_min, k_max
          rec_read = NlatNlonNdep + (k_max-1)*mtle_bnd_nlat(l)*mtle_bnd_nlon(l) + (j_min-1)*mtle_bnd_nlon(l) + (i_max-1) + 1
          ! index bounds
          if (rec_read < 1) rec_read = 1
          if (rec_read > MANTLE_NB*MANTLE_NBZ) rec_read = MANTLE_NB*MANTLE_NBZ

          rho7 = mantle_rho(rec_read)
          coef7(:) = mantle_coef(:,rec_read)

          ! (8) i_max, j_max, k_max
          rec_read = NlatNlonNdep + (k_max-1)*mtle_bnd_nlat(l)*mtle_bnd_nlon(l) + (j_max-1)*mtle_bnd_nlon(l) + (i_max-1) + 1
          ! index bounds
          if (rec_read < 1) rec_read = 1
          if (rec_read > MANTLE_NB*MANTLE_NBZ) rec_read = MANTLE_NB*MANTLE_NBZ

          rho8 = mantle_rho(rec_read)
          coef8(:) = mantle_coef(:,rec_read)

          ! weights
          a = ( lon + 180.d0 - (i_min-1)*dlon ) / dlon
          b = ( lat - mtle_bnd_lat1(l) - (j_min-1)*dlat ) / dlat
          c = ( depth - mtle_bnd_dep1(m) - (k_min-1)*ddep ) / ddep

          mtle_rho = rho1*(1.-a)*(1.-b) + rho2*(1.-a)*b + rho3*a*(1.-b) + rho4*a*b
          mtle_coefs(:) = coef1(:)*(1.-a)*(1.-b) + coef2(:)*(1.-a)*b + coef3(:)*a*(1.-b) + coef4(:)*a*b

          ! all done
          return
        endif

        NlatNlonNdep = NlatNlonNdep + mtle_bnd_nlat(l)*mtle_bnd_nlon(l)*Ndep

      enddo ! do l = 1,Nbnd, with interpolation
    enddo ! do m = 1,Nbndz, with interpolation

  endif ! if (.not. interpolation)

  end subroutine read_mantle_spiral

!
!-------------------------------------------------------------------------------------------------
!

  subroutine add_topography_mantle_spiral(xelm,yelm,zelm)

  use constants
  use shared_parameters, only: R_PLANET,R_PLANET_KM,ELLIPTICITY
  use meshfem_par, only: R220,R400,R670,R771

  implicit none

  double precision, intent(inout) :: xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  ! local parameters
  integer :: ia
  double precision :: topo410,topo660
  double precision :: r,lat,lon
  double precision :: gamma
  double precision :: x,y,z

  ! we loop on all the points of the element
  do ia = 1,NGNOD

    x = xelm(ia)
    y = yelm(ia)
    z = zelm(ia)

    ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
    ! note: at this point, the mesh is still spherical (no need to correct latitude for ellipticity)
    call xyz_2_rlatlon_dble(x,y,z,r,lat,lon,ELLIPTICITY)

    ! The above subroutine produces longitudes with [0,360] convention.
    ! Need to convert to [-180,180] convention expected by subtopo_spiral below
    if (lon > 180.d0) lon = lon - 360.d0

    ! Should be [-90,90], [-180,180] convention
    if (lat > 90.0d0 .or. lat < -90.0d0 .or. lon > 180.0d0 .or. lon < -180.0d0) then
      print *,'Error: in lat/lon:',lat,lon,'in routine add_topography_mantle_spiral()'
      stop 'Error in latitude/longitude range in add_topography_mantle_spiral())'
    endif

    ! stretching occurs between 220 and 770
    if (r > R220/R_PLANET .or. r < R771/R_PLANET) cycle

    ! compute topography on 410 and 660 at current point
    call subtopo_spiral(lat,lon,topo410,topo660)

    ! debug
    !print *,'debug: lat, lon, topo410, topo660',lat,lon,topo410,topo660

    if (topo410 == 0.d0 .and. topo660 == 0.d0) return

    ! non-dimensionalize the topography, which is in km
    ! positive for a depression, so change the sign for a perturbation in radius
    topo410 = -(topo410) / R_PLANET_KM
    topo660 = -(topo660) / R_PLANET_KM

    gamma = 0.d0
    if (r >= R400/R_PLANET .and. r <= R220/R_PLANET) then
      ! stretching between R220 and R400
      gamma = (R220/R_PLANET - r) / (R220/R_PLANET - R400/R_PLANET)
      xelm(ia) = x*(ONE + gamma * topo410 / r)
      yelm(ia) = y*(ONE + gamma * topo410 / r)
      zelm(ia) = z*(ONE + gamma * topo410 / r)
    else if (r >= R771/R_PLANET .and. r <= R670/R_PLANET) then
      ! stretching between R771 and R670
      gamma = (r - R771/R_PLANET) / (R670/R_PLANET - R771/R_PLANET)
      xelm(ia) = x*(ONE + gamma * topo660 / r)
      yelm(ia) = y*(ONE + gamma * topo660 / r)
      zelm(ia) = z*(ONE + gamma * topo660 / r)
    else if (r > R670/R_PLANET .and. r < R400/R_PLANET) then
      ! stretching between R670 and R400
      gamma = (R400/R_PLANET - r) / (R400/R_PLANET - R670/R_PLANET)
      xelm(ia) = x*(ONE + (topo410 + gamma * (topo660 - topo410)) / r)
      yelm(ia) = y*(ONE + (topo410 + gamma * (topo660 - topo410)) / r)
      zelm(ia) = z*(ONE + (topo410 + gamma * (topo660 - topo410)) / r)
    endif
    if (gamma < -0.0001 .or. gamma > 1.0001) call exit_MPI(myrank,'incorrect value of gamma for 410-660 topography')

  enddo

  end subroutine add_topography_mantle_spiral

!
!-------------------------------------------------------------------------------------------------
!
!
! calculate depressions of the 410 and 660km discontinuities in km
! convert depth of 410/660 into preturbations
!

  subroutine subtopo_spiral(lat_in,lon_in,topo410,topo660)

  use constants
  use shared_parameters, only: R_PLANET

  use meshfem_par, only: R400,R670

  use model_spiral_mantle_par

  implicit none

  double precision,intent(in) :: lat_in,lon_in
  double precision,intent(out) :: topo410,topo660

  ! local variables
  integer :: index
  double precision :: dlon,dlat,lat,lon
  ! interpolation variables
  double precision :: a,b  ! weights
  integer :: rec_read  ! position of the record to read in model.dat (direct access file)
  integer :: i_min, i_max, j_min, j_max ! upper and lower bound indices
  double precision :: t410_1,t410_2,t410_3,t410_4   ! t410 corner values in model.dat
  double precision :: t660_1,t660_2,t660_3,t660_4   ! t660 corner values in model.dat

  ! gets position
  lat = lat_in
  lon = lon_in

  ! makes sure lat/lon are within range
  if (lat == 90.0d0) lat = 89.9999d0
  if (lat == -90.0d0) lat = -89.9999d0
  if (lon == 180.0d0) lon = 179.9999d0
  if (lon == -180.0d0) lon = -179.9999d0

  ! initializes
  topo410 = 0.d0
  topo660 = 0.d0

  dlat = 1.d0/TOPO_RES
  dlon = 1.d0/TOPO_RES

  if (.not. interpolate_topo) then
    ! no interpolation
    index = int( (lat + 90.d0)/dlat ) * TOPO_NLO + &
            int( (lon + 180.d0)/dlon ) + 1

    topo410 = mantle_d410(index)
    topo660 = mantle_d660(index)

  else
    ! with interpolation
    ! locate the corners
    i_min = int( (lon + 180.d0)/dlon ) + 1
    i_max = i_min + 1
    j_min = int( (lat + 90.d0)/dlat ) + 1
    j_max = j_min + 1

    ! (1) i_min, j_min
    rec_read = (j_min-1)*TOPO_NLO + (i_min-1) + 1
    ! index bounds
    if (rec_read < 1) rec_read = 1
    if (rec_read > TOPO_NLO*TOPO_NLA) rec_read = TOPO_NLO*TOPO_NLA

    t410_1 = mantle_d410(rec_read)
    t660_1 = mantle_d660(rec_read)

    ! (2) i_min, j_max
    rec_read = (j_max-1)*TOPO_NLO + (i_min-1) + 1
    ! index bounds
    if (rec_read < 1) rec_read = 1
    if (rec_read > TOPO_NLO*TOPO_NLA) rec_read = TOPO_NLO*TOPO_NLA

    t410_2 = mantle_d410(rec_read)
    t660_2 = mantle_d660(rec_read)

    ! (3) i_max, j_min
    rec_read = (j_min-1)*TOPO_NLO + (i_max-1) + 1
    ! index bounds
    if (rec_read < 1) rec_read = 1
    if (rec_read > TOPO_NLO*TOPO_NLA) rec_read = TOPO_NLO*TOPO_NLA

    t410_3 = mantle_d410(rec_read)
    t660_3 = mantle_d660(rec_read)

    ! (4) i_max, j_max
    rec_read = (j_max-1)*TOPO_NLO + (i_max-1) + 1
    ! index bounds
    if (rec_read < 1) rec_read = 1
    if (rec_read > TOPO_NLO*TOPO_NLA) rec_read = TOPO_NLO*TOPO_NLA

    t410_4 = mantle_d410(rec_read)
    t660_4 = mantle_d660(rec_read)

    ! weights
    a = ( lon + 180.d0 - (i_min-1)*dlon ) / dlon
    b = ( lat + 90.d0 - (j_min-1)*dlat ) / dlat

    topo410 = t410_1*(1.d0-a)*(1.d0-b) + t410_2*(1.d0-a)*b + t410_3*a*(1.d0-b) + t410_4*a*b
    topo660 = t660_1*(1.d0-a)*(1.d0-b) + t660_2*(1.d0-a)*b + t660_3*a*(1.d0-b) + t660_4*a*b

  endif ! if (.not. interpolation)

  topo410 = topo410 - (R_PLANET-R400)/1000.d0
  topo660 = topo660 - (R_PLANET-R670)/1000.d0

  end subroutine subtopo_spiral
