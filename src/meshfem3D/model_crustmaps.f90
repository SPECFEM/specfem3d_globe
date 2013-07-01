!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            April 2011
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

!--------------------------------------------------------------------------------------------------
! General Crustmaps
!
! combines Crust2.0 and EUcrust07 for moho depths; the crustal maps are
! interpolating the crustal velocities from Crust2.0 onto the more detailed EUcrust
! crustal depths where ever they are defined.

! current crustmaps (cmaps) take sediment thickness
! and moho depths from EUcrust07 if possible and interpolate corresponding
! velocity/densities given from Crust2.0.
!
! main author: Matthias Meschede (meschede@princeton.edu)
!--------------------------------------------------------------------------------------------------

  subroutine model_crustmaps_broadcast(myrank,GC_V)

! standard routine to setup model

  implicit none

  include "constants.h"
  ! standard include of the MPI library
  include 'mpif.h'

  integer :: myrank

  !model_crustmaps_variables
  type model_crustmaps_variables
    sequence
    double precision, dimension(180*CRUSTMAP_RESOLUTION,&
      360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: thickness
    double precision, dimension(180*CRUSTMAP_RESOLUTION, &
      360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: density
    double precision, dimension(180*CRUSTMAP_RESOLUTION, &
      360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: velocp
    double precision, dimension(180*CRUSTMAP_RESOLUTION, &
      360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: velocs

    double precision thicknessnp(NLAYERS_CRUSTMAP)
    double precision densitynp(NLAYERS_CRUSTMAP)
    double precision velocpnp(NLAYERS_CRUSTMAP)
    double precision velocsnp(NLAYERS_CRUSTMAP)
    double precision thicknesssp(NLAYERS_CRUSTMAP)
    double precision densitysp(NLAYERS_CRUSTMAP)
    double precision velocpsp(NLAYERS_CRUSTMAP)
    double precision velocssp(NLAYERS_CRUSTMAP)

  end type model_crustmaps_variables
  type (model_crustmaps_variables) GC_V
  !model_crustmaps_variables

  ! local parameters
  integer :: ier

  ! master reads in crust maps
  if(myrank == 0) &
    call read_general_crustmap(GC_V)

  ! broadcasts values to all processes
  call MPI_BCAST(GC_V%thickness,180*360*CRUSTMAP_RESOLUTION*CRUSTMAP_RESOLUTION*NLAYERS_CRUSTMAP, &
    MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(GC_V%velocp,180*360*CRUSTMAP_RESOLUTION*CRUSTMAP_RESOLUTION*NLAYERS_CRUSTMAP, &
    MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(GC_V%velocs,180*360*CRUSTMAP_RESOLUTION*CRUSTMAP_RESOLUTION*NLAYERS_CRUSTMAP, &
    MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(GC_V%density,180*360*CRUSTMAP_RESOLUTION*CRUSTMAP_RESOLUTION*NLAYERS_CRUSTMAP, &
    MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

  ! north pole
  call MPI_BCAST(GC_V%thicknessnp,NLAYERS_CRUSTMAP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(GC_V%densitynp,NLAYERS_CRUSTMAP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(GC_V%velocpnp,NLAYERS_CRUSTMAP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(GC_V%velocsnp,NLAYERS_CRUSTMAP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(GC_V%densitynp,NLAYERS_CRUSTMAP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

  ! south pole
  call MPI_BCAST(GC_V%thicknesssp,NLAYERS_CRUSTMAP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(GC_V%densitysp,NLAYERS_CRUSTMAP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(GC_V%velocpsp,NLAYERS_CRUSTMAP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(GC_V%velocssp,NLAYERS_CRUSTMAP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(GC_V%densitysp,NLAYERS_CRUSTMAP,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)


  end subroutine model_crustmaps_broadcast

!
!-------------------------------------------------------------------------------------------------
!

! read general crustmap by Matthias Meschede

  subroutine read_general_crustmap(GC_V)

  implicit none
  include "constants.h"

! Matthias Meschede
 !model_crustmaps_variables
  type model_crustmaps_variables
    sequence
    double precision, dimension(180*CRUSTMAP_RESOLUTION,&
      360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: thickness
    double precision, dimension(180*CRUSTMAP_RESOLUTION, &
      360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: density
    double precision, dimension(180*CRUSTMAP_RESOLUTION, &
      360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: velocp
    double precision, dimension(180*CRUSTMAP_RESOLUTION, &
      360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: velocs

    double precision thicknessnp(NLAYERS_CRUSTMAP)
    double precision densitynp(NLAYERS_CRUSTMAP)
    double precision velocpnp(NLAYERS_CRUSTMAP)
    double precision velocsnp(NLAYERS_CRUSTMAP)
    double precision thicknesssp(NLAYERS_CRUSTMAP)
    double precision densitysp(NLAYERS_CRUSTMAP)
    double precision velocpsp(NLAYERS_CRUSTMAP)
    double precision velocssp(NLAYERS_CRUSTMAP)

  end type model_crustmaps_variables
  type (model_crustmaps_variables) GC_V
  !model_crustmaps_variables

  integer i,l

  do i=3,7
    l = i - 2
    call read_general_crustmap_layer(GC_V%thickness(:,:,l), 't', i)
    call read_general_crustmap_layer(GC_V%density(:,:,l),   'r', i)
    call read_general_crustmap_layer(GC_V%velocp(:,:,l),    'p', i)
    call read_general_crustmap_layer(GC_V%velocs(:,:,l),    's', i)
  enddo

  GC_V%thicknessnp(:) = 0.0
  GC_V%thicknesssp(:) = 0.0
  GC_V%densitynp(:) = 0.0
  GC_V%densitysp(:) = 0.0
  GC_V%velocpnp(:) = 0.0
  GC_V%velocpsp(:) = 0.0
  GC_V%velocsnp(:) = 0.0
  GC_V%velocssp(:) = 0.0

  !compute average values for north and southpole
  do l=1,NLAYERS_CRUSTMAP
    do i=1,360*CRUSTMAP_RESOLUTION
      GC_V%thicknessnp(l) =  GC_V%thicknessnp(l)+GC_V%thickness(1,i,l)
      GC_V%thicknesssp(l) = GC_V%thicknesssp(l)+GC_V%thickness(180*CRUSTMAP_RESOLUTION,i,l)
      GC_V%densitynp(l) = GC_V%densitynp(l)+GC_V%density(1,i,l)
      GC_V%densitysp(l) = GC_V%densitysp(l)+GC_V%density(180*CRUSTMAP_RESOLUTION,i,l)
      GC_V%velocpnp(l) = GC_V%velocpnp(l)+GC_V%velocp(1,i,l)
      GC_V%velocpsp(l) = GC_V%velocpsp(l)+GC_V%velocp(180*CRUSTMAP_RESOLUTION,i,l)
      GC_V%velocsnp(l) = GC_V%velocsnp(l)+GC_V%velocs(1,i,l)
      GC_V%velocssp(l) = GC_V%velocssp(l)+GC_V%velocs(180*CRUSTMAP_RESOLUTION,i,l)
    enddo
    GC_V%thicknessnp(l) = GC_V%thicknessnp(l)/360.0/dble(CRUSTMAP_RESOLUTION)
    GC_V%thicknesssp(l) = GC_V%thicknesssp(l)/360.0/dble(CRUSTMAP_RESOLUTION)
    GC_V%densitynp(l) = GC_V%densitynp(l)/360.0/dble(CRUSTMAP_RESOLUTION)
    GC_V%densitysp(l) = GC_V%densitysp(l)/360.0/dble(CRUSTMAP_RESOLUTION)
    GC_V%velocpnp(l) = GC_V%velocpnp(l)/360.0/dble(CRUSTMAP_RESOLUTION)
    GC_V%velocpsp(l) = GC_V%velocpsp(l)/360.0/dble(CRUSTMAP_RESOLUTION)
    GC_V%velocsnp(l) = GC_V%velocsnp(l)/360.0/dble(CRUSTMAP_RESOLUTION)
    GC_V%velocssp(l) = GC_V%velocssp(l)/360.0/dble(CRUSTMAP_RESOLUTION)

!    print *,'thicknessnp(',l,')',GC_V%thicknessnp(l)
  enddo


  end subroutine read_general_crustmap

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_general_crustmap_layer(var,var_letter,ind)

  implicit none
  include "constants.h"

  double precision, intent(out), &
    dimension(180*CRUSTMAP_RESOLUTION,360*CRUSTMAP_RESOLUTION)&
    :: var
  character(len=1), intent(in) :: var_letter
  integer, intent(in) :: ind

  ! local variables
  character(len=50) :: config_name
  character(len=150) :: default_name
  character(len=150) :: eucrust
  integer :: ier, ila, iln

  write(config_name,'(a,a1,i1)') 'model.eucrust', var_letter, ind
  write(default_name,'(a,a1,i1)') 'DATA/crustmap/eucrust', var_letter, ind

  call get_value_string(eucrust, config_name, default_name)

  open(unit=1,file=eucrust,status='old',action='read')
  if ( ier /= 0 ) then
    write(IMAIN,*) 'error opening "', trim(eucrust), '": ', ier
    call exit_MPI(0, 'error model crustmap')
  endif

  do ila=1,180*CRUSTMAP_RESOLUTION
    read(1,*) (var(ila,iln),iln=1,360*CRUSTMAP_RESOLUTION)
  enddo
  close(1)

  end subroutine read_general_crustmap_layer

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_crustmaps(lat,lon,x,vp,vs,rho,moho,found_crust,GC_V,elem_in_crust)

! Matthias Meschede
! read smooth crust2.0 model (0.25 degree resolution) with eucrust
! based on software routines provided with the crust2.0 model by Bassin et al.
!

  implicit none
  include "constants.h"

!Matthias Meschede
 !model_crustmaps_variables
  type model_crustmaps_variables
    sequence
    double precision, dimension(180*CRUSTMAP_RESOLUTION,&
      360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: thickness
    double precision, dimension(180*CRUSTMAP_RESOLUTION, &
      360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: density
    double precision, dimension(180*CRUSTMAP_RESOLUTION, &
      360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: velocp
    double precision, dimension(180*CRUSTMAP_RESOLUTION, &
      360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: velocs

    double precision thicknessnp(NLAYERS_CRUSTMAP)
    double precision densitynp(NLAYERS_CRUSTMAP)
    double precision velocpnp(NLAYERS_CRUSTMAP)
    double precision velocsnp(NLAYERS_CRUSTMAP)
    double precision thicknesssp(NLAYERS_CRUSTMAP)
    double precision densitysp(NLAYERS_CRUSTMAP)
    double precision velocpsp(NLAYERS_CRUSTMAP)
    double precision velocssp(NLAYERS_CRUSTMAP)

  end type model_crustmaps_variables
  type (model_crustmaps_variables) GC_V
  !model_crustmaps_variables


  double precision :: lat,lon,x,vp,vs,rho,moho
  logical :: found_crust,elem_in_crust
  double precision :: h_sed,h_uc
  double precision :: x3,x4,x5,x6,x7,scaleval
  double precision,dimension(NLAYERS_CRUSTMAP) :: vps,vss,rhos,thicks

  call read_crustmaps(lat,lon,vps,vss,rhos,thicks,GC_V)

  x3 = (R_EARTH-thicks(1)*1000.0d0)/R_EARTH
  h_sed = thicks(1) + thicks(2)
  x4 = (R_EARTH-h_sed*1000.0d0)/R_EARTH
  h_uc = h_sed + thicks(3)
  x5 = (R_EARTH-h_uc*1000.0d0)/R_EARTH
  x6 = (R_EARTH-(h_uc+thicks(4))*1000.0d0)/R_EARTH
  x7 = (R_EARTH-(h_uc+thicks(4)+thicks(5))*1000.0d0)/R_EARTH

  found_crust = .true.
!  if(x > x3 .and. INCLUDE_SEDIMENTS_CRUST &
!   .and. h_sed > MINIMUM_SEDIMENT_THICKNESS) then
  if(x > x3 .and. INCLUDE_SEDIMENTS_CRUST ) then
   vp = vps(1)
   vs = vss(1)
   rho = rhos(1)
!  else if(x > x4 .and. INCLUDE_SEDIMENTS_CRUST &
!   .and. h_sed > MINIMUM_SEDIMENT_THICKNESS) then
  else if(x > x4 .and. INCLUDE_SEDIMENTS_CRUST ) then
   vp = vps(2)
   vs = vss(2)
   rho = rhos(2)
  else if(x > x5) then
   vp = vps(3)
   vs = vss(3)
   rho = rhos(3)
  else if(x > x6) then
   vp = vps(4)
   vs = vss(4)
   rho = rhos(4)
  else if(x > x7 .or. elem_in_crust) then
   vp = vps(5)
   vs = vss(5)
   rho = rhos(5)
  else
   found_crust = .false.
  endif

  if (found_crust) then
  !   non-dimensionalize
    scaleval = dsqrt(PI*GRAV*RHOAV)
    vp = vp*1000.0d0/(R_EARTH*scaleval)
    vs = vs*1000.0d0/(R_EARTH*scaleval)
    rho = rho*1000.0d0/RHOAV
   ! moho = (h_uc+thicks(4)+thicks(5))*1000.0d0/R_EARTH
  else
    scaleval = dsqrt(PI*GRAV*RHOAV)
    vp = 20.0*1000.0d0/(R_EARTH*scaleval)
    vs = 20.0*1000.0d0/(R_EARTH*scaleval)
    rho = 20.0*1000.0d0/RHOAV
  endif

  moho = (h_uc+thicks(4)+thicks(5))*1000.0d0/R_EARTH

  end subroutine model_crustmaps

!
!-------------------------------------------------------------------------------------------------
!


  subroutine read_crustmaps(lat,lon,velp,vels,rhos,thicks,GC_V)

! crustal vp and vs in km/s, layer thickness in km

  implicit none

  include "constants.h"


  ! argument variables
  double precision lat,lon
  double precision rhos(5),thicks(5),velp(5),vels(5)

! Matthias Meschede
 !model_crustmaps_variables
  type model_crustmaps_variables
    sequence
    double precision, dimension(180*CRUSTMAP_RESOLUTION,&
      360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: thickness
    double precision, dimension(180*CRUSTMAP_RESOLUTION, &
      360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: density
    double precision, dimension(180*CRUSTMAP_RESOLUTION, &
      360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: velocp
    double precision, dimension(180*CRUSTMAP_RESOLUTION, &
      360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: velocs

    double precision thicknessnp(NLAYERS_CRUSTMAP)
    double precision densitynp(NLAYERS_CRUSTMAP)
    double precision velocpnp(NLAYERS_CRUSTMAP)
    double precision velocsnp(NLAYERS_CRUSTMAP)
    double precision thicknesssp(NLAYERS_CRUSTMAP)
    double precision densitysp(NLAYERS_CRUSTMAP)
    double precision velocpsp(NLAYERS_CRUSTMAP)
    double precision velocssp(NLAYERS_CRUSTMAP)

  end type model_crustmaps_variables
  type (model_crustmaps_variables) GC_V
  !model_crustmaps_variables

  !-------------------------------
  ! work-around to avoid jacobian problems when stretching mesh elements;
  ! one could also try to slightly change the shape of the doulbing element bricks (which cause the problem)...
  !
  ! defines a "critical" region to have at least a 1-degree smoothing;
  ! critical region can lead to negative jacobians for mesh stretching when CAP smoothing is too small
  double precision,parameter :: LAT_CRITICAL_EUROPE = 50.0d0
  double precision,parameter :: LON_CRITICAL_EUROPE = 22.0d0
  double precision,parameter :: CRITICAL_RANGE_EUROPE = 50.0d0

  ! defines a "critical" region around the andes to have at least a 1-degree smoothing;
  ! critical region can lead to negative jacobians for mesh stretching when CAP smoothing is too small
  double precision,parameter :: LAT_CRITICAL_ANDES = -20.0d0
  double precision,parameter :: LON_CRITICAL_ANDES = -70.0d0
  double precision,parameter :: CRITICAL_RANGE_ANDES = 70.0d0

  ! sampling rate for CAP points
  integer, parameter :: NTHETA = 4
  integer, parameter :: NPHI = 20
  !-------------------------------

! local variables
  double precision weightup,weightleft,weightul,weightur,weightll,weightlr
  double precision xlon(NTHETA*NPHI),xlat(NTHETA*NPHI),weight(NTHETA*NPHI)
  double precision rhol(NLAYERS_CRUSTMAP),thickl(NLAYERS_CRUSTMAP), &
    velpl(NLAYERS_CRUSTMAP),velsl(NLAYERS_CRUSTMAP)
  double precision weightl,cap_degree,dist
  double precision h_sed
  integer num_points
  integer i,ipoin,iupcolat,ileftlng,irightlng

! get integer colatitude and longitude of crustal cap
! -90<lat<90 -180<lon<180
  if(lat > 90.0d0 .or. lat < -90.0d0 .or. lon > 180.0d0 .or. lon < -180.0d0) &
    write(*,*) lat,' ',lon, ' error in latitude/longitude range in crust'
  if(lat==90.0d0) lat=89.9999d0
  if(lat==-90.0d0) lat=-89.9999d0
  if(lon==180.0d0) lon=179.9999d0
  if(lon==-180.0d0) lon=-179.9999d0

  ! by defaults uses only 1 point location
  num_points = 1

  ! checks if inside/outside of critical region for mesh stretching
  if( SMOOTH_CRUST ) then
    dist = dsqrt( (lon-LAT_CRITICAL_EUROPE)**2 + (lat-LAT_CRITICAL_EUROPE )**2 )
    if( dist < CRITICAL_RANGE_EUROPE ) then
      ! sets up smoothing points
      ! by default uses CAP smoothing with crustmap resolution, e.g. 1/4 degree
      cap_degree = 1.d0 / CRUSTMAP_RESOLUTION

      ! increases cap smoothing degree
      ! scales between -1 at center and 0 at border
      dist = dist / CRITICAL_RANGE_EUROPE - 1.0d0
      ! shifts value to 1 at center and 0 to the border with exponential decay
      dist = 1.0d0 - exp( - dist*dist*10.0d0 )
      ! increases smoothing degree inside of critical region
      cap_degree = cap_degree + dist

      ! gets smoothing points and weights
      call CAP_vardegree(lon,lat,xlon,xlat,weight,cap_degree,NTHETA,NPHI)
      num_points = NTHETA*NPHI
    endif
    dist = dsqrt( (lon-LON_CRITICAL_ANDES)**2 + (lat-LAT_CRITICAL_ANDES )**2 )
    if( dist < CRITICAL_RANGE_ANDES ) then
      ! sets up smoothing points
      ! by default uses CAP smoothing with crustmap resolution, e.g. 1/4 degree
      cap_degree = 1.d0 / CRUSTMAP_RESOLUTION

      ! increases cap smoothing degree
      ! scales between -1 at center and 0 at border
      dist = dist / CRITICAL_RANGE_ANDES - 1.0d0
      ! shifts value to 1 at center and 0 to the border with exponential decay
      dist = 1.0d0 - exp( - dist*dist*10.0d0 )
      ! increases smoothing degree inside of critical region
      cap_degree = cap_degree + dist

      ! gets smoothing points and weights
      call CAP_vardegree(lon,lat,xlon,xlat,weight,cap_degree,NTHETA,NPHI)
      num_points = NTHETA*NPHI
    endif
  endif

  ! initializes
  velp(:) = 0.0d0
  vels(:) = 0.0d0
  rhos(:) = 0.0d0
  thicks(:) = 0.0d0

  ! loops over weight points
  do ipoin=1,num_points
    ! checks if more than one weighting points are taken
    if( num_points > 1 ) then
      lat = xlat(ipoin)
      lon = xlon(ipoin)
      ! weighting value
      weightl = weight(ipoin)
    else
      weightl = 1.0d0
    endif

    ! gets crust value indices
    call ibilinearmap(lat,lon,iupcolat,ileftlng,weightup,weightleft)

    ! interpolates location and crust values
    if(iupcolat==0) then
       weightup=weightup*2
    else if(iupcolat==180*CRUSTMAP_RESOLUTION) then
       weightup=2*weightup-1
    endif

    if(ileftlng==360*CRUSTMAP_RESOLUTION) then
      irightlng=1
    else
      irightlng=ileftlng+1
    endif

    weightul=weightup*weightleft
    weightur=weightup*(1.0-weightleft)
    weightll=(1.0-weightup)*weightleft
    weightlr=(1.0-weightup)*(1.0-weightleft)

    if(iupcolat==0) then
      ! north pole
      do i=1,NLAYERS_CRUSTMAP
       thickl(i)=weightul*GC_V%thicknessnp(i)+weightur*GC_V%thicknessnp(i)+&
                 weightll*GC_V%thickness(1,ileftlng,i)+weightlr*GC_V%thickness(1,irightlng,i)

       rhol(i)=weightul*GC_V%densitynp(i)+weightur*GC_V%densitynp(i)+&
               weightll*GC_V%density(1,ileftlng,i)+weightlr*GC_V%density(1,irightlng,i)
       velpl(i)=weightul*GC_V%velocpnp(i)+weightur*GC_V%velocpnp(i)+&
               weightll*GC_V%velocp(1,ileftlng,i)+weightlr*GC_V%velocp(1,irightlng,i)
       velsl(i)=weightul*GC_V%velocsnp(i)+weightur*GC_V%velocsnp(i)+&
               weightll*GC_V%velocs(1,ileftlng,i)+weightlr*GC_V%velocs(1,irightlng,i)
      enddo
    else if(iupcolat==180*CRUSTMAP_RESOLUTION) then
      ! south pole
      do i=1,NLAYERS_CRUSTMAP
       thickl(i)=weightul*GC_V%thickness(iupcolat,ileftlng,i)+weightur*GC_V%thickness(iupcolat,irightlng,i)+&
                 weightll*GC_V%thicknesssp(i)+weightlr*GC_V%thicknesssp(i)
       rhol(i)=weightul*GC_V%density(iupcolat,ileftlng,i)+weightur*GC_V%density(iupcolat,irightlng,i)+&
               weightll*GC_V%densitysp(i)+weightlr*GC_V%densitysp(i)
       velpl(i)=weightul*GC_V%velocp(iupcolat,ileftlng,i)+weightur*GC_V%velocp(iupcolat,irightlng,i)+&
               weightll*GC_V%velocpsp(i)+weightlr*GC_V%velocpsp(i)
       velsl(i)=weightul*GC_V%velocs(iupcolat,ileftlng,i)+weightur*GC_V%velocs(iupcolat,irightlng,i)+&
               weightll*GC_V%velocssp(i)+weightlr*GC_V%velocssp(i)
      enddo
    else
      do i=1,NLAYERS_CRUSTMAP
       thickl(i)=weightul*GC_V%thickness(iupcolat,ileftlng,i)+weightur*GC_V%thickness(iupcolat,irightlng,i)+&
                 weightll*GC_V%thickness(iupcolat+1,ileftlng,i)+weightlr*GC_V%thickness(iupcolat+1,irightlng,i)
       rhol(i)=weightul*GC_V%density(iupcolat,ileftlng,i)+weightur*GC_V%density(iupcolat,irightlng,i)+&
               weightll*GC_V%density(iupcolat+1,ileftlng,i)+weightlr*GC_V%density(iupcolat+1,irightlng,i)
       velpl(i)=weightul*GC_V%velocp(iupcolat,ileftlng,i)+weightur*GC_V%velocp(iupcolat,irightlng,i)+&
               weightll*GC_V%velocp(iupcolat+1,ileftlng,i)+weightlr*GC_V%velocp(iupcolat+1,irightlng,i)
       velsl(i)=weightul*GC_V%velocs(iupcolat,ileftlng,i)+weightur*GC_V%velocs(iupcolat,irightlng,i)+&
               weightll*GC_V%velocs(iupcolat+1,ileftlng,i)+weightlr*GC_V%velocs(iupcolat+1,irightlng,i)
    !   thicks(i)=1.0
    !   rhos(i)=1.0
    !   velp(i)=1.0
    !   vels(i)=1.0i
      enddo
    endif

    ! sediment thickness
    h_sed = thickl(1) + thickl(2)

    ! takes upper crust value if sediment too thin
    if( h_sed < MINIMUM_SEDIMENT_THICKNESS ) then
      velpl(1) = velpl(3)
      velpl(2) = velpl(3)
      velsl(1) = velsl(3)
      velsl(2) = velsl(3)
      rhol(1) = rhol(3)
      rhol(2) = rhol(3)
    endif

    ! total, smoothed values
    rhos(:) = rhos(:) + weightl*rhol(:)
    thicks(:) = thicks(:) + weightl*thickl(:)
    velp(:) = velp(:) + weightl*velpl(:)
    vels(:) = vels(:) + weightl*velsl(:)
  enddo

  end subroutine read_crustmaps

!--------------------------------------------------------------------------------------------

  subroutine ibilinearmap(lat,lng,iupcolat,ileftlng,weightup,weightleft)

  implicit none

  include "constants.h"

  ! argument variables
  double precision weightup,weightleft
  double precision lat,lng, xlng
  double precision buffer
  integer iupcolat
  integer ileftlng

  if(lat > 90.0d0 .or. lat < -90.0d0 .or. lng > 180.0d0 .or. lng < -180.0d0) &
    stop 'error in latitude/longitude range in icolat_ilon'

  if(lng<0) then
    xlng=lng+360.0
  else
    xlng=lng
  endif

  buffer=0.5+((90.0-lat)*CRUSTMAP_RESOLUTION)
  iupcolat=int(buffer)
  weightup=1.0-(buffer-dble(iupcolat))

  if(iupcolat<0) iupcolat=0
  if(iupcolat>180*CRUSTMAP_RESOLUTION)  iupcolat=180*CRUSTMAP_RESOLUTION


  buffer=0.5+(xlng*CRUSTMAP_RESOLUTION)
  ileftlng=int(buffer)
  weightleft=1.0-(buffer-dble(ileftlng))

  if(ileftlng<1) ileftlng=360*CRUSTMAP_RESOLUTION
  if(ileftlng>360*CRUSTMAP_RESOLUTION) ileftlng=1

  end subroutine ibilinearmap
