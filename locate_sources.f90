!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            December 2010
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

!----
!----  locate_sources finds the correct position of the sources
!----

  subroutine locate_sources(NSOURCES,myrank,nspec,nglob,ibool,&
                 xstore,ystore,zstore,xigll,yigll,zigll, &
                 NPROCTOT,ELLIPTICITY,TOPOGRAPHY, &
                 sec,t_cmt,yr,jda,ho,mi,theta_source,phi_source, &
                 NSTEP,DT,hdur,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                 islice_selected_source,ispec_selected_source, &
                 xi_source,eta_source,gamma_source, nu_source, &
                 rspl,espl,espl2,nspl,ibathy_topo,NEX_XI,PRINT_SOURCE_TIME_FUNCTION)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer NPROCTOT
  integer NSTEP,NSOURCES,NEX_XI

  logical ELLIPTICITY,TOPOGRAPHY,PRINT_SOURCE_TIME_FUNCTION

  double precision DT

! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

! for ellipticity
  integer nspl
  double precision rspl(NR),espl(NR),espl2(NR)

  integer nspec,nglob,myrank,isource

  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

! arrays containing coordinates of the points
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore,ystore,zstore

! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX),yigll(NGLLY),zigll(NGLLZ)

  double precision nu_source(NDIM,NDIM,NSOURCES)

  integer yr,jda,ho,mi

  double precision sec
  double precision t_cmt(NSOURCES)
  double precision t0, hdur_gaussian(NSOURCES)

  integer iprocloop

  integer i,j,k,ispec,iglob
  integer ier

  double precision ell
  double precision elevation
  double precision r0,dcost,p20
  double precision theta,phi
  double precision, dimension(NSOURCES) :: theta_source,phi_source
  double precision dist,typical_size
  double precision xi,eta,gamma,dx,dy,dz,dxi,deta

! topology of the control points of the surface element
  integer iax,iay,iaz
  integer iaddx(NGNOD),iaddy(NGNOD),iaddr(NGNOD)

! coordinates of the control points of the surface element
  double precision xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  integer iter_loop

  integer ia
  double precision x,y,z
  double precision xix,xiy,xiz
  double precision etax,etay,etaz
  double precision gammax,gammay,gammaz
  double precision dgamma

  double precision final_distance_source(NSOURCES)
  double precision, dimension(:), allocatable :: final_distance_source_subset

  double precision x_target_source,y_target_source,z_target_source
  double precision r_target_source

  integer islice_selected_source(NSOURCES)

! timer MPI
  double precision time_start,tCPU

  integer isources_already_done,isource_in_this_subset
  integer ispec_selected_source(NSOURCES)
  integer, dimension(:), allocatable :: ispec_selected_source_subset

  integer, dimension(:,:), allocatable :: ispec_selected_source_all
  double precision, dimension(:,:), allocatable :: xi_source_all,eta_source_all,gamma_source_all, &
     final_distance_source_all,x_found_source_all,y_found_source_all,z_found_source_all

  double precision hdur(NSOURCES)

  double precision, dimension(NSOURCES) :: Mxx,Myy,Mzz,Mxy,Mxz,Myz
  double precision, dimension(NSOURCES) :: xi_source,eta_source,gamma_source
  double precision, dimension(:), allocatable :: xi_source_subset,eta_source_subset,gamma_source_subset

  double precision, dimension(NSOURCES) :: lat,long,depth
  double precision scalar_moment
  double precision moment_tensor(6,NSOURCES)
  double precision radius

  character(len=150) OUTPUT_FILES,plot_file

  double precision, dimension(:), allocatable :: x_found_source,y_found_source,z_found_source
  double precision r_found_source
  double precision st,ct,sp,cp
  double precision Mrr,Mtt,Mpp,Mrt,Mrp,Mtp
  double precision colat_source
  double precision distmin

  integer :: ix_initial_guess_source,iy_initial_guess_source,iz_initial_guess_source,NSOURCES_SUBSET_current_size

  logical located_target

! for calculation of source time function and spectrum
  integer it,iom
  double precision time_source,om
  double precision, external :: comp_source_time_function,comp_source_spectrum

! number of points to plot the source time function and spectrum
  integer, parameter :: NSAMP_PLOT_SOURCE = 1000

  integer iorientation
  double precision stazi,stdip,thetan,phin,n(3)

! **************

! make sure we clean the future final array
  ispec_selected_source(:) = 0

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

! read all the sources
  if(myrank == 0) call get_cmt(yr,jda,ho,mi,sec,t_cmt,hdur,lat,long,depth,moment_tensor,DT,NSOURCES)
! broadcast the information read on the master to the nodes
  call MPI_BCAST(yr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(jda,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(ho,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(mi,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  call MPI_BCAST(sec,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

  call MPI_BCAST(t_cmt,NSOURCES,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(hdur,NSOURCES,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(lat,NSOURCES,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(long,NSOURCES,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(depth,NSOURCES,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

  call MPI_BCAST(moment_tensor,6*NSOURCES,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

! define topology of the control element
  call hex_nodes(iaddx,iaddy,iaddr)

! get MPI starting time for all sources
  time_start = MPI_WTIME()

! convert the half duration for triangle STF to the one for gaussian STF
  hdur_gaussian = hdur/SOURCE_DECAY_MIMIC_TRIANGLE

! define t0 as the earliest start time
  t0 = - 1.5d0*minval(t_cmt-hdur)

! loop on all the sources
! gather source information in subsets to reduce memory requirements

! loop over subsets of sources
  do isources_already_done = 0, NSOURCES, NSOURCES_SUBSET_MAX

! the size of the subset can be the maximum size, or less (if we are in the last subset,
! or if there are fewer sources than the maximum size of a subset)
  NSOURCES_SUBSET_current_size = min(NSOURCES_SUBSET_MAX, NSOURCES - isources_already_done)

! allocate arrays specific to each subset
  allocate(final_distance_source_subset(NSOURCES_SUBSET_current_size))

  allocate(ispec_selected_source_subset(NSOURCES_SUBSET_current_size))

  allocate(ispec_selected_source_all(NSOURCES_SUBSET_current_size,0:NPROCTOT-1))

  allocate(xi_source_all(NSOURCES_SUBSET_current_size,0:NPROCTOT-1))
  allocate(eta_source_all(NSOURCES_SUBSET_current_size,0:NPROCTOT-1))
  allocate(gamma_source_all(NSOURCES_SUBSET_current_size,0:NPROCTOT-1))

  allocate(final_distance_source_all(NSOURCES_SUBSET_current_size,0:NPROCTOT-1))

  allocate(x_found_source_all(NSOURCES_SUBSET_current_size,0:NPROCTOT-1))
  allocate(y_found_source_all(NSOURCES_SUBSET_current_size,0:NPROCTOT-1))
  allocate(z_found_source_all(NSOURCES_SUBSET_current_size,0:NPROCTOT-1))

  allocate(xi_source_subset(NSOURCES_SUBSET_current_size))
  allocate(eta_source_subset(NSOURCES_SUBSET_current_size))
  allocate(gamma_source_subset(NSOURCES_SUBSET_current_size))

  allocate(x_found_source(NSOURCES_SUBSET_current_size))
  allocate(y_found_source(NSOURCES_SUBSET_current_size))
  allocate(z_found_source(NSOURCES_SUBSET_current_size))

! make sure we clean the subset array before the gather
  ispec_selected_source_subset(:) = 0

! loop over sources within this subset
  do isource_in_this_subset = 1,NSOURCES_SUBSET_current_size

! mapping from source number in current subset to real source number in all the subsets
  isource = isource_in_this_subset + isources_already_done

! convert geographic latitude lat (degrees) to geocentric colatitude theta (radians)
  if(ASSUME_PERFECT_SPHERE) then
    theta = PI/2.0d0 - lat(isource)*PI/180.0d0
  else
    theta = PI/2.0d0 - atan(0.99329534d0*dtan(lat(isource)*PI/180.0d0))
  endif

  phi = long(isource)*PI/180.0d0
  call reduce(theta,phi)

! get the moment tensor
  Mrr = moment_tensor(1,isource)
  Mtt = moment_tensor(2,isource)
  Mpp = moment_tensor(3,isource)
  Mrt = moment_tensor(4,isource)
  Mrp = moment_tensor(5,isource)
  Mtp = moment_tensor(6,isource)

! convert from a spherical to a Cartesian representation of the moment tensor
  st=dsin(theta)
  ct=dcos(theta)
  sp=dsin(phi)
  cp=dcos(phi)

  Mxx(isource)=st*st*cp*cp*Mrr+ct*ct*cp*cp*Mtt+sp*sp*Mpp &
      +2.0d0*st*ct*cp*cp*Mrt-2.0d0*st*sp*cp*Mrp-2.0d0*ct*sp*cp*Mtp
  Myy(isource)=st*st*sp*sp*Mrr+ct*ct*sp*sp*Mtt+cp*cp*Mpp &
      +2.0d0*st*ct*sp*sp*Mrt+2.0d0*st*sp*cp*Mrp+2.0d0*ct*sp*cp*Mtp
  Mzz(isource)=ct*ct*Mrr+st*st*Mtt-2.0d0*st*ct*Mrt
  Mxy(isource)=st*st*sp*cp*Mrr+ct*ct*sp*cp*Mtt-sp*cp*Mpp &
      +2.0d0*st*ct*sp*cp*Mrt+st*(cp*cp-sp*sp)*Mrp+ct*(cp*cp-sp*sp)*Mtp
  Mxz(isource)=st*ct*cp*Mrr-st*ct*cp*Mtt &
      +(ct*ct-st*st)*cp*Mrt-ct*sp*Mrp+st*sp*Mtp
  Myz(isource)=st*ct*sp*Mrr-st*ct*sp*Mtt &
      +(ct*ct-st*st)*sp*Mrt+ct*cp*Mrp-st*cp*Mtp

! record three components for each station
  do iorientation = 1,3

!   North
    if(iorientation == 1) then
      stazi = 0.d0
      stdip = 0.d0
!   East
    else if(iorientation == 2) then
      stazi = 90.d0
      stdip = 0.d0
!   Vertical
    else if(iorientation == 3) then
      stazi = 0.d0
      stdip = - 90.d0
    else
      call exit_MPI(myrank,'incorrect orientation')
    endif

!   get the orientation of the seismometer
    thetan=(90.0d0+stdip)*PI/180.0d0
    phin=stazi*PI/180.0d0

! we use the same convention as in Harvard normal modes for the orientation

!   vertical component
    n(1) = dcos(thetan)
!   N-S component
    n(2) = - dsin(thetan)*dcos(phin)
!   E-W component
    n(3) = dsin(thetan)*dsin(phin)

!   get the Cartesian components of n in the model: nu
    nu_source(iorientation,1,isource) = n(1)*st*cp+n(2)*ct*cp-n(3)*sp
    nu_source(iorientation,2,isource) = n(1)*st*sp+n(2)*ct*sp+n(3)*cp
    nu_source(iorientation,3,isource) = n(1)*ct-n(2)*st

  enddo

! normalized source radius
  r0 = R_UNIT_SPHERE

  if(ELLIPTICITY) then
    if(TOPOGRAPHY) then
      call get_topo_bathy(lat(isource),long(isource),elevation,ibathy_topo)
      r0 = r0 + elevation/R_EARTH
    endif
    dcost = dcos(theta)
    p20 = 0.5d0*(3.0d0*dcost*dcost-1.0d0)
    radius = r0 - depth(isource)*1000.0d0/R_EARTH
    call spline_evaluation(rspl,espl,espl2,nspl,radius,ell)
    r0 = r0*(1.0d0-(2.0d0/3.0d0)*ell*p20)
  endif

! compute the Cartesian position of the source
  r_target_source = r0 - depth(isource)*1000.0d0/R_EARTH
  x_target_source = r_target_source*dsin(theta)*dcos(phi)
  y_target_source = r_target_source*dsin(theta)*dsin(phi)
  z_target_source = r_target_source*dcos(theta)

  if(myrank == 0) write(IOVTK,*) sngl(x_target_source),sngl(y_target_source),sngl(z_target_source)

! set distance to huge initial value
  distmin = HUGEVAL

! compute typical size of elements at the surface
  typical_size = TWO_PI * R_UNIT_SPHERE / (4.*NEX_XI)

! use 10 times the distance as a criterion for source detection
  typical_size = 10. * typical_size

! flag to check that we located at least one target element
  located_target = .false.

  do ispec = 1,nspec

! exclude elements that are too far from target
  iglob = ibool(1,1,1,ispec)
  dist = dsqrt((x_target_source - dble(xstore(iglob)))**2 &
             + (y_target_source - dble(ystore(iglob)))**2 &
             + (z_target_source - dble(zstore(iglob)))**2)
  if(USE_DISTANCE_CRITERION .and. dist > typical_size) cycle

  located_target = .true.

! loop only on points inside the element
! exclude edges to ensure this point is not shared with other elements
  do k = 2,NGLLZ-1
    do j = 2,NGLLY-1
      do i = 2,NGLLX-1

!       keep this point if it is closer to the receiver
        iglob = ibool(i,j,k,ispec)
        dist = dsqrt((x_target_source - dble(xstore(iglob)))**2 &
                    +(y_target_source - dble(ystore(iglob)))**2 &
                    +(z_target_source - dble(zstore(iglob)))**2)
        if(dist < distmin) then
          distmin = dist
          ispec_selected_source_subset(isource_in_this_subset) = ispec
          ix_initial_guess_source = i
          iy_initial_guess_source = j
          iz_initial_guess_source = k
        endif

      enddo
    enddo
  enddo

! end of loop on all the elements in current slice
  enddo

! *******************************************
! find the best (xi,eta,gamma) for the source
! *******************************************

! if we have not located a target element, the source is not in this slice
! therefore use first element only for fictitious iterative search
  if(.not. located_target) then
    ispec_selected_source_subset(isource_in_this_subset)=1
    ix_initial_guess_source = 2
    iy_initial_guess_source = 2
    iz_initial_guess_source = 2
  endif

! use initial guess in xi, eta and gamma
  xi = xigll(ix_initial_guess_source)
  eta = yigll(iy_initial_guess_source)
  gamma = zigll(iz_initial_guess_source)

! define coordinates of the control points of the element

  do ia=1,NGNOD

    if(iaddx(ia) == 0) then
      iax = 1
    else if(iaddx(ia) == 1) then
      iax = (NGLLX+1)/2
    else if(iaddx(ia) == 2) then
      iax = NGLLX
    else
      call exit_MPI(myrank,'incorrect value of iaddx')
    endif

    if(iaddy(ia) == 0) then
      iay = 1
    else if(iaddy(ia) == 1) then
      iay = (NGLLY+1)/2
    else if(iaddy(ia) == 2) then
      iay = NGLLY
    else
      call exit_MPI(myrank,'incorrect value of iaddy')
    endif

    if(iaddr(ia) == 0) then
      iaz = 1
    else if(iaddr(ia) == 1) then
      iaz = (NGLLZ+1)/2
    else if(iaddr(ia) == 2) then
      iaz = NGLLZ
    else
      call exit_MPI(myrank,'incorrect value of iaddr')
    endif

    iglob = ibool(iax,iay,iaz,ispec_selected_source_subset(isource_in_this_subset))
    xelm(ia) = dble(xstore(iglob))
    yelm(ia) = dble(ystore(iglob))
    zelm(ia) = dble(zstore(iglob))

  enddo

! iterate to solve the non linear system
  do iter_loop = 1,NUM_ITER

! recompute jacobian for the new point
    call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

! compute distance to target location
  dx = - (x - x_target_source)
  dy = - (y - y_target_source)
  dz = - (z - z_target_source)

! compute increments
  dxi  = xix*dx + xiy*dy + xiz*dz
  deta = etax*dx + etay*dy + etaz*dz
  dgamma = gammax*dx + gammay*dy + gammaz*dz

! update values
  xi = xi + dxi
  eta = eta + deta
  gamma = gamma + dgamma

! impose that we stay in that element
! (useful if user gives a source outside the mesh for instance)
  if (xi > 1.d0) xi = 1.d0
  if (xi < -1.d0) xi = -1.d0
  if (eta > 1.d0) eta = 1.d0
  if (eta < -1.d0) eta = -1.d0
  if (gamma > 1.d0) gamma = 1.d0
  if (gamma < -1.d0) gamma = -1.d0

  enddo

! compute final coordinates of point found
  call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z,xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

! store xi,eta,gamma and x,y,z of point found
  xi_source_subset(isource_in_this_subset) = xi
  eta_source_subset(isource_in_this_subset) = eta
  gamma_source_subset(isource_in_this_subset) = gamma
  x_found_source(isource_in_this_subset) = x
  y_found_source(isource_in_this_subset) = y
  z_found_source(isource_in_this_subset) = z

! compute final distance between asked and found (converted to km)
  final_distance_source_subset(isource_in_this_subset) = dsqrt((x_target_source-x_found_source(isource_in_this_subset))**2 + &
    (y_target_source-y_found_source(isource_in_this_subset))**2 + &
    (z_target_source-z_found_source(isource_in_this_subset))**2)*R_EARTH/1000.d0

! end of loop on all the sources
  enddo

! now gather information from all the nodes
! use -1 as a flag to detect if gather fails for some reason
  ispec_selected_source_all(:,:) = -1
  call MPI_GATHER(ispec_selected_source_subset,NSOURCES_SUBSET_current_size,MPI_INTEGER, &
                  ispec_selected_source_all,NSOURCES_SUBSET_current_size,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(xi_source_subset,NSOURCES_SUBSET_current_size,MPI_DOUBLE_PRECISION, &
                  xi_source_all,NSOURCES_SUBSET_current_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(eta_source_subset,NSOURCES_SUBSET_current_size,MPI_DOUBLE_PRECISION, &
                  eta_source_all,NSOURCES_SUBSET_current_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(gamma_source_subset,NSOURCES_SUBSET_current_size,MPI_DOUBLE_PRECISION, &
                  gamma_source_all,NSOURCES_SUBSET_current_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(final_distance_source_subset,NSOURCES_SUBSET_current_size,MPI_DOUBLE_PRECISION, &
    final_distance_source_all,NSOURCES_SUBSET_current_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(x_found_source,NSOURCES_SUBSET_current_size,MPI_DOUBLE_PRECISION, &
    x_found_source_all,NSOURCES_SUBSET_current_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(y_found_source,NSOURCES_SUBSET_current_size,MPI_DOUBLE_PRECISION, &
    y_found_source_all,NSOURCES_SUBSET_current_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(z_found_source,NSOURCES_SUBSET_current_size,MPI_DOUBLE_PRECISION, &
    z_found_source_all,NSOURCES_SUBSET_current_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

! this is executed by main process only
  if(myrank == 0) then

! check that the gather operation went well
  if(minval(ispec_selected_source_all) <= 0) call exit_MPI(myrank,'gather operation failed for source')

! loop on all the sources within subsets
  do isource_in_this_subset = 1,NSOURCES_SUBSET_current_size

! mapping from source number in current subset to real source number in all the subsets
     isource = isources_already_done + isource_in_this_subset

! loop on all the results to determine the best slice
  distmin = HUGEVAL
  do iprocloop = 0,NPROCTOT-1
    if(final_distance_source_all(isource_in_this_subset,iprocloop) < distmin) then
      distmin = final_distance_source_all(isource_in_this_subset,iprocloop)
      islice_selected_source(isource) = iprocloop
      ispec_selected_source(isource) = ispec_selected_source_all(isource_in_this_subset,iprocloop)
      xi_source(isource) = xi_source_all(isource_in_this_subset,iprocloop)
      eta_source(isource) = eta_source_all(isource_in_this_subset,iprocloop)
      gamma_source(isource) = gamma_source_all(isource_in_this_subset,iprocloop)
      x_found_source(isource_in_this_subset) = x_found_source_all(isource_in_this_subset,iprocloop)
      y_found_source(isource_in_this_subset) = y_found_source_all(isource_in_this_subset,iprocloop)
      z_found_source(isource_in_this_subset) = z_found_source_all(isource_in_this_subset,iprocloop)
    endif
  enddo
  final_distance_source(isource) = distmin

    write(IMAIN,*)
    write(IMAIN,*) '*************************************'
    write(IMAIN,*) ' locating source ',isource
    write(IMAIN,*) '*************************************'
    write(IMAIN,*)
    write(IMAIN,*) 'source located in slice ',islice_selected_source(isource_in_this_subset)
    write(IMAIN,*) '               in element ',ispec_selected_source(isource_in_this_subset)
    write(IMAIN,*)
    write(IMAIN,*) '   xi coordinate of source in that element: ',xi_source(isource)
    write(IMAIN,*) '  eta coordinate of source in that element: ',eta_source(isource)
    write(IMAIN,*) 'gamma coordinate of source in that element: ',gamma_source(isource)

! add message if source is a Heaviside
    if(hdur(isource) < 5.*DT) then
      write(IMAIN,*)
      write(IMAIN,*) 'Source time function is a Heaviside, convolve later'
      write(IMAIN,*)
    endif

    write(IMAIN,*)
    write(IMAIN,*) ' half duration: ',hdur(isource),' seconds'
    write(IMAIN,*) '    time shift: ',t_cmt(isource),' seconds'

! get latitude, longitude and depth of the source that will be used
    call xyz_2_rthetaphi_dble(x_found_source(isource_in_this_subset),y_found_source(isource_in_this_subset), &
           z_found_source(isource_in_this_subset),r_found_source,theta_source(isource),phi_source(isource))
    call reduce(theta_source(isource),phi_source(isource))

! convert geocentric to geographic colatitude
    colat_source = PI/2.0d0-datan(1.006760466d0*dcos(theta_source(isource))/dmax1(TINYVAL,dsin(theta_source(isource))))
    if(phi_source(isource)>PI) phi_source(isource)=phi_source(isource)-TWO_PI

    write(IMAIN,*)
    write(IMAIN,*) 'original (requested) position of the source:'
    write(IMAIN,*)
    write(IMAIN,*) '      latitude: ',lat(isource)
    write(IMAIN,*) '     longitude: ',long(isource)
    write(IMAIN,*) '         depth: ',depth(isource),' km'
    write(IMAIN,*)

! compute real position of the source
    write(IMAIN,*) 'position of the source that will be used:'
    write(IMAIN,*)
    write(IMAIN,*) '      latitude: ',(PI/2.0d0-colat_source)*180.0d0/PI
    write(IMAIN,*) '     longitude: ',phi_source(isource)*180.0d0/PI
    write(IMAIN,*) '         depth: ',(r0-r_found_source)*R_EARTH/1000.0d0,' km'
    write(IMAIN,*)

! display error in location estimate
    write(IMAIN,*) 'error in location of the source: ',sngl(final_distance_source(isource)),' km'

! add warning if estimate is poor
! (usually means source outside the mesh given by the user)
    if(final_distance_source(isource) > 50.d0) then
      write(IMAIN,*)
      write(IMAIN,*) '*****************************************************'
      write(IMAIN,*) '*****************************************************'
      write(IMAIN,*) '***** WARNING: source location estimate is poor *****'
      write(IMAIN,*) '*****************************************************'
      write(IMAIN,*) '*****************************************************'
    endif

! print source time function and spectrum
  if(PRINT_SOURCE_TIME_FUNCTION) then

  write(IMAIN,*)
  write(IMAIN,*) 'printing the source-time function'

! print the source-time function
  if(NSOURCES == 1) then
    plot_file = '/plot_source_time_function.txt'
  else
   if(isource < 10) then
      write(plot_file,"('/plot_source_time_function',i1,'.txt')") isource
    elseif(isource < 100) then
      write(plot_file,"('/plot_source_time_function',i2,'.txt')") isource
    else
      write(plot_file,"('/plot_source_time_function',i3,'.txt')") isource
    endif
  endif
  open(unit=27,file=trim(OUTPUT_FILES)//plot_file,status='unknown')

  scalar_moment = 0.
  do i = 1,6
    scalar_moment = scalar_moment + moment_tensor(i,isource)**2
  enddo
  scalar_moment = dsqrt(scalar_moment/2.)

  do it=1,NSTEP
    time_source = dble(it-1)*DT-t0-t_cmt(isource)
    write(27,*) sngl(dble(it-1)*DT-t0),sngl(scalar_moment*comp_source_time_function(time_source,hdur_gaussian(isource)))
  enddo
  close(27)

  write(IMAIN,*)
  write(IMAIN,*) 'printing the source spectrum'

! print the spectrum of the derivative of the source from 0 to 1/8 Hz
  if(NSOURCES == 1) then
   plot_file = '/plot_source_spectrum.txt'
  else
   if(isource < 10) then
      write(plot_file,"('/plot_source_spectrum',i1,'.txt')") isource
    elseif(isource < 100) then
      write(plot_file,"('/plot_source_spectrum',i2,'.txt')") isource
    else
      write(plot_file,"('/plot_source_spectrum',i3,'.txt')") isource
    endif
  endif
  open(unit=27,file=trim(OUTPUT_FILES)//plot_file,status='unknown')

  do iom=1,NSAMP_PLOT_SOURCE
    om=TWO_PI*(1.0d0/8.0d0)*(iom-1)/dble(NSAMP_PLOT_SOURCE-1)
    write(27,*) sngl(om/TWO_PI),sngl(scalar_moment*om*comp_source_spectrum(om,hdur(isource)))
  enddo
  close(27)

  endif

  enddo ! end of loop on all the sources within current source subset

  endif ! end of section executed by main process only

! deallocate arrays specific to each subset
  deallocate(final_distance_source_subset)
  deallocate(ispec_selected_source_subset)
  deallocate(ispec_selected_source_all)
  deallocate(xi_source_all,eta_source_all,gamma_source_all,final_distance_source_all)
  deallocate(x_found_source_all,y_found_source_all,z_found_source_all)
  deallocate(xi_source_subset,eta_source_subset,gamma_source_subset)
  deallocate(x_found_source,y_found_source,z_found_source)

  enddo ! end of loop over all source subsets

! display maximum error in location estimate
  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'maximum error in location of the sources: ',sngl(maxval(final_distance_source)),' km'
    write(IMAIN,*)
  endif


! main process broadcasts the results to all the slices
  call MPI_BCAST(islice_selected_source,NSOURCES,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(ispec_selected_source,NSOURCES,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(xi_source,NSOURCES,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(eta_source,NSOURCES,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(gamma_source,NSOURCES,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

! elapsed time since beginning of source detection
  if(myrank == 0) then
    tCPU = MPI_WTIME() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for detection of sources in seconds = ',tCPU
    write(IMAIN,*)
    write(IMAIN,*) 'End of source detection - done'
    write(IMAIN,*)
  endif

  end subroutine locate_sources

