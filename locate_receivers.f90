!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology September 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

!----
!---- locate_receivers finds the correct position of the receivers
!----

! to locate the receivers we loop on elements located at the surface only

  subroutine locate_receivers(myrank,DT,NSTEP,nspec,nglob,idoubling,ibool, &
                 xstore,ystore,zstore,xigll,yigll, &
                 nrec,islice_selected_rec,ispec_selected_rec, &
                 xi_receiver,eta_receiver,gamma_receiver,station_name,network_name,nu, &
                 yr,jda,ho,mi,sec,NPROCTOT,ELLIPTICITY,TOPOGRAPHY, &
                 theta_source,phi_source, &
                 rspl,espl,espl2,nspl,ibathy_topo)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer NPROCTOT

  logical ELLIPTICITY,TOPOGRAPHY

  integer nspl
  double precision rspl(NR),espl(NR),espl2(NR)

  integer nspec,nglob,nrec,myrank

  integer yr,jda,ho,mi
  double precision sec

  integer idoubling(nspec)
  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)
  integer NSTEP
  double precision DT

! arrays containing coordinates of the points
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore,ystore,zstore

! Gauss-Lobatto-Legendre points of integration
  double precision xigll(NGLLX),yigll(NGLLY)

! use integer array to store values
  integer ibathy_topo(NX_BATHY,NY_BATHY)

  integer, allocatable, dimension(:) :: ix_initial_guess,iy_initial_guess

  integer iorientation
  integer iprocloop
  integer nrec_dummy
  double precision stazi,stdip

  double precision, allocatable, dimension(:) :: x_target,y_target,z_target
  double precision, allocatable, dimension(:) :: epidist
  double precision, allocatable, dimension(:) :: x_found,y_found,z_found
  double precision, allocatable, dimension(:,:) :: x_found_all,y_found_all,z_found_all

  integer irec
  integer i,j,k,ispec,iglob
  integer ier

  double precision ell
  double precision elevation
  double precision n(3)
  double precision thetan,phin
  double precision sint,cost,sinp,cosp
  double precision r0,dcost,p20
  double precision theta,phi
  double precision theta_source,phi_source
  double precision dist
  double precision xi,eta,gamma,dx,dy,dz,dxi,deta,dgamma

! topology of the control points of the surface element
  integer iax,iay,iaz
  integer iaddx(NGNOD),iaddy(NGNOD),iaddz(NGNOD)

! coordinates of the control points of the surface element
  double precision xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  integer iter_loop,ispec_iterate

  integer ia
  double precision x,y,z
  double precision xix,xiy,xiz
  double precision etax,etay,etaz
  double precision gammax,gammay,gammaz

! timer MPI
  double precision time_start,tCPU

! use dynamic allocation
  double precision, dimension(:), allocatable :: final_distance
  double precision, dimension(:,:), allocatable :: final_distance_all
  double precision distmin,final_distance_max

! receiver information
! timing information for the stations
! station information for writing the seismograms
  integer nsamp
  integer, dimension(nrec) :: islice_selected_rec,ispec_selected_rec
  double precision, dimension(nrec) :: xi_receiver,eta_receiver,gamma_receiver
  double precision, dimension(3,3,nrec) :: nu
  character(len=8), dimension(nrec) :: station_name,network_name

  integer, allocatable, dimension(:,:) :: ispec_selected_rec_all
  double precision, allocatable, dimension(:) :: stlat,stlon,stele,stbur
  double precision, allocatable, dimension(:,:) :: xi_receiver_all,eta_receiver_all,gamma_receiver_all

! **************

! get MPI starting time
  time_start = MPI_WTIME()

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' locating receivers'
    write(IMAIN,*) '********************'
    write(IMAIN,*)
  endif

! define topology of the control element
  call usual_hex_nodes(iaddx,iaddy,iaddz)

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '****************************'
    write(IMAIN,*) 'reading receiver information'
    write(IMAIN,*) '****************************'
    write(IMAIN,*)
  endif

! get number of stations from receiver file
  open(unit=1,file='DATA/STATIONS',status='old')
  read(1,*) nrec_dummy

  if(nrec_dummy /= nrec) call exit_MPI(myrank,'problem with number of receivers')

! allocate memory for arrays using number of stations
  allocate(stlat(nrec))
  allocate(stlon(nrec))
  allocate(stele(nrec))
  allocate(stbur(nrec))
  allocate(epidist(nrec))

  allocate(ix_initial_guess(nrec))
  allocate(iy_initial_guess(nrec))
  allocate(x_target(nrec))
  allocate(y_target(nrec))
  allocate(z_target(nrec))
  allocate(x_found(nrec))
  allocate(y_found(nrec))
  allocate(z_found(nrec))
  allocate(final_distance(nrec))

  allocate(ispec_selected_rec_all(nrec,0:NPROCTOT-1))
  allocate(xi_receiver_all(nrec,0:NPROCTOT-1))
  allocate(eta_receiver_all(nrec,0:NPROCTOT-1))
  allocate(gamma_receiver_all(nrec,0:NPROCTOT-1))
  allocate(x_found_all(nrec,0:NPROCTOT-1))
  allocate(y_found_all(nrec,0:NPROCTOT-1))
  allocate(z_found_all(nrec,0:NPROCTOT-1))
  allocate(final_distance_all(nrec,0:NPROCTOT-1))

! loop on all the stations
  do irec=1,nrec

! set distance to huge initial value
  distmin=HUGEVAL

    read(1,*) station_name(irec),network_name(irec),stlat(irec),stlon(irec),stele(irec),stbur(irec)

!
! full Earth, use latitude and longitude
!

!     convert geographic latitude stlat (degrees)
!     to geocentric colatitude theta (radians)
      theta=PI/2.0d0-atan(0.99329534d0*dtan(stlat(irec)*PI/180.0d0))
      phi=stlon(irec)*PI/180.0d0
      call reduce(theta,phi)

! compute epicentral distance
      epidist(irec) = dacos(dcos(theta)*dcos(theta_source) + &
              dsin(theta)*dsin(theta_source)*dcos(phi-phi_source))*180.0d0/PI

! print some information about stations
      if(myrank == 0) &
        write(IMAIN,"('Station #',i5,':  ',a5,'.',a2,'      epicentral distance:  ',f10.3,' degrees')") &
          irec,station_name(irec),network_name(irec),epidist(irec)

! record three components for each station
    do iorientation = 1,3

!     North
      if(iorientation == 1) then
        stazi = 0.d0
        stdip = 0.d0
!     East
      else if(iorientation == 2) then
        stazi = 90.d0
        stdip = 0.d0
!     Vertical
      else if(iorientation == 3) then
        stazi = 0.d0
        stdip = - 90.d0
      else
        call exit_MPI(myrank,'incorrect orientation')
      endif

!     get the orientation of the seismometer
      thetan=(90.0d0+stdip)*PI/180.0d0
      phin=stazi*PI/180.0d0

! we use the same convention as in Harvard normal modes for the orientation

!     vertical component
      n(1) = dcos(thetan)
!     N-S component
      n(2) = - dsin(thetan)*dcos(phin)
!     E-W component
      n(3) = dsin(thetan)*dsin(phin)

!     get the Cartesian components of n in the model: nu
      sint = dsin(theta)
      cost = dcos(theta)
      sinp = dsin(phi)
      cosp = dcos(phi)
      nu(iorientation,1,irec) = n(1)*sint*cosp+n(2)*cost*cosp-n(3)*sinp
      nu(iorientation,2,irec) = n(1)*sint*sinp+n(2)*cost*sinp+n(3)*cosp
      nu(iorientation,3,irec) = n(1)*cost-n(2)*sint

    enddo

!     ellipticity
      r0=1.0d0
      if(ELLIPTICITY) then
        if(TOPOGRAPHY) then
           call get_topo_bathy(stlat(irec),stlon(irec),elevation,ibathy_topo)
           r0 = r0 + elevation/R_EARTH
        endif
        dcost=dcos(theta)
        p20=0.5d0*(3.0d0*dcost*dcost-1.0d0)
        call splint(rspl,espl,espl2,nspl,r0,ell)
        r0=r0*(1.0d0-(2.0d0/3.0d0)*ell*p20)
      endif

! subtract station burial depth (in meters)
      r0 = r0 - stbur(irec)/R_EARTH

! compute the Cartesian position of the receiver
      x_target(irec) = r0*dsin(theta)*dcos(phi)
      y_target(irec) = r0*dsin(theta)*dsin(phi)
      z_target(irec) = r0*dcos(theta)

! examine top of the elements only (receivers always at the surface)
      k = NGLLZ

      do ispec=1,nspec

! examine elements located in the crust only (receivers always at the surface)
      if(idoubling(ispec) == IFLAG_CRUST) then

! loop only on points inside the element
! exclude edges to ensure this point is not shared with other elements
        do j=2,NGLLY-1
          do i=2,NGLLX-1

            iglob = ibool(i,j,k,ispec)
            dist = dsqrt((x_target(irec)-dble(xstore(iglob)))**2 &
                        +(y_target(irec)-dble(ystore(iglob)))**2 &
                        +(z_target(irec)-dble(zstore(iglob)))**2)

!           keep this point if it is closer to the receiver
            if(dist < distmin) then
              distmin = dist
              ispec_selected_rec(irec) = ispec
              ix_initial_guess(irec) = i
              iy_initial_guess(irec) = j
            endif

          enddo
        enddo

      endif

! end of loop on all the spectral elements in current slice
      enddo

! end of loop on all the stations
  enddo

! close receiver file
  close(1)

! create RECORDHEADER file with usual format for normal-mode codes

  if(myrank == 0) then

! create file for QmX Harvard
! Harvard format does not support the network name
! therefore only the station name is included below
! compute total number of samples for normal modes with 1 sample per second
    open(unit=1,file='OUTPUT_FILES/RECORDHEADERS',status='unknown')
    nsamp = nint(dble(NSTEP-1)*DT)
    do irec = 1,nrec
      write(1,"(a8,1x,a3,6x,f8.4,1x,f9.4,1x,f6.1,1x,f6.1,f6.1,1x,f6.1,1x,f12.4,1x,i7,1x,i4,1x,i3,1x,i2,1x,i2,1x,f6.3)") &
        station_name(irec),'LHN',stlat(irec),stlon(irec),stele(irec),stbur(irec), &
        0.,0.,1.,nsamp,yr,jda,ho,mi,sec
      write(1,"(a8,1x,a3,6x,f8.4,1x,f9.4,1x,f6.1,1x,f6.1,f6.1,1x,f6.1,1x,f12.4,1x,i7,1x,i4,1x,i3,1x,i2,1x,i2,1x,f6.3)") &
        station_name(irec),'LHE',stlat(irec),stlon(irec),stele(irec),stbur(irec), &
        90.,0.,1.,nsamp,yr,jda,ho,mi,sec
      write(1,"(a8,1x,a3,6x,f8.4,1x,f9.4,1x,f6.1,1x,f6.1,f6.1,1x,f6.1,1x,f12.4,1x,i7,1x,i4,1x,i3,1x,i2,1x,i2,1x,f6.3)") &
        station_name(irec),'LHZ',stlat(irec),stlon(irec),stele(irec),stbur(irec), &
        0.,-90.,1.,nsamp,yr,jda,ho,mi,sec
    enddo
    close(1)

  endif

! ****************************************
! find the best (xi,eta) for each receiver
! ****************************************

! loop on all the receivers to iterate in that slice
    do irec = 1,nrec

        ispec_iterate = ispec_selected_rec(irec)

! use initial guess in xi and eta
        xi = xigll(ix_initial_guess(irec))
        eta = yigll(iy_initial_guess(irec))

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

    if(iaddz(ia) == 0) then
      iaz = 1
    else if(iaddz(ia) == 1) then
      iaz = (NGLLZ+1)/2
    else if(iaddz(ia) == 2) then
      iaz = NGLLZ
    else
      call exit_MPI(myrank,'incorrect value of iaddz')
    endif

    iglob = ibool(iax,iay,iaz,ispec_iterate)
    xelm(ia) = dble(xstore(iglob))
    yelm(ia) = dble(ystore(iglob))
    zelm(ia) = dble(zstore(iglob))

  enddo

! iterate to solve the non linear system
  do iter_loop = 1,NUM_ITER

! impose receiver exactly at the surface
    if(.not. RECEIVERS_CAN_BE_BURIED) gamma = 1.d0

! recompute jacobian for the new point
    call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
           xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

! compute distance to target location
  dx = - (x - x_target(irec))
  dy = - (y - y_target(irec))
  dz = - (z - z_target(irec))

! compute increments
! gamma does not change since we know the receiver is exactly on the surface
  dxi  = xix*dx + xiy*dy + xiz*dz
  deta = etax*dx + etay*dy + etaz*dz
  if(RECEIVERS_CAN_BE_BURIED) dgamma = gammax*dx + gammay*dy + gammaz*dz

! update values
  xi = xi + dxi
  eta = eta + deta
  if(RECEIVERS_CAN_BE_BURIED) gamma = gamma + dgamma

! impose that we stay in that element
! (useful if user gives a receiver outside the mesh for instance)
! we can go slightly outside the [1,1] segment since with finite elements
! the polynomial solution is defined everywhere
! can be useful for convergence of iterative scheme with distorted elements
  if (xi > 1.10d0) xi = 1.10d0
  if (xi < -1.10d0) xi = -1.10d0
  if (eta > 1.10d0) eta = 1.10d0
  if (eta < -1.10d0) eta = -1.10d0
  if (gamma > 1.10d0) gamma = 1.10d0
  if (gamma < -1.10d0) gamma = -1.10d0

! end of non linear iterations
  enddo

! impose receiver exactly at the surface after final iteration
  if(.not. RECEIVERS_CAN_BE_BURIED) gamma = 1.d0

! compute final coordinates of point found
  call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
         xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

! store xi,eta and x,y,z of point found
  xi_receiver(irec) = xi
  eta_receiver(irec) = eta
  gamma_receiver(irec) = gamma
  x_found(irec) = x
  y_found(irec) = y
  z_found(irec) = z

! compute final distance between asked and found (converted to km)
  final_distance(irec) = dsqrt((x_target(irec)-x_found(irec))**2 + &
    (y_target(irec)-y_found(irec))**2 + (z_target(irec)-z_found(irec))**2)*R_EARTH/1000.d0

    enddo

! for MPI version, gather information from all the nodes
  ispec_selected_rec_all(:,:) = -1
  call MPI_GATHER(ispec_selected_rec,nrec,MPI_INTEGER,ispec_selected_rec_all,nrec,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  call MPI_GATHER(xi_receiver,nrec,MPI_DOUBLE_PRECISION,xi_receiver_all,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(eta_receiver,nrec,MPI_DOUBLE_PRECISION,eta_receiver_all,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(gamma_receiver,nrec,MPI_DOUBLE_PRECISION,gamma_receiver_all,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(final_distance,nrec,MPI_DOUBLE_PRECISION,final_distance_all,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(x_found,nrec,MPI_DOUBLE_PRECISION,x_found_all,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(y_found,nrec,MPI_DOUBLE_PRECISION,y_found_all,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(z_found,nrec,MPI_DOUBLE_PRECISION,z_found_all,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

! this is executed by main process only
  if(myrank == 0) then

! check that the gather operation went well
  if(any(ispec_selected_rec_all(:,:) == -1)) call exit_MPI(myrank,'gather operation failed for receivers')

! MPI loop on all the results to determine the best slice
  islice_selected_rec(:) = -1
  do irec = 1,nrec
  distmin = HUGEVAL
  do iprocloop = 0,NPROCTOT-1
    if(final_distance_all(irec,iprocloop) < distmin) then
      distmin = final_distance_all(irec,iprocloop)
      islice_selected_rec(irec) = iprocloop
      ispec_selected_rec(irec) = ispec_selected_rec_all(irec,iprocloop)
      xi_receiver(irec) = xi_receiver_all(irec,iprocloop)
      eta_receiver(irec) = eta_receiver_all(irec,iprocloop)
      gamma_receiver(irec) = gamma_receiver_all(irec,iprocloop)
      x_found(irec) = x_found_all(irec,iprocloop)
      y_found(irec) = y_found_all(irec,iprocloop)
      z_found(irec) = z_found_all(irec,iprocloop)
    endif
  enddo
  final_distance(irec) = distmin
  enddo

  do irec=1,nrec

    if(final_distance(irec) == HUGEVAL) call exit_MPI(myrank,'error locating receiver')

    if(DISPLAY_DETAILS_STATIONS) then
      write(IMAIN,*)
      write(IMAIN,*) 'station # ',irec,'    ',station_name(irec),network_name(irec)
      write(IMAIN,*) '     original latitude: ',sngl(stlat(irec))
      write(IMAIN,*) '    original longitude: ',sngl(stlon(irec))
      write(IMAIN,*) '   epicentral distance: ',sngl(epidist(irec))
      write(IMAIN,*) 'closest estimate found: ',sngl(final_distance(irec)),' km away'
      write(IMAIN,*) ' in slice ',islice_selected_rec(irec),' in element ',ispec_selected_rec(irec)
      write(IMAIN,*) ' at xi,eta,gamma coordinates = ',xi_receiver(irec),eta_receiver(irec),gamma_receiver(irec)
    endif

! add warning if estimate is poor
! (usually means receiver outside the mesh given by the user)
    if(final_distance(irec) > 50.d0) then
      write(IMAIN,*) 'station # ',irec,'    ',station_name(irec),network_name(irec)
      write(IMAIN,*) '*******************************************************'
      write(IMAIN,*) '***** WARNING: receiver location estimate is poor *****'
      write(IMAIN,*) '*******************************************************'
    endif

  enddo

! compute maximal distance for all the receivers
    final_distance_max = maxval(final_distance(:))

! display maximum error for all the receivers
    write(IMAIN,*)
    write(IMAIN,*) 'maximum error in location of all the receivers: ',sngl(final_distance_max),' km'

! add warning if estimate is poor
! (usually means receiver outside the mesh given by the user)
    if(final_distance_max > 50.d0) then
      write(IMAIN,*)
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '***** WARNING: at least one receiver is poorly located *****'
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '************************************************************'
    endif

! write the list of stations and associated epicentral distance
  open(unit=27,file='OUTPUT_FILES/output_list_stations.txt',status='unknown')
  write(27,*)
  write(27,*) 'total number of stations: ',nrec
  write(27,*)
  do irec=1,nrec
    write(27,"(a5,'.',a2,'   epicentral distance ',f8.3,' deg')") station_name(irec),network_name(irec),epidist(irec)
  enddo
  close(27)

! elapsed time since beginning of mesh generation
  tCPU = MPI_WTIME() - time_start
  write(IMAIN,*)
  write(IMAIN,*) 'Elapsed time for receiver detection in seconds = ',tCPU
  write(IMAIN,*)
  write(IMAIN,*) 'End of receiver detection - done'
  write(IMAIN,*)

  endif    ! end of section executed by main process only

! main process broadcasts the results to all the slices
  call MPI_BCAST(islice_selected_rec,nrec,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(ispec_selected_rec,nrec,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(xi_receiver,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(eta_receiver,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(gamma_receiver,nrec,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

! deallocate arrays
  deallocate(stlat)
  deallocate(stlon)
  deallocate(stele)
  deallocate(stbur)
  deallocate(epidist)
  deallocate(ix_initial_guess)
  deallocate(iy_initial_guess)
  deallocate(x_target)
  deallocate(y_target)
  deallocate(z_target)
  deallocate(x_found)
  deallocate(y_found)
  deallocate(z_found)
  deallocate(final_distance)
  deallocate(ispec_selected_rec_all)
  deallocate(xi_receiver_all)
  deallocate(eta_receiver_all)
  deallocate(gamma_receiver_all)
  deallocate(x_found_all)
  deallocate(y_found_all)
  deallocate(z_found_all)
  deallocate(final_distance_all)

  end subroutine locate_receivers

