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

!----
!---- locate_receivers finds the correct position of the receivers
!----

  subroutine locate_receivers(yr,jda,ho,mi,sec)

  use constants_solver, only: &
    ELLIPTICITY_VAL,NCHUNKS_VAL,NPROCTOT_VAL,NDIM, &
    MAX_LENGTH_STATION_NAME,MAX_LENGTH_NETWORK_NAME, &
    DISPLAY_DETAILS_STATIONS,nrec_SUBSET_MAX, &
    THRESHOLD_EXCLUDE_STATION, &
    HUGEVAL,IMAIN,IOUT,IOUT_VTK, &
    DEGREES_TO_RADIANS,RADIANS_TO_DEGREES,R_UNIT_SPHERE

  use shared_parameters, only: OUTPUT_FILES,R_PLANET

  use specfem_par, only: &
    myrank,DT,NSTEP, &
    nrec,islice_selected_rec,ispec_selected_rec, &
    xi_receiver,eta_receiver,gamma_receiver,station_name,network_name, &
    stlat,stlon,stele,stbur,nu_rec,receiver_final_distance_max, &
    RECEIVERS_CAN_BE_BURIED, &
    ibathy_topo,TOPOGRAPHY

  use specfem_par, only: rspl_ellip,ellipicity_spline,ellipicity_spline2,nspl_ellip

  use specfem_par, only: source_theta_ref,source_phi_ref

  implicit none

  integer,intent(in) :: yr,jda,ho,mi
  double precision,intent(in) :: sec

  ! local parameters
  integer :: iprocloop
  integer :: irec,i
  integer :: ier

  integer, dimension(nrec) :: islice_selected_found,ispec_selected_found
  double precision, dimension(nrec) :: xi_receiver_found,eta_receiver_found,gamma_receiver_found
  double precision, dimension(3,3,nrec) :: nu_found
  integer :: nrec_found

  ! point locations
  double precision, allocatable, dimension(:,:) :: xyz_target
  double precision, allocatable, dimension(:,:) :: xyz_found

  double precision, allocatable, dimension(:,:) :: xyz_found_subset
  double precision, allocatable, dimension(:,:,:) :: xyz_found_all

  double precision, dimension(nrec) :: stlat_found,stlon_found,stele_found,stbur_found,epidist_found

  double precision, allocatable, dimension(:) :: epidist

  integer :: nrec_SUBSET_current_size
  integer :: irec_in_this_subset,irec_already_done

  integer, allocatable, dimension(:) :: ispec_selected_subset
  integer, allocatable, dimension(:,:) :: ispec_selected_all

  double precision, dimension(:), allocatable :: final_distance
  double precision, allocatable, dimension(:) :: final_distance_subset
  double precision, dimension(:,:), allocatable :: final_distance_all

  double precision, allocatable, dimension(:) :: xi_subset,eta_subset,gamma_subset
  double precision, allocatable, dimension(:,:) :: xi_all,eta_all,gamma_all

  double precision :: lat,lon,radius,depth,r_target

  double precision :: theta,phi
  double precision :: sint,cost,sinp,cosp,dist

  double precision :: elevation
  double precision :: r0

  double precision :: distmin_not_squared
  double precision :: x_target,y_target,z_target
  double precision :: x,y,z
  double precision :: xi,eta,gamma

  integer :: ispec_selected

  integer :: iorientation
  double precision :: stazi,stdip
  double precision :: n(3),thetan,phin

! receiver information
! timing information for the stations
! station information for writing the seismograms
  integer :: nsamp

  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name_found
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name_found

  character(len=2) :: bic

  ! sorting order
  integer, allocatable, dimension(:) :: irec_dist_ordered

  ! timer MPI
  double precision :: time_start,tCPU
  double precision, external :: wtime

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' locating receivers'
    write(IMAIN,*) '********************'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! get MPI starting time
  time_start = wtime()

  ! allocate memory for arrays using number of stations
  allocate(epidist(nrec), &
           xyz_target(NDIM,nrec), &
           xyz_found(NDIM,nrec), &
           final_distance(nrec),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating temporary receiver arrays')

  ! read that STATIONS file on the main
  call read_receiver_locations()

  ! loop on all the stations to locate them in the mesh
  do irec = 1,nrec

    ! station lat/lon in degrees
    lat = stlat(irec)
    lon = stlon(irec)

    ! limits longitude to [0.0,360.0]
    if (lon < 0.d0 ) lon = lon + 360.d0
    if (lon > 360.d0 ) lon = lon - 360.d0

    ! converts geographic latitude stlat (degrees) to geocentric colatitude theta (radians)
    call lat_2_geocentric_colat_dble(lat,theta,ELLIPTICITY_VAL)

    ! longitude
    phi = lon*DEGREES_TO_RADIANS

    ! theta to [0,PI] and phi to [0,2PI]
    call reduce(theta,phi)

    sint = sin(theta)
    cost = cos(theta)
    sinp = sin(phi)
    cosp = cos(phi)

    ! compute epicentral distance to reference source position (in radians)
    call get_greatcircle_distance(theta,phi,source_theta_ref,source_phi_ref,dist)

    epidist(irec) = dist * RADIANS_TO_DEGREES

    ! record three components for each station
    do iorientation = 1,3
      !     North
      if (iorientation == 1) then
        stazi = 0.d0
        stdip = 0.d0
      !     East
      else if (iorientation == 2) then
        stazi = 90.d0
        stdip = 0.d0
      !     Vertical
      else if (iorientation == 3) then
        stazi = 0.d0
        stdip = - 90.d0
      else
        call exit_MPI(myrank,'incorrect orientation')
      endif

      ! get the orientation of the seismometer
      thetan = (90.0d0+stdip)*DEGREES_TO_RADIANS
      phin = stazi*DEGREES_TO_RADIANS

      ! we use the same convention as in Harvard normal modes for the orientation

      !     vertical component
      n(1) = cos(thetan)
      !     N-S component
      n(2) = - sin(thetan)*cos(phin)
      !     E-W component
      n(3) = sin(thetan)*sin(phin)

      !     get the Cartesian components of n in the model: nu
      nu_rec(iorientation,1,irec) = n(1)*sint*cosp + n(2)*cost*cosp - n(3)*sinp
      nu_rec(iorientation,2,irec) = n(1)*sint*sinp + n(2)*cost*sinp + n(3)*cosp
      nu_rec(iorientation,3,irec) = n(1)*cost - n(2)*sint
    enddo

    ! point depth (in m)
    depth = stbur(irec)

    ! normalized receiver radius
    r0 = R_UNIT_SPHERE

    ! finds elevation of receiver
    if (TOPOGRAPHY) then
       call get_topo_bathy(lat,lon,elevation,ibathy_topo)
       r0 = r0 + elevation/R_PLANET
    endif

    ! ellipticity
    if (ELLIPTICITY_VAL) then
      ! adds ellipticity factor to radius
      call add_ellipticity_rtheta(r0,theta,nspl_ellip,rspl_ellip,ellipicity_spline,ellipicity_spline2)
    endif

    ! subtract station burial depth (in meters)
    r0 = r0 - depth/R_PLANET

    ! receiver position
    r_target = r0

    ! compute the Cartesian position of the receiver
    x_target = r_target*sint*cosp
    y_target = r_target*sint*sinp
    z_target = r_target*cost

    ! stores Cartesian positions
    xyz_target(1,irec) = x_target
    xyz_target(2,irec) = y_target
    xyz_target(3,irec) = z_target
  enddo

  ! print some information about stations
  if (myrank == 0) then
    ! sorts stations according to epicentral distances
    allocate(irec_dist_ordered(nrec),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating temporary irec_dist_ordered array')

    ! sorts array
    call heap_sort_distances(nrec,epidist,irec_dist_ordered)

    ! outputs info
    write(IMAIN,*) 'Stations sorted by epicentral distance:'
    do i = 1,nrec
      irec = irec_dist_ordered(i)
      write(IMAIN,'(a,i6,a,a24,a,f12.6,a)') ' Station #',irec,': ', &
        trim(network_name(irec))//'.'//trim(station_name(irec)), &
        '    epicentral distance:  ',sngl(epidist(irec)),' degrees'
    enddo

    deallocate(irec_dist_ordered)
  endif


  ! create RECORDHEADERS file with usual format for normal-mode codes
  if (myrank == 0) then

    ! get the base pathname for output files
    call band_instrument_code(DT,bic)

    ! create file for QmX Harvard
    ! Harvard format does not support the network name
    ! therefore only the station name is included below
    ! compute total number of samples for normal modes with 1 sample per second
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/RECORDHEADERS', &
         status='unknown',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file RECORDHEADERS')

    nsamp = nint(dble(NSTEP-1)*DT)

    do irec = 1,nrec

      if (stele(irec) >= -999.9999) then
        write(IOUT,500) station_name(irec),bic(1:2)//'N', &
                     stlat(irec),stlon(irec),stele(irec),stbur(irec), &
                     0.,0.,1.,nsamp,yr,jda,ho,mi,sec
        write(IOUT,500) station_name(irec),bic(1:2)//'E', &
                     stlat(irec),stlon(irec),stele(irec),stbur(irec), &
                     90.,0.,1.,nsamp,yr,jda,ho,mi,sec
        write(IOUT,500) station_name(irec),bic(1:2)//'Z', &
                     stlat(irec),stlon(irec),stele(irec),stbur(irec), &
                     0.,-90.,1.,nsamp,yr,jda,ho,mi,sec
      else
        ! very deep ocean-bottom stations such as H2O are not compatible
        ! with the standard RECORDHEADERS format because of the f6.1 format
        ! therefore suppress decimals for depth in that case
        write(IOUT,600) station_name(irec),bic(1:2)//'N', &
                     stlat(irec),stlon(irec),nint(stele(irec)),stbur(irec), &
                     0.,0.,1.,nsamp,yr,jda,ho,mi,sec
        write(IOUT,600) station_name(irec),bic(1:2)//'E', &
                     stlat(irec),stlon(irec),nint(stele(irec)),stbur(irec), &
                     90.,0.,1.,nsamp,yr,jda,ho,mi,sec
        write(IOUT,600) station_name(irec),bic(1:2)//'Z', &
                     stlat(irec),stlon(irec),nint(stele(irec)),stbur(irec), &
                     0.,-90.,1.,nsamp,yr,jda,ho,mi,sec
      endif
    enddo
    close(IOUT)

  endif

500 format(a8,1x,a3,6x,f9.4,1x,f9.4,1x,f6.1,1x,f9.1,1x,f6.1,1x,f6.1,1x,f12.4,1x,i7,1x,i4.4,1x,i3.3,1x,i2.2,1x,i2.2,1x,f6.3)
600 format(a8,1x,a3,6x,f9.4,1x,f9.4,1x,i6,1x,f9.1,f6.1,1x,f6.1,1x,f12.4,1x,i7,1x,i4.4,1x,i3.3,1x,i2.2,1x,i2.2,1x,f6.3)

  ! make sure we clean the array before the gather
  ispec_selected_rec(:) = 0
  islice_selected_rec(:) = -1

  ! initializes search distances
  final_distance(:) = HUGEVAL
  receiver_final_distance_max = HUGEVAL

  ! loop on all the receivers
  ! gather receiver information in subsets to reduce memory requirements

  ! loop over subsets of receivers
  do irec_already_done = 0, nrec, nrec_SUBSET_MAX

    ! the size of the subset can be the maximum size, or less (if we are in the last subset,
    ! or if there are fewer sources than the maximum size of a subset)
    nrec_SUBSET_current_size = min(nrec_SUBSET_MAX, nrec - irec_already_done)

    ! allocate arrays specific to each subset
    allocate(ispec_selected_subset(nrec_SUBSET_current_size), &
             xi_subset(nrec_SUBSET_current_size), &
             eta_subset(nrec_SUBSET_current_size), &
             gamma_subset(nrec_SUBSET_current_size), &
             xyz_found_subset(NDIM,nrec_SUBSET_current_size), &
             final_distance_subset(nrec_SUBSET_current_size),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary receiver arrays')

    ! gather arrays
    if (myrank == 0) then
      ! only main process needs full arrays allocated
      allocate(ispec_selected_all(nrec_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               xi_all(nrec_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               eta_all(nrec_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               gamma_all(nrec_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               xyz_found_all(NDIM,nrec_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               final_distance_all(nrec_SUBSET_current_size,0:NPROCTOT_VAL-1),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary gather receiver arrays')
    else
      ! dummy arrays
      allocate(ispec_selected_all(1,1), &
               xi_all(1,1), &
               eta_all(1,1), &
               gamma_all(1,1), &
               xyz_found_all(1,1,1), &
               final_distance_all(1,1),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary dummy receiver arrays')
    endif

    ! initializes search results
    ispec_selected_subset(:) = 0
    final_distance_subset(:) = HUGEVAL
    final_distance_all(:,:) = HUGEVAL

    ! find point locations
    ! loop over the stations within this subset
    do irec_in_this_subset = 1,nrec_SUBSET_current_size

      ! mapping from station number in current subset to real station number in all the subsets
      irec = irec_in_this_subset + irec_already_done

      ! station lat/lon in degrees
      lat = stlat(irec)
      lon = stlon(irec)

      ! gets target position
      x_target = xyz_target(1,irec)
      y_target = xyz_target(2,irec)
      z_target = xyz_target(3,irec)

      ! locates best element and xi/eta/gamma interpolation values
      call locate_point(x_target,y_target,z_target,lat,lon,ispec_selected,xi,eta,gamma, &
                        x,y,z,distmin_not_squared,RECEIVERS_CAN_BE_BURIED)

      ! store xi,eta and x,y,z of point found
      xi_subset(irec_in_this_subset) = xi
      eta_subset(irec_in_this_subset) = eta
      gamma_subset(irec_in_this_subset) = gamma

      xyz_found(1,irec) = x
      xyz_found(2,irec) = y
      xyz_found(3,irec) = z

      xyz_found_subset(1,irec_in_this_subset) = x
      xyz_found_subset(2,irec_in_this_subset) = y
      xyz_found_subset(3,irec_in_this_subset) = z

      final_distance(irec) = distmin_not_squared

      final_distance_subset(irec_in_this_subset) = distmin_not_squared
      ispec_selected_subset(irec_in_this_subset) = ispec_selected

    enddo ! end of loop on all stations within current subset

    ! for MPI version, gather information from all the nodes
    ispec_selected_all(:,:) = -1

    call gather_all_i(ispec_selected_subset,nrec_SUBSET_current_size, &
                      ispec_selected_all,nrec_SUBSET_current_size,NPROCTOT_VAL)

    ! this is executed by main process only
    if (myrank == 0) then
      ! check that the gather operation went well
      if (any(ispec_selected_all(:,:) == -1)) then
        print *,'Error ispec all: procs = ',NPROCTOT_VAL,'receivers subset size = ',nrec_SUBSET_current_size
        print *,ispec_selected_all(:,:)
        call exit_MPI(myrank,'gather operation failed for receivers')
      endif
    endif

    call gather_all_dp(xi_subset,nrec_SUBSET_current_size, &
                       xi_all,nrec_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(eta_subset,nrec_SUBSET_current_size, &
                       eta_all,nrec_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(gamma_subset,nrec_SUBSET_current_size, &
                       gamma_all,nrec_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(final_distance_subset,nrec_SUBSET_current_size, &
                       final_distance_all,nrec_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(xyz_found_subset,NDIM*nrec_SUBSET_current_size, &
                       xyz_found_all,NDIM*nrec_SUBSET_current_size,NPROCTOT_VAL)

    ! this is executed by main process only
    if (myrank == 0) then

      ! selects best location in all slices
      ! MPI loop on all the results to determine the best slice
      do irec_in_this_subset = 1,nrec_SUBSET_current_size

        ! mapping from station number in current subset to real station number in all the subsets
        irec = irec_in_this_subset + irec_already_done

        ! loop on all the results to determine the best slice
        distmin_not_squared = HUGEVAL
        do iprocloop = 0,NPROCTOT_VAL-1
          if (final_distance_all(irec_in_this_subset,iprocloop) < distmin_not_squared) then
            ! stores this slice's info
            distmin_not_squared = final_distance_all(irec_in_this_subset,iprocloop)
            islice_selected_rec(irec) = iprocloop
            ispec_selected_rec(irec) = ispec_selected_all(irec_in_this_subset,iprocloop)

            xi_receiver(irec) = xi_all(irec_in_this_subset,iprocloop)
            eta_receiver(irec) = eta_all(irec_in_this_subset,iprocloop)
            gamma_receiver(irec) = gamma_all(irec_in_this_subset,iprocloop)

            xyz_found(:,irec) = xyz_found_all(:,irec_in_this_subset,iprocloop)
          endif
        enddo
        final_distance(irec) = distmin_not_squared
        if (final_distance(irec) == HUGEVAL) call exit_MPI(myrank,'Error locating receiver')

      enddo
    endif ! end of section executed by main process only

    deallocate(ispec_selected_subset)
    deallocate(ispec_selected_all)
    deallocate(xi_subset,eta_subset,gamma_subset)
    deallocate(xi_all,eta_all,gamma_all)
    deallocate(final_distance_all)
    deallocate(final_distance_subset)
    deallocate(xyz_found_subset)
    deallocate(xyz_found_all)

  enddo ! end of loop over all station subsets

  ! deallocate arrays
  deallocate(xyz_target)

  ! this is executed by the main process only
  if (myrank == 0) then

    ! appends receiver locations to sr.vtk file
    open(IOUT_VTK,file=trim(OUTPUT_FILES)//'/sr_tmp.vtk',position='append',status='old',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening and appending receivers to file sr_tmp.vtk')

    ! chooses best receivers locations
    ! if receiver location is too far off, we will exclude the receiver
    islice_selected_found(:) = -1
    nrec_found = 0

    do irec = 1,nrec

      if (DISPLAY_DETAILS_STATIONS .or. final_distance(irec) > 0.01d0 .or. nrec < 50) then

        write(IMAIN,*)
        write(IMAIN,*) 'Station #',irec,': ',trim(network_name(irec))//'.'//trim(station_name(irec))
        write(IMAIN,*) '       original latitude: ',sngl(stlat(irec))
        write(IMAIN,*) '      original longitude: ',sngl(stlon(irec))
        write(IMAIN,*) '     epicentral distance: ',sngl(epidist(irec))
        write(IMAIN,*) '  closest estimate found: ',sngl(final_distance(irec)),' km away'
        write(IMAIN,*) '   in slice ',islice_selected_rec(irec),' in element ',ispec_selected_rec(irec)
        write(IMAIN,*) '   at xi,eta,gamma coordinates = ',xi_receiver(irec),eta_receiver(irec),gamma_receiver(irec)
        write(IMAIN,*) '   at (x,y,z)                  = ',xyz_found(1,irec),xyz_found(2,irec),xyz_found(3,irec)

        ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
        call xyz_2_rlatlon_dble(xyz_found(1,irec),xyz_found(2,irec),xyz_found(3,irec),radius,lat,lon,ELLIPTICITY_VAL)

        ! output same longitude range ([0,360] by default) as input range from stations file stlon(..)
        if (stlon(irec) < 0.d0) lon = lon - 360.d0
        write(IMAIN,*) '   at lat/lon                  = ',sngl(lat),sngl(lon)
      endif

      ! add warning if estimate is poor
      ! (usually means receiver outside the mesh given by the user)
      if (final_distance(irec) > THRESHOLD_EXCLUDE_STATION) then
        write(IMAIN,*) 'Station #',irec,': ',trim(network_name(irec))//'.'//trim(station_name(irec))
        write(IMAIN,*) '*****************************************************************'
        if (NCHUNKS_VAL == 6) then
          write(IMAIN,*) '***** WARNING: receiver location estimate is poor, therefore receiver excluded *****'
        else
          write(IMAIN,*) '***** WARNING: receiver is located outside the mesh, therefore excluded *****'
        endif
        write(IMAIN,*) '*****************************************************************'
      else
        nrec_found = nrec_found + 1

        islice_selected_found(nrec_found) = islice_selected_rec(irec)
        ispec_selected_found(nrec_found) = ispec_selected_rec(irec)

        xi_receiver_found(nrec_found) = xi_receiver(irec)
        eta_receiver_found(nrec_found) = eta_receiver(irec)
        gamma_receiver_found(nrec_found) = gamma_receiver(irec)

        station_name_found(nrec_found) = station_name(irec)
        network_name_found(nrec_found) = network_name(irec)

        stlat_found(nrec_found) = stlat(irec)
        stlon_found(nrec_found) = stlon(irec)
        stele_found(nrec_found) = stele(irec)
        stbur_found(nrec_found) = stbur(irec)

        nu_found(:,:,nrec_found) = nu_rec(:,:,irec)
        epidist_found(nrec_found) = epidist(irec)

        ! writes out actual receiver location to VTK file
        write(IOUT_VTK,'(3e18.6)') sngl(xyz_found(1,irec)), sngl(xyz_found(2,irec)), sngl(xyz_found(3,irec))
      endif
    enddo

    ! finishes sr_tmp.vtk file
    write(IOUT_VTK,*)
    close(IOUT_VTK)

    ! compute maximal distance for all the receivers
    receiver_final_distance_max = maxval(final_distance(:))

    ! display maximum error for all the receivers
    write(IMAIN,*)
    write(IMAIN,*) 'maximum error in location of all the receivers: ',sngl(receiver_final_distance_max),' km'

    ! add warning if estimate is poor
    ! (usually means receiver outside the mesh given by the user)
    if (receiver_final_distance_max > THRESHOLD_EXCLUDE_STATION) then
      write(IMAIN,*)
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '***** WARNING: at least one receiver was excluded from the station list *****'
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '************************************************************'
    endif
    call flush_IMAIN()

    ! replaces list of all receivers with list of only those which were found
    ! (in particular: for a 1-chunk simulation, only stations in this chunk)
    nrec = nrec_found

    islice_selected_rec(1:nrec) = islice_selected_found(1:nrec)
    ispec_selected_rec(1:nrec) = ispec_selected_found(1:nrec)

    xi_receiver(1:nrec) = xi_receiver_found(1:nrec)
    eta_receiver(1:nrec) = eta_receiver_found(1:nrec)
    gamma_receiver(1:nrec) = gamma_receiver_found(1:nrec)

    station_name(1:nrec) = station_name_found(1:nrec)
    network_name(1:nrec) = network_name_found(1:nrec)

    stlat(1:nrec) = stlat_found(1:nrec)
    stlon(1:nrec) = stlon_found(1:nrec)
    stele(1:nrec) = stele_found(1:nrec)
    stbur(1:nrec) = stbur_found(1:nrec)

    nu_rec(:,:,1:nrec) = nu_found(:,:,1:nrec)
    epidist(1:nrec) = epidist_found(1:nrec)

    ! write the list of stations and associated epicentral distance
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/output_list_stations.txt', &
          status='unknown',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file output_list_stations.txt')
    write(IOUT,*)
    write(IOUT,*) 'total number of stations (found): ',nrec
    write(IOUT,*)
    do irec = 1,nrec
      write(IOUT,*) &
        network_name(irec)(1:len_trim(network_name(irec))),'.',station_name(irec)(1:len_trim(station_name(irec))), &
        ' epicentral distance ',sngl(epidist(irec)),' deg'
    enddo
    close(IOUT)

    ! write out a filtered station list
    if (NCHUNKS_VAL /= 6) then
      open(unit=IOUT,file=trim(OUTPUT_FILES)//'/STATIONS_FILTERED', &
            status='unknown',iostat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening file STATIONS_FILTERED')
      ! loop on all the stations to read station information
      do irec = 1,nrec
        write(IOUT,'(a8,1x,a3,6x,f8.4,1x,f9.4,1x,f9.1,1x,f9.1)') trim(station_name(irec)),trim(network_name(irec)), &
          sngl(stlat(irec)),sngl(stlon(irec)),sngl(stele(irec)),sngl(stbur(irec))
      enddo
      ! close receiver file
      close(IOUT)
    endif

  endif    ! end of section executed by main process only

  ! main process broadcasts the results to all the slices
  call bcast_all_singlei(nrec)

  call bcast_all_i(islice_selected_rec,nrec)
  call bcast_all_i(ispec_selected_rec,nrec)

  call bcast_all_dp(xi_receiver,nrec)
  call bcast_all_dp(eta_receiver,nrec)
  call bcast_all_dp(gamma_receiver,nrec)

  call bcast_all_ch_array(station_name,nrec,MAX_LENGTH_STATION_NAME)
  call bcast_all_ch_array(network_name,nrec,MAX_LENGTH_NETWORK_NAME)

  call bcast_all_dp(stlat,nrec)
  call bcast_all_dp(stlon,nrec)
  call bcast_all_dp(stele,nrec)
  call bcast_all_dp(stbur,nrec)
  call bcast_all_dp(nu_rec,nrec*3*3)

  ! deallocate arrays
  deallocate(epidist)
  deallocate(xyz_found)
  deallocate(final_distance)

  ! synchronizes to get right timing
  call synchronize_all()

  ! elapsed time since beginning of mesh generation
  if (myrank == 0) then
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for receiver detection in seconds = ',tCPU
    write(IMAIN,*)
    write(IMAIN,*) 'End of receiver detection - done'
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  end subroutine locate_receivers

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_receiver_locations()

  use constants_solver, only: &
    MAX_LENGTH_STATION_NAME,MAX_LENGTH_NETWORK_NAME,IMAIN,IIN

  use specfem_par, only: &
    myrank,nrec,STATIONS_FILE,RECEIVERS_CAN_BE_BURIED, &
    station_name,network_name,stlat,stlon,stele,stbur

  implicit none

  ! local parameters
  integer :: ier,irec,i
  integer, allocatable, dimension(:) :: station_duplet
  character(len=256) :: line

  ! reads stations files
  if (myrank == 0) then
    ! user output
    write(IMAIN,*) 'reading receiver information...'
    write(IMAIN,*)
    call flush_IMAIN()

    ! opens station file STATIONS or STATIONS_ADJOINT
    open(unit=IIN,file=trim(STATIONS_FILE),status='old',action='read',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening STATIONS file')

    ! loop on all the stations to read station information
    do irec = 1,nrec

      ! old line:
      !read(IIN,*,iostat=ier) station_name(irec),network_name(irec),stlat(irec),stlon(irec),stele(irec),stbur(irec)

      ! reads in line as string
      read(IIN,"(a256)",iostat=ier) line
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading in station ',irec
        call exit_MPI(myrank,'Error reading in station in STATIONS file')
      endif

      ! skips empty lines and comment lines
      do while( len_trim(line) == 0 .or. line(1:1) == '#')
        read(IIN,"(a256)",iostat=ier) line
        if (ier /= 0) then
          write(IMAIN,*) 'Error reading in station ',irec
          call exit_MPI(myrank,'Error reading in station in STATIONS file')
        endif
      enddo

      ! reads in station information
      read(line(1:len_trim(line)),*,iostat=ier) station_name(irec),network_name(irec), &
                                                stlat(irec),stlon(irec),stele(irec),stbur(irec)
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading in station ',irec
        call exit_MPI(myrank,'Error reading in station in STATIONS file')
      endif

      ! checks latitude
      if (stlat(irec) < -90.d0 .or. stlat(irec) > 90.d0) then
        write(IMAIN,*) 'Error station ',trim(station_name(irec)),': latitude ',stlat(irec), &
                       ' is invalid, please check STATIONS record'
        call exit_MPI(myrank,'Error station latitude invalid')
      endif

    enddo
    ! close receiver file
    close(IIN)

    ! In case that the same station and network name appear twice (or more times) in the STATIONS
    ! file, problems occur, as two (or more) seismograms are written (with mode
    ! "append") to a file with same name. The philosophy here is to accept multiple
    ! appearances and to just add a count to the station name in this case.
    allocate(station_duplet(nrec),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating station_duplet array')

    station_duplet(:) = 0
    do irec = 1,nrec
      do i = 1,irec-1
        if ((station_name(irec) == station_name(i)) .and. &
            (network_name(irec) == network_name(i))) then

            station_duplet(i)=station_duplet(i)+1
            if (len_trim(station_name(irec)) <= MAX_LENGTH_STATION_NAME-3) then
              write(station_name(irec),"(a,'_',i2.2)") trim(station_name(irec)),station_duplet(i)+1
            else
              call exit_MPI(myrank,'Please increase MAX_LENGTH_STATION_NAME by at least 3 to name station duplets')
            endif

        endif
      enddo
    enddo
    deallocate(station_duplet)

    ! if receivers can not be buried, sets depth to zero
    if (.not. RECEIVERS_CAN_BE_BURIED ) stbur(:) = 0.d0

  endif

  ! broadcast the information read on the main node to all the nodes
  call bcast_all_ch_array(station_name,nrec,MAX_LENGTH_STATION_NAME)
  call bcast_all_ch_array(network_name,nrec,MAX_LENGTH_NETWORK_NAME)
  call bcast_all_dp(stlat,nrec)
  call bcast_all_dp(stlon,nrec)
  call bcast_all_dp(stele,nrec)
  call bcast_all_dp(stbur,nrec)

  end subroutine read_receiver_locations

!
!-------------------------------------------------------------------------------------------------
!

! sorting routine left here for inlining
!
! Implementation of a Heap Sort Routine
!    Input
!      n = Input
!         Length of arrays
!      X_in = Input
!             Vector to be sorted
!             dimension(n)
!      Y = Output
!         Sorted Indices of vector X
!
!      Example:
!         D = [ 4.0 3.0 1.0 2.0 ] on Input
!         Y = [ 1 2 3 4 ] Computed Internally (in order)
!
!         X = [ 1.0 2.0 3.0 4.0 ] Computed Internally
!         Y = [ 3 4 2 1 ] on Output
!
  subroutine heap_sort_distances(N, X_in, Y)

  implicit none
  integer, intent(in) :: N
  double precision, dimension(N), intent(in) :: X_in
  integer, dimension(N), intent(out) :: Y

  ! local parameters
  double precision, dimension(N) :: X
  double precision :: tmp
  integer :: itmp
  integer :: i

  do i = 1,N
     Y(i) = i
     X(i) = X_in(i)
  enddo

  ! checks if anything to do
  if (N < 2) return

  ! builds heap
  do i = N/2, 1, -1
    call my_heap_sort_siftdown(i, N)
  enddo

  ! sorts array
  do i = N, 2, -1
    ! swaps last and first entry in this section
    tmp = X(1)
    X(1) = X(i)
    X(i) = tmp
    itmp = Y(1)
    Y(1) = Y(i)
    Y(i) = itmp

    call my_heap_sort_siftdown(1, i - 1)
  enddo

  contains

    subroutine my_heap_sort_siftdown(start, bottom)

    implicit none

    integer, intent(in) :: start, bottom

    ! local parameters
    integer :: i, j
    double precision :: xtmp
    integer :: ytmp

    i = start
    xtmp = X(i)
    ytmp = Y(i)

    j = 2 * i
    do while (j <= bottom)
      ! chooses larger value first in this section
      if (j < bottom) then
        if (X(j) <= X(j+1)) j = j + 1
      endif

      ! checks if section already smaller than initial value
      if (X(j) < xtmp) exit

      X(i) = X(j)
      Y(i) = Y(j)
      i = j
      j = 2 * i
    enddo

    X(i) = xtmp
    Y(i) = ytmp

    end subroutine my_heap_sort_siftdown

  end subroutine heap_sort_distances

