!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
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

  subroutine locate_receivers(nspec,nglob,ibool, &
                             xstore,ystore,zstore, &
                             yr,jda,ho,mi,sec, &
                             theta_source,phi_source)

  use constants_solver, only: &
    ELLIPTICITY_VAL,NCHUNKS_VAL,NEX_XI_VAL,NPROCTOT_VAL, &
    CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGNOD,NDIM, &
    MAX_LENGTH_STATION_NAME,MAX_LENGTH_NETWORK_NAME, &
    USE_DISTANCE_CRITERION,DISPLAY_DETAILS_STATIONS,nrec_SUBSET_MAX, &
    THRESHOLD_EXCLUDE_STATION,NUM_ITER, &
    HUGEVAL,IMAIN,IIN,IOUT,IOUT_VTK,MIDX,MIDY,MIDZ, &
    DEGREES_TO_RADIANS,RADIANS_TO_DEGREES,TWO_PI,R_UNIT_SPHERE,R_EARTH,R_EARTH_KM

  use shared_input_parameters, only: OUTPUT_FILES

  use specfem_par, only: &
    myrank,DT,NSTEP, &
    xigll,yigll,zigll, &
    STATIONS_FILE,nrec,islice_selected_rec,ispec_selected_rec, &
    xi_receiver,eta_receiver,gamma_receiver,station_name,network_name, &
    stlat,stlon,stele,stbur,nu, &
    rspl,espl,espl2,nspl,ibathy_topo, &
    TOPOGRAPHY,RECEIVERS_CAN_BE_BURIED

  implicit none

  integer,intent(in) :: nspec,nglob

  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nspec)

  ! arrays containing coordinates of the points
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore,ystore,zstore

  integer,intent(in) :: yr,jda,ho,mi
  double precision,intent(in) :: sec
  double precision,intent(in) :: theta_source,phi_source

  ! local parameters
  integer :: nrec_found
  integer, allocatable, dimension(:) :: ix_initial_guess,iy_initial_guess,iz_initial_guess

  integer :: iorientation
  integer :: iprocloop
  double precision :: stazi,stdip

  double precision, allocatable, dimension(:) :: x_target,y_target,z_target
  double precision, allocatable, dimension(:) :: epidist
  double precision, allocatable, dimension(:) :: x_found,y_found,z_found
  double precision, allocatable, dimension(:,:) :: x_found_all,y_found_all,z_found_all

  integer :: irec
  integer :: i,j,k,ispec,iglob
  integer :: ier

  double precision :: ell
  double precision :: elevation
  double precision :: n(3)
  double precision :: thetan,phin
  double precision :: sint,cost,sinp,cosp
  double precision :: r0,p20
  double precision :: theta,phi
  double precision :: dist_squared
  double precision :: xi,eta,gamma,dx,dy,dz,dxi,deta,dgamma

! topology of the control points of the surface element
  integer :: iax,iay,iaz
  integer :: iaddx(NGNOD),iaddy(NGNOD),iaddr(NGNOD)

! coordinates of the control points of the surface element
  double precision :: xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  integer :: iter_loop,ispec_iterate

  integer :: ia
  double precision :: x,y,z
  double precision :: xix,xiy,xiz
  double precision :: etax,etay,etaz
  double precision :: gammax,gammay,gammaz

  ! timer MPI
  double precision :: time_start,tCPU

  ! use dynamic allocation
  double precision, dimension(:), allocatable :: final_distance
  double precision, dimension(:,:), allocatable :: final_distance_all

  double precision :: distmin_squared,distmin_not_squared,final_distance_max

! receiver information
! timing information for the stations
! station information for writing the seismograms
  integer :: nsamp

  integer, dimension(nrec) :: islice_selected_rec_found,ispec_selected_rec_found
  double precision, dimension(nrec) :: xi_receiver_found,eta_receiver_found,gamma_receiver_found
  double precision, dimension(3,3,nrec) :: nu_found
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name_found
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name_found
  double precision, dimension(nrec) :: stlat_found,stlon_found,stele_found,stbur_found,epidist_found

  integer, allocatable, dimension(:,:) :: ispec_selected_rec_all
  double precision, allocatable, dimension(:,:) :: xi_receiver_all,eta_receiver_all,gamma_receiver_all
  double precision, allocatable, dimension(:) :: xi_receiver_subset,eta_receiver_subset,gamma_receiver_subset

  double precision :: x_target_rec,y_target_rec,z_target_rec

  double precision :: typical_size_squared

  character(len=2) :: bic
  character(len=256) :: string

  integer :: nrec_SUBSET_current_size,irec_in_this_subset,irec_already_done
  double precision, allocatable, dimension(:) :: x_found_subset,y_found_subset,z_found_subset
  double precision, allocatable, dimension(:) :: final_distance_subset
  integer, allocatable, dimension(:) :: ispec_selected_rec_subset

  integer, allocatable, dimension(:) :: station_duplet

  ! timing
  double precision, external :: wtime

  ! distance
  double precision, allocatable, dimension(:,:) :: xyz_midpoints

  ! search range
  double precision :: lat,lon,r
  double precision :: lat_min,lat_max,lon_min,lon_max
  ! search margin in degrees
  double precision,parameter :: LAT_LON_MARGIN = 2.d0
  ! sorting order
  integer, allocatable, dimension(:) :: irec_dist_ordered

  ! get MPI starting time
  time_start = wtime()

  ! make sure we clean the array before the gather
  ispec_selected_rec(:) = 0

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' locating receivers'
    write(IMAIN,*) '********************'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! define topology of the control element
  call hex_nodes(iaddx,iaddy,iaddr)

  ! compute typical size of elements at the surface
  typical_size_squared = TWO_PI * R_UNIT_SPHERE / (4.*NEX_XI_VAL)

  ! use 10 times the distance as a criterion for source detection
  typical_size_squared = (10. * typical_size_squared)**2

  ! allocate memory for arrays using number of stations
  allocate(epidist(nrec), &
           ix_initial_guess(nrec), &
           iy_initial_guess(nrec), &
           iz_initial_guess(nrec), &
           x_target(nrec), &
           y_target(nrec), &
           z_target(nrec), &
           x_found(nrec), &
           y_found(nrec), &
           z_found(nrec), &
           final_distance(nrec),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating temporary receiver arrays')

  ! read that STATIONS file on the master
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
      read(IIN,"(a256)",iostat=ier) string
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading in station ',irec
        call exit_MPI(myrank,'Error reading in station in STATIONS file')
      endif

      ! skips empty lines
      do while( len_trim(string) == 0 )
        read(IIN,"(a256)",iostat=ier) string
        if (ier /= 0) then
          write(IMAIN,*) 'Error reading in station ',irec
          call exit_MPI(myrank,'Error reading in station in STATIONS file')
        endif
      enddo

      ! reads in station information
      read(string(1:len_trim(string)),*,iostat=ier) station_name(irec),network_name(irec), &
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

! BS BS begin
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
! BS BS end

    ! if receivers can not be buried, sets depth to zero
    if (.not. RECEIVERS_CAN_BE_BURIED ) stbur(:) = 0.d0

  endif

  ! broadcast the information read on the master to the nodes
  call bcast_all_ch_array(station_name,nrec,MAX_LENGTH_STATION_NAME)
  call bcast_all_ch_array(network_name,nrec,MAX_LENGTH_NETWORK_NAME)
  call bcast_all_dp(stlat,nrec)
  call bcast_all_dp(stlon,nrec)
  call bcast_all_dp(stele,nrec)
  call bcast_all_dp(stbur,nrec)

  ! limits receiver search
  if (USE_DISTANCE_CRITERION) then
    ! retrieves latitude/longitude range of this slice
    call xyz_2_latlon_minmax(nspec,nglob,ibool,xstore,ystore,zstore,lat_min,lat_max,lon_min,lon_max)

    ! adds search margin
    lat_min = lat_min - LAT_LON_MARGIN
    lat_max = lat_max + LAT_LON_MARGIN

    lon_min = lon_min - LAT_LON_MARGIN
    lon_max = lon_max + LAT_LON_MARGIN

    ! limits latitude to [-90.0,90.0]
    if (lat_min < -90.d0 ) lat_min = -90.d0
    if (lat_max > 90.d0 ) lat_max = 90.d0

    ! limits longitude to [0.0,360.0]
    if (lon_min < 0.d0 ) lon_min = 0.d0
    if (lon_min > 360.d0 ) lon_min = 360.d0
    if (lon_max < 0.d0 ) lon_max = 0.d0
    if (lon_max > 360.d0 ) lon_max = 360.d0

    ! prepares midpoints coordinates
    allocate(xyz_midpoints(NDIM,nspec),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating array xyz_midpoints')
    ! store x/y/z coordinates of center point
    do ispec = 1,nspec
      iglob = ibool(MIDX,MIDY,MIDZ,ispec)
      xyz_midpoints(1,ispec) =  dble(xstore(iglob))
      xyz_midpoints(2,ispec) =  dble(ystore(iglob))
      xyz_midpoints(3,ispec) =  dble(zstore(iglob))
    enddo
  endif

  ! loop on all the stations to locate them in the mesh
  do irec = 1,nrec

    ! set distance to huge initial value
    distmin_squared = HUGEVAL

    ! station lat/lon in degrees
    lat = stlat(irec)
    lon = stlon(irec)

    ! limits longitude to [0.0,360.0]
    if (lon < 0.d0 ) lon = lon + 360.d0
    if (lon > 360.d0 ) lon = lon - 360.d0

    ! converts geographic latitude stlat (degrees) to geocentric colatitude theta (radians)
    call lat_2_geocentric_colat_dble(lat,theta)

    phi = lon*DEGREES_TO_RADIANS
    call reduce(theta,phi)

    ! compute epicentral distance
    epidist(irec) = acos(cos(theta)*cos(theta_source) + &
                         sin(theta)*sin(theta_source)*cos(phi-phi_source))*RADIANS_TO_DEGREES

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

      !     get the orientation of the seismometer
      thetan=(90.0d0+stdip)*DEGREES_TO_RADIANS
      phin=stazi*DEGREES_TO_RADIANS

      ! we use the same convention as in Harvard normal modes for the orientation

      !     vertical component
      n(1) = cos(thetan)
      !     N-S component
      n(2) = - sin(thetan)*cos(phin)
      !     E-W component
      n(3) = sin(thetan)*sin(phin)

      !     get the Cartesian components of n in the model: nu
      sint = sin(theta)
      cost = cos(theta)
      sinp = sin(phi)
      cosp = cos(phi)
      nu(iorientation,1,irec) = n(1)*sint*cosp+n(2)*cost*cosp-n(3)*sinp
      nu(iorientation,2,irec) = n(1)*sint*sinp+n(2)*cost*sinp+n(3)*cosp
      nu(iorientation,3,irec) = n(1)*cost-n(2)*sint
    enddo

    ! normalized receiver radius
    r0 = R_UNIT_SPHERE

    ! finds elevation of receiver
    if (TOPOGRAPHY) then
       call get_topo_bathy(lat,lon,elevation,ibathy_topo)
       r0 = r0 + elevation/R_EARTH
    endif
    ! ellipticity
    if (ELLIPTICITY_VAL) then
      cost=cos(theta)
      ! this is the Legendre polynomial of degree two, P2(cos(theta)),
      ! see the discussion above eq (14.4) in Dahlen and Tromp (1998)
      p20=0.5d0*(3.0d0*cost*cost-1.0d0)
      ! get ellipticity using spline evaluation
      call spline_evaluation(rspl,espl,espl2,nspl,r0,ell)
      ! this is eq (14.4) in Dahlen and Tromp (1998)
      r0=r0*(1.0d0-(2.0d0/3.0d0)*ell*p20)
    endif

    ! subtract station burial depth (in meters)
    r0 = r0 - stbur(irec)/R_EARTH

    ! compute the Cartesian position of the receiver
    x_target_rec = r0*sin(theta)*cos(phi)
    y_target_rec = r0*sin(theta)*sin(phi)
    z_target_rec = r0*cos(theta)

    x_target(irec) = x_target_rec
    y_target(irec) = y_target_rec
    z_target(irec) = z_target_rec

    ! would write out desired target locations of receivers
    !if (myrank == 0) write(IOUT_VTK,*) sngl(x_target(irec)), sngl(y_target(irec)), sngl(z_target(irec))

    ! initializes located target
    ! if we have not located a target element, the receiver is not in this slice
    ! therefore use first element only for fictitious iterative search
    ispec_selected_rec(irec) = 1
    ix_initial_guess(irec) = MIDX
    iy_initial_guess(irec) = MIDY
    iz_initial_guess(irec) = MIDZ
    distmin_squared = HUGEVAL

    ! searches closest GLL point
    if (USE_DISTANCE_CRITERION) then
      ! checks if receiver in this slice
      if (lat >= lat_min .and. lat <= lat_max .and. &
          lon >= lon_min .and. lon <= lon_max) then
        ! loops over all elements
        do ispec = 1,nspec
          ! exclude elements that are too far from target
          dist_squared = (x_target_rec - xyz_midpoints(1,ispec))*(x_target_rec - xyz_midpoints(1,ispec)) &
                       + (y_target_rec - xyz_midpoints(2,ispec))*(y_target_rec - xyz_midpoints(2,ispec)) &
                       + (z_target_rec - xyz_midpoints(3,ispec))*(z_target_rec - xyz_midpoints(3,ispec))
          !  we compare squared distances instead of distances themselves to significantly speed up calculations
          if (dist_squared > typical_size_squared) cycle
          ! loop only on points inside the element
          ! exclude edges to ensure this point is not shared with other elements
          do k=2,NGLLZ-1
            do j=2,NGLLY-1
              do i=2,NGLLX-1
                iglob = ibool(i,j,k,ispec)
                dist_squared = (x_target_rec - dble(xstore(iglob)))*(x_target_rec - dble(xstore(iglob))) &
                             + (y_target_rec - dble(ystore(iglob)))*(y_target_rec - dble(ystore(iglob))) &
                             + (z_target_rec - dble(zstore(iglob)))*(z_target_rec - dble(zstore(iglob)))
                !  keep this point if it is closer to the receiver
                !  we compare squared distances instead of distances themselves to significantly speed up calculations
                if (dist_squared < distmin_squared) then
                  distmin_squared = dist_squared
                  ispec_selected_rec(irec) = ispec
                  ix_initial_guess(irec) = i
                  iy_initial_guess(irec) = j
                  iz_initial_guess(irec) = k
                endif
              enddo
            enddo
          enddo
        ! end of loop on all the spectral elements in current slice
        enddo
      endif
    else
      ! searches through all elements
      do ispec = 1,nspec
        ! loop only on points inside the element
        ! exclude edges to ensure this point is not shared with other elements
        do k=2,NGLLZ-1
          do j=2,NGLLY-1
            do i=2,NGLLX-1
              iglob = ibool(i,j,k,ispec)
              dist_squared = (x_target_rec - dble(xstore(iglob)))*(x_target_rec - dble(xstore(iglob))) &
                           + (y_target_rec - dble(ystore(iglob)))*(y_target_rec - dble(ystore(iglob))) &
                           + (z_target_rec - dble(zstore(iglob)))*(z_target_rec - dble(zstore(iglob)))
              !  keep this point if it is closer to the receiver
              !  we compare squared distances instead of distances themselves to significantly speed up calculations
              if (dist_squared < distmin_squared) then
                distmin_squared = dist_squared
                ispec_selected_rec(irec) = ispec
                ix_initial_guess(irec) = i
                iy_initial_guess(irec) = j
                iz_initial_guess(irec) = k
              endif
            enddo
          enddo
        enddo
      ! end of loop on all the spectral elements in current slice
      enddo
    endif ! USE_DISTANCE_CRITERION

  ! end of loop on all the stations
  enddo

  if (USE_DISTANCE_CRITERION ) deallocate(xyz_midpoints)

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

500 format(a8,1x,a3,6x,f9.4,1x,f9.4,1x,f6.1,1x,f6.1,1x,f6.1,1x,f6.1,1x,f12.4,1x,i7,1x,i4.4,1x,i3.3,1x,i2.2,1x,i2.2,1x,f6.3)
600 format(a8,1x,a3,6x,f9.4,1x,f9.4,1x,i6,1x,f6.1,f6.1,1x,f6.1,1x,f12.4,1x,i7,1x,i4.4,1x,i3.3,1x,i2.2,1x,i2.2,1x,f6.3)


! ****************************************
! find the best (xi,eta) for each receiver
! ****************************************

  ! initializes search distances
  final_distance(:) = HUGEVAL

  ! loop on all the receivers
  ! gather source information in subsets to reduce memory requirements
  islice_selected_rec(:) = -1

  ! loop over subsets of receivers
  do irec_already_done = 0, nrec, nrec_SUBSET_MAX

    ! the size of the subset can be the maximum size, or less (if we are in the last subset,
    ! or if there are fewer sources than the maximum size of a subset)
    nrec_SUBSET_current_size = min(nrec_SUBSET_MAX, nrec - irec_already_done)

    ! allocate arrays specific to each subset
    allocate(ispec_selected_rec_subset(nrec_SUBSET_current_size), &
             xi_receiver_subset(nrec_SUBSET_current_size), &
             eta_receiver_subset(nrec_SUBSET_current_size), &
             gamma_receiver_subset(nrec_SUBSET_current_size), &
             x_found_subset(nrec_SUBSET_current_size), &
             y_found_subset(nrec_SUBSET_current_size), &
             z_found_subset(nrec_SUBSET_current_size), &
             final_distance_subset(nrec_SUBSET_current_size),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary receiver arrays')

    ! gather arrays
    if (myrank == 0) then
      ! only master process needs full arrays allocated
      allocate(ispec_selected_rec_all(nrec_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               xi_receiver_all(nrec_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               eta_receiver_all(nrec_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               gamma_receiver_all(nrec_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               x_found_all(nrec_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               y_found_all(nrec_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               z_found_all(nrec_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               final_distance_all(nrec_SUBSET_current_size,0:NPROCTOT_VAL-1),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary gather receiver arrays')
    else
      ! dummy arrays
      allocate(ispec_selected_rec_all(1,1), &
               xi_receiver_all(1,1), &
               eta_receiver_all(1,1), &
               gamma_receiver_all(1,1), &
               x_found_all(1,1), &
               y_found_all(1,1), &
               z_found_all(1,1), &
               final_distance_all(1,1),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary dummy receiver arrays')
    endif

    ! initializes search results
    final_distance_all(:,:) = HUGEVAL
    final_distance_subset(:) = HUGEVAL

    ! loop over the stations within this subset
    do irec_in_this_subset = 1,nrec_SUBSET_current_size

      ! mapping from station number in current subset to real station number in all the subsets
      irec = irec_in_this_subset + irec_already_done

      ispec_iterate = ispec_selected_rec(irec)
      ispec_selected_rec_subset(irec_in_this_subset) = ispec_iterate

      ! define coordinates of the control points of the element
      do ia = 1,NGNOD

        iax = 0
        if (iaddx(ia) == 0) then
          iax = 1
        else if (iaddx(ia) == 1) then
          iax = MIDX
        else if (iaddx(ia) == 2) then
          iax = NGLLX
        else
          call exit_MPI(myrank,'incorrect value of iaddx')
        endif

        iay = 0
        if (iaddy(ia) == 0) then
          iay = 1
        else if (iaddy(ia) == 1) then
          iay = MIDY
        else if (iaddy(ia) == 2) then
          iay = NGLLY
        else
          call exit_MPI(myrank,'incorrect value of iaddy')
        endif

        iaz = 0
        if (iaddr(ia) == 0) then
          iaz = 1
        else if (iaddr(ia) == 1) then
          iaz = MIDZ
        else if (iaddr(ia) == 2) then
          iaz = NGLLZ
        else
          call exit_MPI(myrank,'incorrect value of iaddr')
        endif

        iglob = ibool(iax,iay,iaz,ispec_iterate)
        xelm(ia) = dble(xstore(iglob))
        yelm(ia) = dble(ystore(iglob))
        zelm(ia) = dble(zstore(iglob))

      enddo

      ! use initial guess in xi and eta
      xi = xigll(ix_initial_guess(irec))
      eta = yigll(iy_initial_guess(irec))
      gamma = zigll(iz_initial_guess(irec))

      ! impose receiver exactly at the surface
      if (.not. RECEIVERS_CAN_BE_BURIED) gamma = 1.d0

      x_target_rec = x_target(irec)
      y_target_rec = y_target(irec)
      z_target_rec = z_target(irec)

      ! iterate to solve the non linear system
      do iter_loop = 1,NUM_ITER

        ! recompute Jacobian for the new point
        call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                                xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

        ! compute distance to target location
        dx = - (x - x_target_rec)
        dy = - (y - y_target_rec)
        dz = - (z - z_target_rec)

        ! compute increments
        ! gamma does not change since we know the receiver is exactly on the surface
        dxi  = xix*dx + xiy*dy + xiz*dz
        deta = etax*dx + etay*dy + etaz*dz

        ! update values
        xi = xi + dxi
        eta = eta + deta

        ! buried receivers vary in z depth
        if (RECEIVERS_CAN_BE_BURIED) then
          dgamma = gammax*dx + gammay*dy + gammaz*dz
          gamma = gamma + dgamma
        endif

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
      if (.not. RECEIVERS_CAN_BE_BURIED) gamma = 1.d0

      ! compute final coordinates of point found
      call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                              xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

      ! store xi,eta and x,y,z of point found
      xi_receiver_subset(irec_in_this_subset) = xi
      eta_receiver_subset(irec_in_this_subset) = eta
      gamma_receiver_subset(irec_in_this_subset) = gamma

      x_found(irec) = x
      y_found(irec) = y
      z_found(irec) = z

      x_found_subset(irec_in_this_subset) = x_found(irec)
      y_found_subset(irec_in_this_subset) = y_found(irec)
      z_found_subset(irec_in_this_subset) = z_found(irec)

      ! compute final distance between asked and found (converted to km)
      distmin_not_squared = dsqrt((x_target_rec-x)**2 + &
                                  (y_target_rec-y)**2 + &
                                  (z_target_rec-z)**2)*R_EARTH_KM

      final_distance(irec) = distmin_not_squared
      final_distance_subset(irec_in_this_subset) = distmin_not_squared

    enddo ! end of loop on all stations within current subset

    ! for MPI version, gather information from all the nodes
    ispec_selected_rec_all(:,:) = -1

    call gather_all_i(ispec_selected_rec_subset,nrec_SUBSET_current_size, &
                      ispec_selected_rec_all,nrec_SUBSET_current_size,NPROCTOT_VAL)

    ! this is executed by main process only
    if (myrank == 0) then
      ! check that the gather operation went well
      if (any(ispec_selected_rec_all(:,:) == -1)) call exit_MPI(myrank,'gather operation failed for receivers')
    endif

    call gather_all_dp(xi_receiver_subset,nrec_SUBSET_current_size,xi_receiver_all,nrec_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(eta_receiver_subset,nrec_SUBSET_current_size,eta_receiver_all,nrec_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(gamma_receiver_subset,nrec_SUBSET_current_size,gamma_receiver_all,nrec_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(final_distance_subset,nrec_SUBSET_current_size,final_distance_all,nrec_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(x_found_subset,nrec_SUBSET_current_size,x_found_all,nrec_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(y_found_subset,nrec_SUBSET_current_size,y_found_all,nrec_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(z_found_subset,nrec_SUBSET_current_size,z_found_all,nrec_SUBSET_current_size,NPROCTOT_VAL)

    ! this is executed by main process only
    if (myrank == 0) then

      ! MPI loop on all the results to determine the best slice
      do irec_in_this_subset = 1,nrec_SUBSET_current_size

        ! mapping from station number in current subset to real station number in all the subsets
        irec = irec_in_this_subset + irec_already_done

        distmin_not_squared = HUGEVAL
        do iprocloop = 0,NPROCTOT_VAL-1
          if (final_distance_all(irec_in_this_subset,iprocloop) < distmin_not_squared) then
            distmin_not_squared = final_distance_all(irec_in_this_subset,iprocloop)

            islice_selected_rec(irec) = iprocloop
            ispec_selected_rec(irec) = ispec_selected_rec_all(irec_in_this_subset,iprocloop)
            xi_receiver(irec) = xi_receiver_all(irec_in_this_subset,iprocloop)
            eta_receiver(irec) = eta_receiver_all(irec_in_this_subset,iprocloop)
            gamma_receiver(irec) = gamma_receiver_all(irec_in_this_subset,iprocloop)
            x_found(irec) = x_found_all(irec_in_this_subset,iprocloop)
            y_found(irec) = y_found_all(irec_in_this_subset,iprocloop)
            z_found(irec) = z_found_all(irec_in_this_subset,iprocloop)
          endif
        enddo
        final_distance(irec) = distmin_not_squared
      enddo
    endif ! end of section executed by main process only

    deallocate(ispec_selected_rec_subset)
    deallocate(ispec_selected_rec_all)
    deallocate(xi_receiver_subset)
    deallocate(eta_receiver_subset)
    deallocate(gamma_receiver_subset)
    deallocate(xi_receiver_all)
    deallocate(eta_receiver_all)
    deallocate(gamma_receiver_all)
    deallocate(x_found_subset)
    deallocate(y_found_subset)
    deallocate(z_found_subset)
    deallocate(x_found_all)
    deallocate(y_found_all)
    deallocate(z_found_all)
    deallocate(final_distance_subset)
    deallocate(final_distance_all)

  enddo ! end of loop over all station subsets

  ! this is executed by the main process only
  if (myrank == 0) then

    ! appends receiver locations to sr.vtk file
    open(IOUT_VTK,file=trim(OUTPUT_FILES)//'/sr_tmp.vtk',position='append',status='old',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening and appending receivers to file sr_tmp.vtk')

    ! chooses best receivers locations
    islice_selected_rec_found(:) = -1
    nrec_found = 0
    do irec = 1,nrec

      if (final_distance(irec) == HUGEVAL) call exit_MPI(myrank,'Error locating receiver')

      if (DISPLAY_DETAILS_STATIONS .or. final_distance(irec) > 0.01d0) then
        write(IMAIN,*)
        write(IMAIN,*) 'Station #',irec,': ',trim(network_name(irec))//'.'//trim(station_name(irec))
        write(IMAIN,*) '       original latitude: ',sngl(stlat(irec))
        write(IMAIN,*) '      original longitude: ',sngl(stlon(irec))
        write(IMAIN,*) '     epicentral distance: ',sngl(epidist(irec))
        write(IMAIN,*) '  closest estimate found: ',sngl(final_distance(irec)),' km away'
        write(IMAIN,*) '   in slice ',islice_selected_rec(irec),' in element ',ispec_selected_rec(irec)
        write(IMAIN,*) '   at xi,eta,gamma coordinates = ',xi_receiver(irec),eta_receiver(irec),gamma_receiver(irec)
        ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
        call xyz_2_rlatlon_dble(x_found(irec),y_found(irec),z_found(irec),r,lat,lon)
        write(IMAIN,*) '   at lat/lon = ',sngl(lat),sngl(lon)
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

        islice_selected_rec_found(nrec_found) = islice_selected_rec(irec)
        ispec_selected_rec_found(nrec_found) = ispec_selected_rec(irec)
        xi_receiver_found(nrec_found) = xi_receiver(irec)
        eta_receiver_found(nrec_found) = eta_receiver(irec)
        gamma_receiver_found(nrec_found) = gamma_receiver(irec)
        station_name_found(nrec_found) = station_name(irec)
        network_name_found(nrec_found) = network_name(irec)
        stlat_found(nrec_found) = stlat(irec)
        stlon_found(nrec_found) = stlon(irec)
        stele_found(nrec_found) = stele(irec)
        stbur_found(nrec_found) = stbur(irec)
        nu_found(:,:,nrec_found) = nu(:,:,irec)
        epidist_found(nrec_found) = epidist(irec)

        ! writes out actual receiver location to VTK file
        write(IOUT_VTK,'(3e18.6)') sngl(x_found(irec)), sngl(y_found(irec)), sngl(z_found(irec))
      endif
    enddo

    ! finishes sr_tmp.vtk file
    write(IOUT_VTK,*)
    close(IOUT_VTK)

    ! compute maximal distance for all the receivers
    final_distance_max = maxval(final_distance(:))

    ! display maximum error for all the receivers
    write(IMAIN,*)
    write(IMAIN,*) 'maximum error in location of all the receivers: ',sngl(final_distance_max),' km'

    ! add warning if estimate is poor
    ! (usually means receiver outside the mesh given by the user)
    if (final_distance_max > THRESHOLD_EXCLUDE_STATION) then
      write(IMAIN,*)
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '***** WARNING: at least one receiver was excluded from the station list *****'
      write(IMAIN,*) '************************************************************'
      write(IMAIN,*) '************************************************************'
    endif
    call flush_IMAIN()

    nrec = nrec_found
    islice_selected_rec(1:nrec) = islice_selected_rec_found(1:nrec)
    ispec_selected_rec(1:nrec) = ispec_selected_rec_found(1:nrec)
    xi_receiver(1:nrec) = xi_receiver_found(1:nrec)
    eta_receiver(1:nrec) = eta_receiver_found(1:nrec)
    gamma_receiver(1:nrec) = gamma_receiver_found(1:nrec)
    station_name(1:nrec) = station_name_found(1:nrec)
    network_name(1:nrec) = network_name_found(1:nrec)
    stlat(1:nrec) = stlat_found(1:nrec)
    stlon(1:nrec) = stlon_found(1:nrec)
    stele(1:nrec) = stele_found(1:nrec)
    stbur(1:nrec) = stbur_found(1:nrec)
    nu(:,:,1:nrec) = nu_found(:,:,1:nrec)
    epidist(1:nrec) = epidist_found(1:nrec)

    ! write the list of stations and associated epicentral distance
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/output_list_stations.txt', &
          status='unknown',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file output_list_stations.txt')
    write(IOUT,*)
    write(IOUT,*) 'total number of stations: ',nrec
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
        write(IOUT,'(a8,1x,a3,6x,f8.4,1x,f9.4,1x,f6.1,1x,f6.1)') trim(station_name(irec)),trim(network_name(irec)), &
          sngl(stlat(irec)),sngl(stlon(irec)),sngl(stele(irec)),sngl(stbur(irec))
      enddo
      ! close receiver file
      close(IOUT)
    endif

  endif    ! end of section executed by main process only

  call synchronize_all()

  ! deallocate arrays
  deallocate(epidist)
  deallocate(ix_initial_guess,iy_initial_guess,iz_initial_guess)
  deallocate(x_target,y_target,z_target)
  deallocate(x_found,y_found,z_found)
  deallocate(final_distance)

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
  call bcast_all_dp(nu,nrec*3*3)

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
!--------------------
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
    call heap_sort_siftdown(i, N)
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

    call heap_sort_siftdown(1, i - 1)
  enddo

  contains

    subroutine heap_sort_siftdown(start, bottom)

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

    end subroutine heap_sort_siftdown

  end subroutine heap_sort_distances

