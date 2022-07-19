!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  8 . 0
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
!---- locate_green_locations finds the locations to be s
!----

  subroutine locate_green_locations()

  use constants_solver, only: &
    ELLIPTICITY_VAL,NCHUNKS_VAL,NPROCTOT_VAL,NDIM, &
    MAX_LENGTH_STATION_NAME,MAX_LENGTH_NETWORK_NAME, &
    DISPLAY_DETAILS_GREEN_FUNCTIONS,NGF_SUBSET_MAX, &
    THRESHOLD_EXCLUDE_STATION, &
    HUGEVAL,IMAIN,IOUT,IOUT_VTK, &
    DEGREES_TO_RADIANS,RADIANS_TO_DEGREES,R_UNIT_SPHERE

  use shared_parameters, only: OUTPUT_FILES,R_PLANET

  use specfem_par, only: &
    myrank,DT,NSTEP, &
    ngf,islice_selected_gf_loc,ispec_selected_gf_loc, &
    xi_gf_loc,eta_gf_loc,gamma_gf_loc, &
    gf_loc_lat,gf_loc_lon,gf_loc_depth,gf_final_distance_max, nu_gf_loc, &
    rspl,ellipicity_spline,ellipicity_spline2,nspl,ibathy_topo, &
    TOPOGRAPHY,RECEIVERS_CAN_BE_BURIED

  implicit none

  ! local parameters
  integer :: iprocloop
  integer :: igf,i
  integer :: ier

  integer, dimension(ngf) :: islice_selected_found,ispec_selected_found
  double precision, dimension(ngf) :: xi_gf_loc_found,eta_gf_loc_found,gamma_gf_loc_found
  double precision, dimension(3,3,ngf) :: nu_gf_loc_found
  integer :: ngf_found

  ! point locations
  double precision, allocatable, dimension(:,:) :: xyz_target
  double precision, allocatable, dimension(:,:) :: xyz_found

  double precision, allocatable, dimension(:,:) :: xyz_found_subset
  double precision, allocatable, dimension(:,:,:) :: xyz_found_all

  double precision, dimension(ngf) :: gf_loc_lat_found,gf_loc_lon_found,gf_loc_depth_found

  integer :: ngf_SUBSET_current_size
  integer :: igf_in_this_subset,igf_already_done

  integer, allocatable, dimension(:) :: ispec_selected_subset
  integer, allocatable, dimension(:,:) :: ispec_selected_all

  double precision, dimension(:), allocatable :: final_distance
  double precision, allocatable, dimension(:) :: final_distance_subset
  double precision, dimension(:,:), allocatable :: final_distance_all

  double precision, allocatable, dimension(:) :: xi_subset,eta_subset,gamma_subset
  double precision, allocatable, dimension(:,:) :: xi_all,eta_all,gamma_all

  double precision :: lat,lon,radius,depth,r_target

  double precision :: theta,phi
  double precision :: sint,cost,sinp,cosp

  double precision :: ell
  double precision :: elevation
  double precision :: r0,p20

  double precision :: distmin_not_squared
  double precision :: x_target,y_target,z_target
  double precision :: x,y,z
  double precision :: xi,eta,gamma

  integer :: ispec_selected

  ! For computing rotation matrices for the orientation of the particle motion.
  integer :: iorientation
  double precision :: stazi,stdip
  double precision :: n(3),thetan,phin

  character(len=2) :: bic

  ! timer MPI
  double precision :: time_start,tCPU
  double precision, external :: wtime

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '**************************'
    write(IMAIN,*) ' Green Function Locations '
    write(IMAIN,*) '**************************'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! get MPI starting time
  time_start = wtime()

  ! allocate memory for arrays using number of stations
  allocate(xyz_target(NDIM,ngf), &
           xyz_found(NDIM,ngf), &
           final_distance(ngf),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating temporary green function arrays')

  ! read that STATIONS file on the main
  call read_green_locations()

  ! loop on all the stations to locate them in the mesh
  do igf = 1,ngf

    ! station lat/lon in degrees
    lat = gf_loc_lat(igf)
    lon = gf_loc_lon(igf)

    ! limits longitude to [0.0,360.0]
    if (lon < 0.d0 ) lon = lon + 360.d0
    if (lon > 360.d0 ) lon = lon - 360.d0

    ! converts geographic latitude stlat (degrees) to geocentric colatitude theta (radians)
    call lat_2_geocentric_colat_dble(lat,theta)

    phi = lon*DEGREES_TO_RADIANS
    call reduce(theta,phi)

    sint = sin(theta)
    cost = cos(theta)
    sinp = sin(phi)
    cosp = cos(phi)

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
      nu_gf_loc(iorientation,1,igf) = n(1)*sint*cosp + n(2)*cost*cosp - n(3)*sinp
      nu_gf_loc(iorientation,2,igf) = n(1)*sint*sinp + n(2)*cost*sinp + n(3)*cosp
      nu_gf_loc(iorientation,3,igf) = n(1)*cost - n(2)*sint
    enddo

    ! point depth (in m)
    depth = gf_loc_depth(igf)

    ! normalized receiver radius
    r0 = R_UNIT_SPHERE

    ! finds elevation of receiver
    if (TOPOGRAPHY) then
       call get_topo_bathy(lat,lon,elevation,ibathy_topo)
       r0 = r0 + elevation/R_PLANET
    endif

    ! ellipticity
    if (ELLIPTICITY_VAL) then
      ! this is the Legendre polynomial of degree two, P2(cos(theta)),
      ! see the discussion above eq (14.4) in Dahlen and Tromp (1998)
      p20 = 0.5d0*(3.0d0*cost*cost-1.0d0)

      ! get ellipticity using spline evaluation
      call spline_evaluation(rspl,ellipicity_spline,ellipicity_spline2,nspl,r0,ell)

      ! this is eq (14.4) in Dahlen and Tromp (1998)
      r0 = r0*(1.0d0-(2.0d0/3.0d0)*ell*p20)
    endif

    ! subtract location depth (in meters)
    r_target = r0 - depth/R_PLANET

    ! compute the Cartesian position of the receiver
    x_target = r_target*sint*cosp
    y_target = r_target*sint*sinp
    z_target = r_target*cost

    ! stores Cartesian positions
    xyz_target(1,igf) = x_target
    xyz_target(2,igf) = y_target
    xyz_target(3,igf) = z_target
  enddo

  ! make sure we clean the array before the gather
  ispec_selected_gf_loc(:) = 0
  islice_selected_gf_loc(:) = -1

  ! initializes search distances
  final_distance(:) = HUGEVAL
  gf_final_distance_max = HUGEVAL

  ! loop on all the receivers
  ! gather receiver information in subsets to reduce memory requirements

  ! loop over subsets of receivers
  do igf_already_done = 0, ngf, NGF_SUBSET_MAX

    ! the size of the subset can be the maximum size, or less (if we are in the last subset,
    ! or if there are fewer sources than the maximum size of a subset)
    ngf_SUBSET_current_size = min(NGF_SUBSET_MAX, ngf - igf_already_done)

    ! allocate arrays specific to each subset
    allocate(ispec_selected_subset(ngf_SUBSET_current_size), &
             xi_subset(ngf_SUBSET_current_size), &
             eta_subset(ngf_SUBSET_current_size), &
             gamma_subset(ngf_SUBSET_current_size), &
             xyz_found_subset(NDIM,ngf_SUBSET_current_size), &
             final_distance_subset(ngf_SUBSET_current_size),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary green function arrays')

    ! gather arrays
    if (myrank == 0) then
      ! only main process needs full arrays allocated
      allocate(ispec_selected_all(ngf_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               xi_all(ngf_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               eta_all(ngf_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               gamma_all(ngf_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               xyz_found_all(NDIM,ngf_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               final_distance_all(ngf_SUBSET_current_size,0:NPROCTOT_VAL-1),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary gather green function arrays')
    else
      ! dummy arrays
      allocate(ispec_selected_all(1,1), &
               xi_all(1,1), &
               eta_all(1,1), &
               gamma_all(1,1), &
               xyz_found_all(1,1,1), &
               final_distance_all(1,1),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary dummy green function arrays')
    endif

    ! initializes search results
    ispec_selected_subset(:) = 0
    final_distance_subset(:) = HUGEVAL
    final_distance_all(:,:) = HUGEVAL

    ! find point locations
    ! loop over the stations within this subset
    do igf_in_this_subset = 1,ngf_SUBSET_current_size

      ! mapping from station number in current subset to real station number in all the subsets
      igf = igf_in_this_subset + igf_already_done

      ! station lat/lon in degrees
      lat = gf_loc_lat(igf)
      lon = gf_loc_lon(igf)

      ! gets target position
      x_target = xyz_target(1,igf)
      y_target = xyz_target(2,igf)
      z_target = xyz_target(3,igf)

      ! locates best element and xi/eta/gamma interpolation values
      call locate_point(x_target,y_target,z_target,lat,lon,ispec_selected,xi,eta,gamma, &
                        x,y,z,distmin_not_squared,RECEIVERS_CAN_BE_BURIED)

      ! store xi,eta and x,y,z of point found
      xi_subset(igf_in_this_subset) = xi
      eta_subset(igf_in_this_subset) = eta
      gamma_subset(igf_in_this_subset) = gamma

      xyz_found(1,igf) = x
      xyz_found(2,igf) = y
      xyz_found(3,igf) = z

      xyz_found_subset(1,igf_in_this_subset) = x
      xyz_found_subset(2,igf_in_this_subset) = y
      xyz_found_subset(3,igf_in_this_subset) = z

      final_distance(igf) = distmin_not_squared

      final_distance_subset(igf_in_this_subset) = distmin_not_squared
      ispec_selected_subset(igf_in_this_subset) = ispec_selected

    enddo ! end of loop on all stations within current subset

    ! for MPI version, gather information from all the nodes
    ispec_selected_all(:,:) = -1

    call gather_all_i(ispec_selected_subset,ngf_SUBSET_current_size, &
                      ispec_selected_all,ngf_SUBSET_current_size,NPROCTOT_VAL)

    ! this is executed by main process only
    if (myrank == 0) then
      ! check that the gather operation went well
      if (any(ispec_selected_all(:,:) == -1)) then
        print *,'Error ispec all: procs = ',NPROCTOT_VAL,'green function location subset size = ',ngf_SUBSET_current_size
        print *,ispec_selected_all(:,:)
        call exit_MPI(myrank,'gather operation failed for green function locations')
      endif
    endif

    call gather_all_dp(xi_subset,ngf_SUBSET_current_size, &
                       xi_all,ngf_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(eta_subset,ngf_SUBSET_current_size, &
                       eta_all,ngf_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(gamma_subset,ngf_SUBSET_current_size, &
                       gamma_all,ngf_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(final_distance_subset,ngf_SUBSET_current_size, &
                       final_distance_all,ngf_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(xyz_found_subset,NDIM*ngf_SUBSET_current_size, &
                       xyz_found_all,NDIM*ngf_SUBSET_current_size,NPROCTOT_VAL)

    ! this is executed by main process only
    if (myrank == 0) then

      ! selects best location in all slices
      ! MPI loop on all the results to determine the best slice
      do igf_in_this_subset = 1,ngf_SUBSET_current_size

        ! mapping from station number in current subset to real station number in all the subsets
        igf = igf_in_this_subset + igf_already_done

        ! loop on all the results to determine the best slice
        distmin_not_squared = HUGEVAL
        do iprocloop = 0,NPROCTOT_VAL-1
          if (final_distance_all(igf_in_this_subset,iprocloop) < distmin_not_squared) then
            ! stores this slice's info
            distmin_not_squared = final_distance_all(igf_in_this_subset,iprocloop)
            islice_selected_gf_loc(igf) = iprocloop
            ispec_selected_gf_loc(igf) = ispec_selected_all(igf_in_this_subset,iprocloop)

            xi_gf_loc(igf) = xi_all(igf_in_this_subset,iprocloop)
            eta_gf_loc(igf) = eta_all(igf_in_this_subset,iprocloop)
            gamma_gf_loc(igf) = gamma_all(igf_in_this_subset,iprocloop)

            xyz_found(:,igf) = xyz_found_all(:,igf_in_this_subset,iprocloop)
          endif
        enddo
        final_distance(igf) = distmin_not_squared
        if (final_distance(igf) == HUGEVAL) call exit_MPI(myrank,'Error locating green function location')

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
    open(IOUT_VTK,file=trim(OUTPUT_FILES)//'/gf_tmp.vtk',position='append',status='old',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening and appending green function locations to file gf_tmp.vtk')

    ! chooses best receivers locations
    ! if receiver location is too far off, we will exclude the receiver
    islice_selected_found(:) = -1
    ngf_found = 0

    do igf = 1,ngf

      if (DISPLAY_DETAILS_GREEN_FUNCTIONS .or. final_distance(igf) > 0.01d0 .or. ngf < 50) then

        ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
        call xyz_2_rlatlon_dble(xyz_found(1,igf),xyz_found(2,igf),xyz_found(3,igf),radius,lat,lon)

        ! Write table header 
        if (igf == 1) then
          write(IMAIN,*)
          write(IMAIN,*) 
          write(IMAIN,'(A6,A14,A14,A14,A14,A14,A14,A14,A14,A14,A14)') &
            'GF #', 'orig. lat','orig. lon','dist (km)', &
            'slice', 'element', &
            'xi','eta','gamma', &
            'new lat','new lon'
        endif
        
        ! Write table entry
        write(IMAIN,'(I6,F14.4,F14.4,F14.4,I14,I14,F14.4,F14.4,F14.4,F14.4,F14.4)') &
          igf, sngl(gf_loc_lat(igf)), sngl(gf_loc_lon(igf)), sngl(final_distance(igf)), &
          islice_selected_gf_loc(igf), ispec_selected_gf_loc(igf), &
          xi_gf_loc(igf),eta_gf_loc(igf),gamma_gf_loc(igf), &
          sngl(lat),sngl(lon)
      endif

      ngf_found = ngf_found + 1

      islice_selected_found(ngf_found) = islice_selected_gf_loc(igf)
      ispec_selected_found(ngf_found) = ispec_selected_gf_loc(igf)

      xi_gf_loc_found(ngf_found) = xi_gf_loc(igf)
      eta_gf_loc_found(ngf_found) = eta_gf_loc(igf)
      gamma_gf_loc_found(ngf_found) = gamma_gf_loc(igf)

      gf_loc_lat_found(ngf_found) = gf_loc_lat(igf)
      gf_loc_lon_found(ngf_found) = gf_loc_lon(igf)
      gf_loc_depth_found(ngf_found) = gf_loc_depth(igf)

      nu_gf_loc_found(:,:,ngf_found) = nu_gf_loc(:,:,igf)

      ! writes out actual receiver location to VTK file
      write(IOUT_VTK,'(3e18.6)') sngl(xyz_found(1,igf)), sngl(xyz_found(2,igf)), sngl(xyz_found(3,igf))
    enddo

    ! finishes sr_tmp.vtk file
    write(IOUT_VTK,*)
    close(IOUT_VTK)

    ! compute maximal distance for all the receivers
    gf_final_distance_max = maxval(final_distance(:))

    ! display maximum error for all the receivers
    write(IMAIN,*)
    write(IMAIN,*) 'maximum error in location of all the receivers: ',sngl(gf_final_distance_max),' km'

    ! replaces list of all receivers with list of only those which were found
    ! (in particular: for a 1-chunk simulation, only stations in this chunk)
    ngf = ngf_found

    islice_selected_gf_loc(1:ngf) = islice_selected_found(1:ngf)
    ispec_selected_gf_loc(1:ngf) = ispec_selected_found(1:ngf)

    xi_gf_loc(1:ngf) = xi_gf_loc_found(1:ngf)
    eta_gf_loc(1:ngf) = eta_gf_loc_found(1:ngf)
    gamma_gf_loc(1:ngf) = gamma_gf_loc_found(1:ngf)

    gf_loc_lat(1:ngf) = gf_loc_lat_found(1:ngf)
    gf_loc_lon(1:ngf) = gf_loc_lon_found(1:ngf)
    gf_loc_depth(1:ngf) = gf_loc_depth_found(1:ngf)  

  endif    ! end of section executed by main process only

  ! main process broadcasts the results to all the slices
  call bcast_all_singlei(ngf)

  call bcast_all_i(islice_selected_gf_loc,ngf)
  call bcast_all_i(ispec_selected_gf_loc,ngf)

  call bcast_all_dp(xi_gf_loc,ngf)
  call bcast_all_dp(eta_gf_loc,ngf)
  call bcast_all_dp(gamma_gf_loc,ngf)

  call bcast_all_dp(gf_loc_lat,ngf)
  call bcast_all_dp(gf_loc_lon,ngf)
  call bcast_all_dp(gf_loc_depth,ngf)
  call bcast_all_dp(nu_gf_loc,ngf*3*3)

  ! deallocate arrays
  deallocate(xyz_found)
  deallocate(final_distance)

  ! synchronizes to get right timing
  call synchronize_all()

  ! elapsed time since beginning of mesh generation
  if (myrank == 0) then
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for GF Location detection in seconds = ',tCPU
    write(IMAIN,*)
    write(IMAIN,*) 'End of receiver detection - done'
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  end subroutine locate_green_locations

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_green_locations()

  use constants_solver, only: &
    IMAIN,IIN

  use specfem_par, only: &
    myrank,ngf,GF_LOCATIONS_FILE,RECEIVERS_CAN_BE_BURIED, &
    gf_loc_lat,gf_loc_lon,gf_loc_depth

  implicit none

  ! local parameters
  integer :: ier,igf,i
  character(len=256) :: line

  ! reads stations files
  if (myrank == 0) then
    ! user output
    write(IMAIN,*) 'reading green function location information...'
    write(IMAIN,*)
    call flush_IMAIN()

    ! opens station file GF_LOCATIONS FILE
    open(unit=IIN,file=trim(GF_LOCATIONS_FILE),status='old',action='read',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening GF_LOCATIONS file')

    ! loop on all the stations to read station information
    do igf = 1,ngf

      ! old line:
      !read(IIN,*,iostat=ier) station_name(igf),network_name(igf),stlat(igf),stlon(igf),stele(igf),stbur(igf)

      ! reads in line as string
      read(IIN,"(a256)",iostat=ier) line
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading in station ',igf
        call exit_MPI(myrank,'Error reading in station in STATIONS file')
      endif

      ! skips empty lines and comment lines
      do while( len_trim(line) == 0 .or. line(1:1) == '#')
        read(IIN,"(a256)",iostat=ier) line
        if (ier /= 0) then
          write(IMAIN,*) 'Error reading in station ',igf
          call exit_MPI(myrank,'Error reading in station in STATIONS file')
        endif
      enddo

      ! reads in station information
      read(line(1:len_trim(line)),*,iostat=ier) gf_loc_lat(igf),gf_loc_lon(igf),gf_loc_depth(igf)
      if (ier /= 0) then
        write(IMAIN,*) 'Error reading in station ',igf
        call exit_MPI(myrank,'Error reading in station in STATIONS file')
      endif

      ! checks latitude
      if (gf_loc_lat(igf) < -90.d0 .or. gf_loc_lat(igf) > 90.d0) then
        write(IMAIN,*) 'Error in line ', igf,': latitude ', gf_loc_lat(igf), &
                       ' is invalid, please check STATIONS record'
        call exit_MPI(myrank,'Error station latitude invalid')
      endif

    enddo
    ! close receiver file
    close(IIN)

  endif

  ! broadcast the information read on the main node to all the nodes
  call bcast_all_dp(gf_loc_lat,ngf)
  call bcast_all_dp(gf_loc_lon,ngf)
  call bcast_all_dp(gf_loc_depth,ngf)

  end subroutine read_green_locations