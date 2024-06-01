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
!----  locate_sources finds the correct position of the sources
!----

  subroutine locate_sources()

  use constants_solver

  use shared_parameters, only: OUTPUT_FILES,R_PLANET,RHOAV

  use specfem_par, only: NSOURCES,DT, &
    ibathy_topo,TOPOGRAPHY, &
    LOCAL_TMP_PATH,SIMULATION_TYPE, &
    SAVE_SOURCE_MASK

  use specfem_par, only: rspl_ellip,ellipicity_spline,ellipicity_spline2,nspl_ellip

  ! sources
  use specfem_par, only: &
    hdur,Mxx,Myy,Mzz,Mxy,Mxz,Myz,Mw,M0, &
    xi_source,eta_source,gamma_source,nu_source, &
    islice_selected_source,ispec_selected_source, &
    tshift_src,source_theta_ref,source_phi_ref,source_final_distance_max

  ! forces
  use specfem_par, only: &
    USE_FORCE_POINT_SOURCE, USE_MONOCHROMATIC_CMT_SOURCE,force_stf,factor_force_source, &
    comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP

  use specfem_par, only: &
    nspec => NSPEC_CRUST_MANTLE

  use specfem_par_movie, only: vtkdata_source_x,vtkdata_source_y,vtkdata_source_z

  implicit none

  ! local parameters
  integer :: iprocloop
  integer :: isource
  integer :: ier

  double precision, dimension(NSOURCES) :: srclat,srclon,srcdepth
  double precision, dimension(6,NSOURCES) :: moment_tensor

  double precision, dimension(NSOURCES) :: final_distance
  double precision, dimension(NSOURCES) :: rsurface_source

  ! point locations
  double precision, allocatable, dimension(:,:) :: xyz_target
  double precision, allocatable, dimension(:,:) :: xyz_found_subset
  double precision, allocatable, dimension(:,:,:) :: xyz_found_all

  integer :: NSOURCES_SUBSET_current_size
  integer :: isource_in_this_subset,isources_already_done

  integer, dimension(:), allocatable :: ispec_selected_subset
  integer, dimension(:,:), allocatable :: ispec_selected_all

  double precision, dimension(:), allocatable :: final_distance_subset
  double precision, dimension(:,:), allocatable :: final_distance_all

  double precision, dimension(:), allocatable :: xi_subset,eta_subset,gamma_subset
  double precision, dimension(:,:), allocatable :: xi_all,eta_all,gamma_all

  double precision :: lat,lon,depth,r_target

  double precision :: theta,phi
  double precision :: sint,cost,sinp,cosp

  double precision :: elevation
  double precision :: r0

  double precision :: r_found
  double precision :: Mrr,Mtt,Mpp,Mrt,Mrp,Mtp

  double precision :: distmin_not_squared
  double precision :: x_target,y_target,z_target
  double precision :: x,y,z
  double precision :: xi,eta,gamma

  integer :: ispec_selected

  integer :: iorientation
  double precision :: stazi,stdip
  double precision :: n(3),thetan,phin

  double precision :: f0,t0_ricker,scaleF,force_N
  double precision :: total_M0,total_Mw,total_force_N

  double precision, external :: get_cmt_scalar_moment
  double precision, external :: get_cmt_moment_magnitude
  double precision, external :: get_cmt_moment_magnitude_from_M0

  ! mask source region (mask values are between 0 and 1, with 0 around sources)
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: mask_source

  ! timer MPI
  double precision :: time_start,tCPU
  double precision, external :: wtime

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '********************'
    write(IMAIN,*) ' locating sources'
    write(IMAIN,*) '********************'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! get MPI starting time for all sources
  time_start = wtime()

  ! make sure we clean the future final array
  ispec_selected_source(:) = 0
  islice_selected_source(:) = -1

  final_distance(:) = HUGEVAL
  source_final_distance_max = HUGEVAL

  total_M0 = 0.d0
  total_Mw = 0.d0
  total_force_N = 0.d0

  source_phi_ref = 0.d0
  source_theta_ref = 0.d0

  ! normalized source radius
  r0 = R_UNIT_SPHERE

  ! read all the sources
  call read_source_locations(srclat,srclon,srcdepth,moment_tensor)

  ! initializes source mask
  if (SAVE_SOURCE_MASK .and. SIMULATION_TYPE == 3) then
    allocate(mask_source(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating mask source array')
    mask_source(:,:,:,:) = 1.0_CUSTOM_REAL
  endif

  ! appends locations to sr.vtk file
  if (myrank == 0) then
    open(IOUT_VTK,file=trim(OUTPUT_FILES)//'/sr_tmp.vtk', &
          position='append',status='old',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening and appending sources to file sr_tmp.vtk')
  endif

  ! loop on all the sources
  ! gather source information in subsets to reduce memory requirements

  ! loop over subsets of sources
  do isources_already_done = 0, NSOURCES, NSOURCES_SUBSET_MAX

    ! the size of the subset can be the maximum size, or less (if we are in the last subset,
    ! or if there are fewer sources than the maximum size of a subset)
    NSOURCES_SUBSET_current_size = min(NSOURCES_SUBSET_MAX, NSOURCES - isources_already_done)

    ! allocate arrays specific to each subset
    allocate(ispec_selected_subset(NSOURCES_SUBSET_current_size), &
             xi_subset(NSOURCES_SUBSET_current_size), &
             eta_subset(NSOURCES_SUBSET_current_size), &
             gamma_subset(NSOURCES_SUBSET_current_size), &
             xyz_found_subset(NDIM,NSOURCES_SUBSET_current_size), &
             final_distance_subset(NSOURCES_SUBSET_current_size), &
             xyz_target(NDIM,NSOURCES_SUBSET_current_size),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary source arrays')

    ! initializes
    ispec_selected_subset(:) = 0
    xi_subset(:) = 0.d0; eta_subset(:) = 0.d0; gamma_subset(:) = 0.d0
    xyz_found_subset(:,:) = 0.d0; xyz_target(:,:) = 0.d0
    final_distance_subset(:) = HUGEVAL

    ! arrays to collect data
    if (myrank == 0) then
      allocate(ispec_selected_all(NSOURCES_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               xi_all(NSOURCES_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               eta_all(NSOURCES_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               gamma_all(NSOURCES_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               xyz_found_all(NDIM,NSOURCES_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               final_distance_all(NSOURCES_SUBSET_current_size,0:NPROCTOT_VAL-1),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary source arrays for gather')
    else
      ! dummy arrays
      allocate(ispec_selected_all(1,1), &
               xi_all(1,1), &
               eta_all(1,1), &
               gamma_all(1,1), &
               xyz_found_all(1,1,1), &
               final_distance_all(1,1),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary source dummy arrays for gather')
    endif

    ! initializes
    ispec_selected_all(:,:) = 0
    xi_all(:,:) = 0.d0; eta_all(:,:) = 0.d0; gamma_all(:,:) = 0.d0
    xyz_found_all(:,:,:) = 0.d0; final_distance_all(:,:) = HUGEVAL

    ! loop over sources within this subset
    do isource_in_this_subset = 1,NSOURCES_SUBSET_current_size

      ! mapping from source number in current subset to real source number in all the subsets
      isource = isource_in_this_subset + isources_already_done

      ! source lat/lon in degrees
      lat = srclat(isource)
      lon = srclon(isource)

      ! limits longitude to [0.0,360.0]
      if (lon < 0.d0 ) lon = lon + 360.d0
      if (lon > 360.d0 ) lon = lon - 360.d0

      ! convert geographic latitude lat (degrees) to geocentric colatitude theta (radians)
      call lat_2_geocentric_colat_dble(lat,theta,ELLIPTICITY_VAL)

      ! longitude
      phi = lon*DEGREES_TO_RADIANS

      ! theta to [0,PI] and phi to [0,2PI]
      call reduce(theta,phi)

      ! stores first source theta/phi position (geocentric) as reference
      ! for epicentral distance calculations of receivers
      if (isource == 1) then
        source_phi_ref = phi
        source_theta_ref = theta
      endif

      sint = sin(theta)
      cost = cos(theta)
      sinp = sin(phi)
      cosp = cos(phi)

      ! get the moment tensor
      Mrr = moment_tensor(1,isource)
      Mtt = moment_tensor(2,isource)
      Mpp = moment_tensor(3,isource)
      Mrt = moment_tensor(4,isource)
      Mrp = moment_tensor(5,isource)
      Mtp = moment_tensor(6,isource)

      ! convert from a spherical to a Cartesian representation of the moment tensor
      Mxx(isource) = sint*sint*cosp*cosp*Mrr + cost*cost*cosp*cosp*Mtt + sinp*sinp*Mpp &
          + 2.0d0*sint*cost*cosp*cosp*Mrt - 2.0d0*sint*sinp*cosp*Mrp - 2.0d0*cost*sinp*cosp*Mtp

      Myy(isource) = sint*sint*sinp*sinp*Mrr + cost*cost*sinp*sinp*Mtt + cosp*cosp*Mpp &
          + 2.0d0*sint*cost*sinp*sinp*Mrt + 2.0d0*sint*sinp*cosp*Mrp + 2.0d0*cost*sinp*cosp*Mtp

      Mzz(isource) = cost*cost*Mrr + sint*sint*Mtt - 2.0d0*sint*cost*Mrt

      Mxy(isource) = sint*sint*sinp*cosp*Mrr + cost*cost*sinp*cosp*Mtt - sinp*cosp*Mpp &
          + 2.0d0*sint*cost*sinp*cosp*Mrt + sint*(cosp*cosp-sinp*sinp)*Mrp + cost*(cosp*cosp-sinp*sinp)*Mtp

      Mxz(isource) = sint*cost*cosp*Mrr - sint*cost*cosp*Mtt &
          + (cost*cost-sint*sint)*cosp*Mrt - cost*sinp*Mrp + sint*sinp*Mtp

      Myz(isource) = sint*cost*sinp*Mrr - sint*cost*sinp*Mtt &
          + (cost*cost-sint*sint)*sinp*Mrt + cost*cosp*Mrp - sint*cosp*Mtp

      ! record three components for each station
      do iorientation = 1,3
        !   North
        if (iorientation == 1) then
          stazi = 0.d0
          stdip = 0.d0
        !   East
        else if (iorientation == 2) then
          stazi = 90.d0
          stdip = 0.d0
        !   Vertical
        else if (iorientation == 3) then
          stazi = 0.d0
          stdip = - 90.d0
        else
          call exit_MPI(myrank,'incorrect orientation')
        endif

        ! get the orientation of the position
        thetan = (90.0d0+stdip)*DEGREES_TO_RADIANS
        phin = stazi*DEGREES_TO_RADIANS

        ! we use the same convention as in Harvard normal modes for the orientation

        !   vertical component
        n(1) = cos(thetan)
        !   N-S component
        n(2) = - sin(thetan)*cos(phin)
        !   E-W component
        n(3) = sin(thetan)*sin(phin)

        !   get the Cartesian components of n in the model: nu
        nu_source(iorientation,1,isource) = n(1)*sint*cosp + n(2)*cost*cosp - n(3)*sinp
        nu_source(iorientation,2,isource) = n(1)*sint*sinp + n(2)*cost*sinp + n(3)*cosp
        nu_source(iorientation,3,isource) = n(1)*cost - n(2)*sint
      enddo

      ! point depth (in m)
      depth = srcdepth(isource)*1000.0d0

      ! normalized source radius
      r0 = R_UNIT_SPHERE

      ! finds elevation of position
      if (TOPOGRAPHY) then
        call get_topo_bathy(lat,lon,elevation,ibathy_topo)
        r0 = r0 + elevation/R_PLANET
      endif

      ! ellipticity
      if (ELLIPTICITY_VAL) then
        ! adds ellipticity factor to radius
        call add_ellipticity_rtheta(r0,theta,nspl_ellip,rspl_ellip,ellipicity_spline,ellipicity_spline2)
      endif

      ! stores surface radius for info output
      rsurface_source(isource) = r0

      ! subtracts source depth (given in m)
      r0 = r0 - depth/R_PLANET

      ! source position
      r_target = r0

      ! compute the Cartesian position of the source
      x_target = r_target*sint*cosp
      y_target = r_target*sint*sinp
      z_target = r_target*cost

      ! stores Cartesian positions
      xyz_target(1,isource_in_this_subset) = x_target
      xyz_target(2,isource_in_this_subset) = y_target
      xyz_target(3,isource_in_this_subset) = z_target
    enddo

    ! make sure we clean the subset array before the gather
    ispec_selected_subset(:) = 0
    final_distance_subset(:) = HUGEVAL
    final_distance_all(:,:) = HUGEVAL

    ! find point locations
    do isource_in_this_subset = 1,NSOURCES_SUBSET_current_size

      ! mapping from source number in current subset to real source number in all the subsets
      isource = isource_in_this_subset + isources_already_done

      ! source lat/lon in degrees
      lat = srclat(isource)
      lon = srclon(isource)

      ! gets target position
      x_target = xyz_target(1,isource_in_this_subset)
      y_target = xyz_target(2,isource_in_this_subset)
      z_target = xyz_target(3,isource_in_this_subset)

      ! locates best element and xi/eta/gamma interpolation values
      call locate_point(x_target,y_target,z_target,lat,lon,ispec_selected,xi,eta,gamma, &
                        x,y,z,distmin_not_squared,.true.)

      ! store xi,eta,gamma and x,y,z of point found
      xi_subset(isource_in_this_subset) = xi
      eta_subset(isource_in_this_subset) = eta
      gamma_subset(isource_in_this_subset) = gamma

      xyz_found_subset(1,isource_in_this_subset) = x
      xyz_found_subset(2,isource_in_this_subset) = y
      xyz_found_subset(3,isource_in_this_subset) = z

      final_distance_subset(isource_in_this_subset) = distmin_not_squared
      ispec_selected_subset(isource_in_this_subset) = ispec_selected

      ! calculates a Gaussian mask around this source point
      if (SAVE_SOURCE_MASK .and. SIMULATION_TYPE == 3) then
        call calc_mask_source(mask_source,x_target,y_target,z_target)
      endif

    ! end of loop on all the sources
    enddo

    ! for MPI version, gather information from all the nodes
    ! use -1 as a flag to detect if gather fails for some reason
    ispec_selected_all(:,:) = -1

    ! now gather information from all the nodes
    call gather_all_i(ispec_selected_subset,NSOURCES_SUBSET_current_size, &
                      ispec_selected_all,NSOURCES_SUBSET_current_size,NPROCTOT_VAL)

    ! this is executed by main process only
    if (myrank == 0) then
      ! check that the gather operation went well
      if (any(ispec_selected_all(:,:) == -1)) then
        print *,'Error ispec all: procs = ',NPROCTOT_VAL,'sources subset size = ',NSOURCES_SUBSET_current_size
        print *,ispec_selected_all(:,:)
        call exit_MPI(myrank,'gather operation failed for source')
      endif
    endif

    call gather_all_dp(xi_subset,NSOURCES_SUBSET_current_size, &
                       xi_all,NSOURCES_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(eta_subset,NSOURCES_SUBSET_current_size, &
                       eta_all,NSOURCES_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(gamma_subset,NSOURCES_SUBSET_current_size, &
                       gamma_all,NSOURCES_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(final_distance_subset,NSOURCES_SUBSET_current_size, &
                       final_distance_all,NSOURCES_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(xyz_found_subset,NDIM*NSOURCES_SUBSET_current_size, &
                       xyz_found_all,NDIM*NSOURCES_SUBSET_current_size,NPROCTOT_VAL)

    ! this is executed by main process only
    if (myrank == 0) then

      ! selects best location in all slices
      ! loop on all the sources within subsets
      do isource_in_this_subset = 1,NSOURCES_SUBSET_current_size

        ! mapping from source number in current subset to real source number in all the subsets
        isource = isources_already_done + isource_in_this_subset

        ! loop on all the results to determine the best slice
        distmin_not_squared = HUGEVAL
        do iprocloop = 0,NPROCTOT_VAL-1
          if (final_distance_all(isource_in_this_subset,iprocloop) < distmin_not_squared) then
            ! stores this slice's info
            distmin_not_squared = final_distance_all(isource_in_this_subset,iprocloop)
            islice_selected_source(isource) = iprocloop
            ispec_selected_source(isource) = ispec_selected_all(isource_in_this_subset,iprocloop)

            xi_source(isource) = xi_all(isource_in_this_subset,iprocloop)
            eta_source(isource) = eta_all(isource_in_this_subset,iprocloop)
            gamma_source(isource) = gamma_all(isource_in_this_subset,iprocloop)

            xyz_found_subset(:,isource_in_this_subset) = xyz_found_all(:,isource_in_this_subset,iprocloop)
          endif
        enddo
        final_distance(isource) = distmin_not_squared
      enddo
    endif ! end of section executed by main process only

    ! output source infos
    if (myrank == 0) then

      ! loop on all the sources within subsets
      do isource_in_this_subset = 1,NSOURCES_SUBSET_current_size
        ! mapping from source number in current subset to real source number in all the subsets
        isource = isources_already_done + isource_in_this_subset

        ! source info
        write(IMAIN,*)
        write(IMAIN,*) 'source # ',isource
        write(IMAIN,*)
        write(IMAIN,*) '  source located in slice ',islice_selected_source(isource_in_this_subset)
        write(IMAIN,*) '                 in element ',ispec_selected_source(isource_in_this_subset)
        write(IMAIN,*)
        ! different output for force point sources
        if (USE_FORCE_POINT_SOURCE) then
          write(IMAIN,*) '  using force point source:'
          write(IMAIN,*) '    xi coordinate of source in that element: ',xi_source(isource)
          write(IMAIN,*) '    eta coordinate of source in that element: ',eta_source(isource)
          write(IMAIN,*) '    gamma coordinate of source in that element: ',gamma_source(isource)

          write(IMAIN,*)
          write(IMAIN,*) '    component of direction vector in East direction: ',comp_dir_vect_source_E(isource)
          write(IMAIN,*) '    component of direction vector in North direction: ',comp_dir_vect_source_N(isource)
          write(IMAIN,*) '    component of direction vector in Vertical direction: ',comp_dir_vect_source_Z_UP(isource)

          !write(IMAIN,*) '  i index of source in that element: ',nint(xi_source(isource))
          !write(IMAIN,*) '  j index of source in that element: ',nint(eta_source(isource))
          !write(IMAIN,*) '  k index of source in that element: ',nint(gamma_source(isource))
          !write(IMAIN,*)
          !write(IMAIN,*) '  component direction: ',COMPONENT_FORCE_SOURCE
          write(IMAIN,*)
          write(IMAIN,*) '    nu1 = ',nu_source(1,:,isource),'North'
          write(IMAIN,*) '    nu2 = ',nu_source(2,:,isource),'East'
          write(IMAIN,*) '    nu3 = ',nu_source(3,:,isource),'Vertical'
          write(IMAIN,*)
          write(IMAIN,*) '    at (x,y,z) coordinates = ',xyz_found_subset(1,isource_in_this_subset), &
            xyz_found_subset(2,isource_in_this_subset),xyz_found_subset(3,isource_in_this_subset)
        else
          ! moment tensor
          write(IMAIN,*) '  using moment tensor source:'
          write(IMAIN,*) '    xi coordinate of source in that element: ',xi_source(isource)
          write(IMAIN,*) '    eta coordinate of source in that element: ',eta_source(isource)
          write(IMAIN,*) '    gamma coordinate of source in that element: ',gamma_source(isource)
          write(IMAIN,*)
          write(IMAIN,*) '    at (x,y,z) coordinates = ',xyz_found_subset(1,isource_in_this_subset), &
            xyz_found_subset(2,isource_in_this_subset),xyz_found_subset(3,isource_in_this_subset)
        endif
        write(IMAIN,*)

        ! source time function info
        write(IMAIN,*) '  source time function:'
        if (EXTERNAL_SOURCE_TIME_FUNCTION) then
          ! external STF
          write(IMAIN,*) '    using external source time function'
          write(IMAIN,*)
        else
          ! frequency/half-duration
          if (USE_FORCE_POINT_SOURCE) then
            ! single point force
            ! prints frequency content for point forces
            select case(force_stf(isource))
            case (0)
              ! Gaussian
              write(IMAIN,*) '    using Gaussian source time function'
              write(IMAIN,*) '             half duration: ',hdur(isource),' seconds'
              write(IMAIN,*) '    Gaussian half duration: ',hdur(isource)/SOURCE_DECAY_MIMIC_TRIANGLE,' seconds'
            case (1)
              ! Ricker
              write(IMAIN,*) '    using Ricker source time function'
              ! prints frequency content for point forces
              f0 = hdur(isource)
              t0_ricker = 1.2d0/f0
              write(IMAIN,*)
              write(IMAIN,*) '    using a source of dominant frequency ',f0
              write(IMAIN,*) '    t0_ricker = ',t0_ricker,'tshift_src = ',tshift_src(isource)
              write(IMAIN,*)
              write(IMAIN,*) '    lambda_S at dominant frequency = ',3000./sqrt(3.)/f0
              write(IMAIN,*) '    lambda_S at highest significant frequency = ',3000./sqrt(3.)/(2.5*f0)
              write(IMAIN,*)
              write(IMAIN,*) '    half duration in frequency: ',hdur(isource),' seconds**(-1)'
            case (2)
              ! Heaviside
              write(IMAIN,*) '    using (quasi) Heaviside source time function'
              write(IMAIN,*) '             half duration: ',hdur(isource),' seconds'
            case (3)
              ! Monochromatic
              write(IMAIN,*) '    using monochromatic source time function'
              ! prints frequency content for point forces
              f0 = hdur(isource)
              write(IMAIN,*)
              write(IMAIN,*) '    using a source of period ',f0
              write(IMAIN,*)
              write(IMAIN,*) '    half duration in period: ',hdur(isource),' seconds'
            case (4)
              ! Gaussian by Meschede et al. (2011)
              write(IMAIN,*) '    using Gaussian source time function by Meschede et al. (2011), eq.(2)'
              write(IMAIN,*) '             tau: ',hdur(isource),' seconds'
            case default
              stop 'unsupported force_stf value!'
            end select
          else if (USE_MONOCHROMATIC_CMT_SOURCE) then
            ! moment tensor
            write(IMAIN,*) '    using monochromatic source time function'
            ! add message if source is monochromatic
            write(IMAIN,*)
            write(IMAIN,*) '    period: ',hdur(isource),' seconds'
          else
            ! moment tensor
            write(IMAIN,*) '    using (quasi) Heaviside source time function'
            ! add message if source is a Heaviside
            if (hdur(isource) <= 5.0*DT) then
              write(IMAIN,*)
              write(IMAIN,*) '    Source time function is a Heaviside, convolve later'
              write(IMAIN,*)
            endif
            write(IMAIN,*)
            write(IMAIN,*) '    half duration: ',hdur(isource),' seconds'
          endif
        endif
        write(IMAIN,*) '    time shift: ',tshift_src(isource),' seconds'
        write(IMAIN,*)

        ! magnitude
        write(IMAIN,*) '  magnitude of the source:'
        if (USE_FORCE_POINT_SOURCE) then
          ! single point force
          ! scale and non-dimensionalize the factor_force_source
          ! factor_force_source in FORCESOLUTION file is in Newton
          ! 1 Newton is 1 kg * 1 m / (1 second)^2
          scaleF = RHOAV * (R_PLANET**4) * PI*GRAV*RHOAV
          ! force in Newton
          force_N = factor_force_source(isource) * scaleF
          ! adds to total force applied (sum over all force point sources)
          total_force_N = total_force_N + force_N

          write(IMAIN,*) '    force = ', sngl(force_N),'(Newton)' ! dimensionalized
        else
          ! moment-tensor
          M0 = get_cmt_scalar_moment(Mxx(isource),Myy(isource),Mzz(isource),Mxy(isource),Mxz(isource),Myz(isource))
          Mw =  get_cmt_moment_magnitude(Mxx(isource),Myy(isource),Mzz(isource),Mxy(isource),Mxz(isource),Myz(isource))
          ! adds to total moment
          total_M0 = total_M0 + M0
          total_Mw = get_cmt_moment_magnitude_from_M0(total_M0)

          write(IMAIN,*) '       scalar moment M0 = ', M0,' dyne-cm'
          write(IMAIN,*) '    moment magnitude Mw = ', Mw
        endif
        write(IMAIN,*)

        ! writes out actual source position to VTK file
        write(IOUT_VTK,'(3e18.6)') sngl(xyz_found_subset(1,isource_in_this_subset)), &
                                   sngl(xyz_found_subset(2,isource_in_this_subset)), &
                                   sngl(xyz_found_subset(3,isource_in_this_subset))

        ! get latitude, longitude and depth of the source that will be used
        call xyz_2_rlatlon_dble(xyz_found_subset(1,isource_in_this_subset), &
                                xyz_found_subset(2,isource_in_this_subset), &
                                xyz_found_subset(3,isource_in_this_subset), &
                                r_found,lat,lon,ELLIPTICITY_VAL)

        ! brings longitude range ([0,360] by default) to source longitude range for display
        if (srclon(isource) < 0.d0) lon = lon - 360.d0

        write(IMAIN,*)
        write(IMAIN,*) '  original (requested) position of the source:'
        write(IMAIN,*)
        write(IMAIN,*) '        latitude: ',srclat(isource)
        write(IMAIN,*) '       longitude: ',srclon(isource)
        write(IMAIN,*) '           depth: ',srcdepth(isource),' km'
        write(IMAIN,*)

        ! compute real position of the source
        write(IMAIN,*) '  position of the source that will be used:'
        write(IMAIN,*)
        write(IMAIN,*) '        latitude: ',lat
        write(IMAIN,*) '       longitude: ',lon
        write(IMAIN,*) '           depth: ',(rsurface_source(isource)-r_found)*R_PLANET/1000.0d0,' km'
        write(IMAIN,*)

        ! display error in location estimate
        write(IMAIN,*) '  Error in location of the source: ',sngl(final_distance(isource)),' km'

        ! add warning if estimate is poor
        ! (usually means source outside the mesh given by the user)
        if (final_distance(isource) > 5.d0) then
          write(IMAIN,*)
          write(IMAIN,*) '*****************************************************'
          write(IMAIN,*) '*****************************************************'
          write(IMAIN,*) '***** WARNING: source location estimate is poor *****'
          write(IMAIN,*) '*****************************************************'
          write(IMAIN,*) '*****************************************************'
        endif
        call flush_IMAIN()

        ! stores location for VTK visualization
        if (isource == 1) then
          vtkdata_source_x = sngl(xyz_found_subset(1,isource_in_this_subset))
          vtkdata_source_y = sngl(xyz_found_subset(2,isource_in_this_subset))
          vtkdata_source_z = sngl(xyz_found_subset(3,isource_in_this_subset))
        endif

      enddo ! end of loop on all the sources within current source subset

    endif ! end of section executed by main process only

    ! deallocate arrays specific to each subset
    deallocate(ispec_selected_subset)
    deallocate(ispec_selected_all)
    deallocate(xi_subset,eta_subset,gamma_subset)
    deallocate(xi_all,eta_all,gamma_all)
    deallocate(final_distance_all)
    deallocate(final_distance_subset)
    deallocate(xyz_target)
    deallocate(xyz_found_subset)
    deallocate(xyz_found_all)

  enddo ! end of loop over all source subsets

  ! display maximum error in location estimate
  if (myrank == 0) then
    ! sets total magnitude (for finite sources)
    M0 = total_M0
    Mw = total_Mw
    force_N = total_force_N

    if (NSOURCES > 1) then
      write(IMAIN,*)
      write(IMAIN,*) '********************'
      write(IMAIN,*) 'finite source combined over all ',NSOURCES,' sources applied:'
      if (USE_FORCE_POINT_SOURCE) then
        ! total force in Newton
        write(IMAIN,*) '  total force = ', sngl(force_N),'(Newton)' ! dimensionalized
      else
        ! moment-tensor
        write(IMAIN,*) '     total scalar moment M0 = ', M0,' dyne-cm'
        write(IMAIN,*) '  total moment magnitude Mw = ', Mw
      endif
      write(IMAIN,*) '********************'
    endif

    ! compute maximal distance for all the sources
    source_final_distance_max = maxval(final_distance(:))

    write(IMAIN,*)
    write(IMAIN,*) 'maximum error in location of the sources: ',sngl(source_final_distance_max),' km'
    write(IMAIN,*)
    call flush_IMAIN()

    ! closing sr_tmp.vtk
    close(IOUT_VTK)
  endif

  ! main process broadcasts the results to all the slices
  call bcast_all_i(islice_selected_source,NSOURCES)
  call bcast_all_i(ispec_selected_source,NSOURCES)

  call bcast_all_dp(xi_source,NSOURCES)
  call bcast_all_dp(eta_source,NSOURCES)
  call bcast_all_dp(gamma_source,NSOURCES)

  ! Broadcast magnitude and scalar moment to all processers
  call bcast_all_singledp(M0)
  call bcast_all_singledp(Mw)

  ! stores source mask
  if (SAVE_SOURCE_MASK .and. SIMULATION_TYPE == 3) then
    call save_mask_source(mask_source,nspec,LOCAL_TMP_PATH)
    deallocate(mask_source)
  endif

  ! synchronizes to get right timing
  call synchronize_all()

  ! elapsed time since beginning of source detection
  if (myrank == 0) then
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for detection of sources in seconds = ',tCPU
    write(IMAIN,*)
    write(IMAIN,*) 'End of source detection - done'
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  call synchronize_all()

  end subroutine locate_sources


!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_source_locations(srclat,srclon,srcdepth,moment_tensor)

  use constants, only: &
    myrank

  use specfem_par, only: &
    DT,NSOURCES,tshift_src,hdur,min_tshift_src_original

  ! forces
  use specfem_par, only: &
    USE_FORCE_POINT_SOURCE,force_stf,factor_force_source, &
    comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP

  implicit none

  ! (uses these as explicit routine arguments to avoid compiler warnings)
  double precision, dimension(NSOURCES),intent(out) :: srclat,srclon,srcdepth
  double precision, dimension(6,NSOURCES),intent(out) :: moment_tensor

  ! local parameters
  ! event time
  integer :: yr,jda,mo,da,ho,mi
  double precision :: sec

  ! initializes
  srclat(:) = 0.d0
  srclon(:) = 0.d0
  srcdepth(:) = 0.d0
  moment_tensor(:,:) = 0.d0

  tshift_src(:) = 0.d0
  hdur(:) = 0.d0
  min_tshift_src_original = 0.d0

  ! reads in source descriptions
  if (USE_FORCE_POINT_SOURCE) then
    ! point forces
    if (myrank == 0) then
      ! only main process reads in FORCESOLUTION file
      call get_force(tshift_src,hdur, &
                     srclat,srclon,srcdepth,DT,NSOURCES, &
                     min_tshift_src_original,force_stf,factor_force_source, &
                     comp_dir_vect_source_E,comp_dir_vect_source_N,comp_dir_vect_source_Z_UP)
    endif
    ! broadcasts specific point force infos
    call bcast_all_i(force_stf,NSOURCES)
    call bcast_all_dp(factor_force_source,NSOURCES)
    call bcast_all_dp(comp_dir_vect_source_E,NSOURCES)
    call bcast_all_dp(comp_dir_vect_source_N,NSOURCES)
    call bcast_all_dp(comp_dir_vect_source_Z_UP,NSOURCES)
  else
    ! CMT moment tensors
    if (myrank == 0) then
      ! only main process reads in CMTSOLUTION file
      call get_cmt(yr,jda,mo,da,ho,mi,sec, &
                   tshift_src,hdur, &
                   srclat,srclon,srcdepth,moment_tensor, &
                   DT,NSOURCES,min_tshift_src_original)
    endif
    ! broadcast ispecific moment tensor infos
    call bcast_all_dp(moment_tensor,6*NSOURCES)
  endif

  ! broadcast the information read on the main node to all the nodes
  call bcast_all_dp(tshift_src,NSOURCES)
  call bcast_all_dp(hdur,NSOURCES)
  call bcast_all_dp(srclat,NSOURCES)
  call bcast_all_dp(srclon,NSOURCES)
  call bcast_all_dp(srcdepth,NSOURCES)
  call bcast_all_singledp(min_tshift_src_original)

  end subroutine read_source_locations

!
!-------------------------------------------------------------------------------------------------
!

  subroutine calc_mask_source(mask_source,x_target,y_target,z_target)

! calculate a Gaussian function mask in the crust_mantle region
! which is 0 around the source locations and 1 everywhere else

  use constants

  use specfem_par, only: &
    nspec => NSPEC_CRUST_MANTLE

  use specfem_par_crustmantle, only: &
    ibool => ibool_crust_mantle, &
    xstore => xstore_crust_mantle,ystore => ystore_crust_mantle,zstore => zstore_crust_mantle

  ! for point search
  use specfem_par, only: &
    typical_size_squared

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout) :: mask_source

  double precision,intent(in) :: x_target,y_target,z_target

  ! local parameters
  integer :: i,j,k,iglob,ispec
  double precision :: dist_squared,sigma_squared

  ! standard deviation for Gaussian
  ! (removes factor of 100 added for search radius from typical_size_squared)
  sigma_squared = typical_size_squared / 100.d0

  ! searches through all elements
  do ispec = 1,nspec
    ! loops over GLL points within this ispec element
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ! gets distance (squared) to source
          iglob = ibool(i,j,k,ispec)
          dist_squared = (x_target - dble(xstore(iglob)))**2 &
                       + (y_target - dble(ystore(iglob)))**2 &
                       + (z_target - dble(zstore(iglob)))**2

          ! adds Gaussian function value to mask
          ! (mask value becomes 0 closer to source location, 1 everywhere else )
          mask_source(i,j,k,ispec) = mask_source(i,j,k,ispec) &
                    * ( 1.0_CUSTOM_REAL - real(exp( - dist_squared / sigma_squared ),kind=CUSTOM_REAL) )
        enddo
      enddo
    enddo
  enddo

  end subroutine calc_mask_source

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_mask_source(mask_source,nspec,LOCAL_TMP_PATH)

! saves a mask in the crust_mantle region which is 0 around the source locations
! and 1 everywhere else

  use constants

  implicit none

  integer,intent(in) :: nspec

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: mask_source
  character(len=MAX_STRING_LEN) :: LOCAL_TMP_PATH

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: prname

  ! stores into file
  call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_TMP_PATH)

  open(unit=IOUT,file=trim(prname)//'mask_source.bin', &
       status='unknown',form='unformatted',action='write',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening mask_source.bin file')

  write(IOUT) mask_source
  close(IOUT)

  end subroutine save_mask_source

