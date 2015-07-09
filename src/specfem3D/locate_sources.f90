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

  subroutine locate_sources(nspec,nglob,ibool, &
                            xstore,ystore,zstore, &
                            ELLIPTICITY,min_tshift_cmt_original)

  use constants_solver

  use shared_input_parameters, only: OUTPUT_FILES

  use specfem_par,only: &
    NSOURCES,myrank, &
    tshift_cmt,theta_source,phi_source, &
    DT,hdur,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
    rspl,espl,espl2,nspl,ibathy_topo, &
    LOCAL_TMP_PATH,SIMULATION_TYPE,TOPOGRAPHY, &
    xigll,yigll,zigll, &
    xi_source,eta_source,gamma_source,nu_source, &
    islice_selected_source,ispec_selected_source, &
    SAVE_SOURCE_MASK

  use specfem_par_movie,only: vtkdata_source_x,vtkdata_source_y,vtkdata_source_z

  implicit none

  integer,intent(in) :: nspec,nglob
  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nspec)

  ! arrays containing coordinates of the points
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: xstore,ystore,zstore

  logical,intent(in) :: ELLIPTICITY

  double precision,intent(out) :: min_tshift_cmt_original

  ! local parameters
  integer :: isource
  integer :: iprocloop
  integer :: i,j,k,ispec,iglob
  integer :: ier

  double precision :: ell
  double precision :: elevation
  double precision :: r0,dcost,p20
  double precision :: theta,phi
  double precision :: dist,typical_size
  double precision :: xi,eta,gamma,dx,dy,dz,dxi,deta

  ! topology of the control points of the surface element
  integer :: iax,iay,iaz
  integer :: iaddx(NGNOD),iaddy(NGNOD),iaddr(NGNOD)

  ! coordinates of the control points of the surface element
  double precision :: xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  integer :: iter_loop
  integer :: ia
  double precision :: x,y,z
  double precision :: xix,xiy,xiz
  double precision :: etax,etay,etaz
  double precision :: gammax,gammay,gammaz
  double precision :: dgamma

  double precision, dimension(NSOURCES) :: final_distance_source
  double precision, dimension(:), allocatable :: final_distance_source_subset

  double precision :: x_target_source,y_target_source,z_target_source
  double precision :: r_target_source

  integer :: isources_already_done,isource_in_this_subset
  integer, dimension(:), allocatable :: ispec_selected_source_subset

  integer, dimension(:,:), allocatable :: ispec_selected_source_all
  double precision, dimension(:,:), allocatable :: xi_source_all,eta_source_all,gamma_source_all, &
     final_distance_source_all,x_found_source_all,y_found_source_all,z_found_source_all

  double precision, dimension(:), allocatable :: xi_source_subset,eta_source_subset,gamma_source_subset

  double precision, dimension(NSOURCES) :: lat,long,depth
  double precision, dimension(6,NSOURCES) :: moment_tensor
  double precision :: radius

  double precision, dimension(:), allocatable :: x_found_source,y_found_source,z_found_source
  double precision :: r_found_source
  double precision :: st,ct,sp,cp
  double precision :: Mrr,Mtt,Mpp,Mrt,Mrp,Mtp
  double precision :: colat_source
  double precision :: distmin

  integer :: ix_initial_guess_source,iy_initial_guess_source,iz_initial_guess_source
  integer :: NSOURCES_SUBSET_current_size

  logical :: located_target

  integer :: iorientation
  double precision :: stazi,stdip,thetan,phin,n(3)
  integer :: imin,imax,jmin,jmax,kmin,kmax
  double precision :: f0,t0_ricker

  double precision, external :: get_cmt_scalar_moment
  double precision, external :: get_cmt_moment_magnitude

  ! mask source region (mask values are between 0 and 1, with 0 around sources)
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: mask_source

  ! event time
  integer :: yr,jda,ho,mi
  double precision :: sec

  ! timer MPI
  double precision :: time_start,tCPU
  double precision, external :: wtime

  ! get MPI starting time for all sources
  time_start = wtime()

  ! make sure we clean the future final array
  ispec_selected_source(:) = 0
  final_distance_source(:) = HUGEVAL

  ! read all the sources
  if (myrank == 0) then
    ! only master process reads in CMTSOLUTION file
    call get_cmt(yr,jda,ho,mi,sec,tshift_cmt,hdur,lat,long,depth,moment_tensor, &
                 DT,NSOURCES,min_tshift_cmt_original)
  endif

  ! broadcast the information read on the master to the nodes
  call bcast_all_dp(tshift_cmt,NSOURCES)
  call bcast_all_dp(hdur,NSOURCES)
  call bcast_all_dp(lat,NSOURCES)
  call bcast_all_dp(long,NSOURCES)
  call bcast_all_dp(depth,NSOURCES)

  call bcast_all_dp(moment_tensor,6*NSOURCES)
  call bcast_all_singledp(min_tshift_cmt_original)

  ! define topology of the control element
  call hex_nodes(iaddx,iaddy,iaddr)

  ! compute typical size of elements at the surface
  typical_size = TWO_PI * R_UNIT_SPHERE / (4.0 * NEX_XI_VAL)

  ! use 10 times the distance as a criterion for source detection
  typical_size = 10.0 * typical_size

  ! initializes source mask
  if (SAVE_SOURCE_MASK .and. SIMULATION_TYPE == 3) then
    allocate(mask_source(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating mask source array')
    mask_source(:,:,:,:) = 1.0_CUSTOM_REAL
  endif

  ! appends receiver locations to sr.vtk file
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
    allocate(final_distance_source_subset(NSOURCES_SUBSET_current_size), &
             ispec_selected_source_subset(NSOURCES_SUBSET_current_size), &
             xi_source_subset(NSOURCES_SUBSET_current_size), &
             eta_source_subset(NSOURCES_SUBSET_current_size), &
             gamma_source_subset(NSOURCES_SUBSET_current_size), &
             x_found_source(NSOURCES_SUBSET_current_size), &
             y_found_source(NSOURCES_SUBSET_current_size), &
             z_found_source(NSOURCES_SUBSET_current_size),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary source arrays')

    ! arrays to collect data
    if (myrank == 0) then
      allocate(ispec_selected_source_all(NSOURCES_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               xi_source_all(NSOURCES_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               eta_source_all(NSOURCES_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               gamma_source_all(NSOURCES_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               final_distance_source_all(NSOURCES_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               x_found_source_all(NSOURCES_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               y_found_source_all(NSOURCES_SUBSET_current_size,0:NPROCTOT_VAL-1), &
               z_found_source_all(NSOURCES_SUBSET_current_size,0:NPROCTOT_VAL-1),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary source arrays for gather')
    else
      ! dummy arrays
      allocate(ispec_selected_source_all(1,1), &
               xi_source_all(1,1), &
               eta_source_all(1,1), &
               gamma_source_all(1,1), &
               final_distance_source_all(1,1), &
               x_found_source_all(1,1), &
               y_found_source_all(1,1), &
               z_found_source_all(1,1),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary source dummy arrays for gather')
    endif
    ! use -1 as a flag to detect if gather fails for some reason
    ispec_selected_source_all(:,:) = -1

    ! make sure we clean the subset array before the gather
    ispec_selected_source_subset(:) = 0
    final_distance_source_subset(:) = HUGEVAL
    final_distance_source_all(:,:) = HUGEVAL

    ! loop over sources within this subset
    do isource_in_this_subset = 1,NSOURCES_SUBSET_current_size

      ! mapping from source number in current subset to real source number in all the subsets
      isource = isource_in_this_subset + isources_already_done

      ! convert geographic latitude lat (degrees) to geocentric colatitude theta (radians)
      call lat_2_geocentric_colat_dble(lat(isource),theta)

      phi = long(isource)*DEGREES_TO_RADIANS
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

        !   get the orientation of the seismometer
        thetan=(90.0d0+stdip)*DEGREES_TO_RADIANS
        phin=stazi*DEGREES_TO_RADIANS

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

      ! finds elevation of position
      if (TOPOGRAPHY) then
        call get_topo_bathy(lat(isource),long(isource),elevation,ibathy_topo)
        r0 = r0 + elevation/R_EARTH
      endif
      if (ELLIPTICITY) then
        dcost = dcos(theta)
! this is the Legendre polynomial of degree two, P2(cos(theta)), see the discussion above eq (14.4) in Dahlen and Tromp (1998)
        p20 = 0.5d0*(3.0d0*dcost*dcost-1.0d0)
        radius = r0 - depth(isource)*1000.0d0/R_EARTH
! get ellipticity using spline evaluation
        call spline_evaluation(rspl,espl,espl2,nspl,radius,ell)
! this is eq (14.4) in Dahlen and Tromp (1998)
        r0 = r0*(1.0d0-(2.0d0/3.0d0)*ell*p20)
      endif

      ! subtracts source depth (given in km)
      r_target_source = r0 - depth(isource)*1000.0d0/R_EARTH

      ! compute the Cartesian position of the source
      x_target_source = r_target_source*dsin(theta)*dcos(phi)
      y_target_source = r_target_source*dsin(theta)*dsin(phi)
      z_target_source = r_target_source*dcos(theta)

      ! set distance to huge initial value
      distmin = HUGEVAL

      ! flag to check that we located at least one target element
      located_target = .false.
      ix_initial_guess_source = 0
      iy_initial_guess_source = 0
      iz_initial_guess_source = 0

      do ispec = 1,nspec

        ! exclude elements that are too far from target
        if (USE_DISTANCE_CRITERION) then
          iglob = ibool(MIDX,MIDY,MIDZ,ispec)
          dist = dsqrt((x_target_source - dble(xstore(iglob)))**2 &
                     + (y_target_source - dble(ystore(iglob)))**2 &
                     + (z_target_source - dble(zstore(iglob)))**2)
          if (dist > typical_size) cycle
        endif

        ! define the interval in which we look for points
        if (USE_FORCE_POINT_SOURCE) then
          ! force sources will be put on an exact GLL point
          imin = 1
          imax = NGLLX

          jmin = 1
          jmax = NGLLY

          kmin = 1
          kmax = NGLLZ

        else
          ! double-couple CMTSOLUTION
          ! loop only on points inside the element
          ! exclude edges to ensure this point is not shared with other elements
          imin = 2
          imax = NGLLX - 1

          jmin = 2
          jmax = NGLLY - 1

          kmin = 2
          kmax = NGLLZ - 1
        endif

        do k = kmin,kmax
          do j = jmin,jmax
            do i = imin,imax

              ! keep this point if it is closer to the receiver
              iglob = ibool(i,j,k,ispec)
              dist = dsqrt((x_target_source - dble(xstore(iglob)))**2 &
                          +(y_target_source - dble(ystore(iglob)))**2 &
                          +(z_target_source - dble(zstore(iglob)))**2)
              if (dist < distmin) then
                distmin = dist
                ispec_selected_source_subset(isource_in_this_subset) = ispec
                ix_initial_guess_source = i
                iy_initial_guess_source = j
                iz_initial_guess_source = k
                located_target = .true.
              endif

            enddo
          enddo
        enddo

        ! calculates a Gaussian mask around source point
        if (SAVE_SOURCE_MASK .and. SIMULATION_TYPE == 3) then
          call calc_mask_source(mask_source,ispec,NSPEC,typical_size, &
                                x_target_source,y_target_source,z_target_source, &
                                ibool,xstore,ystore,zstore,NGLOB)
        endif

      ! end of loop on all the elements in current slice
      enddo

      ! *******************************************
      ! find the best (xi,eta,gamma) for the source
      ! *******************************************

      ! if we have not located a target element, the source is not in this slice
      ! therefore use first element only for fictitious iterative search
      if (.not. located_target) then
        ispec_selected_source_subset(isource_in_this_subset) = 1
        ix_initial_guess_source = MIDX
        iy_initial_guess_source = MIDY
        iz_initial_guess_source = MIDZ
      endif

      ! for point sources, the location will be exactly at a GLL point
      ! otherwise this tries to find best location
      if (USE_FORCE_POINT_SOURCE) then
        ! store xi,eta,gamma and x,y,z of point found
        ! note: they have range [1.0d0,NGLLX/Y/Z], used for point sources
        !          see e.g. in compute_add_sources.f90
        xi_source_subset(isource_in_this_subset) = dble(ix_initial_guess_source)
        eta_source_subset(isource_in_this_subset) = dble(iy_initial_guess_source)
        gamma_source_subset(isource_in_this_subset) = dble(iz_initial_guess_source)

        iglob = ibool(ix_initial_guess_source,iy_initial_guess_source, &
            iz_initial_guess_source,ispec_selected_source_subset(isource_in_this_subset))
        x_found_source(isource_in_this_subset) = xstore(iglob)
        y_found_source(isource_in_this_subset) = ystore(iglob)
        z_found_source(isource_in_this_subset) = zstore(iglob)

        ! compute final distance between asked and found (converted to km)
        final_distance_source_subset(isource_in_this_subset) = &
          dsqrt((x_target_source-x_found_source(isource_in_this_subset))**2 + &
                (y_target_source-y_found_source(isource_in_this_subset))**2 + &
                (z_target_source-z_found_source(isource_in_this_subset))**2)*R_EARTH/1000.d0

      else

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

          iglob = ibool(iax,iay,iaz,ispec_selected_source_subset(isource_in_this_subset))
          xelm(ia) = dble(xstore(iglob))
          yelm(ia) = dble(ystore(iglob))
          zelm(ia) = dble(zstore(iglob))

        enddo

        ! use initial guess in xi, eta and gamma
        xi = xigll(ix_initial_guess_source)
        eta = yigll(iy_initial_guess_source)
        gamma = zigll(iz_initial_guess_source)

        ! iterate to solve the non linear system
        do iter_loop = 1,NUM_ITER

          ! recompute Jacobian for the new point
          call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                                  xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

          ! compute distance to target location
          dx = - (x - x_target_source)
          dy = - (y - y_target_source)
          dz = - (z - z_target_source)

          ! compute increments
          dxi  = xix*dx + xiy*dy + xiz*dz
          deta = etax*dx + etay*dy + etaz*dz
          dgamma =  gammax*dx + gammay*dy + gammaz*dz

          ! impose limit on increments
          if (abs(dxi) > 0.3d0 ) dxi = sign(1.0d0,dxi)*0.3d0
          if (abs(deta) > 0.3d0 ) deta = sign(1.0d0,deta)*0.3d0
          if (abs(dgamma) > 0.3d0 ) dgamma = sign(1.0d0,dgamma)*0.3d0

          ! update values
          xi = xi + dxi
          eta = eta + deta
          gamma = gamma + dgamma

          ! impose that we stay in that element
          ! (useful if user gives a source outside the mesh for instance)
          ! we can go slightly outside the [1,1] segment since with finite elements
          ! the polynomial solution is defined everywhere
          ! can be useful for convergence of iterative scheme with distorted elements
          if (xi > 1.10d0) xi = 1.10d0
          if (xi < -1.10d0) xi = -1.10d0
          if (eta > 1.10d0) eta = 1.10d0
          if (eta < -1.10d0) eta = -1.10d0
          if (gamma > 1.10d0) gamma = 1.10d0
          if (gamma < -1.10d0) gamma = -1.10d0

        enddo

        ! compute final coordinates of point found
        call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                                xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

        ! store xi,eta,gamma and x,y,z of point found
        xi_source_subset(isource_in_this_subset) = xi
        eta_source_subset(isource_in_this_subset) = eta
        gamma_source_subset(isource_in_this_subset) = gamma
        x_found_source(isource_in_this_subset) = x
        y_found_source(isource_in_this_subset) = y
        z_found_source(isource_in_this_subset) = z

        ! compute final distance between asked and found (converted to km)
        final_distance_source_subset(isource_in_this_subset) = &
          dsqrt((x_target_source-x)**2 + &
                (y_target_source-y)**2 + &
                (z_target_source-z)**2)*R_EARTH/1000.d0

      endif ! USE_FORCE_POINT_SOURCE

    ! end of loop on all the sources
    enddo
    ! synchronizes processes
    call synchronize_all()

    ! now gather information from all the nodes
    call gather_all_i(ispec_selected_source_subset,NSOURCES_SUBSET_current_size, &
      ispec_selected_source_all,NSOURCES_SUBSET_current_size,NPROCTOT_VAL)

    ! checks that the gather operation went well
    if (myrank == 0) then
      if (minval(ispec_selected_source_all(:,:)) <= 0) then
        print*,'Error ispec all: procs = ',NPROCTOT_VAL,'sources subset size = ',NSOURCES_SUBSET_current_size
        print*,ispec_selected_source_all(:,:)
        call exit_MPI(myrank,'gather operation failed for source')
      endif
    endif

    call gather_all_dp(xi_source_subset,NSOURCES_SUBSET_current_size, &
      xi_source_all,NSOURCES_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(eta_source_subset,NSOURCES_SUBSET_current_size, &
      eta_source_all,NSOURCES_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(gamma_source_subset,NSOURCES_SUBSET_current_size, &
      gamma_source_all,NSOURCES_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(final_distance_source_subset,NSOURCES_SUBSET_current_size, &
      final_distance_source_all,NSOURCES_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(x_found_source,NSOURCES_SUBSET_current_size, &
      x_found_source_all,NSOURCES_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(y_found_source,NSOURCES_SUBSET_current_size, &
      y_found_source_all,NSOURCES_SUBSET_current_size,NPROCTOT_VAL)
    call gather_all_dp(z_found_source,NSOURCES_SUBSET_current_size, &
      z_found_source_all,NSOURCES_SUBSET_current_size,NPROCTOT_VAL)

    ! this is executed by main process only
    if (myrank == 0) then

      ! loop on all the sources within subsets
      do isource_in_this_subset = 1,NSOURCES_SUBSET_current_size

        ! mapping from source number in current subset to real source number in all the subsets
        isource = isources_already_done + isource_in_this_subset

        ! loop on all the results to determine the best slice
        distmin = HUGEVAL
        do iprocloop = 0,NPROCTOT_VAL-1
          if (final_distance_source_all(isource_in_this_subset,iprocloop) < distmin) then
            ! stores this slice's info
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
        ! different output for force point sources
        if (USE_FORCE_POINT_SOURCE) then
          write(IMAIN,*) '  i index of source in that element: ',nint(xi_source(isource))
          write(IMAIN,*) '  j index of source in that element: ',nint(eta_source(isource))
          write(IMAIN,*) '  k index of source in that element: ',nint(gamma_source(isource))
          write(IMAIN,*)
          write(IMAIN,*) '  component direction: ',COMPONENT_FORCE_SOURCE
          write(IMAIN,*)
          write(IMAIN,*) '  nu1 = ',nu_source(1,:,isource)
          write(IMAIN,*) '  nu2 = ',nu_source(2,:,isource)
          write(IMAIN,*) '  nu3 = ',nu_source(3,:,isource)
          write(IMAIN,*)
          write(IMAIN,*) '  at (x,y,z) coordinates = ',x_found_source(isource_in_this_subset),&
            y_found_source(isource_in_this_subset),z_found_source(isource_in_this_subset)

          ! prints frequency content for point forces
          f0 = hdur(isource)
          t0_ricker = 1.2d0/f0
          write(IMAIN,*) '  using a source of dominant frequency ',f0
          write(IMAIN,*) '  lambda_S at dominant frequency = ',3000./sqrt(3.)/f0
          write(IMAIN,*) '  lambda_S at highest significant frequency = ',3000./sqrt(3.)/(2.5*f0)
          write(IMAIN,*) '  t0_ricker = ',t0_ricker,'tshift_cmt = ',tshift_cmt(isource)
          write(IMAIN,*)
          write(IMAIN,*) '  half duration -> frequency: ',hdur(isource),' seconds**(-1)'
        else
          write(IMAIN,*) '   xi coordinate of source in that element: ',xi_source(isource)
          write(IMAIN,*) '  eta coordinate of source in that element: ',eta_source(isource)
          write(IMAIN,*) 'gamma coordinate of source in that element: ',gamma_source(isource)
          ! add message if source is a Heaviside
          if (hdur(isource) <= 5.*DT) then
            write(IMAIN,*)
            write(IMAIN,*) 'Source time function is a Heaviside, convolve later'
            write(IMAIN,*)
          endif
          write(IMAIN,*)
          write(IMAIN,*) ' half duration: ',hdur(isource),' seconds'
        endif
        write(IMAIN,*) '    time shift: ',tshift_cmt(isource),' seconds'
        write(IMAIN,*)
        write(IMAIN,*) 'magnitude of the source:'
        write(IMAIN,*) '     scalar moment M0 = ', &
          get_cmt_scalar_moment(Mxx(isource),Myy(isource),Mzz(isource),Mxy(isource),Mxz(isource),Myz(isource)),' dyne-cm'
        write(IMAIN,*) '  moment magnitude Mw = ', &
          get_cmt_moment_magnitude(Mxx(isource),Myy(isource),Mzz(isource),Mxy(isource),Mxz(isource),Myz(isource))
        write(IMAIN,*)

        ! writes out actual source position to VTK file
        write(IOUT_VTK,'(3e18.6)') sngl(x_found_source(isource_in_this_subset)), &
                                   sngl(y_found_source(isource_in_this_subset)), &
                                   sngl(z_found_source(isource_in_this_subset))

        ! get latitude, longitude and depth of the source that will be used
        call xyz_2_rthetaphi_dble(x_found_source(isource_in_this_subset), &
                                  y_found_source(isource_in_this_subset), &
                                  z_found_source(isource_in_this_subset), &
                                  r_found_source,theta_source(isource),phi_source(isource))
        call reduce(theta_source(isource),phi_source(isource))

        ! converts geocentric to geographic colatitude
        call geocentric_2_geographic_dble(theta_source(isource),colat_source)

        ! brings longitude between -PI and PI
        if (phi_source(isource)>PI) phi_source(isource)=phi_source(isource)-TWO_PI

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
        write(IMAIN,*) '      latitude: ',(PI_OVER_TWO-colat_source)*RADIANS_TO_DEGREES
        write(IMAIN,*) '     longitude: ',phi_source(isource)*RADIANS_TO_DEGREES
        write(IMAIN,*) '         depth: ',(r0-r_found_source)*R_EARTH/1000.0d0,' km'
        write(IMAIN,*)

        ! display error in location estimate
        write(IMAIN,*) 'Error in location of the source: ',sngl(final_distance_source(isource)),' km'

        ! add warning if estimate is poor
        ! (usually means source outside the mesh given by the user)
        if (final_distance_source(isource) > 50.d0) then
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
          vtkdata_source_x = sngl(x_found_source(isource_in_this_subset))
          vtkdata_source_y = sngl(y_found_source(isource_in_this_subset))
          vtkdata_source_z = sngl(z_found_source(isource_in_this_subset))
        endif

      enddo ! end of loop on all the sources within current source subset

    endif ! end of section executed by main process only

    ! deallocate arrays specific to each subset
    deallocate(final_distance_source_subset)
    deallocate(ispec_selected_source_subset)
    deallocate(xi_source_subset,eta_source_subset,gamma_source_subset)
    deallocate(x_found_source,y_found_source,z_found_source)
    deallocate(ispec_selected_source_all)
    deallocate(xi_source_all,eta_source_all,gamma_source_all,final_distance_source_all)
    deallocate(x_found_source_all,y_found_source_all,z_found_source_all)

  enddo ! end of loop over all source subsets

  ! display maximum error in location estimate
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'maximum error in location of the sources: ',sngl(maxval(final_distance_source)),' km'
    write(IMAIN,*)
    call flush_IMAIN()

    ! closing sr_tmp.vtk
    close(IOUT_VTK)
  endif
  call synchronize_all()

  ! main process broadcasts the results to all the slices
  call bcast_all_i(islice_selected_source,NSOURCES)
  call bcast_all_i(ispec_selected_source,NSOURCES)

  call bcast_all_dp(xi_source,NSOURCES)
  call bcast_all_dp(eta_source,NSOURCES)
  call bcast_all_dp(gamma_source,NSOURCES)

  ! stores source mask
  if (SAVE_SOURCE_MASK .and. SIMULATION_TYPE == 3) then
    call save_mask_source(myrank,mask_source,NSPEC,LOCAL_TMP_PATH)
    deallocate( mask_source )
  endif

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

  subroutine calc_mask_source(mask_source,ispec,NSPEC,typical_size, &
                            x_target_source,y_target_source,z_target_source, &
                            ibool,xstore,ystore,zstore,NGLOB)

! calculate a Gaussian function mask in the crust_mantle region
! which is 0 around the source locations and 1 everywhere else

  use constants

  implicit none

  integer :: ispec,NSPEC,NGLOB

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: mask_source
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: xstore,ystore,zstore
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: ibool

  double precision :: typical_size
  double precision :: x_target_source,y_target_source,z_target_source

  ! local parameters
  integer i,j,k,iglob
  double precision dist_sq,sigma_sq

  ! standard deviation for Gaussian
  ! (removes factor 10 added for search radius from typical_size)
  sigma_sq = typical_size * typical_size / 100.0

  ! loops over GLL points within this ispec element
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        ! gets distance (squared) to source
        iglob = ibool(i,j,k,ispec)
        dist_sq = (x_target_source - dble(xstore(iglob)))**2 &
                  +(y_target_source - dble(ystore(iglob)))**2 &
                  +(z_target_source - dble(zstore(iglob)))**2

        ! adds Gaussian function value to mask
        ! (mask value becomes 0 closer to source location, 1 everywhere else )
        mask_source(i,j,k,ispec) = mask_source(i,j,k,ispec) &
                  * ( 1.0_CUSTOM_REAL - exp( - dist_sq / sigma_sq ) )

      enddo
    enddo
  enddo

  end subroutine calc_mask_source

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_mask_source(myrank,mask_source,NSPEC,LOCAL_TMP_PATH)

! saves a mask in the crust_mantle region which is 0 around the source locations
! and 1 everywhere else

  use constants

  implicit none

  integer :: myrank,NSPEC

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: mask_source
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

!
!-------------------------------------------------------------------------------------------------
!

  subroutine print_stf(NSOURCES,isource,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                      tshift_cmt,hdur,min_tshift_cmt_original,NSTEP,DT)

! prints source time function

  use constants
  use shared_input_parameters

  implicit none

  integer :: NSOURCES,isource

  double precision,dimension(NSOURCES) :: Mxx,Myy,Mzz,Mxy,Mxz,Myz
  double precision,dimension(NSOURCES) :: tshift_cmt,hdur

  double precision :: min_tshift_cmt_original
  integer :: NSTEP
  double precision :: DT


  ! local parameters
  integer :: it,iom,ier
  double precision :: scalar_moment
  double precision :: t0, hdur_gaussian(NSOURCES)
  double precision :: t_cmt_used(NSOURCES)
  double precision time_source,om
  double precision :: f0

  double precision, external :: comp_source_time_function,comp_source_spectrum
  double precision, external :: comp_source_time_function_rickr
  double precision, external :: get_cmt_scalar_moment

  character(len=MAX_STRING_LEN) :: plot_file

  ! number of points to plot the source time function and spectrum
  integer, parameter :: NSAMP_PLOT_SOURCE = 1000

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) 'printing the source-time function'
  call flush_IMAIN()

  ! print the source-time function
  if (NSOURCES == 1) then
    plot_file = '/plot_source_time_function.txt'
  else
   if (isource < 10) then
      write(plot_file,"('/plot_source_time_function',i1,'.txt')") isource
    else if (isource < 100) then
      write(plot_file,"('/plot_source_time_function',i2,'.txt')") isource
    else
      write(plot_file,"('/plot_source_time_function',i3,'.txt')") isource
    endif
  endif

  ! output file
  open(unit=IOUT,file=trim(OUTPUT_FILES)//plot_file, &
        status='unknown',iostat=ier)
  if (ier /= 0 ) call exit_mpi(0,'Error opening plot_source_time_function file')

  ! calculates scalar moment M0
  scalar_moment = get_cmt_scalar_moment(Mxx(isource),Myy(isource),Mzz(isource),Mxy(isource),Mxz(isource),Myz(isource))

  ! define t0 as the earliest start time
  ! note: this calculation here is only used for outputting the plot_source_time_function file
  !          (see setup_sources_receivers.f90)
  t0 = - 1.5d0*minval( tshift_cmt(:) - hdur(:) )
  if (USE_FORCE_POINT_SOURCE ) t0 = - 1.2d0 * minval(tshift_cmt(:) - 1.0d0/hdur(:))

  t_cmt_used(:) = tshift_cmt(:)
  if (USER_T0 > 0.d0) then
    if (t0 <= USER_T0 + min_tshift_cmt_original) then
      t_cmt_used(:) = tshift_cmt(:) + min_tshift_cmt_original
      t0 = USER_T0
    endif
  endif
  ! convert the half duration for triangle STF to the one for Gaussian STF
  ! note: this calculation here is only used for outputting the plot_source_time_function file
  !          (see setup_sources_receivers.f90)
  hdur_gaussian(:) = hdur(:)/SOURCE_DECAY_MIMIC_TRIANGLE

  ! writes out source time function to file
  do it = 1,NSTEP
    time_source = dble(it-1)*DT-t0-t_cmt_used(isource)
    if (USE_FORCE_POINT_SOURCE) then
      ! Ricker source time function
      f0 = hdur(isource)
      write(IOUT,*) sngl(dble(it-1)*DT-t0), &
        sngl(FACTOR_FORCE_SOURCE*comp_source_time_function_rickr(time_source,f0))
    else
      ! Gaussian source time function
      write(IOUT,*) sngl(dble(it-1)*DT-t0), &
        sngl(scalar_moment*comp_source_time_function(time_source,hdur_gaussian(isource)))
    endif
  enddo
  close(IOUT)

  write(IMAIN,*)
  write(IMAIN,*) 'printing the source spectrum'
  call flush_IMAIN()

  ! print the spectrum of the derivative of the source from 0 to 1/8 Hz
  if (NSOURCES == 1) then
    plot_file = '/plot_source_spectrum.txt'
  else
   if (isource < 10) then
      write(plot_file,"('/plot_source_spectrum',i1,'.txt')") isource
    else if (isource < 100) then
      write(plot_file,"('/plot_source_spectrum',i2,'.txt')") isource
    else
      write(plot_file,"('/plot_source_spectrum',i3,'.txt')") isource
    endif
  endif

  open(unit=IOUT,file=trim(OUTPUT_FILES)//plot_file, &
        status='unknown',iostat=ier)
  if (ier /= 0 ) call exit_mpi(0,'Error opening plot_source_spectrum file')

  do iom = 1,NSAMP_PLOT_SOURCE
    om=TWO_PI*(1.0d0/8.0d0)*(iom-1)/dble(NSAMP_PLOT_SOURCE-1)
    write(IOUT,*) sngl(om/TWO_PI), &
      sngl(scalar_moment*om*comp_source_spectrum(om,hdur(isource)))
  enddo
  close(IOUT)

  end subroutine print_stf
