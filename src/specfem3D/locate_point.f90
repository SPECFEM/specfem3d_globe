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


  subroutine locate_point(x_target,y_target,z_target,lat_target,lon_target,ispec_selected,xi,eta,gamma, &
                          x,y,z,distmin_not_squared,CAN_BE_BURIED)

! locates target point inside best mesh element
  use constants_solver, only: &
    NGLLX,NGLLY,NGLLZ,MIDX,MIDY,MIDZ,NGNOD,R_EARTH_KM,HUGEVAL, &
    NUM_ITER,USE_DISTANCE_CRITERION

  use specfem_par, only: &
    nspec => NSPEC_CRUST_MANTLE

  use specfem_par_crustmantle, only: &
    ibool => ibool_crust_mantle, &
    xstore => xstore_crust_mantle,ystore => ystore_crust_mantle,zstore => zstore_crust_mantle

  ! for point search
  use specfem_par,only: &
    typical_size_squared, &
    anchor_iax,anchor_iay,anchor_iaz, &
    lat_min,lat_max,lon_min,lon_max,xyz_midpoints, &
    xigll,yigll,zigll

  implicit none

  double precision,intent(in) :: x_target,y_target,z_target
  double precision,intent(in) :: lat_target,lon_target

  integer,intent(out) :: ispec_selected
  double precision,intent(out) :: xi,eta,gamma
  double precision,intent(out) :: x,y,z
  double precision,intent(out) :: distmin_not_squared

  logical,intent(in) :: CAN_BE_BURIED

  ! local parameters
  integer :: ix_initial_guess,iy_initial_guess,iz_initial_guess
  integer :: ispec,i,j,k,iglob
  integer :: ia,iter_loop

  ! coordinates of the control points of the surface element
  double precision :: xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)
  double precision :: lat,lon
  double precision :: distmin_squared,dist_squared
  double precision :: dx,dy,dz,dx_min,dy_min,dz_min,d_min_sq
  double precision :: dxi,deta,dgamma
  double precision :: xix,xiy,xiz
  double precision :: etax,etay,etaz
  double precision :: gammax,gammay,gammaz
  logical :: target_located

  ! set distance to huge initial value
  distmin_squared = HUGEVAL

  ! initializes located target
  ! if we have not located a target element, the receiver is not in this slice
  ! therefore use first element only for fictitious iterative search
  ispec_selected = 1
  ix_initial_guess = MIDX
  iy_initial_guess = MIDY
  iz_initial_guess = MIDZ

  ! limits latitude to [-90.0,90.0]
  lat = lat_target
  if (lat < -90.d0) lat = -90.d0
  if (lat > 90.d0) lat = 90.d0

  ! limits longitude to [0.0,360.0]
  lon = lon_target
  if (lon < 0.d0) lon = lon + 360.d0
  if (lon > 360.d0) lon = lon - 360.d0

  ! checks if receiver in this slice
  if (lat >= lat_min .and. lat <= lat_max .and. &
      lon >= lon_min .and. lon <= lon_max) then
    target_located = .true.
  else
    target_located = .false.
  endif

  ! debug
  !print*,'target located:',target_located,'lat',sngl(lat),sngl(lat_min),sngl(lat_max),'lon',sngl(lon),sngl(lon_min),sngl(lon_max)

  if (target_located) then
    ! point in this slice

    ! loops over elements to find nearest location
    ! searches closest GLL point
    if (USE_DISTANCE_CRITERION) then
      ! loops over all elements
      do ispec = 1,nspec
        ! exclude elements that are too far from target
        dist_squared = (x_target - xyz_midpoints(1,ispec))*(x_target - xyz_midpoints(1,ispec)) &
                     + (y_target - xyz_midpoints(2,ispec))*(y_target - xyz_midpoints(2,ispec)) &
                     + (z_target - xyz_midpoints(3,ispec))*(z_target - xyz_midpoints(3,ispec))
        !  we compare squared distances instead of distances themselves to significantly speed up calculations
        if (dist_squared > typical_size_squared) cycle

        ! loop only on points inside the element
        ! exclude edges to ensure this point is not shared with other elements
        do k = 2,NGLLZ-1
          do j = 2,NGLLY-1
            do i = 2,NGLLX-1
              iglob = ibool(i,j,k,ispec)
              dist_squared = (x_target - dble(xstore(iglob)))*(x_target - dble(xstore(iglob))) &
                           + (y_target - dble(ystore(iglob)))*(y_target - dble(ystore(iglob))) &
                           + (z_target - dble(zstore(iglob)))*(z_target - dble(zstore(iglob)))

              !  keep this point if it is closer to the receiver
              !  we compare squared distances instead of distances themselves to significantly speed up calculations
              if (dist_squared < distmin_squared) then
                distmin_squared = dist_squared
                ispec_selected = ispec
                ix_initial_guess = i
                iy_initial_guess = j
                iz_initial_guess = k
              endif

            enddo
          enddo
        enddo

      ! end of loop on all the spectral elements in current slice
      enddo
    else
      ! searches through all elements
      do ispec = 1,nspec
        ! loop only on points inside the element
        ! exclude edges to ensure this point is not shared with other elements
        do k = 2,NGLLZ-1
          do j = 2,NGLLY-1
            do i = 2,NGLLX-1
              iglob = ibool(i,j,k,ispec)
              dist_squared = (x_target - dble(xstore(iglob)))*(x_target - dble(xstore(iglob))) &
                           + (y_target - dble(ystore(iglob)))*(y_target - dble(ystore(iglob))) &
                           + (z_target - dble(zstore(iglob)))*(z_target - dble(zstore(iglob)))

              !  keep this point if it is closer to the receiver
              !  we compare squared distances instead of distances themselves to significantly speed up calculations
              if (dist_squared < distmin_squared) then
                distmin_squared = dist_squared
                ispec_selected = ispec
                ix_initial_guess = i
                iy_initial_guess = j
                iz_initial_guess = k
                target_located = .true.
              endif
            enddo
          enddo
        enddo
      ! end of loop on all the spectral elements in current slice
      enddo
    endif ! USE_DISTANCE_CRITERION

  endif ! target_located


  ! ****************************************
  ! find the best (xi,eta,gamma)
  ! ****************************************
  if (target_located) then

!      ! for point sources, the location will be exactly at a GLL point
!      ! otherwise this tries to find best location
!
!      if (USE_FORCE_POINT_SOURCE) then
!        ! store xi,eta,gamma and x,y,z of point found
!        ! note: they have range [1.0d0,NGLLX/Y/Z], used for point sources
!        !          see e.g. in compute_add_sources.f90
!        xi_subset(isource_in_this_subset) = dble(ix_initial_guess)
!        eta_subset(isource_in_this_subset) = dble(iy_initial_guess)
!        gamma_subset(isource_in_this_subset) = dble(iz_initial_guess)
!
!        iglob = ibool(ix_initial_guess,iy_initial_guess, &
!            iz_initial_guess,ispec_selected)
!        xyz_found_subset(1,isource_in_this_subset) = xstore(iglob)
!        xyz_found_subset(2,isource_in_this_subset) = ystore(iglob)
!        xyz_found_subset(3,isource_in_this_subset) = zstore(iglob)
!
!        ! compute final distance between asked and found (converted to km)
!        final_distance_subset(isource_in_this_subset) = &
!          dsqrt((x_target-xyz_found_subset(1,isource_in_this_subset))**2 + &
!                (y_target-xyz_found_subset(2,isource_in_this_subset))**2 + &
!                (z_target-xyz_found_subset(3,isource_in_this_subset))**2)*R_EARTH/1000.d0
!
!      else

    ! define coordinates of the control points of the element
    do ia = 1,NGNOD
      iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec_selected)
      xelm(ia) = dble(xstore(iglob))
      yelm(ia) = dble(ystore(iglob))
      zelm(ia) = dble(zstore(iglob))
    enddo

    ! use initial guess in xi and eta
    xi = xigll(ix_initial_guess)
    eta = yigll(iy_initial_guess)
    gamma = zigll(iz_initial_guess)

    ! impose receiver exactly at the surface
    if (.not. CAN_BE_BURIED) gamma = 1.d0

    d_min_sq = HUGEVAL
    dx_min = HUGEVAL
    dy_min = HUGEVAL
    dz_min = HUGEVAL

    ! iterate to solve the non linear system
    do iter_loop = 1,NUM_ITER

      ! recompute Jacobian for the new point
      call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                              xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

      ! compute distance to target location
      dx = - (x - x_target)
      dy = - (y - y_target)
      dz = - (z - z_target)

      ! compute increments
      dxi = xix*dx + xiy*dy + xiz*dz
      deta = etax*dx + etay*dy + etaz*dz
      dgamma = gammax*dx + gammay*dy + gammaz*dz

      ! impose limit on increments
      if (abs(dxi) > 0.3d0 ) dxi = sign(1.0d0,dxi)*0.3d0
      if (abs(deta) > 0.3d0 ) deta = sign(1.0d0,deta)*0.3d0
      if (abs(dgamma) > 0.3d0 ) dgamma = sign(1.0d0,dgamma)*0.3d0

      ! update values
      xi = xi + dxi
      eta = eta + deta
      gamma = gamma + dgamma

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

    ! impose receiver exactly at the surface
    if (.not. CAN_BE_BURIED) gamma = 1.d0

    ! compute final coordinates of point found
    call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                            xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

  else
    ! point not found in this slice
    iglob = ibool(ix_initial_guess,iy_initial_guess,iz_initial_guess,ispec_selected)
    x = dble(xstore(iglob))
    y = dble(ystore(iglob))
    z = dble(zstore(iglob))
    xi = 0.d0
    eta = 0.d0
    gamma = 0.d0
  endif

  ! compute final distance between asked and found (converted to km)
  distmin_not_squared = dsqrt((x_target-x)**2 + &
                              (y_target-y)**2 + &
                              (z_target-z)**2)*R_EARTH_KM

  end subroutine locate_point
