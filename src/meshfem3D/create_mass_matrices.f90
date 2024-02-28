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

  subroutine create_mass_matrices(idoubling,ibool, &
                                  iregion_code,xstore,ystore,zstore, &
                                  NSPEC2D_TOP)

  ! creates rmassx, rmassy, rmassz and rmass_ocean_load

!
! ****************************************************************************************************
! IMPORTANT: this routine must *NOT* use flag SIMULATION_TYPE (nor SAVE_FORWARD), i.e. none of the parameters it computes
! should depend on SIMULATION_TYPE, because most users do not recompile the code nor rerun the mesher
! when switching from SIMULATION_TYPE == 1 to SIMULATION_TYPE == 3 and thus the header file created
! by this routine would become wrong in the case of a run with SIMULATION_TYPE == 3 if the code
! was compiled with SIMULATION_TYPE == 1
! ****************************************************************************************************
!

  use constants

  use meshfem_models_par, only: &
    OCEANS

  use meshfem_par, only: &
    myrank,nspec,nglob,NCHUNKS,ABSORBING_CONDITIONS, &
    ROTATION,EXACT_MASS_MATRIX_FOR_ROTATION,INCLUDE_CENTRAL_CUBE

  use regions_mesh_par, only: &
    wxgll,wygll,wzgll

  use regions_mesh_par2, only: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
    gammaxstore,gammaystore,gammazstore,rhostore,kappavstore, &
    rmassx,rmassy,rmassz,b_rmassx,b_rmassy, &
    nglob_xy

  implicit none

  integer,dimension(nspec),intent(in) :: idoubling
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool

  integer,intent(in) :: iregion_code

  ! arrays with the mesh in double precision
  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: xstore,ystore,zstore

  ! Stacey conditions put back
  integer,intent(in) :: NSPEC2D_TOP

  ! local parameters
  double precision :: weight
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl

  integer :: ispec,i,j,k,iglob

  ! initializes matrices
  !
  ! in the case of Stacey boundary conditions, add C*delta/2 contribution to the mass matrix
  ! on the Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be fictitious / unused
  !
  ! Now also handle EXACT_MASS_MATRIX_FOR_ROTATION, which requires similar corrections

  rmassx(:) = 0._CUSTOM_REAL
  rmassy(:) = 0._CUSTOM_REAL
  rmassz(:) = 0._CUSTOM_REAL

  b_rmassx(:) = 0._CUSTOM_REAL
  b_rmassy(:) = 0._CUSTOM_REAL

  ! checks if anything to do
  if (iregion_code == IREGION_TRINFINITE .or. iregion_code == IREGION_INFINITE) return

  ! first create the main standard mass matrix with no corrections

! openmp mesher
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ispec,i,j,k,iglob,weight, &
!$OMP xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl)
!$OMP DO
  do ispec = 1,nspec

    ! suppress fictitious elements in central cube
    if (idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          weight = wxgll(i)*wygll(j)*wzgll(k)
          iglob = ibool(i,j,k,ispec)

          ! compute the Jacobian
          xixl = xixstore(i,j,k,ispec)
          xiyl = xiystore(i,j,k,ispec)
          xizl = xizstore(i,j,k,ispec)
          etaxl = etaxstore(i,j,k,ispec)
          etayl = etaystore(i,j,k,ispec)
          etazl = etazstore(i,j,k,ispec)
          gammaxl = gammaxstore(i,j,k,ispec)
          gammayl = gammaystore(i,j,k,ispec)
          gammazl = gammazstore(i,j,k,ispec)

          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                          - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                          + xizl*(etaxl*gammayl-etayl*gammaxl))

          ! definition depends if region is fluid or solid
          select case (iregion_code)
          case (IREGION_CRUST_MANTLE, IREGION_INNER_CORE)
            ! distinguish between single and double precision for reals
!$OMP ATOMIC
            rmassz(iglob) = rmassz(iglob) + &
                   real(dble(rhostore(i,j,k,ispec)) * dble(jacobianl) * weight, kind=CUSTOM_REAL)

          ! fluid in outer core
          case (IREGION_OUTER_CORE)
            !checks division
            if (kappavstore(i,j,k,ispec) <= 0.0_CUSTOM_REAL) stop 'Error invalid kappav in outer core mass matrix'

            ! no anisotropy in the fluid, use kappav

            ! distinguish between single and double precision for reals
!$OMP ATOMIC
            rmassz(iglob) = rmassz(iglob) + &
                   real(dble(jacobianl) * weight * dble(rhostore(i,j,k,ispec)) / dble(kappavstore(i,j,k,ispec)), &
                        kind=CUSTOM_REAL)

          case default
            call exit_MPI(myrank,'wrong region code in create_mass_matrix')

          end select

        enddo
      enddo
    enddo
  enddo ! of loop on ispec
!$OMP ENDDO
!$OMP END PARALLEL

  ! copy the initial mass matrix if needed
  if (nglob_xy == nglob) then
    rmassx(:) = rmassz(:)
    rmassy(:) = rmassz(:)

    b_rmassx(:) = rmassz(:)
    b_rmassy(:) = rmassz(:)
  endif

  ! then make the corrections to the copied mass matrices if needed
  if (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
    call create_mass_matrices_rotation(ibool,idoubling,iregion_code)
  endif

  ! absorbing boundaries
  ! add C*deltat/2 contribution to the mass matrices on the Stacey edges
  if (NCHUNKS /= 6 .and. ABSORBING_CONDITIONS) then
    call create_mass_matrices_Stacey(ibool,iregion_code)
  endif

  ! check that mass matrix is positive
  if (iregion_code == IREGION_INNER_CORE .and. INCLUDE_CENTRAL_CUBE) then
    ! note: in fictitious elements mass matrix is still zero
    if (minval(rmassz(:)) < 0._CUSTOM_REAL) call exit_MPI(myrank,'negative rmassz matrix term')
    ! check that the additional mass matrices are positive, if they exist
    if (nglob_xy == nglob) then
      if (minval(rmassx) < 0._CUSTOM_REAL) call exit_MPI(myrank,'negative rmassx matrix term')
      if (minval(rmassy) < 0._CUSTOM_REAL) call exit_MPI(myrank,'negative rmassy matrix term')
      if (minval(b_rmassx) < 0._CUSTOM_REAL) call exit_MPI(myrank,'negative b_rmassx matrix term')
      if (minval(b_rmassy) < 0._CUSTOM_REAL) call exit_MPI(myrank,'negative b_rmassy matrix term')
    endif
  else
    ! no fictitious elements, mass matrix must be strictly positive
    if (minval(rmassz(:)) <= 0._CUSTOM_REAL) call exit_MPI(myrank,'negative rmassz matrix term')
    ! check that the additional mass matrices are strictly positive, if they exist
    if (nglob_xy == nglob) then
      if (minval(rmassx) <= 0._CUSTOM_REAL) call exit_MPI(myrank,'negative rmassx matrix term')
      if (minval(rmassy) <= 0._CUSTOM_REAL) call exit_MPI(myrank,'negative rmassy matrix term')
      if (minval(b_rmassx) <= 0._CUSTOM_REAL) call exit_MPI(myrank,'negative b_rmassx matrix term')
      if (minval(b_rmassy) <= 0._CUSTOM_REAL) call exit_MPI(myrank,'negative b_rmassy matrix term')
    endif
  endif

  ! save ocean load mass matrix as well if oceans
  if (OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) then
    call create_mass_matrices_ocean_load(ibool,xstore,ystore,zstore,NSPEC2D_TOP)
  endif

  end subroutine create_mass_matrices


!
!-------------------------------------------------------------------------------------------------
!

  subroutine create_mass_matrices_rotation(ibool,idoubling,iregion_code)

! in the case of rotation, add C*deltat/2 contribution to the mass matrix
! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
!
! only called in case of (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION)
  use constants

  use meshfem_par, only: &
    myrank,DT,nspec

  use regions_mesh_par, only: &
    wxgll,wygll,wzgll

  use regions_mesh_par2, only: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
    gammaxstore,gammaystore,gammazstore, &
    rmassx,rmassy,b_rmassx,b_rmassy

  use shared_parameters, only: UNDO_ATTENUATION,HOURS_PER_DAY,SECONDS_PER_HOUR,RHOAV

  implicit none

  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  integer,dimension(nspec) :: idoubling

  integer :: iregion_code

  ! local parameters
  double precision :: weight

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: deltat,two_omega_earth_dt,b_two_omega_earth_dt,scale_t_inv

  integer :: ispec,i,j,k,iglob

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    creates exact mass matrix for rotation'
    call flush_IMAIN()
  endif

  ! use the non-dimensional time step to make the mass matrix correction
  deltat = real(DT*dsqrt(PI*GRAV*RHOAV), kind=CUSTOM_REAL)

  scale_t_inv = real(dsqrt(PI*GRAV*RHOAV), kind=CUSTOM_REAL)

  ! distinguish between single and double precision for reals
  two_omega_earth_dt = real(2.d0 * TWO_PI / (HOURS_PER_DAY * SECONDS_PER_HOUR * scale_t_inv) * deltat, kind=CUSTOM_REAL)

  ! reconstructed wavefield
  if (UNDO_ATTENUATION) then
    ! spinning forward
    b_two_omega_earth_dt = two_omega_earth_dt
  else
    ! spinning backward to reconstruct wavefield
    b_two_omega_earth_dt = - two_omega_earth_dt
  endif

  ! definition depends if region is fluid or solid
  select case (iregion_code)
  case (IREGION_CRUST_MANTLE, IREGION_INNER_CORE)
! openmp mesher
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ispec,i,j,k,iglob,weight, &
!$OMP xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl)
!$OMP DO
    do ispec = 1,nspec

      ! suppress fictitious elements in central cube
      if (idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX

            weight = wxgll(i)*wygll(j)*wzgll(k)
            iglob = ibool(i,j,k,ispec)

            ! compute the Jacobian
            xixl = xixstore(i,j,k,ispec)
            xiyl = xiystore(i,j,k,ispec)
            xizl = xizstore(i,j,k,ispec)
            etaxl = etaxstore(i,j,k,ispec)
            etayl = etaystore(i,j,k,ispec)
            etazl = etazstore(i,j,k,ispec)
            gammaxl = gammaxstore(i,j,k,ispec)
            gammayl = gammaystore(i,j,k,ispec)
            gammazl = gammazstore(i,j,k,ispec)

            jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                            - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                            + xizl*(etaxl*gammayl-etayl*gammaxl))

            ! distinguish between single and double precision for reals
!$OMP ATOMIC
            rmassx(iglob) = rmassx(iglob) &
                - two_omega_earth_dt * 0.5_CUSTOM_REAL*real(dble(jacobianl) * weight, kind=CUSTOM_REAL)
!$OMP ATOMIC
            rmassy(iglob) = rmassy(iglob) &
                + two_omega_earth_dt * 0.5_CUSTOM_REAL*real(dble(jacobianl) * weight, kind=CUSTOM_REAL)

!$OMP ATOMIC
            b_rmassx(iglob) = b_rmassx(iglob) &
                - b_two_omega_earth_dt * 0.5_CUSTOM_REAL*real(dble(jacobianl) * weight, kind=CUSTOM_REAL)
!$OMP ATOMIC
            b_rmassy(iglob) = b_rmassy(iglob) &
                + b_two_omega_earth_dt * 0.5_CUSTOM_REAL*real(dble(jacobianl) * weight, kind=CUSTOM_REAL)
          enddo
        enddo
      enddo
    enddo ! of loop on ispec
!$OMP ENDDO
!$OMP END PARALLEL

  case (IREGION_OUTER_CORE)
    ! nothing to do
    continue

  case default
    call exit_MPI(myrank,'wrong region code in create_mass_matrices_rotation')
  end select

  end subroutine create_mass_matrices_rotation


!
!-------------------------------------------------------------------------------------------------
!

  subroutine create_mass_matrices_Stacey(ibool,iregion_code)

! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
! on Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix

  use constants

  use meshfem_par, only: &
    myrank,DT,nspec,nglob, &
    ROTATION,EXACT_MASS_MATRIX_FOR_ROTATION,ABSORBING_CONDITIONS

  use regions_mesh_par2, only: &
    rmassx,rmassy,rmassz,b_rmassx,b_rmassy, &
    rho_vp,rho_vs, &
    nglob_xy

  ! absorb
  use regions_mesh_par2, only: num_abs_boundary_faces, &
    abs_boundary_ispec,abs_boundary_npoin, &
    abs_boundary_ijk,abs_boundary_normal,abs_boundary_jacobian2Dw

  use shared_parameters, only: RHOAV,USE_LDDRK

  implicit none

  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool

  integer,intent(in) :: iregion_code

  ! local parameters
  double precision :: weight
  real(kind=CUSTOM_REAL) :: deltat,deltatover2
  real(kind=CUSTOM_REAL) :: tx,ty,tz,sn
  real(kind=CUSTOM_REAL) :: nx,ny,nz,vn
  integer :: ispec,i,j,k,iglob
  integer :: igll,iface

!--- DK and Zhinan Xie: add C Delta_t / 2 contribution to the mass matrix
!--- DK and Zhinan Xie: in the case of Clayton-Engquist absorbing boundaries;
!--- DK and Zhinan Xie: see for instance the book of Hughes (1987) chapter 9.
!--- DK and Zhinan Xie: IMPORTANT: note that this implies that we must have two different mass matrices,
!--- DK and Zhinan Xie: one per component of the wave field i.e. one per spatial dimension.
!--- DK and Zhinan Xie: This was also suggested by Jean-Paul Ampuero in 2003.

  ! checks if anything to do
  if (.not. ABSORBING_CONDITIONS) return

  ! no Stacey boundary on inner core, only crust/mantle and outer core
  if (iregion_code == IREGION_INNER_CORE) return

  ! saftey check
  ! crust/mantle region must have nglob_xy set to have 3 different rmassx/rmassy/rmassz
  if (iregion_code == IREGION_CRUST_MANTLE .and. (nglob_xy /= nglob)) &
    stop 'Invalid nglob_xy for crust/mantle Stacey boundary'

  ! note: for LDDRK, the time scheme needs no mass matrix contribution due to the absorbing boundary term.
  !       the additional contribution comes from the Newmark formulation and only needs to be added in those cases.
  if (USE_LDDRK) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    updates mass matrix with Stacey boundary corrections'
    call flush_IMAIN()
  endif

  ! checks if we have absorbing boundary arrays
  if (.not. allocated(abs_boundary_ispec) ) call exit_MPI(myrank,'Error Stacey array not allocated')

  ! use the non-dimensional time step to make the mass matrix correction
  deltat = real(DT*dsqrt(PI*GRAV*RHOAV), kind=CUSTOM_REAL)
  deltatover2 = real(0.5d0*deltat, kind=CUSTOM_REAL)

  ! adds contributions to mass matrix to stabilize Stacey conditions
  select case (iregion_code)
  case (IREGION_CRUST_MANTLE)
    do iface = 1,num_abs_boundary_faces

      ispec = abs_boundary_ispec(iface)

      do igll = 1,abs_boundary_npoin(iface)
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        nx = abs_boundary_normal(1,igll,iface)
        ny = abs_boundary_normal(2,igll,iface)
        nz = abs_boundary_normal(3,igll,iface)

        weight = abs_boundary_jacobian2Dw(igll,iface)

        iglob = ibool(i,j,k,ispec)

        vn = deltatover2 * (nx+ny+nz)

        tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(deltatover2-vn*nx)
        ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(deltatover2-vn*ny)
        tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(deltatover2-vn*nz)

        rmassx(iglob) = rmassx(iglob) + real(tx*weight, kind=CUSTOM_REAL)
        rmassy(iglob) = rmassy(iglob) + real(ty*weight, kind=CUSTOM_REAL)
        rmassz(iglob) = rmassz(iglob) + real(tz*weight, kind=CUSTOM_REAL)

        if (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
          b_rmassx(iglob) = b_rmassx(iglob) + real(tx*weight, kind=CUSTOM_REAL)
          b_rmassy(iglob) = b_rmassy(iglob) + real(ty*weight, kind=CUSTOM_REAL)
        endif
      enddo
    enddo

    ! check that mass matrix is positive
    if (minval(rmassx(:)) <= 0.) call exit_MPI(myrank,'negative rmassx matrix term')
    if (minval(rmassy(:)) <= 0.) call exit_MPI(myrank,'negative rmassy matrix term')

  case (IREGION_OUTER_CORE)
    do iface = 1,num_abs_boundary_faces

      ispec = abs_boundary_ispec(iface)

      do igll = 1,abs_boundary_npoin(iface)
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)

        weight = abs_boundary_jacobian2Dw(igll,iface)

        iglob = ibool(i,j,k,ispec)

        !checks division
        if (rho_vp(i,j,k,ispec) <= 0.0_CUSTOM_REAL) stop 'Error invalid rho_vp in outer core stacey xmin'

        sn = deltatover2/rho_vp(i,j,k,ispec)

        rmassz(iglob) = rmassz(iglob) + real(weight*sn, kind=CUSTOM_REAL)
      enddo
    enddo

  case (IREGION_INNER_CORE)
    continue

  case default
    call exit_MPI(myrank,'wrong region code in create_mass_matrices_Stacey')

  end select

  end subroutine create_mass_matrices_Stacey

!
!-------------------------------------------------------------------------------------------------
!

  subroutine create_mass_matrices_ocean_load(ibool,xstore,ystore,zstore,NSPEC2D_TOP)

  use constants

  use meshfem_models_par, only: &
    OCEANS,TOPOGRAPHY,ELLIPTICITY, &
    ibathy_topo,CASE_3D

  use meshfem_par, only: &
    myrank,RHO_OCEANS,nspec

  use regions_mesh_par, only: &
    wxgll,wygll

  use regions_mesh_par2, only: &
    rmassz,rmass_ocean_load, &
    ibelm_top,jacobian2D_top

  use shared_parameters, only: R_PLANET

  implicit none

  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool

  ! arrays with the mesh in double precision
  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: xstore,ystore,zstore

  integer,intent(in) :: NSPEC2D_TOP

  ! local parameters
  double precision :: x,y,z,weight
  double precision :: r,lat,lon
  double precision :: elevation,height_oceans

  integer :: ispec,i,j,k,iglob,ispec2D

  logical :: do_ocean_load

  ! checks if anything to do
  if (.not. OCEANS) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    updates mass matrix with ocean load'
    call flush_IMAIN()
  endif

  ! initializes
  do_ocean_load = .false.

  ! note: old version (5.1.5)
  ! only for models where 3D crustal stretching was used (even without topography?)
  if (USE_OLD_VERSION_5_1_5_FORMAT) then
    if (CASE_3D) then
      do_ocean_load = .true.
    endif
  else
    ! note: new version:
    ! for 3D Earth with topography, compute local height of oceans
    if (TOPOGRAPHY) then
      do_ocean_load = .true.
    endif
  endif

  ! create ocean load mass matrix for degrees of freedom at ocean bottom
  rmass_ocean_load(:) = 0._CUSTOM_REAL

  ! add contribution of the oceans
  ! for surface elements exactly at the top of the crust (ocean bottom)

! openmp mesher
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ispec2D,ispec,i,j,iglob,x,y,z,r,lat,lon,elevation,height_oceans,weight)
!$OMP DO
  do ispec2D = 1,NSPEC2D_TOP

    ! gets spectral element index
    ispec = ibelm_top(ispec2D)

    ! assumes elements are ordered such that k == NGLLZ is the top surface
    k = NGLLZ

    ! loops over surface points
    do j = 1,NGLLY
      do i = 1,NGLLX

        ! for 3D Earth with topography, compute local height of oceans
        if (do_ocean_load) then

          ! get coordinates of current point
          x = xstore(i,j,k,ispec)
          y = ystore(i,j,k,ispec)
          z = zstore(i,j,k,ispec)

          ! map to geographic latitude and longitude (in degrees) for bathymetry routine
          ! note: at this point, the mesh can be elliptical, depending on the Par_file flag choosen
          call xyz_2_rlatlon_dble(x,y,z,r,lat,lon,ELLIPTICITY)

          ! compute elevation at current point
          call get_topo_bathy(lat,lon,elevation,ibathy_topo)

          ! non-dimensionalize the elevation, which is in meters
          ! and suppress positive elevation, which means no oceans
          if (elevation >= - MINIMUM_THICKNESS_3D_OCEANS) then
            height_oceans = 0.d0
          else
            height_oceans = dabs(elevation) / R_PLANET
          endif

        else
          ! if 1D Earth, use oceans of constant thickness everywhere
          height_oceans = THICKNESS_OCEANS_PREM
        endif

        ! take into account inertia of water column
        weight = wxgll(i) * wygll(j) &
                  * dble(jacobian2D_top(i,j,ispec2D)) &
                  * dble(RHO_OCEANS) * height_oceans

        ! gets global point index
        iglob = ibool(i,j,k,ispec)

        ! distinguish between single and double precision for reals
!$OMP ATOMIC
        rmass_ocean_load(iglob) = rmass_ocean_load(iglob) + real(weight, kind=CUSTOM_REAL)

      enddo
    enddo
  enddo
!$OMP ENDDO
!$OMP END PARALLEL

  ! add regular mass matrix to ocean load contribution
  rmass_ocean_load(:) = rmass_ocean_load(:) + rmassz(:)

  end subroutine create_mass_matrices_ocean_load

