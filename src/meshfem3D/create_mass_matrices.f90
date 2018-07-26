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

  subroutine create_mass_matrices(idoubling,ibool, &
                                  iregion_code,xstore,ystore,zstore, &
                                  NSPEC2D_TOP,NSPEC2D_BOTTOM)

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

  use meshfem3D_models_par, only: &
    OCEANS

  use meshfem3D_par, only: &
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
  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  ! Stacey conditions put back
  integer :: NSPEC2D_TOP,NSPEC2D_BOTTOM

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

!----------------------------------------------------------------

! first create the main standard mass matrix with no corrections
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
            rmassz(iglob) = rmassz(iglob) + &
                   real(dble(rhostore(i,j,k,ispec)) * dble(jacobianl) * weight, kind=CUSTOM_REAL)

          ! fluid in outer core
          case (IREGION_OUTER_CORE)

            ! no anisotropy in the fluid, use kappav

            ! distinguish between single and double precision for reals
            rmassz(iglob) = rmassz(iglob) + &
                   real(dble(jacobianl) * weight * dble(rhostore(i,j,k,ispec)) / dble(kappavstore(i,j,k,ispec)), &
                        kind=CUSTOM_REAL)

          case default
            call exit_MPI(myrank,'wrong region code')

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
    call create_mass_matrices_Stacey(ibool,iregion_code,NSPEC2D_BOTTOM)
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

  use meshfem3D_par, only: &
    myrank,DT,nspec

  use regions_mesh_par, only: &
    wxgll,wygll,wzgll

  use regions_mesh_par2, only: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
    gammaxstore,gammaystore,gammazstore, &
    rmassx,rmassy,b_rmassx,b_rmassy

  use shared_parameters, only: UNDO_ATTENUATION

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
            rmassx(iglob) = rmassx(iglob) &
                - two_omega_earth_dt * 0.5_CUSTOM_REAL*real(dble(jacobianl) * weight, kind=CUSTOM_REAL)
            rmassy(iglob) = rmassy(iglob) &
                + two_omega_earth_dt * 0.5_CUSTOM_REAL*real(dble(jacobianl) * weight, kind=CUSTOM_REAL)

            b_rmassx(iglob) = b_rmassx(iglob) &
                - b_two_omega_earth_dt * 0.5_CUSTOM_REAL*real(dble(jacobianl) * weight, kind=CUSTOM_REAL)
            b_rmassy(iglob) = b_rmassy(iglob) &
                + b_two_omega_earth_dt * 0.5_CUSTOM_REAL*real(dble(jacobianl) * weight, kind=CUSTOM_REAL)
          enddo
        enddo
      enddo
    enddo ! of loop on ispec
!$OMP ENDDO
!$OMP END PARALLEL

  end select


  end subroutine create_mass_matrices_rotation


!
!-------------------------------------------------------------------------------------------------
!

  subroutine create_mass_matrices_Stacey(ibool,iregion_code,NSPEC2D_BOTTOM)

! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
! on Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix

  use constants

  use meshfem3D_par, only: &
    myrank,DT,NCHUNKS,ichunk,nspec, &
    ROTATION,EXACT_MASS_MATRIX_FOR_ROTATION

  use regions_mesh_par, only: &
    wxgll,wygll,wzgll

  use regions_mesh_par2, only: &
    rmassx,rmassy,rmassz,b_rmassx,b_rmassy, &
    ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom, &
    normal_xmin,normal_xmax,normal_ymin,normal_ymax, &
    jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax, &
    jacobian2D_bottom, &
    rho_vp,rho_vs, &
    nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
    nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta

  implicit none

  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  integer :: iregion_code

  ! Stacey conditions
  integer :: NSPEC2D_BOTTOM

  ! local parameters
  double precision :: weight
  double precision, dimension(NGLLX,NGLLY) :: wgllwgll_xy
  double precision, dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  double precision, dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  real(kind=CUSTOM_REAL) :: deltat,deltatover2
  real(kind=CUSTOM_REAL) :: tx,ty,tz,sn
  real(kind=CUSTOM_REAL) :: nx,ny,nz,vn

  integer :: ispec,i,j,k,iglob
  integer :: ispec2D

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    updates mass matrix with Stacey boundary corrections'
    call flush_IMAIN()
  endif

  ! checks if we have absorbing boundary arrays
  if (.not. allocated(nimin) ) call exit_MPI(myrank,'Error Stacey array not allocated')

  ! use the non-dimensional time step to make the mass matrix correction
  deltat = real(DT*dsqrt(PI*GRAV*RHOAV), kind=CUSTOM_REAL)
  deltatover2 = real(0.5d0*deltat, kind=CUSTOM_REAL)

  ! weights on surfaces
  do i = 1,NGLLX
    do j = 1,NGLLY
       wgllwgll_xy(i,j) = wxgll(i)*wygll(j)
    enddo
  enddo
  do i = 1,NGLLX
    do k = 1,NGLLZ
       wgllwgll_xz(i,k) = wxgll(i)*wzgll(k)
    enddo
  enddo
  do j = 1,NGLLY
    do k = 1,NGLLZ
       wgllwgll_yz(j,k) = wygll(j)*wzgll(k)
    enddo
  enddo

  ! adds contributions to mass matrix to stabilize Stacey conditions
  select case (iregion_code)
  case (IREGION_CRUST_MANTLE)

    rmassx(:) = rmassz(:)
    rmassy(:) = rmassz(:)

    !   xmin
    ! if two chunks exclude this face for one of them
    if (NCHUNKS == 1 .or. ichunk == CHUNK_AC) then

       do ispec2D = 1,nspec2D_xmin

          ispec=ibelm_xmin(ispec2D)

          ! exclude elements that are not on absorbing edges
          if (nkmin_xi(1,ispec2D) == 0 .or. njmin(1,ispec2D) == 0) cycle

          i = 1
          do k = nkmin_xi(1,ispec2D),NGLLZ
             do j = njmin(1,ispec2D),njmax(1,ispec2D)
                iglob=ibool(i,j,k,ispec)

                nx = normal_xmin(1,j,k,ispec2D)
                ny = normal_xmin(2,j,k,ispec2D)
                nz = normal_xmin(3,j,k,ispec2D)

                vn = deltatover2*(nx+ny+nz)

                tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(deltatover2-vn*nx)
                ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(deltatover2-vn*ny)
                tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(deltatover2-vn*nz)

                weight = jacobian2D_xmin(j,k,ispec2D)*wgllwgll_yz(j,k)

                rmassx(iglob) = rmassx(iglob) + real(tx*weight, kind=CUSTOM_REAL)
                rmassy(iglob) = rmassy(iglob) + real(ty*weight, kind=CUSTOM_REAL)
                rmassz(iglob) = rmassz(iglob) + real(tz*weight, kind=CUSTOM_REAL)
                if (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
                  b_rmassx(iglob) = b_rmassx(iglob) + real(tx*weight, kind=CUSTOM_REAL)
                  b_rmassy(iglob) = b_rmassy(iglob) + real(ty*weight, kind=CUSTOM_REAL)
                endif
             enddo
          enddo
       enddo

    endif ! NCHUNKS == 1 .or. ichunk == CHUNK_AC

    !   xmax
    ! if two chunks exclude this face for one of them
    if (NCHUNKS == 1 .or. ichunk == CHUNK_AB) then

       do ispec2D = 1,nspec2D_xmax

          ispec=ibelm_xmax(ispec2D)

          ! exclude elements that are not on absorbing edges
          if (nkmin_xi(2,ispec2D) == 0 .or. njmin(2,ispec2D) == 0) cycle

          i=NGLLX
          do k=nkmin_xi(2,ispec2D),NGLLZ
             do j=njmin(2,ispec2D),njmax(2,ispec2D)
                iglob=ibool(i,j,k,ispec)

                nx = normal_xmax(1,j,k,ispec2D)
                ny = normal_xmax(2,j,k,ispec2D)
                nz = normal_xmax(3,j,k,ispec2D)

                vn = deltatover2*(nx+ny+nz)

                tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(deltatover2-vn*nx)
                ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(deltatover2-vn*ny)
                tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(deltatover2-vn*nz)

                weight = jacobian2D_xmax(j,k,ispec2D)*wgllwgll_yz(j,k)

                rmassx(iglob) = rmassx(iglob) + real(tx*weight, kind=CUSTOM_REAL)
                rmassy(iglob) = rmassy(iglob) + real(ty*weight, kind=CUSTOM_REAL)
                rmassz(iglob) = rmassz(iglob) + real(tz*weight, kind=CUSTOM_REAL)
                if (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
                  b_rmassx(iglob) = b_rmassx(iglob) + real(tx*weight, kind=CUSTOM_REAL)
                  b_rmassy(iglob) = b_rmassy(iglob) + real(ty*weight, kind=CUSTOM_REAL)
                endif
             enddo
          enddo
       enddo

    endif ! NCHUNKS == 1 .or. ichunk == CHUNK_AB

    !   ymin
    do ispec2D = 1,nspec2D_ymin

       ispec=ibelm_ymin(ispec2D)

       ! exclude elements that are not on absorbing edges
       if (nkmin_eta(1,ispec2D) == 0 .or. nimin(1,ispec2D) == 0) cycle

       j = 1
       do k = nkmin_eta(1,ispec2D),NGLLZ
          do i = nimin(1,ispec2D),nimax(1,ispec2D)
            iglob=ibool(i,j,k,ispec)

             nx = normal_ymin(1,i,k,ispec2D)
             ny = normal_ymin(2,i,k,ispec2D)
             nz = normal_ymin(3,i,k,ispec2D)

             vn = deltatover2*(nx+ny+nz)

             tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(deltatover2-vn*nx)
             ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(deltatover2-vn*ny)
             tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(deltatover2-vn*nz)

             weight = jacobian2D_ymin(i,k,ispec2D)*wgllwgll_xz(i,k)

             rmassx(iglob) = rmassx(iglob) + real(tx*weight, kind=CUSTOM_REAL)
             rmassy(iglob) = rmassy(iglob) + real(ty*weight, kind=CUSTOM_REAL)
             rmassz(iglob) = rmassz(iglob) + real(tz*weight, kind=CUSTOM_REAL)
             if (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
               b_rmassx(iglob) = b_rmassx(iglob) + real(tx*weight, kind=CUSTOM_REAL)
               b_rmassy(iglob) = b_rmassy(iglob) + real(ty*weight, kind=CUSTOM_REAL)
             endif
          enddo
       enddo
    enddo

    !   ymax
    do ispec2D = 1,nspec2D_ymax

       ispec=ibelm_ymax(ispec2D)

       ! exclude elements that are not on absorbing edges
       if (nkmin_eta(2,ispec2D) == 0 .or. nimin(2,ispec2D) == 0) cycle

       j=NGLLY
       do k=nkmin_eta(2,ispec2D),NGLLZ
          do i=nimin(2,ispec2D),nimax(2,ispec2D)
             iglob=ibool(i,j,k,ispec)

             nx = normal_ymax(1,i,k,ispec2D)
             ny = normal_ymax(2,i,k,ispec2D)
             nz = normal_ymax(3,i,k,ispec2D)

             vn = deltatover2*(nx+ny+nz)

             tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(deltatover2-vn*nx)
             ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(deltatover2-vn*ny)
             tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(deltatover2-vn*nz)

             weight = jacobian2D_ymax(i,k,ispec2D)*wgllwgll_xz(i,k)

             rmassx(iglob) = rmassx(iglob) + real(tx*weight, kind=CUSTOM_REAL)
             rmassy(iglob) = rmassy(iglob) + real(ty*weight, kind=CUSTOM_REAL)
             rmassz(iglob) = rmassz(iglob) + real(tz*weight, kind=CUSTOM_REAL)
             if (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
               b_rmassx(iglob) = b_rmassx(iglob) + real(tx*weight, kind=CUSTOM_REAL)
               b_rmassy(iglob) = b_rmassy(iglob) + real(ty*weight, kind=CUSTOM_REAL)
             endif
          enddo
       enddo
    enddo

    ! check that mass matrix is positive
    if (minval(rmassx(:)) <= 0.) call exit_MPI(myrank,'negative rmassx matrix term')
    if (minval(rmassy(:)) <= 0.) call exit_MPI(myrank,'negative rmassy matrix term')

  case (IREGION_OUTER_CORE)

    !   xmin
    ! if two chunks exclude this face for one of them
    if (NCHUNKS == 1 .or. ichunk == CHUNK_AC) then

       do ispec2D = 1,nspec2D_xmin

          ispec=ibelm_xmin(ispec2D)

          ! exclude elements that are not on absorbing edges
          if (nkmin_xi(1,ispec2D) == 0 .or. njmin(1,ispec2D) == 0) cycle

          i = 1
          do k = nkmin_xi(1,ispec2D),NGLLZ
             do j = njmin(1,ispec2D),njmax(1,ispec2D)
                iglob=ibool(i,j,k,ispec)

                sn = deltatover2/rho_vp(i,j,k,ispec)

                weight = jacobian2D_xmin(j,k,ispec2D)*wgllwgll_yz(j,k)

                rmassz(iglob) = rmassz(iglob) + real(weight*sn, kind=CUSTOM_REAL)
             enddo
          enddo
       enddo

    endif ! NCHUNKS == 1 .or. ichunk == CHUNK_AC

    !   xmax
    ! if two chunks exclude this face for one of them
    if (NCHUNKS == 1 .or. ichunk == CHUNK_AB) then

       do ispec2D = 1,nspec2D_xmax

          ispec=ibelm_xmax(ispec2D)

          ! exclude elements that are not on absorbing edges
          if (nkmin_xi(2,ispec2D) == 0 .or. njmin(2,ispec2D) == 0) cycle

          i=NGLLX
          do k=nkmin_xi(2,ispec2D),NGLLZ
             do j=njmin(2,ispec2D),njmax(2,ispec2D)
                iglob=ibool(i,j,k,ispec)

                sn = deltatover2/rho_vp(i,j,k,ispec)

                weight = jacobian2D_xmax(j,k,ispec2D)*wgllwgll_yz(j,k)

                rmassz(iglob) = rmassz(iglob) + real(weight*sn, kind=CUSTOM_REAL)
             enddo
          enddo
       enddo

    endif ! NCHUNKS == 1 .or. ichunk == CHUNK_AB

    !   ymin
    do ispec2D = 1,nspec2D_ymin

       ispec=ibelm_ymin(ispec2D)

       ! exclude elements that are not on absorbing edges
       if (nkmin_eta(1,ispec2D) == 0 .or. nimin(1,ispec2D) == 0) cycle

       j = 1
       do k = nkmin_eta(1,ispec2D),NGLLZ
          do i = nimin(1,ispec2D),nimax(1,ispec2D)
             iglob=ibool(i,j,k,ispec)

             sn = deltatover2/rho_vp(i,j,k,ispec)

             weight = jacobian2D_ymin(i,k,ispec2D)*wgllwgll_xz(i,k)

             rmassz(iglob) = rmassz(iglob) + real(weight*sn, kind=CUSTOM_REAL)
          enddo
       enddo
    enddo

    !   ymax
    do ispec2D = 1,nspec2D_ymax

       ispec=ibelm_ymax(ispec2D)

       ! exclude elements that are not on absorbing edges
       if (nkmin_eta(2,ispec2D) == 0 .or. nimin(2,ispec2D) == 0) cycle

       j=NGLLY
       do k=nkmin_eta(2,ispec2D),NGLLZ
          do i=nimin(2,ispec2D),nimax(2,ispec2D)
             iglob=ibool(i,j,k,ispec)

             sn = deltatover2/rho_vp(i,j,k,ispec)

             weight = jacobian2D_ymax(i,k,ispec2D)*wgllwgll_xz(i,k)

             rmassz(iglob) = rmassz(iglob) + real(weight*sn, kind=CUSTOM_REAL)
          enddo
       enddo
    enddo

    !   bottom (zmin)
    do ispec2D = 1,NSPEC2D_BOTTOM

       ispec=ibelm_bottom(ispec2D)

       k = 1
       do j = 1,NGLLY
          do i = 1,NGLLX
             iglob=ibool(i,j,k,ispec)

             sn = deltatover2/rho_vp(i,j,k,ispec)

             weight = jacobian2D_bottom(i,j,ispec2D)*wgllwgll_xy(i,j)

             rmassz(iglob) = rmassz(iglob) + real(weight*sn, kind=CUSTOM_REAL)
          enddo
       enddo
    enddo

  case (IREGION_INNER_CORE)
    continue

  case default
    call exit_MPI(myrank,'wrong region code')

  end select

  end subroutine create_mass_matrices_Stacey

!
!-------------------------------------------------------------------------------------------------
!

  subroutine create_mass_matrices_ocean_load(ibool,xstore,ystore,zstore,NSPEC2D_TOP)

  use constants

  use meshfem3D_models_par, only: &
    OCEANS,TOPOGRAPHY,ibathy_topo,CASE_3D

  use meshfem3D_par, only: &
    myrank,RHO_OCEANS,nspec

  use regions_mesh_par, only: &
    wxgll,wygll

  use regions_mesh_par2, only: &
    rmassz,rmass_ocean_load, &
    ibelm_top,jacobian2D_top

  implicit none

  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  ! arrays with the mesh in double precision
  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  integer :: NSPEC2D_TOP

  ! local parameters
  double precision :: x,y,z,r,theta,phi,weight
  double precision :: lat,lon
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
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ispec2D,ispec,i,j,iglob,x,y,z,r,theta,phi,lat,lon,elevation,height_oceans,weight)
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

          ! map to latitude and longitude for bathymetry routine
          ! slightly move points to avoid roundoff problem when exactly on the polar axis
          call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)

          if (.not. USE_OLD_VERSION_5_1_5_FORMAT) then
            ! adds small margins
  !! DK DK: added a test to only do this if we are on the axis
            if (abs(theta) > 89.99d0) then
              theta = theta + 0.0000001d0
              phi = phi + 0.0000001d0
            endif
          endif

          call reduce(theta,phi)

          ! converts the geocentric colatitude to a geographic colatitude
          ! note: bathymetry is given in geographic lat/lon
          !       (i.e., latitude with respect to reference ellipsoid)
          !       we will need convert the geocentric positions here to geographic ones
          if (USE_OLD_VERSION_5_1_5_FORMAT) then
            ! always converts
            theta = PI_OVER_TWO - datan(1.006760466d0*dcos(theta)/dmax1(TINYVAL,dsin(theta)))
          else
            ! will take flag ASSUME_PERFECT_SPHERE into account
            call geocentric_2_geographic_dble(theta,theta)
          endif

          ! get geographic latitude and longitude in degrees
          lat = (PI_OVER_TWO-theta)*RADIANS_TO_DEGREES
          lon = phi * RADIANS_TO_DEGREES

          ! compute elevation at current point
          call get_topo_bathy(lat,lon,elevation,ibathy_topo)

          ! non-dimensionalize the elevation, which is in meters
          ! and suppress positive elevation, which means no oceans
          if (elevation >= - MINIMUM_THICKNESS_3D_OCEANS) then
            height_oceans = 0.d0
          else
            height_oceans = dabs(elevation) / R_EARTH
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
        rmass_ocean_load(iglob) = rmass_ocean_load(iglob) + real(weight, kind=CUSTOM_REAL)

      enddo
    enddo
  enddo
!$OMP ENDDO
!$OMP END PARALLEL

  ! add regular mass matrix to ocean load contribution
  rmass_ocean_load(:) = rmass_ocean_load(:) + rmassz(:)

  end subroutine create_mass_matrices_ocean_load

