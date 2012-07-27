!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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

  subroutine create_mass_matrices(myrank,nspec,idoubling,wxgll,wygll,wzgll,ibool, &
                          nspec_actually,xixstore,xiystore,xizstore, &
                          etaxstore,etaystore,etazstore, &
                          gammaxstore,gammaystore,gammazstore, &
                          iregion_code,rhostore,kappavstore, &
                          nglob_xy,nglob,prname, &
                          rmassx,rmassy,rmassz, &
                          nglob_oceans,rmass_ocean_load, &
                          xstore,ystore,zstore,RHO_OCEANS, &
                          NSPEC2D_TOP,NSPEC2D_BOTTOM,NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                          ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                          nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                          normal_xmin,normal_xmax,normal_ymin,normal_ymax, &
                          rho_vp,rho_vs,nspec_stacey, &
                          jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax, &
                          jacobian2D_bottom,jacobian2D_top)

  ! creates rmassx, rmassy, rmassz and rmass_ocean_load

  use meshfem3D_models_par
  use meshfem3D_par,only: DT, NCHUNKS, ABSORBING_CONDITIONS, ichunk

  implicit none

  integer :: myrank,nspec
  integer :: idoubling(nspec)
  integer :: ibool(NGLLX,NGLLY,NGLLZ,nspec)
  integer :: nspec_actually

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_actually) ::  &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

  integer :: iregion_code,nglob_xy,nglob

  ! mass matrices
  ! add C*deltat/2 contribution to the mass matrices on Stacey edges
  real(kind=CUSTOM_REAL), dimension(nglob_xy) :: rmassx,rmassy
  real(kind=CUSTOM_REAL), dimension(nglob)    :: rmassz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: rhostore,kappavstore

  ! ocean mass matrix
  integer :: nglob_oceans
  real(kind=CUSTOM_REAL), dimension(nglob_oceans) :: rmass_ocean_load

  ! arrays with the mesh in double precision
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ystore
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: zstore

  double precision :: RHO_OCEANS

  ! processor identification
  character(len=150) prname

  ! time scheme
  real(kind=CUSTOM_REAL) :: deltat,deltatover2

  ! Stacey conditions put back
  integer :: NSPEC2D_TOP,NSPEC2D_BOTTOM,NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX
  integer :: nspec_stacey

  double precision, dimension(NGLLX) :: wxgll
  double precision, dimension(NGLLY) :: wygll
  double precision, dimension(NGLLZ) :: wzgll

  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX) :: jacobian2D_xmin,jacobian2D_xmax
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX) :: jacobian2D_ymin,jacobian2D_ymax
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM) :: jacobian2D_bottom
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP) :: jacobian2D_top

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX) :: normal_xmin,normal_xmax
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX) :: normal_ymin,normal_ymax

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_stacey) :: rho_vp,rho_vs

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  real(kind=CUSTOM_REAL) :: tx,ty,tz,sn
  real(kind=CUSTOM_REAL) :: nx,ny,nz,vn

  integer, dimension(NSPEC2D_TOP) :: ibelm_top
  integer, dimension(NSPEC2D_BOTTOM) :: ibelm_bottom
  integer, dimension(NSPEC2DMAX_XMIN_XMAX) :: ibelm_xmin,ibelm_xmax
  integer, dimension(NSPEC2DMAX_YMIN_YMAX) :: ibelm_ymin,ibelm_ymax

  integer, dimension(2,NSPEC2DMAX_YMIN_YMAX) :: nimin,nimax,nkmin_eta
  integer, dimension(2,NSPEC2DMAX_XMIN_XMAX) :: njmin,njmax,nkmin_xi

  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,ispec2D

  ! local parameters
  double precision :: xval,yval,zval,rval,thetaval,phival,weight
  double precision :: lat,lon,colat
  double precision :: elevation,height_oceans
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl

  integer :: ispec,i,j,k,iglob,ier
  integer :: ix_oceans,iy_oceans,iz_oceans,ispec_oceans,ispec2D_top_crust

  ! initializes matrices
  !
  ! in the case of stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
  ! on Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be obsolete

  rmassx(:) = 0._CUSTOM_REAL
  rmassy(:) = 0._CUSTOM_REAL
  rmassz(:) = 0._CUSTOM_REAL

  ! use the non-dimensional time step to make the mass matrix correction
  if(CUSTOM_REAL == SIZE_REAL) then
    deltat = sngl(DT*dsqrt(PI*GRAV*RHOAV))
  else
    deltat = DT*dsqrt(PI*GRAV*RHOAV)
  endif

  ! distinguish between single and double precision for reals
  if(CUSTOM_REAL == SIZE_REAL) then
     do i=1,NGLLX
        do j=1,NGLLY
           wgllwgll_xy(i,j) = sngl(wxgll(i)*wygll(j))
        enddo
     enddo

     do i=1,NGLLX
        do k=1,NGLLZ
           wgllwgll_xz(i,k) = sngl(wxgll(i)*wzgll(k))
        enddo
     enddo

     do j=1,NGLLY
        do k=1,NGLLZ
           wgllwgll_yz(j,k) = sngl(wygll(j)*wzgll(k))
        enddo
     enddo

     deltatover2 = sngl(0.5d0*deltat)

  else  ! double precision version
     do i=1,NGLLX
        do j=1,NGLLY
           wgllwgll_xy(i,j) = wxgll(i)*wygll(j)
        enddo
     enddo

     do i=1,NGLLX
        do k=1,NGLLZ
           wgllwgll_xz(i,k) = wxgll(i)*wzgll(k)
        enddo
     enddo

     do j=1,NGLLY
        do k=1,NGLLZ
           wgllwgll_yz(j,k) = wygll(j)*wzgll(k)
        enddo
     enddo

     deltatover2 = 0.5d0*deltat

  endif

  do ispec=1,nspec

    ! suppress fictitious elements in central cube
    if(idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          weight = wxgll(i)*wygll(j)*wzgll(k)
          iglob = ibool(i,j,k,ispec)

          ! compute the jacobian
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
          select case( iregion_code)

          case( IREGION_CRUST_MANTLE, IREGION_INNER_CORE )
            ! distinguish between single and double precision for reals
            if(CUSTOM_REAL == SIZE_REAL) then
              rmassz(iglob) = rmassz(iglob) + &
                     sngl(dble(rhostore(i,j,k,ispec)) * dble(jacobianl) * weight)
            else
              rmassz(iglob) = rmassz(iglob) + rhostore(i,j,k,ispec) * jacobianl * weight
            endif

          ! fluid in outer core
          case( IREGION_OUTER_CORE )

            ! no anisotropy in the fluid, use kappav

            ! distinguish between single and double precision for reals
            if(CUSTOM_REAL == SIZE_REAL) then
              rmassz(iglob) = rmassz(iglob) + &
                     sngl(dble(jacobianl) * weight * dble(rhostore(i,j,k,ispec)) / dble(kappavstore(i,j,k,ispec)))
            else
              rmassz(iglob) = rmassz(iglob) + &
                     jacobianl * weight * rhostore(i,j,k,ispec) / kappavstore(i,j,k,ispec)
            endif

          case default
            call exit_MPI(myrank,'wrong region code')

          end select

        enddo
      enddo
    enddo
  enddo

  ! save ocean load mass matrix as well if oceans
  if(OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) then

    ! create ocean load mass matrix for degrees of freedom at ocean bottom
    rmass_ocean_load(:) = 0._CUSTOM_REAL

    ! add contribution of the oceans
    ! for surface elements exactly at the top of the crust (ocean bottom)
    do ispec2D_top_crust = 1,NSPEC2D_TOP

      ispec_oceans = ibelm_top(ispec2D_top_crust)

      iz_oceans = NGLLZ

      do ix_oceans = 1,NGLLX
        do iy_oceans = 1,NGLLY

          iglob=ibool(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)

          ! if 3D Earth, compute local height of oceans
          if(CASE_3D) then

            ! get coordinates of current point
            xval = xstore(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)
            yval = ystore(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)
            zval = zstore(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)

            ! map to latitude and longitude for bathymetry routine
            call xyz_2_rthetaphi_dble(xval,yval,zval,rval,thetaval,phival)
            call reduce(thetaval,phival)

            ! convert the geocentric colatitude to a geographic colatitude
            colat = PI/2.0d0 - datan(1.006760466d0*dcos(thetaval)/dmax1(TINYVAL,dsin(thetaval)))

            ! get geographic latitude and longitude in degrees
            lat = 90.0d0 - colat*180.0d0/PI
            lon = phival*180.0d0/PI

            ! compute elevation at current point
            elevation = 0.d0
            if( TOPOGRAPHY ) then
              call get_topo_bathy(lat,lon,elevation,ibathy_topo)
            endif

            ! non-dimensionalize the elevation, which is in meters
            ! and suppress positive elevation, which means no oceans
            if(elevation >= - MINIMUM_THICKNESS_3D_OCEANS) then
              height_oceans = 0.d0
            else
              height_oceans = dabs(elevation) / R_EARTH
            endif

          else
            ! if 1D Earth, use oceans of constant thickness everywhere
            height_oceans = THICKNESS_OCEANS_PREM
          endif

          ! take into account inertia of water column
          weight = wxgll(ix_oceans) * wygll(iy_oceans) &
                    * dble(jacobian2D_top(ix_oceans,iy_oceans,ispec2D_top_crust)) &
                    * dble(RHO_OCEANS) * height_oceans

          ! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            rmass_ocean_load(iglob) = rmass_ocean_load(iglob) + sngl(weight)
          else
            rmass_ocean_load(iglob) = rmass_ocean_load(iglob) + weight
          endif

        enddo
      enddo

    enddo

    ! add regular mass matrix to ocean load contribution
    rmass_ocean_load(:) = rmass_ocean_load(:) + rmassz(:)

  endif

  ! add C*deltat/2 contribution to the mass matrices on Stacey edges
  if(NCHUNKS /= 6 .and. ABSORBING_CONDITIONS) then

     ! read arrays for Stacey conditions
     open(unit=27,file=prname(1:len_trim(prname))//'stacey.bin', &
          status='old',form='unformatted',action='read',iostat=ier)
     if( ier /= 0 ) call exit_mpi(myrank,'error opening stacey.bin in create_mass_matrices')
     read(27) nimin
     read(27) nimax
     read(27) njmin
     read(27) njmax
     read(27) nkmin_xi
     read(27) nkmin_eta
     close(27)

     select case(iregion_code)

     case(IREGION_CRUST_MANTLE)

        rmassx(:) = rmassz(:)
        rmassy(:) = rmassz(:)

        !   xmin
        ! if two chunks exclude this face for one of them
        if(NCHUNKS == 1 .or. ichunk == CHUNK_AC) then

           do ispec2D=1,nspec2D_xmin

              ispec=ibelm_xmin(ispec2D)

              ! exclude elements that are not on absorbing edges
              if(nkmin_xi(1,ispec2D) == 0 .or. njmin(1,ispec2D) == 0) cycle

              i=1
              do k=nkmin_xi(1,ispec2D),NGLLZ
                 do j=njmin(1,ispec2D),njmax(1,ispec2D)
                    iglob=ibool(i,j,k,ispec)

                    nx = normal_xmin(1,j,k,ispec2D)
                    ny = normal_xmin(2,j,k,ispec2D)
                    nz = normal_xmin(3,j,k,ispec2D)

                    vn = deltatover2*(nx+ny+nz)

                    tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(deltatover2-vn*nx)
                    ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(deltatover2-vn*ny)
                    tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(deltatover2-vn*nz)

                    weight = jacobian2D_xmin(j,k,ispec2D)*wgllwgll_yz(j,k)

                    if(CUSTOM_REAL == SIZE_REAL) then
                       rmassx(iglob) = rmassx(iglob) + sngl(tx*weight)
                       rmassy(iglob) = rmassy(iglob) + sngl(ty*weight)
                       rmassz(iglob) = rmassz(iglob) + sngl(tz*weight)
                    else
                       rmassx(iglob) = rmassx(iglob) + tx*weight
                       rmassy(iglob) = rmassy(iglob) + ty*weight
                       rmassz(iglob) = rmassz(iglob) + tz*weight
                    endif
                 enddo
              enddo
           enddo

        endif ! NCHUNKS == 1 .or. ichunk == CHUNK_AC

        !   xmax
        ! if two chunks exclude this face for one of them
        if(NCHUNKS == 1 .or. ichunk == CHUNK_AB) then

           do ispec2D=1,nspec2D_xmax

              ispec=ibelm_xmax(ispec2D)

              ! exclude elements that are not on absorbing edges
              if(nkmin_xi(2,ispec2D) == 0 .or. njmin(2,ispec2D) == 0) cycle

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

                    if(CUSTOM_REAL == SIZE_REAL) then
                       rmassx(iglob) = rmassx(iglob) + sngl(tx*weight)
                       rmassy(iglob) = rmassy(iglob) + sngl(ty*weight)
                       rmassz(iglob) = rmassz(iglob) + sngl(tz*weight)
                    else
                       rmassx(iglob) = rmassx(iglob) + tx*weight
                       rmassy(iglob) = rmassy(iglob) + ty*weight
                       rmassz(iglob) = rmassz(iglob) + tz*weight
                    endif
                 enddo
              enddo
           enddo

        endif ! NCHUNKS == 1 .or. ichunk == CHUNK_AB

        !   ymin
        do ispec2D=1,nspec2D_ymin

           ispec=ibelm_ymin(ispec2D)

           ! exclude elements that are not on absorbing edges
           if(nkmin_eta(1,ispec2D) == 0 .or. nimin(1,ispec2D) == 0) cycle

           j=1
           do k=nkmin_eta(1,ispec2D),NGLLZ
              do i=nimin(1,ispec2D),nimax(1,ispec2D)
                iglob=ibool(i,j,k,ispec)

                 nx = normal_ymin(1,i,k,ispec2D)
                 ny = normal_ymin(2,i,k,ispec2D)
                 nz = normal_ymin(3,i,k,ispec2D)

                 vn = deltatover2*(nx+ny+nz)

                 tx = rho_vp(i,j,k,ispec)*vn*nx + rho_vs(i,j,k,ispec)*(deltatover2-vn*nx)
                 ty = rho_vp(i,j,k,ispec)*vn*ny + rho_vs(i,j,k,ispec)*(deltatover2-vn*ny)
                 tz = rho_vp(i,j,k,ispec)*vn*nz + rho_vs(i,j,k,ispec)*(deltatover2-vn*nz)

                 weight = jacobian2D_ymin(i,k,ispec2D)*wgllwgll_xz(i,k)

                 if(CUSTOM_REAL == SIZE_REAL) then
                    rmassx(iglob) = rmassx(iglob) + sngl(tx*weight)
                    rmassy(iglob) = rmassy(iglob) + sngl(ty*weight)
                    rmassz(iglob) = rmassz(iglob) + sngl(tz*weight)
                 else
                    rmassx(iglob) = rmassx(iglob) + tx*weight
                    rmassy(iglob) = rmassy(iglob) + ty*weight
                    rmassz(iglob) = rmassz(iglob) + tz*weight
                 endif
              enddo
           enddo
        enddo

        !   ymax
        do ispec2D=1,nspec2D_ymax

           ispec=ibelm_ymax(ispec2D)

           ! exclude elements that are not on absorbing edges
           if(nkmin_eta(2,ispec2D) == 0 .or. nimin(2,ispec2D) == 0) cycle

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

                 if(CUSTOM_REAL == SIZE_REAL) then
                    rmassx(iglob) = rmassx(iglob) + sngl(tx*weight)
                    rmassy(iglob) = rmassy(iglob) + sngl(ty*weight)
                    rmassz(iglob) = rmassz(iglob) + sngl(tz*weight)
                 else
                    rmassx(iglob) = rmassx(iglob) + tx*weight
                    rmassy(iglob) = rmassy(iglob) + ty*weight
                    rmassz(iglob) = rmassz(iglob) + tz*weight
                 endif
              enddo
           enddo
        enddo

        ! check that mass matrix is positive
        if(minval(rmassx(:)) <= 0.) call exit_MPI(myrank,'negative rmassx matrix term')
        if(minval(rmassy(:)) <= 0.) call exit_MPI(myrank,'negative rmassy matrix term')

     case(IREGION_OUTER_CORE)

        !   xmin
        ! if two chunks exclude this face for one of them
        if(NCHUNKS == 1 .or. ichunk == CHUNK_AC) then

           do ispec2D=1,nspec2D_xmin

              ispec=ibelm_xmin(ispec2D)

              ! exclude elements that are not on absorbing edges
              if(nkmin_xi(1,ispec2D) == 0 .or. njmin(1,ispec2D) == 0) cycle

              i=1
              do k=nkmin_xi(1,ispec2D),NGLLZ
                 do j=njmin(1,ispec2D),njmax(1,ispec2D)
                    iglob=ibool(i,j,k,ispec)

                    sn = deltatover2/rho_vp(i,j,k,ispec)

                    weight = jacobian2D_xmin(j,k,ispec2D)*wgllwgll_yz(j,k)

                    if(CUSTOM_REAL == SIZE_REAL) then
                       rmassz(iglob) = rmassz(iglob) + sngl(weight*sn)
                    else
                       rmassz(iglob) = rmassz(iglob) + weight*sn
                    endif
                 enddo
              enddo
           enddo

        endif ! NCHUNKS == 1 .or. ichunk == CHUNK_AC

        !   xmax
        ! if two chunks exclude this face for one of them
        if(NCHUNKS == 1 .or. ichunk == CHUNK_AB) then

           do ispec2D=1,nspec2D_xmax

              ispec=ibelm_xmax(ispec2D)

              ! exclude elements that are not on absorbing edges
              if(nkmin_xi(2,ispec2D) == 0 .or. njmin(2,ispec2D) == 0) cycle

              i=NGLLX
              do k=nkmin_xi(2,ispec2D),NGLLZ
                 do j=njmin(2,ispec2D),njmax(2,ispec2D)
                    iglob=ibool(i,j,k,ispec)

                    sn = deltatover2/rho_vp(i,j,k,ispec)

                    weight = jacobian2D_xmax(j,k,ispec2D)*wgllwgll_yz(j,k)

                    if(CUSTOM_REAL == SIZE_REAL) then
                       rmassz(iglob) = rmassz(iglob) + sngl(weight*sn)
                    else
                       rmassz(iglob) = rmassz(iglob) + weight*sn
                    endif
                 enddo
              enddo
           enddo

        endif ! NCHUNKS == 1 .or. ichunk == CHUNK_AB

        !   ymin
        do ispec2D=1,nspec2D_ymin

           ispec=ibelm_ymin(ispec2D)

           ! exclude elements that are not on absorbing edges
           if(nkmin_eta(1,ispec2D) == 0 .or. nimin(1,ispec2D) == 0) cycle

           j=1
           do k=nkmin_eta(1,ispec2D),NGLLZ
              do i=nimin(1,ispec2D),nimax(1,ispec2D)
                 iglob=ibool(i,j,k,ispec)

                 sn = deltatover2/rho_vp(i,j,k,ispec)

                 weight = jacobian2D_ymin(i,k,ispec2D)*wgllwgll_xz(i,k)

                 if(CUSTOM_REAL == SIZE_REAL) then
                    rmassz(iglob) = rmassz(iglob) + sngl(weight*sn)
                 else
                    rmassz(iglob) = rmassz(iglob) + weight*sn
                 endif
              enddo
           enddo
        enddo

        !   ymax
        do ispec2D=1,nspec2D_ymax

           ispec=ibelm_ymax(ispec2D)

           ! exclude elements that are not on absorbing edges
           if(nkmin_eta(2,ispec2D) == 0 .or. nimin(2,ispec2D) == 0) cycle

           j=NGLLY
           do k=nkmin_eta(2,ispec2D),NGLLZ
              do i=nimin(2,ispec2D),nimax(2,ispec2D)
                 iglob=ibool(i,j,k,ispec)

                 sn = deltatover2/rho_vp(i,j,k,ispec)

                 weight = jacobian2D_ymax(i,k,ispec2D)*wgllwgll_xz(i,k)

                 if(CUSTOM_REAL == SIZE_REAL) then
                    rmassz(iglob) = rmassz(iglob) + sngl(weight*sn)
                 else
                    rmassz(iglob) = rmassz(iglob) + weight*sn
                 endif
              enddo
           enddo
        enddo

        !   bottom (zmin)
        do ispec2D=1,NSPEC2D_BOTTOM

           ispec=ibelm_bottom(ispec2D)

           k=1
           do j=1,NGLLY
              do i=1,NGLLX
                 iglob=ibool(i,j,k,ispec)

                 sn = deltatover2/rho_vp(i,j,k,ispec)

                 weight = jacobian2D_bottom(i,j,ispec2D)*wgllwgll_xy(i,j)

                 if(CUSTOM_REAL == SIZE_REAL) then
                    rmassz(iglob) = rmassz(iglob) + sngl(weight*sn)
                 else
                    rmassz(iglob) = rmassz(iglob) + weight*sn
                 endif
              enddo
           enddo
        enddo

     case( IREGION_INNER_CORE )

     case default
        call exit_MPI(myrank,'wrong region code')

     end select

  endif

  ! check that mass matrix is positive
  ! note: in fictitious elements it is still zero
  if(minval(rmassz(:)) < 0._CUSTOM_REAL) call exit_MPI(myrank,'negative rmassz matrix term')

  end subroutine create_mass_matrices
