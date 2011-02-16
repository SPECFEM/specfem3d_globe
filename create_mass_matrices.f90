!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2011
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
                          iregion_code,nglob,rmass,rhostore,kappavstore, &
                          nglob_oceans,rmass_ocean_load,NSPEC2D_TOP,ibelm_top,jacobian2D_top, &
                          xstore,ystore,zstore,RHO_OCEANS)

! creates rmass and rmass_ocean_load

  use meshfem3D_models_par

  implicit none

  integer myrank,nspec

  integer idoubling(nspec)

  double precision wxgll(NGLLX),wygll(NGLLY),wzgll(NGLLZ)

  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

  integer nspec_actually
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_actually) ::  &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

  integer iregion_code

  ! mass matrix
  integer nglob
  real(kind=CUSTOM_REAL), dimension(nglob) :: rmass

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: rhostore,kappavstore

  ! ocean mass matrix
  integer nglob_oceans
  real(kind=CUSTOM_REAL), dimension(nglob_oceans) :: rmass_ocean_load

  integer NSPEC2D_TOP
  integer, dimension(NSPEC2D_TOP) :: ibelm_top
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP) :: jacobian2D_top

  ! arrays with the mesh in double precision
  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

  double precision RHO_OCEANS

  ! local parameters
  double precision weight
  double precision xval,yval,zval,rval,thetaval,phival
  double precision lat,lon,colat
  double precision elevation,height_oceans
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl

  integer :: ispec,i,j,k,iglobnum
  integer :: ix_oceans,iy_oceans,iz_oceans,ispec_oceans,ispec2D_top_crust


  ! initializes
  rmass(:) = 0._CUSTOM_REAL

  do ispec=1,nspec

    ! suppress fictitious elements in central cube
    if(idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          weight = wxgll(i)*wygll(j)*wzgll(k)
          iglobnum = ibool(i,j,k,ispec)

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
          if(iregion_code == IREGION_CRUST_MANTLE .or. iregion_code == IREGION_INNER_CORE) then

            ! distinguish between single and double precision for reals
            if(CUSTOM_REAL == SIZE_REAL) then
              rmass(iglobnum) = rmass(iglobnum) + &
                     sngl(dble(rhostore(i,j,k,ispec)) * dble(jacobianl) * weight)
            else
              rmass(iglobnum) = rmass(iglobnum) + rhostore(i,j,k,ispec) * jacobianl * weight
            endif

          ! fluid in outer core
          else if(iregion_code == IREGION_OUTER_CORE) then

            ! no anisotropy in the fluid, use kappav

            ! distinguish between single and double precision for reals
            if(CUSTOM_REAL == SIZE_REAL) then
              rmass(iglobnum) = rmass(iglobnum) + &
                     sngl(dble(jacobianl) * weight * dble(rhostore(i,j,k,ispec)) / dble(kappavstore(i,j,k,ispec)))
            else
              rmass(iglobnum) = rmass(iglobnum) + &
                     jacobianl * weight * rhostore(i,j,k,ispec) / kappavstore(i,j,k,ispec)
            endif

          else
            call exit_MPI(myrank,'wrong region code')
          endif

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

          iglobnum=ibool(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)

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
            elevation = 0.d0

            ! compute elevation at current point
            call get_topo_bathy(lat,lon,elevation,ibathy_topo)

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
          weight = wxgll(ix_oceans)*wygll(iy_oceans)*dble(jacobian2D_top(ix_oceans,iy_oceans,ispec2D_top_crust)) &
                 * dble(RHO_OCEANS) * height_oceans

          ! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            rmass_ocean_load(iglobnum) = rmass_ocean_load(iglobnum) + sngl(weight)
          else
            rmass_ocean_load(iglobnum) = rmass_ocean_load(iglobnum) + weight
          endif

        enddo
      enddo

    enddo

    ! add regular mass matrix to ocean load contribution
    rmass_ocean_load(:) = rmass_ocean_load(:) + rmass(:)

  endif

  end subroutine create_mass_matrices
