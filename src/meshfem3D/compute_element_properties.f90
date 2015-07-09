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

! compute several rheological and geometrical properties for a given spectral element
  subroutine compute_element_properties(ispec,iregion_code,idoubling,ipass, &
                         xstore,ystore,zstore,nspec,myrank, &
                         xelm,yelm,zelm,shape3D, &
                         rmin,rmax, &
                         rhostore,dvpstore, &
                         kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                         xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                         gammaxstore,gammaystore,gammazstore,nspec_actually, &
                         c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                         c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                         c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                         nspec_ani,nspec_stacey, &
                         rho_vp,rho_vs,&
                         xigll,yigll,zigll,ispec_is_tiso)

  use meshfem3D_models_par

  implicit none

  ! correct number of spectral elements in each block depending on chunk type
  integer ispec,nspec,nspec_stacey

! arrays with the mesh in double precision
  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

! code for the four regions of the mesh
  integer :: iregion_code

! meshing phase
  integer :: ipass

! 3D shape functions and their derivatives
  double precision, dimension(NGNOD,NGLLX,NGLLY,NGLLZ) :: shape3D

  double precision, dimension(NGNOD) :: xelm,yelm,zelm

! parameters needed to store the radii of the grid points
! in the spherically symmetric Earth
  integer,dimension(nspec) :: idoubling
  double precision :: rmin,rmax

! for model density and anisotropy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: rhostore,dvpstore,kappavstore, &
    kappahstore,muvstore,muhstore,eta_anisostore

! the 21 coefficients for an anisotropic medium in reduced notation
  integer :: nspec_ani
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_ani) :: &
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store

! arrays with mesh parameters
  integer :: nspec_actually
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_actually) :: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

! proc numbers for MPI
  integer :: myrank

! Stacey, indices for Clayton-Engquist absorbing conditions
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_stacey) :: rho_vp,rho_vs

  ! Parameters used to calculate Jacobian based upon 125 GLL points
  double precision :: xigll(NGLLX)
  double precision :: yigll(NGLLY)
  double precision :: zigll(NGLLZ)

  logical, dimension(nspec) :: ispec_is_tiso

  ! Parameter used to decide whether this element is in the crust or not
  logical:: elem_in_crust,elem_in_mantle
  ! flag for transverse isotropic elements
  logical:: elem_is_tiso

  ! note: at this point, the mesh is still perfectly spherical

  ! add topography of the Moho *before* adding the 3D crustal velocity model so that the stretched
  ! mesh gets assigned the right model values
  elem_in_crust = .false.
  elem_in_mantle = .false.
  elem_is_tiso = .false.
  if (iregion_code == IREGION_CRUST_MANTLE) then
    if (CRUSTAL .and. CASE_3D) then
      ! 3D crustal models
      if (idoubling(ispec) == IFLAG_CRUST &
        .or. idoubling(ispec) == IFLAG_220_80 &
        .or. idoubling(ispec) == IFLAG_80_MOHO) then
        ! Stretch mesh to honor smoothed moho thickness from crust2.0
        ! mesh is stretched between surface and 220
        !
        ! differentiate between regional and global meshing
        if (REGIONAL_MOHO_MESH) then
          call moho_stretching_honor_crust_reg(myrank,xelm,yelm,zelm, &
                                               elem_in_crust,elem_in_mantle)
        else
          call moho_stretching_honor_crust(myrank,xelm,yelm,zelm, &
                                           elem_in_crust,elem_in_mantle)
        endif
      else
        ! element below 220km
        ! sets element flag for mantle
        elem_in_mantle = .true.
      endif
    else
      ! 1D crust, no stretching
      ! sets element flags
      if (idoubling(ispec) == IFLAG_CRUST) then
        elem_in_crust = .true.
      else
        elem_in_mantle = .true.
      endif
    endif

    ! sets transverse isotropic flag for elements in mantle
    if (TRANSVERSE_ISOTROPY) then
      ! modifies tiso to have it for all mantle elements
      ! preferred for example, when using 1Dref (STW model)
      if (USE_FULL_TISO_MANTLE) then
        ! all elements below the actual moho will be used for transverse isotropy
        ! note: this will increase the computation time by ~ 45 %
        if (elem_in_mantle) then
          elem_is_tiso = .true.
        endif
      else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1DREF) then
        ! transverse isotropic mantle between fictitious moho to 670km depth
        ! preferred for Harvard (Kustowski's) models using STW 1D reference, i.e.
        ! THREE_D_MODEL_S362ANI
        ! THREE_D_MODEL_S362WMANI
        ! THREE_D_MODEL_S29EA
        ! THREE_D_MODEL_GLL
        ! which show significant transverse isotropy also below 220km depth
        if (USE_OLD_VERSION_5_1_5_FORMAT) then
          ! assigns TI only to elements below (2-layer) fictitious moho down to 670
          if (idoubling(ispec)==IFLAG_220_80 &
            .or. idoubling(ispec)==IFLAG_80_MOHO &
            .or. idoubling(ispec)==IFLAG_670_220) then
            elem_is_tiso = .true.
          endif
        else
          ! assigns TI to elements in mantle elements just below actual moho down to 670
          if (idoubling(ispec)==IFLAG_220_80 &
            .or. idoubling(ispec)==IFLAG_80_MOHO &
            .or. idoubling(ispec)==IFLAG_670_220 &
            .or. (idoubling(ispec)==IFLAG_CRUST .and. elem_in_mantle )) then
            elem_is_tiso = .true.
          endif
        endif
      else if (idoubling(ispec)==IFLAG_220_80 .or. idoubling(ispec)==IFLAG_80_MOHO) then
        ! default case for PREM reference models:
        ! models use only transverse isotropy between moho and 220 km depth
        elem_is_tiso = .true.
        ! checks mantle flag to be sure
        !if (elem_in_mantle .eqv. .false. ) stop 'Error mantle flag confused between moho and 220'
      endif
    endif

  endif ! IREGION_CRUST_MANTLE

  ! sets element tiso flag
  ispec_is_tiso(ispec) = elem_is_tiso

  ! interpolates and stores GLL point locations
  call compute_element_GLL_locations(xelm,yelm,zelm,ispec,nspec, &
                                     xstore,ystore,zstore,shape3D)

  ! computes velocity/density/... values for the chosen Earth model
  ! (only needed for second meshing phase)
  if (ipass == 2) then
    call get_model(myrank,iregion_code,ispec,nspec,idoubling(ispec), &
                   kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                   rhostore,dvpstore,nspec_ani, &
                   c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                   c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                   c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                   nspec_stacey,rho_vp,rho_vs, &
                   xstore,ystore,zstore, &
                   rmin,rmax, &
                   elem_in_crust,elem_in_mantle)

  endif

  ! either use GLL points or anchor points to capture TOPOGRAPHY and ELLIPTICITY
  !
  ! note:  using GLL points to capture them results in a slightly more accurate mesh.
  !           however, it introduces more deformations to the elements which might lead to
  !           problems with the Jacobian. using the anchors is therefore more robust.

  ! adds surface topography
  if (TOPOGRAPHY) then
    if (idoubling(ispec) == IFLAG_CRUST .or. &
        idoubling(ispec) == IFLAG_220_80 .or. &
        idoubling(ispec) == IFLAG_80_MOHO) then
      ! stretches mesh between surface and R220 accordingly
      if (USE_GLL) then
        ! stretches every GLL point accordingly
        call add_topography_gll(myrank,xstore,ystore,zstore,ispec,nspec,ibathy_topo)
      else
        ! stretches anchor points only, interpolates GLL points later on
        call add_topography(myrank,xelm,yelm,zelm,ibathy_topo)
      endif
    endif
  endif

  ! adds topography on 410 km and 650 km discontinuity in model S362ANI
  if (.not. SUPPRESS_INTERNAL_TOPOGRAPHY) then
    if (THREE_D_MODEL == THREE_D_MODEL_S362ANI .or. THREE_D_MODEL == THREE_D_MODEL_S362WMANI &
      .or. THREE_D_MODEL == THREE_D_MODEL_S362ANI_PREM .or. THREE_D_MODEL == THREE_D_MODEL_S29EA) then
      ! stretching between 220 and 770
      if (idoubling(ispec) == IFLAG_670_220 .or. &
          idoubling(ispec) == IFLAG_MANTLE_NORMAL) then
        if (USE_GLL) then
          ! stretches every GLL point accordingly
          call add_topography_410_650_gll(myrank,xstore,ystore,zstore,ispec,nspec)
        else
          ! stretches anchor points only, interpolates GLL points later on
          call add_topography_410_650(myrank,xelm,yelm,zelm)
        endif
      endif
    endif

    ! these are placeholders:
    ! their corresponding subroutines subtopo_cmb() and subtopo_icb() are not implemented yet....
    ! must be done/supplied by the user; uncomment in case
    ! CMB topography
    !  if (THREE_D_MODEL == THREE_D_MODEL_S362ANI .and. (idoubling(ispec)==IFLAG_MANTLE_NORMAL &
    !     .or. idoubling(ispec)==IFLAG_OUTER_CORE_NORMAL)) &
    !           call add_topography_cmb(myrank,xelm,yelm,zelm,RTOPDDOUBLEPRIME,RCMB)

    ! ICB topography
    !  if (THREE_D_MODEL == THREE_D_MODEL_S362ANI .and. (idoubling(ispec)==IFLAG_OUTER_CORE_NORMAL &
    !     .or. idoubling(ispec)==IFLAG_INNER_CORE_NORMAL .or. idoubling(ispec)==IFLAG_MIDDLE_CENTRAL_CUBE &
    !     .or. idoubling(ispec)==IFLAG_BOTTOM_CENTRAL_CUBE .or. idoubling(ispec)==IFLAG_TOP_CENTRAL_CUBE &
    !     .or. idoubling(ispec)==IFLAG_IN_FICTITIOUS_CUBE)) &
    !           call add_topography_icb(myrank,xelm,yelm,zelm,RICB,RCMB)
  endif


  ! make the Earth elliptical
  if (ELLIPTICITY) then
    ! note: after adding ellipticity, the mesh becomes elliptical and geocentric and geodetic/geographic colatitudes differ.
    if (USE_GLL) then
      ! make the Earth's ellipticity, use GLL points
      call get_ellipticity_gll(xstore,ystore,zstore,ispec,nspec,nspl,rspl,espl,espl2)
    else
      ! make the Earth's ellipticity, use element anchor points
      call get_ellipticity(xelm,yelm,zelm,nspl,rspl,espl,espl2)
    endif
  endif

  ! re-interpolates and creates the GLL point locations since the anchor points might have moved
  !
  ! note: velocity values associated for each GLL point will "move" along together with
  !          their associated points. however, we don't re-calculate the velocity model values since the
  !          models are/should be referenced with respect to a spherical Earth.
  if (.not. USE_GLL) then
    call compute_element_GLL_locations(xelm,yelm,zelm,ispec,nspec, &
                                      xstore,ystore,zstore,shape3D)
  endif

  ! updates Jacobian
  ! (only needed for second meshing phase)
  if (ipass == 2) then
    call recalc_jacobian_gll3D(myrank,xstore,ystore,zstore,xigll,yigll,zigll,&
                                ispec,nspec,&
                                xixstore,xiystore,xizstore,&
                                etaxstore,etaystore,etazstore,&
                                gammaxstore,gammaystore,gammazstore)
  endif

  end subroutine compute_element_properties

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_element_GLL_locations(xelm,yelm,zelm,ispec,nspec, &
                                      xstore,ystore,zstore,shape3D)

  use constants

  implicit none

  integer :: ispec,nspec

  double precision,dimension(NGNOD) :: xelm,yelm,zelm

  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  double precision,dimension(NGNOD,NGLLX,NGLLY,NGLLZ) :: shape3D

  ! local parameters
  double precision :: xmesh,ymesh,zmesh
  integer :: i,j,k,ia

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        xmesh = ZERO
        ymesh = ZERO
        zmesh = ZERO

        ! interpolates the location using 3D shape functions
        do ia = 1,NGNOD

          xmesh = xmesh + shape3D(ia,i,j,k)*xelm(ia)
          ymesh = ymesh + shape3D(ia,i,j,k)*yelm(ia)
          zmesh = zmesh + shape3D(ia,i,j,k)*zelm(ia)

        enddo

        ! stores mesh coordinates
        xstore(i,j,k,ispec) = xmesh
        ystore(i,j,k,ispec) = ymesh
        zstore(i,j,k,ispec) = zmesh

      enddo
    enddo
  enddo

  end subroutine compute_element_GLL_locations

