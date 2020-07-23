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

! compute several rheological and geometrical properties for a given spectral element
  subroutine compute_element_properties(ispec,iregion_code,idoubling,ipass, &
                                        xstore,ystore,zstore,nspec, &
                                        xelm,yelm,zelm,shape3D, &
                                        rmin,rmax, &
                                        xigll,yigll,zigll,ispec_is_tiso)

  use constants, only: myrank,NGLLX,NGLLY,NGLLZ,NGNOD,CUSTOM_REAL, &
    IFLAG_220_80,IFLAG_670_220,IFLAG_80_MOHO,IFLAG_MANTLE_NORMAL,IFLAG_CRUST, &
    IFLAG_OUTER_CORE_NORMAL,IFLAG_IN_FICTITIOUS_CUBE, &
    IREGION_CRUST_MANTLE,SUPPRESS_INTERNAL_TOPOGRAPHY,USE_GLL

  use meshfem3D_models_par, only: &
    TOPOGRAPHY,ELLIPTICITY,CRUSTAL,CASE_3D, &
    THREE_D_MODEL,THREE_D_MODEL_MANTLE_SH,THREE_D_MODEL_S29EA, &
    THREE_D_MODEL_S362ANI,THREE_D_MODEL_S362WMANI,THREE_D_MODEL_S362ANI_PREM, &
    ibathy_topo,nspl,rspl,ellipicity_spline,ellipicity_spline2, &
    REGIONAL_MOHO_MESH


  use regions_mesh_par2, only: &
    xixstore,xiystore,xizstore, &
    etaxstore,etaystore,etazstore, &
    gammaxstore,gammaystore,gammazstore

  implicit none

  ! correct number of spectral elements in each block depending on chunk type
  integer,intent(in) :: ispec,nspec

! arrays with the mesh in double precision
  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout) :: xstore,ystore,zstore

! code for the four regions of the mesh
  integer,intent(in) :: iregion_code

! meshing phase
  integer,intent(in) :: ipass

! 3D shape functions and their derivatives
  double precision, dimension(NGNOD,NGLLX,NGLLY,NGLLZ),intent(in) :: shape3D

  double precision, dimension(NGNOD),intent(inout) :: xelm,yelm,zelm

! parameters needed to store the radii of the grid points
! in the spherically symmetric Earth
  integer,dimension(nspec),intent(in) :: idoubling
  double precision,intent(in) :: rmin,rmax

  ! Parameters used to calculate Jacobian based upon 125 GLL points
  double precision,intent(in) :: xigll(NGLLX)
  double precision,intent(in) :: yigll(NGLLY)
  double precision,intent(in) :: zigll(NGLLZ)

  logical, dimension(nspec),intent(out) :: ispec_is_tiso

  ! local parameters
  ! Parameter used to decide whether this element is in the crust or not
  logical :: elem_in_crust,elem_in_mantle
  ! flag for transverse isotropic elements
  logical :: elem_is_tiso

  !debug
  logical, parameter :: DEBUG_OUTPUT = .false.

! note: at this point, the mesh is still perfectly spherical

  ! flag if element completely in crust (all corners above moho)
  elem_in_crust = .false.
  ! flag if element completely in mantle (all corners below moho)
  elem_in_mantle = .false.

  !debug
  if (DEBUG_OUTPUT) then
    if (myrank == 0) print *,'element ',ispec,' properties:'
  endif

  ! add topography of the Moho *before* adding the 3D crustal velocity model so that the stretched
  ! mesh gets assigned the right model values
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
          call moho_stretching_honor_crust_reg(xelm,yelm,zelm, &
                                               elem_in_crust,elem_in_mantle)
        else
          call moho_stretching_honor_crust(xelm,yelm,zelm, &
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
    !debug
    if (DEBUG_OUTPUT) then
      if (myrank == 0) print *,'  in crust',elem_in_crust,' in mantle ',elem_in_mantle
    endif
  endif ! IREGION_CRUST_MANTLE

  ! sets element tiso flag
  call compute_element_tiso_flag(elem_is_tiso,elem_in_mantle,iregion_code,ispec,nspec,idoubling)

  ! stores as element flags
  ispec_is_tiso(ispec) = elem_is_tiso

  !debug
  if (DEBUG_OUTPUT) then
    if (myrank == 0) print *,'  tiso flag ',elem_is_tiso
  endif

  ! interpolates and stores GLL point locations
  call compute_element_GLL_locations(xelm,yelm,zelm,ispec,nspec, &
                                     xstore,ystore,zstore,shape3D)

  !debug
  if (DEBUG_OUTPUT) then
    if (myrank == 0) print *,'  locations done'
  endif

  ! computes velocity/density/... values for the chosen Earth model
  ! (only needed for second meshing phase)
  if (ipass == 2) then
    call get_model(iregion_code,ispec,nspec,idoubling(ispec), &
                   xstore,ystore,zstore, &
                   rmin,rmax, &
                   elem_in_crust,elem_in_mantle)

    !debug
    if (DEBUG_OUTPUT) then
      if (myrank == 0) print *,'  model properties done'
    endif
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
        call add_topography_gll(xstore,ystore,zstore,ispec,nspec,ibathy_topo)
      else
        ! stretches anchor points only, interpolates GLL points later on
        call add_topography(xelm,yelm,zelm,ibathy_topo)
      endif
    endif

    !debug
    if (DEBUG_OUTPUT) then
      if (myrank == 0) print *,'  topography done'
    endif
  endif


  ! adds topography on 410 km and 650 km discontinuity in model S362ANI
  if (.not. SUPPRESS_INTERNAL_TOPOGRAPHY) then
    ! s362ani internal topography
    if (THREE_D_MODEL == THREE_D_MODEL_S362ANI .or. THREE_D_MODEL == THREE_D_MODEL_S362WMANI &
      .or. THREE_D_MODEL == THREE_D_MODEL_S362ANI_PREM .or. THREE_D_MODEL == THREE_D_MODEL_S29EA) then
      ! stretching between 220 and 770
      if (idoubling(ispec) == IFLAG_670_220 .or. &
          idoubling(ispec) == IFLAG_MANTLE_NORMAL) then
        if (USE_GLL) then
          ! stretches every GLL point accordingly
          call add_topography_410_650_gll(xstore,ystore,zstore,ispec,nspec)
        else
          ! stretches anchor points only, interpolates GLL points later on
          call add_topography_410_650(xelm,yelm,zelm)
        endif
      endif
    endif

    ! full_sh model
    if (THREE_D_MODEL == THREE_D_MODEL_MANTLE_SH) then
      ! 410/650 topography
      ! stretching between 220 and 770
      if (idoubling(ispec) == IFLAG_670_220 .or. &
          idoubling(ispec) == IFLAG_MANTLE_NORMAL) then
        call add_topography_sh_mantle(xelm,yelm,zelm)
      endif

      ! CMB topography
      ! stretching lower mantle/outer core
      if (idoubling(ispec) == IFLAG_MANTLE_NORMAL .or. &
          idoubling(ispec) == IFLAG_OUTER_CORE_NORMAL) then
        call add_topography_sh_cmb(xelm,yelm,zelm)
      endif
    endif

    ! these are placeholders:
    ! their corresponding subroutines subtopo_cmb() and subtopo_icb() are not implemented yet....
    ! must be done/supplied by the user; uncomment in case
    ! CMB topography
    !  if (THREE_D_MODEL == THREE_D_MODEL_S362ANI .and. (idoubling(ispec)==IFLAG_MANTLE_NORMAL &
    ! .or. idoubling(ispec)==IFLAG_OUTER_CORE_NORMAL)) &
    !           call add_topography_cmb(xelm,yelm,zelm)

    ! ICB topography
    !  if (THREE_D_MODEL == THREE_D_MODEL_S362ANI .and. (idoubling(ispec)==IFLAG_OUTER_CORE_NORMAL &
    ! .or. idoubling(ispec)==IFLAG_INNER_CORE_NORMAL .or. idoubling(ispec)==IFLAG_MIDDLE_CENTRAL_CUBE &
    ! .or. idoubling(ispec)==IFLAG_BOTTOM_CENTRAL_CUBE .or. idoubling(ispec)==IFLAG_TOP_CENTRAL_CUBE &
    ! .or. idoubling(ispec)==IFLAG_IN_FICTITIOUS_CUBE)) &
    !           call add_topography_icb(xelm,yelm,zelm)

    !debug
    if (DEBUG_OUTPUT) then
      if (myrank == 0) print *,'  internal topography done'
    endif
  endif

  ! make the Earth elliptical
  if (ELLIPTICITY) then
    ! note: after adding ellipticity, the mesh becomes elliptical and geocentric and geodetic/geographic colatitudes differ.
    if (USE_GLL) then
      ! make the Earth's ellipticity, use GLL points
      call get_ellipticity_gll(xstore,ystore,zstore,ispec,nspec,nspl,rspl,ellipicity_spline,ellipicity_spline2)
    else
      ! make the Earth's ellipticity, use element anchor points
      call get_ellipticity(xelm,yelm,zelm,nspl,rspl,ellipicity_spline,ellipicity_spline2)
    endif

    !debug
    if (DEBUG_OUTPUT) then
      if (myrank == 0) print *,'  ellipticity done'
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

    !debug
    if (DEBUG_OUTPUT) then
      if (myrank == 0) print *,'  GLL locations done'
    endif
  endif

  ! updates Jacobian
  ! (only needed for second meshing phase)
  if (ipass == 2) then
    if (idoubling(ispec) /= IFLAG_IN_FICTITIOUS_CUBE) then
      call recalc_jacobian_gll3D(xstore,ystore,zstore,xigll,yigll,zigll, &
                                 ispec,nspec, &
                                 xixstore,xiystore,xizstore, &
                                 etaxstore,etaystore,etazstore, &
                                 gammaxstore,gammaystore,gammazstore)
    endif

    !debug
    if (DEBUG_OUTPUT) then
      if (myrank == 0) print *,'  jacobian done'
    endif
  endif

  end subroutine compute_element_properties

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_element_GLL_locations(xelm,yelm,zelm,ispec,nspec, &
                                           xstore,ystore,zstore,shape3D)

  use constants, only: NGLLX,NGLLY,NGLLZ,NGNOD,ZERO

  implicit none

  integer,intent(in) :: ispec,nspec

  double precision,dimension(NGNOD),intent(in) :: xelm,yelm,zelm

  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout) :: xstore,ystore,zstore

  double precision,dimension(NGNOD,NGLLX,NGLLY,NGLLZ),intent(in) :: shape3D

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


!
!-------------------------------------------------------------------------------------------------
!


  subroutine compute_element_tiso_flag(elem_is_tiso,elem_in_mantle,iregion_code,ispec,nspec,idoubling)

! sets transverse isotropic flag for elements in crust/mantle

  use constants, only: IMAIN,myrank,USE_OLD_VERSION_7_0_0_FORMAT,USE_OLD_VERSION_5_1_5_FORMAT, &
    IFLAG_CRUST,IFLAG_220_80,IFLAG_80_MOHO,IFLAG_670_220,IFLAG_MANTLE_NORMAL,IREGION_CRUST_MANTLE, &
    REFERENCE_MODEL_1DREF,REFERENCE_MODEL_1DREF, &
    THREE_D_MODEL_S362WMANI,THREE_D_MODEL_SGLOBE

  use meshfem3D_models_par, only: &
    TRANSVERSE_ISOTROPY,USE_FULL_TISO_MANTLE,REFERENCE_1D_MODEL,THREE_D_MODEL

  implicit none

  logical,intent(out) :: elem_is_tiso
  logical,intent(in) :: elem_in_mantle
  integer,intent(in) :: iregion_code,ispec
  integer,intent(in) :: nspec
  integer,dimension(nspec),intent(in) :: idoubling

  ! local parameters
  logical, save :: is_first_call = .true.

  ! initializes
  elem_is_tiso = .false.

  ! checks if anything to do
  if (.not. TRANSVERSE_ISOTROPY) return
  if (iregion_code /= IREGION_CRUST_MANTLE) return

  ! user output
  if (is_first_call) then
    is_first_call = .false.
    if (myrank == 0) then
      ! only output once
      write(IMAIN,*) '  setting tiso flags in mantle model'
      if (USE_FULL_TISO_MANTLE) &
        write(IMAIN,*) '    using fully transverse isotopic mantle'
      if (USE_OLD_VERSION_7_0_0_FORMAT) &
        write(IMAIN,*) '    using formatting from version 7.0.0'
      if (USE_OLD_VERSION_5_1_5_FORMAT) &
        write(IMAIN,*) '    using formatting from version 5.1.5'
      call flush_IMAIN()
    endif
  endif

  ! transverse isotropic models
  ! modifies tiso to have it for all mantle elements
  ! preferred for example, when using 1Dref (STW model)
  if (USE_FULL_TISO_MANTLE) then
    ! all elements below the actual moho will be used for transverse isotropy
    ! note: this will increase the computation time by ~ 45 %
    if (USE_OLD_VERSION_7_0_0_FORMAT) then
      if (elem_in_mantle) elem_is_tiso = .true.
    else
      if (idoubling(ispec) == IFLAG_MANTLE_NORMAL &
        .or. idoubling(ispec) == IFLAG_670_220 &
        .or. idoubling(ispec) == IFLAG_220_80 &
        .or. idoubling(ispec) == IFLAG_80_MOHO &
        .or. idoubling(ispec) == IFLAG_CRUST) then
        elem_is_tiso = .true.
      endif
    endif

    ! all done
    return
  endif

  ! TISO flags distributions

  ! reference model specifics
  select case (REFERENCE_1D_MODEL)
  case (REFERENCE_MODEL_1DREF)
    ! transverse isotropic mantle between fictitious moho to 670km depth
    ! preferred for Harvard (Kustowski's) models using STW 1D reference, i.e.
    ! THREE_D_MODEL_S362ANI
    ! THREE_D_MODEL_S362WMANI
    ! THREE_D_MODEL_S29EA
    ! THREE_D_MODEL_GLL
    ! which show significant transverse isotropy also below 220km depth
    if (USE_OLD_VERSION_5_1_5_FORMAT) then
      ! assigns TI only to elements below (2-layer) fictitious moho down to 670
      if (idoubling(ispec) == IFLAG_220_80 &
        .or. idoubling(ispec) == IFLAG_80_MOHO &
        .or. idoubling(ispec) == IFLAG_670_220) then
        elem_is_tiso = .true.
      endif
    else if (USE_OLD_VERSION_7_0_0_FORMAT) then
      ! assigns TI to elements in mantle elements just below actual moho down to 670
      if (idoubling(ispec) == IFLAG_670_220 &
          .or. idoubling(ispec) == IFLAG_220_80 &
          .or. idoubling(ispec) == IFLAG_80_MOHO &
          .or. (idoubling(ispec) == IFLAG_CRUST .and. elem_in_mantle) ) then
        elem_is_tiso = .true.
      endif
    else
      ! assigns TI to elements in crust and mantle down to 670
      if (idoubling(ispec) == IFLAG_670_220 &
          .or. idoubling(ispec) == IFLAG_220_80 &
          .or. idoubling(ispec) == IFLAG_80_MOHO &
          .or. idoubling(ispec) == IFLAG_CRUST) then
        elem_is_tiso = .true.
      endif
    endif

  case default
    ! default reference models
    ! for example, PREM assigns transverse isotropy between Moho and 220km
    if (USE_OLD_VERSION_5_1_5_FORMAT) then
      ! assigns TI only to elements below (2-layer) fictitious moho
      if (idoubling(ispec) == IFLAG_220_80 &
          .or. idoubling(ispec) == IFLAG_80_MOHO) then
        ! default case for PREM reference models:
        ! models use only transverse isotropy between moho and 220 km depth
        elem_is_tiso = .true.
        ! checks mantle flag to be sure
        if (elem_in_mantle .eqv. .false. ) stop 'Error mantle flag confused between moho and 220'
      endif
    else if (USE_OLD_VERSION_7_0_0_FORMAT) then
      ! assigns TI to elements in mantle elements just below actual moho
      if (idoubling(ispec) == IFLAG_220_80 &
          .or. idoubling(ispec) == IFLAG_80_MOHO &
          .or. (idoubling(ispec) == IFLAG_CRUST .and. elem_in_mantle )) then
        ! default case for PREM reference models:
        ! models use only transverse isotropy between moho and 220 km depth
        elem_is_tiso = .true.
      endif
    else
      ! assigns TI to elements in crust and mantle elements (down to 220),
      ! to allow for tiso in crust and below actual moho (especially for oceanic crusts);
      ! the crustal models will decide if model parameters are tiso or iso
      if (idoubling(ispec) == IFLAG_220_80 &
          .or. idoubling(ispec) == IFLAG_80_MOHO &
          .or. idoubling(ispec) == IFLAG_CRUST) then
        ! default case for PREM reference models:
        ! models use only transverse isotropy between moho and 220 km depth
        elem_is_tiso = .true.
      endif
    endif
  end select

  ! note: for this new implementation, we will allow to have crustal elements computed as tiso.
  !       the crustal model will decide if the vpv,vph,vsv,vsh,eta values are actually isotropic or not.
  !       the additional costs for treating these elements as tiso have become quite small.
  !
  !       older implementations used also the flag elem_in_mantle to distinguish elements in the IFLAG_CRUST layer
  !       which were below the actual moho (needed for the case of 3D crustal models and thin oceanic moho).
  !       this flag however becomes obsolete in this newer implementation. it is kept for backward compatibility,
  !       in case one wants to run a simulation with USE_OLD_VERSION_5_1_5_FORMAT or USE_OLD_VERSION_7_0_0_FORMAT turned on.

  ! 3D mantle model specific additions
  select case (THREE_D_MODEL)

  case (THREE_D_MODEL_S362WMANI)
    ! S362wmani allows radial anisotropy throughout the mantle
    ! allows tiso down to CMB
    if (idoubling(ispec) == IFLAG_MANTLE_NORMAL &
        .or. elem_in_mantle) then
      elem_is_tiso = .true.
    endif

  case (THREE_D_MODEL_SGLOBE)
    ! sgloberani_aniso model allows for additional tiso throughout the mantle
    ! (from above based on PREM, it will have tiso already set from crust down to 220)
    if (idoubling(ispec) == IFLAG_MANTLE_NORMAL &
        .or. idoubling(ispec) == IFLAG_670_220 &
        .or. elem_in_mantle) then
      elem_is_tiso = .true.
    endif

    ! note: THREE_D_MODEL_SGLOBE_ISO
    !       sgloberani_iso model based on PREM, it will have tiso already set from crust down to 220

  case default
    ! nothing special to add
    continue

  end select

  !debug
  !if (myrank == 0) print *,'  element ',ispec,' tiso flag: ',elem_is_tiso

  end subroutine compute_element_tiso_flag
