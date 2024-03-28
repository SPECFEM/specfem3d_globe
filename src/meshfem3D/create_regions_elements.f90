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

  subroutine create_regions_elements(ipass, &
                                     NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                                     offset_proc_xi,offset_proc_eta)

! creates all elements belonging to different regions of the mesh
! this will create a cubed-sphere mesh that accommodates Moho variations, topography and ellipticity.

  use constants, only: &
    myrank,IMAIN, &
    IREGION_CRUST_MANTLE,IREGION_INNER_CORE,IREGION_OUTER_CORE, &
    IREGION_TRINFINITE,IREGION_INFINITE, &
    IFLAG_IN_FICTITIOUS_CUBE

  use meshfem_par, only: &
    nspec,iregion_code, &
    idoubling,is_on_a_slice_edge, &
    NPROC_XI,NPROC_ETA,NCHUNKS, &
    R_CENTRAL_CUBE,INCLUDE_CENTRAL_CUBE, &
    rmins,rmaxs,iproc_xi,iproc_eta,ichunk,NEX_XI, &
    rotation_matrix,ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD, &
    this_region_has_a_doubling, &
    CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, &
    ner_mesh_layers

  use meshfem_models_par, only: &
    TRANSVERSE_ISOTROPY

  use regions_mesh_par2, only: &
    USE_ONE_LAYER_SB,ispec_is_tiso,ifirst_region,ilast_region,perm_layer,NUMBER_OF_MESH_LAYERS, &
    iMPIcut_xi,iMPIcut_eta, &
    stretch_tab, &
    iboun

! boundary mesh
  use regions_mesh_par2, only: &
    r_moho,r_400,r_670,NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,nex_eta_moho, &
    ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot,ibelm_670_top,ibelm_670_bot, &
    normal_moho,normal_400,normal_670,jacobian2D_moho,jacobian2D_400,jacobian2D_670, &
    ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot

  implicit none

  integer,intent(in) :: ipass

  integer,intent(in) :: NEX_PER_PROC_XI,NEX_PER_PROC_ETA

  integer,intent(in) :: offset_proc_xi,offset_proc_eta

  ! local parameters
  integer :: ispec_count,ispec_count_all,nspec_tiso

  ! parameters needed to store the radii of the grid points in the spherically symmetric Earth
  double precision :: rmin,rmax
  integer :: ner_without_doubling,ilayer,ilayer_loop

  ! timing
  double precision, external :: wtime
  !double precision :: time_start,tCPU
  integer,dimension(8) :: tval

  ! tolerance value to zero layer thickness
  double precision, parameter :: TOLERANCE_LAYER_THICKNESS = 1.d-9

  ! initializes flags for transverse isotropic elements
  ispec_is_tiso(:) = .false.
  nspec_tiso = 0

  ! get MPI starting time
  !time_start = wtime()

  ! counts all the elements in this region of the mesh
  ispec_count = 0

  ! checks if anything to do
  if (ifirst_region == 0 .and. ilast_region == 0) return

  ! loop on all the layers in this region of the mesh
  do ilayer_loop = ifirst_region,ilast_region

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  creating layer ',ilayer_loop-ifirst_region+1, &
                                   'out of ',ilast_region-ifirst_region+1
      call flush_IMAIN()
    endif

    ! creates all elements per layer
    select case(iregion_code)
    case (IREGION_INNER_CORE,IREGION_OUTER_CORE,IREGION_CRUST_MANTLE)
      ! default regions for global mesh (inner core, outer core, crust/mantle)
      ! current layer
      ilayer = perm_layer(ilayer_loop)

      ! determine the radii that define the shell
      rmin = rmins(ilayer)
      rmax = rmaxs(ilayer)

      ! skips layers with zero thickness
      if (abs(rmax - rmin) < TOLERANCE_LAYER_THICKNESS) cycle

      ner_without_doubling = ner_mesh_layers(ilayer)

      ! checks if anything to do
      if (ner_without_doubling == 0) cycle

      ! if there is a doubling at the top of this region, we implement it in the last two layers of elements
      ! and therefore we suppress two layers of regular elements here
      USE_ONE_LAYER_SB = .false.
      if (this_region_has_a_doubling(ilayer)) then
        if (ner_mesh_layers(ilayer) == 1) then
          ner_without_doubling = ner_without_doubling - 1
          USE_ONE_LAYER_SB = .true.
        else
          ner_without_doubling = ner_without_doubling - 2
          USE_ONE_LAYER_SB = .false.
        endif
      endif

      ! regular mesh elements
      call create_regular_elements(ilayer,ichunk,ispec_count,ipass, &
                                   ifirst_region,ilast_region,iregion_code, &
                                   nspec,NCHUNKS,NUMBER_OF_MESH_LAYERS, &
                                   NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                                   ner_without_doubling, &
                                   INCLUDE_CENTRAL_CUBE, &
                                   rmin,rmax,r_moho,r_400,r_670, &
                                   iMPIcut_xi,iMPIcut_eta, &
                                   ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
                                   rotation_matrix,idoubling,USE_ONE_LAYER_SB, &
                                   stretch_tab, &
                                   NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,nex_eta_moho, &
                                   ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot,ibelm_670_top,ibelm_670_bot, &
                                   normal_moho,normal_400,normal_670,jacobian2D_moho,jacobian2D_400,jacobian2D_670, &
                                   ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top, &
                                   ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot, &
                                   ispec_is_tiso)


      ! mesh doubling elements
      if (this_region_has_a_doubling(ilayer) ) then
        call create_doubling_elements(ilayer,ichunk,ispec_count,ipass, &
                                      ifirst_region,ilast_region,iregion_code, &
                                      nspec,NCHUNKS,NUMBER_OF_MESH_LAYERS, &
                                      NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                                      INCLUDE_CENTRAL_CUBE, &
                                      rmin,rmax,r_moho,r_400,r_670, &
                                      iMPIcut_xi,iMPIcut_eta, &
                                      ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
                                      rotation_matrix,idoubling,USE_ONE_LAYER_SB, &
                                      NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,nex_eta_moho, &
                                      ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot,ibelm_670_top,ibelm_670_bot, &
                                      normal_moho,normal_400,normal_670,jacobian2D_moho,jacobian2D_400,jacobian2D_670, &
                                      ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top, &
                                      ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot, &
                                      CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,offset_proc_xi,offset_proc_eta, &
                                      ispec_is_tiso)
      endif

    case (IREGION_TRINFINITE,IREGION_INFINITE)
      ! infinite-element mesh regions
      ! current layer - no need for layer permutations
      ilayer = ilayer_loop
      ! mesh spectral-infinite elements
      call SIEM_mesh_create_elements(ilayer,ispec_count,ipass,iregion_code)

    case default
      call exit_MPI(myrank,'Invalid region code for creating regions elements')
    end select

    ! user output
    if (myrank == 0) then
      ! time estimate
      !tCPU = wtime() - time_start

      ! outputs current time on system
      call date_and_time(VALUES=tval)

      ! debug: outputs remaining time (poor estimation)
      !tCPU = (1.0-(ilayer_loop-ifirst_region+1.0)/(ilast_region-ifirst_region+1.0)) &
      !          /(ilayer_loop-ifirst_region+1.0)/(ilast_region-ifirst_region+1.0)*tCPU*10.0

      ! user output
      write(IMAIN,'(a,f5.1,a,a,i2.2,a,i2.2,a,i2.2,a)') &
        "    ",(ilayer_loop-ifirst_region+1.0)/(ilast_region-ifirst_region+1.0) * 100.0,"%", &
        "    current clock (NOT elapsed) time is: ",tval(5),"h ",tval(6),"min ",tval(7),"sec"

      ! flushes I/O buffer
      call flush_IMAIN()
    endif

  enddo ! of ilayer_loop

  ! stats
  call sum_all_i(ispec_count,ispec_count_all)

  ! user output
  if (myrank == 0 ) then
    write(IMAIN,*) '  layers done'
    write(IMAIN,*)
    write(IMAIN,*) '  number of elements (per slice)        = ',ispec_count
    write(IMAIN,*) '  total number of elements (all slices) = ',ispec_count_all
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! frees temporary arrays
  deallocate(stretch_tab)
  deallocate(perm_layer)
  deallocate(jacobian2D_moho,jacobian2D_400,jacobian2D_670)

  ! define central cube in inner core
  if (INCLUDE_CENTRAL_CUBE .and. iregion_code == IREGION_INNER_CORE) then
    ! user output
    if (myrank == 0 ) then
      write(IMAIN,*) '  creating central cube'
      call flush_IMAIN()
    endif

    call create_central_cube(ichunk,ispec_count,ipass, &
                             nspec,NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                             iproc_xi,iproc_eta,NPROC_XI,NPROC_ETA, &
                             iMPIcut_xi,iMPIcut_eta, &
                             idoubling,iregion_code, &
                             rmin,rmax,R_CENTRAL_CUBE, &
                             ispec_is_tiso)

    ! user output
    if (myrank == 0 ) then
      write(IMAIN,*) '  central cube done'
      write(IMAIN,*) '    total number of elements = ',ispec_count
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! check total number of spectral elements created
  if (ispec_count /= nspec) then
    print *,'Error: rank ',myrank,' has created invalid total number of elements ',ispec_count,'; should be ',nspec
    call exit_MPI(myrank,'ispec_count should equal nspec')
  endif

  ! if any of these flags is true, the element is on a communication edge
  ! this is not enough because it can also be in contact by an edge or a corner but not a full face
  ! therefore we will have to fix array "is_on_a_slice_edge" later in the solver to take this into account
  is_on_a_slice_edge(:) = &
      iMPIcut_xi(1,:) .or. iMPIcut_xi(2,:) .or. &
      iMPIcut_eta(1,:) .or. iMPIcut_eta(2,:) .or. &
      iboun(1,:) .or. iboun(2,:) .or. &
      iboun(3,:) .or. iboun(4,:) .or. &
      iboun(5,:) .or. iboun(6,:)

  ! no need to count fictitious elements on the edges
  ! for which communications cannot be overlapped with calculations
  where(idoubling == IFLAG_IN_FICTITIOUS_CUBE) is_on_a_slice_edge = .false.

  ! checks transverse isotropic elements
  if (ipass == 2) then
    ! count number of anisotropic elements in current region
    ! should be zero in all the regions except in the mantle
    nspec_tiso = count(ispec_is_tiso(:))

    ! checks number of anisotropic elements found in the mantle
    if (iregion_code /= IREGION_CRUST_MANTLE .and. nspec_tiso /= 0 ) &
      call exit_MPI(myrank,'found anisotropic elements outside of the mantle')
    if (TRANSVERSE_ISOTROPY) then
      if (iregion_code == IREGION_CRUST_MANTLE .and. nspec_tiso == 0) &
        call exit_MPI(myrank,'found no anisotropic elements in the mantle')
    endif
  endif

  end subroutine create_regions_elements

