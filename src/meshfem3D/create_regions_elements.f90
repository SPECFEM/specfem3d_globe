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

  subroutine create_regions_elements(iregion_code,ipass, &
                                     NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                                     offset_proc_xi,offset_proc_eta)

! creates all elements belonging to different regions of the mesh

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN

  use meshfem3D_par, only: &
    nspec,idoubling,is_on_a_slice_edge, &
    xstore,ystore,zstore, &
    IMAIN,myrank, &
    IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE,IFLAG_IN_FICTITIOUS_CUBE, &
    NPROC_XI,NPROC_ETA,NCHUNKS, &
    INCLUDE_CENTRAL_CUBE,R_CENTRAL_CUBE, &
    MAX_NUMBER_OF_MESH_LAYERS,MAX_NUM_REGIONS,NB_SQUARE_CORNERS, &
    rmins,rmaxs,iproc_xi,iproc_eta,ichunk,NEX_XI, &
    rotation_matrix,ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD, &
    ratio_sampling_array,doubling_index,this_region_has_a_doubling, &
    ratio_divide_central_cube,CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, &
    ner,r_top,r_bottom

  use meshfem3D_models_par, only: &
    SAVE_BOUNDARY_MESH,SUPPRESS_CRUSTAL_MESH,REGIONAL_MOHO_MESH, &
    TRANSVERSE_ISOTROPY

  use regions_mesh_par, only: &
    xigll,yigll,zigll,iaddx,iaddy,iaddz,shape3D,dershape2D_bottom

  use regions_mesh_par2, only: &
    USE_ONE_LAYER_SB,ispec_is_tiso,ifirst_region,ilast_region,perm_layer,NUMBER_OF_MESH_LAYERS, &
    nspec_ani,nspec_actually,nspec_stacey,iMPIcut_xi,iMPIcut_eta, &
    stretch_tab, &
    rhostore,dvpstore,rho_vp,rho_vs, &
    kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
    xixstore,xiystore,xizstore, &
    etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,iboun

! boundary mesh
  use regions_mesh_par2, only: &
    r_moho,r_400,r_670,NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,nex_eta_moho, &
    ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot,ibelm_670_top,ibelm_670_bot, &
    normal_moho,normal_400,normal_670,jacobian2D_moho,jacobian2D_400,jacobian2D_670, &
    ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot

  implicit none

  integer,intent(in) :: iregion_code
  integer,intent(in) :: ipass

  integer,intent(in) :: NEX_PER_PROC_XI,NEX_PER_PROC_ETA

  integer,intent(in) :: offset_proc_xi,offset_proc_eta

  ! local parameters
  integer :: ispec,nspec_tiso
  ! parameters needed to store the radii of the grid points in the spherically symmetric Earth
  double precision :: rmin,rmax
  integer :: ner_without_doubling,ilayer,ilayer_loop
  ! timing
  double precision, external :: wtime
  !double precision :: time_start,tCPU
  integer,dimension(8) :: tval

  ! initializes flags for transverse isotropic elements
  ispec_is_tiso(:) = .false.

  ! get MPI starting time
  !time_start = wtime()

  ! loop on all the layers in this region of the mesh
  ispec = 0 ! counts all the elements in this region of the mesh
  do ilayer_loop = ifirst_region,ilast_region

    ilayer = perm_layer(ilayer_loop)

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  creating layer ',ilayer_loop-ifirst_region+1, &
                                   'out of ',ilast_region-ifirst_region+1
      call flush_IMAIN()
    endif

    ! determine the radii that define the shell
    rmin = rmins(ilayer)
    rmax = rmaxs(ilayer)

    ner_without_doubling = ner(ilayer)

    ! if there is a doubling at the top of this region, we implement it in the last two layers of elements
    ! and therefore we suppress two layers of regular elements here
    USE_ONE_LAYER_SB = .false.
    if (this_region_has_a_doubling(ilayer)) then
      if (ner(ilayer) == 1) then
        ner_without_doubling = ner_without_doubling - 1
        USE_ONE_LAYER_SB = .true.
      else
        ner_without_doubling = ner_without_doubling - 2
        USE_ONE_LAYER_SB = .false.
      endif
    endif

    ! regular mesh elements
    call create_regular_elements(ilayer,ichunk,ispec,ipass, &
                    ifirst_region,ilast_region,iregion_code, &
                    nspec,NCHUNKS,NUMBER_OF_MESH_LAYERS, &
                    NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                    ner_without_doubling,ner,ratio_sampling_array,r_top,r_bottom, &
                    xstore,ystore,zstore, &
                    iaddx,iaddy,iaddz,xigll,yigll,zigll, &
                    shape3D,dershape2D_bottom, &
                    INCLUDE_CENTRAL_CUBE, &
                    rmin,rmax,r_moho,r_400,r_670, &
                    rhostore,dvpstore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                    nspec_ani,c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                    nspec_actually,xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                    gammaxstore,gammaystore,gammazstore, &
                    nspec_stacey,rho_vp,rho_vs,iboun,iMPIcut_xi,iMPIcut_eta, &
                    ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
                    rotation_matrix,idoubling,doubling_index,USE_ONE_LAYER_SB, &
                    stretch_tab, &
                    NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,nex_eta_moho, &
                    ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot,ibelm_670_top,ibelm_670_bot, &
                    normal_moho,normal_400,normal_670,jacobian2D_moho,jacobian2D_400,jacobian2D_670, &
                    ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top, &
                    ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot, &
                    ispec_is_tiso)


    ! mesh doubling elements
    if (this_region_has_a_doubling(ilayer) ) &
      call create_doubling_elements(ilayer,ichunk,ispec,ipass, &
                    ifirst_region,ilast_region,iregion_code, &
                    nspec,NCHUNKS,NUMBER_OF_MESH_LAYERS, &
                    NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                    ner,ratio_sampling_array,r_top,r_bottom, &
                    xstore,ystore,zstore,xigll,yigll,zigll, &
                    shape3D,dershape2D_bottom, &
                    INCLUDE_CENTRAL_CUBE, &
                    rmin,rmax,r_moho,r_400,r_670, &
                    rhostore,dvpstore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                    nspec_ani,c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                    nspec_actually,xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                    gammaxstore,gammaystore,gammazstore, &
                    nspec_stacey,rho_vp,rho_vs,iboun,iMPIcut_xi,iMPIcut_eta, &
                    ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
                    rotation_matrix,idoubling,doubling_index,USE_ONE_LAYER_SB, &
                    NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,nex_eta_moho, &
                    ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot,ibelm_670_top,ibelm_670_bot, &
                    normal_moho,normal_400,normal_670,jacobian2D_moho,jacobian2D_400,jacobian2D_670, &
                    ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top, &
                    ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot, &
                    CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,offset_proc_xi,offset_proc_eta, &
                    ispec_is_tiso)

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

  deallocate(stretch_tab)
  deallocate(perm_layer)
  deallocate(jacobian2D_moho,jacobian2D_400,jacobian2D_670)

  if (myrank == 0 ) write(IMAIN,*)

  ! define central cube in inner core
  if (INCLUDE_CENTRAL_CUBE .and. iregion_code == IREGION_INNER_CORE) then
    ! user output
    if (myrank == 0 ) write(IMAIN,*) '  creating central cube'

    call create_central_cube(ichunk,ispec,iaddx,iaddy,iaddz,ipass, &
                        nspec,NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,R_CENTRAL_CUBE, &
                        iproc_xi,iproc_eta,NPROC_XI,NPROC_ETA,ratio_divide_central_cube, &
                        iMPIcut_xi,iMPIcut_eta,iboun, &
                        idoubling,iregion_code,xstore,ystore,zstore, &
                        shape3D,rmin,rmax,rhostore,dvpstore, &
                        kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                        xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                        gammaxstore,gammaystore,gammazstore,nspec_actually, &
                        c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                        c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                        c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                        nspec_ani,nspec_stacey, &
                        rho_vp,rho_vs,xigll,yigll,zigll, &
                        ispec_is_tiso)
  endif

  ! check total number of spectral elements created
  if (ispec /= nspec) call exit_MPI(myrank,'ispec should equal nspec')

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

