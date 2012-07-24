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



  subroutine create_regions_mesh(iregion_code,ibool,idoubling,is_on_a_slice_edge, &
                          xstore,ystore,zstore, &
                          nspec, &
                          nglob_theor,npointot, &
                          NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                          NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                          NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
                          NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                          offset_proc_xi,offset_proc_eta, &
                          ipass)

! creates the different regions of the mesh

  use meshfem3D_par,only: &
    IMAIN,volume_total,addressing,ichunk_slice,iproc_xi_slice,iproc_eta_slice, &
    myrank,LOCAL_PATH, &
    IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE,IFLAG_IN_FICTITIOUS_CUBE, &
    NPROC,NPROCTOT,NPROC_XI,NPROC_ETA,NCHUNKS, &
    SAVE_MESH_FILES,INCLUDE_CENTRAL_CUBE,ABSORBING_CONDITIONS, &
    R_CENTRAL_CUBE,RICB,RHO_OCEANS,RCMB,R670,RMOHO,RMOHO_FICTITIOUS_IN_MESHER,&
    RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
    MAX_NUMBER_OF_MESH_LAYERS,MAX_NUM_REGIONS,NB_SQUARE_EDGES_ONEDIR,NB_SQUARE_CORNERS, &
    rmins,rmaxs,iproc_xi,iproc_eta,ichunk,NEX_XI, &
    rotation_matrix,ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD, &
    ratio_sampling_array,doubling_index,this_region_has_a_doubling, &
    ratio_divide_central_cube,CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, &
    NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER, &
    ner,r_bottom,r_top

  use meshfem3D_models_par,only: &
    ATTENUATION,ANISOTROPIC_INNER_CORE,ANISOTROPIC_3D_MANTLE, &
    SAVE_BOUNDARY_MESH,ONE_CRUST,CASE_3D,SUPPRESS_CRUSTAL_MESH,REGIONAL_MOHO_MESH, &
    OCEANS,TRANSVERSE_ISOTROPY,HETEROGEN_3D_MANTLE,HONOR_1D_SPHERICAL_MOHO, &
    ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
    AM_V,nspl,rspl,espl,espl2

  use create_regions_mesh_par
  use create_regions_mesh_par2

  implicit none

  ! code for the four regions of the mesh
  integer :: iregion_code

  ! correct number of spectral elements in each block depending on chunk type
  integer :: nspec

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: idoubling

  ! this for non blocking MPI
  logical, dimension(nspec) :: is_on_a_slice_edge

  ! arrays with the mesh in double precision
  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

  integer :: nglob_theor,npointot

  integer :: NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer :: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX
  integer :: NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX
  integer :: NSPEC2D_BOTTOM,NSPEC2D_TOP

  integer :: offset_proc_xi,offset_proc_eta

  ! now perform two passes in this part to be able to save memory
  integer,intent(in) :: ipass

  ! local parameters
  integer :: ier
  integer :: nglob
  ! saved for second call
  integer, save :: npoin2D_xi,npoin2D_eta
  ! check area and volume of the final mesh
  double precision :: area_local_bottom,area_local_top
  double precision :: volume_local

  ! user output
  if(myrank == 0 ) then
    select case(ipass)
    case(1)
      write(IMAIN,*) 'first pass'
    case(2)
      write(IMAIN,*) 'second pass'
    case default
      call exit_MPI(myrank,'error ipass value in create_regions_mesh')
    end select
  endif

  ! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,iregion_code,LOCAL_PATH)

  ! initializes arrays
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...allocating arrays '
  endif
  call crm_allocate_arrays(iregion_code,nspec, &
                                NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                                NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                                ipass)


  ! initialize number of layers
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...setting up layers '
  endif
  call crm_setup_layers(iregion_code,nspec,ipass, &
                              xstore,ystore,zstore,ibool,idoubling, &
                              NEX_PER_PROC_ETA,is_on_a_slice_edge)

  !  creates mesh elements
  call sync_all()
  if( myrank == 0) then
    write(IMAIN,*) '  ...creating mesh elements '
  endif
  call crm_create_elements(iregion_code,nspec,ipass, &
                              xstore,ystore,zstore,idoubling, &
                              NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                              is_on_a_slice_edge, &
                              offset_proc_xi,offset_proc_eta)


  ! only create global addressing and the MPI buffers in the first pass
  select case(ipass)
  case( 1 )
    ! creates ibool index array for projection from local to global points
    call sync_all()
    if( myrank == 0) then
      write(IMAIN,*) '  ...creating global addressing'
    endif
    call crm_setup_indexing(ibool,xstore,ystore,zstore, &
                                nspec,nglob_theor,npointot)


    ! create MPI buffers
    call sync_all()
    if( myrank == 0) then
      write(IMAIN,*) '  ...creating MPI buffers'
    endif
    call crm_setup_mpi_buffers(npointot,nspec,ibool,idoubling, &
                                  xstore,ystore,zstore,iregion_code, &
                                  npoin2D_xi,npoin2D_eta)


    ! Stacey
    if(NCHUNKS /= 6) then
      call get_absorb(myrank,prname,iboun,nspec,nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta, &
                         NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM)

      deallocate(nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta)
    endif
    
  ! only create mass matrix and save all the final arrays in the second pass
  case( 2 )
    ! precomputes jacobian for 2d absorbing boundary surfaces
    call sync_all()
    if( myrank == 0) then
      write(IMAIN,*) '  ...precomputing jacobian'
    endif
    call get_jacobian_boundaries(myrank,iboun,nspec,xstore,ystore,zstore, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
              nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
              jacobian2D_xmin,jacobian2D_xmax, &
              jacobian2D_ymin,jacobian2D_ymax, &
              jacobian2D_bottom,jacobian2D_top, &
              normal_xmin,normal_xmax, &
              normal_ymin,normal_ymax, &
              normal_bottom,normal_top, &
              NSPEC2D_BOTTOM,NSPEC2D_TOP, &
              NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,&
              xigll,yigll,zigll)

    ! create chunk buffers if more than one chunk
    call sync_all()
    if( myrank == 0) then
      write(IMAIN,*) '  ...creating chunk buffers'
    endif
    if(NCHUNKS > 1) then
      call create_chunk_buffers(iregion_code,nspec,ibool,idoubling, &
                              xstore,ystore,zstore, &
                              nglob_theor, &
                              NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                              NPROC_XI,NPROC_ETA, &
                              NPROC,NPROCTOT, &
                              NGLOB1D_RADIAL_CORNER,maxval(NGLOB1D_RADIAL_CORNER(iregion_code,:)), &
                              NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
                              myrank,LOCAL_PATH,addressing, &
                              ichunk_slice,iproc_xi_slice,iproc_eta_slice,NCHUNKS, &
                              nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                              ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax)
    else
      if(myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'only one chunk, no need to create chunk buffers'
        write(IMAIN,*)
      endif
    endif

    !uncomment: adds model smoothing for point profile models
    !    if( THREE_D_MODEL == THREE_D_MODEL_PPM ) then
    !     call smooth_model(myrank, nproc_xi,nproc_eta,&
    !        rho_vp,rho_vs,nspec_stacey, &
    !        iregion_code,xixstore,xiystore,xizstore, &
    !        etaxstore,etaystore,etazstore, &
    !        gammaxstore,gammaystore,gammazstore, &
    !        xstore,ystore,zstore,rhostore,dvpstore, &
    !        kappavstore,kappahstore,muvstore,muhstore,eta_anisostore,&
    !        nspec,HETEROGEN_3D_MANTLE, &
    !        NEX_XI,NCHUNKS,ABSORBING_CONDITIONS,PPM_V )

    ! creates mass matrix
    call sync_all()
    if( myrank == 0) then
      write(IMAIN,*) '  ...creating mass matrix'
    endif

    ! allocates mass matrices in this slice (will be fully assembled in the solver)
    !
    ! in the case of stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
    ! on Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
    ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
    !
    ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
    ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be obsolete

    ! copy the theoretical number of points for the second pass
    nglob = nglob_theor

    if(NCHUNKS /= 6 .and. ABSORBING_CONDITIONS) then
       select case(iregion_code)
       case( IREGION_CRUST_MANTLE )
          nglob_xy = nglob
       case( IREGION_INNER_CORE, IREGION_OUTER_CORE )
          nglob_xy = 1
       endselect
    else
       nglob_xy = 1
    endif

    allocate(rmassx(nglob_xy),stat=ier)
    if(ier /= 0) stop 'error in allocate 21'
    allocate(rmassy(nglob_xy),stat=ier)
    if(ier /= 0) stop 'error in allocate 21'
    allocate(rmassz(nglob),stat=ier)
    if(ier /= 0) stop 'error in allocate 21'

    ! allocates ocean load mass matrix as well if oceans
    if(OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) then
      nglob_oceans = nglob
    else
      ! allocate dummy array if no oceans
      nglob_oceans = 1
    endif
    allocate(rmass_ocean_load(nglob_oceans),stat=ier)
    if(ier /= 0) stop 'error in allocate 22'

    ! creating mass matrices in this slice (will be fully assembled in the solver)
    call create_mass_matrices(myrank,nspec,idoubling,wxgll,wygll,wzgll,ibool, &
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

    ! save the binary files
    call sync_all()
    if( myrank == 0) then
      write(IMAIN,*) '  ...saving binary files'
    endif
    call save_arrays_solver(myrank,rho_vp,rho_vs,nspec_stacey, &
                  prname,iregion_code,xixstore,xiystore,xizstore, &
                  etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
                  xstore,ystore,zstore,rhostore,dvpstore, &
                  kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                  nspec_ani,c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                  c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                  c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                  ibool,idoubling,is_on_a_slice_edge,nglob_xy,nglob, &
                  rmassx,rmassy,rmassz,rmass_ocean_load,nglob_oceans, &
                  ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                  nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                  normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top, &
                  jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax, &
                  jacobian2D_bottom,jacobian2D_top,nspec, &
                  NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                  TRANSVERSE_ISOTROPY,HETEROGEN_3D_MANTLE,ANISOTROPIC_3D_MANTLE, &
                  ANISOTROPIC_INNER_CORE,OCEANS, &
                  tau_s,tau_e_store,Qmu_store,T_c_source,ATTENUATION, &
                  size(tau_e_store,2),size(tau_e_store,3),size(tau_e_store,4),size(tau_e_store,5),&
                  ABSORBING_CONDITIONS,SAVE_MESH_FILES,ispec_is_tiso)

    deallocate(rmassx,rmassy,rmassz)
    deallocate(rmass_ocean_load)

    ! setup mpi communication interfaces
    call sync_all()
    if( myrank == 0) then
      write(IMAIN,*) '  ...preparing MPI interfaces'
    endif    
    call create_MPI_interfaces(iregion_code)

    ! sets up inner/outer element arrays
    call sync_all()
    if( myrank == 0) then
      write(IMAIN,*) '  ...element inner/outer separation '
    endif
    call setup_inner_outer(iregion_code)
    
    ! sets up mesh coloring
    call sync_all()
    if( myrank == 0) then
      write(IMAIN,*) '  ...element mesh coloring '
    endif
    call setup_color_perm(iregion_code)

    ! saves MPI interface infos
    call cmi_save_MPI_interfaces(iregion_code)

    ! frees memory
    call cmi_free_MPI_arrays(iregion_code)


    ! boundary mesh
    if (SAVE_BOUNDARY_MESH .and. iregion_code == IREGION_CRUST_MANTLE) then
      ! user output
      call sync_all()
      if( myrank == 0) then
        write(IMAIN,*) '  ...saving boundary mesh files'
      endif
      call crm_save_boundary_mesh()
    endif

    ! compute volume, bottom and top area of that part of the slice
    call compute_volumes(volume_local,area_local_bottom,area_local_top, &
                            nspec,wxgll,wygll,wzgll,xixstore,xiystore,xizstore, &
                            etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
                            NSPEC2D_BOTTOM,jacobian2D_bottom,NSPEC2D_TOP,jacobian2D_top)

    ! computes total area and volume
    call compute_area(myrank,NCHUNKS,iregion_code, &
                              area_local_bottom,area_local_top,&
                              volume_local,volume_total, &
                              RCMB,RICB,R_CENTRAL_CUBE)

    ! create AVS or DX mesh data for the slices
    if(SAVE_MESH_FILES) then
      ! user output
      call sync_all()
      if( myrank == 0) then
        write(IMAIN,*) '  ...saving AVS mesh files'
      endif
      call crm_save_mesh_files(nspec,npointot,iregion_code,ibool,idoubling, &
                                  xstore,ystore,zstore)
    endif

  case default
    stop 'there cannot be more than two passes in mesh creation'

  end select  ! end of test if first or second pass


  ! deallocate these arrays after each pass
  ! because they have a different size in each pass to save memory
  deallocate(xixstore,xiystore,xizstore)
  deallocate(etaxstore,etaystore,etazstore)
  deallocate(gammaxstore,gammaystore,gammazstore)

  ! deallocate arrays
  deallocate(rhostore,dvpstore,kappavstore,kappahstore)
  deallocate(muvstore,muhstore)
  deallocate(eta_anisostore)
  deallocate(ispec_is_tiso)
  deallocate(c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
            c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
            c36store,c44store,c45store,c46store,c55store,c56store,c66store)
  deallocate(iboun)
  deallocate(ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax)
  deallocate(ibelm_bottom,ibelm_top)
  deallocate(jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax)
  deallocate(jacobian2D_bottom,jacobian2D_top)
  deallocate(normal_xmin,normal_xmax,normal_ymin,normal_ymax)
  deallocate(normal_bottom,normal_top)
  deallocate(iMPIcut_xi,iMPIcut_eta)

  deallocate(rho_vp,rho_vs)
  deallocate(Qmu_store)
  deallocate(tau_e_store)
  
  deallocate(ibelm_moho_top,ibelm_moho_bot)
  deallocate(ibelm_400_top,ibelm_400_bot)
  deallocate(ibelm_670_top,ibelm_670_bot)
  deallocate(normal_moho,normal_400,normal_670)


  ! user output
  if(myrank == 0 ) write(IMAIN,*)

  end subroutine create_regions_mesh

!
!-------------------------------------------------------------------------------------------------
!

  subroutine crm_allocate_arrays(iregion_code,nspec, &
                                NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                                NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                                ipass)

  use meshfem3D_par,only: &
    myrank,NGLLX,NGLLY,NGLLZ,NDIM, &
    IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE, &
    NCHUNKS

  use meshfem3D_models_par,only: &
    ATTENUATION,ATTENUATION_3D,ANISOTROPIC_INNER_CORE,ANISOTROPIC_3D_MANTLE, &
    SAVE_BOUNDARY_MESH,AM_V

  use create_regions_mesh_par2

  implicit none

  integer,intent(in) :: iregion_code,nspec
  integer,intent(in) :: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX
  integer,intent(in) :: NSPEC2D_BOTTOM,NSPEC2D_TOP
  integer,intent(in) :: ipass

  ! local parameters
  integer :: ier

  ! New Attenuation definition on all GLL points
  ! Attenuation
  if (ATTENUATION) then
    T_c_source = AM_V%QT_c_source
    tau_s(:)   = AM_V%Qtau_s(:)
    nspec_att = nspec

    if (ATTENUATION_3D) then
      ! attenuation arrays are fully 3D
      allocate(Qmu_store(NGLLX,NGLLY,NGLLZ,nspec_att), &
              tau_e_store(N_SLS,NGLLX,NGLLY,NGLLZ,nspec_att),stat=ier)
    else if (.not. ATTENUATION_3D) then
      ! save some memory in case of 1D attenuation models
      allocate(Qmu_store(1,1,1,nspec_att), &
              tau_e_store(N_SLS,1,1,1,nspec_att),stat=ier)
    endif
  else
    ! allocates dummy size arrays
    nspec_att = 1
    allocate(Qmu_store(NGLLX,NGLLY,NGLLZ,nspec_att), &
            tau_e_store(N_SLS,NGLLX,NGLLY,NGLLZ,nspec_att),stat=ier)
  end if
  if(ier /= 0) stop 'error in allocate 1'

  ! array with model density
  allocate(rhostore(NGLLX,NGLLY,NGLLZ,nspec), &
          dvpstore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if(ier /= 0) stop 'error in allocate 6'

  ! for anisotropy
  allocate(kappavstore(NGLLX,NGLLY,NGLLZ,nspec), &
          muvstore(NGLLX,NGLLY,NGLLZ,nspec), &
          kappahstore(NGLLX,NGLLY,NGLLZ,nspec), &
          muhstore(NGLLX,NGLLY,NGLLZ,nspec), &
          eta_anisostore(NGLLX,NGLLY,NGLLZ,nspec), &
          ispec_is_tiso(nspec),stat=ier)
  if(ier /= 0) stop 'error in allocate 7'

  ! Stacey absorbing boundaries
  if(NCHUNKS /= 6) then
    nspec_stacey = nspec
  else
    nspec_stacey = 1
  endif
  allocate(rho_vp(NGLLX,NGLLY,NGLLZ,nspec_stacey), &
          rho_vs(NGLLX,NGLLY,NGLLZ,nspec_stacey),stat=ier)
  if(ier /= 0) stop 'error in allocate 8'

  ! anisotropy
  if((ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) .or. &
     (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE)) then
    nspec_ani = nspec
  else
    nspec_ani = 1
  endif
  allocate(c11store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c12store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c13store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c14store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c15store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c16store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c22store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c23store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c24store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c25store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c26store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c33store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c34store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c35store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c36store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c44store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c45store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c46store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c55store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c56store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c66store(NGLLX,NGLLY,NGLLZ,nspec_ani),stat=ier)
  if(ier /= 0) stop 'error in allocate 9'

  ! boundary locator
  allocate(iboun(6,nspec),stat=ier)
  if(ier /= 0) stop 'error in allocate 10'

  ! boundary parameters locator
  allocate(ibelm_xmin(NSPEC2DMAX_XMIN_XMAX), &
          ibelm_xmax(NSPEC2DMAX_XMIN_XMAX), &
          ibelm_ymin(NSPEC2DMAX_YMIN_YMAX), &
          ibelm_ymax(NSPEC2DMAX_YMIN_YMAX), &
          ibelm_bottom(NSPEC2D_BOTTOM), &
          ibelm_top(NSPEC2D_TOP),stat=ier)
  if(ier /= 0) stop 'error in allocate 11'

  ! 2-D jacobians and normals
  allocate(jacobian2D_xmin(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX), &
          jacobian2D_xmax(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX), &
          jacobian2D_ymin(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX), &
          jacobian2D_ymax(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX), &
          jacobian2D_bottom(NGLLX,NGLLY,NSPEC2D_BOTTOM), &
          jacobian2D_top(NGLLX,NGLLY,NSPEC2D_TOP),stat=ier)
  if(ier /= 0) stop 'error in allocate 12'

  allocate(normal_xmin(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX), &
          normal_xmax(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX), &
          normal_ymin(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX), &
          normal_ymax(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX), &
          normal_bottom(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM), &
          normal_top(NDIM,NGLLX,NGLLY,NSPEC2D_TOP),stat=ier)
  if(ier /= 0) stop 'error in allocate 13'

  ! Stacey
  if( ipass == 1 .and. NCHUNKS /= 6 ) then
    allocate(nimin(2,NSPEC2DMAX_YMIN_YMAX), &
            nimax(2,NSPEC2DMAX_YMIN_YMAX), &
            njmin(2,NSPEC2DMAX_XMIN_XMAX), &
            njmax(2,NSPEC2DMAX_XMIN_XMAX), &
            nkmin_xi(2,NSPEC2DMAX_XMIN_XMAX), &
            nkmin_eta(2,NSPEC2DMAX_YMIN_YMAX),stat=ier)
    if(ier /= 0) stop 'error in allocate 14'
  endif
  
  ! MPI cut-planes parameters along xi and along eta
  allocate(iMPIcut_xi(2,nspec), &
          iMPIcut_eta(2,nspec),stat=ier)
  if(ier /= 0) stop 'error in allocate 15'

  ! store and save the final arrays only in the second pass
  ! therefore in the first pass some arrays can be allocated with a dummy size
  if(ipass == 1) then
    ACTUALLY_STORE_ARRAYS = .false.
    nspec_actually = 1
  else
    ACTUALLY_STORE_ARRAYS = .true.
    nspec_actually = nspec
  endif
  allocate(xixstore(NGLLX,NGLLY,NGLLZ,nspec_actually), &
          xiystore(NGLLX,NGLLY,NGLLZ,nspec_actually), &
          xizstore(NGLLX,NGLLY,NGLLZ,nspec_actually), &
          etaxstore(NGLLX,NGLLY,NGLLZ,nspec_actually), &
          etaystore(NGLLX,NGLLY,NGLLZ,nspec_actually), &
          etazstore(NGLLX,NGLLY,NGLLZ,nspec_actually), &
          gammaxstore(NGLLX,NGLLY,NGLLZ,nspec_actually), &
          gammaystore(NGLLX,NGLLY,NGLLZ,nspec_actually), &
          gammazstore(NGLLX,NGLLY,NGLLZ,nspec_actually),stat=ier)
  if(ier /= 0) stop 'error in allocate 16'

  ! boundary mesh
  if (ipass == 2 .and. SAVE_BOUNDARY_MESH .and. iregion_code == IREGION_CRUST_MANTLE) then
    NSPEC2D_MOHO = NSPEC2D_TOP
    NSPEC2D_400 = NSPEC2D_MOHO / 4
    NSPEC2D_670 = NSPEC2D_400
  else
    NSPEC2D_MOHO = 1
    NSPEC2D_400 = 1
    NSPEC2D_670 = 1
  endif
  allocate(ibelm_moho_top(NSPEC2D_MOHO),ibelm_moho_bot(NSPEC2D_MOHO), &
          ibelm_400_top(NSPEC2D_400),ibelm_400_bot(NSPEC2D_400), &
          ibelm_670_top(NSPEC2D_670),ibelm_670_bot(NSPEC2D_670), &
          normal_moho(NDIM,NGLLX,NGLLY,NSPEC2D_MOHO), &
          normal_400(NDIM,NGLLX,NGLLY,NSPEC2D_400), &
          normal_670(NDIM,NGLLX,NGLLY,NSPEC2D_670), &
          jacobian2D_moho(NGLLX,NGLLY,NSPEC2D_MOHO), &
          jacobian2D_400(NGLLX,NGLLY,NSPEC2D_400), &
          jacobian2D_670(NGLLX,NGLLY,NSPEC2D_670),stat=ier)
  if(ier /= 0) stop 'error in allocate 17'

  end subroutine crm_allocate_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine crm_setup_layers(iregion_code,nspec,ipass, &
                              xstore,ystore,zstore,ibool,idoubling, &
                              NEX_PER_PROC_ETA,is_on_a_slice_edge)

  use meshfem3D_par,only: &
    myrank,NGLLX,NGLLY,NGLLZ, &
    IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE, &
    R670,RMOHO,R400,RMIDDLE_CRUST,MAX_NUMBER_OF_MESH_LAYERS, &
    ner,r_top,r_bottom

  use meshfem3D_models_par,only: &
    CASE_3D,SUPPRESS_CRUSTAL_MESH,ONE_CRUST,REGIONAL_MOHO_MESH

  use create_regions_mesh_par
  use create_regions_mesh_par2

  implicit none

  integer,intent(in) :: iregion_code,nspec
  integer,intent(in) :: ipass

  ! arrays with the mesh in double precision
  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: idoubling

  ! this for non blocking MPI
  logical, dimension(nspec) :: is_on_a_slice_edge

  integer :: NEX_PER_PROC_ETA

  ! local parameters
  integer :: i,ier

  ! initializes element layers
  call initialize_layers(myrank,ipass,xigll,yigll,zigll,wxgll,wygll,wzgll, &
                        shape3D,dershape3D,shape2D_x,shape2D_y,shape2D_bottom,shape2D_top, &
                        dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                        iaddx,iaddy,iaddz,nspec,xstore,ystore,zstore,ibool,idoubling, &
                        iboun,iMPIcut_xi,iMPIcut_eta,ispec2D_moho_top,ispec2D_moho_bot, &
                        ispec2D_400_top,ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot, &
                        NEX_PER_PROC_ETA,nex_eta_moho,RMOHO,R400,R670,r_moho,r_400,r_670, &
                        ONE_CRUST,NUMBER_OF_MESH_LAYERS,layer_shift, &
                        iregion_code,ifirst_region,ilast_region, &
                        first_layer_aniso,last_layer_aniso,nb_layer_above_aniso,is_on_a_slice_edge)

  ! to consider anisotropic elements first and to build the mesh from the bottom to the top of the region
  allocate (perm_layer(ifirst_region:ilast_region),stat=ier)
  if(ier /= 0) stop 'error in allocate 18'
  perm_layer = (/ (i, i=ilast_region,ifirst_region,-1) /)

  if(iregion_code == IREGION_CRUST_MANTLE) then
    cpt=3
    perm_layer(1)=first_layer_aniso
    perm_layer(2)=last_layer_aniso
    do i = ilast_region,ifirst_region,-1
      if (i/=first_layer_aniso .and. i/=last_layer_aniso) then
        perm_layer(cpt) = i
        cpt=cpt+1
      endif
    enddo
  endif

  ! crustal layer stretching: element layer's top and bottom radii will get stretched when in crust
  ! (number of element layers in crust can vary for different resolutions and 1chunk simulations)
  allocate(stretch_tab(2,ner(1)),stat=ier)
  if(ier /= 0) stop 'error in allocate 19'
  if (CASE_3D .and. iregion_code == IREGION_CRUST_MANTLE .and. .not. SUPPRESS_CRUSTAL_MESH) then
    ! stretching function determines top and bottom of each element layer in the
    ! crust region (between r_top(1) and r_bottom(1)), where ner(1) is the
    ! number of element layers in this crust region

    ! differentiate between regional meshes or global meshes
    if( REGIONAL_MOHO_MESH ) then
      call stretching_function_regional(r_top(1),r_bottom(1),ner(1),stretch_tab)
    else
      call stretching_function(r_top(1),r_bottom(1),ner(1),stretch_tab)
    endif

    ! RMIDDLE_CRUST so far is only used for 1D - models with two layers in the crust
    ! (i.e. ONE_CRUST is set to .false.), those models do not use CASE_3D

    ! all 3D models use this stretching function to honor a 3D crustal model
    ! for those models, we set RMIDDLE_CRUST to the bottom of the first element layer
    ! this value will be used in moho_stretching.f90 to decide whether or not elements
    ! have to be stretched under oceanic crust.
    !
    ! note: stretch_tab uses (dimensionalized) radii from r_top and r_bottom
    !(with stretch_tab( index_radius(1=top,2=bottom), index_layer( 1=first layer, 2=second layer, 3= ...) )
    RMIDDLE_CRUST = stretch_tab(2,1)
  endif

  end subroutine crm_setup_layers

!
!-------------------------------------------------------------------------------------------------
!

  subroutine crm_create_elements(iregion_code,nspec,ipass, &
                              xstore,ystore,zstore,idoubling, &
                              NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                              is_on_a_slice_edge, &
                              offset_proc_xi,offset_proc_eta)

! creates the different regions of the mesh

  use meshfem3D_par,only: &
    IMAIN,volume_total,addressing,ichunk_slice,iproc_xi_slice,iproc_eta_slice, &
    myrank,LOCAL_PATH, &
    IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE,IFLAG_IN_FICTITIOUS_CUBE, &
    NPROC,NPROCTOT,NPROC_XI,NPROC_ETA,NCHUNKS, &
    SAVE_MESH_FILES,INCLUDE_CENTRAL_CUBE,ABSORBING_CONDITIONS, &
    R_CENTRAL_CUBE,RICB,RHO_OCEANS,RCMB,R670,RMOHO,RMOHO_FICTITIOUS_IN_MESHER,&
    RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
    MAX_NUMBER_OF_MESH_LAYERS,MAX_NUM_REGIONS,NB_SQUARE_EDGES_ONEDIR,NB_SQUARE_CORNERS, &
    rmins,rmaxs,iproc_xi,iproc_eta,ichunk,NEX_XI, &
    rotation_matrix,ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD, &
    ratio_sampling_array,doubling_index,this_region_has_a_doubling, &
    ratio_divide_central_cube,CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, &
    ner,r_top,r_bottom

  use meshfem3D_models_par,only: &
    ATTENUATION,ANISOTROPIC_INNER_CORE,ANISOTROPIC_3D_MANTLE, &
    SAVE_BOUNDARY_MESH,ONE_CRUST,CASE_3D,SUPPRESS_CRUSTAL_MESH,REGIONAL_MOHO_MESH, &
    OCEANS,TRANSVERSE_ISOTROPY,HETEROGEN_3D_MANTLE,HONOR_1D_SPHERICAL_MOHO, &
    ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
    AM_V,nspl,rspl,espl,espl2

  use create_regions_mesh_par
  use create_regions_mesh_par2

  implicit none

  integer,intent(in) :: iregion_code,nspec
  integer,intent(in) :: ipass

  ! arrays with the mesh in double precision
  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

  integer, dimension(nspec) :: idoubling

  ! this for non blocking MPI
  logical, dimension(nspec) :: is_on_a_slice_edge

  integer :: NEX_PER_PROC_XI,NEX_PER_PROC_ETA

  integer :: offset_proc_xi,offset_proc_eta

  ! local parameters
  integer :: ispec,nspec_tiso
  ! parameters needed to store the radii of the grid points in the spherically symmetric Earth
  double precision :: rmin,rmax
  integer :: ner_without_doubling,ilayer,ilayer_loop
  ! timing
  double precision, external :: wtime
  double precision :: time_start,tCPU
  integer,dimension(8) :: tval

  ! initializes flags for transverse isotropic elements
  ispec_is_tiso(:) = .false.

  ! get MPI starting time
  time_start = wtime()

  ! loop on all the layers in this region of the mesh
  ispec = 0 ! counts all the elements in this region of the mesh
  do ilayer_loop = ifirst_region,ilast_region

    ilayer = perm_layer(ilayer_loop)

    ! user output
    if(myrank == 0 ) then
      write(IMAIN,*) '  creating layer ',ilayer_loop-ifirst_region+1, &
                                   'out of ',ilast_region-ifirst_region+1
    endif

    ! determine the radii that define the shell
    rmin = rmins(ilayer)
    rmax = rmaxs(ilayer)

    if(iregion_code == IREGION_CRUST_MANTLE .and. ilayer_loop==3) then
      FIRST_ELT_NON_ANISO = ispec+1
    endif
    if(iregion_code == IREGION_CRUST_MANTLE &
      .and. ilayer_loop==(ilast_region-nb_layer_above_aniso+1)) then
      FIRST_ELT_ABOVE_ANISO = ispec+1
    endif

    ner_without_doubling = ner(ilayer)

    ! if there is a doubling at the top of this region, we implement it in the last two layers of elements
    ! and therefore we suppress two layers of regular elements here
    USE_ONE_LAYER_SB = .false.
    if(this_region_has_a_doubling(ilayer)) then
      if (ner(ilayer) == 1) then
        ner_without_doubling = ner_without_doubling - 1
        USE_ONE_LAYER_SB = .true.
      else
        ner_without_doubling = ner_without_doubling - 2
        USE_ONE_LAYER_SB = .false.
      endif
    endif

    ! regular mesh elements
    call create_regular_elements(myrank,ilayer,ichunk,ispec,ipass, &
                    ifirst_region,ilast_region,iregion_code, &
                    nspec,NCHUNKS,NUMBER_OF_MESH_LAYERS, &
                    NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                    ner_without_doubling,ner,ratio_sampling_array,r_top,r_bottom, &
                    xstore,ystore,zstore, &
                    iaddx,iaddy,iaddz,xigll,yigll,zigll, &
                    shape3D,dershape2D_bottom, &
                    INCLUDE_CENTRAL_CUBE,ABSORBING_CONDITIONS, &
                    RICB,RCMB,R670,RMOHO,RMOHO_FICTITIOUS_IN_MESHER,RTOPDDOUBLEPRIME, &
                    R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
                    rmin,rmax,r_moho,r_400,r_670, &
                    rhostore,dvpstore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                    nspec_ani,c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                    nspec_actually,xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,&
                    gammaxstore,gammaystore,gammazstore,&
                    nspec_stacey,rho_vp,rho_vs,iboun,iMPIcut_xi,iMPIcut_eta, &
                    ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
                    nspec_att,Qmu_store,tau_e_store,tau_s,T_c_source, &
                    size(tau_e_store,2),size(tau_e_store,3),size(tau_e_store,4), &
                    rotation_matrix,idoubling,doubling_index,USE_ONE_LAYER_SB, &
                    stretch_tab,ACTUALLY_STORE_ARRAYS, &
                    NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,nex_eta_moho, &
                    ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot,ibelm_670_top,ibelm_670_bot, &
                    normal_moho,normal_400,normal_670,jacobian2D_moho,jacobian2D_400,jacobian2D_670, &
                    ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,&
                    ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot,&
                    ispec_is_tiso)


    ! mesh doubling elements
    if( this_region_has_a_doubling(ilayer) ) &
      call create_doubling_elements(myrank,ilayer,ichunk,ispec,ipass, &
                    ifirst_region,ilast_region,iregion_code, &
                    nspec,NCHUNKS,NUMBER_OF_MESH_LAYERS, &
                    NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                    ner,ratio_sampling_array,r_top,r_bottom, &
                    xstore,ystore,zstore,xigll,yigll,zigll, &
                    shape3D,dershape2D_bottom, &
                    INCLUDE_CENTRAL_CUBE,ABSORBING_CONDITIONS, &
                    RICB,RCMB,R670,RMOHO,RMOHO_FICTITIOUS_IN_MESHER,RTOPDDOUBLEPRIME, &
                    R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
                    rmin,rmax,r_moho,r_400,r_670, &
                    rhostore,dvpstore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                    nspec_ani,c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                    nspec_actually,xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,&
                    gammaxstore,gammaystore,gammazstore,&
                    nspec_stacey,rho_vp,rho_vs,iboun,iMPIcut_xi,iMPIcut_eta, &
                    ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
                    nspec_att,Qmu_store,tau_e_store,tau_s,T_c_source, &
                    size(tau_e_store,2),size(tau_e_store,3),size(tau_e_store,4), &
                    rotation_matrix,idoubling,doubling_index,USE_ONE_LAYER_SB,ACTUALLY_STORE_ARRAYS, &
                    NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,nex_eta_moho, &
                    ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot,ibelm_670_top,ibelm_670_bot, &
                    normal_moho,normal_400,normal_670,jacobian2D_moho,jacobian2D_400,jacobian2D_670, &
                    ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,&
                    ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot, &
                    CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,offset_proc_xi,offset_proc_eta, &
                    ispec_is_tiso)

    ! user output
    if(myrank == 0 ) then
      ! time estimate
      tCPU = wtime() - time_start
      ! remaining
      tCPU = (1.0-(ilayer_loop-ifirst_region+1.0)/(ilast_region-ifirst_region+1.0)) &
                /(ilayer_loop-ifirst_region+1.0)/(ilast_region-ifirst_region+1.0)*tCPU*10.0
      write(IMAIN,*) "    ",(ilayer_loop-ifirst_region+1.0)/(ilast_region-ifirst_region+1.0) * 100.0," %", &
                    " time remaining:", tCPU,"s"
      call date_and_time(VALUES=tval)
      write(IMAIN,*) "    ",tval(5),"h",tval(6),"min",tval(7),".",tval(8),"sec"
    endif

  enddo !ilayer_loop

  deallocate(stretch_tab)
  deallocate(perm_layer)
  deallocate(jacobian2D_moho,jacobian2D_400,jacobian2D_670)
      
  if(myrank == 0 ) write(IMAIN,*)

  ! define central cube in inner core
  if(INCLUDE_CENTRAL_CUBE .and. iregion_code == IREGION_INNER_CORE) then
    ! user output
    if(myrank == 0 ) write(IMAIN,*) '  creating central cube'

    call create_central_cube(myrank,ichunk,ispec,iaddx,iaddy,iaddz,ipass, &
                        nspec,NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,R_CENTRAL_CUBE, &
                        iproc_xi,iproc_eta,NPROC_XI,NPROC_ETA,ratio_divide_central_cube, &
                        iMPIcut_xi,iMPIcut_eta,iboun, &
                        idoubling,iregion_code,xstore,ystore,zstore, &
                        RICB,RCMB,R670,RMOHO,RMOHO_FICTITIOUS_IN_MESHER,RTOPDDOUBLEPRIME,&
                        R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
                        shape3D,rmin,rmax,rhostore,dvpstore,&
                        kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                        xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,&
                        gammaxstore,gammaystore,gammazstore,nspec_actually, &
                        c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                        c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                        c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                        nspec_ani,nspec_stacey,nspec_att,Qmu_store,tau_e_store,tau_s,T_c_source,&
                        size(tau_e_store,2),size(tau_e_store,3),size(tau_e_store,4), &
                        rho_vp,rho_vs,ABSORBING_CONDITIONS,ACTUALLY_STORE_ARRAYS,xigll,yigll,zigll, &
                        ispec_is_tiso)
  endif

  ! check total number of spectral elements created
  if(ispec /= nspec) call exit_MPI(myrank,'ispec should equal nspec')

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
  if( ipass == 2 ) then
    ! count number of anisotropic elements in current region
    ! should be zero in all the regions except in the mantle
    nspec_tiso = count(ispec_is_tiso(:))

    ! checks number of anisotropic elements found in the mantle
    if(iregion_code /= IREGION_CRUST_MANTLE .and. nspec_tiso /= 0 ) &
      call exit_MPI(myrank,'found anisotropic elements outside of the mantle')
    if( TRANSVERSE_ISOTROPY ) then
      if(iregion_code == IREGION_CRUST_MANTLE .and. nspec_tiso == 0) &
        call exit_MPI(myrank,'found no anisotropic elements in the mantle')
    endif
  endif

  end subroutine crm_create_elements

!
!-------------------------------------------------------------------------------------------------
!

  subroutine crm_setup_indexing(ibool,xstore,ystore,zstore, &
                                nspec,nglob_theor,npointot)

! creates global indexing array ibool

  use meshfem3d_par,only: &
    myrank,NGLLX,NGLLY,NGLLZ

  use create_regions_mesh_par2
  implicit none

  ! number of spectral elements in each block
  integer,intent(in) :: nspec,npointot,nglob_theor

  ! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  ! local parameters
  ! variables for creating array ibool
  double precision, dimension(:), allocatable :: xp,yp,zp
  integer, dimension(:), allocatable :: locval
  logical, dimension(:), allocatable :: ifseg

  integer :: nglob
  integer :: ieoff,ilocnum,ier
  integer :: i,j,k,ispec
  character(len=150) :: errmsg

  ! allocate memory for arrays
  allocate(locval(npointot), &
          ifseg(npointot), &
          xp(npointot), &
          yp(npointot), &
          zp(npointot),stat=ier)
  if(ier /= 0) stop 'error in allocate 20'

  locval = 0
  ifseg = .false.
  xp = 0.d0
  yp = 0.d0
  zp = 0.d0

  ! we need to create a copy of the x, y and z arrays because sorting in get_global will swap
  ! these arrays and therefore destroy them
  do ispec=1,nspec
    ieoff = NGLLX * NGLLY * NGLLZ * (ispec-1)
    ilocnum = 0
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          ilocnum = ilocnum + 1
          xp(ilocnum+ieoff) = xstore(i,j,k,ispec)
          yp(ilocnum+ieoff) = ystore(i,j,k,ispec)
          zp(ilocnum+ieoff) = zstore(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo

  call get_global(nspec,xp,yp,zp,ibool,locval,ifseg,nglob,npointot)

  deallocate(xp,yp,zp)
  deallocate(locval,ifseg)

  ! check that number of points found equals theoretical value
  if(nglob /= nglob_theor) then
    write(errmsg,*) 'incorrect total number of points found: myrank,nglob,nglob_theor = ',&
      myrank,nglob,nglob_theor
    call exit_MPI(myrank,errmsg)
  endif
  if(minval(ibool) /= 1 .or. maxval(ibool) /= nglob_theor) call exit_MPI(myrank,'incorrect global numbering')

  ! creates a new indirect addressing to reduce cache misses in memory access in the solver
  ! this is *critical* to improve performance in the solver
  call get_global_indirect_addressing(nspec,nglob_theor,ibool)

  ! checks again
  if(minval(ibool) /= 1 .or. maxval(ibool) /= nglob_theor) call exit_MPI(myrank,'incorrect global numbering after sorting')

  end subroutine crm_setup_indexing

!
!-------------------------------------------------------------------------------------------------
!

  subroutine crm_setup_mpi_buffers(npointot,nspec,ibool,idoubling, &
                                  xstore,ystore,zstore,iregion_code, &
                                  npoin2D_xi,npoin2D_eta)

! sets up MPI cutplane arrays

  use meshfem3d_par,only: &
    myrank,NGLLX,NGLLY,NGLLZ, &
    NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER, &
    NSPEC2D_XI_FACE,NSPEC2D_ETA_FACE

  use create_regions_mesh_par2
  implicit none

  ! number of spectral elements in each block
  integer,intent(in) :: nspec,npointot

  ! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: idoubling

  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  integer,intent(in) :: iregion_code
  integer :: npoin2D_xi,npoin2D_eta

  ! local parameters
  logical, dimension(:), allocatable :: mask_ibool
  integer :: ier

  ! arrays mask_ibool(npointot) used to save memory
  ! allocate memory for arrays
  allocate(mask_ibool(npointot), &
          stat=ier)
  if(ier /= 0) stop 'error in allocate 20b'

  call get_MPI_cutplanes_xi(myrank,prname,nspec,iMPIcut_xi,ibool, &
                  xstore,ystore,zstore,mask_ibool,npointot, &
                  NSPEC2D_ETA_FACE,iregion_code,npoin2D_xi)

  call get_MPI_cutplanes_eta(myrank,prname,nspec,iMPIcut_eta,ibool, &
                  xstore,ystore,zstore,mask_ibool,npointot, &
                  NSPEC2D_XI_FACE,iregion_code,npoin2D_eta)

  call get_MPI_1D_buffers(myrank,prname,nspec,iMPIcut_xi,iMPIcut_eta,ibool,idoubling, &
                  xstore,ystore,zstore,mask_ibool,npointot, &
                  NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER,iregion_code)

  deallocate(mask_ibool)

  end subroutine crm_setup_mpi_buffers

!
!-------------------------------------------------------------------------------------------------
!

  subroutine crm_save_boundary_mesh()

! creates global indexing array ibool

  use meshfem3d_par,only: &
    myrank

  use meshfem3D_models_par,only: &
    SAVE_BOUNDARY_MESH,HONOR_1D_SPHERICAL_MOHO,SUPPRESS_CRUSTAL_MESH

  use create_regions_mesh_par2
  implicit none

  ! local parameters
  integer :: ier

  ! first check the number of surface elements are the same for Moho, 400, 670
  if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
    if (ispec2D_moho_top /= NSPEC2D_MOHO .or. ispec2D_moho_bot /= NSPEC2D_MOHO) &
           call exit_mpi(myrank, 'Not the same number of Moho surface elements')
  endif
  if (ispec2D_400_top /= NSPEC2D_400 .or. ispec2D_400_bot /= NSPEC2D_400) &
           call exit_mpi(myrank,'Not the same number of 400 surface elements')
  if (ispec2D_670_top /= NSPEC2D_670 .or. ispec2D_670_bot /= NSPEC2D_670) &
           call exit_mpi(myrank,'Not the same number of 670 surface elements')

  ! writing surface topology databases
  open(unit=27,file=prname(1:len_trim(prname))//'boundary_disc.bin', &
       status='unknown',form='unformatted',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error opening boundary_disc.bin file')

  write(27) NSPEC2D_MOHO, NSPEC2D_400, NSPEC2D_670
  write(27) ibelm_moho_top
  write(27) ibelm_moho_bot
  write(27) ibelm_400_top
  write(27) ibelm_400_bot
  write(27) ibelm_670_top
  write(27) ibelm_670_bot
  write(27) normal_moho
  write(27) normal_400
  write(27) normal_670
  close(27)

  end subroutine crm_save_boundary_mesh

!
!-------------------------------------------------------------------------------------------------
!

  subroutine crm_save_mesh_files(nspec,npointot,iregion_code,ibool,idoubling, &
                                  xstore,ystore,zstore)

  use meshfem3d_par,only: &
    myrank,NGLLX,NGLLY,NGLLZ, &
    RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
    RMIDDLE_CRUST,ROCEAN

  use meshfem3D_models_par,only: &
    ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
    nspl,rspl,espl,espl2

  use create_regions_mesh_par2
  implicit none

  ! number of spectral elements in each block
  integer,intent(in) :: nspec,npointot,iregion_code

  ! arrays with the mesh
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: idoubling

  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  ! local parameters
  ! arrays used for AVS or DX files
  integer, dimension(:), allocatable :: num_ibool_AVS_DX
  logical, dimension(:), allocatable :: mask_ibool
  integer :: ier

  ! arrays num_ibool_AVS_DX and mask_ibool used to save memory
  ! allocate memory for arrays
  allocate(num_ibool_AVS_DX(npointot), &
          mask_ibool(npointot), &
          stat=ier)
  if(ier /= 0) stop 'error in allocate 21'

  call write_AVS_DX_global_data(myrank,prname,nspec,ibool,idoubling,xstore,ystore,zstore, &
                                num_ibool_AVS_DX,mask_ibool,npointot)

  call write_AVS_DX_global_faces_data(myrank,prname,nspec,iMPIcut_xi,iMPIcut_eta,ibool, &
          idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool,npointot, &
          rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
          ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
          RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
          RMIDDLE_CRUST,ROCEAN,iregion_code)

  call write_AVS_DX_global_chunks_data(myrank,prname,nspec,iboun,ibool, &
          idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool,npointot, &
          rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
          ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
          RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
          RMIDDLE_CRUST,ROCEAN,iregion_code)

  call write_AVS_DX_surface_data(myrank,prname,nspec,iboun,ibool, &
          idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool,npointot, &
          rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
          ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
          RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
          RMIDDLE_CRUST,ROCEAN,iregion_code)

  ! Output material information for all GLL points
  ! Can be use to check the mesh
  !    call write_AVS_DX_global_data_gll(prname,nspec,xstore,ystore,zstore,&
  !                rhostore,kappavstore,muvstore,Qmu_store,ATTENUATION)
  deallocate(num_ibool_AVS_DX,mask_ibool)

  end subroutine crm_save_mesh_files
