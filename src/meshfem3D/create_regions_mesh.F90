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

  subroutine create_regions_mesh(iregion_code,npointot, &
                                 NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                                 NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                                 NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                                 offset_proc_xi,offset_proc_eta, &
                                 ipass)

! creates the different regions of the mesh

!
! ****************************************************************************************************
! IMPORTANT: this routine must *NOT* use flag SIMULATION_TYPE (nor SAVE_FORWARD), i.e. none of the parameters it computes
! should depend on SIMULATION_TYPE, because most users do not recompile the code nor rerun the mesher
! when switching from SIMULATION_TYPE == 1 to SIMULATION_TYPE == 3 and thus the header file created
! by this routine would become wrong in the case of a run with SIMULATION_TYPE == 3 if the code
! was compiled with SIMULATION_TYPE == 1
! ****************************************************************************************************
!

  use meshfem3D_par, only: &
    nspec,nglob, &
    ibool,idoubling,xstore,ystore,zstore, &
    xstore_glob,ystore_glob,zstore_glob, &
    IMAIN,volume_total,Earth_mass_total,Earth_center_of_mass_x_total,Earth_center_of_mass_y_total,Earth_center_of_mass_z_total, &
    myrank,LOCAL_PATH,IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE, &
    IFLAG_IN_FICTITIOUS_CUBE, &
    NCHUNKS,SAVE_MESH_FILES,ABSORBING_CONDITIONS, &
    R_CENTRAL_CUBE,RICB,RCMB, &
    MAX_NUMBER_OF_MESH_LAYERS,MAX_NUM_REGIONS,NB_SQUARE_CORNERS, &
    NGLOB1D_RADIAL_CORNER, &
    NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
    ADIOS_FOR_ARRAYS_SOLVER, &
    ROTATION,EXACT_MASS_MATRIX_FOR_ROTATION,GRAVITY_INTEGRALS

  use meshfem3D_models_par, only: &
    SAVE_BOUNDARY_MESH,SUPPRESS_CRUSTAL_MESH,REGIONAL_MOHO_MESH, &
    OCEANS
#ifdef CEM
  use meshfem3D_models_par, only: CEM_REQUEST
#endif

  use MPI_interfaces_par, only: &
    NGLOB1D_RADIAL_MAX,iboolcorner,iboolfaces, &
    iboolleft_xi,iboolright_xi,iboolleft_eta,iboolright_eta, &
    ibool1D_leftxi_lefteta,ibool1D_rightxi_lefteta, &
    ibool1D_leftxi_righteta,ibool1D_rightxi_righteta, &
    xyz1D_leftxi_lefteta,xyz1D_rightxi_lefteta, &
    xyz1D_leftxi_righteta,xyz1D_rightxi_righteta

  use regions_mesh_par

  use regions_mesh_par2

  implicit none

  ! code for the four regions of the mesh
  integer,intent(in) :: iregion_code

  ! correct number of spectral elements in each block depending on chunk type
  integer,intent(in) :: npointot

  integer :: NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer :: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX
  integer :: NSPEC2D_BOTTOM,NSPEC2D_TOP

  integer :: offset_proc_xi,offset_proc_eta

  ! now perform two passes in this part to be able to save memory
  integer,intent(in) :: ipass

  ! local parameters
  integer :: ier

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    select case (ipass)
    case (1)
      write(IMAIN,*) 'first pass'
    case (2)
      write(IMAIN,*) 'second pass'
    case default
      call exit_MPI(myrank,'Error ipass value in create_regions_mesh')
    end select
    call flush_IMAIN()
  endif

  ! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,iregion_code,LOCAL_PATH)

  ! initializes arrays
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...allocating arrays '
    call flush_IMAIN()
  endif
  call crm_allocate_arrays(iregion_code,ipass, &
                           NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                           NSPEC2D_BOTTOM,NSPEC2D_TOP)


  ! initialize number of layers
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...setting up layers '
    call flush_IMAIN()
  endif
  call crm_setup_layers(iregion_code,ipass,NEX_PER_PROC_ETA)

  !  creates mesh elements
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  ...creating mesh elements '
    call flush_IMAIN()
  endif
  call create_regions_elements(iregion_code,ipass, &
                               NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                               offset_proc_xi,offset_proc_eta)

  select case (ipass)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case (1) !!!!!!!!!!! first pass of the mesher
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! only create global addressing and the MPI buffers in the first pass

    ! creates ibool index array for projection from local to global points
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...creating global addressing'
      call flush_IMAIN()
    endif
    call crm_setup_indexing(npointot)

    ! create MPI buffers
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...creating MPI buffers'
      call flush_IMAIN()
    endif
    call crm_setup_mpi_buffers(npointot,iregion_code)


    ! sets up Stacey absorbing boundary indices (nimin,nimax,..)
    if (NCHUNKS /= 6) call get_absorb(prname,iregion_code, &
                                      iboun,nimin,nimax, &
                                      njmin,njmax, nkmin_xi,nkmin_eta, NSPEC2DMAX_XMIN_XMAX, &
                                      NSPEC2DMAX_YMIN_YMAX, NSPEC2D_BOTTOM)

    ! Only go into here if we're requesting xyz files for CEM
#ifdef CEM
    if (CEM_REQUEST) then

      call build_global_coordinates(iregion_code)
      call write_cem_request(iregion_code)
      call synchronize_all()

      deallocate(ibool1D_leftxi_lefteta,ibool1D_rightxi_lefteta, &
                 ibool1D_leftxi_righteta,ibool1D_rightxi_righteta, &
                 xyz1D_leftxi_lefteta,xyz1D_rightxi_lefteta, &
                 xyz1D_leftxi_righteta,xyz1D_rightxi_righteta,iboolleft_xi, &
                 iboolright_xi,iboolleft_eta,iboolright_eta,nimin,nimax, &
                 njmin, njmax,nkmin_xi,nkmin_eta,iboolfaces,iboolcorner)

    endif
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  case (2) !!!!!!!!!!! second pass of the mesher
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! create mass matrix and save all the final arrays in the second pass

    ! allocates arrays for global mesh coordinates (used for MPI interfaces setup and saving mesh outputs)
    call crm_fill_global_meshes()

    ! checks time step size since now all elements have material properties assigned
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...checking mesh resolution and time step'
      call flush_IMAIN()
    endif
    call check_mesh_resolution(iregion_code,xstore,ystore,zstore, &
                               kappavstore,kappahstore,muvstore,muhstore,rhostore)

    ! precomputes Jacobian for 2D absorbing boundary surfaces
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...precomputing Jacobian'
      call flush_IMAIN()
    endif
    call get_jacobian_boundaries(iboun,xstore,ystore,zstore, &
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
                                 NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                                 xigll,yigll,zigll)

!! DK DK for gravity integrals
    ! creation of the top observation surface if region is the crust_mantle
    if (GRAVITY_INTEGRALS) then
      if (iregion_code == IREGION_CRUST_MANTLE) call gravity_observation_surface()
    endif

    ! create chunk buffers if more than one chunk
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...creating chunk buffers'
      call flush_IMAIN()
    endif
    call create_chunk_buffers(iregion_code,ibool,idoubling, &
                              xstore,ystore,zstore, &
                              NGLOB1D_RADIAL_CORNER,NGLOB1D_RADIAL_MAX, &
                              NGLOB2DMAX_XMIN_XMAX(iregion_code),NGLOB2DMAX_YMIN_YMAX(iregion_code))

    ! only deallocates after second pass
    deallocate(ibool1D_leftxi_lefteta,ibool1D_rightxi_lefteta, &
               ibool1D_leftxi_righteta,ibool1D_rightxi_righteta, &
               xyz1D_leftxi_lefteta,xyz1D_rightxi_lefteta, &
               xyz1D_leftxi_righteta,xyz1D_rightxi_righteta)

    ! setup MPI communication interfaces
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...preparing MPI interfaces'
      call flush_IMAIN()
    endif
    ! creates MPI interface arrays
    call create_MPI_interfaces(iregion_code)

    ! sets up MPI interface arrays
    call setup_MPI_interfaces(iregion_code)

    ! only deallocates after second pass
    deallocate(iboolleft_xi,iboolright_xi,iboolleft_eta,iboolright_eta)
    deallocate(iboolfaces)
    deallocate(iboolcorner)

    ! sets up inner/outer element arrays
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...element inner/outer separation '
      call flush_IMAIN()
    endif
    call setup_inner_outer(iregion_code)

    ! sets up mesh coloring
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...element mesh coloring '
      call flush_IMAIN()
    endif
    call setup_color_perm(iregion_code)

    !uncomment: adds model smoothing for point profile models
    !    if (THREE_D_MODEL == THREE_D_MODEL_PPM) then
    !     call smooth_model(nproc_xi,nproc_eta, &
    !        rho_vp,rho_vs,nspec_stacey, &
    !        iregion_code,xixstore,xiystore,xizstore, &
    !        etaxstore,etaystore,etazstore, &
    !        gammaxstore,gammaystore,gammazstore, &
    !        xstore,ystore,zstore,rhostore,dvpstore, &
    !        kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
    !        nspec,HETEROGEN_3D_MANTLE, &
    !        NEX_XI,NCHUNKS,ABSORBING_CONDITIONS )

    ! creates mass matrix
    call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  ...creating mass matrix'
      call flush_IMAIN()
    endif

    ! allocates mass matrices in this slice (will be fully assembled in the solver)
    !
    ! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
    ! on Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
    ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
    !
    ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
    ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be obsolete
    if (NCHUNKS /= 6 .and. ABSORBING_CONDITIONS) then
      select case (iregion_code)
      case (IREGION_CRUST_MANTLE)
        nglob_xy = nglob
      case (IREGION_INNER_CORE, IREGION_OUTER_CORE)
        nglob_xy = 1
      end select
    else
       nglob_xy = 1
    endif

    if (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
      select case (iregion_code)
      case (IREGION_CRUST_MANTLE,IREGION_INNER_CORE)
         nglob_xy = nglob
      case (IREGION_OUTER_CORE)
         nglob_xy = 1
      end select
    endif

    allocate(rmassx(nglob_xy), &
             rmassy(nglob_xy), &
             stat=ier)
    if (ier /= 0) stop 'Error in allocate 21'

    allocate(b_rmassx(nglob_xy), &
             b_rmassy(nglob_xy),stat=ier)
    if (ier /= 0) stop 'Error in allocate b_21'

    allocate(rmassz(nglob),stat=ier)
    if (ier /= 0) stop 'Error in allocate 22'

    ! allocates ocean load mass matrix as well if oceans
    if (OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) then
      nglob_oceans = nglob
    else
      ! allocate dummy array if no oceans
      nglob_oceans = 1
    endif
    allocate(rmass_ocean_load(nglob_oceans),stat=ier)
    if (ier /= 0) stop 'Error in allocate 22'

    ! creating mass matrices in this slice (will be fully assembled in the solver)
    ! note: for Stacey boundaries, needs indexing nimin,.. filled in the first pass
    call create_mass_matrices(idoubling,ibool, &
                              iregion_code,xstore,ystore,zstore, &
                              NSPEC2D_TOP,NSPEC2D_BOTTOM)

    ! Stacey
    deallocate(nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta)

    ! save the binary files
    call synchronize_all()

!! DK DK for gravity integrals
    ! gravity integrals computation won't need mesh nor solver runs
    if (.not. GRAVITY_INTEGRALS) then
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) '  ...saving binary files'
        call flush_IMAIN()
      endif
      ! saves mesh and model parameters
      if (ADIOS_FOR_ARRAYS_SOLVER) then
        if (myrank == 0) write(IMAIN,*) '    in ADIOS file format'
        call save_arrays_solver_adios(idoubling,ibool, &
                                      iregion_code,xstore,ystore,zstore, &
                                      NSPEC2DMAX_XMIN_XMAX, NSPEC2DMAX_YMIN_YMAX, &
                                      NSPEC2D_TOP,NSPEC2D_BOTTOM)
      else
        call save_arrays_solver(idoubling,ibool, &
                                iregion_code,xstore,ystore,zstore, &
                                NSPEC2D_TOP,NSPEC2D_BOTTOM)
      endif

      ! saves MPI interface info
      call save_arrays_solver_MPI(iregion_code)

      ! boundary mesh for MOHO, 400 and 670 discontinuities
      if (SAVE_BOUNDARY_MESH .and. iregion_code == IREGION_CRUST_MANTLE) then
        ! user output
        call synchronize_all()
        if (myrank == 0) then
          write(IMAIN,*)
          write(IMAIN,*) '  ...saving boundary mesh files'
          call flush_IMAIN()
        endif
        ! saves boundary file
        if (ADIOS_FOR_ARRAYS_SOLVER) then
          if (myrank == 0) write(IMAIN,*) '    in ADIOS file format'
          call save_arrays_boundary_adios()
        else
          call save_arrays_boundary()
        endif
      endif

      ! create AVS or DX mesh data for the slices
      if (SAVE_MESH_FILES) then
        ! user output
        call synchronize_all()
        if (myrank == 0) then
          write(IMAIN,*)
          write(IMAIN,*) '  ...saving AVS or DX mesh files'
          call flush_IMAIN()
        endif
        call write_AVS_DX_output(npointot,iregion_code)
      endif

    endif ! .not. GRAVITY_INTEGRALS

    ! frees memory
    deallocate(rmassx,rmassy,rmassz)
    deallocate(b_rmassx,b_rmassy)
    deallocate(rmass_ocean_load)
    ! frees allocated mesh memory
    deallocate(xstore_glob,ystore_glob,zstore_glob)
    ! frees MPI arrays memory
    call crm_free_MPI_arrays(iregion_code)

    ! compute volume, bottom and top area of that part of the slice, and then the total
    call compute_volumes_and_areas(NCHUNKS,iregion_code,nspec,wxgll,wygll,wzgll, &
                                   xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
                                   NSPEC2D_BOTTOM,jacobian2D_bottom,NSPEC2D_TOP,jacobian2D_top,idoubling, &
                                   volume_total,RCMB,RICB,R_CENTRAL_CUBE)

    ! compute Earth mass of that part of the slice, and then total Earth mass
    call compute_Earth_mass(Earth_mass_total, &
                            Earth_center_of_mass_x_total,Earth_center_of_mass_y_total,Earth_center_of_mass_z_total, &
                            nspec,wxgll,wygll,wzgll,xstore,ystore,zstore,xixstore,xiystore,xizstore, &
                            etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,rhostore,idoubling)

!! DK DK for gravity integrals
    if (GRAVITY_INTEGRALS) then
      call synchronize_all()
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) '  ...computing gravity integrals'
        write(IMAIN,*)
        call flush_IMAIN()
      endif
      ! compute gravity integrals of that part of the slice, and then total integrals for the whole Earth
      call gravity_compute_integrals(iregion_code,nspec,wxgll,wygll,wzgll,xstore,ystore,zstore, &
                                     xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
                                     rhostore,idoubling)

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

  end subroutine create_regions_mesh

!
!-------------------------------------------------------------------------------------------------
!

  subroutine crm_allocate_arrays(iregion_code,ipass, &
                                 NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                                 NSPEC2D_BOTTOM,NSPEC2D_TOP)

  use constants

  use meshfem3D_par, only: &
    nspec, &
    NCHUNKS,NUMCORNERS_SHARED,NUMFACES_SHARED, &
    NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
    NGLOB1D_RADIAL,NGLOB1D_RADIAL_CORNER, &
    ATT1,ATT2,ATT3

  use meshfem3D_models_par, only: &
    ATTENUATION,ANISOTROPIC_INNER_CORE,ANISOTROPIC_3D_MANTLE, &
    SAVE_BOUNDARY_MESH,AM_V

  use regions_mesh_par2

  use MPI_interfaces_par

  implicit none

  integer,intent(in) :: iregion_code
  integer,intent(in) :: ipass

  integer,intent(in) :: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX
  integer,intent(in) :: NSPEC2D_BOTTOM,NSPEC2D_TOP

  ! local parameters
  integer :: ier

  ! adios needs properly initialized arrays, otherwise its intrinsic check procedures will cause undefined operations

  ! new attenuation definition on all GLL points
  ! attenuation
  if (ATTENUATION) then
    T_c_source = AM_V%QT_c_source
    tau_s(:)   = AM_V%Qtau_s(:)
    nspec_att = nspec
  else
    ! allocates dummy size arrays
    nspec_att = 1
  endif
  allocate(Qmu_store(ATT1,ATT2,ATT3,nspec_att), &
           tau_e_store(ATT1,ATT2,ATT3,N_SLS,nspec_att),stat=ier)
  if (ier /= 0) stop 'Error in allocate 1'

  Qmu_store(:,:,:,:) = 0.0; tau_e_store(:,:,:,:,:) = 0.0

  ! array with model density
  allocate(rhostore(NGLLX,NGLLY,NGLLZ,nspec), &
           dvpstore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) stop 'Error in allocate 6'

  rhostore(:,:,:,:) = 0.0; dvpstore(:,:,:,:) = 0.0

  ! for anisotropy
  allocate(kappavstore(NGLLX,NGLLY,NGLLZ,nspec), &
           muvstore(NGLLX,NGLLY,NGLLZ,nspec), &
           kappahstore(NGLLX,NGLLY,NGLLZ,nspec), &
           muhstore(NGLLX,NGLLY,NGLLZ,nspec), &
           eta_anisostore(NGLLX,NGLLY,NGLLZ,nspec), &
           ispec_is_tiso(nspec),stat=ier)
  if (ier /= 0) stop 'Error in allocate 7'

  kappavstore(:,:,:,:) = 0.0; kappahstore(:,:,:,:) = 0.0
  muvstore(:,:,:,:) = 0.0; muhstore(:,:,:,:) = 0.0
  eta_anisostore(:,:,:,:) = 0.0
  ispec_is_tiso(:) = .false.

  ! Stacey absorbing boundaries
  if (NCHUNKS /= 6) then
    nspec_stacey = nspec
  else
    nspec_stacey = 1
  endif
  allocate(rho_vp(NGLLX,NGLLY,NGLLZ,nspec_stacey), &
           rho_vs(NGLLX,NGLLY,NGLLZ,nspec_stacey),stat=ier)
  if (ier /= 0) stop 'Error in allocate 8'

  rho_vp(:,:,:,:) = 0.0; rho_vs(:,:,:,:) = 0.0

  ! anisotropy
  if ((ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) .or. &
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
  if (ier /= 0) stop 'Error in allocate 9'

  c11store(:,:,:,:) = 0.0; c12store(:,:,:,:) = 0.0; c13store(:,:,:,:) = 0.0
  c14store(:,:,:,:) = 0.0; c15store(:,:,:,:) = 0.0; c16store(:,:,:,:) = 0.0
  c22store(:,:,:,:) = 0.0; c23store(:,:,:,:) = 0.0; c24store(:,:,:,:) = 0.0
  c25store(:,:,:,:) = 0.0; c26store(:,:,:,:) = 0.0; c33store(:,:,:,:) = 0.0
  c34store(:,:,:,:) = 0.0; c35store(:,:,:,:) = 0.0; c36store(:,:,:,:) = 0.0
  c44store(:,:,:,:) = 0.0; c45store(:,:,:,:) = 0.0; c46store(:,:,:,:) = 0.0
  c55store(:,:,:,:) = 0.0; c56store(:,:,:,:) = 0.0; c66store(:,:,:,:) = 0.0

  ! boundary locator
  allocate(iboun(6,nspec),stat=ier)
  if (ier /= 0) stop 'Error in allocate 10'

  ! boundary parameters locator
  allocate(ibelm_xmin(NSPEC2DMAX_XMIN_XMAX), &
           ibelm_xmax(NSPEC2DMAX_XMIN_XMAX), &
           ibelm_ymin(NSPEC2DMAX_YMIN_YMAX), &
           ibelm_ymax(NSPEC2DMAX_YMIN_YMAX), &
           ibelm_bottom(NSPEC2D_BOTTOM), &
           ibelm_top(NSPEC2D_TOP),stat=ier)
  if (ier /= 0) stop 'Error in allocate 11'

  ibelm_xmin(:) = 0; ibelm_xmax(:) = 0
  ibelm_ymin(:) = 0; ibelm_ymax(:) = 0
  ibelm_bottom(:) = 0; ibelm_top(:) = 0

  ! 2-D Jacobians and normals
  allocate(jacobian2D_xmin(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX), &
           jacobian2D_xmax(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX), &
           jacobian2D_ymin(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX), &
           jacobian2D_ymax(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX), &
           jacobian2D_bottom(NGLLX,NGLLY,NSPEC2D_BOTTOM), &
           jacobian2D_top(NGLLX,NGLLY,NSPEC2D_TOP),stat=ier)
  if (ier /= 0) stop 'Error in allocate 12'

  jacobian2D_xmin(:,:,:) = 0.0; jacobian2D_xmax(:,:,:) = 0.0
  jacobian2D_ymin(:,:,:) = 0.0; jacobian2D_ymax(:,:,:) = 0.0
  jacobian2D_bottom(:,:,:) = 0.0; jacobian2D_top(:,:,:) = 0.0

  allocate(normal_xmin(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX), &
           normal_xmax(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX), &
           normal_ymin(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX), &
           normal_ymax(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX), &
           normal_bottom(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM), &
           normal_top(NDIM,NGLLX,NGLLY,NSPEC2D_TOP),stat=ier)
  if (ier /= 0) stop 'Error in allocate 13'

  normal_xmin(:,:,:,:) = 0.0; normal_xmax(:,:,:,:) = 0.0
  normal_ymin(:,:,:,:) = 0.0; normal_ymax(:,:,:,:) = 0.0
  normal_bottom(:,:,:,:) = 0.0; normal_top(:,:,:,:) = 0.0

  ! Stacey
  if (ipass == 1) then
    if (NCHUNKS /= 6) then
      allocate(nimin(2,NSPEC2DMAX_YMIN_YMAX), &
               nimax(2,NSPEC2DMAX_YMIN_YMAX), &
               njmin(2,NSPEC2DMAX_XMIN_XMAX), &
               njmax(2,NSPEC2DMAX_XMIN_XMAX), &
               nkmin_xi(2,NSPEC2DMAX_XMIN_XMAX), &
               nkmin_eta(2,NSPEC2DMAX_YMIN_YMAX),stat=ier)
      if (ier /= 0) stop 'Error in allocate 14'
    else
      allocate(nimin(1,1), &
               nimax(1,1), &
               njmin(1,1), &
               njmax(1,1), &
               nkmin_xi(1,1), &
               nkmin_eta(1,1),stat=ier)
      if (ier /= 0) stop 'Error in allocate 14'
    endif

    ! initializes boundary indices only during first pass, we need then the stored index values
    ! for creating mass matrices for Stacey in second pass
    nimin(:,:) = 0; nimax(:,:) = 0
    njmin(:,:) = 0; njmax(:,:) = 0
    nkmin_xi(:,:) = 0; nkmin_eta(:,:) = 0
  endif


  ! MPI cut-planes parameters along xi and along eta
  allocate(iMPIcut_xi(2,nspec), &
           iMPIcut_eta(2,nspec),stat=ier)
  if (ier /= 0) stop 'Error in allocate 15'

  iMPIcut_xi(:,:) = .false.; iMPIcut_eta(:,:) = .false.

  ! MPI buffer indices
  !
  ! define maximum size for message buffers
  ! use number of elements found in the mantle since it is the largest region
  NGLOB2DMAX_XY = max(NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE))
  ! 1-D buffers
  NGLOB1D_RADIAL_MAX = maxval(NGLOB1D_RADIAL_CORNER(iregion_code,:))

  if (ipass == 1) then
    allocate(iboolleft_xi(NGLOB2DMAX_XMIN_XMAX(iregion_code)), &
             iboolright_xi(NGLOB2DMAX_XMIN_XMAX(iregion_code)), &
             iboolleft_eta(NGLOB2DMAX_YMIN_YMAX(iregion_code)), &
             iboolright_eta(NGLOB2DMAX_YMIN_YMAX(iregion_code)), &
             stat=ier)
    if (ier /= 0) stop 'Error in allocate 15b'

    iboolleft_xi(:) = 0; iboolright_xi(:) = 0
    iboolleft_eta(:) = 0; iboolright_eta(:) = 0

    allocate(ibool1D_leftxi_lefteta(NGLOB1D_RADIAL_MAX), &
             ibool1D_rightxi_lefteta(NGLOB1D_RADIAL_MAX), &
             ibool1D_leftxi_righteta(NGLOB1D_RADIAL_MAX), &
             ibool1D_rightxi_righteta(NGLOB1D_RADIAL_MAX), &
             stat=ier)
    if (ier /= 0) stop 'Error in allocate 15c'

    ibool1D_leftxi_lefteta(:) = 0; ibool1D_rightxi_lefteta(:) = 0
    ibool1D_leftxi_righteta(:) = 0; ibool1D_rightxi_righteta(:) = 0

    allocate(xyz1D_leftxi_lefteta(NGLOB1D_RADIAL_MAX,NDIM), &
             xyz1D_rightxi_lefteta(NGLOB1D_RADIAL_MAX,NDIM), &
             xyz1D_leftxi_righteta(NGLOB1D_RADIAL_MAX,NDIM), &
             xyz1D_rightxi_righteta(NGLOB1D_RADIAL_MAX,NDIM), &
             stat=ier)
    if (ier /= 0) stop 'Error in allocate 15c'

    xyz1D_leftxi_lefteta(:,:) = 0.0; xyz1D_rightxi_lefteta(:,:) = 0.0
    xyz1D_leftxi_righteta(:,:) = 0.0; xyz1D_rightxi_righteta(:,:) = 0.0

    allocate(iboolcorner(NGLOB1D_RADIAL(iregion_code),NUMCORNERS_SHARED), &
             iboolfaces(NGLOB2DMAX_XY,NUMFACES_SHARED), &
             stat=ier)
    if (ier /= 0) stop 'Error in allocate 15b'

    iboolcorner(:,:) = 0; iboolfaces(:,:) = 0
  endif

  ! store and save the final arrays only in the second pass
  ! therefore in the first pass some arrays can be allocated with a dummy size
  if (ipass == 1) then
    nspec_actually = 1
  else
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
  if (ier /= 0) stop 'Error in allocate 16'

  xixstore(:,:,:,:) = 0.0; xiystore(:,:,:,:) = 0.0; xizstore(:,:,:,:) = 0.0
  etaxstore(:,:,:,:) = 0.0; etaystore(:,:,:,:) = 0.0; etazstore(:,:,:,:) = 0.0
  gammaxstore(:,:,:,:) = 0.0; gammaystore(:,:,:,:) = 0.0; gammazstore(:,:,:,:) = 0.0

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
  if (ier /= 0) stop 'Error in allocate 17'

  ibelm_moho_top(:) = 0; ibelm_moho_bot(:) = 0
  ibelm_400_top(:) = 0; ibelm_400_bot(:) = 0
  ibelm_670_top(:) = 0; ibelm_670_bot(:) = 0
  normal_moho(:,:,:,:) = 0.0; normal_400(:,:,:,:) = 0.0; normal_670(:,:,:,:) = 0.0
  jacobian2D_moho(:,:,:) = 0.0; jacobian2D_400(:,:,:) = 0.0; jacobian2D_670(:,:,:) = 0.0

  end subroutine crm_allocate_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine crm_setup_layers(iregion_code,ipass, &
                              NEX_PER_PROC_ETA)

  use meshfem3D_par, only: &
    nspec, &
    ibool,idoubling,is_on_a_slice_edge, &
    xstore,ystore,zstore, &
    NGLLX,NGLLY,NGLLZ, &
    IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE, &
    R670,RMOHO,R400,RMIDDLE_CRUST,MAX_NUMBER_OF_MESH_LAYERS, &
    ner,r_top,r_bottom

  use meshfem3D_models_par, only: &
    CASE_3D,SUPPRESS_CRUSTAL_MESH,ONE_CRUST,REGIONAL_MOHO_MESH

  use regions_mesh_par

  use regions_mesh_par2

  implicit none

  integer,intent(in) :: iregion_code
  integer,intent(in) :: ipass
  integer :: NEX_PER_PROC_ETA

  ! local parameters
  integer :: cpt
  integer :: i,ier

  ! initializes element layers
  call initialize_layers(ipass,xigll,yigll,zigll,wxgll,wygll,wzgll, &
                         shape3D,dershape3D,shape2D_x,shape2D_y,shape2D_bottom,shape2D_top, &
                         dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                         iaddx,iaddy,iaddz, &
                         nspec,xstore,ystore,zstore,ibool,idoubling, &
                         iboun,iMPIcut_xi,iMPIcut_eta,ispec2D_moho_top,ispec2D_moho_bot, &
                         ispec2D_400_top,ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot, &
                         NEX_PER_PROC_ETA,nex_eta_moho,RMOHO,R400,R670,r_moho,r_400,r_670, &
                         ONE_CRUST,NUMBER_OF_MESH_LAYERS,layer_shift, &
                         iregion_code,ifirst_region,ilast_region, &
                         first_layer_aniso,last_layer_aniso,is_on_a_slice_edge)

  ! to consider anisotropic elements first and to build the mesh from the bottom to the top of the region
  allocate (perm_layer(ifirst_region:ilast_region),stat=ier)
  if (ier /= 0) stop 'Error in allocate 18'
  perm_layer = (/ (i, i=ilast_region,ifirst_region,-1) /)

  if (iregion_code == IREGION_CRUST_MANTLE) then
    cpt = 3
    perm_layer(1) = first_layer_aniso
    perm_layer(2) = last_layer_aniso
    do i = ilast_region,ifirst_region,-1
      if (i /= first_layer_aniso .and. i /= last_layer_aniso) then
        perm_layer(cpt) = i
        cpt = cpt + 1
      endif
    enddo
  endif

  ! crustal layer stretching: element layer's top and bottom radii will get stretched when in crust
  ! (number of element layers in crust can vary for different resolutions and 1chunk simulations)
  allocate(stretch_tab(2,ner(1)),stat=ier)
  if (ier /= 0) stop 'Error in allocate 19'
  if (CASE_3D .and. iregion_code == IREGION_CRUST_MANTLE .and. .not. SUPPRESS_CRUSTAL_MESH) then
    ! stretching function determines top and bottom of each element layer in the
    ! crust region (between r_top(1) and r_bottom(1)), where ner(1) is the
    ! number of element layers in this crust region

    ! differentiate between regional meshes or global meshes
    if (REGIONAL_MOHO_MESH) then
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

  subroutine crm_setup_indexing(npointot)

! creates global indexing array ibool

  use constants, only: NGLLX,NGLLY,NGLLZ,ZERO

  use meshfem3d_par, only: &
    nspec,nglob, &
    ibool,xstore,ystore,zstore, &
    myrank

  use regions_mesh_par2

  implicit none

  ! number of total points
  integer,intent(in) :: npointot

  ! local parameters
  ! variables for creating array ibool
  double precision, dimension(:), allocatable :: xp,yp,zp
  integer, dimension(:), allocatable :: locval
  logical, dimension(:), allocatable :: ifseg

  integer :: nglob_new
  integer :: ieoff,ilocnum,ier
  integer :: i,j,k,ispec
  character(len=MAX_STRING_LEN) :: errmsg

  ! allocate memory for arrays
  allocate(locval(npointot), &
          ifseg(npointot), &
          xp(npointot), &
          yp(npointot), &
          zp(npointot),stat=ier)
  if (ier /= 0) stop 'Error in allocate 20'

  locval(:) = 0
  ifseg(:) = .false.
  xp(:) = ZERO
  yp(:) = ZERO
  zp(:) = ZERO

  ! we need to create a copy of the x, y and z arrays because sorting in get_global will swap
  ! these arrays and therefore destroy them
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ispec,ieoff,ilocnum,i,j,k)
!$OMP DO
  do ispec = 1,nspec
    ieoff = NGLLX * NGLLY * NGLLZ * (ispec-1)
    ilocnum = 0
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ! increases point counter
          !ilocnum = ilocnum + 1
          ! without dependency
          ilocnum = i + ((j-1) + (k-1)*NGLLY) * NGLLX
          ! fills 1D arrays
          xp(ilocnum+ieoff) = xstore(i,j,k,ispec)
          yp(ilocnum+ieoff) = ystore(i,j,k,ispec)
          zp(ilocnum+ieoff) = zstore(i,j,k,ispec)
        enddo
      enddo
    enddo
  enddo
!$OMP ENDDO
!$OMP END PARALLEL

  call get_global(npointot,xp,yp,zp,ibool,locval,ifseg,nglob_new)

  deallocate(xp,yp,zp)
  deallocate(locval,ifseg)

  ! check that number of points found equals theoretical value
  if (nglob_new /= nglob) then
    write(errmsg,*) 'incorrect total number of points found: myrank,nglob_new,nglob = ',myrank,nglob_new,nglob
    call exit_MPI(myrank,errmsg)
  endif
  if (minval(ibool) /= 1 .or. maxval(ibool) /= nglob) &
    call exit_MPI(myrank,'incorrect global numbering')

  ! creates a new indirect addressing to reduce cache misses in memory access in the solver
  ! this is *critical* to improve performance in the solver
  call get_global_indirect_addressing(nspec,nglob,ibool)

  ! checks again
  if (minval(ibool) /= 1 .or. maxval(ibool) /= nglob) &
    call exit_MPI(myrank,'incorrect global numbering after sorting')

  end subroutine crm_setup_indexing

!
!-------------------------------------------------------------------------------------------------
!

  subroutine crm_setup_mpi_buffers(npointot,iregion_code)

! sets up MPI cutplane arrays

  use meshfem3d_par, only: &
    nspec,ibool,idoubling, &
    xstore,ystore,zstore, &
    NGLLX,NGLLY,NGLLZ, &
    NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER, &
    NSPEC2D_XI_FACE,NSPEC2D_ETA_FACE, &
    NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX

  use regions_mesh_par2

  use MPI_interfaces_par

  implicit none

  ! number of total points
  integer,intent(in) :: npointot

  integer,intent(in) :: iregion_code

  ! local parameters
  logical, dimension(:), allocatable :: mask_ibool
  integer :: ier

  ! arrays mask_ibool(npointot) used to save memory
  ! allocate memory for arrays
  allocate(mask_ibool(npointot), &
          stat=ier)
  if (ier /= 0) stop 'Error in allocate 20b'

  ! initializes
  npoin2D_xi_all(:) = 0
  npoin2D_eta_all(:) = 0
  iboolleft_xi(:) = 0
  iboolleft_eta(:) = 0
  iboolright_xi(:) = 0
  iboolright_eta(:) = 0

  ! gets MPI buffer indices
  call get_MPI_cutplanes_xi(prname,nspec,iMPIcut_xi,ibool, &
                  xstore,ystore,zstore,mask_ibool,npointot, &
                  NSPEC2D_ETA_FACE,iregion_code,npoin2D_xi, &
                  iboolleft_xi,iboolright_xi, &
                  npoin2D_xi_all,NGLOB2DMAX_XMIN_XMAX(iregion_code))

  call get_MPI_cutplanes_eta(prname,nspec,iMPIcut_eta,ibool, &
                  xstore,ystore,zstore,mask_ibool,npointot, &
                  NSPEC2D_XI_FACE,iregion_code,npoin2D_eta, &
                  iboolleft_eta,iboolright_eta, &
                  npoin2D_eta_all,NGLOB2DMAX_YMIN_YMAX(iregion_code))

  call get_MPI_1D_buffers(prname,nspec,iMPIcut_xi,iMPIcut_eta, &
                  ibool,idoubling, &
                  xstore,ystore,zstore,mask_ibool,npointot, &
                  NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER,iregion_code, &
                  ibool1D_leftxi_lefteta,ibool1D_rightxi_lefteta, &
                  ibool1D_leftxi_righteta,ibool1D_rightxi_righteta, &
                  xyz1D_leftxi_lefteta,xyz1D_rightxi_lefteta, &
                  xyz1D_leftxi_righteta,xyz1D_rightxi_righteta, &
                  NGLOB1D_RADIAL_MAX)

  deallocate(mask_ibool)

  end subroutine crm_setup_mpi_buffers


!
!-------------------------------------------------------------------------------
!


  subroutine crm_free_MPI_arrays(iregion_code)

  use MPI_interfaces_par

  use MPI_crust_mantle_par
  use MPI_outer_core_par
  use MPI_inner_core_par

  implicit none

  integer,intent(in):: iregion_code

  ! free memory
  select case (iregion_code)
  case (IREGION_CRUST_MANTLE)
    ! crust mantle
    deallocate(phase_ispec_inner_crust_mantle)
    deallocate(num_elem_colors_crust_mantle)
  case (IREGION_OUTER_CORE)
    ! outer core
    deallocate(phase_ispec_inner_outer_core)
    deallocate(num_elem_colors_outer_core)
  case (IREGION_INNER_CORE)
    ! inner core
    deallocate(phase_ispec_inner_inner_core)
    deallocate(num_elem_colors_inner_core)
  end select

  end subroutine crm_free_MPI_arrays


!
!-------------------------------------------------------------------------------------------------
!

  subroutine crm_fill_global_meshes()

! fills in global mesh arrays
! (defines on global points xstore_glob(iglob) rather than locally as in xstore(i,j,k,ispec))

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  use meshfem3D_par, only: &
    nspec,nglob, &
    myrank,ibool,xstore,ystore,zstore, &
    xstore_glob,ystore_glob,zstore_glob

  implicit none

  ! local parameters
  integer :: i,j,k,ispec,iglob,ier

  ! safety check
  if (allocated(xstore_glob)) stop 'Error global mesh already allocated'
  if (nglob < 1) stop 'Invalid number of global mesh points'
  if (nspec < 1) stop 'Invalid number of elements'

  ! allocates temporary global mesh
  allocate(xstore_glob(nglob),ystore_glob(nglob),zstore_glob(nglob), &
           stat=ier)
  if (ier /= 0) call exit_mpi(myrank,'Error allocating temporary global mesh arrays')

  ! fill CUSTOM_REAL arrays
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! distinguish between single and double precision for reals
          xstore_glob(iglob) = real(xstore(i,j,k,ispec), kind=CUSTOM_REAL)
          ystore_glob(iglob) = real(ystore(i,j,k,ispec), kind=CUSTOM_REAL)
          zstore_glob(iglob) = real(zstore(i,j,k,ispec), kind=CUSTOM_REAL)
        enddo
      enddo
    enddo
  enddo

  ! synchronize processes
  call synchronize_all()

  end subroutine crm_fill_global_meshes
