!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            August 2013
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

!-------------------------------------------------------------------------------
!> \file get_absorb_adios.f90
!! \brief Function to write stacey boundary condition to disk with ADIOS.
!! \author MPBL
!-------------------------------------------------------------------------------

!===============================================================================
!> \brief Main routine to save the arrays from the mesher to the solver with the
!!        help of ADIOS
!! \param myrank The MPI rank of the current process
!! \param nspec Number of GLL points per element
!! \param nglob Number of mesh points
!! \param idoubling Array of information on every mesh point
!! \param ibool Array of information on every mesh point
!! \param iregion_code The region the absorbing conditon is written for. Check
!!                     constant.h files to see what these regions are.
!! \param xstore Array with the x coordinates of the mesh points
!! \param ystore Array with the y coordinates of the mesh points
!! \param zstore Array with the z coordinates of the mesh points
!! \param NSPEC2DMAX_XMIN_XMAX Integer to compute the size of the arrays
!!                             in argument
!! \param NSPEC2DMAX_YMIN_YMAX Integer to compute the size of the arrays
!!                             in argument
!! \param NSPEC2D_TOP Integer to compute the size of the arrays
!!                    in argument
!! \param NSPEC2D_BOTTOM Integer to compute the size of the arrays
!!                       in argument
subroutine save_arrays_solver_adios(myrank,nspec,nglob,idoubling,ibool, &
                    iregion_code,xstore,ystore,zstore, &
                    NSPEC2DMAX_XMIN_XMAX, NSPEC2DMAX_YMIN_YMAX, &
                    NSPEC2D_TOP,NSPEC2D_BOTTOM)

  use mpi
  use adios_write_mod

  use constants

  use meshfem3D_models_par,only: &
    OCEANS,TRANSVERSE_ISOTROPY,HETEROGEN_3D_MANTLE,ANISOTROPIC_3D_MANTLE, &
    ANISOTROPIC_INNER_CORE,ATTENUATION

  use meshfem3D_par,only: &
    NCHUNKS,ABSORBING_CONDITIONS,SAVE_MESH_FILES, LOCAL_PATH, &
    ADIOS_FOR_SOLVER_MESHFILES

  use create_regions_mesh_par2,only: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
    gammaxstore,gammaystore,gammazstore, &
    rhostore,dvpstore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore,&
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
    rmassx,rmassy,rmassz,rmass_ocean_load, &
    ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
    normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top, &
    jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax, &
    jacobian2D_bottom,jacobian2D_top, &
    rho_vp,rho_vs, &
    nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
    ispec_is_tiso,tau_s,T_c_source,tau_e_store,Qmu_store, &
    prname, nspec_actually, nspec_ani, nspec_stacey, nglob_xy, nglob_oceans

  implicit none

  integer :: myrank
  integer :: nspec,nglob

  ! doubling mesh flag
  integer, dimension(nspec) :: idoubling
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  integer :: iregion_code

  ! arrays with the mesh in double precision
  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  ! boundary parameters locator
  integer :: NSPEC2D_TOP,NSPEC2D_BOTTOM, &
      NSPEC2DMAX_XMIN_XMAX, NSPEC2DMAX_YMIN_YMAX

  ! local parameters
  integer :: i,j,k,ispec,iglob,ier
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: tmp_array_x, &
      tmp_array_y, tmp_array_z

  ! local parameters
  character(len=150) :: reg_name, outputname, group_name
  integer :: ierr, sizeprocs, comm, local_dim
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, adios_handle, varid
  integer(kind=8)         :: adios_groupsize, adios_totalsize

  ! create a prefix for the file name such as LOCAL_PATH/regX_
  call create_name_database_adios(reg_name,iregion_code,LOCAL_PATH)

  !---------------------------------------------------------
  !--- Solver data arrays ----------------------------------
  !---------------------------------------------------------

  ! create the name for the database of the current slide and region
  outputname = trim(reg_name) // "solver_data.bp"

  ! save arrays for the solver to run.
  write(group_name,"('SPECFEM3D_GLOBE_ARRAYS_SOLVER_reg',i1)") iregion_code
  call world_size(sizeprocs) ! TODO keep it in parameters
  ! Alias COMM_WORLD to use ADIOS
  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
  ! set the adios group size to 0 before incremented by calls to
  ! helpers functions.
  group_size_inc = 0
  call adios_declare_group(adios_group, group_name, &
      "", 0, adios_err)
  ! We set the transport method to 'MPI'. This seems to be the correct choice
  ! for now. We might want to move this to the constant.h file later on.
  call adios_select_method(adios_group, "MPI", "", "", adios_err)

  !--- Define ADIOS variables -----------------------------
  ! save nspec and nglob, to be used in combine_paraview_data
  call define_adios_integer_scalar (adios_group, "nspec", "", &
      group_size_inc)
  call define_adios_integer_scalar (adios_group, "nglob", "", &
      group_size_inc)

  local_dim = nglob
  call define_adios_global_real_1d_array(adios_group, "xstore", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "ystore", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "zstore", &
      local_dim, group_size_inc)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec
  call define_adios_global_real_1d_array(adios_group, "rhostore", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "kappavstore", &
      local_dim, group_size_inc)
  call define_adios_global_integer_1d_array(adios_group, "ibool", &
      local_dim, group_size_inc)
  if(iregion_code /= IREGION_OUTER_CORE) then
    if(.not. (ANISOTROPIC_3D_MANTLE .and. &
        iregion_code == IREGION_CRUST_MANTLE)) then
      call define_adios_global_real_1d_array(adios_group, "muvstore", &
          local_dim, group_size_inc)
    endif
    if(TRANSVERSE_ISOTROPY) then
      if(iregion_code == IREGION_CRUST_MANTLE .and. &
          .not. ANISOTROPIC_3D_MANTLE) then
        call define_adios_global_real_1d_array(adios_group, "kappahstore", &
            local_dim, group_size_inc)
        call define_adios_global_real_1d_array(adios_group, "muhstore", &
            local_dim, group_size_inc)
        call define_adios_global_real_1d_array(adios_group, "eta_anisostore", &
            local_dim, group_size_inc)
      endif
    endif
  endif

  local_dim = nspec
  call define_adios_global_integer_1d_array(adios_group, "idoubling", &
      local_dim, group_size_inc)
  call define_adios_global_integer_1d_array(adios_group, "ispec_is_tiso", &
      local_dim, group_size_inc)
  local_dim = NGLLX * NGLLY * NGLLZ * nspec_actually
  call define_adios_global_real_1d_array(adios_group, "xixstore", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "xiystore", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "xizstore", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "etaxstore", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "etaystore", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "etazstore", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "gammaxstore", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "gammaystore", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "gammazstore", &
      local_dim, group_size_inc)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_ani
  if(iregion_code /= IREGION_OUTER_CORE) then
    !   save anisotropy in the inner core only
    if(ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) then
      call define_adios_global_real_1d_array(adios_group, "c11store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c33store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c12store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c13store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c44store", &
          local_dim, group_size_inc)
    endif
    if(ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
      call define_adios_global_real_1d_array(adios_group, "c11store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c12store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c13store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c14store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c15store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c16store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c22store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c23store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c24store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c25store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c26store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c33store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c34store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c35store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c36store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c44store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c45store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c46store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c55store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c56store", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "c66store", &
          local_dim, group_size_inc)
    endif
  endif

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_stacey
  if(ABSORBING_CONDITIONS) then
    if(iregion_code == IREGION_CRUST_MANTLE) then
      call define_adios_global_real_1d_array(adios_group, "rho_vp", &
          local_dim, group_size_inc)
      call define_adios_global_real_1d_array(adios_group, "rho_vs", &
          local_dim, group_size_inc)
    else if(iregion_code == IREGION_OUTER_CORE) then
      call define_adios_global_real_1d_array(adios_group, "rho_vp", &
          local_dim, group_size_inc)
    endif
  endif

  local_dim = nglob_xy
  if(NCHUNKS /= 6 .and. ABSORBING_CONDITIONS .and. &
      iregion_code == IREGION_CRUST_MANTLE) then
    call define_adios_global_real_1d_array(adios_group, "rmassx", &
        local_dim, group_size_inc)
    call define_adios_global_real_1d_array(adios_group, "rmassy", &
        local_dim, group_size_inc)
  endif
  local_dim = nglob
  call define_adios_global_real_1d_array(adios_group, "rmassz", &
      local_dim, group_size_inc)

  local_dim = nglob_oceans
  if(OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) then
    call define_adios_global_real_1d_array(adios_group, "rmass_ocean_load", &
        local_dim, group_size_inc)
  endif

  !--- Open an ADIOS handler to the restart file. ---------
  call adios_open (adios_handle, group_name, &
      outputname, "w", comm, adios_err);
  call adios_group_size (adios_handle, group_size_inc, &
                         adios_totalsize, adios_err)

  ! mesh topology

  ! mesh arrays used in the solver to locate source and receivers
  ! and for anisotropy and gravity, save in single precision
  ! use tmp_array for temporary storage to perform conversion
  allocate(tmp_array_x(nglob),stat=ier)
  if( ier /=0 ) call exit_MPI(myrank,&
      'error allocating temporary array for mesh topology')
  allocate(tmp_array_y(nglob),stat=ier)
  if( ier /=0 ) call exit_MPI(myrank,&
      'error allocating temporary array for mesh topology')
  allocate(tmp_array_z(nglob),stat=ier)
  if( ier /=0 ) call exit_MPI(myrank,&
      'error allocating temporary array for mesh topology')

  !--- x coordinate
  tmp_array_x(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            tmp_array_x(iglob) = sngl(xstore(i,j,k,ispec))
          else
            tmp_array_x(iglob) = xstore(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo
  !--- y coordinate
  tmp_array_y(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            tmp_array_y(iglob) = sngl(ystore(i,j,k,ispec))
          else
            tmp_array_y(iglob) = ystore(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo
  !--- z coordinate
  tmp_array_z(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            tmp_array_z(iglob) = sngl(zstore(i,j,k,ispec))
          else
            tmp_array_z(iglob) = zstore(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo

  !--- Schedule writes for the previously defined ADIOS variables
  ! save nspec and nglob, to be used in combine_paraview_data
  call adios_write(adios_handle, "nspec", nspec, adios_err)
  call adios_write(adios_handle, "nglob", nglob, adios_err)

  local_dim = nglob
  call adios_set_path (adios_handle, "xstore", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", tmp_array_x, adios_err)

  call adios_set_path (adios_handle, "ystore", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", tmp_array_y, adios_err)

  call adios_set_path (adios_handle, "zstore", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", tmp_array_z, adios_err)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec
  call adios_set_path (adios_handle, "rhostore", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", rhostore, adios_err)

  call adios_set_path (adios_handle, "kappavstore", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", kappavstore, adios_err)

  call adios_set_path (adios_handle, "ibool", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", ibool, adios_err)

  if(iregion_code /= IREGION_OUTER_CORE) then
    if(.not. (ANISOTROPIC_3D_MANTLE .and. &
        iregion_code == IREGION_CRUST_MANTLE)) then
      call adios_set_path (adios_handle, "muvstore", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", muvstore, adios_err)
    endif
    if(TRANSVERSE_ISOTROPY) then
      if(iregion_code == IREGION_CRUST_MANTLE .and. &
          .not. ANISOTROPIC_3D_MANTLE) then
        call adios_set_path (adios_handle, "kappahstore", adios_err)
        call write_1D_global_array_adios_dims(adios_handle, myrank, &
            local_dim, sizeprocs)
        call adios_write(adios_handle, "array", kappahstore, adios_err)

        call adios_set_path (adios_handle, "muhstore", adios_err)
        call write_1D_global_array_adios_dims(adios_handle, myrank, &
            local_dim, sizeprocs)
        call adios_write(adios_handle, "array", muhstore, adios_err)

        call adios_set_path (adios_handle, "eta_anisostore", adios_err)
        call write_1D_global_array_adios_dims(adios_handle, myrank, &
            local_dim, sizeprocs)
        call adios_write(adios_handle, "array", eta_anisostore, adios_err)
      endif
    endif
  endif

  local_dim = nspec
  call adios_set_path (adios_handle, "idoubling", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", idoubling, adios_err)

  call adios_set_path (adios_handle, "ispec_is_tiso", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", ispec_is_tiso, adios_err)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_actually
  call adios_set_path (adios_handle, "xixstore", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", xixstore, adios_err)

  call adios_set_path (adios_handle, "xiystore", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", xiystore, adios_err)

  call adios_set_path (adios_handle, "xizstore", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", xizstore, adios_err)

  call adios_set_path (adios_handle, "etaxstore", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", etaxstore, adios_err)

  call adios_set_path (adios_handle, "etaystore", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", etaystore, adios_err)

  call adios_set_path (adios_handle, "etazstore", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", etazstore, adios_err)

  call adios_set_path (adios_handle, "gammaxstore", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", gammaxstore, adios_err)

  call adios_set_path (adios_handle, "gammaystore", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", gammaystore, adios_err)

  call adios_set_path (adios_handle, "gammazstore", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", gammazstore, adios_err)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_ani
  if(iregion_code /= IREGION_OUTER_CORE) then
    !   save anisotropy in the inner core only
    if(ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) then
      call adios_set_path (adios_handle, "c11store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c11store, adios_err)

      call adios_set_path (adios_handle, "c33store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c33store, adios_err)

      call adios_set_path (adios_handle, "c12store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c12store, adios_err)

      call adios_set_path (adios_handle, "c13store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c13store, adios_err)

      call adios_set_path (adios_handle, "c44store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c44store, adios_err)
    endif
    if(ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
      call adios_set_path (adios_handle, "c11store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c11store, adios_err)

      call adios_set_path (adios_handle, "c12store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c12store, adios_err)

      call adios_set_path (adios_handle, "c13store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c13store, adios_err)

      call adios_set_path (adios_handle, "c14store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c14store, adios_err)

      call adios_set_path (adios_handle, "c15store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c15store, adios_err)

      call adios_set_path (adios_handle, "c16store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c16store, adios_err)

      call adios_set_path (adios_handle, "c22store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c22store, adios_err)

      call adios_set_path (adios_handle, "c23store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c23store, adios_err)

      call adios_set_path (adios_handle, "c24store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c24store, adios_err)

      call adios_set_path (adios_handle, "c25store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c25store, adios_err)

      call adios_set_path (adios_handle, "c26store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c26store, adios_err)

      call adios_set_path (adios_handle, "c33store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c33store, adios_err)

      call adios_set_path (adios_handle, "c34store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c34store, adios_err)

      call adios_set_path (adios_handle, "c35store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c35store, adios_err)

      call adios_set_path (adios_handle, "c36store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c36store, adios_err)

      call adios_set_path (adios_handle, "c44store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c44store, adios_err)

      call adios_set_path (adios_handle, "c45store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c45store, adios_err)

      call adios_set_path (adios_handle, "c46store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c46store, adios_err)

      call adios_set_path (adios_handle, "c55store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c55store, adios_err)

      call adios_set_path (adios_handle, "c56store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c56store, adios_err)

      call adios_set_path (adios_handle, "c66store", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", c66store, adios_err)
    endif
  endif

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_stacey
  if(ABSORBING_CONDITIONS) then
    if(iregion_code == IREGION_CRUST_MANTLE) then
      call adios_set_path (adios_handle, "rho_vp", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", rho_vp, adios_err)

      call adios_set_path (adios_handle, "rho_vs", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", rho_vs, adios_err)

    else if(iregion_code == IREGION_OUTER_CORE) then
      call adios_set_path (adios_handle, "rho_vp", adios_err)
      call write_1D_global_array_adios_dims(adios_handle, myrank, &
          local_dim, sizeprocs)
      call adios_write(adios_handle, "array", rho_vp, adios_err)
    endif
  endif

  local_dim = nglob_xy
  if(NCHUNKS /= 6 .and. ABSORBING_CONDITIONS .and. &
      iregion_code == IREGION_CRUST_MANTLE) then
    call adios_set_path (adios_handle, "rmassx", adios_err)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        local_dim, sizeprocs)
    call adios_write(adios_handle, "array", rmassx, adios_err)

    call adios_set_path (adios_handle, "rmassy", adios_err)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        local_dim, sizeprocs)
    call adios_write(adios_handle, "array", rmassy, adios_err)
  endif

  local_dim = nglob
  call adios_set_path (adios_handle, "rmassz", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", rmassz, adios_err)

  local_dim = nglob_oceans
  if(OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) then
    call adios_set_path (adios_handle, "rmass_ocean_load", adios_err)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        local_dim, sizeprocs)
    call adios_write(adios_handle, "array", rmass_ocean_load, adios_err)
    if(minval(rmass_ocean_load) <= 0._CUSTOM_REAL) &
        call exit_MPI(myrank,'negative mass matrix term for the oceans')
  endif

  !--- Reset the path to zero and perform the actual write to disk
  call adios_set_path (adios_handle, "", adios_err)
  call adios_close(adios_handle, adios_err)

  ! Clean the temporary arrays containing the node information
  deallocate(tmp_array_x)
  deallocate(tmp_array_y)
  deallocate(tmp_array_z)

  !---------------------------------------------------------
  !--- Boundary arrays -------------------------------------
  !---------------------------------------------------------

  ! Postpend the actual file name.
  outputname = trim(reg_name) // "boundary.bp"

  ! save boundary arrays in ADIOS files
  write(group_name,"('SPECFEM3D_GLOBE_BOUNDARY_reg',i1)") iregion_code
  ! set the adios group size to 0 before incremented by calls to
  ! helpers functions.
  group_size_inc = 0
  call adios_declare_group(adios_group, group_name, &
      "", 0, adios_err)
  call adios_select_method(adios_group, "MPI", "", "", adios_err)

  !--- Define ADIOS variables -----------------------------
  call define_adios_integer_scalar (adios_group, "nspec2D_xmin", "", &
      group_size_inc)
  call define_adios_integer_scalar (adios_group, "nspec2D_xmax", "", &
      group_size_inc)
  call define_adios_integer_scalar (adios_group, "nspec2D_ymin", "", &
      group_size_inc)
  call define_adios_integer_scalar (adios_group, "nspec2D_ymax", "", &
      group_size_inc)
  call define_adios_integer_scalar (adios_group, "NSPEC2D_BOTTOM", "", &
      group_size_inc)
  call define_adios_integer_scalar (adios_group, "NSPEC2D_TOP", "", &
      group_size_inc)

  !local_dim = NSPEC2DMAX_XMIN_YMAX
  local_dim = size (ibelm_xmin)
  call define_adios_global_integer_1d_array(adios_group, "ibelm_xmin", &
      local_dim, group_size_inc)
  call define_adios_global_integer_1d_array(adios_group, "ibelm_xmax", &
      local_dim, group_size_inc)

  !local_dim = NSPEC2DMAX_YMIN_YMAX
  local_dim = size (ibelm_ymin)
  call define_adios_global_integer_1d_array(adios_group, "ibelm_ymin", &
      local_dim, group_size_inc)
  call define_adios_global_integer_1d_array(adios_group, "ibelm_ymax", &
      local_dim, group_size_inc)

  !local_dim = NSPEC2D_BOTTOM
  local_dim = size (ibelm_bottom)
  call define_adios_global_integer_1d_array(adios_group, "ibelm_bottom", &
      local_dim, group_size_inc)

  !local_dim = NSPEC2D_TOP
  local_dim = size (ibelm_top)
  call define_adios_global_integer_1d_array(adios_group, "ibelm_top", &
      local_dim, group_size_inc)

  !local_dim = NDIM*NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX
  local_dim = size (normal_xmin)
  call define_adios_global_real_1d_array(adios_group, "normal_xmin", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "normal_xmax", &
      local_dim, group_size_inc)

  !local_dim = NDIM*NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX
  local_dim = size (normal_ymin)
  call define_adios_global_real_1d_array(adios_group, "normal_ymin", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "normal_ymax", &
      local_dim, group_size_inc)

  !local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_BOTTOM
  local_dim = size (normal_bottom)
  call define_adios_global_real_1d_array(adios_group, "normal_bottom", &
      local_dim, group_size_inc)

  !local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_TOP
  local_dim = size (normal_top)
  call define_adios_global_real_1d_array(adios_group, "normal_top", &
      local_dim, group_size_inc)

  !local_dim = NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX
  local_dim = size (jacobian2D_xmin)
  call define_adios_global_real_1d_array(adios_group, "jacobian2D_xmin", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "jacobian2D_xmax", &
      local_dim, group_size_inc)
  !local_dim = NDIM*NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX
  local_dim = size (jacobian2D_ymin)
  call define_adios_global_real_1d_array(adios_group, "jacobian2D_ymin", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "jacobian2D_ymax", &
      local_dim, group_size_inc)
  !local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_BOTTOM
  local_dim = size (jacobian2D_bottom)
  call define_adios_global_real_1d_array(adios_group, "jacobian2D_bottom", &
      local_dim, group_size_inc)
  !local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_TOP
  local_dim = size (jacobian2D_top)
  call define_adios_global_real_1d_array(adios_group, "jacobian2D_top", &
      local_dim, group_size_inc)

  !--- Open an ADIOS handler to the restart file. ---------
  call adios_open (adios_handle, group_name, &
      outputname, "w", comm, adios_err);
  call adios_group_size (adios_handle, group_size_inc, &
                         adios_totalsize, adios_err)

  !--- Schedule writes for the previously defined ADIOS variables
  call adios_write(adios_handle, "nspec2D_xmin", nspec2D_xmin, adios_err)
  call adios_write(adios_handle, "nspec2D_xmax", nspec2D_xmax, adios_err)
  call adios_write(adios_handle, "nspec2D_ymin", nspec2D_ymin, adios_err)
  call adios_write(adios_handle, "nspec2D_ymax", nspec2D_ymax, adios_err)
  call adios_write(adios_handle, "NSPEC2D_BOTTOM", NSPEC2D_BOTTOM, adios_err)
  call adios_write(adios_handle, "NSPEC2D_TOP", NSPEC2D_TOP, adios_err)

  !local_dim = NSPEC2DMAX_XMIN_XMAX
  local_dim = size (ibelm_xmin)
  call adios_set_path (adios_handle, "ibelm_xmin", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", ibelm_xmin, adios_err)
  call adios_set_path (adios_handle, "ibelm_xmax", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", ibelm_xmax, adios_err)

  !local_dim = NSPEC2DMAX_YMIN_YMAX
  local_dim = size (ibelm_ymin)
  call adios_set_path (adios_handle, "ibelm_ymin", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", ibelm_ymin, adios_err)
  call adios_set_path (adios_handle, "ibelm_ymax", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", ibelm_ymax, adios_err)

  !local_dim = NSPEC2D_BOTTOM
  local_dim = size (ibelm_bottom)
  call adios_set_path (adios_handle, "ibelm_bottom", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", ibelm_bottom, adios_err)

  !local_dim = NSPEC2D_TOP
  local_dim = size (ibelm_top)
  call adios_set_path (adios_handle, "ibelm_top", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", ibelm_top, adios_err)

  !local_dim = NDIM*NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX
  local_dim = size (normal_xmin)
  call adios_set_path (adios_handle, "normal_xmin", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", normal_xmin, adios_err)
  call adios_set_path (adios_handle, "normal_xmax", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", normal_xmax, adios_err)

  !local_dim = NDIM*NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX
  local_dim = size (normal_ymin)
  call adios_set_path (adios_handle, "normal_ymin", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", normal_ymin, adios_err)
  call adios_set_path (adios_handle, "normal_ymax", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", normal_ymax, adios_err)

  !local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_BOTTOM
  local_dim = size (normal_bottom)
  call adios_set_path (adios_handle, "normal_bottom", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", normal_bottom, adios_err)

  !local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_TOP
  local_dim = size (normal_top)
  call adios_set_path (adios_handle, "normal_top", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", normal_top, adios_err)

  !local_dim = NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX
  local_dim = size (jacobian2D_xmin)
  call adios_set_path (adios_handle, "jacobian2D_xmin", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", jacobian2D_xmin, adios_err)
  call adios_set_path (adios_handle, "jacobian2D_xmax", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", jacobian2D_xmax, adios_err)

  !local_dim = NDIM*NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX
  local_dim = size (jacobian2D_ymin)
  call adios_set_path (adios_handle, "jacobian2D_ymin", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", jacobian2D_ymin, adios_err)
  call adios_set_path (adios_handle, "jacobian2D_ymax", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", jacobian2D_ymax, adios_err)

  !local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_BOTTOM
  local_dim = size (jacobian2D_bottom)
  call adios_set_path (adios_handle, "jacobian2D_bottom", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", jacobian2D_bottom, adios_err)

  !local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_TOP
  local_dim = size (jacobian2D_top)
  call adios_set_path (adios_handle, "jacobian2D_top", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", jacobian2D_top, adios_err)

  !--- Reset the path to zero and perform the actual write to disk
  call adios_set_path (adios_handle, "", adios_err)
  call adios_close(adios_handle, adios_err)

  !---------------------------------------------------------
  !--- Attenuation arrays ----------------------------------
  !---------------------------------------------------------
  if(ATTENUATION) then
    outputname = trim(reg_name) // "attenuation.bp"
    write(group_name,"('SPECFEM3D_GLOBE_ATTENUATION_reg',i1)") iregion_code
    group_size_inc = 0
    call adios_declare_group(adios_group, group_name, &
        "", 0, adios_err)
    call adios_select_method(adios_group, "MPI", "", "", adios_err)

    !--- Define ADIOS variables -----------------------------
    call define_adios_double_scalar(adios_group, "T_c_source", "", &
        group_size_inc)

    local_dim = size(tau_s)
    call define_adios_global_double_1d_array(adios_group, "tau_s", &
        local_dim, group_size_inc)
    local_dim = size(tau_e_store)
    call define_adios_global_double_1d_array(adios_group, "tau_e_store", &
        local_dim, group_size_inc)
    local_dim = size(Qmu_store)
    call define_adios_global_double_1d_array(adios_group, "Qmu_store", &
        local_dim, group_size_inc)

    !--- Open an ADIOS handler to the restart file. ---------
    call adios_open (adios_handle, group_name, &
        outputname, "w", comm, adios_err);
    call adios_group_size (adios_handle, group_size_inc, &
                           adios_totalsize, adios_err)

    !--- Schedule writes for the previously defined ADIOS variables
    call adios_write(adios_handle, "T_c_source", T_c_source, adios_err)

    local_dim = size (tau_s)
    call adios_set_path (adios_handle, "tau_s", adios_err)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        local_dim, sizeprocs)
    call adios_write(adios_handle, "array", tau_s, adios_err)
    local_dim = size (tau_e_store)
    call adios_set_path (adios_handle, "tau_e_store", adios_err)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        local_dim, sizeprocs)
    call adios_write(adios_handle, "array", tau_e_store, adios_err)
    local_dim = size (Qmu_store)
    call adios_set_path (adios_handle, "Qmu_store", adios_err)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        local_dim, sizeprocs)
    call adios_write(adios_handle, "array", Qmu_store, adios_err)

    !--- Reset the path to zero and perform the actual write to disk
    call adios_set_path (adios_handle, "", adios_err)
    call adios_close(adios_handle, adios_err)
  endif

  !---------------------------------------------------------
  !--- dvp arrays ------------------------------------------
  !---------------------------------------------------------
  if(HETEROGEN_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
    outputname = trim(reg_name) // "dvp.bp"
    write(group_name,"('SPECFEM3D_GLOBE_DVP_reg',i1)") iregion_code
    group_size_inc = 0
    call adios_declare_group(adios_group, group_name, &
        "", 0, adios_err)
    call adios_select_method(adios_group, "MPI", "", "", adios_err)

    !--- Define ADIOS variables -----------------------------
    local_dim = size (dvpstore)
    call define_adios_global_real_1d_array(adios_group, "dvp", &
        local_dim, group_size_inc)
    !--- Open an ADIOS handler to the restart file. ---------
    call adios_open (adios_handle, group_name, &
        outputname, "w", comm, adios_err);
    call adios_group_size (adios_handle, group_size_inc, &
                           adios_totalsize, adios_err)
    call adios_set_path (adios_handle, "dvp", adios_err)
    !--- Schedule writes for the previously defined ADIOS variables
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        local_dim, sizeprocs)
    call adios_write(adios_handle, "array", dvpstore, adios_err)

    !--- Reset the path to zero and perform the actual write to disk
    call adios_set_path (adios_handle, "", adios_err)
    call adios_close(adios_handle, adios_err)
  endif

  !---------------------------------------------------------
  !--- mehsfiles arrays ------------------------------------
  !---------------------------------------------------------
  ! uncomment for vp & vs model storage
  if( SAVE_MESH_FILES ) then
    ! outputs model files in binary format
    if (ADIOS_FOR_SOLVER_MESHFILES) then
      call save_arrays_solver_meshfiles_adios(myrank,iregion_code, &
          reg_name, nspec)
    else
      call save_arrays_solver_meshfiles(myrank,nspec)
    endif
  endif

end subroutine save_arrays_solver_adios


!===============================================================================
!> \brief Save the meshfiles that will be used by the solver in an ADIOS format.
!!
!! \param myrank The MPI rank of the current process.
!! \param iregion_code Code of the region considered. See constant.h for details
!! \param reg_name Output file prefix with the name of the region included
!! \param nspec Number of GLL points per spectral elements
subroutine save_arrays_solver_meshfiles_adios(myrank, iregion_code, &
    reg_name, nspec)

  ! outputs model files in binary format
  use mpi
  use adios_write_mod
  use constants

  use meshfem3D_models_par,only: &
    TRANSVERSE_ISOTROPY,ATTENUATION,ATTENUATION_3D,ATTENUATION_1D_WITH_3D_STORAGE

  use create_regions_mesh_par2,only: &
    rhostore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
    Qmu_store, &
    prname

  implicit none

  integer :: myrank, nspec, iregion_code
  character(len=150) :: reg_name

  ! local parameters
  integer :: i,j,k,ispec,ier
  real(kind=CUSTOM_REAL) :: scaleval1,scaleval2
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: temp_store

  ! local parameters
  character(len=150) :: outputname, group_name
  integer :: ierr, sizeprocs, comm, local_dim
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, adios_handle, varid
  integer(kind=8)         :: adios_groupsize, adios_totalsize

  ! scaling factors to re-dimensionalize units
  scaleval1 = sngl( sqrt(PI*GRAV*RHOAV)*(R_EARTH/1000.0d0) )
  scaleval2 = sngl( RHOAV/1000.0d0 )

  call world_size(sizeprocs) ! TODO keep it in parameters
  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)

  ! isotropic model
  outputname = trim(reg_name) // "solver_meshfiles.bp"
  write(group_name,"('SPECFEM3D_GLOBE_solver_meshfiles_reg',i1)") iregion_code

  group_size_inc = 0

  call adios_declare_group(adios_group, group_name, &
      "", 0, adios_err)
  call adios_select_method(adios_group, "MPI", "", "", adios_err)

  !--- Define ADIOS variables -----------------------------
  !--- vp arrays -------------------------------------------
  local_dim = size (kappavstore)
  call define_adios_global_real_1d_array(adios_group, "vp", &
      local_dim, group_size_inc)
  !--- vs arrays -------------------------------------------
  local_dim = size (rhostore)
  call define_adios_global_real_1d_array(adios_group, "vs", &
      local_dim, group_size_inc)
  !--- rho arrays ------------------------------------------
  local_dim = size (rhostore)
  call define_adios_global_real_1d_array(adios_group, "rho", &
      local_dim, group_size_inc)
  ! transverse isotropic model
  if( TRANSVERSE_ISOTROPY ) then
    !--- vpv arrays ----------------------------------------
    local_dim = size (kappavstore)
    call define_adios_global_real_1d_array(adios_group, "vpv", &
        local_dim, group_size_inc)
    !--- vph arrays ----------------------------------------
    local_dim = size (kappavstore)
    call define_adios_global_real_1d_array(adios_group, "vph", &
        local_dim, group_size_inc)
    !--- vsv arrays ----------------------------------------
    local_dim = size (rhostore)
    call define_adios_global_real_1d_array(adios_group, "vsv", &
        local_dim, group_size_inc)
    !--- vsh arrays ----------------------------------------
    local_dim = size (rhostore)
    call define_adios_global_real_1d_array(adios_group, "vsh", &
        local_dim, group_size_inc)
    !--- eta arrays ----------------------------------------
    local_dim = size (eta_anisostore)
    call define_adios_global_real_1d_array(adios_group, "eta", &
        local_dim, group_size_inc)
  endif
  if( ATTENUATION ) then
    !--- Qmu arrays ----------------------------------------
    local_dim = NGLLX * NGLLY * NGLLZ * nspec
    call define_adios_global_real_1d_array(adios_group, "qmu", &
        local_dim, group_size_inc)
  endif

  !--- Open an ADIOS handler to the restart file. ---------
  call adios_open (adios_handle, group_name, &
      outputname, "w", comm, adios_err);
  call adios_group_size (adios_handle, group_size_inc, &
                         adios_totalsize, adios_err)

  !--- Schedule writes for the previously defined ADIOS variables
  !--- vp arrays -------------------------------------------
  local_dim = size (kappavstore)
  call adios_set_path (adios_handle, "vp", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", &
      sqrt( (kappavstore+4.*muvstore/3.)/rhostore )*scaleval1, &
      adios_err)
  !--- vs arrays -------------------------------------------
  local_dim = size (rhostore)
  call adios_set_path (adios_handle, "vs", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", &
      sqrt( muvstore/rhostore )*scaleval1, &
      adios_err)
  !--- rho arrays ------------------------------------------
  local_dim = size (rhostore)
  call adios_set_path (adios_handle, "rho", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", &
      rhostore *scaleval2, &
      adios_err)

  ! transverse isotropic model
  if( TRANSVERSE_ISOTROPY ) then
    !--- vpv arrays ----------------------------------------
    local_dim = size (kappavstore)
    call adios_set_path (adios_handle, "vpv", adios_err)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        local_dim, sizeprocs)
    call adios_write(adios_handle, "array", &
        sqrt( (kappavstore+4.*muvstore/3.)/rhostore )*scaleval1, &
        adios_err)
    !--- vph arrays ----------------------------------------
    local_dim = size (kappavstore)
    call adios_set_path (adios_handle, "vph", adios_err)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        local_dim, sizeprocs)
    call adios_write(adios_handle, "array", &
        sqrt( (kappahstore+4.*muhstore/3.)/rhostore )*scaleval1, &
        adios_err)
    !--- vsv arrays ----------------------------------------
    local_dim = size (rhostore)
    call adios_set_path (adios_handle, "vsv", adios_err)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        local_dim, sizeprocs)
    call adios_write(adios_handle, "array", &
        sqrt( muvstore/rhostore )*scaleval1, &
        adios_err)
    !--- vsh arrays ----------------------------------------
    local_dim = size (rhostore)
    call adios_set_path (adios_handle, "vsh", adios_err)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        local_dim, sizeprocs)
    call adios_write(adios_handle, "array", &
        sqrt( muhstore/rhostore )*scaleval1, &
        adios_err)
    !--- eta arrays ----------------------------------------
    local_dim = size (eta_anisostore)
    call adios_set_path (adios_handle, "eta", adios_err)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        local_dim, sizeprocs)
    call adios_write(adios_handle, "array", &
        eta_anisostore, &
        adios_err)
  endif ! TRANSVERSE_ISOTROPY

  ! shear attenuation
  if( ATTENUATION ) then
    !-------------------------------------------------------
    !--- Qmu arrays ----------------------------------------
    !-------------------------------------------------------
    ! saves Qmu_store to full custom_real array
    ! uses temporary array
    allocate(temp_store(NGLLX,NGLLY,NGLLZ,nspec))
    if (ATTENUATION_3D .or. ATTENUATION_1D_WITH_3D_STORAGE) then
      ! attenuation arrays are fully 3D
      if(CUSTOM_REAL == SIZE_REAL) then
        temp_store(:,:,:,:) = sngl(Qmu_store(:,:,:,:))
      else
        temp_store(:,:,:,:) = Qmu_store(:,:,:,:)
      endif
    else
      ! attenuation array dimensions: Q_mustore(1,1,1,nspec)
      do ispec = 1,nspec
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              ! distinguish between single and double precision for reals
              if(CUSTOM_REAL == SIZE_REAL) then
                temp_store(i,j,k,ispec) = sngl(Qmu_store(1,1,1,ispec))
              else
                temp_store(i,j,k,ispec) = Qmu_store(1,1,1,ispec)
              endif
            enddo
          enddo
        enddo
      enddo
    endif

    local_dim = NGLLX * NGLLY * NGLLZ * nspec
    call adios_set_path (adios_handle, "qmu", adios_err)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        local_dim, sizeprocs)
    call adios_write(adios_handle, "array", &
        temp_store, &
        adios_err)

    ! frees temporary memory
    deallocate(temp_store)
  endif ! ATTENUATION

  !--- Reset the path to zero and perform the actual write to disk
  call adios_set_path (adios_handle, "", adios_err)
  call adios_close(adios_handle, adios_err)

end subroutine save_arrays_solver_meshfiles_adios


!===============================================================================
!> \brief Save the arrays use by the solver for MPI communications.
!!
!! \param myrank The MPI rank of the current process.
!! \param iregion_code Code of the region considered. See constant.h for details
!! \param LOCAL_PATH The full path to the output directory
!! \param num_interfaces The number of interfaces between processors
!! \param max_nibool_interfaces
!! \param my_neighbours
!! \param nibool_interfaces
!! \param ibool_interfaces
!! \param nspec_inner Number of spectral elements in the inner core
!! \param nspec_outer Number of spectral elemetns in the outer core
!! \param num_phase_ispec
!! \param phase_ispec_inner
!! \param num_colors_inner Number of colors for GPU computing in the inner core.
!! \param num_colors_outer Number of colors for GPU computing in the outer core.
subroutine save_MPI_arrays_adios(myrank,iregion_code,LOCAL_PATH, &
   num_interfaces,max_nibool_interfaces, my_neighbours,nibool_interfaces, &
   ibool_interfaces, nspec_inner,nspec_outer, num_phase_ispec, &
   phase_ispec_inner, num_colors_outer,num_colors_inner, num_elem_colors)

  use mpi
  use adios_write_mod
  implicit none

  include "constants.h"

  integer :: iregion_code,myrank
  character(len=150) :: LOCAL_PATH
  ! MPI interfaces
  integer :: num_interfaces,max_nibool_interfaces
  integer, dimension(num_interfaces) :: my_neighbours
  integer, dimension(num_interfaces) :: nibool_interfaces
  integer, dimension(max_nibool_interfaces,num_interfaces) :: &
      ibool_interfaces
  ! inner/outer elements
  integer :: nspec_inner,nspec_outer
  integer :: num_phase_ispec
  integer,dimension(num_phase_ispec,2) :: phase_ispec_inner
  ! mesh coloring
  integer :: num_colors_outer,num_colors_inner
  integer, dimension(num_colors_outer + num_colors_inner) :: &
    num_elem_colors

  ! local parameters
  character(len=150) :: prname, outputname, group_name
  integer :: ierr, sizeprocs, comm, local_dim
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, adios_handle, varid
  integer(kind=8)         :: adios_groupsize, adios_totalsize

  ! create the name for the database of the current slide and region
  call create_name_database_adios(prname,iregion_code,LOCAL_PATH)

  outputname = trim(prname) // "solver_data_mpi.bp"
  write(group_name,"('SPECFEM3D_GLOBE_MPI_ARRAYS_reg',i1)") iregion_code
  call world_size(sizeprocs) ! TODO keep it in parameters
  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
  group_size_inc = 0
  call adios_declare_group(adios_group, group_name, &
      "", 0, adios_err)
  call adios_select_method(adios_group, "MPI", "", "", adios_err)

  !--- Define ADIOS variables -----------------------------
  !! MPI interfaces
  call define_adios_integer_scalar (adios_group, "num_interfaces", "", &
      group_size_inc)
  if( num_interfaces > 0 ) then
    call define_adios_integer_scalar(adios_group, "max_nibool_interfaces", &
        "", group_size_inc)
    call define_adios_global_integer_1d_array(adios_group, "my_neighbours", &
        num_interfaces, group_size_inc)
    call define_adios_global_integer_1d_array(adios_group, "nibool_interfaces",&
        num_interfaces, group_size_inc)
    local_dim = max_nibool_interfaces*num_interfaces
    call define_adios_global_integer_1d_array(adios_group, "ibool_interfaces", &
        local_dim, group_size_inc)
  endif

  ! inner/outer elements
  call define_adios_integer_scalar (adios_group, "nspec_inner", "", &
      group_size_inc)
  call define_adios_integer_scalar (adios_group, "nspec_outer", "", &
      group_size_inc)
  call define_adios_integer_scalar (adios_group, "num_phase_ispec", "", &
      group_size_inc)
  if(num_phase_ispec > 0 ) then
    local_dim = num_phase_ispec * 2
    call define_adios_global_integer_1d_array(adios_group, "phase_ispec_inner", &
        local_dim, group_size_inc)
  endif

  ! mesh coloring
  if( USE_MESH_COLORING_GPU ) then
    call define_adios_integer_scalar (adios_group, "num_colors_outer", "", &
        group_size_inc)
    call define_adios_integer_scalar (adios_group, "num_colors_inner", "", &
        group_size_inc)
    call define_adios_global_integer_1d_array(adios_group, "num_elem_colors", &
        num_colors_outer + num_colors_inner, group_size_inc)
  endif

  !--- Open an ADIOS handler to the restart file. ---------
  call adios_open (adios_handle, group_name, &
      outputname, "w", comm, adios_err);
  call adios_group_size (adios_handle, group_size_inc, &
                         adios_totalsize, adios_err)

  !--- Schedule writes for the previously defined ADIOS variables
  ! MPI interfaces
  call adios_write(adios_handle, "num_interfaces", num_interfaces, adios_err)
  if( num_interfaces > 0 ) then
    call adios_write(adios_handle, "max_nibool_interfaces", &
        max_nibool_interfaces, adios_err)

    local_dim = num_interfaces

    call adios_set_path (adios_handle, "my_neighbours", adios_err)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        local_dim, sizeprocs)
    call adios_write(adios_handle, "array", my_neighbours, adios_err)

    call adios_set_path (adios_handle, "nibool_interfaces", adios_err)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        local_dim, sizeprocs)
    call adios_write(adios_handle, "array", nibool_interfaces, adios_err)

    local_dim = max_nibool_interfaces * num_interfaces

    call adios_set_path (adios_handle, "ibool_interfaces", adios_err)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        local_dim, sizeprocs)
    call adios_write(adios_handle, "array", &
        ibool_interfaces, adios_err)
    call adios_set_path (adios_handle, "", adios_err)
  endif

  ! inner/outer elements
  call adios_write(adios_handle, "nspec_inner", nspec_inner, adios_err)
  call adios_write(adios_handle, "nspec_outer", nspec_outer, adios_err)
  call adios_write(adios_handle, "num_phase_ispec", num_phase_ispec, adios_err)

  if(num_phase_ispec > 0 ) then
    local_dim = num_phase_ispec * 2
    call adios_set_path (adios_handle, "phase_ispec_inner", adios_err)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        local_dim, sizeprocs)
    call adios_write(adios_handle, "array", &
        phase_ispec_inner, adios_err)
    call adios_set_path (adios_handle, "", adios_err)
  endif

  ! mesh coloring
  if( USE_MESH_COLORING_GPU ) then
    call adios_write(adios_handle, "num_colors_outer", nspec_inner, adios_err)
    call adios_write(adios_handle, "num_colors_inner", nspec_inner, adios_err)
    local_dim = num_colors_outer + num_colors_inner
    call adios_set_path (adios_handle, "num_elem_colors", adios_err)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        local_dim, sizeprocs)
    call adios_write(adios_handle, "array", &
        num_elem_colors, adios_err)
    call adios_set_path (adios_handle, "", adios_err)
  endif

  !--- Reset the path to zero and perform the actual write to disk
  call adios_close(adios_handle, adios_err)

end subroutine save_MPI_arrays_adios


!===============================================================================
!> \brief Write boundary conditions (MOHO, 400, 600) to a single ADIOS file.
subroutine save_arrays_solver_boundary_adios()

! saves arrays for boundaries such as MOHO, 400 and 670 discontinuities
  use mpi

  use meshfem3d_par,only: &
    myrank, LOCAL_PATH

  use meshfem3D_models_par,only: &
    HONOR_1D_SPHERICAL_MOHO
    !SAVE_BOUNDARY_MESH,HONOR_1D_SPHERICAL_MOHO,SUPPRESS_CRUSTAL_MESH

  use create_regions_mesh_par2, only: &
    NSPEC2D_MOHO, NSPEC2D_400, NSPEC2D_670, &
    ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot, &
    ibelm_670_top,ibelm_670_bot,normal_moho,normal_400,normal_670, &
    ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,ispec2D_400_bot, &
    ispec2D_670_top,ispec2D_670_bot, &
    prname

  implicit none
  include "constants.h"

  ! local parameters
  ! local parameters
  character(len=150) :: outputname, group_name
  integer :: ierr, sizeprocs, comm, local_dim
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, adios_handle, varid
  integer(kind=8)         :: adios_groupsize, adios_totalsize

  ! first check the number of surface elements are the same for Moho, 400, 670
  if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
    if (ispec2D_moho_top /= NSPEC2D_MOHO .or. ispec2D_moho_bot /= NSPEC2D_MOHO) &
           call exit_mpi(myrank, 'Not the same number of Moho surface elements')
  endif
  if (ispec2D_400_top /= NSPEC2D_400 .or. ispec2D_400_bot /= NSPEC2D_400) &
           call exit_mpi(myrank,'Not the same number of 400 surface elements')
  if (ispec2D_670_top /= NSPEC2D_670 .or. ispec2D_670_bot /= NSPEC2D_670) &
           call exit_mpi(myrank,'Not the same number of 670 surface elements')

  outputname = trim(LOCAL_PATH) // "/boundary_disc.bp"
  group_name = "SPECFEM3D_GLOBE_BOUNDARY_DISC"
  call world_size(sizeprocs) ! TODO keep it in parameters
  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
  group_size_inc = 0
  call adios_declare_group(adios_group, group_name, &
      "", 0, adios_err)
  call adios_select_method(adios_group, "MPI", "", "", adios_err)

  !--- Define ADIOS variables -----------------------------
  call define_adios_integer_scalar (adios_group, "NSPEC2D_MOHO", "", &
      group_size_inc)
  call define_adios_integer_scalar (adios_group, "NSPEC2D_400", "", &
      group_size_inc)
  call define_adios_integer_scalar (adios_group, "NSPEC2D_670", "", &
      group_size_inc)

  local_dim = NSPEC2D_MOHO
  call define_adios_global_integer_1d_array(adios_group, "ibelm_moho_top", &
        local_dim, group_size_inc)
  call define_adios_global_integer_1d_array(adios_group, "ibelm_moho_bot", &
        local_dim, group_size_inc)
  local_dim = NSPEC2D_400
  call define_adios_global_integer_1d_array(adios_group, "ibelm_400_top", &
        local_dim, group_size_inc)
  call define_adios_global_integer_1d_array(adios_group, "ibelm_400_bot", &
        local_dim, group_size_inc)
  local_dim = NSPEC2D_670
  call define_adios_global_integer_1d_array(adios_group, "ibelm_670_top", &
        local_dim, group_size_inc)
  call define_adios_global_integer_1d_array(adios_group, "ibelm_670_bot", &
        local_dim, group_size_inc)
  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_MOHO
  call define_adios_global_real_1d_array(adios_group, "normal_moho", &
        local_dim, group_size_inc)
  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_400
  call define_adios_global_real_1d_array(adios_group, "normal_400", &
        local_dim, group_size_inc)
  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_670
  call define_adios_global_real_1d_array(adios_group, "normal_670", &
        local_dim, group_size_inc)

  !--- Open an ADIOS handler to the restart file. ---------
  call adios_open (adios_handle, group_name, &
      outputname, "w", comm, adios_err);
  call adios_group_size (adios_handle, group_size_inc, &
                         adios_totalsize, adios_err)

  !--- Schedule writes for the previously defined ADIOS variables
  call adios_write(adios_handle, "NSPEC2D_MOHO", NSPEC2D_MOHO, adios_err)
  call adios_write(adios_handle, "NSPEC2D_400", NSPEC2D_400, adios_err)
  call adios_write(adios_handle, "NSPEC2D_670", NSPEC2D_670, adios_err)

  local_dim = NSPEC2D_MOHO
  call adios_set_path (adios_handle, "ibelm_moho_top", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", ibelm_moho_top, adios_err)
  call adios_set_path (adios_handle, "ibelm_moho_bot", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", ibelm_moho_bot, adios_err)

  local_dim = NSPEC2D_400
  call adios_set_path (adios_handle, "ibelm_400_top", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", ibelm_400_top, adios_err)
  call adios_set_path (adios_handle, "ibelm_400_bot", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", ibelm_400_bot, adios_err)

  local_dim = NSPEC2D_670
  call adios_set_path (adios_handle, "ibelm_670_top", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", ibelm_670_top, adios_err)
  call adios_set_path (adios_handle, "ibelm_670_bot", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", ibelm_670_bot, adios_err)

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_MOHO
  call adios_set_path (adios_handle, "normal_moho", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", normal_moho, adios_err)

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_400
  call adios_set_path (adios_handle, "normal_400", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", normal_400, adios_err)

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_670
  call adios_set_path (adios_handle, "normal_670", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", normal_670, adios_err)

  !--- Reset the path to zero and perform the actual write to disk
  call adios_set_path (adios_handle, "", adios_err)
  call adios_close(adios_handle, adios_err)
end subroutine save_arrays_solver_boundary_adios

!-------------------------------------------------------------------------------
!> Write local, global and offset dimensions to ADIOS
!! \param adios_handle Handle to the adios file
!! \param local_dim Number of elements to be written by one process
!! \param sizeprocs Number of MPI processes
subroutine write_1D_global_array_adios_dims(adios_handle, myrank, &
    local_dim, sizeprocs)
  use adios_write_mod

  implicit none

  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: sizeprocs, local_dim, myrank

  integer :: adios_err

  call adios_write(adios_handle, "local_dim", local_dim, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_write(adios_handle, "global_dim", local_dim*sizeprocs, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_write(adios_handle, "offset", local_dim*myrank, adios_err)
  call check_adios_err(myrank,adios_err)
end subroutine write_1D_global_array_adios_dims

