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

#include "config.fh"

!===============================================================================
!> \brief Main routine to save the arrays from the mesher to the solver with the
!!        help of ADIOS
!! \param myrank The MPI rank of the current process
!! \param nspec Number of GLL points per element
!! \param nglob Number of mesh points
!! \param idoubling Array of information on every mesh point
!! \param ibool Array of information on every mesh point
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
  subroutine save_arrays_solver_adios(idoubling,ibool,xstore,ystore,zstore, &
                                      NSPEC2DMAX_XMIN_XMAX, NSPEC2DMAX_YMIN_YMAX, &
                                      NSPEC2D_TOP,NSPEC2D_BOTTOM)

  use constants

  use shared_parameters, only: ATT_F_C_SOURCE

  use meshfem_models_par, only: &
    OCEANS,TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
    ANISOTROPIC_INNER_CORE,ATTENUATION

  use meshfem_par, only: &
    NCHUNKS,ABSORBING_CONDITIONS,LOCAL_PATH, &
    ROTATION,EXACT_MASS_MATRIX_FOR_ROTATION, &
    nspec,nglob,iregion_code

  use regions_mesh_par2, only: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
    gammaxstore,gammaystore,gammazstore, &
    rhostore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
    mu0store, &
    rmassx,rmassy,rmassz,rmass_ocean_load, &
    b_rmassx,b_rmassy, &
    ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
    normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top, &
    jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax, &
    jacobian2D_bottom,jacobian2D_top, &
    rho_vp,rho_vs, &
    nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
    ispec_is_tiso, &
    tau_s_store,tau_e_store,Qmu_store, &
    nspec_ani, nspec_stacey, nglob_xy, nglob_oceans ! prname,

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! doubling mesh flag
  integer, dimension(nspec) :: idoubling
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  ! arrays with the mesh in double precision
  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  ! boundary parameters locator
  integer,intent(in) :: NSPEC2D_TOP,NSPEC2D_BOTTOM, NSPEC2DMAX_XMIN_XMAX, NSPEC2DMAX_YMIN_YMAX

  ! local parameters
  integer :: NSPEC2D_TOP_wmax,NSPEC2D_BOTTOM_wmax, NSPEC2DMAX_XMIN_XMAX_wmax, NSPEC2DMAX_YMIN_YMAX_wmax
  integer, parameter :: num_ints_to_reduce = 4
  integer, dimension(num_ints_to_reduce) :: ints_to_reduce

  integer :: i,j,k,ispec,iglob,ier
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: tmp_array_x, tmp_array_y, tmp_array_z
  double precision :: f_c_source

  ! local parameters
  character(len=MAX_STRING_LEN) :: outputname, group_name ! reg_name
  integer(kind=8) :: local_dim
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  character(len=128)      :: region_name, region_name_scalar

  !--- Save the number of region written. Open the file in "w" mode if 0, else
  !    in "a"  mode
  integer, save :: num_regions_written = 0

  ! user output
  if (myrank == 0) then
#if defined(USE_ADIOS)
    write(IMAIN,*) '    solver   in ADIOS1 file format'
#elif defined(USE_ADIOS2)
    write(IMAIN,*) '    solver   in ADIOS2 file format'
#endif
    call flush_IMAIN()
  endif

  ! mesh topology

  ! mesh arrays used in the solver to locate source and receivers
  ! and for anisotropy and gravity, save in single precision
  ! use tmp_array for temporary storage to perform conversion
  allocate(tmp_array_x(nglob),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary x array for mesh topology')
  allocate(tmp_array_y(nglob),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary y array for mesh topology')
  allocate(tmp_array_z(nglob),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating temporary z array for mesh topology')

  ! note: we keep these temporary arrays until after the perform/close/end_step call is done,
  !       since the write_adios_** calls might be in deferred mode.

  !--- x coordinate
  tmp_array_x(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! distinguish between single and double precision for reals
          tmp_array_x(iglob) = real(xstore(i,j,k,ispec), kind=CUSTOM_REAL)
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
          tmp_array_y(iglob) = real(ystore(i,j,k,ispec), kind=CUSTOM_REAL)
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
          tmp_array_z(iglob) = real(zstore(i,j,k,ispec), kind=CUSTOM_REAL)
        enddo
      enddo
    enddo
  enddo

  ! create a prefix for the file name such as LOCAL_PATH/regX_
  !call create_name_database_adios(reg_name,iregion_code,LOCAL_PATH)

  write(region_name,"('reg',i1, '/')") iregion_code
  write(region_name_scalar,"('reg',i1)") iregion_code

  !---------------------------------------------------------
  !--- Solver data arrays ----------------------------------
  !---------------------------------------------------------

  ! save arrays for the solver to run.
  write(group_name,"('SPECFEM3D_GLOBE_ARRAYS_SOLVER_reg',i1)") iregion_code

  ! set the adios group size to 0 before incremented by calls to
  ! helpers functions.
  call init_adios_group(myadios_group,group_name)

  !--- Define ADIOS variables -----------------------------
  ! save nspec and nglob, to be used in combine_paraview_data
  group_size_inc = 0
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec))
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nglob))

  local_dim = nglob
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, "x_global", tmp_array_x)
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, "y_global", tmp_array_y)
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, "z_global", tmp_array_z)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(xstore))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ystore))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(zstore))

  local_dim = NGLLX * NGLLY * NGLLZ * nspec
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibool))

  local_dim = nspec
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(idoubling))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ispec_is_tiso))

  ! local GLL points
  local_dim = NGLLX * NGLLY * NGLLZ * nspec
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(xixstore))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(xiystore))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(xizstore))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(etaxstore))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(etaystore))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(etazstore))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(gammaxstore))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(gammaystore))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(gammazstore))

  local_dim = NGLLX * NGLLY * NGLLZ * nspec
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(rhostore))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(kappavstore))
  if (iregion_code /= IREGION_OUTER_CORE) then
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(muvstore))
  endif

  ! solid regions
  select case(iregion_code)
  case (IREGION_CRUST_MANTLE)
    ! save anisotropy in the mantle only
    if (ANISOTROPIC_3D_MANTLE) then
      local_dim = NGLLX * NGLLY * NGLLZ * nspec_ani
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c11store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c12store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c13store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c14store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c15store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c16store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c22store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c23store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c24store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c25store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c26store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c33store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c34store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c35store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c36store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c44store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c45store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c46store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c55store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c56store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c66store))
    else
      if (TRANSVERSE_ISOTROPY) then
        local_dim = NGLLX * NGLLY * NGLLZ * nspec
        call define_adios_global_array1D(myadios_group,group_size_inc,local_dim,region_name,STRINGIFY_VAR(kappahstore))
        call define_adios_global_array1D(myadios_group,group_size_inc,local_dim,region_name,STRINGIFY_VAR(muhstore))
        call define_adios_global_array1D(myadios_group,group_size_inc,local_dim,region_name,STRINGIFY_VAR(eta_anisostore))
      endif
    endif
    ! for azimuthal aniso kernels
    local_dim = NGLLX * NGLLY * NGLLZ * nspec
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(mu0store))

  case (IREGION_INNER_CORE)
    !   save anisotropy in the inner core only
    if (ANISOTROPIC_INNER_CORE) then
      local_dim = NGLLX * NGLLY * NGLLZ * nspec_ani
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c11store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c12store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c13store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c33store))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c44store))
    endif
  end select

  ! Stacey
  if (ABSORBING_CONDITIONS) then
    local_dim = NGLLX * NGLLY * NGLLZ * nspec_stacey
    if (iregion_code == IREGION_CRUST_MANTLE) then
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(rho_vp))
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(rho_vs))
    else if (iregion_code == IREGION_OUTER_CORE) then
      call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(rho_vp))
    endif
  endif

  ! mass matrices
  if ((NCHUNKS /= 6 .and. ABSORBING_CONDITIONS .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
     (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
     (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION .and. iregion_code == IREGION_INNER_CORE)) then
    local_dim = nglob_xy
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(rmassx))
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(rmassy))
  endif
  local_dim = nglob
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(rmassz))

 ! mass matrices for backward simulation when ROTATION is .true.
  if (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
    if (iregion_code == IREGION_CRUST_MANTLE .or. iregion_code == IREGION_INNER_CORE) then
      local_dim = nglob_xy
      call define_adios_global_array1D(myadios_group,group_size_inc, local_dim,region_name, STRINGIFY_VAR(b_rmassx))
      call define_adios_global_array1D(myadios_group,group_size_inc, local_dim,region_name, STRINGIFY_VAR(b_rmassy))
    endif
  endif

  ! additional ocean load mass matrix if oceans and if we are in the crust
  if (OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) then
    local_dim = nglob_oceans
    call define_adios_global_array1D(myadios_group,group_size_inc,local_dim,region_name,STRINGIFY_VAR(rmass_ocean_load))
  endif


  !--- Open an ADIOS handler to the restart file. ---------
  outputname = get_adios_filename(trim(LOCAL_PATH) // "/solver_data")
  ! user output
  if (myrank == 0) write(IMAIN,*) '    saving arrays in ADIOS file: ',trim(outputname)

  if (num_regions_written == 0) then
    ! opens file for writing
    call open_file_adios_write(myadios_file,myadios_group,outputname,group_name)
  else
    ! opens file for writing in append mode
    call open_file_adios_write_append(myadios_file,myadios_group,outputname,group_name)
  endif

  ! note: adios2 increases step numbers on variables when appending to a file.
  !       this can lead to issues when reading back values for the next regions, for example, reg2/nspec
  !       to work-around this, we explicitly call begin_step() and end_step() for writing out region1/2/3 data
  call write_adios_begin_step(myadios_file)

  ! sets group size
  call set_adios_group_size(myadios_file,group_size_inc)

  !--- Schedule writes for the previously defined ADIOS variables
  ! save nspec and nglob, to be used in combine_paraview_data
  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "nspec",nspec)
  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "nglob",nglob)

  local_dim = nglob
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // "x_global", tmp_array_x)
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // "y_global", tmp_array_y)
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // "z_global", tmp_array_z)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(xstore))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(ystore))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(zstore))

  local_dim = NGLLX * NGLLY * NGLLZ * nspec
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(ibool))

  local_dim = nspec
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                  trim(region_name) // STRINGIFY_VAR(idoubling))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(ispec_is_tiso))

  local_dim = NGLLX * NGLLY * NGLLZ * nspec
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(xixstore))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(xiystore))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(xizstore))

  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(etaxstore))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(etaystore))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(etazstore))

  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(gammaxstore))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(gammaystore))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(gammazstore))

  local_dim = NGLLX * NGLLY * NGLLZ * nspec
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(rhostore))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(kappavstore))
  if (iregion_code /= IREGION_OUTER_CORE) then
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // STRINGIFY_VAR(muvstore))
  endif

  ! solid regions
  select case(iregion_code)
  case (IREGION_CRUST_MANTLE)
    ! save anisotropy in the mantle only
    if (ANISOTROPIC_3D_MANTLE) then
      local_dim = NGLLX * NGLLY * NGLLZ * nspec_ani
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c11store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c12store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c13store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c14store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c15store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c16store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c22store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c23store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c24store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c25store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c26store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c33store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c34store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c35store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c36store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c44store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c45store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c46store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c55store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c56store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c66store))
    else
      if (TRANSVERSE_ISOTROPY) then
        local_dim = NGLLX * NGLLY * NGLLZ * nspec
        call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                    trim(region_name) // STRINGIFY_VAR(kappahstore))
        call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                    trim(region_name) // STRINGIFY_VAR(muhstore))
        call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                    trim(region_name) // STRINGIFY_VAR(eta_anisostore))
      endif
    endif
    local_dim = NGLLX * NGLLY * NGLLZ * nspec
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                 trim(region_name) // STRINGIFY_VAR(mu0store))

  case (IREGION_INNER_CORE)
    ! save anisotropy in the inner core only
    if (ANISOTROPIC_INNER_CORE) then
      local_dim = NGLLX * NGLLY * NGLLZ * nspec_ani
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c11store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c12store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c13store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c33store))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c44store))
    endif
  end select

  if (ABSORBING_CONDITIONS) then
    local_dim = NGLLX * NGLLY * NGLLZ * nspec_stacey
    if (iregion_code == IREGION_CRUST_MANTLE) then
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // STRINGIFY_VAR(rho_vp))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // STRINGIFY_VAR(rho_vs))

    else if (iregion_code == IREGION_OUTER_CORE) then
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // STRINGIFY_VAR(rho_vp))
    endif
  endif

  ! mass matrices
  if ((NCHUNKS /= 6 .and. ABSORBING_CONDITIONS .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
     (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
     (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION .and. iregion_code == IREGION_INNER_CORE)) then
    local_dim = nglob_xy
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // STRINGIFY_VAR(rmassx))
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // STRINGIFY_VAR(rmassy))
  endif

  local_dim = nglob
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(rmassz))

  ! mass matrices for backward simulation when ROTATION is .true.
  if (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
    if (iregion_code == IREGION_CRUST_MANTLE .or. iregion_code == IREGION_INNER_CORE) then
      local_dim = nglob_xy
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                       trim(region_name) // STRINGIFY_VAR(b_rmassx))
      call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                       trim(region_name) // STRINGIFY_VAR(b_rmassy))
    endif
  endif

  ! additional ocean load mass matrix if oceans and if we are in the crust
  if (OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) then
    local_dim = nglob_oceans
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // STRINGIFY_VAR(rmass_ocean_load))
    ! checks
    if (minval(rmass_ocean_load) <= 0._CUSTOM_REAL) &
        call exit_MPI(myrank,'negative mass matrix term for the oceans')
  endif

  ! end step to indicate output is completed. ADIOS2 can do I/O
  call write_adios_end_step(myadios_file)

  !--- perform the actual write to disk
  ! Reset the path to its original value to avoid bugs.
  call write_adios_perform(myadios_file)
  ! closes file
  call close_file_adios(myadios_file)

  ! makes sure all done
  call synchronize_all()

  ! Clean the temporary arrays containing the node information
  deallocate(tmp_array_x)
  deallocate(tmp_array_y)
  deallocate(tmp_array_z)

  !---------------------------------------------------------
  !--- Boundary arrays -------------------------------------
  !---------------------------------------------------------

  ints_to_reduce(1) = NSPEC2D_TOP
  ints_to_reduce(2) = NSPEC2D_BOTTOM
  ints_to_reduce(3) = NSPEC2DMAX_XMIN_XMAX
  ints_to_reduce(4) = NSPEC2DMAX_YMIN_YMAX

  call max_all_all_veci(ints_to_reduce,num_ints_to_reduce)

  NSPEC2D_TOP_wmax          = ints_to_reduce(1)
  NSPEC2D_BOTTOM_wmax       = ints_to_reduce(2)
  NSPEC2DMAX_XMIN_XMAX_wmax = ints_to_reduce(3)
  NSPEC2DMAX_YMIN_YMAX_wmax = ints_to_reduce(4)

  ! debug daniel
  !call synchronize_all()
  !print *,'debug: ',myrank,'nspec2d top      :',NSPEC2D_TOP,NSPEC2D_TOP_wmax
  !print *,'debug: ',myrank,'nspec2d bottom   :',NSPEC2D_BOTTOM,NSPEC2D_BOTTOM_wmax
  !print *,'debug: ',myrank,'nspec2d xmin_xmax:',NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_XMIN_XMAX_wmax
  !print *,'debug: ',myrank,'nspec2d ymin_ymax:',NSPEC2DMAX_YMIN_YMAX,NSPEC2DMAX_YMIN_YMAX_wmax
  !call synchronize_all()

  ! checks
  if (NSPEC2D_TOP /= NSPEC2D_TOP_wmax .or. &
      NSPEC2D_BOTTOM /= NSPEC2D_BOTTOM_wmax .or. &
      NSPEC2DMAX_XMIN_XMAX /= NSPEC2DMAX_XMIN_XMAX_wmax .or. &
      NSPEC2DMAX_YMIN_YMAX /= NSPEC2DMAX_YMIN_YMAX_wmax) then
    print *,myrank,'Error nspec2d for coupling surfaces'
    call exit_mpi(myrank,'Error nspec2d for coupling surfaces in adios saved file')
  endif

  ! save boundary arrays in ADIOS files
  write(group_name,"('SPECFEM3D_GLOBE_BOUNDARY_reg',i1)") iregion_code
  ! set the adios group size to 0 before incremented by calls to
  ! helpers functions.
  call init_adios_group(myadios_group,group_name)

  !--- Define ADIOS variables -----------------------------
  group_size_inc = 0
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec2D_xmin))
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec2D_xmax))
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec2D_ymin))
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec2D_ymax))
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(NSPEC2D_BOTTOM))
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(NSPEC2D_TOP))

  ! boundary elements
  local_dim = NSPEC2DMAX_XMIN_XMAX
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_xmin))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_xmax))

  local_dim = NSPEC2DMAX_YMIN_YMAX
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_ymin))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_ymax))

  local_dim = NSPEC2D_BOTTOM
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_bottom))

  local_dim = NSPEC2D_TOP
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_top))

  ! normals
  local_dim = NDIM*NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(normal_xmin))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(normal_xmax))

  local_dim = NDIM*NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(normal_ymin))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(normal_ymax))

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_BOTTOM
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(normal_bottom))

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_TOP
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(normal_top))

  ! Jacobians
  local_dim = NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(jacobian2D_xmin))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(jacobian2D_xmax))

  local_dim = NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(jacobian2D_ymin))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(jacobian2D_ymax))

  local_dim = NGLLX*NGLLY*NSPEC2D_BOTTOM
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(jacobian2D_bottom))

  local_dim = NGLLX*NGLLY*NSPEC2D_TOP
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(jacobian2D_top))

  !--- Open an ADIOS handler to the restart file. ---------
  outputname = get_adios_filename(trim(LOCAL_PATH) // "/boundary")
  ! user output
  if (myrank == 0) write(IMAIN,*) '    saving arrays in ADIOS file: ',trim(outputname)

  if (num_regions_written == 0) then
    ! opens file for writing
    call open_file_adios_write(myadios_file,myadios_group,outputname,group_name)
  else
    ! opens file for writing in append mode
    call open_file_adios_write_append(myadios_file,myadios_group,outputname,group_name)
  endif

  call set_adios_group_size(myadios_file,group_size_inc)

  !--- Schedule writes for the previously defined ADIOS variables
  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "nspec2D_xmin", nspec2D_xmin)
  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "nspec2D_xmax", nspec2D_xmax)

  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "nspec2D_ymin", nspec2D_ymin)
  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "nspec2D_ymax", nspec2D_ymax)

  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "NSPEC2D_BOTTOM", NSPEC2D_BOTTOM)
  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "NSPEC2D_TOP", NSPEC2D_TOP)

  ! boundary elements
  local_dim = NSPEC2DMAX_XMIN_XMAX
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                 trim(region_name) // STRINGIFY_VAR(ibelm_xmin))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                 trim(region_name) // STRINGIFY_VAR(ibelm_xmax))

  local_dim = NSPEC2DMAX_YMIN_YMAX
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                 trim(region_name) // STRINGIFY_VAR(ibelm_ymin))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                 trim(region_name) // STRINGIFY_VAR(ibelm_ymax))

  local_dim = NSPEC2D_BOTTOM
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                               trim(region_name) // STRINGIFY_VAR(ibelm_bottom))

  local_dim = NSPEC2D_TOP
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                  trim(region_name) // STRINGIFY_VAR(ibelm_top))

  ! normals
  local_dim = NDIM*NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                trim(region_name) // STRINGIFY_VAR(normal_xmin))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                trim(region_name) // STRINGIFY_VAR(normal_xmax))

  local_dim = NDIM*NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                trim(region_name) // STRINGIFY_VAR(normal_ymin))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                trim(region_name) // STRINGIFY_VAR(normal_ymax))

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_BOTTOM
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                              trim(region_name) // STRINGIFY_VAR(normal_bottom))

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_TOP
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                 trim(region_name) // STRINGIFY_VAR(normal_top))

  ! Jacobians
  local_dim = NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                            trim(region_name) // STRINGIFY_VAR(jacobian2D_xmin))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                            trim(region_name) // STRINGIFY_VAR(jacobian2D_xmax))

  local_dim = NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                            trim(region_name) // STRINGIFY_VAR(jacobian2D_ymin))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                            trim(region_name) // STRINGIFY_VAR(jacobian2D_ymax))

  local_dim = NGLLX*NGLLY*NSPEC2D_BOTTOM
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                          trim(region_name) // STRINGIFY_VAR(jacobian2D_bottom))

  local_dim = NGLLX*NGLLY*NSPEC2D_TOP
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                             trim(region_name) // STRINGIFY_VAR(jacobian2D_top))

  !--- Reset the path to zero and perform the actual write to disk
  call write_adios_perform(myadios_file)
  ! closes file
  call close_file_adios(myadios_file)

  !---------------------------------------------------------
  !--- Attenuation arrays ----------------------------------
  !---------------------------------------------------------
  if (ATTENUATION) then
    ! attenuation center frequency
    f_c_source = ATT_F_C_SOURCE

    write(group_name,"('SPECFEM3D_GLOBE_ATTENUATION_reg',i1)") iregion_code

    ! set the adios group size to 0 before incremented by calls to helpers functions.
    group_size_inc = 0
    call init_adios_group(myadios_group,group_name)

    !--- Define ADIOS variables -----------------------------
    call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(f_c_source))

    local_dim = size(tau_s_store)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(tau_s_store))
    local_dim = size(tau_e_store)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(tau_e_store))
    local_dim = size(Qmu_store)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(Qmu_store))

    !--- Open an ADIOS handler to the restart file. ---------
    outputname = get_adios_filename(trim(LOCAL_PATH) // "/attenuation")
    ! user output
    if (myrank == 0) write(IMAIN,*) '    saving arrays in ADIOS file: ',trim(outputname)

    if (num_regions_written == 0) then
      ! opens file for writing
      call open_file_adios_write(myadios_file,myadios_group,outputname,group_name)
    else
      ! opens file for writing in append mode
      call open_file_adios_write_append(myadios_file,myadios_group,outputname,group_name)
    endif

    call set_adios_group_size(myadios_file,group_size_inc)

    !--- Schedule writes for the previously defined ADIOS variables
    call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // STRINGIFY_VAR(f_c_source))

    local_dim = size (tau_s_store)
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // STRINGIFY_VAR(tau_s_store))
    local_dim = size (tau_e_store)
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                trim(region_name) // STRINGIFY_VAR(tau_e_store))
    local_dim = size (Qmu_store)
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                  trim(region_name) // STRINGIFY_VAR(Qmu_store))

    !--- Reset the path to zero and perform the actual write to disk
    call write_adios_perform(myadios_file)
    ! closes file
    call close_file_adios(myadios_file)

  endif

  num_regions_written = num_regions_written + 1

  end subroutine save_arrays_solver_adios


!===============================================================================
!> \brief Save the arrays use by the solver for MPI communications.
!!
!! \param myrank The MPI rank of the current process.
!! \param iregion_code Code of the region considered. See constant.h for details
!! \param LOCAL_PATH The full path to the output directory
!! \param num_interfaces The number of interfaces between processors
!! \param max_nibool_interfaces
!! \param my_neighbors
!! \param nibool_interfaces
!! \param ibool_interfaces
!! \param nspec_inner Number of spectral elements in the inner core
!! \param nspec_outer Number of spectral elements in the outer core
!! \param num_phase_ispec
!! \param phase_ispec_inner
!! \param num_colors_inner Number of colors for GPU computing in the inner core.
!! \param num_colors_outer Number of colors for GPU computing in the outer core.

  subroutine save_mpi_arrays_adios(iregion_code,LOCAL_PATH, &
                                   num_interfaces,max_nibool_interfaces, &
                                   my_neighbors,nibool_interfaces,ibool_interfaces, &
                                   nspec_inner,nspec_outer, &
                                   num_phase_ispec,phase_ispec_inner, &
                                   num_colors_outer,num_colors_inner,num_elem_colors)

  use constants

  use adios_helpers_mod
  use manager_adios

  implicit none

  integer, intent(in) :: iregion_code
  character(len=MAX_STRING_LEN), intent(in) :: LOCAL_PATH
  ! MPI interfaces
  integer, intent(in) :: num_interfaces,max_nibool_interfaces
  integer, dimension(num_interfaces), intent(in) :: my_neighbors
  integer, dimension(num_interfaces), intent(in) :: nibool_interfaces
  integer, dimension(max_nibool_interfaces,num_interfaces), intent(in) :: ibool_interfaces
  ! inner/outer elements
  integer, intent(in) :: nspec_inner,nspec_outer
  integer, intent(in) :: num_phase_ispec
  integer,dimension(num_phase_ispec,2), intent(in) :: phase_ispec_inner
  ! mesh coloring
  integer, intent(in) :: num_colors_outer,num_colors_inner
  integer, dimension(num_colors_outer + num_colors_inner), intent(in) :: num_elem_colors

  ! local parameters
  character(len=MAX_STRING_LEN) :: outputname, group_name ! prname,
  integer(kind=8) :: local_dim
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  character(len=128) :: region_name, region_name_scalar
  !--- Save the number of region written. Open the file in "w" mode if 0, else
  !    in "a"  mode
  integer, save :: num_regions_written = 0

  integer, parameter :: num_ints_to_reduce = 5
  integer, dimension(num_ints_to_reduce) :: ints_to_reduce

  ! wmax = world_max variables to have constant strides in adios file
  integer :: num_interfaces_wmax, max_nibool_interfaces_wmax, &
             num_phase_ispec_wmax, num_colors_outer_wmax, num_colors_inner_wmax

  ! user output
  if (myrank == 0) then
#if defined(USE_ADIOS)
    write(IMAIN,*) '    MPI      in ADIOS 1 file format'
#elif defined(USE_ADIOS2)
    write(IMAIN,*) '    MPI      in ADIOS 2 file format'
#endif
    call flush_IMAIN()
  endif

  ints_to_reduce(1) = num_interfaces
  ints_to_reduce(2) = max_nibool_interfaces
  ints_to_reduce(3) = num_phase_ispec
  ints_to_reduce(4) = num_colors_outer
  ints_to_reduce(5) = num_colors_inner

  call max_all_all_veci(ints_to_reduce,num_ints_to_reduce)

  num_interfaces_wmax        = ints_to_reduce(1)
  max_nibool_interfaces_wmax = ints_to_reduce(2)
  num_phase_ispec_wmax       = ints_to_reduce(3)
  num_colors_outer_wmax      = ints_to_reduce(4)
  num_colors_inner_wmax      = ints_to_reduce(5)

  ! note: the number of interfaces could be different for different MPI processes
  !       (e.g., depending if a slice is located at a chunk edge or inside)
  !       determining the maximum values helps to assign a common rule for array offsets for all processes.
  !
  ! checks
  if (num_interfaces > num_interfaces_wmax .or. &
      max_nibool_interfaces > max_nibool_interfaces_wmax .or. &
      num_phase_ispec > num_phase_ispec_wmax .or. &
      num_colors_inner > num_colors_inner_wmax .or. &
      num_colors_outer > num_colors_outer_wmax) then
    print *,'Error: rank ',myrank,' num_interfaces ', &
            num_interfaces,max_nibool_interfaces,num_phase_ispec,num_colors_inner,num_colors_outer, &
            ' for MPI arrays differ from max values: ', &
            num_interfaces_wmax,max_nibool_interfaces_wmax,num_phase_ispec_wmax,num_colors_inner_wmax,num_colors_outer_wmax
    call exit_mpi(myrank,'Error num_interfaces for MPI arrays in adios saved file')
  endif

  ! create the name for the database of the current slide and region
  !call create_name_database_adios(prname,iregion_code,LOCAL_PATH)

  write(region_name,"('reg',i1, '/')") iregion_code
  write(region_name_scalar,"('reg',i1)") iregion_code
  write(group_name,"('SPECFEM3D_GLOBE_MPI_ARRAYS_reg',i1)") iregion_code

  call init_adios_group(myadios_group,group_name)

  !--- Define ADIOS variables -----------------------------
  !! MPI interfaces
  group_size_inc = 0

  ! used only for file checking
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, "myrank", myrank)

  ! interfaces
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(num_interfaces))
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(max_nibool_interfaces))

  if (num_interfaces_wmax > 0) then
    local_dim = num_interfaces_wmax
    call define_adios_global_array1D(myadios_group,group_size_inc,local_dim,region_name,STRINGIFY_VAR(my_neighbors))
    call define_adios_global_array1D(myadios_group,group_size_inc,local_dim,region_name,STRINGIFY_VAR(nibool_interfaces))

    local_dim = max_nibool_interfaces_wmax * num_interfaces_wmax
    call define_adios_global_array1D(myadios_group,group_size_inc,local_dim,region_name,STRINGIFY_VAR(ibool_interfaces))
  endif

  ! inner/outer elements
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec_inner))
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec_outer))
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(num_phase_ispec))

  if (num_phase_ispec_wmax > 0) then
    local_dim = num_phase_ispec_wmax * 2
    call define_adios_global_array1D(myadios_group,group_size_inc,local_dim,region_name,STRINGIFY_VAR(phase_ispec_inner))
  endif

  ! mesh coloring
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(num_colors_outer))
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(num_colors_inner))

  if (USE_MESH_COLORING_GPU) then
    local_dim = num_colors_outer_wmax + num_colors_inner_wmax
    call define_adios_global_array1D(myadios_group,group_size_inc,local_dim,region_name,STRINGIFY_VAR(num_elem_colors))
  endif

  !--- Open an ADIOS handler to the restart file. ---------
  outputname = get_adios_filename(trim(LOCAL_PATH) // "/solver_data_mpi")
  ! user output
  if (myrank == 0) write(IMAIN,*) '    saving arrays in ADIOS file: ',trim(outputname)

  if (num_regions_written == 0) then
    ! opens file for writing
    call open_file_adios_write(myadios_file,myadios_group,outputname,group_name)
  else
    ! opens file for writing in append mode
    call open_file_adios_write_append(myadios_file,myadios_group,outputname,group_name)
  endif

  call set_adios_group_size(myadios_file,group_size_inc)

  !--- Schedule writes for the previously defined ADIOS variables

  ! used only for file checking
  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "myrank", myrank)

  ! MPI interfaces
  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "num_interfaces", num_interfaces)
  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "max_nibool_interfaces", max_nibool_interfaces)

  if (num_interfaces_wmax > 0) then
    ! note: using the *_wmax array sizes for local_dim is providing the same local_dim/global_dim/offset values
    !       in the adios file for all rank processes. this mimicks the same chunk sizes for all processes in ADIOS files.
    !       this helps when reading back arrays using offsets based on the local_dim value.
    !
    ! we thus use num_interfaces_wmax here rather than num_interfaces
    local_dim = num_interfaces_wmax
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // "my_neighbors", my_neighbors)
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // "nibool_interfaces", nibool_interfaces)

    local_dim = max_nibool_interfaces_wmax * num_interfaces_wmax
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // "ibool_interfaces", ibool_interfaces)
  endif

  ! inner/outer elements
  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "nspec_inner", nspec_inner)
  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "nspec_outer", nspec_outer)
  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "num_phase_ispec", num_phase_ispec)

  if (num_phase_ispec_wmax > 0) then
    local_dim = num_phase_ispec_wmax * 2
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // "phase_ispec_inner", phase_ispec_inner)
  endif

  ! mesh coloring
  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "num_colors_outer", num_colors_outer)
  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "num_colors_inner", num_colors_inner)

  if (USE_MESH_COLORING_GPU) then
    local_dim = num_colors_outer_wmax + num_colors_inner_wmax
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // "num_elem_colors", num_elem_colors)
  endif

  !--- Reset the path to zero and perform the actual write to disk
  call write_adios_perform(myadios_file)
  ! closes file
  call close_file_adios(myadios_file)

  num_regions_written = num_regions_written + 1

  end subroutine save_mpi_arrays_adios


!===============================================================================
!> \brief Write boundary conditions (MOHO, 400, 600) to a single ADIOS file.

  subroutine save_arrays_boundary_adios()

! saves arrays for boundaries such as MOHO, 400 and 670 discontinuities

  use constants

  use meshfem_par, only: &
    myrank, LOCAL_PATH

  use meshfem_models_par, only: &
    HONOR_1D_SPHERICAL_MOHO
    !SAVE_BOUNDARY_MESH,HONOR_1D_SPHERICAL_MOHO,SUPPRESS_CRUSTAL_MESH

  use regions_mesh_par2, only: &
    NSPEC2D_MOHO, NSPEC2D_400, NSPEC2D_670, &
    ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot, &
    ibelm_670_top,ibelm_670_bot,normal_moho,normal_400,normal_670, &
    ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,ispec2D_400_bot, &
    ispec2D_670_top,ispec2D_670_bot ! prname

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! local parameters
  ! local parameters
  character(len=MAX_STRING_LEN) :: outputname, group_name
  integer(kind=8) :: local_dim
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  character(len=128) :: region_name, region_name_scalar
  integer, parameter :: iregion_code = IREGION_CRUST_MANTLE
  !--- Save the number of region written. Open the file in "w" mode if 0, else
  !    in "a"  mode
  integer, save :: num_regions_written = 0

  ! user output
  if (myrank == 0) then
#if defined(USE_ADIOS)
    write(IMAIN,*) '    boundary in ADIOS 1 file format'
#elif defined(USE_ADIOS2)
    write(IMAIN,*) '    boundary in ADIOS 2 file format'
#endif
    call flush_IMAIN()
  endif

  ! first check the number of surface elements are the same for Moho, 400, 670
  if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
    if (ispec2D_moho_top /= NSPEC2D_MOHO .or. ispec2D_moho_bot /= NSPEC2D_MOHO)&
           call exit_mpi(myrank, 'Not the same number of Moho surface elements')
  endif
  if (ispec2D_400_top /= NSPEC2D_400 .or. ispec2D_400_bot /= NSPEC2D_400) &
           call exit_mpi(myrank,'Not the same number of 400 surface elements')
  if (ispec2D_670_top /= NSPEC2D_670 .or. ispec2D_670_bot /= NSPEC2D_670) &
           call exit_mpi(myrank,'Not the same number of 670 surface elements')

  write(region_name,"('reg',i1, '/')") iregion_code
  write(region_name_scalar,"('reg',i1)") iregion_code

  group_name = "SPECFEM3D_GLOBE_BOUNDARY_DISC"
  call init_adios_group(myadios_group,group_name)

  !--- Define ADIOS variables -----------------------------
  group_size_inc = 0
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec2d_moho))
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec2d_400))
  call define_adios_scalar(myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec2d_670))

  local_dim = NSPEC2D_MOHO
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_moho_top))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_moho_bot))

  local_dim = NSPEC2D_400
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_400_top))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_400_bot))

  local_dim = NSPEC2D_670
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_670_top))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_670_bot))

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_MOHO
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(normal_moho))

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_400
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(normal_400))

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_670
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(normal_670))

  !--- Open an ADIOS handler to the restart file. ---------
  outputname = get_adios_filename(trim(LOCAL_PATH) // "/boundary_disc")
  ! user output
  if (myrank == 0) write(IMAIN,*) '    saving arrays in ADIOS file: ',trim(outputname)

  if (num_regions_written == 0) then
    ! opens file for writing
    call open_file_adios_write(myadios_file,myadios_group,outputname,group_name)
  else
    ! opens file for writing in append mode
    call open_file_adios_write_append(myadios_file,myadios_group,outputname,group_name)
  endif

  call set_adios_group_size(myadios_file,group_size_inc)

  !--- Schedule writes for the previously defined ADIOS variables
  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "NSPEC2D_MOHO", NSPEC2D_MOHO)
  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "NSPEC2D_400", NSPEC2D_400)
  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // "NSPEC2D_670", NSPEC2D_670)

  local_dim = NSPEC2D_MOHO
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                             trim(region_name) // STRINGIFY_VAR(ibelm_moho_top))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                             trim(region_name) // STRINGIFY_VAR(ibelm_moho_bot))

  local_dim = NSPEC2D_400
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                              trim(region_name) // STRINGIFY_VAR(ibelm_400_top))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                              trim(region_name) // STRINGIFY_VAR(ibelm_400_bot))

  local_dim = NSPEC2D_670
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                              trim(region_name) // STRINGIFY_VAR(ibelm_670_top))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                              trim(region_name) // STRINGIFY_VAR(ibelm_670_bot))

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_MOHO
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                trim(region_name) // STRINGIFY_VAR(normal_moho))

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_400
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                 trim(region_name) // STRINGIFY_VAR(normal_400))

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_670
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                 trim(region_name) // STRINGIFY_VAR(normal_670))

  !--- Reset the path to zero and perform the actual write to disk
  call write_adios_perform(myadios_file)
  ! closes file
  call close_file_adios(myadios_file)

  num_regions_written = num_regions_written + 1

  end subroutine save_arrays_boundary_adios

