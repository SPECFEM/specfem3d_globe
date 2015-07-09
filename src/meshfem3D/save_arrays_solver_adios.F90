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

#include "config.fh"

!===============================================================================
!> \brief Main routine to save the arrays from the mesher to the solver with the
!!        help of ADIOS
!! \param myrank The MPI rank of the current process
!! \param nspec Number of GLL points per element
!! \param nglob Number of mesh points
!! \param idoubling Array of information on every mesh point
!! \param ibool Array of information on every mesh point
!! \param iregion_code The region the absorbing condition is written for. Check
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


  use adios_write_mod,only: adios_declare_group,adios_select_method,adios_open,adios_group_size

  use adios_helpers_mod,only: define_adios_global_array1D,write_adios_global_1d_array, &
    define_adios_global_real_1d_array,define_adios_scalar, &
    check_adios_err

  use constants

  use meshfem3D_models_par,only: &
    OCEANS,TRANSVERSE_ISOTROPY,HETEROGEN_3D_MANTLE,ANISOTROPIC_3D_MANTLE, &
    ANISOTROPIC_INNER_CORE,ATTENUATION

  use meshfem3D_par,only: &
    NCHUNKS,ABSORBING_CONDITIONS,SAVE_MESH_FILES, LOCAL_PATH, &
    ADIOS_FOR_SOLVER_MESHFILES, &
    ROTATION,EXACT_MASS_MATRIX_FOR_ROTATION

  use create_regions_mesh_par2,only: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
    gammaxstore,gammaystore,gammazstore, &
    rhostore,dvpstore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore,&
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
    rmassx,rmassy,rmassz,rmass_ocean_load, &
    b_rmassx,b_rmassy, &
    ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
    normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top, &
    jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax, &
    jacobian2D_bottom,jacobian2D_top, &
    rho_vp,rho_vs, &
    nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
    ispec_is_tiso,tau_s,T_c_source,tau_e_store,Qmu_store, &
    nspec_actually, nspec_ani, nspec_stacey, nglob_xy, nglob_oceans ! prname,

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
  integer,intent(in) :: NSPEC2D_TOP,NSPEC2D_BOTTOM, &
      NSPEC2DMAX_XMIN_XMAX, NSPEC2DMAX_YMIN_YMAX

  ! local parameters
  integer :: NSPEC2D_TOP_wmax,NSPEC2D_BOTTOM_wmax, &
      NSPEC2DMAX_XMIN_XMAX_wmax, NSPEC2DMAX_YMIN_YMAX_wmax
  integer, parameter :: num_ints_to_reduce = 4
  integer, dimension(num_ints_to_reduce) :: ints_to_reduce

  integer :: i,j,k,ispec,iglob,ier
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: tmp_array_x, &
      tmp_array_y, tmp_array_z

  ! local parameters
  character(len=MAX_STRING_LEN) :: reg_name, outputname, group_name
  integer :: comm, local_dim
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, handle
  integer(kind=8)         :: adios_totalsize
  character(len=128)      :: region_name, region_name_scalar

  !--- Save the number of region written. Open the file in "w" mode if 0, else
  !    in "a"  mode
  integer, save :: num_regions_written = 0

  integer :: sizeprocs

  ! number of MPI processes
  call world_size(sizeprocs)

  ! create a prefix for the file name such as LOCAL_PATH/regX_
  call create_name_database_adios(reg_name,iregion_code,LOCAL_PATH)
  write(region_name,"('reg',i1, '/')") iregion_code
  write(region_name_scalar,"('reg',i1)") iregion_code

  !---------------------------------------------------------
  !--- Solver data arrays ----------------------------------
  !---------------------------------------------------------

  ! save arrays for the solver to run.
  write(group_name,"('SPECFEM3D_GLOBE_ARRAYS_SOLVER_reg',i1)") iregion_code

  ! Alias COMM_WORLD to use ADIOS
  call world_duplicate(comm)

  ! set the adios group size to 0 before incremented by calls to
  ! helpers functions.
  group_size_inc = 0

  call adios_declare_group(adios_group, group_name, "", 1, adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  ! We set the transport method to 'MPI'. This seems to be the correct choice
  ! for now. We might want to move this to the constant.h file later on.
  call adios_select_method(adios_group, ADIOS_TRANSPORT_METHOD, "", "", adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  !--- Define ADIOS variables -----------------------------
  ! save nspec and nglob, to be used in combine_paraview_data
  call define_adios_scalar (adios_group, group_size_inc, &
                            region_name_scalar, STRINGIFY_VAR(nspec))
  call define_adios_scalar (adios_group, group_size_inc, &
                            region_name_scalar, STRINGIFY_VAR(nglob))

  local_dim = nglob
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "x_global", tmp_array_x)
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "y_global", tmp_array_y)
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "z_global", tmp_array_z)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(xstore))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ystore))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(zstore))

  local_dim = NGLLX * NGLLY * NGLLZ * nspec
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(rhostore))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(kappavstore))

  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibool))

  if (iregion_code /= IREGION_OUTER_CORE) then

    if (.not. (ANISOTROPIC_3D_MANTLE .and. &
               iregion_code == IREGION_CRUST_MANTLE)) then
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(muvstore))
    endif

    if (TRANSVERSE_ISOTROPY) then
      if (iregion_code == IREGION_CRUST_MANTLE .and. &
          .not. ANISOTROPIC_3D_MANTLE) then
        call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(kappahstore))
        call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(muhstore))
        call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(eta_anisostore))
      endif
    endif
  endif

  local_dim = nspec
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(idoubling))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ispec_is_tiso))

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_actually
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(xixstore))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(xiystore))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(xizstore))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(etaxstore))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(etaystore))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(etazstore))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(gammaxstore))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(gammaystore))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(gammazstore))

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_ani
  if (iregion_code /= IREGION_OUTER_CORE) then
    !   save anisotropy in the inner core only
    if (ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) then
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c11store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c33store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c12store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c13store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c44store))
    endif
    if (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c11store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c12store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c13store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c14store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c15store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c16store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c22store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c23store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c24store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c25store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c26store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c33store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c34store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c35store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c36store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c44store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c45store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c46store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c55store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c56store))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(c66store))
    endif
  endif

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_stacey
  if (ABSORBING_CONDITIONS) then
    if (iregion_code == IREGION_CRUST_MANTLE) then
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(rho_vp))
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(rho_vs))
    else if (iregion_code == IREGION_OUTER_CORE) then
      call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(rho_vp))
    endif
  endif

  local_dim = nglob_xy
  if ((NCHUNKS /= 6 .and. ABSORBING_CONDITIONS .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
     (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
     (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION .and. iregion_code == IREGION_INNER_CORE)) then
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(rmassx))
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(rmassy))
  endif
  local_dim = nglob
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(rmassz))

 ! mass matrices for backward simulation when ROTATION is .true.
  local_dim = nglob_xy
  if (EXACT_MASS_MATRIX_FOR_ROTATION) then
    if ((ROTATION .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
       (ROTATION .and. iregion_code == IREGION_INNER_CORE)) then
      call define_adios_global_real_1d_array(adios_group,group_size_inc, local_dim,region_name, STRINGIFY_VAR(b_rmassx) )
      call define_adios_global_real_1d_array(adios_group,group_size_inc, local_dim,region_name, STRINGIFY_VAR(b_rmassy) )
    endif
  endif

  local_dim = nglob_oceans
  if (OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) then
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(rmass_ocean_load))
  endif


  !--- Open an ADIOS handler to the restart file. ---------
  outputname = trim(LOCAL_PATH) // "/solver_data.bp"
  ! user output
  if (myrank == 0) write(IMAIN,*) '    saving arrays in ADIOS file: ',trim(outputname)

  if (num_regions_written == 0) then
    call adios_open (handle, group_name, outputname, "w", comm, adios_err)
    if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed'
  else
    call adios_open (handle, group_name, outputname, "a", comm, adios_err)
    if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed'
  endif

  call adios_group_size (handle, group_size_inc, adios_totalsize, adios_err)
  if (adios_err /= 0 ) stop 'Error calling adios_group_size() routine failed'

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

  !--- Schedule writes for the previously defined ADIOS variables
  ! save nspec and nglob, to be used in combine_paraview_data
  call adios_write(handle, trim(region_name) // "nspec", nspec, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_write(handle, trim(region_name) // "nglob", nglob, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = nglob
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // "x_global", tmp_array_x)
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // "y_global", tmp_array_y)
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // "z_global", tmp_array_z)
  local_dim = NGLLX * NGLLY * NGLLZ * nspec
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(xstore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(ystore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(zstore))

  local_dim = NGLLX * NGLLY * NGLLZ * nspec
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(rhostore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                trim(region_name) // STRINGIFY_VAR(kappavstore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(ibool))

  if (iregion_code /= IREGION_OUTER_CORE) then
    if (.not. (ANISOTROPIC_3D_MANTLE .and. &
               iregion_code == IREGION_CRUST_MANTLE)) then
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(muvstore))
    endif
    if (TRANSVERSE_ISOTROPY) then
      if (iregion_code == IREGION_CRUST_MANTLE .and. &
          .not. ANISOTROPIC_3D_MANTLE) then
        call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                trim(region_name) // STRINGIFY_VAR(kappahstore))
        call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(muhstore))
        call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                             trim(region_name) // STRINGIFY_VAR(eta_anisostore))
      endif
    endif
  endif

  local_dim = nspec
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                  trim(region_name) // STRINGIFY_VAR(idoubling))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                              trim(region_name) // STRINGIFY_VAR(ispec_is_tiso))

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_actually
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(xixstore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(xiystore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(xizstore))

  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                  trim(region_name) // STRINGIFY_VAR(etaxstore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                  trim(region_name) // STRINGIFY_VAR(etaystore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                  trim(region_name) // STRINGIFY_VAR(etazstore))

  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                trim(region_name) // STRINGIFY_VAR(gammaxstore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                trim(region_name) // STRINGIFY_VAR(gammaystore))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                trim(region_name) // STRINGIFY_VAR(gammazstore))

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_ani
  if (iregion_code /= IREGION_OUTER_CORE) then

    !   save anisotropy in the inner core only
    if (ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) then
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c11store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c33store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c12store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c13store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c44store))
    endif

    if (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c11store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c12store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c13store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c14store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c15store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c16store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c22store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c23store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c24store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c25store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c33store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c34store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c35store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c36store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c44store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c45store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c46store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c55store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c56store))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(c66store))
    endif
  endif

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_stacey
  if (ABSORBING_CONDITIONS) then
    if (iregion_code == IREGION_CRUST_MANTLE) then
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     trim(region_name) // STRINGIFY_VAR(rho_vp))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     trim(region_name) // STRINGIFY_VAR(rho_vs))

    else if (iregion_code == IREGION_OUTER_CORE) then
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     trim(region_name) // STRINGIFY_VAR(rho_vp))
    endif
  endif

  ! mass matrices
  local_dim = nglob_xy
  if ((NCHUNKS /= 6 .and. ABSORBING_CONDITIONS .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
     (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
     (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION .and. iregion_code == IREGION_INNER_CORE)) then
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     trim(region_name) // STRINGIFY_VAR(rmassx))
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     trim(region_name) // STRINGIFY_VAR(rmassy))
  endif

  local_dim = nglob
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                   trim(region_name) // STRINGIFY_VAR(rmassz))

  ! mass matrices for backward simulation when ROTATION is .true.
  local_dim = nglob_xy
  if (EXACT_MASS_MATRIX_FOR_ROTATION) then
    if ((ROTATION .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
       (ROTATION .and. iregion_code == IREGION_INNER_CORE)) then
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       trim(region_name) // STRINGIFY_VAR(b_rmassx))
      call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                       trim(region_name) // STRINGIFY_VAR(b_rmassy))
    endif
  endif

  ! additional ocean load mass matrix if oceans and if we are in the crust
  local_dim = nglob_oceans
  if (OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) then
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                           trim(region_name) // STRINGIFY_VAR(rmass_ocean_load))
    if (minval(rmass_ocean_load) <= 0._CUSTOM_REAL) &
        call exit_MPI(myrank,'negative mass matrix term for the oceans')
  endif

  !--- perform the actual write to disk
  call adios_close(handle, adios_err)

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

  call max_allreduce_i(ints_to_reduce,num_ints_to_reduce)

  NSPEC2D_TOP_wmax          = ints_to_reduce(1)
  NSPEC2D_BOTTOM_wmax       = ints_to_reduce(2)
  NSPEC2DMAX_XMIN_XMAX_wmax = ints_to_reduce(3)
  NSPEC2DMAX_YMIN_YMAX_wmax = ints_to_reduce(4)

  ! debug daniel
  !call synchronize_all()
  !print*,myrank,'nspec2d top      :',NSPEC2D_TOP,NSPEC2D_TOP_wmax
  !print*,myrank,'nspec2d bottom   :',NSPEC2D_BOTTOM,NSPEC2D_BOTTOM_wmax
  !print*,myrank,'nspec2d xmin_xmax:',NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_XMIN_XMAX_wmax
  !print*,myrank,'nspec2d ymin_ymax:',NSPEC2DMAX_YMIN_YMAX,NSPEC2DMAX_YMIN_YMAX_wmax
  !call synchronize_all()

  ! checks
  if (NSPEC2D_TOP /= NSPEC2D_TOP_wmax .or. &
      NSPEC2D_BOTTOM /= NSPEC2D_BOTTOM_wmax .or. &
      NSPEC2DMAX_XMIN_XMAX /= NSPEC2DMAX_XMIN_XMAX_wmax .or. &
      NSPEC2DMAX_YMIN_YMAX /= NSPEC2DMAX_YMIN_YMAX_wmax) then
    print*,myrank,'Error nspec2d for coupling surfaces'
    call exit_mpi(myrank,'Error nspec2d for coupling surfaces in adios saved file')
  endif

  ! save boundary arrays in ADIOS files
  write(group_name,"('SPECFEM3D_GLOBE_BOUNDARY_reg',i1)") iregion_code
  ! set the adios group size to 0 before incremented by calls to
  ! helpers functions.
  group_size_inc = 0
  call adios_declare_group(adios_group, group_name, "", 1, adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  call adios_select_method(adios_group, ADIOS_TRANSPORT_METHOD, "", "", adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  !--- Define ADIOS variables -----------------------------
  call define_adios_scalar (adios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec2D_xmin))
  call define_adios_scalar (adios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec2D_xmax))
  call define_adios_scalar (adios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec2D_ymin))
  call define_adios_scalar (adios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec2D_ymax))
  call define_adios_scalar (adios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(NSPEC2D_BOTTOM))
  call define_adios_scalar (adios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(NSPEC2D_TOP))

  ! boundary elements
  local_dim = NSPEC2DMAX_XMIN_XMAX
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_xmin))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_xmax))

  local_dim = NSPEC2DMAX_YMIN_YMAX
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_ymin))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_ymax))

  local_dim = NSPEC2D_BOTTOM
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_bottom))

  local_dim = NSPEC2D_TOP
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_top))

  ! normals
  local_dim = NDIM*NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(normal_xmin))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(normal_xmax))

  local_dim = NDIM*NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(normal_ymin))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(normal_ymax))

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_BOTTOM
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(normal_bottom))

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_TOP
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(normal_top))

  ! Jacobians
  local_dim = NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(jacobian2D_xmin))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(jacobian2D_xmax))

  local_dim = NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(jacobian2D_ymin))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(jacobian2D_ymax))

  local_dim = NGLLX*NGLLY*NSPEC2D_BOTTOM
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(jacobian2D_bottom))

  local_dim = NGLLX*NGLLY*NSPEC2D_TOP
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(jacobian2D_top))

  !--- Open an ADIOS handler to the restart file. ---------
  outputname = trim(LOCAL_PATH) // "/boundary.bp"
  ! user output
  if (myrank == 0) write(IMAIN,*) '    saving arrays in ADIOS file: ',trim(outputname)

  if (num_regions_written == 0) then
    call adios_open (handle, group_name, outputname, "w", comm, adios_err)
    if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed'
  else
    call adios_open (handle, group_name, outputname, "a", comm, adios_err)
    if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed'
  endif

  call adios_group_size (handle, group_size_inc, adios_totalsize, adios_err)
  if (adios_err /= 0 ) stop 'Error calling adios_group_size() routine failed'

  !--- Schedule writes for the previously defined ADIOS variables
  call adios_write(handle, trim(region_name) // "nspec2D_xmin", nspec2D_xmin, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_write(handle, trim(region_name) // "nspec2D_xmax", nspec2D_xmax, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_write(handle, trim(region_name) // "nspec2D_ymin", nspec2D_ymin, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_write(handle, trim(region_name) // "nspec2D_ymax", nspec2D_ymax, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_write(handle, trim(region_name) // "NSPEC2D_BOTTOM", NSPEC2D_BOTTOM, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_write(handle, trim(region_name) // "NSPEC2D_TOP", NSPEC2D_TOP, adios_err)
  call check_adios_err(myrank,adios_err)

  ! boundary elements
  local_dim = NSPEC2DMAX_XMIN_XMAX
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                 trim(region_name) // STRINGIFY_VAR(ibelm_xmin))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                 trim(region_name) // STRINGIFY_VAR(ibelm_xmax))

  local_dim = NSPEC2DMAX_YMIN_YMAX
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                 trim(region_name) // STRINGIFY_VAR(ibelm_ymin))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                 trim(region_name) // STRINGIFY_VAR(ibelm_ymax))

  local_dim = NSPEC2D_BOTTOM
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                               trim(region_name) // STRINGIFY_VAR(ibelm_bottom))

  local_dim = NSPEC2D_TOP
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                  trim(region_name) // STRINGIFY_VAR(ibelm_top))

  ! normals
  local_dim = NDIM*NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                trim(region_name) // STRINGIFY_VAR(normal_xmin))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                trim(region_name) // STRINGIFY_VAR(normal_xmax))

  local_dim = NDIM*NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                trim(region_name) // STRINGIFY_VAR(normal_ymin))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                trim(region_name) // STRINGIFY_VAR(normal_ymax))

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_BOTTOM
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                              trim(region_name) // STRINGIFY_VAR(normal_bottom))

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_TOP
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                 trim(region_name) // STRINGIFY_VAR(normal_top))

  ! Jacobians
  local_dim = NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                            trim(region_name) // STRINGIFY_VAR(jacobian2D_xmin))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                            trim(region_name) // STRINGIFY_VAR(jacobian2D_xmax))

  local_dim = NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                            trim(region_name) // STRINGIFY_VAR(jacobian2D_ymin))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                            trim(region_name) // STRINGIFY_VAR(jacobian2D_ymax))

  local_dim = NGLLX*NGLLY*NSPEC2D_BOTTOM
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                          trim(region_name) // STRINGIFY_VAR(jacobian2D_bottom))

  local_dim = NGLLX*NGLLY*NSPEC2D_TOP
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                             trim(region_name) // STRINGIFY_VAR(jacobian2D_top))

  !--- Reset the path to zero and perform the actual write to disk
  call adios_close(handle, adios_err)

  !---------------------------------------------------------
  !--- Attenuation arrays ----------------------------------
  !---------------------------------------------------------
  if (ATTENUATION) then
    write(group_name,"('SPECFEM3D_GLOBE_ATTENUATION_reg',i1)") iregion_code
    group_size_inc = 0
    call adios_declare_group(adios_group, group_name, "", 0, adios_err)
    ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
    !call check_adios_err(myrank,adios_err)

    call adios_select_method(adios_group, ADIOS_TRANSPORT_METHOD, "", "", adios_err)
    ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
    !call check_adios_err(myrank,adios_err)

    !--- Define ADIOS variables -----------------------------
    call define_adios_scalar (adios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(t_c_source))

    local_dim = size(tau_s)
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(tau_s))
    local_dim = size(tau_e_store)
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(tau_e_store))
    local_dim = size(Qmu_store)
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(Qmu_store))

    !--- Open an ADIOS handler to the restart file. ---------
    outputname = trim(LOCAL_PATH) // "/attenuation.bp"
    ! user output
    if (myrank == 0) write(IMAIN,*) '    saving arrays in ADIOS file: ',trim(outputname)

    if (num_regions_written == 0) then
      call adios_open (handle, group_name, outputname, "w", comm, adios_err)
      if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed'
    else
      call adios_open (handle, group_name, outputname, "a", comm, adios_err)
      if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed'
    endif

    call adios_group_size (handle, group_size_inc, adios_totalsize, adios_err)
    if (adios_err /= 0 ) stop 'Error calling adios_group_size() routine failed'

    !--- Schedule writes for the previously defined ADIOS variables
    call adios_write(handle, trim(region_name) // "t_c_source", T_c_source, adios_err)

    local_dim = size (tau_s)
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                     trim(region_name) // STRINGIFY_VAR(tau_s))
    local_dim = size (tau_e_store)
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                trim(region_name) // STRINGIFY_VAR(tau_e_store))
    local_dim = size (Qmu_store)
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                  trim(region_name) // STRINGIFY_VAR(Qmu_store))

    !--- Reset the path to zero and perform the actual write to disk
    call adios_close(handle, adios_err)

  endif

  !---------------------------------------------------------
  !--- dvp arrays ------------------------------------------
  !---------------------------------------------------------
  if (HETEROGEN_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
    write(group_name,"('SPECFEM3D_GLOBE_DVP_reg',i1)") iregion_code
    group_size_inc = 0
    call adios_declare_group(adios_group, group_name, "", 0, adios_err)
    ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
    !call check_adios_err(myrank,adios_err)

    call adios_select_method(adios_group, ADIOS_TRANSPORT_METHOD, "", "", adios_err)
    ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
    !call check_adios_err(myrank,adios_err)

    !--- Define ADIOS variables -----------------------------
    local_dim = size (dvpstore)
    call define_adios_global_array1D(adios_group, group_size_inc, &
                                     local_dim, region_name, &
                                     "dvp", dvpstore)
    !--- Open an ADIOS handler to the restart file. ---------
    outputname = trim(LOCAL_PATH) // "/dvp.bp"
    ! user output
    if (myrank == 0) write(IMAIN,*) '    saving arrays in ADIOS file: ',trim(outputname)

    if (num_regions_written == 0) then
      call adios_open (handle, group_name, outputname, "w", comm, adios_err)
      if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed'
    else
      call adios_open (handle, group_name, outputname, "a", comm, adios_err)
      if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed'
    endif

    call adios_group_size (handle, group_size_inc, adios_totalsize, adios_err)
    if (adios_err /= 0 ) stop 'Error calling adios_group_size() routine failed'

    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, trim(region_name) // "dvp", dvpstore)

    !--- Reset the path to zero and perform the actual write to disk
    call adios_close(handle, adios_err)
  endif

  !---------------------------------------------------------
  !--- meshfiles arrays ------------------------------------
  !---------------------------------------------------------
  ! uncomment for vp & vs model storage
  if (SAVE_MESH_FILES) then
    ! outputs model files in binary format
    if (ADIOS_FOR_SOLVER_MESHFILES) then
      call save_arrays_solver_meshfiles_adios(myrank,iregion_code,nspec)
    else
      call save_arrays_solver_meshfiles(myrank,nspec)
    endif
  endif

  num_regions_written = num_regions_written + 1

end subroutine save_arrays_solver_adios


!===============================================================================
!> \brief Save the meshfiles that will be used by the solver in an ADIOS format.
!!
!! \param myrank The MPI rank of the current process.
!! \param iregion_code Code of the region considered. See constant.h for details
!! \param reg_name Output file prefix with the name of the region included
!! \param nspec Number of GLL points per spectral elements
subroutine save_arrays_solver_meshfiles_adios(myrank, iregion_code, nspec)

  ! outputs model files in binary format
  use adios_write_mod,only: adios_declare_group,adios_select_method,adios_open,adios_group_size

  use adios_helpers_mod,only: define_adios_global_array1D,write_adios_global_1d_array,check_adios_err

  use constants

  use meshfem3D_par, only: &
    LOCAL_PATH

  use meshfem3D_models_par,only: &
    TRANSVERSE_ISOTROPY,ATTENUATION, &
    ATTENUATION_3D,ATTENUATION_1D_WITH_3D_STORAGE

  use create_regions_mesh_par2,only: &
    rhostore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
    Qmu_store

  implicit none

  integer :: myrank, nspec, iregion_code

  ! local parameters
  integer :: i,j,k,ispec
  real(kind=CUSTOM_REAL) :: scaleval1,scaleval2
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: temp_store
  real(kind=CUSTOM_REAL), dimension(1,1,1,1) :: dummy_ijke

  ! local parameters
  character(len=MAX_STRING_LEN) :: outputname, group_name
  integer :: comm, local_dim
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, handle
  integer(kind=8)         :: adios_totalsize
  character(len=128)      :: region_name, region_name_scalar
  !--- Save the number of region written. Open the file in "w" mode if 0, else
  !    in "a"  mode
  integer, save :: num_regions_written = 0

  integer :: sizeprocs

  ! number of MPI processes
  call world_size(sizeprocs)

  ! scaling factors to re-dimensionalize units
  scaleval1 = sngl( sqrt(PI*GRAV*RHOAV)*(R_EARTH/1000.0d0) )
  scaleval2 = sngl( RHOAV/1000.0d0 )

  call world_duplicate(comm)

  ! isotropic model
  write(region_name,"('reg',i1, '/')") iregion_code
  write(region_name_scalar,"('reg',i1)") iregion_code
  write(group_name,"('SPECFEM3D_GLOBE_solver_meshfiles_reg',i1)") iregion_code

  group_size_inc = 0

  call adios_declare_group(adios_group, group_name, "", 1, adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  call adios_select_method(adios_group, ADIOS_TRANSPORT_METHOD, "", "", adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  !--- Define ADIOS variables -----------------------------
  !--- vp arrays -------------------------------------------
  local_dim = size (kappavstore)
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "vp", dummy_ijke)
  !--- vs arrays -------------------------------------------
  local_dim = size (rhostore)
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "vs", dummy_ijke)
  !--- rho arrays ------------------------------------------
  local_dim = size (rhostore)
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "rho", dummy_ijke)
  ! transverse isotropic model
  if (TRANSVERSE_ISOTROPY) then
    !--- vpv arrays ----------------------------------------
    local_dim = size (kappavstore)
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "vpv", dummy_ijke)
    !--- vph arrays ----------------------------------------
    local_dim = size (kappavstore)
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "vph", dummy_ijke)
    !--- vsv arrays ----------------------------------------
    local_dim = size (rhostore)
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "vsv", dummy_ijke)
    !--- vsh arrays ----------------------------------------
    local_dim = size (rhostore)
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "vsh", dummy_ijke)
    !--- eta arrays ----------------------------------------
    local_dim = size (eta_anisostore)
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "eta", eta_anisostore)
  endif

  if (ATTENUATION) then
    !--- Qmu arrays ----------------------------------------
    local_dim = NGLLX * NGLLY * NGLLZ * nspec
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, "qmu", dummy_ijke)
  endif

  !--- Open an ADIOS handler to the restart file. ---------
  outputname = trim(LOCAL_PATH) // "/solver_meshfiles.bp"
  ! user output
  if (myrank == 0) write(IMAIN,*) '    saving arrays in ADIOS file: ',trim(outputname)

  if (num_regions_written == 0) then
    call adios_open (handle, group_name, outputname, "w", comm, adios_err)
    if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed'
  else
    call adios_open (handle, group_name, outputname, "a", comm, adios_err)
    if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed'
  endif

  call adios_group_size (handle, group_size_inc, adios_totalsize, adios_err)
  if (adios_err /= 0 ) stop 'Error calling adios_group_size() routine failed'


  !--- Schedule writes for the previously defined ADIOS variables
  ! TODO Try the new write helpers routines
  !--- vp arrays -------------------------------------------
  local_dim = size (kappavstore)
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, trim(region_name) // "vp", &
                                   sqrt( (kappavstore+4.*muvstore/3.)/rhostore )*scaleval1)

  !--- vs arrays -------------------------------------------
  local_dim = size (rhostore)
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, trim(region_name) // "vs", &
                                   sqrt( muvstore/rhostore )*scaleval1 )
  !--- rho arrays ------------------------------------------
  local_dim = size (rhostore)
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, trim(region_name) // "rho", &
                                   rhostore *scaleval2)
  ! transverse isotropic model
  if (TRANSVERSE_ISOTROPY) then
    !--- vps arrays ----------------------------------------
    local_dim = size (kappavstore)
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, trim(region_name) // "vpv", &
                                     sqrt( (kappavstore+4.*muvstore/3.)/rhostore )*scaleval1)

    local_dim = size (kappavstore)
    !--- vph arrays ----------------------------------------
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, trim(region_name) // "vph", &
                                     sqrt( (kappahstore+4.*muhstore/3.)/rhostore )*scaleval1)
    !--- vsv arrays ----------------------------------------
    local_dim = size (rhostore)
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, trim(region_name) // "vsv", &
                                     sqrt( muvstore/rhostore )*scaleval1)
    !--- vsh arrays ----------------------------------------
    local_dim = size (rhostore)
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, trim(region_name) // "vsh", &
                                     sqrt( muhstore/rhostore )*scaleval1)
    !--- eta arrays ----------------------------------------
    local_dim = size (eta_anisostore)
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, trim(region_name) // "eta", eta_anisostore)
  endif ! TRANSVERSE_ISOTROPY
  ! shear attenuation
  if (ATTENUATION) then
    !-------------------------------------------------------
    !--- Qmu arrays ----------------------------------------
    !-------------------------------------------------------
    ! saves Qmu_store to full custom_real array
    ! uses temporary array
    allocate(temp_store(NGLLX,NGLLY,NGLLZ,nspec))
    if (ATTENUATION_3D .or. ATTENUATION_1D_WITH_3D_STORAGE) then
      ! attenuation arrays are fully 3D
      temp_store(:,:,:,:) = Qmu_store(:,:,:,:)
    else
      ! attenuation array dimensions: Q_mustore(1,1,1,nspec)
      do ispec = 1,nspec
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              temp_store(i,j,k,ispec) = Qmu_store(1,1,1,ispec)
            enddo
          enddo
        enddo
      enddo
    endif

    local_dim = NGLLX * NGLLY * NGLLZ * nspec
    call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, trim(region_name) // "qmu", temp_store)

    ! frees temporary memory
    deallocate(temp_store)
  endif ! ATTENUATION

  call adios_close(handle, adios_err)

  num_regions_written = num_regions_written + 1

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
!! \param nspec_outer Number of spectral elements in the outer core
!! \param num_phase_ispec
!! \param phase_ispec_inner
!! \param num_colors_inner Number of colors for GPU computing in the inner core.
!! \param num_colors_outer Number of colors for GPU computing in the outer core.

subroutine save_mpi_arrays_adios(myrank,iregion_code,LOCAL_PATH, &
                                 num_interfaces,max_nibool_interfaces, my_neighbours,nibool_interfaces, &
                                 ibool_interfaces, nspec_inner,nspec_outer, num_phase_ispec, &
                                 phase_ispec_inner, num_colors_outer,num_colors_inner, num_elem_colors)

  use adios_write_mod,only: adios_declare_group,adios_select_method,adios_open,adios_group_size

  use adios_helpers_mod,only: define_adios_global_array1D,write_adios_global_1d_array, &
    define_adios_scalar,check_adios_err

  use constants

  implicit none

  integer :: iregion_code,myrank
  character(len=MAX_STRING_LEN) :: LOCAL_PATH
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
  character(len=MAX_STRING_LEN) :: prname, outputname, group_name
  integer :: comm, local_dim
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, handle
  integer(kind=8)         :: adios_totalsize
  character(len=128)      :: region_name, region_name_scalar
  !--- Save the number of region written. Open the file in "w" mode if 0, else
  !    in "a"  mode
  integer, save :: num_regions_written = 0

  integer, parameter :: num_ints_to_reduce = 5
  integer, dimension(num_ints_to_reduce) :: ints_to_reduce
  ! wmax = world_max variables to have constant strides in adios file
  integer :: num_interfaces_wmax, max_nibool_interfaces_wmax, &
             num_phase_ispec_wmax, num_colors_outer_wmax, num_colors_inner_wmax

  integer :: sizeprocs

  ! number of MPI processes
  call world_size(sizeprocs)

  ints_to_reduce(1) = num_interfaces
  ints_to_reduce(2) = max_nibool_interfaces
  ints_to_reduce(3) = num_phase_ispec
  ints_to_reduce(4) = num_colors_outer
  ints_to_reduce(5) = num_colors_inner

  call max_allreduce_i(ints_to_reduce,num_ints_to_reduce)

  num_interfaces_wmax        = ints_to_reduce(1)
  max_nibool_interfaces_wmax = ints_to_reduce(2)
  num_phase_ispec_wmax       = ints_to_reduce(3)
  num_colors_outer_wmax      = ints_to_reduce(4)
  num_colors_inner_wmax      = ints_to_reduce(5)

  ! create the name for the database of the current slide and region
  call create_name_database_adios(prname,iregion_code,LOCAL_PATH)

  write(region_name,"('reg',i1, '/')") iregion_code
  write(region_name_scalar,"('reg',i1)") iregion_code
  write(group_name,"('SPECFEM3D_GLOBE_MPI_ARRAYS_reg',i1)") iregion_code

  call world_duplicate(comm)

  group_size_inc = 0
  call adios_declare_group(adios_group, group_name, "", 1, adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  call adios_select_method(adios_group, ADIOS_TRANSPORT_METHOD, "", "", adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  !--- Define ADIOS variables -----------------------------
  !! MPI interfaces
  call define_adios_scalar (adios_group, group_size_inc, &
                            region_name_scalar, STRINGIFY_VAR(num_interfaces))
  if (num_interfaces > 0) then
    call define_adios_scalar (adios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(max_nibool_interfaces))
    local_dim = num_interfaces_wmax
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(my_neighbours))
    local_dim = num_interfaces_wmax
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(nibool_interfaces))
    local_dim = max_nibool_interfaces_wmax * num_interfaces_wmax
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibool_interfaces))
  endif

  ! inner/outer elements
  call define_adios_scalar (adios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec_inner))
  call define_adios_scalar (adios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec_outer))
  call define_adios_scalar (adios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(num_phase_ispec))
  if (num_phase_ispec > 0) then
    local_dim = num_phase_ispec_wmax * 2
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(phase_ispec_inner))
  endif

  ! mesh coloring
  if (USE_MESH_COLORING_GPU) then
    call define_adios_scalar (adios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(num_colors_outer))
    call define_adios_scalar (adios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(num_colors_inner))
    local_dim = num_colors_outer_wmax + num_colors_inner_wmax
    call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(num_elem_colors))
  endif

  !--- Open an ADIOS handler to the restart file. ---------
  outputname = trim(LOCAL_PATH) // "/solver_data_mpi.bp"
  ! user output
  if (myrank == 0) write(IMAIN,*) '    saving arrays in ADIOS file: ',trim(outputname)

  if (num_regions_written == 0) then
    call adios_open (handle, group_name, outputname, "w", comm, adios_err)
    if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed'
  else
    call adios_open (handle, group_name, outputname, "a", comm, adios_err)
    if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed'
  endif

  call adios_group_size (handle, group_size_inc, adios_totalsize, adios_err)
  if (adios_err /= 0 ) stop 'Error calling adios_group_size() routine failed'

  !--- Schedule writes for the previously defined ADIOS variables
  ! MPI interfaces
  call adios_write(handle, trim(region_name) // "num_interfaces", &
                   num_interfaces, adios_err)

  if (num_interfaces > 0) then
    call adios_write(handle, trim(region_name) // "max_nibool_interfaces", &
        max_nibool_interfaces, adios_err)
    !local_dim = num_interfaces_wmax
    !call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                     !local_dim, trim(region_name) // &
                                     !STRINGIFY_VAR(my_neighbours))
    !call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                     !local_dim,  trim(region_name) // &
                                     !STRINGIFY_VAR(nibool_interfaces))

    !local_dim = max_nibool_interfaces_wmax * num_interfaces_wmax
    !call write_adios_global_1d_array(handle, myrank, sizeprocs, &
                                     !local_dim, trim(region_name) // &
                                     !STRINGIFY_VAR(ibool_interfaces))

    call adios_write(handle, trim(region_name) // "my_neighbours/local_dim", num_interfaces, adios_err)
    call adios_write(handle, trim(region_name) // "my_neighbours/global_dim", num_interfaces_wmax*sizeprocs, adios_err)
    call adios_write(handle, trim(region_name) // "my_neighbours/offset", num_interfaces_wmax*myrank, adios_err)
    call adios_write(handle, trim(region_name) // "my_neighbours/array", my_neighbours, adios_err)

    call adios_write(handle, trim(region_name) // "nibool_interfaces/local_dim", num_interfaces, adios_err)
    call adios_write(handle, trim(region_name) // "nibool_interfaces/global_dim", num_interfaces_wmax*sizeprocs, adios_err)
    call adios_write(handle, trim(region_name) // "nibool_interfaces/offset", num_interfaces_wmax*myrank, adios_err)
    call adios_write(handle, trim(region_name) // "nibool_interfaces/array", nibool_interfaces, adios_err)

    call adios_write(handle, trim(region_name) // "ibool_interfaces/local_dim", &
                     max_nibool_interfaces * num_interfaces, adios_err)
    call adios_write(handle, trim(region_name) // "ibool_interfaces/global_dim", &
                     max_nibool_interfaces_wmax * num_interfaces_wmax*sizeprocs, adios_err)
    call adios_write(handle, trim(region_name) // "ibool_interfaces/offset", &
                     max_nibool_interfaces_wmax * num_interfaces_wmax*myrank, adios_err)
    call adios_write(handle, trim(region_name) // "ibool_interfaces/array", &
                     ibool_interfaces, adios_err)
  endif

  ! inner/outer elements
  call adios_write(handle, trim(region_name) // "nspec_inner", nspec_inner, adios_err)
  call adios_write(handle, trim(region_name) // "nspec_outer", nspec_outer, adios_err)
  call adios_write(handle, trim(region_name) // "num_phase_ispec", num_phase_ispec, adios_err)

  if (num_phase_ispec > 0) then
    !local_dim = num_phase_ispec_wmax * 2
    !call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                          !trim(region_name) // STRINGIFY_VAR(phase_ispec_inner))
    call adios_write(handle, trim(region_name) // "phase_ispec_inner/local_dim", 2 * num_phase_ispec, adios_err)
    call adios_write(handle, trim(region_name) // "phase_ispec_inner/global_dim", 2 * num_phase_ispec_wmax * sizeprocs, adios_err)
    call adios_write(handle, trim(region_name) // "phase_ispec_inner/offset", 2 * num_phase_ispec_wmax * myrank, adios_err)
    call adios_write(handle, trim(region_name) // "phase_ispec_inner/array", phase_ispec_inner, adios_err)
  endif

  ! mesh coloring
  if (USE_MESH_COLORING_GPU) then
    call adios_write(handle, trim(region_name) // "num_colors_outer", nspec_inner, adios_err)
    call adios_write(handle, trim(region_name) // "num_colors_inner", nspec_inner, adios_err)

    !local_dim = num_colors_outer_wmax + num_colors_inner_wmax
    !call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                            !trim(region_name) // STRINGIFY_VAR(num_elem_colors))

    call adios_write(handle, trim(region_name) // "num_elem_colors/local_dim", &
                     (num_colors_inner + num_colors_outer), adios_err)
     call adios_write(handle, trim(region_name) // "num_elem_colors/global_dim", &
                      (num_colors_inner_wmax + num_colors_outer_wmax)* sizeprocs, adios_err)
     call adios_write(handle, trim(region_name) // "num_elem_colors/offset", &
                      (num_colors_inner_wmax + num_colors_outer_wmax)* myrank, adios_err)
     call adios_write(handle, trim(region_name) // "num_elem_colors/array", &
                      num_elem_colors, adios_err)
  endif

  !--- Reset the path to zero and perform the actual write to disk
  call adios_close(handle, adios_err)

  num_regions_written = num_regions_written + 1

end subroutine save_mpi_arrays_adios


!===============================================================================
!> \brief Write boundary conditions (MOHO, 400, 600) to a single ADIOS file.

subroutine save_arrays_boundary_adios()

! saves arrays for boundaries such as MOHO, 400 and 670 discontinuities

  use constants

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
    ispec2D_670_top,ispec2D_670_bot ! prname

  use adios_helpers_mod,only: define_adios_global_array1D,write_adios_global_1d_array, &
    define_adios_scalar,check_adios_err

  implicit none

  ! local parameters
  ! local parameters
  character(len=MAX_STRING_LEN) :: outputname, group_name
  integer :: comm, local_dim
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, handle !, varid
  integer(kind=8)         :: adios_totalsize
  character(len=128)      :: region_name, region_name_scalar
  integer, parameter :: iregion_code = IREGION_CRUST_MANTLE
  !--- Save the number of region written. Open the file in "w" mode if 0, else
  !    in "a"  mode
  integer, save :: num_regions_written = 0

  integer :: sizeprocs

  ! number of MPI processes
  call world_size(sizeprocs)

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

  call world_duplicate(comm)

  group_size_inc = 0
  call adios_declare_group(adios_group, group_name, "", 1, adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  call adios_select_method(adios_group, ADIOS_TRANSPORT_METHOD, "", "", adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  !--- Define ADIOS variables -----------------------------
  call define_adios_scalar (adios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec2d_moho))
  call define_adios_scalar (adios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec2d_400))
  call define_adios_scalar (adios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(nspec2d_670))

  local_dim = NSPEC2D_MOHO
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_moho_top))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_moho_bot))
  local_dim = NSPEC2D_400
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_400_top))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_400_bot))
  local_dim = NSPEC2D_670
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_670_top))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(ibelm_670_bot))
  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_MOHO
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(normal_moho))
  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_400
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(normal_400))
  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_670
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(normal_670))

  !--- Open an ADIOS handler to the restart file. ---------
  outputname = trim(LOCAL_PATH) // "/boundary_disc.bp"
  ! user output
  if (myrank == 0) write(IMAIN,*) '    saving arrays in ADIOS file: ',trim(outputname)

  if (num_regions_written == 0) then
    call adios_open (handle, group_name, outputname, "w", comm, adios_err)
    if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed'
  else
    call adios_open (handle, group_name, outputname, "a", comm, adios_err)
    if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed'
  endif

  call adios_group_size (handle, group_size_inc, adios_totalsize, adios_err)
  if (adios_err /= 0 ) stop 'Error calling adios_group_size() routine failed'

  !--- Schedule writes for the previously defined ADIOS variables
  call adios_write(handle, trim(region_name) // "NSPEC2D_MOHO", NSPEC2D_MOHO, adios_err)
  call adios_write(handle, trim(region_name) // "NSPEC2D_400", NSPEC2D_400, adios_err)
  call adios_write(handle, trim(region_name) // "NSPEC2D_670", NSPEC2D_670, adios_err)

  local_dim = NSPEC2D_MOHO
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                             trim(region_name) // STRINGIFY_VAR(ibelm_moho_top))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                             trim(region_name) // STRINGIFY_VAR(ibelm_moho_bot))

  local_dim = NSPEC2D_400
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                              trim(region_name) // STRINGIFY_VAR(ibelm_400_top))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                              trim(region_name) // STRINGIFY_VAR(ibelm_400_bot))

  local_dim = NSPEC2D_670
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                              trim(region_name) // STRINGIFY_VAR(ibelm_670_top))
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                              trim(region_name) // STRINGIFY_VAR(ibelm_670_bot))

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_MOHO
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                trim(region_name) // STRINGIFY_VAR(normal_moho))

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_400
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                 trim(region_name) // STRINGIFY_VAR(normal_400))

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_670
  call write_adios_global_1d_array(handle, myrank, sizeprocs, local_dim, &
                                 trim(region_name) // STRINGIFY_VAR(normal_670))

  !--- Reset the path to zero and perform the actual write to disk
  call adios_close(handle, adios_err)

  num_regions_written = num_regions_written + 1

end subroutine save_arrays_boundary_adios

!!------------------------------------------------------------------------------
!!> Write local, global and offset dimensions to ADIOS
!!! \param handle Handle to the adios file
!!! \param local_dim Number of elements to be written by one process
!!! \param sizeprocs Number of MPI processes
!subroutine write_1D_global_array_adios_dims(handle, myrank, &
    !local_dim, sizeprocs)
  !use adios_write_mod

  !implicit none

  !integer(kind=8), intent(in) :: handle
  !integer, intent(in) :: sizeprocs, local_dim, myrank

  !integer :: adios_err

  !call adios_write(handle, "local_dim", local_dim, adios_err)
  !call check_adios_err(myrank,adios_err)
  !call adios_write(handle, "global_dim", local_dim*sizeprocs, adios_err)
  !call check_adios_err(myrank,adios_err)
  !call adios_write(handle, "offset", local_dim*myrank, adios_err)
  !call check_adios_err(myrank,adios_err)
!end subroutine write_1D_global_array_adios_dims

