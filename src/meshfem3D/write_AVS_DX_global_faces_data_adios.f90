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

!-------------------------------------------------------------------------------
!> \file write_AVS_DX_global_faces_data_adios.f90
!! \brief create AVS or DX 2D data for the faces of the slice,
!!        to be recombined in postprocessing
!! \author MPBL
!-------------------------------------------------------------------------------

!===============================================================================
module AVS_DX_global_faces_mod

  implicit none

  type avs_dx_global_faces_t
    integer(kind=4) :: npoin, nspecface
    real(kind=4), dimension(:), allocatable :: x_adios, y_adios, z_adios
    integer(kind=4), dimension(:), allocatable :: idoubling, iglob1, iglob2, &
        iglob3, iglob4
    real(kind=4), dimension(:), allocatable :: dvp, dvs
  endtype

contains

!===============================================================================
subroutine define_AVS_DX_global_faces_data_adios (nspec, iMPIcut_xi,iMPIcut_eta, &
                                                  ibool,mask_ibool, &
                                                  npointot, &
                                                  MODEL_3D_MANTLE_PERTUBATIONS, &
                                                  group_size_inc, avs_dx_adios)

  use constants
  use adios_helpers_mod
  use manager_adios

  implicit none

  integer :: nspec
  integer :: ibool(NGLLX,NGLLY,NGLLZ,nspec)

  logical :: MODEL_3D_MANTLE_PERTUBATIONS

  logical :: iMPIcut_xi(2,nspec)
  logical :: iMPIcut_eta(2,nspec)

  ! logical mask used to output global points only once
  integer :: npointot
  logical :: mask_ibool(npointot)

  integer(kind=8), intent(inout) :: group_size_inc
  type(avs_dx_global_faces_t), intent(inout) :: avs_dx_adios

  ! local parameters
  integer :: ispec
  integer :: iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
  integer :: npoin,nspecface
  integer :: ierr
  integer(kind=8) :: local_dim

  ! Dummy arrays for type inference inside adios helpers
  real(kind=4), dimension(:), allocatable :: dummy_real1d
  integer(kind=4), dimension(:), allocatable :: dummy_int1d

  ! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  nspecface = 0

  ! mark global AVS or DX points
  do ispec = 1,nspec
    ! only if on face
    if (iMPIcut_xi(1,ispec) .or. iMPIcut_xi(2,ispec) .or. &
        iMPIcut_eta(1,ispec) .or. iMPIcut_eta(2,ispec)) then
      iglob1 = ibool(1,1,1,ispec)
      iglob2 = ibool(NGLLX,1,1,ispec)
      iglob3 = ibool(NGLLX,NGLLY,1,ispec)
      iglob4 = ibool(1,NGLLY,1,ispec)
      iglob5 = ibool(1,1,NGLLZ,ispec)
      iglob6 = ibool(NGLLX,1,NGLLZ,ispec)
      iglob7 = ibool(NGLLX,NGLLY,NGLLZ,ispec)
      iglob8 = ibool(1,NGLLY,NGLLZ,ispec)

      ! face xi = xi_min
      if (iMPIcut_xi(1,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglob1) = .true.
        mask_ibool(iglob4) = .true.
        mask_ibool(iglob8) = .true.
        mask_ibool(iglob5) = .true.
      endif

    ! face xi = xi_max
      if (iMPIcut_xi(2,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglob2) = .true.
        mask_ibool(iglob3) = .true.
        mask_ibool(iglob7) = .true.
        mask_ibool(iglob6) = .true.
      endif

    ! face eta = eta_min
      if (iMPIcut_eta(1,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglob1) = .true.
        mask_ibool(iglob2) = .true.
        mask_ibool(iglob6) = .true.
        mask_ibool(iglob5) = .true.
      endif

    ! face eta = eta_max
      if (iMPIcut_eta(2,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglob4) = .true.
        mask_ibool(iglob3) = .true.
        mask_ibool(iglob7) = .true.
        mask_ibool(iglob8) = .true.
      endif
    endif
  enddo
  ! count global number of AVS or DX points
  npoin = count(mask_ibool(:))

  avs_dx_adios%npoin = npoin
  avs_dx_adios%nspecface = nspecface

  allocate(avs_dx_adios%x_adios(npoin), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating x_adios.")
  allocate(avs_dx_adios%y_adios(npoin), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating y_adios.")
  allocate(avs_dx_adios%z_adios(npoin), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating z_adios.")

  ! Allocate temporary arrays for AVS/DX elements.
  allocate(avs_dx_adios%idoubling(nspecface), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating idoubling.")
  allocate(avs_dx_adios%iglob1(nspecface), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating iglob1.")
  allocate(avs_dx_adios%iglob2(nspecface), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating iglob2.")
  allocate(avs_dx_adios%iglob3(nspecface), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating iglob3.")
  allocate(avs_dx_adios%iglob4(nspecface), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating iglob4.")

  allocate(dummy_real1d(npoin), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating dummy_real.")
  dummy_real1d(:) = 0.0

  allocate(dummy_int1d(nspecface), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating dummy_real.")
  dummy_int1d(:) = 0

  !--- Variables for '...AVS_DXpointsfaces.txt'
  local_dim = npoin
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                   '', "points_faces/x_value", dummy_real1d)
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                   '', "points_faces/y_value", dummy_real1d)
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                   '', "points_faces/z_value", dummy_real1d)

  !--- Variables for AVS_DXelementsfaces.txt
  local_dim = nspecface
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                  '', "elements_faces/idoubling", dummy_int1d)

  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                 '', "elements_faces/num_ibool_AVS_DX_iglob1", dummy_int1d)
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                 '', "elements_faces/num_ibool_AVS_DX_iglob2", dummy_int1d)
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                 '', "elements_faces/num_ibool_AVS_DX_iglob3", dummy_int1d)
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                 '', "elements_faces/num_ibool_AVS_DX_iglob4", dummy_int1d)

  if (MODEL_3D_MANTLE_PERTUBATIONS) then
    allocate(avs_dx_adios%dvp(nspecface), stat=ierr)
    if (ierr /= 0) call exit_MPI(myrank, "Error allocating dvp.")
    allocate(avs_dx_adios%dvs(nspecface), stat=ierr)
    if (ierr /= 0) call exit_MPI(myrank, "Error allocating dvs.")

    ! re-allocate dummy with correct size
    deallocate(dummy_real1d)
    allocate(dummy_real1d(nspecface), stat=ierr)
    if (ierr /= 0) call exit_MPI(myrank, "Error allocating dummy_real.")
    dummy_real1d(:) = 0.0

    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                     '', "elements_faces/dvp", dummy_real1d)
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, &
                                     '', "elements_faces/dvs", dummy_real1d)
  endif

  deallocate(dummy_real1d,dummy_int1d)

end subroutine define_AVS_DX_global_faces_data_adios

!===============================================================================

subroutine prepare_AVS_DX_global_faces_data_adios(nspec, &
                                                  iMPIcut_xi,iMPIcut_eta, &
                                                  ibool,idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool, &
                                                  npointot,rhostore,kappavstore,muvstore, &
                                                  nspl,rspl,ellipicity_spline,ellipicity_spline2,ELLIPTICITY, &
                                                  MODEL_3D_MANTLE_PERTUBATIONS, &
                                                  RICB,RCMB,RTOPDDOUBLEPRIME,R670,R220,R771,R400,R120,R80,RMOHO, &
                                                  RMIDDLE_CRUST,iregion_code, &
                                                  avs_dx_adios)

  use constants
  use shared_parameters, only: R_PLANET

  implicit none

  integer nspec
  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

  integer idoubling(nspec)

  logical ELLIPTICITY,MODEL_3D_MANTLE_PERTUBATIONS

  logical iMPIcut_xi(2,nspec)
  logical iMPIcut_eta(2,nspec)

  double precision RICB,RCMB,RTOPDDOUBLEPRIME,R670,R220,R771, &
    R400,R120,R80,RMOHO,RMIDDLE_CRUST

  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

  real(kind=CUSTOM_REAL) kappavstore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) muvstore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) rhostore(NGLLX,NGLLY,NGLLZ,nspec)

! logical mask used to output global points only once
  integer npointot
  logical mask_ibool(npointot)

! numbering of global AVS or DX points
  integer num_ibool_AVS_DX(npointot)

  integer ispec
  integer i,j,k,np
  integer iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
  integer npoin,numpoin,nspecface,ispecface

  double precision r,rho,vp,vs,Qkappa,Qmu
  double precision vpv,vph,vsv,vsh,eta_aniso
  double precision x,y,z
  real(kind=CUSTOM_REAL) dvp,dvs

! for ellipticity
  integer nspl
  double precision rspl(NR_DENSITY),ellipicity_spline(NR_DENSITY),ellipicity_spline2(NR_DENSITY)

  integer iregion_code

  type(avs_dx_global_faces_t), intent(inout) :: avs_dx_adios

  ! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  nspecface = 0
  dvp = 0.0
  dvs = 0.0

! mark global AVS or DX points
  do ispec = 1,nspec
  ! only if on face
    if (iMPIcut_xi(1,ispec) .or. iMPIcut_xi(2,ispec) .or. &
                iMPIcut_eta(1,ispec) .or. iMPIcut_eta(2,ispec)) then
      iglob1 = ibool(1,1,1,ispec)
      iglob2 = ibool(NGLLX,1,1,ispec)
      iglob3 = ibool(NGLLX,NGLLY,1,ispec)
      iglob4 = ibool(1,NGLLY,1,ispec)
      iglob5 = ibool(1,1,NGLLZ,ispec)
      iglob6 = ibool(NGLLX,1,NGLLZ,ispec)
      iglob7 = ibool(NGLLX,NGLLY,NGLLZ,ispec)
      iglob8 = ibool(1,NGLLY,NGLLZ,ispec)

      ! face xi = xi_min
      if (iMPIcut_xi(1,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglob1) = .true.
        mask_ibool(iglob4) = .true.
        mask_ibool(iglob8) = .true.
        mask_ibool(iglob5) = .true.
      endif

      ! face xi = xi_max
      if (iMPIcut_xi(2,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglob2) = .true.
        mask_ibool(iglob3) = .true.
        mask_ibool(iglob7) = .true.
        mask_ibool(iglob6) = .true.
      endif

      ! face eta = eta_min
      if (iMPIcut_eta(1,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglob1) = .true.
        mask_ibool(iglob2) = .true.
        mask_ibool(iglob6) = .true.
        mask_ibool(iglob5) = .true.
      endif

      ! face eta = eta_max
      if (iMPIcut_eta(2,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglob4) = .true.
        mask_ibool(iglob3) = .true.
        mask_ibool(iglob7) = .true.
        mask_ibool(iglob8) = .true.
      endif

    endif
  enddo

  ! count global number of AVS or DX points
  npoin = count(mask_ibool(:))

  ! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  ! output global AVS or DX points
  numpoin = 0
  do ispec = 1,nspec
    ! only if on face
    if (iMPIcut_xi(1,ispec) .or. iMPIcut_xi(2,ispec) .or. &
                iMPIcut_eta(1,ispec) .or. iMPIcut_eta(2,ispec)) then
      iglob1 = ibool(1,1,1,ispec)
      iglob2 = ibool(NGLLX,1,1,ispec)
      iglob3 = ibool(NGLLX,NGLLY,1,ispec)
      iglob4 = ibool(1,NGLLY,1,ispec)
      iglob5 = ibool(1,1,NGLLZ,ispec)
      iglob6 = ibool(NGLLX,1,NGLLZ,ispec)
      iglob7 = ibool(NGLLX,NGLLY,NGLLZ,ispec)
      iglob8 = ibool(1,NGLLY,NGLLZ,ispec)

      ! face xi = xi_min
      if (iMPIcut_xi(1,ispec)) then
        if (.not. mask_ibool(iglob1)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob1) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,1,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,1,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,1,1,ispec))
        endif
        if (.not. mask_ibool(iglob4)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob4) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,NGLLY,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,NGLLY,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,NGLLY,1,ispec))
        endif
        if (.not. mask_ibool(iglob8)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob8) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,NGLLY,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,NGLLY,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,NGLLY,NGLLZ,ispec))
        endif
        if (.not. mask_ibool(iglob5)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob5) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,1,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,1,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,1,NGLLZ,ispec))
        endif
        mask_ibool(iglob1) = .true.
        mask_ibool(iglob4) = .true.
        mask_ibool(iglob8) = .true.
        mask_ibool(iglob5) = .true.
      endif

      ! face xi = xi_max
      if (iMPIcut_xi(2,ispec)) then
        if (.not. mask_ibool(iglob2)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob2) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,1,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,1,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,1,1,ispec))
        endif
        if (.not. mask_ibool(iglob3)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob3) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,NGLLY,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,NGLLY,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,NGLLY,1,ispec))
        endif
        if (.not. mask_ibool(iglob7)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob7) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,NGLLY,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,NGLLY,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,NGLLY,NGLLZ,ispec))
        endif
        if (.not. mask_ibool(iglob6)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob6) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,1,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,1,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,1,NGLLZ,ispec))
        endif
        mask_ibool(iglob2) = .true.
        mask_ibool(iglob3) = .true.
        mask_ibool(iglob7) = .true.
        mask_ibool(iglob6) = .true.
      endif

      ! face eta = eta_min
      if (iMPIcut_eta(1,ispec)) then
        if (.not. mask_ibool(iglob1)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob1) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,1,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,1,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,1,1,ispec))
        endif
        if (.not. mask_ibool(iglob2)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob2) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,1,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,1,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,1,1,ispec))
        endif
        if (.not. mask_ibool(iglob6)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob6) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,1,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,1,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,1,NGLLZ,ispec))
        endif
        if (.not. mask_ibool(iglob5)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob5) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,1,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,1,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,1,NGLLZ,ispec))
        endif
        mask_ibool(iglob1) = .true.
        mask_ibool(iglob2) = .true.
        mask_ibool(iglob6) = .true.
        mask_ibool(iglob5) = .true.
      endif

      ! face eta = eta_max
      if (iMPIcut_eta(2,ispec)) then
        if (.not. mask_ibool(iglob4)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob4) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,NGLLY,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,NGLLY,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,NGLLY,1,ispec))
        endif
        if (.not. mask_ibool(iglob3)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob3) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,NGLLY,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,NGLLY,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,NGLLY,1,ispec))
        endif
        if (.not. mask_ibool(iglob7)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob7) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,NGLLY,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,NGLLY,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,NGLLY,NGLLZ,ispec))
        endif
        if (.not. mask_ibool(iglob8)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob8) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,NGLLY,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,NGLLY,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,NGLLY,NGLLZ,ispec))
        endif
        mask_ibool(iglob4) = .true.
        mask_ibool(iglob3) = .true.
        mask_ibool(iglob7) = .true.
        mask_ibool(iglob8) = .true.
      endif

    endif
  enddo

  ! check that number of global points output is okay
  if (numpoin /= npoin) &
    call exit_MPI(myrank, &
        'incorrect number of global points in AVS or DX file creation')

  ! output global AVS or DX elements

  ispecface = 0
  do ispec = 1,nspec
!    print *, ispecface, nspecface
    ! only if on face
    if (iMPIcut_xi(1,ispec) .or. iMPIcut_xi(2,ispec) .or. &
                iMPIcut_eta(1,ispec) .or. iMPIcut_eta(2,ispec)) then
      iglob1 = ibool(1,1,1,ispec)
      iglob2 = ibool(NGLLX,1,1,ispec)
      iglob3 = ibool(NGLLX,NGLLY,1,ispec)
      iglob4 = ibool(1,NGLLY,1,ispec)
      iglob5 = ibool(1,1,NGLLZ,ispec)
      iglob6 = ibool(NGLLX,1,NGLLZ,ispec)
      iglob7 = ibool(NGLLX,NGLLY,NGLLZ,ispec)
      iglob8 = ibool(1,NGLLY,NGLLZ,ispec)

      ! include lateral variations if needed
      if (MODEL_3D_MANTLE_PERTUBATIONS) then
        ! pick a point within the element and get its radius
        r = dsqrt(xstore(2,2,2,ispec)**2 &
            + ystore(2,2,2,ispec)**2 &
            + zstore(2,2,2,ispec)**2)

        if (r > RCMB/R_PLANET .and. r < R_UNIT_SPHERE) then
          ! average over the element
          dvp = 0.0
          dvs = 0.0
          np  = 0
          do k = 2,NGLLZ-1
            do j = 2,NGLLY-1
              do i = 2,NGLLX-1
                np = np+1
                x = xstore(i,j,k,ispec)
                y = ystore(i,j,k,ispec)
                z = zstore(i,j,k,ispec)

                ! take out ellipticity
                if (ELLIPTICITY) then
                  ! removes ellipticity stretch from position x/y/z
                  call revert_ellipticity(x,y,z,nspl,rspl,ellipicity_spline,ellipicity_spline2)
                endif

                ! updates radius
                r = dsqrt(x*x+y*y+z*z)

                ! gets reference model values:rho,vpv,vph,vsv,vsh and eta_aniso
                call meshfem3D_models_get1D_val(iregion_code, &
                                                idoubling(ispec), r, rho, vpv, vph, vsv, vsh, eta_aniso, &
                                                Qkappa, Qmu, RICB, RCMB, RTOPDDOUBLEPRIME, R80, R120, &
                                                R220, R400, R670, R771, RMOHO, RMIDDLE_CRUST)

                 ! calculates isotropic values
                vp = sqrt(((8.d0+4.d0*eta_aniso)*vph*vph + 3.d0*vpv*vpv &
                        + (8.d0 - 8.d0*eta_aniso)*vsv*vsv)/15.d0)
                vs = sqrt(((1.d0-2.d0*eta_aniso)*vph*vph + vpv*vpv &
                        + 5.d0*vsh*vsh + (6.d0+4.d0*eta_aniso)*vsv*vsv)/15.d0)

                if (abs(rhostore(i,j,k,ispec)) < 1.e-20) then
                  print *,'attention: rhostore close to zero', &
                      rhostore(i,j,k,ispec),r,i,j,k,ispec
                  dvp = 0.0
                  dvs = 0.0
                else if (abs(sngl(vp)) < 1.e-20) then
                  print *,' attention: vp close to zero',sngl(vp),r,i,j,k,ispec
                  dvp = 0.0
                else if (abs(sngl(vs)) < 1.e-20) then
                  print *,' attention: vs close to zero',sngl(vs),r,i,j,k,ispec
                  dvs = 0.0
                else
                  dvp = dvp + (sqrt((kappavstore(i,j,k,ispec) &
                      +4.*muvstore(i,j,k,ispec)/3.)/rhostore(i,j,k,ispec)) &
                      - sngl(vp))/sngl(vp)
                  dvs = dvs &
                      + (sqrt(muvstore(i,j,k,ispec)/rhostore(i,j,k,ispec)) &
                      - sngl(vs))/sngl(vs)
                endif

              enddo
            enddo
          enddo
          dvp = dvp / np
          dvs = dvs / np
        else
          dvp = 0.0
          dvs = 0.0
        endif
      endif

      ! face xi = xi_min
      if (iMPIcut_xi(1,ispec)) then
        ispecface = ispecface + 1
        avs_dx_adios%idoubling(ispecface) = idoubling(ispec)
        avs_dx_adios%iglob1(ispecface) = num_ibool_AVS_DX(iglob1)
        avs_dx_adios%iglob2(ispecface) = num_ibool_AVS_DX(iglob4)
        avs_dx_adios%iglob3(ispecface) = num_ibool_AVS_DX(iglob8)
        avs_dx_adios%iglob4(ispecface) = num_ibool_AVS_DX(iglob5)
        if (MODEL_3D_MANTLE_PERTUBATIONS) then
          avs_dx_adios%dvp(ispecface) = dvp
          avs_dx_adios%dvs(ispecface) = dvs
        endif
      endif

      ! face xi = xi_max
      if (iMPIcut_xi(2,ispec)) then
        ispecface = ispecface + 1
        avs_dx_adios%idoubling(ispecface) = idoubling(ispec)
        avs_dx_adios%iglob1(ispecface)= num_ibool_AVS_DX(iglob2)
        avs_dx_adios%iglob2(ispecface) = num_ibool_AVS_DX(iglob3)
        avs_dx_adios%iglob3(ispecface) = num_ibool_AVS_DX(iglob7)
        avs_dx_adios%iglob4(ispecface) = num_ibool_AVS_DX(iglob6)
        if (MODEL_3D_MANTLE_PERTUBATIONS) then
          avs_dx_adios%dvp(ispecface) = dvp
          avs_dx_adios%dvs(ispecface) = dvs
        endif
      endif

      ! face eta = eta_min
      if (iMPIcut_eta(1,ispec)) then
        ispecface = ispecface + 1
        avs_dx_adios%idoubling(ispecface) = idoubling(ispec)
        avs_dx_adios%iglob1(ispecface) = num_ibool_AVS_DX(iglob1)
        avs_dx_adios%iglob2(ispecface) = num_ibool_AVS_DX(iglob2)
        avs_dx_adios%iglob3(ispecface) = num_ibool_AVS_DX(iglob6)
        avs_dx_adios%iglob4(ispecface) = num_ibool_AVS_DX(iglob5)
        if (MODEL_3D_MANTLE_PERTUBATIONS) then
          avs_dx_adios%dvp(ispecface) = dvp
          avs_dx_adios%dvs(ispecface) = dvs
        endif
      endif

      ! face eta = eta_max
      if (iMPIcut_eta(2,ispec)) then
        ispecface = ispecface + 1
        avs_dx_adios%idoubling(ispecface) = idoubling(ispec)
        avs_dx_adios%iglob1(ispecface) = num_ibool_AVS_DX(iglob4)
        avs_dx_adios%iglob2(ispecface) = num_ibool_AVS_DX(iglob3)
        avs_dx_adios%iglob3(ispecface) = num_ibool_AVS_DX(iglob7)
        avs_dx_adios%iglob4(ispecface) = num_ibool_AVS_DX(iglob8)
        if (MODEL_3D_MANTLE_PERTUBATIONS) then
          avs_dx_adios%dvp(ispecface) = dvp
          avs_dx_adios%dvs(ispecface) = dvs
        endif
      endif

    endif
  enddo

  ! check that number of surface elements output is okay
  if (ispecface /= nspecface) &
    call exit_MPI(myrank, &
        'incorrect number of surface elements in AVS or DX file creation')

end subroutine prepare_AVS_DX_global_faces_data_adios

!===============================================================================

subroutine write_AVS_DX_global_faces_data_adios(myrank, sizeprocs, avs_dx_adios, MODEL_3D_MANTLE_PERTUBATIONS)

  use adios_helpers_mod
  use manager_adios

  implicit none
  !--- Arguments
  integer, intent(in) :: myrank, sizeprocs
  type(avs_dx_global_faces_t), intent(inout) :: avs_dx_adios ! out for adios_write
  logical MODEL_3D_MANTLE_PERTUBATIONS
  !--- Variables
  integer :: npoin, nspec
  integer(kind=8) :: local_dim

  npoin = avs_dx_adios%npoin
  nspec = avs_dx_adios%nspecface

  local_dim = npoin
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                   "points_faces/x_value", avs_dx_adios%x_adios)
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                   "points_faces/y_value", avs_dx_adios%y_adios)
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                   "points_faces/z_value", avs_dx_adios%z_adios)

  local_dim = nspec
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                   "elements_faces/idoubling", avs_dx_adios%idoubling)

  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                   "elements_faces/num_ibool_AVS_DX_iglob1", avs_dx_adios%iglob1)
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                   "elements_faces/num_ibool_AVS_DX_iglob2", avs_dx_adios%iglob2)
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                   "elements_faces/num_ibool_AVS_DX_iglob3", avs_dx_adios%iglob3)
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                   "elements_faces/num_ibool_AVS_DX_iglob4", avs_dx_adios%iglob4)


  if (MODEL_3D_MANTLE_PERTUBATIONS) then
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                     "elements_faces/dvp", avs_dx_adios%dvp)
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs, local_dim, &
                                     "elements_faces/dvs", avs_dx_adios%dvs)
  endif

end subroutine write_AVS_DX_global_faces_data_adios

!===============================================================================

subroutine free_AVS_DX_global_faces_data_adios(avs_dx_adios, MODEL_3D_MANTLE_PERTUBATIONS)

  implicit none
  !--- Arguments
  type(avs_dx_global_faces_t), intent(inout) :: avs_dx_adios
  logical MODEL_3D_MANTLE_PERTUBATIONS

  deallocate(avs_dx_adios%x_adios)
  deallocate(avs_dx_adios%y_adios)
  deallocate(avs_dx_adios%z_adios)

  deallocate(avs_dx_adios%idoubling)
  deallocate(avs_dx_adios%iglob1)
  deallocate(avs_dx_adios%iglob2)
  deallocate(avs_dx_adios%iglob3)
  deallocate(avs_dx_adios%iglob4)

  if (MODEL_3D_MANTLE_PERTUBATIONS) then
    deallocate(avs_dx_adios%dvp)
    deallocate(avs_dx_adios%dvs)
  endif

  avs_dx_adios%npoin = 0
  avs_dx_adios%nspecface = 0

end subroutine free_AVS_DX_global_faces_data_adios

end module
