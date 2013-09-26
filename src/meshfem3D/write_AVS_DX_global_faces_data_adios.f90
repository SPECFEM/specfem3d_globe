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
    real, dimension(:), allocatable :: dvp, dvs
  endtype

contains

!===============================================================================
subroutine define_AVS_DX_global_faces_data_adios (adios_group, &
    myrank, prname, nspec, iMPIcut_xi,iMPIcut_eta, &
    ibool,idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool, &
    npointot,rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
    ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
    RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
    RMIDDLE_CRUST,ROCEAN,iregion_code, &
    group_size_inc, avs_dx_adios)

  use constants
  use adios_write_mod

  implicit none

  integer(kind=8), intent(in) :: adios_group
  integer(kind=8), intent(inout) :: group_size_inc

  integer nspec,myrank
  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

  integer idoubling(nspec)

  logical ELLIPTICITY,ISOTROPIC_3D_MANTLE

  logical iMPIcut_xi(2,nspec)
  logical iMPIcut_eta(2,nspec)

  double precision RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771, &
      R400,R120,R80,RMOHO,RMIDDLE_CRUST,ROCEAN

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
  double precision x,y,z,theta,phi_dummy,cost,p20,ell,factor
  real(kind=CUSTOM_REAL) dvp,dvs

! for ellipticity
  integer nspl
  double precision rspl(NR),espl(NR),espl2(NR)

! processor identification
  character(len=150) prname

  integer iregion_code

  type(avs_dx_global_faces_t), intent(inout) :: avs_dx_adios

  integer :: ierr

  ! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  nspecface = 0

  ! mark global AVS or DX points
  do ispec=1,nspec
    ! only if on face
    if(iMPIcut_xi(1,ispec) .or. iMPIcut_xi(2,ispec) .or. &
        iMPIcut_eta(1,ispec) .or. iMPIcut_eta(2,ispec)) then
      iglob1=ibool(1,1,1,ispec)
      iglob2=ibool(NGLLX,1,1,ispec)
      iglob3=ibool(NGLLX,NGLLY,1,ispec)
      iglob4=ibool(1,NGLLY,1,ispec)
      iglob5=ibool(1,1,NGLLZ,ispec)
      iglob6=ibool(NGLLX,1,NGLLZ,ispec)
      iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
      iglob8=ibool(1,NGLLY,NGLLZ,ispec)

      ! face xi = xi_min
      if(iMPIcut_xi(1,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglob1) = .true.
        mask_ibool(iglob4) = .true.
        mask_ibool(iglob8) = .true.
        mask_ibool(iglob5) = .true.
      endif

    ! face xi = xi_max
      if(iMPIcut_xi(2,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglob2) = .true.
        mask_ibool(iglob3) = .true.
        mask_ibool(iglob7) = .true.
        mask_ibool(iglob6) = .true.
      endif

    ! face eta = eta_min
      if(iMPIcut_eta(1,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglob1) = .true.
        mask_ibool(iglob2) = .true.
        mask_ibool(iglob6) = .true.
        mask_ibool(iglob5) = .true.
      endif

    ! face eta = eta_max
      if(iMPIcut_eta(2,ispec)) then
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

  !--- Variables for '...AVS_DXpointsfaces.txt'
  call define_adios_global_real_1d_array(adios_group, "points_faces/x_value", &
      npoin, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "points_faces/y_value", &
      npoin, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "points_faces/z_value", &
      npoin, group_size_inc)
  !--- Variables for AVS_DXelementsfaces.txt
  call define_adios_global_real_1d_array(adios_group, &
      "elements_faces/idoubling", nspecface, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "elements_faces/num_ibool_AVS_DX_iglob1", nspecface, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "elements_faces/num_ibool_AVS_DX_iglob2", nspecface, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "elements_faces/num_ibool_AVS_DX_iglob3", nspecface, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "elements_faces/num_ibool_AVS_DX_iglob4", nspecface, group_size_inc)

  if(ISOTROPIC_3D_MANTLE) then
    allocate(avs_dx_adios%dvp(nspecface), stat=ierr)
    if (ierr /= 0) call exit_MPI(myrank, "Error allocating dvp.")
    allocate(avs_dx_adios%dvs(nspecface), stat=ierr)
    if (ierr /= 0) call exit_MPI(myrank, "Error allocating dvs.")
    call define_adios_global_real_1d_array(adios_group, &
        "elements_faces/dvp", dvp, group_size_inc)
    call define_adios_global_real_1d_array(adios_group, &
        "elements_faces/dvp", dvs, group_size_inc)
  endif

end subroutine define_AVS_DX_global_faces_data_adios

!===============================================================================
subroutine prepare_AVS_DX_global_faces_data_adios (myrank, prname, nspec, &
        iMPIcut_xi,iMPIcut_eta, &
        ibool,idoubling,xstore,ystore,zstore,num_ibool_AVS_DX,mask_ibool, &
        npointot,rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
        ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
        RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
        RMIDDLE_CRUST,ROCEAN,iregion_code, &
        avs_dx_adios)

  use constants

  implicit none

  integer nspec,myrank
  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

  integer idoubling(nspec)

  logical ELLIPTICITY,ISOTROPIC_3D_MANTLE

  logical iMPIcut_xi(2,nspec)
  logical iMPIcut_eta(2,nspec)

  double precision RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771, &
    R400,R120,R80,RMOHO,RMIDDLE_CRUST,ROCEAN

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
  double precision x,y,z,theta,phi_dummy,cost,p20,ell,factor
  real(kind=CUSTOM_REAL) dvp,dvs

! for ellipticity
  integer nspl
  double precision rspl(NR),espl(NR),espl2(NR)

! processor identification
  character(len=150) prname

  integer iregion_code

  type(avs_dx_global_faces_t), intent(inout) :: avs_dx_adios

  ! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  nspecface = 0

! mark global AVS or DX points
  do ispec=1,nspec
  ! only if on face
    if(iMPIcut_xi(1,ispec) .or. iMPIcut_xi(2,ispec) .or. &
                iMPIcut_eta(1,ispec) .or. iMPIcut_eta(2,ispec)) then
      iglob1=ibool(1,1,1,ispec)
      iglob2=ibool(NGLLX,1,1,ispec)
      iglob3=ibool(NGLLX,NGLLY,1,ispec)
      iglob4=ibool(1,NGLLY,1,ispec)
      iglob5=ibool(1,1,NGLLZ,ispec)
      iglob6=ibool(NGLLX,1,NGLLZ,ispec)
      iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
      iglob8=ibool(1,NGLLY,NGLLZ,ispec)

      ! face xi = xi_min
      if(iMPIcut_xi(1,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglob1) = .true.
        mask_ibool(iglob4) = .true.
        mask_ibool(iglob8) = .true.
        mask_ibool(iglob5) = .true.
      endif

      ! face xi = xi_max
      if(iMPIcut_xi(2,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglob2) = .true.
        mask_ibool(iglob3) = .true.
        mask_ibool(iglob7) = .true.
        mask_ibool(iglob6) = .true.
      endif

      ! face eta = eta_min
      if(iMPIcut_eta(1,ispec)) then
        nspecface = nspecface + 1
        mask_ibool(iglob1) = .true.
        mask_ibool(iglob2) = .true.
        mask_ibool(iglob6) = .true.
        mask_ibool(iglob5) = .true.
      endif

      ! face eta = eta_max
      if(iMPIcut_eta(2,ispec)) then
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
  do ispec=1,nspec
    ! only if on face
    if(iMPIcut_xi(1,ispec) .or. iMPIcut_xi(2,ispec) .or. &
                iMPIcut_eta(1,ispec) .or. iMPIcut_eta(2,ispec)) then
      iglob1=ibool(1,1,1,ispec)
      iglob2=ibool(NGLLX,1,1,ispec)
      iglob3=ibool(NGLLX,NGLLY,1,ispec)
      iglob4=ibool(1,NGLLY,1,ispec)
      iglob5=ibool(1,1,NGLLZ,ispec)
      iglob6=ibool(NGLLX,1,NGLLZ,ispec)
      iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
      iglob8=ibool(1,NGLLY,NGLLZ,ispec)

      ! face xi = xi_min
      if(iMPIcut_xi(1,ispec)) then
        if(.not. mask_ibool(iglob1)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob1) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,1,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,1,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,1,1,ispec))
        endif
        if(.not. mask_ibool(iglob4)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob4) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,NGLLY,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,NGLLY,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,NGLLY,1,ispec))
        endif
        if(.not. mask_ibool(iglob8)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob8) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,NGLLY,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,NGLLY,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,NGLLY,NGLLZ,ispec))
        endif
        if(.not. mask_ibool(iglob5)) then
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
      if(iMPIcut_xi(2,ispec)) then
        if(.not. mask_ibool(iglob2)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob2) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,1,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,1,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,1,1,ispec))
        endif
        if(.not. mask_ibool(iglob3)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob3) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,NGLLY,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,NGLLY,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,NGLLY,1,ispec))
        endif
        if(.not. mask_ibool(iglob7)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob7) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,NGLLY,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,NGLLY,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,NGLLY,NGLLZ,ispec))
        endif
        if(.not. mask_ibool(iglob6)) then
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
      if(iMPIcut_eta(1,ispec)) then
        if(.not. mask_ibool(iglob1)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob1) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,1,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,1,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,1,1,ispec))
        endif
        if(.not. mask_ibool(iglob2)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob2) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,1,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,1,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,1,1,ispec))
        endif
        if(.not. mask_ibool(iglob6)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob6) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,1,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,1,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,1,NGLLZ,ispec))
        endif
        if(.not. mask_ibool(iglob5)) then
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
      if(iMPIcut_eta(2,ispec)) then
        if(.not. mask_ibool(iglob4)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob4) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,NGLLY,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,NGLLY,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,NGLLY,1,ispec))
        endif
        if(.not. mask_ibool(iglob3)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob3) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,NGLLY,1,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,NGLLY,1,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,NGLLY,1,ispec))
        endif
        if(.not. mask_ibool(iglob7)) then
          numpoin = numpoin + 1
          num_ibool_AVS_DX(iglob7) = numpoin
          avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,NGLLY,NGLLZ,ispec))
          avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,NGLLY,NGLLZ,ispec))
          avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,NGLLY,NGLLZ,ispec))
        endif
        if(.not. mask_ibool(iglob8)) then
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
  if(numpoin /= npoin) &
    call exit_MPI(myrank, &
        'incorrect number of global points in AVS or DX file creation')

  ! output global AVS or DX elements

  ispecface = 0
  do ispec=1,nspec
!    print *, ispecface, nspecface
    ! only if on face
    if(iMPIcut_xi(1,ispec) .or. iMPIcut_xi(2,ispec) .or. &
                iMPIcut_eta(1,ispec) .or. iMPIcut_eta(2,ispec)) then
      iglob1=ibool(1,1,1,ispec)
      iglob2=ibool(NGLLX,1,1,ispec)
      iglob3=ibool(NGLLX,NGLLY,1,ispec)
      iglob4=ibool(1,NGLLY,1,ispec)
      iglob5=ibool(1,1,NGLLZ,ispec)
      iglob6=ibool(NGLLX,1,NGLLZ,ispec)
      iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
      iglob8=ibool(1,NGLLY,NGLLZ,ispec)
if (iglob1 > npointot) print *, myrank, "problem with iglob1", iglob1, npointot
if (iglob2 > npointot) print *, myrank, "problem with iglob2", iglob2, npointot
if (iglob3 > npointot) print *, myrank, "problem with iglob3", iglob3, npointot
if (iglob4 > npointot) print *, myrank, "problem with iglob4", iglob4, npointot
if (iglob5 > npointot) print *, myrank, "problem with iglob5", iglob5, npointot
if (iglob6 > npointot) print *, myrank, "problem with iglob6", iglob6, npointot
if (iglob7 > npointot) print *, myrank, "problem with iglob7", iglob7, npointot
if (iglob8 > npointot) print *, myrank, "problem with iglob8", iglob8, npointot

if (iglob1 < 0) print *, myrank, "problem with iglob1", iglob1, npointot
if (iglob2 < 0) print *, myrank, "problem with iglob2", iglob2, npointot
if (iglob3 < 0) print *, myrank, "problem with iglob3", iglob3, npointot
if (iglob4 < 0) print *, myrank, "problem with iglob4", iglob4, npointot
if (iglob5 < 0) print *, myrank, "problem with iglob5", iglob5, npointot
if (iglob6 < 0) print *, myrank, "problem with iglob6", iglob6, npointot
if (iglob7 < 0) print *, myrank, "problem with iglob7", iglob7, npointot
if (iglob8 < 0) print *, myrank, "problem with iglob8", iglob8, npointot

      ! include lateral variations if needed
      if(ISOTROPIC_3D_MANTLE) then
        ! pick a point within the element and get its radius
        r = dsqrt(xstore(2,2,2,ispec)**2 &
            + ystore(2,2,2,ispec)**2 &
            + zstore(2,2,2,ispec)**2)

        if(r > RCMB/R_EARTH .and. r < R_UNIT_SPHERE) then
          ! average over the element
          dvp = 0.0
          dvs = 0.0
          np =0
          do k=2,NGLLZ-1
            do j=2,NGLLY-1
              do i=2,NGLLX-1
                np=np+1
                x=xstore(i,j,k,ispec)
                y=ystore(i,j,k,ispec)
                z=zstore(i,j,k,ispec)
                r=dsqrt(x*x+y*y+z*z)
                ! take out ellipticity
                if(ELLIPTICITY) then
                  call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi_dummy)
                  cost=dcos(theta)
                  p20=0.5d0*(3.0d0*cost*cost-1.0d0)
                  call spline_evaluation(rspl,espl,espl2,nspl,r,ell)
                  factor=ONE-(TWO/3.0d0)*ell*p20
                  r=r/factor
                endif

                ! gets reference model values:rho,vpv,vph,vsv,vsh and eta_aniso
                call meshfem3D_models_get1D_val(myrank, iregion_code, &
                    idoubling(ispec), r, rho, vpv, vph, vsv, vsh, eta_aniso, &
                    Qkappa, Qmu, RICB, RCMB, RTOPDDOUBLEPRIME, R80, R120, &
                    R220, R400, R600, R670, R771, RMOHO, RMIDDLE_CRUST, ROCEAN)

                 ! calculates isotropic values
                vp = sqrt(((8.d0+4.d0*eta_aniso)*vph*vph + 3.d0*vpv*vpv &
                        + (8.d0 - 8.d0*eta_aniso)*vsv*vsv)/15.d0)
                vs = sqrt(((1.d0-2.d0*eta_aniso)*vph*vph + vpv*vpv &
                        + 5.d0*vsh*vsh + (6.d0+4.d0*eta_aniso)*vsv*vsv)/15.d0)

                if( abs(rhostore(i,j,k,ispec))< 1.e-20 ) then
                  print*,'attention: rhostore close to zero', &
                      rhostore(i,j,k,ispec),r,i,j,k,ispec
                  dvp = 0.0
                  dvs = 0.0
                else if( abs(sngl(vp))< 1.e-20 ) then
                  print*,' attention: vp close to zero',sngl(vp),r,i,j,k,ispec
                  dvp = 0.0
                else if( abs(sngl(vs))< 1.e-20 ) then
                  print*,' attention: vs close to zero',sngl(vs),r,i,j,k,ispec
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
      if(iMPIcut_xi(1,ispec)) then
        ispecface = ispecface + 1
        avs_dx_adios%idoubling(ispecface) = idoubling(ispec)
        avs_dx_adios%iglob1(ispecface) = num_ibool_AVS_DX(iglob1)
        avs_dx_adios%iglob2(ispecface) = num_ibool_AVS_DX(iglob4)
        avs_dx_adios%iglob3(ispecface) = num_ibool_AVS_DX(iglob8)
        avs_dx_adios%iglob4(ispecface) = num_ibool_AVS_DX(iglob5)
        if(ISOTROPIC_3D_MANTLE)  then
          avs_dx_adios%dvp(ispecface) = dvp
          avs_dx_adios%dvs(ispecface) = dvs
        endif
      endif

      ! face xi = xi_max
      if(iMPIcut_xi(2,ispec)) then
        ispecface = ispecface + 1
        avs_dx_adios%idoubling(ispecface) = idoubling(ispec)
        avs_dx_adios%iglob1(ispecface)= num_ibool_AVS_DX(iglob2)
        avs_dx_adios%iglob2(ispecface) = num_ibool_AVS_DX(iglob3)
        avs_dx_adios%iglob3(ispecface) = num_ibool_AVS_DX(iglob7)
        avs_dx_adios%iglob4(ispecface) = num_ibool_AVS_DX(iglob6)
        if(ISOTROPIC_3D_MANTLE) then
          avs_dx_adios%dvp(ispecface) = dvp
          avs_dx_adios%dvs(ispecface) = dvs
        endif
      endif

      ! face eta = eta_min
      if(iMPIcut_eta(1,ispec)) then
        ispecface = ispecface + 1
        avs_dx_adios%idoubling(ispecface) = idoubling(ispec)
        avs_dx_adios%iglob1(ispecface) = num_ibool_AVS_DX(iglob1)
        avs_dx_adios%iglob2(ispecface) = num_ibool_AVS_DX(iglob2)
        avs_dx_adios%iglob3(ispecface) = num_ibool_AVS_DX(iglob6)
        avs_dx_adios%iglob4(ispecface) = num_ibool_AVS_DX(iglob5)
        if(ISOTROPIC_3D_MANTLE) then
          avs_dx_adios%dvp(ispecface) = dvp
          avs_dx_adios%dvs(ispecface) = dvs
        endif
      endif

      ! face eta = eta_max
      if(iMPIcut_eta(2,ispec)) then
        ispecface = ispecface + 1
        avs_dx_adios%idoubling(ispecface) = idoubling(ispec)
        avs_dx_adios%iglob1(ispecface) = num_ibool_AVS_DX(iglob4)
        avs_dx_adios%iglob2(ispecface) = num_ibool_AVS_DX(iglob3)
        avs_dx_adios%iglob3(ispecface) = num_ibool_AVS_DX(iglob7)
        avs_dx_adios%iglob4(ispecface) = num_ibool_AVS_DX(iglob8)
        if(ISOTROPIC_3D_MANTLE) then
          avs_dx_adios%dvp(ispecface) = dvp
          avs_dx_adios%dvs(ispecface) = dvs
        endif
      endif

    endif
  enddo

  ! check that number of surface elements output is okay
  if(ispecface /= nspecface) &
    call exit_MPI(myrank,&
        'incorrect number of surface elements in AVS or DX file creation')

end subroutine prepare_AVS_DX_global_faces_data_adios

!===============================================================================
subroutine write_AVS_DX_global_faces_data_adios(adios_handle, myrank, &
    sizeprocs, avs_dx_adios, ISOTROPIC_3D_MANTLE)

  use adios_write_mod
  implicit none
  !--- Arguments
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  type(avs_dx_global_faces_t), intent(inout) :: avs_dx_adios ! out for adios_write
  logical ISOTROPIC_3D_MANTLE
  !--- Variables
  integer :: npoin, nspec
  integer :: ierr

  npoin = avs_dx_adios%npoin
  nspec = avs_dx_adios%nspecface

  call adios_set_path(adios_handle, "points_faces/x_value", ierr)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      npoin, sizeprocs)
  call adios_write(adios_handle, "array", avs_dx_adios%x_adios, ierr)

  call adios_set_path(adios_handle, "points_faces/y_value", ierr)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      npoin, sizeprocs)
  call adios_write(adios_handle, "array", avs_dx_adios%y_adios, ierr)

  call adios_set_path(adios_handle, "points_faces/z_value", ierr)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      npoin, sizeprocs)
  call adios_write(adios_handle, "array", avs_dx_adios%z_adios, ierr)


  call adios_set_path(adios_handle, "elements_faces/idoubling", ierr)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      nspec, sizeprocs)
  call adios_write(adios_handle, "array", avs_dx_adios%idoubling, ierr)


  call adios_set_path(adios_handle, &
      "elements_faces/num_ibool_AVS_DX_iglob1", ierr)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      nspec, sizeprocs)
  call adios_write(adios_handle, "array", avs_dx_adios%iglob1, ierr)

  call adios_set_path(adios_handle, &
      "elements_faces/num_ibool_AVS_DX_iglob2", ierr)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      nspec, sizeprocs)
  call adios_write(adios_handle, "array", avs_dx_adios%iglob2, ierr)

  call adios_set_path(adios_handle, &
      "elements_faces/num_ibool_AVS_DX_iglob3", ierr)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      nspec, sizeprocs)
  call adios_write(adios_handle, "array", avs_dx_adios%iglob3, ierr)

  call adios_set_path(adios_handle, &
      "elements_faces/num_ibool_AVS_DX_iglob4", ierr)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      nspec, sizeprocs)
  call adios_write(adios_handle, "array", avs_dx_adios%iglob4, ierr)


  if(ISOTROPIC_3D_MANTLE) then
    call adios_set_path(adios_handle, "elements_faces/dvp", ierr)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        nspec, sizeprocs)
    call adios_write(adios_handle, "array", avs_dx_adios%dvp, ierr)
    call adios_set_path(adios_handle, "elements_faces/dvs", ierr)
    call write_1D_global_array_adios_dims(adios_handle, myrank, &
        nspec, sizeprocs)
    call adios_write(adios_handle, "array", avs_dx_adios%dvs, ierr)
  endif

end subroutine write_AVS_DX_global_faces_data_adios

!===============================================================================
subroutine free_AVS_DX_global_faces_data_adios(myrank, avs_dx_adios, &
    ISOTROPIC_3D_MANTLE)
  implicit none
  !--- Arguments
  integer, intent(in) :: myrank
  type(avs_dx_global_faces_t), intent(inout) :: avs_dx_adios
  logical ISOTROPIC_3D_MANTLE
  !--- Variables
  !--- Variables
  integer :: ierr

  deallocate(avs_dx_adios%x_adios, stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error deallocating x_adios.")
  deallocate(avs_dx_adios%y_adios, stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error deallocating y_adios.")
  deallocate(avs_dx_adios%z_adios, stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error deallocating z_adios.")

  deallocate(avs_dx_adios%idoubling, stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, &
      "Error deallocating num_ibool_AVS_DX_iglob1.")
  deallocate(avs_dx_adios%iglob1, stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, &
      "Error deallocating num_ibool_AVS_DX_iglob1.")
  deallocate(avs_dx_adios%iglob2, stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, &
      "Error deallocating num_ibool_AVS_DX_iglob2.")
  deallocate(avs_dx_adios%iglob3, stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, &
      "Error deallocating num_ibool_AVS_DX_iglob3.")
  deallocate(avs_dx_adios%iglob4, stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, &
      "Error deallocating num_ibool_AVS_DX_iglob4.")

  if(ISOTROPIC_3D_MANTLE) then
    deallocate(avs_dx_adios%dvp, stat=ierr)
    if (ierr /= 0) call exit_MPI(myrank, &
        "Error deallocating dvp.")
    deallocate(avs_dx_adios%dvs, stat=ierr)
    if (ierr /= 0) call exit_MPI(myrank, &
        "Error deallocating dvs.")
  endif

  avs_dx_adios%npoin = 0
  avs_dx_adios%nspecface = 0
end subroutine free_AVS_DX_global_faces_data_adios

end module
